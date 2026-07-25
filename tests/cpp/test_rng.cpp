// Tests for the deterministic RNG.
//
// These exist because reproducibility is a load-bearing property of this
// project, not a nicety. Benchmarks are reported with confidence intervals over
// seeds, and golden trajectory fixtures are compared across Linux and Windows.
// If the RNG drifts between platforms every one of those comparisons is void.

#include <gtest/gtest.h>

#include <cmath>
#include <random>
#include <vector>

#include "core/types.hpp"

using namespace aal;

TEST(Rng, Mt19937_64MatchesTheStandardMandatedValue) {
  // The C++ standard pins this exactly: the 10000th consecutive invocation of a
  // default-constructed mt19937_64 shall produce 9981545732273789042. If this
  // ever fails, the standard library is non-conforming and no cross-platform
  // reproducibility claim in this repository holds.
  std::mt19937_64 eng;
  std::uint64_t v = 0;
  for (int i = 0; i < 10000; ++i) v = eng();
  EXPECT_EQ(v, 9981545732273789042ULL);
}

TEST(Rng, IsReproducibleAcrossInstances) {
  Rng a(12345), b(12345);
  for (int i = 0; i < 1000; ++i) {
    EXPECT_EQ(a.nextBits(), b.nextBits());
  }
  Rng c(12345), d(12345);
  for (int i = 0; i < 1000; ++i) {
    EXPECT_DOUBLE_EQ(c.uniform(), d.uniform());
  }
}

TEST(Rng, UniformStaysInRange) {
  Rng r(99);
  for (int i = 0; i < 200000; ++i) {
    const Scalar u = r.uniform();
    // Half-open [0,1): exactly 1.0 must never appear, or mapping into an array
    // index by multiplication overflows by one.
    EXPECT_GE(u, 0.0);
    EXPECT_LT(u, 1.0);
  }
}

TEST(Rng, UniformIsApproximatelyFlat) {
  Rng r(4);
  constexpr int kBins = 20;
  constexpr int kDraws = 400000;
  std::vector<int> hist(kBins, 0);
  for (int i = 0; i < kDraws; ++i) hist[static_cast<int>(r.uniform() * kBins)]++;
  const double expected = double(kDraws) / kBins;
  // Chi-square with 19 degrees of freedom: the 99.9th percentile is 43.8.
  double chi2 = 0;
  for (int c : hist) chi2 += (c - expected) * (c - expected) / expected;
  EXPECT_LT(chi2, 43.8) << "uniform draws are not flat, chi2 = " << chi2;
}

TEST(Rng, NormalHasCorrectMomentsAndIsReproducible) {
  Rng r(2024);
  constexpr int kN = 400000;
  double sum = 0, sumsq = 0;
  for (int i = 0; i < kN; ++i) {
    const double x = r.normal();
    sum += x;
    sumsq += x * x;
  }
  const double mean = sum / kN;
  const double var = sumsq / kN - mean * mean;
  // Standard error of the mean is 1/sqrt(N) = 0.00158, so 5 sigma is 0.0079.
  EXPECT_NEAR(mean, 0.0, 0.008);
  EXPECT_NEAR(var, 1.0, 0.02);

  Rng a(77), b(77);
  for (int i = 0; i < 500; ++i) EXPECT_DOUBLE_EQ(a.normal(), b.normal());
}

TEST(Rng, UnitBallSamplesAreInsideAndRadiallyUniform) {
  Rng r(5);
  constexpr int kN = 200000;
  int inner_octant = 0;
  for (int i = 0; i < kN; ++i) {
    const Vec3 p = r.uniformInUnitBall();
    ASSERT_LE(p.norm(), 1.0 + 1e-12);
    // For a uniform ball the fraction within radius 0.5 is (0.5)^3 = 0.125.
    if (p.norm() <= 0.5) ++inner_octant;
  }
  const double frac = double(inner_octant) / kN;
  EXPECT_NEAR(frac, 0.125, 0.005);
}
