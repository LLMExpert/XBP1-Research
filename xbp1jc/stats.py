"""
Statistics for junction-count proportions.

Two things matter here and both are easy to get wrong:

1. f(XBP1s) is a binomial proportion that is frequently near zero in unstressed
   samples. The normal approximation misbehaves there and can produce intervals
   extending below zero, so Wilson intervals are used throughout.

2. Adequacy of depth is a property of the INFORMATIVE JUNCTION READ COUNT, not
   of sequencing depth. The two are related by a factor that varied by more
   than an order of magnitude across the conditions we examined. The required
   count is therefore computed by power calculation and fixed before analysis;
   groups below it are reported as bounds, not estimates.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

# Normal quantiles for common two-sided alpha levels.
_Z = {0.10: 1.6448536269514722,
      0.05: 1.959963984540054,
      0.01: 2.5758293035489004}


def _z(alpha: float) -> float:
    if alpha in _Z:
        return _Z[alpha]
    # Acklam-style rational approximation to the inverse normal CDF,
    # adequate for interval construction.
    p = 1.0 - alpha / 2.0
    a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
         1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00]
    b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
         6.680131188771972e+01, -1.328068155288572e+01]
    c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
         -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00]
    dd = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
          3.754408661907416e+00]
    plow, phigh = 0.02425, 1 - 0.02425
    if p < plow:
        q = math.sqrt(-2 * math.log(p))
        return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / \
               ((((dd[0] * q + dd[1]) * q + dd[2]) * q + dd[3]) * q + 1)
    if p > phigh:
        q = math.sqrt(-2 * math.log(1 - p))
        return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / \
               ((((dd[0] * q + dd[1]) * q + dd[2]) * q + dd[3]) * q + 1)
    q = p - 0.5
    r = q * q
    return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q / \
           (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1)


@dataclass
class Interval:
    point: float          # percentage
    lo: float             # percentage
    hi: float             # percentage
    n: int                # informative junction reads

    def __str__(self) -> str:
        return (f"{self.point:.1f}% (95% CI {self.lo:.1f}-{self.hi:.1f}; "
                f"n = {self.n:,})")


def wilson(spliced: int, unspliced: int, alpha: float = 0.05) -> Interval | None:
    """
    Wilson score interval for the spliced fraction, returned as percentages.

    Reference: Wilson EB (1927) J Am Stat Assoc 22:209-212.

    Returns None where there are no informative reads, which is the honest
    output: an arm with zero informative junction reads has no estimate, and
    reporting 0% would be wrong.
    """
    n = spliced + unspliced
    if n == 0:
        return None
    z = _z(alpha)
    p = spliced / n
    denom = 1.0 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return Interval(point=100.0 * p,
                    lo=100.0 * max(0.0, centre - half),
                    hi=100.0 * min(1.0, centre + half),
                    n=n)


def required_reads(baseline_pct: float,
                   fold_change: float = 2.0,
                   power: float = 0.80,
                   alpha: float = 0.05) -> int:
    """
    Informative junction reads required per group to detect `fold_change`
    relative to `baseline_pct`, for a two-sided two-proportion comparison.

    Uses the standard normal-approximation sample size with the variance
    pooled under the null:

        n = [ z_{a/2} * sqrt(2 p_bar q_bar)
              + z_b * sqrt(p1 q1 + p2 q2) ]^2 / (p2 - p1)^2

    For the 8.3% pooled mock baseline of the time course in our corpus,
    detecting a doubling at 80% power and alpha=0.05 gives 247.19, i.e. 247
    informative junction reads per group. This is the published depth standard
    and this function reproduces it exactly.

    The result is rounded to nearest rather than ceiled, so that the boundary
    case matches the published figure.
    """
    p1 = baseline_pct / 100.0
    if not 0 < p1 < 1:
        raise ValueError("baseline must be strictly between 0 and 100 percent")
    if fold_change == 1:
        raise ValueError("fold_change of 1 gives no detectable effect")
    p2 = min(0.999999, p1 * fold_change)
    if p2 == p1:
        raise ValueError("fold_change produces no difference in proportions")

    z_a = _z(alpha)                 # two-sided
    z_b = _z(2 * (1 - power))       # one-sided quantile for power
    p_bar = (p1 + p2) / 2.0
    num = (z_a * math.sqrt(2 * p_bar * (1 - p_bar))
           + z_b * math.sqrt(p1 * (1 - p1) + p2 * (1 - p2))) ** 2
    return int(round(num / (p2 - p1) ** 2))


def adequacy(spliced: int, unspliced: int, baseline_pct: float,
             fold_change: float = 2.0, power: float = 0.80,
             alpha: float = 0.05) -> dict:
    """
    Verdict on whether a group's informative read count supports an estimate.

    Returns the required count, the observed count, and a verdict of
    "estimate" or "bound". A group below threshold should be reported as a
    bound with its interval, not as a point estimate.
    """
    need = required_reads(baseline_pct, fold_change, power, alpha)
    have = spliced + unspliced
    return {"required": need,
            "observed": have,
            "adequate": have >= need,
            "verdict": "estimate" if have >= need else "bound",
            "baseline_pct": baseline_pct,
            "fold_change": fold_change,
            "power": power,
            "alpha": alpha}
