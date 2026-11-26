"""Why skiver's Illumina hazard fit is distorted: simulate the quantity skiver fits.

Skiver observes h(t) only for t in (k, k+v] -- and v0.3.2 narrows that to [k+3, k+v-2]
via --ignore-{smallest,largest}-t 2 -- then extrapolates a discrete Weibull back to
t = 1 to report the per-base error rate.  A (k,v)-mer spans k+v = 34 bases, a third of
a 101 bp Illumina read, so position-within-value aliases onto position-within-read:
larger t sits systematically later in the read, where the error rate is higher.

This script computes h(t) exactly as skiver does -- conditioned on a clean k-base key,
aggregated over every (k,v)-mer offset in the read -- for a read set with
  (a) read-level frailty (heterogeneous per-read rate; what beta < 1 reflects), and
  (b) a positional error ramp,
then runs skiver's clog-log fit on the window and compares the extrapolated h(1) with
the true mean per-base rate.

Run with the measured ramp (default; taken from
`skiver analyze` -> summary_read_position.csv on ERR2935851, which is derived from
consensus disagreement, never from reported quality):

    python3 -I docs/skiver_illumina_sim.py [ramp.npy]
"""

import sys

import numpy as np

K = V = 17
KV = K + V
IGNORE_LO, IGNORE_HI = 2, 2  # skiver v0.3.2 defaults

# Measured per-cycle error ramp, ERR2935851 R1 (B. subtilis, HiSeq 1500, 101 bp), as
# (relative rate, number of cycles), normalised to mean 1.0 by ramp_profile().  Taken from
# `skiver analyze` -> summary_read_position.csv, i.e. consensus disagreement, never from the
# reported quality scores.  6.4x from read start to read end.
MEASURED_RAMP = np.repeat([0.16, 0.32, 0.55, 1.54, 1.68, 1.61], [10, 15, 25, 25, 15, 11])


def ramp_profile(length, shape, strength):
    """Per-position multiplier, mean 1.0, stretched to `length`."""
    if shape is None:
        return np.ones(length)
    base = np.interp(np.linspace(0, 1, length), np.linspace(0, 1, len(shape)), shape)
    base = base / base.mean()
    # strength=1 reproduces the profile as measured; <1 flattens it toward uniform
    return 1.0 + strength * (base - 1.0)


def h_of_t_as_skiver_sees_it(length, mean_rate, profile, frailty, weights):
    """Exact h(t), t = 1..KV, conditioned on a clean key and pooled over offsets."""
    base = mean_rate * profile
    offsets = np.arange(0, length - KV + 1)
    surv = np.zeros(KV + 1)  # sum of P(KV-mer positions 1..t all clean)
    num = np.zeros(KV + 1)  # sum of P(first mismatch exactly at t)
    total_w = 0.0
    for f, w in zip(frailty, weights, strict=True):
        e_read = np.clip(base * f, 0.0, 1.0)
        for s in offsets:
            e = e_read[s : s + KV]
            clean = np.cumprod(1.0 - e)
            pre = np.concatenate(([1.0], clean[:-1]))
            surv[1:] += w * clean
            num[1:] += w * pre * e
            total_w += w
    h = np.zeros(KV + 1)
    for t in range(1, KV + 1):
        denom = surv[t - 1] if t >= 2 else total_w
        h[t] = num[t] / denom
    true_mean = sum(w * np.mean(np.clip(base * f, 0.0, 1.0)) for f, w in zip(frailty, weights, strict=True)) / sum(
        weights
    )
    return h, true_mean


def fit_window(h, ignore_lo=IGNORE_LO, ignore_hi=IGNORE_HI):
    """Skiver's clog-log regression over its (narrowed) observable window."""
    t = np.arange(K + 1 + ignore_lo, KV + 1 - ignore_hi)
    y = np.log(-np.log(1.0 - h[t]))
    slope, intercept = np.polyfit(np.log(t), y, 1)
    beta = slope + 1.0
    lam = np.exp(intercept) / beta
    return lam, beta, 1.0 - np.exp(-lam)


def main(argv):
    shape = MEASURED_RAMP
    label = "measured (ERR2935851, 6.4x end/start)"
    if len(argv) > 1:
        pos_err = np.load(argv[1])
        shape = pos_err[1]
        label = f"measured from {argv[1]}"
    print(f"ramp profile: {label}")
    print("frailty: 3-class mixture, relative rates 0.3 / 1.0 / 3.0 at 0.5 / 0.4 / 0.1\n")

    frailty, weights = [0.3, 1.0, 3.0], [0.5, 0.4, 0.1]
    flat = ([1.0], [1.0])

    scenarios = [
        ("long read 10kb, flat rate, no frailty", 10000, 0.006, None, 0.0, flat),
        ("long read 10kb, frailty, no ramp", 10000, 0.006, None, 0.0, (frailty, weights)),
        ("101bp, frailty, no ramp", 101, 0.006, None, 0.0, (frailty, weights)),
        ("101bp, frailty, ramp x0.25 strength", 101, 0.006, shape, 0.25, (frailty, weights)),
        ("101bp, frailty, ramp x0.5 strength", 101, 0.006, shape, 0.5, (frailty, weights)),
        ("101bp, frailty, MEASURED ramp", 101, 0.006, shape, 1.0, (frailty, weights)),
        ("151bp, frailty, MEASURED ramp", 151, 0.006, shape, 1.0, (frailty, weights)),
        ("301bp, frailty, MEASURED ramp", 301, 0.006, shape, 1.0, (frailty, weights)),
        ("long read 10kb, frailty, MEASURED ramp", 10000, 0.006, shape, 1.0, (frailty, weights)),
    ]

    print(f"{'scenario':<40}{'true h(1)':>11}{'skiver h1':>11}{'ratio':>8}{'beta_hat':>10}")
    for label, length, rate, shp, strength, (fr, w) in scenarios:
        profile = ramp_profile(length, shp, strength)
        h, true_mean = h_of_t_as_skiver_sees_it(length, rate, profile, fr, w)
        _, beta, h1 = fit_window(h)
        print(f"{label:<40}{true_mean:>11.5f}{h1:>11.5f}{h1 / true_mean:>7.2f}x{beta:>10.3f}")

    # Sensitivity of the backward extrapolation to beta, on skiver's real ERR2935851 window.
    # Sensitivity of the backward extrapolation to beta, on skiver v0.3.2's real window
    # from the FULL ERR2935851 R1 (23.16 M reads, 557x; k = v = 17, so t = 20..32 after the
    # default --ignore-{smallest,largest}-t 2).  skiver reports beta = 0.938, h(1) = 0.004854.
    print("\nSensitivity of h(1) to beta, real window from full ERR2935851 R1 (557x):")
    t = np.arange(20, 33)
    h_real = np.array(
        [
            0.003979,
            0.003895,
            0.003837,
            0.003785,
            0.003762,
            0.003775,
            0.003684,
            0.003733,
            0.003700,
            0.003676,
            0.003733,
            0.003618,
            0.003620,
        ]
    )
    y = np.log(-np.log(1 - h_real))
    slope, intercept = np.polyfit(np.log(t), y, 1)
    beta_ols = slope + 1.0
    h1_ols = 1 - np.exp(-np.exp(intercept) / beta_ols)
    print("  skiver (ridge + Huber): beta=0.938  h(1)=0.004854")
    print(
        f"  plain OLS, same points: beta={beta_ols:.3f}  h(1)={h1_ols:.5f}"
        f"   <- {h1_ols / 0.004854:.2f}x skiver, and ~BEST's 0.0072"
    )
    for beta_fixed in (1.0, 0.938, 0.90, 0.85, 0.80, 0.75, 0.70):
        intercept_f = np.mean(y - (beta_fixed - 1.0) * np.log(t))
        lam = np.exp(intercept_f) / beta_fixed
        print(f"  beta fixed {beta_fixed:.3f}:      h(1)={1 - np.exp(-lam):.5f}")


if __name__ == "__main__":
    main(sys.argv)
