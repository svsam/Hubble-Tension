"""Summarise the marginalized H0 posterior from a CosmoSIS text chain.

The reported interval is a sample percentile interval after a chosen burn-in
cut. It is descriptive of this run and is not a convergence diagnostic.
"""
import argparse
from pathlib import Path

import numpy as np
from scipy.stats import gaussian_kde

from chain_io import load_chain


def summarize(path: str | Path, burn_fraction: float = 0.40) -> tuple:
    """Load samples, trim burn-in, estimate the H0 density, and summarize it."""
    if not 0 <= burn_fraction < 1:
        raise ValueError("burn_fraction must be in [0, 1)")
    frame = load_chain(path)
    chain = frame.loc[np.isfinite(frame["H0"])].copy()
    # The chain is written iteration-by-iteration across walkers, so cutting
    # this fraction of rows removes the same initial fraction of each walker.
    cut = int(len(chain) * burn_fraction)
    samples = chain["H0"].to_numpy()[cut:]
    if len(samples) < 20:
        raise ValueError(f"Only {len(samples)} finite post-burn-in samples; run a longer chain")
    # The 16th/84th percentiles give a central 68% interval; retain its
    # asymmetry rather than assuming the finite chain is exactly Gaussian.
    q16, median, q84 = np.quantile(samples, [0.16, 0.50, 0.84])
    mean = float(np.mean(samples))
    # KDE is used for a legible curve in the plots; quantiles above remain the
    # reported interval so conclusions do not depend on KDE bandwidth.
    kde = gaussian_kde(samples)
    grid = np.linspace(float(np.min(samples)), float(np.max(samples)), 600)
    density = kde(grid)
    # Published literature references stay distinct from a separate Planck-lite
    # chain. SH0ES calibration does enter the Pantheon+SH0ES likelihood.
    local_mu, local_sigma = 73.04, 1.04
    planck_mu, planck_sigma = 67.4, 0.5
    # Approximate Gaussian comparison: use half the project's 68% width as a
    # single sigma. This is illustrative, not a complete consistency test.
    tension = abs(median - planck_mu) / np.sqrt(((q84 - q16) / 2) ** 2 + planck_sigma**2)
    # Keep published reference numbers beside (but distinguishable from) the
    # chain-derived posterior values returned below.
    summary = {
        "raw_rows": len(frame), "post_burn_rows": len(samples), "burn_fraction": burn_fraction,
        "mean": mean, "median": float(median), "q16": float(q16), "q84": float(q84),
        "local_published_h0": local_mu, "local_published_sigma": local_sigma,
        "planck_h0": planck_mu, "planck_sigma": planck_sigma,
        "illustrative_planck_tension_sigma": float(tension),
    }
    return frame, chain, samples, grid, density, summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("chain", nargs="?", default="output/pantheon_plus_shoes_chain.txt")
    parser.add_argument("--burn-fraction", type=float, default=0.40)
    args = parser.parse_args()
    *_, summary = summarize(args.chain, args.burn_fraction)
    for key, value in summary.items():
        print(f"{key}: {value}")
