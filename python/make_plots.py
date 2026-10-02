"""Create compact, journal-style figures from the directly read chain.

The trace uses muted hues only to distinguish walkers. Posterior and
comparison figures are monochrome and show the retained samples themselves.
Each figure is saved as a 300-dpi PNG and a vector PDF.
"""
import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from analyze_chain import summarize


H0_LABEL = r"$H_0$ [km s$^{-1}$ Mpc$^{-1}$]"


def apply_journal_style():
    """Set compact, consistent typography and line weights for print figures."""
    mpl.rcParams.update({
        "font.family": "serif",
        "font.serif": ["DejaVu Serif"],
        "mathtext.fontset": "stix",
        "font.size": 8,
        "axes.labelsize": 9,
        "axes.titlesize": 9,
        "axes.linewidth": 0.7,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
        "legend.fontsize": 7.5,
        "lines.linewidth": 1.0,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.facecolor": "white",
    })


def save_figure(fig, out: Path, name: str, dpi: int = 300):
    """Save one figure in a crisp raster format and an editable vector format."""
    fig.savefig(out / f"{name}.png", dpi=dpi, bbox_inches="tight", facecolor="white")
    fig.savefig(out / f"{name}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def clean_axes(*axes):
    """Use restrained framing so the plotted data carry visual emphasis."""
    for axis in axes:
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
        axis.tick_params(top=False, right=False)


def make_plots(chain_path="output/pantheon_plus_shoes_chain.txt", out_dir="plots", walkers=8,
               burn_fraction=0.40, planck_chain_path="output/planck_lite_chain.txt"):
    """Build figures comparing the local chain, a Planck-lite chain, and references."""
    apply_journal_style()
    _, chain, samples, grid, density, result = summarize(chain_path, burn_fraction)
    planck_result = None
    if planck_chain_path and Path(planck_chain_path).is_file():
        from analyze_chain import summarize as summarize_chain
        _, _, _, _, _, planck_result = summarize_chain(
            planck_chain_path, burn_fraction
        )
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    h0 = chain["H0"].to_numpy()

    # CosmoSIS writes rows in walker groups for each iteration. Refuse a
    # partial group rather than silently truncating or mislabelling a trace.
    if walkers <= 0 or len(h0) % walkers:
        raise ValueError(f"Chain has {len(h0)} rows, which cannot be split into {walkers} walkers")
    steps = len(h0) // walkers
    traces = h0.reshape(steps, walkers)

    # Trace: every stored point is marked; muted colors help track walkers.
    fig, ax = plt.subplots(figsize=(7.2, 3.3), constrained_layout=True)
    muted_palette = ["#71869A", "#9A8876", "#789184", "#9A7F86",
                     "#858395", "#A29266", "#6F9293", "#9B8D9E"]
    line_styles = ["-", "--", ":", "-."]
    iteration = np.arange(steps)
    for i in range(walkers):
        ax.plot(iteration, traces[:, i], color=muted_palette[i % len(muted_palette)],
                ls=line_styles[i % len(line_styles)], lw=0.65, marker=".", markersize=1.7,
                alpha=0.82)
    ax.axvline(steps * burn_fraction, color="black", ls=":", lw=0.9,
               label=f"burn-in cut ({burn_fraction:.0%})")
    ax.axhline(result["median"], color="black", ls="--", lw=0.9,
               label="post-burn-in median")
    ax.set(xlabel="MCMC iteration", ylabel=H0_LABEL,
           xlim=(0, steps - 1), title="Pantheon+SH0ES sampling trace")
    ax.minorticks_on()
    ax.tick_params(which="minor", length=2.2)
    ax.legend(frameon=False, loc="upper left", ncol=2, handlelength=2.0)
    ax.text(0.99, 0.03, f"{walkers} walkers · {steps} iterations",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=7.5)
    clean_axes(ax)
    save_figure(fig, out, "h0_trace")

    # Posterior: the histogram preserves the empirical distribution; the KDE
    # is a visual guide. Rug marks show every retained sample, not synthetic data.
    fig, ax = plt.subplots(figsize=(3.55, 2.85), constrained_layout=True)
    ax.hist(samples, bins="fd", density=True, color="0.88", edgecolor="black",
            linewidth=0.55, label="MCMC samples")
    ax.plot(grid, density, color="black", lw=1.5, label="Gaussian KDE")
    ax.axvspan(result["q16"], result["q84"], color="0.70", alpha=0.25,
               label="central 68% interval")
    ax.axvline(result["median"], color="black", ls="--", lw=0.9,
               label=f"median = {result['median']:.2f}")
    rug_y = -0.035 * float(np.max(density))
    ax.plot(samples, np.full_like(samples, rug_y), "|", color="0.25", markersize=3.3,
            markeredgewidth=0.35, alpha=0.35)
    ax.set_ylim(bottom=3.0 * rug_y)
    ax.set(xlabel=H0_LABEL, ylabel=r"Density [Mpc s km$^{-1}$]",
           title=fr"Marginalized $H_0$ posterior ($N={len(samples):,}$)")
    ax.minorticks_on()
    ax.tick_params(which="minor", length=2.2)
    ax.legend(frameon=False, loc="upper right", fontsize=6.8, handlelength=1.6)
    clean_axes(ax)
    save_figure(fig, out, "h0_posterior")

    # Comparison: the upper panel shows this chain's density and sample rug;
    # the lower forest plot puts published and chain-derived intervals on one
    # physical scale. SH0ES is contextual and is included in this likelihood.
    fig, (ax, ax_ref) = plt.subplots(
        2, 1, figsize=(7.2, 4.8), sharex=True, constrained_layout=True,
        gridspec_kw={"height_ratios": [2.5, 1.55], "hspace": 0.04},
    )
    ax.plot(grid, density, color="black", lw=1.5, label="Pantheon+SH0ES posterior")
    ax.axvspan(67.4 - 0.5, 67.4 + 0.5, color="0.75", alpha=0.32,
               label="Planck 2018 published 68% interval")
    ax.axvline(67.4, color="0.35", ls="--", lw=0.9)
    ax.axvspan(result["q16"], result["q84"], color="0.78", alpha=0.32)
    ax.axvline(result["median"], color="black", ls="--", lw=0.9)
    rug_y = -0.035 * float(np.max(density))
    ax.plot(samples, np.full_like(samples, rug_y), "|", color="0.25", markersize=3.1,
            markeredgewidth=0.35, alpha=0.35)
    ax.set_ylim(bottom=3.0 * rug_y)
    ax.set_ylabel(r"Density [Mpc s km$^{-1}$]")
    ax.set_title("Local distance-ladder and Planck CMB constraints")
    ax.legend(frameon=False, loc="upper right", fontsize=6.8, handlelength=1.6)

    # Markers separate provenance; all intervals are plotted in monochrome.
    comparisons = [
        (3, 67.4, 0.5, 0.5, "s", "Planck 2018 published", "$67.4\\pm0.5$"),
        (2, 73.04, 1.04, 1.04, "D", "SH0ES 2022 published", "$73.04\\pm1.04$"),
        (0, result["median"], result["median"] - result["q16"],
         result["q84"] - result["median"], "o", "This analysis",
         f"${result['median']:.2f}^{{+{result['q84']-result['median']:.2f}}}_{{-{result['median']-result['q16']:.2f}}}$"),
    ]
    if planck_result is not None:
        # This run is too short for a defensible Planck-lite credible interval;
        # show its exploratory median only and keep the published error bar.
        comparisons.insert(1, (1, planck_result["median"], 0.0, 0.0, "^",
                               "Planck-lite chain median†", f"{planck_result['median']:.2f}"))
    for y, center, lower, upper, marker, _, _ in comparisons:
        if marker == "^":
            ax_ref.plot(center, y, marker=marker, color="0.4", markersize=4.0, zorder=3)
        else:
            ax_ref.errorbar(center, y, xerr=np.array([[lower], [upper]]), fmt=marker,
                            color="black", ecolor="black", elinewidth=1.0, capsize=2.8,
                            capthick=0.85, markersize=4.0, zorder=3)
    y_positions = [item[0] for item in comparisons]
    ax_ref.set_yticks(y_positions, [item[5] for item in comparisons])
    ax_ref.set_ylim(-0.55, max(y_positions) + 0.55)
    ax_ref.set_xlabel(H0_LABEL)
    ax_ref.set_xlim(66.2, 77.3)
    ax_ref.grid(axis="x", color="0.86", lw=0.45)
    ax_ref.set_axisbelow(True)
    # Put the numeric interval beside each bar, making the error bars readable
    # in print without relying on a legend or color coding.
    for y, center, lower, upper, marker, label, value in comparisons:
        if marker == "^":
            ax_ref.annotate(f"{value} (median only)", (center, y), xytext=(5, 0),
                            textcoords="offset points", ha="left", va="center", fontsize=7)
            continue
        ax_ref.annotate(value, (center + upper, y), xytext=(4, 0),
                        textcoords="offset points", ha="left", va="center", fontsize=7)
    for axis in (ax, ax_ref):
        axis.minorticks_on()
        axis.tick_params(which="minor", length=2.0)
    clean_axes(ax, ax_ref)
    save_figure(fig, out, "h0_planck_comparison")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("chain", nargs="?", default="output/pantheon_plus_shoes_chain.txt")
    parser.add_argument("--out-dir", default="plots")
    parser.add_argument("--walkers", type=int, default=8)
    parser.add_argument("--burn-fraction", type=float, default=0.40)
    parser.add_argument("--planck-chain", default="output/planck_lite_chain.txt",
                        help="Optional Planck-lite CosmoSIS chain for a direct posterior comparison")
    args = parser.parse_args()
    print(make_plots(args.chain, args.out_dir, args.walkers, args.burn_fraction,
                     args.planck_chain))
