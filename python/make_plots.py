"""Create compact, journal-style figures from the directly read chain.

The trace uses muted hues only to distinguish walkers. Posterior and
comparison figures are monochrome and show the retained samples themselves.
Each figure is saved as a 300-dpi PNG.
"""
import argparse
import importlib.util
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from analyze_chain import summarize
from analyze_maxlike import summarize_maxlike


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
    """Save one figure as a high-resolution PNG."""
    fig.savefig(out / f"{name}.png", dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def clean_axes(*axes):
    """Use restrained framing so the plotted data carry visual emphasis."""
    for axis in axes:
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
        axis.tick_params(top=False, right=False)


def make_pantheon_hubble_diagram(maxlike, out: Path):
    """Plot every likelihood-selected SN with its covariance diagonal error.

    Cepheid calibrators are plotted as measured distance moduli; the cosmology
    curve is compared only with Hubble-flow SNe, matching the likelihood logic.
    """
    project_root = Path(__file__).resolve().parents[1]
    library_root = project_root.parent
    likelihood_dir = library_root / "likelihood" / "pantheon_plus"
    data = pd.read_csv(likelihood_dir / "Pantheon+SH0ES.dat", sep=r"\s+")
    calibrator_all = data["IS_CALIBRATOR"].astype(bool).to_numpy()
    selected = (data["zHD"].to_numpy() > 0.01) | calibrator_all

    # Read the same compressed covariance that the CosmoSIS likelihood uses.
    covariance_module = likelihood_dir / "pantheon_covariance_io.py"
    spec = importlib.util.spec_from_file_location("pantheon_covariance_io", covariance_module)
    covio = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(covio)
    covariance = covio.read_compressed(
        str(likelihood_dir / "Pantheon+SH0ES_STAT+SYS.cov_compressed.gz")
    )
    covariance = covariance[np.ix_(selected, selected)]
    sigma = np.sqrt(np.diag(covariance))

    selected_data = data.loc[selected].reset_index(drop=True)
    is_calibrator = selected_data["IS_CALIBRATOR"].astype(bool).to_numpy()
    z_hd = selected_data["zHD"].to_numpy()
    z_hel = selected_data["zHEL"].to_numpy()
    observed_mu = np.empty(len(selected_data))
    observed_mu[is_calibrator] = selected_data.loc[is_calibrator, "CEPH_DIST"]
    observed_mu[~is_calibrator] = (
        selected_data.loc[~is_calibrator, "m_b_corr"].to_numpy() - maxlike["M"]
    )

    # Match the Pantheon+ module's luminosity-distance redshift factors.
    z_grid = np.linspace(0.0, max(2.4, float(z_hd.max()) * 1.01), 5000)
    e_z = np.sqrt(maxlike["omega_m"] * (1.0 + z_grid) ** 3 + (1.0 - maxlike["omega_m"]))
    from scipy.integrate import cumulative_trapezoid
    comoving_distance = (299792.458 / maxlike["H0"]) * cumulative_trapezoid(
        1.0 / e_z, z_grid, initial=0.0
    )
    hf = ~is_calibrator
    dc_hf = np.interp(z_hd[hf], z_grid, comoving_distance)
    # The likelihood multiplies D_A by (1+zHD)(1+zHEL); in flat LCDM,
    # D_A=D_C/(1+zHD), leaving (1+zHEL)D_C for its luminosity-distance term.
    model_mu = 5.0 * np.log10((1.0 + z_hel[hf]) * dc_hf) + 25.0
    residual = observed_mu[hf] - model_mu

    fig, (ax, ax_res) = plt.subplots(
        2, 1, figsize=(7.2, 5.3), sharex=True, constrained_layout=True,
        gridspec_kw={"height_ratios": [2.5, 1.0], "hspace": 0.04},
    )
    ax.errorbar(z_hd[hf], observed_mu[hf], yerr=sigma[hf], fmt=".", ms=2.0,
                color="0.12", ecolor="0.45", elinewidth=0.28, capsize=0,
                alpha=0.30, label="Hubble-flow SNe")
    ax.errorbar(z_hd[is_calibrator], observed_mu[is_calibrator], yerr=sigma[is_calibrator],
                fmt="^", ms=3.2, color="0.55", ecolor="0.65", elinewidth=0.45,
                capsize=1.3, alpha=0.8, label="Cepheid calibrators")
    z_curve = np.geomspace(max(0.01, float(z_hd[hf].min())), float(z_hd[hf].max()), 900)
    dc_curve = np.interp(z_curve, z_grid, comoving_distance)
    # Use zHEL≈zHD for the smooth display curve; the likelihood evaluates the
    # exact per-supernova pair above, and the difference is small in this range.
    mu_curve = 5.0 * np.log10((1.0 + z_curve) * dc_curve) + 25.0
    ax.plot(z_curve, mu_curve, color="black", lw=1.35,
            label=fr"Flat $\Lambda$CDM maximum likelihood ($H_0={maxlike['H0']:.2f}$)")
    ax.set_xscale("log")
    ax.set_ylabel(r"Distance modulus $\mu$ [mag]")
    ax.set_title("Pantheon+SH0ES data and the fitted Hubble diagram")
    ax.legend(frameon=False, loc="lower right", fontsize=7.0)
    ax.text(0.015, 0.98,
            f"Selected: {hf.sum():,} Hubble-flow SNe + {is_calibrator.sum():,} Cepheid calibrators\n"
            f"Maximum likelihood: H₀={maxlike['H0']:.2f}, Ωₘ={maxlike['omega_m']:.3f}, M={maxlike['M']:.3f}",
            transform=ax.transAxes, ha="left", va="top", fontsize=7.0)

    ax_res.errorbar(z_hd[hf], residual, yerr=sigma[hf], fmt=".", ms=2.0,
                    color="0.12", ecolor="0.50", elinewidth=0.28, capsize=0,
                    alpha=0.30)
    ax_res.axhline(0.0, color="black", lw=0.8, ls="--")
    ax_res.set_xscale("log")
    ax_res.set_xlabel("CMB-frame redshift $z_{HD}$")
    ax_res.set_ylabel(r"Residual $\Delta\mu$ [mag]")
    ax_res.set_xlim(0.008, 2.6)
    ax_res.text(0.99, 0.96, "Hubble-flow SNe only", transform=ax_res.transAxes,
                ha="right", va="top", fontsize=7.0, color="0.3")
    clean_axes(ax, ax_res)
    save_figure(fig, out, "pantheon_hubble_diagram")


def make_plots(chain_path="output/pantheon_plus_shoes_chain.txt", out_dir="plots", walkers=8,
               burn_fraction=0.40,
               maxlike_chain_path="output/pantheon_plus_shoes_maxlike.txt",
               covariance_path="output/pantheon_plus_shoes_maxlike_covmat.txt"):
    """Build SN data, posterior, and method/reference comparison figures."""
    apply_journal_style()
    _, chain, samples, grid, density, result = summarize(chain_path, burn_fraction)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    maxlike = None
    if Path(maxlike_chain_path).is_file() and Path(covariance_path).is_file():
        maxlike = summarize_maxlike(maxlike_chain_path, covariance_path)
        make_pantheon_hubble_diagram(maxlike, out)
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
    fig, ax = plt.subplots(figsize=(4.35, 2.85), constrained_layout=True)
    ax.hist(samples, bins="fd", density=True, color="0.88", edgecolor="black",
            linewidth=0.55, label="MCMC samples")
    ax.plot(grid, density, color="black", lw=1.5, label="Gaussian KDE")
    ax.axvspan(result["q16"], result["q84"], color="0.70", alpha=0.25,
               label="central 68% interval")
    ax.axvline(result["median"], color="black", ls="--", lw=0.9,
               label=f"median = {result['median']:.2f}")
    if maxlike is not None:
        ax.axvline(maxlike["H0"], color="0.25", ls="-.", lw=0.9,
                   label="maximum-likelihood fit")
    rug_y = -0.035 * float(np.max(density))
    ax.plot(samples, np.full_like(samples, rug_y), "|", color="0.25", markersize=3.3,
            markeredgewidth=0.35, alpha=0.35)
    ax.set_ylim(bottom=3.0 * rug_y)
    ax.set(xlabel=H0_LABEL, ylabel=r"Density [Mpc s km$^{-1}$]",
           title=fr"Marginalized $H_0$ posterior ($N={len(samples):,}$)")
    ax.minorticks_on()
    ax.tick_params(which="minor", length=2.2)
    ax.legend(frameon=False, loc="center left", bbox_to_anchor=(1.01, 0.5),
              fontsize=6.5, handlelength=1.6)
    clean_axes(ax)
    save_figure(fig, out, "h0_posterior")

    # Compare the MCMC posterior, deterministic fit, and published references.
    # SH0ES is contextual because its calibration enters the fitted likelihood.
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
    ax.legend(frameon=False, loc="upper left", fontsize=6.8, handlelength=1.6)

    # Markers separate provenance; all intervals are plotted in monochrome.
    comparisons = [
        (3, 67.4, 0.5, 0.5, "s", "Planck 2018 published", "$67.4\\pm0.5$"),
        (2, 73.04, 1.04, 1.04, "D", "SH0ES 2022 published", "$73.04\\pm1.04$"),
        (0, result["median"], result["median"] - result["q16"],
         result["q84"] - result["median"], "o", "This analysis",
         f"${result['median']:.2f}^{{+{result['q84']-result['median']:.2f}}}_{{-{result['median']-result['q16']:.2f}}}$"),
    ]
    if maxlike is not None:
        # BFGS's inverse Hessian estimates local covariance near the optimum.
        comparisons.insert(1, (1, maxlike["H0"], maxlike["H0_sigma"], maxlike["H0_sigma"], "^",
                               "CosmoSIS maximum likelihood (Hessian)",
                               f"${maxlike['H0']:.2f}\\pm{maxlike['H0_sigma']:.2f}$"))
    for y, center, lower, upper, marker, _, _ in comparisons:
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
    parser.add_argument("--maxlike-chain", default="output/pantheon_plus_shoes_maxlike.txt")
    parser.add_argument("--covariance", default="output/pantheon_plus_shoes_maxlike_covmat.txt")
    args = parser.parse_args()
    print(make_plots(args.chain, args.out_dir, args.walkers, args.burn_fraction,
                     args.maxlike_chain, args.covariance))
