"""Read a CosmoSIS maximum-likelihood fit and estimate local H0 errors.

The optimizer's inverse Hessian approximates the covariance near the best fit.
Unlike MCMC quantiles, these symmetric errors rely on a locally Gaussian fit.
"""
import numpy as np

from chain_io import load_chain


def summarize_maxlike(chain_path, covariance_path):
    """Return fitted parameters and the Hessian-based H0 estimate/error."""
    frame = load_chain(chain_path)
    if len(frame) != 1:
        raise ValueError(f"Expected one maximum-likelihood row, found {len(frame)}")
    row = frame.iloc[0]

    varied_names = [
        name for name in frame.columns
        if name.lower() in {
            "cosmological_parameters--omega_m",
            "cosmological_parameters--h0",
            "supernova_params--m",
        }
    ]
    covariance = np.loadtxt(covariance_path, ndmin=2)
    if covariance.shape != (len(varied_names), len(varied_names)):
        raise ValueError(
            f"Covariance shape {covariance.shape} does not match {len(varied_names)} fitted parameters"
        )
    h_name = next(name for name in varied_names if name.lower().endswith("--h0"))
    h_index = varied_names.index(h_name)
    h_variance = float(covariance[h_index, h_index])
    if not np.isfinite(h_variance) or h_variance <= 0:
        raise ValueError(f"Invalid inverse-Hessian variance for h: {h_variance}")

    m_name = next(name for name in varied_names if name.lower().endswith("--m"))
    omega_name = next(name for name in varied_names if name.lower().endswith("--omega_m"))
    return {
        "H0": float(row["H0"]),
        "H0_sigma": 100.0 * np.sqrt(h_variance),
        "omega_m": float(row[omega_name]),
        "M": float(row[m_name]),
        "loglike": float(row["like"]),
        "varied_names": varied_names,
        "covariance": covariance,
    }


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("chain", nargs="?", default="output/pantheon_plus_shoes_maxlike.txt")
    parser.add_argument("--covariance", default="output/pantheon_plus_shoes_maxlike_covmat.txt")
    args = parser.parse_args()
    for key, value in summarize_maxlike(args.chain, args.covariance).items():
        if key not in {"varied_names", "covariance"}:
            print(f"{key}: {value}")
