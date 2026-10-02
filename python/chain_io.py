"""Read and prepare a CosmoSIS text chain for the analysis scripts.

This module isolates file parsing and the dimensionless-h to physical-H0
conversion from posterior calculations and plotting.
"""
from pathlib import Path

import pandas as pd


def load_chain(path: str | Path) -> pd.DataFrame:
    """Load the whitespace chain, retain its names, and add H0 in km/s/Mpc.

    CosmoSIS writes parameter names on the first line with a leading ``#``;
    pandas' normal comment handling would discard that header, so it is read
    separately before the remaining commented metadata and numeric rows.
    """
    with Path(path).open(encoding="utf-8") as stream:
        header = stream.readline().strip().lstrip("#").split()
    frame = pd.read_csv(path, sep=r"\s+", comment="#", header=None,
                        names=header, engine="python")
    h_column = next((name for name in frame.columns if name.lower().endswith("--h0")), None)
    if h_column is None:
        raise ValueError(f"No cosmological_parameters--h0 column in {path}; found {list(frame.columns)}")
    # CosmoSIS stores h = H0 / 100, with H0 measured in km s^-1 Mpc^-1.
    frame["H0"] = 100.0 * pd.to_numeric(frame[h_column], errors="coerce")
    return frame
