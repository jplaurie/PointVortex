# PointVortex analysis notebook

`point_vortex_analysis.ipynb` reads a trajectory CSV and diagnostics CSV produced by any solver
backend. It writes two vector PDF figures:

- `diagnostics_evolution.pdf` shows the Hamiltonian, geometry-relevant impulses, conservation
  drift, and cumulative dipole events;
- `vortex_configurations.pdf` shows evenly spaced physical-space snapshots, including the
  periodic square or disk boundary when applicable.

Launch Jupyter from the repository root:

```bash
jupyter lab scripts/analysis/point_vortex_analysis.ipynb
```

Edit the settings in the first code cell to select `infinite`, `periodic`, or `disk`, specify the
domain size, input CSV paths, number of snapshots, and output directory. The defaults read
`data/vortices.csv` and `data/diagnostics.csv` and write PDFs under `data/figures/`.

The notebook requires Python 3, NumPy, Matplotlib, and Jupyter. It does not require pandas or a
LaTeX installation. It searches upward for the repository root, so it can also be launched from
the notebook directory.
