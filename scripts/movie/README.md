# Point-vortex movie tool

`make_vortex_movie.py` reads the `vortices.csv` trajectory written by the C++
solver and streams it into an MP4. It groups rows by `time`, so the entire
trajectory does not need to fit in memory. The extra `index`, `u`, and `v`
columns may be present and are ignored. The number of vortices may vary between
saved frames when dipole removal is enabled.

## Requirements

- Python 3.10 or newer
- Matplotlib
- FFmpeg available as `ffmpeg`
- A LaTeX installation available as `latex`

For example, install Matplotlib with your Python environment's package manager:

```bash
python3 -m pip install matplotlib
```

## Examples

Run these commands from the repository root:

```bash
# Infinite plane; the script determines fixed limits from the selected frames.
python3 scripts/movie/make_vortex_movie.py data/vortices.csv \
    --geometry infinite --output data/infinite.mp4

# Square periodic box. Coordinates are wrapped for display only.
python3 scripts/movie/make_vortex_movie.py data/vortices.csv \
    --geometry periodic --box-length 2.0 --output data/periodic.mp4

# Circular disk.
python3 scripts/movie/make_vortex_movie.py data/vortices.csv \
    --geometry disk --radius 1.0 --output data/disk.mp4
```

`--start`, `--stop`, and `--stride` select output-frame indices. For example,
this renders frames 100 through 499, taking every second frame:

```bash
python3 scripts/movie/make_vortex_movie.py data/vortices.csv \
    --geometry disk --start 100 --stop 500 --stride 2 --output data/excerpt.mp4
```

For an infinite-plane movie, `--xlim MIN MAX` and `--ylim MIN MAX` override the
automatic limits. Other useful controls include `--fps`, `--dpi`,
`--marker-size`, and `--title`. Run the script with `--help` for the complete
command-line reference.

The script validates numeric CSV data and axis ranges before rendering, and
reports missing LaTeX or FFmpeg executables explicitly.

The trajectory must be ordered by output time, as produced by PointVortex.
Positive and negative vortices are drawn in red and blue respectively. Rows
with zero circulation are not displayed. Plot text is rendered with LaTeX and
Computer Modern serif fonts.
