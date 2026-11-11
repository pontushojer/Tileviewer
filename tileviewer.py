"""
Tileviewer: A Streamlit app for visualizing genomic data.
"""

import base64
import io
import zipfile

import numpy as np
import plotly.graph_objects as go
import streamlit as st
from matplotlib import colormaps
from PIL import Image
from scipy import sparse
from scipy.ndimage import maximum_filter


RASTER_SIZE = 900
HOVER_LIMIT = 100_000
MAX_EXACT_INT = 2**53 - 1
FOCUS_RADIUS_BP = 1_000_000


@st.cache_resource(show_spinner="Loading sparse matrix...")
def load_matrix(content):
    matrix = sparse.load_npz(io.BytesIO(content))
    if not sparse.isspmatrix_csr(matrix):
        raise ValueError("The NPZ must contain a SciPy CSR matrix.")
    if not np.issubdtype(matrix.dtype, np.integer):
        raise ValueError("Matrix values must have an integer dtype.")
    if not all(matrix.shape):
        raise ValueError("Both matrix dimensions must be nonzero.")
    matrix.check_format(full_check=True)
    matrix.sort_indices()
    if not matrix.has_canonical_format:
        raise ValueError("Duplicate matrix entries are not supported.")
    if matrix.nnz and np.any(matrix.data < 0):
        raise ValueError("Mapped coordinates must be positive.")
    return matrix


def parse_fai(content):
    names, ends = [], []
    total = 0
    for line_number, line in enumerate(content.decode("utf-8-sig").splitlines(), 1):
        fields = line.split("\t")
        if len(fields) < 2 or not fields[0]:
            raise ValueError(f"Invalid FAI record on line {line_number}.")
        if fields[0] in names:
            raise ValueError(f"Duplicate reference name on line {line_number}.")
        try:
            length = int(fields[1])
        except ValueError as exc:
            raise ValueError(f"Invalid reference length on line {line_number}.") from exc
        if length <= 0:
            raise ValueError(f"Reference length must be positive on line {line_number}.")
        total += length
        if total > MAX_EXACT_INT:
            raise ValueError("Total reference length exceeds the supported exact integer range.")
        names.append(fields[0])
        ends.append(total)
    if not names:
        raise ValueError("The FAI file is empty.")
    return names, np.asarray(ends, dtype=np.int64)


def bounds_controls(axis, size):
    st.sidebar.caption(f"{axis} bounds: zero-based, end exclusive")
    if st.session_state.manual_bounds:
        start = st.sidebar.number_input(f"{axis} min", 0, size - 1, 0, 1)
        end = st.sidebar.number_input(f"{axis} max (exclusive)", 1, size, size, 1)
    else:
        start, end = st.sidebar.slider(
            f"{axis} range", 0, size, (0, size),
            help="The upper bound is excluded.",
        )
    if start >= end:
        st.sidebar.error(f"{axis} min must be less than {axis} max.")
        st.stop()
    return int(start), int(end)


def make_raster(x, y, values, ux, uy, x0, x1, y0, y1, total, cmap,
                background, unmapped_color, scale, smoothing):
    nx, ny = x1 - x0, y1 - y0
    width, height = min(nx, RASTER_SIZE), min(ny, RASTER_SIZE)
    # Each output bin represents the largest genomic coordinate within it.
    bins = np.zeros((height, width), dtype=np.int64)
    bx = (x - x0) * width // nx
    by = (y - y0) * height // ny
    np.maximum.at(bins, (by, bx), values)
    unmapped_bins = np.zeros((height, width), dtype=np.uint8)
    unmapped_bins[(uy - y0) * height // ny, (ux - x0) * width // nx] = 1
    if scale >= 1:
        bins = maximum_filter(bins, size=scale, mode="constant", cval=0)
        unmapped_bins = maximum_filter(unmapped_bins, size=scale, mode="constant", cval=0)

    colors = (colormaps[cmap](np.linspace(0, 1, 256)) * 255).astype(np.uint8)
    rgba = np.empty((height, width, 4), dtype=np.uint8)
    rgba[:] = background
    occupied = bins != 0
    indices = ((bins[occupied] - 1) / max(total - 1, 1) * 255).astype(np.uint8)
    rgba[occupied] = colors[indices]
    rgba[unmapped_bins != 0] = unmapped_color
    image = Image.fromarray(rgba, "RGBA")
    # Explicit interpolation makes the display mode independent of the browser.
    target_width, target_height = width * 2, height * 2
    resampling = Image.Resampling.NEAREST if smoothing == "nearest" else Image.Resampling.BILINEAR
    image = image.resize((target_width, target_height), resampling)
    buffer = io.BytesIO()
    image.save(buffer, format="PNG")
    return "data:image/png;base64," + base64.b64encode(buffer.getvalue()).decode("ascii"), target_width, target_height


st.set_page_config(page_title="Tileviewer", layout="wide")
st.sidebar.header("Tileviewer inputs")
matrix_file = st.sidebar.file_uploader(
    "Sparse CSR matrix (.npz) storing (x, y) coordinates with absolute genome positions for each cluster. "
    "Explicit stored zeros are unmapped clusters.", 
    type="npz"
)
fai_file = st.sidebar.file_uploader("Reference index (.fai)", type="fai")
if not matrix_file or not fai_file:
    st.info("Upload both a sparse matrix and its reference FAI index.")
    st.stop()

try:
    matrix = load_matrix(matrix_file.getvalue())
    names, ends = parse_fai(fai_file.getvalue())
    if matrix.nnz and int(matrix.data.max()) > int(ends[-1]):
        raise ValueError("A matrix coordinate exceeds the total reference length.")
except (ValueError, TypeError, OSError, UnicodeError, EOFError, KeyError, zipfile.BadZipFile) as exc:
    st.error(f"Could not load inputs: {exc}")
    st.stop()

st.sidebar.header("Crop")
st.sidebar.checkbox("Enter exact bounds", key="manual_bounds")
x0, x1 = bounds_controls("X (matrix rows)", matrix.shape[0])
y0, y1 = bounds_controls("Y (matrix columns)", matrix.shape[1])

st.sidebar.header("Display")
cmap = st.sidebar.selectbox("Colormap", ["viridis", "plasma", "turbo", "rainbow"])
background_name = st.sidebar.selectbox("Empty-cell background", ["Transparent", "White", "Dark grey", "Black"])
backgrounds = {
    "Transparent": (0, 0, 0, 0),
    "White": (255, 255, 255, 255),
    "Dark grey": (48, 48, 48, 255),
    "Black": (0, 0, 0, 255),
}
smoothing = st.sidebar.radio("Pixel rendering", ["nearest", "bilinear"], horizontal=True)
aspect = st.sidebar.radio("Aspect ratio", ["Equal", "Auto"], horizontal=True)
scale = st.sidebar.slider("Overview spot size (raster pixels)", 1, 10, 1, step=1)
dot_size = st.sidebar.slider("Cluster diameter (screen pixels)", 1, 20, 5)
unmapped_color = st.sidebar.color_picker("Unmapped cluster color", "#ff00ff")
image_rendering = "pixelated" if smoothing == "nearest" else "auto"
st.markdown(
    f"<style>.js-plotly-plot image {{ image-rendering: {image_rendering}; }}</style>",
    unsafe_allow_html=True,
)

# Slice before creating any dense raster; no dense matrix-sized array is made.
crop = matrix[x0:x1, y0:y1].tocoo(copy=False)
mapped = crop.data != 0
x = crop.row[mapped].astype(np.int64) + x0
y = crop.col[mapped].astype(np.int64) + y0
ux = crop.row[~mapped].astype(np.int64) + x0
uy = crop.col[~mapped].astype(np.int64) + y0
values = crop.data[mapped].astype(np.int64, copy=False)
area = (x1 - x0) * (y1 - y0)
count = len(values)
unmapped_count = len(ux)
cluster_count = count + unmapped_count
chromosome_indices = np.searchsorted(ends, values, side="left")
counts = np.bincount(chromosome_indices, minlength=len(names))

st.subheader("Active crop")
a, b = st.columns(2)
a.metric("Mapped clusters", f"{count:,}")
b.metric("Unmapped clusters", f"{unmapped_count:,}")

def select_cluster():
    points = st.session_state["tile_plot"].selection.points
    if points and points[-1].get("customdata"):
        st.session_state["pending_focus"] = points[-1]["customdata"]

dot_mode = cluster_count <= HOVER_LIMIT
pending_focus = st.session_state.pop("pending_focus", None)
if dot_mode and count and pending_focus:
    _, clicked_name, clicked_position = pending_focus
    if clicked_name in names:
        clicked_index = names.index(clicked_name)
        clicked_position = int(clicked_position)
        clicked_length = int(ends[clicked_index] - (ends[clicked_index - 1] if clicked_index else 0))
        if 1 <= clicked_position <= clicked_length:
            st.session_state["color_mode"] = "Chromosome + interval"
            st.session_state["focus_chromosome"] = clicked_name
            st.session_state[f"focus_min_{clicked_name}"] = max(1, clicked_position - FOCUS_RADIUS_BP)
            st.session_state[f"focus_max_{clicked_name}"] = min(clicked_length, clicked_position + FOCUS_RADIUS_BP)

color_mode = "Whole genome"
if dot_mode and count:
    color_mode = st.sidebar.radio(
        "Color mapped clusters by", ["Whole genome", "Chromosome + interval"], key="color_mode"
    )

starts = np.concatenate((np.array([0], dtype=np.int64), ends[:-1]))
relative_all = values - starts[chromosome_indices]
focus = None
if color_mode == "Chromosome + interval":
    selected_name = st.sidebar.selectbox(
        "Chromosome", names, index=int(np.argmax(counts)), key="focus_chromosome"
    )
    selected_index = names.index(selected_name)
    chr_length = int(ends[selected_index] - starts[selected_index])
    in_chr = chromosome_indices == selected_index
    chr_positions = relative_all[in_chr]
    default_min = int(chr_positions.min()) if chr_positions.size else 1
    default_max = int(chr_positions.max()) if chr_positions.size else chr_length
    focus_min = st.sidebar.number_input(
        "Interval start (1-based, inclusive)", 1, chr_length, default_min,
        key=f"focus_min_{selected_name}",
    )
    focus_max = st.sidebar.number_input(
        "Interval end (1-based, inclusive)", 1, chr_length, default_max,
        key=f"focus_max_{selected_name}",
    )
    if focus_min > focus_max:
        st.sidebar.error("Interval start must not exceed interval end.")
        st.stop()
    focus = in_chr & (relative_all >= focus_min) & (relative_all <= focus_max)

fig = go.Figure()
if not dot_mode:
    source, width, height = make_raster(
        x, y, values, ux, uy, x0, x1, y0, y1, int(ends[-1]), cmap,
        backgrounds[background_name],
        tuple(bytes.fromhex(unmapped_color[1:])) + (255,), scale, smoothing,
    )
    dx, dy = (x1 - x0) / width, (y1 - y0) / height
    fig.add_trace(go.Image(
        source=source, x0=x0 - 0.5 + dx / 2, y0=y0 - 0.5 + dy / 2,
        dx=dx, dy=dy, hoverinfo="skip",
    ))

def sample_indices(n, limit):
    if n <= limit:
        return np.arange(n)
    return np.linspace(0, n - 1, limit, dtype=np.int64)


if dot_mode:
    mapped_limit, unmapped_limit = count, unmapped_count
else:
    mapped_limit = round(HOVER_LIMIT * count / cluster_count)
    if count and unmapped_count:
        mapped_limit = min(HOVER_LIMIT - 1, max(1, mapped_limit))
    unmapped_limit = HOVER_LIMIT - mapped_limit

if unmapped_count:
    shown_unmapped = sample_indices(unmapped_count, unmapped_limit)
    fig.add_trace(go.Scattergl(
        x=ux[shown_unmapped], y=uy[shown_unmapped], mode="markers",
        marker={
            "size": dot_size*1.5 if dot_mode else scale,
            "symbol": "x",
            "color": unmapped_color if dot_mode else "rgba(0,0,0,0.01)",
        },
        name="Unmapped",
        hovertemplate="Grid: (%{x}, %{y})<br>Unmapped cluster<extra></extra>",
        showlegend=dot_mode,
    ))

if count:
    shown = sample_indices(count, mapped_limit)
    chrom = chromosome_indices[shown]
    relative = relative_all[shown]
    custom = [
        [str(int(value)), names[int(ci)], str(int(pos))]
        for value, ci, pos in zip(values[shown], chrom, relative)
    ]
    if dot_mode:
        colorscale = [
            [float(t), f"rgb({r},{g},{b})"]
            for t, (r, g, b, _) in zip(
                np.linspace(0, 1, 256),
                (colormaps[cmap](np.linspace(0, 1, 256)) * 255).astype(np.uint8),
            )
        ]
        if focus is None:
            groups = [(np.ones(len(shown), dtype=bool), {
                "size": dot_size*1.5, "color": values[shown], "colorscale": colorscale,
                "cmin": 1, "cmax": max(int(ends[-1]), 2), "showscale": True,
                "colorbar": {"title": "Absolute bp"},
            }, "Mapped")]
        else:
            highlighted = focus[shown]
            groups = [
                (~highlighted, {"size": dot_size*1.5, "color": "#888888"}, "Outside interval"),
                (highlighted, {
                    "size": dot_size*1.5, "color": relative[highlighted], "colorscale": colorscale,
                    "cmin": int(focus_min), "cmax": max(int(focus_max), int(focus_min) + 1),
                    "showscale": True, "colorbar": {"title": f"{selected_name} bp"},
                }, "In interval"),
            ]
    else:
        groups = [(np.ones(len(shown), dtype=bool),
                   {"size": 6, "color": "rgba(0,0,0,0.01)"}, "Mapped")]
    custom = np.asarray(custom, dtype=object)
    for mask, marker, label in groups:
        if not np.any(mask):
            continue
        fig.add_trace(go.Scattergl(
            x=x[shown][mask], y=y[shown][mask], mode="markers",
            marker=marker,
            customdata=custom[mask],
            hovertemplate=(
                "Grid: (%{x}, %{y})<br>Absolute: %{customdata[0]} bp"
                "<br>%{customdata[1]}:%{customdata[2]:,}<extra></extra>"
            ),
            name=label, showlegend=dot_mode and focus is not None,
        ))

fig.update_layout(
    height=700, margin=dict(l=20, r=20, t=80, b=20), hovermode="closest",
    uirevision=f"{x0}:{x1}:{y0}:{y1}:{aspect}",
    legend=dict(orientation="v", x=1.02, xanchor="left", y=1, yanchor="top"),
    plot_bgcolor=(
        "rgba(0,0,0,0)" if background_name == "Transparent"
        else {"White": "white", "Dark grey": "#303030", "Black": "black"}[background_name]
    ),
    dragmode="zoom",
    xaxis=dict(title="X (matrix row)", range=[x0 - 0.5, x1 - 0.5], constrain="domain"),
    yaxis=dict(
        title="Y (matrix column)", range=[y1 - 0.5, y0 - 0.5],
        scaleanchor="x" if aspect == "Equal" else False, scaleratio=1,
        constrain="domain",
    ),
)
for trace in fig.data:
    if trace.type == "scattergl" and trace.marker.showscale is True:
        trace.marker.colorbar.update(len=0.80, y=0, yanchor="bottom")
if dot_mode:
    st.plotly_chart(
        fig, use_container_width=True, config={"scrollZoom": True}, key="tile_plot",
        on_select=select_cluster, selection_mode="points",
    )
else:
    st.plotly_chart(fig, use_container_width=True, config={"scrollZoom": True})
st.caption(
    "Drag to zoom; use the toolbar to pan or reset. Statistics describe the sidebar crop, "
    "not the Plotly zoom. Cluster diameter is in screen pixels and does not change the data."
)
if dot_mode and count:
    st.caption("Click a mapped cluster to color its chromosome interval ±1 Mb. Unmapped clusters cannot be focused.")
if not dot_mode:
    st.warning(
        f"Showing a raster overview because this crop contains {cluster_count:,} clusters. "
        f"Crop to at most {HOVER_LIMIT:,} clusters for individual dots, interval coloring, "
        "and complete hover coverage. Unmapped clusters take priority in overlapping raster "
        "bins; otherwise each bin shows its largest genomic coordinate."
    )

st.subheader("Chromosome distribution")
present = counts > 0
if present.any():
    chart = go.Figure(go.Bar(
        x=np.asarray(names)[present], y=counts[present],
        hovertemplate="%{x}: %{y:,} cells<extra></extra>",
    ))
    chart.update_layout(
        height=300, margin=dict(l=20, r=20, t=10, b=20),
        xaxis_title="Reference sequence", yaxis_title="Mapped cells",
    )
    st.plotly_chart(chart, use_container_width=True)
else:
    st.info("No mapped cells in this crop.")
