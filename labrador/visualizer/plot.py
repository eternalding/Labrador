import plotly.express as px
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objs as go
import plotly.offline as py
from plotly.subplots import make_subplots
from ipywidgets import interactive, HBox, VBox, widgets, interact, fixed
from IPython.display import display

MAX_DISPLAY_PER_ROW = 3


def plot_interactive_matrix(arr: np.ndarray, x_labels: list, y_labels: list, *args, **kwargs):
    fig = px.imshow(arr, x=x_labels, y=y_labels,
                    color_continuous_scale='YlOrRd', aspect="equal", *args, **kwargs)

    if 'range_color' not in kwargs:
        c_min = arr.min()
        c_max = arr.max()
    else:
        c_min, c_max = kwargs['range_color']

    vb = set_ipython_color_slider(fig, c_min, c_max)
    return vb


def get_slider(c_min: float, c_max: float, description: str = "ColorBar"):
    slider = widgets.FloatRangeSlider(max=c_max,
                                      min=c_min,
                                      step=(c_max-c_min)/1000,
                                      readout=False,
                                      readout_format='.2f',
                                      description=description)
    return slider


def update_color_range(figure, range_color, coloraxis: str = "coloraxis"):
    display(range_color)
    figure.layout[coloraxis]["cmin"] = range_color[0]
    figure.layout[coloraxis]["cmax"] = range_color[1]


def set_ipython_color_slider(fig, c_min: float, c_max: float):
    f = go.FigureWidget(fig)
    slider = get_slider(c_min, c_max)
    vb = VBox([f, interactive(update_color_range,
                              figure=fixed(f),
                              range_color=slider,
                              coloraxis=fixed("coloraxis"))])
    vb.layout.align_items = 'center'
    return vb


def gen_heatmap(idx: int, arr: np.ndarray, start: int, end: int, name: str, resolution: int = 5000):
    HOVERTEMPLATE = "<b>Bin1:</b> %{x}<br>" \
                    "<b>Bin2:</b> %{y}<br>" \
                    "<b>Value</b>: %{z}<br>"
    heatmap = go.Heatmap(x=np.arange(start, end, resolution),
                         y=np.arange(start, end, resolution),
                         z=arr,
                         coloraxis=f"coloraxis{idx}",
                         name=name,
                         colorscale="YlOrRd",
                         hovertemplate=HOVERTEMPLATE,
                         zauto=True,
                         showscale=False,
                         )
    return heatmap


def plot_multiple_arrays(arrs: dict, title: str = None):
    max_rows = ((len(arrs)-1) // MAX_DISPLAY_PER_ROW) + 1

    # All arrs were plotted by this order
    fields = list(arrs.keys())

    # Define figure
    fig = make_subplots(rows=max_rows,
                        cols=MAX_DISPLAY_PER_ROW,
                        subplot_titles=fields)

    fig.update_coloraxes(colorscale='YlOrRd', showscale=False)

    for idx, field in enumerate(fields, 1):
        row = ((idx-1) // MAX_DISPLAY_PER_ROW) + 1
        col = ((idx-1) % MAX_DISPLAY_PER_ROW) + 1
        fig.add_trace(arrs[field], row, col)
        fig.layout[f'coloraxis{idx}'] = fig.layout['coloraxis']
        fig.layout[f'coloraxis{idx}'].cmin = arrs[field]['z'].min()
        fig.layout[f'coloraxis{idx}'].cmax = arrs[field]['z'].max()

    for key in fig.layout:
        if 'xaxis' in key or 'yaxis' in key:
            fig.layout[key]['tickangle'] = -45
            if 'yaxis' in key:
                # Reverse y-axis to make origin on upper-left
                fig.layout[key]['autorange'] = 'reversed'

    fig.update_layout(height=max_rows*300 + 100,
                      width=MAX_DISPLAY_PER_ROW*300,
                      title_text=title)

    # Create multiple slide-bars
    f = go.FigureWidget(fig)

    sliders = []
    for idx, field in enumerate(fields, 1):
        c_min = arrs[field]['z'].min()
        c_max = arrs[field]['z'].max()
        slider = get_slider(c_min, c_max, description=field)
        widget = interactive(update_color_range,
                             figure=fixed(f),
                             range_color=slider,
                             coloraxis=fixed(f"coloraxis{idx}"))
        sliders.append(widget)

    slider_vb = VBox(sliders)
    total_vb = VBox([slider_vb, f])
    total_vb.layout.align_items = 'center'
    return total_vb
