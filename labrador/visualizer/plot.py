import plotly.express as px
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objs as go
import plotly.offline as py
from ipywidgets import interactive, HBox, VBox, widgets, interact


def plot_interactive_matrix(arr: np.ndarray, x_labels: list, y_labels: list, resolution: int, *args, **kwargs):
    fig = px.imshow(arr, x=x_labels, y=y_labels,
                    color_continuous_scale='YlOrRd', aspect="equal", *args, **kwargs)

    if 'range_color' not in kwargs:
        c_min = arr.min()
        c_max = arr.max()
    else:
        c_min, c_max = kwargs['range_color']

    vb = set_ipython_color_slider(fig, c_min, c_max)
    return vb


def set_ipython_color_slider(fig, c_min: float, c_max: float):
    f = go.FigureWidget(fig)
    slider = widgets.FloatRangeSlider(max=c_max,
                                      min=c_min,
                                      step=0.1,
                                      readout=False,
                                      readout_format='.2f',
                                      description='ColorBar')

    def update_range(range_color):
        print(range_color)
        f.layout.coloraxis.cmin = range_color[0]
        f.layout.coloraxis.cmax = range_color[1]

    vb = VBox((f, interactive(update_range, range_color=slider)))
    vb.layout.align_items = 'center'
    return vb
