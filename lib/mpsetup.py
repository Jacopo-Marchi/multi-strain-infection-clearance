import matplotlib.pyplot as plt
import matplotlib as mpl
import string
import itertools


def despine(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')


def label_axes(fig_or_axes, labels,
               labelstyle=r'{\sf \textbf{%s}}',
               xy=(-0.05, 0.95), xycoords='axes fraction', **kwargs):
    """
    Walks through axes and labels each.
    kwargs are collected and passed to `annotate`

    Parameters
    ----------
    fig : Figure or Axes to work on
    labels : iterable or None
        iterable of strings to use to label the axes.
        If None, lower case letters are used.

    loc : Where to put the label units (len=2 tuple of floats)
    xycoords : loc relative to axes, figure, etc.
    kwargs : to be passed to annotate
    """
    # re-use labels rather than stop labeling
    labels = itertools.cycle(labels.upper())
    axes = fig_or_axes.axes if isinstance(fig_or_axes, plt.Figure) else fig_or_axes
    for ax, label in zip(axes, labels):
        ax.annotate(labelstyle % label, xy=xy, xycoords=xycoords,
                    **kwargs)


class OffsetHandlerTuple(mpl.legend_handler.HandlerTuple):
    """
    Legend Handler for tuple plotting markers on top of each other
    """
    def __init__(self, **kwargs):
        mpl.legend_handler.HandlerTuple.__init__(self, **kwargs)

    def create_artists(self, legend, orig_handle,
                       xdescent, ydescent, width, height, fontsize,
                       trans):
        nhandles = len(orig_handle)
        perside = (nhandles - 1) / 2
        offset = height / nhandles
        handler_map = legend.get_legend_handler_map()
        a_list = []
        for i, handle1 in enumerate(orig_handle):
            handler = legend.get_legend_handler(handler_map, handle1)
            _a_list = handler.create_artists(legend, handle1,
                                             xdescent,
                                             offset*i+ydescent-offset*perside,
                                             width, height,
                                             fontsize,
                                             trans)
            a_list.extend(_a_list)
        return a_list



def get_cmap(N):
    ''' Returns a function that maps each index in 0, 1, ...
        N-1 to a distinct RGB color.'''
    color_norm  = colors.Normalize(vmin=0, vmax=N-1)
    scalar_map = cm.ScalarMappable(norm=color_norm, cmap='jet')
    def map_index_to_rgb_color(index):
        return scalar_map.to_rgba(index)
    return map_index_to_rgb_color
    
    
def gridint(x, y, z, resX=1000, resY=1000):
    "Convert 3 column data to matplotlib grid"
    xi = np.geomspace(min(x), max(x), resX)
    yi = np.geomspace(min(y), max(y-0.01), resY)
    
    print((min(xi), min(yi)))
    print((max(xi), max(yi), max(y)))
    Z = interpn((x, y), z, (xi[None,:], yi[:,None]), method='nearest')
    X, Y = np.meshgrid(xi, yi)
    return X, Y, Z
    
