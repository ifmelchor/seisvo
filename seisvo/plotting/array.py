
import matplotlib.pyplot as plt
import matplotlib.colors as mcolor
import matplotlib.ticker as mtick
import matplotlib.dates as mdates


def location_map(arr, exclude_locs=[], show=True):

    fig, ax = plt.subplots(1,1)

    for loc, utm in arr.utm.items():
        if loc not in exclude_locs:
            x, y = utm["easting"], utm["northing"]
            xdisp = x+0.01*x
            ydisp = y+0.01*y
            ax.scatter(x, y, marker="^", color="k")
            ax.annotate(loc, xy=(xdisp, ydisp))
    
    ax.set_ylabel("UTM northing")
    ax.set_xlabel("UTM easting")

    if show:
        plt.show()

    return fig

    