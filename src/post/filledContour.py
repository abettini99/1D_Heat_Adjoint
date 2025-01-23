import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
from typing import Any
plt.style.use('src/post/plotting.style')

## @brief Creates a contour plot.
#
#  @param x      2D numpy array of data x-positions
#  @param y      2D numpy array of data y-positions
#  @param f      2D numpy array of data values
#  @param levels the number of visible contours, default 10
#  @param title  the title of the plot, default none
#  @param xlabel the x-label of the plot, default none
#  @param ylabel the y-label of the plot, default none
#
#  @return fig   container of all plot elements
#  @return ax    (list of) container(s) of all elements of individual (sub-)plot in the \p fig.
def plot(x : npt.NDArray[np.float64], 
         y : npt.NDArray[np.float64], 
         f : npt.NDArray[np.float64], 
         levels : int = 10,
         title  : str = None,
         xlabel : str = None,
         ylabel : str = None
         ) -> tuple[plt.Figure, Any]:

    # Create figure and axes
    fig, ax = plt.subplots(1,1, squeeze=False, figsize=(6.4,4.8))

    # Plot filled contour
    axc = ax[0,0].contourf(x,y,f, levels=levels)

    # Set limits
    ax[0,0].set_xlim(x[0,:][0], x[0,:][-1])
    ax[0,0].set_ylim(y[:,0][0], y[:,0][-1])

    # Set title
    if title:
        ax[0,0].set_title(title)
    # Set xlabel
    if xlabel:
        ax[0,0].set_xlabel(xlabel)
    # Set ylabel
    if ylabel:
        ax[0,0].set_ylabel(ylabel)

    # Set grid and ticks
    ax[0,0].grid(True,which="major",color="#999999",alpha=0.50)
    ax[0,0].grid(True,which="minor",color="#DDDDDD",ls="--",alpha=0.25)
    ax[0,0].minorticks_on()
    ax[0,0].tick_params(which='major', length=3.5, width=2, direction='out',
                        bottom=True, top=True, left=True, right=True)
    ax[0,0].tick_params(which='minor', length=1.75, width=2, direction='out',
                        bottom=True, top=True, left=True, right=True)
    ax[0,0].set_aspect('equal')
    fig.colorbar(axc)
    
    return fig, ax