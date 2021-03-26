import numpy as np
import plot_grid as mf
import pyqtgraph as pg
import pyqtgraph.exporters
import matplotlib.pyplot as plt
from PyQt5 import QtCore, QtGui
from functions_2d import generate_xy_coord

def freqresponse(freq_range, delta, obj_func, func_name, save=None):
    """ Plot the frequency response.

    Args:
        freq_range (:obj:`list`): Range of frequencies analyzed.
            First value is the minimum frequency.
            Second value is the maximum frequency.
        delta (:obj:`int`): Step between each calculation of the function. 
        obj_func (:obj:`list`): Objective function values.
        func_name (:obj:`str`): Objective function name.
            It can be: 'Compliance', 'Input Power', 'Elastic Potential Energy', 'Kinetic Energy' or 'R Ratio'.
        save (:obj:`bool`, optional): True for save the graphic in PNG. Defaults to False.
    """
    fig, ax = plt.subplots()
    freq = np.arange(freq_range[0], freq_range[1] + 1, delta)
    ax.plot(freq, obj_func.real)
    ax.set_xlabel('frequency [Hz]', fontsize=16)
    ax.set_ylabel(func_name.lower(), fontsize=16)
    ax.set_yscale('log')
    if save is not None:
        plt.savefig(save + ".eps")
    plt.show()

def compare_freqresponse(freq_range, delta, newf, oldf, func_name, save):
    """ Plot the frequency response of the original and the optimized function.

    Args:
        freq_range (:obj:`list`): Range of frequencies analyzed.
            First value is the minimum frequency.
            Second value is the maximum frequency.
        delta (:obj:`int`): Step between each calculation of the function. 
        newf (array): Optimized function.
        oldf (array): Original function.
        func_name (:obj:`str`): Objective function name.
            It can be: 'Compliance', 'Input Power', 'Elastic Potential Energy', 'Kinetic Energy' or 'R Ratio'.
        save (:obj:`bool`, optional): True for save the graphic in PNG. Defaults to False.

    Returns:
        A single Axes object from matplotlib.pyplot.
    """
    freq = np.arange(freq_range[0], freq_range[1] + 1, delta)
    fig, ax = plt.subplots()
    ax.plot(freq, oldf, label = 'original')
    ax.plot(freq, newf, label = 'optimized')
    ax.set_xlabel("frequency [Hz]")
    ax.set_ylabel(func_name.lower())
    ax.legend(loc='best', frameon=False)
    ax.set_yscale('log')
    if save:
        plt.savefig("compare_freq_response.png")

    return ax

def window_each_iter(constr_func, list_iter, list_f0val, list_fvals, func_name, xval, lx, ly, nelx, nely):
    """ Generate a window to plot the optimized mesh and the convergence graph in the same window.

    Args:
        constr_func (:obj:`list`): Restriction functions applied.
        list_iter (:obj:`list`): All iteration values.
        list_f0val (:obj:`list`): All objective function values.
        list_fvals (:obj:`list`): All constraint function values.
        func_name (:obj:`str`): Objective function name.

    Returns:
        Principal window, convergece graph, optimezed part.
    """
    win = pg.GraphicsLayoutWidget(show=True)
    win.resize(900,600)
    win.setWindowTitle('MMA')
    #
    grid = mf.PColorMeshItem(cmap='grey')
    x_plot, y_plot = set_coord_grid(lx, ly, nelx, nely)
    set_grid_data(grid, xval, x_plot, y_plot, nelx, nely)
    plot = win.addPlot()
    plot.setAspectLocked(True)
    plot.hideAxis('bottom')
    plot.hideAxis('left')
    plot.addItem(grid)
    # Plot Objective function and area
    win.nextRow()
    p2 = win.addPlot()
    p2.addLegend(labelTextColor=(0,0,0), offset=(700,10))
    p2.plot(list_iter, list_f0val, pen={'color': (0,0,0), 'width': 2}, name=func_name.lower())
    p2.setLabel('left', func_name.lower())
    p2.setLabel('bottom', "iteration") 

    colors = [(43,174,179), (64,66,114), (255,110,60), (255,215,75)]
    if 'Area' in constr_func:
        ind = constr_func.index("Area")
        p2.plot(list_iter, list_fvals[ind], pen={'color': colors[0], 'width': 2}, name='constraint - area')
        
    if 'R Ratio' in constr_func:
        ind = constr_func.index("R Ratio")
        p2.plot(list_iter, list_fvals[ind], pen={'color': colors[1], 'width': 2}, name='constraint - r ratio')
    pg.QtGui.QApplication.processEvents()
    return win, p2, grid

def simple_window():
    """ Generate a window to plot the optimized mesh.

    Returns:
        Principal window, optimized part.
    """
    win = pg.GraphicsLayoutWidget(show=True)
    win.setWindowTitle('MMA')

    grid = mf.PColorMeshItem(cmap='grey')
    plot = win.addPlot()
    plot.setAspectLocked(True)
    plot.hideAxis('bottom')
    plot.hideAxis('left')
    plot.addItem(grid)
    pg.QtGui.QApplication.processEvents()

    return win, grid

def win_convergence(constr_func, list_iter, list_f0val, list_fvals, func_name):
    """ Generate a window to plot the convergence graph.

    Args:
        constr_func (:obj:`list`): Restriction functions applied.
        list_iter (:obj:`list`): All iteration values.
        list_f0val (:obj:`list`): All objective function values.
        list_fvals (:obj:`list`): All constraint function values.
        func_name (:obj:`str`): Objective function name.

    Returns:
        Principal window, convergece graph.
    """
    win = pg.GraphicsLayoutWidget(show=True, title="MMA")
    win.resize(1000,600)
    #p = win.addPlot(title="Convergence")
    p = win.addPlot()
    p.addLegend(labelTextColor=(0,0,0), offset=(800,10))
    p.plot(list_iter, list_f0val, pen={'color': (0,0,0), 'width': 2}, name=func_name.lower())
    p.setLabel('left', func_name.lower())
    p.setLabel('bottom', "iteration") 
    colors = [(43,174,179), (64,66,114), (255,110,60), (255,215,75)]
    if 'Area' in constr_func:
        ind = constr_func.index("Area")
        p.plot(list_iter, list_fvals[ind], pen={'color': colors[0], 'width': 2}, name='constraint - area')
        
    if 'R Ratio' in constr_func:
        ind = constr_func.index("R Ratio")
        p.plot(list_iter, list_fvals[ind], pen={'color': colors[1], 'width': 2}, name='constraint - r ratio')

    return win, p

def set_coord_grid(lx, ly, nelx, nely):
    x, y   = generate_xy_coord(lx, ly, nelx, nely)
    x_plot = np.repeat(x, nely+1).reshape(nelx+1, nely+1)
    y_plot = np.tile(y, nelx+1).reshape(nelx+1, nely+1)
    return x_plot, y_plot

def set_grid_data(grid, xval, x_plot, y_plot, nelx, nely):
    grid.setData(x_plot, y_plot, xval.reshape(nelx, nely, order='F'))
    
def convergence(constr_func, p, list_iter, list_f0val, list_fvals):
    """ Update the values of the objective function and the constraint function to plot the convergence graph.

    Args:
        constr_func (:obj:`list`): Restriction functions applied.
        p (pyqtgraph.graphicsItems.PlotItem): Convergence graph window
        list_iter (:obj:`list`): All iteration values.
        list_f0val (:obj:`list`): All objective function values.
        list_fvals (:obj:`list`): All constraint function values.

    Returns:
        Convergence graph window.
    """
    p.plot(list_iter, list_f0val, pen={'color': (0,0,0), 'width': 2})
    colors = [(43,174,179), (64,66,114), (255,110,60), (255,215,75)]
    if 'Area' in constr_func:
        ind = constr_func.index("Area")
        p.plot(list_iter, list_fvals[ind], pen={'color': colors[0], 'width': 2})
        
    if 'R Ratio' in constr_func:
        ind = constr_func.index("R Ratio")
        p.plot(list_iter, list_fvals[ind], pen={'color': colors[1], 'width': 2})

    return p

def save_fig(opt_part, convergence):
    """ Save the optimized part and the convergence graph.

    Args:
        opt_part (pyqtgraph.graphicsItems.PlotItem): Optimized part window. 
        convergence (pyqtgraph.graphicsItems.PlotItem): Convergence graph window.
    """   
    exporter = pg.exporters.ImageExporter(opt_part)
    exporter.export('optimization.png')

    exporter2 = pg.exporters.ImageExporter(convergence)
    exporter2.export('convergence.png')