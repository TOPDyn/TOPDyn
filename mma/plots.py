import numpy as np
import pyqtgraph as pg
import pyqtgraph.exporters
from PyQt5 import QtCore
import matplotlib.pyplot as plt

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
    #plt.show()

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

def window_each_iter(constr_func, list_iter, list_f0val, list_fvals, xval, nelx, nely, lx, ly, func_name):
    """ Generate a window to plot the optimized mesh and the convergence graph in the same window.

    Args:
        constr_func (:obj:`list`): Restriction functions applied.
        list_iter (:obj:`list`): All iteration values.
        list_f0val (:obj:`list`): All objective function values.
        list_fvals (:obj:`list`): All constraint function values.
        xval (:obj:`numpy.array`): Indicates where there is mass.
        nelx (:obj:`int`): Number of elements on the X-axis.
        nely (:obj:`int`): Number of elements on the Y-axis.
        func_name (:obj:`str`): Objective function name.

    Returns:
        Principal window, convergece graph, optimezed part.
    """
    win = pg.GraphicsLayoutWidget(show=True)
    win.resize(900,600)
    win.setWindowTitle('MMA')
    # Enable antialiasing for prettier plots
    pg.setConfigOptions(antialias=True)
    p1 = win.addPlot()
    p1.hideAxis('bottom')
    p1.hideAxis('left')
    win.setAspectLocked(False)
    p1.setAspectLocked(False)
    image = pg.ImageItem()
    image.setImage(-xval.real.reshape((nelx, nely), order='F'), levels=[-1,0], autoRange=False)
    rect = QtCore.QRectF(0, 0, lx, ly)
    image.setRect(rect)
    p1.addItem(image)
    win.nextRow()
    # Plot Objective function and area
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

    return win, p2, image

def simple_window(xval, nelx, nely, lx, ly):
    """ Generate a window to plot the optimized mesh.

    Args:
        xval (:obj:`numpy.array`): Indicates where there is mass.
        nelx (:obj:`int`): Number of elements on the X-axis.
        nely (:obj:`int`): Number of elements on the Y-axis.

    Returns:
        Principal window, optimezed part.
    """
    win = pg.GraphicsLayoutWidget(show=True)
    win.resize(900,600)
    win.setWindowTitle('MMA')
    # Enable antialiasing for prettier plots
    pg.setConfigOptions(antialias=True)
    p1 = win.addPlot()
    p1.hideAxis('bottom')
    p1.hideAxis('left')
    win.setAspectLocked(False)
    p1.setAspectLocked(False)
    image = pg.ImageItem()
    image.setImage(-xval.real.reshape((nelx, nely), order='F'), levels=[-1,0], autoRange=False)
    rect = QtCore.QRectF(0, 0, lx, ly)
    image.setRect(rect)
    p1.addItem(image)

    return win, image

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