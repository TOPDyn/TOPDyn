import numpy as np
import plot_grid as mf
import pyqtgraph as pg
import pyqtgraph.exporters
import matplotlib.pyplot as plt
from PyQt5 import QtCore, QtGui
from functions_2d import generate_xy_coord

def legend_constr(constr_func):
    """ Defines the labels  of the constraint functions.

    Args:
        constr_func (:obj:`list`): Restriction functions applied.

    Returns:
        Numpy array with the labels.
    """
    label = np.empty(len(constr_func), dtype=object)
    func, index = np.unique(constr_func, return_index=True)
    aux = np.arange(0, len(constr_func))
    aux2 = np.setdiff1d(aux, index)
    label[aux2] = None
    i = 0
    for f in func:
        if f =="area":
            label[index[i]] = "constraint - area"
        elif f == "r_ratio":
            label[index[i]] = "constraint - r ratio"
        elif f == "compliance":
            label[index[i]] = "constraint - compliance"
        elif f == "local_ep":
            label[index[i]] = "constraint - local ep"
        elif f == "local_ki":
            label[index[i]] = "constraint - local ki"
        elif f == "local_r":
            label[index[i]] = "constraint - local r"
        i += 1
    return label

def set_pen(f):
    colors = [(43,174,179), (64,66,114), (255,110,60), (255,215,75), (255,102,0), (255,128,128), (0,51,0)]
    
    if f == "area":
        pen_set = {'color': colors[0], 'width': 2}
    elif f == "r_ratio":
        pen_set = {'color': colors[1], 'width': 2}
    elif f == "compliance":
        pen_set = {'color': colors[2], 'width': 2}
    elif f == "local_ep":
        pen_set = {'color': colors[3], 'width': 2}
    elif f == "local_ki":
        pen_set = {'color': colors[4], 'width': 2}
    elif f == "local_r":
        pen_set = {'color': colors[5], 'width': 2}
    return pen_set

def window_each_iter(constr_func, func_name, label):
    """ Generates a window to plot the optimized mesh and the convergence graph in the same window.

    Args:
        constr_func (:obj:`list`): Restriction functions applied.
        func_name (:obj:`str`): Objective function name.
        label (:obj:`numpy.array`): Label of the constraint functions.

    Returns:
        Principal window, convergece graph, optimized part.
    """
    win = pg.GraphicsLayoutWidget(show=True)
    win.resize(900,600)
    win.setWindowTitle('MMA')
    
    grid = mf.PColorMeshItem(cmap='grey')
    plot = win.addPlot()
    plot.setAspectLocked(True)
    plot.hideAxis('bottom')
    plot.hideAxis('left')
    plot.addItem(grid)
    # Plot Objective function and area
    win.nextRow()
    p2 = win.addPlot()
    p2.addLegend(labelTextColor=(0,0,0), offset=(700,10))
    p2.setLabel('left', func_name.lower())
    p2.setLabel('bottom', "iteration") 
    
    curves_p2 = []
    curves_p2.append(p2.plot(pen={'color': (0,0,0), 'width': 2}, name=func_name.lower()))
    for ind, f in enumerate(constr_func):
        pen_set = set_pen(f)
        curves_p2.append(p2.plot(name=label[ind], pen=pen_set))
    return win, p2, curves_p2, grid

def simple_window():
    """ Generates a window to plot the optimized mesh.

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

def win_convergence(constr_func, list_iter, list_f0val, list_fvals, func_name, label):
    """ Generates a window to plot the convergence graph.

    Args:
        constr_func (:obj:`list`): Restriction functions applied.
        list_iter (:obj:`list`): All iteration values.
        list_f0val (:obj:`list`): All objective function values.
        list_fvals (:obj:`list`): All constraint function values.
        func_name (:obj:`str`): Objective function name.
        label (:obj:`numpy.array`): Label of the constraint functions.

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
   
    for ind, f in enumerate(constr_func):
        pen_set = set_pen(f)
        p.plot(list_iter, list_fvals[ind], pen=pen_set, name=label[ind])
    return win, p

def set_coord_grid(lx, ly, nelx, nely):
    """ Defines dimensions of the optimized part.   """
    x, y   = generate_xy_coord(lx, ly, nelx, nely)
    x_plot = np.repeat(x, nely+1).reshape(nelx+1, nely+1)
    y_plot = np.tile(y, nelx+1).reshape(nelx+1, nely+1)
    return x_plot, y_plot

def set_grid_data(grid, xval, x_plot, y_plot, nelx, nely):
    grid.setData(x_plot, y_plot, xval.reshape(nelx, nely, order='F'))

def set_conv_data(outeriter, curves_funcs, list_iter, list_f0val, list_fvals, constr_func):
    """ Updates values of the objective function and the constraint function to plot the convergence graph.   """
    curves_funcs[0].setData(list_iter[:outeriter+1], list_f0val[:outeriter+1])

    for ind in range(len(constr_func)):
        curves_funcs[ind+1].setData(list_iter[:outeriter+1], list_fvals[:outeriter+1, ind])
   
def update_conv(constr_func, p, list_iter, list_f0val, list_fvals):
    """ Updates values of the objective function and the constraint function to plot the convergence graph.

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
    for ind, f in enumerate(constr_func):
        pen_set = set_pen(f)
        p.plot(list_iter, list_fvals[ind], pen=pen_set)
    return p

def save_fig(fig, path, pg_graph):
    """ Saves the graphics: optimized part, convergence graph, frequency response graph and deformed mesh.

    Args:
        fig: Object with the graph. It can be: - pyqtgraph.graphicsItems.PlotItem or - matplotlib.pyplot.figure.
        path: Directory to save the graph.
        pg_graph (:obj:`bool`): Specifies the type of figure.
    """   
    if pg_graph:
        exporter = pg.exporters.ImageExporter(fig)
        exporter.export(path)
    else:
        fig.savefig(path)

def compare_freqresponse(freq_range, newf, oldf, func_name):
    """ Plot the frequency response of the original and the optimized function.

    Args:
        freq_range (:obj:`list`): Range of analyzed frequencies.
            First value is the minimum frequency.
            Second value is the maximum frequency.
            Third value is the step between each calculation of the objective function. 
        newf (array): Optimized function.
        oldf (array): Original function.
        func_name (:obj:`str`): Objective function name.
            It can be: "compliance", "input_power", "elastic_potential_energy", "kinetic_energy" or "r_ratio".

    Returns:
        A figure object and a single Axes object from matplotlib.pyplot.
    """
    freq = np.arange(freq_range[0], freq_range[1] + 1, freq_range[2])
    fig, ax = plt.subplots()
    ax.plot(freq, oldf, label = 'original')
    ax.plot(freq, newf, label = 'optimized')
    ax.set_xlabel("frequency [Hz]")
    ax.set_ylabel(func_name.lower())
    ax.legend(loc='best', frameon=False)
    ax.set_yscale('log')
    return fig, ax

def freqrsp_modes(freq_range, newf, oldf, modes, func_name, save):
    """ Plot the frequency response of the function with multiple modes.

    Args:
        freq_range (:obj:`list`): Range of frequencies analyzed.
            First value is the minimum frequency.
            Second value is the maximum frequency.
            Third value is the step between each calculation of the objective function. 
        newf (:obj:`numpy.array`): Optimized function.
        oldf (:obj:`numpy.array`): Original function.
        modes (:obj:`list`): Analyzed modes. 
        func_name (:obj:`str`): Objective function name.
            It can be: "compliance", "input_power", "elastic_potential_energy", "kinetic_energy" or "r_ratio".
        save (:obj:`bool`): True for save the graphic in PNG.
    """
    freq = np.arange(freq_range[0], freq_range[1] + 1, freq_range[2])
    
    for i, mode in enumerate(modes):
        plt.figure(i+1)
        plt.plot(freq, oldf.real, label='original')
        plt.plot(freq, newf[:, i].real, label=str(mode) + ' mode')
        plt.title(str(mode) + ' mode')
        plt.xlabel('frequency [Hz]', fontsize=16)
        plt.ylabel(func_name.lower(), fontsize=16)
        plt.yscale('log')
        plt.legend()
        if save:
            plt.savefig(str(mode) + ".eps")

    plt.figure(i+1)
    #plt.plot(freq, oldf.real, label='original')
    for i, mode in enumerate(modes):
        plt.plot(freq, newf[:, i].real, label=str(mode) + ' mode')
        plt.xlabel('frequency [Hz]', fontsize=16)
        plt.ylabel(func_name.lower(), fontsize=16)
        plt.yscale('log')
    plt.title('All modes')
    plt.legend()
    if save:
        plt.savefig('all' + ".eps")   
    plt.show()

def compare_deriv(nodes, delta_d, dw, dw_orig):
    for ind2, node in enumerate(nodes):
        plt.figure(ind2+1)
        plt.plot(delta_d, dw[:, ind2], marker='o', label='FDM')
        plt.plot(delta_d, np.repeat(dw_orig[ind2], len(delta_d)), marker='o', label='Analytical Method')
        plt.title('Node ' + str(node), fontsize=18)
        plt.xlabel(r'$\Delta D_i$', fontsize=18)
        plt.ylabel(r'$\frac{d \alpha}{d D_i}$', fontsize=18) #, rotation=0
        plt.legend()
    plt.show(block=True)