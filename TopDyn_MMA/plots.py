import numpy as np
import pyqtgraph as pg
import pyqtgraph.exporters
import matplotlib.pyplot as plt

def freqresponse(freq_range, delta, obj_func, func_name, save=False):
    """ Plot the frequency response.

    Args:
        freq_range (list): Range of frequencies analyzed.
        delta (int): Steps between 1 and the parameter freq_range.
        obj_func (list): Objective function.
        
    """
    freq = np.arange(freq_range[0], freq_range[1] + 1, delta)
    plt.plot(freq, obj_func.real)
    plt.xlabel('frequency [Hz]', fontsize=16)
    plt.ylabel(func_name.lower(), fontsize=16)
    plt.yscale('log')
    if save:
        plt.savefig("freq_response1.png")
    plt.show()

def compare_freqresponse(freq_range, delta, newf, oldf, func_name, save):
    """ Plot the frequency response of the original and the optimized function.

    Args:
        freq_range (list): Range of frequencies analyzed.
        delta (int): Steps between 1 and the parameter freq_range.
        newf (array): Optimized function.
        oldf (array): Original function.
        func_name (str): Objective function used.
            It can be: 'Compliance', 'Input Power', 'Elastic Potential Energy', 'Kinetic Energy' or 'R Ratio'.
        save (bool, optional): True for save the graphic in PNG.

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

def window_each_iter(constr_func, list_iter, list_f0val, list_fvals, xval, nelx, nely, func_name):

    win = pg.GraphicsLayoutWidget(show=True)
    win.resize(900,600)
    win.setWindowTitle('MMA')
    # Enable antialiasing for prettier plots
    pg.setConfigOptions(antialias=True)
    vb = pg.ViewBox()
    p1 = win.addItem(vb)
    # Configure view for images
    vb.setAspectLocked()
    vb.invertY()
    vb.show()
    # Create image
    image = pg.ImageItem(-xval.real.reshape((nely, nelx)), axisOrder='row-major')
    # Display image
    vb.addItem(image)
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

def simple_window(xval, nelx, nely):
    gv = pg.GraphicsView()
    gv.setWindowTitle('MMA')
    gv.resize(1000,600)
    vb = pg.ViewBox()
    gv.setCentralItem(vb)
    gv.show()
    # Configure view for images
    vb.setAspectLocked()
    vb.invertY()
    vb.show()
    # Create image
    image = pg.ImageItem(-xval.real.reshape((nely, nelx)), axisOrder='row-major')
    # Display image
    vb.addItem(image)

    return gv, image

def win_convergence(constr_func, list_iter, list_f0val, list_fvals, func_name):

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
    p.plot(list_iter, list_f0val, pen={'color': (0,0,0), 'width': 2})
    colors = [(43,174,179), (64,66,114), (255,110,60), (255,215,75)]
    if 'Area' in constr_func:
        ind = constr_func.index("Area")
        p.plot(list_iter, list_fvals[ind], pen={'color': colors[0], 'width': 2})
        
    if 'R Ratio' in constr_func:
        ind = constr_func.index("R Ratio")
        p.plot(list_iter, list_fvals[ind], pen={'color': colors[1], 'width': 2})

    return p

def save_fig(imagewindow, convergence):
    
    exporter = pg.exporters.ImageExporter(imagewindow)
    exporter.export('optimization.png')

    exporter2 = pg.exporters.ImageExporter(convergence)
    exporter2.export('convergence.png')