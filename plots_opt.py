import os
import numpy as np
import pyqtgraph as pg
import pyqtgraph.exporters # to save

class PlotOpt():
    def __init__(self, lx, ly, nelx, nely, constr_func, constr_values, max_iter):
        self.nelx = nelx
        self.nely = nely

        self.lx = lx
        self.ly = ly

        self.directory = os.path.join(os.path.dirname(__file__), 'temp')

        self.constr_func = constr_func
        if constr_values.size == 1:
            self.constr_values = np.array([constr_values])
        else:
            self.constr_values = constr_values
        self.labels_constr = self._legend_constr()

        self.x_plot, self.y_plot = self._set_coord_grid(self.lx, self.ly, self.nelx, self.nely)

        self.list_iter  = np.empty(max_iter + 1)
        self.list_f0val = np.empty(max_iter + 1)
        self.list_fvals = np.empty((max_iter + 1, len(self.constr_func)))

    def update_lists(self, outit, fval, f0val):
        """ Adds new values to list of functions to plot convergence.

        Args:
            outit (:obj:`int`): Iteration.
            fval (:obj:`numpy.array`): Constraint function.
            f0val (:obj:`numpy.array`): Objective function.

        Returns:
            A tuple of lists with iterations, objective and constraint function.
        """
        if fval.size == 1:
            fval = np.array([fval])
        self.list_iter[outit] = outit+1
        self.list_f0val[outit] = f0val

        for i in range(len(self.constr_values)):
            if self.constr_values[i] > 0:
                self.list_fvals[outit, i] = (fval[i] + self.constr_values[i])
            else:
                self.list_fvals[outit, i] = (fval[i] - self.constr_values[i])

    def _legend_constr(self):
        """ Defines the labels  of the constraint functions.

        Args:
            constr_func (:obj:`list`): Restriction functions applied.

        Returns:
            Numpy array with the labels.
        """
        label = np.empty(len(self.constr_func), dtype=object)
        func, index = np.unique(self.constr_func, return_index=True)
        aux = np.arange(0, len(self.constr_func))
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

    def set_pen(self, f):
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

    def _set_coord_grid(self, lx, ly, nelx, nely):
        """ Defines dimensions of the optimized part. """
        x = lx/nelx * np.arange(nelx + 1) 
        y = ly/nely * np.arange(nely + 1)
        x_plot = np.repeat(x, nely + 1).reshape(nelx + 1, nely + 1)
        y_plot = np.tile(y, nelx + 1).reshape(nelx + 1, nely + 1)
        return x_plot, y_plot

    def set_grid_data(self, xval):
        self.grid.setData(self.x_plot, self.y_plot, xval.reshape(self.nelx, self.nely, order='F'))

    def set_conv_data(self, outit):
        """ Updates values of the objective function and the constraint function to plot the convergence graph.   """
        self.curves_funcs[0].setData(self.list_iter[:outit + 1], self.list_f0val[:outit + 1])

        for ind in range(len(self.constr_func)):
            self.curves_funcs[ind + 1].setData(self.list_iter[:outit+  1], self.list_fvals[:outit + 1, ind])

    def save_figs(self, graph_grid, graph_conv, graph_freq, canvas_mesh):
        """ Saves figure.
        Args:
            fig: Object with the graph. 
        """ 
        folder_name = 'images'
        directory = os.path.join(os.path.join(os.path.dirname(__file__), 'data'), folder_name)
        os.makedirs(directory, exist_ok=True)

        self.save_pg_fig(graph_grid, os.path.join(directory, 'part.png'))
        self.save_pg_fig(graph_conv.plotItem, os.path.join(directory, 'graph_conv.png'))

        if graph_freq is not None:
            self.save_pg_fig(graph_freq.plotItem, os.path.join(directory, 'graph_freq.png'))

        if canvas_mesh is not None:
            canvas_mesh.figure.savefig(os.path.join(directory, 'mesh.png'))

    def save_pg_fig(self, fig, path):
        """ Saves the graphics: optimized part, convergence graph, frequency response graph and deformed mesh.
        Args:
            fig (:obj:`pyqtgraph.graphicsItems.PlotItem`): Object with the graph.
            path: Directory to save the graph.
        """   
        exporter = pg.exporters.ImageExporter(fig)
        exporter.export(path)

    # def freqrsp_modes(freq_range, newf, oldf, modes, func_name, save):
    #     """ Plot the frequency response of the function with multiple modes.
    #     Args:
    #         freq_range (:obj:`list`): Range of frequencies analyzed.
    #             First value is the minimum frequency.
    #             Second value is the maximum frequency.
    #             Third value is the step between each calculation of the objective function. 
    #         newf (:obj:`numpy.array`): Optimized function.
    #         oldf (:obj:`numpy.array`): Original function.
    #         modes (:obj:`list`): Analyzed modes. 
    #         func_name (:obj:`str`): Objective function name.
    #             It can be: "compliance", "input_power", "elastic_potential_energy", "kinetic_energy" or "r_ratio".
    #         save (:obj:`bool`): True for save the graphic in PNG.
    #     """
    #     freq = np.arange(freq_range[0], freq_range[1] + 1, freq_range[2])
        
    #     for i, mode in enumerate(modes):
    #         plt.figure(i+1)
    #         plt.plot(freq, oldf.real, label='original')
    #         plt.plot(freq, newf[:, i].real, label=str(mode) + ' mode')
    #         plt.title(str(mode) + ' mode')
    #         plt.xlabel('frequency [Hz]', fontsize=16)
    #         plt.ylabel(func_name.lower(), fontsize=16)
    #         plt.yscale('log')
    #         plt.legend()
    #         if save:
    #             plt.savefig(str(mode) + ".eps")

    #     plt.figure(i+1)
    #     #plt.plot(freq, oldf.real, label='original')
    #     for i, mode in enumerate(modes):
    #         plt.plot(freq, newf[:, i].real, label=str(mode) + ' mode')
    #         plt.xlabel('frequency [Hz]', fontsize=16)
    #         plt.ylabel(func_name.lower(), fontsize=16)
    #         plt.yscale('log')
    #     plt.title('All modes')
    #     plt.legend()
    #     if save:
    #         plt.savefig('all' + ".eps")   
    #     plt.show()

    # def compare_deriv(nodes, delta_d, dw, dw_orig):
    #     for ind2, node in enumerate(nodes):
    #         plt.figure(ind2+1)
    #         plt.plot(delta_d, dw[:, ind2], marker='o', label='FDM')
    #         plt.plot(delta_d, np.repeat(dw_orig[ind2], len(delta_d)), marker='o', label='Analytical Method')
    #         plt.title('Node ' + str(node), fontsize=18)
    #         plt.xlabel(r'$\Delta D_i$', fontsize=18)
    #         plt.ylabel(r'$\frac{d \alpha}{d D_i}$', fontsize=18) #, rotation=0
    #         plt.legend()
    #     plt.show(block=True)
    
        # update_lists
        # create_header VAI NAS FUNÇÕES DA SCREEN MESMO OU PLOT, SE FOR NO PLOT TER DUAS LISTAS DE CADA, UMA PARA SALVAR E OUTRA PARA PLOTAR
        # printProgressBar
        # TODO: ULTIMAS FUNÇÕES PERGUNTAR ONDE O OLAVO ACHA MELHOR ENCAIXAR
        # TODO: PERGUNTAR SE ELE CHECOU COMO É COLOCADA AS CONDIÇÕES DE CONTORNO
        # TODO: SE ELE TEM UMA MANEIRA MELHOR E O QUE ACHA QUE PRECISA MUDAR