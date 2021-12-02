import numpy as np

class Xval:
    def __init__(self, nelx, nely, constr_func, constr_values, dens_filter, H, neighbors) -> None:
       
        value_xval = self.__set_value_initxval(constr_func, constr_values)
        self.xval = value_xval * np.ones((nelx * nely, 1))
        self.initial_xval = self.xval.copy()
        self.xnew = self.xval.copy()

        self.dens_filter = dens_filter

        self.H = H
        self.neighbors = neighbors

    def update_xval(self, new_xval):
        self.xval = new_xval
    
    def calc_xnew(self, xval=None):
        """ Recalculates xval. """
        if xval is None:
            xval = self.xval

        if self.dens_filter:
            a = 1/np.sum(self.H, axis=1)
            b = np.sum(self.H * xval[self.neighbors.flatten()].reshape(self.H.shape), axis=1)
            xe = a * b
            aux = xe.reshape(-1, 1)
            self.xnew = aux
        else:
            self.xnew = xval

    def __set_value_initxval(self, constr_func, constr_values):
        """ Calculates the initial value of xval.
        
        Args:
            constr_func (:obj:`list`): constraint functions applied.
            constr_values (:obj:`list`): Values of constraint functions applied.
        
        Returns:
            initial_xval (:obj:`float`): First value of xval
        """
        idx = [i for i, e in enumerate(constr_func) if e == constr_func[0]]
        if len(idx) > 1:
            initial_xval = ((abs(constr_values[idx[0]]) + abs(constr_values[idx[1]]))/2)/100
        else:
            initial_xval = abs(constr_values[idx[0]])/100
        return initial_xval
    
    def set_passive_el_xval(self, passive_el):
        """ Sets the values of passive elements.

        Args:
            passive_el (:obj:`numpy.array`): Passive element nodes.

        Returns:
            Updated xval.
        """
        self.xval[passive_el] = 1
      
    def set_passive_el_xmin(self, xmin, passive_el):
        """ Sets the values of passive elements.

        Args:
            xmin (:obj:`numpy.array`): 
            passive_el (:obj:`numpy.array`): Passive element nodes.

        Returns:
            Updated xmin.
        """
        xmin[passive_el] = 0.99
        return xmin