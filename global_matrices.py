
import numpy as np

class GlobalMatrices:
    def __init__(self, coord, connect, E, v, rho):
        self.coord = coord
        self.connect = connect
        self.E = E
        self.v = v
        self.rho = rho

    def shape_Q4(self, ssx, ttx):
        """ Linear Shape Functions and Derivatives.

        Args:
            ssx (:obj:`float`): Local coordinate of the element on the X-axis.
            ttx (:obj:`float`): Local coordinate of the element on the Y-axis.
        
        Returns:
            Tuple with the shape function and its derivative.
        """
        numerator = 4
        #shape functions
        phi = np.array([(1-ssx)*(1-ttx), (1+ssx)*(1-ttx), (1+ssx)*(1+ttx), (1-ssx)*(1+ttx)], dtype=float)/numerator
        #derivatives
        dphi = np.array([[-(1-ttx),  (1-ttx), (1+ttx), -(1+ttx)],
                        [-(1-ssx), -(1+ssx), (1+ssx),  (1-ssx)]], dtype=float)/numerator
        return phi, dphi

    def matrices_Q4(self, ee):
        """ Q4 stiffness and mass matrices.

        Args:
            ee (:obj:`int`): Element.
            coord (:obj:`numpy.array`): Coordinates of the element.
            connect (:obj:`numpy.array`): Element connectivity.
            E (:obj:`float`): Elastic modulus.
            v (:obj:`float`): Poisson's ratio.
            rho (:obj:`float`): Density.

        Returns:
            Tuple with stiffness and mass matrices.
        """
        # constitutive info
        tempc = self.E/(1 - self.v**2)
        tempn = tempc*(1 - self.v)/2
        tempt = tempc*self.v
        cttv = np.array([[tempc, tempt, 0], [tempt, tempc, 0], [0, 0, tempn]], dtype=float)
        # integration points
        nint, con, wps = 4, 1/np.sqrt(3), 1
        pint = np.array([[-con, -con], [con, -con], [con, con], [-con, con]], dtype=float)
        # preallocating elementary matrices
        Ke, Me = 0, 0
        # integration
        for i in range(nint):
            ssx, ttx = pint[i, 0], pint[i, 1]
            phi, dphi = self.shape_Q4(ssx,ttx)
            ie = self.connect[ee,1:]-1
            dxdy = dphi@self.coord[ie, 1:3] # note: dxds, dyds, dxdt, dydt = dxdy[0,0], dxdy[0,1], dxdy[1,0],dxdy[1,1] 
            JAC = np.array([[dxdy[0,0], dxdy[0,1]],[dxdy[1,0], dxdy[1,1]]], dtype=float) # JAC = np.array([[dxds, dyds],[dxdt, dydt]], dtype=float)
            detJAC = JAC[0,0]*JAC[1,1] - JAC[0,1]*JAC[1,0]
            iJAC = (1/detJAC)*np.array([[JAC[1,1], -JAC[0,1]], [-JAC[1,0], JAC[0,0]]], dtype=float)
            dphi_t = iJAC @ dphi
            
            B = np.array([[dphi_t[0,0], 0, dphi_t[0,1], 0, dphi_t[0,2], 0, dphi_t[0,3], 0],
                        [0, dphi_t[1,0], 0, dphi_t[1,1], 0, dphi_t[1,2], 0, dphi_t[1,3]],
                        [dphi_t[1,0], dphi_t[0,0], dphi_t[1,1], dphi_t[0,1],dphi_t[1,2], dphi_t[0,2], dphi_t[1,3], dphi_t[0,3]]], dtype=float)
            
            N = np.array([[phi[0], 0, phi[1], 0, phi[2], 0, phi[3], 0],
                        [0, phi[0], 0, phi[1], 0, phi[2], 0, phi[1]]], dtype=float)
            
            Ke += B.T@(cttv@B)*(detJAC*wps)
            Me += self.rho*N.T@N*(detJAC*wps)
        return Ke.real, Me.real