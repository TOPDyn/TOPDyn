from PyQt5 import QtCore, QtWidgets, QtGui
import os
import pickle 
import ast
import numpy as np

class Parameters():
    def update_params(self): 
        self.nelx = ast.literal_eval(self.nelx_spin.text())
        self.nely = ast.literal_eval(self.nely_spin.text())
        self.lx   = ast.literal_eval(self.lx_spin.text())
        self.ly   = ast.literal_eval(self.ly_spin.text())
        
        self.E = ast.literal_eval(self.E_spin.text())
        self.v = ast.literal_eval(self.v_spin.text())
        self.rho   = ast.literal_eval(self.rho_spin.text())
        
        self.alpha = ast.literal_eval(self.alpha_spin.text())
        self.beta  = ast.literal_eval(self.beta_spin.text())
        self.eta   = ast.literal_eval(self.eta_spin.text())
        
        self.factor = ast.literal_eval(self.factor_spin.text())
        self.freq = ast.literal_eval(self.freq_spin.text())

        self.freqrsp = self.freqrsp_check.checkState()
        self.freq_range = ast.literal_eval(self.freq_range_spin.text())
        
        self.load_matrix = ast.literal_eval(self.load_matrix_spin.text())
        self.constr_matrix = ast.literal_eval(self.constrain_spin.text())
        
        self.save = self.save_check.checkState()

    def check_param(self, input_val, types, war):
        try:
            for type in types:
                isinstance(ast.literal_eval(input_val), type)
        except ValueError:
            warning = QtWidgets.QLabel(war)
            self.warnings.append(warning)
        
    def check_params(self):
        self.warnings = []

        self.check_param(self.nelx_spin.text(), [int], 'nelx must be an intereger')

        self.check_param(self.nely_spin.text(), [int], 'nely must be an intereger')

        self.check_param(self.lx_spin.text(), [int, float], 'lx must be an intereger or float')

        self.check_param(self.ly_spin.text(), [int, float], 'ly must be an intereger or float')

        self.check_param(self.E_spin.text(), [int, float], "E must be an intereger or float")

        self.check_param(self.v_spin.text(), [int, float], 'v must be an intereger or float')

        self.check_param(self.rho_spin.text(), [int, float], 'rho must be an intereger or float')

        self.check_param(self.alpha_spin.text(), [int, float], 'alpha must be an intereger or float')

        self.check_param(self.beta_spin.text(), [int, float], 'beta must be an intereger or float')

        self.check_param(self.eta_spin.text(), [int, float], 'eta must be an intereger or float')

        self.check_param(self.factor_spin.text(), [int, float], 'Factor must be an intereger or float')

        self.check_param(self.freq_spin.text(), [int], 'Frequency must be an intereger')

        self.check_param(self.freq_range_spin.text(), [list], 'Frequency range must be a list')

class ParametersFEM2D(Parameters):
    def __init__(self):
        # QLineEdit
        self.nelx_spin = QtWidgets.QLineEdit()
        self.nely_spin = QtWidgets.QLineEdit()
        self.lx_spin = QtWidgets.QLineEdit()
        self.ly_spin = QtWidgets.QLineEdit()
        
        self.E_spin = QtWidgets.QLineEdit()
        self.v_spin = QtWidgets.QLineEdit()
        self.rho_spin = QtWidgets.QLineEdit()
        
        self.alpha_spin = QtWidgets.QLineEdit()
        self.beta_spin = QtWidgets.QLineEdit()
        self.eta_spin = QtWidgets.QLineEdit()
        
        self.freq_spin = QtWidgets.QLineEdit()
        self.load_matrix_spin = QtWidgets.QLineEdit()
        self.constrain_spin = QtWidgets.QLineEdit()

        self.factor_spin = QtWidgets.QLineEdit()

        self.freqrsp_check = QtWidgets.QCheckBox("Plot freq rsp")  
        self.freq_range_spin = QtWidgets.QLineEdit()
        self.node_plot_spin = QtWidgets.QLineEdit()

        self.freq_range_spin.setDisabled(True)
        self.node_plot_spin.setDisabled(True)

        self.save_check = QtWidgets.QCheckBox("Save data")

        self.warnings = []

        self._set_default()
        self.update_params()

        self.text = """ 
                    <p> <b> <font size="+7"> Parameters: </font> </b>
                    <hr>
                    <p> <b> <font size="+1"> nelx </font> </b> <font size="+1"> (<i>int</i>): Number of elements on the x-axis.  </font> </p>
                    <p> <b> <font size="+1"> nely </font> </b> <font size="+1"> (<i>int</i>): Number of elements on the y-axis.</font> </p>
                    <p> <b> <font size="+1"> lx </font> </b> <font size="+1"> (<i>float</i>): x-axis length.</font> </p>
                    <p> <b> <font size="+1"> ly </font> </b> <font size="+1"> (<i>float</i>): y-axis length.</font> </p>
                    <p> <b> <font size="+1"> load_matrix </font> </b> <font size="+1"> (<i>numpy.array</i>): It's a list of lists. </font> </p>
                    <p style="margin-left:2em"> <font size="+1"> The list can be:
                        <UL>
                        <LI>[x_coordinate, y_coordinate, force_applied_x, force_applied_y, force_value]</LI>
                        <LI>[value_coordinate, column_to_compare, force_applied_x, force_applied_y, force_value, error_margin]</LI>
                        </UL> </font> </p>
                    <p style="margin-left:2em"> <font size="+1"> It is possible to merge the two options. Examples:
                        <UL>
                        <LI>force_matrix = [[1, 1, 0, -1, 100]] -> Apply a negative force of modulus 100 N in the Y direction to the node at coordinate (1,1)</LI>
                        <LI>force_matrix = [[0, 1, 1, 0, 200, 0.001]] -> Apply a positive force of modulus 200 N in X direction to all the nodes with x=0</LI>
                        <LI>force_matrix = [[1, 1, 0, -1, 100], [0, 1, -1, 0, 200, 0.001]] -> Apply the two options above.</LI>
                        </UL> </font> </p>
                    <p> <b> <font size="+1"> constr_matrix </font> </b> <font size="+1"> (<i>numpy.array</i>): It's a list of lists. </font> </p>
                    <p style="margin-left:2em"> <font size="+1"> The list can be:
                        <UL>
                        <LI>[x_coordinate, y_coordinate, constrain_disp_x, constrain_disp_y]</LI>
                        <LI>[value_coordinate, column_to_compare, constrain_disp_x, constrain_disp_y, error_margin]</LI>
                        </UL> </font> </p>                    
                    <p> <b> <font size="+1"> E </font> </b> <font size="+1"> (<i>float</i>): Elastic modulus. </font> </p>
                    <p> <b> <font size="+1"> v </font> </b> <font size="+1"> (<i>float</i>): Poisson's ratio. </font> </p>
                    <p> <b> <font size="+1"> rho </font> </b> <font size="+1"> (<i>float</i>): Density. </font> </p>
                    <p> <b> <font size="+1"> alpha </font> </b> <font size="+1"> (<i>float</i>): Damping coefficient proportional to mass. </font> </p>
                    <p> <b> <font size="+1"> beta </font> </b> <font size="+1"> (<i>float</i>): Damping coefficient proportional to stiffness. </font> </p> 
                    <p> <b> <font size="+1"> eta </font> </b> <font size="+1"> (<i>float</i>): Damping coefficient. </font> </p>
                    <p> <b> <font size="+1"> factor </font> </b> <font size="+1"> (<i>float</i>): Factor to deform the mesh. </font> </p>
                    <p> <b> <font size="+1"> freq </font> </b> <font size="+1"> (<i>int</i>): Optimized frequency. </font> </p>
                    <p> <b> <font size="+1"> node_plot </font> </b> <font size="+1"> (<i>numpy.array</i>): The node to plot the frequency response graph.  </font> </p>
                        <p style="margin-left:2em"> <font size="+1"> The columns are respectively node, x direction, y direction. </p>
                        <p style="margin-left:2em"> <font size="+1"> It is a 1x3 matrix. </font> </p>
                    <p> <b> <font size="+1"> freq_rsp </font> </b> <font size="+1"> (<i>list</i>): If len is 3, a frequency response graph of the original and optimized structure is generated.  </font> </p>
                        <p style="margin-left:2em"> <font size="+1"> First value is the minimum frequency of the graph. </font> </p>
                        <p style="margin-left:2em"> <font size="+1"> Second value is the maximum frequency of the graph. </font> </p>
                        <p style="margin-left:2em"> <font size="+1"> Third value is the step between each calculation of the objective function. </font> <</p>
                    <p> <b> <font size="+1"> save </font> </b> <font size="+1"> (<i>bool</i>): if True save the optimization and frequency response graphs as PNG. </font> </p>
                    """

    def add_widgtes(self, layout):
        layout.addWidget(QtWidgets.QLabel('Nelx'))
        layout.addWidget(self.nelx_spin)

        layout.addWidget(QtWidgets.QLabel('Nely'))
        layout.addWidget(self.nely_spin)

        layout.addWidget(QtWidgets.QLabel('lx'))
        layout.addWidget(self.lx_spin)

        layout.addWidget(QtWidgets.QLabel('ly'))
        layout.addWidget(self.ly_spin)

        layout.addWidget(QtWidgets.QLabel('E'))
        layout.addWidget(self.E_spin)

        layout.addWidget(QtWidgets.QLabel('v'))
        layout.addWidget(self.v_spin)

        layout.addWidget(QtWidgets.QLabel('rho'))
        layout.addWidget(self.rho_spin)

        layout.addWidget(QtWidgets.QLabel('Alpha'))
        layout.addWidget(self.alpha_spin)

        layout.addWidget(QtWidgets.QLabel('Beta'))
        layout.addWidget(self.beta_spin)

        layout.addWidget(QtWidgets.QLabel('Eta'))
        layout.addWidget(self.eta_spin)

        layout.addWidget(QtWidgets.QLabel('Factor'))
        layout.addWidget(self.factor_spin)
        
        layout.addWidget(QtWidgets.QLabel('Frequency'))
        layout.addWidget(self.freq_spin)

        layout.addWidget(QtWidgets.QLabel('Force matrix'))
        layout.addWidget(self.load_matrix_spin)

        layout.addWidget(QtWidgets.QLabel('Contrain nodes disp.'))
        layout.addWidget(self.constrain_spin)

        layout.addWidget(self.freqrsp_check)

        layout.addWidget(QtWidgets.QLabel('Frequency range'))
        layout.addWidget(self.freq_range_spin)

        layout.addWidget(QtWidgets.QLabel('Node plot'))
        layout.addWidget(self.node_plot_spin)
        layout.addWidget(self.save_check)
    
    def _set_default(self):       
        self.nelx_spin.setText('40')
        self.nely_spin.setText('20')
        self.lx_spin.setText('1')
        self.ly_spin.setText('0.5')

        self.E_spin.setText('210e9')
        self.v_spin.setText('0.3')
        self.rho_spin.setText('7860')

        self.alpha_spin.setText('0')
        self.beta_spin.setText('0')
        self.eta_spin.setText('0')

        self.factor_spin.setText('10000')

        self.freq_spin.setText('200')

        self.freq_range_spin.setText("[0,25,5]")
        self.node_plot_spin.setText("[1, 0.25, 0, 1]")

        self.load_matrix_spin.setText('[{"x_coord":1, "y_coord":0.25, "x_direc":0, "y_direc":-1, "force":1}]')
        self.constrain_spin.setText('[{"coord":0, "axis":1, "eps":0.001, "constrain_disp_x":1, "constrain_disp_y":1}]')

    def update_params(self): 
        super().update_params()
        self.node_plot = ast.literal_eval(self.node_plot_spin.text())

    def export_param(self):
        param = {"nelx":self.nelx, "nely":self.nely, "lx":self.lx, "ly":self.ly, "E":self.E, "v":self.v, "rho":self.rho,
                "alpha":self.alpha, "beta":self.beta, "eta":self.eta, "factor":self.factor, "freq":self.freq, 
                "freqrsp":self.freqrsp, "freq_range":self.freq_range, "load_matrix":self.load_matrix, 
                "constr_matrix":self.constr_matrix, "node_plot":self.node_plot, "save":self.save, "mesh_file":None}
        
        folder_name = 'temp'
        directory = os.path.join(os.path.dirname(__file__), folder_name)
        os.makedirs(directory, exist_ok=True)

        try: 
            # WRITE
            dir_file = os.path.join(directory, 'param_file.txt')
            # open file for writing
            f = open(dir_file,"w")
            # write file
            f.write(str(param))
            # close file
            f.close()

        except: 
            print("Something went wrong")

    def check_params(self):
        super().check_params()
        self.check_param(self.node_plot_spin.text(), [list], 'Node must be a list')

class ParametersOpt(Parameters):
    def __init__(self):
        self.mma_radio = QtWidgets.QRadioButton("MMA")
        self.gcmma_radio = QtWidgets.QRadioButton("GMMA")

        # QLineEdit
        self.nelx_spin = QtWidgets.QLineEdit()
        self.nely_spin = QtWidgets.QLineEdit()
        self.lx_spin = QtWidgets.QLineEdit()
        self.ly_spin = QtWidgets.QLineEdit()
        
        self.E_spin = QtWidgets.QLineEdit()
        self.v_spin = QtWidgets.QLineEdit()
        self.rho_spin = QtWidgets.QLineEdit()
        self.fac_ratio_spin = QtWidgets.QLineEdit()

        self.x_min_m_spin = QtWidgets.QLineEdit()
        self.x_min_k_spin = QtWidgets.QLineEdit()
        self.penal_k_spin   = QtWidgets.QLineEdit() #p_par
        self.penal_m_spin   = QtWidgets.QLineEdit() #q_par
        
        self.alpha_spin = QtWidgets.QLineEdit()
        self.beta_spin = QtWidgets.QLineEdit()
        self.eta_spin = QtWidgets.QLineEdit()

        self.load_matrix_spin = QtWidgets.QLineEdit()
        self.constrain_spin = QtWidgets.QLineEdit()

        self.constr_func_spin = QtWidgets.QLineEdit()
        self.constr_values_spin = QtWidgets.QLineEdit()
        self.passive_coord_spin = QtWidgets.QLineEdit()    

        self.modes_spin = QtWidgets.QLineEdit()
        self.const_func_spin = QtWidgets.QLineEdit()
        self.n1_spin = QtWidgets.QLineEdit()
        self.freq_spin = QtWidgets.QLineEdit()
        self.func_name_box = QtWidgets.QComboBox()
        self.func_name_box.addItem("compliance")
        self.func_name_box.addItem("input_power")
        self.func_name_box.addItem("elastic_potential_energy")
        self.func_name_box.addItem("kinetic_energy")
        self.func_name_box.addItem("r_ratio")
        self.func_name_box.addItem("local_ep")
        self.func_name_box.addItem("local_ki")
        self.func_name_box.addItem("local_r")


        self.func_name2_box = QtWidgets.QComboBox()
        self.func_name2_box.addItem("None")
        self.func_name2_box.addItem("compliance")
        self.func_name2_box.addItem("input_power")
        self.func_name2_box.addItem("elastic_potential_energy")
        self.func_name2_box.addItem("kinetic_energy")
        self.func_name2_box.addItem("r_ratio")
        self.func_name2_box.addItem("local_ep")
        self.func_name2_box.addItem("local_ki")
        self.func_name2_box.addItem("local_r")
        self.freq2_spin = QtWidgets.QLineEdit()
        self.freq2_spin.setDisabled(True)

        self.freqrsp_check = QtWidgets.QCheckBox("Plot freq rsp")  
        self.freq_range_spin = QtWidgets.QLineEdit()
        self.alpha_plot_spin = QtWidgets.QLineEdit()
        self.beta_plot_spin  = QtWidgets.QLineEdit()
        self.eta_plot_spin   = QtWidgets.QLineEdit()

        self.freq_range_spin.setDisabled(True)
        self.alpha_plot_spin.setDisabled(True)
        self.beta_plot_spin.setDisabled(True)
        self.eta_plot_spin.setDisabled(True)

        self.max_iter_spin = QtWidgets.QLineEdit() 
        self.save_check = QtWidgets.QCheckBox("Save data")
        self.dens_filter_check = QtWidgets.QCheckBox("Density Filter")      
        self.mesh_deform_check = QtWidgets.QCheckBox("Plot deformed mesh")
        self.factor_spin = QtWidgets.QLineEdit()
        self.factor_spin.setDisabled(True)

        self.warnings = []

        self.text = """ <p> <b> <font size="+7"> Parameters: </font> </b>
                    <hr>
                    <p> <b> <font size="+1"> nelx </b> <font size="+1"> (<i>int</i>): Number of elements on the x-axis.  </font> </p>
                    <p> <b> <font size="+1"> nely </b> <font size="+1"> (<i>int</i>): Number of elements on the y-axis.  </font> </p>
                    <p> <b> <font size="+1"> lx </b> <font size="+1"> (<i>float</i>): x-axis length.  </font> </p>
                    <p> <b> <font size="+1"> ly </b> <font size="+1"> (<i>float</i>): x-axis length.  </font> </p>
                    <p> <b> <font size="+1"> func_name </b> <font size="+1"> (<i>str</i>): Objective function used.  </font> </p>
                        <p style="margin-left:2em"> <font size="+1"> It can be: 'Compliance', 'Input Power', 'Elastic Potential Energy', 'Kinetic Energy' or 'R Ratio'.  </font> </p>
                        <p style="margin-left:2em"> <font size="+1"> If the multiobjective function is being calculated, weight n1 is assigned.  </font> </p>
                        <p style="margin-left:2em"> <font size="+1"> The list can be: </font> </p>
                       
                            <p style="margin-left:4em"> <font size="+1"> [x_coordinate, y_coordinate, force_applied_x, force_applied_y, force_value] </font> </p>
                            <p style="margin-left:4em"> <font size="+1"> [value_coordinate, column_to_compare, force_applied_x, force_applied_y, force_value, error_margin] </font> </p>
                        
                        <p style="margin-left:2em"> <font size="+1"> It is possible to merge the two options. Examples: </font> </p>
                      
                            <p style="margin-left:4em"> <font size="+1"> force_matrix = [[1, 1, 0, -1, 100]] -> Apply a negative force of modulus 100 N in the Y direction to the node at coordinate (1,1) </font> </p>
                            <p style="margin-left:4em"> <font size="+1"> force_matrix = [[0, 1, 1, 0, 200, 0.001]] -> Apply a positive force of modulus 200 N in X direction to all the nodes with x=0 </font> </p>
                            <p style="margin-left:4em"> <font size="+1"> force_matrix = [[1, 1, 0, -1, 100], [0, 1, -1, 0, 200, 0.001]] -> Apply the two options above. </font> </p>
                        
                    <p> <b> <font size="+1"> constr_matrix </font> </b> <font size="+1"> (<i>numpy.array</i>): It's a list of lists. </font> </p>
                        <p style="margin-left:2em"> <font size="+1"> The list can be: </font> </p> 
        
                            <p style="margin-left:4em"> <font size="+1"> [x_coordinate, y_coordinate, constrain_disp_x, constrain_disp_y] </font> </p>
                            <p style="margin-left:4em"> <font size="+1"> [value_coordinate, column_to_compare, constrain_disp_x, constrain_disp_y, error_margin] </font> </p>
                        
                    <p> <b> <font size="+1"> freq1 </b> <font size="+1"> (<i>int</i>): Optimized frequency. </font> </p>
                    <p> <b> <font size="+1"> constr_func </b> <font size="+1"> (<i>list</i>): Constraint functions applied.  </font> </p>
                            <p style="margin-left:2em"> <font size="+1"> It can be: 'Area', 'R Ratio' or 'Compliance.  </font> </p>
                            <p style="margin-left:2em"> <font size="+1"> The first function in the list will be used to define the initial value of xval.  </font> </p>
                            <p style="margin-left:2em"> <font size="+1"> If the same function is passed 2x,the box constraint is used. Negative values indicate the lower constraint.  </font> </p>
                            <p style="margin-left:2em"> <font size="+1"> Example:  </font> </p>
                                <p style="margin-left:4em"> <font size="+1"> constr_func   = ['Area', 'Area']  </font> </p>
                                <p style="margin-left:4em"> <font size="+1"> constr_values = [50, -20]  </font> </p>
                    <p> <b> <font size="+1"> constr_values </b> <font size="+1"> (<i>list</i>): Values of constraint functions applied.   </font> </p>
                                <p style="margin-left:2em"> <font size="+1"> Value in position i relates to the function in position i of the list constr_func.  </font> </p>
                                <p style="margin-left:2em"> <font size="+1"> It can be a maximum of 6 values.  </font> </p>
                                <p style="margin-left:2em"> <font size="+1"> constr_values[i] &#60; 0 = lower constraint  </font> </p>
                                <p style="margin-left:2em"> <font size="+1"> constr_values[i] > 0 = upper constraint  </font> </p>
                                <p style="margin-left:2em"> <font size="+1"> If 'Compliance' is passed a tuple with constraint value and frequency respectively.  </font> </p>
                                <p style="margin-left:2em"> <font size="+1"> Example:   </font> </p>
                                    <p style="margin-left:4em"> <font size="+1"> constr_func   = ['Area', 'Area', 'Compliance, 'R Ratio]  </font> </p>
                                    <p style="margin-left:4em"> <font size="+1"> constr_values = [50, -20, (50, 1000), 10]  </font> </p>
                    <p> <b> <font size="+1"> n1 </b> <font size="+1"> (<i>float</i>): Weight n1 used in func_name. </font> </p>
                            <p style="margin-left:2em"> <font size="+1"> If n1 &#60; 0: Maximize objective function  </font> </p>
                            <p style="margin-left:2em"> <font size="+1"> If n1 > 0: Minimize objective function  </font> </p>
                    <p> <b> <font size="+1"> multiobjective </b> <font size="+1"> (<i>tuple</i>): Second function and frequency in the multiobjective.  </font> </p>
                            <p style="margin-left:2em"> <font size="+1"> First value is the second function of the multiobjective function. The assigned weight is (1 - n1).  </font> </p>
                            <p style="margin-left:2em"> <font size="+1"> Second value is the frequency that func_name2 is being optimized.   </font> </p>
                    <p> <b> <font size="+1"> const_func </b> <font size="+1"> (<i>float</i>):  </font> </p>
                    <p> <b> <font size="+1"> fac_ratio </b> <font size="+1"> (<i>float</i>): Factor applied in the radius to get elements in the vicinity of each element. </font> </p>
                    <p> <b> <font size="+1"> modes </b> <font size="+1"> (<i>int</i>): If not None is used the Mode Superposition Method to calculate the displacement. </font> </p>
                    <p> <b> <font size="+1"> rho </b> <font size="+1"> (<i>float</i>): Density.    </font> </p>
                    <p> <b> <font size="+1"> E </b> <font size="+1"> (<i>float</i>): Elastic modulus. </font> </p>
                    <p> <b> <font size="+1"> v </b> <font size="+1"> (<i>float</i>): Poisson's ratio. </font> </p>
                    <p> <b> <font size="+1"> x_min_m </b> <font size="+1"> (<i>float</i>): Minimum relative densities to mass.  </font> </p>
                    <p> <b> <font size="+1"> x_min_k </b> <font size="+1"> (<i>float</i>): Minimum relative densities to stiffness.  </font> </p>
                    <p> <b> <font size="+1"> alpha_par </b> <font size="+1"> (<i>float</i>): Damping coefficient proportional to mass.   </font> </p>
                    <p> <b> <font size="+1"> beta_par </b> <font size="+1"> (<i>float</i>): Damping coefficient proportional to stiffness.   </font> </p>
                    <p> <b> <font size="+1"> eta_par </b> <font size="+1"> (<i>float</i>): Damping coefficient.   </font> </p>
                    <p> <b> <font size="+1"> alpha_plot </b> <font size="+1"> (<i>float</i>): Damping coefficient proportional to mass to generate the graph. </font> </p>
                    <p> <b> <font size="+1"> beta_plot </b> <font size="+1"> (<i>float</i>): Damping coefficient proportional to stiffness to generate the graph.  </font> </p>
                    <p> <b> <font size="+1"> eta_plot </b> <font size="+1"> (<i>float</i>): Damping coefficient to generate the graph. </font> </p>
                    <p> <b> <font size="+1"> p_par </b> <font size="+1"> (<i>int</i>): Penalization power to stiffness.  </font> </p>
                    <p> <b> <font size="+1"> q_par </b> <font size="+1"> (<i>int</i>): Penalization power to mass.  </font> </p>
                    <p> <b> <font size="+1"> passive_coord </b> <font size="+1"> (<i>tuple</i>): Region that the shape will not be changed.   </font> </p>
                        <p style="margin-left:2em"> <font size="+1"> Example: </font> </p>
                        <p style="margin-left:4em"> <font size="+1"> ((0.5, 1), (0.3, 0.6)) = ((x_initial, x_final), (y_initial, y_final))  </font> </p>
                    <p> <b> <font size="+1"> freq_rsp </b> <font size="+1"> (<i>list</i>): If len is 3, a frequency response graph of the original and optimized structure is generated. </font> </p>
                        <p style="margin-left:2em"> <font size="+1"> First value is the minimum frequency of the graph.  </font> </p>
                        <p style="margin-left:2em"> <font size="+1"> Second value is the maximum frequency of the graph.  </font> </p>
                        <p style="margin-left:2em"> <font size="+1"> Third value is the step between each calculation of the objective function.   </font> </p>
                    <p> <b> <font size="+1"> dens_filter </b> <font size="+1"> (<i>bool</i>): If True use density filter and False use sensitivity filter.  </font> </p>
                    <p> <b> <font size="+1"> each_iter </b> <font size="+1"> (<i>bool</i>): If True plots the convergence graph for each iteration of the optimization.   </font> </p>
                    <p> <b> <font size="+1"> max_iter </b> <font size="+1"> (<i>int</i>): Number of iterations. </font> </p>
                    <p> <b> <font size="+1"> mesh_deform </b> <font size="+1"> (<i>bool</i>): If True plots the mesh deformation of the dynamic function.  </font> </p>
                    <p> <b> <font size="+1"> factor </b> <font size="+1"> (<i>float</i>): Factor to deform the mesh.  </font> </p>
                    <p> <b> <font size="+1"> save </b> <font size="+1"> (<i>bool</i>): if True save the optimization and frequency response graphs as PNG.  </font> </p>
                    """

        self._set_default()
        self.update_params()

    def export_param(self):
        param = {"nelx":self.nelx, "nely":self.nely, "lx":self.lx, "ly":self.ly, "E":self.E, "v":self.v, "rho":self.rho,
                "alpha_par":self.alpha, "beta_par":self.beta, "eta_par":self.eta, "factor":self.factor, "freq":self.freq, 
                "freqrsp":self.freqrsp, "freq_range":self.freq_range, "load_matrix":self.load_matrix,
                "constr_matrix":self.constr_matrix, "save":self.save, "mesh_file":None, "mma":self.mma,
                "fac_ratio":self.fac_ratio, "x_min_m":self.x_min_m, "x_min_k":self.x_min_k, "penal_k":self.penal_k,
                "penal_m":self.penal_m, "constr_func":self.constr_func, "constr_values":self.constr_values,
                "passive_coord":self.passive_coord, "modes":self.modes,"const_func":self.const_func,"n1":self.n1,
                "func_name":self.func_name,"func_name2":self.func_name2,"freq2":self.freq2,"alpha_plot":self.alpha_plot,
                "beta_plot":self.beta_plot,"eta_plot":self.eta_plot,"max_iter":self.max_iter,
                "save":self.const_func,"save":self.save,"dens_filter":self.dens_filter,"mesh_deform":self.mesh_deform}
        
        folder_name = 'temp'
        directory = os.path.join(os.path.dirname(__file__), folder_name)
        os.makedirs(directory, exist_ok=True)

        try: 
            # WRITE
            dir_file = os.path.join(directory, 'param_file.txt')
            # open file for writing
            f = open(dir_file,"w")
            # write file
            f.write(str(param))
            # close file
            f.close()

        except: 
            print("Something went wrong")

    def add_widgtes(self, layout):
        layout.addWidget(QtWidgets.QLabel('Optimization Method'))
        layout.addWidget(self.mma_radio)
        layout.addWidget(self.gcmma_radio)

        layout.addWidget(QtWidgets.QLabel('Nelx'))
        layout.addWidget(self.nelx_spin)
        
        layout.addWidget(QtWidgets.QLabel('Nely'))
        layout.addWidget(self.nely_spin)

        layout.addWidget(QtWidgets.QLabel('lx'))
        layout.addWidget(self.lx_spin)

        layout.addWidget(QtWidgets.QLabel('ly'))
        layout.addWidget(self.ly_spin)

        layout.addWidget(QtWidgets.QLabel('E'))
        layout.addWidget(self.E_spin)

        layout.addWidget(QtWidgets.QLabel('v'))
        layout.addWidget(self.v_spin)

        layout.addWidget(QtWidgets.QLabel('rho'))
        layout.addWidget(self.rho_spin)

        layout.addWidget(QtWidgets.QLabel('Factor ratio'))
        layout.addWidget(self.fac_ratio_spin)

        layout.addWidget(QtWidgets.QLabel('x_min_mass'))
        layout.addWidget(self.x_min_m_spin)

        layout.addWidget(QtWidgets.QLabel('x_min_stif'))
        layout.addWidget(self.x_min_k_spin)

        layout.addWidget(QtWidgets.QLabel('penal stif'))
        layout.addWidget(self.penal_k_spin)

        layout.addWidget(QtWidgets.QLabel('penal mass'))
        layout.addWidget(self.penal_m_spin)        

        layout.addWidget(QtWidgets.QLabel('Alpha'))
        layout.addWidget(self.alpha_spin)

        layout.addWidget(QtWidgets.QLabel('Beta'))
        layout.addWidget(self.beta_spin)

        layout.addWidget(QtWidgets.QLabel('Eta'))
        layout.addWidget(self.eta_spin)

        layout.addWidget(QtWidgets.QLabel('Load matrix'))
        layout.addWidget(self.load_matrix_spin)

        layout.addWidget(QtWidgets.QLabel('Constraint nodes'))
        layout.addWidget(self.constrain_spin)

        layout.addWidget(QtWidgets.QLabel('Constrain funcs.'))
        layout.addWidget(self.constr_func_spin)

        layout.addWidget(QtWidgets.QLabel('Constrain values'))
        layout.addWidget(self.constr_values_spin)

        layout.addWidget(QtWidgets.QLabel('Passive coords'))
        layout.addWidget(self.passive_coord_spin)

        layout.addWidget(QtWidgets.QLabel('Modes'))
        layout.addWidget(self.modes_spin)

        layout.addWidget(QtWidgets.QLabel('Constant Function'))
        layout.addWidget(self.const_func_spin)

        layout.addWidget(QtWidgets.QLabel('n1'))
        layout.addWidget(self.n1_spin)

        layout.addWidget(QtWidgets.QLabel('Frequency'))
        layout.addWidget(self.freq_spin)

        layout.addWidget(QtWidgets.QLabel('Objective Function'))
        layout.addWidget(self.func_name_box)

        layout.addWidget(QtWidgets.QLabel('Multiobjective Function'))
        layout.addWidget(self.func_name2_box)
        layout.addWidget(QtWidgets.QLabel('Multiobjective Freq'))
        layout.addWidget(self.freq2_spin)

        layout.addWidget(QtWidgets.QLabel('Max iterations'))
        layout.addWidget(self.max_iter_spin)

        layout.addWidget(self.save_check)
        layout.addWidget(self.dens_filter_check)

        layout.addWidget(self.mesh_deform_check)
        layout.addWidget(QtWidgets.QLabel('Factor'))
        layout.addWidget(self.factor_spin)
      
        layout.addWidget(self.freqrsp_check)
        layout.addWidget(QtWidgets.QLabel('Freq range'))
        layout.addWidget(self.freq_range_spin)

        layout.addWidget(QtWidgets.QLabel('Alpha plot'))
        layout.addWidget(self.alpha_plot_spin)

        layout.addWidget(QtWidgets.QLabel('Beta plot'))
        layout.addWidget(self.beta_plot_spin)

        layout.addWidget(QtWidgets.QLabel('Eta plot'))
        layout.addWidget(self.eta_plot_spin)

    def update_params(self): 
        super().update_params()

        self.mma = self.check_method()
        self.fac_ratio = ast.literal_eval(self.fac_ratio_spin.text())

        self.x_min_m = ast.literal_eval(self.x_min_m_spin.text())
        self.x_min_k = ast.literal_eval(self.x_min_k_spin.text())
        self.penal_k = ast.literal_eval(self.penal_k_spin.text()) #p_par
        self.penal_m = ast.literal_eval(self.penal_m_spin.text()) #q_par

        self.constr_func = ast.literal_eval(self.constr_func_spin.text())
        self.constr_values = ast.literal_eval(self.constr_values_spin.text())
        self.passive_coord = ast.literal_eval(self.passive_coord_spin.text())

        self.modes = ast.literal_eval(self.modes_spin.text())
        self.const_func = ast.literal_eval(self.const_func_spin.text())
        self.n1 = ast.literal_eval(self.n1_spin.text())
        self.func_name = str(self.func_name_box.currentText())

        self.func_name2 = str(self.func_name2_box.currentText())
        self.freq2 = ast.literal_eval(self.freq2_spin.text())

        self.alpha_plot = ast.literal_eval(self.alpha_plot_spin.text())
        self.beta_plot  = ast.literal_eval(self.beta_plot_spin.text())
        self.eta_plot   = ast.literal_eval(self.eta_plot_spin.text())

        self.max_iter = ast.literal_eval(self.max_iter_spin.text())
        self.save = self.save_check.checkState()

        self.dens_filter  = self.dens_filter_check.checkState()
        self.mesh_deform  = self.mesh_deform_check.checkState()
       
    def check_method(self):
        if self.mma_radio.isChecked():
            return True
        else:
            return False

    def _set_default(self):
        self.mma_radio.setChecked(True)

        self.nelx_spin.setText('50')
        self.nely_spin.setText('100')
        self.lx_spin.setText('0.5')
        self.ly_spin.setText('1')

        self.E_spin.setText('210e9')
        self.v_spin.setText('0.3')
        self.rho_spin.setText('7860')
        self.fac_ratio_spin.setText('2.2')

        self.x_min_m_spin.setText('1e-12')
        self.x_min_k_spin.setText('1e-9')
        self.penal_k_spin.setText('3')
        self.penal_m_spin.setText('1')

        self.alpha_spin.setText('0')
        self.beta_spin.setText('1e-5')
        self.eta_spin.setText('0')

        self.load_matrix_spin.setText('[{"coord":1, "axis":2, "x_direc":0, "y_direc":-1, "force":10000, "eps":0.001}]')
        self.constrain_spin.setText('[{"coord":0, "axis":2, "eps":0.001, "constrain_disp_x":1, "constrain_disp_y":1}]')

        self.constr_func_spin.setText("['area']")
        self.constr_values_spin.setText('[30]')
        self.passive_coord_spin.setText('((0, 0.5), (0.95, 1))')

        self.modes_spin.setText('None')
        self.const_func_spin.setText('100')
        self.n1_spin.setText('1')
        self.freq_spin.setText('10')
        
        self.freq2_spin.setText('0')
        
        self.freq_range_spin.setText("[5, 500, 5]")
        self.alpha_plot_spin.setText('0')
        self.beta_plot_spin.setText('1e-6')
        self.eta_plot_spin.setText('0')

        self.max_iter_spin.setText('100')
        self.factor_spin.setText('10000')       

    def check_params(self):
        super().check_params()

        if str(self.func_name2_box.currentText()) is not None:
            self.check_param(self.freq2_spin.text(), [int], 'Frequency must be an intereger')

        self.check_param(self.x_min_m_spin.text(), [int, float], 'x_min_m must be an intereger or float')

        self.check_param(self.x_min_k_spin.text(), [int, float], 'x_min_k must be an intereger or float')

        self.check_param(self.penal_k_spin.text(), [int], 'penal_k must be an intereger')

        self.check_param(self.penal_m_spin.text(), [int], 'penal_m must be an intereger')

        self.check_param(self.const_func_spin.text(), [int, float], 'const_func must be an intereger or float')

        self.check_param(self.n1_spin.text(), [int, float], 'n1 must be an intereger or float')

        self.check_param(self.fac_ratio_spin.text(), [int, float], 'fac_ratio must be an intereger or float')

        if not isinstance(ast.literal_eval(self.modes_spin.text()), int) and not (ast.literal_eval(self.modes_spin.text()) is None): 
            warning = QtWidgets.QLabel('modes must be an intereger or None')
            self.warnings.append(warning)

        self.check_param(self.max_iter_spin.text(), [int, float], 'max_iter must be an intereger or float')

        self.check_param(self.alpha_plot_spin.text(), [int, float], 'alpha_plot must be an intereger or float')

        self.check_param(self.beta_plot_spin.text(), [int, float], 'beta_plot must be an intereger or float')

        self.check_param(self.eta_plot_spin.text(), [int, float], 'eta_plot must be an intereger or float')

        # TODO: FALTA COLOCAR ELES
        # self.constr_func_par = ast.literal_eval(self.constr_func_spin.text())
        # self.constr_values_par = ast.literal_eval(self.constr_values_spin.text())
        # self.passive_coord_par = ast.literal_eval(self.passive_coord_spin.text())

class ParametersText():
    def __init__(self):
        self.editor = QtWidgets.QLabel()
        self.editor.setTextInteractionFlags(QtCore.Qt.LinksAccessibleByMouse)
        self.editor.setTextFormat(QtCore.Qt.RichText)
        self.editor.setTextInteractionFlags(QtCore.Qt.TextBrowserInteraction)
        self.editor.setOpenExternalLinks(True)

    def set_text(self, text):
        self.editor.setText(text)