from PyQt5 import QtCore, QtWidgets, QtGui
import os
import ast
import numpy as np

class Parameters():
    def __init__(self) -> None:
        # btn load
        self.load_type_btn = [] # by_coord ou by_col
        self.load_coord_btn = [] # uma tupla com (x,y) ou um unico
        self.load_col_btn = [] # valores None
        self.load_error_btn = [] # valores None
        self.load_x_dir_btn = [] # tuplas ou valores sozinhos
        self.load_y_dir_btn = [] # tuplas ou valores sozinhos
        self.load_value_btn = []

        # values load
        self.by_coord = [] # by_coord ou by_col
        self.load_coord = [] # uma tupla com (x,y) ou um unico
        self.x_load_col = [] # valores None
        self.load_error = [] # valores None
        self.load_x_dir = [] # tuplas ou valores sozinhos
        self.load_y_dir = [] # tuplas ou valores sozinhos
        self.load_value = []

        self.warnings = []
        self.warnings_load = []

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
                
        self.save = self.save_check.checkState()

    def check_param(self, warnings, input_val, types, war):
        try:
            if input_val:
                for type in types:
                    isinstance(ast.literal_eval(input_val), type)
            else:
                warning = QtWidgets.QLabel(war)
                warnings.append(warning)

        except ValueError:
            warning = QtWidgets.QLabel(war)
            warnings.append(warning)
        
    def check_params(self):
        self.warnings = []

        self.check_param(self.warnings, self.nelx_spin.text(), [int], 'nelx must be an integer')

        self.check_param(self.warnings, self.nely_spin.text(), [int], 'nely must be an integer')

        self.check_param(self.warnings, self.lx_spin.text(), [int, float], 'lx must be an integer or float')

        self.check_param(self.warnings, self.ly_spin.text(), [int, float], 'ly must be an integer or float')

        self.check_param(self.warnings, self.E_spin.text(), [int, float], "E must be an integer or float")

        self.check_param(self.warnings, self.v_spin.text(), [int, float], 'v must be an integer or float')

        self.check_param(self.warnings, self.rho_spin.text(), [int, float], 'rho must be an integer or float')

        self.check_param(self.warnings, self.alpha_spin.text(), [int, float], 'alpha must be an integer or float')

        self.check_param(self.warnings, self.beta_spin.text(), [int, float], 'beta must be an integer or float')

        self.check_param(self.warnings, self.eta_spin.text(), [int, float], 'eta must be an integer or float')

        self.check_param(self.warnings, self.factor_spin.text(), [int, float], 'Factor must be an integer or float')

        self.check_param(self.warnings, self.freq_spin.text(), [int], 'Frequency must be an integer')

        self.check_param(self.warnings, self.freq_range_spin.text(), [list], 'Frequency range must be a list')

class ParametersOpt(Parameters):
    def __init__(self):
        super().__init__()

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

        self.create_param_btns()
        self.set_default()
        #self.update_params()

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

# Optimization param
    def create_param_btns(self):
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
        #self.freq2_spin.setDisabled(True)

        self.freqrsp_check = QtWidgets.QCheckBox("Plot freq rsp")  
        self.freq_range_spin = QtWidgets.QLineEdit()
        self.alpha_plot_spin = QtWidgets.QLineEdit()
        self.beta_plot_spin  = QtWidgets.QLineEdit()
        self.eta_plot_spin   = QtWidgets.QLineEdit()

        self.max_iter_spin = QtWidgets.QLineEdit() 
        self.save_check = QtWidgets.QCheckBox("Save data")
        self.dens_filter_check = QtWidgets.QCheckBox("Density Filter")      
        self.mesh_deform_check = QtWidgets.QCheckBox("Plot deformed mesh")
        self.factor_spin = QtWidgets.QLineEdit()
    
    def add_btns(self, layout):
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

        self.mma = True if self.mma_radio.isChecked() else False
        self.fac_ratio = ast.literal_eval(self.fac_ratio_spin.text())

        self.x_min_m = ast.literal_eval(self.x_min_m_spin.text())
        self.x_min_k = ast.literal_eval(self.x_min_k_spin.text())
        self.penal_k = ast.literal_eval(self.penal_k_spin.text()) #p_par
        self.penal_m = ast.literal_eval(self.penal_m_spin.text()) #q_par

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

    def set_default(self):
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
        
        self.freq2_spin.setDisabled(True)
        self.freq_range_spin.setDisabled(True)
        self.alpha_plot_spin.setDisabled(True)
        self.beta_plot_spin.setDisabled(True)
        self.eta_plot_spin.setDisabled(True)
        self.factor_spin.setDisabled(True)

        self.mesh_deform_check.toggled.connect(self.factor_spin.setEnabled)   
        self.func_name2_box.activated.connect(self.freq2_spin.setEnabled) 
        self.freqrsp_check.toggled.connect(self.freq_range_spin.setEnabled)
        self.freqrsp_check.toggled.connect(self.alpha_plot_spin.setEnabled)
        self.freqrsp_check.toggled.connect(self.beta_plot_spin.setEnabled)
        self.freqrsp_check.toggled.connect(self.eta_plot_spin.setEnabled)     

    def update_default(self):
        self.mma_radio.setChecked(self.mma)
        self.gcmma_radio.setChecked(not self.mma)

        self.nelx_spin.setText(str(self.nelx))
        self.nely_spin.setText(str(self.nely))
        self.lx_spin.setText(str(self.lx))
        self.ly_spin.setText(str(self.ly))

        self.E_spin.setText(str(self.E))
        self.v_spin.setText(str(self.v))
        self.rho_spin.setText(str(self.rho))
        self.fac_ratio_spin.setText(str(self.fac_ratio))

        self.x_min_m_spin.setText(str(self.x_min_m))
        self.x_min_k_spin.setText(str(self.x_min_k))
        self.penal_k_spin.setText(str(self.penal_k))
        self.penal_m_spin.setText(str(self.penal_m))

        self.alpha_spin.setText(str(self.alpha))
        self.beta_spin.setText(str(self.beta))
        self.eta_spin.setText(str(self.eta))

        self.passive_coord_spin.setText(str(self.passive_coord))
        self.modes_spin.setText(str(self.modes))
        self.const_func_spin.setText(str(self.const_func))
        self.n1_spin.setText(str(self.n1))
        self.freq_spin.setText(str(self.freq))
        self.freq2_spin.setText(str(self.freq2))
        
        self.freq_range_spin.setText(str(self.freq_range))
        self.alpha_plot_spin.setText(str(self.alpha_plot))
        self.beta_plot_spin.setText(str(self.beta_plot))
        self.eta_plot_spin.setText(str(self.eta_plot))

        self.max_iter_spin.setText(str(self.max_iter))
        self.factor_spin.setText(str(self.factor))

        self.func_name_box.setCurrentText(self.func_name)  

        if self.freqrsp:
            self.freqrsp_check.setChecked(True) 
            self.freq_range_spin.setEnabled(True)
            self.alpha_plot_spin.setEnabled(True)
            self.beta_plot_spin.setEnabled(True)
            self.eta_plot_spin.setEnabled(True)
        else:
            self.freq_range_spin.setDisabled(True)
            self.alpha_plot_spin.setDisabled(True)
            self.beta_plot_spin.setDisabled(True)
            self.eta_plot_spin.setDisabled(True)

        if self.mesh_deform:
            self.mesh_deform_check.setChecked(True)
            self.factor_spin.setEnabled(True)

        if self.save:
            self.save_check.setChecked(True)

        if self.dens_filter:
            self.dens_filter_check.setChecked(True)

        if self.func_name2 != 'None':
            self.freq2_spin.setEnabled(True)
            self.func_name2_box.setCurrentText(self.func_name2)        

        self.mesh_deform_check.toggled.connect(self.factor_spin.setEnabled)   
        self.func_name2_box.activated.connect(self.freq2_spin.setEnabled) 
        self.freqrsp_check.toggled.connect(self.freq_range_spin.setEnabled)
        self.freqrsp_check.toggled.connect(self.alpha_plot_spin.setEnabled)
        self.freqrsp_check.toggled.connect(self.beta_plot_spin.setEnabled)
        self.freqrsp_check.toggled.connect(self.eta_plot_spin.setEnabled)  

    def check_params(self):
        super().check_params()

        if str(self.func_name2_box.currentText()) is not None:
            self.check_param(self.warnings, self.freq2_spin.text(), [int], 'Frequency must be an integer')

        self.check_param(self.warnings, self.x_min_m_spin.text(), [int, float], 'x_min_m must be an integer or float')

        self.check_param(self.warnings, self.x_min_k_spin.text(), [int, float], 'x_min_k must be an integer or float')

        self.check_param(self.warnings, self.penal_k_spin.text(), [int], 'penal_k must be an integer')

        self.check_param(self.warnings, self.penal_m_spin.text(), [int], 'penal_m must be an integer')

        self.check_param(self.warnings, self.const_func_spin.text(), [int, float], 'const_func must be an integer or float')

        self.check_param(self.warnings, self.n1_spin.text(), [int, float], 'n1 must be an integer or float')

        self.check_param(self.warnings, self.fac_ratio_spin.text(), [int, float], 'fac_ratio must be an integer or float')

        if not isinstance(ast.literal_eval(self.modes_spin.text()), int) and not (ast.literal_eval(self.modes_spin.text()) is None): 
            warning = QtWidgets.QLabel('modes must be an integer or None')
            self.warnings.append(warning)

        self.check_param(self.warnings, self.max_iter_spin.text(), [int, float], 'max_iter must be an integer or float')

        self.check_param(self.warnings, self.alpha_plot_spin.text(), [int, float], 'alpha_plot must be an integer or float')

        self.check_param(self.warnings, self.beta_plot_spin.text(), [int, float], 'beta_plot must be an integer or float')

        self.check_param(self.warnings, self.eta_plot_spin.text(), [int, float], 'eta_plot must be an integer or float')

# Load param
    def create_param_load(self):
        by_coord_load_btn = QtWidgets.QRadioButton("Add load by coordinate")
        by_column_load_btn = QtWidgets.QRadioButton("Add load by column")
        self.load_type_btn.append([by_coord_load_btn, by_column_load_btn])

        coord_load_x_line  = QtWidgets.QLineEdit()
        coord_load_y_line  = QtWidgets.QLineEdit()
        coord_load_line = QtWidgets.QLineEdit()
        self.load_coord_btn.append([coord_load_x_line, coord_load_y_line, coord_load_line])
        
        x_by_col_load_btn = QtWidgets.QRadioButton("X")
        y_by_col_load_btn = QtWidgets.QRadioButton("Y")
        self.load_col_btn.append([x_by_col_load_btn, y_by_col_load_btn])

        erro_margin_load_line  = QtWidgets.QLineEdit()
        self.load_error_btn.append(erro_margin_load_line)

        load_x_dir_line_pos = QtWidgets.QRadioButton("Positive") 
        load_x_dir_line_neg = QtWidgets.QRadioButton("Negative")
        load_x_dir_line_nan = QtWidgets.QRadioButton("None")
        self.load_x_dir_btn.append([load_x_dir_line_pos, load_x_dir_line_neg, load_x_dir_line_nan])
   
        load_y_dir_line_pos  = QtWidgets.QRadioButton("Positive") 
        load_y_dir_line_neg  = QtWidgets.QRadioButton("Negative")         
        load_y_dir_line_nan = QtWidgets.QRadioButton("None")
        self.load_y_dir_btn.append([load_y_dir_line_pos, load_y_dir_line_neg, load_y_dir_line_nan])

        load_value_line = QtWidgets.QLineEdit()
        self.load_value_btn.append(load_value_line)

    def add_param_load(self, layout):
        load_type = QtWidgets.QButtonGroup(layout)
        load_label = QtWidgets.QLabel('---- Load ----')
        load_label.setAlignment(QtCore.Qt.AlignCenter)
        layout.addRow(load_label)
        
        load_type.addButton(self.load_type_btn[-1][0])
        layout.addRow(self.load_type_btn[-1][0])
       
        layout.addRow(QtWidgets.QLabel('X-coor'), self.load_coord_btn[-1][0])
        layout.addRow(QtWidgets.QLabel('Y-coord'), self.load_coord_btn[-1][1])
       
        load_type.addButton(self.load_type_btn[-1][1])
        layout.addRow(self.load_type_btn[-1][1])
       
        layout.addRow(QtWidgets.QLabel('Coord'), self.load_coord_btn[-1][2])
     
        layout.addRow(QtWidgets.QLabel('Column'))
        column_load = QtWidgets.QButtonGroup(layout)
        column_load.addButton(self.load_col_btn[-1][0])
        column_load.addButton(self.load_col_btn[-1][1])
        layout.addRow(self.load_col_btn[-1][0], self.load_col_btn[-1][1])
      
        layout.addRow(QtWidgets.QLabel('Error margin'), self.load_error_btn[-1])
       
        x_dir = QtWidgets.QButtonGroup(layout)
        x_dir.addButton(self.load_x_dir_btn[-1][0])
        x_dir.addButton(self.load_x_dir_btn[-1][1])
        x_dir.addButton(self.load_x_dir_btn[-1][2])
        layout.addRow(QtWidgets.QLabel('Add load in X direction'))
        lay = QtWidgets.QHBoxLayout()
        lay.addWidget(self.load_x_dir_btn[-1][0])
        lay.addWidget(self.load_x_dir_btn[-1][1])
        lay.addWidget(self.load_x_dir_btn[-1][2])
        layout.addRow(lay)

        y_dir = QtWidgets.QButtonGroup(layout)
        y_dir.addButton(self.load_y_dir_btn[-1][0])
        y_dir.addButton(self.load_y_dir_btn[-1][1])
        y_dir.addButton(self.load_y_dir_btn[-1][2])
        layout.addRow(QtWidgets.QLabel('Add load in Y direction'))
        lay = QtWidgets.QHBoxLayout()
        lay.addWidget(self.load_y_dir_btn[-1][0])
        lay.addWidget(self.load_y_dir_btn[-1][1])
        lay.addWidget(self.load_y_dir_btn[-1][2])
        layout.addRow(lay)

        layout.addRow(QtWidgets.QLabel('Force value'), self.load_value_btn[-1])

    def set_default_load(self):
        self.load_type_btn[-1][0].setChecked(True)
        self.load_type_btn[-1][1].setChecked(False)

        self.load_coord_btn[-1][0].setText(str(self.lx))
        self.load_coord_btn[-1][1].setText(str(self.ly/2))
        
        self.load_coord_btn[-1][2].setDisabled(True)
        self.load_col_btn[-1][0].setDisabled(True)
        self.load_col_btn[-1][1].setDisabled(True)
        self.load_error_btn[-1].setDisabled(True)

        self.load_x_dir_btn[-1][0].setChecked(True)
        self.load_y_dir_btn[-1][0].setChecked(True)
        self.load_value_btn[-1].setText('100')

        self.load_type_btn[-1][0].toggled.connect(self.load_coord_btn[-1][0].setEnabled)  
        self.load_type_btn[-1][0].toggled.connect(self.load_coord_btn[-1][1].setEnabled)
        self.load_type_btn[-1][0].toggled.connect(self.load_col_btn[-1][0].setDisabled)
        self.load_type_btn[-1][0].toggled.connect(self.load_col_btn[-1][1].setDisabled)
        self.load_type_btn[-1][0].toggled.connect(self.load_error_btn[-1].setDisabled)
        self.load_type_btn[-1][0].toggled.connect(self.load_coord_btn[-1][2].setDisabled)

        self.load_type_btn[-1][1].toggled.connect(self.load_coord_btn[-1][0].setDisabled)  
        self.load_type_btn[-1][1].toggled.connect(self.load_coord_btn[-1][1].setDisabled)
        self.load_type_btn[-1][1].toggled.connect(self.load_col_btn[-1][0].setEnabled)
        self.load_type_btn[-1][1].toggled.connect(self.load_col_btn[-1][1].setEnabled)
        self.load_type_btn[-1][1].toggled.connect(self.load_error_btn[-1].setEnabled)
        self.load_type_btn[-1][1].toggled.connect(self.load_coord_btn[-1][2].setEnabled)

    def reset_load_list(self):
        self.load_type_btn = [] # by_coord ou by_col
        self.load_coord_btn = [] # uma tupla com (x,y) ou um unico
        self.load_col_btn = [] # valores None
        self.load_error_btn = [] # valores None
        self.load_x_dir_btn = [] # tuplas ou valores sozinhos
        self.load_y_dir_btn = [] # tuplas ou valores sozinhos
        self.load_value_btn = []

        # values load
        self.by_coord = [] # by_coord ou by_col
        self.load_coord = [] # uma tupla com (x,y) ou um unico
        self.x_load_col = [] # valores None
        self.load_error = [] # valores None
        self.load_x_dir = [] # tuplas ou valores sozinhos
        self.load_y_dir = [] # tuplas ou valores sozinhos
        self.load_value = []

    def check_load_param(self):
        self.warnings_load = []
        for i in range(len(self.load_type_btn)):
            if self.load_type_btn[i][0].isChecked():
                self.check_param(self.warnings_load, self.load_coord_btn[i][0].text(), [int, float], 'x coordinate must be an integer or float')
                self.check_param(self.warnings_load, self.load_coord_btn[i][1].text(), [int, float], 'y coordinate must be an integer or float')
            else:
                if not (self.load_col_btn[-1][0].isChecked()) and not (self.load_col_btn[-1][1].isChecked()):
                    warning = QtWidgets.QLabel("select column")
                    self.warnings_load.append(warning)
                self.check_param(self.warnings_load, self.load_coord_btn[i][2].text(), [int, float], 'coordinate must be an integer or float')
                self.check_param(self.warnings_load, self.load_error_btn[i].text(), [int, float], 'error margin must be an integer or float')
            
            self.check_param(self.warnings_load, self.load_value_btn[i].text(), [int, float], 'load value must be an integer or float')

    def check_load_values(self):
        self.warnings_load = []
        for i in range(len(self.load_type_btn)):
            if self.load_type_btn[i][0].isChecked():
                if ast.literal_eval(self.load_coord_btn[i][0].text()) > self.lx:
                    warning = QtWidgets.QLabel("x coordinate exceeds mesh boundaries")
                    self.warnings_load.append(warning)
                if ast.literal_eval(self.load_coord_btn[i][1].text()) > self.ly:
                    warning = QtWidgets.QLabel("y coordinate exceeds mesh boundaries")
                    self.warnings_load.append(warning)
            else:
                if self.load_col_btn[-1][0].isChecked():
                    if ast.literal_eval(self.load_coord_btn[i][2].text()) > self.lx:
                        warning = QtWidgets.QLabel("coordinate exceeds mesh boundaries")
                        self.warnings_load.append(warning)
                else:
                    if ast.literal_eval(self.load_coord_btn[i][2].text()) > self.ly:
                        warning = QtWidgets.QLabel("coordinate exceeds mesh boundaries")
                        self.warnings_load.append(warning)

    def update_param_load(self):
        for i in range(len(self.load_type_btn)):
            by_coord = True if self.load_type_btn[i][0].isChecked() else False
            self.by_coord.append(by_coord)
            if by_coord:
                x = ast.literal_eval(self.load_coord_btn[i][0].text())
                y = ast.literal_eval(self.load_coord_btn[i][1].text())
                self.load_coord.append([x,y])
                self.x_load_col.append(None)
                self.load_error.append(None)
            else:
                coord = ast.literal_eval(self.load_coord_btn[i][2].text())
                self.load_coord.append(coord)
                col = 1 if self.load_col_btn[i][0].isChecked() else 2
                self.x_load_col.append(col)
                error = ast.literal_eval(self.load_error_btn[i].text())
                self.load_error.append(error)
            
            x = 1 if self.load_x_dir_btn[i][0].isChecked() else -1 if self.load_x_dir_btn[i][1].isChecked() else 0
            self.load_x_dir.append(x)

            y = 1 if self.load_y_dir_btn[i][0].isChecked() else -1 if self.load_x_dir_btn[i][1].isChecked() else 0
            self.load_y_dir.append(y)

            self.load_value = ast.literal_eval(self.load_value_btn[i].text())

    def convert_load_to_dict(self):
        # list of dicts
        self.load = []
        for i in range(len(self.load_type_btn)):
            if self.by_coord[i]:
                aux_key = ["x_coord", "y_coord", "x_direc", "y_direc", "force"]
                aux_val = [self.load_coord[i][0], self.load_coord[i][1], self.load_x_dir[i], self.load_y_dir[i], self.load_value[i]]
                #dicti = {"x_coord":self.load_coord[i][0], "y_coord":self.load_coord[i][1], "x_direc":self.load_x_dir[i], "y_direc":self.load_y_dir[i], "force":self.load_value[i]}
            else:
                aux_key = ["coord", "axis", "eps", "x_direc", "y_direc", "force"]
                aux_val = [self.load_coord[i], self.load_coord[i], self.load_error[i], self.load_x_dir[i], self.load_y_dir[i], self.load_value[i]]
                #dicti = {"coord":self.load_coord[i], "axis":self.load_coord[i], "eps":self.load_error[i], "x_direc":self.load_x_dir[i], "y_direc":self.load_y_dir[i], "force":self.load_value[i]}

            dicti = dict(zip(aux_key, aux_val))
            self.load.append(dicti)

# Node contrain param
    def add_node_constrain(self):
        pass

# Node constraint param
    def create_constraint(self):
        self.area_check = QtWidgets.QCheckBox("area")
        self.min_area_line  = QtWidgets.QLineEdit()
        self.max_area_line  = QtWidgets.QLineEdit()
 
        self.r_ratio_check = QtWidgets.QCheckBox("R ratio")
        self.min_r_ratio_line  = QtWidgets.QLineEdit()
        self.max_r_ratio_line  = QtWidgets.QLineEdit()

        self.compliance_check = QtWidgets.QCheckBox("compliance")
        self.min_compliance_line  = QtWidgets.QLineEdit()
        self.max_compliance_line  = QtWidgets.QLineEdit()
        self.freq_compliance_line  = QtWidgets.QLineEdit()

        self.local_ep_check = QtWidgets.QCheckBox("local_ep")
        self.min_local_ep_line  = QtWidgets.QLineEdit()
        self.max_local_ep_line  = QtWidgets.QLineEdit()
        self.freq_local_ep_line  = QtWidgets.QLineEdit()

        self.local_ki_check = QtWidgets.QCheckBox("local_ki")
        self.min_local_ki_line  = QtWidgets.QLineEdit()
        self.max_local_ki_line  = QtWidgets.QLineEdit()
        self.freq_local_ki_line  = QtWidgets.QLineEdit()

        self.local_r_check = QtWidgets.QCheckBox("local_r")
        self.min_local_r_line  = QtWidgets.QLineEdit()
        self.max_local_r_line  = QtWidgets.QLineEdit()
        self.freq_local_r_line  = QtWidgets.QLineEdit()

    def add_constraint_param(self, layout):
        label = QtWidgets.QLabel('---- Constraint ----')
        layout.addRow(label)
        layout.addRow(self.area_check)
        layout.addRow(QtWidgets.QLabel('min'), self.min_area_line)
        layout.addRow(QtWidgets.QLabel('max'), self.max_area_line)

        layout.addRow(self.r_ratio_check)
        layout.addRow(QtWidgets.QLabel('min'), self.min_r_ratio_line)
        layout.addRow(QtWidgets.QLabel('max'), self.max_r_ratio_line)

        layout.addRow(self.compliance_check)
        layout.addRow(QtWidgets.QLabel('min'), self.min_compliance_line)
        layout.addRow(QtWidgets.QLabel('max'), self.max_compliance_line)
        layout.addRow(QtWidgets.QLabel('freq'), self.freq_compliance_line)

        layout.addRow(self.local_ep_check)
        layout.addRow(QtWidgets.QLabel('min'), self.min_local_ep_line)
        layout.addRow(QtWidgets.QLabel('max'), self.max_local_ep_line)
        layout.addRow(QtWidgets.QLabel('freq'), self.freq_local_ep_line)

        layout.addRow(self.local_ki_check)
        layout.addRow(QtWidgets.QLabel('min'), self.min_local_ki_line)
        layout.addRow(QtWidgets.QLabel('max'), self.max_local_ki_line)
        layout.addRow(QtWidgets.QLabel('freq'), self.freq_local_ki_line)

        layout.addRow(self.local_r_check)
        layout.addRow(QtWidgets.QLabel('min'), self.min_local_r_line)
        layout.addRow(QtWidgets.QLabel('max'), self.max_local_r_line)
        layout.addRow(QtWidgets.QLabel('freq'), self.freq_local_r_line)
       
    def set_default_constraint(self):
        self.area_check.setChecked(True)
        self.min_area_line.setText('30')
 
        self.min_r_ratio_line.setDisabled(True)
        self.max_r_ratio_line.setDisabled(True)

        self.min_compliance_line.setDisabled(True)
        self.max_compliance_line.setDisabled(True)
        self.freq_compliance_line.setDisabled(True)

        self.min_local_ep_line.setDisabled(True)
        self.max_local_ep_line.setDisabled(True)
        self.freq_local_ep_line.setDisabled(True)

        self.min_local_ki_line.setDisabled(True)
        self.max_local_ki_line.setDisabled(True)
        self.freq_local_ki_line.setDisabled(True)

        self.min_local_r_line.setDisabled(True)
        self.max_local_r_line.setDisabled(True)
        self.freq_local_r_line.setDisabled(True)

        self.area_check.toggled.connect(self.min_area_line.setEnabled)
        self.area_check.toggled.connect(self.max_area_line.setEnabled)

        self.r_ratio_check.toggled.connect(self.min_r_ratio_line.setEnabled)
        self.r_ratio_check.toggled.connect(self.max_r_ratio_line.setEnabled)

        self.local_r_check.toggled.connect(self.min_local_r_line.setEnabled)
        self.local_r_check.toggled.connect(self.max_local_r_line.setEnabled)
        self.local_r_check.toggled.connect(self.freq_local_r_line.setEnabled)

        self.compliance_check.toggled.connect(self.min_compliance_line.setEnabled)
        self.compliance_check.toggled.connect(self.max_compliance_line.setEnabled)
        self.compliance_check.toggled.connect(self.freq_compliance_line.setEnabled)

        self.local_ep_check.toggled.connect(self.min_local_ep_line.setEnabled)
        self.local_ep_check.toggled.connect(self.max_local_ep_line.setEnabled)
        self.local_ep_check.toggled.connect(self.freq_local_ep_line.setEnabled)

        self.local_ki_check.toggled.connect(self.min_local_ki_line.setEnabled)
        self.local_ki_check.toggled.connect(self.max_local_ki_line.setEnabled)
        self.local_ki_check.toggled.connect(self.freq_local_ki_line.setEnabled)

         


class ParametersText():
    def __init__(self):
        self.editor = QtWidgets.QLabel()
        self.editor.setTextInteractionFlags(QtCore.Qt.LinksAccessibleByMouse)
        self.editor.setTextFormat(QtCore.Qt.RichText)
        self.editor.setTextInteractionFlags(QtCore.Qt.TextBrowserInteraction)
        self.editor.setOpenExternalLinks(True)

    def set_text(self, text):
        self.editor.setText(text)