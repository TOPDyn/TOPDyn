from PyQt5 import QtCore, QtWidgets
import numpy as np
import os
import ast
from param_bond_cond import Verification, ParamLoad, ParamNodeConstrain

class Parameters(Verification):
    def __init__(self) -> None:
        self.load = ParamLoad()
        self.node_constrain =  ParamNodeConstrain()
        self.warnings = []

    def create_btns(self):
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
        self.factor_spin = QtWidgets.QLineEdit()

        self.save_check = QtWidgets.QCheckBox("Save Data")

        self.freqrsp_check = QtWidgets.QCheckBox("Frequency Response")  
        self.freq_range_spin = QtWidgets.QLineEdit()

    def update_params(self):
        self.nelx = ast.literal_eval(self.nelx_spin.text())
        self.nely = ast.literal_eval(self.nely_spin.text())
        self.lx   = ast.literal_eval(self.lx_spin.text())
        self.ly   = ast.literal_eval(self.ly_spin.text())
        self.load.set_size(self.lx, self.ly)
        self.node_constrain.set_size(self.lx, self.ly)

        self.E = ast.literal_eval(self.E_spin.text())
        self.v = ast.literal_eval(self.v_spin.text())
        self.rho   = ast.literal_eval(self.rho_spin.text())
        
        self.alpha = ast.literal_eval(self.alpha_spin.text())
        self.beta  = ast.literal_eval(self.beta_spin.text())
        self.eta   = ast.literal_eval(self.eta_spin.text())
        
        self.factor = ast.literal_eval(self.factor_spin.text())
        self.freq = ast.literal_eval(self.freq_spin.text())

        self.freqrsp = self.freqrsp_check.isChecked()
        if self.freqrsp:
            self.freq_range = ast.literal_eval(self.freq_range_spin.text())
        else:
            self.freq_range = None
                
        self.save = self.save_check.isChecked()

    def update_default(self):
        self.nelx_spin.setText(str(self.nelx))
        self.nely_spin.setText(str(self.nely))
        self.lx_spin.setText(str(self.lx))
        self.ly_spin.setText(str(self.ly))

        self.E_spin.setText(str(self.E))
        self.v_spin.setText(str(self.v))
        self.rho_spin.setText(str(self.rho))
        
        self.alpha_spin.setText(str(self.alpha))
        self.beta_spin.setText(str(self.beta))
        self.eta_spin.setText(str(self.eta))
        
        self.freq_spin.setText(str(self.freq))
        self.factor_spin.setText(str(self.factor))

        self.freqrsp_check.setChecked(self.freqrsp) 
        
        self.save_check.setChecked(self.save)
        
    def check_params(self):
        self.warnings = []

        self.check_param(self.warnings, self.nelx_spin.text(), [int], 'Nelx must be an integer')

        self.check_param(self.warnings, self.nely_spin.text(), [int], 'Nely must be an integer')

        self.check_param(self.warnings, self.lx_spin.text(), [int, float], 'Lx must be an integer or float')

        self.check_param(self.warnings, self.ly_spin.text(), [int, float], 'Ly must be an integer or float')

        self.check_param(self.warnings, self.E_spin.text(), [int, float], "E must be an integer or float")

        self.check_param(self.warnings, self.v_spin.text(), [int, float], 'v must be an integer or float')

        self.check_param(self.warnings, self.rho_spin.text(), [int, float], 'rho must be an integer or float')

        self.check_param(self.warnings, self.alpha_spin.text(), [int, float], 'Alpha must be an integer or float')

        self.check_param(self.warnings, self.beta_spin.text(), [int, float], 'Beta must be an integer or float')

        self.check_param(self.warnings, self.eta_spin.text(), [int, float], 'Eta must be an integer or float')

        self.check_param(self.warnings, self.factor_spin.text(), [int, float], 'Factor must be an integer or float')

        self.check_param(self.warnings, self.freq_spin.text(), [int], 'Frequency must be an integer')
        
        if self.freqrsp_check.isChecked():
            self.check_param(self.warnings, self.freq_range_spin.text(), [list], 'Frequency range must be a list')

    def create_dict_param(self):
        self.node_constrain.input_to_np()
        self.load.input_to_np()

    def export_param(self):  
        folder_name = 'temp'
        directory = os.path.join(os.path.dirname(__file__), folder_name)
        os.makedirs(directory, exist_ok=True)
        try: 
            # WRITE PARAMS
            dir_file = os.path.join(directory, 'param_file.txt')
            # open file for writing
            f = open(dir_file,"w")
            # write file
            f.write(str(self.dict_param))
            # close file
            f.close()

            # WRITE LOAD MATRIX
            dir_file = os.path.join(directory, 'param_load_matrix.txt')
            np.savetxt(dir_file, self.load.np_array)

            # WRITE NODE CONSTRAIN MATRIX
            dir_file = os.path.join(directory, 'param_node_constrain_matrix.txt')
            np.savetxt(dir_file, self.node_constrain.np_array)     
        except: 
            print("Something went wrong")

class ParametersFEM2D(Parameters):
    def __init__(self):
        super().__init__()
        self.create_btns()
        self.set_default()
        #self.update_params()

    def create_btns(self):
        super().create_btns()
        self.x_coord_plot_btn = QtWidgets.QLineEdit()
        self.y_coord_plot_btn = QtWidgets.QLineEdit()
        self.x_dir_plot_btn = QtWidgets.QRadioButton("X")
        self.y_dir_plot_btn = QtWidgets.QRadioButton("Y")

    def add_white_space(self, layout, updt):
        if updt:
            layout.addRow(QtWidgets.QLabel(''))

    def add_btns(self, layout, updt=True):
        layout.addRow(QtWidgets.QLabel('Nelx'))
        layout.addRow(self.nelx_spin)
        #self.add_white_space(layout, updt)
        
        layout.addRow(QtWidgets.QLabel('Nely'))
        layout.addRow(self.nely_spin)
        
        layout.addRow(QtWidgets.QLabel('Lx'))
        layout.addRow(self.lx_spin)

        layout.addRow(QtWidgets.QLabel('Ly'))
        layout.addRow(self.ly_spin)

        layout.addRow(QtWidgets.QLabel('E'))
        layout.addRow(self.E_spin)

        layout.addRow(QtWidgets.QLabel('v'))
        layout.addRow(self.v_spin)

        layout.addRow(QtWidgets.QLabel('rho'))
        layout.addRow(self.rho_spin)

        layout.addRow(QtWidgets.QLabel('Alpha'))
        layout.addRow(self.alpha_spin)

        layout.addRow(QtWidgets.QLabel('Beta'))
        layout.addRow(self.beta_spin)

        layout.addRow(QtWidgets.QLabel('Eta'))
        layout.addRow(self.eta_spin)

        layout.addRow(QtWidgets.QLabel('Factor'))
        layout.addRow(self.factor_spin)
        
        layout.addRow(QtWidgets.QLabel('Frequency'))
        layout.addRow(self.freq_spin)

        layout.addRow(self.save_check)
        layout.addRow(self.freqrsp_check)
        layout.addRow(QtWidgets.QLabel('Frequency Range'))
        layout.addRow(self.freq_range_spin)

        layout.addRow(QtWidgets.QLabel('Node to plot:'))
        layout.addRow(QtWidgets.QLabel('X-coord'), self.x_coord_plot_btn)
        layout.addRow(QtWidgets.QLabel('Y-coord'), self.y_coord_plot_btn)

        directions = QtWidgets.QButtonGroup(layout)
        directions.addButton(self.x_dir_plot_btn)
        directions.addButton(self.y_dir_plot_btn)
        layout.addRow(self.x_dir_plot_btn)
        layout.addRow(self.y_dir_plot_btn)     

    def toggled_fem2d(self):
        self.freqrsp_check.toggled.connect(self.freq_range_spin.setEnabled)
        self.freqrsp_check.toggled.connect(self.x_coord_plot_btn.setEnabled)
        self.freqrsp_check.toggled.connect(self.y_coord_plot_btn.setEnabled)
        self.freqrsp_check.toggled.connect(self.x_dir_plot_btn.setEnabled)
        self.freqrsp_check.toggled.connect(self.y_dir_plot_btn.setEnabled)

    def set_default(self):       
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

        self.freq_range_spin.setText("[0, 25, 5]")    
        self.freq_range_spin.setDisabled(True)
        self.x_coord_plot_btn.setDisabled(True)
        self.y_coord_plot_btn.setDisabled(True)
        self.x_dir_plot_btn.setDisabled(True)
        self.y_dir_plot_btn.setDisabled(True)

        self.toggled_fem2d()

    def update_params(self): 
        super().update_params()
        if self.freqrsp:
            self.x_coord_plot = ast.literal_eval(self.x_coord_plot_btn.text())
            self.y_coord_plot = ast.literal_eval(self.y_coord_plot_btn.text())
            self.x_dir_plot = 1 if self.x_dir_plot_btn.isChecked() else 0
            self.y_dir_plot = 1 if self.y_dir_plot_btn.isChecked() else 0

    def update_default(self):
        super().update_default()
        if self.freqrsp:
            self.freq_range_spin.setText(str(self.freq_range))
            self.x_coord_plot_btn.setText(str(self.x_coord_plot))
            self.y_coord_plot_btn.setText(str(self.y_coord_plot))
            x = True if self.x_dir_plot == 1 else False
            self.x_dir_plot_btn.setChecked(x)
            self.y_dir_plot_btn.setChecked(not x)
        else:
            self.freq_range_spin.setText('[5,500,5]')
            self.x_coord_plot_btn.setText('0')
            self.y_coord_plot_btn.setText('0')
            self.x_dir_plot_btn.setChecked(True)
            self.freq_range_spin.setDisabled(True)
            self.x_coord_plot_btn.setDisabled(True)
            self.y_coord_plot_btn.setDisabled(True)
            self.x_dir_plot_btn.setDisabled(True)
            self.y_dir_plot_btn.setDisabled(True)

        self.toggled_fem2d()

    def set_node_plot(self):
        if self.freqrsp:
            self.node_plot = [self.x_coord_plot, self.y_coord_plot, self.x_dir_plot, self.y_dir_plot]
        else:
            self.node_plot = None

    def create_dict_param(self):
        super().create_dict_param()
        self.set_node_plot()

        self.dict_param = {"nelx":self.nelx, "nely":self.nely, "lx":self.lx, "ly":self.ly, "E":self.E, "v":self.v,
                            "alpha":self.alpha, "beta":self.beta, "eta":self.eta, "freq":self.freq, "rho":self.rho,
                            "freqrsp":self.freqrsp, "freq_range":self.freq_range, "factor":self.factor,
                            "node_plot":self.node_plot, "save":self.save, "mesh_file":None}

    def check_params(self):
        super().check_params()
        if self.freqrsp_check.isChecked():
            self.check_param(self.warnings, self.x_coord_plot_btn.text(), [int, float], 'Node to plot: node coordinate must be an integer or float')
            self.check_param(self.warnings, self.y_coord_plot_btn.text(), [int, float], 'Node to plot: node coordinate must be an integer or float')
        if len(self.warnings) == 0:
            self.check_node_plot()

    def check_node_plot(self):
        if self.freqrsp_check.isChecked():
            if ast.literal_eval(self.x_coord_plot_btn.text()) > ast.literal_eval(self.lx_spin.text()):
                warning = QtWidgets.QLabel("Node to plot: x coordinate exceeds mesh boundaries")
                self.warnings.append(warning)
            if ast.literal_eval(self.y_coord_plot_btn.text()) > ast.literal_eval(self.ly_spin.text()):
                warning = QtWidgets.QLabel("Node to plot: y coordinate exceeds mesh boundaries")
                self.warnings.append(warning)

class ParametersOpt(Parameters):
    def __init__(self):
        super().__init__()

        self.create_btns()
        self.set_default()
        #self.update_params()
       
    def create_dict_param(self):
        super().create_dict_param()
        self.dict_param = {"nelx":self.nelx, "nely":self.nely, "lx":self.lx, "ly":self.ly, "E":self.E, "v":self.v, "rho":self.rho,
                            "alpha_par":self.alpha, "beta_par":self.beta, "eta_par":self.eta, "factor":self.factor, "freq":self.freq, 
                            "freqrsp":self.freqrsp, "freq_range":self.freq_range, "save":self.save, "mesh_file":None, "mma":self.mma,
                            "fac_ratio":self.fac_ratio, "x_min_m":self.x_min_m, "x_min_k":self.x_min_k, "penal_k":self.penal_k,
                            "penal_m":self.penal_m, "constr_func":self.constr_func, "constr_values":self.constr_values,
                            "passive_coord":self.passive_coord, "modes":self.modes,"const_func":self.const_func,"n1":self.n1,
                            "func_name":self.func_name,"func_name2":self.func_name2,"freq2":self.freq2,"alpha_plot":self.alpha_plot,
                            "beta_plot":self.beta_plot,"eta_plot":self.eta_plot,"max_iter":self.max_iter,
                            "save":self.const_func,"save":self.save,"dens_filter":self.dens_filter,"mesh_deform":self.mesh_deform}

    def create_btns(self):
        super().create_btns()
        self.mma_radio = QtWidgets.QRadioButton("MMA")
        self.gcmma_radio = QtWidgets.QRadioButton("GMMA")
        self.fac_ratio_spin = QtWidgets.QLineEdit()

        self.x_min_m_spin = QtWidgets.QLineEdit()
        self.x_min_k_spin = QtWidgets.QLineEdit()
        self.penal_k_spin   = QtWidgets.QLineEdit() #p_par
        self.penal_m_spin   = QtWidgets.QLineEdit() #q_par

        self.passive_coord_spin = QtWidgets.QLineEdit()    
        self.modes_spin = QtWidgets.QLineEdit()
        self.const_func_spin = QtWidgets.QLineEdit()
        self.n1_spin = QtWidgets.QLineEdit()

        self.func_name_box = QtWidgets.QComboBox()
        self.func_name_box.addItem("Compliance")
        self.func_name_box.addItem("Input power")
        self.func_name_box.addItem("Elastic Potential Energy")
        self.func_name_box.addItem("Kinetic Energy")
        self.func_name_box.addItem("Strain-to-kinetic energy ratio")
        self.func_name_box.addItem("Local elastic potential energy")
        self.func_name_box.addItem("Local kinetic energy")
        self.func_name_box.addItem("Local strain-to-kinetic energy ratio")

        self.func_name2_box = QtWidgets.QComboBox()
        self.func_name2_box.addItem("None")
        self.func_name2_box.addItem("Compliance")
        self.func_name2_box.addItem("Input power")
        self.func_name2_box.addItem("Elastic Potential Energy")
        self.func_name2_box.addItem("Kinetic Energy")
        self.func_name2_box.addItem("Strain-to-kinetic energy ratio")
        self.func_name2_box.addItem("Local elastic potential energy")
        self.func_name2_box.addItem("Local kinetic energy")
        self.func_name2_box.addItem("Local strain-to-kinetic energy ratio")

        self.freq2_spin = QtWidgets.QLineEdit()

        self.alpha_plot_spin = QtWidgets.QLineEdit()
        self.beta_plot_spin  = QtWidgets.QLineEdit()
        self.eta_plot_spin   = QtWidgets.QLineEdit()

        self.max_iter_spin = QtWidgets.QLineEdit() 
        self.dens_filter_check = QtWidgets.QCheckBox("Density Filter")      
        self.mesh_deform_check = QtWidgets.QCheckBox("Plot Deformed Mesh")
    
    def add_btns(self, layout):
        layout.addRow(QtWidgets.QLabel('Optimization Method'))
        layout.addRow(self.mma_radio)
        layout.addRow(self.gcmma_radio)

        layout.addRow(QtWidgets.QLabel('Nelx'))
        layout.addRow(self.nelx_spin)
        
        layout.addRow(QtWidgets.QLabel('Nely'))
        layout.addRow(self.nely_spin)

        layout.addRow(QtWidgets.QLabel('Lx'))
        layout.addRow(self.lx_spin)

        layout.addRow(QtWidgets.QLabel('Ly'))
        layout.addRow(self.ly_spin)

        layout.addRow(QtWidgets.QLabel('E'))
        layout.addRow(self.E_spin)

        layout.addRow(QtWidgets.QLabel('v'))
        layout.addRow(self.v_spin)

        layout.addRow(QtWidgets.QLabel('rho'))
        layout.addRow(self.rho_spin)

        layout.addRow(QtWidgets.QLabel('Ratio Factor'))
        layout.addRow(self.fac_ratio_spin)

        layout.addRow(QtWidgets.QLabel('x_min_mass'))
        layout.addRow(self.x_min_m_spin)

        layout.addRow(QtWidgets.QLabel('x_min_stif'))
        layout.addRow(self.x_min_k_spin)

        layout.addRow(QtWidgets.QLabel('Penal. Stiffness'))
        layout.addRow(self.penal_k_spin)

        layout.addRow(QtWidgets.QLabel('Penal. Mass'))
        layout.addRow(self.penal_m_spin)        

        layout.addRow(QtWidgets.QLabel('Alpha'))
        layout.addRow(self.alpha_spin)

        layout.addRow(QtWidgets.QLabel('Beta'))
        layout.addRow(self.beta_spin)

        layout.addRow(QtWidgets.QLabel('Eta'))
        layout.addRow(self.eta_spin)

        layout.addRow(QtWidgets.QLabel('Passive Coordinates'))
        layout.addRow(self.passive_coord_spin)

        layout.addRow(QtWidgets.QLabel('Modes'))
        layout.addRow(self.modes_spin)

        layout.addRow(QtWidgets.QLabel('Constant Function'))
        layout.addRow(self.const_func_spin)

        layout.addRow(QtWidgets.QLabel('Objective Function Weight'))
        layout.addRow(self.n1_spin)

        layout.addRow(QtWidgets.QLabel('Frequency'))
        layout.addRow(self.freq_spin)

        layout.addRow(QtWidgets.QLabel('Objective Function'))
        layout.addRow(self.func_name_box)

        layout.addRow(QtWidgets.QLabel('Multiobjective Function'))
        layout.addRow(self.func_name2_box)
        layout.addRow(QtWidgets.QLabel('Multiobjective Frequency'))
        layout.addRow(self.freq2_spin)

        layout.addRow(QtWidgets.QLabel('Max. Iterations'))
        layout.addRow(self.max_iter_spin)

        layout.addRow(self.save_check)
        layout.addRow(self.dens_filter_check)

        layout.addRow(self.mesh_deform_check)
        layout.addRow(QtWidgets.QLabel('Factor'))
        layout.addRow(self.factor_spin)
      
        layout.addRow(self.freqrsp_check)
        layout.addRow(QtWidgets.QLabel('Frequency Range'))
        layout.addRow(self.freq_range_spin)

        layout.addRow(QtWidgets.QLabel('Alpha Plot'))
        layout.addRow(self.alpha_plot_spin)

        layout.addRow(QtWidgets.QLabel('Beta Plot'))
        layout.addRow(self.beta_plot_spin)

        layout.addRow(QtWidgets.QLabel('Eta Plot'))
        layout.addRow(self.eta_plot_spin)

    def update_params(self): 
        super().update_params()
        self.mma = True if self.mma_radio.isChecked() else False
        self.fac_ratio = ast.literal_eval(self.fac_ratio_spin.text())

        self.x_min_m = ast.literal_eval(self.x_min_m_spin.text())
        self.x_min_k = ast.literal_eval(self.x_min_k_spin.text())
        self.penal_k = ast.literal_eval(self.penal_k_spin.text()) #p_par
        self.penal_m = ast.literal_eval(self.penal_m_spin.text()) #q_par

        if self.passive_coord_spin.text():
            self.passive_coord = ast.literal_eval(self.passive_coord_spin.text())
        else: 
            self.passive_coord = None
        if self.modes_spin.text():
            self.modes = ast.literal_eval(self.modes_spin.text())
        else:
            self.modes = None

        self.const_func = ast.literal_eval(self.const_func_spin.text())
        self.n1 = ast.literal_eval(self.n1_spin.text())
        self.func_name = self.get_func_name(self.func_name_box)
        
        self.func_name2 = self.get_func_name(self.func_name2_box)
        if self.func_name2 is not None:
            self.freq2 = ast.literal_eval(self.freq2_spin.text())
        else:
            self.freq2 = None
        
        if self.freqrsp:
            self.alpha_plot = ast.literal_eval(self.alpha_plot_spin.text())
            self.beta_plot  = ast.literal_eval(self.beta_plot_spin.text())
            self.eta_plot   = ast.literal_eval(self.eta_plot_spin.text())
        else:
            self.alpha_plot = None
            self.beta_plot  = None
            self.eta_plot   = None

        self.max_iter = ast.literal_eval(self.max_iter_spin.text())
        self.save = self.save_check.isChecked()

        self.dens_filter  = self.dens_filter_check.isChecked()
        self.mesh_deform  = self.mesh_deform_check.isChecked()

    def get_func_name(self, btn):       
        if str(btn.currentText()) == "Compliance":
            func_name = "compliance"
        elif str(btn.currentText()) == "Input power":
            func_name = "input_power"
        elif str(btn.currentText()) == "Kinetic Energy":
            func_name = "kinetic_energy"
        elif str(btn.currentText()) == "Elastic Potential Energy":
            func_name = "elastic_potential_energy"
        elif str(btn.currentText()) == "Strain-to-kinetic energy ratio":
            func_name = "r_ratio"
        elif str(btn.currentText()) == "Local elastic potential energy":
            func_name = "local_ep"
        elif str(btn.currentText()) == "Local kinetic energy":
            func_name = "local_ki"
        elif str(btn.currentText()) == "Local strain-to-kinetic energy ratio":
            func_name = "local_r"
        else:
            func_name = None
        return func_name

    def get_box_name(self, func_name):
        if func_name == "compliance":
            box = "Compliance"
        elif func_name == "input_power":
            box = "Input power"
        elif func_name == "kinetic_energy":
            box = "Kinetic Energy"
        elif func_name == "elastic_potential_energy":
            box = "Elastic Potential Energy"
        elif func_name == "r_ratio":
            box = "Strain-to-kinetic energy ratio"
        elif func_name == "local_ep":
            box = "Local elastic potential energy"
        elif func_name == "local_ki":
            box = "Local kinetic energy"
        elif func_name == "local_r":
            box = "Local strain-to-kinetic energy ratio"
        else:
            box = "None"
        return box

    def toggled_opt(self):
        self.mesh_deform_check.toggled.connect(self.factor_spin.setEnabled)   
        self.func_name2_box.activated.connect(self.freq2_spin.setEnabled) 
        self.freqrsp_check.toggled.connect(self.freq_range_spin.setEnabled)
        self.freqrsp_check.toggled.connect(self.alpha_plot_spin.setEnabled)
        self.freqrsp_check.toggled.connect(self.beta_plot_spin.setEnabled)
        self.freqrsp_check.toggled.connect(self.eta_plot_spin.setEnabled)

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
        self.dens_filter_check.setChecked(True)
        self.toggled_opt()

    def update_default(self):
        super().update_default()
        self.mma_radio.setChecked(self.mma)
        self.gcmma_radio.setChecked(not self.mma)
        self.fac_ratio_spin.setText(str(self.fac_ratio))

        self.x_min_m_spin.setText(str(self.x_min_m))
        self.x_min_k_spin.setText(str(self.x_min_k))
        self.penal_k_spin.setText(str(self.penal_k))
        self.penal_m_spin.setText(str(self.penal_m))

        if self.passive_coord is not None:
            self.passive_coord_spin.setText(str(self.passive_coord))
        if self.modes is not None:
            self.modes_spin.setText(str(self.modes))
        self.const_func_spin.setText(str(self.const_func))
        self.n1_spin.setText(str(self.n1))    
        
        self.max_iter_spin.setText(str(self.max_iter))

        name = self.get_box_name(self.func_name)
        index = self.func_name_box.findText(name, QtCore.Qt.MatchFixedString)
        if index > -1:
            self.func_name_box.setCurrentIndex(index) 

        if self.freqrsp:
            self.freq_range_spin.setText(str(self.freq_range))
            self.alpha_plot_spin.setText(str(self.alpha_plot))
            self.beta_plot_spin.setText(str(self.beta_plot))
            self.eta_plot_spin.setText(str(self.eta_plot))
            self.freqrsp_check.setChecked(True) 
            self.freq_range_spin.setEnabled(True)
            self.alpha_plot_spin.setEnabled(True)
            self.beta_plot_spin.setEnabled(True)
            self.eta_plot_spin.setEnabled(True)
        else:
            self.freq_range_spin.setText('[5,500,5]')
            self.alpha_plot_spin.setText('0')
            self.beta_plot_spin.setText('1e-6')
            self.eta_plot_spin.setText('0')
            self.freq_range_spin.setDisabled(True)
            self.alpha_plot_spin.setDisabled(True)
            self.beta_plot_spin.setDisabled(True)
            self.eta_plot_spin.setDisabled(True)

        if self.mesh_deform:
            self.mesh_deform_check.setChecked(True)
            self.factor_spin.setEnabled(True)
        else:
            self.factor_spin.setDisabled(True)

        if self.dens_filter:
            self.dens_filter_check.setChecked(True)

        if self.func_name2 != None:
            self.freq2_spin.setText(str(self.freq2))
            self.freq2_spin.setEnabled(True)
            name = self.get_box_name(self.func_name2)
            index = self.func_name_box.findText(name, QtCore.Qt.MatchFixedString)
            if index > -1:
                self.func_name2_box.setCurrentIndex(index+1) 
        else:
            self.freq2_spin.setDisabled(True)
            self.freq2_spin.setText('0')

        self.toggled_opt()

    def check_params(self):
        super().check_params()

        if str(self.func_name2_box.currentText()) is not None:
            self.check_param(self.warnings, self.freq2_spin.text(), [int], 'Multiobjective frequency function must be an integer')

        self.check_param(self.warnings, self.x_min_m_spin.text(), [int, float], 'x_min_mass must be an integer or float')

        self.check_param(self.warnings, self.x_min_k_spin.text(), [int, float], 'x_min_stif must be an integer or float')

        self.check_param(self.warnings, self.penal_k_spin.text(), [int], 'Penal. stiffness must be an integer')

        self.check_param(self.warnings, self.penal_m_spin.text(), [int], 'Penal. mass must be an integer')

        self.check_param(self.warnings, self.const_func_spin.text(), [int, float], 'Constant Function must be an integer or float')

        self.check_param(self.warnings, self.n1_spin.text(), [int, float], 'Objective function weight must be an integer or float')

        self.check_param(self.warnings, self.fac_ratio_spin.text(), [int, float], 'Ratio Factor must be an integer or float')

        if self.modes_spin.text():
            self.check_param(self.warnings, self.modes_spin.text(), [int], 'Modes must be an integer.')
            #TODO: PRECISA VERIFICAR SE O NUMERO DE MODES É VALIDO.
  
        if self.passive_coord_spin.text():
            aux_war = 'Passive Coordinates must be a tuple of tuples'
            try:
                if isinstance(ast.literal_eval(self.passive_coord_spin.text()), tuple):
                    for item in ast.literal_eval(self.passive_coord_spin.text()):
                        if not isinstance(item, tuple):
                            self.warnings.append(QtWidgets.QLabel(aux_war))
                            break
            except:
                warning = QtWidgets.QLabel(aux_war)
                self.warnings.append(warning)

            if len(self.warnings) == 0:
                aux_passive = ast.literal_eval(self.passive_coord_spin.text())

                self.check_passive_coord(aux_passive, ast.literal_eval(self.lx_spin.text()), 'x')
                self.check_passive_coord(aux_passive, ast.literal_eval(self.ly_spin.text()), 'y')
            
        self.check_param(self.warnings, self.max_iter_spin.text(), [int, float], 'Max. Iterations must be an integer or float')

        if self.freqrsp_check.isChecked():
            self.check_param(self.warnings, self.alpha_plot_spin.text(), [int, float], 'Alpha Plot must be an integer or float')
            self.check_param(self.warnings, self.beta_plot_spin.text(), [int, float], 'Beta Plot must be an integer or float')
            self.check_param(self.warnings, self.eta_plot_spin.text(), [int, float], 'Eta Plot must be an integer or float')

    def check_passive_coord(self, passive_coord, coord, axis):
        if axis == 'x':
            war1 = "Passive Coordinates on X-Axis exceeds Lx."
            war2 = "Passive Coordinates on X-Axis are invalid:"
        else:
            war1 = "Passive Coordinates on Y-Axis exceeds Ly."
            war2 = "Passive Coordinates on Y-Axis are invalid:"
          
        if passive_coord[0][0] < passive_coord[0][1]:
            if passive_coord[0][0] > coord or passive_coord[0][1] > coord:
                self.warnings.append(QtWidgets.QLabel(war1))
        else:
            self.warnings.append(QtWidgets.QLabel(war2 + " The first value must be less than the second value."))

# Constraint
    def create_constraint(self):
        self.area_check = QtWidgets.QCheckBox("Area")
        self.min_area_line  = QtWidgets.QLineEdit()
        self.max_area_line  = QtWidgets.QLineEdit()
 
        self.r_ratio_check = QtWidgets.QCheckBox("Strain-to-kinetic energy ratio")
        self.min_r_ratio_line  = QtWidgets.QLineEdit()
        self.max_r_ratio_line  = QtWidgets.QLineEdit()

        self.compliance_check = QtWidgets.QCheckBox("Compliance")
        self.min_compliance_line  = QtWidgets.QLineEdit()
        self.max_compliance_line  = QtWidgets.QLineEdit()
        self.freq_compliance_line  = QtWidgets.QLineEdit()

        self.local_ep_check = QtWidgets.QCheckBox("Local elastic potential energy")
        self.min_local_ep_line  = QtWidgets.QLineEdit()
        self.max_local_ep_line  = QtWidgets.QLineEdit()
        self.freq_local_ep_line  = QtWidgets.QLineEdit()

        self.local_ki_check = QtWidgets.QCheckBox("Local kinetic energy")
        self.min_local_ki_line  = QtWidgets.QLineEdit()
        self.max_local_ki_line  = QtWidgets.QLineEdit()
        self.freq_local_ki_line  = QtWidgets.QLineEdit()

        self.local_r_check = QtWidgets.QCheckBox("Local strain-to-kinetic energy ratio")
        self.min_local_r_line  = QtWidgets.QLineEdit()
        self.max_local_r_line  = QtWidgets.QLineEdit()
        self.freq_local_r_line  = QtWidgets.QLineEdit()

    def add_constraint_param(self, layout):
        label = QtWidgets.QLabel('---- Constraint ----')
        label.setAlignment(QtCore.Qt.AlignCenter)
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

    def check_min_max_constraint(self, maxi, mini, func):
        if float(maxi) <= float(mini):
            self.warnings_constraint.append(QtWidgets.QLabel(func + " - The maximum value must be greater than the minimum value."))

    def check_constraint_values(self):
        self.warnings_constraint = []
        if self.area_check.isChecked() and self.max_area_line.text():
            self.check_min_max_constraint(self.max_area_line.text(), self.min_area_line.text(), 'Area')
                   
        if self.r_ratio_check.isChecked() and self.max_r_ratio_line.text():
            self.check_min_max_constraint(self.max_r_ratio_line.text(), self.min_r_ratio_line.text(), 'Strain-to-kinetic energy ratio')

        if self.compliance_check.isChecked() and self.max_compliance_line.text():
            self.check_min_max_constraint(self.max_compliance_line.text(), self.min_compliance_line.text(), 'Compliance')

        if self.local_ep_check.isChecked() and self.max_local_ep_line.text():
            self.check_min_max_constraint(self.max_local_ep_line.text(), self.min_local_ep_line.text(), 'Local elastic potential energy')

        if self.local_ki_check.isChecked() and self.max_local_ki_line.text():
            self.check_min_max_constraint(self.max_local_ki_line.text(), self.min_local_ki_line.text(), 'Local kinetic energy')

        if self.local_r_check.isChecked() and self.max_local_r_line.text():
            self.check_min_max_constraint(self.max_local_r_line.text(), self.min_local_r_line.text(), 'Local strain-to-kinetic energy ratio')

    def check_constraint(self):
        self.warnings_constraint = []

        if self.area_check.isChecked():
            self.check_param(self.warnings_constraint, self.min_area_line.text(), [int, float], 'Area - min value must be an integer or float')
            if self.max_area_line.text():
                self.check_param(self.warnings_constraint, self.max_area_line.text(), [int, float], 'Area - max value must be an integer or float')
                   
        if self.r_ratio_check.isChecked():
            self.check_param(self.warnings_constraint, self.min_r_ratio_line.text(), [int, float], 'Strain-to-kinetic energy ratio - min value must be an integer or float')
            if self.max_r_ratio_line.text():
                self.check_param(self.warnings_constraint, self.max_r_ratio_line.text(), [int, float], 'Strain-to-kinetic energy ratio - max value must be an integer or float')

        if self.compliance_check.isChecked():
            self.check_param(self.warnings_constraint, self.min_compliance_line.text(), [int, float], 'Compliance - minimum value must be an integer or float')
            self.check_param(self.warnings_constraint, self.freq_compliance_line.text(), [int, float], 'Compliance - frequency must be an integer or float')
            if self.max_compliance_line.text():
                self.check_param(self.warnings_constraint, self.max_compliance_line.text(), [int, float], 'Compliance - maximum value must be an integer or float')

        if self.local_ep_check.isChecked():
            self.check_param(self.warnings_constraint, self.min_local_ep_line.text(), [int, float], 'Local elastic potential energy - minimum value must be an integer or float')
            self.check_param(self.warnings_constraint, self.freq_local_ep_line.text(), [int, float], 'Local elastic potential energy - frequency must be an integer or float')
            if self.max_local_ep_line.text():
                self.check_param(self.warnings_constraint, self.max_local_ep_line.text(), [int, float], 'Local elastic potential energy - maximum value must be an integer or float')

        if self.local_ki_check.isChecked():
            self.check_param(self.warnings_constraint, self.min_local_ki_line.text(), [int, float], 'Local kinetic energy - minimum value must be an integer or float')
            self.check_param(self.warnings_constraint, self.freq_local_ki_line.text(), [int, float], 'Local kinetic energy - frequency must be an integer or float')
            if self.max_local_ki_line.text():
                self.check_param(self.warnings_constraint, self.max_local_ki_line.text(), [int, float], 'Local kinetic energy - maximum value must be an integer or float')

        if self.local_r_check.isChecked():
            self.check_param(self.warnings_constraint, self.min_local_r_line.text(), [int, float], 'Local strain-to-kinetic energy ratio - minimum value must be an integer or float')
            self.check_param(self.warnings_constraint, self.freq_local_r_line.text(), [int, float], 'Local strain-to-kinetic energy ratio - frequency must be an integer or float')
            if self.max_local_r_line.text():
                self.check_param(self.warnings_constraint, self.max_local_r_line.text(), [int, float], 'Local strain-to-kinetic energy ratio - maximum value must be an integer or float')

    def update_constraint(self):
        self.area = self.area_check.isChecked()
        if self.area:
            self.min_area = ast.literal_eval(self.min_area_line.text())
            self.max_area = ast.literal_eval(self.max_area_line.text()) if self.max_area_line.text() else None
        
        self.r_ratio = self.r_ratio_check.isChecked()
        if self.r_ratio:
            self.min_r_ratio = ast.literal_eval(self.min_r_ratio_line.text())
            self.max_r_ratio = ast.literal_eval(self.max_r_ratio_line.text()) if self.max_r_ratio_line.text() else None

        self.compliance = self.compliance_check.isChecked()
        if self.compliance:
            self.min_compliance = ast.literal_eval(self.min_compliance_line.text())
            self.max_compliance = ast.literal_eval(self.max_compliance_line.text()) if self.max_compliance_line.text() else None
            self.freq_compliance = ast.literal_eval(self.freq_compliance_line.text())

        self.local_ep = self.local_ep_check.isChecked()
        if self.local_ep:
            self.min_local_ep = ast.literal_eval(self.min_local_ep_line.text())
            self.max_local_ep = ast.literal_eval(self.max_local_ep_line.text()) if self.max_local_ep_line.text() else None
            self.freq_local_ep = ast.literal_eval(self.freq_local_ep_line.text()) if self.freq_local_ep_line.text() else None

        self.local_ki = self.local_ki_check.isChecked()
        if self.local_ki:
            self.min_local_ki = ast.literal_eval(self.min_local_ki_line.text())
            self.max_local_ki = ast.literal_eval(self.max_local_ki_line.text()) if self.max_local_ki_line.text() else None
            self.freq_local_ki = ast.literal_eval(self.freq_local_ki_line.text())

        self.local_r = self.local_r_check.isChecked()
        if self.local_r:
            self.min_local_r = ast.literal_eval(self.min_local_r_line.text())
            self.max_local_r = ast.literal_eval(self.max_local_r_line.text()) if self.max_local_r_line.text() else None
            self.freq_local_r = ast.literal_eval(self.freq_local_r_line.text())

    def upd_constr_func(self, func, mini, maxi, freq=None):
        self.constr_func.append(func)
        if freq is not None:
            self.constr_values.append([mini, freq])
        else:
            self.constr_values.append(mini)
        if maxi is not None:
            self.constr_func.append(func)
            if freq is not None:
                self.constr_values[-1][0] *= -1
                self.constr_values.append([maxi, freq])
            else:
                self.constr_values[-1] *= -1
                self.constr_values.append(maxi)

    def constraint_to_list(self):
        self.constr_func = []
        self.constr_values = []

        if self.area:
            self.upd_constr_func("area", self.min_area, self.max_area)
         
        if self.r_ratio:
            self.upd_constr_func("r_ratio", self.min_r_ratio, self.max_r_ratio)

        if self.compliance:
            self.upd_constr_func("compliance", self.min_compliance, self.max_compliance, self.freq_compliance)

        if self.local_ep:
            self.upd_constr_func("local_ep", self.min_local_ep, self.max_local_ep, self.freq_local_ep)

        if self.local_ki:
            self.upd_constr_func("local_ki", self.min_local_ki, self.max_local_ki, self.freq_local_ki)

        if self.local_r:
            self.upd_constr_func("local_r", self.min_local_r, self.max_local_r, self.freq_local_r)