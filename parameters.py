from PyQt5 import QtCore, QtWidgets
import os
import ast

class Parameters():
    def __init__(self) -> None:
        # btn load
        self.load_type_btn = []
        self.load_coord_btn = []
        self.load_col_btn = []
        self.load_error_btn = []
        self.load_x_dir_btn = []
        self.load_y_dir_btn = []
        self.load_value_btn = []

        # values load
        self.load_by_coord = []
        self.load_coord = []
        self.x_load_col = []
        self.load_error = []
        self.load_x_dir = []
        self.load_y_dir = []
        self.load_value = []

        # nodes displacements constrained
        self.node_constrain_type_btn = []
        self.node_constrain_coord_btn = []
        self.node_constrain_col_btn = []
        self.node_constrain_error_btn = []
        self.node_constrain_x_dir_btn = []
        self.node_constrain_y_dir_btn = []

        # values nodes displacements constrained
        self.node_constrain_by_coord = []
        self.node_constrain_coord = []
        self.x_node_constrain_col = []
        self.node_constrain_error = []
        self.node_constrain_x_dir = []
        self.node_constrain_y_dir = []

        self.warnings = []
        self.warnings_load = []
        self.warnings_node_constrain = []

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
        if self.freqrsp:
            self.freq_range = ast.literal_eval(self.freq_range_spin.text())
        else:
            self.freq_range = None
                
        self.save = self.save_check.checkState()

    def check_param(self, warnings, input_val, types, war):
        try:
            if input_val:
                for type in types:
                    isinstance(ast.literal_eval(input_val), type)
            else:
                warning = QtWidgets.QLabel(war)
                warnings.append(warning)
        except:
            warning = QtWidgets.QLabel(war)
            warnings.append(warning)
        
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

# Load
    def create_load_btn(self):
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

    def add_load_btn(self, layout):
        load_type = QtWidgets.QButtonGroup(layout)
        load_label = QtWidgets.QLabel('---- Load ----')
        load_label.setAlignment(QtCore.Qt.AlignCenter)
        layout.addRow(load_label)
        
        load_type.addButton(self.load_type_btn[-1][0])
        layout.addRow(self.load_type_btn[-1][0])
       
        layout.addRow(QtWidgets.QLabel('X-coord'), self.load_coord_btn[-1][0])
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

        layout.addRow(QtWidgets.QLabel('Load value'), self.load_value_btn[-1])

    def toggled_load(self):
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

        self.toggled_load()

    def update_default_load(self, ind):
        self.load_type_btn[ind][0].setChecked(self.load_by_coord[ind])
        self.load_type_btn[ind][1].setChecked(not self.load_by_coord[ind])

        if self.load_by_coord[ind]:
            self.load_coord_btn[ind][0].setText(str(self.load_coord[ind][0]))
            self.load_coord_btn[ind][1].setText(str(self.load_coord[ind][1]))
        
            self.load_coord_btn[ind][2].setDisabled(True)
            self.load_col_btn[ind][0].setDisabled(True)
            self.load_col_btn[ind][1].setDisabled(True)
            self.load_error_btn[ind].setDisabled(True)
        else:
            self.load_coord_btn[ind][0].setDisabled(True)
            self.load_coord_btn[ind][1].setDisabled(True)

            self.load_coord_btn[ind][2].setText(str(self.load_coord[ind]))
            col = True if self.x_load_col[ind] else False
            self.load_col_btn[ind][0].setChecked(col)
            self.load_col_btn[ind][1].setChecked(not col)
            self.load_error_btn[ind].setText(str(self.load_error[ind]))

        ind_x = 0 if self.load_x_dir[ind] == 1 else 1 if self.load_x_dir[ind] == -1 else 2
        self.load_x_dir_btn[ind][ind_x].setChecked(True)

        ind_y = 0 if self.load_y_dir[ind] == 1 else 1 if self.load_y_dir[ind] == -1 else 2
        self.load_y_dir_btn[ind][ind_y].setChecked(True)
    
        self.load_value_btn[ind].setText(str(self.load_value[ind]))

        self.toggled_load()

    def rewrite_load(self, layout):
        self.reset_load_list(vals=False)
        for ind in range(len(self.load_by_coord)):
            self.create_load_btn()
            self.update_default_load(ind)
            self.add_load_btn(layout)
        
    def reset_load_list(self, btn=True, vals=True):
        if btn:
            self.load_type_btn = []
            self.load_coord_btn = []
            self.load_col_btn = []
            self.load_error_btn = []
            self.load_x_dir_btn = []
            self.load_y_dir_btn = []
            self.load_value_btn = []
        if vals:
            self.load_by_coord = []
            self.load_coord = []
            self.x_load_col = []
            self.load_error = []
            self.load_x_dir = []
            self.load_y_dir = []
            self.load_value = []

    def check_load_btn(self):
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

    def update_load(self):
        for i in range(len(self.load_type_btn)):
            load_by_coord = True if self.load_type_btn[i][0].isChecked() else False
            self.load_by_coord.append(load_by_coord)
            if load_by_coord:
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

            y = 1 if self.load_y_dir_btn[i][0].isChecked() else -1 if self.load_y_dir_btn[i][1].isChecked() else 0
            self.load_y_dir.append(y)

            self.load_value.append(ast.literal_eval(self.load_value_btn[i].text()))

    def convert_load_to_dict(self):
        # list of dicts
        self.load = []
        for i in range(len(self.load_type_btn)):
            if self.load_by_coord[i]:
                aux_key = ["x_coord", "y_coord", "x_direc", "y_direc", "force"]
                aux_val = [self.load_coord[i][0], self.load_coord[i][1], self.load_x_dir[i], self.load_y_dir[i], self.load_value[i]]
            else:
                aux_key = ["coord", "axis", "eps", "x_direc", "y_direc", "force"]
                aux_val = [self.load_coord[i], self.x_load_col[i], self.load_error[i], self.load_x_dir[i], self.load_y_dir[i], self.load_value[i]]
            dicti = dict(zip(aux_key, aux_val))
            self.load.append(dicti)

# Nodes displacements constrained
    def create_node_constrain_btn(self):
        by_coord_node_constrain_btn = QtWidgets.QRadioButton("Constrain node displacement by coordinate")
        by_column_node_constrain_btn = QtWidgets.QRadioButton("Constrain node displacement by column")
        self.node_constrain_type_btn.append([by_coord_node_constrain_btn, by_column_node_constrain_btn])

        coord_node_constrain_x_line  = QtWidgets.QLineEdit()
        coord_node_constrain_y_line  = QtWidgets.QLineEdit()
        coord_node_constrain_line = QtWidgets.QLineEdit()
        self.node_constrain_coord_btn.append([coord_node_constrain_x_line, coord_node_constrain_y_line, coord_node_constrain_line])
        
        x_by_col_node_constrain_btn = QtWidgets.QRadioButton("X")
        y_by_col_node_constrain_btn = QtWidgets.QRadioButton("Y")
        self.node_constrain_col_btn.append([x_by_col_node_constrain_btn, y_by_col_node_constrain_btn])

        erro_margin_node_constrain_line  = QtWidgets.QLineEdit()
        self.node_constrain_error_btn.append(erro_margin_node_constrain_line)

        node_constrain_x_dir_line_yes = QtWidgets.QRadioButton("Yes") 
        node_constrain_x_dir_line_no = QtWidgets.QRadioButton("No")
        self.node_constrain_x_dir_btn.append([node_constrain_x_dir_line_yes, node_constrain_x_dir_line_no])
   
        node_constrain_y_dir_line_yes  = QtWidgets.QRadioButton("Yes") 
        node_constrain_y_dir_line_no  = QtWidgets.QRadioButton("No")         
        self.node_constrain_y_dir_btn.append([node_constrain_y_dir_line_yes, node_constrain_y_dir_line_no])

    def add_node_constrain_btn(self, layout):
        node_constrain_type = QtWidgets.QButtonGroup(layout)
        node_constrain_label = QtWidgets.QLabel('---- Constrain node displacement ----')
        node_constrain_label.setAlignment(QtCore.Qt.AlignCenter)
        layout.addRow(node_constrain_label)
        
        node_constrain_type.addButton(self.node_constrain_type_btn[-1][0])
        layout.addRow(self.node_constrain_type_btn[-1][0])
       
        layout.addRow(QtWidgets.QLabel('X-coor'), self.node_constrain_coord_btn[-1][0])
        layout.addRow(QtWidgets.QLabel('Y-coord'), self.node_constrain_coord_btn[-1][1])
       
        node_constrain_type.addButton(self.node_constrain_type_btn[-1][1])
        layout.addRow(self.node_constrain_type_btn[-1][1])
       
        layout.addRow(QtWidgets.QLabel('Coord'), self.node_constrain_coord_btn[-1][2])
     
        layout.addRow(QtWidgets.QLabel('Column'))
        column_node_constrain = QtWidgets.QButtonGroup(layout)
        column_node_constrain.addButton(self.node_constrain_col_btn[-1][0])
        column_node_constrain.addButton(self.node_constrain_col_btn[-1][1])
        layout.addRow(self.node_constrain_col_btn[-1][0], self.node_constrain_col_btn[-1][1])
      
        layout.addRow(QtWidgets.QLabel('Error margin'), self.node_constrain_error_btn[-1])
       
        x_dir = QtWidgets.QButtonGroup(layout)
        x_dir.addButton(self.node_constrain_x_dir_btn[-1][0])
        x_dir.addButton(self.node_constrain_x_dir_btn[-1][1])
        layout.addRow(QtWidgets.QLabel('Contrain node displacement in X direction'))
        layout.addRow(self.node_constrain_x_dir_btn[-1][0], self.node_constrain_x_dir_btn[-1][1])

        y_dir = QtWidgets.QButtonGroup(layout)
        y_dir.addButton(self.node_constrain_y_dir_btn[-1][0])
        y_dir.addButton(self.node_constrain_y_dir_btn[-1][1])
        layout.addRow(QtWidgets.QLabel('Contrain node displacement in Y direction'))
        layout.addRow(self.node_constrain_y_dir_btn[-1][0], self.node_constrain_y_dir_btn[-1][1])

    def toggled_node_constrain(self):
        self.node_constrain_type_btn[-1][0].toggled.connect(self.node_constrain_coord_btn[-1][0].setEnabled)  
        self.node_constrain_type_btn[-1][0].toggled.connect(self.node_constrain_coord_btn[-1][1].setEnabled)
        self.node_constrain_type_btn[-1][0].toggled.connect(self.node_constrain_col_btn[-1][0].setDisabled)
        self.node_constrain_type_btn[-1][0].toggled.connect(self.node_constrain_col_btn[-1][1].setDisabled)
        self.node_constrain_type_btn[-1][0].toggled.connect(self.node_constrain_error_btn[-1].setDisabled)
        self.node_constrain_type_btn[-1][0].toggled.connect(self.node_constrain_coord_btn[-1][2].setDisabled)

        self.node_constrain_type_btn[-1][1].toggled.connect(self.node_constrain_coord_btn[-1][0].setDisabled)  
        self.node_constrain_type_btn[-1][1].toggled.connect(self.node_constrain_coord_btn[-1][1].setDisabled)
        self.node_constrain_type_btn[-1][1].toggled.connect(self.node_constrain_col_btn[-1][0].setEnabled)
        self.node_constrain_type_btn[-1][1].toggled.connect(self.node_constrain_col_btn[-1][1].setEnabled)
        self.node_constrain_type_btn[-1][1].toggled.connect(self.node_constrain_error_btn[-1].setEnabled)
        self.node_constrain_type_btn[-1][1].toggled.connect(self.node_constrain_coord_btn[-1][2].setEnabled)

    def set_default_node_constrain(self):
        self.node_constrain_type_btn[-1][0].setChecked(False)
        self.node_constrain_type_btn[-1][1].setChecked(True)

        self.node_constrain_coord_btn[-1][0].setDisabled(True)
        self.node_constrain_coord_btn[-1][1].setDisabled(True)
        
        self.node_constrain_coord_btn[-1][2].setEnabled(True)
        self.node_constrain_col_btn[-1][0].setEnabled(True)
        self.node_constrain_col_btn[-1][1].setEnabled(True)
        self.node_constrain_error_btn[-1].setEnabled(True)

        self.node_constrain_coord_btn[-1][2].setText('0')
        self.node_constrain_col_btn[-1][0].setChecked(True)
        self.node_constrain_error_btn[-1].setText('1e-6')

        self.node_constrain_x_dir_btn[-1][0].setChecked(True)
        self.node_constrain_y_dir_btn[-1][0].setChecked(True)

        self.toggled_node_constrain()

    def update_default_node_constrain(self, ind):
        node_constrain_by_coord = True if self.node_constrain_by_coord[ind] else False
        self.node_constrain_type_btn[-1][0].setChecked(node_constrain_by_coord)
        self.node_constrain_type_btn[-1][1].setChecked(not node_constrain_by_coord)

        if node_constrain_by_coord:
            self.node_constrain_coord_btn[-1][0].setText(str(self.node_constrain_coord[ind][0]))
            self.node_constrain_coord_btn[-1][1].setText(str(self.node_constrain_coord[ind][1]))

            self.node_constrain_coord_btn[-1][2].setDisabled(True)
            self.node_constrain_col_btn[-1][0].setDisabled(True)
            self.node_constrain_col_btn[-1][1].setDisabled(True)
            self.node_constrain_error_btn[-1].setDisabled(True)
        else:
            self.node_constrain_coord_btn[-1][0].setDisabled(True)
            self.node_constrain_coord_btn[-1][1].setDisabled(True)
        
            self.node_constrain_coord_btn[-1][2].setText(str(self.node_constrain_coord[ind]))
            col = True if self.x_node_constrain_col[ind] else False
            self.node_constrain_col_btn[-1][0].setChecked(col)
            self.node_constrain_col_btn[-1][1].setChecked(not col)
            self.node_constrain_error_btn[-1].setText(str(self.node_constrain_error[ind]))

        x_ind = 0 if self.node_constrain_x_dir[ind] else 1
        self.node_constrain_x_dir_btn[-1][x_ind].setChecked(True)

        y_ind = 0 if self.node_constrain_y_dir[ind] else 1  
        self.node_constrain_y_dir_btn[-1][y_ind].setChecked(True)

        self.toggled_node_constrain()

    def rewrite_node_constrain(self, layout):
        self.reset_node_constrain_list(vals=False)
        for ind in range(len(self.node_constrain_by_coord)):
            self.create_node_constrain_btn()
            self.update_default_node_constrain(ind)
            self.add_node_constrain_btn(layout)

    def reset_node_constrain_list(self, btn=True, vals=True):
        if btn:
            self.node_constrain_type_btn = []
            self.node_constrain_coord_btn = []
            self.node_constrain_col_btn = []
            self.node_constrain_error_btn = []
            self.node_constrain_x_dir_btn = []
            self.node_constrain_y_dir_btn = []
        if vals:
            self.node_constrain_by_coord = [] 
            self.node_constrain_coord = []
            self.x_node_constrain_col = []
            self.node_constrain_error = []
            self.node_constrain_x_dir = []
            self.node_constrain_y_dir = []

    def check_node_constrain_btn(self):
        self.warnings_node_constrain = []
        for i in range(len(self.node_constrain_type_btn)):
            if self.node_constrain_type_btn[i][0].isChecked():
                self.check_param(self.warnings_node_constrain, self.node_constrain_coord_btn[i][0].text(), [int, float], 'x coordinate must be an integer or float')
                self.check_param(self.warnings_node_constrain, self.node_constrain_coord_btn[i][1].text(), [int, float], 'y coordinate must be an integer or float')
            else:
                if not (self.node_constrain_col_btn[-1][0].isChecked()) and not (self.node_constrain_col_btn[-1][1].isChecked()):
                    warning = QtWidgets.QLabel("select column")
                    self.warnings_node_constrain.append(warning)
                self.check_param(self.warnings_node_constrain, self.node_constrain_coord_btn[i][2].text(), [int, float], 'coordinate must be an integer or float')
                self.check_param(self.warnings_node_constrain, self.node_constrain_error_btn[i].text(), [int, float], 'error margin must be an integer or float')

    def check_node_constrain_values(self):
        self.warnings_node_constrain = []
        for i in range(len(self.node_constrain_type_btn)):
            if self.node_constrain_type_btn[i][0].isChecked():
                if ast.literal_eval(self.node_constrain_coord_btn[i][0].text()) > self.lx:
                    warning = QtWidgets.QLabel("x coordinate exceeds mesh boundaries")
                    self.warnings_node_constrain.append(warning)
                if ast.literal_eval(self.node_constrain_coord_btn[i][1].text()) > self.ly:
                    warning = QtWidgets.QLabel("y coordinate exceeds mesh boundaries")
                    self.warnings_node_constrain.append(warning)
            else:
                if self.node_constrain_col_btn[-1][0].isChecked():
                    if ast.literal_eval(self.node_constrain_coord_btn[i][2].text()) > self.lx:
                        warning = QtWidgets.QLabel("coordinate exceeds mesh boundaries")
                        self.warnings_node_constrain.append(warning)
                else:
                    if ast.literal_eval(self.node_constrain_coord_btn[i][2].text()) > self.ly:
                        warning = QtWidgets.QLabel("coordinate exceeds mesh boundaries")
                        self.warnings_node_constrain.append(warning)

    def update_node_constrain(self):
        for i in range(len(self.node_constrain_type_btn)):
            node_constrain_by_coord = True if self.node_constrain_type_btn[i][0].isChecked() else False
            self.node_constrain_by_coord.append(node_constrain_by_coord)
            if node_constrain_by_coord:
                x = ast.literal_eval(self.node_constrain_coord_btn[i][0].text())
                y = ast.literal_eval(self.node_constrain_coord_btn[i][1].text())
                self.node_constrain_coord.append([x,y])
                self.x_node_constrain_col.append(None)
                self.node_constrain_error.append(None)
            else:
                coord = ast.literal_eval(self.node_constrain_coord_btn[i][2].text())
                self.node_constrain_coord.append(coord)
                col = 1 if self.node_constrain_col_btn[i][0].isChecked() else 2
                self.x_node_constrain_col.append(col)
                error = ast.literal_eval(self.node_constrain_error_btn[i].text())
                self.node_constrain_error.append(error)
            
            x = 1 if self.node_constrain_x_dir_btn[i][0].isChecked() else 0
            self.node_constrain_x_dir.append(x)

            y = 1 if self.node_constrain_y_dir_btn[i][0].isChecked() else 0
            self.node_constrain_y_dir.append(y)

    def convert_node_constrain_to_dict(self):
        # list of dicts
        self.node_constrain = []
        for i in range(len(self.node_constrain_type_btn)):
            if self.node_constrain_by_coord[i]:
                aux_key = ["x_coord", "y_coord", "constrain_disp_x", "constrain_disp_y"]
                aux_val = [self.node_constrain_coord[i][0], self.node_constrain_coord[i][1], self.node_constrain_x_dir[i], self.node_constrain_y_dir[i]]
            else:
                aux_key = ["coord", "axis", "eps", "constrain_disp_x", "constrain_disp_y"]
                aux_val = [self.node_constrain_coord[i], self.x_node_constrain_col[i], self.node_constrain_error[i], self.node_constrain_x_dir[i], self.node_constrain_y_dir[i]]
            dicti = dict(zip(aux_key, aux_val))
            self.node_constrain.append(dicti)

class ParametersFEM2D(Parameters):
    def __init__(self):
        self.create_btns()
        self.set_default()
        self.update_params()

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
        self.x_coord_plot_btn = QtWidgets.QLineEdit()
        self.y_coord_plot_btn = QtWidgets.QLineEdit()
        self.x_dir_plot_btn = QtWidgets.QRadioButton("X")
        self.y_dir_plot_btn = QtWidgets.QRadioButton("Y")       

    def add_btns(self, layout):
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

        self.save_check.setChecked(self.save)

        self.toggled_fem2d()

    def set_node_plot(self):
        if self.freqrsp:
            self.node_plot = [self.x_coord_plot, self.y_coord_plot, self.x_dir_plot, self.y_dir_plot]
        else:
            self.node_plot = None

    def export_param(self):

        self.set_node_plot()

        param = {"nelx":self.nelx, "nely":self.nely, "lx":self.lx, "ly":self.ly, "E":self.E, "v":self.v, "rho":self.rho,
                "alpha":self.alpha, "beta":self.beta, "eta":self.eta, "factor":self.factor, "freq":self.freq, 
                "freqrsp":self.freqrsp, "freq_range":self.freq_range, "load_matrix":self.load, 
                "constr_matrix":self.node_constrain, "node_plot":self.node_plot, "save":self.save, "mesh_file":None}
        
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
        if self.freqrsp_check.isChecked():
            self.check_param(self.warnings, self.x_coord_plot_btn.text(), [int, float], 'Node to plot: node coordinate must be an integer or float')
            self.check_param(self.warnings, self.y_coord_plot_btn.text(), [int, float], 'Node to plot: node coordinate must be an integer or float')
        if len(self.warnings) == 0:
            self.check_node_plot()

    def check_node_plot(self):
        if self.freqrsp_check.isChecked():
            if ast.literal_eval(self.x_coord_plot_btn.text()) > self.lx:
                warning = QtWidgets.QLabel("Node to plot: x coordinate exceeds mesh boundaries")
                self.warnings.append(warning)
            if ast.literal_eval(self.y_coord_plot_btn.text()) > self.ly:
                warning = QtWidgets.QLabel("Node to plot: y coordinate exceeds mesh boundaries")
                self.warnings.append(warning)

class ParametersOpt(Parameters):
    def __init__(self):
        super().__init__()

        self.create_btns()
        self.set_default()
        #self.update_params()

    def export_param(self):
        param = {"nelx":self.nelx, "nely":self.nely, "lx":self.lx, "ly":self.ly, "E":self.E, "v":self.v, "rho":self.rho,
                "alpha_par":self.alpha, "beta_par":self.beta, "eta_par":self.eta, "factor":self.factor, "freq":self.freq, 
                "freqrsp":self.freqrsp, "freq_range":self.freq_range, "load_matrix":self.load,
                "constr_matrix":self.node_constrain, "save":self.save, "mesh_file":None, "mma":self.mma,
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
    def create_btns(self):
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

        self.freqrsp_check = QtWidgets.QCheckBox("Frequency Response")  
        self.freq_range_spin = QtWidgets.QLineEdit()
        self.alpha_plot_spin = QtWidgets.QLineEdit()
        self.beta_plot_spin  = QtWidgets.QLineEdit()
        self.eta_plot_spin   = QtWidgets.QLineEdit()

        self.max_iter_spin = QtWidgets.QLineEdit() 
        self.save_check = QtWidgets.QCheckBox("Save Data")
        self.dens_filter_check = QtWidgets.QCheckBox("Density Filter")      
        self.mesh_deform_check = QtWidgets.QCheckBox("Plot Deformed Mesh")
        self.factor_spin = QtWidgets.QLineEdit()
    
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
        self.save = self.save_check.checkState()

        self.dens_filter  = self.dens_filter_check.checkState()
        self.mesh_deform  = self.mesh_deform_check.checkState()

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

        if self.passive_coord is not None:
            self.passive_coord_spin.setText(str(self.passive_coord))
        if self.modes is not None:
            self.modes_spin.setText(str(self.modes))
        self.const_func_spin.setText(str(self.const_func))
        self.n1_spin.setText(str(self.n1))
        self.freq_spin.setText(str(self.freq))       
        
        self.max_iter_spin.setText(str(self.max_iter))
        self.factor_spin.setText(str(self.factor))

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

        if self.save:
            self.save_check.setChecked(True)

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

class TextFem2d():
    def __init__(self):
        self.editor = QtWidgets.QLabel()
        self.editor.setTextInteractionFlags(QtCore.Qt.LinksAccessibleByMouse)
        self.editor.setTextFormat(QtCore.Qt.RichText)
        self.editor.setTextInteractionFlags(QtCore.Qt.TextBrowserInteraction)
        self.editor.setOpenExternalLinks(True)
        self.text_fem = """ <p> <b> <font size="+7"> Parameters: </font> </b>
                <hr>
                <p> <b> <font size="+1"> Nelx </font> </b> <font size="+1"> (int): Number of elements on the x-axis.  </font> </p>
                <p> <b> <font size="+1"> Nely </font> </b> <font size="+1"> (int): Number of elements on the y-axis.</font> </p>
                <p> <b> <font size="+1"> Lx </font> </b> <font size="+1"> (float): X-axis length.</font> </p>
                <p> <b> <font size="+1"> Ly </font> </b> <font size="+1"> (float): Y-axis length.</font> </p>           
                <p> <b> <font size="+1"> E </font> </b> <font size="+1"> (float): Elastic modulus. </font> </p>
                <p> <b> <font size="+1"> v </font> </b> <font size="+1"> (float): Poisson's ratio. </font> </p>
                <p> <b> <font size="+1"> rho </font> </b> <font size="+1"> (float): Density. </font> </p>
                <p> <b> <font size="+1"> Alpha </font> </b> <font size="+1"> (float): Damping coefficient proportional to mass. </font> </p>
                <p> <b> <font size="+1"> Beta </font> </b> <font size="+1"> (float): Damping coefficient proportional to stiffness. </font> </p> 
                <p> <b> <font size="+1"> Eta </font> </b> <font size="+1"> (float): Damping coefficient. </font> </p>
                <p> <b> <font size="+1"> Factor </font> </b> <font size="+1"> (float): Factor to deform the mesh. </font> </p>
                <p> <b> <font size="+1"> Frequency </font> </b> <font size="+1"> (int): Optimized frequency. </font> </p>
                <p> <b> <font size="+1"> Save Data </font> </b> <font size="+1"> (bool): if True save the optimization and frequency response graphs as PNG. </font> </p>
                <p> <b> <font size="+1"> Frequency Response </font> </b>  <font size="+1"> It's necessary to pass the frequency range and the node that will be plotted.  </font> </p>
                    <p style="margin-left:2em"> <font size="+1"> <b> Frequency range </b> (list): It's a list with three values. </font> </p>
                        <p style="margin-left:4em"> <font size="+1"> - First value is the minimum frequency of the graph. </font> </p>
                        <p style="margin-left:4em"> <font size="+1"> - Second value is the maximum frequency of the graph. </font> </p>
                        <p style="margin-left:4em"> <font size="+1"> - Third value is the step between each calculation of the objective function. </font> <</p>
                
                    <p style="margin-left:2em"> <font size="+1"> <b> Node to plot </b> : Node that will be calculated the frequency response. </font> </p>
                        <p style="margin-left:4em"> <font size="+1"> - X-coord (float): X-axis coordinate of the node. </font> </p>
                        <p style="margin-left:4em"> <font size="+1"> - Y-coord (float): Y-axis coordinate of the node. </font> </p>   
                        <p style="margin-left:4em"> <font size="+1"> - X or Y: The direction. </font> </p>   
                """
        self.set_text()

    def set_text(self):
        self.editor.setText(self.text_fem)

class TextOpt():
    def __init__(self):
        self.editor = QtWidgets.QLabel()
        self.editor.setTextInteractionFlags(QtCore.Qt.LinksAccessibleByMouse)
        self.editor.setTextFormat(QtCore.Qt.RichText)
        self.editor.setTextInteractionFlags(QtCore.Qt.TextBrowserInteraction)
        self.editor.setOpenExternalLinks(True)

        self.text_opt = """ <p> <b> <font size="+7"> Parameters: </font> </b>
            <hr>
            <p> <b> <font size="+1"> Nelx </b> <font size="+1"> (int): Number of elements on the X-axis.  </font> </p>
            <p> <b> <font size="+1"> Nely </b> <font size="+1"> (int): Number of elements on the Y-axis.  </font> </p>
            <p> <b> <font size="+1"> Lx </b> <font size="+1"> (float): X-axis length.  </font> </p>
            <p> <b> <font size="+1"> Ly </b> <font size="+1"> (float): Y-axis length.  </font> </p>
            <p> <b> <font size="+1"> E </b> <font size="+1"> (float): Elastic modulus. </font> </p>
            <p> <b> <font size="+1"> v </b> <font size="+1"> (float): Poisson's ratio. </font> </p>
            <p> <b> <font size="+1"> rho </b> <font size="+1"> (float): Density.
            <p> <b> <font size="+1"> Ratio Factor </b> <font size="+1"> (float): Factor applied in the radius to get elements in the vicinity of each element. </font> </p> 
            <p> <b> <font size="+1"> x_min_mass </b> <font size="+1"> (float): Minimum relative densities to mass.  </font> </p>
            <p> <b> <font size="+1"> x_min_stif </b> <font size="+1"> (float): Minimum relative densities to stiffness.  </font> </p>
            <font size="+1"> Penal. Stiffness </b> <font size="+1"> (int): Penalization power to stiffness.  </font> </p>
            <p> <b> <font size="+1"> Penal. Mass </b> <font size="+1"> (int): Penalization power to mass.  </font> </p>
            <p> <b> <font size="+1"> Alpha </b> <font size="+1"> (float): Damping coefficient proportional to mass.   </font> </p>
            <p> <b> <font size="+1"> Beta </b> <font size="+1"> (float): Damping coefficient proportional to stiffness.   </font> </p>
            <p> <b> <font size="+1"> Eta </b> <font size="+1"> (float): Damping coefficient.   </font> </p>
            <p> <b> <font size="+1"> Passive Coordinates </b> <font size="+1"> (tuple): Region that the shape will not be changed.   </font> </p>
                <p style="margin-left:2em"> <font size="+1"> Example: </font> </p>
                <p style="margin-left:4em"> <font size="+1"> ((0.5, 1), (0.3, 0.6)) = ((x_initial, x_final), (y_initial, y_final))  </font> </p>
            <p> <b> <font size="+1"> Modes </b> <font size="+1"> (int): If not None is used the Mode Superposition Method to calculate the displacement. </font> </p>
            <p> <b> <font size="+1"> Constant Function </b> <font size="+1"> (float):  </font> </p>
            <p> <b> <font size="+1"> Objective Function Weight </b> <font size="+1"> (float): Weight associated with the objective function.. </font> </p>
                    <p style="margin-left:2em"> <font size="+1"> - If n1 &#60; 0: Maximize objective function.  </font> </p>
                    <p style="margin-left:2em"> <font size="+1"> - If n1 > 0: Minimize objective function.  </font> </p>
            <p> <b> <font size="+1"> Frequency </b> <font size="+1"> (int): Optimized frequency. </font> </p>
            <p> <b> <font size="+1"> Objective Function </b> <font size="+1"> (str): Objective function used.  </font> </p>
                <p style="margin-left:2em"> <font size="+1"> If the multiobjective function is being calculated, weight n1 is assigned.  </font> </p>
            <p> <b> <font size="+1"> Multiobjective Function </b> <font size="+1"> (tuple): Second function calculated.  </font> </p>
                    <p style="margin-left:2em"> <font size="+1"> The assigned weight is (1 - n1).  </font> </p>            
            <p> <b> <font size="+1"> Multiobjective Frequency </b> <font size="+1"> (tuple): frequency that the multiobjective function is being optimized.  </font> </p>
            <p> <b> <font size="+1"> Max. Iterations </b> <font size="+1"> (int): Number of iterations. </font> </p>
            <p> <b> <font size="+1"> Save Data </b> <font size="+1"> (bool): if checked saves: </font> </p>
                <p style="margin-left:2em"> <font size="+1"> - Xval, objective function, constraint function and frequency response values.  </font> </p>
                <p style="margin-left:2em"> <font size="+1"> - Optimization part, deformed mesh and frequency response graphs as PNG.  </font> </p>
            <p> <b> <font size="+1"> Density Filter </b> <font size="+1"> (bool): If checked uses density filter. Otherwise uses sensitivity filter.  </font> </p>
            <p> <b> <font size="+1"> Plot Deformed Mesh </b> <font size="+1"> (bool): If checked plots the mesh deformation of the dynamic function.  </font> </p>
                <p style="margin-left:2em"> <font size="+1"> <b> Factor </b> (float): Factor to deform the mesh.  </font> </p>
                    
            <p> <b> <font size="+1"> Frequency Response </b> <font size="+1">: It's necessary to pass the frequency range and the coefficient to plot.
                <p style="margin-left:2em"> <font size="+1"> <b> Frequency Range </b> (list): It's a list with three values. </font> </p>
                    <p style="margin-left:4em"> <font size="+1"> - First value is the minimum frequency of the graph. </font> </p>
                    <p style="margin-left:4em"> <font size="+1"> - Second value is the maximum frequency of the graph. </font> </p>
                    <p style="margin-left:4em"> <font size="+1"> - Third value is the step between each calculation of the objective function. </font> <</p>
            <p style="margin-left:2em"> <font size="+1"> <b> Alpha Plot </b> <font size="+1"> (float): Damping coefficient proportional to mass.   </font> </p>
            <p style="margin-left:2em"> <font size="+1"> <b> Beta Plot </b> <font size="+1"> (float): Damping coefficient proportional to stiffness.   </font> </p>
            <p style="margin-left:2em"> <font size="+1"> <b> Eta Plot </b> <font size="+1"> (float): Damping coefficient.   </font> </p>
            """
        self.set_text()

    def set_text(self):
        self.editor.setText(self.text_opt)

class TextBc():
    def __init__(self):
        self.text_load = """<p> <b> <font size="+7"> Add loads </font> </b>
            <hr>
            <p style="margin-left:2em"> <font size="+3"> Add loads by coordinate: </font> </p>
                    
                <p style="margin-left:4em"> <font size="+1"> <b> X-coord </b> (float): coordinate in X-axis. </font> </p>
                <p style="margin-left:4em"> <font size="+1"> <b> Y-coord </b> (float): coordinate in Y-axis. </font> </p>
            
            <p style="margin-left:2em"> <font size="+3"> Add loads by column: </font> </p>
                    
                <p style="margin-left:4em"> <font size="+1"> <b> Coord </b> (float): Coordinate. </font> </p>
                <p style="margin-left:4em"> <font size="+1"> <b> Column </b>: Direction that the load is applied. It can be: X-Axis or Y-Axis. </font> </p>
                <p style="margin-left:4em"> <font size="+1"> <b> Error margin </b> (float): Margin of error. The more discretized the mesh, the smaller the error margin </font> </p>

            <p style="margin-left:2em"> <font size="+1"> <b> Load direction </b>: </font> </p>

                <p style="margin-left:4em"> <font size="+1"> Positive: Indicates that the load is applied in the positive direction. </font> </p>
                <p style="margin-left:4em"> <font size="+1"> Negative: Indicates that the load is applied in the negative direction. </font> </p>
                <p style="margin-left:4em"> <font size="+1"> None: Indicates that no load is applied. </font> </p>
            
            <p style="margin-left:2em"> <font size="+1"> <b> Load value </b>: The module of the load in Newton. </font> </p>            
            """
        self.editor_load = QtWidgets.QLabel()
        self.editor_load.setTextInteractionFlags(QtCore.Qt.LinksAccessibleByMouse)
        self.editor_load.setTextFormat(QtCore.Qt.RichText)
        self.editor_load.setTextInteractionFlags(QtCore.Qt.TextBrowserInteraction)
        self.editor_load.setOpenExternalLinks(True)
        self.set_load_text()

        self.text_node_constrain = """<p> <b> <font size="+7"> Constrain nodes displacements </font> </b>
            <hr>
            <p style="margin-left:2em"> <font size="+3"> Constrain nodes displacements by coordinate: </font> </p>
                    
                <p style="margin-left:4em"> <font size="+1"> <b> X-coord </b> (float): coordinate in X-axis. </font> </p>
                <p style="margin-left:4em"> <font size="+1"> <b> Y-coord </b> (float): coordinate in Y-axis. </font> </p>
            
            <p style="margin-left:2em"> <font size="+3"> Constrain nodes displacements by column: </font> </p>
                    
                <p style="margin-left:4em"> <font size="+1"> <b> Coord </b> (float): Coordinate. </font> </p>
                <p style="margin-left:4em"> <font size="+1"> <b> Column </b>: Direction that the load is applied. It can be: X-Axis or Y-Axis. </font> </p>
                <p style="margin-left:4em"> <font size="+1"> <b> Error margin </b> (float): Margin of error. The more discretized the mesh, the smaller the error margin. </font> </p>

            <p style="margin-left:2em"> <font size="+1"> <b> Direction: </b> </font> </p>

                <p style="margin-left:4em"> <font size="+1"> Yes: Indicates that the node displacement is constrained. </font> </p>
                <p style="margin-left:4em"> <font size="+1"> No: Indicates that no node displacement is constrained. </font> </p>
            """
        self.editor_node_constrain = QtWidgets.QLabel()
        self.editor_node_constrain.setTextInteractionFlags(QtCore.Qt.LinksAccessibleByMouse)
        self.editor_node_constrain.setTextFormat(QtCore.Qt.RichText)
        self.editor_node_constrain.setTextInteractionFlags(QtCore.Qt.TextBrowserInteraction)
        self.editor_node_constrain.setOpenExternalLinks(True)
        self.set_node_constrain_text()

    def set_load_text(self):
        self.editor_load.setText(self.text_load)

    def set_node_constrain_text(self):
        self.editor_node_constrain.setText(self.text_node_constrain)

class TextConstraint():
    def __init__(self):
        self.editor = QtWidgets.QLabel()
        self.editor.setTextInteractionFlags(QtCore.Qt.LinksAccessibleByMouse)
        self.editor.setTextFormat(QtCore.Qt.RichText)
        self.editor.setTextInteractionFlags(QtCore.Qt.TextBrowserInteraction)
        self.editor.setOpenExternalLinks(True)

        self.text_constraint = """  """
        self.set_text()
        
    def set_text(self):
        self.editor.setText(self.text_constraint)