from PyQt5 import QtCore, QtWidgets
import numpy as np
import ast

class Verification():

    def check_param(self, warnings, input_val, types, war):
        try:
            if input_val:
                isinstance(ast.literal_eval(input_val), types)
            else:
                warning = QtWidgets.QLabel(war)
                warnings.append(warning)
        except:
            warning = QtWidgets.QLabel(war)
            warnings.append(warning)

class ParamLoad(Verification):
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

        self.warnings = []

    # Load
    def create_btn(self):
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

    def add_btn(self, layout):
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

    def toggled(self):
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

    def set_size(self, lx, ly):
        self.lx = lx
        self.ly = ly

    def set_default(self):
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

        self.toggled()

    def update_default(self, ind):
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

        self.toggled()

    def rewrite(self, layout):
        self.reset_list(vals=False)
        for ind in range(len(self.load_by_coord)):
            self.create_btn()
            self.update_default(ind)
            self.add_btn(layout)
        
    def reset_list(self, btn=True, vals=True):
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

    def check_btn(self):
        self.warnings = []
        for i in range(len(self.load_type_btn)):
            if self.load_type_btn[i][0].isChecked():
                self.check_param(self.warnings, self.load_coord_btn[i][0].text(), (int, float), 'x coordinate must be an integer or float')
                self.check_param(self.warnings, self.load_coord_btn[i][1].text(), (int, float), 'y coordinate must be an integer or float')
            else:
                if not (self.load_col_btn[-1][0].isChecked()) and not (self.load_col_btn[-1][1].isChecked()):
                    warning = QtWidgets.QLabel("select column")
                    self.warnings.append(warning)
                self.check_param(self.warnings, self.load_coord_btn[i][2].text(), (int, float), 'coordinate must be an integer or float')
                self.check_param(self.warnings, self.load_error_btn[i].text(), (int, float), 'error margin must be an integer or float')
            
            self.check_param(self.warnings, self.load_value_btn[i].text(), (int, float), 'load value must be an integer or float')

    def check_values(self):
        self.warnings = []
        for i in range(len(self.load_type_btn)):
            if self.load_type_btn[i][0].isChecked():
                if ast.literal_eval(self.load_coord_btn[i][0].text()) > self.lx:
                    warning = QtWidgets.QLabel("x coordinate exceeds mesh boundaries")
                    self.warnings.append(warning)
                if ast.literal_eval(self.load_coord_btn[i][1].text()) > self.ly:
                    warning = QtWidgets.QLabel("y coordinate exceeds mesh boundaries")
                    self.warnings.append(warning)
            else:
                if self.load_col_btn[-1][0].isChecked():
                    if ast.literal_eval(self.load_coord_btn[i][2].text()) > self.lx:
                        warning = QtWidgets.QLabel("coordinate exceeds mesh boundaries")
                        self.warnings.append(warning)
                else:
                    if ast.literal_eval(self.load_coord_btn[i][2].text()) > self.ly:
                        warning = QtWidgets.QLabel("coordinate exceeds mesh boundaries")
                        self.warnings.append(warning)

    def update(self):
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

    def input_to_np(self):
        length = len(self.load_type_btn)
        self.np_array = np.empty((length, 6))

        for i in range(length):
            if self.load_by_coord[i]:
                aux = [self.load_coord[i][0], self.load_coord[i][1], self.load_x_dir[i], self.load_y_dir[i], self.load_value[i], np.nan]
            else:
                aux = [self.load_coord[i], self.x_load_col[i], self.load_x_dir[i], self.load_y_dir[i], self.load_value[i], self.load_error[i]]
            self.np_array[i, :] = aux

class ParamNodeConstrain(Verification):
    def __init__(self) -> None:
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

    # Nodes displacements constrained
    def create_btn(self):
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

    def add_btn(self, layout):
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

    def toggled(self):
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

    def set_default(self):
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

        self.toggled()

    def update_default(self, ind):
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

        self.toggled()

    def rewrite(self, layout):
        self.reset_list(vals=False)
        for ind in range(len(self.node_constrain_by_coord)):
            self.create_btn()
            self.update_default(ind)
            self.add_btn(layout)

    def reset_list(self, btn=True, vals=True):
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

    def check_btn(self):
        self.warnings = []
        for i in range(len(self.node_constrain_type_btn)):
            if self.node_constrain_type_btn[i][0].isChecked():
                self.check_param(self.warnings, self.node_constrain_coord_btn[i][0].text(), (int, float), 'x coordinate must be an integer or float')
                self.check_param(self.warnings, self.node_constrain_coord_btn[i][1].text(), (int, float), 'y coordinate must be an integer or float')
            else:
                if not (self.node_constrain_col_btn[-1][0].isChecked()) and not (self.node_constrain_col_btn[-1][1].isChecked()):
                    warning = QtWidgets.QLabel("select column")
                    self.warnings.append(warning)
                self.check_param(self.warnings, self.node_constrain_coord_btn[i][2].text(), (int, float), 'coordinate must be an integer or float')
                self.check_param(self.warnings, self.node_constrain_error_btn[i].text(), (int, float), 'error margin must be an integer or float')

    def set_size(self, lx, ly):
        self.lx = lx
        self.ly = ly

    def check_values(self):
        self.warnings = []
        for i in range(len(self.node_constrain_type_btn)):
            if self.node_constrain_type_btn[i][0].isChecked():
                if ast.literal_eval(self.node_constrain_coord_btn[i][0].text()) > self.lx:
                    warning = QtWidgets.QLabel("x coordinate exceeds mesh boundaries")
                    self.warnings.append(warning)
                if ast.literal_eval(self.node_constrain_coord_btn[i][1].text()) > self.ly:
                    warning = QtWidgets.QLabel("y coordinate exceeds mesh boundaries")
                    self.warnings.append(warning)
            else:
                if self.node_constrain_col_btn[-1][0].isChecked():
                    if ast.literal_eval(self.node_constrain_coord_btn[i][2].text()) > self.lx:
                        warning = QtWidgets.QLabel("coordinate exceeds mesh boundaries")
                        self.warnings.append(warning)
                else:
                    if ast.literal_eval(self.node_constrain_coord_btn[i][2].text()) > self.ly:
                        warning = QtWidgets.QLabel("coordinate exceeds mesh boundaries")
                        self.warnings.append(warning)

    def update(self):
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

    def input_to_np(self):
        length = len(self.node_constrain_type_btn)
        self.np_array = np.empty((length, 5))

        for i in range(length):
            if self.node_constrain_by_coord[i]:
                aux = [self.node_constrain_coord[i][0], self.node_constrain_coord[i][1], self.node_constrain_x_dir[i], self.node_constrain_y_dir[i], np.nan]
            else:
                aux = [self.node_constrain_coord[i], self.x_node_constrain_col[i], self.node_constrain_x_dir[i], self.node_constrain_y_dir[i], self.node_constrain_error[i]]
            self.np_array[i, :] = aux