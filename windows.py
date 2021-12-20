from PyQt5 import QtWidgets, QtCore, QtGui
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

import os
import sys
import re
import ast
import shutil

import matplotlib.pyplot as plt
plt.switch_backend('Agg')
import numpy as np
import pyqtgraph as pg

import plot_grid as mf
from plots_2d import PlotsFem2d
from plots_opt import PlotOpt
from parameters import ParametersOpt, ParametersText

class PopUp(QtWidgets.QDialog):
    def __init__(self, warnings):
        super().__init__()

        self.warnings = warnings

        self.setWindowTitle("Warning")
        self.layout = QtWidgets.QVBoxLayout()

        for warning in self.warnings:
            self.layout.addWidget(warning)
        self.setLayout(self.layout)

class MainWindow(QtWidgets.QDialog):
    def __init__(self):
        super(MainWindow, self).__init__()
        layout = QtWidgets.QHBoxLayout()

        button_fem2d = QtWidgets.QPushButton('FEM 2D', self)
        button_fem2d.clicked.connect(self.go_to_fem2d)
        layout.addWidget(button_fem2d)

        button_opt = QtWidgets.QPushButton('Optimization', self)
        button_opt.clicked.connect(self.go_to_opt)
        layout.addWidget(button_opt)

        #button_fem3d = QtWidgets.QPushButton('FEM 3D', self)
        #button_fem3d.clicked.connect(self.go_to_fem2d)
        #layout.addWidget(button_fem3d)
        self.setLayout(layout)
       
    def go_to_fem2d(self):
        param_2d = ParametersFEM2D()
        fem2d = WindowsFem2d(param_2d)
        widget.addWidget(fem2d)
        widget.setCurrentIndex(widget.currentIndex()+1)

    def go_to_opt(self):
        param_opt = ParametersOpt()
        opt = WindowsOptimization(param_opt)
        widget.addWidget(opt)
        widget.setCurrentIndex(widget.currentIndex()+1)

class SecondWindow(QtWidgets.QDialog):
    def __init__(self, param):
        super(SecondWindow, self).__init__()

        self.pbar = QtWidgets.QProgressBar()
        
        # lateral menu
        left_layout = QtWidgets.QVBoxLayout()
        self.btn_back_menu = QtWidgets.QPushButton('Back to Menu')
        self.btn_back_menu.clicked.connect(self.go_to_screen1)
        left_layout.addWidget(self.btn_back_menu)

        self.param = param
        #self.param.create_param_btns()
        #self.param.set_default()
        self.param.add_btns(left_layout)

        self.btn_to_load = QtWidgets.QPushButton('Next')
        self.btn_to_load.clicked.connect(self.go_to_load)
        left_layout.addWidget(self.btn_to_load)

        self._counter_button_run = 0
    
        left_layout.addStretch(5)
        left_layout.setSpacing(20)
        self.left_widget = QtWidgets.QWidget()
        self.left_widget.setLayout(left_layout)
        scroll = QtWidgets.QScrollArea()
        scroll.setWidget(self.left_widget)
        scroll.setWidgetResizable(True)

        # window 
        self.right_widget = QtWidgets.QTabWidget()
        
        self.param_text = ParametersText()
        self.param_text.set_text(self.param.text)
        
        self.tab1 = self.ui1()
        self.right_widget.addTab(self.tab1, '')
        
        self.tab2 = self.ui2()
        self.right_widget.addTab(self.tab2, '')
        
        self.right_widget.setCurrentIndex(0)
        self.right_widget.setStyleSheet('''QTabBar::tab{width: 0; \
            height: 0; margin: 0; padding: 0; border: none;}''')

        scroll_right = QtWidgets.QScrollArea()
        scroll_right.setWidget(self.right_widget)
        scroll_right.setWidgetResizable(True)

        # main window 
        main_layout = QtWidgets.QHBoxLayout()
        main_layout.addWidget(scroll)
        main_layout.addWidget(scroll_right)
        main_layout.setStretch(0, 56) #aqui aumenta a barra lateral
        main_layout.setStretch(1, 300)
        self.setLayout(main_layout) 

    # ------- pages -------
    def ui1(self):
        layout_ui1 = QtWidgets.QVBoxLayout()
        layout_ui1.addWidget(self.param_text.editor)
        layout_ui1.addStretch(20)
        
        main = QtWidgets.QWidget()
        main.setLayout(layout_ui1)
        return main
    
    def ui2(self):
        self.init_ui2_layout = QtWidgets.QVBoxLayout()
        main = QtWidgets.QWidget()
        main.setLayout(self.init_ui2_layout)
        return main

    def upd_ui2_layout(self, new_layout): #TODO:LEMBRAR DE MUDAR PARA O FEM2D
        if self.tab2.layout() is not None:
            QtWidgets.QWidget().setLayout(self.tab2.layout())
        self.tab2.setLayout(new_layout)

    def upd_left_layout(self, new_layout):
        if self.left_widget.layout() is not None:
            QtWidgets.QWidget().setLayout(self.left_widget.layout())
        self.left_widget.setLayout(new_layout)

    def ui_add_load(self, update=False):
        ui_layout = QtWidgets.QFormLayout()
        self.btn_back_param = QtWidgets.QPushButton('Back')
        self.btn_back_param.clicked.connect(self.back_to_param)
        ui_layout.addRow(self.btn_back_param)

        self.btn_to_node_constrain = QtWidgets.QPushButton('Next')
        self.btn_to_node_constrain.clicked.connect(self.go_to_node_constrain)
        ui_layout.addRow(self.btn_to_node_constrain)

        self.btn_add_new_load = QtWidgets.QPushButton('Add new load')
        self.btn_add_new_load.clicked.connect(self.add_new_load)
        ui_layout.addRow(self.btn_add_new_load)

        self.btn_reset_load = QtWidgets.QPushButton('Reset Load')
        self.btn_reset_load.clicked.connect(self.reset_load)
        ui_layout.addRow(self.btn_reset_load)
        if update:
            self.param.rewrite_load(ui_layout)
        else:
            self.param.reset_load_list()
            self.param.create_load_btn()
            self.param.set_default_load()
            self.param.add_load_btn(ui_layout)
        self.upd_left_layout(ui_layout)

    def ui_add_node_constrain(self):
        ui_layout = QtWidgets.QFormLayout()
        self.btn_back_load = QtWidgets.QPushButton('Back')
        self.btn_back_load.clicked.connect(self.back_to_load)
        ui_layout.addRow(self.btn_back_load)

        self.btn_to_constraint = QtWidgets.QPushButton('Next')
        self.btn_to_constraint.clicked.connect(self.go_to_constraint)
        ui_layout.addRow(self.btn_to_constraint)

        self.btn_add_new_node_constrain = QtWidgets.QPushButton('Constrain new node displacement')
        self.btn_add_new_node_constrain.clicked.connect(self.add_new_node_constrain)
        ui_layout.addRow(self.btn_add_new_node_constrain)

        self.btn_reset_node_constrain = QtWidgets.QPushButton('Reset constraint node displacement')
        self.btn_reset_node_constrain.clicked.connect(self.reset_node_constrain)
        ui_layout.addRow(self.btn_reset_node_constrain)

        self.param.reset_node_constrain_list()
        self.param.create_node_constrain_btn()
        self.param.set_default_node_constrain()
        self.param.add_node_constrain_btn(ui_layout)
        self.upd_left_layout(ui_layout)

    def ui_add_constraint(self):
        ui_layout = QtWidgets.QFormLayout()
        self.btn_back_constrain = QtWidgets.QPushButton('Back')
        self.btn_back_constrain.clicked.connect(self.back_to_node_constrain)
        ui_layout.addRow(self.btn_back_constrain)

        self.param.create_constraint()
        self.param.add_constraint_param(ui_layout)
        self.param.set_default_constraint()

        self.button_run = QtWidgets.QPushButton('Run')
        self.button_run.clicked.connect(lambda: self.run())
        ui_layout.addRow(self.button_run)

        self.button_stop = QtWidgets.QPushButton('Stop')
        self.button_stop.setEnabled(False)
        ui_layout.addRow(self.button_stop)
        self.upd_left_layout(ui_layout)

    # ------- buttons ------- 
    def go_to_screen1(self):
        mainwindow = MainWindow()
        widget.addWidget(mainwindow)
        widget.setCurrentIndex(widget.currentIndex()+1)

    def go_to_load(self):
        self.param.check_params()
        if len(self.param.warnings) != 0:
            dlg = PopUp(self.param.warnings)      
            dlg.exec_()
        else:
            self.param.update_params()
            # TODO:save parameters
            self.ui_add_load()

    def go_to_node_constrain(self):
        self.param.check_load_btn()
        if len(self.param.warnings_load) != 0:
            dlg = PopUp(self.param.warnings_load)      
            dlg.exec_()
        else:
            self.param.check_load_values()
            if len(self.param.warnings_load) != 0:
                dlg = PopUp(self.param.warnings_load)      
                dlg.exec_()
            else:
                self.param.update_load()
                self.ui_add_node_constrain()

    def go_to_constraint(self):
        self.param.check_node_constrain_btn()
        if len(self.param.warnings_node_constrain) != 0:
            dlg = PopUp(self.param.warnings_node_constrain)      
            dlg.exec_()
        else:
            self.param.check_node_constrain_values()
            if len(self.param.warnings_node_constrain) != 0:
                dlg = PopUp(self.param.warnings_node_constrain)      
                dlg.exec_()
            else:
                self.param.update_node_constrain()
                self.ui_add_constraint()

    def back_to_param(self):
        left_layout = QtWidgets.QVBoxLayout()
        self.btn_back_menu = QtWidgets.QPushButton('Back to Menu')
        self.btn_back_menu.clicked.connect(self.go_to_screen1)
        left_layout.addWidget(self.btn_back_menu)
        
        self.param.create_param_btns()
        self.param.update_default()
        self.param.add_btns(left_layout)

        self.btn_to_load = QtWidgets.QPushButton('Next')
        self.btn_to_load.clicked.connect(self.go_to_load)
        left_layout.addWidget(self.btn_to_load)
        self.upd_left_layout(left_layout)

    def back_to_load(self):
        self.ui_add_load(update=True) # TODO: Precisa mostrar os valores salvos
    
    def back_to_node_constrain(self):
        self.ui_add_node_constrain()

    def add_new_load(self):
        self.param.create_load_btn()
        self.param.set_default_load()
        self.param.add_load_btn(self.left_widget.layout()) #TODO: Adiciona um novo esquema

    def add_new_node_constrain(self):
        self.param.create_node_constrain_btn()
        self.param.set_default_node_constrain()
        self.param.add_node_constrain_btn(self.left_widget.layout()) #TODO: Adiciona um novo esquema

    def reset_load(self):
        self.ui_add_load()

    def reset_node_constrain(self):
        self.ui_add_node_constrain()

    def run(self):
        self.param.check_constraint()
        if len(self.param.warnings_constraint) != 0:
            dlg = PopUp(self.param.warnings_constraint)      
            dlg.exec_()
        else:
            self.button_stop.setEnabled(True)
            self.button_run.setEnabled(False)
    
            self.param.update_constraint()
            self.param.constraint_to_list()
            self.param.convert_load_to_dict()
            self.param.convert_node_constrain_to_dict()
            
            self._counter_button_run += 1
            if self._counter_button_run == 1:
                self.right_widget.setCurrentIndex(1)
                print("entrou")

    # ------- functions -------     
    def delete_temp(self):
        folder_name = 'temp'
        directory = os.path.join(os.path.dirname(__file__), folder_name)
        if os.path.exists(directory):
            shutil.rmtree(directory)

# Regular expressions to extract values from QProcess
progress_re = re.compile("Total complete: (\d+)%")

update_data_re = re.compile("Update data")

plot_freqresp = re.compile("Plot frequency response")

create_plot_re = re.compile("Create plot")

mesh_deform = re.compile("Plot deformed mesh")

save_figs = re.compile("Save figs")

def simple_percent_parser(output):
    """
    Matches lines using the update_data_re regex,
    returning a single integer for the % progress.
    """
    m = progress_re.search(output)
    if m:
        pc_complete = m.group(1)
        return int(pc_complete)

def simple_update_parser(output):
    """
    Matches lines using the progress_re regex,
    returning a bool.
    """
    m = update_data_re.search(output)
    if m:
        return True

def simple_freq_parser(output):
    """
    Matches lines using the progress_re regex,
    returning a bool.
    """
    m = plot_freqresp.search(output)
    if m:
        return True

def simple_save_parser(output):
    """
    Matches lines using the progress_re regex,
    returning a bool.
    """
    m = save_figs.search(output)
    if m:
        return True

def simple_create_parser(output):
    """
    Matches lines using the create_plot_re regex,
    returning a bool.
    """
    m = create_plot_re.search(output)
    if m:
        return True

def simple_deform_parser(output):
    """
    Matches lines using the create_plot_re regex,
    returning a bool.
    """
    m = mesh_deform.search(output)
    if m:
        return True

class WindowsOptimization(SecondWindow):
    def __init__(self, param):
        super().__init__(param)

        self.stop_thread = False
        
        #self.button_run.clicked.connect(lambda: self.run())
        self.param.mesh_deform_check.toggled.connect(self.param.factor_spin.setEnabled)   
        self.param.func_name2_box.activated.connect(self.param.freq2_spin.setEnabled) 
        self.param.freqrsp_check.toggled.connect(self.param.freq_range_spin.setEnabled)
        self.param.freqrsp_check.toggled.connect(self.param.alpha_plot_spin.setEnabled)
        self.param.freqrsp_check.toggled.connect(self.param.beta_plot_spin.setEnabled)
        self.param.freqrsp_check.toggled.connect(self.param.eta_plot_spin.setEnabled)

        # Boundary Conditions
        #TODO: Aqui preciso pegar de cada uma da lista
        
        self.directory = os.path.join(os.path.dirname(__file__))
        self.dir_temp = os.path.join(self.directory, 'temp')
    
    # ------- pages -------
    def add_canvas_ui2(self): 
        ui2_layout = QtWidgets.QVBoxLayout()
        ui2_layout.addWidget(self.pbar)

        self.graph_grid = pg.PlotWidget()
        self.grid = mf.PColorMeshItem(cmap='grey')
        self.graph_grid.setAspectLocked(True)
        self.graph_grid.hideAxis('bottom')
        self.graph_grid.hideAxis('left')
        self.graph_grid.addItem(self.grid)
        ui2_layout.addWidget(self.graph_grid)

        self.graph_conv = pg.PlotWidget()
        self.graph_conv.addLegend(labelTextColor=(0,0,0), offset=(800,10))
        self.graph_conv.setLabel('left', self.param.func_name.lower())
        self.graph_conv.setLabel('bottom', "iteration")
        #self.graph_conv.setFixedWidth(self.graph_grid.width())
        #self.graph_conv.setFixedHeight(self.graph_grid.height())
        ui2_layout.addWidget(self.graph_conv)

        self.text = QtWidgets.QPlainTextEdit()
        self.text.setReadOnly(True)

        ui2_layout.addWidget(self.text)

        self.upd_ui2_layout(ui2_layout)

    def add_freq_to_canvas_ui2(self): 
        ui2_layout = QtWidgets.QVBoxLayout()
        ui2_layout.addWidget(self.pbar)
        ui2_layout.addWidget(self.graph_grid)
        ui2_layout.addWidget(self.graph_conv)
        ui2_layout.addWidget(self.graph_freq)     
        ui2_layout.addWidget(self.text)

        self.upd_ui2_layout(ui2_layout)

    def add_deformed_mesh_to_canvas_ui2(self): 
        ui2_layout = QtWidgets.QVBoxLayout()
        ui2_layout.addWidget(self.pbar)
        ui2_layout.addWidget(self.graph_grid)
        ui2_layout.addWidget(self.graph_conv)
        if self.param.freqrsp:
            ui2_layout.addWidget(self.graph_freq)

        self.canvas_mesh.setFixedWidth(self.graph_grid.width())
        self.canvas_mesh.setFixedHeight(self.graph_grid.height())
        ui2_layout.addWidget(self.canvas_mesh)   
        ui2_layout.addWidget(self.text)

        self.upd_ui2_layout(ui2_layout)

    # ------- buttons -------
    def run(self):
        super().run()

        if len(self.param.warnings) == 0:
            pg.setConfigOption('background', 'w')
            pg.setConfigOption('foreground', 'k')

            self.add_canvas_ui2()
            self.param.export_param()

            self.process_opt = QtCore.QProcess()
            self.process_opt.readyReadStandardOutput.connect(self.handle_stdout)
            self.process_opt.readyReadStandardError.connect(self.handle_stderr)
            self.process_opt.stateChanged.connect(self.handle_state)
            self.process_opt.finished.connect(lambda: self.finish_process())
            if self.param.mma:
                self.process_opt.start("python", ['process_mma.py'])
            else:
                self.process_opt.start("python", ['process_gcmma.py'])
            
            self.button_stop.clicked.connect(lambda: self.stop_execution())

    def handle_stderr(self):
        data = self.process_opt.readAllStandardError()
        stderr = bytes(data).decode("utf8")
        # Extract progress if it is in the data.
        progress = simple_percent_parser(stderr)
        if progress:
            self.evt_update_progress(progress)
        
        updata_data = simple_update_parser(stderr)
        if updata_data:
            # Load values
            file = open(os.path.join(self.dir_temp, 'param_plot.txt'), "r")
            contents = file.read()
            param = ast.literal_eval(contents)
            file.close()  

            outeriter = param["outeriter"]
            outit = param["outit"]
            f0val = param["f0val"]
            fval = np.loadtxt(os.path.join(self.dir_temp, 'fval.txt'))
            xval = np.loadtxt(os.path.join(self.dir_temp, 'xval_log.txt'))
            xnew = np.loadtxt(os.path.join(self.dir_temp, 'xnew_plot.txt'))
            if os.path.exists(os.path.join(self.dir_temp, 'natural_freqs.txt')):
                natural_freqs = np.loadtxt(os.path.join(self.dir_temp, 'natural_freqs.txt'))
            else:
                natural_freqs = None

            # Update lists
            self.plot_opt.update_lists(outit, fval, f0val)

            # Plot values
            self.grid.setData(self.plot_opt.x_plot, self.plot_opt.y_plot, xnew.reshape(self.plot_opt.nelx, self.plot_opt.nely, order='F'))

            self.curves_conv[0].setData(self.plot_opt.list_iter[:outit + 1], self.plot_opt.list_f0val[:outit + 1])
            for ind in range(len(self.plot_opt.constr_func)):
                self.curves_conv[ind + 1].setData(self.plot_opt.list_iter[:outit + 1], self.plot_opt.list_fvals[:outit + 1, ind])

            out1, out2, out3 = self.evt_set_logger(outeriter, f0val, fval, xval, natural_freqs)
            self.set_message("Functions value: " + out1 + "\n")
            self.set_message("Xval vector: " + out2 + "\n")
            if out3 is not None:
                self.set_message("Natural frequencies: " + out3)

        plot_freq = simple_freq_parser(stderr)
        if plot_freq:
            orig = np.loadtxt(os.path.join(self.dir_temp, 'f_original.txt'))
            opt = np.loadtxt(os.path.join(self.dir_temp, 'f_optimized.txt'))
            self.graph_freq = pg.PlotWidget()
            self.graph_freq.setLabel('left', self.param.func_name.lower())
            self.graph_freq.setLabel('bottom', "frequency [Hz]")
            self.graph_freq.addLegend()
            self.graph_freq.plotItem.setLogMode(False, True)
            interval = np.arange(self.param.freq_range[0], self.param.freq_range[1] + 1, self.param.freq_range[2])
            self.graph_freq.plot(interval, orig, name='original', pen=pg.mkPen(color=(43,174,179), width=3, style=QtCore.Qt.DashLine))
            self.graph_freq.plot(interval, opt, name='optimized', pen=pg.mkPen(color=(64,66,114), width=2))
            self.add_freq_to_canvas_ui2()
        
        mesh_deform = simple_deform_parser(stderr)
        if mesh_deform:
            load_matrix = np.loadtxt(os.path.join(self.dir_temp, 'load_matrix.txt')).reshape(-1, 4)
            constr_matrix = np.loadtxt(os.path.join(self.dir_temp, 'constr_matrix.txt'), dtype=int)
            coord = np.loadtxt(os.path.join(self.dir_temp, 'coord.txt'))
            connect = np.loadtxt(os.path.join(self.dir_temp, 'connect.txt'), dtype=int)
            disp_vector = np.loadtxt(os.path.join(self.dir_temp, 'disp_vector.txt'))

            plot_mesh = PlotsFem2d(coord, connect)
            disp_vector = plot_mesh.change_disp_shape(disp_vector)
            coord_U = plot_mesh.apply_disp(disp_vector, self.param.factor)
            collection = plot_mesh.build_collection(coord_U)
        
            self.canvas_mesh = FigureCanvas(plt.Figure())
            self.canvas_mesh.figure.tight_layout()
            ax_mesh = self.canvas_mesh.figure.subplots()
            plot_mesh.plot_collection(ax_mesh, self.plot_opt.lx, self.plot_opt.ly, coord_U, collection, \
                                        load_matrix, constr_matrix)

            self.add_deformed_mesh_to_canvas_ui2()

        save = simple_save_parser(stderr)
        if save: 
            if not self.param.freqrsp:
                self.graph_freq = None
            if not self.param.mesh_deform:
                self.canvas_mesh = None
            self.plot_opt.save_figs(self.grid, self.graph_conv, self.graph_freq, self.canvas_mesh)
            
        create_plot = simple_create_parser(stderr)
        if create_plot:
            constr_values = np.loadtxt(os.path.join(self.dir_temp, 'constr_values.txt'))
            file = open(os.path.join(self.dir_temp, 'param_plot.txt'), "r")
            contents = file.read()
            param_plot = ast.literal_eval(contents)
            file.close()

            self.plot_opt = PlotOpt(param_plot["lx"], param_plot["ly"], param_plot["nelx"], param_plot["nely"], \
                                    param_plot["constr_func"], constr_values, self.param.max_iter)
            self.curves_conv = []
            self.curves_conv.append(self.graph_conv.plot(pen={'color': (0,0,0), 'width': 2}, name=self.param.func_name.lower()))
            for ind, f in enumerate(self.plot_opt.constr_func):
                pen_set = self.plot_opt.set_pen(f)
                self.curves_conv.append(self.graph_conv.plot(name=self.plot_opt.labels_constr[ind], pen=pen_set))
        
        self.set_message(stderr) # TODO: APAGAR DEPOIS

    def handle_stdout(self): #TODO: APAGAR DEPOIS
        data = self.process_opt.readAllStandardOutput()
        stdout = bytes(data).decode("utf8")
        self.set_message(stdout)

    def handle_state(self, state): #TODO: APAGAR DEPOIS
        states = {
            QtCore.QProcess.NotRunning: 'Not running',
            QtCore.QProcess.Starting: 'Starting',
            QtCore.QProcess.Running: 'Running',
        }
        state_name = states[state]
        print(state_name)

    # ------- functions -------
    def evt_update_progress(self, val):
        self.pbar.setValue(val)

    def set_message(self, s):
        self.text.appendPlainText(s)

    def evt_set_logger(self, outeriter, f0val, fval, xval, natural_freqs):
        outvector1 = str([outeriter, round(f0val, 6), np.round(fval.flatten(),6)])
        outvector2 = str(xval.flatten())
        if natural_freqs is not None:
            outvector3 = str(natural_freqs)
        else:
            outvector3 = None
        return outvector1, outvector2, outvector3

    def stop_execution(self):   
        if self.process_opt is not None:
            self.process_opt.kill()
            self.finish_process()

    def finish_process(self):
        self.process_opt = None
        #self.delete_temp()            
        self.evt_update_progress(100)
        self.button_stop.setEnabled(False)
        self.button_run.setEnabled(True)

if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    widget = QtWidgets.QStackedWidget()
    mainwindow = MainWindow()
    widget.addWidget(mainwindow)
    widget.showMaximized()
    sys.exit(app.exec_())