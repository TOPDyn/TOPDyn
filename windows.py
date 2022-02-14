from PyQt5 import QtWidgets, QtCore, QtGui
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.pyplot as plt
plt.switch_backend('Agg')
import pyqtgraph as pg

import os
import sys
import re
import ast
import shutil
import numpy as np

import plot_grid as mf
from plots_2d import PlotsFem2d
from plots_opt import PlotOpt
from parameters import ParametersOpt, ParametersFEM2D
from param_texts import TextBc, TextConstraint, TextOpt, TextFem2d

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
        self.setLayout(layout)
       
    def go_to_fem2d(self):
        param_2d = ParametersFEM2D()
        text_main = TextFem2d()
        fem2d = WindowsFem2d(param_2d, text_main)
        widget.addWidget(fem2d)
        widget.setCurrentIndex(widget.currentIndex()+1)

    def go_to_opt(self):
        param_opt = ParametersOpt()
        text_main = TextOpt()
        opt = WindowsOptimization(param_opt, text_main)
        widget.addWidget(opt)
        widget.setCurrentIndex(widget.currentIndex()+1)

class SecondWindow(QtWidgets.QDialog):
    def __init__(self, param, text_main):
        super(SecondWindow, self).__init__()
        self._counter_button_run = 0

        self.pbar = QtWidgets.QProgressBar()
        # lateral menu
        left_layout = QtWidgets.QFormLayout()
        self.btn_back_menu = QtWidgets.QPushButton('Back to Menu')
        self.btn_back_menu.clicked.connect(self.go_to_screen1)
        left_layout.addRow(self.btn_back_menu)
        self.param = param
        self.param.add_btns(left_layout)
        self.btn_to_load = QtWidgets.QPushButton('Next')
        self.btn_to_load.clicked.connect(self.go_to_load)
        left_layout.addRow(self.btn_to_load)
    
        left_layout.setSpacing(20)
        self.left_widget = QtWidgets.QWidget()
        self.left_widget.setLayout(left_layout)
        scroll = QtWidgets.QScrollArea()
        #scroll.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)
        #scroll.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)
        scroll.setWidget(self.left_widget)
        scroll.setWidgetResizable(True)

        # window 
        self.right_widget = QtWidgets.QTabWidget()       
        self.tab1 = self.ui1(text_main)
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
    def ui1(self, text_main):
        layout_ui1 = QtWidgets.QVBoxLayout()
        layout_ui1.addWidget(text_main.editor)
        layout_ui1.addStretch(20)
        
        main = QtWidgets.QWidget()
        main.setLayout(layout_ui1)
        return main

    def upd_ui1_layout(self, new_layout):
        if self.tab1.layout() is not None:
            QtWidgets.QWidget().setLayout(self.tab1.layout())
        self.tab1.setLayout(new_layout)
    
    def ui2(self):
        self.init_ui2_layout = QtWidgets.QFormLayout()
        main = QtWidgets.QWidget()
        main.setLayout(self.init_ui2_layout)
        return main

    def upd_ui2_layout(self, new_layout):
        if self.tab2.layout() is not None:
            QtWidgets.QWidget().setLayout(self.tab2.layout())
        self.tab2.setLayout(new_layout)

    def upd_left_layout(self, new_layout):
        if self.left_widget.layout() is not None:
            QtWidgets.QWidget().setLayout(self.left_widget.layout())
        self.left_widget.setLayout(new_layout)

    def ui_add_load(self, update=False):
        self.add_load_text()
        ui_layout = QtWidgets.QFormLayout()
        ui_layout.setVerticalSpacing(10)
        self.btn_add_new_load = QtWidgets.QPushButton('Add new load')
        self.btn_add_new_load.clicked.connect(self.add_new_load)
        ui_layout.addRow(self.btn_add_new_load)

        self.btn_reset_load = QtWidgets.QPushButton('Reset Load')
        self.btn_reset_load.clicked.connect(self.reset_load)
        ui_layout.addRow(self.btn_reset_load)

        self.btn_back_param = QtWidgets.QPushButton('Back')
        self.btn_back_param.clicked.connect(self.back_to_param)
        ui_layout.addRow(self.btn_back_param)

        self.btn_to_node_constrain = QtWidgets.QPushButton('Next')
        self.btn_to_node_constrain.clicked.connect(self.go_to_node_constrain)
        ui_layout.addRow(self.btn_to_node_constrain)

        if update:
            self.param.load.rewrite(ui_layout)
        else:
            self.param.load.reset_list()
            self.param.load.create_btn()
            self.param.load.set_default()
            self.param.load.add_btn(ui_layout)
        self.upd_left_layout(ui_layout)

    def ui_add_node_constrain(self):
        self.add_node_constrain_text()

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
            self.ui_add_load()

    def go_to_node_constrain(self):
        self.param.load.check_btn()
        if len(self.param.load.warnings) != 0:
            dlg = PopUp(self.param.load.warnings)      
            dlg.exec_()
        else:
            self.param.load.check_values()
            if len(self.param.load.warnings) != 0:
                dlg = PopUp(self.param.load.warnings)      
                dlg.exec_()
            else:
                self.param.load.update()
                self.ui_add_node_constrain()

    def back_to_param(self):
        left_layout = QtWidgets.QFormLayout()
        left_layout.setVerticalSpacing(20)
        self.btn_back_menu = QtWidgets.QPushButton('Back to Menu')
        self.btn_back_menu.clicked.connect(self.go_to_screen1)
        left_layout.addRow(self.btn_back_menu)
        
        self.param.create_btns()
        self.param.update_default()
        self.param.add_btns(left_layout)

        self.btn_to_load = QtWidgets.QPushButton('Next')
        self.btn_to_load.clicked.connect(self.go_to_load)
        left_layout.addRow(self.btn_to_load)
        self.upd_left_layout(left_layout)

    def back_to_load(self):
        self.ui_add_load(update=True)

    def add_new_load(self):
        self.param.load.create_btn()
        self.param.load.set_default()
        self.param.load.add_btn(self.left_widget.layout())

    def add_new_node_constrain(self):
        self.param.node_constrain.create_btn()
        self.param.node_constrain.set_default()
        self.param.node_constrain.add_btn(self.left_widget.layout())

    def add_load_text(self):
        ui_layout = QtWidgets.QFormLayout()
        self.bc_text = TextBc()
        ui_layout.addRow(self.bc_text.editor_load)
        self.upd_ui1_layout(ui_layout)

    def add_node_constrain_text(self):
        ui_layout = QtWidgets.QFormLayout()
        self.bc_text = TextBc()
        ui_layout.addRow(self.bc_text.editor_node_constrain)
        self.upd_ui1_layout(ui_layout)

    def reset_load(self):
        self.ui_add_load()

    def reset_node_constrain(self):
        self.ui_add_node_constrain()

    # ------- functions -------     
    def delete_temp(self):
        folder_name = 'temp'
        directory = os.path.join(os.path.dirname(__file__), folder_name)
        if os.path.exists(directory):
            shutil.rmtree(directory)

# Regular expressions to extract values from QProcess
progress_re = re.compile("Total complete: (\d+)%")

update_data_re = re.compile("Updating data")

plot_freqresp_re = re.compile("Ploting frequency response")

create_plot_re = re.compile("Creating plot")

mesh_deform_re = re.compile("Ploting deformed mesh")

save_figs_re = re.compile("Saving figures")

def simple_percent_parser(output):
    """ Matches lines using the progress_re regex,
    returning a single integer for the % progress. """
    m = progress_re.search(output)
    if m:
        pc_complete = m.group(1)
        return int(pc_complete)

def simple_update_parser(output):
    """ Matches lines using the update_data_re regex,
    returning a bool. """
    m = update_data_re.search(output)
    if m:
        return True

def simple_freq_parser(output):
    """ Matches lines using the plot_freqresp_re regex,
    returning a bool. """
    m = plot_freqresp_re.search(output)
    if m:
        return True

def simple_save_parser(output):
    """ Matches lines using the save_figs_re regex,
    returning a bool. """
    m = save_figs_re.search(output)
    if m:
        return True

def simple_create_parser(output):
    """ Matches lines using the create_plot_re regex,
    returning a bool. """
    m = create_plot_re.search(output)
    if m:
        return True

def simple_deform_parser(output):
    """ Matches lines using the mesh_deform_re regex,
    returning a bool. """
    m = mesh_deform_re.search(output)
    if m:
        return True

class WindowsFem2d(SecondWindow):
    def __init__(self, param, text_main):
        super().__init__(param, text_main)

        self.stop_thread = False

    # ------- pages -------
    def add_canvas_ui2(self): 
        ui2_layout = QtWidgets.QVBoxLayout()
        ui2_layout.addWidget(self.pbar)

        self.canvas_mesh = FigureCanvas(plt.Figure())
        self.canvas_mesh.figure.tight_layout()
        ui2_layout.addWidget(self.canvas_mesh)

        if self.param.freqrsp:
            self.canvas_freq = FigureCanvas(plt.Figure())
            self.canvas_freq.figure.tight_layout()
            ui2_layout.addWidget(self.canvas_freq)
        self.upd_ui2_layout(ui2_layout)

    def add_fem_text(self):
        ui_layout = QtWidgets.QFormLayout()
        self.fem_text = TextFem2d()
        ui_layout.addRow(self.fem_text.editor)
        self.upd_ui1_layout(ui_layout)

    def ui_add_node_constrain(self, update=False):
        super().ui_add_node_constrain()

        ui_layout = QtWidgets.QFormLayout()
        self.btn_back_load = QtWidgets.QPushButton('Back')
        self.btn_back_load.clicked.connect(self.back_to_load)
        ui_layout.addRow(self.btn_back_load)

        self.button_run = QtWidgets.QPushButton('Run')
        self.button_run.clicked.connect(lambda: self.run())
        ui_layout.addRow(self.button_run)

        self.button_stop = QtWidgets.QPushButton('Stop')
        self.button_stop.setEnabled(False)
        ui_layout.addRow(self.button_stop)

        self.btn_add_new_node_constrain = QtWidgets.QPushButton('Constrain new node displacement')
        self.btn_add_new_node_constrain.clicked.connect(self.add_new_node_constrain)
        ui_layout.addRow(self.btn_add_new_node_constrain)

        self.btn_reset_node_constrain = QtWidgets.QPushButton('Reset constraint node displacement')
        self.btn_reset_node_constrain.clicked.connect(self.reset_node_constrain)
        ui_layout.addRow(self.btn_reset_node_constrain)
        if update:
            self.param.node_constrain.rewrite(ui_layout)
        else:
            self.param.node_constrain.reset_list()
            self.param.node_constrain.create_btn()
            self.param.node_constrain.set_default()
            self.param.node_constrain.add_btn(ui_layout)
        self.upd_left_layout(ui_layout)

    # ------- buttons -------
    def back_to_param(self):
        super().back_to_param()
        self.add_fem_text()
     
    def run(self):
        self.param.node_constrain.check_btn()
        if len(self.param.node_constrain.warnings) != 0:
            dlg = PopUp(self.param.node_constrain.warnings)      
            dlg.exec_()
        else:
            self.param.node_constrain.check_values()
            if len(self.param.node_constrain.warnings) != 0:
                dlg = PopUp(self.param.node_constrain.warnings)      
                dlg.exec_()
            else:
                self.button_stop.setEnabled(True)
                self.button_run.setEnabled(False)
                self.btn_back_load.setEnabled(False)
        
                self.param.node_constrain.update()
                
                self._counter_button_run += 1
                if self._counter_button_run == 1:
                    self.right_widget.setCurrentIndex(1)

                self.add_canvas_ui2()
                self.param.create_dict_param()
                self.param.export_param() 
                self.process_fem = QtCore.QProcess()
                self.process_fem.readyReadStandardError.connect(self.handle_stderr)
                self.process_fem.stateChanged.connect(self.handle_state)
                self.process_fem.finished.connect(self.plot_graphs)
                self.process_fem.start("python", ['process_fem2d.py'])
                self.button_stop.clicked.connect(lambda: self.stop_execution())

    def handle_stderr(self):
        data = self.process_fem.readAllStandardError()
        stderr = bytes(data).decode("utf8")
        # Extract progress if it is in the data.
        progress = simple_percent_parser(stderr)
        if progress:
            self.evt_update_progress(progress)

    def handle_state(self, state):
        states = {
            QtCore.QProcess.NotRunning: 'Not running',
            QtCore.QProcess.Starting: 'Starting',
            QtCore.QProcess.Running: 'Running',
        }
        state_name = states[state]
        print(state_name)

    # ------- functions -------
    def evt_update_mesh(self):
        self.canvas_mesh.draw()

    def evt_update_freq(self):
        self.canvas_freq.draw()

    def evt_update_progress(self, val):
        self.pbar.setValue(val)

    def stop_execution(self):
        self.evt_update_progress(0)
        if self.process_fem is not None:
            self.process_fem.kill()
            self.process_fem = None
            self.delete_temp()
            self.button_stop.setEnabled(False)
            self.button_run.setEnabled(True)
            self.btn_back_load.setEnabled(True)
        else:
            self.stop_thread = True

    def plot_graphs(self):
        if self.process_fem is not None:
            self.worker_plot = WorkerPlot(parent=self)
            self.worker_plot.start()
            self.worker_plot.complete_worker.connect(lambda: self.worker_plot.terminate())
            self.worker_plot.complete_worker.connect(lambda: self.worker_plot.deleteLater())
            self.worker_plot.complete_worker.connect(lambda: self.button_stop.setEnabled(False))
            self.worker_plot.complete_worker.connect(lambda: self.button_run.setEnabled(True))
            self.worker_plot.complete_worker.connect(lambda: self.btn_back_load.setEnabled(True))
            self.worker_plot.complete_worker.connect(self.delete_temp)
            self.worker_plot.update_mesh.connect(self.evt_update_mesh)
            self.worker_plot.update_freq.connect(self.evt_update_freq)
            self.worker_plot.update_progess.connect(self.evt_update_progress)
            self.process_fem = None

class WorkerPlot(QtCore.QThread):
    update_mesh = QtCore.pyqtSignal(bool)
    update_freq = QtCore.pyqtSignal(bool)
    update_progess = QtCore.pyqtSignal(int)
    complete_worker = QtCore.pyqtSignal(bool)

    def __init__(self, parent):    
        super().__init__()
        self.parent = parent
               
    @QtCore.pyqtSlot()
    def run(self):
        while True:
            # READ
            folder_name = 'temp'
            directory = os.path.join(os.path.dirname(__file__), folder_name) 

            file = open(os.path.join(directory, 'mesh_data.txt'), "r")
            contents = file.read()
            mesh_data = ast.literal_eval(contents)
            file.close()

            coord = np.loadtxt(os.path.join(directory, 'coord.txt'))
            connect = np.loadtxt(os.path.join(directory, 'connect.txt'), dtype=int)
            disp_vector_org = np.loadtxt(os.path.join(directory, 'disp_vector.txt'), dtype=np.complex)

            load_matrix = np.loadtxt(os.path.join(directory, 'load_matrix.txt')).reshape(-1, 4) #TODO: VER A QUESTAO DO VALOR DA FORÇA NEGATIVO
            constr_matrix = np.loadtxt(os.path.join(directory, 'constr_matrix.txt'), dtype=int)

            # Deformed mesh
            plot_2d = PlotsFem2d(coord, connect)
            disp_vector = plot_2d.change_disp_shape(disp_vector_org.real)

            plot_2d.apply_disp(disp_vector, self.parent.param.factor)
            plot_2d.build_collection()

            if self.parent.stop_thread:
                self.parent.stop_thread = False
                break

            self.parent.ax_mesh = self.parent.canvas_mesh.figure.subplots()
            plot_2d.plot_collection(self.parent.ax_mesh, mesh_data["lx"], mesh_data["ly"], load_matrix, constr_matrix)

            if self.parent.stop_thread:
                self.parent.stop_thread = False
                break
            
            self.update_mesh.emit(True)
            self.update_progess.emit(90)

            # Frequency response
            vector_U = None
            if self.parent.param.freqrsp:
                node_plot = np.loadtxt(os.path.join(directory, 'node_plot.txt'), dtype=int).reshape(1,3)
                vector_U = np.loadtxt(os.path.join(directory, 'vector_U.txt'), dtype=np.complex)
        
                self.parent.ax_freq = self.parent.canvas_freq.figure.subplots()
                plot_2d.plot_freq_rsp(self.parent.ax_freq, node_plot, self.parent.param.freq_range, vector_U)

                if self.parent.stop_thread:
                    self.parent.stop_thread = False
                    break

                self.update_freq.emit(True)
                self.update_progess.emit(95)
            
            if self.parent.param.save:
                data = os.path.join(os.path.join(os.path.dirname(__file__)), 'data_2d')
                os.makedirs(data, exist_ok=True)

                np.savetxt(os.path.join(data, 'disp_vector.txt'), disp_vector_org)
                if vector_U is not None:
                    np.savetxt(os.path.join(data, 'freqrsp_vector.txt'), vector_U)

                self.parent.canvas_mesh.figure.savefig(os.path.join(data, 'mesh.png'))
                
                if self.parent.param.freqrsp:
                    self.parent.canvas_freq.figure.savefig(os.path.join(data, 'frequency_response.png'))

            self.update_progess.emit(100)
            break

        self.complete_worker.emit(True)

class WindowsOptimization(SecondWindow):
    def __init__(self, param, text_main):
        super().__init__(param, text_main)
        self.stop_thread = False
        self.directory = os.path.join(os.path.dirname(__file__))
        self.dir_temp = os.path.join(self.directory, 'temp')

    # ------- param windows -------
    def add_opt_text(self):
        ui_layout = QtWidgets.QFormLayout()
        self.text_main = TextOpt()
        ui_layout.addRow(self.text_main.editor)
        self.upd_ui1_layout(ui_layout)

    def add_contraint_text(self):
        ui_layout = QtWidgets.QFormLayout()
        constraint_text = TextConstraint()
        ui_layout.addRow(constraint_text.editor)
        self.upd_ui1_layout(ui_layout)

    def ui_add_node_constrain(self, update=False):
        super().ui_add_node_constrain()

        ui_layout = QtWidgets.QFormLayout()
        ui_layout.setVerticalSpacing(10)
        self.btn_add_new_node_constrain = QtWidgets.QPushButton('Constrain new node displacement')
        self.btn_add_new_node_constrain.clicked.connect(self.add_new_node_constrain)
        ui_layout.addRow(self.btn_add_new_node_constrain)

        self.btn_reset_node_constrain = QtWidgets.QPushButton('Reset constraint node displacement')
        self.btn_reset_node_constrain.clicked.connect(self.reset_node_constrain)
        ui_layout.addRow(self.btn_reset_node_constrain)

        self.btn_back_load = QtWidgets.QPushButton('Back')
        self.btn_back_load.clicked.connect(self.back_to_load)
        ui_layout.addRow(self.btn_back_load)

        self.btn_to_constraint = QtWidgets.QPushButton('Next')
        self.btn_to_constraint.clicked.connect(self.go_to_constraint)
        ui_layout.addRow(self.btn_to_constraint)
        
        if update:
            self.param.node_constrain.rewrite(ui_layout)
        else:
            self.param.node_constrain.reset_list()
            self.param.node_constrain.create_btn()
            self.param.node_constrain.set_default()
            self.param.node_constrain.add_btn(ui_layout)
        self.upd_left_layout(ui_layout)

    def ui_add_constraint(self):
        self.add_contraint_text()

        ui_layout = QtWidgets.QFormLayout()
        ui_layout.setVerticalSpacing(15)
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
    def go_to_constraint(self):
        self.param.node_constrain.check_btn()
        if len(self.param.node_constrain.warnings) != 0:
            dlg = PopUp(self.param.node_constrain.warnings)      
            dlg.exec_()
        else:
            self.param.node_constrain.check_values()
            if len(self.param.node_constrain.warnings) != 0:
                dlg = PopUp(self.param.node_constrain.warnings)      
                dlg.exec_()
            else:
                self.param.node_constrain.update()
                self.ui_add_constraint()
    
    def back_to_param(self):
        super().back_to_param()
        self.add_opt_text()

    def back_to_node_constrain(self):
        self.ui_add_node_constrain(update=True)

    def run(self):
        self.param.check_constraint()
        if len(self.param.warnings_constraint) != 0:
            dlg = PopUp(self.param.warnings_constraint)      
            dlg.exec_()
        else:
            self.param.check_constraint_values()
            if len(self.param.warnings_constraint) != 0:
                dlg = PopUp(self.param.warnings_constraint)      
                dlg.exec_()
            else:
                self.button_stop.setEnabled(True)
                self.button_run.setEnabled(False)
                self.btn_back_constrain.setEnabled(False)
        
                self.param.update_constraint()
                self.param.constraint_to_list()
                        
                self._counter_button_run += 1
                if self._counter_button_run == 1:
                    self.right_widget.setCurrentIndex(1)
                    print("entrou")

                pg.setConfigOption('background', 'w')
                pg.setConfigOption('foreground', 'k')

                self.add_canvas_ui2()
                self.param.create_dict_param()
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

    # ------- plot windows -------
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
        self.graph_conv.setLabel('left', "normalized functions")
        self.graph_conv.setLabel('bottom', "iteration")
        ui2_layout.addWidget(self.graph_conv)

        self.text = QtWidgets.QPlainTextEdit()
        self.text.setReadOnly(True)
        ui2_layout.addWidget(self.text)
        self.upd_ui2_layout(ui2_layout)

    def add_freq_to_canvas_ui2(self): 
        ui2_layout = QtWidgets.QVBoxLayout()
        ui2_layout.addWidget(self.pbar)
        ui2_layout.addWidget(self.graph_grid)
        self.graph_grid.setFixedWidth(self.grid_width)
        self.graph_grid.setFixedHeight(self.grid_height)
        ui2_layout.addWidget(self.graph_conv)
        self.graph_conv.setFixedWidth(self.conv_width)
        self.graph_conv.setFixedHeight(self.conv_height)
        ui2_layout.addWidget(self.graph_freq)  
        self.graph_freq.setFixedWidth(self.conv_width)
        self.graph_freq.setFixedHeight(self.conv_height)   
        ui2_layout.addWidget(self.text)
        self.upd_ui2_layout(ui2_layout)

    def add_deformed_mesh_to_canvas_ui2(self): 
        ui2_layout = QtWidgets.QVBoxLayout()
        ui2_layout.addWidget(self.pbar)
        ui2_layout.addWidget(self.graph_grid)
        self.graph_grid.setFixedWidth(self.grid_width)
        self.graph_grid.setFixedHeight(self.grid_height)
        ui2_layout.addWidget(self.graph_conv)
        self.graph_conv.setFixedWidth(self.conv_width)
        self.graph_conv.setFixedHeight(self.conv_height)
        if self.param.freqrsp:
            ui2_layout.addWidget(self.graph_freq)
            self.graph_freq.setFixedWidth(self.conv_width)
            self.graph_freq.setFixedHeight(self.conv_height)  
        ui2_layout.addWidget(self.canvas_mesh)
        self.canvas_mesh.setFixedWidth(self.graph_grid.width())
        self.canvas_mesh.setFixedHeight(self.graph_grid.height())
        ui2_layout.addWidget(self.text)
        self.text.setFixedWidth(self.text_width)
        self.text.setFixedHeight(self.text_height)
        self.upd_ui2_layout(ui2_layout)

    # ------- handle process output -------
    def handle_stderr(self):
        data = self.process_opt.readAllStandardError()
        stderr = bytes(data).decode("utf8")
        # Extract progress if it is in the data.
        progress = simple_percent_parser(stderr)
        if progress:
            self.evt_update_progress(progress)
            self.set_message(stderr)
        
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

            # keep standard size
            if outeriter == 1:
                self.conv_height = self.graph_conv.height()
                self.conv_width = self.graph_conv.width()

                self.grid_height = self.graph_grid.height()
                self.grid_width = self.graph_grid.width()

                self.text_height = self.text.height()
                self.text_width = self.text.width() 

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
            plot_mesh.apply_disp(disp_vector, self.param.factor)

            if (self.param.passive_type is not None) and (self.param.passive_type==0):
                passive_el = True
                passive_coord = self.param.passive_coord
            else:
                passive_el = False
                passive_coord = None
            plot_mesh.build_collection(passive_coord)
        
            self.canvas_mesh = FigureCanvas(plt.Figure())
            self.canvas_mesh.figure.tight_layout()
            ax_mesh = self.canvas_mesh.figure.subplots()
            plot_mesh.plot_collection(ax_mesh, self.plot_opt.lx, self.plot_opt.ly, 
                                        load_matrix, constr_matrix, passive_el)

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

    def handle_stdout(self):
        data = self.process_opt.readAllStandardOutput()
        stdout = bytes(data).decode("utf8")
        self.set_message(stdout)

    def handle_state(self, state):
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
        self.delete_temp()            
        self.evt_update_progress(100)
        self.button_stop.setEnabled(False)
        self.button_run.setEnabled(True)
        self.btn_back_constrain.setEnabled(True)

if __name__ == '__main__':
    app = QtWidgets.QApplication(sys.argv)
    widget = QtWidgets.QStackedWidget()
    mainwindow = MainWindow()
    widget.addWidget(mainwindow)
    widget.showMaximized()
    sys.exit(app.exec_())