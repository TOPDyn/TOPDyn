from PyQt5 import QtCore, QtWidgets

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
            <p> <b> <font size="+1"> Penal. Stiffness </b> <font size="+1"> (int): Penalization power to stiffness.  </font> </p>
            <p> <b> <font size="+1"> Penal. Mass </b> <font size="+1"> (int): Penalization power to mass.  </font> </p>
            <p> <b> <font size="+1"> Alpha </b> <font size="+1"> (float): Damping coefficient proportional to mass.   </font> </p>
            <p> <b> <font size="+1"> Beta </b> <font size="+1"> (float): Damping coefficient proportional to stiffness.   </font> </p>
            <p> <b> <font size="+1"> Eta </b> <font size="+1"> (float): Damping coefficient.   </font> </p>
            <p> <b> <font size="+1"> Passive Elements </b> <font size="+1"> (tuple): Region that the shape will not be changed.   </font> </p>
                <p style="margin-left:2em"> <font size="+1"> <b> Passive coordinates </b>: Passive element coordinates. </font> </p>
                    <p style="margin-left:4em"> <font size="+1"> Example: ((0.5, 1), (0.3, 0.6)) = ((x_initial, x_final), (y_initial, y_final))</font> </p> 
                <p style="margin-left:2em"> <font size="+1"> - If 0 is selected the elements in the region are disregarded. </font> </p>
                <p style="margin-left:2em"> <font size="+1"> - If 1 is selected the elements in the region are static. </font> </p>
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
                <p style="margin-left:2em"> <font size="+1"> - Relative density vector, objective functions, constraint functions and frequency response values.  </font> </p>
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