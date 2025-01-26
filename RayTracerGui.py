#!/usr/bin/env python
import sys
import numpy as np
import pickle  # Use 'pickle' instead of 'cPickle' in Python 3
from PyQt5.QtWidgets import (
    QWidget, QLabel, QLineEdit, QTextEdit, QGridLayout, QApplication,
    QPushButton, QComboBox, QFileDialog, QSpinBox, QCheckBox
)
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

import RayTrace
import copy as cp
import ExternalScript as ExScr
import RTUtils as RTU

class Example(QWidget):

    def __init__(self):
        super(Example, self).__init__()

        self.initUI()
        self.ClearAll()

        self.colors = ['blue', 'yellow', 'red', 'black']
        self.alphas = [0.05, 0.1, 0.1, 0.2]

    def initUI(self):
        ellabel = QLabel('Elements')
        title = QLabel('Add Element')
        Type = QLabel('Element Type')
        Geo = QLabel('Element Geometry')
        Position = QLabel('Midpoint of Element (x,y,z)')
        Direction = QLabel('Direction of Element Face (x,y,z)')
        Size = QLabel('Element Size (a,b) along axis or (a,b,[x,y,z]) to specify axis')
        Points = QLabel('Specify points (x1,y1,z1,x2,y2...) triangle only')
        Material1 = QLabel('Material 1 (Interface)')
        Material2 = QLabel('Material 2 (Interface)')
        Focal = QLabel('Focal Length')
        DetSpec = QLabel('Spectrum for EQE')

        self.ElementCombo = QComboBox()
        self.ElementCombo.currentTextChanged.connect(self.FillSelectedElement)
        self.TypeEdit = QComboBox()
        self.TypeEdit.addItem("Mirror")
        self.TypeEdit.addItem("Interface")
        self.TypeEdit.addItem("Thin Lens")
        self.TypeEdit.addItem("Dump")
        self.TypeEdit.currentTextChanged.connect(self.TypeEditDisabler)

        self.GeoEdit = QComboBox()
        self.GeoEdit.addItem("Plane")
        self.GeoEdit.addItem("Circle")
        self.GeoEdit.addItem("Cylinder")
        self.GeoEdit.addItem("Triangle")
        self.GeoEdit.addItem("Box")
        self.GeoEdit.currentTextChanged.connect(self.GeoEditDisabler)

        self.PositionEdit = QLineEdit()
        self.DirectionEdit = QLineEdit()
        self.SizeEdit = QLineEdit()
        self.PointsEdit = QLineEdit()
        self.MaterialEdit1 = QComboBox()
        self.MaterialEdit2 = QComboBox()
        self.FocalLengthEdit = QLineEdit()
        self.DetSpecEdit = QLineEdit()

        self.AddBtn = QPushButton("Accept")
        self.AddBtn.clicked.connect(self.AcceptElement)
        NewElBtn = QPushButton("New")
        NewElBtn.clicked.connect(self.EnableNew)
        DltElBtn = QPushButton("Delete")
        DltElBtn.clicked.connect(self.DeleteElement)
        EditBtn = QPushButton("Edit")
        EditBtn.clicked.connect(self.EnableEdit)
        SaveBtn = QPushButton("Save")
        SaveBtn.clicked.connect(self.SaveElements)
        LoadBtn = QPushButton("Load")
        LoadBtn.clicked.connect(self.LoadElements)
        ActionBtn = QPushButton("Execute External Script")
        ActionBtn.clicked.connect(self.ExtScript)
        # Disable editing elements at start
        self.DisableEdit()

        self.ElList = QTextEdit()
        self.ElPropList = QTextEdit()
        self.ElSizeList = QTextEdit()
        self.ElGeoList = QTextEdit()

        # Lay out the element editing widgets
        grid = QGridLayout()
        grid.setSpacing(2)

        grid.addWidget(ellabel, 0, 0)
        grid.addWidget(self.ElementCombo, 0, 1)
        grid.addWidget(NewElBtn, 1, 1)
        grid.addWidget(EditBtn, 0, 2)
        grid.addWidget(DltElBtn, 0, 3)

        grid.addWidget(self.ElList, 1, 2, 10, 1)
        grid.addWidget(self.ElPropList, 1, 3, 10, 1)
        grid.addWidget(self.ElSizeList, 1, 4, 10, 1)
        grid.addWidget(self.ElGeoList, 1, 5, 10, 1)

        grid.addWidget(title, 3, 0)

        grid.addWidget(Type, 4, 0)
        grid.addWidget(self.TypeEdit, 4, 1)

        grid.addWidget(Geo, 5, 0)
        grid.addWidget(self.GeoEdit, 5, 1)

        grid.addWidget(Position, 6, 0)
        grid.addWidget(self.PositionEdit, 6, 1)

        grid.addWidget(Direction, 7, 0)
        grid.addWidget(self.DirectionEdit, 7, 1)

        grid.addWidget(Size, 8, 0)
        grid.addWidget(self.SizeEdit, 8, 1)

        grid.addWidget(Points, 9, 0)
        grid.addWidget(self.PointsEdit, 9, 1)

        grid.addWidget(Material1, 10, 0)
        grid.addWidget(self.MaterialEdit1, 10, 1)

        grid.addWidget(Material2, 11, 0)
        grid.addWidget(self.MaterialEdit2, 11, 1)

        grid.addWidget(Focal, 12, 0)
        grid.addWidget(self.FocalLengthEdit, 12, 1)

        grid.addWidget(DetSpec, 13, 0)
        grid.addWidget(self.DetSpecEdit, 13, 1)

        grid.addWidget(self.AddBtn, 14, 0)
        grid.addWidget(SaveBtn, 14, 1)
        grid.addWidget(LoadBtn, 14, 2)
        grid.addWidget(ActionBtn, 14, 3)

        # ---------- Material section ----------
        MaterialTitle = QLabel('Add Material')
        Name = QLabel('Material Name')
        Ref = QLabel('Refractive Index')
        Absorb = QLabel('Is Absorber ?')
        ConcD = QLabel('Concentration Donor')
        ConcA = QLabel('Concentration Acceptor')
        color = QLabel('Color')
        Aligned = QLabel('Is Aligned ?')
        AlignV = QLabel('Direction of alignment (x,y,z)')
        AlignS = QLabel('Std. Dev. of alignment (rad)')
        AbsEmAngle = QLabel('Angle between Abs/Emission Dipole (°)')

        FixDist = QLabel('Fixed Distance D/A?')
        FixDistR = QLabel('Distance (nm)')
        FixDistNeighbors = QLabel('Number of Donors per Acceptor')

        r_0 = QLabel('FRET Radius R0 (nm)')
        Q_D = QLabel('Fl. Quantum yield (donor)')
        Q_A = QLabel('Fl. Quantum yield (acceptor)')
        tau_D = QLabel('Lifetime (donor) ns')
        tau_A = QLabel('Lifetime (acceptor) ns')
        emax_D = QLabel('Max extinction (donor) (M^-1 cm^-1)')
        emax_A = QLabel('Max extinction (acceptor) (M^-1 cm^-1)')
        sp_abs_D = QLabel('Absorption spectrum Donor path')
        sp_em_D = QLabel('Emission spectrum Donor path')
        sp_abs_A = QLabel('Absorption spectrum Acceptor path')
        sp_em_A = QLabel('Emission spectrum Acceptor path')

        self.MaterialsCombo = QComboBox()
        self.MaterialsCombo.currentTextChanged.connect(self.FillSelectedMaterial)
        self.NameEdit = QLineEdit()
        self.RefEdit = QLineEdit()
        self.AbsorbEdit = QComboBox()
        self.AbsorbEdit.addItem("False")
        self.AbsorbEdit.addItem("True")

        self.ConcAEdit = QLineEdit()
        self.ConcDEdit = QLineEdit()
        self.colorEdit = QLineEdit()
        self.AlignedEdit = QComboBox()
        self.AlignedEdit.addItem("False")
        self.AlignedEdit.addItem("True")
        self.AlignVEdit = QLineEdit()
        self.AlignSEdit = QLineEdit()

        self.AbsEmAngleEdit = QSpinBox()
        self.AbsEmAngleEdit.setRange(0, 90)
        self.AbsEmAngleEdit.setValue(0)

        self.FixDistEdit = QComboBox()
        self.FixDistEdit.addItem("False")
        self.FixDistEdit.addItem("True")
        self.FixDistREdit = QLineEdit()
        self.FixDistNeighborsEdit = QSpinBox()
        self.FixDistNeighborsEdit.setRange(1, 10)
        self.FixDistNeighborsEdit.setValue(1)

        self.r_0Edit = QLineEdit()
        self.Q_DEdit = QLineEdit()
        self.Q_AEdit = QLineEdit()
        self.tau_DEdit = QLineEdit()
        self.tau_AEdit = QLineEdit()
        self.emax_DEdit = QLineEdit()
        self.emax_AEdit = QLineEdit()
        self.sp_abs_DEdit = QLineEdit()
        self.sp_em_DEdit = QLineEdit()
        self.sp_abs_AEdit = QLineEdit()
        self.sp_em_AEdit = QLineEdit()

        # Connect combos to disablers
        self.AbsorbEdit.currentTextChanged.connect(self.AbsorberDisabler)
        self.AlignedEdit.currentTextChanged.connect(self.AlignedDisabler)
        self.FixDistEdit.currentTextChanged.connect(self.FixDistDisabler)

        self.AddMatBtn = QPushButton("Accept")
        self.AddMatBtn.clicked.connect(self.AcceptMaterial)
        NewMatBtn = QPushButton("New")
        NewMatBtn.clicked.connect(self.NewMatEnabler)
        DltMatBtn = QPushButton("Delete")
        DltMatBtn.clicked.connect(self.DeleteMaterial)
        EditMatBtn = QPushButton("Edit")
        EditMatBtn.clicked.connect(self.EditMatFun)
        self.MatList = QTextEdit()

        grid.addWidget(MaterialTitle, 15, 0)
        grid.addWidget(self.MaterialsCombo, 15, 1)

        grid.addWidget(NewMatBtn, 16, 1)
        grid.addWidget(EditMatBtn, 15, 2)
        grid.addWidget(DltMatBtn, 15, 3)

        grid.addWidget(Name, 17, 0)
        grid.addWidget(self.NameEdit, 17, 1)

        grid.addWidget(Ref, 18, 0)
        grid.addWidget(self.RefEdit, 18, 1)

        grid.addWidget(color, 19, 0)
        grid.addWidget(self.colorEdit, 19, 1)

        grid.addWidget(Absorb, 20, 0)
        grid.addWidget(self.AbsorbEdit, 20, 1)

        grid.addWidget(ConcD, 21, 0)
        grid.addWidget(self.ConcDEdit, 21, 1)

        grid.addWidget(ConcA, 22, 0)
        grid.addWidget(self.ConcAEdit, 22, 1)

        grid.addWidget(Aligned, 24, 0)
        grid.addWidget(self.AlignedEdit, 24, 1)

        grid.addWidget(AlignV, 25, 0)
        grid.addWidget(self.AlignVEdit, 25, 1)

        grid.addWidget(AlignS, 26, 0)
        grid.addWidget(self.AlignSEdit, 26, 1)

        grid.addWidget(AbsEmAngle, 27, 0)
        grid.addWidget(self.AbsEmAngleEdit, 27, 1)

        grid.addWidget(FixDist, 28, 0)
        grid.addWidget(self.FixDistEdit, 28, 1)

        grid.addWidget(FixDistR, 29, 0)
        grid.addWidget(self.FixDistREdit, 29, 1)

        grid.addWidget(FixDistNeighbors, 30, 0)
        grid.addWidget(self.FixDistNeighborsEdit, 30, 1)

        grid.addWidget(r_0, 31, 0)
        grid.addWidget(self.r_0Edit, 31, 1)

        grid.addWidget(Q_D, 32, 0)
        grid.addWidget(self.Q_DEdit, 32, 1)

        grid.addWidget(Q_A, 33, 0)
        grid.addWidget(self.Q_AEdit, 33, 1)

        grid.addWidget(tau_D, 34, 0)
        grid.addWidget(self.tau_DEdit, 34, 1)

        grid.addWidget(tau_A, 35, 0)
        grid.addWidget(self.tau_AEdit, 35, 1)

        grid.addWidget(emax_D, 36, 0)
        grid.addWidget(self.emax_DEdit, 36, 1)

        grid.addWidget(emax_A, 37, 0)
        grid.addWidget(self.emax_AEdit, 37, 1)

        grid.addWidget(sp_abs_D, 38, 0)
        grid.addWidget(self.sp_abs_DEdit, 38, 1)

        grid.addWidget(sp_em_D, 39, 0)
        grid.addWidget(self.sp_em_DEdit, 39, 1)

        grid.addWidget(sp_abs_A, 40, 0)
        grid.addWidget(self.sp_abs_AEdit, 40, 1)

        grid.addWidget(sp_em_A, 41, 0)
        grid.addWidget(self.sp_em_AEdit, 41, 1)

        grid.addWidget(self.AddMatBtn, 42, 0)
        self.AddMatBtn.setEnabled(False)
        grid.addWidget(self.MatList, 16, 2, 12, 1)

        ShowBtn = QPushButton("Show")
        ShowBtn.clicked.connect(self.DrawFigure)
        grid.addWidget(ShowBtn, 44, 0)

        RunBtn = QPushButton("Run")
        RunBtn.clicked.connect(self.RunRayTracing)
        grid.addWidget(RunBtn, 45, 0)

        self.MatEnabler(False)
        self.ClearMaterials("0")

        # ---------- Photon generation ----------
        PhLabel = QLabel("Photon Generation")
        NumPh = QLabel("Number of Photons")
        StartPh = QLabel("Starting Point Photon (x,y,z)")
        DirectionPh = QLabel("Direction of Emission")
        EmDPh = QLabel("Emission Direction Distribution")
        SpectWavel = QLabel("Emission Wavelength Type")
        WavelPh = QLabel("Wavelength Photon (nm)")
        StMatPh = QLabel("Starting Material")
        EmPol = QLabel('Polarize emitted light')
        ScaleFactor = QLabel("Scale Factor")

        self.NumPhEdit = QSpinBox()
        self.NumPhEdit.setRange(1, 100000000)
        self.NumPhEdit.setValue(100)
        self.StartPhEdit = QLineEdit("0,0,0")
        self.EmDPhEdit = QComboBox()
        self.EmDPhEdit.addItem("Directional Isotropic")
        self.EmDPhEdit.addItem("Directional Strict")

        self.DirectionPhEdit = QLineEdit("0,0,-1")
        self.SpectWavel = QComboBox()
        self.SpectWavel.addItem("Sun spectrum AG15g")
        self.SpectWavel.addItem("Specified Wavelength...")
        self.WavelPhEdit = QSpinBox()
        self.WavelPhEdit.setRange(200, 2400)
        self.WavelPhEdit.setValue(450)
        self.StMatPhEdit = QComboBox()
        self.EmPolEdit = QCheckBox()
        self.EmPolEdit.setChecked(True)
        self.ScaleFactorEdit = QLineEdit("1.0")

        grid.addWidget(PhLabel, 16, 4)

        grid.addWidget(NumPh, 17, 4)
        grid.addWidget(self.NumPhEdit, 17, 5)

        grid.addWidget(StartPh, 18, 4)
        grid.addWidget(self.StartPhEdit, 18, 5)

        grid.addWidget(EmDPh, 19, 4)
        grid.addWidget(self.EmDPhEdit, 19, 5)

        grid.addWidget(DirectionPh, 20, 4)
        grid.addWidget(self.DirectionPhEdit, 20, 5)

        grid.addWidget(SpectWavel, 21, 4)
        grid.addWidget(self.SpectWavel, 21, 5)

        grid.addWidget(WavelPh, 22, 4)
        grid.addWidget(self.WavelPhEdit, 22, 5)

        grid.addWidget(StMatPh, 23, 4)
        grid.addWidget(self.StMatPhEdit, 23, 5)

        grid.addWidget(EmPol, 24, 4)
        grid.addWidget(self.EmPolEdit, 24, 5)

        grid.addWidget(ScaleFactor, 25, 4)
        grid.addWidget(self.ScaleFactorEdit, 25, 5)

        self.setLayout(grid)
        self.setGeometry(300, 300, 1100, 700)
        self.setWindowTitle('RayTracerGui')
        self.show()

    def ExtScript(self):
        # Execute external script code
        ExScr.main(self)

    def ClearAll(self):
        self.FirstElement = True
        self.Elements = np.zeros((1, 2, 3))
        self.ElementProperties = []
        self.ElementGeometry = []
        self.ElementSize = []
        self.Materials = []
        self.n_els = 0
        self.EditMode = False
        self.EditMatMode = False

        self.ElementCombo.clear()

        self.MaterialEdit1.clear()
        self.MaterialEdit2.clear()
        self.MaterialsCombo.clear()
        self.StMatPhEdit.clear()

    def DisableEdit(self):
        self.TypeEdit.setEnabled(False)
        self.GeoEdit.setEnabled(False)
        self.PositionEdit.setEnabled(False)
        self.DirectionEdit.setEnabled(False)
        self.SizeEdit.setEnabled(False)
        self.PointsEdit.setEnabled(False)
        self.MaterialEdit1.setEnabled(False)
        self.MaterialEdit2.setEnabled(False)
        self.FocalLengthEdit.setEnabled(False)
        self.DetSpecEdit.setEnabled(False)

    def EnableNew(self):
        self.EditMode = False
        self.TypeEdit.setEnabled(True)
        self.GeoEdit.setEnabled(True)
        self.PositionEdit.setEnabled(True)
        self.DirectionEdit.setEnabled(True)
        self.SizeEdit.setEnabled(True)
        self.AddBtn.setEnabled(True)
        self.GeoEditDisabler()
        self.TypeEditDisabler()

    def EnableEdit(self):
        self.EditMode = True
        self.TypeEdit.setEnabled(True)
        self.GeoEdit.setEnabled(True)
        self.PositionEdit.setEnabled(True)
        self.DirectionEdit.setEnabled(True)
        self.SizeEdit.setEnabled(True)
        self.AddBtn.setEnabled(True)

        self.GeoEditDisabler()
        self.TypeEditDisabler()

    def RunRayTracing(self):
        foldername = str(QFileDialog.getExistingDirectory(self, "Select Directory"))
        if not foldername:
            print("No directory selected.")
            return
        LSC = RayTrace.RayTrace()
        LSC.IOOperator.ParsePath(foldername)
        LSC.runSimulation()

    def SaveElements(self):
        foldername = str(QFileDialog.getExistingDirectory(self, "Select Directory"))
        if not foldername:
            print("No directory selected.")
            return

        with open(foldername + "/Elements.prt", 'wb') as f:
            pickle.dump(self.Elements, f)
        with open(foldername + "/ElementProperties.prt", 'wb') as f:
            pickle.dump(self.ElementProperties, f)
        with open(foldername + "/ElementGeometry.prt", 'wb') as f:
            pickle.dump(self.ElementGeometry, f)
        with open(foldername + "/ElementSize.prt", 'wb') as f:
            pickle.dump(self.ElementSize, f)
        with open(foldername + "/Materials.prt", 'wb') as f:
            pickle.dump(self.Materials, f)

        # Create and save array for photon generation
        try:
            start_ph = [float(x) for x in self.StartPhEdit.text().split(",")]
            dir_ph = [float(x) for x in self.DirectionPhEdit.text().split(",")]
            self.PhArray = [
                self.NumPhEdit.value(),
                start_ph,
                self.EmDPhEdit.currentText(),
                self.SpectWavel.currentText(),
                self.WavelPhEdit.value(),
                self.StMatPhEdit.currentText(),
                dir_ph,
                self.EmPolEdit.isChecked(),
                float(self.ScaleFactorEdit.text())
            ]
        except:
            print("Please enter all Photon Parameters properly.")
            return

        with open(foldername + "/Photon.prt", 'wb') as f:
            pickle.dump(self.PhArray, f)

    def LoadElements(self):
        self.ClearAll()

        foldername = str(QFileDialog.getExistingDirectory(self, "Select Directory"))
        if not foldername:
            print("No directory selected.")
            return

        types = {1: "Mirror", 2: "Interface", 3: "ThinLens", 4: "Dump"}
        geos = {1: "Plane", 2: "Circle", 3: "Cylinder", 4: "Triangle"}

        with open(foldername + "/Elements.prt", 'rb') as f:
            self.Elements = np.array(pickle.load(f))
        with open(foldername + "/ElementProperties.prt", 'rb') as f:
            self.ElementProperties = pickle.load(f)
        with open(foldername + "/ElementGeometry.prt", 'rb') as f:
            self.ElementGeometry = pickle.load(f)
        with open(foldername + "/ElementSize.prt", 'rb') as f:
            self.ElementSize = pickle.load(f)
        with open(foldername + "/Materials.prt", 'rb') as f:
            self.Materials = pickle.load(f)
        with open(foldername + "/Photon.prt", 'rb') as f:
            self.PhArray = pickle.load(f)

        for i in range(len(self.Elements)):
            element_str = (
                str(i + 1) + " " +
                types.get(self.ElementProperties[i][0], "Unknown") + " " +
                geos.get(self.ElementGeometry[i], "Unknown")
            )
            self.ElementCombo.addItem(element_str)

        for mat in self.Materials:
            self.MaterialEdit1.addItem(mat[0])
            self.MaterialEdit2.addItem(mat[0])
            self.MaterialsCombo.addItem(mat[0])
            self.StMatPhEdit.addItem(mat[0])

        self.n_els = len(self.Elements) + 1
        self.FirstElement = False
        self.FillSelectedElement()
        self.FillSelectedMaterial()
        self.FillPhoton()
        self.UpdateTextBoxes()

    def FillPhoton(self):
        """
        [self.NumPhEdit.value(), [start_x, start_y, start_z], emission_mode,
         wavelength_type, wavelength_val, start_material, [dir_x, dir_y, dir_z],
         pol_bool, scale_factor]
        """
        self.NumPhEdit.setValue(self.PhArray[0])
        sp = self.PhArray[1]
        self.StartPhEdit.setText(f"{sp[0]},{sp[1]},{sp[2]}")
        em_text = self.PhArray[2]
        idx_em = self.EmDPhEdit.findText(em_text)
        self.EmDPhEdit.setCurrentIndex(idx_em)
        w_text = self.PhArray[3]
        idx_sw = self.SpectWavel.findText(w_text)
        self.SpectWavel.setCurrentIndex(idx_sw)
        self.WavelPhEdit.setValue(self.PhArray[4])
        idx_mat = self.StMatPhEdit.findText(self.PhArray[5])
        self.StMatPhEdit.setCurrentIndex(idx_mat if idx_mat >= 0 else 0)
        dir_ph = self.PhArray[6]
        self.DirectionPhEdit.setText(f"{dir_ph[0]},{dir_ph[1]},{dir_ph[2]}")
        self.EmPolEdit.setChecked(self.PhArray[7])

        if len(self.PhArray) == 9:
            self.ScaleFactorEdit.setText(str(self.PhArray[8]))
        else:
            self.ScaleFactorEdit.setText("1.0")

    def AcceptElement(self):
        self.DisableEdit()
        self.AddBtn.setEnabled(False)

        if not self.EditMode:
            self.AddElement()
        else:
            self.EditElement()

        self.UpdateTextBoxes()

    def DeleteElement(self):
        if len(self.ElementGeometry) > 0:
            ind = self.ElementCombo.currentIndex()
            if ind < 0:
                return
            self.Elements = np.delete(self.Elements, ind, axis=0)
            self.n_els -= 1
            del self.ElementProperties[ind]
            del self.ElementGeometry[ind]
            del self.ElementSize[ind]
            self.ElementCombo.removeItem(ind)
            self.UpdateTextBoxes()

    def EditElement(self):
        ind = self.ElementCombo.currentIndex()
        t = self.TypeEdit.currentText()
        if t == "Mirror":
            proparray = [1]
        elif t == "Interface":
            proparray = [
                2,
                self.MaterialEdit1.currentText(),
                self.MaterialEdit2.currentText()
            ]
        elif t == "Thin Lens":
            proparray = [3, float(self.FocalLengthEdit.text())]
        elif t == "Dump":
            proparray = [4, 0, 0, self.DetSpecEdit.text()]
        else:
            proparray = [1]

        try:
            pos = [float(x) for x in self.PositionEdit.text().split(",")]
            direc = [float(x) for x in self.DirectionEdit.text().split(",")]
        except:
            print("Please enter valid (x,y,z) in Position or Direction.")
            return

        self.Elements[ind] = [pos, direc]
        self.ElementProperties[ind] = cp.deepcopy(proparray)
        self.ElementGeometry[ind] = self.GeoEdit.currentIndex() + 1

        shape_text = self.GeoEdit.currentText()
        if shape_text == "Plane":
            parsesize = self.SizeEdit.text().split(",")
            try:
                if len(parsesize) == 2:
                    self.ElementSize[ind] = [
                        float(parsesize[0]),
                        float(parsesize[1])
                    ]
                else:
                    # Example: "1, 2, [0, 1, 2]"
                    # The '...' might need parsing or cleaning
                    # We'll assume something like 1,2,[0,1,2]
                    # or 1,2,(0,1,2)
                    a0 = float(parsesize[0])
                    a1 = float(parsesize[1])
                    # next part is the vector
                    vec_str = ",".join(parsesize[2:])
                    vec_str = vec_str.strip("[]() ")
                    v_nums = [float(x) for x in vec_str.split(",")]
                    self.ElementSize[ind] = [a0, a1, v_nums]
            except:
                print("Please provide correct Plane size input.")
                return

        elif shape_text in ["Circle", "Cylinder"]:
            try:
                self.ElementSize[ind] = [float(self.SizeEdit.text())]
            except:
                print("Please provide correct radius for Circle/Cylinder.")
                return
        elif shape_text == "Triangle":
            # parse 9 floats
            try:
                parsesize = [float(x) for x in self.PointsEdit.text().split(",")]
                if len(parsesize) != 9:
                    print("Need 9 values (3 points) for Triangle.")
                    return
                self.ElementSize[ind] = [
                    [parsesize[0], parsesize[1], parsesize[2]],
                    [parsesize[3], parsesize[4], parsesize[5]],
                    [parsesize[6], parsesize[7], parsesize[8]]
                ]
                self.Elements[ind] = RTU.MakeTriangle(
                    [parsesize[0], parsesize[1], parsesize[2]],
                    [parsesize[3], parsesize[4], parsesize[5]],
                    [parsesize[6], parsesize[7], parsesize[8]]
                )
            except:
                print("Please provide correct coordinates for Triangle.")
                return

        self.ElementCombo.setItemText(
            ind,
            str(ind + 1) + " " + self.TypeEdit.currentText() + " " + shape_text
        )

    def AddElement(self):
        t = self.TypeEdit.currentText()
        if t == "Mirror":
            proparray = [1]
        elif t == "Interface":
            proparray = [
                2,
                self.MaterialEdit1.currentText(),
                self.MaterialEdit2.currentText()
            ]
        elif t == "Thin Lens":
            proparray = [3, float(self.FocalLengthEdit.text())]
        elif t == "Dump":
            proparray = [4, 0, 0, self.DetSpecEdit.text()]
        else:
            proparray = [1]

        shape_text = self.GeoEdit.currentText()

        try:
            pos = [float(x) for x in self.PositionEdit.text().split(",")]
            direc = [float(x) for x in self.DirectionEdit.text().split(",")]
        except:
            print("Please enter valid (x,y,z) for Position or Direction.")
            return

        # For the first element
        if self.FirstElement:
            self.Elements = np.array([ [pos, direc] ])
            self.ElementProperties = [cp.deepcopy(proparray)]
            self.FirstElement = False
        elif shape_text != "Triangle":
            self.Elements = np.append(self.Elements, [[pos, direc]], axis=0)
            self.ElementProperties.append(cp.deepcopy(proparray))

        # geometry
        g_idx = self.GeoEdit.currentIndex() + 1
        if shape_text != "Box":
            self.ElementGeometry.append(g_idx)

        if shape_text == "Plane":
            parsesize = self.SizeEdit.text().split(",")
            try:
                if len(parsesize) == 2:
                    self.ElementSize.append([
                        float(parsesize[0]),
                        float(parsesize[1])
                    ])
                else:
                    # parse something like "1,2,[0,1,2]"
                    a0 = float(parsesize[0])
                    a1 = float(parsesize[1])
                    vec_str = ",".join(parsesize[2:])
                    vec_str = vec_str.strip("[]() ")
                    v_nums = [float(x) for x in vec_str.split(",")]
                    self.ElementSize.append([a0, a1, v_nums])
            except:
                print("Error parsing plane size.")
                return

        elif shape_text in ["Circle", "Cylinder"]:
            try:
                r_val = float(self.SizeEdit.text())
                self.ElementSize.append([r_val])
            except:
                print("Error parsing radius.")
                return

        elif shape_text == "Triangle":
            try:
                parsesize = [float(x) for x in self.PointsEdit.text().split(",")]
                if len(parsesize) != 9:
                    print("Need 9 values for Triangle.")
                    return
                self.ElementSize.append([
                    [parsesize[0], parsesize[1], parsesize[2]],
                    [parsesize[3], parsesize[4], parsesize[5]],
                    [parsesize[6], parsesize[7], parsesize[8]]
                ])
                self.Elements = np.append(self.Elements, [
                    RTU.MakeTriangle(
                        [parsesize[0], parsesize[1], parsesize[2]],
                        [parsesize[3], parsesize[4], parsesize[5]],
                        [parsesize[6], parsesize[7], parsesize[8]]
                    )
                ], axis=0)
                self.ElementProperties.append(cp.deepcopy(proparray))
            except:
                print("Error parsing triangle coordinates.")
                return

        elif shape_text == "Box":
            # parse box
            szed = [float(x) for x in self.SizeEdit.text().split(",")]
            # 'direc' is the direction
            if len(szed) != 3 or len(direc) != 3:
                print("Please enter valid box size or direction.")
                return
            bx = RTU.MakeBox(pos, szed, direc)

            # For the box, we are effectively adding 6 planes
            if len(self.ElementProperties) > 1:
                # remove the last element if it was appended
                self.Elements = np.append(self.Elements[:-1], bx, axis=0)
                self.ElementProperties = self.ElementProperties[:-1]
            else:
                self.Elements = bx
                self.ElementProperties = []

            # each face is geometry=1(plane)
            for _ in range(6):
                self.ElementGeometry.append(1)
                self.ElementProperties.append(cp.deepcopy(proparray))

            # figure out a perpendicular vector if needed
            drv = direc
            if np.count_nonzero(drv) >= 2:
                if drv[2] == 0:
                    drv1 = [drv[1], -drv[0], drv[2]]
                elif drv[1] == 0:
                    drv1 = [-drv[2], drv[1], drv[0]]
                elif drv[0] == 0:
                    drv1 = [drv[0], -drv[2], drv[1]]
                    drv1 = np.cross(drv, drv1)
            else:
                drv1 = [drv[1], -drv[0], drv[2]]

            drv2 = drv
            # extend self.ElementSize with 6 planes
            self.ElementSize.extend([
                [szed[2], szed[1], drv1],
                [szed[2], szed[1], drv1],
                [szed[0], szed[2], drv2],
                [szed[0], szed[2], drv2],
                [szed[0], szed[1], drv],
                [szed[0], szed[1], drv]
            ])

        # Now handle element combo labeling
        if shape_text != "Box":
            self.n_els += 1
            self.ElementCombo.addItem(str(self.n_els) + " " + t + " " + shape_text)
        else:
            for _ in range(6):
                self.n_els += 1
                self.ElementCombo.addItem(str(self.n_els) + " " + t + " Plane")

    def UpdateTextBoxes(self):
        self.ElSizeList.setText(str(self.ElementSize))
        self.ElList.setText(str(self.Elements))
        self.ElPropList.setText(str(self.ElementProperties))
        self.ElGeoList.setText(str(self.ElementGeometry))
        self.MatList.setText(str(self.Materials))

    def FillSelectedElement(self):
        if len(self.ElementGeometry) > 0:
            ind = min(self.ElementCombo.currentIndex(), len(self.ElementGeometry) - 1)

            # Type
            t_idx = self.ElementProperties[ind][0] - 1
            self.TypeEdit.setCurrentIndex(t_idx)

            # Geometry
            g_idx = self.ElementGeometry[ind] - 1
            self.GeoEdit.setCurrentIndex(g_idx)

            # If interface, lens, dump:
            shape_text = self.TypeEdit.currentText()
            if shape_text == "Interface":
                mat1 = self.ElementProperties[ind][1]
                mat2 = self.ElementProperties[ind][2]
                idx_m1 = self.MaterialEdit1.findText(mat1)
                idx_m2 = self.MaterialEdit2.findText(mat2)
                if idx_m1 >= 0:
                    self.MaterialEdit1.setCurrentIndex(idx_m1)
                if idx_m2 >= 0:
                    self.MaterialEdit2.setCurrentIndex(idx_m2)
            elif shape_text == "Thin Lens":
                try:
                    self.FocalLengthEdit.setText(str(self.ElementProperties[ind][1]))
                except:
                    print("Could not parse focal length.")
            elif shape_text == "Dump":
                try:
                    self.DetSpecEdit.setText(str(self.ElementProperties[ind][3]))
                except:
                    print("Loading old file, Dump specification missing.")

            # Now geometry details
            geo_text = self.GeoEdit.currentText()
            elem = self.Elements[ind]
            if geo_text == "Plane":
                # plane
                self.PositionEdit.setText(str(elem[0])[1:-1])
                self.DirectionEdit.setText(str(elem[1])[1:-1])
                if len(self.ElementSize[ind]) > 2:
                    # e.g. [a0,a1,[x,y,z]]
                    s = self.ElementSize[ind]
                    self.SizeEdit.setText(
                        f"{s[0]},{s[1]},{s[2]}"
                    )
                else:
                    s = self.ElementSize[ind]
                    self.SizeEdit.setText(f"{s[0]},{s[1]}")
            elif geo_text in ["Circle", "Cylinder"]:
                self.PositionEdit.setText(str(elem[0])[1:-1])
                self.DirectionEdit.setText(str(elem[1])[1:-1])
                self.SizeEdit.setText(str(self.ElementSize[ind][0]))
            elif geo_text == "Triangle":
                # 3 points
                ps = self.ElementSize[ind]
                self.PointsEdit.setText(
                    str(ps[0])[1:-1] + "," +
                    str(ps[1])[1:-1] + "," +
                    str(ps[2])[1:-1]
                )

            self.DisableEdit()

    def FillSelectedMaterial(self):
        if len(self.Materials) > 0:
            ind = self.MaterialsCombo.currentIndex()
            mat = self.Materials[ind]
            # mat is [name, n, isAbs, conc_d, conc_a, color, isAligned, ...]
            self.NameEdit.setText(str(mat[0]))
            self.RefEdit.setText(str(mat[1]))
            self.AbsorbEdit.setCurrentIndex(int(mat[2]))
            self.ConcDEdit.setText(str(mat[3]))
            self.ConcAEdit.setText(str(mat[4]))
            self.colorEdit.setText(str(mat[5]))
            self.AlignedEdit.setCurrentIndex(int(mat[6]))
            self.AlignVEdit.setText(
                f"{mat[7][0]},{mat[7][1]},{mat[7][2]}"
            )
            self.AlignSEdit.setText(str(mat[8]))
            self.r_0Edit.setText(str(mat[9]))
            self.Q_DEdit.setText(str(mat[10]))
            self.Q_AEdit.setText(str(mat[11]))
            self.tau_DEdit.setText(str(mat[12]))
            self.tau_AEdit.setText(str(mat[13]))
            self.emax_DEdit.setText(str(mat[14]))
            self.emax_AEdit.setText(str(mat[15]))
            self.sp_abs_DEdit.setText(str(mat[16]))
            self.sp_em_DEdit.setText(str(mat[17]))
            self.sp_abs_AEdit.setText(str(mat[18]))
            self.sp_em_AEdit.setText(str(mat[19]))
            # Abs/Em angle
            try:
                self.AbsEmAngleEdit.setValue(mat[20])
            except:
                self.AbsEmAngleEdit.setValue(0)
            # FixDist
            try:
                self.FixDistEdit.setCurrentIndex(int(mat[21]))
                self.FixDistREdit.setText(str(mat[22]))
                self.FixDistNeighborsEdit.setValue(mat[23])
            except:
                self.FixDistEdit.setCurrentIndex(0)
                self.FixDistREdit.setText("0")
                self.FixDistNeighborsEdit.setValue(1)

            self.MatEnabler(False)

    def AcceptMaterial(self):
        self.AddMatBtn.setEnabled(False)
        if not self.EditMatMode:
            self.AddMaterial()
        else:
            self.EditMaterial()
        self.MatEnabler(False)
        self.UpdateTextBoxes()

    def DeleteMaterial(self):
        ind = self.MaterialsCombo.currentIndex()
        if ind < 0 or ind >= len(self.Materials):
            return

        self.MaterialsCombo.removeItem(ind)
        self.MaterialEdit1.removeItem(ind)
        self.MaterialEdit2.removeItem(ind)
        self.StMatPhEdit.removeItem(ind)
        del self.Materials[ind]
        self.UpdateTextBoxes()

    def EditMaterial(self):
        ind = self.MaterialsCombo.currentIndex()
        if ind < 0 or ind >= len(self.Materials):
            return
        try:
            v_al = [float(x) for x in self.AlignVEdit.text().split(",")]
        except:
            print("Invalid alignment vector.")
            return

        mat = [
            self.NameEdit.text(),
            float(self.RefEdit.text()),
            bool(self.AbsorbEdit.currentIndex()),
            float(self.ConcDEdit.text()),
            float(self.ConcAEdit.text()),
            self.colorEdit.text(),
            bool(self.AlignedEdit.currentIndex()),
            v_al,
            float(self.AlignSEdit.text()),
            float(self.r_0Edit.text()),
            float(self.Q_DEdit.text()),
            float(self.Q_AEdit.text()),
            float(self.tau_DEdit.text()),
            float(self.tau_AEdit.text()),
            float(self.emax_DEdit.text()),
            float(self.emax_AEdit.text()),
            self.sp_abs_DEdit.text(),
            self.sp_em_DEdit.text(),
            self.sp_abs_AEdit.text(),
            self.sp_em_AEdit.text(),
            self.AbsEmAngleEdit.value(),
            bool(self.FixDistEdit.currentIndex()),
            float(self.FixDistREdit.text()),
            self.FixDistNeighborsEdit.value()
        ]
        self.Materials[ind] = mat

        # Update combos
        name = mat[0]
        self.MaterialEdit1.setItemText(ind, name)
        self.MaterialEdit2.setItemText(ind, name)
        self.MaterialsCombo.setItemText(ind, name)
        self.StMatPhEdit.setItemText(ind, name)

    def AddMaterial(self):
        try:
            v_al = [float(x) for x in self.AlignVEdit.text().split(",")]
        except:
            print("Invalid alignment vector.")
            return
        mat = [
            self.NameEdit.text(),
            float(self.RefEdit.text()),
            bool(self.AbsorbEdit.currentIndex()),
            float(self.ConcDEdit.text()),
            float(self.ConcAEdit.text()),
            self.colorEdit.text(),
            bool(self.AlignedEdit.currentIndex()),
            v_al,
            float(self.AlignSEdit.text()),
            float(self.r_0Edit.text()),
            float(self.Q_DEdit.text()),
            float(self.Q_AEdit.text()),
            float(self.tau_DEdit.text()),
            float(self.tau_AEdit.text()),
            float(self.emax_DEdit.text()),
            float(self.emax_AEdit.text()),
            self.sp_abs_DEdit.text(),
            self.sp_em_DEdit.text(),
            self.sp_abs_AEdit.text(),
            self.sp_em_AEdit.text(),
            self.AbsEmAngleEdit.value(),
            bool(self.FixDistEdit.currentIndex()),
            float(self.FixDistREdit.text()),
            self.FixDistNeighborsEdit.value()
        ]

        self.Materials.append(mat)
        self.MatList.setText(str(self.Materials))

        # Add to combos
        name = mat[0]
        self.MaterialEdit1.addItem(name)
        self.MaterialEdit2.addItem(name)
        self.MaterialsCombo.addItem(name)
        self.StMatPhEdit.addItem(name)

        self.ClearMaterials('0')

    def EditMatFun(self):
        self.FillSelectedMaterial()
        self.EditMatEnabler()

    def EditMatEnabler(self):
        self.MatEnabler(True)
        self.AddMatBtn.setEnabled(True)
        self.AbsorberDisabler()
        self.AlignedDisabler()
        self.FixDistDisabler()
        self.EditMatMode = True

    def NewMatEnabler(self):
        self.AddMatBtn.setEnabled(True)
        self.MatEnabler(True)
        self.AbsorberDisabler()
        self.AlignedDisabler()
        self.FixDistDisabler()
        self.EditMatMode = False

    def MatEnabler(self, state):
        self.NameEdit.setEnabled(state)
        self.RefEdit.setEnabled(state)
        self.AbsorbEdit.setEnabled(state)
        self.ConcAEdit.setEnabled(state)
        self.ConcDEdit.setEnabled(state)
        self.colorEdit.setEnabled(state)
        self.AlignedEdit.setEnabled(state)
        self.AlignVEdit.setEnabled(state)
        self.AlignSEdit.setEnabled(state)
        self.FixDistEdit.setEnabled(state)
        self.FixDistREdit.setEnabled(state)
        self.FixDistNeighborsEdit.setEnabled(state)
        self.r_0Edit.setEnabled(state)
        self.Q_DEdit.setEnabled(state)
        self.Q_AEdit.setEnabled(state)
        self.tau_DEdit.setEnabled(state)
        self.tau_AEdit.setEnabled(state)
        self.emax_DEdit.setEnabled(state)
        self.emax_AEdit.setEnabled(state)
        self.sp_abs_DEdit.setEnabled(state)
        self.sp_em_DEdit.setEnabled(state)
        self.sp_abs_AEdit.setEnabled(state)
        self.sp_em_AEdit.setEnabled(state)
        self.AbsEmAngleEdit.setEnabled(state)

    def AbsorberDisabler(self):
        state = bool(self.AbsorbEdit.currentIndex())
        self.ConcAEdit.setEnabled(state)
        self.ConcDEdit.setEnabled(state)
        self.AlignedEdit.setEnabled(state)
        self.AlignVEdit.setEnabled(state)
        self.AlignSEdit.setEnabled(state)
        self.FixDistEdit.setEnabled(state)
        self.FixDistREdit.setEnabled(state)
        self.FixDistNeighborsEdit.setEnabled(state)
        self.r_0Edit.setEnabled(state)
        self.Q_DEdit.setEnabled(state)
        self.Q_AEdit.setEnabled(state)
        self.tau_DEdit.setEnabled(state)
        self.tau_AEdit.setEnabled(state)
        self.emax_DEdit.setEnabled(state)
        self.emax_AEdit.setEnabled(state)
        self.sp_abs_DEdit.setEnabled(state)
        self.sp_em_DEdit.setEnabled(state)
        self.sp_abs_AEdit.setEnabled(state)
        self.sp_em_AEdit.setEnabled(state)
        self.AlignedDisabler()
        self.FixDistDisabler()

    def AlignedDisabler(self):
        state = bool(self.AlignedEdit.currentIndex())
        self.AlignVEdit.setEnabled(state)
        self.AlignSEdit.setEnabled(state)
        self.AbsEmAngleEdit.setEnabled(state)

    def FixDistDisabler(self):
        state = bool(self.FixDistEdit.currentIndex())
        self.FixDistREdit.setEnabled(state)
        self.FixDistNeighborsEdit.setEnabled(state)

    def ClearMaterials(self, text=''):
        self.NameEdit.setText(text)
        self.RefEdit.setText(text)
        self.ConcAEdit.setText(text)
        self.ConcDEdit.setText(text)
        self.colorEdit.setText(text)
        self.AlignVEdit.setText(f"{text},{text},{text}")
        self.AlignSEdit.setText(text)
        self.FixDistREdit.setText(text)
        self.FixDistNeighborsEdit.setValue(0)
        self.r_0Edit.setText(text)
        self.Q_DEdit.setText(text)
        self.Q_AEdit.setText(text)
        self.tau_DEdit.setText(text)
        self.tau_AEdit.setText(text)
        self.emax_DEdit.setText(text)
        self.emax_AEdit.setText(text)
        self.sp_abs_DEdit.setText(text)
        self.sp_em_DEdit.setText(text)
        self.sp_abs_AEdit.setText(text)
        self.sp_em_AEdit.setText(text)
        self.AbsEmAngleEdit.setValue(0)

    def TypeEditDisabler(self):
        t = self.TypeEdit.currentText()
        if t == "Mirror":
            self.MaterialEdit1.setEnabled(False)
            self.MaterialEdit2.setEnabled(False)
            self.FocalLengthEdit.setEnabled(False)
            self.DetSpecEdit.setEnabled(False)
        elif t == "Interface":
            self.MaterialEdit1.setEnabled(True)
            self.MaterialEdit2.setEnabled(True)
            self.FocalLengthEdit.setEnabled(False)
            self.DetSpecEdit.setEnabled(False)
        elif t == "Thin Lens":
            self.MaterialEdit1.setEnabled(False)
            self.MaterialEdit2.setEnabled(False)
            self.FocalLengthEdit.setEnabled(True)
            self.DetSpecEdit.setEnabled(False)
        elif t == "Dump":
            self.MaterialEdit1.setEnabled(False)
            self.MaterialEdit2.setEnabled(False)
            self.FocalLengthEdit.setEnabled(False)
            self.DetSpecEdit.setEnabled(True)

    def GeoEditDisabler(self):
        geo = self.GeoEdit.currentText()
        if geo == "Plane":
            self.PositionEdit.setEnabled(True)
            self.DirectionEdit.setEnabled(True)
            self.SizeEdit.setEnabled(True)
            self.PointsEdit.setEnabled(False)
        elif geo == "Circle":
            self.PositionEdit.setEnabled(True)
            self.DirectionEdit.setEnabled(True)
            self.SizeEdit.setEnabled(True)
            self.PointsEdit.setEnabled(False)
        elif geo == "Cylinder":
            self.PositionEdit.setEnabled(True)
            self.DirectionEdit.setEnabled(True)
            self.SizeEdit.setEnabled(True)
            self.PointsEdit.setEnabled(False)
        elif geo == "Triangle":
            self.PositionEdit.setEnabled(False)
            self.DirectionEdit.setEnabled(False)
            self.SizeEdit.setEnabled(False)
            self.PointsEdit.setEnabled(True)
        elif geo == "Box":
            self.PositionEdit.setEnabled(True)
            self.DirectionEdit.setEnabled(True)
            self.SizeEdit.setEnabled(True)
            self.PointsEdit.setEnabled(False)

    def DrawCylinder(self, p1, p2, r, i):
        """
        Draw a Cylinder from p1 to p2 with a radius r
        Adapted from:
        https://stackoverflow.com/questions/32317247/
               how-to-draw-a-cylinder-using-matplotlib-along-length-of-point-x1-y1-and-x2-y2
        """
        clr = i
        a = p2 - p1
        abs_a = np.linalg.norm(a)
        if abs_a == 0:
            return
        a /= abs_a

        not_a = np.array([1, 0, 0])
        if np.allclose(a, not_a):
            not_a = np.array([0, 1, 0])

        n1 = np.cross(a, not_a)
        n1 /= np.linalg.norm(n1)
        n2 = np.cross(a, n1)

        t = np.linspace(0, abs_a, 6)
        theta = np.linspace(0, 2 * np.pi, 6)
        t, theta = np.meshgrid(t, theta)

        X, Y, Z = [
            p1[j] + a[j] * t + r * np.sin(theta) * n1[j] + r * np.cos(theta) * n2[j]
            for j in [0, 1, 2]
        ]
        self.ax.plot_surface(
            X, Y, Z,
            color=self.colors[self.ElementProperties[clr][0] - 1],
            alpha=self.alphas[self.ElementProperties[clr][0] - 1]
        )

    def DrawCircle(self, M, n, r, i):
        """
        Draw Circle with radius r at center M with normal n.
        """
        not_n = np.array([1, 0, 0])
        if np.allclose(n, not_n):
            not_n = np.array([0, 1, 0])

        n1 = np.cross(n, not_n)
        n1 /= np.linalg.norm(n1)
        n2 = np.cross(n, n1)

        theta = np.linspace(0, 2 * np.pi, 100)
        X, Y, Z = [
            M[j] + r * np.sin(theta) * n1[j] + r * np.cos(theta) * n2[j]
            for j in [0, 1, 2]
        ]
        verts = [list(zip(X, Y, Z))]
        Collection = Poly3DCollection(
            verts, linewidths=1,
            alpha=self.alphas[self.ElementProperties[i][0] - 1]
        )
        FaceColor = self.colors[self.ElementProperties[i][0] - 1]
        Collection.set_facecolor(FaceColor)
        Collection.set_edgecolor("black")
        self.ax.add_collection3d(Collection)

    def DrawFigure(self):
        """
        Draw the 3D elements in a separate matplotlib figure.
        """
        plt.ion()
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111, projection='3d')
        self.ax.set_xlabel('x')
        self.ax.set_ylabel('y')
        self.ax.set_zlabel('z')

        i = 0
        for el in self.Elements:
            geom = self.ElementGeometry[i]
            # Plane
            if geom == 1:
                # e1,e2 for bounding box
                if len(self.ElementSize[i]) == 2:
                    # normal plane
                    e1 = np.array([-el[1][1], el[1][0], 0], dtype=float)
                    e2 = np.array([0, -el[1][2], el[1][1]], dtype=float)
                    ealt = np.array([-el[1][2], 0, el[1][0]], dtype=float)
                    if np.allclose(e1, [0, 0, 0]):
                        e1 = ealt
                    elif np.allclose(e2, [0, 0, 0]):
                        e2 = ealt
                    elif np.allclose(e1, e2) or np.allclose(e1, -e2):
                        e2 = ealt
                else:
                    # if rotated around axis
                    e1 = np.array(self.ElementSize[i][2], dtype=float)
                    e2 = np.cross(e1, el[1])

                e1 /= np.linalg.norm(e1)
                e2 /= np.linalg.norm(e2)

                h = self.ElementSize[i][0]
                w = self.ElementSize[i][1]
                a = el[0] + e1 * h + e2 * w
                b = el[0] + e1 * h - e2 * w
                c = el[0] - e1 * h - e2 * w
                d = el[0] - e1 * h + e2 * w
                x = [a[0], b[0], c[0], d[0]]
                y = [a[1], b[1], c[1], d[1]]
                z = [a[2], b[2], c[2], d[2]]

                verts = [list(zip(x, y, z))]
                Collection = Poly3DCollection(
                    verts, linewidths=1,
                    alpha=self.alphas[self.ElementProperties[i][0] - 1]
                )
                FaceColor = self.colors[self.ElementProperties[i][0] - 1]
                Collection.set_facecolor(FaceColor)
                Collection.set_edgecolor("black")
                self.ax.add_collection3d(Collection)

            elif geom == 2:
                # Circle
                self.DrawCircle(
                    np.array(el[0], dtype=float),
                    np.array(el[1], dtype=float),
                    self.ElementSize[i][0],
                    i
                )
            elif geom == 3:
                # Cylinder
                p1 = np.array(el[0], dtype=float)
                p2 = np.array(el[1], dtype=float)
                r_val = self.ElementSize[i][0]
                self.DrawCylinder(p1, p2, r_val, i)
            elif geom == 4:
                # Triangle
                p1 = self.ElementSize[i][0]
                p2 = self.ElementSize[i][1]
                p3 = self.ElementSize[i][2]
                x = [p1[0], p2[0], p3[0]]
                y = [p1[1], p2[1], p3[1]]
                z = [p1[2], p2[2], p3[2]]
                verts = [list(zip(x, y, z))]
                Collection = Poly3DCollection(
                    verts, linewidths=1,
                    alpha=self.alphas[self.ElementProperties[i][0] - 1]
                )
                FaceColor = self.colors[self.ElementProperties[i][0] - 1]
                Collection.set_facecolor(FaceColor)
                Collection.set_edgecolor("black")
                self.ax.add_collection3d(Collection)

            i += 1

        plt.show()


if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = Example()
    sys.exit(app.exec_())
