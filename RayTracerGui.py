#!/usr/bin/env python
import sys
import numpy as np
import cPickle
from PyQt5.QtWidgets import (QWidget, QLabel, QLineEdit,
    QTextEdit, QGridLayout, QApplication, QPushButton,QComboBox, QFileDialog, QSpinBox, QCheckBox)
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
#import matplotlib.pyplot as plt
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
        #Disable stuff for the start
        self.DisableEdit()

        self.ElList = QTextEdit()
        self.ElPropList = QTextEdit()
        self.ElSizeList = QTextEdit()
        self.ElGeoList = QTextEdit()


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
        grid.addWidget(ActionBtn,14, 3)

        '''def __init__(self, n, IsAbsorber = False, conc_a=0, conc_b=0, color = 'blue', IsAligned = False , AlignVect = [1,0,0], AlignS = 0.1, R_0= 5, Q_D = 0.9, Q_A = 0.9, tau_D = 8, tau_A = 8, emax_D = 23500, emax_A = 54000, sp_abs_D = "Coumarin1_Abs.txt", sp_em_D = "Coumarin1_Em.txt", sp_abs_A = "Coumarin6_Abs.txt", sp_em_A = "Coumarin6_Em.txt"):'''

        MaterialTitle = QLabel('Add Material')
        Name = QLabel('Material Name')
        Ref = QLabel('Refractive Index')
        Absorb =  QLabel('Is Absorber ?')
        ConcD = QLabel('Concentration Donor')
        ConcA = QLabel('Concentration Acceptor')
        color = QLabel('Color')
        Aligned = QLabel('Is Aligned ?')
        AlignV = QLabel('Direction of alignment (x,y,z)')
        AlignS = QLabel('Standard Deviation of alignment in rad')
        AbsEmAngle = QLabel('Angle between absorption and Emission Dipole moment')

        FixDist = QLabel('Fixed Distance between donor & Acceptor')
        FixDistR = QLabel('Distance')
        FixDistNeighbors = QLabel('Number of Donors per Acceptor')

        r_0 = QLabel('FRET Radius R0 in nm')
        Q_D = QLabel('Fl. Quantum yiel of donor')
        Q_A = QLabel('Fl. Quantum yiel of acceptor')
        tau_D = QLabel('Lifetime donor')
        tau_A = QLabel('Lifetime acceptor')
        emax_D = QLabel('Maximum extincion donor M^-1 cm^-1')
        emax_A = QLabel('Maximum extincion acceptor M^-1 cm^-1')
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
        grid.addWidget(self.MaterialsCombo,15,1)

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


        ShowBtn =  QPushButton("Show")
        ShowBtn.clicked.connect(self.DrawFigure)
        grid.addWidget(ShowBtn,44,0)

        RunBtn =  QPushButton("Run")
        RunBtn.clicked.connect(self.RunRayTracing)
        grid.addWidget(RunBtn,45,0)

        self.MatEnabler(False)
        self.ClearMaterials("0")


        '''Photon generation        '''
        PhLabel = QLabel("Photon Generation")
        NumPh = QLabel("Number of Photons")
        StartPh = QLabel("Starting Point Photon (x,y,z)")
        DirectionPh = QLabel("Direction of Emission")
        EmDPh = QLabel("Emission Direction Distribution (x,y,z)")
        SpectWavel = QLabel("Emission Wavelength")
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
        self.WavelPhEdit.setRange(200,2400)
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

        grid.addWidget(StMatPh, 23,4)
        grid.addWidget(self.StMatPhEdit,23,5)

        grid.addWidget(EmPol,24,4)
        grid.addWidget(self.EmPolEdit,24,5)

        grid.addWidget(ScaleFactor,25,4)
        grid.addWidget(self.ScaleFactorEdit,25,5)

        self.setLayout(grid)


        self.setGeometry(300, 300, 350, 300)
        self.setWindowTitle('RayTracerGui')
        self.show()


    def ExtScript(self):
        ExScr.main(self)
        # execfile('ExternalScript.py')

    def ClearAll(self):
        self.FirstElement = True
        self.Elements = np.zeros((1,2,3))
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
        LSC = RayTrace.RayTrace()        
        LSC.IOOperator.ParsePath(foldername)
        #LSC.ParsePath(foldername)
        LSC.runSimulation()

    def SaveElements(self):
        foldername = str(QFileDialog.getExistingDirectory(self, "Select Directory"))

        f = open(foldername + "/Elements.prt",'wb')
        cPickle.dump(self.Elements,f)
        f = open(foldername + "/ElementProperties.prt",'wb')
        cPickle.dump(self.ElementProperties,f)
        f = open(foldername + "/ElementGeometry.prt",'wb')
        cPickle.dump(self.ElementGeometry,f)
        f = open(foldername + "/ElementSize.prt",'wb')
        cPickle.dump(self.ElementSize,f)
        f = open(foldername + "/Materials.prt",'wb')
        cPickle.dump(self.Materials,f)

        '''Create and save array for photon generation'''
        try:
            self.PhArray = [self.NumPhEdit.value(), map(float,self.StartPhEdit.text().split(",")), self.EmDPhEdit.currentText(),self.SpectWavel.currentText(),self.WavelPhEdit.value(),self.StMatPhEdit.currentText(), map(float,self.DirectionPhEdit.text().split(",")), self.EmPolEdit.isChecked(),float(self.ScaleFactorEdit.text())]
        except:
            print "Please enter all Photon Parameters"

        f = open(foldername + "/Photon.prt",'wb')
        cPickle.dump(self.PhArray,f)


        # self.ElementProperties.tofile(foldername + "/ElementProperties.prt")
        # self.ElementGeometry.tofile(foldername + "/ElementGeometry.prt")
        # self.ElementSize.tofile(foldername + "/ElementSize.prt")

    def LoadElements(self):
        self.ClearAll()

        foldername = str(QFileDialog.getExistingDirectory(self, "Select Directory"))
        types = {1:"Mirror",2:"Interface",3:"ThinLens",4:"Dump"}
        geos = {1:"Plane", 2:"Circle",3:"Cylinder",4:"Triangle"}

        f = open(foldername + "/Elements.prt",'rb')
        self.Elements = np.array(cPickle.load(f))
        f = open(foldername + "/ElementProperties.prt",'rb')
        self.ElementProperties = cPickle.load(f)
        f = open(foldername + "/ElementGeometry.prt",'rb')
        self.ElementGeometry = cPickle.load(f)
        f = open(foldername + "/ElementSize.prt",'rb')
        self.ElementSize = cPickle.load(f)
        f = open(foldername + "/Materials.prt",'rb')
        self.Materials = cPickle.load(f)
        f = open(foldername + "/Photon.prt",'rb')
        self.PhArray = cPickle.load(f)


        for i in range(len(self.Elements)):
            self.ElementCombo.addItem(str(i+1) + " " + types[self.ElementProperties[i][0]] + " " + geos[self.ElementGeometry[i]])
        for i in range(len(self.Materials)):
            self.MaterialEdit1.addItem(self.Materials[i][0])
            self.MaterialEdit2.addItem(self.Materials[i][0])
            self.MaterialsCombo.addItem(self.Materials[i][0])
            self.StMatPhEdit.addItem(self.Materials[i][0])

        self.n_els = len(self.Elements) + 1
        self.FirstElement = False
        self.FillSelectedElement()
        self.FillSelectedMaterial()
        self.FillPhoton()
        self.UpdateTextBoxes()

    def FillPhoton(self):
        '''[self.NumPhEdit.value(), map(float,self.StartPhEdit.text().split(",")), self.EmDPhEdit.currentText(),self.WavelPhEdit.value(),self.StMatPhEdit.currentText()]'''
        self.NumPhEdit.setValue(self.PhArray[0])
        self.StartPhEdit.setText(str(self.PhArray[1][0]) + "," + str(self.PhArray[1][1]) + "," + str(self.PhArray[1][2]))
        index = self.EmDPhEdit.findText(self.PhArray[2])
        self.EmDPhEdit.setCurrentIndex(index)
        index = self.SpectWavel.findText(self.PhArray[3])
        self.SpectWavel.setCurrentIndex(index)
        self.WavelPhEdit.setValue(self.PhArray[4])
        index = self.StMatPhEdit.findText(self.PhArray[5])
        self.StMatPhEdit.setCurrentIndex(index)
        self.DirectionPhEdit.setText(str(self.PhArray[6][0]) + "," + str(self.PhArray[6][1]) + "," + str(self.PhArray[6][2]))
        self.EmPolEdit.setChecked(self.PhArray[7])

        if len(self.PhArray) == 9:
            self.ScaleFactorEdit.setText(str(self.PhArray[8]))
        else:
            self.ScaleFactorEdit.setText("1.0")

    def AcceptElement(self):
        self.DisableEdit()
        self.AddBtn.setEnabled(False)

        if self.EditMode == False:
            self.AddElement()
        else:
            self.EditElement()

        self.UpdateTextBoxes()


    def DeleteElement(self):
        if len(self.ElementGeometry) > 0:
            ind = self.ElementCombo.currentIndex()
            self.Elements =  np.delete(self.Elements,ind,axis=0)
            self.n_els -= 1
            del self.ElementProperties[ind]
            del self.ElementGeometry[ind]
            del self.ElementSize[ind]
            self.ElementCombo.removeItem(ind)
            self.UpdateTextBoxes()

    def EditElement(self):
        ind = self.ElementCombo.currentIndex()

        if self.TypeEdit.currentText() == "Mirror":
            proparray = [1]
        elif self.TypeEdit.currentText() == "Interface":
            proparray = [2,self.MaterialEdit1.currentText(),self.MaterialEdit2.currentText()]
        elif self.TypeEdit.currentText() == "Thin Lens":
            proparray = [3,float(self.FocalLengthEdit.text())]
        elif self.TypeEdit.currentText() == "Dump":
            proparray = [4,0,0,self.DetSpecEdit.text()]

        self.Elements[ind] = [map(float,self.PositionEdit.text().split(",")),map(float,self.DirectionEdit.text().split(","))]
        self.ElementProperties[ind] = cp.deepcopy(proparray)
        self.ElementGeometry[ind] = self.GeoEdit.currentIndex() + 1

        if self.GeoEdit.currentText() == "Plane":
            parsesize = self.SizeEdit.text().split(",")
            if len(parsesize) == 2:
                self.ElementSize[ind] = [float(parsesize[0]),float(parsesize[1])]
            else:
                self.ElementSize[ind] = [float(parsesize[0]),float(parsesize[1]),[float(parsesize[2][1:]),float(parsesize[3]),float(parsesize[4][:-1])]]
        elif self.GeoEdit.currentText() == "Circle" or self.GeoEdit.currentText() == "Cylinder":
            self.ElementSize[ind] = [float(self.SizeEdit.text())]
        elif self.GeoEdit.currentText() == "Triangle":
            parsesize = map(float,self.PointsEdit.text().split(","))
            self.ElementSize[ind] = [[parsesize[0],parsesize[1],parsesize[2]],[parsesize[3],parsesize[4],parsesize[5]],[parsesize[6],parsesize[7],parsesize[8]]]
            self.Elements[ind] = RTU.MakeTriangle([parsesize[0],parsesize[1],parsesize[2]],[parsesize[3],parsesize[4],parsesize[5]],[parsesize[6],parsesize[7],parsesize[8]])
         

        self.ElementCombo.setItemText(ind, str(ind+1) + " " + self.TypeEdit.currentText() + " " + self.GeoEdit.currentText())


    def AddElement(self):


        if self.TypeEdit.currentText() == "Mirror":
            proparray = [1]
        elif self.TypeEdit.currentText() == "Interface":
            proparray = [2,self.MaterialEdit1.currentText(),self.MaterialEdit2.currentText()]
        elif self.TypeEdit.currentText() == "Thin Lens":
            proparray = [3,float(self.FocalLengthEdit.text())]
        elif self.TypeEdit.currentText() == "Dump":
            proparray = [4,0,0,self.DetSpecEdit.text()]

        if self.FirstElement == True:
            self.Elements = np.array([[map(float,self.PositionEdit.text().split(",")),map(float,self.DirectionEdit.text().split(","))]])
            self.ElementProperties = cp.deepcopy([proparray])
            self.FirstElement = False
        elif self.GeoEdit.currentText() != "Triangle":
            self.Elements = np.append(self.Elements,[[map(float,self.PositionEdit.text().split(",")),map(float,self.DirectionEdit.text().split(","))]],axis=0)
            self.ElementProperties.append(proparray)


        if self.GeoEdit.currentIndex() < 4:
            self.ElementGeometry.append(self.GeoEdit.currentIndex() + 1)
        elif self.GeoEdit.currentIndex() == 4:
            szed = map(float,self.SizeEdit.text().split(","))
            drv = map(float,self.DirectionEdit.text().split(","))

            if len(szed) != 3 or len(drv) != 3:
                print "Please enter valid coordinates"
            else:
                bx = RTU.MakeBox(map(float,self.PositionEdit.text().split(",")), szed , drv)

                if len(self.ElementProperties)>1:
                    self.Elements = np.append(self.Elements[:-1], bx,axis=0)
                    self.ElementProperties = self.ElementProperties[:-1]
                else:
                    self.Elements = bx
                    self.ElementProperties = []

                for i in range(6):
                    self.ElementGeometry.append(1)
                    self.ElementProperties.append(proparray)

                if np.count_nonzero(drv)>=2:
                    if drv[2] == 0:
                        drv1 = [drv[1],-drv[0],drv[2]]
                    elif drv[1] == 0:
                        drv1 = [-drv[2],drv[1],drv[0]]
                    elif drv[0] == 0:
                        drv1 = [drv[0],-drv[2],drv[1]]
                        drv1 = np.cross(drv,drv1)
                else:
                    drv1 = [drv[1],-drv[0],drv[2]]

                drv2 = drv

                self.ElementSize.extend([[szed[2],szed[1],drv1],[szed[2],szed[1],drv1],[szed[0],szed[2],drv2],[szed[0],szed[2],drv2],[szed[0],szed[1],drv],[szed[0],szed[1],drv]])


        if self.GeoEdit.currentText() == "Plane":
            parsesize = self.SizeEdit.text().split(",")
            
            if len(parsesize) == 2:
                self.ElementSize.append([float(parsesize[0]),float(parsesize[1])])
            else:
                self.ElementSize.append([float(parsesize[0]),float(parsesize[1]),[float(parsesize[2][1:]),float(parsesize[3]),float(parsesize[4][:-1])]])
        elif self.GeoEdit.currentText() == "Circle" or self.GeoEdit.currentText() == "Cylinder":
            self.ElementSize.append([float(self.SizeEdit.text())])
        elif self.GeoEdit.currentText() == "Triangle":
            parsesize = map(float,self.PointsEdit.text().split(","))
            self.ElementSize.append([[parsesize[0],parsesize[1],parsesize[2]],[parsesize[3],parsesize[4],parsesize[5]],[parsesize[6],parsesize[7],parsesize[8]]])
            self.Elements = np.append(self.Elements,[RTU.MakeTriangle([parsesize[0],parsesize[1],parsesize[2]],[parsesize[3],parsesize[4],parsesize[5]],[parsesize[6],parsesize[7],parsesize[8]])],axis=0)
            self.ElementProperties.append(proparray)


        if self.GeoEdit.currentText() !="Box":
            self.n_els += 1
            self.ElementCombo.addItem(str(self.n_els) + " " + self.TypeEdit.currentText() + " " + self.GeoEdit.currentText())
        else:
            for i in range(6):
                self.n_els += 1
                self.ElementCombo.addItem(str(self.n_els) + " " + self.TypeEdit.currentText() + " " + "Plane")




    def UpdateTextBoxes(self):
        self.ElSizeList.setText(str(self.ElementSize))
        self.ElList.setText(str(self.Elements))
        self.ElPropList.setText(str(self.ElementProperties))
        self.ElGeoList.setText(str(self.ElementGeometry))
        self.MatList.setText(str(self.Materials))

    def FillSelectedElement(self):
        if len(self.ElementGeometry) > 0:

            
            ind = np.minimum(self.ElementCombo.currentIndex(),len(self.ElementGeometry)-1)



            self.TypeEdit.setCurrentIndex(self.ElementProperties[ind][0]-1)
            self.GeoEdit.setCurrentIndex(self.ElementGeometry[ind] - 1)

            if self.TypeEdit.currentText() == "Interface":
                index = self.MaterialEdit1.findText(self.ElementProperties[ind][1])
                self.MaterialEdit1.setCurrentIndex(index)
                index = self.MaterialEdit2.findText(self.ElementProperties[ind][2])
                self.MaterialEdit2.setCurrentIndex(index)
            elif self.TypeEdit.currentText() == "Thin Lens":
                self.FocalLengthEdit.setText(str(self.ElementProperties[ind][1]))
            elif self.TypeEdit.currentText() == "Dump":
                try:
                    self.DetSpecEdit.setText(str(self.ElementProperties[ind][3]))
                except:
                    print "Loading Old File"


            if self.GeoEdit.currentText() == "Plane":
                self.PositionEdit.setText(str(self.Elements[ind][0].tolist())[1:-1])
                self.DirectionEdit.setText(str(self.Elements[ind][1].tolist())[1:-1])

                if len(self.ElementSize[ind]) > 2:
                    self.SizeEdit.setText(str(self.ElementSize[ind][0]) + "," + str(self.ElementSize[ind][1]) + "," + str(self.ElementSize[ind][2]))
                else:
                    self.SizeEdit.setText(str(self.ElementSize[ind][0]) + "," + str(self.ElementSize[ind][1]))

            elif self.GeoEdit.currentText() == "Circle" or self.GeoEdit.currentText() == "Cylinder":
                self.PositionEdit.setText(str(self.Elements[ind][0].tolist())[1:-1])
                self.DirectionEdit.setText(str(self.Elements[ind][1].tolist())[1:-1])
                self.SizeEdit.setText(str(self.ElementSize[ind][0]))

            elif self.GeoEdit.currentText() == "Triangle":
                self.PointsEdit.setText(str(self.ElementSize[ind][0])[1:-1] + "," + str(self.ElementSize[ind][1])[1:-1] + "," + str(self.ElementSize[ind][2])[1:-1])

            self.DisableEdit()



    def FillSelectedMaterial(self):
        if len(self.Materials) > 0:
            ind = self.MaterialsCombo.currentIndex()
            self.NameEdit.setText(str(self.Materials[ind][0]))
            self.RefEdit.setText(str(self.Materials[ind][1]))
            self.AbsorbEdit.setCurrentIndex(int(self.Materials[ind][2]))
            self.ConcDEdit.setText(str(self.Materials[ind][3]))
            self.ConcAEdit.setText(str(self.Materials[ind][4]))
            self.colorEdit.setText(str(self.Materials[ind][5]))
            self.AlignedEdit.setCurrentIndex(int(self.Materials[ind][6]))
            self.AlignVEdit.setText(str(self.Materials[ind][7][0]) + "," + str(self.Materials[ind][7][1]) + "," + str(self.Materials[ind][7][2]) )
            self.AlignSEdit.setText(str(self.Materials[ind][8]))
            self.r_0Edit.setText(str(self.Materials[ind][9]))
            self.Q_DEdit.setText(str(self.Materials[ind][10]))
            self.Q_AEdit.setText(str(self.Materials[ind][11]))
            self.tau_DEdit.setText(str(self.Materials[ind][12]))
            self.tau_AEdit.setText(str(self.Materials[ind][13]))
            self.emax_DEdit.setText(str(self.Materials[ind][14]))
            self.emax_AEdit.setText(str(self.Materials[ind][15]))
            self.sp_abs_DEdit.setText(str(self.Materials[ind][16]))
            self.sp_em_DEdit.setText(str(self.Materials[ind][17]))
            self.sp_abs_AEdit.setText(str(self.Materials[ind][18]))
            self.sp_em_AEdit.setText(str(self.Materials[ind][19]))
            '''Only temporarily necessary for old files'''
            try:
                self.AbsEmAngleEdit.setValue(self.Materials[ind][20])
            except:
                self.AbsEmAngleEdit.setValue(0)
            try:
                self.FixDistEdit.setCurrentIndex(int(self.Materials[ind][21]))
                self.FixDistREdit.setText(str(self.Materials[ind][22]))
                self.FixDistNeighborsEdit.setValue(self.Materials[ind][23])
            except:
                self.FixDistEdit.setCurrentIndex(int(0))
                self.FixDistREdit.setText(str(0))
                self.FixDistNeighborsEdit.setValue(1)

            self.MatEnabler(False)

    def AcceptMaterial(self):
        self.AddMatBtn.setEnabled(False)
        if self.EditMatMode == False:
            self.AddMaterial()
        else:
            self.EditMaterial()
        self.MatEnabler(False)
        self.UpdateTextBoxes()

    def DeleteMaterial(self):
        ind = self.MaterialsCombo.currentIndex()
        self.MaterialsCombo.removeItem(ind)
        self.MaterialEdit1.removeItem(ind)
        self.MaterialEdit2.removeItem(ind)
        self.StMatPhEdit.removeItem(ind)
        del self.Materials[ind]
        self.UpdateTextBoxes()

    def EditMaterial(self):
        ind = self.MaterialsCombo.currentIndex()
        v_al = map(float,self.AlignVEdit.text().split(","))

        indexm1 = self.MaterialEdit1.findText(self.Materials[ind][0])
        indexm2 = self.MaterialEdit2.findText(self.Materials[ind][0])
        indexPh = self.StMatPhEdit.findText(self.Materials[ind][0])

        self.Materials[ind] = [self.NameEdit.text(),float(self.RefEdit.text()),bool(self.AbsorbEdit.currentIndex()),float(self.ConcDEdit.text()),float(self.ConcAEdit.text()),self.colorEdit.text(),bool(self.AlignedEdit.currentIndex()),v_al,float(self.AlignSEdit.text()),float(self.r_0Edit.text()),float(self.Q_DEdit.text()),float(self.Q_AEdit.text()),float(self.tau_DEdit.text()),float(self.tau_AEdit.text()),float(self.emax_DEdit.text()),float(self.emax_AEdit.text()),self.sp_abs_DEdit.text(),self.sp_em_DEdit.text(),self.sp_abs_AEdit.text(),self.sp_em_AEdit.text(), self.AbsEmAngleEdit.value(), bool(self.FixDistEdit.currentIndex()),float(self.FixDistREdit.text()),self.FixDistNeighborsEdit.value()]

        self.MaterialEdit1.setItemText(indexm1,self.Materials[ind][0])
        self.MaterialEdit2.setItemText(indexm2,self.Materials[ind][0])
        self.MaterialsCombo.setItemText(ind,self.Materials[ind][0])
        self.StMatPhEdit.setItemText(indexPh,self.Materials[ind][0])


    def AddMaterial(self):
        # self.Materials = np.append(self.Materials,[self.NameEdit.text(),bool(self.AbsorbEdit.currentIndex()),float(self.RefEdit.text())])
        v_al = map(float,self.AlignVEdit.text().split(","))

        self.Materials.append([self.NameEdit.text(),float(self.RefEdit.text()),bool(self.AbsorbEdit.currentIndex()),float(self.ConcDEdit.text()),float(self.ConcAEdit.text()),self.colorEdit.text(),bool(self.AlignedEdit.currentIndex()),v_al,float(self.AlignSEdit.text()),float(self.r_0Edit.text()),float(self.Q_DEdit.text()),float(self.Q_AEdit.text()),float(self.tau_DEdit.text()),float(self.tau_AEdit.text()),float(self.emax_DEdit.text()),float(self.emax_AEdit.text()),self.sp_abs_DEdit.text(),self.sp_em_DEdit.text(),self.sp_abs_AEdit.text(),self.sp_em_AEdit.text(),self.AbsEmAngleEdit.value(), bool(self.FixDistEdit.currentIndex()),float(self.FixDistREdit.text()),self.FixDistNeighborsEdit.value()])

        self.MatList.setText(str(self.Materials))
        self.MaterialEdit1.addItem(self.NameEdit.text())
        self.MaterialEdit2.addItem(self.NameEdit.text())
        self.MaterialsCombo.addItem(self.NameEdit.text())
        self.StMatPhEdit.addItem(self.NameEdit.text())

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

    def MatEnabler(self,state):

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
        # self.AbsorberDisabler()
        # self.AlignedDisabler()



    def AbsorberDisabler(self):
        state = bool(self.AbsorbEdit.currentIndex())

        self.ConcAEdit.setEnabled(state)
        self.ConcDEdit.setEnabled(state)
        # self.colorEdit.setEnabled(state)
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
        # self.AbsorberDisabler()

        self.AlignVEdit.setEnabled(state)
        self.AlignSEdit.setEnabled(state)
        self.AbsEmAngleEdit.setEnabled(state)

    def FixDistDisabler(self):
        state = bool(self.FixDistEdit.currentIndex())

        self.FixDistREdit.setEnabled(state)
        self.FixDistNeighborsEdit.setEnabled(state)


    def ClearMaterials(self,text = ''):
        self.NameEdit.setText(text)
        self.RefEdit.setText(text)
        self.ConcAEdit.setText(text)
        self.ConcDEdit.setText(text)
        self.colorEdit.setText(text)
        self.AlignVEdit.setText(str(text) + "," + str(text) + "," + str(text))
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
        # print self.TypeEdit.currentText()
        if self.TypeEdit.currentText() == "Mirror":
            self.MaterialEdit1.setEnabled(False)
            self.MaterialEdit2.setEnabled(False)
            self.FocalLengthEdit.setEnabled(False)
            self.DetSpecEdit.setEnabled(False)

        elif self.TypeEdit.currentText() == "Interface":
            self.MaterialEdit1.setEnabled(True)
            self.MaterialEdit2.setEnabled(True)
            self.FocalLengthEdit.setEnabled(False)
            self.DetSpecEdit.setEnabled(False)
        elif self.TypeEdit.currentText() == "Thin Lens":
            self.MaterialEdit1.setEnabled(False)
            self.MaterialEdit2.setEnabled(False)
            self.FocalLengthEdit.setEnabled(True)
            self.DetSpecEdit.setEnabled(False)
        elif self.TypeEdit.currentText() == "Dump":
            self.MaterialEdit1.setEnabled(False)
            self.MaterialEdit2.setEnabled(False)
            self.FocalLengthEdit.setEnabled(False)
            self.DetSpecEdit.setEnabled(True)

    def GeoEditDisabler(self):
        if self.GeoEdit.currentText() == "Plane":
            self.PositionEdit.setEnabled(True)
            self.DirectionEdit.setEnabled(True)
            self.SizeEdit.setEnabled(True)
            self.PointsEdit.setEnabled(False)
        elif self.GeoEdit.currentText() == "Circle":
            self.PositionEdit.setEnabled(True)
            self.DirectionEdit.setEnabled(True)
            self.SizeEdit.setEnabled(True)
            self.PointsEdit.setEnabled(False)
        elif self.GeoEdit.currentText() == "Cylinder":
            self.PositionEdit.setEnabled(True)
            self.DirectionEdit.setEnabled(True)
            self.SizeEdit.setEnabled(True)
            self.PointsEdit.setEnabled(False)
        elif self.GeoEdit.currentText() == "Triangle":
            self.PositionEdit.setEnabled(False)
            self.DirectionEdit.setEnabled(False)
            self.SizeEdit.setEnabled(False)
            self.PointsEdit.setEnabled(True)
        elif self.GeoEdit.currentText() == "Box":
            self.PositionEdit.setEnabled(True)
            self.DirectionEdit.setEnabled(True)
            self.SizeEdit.setEnabled(True)
            self.PointsEdit.setEnabled(False)



    def DrawCylinder(self,p1,p2,r,i):
        '''Draw a Cylinder from p1 to p2 with a radius r
        Adapted from here: http://stackoverflow.com/questions/32317247/how-to-draw-a-cylinder-using-matplotlib-along-length-of-point-x1-y1-and-x2-y2
        '''
        clr = i
        #find the vector from point p1 to point p2
        a = p2 -p1
        abs_a = np.linalg.norm(a)
        a /= abs_a

        #make some vector not in the same direction as v
        not_a = np.array([1, 0, 0])
        if (a == not_a).all():
            not_a = np.array([0, 1, 0])
        #make vector perpendicular to v
        n1 = np.cross(a, not_a)
        n1 /= np.linalg.norm(n1)
        #make unit vector perpendicular to v and n1
        n2 = np.cross(a, n1)

        t = np.linspace(0, abs_a, 6)
        theta = np.linspace(0, 2 * np.pi, 6)
        t, theta = np.meshgrid(t, theta)
        #generate coordinates for surface
        X, Y, Z = [p1[j] + a[j] * t + r * np.sin(theta) * n1[j] + r * np.cos(theta) * n2[j] for j in [0, 1, 2]]
        self.ax.plot_surface(X, Y, Z,color=self.colors[self.ElementProperties[clr][0]-1],alpha=self.alphas[self.ElementProperties[clr][0]-1])

    def DrawCircle(self,M,n,r,i):
        '''Draw Circle with radius r'''
        #make some vector not in the same direction as v
        not_n = np.array([1, 0, 0])
        if (n == not_n).all():
            not_n = np.array([0, 1, 0])
        #make vector perpendicular to v
        n1 = np.cross(n, not_n)
        n1 /= np.linalg.norm(n1)
        #make unit vector perpendicular to v and n1
        n2 = np.cross(n, n1)

        theta = np.linspace(0, 2 * np.pi, 100)
        X, Y, Z = [M[j] +  r * np.sin(theta) * n1[j] + r * np.cos(theta) * n2[j] for j in [0, 1, 2]]
        verts = [zip(X,Y,Z)]
        #color=self.colors[self.ElementProperties[i][0]-1],

        Collection = Poly3DCollection(verts,linewidths=1,alpha = self.alphas[self.ElementProperties[i][0]-1])
        FaceColor = self.colors[self.ElementProperties[i][0]-1]
        Collection.set_facecolor(FaceColor)
        Collection.set_edgecolor("black")
        self.ax.add_collection3d(Collection)
        ##self.ax.add_collection3d(Poly3DCollection(verts,facecolors=self.colors[self.ElementProperties[clr][0]-1],alpha=self.alphas[self.ElementProperties[clr][0]-1]))

    def DrawFigure(self):
        '''Draw Stuff in 3d,using polys, works way better'''
        plt.ion()
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111, projection='3d')
        # self.ax.set_autoscale_on(False)
        # self.ax.set_xlim([-1,1])
        # self.ax.set_ylim([-1,1])
        # self.ax.set_zlim([-1,1])
        self.ax.set_xlabel('x')
        self.ax.set_ylabel('y')
        self.ax.set_zlabel('z')
        i=0
        for el in self.Elements:
            #if element is plane
            if self.ElementGeometry[i] == 1:
                #If plane is not rotated around axis, find vectors in plane along axis
                if len(self.ElementSize[i])==2:
                    e1 = [-el[1][1],el[1][0],0]
                    e2 = [0,-el[1][2],el[1][1]]
                    ealt = [-el[1][2],0,el[1][0]]
                    if e1 == [0,0,0]:
                        e1 = ealt
                    elif e2 == [0,0,0]:
                        e2 = ealt
                    elif e1 == e2 or np.all(np.array(e1) == -np.array(e2)):
                        e2 = ealt
                else:
                #if rotated around axis
                    e1 = self.ElementSize[i][2]
                    e2 = np.cross(e1,el[1])

                e1 /= np.linalg.norm(e1)
                e2 /= np.linalg.norm(e2)


                a = el[0]+e1*self.ElementSize[i][0]+e2*self.ElementSize[i][1]
                b = el[0]+e1*self.ElementSize[i][0]-e2*self.ElementSize[i][1]
                c = el[0]-e1*self.ElementSize[i][0]-e2*self.ElementSize[i][1]
                d = el[0]-e1*self.ElementSize[i][0]+e2*self.ElementSize[i][1]
                x = [a[0],b[0],c[0],d[0]]
                y = [a[1],b[1],c[1],d[1]]
                z = [a[2],b[2],c[2],d[2]]

                verts = [zip(x,y,z)]
                #color=self.colors[self.ElementProperties[i][0]-1],
                print(verts)
                Collection = Poly3DCollection(verts,linewidths=1,alpha = self.alphas[self.ElementProperties[i][0]-1])
                FaceColor = self.colors[self.ElementProperties[i][0]-1]
                Collection.set_facecolor(FaceColor)
                Collection.set_edgecolor("black")
                self.ax.add_collection3d(Collection)
                #This Code works as well, Some versions of matplotlib do have a bug here though
                ##self.ax.add_collection3d(Poly3DCollection(verts,facecolors=self.colors[self.ElementProperties[i][0]-1],alpha=self.alphas[self.ElementProperties[i][0]-1]))

            if self.ElementGeometry[i] == 2:
                self.DrawCircle(el[0],el[1],self.ElementSize[i][0],i)
            if self.ElementGeometry[i] == 3:
                self.DrawCylinder(el[0],el[1],self.ElementSize[i][0],i)
            if self.ElementGeometry[i] == 4:
                x = [self.ElementSize[i][0][0],self.ElementSize[i][1][0],self.ElementSize[i][2][0]]
                y = [self.ElementSize[i][0][1],self.ElementSize[i][1][1],self.ElementSize[i][2][1]]
                z = [self.ElementSize[i][0][2],self.ElementSize[i][1][2],self.ElementSize[i][2][2]]
                verts = [zip(x,y,z)]
                Collection = Poly3DCollection(verts,linewidths=1,alpha = self.alphas[self.ElementProperties[i][0]-1])
                FaceColor = self.colors[self.ElementProperties[i][0]-1]
                Collection.set_facecolor(FaceColor)
                Collection.set_edgecolor("black")
                self.ax.add_collection3d(Collection)
                #This Code works as well, Some versions of matplotlib do have a bug here though
                ##self.ax.add_collection3d(Poly3DCollection(verts,facecolors=self.colors[self.ElementProperties[i][0]-1],alpha=self.alphas[self.ElementProperties[i][0]-1]))


            i+=1



if __name__ == '__main__':

    app = QApplication(sys.argv)
    ex = Example()
    sys.exit(app.exec_())
    main()
