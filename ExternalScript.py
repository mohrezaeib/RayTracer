import numpy as np

def main(Obj):

    print "External Script Run"
    #ZigZag for 2x2 Structure
    #Old
    #AddZigZag(Obj,-1, 0,1,-0.2,0.005, -0.04)
    #New
    #AddZigZag(Obj,-1, 0,1,-0.2,0.005, -(0.2/11))

    #AddZigZag(Obj,-1, 0.0,1,-0.2,0.005, -(0.2/11))

    AddZigZag(Obj,-10, -0.6,10,-0.8,0.005, -(0.2/11))

    #AddZigZagPlanes(Obj,-0.5,0.0,0.5,-0.05)
    #AddZigZagPlanes(Obj,-0.5,-0.0505,0.5,-0.1005)
    #AddZigZagPlanes(Obj,-0.5,-0.1010,0.5,-0.1510)
    #AddZigZagPlanes(Obj,-0.5,-0.1515,0.5,-0.2015)

    #Vert Planes for 2x2 Structure
    #AddVerticalPlanes(Obj,-1,1)

    #Vert Planes for 20x20 Structure
    #AddVerticalPlanes(Obj,-10,10)

    #Vert Planes for Hexagon
    #AddVerticalPlanes(Obj,-1.5,1.5)

    #AddCylinders(Obj, -0.0001, -0.1999,0.028,-1,1,-1,1,0.0025)

def AddCylinders(Obj,Top,Bottom,Distance,FromX,ToX,FromY,ToY,r=0.005):
    CurrentX = FromX + 0.5*Distance
    CurrentY = FromY + 0.5*Distance

    Obj.GeoEdit.setCurrentIndex(2)
    Obj.SizeEdit.setText(str(r))


    for i in range(2):
        CurrentX = FromX + (0.5+0.5*i)*Distance
        while CurrentX <= ToX-0.5*Distance:
            CurrentY = FromY + (0.5+0.5*i)*Distance
            while CurrentY <= ToY-0.5*Distance:
                #Obj.SizeEdit.setText(str(r))
                Obj.GeoEdit.setCurrentIndex(2)
                Obj.PositionEdit.setText(str(CurrentX)+","+str(CurrentY)+","+str(Top))
                Obj.DirectionEdit.setText(str(CurrentX)+","+str(CurrentY)+","+str(Bottom))
                Obj.AddElement()

                #Obj.SizeEdit.setText(str(r/2.0))
                Obj.GeoEdit.setCurrentIndex(1)
                Obj.PositionEdit.setText(str(CurrentX)+","+str(CurrentY)+","+str(Top))
                Obj.DirectionEdit.setText("0,0,1")
                Obj.AddElement()
                Obj.PositionEdit.setText(str(CurrentX)+","+str(CurrentY)+","+str(Bottom))
                Obj.DirectionEdit.setText("0,0,1")
                Obj.AddElement()

                CurrentY += Distance

            CurrentX += Distance

def AddVerticalPlanes(Obj,VariableFromX, VariabletoX, Step = 0.05 , Thickness = 0.005):
    CurrentPointX = VariableFromX + Step

    Width = Thickness/2.0

    Obj.GeoEdit.setCurrentIndex(4)
    Obj.SizeEdit.setText(str(Width) + ", 0.099, 0.99")
    # Obj.SizeEdit.setText(str(Width) + ", 0.0499, 9.99")

    while CurrentPointX <= VariabletoX:
        Obj.DirectionEdit.setText("1, 0, 0")
        Obj.PositionEdit.setText(str(CurrentPointX) + ",0,-0.7")
        Obj.AddElement()

        CurrentPointX += Step



def AddZigZagPlanes(Obj,VariableFromX,VariableFromZ,VariabletoX,VariableToZ, Step = 0.045 , Thickness = 0.005):

    CurrentPointX = VariableFromX + Step/2.0
    CurrentPointZ = VariableFromZ - Step/2.0
    sgn = 1

    Width = np.sqrt(2*(Step**2))

    Obj.GeoEdit.setCurrentIndex(0)
    Obj.SizeEdit.setText("0.5," + str(round(Width/2.0,4)))

    while CurrentPointX < (VariabletoX+Step):
        Obj.DirectionEdit.setText(str(sgn)+",0,1")
        Obj.PositionEdit.setText(str(CurrentPointX) + "," + str(0.0) + "," + str(CurrentPointZ))
        Obj.AddElement()
        CurrentPointZ -= Thickness
        Obj.PositionEdit.setText(str(CurrentPointX) + "," + str(0.0) + "," + str(CurrentPointZ))
        Obj.AddElement()
        CurrentPointZ += Thickness
        sgn *= -1

        CurrentPointX += Step



def AddZigZag(Obj,VariableFromX, VariableFromZ, VariableToX, VariableToZ, StepX, StepZ):
    '''            szed = map(float,self.SizeEdit.text().split(","))
                drv = map(float,self.DirectionEdit.text().split(","))
                bx = self.MakeBox(map(float,self.PositionEdit.text().split(",")), szed , drv)

self.Elements = np.append(self.Elements,[[map(float,self.PositionEdit.text().split(",")),map(float,self.DirectionEdit.text().split(","))]],axis=0)

                '''
    CurrentPointX = VariableFromX+StepX/2.0
    CurrentPointZ = VariableFromZ+StepZ/2.0
    sgn = -1

    Obj.DirectionEdit.setText("1, 0, 0")
    Obj.SizeEdit.setText("0.0025,0.00909,9.99")

    while CurrentPointX < VariableToX:
        Obj.GeoEdit.setCurrentIndex(4)
        # print CurrentPointX, CurrentPointZ,Step
        print str(CurrentPointX) + "," + str(0.0) + "," + str(CurrentPointZ)
        Obj.PositionEdit.setText(str(CurrentPointX) + "," + str(0.0) + "," + str(CurrentPointZ))
        Obj.AddElement()

        if CurrentPointZ > VariableFromZ + StepZ:
            sgn *= -1
            # print "sgnchng",CurrentPointZ , VariableFromZ
        elif CurrentPointZ < VariableToZ - StepZ:
            sgn *= -1
            # print "sgnchng",CurrentPointZ , VariableToZ

        CurrentPointX += StepX
        CurrentPointZ += sgn*StepZ



#
