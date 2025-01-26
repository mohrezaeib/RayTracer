import numpy as np
import scipy as sp
import copy as cp


def CheckSides(v1, v2, Plane):
    '''Checks if two vectors v1/v2 are on the same side of a plane'''
    if np.sign(Plane[1][0]*v1[0] + Plane[1][1]*v1[1] + Plane[1][2]*v1[2] - (Plane[0][0]*Plane[1][0] + Plane[0][1]*Plane[1][1] + Plane[0][2]*Plane[1][2])) == np.sign(Plane[1][0]*v2[0] + Plane[1][1]*v2[1] + Plane[1][2]*v2[2] - (Plane[0][0]*Plane[1][0] + Plane[0][1]*Plane[1][1] + Plane[0][2]*Plane[1][2])):
        return True
    else:
        return False



def FindClosestListElement(A, target):
    #A must be sorted
    idx = A.searchsorted(target)
    idx = np.clip(idx, 1, len(A)-1)
    left = A[idx-1]
    right = A[idx]
    idx -= target - left < right - target
    return idx


def MakeVectAxis(vector):
    '''Return transformation matrix that makes the given vector the new x-axis'''
    #make some vector not in the same direction as v
    not_v = np.array([1, 0, 0.0])
    if (vector == not_v).all():
        not_v = np.array([0, 1, 0.0])
    #make vector perpendicular to v
    n1 = np.cross(vector, not_v)
    n1 /= np.linalg.norm(n1)
    #make unit vector perpendicular to v and n1
    n2 = np.cross(vector, n1)
    n2 /= np.linalg.norm(n2)
    #Find Transformation matrix to make central line of cylinder one axis
    tm = np.zeros((3,3))
    tm[:,0] = vector
    tm[:,1] = n1
    tm[:,2] = n2
    return tm


def GenerateDirectionGauss(Photon):
    '''Generate vector with gaussian distributed angles from the z axis'''
    theta = Photon.CurrentMaterial.AlignS * np.random.randn()
    phi = np.random.rand()*2.0*np.pi
    v_m  = [np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi),np.cos(theta)]
    #make the alignment vector the z-vector [0,0,1]
    tm = np.roll(MakeVectAxis(Photon.CurrentMaterial.AlignVect),2,axis=1)
    v_m = np.dot(tm,v_m)
    v_m /= np.linalg.norm(v_m)
    return v_m

def GenerateDirectionSineSq():
    '''Sample from a sine squared distribution See: https://math.stackexchange.com/questions/498651/inverse-quantile-function-for-sin2x/499356'''
    # unirnd_x = np.random.rand() * (np.pi/2.0)
    # unirnd_y = np.random.rand()
    # sinsqsample = np.where(unirnd_y > np.sin(unirnd_x)**2,unirnd_x+np.pi/2.0,unirnd_x)

    '''For three dimensions it is needed to sample from the sin^3 distribution to weigh for the Spherical surface just life integration over spherical coordinates. This is hard thus here I used rejection sampling which is slow but works just fine'''

    unirnd_x = np.random.rand() * (np.pi)
    unirnd_y = np.random.rand()
    while unirnd_y > np.sin(unirnd_x)**3:
        unirnd_x = np.random.rand() * (np.pi)
        unirnd_y = np.random.rand()
    sinsqsample = unirnd_x

#_____________________________________________________
    # buff = 0
    # while buff == 0:
    #     unirnd = np.random.rand() + 0j
    #
    #     thrdrt1 = (2*unirnd+2*np.sqrt(unirnd-1)*np.sqrt(unirnd)-1)**(1/3.0)
    #     r3 = np.absolute(thrdrt1)
    #     th = np.angle(thrdrt1)
    #     thrdrt2 = r3 * np.exp(th*1j + 2j*np.pi/3)
    #     thrdrt3 = r3 * np.exp(th*1j - 2j*np.pi/3)
    #
    #
    #     term1 = np.arccos(thrdrt1+1/thrdrt1)
    #     term2 = np.arccos(thrdrt2+1/thrdrt2)
    #     term3 = np.arccos(thrdrt3+1/thrdrt3)
    #
    #
    #     a = np.array([term1,term2,term3])
    #     # print a,  a[np.absolute(np.imag(a))<1e-5]
    #     if len(a[np.absolute(np.imag(a))<1e-5])>0:
    #         buff = a[np.absolute(np.imag(a))<1e-5][0]
    #
    # sinsqsample = np.real(buff)
#_____________________________________________________

    psi= np.random.rand()*2*np.pi
    a = np.sin(sinsqsample)*np.cos(psi)
    b = np.sin(sinsqsample)*np.sin(psi)
    c = np.cos(sinsqsample)
    return [a,b,c]

def GenerateDirectionSine():
    '''generate random spherical directions See https://www.particleincell.com/2015/cosine-distribution/
    https://stackoverflow.com/questions/9750908/how-to-generate-a-unit-vector-pointing-in-a-random-direction-with-isotropic-dist'''

    v = np.random.randn(3)
    v /= np.linalg.norm(v)
    return v

def RotationMatrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians. https://stackoverflow.com/questions/6802577/python-rotation-of-3d-vector
    """
    axis = np.asarray(axis)
    axis = axis/np.sqrt(np.dot(axis, axis))
    a = np.cos(theta/2.0)
    b, c, d = -axis*np.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])


def MakeTriangle(p1,p2,p3):
    '''Create a triangle from three points, returns the center of mass as support vector and the normal vector'''
    p1,p2,p3 = np.array([p1,p2,p3])
    s = (p1+p2+p3)/3.0
    d = -np.cross(p3-p1,p2-p1)
    d /= np.linalg.norm(d)

    return [s,d]

def MakeBox(M=[0,0,0],Size=[0.2,0.2,0.2], n_x = [1,0,0]):
    '''Returns an array of arrays in the form [SupportVector,NormalVector] defining planes that define a Box'''
    #Define normal box
    a = np.array([ [[0.0,0.0,0],[1,0,0]] , [[0,0,0],[-1,0,0]] , [[0,0,0],[0,1,0]] , [[0,0,0],[0,-1,0]] , [[0,0,0],[0,0,1]] , [[0,0,0],[0,0,-1]] ])
    M = np.array(M)
    size = np.array(Size)
    n_x /= np.linalg.norm(n_x)
    tm = MakeVectAxis(n_x)
    #Transform vectors
    for i in range(np.shape(a)[0]):
        a[i][1] = np.dot(tm,a[i][1])
        a[i][1] = a[i][1]/np.linalg.norm(a[i][1])
        # print M + a[i][1] *size
        a[i][0] = M + a[i][1] *size[i/2]
    return a


def RectifySpectrum(Spectrum):
    '''Adjust spectrum for sampling. Interval widths have to be taken into account, also it is normalized to unity'''
    WorkingSpectrum = cp.deepcopy(Spectrum)
    Intervalwidths = np.ediff1d(Spectrum[:, 0])
    WorkingSpectrum[:, 1][1:] *= Intervalwidths
    WorkingSpectrum[:, 1] /= np.sum(WorkingSpectrum[:, 1])
    return WorkingSpectrum







##EOF
