import numpy as np
import scipy as sp
import copy as cp


def CheckSides(v1, v2, Plane):
    """
    Checks if two points/vectors (v1, v2) lie on the same side of a given plane.

    Plane is in the form:
      Plane[0] = a point on the plane, e.g. (px, py, pz)
      Plane[1] = a normal vector to the plane, e.g. (nx, ny, nz)
    """
    # Evaluate plane equation for v1 and v2, check signs
    lhs1 = (Plane[1][0]*v1[0]
            + Plane[1][1]*v1[1]
            + Plane[1][2]*v1[2]
            - (Plane[0][0]*Plane[1][0]
               + Plane[0][1]*Plane[1][1]
               + Plane[0][2]*Plane[1][2]))

    lhs2 = (Plane[1][0]*v2[0]
            + Plane[1][1]*v2[1]
            + Plane[1][2]*v2[2]
            - (Plane[0][0]*Plane[1][0]
               + Plane[0][1]*Plane[1][1]
               + Plane[0][2]*Plane[1][2]))

    return np.sign(lhs1) == np.sign(lhs2)


def FindClosestListElement(A, target):
    """
    Finds the index of the closest value in array A to the given target.
    The array A must be sorted. Uses np.searchsorted for efficiency.
    """
    idx = A.searchsorted(target)
    idx = np.clip(idx, 1, len(A) - 1)
    left = A[idx - 1]
    right = A[idx]
    # If the target is closer to left than right, decrement idx
    idx -= (target - left) < (right - target)
    return idx


def MakeVectAxis(vector):
    """
    Return a 3x3 transformation matrix whose first column is `vector`
    (normalized) and whose other two columns are perpendicular directions.

    The matrix, when multiplied by a vector in local coordinates, aligns
    the local x-axis with `vector` in global coordinates.
    """
    vector = np.array(vector, dtype=float)
    vector /= np.linalg.norm(vector)

    not_v = np.array([1.0, 0.0, 0.0])
    if (vector == not_v).all():
        not_v = np.array([0.0, 1.0, 0.0])

    n1 = np.cross(vector, not_v)
    n1 /= np.linalg.norm(n1)
    n2 = np.cross(vector, n1)
    n2 /= np.linalg.norm(n2)

    tm = np.zeros((3, 3), dtype=float)
    # First column is our main axis
    tm[:, 0] = vector
    # Second column is perpendicular axis
    tm[:, 1] = n1
    # Third column is final perpendicular axis
    tm[:, 2] = n2

    return tm


def GenerateDirectionGauss(Photon):
    """
    Generate a random direction vector with an angular spread given by
    Photon.CurrentMaterial.AlignS (sigma). This is effectively a small-angle
    Gaussian distribution around the alignment vector.

    Returns a unit 3D vector.
    """
    theta = Photon.CurrentMaterial.AlignS * np.random.randn()
    phi = np.random.rand() * 2.0 * np.pi
    v_m = [
        np.sin(theta) * np.cos(phi),
        np.sin(theta) * np.sin(phi),
        np.cos(theta)
    ]
    # Make the alignment vector the z-axis
    tm = np.roll(MakeVectAxis(Photon.CurrentMaterial.AlignVect), 2, axis=1)
    v_m = np.dot(tm, v_m)
    v_m /= np.linalg.norm(v_m)
    return v_m


def GenerateDirectionSineSq():
    """
    Generates a random 3D direction with a sin^3 weighting (approx),
    using a simple rejection sampling approach. Then orients it randomly
    in phi.

    This is a custom distribution, often referred to for certain emission
    patterns. For truly isotropic, see GenerateDirectionSine.
    """
    while True:
        unirnd_x = np.random.rand() * np.pi
        unirnd_y = np.random.rand()
        # Rejection test with sin^3
        if unirnd_y <= (np.sin(unirnd_x)**3):
            break

    sinsqsample = unirnd_x
    psi = np.random.rand() * 2.0 * np.pi
    a = np.sin(sinsqsample) * np.cos(psi)
    b = np.sin(sinsqsample) * np.sin(psi)
    c = np.cos(sinsqsample)
    return [a, b, c]


def GenerateDirectionSine():
    """
    Generates a random 3D unit vector with an isotropic (spherical) distribution.
    Uses the standard normal method and normalizes the result.
    """
    v = np.random.randn(3)
    v /= np.linalg.norm(v)
    return v


def RotationMatrix(axis, theta):
    """
    Return the rotation matrix associated with a rotation about
    the given axis by theta radians (right-hand rule).
    Source:
    https://stackoverflow.com/questions/6802577/python-rotation-of-3d-vector
    """
    axis = np.asarray(axis, dtype=float)
    axis /= np.sqrt(np.dot(axis, axis))
    a = np.cos(theta / 2.0)
    b, c, d = -axis * np.sin(theta / 2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad = b*c, a*d
    ac, ab = a*c, a*b
    bd, cd = b*d, c*d
    return np.array([
        [aa + bb - cc - dd, 2 * (bc + ad),     2 * (bd - ac)],
        [2 * (bc - ad),     aa + cc - bb - dd, 2 * (cd + ab)],
        [2 * (bd + ac),     2 * (cd - ab),     aa + dd - bb - cc]
    ])


def MakeTriangle(p1, p2, p3):
    """
    Given three points p1, p2, p3 in 3D, returns:
      [ center_of_mass, normal_vector ]
    where center_of_mass is the average of p1, p2, p3
    and normal_vector is a unit normal to the plane of the triangle.
    """
    p1 = np.array(p1, dtype=float)
    p2 = np.array(p2, dtype=float)
    p3 = np.array(p3, dtype=float)
    s = (p1 + p2 + p3) / 3.0
    d = -np.cross(p3 - p1, p2 - p1)
    d /= np.linalg.norm(d)
    return [s, d]


def MakeBox(M=[0, 0, 0], Size=[0.2, 0.2, 0.2], n_x=[1, 0, 0]):
    """
    Returns an array of 6 planes in the form [SupportVector, NormalVector]
    that define a rectangular box oriented by n_x, with 'Size' as half-edges
    in each dimension. The planes are arranged in pairs around M.
    """
    # Basic 6-plane box, oriented along standard axes
    a = np.array([
        [[0.0, 0.0, 0.0], [1, 0, 0]],
        [[0.0, 0.0, 0.0], [-1, 0, 0]],
        [[0.0, 0.0, 0.0], [0, 1, 0]],
        [[0.0, 0.0, 0.0], [0, -1, 0]],
        [[0.0, 0.0, 0.0], [0, 0, 1]],
        [[0.0, 0.0, 0.0], [0, 0, -1]]
    ], dtype=float)
    M = np.array(M, dtype=float)
    size = np.array(Size, dtype=float)
    n_x = np.array(n_x, dtype=float)
    n_x /= np.linalg.norm(n_x)
    tm = MakeVectAxis(n_x)

    for i in range(a.shape[0]):
        # Rotate the normal by tm
        a[i][1] = np.dot(tm, a[i][1])
        a[i][1] /= np.linalg.norm(a[i][1])
        # Adjust support vector by half-size along that normal
        # In Python 3, i//2 ensures integer division for half-plane indexing
        a[i][0] = M + a[i][1] * size[i // 2]

    return a


def RectifySpectrum(Spectrum):
    """
    Takes a 2D array Spectrum of shape (N,2): [ [wavelength, intensity], ... ]
    and returns a new array where intensity is adjusted by the interval width
    and normalized so the sum of intensities = 1.0.

    Assumes Spectrum[:, 0] is strictly increasing in wavelength.
    """
    WorkingSpectrum = cp.deepcopy(Spectrum)
    # Interval widths for each step in wavelengths
    Intervalwidths = np.ediff1d(Spectrum[:, 0])
    # Multiply intensities by interval widths (except the first intensity)
    WorkingSpectrum[1:, 1] *= Intervalwidths

    # Normalize so total area = 1
    total = np.sum(WorkingSpectrum[:, 1])
    if total > 0:
        WorkingSpectrum[:, 1] /= total

    return WorkingSpectrum
