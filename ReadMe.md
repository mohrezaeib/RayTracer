# RayTracer: A Comprehensive Overview

This repository contains a Python-based **Ray Tracing** and **FRET (Förster Resonance Energy Transfer)** simulation environment. It is designed to:

1. **Parse** user-defined geometry and material properties from external files.
2. **Simulate** the propagation of photons through these geometries, accounting for reflection, refraction, and absorption.
3. **Handle** specialized processes such as FRET, thin-lens focusing, aligned dipoles, and partial polarization.

Below, you will find an extensive explanation of the underlying physics and math, details about how the code is structured, and examples demonstrating typical use cases.

---

## 1. Introduction to the Software

**RayTracer** is both a general-purpose optical ray tracer and a specialized simulator for **FRET** phenomena. Users can:

- **Define optical elements** (mirrors, lenses, interfaces, dumps/detectors, etc.) by specifying their geometry (plane, circle, cylinder, triangle, box) and their properties (material indices, absorption spectra, reflectivity, etc.).
- **Create or load materials** with given refractive indices, absorption/emission spectra (for dyes or fluorophores), alignment parameters, and other properties relevant to FRET.
- **Generate photons** from user-defined emission distributions (directional isotropic, solar spectrum, a specified wavelength, etc.).
- **Run a simulation** that tracks how each photon interacts with the geometry:  
  - **Intersections** with surfaces  
  - **Reflections** (mirror-like)  
  - **Refractions** (Snell’s law at interfaces)  
  - **Absorption** (Lambert-Beer law)  
  - **FRET events** (non-radiative energy transfer between donor and acceptor dyes)  

After running, the software generates logs, including details about absorbed photons, detectors (dumps), or any custom user-specified output.

---

## 2. Physics and Mathematics in RayTracer

### 2.1 Ray Tracing Basics

1. **Ray Definition**:  
   A photon is treated as a classical ray characterized by a `Location` (3D coordinates) and a `Direction` (unit vector).  
   \[
   \mathbf{r}(t) = \mathbf{r}_0 + t \, \mathbf{d}, 
   \]
   where \(\mathbf{r}_0\) is the starting point and \(\mathbf{d}\) is the direction vector (normalized).

2. **Intersection with Geometry**:  
   - **Plane**: A plane can be defined by a point \(\mathbf{p}\) and a normal \(\mathbf{n}\). The intersection parameter \( t \) solves  
     \[
     (\mathbf{p} - \mathbf{r}_0)\cdot\mathbf{n} = t \, (\mathbf{d}\cdot \mathbf{n}).
     \]
     If \( \mathbf{d}\cdot \mathbf{n} = 0 \), the ray is parallel and either does not intersect or lies in the plane.  
   - **Circle**, **Rectangle**, **Triangle**: These are “bounded planes.” A solution for the plane intersection is found, and then a check ensures the intersection is within the circle’s radius, rectangle’s edges, or triangle’s boundaries.  
   - **Cylinder**: This involves solving a quadratic for the intersection along the cylinder axis, then clipping by cylinder height.  

3. **Reflection and Refraction**:  
   - **Reflection** uses the standard formula:  
     \[
     \mathbf{d}_{\text{reflected}} = \mathbf{d}_{\text{inc}} - 2\,(\mathbf{d}_{\text{inc}}\cdot \mathbf{n})\,\mathbf{n},
     \]
     assuming a perfectly specular mirror.  
   - **Refraction** uses Snell’s law:  
     \[
     n_1 \sin \theta_i = n_2 \sin \theta_t,
     \]
     implemented via vector projections, checking for total internal reflection when the expression under the square root becomes negative.  

4. **Lambert-Beer Absorption**:  
   - If the material is an “absorber,” the photon has a finite “free path” according to the exponential distribution. The code calculates a random absorption distance:  
     \[
     d_{\text{abs}} = -\frac{\ln(\zeta)}{\alpha},
     \]
     where \(\zeta\) is a uniform random number in \([0,1]\) and \(\alpha\) is the absorption coefficient derived from the material’s spectra, concentration, and the Beer-Lambert law.  

### 2.2 FRET (Förster Resonance Energy Transfer)

1. **Förster Distance \( R_0 \)**:  
   Defines the distance at which there is a 50% energy transfer efficiency from donor to acceptor. Computed from donor emission, acceptor absorption spectra, quantum yield, etc. 
   \[
   R_0^6 \propto \frac{9 \ln(10)\,Q_D}{128 \pi^5 \,N_A \, n^4}\,\int F_D(\lambda)\,\epsilon_A(\lambda)\,\lambda^4\, d\lambda,
   \]
   with \(\epsilon_A\) the acceptor’s molar absorptivity, \(F_D\) the donor’s normalized emission, \(N_A\) Avogadro’s number, \(n\) refractive index, etc.

2. **Probability of Energy Transfer**:  
   - If the donor and acceptor are at distance \( r \), the probability of FRET can be computed as  
     \[
     P_{\text{FRET}} = \frac{1}{1 + (r/R_0)^6}.
     \]
   - Code also supports **“fixed distance”** or **“variable concentrations”** approaches, allowing Monte Carlo selection of donor/acceptor collisions.  

3. **Alignment (Polarization)**:  
   - Some dyes can be aligned (e.g., polymer films). The code can incorporate a distribution for dipole orientations, such as a Gaussian distribution around a preferred orientation.  
   - Emission can be partially polarized.  

### 2.3 Emission Spectra and Polarization

- **Emission Sine** and **Sine-Squared**: In certain calculations, random directions are generated by sampling from \( \sin(\theta) \) or \( \sin^2(\theta) \) distributions.  
- **Rectified Spectra**: The software can read text files describing absorption/emission spectra. These are integrated, scaled, and used to stochastically pick emission or absorption wavelengths.

---

## 3. Implementation Details

### 3.1 Code Structure

- **`RayTrace.py`**  
  Main class orchestrating the simulation. It loads geometry, materials, and photons from saved “.prt” files. Then it runs the simulation loop:
  1. **Check Elements For Crossing Point**: Finds the nearest intersection among all geometry.  
  2. **Check For Absorption**: If the photon is within an absorbing medium, might be absorbed (Lambert-Beer).  
  3. **Execute Operation**: If intersection is valid, decide if we reflect, refract, pass, or do lens focusing.  

- **`PhotonObject.py`**  
  Represents a photon with attributes:
  - `Location`, `Direction`, `Wavelength`, `CurrentMaterial`  
  - `ReflectOnPlane()`, `RefractOnPlane()`, `FocusThinLens()`, etc.  
  - Tracks distances left to absorption and polarization vectors.

- **`Elements.py`**  
  Contains classes for optical elements:  
  - `Mirror()`, `Interface()`, `Lens()`, `Dump()`, etc.  
  Each has an `ExecuteOperation()` method describing how it modifies the photon.  

- **`Geometries.py`**  
  Contains geometry classes:  
  - `Rectangle()`, `Circle()`, `Cylinder()`, `Triangle()`  
  Each has a `FindIntersections()` method to compute intersection points with the current photon.  

- **`MonteCarloFRET.py`**  
  Describes how FRET probabilities are computed, how the photon might be re-emitted at a new wavelength, or how energy is lost to internal conversion.

- **`DrawGraph.py`**  
  Utility class for 3D plotting using `matplotlib`. It can visualize cylinders, circles, rectangles, and the photon’s trajectory lines.

- **`RTUtils.py`**  
  Collection of helper functions: generating random directions, rotation matrices, checking sides of planes, etc.

- **`IOOperation.py`**  
  Handles saving and loading `.prt` files (pickle) for elements, materials, and photon generation parameters.

### 3.2 Typical Workflow

1. **Configure Materials and Geometry**  
   - Create or edit `.prt` files specifying arrays for “Elements,” “ElementProperties,” “ElementGeometry,” “ElementSize,” and “Materials.”  
   - The parser will read these into arrays.

2. **Photon Setup**  
   - Set the number of photons to simulate, starting positions, directional distributions, emission wavelength type, etc.

3. **Run the Simulation**  
   ```python
   myRayTracer = RayTrace.RayTrace()
   myRayTracer.ParsePath("./MyJob/")  # Where the .prt files are stored
   myRayTracer.runSimulation()
   ```
   - The simulation loop tracks each photon until it gets dumped, internally converted, or is lost after many reflections.

4. **Logs and Output**  
   - After finishing, `IOOperation.SaveLog()` is called, writing results such as `Logfile.txt`, spectra for detectors, `AbsorbedPhotons.txt`, etc.

---

## 4. Usage Examples

### 4.1 Minimal Example for Reflection

**Objective**: Create a mirror plane and a single photon hitting it.

1. **Geometry**:  
   - One plane with type=Mirror and normal along z-axis at `z=0`.
2. **Photon**:  
   - Starting at `[0, 0, 1]`, pointing downward `[0, 0, -1]`.

**Code Sketch**:
```python
import RayTrace, Materials, PhotonObject

# Create a RayTrace instance
rt = RayTrace.RayTrace()

# Manually define a single mirror plane (or load from .prt)
rt.Elements = [
    [[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]]  # plane support and normal
]
rt.ElementProperties = [[1]]  # '1' => Mirror
rt.ElementGeometry = [1]      # '1' => Rectangle (infinite for test)
rt.ElementSize = [[1.0, 1.0]] # bounding rectangle size (±1 in each direction)
# Simple material
air = Materials.MaterialObj("air", 1.0)
rt.MaterialDict["air"] = air

# A single photon
rt.NumberOfPhotons = 1
rt.Photon = PhotonObject.Photon(Location=[0,0,1], EmissionMode="Directional Strict",
                                EmissionDirection=[0,0,-1], EmissionWavelengthType="Specified Wavelength...",
                                StartingWavelength=500.0, StartingMaterial=air)

# Start simulation
rt.CreateOpticalElementList()
rt.runSimulation()
```
You’ll see logs for how the photon intersects `z=0` plane and is reflected.

### 4.2 FRET Example with Two Dyes

- We define two materials: `donorMat` and `acceptorMat`, specifying their absorption/emission spectra, FRET parameters, etc.
- The geometry might be a block of “mixed doping,” or we can fix a distance for testing.

**Key Steps**:

1. In `Materials.prt`, set:  
   ```python
   # Format: [name, refr_index, IsAbsorber, conc_b, conc_a, color, IsAligned, ... , sp_abs_D, sp_em_D, sp_abs_A, sp_em_A, ...]
   donorMat = [
       "donor", 1.5, True, 5e-3, 0.0, "blue", False, [0,0,1], 0.1, 6.0, 0.8, 0.0, 10.0, 0.0,
       23500, 0, "DonorAbs.txt", "DonorEm.txt", "None.txt", "None.txt", 0, False, 1.0, 0
   ]
   acceptorMat = [
       "acceptor", 1.5, True, 0.0, 2e-3, "red", False, [0,0,1], 0.1, 5.0, 0.0, 0.9, 0.0, 8.0,
       0, 54000, "None.txt", "None.txt", "AcceptorAbs.txt", "AcceptorEm.txt", 0, False, 1.0, 0
   ]
   ```
2. In your main Python script, load these materials. The code automatically sets up the FRET logic from `MonteCarloFRET.py`.
3. Run the simulation, track how many photons get “transferred” from donor to acceptor, or lost to internal conversion.

---

## 5. Further Reading and References

- **For basic ray tracing**:  
  - [“Physically Based Rendering” (Pharr et al.)](https://www.pbrt.org/) for an in-depth, physically-based approach.
- **For FRET**:  
  - The original works of **Theodor Förster** (1959).  
  - P. Wu and L. Brand, *Resonance Energy Transfer: Methods and Applications*, Analytical Biochemistry, 200, 1989.

---

## 6. Conclusion

**RayTracer** provides a flexible framework for simulating optical phenomena combined with **FRET** processes. By mixing ray optics (geometrical intersection, reflection, refraction) with the advanced physics of excited-state energy transfer, it is well-suited for:

- **Polymer waveguides** doping analysis  
- **Biophotonics** applications (fluorescent molecules, quantum yield studies)  
- **Optical design** of mirrors, lenses, and high-level system simulations  

The modular structure (distinct classes for geometry, elements, photons, materials) makes it easy to **extend** or **customize** for specific needs, whether you want to add a new geometry type or incorporate specialized spectroscopic data.

For **getting started**, refer to:
1. The example script in `Main.py` or `RayTrace.py`.
2. The `.prt` files in the `Job` folder, showing how to define geometry, materials, and photons.
3. The docstrings throughout the code for details on each method’s usage.

Feel free to reach out via issues or pull requests if you have questions, find bugs, or want to contribute new features!