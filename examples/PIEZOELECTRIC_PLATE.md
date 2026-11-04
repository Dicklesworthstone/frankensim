# Three-dimensional piezoelectric plate

Run the production source, dielectric, solid and coupled solver together:

```sh
cargo run -p fs-couple --example piezoelectric_plate
```

The example builds a 10 × 10 × 2 mm plate with 48 P1 tetrahedra, uniform
electrodes on its broad faces and six supports that remove rigid motion while
allowing free transverse and thickness strains. A 100 V voltage difference
produces displacement and electrode charge through one coupled solve. A second
solve replaces the driven electrode by an isolated conductor with the first
solve's total charge: its voltage is then an unknown shared by the entire
electrode. Both solves include electrical feedback on mechanical deformation.

The example then starts the same finite plate at rest and ramps the voltage
difference from 0 to 100 V over 400 ns, in 20 ns steps. The consistent native
mass matrix now evolves its nodal velocities and displacement while the
dielectric equation solves the simultaneous electric field. After the ramp,
the driven electrode is disconnected at its **computed charge**. The plate
keeps its accepted displacement and velocity, and the isolated sensor voltage
changes with its subsequent free motion. The printed samples include thickness
change, voltage and accumulated energy drift. The timestep and coarse mesh are
illustrative numerical inputs, not a converged prediction of a physical device.
A second trajectory starts from that identical accepted state with a declared
1 kΩ external resistor. Its current changes the electrode charge inside the
same coupled solve, changing the mechanical motion as well. The final report
shows the resulting voltage, deformation, resistor heat and complete energy
defect. The resistor is an explicitly authored Ohmic load, with resistivity
`1e-3 Ω m`, length `10 mm` and area `0.01 mm²`; it is not a measured PZT loss
or a temperature-dependent resistor model.

## Data and conventions

The source is PIEZO.COM's
[Materials Technical Data, typical values, 2020-11-12](https://info.piezo.com/hubfs/Data-Sheets/piezo-material-properties-data-sheet-20201112.pdf),
PZT-5A column. This example uses ten independent compliance/charge/dielectric
coefficients and density, and explicitly **assumes an ideal transversely
isotropic poled ceramic** to fill the tensor. The completed tensor is an
authored model derived from nominal published values, not 64 independent
measurements or a calibrated specimen model.

| Input | Value | Convention |
| --- | --- | --- |
| `s11`, `s12`, `s13`, `s33`, `s55` | 15.1, −4.8, −7.6, 18.6, 40.0 × 10⁻¹² Pa⁻¹ | Fixed electric field; engineering shear |
| `d31`, `d33`, `d15` | −190, 390, 460 × 10⁻¹² m/V | Reciprocal strain/charge coefficients |
| Relative `epsilon11^T`, `epsilon33^T` | 1600, 1800 | Fixed stress; nominal dielectric data at 1 kHz |
| Density | 7950 kg/m³ | Reference solid |

The code converts engineering shear to Mandel coordinates explicitly, converts
relative to absolute permittivity, and resolves a full same-card source bundle.
The material owner then computes

```text
C^E = inverse(s^E)
e = d C^E
epsilon^S = epsilon^T - d C^E d^T
```

Admission requires the full coupled storage to remain stable. The independent
source domain and provenance remain attached to the derived law. This example
uses a frozen quasistatic approximation to the nominal data; it does not invent
a measured temperature, frequency or field validity interval.

## Expected physical behavior

For a stress-free uniform plate the thickness extension is `d33 × 100 V`,
or **39 nm**. The lateral contraction follows `d31` from the same tensor.
Free capacitance is `epsilon33^T × area / thickness`, about **0.797 nF**;
using the fully clamped permittivity here would give the wrong free response.
Returning the computed total charge to an isolated electrode recovers its
original −100 V potential. The printed stored energy must close against
mechanical and electrical work within the declared tolerance.

Local acoustic phase speeds are also printed for fixed electric field and
charge-free insulating bulk. For the latter, Gauss's law removes only the
normal electric displacement: `n·D=0` with `E` parallel to `n`. Using the full
local `D=0` stiffness instead would incorrectly stiffen some shear modes.

## Scope

The reusable API is `fs_couple::piezoelectric::spatial::TetPiezoelectricModel`.
It accepts one tensor or one tensor per tetrahedron, explicit material/mesh
frames, nonzero prescribed displacements, nodal force and volume-charge loads,
and prescribed-voltage or total-charge electrodes. The sparse operators are
prepared once for repeated load cases. The full native consistent mass matrix
also drives `spatial::transient::TetPiezoelectricTransient`. Its implicit
midpoint solve uses `K + 4 M / dt²`, the same reciprocal coupling matrix and
the original electrostatic operator. Electrical potential has no artificial
inertia. Fixed supports and electrode topology are prepared once; force is a
midpoint input and prescribed charge or voltage varies between endpoints.

Every accepted step independently checks momentum and charge balance, as well
as the change in kinetic, strain and electric storage against actual mechanical
and electrical work. The closed isolated response has no energy source or loss.
Initialization solves the electric field at the supplied displacement; it does
not reset the mechanical state to static equilibrium. Candidate construction
leaves the original state intact if a solve, material envelope or cancellation
check fails.

For loaded devices, `transient::load::step_with_loads_candidate` accepts
source-bound resistors and voltage or current drives on the retained electrode
types. The additional conductance enters the original sparse solve; there is
no separate load iteration or substituted mechanical damping. A load report
separates external electrical supply from native Ohmic heat while keeping the
original conductor source attached. The ordinary candidate's work report has
the smaller accounting boundary of the body alone. Resistors must connect
electrodes in the same retained electrical reference component.

The same full tensor also supplies the existing reduced rod and circuit models
through `TensorPiezoelectric::longitudinal(direction)`. This explicitly imposes
uniaxial stress and a longitudinal electric field while retaining transverse
mechanical relaxation. It keeps the original source claims and a receipt naming
the tensor, frame and oriented cut. It does not impose insulating side faces or
solve fringing fields. The existing scalar thermal-coefficient model refuses
tensor reductions until a tensor-aware thermal law is supplied.

The model assumes infinitesimal strain, linear reversible constitutive laws,
fixed geometry and electroquasistatic fields inside the meshed dielectric.
The spatial transient resolves the modes retained by its mesh, with temporal
dispersion from the chosen timestep. It supports the described resistive loads,
but does not add general RLC circuits, exterior fringing, acoustic radiation,
internal dielectric/structural loss, hysteresis, depoling, breakdown or thermal
coefficient evolution. Resistor dissipation is reported at the original load;
the example does not evolve resistor temperature or add its heat to the ceramic.
The existing longitudinal circuit/thermal transducer remains available for
those already-supported reduced coupling paths. Synthetic G1/G3 regressions
test the numerical equations and frame/boundary conventions; they do not
establish experimental accuracy for PZT-5A or any manufactured device.
