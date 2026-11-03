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
is available, but the spatial solver in this increment is static.

The model assumes infinitesimal strain, linear reversible constitutive laws,
fixed geometry and electroquasistatic fields inside the meshed dielectric.
It does not add finite-device resonances, transient circuit loading, exterior
fringing, loss, hysteresis, depoling, breakdown or thermal coefficient evolution.
The existing longitudinal circuit/thermal transducer remains available for
those already-supported reduced coupling paths. Synthetic G1/G3 regressions
test the numerical equations and frame/boundary conventions; they do not
establish experimental accuracy for PZT-5A or any manufactured device.
