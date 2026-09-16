# Directional and temperature-dependent cooling materials

`cooling-network` accepts explicit material laws in the existing
`solid.materials` table. Every tetrahedron still needs exactly one entry in
`solid.element_materials`, referencing a declared material name. The legacy
uniform `solid.conductivity_w_m_k` mode is unchanged; do not combine it with a
material table.

Run the fan-cooled hotspot with an oriented substrate and an isotropic spreader:

```bash
frankensim --json cooling-network \
  examples/cooling-network/orthotropic-hotspot.json
```

This example uses illustrative geometry, conductivity, heating, fan and channel
inputs. It is not a measured PCB, a laminate homogenization, or a validated
hardware design.

## Constitutive choices

Each named material has `name`, `source`, and exactly ONE conductivity form.
All conductivity values are in W/(m K); coordinates are the mesh's Cartesian
coordinates. No material frame is inferred from geometry or from a name.

An existing isotropic material remains:

```json
{"name":"metal","conductivity_w_m_k":20,"source":"caller declaration"}
```

A full constant tensor must be finite, symmetric and positive definite under
the existing `fs-conduction` admission rules:

```json
{
  "name":"substrate",
  "conductivity_tensor_w_m_k":[[8.48,8.64,0],[8.64,13.52,0],[0,0,1]],
  "source":"caller-declared tensor in mesh coordinates"
}
```

Alternatively declare the three principal directions and conductivities:

```json
{
  "name":"substrate",
  "orthotropic":{
    "principal_axes":[[0.6,0.8,0],[-0.8,0.6,0],[0,0,1]],
    "conductivity_w_m_k":[20,2,1]
  },
  "source":"caller-declared principal axes in mesh coordinates"
}
```

The matrix ROWS are unit principal directions in the mesh frame, mutually
orthogonal. Entry `i` in the conductivity vector belongs to axis row `i`.
`fs-conduction` constructs `K = sum_i k_i e_i e_i^T`; off-diagonal terms are
retained. It does not average the three conductivities into a scalar. The two
directional examples above describe the same tensor, up to floating-point
rounding. Results preserve the declaration, coordinate convention and, for an
orthotropic declaration, the resolved tensor.

Constant directional materials use the existing steady and backward-Euler
assembly paths, including element assignments and matching-face contacts.
The existing inlet/effective-convection-coefficient adjoints use that same
operator. No material-parameter, orientation or geometry adjoint is added.

## Bounded nonlinear scalar conductivity

A steady material can instead declare a piecewise-linear `k(T)` law:

```json
{
  "name":"temperature-sensitive-solid",
  "conductivity_curve":{
    "temperature_k":[250,400],
    "conductivity_w_m_k":[17,2]
  },
  "source":"illustrative bounded k(T), not a measured material card"
}
```

The parallel arrays need two to 256 entries. Temperatures are positive kelvin
in strictly increasing order; conductivities are finite and strictly positive.
The existing `ConductivityTable` owns interpolation and refuses extrapolation.
The initial guess and subsequent solver evaluations must remain within the
material's supported temperature span; no endpoint clamping is performed.

The full temperature-dependent law reaches the production nonlinear solver.
Its implicit tangent and adjoint include `K'(T)` through the existing
`RobinLinearization`; this is not an adjoint of a frozen representative
conductivity. Existing derivative refusals at slope discontinuities or validity
endpoints remain in force. Multiple material laws may be assigned to different
elements, with source and contact terms retained.

**Nonlinear conductivity is steady-only in this product.** The existing
backward-Euler producer rejects temperature-dependent conductivity. A transient
request is not silently evaluated using a frozen curve. Fluid properties and
convection coefficients remain frozen during each thermal solve in either mode.

## Material uncertainty through real cooling solves

Run the corresponding material-uncertainty example:

```bash
frankensim --json cooling-network-uq \
  examples/cooling-network/orthotropic-hotspot.json \
  examples/cooling-network/uq-orthotropic-hotspot.json \
  --checkpoint material-study.uqcp
```

The destination must be new. Existing checkpoint/resume, sample-chunk and
sequential-compliance controls described in `COOLING_UQ.md` apply unchanged.
Each sample modifies the real material declaration and invokes the normal
cooling producer. The new target forms are:

```json
{"kind":"solid-conductivity"}
{"kind":"material-conductivity","material":"spreader"}
{"kind":"material-principal-conductivity","material":"substrate","axis":0}
```

The first requires a uniform scalar solid; the second a named isotropic scalar
material. The third changes only one explicitly declared orthotropic principal
conductivity, with zero-based axis 0, 1 or 2. The axes, other principal values
and element assignments remain fixed. The units are W/(m K).

Positive uniform supports and the existing explicit dependence models are
available. Gaussian samples leaving the positive physical domain refuse the
entire execution; they are not clipped, redrawn or skipped. A scalar uncertainty
target cannot overwrite an anisotropic tensor or a temperature-dependent law.
Independent tensor-entry sampling, orientation uncertainty and uncertain curve
knots are not implemented. A fixed `k(T)` law may coexist with other supported
uncertain parameters in a steady base model.

## Focused checks

```bash
cargo test -p fs-cli --bin frankensim network_command::solid_data
cargo test -p fs-cli --bin frankensim uq_command::model::material
cargo test -p fs-cli --test cooling_materials
```

The tests compare directional heat flow against an independent slab resistance
solution, rotate material and geometry together, compare implicit derivatives
with perturbed coupled solves, and use a manufactured nonlinear field with
known source and boundary heat rates. Actual-binary regressions also cover
material UQ and checkpoint replay. Source availability is not an assertion that
these tests have been executed on a particular checkout.

All these constitutive inputs are caller declarations, not material-database
receipts or physical uncertainty certificates. Reported temperatures remain
nominal discrete-model results. Sampling confidence concerns the declared
model and probability law, not missing material data or experimental validation.
