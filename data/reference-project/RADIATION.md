# Surface radiation in the cooling project

`cooling-radiation.fsim` runs the canonical seven-stage cooling pipeline with
an additional gray exterior surface. It reuses `plate.stl` and the bulk
conductivity card `aa6061.fsmcdpk`. The additional `gray-surface.fsmcdpk`
contains the exact immutable emissivity card named by the project.

This is a manufactured reference case: 5 W in the existing tetrahedron,
convection at 10 W/(m² K) to 293.15 K, and emissivity 0.85 facing a black
isothermal reservoir at 293.15 K. Its material data are fixture declarations,
not experimental validation or recommended material properties.

From the repository root, using a built `frankensim` binary on `PATH`:

```bash
frankensim --json validate data/reference-project/cooling-radiation.fsim
frankensim --json import data/reference-project/cooling-radiation.fsim \
  data/reference-project/plate.stl radiation.db --unit m --max-hole-edges 0
frankensim --json solve data/reference-project/cooling-radiation.fsim radiation.db \
  --materials data/reference-project/aa6061.fsmcdpk \
  --materials data/reference-project/gray-surface.fsmcdpk
```

The solve result supplies a run ID. Export its retained report and package:

```bash
frankensim --json report <run-id> radiation.db
frankensim --json package <run-id> radiation.db
```

The conduction receipt separates `radiative_out_w` and `convective_out_w`.
When a surface uses `airflow-convection`, only the latter enters the air's
energy equation. Each patch retains its card identity, exact emissivity query
receipt, surface temperature, applied and nonlinear heat, and closure gates.
The HTML and JSON reports show the paired radiation-on/off temperature
comparison under model-form sensitivity. That difference is an Estimated
comparison of two physical models; it does not fill an unknown model-error
bound or certify compliance.

The admitted law is `epsilon sigma A (T_mean^4 - T_reservoir^4)` for each
named exterior convection patch. A positive secant coefficient acts on the
pointwise P1 Robin trace and is iterated until both temperature and heat
criteria close. Emissivity stays fixed at its declared card query temperature.
Base and uniform `ladder` fidelity support this law, including heterogeneous
conductivity and fixed contact resistance. `adaptive` fidelity explicitly
refuses it until the complete radiative tangent enters the goal calculation.
This first product port does not model enclosure exchange, occlusion or a
participating medium.

The existing test-fixture mechanism can regenerate exactly the two radiation
files without rewriting the original reference:

```bash
cargo test -p fs-cli --test solve -- \
  --ignored radiation_product::dump_radiation_reference_project_fixture
```
