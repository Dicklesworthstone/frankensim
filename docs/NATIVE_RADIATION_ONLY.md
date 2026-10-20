# Radiation-only native cooling

`data/reference-project/cooling-radiation-only.fsim` uses the retained tetrahedron,
conductivity card and gray-surface fixture card without an air network or a
convection coefficient. The inherited region name `air` is just an assignment
identity; its material is the same solid aluminium fixture as the reference.

```bash
frankensim --json validate data/reference-project/cooling-radiation-only.fsim
frankensim --json import data/reference-project/cooling-radiation-only.fsim \
  data/reference-project/plate.stl radiation-only.db --unit m --max-hole-edges 0
frankensim --json run data/reference-project/cooling-radiation-only.fsim \
  radiation-only.db --materials data/reference-project/aa6061.fsmcdpk \
  --materials data/reference-project/gray-surface.fsmcdpk
```

Use the returned run identity with the ordinary `report` and `package` verbs.
The example retains the original version-5 wire representation and passes
through the existing versioned migration path; zero-flux syntax is not new.

## Declare no nonradiative exchange

A radiation-only surface has an explicit physical boundary:

```lisp
(heat-flux :target "air" :outward-flux 0W/m^2)
```

The separately declared radiation surface must match that target and supply its
immutable surface card, optional exact claim pin, query temperature, independent
black-reservoir temperature and bounded outer iteration controls. Merely adding
an empty or unmatched radiation block does not anchor a steady solve. Neither
nonzero prescribed flux nor fixed temperature is admitted on the same radiating
trace. Existing convection/radiation combinations keep their semantics.

The actual positive radiative secant supplies the Robin operator. No artificial
convection coefficient or ambient-fluid temperature is introduced. Radiative
heat remains separate from air heat; a wholly radiation-only body reports zero
convective watts. The steady radiation-off comparison is `no-data` when turning
radiation off removes every temperature anchor: there is no unique equilibrium
to subtract. A declared finite-time capacity model can still have an ordinary
radiation-off endpoint, and that comparison continues through the existing solver.

## Supported path and limits

Native base steady solves and fixed-capacity transient solves use the existing
card-backed mean-temperature radiation kernel. The three native regressions in
`solve/radiation_product/vacuum_product.rs` cover steady global balance, invalid
boundary declarations and both cooling/heating transient endpoints. They also
check retained reports and the distinction between steady and finite-time
radiation-off comparisons.

The lower-level ambient-radiation adjoint now supports explicit zero-flux
radiating traces, including reservoir, emissivity and nodal-load derivatives.
The native adaptive and native adjoint selectors are not extended by this
integration; do not infer those capabilities from base-solve admission.

This remains an Estimated area-mean gray-surface model with view factor one to a
black isothermal reservoir. It is not enclosure radiosity, a pointwise fourth-
power integral, temperature-dependent emissivity or experimental validation.
The fixture cards do not authenticate a real coating or qualify a vacuum design.
