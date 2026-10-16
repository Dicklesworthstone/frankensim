# Native contact-resistance controls

The native adjoint can include one resistance multiplier per declared matching
contact, using the existing imported geometry, interface cards and complete
thermal adjoint. Copy `contact-pair.fsim` to a separate project file and extend
its output clause:

```lisp
(outputs
  (qoi :name "temperature-max" :kind "scalar")
  (qoi :name "temperature-max-contact-adjoint" :kind "report"))
```

Import and solve that exact project with the same two mesh sources and card
packs as the ordinary contact-pair example. The derivative request changes
the project/run identity; it does not change the physical boundary or load.
Choose this report OR `temperature-max-adjoint`, not both. The original report
continues to omit contact controls and performs no additional contact work.

The conduction receipt's `nominal_adjoint` object identifies its `output` as
`temperature-max-contact-adjoint`. It retains the existing power, boundary,
material and admitted fan controls and adds `contact-resistance-multiplier`
rows. Each contact row binds its interface name, original card identity,
reference value 1, paired-face count and mapped/uniform status.

The coordinate is `R''_face -> s R''_face`, evaluated at `s=1`: `derivative`
is kelvin per dimensionless multiplier. For a uniform contact, dividing by its
nominal area-specific resistance converts this to the absolute resistance
partial. A mapped contact has no single absolute resistance to divide by;
the shared multiplier preserves its complete spatial resistance pattern.

The contraction uses both independent P1 jumps and each face's actual
resistance. It is not the product of patch-mean jumps, does not weld nodes,
and does not run another primal or adjoint. The supplied native dual already
includes the accepted material, air and radiation feedback. Prescribed
values remain in the full primal, with zero prescribed-node dual entries.

Results remain Estimated local derivatives of the selected hottest node, not
unique maximum derivatives at ties, material-uncertainty propagation, shape
or pressure/finish derivatives, temperature-dependent contact laws, or
physical-validation/error certificates. Native matching contact admission is
unchanged. The reusable matching-only API refuses nonmatching sets rather
than silently returning a delegated exact subpatch as the whole derivative.

`fs-conduction --test contact_controls` contains independent named/mapped
contact FEM comparisons and refusal checks. `fs-cli --test
native_contact_adjoint` exercises real imports, normalized interface-card
changes, physical re-solves and unchanged existing derivative rows. Test
presence is not a claim that a particular build has passed either target.
