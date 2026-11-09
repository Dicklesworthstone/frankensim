


## Implicit laminar enthalpy with external RCR storage

With `laminar-enthalpy`, `LaminarEnthalpy::advance_windkessel` resolves the
original nonlinear BDM momentum/rheology and RCR/Windkessel ports, including
optional series inertance, inside each native backward-Euler caloric residual.
The original material is queried at the actual trial temperature and shear;
its accepted source snapshot remains owned by the returned flow. Original
viscosity, convection, pressure response, final momentum equations and both
unreduced circuit laws retain their native admission and residual gates.
Entirely circuit-connected components require no artificial traction datum.

Every thermal/Newton probe starts from the SAME accepted velocity and compliance
pressures. Previous port Q is integrated from that native velocity, never
supplied independently or advanced by a provisional thermal trial. Compliance
coordinates use the same reservoir references and physical clock; geometry,
port ordering, density and C/L/R remain fixed, and reservoir increments are held
through the interval. BDF1 is used consistently for momentum, external storage
and caloric transport. This entry adds no adaptive circuit window, moving mesh,
variable density or Boussinesq/load combination.

The heat callback borrows the complete native circuit/flow result. Exact
sign-split BDM rates drive the original enthalpy donor selection, with explicit
external inlet caloric states required on every incoming branch, including
reversal. Bulk strain and explicitly routed relative-slip heat are deposited
once. External resistor losses, compliance/inertance storage, reservoir/port
work and numerical BDF defects remain separate; no thermal receiver is inferred
for an external load. Original additional heat links and sources remain usable.

Acceptance moves the original resolved field into the standard coupled
candidate and native pressure/Q/storage/power receipts into its typed exchange
payload without cloning or reconstructing the field. The receipt exposes
accepted compliance coordinates for continuation with that same candidate's
velocity and thermal clock. Native energies/powers remain J/m and W/m; declared
depth and actual interval convert them to joules. Port power covers added
pressure increments relative to original traction, not an inferred full-system
energy or passivity certificate. Source, stream, numerical-budget, cancellation
or final thermal refusal publishes none of the three physical histories.

Focused native regression definitions cover independently eliminated two-port
RCR/inertance momentum and storage equations, signed reversal with two-cell
caloric transport, external energy versus BDF damping, original Arrhenius/
Carreau viscosity with a manufactured finite-slip thermal endpoint, and late
source/stream/budget/cancellation rollback. All four changed Rust files passed
the pinned parser, formatter and `rustfmt --check`. The full native FEM/thermal
dependency build and these integration test executions remain unverified in
the current partial restored environment; no numerical pass is inferred from
formatting or static review.
