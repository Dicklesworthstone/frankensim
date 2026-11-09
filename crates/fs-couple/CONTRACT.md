


### Adaptive laminar enthalpy and external hydraulic continuation

`LaminarEnthalpy::advance_windkessel_adaptive` reuses the original coupled
BDF1 coarse/two-fine controller. Coarse and first-half solves borrow the same
accepted flow/compliance history; the second half borrows the actual first-half
pair. Rejected comparisons advance no external state. The existing ordinary,
Boussinesq and driven entries retain their original semantics and leave the
new optional hydraulic observations empty. No second time controller or
extrapolated physical field is introduced.

An explicit hydraulic norm joins the original caloric, velocity and work
criteria. It checks each port's compliance and added pressure, actual integrated
native Q, separate C/L storage, and signed plus resolved-magnitude resistor,
reservoir and port work. Signed and magnitude Q*dt*depth throughput have their
own cubic-metre tolerance, so equal endpoint flow cannot conceal different
transported volume. Pressure scaling removes the held reservoir offset and
uses the start compliance/Ohmic state; it invents no initial acceleration.
Storage changes and errors use difference of squares before dimensional
multiplication. All spans and local corrections come from the original
first-order controller, including exact target clipping and one-ULP cap tails.

Only accepted fine receipts contribute compensated per-port work and signed/
absolute transported volume totals. Initial Q is independently obtained from
the original accepted BDM edge transport; initial/final C*x²/2 and L*Q²/2
changes use the paired physical states and declared depth. Each port and their
aggregate must satisfy port work minus reservoir work, both resistor losses,
both storage changes and the separate native BDF defect within the explicit
window allowance. This external-circuit identity does not add heat, replace
the independent native caloric balance, or certify whole-system passivity.

All coarse/rejected work consumes the shared native-step and thermal-endpoint
ceilings. Source, material, momentum, stream and thermal refusals remain
terminal. Cancellation or a failed final circuit account rejects the complete
window, including any internally accepted prefix. The final step's retained
compliance coordinates continue with the returned original velocity and clock.
Circuit parameters/references and physical boundary data remain held; this
entry adds no scheduled circuit events or Boussinesq/load combination.

Native regression definitions cover symmetric compliance discharge with zero
fluid/thermal change, hydraulic-only refinement and an independent accepted BE
recurrence; signed equilibrium all-port plug transport and exact work/storage;
nonzero resistor losses excluded from heat; decimal target completion; and
late source, work-budget and cancellation rollback after an accepted prefix.
All seven changed Rust files passed the pinned parser, formatter and
`rustfmt --check`. Bounded independent source review found no paired-history
or accounting blocker. Full native FEM/thermal compilation and these native
integration test executions remain unverified in the current partial restored
environment; no runtime pass is inferred from those checks.
