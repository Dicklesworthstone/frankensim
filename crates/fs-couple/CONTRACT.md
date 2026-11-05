

## Periodic driving of nonlinear physical-port feedback

`phs_orbit::solve_forced_feedback_harmonic_balance_*scaled_seeded` combines an
authored `PeriodicPortInput` with the original instantaneous
`AutonomousPortFeedback` law in the same full pHS ODE. Each residual uses
`(J-R) grad(H) + G [u_periodic(t) + u_feedback(x,G^T grad(H))]`. The externally
prescribed frequency and phase remain fixed, with no autonomous phase row or
unknown frequency. Original storage, native resistance, complete port order
and every explicit material-memory coordinate remain present.

The existing feedback tangent owner supplies `(J-R) H'' + G du_total/dx` to
the shared scaled/adaptive AFT solver. Its declared total feedback derivative
includes the physical flow/Hessian chain once; the periodic drive is
state-independent. Native complete storage tangents are discovered, with
explicit complete Hessian actions available for other laws. A refusing or
incomplete declared action never falls back to finite differences. Forcing
bandwidth is admitted before source preparation, including arbitrarily small
nonzero lines; zero padding does not consume resolution. All callbacks are
deterministic and read-only, with physical memory in explicit states.

`audit_forced_feedback_orbit` evaluates the original storage gradient once per
quadrature sample and retains separate signed external-drive and feedback
work through the original compensated physical input map. Negative feedback
work is net absorption by that authored law, separate from native R loss; it
is not automatically material heat. The combined original input supplies the
cycle power residual. Separately rounded large opposing component averages
need not sum bit-for-bit to that combined account.

`solve_forced_feedback_shooting_variational` and its `_analytic` variant expose
that same complete ODE and total derivative to the original discrete RK4
variational owner. Period is exactly the external source period, without an
autonomous phase row or period column. The complete state, monodromy and Floquet
spectrum retain material memory. Adapter and native solver resource costs share
the caller's byte/work caps; original forcing bandwidth is checked before
source preparation, and cancellation precedes preparation and is polled at
native stage/product boundaries. A declared tangent refusal never changes to
finite differences. No second integrator or constitutive law is introduced.

The bounded finite-order HB lane has no cooperative Cx cancellation. Neither
its sampled balance nor one discrete shooting map establishes an infinite
harmonic tail, physical-law validity, global stability or continuous-orbit
accuracy. Audit quadrature must resolve drive and retained state orders;
nonlinear work additionally needs sampling refinement. Focused G1/G3 tests
cover exact phased multiport quadratic/cubic harmonic conversion under unit
rescaling, retained relaxation memory, passive asymmetric drag versus native
time-domain shooting at different phases, the analytic discrete scalar Floquet
multiplier, nonlinear storage and signed power splits, and strict bandwidth,
source, tangent, resource and cancellation refusals.

