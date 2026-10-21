# Binary endpoints with response and physical-stress requirements

A feasible gray design is not a feasible thresholded design. The binary path
re-solves every independent force/motion experiment at the actual two-phase
stiffness and checks its original response bounds AND physical material stress.

```sh
cargo run -p fs-topopt --features cutfem-marquee --release \
  --example response_stress_sdf3 -- 80 500000 0.0012 0.008 0.1 \
  --binary 0.25,0.5,0.75 0.001
```

Use the repository's DSR/RCH lane where available. The last argument is explicit
void stiffness, not a density floor inferred from the gray problem. Thresholds
are an ordered comma-separated schedule of at most 64 values in (0,1]; duplicates
are permitted and solved independently. The original 500000-iteration budget
covers optimization AND all endpoint trials; no work allowance is renewed.
Without `--binary`, the existing continuous command behavior is unchanged.

## Physical endpoint, not gray stress

`evaluate_binary_response_stress_limits` accepts `BinaryResponseStressOptions3`
from `sdf3::response::limits::binary::strength`. The existing endpoint owner
filters/projects the source, thresholds once, respects protected solid/void
labels and constructs stiffness scales of exactly 1 or the supplied positive
void stiffness. It does not refilter the binary mask. All original loads,
support patches and retained quadrature remain, including those on weak cells.
Prescribed-motion liftings and reaction functionals use the endpoint material.

Stress uses the fresh endpoint displacements and actual `scale*C_ref*eps(u)`
from the original bulk-stress integration. No gray displacement, qp-relaxed
stress, adjoint or derivative is reused. Every declared case contributes; weak
cells are NOT suppressed. Per-case solid and void maxima are reported separately.
There is no public cached-equilibrium entry point. Incoming study scales are
restored on every returned result, including cancellation after the response
solve but before completion of stress sampling. Callback panics are not caught.

The physical-stress aggregate is the unweighted p-norm of all retained von Mises
samples, with fixed p in [2,64]. Streaming rescaling avoids overflowing merely
because individual powers are large; no averaging can hide a failed sample.
Only one case's stress-point vector is retained at a time. The complete family
shares `max_points`; every candidate also uses the existing cumulative solver
control. Stress sampling adds no preparation, primal, adjoint or geometry solve.

This physical-only aggregate is deliberately distinct from the differentiable
gray problem's combined relaxed/physical stress measure. The numerical stress
cap and actual maximum are checked separately, alongside every original signed
response residual and the same unnormalized joint log-sum-exp upper aggregate.
The resulting conservative finite-set constraint can exclude a tight feasible
window. Its value depends on the retained sample count; it is not a continuum
maximum, error bound or strength certificate.

## Explicit threshold selection

`scan_binary_response_stress_limits` takes a schedule and
`BinaryStressScanOptions3`. Its `stress.endpoint.threshold` is replaced by each
schedule entry; the other material/response/stress policies stay fixed. It
projects the source once and prepares/re-solves each candidate. Reaction windows
and physical stress need not be monotone in threshold, so the scanner does not
bisect or interpolate feasibility. The least-volume completed joint-feasible
candidate under the explicit numerical volume cap wins; equal-volume ties keep
the first. This is not a global binary topology optimizer.

Bad schedules or physical policies refuse before law callbacks or solves.
Later errors are retained in `outcome` together with completed summaries, their
best feasible fields and cumulative work. `Ok` with no incumbent means that the
schedule found no feasible endpoint. Cancellation before retaining a candidate
omits it; cancellation afterward keeps its full evidence. Earlier incumbents
survive failure of a later trial. Only the incumbent/current field families
are retained, not displacement vectors for the entire threshold schedule.

The example prints candidate response/stress acceptance separately, then every
original value/bound/residual and the selected stiffness mask. It prints complete
prefix evidence before returning a late scan error. No feasible binary endpoint
or an unfinished scan returns nonzero; the previously selected continuous design
is not overwritten. The example imposes no additional binary-volume restriction
beyond a fraction of 1; the library exposes that cap explicitly.

Positive void stiffness, unchanged loads on weak cells and retained geometry
are declared modeling assumptions, not a manufactured solid. There is no
zero-void limit, surface export, connectivity/minimum-feature proof, follower
load, binary threshold derivative or continuum safety claim.

Seven endpoint regressions, four scanner regressions and two command argument
tests are included. Native execution is still unverified: DSR, its documented
absolute fallback, RCH and Cargo are absent in the editing environment. Source
identity and lexical checks do not replace native compilation or physics tests.
