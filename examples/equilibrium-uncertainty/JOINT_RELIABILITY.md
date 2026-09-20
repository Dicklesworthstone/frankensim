# Joint reliability of authored mechanical requirements

`equilibrium_uq --all-constraints` estimates the probability that **every
response requirement holds on the same parameter draw**, across all declared
load cases. It does not multiply marginal probabilities or assume independent
responses merely because the input variables are independent.

The existing `.fit` constraint records define displacement, signed left-side
spring force, compressive contact force, and positive penetration requirements.
Upper and lower inequalities are compared in their physical units. Each equality
needs an explicit `--equality-tolerance NAME VALUE`, also in its physical units:
`abs(response - bound) <= tolerance`. Zero is an explicitly requested exact
comparison. Normalization scales, solver tolerances and confidence parameters
never supply or widen that physical acceptance band.

## Example with an analytical joint event

From the repository root:

```sh
cargo run -p fs-cli --bin equilibrium_uq -- \
  examples/equilibrium-uncertainty/joint-reliability.model \
  examples/equilibrium-uncertainty/joint-reliability.fit \
  --method rqmc --replicates 4 --samples 512 --seed 73 \
  --all-constraints --independent --uniform-x force-N -0.5 0.5 \
  --equality-tolerance settled-band 0.00390625
```

The example has a 0.25 kg free mass restrained by a 64 N/m spring and a
128 N/m unilateral contact at a 1/64 m gap. Two independent stationary
experiments share the same declared uncertain physical force, uniform on
[0.5, 1.5] N. Their complete physical displacement is

- `q = F/64` below contact onset, `F <= 1 N`;
- `q = (F+2)/192` above contact onset.

The force ceiling admits `F <= 1.25 N` and the minimum displacement admits
`F >= 0.75 N`. Each of these two constraints has probability 0.75, but their
joint probability is **0.5**, not `0.75*0.75 = 0.5625`. The penetration ceiling
has the same upper-force event, the spring-reaction lower bound is inactive,
and the explicitly declared equality band has the same lower-force event
within this domain. All five requirements therefore jointly pass with
probability 0.5 under the analytical continuous-input model. Numerical sampling
and the original preload tolerances are not an exact probability certificate.

For Monte Carlo, replace `--method rqmc --replicates 4` with `--method mc`.
The existing paired `--require-probability P --confidence-alpha A` options can
stop MC when the joint indicator confidence sequence resolves that fixed
requirement. They remain forbidden for RQMC. Equality bands are fixed before
sampling and retained in the output. Searching among bands, thresholds, seeds
or models after looking at results requires separate multiplicity control.

## Output and boundaries

Joint mode emits `frankensim-equilibrium-reliability-v1`, with the joint
compliance estimate, probability of any failure, actual sample/case work,
source hashes, laws, coordinate scaling, and the complete requirements.
Per-constraint failure counts, empirical compliance rates and observed physical
ranges help identify which requirements fail; these are descriptive marginals,
**not simultaneous per-row confidence bounds**. MC confidence, when requested,
concerns the one joint event. RQMC standard errors come from complete independent
scramble means; they are not confidence sequences or iid pointwise errors.

All cases are solved once per draw using the existing primal owners. Force,
penetration, energy and equilibrium-residual gates remain mandatory. A failed
case or response arithmetic stops the study without partial statistics, dropped
samples or replacement draws. Contact onset can be observed without a derivative
margin: this grants no derivative at the kink. No new physics, optimizer,
sampler, independence assumption between responses, or physical validation is
introduced. Separate physical tolerances are requirements, not numerical-error
bounds. Evidence remains Estimated under the fixed reduced model.

An empty constraint family refuses rather than producing a vacuous pass.
Unknown, duplicate, negative or missing equality bands refuse. Joint mode
forbids `--case`, `--target` and `--limit-m`; the unchanged displacement-only
mode retains its v2 output and explicit unassessed-constraint count.

Focused native checks:

```sh
rch exec -- cargo test -p fs-couple --test equilibrium_forward_responses
rch exec -- cargo test -p fs-cli --bin equilibrium_uq --test equilibrium_joint_reliability --test equilibrium_uq
```

These tests were authored but not executed in the implementation environment,
which had no `rch`, `cargo` or `rustc`. Independent numerical reference checks
are not execution of Rust, its deterministic math, Philox or Owen scrambling.
