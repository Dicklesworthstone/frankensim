# fs-g1-train contract

`fs-g1-train` owns the GQA transformer, its manual backward pass, PPO/GAE,
Muon/Adam optimizers, deterministic weight codec, and training receipts. It
does not own Unitree G1 rigid-body or contact dynamics.

## Purpose and layer

Layer L6 standalone policy-training support for the G1 browser product; it
owns the transformer, PPO/GAE, optimizers, codec, and training receipts.

## Public types and semantics

The public surface re-exports `GaitTransformer`, `G1Env`, `PpoConfig`,
`RunningNorm`, `Trajectory`, `MuonParam`, and `AdamParam`; the stand-in
environment and checkpoint rules below define their product boundary.

## Invariants

Published stand-in checkpoints require a post-update deterministic greedy
evaluation with nonzero policy-head norm, forward distance, and actuator work.

## Error model

The finalizer refuses legacy checkpoints lacking the required identity or any
of those three publication facts; it does not relabel them.

## Determinism class

Training and export use explicit seeds and a length-prefixed little-endian f32
codec; the documented focused gates check this boundary.

## Cancellation behavior

No cancellation behavior is documented for this standalone training contract.

## Unsafe boundary

No unsafe boundary is claimed by this contract.

## Feature flags

None are declared in the current manifest; `g1-learned` belongs to the sibling
`fs-cmaes-viz-wasm` integration feature.

## Conformance tests

The standalone and owner-seam focused commands listed below are the current
conformance gates. The owner seam's declared training environment is the
source-bound `G1Task::Walking` / `G1Challenge::Flat` owner rollout at 480 Hz
for 0.5 s, initialized from the disclosed walking curriculum. It records PPO
updates and the deterministic zero-head/base-policy counterfactual. It is a
training seam only: no held-out validation split, seed-disjoint terrain/push
campaign, or publishable locomotion artifact is implemented by this gate.

## No-claim boundaries

A stand-in checkpoint is not learned humanoid locomotion, robust performance,
hardware transfer, or superiority evidence without the additional receipts
already specified below.

## Environment boundaries

- `standin_env::StandinEnv` is the disclosed `action-causal-standin-v2`
  explanatory environment shared with `cmaes_explainer`. Its 42 inputs are
  joint positions, joint velocities, base roll/pitch/yaw, base angular
  velocity, base linear acceleration, gait phase sine/cosine, and target
  speed. Its 15 bounded actions are reduced lower-body residual commands.
- Stand-in forward motion is a kinematic proxy derived from lower-body joint
  velocity and bilateral hip opposition. Zero or non-finite action produces
  no commanded progress. The stand-in is not FrankenSim SE(3), contact,
  terrain, or push authority.
- Real 29-DoF owner training is composed by `fs-cmaes-viz-wasm` behind its
  `g1-learned` feature. That adapter consumes the owner observation and action
  layouts directly; `fs-g1-train` must not duplicate them. The focused seam
  uses the owner's 42-element observation vector and 15 residual-action
  coordinates, as defined by the owner, rather than the stand-in units.

## Artifact publication

The stand-in ablation exporter may publish a checkpoint only after a PPO
update and only when deterministic greedy evaluation under
`action-causal-standin-v2` proves all of the following:

1. the policy head has nonzero norm;
2. the policy causes nonzero forward distance;
3. the policy performs nonzero actuator work.

The finalizer independently rechecks the contract identity and those three
conditions. Legacy checkpoints without the identity, including the historical
iteration-zero all-zero-head artifact, are refused rather than relabeled.

A stand-in checkpoint is an optimizer/export-parity artifact, not learned
humanoid locomotion. A downstream UI may use a learned-locomotion label only
for an owner-coupled artifact with held-out terrain and push receipts, action
counterfactuals, source and environment identities, a content digest, sample
and wall-clock counts, and deterministic import/forward parity.

## Determinism and gates

Training and export use explicit seeds and length-prefixed little-endian f32
arrays. The focused standalone gate is:

```text
cargo test --manifest-path crates/fs-g1-train/Cargo.toml
```

The real-owner Walking/Flat training-seam gate is:

```text
cargo test --manifest-path crates/fs-cmaes-viz-wasm/Cargo.toml \
  --features g1-learned --test g1_learned_ppo
```

Neither gate by itself proves learned locomotion, athletic performance,
held-out terrain/push robustness, sample efficiency, hardware transfer, or
superiority to the phase-basis controller. It also does not make the walking
curriculum a from-scratch baseline: the curriculum remains an explicitly
disclosed owner-policy initialization.
