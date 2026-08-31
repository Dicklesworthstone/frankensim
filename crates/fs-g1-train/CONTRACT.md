# fs-g1-train contract

`fs-g1-train` owns the GQA transformer, its manual backward pass, PPO/GAE,
Muon/Adam optimizers, deterministic weight codec, and training receipts. It
does not own Unitree G1 rigid-body or contact dynamics.

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
  layouts directly; `fs-g1-train` must not duplicate them.

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

The real-owner seam gate is:

```text
cargo test --manifest-path crates/fs-cmaes-viz-wasm/Cargo.toml \
  --features g1-learned --test g1_learned_ppo
```

Neither gate by itself proves athletic performance, robustness, sample
efficiency, hardware transfer, or superiority to the phase-basis controller.
