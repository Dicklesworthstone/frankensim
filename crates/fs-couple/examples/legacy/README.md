# fs-couple: legacy examples

The examples in this folder run modules from `src/legacy/`: cavitation, rocket_blowdown, rod_motion,
superfluid, porous_heating and molecular_exchange. No app, script or CI job runs them.

Each keeps its name through an `[[example]]` entry in Cargo.toml, with the same required features
as before. For example, `cargo run -p fs-couple --features superfluid --example superfluid_piston`
runs as it did.
