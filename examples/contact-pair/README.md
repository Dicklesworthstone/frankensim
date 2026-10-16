# Worked example: two bodies and a declared contact joint

Two aluminium blocks meet on a 1 m² face. The hot block dissipates 5 W and
its outer faces carry zero flux. The cold block's outer faces are held at
293.15 K. Between them sits a declared, card-backed dry contact
(`cold-hot.fsintpk`, R'' = 0.1 m² K/W at the declared 1 MPa). It is the
product-path form of the contact oracles in `crates/fs-cli/tests/solve.rs`.
Those oracles generate this directory (`dump_contact_pair_example`), and a
G0 test keeps the tracked bytes equal to the generator's output.

```bash
cargo run -p fs-cli --bin frankensim -- --json import \
  examples/contact-pair/contact-pair.fsim \
  examples/contact-pair/cold-body.stl examples/contact-pair/hot-body.stl \
  "${WORK}/pair.db" --unit m --max-hole-edges 0
cargo run -p fs-cli --bin frankensim -- --json solve \
  examples/contact-pair/contact-pair.fsim "${WORK}/pair.db" \
  --materials data/reference-project/aa6061.fsmcdpk \
  --interfaces examples/contact-pair/cold-hot.fsintpk
```

`import` takes one source per geometry row, in declaration order. Swapping
the two STLs refuses on the pinned source hash, and passing one refuses
`cli-import-source-count`.

## What conservation fixes

All 5 W have one exit, the joint, so the contact heat is exactly −5 W (hot to
cold is negative in the receipt's a→b orientation). The mean temperature jump
is exactly −0.5 K = −5 W × 0.1 m² K/W / 1 m², whatever the mesh. MEASURED
2026-09-29: `heat_rate_a_to_b_w` −5.0000000037, `mean_jump_k` −0.5000000004,
energy closure 7.5e-10, seven stages. The G1 test
`g1_the_contact_pair_imports_two_sources_and_conducts_through_the_declared_joint`
checks these numbers. The interface card is fixture data (its provenance says
so), not a measured material pair.
