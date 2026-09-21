# Snare wires and distributed air in one mechanical solve

`--cavity-modes` now works with `snare`, `snare-off`, and `drum-modal`, including
all their existing `-wav` and `-mic` variants. It extends the original nonlinear
path in `CAVITY.md`: the old restriction to `drum` and `drum-stretch` no longer
applies to sealed prepared cavities. Commands without the option are unchanged.

```sh
# Mechanical onset: 20 distinct strands, two heads, spatial enclosed air.
(set -C; cargo run --release -p fs-couple --example percussion -- \
  snare 4096 --cavity-modes --strike-speed-m-s 4 \
  --strike-position-m 0.06 0.01 > snare-cavity.csv)

# The same physical system observed through the existing exterior BEM path.
(set -C; cargo run --release -p fs-couple --example percussion -- \
  snare-mic 4800 2 0.08 0.05 0.35 --cavity-modes \
  --drum-spec crates/fs-couple/examples/percussion/estimated_drum.fsd \
  --strike-position-m 0.06 0.01 > snare-cavity.wav)
```

Do not add `--prepared-nonlinear` to these commands: they already use the
prepared modal/contact image. `drum-stretch --cavity-modes --prepared-nonlinear`
remains the separate nonlinear-head image. Stretching heads combined with the
full snare bank, necks on the prepared snare, and vented exterior audio still
refuse. No nonlinear potential or declared loss is silently dropped.

## Actual coupling, not a pressure overlay

The original film FEM, installed tensions, damping, striker, wire geometries,
and contact laws are unchanged. The reference snare retains all 20 strands,
160 wire coordinates and 240 distributed wire/head contact points, plus the
stick contact. Each strand still moves independently. The cavity basis and head
surface quadrature are the same ones used by the nonlinear drum. A shared
construction supplies the uniform compression mode and nonuniform pressure
modes; it does not install a second copy of the compact gas spring.

Each acoustic spring acts through its complete signed column: both heads plus
its air inertia. That column cannot be replaced by independent pairwise springs
without changing the Hamiltonian's cross terms. The new reusable
`CavityCoupling::build_linear` compiles a direct sum of the unchanged diagonal
coordinates into the existing prepared connection/contact solver. No eigenbasis
is changed, no wire is homogenized, and no new time integrator is introduced.
All cavity and contact reactions are simultaneous, not previous-sample forcing.

Striker, heads and wires keep their original state addresses; acoustic inertia
is appended last. Wires and striker have zero direct enclosed-volume coupling.
Wires affect air through the resonant head. External audio still observes actual
head motion through the original closed exterior boundary and BEM transfers.
Neither interior pressure nor an extra wire/gas oscillator is summed into the
microphone signal. Direct wire radiation, radiation backreaction and room
scattering remain outside this model.

## Budgets and evidence

The prepared snare keeps its existing limits of 256 total coordinates, 512
contact points and 50 million contact setup terms. Air coordinates count toward
that same total; overflow refuses rather than truncating the wire bank. The
nonlinear reference's 64-mode limit is unchanged. The prepared cavity requires
zero acoustic momentum drag because its free-coordinate owner does not provide
that damping law. This example already declared zero acoustic drag; no loss
coefficient was changed to make the composition work.

The five library tests cover uniform-volume equivalence, independent coupled
eigen-dynamics under time refinement, basis-rescaling invariance, many-body
contact and cancellation/retry, and unchanged loss/budget refusals. Five example
tests cover real 20-strand assembly, changed head motion and spatial pressure,
prepared/reference onset comparison, supplied-geometry audio construction, and
command/physics restrictions. The audio-construction test does not execute BEM
fitting or certify a waveform. Focused native commands are:

```sh
cargo test --release -p fs-couple --lib render::plate::impact::cavity::prepared::tests
cargo test --release -p fs-couple --example percussion -- --test-threads=1
```

The conservative equations were also evaluated independently in Python: halving
the test step reduced the trajectory error by about fourfold. That is not a
native Rust result. Native builds, full-band/modal convergence, specimen
calibration and real-time deadlines require their own execution and evidence.
