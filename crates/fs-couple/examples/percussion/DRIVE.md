# Physical stick-force performances

`--stick-force-file player.csv` drives the **existing stick**, through its
geometry-derived effective mass and contact point, throughout a render. It
works with every existing splash/drum/snare command and numerical image,
including the admitted cavity and stretching combinations. The same stepped
mechanics feeds CSV observations and the existing BEM pressure renderer.

The file has no header. Each row is `time_s,force_n` in seconds and newtons;
blank lines and `#` comments are allowed. Times are finite, nonnegative and
strictly increasing. At least two rows are required; the first and last force
must be zero. Force is zero outside these endpoints. Positive force pushes
toward the head; negative force pulls the stick away. For example:

```csv
# A delayed push from rest; this is player input, not a contact-force fit.
0,0
0.01,0
0.012,2
0.014,0
```

```bash
cargo run --release -p fs-couple --example percussion -- \
  drum-modal 20000 --strike-speed-m-s 0 --stick-force-file player.csv > motion.csv

cargo run --release -p fs-couple --example percussion -- \
  drum-modal-mic 48000 20 --strike-speed-m-s 0 \
  --stick-force-file player.csv > pressure.wav
```

The initial `--strike-speed-m-s` remains an independent physical input. It is
**not automatically set to zero**: the legacy default is still 0.8 m/s. Supply
zero for a performance that starts with a stationary stick. The strike location,
materials, geometry, head basis, snare contacts, cavity and acoustic boundary
are unchanged. Further pushes and lifts act on the moving stick; they do not
teleport it, prescribe an impact velocity, restart a one-shot sound or clear
head/air/felt/snare history. Impact timing and rebound follow from mechanics.
This is a fixed-axis, single-stick player port, not a two-hand motion capture or
feedback controller. A row does not guarantee a strike at its timestamp.

Linear interpolation is integrated over each mechanical interval, including
knots between ticks. Its mean force is supplied as that interval's external
load, as `F/sqrt(m)` in the mass-normalized stick coordinate. This preserves
input impulse under time subdivision; it does not claim substep trajectory or
sound convergence for unresolved rapid force changes. The existing owner still
checks external work, energy, contact and state limits. Refusal does not advance
the playing clock. Drive preparation allocates once; it does not make the
reference Newton solver allocation-free or certify real-time performance.

The complete file must fit the render duration; it is never silently truncated.
Limits are 4 MiB, 65,536 knots and the existing generalized-force envelope.
Driven mechanics CSV appends `player_work_j`, the accepted step's external work,
not a cumulative total. Unforced CSV columns and all PCM scaling remain unchanged.
The WAV observer still has its documented causal filtering and propagation delay.

G0/G3 regressions cover signed pulse integration, subdivision, malformed input,
refusal/retry, physical mass/work scaling and a delayed drive through the actual
reference and prepared drum contact assemblies. Run the example test target via
the repository's DSR/RCH development lane; native test execution is required.
