# render-gpu: the drum kit on an NVIDIA GPU, or on the CPU

`render-gpu` renders fs-vkplate kit pieces (kick, snare, toms, hi-hat, crash, ride). When a CUDA driver and an NVIDIA GPU are present, it steps the simulation on the GPU. Otherwise it runs the ordinary CPU solver, `fs_vkplate::kit::render`.

The GPU path is an optional accelerator, not a second solver:

- It steps the instrument that fs-vkplate exports (`Instrument::image`).
- It is checked against the CPU step (`--check`).
- The CPU step stays the reference.

Normal-contact damping reads current node positions rather than a cached intermediate
compression. Sliding friction solves the same anisotropic Coulomb disk problem as the
CPU, so its force opposes final slip. These changes require a fresh CPU/GPU comparison;
the historical measurements below do not qualify the corrected equations or their timbre.

The tool is its own cargo workspace. Nothing in the main workspace builds it or links its GPU bindings (the engine-kit pattern).

## Build and run

```sh
cargo build --release              # no CUDA toolkit needed to build; on macOS this is the CPU path only
render-gpu PIECE HITS.csv OUT.wav SECONDS [--backend auto|cpu|cuda] [--airy f64|f32] [--device N|all] [--preset-heads]
render-gpu PIECE HITS.csv OUT.wav SECONDS --chunks N [--preroll S]   # N overlapping chunks at once
render-gpu --batch PIECE SECONDS OUT_PREFIX HITS1.csv HITS2.csv ...  # one simulation per list, in lockstep
render-gpu --check PIECE SECONDS [--airy f64|f32]                    # CPU vs GPU from the same strokes
render-gpu --precompile sm_89                                        # fill the kernel cache, no GPU needed
```

`PIECE` and `HITS.csv` are those of fs-vkplate's `perform` example.

**Drums need `--preset-heads`.** Drum pieces (kick, snare, snare-off, tom-12,
tom-14) refuse by default. Their only card-free heads are the presets, whose T60
was fitted to recordings, so no default render path reads them (bead .370).
The sourced heads carry a head-radiation load, and the film loss curve's heads a
film law on their bending operator, that this GPU path does not model (the plan
refuses either), so render those natively with fs-vkplate's examples. `--preset-heads`
renders the fitted heads on purpose, and the notes name them as fitted. The
cymbals and the hi-hat are unchanged.

**Options.**

- `--device all` splits a batch, or a chunked render, over every GPU (one host thread each).
- `--exact-friction-order` solves frictions in the CPU's list order instead of a coloured sweep. It is slower; `--check` always uses it.
- `--grid-scale S` multiplies the hi-hat/crash/ride mesh spacing (default `1`, resource-bounded to `0.5..=8`). Values below one refine the same physical model for convergence studies; the library still chooses and checks the integration clock and can refuse a mesh within that range. Material, contact, support and microphone parameters stay fixed. This applies to CPU, CUDA, batches, chunks and `--check`. Drum-head pieces refuse a non-default scale. A coarser grid is not evidence of equivalent sound quality.

For an isolated time-refinement diagnostic, run `cargo run --release --example shell_energy_probe` from this workspace. It compares a lossless curved annulus at four time resolutions, with both initial velocity and a force pulse from rest. It reports energy drift and mass-weighted displacement error against the finest clock for linear and von Karman mechanics. This separates initialization and integration effects; it is not an instrument sound-quality test.

**Run-time requirements.** cudarc loads the CUDA driver and NVRTC at run time:

- A machine without them renders on the CPU, saying so (`--backend cuda` refuses instead).
- The kernel is compiled for the GPU found, straight to its machine code (CUBIN). An NVRTC newer than the driver therefore still works.
- The compiled kernel is cached in `~/.cache/fs-render-gpu/`.
- A compile takes about 5 s on an idle host. On a heavily loaded build host it can take 15 minutes. In that case, run `--precompile` elsewhere (only NVRTC is needed, no GPU) and copy the cache directory over.

**Without a GPU,** `--chunks` and `--batch` run one CPU thread per simulation. On macOS (no CUDA) that is the fast path for a long piece: the Strandbar hi-hat (17.4 s, 32 chunks) took 793 s on an M5 Ultra's 36 cores, on a busy host, against 2551 s rendered continuously on one core.

## How it works

`fs_vkplate::Instrument::image()` exports every operator, coefficient and piece of state one time step reads. It changes nothing; a test checks that the output is bit-identical with and without an export.

`plan.rs` flattens that export:

- one node numbering over all bodies;
- fixed-width column-major stencil tables;
- the nonlinear gradient as a gather, in the CPU scatter's order;
- contacts in node-disjoint groups (which solve exactly as one-at-a-time in list order);
- frictions in coloured groups;
- microphone taps grouped by delay, in balanced segments;
- the dense inverse of each cymbal's Airy operator, built from the CPU's own band factor.

`kernel.cu` is one persistent cooperative kernel. All blocks stay resident and meet at a grid barrier between the dependent phases of a step:

1. linear forces and link forces;
2. the Airy solve as a dense product (plates only);
3. cavity, port, radiation, links, drives and the nonlinear gradient;
4. the SAV update;
5. the Woodbury correction, strikes, contacts and frictions;
6. energy sums and microphone taps;
7. outputs.

That is five barriers per step for a drum and six for a cymbal. A large contact or friction group is spread over every multiprocessor, with a barrier in front of it; a small one runs on its simulation's block.

`--batch` and `--chunks` step many simulations of one piece in lockstep. They share every operator, and each pass over a cymbal's Airy inverse serves up to 16 of them.

**Numbers.**

- **Precision.** All state and arithmetic is FP64, except that `--airy f32` stores the Airy inverse in single precision.
- **Determinism.** Sums use fixed trees, so a run repeats bit for bit.
- **Difference from the CPU.** It rounds differently from the CPU's sequential sums; the measured differences are below.
- **Frictions.** By default they are solved as a coloured Gauss-Seidel sweep: 6 groups on a closed hi-hat, where the CPU's list order needs 603.
  - Each projection reduces free-mass kinetic energy for the state it sees; this is not a
    passivity certificate for the complete coupled oscillator discretization.
  - The sweeps can differ in waveform and spectrum; the comparisons below describe only
    the measured historical inputs and implementation.

**Chunks** follow the method validated for the Strandbar studio hi-hat (stitched within 0.06–0.3 dB per octave of a continuous render):

- each chunk starts with a pre-roll of the real strokes;
- each cut sits 10 ms before a stroke;
- the windows meet in a 4 ms raised-cosine crossfade.

A hi-hat adds a stick contact per stroke, so chunks are padded to the same stroke count with strokes after the end. To the CPU model, such a stroke is no different from a later real one.

## Measured (RTX 4090 on trj; CPU per step from one M4 core)

### Agreement with the CPU

`--check`, 2 ms windows, compared with the signal RMS. The CPU friction order and f64 Airy storage are used throughout.
The command now exits nonzero if any window exceeds the existing `1e-3` relative-RMS threshold.
Missing microphones/samples, nonfinite trajectories, and an entirely silent CPU reference also refuse qualification.
A successful process exit therefore requires a complete numerical comparison; it does not certify acoustic quality or real-time throughput.
The comparison constructs the instrument once, then copies its complete initial state into two independent trajectories. This avoids repeating mesh factorization and the builders' stability checks. Construction time is reported separately from CPU stepping and GPU setup/stepping; use the complete process duration when assessing job latency.

| piece | largest window difference | octave bands |
|---|---|---|
| kick, 0.3 s | 2.3e-9 | 0.00 dB |
| tom-12, 0.3 s | 1.1e-9 | 0.00 dB |
| closed hi-hat, 0.65 s, two strokes | 3.3e-8 (never above 1e-3) | 0.00 dB |
| snare | 7e-8 until the second stroke, about 300 ms in | within 0.22 dB |
| crash, 0.15 s, two strokes (GPU render against an M4 CPU render, 48 kHz float WAVs) | 1.2e-6 | 0.00 dB |
| ride, likewise | 7.5e-6 | 0.00 dB |

The whole closed hi-hat step agrees with the CPU to 3e-8 over the run: von Kármán shells, 1204 rim contacts with friction, clutch links and the pedal drive. The snare's buzzing wires are chaotic, so after the second stroke the two renders differ sample by sample but not in their spectrum.

### Default settings, against a CPU render on an M4

Closed hi-hat, two strokes:

- **Coloured frictions:** 0.5% relative RMS difference from 12 ms on; every octave band within 0.01 dB, with the same peak.
- **`--airy f32` on top:** no further change (5.11e-3 against 5.09e-3).

### Studio renders

**Strandbar kick** (17.368 s, 32 strokes):

- Continuously on one GPU: within 7e-9 relative RMS of the CPU studio stem, and 4.8e-7 Pa at most on an 8.86 Pa peak. 27.6 s on the GPU, against 304.6 s on the CPU of the same busy host.
- As 16 lockstep chunks: 0.00 dB in every octave band.

**Strandbar hi-hat** (17.368 s, 128 strokes; the CPU studio render took 2188.9 s on one core):

- Run as 32 chunks, 0.5 s pre-roll, `--airy f32`, on both GPUs.
- **287 s from start to end, 7.6x faster;** the stepping itself took 208 s on each GPU in parallel.
- Against the continuous CPU stem: within 0.24 dB in every octave band (+0.18 to +0.24 dB from 4 to 16 kHz); RMS +0.05 dB; peak 2.555 against 2.567 Pa.

### Speed

µs per step, and per simulation-step for a batch:

| piece | CPU (one M4 core) | GPU, one simulation | GPU, batch |
|---|---|---|---|
| kick (2 × 7589 nodes, 48 kHz) | 172 | 33-38 (4.5-5x) | 12.7 per simulation at 32 (13.5x) |
| closed hi-hat (2 × 3292 nodes, 288 kHz, 1204 rim contacts + frictions), `--airy f32` | 325-356 | 156 (2.1-2.3x) | 30.6 per simulation at 16 (11x) |
| closed hi-hat, `--airy f64` | 325-356 | 282 | 78 per simulation at 4 |

With both GPUs, a batch or chunked render runs at twice the rate in that last column.

### Where a step's time goes

Cycle counters per phase, with `RENDER_GPU_PROF=1`:

- **A drum step** is about ten dependent memory round trips plus five grid barriers of about 1 µs each.
- **A hi-hat step** is dominated by:
  - its 1204 implicit rim contacts and frictions: nine sequential groups of FP64 Newton solves, about 100 µs;
  - the Airy product: 73 MB in f32, just over the 72 MB L2, about 35-40 µs.

## Limits, and where more speed is

**FP64 on a GeForce card.** An RTX 4090 runs double precision at 1/64 of its single-precision rate: two FP64 units per multiprocessor.

- One multiprocessor does FP64 work about eight times slower than one CPU core.
- So one simulation on the GPU is latency-bound. Each step's phases are a handful of dependent memory trips plus about 1 µs per grid barrier.
- The GPU wins large through batches and chunks, not through a single simulation.
- A 64-core AVX2 CPU running 64 renders at once has more peak FP64 throughput than one 4090 (about 3.6 against 1.3 TFLOPS). The GPUs pay off as extra capacity, and for one long render that cannot be split.

**The next large step is single precision** where it is safe:

- a difference-form leapfrog (store w_n - w_{n-1}, and keep w as a compensated pair);
- the Airy right-hand side expanded bilinearly about the rest shape, so there is no cancellation between w and z0;
- FP32 Newton iterations polished in FP64 for the contacts.

That turns most of the arithmetic into FP32, 64 times the throughput, and it is also what an Apple GPU needs, since Metal has no FP64 at all (bead frankensim-7z7ew).

**The cymbals' Airy solve.** On the CPU this is a banded Cholesky, sequential in its rows and 70% of a hi-hat step; a GPU cannot run it in parallel. Here it is a dense inverse:

- the matrix is read once per step;
- f32 storage halves the bytes;
- a batch shares the reads.

A solver with far fewer operations per step would speed up both paths (bead frankensim-b1y5x). Candidates are a hierarchical low-rank inverse, or a fast-sine-transform solve with a capacitance correction for the disc boundary.
