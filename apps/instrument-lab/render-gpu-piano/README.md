# render-gpu-piano: the grand piano on a GPU

`render-gpu-piano` is the `grand_piano` renderer from fs-couple's examples with one addition: stretches of the performance in which no hammer can touch a string are stepped on a GPU.

- On Apple silicon it uses Metal.
- On machines with an NVIDIA driver it uses CUDA.
- Otherwise it uses a CPU emulator of the same phases.

It compiles the example's own source files (`#[path]`, nothing copied) and takes the same command line, so it writes the same outputs. The CPU step stays the reference: the GPU path is checked against it step by step (`FS_PIANO_GPU_CHECK`), and every stretch must close the instrument's energy ledger.

## Build and run

```sh
cargo build --release      # Metal on macOS; CUDA loaded at run time elsewhere (no toolkit needed to build)
render-gpu-piano --preset steinway-d --performance score.csv --duration 17.4 --fast-dampers \
  --board-band-hz 2400 --board-reduction 128,64,200,400,800,1600,2300 --render piano.wav
render-gpu-piano --precompile sm_89   # fill the CUDA kernel cache (NVRTC only, no GPU needed)
```

The tool is its own cargo workspace. Nothing in the main workspace builds it or links its GPU bindings (the render-gpu and engine-piano pattern).

**Environment.**

| Variable | Effect |
|---|---|
| `FS_PIANO_GPU` | `auto` (default; falls back to `cores` when no GPU backend is usable), `metal`, `cuda`, `cores` (the phase emulator on the CPU's performance cores), `cpu` (the same on one core) or `off` (plain CPU render) |
| `FS_PIANO_CPU_THREADS` | Threads for `cores` (default: performance cores on macOS, available parallelism elsewhere) |
| `FS_PIANO_GPU_CHECK=N` | Limits stretches to N steps and compares every step with `Instrument::step_linear_reference` from the same state. Reports the worst relative state, board-trace, contact-position and loss differences. |
| `FS_PIANO_GPU_REAL` | CUDA only: `double` (default) or `df64` |
| `FS_PIANO_GPU_ADD` | Metal: `sloppy` (default) or `ieee` double-float addition |
| `FS_PIANO_GPU_BENCH=0,1,...` | Times stretches from rest for each `P_SKIP` mask in kernel.h (0 is the real kernel), then exits |
| `FS_PIANO_GPU_BLOCKS`, `FS_PIANO_GPU_WARPS`, `FS_PIANO_GPU_DEVICE` | Overrides for the block budget, warps per block (default 8) and CUDA device |

**Configurations it models.**

- The contact-free step of the linear string and board bank, with or without string loss and board damping.
- Spatial viscous dampers (`--dampers estimated`, pad files, the preset), including the pedal and lifted keys.
- One substep per sample.

Felt dampers, point-drag dampers, exterior radiation, the cavity, string stretching and hammer footprints are not modeled; the tool says so and renders on the CPU. The CPU also renders without a physical microphone (`--diagnostic-volume`) and when the instrument is observed between steps (`--bridge-trace-csv`). With the GPU path the CPU's own steps carry the dampers on reduced coordinates, as `--fast-dampers` does (the same flows; rounding differs at about 1e-16). If the GPU fails during a render, the remaining stretches run on the CPU phase emulator.

## How it works

**Split.** The hammers, the felt and every contact stay on the CPU (`grand_piano::audio::Accelerator`).

- The planner (`Instrument::contact_free_steps`) flies the active hammers forward against each string's amplitude bound. The stretch it allows runs on the GPU, up to the next control event and at most 2048 steps.
- The GPU reports every step's contact positions. The CPU then replays the hammers' free flight against those positions (`Instrument::advance_hammers_external`), which refuses any step where a contact could engage.
- On a refusal, the stretch is rerun up to that step and the CPU takes over. Contacts, events and the frames around them always step on the CPU.

**Step as phases (`plan.rs`).** One contact-free step (damp, predict, finish with zero force, commit, damp) is linear in the bank state. The work splits into one item per bridge group (its string modes, damper pads and contacts) and one per board coordinate. Items meet only through small vectors, and every coupling is a product with a matrix composed on the host from the bank's own operators. The kernel needs 7 grid barriers per step with dampers engaged and 3 without. In exact arithmetic it is the CPU's step; only the rounding differs.

**Dampers in station space.** Each pad's chain of rank-one flows is written as three products, precomposed per pedal state (`plan::ChainOp`): u = S'v, c = T u, v += S c. The superposed form (`dampers::Prepared::apply_superposed`) couples pads through the board. A damper call therefore needs only each station's speed u. The changes accumulate as S Δ, which reaches the string velocities once per step.

**Precision.**

- Apple GPUs have no FP64, so Metal uses double-float pairs (about 48 bits), with fast math off.
- CUDA uses `double` by default.

**Grid barrier and residency.** All blocks must be resident at once.

- CUDA guarantees residency with a cooperative launch.
- Metal does not guarantee it, so the kernel launches at most one threadgroup per GPU core (`gpu-core-count`; more hung).
- The Metal barrier also aborts when a wait lasts far longer than any phase. The host then reruns the stretch from its saved state. The end-of-run line counts these reruns.

**On the CPU (`cores`).** The same phases run on the CPU cores (parallel.rs). Each thread owns a run of bridge groups and board rows, and the threads meet at a spin barrier between phases. Every item's arithmetic is the single-core emulator's, so states and traces are bit-identical to it.

## Measured (Apple M5 Ultra, 80-core GPU)

Strandbar piano, `--fast-dampers`:

| | 4 s excerpt (192,000 steps, 46 notes) | Full stem (17.4 s, 833,664 steps, 260 notes) |
|---|---|---|
| CPU (`FS_PIANO_GPU=off`, one core) | 208.8 s | 929.1 s |
| Metal | 22.1 s (97.5 % of steps on the GPU) | 94.5 s (97.8 % on the GPU) |

**Contact-free step cost (2048-step stretches from rest).**

| Backend | Dampers engaged | Dampers lifted |
|---|---|---|
| Metal | 80 µs | 28 µs |
| `cores`, 12 performance cores | 103–115 µs | 33–35 µs |
| `cpu`, one core | 613 µs | |

The `cores` figures were taken on a busy host.

**Agreement with the CPU render.**

- Output, 4 s excerpt: 1.8e-10 relative RMS difference (max 3e-8 Pa at a 3.8 Pa peak, below float32 resolution); 0.0000 dB in every octave.
- Output, full stem: 2.3e-10 relative RMS (max 6e-8 Pa); 0.0000 dB in every octave.
- Energy: same stored energy and damper loss to 9 digits; ledger closure 1.3e-14 J (excerpt), 1.5e-13 J (stem).
- Step by step (`FS_PIANO_GPU_CHECK=64`, 35,000 steps): state 1.7e-13 relative, board trace 7e-13, contact positions 7e-17 m, losses 2e-11.

The rest of the GPU run's wall time is the CPU's contact steps (about 1.1 ms each), the hammer replay and the microphone.

## Measured (RTX 4090, CUDA, `double`)

The 4 s excerpt on trj's RTX 4090 matches the CPU render made on the M5 Ultra:

- 1.2e-14 relative RMS (max 3.6e-12 Pa at a 3.8 Pa peak) and 0.0000 dB in every octave;
- the same stored energy to 9 digits; ledger closure 2e-14 J.

Speed was not measurable there. The host ran at load 160, so its CPU work took about ten times longer than usual: damper rebuilds took 49 s instead of under 0.3 s. The GPU stepped at 325 µs per step. Measuring on a quiet host is bead frankensim-y3yoo.

