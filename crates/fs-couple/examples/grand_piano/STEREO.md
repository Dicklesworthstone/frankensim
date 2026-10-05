# Two spatial microphones, one physical piano

`--microphone-right x_m,y_m,z_m` adds a second Rayleigh receiver and writes
frame-interleaved left/right PCM16. The original `--microphone` position, or
its unchanged default `(0.675,1,1)` m, is the left channel. Omitting the new
option preserves mono output. Coordinates use the supplied board's acoustic
frame, including the existing flattened projection for crowned boards.

```sh
cargo run --release -p fs-couple --example grand_piano -- \
  --preset steinway-d --midi performance.mid --dampers estimated \
  --microphone 0.2,0.8,1.0 --microphone-right 1.2,0.8,1.0 \
  --duration 30 --render stereo.wav
```

Each microphone uses its own distributed signed surface projection, geometric
travel times, causal filter and pressure history. Both receive the SAME loaded
soundboard velocity trace from every mechanical substep. The engine and score
advance once per audio frame, not once per receiver. No second instrument,
voice stealing, panning law, copied mono signal or decorrelation effect is used.
Coincident positions deliberately produce identical channels; swapping the
positions swaps the observations without changing the physical performance.

Both channels use the same PCM full-scale pressure, defaulting to 2 Pa. Set
`--pcm-full-scale-pa` above the reported pressure peak when the default range
is too small; this changes only the PCM conversion, not the physical pressure
or microphone position. An over-range render refuses before writing a clipped
WAV. No independent channel normalization is performed. The reported peak is
across both channels, and clips count scalar channel samples. Duration, score
event indices, sample rate and block-error progress remain frame counts,
regardless of channel count.
Input scale, geometry, hammer cards, damper cards, MIDI/CSV controls, tuning and
mechanical limits are unchanged. Stereo requires a geometric pressure render;
`--diagnostic-volume` is not silently promoted to physical spatial sound.

For hosts, `AudioStream::new_stereo` prepares both receivers and
`render_interleaved_block` writes `[L0,R0,L1,R1,...]`. Arbitrary complete-frame
block sizes retain both histories. An odd-sized stereo buffer refuses before
consuming controls or mechanics and can be retried with a complete buffer.
The old `render_block` remains mono-only and refuses an implicit downmix.
A physics/observer failure preserves only the completed frame prefix, silences
the entire remainder and latches the stream. This includes an observer failure
after mechanical acceptance; resuming from that state would desynchronize clocks.
The successful host loop adds no allocation. This is not a measured real-time
or audio-device-backend claim.

The original Rayleigh assumptions remain: prescribed linear radiation into an
infinite baffled half-space; no room/lid/cabinet scattering, binaural HRTF,
radiation backreaction or measured SPL/realism claim. Stereo does not enlarge
the board's structural bandwidth or validate its mesh. Native regressions
compare both channels to independent mono receivers, mechanical/work identity,
block partitioning, coincident receivers, fault framing and PCM layout.

`--acoustic-refinement-levels 0..3` independently refines the flat P1 board's
Rayleigh integration. Each level splits a structural triangle into four
equal-area acoustic cells, retaining the positive three-point rule and evaluating
the same signed P1 modal displacement on each cell. Structural nodes, mass,
eigenpairs, bridge projections, source dynamics and total radiating area stay
fixed. Level zero preserves the original points and floating-point path.
The existing 120,000-point receiver budget still refuses oversized requests
before eigenanalysis. The flag requires geometric pressure output and excludes
diagnostic volume and edge-cubic controls, including explicit zero. Crowned
fields reject nonzero levels; zero preserves their original integration path.
It provides a spatial integration convergence
test; it does not extend structural bandwidth or establish a measured piano
match. Compare successive levels at fixed source controls and microphone
positions before interpreting pressure changes as improved accuracy.

`--receiver-pressure-csv receivers.csv` exports total unquantized pressure
with `sample,time_s,pressure_left_pa` and, for stereo, `pressure_right_pa`.
It reads the same completed pressure buffer encoded into the WAV, so it adds
no observer or mechanical step and allocates no per-mode history. It also
supports music. Use `--modal-pressure-csv` when signed mode attribution is
needed; use this smaller export for repeated receiver and PCM comparisons.
