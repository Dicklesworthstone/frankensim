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

Both channels retain the existing 2 Pa PCM full-scale. No independent channel
normalization is performed. The reported peak is across both channels, and
clips count scalar channel samples. Duration, score event indices, sample rate
and block-error progress remain frame counts, regardless of channel count.
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
