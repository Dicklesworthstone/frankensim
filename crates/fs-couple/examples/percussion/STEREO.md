# Two spatial microphones on one drum or cymbal

Add `--microphone-right x_m,y_m,z_m` to any existing `-mic` command. The
original positional receiver, or its unchanged default `(0.08,0.05,0.35)` m,
is the left channel. The new comma-separated position is the right channel.
Output is one frame-interleaved left/right PCM16 WAV at 48 kHz per channel.
Omitting the option preserves mono; it does not change CSV or far-field WAV.

```sh
cargo run --release -p fs-couple --example percussion -- \
  snare-mic 48000 20 0.08 0.05 0.35 \
  --microphone-right -0.12,0.05,0.4 --cavity-modes \
  --second-stick-position-m -0.05 0.02 --second-stick-speed-m-s 1.6 \
  > stereo-snare.wav
```

Positions use the same acoustic geometry frame as the existing microphone.
For drums the origin is at mid-depth, with the batter at positive z. For a
cymbal it is the supplied shell frame. Each receiver must lie outside the
same enclosing sphere and leave at least two output samples of propagation.
Both receivers are admitted before any BEM source solve or mechanical step.
Stereo is deliberately refused on mechanics CSV and far-field `-wav` commands
rather than silently changing their meanings. Use `-mic` for spatial stereo.

The triangle boundary and modal source velocities are unchanged. At each
sampled frequency the source BEM batch is solved ONCE, then the existing
exterior Green representation evaluates both observation points. Each point
retains its own complex response, proper causal fit, held-out error check,
propagation delay and runtime history. A failed fit in either channel refuses
the entire candidate. No shared gain curve, panning law or dual-mono copy is
substituted. The frequency window and fit tolerances remain unchanged.

During rendering, one mechanical block supplies one interval-acceleration
trace and one causal decimator. BOTH receiver banks consume that same result;
adding a microphone never steps the stick, snare, head, felt or gas twice.
Both independent player-force files, two sticks, supplied geometry, mufflers,
sealed cavity drag and admitted analytic/nonlinear images remain available.
The existing full-scale pressure is applied identically to both channels;
there is no channel normalization. Frame count and duration are unchanged,
while the WAV payload doubles. Clip counts count scalar channel samples.

All fitting precedes playback. The offline front door publishes the complete
WAV only after successful mechanics, both observers and PCM admission. A
failure may follow accepted mechanics, so discard the candidate experiment
rather than retrying its acoustic state. It is not a transactional streaming
or real-time audio-device interface.

The original model limitations remain: one-way, linear, undeformed exterior
acoustics, the declared 40..1640 Hz fit window, no radiation feedback, room,
headphones/HRTF, calibrated SPL or full-band adequacy claim. **Vented exterior
audio is still refused**; stereo does not supply the missing aperture radiation.
Native regressions compare shared BEM bakes with independent receiver bakes,
actual drum state/time/work with mono execution, interleaved PCM, receiver and
CLI admission. Test-only authored filters isolate frame plumbing; they are
never selectable in the production render path.
