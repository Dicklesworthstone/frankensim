# MIDI scores through the physical grand piano

`--midi score.mid` imports a Standard MIDI File into the existing sample-accurate
performance and render path. It does not add a synthesizer, sample bank, voice
allocator, authored oscillator, gain envelope or alternative acoustic model.

```bash
cargo run --release -p fs-couple --example grand_piano -- \
  --preset steinway-d --midi performance.mid --midi-channel 1 \
  --midi-velocity-max-m-s 4.5 --duration 30 --render piano.wav
```

Select the actual piano channel with `--midi-channel 1..16` (default 1). Tracks
are merged in time, but only the selected channel moves this piano's keys or
pedals. Other-channel messages are counted and ignored: a drum channel cannot
accidentally strike a nonexistent piano string or release this piano's sustain.
The supplied scale, board, felt cards, hammer masses, flexible preset shanks and
microphone options remain independently selectable. `--midi` requires `--render`
and excludes `--performance`, `--note` and `--velocity`; MIDI-only mapping flags
without a MIDI file are errors. The score path is protected by the same direct
input/output path-collision checks as the existing physical input files.

## Explicit physical mapping

MIDI velocity is dimensionless. This importer deliberately declares a simple,
**uncalibrated** mapping: velocity 1..127 becomes
`maximum_velocity_m_s * velocity / 127` in post-escapement hammer speed. The
default maximum is 4.5 m/s; `--midi-velocity-max-m-s` accepts finite values in
(0,8] m/s. It changes mechanical input energy, not pressure gain. This is not
a measured controller curve, a key-velocity measurement or a full grand action.
Use the existing CSV `jack_staccato`/`jack_legato` controls for physical jack
force input instead of inventing a MIDI-byte-to-newton equivalence.

Note-on with zero velocity is note-off. Note-off releases the physical key;
release-velocity bytes are not modeled. Overlapping note-ons of the same key
refuse: one physical key cannot allocate independent voices. The engine also
retains its actual hammer-rearm and contact-validity checks for fast repetitions.
Missing keys refuse rather than transposing, substituting a default scale or
quietly dropping notes.

CC64 operates sustain with the MIDI switch threshold of 64 by default. Explicit
`--midi-half-pedal` instead interprets CC64/127 as this engine's physical pedal
travel; this is an additional mapping assumption, not an automatic claim that
the controller recorded continuous pedal depth. CC66 and CC67 use the same
switch threshold for sostenuto and una corda. Same-time events preserve source
order, including whether sostenuto captures a just-pressed key.

CC123 releases held keys but respects the current pedals. CC121 resets the
three mapped pedal controls without erasing resonator state. At the last
end-of-track timestamp the importer releases any still-held keys and pedals,
so the remaining render is physical ringdown. It never mutes PCM or clears
string, soundboard, felt or contact history. These final releases are reported.

Selected-channel non-center pitch bend and CC120 All Sound Off **refuse** rather
than retuning string frequencies or zeroing mechanical state. Program changes,
pressure, volume/expression and other unmapped channel messages are counted and
ignored. SysEx is length-checked, counted and ignored; no device-specific setup,
General MIDI timbre selection or synthesizer effect is inferred. The render
prints the selected channel, mapping, ignored-message counts and score endpoint.

## Timing, bounds and evidence

The reader accepts format 0 and simultaneous format 1 with up to 256 tracks.
PPQN timing includes the default 500,000 microseconds per quarter and the full
tempo map; format-1 tempo events must be in the first, conductor track. SMPTE
-24/-25/-30 clocks and -29 (30000/1001 fps) are supported. SMPTE timing is not
retimed by tempo events. The timeline starts at file tick zero; absolute SMPTE
Offset metadata does not position the file on an external transport.

Track merging is stable by tick, then file track order, then in-track order.
Integer rational time is accumulated before converting each absolute timestamp
to the first output sample at or after that time. There is no accumulated
rounding of individual deltas or tempo segments. This is output-sample-accurate
control dispatch, not fractional-sample reconstruction of a key action.

Length-delimited metadata, header extensions and unknown chunks are skipped;
channel running status is implemented and is cleared by meta/SysEx events.
Malformed/truncated input, missing end-of-track, invalid data bytes, overlong
variable-length quantities, zero tempo, format 2 and inconsistent track counts
refuse before the performance is published. Cold import is limited to 32 MiB,
4,096 chunks and one million raw events **and** one million expanded physical
controls. File reads enforce the byte limit before parser allocation growth.

The **entire file**, including its final end-of-track, must fit `--duration`;
notes and releases are never silently truncated. Add time after that endpoint
for ringdown and microphone/filter delay. The existing maximum duration remains
120 seconds; this change is an offline score importer, not a live MIDI device,
DAW transport, MIDI 2.0 implementation or real-time-performance qualification.

G0/G3 tests exercise tempo maps, fractional timing, SMPTE drop-frame, running
status, malformed prefixes, channel isolation, pedal/source ordering, source
end releases and physical-control admission. An integration test compares MIDI
against equivalent SI CSV controls through the actual preset hammer/contact/
string/board engine, including pedal release and its energy balance. These are
native Rust regressions; they require execution through the repository's
DSR/RCH development lane before a passing native-test claim is made.

Format and controller references: MIDI Association Standard MIDI Files v1.0
(RP-001), and the MIDI 1.0 Control Change Messages table. These establish file
and message semantics, not the uncalibrated physical velocity mapping above.
