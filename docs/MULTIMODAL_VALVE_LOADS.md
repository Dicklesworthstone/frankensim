# Referenced loads for multimodal valve paths

`MmLineBank::new` and `BrassVoice::new` accept `MmBankLoad::PerCombo(&loads)`.
The array has exactly one `MmLoad` per original duct, in the same order as the
ducts and labels. Existing `MouthLoad`, `MmLoad` and `&MmLoad` calls remain
shared-load calls with unchanged realization arithmetic.

A strict `MmLoad::Tabulated` retains its own declared mouth position, radius,
gas and DC impedance. Every combination is checked before frequency realization,
including combinations inactive at startup. Differing path lengths therefore
need independently referenced tables; assigning one table twice does not move
its physical plane. The original table admission, complete ladder bandwidth,
modal recurrence, shared FFT-size ladder and passivity tests remain in force.

```rust,ignore
let loads = [MmLoad::Tabulated(&open_table), MmLoad::Tabulated(&valved_table)];
let voice = BrassVoice::new(
    &ducts, &labels, &gas, MmBankLoad::PerCombo(&loads), &config, lip, cup_volume,
)?;
```

Valves use the existing outgoing-history lift, optionally with a crossfade.
Lip position/momentum and cup pressure/flow are not reset by a valve change.
The focused source regression compares each returned incoming sample with a
direct convolution of the actual outgoing history through both valve changes,
including waves sent before a combination was selected.

This supersedes the earlier fs-couple contract statement that per-combination
tables are outside the bank interface. It does not change the one-bore lip-file
format or its plane-wave mouth-flow observer, provide a source radiation bake,
qualify arbitrary tables as causal, or model physical valve travel. The broader
instrument-realism bead remains open. Native Rust execution of these new tests
has not been completed in the implementation environment.

## Sample-indexed performances

`brass_sequence::BrassSequence` binds these same paths and loads to a finite
sequence of the original `BrassControl` assignments. `BrassEvent::sample` is
the zero-based mechanics sample immediately before which the control applies.
Events at zero and at the last sample are valid; the terminal horizon is not
an additional sample. The callback splits at events without pre-rendering,
interpolating control times, or resetting lip, cup or acoustic history.

The wrapper implements the existing `PressureRenderer`. It emits **mouthpiece
pressure in pascals**, not a radiated microphone signal. `block_outgoing` also
returns the actual outgoing inlet wave for a separately specified acoustic
observer. Existing PCM and integer-decimation consumers can observe the source;
decimation does not retime its controls. The original voice's diagnostic block
statistics still depend on the caller's partition, unlike the physical samples.

Initial pressure, tension, active path, sample horizon and callback size are
explicit. Admission caps are 16 paths, 16384 events, 600 seconds and 65536
samples per callback. It rejects unordered events, duplicate assignments of
one control at one sample, invalid physical values and overlapping valve fades.
Different simultaneous controls retain authored order. A same-path valve
assignment remains the original bank's no-op; an actual fade must finish within
the horizon. Schedule validation precedes acoustic realization. Size failures
leave the trajectory untouched; a physical failure poisons continuation because
an earlier callback prefix may have advanced. There is no rollback, live gas
rebinding, ramp interpolation or whole-system passivity claim here.

A runnable manufactured two-path example writes a new diagnostic WAV without
overwriting an existing file or peak-normalizing its pressure:

```sh
cargo run -p fs-couple --example brass_valve_phrase -- mouthpiece.wav
cargo test -p fs-couple --lib brass_sequence::tests
```

Its geometry, lip and resistive terminal laws are analytic/manufactured fixtures,
not measured instrument or exterior-radiation data. Five new Rust tests compare
manual voice control with multiple callback partitions, validate schedule and
callback refusals, exercise poisoning, and reuse the original decimator. They
remain natively unexecuted here. An independent Python schedule-arithmetic
check passed 2004 partitions; that does not execute or qualify Rust acoustics.
