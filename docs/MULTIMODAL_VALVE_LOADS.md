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
