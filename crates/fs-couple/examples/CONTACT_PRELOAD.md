# Settle contact-loaded structures, then release them

```bash
cargo run -p fs-couple --bin music_render -- \
  modal crates/fs-couple/examples/modal-contact-preload.performance \
  /tmp/contact-preload.wav --block 37
```

Normal-contact performance versions 3 and 4 now accept `static-preload` on every
voice. Supply zero Q/V fields and the actual held physical port forces. The
loader solves the **complete nonlinear contact-loaded equilibrium**, retaining
every modal basis and bilateral spring. It does not settle the components
independently or simulate a fictitious warm-up. Force events later change those
loads without resetting the resulting displacement or contact state.

The example has three bodies, two nonlinear contacts, and an unforced center
body that is the only pressure observer. It starts settled under opposing loads.
Both outside loads release before sample 37, then the left body is driven again
at sample 300 and released at sample 340. The initial held interval should have
only numerical-roundoff motion, rather than an unwanted startup impact. The
4801-sample output retains its short final callback and declared 0.2 Pa scale.
These are authored reduced-model parameters, not a calibrated instrument.

## Native API and limits

`CoupledModalSystem::initialize_contact_equilibrium` accepts all contacts and
held generalized forces, solves stationary balance, and returns total stored
energy including the contact potentials. It consumes no sample. Attach those
same contacts using `ContactModalSystem::new` or `MultiContactModalSystem::new`
and continue with the same held forces until release.

For file-driven v4 performances, `multi_contact_limits` bounds both the static
setup/sweeps and subsequent dynamic work. Every contact keeps its declared root,
force, penetration and residual limits. V3 has one contact and needs one joint
coordinate sweep; it uses `coupling_limits`' setup-term ceiling and the existing
per-contact root cap. Budgets are never raised automatically. Mixed retained and
preloaded components, nonzero authored Q/V on a preload, and unconverged solves
refuse before output creation. The friction-enabled v5 path still requires
`retain-state`; no tangential static sticking/preload has been inferred.

The solution uses the static inverse stiffness, not the audio-step compliance.
Bilateral rest offsets enter the free equilibrium but not its derivative.
Contact activity follows the actual solved gaps. Physical limits apply to the
final restrained state; a large unrestrained prediction is not falsely published
or used to reject a smaller admissible contacted solution. Contact-inclusive
energy, every contact law and the actual modal force balance are checked before
publication. Cancellation/refusal preserves the accepted network and clock.

This is an explicitly pre-window settled load, not a measured work history or
an in-run equilibrium reset. Fixed mass-normalized linear components and the
existing normal power-law potentials remain the model assumptions.

## Stable near-equilibrium contact evaluation

The previous scalar contact adapter subtracted nearly equal potential energies.
At settled states that amplified roundoff enough to stall the force solve even
when the actual motion was negligible. `fs_dcontact::OpeningContactStep` now
computes the same one-coordinate potential secant through a stable ratio and
`fs-math`'s deterministic `expm1`; contact entry/exit still accounts for the full
opening displacement. The modal contact paths use this owner. The older reed
adapter and general pHS discrete-gradient implementation are not changed.

No force or energy tolerance is loosened. This numerical correction changes
rounding and can change PCM bytes of older contact renders; byte equality with
pre-correction builds is not claimed. Input schema/units, physical contact laws,
force scheduling and the WAV encoder are unchanged.

Focused checks (native Rust execution was unavailable in the authoring environment):

```bash
cargo test -p fs-dcontact --lib opening
cargo test -p fs-couple --test modal_contact_preload --test music_render_preload
cargo test -p fs-couple --test modal_contact --test modal_multi_contact --test modal_contact_friction
```
