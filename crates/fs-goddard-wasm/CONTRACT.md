# fs-goddard-wasm — CONTRACT

Layer: L6 (boundary crate). Standalone `[workspace]` following the
`fs-wasm` / `fs-flyer-wasm` pattern: browser builds never depend on the
state of unrelated in-progress native crates.

## Purpose and layer

The boundary between the Goddard Rocket museum page (browser JS) and the
simulation plane (`fs-mbd` rigid-body physics), compiled to
`wasm32-unknown-unknown`. Its primary export is a source-bounded apparatus
step for US 1,102,653: torque-free primary/gyroscope poses, the printed firing
sequence, the Claim 2 `L/D` limit, and the Claim 7 ideal isolation probe. The
older liquid-nozzle export is retained only as an explicitly adjacent
interpretive calculation for other presentation contexts.
Native builds compile the same pure function for tests; `wasm-bindgen` is
confined to wasm32 via a target-specific dependency so native builds stay
dependency-clean.

## Public types and semantics

| Entry | Kind | Contract |
|---|---|---|
| `goddard_apparatus_step(elapsed_seconds, primary_spin_rpm, gyro_spin_rpm, tube_length_ratio, auxiliary_release_fraction, primary_charge_substantially_consumed, gyro_enabled) -> String` | wasm + native | Admits one bounded finite source-apparatus tuple, executes normalized torque-free poses through `fs_mbd::goddard::step_goddard_apparatus`, and returns quaternions, SI angular rates, and explicit Claim 1/2/7 states or a typed refusal. Ideal Claim 7 isolation applies only when the declared gyroscope is present and spinning; a stopped gyroscope shares the primary world rate. |
| `goddard_rocket_step(chamber_pressure_psi, chamber_temperature_k, specific_gas_constant_j_kg_k, heat_capacity_ratio, throat_area_cm2, expansion_ratio, ambient_pressure_pa) -> String` | wasm + native | Resolves an ideal calorically-perfect nozzle through `fs_mbd::goddard::step_goddard_rocket`, with explicit chamber thermodynamics and ambient pressure. The sonic throat determines mass flow. Returns exit T/p/Mach/speed, choked mass flow, momentum and pressure thrust, total thrust, specific impulse and characteristic velocity inside `{"ok":{...}}`, or a typed refusal. |

The previous four-argument rocket signature is replaced. No repository browser
caller uses it; external callers must provide the seven explicit inputs and
consume `mass_flow_kg_per_sec` as a result. The source-bounded apparatus export
keeps its existing signature and semantics. For example,
`goddard_rocket_step(350, 2800, 365, 1.24, 4.2, 3.5, 0)` declares an assumed
calorically-perfect chamber in vacuum; those gas properties are caller inputs,
not historical or empirically validated propellant data.

## Invariants

1. **Shape stability.** Each emitted record keeps exactly the fields of its
   owning result under their struct names inside an `"ok"` object; no field is
   silently repurposed by the renderer.
2. **Dependency purity.** Runtime dependencies are the workspace crate
   `fs-mbd` only (Decalogue P1); serialization is explicit formatting, not
   a serde edge.
3. **Deterministic content.** Identical inputs produce byte-identical
   output on the same platform; no wall clock, randomness, host identity,
   or allocation-order dependence reaches the string.
4. **Fail-closed admission.** Non-finite or above-cap inputs refuse before the
   owning kernel runs. The source apparatus deliberately admits zero spin and
   below-Claim-2 ratios inside a bounded break probe, returning explicit claim
   failure rather than laundering the state as valid. Non-finite outputs refuse
   rather than being rendered as JSON `null`.

## Error model

The boundary returns stable refusal codes with a human-readable message and
ranked repairs. `non-finite-input` covers NaN and infinity;
`input-outside-domain` covers values outside the documented bounds;
`non-finite-output` catches a broken downstream result; and
`rigid-body-refusal` preserves a refusal from the generic owner. The source
apparatus admits 0–600 s, 0–1,200 primary RPM, 0–60,000 gyro RPM, `L/D` from
1–12, and release fraction from 0–1. The ideal-nozzle boundary admits positive
chamber pressure up to 10,000 psi, temperature up to 10,000 K, gas constant up
to 10,000 J/(kg K), throat area up to 100,000 cm², gamma in `(1,5/3]`, area
ratio in `[1,1000]`, and nonnegative ambient pressure up to 100 MPa. These are
input/work bounds, not material-validity claims. The kernel additionally
requires ambient pressure no higher than its isentropic exit pressure;
`unsupported-back-pressure` explicitly refuses overexpanded operation instead
of ignoring shocks or separation. Vacuum is admitted without evaluating gas
at zero pressure. `output-outside-domain` covers a required positive quantity
lost to overflow or underflow; `nozzle-refusal` preserves other kernel failures.

## Determinism class

Bit-identical replay on the same ISA for identical inputs. Cross-ISA
identity is not claimed (floating-point rendering follows the platform's
shortest-roundtrip display).

## Cancellation behavior

Synchronous single-step entry. The rocket area-Mach inversion has a hard cap
of 80 allocation-free iterations. A JS caller can discard the returned string
if the result is no longer needed; this is not an asynchronous solver.

## Unsafe boundary

None: `#![forbid(unsafe_code)]` is declared in Cargo.toml lints and the
crate contains no unsafe block.

## Feature flags

None. The wasm32 binding is selected by target architecture
(`cfg_attr(target_arch = "wasm32", wasm_bindgen)`), not by a feature.

## Conformance tests

Native unit tests exercise a valid source apparatus state, both claim-failure
probes, non-finite and bounded-domain refusals, an ideal-nozzle result, its
lower/upper admission boundaries, underflowed output refusal, sonic/vacuum
operation, unsupported back pressure and exact replay. The owning fs-mbd
physics tests independently check the analytic Mach-2 nozzle, near-sonic and
gamma-to-one limits, mass/enthalpy balance, area and pressure scaling, and
pressure thrust. The dependency-policy/contract gates remain complementary
evidence.

## No-claim boundaries

This crate is a presentation seam. It makes no claim that either embedded
Goddard model is validated, that visitor-entered spin speeds are source values,
that normalized poses establish mass properties or loads, or that the rendered
numbers carry uncertainty bounds. The source apparatus returns no thrust,
Mach, liquid propellant, de Laval, or trajectory quantity. The older liquid
export must never be described as a mechanism claimed or disclosed by US
1,102,653. Its declared constant heat capacities are a reduced-model
assumption, not a chemistry or high-temperature real-gas solution; use the
source-backed `fs-material` thermochemical nozzle for temperature-dependent
properties. No tank depletion, shock/separation, viscous loss or vehicle
trajectory is inferred. Physics claims belong to the owning kernels and
their evidence, never to this renderer.
