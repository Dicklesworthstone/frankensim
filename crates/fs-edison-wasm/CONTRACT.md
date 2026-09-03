# fs-edison-wasm — CONTRACT

Layer: L6 browser boundary. This standalone workspace is a narrow patent
binding around the generic `fs-conduction` incandescent radiative balance.

## Purpose and layer

Admit the declared operating-point inputs used by the US 223,898 exhibit and
serialize `fs_conduction::incandescent::solve_incandescent_radiative_balance`
for browsers. The generic crate owns the law; this crate owns only admission
and stable envelope serialization.

## Public types and semantics

| Entry | Kind | Contract |
|---|---|---|
| `edison_radiative_step(...) -> String` | wasm + native | Accepts voltage, hot resistance, radiating area, emissivity, and ambient temperature in SI; returns current, Joule power, equilibrium temperature, rerun radiative power, and closure in `{"ok":{...}}`; every invalid domain returns `{"refusal":{...}}`. |

## Invariants

1. `V²/R` is balanced against `εσA(T⁴ - T_ambient⁴)` using the generic
   `fs-conduction` Stefan-Boltzmann constant.
2. The boundary never invents resistance, filament dimensions, emissivity, or
   ambient temperature. Every one crosses the call explicitly.
3. Returned energy closure is re-evaluated from the returned temperature.
4. Input refusal never returns an `ok` payload or silently clamps a quantity.

## Error model

Refusals distinguish non-finite input, negative voltage, non-positive
resistance/area/ambient temperature, and emissivity outside `(0,1]`. Each
refusal carries a repair.

## Determinism and cancellation

The bounded synchronous query is pure and deterministic and has no useful
cancellation boundary.

## Unsafe boundary and features

No unsafe code. The wasm binding is selected only by target architecture.

## No-claim boundaries

This rung omits lead conduction, residual-gas conduction, convection,
temperature-dependent resistance, material aging, useful visible-light
efficacy, and lifetime. The US 223,898 source prints resistance and geometry
examples but not a calibrated commercial operating point; exhibit choices are
identified as declared or illustrative inputs.
