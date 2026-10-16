# fs-soa-derive CONTRACT

## Purpose and layer

Layer: **UTIL** (proc-macro crate; zero dependencies — std `proc_macro`
only). The in-house `#[derive(Soa)]` macro (plan §5.3): turns a named-
field POD struct into a structure-of-arrays container built on fs-soa's
runtime. In-house because the Franken-only dependency law covers macro
dependencies too: no syn, no quote, no proc-macro2 — a hand-written
token walker and a text generator.

## Public types and semantics

- `#[derive(Soa)]` on `struct Name { … }` generates `NameSoa` with:
  one `fs_soa::FieldBuf` per leaf field (a nested container per
  `#[soa(nested)]` field), `new`/`with_capacity`/`len`/`is_empty`/
  `capacity`/`clear`/`reserve`, `push`/`get`/`set` (AoS
  scatter/gather, bit-exact), per-field accessors `{field}()` /
  `{field}_mut()` (dense aligned slices; nested: the inner container),
  `iter()` (gathered values), `field_views()` (dotted-path
  `RawView`s), `layout_descr()` (address-free JSON lines), plus
  `Default`, `fs_soa::SoaAble` and `fs_soa::SoaContainer` impls so
  containers nest.
- Helper attribute `#[soa(nested)]`: the field type must itself derive
  `Soa`; storage recurses at the type level via `SoaAble::Soa`.
- Supported shapes: named-field structs; nested POD structs (via the
  attribute); generics with bounds and const parameters (re-emitted
  verbatim); `fs_qty::Qty`-typed fields (plain 8-byte leaves — columns
  stay dimensionally typed).
- Generated accessors adopt the struct's visibility; leaf fields must
  be `Copy` (enforced post-expansion by `FieldBuf<T: Copy>`).

## Invariants

- Output is a pure function of the input token stream (no
  environment, filesystem, or randomness) — deterministic expansion.
- Generated code compiles under the workspace's clippy
  all+pedantic wall and `missing_docs` (every public item carries
  `#[doc]`; impls are `#[automatically_derived]`).
- Nested storage is driven in UFCS form (`<Storage as
  SoaContainer<T>>::…`) so callers never need trait imports and
  inference never ambiguates.
- Token rendering respects `Spacing::Joint` (`::`, `->` are never
  split — the bug class the first build caught).

## Error model

Unsupported shapes produce `compile_error!` with a structured message
— no silent fallbacks: tuple structs, unit structs, enums, unions,
lifetime parameters (containers own their storage), generic parameter
defaults, zero fields, field names colliding with the generated API
(reserved list in the message), unknown `#[soa(…)]` arguments.

## Determinism class

Bit-irrelevant (compile-time only); expansion is deterministic per the
invariant above. No floating point.

## Cancellation behavior

Not applicable: runs inside rustc, bounded by input size (single pass
over tokens plus string assembly).

## Unsafe boundary

None. `unsafe_code = "deny"` via workspace lints.

## Feature flags

None.

## Conformance tests

Exercised through fs-soa (a proc-macro crate cannot unit-test its own
expansion): `crates/fs-soa/tests/soa_battery.rs` covers plain, nested,
generic-with-bounds+const, and Qty-typed fixtures end to end;
`crates/fs-soa/tests/compile_fail.rs` pins all eight diagnostic paths
via the in-house offline cargo-check harness.

## No-claim boundaries

- No span-precise diagnostics: errors point at the derive site, not
  the offending token (string-based generation; revisit only if the
  ergonomics cost shows up in practice).
- No tuple-struct support (no field names to become accessors), no
  enums/unions (not SoA-decomposable), no lifetime-parameterized
  structs (POD only), no generic parameter defaults.
- No `#[soa(skip)]`, per-field alignment overrides, or rename
  attributes — add only with a consumer that needs them.
- Automatic recursion into un-annotated nested structs is impossible
  in a derive macro (no cross-item visibility) — `#[soa(nested)]` is
  the explicit, honest marker.
