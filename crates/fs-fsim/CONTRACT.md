# fs-fsim Contract

## Purpose and layer
`fs-fsim` (layer `L6`, HELM) owns the wasm32-clean FrankenScript AST (`ast`), isomorphic canonical s-expression (`sexpr`) and JSON (`json`) syntaxes, and structured `IrError` shared by `fs-ir` and downstream runtimes.

## Public types and semantics
- `ast::{CountUnit, CountValue, Node, NodeKind, Span}`: span-bearing typed AST over `fs-qty` dimensions.
- `sexpr` and `json`: deterministic total parsers and canonical checked printers with lossless round-trip isomorphism.
- `IrError`, `IrErrorKind`: structured syntax and shape refusals carrying byte spans and fix hints.

## Invariants
- Pure Rust, zero filesystem, zero `fs-ledger`/SQLite dependencies; compiles natively and for `wasm32-unknown-unknown`.
- `fs-ir` re-exports `ast`, `json`, `sexpr`, `IrError`, and `IrErrorKind` verbatim from `fs-fsim`.

## Error model
All parse and shape failures return `IrError { span, kind, detail, hint }`. Depth caps (`sexpr::MAX_DEPTH = 256`) prevent stack overflow on hostile input.

## Determinism class
Pure deterministic parsing and canonical formatting; identical input bytes produce identical ASTs and output bytes.

## Cancellation behavior
Synchronous parsing and printing bounded by input length and `sexpr::MAX_DEPTH` (`256`); no background tasks or `Cx` cancellation polling.

## Unsafe boundary
`#![forbid(unsafe_code)]`. Zero `unsafe` blocks.

## Feature flags
None.

## Conformance tests
- `fs-fsim` `tests::sexpr_and_json_round_trip` (`src/lib.rs`).
- `fs-ir` unit and integration test suites exercising the re-exported `ast`, `sexpr`, `json`, and `IrError` API.

## No-claim boundaries
Syntax parsing and canonical formatting establish structural validity only, not scientific authority or material provenance.

## Planned (R2b)
- Cut R2b: move `fs-cli`'s `.fsim` assembly codec and unified `UNIT_DIMS`/`dims_to_unit` into `fs-fsim`.
