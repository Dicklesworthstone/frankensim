# fs-fsim Contract

## Purpose and layer
`fs-fsim` (layer `L6`, HELM) owns the wasm32-clean FrankenScript AST (`ast`), isomorphic canonical s-expression (`sexpr`) and JSON (`json`) syntaxes, structured `IrError`, the single `.fsim` assembly codec (`codec`), and the unified SI display unit table (`UNIT_DIMS`, `dims_to_unit`) shared by `fs-cli`, `fs-ir`, and `fs-runtime`.

## Public types and semantics
- `ast::{CountUnit, CountValue, Node, NodeKind, Span}`: span-bearing typed AST over `fs-qty` dimensions.
- `sexpr` and `json`: deterministic total parsers and canonical checked printers with lossless round-trip isomorphism.
- `IrError`, `IrErrorKind`: structured syntax and shape refusals carrying byte spans and fix hints.
- `codec::{AdmittedAssembly, AssemblyAnalysisMode, CodecRefusal, MAX_PACKS, PackCatalog, PackEntry, PackKind, admit_assembly_source}`: single shared `.fsim` parser, Five Explicits admission, and `fs-matdb` card-pack store builder used by both `fs-cli` and `fs-runtime`.
- `UNIT_DIMS`, `dims_to_unit`: canonical 17-unit SI dimension-to-unit mapping (exact union of `fs-cli` and `fs-runtime` tables) re-exported by both crates.

## Invariants
- Pure Rust, zero filesystem, zero `fs-ledger`/SQLite dependencies; compiles natively and for `wasm32-unknown-unknown`.
- `fs-ir` re-exports `ast`, `json`, `sexpr`, `IrError`, and `IrErrorKind` verbatim from `fs-fsim`.
- `PackCatalog::from_raw_bytes` enforces `MAX_PACKS` (`128`), deduplicates identical packs by content hash, and refuses conflicting cards or unknown magic bytes.

## Error model
Syntax and shape failures return `IrError { span, kind, detail, hint }`; `.fsim` admission and card-pack failures return `CodecRefusal { code, what, fix }`. Depth caps (`sexpr::MAX_DEPTH = 256`) and `MAX_PACKS` (`128`) guard hostile input.

## Determinism class
Pure deterministic parsing, admission, and canonical formatting; identical input bytes and card packs produce identical ASTs, `DraftAssembly` structures, and canonical output bytes.

## Cancellation behavior
Synchronous parsing and admission bounded by input length, `sexpr::MAX_DEPTH` (`256`), and `MAX_PACKS` (`128`); no background tasks or `Cx` cancellation polling.

## Unsafe boundary
`#![forbid(unsafe_code)]`. Zero `unsafe` blocks.

## Feature flags
None.

## Conformance tests
- `fs-fsim` `tests::sexpr_and_json_round_trip_and_dims_to_unit_match` (`src/lib.rs`) and `codec::tests` (`.fsim` admission, pack dedup/cap, and round-trip unit tests).
- `fs-ir` unit tests exercising the re-exported `ast`, `sexpr`, `json`, and `IrError` API.
- `fs-cli` `assembly_run` lib and integration tests exercising the shared `fs_fsim::admit_assembly_source` path.

## No-claim boundaries
Syntax parsing, card-pack decoding, and structural admission do not mint scientific authority or run physics stepping; those belong to `fs-matdb` cards and `fs-system`.
