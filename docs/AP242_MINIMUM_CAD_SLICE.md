# AP242 Minimum CAD Slice for the Cooling Vertical

Status: scope contract for
`frankensim-extreal-program-f85xj.11.2`; implementation remains staged under
the existing AP242 program `frankensim-leapfrog-2026-program-i94v.2.2`.

## Decision

Cooling 0.1/0.2 is the first bounded consumer of the repository's AP242
program. It does not own a second schema kernel or a vertical-specific STEP
interpreter. The minimum semantic slice is:

1. one root assembly with explicit component occurrences;
2. a stable product/part identity and retained source entity identity for each
   occurrence;
3. an explicit proper-rigid placement for every occurrence, including parent
   and child frame identity, transform convention, and handedness;
4. one unambiguous source length-unit context normalized through `fs-qty`
   before geometry or placement admission;
5. part names and bounded attributes that can be mapped deliberately to
   `fs-project` region bindings and geometric selectors; and
6. finite faceted or tessellated geometry for every admitted part, with the
   source representation, tolerance, conversion receipt, and provenance
   retained separately from the assembly graph.

The first reference assembly is an enclosure, board, and heat sink represented
as three distinct occurrences. Coincident or repeated geometry must not merge
occurrence identity.

## Authority and sequencing

The existing leapfrog program remains authoritative:

| Concern | Authority | Cooling dependency |
| --- | --- | --- |
| Bounded versioned Part-21/AP242 entity decoding | `frankensim-leapfrog-2026-program-i94v.2.2.1` | The vertical consumes the supported schema subset; it does not fork it. |
| Product, assembly, occurrence, configuration, and lineage semantics | `frankensim-leapfrog-2026-program-i94v.2.2.2` | The vertical requires enclosure, board, and sink occurrences with distinct identities. |
| Exact/tessellated geometry and validation-property binding | `frankensim-leapfrog-2026-program-i94v.2.2.3` | The minimum slice admits per-part tessellation; it does not call a mesh exact geometry. |
| Unsupported entities and semantic-difference receipts | `frankensim-leapfrog-2026-program-i94v.2.2.7` | Unsupported AP242 semantics refuse or remain opaque with an explicit loss record. |
| Existing strict faceted STEP bridge | `frankensim-ext-cad-step-interchange-8uxb` | This is the current fallback and lower-layer implementation base, not full AP242 support. |
| Persistent project assignment and named-group path | `frankensim-extreal-program-f85xj.6.3` | Names and attributes must be bound explicitly; source labels are not trusted region identity. |
| Product import acceptance lane | `frankensim-extreal-program-f85xj.6.11` | Its future golden is the integration bar; current component tests do not satisfy it. |
| Permissioned supplier-pair corpus | `frankensim-extreal-program-f85xj.11.6` | The corpus must falsify this scope before any supplier-coverage claim. |

Related dependency edges in Beads connect this scope record to each authority.
The edge records relevance, not completion or proof.

## Gap children

The existing AP242 program covered the broad semantic areas but did not assign
two import-boundary seams precisely enough for the vertical. They therefore
live under the leapfrog tree:

- `frankensim-leapfrog-2026-program-i94v.2.2.10` binds AP242 occurrence
  transforms into explicit `fs-geom` placements. It depends on the existing
  assembly and geometry owners and must refuse ambiguous, singular, non-rigid,
  or unsupported frame mappings.
- `frankensim-leapfrog-2026-program-i94v.2.2.11` binds supported STEP length
  units into `fs-qty` at import admission. It depends on the schema kernel and
  geometry binding and must refuse missing, conflicting, ambiguous, or
  unsupported unit contexts.

Neither gap is refiled under EXTREAL. Closing this scope task does not close or
implement either child.

## Current executable fallback

The current product import surface is narrower than the minimum AP242 slice:

| Input surface | Current executable support | Evidence and boundary |
| --- | --- | --- |
| STL | CLI/library quarantine, bounded repair, promotion, persistent assignment, and ledger lineage | `crates/fs-cli/tests/import.rs::g0_dirty_stl_promotes_assigns_and_retains_complete_lineage` exercises this route. STL carries positions only. |
| OBJ | CLI/library `v`/`f` subset through the mesh quarantine path | Implemented in `fs-io`; not exercised by the current product import e2e fixture. Texture coordinates, normals, and materials are not retained as semantics. |
| PLY | CLI/library ASCII and little-endian mesh subset through the quarantine path | Implemented in `fs-io`; not exercised by the current product import e2e fixture. Unsupported properties are not product attributes. |
| Strict faceted STEP | Caller-selected triangular `FACETED_BREP` closure, topology quarantine, estimated SDF handoff, assignment, and separate receipts | `crates/fs-cli/tests/import.rs::g0_faceted_step_import_retains_both_receipts_and_repaired_mesh` exercises this route. The caller supplies the root entity, length unit, and sampling spacing. |
| AP242 semantic assembly | Not implemented | No discovery of products, occurrences, transforms, unit contexts, or part attributes is claimed. |
| 3MF/GLB/VTK | Export-only surfaces in `fs-io` | They are not current product import formats. |

The CLI accepts mesh sources only when the canonical project row declares
`stl`, `obj`, or `ply`; it does not infer a format from a suggestive file name.
The strict STEP route recognizes a pinned triangular resource subset rather
than a complete AP242 representation. Caller/project unit declarations and
caller-supplied named face groups are retained policy, not decoded CAD
semantics.

This fallback is recorded as maturity L1. The two named CLI tests demonstrate
that the retained routes execute; they do not raise the capability to an
integrated AP242 workflow. The future `f85xj.6.11` lane and permissioned
`f85xj.11.6` corpus remain required before stronger integration or
supplier-coverage language.

## Explicitly staged semantics

The following are outside the minimum slice and remain with the broader
leapfrog program:

- trimmed NURBS and exact B-rep interpretation;
- PMI, GD&T, datum, and surface-texture semantics;
- material, process, lot, harness, and kinematic semantics;
- persistent topological naming across revisions;
- inferred joints or motion constraints from static placements;
- certified physical/CAD sameness or continuum coverage; and
- export/import semantic round-trip claims beyond the admitted subset.

Unsupported or absent semantics must be represented as a structured refusal,
opaque extension, or semantic-difference receipt. Visual similarity, matching
part names, and a clean triangle soup are not substitutes for those claims.

## Acceptance and falsifiers

The minimum slice is ready for promotion only when all of the following are
retained:

1. an enclosure/board/sink supplier fixture whose assembly occurrences,
   transforms, unit context, names/attributes, and per-part tessellation are
   independently enumerated;
2. deterministic import receipts that reproduce those identities and
   placements under entity-order permutation;
3. unit-rescaling tests proving that equivalent source units produce the same
   canonical `fs-qty` geometry;
4. ambiguity, missing-unit, unsupported-entity, singular/non-rigid transform,
   repeated-part, and cancellation refusals;
5. a full `f85xj.6.11` product lane binding the admitted parts to project
   regions and preserving all lower receipts; and
6. DSR evidence and retained artifact hashes on a stable source snapshot.

The supplier corpus in `f85xj.11.6` is not yet retained, and `f85xj.6.11` is
not currently green. This document therefore fixes the target and the honest
fallback; it does not claim that the minimum AP242 slice has landed or that it
has been validated against real supplier files.
