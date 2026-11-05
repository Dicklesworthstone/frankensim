# SAFETY: fs-feec highorder fma capsule

## What is unsafe here

`unsafe fn apply_mono_x86`, `unsafe fn apply_tensor_mass_x86`, and
`unsafe fn apply_vector_mass_x86`
(only because `#[target_feature]` requires it), each with one
`unsafe { }` call site in its dispatcher, compiled only for x86-64.

## Why it is sound

- The dispatcher verifies `avx2` and `fma` via
  `is_x86_feature_detected!` immediately before the call — the only
  precondition `target_feature` imposes.
- The bodies are `TensorSpace::apply_mono_body` and
  `TensorMass::apply_body` with its rectangular contractions, and
  `VectorMaterialMass::apply_body` with its polynomial contractions, ordinary safe
  Rust (slice indexing, `mul_add`), marked `#[inline(always)]` so its
  code is GENERATED inside this function under the enabled features;
  no intrinsics, no pointers, no transmutes anywhere.
- Bit-identity: `f64::mul_add` is a single IEEE-754 fused operation
  under both codegens (libm call on baseline, `vfmadd` here) — same
  value per element. The stiffness path is gated by the crate's sf-kron
  golden hash; the material mass paths have portable-versus-dispatched bit comparisons.
  Cancellation, checked slice admission and atomic output publication use the
  identical safe mass body in both codegens.

## Blast radius

None beyond the capsule contract: if feature detection were wrong the
CPU would fault on the first FMA instruction (no silent corruption);
results cannot differ because the arithmetic is operation-for-
operation identical.
