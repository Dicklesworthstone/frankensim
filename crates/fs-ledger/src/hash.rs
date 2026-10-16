//! In-house BLAKE3 content hashing (plan §11.2, Bet 10; Decalogue P1/P9).
//!
//! The implementation lives in the dependency-light UTIL crate
//! [`fs_blake3`] (extracted verbatim from here by bead 7uq9 so exactly one
//! BLAKE3 exists in-tree); this module re-exports it unchanged, so every
//! historical `fs_ledger::hash::…` and `fs_ledger::ContentHash` path —
//! including the `ledger_001` conformance vectors — keeps working with
//! identical semantics and identical bits.
//!
//! Output is the fixed 32-byte root hash used as artifact identity
//! everywhere in the Design Ledger. Only the plain hash mode is
//! implemented; keyed hashing, key derivation, and extended (XOF) output
//! are no-claim (CONTRACT.md).
//!
//! Determinism class: pure function of the input bytes — bit-stable across
//! runs, thread counts, and ISAs.

pub use fs_blake3::{Blake3, ContentHash, hash_bytes};
