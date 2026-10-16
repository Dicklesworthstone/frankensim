//! Raw byte access for the opt-in reclaimed-chunk poison detector.
//!
//! This capsule is deliberately separate from the nearly-full bump-pointer
//! capsule. Its crate-private methods are exposed only through `ArenaPool`'s
//! safe, lock-serialized diagnostic facade.
#![allow(unsafe_code)]
// registered capsule — see SAFETY.md beside this file

use super::Chunk;

impl Chunk {
    /// Overwrite the whole reclaimed block with the seed-derived poison byte.
    pub(crate) fn poison_for_reclaim(&mut self, seed: u64) {
        // SAFETY: `Chunk` exclusively owns a live allocation of `self.len()`
        // bytes. The caller either extracted it through `&mut Arena` after all
        // arena references died, or is rolling back a pending acquisition that
        // was never installed, so no reference into the block is alive.
        unsafe { core::ptr::write_bytes(self.base.as_ptr(), poison_byte(seed), self.len()) }
    }

    /// Return the first byte that differs from the expected reclaimed poison.
    pub(crate) fn reclaimed_poison_mismatch(&self, seed: u64) -> Option<(usize, u8, u8)> {
        let expected = poison_byte(seed);
        // SAFETY: the block is live for `self.len()` bytes and is read while
        // the free-list mutex prevents reuse or deallocation.
        let bytes = unsafe { core::slice::from_raw_parts(self.base.as_ptr(), self.len()) };
        bytes
            .iter()
            .position(|&actual| actual != expected)
            .map(|offset| (offset, expected, bytes[offset]))
    }

    /// Deterministically alter one poisoned byte for the safe G4 fault hook.
    pub(crate) fn inject_reclaimed_corruption(&mut self, seed: u64) -> (usize, u8, u8) {
        debug_assert!(self.len() > 0, "allocated chunks are never empty");
        let expected = poison_byte(seed);
        let actual = expected ^ 0xff;
        let offset = corruption_offset(seed, self.len());
        // SAFETY: the chosen offset is in-bounds, the block is live and
        // exclusively free-listed, and the free-list mutex excludes reuse.
        unsafe { self.base.as_ptr().add(offset).write(actual) };
        (offset, expected, actual)
    }
}

fn poison_byte(seed: u64) -> u8 {
    let mixed = mix(seed ^ 0xA110_C0DE_5EED_0004);
    (mixed ^ (mixed >> 8) ^ (mixed >> 16) ^ (mixed >> 24)) as u8
}

fn corruption_offset(seed: u64, len: usize) -> usize {
    (mix(seed ^ len as u64) % len as u64) as usize
}

fn mix(mut value: u64) -> u64 {
    value ^= value >> 30;
    value = value.wrapping_mul(0xbf58_476d_1ce4_e5b9);
    value ^= value >> 27;
    value = value.wrapping_mul(0x94d0_49bb_1331_11eb);
    value ^ (value >> 31)
}
