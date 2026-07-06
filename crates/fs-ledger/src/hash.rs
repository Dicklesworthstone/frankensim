//! In-house BLAKE3 content hashing (plan §11.2, Bet 10; Decalogue P1/P9).
//!
//! A pure safe-Rust implementation of the BLAKE3 hash function following the
//! reference algorithm: 1024-byte chunks compressed in 64-byte blocks by a
//! 7-round ChaCha-derived compression function, merged into a binary tree of
//! parent nodes. Output is the fixed 32-byte root hash used as artifact
//! identity everywhere in the Design Ledger.
//!
//! Verified against oracle-generated official-pattern test vectors in the
//! conformance suite (`ledger_001`), including multi-chunk and multi-level
//! tree inputs. Only the plain hash mode is implemented; keyed hashing,
//! key derivation, and extended (XOF) output are no-claim (CONTRACT.md).
//!
//! Determinism class: pure function of the input bytes — bit-stable across
//! runs, thread counts, and ISAs.

use core::fmt;

/// Initialization vector (the SHA-256 IV words, per the BLAKE3 spec).
const IV: [u32; 8] = [
    0x6A09_E667,
    0xBB67_AE85,
    0x3C6E_F372,
    0xA54F_F53A,
    0x510E_527F,
    0x9B05_688C,
    0x1F83_D9AB,
    0x5BE0_CD19,
];

/// The fixed message-word permutation applied between rounds.
const MSG_PERMUTATION: [usize; 16] = [2, 6, 3, 10, 7, 0, 4, 13, 1, 11, 12, 5, 9, 14, 15, 8];

const BLOCK_LEN: usize = 64;
const CHUNK_LEN: usize = 1024;

const CHUNK_START: u32 = 1 << 0;
const CHUNK_END: u32 = 1 << 1;
const PARENT: u32 = 1 << 2;
const ROOT: u32 = 1 << 3;

/// The quarter-round mixing function.
fn g(state: &mut [u32; 16], a: usize, b: usize, c: usize, d: usize, mx: u32, my: u32) {
    state[a] = state[a].wrapping_add(state[b]).wrapping_add(mx);
    state[d] = (state[d] ^ state[a]).rotate_right(16);
    state[c] = state[c].wrapping_add(state[d]);
    state[b] = (state[b] ^ state[c]).rotate_right(12);
    state[a] = state[a].wrapping_add(state[b]).wrapping_add(my);
    state[d] = (state[d] ^ state[a]).rotate_right(8);
    state[c] = state[c].wrapping_add(state[d]);
    state[b] = (state[b] ^ state[c]).rotate_right(7);
}

/// One full round: four column mixes then four diagonal mixes.
fn round(state: &mut [u32; 16], m: &[u32; 16]) {
    g(state, 0, 4, 8, 12, m[0], m[1]);
    g(state, 1, 5, 9, 13, m[2], m[3]);
    g(state, 2, 6, 10, 14, m[4], m[5]);
    g(state, 3, 7, 11, 15, m[6], m[7]);
    g(state, 0, 5, 10, 15, m[8], m[9]);
    g(state, 1, 6, 11, 12, m[10], m[11]);
    g(state, 2, 7, 8, 13, m[12], m[13]);
    g(state, 3, 4, 9, 14, m[14], m[15]);
}

fn permute(m: &mut [u32; 16]) {
    let mut permuted = [0u32; 16];
    for (dst, &src) in permuted.iter_mut().zip(MSG_PERMUTATION.iter()) {
        *dst = m[src];
    }
    *m = permuted;
}

/// The BLAKE3 compression function (7 rounds + feed-forward).
fn compress(
    chaining_value: &[u32; 8],
    block_words: &[u32; 16],
    counter: u64,
    block_len: u32,
    flags: u32,
) -> [u32; 16] {
    let mut state = [
        chaining_value[0],
        chaining_value[1],
        chaining_value[2],
        chaining_value[3],
        chaining_value[4],
        chaining_value[5],
        chaining_value[6],
        chaining_value[7],
        IV[0],
        IV[1],
        IV[2],
        IV[3],
        counter as u32,
        (counter >> 32) as u32,
        block_len,
        flags,
    ];
    let mut block = *block_words;
    round(&mut state, &block); // round 1
    permute(&mut block);
    round(&mut state, &block); // round 2
    permute(&mut block);
    round(&mut state, &block); // round 3
    permute(&mut block);
    round(&mut state, &block); // round 4
    permute(&mut block);
    round(&mut state, &block); // round 5
    permute(&mut block);
    round(&mut state, &block); // round 6
    permute(&mut block);
    round(&mut state, &block); // round 7
    for i in 0..8 {
        state[i] ^= state[i + 8];
        state[i + 8] ^= chaining_value[i];
    }
    state
}

/// Little-endian words from a (zero-padded) 64-byte block.
fn words_from_block(block: &[u8; BLOCK_LEN]) -> [u32; 16] {
    let mut words = [0u32; 16];
    let (quads, _) = block.as_chunks::<4>();
    for (word, bytes) in words.iter_mut().zip(quads) {
        *word = u32::from_le_bytes(*bytes);
    }
    words
}

fn first_8_words(compression_output: [u32; 16]) -> [u32; 8] {
    let mut cv = [0u32; 8];
    cv.copy_from_slice(&compression_output[..8]);
    cv
}

/// A pending compression whose ROOT flag decision is deferred until we know
/// whether it is the tree root (the reference implementation's `Output`).
struct Output {
    input_chaining_value: [u32; 8],
    block_words: [u32; 16],
    counter: u64,
    block_len: u32,
    flags: u32,
}

impl Output {
    fn chaining_value(&self) -> [u32; 8] {
        first_8_words(compress(
            &self.input_chaining_value,
            &self.block_words,
            self.counter,
            self.block_len,
            self.flags,
        ))
    }

    fn root_hash(&self) -> [u8; 32] {
        // Fixed 32-byte output: the first output block (counter 0) suffices.
        let words = compress(
            &self.input_chaining_value,
            &self.block_words,
            0,
            self.block_len,
            self.flags | ROOT,
        );
        let mut bytes = [0u8; 32];
        let (quads, _) = bytes.as_chunks_mut::<4>();
        for (chunk, word) in quads.iter_mut().zip(words.iter()) {
            *chunk = word.to_le_bytes();
        }
        bytes
    }
}

fn parent_output(left_child_cv: [u32; 8], right_child_cv: [u32; 8]) -> Output {
    let mut block_words = [0u32; 16];
    block_words[..8].copy_from_slice(&left_child_cv);
    block_words[8..].copy_from_slice(&right_child_cv);
    Output {
        input_chaining_value: IV,
        block_words,
        counter: 0, // parent nodes always use counter 0
        block_len: BLOCK_LEN as u32,
        flags: PARENT,
    }
}

/// Incremental state for one 1024-byte chunk.
struct ChunkState {
    chaining_value: [u32; 8],
    chunk_counter: u64,
    block: [u8; BLOCK_LEN],
    block_len: u8,
    blocks_compressed: u8,
}

impl ChunkState {
    fn new(chunk_counter: u64) -> Self {
        ChunkState {
            chaining_value: IV,
            chunk_counter,
            block: [0; BLOCK_LEN],
            block_len: 0,
            blocks_compressed: 0,
        }
    }

    fn len(&self) -> usize {
        BLOCK_LEN * self.blocks_compressed as usize + self.block_len as usize
    }

    fn start_flag(&self) -> u32 {
        if self.blocks_compressed == 0 {
            CHUNK_START
        } else {
            0
        }
    }

    fn update(&mut self, mut input: &[u8]) {
        while !input.is_empty() {
            // Compress a full buffered block only once MORE input arrives, so
            // the final block (which needs CHUNK_END) is always identifiable.
            if self.block_len as usize == BLOCK_LEN {
                let block_words = words_from_block(&self.block);
                self.chaining_value = first_8_words(compress(
                    &self.chaining_value,
                    &block_words,
                    self.chunk_counter,
                    BLOCK_LEN as u32,
                    self.start_flag(),
                ));
                self.blocks_compressed += 1;
                self.block = [0; BLOCK_LEN];
                self.block_len = 0;
            }
            let want = BLOCK_LEN - self.block_len as usize;
            let take = want.min(input.len());
            self.block[self.block_len as usize..self.block_len as usize + take]
                .copy_from_slice(&input[..take]);
            self.block_len += take as u8;
            input = &input[take..];
        }
    }

    fn output(&self) -> Output {
        Output {
            input_chaining_value: self.chaining_value,
            block_words: words_from_block(&self.block),
            counter: self.chunk_counter,
            block_len: u32::from(self.block_len),
            flags: self.start_flag() | CHUNK_END,
        }
    }
}

/// Maximum tree depth: 2^54 chunks covers the largest input BLAKE3 defines
/// (and vastly exceeds fsqlite's 1 GB value limit).
const MAX_DEPTH: usize = 54;

/// Streaming BLAKE3 hasher (plain hash mode, 32-byte output).
///
/// Feed bytes with [`Blake3::update`] in any split pattern; the digest is
/// identical to hashing the concatenation in one call (property-tested).
pub struct Blake3 {
    chunk_state: ChunkState,
    cv_stack: [[u32; 8]; MAX_DEPTH],
    cv_stack_len: u8,
}

impl Default for Blake3 {
    fn default() -> Self {
        Blake3::new()
    }
}

impl Blake3 {
    /// A fresh hasher.
    #[must_use]
    pub fn new() -> Self {
        Blake3 {
            chunk_state: ChunkState::new(0),
            cv_stack: [[0u32; 8]; MAX_DEPTH],
            cv_stack_len: 0,
        }
    }

    fn push_stack(&mut self, cv: [u32; 8]) {
        self.cv_stack[self.cv_stack_len as usize] = cv;
        self.cv_stack_len += 1;
    }

    fn pop_stack(&mut self) -> [u32; 8] {
        self.cv_stack_len -= 1;
        self.cv_stack[self.cv_stack_len as usize]
    }

    /// Fold completed subtrees: a chunk count with k trailing zero bits has
    /// k completed subtrees to merge before its chaining value is pushed.
    fn add_chunk_chaining_value(&mut self, mut new_cv: [u32; 8], mut total_chunks: u64) {
        while total_chunks & 1 == 0 {
            new_cv = parent_output(self.pop_stack(), new_cv).chaining_value();
            total_chunks >>= 1;
        }
        self.push_stack(new_cv);
    }

    /// Absorb `input`.
    pub fn update(&mut self, mut input: &[u8]) {
        while !input.is_empty() {
            // Close a full chunk only once more input arrives, so the final
            // chunk (whose output may become the root) stays open.
            if self.chunk_state.len() == CHUNK_LEN {
                let chunk_cv = self.chunk_state.output().chaining_value();
                let total_chunks = self.chunk_state.chunk_counter + 1;
                self.add_chunk_chaining_value(chunk_cv, total_chunks);
                self.chunk_state = ChunkState::new(total_chunks);
            }
            let want = CHUNK_LEN - self.chunk_state.len();
            let take = want.min(input.len());
            self.chunk_state.update(&input[..take]);
            input = &input[take..];
        }
    }

    /// The 32-byte root hash of everything absorbed so far.
    #[must_use]
    pub fn finalize(&self) -> ContentHash {
        let mut output = self.chunk_state.output();
        let mut remaining = self.cv_stack_len as usize;
        while remaining > 0 {
            remaining -= 1;
            output = parent_output(self.cv_stack[remaining], output.chaining_value());
        }
        ContentHash(output.root_hash())
    }
}

/// Hash a byte slice in one call.
#[must_use]
pub fn hash_bytes(bytes: &[u8]) -> ContentHash {
    let mut hasher = Blake3::new();
    hasher.update(bytes);
    hasher.finalize()
}

/// A 32-byte BLAKE3 content hash: the identity of every ledger artifact.
///
/// Stored as a 32-byte BLOB primary key; rendered as 64 lowercase hex chars.
#[derive(Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct ContentHash(pub [u8; 32]);

impl ContentHash {
    /// The raw 32 bytes (the database key).
    #[must_use]
    pub fn as_bytes(&self) -> &[u8; 32] {
        &self.0
    }

    /// Lowercase hex rendering (64 chars).
    #[must_use]
    pub fn to_hex(&self) -> String {
        use core::fmt::Write as _;
        let mut s = String::with_capacity(64);
        for b in &self.0 {
            let _ = write!(s, "{b:02x}");
        }
        s
    }

    /// Parse 64 hex chars (either case). Returns `None` on any malformation.
    #[must_use]
    pub fn from_hex(s: &str) -> Option<ContentHash> {
        let raw = s.as_bytes();
        if raw.len() != 64 {
            return None;
        }
        let mut out = [0u8; 32];
        let (pairs, _) = raw.as_chunks::<2>();
        for (dst, pair) in out.iter_mut().zip(pairs) {
            let hi = (pair[0] as char).to_digit(16)?;
            let lo = (pair[1] as char).to_digit(16)?;
            *dst = (hi * 16 + lo) as u8;
        }
        Some(ContentHash(out))
    }

    /// Construct from a raw 32-byte slice (e.g. a database BLOB column).
    #[must_use]
    pub fn from_slice(bytes: &[u8]) -> Option<ContentHash> {
        let arr: [u8; 32] = bytes.try_into().ok()?;
        Some(ContentHash(arr))
    }
}

impl fmt::Display for ContentHash {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.to_hex())
    }
}

impl fmt::Debug for ContentHash {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "ContentHash({})", self.to_hex())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The official BLAKE3 empty-input digest (spec test vector).
    #[test]
    fn empty_input_matches_spec() {
        assert_eq!(
            hash_bytes(b"").to_hex(),
            "af1349b9f5f9a1a6a0404dea36dcc9499bcb25c9adc112b7cc9a93cae41f3262"
        );
    }

    #[test]
    fn abc_matches_oracle() {
        assert_eq!(
            hash_bytes(b"abc").to_hex(),
            "6437b3ac38465133ffb63b75273a8db548c558465d79db03fd359c6cd5bd9d85"
        );
    }

    #[test]
    fn hex_round_trip_and_rejection() {
        let h = hash_bytes(b"round trip");
        assert_eq!(ContentHash::from_hex(&h.to_hex()), Some(h));
        assert_eq!(ContentHash::from_hex("zz"), None);
        assert_eq!(ContentHash::from_hex(&"0".repeat(63)), None);
        assert_eq!(ContentHash::from_hex(&"g".repeat(64)), None);
        assert_eq!(ContentHash::from_slice(&[0u8; 31]), None);
        assert!(ContentHash::from_slice(&[0u8; 32]).is_some());
    }

    #[test]
    fn streaming_equals_one_shot_across_block_and_chunk_edges() {
        // Cover splits straddling block (64) and chunk (1024) boundaries.
        let data: Vec<u8> = (0..4099u32).map(|i| (i % 251) as u8).collect();
        let whole = hash_bytes(&data);
        for split in [0usize, 1, 63, 64, 65, 1023, 1024, 1025, 2048, 4098] {
            let mut h = Blake3::new();
            h.update(&data[..split]);
            h.update(&data[split..]);
            assert_eq!(h.finalize(), whole, "split at {split}");
        }
    }
}
