//! fs-soa — structure-of-arrays runtime (plan §5.3). Layer: L0
//! SUBSTRATE.
//!
//! The machinery `#[derive(Soa)]` targets: per-field growable buffers
//! aligned to [`SOA_ALIGN`] WITHOUT unsafe (over-allocate + whole-
//! element `align_offset`, the pattern proven by fs-substrate's
//! `AlignedBuf`), AoS gather/scatter, chunked SIMD-friendly slice
//! access with explicit masked-tail handling, zero-copy strided view
//! descriptors for the FrankenNumpy membrane (§12), and chunk-quantum
//! grouping as the tile-identity hook.
//!
//! "SoA everywhere hot" is a memory-discipline pillar: batched small
//! dense LA and LBM lattices need SIMD lanes running ACROSS elements,
//! which is only natural in this layout.

pub use fs_soa_derive::Soa;

/// Target byte alignment for every field buffer: 128 unconditionally
/// (superset of Apple 128B and x86 64B lines; matches
/// `fs_alloc::ALLOC_ALIGN`).
pub const SOA_ALIGN: usize = 128;

/// Default chunk quantum for tile-identity grouping: 512 elements =
/// the E8 tile volume (8³) from fs-substrate's tile doctrine, so SoA
/// chunk indices can serve as executor tile identities for free.
pub const DEFAULT_CHUNK_QUANTUM: usize = 512;

/// Slack elements prepended to each allocation so a 128-byte-aligned
/// start is reachable by whole-element steps: the smallest solution k
/// of (base + k·size) ≡ 0 (mod 128) is < 128 whenever one exists.
const PAD_ELEMS: usize = 128;

/// One field's growable aligned buffer. NO unsafe: the backing `Vec`
/// is over-allocated and the payload starts at the first 128-byte-
/// aligned whole-element offset. Slots outside the payload hold copies
/// of pushed values and are never exposed.
#[derive(Debug, Clone)]
pub struct FieldBuf<T: Copy> {
    data: Vec<T>,
    offset: usize,
    len: usize,
    hint: usize,
}

impl<T: Copy> Default for FieldBuf<T> {
    fn default() -> FieldBuf<T> {
        FieldBuf::new()
    }
}

impl<T: Copy> FieldBuf<T> {
    /// Empty buffer (no allocation until the first push).
    #[must_use]
    pub const fn new() -> FieldBuf<T> {
        FieldBuf {
            data: Vec::new(),
            offset: 0,
            len: 0,
            hint: 0,
        }
    }

    /// Empty buffer that will allocate for at least `cap` elements on
    /// first push.
    #[must_use]
    pub const fn with_capacity(cap: usize) -> FieldBuf<T> {
        FieldBuf {
            data: Vec::new(),
            offset: 0,
            len: 0,
            hint: cap,
        }
    }

    /// Elements stored.
    #[must_use]
    pub const fn len(&self) -> usize {
        self.len
    }

    /// True when no elements are stored.
    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Elements storable without reallocation.
    #[must_use]
    pub fn capacity(&self) -> usize {
        if self.data.is_empty() {
            0
        } else {
            self.data.len() - self.offset
        }
    }

    /// The payload as a slice.
    #[must_use]
    pub fn as_slice(&self) -> &[T] {
        &self.data[self.offset..self.offset + self.len]
    }

    /// The payload as a mutable slice.
    pub fn as_mut_slice(&mut self) -> &mut [T] {
        &mut self.data[self.offset..self.offset + self.len]
    }

    /// Drop all elements (keeps the allocation).
    pub fn clear(&mut self) {
        self.len = 0;
    }

    /// Ensure capacity for `additional` more elements.
    pub fn reserve(&mut self, additional: usize) {
        if self.data.is_empty() {
            self.hint = self.hint.max(self.len + additional);
        } else if self.len + additional > self.capacity() {
            self.grow(self.len + additional);
        }
    }

    /// Append one element.
    pub fn push(&mut self, value: T) {
        if self.data.is_empty() {
            let cap = self.hint.max(8);
            let data = vec![value; cap + PAD_ELEMS];
            self.offset = aligned_offset(data.as_ptr());
            self.data = data;
            self.len = 1;
            return;
        }
        if self.len == self.capacity() {
            self.grow(self.len * 2);
        }
        self.data[self.offset + self.len] = value;
        self.len += 1;
    }

    fn grow(&mut self, new_cap: usize) {
        let seed = self.data[self.offset];
        let mut data = vec![seed; new_cap.max(8) + PAD_ELEMS];
        let offset = aligned_offset(data.as_ptr());
        data[offset..offset + self.len]
            .copy_from_slice(&self.data[self.offset..self.offset + self.len]);
        self.data = data;
        self.offset = offset;
    }

    /// Zero-copy view descriptor for the payload (the FrankenNumpy
    /// membrane shape: address + dense stride; `addr` is 0 for an
    /// unallocated buffer).
    #[must_use]
    pub fn view(&self, name: &str) -> RawView {
        let addr = if self.data.is_empty() {
            0
        } else {
            self.as_slice().as_ptr().addr()
        };
        RawView {
            name: name.to_string(),
            addr,
            len: self.len,
            elem_bytes: size_of::<T>(),
            stride_bytes: size_of::<T>(),
            achieved_align: achieved_align(addr),
            dtype: std::any::type_name::<T>(),
        }
    }
}

/// First whole-element offset reaching [`SOA_ALIGN`], degrading to the
/// best reachable power-of-two alignment when 128 is unreachable by
/// element-sized steps (possible for exotic element sizes; the view
/// descriptor always reports what was ACHIEVED).
fn aligned_offset<T>(base: *const T) -> usize {
    let mut align = SOA_ALIGN;
    while align > 1 {
        let off = base.align_offset(align);
        if off != usize::MAX && off < PAD_ELEMS {
            return off;
        }
        align /= 2;
    }
    0
}

const fn achieved_align(addr: usize) -> usize {
    if addr == 0 {
        SOA_ALIGN // unallocated: alignment is vacuous, report the target
    } else {
        1 << addr.trailing_zeros()
    }
}

/// Zero-copy strided view descriptor of one leaf field buffer — the
/// §12 membrane contract's shape (pointer address, element count,
/// element/stride bytes, achieved alignment, element type name).
/// Live FrankenNumpy wiring is deliberately NOT here (see CONTRACT
/// no-claims); this is the stable descriptor it will consume.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RawView {
    /// Dotted field path ("inner.pos" for nested containers).
    pub name: String,
    /// Payload start address (0 when unallocated).
    pub addr: usize,
    /// Element count.
    pub len: usize,
    /// Bytes per element.
    pub elem_bytes: usize,
    /// Bytes between consecutive elements (dense: == `elem_bytes`).
    pub stride_bytes: usize,
    /// Largest power of two dividing `addr` (capped conceptually by
    /// the allocation; [`SOA_ALIGN`] when unallocated).
    pub achieved_align: usize,
    /// Element type name (auditability, not ABI).
    pub dtype: &'static str,
}

impl RawView {
    /// Address-free JSON description (stable across runs — addresses
    /// are excluded on purpose so logs and goldens stay deterministic).
    #[must_use]
    pub fn descr(&self) -> String {
        format!(
            "{{\"field\":\"{}\",\"len\":{},\"elem_bytes\":{},\"stride_bytes\":{},\"dtype\":\"{}\"}}",
            self.name, self.len, self.elem_bytes, self.stride_bytes, self.dtype
        )
    }
}

/// Types with a generated SoA container (implemented by
/// `#[derive(Soa)]`, consumed by `#[soa(nested)]` fields).
pub trait SoaAble: Sized {
    /// The generated container type.
    type Soa: SoaContainer<Self>;
}

/// The container operations `#[derive(Soa)]` generates, as a trait so
/// nested containers compose (an outer container's field of type
/// `Inner` stores an `InnerSoa` and drives it through this interface).
pub trait SoaContainer<T> {
    /// Empty container.
    fn c_new() -> Self;
    /// Empty container with a capacity hint.
    fn c_with_capacity(cap: usize) -> Self;
    /// Elements stored.
    fn c_len(&self) -> usize;
    /// Append one value (scattered across field buffers).
    fn c_push(&mut self, value: T);
    /// Gather element `i` back into a value.
    fn c_get(&self, i: usize) -> T;
    /// Scatter `value` into slot `i`.
    fn c_set(&mut self, i: usize, value: T);
    /// Drop all elements (keep allocations).
    fn c_clear(&mut self);
    /// Ensure room for `additional` more elements.
    fn c_reserve(&mut self, additional: usize);
    /// Append this container's leaf views under `prefix`.
    fn c_views(&self, prefix: &str, out: &mut Vec<RawView>);
    /// Append this container's leaf layout lines under `prefix`.
    fn c_layout(prefix: &str, out: &mut Vec<String>);
}

/// Join a view-path prefix and a field name with '.'.
#[must_use]
pub fn view_name(prefix: &str, field: &str) -> String {
    if prefix.is_empty() {
        field.to_string()
    } else {
        format!("{prefix}.{field}")
    }
}

/// Address-free layout line for one leaf field of type `T` (used by
/// generated `layout_descr`; logged in tests for auditability).
#[must_use]
pub fn leaf_layout<T>(name: &str) -> String {
    format!(
        "{{\"field\":\"{name}\",\"elem_bytes\":{},\"elem_align\":{},\"dtype\":\"{}\"}}",
        size_of::<T>(),
        align_of::<T>(),
        std::any::type_name::<T>()
    )
}

/// Split a field slice into SIMD-width chunks plus the masked tail:
/// the explicit (full chunks, remainder) shape kernels want. Panics if
/// `width == 0`.
pub fn chunks_with_tail<T>(s: &[T], width: usize) -> (std::slice::ChunksExact<'_, T>, &[T]) {
    let it = s.chunks_exact(width);
    let tail = it.remainder();
    (it, tail)
}

/// Mutable variant of [`chunks_with_tail`]. The tail is returned by
/// the iterator's `into_remainder` after iteration:
/// `let mut it = chunks_with_tail_mut(s, w); for c in &mut it { … }
/// let tail = it.into_remainder();`.
pub fn chunks_with_tail_mut<T>(s: &mut [T], width: usize) -> std::slice::ChunksExactMut<'_, T> {
    s.chunks_exact_mut(width)
}

/// Number of chunk-quantum groups covering `len` elements (the tile-
/// identity hook: group index = element index / quantum, stable under
/// growth). Panics if `quantum == 0`.
#[must_use]
pub const fn chunk_count(len: usize, quantum: usize) -> usize {
    len.div_ceil(quantum)
}

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
