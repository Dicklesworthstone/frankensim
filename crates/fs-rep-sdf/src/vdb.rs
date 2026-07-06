//! FrankenVDB (plan §7.2): the in-house sparse hierarchical tile tree —
//! root map → 32³ internal nodes → 8³ bitmasked leaves. THE shared spatial
//! substrate: LBM lattices (§8.3), CutFEM background grids (§8.1), voxel
//! charts, and the narrow band all sit on it.
//!
//! Determinism discipline: the root and internal children are BTreeMaps
//! (sorted iteration — HashMap iteration order is banned from results by
//! the fs-exec determinism contract), and in-leaf order is linear index
//! order, so `iter_active` is a total deterministic order.

use std::collections::BTreeMap;
use std::fmt::Write as _;

const LEAF_LOG2: i32 = 3;
const LEAF_DIM: i32 = 1 << LEAF_LOG2; // 8
const LEAF_VOLUME: usize = 512;
const INTERNAL_LOG2: i32 = 5;
const INTERNAL_DIM: i32 = 1 << INTERNAL_LOG2; // 32

/// An 8³ bitmasked leaf: full value array + active mask.
struct Leaf<T> {
    mask: [u64; LEAF_VOLUME / 64],
    values: [T; LEAF_VOLUME],
}

impl<T: Copy> Leaf<T> {
    fn new(background: T) -> Self {
        Leaf {
            mask: [0; LEAF_VOLUME / 64],
            values: [background; LEAF_VOLUME],
        }
    }

    fn bit(idx: usize) -> (usize, u64) {
        (idx / 64, 1u64 << (idx % 64))
    }

    fn is_active(&self, idx: usize) -> bool {
        let (w, b) = Self::bit(idx);
        self.mask[w] & b != 0
    }

    fn set(&mut self, idx: usize, v: T) {
        let (w, b) = Self::bit(idx);
        self.mask[w] |= b;
        self.values[idx] = v;
    }

    fn deactivate(&mut self, idx: usize) {
        let (w, b) = Self::bit(idx);
        self.mask[w] &= !b;
    }

    fn active_count(&self) -> u64 {
        self.mask.iter().map(|w| u64::from(w.count_ones())).sum()
    }
}

/// One 32³ internal node: sparse children keyed by packed in-node index
/// (BTreeMap keeps memory proportional to occupancy AND iteration
/// deterministic).
struct InternalNode<T> {
    children: BTreeMap<u32, Box<Leaf<T>>>,
}

impl<T> InternalNode<T> {
    fn new() -> Self {
        InternalNode {
            children: BTreeMap::new(),
        }
    }
}

/// The sparse grid. `T` is the voxel payload (f32 for SDF/level-set use).
pub struct VdbGrid<T: Copy> {
    background: T,
    root: BTreeMap<[i32; 3], InternalNode<T>>,
}

/// Memory/topology statistics (the footprint evidence the acceptance
/// criteria ask for — measured, ledgered, deterministic).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VdbStats {
    /// Active voxels.
    pub active: u64,
    /// Allocated leaves.
    pub leaves: u64,
    /// Allocated internal nodes.
    pub internals: u64,
    /// Resident payload bytes (leaves + masks + node maps, estimated from
    /// layout — not an allocator measurement).
    pub resident_bytes: u64,
    /// Resident bytes per active voxel.
    pub bytes_per_active: f64,
}

impl VdbStats {
    /// Canonical JSON.
    #[must_use]
    pub fn to_json(&self) -> String {
        let mut s = String::with_capacity(96);
        let _ = write!(
            s,
            "{{\"active\":{},\"leaves\":{},\"internals\":{},\"resident_bytes\":{},\
             \"bytes_per_active\":{:.2}}}",
            self.active, self.leaves, self.internals, self.resident_bytes, self.bytes_per_active
        );
        s
    }
}

/// Split a global voxel coordinate into (root key, packed internal index,
/// packed leaf index). Arithmetic shifts give floor semantics for
/// negatives; low-bit masks are non-negative in two's complement.
fn split(coord: [i32; 3]) -> ([i32; 3], u32, usize) {
    let leaf = [
        coord[0] >> LEAF_LOG2,
        coord[1] >> LEAF_LOG2,
        coord[2] >> LEAF_LOG2,
    ];
    let root = [
        leaf[0] >> INTERNAL_LOG2,
        leaf[1] >> INTERNAL_LOG2,
        leaf[2] >> INTERNAL_LOG2,
    ];
    let li = [
        (leaf[0] & (INTERNAL_DIM - 1)) as u32,
        (leaf[1] & (INTERNAL_DIM - 1)) as u32,
        (leaf[2] & (INTERNAL_DIM - 1)) as u32,
    ];
    let internal_idx = (li[2] << (2 * INTERNAL_LOG2)) | (li[1] << INTERNAL_LOG2) | li[0];
    let vi = [
        (coord[0] & (LEAF_DIM - 1)) as usize,
        (coord[1] & (LEAF_DIM - 1)) as usize,
        (coord[2] & (LEAF_DIM - 1)) as usize,
    ];
    let leaf_idx = (vi[2] << (2 * LEAF_LOG2 as usize)) | (vi[1] << (LEAF_LOG2 as usize)) | vi[0];
    (root, internal_idx, leaf_idx)
}

/// Inverse of [`split`]'s leaf-local part: rebuild the global coordinate.
fn unsplit(root: [i32; 3], internal_idx: u32, leaf_idx: usize) -> [i32; 3] {
    let li = [
        (internal_idx & (INTERNAL_DIM as u32 - 1)) as i32,
        ((internal_idx >> INTERNAL_LOG2) & (INTERNAL_DIM as u32 - 1)) as i32,
        (internal_idx >> (2 * INTERNAL_LOG2)) as i32,
    ];
    let vi = [
        (leaf_idx & (LEAF_DIM as usize - 1)) as i32,
        ((leaf_idx >> LEAF_LOG2) & (LEAF_DIM as usize - 1)) as i32,
        (leaf_idx >> (2 * LEAF_LOG2)) as i32,
    ];
    [
        ((root[0] << INTERNAL_LOG2) + li[0]) * LEAF_DIM + vi[0],
        ((root[1] << INTERNAL_LOG2) + li[1]) * LEAF_DIM + vi[1],
        ((root[2] << INTERNAL_LOG2) + li[2]) * LEAF_DIM + vi[2],
    ]
}

impl<T: Copy> VdbGrid<T> {
    /// Empty grid returning `background` everywhere.
    #[must_use]
    pub fn new(background: T) -> Self {
        VdbGrid {
            background,
            root: BTreeMap::new(),
        }
    }

    /// The background value.
    #[must_use]
    pub fn background(&self) -> T {
        self.background
    }

    /// Activate a voxel with a value (O(log root) + O(log node) — the
    /// "O(1)-ish" of the acceptance criteria, verified by the scaling
    /// probe in the conformance suite).
    pub fn set(&mut self, coord: [i32; 3], value: T) {
        let (root, ii, li) = split(coord);
        self.root
            .entry(root)
            .or_insert_with(InternalNode::new)
            .children
            .entry(ii)
            .or_insert_with(|| Box::new(Leaf::new(self.background)))
            .set(li, value);
    }

    /// Read a voxel (background when inactive).
    #[must_use]
    pub fn get(&self, coord: [i32; 3]) -> T {
        let (root, ii, li) = split(coord);
        self.root
            .get(&root)
            .and_then(|n| n.children.get(&ii))
            .and_then(|leaf| leaf.is_active(li).then(|| leaf.values[li]))
            .unwrap_or(self.background)
    }

    /// True when the voxel is active.
    #[must_use]
    pub fn is_active(&self, coord: [i32; 3]) -> bool {
        let (root, ii, li) = split(coord);
        self.root
            .get(&root)
            .and_then(|n| n.children.get(&ii))
            .is_some_and(|leaf| leaf.is_active(li))
    }

    /// Deactivate a voxel (value reverts to background on read).
    pub fn deactivate(&mut self, coord: [i32; 3]) {
        let (root, ii, li) = split(coord);
        if let Some(node) = self.root.get_mut(&root)
            && let Some(leaf) = node.children.get_mut(&ii)
        {
            leaf.deactivate(li);
        }
    }

    /// Active voxel count.
    #[must_use]
    pub fn active_count(&self) -> u64 {
        self.root
            .values()
            .flat_map(|n| n.children.values())
            .map(|l| l.active_count())
            .sum()
    }

    /// Deterministic iteration over active voxels (root order → internal
    /// order → leaf linear order).
    pub fn iter_active(&self) -> impl Iterator<Item = ([i32; 3], T)> + '_ {
        self.root.iter().flat_map(|(rk, node)| {
            node.children.iter().flat_map(move |(ii, leaf)| {
                (0..LEAF_VOLUME).filter_map(move |li| {
                    leaf.is_active(li)
                        .then(|| (unsplit(*rk, *ii, li), leaf.values[li]))
                })
            })
        })
    }

    /// Dilate the active set by one voxel across the 6 face directions.
    /// New voxels take the activating neighbor's value; already-active
    /// voxels keep theirs. Deterministic (iteration order + fixed
    /// direction order).
    pub fn dilate(&mut self) {
        const DIRS: [[i32; 3]; 6] = [
            [-1, 0, 0],
            [1, 0, 0],
            [0, -1, 0],
            [0, 1, 0],
            [0, 0, -1],
            [0, 0, 1],
        ];
        let snapshot: Vec<([i32; 3], T)> = self.iter_active().collect();
        for (c, v) in snapshot {
            for d in DIRS {
                let n = [c[0] + d[0], c[1] + d[1], c[2] + d[2]];
                if !self.is_active(n) {
                    self.set(n, v);
                }
            }
        }
    }

    /// Erode the active set by one voxel: deactivate actives with any
    /// inactive face neighbor. Deterministic.
    pub fn erode(&mut self) {
        const DIRS: [[i32; 3]; 6] = [
            [-1, 0, 0],
            [1, 0, 0],
            [0, -1, 0],
            [0, 1, 0],
            [0, 0, -1],
            [0, 0, 1],
        ];
        let doomed: Vec<[i32; 3]> = self
            .iter_active()
            .filter(|&(c, _)| {
                DIRS.iter().any(|d| {
                    !self.is_active([c[0] + d[0], c[1] + d[1], c[2] + d[2]])
                })
            })
            .map(|(c, _)| c)
            .collect();
        for c in doomed {
            self.deactivate(c);
        }
    }

    /// Layout-derived memory statistics.
    #[must_use]
    pub fn memory_stats(&self) -> VdbStats {
        let leaves: u64 = self.root.values().map(|n| n.children.len() as u64).sum();
        let internals = self.root.len() as u64;
        let leaf_bytes = (size_of::<Leaf<T>>() + size_of::<u32>() + size_of::<usize>()) as u64;
        let resident = leaves * leaf_bytes + internals * 64;
        let active = self.active_count();
        VdbStats {
            active,
            leaves,
            internals,
            resident_bytes: resident,
            bytes_per_active: resident as f64 / active.max(1) as f64,
        }
    }
}
