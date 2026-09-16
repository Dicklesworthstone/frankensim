//! Uniform refinement of an existing tetrahedral mesh, with exact-topology
//! parentage for P1 fields and surface triangles. This adapts the existing
//! red-refinement kernel; it is not a second tetrahedralization algorithm.
//!
//! Vertex identity, not coordinate equality, defines an edge. In particular,
//! separately numbered coincident contact traces are never welded together.
//! Output-count limits are not a total allocator/workspace memory guarantee.

use std::collections::{BTreeMap, BTreeSet};
use fs_exec::Cx;
use fs_ivl::{Sign, orient3d};

/// Hard output-count limits checked before allocating the refined mesh.
#[derive(Debug, Clone, Copy)]
pub struct TetRefinementLimits {
    /// Maximum number of output vertices, including retained original nodes.
    pub max_vertices: usize,
    /// Maximum number of output tetrahedra.
    pub max_tetrahedra: usize,
}

/// A rejected input, exceeded resource limit, or cancelled refinement.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TetRefinementError {
    /// Invalid vertex indices, non-finite coordinates, or degenerate cells.
    InvalidMesh,
    /// The next complete refinement exceeds an output-count limit.
    OutputLimit,
    /// A midpoint or child cell is not representable without degeneration.
    UnrepresentableSplit,
    /// Field arity or finiteness does not match the parent mesh.
    InvalidField,
    /// The requested triangle does not belong to a parent tetrahedron.
    UnknownFace,
    /// Cancellation was observed at a checkpoint.
    Cancelled,
}
impl core::fmt::Display for TetRefinementError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "tet refinement: {self:?}")
    }
}
impl std::error::Error for TetRefinementError {}

/// One conforming 1-to-8 refinement. Children of cell `i` occupy `8*i..8*i+8`.
/// Original vertices keep their indices; every new vertex is an edge midpoint.
#[derive(Debug)]
pub struct TetRefinement {
    original_vertices: usize,
    positions: Vec<[f64; 3]>,
    tetrahedra: Vec<[u32; 4]>,
    midpoint_parents: Vec<[u32; 2]>,
    edges: BTreeMap<[u32; 2], u32>,
    faces: BTreeSet<[u32; 3]>,
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), TetRefinementError> {
    cx.checkpoint().map_err(|_| TetRefinementError::Cancelled)
}
fn edge(a: u32, b: u32) -> [u32; 2] { if a < b { [a,b] } else { [b,a] } }

impl TetRefinement {
    /// Refine an admitted array mesh with the existing shortest-diagonal kernel.
    ///
    /// Cancellation is checked during admission and output verification, and
    /// around the single bounded split-kernel call. This does not repair an
    /// overlapping or otherwise invalid parent mesh or improve its quality.
    ///
    /// # Errors
    /// Refuses invalid/degenerate cells, exceeded output limits, cancellation,
    /// and floating-point geometry that cannot represent the requested split.
    pub fn build(cx: &Cx<'_>, positions: &[[f64; 3]], tetrahedra: &[[u32; 4]],
        limits: TetRefinementLimits) -> Result<Self, TetRefinementError> {
        checkpoint(cx)?;
        if positions.len() < 4 || tetrahedra.is_empty() {
            return Err(TetRefinementError::InvalidMesh);
        }
        let children = tetrahedra.len().checked_mul(8)
            .filter(|&n| n <= limits.max_tetrahedra)
            .ok_or(TetRefinementError::OutputLimit)?;
        if positions.len() > limits.max_vertices || positions.len() > u32::MAX as usize {
            return Err(TetRefinementError::OutputLimit);
        }
        for (i, position) in positions.iter().enumerate() {
            if i % 512 == 0 { checkpoint(cx)?; }
            if position.iter().any(|v| !v.is_finite()) { return Err(TetRefinementError::InvalidMesh); }
        }
        let mut edges = BTreeMap::new();
        let mut midpoint_parents = Vec::new();
        let mut faces = BTreeSet::new();
        for (i, cell) in tetrahedra.iter().enumerate() {
            if i % 256 == 0 { checkpoint(cx)?; }
            let mut unique = *cell;
            unique.sort_unstable();
            if unique.windows(2).any(|p|p[0]==p[1])
                || cell.iter().any(|&v|v as usize >= positions.len()) {
                return Err(TetRefinementError::InvalidMesh);
            }
            let p = |v:u32| positions[v as usize];
            if orient3d(p(cell[0]),p(cell[1]),p(cell[2]),p(cell[3])) == Sign::Zero {
                return Err(TetRefinementError::InvalidMesh);
            }
            // This is lineage enumeration, in the kernel's edge-visit order.
            // The returned coordinates are checked below, so a future change
            // to that order cannot silently transfer fields to the wrong node.
            for (a,b) in [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)] {
                let key = edge(cell[a],cell[b]);
                if !edges.contains_key(&key) {
                    let index = positions.len().checked_add(midpoint_parents.len())
                        .filter(|&n| n < limits.max_vertices && n < u32::MAX as usize)
                        .ok_or(TetRefinementError::OutputLimit)?;
                    edges.insert(key,index as u32);
                    midpoint_parents.push(key);
                }
            }
            for omitted in 0..4 {
                let mut face = [0;3];
                let mut slot = 0;
                for (j,&v) in cell.iter().enumerate() {
                    if j != omitted { face[slot]=v; slot+=1; }
                }
                face.sort_unstable();
                faces.insert(face);
            }
        }
        checkpoint(cx)?;
        let split = crate::uniform::split_uniform(positions,tetrahedra,&[]);
        checkpoint(cx)?;
        if split.tets.len() != children || split.positions.len() != positions.len()+midpoint_parents.len() {
            return Err(TetRefinementError::UnrepresentableSplit);
        }
        for (i,&[a,b]) in midpoint_parents.iter().enumerate() {
            if i % 256 == 0 { checkpoint(cx)?; }
            let midpoint = positions[a as usize].map_with(positions[b as usize], |x,y| f64::midpoint(x,y));
            let actual = split.positions[positions.len()+i];
            if actual != midpoint || actual == positions[a as usize] || actual == positions[b as usize] {
                return Err(TetRefinementError::UnrepresentableSplit);
            }
        }
        for (i, cell) in split.tets.iter().enumerate() {
            if i % 256 == 0 { checkpoint(cx)?; }
            let p = |v:u32| split.positions[v as usize];
            if orient3d(p(cell[0]),p(cell[1]),p(cell[2]),p(cell[3])) == Sign::Zero {
                return Err(TetRefinementError::UnrepresentableSplit);
            }
        }
        Ok(Self { original_vertices: positions.len(), positions: split.positions,
            tetrahedra: split.tets, midpoint_parents, edges, faces })
    }

    /// Refined positions, with original nodes first and unchanged.
    #[must_use]
    pub fn positions(&self) -> &[[f64; 3]] { &self.positions }
    /// Refined tetrahedra, grouped eight at a time by parent cell.
    #[must_use]
    pub fn tetrahedra(&self) -> &[[u32; 4]] { &self.tetrahedra }
    /// Parent edge of each new vertex, in appended-vertex order.
    #[must_use]
    pub fn midpoint_parents(&self) -> &[[u32; 2]] { &self.midpoint_parents }

    /// Inject the parent P1 scalar function into the refined P1 space.
    /// This preserves its physical support and integrals in exact arithmetic;
    /// it does not renormalize a nodal source or reinterpret its footprint.
    ///
    /// # Errors
    /// Refuses a non-finite/wrong-length parent field or cancellation.
    pub fn prolongate(&self, cx: &Cx<'_>, field: &[f64]) -> Result<Vec<f64>, TetRefinementError> {
        checkpoint(cx)?;
        if field.len()!=self.original_vertices || field.iter().any(|v|!v.is_finite()) {
            return Err(TetRefinementError::InvalidField);
        }
        let mut result=Vec::with_capacity(self.positions.len());
        result.extend_from_slice(field);
        for (i,&[a,b]) in self.midpoint_parents.iter().enumerate() {
            if i%512==0 {checkpoint(cx)?;}
            result.push(f64::midpoint(field[a as usize],field[b as usize]));
        }
        Ok(result)
    }

    /// Four child triangles in the supplied local vertex order. Matching
    /// contact sides should first be ordered by their geometric correspondence.
    ///
    /// # Errors
    /// Refuses a triangle not present in the parent tetrahedral complex.
    pub fn face_children(&self, face:[u32;3]) -> Result<[[u32;3];4],TetRefinementError> {
        let mut key=face; key.sort_unstable();
        if !self.faces.contains(&key) {return Err(TetRefinementError::UnknownFace);}
        let [a,b,c]=face;
        let ab=self.edges[&edge(a,b)]; let bc=self.edges[&edge(b,c)]; let ca=self.edges[&edge(c,a)];
        Ok([[a,ab,ca],[ab,b,bc],[ca,bc,c],[ab,bc,ca]])
    }
}
