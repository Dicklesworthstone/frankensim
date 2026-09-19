//! Bounded, immutable dyadic refinement for the 3-D CutFEM background.
//!
//! Leaves cover the whole design box. Refinement closes 2:1 balance across
//! faces, edges AND vertices, without coarsening or rebuilding unchanged leaves.
//! Quarter-lattice queries find neighbors through bounded ancestor lookups;
//! there is no all-pairs face search. Physical geometry is supplied separately.
use std::collections::{BTreeMap, BTreeSet};
use std::ops::ControlFlow;
use crate::HexCell;

/// A dyadic box: coordinates are integer cell indices at `level`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct Octant3 { level: u8, index: [u32; 3] }
impl Octant3 {
    /// Refinement level; zero is the entire background.
    #[must_use] pub const fn level(self) -> u8 { self.level }
    /// Cell indices at this level.
    #[must_use] pub const fn index(self) -> [u32; 3] { self.index }
    fn children(self) -> [Self; 8] {
        std::array::from_fn(|child| Self { level: self.level + 1,
            index: std::array::from_fn(|a| 2*self.index[a] + ((child >> a)&1) as u32) })
    }
}
/// Exact integer vertex coordinate on the tree's finest permitted lattice.
pub type OctreeNode3 = [u32; 3];
/// A positive-area shared face, emitted once in positive-axis orientation.
#[derive(Debug, Clone, Copy)]
pub struct OctreeFace3 {
    /// Lower-coordinate cell along `axis`.
    pub lower: Octant3,
    /// Higher-coordinate cell along `axis`.
    pub upper: Octant3,
    /// Cartesian axis normal to the face.
    pub axis: usize,
    /// Finest-lattice patch bounds; the normal coordinates coincide.
    pub lo: OctreeNode3,
    /// Upper patch corner.
    pub hi: OctreeNode3,
}
/// A refusal never changes the input tree.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum OctreeError3 { Invalid(&'static str), LeafBudget, LevelBudget, Cancelled }
impl std::fmt::Display for OctreeError3 {
    fn fmt(&self,f:&mut std::fmt::Formatter<'_>)->std::fmt::Result { write!(f,"3-D octree refused: {self:?}") }
}
impl std::error::Error for OctreeError3 {}
fn poll(c:&mut impl FnMut()->ControlFlow<()>)->Result<(),OctreeError3> {
    if c().is_break() { Err(OctreeError3::Cancelled) } else { Ok(()) }
}
/// Complete, strongly 2:1-balanced dyadic partition with retained resource caps.
#[derive(Debug, Clone)]
pub struct Octree3 { leaves: BTreeSet<Octant3>, max_level: u8, max_leaves: usize }
impl Octree3 {
    /// Build a uniform partition. At most 20 levels are supported; admission
    /// checks the exact leaf count before allocating it.
    pub fn uniform(level:u8,max_level:u8,max_leaves:usize)->Result<Self,OctreeError3> {
        if level>max_level || max_level>20 { return Err(OctreeError3::Invalid("invalid levels")); }
        let count=1usize.checked_shl(3*u32::from(level)).ok_or(OctreeError3::LeafBudget)?;
        if count>max_leaves { return Err(OctreeError3::LeafBudget); }
        let side=1u32<<level;
        let mut leaves=BTreeSet::new();
        for z in 0..side { for y in 0..side { for x in 0..side {
            leaves.insert(Octant3 {level,index:[x,y,z]});
        } } }
        Ok(Self {leaves,max_level,max_leaves})
    }
    /// Stable level/index order, independent of marking order.
    #[must_use] pub fn leaves(&self)->&BTreeSet<Octant3> { &self.leaves }
    /// Finest-lattice extent in each coordinate.
    #[must_use] pub const fn extent(&self)->u32 { 1u32<<self.max_level }
    /// Lower corner and edge length in finest-lattice coordinates.
    #[must_use] pub fn lattice_box(&self,c:Octant3)->(OctreeNode3,u32) {
        assert!(c.level<=self.max_level,"octant exceeds this tree's lattice");
        let size=1u32<<(self.max_level-c.level);
        (c.index.map(|i|i*size),size)
    }
    /// Eight vertices, with x/y/z chosen by the corner bits.
    #[must_use] pub fn corners(&self,c:Octant3)->[OctreeNode3;8] {
        let (lo,s)=self.lattice_box(c);
        std::array::from_fn(|v|std::array::from_fn(|a|lo[a]+s*((v>>a)&1) as u32))
    }
    /// Map an exact lattice point to a caller's physical design box.
    #[must_use] pub fn position(&self,n:OctreeNode3,domain:HexCell)->[f64;3] {
        let (lo,hi)=(domain.lo(),domain.hi());
        std::array::from_fn(|a| if n[a]==0 {lo[a]} else if n[a]==self.extent() {hi[a]}
            else {lo[a]+(hi[a]-lo[a])*(f64::from(n[a])/f64::from(self.extent()))})
    }
    /// Physical bounds, refusing collapsed floating-point coordinates.
    pub fn bounds(&self,c:Octant3,domain:HexCell)->Result<HexCell,OctreeError3> {
        let v=self.corners(c);
        HexCell::try_new(self.position(v[0],domain),self.position(v[7],domain))
            .map_err(|_|OctreeError3::Invalid("unrepresentable physical octant"))
    }
    fn covering(&self,p:[i64;3])->Option<Octant3> {
        if p.iter().any(|&x| x<0 || x>=4*i64::from(self.extent())) { return None; }
        for level in (0..=self.max_level).rev() {
            let divisor=4i64<<(self.max_level-level);
            let c=Octant3 {level,index:p.map(|x|(x/divisor) as u32)};
            if self.leaves.contains(&c) { return Some(c); }
        }
        None
    }
    fn split(&mut self,c:Octant3)->Result<(),OctreeError3> {
        if c.level==self.max_level { return Err(OctreeError3::LevelBudget); }
        if self.leaves.len().checked_add(7).is_none_or(|n|n>self.max_leaves) { return Err(OctreeError3::LeafBudget); }
        self.leaves.remove(&c); self.leaves.extend(c.children()); Ok(())
    }
    /// Refine a set of CURRENT leaves once, then close strong 2:1 balance.
    /// Duplicate marks have no extra effect. Unknown marks, caps and cancellation
    /// return no partial tree; the original remains usable for retry.
    pub fn refined(&self,marks:&[Octant3],mut checkpoint:impl FnMut()->ControlFlow<()>)
        ->Result<Self,OctreeError3> {
        poll(&mut checkpoint)?;
        let marked:BTreeSet<_>=marks.iter().copied().collect();
        if marked.iter().any(|c|!self.leaves.contains(c)) { return Err(OctreeError3::Invalid("mark is not a current leaf")); }
        let mut next=self.clone();
        for c in marked { poll(&mut checkpoint)?; next.split(c)?; }
        loop {
            let mut coarse=BTreeSet::new();
            for &c in &next.leaves {
                poll(&mut checkpoint)?;
                let (lo,s)=next.lattice_box(c);
                for z in -1..=1 { for y in -1..=1 { for x in -1..=1 {
                    let d=[x,y,z]; if d==[0,0,0] {continue;}
                    let p=std::array::from_fn(|a|4*i64::from(lo[a])+match d[a] {
                        -1=>-1,0=>2*i64::from(s),_=>4*i64::from(s)+1,
                    });
                    if let Some(n)=next.covering(p) { if c.level>n.level+1 {coarse.insert(n);} }
                } } }
            }
            if coarse.is_empty() {break;}
            for c in coarse {poll(&mut checkpoint)?;next.split(c)?;}
        }
        poll(&mut checkpoint)?; Ok(next)
    }
    /// All shared face patches with exact tangential overlap. A coarse/fine
    /// interface produces four non-overlapping patches, not one full coarse face.
    pub fn faces(&self,mut checkpoint:impl FnMut()->ControlFlow<()>)->Result<Vec<OctreeFace3>,OctreeError3> {
        let mut result=Vec::new();
        for &c in &self.leaves {
            poll(&mut checkpoint)?;
            let (lo,s)=self.lattice_box(c);
            for axis in 0..3 {
                let (a,b)=((axis+1)%3,(axis+2)%3);let mut neighbors=BTreeSet::new();
                for i in [1,3] {for j in [1,3] {
                    let mut p=lo.map(|v|4*i64::from(v));p[axis]+=4*i64::from(s)+1;
                    p[a]+=i*i64::from(s);p[b]+=j*i64::from(s);
                    if let Some(n)=self.covering(p) {neighbors.insert(n);}
                }}
                for n in neighbors {
                    let (nl,ns)=self.lattice_box(n);
                    let lower=std::array::from_fn(|d|lo[d].max(nl[d]));
                    let upper=std::array::from_fn(|d|(lo[d]+s).min(nl[d]+ns));
                    if lower[axis]!=upper[axis] || lower[a]>=upper[a] || lower[b]>=upper[b] {
                        return Err(OctreeError3::Invalid("invalid face overlap"));
                    }
                    result.push(OctreeFace3 {lower:c,upper:n,axis,lo:lower,hi:upper});
                }
            }
        }
        poll(&mut checkpoint)?;Ok(result)
    }
    /// Direct Q1 hanging-node constraints. Each row is interpolation on a
    /// coarser face; callers recursively substitute rows before elimination.
    /// Dyadic half/quarter coefficients reproduce all trilinear polynomials.
    pub fn constraints(&self,mut checkpoint:impl FnMut()->ControlFlow<()>)->Result<BTreeMap<OctreeNode3,Vec<(OctreeNode3,f64)>>,OctreeError3> {
        let faces=self.faces(&mut checkpoint)?;
        let mut rows:BTreeMap<_,(u8,Vec<_>)>=BTreeMap::new();
        for face in faces {
            poll(&mut checkpoint)?;
            if face.lower.level==face.upper.level {continue;}
            let (coarse,fine)=if face.lower.level<face.upper.level {(face.lower,face.upper)}else{(face.upper,face.lower)};
            let vertices=self.corners(coarse);let (lo,s)=self.lattice_box(coarse);
            for node in self.corners(fine) {
                if node[face.axis]!=face.lo[face.axis] || vertices.contains(&node) {continue;}
                let t:[f64;3]=std::array::from_fn(|a|f64::from(node[a]-lo[a])/f64::from(s));
                let mut row=Vec::new();
                for (i,&v) in vertices.iter().enumerate() {
                    let w=(0..3).map(|a|if i&(1<<a)==0 {1.0-t[a]}else{t[a]}).product::<f64>();
                    if w>0.0 {row.push((v,w));}
                }
                row.sort_by_key(|&(n,_)|n);
                match rows.get(&node) {
                    Some((level,old)) if *level==coarse.level && *old!=row=>return Err(OctreeError3::Invalid("inconsistent face constraints")),
                    Some((level,_)) if *level<=coarse.level=>{},
                    _=>{rows.insert(node,(coarse.level,row));}
                }
            }
        }
        poll(&mut checkpoint)?;Ok(rows.into_iter().map(|(n,(_,r))|(n,r)).collect())
    }
}
