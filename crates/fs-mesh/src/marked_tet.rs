//! Local conforming refinement by complete edge-star bisection.
//!
//! Marked cells nominate their longest original edge. Every cell incident to
//! a nominated edge is split, so a face cannot acquire a hanging midpoint.
//! Explicit coincident-edge pairs close contact constraints WITHOUT identifying
//! their vertices. Geometric edge ordering makes both contact triangulations
//! agree even when their vertex IDs/local orders differ.
//!
//! This is not Delaunay refinement or a shape-regularity theorem. It preserves
//! an admitted parent's topology and P1 fields, subject to floating-point
//! representability checks. Output caps are not a total workspace memory cap.

use std::collections::{BTreeMap, BTreeSet};
use fs_exec::Cx;
use fs_ivl::{orient3d, Sign};
use crate::{TetRefinementError as Error, TetRefinementLimits};

type Edge = [u32; 2];
const PAIRS: [(usize, usize); 6] = [(0,1),(0,2),(0,3),(1,2),(1,3),(2,3)];
fn edge(a:u32,b:u32)->Edge { if a < b { [a,b] } else { [b,a] } }
fn edges(t:[u32;4])->[Edge;6] { PAIRS.map(|(a,b)|edge(t[a],t[b])) }
fn check(cx:&Cx<'_>)->Result<(),Error> { cx.checkpoint().map_err(|_|Error::Cancelled) }
fn point_order(a:[f64;3],b:[f64;3])->std::cmp::Ordering {
    for (a,b) in a.into_iter().zip(b) {
        // Coordinate equality admits signed zero; use the same equality here.
        if a != b { return a.total_cmp(&b); }
    }
    std::cmp::Ordering::Equal
}
fn endpoints(p:&[[f64;3]],e:Edge)->[[f64;3];2] {
    let a=p[e[0] as usize]; let b=p[e[1] as usize];
    if point_order(a,b).is_gt() { [b,a] } else { [a,b] }
}
fn sign(p:&[[f64;3]],t:[u32;4])->Sign {
    orient3d(p[t[0] as usize],p[t[1] as usize],p[t[2] as usize],p[t[3] as usize])
}
fn longest(p:&[[f64;3]],t:[u32;4])->Result<Edge,Error> {
    let es=edges(t);
    let mut differences=[[0.0;3];6]; let mut scale=0.0_f64;
    for (d,e) in differences.iter_mut().zip(es) {
        for axis in 0..3 {
            d[axis]=p[e[0] as usize][axis]-p[e[1] as usize][axis];
            if !d[axis].is_finite() { return Err(Error::InvalidMesh); }
            scale=scale.max(d[axis].abs());
        }
    }
    if scale==0.0 { return Err(Error::InvalidMesh); }
    let mut best=es[0]; let mut length=-1.0;
    for (d,e) in differences.into_iter().zip(es) {
        let n=d.into_iter().map(|x|(x/scale)*(x/scale)).sum::<f64>();
        if n>length || (n==length && e<best) { best=e; length=n; }
    }
    Ok(best)
}

/// An immutable, complete local refinement with explicit cell and P1 lineage.
#[derive(Debug)]
pub struct MarkedTetRefinement {
    original_vertices:usize,
    positions:Vec<[f64;3]>,
    tetrahedra:Vec<[u32;4]>,
    parent_elements:Vec<usize>,
    midpoint_parents:Vec<Edge>,
    split_edges:BTreeMap<Edge,(usize,u32)>,
    faces:BTreeSet<[u32;3]>,
}
impl MarkedTetRefinement {
    /// Refine the marked cells' longest edges and their full incident stars.
    /// `paired_edges` declares separately numbered, exactly coincident contact
    /// edges that must be split together. It does not weld either trace.
    /// Mark order and duplicate marks do not affect the result.
    ///
    /// # Errors
    /// Refuses invalid/empty marks, malformed input or edge constraints,
    /// exceeded output limits, cancellation, or unrepresentable children.
    pub fn build(cx:&Cx<'_>,positions:&[[f64;3]],tetrahedra:&[[u32;4]],
        marked:&[usize],paired_edges:&[(Edge,Edge)],limits:TetRefinementLimits)->Result<Self,Error> {
        check(cx)?;
        if positions.len()<4 || tetrahedra.is_empty() || marked.is_empty() {
            return Err(Error::InvalidMesh);
        }
        if positions.len()>limits.max_vertices || positions.len()>u32::MAX as usize
            || tetrahedra.len()>limits.max_tetrahedra { return Err(Error::OutputLimit); }
        for (i,p) in positions.iter().enumerate() {
            if i%256==0 { check(cx)?; }
            if p.iter().any(|x|!x.is_finite()) { return Err(Error::InvalidMesh); }
        }
        let mut stars:BTreeMap<Edge,BTreeSet<usize>>=BTreeMap::new();
        let mut faces=BTreeMap::new(); let mut unique=BTreeSet::new();
        for (i,&t) in tetrahedra.iter().enumerate() {
            if i%256==0 { check(cx)?; }
            let mut key=t; key.sort_unstable();
            if key.windows(2).any(|w|w[0]==w[1]) || key[3] as usize>=positions.len()
                || !unique.insert(key) || sign(positions,t)==Sign::Zero { return Err(Error::InvalidMesh); }
            for e in edges(t) { stars.entry(e).or_default().insert(i); }
            for omitted in 0..4 {
                let mut face=[0;3]; let mut slot=0;
                for (j,&v) in key.iter().enumerate() { if j!=omitted { face[slot]=v; slot+=1; } }
                let count=faces.entry(face).or_insert(0_usize); *count+=1;
                if *count>2 { return Err(Error::InvalidMesh); }
            }
        }
        let mut selected=BTreeSet::new();
        for (i,&index) in marked.iter().enumerate() {
            if i%256==0 { check(cx)?; }
            selected.insert(longest(positions,*tetrahedra.get(index).ok_or(Error::InvalidMesh)?)?);
        }
        let mut links:BTreeMap<Edge,BTreeSet<Edge>>=BTreeMap::new();
        for (i,&(a,b)) in paired_edges.iter().enumerate() {
            if i%256==0 { check(cx)?; }
            let a=edge(a[0],a[1]); let b=edge(b[0],b[1]);
            if !stars.contains_key(&a) || !stars.contains_key(&b)
                || endpoints(positions,a)!=endpoints(positions,b) { return Err(Error::InvalidMesh); }
            links.entry(a).or_default().insert(b); links.entry(b).or_default().insert(a);
        }
        let mut pending:Vec<_>=selected.iter().copied().collect();
        while let Some(e)=pending.pop() {
            check(cx)?;
            if let Some(others)=links.get(&e) {
                for &other in others { if selected.insert(other) { pending.push(other); } }
            }
        }
        let vertices=positions.len().checked_add(selected.len())
            .filter(|&n|n<=limits.max_vertices && n<=u32::MAX as usize).ok_or(Error::OutputLimit)?;
        let mut selected:Vec<_>=selected.into_iter().collect();
        selected.sort_by(|&a,&b| {
            let aa=endpoints(positions,a); let bb=endpoints(positions,b);
            point_order(aa[0],bb[0]).then(point_order(aa[1],bb[1])).then(a.cmp(&b))
        });
        let mut result=Self { original_vertices:positions.len(),positions:positions.to_vec(),
            tetrahedra:tetrahedra.to_vec(),parent_elements:(0..tetrahedra.len()).collect(),
            midpoint_parents:Vec::new(),split_edges:BTreeMap::new(),faces:faces.into_keys().collect() };
        result.positions.try_reserve_exact(vertices-positions.len()).map_err(|_|Error::OutputLimit)?;
        for (ordinal,e) in selected.into_iter().enumerate() {
            check(cx)?;
            let incident:Vec<_>=stars.get(&e).ok_or(Error::InvalidMesh)?.iter().copied().collect();
            if incident.is_empty() { return Err(Error::InvalidMesh); }
            result.tetrahedra.len().checked_add(incident.len())
                .filter(|&n|n<=limits.max_tetrahedra).ok_or(Error::OutputLimit)?;
            result.tetrahedra.try_reserve(incident.len()).map_err(|_|Error::OutputLimit)?;
            result.parent_elements.try_reserve(incident.len()).map_err(|_|Error::OutputLimit)?;
            let m=result.positions.len() as u32;
            let midpoint=std::array::from_fn(|axis|
                f64::midpoint(positions[e[0] as usize][axis],positions[e[1] as usize][axis]));
            if midpoint==positions[e[0] as usize] || midpoint==positions[e[1] as usize] {
                return Err(Error::UnrepresentableSplit);
            }
            result.positions.push(midpoint); result.midpoint_parents.push(e);
            result.split_edges.insert(e,(ordinal,m));
            for (at,index) in incident.into_iter().enumerate() {
                if at%256==0 { check(cx)?; }
                let old=result.tetrahedra[index]; let parent=result.parent_elements[index];
                let mut left=old; let mut right=old;
                for v in &mut left { if *v==e[1] { *v=m; } }
                for v in &mut right { if *v==e[0] { *v=m; } }
                let want=sign(&result.positions,old);
                if sign(&result.positions,left)!=want || sign(&result.positions,right)!=want {
                    return Err(Error::UnrepresentableSplit);
                }
                for key in edges(old) { stars.get_mut(&key).ok_or(Error::InvalidMesh)?.remove(&index); }
                let appended=result.tetrahedra.len();
                result.tetrahedra[index]=left;
                result.tetrahedra.push(right); result.parent_elements.push(parent);
                for (id,t) in [(index,left),(appended,right)] {
                    for key in edges(t) { stars.entry(key).or_default().insert(id); }
                }
            }
        }
        check(cx)?; Ok(result)
    }
    /// Refined positions; every original vertex and its number are retained.
    #[must_use]
    pub fn positions(&self)->&[[f64;3]] { &self.positions }
    /// Refined cells, including cells changed only by conformity closure.
    #[must_use]
    pub fn tetrahedra(&self)->&[[u32;4]] { &self.tetrahedra }
    /// Original parent of each returned cell; use this for material transfer.
    #[must_use]
    pub fn parent_elements(&self)->&[usize] { &self.parent_elements }
    /// Parent edges of appended midpoint vertices, in output order.
    #[must_use]
    pub fn midpoint_parents(&self)->&[Edge] { &self.midpoint_parents }
    /// Prolongate the exact parent P1 function without changing its support.
    ///
    /// # Errors
    /// Refuses wrong arity, non-finite values, or cancellation.
    pub fn prolongate(&self,cx:&Cx<'_>,field:&[f64])->Result<Vec<f64>,Error> {
        check(cx)?;
        if field.len()!=self.original_vertices || field.iter().any(|v|!v.is_finite()) {
            return Err(Error::InvalidField);
        }
        let mut out=field.to_vec(); out.try_reserve(self.midpoint_parents.len()).map_err(|_|Error::OutputLimit)?;
        for (i,&[a,b]) in self.midpoint_parents.iter().enumerate() {
            if i%256==0 { check(cx)?; }
            out.push(f64::midpoint(field[a as usize],field[b as usize]));
        }
        Ok(out)
    }
    /// Subdivide one original triangle using the exact same edge order as the
    /// volume. Unchanged faces return one triangle, not an invented refinement.
    ///
    /// # Errors
    /// Refuses a face absent from the original tetrahedral complex.
    pub fn face_children(&self,face:[u32;3])->Result<Vec<[u32;3]>,Error> {
        let mut key=face; key.sort_unstable();
        if !self.faces.contains(&key) { return Err(Error::UnknownFace); }
        let mut operations=Vec::new();
        for (a,b) in [(0,1),(0,2),(1,2)] {
            let e=edge(face[a],face[b]);
            if let Some(&(order,m))=self.split_edges.get(&e) { operations.push((order,e,m)); }
        }
        operations.sort_unstable(); let mut out=vec![face];
        for (_,e,m) in operations {
            let mut next=Vec::new();
            for t in out {
                if t.contains(&e[0]) && t.contains(&e[1]) {
                    let mut a=t; let mut b=t;
                    for v in &mut a { if *v==e[1] { *v=m; } }
                    for v in &mut b { if *v==e[0] { *v=m; } }
                    next.push(a); next.push(b);
                } else { next.push(t); }
            }
            out=next;
        }
        Ok(out)
    }
}
