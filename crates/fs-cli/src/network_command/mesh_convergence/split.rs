//! Keep one problem-transfer path for uniform and marked meshes.
use super::*;
use fs_mesh::MarkedTetRefinement;

pub(super) enum Split { Uniform(TetRefinement), Marked(MarkedTetRefinement) }
impl Split {
    pub(super) fn build(cx:&Cx<'_>,root:&J,request:&Request,tets:&[[u32;4]],
        limits:TetRefinementLimits,marked:Option<&[usize]>)->Result<Self> {
        let Some(marked)=marked else {
            return TetRefinement::build(cx,request.mesh.positions(),tets,limits)
                .map(Self::Uniform).map_err(refinement_error);
        };
        let mut links=Vec::new();
        if let Some(contacts)=get(root,"solid")?.get("contacts") {
            for contact in array(contacts,"contacts",4096)? {
                // Independent planar traces need no synchronized edges. Their
                // overlap integration is rebuilt after transferring each side.
                // Exact subpatches may cease to match without changing the law.
                if contact.get("nonmatching").is_some() { continue; }
                for pair in array(get(contact,"face_pairs")?,"face_pairs",200_000)? {
                    poll(cx)?;
                    let a=indices::<3>(get(pair,"side_a")?,"side_a",request.mesh.vertex_count())?;
                    let b=indices::<3>(get(pair,"side_b")?,"side_b",request.mesh.vertex_count())?;
                    let b=match_vertices(request.mesh.positions(),a,b)?;
                    for (i,j) in [(0,1),(0,2),(1,2)] { links.push(([a[i],a[j]],[b[i],b[j]])); }
                }
            }
        }
        MarkedTetRefinement::build(cx,request.mesh.positions(),tets,marked,&links,limits)
            .map(Self::Marked).map_err(refinement_error)
    }
    pub(super) fn positions(&self)->&[[f64;3]] {
        match self { Self::Uniform(s)=>s.positions(),Self::Marked(s)=>s.positions() }
    }
    pub(super) fn tetrahedra(&self)->&[[u32;4]] {
        match self { Self::Uniform(s)=>s.tetrahedra(),Self::Marked(s)=>s.tetrahedra() }
    }
    pub(super) fn parent(&self,index:usize)->usize {
        match self { Self::Uniform(_)=>index/8,Self::Marked(s)=>s.parent_elements()[index] }
    }
    pub(super) fn prolongate(&self,cx:&Cx<'_>,field:&[f64])->Result<Vec<f64>> {
        match self { Self::Uniform(s)=>s.prolongate(cx,field),Self::Marked(s)=>s.prolongate(cx,field) }
            .map_err(refinement_error)
    }
    pub(super) fn face_children(&self,face:[u32;3])->Result<Vec<[u32;3]>> {
        match self {
            Self::Uniform(s)=>s.face_children(face).map(Vec::from),
            Self::Marked(s)=>s.face_children(face),
        }.map_err(refinement_error)
    }
    pub(super) fn contact_children(&self,a:[u32;3],b:[u32;3])->Result<Vec<([u32;3],[u32;3])>> {
        let b=match_vertices(self.positions(),a,b)?;
        let aa=self.face_children(a)?; let bb=self.face_children(b)?;
        if aa.len()!=bb.len() { return Err(producer("local refinement left unmatched contact subdivisions")); }
        let mut used=BTreeSet::new(); let mut pairs=Vec::new();
        for a in aa {
            let found=bb.iter().enumerate().find(|(index,b)|!used.contains(index)
                && a.iter().all(|&v|b.iter().any(|&w|self.positions()[v as usize]==self.positions()[w as usize])))
                .ok_or_else(||producer("refined contact triangles do not match geometrically"))?;
            used.insert(found.0); pairs.push((a,*found.1));
        }
        Ok(pairs)
    }
}
fn match_vertices(positions:&[[f64;3]],a:[u32;3],b:[u32;3])->Result<[u32;3]> {
    let mut matched=[0;3];
    for (i,&vertex) in a.iter().enumerate() {
        matched[i]=b.iter().copied().find(|&other|positions[vertex as usize]==positions[other as usize])
            .ok_or_else(||bad("mesh refinement requires exactly coincident matching contact vertices"))?;
    }
    let mut unique=matched; unique.sort_unstable();
    if unique.windows(2).any(|v|v[0]==v[1]) { return Err(bad("ambiguous contact vertex correspondence")); }
    Ok(matched)
}
