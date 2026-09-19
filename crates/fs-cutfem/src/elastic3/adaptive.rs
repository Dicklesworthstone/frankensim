//! Locally refined 3-D CutFEM with exact Q1 hanging-node elimination.
//!
//! Reuses the Cartesian kernel's cell blocks, bulk rules, scale-energy and
//! ghost representations. The physical apply is T^T K T; body loads use T^T b.
//! Ghost jumps are integrated once over each actual coarse/fine face patch.
//! T is independent of density, so the existing bulk-plus-ghost pullback stays
//! valid after prolongation. No coarsening, DWR estimate or continuum certificate.
use super::*;
use crate::octree3::{Octant3,Octree3,OctreeNode3,OctreeError3};
type Row=Vec<(OctreeNode3,f64)>;

/// A reduced operator on independent master nodes of a balanced octree.
pub struct AdaptiveElasticity3 {
    raw: CutElasticity3,
    rows: Vec<Vec<(usize,f64)>>,
    nodes: Vec<[f64;3]>,
    fixed: Vec<bool>,
    leaves: Vec<Octant3>,
    edges: Vec<(usize,usize,f64)>,
}
fn tree_error(error:OctreeError3)->ElasticityError3 {
    match error {
        OctreeError3::Cancelled=>ElasticityError3::Quadrature(QuadratureError3::Cancelled),
        _=>ElasticityError3::Invalid("octree topology/constraint construction refused"),
    }
}
fn terminal_row(node:OctreeNode3,direct:&BTreeMap<OctreeNode3,Row>,memo:&mut BTreeMap<OctreeNode3,Row>,
    path:&mut BTreeSet<OctreeNode3>,control:&mut QuadratureControl3<'_>)->Result<Row,ElasticityError3> {
    control.poll()?;
    if let Some(row)=memo.get(&node) {return Ok(row.clone());}
    let Some(row)=direct.get(&node) else {return Ok(vec![(node,1.0)]);};
    if path.len()>=64 || !path.insert(node) {return Err(ElasticityError3::Invalid("hanging-node cycle/depth budget"));}
    let mut out=BTreeMap::new();
    for &(parent,weight) in row {
        for (master,w) in terminal_row(parent,direct,memo,path,control)? {
            *out.entry(master).or_insert(0.0)+=weight*w;
            if out.len()>64 {return Err(ElasticityError3::Invalid("hanging row exceeds 64 terminal masters"));}
        }
    }
    path.remove(&node);
    let result:Row=out.into_iter().collect();
    if result.iter().map(|(_,w)|w).sum::<f64>()!=1.0 {return Err(ElasticityError3::Invalid("constraint partition of unity lost"));}
    memo.insert(node,result.clone());Ok(result)
}
impl AdaptiveElasticity3 {
    /// Integrate an admitted octree on a fixed implicit domain. `max_cells`
    /// bounds the full background and `max_dofs` bounds BOTH physical and master
    /// displacement spaces. At most 64 terminal masters per hanging row are
    /// admitted. Homogeneous box clamps select independent master nodes; slave
    /// values inherit those conditions through T. Every component still needs
    /// sufficient supports; a small residual is not a uniqueness certificate.
    #[allow(clippy::too_many_arguments)]
    pub fn build(domain:HexCell,tree:&Octree3,sdf:&dyn CutSdf3,material:&IsotropicElastic,
        clamp:&dyn Fn([f64;3])->bool,options:ElasticityOptions3,control:&mut QuadratureControl3<'_>)
        ->Result<Self,ElasticityError3> {
        control.poll()?;
        if tree.leaves().len()>options.max_cells || !options.ghost_gamma.is_finite() || options.ghost_gamma<0.0 {
            return Err(ElasticityError3::Invalid("invalid adaptive assembly allowance"));
        }
        let card=IsotropicElastic::new(material.youngs,material.poisson,material.strain_limit)
            .map_err(|_|ElasticityError3::Invalid("invalid isotropic material"))?;
        let (lambda,mu)=card.lame();
        if !lambda.is_finite()||!mu.is_finite()||mu<=0.0||!card.strain_limit.is_finite()||card.strain_limit<=0.0
            ||!((lambda+2.0*mu)/mu<=4.0) {return Err(ElasticityError3::Invalid("material outside compressible regime"));}
        let faces=tree.faces(||if control.poll().is_ok(){ControlFlow::Continue(())}else{ControlFlow::Break(())}).map_err(tree_error)?;
        let direct=tree.constraints(||if control.poll().is_ok(){ControlFlow::Continue(())}else{ControlFlow::Break(())}).map_err(tree_error)?;
        let mut cells=Vec::new();let mut leaves=Vec::new();let mut keys=BTreeSet::new();
        let mut volume_bounds=Interval::new(0.0,0.0);
        for &leaf in tree.leaves() {
            control.poll()?;
            let bounds=tree.bounds(leaf,domain).map_err(tree_error)?;
            let rules=cut_cell_rules3(sdf,bounds,control)?;
            volume_bounds=volume_bounds+rules.volume_bounds();
            if rules.bulk().is_empty() {
                if rules.cut_boxes()>0 {return Err(ElasticityError3::UnresolvedSupport(bounds));}
                continue;
            }
            let mut stiffness=[[0.0;24];24];
            for &(p,w) in rules.bulk() {
                control.poll()?;let (_,g)=q1(bounds,p);
                for i in 0..24 {for j in i..24 {
                    let (a,b,ci,cj)=(i/3,j/3,i%3,j%3);
                    let dot=g[a].iter().zip(g[b]).map(|(u,v)|u*v).sum::<f64>();
                    stiffness[i][j]+=w*(lambda*g[a][ci]*g[b][cj]+mu*g[a][cj]*g[b][ci]+if ci==cj {mu*dot}else{0.0});
                }}
            }
            for i in 0..24 {for j in i..24 {
                if !stiffness[i][j].is_finite(){return Err(ElasticityError3::Invalid("adaptive stiffness overflow"));}
                stiffness[j][i]=stiffness[i][j];
            }}
            control.poll()?;let enclosure=sdf.enclose(bounds.lo(),bounds.hi());control.poll()?;
            if !enclosure.lo().is_finite()||!enclosure.hi().is_finite()||enclosure.lo()>enclosure.hi() {
                return Err(ElasticityError3::Invalid("invalid assembly enclosure"));
            }
            keys.extend(tree.corners(leaf));
            if keys.len()>options.max_dofs/3 {return Err(ElasticityError3::Invalid("physical dof allowance"));}
            leaves.push(leaf);
            cells.push(Cell3 {key:tree.lattice_box(leaf).0.map(|v|v as usize),bounds,nodes:[0;8],stiffness,rules,cut:enclosure.hi()>=0.0});
        }
        if cells.is_empty(){return Err(ElasticityError3::EmptyDomain);}
        if !volume_bounds.hi().is_finite(){return Err(ElasticityError3::Invalid("volume overflow"));}
        let mut memo=BTreeMap::new();let mut masters=BTreeSet::new();let mut node_rows=Vec::new();
        for &node in &keys {
            let row=terminal_row(node,&direct,&mut memo,&mut BTreeSet::new(),control)?;
            masters.extend(row.iter().map(|&(n,_)|n));node_rows.push(row);
            if masters.len()>options.max_dofs/3 {return Err(ElasticityError3::Invalid("master dof allowance"));}
        }
        let master_ids:BTreeMap<_,_>=masters.iter().enumerate().map(|(i,&n)|(n,i)).collect();
        let raw_ids:BTreeMap<_,_>=keys.iter().enumerate().map(|(i,&n)|(n,i)).collect();
        let rows=node_rows.into_iter().map(|row|row.into_iter().map(|(n,w)|(master_ids[&n],w)).collect()).collect();
        let nodes:Vec<_>=masters.iter().map(|&n|tree.position(n,domain)).collect();
        let mut fixed=Vec::new();
        for (&key,&p) in masters.iter().zip(&nodes) {
            control.poll()?;let selected=clamp(p);control.poll()?;
            if selected && !(0..3).any(|a|key[a]==0||key[a]==tree.extent()) {return Err(ElasticityError3::Invalid("clamp is not on box boundary"));}
            fixed.push(selected);
        }
        if !fixed.iter().any(|v|*v){return Err(ElasticityError3::Invalid("no displacement support selected"));}
        for (cell,&leaf) in cells.iter_mut().zip(&leaves) {cell.nodes=tree.corners(leaf).map(|n|raw_ids[&n]);}
        let cell_ids:BTreeMap<_,_>=leaves.iter().enumerate().map(|(i,&c)|(c,i)).collect();
        let mut ghosts=Vec::new();let mut edges=Vec::new();
        for face in faces {
            control.poll()?;
            let (Some(&l),Some(&r))=(cell_ids.get(&face.lower),cell_ids.get(&face.upper)) else {continue;};
            let lo=tree.position(face.lo,domain);let hi=tree.position(face.hi,domain);
            let distance=f64::midpoint(cells[r].bounds.lo()[face.axis],cells[r].bounds.hi()[face.axis])
                -f64::midpoint(cells[l].bounds.lo()[face.axis],cells[l].bounds.hi()[face.axis]);
            if !distance.is_finite()||distance<=0.0 {return Err(ElasticityError3::Invalid("invalid neighbor distance"));}
            edges.push((l,r,distance));
            if options.ghost_gamma>0.0&&(cells[l].cut||cells[r].cut) {
                patch_ghosts(&cells,l,r,face.axis,lo,hi,options.ghost_gamma*mu,control,&mut ghosts)?;
            }
        }
        let raw_nodes:Vec<_>=keys.iter().map(|&n|tree.position(n,domain)).collect();
        let raw=CutElasticity3 {fixed:vec![false;raw_nodes.len()],nodes:raw_nodes,scales:vec![1.0;cells.len()],cells,ghosts,
            volume_bounds:Interval::new(volume_bounds.lo().max(0.0),volume_bounds.hi())};
        control.poll()?;Ok(Self {raw,rows,nodes,fixed,leaves,edges})
    }
    /// Independent master positions; reduced displacement coefficients use this order.
    #[must_use] pub fn nodes(&self)->&[[f64;3]] {&self.nodes}
    /// Homogeneous master-node clamp mask.
    #[must_use] pub fn fixed(&self)->&[bool] {&self.fixed}
    /// Physical active-node positions, including eliminated hanging nodes.
    #[must_use] pub fn physical_nodes(&self)->&[[f64;3]] {self.raw.nodes()}
    /// Physical cell connectivity, indexing `physical_nodes()`.
    #[must_use] pub fn cell_nodes(&self)->Vec<[usize;8]> {self.raw.cell_nodes()}
    /// Active leaf identities, in scale/volume order.
    #[must_use] pub fn leaves(&self)->&[Octant3] {&self.leaves}
    /// Active face neighbors and their positive normal center separation.
    #[must_use] pub fn filter_edges(&self)->&[(usize,usize,f64)] {&self.edges}
    /// Active design dimension.
    #[must_use] pub fn cells(&self)->usize {self.raw.cells()}
    /// Retained numerical cut measures.
    #[must_use] pub fn volumes(&self)->Vec<f64> {self.raw.volumes()}
    /// Conservative geometry-volume enclosure, not a solution error bound.
    #[must_use] pub fn volume_bounds(&self)->Interval {self.raw.volume_bounds()}
    /// Positive relative cell stiffnesses.
    #[must_use] pub fn scales(&self)->&[f64] {self.raw.scales()}
    /// Update material distribution without reintegrating or changing constraints.
    pub fn set_scales(&mut self,s:&[f64])->Result<(),ElasticityError3> {self.raw.set_scales(s)}
    fn expand(&self,x:&[f64])->Vec<f64> {
        self.rows.iter().flat_map(|row|(0..3).map(move |c|row.iter().map(|&(m,w)|if self.fixed[m]{0.0}else{w*x[3*m+c]}).sum())).collect()
    }
    fn restrict(&self,x:&[f64],y:&mut[f64]) {
        y.fill(0.0);
        for (i,row) in self.rows.iter().enumerate() {for &(m,w) in row {if !self.fixed[m] {for c in 0..3 {y[3*m+c]+=w*x[3*i+c];}}}}
    }
    /// Reconstruct the conforming physical nodal field for observation/export.
    pub fn physical_displacements(&self,x:&[f64])->Result<Vec<f64>,ElasticityError3> {
        if x.len()!=self.n()||!x.iter().all(|v|v.is_finite()){return Err(ElasticityError3::Invalid("invalid master field"));}
        Ok(self.expand(x))
    }
    /// Integrate a dead body load and apply exactly the transpose of reconstruction.
    pub fn body_load(&self,f:&dyn Fn([f64;3])->[f64;3],mut poll:impl FnMut()->ControlFlow<()>)
        ->Result<Vec<f64>,ElasticityError3> {
        let raw=self.raw.body_load(f,&mut poll)?;let mut rhs=vec![0.0;self.n()];self.restrict(&raw,&mut rhs);
        if poll().is_break(){return Err(ElasticityError3::Cancelled);}
        if !rhs.iter().all(|v|v.is_finite()){return Err(ElasticityError3::Invalid("reduced load overflow"));}Ok(rhs)
    }
    /// Full discrete density contraction after Q1 prolongation, including ghost terms.
    pub fn scale_quadratic_forms(&self,x:&[f64])->Result<Vec<f64>,ElasticityError3> {
        self.raw.scale_quadratic_forms(&self.physical_displacements(x)?)
    }
    /// Solve the reduced operator with the canonical CG recurrence and an explicit
    /// true-residual gate. No physical or master field is returned on a stop.
    /// Identity preconditioning is the initial adaptive path; there is no multigrid
    /// or mesh-independent iteration claim. Budget is shared by all CG batches.
    pub fn solve_controlled(&self,force:&[f64],tol:f64,max_iters:usize,poll_iters:usize,
        mut checkpoint:impl FnMut(usize)->ControlFlow<()>)->Result<ElasticitySolution3,ElasticityError3> {
        if force.len()!=self.n()||!force.iter().all(|v|v.is_finite())||!tol.is_finite()||tol<=0.0||tol>=1.0||poll_iters==0 {
            return Err(ElasticityError3::Invalid("invalid adaptive solve controls"));
        }
        if checkpoint(0).is_break(){return Err(ElasticityError3::Cancelled);}
        let rhs:Vec<f64>=force.iter().enumerate().map(|(i,&v)|if self.fixed[i/3]{0.0}else{v}).collect();
        if rhs.iter().all(|v|*v==0.0) {
            let x=vec![0.0;self.n()];let residual=recomputed_euclidean_residual_claim(self,&x,&rhs);
            if checkpoint(0).is_break(){return Err(ElasticityError3::Cancelled);}
            return Ok(ElasticitySolution3 {coefficients:x,compliance:0.0,iterations:0,residual});
        }
        let precond=fs_sparse::precond::IdentityPrecond;
        let mut state=CgState::new(self,&precond,&rhs);
        loop {
            if checkpoint(state.iters).is_break(){return Err(ElasticityError3::Cancelled);}
            if !state.rel_residual().is_finite(){return Err(ElasticityError3::NotConverged{iterations:state.iters,relative_residual:state.rel_residual()});}
            if state.rel_residual()<0.1*tol || state.iters>=max_iters {
                let residual=recomputed_euclidean_residual_claim(self,&state.x,&rhs);
                if checkpoint(state.iters).is_break(){return Err(ElasticityError3::Cancelled);}
                let r=residual.euclidean().expect("explicit Euclidean residual");
                if !(r<tol) {return Err(ElasticityError3::NotConverged{iterations:state.iters,relative_residual:r});}
                let compliance:f64=rhs.iter().zip(&state.x).map(|(b,x)|b*x).sum();
                if !compliance.is_finite(){return Err(ElasticityError3::Invalid("adaptive compliance overflow"));}
                if checkpoint(state.iters).is_break(){return Err(ElasticityError3::Cancelled);}
                return Ok(ElasticitySolution3 {coefficients:state.x,compliance,iterations:state.iters,residual});
            }
            let before=state.iters;let _=state.run(self,&precond,0.1*tol,poll_iters.min(max_iters-before));
            state.history.clear();
            if state.iters==before {return Err(ElasticityError3::NotConverged{iterations:before,relative_residual:state.rel_residual()});}
        }
    }
}
impl LinearOp for AdaptiveElasticity3 {
    fn n(&self)->usize {3*self.nodes.len()}
    fn apply(&self,x:&[f64],y:&mut[f64]) {
        assert_eq!(x.len(),self.n());assert_eq!(y.len(),self.n());
        let full=self.expand(x);let mut result=vec![0.0;full.len()];self.raw.apply(&full,&mut result);self.restrict(&result,y);
        for (m,&fixed) in self.fixed.iter().enumerate(){if fixed {y[3*m..3*m+3].copy_from_slice(&x[3*m..3*m+3]);}}
    }
}
#[allow(clippy::too_many_arguments)]
fn patch_ghosts(cells:&[Cell3],left:usize,right:usize,axis:usize,lo:[f64;3],hi:[f64;3],coefficient:f64,
    control:&mut QuadratureControl3<'_>,out:&mut Vec<GhostPoint3>)->Result<(),ElasticityError3> {
    let l=&cells[left];let r=&cells[right];let (a,b)=((axis+1)%3,(axis+2)%3);
    let h=(l.bounds.hi()[axis]-l.bounds.lo()[axis]).min(r.bounds.hi()[axis]-r.bounds.lo()[axis]);
    let gauss=[(-0.774_596_669_241_483_4,5.0/9.0),(0.0,8.0/9.0),(0.774_596_669_241_483_4,5.0/9.0)];
    for (u,wu) in gauss {for (v,wv) in gauss {
        control.poll()?;let mut p=lo;p[a]=f64::midpoint(lo[a],hi[a])+0.5*(hi[a]-lo[a])*u;
        p[b]=f64::midpoint(lo[b],hi[b])+0.5*(hi[b]-lo[b])*v;
        let (_,gl)=q1(l.bounds,p);let (_,gr)=q1(r.bounds,p);let mut jump=BTreeMap::new();
        for i in 0..8 {*jump.entry(l.nodes[i]).or_insert(0.0)+=gl[i][axis];*jump.entry(r.nodes[i]).or_insert(0.0)-=gr[i][axis];}
        let weight=coefficient*h*0.25*(hi[a]-lo[a])*(hi[b]-lo[b])*wu*wv;
        if !weight.is_finite()||weight<=0.0||!jump.values().all(|j:&f64|j.is_finite()){return Err(ElasticityError3::Invalid("adaptive ghost overflow"));}
        let (nodes,jump)=jump.into_iter().unzip();out.push(GhostPoint3{cells:[left,right],nodes,jump,weight});
    }}Ok(())
}
