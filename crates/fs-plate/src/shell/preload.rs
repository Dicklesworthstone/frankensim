//! Full-coordinate static loading and small vibrations about equilibrium.
//!
//! Uses the SAME Green--Lagrange P1 membrane energy as `shell::reduction`,
//! together with the existing linear DKT, drilling and bonded-beam stiffness.
//! The original membrane is replaced, not counted twice. Newton/globalization
//! is injected by the caller (the piano uses fs-solver); factorization and
//! stability checks belong to fs-sparse. This avoids a plate/FEEC dependency cycle.
//! No reduced modal basis is used to find static equilibrium: in-plane response
//! must not disappear merely because its modes lie above the audible band.
//!
//! Dead forces and optional conservative anchored tethers, fixed supports,
//! moderate rotations. Beam bending/axial laws and offsets remain linear in
//! their reference frames. No nonconservative follower loads, beam buckling,
//! glue creep, manufacturing stress, large-rotation bending or unstable branch.
use super::{FacetGeometry, ShellMesh, ShellModel, ShellSupport};
use super::stiffened::{ShellBeam, assemble_stiffened_shell};
use crate::PlateSection;
use fs_sparse::{Coo, Csr, DirectOrdering, LdltFactor, LdltOptions, SymbolicLdlt};
use std::cell::Cell;

/// Direction-updating tensile connections at existing shell point ports.
pub mod tethers;
use tethers::{PreparedTether, ShellTether, TetherResponse};

/// Explicit cold work and moderate-deformation admission limits.
#[derive(Clone, Copy, Debug)]
pub struct PreloadOptions {
    /// Monotone equal load increments; each must converge before the next.
    pub load_steps: usize,
    /// Maximum Newton iterations per increment.
    pub max_iterations: usize,
    /// Combined residual/Jacobian evaluations across the entire solve.
    pub max_evaluations: usize,
    /// Maximum norm of either midsurface displacement derivative.
    pub maximum_gradient: f64,
    /// Maximum norm of a nodal physical rotation, in radians.
    pub maximum_rotation_rad: f64,
    /// Final physical force residual relative to the applied equivalent-force norm.
    pub relative_force_tolerance: f64,
}
impl Default for PreloadOptions {
    fn default() -> Self {
        Self { load_steps: 8, max_iterations: 32, max_evaluations: 20_000,
            maximum_gradient: 0.1, maximum_rotation_rad: 0.1,
            relative_force_tolerance: 1e-7 }
    }
}

/// Stable incremental pencil and the reference-to-equilibrium displacement.
#[derive(Debug)]
pub struct PreloadedShell {
    /// Consistent equilibrium tangent INCLUDING any massless tethers; reference M.
    /// Do not add these same tethers again in a downstream dynamic model.
    pub model: ShellModel,
    /// Full nodal u,v,w [m], theta_x,y,z [rad], including zero fixed coordinates.
    pub displacement: Vec<f64>,
    /// Norm of unbalanced force; moments are divided by the declared length.
    pub residual_force_n: f64,
    /// Length used only to scale rotations/moments in solver norms [m].
    pub norm_length_m: f64,
    /// Shell/beam stored energy, excluding external/tether potential [J].
    pub stored_energy_j: f64,
    /// Final forces, tensions, lengths and potential changes in input order.
    pub tether_responses: Vec<TetherResponse>,
    /// Accepted Newton attempts across all increments.
    pub iterations: usize,
    /// Actual cold residual/Jacobian evaluations.
    pub evaluations: usize,
}

fn dot(a: [f64;3], b: [f64;3]) -> f64 { a.iter().zip(b).map(|(a,b)|a*b).sum() }
fn mul(a: &[f64;9], b: [f64;3]) -> [f64;3] {
    std::array::from_fn(|i| (0..3).map(|j| a[3*i+j]*b[j]).sum())
}
fn norm(x: &[f64]) -> f64 { x.iter().fold(0.0_f64, |a,b| a.hypot(*b)) }

struct Element {
    geometry: FacetGeometry,
    membrane: [f64;9], // area-integrated resultant modulus, N m
    dofs: [Option<usize>;9],
    linear: [[f64;3];9],
}
struct Correction {
    energy: f64,
    force: [f64;9],
    tangent: [f64;81],
    gradient: f64,
}
impl Element {
    fn new(mesh: &ShellMesh, model: &ShellModel, section: &PlateSection, e: usize)
        -> Result<Self,String> {
        let geometry = mesh.facet(e).map_err(|e|e.to_string())?;
        let membrane = section.d.map(|d| d*geometry.area_m2*12.0/section.thickness.powi(2));
        let dofs = std::array::from_fn(|i| model.dof_map[6*mesh.tris[e][i/3]+i%3]);
        let linear = std::array::from_fn(|i| {
            let [x,y] = geometry.gradient[i/3]; let c = i%3;
            [x*geometry.frame[0][c], y*geometry.frame[1][c],
                y*geometry.frame[0][c]+x*geometry.frame[1][c]]
        });
        if membrane.iter().any(|x|!x.is_finite()) { return Err("shell membrane modulus overflow".into()); }
        Ok(Self { geometry, membrane, dofs, linear })
    }
    // Polynomial DIFFERENCES avoid subtracting two nearly equal large membrane
    // forces/stiffnesses. At zero displacement every correction is exactly zero.
    fn correction(&self, u: &[f64]) -> Correction {
        let at = |a:usize,c:usize| self.dofs[3*a+c].map_or(0.0,|i|u[i]);
        let mut dx = [0.;3]; let mut dy = [0.;3];
        for a in 1..3 { for c in 0..3 {
            let v = at(a,c)-at(0,c);
            dx[c] += self.geometry.gradient[a][0]*v;
            dy[c] += self.geometry.gradient[a][1]*v;
        }}
        let linear = [dot(self.geometry.frame[0],dx),dot(self.geometry.frame[1],dy),
            dot(self.geometry.frame[0],dy)+dot(self.geometry.frame[1],dx)];
        let quadratic = [0.5*dot(dx,dx),0.5*dot(dy,dy),dot(dx,dy)];
        let sl = mul(&self.membrane,linear);
        let sq = mul(&self.membrane,quadratic);
        let stress = std::array::from_fn(|i|sl[i]+sq[i]);
        let nonlinear: [[f64;3];9] = std::array::from_fn(|i| {
            let [x,y] = self.geometry.gradient[i/3]; let c = i%3;
            [x*dx[c],y*dy[c],y*dx[c]+x*dy[c]]
        });
        let mut out = Correction { energy: dot(linear,sq)+0.5*dot(quadratic,sq),
            force:[0.;9],tangent:[0.;81],gradient:norm(&dx).max(norm(&dy)) };
        for i in 0..9 {
            out.force[i] = dot(self.linear[i],sq)+dot(nonlinear[i],stress);
            for j in 0..=i {
                let mut v = dot(self.linear[i],mul(&self.membrane,nonlinear[j]))
                    +dot(nonlinear[i],mul(&self.membrane,self.linear[j]))
                    +dot(nonlinear[i],mul(&self.membrane,nonlinear[j]));
                if i%3 == j%3 {
                    let [ix,iy] = self.geometry.gradient[i/3];
                    let [jx,jy] = self.geometry.gradient[j/3];
                    v += ix*jx*stress[0]+iy*jy*stress[1]+(ix*jy+iy*jx)*stress[2];
                }
                out.tangent[9*i+j] = v; out.tangent[9*j+i] = v;
            }
        }
        out
    }
}

/// Cold nonlinear residual used by an injected solver. Coordinates scale physical
/// rotations by the reference span [m]; residuals are reference-K-inverse force
/// imbalance in those same coordinates. Final admission separately checks forces.
/// Invalid trials/evaluation-budget exhaustion return NaNs for solver rejection.
pub struct PreloadProblem<'a> {
    model: &'a ShellModel,
    elements: Vec<Element>,
    tethers: Vec<PreparedTether>,
    base: LdltFactor,
    scale: Vec<f64>,
    load: Vec<f64>,
    fraction: f64,
    options: PreloadOptions,
    evaluations: Cell<usize>,
}
impl PreloadProblem<'_> {
    fn physical(&self, x:&[f64]) -> Vec<f64> {
        x.iter().zip(&self.scale).map(|(x,s)|x/s).collect()
    }
    fn admitted(&self, u:&[f64]) -> bool {
        if u.len()!=self.model.free || u.iter().any(|v|!v.is_finite()) { return false; }
        for node in self.model.dof_map.chunks_exact(6) {
            let rotation: [f64;3] = std::array::from_fn(|i|node[i+3].map_or(0.,|d|u[d]));
            if norm(&rotation)>self.options.maximum_rotation_rad { return false; }
        }
        true
    }
    fn reserve(&self) -> bool {
        let count=self.evaluations.get();
        if count>=self.options.max_evaluations { return false; }
        self.evaluations.set(count+1); true
    }
    fn force_energy(&self,u:&[f64],force:&mut[f64]) -> Option<f64> {
        if !self.admitted(u) || force.len()!=self.model.free { return None; }
        self.model.k.spmv(u,force);
        let mut energy=0.5*u.iter().zip(&*force).map(|(a,b)|a*b).sum::<f64>();
        for element in &self.elements {
            let c=element.correction(u);
            if !c.gradient.is_finite() || c.gradient>self.options.maximum_gradient {return None;}
            energy+=c.energy;
            for i in 0..9 {if let Some(d)=element.dofs[i] {force[d]+=c.force[i];}}
        }
        (energy.is_finite() && force.iter().all(|v|v.is_finite())).then_some(energy)
    }
    fn applied_load(&self,u:&[f64])->Result<Vec<f64>,String> {
        let mut load=self.load.clone();
        for tether in &self.tethers {tether.add_gradient(u,&mut load,-1.)?;}
        for p in &mut load {*p*=self.fraction;}
        if load.iter().any(|v|!v.is_finite()) {return Err("nonfinite combined shell load".into());}
        Ok(load)
    }
    fn check_balance(&self,u:&[f64])->Result<f64,String> {
        let mut force=vec![0.;self.model.free];
        self.force_energy(u,&mut force).ok_or("equilibrium exceeds finite/moderate-deformation limits")?;
        let load=self.applied_load(u)?;
        let residual:Vec<_>=force.iter().zip(&load).zip(&self.scale)
            .map(|((f,p),s)|(f-p)/s).collect();
        let load_norm=norm(&load.iter().zip(&self.scale).map(|(p,s)|p/s).collect::<Vec<_>>());
        let residual=norm(&residual);
        if !residual.is_finite() || residual>1e-8+self.options.relative_force_tolerance*load_norm {
            return Err(format!("shell preload failed physical force balance: {residual} N for {load_norm} N load"));
        }
        Ok(residual)
    }
    fn tangent(&self,u:&[f64]) -> Result<Csr,String> {
        let mut k=Coo::new(self.model.free,self.model.free);
        for row in 0..self.model.free {
            let (columns,values)=self.model.k.row(row);
            for (&col,&value) in columns.iter().zip(values) {k.push(row,col,value);}
        }
        for element in &self.elements {
            let c=element.correction(u);
            if !c.gradient.is_finite() || c.gradient>self.options.maximum_gradient
                || c.tangent.iter().any(|v|!v.is_finite()) {return Err("invalid equilibrium tangent".into());}
            for i in 0..9 {for j in 0..9 {
                if let (Some(a),Some(b))=(element.dofs[i],element.dofs[j]) {k.push(a,b,c.tangent[9*i+j]);}
            }}
        }
        for tether in &self.tethers {tether.add_tangent(u,&mut k,self.fraction)?;}
        Ok(k.assemble())
    }
}
impl PreloadProblem<'_> {
    /// Number of unconstrained coordinates.
    pub fn dimension(&self)->usize {self.model.free}
    /// Overwrite the preconditioned shell minus current applied-load residual.
    pub fn residual(&self,x:&[f64],out:&mut[f64]) {
        if !self.reserve() || x.len()!=self.model.free || out.len()!=self.model.free {out.fill(f64::NAN);return;}
        let u=self.physical(x);
        if self.force_energy(&u,out).is_none() {out.fill(f64::NAN);return;}
        let Ok(load)=self.applied_load(&u) else {out.fill(f64::NAN);return;};
        for (f,load) in out.iter_mut().zip(load) {*f-=load;}
        let value=self.base.solve(out);
        for ((out,v),scale) in out.iter_mut().zip(value).zip(&self.scale) {*out=v*scale;}
    }
    /// Overwrite the analytic Jacobian action in the scaled solver coordinates.
    pub fn jacobian_apply(&self,x:&[f64],direction:&[f64],out:&mut[f64]) {
        if !self.reserve() || x.len()!=self.model.free || direction.len()!=self.model.free
            || out.len()!=self.model.free {out.fill(f64::NAN);return;}
        let u=self.physical(x);let v=self.physical(direction);
        if !self.admitted(&u) {out.fill(f64::NAN);return;}
        self.model.k.spmv(&v,out);
        for element in &self.elements {
            let c=element.correction(&u);
            if !c.gradient.is_finite() || c.gradient>self.options.maximum_gradient {out.fill(f64::NAN);return;}
            for i in 0..9 {if let Some(a)=element.dofs[i] {
                for j in 0..9 {if let Some(b)=element.dofs[j] {out[a]+=c.tangent[9*i+j]*v[b];}}
            }}
        }
        for tether in &self.tethers {
            if tether.add_action(&u,&v,out,self.fraction).is_err() {out.fill(f64::NAN);return;}
        }
        let value=self.base.solve(out);
        for ((out,v),scale) in out.iter_mut().zip(value).zip(&self.scale) {*out=v*scale;}
    }
}
fn positive_factor(k:&Csr)->Result<LdltFactor,String> {
    let factor=SymbolicLdlt::analyze(k,DirectOrdering::Amd).map_err(|e|e.to_string())?
        .factor(k,&LdltOptions::default()).map_err(|e|e.to_string())?;
    if factor.inertia().negative!=0 {return Err("static shell equilibrium is unstable (negative tangent inertia)".into());}
    Ok(factor)
}

/// Solve full-coordinate nonlinear static equilibrium, then return its stable
/// tangent pencil for the existing modal solver. `loads` is 6*node_count:
/// Cartesian forces [N] followed by moments [N m] at each reference node.
/// Clamped/pinned loads at fixed coordinates become support reactions, not
/// a fictitious displacement. Inputs and caller-owned geometry never mutate.
///
/// This is an incremental small-vibration model about a solved dead-load
/// equilibrium. `solve` returns coordinates and iteration count; a claimed
/// convergence is NOT trusted without force balance, deformation admission
/// and stable inertia.
/// # Errors
/// Invalid geometry/budgets, missing supports, excessive deformation, exhausted
/// Newton/Krylov work, unresolved force balance or nonpositive tangent inertia.
#[allow(clippy::too_many_arguments)]
pub fn equilibrate_stiffened_shell(mesh:&ShellMesh,sections:&[PlateSection],
    boundary:&[usize],support:ShellSupport,beams:&[ShellBeam],loads:&[f64],options:PreloadOptions,
    solve: &mut impl FnMut(&PreloadProblem<'_>, Vec<f64>, usize)->Result<(Vec<f64>,usize),String>)
    ->Result<PreloadedShell,String> {
    equilibrate_tethered_shell(mesh,sections,boundary,support,beams,loads,&[],options,solve)
}

/// As [`equilibrate_stiffened_shell`], with optional direction-updating tensile
/// connections. The continuation fraction scales their whole potential, so
/// residual, analytic Jacobian and accepted stability tests describe ONE system.
/// Final K includes the connections' material AND geometric stiffness. They
/// add no mass. At most 4096 connections are admitted in this cold solve.
/// # Errors
/// The dead-load refusals plus invalid point ports, anchors or tensile laws.
#[allow(clippy::too_many_arguments)]
pub fn equilibrate_tethered_shell(mesh:&ShellMesh,sections:&[PlateSection],
    boundary:&[usize],support:ShellSupport,beams:&[ShellBeam],loads:&[f64],tethers:&[ShellTether],
    options:PreloadOptions,
    solve: &mut impl FnMut(&PreloadProblem<'_>, Vec<f64>, usize)->Result<(Vec<f64>,usize),String>)
    ->Result<PreloadedShell,String> {
    if mesh.nodes.len()>10_000 || mesh.tris.len()>40_000 || loads.len()!=6*mesh.nodes.len()
        || tethers.len()>4096
        || loads.iter().any(|v|!v.is_finite()) || support==ShellSupport::Free || boundary.is_empty()
        || !(1..=64).contains(&options.load_steps) || !(1..=128).contains(&options.max_iterations)
        || options.max_evaluations==0 || options.max_evaluations>100_000
        || !options.maximum_gradient.is_finite() || !(0.0..=0.25).contains(&options.maximum_gradient)
        || options.maximum_gradient==0. || !options.maximum_rotation_rad.is_finite()
        || !(0.0..=0.25).contains(&options.maximum_rotation_rad) || options.maximum_rotation_rad==0.
        || !options.relative_force_tolerance.is_finite()
        || !(1e-12..=1e-4).contains(&options.relative_force_tolerance) {
        return Err("invalid supported-shell preload, load cardinality or explicit budget".into());
    }
    let model=assemble_stiffened_shell(mesh,sections,boundary,support,beams).map_err(|e|e.to_string())?;
    let base=positive_factor(&model.k)?;
    let length=(0..3).map(|c| {
        let lo=mesh.nodes.iter().map(|p|p[c]).fold(f64::INFINITY,f64::min);
        let hi=mesh.nodes.iter().map(|p|p[c]).fold(f64::NEG_INFINITY,f64::max);hi-lo
    }).fold(0.0_f64,f64::max);
    if !length.is_finite() || length<=0. {return Err("invalid shell norm length".into());}
    let mut scale=vec![1.;model.free];let mut load=vec![0.;model.free];
    for (full,free) in model.dof_map.iter().enumerate() {if let Some(i)=free {
        scale[*i]=if full%6<3 {1.}else{length};load[*i]=loads[full];
    }}
    let elements=sections.iter().enumerate().map(|(e,s)|Element::new(mesh,&model,s,e))
        .collect::<Result<Vec<_>,_>>()?;
    let tethers=tethers.iter().map(|t|PreparedTether::new(t,mesh,&model)).collect::<Result<Vec<_>,_>>()?;
    let mut problem=PreloadProblem {model:&model,elements,tethers,base,scale,load,fraction:0.,options,evaluations:Cell::new(0)};
    let mut x=vec![0.;model.free];let mut iterations=0;let mut final_k=None;
    for step in 1..=options.load_steps {
        problem.fraction=step as f64/options.load_steps as f64;
        let (candidate,count)=solve(&problem,x,options.max_iterations)
            .map_err(|e|format!("shell preload increment {step}/{}: {e}",options.load_steps))?;
        if candidate.len()!=model.free || candidate.iter().any(|v|!v.is_finite())
            || count>options.max_iterations {return Err("invalid injected equilibrium solution".into());}
        x=candidate;iterations+=count;
        let u=problem.physical(&x);
        problem.check_balance(&u)?;
        let k=problem.tangent(&u)?;positive_factor(&k)?;final_k=Some(k);
    }
    let u=problem.physical(&x);let mut force=vec![0.;model.free];
    let energy=problem.force_energy(&u,&mut force).ok_or("invalid converged shell state")?;
    let residual_force_n=problem.check_balance(&u)?;
    let k=final_k.ok_or("no completed preload increment")?;
    let evaluations=problem.evaluations.get();
    let tether_responses=problem.tethers.iter().map(|t|t.response(&u)).collect::<Result<Vec<_>,_>>()?;
    let displacement=model.dof_map.iter().map(|i|i.map_or(0.,|i|u[i])).collect();
    drop(problem);
    Ok(PreloadedShell {model:ShellModel {k,..model},displacement,residual_force_n,
        norm_length_m:length,stored_energy_j:energy,tether_responses,iterations,evaluations})
}

#[cfg(test)]
mod tests {
    use super::*;
    // Independent sparse-direct Newton oracle for small deterministic fixtures.
    // The product adapter supplies the existing fs-solver globalized driver.
    fn reference_newton(p:&PreloadProblem<'_>,mut x:Vec<f64>,limit:usize)
        ->Result<(Vec<f64>,usize),String> {
        for iteration in 0..limit {
            if !p.reserve() {return Err("test solve budget exhausted".into());}
            let u=p.physical(&x);
            if p.check_balance(&u).is_ok() {return Ok((x,iteration));}
            let mut force=vec![0.;p.dimension()];
            p.force_energy(&u,&mut force).ok_or("test trial refused")?;
            for (f,load) in force.iter_mut().zip(p.applied_load(&u)?) {*f-=load;}
            let factor=positive_factor(&p.tangent(&u)?)?;let delta=factor.solve(&force);
            for ((x,d),s) in x.iter_mut().zip(delta).zip(&p.scale) {*x-=d*s;}
        }
        Err("test Newton did not converge".into())
    }
    fn fixture(crown:f64)->(ShellMesh,Vec<PlateSection>,Vec<usize>) {
        let mut nodes=Vec::new();let mut tris=Vec::new();let mut fixed=Vec::new();
        for j in 0..3 {for i in 0..3 {
            nodes.push([i as f64*0.25,j as f64*0.25,if i==1&&j==1 {crown}else{0.}]);
            if i!=1||j!=1 {fixed.push(3*j+i);}
        }}
        for j in 0..2 {for i in 0..2 {let a=3*j+i;tris.push([a,a+1,a+4]);tris.push([a,a+4,a+3]);}}
        let mesh=ShellMesh::new(nodes,tris).unwrap();
        let sections=mesh.tris.iter().map(|_|PlateSection::orthotropic_plane_stress_at_angle(
            9e9,0.6e9,0.3,0.4e9,0.004,430.,0.).unwrap()).collect();
        (mesh,sections,fixed)
    }
    #[test]
    fn correction_gradient_and_consistent_tangent_match_independent_differences() {
        let (mesh,sections,fixed)=fixture(0.008);
        let model=assemble_stiffened_shell(&mesh,&sections,&fixed,ShellSupport::Clamped,&[]).unwrap();
        let element=Element::new(&mesh,&model,&sections[0],0).unwrap();
        let mut u=vec![0.;model.free];u[0]=0.0002;u[1]=-0.0001;u[2]=0.001;
        let c=element.correction(&u);let eps=1e-7;
        for i in 0..9 {if let Some(d)=element.dofs[i] {
            let mut a=u.clone();let mut b=u.clone();a[d]+=eps;b[d]-=eps;
            let ca=element.correction(&a);let cb=element.correction(&b);
            let fd=(ca.energy-cb.energy)/(2.*eps);
            assert!((fd-c.force[i]).abs()<1e-6*(1.+c.force[i].abs()));
            for j in 0..9 {if element.dofs[j].is_some() {
                let fd=(ca.force[j]-cb.force[j])/(2.*eps);
                assert!((fd-c.tangent[9*j+i]).abs()<1e-6*(1.+c.tangent[9*j+i].abs()));
            }}
        }}
    }
    #[test]
    fn zero_load_preserves_the_original_mass_stiffness_and_geometry() {
        let (mesh,sections,fixed)=fixture(0.008);let original=mesh.clone();
        let base=assemble_stiffened_shell(&mesh,&sections,&fixed,ShellSupport::Clamped,&[]).unwrap();
        let loaded=equilibrate_stiffened_shell(&mesh,&sections,&fixed,ShellSupport::Clamped,&[],
            &vec![0.;6*mesh.nodes.len()],PreloadOptions::default(),&mut reference_newton).unwrap();
        assert_eq!(mesh,original);assert!(loaded.displacement.iter().all(|u|*u==0.));
        assert_eq!(loaded.residual_force_n,0.);assert_eq!(loaded.iterations,0);
        assert!(loaded.tether_responses.is_empty());
        for i in 0..base.free {
            assert_eq!(loaded.model.k.row(i),base.k.row(i));assert_eq!(loaded.model.m.row(i),base.m.row(i));
        }
    }
    #[test]
    fn loaded_flat_board_has_equilibrium_and_a_changed_small_signal_pencil() {
        let (mesh,sections,fixed)=fixture(0.);
        let base=assemble_stiffened_shell(&mesh,&sections,&fixed,ShellSupport::Clamped,&[]).unwrap();
        let mut load=vec![0.;54];load[6*4+2]=-10.;
        let loaded=equilibrate_stiffened_shell(&mesh,&sections,&fixed,ShellSupport::Clamped,&[],&load,
            PreloadOptions::default(),&mut reference_newton).unwrap();
        assert!(loaded.displacement[26]<0.);assert!(loaded.residual_force_n<1e-6);
        assert!(loaded.stored_energy_j>0.);assert!(loaded.iterations>0);
        let mut v=vec![0.;base.free];v[base.dof_map[26].unwrap()]=1.;
        let mut a=v.clone();let mut b=v.clone();base.k.spmv(&v,&mut a);loaded.model.k.spmv(&v,&mut b);
        assert!(b[2]>a[2]*1.0001,"stretching must change physical tangent, not just a logged load");
        for node in fixed {assert!(loaded.displacement[6*node..6*node+6].iter().all(|u|*u==0.));}
    }
    #[test]
    fn invalid_load_budget_and_excessive_deformation_do_not_return_a_pencil() {
        let (mesh,sections,fixed)=fixture(0.008);let load=vec![0.;54];
        let solve=|load:&[f64],options|equilibrate_stiffened_shell(&mesh,&sections,&fixed,
            ShellSupport::Clamped,&[],load,options,&mut reference_newton);
        assert!(solve(&load[..53],PreloadOptions::default()).is_err());
        let mut bad=load.clone();bad[26]=f64::NAN;assert!(solve(&bad,PreloadOptions::default()).is_err());
        assert!(solve(&load,PreloadOptions {load_steps:0,..PreloadOptions::default()}).is_err());
        bad[26]=-1e8;assert!(solve(&bad,PreloadOptions {max_evaluations:64,..PreloadOptions::default()}).is_err());
    }
    #[test]
    fn coupled_tethers_change_equilibrium_tension_and_the_actual_tangent() {
        let (mesh,sections,fixed)=fixture(0.);
        let connections:Vec<_>=[-0.25,0.75].into_iter().map(|x|ShellTether {
            nodes:mesh.tris[0],weights:[0.,0.,1.],arm_m:[0.;3],anchor_m:[x,0.25,-0.04],
            reference_tension_n:100.,axial_stiffness_n_per_m:20_000.,
        }).collect();
        let load=vec![0.;54];
        let coupled=equilibrate_tethered_shell(&mesh,&sections,&fixed,ShellSupport::Clamped,&[],
            &load,&connections,PreloadOptions::default(),&mut reference_newton).unwrap();
        assert!(coupled.displacement[26]<0.);assert!(coupled.residual_force_n<1e-6);
        assert_eq!(coupled.tether_responses.len(),2);
        for response in &coupled.tether_responses {
            assert!(response.tension_n>0. && response.tension_n<100.);
            assert!(response.length_change_m<0.);
        }
        // Freeze the converged loads, not the unloaded-direction loads. This
        // recovers the SAME structural equilibrium but not the tether tangent.
        let mut frozen=load;
        for response in &coupled.tether_responses {
            for c in 0..3 {frozen[24+c]+=response.force_n[c];}
        }
        let bare=equilibrate_stiffened_shell(&mesh,&sections,&fixed,ShellSupport::Clamped,&[],
            &frozen,PreloadOptions::default(),&mut reference_newton).unwrap();
        for (a,b) in bare.displacement.iter().zip(&coupled.displacement) {assert!((a-b).abs()<1e-9);}
        let z=coupled.model.dof_map[26].unwrap();
        assert!(coupled.model.k.get(z,z)>bare.model.k.get(z,z));
        let mut bad=connections;bad[0].anchor_m=mesh.nodes[4];
        assert!(equilibrate_tethered_shell(&mesh,&sections,&fixed,ShellSupport::Clamped,&[],
            &vec![0.;54],&bad,PreloadOptions::default(),&mut reference_newton).is_err());
    }
}
