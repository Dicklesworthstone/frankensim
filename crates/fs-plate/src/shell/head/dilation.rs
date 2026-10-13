//! Uniform edge actuation of the existing homogeneous nonlinear film.
//!
//! q retains the original mass-normalized transverse basis; the final coordinate
//! Q=sqrt(M)*s represents outward edge dilation s [m]. Incremental uniform
//! tension is E*h/(1-nu)*s/R. With G_ij=integral grad(phi_i).grad(phi_j) dA,
//! H=H_membrane(q) + (E*h/(1-nu))*s/(2R)*q^T G q + (K_film+K_hardware)*s^2/2.
//! K_film=2*E*h/(1-nu)*area/R^2 uses the ACTUAL polygonal area.
//!
//! This is the SAME positive facet energy with affine strain s/R*[1,1,0]
//! added before squaring, not a Berger replacement for the original nonlinear
//! reduction. Its cross term reduces exactly by integration of div(u): relaxed
//! in-plane displacement vanishes on the complete rim, so its area integral is
//! zero. Uniform affine dilation needs no second static solve. Mixed derivatives
//! retain both changing head stress and the reciprocal actuator reaction.
//!
//! Geometry/acoustic charts remain fixed-reference; only small edge strain is
//! appropriate. M is the explicit TOTAL effective inertia of this coordinate.
//! Quasi-static transverse-induced in-plane relaxation remains unchanged. Static
//! installed preload is balanced at s=0; forces/work here are increments about
//! that equilibrium, not total pedal work. No particular pedal/hoop geometry,
//! moving acoustic boundary or measured instrument is inferred.
use super::{TensionedDisk, nonlinear::{MembraneReduction,MembraneReductionBudget}};
use crate::{ModePair,PlateError};
fn bad(what:&'static str)->PlateError {PlateError::BadSection {what}}

#[derive(Debug,Clone)]
pub struct DilationReduction {
    membrane:MembraneReduction,
    geometric:Vec<f64>,
    inverse_sqrt_mass:f64,
    radius_m:f64,
    tension_per_strain:f64,
    initial_tension:f64,
    tension_per_coordinate:f64,
    rim_per_mass:f64,
    rim_omega:f64,
}
impl DilationReduction {
    /// Reduce the actual uniform disk with its complete in-plane rim constrained.
    /// `effective_mass_kg` includes all represented inertia; hardware stiffness
    /// is ADDITIONAL to the film's own radial stiffness. The original cold solve
    /// callback is invoked once, as in MembraneReduction::new.
    #[allow(clippy::too_many_arguments)]
    pub fn new<F>(film:&TensionedDisk,modes:&[ModePair],effective_mass_kg:f64,
        hardware_stiffness_n_m:f64,budget:MembraneReductionBudget,solve:F)->Result<Self,PlateError>
    where F:FnOnce(usize,&[f64],&[Vec<f64>])->Result<Vec<Vec<f64>>,PlateError> {
        let spec=film.spec;
        if !film.tension_variation.is_zero() || !effective_mass_kg.is_finite() || effective_mass_kg<=0.0
            || !hardware_stiffness_n_m.is_finite() || hardware_stiffness_n_m<0.0
            || !spec.radius_m.is_finite() || spec.radius_m<=0.0
            || !spec.tension_n_m.is_finite() || spec.tension_n_m<=0.0 {
            return Err(bad("edge dilation needs a uniform film, positive total effective mass and passive extra stiffness"));
        }
        // Public disk fields must not substitute a different section for the
        // material that identifies the uniform biaxial tension increment.
        let section=crate::PlateSection::isotropic(spec.young_pa,spec.poisson,
            spec.thickness_m,spec.density_kg_m3)?;
        if section.d!=film.section.d || section.thickness!=film.section.thickness {
            return Err(bad("edge dilation requires the disk's original isotropic section"));
        }
        let mut uses=std::collections::BTreeMap::<(usize,usize),usize>::new();
        for tri in &film.mesh.tris {for i in 0..3 {
            let (a,b)=(tri[i],tri[(i+1)%3]);*uses.entry((a.min(b),a.max(b))).or_default()+=1;
        }}
        let mut rim:Vec<_>=uses.iter().filter(|(_,n)|**n==1).flat_map(|(&(a,b),_)|[a,b]).collect();
        rim.sort_unstable();rim.dedup();
        if rim.len()<3 || film.mesh.nodes.len().checked_mul(3)!=Some(film.model.dof_map.len())
            || uses.values().any(|n|*n>2)
            || rim.iter().any(|&i|i>=film.mesh.nodes.len()
                || (film.mesh.nodes[i].0.hypot(film.mesh.nodes[i].1)-spec.radius_m).abs()>1e-10*spec.radius_m
                || film.model.dof_map[3*i].is_some()) {
            return Err(bad("edge dilation requires the complete fixed-transverse circular film rim"));
        }
        let membrane=MembraneReduction::new(&film.mesh,&film.section,&film.model,modes,&rim,budget,solve)?;
        let inverse_sqrt_mass=1.0/effective_mass_kg.sqrt();
        let tension_per_strain=spec.young_pa*spec.thickness_m/(1.0-spec.poisson);
        let tension_per_coordinate=tension_per_strain/spec.radius_m*inverse_sqrt_mass;
        let n=modes.len();let mut geometric=vec![0.0;n*n];let mut area=0.0;
        for t in &film.mesh.tris {
            let p=t.map(|i|film.mesh.nodes[i]);
            let a2=(p[1].0-p[0].0)*(p[2].1-p[0].1)-(p[2].0-p[0].0)*(p[1].1-p[0].1);
            let grad:[[f64;2];3]=core::array::from_fn(|i|[
                (p[(i+1)%3].1-p[(i+2)%3].1)/a2,(p[(i+2)%3].0-p[(i+1)%3].0)/a2]);
            let g:Vec<[f64;2]>=modes.iter().map(|m|core::array::from_fn(|c|(0..3).map(|i|
                grad[i][c]*film.model.dof_map[3*t[i]].map_or(0.0,|j|m.phi[j])).sum())).collect();
            area+=0.5*a2;
            for i in 0..n {for j in 0..=i {
                let v=0.5*a2*(g[i][0]*g[j][0]+g[i][1]*g[j][1]);
                geometric[i*n+j]+=v;if i!=j {geometric[j*n+i]+=v;}
            }}
        }
        let rim_per_mass=(2.0*tension_per_strain*area/spec.radius_m.powi(2)+hardware_stiffness_n_m)/effective_mass_kg;
        let rim_omega=rim_per_mass.sqrt();
        if ![inverse_sqrt_mass,tension_per_strain,tension_per_coordinate,rim_per_mass,rim_omega]
            .iter().all(|v|v.is_finite() && *v>0.0) || geometric.iter().any(|v|!v.is_finite()) {
            return Err(bad("edge dilation coefficients are not representable"));
        }
        Ok(Self {membrane,geometric,inverse_sqrt_mass,radius_m:spec.radius_m,
            tension_per_strain,initial_tension:spec.tension_n_m,tension_per_coordinate,rim_per_mass,rim_omega})
    }
    pub fn membrane(&self)->&MembraneReduction {&self.membrane}
    pub fn mode_count(&self)->usize {self.membrane.mode_count()+1}
    pub fn rim_coordinate(&self)->usize {self.membrane.mode_count()}
    /// Physical radial force [N] times this weight is the canonical input effort.
    pub fn force_weight(&self)->f64 {self.inverse_sqrt_mass}
    pub fn rim_displacement_m(&self,q:&[f64])->f64 {
        if q.len()!=self.mode_count() {return f64::NAN;}q[self.rim_coordinate()]*self.inverse_sqrt_mass
    }
    pub fn edge_strain(&self,q:&[f64])->f64 {self.rim_displacement_m(q)/self.radius_m}
    /// Uniform prestress, EXCLUDING local vibration-induced stresses.
    pub fn baseline_tension_n_m(&self,q:&[f64])->f64 {self.initial_tension+self.tension_per_strain*self.edge_strain(q)}
    /// Bound the tunable LINEAR head pencil over an explicit strain cap. This
    /// is not a bound on all finite-amplitude frequencies or integration error.
    pub fn clock_bounds(&self,maximum_abs_strain:f64)->Result<Vec<f64>,PlateError> {
        if !maximum_abs_strain.is_finite() || maximum_abs_strain<=0.0 || maximum_abs_strain>=0.05 {
            return Err(bad("edge dilation needs an explicit small strain cap in (0,.05)"));
        }
        let ratio=(1.0+self.tension_per_strain*maximum_abs_strain/self.initial_tension).sqrt();
        let mut out:Vec<_>=self.membrane.omegas().iter().map(|w|w*ratio).collect();out.push(self.rim_omega);
        if out.iter().any(|w|!w.is_finite()) {return Err(bad("tunable film clock bound overflow"));}Ok(out)
    }
    fn gq(&self,q:&[f64],i:usize)->f64 {
        let n=self.rim_coordinate();(0..n).map(|j|self.geometric[i*n+j]*q[j]).sum()
    }
    fn quadratic(&self,q:&[f64])->f64 {(0..self.rim_coordinate()).map(|i|q[i]*self.gq(q,i)).sum()}
    /// Original nonlinear membrane energy plus the computed actuator coupling.
    /// Trial states outside positive-tension validity may return negative values;
    /// the time host must admit strain, tension and energy before publication.
    pub fn potential(&self,q:&[f64])->f64 {
        let n=self.rim_coordinate();
        if q.len()!=n+1 || q.iter().any(|v|!v.is_finite()) {return f64::NAN;}
        self.membrane.potential(&q[..n])+0.5*self.tension_per_coordinate*q[n]*self.quadratic(q)
            +0.5*self.rim_per_mass*q[n]*q[n]
    }
    pub fn gradient(&self,q:&[f64],out:&mut[f64]) {
        let n=self.rim_coordinate();
        if q.len()!=n+1 || out.len()!=n+1 || q.iter().any(|v|!v.is_finite()) {out.fill(f64::NAN);return;}
        self.membrane.gradient(&q[..n],&mut out[..n]);let scale=self.tension_per_coordinate;
        for i in 0..n {out[i]+=scale*q[n]*self.gq(q,i);}
        out[n]=0.5*scale*self.quadratic(q)+self.rim_per_mass*q[n];
    }
    pub fn hessian_vector(&self,q:&[f64],d:&[f64],out:&mut[f64]) {
        let n=self.rim_coordinate();
        if q.len()!=n+1 || d.len()!=n+1 || out.len()!=n+1 || q.iter().chain(d).any(|v|!v.is_finite()) {out.fill(f64::NAN);return;}
        self.membrane.hessian_vector(&q[..n],&d[..n],&mut out[..n]);out[n]=self.rim_per_mass*d[n];
        for i in 0..n {
            let cross=self.tension_per_coordinate*self.gq(q,i);
            out[i]+=self.tension_per_coordinate*q[n]*self.gq(d,i)+cross*d[n];out[n]+=cross*d[i];
        }
    }
    /// Exact symmetric tangent. Reuses the original batched membrane tangent,
    /// plus O(modes^2) work; no modal re-solve or hot-path allocation.
    pub fn hessian_matrix(&self,q:&[f64],out:&mut[f64],scratch:&mut[[f64;3]]) {
        let n=self.rim_coordinate();let size=n+1;
        if q.len()!=size || out.len()!=size*size || scratch.len()!=size || q.iter().any(|v|!v.is_finite()) {out.fill(f64::NAN);return;}
        self.membrane.hessian_matrix(&q[..n],&mut out[..n*n],&mut scratch[..n]);
        // Expand the original row stride backwards, preserving unconsumed rows.
        for i in (0..n).rev() {out.copy_within(i*n..(i+1)*n,i*size);}
        for i in 0..n {
            for j in 0..n {out[i*size+j]+=self.tension_per_coordinate*q[n]*self.geometric[i*n+j];}
            let cross=self.tension_per_coordinate*self.gq(q,i);out[i*size+n]=cross;out[n*size+i]=cross;
        }
        out[n*size+n]=self.rim_per_mass;
    }
}
#[cfg(test)]
mod tests;
