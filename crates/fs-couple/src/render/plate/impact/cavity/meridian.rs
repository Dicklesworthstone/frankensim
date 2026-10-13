//! Axisymmetric enclosed air from a supplied radius/height profile.
//!
//! Continuous P1 functions in (r,z), times cos(m theta)/sin(m theta), solve
//! integral [r grad f.grad g + m^2 f*g/r] dr dz = k^2 integral r*f*g dr dz.
//! Pressure is natural Neumann on the profile and end disks. For m>0, axis
//! vertices are fixed to zero; for m=0 the constant mode is retained exactly.
//! Four-by-four positive Gauss/Duffy quadrature supplies the angular term;
//! the same mass matrix supplies norms. fs-modal owns the generalized solve.
//! No second time integrator or cylinder-equivalent volume approximation.
//!
//! This is a radial-graph chart, not arbitrary 3-D geometry. Interior stations
//! have positive radius; either end may close on axis. Profile, meridional mesh,
//! angular families and integration need independent convergence. A rigid-wall
//! pressure basis can couple to moving surfaces through the existing overlaps.
use super::{ImpactError, invalid};
use crate::vibroacoustic::{AcousticMedium, CavityModes};
use fs_exec::CancelGate;
use core::f64::consts::{PI, TAU};

/// One physical profile station. Straight segments between stations are the
/// actual supplied geometry, not samples silently fitted to a hemisphere.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Station {
    /// Axial coordinate [m], strictly increasing in input order.
    pub z_m: f64,
    /// Nonnegative radial extent [m]. Zero is permitted only at an end.
    pub radius_m: f64,
}
/// Independent cold mesh, spectral and memory ceilings.
#[derive(Clone, Copy, Debug)]
pub struct Budget {
    /// Radial intervals at each positive-radius station (2..=32).
    pub radial_intervals: usize,
    /// Inclusive angular Fourier order (0..=8); complete pairs only.
    pub maximum_azimuthal_order: usize,
    /// Retain every mode below this ceiling in the declared angular space [Hz].
    pub maximum_frequency_hz: f64,
    /// Refuse, never truncate, a larger retained set (1..=256).
    pub maximum_modes: usize,
    /// P1 nodes, before removing axis vertices (at most 512).
    pub maximum_nodes: usize,
    /// Bound the sum of squared free dimensions across all angular pencils.
    pub maximum_matrix_entries: usize,
    /// Algebraic eigensolve tolerance, not a mesh-error certificate.
    pub relative_tolerance: f64,
}
/// One retained mode with its immutable nodal pressure coefficients.
#[derive(Clone, Debug)]
pub struct Mode {
    pub azimuthal_order: usize,
    pub family_index: usize,
    pub sine: bool,
    pub omega_rad_s: f64,
    pub norm_m3: f64,
    pub relative_residual: f64,
    coefficients: Vec<f64>,
}
/// One pressure basis and geometry used for volume, overlaps and point samples.
#[derive(Clone, Debug)]
pub struct MeridianCavity {
    profile: Vec<Station>,
    nodes: Vec<[f64;2]>,
    triangles: Vec<[usize;3]>,
    modes: Vec<Mode>,
    medium: AcousticMedium,
    scale: f64,
    volume_m3: f64,
}
const X: [f64;4] = [0.06943184420297371,0.33000947820757187,0.6699905217924281,0.9305681557970262];
const W: [f64;4] = [0.17392742256872692,0.32607257743127307,0.32607257743127307,0.17392742256872692];
fn cross(a:[f64;2],b:[f64;2])->f64 {a[0]*b[1]-a[1]*b[0]}
fn sub(a:[f64;2],b:[f64;2])->[f64;2] {[a[0]-b[0],a[1]-b[1]]}
fn cancelled(gate:&CancelGate)->Result<(),ImpactError> {
    if gate.is_requested() {Err(ImpactError::Cancelled)} else {Ok(())}
}

impl MeridianCavity {
    /// Assemble a connected profile-defined chamber. End disks are pressure
    /// boundaries, including a moving head when connected by surface overlaps.
    /// No initial pressure, material damping or external radiation is inferred.
    pub fn new(profile:Vec<Station>,medium:AcousticMedium,budget:Budget,gate:&CancelGate)
        ->Result<Self,ImpactError> {
        cancelled(gate)?;
        if !(2..=128).contains(&profile.len()) || !(2..=32).contains(&budget.radial_intervals)
            || !(1..=512).contains(&budget.maximum_nodes) || !(1..=256).contains(&budget.maximum_modes)
            || budget.maximum_azimuthal_order>8 || budget.maximum_matrix_entries==0
            || ![medium.rho0,medium.c0,budget.maximum_frequency_hz,budget.relative_tolerance]
                .iter().all(|v|v.is_finite() && *v>0.0) || budget.relative_tolerance>1e-3
            || profile.iter().enumerate().any(|(i,s)| !s.z_m.is_finite() || !s.radius_m.is_finite()
                || s.radius_m<0.0 || (i>0 && i+1<profile.len() && s.radius_m==0.0))
            || profile.windows(2).any(|s|s[1].z_m<=s[0].z_m || (s[0].radius_m==0.0 && s[1].radius_m==0.0)) {
            return Err(invalid("meridian needs ordered finite geometry and bounded mesh/spectral controls"));
        }
        let scale=profile.iter().map(|s|s.radius_m).fold(0.0_f64,f64::max);
        let nr=budget.radial_intervals;
        let count=profile.iter().map(|s|if s.radius_m==0.0 {1}else{nr+1}).sum::<usize>();
        let free=count-profile.len();
        let entries=count*count+budget.maximum_azimuthal_order*free*free;
        let cutoff=TAU*budget.maximum_frequency_hz;
        let mut volume=0.0;
        for s in profile.windows(2) {
            let a=s[0].radius_m;let b=s[1].radius_m;
            volume+=PI/3.0*(s[1].z_m-s[0].z_m)*(a*a+a*b+b*b);
        }
        if count>budget.maximum_nodes || entries>budget.maximum_matrix_entries
            || !scale.is_finite() || scale<=0.0 || !volume.is_finite() || volume<=0.0
            || !cutoff.is_finite() || !(scale*scale*scale).is_finite() || scale*scale*scale==0.0 {
            return Err(invalid("meridian exceeds geometry/pencil budget or representable volume"));
        }
        let mut nodes=Vec::with_capacity(count);let mut rows=Vec::with_capacity(profile.len());
        for s in &profile {
            let mut row=Vec::new();
            for i in 0..if s.radius_m==0.0 {1}else{nr+1} {
                row.push(nodes.len());nodes.push([s.radius_m/scale*i as f64/nr as f64,
                    (s.z_m-profile[0].z_m)/scale]);
            }
            rows.push(row);
        }
        if nodes.iter().flatten().any(|v|!v.is_finite()) {return Err(invalid("meridian normalized coordinate overflow"));}
        let mut triangles=Vec::new();
        for pair in rows.windows(2) {
            let (a,b)=(&pair[0],&pair[1]);
            for i in 0..nr {
                if a.len()==1 {triangles.push([a[0],b[i+1],b[i]]);}
                else if b.len()==1 {triangles.push([a[i],a[i+1],b[0]]);}
                else {triangles.extend([[a[i],a[i+1],b[i+1]],[a[i],b[i+1],b[i]]]);}
            }
        }
        let mut result=Self {profile,nodes,triangles,modes:Vec::new(),medium,scale,volume_m3:volume};
        for order in 0..=budget.maximum_azimuthal_order {
            cancelled(gate)?;
            let (k,m,indices)=result.pencil(order,gate)?;let n=indices.len();
            let pairs=fs_modal::eigh_gen_dense(&k,&m,n).map_err(|e|ImpactError::Owner(e.to_string()))?;
            for (family,pair) in pairs.into_iter().enumerate() {
                let constant=order==0 && family==0;
                let lambda=if constant {0.0}else{pair.lambda};
                if !lambda.is_finite() || lambda<0.0 {return Err(invalid("meridian eigenvalue outside finite positive spectrum"));}
                let omega=medium.c0*lambda.sqrt()/scale;
                if !omega.is_finite() {return Err(invalid("meridian frequency overflow"));}
                if omega>cutoff {continue;}
                let phi=if constant {vec![1.0;n]}else{pair.phi};
                let mut peak=0.0_f64;let mut sign=1.0;
                for &p in &phi {if p.abs()>peak {peak=p.abs();sign=if p<0.0 {-1.0}else{1.0};}}
                if !peak.is_finite() || peak<=0.0 {return Err(invalid("meridian eigenvector has no finite scale"));}
                let mut residual=0.0_f64;let mut denominator=0.0_f64;let mut norm=0.0;
                for i in 0..n {
                    let mut value=0.0;let mut row=0.0;
                    for j in 0..n {
                        value+=(k[i*n+j]-lambda*m[i*n+j])*phi[j];
                        row+=k[i*n+j].abs()+lambda*m[i*n+j].abs();
                        norm+=(phi[i]/peak)*m[i*n+j]*(phi[j]/peak);
                    }
                    residual=residual.max(value.abs());denominator=denominator.max(row*peak);
                }
                let error=residual/denominator.max(f64::MIN_POSITIVE);
                if !error.is_finite() || error>budget.relative_tolerance {
                    return Err(invalid("meridian eigenvector failed the algebraic residual gate"));
                }
                let norm=if constant {volume}else{norm*scale*scale*scale*if order==0 {TAU}else{PI}};
                if !norm.is_finite() || norm<=0.0 {return Err(invalid("meridian mode norm is unrepresentable"));}
                let multiplicity=if order==0 {1}else{2};
                if result.modes.len()+multiplicity>budget.maximum_modes {
                    return Err(invalid("meridian window exceeds the complete angular-pair mode budget"));
                }
                let mut coefficients=vec![0.0;count];
                for (&index,&p) in indices.iter().zip(&phi) {coefficients[index]=if constant {1.0}else{sign*p/peak};}
                for member in 0..multiplicity {result.modes.push(Mode {azimuthal_order:order,
                    family_index:family,sine:member==1,omega_rad_s:omega,norm_m3:norm,
                    relative_residual:error,coefficients:coefficients.clone()});}
            }
        }
        result.modes.sort_by(|a,b|a.omega_rad_s.total_cmp(&b.omega_rad_s)
            .then(a.azimuthal_order.cmp(&b.azimuthal_order)).then(a.family_index.cmp(&b.family_index))
            .then(a.sine.cmp(&b.sine)));
        cancelled(gate)?;Ok(result)
    }
    fn pencil(&self,order:usize,gate:&CancelGate)->Result<(Vec<f64>,Vec<f64>,Vec<usize>),ImpactError> {
        let indices:Vec<_>=(0..self.nodes.len()).filter(|&i|order==0 || self.nodes[i][0]>0.0).collect();
        let mut map=vec![None;self.nodes.len()];for (j,&i) in indices.iter().enumerate() {map[i]=Some(j);}
        let n=indices.len();let (mut k,mut m)=(vec![0.0;n*n],vec![0.0;n*n]);
        for &tri in &self.triangles {
            cancelled(gate)?;
            let p=tri.map(|i|self.nodes[i]);let area=cross(sub(p[1],p[0]),sub(p[2],p[0]));
            if !area.is_finite() || area<=0.0 {return Err(invalid("meridian has degenerate or reversed triangles"));}
            let grad:[[f64;2];3]=core::array::from_fn(|i| {
                let (a,b)=(p[(i+1)%3],p[(i+2)%3]);[(a[1]-b[1])/area,(b[0]-a[0])/area]
            });
            for (u,&wu) in X.iter().zip(&W) {for (v,&wv) in X.iter().zip(&W) {
                let shape=[(1.0-u)*(1.0-v),*u,(1.0-u)*v];
                let r=(0..3).map(|i|p[i][0]*shape[i]).sum::<f64>();
                let weight=area*(1.0-u)*wu*wv;
                if !r.is_finite() || r<=0.0 {return Err(invalid("meridian quadrature reached invalid axis radius"));}
                for i in 0..3 {if let Some(a)=map[tri[i]] {for j in 0..3 {if let Some(b)=map[tri[j]] {
                    m[a*n+b]+=weight*r*shape[i]*shape[j];
                    k[a*n+b]+=weight*(r*(grad[i][0]*grad[j][0]+grad[i][1]*grad[j][1])
                        +(order*order) as f64*shape[i]*shape[j]/r);
                }}}}
            }}
        }
        if k.iter().chain(&m).any(|v|!v.is_finite()) {return Err(invalid("nonfinite meridian pressure pencil"));}
        Ok((k,m,indices))
    }
    /// The original piecewise-linear physical wall, unchanged by discretization.
    pub fn profile(&self)->&[Station] {&self.profile}
    pub fn modes(&self)->&[Mode] {&self.modes}
    pub fn volume_m3(&self)->f64 {self.volume_m3}
    /// Interpolate the P1 field in its actual meridian cell, with exact angular
    /// factors. Outside points refuse, including points in the enclosing box
    /// but outside the bowl. This is cold receiver/overlap preparation.
    pub fn values_at(&self,point:[f64;3])->Result<Vec<f64>,ImpactError> {
        if point.iter().any(|v|!v.is_finite()) {return Err(invalid("nonfinite meridian receiver"));}
        let r=point[0].hypot(point[1]);let p=[r/self.scale,(point[2]-self.profile[0].z_m)/self.scale];
        if !r.is_finite() || p.iter().any(|v|!v.is_finite()) {return Err(invalid("meridian receiver scaling overflow"));}
        let mut found=None;
        for &tri in &self.triangles {
            let v=tri.map(|i|self.nodes[i]);let a=sub(v[1],v[0]);let b=sub(v[2],v[0]);let x=sub(p,v[0]);
            let d=cross(a,b);let u=cross(x,b)/d;let w=cross(a,x)/d;let bary=[1.0-u-w,u,w];
            if bary.iter().all(|b|*b>=-64.0*f64::EPSILON && *b<=1.0+64.0*f64::EPSILON) {found=Some((tri,bary));break;}
        }
        let (tri,bary)=found.ok_or_else(||invalid("point is outside the actual meridian cavity"))?;
        let (cosine,sine)=if r==0.0 {(1.0,0.0)}else{(point[0]/r,point[1]/r)};
        let mut out=Vec::with_capacity(self.modes.len());
        for mode in &self.modes {
            if mode.omega_rad_s==0.0 {out.push(1.0);continue;}
            let amplitude=(0..3).map(|i|bary[i]*mode.coefficients[tri[i]]).sum::<f64>();
            let (mut c,mut s)=(1.0,0.0);
            for _ in 0..mode.azimuthal_order {(c,s)=(c*cosine-s*sine,s*cosine+c*sine);}
            let value=amplitude*if mode.sine {s}else{c};
            if !value.is_finite() {return Err(invalid("meridian pressure shape overflow"));}out.push(value);
        }
        Ok(out)
    }
    /// Existing generic carrier for the SAME basis, volume norms and interface.
    pub fn sample(&self,points:&[[f64;3]],maximum_terms:usize,gate:&CancelGate)->Result<CavityModes,ImpactError> {
        cancelled(gate)?;
        if points.is_empty() || points.len().checked_mul(self.modes.len()).is_none_or(|n|n>maximum_terms) {
            return Err(invalid("meridian interface sampling exceeds its work budget"));
        }
        let mut interface=vec![vec![0.0;points.len()];self.modes.len()];
        for (i,&p) in points.iter().enumerate() {cancelled(gate)?;
            for (j,value) in self.values_at(p)?.into_iter().enumerate() {interface[j][i]=value;}
        }
        Ok(CavityModes {omegas:self.modes.iter().map(|m|m.omega_rad_s).collect(),
            lambdas:self.modes.iter().map(|m|m.norm_m3).collect(),interface,loss_factor:0.0,
            rho0:self.medium.rho0,c0:self.medium.c0})
    }
}

#[cfg(test)]
mod tests;
