//! Rigid cylindrical cavity: radial P1 Galerkin x azimuthal Fourier x axial cosine.
//!
//! Separation of -Laplacian gives integral(r f'g' + m^2 f g/r) dr =
//! k_r^2 integral(r f g) dr. The axis is regular (f(0)=0 for m>0);
//! the wall is natural Neumann, not a pressure-release Dirichlet drumhead.
//! fs-modal owns the generalized eigensolve. No Bessel approximation or
//! hand-authored acoustic frequency is introduced. Radial quadrature uses three
//! Gauss points; the basis must be refined independently of the head mesh.
//! The analytic reference is k_r R = zeros of J_m', not zeros of J_m.
//! https://dlmf.nist.gov/10.21 and https://dlmf.nist.gov/10.6
use super::{ImpactError,invalid};
use crate::vibroacoustic::{AcousticMedium,CavityModes};
use fs_exec::CancelGate;
use fs_math::det;

/// Geometry and a bounded, explicit spectral/discretization window.
#[derive(Debug,Clone,Copy)]
pub struct CylinderSpec {
    pub radius_m:f64,
    pub depth_m:f64,
    /// Uniform radial intervals (4..=96); refine this for eigenfrequency accuracy.
    pub radial_intervals:usize,
    /// Both sine and cosine members are kept for every m>0, up to this order.
    pub maximum_azimuthal_order:usize,
    /// Axial cosine index 0..=this value. End planes are z=0 and z=depth.
    pub maximum_axial_order:usize,
    /// Retain every mode in the declared tensor window below this frequency.
    /// This is not an assertion that missing radial/angular families are converged.
    pub maximum_frequency_hz:f64,
    /// Refuse an overfull window; never split a degenerate angular pair to fit.
    pub maximum_modes:usize,
    /// Relative algebraic radial eigensolve residual ceiling, not a mesh error.
    pub eigen_residual_tolerance:f64,
}
/// One retained pressure basis member; radial coefficients are immutable.
#[derive(Debug,Clone)]
pub struct CylinderMode {
    pub azimuthal_order:usize,
    pub radial_index:usize,
    pub axial_order:usize,
    pub sine:bool,
    pub omega_rad_s:f64,
    pub norm_m3:f64,
    pub relative_radial_residual:f64,
    radial:Vec<f64>,
}
/// Reusable cold geometric basis for interface integration and interior probes.
#[derive(Debug,Clone)]
pub struct CylindricalCavity {spec:CylinderSpec,medium:AcousticMedium,modes:Vec<CylinderMode>}

fn radial_pencil(intervals:usize,order:usize)->(Vec<f64>,Vec<f64>) {
    let skip=usize::from(order>0);let n=intervals+1-skip;
    let(mut k,mut mass)=(vec![0.0;n*n],vec![0.0;n*n]);
    let h=1.0/intervals as f64;let root=(3.0_f64/5.0).sqrt();
    for e in 0..intervals {
        for (p,w) in [(-root,5.0/9.0),(0.0,8.0/9.0),(root,5.0/9.0)] {
            let t=0.5*(p+1.0);let r=(e as f64+t)*h;let measure=0.5*h*w;
            let shape=[1.0-t,t];let derivative=[-1.0/h,1.0/h];
            for a in 0..2 {for b in 0..2 {
                if e+a<skip || e+b<skip {continue;}
                let i=e+a-skip;let j=e+b-skip;
                mass[i*n+j]+=measure*r*shape[a]*shape[b];
                k[i*n+j]+=measure*(r*derivative[a]*derivative[b]
                    +(order*order) as f64*shape[a]*shape[b]/r);
            }}
        }
    }
    (k,mass)
}
fn residual(k:&[f64],mass:&[f64],phi:&[f64],lambda:f64)->f64 {
    let n=phi.len();let mut error=0.0_f64;let mut scale=0.0_f64;
    let amplitude=phi.iter().fold(0.0_f64,|a,p|a.max(p.abs()));
    for i in 0..n {
        let mut value=0.0;let mut row=0.0;
        for j in 0..n {value+=(k[i*n+j]-lambda*mass[i*n+j])*phi[j];
            row+=k[i*n+j].abs()+lambda.abs()*mass[i*n+j].abs();}
        error=error.max(value.abs());scale=scale.max(row*amplitude);
    }
    error/scale.max(f64::MIN_POSITIVE)
}
impl CylindricalCavity {
    /// Assemble the declared circular cylinder. Rigid-wall basis modes do not
    /// imply rigid moving heads: interface coupling supplies their actual motion.
    /// The zero radial/axial mode is the analytically known constant, exactly.
    pub fn new(spec:CylinderSpec,medium:AcousticMedium,gate:&CancelGate)->Result<Self,ImpactError> {
        if gate.is_requested() {return Err(ImpactError::Cancelled);}
        if ![spec.radius_m,spec.depth_m,spec.maximum_frequency_hz,spec.eigen_residual_tolerance,
            medium.rho0,medium.c0].iter().all(|v|v.is_finite() && *v>0.0)
            || !(4..=96).contains(&spec.radial_intervals) || spec.maximum_azimuthal_order>8
            || spec.maximum_axial_order>16 || !(1..=256).contains(&spec.maximum_modes)
            || spec.eigen_residual_tolerance>1e-3
        {return Err(invalid("cylinder needs finite positive geometry and bounded spectral/eigen budgets"));}
        let cutoff=core::f64::consts::TAU*spec.maximum_frequency_hz;
        if !cutoff.is_finite() {return Err(invalid("cylinder frequency ceiling overflow"));}
        let mut modes=Vec::new();
        for order in 0..=spec.maximum_azimuthal_order {
            if gate.is_requested() {return Err(ImpactError::Cancelled);}
            let (k,mass)=radial_pencil(spec.radial_intervals,order);
            let skip=usize::from(order>0);let n=spec.radial_intervals+1-skip;
            let radial_modes=fs_modal::eigh_gen_dense(&k,&mass,n).map_err(|e|ImpactError::Owner(e.to_string()))?;
            for (radial_index,pair) in radial_modes.into_iter().enumerate() {
                let constant=order==0 && radial_index==0;
                let lambda=if constant {0.0}else{pair.lambda};
                let phi=if constant {vec![1.0;n]}else{pair.phi};
                let error=residual(&k,&mass,&phi,lambda);
                if !lambda.is_finite() || lambda<0.0 || !error.is_finite() || error>spec.eigen_residual_tolerance {
                    return Err(invalid("cylinder radial eigenpair failed its algebraic residual gate"));
                }
                if medium.c0*lambda.sqrt()/spec.radius_m>cutoff {continue;}
                let peak=phi.iter().fold(0.0_f64,|a,p|a.max(p.abs()));
                if !peak.is_finite() || peak<=0.0 {return Err(invalid("cylinder radial shape has no finite amplitude"));}
                let sign=if phi[n-1]<0.0 {-1.0}else{1.0};
                let mut radial=vec![0.0;spec.radial_intervals+1];
                for i in 0..n {radial[i+skip]=sign*phi[i]/peak;}
                // Physical norm from the SAME radial mass matrix as the solve.
                let mut norm=0.0;
                for i in 0..n {for j in 0..n {norm+=radial[i+skip]*mass[i*n+j]*radial[j+skip];}}
                // Preserve the analytic constant exactly, including the volume.
                if constant {norm=0.5;radial.fill(1.0);}
                for axial in 0..=spec.maximum_axial_order {
                    let kz=core::f64::consts::PI*axial as f64/spec.depth_m;
                    let omega=medium.c0*(lambda.sqrt()/spec.radius_m).hypot(kz);
                    if !omega.is_finite() {return Err(invalid("cylinder derived frequency overflow"));}
                    if omega>cutoff {continue;}
                    let angular=if order==0 {core::f64::consts::TAU}else{core::f64::consts::PI};
                    let depth=spec.depth_m/if axial==0 {1.0}else{2.0};
                    let norm_m3=norm*spec.radius_m*spec.radius_m*angular*depth;
                    if !norm_m3.is_finite() || norm_m3<=0.0 {return Err(invalid("cylinder pressure norm overflow"));}
                    let multiplicity=if order==0 {1}else{2};
                    if modes.len()+multiplicity>spec.maximum_modes {return Err(invalid("cylinder frequency window exceeds complete-pair mode budget"));}
                    for member in 0..multiplicity {modes.push(CylinderMode {azimuthal_order:order,
                        radial_index,axial_order:axial,sine:member==1,omega_rad_s:omega,norm_m3,
                        relative_radial_residual:error,radial:radial.clone()});}
                }
            }
        }
        modes.sort_by(|a,b|a.omega_rad_s.total_cmp(&b.omega_rad_s)
            .then(a.azimuthal_order.cmp(&b.azimuthal_order)).then(a.radial_index.cmp(&b.radial_index))
            .then(a.axial_order.cmp(&b.axial_order)).then(a.sine.cmp(&b.sine)));
        if gate.is_requested() {return Err(ImpactError::Cancelled);}
        Ok(Self {spec,medium,modes})
    }
    #[must_use]
    pub fn modes(&self)->&[CylinderMode] {&self.modes}
    #[must_use]
    pub fn spec(&self)->CylinderSpec {self.spec}

    /// Dimensionless pressure mode values inside the cylinder, in retained order.
    /// Geometry origin is the centre of the z=0 end plane, not the cavity centre.
    pub fn values_at(&self,point:[f64;3])->Result<Vec<f64>,ImpactError> {
        let r=point[0].hypot(point[1]);
        if point.iter().any(|x|!x.is_finite()) || !r.is_finite() || r>self.spec.radius_m*(1.0+32.0*f64::EPSILON)
            || point[2]<0.0 || point[2]>self.spec.depth_m {
            return Err(invalid("cavity observation point is outside the declared cylinder"));
        }
        let station=(r/self.spec.radius_m).min(1.0)*self.spec.radial_intervals as f64;
        let cell=(station as usize).min(self.spec.radial_intervals-1);let t=station-cell as f64;
        let (cosine,sine)=if r==0.0 {(1.0,0.0)}else{(point[0]/r,point[1]/r)};
        let mut values=Vec::with_capacity(self.modes.len());
        for mode in &self.modes {
            let(mut c,mut s)=(1.0,0.0);
            for _ in 0..mode.azimuthal_order {(c,s)=(c*cosine-s*sine,s*cosine+c*sine);}
            let radial=(1.0-t)*mode.radial[cell]+t*mode.radial[cell+1];
            let axial=det::cos(core::f64::consts::PI*mode.axial_order as f64*point[2]/self.spec.depth_m);
            let value=radial*axial*if mode.sine {s}else{c};
            if !value.is_finite() {return Err(invalid("cylinder shape evaluation overflow"));}values.push(value);
        }
        Ok(values)
    }
    /// Sample the existing CavityModes carrier on any admitted interface points.
    /// Damping is not inferred from geometry. The caller supplies it to the time
    /// adapter or explicitly sets frequency-domain losses on its separate image.
    pub fn sample(&self,points:&[[f64;3]],maximum_terms:usize)->Result<CavityModes,ImpactError> {
        if points.is_empty() || self.modes.len().checked_mul(points.len()).is_none_or(|n|n>maximum_terms) {
            return Err(invalid("cylinder interface sample budget exceeded"));
        }
        let mut interface=vec![vec![0.0;points.len()];self.modes.len()];
        for (i,&point) in points.iter().enumerate() {for (j,value) in self.values_at(point)?.into_iter().enumerate() {
            interface[j][i]=value;
        }}
        Ok(CavityModes {omegas:self.modes.iter().map(|m|m.omega_rad_s).collect(),
            lambdas:self.modes.iter().map(|m|m.norm_m3).collect(),interface,loss_factor:0.0,
            rho0:self.medium.rho0,c0:self.medium.c0})
    }
}
