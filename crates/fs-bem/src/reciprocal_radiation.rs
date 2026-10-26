//! Symmetric P1/P0 Galerkin radiation for reciprocal mechanical feedback.
//!
//! Unlike centroid collocation this is a variational, work-conjugate boundary
//! image. With outward normals, e^{-i wt}, pressure p in continuous P1 and an
//! auxiliary single-layer density lambda in P0, solve
//!     [ N       (K-M/2)^T ] [ p      ] = [ M^T i*w*rho*v ]
//!     [ K-M/2    V        ] [ lambda ]   [ 0             ].
//! Eliminating lambda gives the symmetric exterior DtN operator
//! N-(K-M/2)^T V^{-1}(K-M/2). The exterior field is D*p+S*lambda;
//! lambda is NOT the prescribed Neumann datum and is never labelled as velocity.
//!
//! Real hypersingular terms use the surface-curl weak form. Static singularities
//! have analytic triangle moments. The imaginary part of the ENTIRE Calderon
//! block is integrated as the positive plane-wave identity for Im G. Its Gram
//! structure gives the discrete outgoing-flux identity before any mechanical
//! projection. No eigenvalue clipping, response symmetrization, added resistance,
//! or replacement of a refused collocation result occurs.
//!
//! Reference: Chen, Hofstetter & Mang, "A Galerkin-type BE-FE formulation for
//! elasto-acoustic coupling", CMAME 152 (1998), 147-155,
//! doi:10.1016/S0045-7825(97)00187-4.
//!
//! This bounded dense image admits only a band below the closed body's slab
//! Dirichlet lower bound. It is NOT a resonance-safe Burton-Miller replacement.
//! P1/P0 spaces, 4x4 triangle quadrature and angular refinement must converge;
//! a passive discrete load is not an acoustic-accuracy or real-time certificate.
use crate::helmholtz::{HelmholtzError,Medium,MIN_PANELS_PER_WAVELENGTH,MAX_RADIATION_FIELDS_PER_BATCH};
use crate::panel3d::SpherePanels;
use fs_la::eigen_complex::{LuComplexError, lu_complex_with_cancel};
use fs_math::{c64::C64,det};
use core::f64::consts::{PI,TAU};
mod geometry;
use geometry::*;
/// Open P1 pressure-jump screens using the same singular triangle moments.
pub mod open_screen;
/// Reflection at the open end of a thin unflanged pipe, on the open screen.
pub mod pipe_end;
mod field;
pub use field::Solution;
mod ports;
pub use ports::VelocityPortSolution;
mod receivers;
pub use receivers::{ReceiverGeometry,PreparedReceivers};

/// Largest admitted P1 vertex plus P0 panel count. It also bounds each
/// frequency's dense system, 16 n^2 bytes beside the retained cache:
/// 1.07 GB at this ceiling.
pub const MAX_RECIPROCAL_UNKNOWNS:usize=8192;
/// The panel count whose closed genus-0 surface has MAX_RECIPROCAL_UNKNOWNS
/// unknowns (3F/2 + 2).
pub const MAX_RECIPROCAL_PANELS:usize=5460;
/// Largest angular refinement of the plane-wave identity for Im G. Each level
/// doubles the admitted k*diameter and quadruples the direction count.
pub const MAX_ANGULAR_REFINEMENT:u32=3;

/// Explicit cold storage/work ceilings; no mechanical modes are truncated.
#[derive(Debug,Clone,Copy)]
pub struct Budget {
    /// P1 vertices plus P0 panels, at most [`MAX_RECIPROCAL_UNKNOWNS`].
    pub maximum_unknowns:usize,
    /// Retained geometry/static moments, not allocator overhead or per-solve LU.
    pub maximum_cache_bytes:usize,
    /// Conservative scalar-work units for ONE frequency including its batch.
    /// Not measured flops, wall time, or a cumulative sweep budget.
    pub maximum_work:u64,
    /// 0..=[`MAX_ANGULAR_REFINEMENT`]: subdivides the 8-point polar rule and periodic azimuthal rule.
    /// Antipodal partners share weights exactly. Refine independently of mesh.
    pub angular_refinement:u32,
}

/// Borrowed immutable boundary: solutions cannot be evaluated on another mesh.
pub struct Prepared<'s> {
    surface:&'s SpherePanels, facets:Vec<Facet>, vertices:usize, length:f64,
    highest_k:f64, ppw_at_highest:f64, budget:Budget,
    pairs:Vec<Pair>, k0:Vec<f64>, bytes:usize,
}
fn bad(what:&'static str)->HelmholtzError {HelmholtzError::BadParameter{what}}
fn finite(v:C64)->bool {v.re.is_finite() && v.im.is_finite()}
fn poll(cancelled:&mut impl FnMut()->bool)->Result<(),HelmholtzError> {
    if cancelled(){Err(HelmholtzError::Cancelled)}else{Ok(())}
}
fn lu_error(error:LuComplexError)->HelmholtzError {
    match error {
        LuComplexError::Cancelled=>HelmholtzError::Cancelled,
        LuComplexError::Singular{..}=>HelmholtzError::Singular,
    }
}
fn zeroed<T:Clone>(n:usize,value:T)->Result<Vec<T>,HelmholtzError> {
    let mut v=Vec::new();v.try_reserve_exact(n).map_err(|_|HelmholtzError::AllocationFailed{what:"reciprocal acoustic workspace"})?;
    v.resize(n,value);Ok(v)
}
fn add_symmetric(h:&mut [C64],n:usize,i:usize,j:usize,value:f64) {
    let (a,b)=if i<=j {(i,j)}else{(j,i)};
    h[a*n+b].re+=value;h[b*n+a]=h[a*n+b];
}

impl<'s> Prepared<'s> {
    /// Prepare the same supplied closed triangles, without changing geometry.
    /// Cancellation is checked before allocation and between triangle pairs.
    /// # Errors
    /// Open/nonconforming surfaces, unresolved/high bands, exceeded budgets,
    /// nonfinite moments, allocation failure or cooperative cancellation.
    pub fn new(surface:&'s SpherePanels,highest_k:f64,budget:Budget,
        mut cancelled:impl FnMut()->bool)->Result<Self,HelmholtzError> {
        poll(&mut cancelled)?;
        let nf=surface.areas().len();
        if !(4..=MAX_RECIPROCAL_PANELS).contains(&nf) || !(1..=MAX_RECIPROCAL_UNKNOWNS).contains(&budget.maximum_unknowns)
            || budget.angular_refinement>MAX_ANGULAR_REFINEMENT || budget.maximum_work==0
            || !highest_k.is_finite() || highest_k<=0.0 {
            return Err(bad("invalid reciprocal geometry/frequency/work envelope"));
        }
        // Same area-based wavelength floor as the collocation owner.
        let ppw=TAU/highest_k/surface.areas().iter().copied().fold(0.0_f64,f64::max).sqrt();
        if !ppw.is_finite() || ppw<MIN_PANELS_PER_WAVELENGTH {
            return Err(HelmholtzError::TooCoarse{panels_per_wavelength:ppw});
        }
        let (facets,vertices,widths,length)=facets(surface)?;
        let n=vertices+nf;
        if n>budget.maximum_unknowns {return Err(HelmholtzError::WorkCap{panels:n});}
        // Reuse the geometry owner's component/orientation admission. This
        // point is beyond the bounding box; no observer is retained or solved.
        let mut outside=facets[0].triangle[0];outside[0]+=2.0*length;
        let _closed=crate::near_field::Geometry::new(surface,&[outside],64.0*f64::EPSILON*length)?;
        let slab=widths.into_iter().fold(f64::INFINITY,f64::min);
        if slab<=0.0 || (highest_k*slab).next_up()>=PI.next_down()
            || highest_k*length>4.0*f64::from(1_u32<<budget.angular_refinement) {
            return Err(bad("reciprocal band exceeds its resonance-exclusion or angular-resolution chart"));
        }
        let pair_count=nf*(nf+1)/2;
        let bytes=pair_count*core::mem::size_of::<Pair>()+nf*vertices*8
            +nf*core::mem::size_of::<Facet>();
        if bytes>budget.maximum_cache_bytes {
            return Err(HelmholtzError::CacheBudget{required_bytes:bytes,maximum_bytes:budget.maximum_cache_bytes});
        }
        let mut pairs=Vec::new();pairs.try_reserve_exact(pair_count).map_err(|_|
            HelmholtzError::AllocationFailed{what:"reciprocal static triangle pairs"})?;
        let mut k0=zeroed(nf*vertices,0.0)?;
        // Preflight frequency work before integrating a static table. A later
        // batch repeats this check with its actual source count.
        let mut out=Self {surface,facets,vertices,length,highest_k,ppw_at_highest:ppw,
            budget,pairs:Vec::new(),k0:Vec::new(),bytes};
        out.admit_work(1)?;
        for i in 0..nf {
            for a in out.facets[i].ids {k0[i*vertices+a]-=out.facets[i].area/6.0;}
            for j in i..nf {
                poll(&mut cancelled)?;
                let (pair,forward,back)=pair(&out.facets[i],&out.facets[j],i==j);
                if !pair.single.is_finite() || pair.shape.iter().flatten().chain(&forward).chain(&back).any(|v|!v.is_finite()) {
                    return Err(bad("nonfinite reciprocal static triangle moments"));
                }
                for a in 0..3 {
                    k0[i*vertices+out.facets[j].ids[a]]+=forward[a];
                    if i!=j {k0[j*vertices+out.facets[i].ids[a]]+=back[a];}
                }
                pairs.push(pair);
            }
        }
        out.pairs=pairs;out.k0=k0;poll(&mut cancelled)?;Ok(out)
    }
    /// P1 plus P0 unknown count, not the number of mechanical sources.
    pub fn unknowns(&self)->usize {self.vertices+self.facets.len()}
    /// Logical retained bytes; LU/solutions allocate additional bounded storage.
    pub fn cache_bytes(&self)->usize {self.bytes}
    /// Conservative assembly/factorization units for one frequency.
    pub fn frequency_work(&self,fields:usize)->Result<u64,HelmholtzError> {
        if !(1..=MAX_RADIATION_FIELDS_PER_BATCH).contains(&fields) {return Err(bad("reciprocal source-count budget"));}
        let n=self.unknowns() as u128;let f=self.facets.len() as u128;
        let directions=64_u128*(1_u128<<(2*self.budget.angular_refinement));
        // Includes singular-moment preparation conservatively on every call,
        // tensor regular-kernel work, angular Gram work, LU and all RHS solves.
        let work=n*n*n+2*n*n*fields as u128+f*f*16*16*96+directions*(n*n+f*16*24);
        u64::try_from(work).map_err(|_|bad("reciprocal work overflow"))
    }
    fn admit_work(&self,fields:usize)->Result<(),HelmholtzError> {
        if self.frequency_work(fields)?>self.budget.maximum_work {Err(bad("reciprocal frequency exceeds declared work budget"))}else{Ok(())}
    }

    fn assemble(&self,k:f64,cancelled:&mut impl FnMut()->bool)->Result<(Vec<C64>,Vec<Wave>),HelmholtzError> {
        let nf=self.facets.len();let nv=self.vertices;let n=self.unknowns();
        let mut h=zeroed(n*n,C64::ZERO)?;let mut at=0;
        for i in 0..nf {
            poll(cancelled)?;
            for a in 0..nv {
                h[(nv+i)*n+a]=C64::from_re(self.k0[i*nv+a]);h[a*n+nv+i]=h[(nv+i)*n+a];
            }
        }
        for i in 0..nf {for j in i..nf {
            poll(cancelled)?;
            let a=&self.facets[i];let b=&self.facets[j];
            let mut pair=self.pairs[at];at+=1;
            let mut kij=[0.0;3];let mut kji=[0.0;3];
            for x in a.samples {for y in b.samples {
                let d=sub(x.point,y.point);let r=norm(d);let weight=x.weight*y.weight;
                let (g,derivative)=regular_real(k,r);
                pair.single+=weight*g;
                for u in 0..3 {for v in 0..3 {pair.shape[u][v]+=weight*x.shape[u]*y.shape[v]*g;}}
                if i!=j {for u in 0..3 {
                    kij[u]+=weight*derivative*dot(b.normal,d)*y.shape[u];
                    kji[u]-=weight*derivative*dot(a.normal,d)*x.shape[u];
                }}
            }}
            h[(nv+i)*n+nv+j]=C64::from_re(pair.single);h[(nv+j)*n+nv+i]=h[(nv+i)*n+nv+j];
            for u in 0..3 {
                h[(nv+i)*n+b.ids[u]].re+=kij[u];h[b.ids[u]*n+nv+i]=h[(nv+i)*n+b.ids[u]];
                if i!=j {h[(nv+j)*n+a.ids[u]].re+=kji[u];h[a.ids[u]*n+nv+j]=h[(nv+j)*n+a.ids[u]];}
                for v in 0..3 {
                    if i==j && v<u {continue;}
                    let moment=if i==j {0.5*pair.shape[u][v]+0.5*pair.shape[v][u]}else{pair.shape[u][v]};
                    let value=k*k*dot(a.normal,b.normal)*moment-dot(a.curl[u],b.curl[v])*pair.single;
                    let multiplicity=if i!=j && a.ids[u]==b.ids[v] {2.0}else{1.0};
                    add_symmetric(&mut h,n,a.ids[u],b.ids[v],multiplicity*value);
                }
            }
        }}
        let waves=self.waves(k,cancelled)?;
        for wave in &waves {
            poll(cancelled)?;
            let scale=k*wave.weight/(8.0*PI*PI); // antipodal pair
            for i in 0..n {
                poll(cancelled)?;
                for j in i..n {
                    h[i*n+j].im+=scale*(wave.coefficients[i].conj()*wave.coefficients[j]).re;
                    h[j*n+i]=h[i*n+j];
                }
            }
        }
        for row in h.chunks_exact(n) {
            poll(cancelled)?;
            if row.iter().any(|v|!finite(*v)) {return Err(bad("nonfinite reciprocal frequency operator"));}
        }
        Ok((h,waves))
    }

    /// Solve a whole source batch with one complex symmetric operator/factor.
    /// All sources are validated before work; no partial batch is published.
    /// Cancellation reaches assembly, each dense LU elimination/substitution
    /// row, scaling/residual rows and outgoing-wave sums, not just frequencies.
    /// The immutable geometry remains reusable after any cancelled batch.
    /// # Errors
    /// Input/band/budget errors, cancellation, singular LU, nonfinite output,
    /// or failed backward-residual/physical outgoing-flux balance checks.
    pub fn solve_batch<'p>(&'p self,k:f64,medium:Medium,fields:&[&[C64]],
        mut cancelled:impl FnMut()->bool)->Result<Vec<Solution<'p,'s>>,HelmholtzError> {
        poll(&mut cancelled)?;self.admit_work(fields.len())?;
        if !k.is_finite() || k<=0.0 || k>self.highest_k
            || !medium.density.is_finite() || medium.density<=0.0
            || !medium.sound_speed.is_finite() || medium.sound_speed<=0.0
            || fields.iter().any(|v|v.len()!=self.facets.len() || v.iter().any(|x|!finite(*x))) {
            return Err(bad("invalid reciprocal medium, source field or prepared band"));
        }
        let w=k*medium.sound_speed;let omega_rho=w*medium.density;
        if !omega_rho.is_finite() || omega_rho<=0.0 {return Err(bad("reciprocal Neumann scaling is unrepresentable"));}
        let (mut h,waves)=self.assemble(k,&mut cancelled)?;let n=self.unknowns();let nv=self.vertices;
        // Congruence scaling puts auxiliary and pressure equations in comparable
        // units, without changing source work, the physical neck mass or geometry.
        for i in 0..n {
            poll(&mut cancelled)?;
            for j in 0..n {
                let scale=(if i>=nv {1.0/self.length}else{1.0})
                    *(if j>=nv {1.0/self.length}else{1.0});
                h[i*n+j]=h[i*n+j].scale(scale);
            }
        }
        let lu=lu_complex_with_cancel(&h,n,&mut cancelled).map_err(lu_error)?;
        poll(&mut cancelled)?;
        let mut output=Vec::with_capacity(fields.len());
        for &velocity in fields {
            poll(&mut cancelled)?;
            let mut rhs=zeroed(n,C64::ZERO)?;
            for (f,&v) in self.facets.iter().zip(velocity) {
                for a in f.ids {rhs[a]=rhs[a]+v*C64::new(0.0,omega_rho*f.area/3.0);}
            }
            let mut x=rhs.clone();
            lu.solve_with_cancel(&mut x,&mut cancelled).map_err(lu_error)?;
            if x.iter().any(|v|!finite(*v)) {return Err(bad("nonfinite reciprocal LU solution"));}
            let mut residual=0.0_f64;let mut denominator=0.0_f64;
            for i in 0..n {
                poll(&mut cancelled)?;
                let mut value=C64::ZERO;let mut scale=rhs[i].abs();
                for j in 0..n {value=value+h[i*n+j]*x[j];scale+=h[i*n+j].abs()*x[j].abs();}
                residual=residual.max((value-rhs[i]).abs());denominator=denominator.max(scale);
            }
            let backward_residual=residual/denominator.max(f64::MIN_POSITIVE);
            if !backward_residual.is_finite() || backward_residual>1e-10 {return Err(bad("reciprocal LU backward residual"));}
            for value in &mut x[nv..] {*value=value.scale(1.0/self.length);}
            let pressure:Vec<_>=self.facets.iter().map(|f|
                (x[f.ids[0]]+x[f.ids[1]]+x[f.ids[2]]).scale(1.0/3.0)).collect();
            let mut surface_power=C64::ZERO;let mut apparent=0.0;
            for ((p,v),f) in pressure.iter().zip(velocity).zip(&self.facets) {
                surface_power=surface_power+*p*v.conj().scale(0.5*f.area);
                apparent+=0.5*f.area*p.abs()*v.abs();
            }
            let mut outgoing=0.0;
            for wave in &waves {
                poll(&mut cancelled)?;
                let mut forward=C64::ZERO;let mut opposite=C64::ZERO;
                for (a,c) in wave.coefficients.iter().zip(&x) {
                    forward=forward+*a * *c;opposite=opposite+a.conj() * *c;
                }
                outgoing+=wave.weight*(forward.abs().powi(2)+opposite.abs().powi(2))
                    /(32.0*PI*PI*medium.density*medium.sound_speed);
            }
            let flux_balance=(surface_power.re-outgoing).abs()/apparent.max(outgoing).max(f64::MIN_POSITIVE);
            if !finite(surface_power) || !outgoing.is_finite() || !flux_balance.is_finite() || flux_balance>1e-8 {
                return Err(bad("reciprocal surface work and outgoing flux disagree"));
            }
            output.push(Solution {owner:self,k,medium,pressure,velocity:velocity.to_vec(),coefficients:x,
                surface_power,outgoing_power_w:outgoing,flux_balance,backward_residual,
                panels_per_wavelength:self.ppw_at_highest*self.highest_k/k});
        }
        poll(&mut cancelled)?;Ok(output)
    }
}

// Real Helmholtz remainders after analytic Laplace subtraction. The limits at
// a coincident quadrature node are finite; coplanar double layers are zero.
fn regular_real(k:f64,r:f64)->(f64,f64) {
    if r==0.0 {return (0.0,0.0);}
    let z=k*r;let sine=det::sin(0.5*z);
    let numerator=if z.abs()<1e-3 {
        let z2=z*z;z2*(0.5+z2*(-0.125+z2/144.0))
    }else{det::cos(z)+z*det::sin(z)-1.0};
    (-2.0*sine*sine/(4.0*PI*r),numerator/(4.0*PI*r*r*r))
}
struct Wave {weight:f64,coefficients:Vec<C64>}
impl Prepared<'_> {
    fn waves(&self,k:f64,cancelled:&mut impl FnMut()->bool)->Result<Vec<Wave>,HelmholtzError> {
        let mut waves=Vec::new();
        for (direction,weight) in directions(self.budget.angular_refinement) {
            poll(cancelled)?;
            let coefficients=self.wave_coefficients(k,direction)?;
            waves.push(Wave {weight,coefficients});
        }
        Ok(waves)
    }
    fn wave_coefficients(&self,k:f64,direction:Point)->Result<Vec<C64>,HelmholtzError> {
        let mut out=zeroed(self.unknowns(),C64::ZERO)?;
        for (i,f) in self.facets.iter().enumerate() {
            let derivative=C64::new(0.0,-k*dot(f.normal,direction));
            for sample in f.samples {
                let phase=-k*dot(direction,sample.point);
                let value=C64::new(det::cos(phase),det::sin(phase)).scale(sample.weight);
                out[self.vertices+i]=out[self.vertices+i]+value;
                for a in 0..3 {out[f.ids[a]]=out[f.ids[a]]+(value*derivative).scale(sample.shape[a]);}
            }
        }
        if out.iter().any(|v|!finite(*v)) {return Err(bad("reciprocal plane-wave integration overflow"));}
        Ok(out)
    }
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod cancellation_tests;
