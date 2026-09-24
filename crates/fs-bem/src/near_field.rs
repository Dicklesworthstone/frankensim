//! Finite receivers near a retained closed triangle surface.
//!
//! The boundary solve remains `helmholtz`. This observer integrates its SAME
//! constant pressure/normal-velocity traces with adaptive triangle quadrature,
//! rather than treating a nearby panel as a monopole at its centroid. Admission
//! uses exact indexed closure after exact-coordinate welding, component-wise
//! solid-angle winding, and point-to-triangle clearance. No enclosing-sphere
//! exclusion, surface repair, tolerance welding or new boundary unknowns.
//!
//! Inputs must be disjoint, outward, non-self-intersecting closed components.
//! Closure/winding checks are not a global intersection certificate. Quadrature
//! discrepancies below are estimates, not bounds on the BEM or model error.
use crate::{helmholtz::{self, HelmholtzError, Medium, RadiationSolution}, panel3d::SpherePanels};
use fs_math::{c64::C64, det};
use std::collections::{BTreeMap, BTreeSet};

#[path="near_field_velocity.rs"]
mod velocity;
pub use velocity::{FirstOrder, VectorObservation};

const MAX_PANELS: usize = helmholtz::MAX_DENSE_PANELS;
const MAX_RECEIVERS: usize = 64;
const FOUR_PI: f64 = 4. * std::f64::consts::PI;
type Point = [f64; 3];
type Triangle = [Point; 3];
fn bad(what: &'static str) -> HelmholtzError { HelmholtzError::BadParameter { what } }
fn sub(a: Point, b: Point) -> Point { std::array::from_fn(|i| a[i]-b[i]) }
fn dot(a: Point, b: Point) -> f64 { a.iter().zip(b).map(|(a,b)| a*b).sum() }
fn cross(a: Point, b: Point) -> Point { [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]] }
fn norm(a: Point) -> f64 { dot(a,a).sqrt() }
fn midpoint(a: Point, b: Point) -> Point { std::array::from_fn(|i| a[i]+0.5*(b[i]-a[i])) }

// The input coordinate envelope prevents products here from overflowing. Work
// near a far-away tiny body is additionally refused at its roundoff clearance.
fn distance(x: Point, t: Triangle) -> f64 {
    let [a,b,c]=t; let u=sub(b,a); let v=sub(c,a); let w=sub(x,a);
    let n=cross(u,v); let nn=dot(n,n);
    let beta=dot(cross(w,v),n)/nn; let gamma=dot(cross(u,w),n)/nn;
    if beta>=0. && gamma>=0. && beta+gamma<=1. {
        return dot(w,n).abs()/nn.sqrt();
    }
    [(a,b),(b,c),(c,a)].into_iter().map(|(a,b)| {
        let e=sub(b,a); let d=sub(x,a); let s=(dot(d,e)/dot(e,e)).clamp(0.,1.);
        norm(std::array::from_fn(|i| d[i]-s*e[i]))
    }).fold(f64::INFINITY,f64::min)
}
// Signed triangle solid angle, normalized before evaluation to avoid scale
// products. Van Oosterom-Strackee / generalized winding-number formulation.
fn angle(x: Point, t: Triangle) -> f64 {
    let [a,b,c]=t.map(|p| { let v=sub(p,x); let r=norm(v); v.map(|q|q/r) });
    2.*det::atan2(dot(a,cross(b,c)),1.+dot(a,b)+dot(b,c)+dot(c,a))
}
fn compensated(sum: &mut f64, correction: &mut f64, value: f64) {
    let y=value-*correction; let next=*sum+y; *correction=(next-*sum)-y; *sum=next;
}

/// Cold work and integration accuracy controls. A failed request returns no
/// prepared observer; it never substitutes centroid evaluation.
#[derive(Clone, Copy, Debug)]
pub struct Options {
    /// Relative local discrepancy against the integral of absolute kernel size.
    pub relative_tolerance: f64,
    /// Triangle subdivision depth, at most 16.
    pub maximum_depth: usize,
    /// Kernel evaluations across ALL receivers/panels at this frequency.
    pub maximum_kernel_evaluations: usize,
}
impl Default for Options {
    fn default() -> Self { Self { relative_tolerance:1e-7, maximum_depth:14, maximum_kernel_evaluations:2_000_000 } }
}
impl Options {
    fn validate(self) -> Result<(),HelmholtzError> {
        if !self.relative_tolerance.is_finite() || !(1e-12..=1e-3).contains(&self.relative_tolerance)
            || self.maximum_depth>16 || self.maximum_kernel_evaluations<80
            || self.maximum_kernel_evaluations>20_000_000 {
            return Err(bad("invalid near-field quadrature accuracy or work budget"));
        }
        Ok(())
    }
}

/// Admitted receiver positions and actual distance to the nearest panel.
/// Retains a borrow of the immutable source surface. Reuse for a frequency grid.
pub struct Geometry<'a> {
    surface: &'a SpherePanels,
    points: Vec<Point>,
    clearances: Vec<f64>,
}
impl<'a> Geometry<'a> {
    /// Require exact triangle geometry, a positive clearance, and outside every
    /// consistently outward closed component. Distances use faces AND edges,
    /// not centroids. This does not certify global mesh nonintersection.
    /// # Errors
    /// Nonfinite/unresolved coordinates, open/inward surfaces, an interior or
    /// near-surface point, missing triangle geometry, or cold resource bounds.
    pub fn new(surface: &'a SpherePanels, points: &[Point], minimum_clearance_m: f64) -> Result<Self,HelmholtzError> {
        let triangles=surface.triangles().ok_or_else(||bad("near-field receivers require retained triangles"))?;
        if triangles.is_empty() || triangles.len()>MAX_PANELS || points.is_empty() || points.len()>MAX_RECEIVERS
            || !minimum_clearance_m.is_finite() || minimum_clearance_m<=0.
            || triangles.iter().flatten().flatten().chain(points.iter().flatten()).any(|v|!v.is_finite()||v.abs()>1e9) {
            return Err(bad("invalid bounded near-field surface, receiver or clearance"));
        }
        let components=components(triangles)?;
        let scale=triangles.iter().flatten().flatten().chain(points.iter().flatten())
            .fold(0.0_f64,|s,v|s.max(v.abs()));
        let guard=minimum_clearance_m.max(128.*f64::EPSILON*scale);
        let mut clearances=Vec::with_capacity(points.len());
        for &point in points {
            let clearance=triangles.iter().map(|&t|distance(point,t)).fold(f64::INFINITY,f64::min);
            if !clearance.is_finite() || clearance<guard { return Err(bad("receiver touches or is too close to a source triangle")); }
            for component in &components {
                let(mut sum,mut correction)=(0.,0.);
                for &i in component { compensated(&mut sum,&mut correction,angle(point,triangles[i])); }
                // Refuse ambiguous values instead of guessing from a half-space
                // or approximate centroid solid-angle sign. No cancellation
                // between different components can hide an interior receiver.
                if !sum.is_finite() || (sum/FOUR_PI).abs()>1e-7 {
                    return Err(bad("receiver is interior or its exterior winding is unresolved"));
                }
            }
            clearances.push(clearance);
        }
        Ok(Self {surface,points:points.to_vec(),clearances})
    }
    /// Minimum Euclidean point-to-triangle distance [m], in receiver order.
    pub fn clearances_m(&self) -> &[f64] { &self.clearances }
    /// Original, unshifted physical receiver coordinates [m].
    pub fn points_m(&self) -> &[Point] { &self.points }
    /// Integrate a shared pair of Green-operator rows per receiver. Apply to
    /// all modal solutions at this wavenumber without repeating quadrature.
    /// # Errors
    /// Invalid k/medium, exhausted refinement/work, or nonfinite quadrature.
    pub fn prepare(&self, k:f64, medium:Medium, options:Options) -> Result<Prepared<'a>,HelmholtzError> {
        self.prepare_fields(k,medium,options,false)
    }
    fn prepare_fields(&self,k:f64,medium:Medium,options:Options,with_velocity:bool) -> Result<Prepared<'a>,HelmholtzError> {
        options.validate()?;
        let omega_rho=k*medium.sound_speed*medium.density;
        if !k.is_finite() || k<=0. || !medium.density.is_finite() || medium.density<=0.
            || !medium.sound_speed.is_finite() || medium.sound_speed<=0. || !omega_rho.is_finite() || omega_rho<=0. {
            return Err(bad("near-field evaluation needs positive finite k and medium"));
        }
        let triangles=self.surface.triangles().ok_or_else(||bad("missing retained receiver triangles"))?;
        let mut work=Work {options,evaluations:0};
        let mut rows=Vec::with_capacity(self.points.len());
        for &x in &self.points {
            let mut row=Vec::with_capacity(triangles.len());
            for (i,&t) in triangles.iter().enumerate() {
                row.push(integrate_fields(k,x,t,self.surface.normals()[i],0,&mut work,with_velocity)?);
            }
            rows.push(row);
        }
        Ok(Prepared {surface:self.surface,k,medium,rows,with_velocity,kernel_evaluations:work.evaluations})
    }
}

// Exact-coordinate seam welding; no tolerance repair. Component-wise winding
// needs oriented closed components, not the approximate centroid flux test.
pub(crate) fn components(triangles:&[Triangle]) -> Result<Vec<Vec<usize>>,HelmholtzError> {
    let mut vertices=BTreeMap::<[u64;3],usize>::new();
    let mut edges=BTreeMap::<(usize,usize),Vec<(usize,usize,usize)>>::new(); let mut faces=BTreeSet::new();
    for (face,t) in triangles.iter().enumerate() {
        let mut ids=[0;3];
        for (i,p) in t.iter().enumerate() {
            let key=p.map(|v|if v==0. {0}else{v.to_bits()}); let next=vertices.len();
            ids[i]=*vertices.entry(key).or_insert(next);
        }
        let mut sorted=ids; sorted.sort_unstable();
        let area=norm(cross(sub(t[1],t[0]),sub(t[2],t[0])));
        if sorted[0]==sorted[1] || sorted[1]==sorted[2] || !faces.insert(sorted) || !area.is_finite() || area<=0. {
            return Err(bad("duplicate or unresolved near-field triangle"));
        }
        for i in 0..3 {let(a,b)=(ids[i],ids[(i+1)%3]); edges.entry((a.min(b),a.max(b))).or_default().push((face,a,b));}
    }
    let mut adjacency=vec![Vec::new();triangles.len()];
    for uses in edges.values() {
        if uses.len()!=2 || uses[0].1!=uses[1].2 || uses[0].2!=uses[1].1 {
            return Err(bad("near-field surface must have two opposite uses of each edge"));
        }
        adjacency[uses[0].0].push(uses[1].0); adjacency[uses[1].0].push(uses[0].0);
    }
    let mut visited=vec![false;triangles.len()]; let mut result=Vec::new();
    for start in 0..triangles.len() {
        if visited[start] {continue;}
        let origin=triangles[start][0]; let(mut volume,mut correction)=(0.,0.);
        let mut pending=vec![start]; let mut component=Vec::new(); visited[start]=true;
        while let Some(i)=pending.pop() {
            component.push(i); let [a,b,c]=triangles[i].map(|p|sub(p,origin));
            compensated(&mut volume,&mut correction,dot(a,cross(b,c)));
            for &j in &adjacency[i] {if !visited[j] {visited[j]=true;pending.push(j);}}
        }
        if !volume.is_finite() || volume<=0. {return Err(bad("near-field components must enclose positive outward volume"));}
        result.push(component);
    }
    Ok(result)
}

#[derive(Clone,Copy,Default)]
struct Integral { s:C64, d:C64, abs_s:f64, abs_d:f64, error_s:f64, error_d:f64, gradient:velocity::Derivative }
impl Integral {
    fn add(&mut self,b:Self) {
        self.s=self.s+b.s;self.d=self.d+b.d;self.abs_s+=b.abs_s;self.abs_d+=b.abs_d;
        self.error_s+=b.error_s;self.error_d+=b.error_d;
        self.gradient.add(b.gradient);
    }
}
struct Work { options:Options, evaluations:usize }
const X4:[f64;4]=[0.06943184420297371,0.33000947820757187,0.6699905217924281,0.9305681557970262];
const W4:[f64;4]=[0.17392742256872693,0.32607257743127307,0.32607257743127307,0.17392742256872693];
const X8:[f64;8]=[0.019855071751231884,0.10166676129318664,0.2372337950418355,0.4082826787521751,
    0.5917173212478249,0.7627662049581645,0.8983332387068134,0.9801449282487681];
const W8:[f64;8]=[0.05061426814518813,0.11119051722668724,0.15685332293894366,0.181341891689181,
    0.181341891689181,0.15685332293894366,0.11119051722668724,0.05061426814518813];
fn rule(k:f64,x:Point,t:Triangle,normal:Point,nodes:&[f64],weights:&[f64],work:&mut Work,with_velocity:bool) -> Result<Integral,HelmholtzError> {
    let count=nodes.len()*nodes.len();
    work.evaluations=work.evaluations.checked_add(count).filter(|n|*n<=work.options.maximum_kernel_evaluations)
        .ok_or_else(||bad("near-field kernel-evaluation budget exhausted"))?;
    let u=sub(t[1],t[0]);let v=sub(t[2],t[1]);let area2=norm(cross(u,sub(t[2],t[0])));
    let mut out=Integral::default();
    for (&a,&wa) in nodes.iter().zip(weights) {for (&b,&wb) in nodes.iter().zip(weights) {
        let y=std::array::from_fn(|c|t[0][c]+a*u[c]+a*b*v[c]);let delta=sub(x,y);let r=norm(delta);
        if !r.is_finite() || r<=0. {return Err(bad("unresolved near-field integration point"));}
        let phase=k*r;let g=C64::new(det::cos(phase),det::sin(phase)).scale(1./(FOUR_PI*r));
        let projection=dot(normal,delta)/r;let dg=g*C64::new(1./r,-k).scale(projection);
        let jac=wa*wb*area2*a;
        out.s=out.s+g.scale(jac);out.d=out.d+dg.scale(jac);
        if with_velocity {out.gradient.accumulate(k,g,delta,r,normal,jac)?;}
        out.abs_s+=jac/(FOUR_PI*r);out.abs_d+=jac/(FOUR_PI*r)*(1./r).hypot(k)*projection.abs();
    }}
    if ![out.s.re,out.s.im,out.d.re,out.d.im,out.abs_s,out.abs_d].iter().all(|v|v.is_finite()) {
        return Err(bad("nonfinite near-field quadrature"));
    }
    Ok(out)
}
fn children([a,b,c]:Triangle) -> [Triangle;4] {
    let ab=midpoint(a,b);let bc=midpoint(b,c);let ca=midpoint(c,a);
    [[a,ab,ca],[ab,b,bc],[ca,bc,c],[ab,bc,ca]]
}
fn integrate(k:f64,x:Point,t:Triangle,normal:Point,depth:usize,work:&mut Work) -> Result<Integral,HelmholtzError> {
    integrate_fields(k,x,t,normal,depth,work,false)
}
fn integrate_fields(k:f64,x:Point,t:Triangle,normal:Point,depth:usize,work:&mut Work,with_velocity:bool) -> Result<Integral,HelmholtzError> {
    let coarse=rule(k,x,t,normal,&X4,&W4,work,with_velocity)?;let mut fine=rule(k,x,t,normal,&X8,&W8,work,with_velocity)?;
    fine.error_s=(fine.s-coarse.s).abs();fine.error_d=(fine.d-coarse.d).abs();
    let size=[norm(sub(t[1],t[0])),norm(sub(t[2],t[1])),norm(sub(t[0],t[2]))].into_iter().fold(0.0_f64,f64::max);
    let separated=size<=2.*distance(x,t) && k*size<=3.;
    let tol=work.options.relative_tolerance;
    // The geometric/phase condition prevents two rules from jointly missing
    // a sharp nearby peak or an under-sampled oscillation.
    let gradient_ok=!with_velocity || fine.gradient.accept(&coarse.gradient,tol);
    if separated && fine.error_s<=tol*fine.abs_s && fine.error_d<=tol*fine.abs_d && gradient_ok {return Ok(fine);}
    if depth==work.options.maximum_depth {return Err(bad("near-field triangle refinement depth exhausted"));}
    let mut sum=Integral::default();
    for child in children(t) {sum.add(integrate_fields(k,x,child,normal,depth+1,work,with_velocity)?);}
    Ok(sum)
}

/// Frequency-specific Green representation rows. No new solve is performed.
pub struct Prepared<'a> {
    surface:&'a SpherePanels,k:f64,medium:Medium,rows:Vec<Vec<Integral>>,
    with_velocity:bool,
    /// Actual quadrature work during preparation, shared by every modal field.
    pub kernel_evaluations:usize,
}
/// Evaluated pressure and propagated local quadrature discrepancy estimates.
pub struct Observation {
    /// Complex pressure [Pa], receiver order, exp(-i omega t).
    pub pressure:Vec<C64>,
    /// Sum of |pressure_j| error(D_j) + |i omega rho velocity_j| error(S_j).
    /// Does not cover solve, spatial discretization or model error.
    pub quadrature_error_estimate_pa:Vec<f64>,
}
impl Prepared<'_> {
    /// Apply the shared rows to an existing solution at the identical k, medium
    /// and surface. Public field data are rechecked; cross-wiring is refused.
    /// # Errors
    /// A mismatched or nonfinite solution, or nonfinite evaluated pressure.
    pub fn evaluate(&self,solution:&RadiationSolution) -> Result<Observation,HelmholtzError> {
        // Reuse the original solver's fingerprint, shape and medium admission.
        // Empty points perform no old centroid receiver evaluation.
        helmholtz::exterior_pressure_at_points(self.surface,solution,self.medium,&[])?;
        if solution.k.to_bits()!=self.k.to_bits() {return Err(bad("receiver rows and solution wavenumbers differ"));}
        let omega_rho=self.k*self.medium.sound_speed*self.medium.density;
        let mut pressure=Vec::with_capacity(self.rows.len());let mut errors=Vec::with_capacity(self.rows.len());
        for row in &self.rows {
            let(mut re,mut im,mut cr,mut ci,mut error)=(0.,0.,0.,0.,0.);
            for ((cell,&p),&v) in row.iter().zip(&solution.pressure).zip(&solution.velocity) {
                let q=v*C64::new(0.,omega_rho);let term=cell.d*p-cell.s*q;
                compensated(&mut re,&mut cr,term.re);compensated(&mut im,&mut ci,term.im);
                error+=p.abs()*cell.error_d+q.abs()*cell.error_s;
            }
            if ![re,im,error].iter().all(|v|v.is_finite()) {return Err(bad("nonfinite near-field pressure or error estimate"));}
            pressure.push(C64::new(re,im));errors.push(error);
        }
        Ok(Observation {pressure,quadrature_error_estimate_pa:errors})
    }
}

#[cfg(test)]
#[path="near_field_tests.rs"]
mod tests;

#[cfg(test)]
#[path="near_field_velocity_tests.rs"]
mod velocity_tests;
