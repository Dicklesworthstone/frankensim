//! Cold planar flexure of a tapered round beam. Extends the Hermite-beam owner.
//! Consistent Rayleigh inertia retains both transverse motion and section rotation.
//! The pivot fixes translation, not angle; its exact rigid rotation is separated
//! before the existing generalized eigensolve. No authored elastic frequencies.
use super::{PlateError, bad};

/// Circular section at an axial station; radius is linear between stations.
#[derive(Debug, Clone, Copy)]
pub struct RoundStation { pub x_m: f64, pub radius_m: f64 }
/// Physical and cold-work inputs. Coordinates are measured along the same shaft.
#[derive(Debug, Clone, Copy)]
pub struct RoundBeamSpec {
    pub young_pa: f64, pub density_kg_m3: f64,
    pub pivot_m: f64, pub contact_m: f64, pub hand_m: f64,
    pub subdivisions: usize, pub maximum_hz: f64, pub maximum_modes: usize,
}
/// A body rigidly attached at a shaft station, with its centre of mass ON that
/// station. Rotary inertia is about its centre, normal to the bending plane.
/// This is actual added mass, never a tip-effective mass already including the
/// shaft. Off-axis centres require a different, coupled inertial operator.
#[derive(Debug, Clone, Copy)]
pub struct RoundInertia {
    pub x_m: f64,
    pub mass_kg: f64,
    pub rotary_kg_m2: f64,
}

/// Mass-normalized eigenbasis of a pin-supported shaft, with rigid mode first.
#[derive(Debug, Clone)]
pub struct RoundBeamModes {
    pub omega: Vec<f64>,
    /// Point displacement / modal displacement, and its conjugate force map.
    pub tip: Vec<f64>, pub hand: Vec<f64>,
    pub pivot_inertia_kg_m2: f64,
    pub nodes_m: Vec<f64>,
    /// Per mode [w, theta] at each source node. theta is physical slope.
    pub shapes: Vec<Vec<f64>>,
}

/// Cubic Hermite displacement, first derivative and second derivative in x.
/// Node ordering is [w0, theta0, w1, theta1], theta=dw/dx.
pub fn hermite(t: f64, l: f64) -> ([f64;4], [f64;4], [f64;4]) {
    ([1.-3.*t*t+2.*t*t*t, l*(t-2.*t*t+t*t*t), 3.*t*t-2.*t*t*t, l*(-t*t+t*t*t)],
     [(-6.*t+6.*t*t)/l, 1.-4.*t+3.*t*t, (6.*t-6.*t*t)/l, -2.*t+3.*t*t],
     [(-6.+12.*t)/(l*l), (-4.+6.*t)/l, (6.-12.*t)/(l*l), (-2.+6.*t)/l])
}

/// Exact polynomial integration for a linearly tapered radius, constant E/rho.
/// K = integral EI N'' N''^T; M = integral rho(A N N^T + I N' N'^T).
/// Five-point Gauss integrates every term (degree at most eight), up to rounding.
pub fn round_element(l: f64, radii: [f64;2], e: f64, rho: f64)
    -> Result<([f64;16], [f64;16]), PlateError> {
    if [l,e,rho].iter().any(|v| !v.is_finite() || *v<=0.)
        || radii.iter().any(|v| !v.is_finite() || *v<0.) || radii==[0.,0.] {
        return Err(bad("round beam element needs positive length/E/rho and a nonempty radius profile"));
    }
    let mut k=[0.;16]; let mut m=[0.;16];
    for (t,w) in [(-0.906179845938664,0.2369268850561891),
        (-0.5384693101056831,0.4786286704993665),(0.,0.5688888888888889),
        (0.5384693101056831,0.4786286704993665),(0.906179845938664,0.2369268850561891)] {
        let t=0.5*(1.+t); let w=0.5*l*w;
        let radius=(1.-t)*radii[0]+t*radii[1];
        let a=std::f64::consts::PI*radius*radius; let inertia=a*radius*radius/4.;
        let (n,d,dd)=hermite(t,l);
        for i in 0..4 { for j in 0..4 {
            k[4*i+j]+=w*e*inertia*dd[i]*dd[j];
            m[4*i+j]+=w*rho*(a*n[i]*n[j]+inertia*d[i]*d[j]);
        }}
    }
    if k.iter().chain(&m).any(|v| !v.is_finite()) { return Err(bad("round beam element overflow")); }
    Ok((k,m))
}

impl RoundBeamModes {
    /// Retain every elastic eigenpair up to the explicit upper frequency, or
    /// refuse the budget. Includes the butt behind the pivot and finite tip.
    pub fn new(stations: &[RoundStation], spec: RoundBeamSpec) -> Result<Self, PlateError> {
        Self::with_inertias(stations, spec, &[])
    }

    /// Assemble attached bodies into the ORIGINAL pencil before separating the
    /// rigid mode or truncating the spectrum: M += m N N^T + J N' N'^T.
    /// All retained frequencies, force rows and rigid inertia therefore use
    /// the same loaded mass metric. Empty input preserves `new` exactly.
    pub fn with_inertias(stations: &[RoundStation], spec: RoundBeamSpec,
        attached: &[RoundInertia]) -> Result<Self, PlateError> {
        let s=spec;
        if stations.len()<2 || stations.len()>33 || s.subdivisions==0 || s.subdivisions>16
            || !(2..=17).contains(&s.maximum_modes)
            || [s.young_pa,s.density_kg_m3,s.maximum_hz].iter().any(|x| !x.is_finite() || *x<=0.)
            || [s.pivot_m,s.contact_m,s.hand_m].iter().any(|x| !x.is_finite())
            || stations.iter().any(|v| !v.x_m.is_finite() || !v.radius_m.is_finite() || v.radius_m<0.)
            || stations.windows(2).any(|w| w[0].x_m>=w[1].x_m || w[0].radius_m==0. && w[1].radius_m==0.) {
            return Err(bad("invalid bounded round beam specification"));
        }
        let start=stations[0].x_m; let end=stations[stations.len()-1].x_m; let length=end-start;
        if !length.is_finite() || length<=0. || !(start..end).contains(&s.pivot_m)
            || s.contact_m<=s.pivot_m || s.contact_m>end || s.hand_m<=s.pivot_m || s.hand_m>end {
            return Err(bad("beam pivot and positive-lever force stations must lie on the shaft"));
        }
        if attached.len()>8 || attached.iter().any(|a|
            !a.x_m.is_finite() || !(start..=end).contains(&a.x_m)
            || !a.mass_kg.is_finite() || a.mass_kg<0.
            || !a.rotary_kg_m2.is_finite() || a.rotary_kg_m2<0.
            || (a.mass_kg==0. && a.rotary_kg_m2==0.)) {
            return Err(bad("attached beam inertia needs at most eight in-profile stations and nonnegative, nonzero physical mass/inertia"));
        }
        let count=(stations.len()-1).checked_mul(s.subdivisions).and_then(|n|n.checked_add(2))
            .ok_or_else(||bad("round beam mesh overflow"))?;
        if count>66 { return Err(bad("round beam cold mesh exceeds 66 nodes")); }
        let mut nodes=Vec::with_capacity(count);
        for w in stations.windows(2) { for i in 0..s.subdivisions {
            nodes.push(w[0].x_m+(w[1].x_m-w[0].x_m)*i as f64/s.subdivisions as f64);
        }}
        nodes.push(end);
        if !nodes.contains(&s.pivot_m) { nodes.push(s.pivot_m); }
        nodes.sort_by(f64::total_cmp);
        if nodes.windows(2).any(|w|w[0]>=w[1]) {return Err(bad("round beam mesh lost axial resolution"));}
        let pivot=nodes.iter().position(|x|*x==s.pivot_m).expect("inserted pivot");
        let fixed=2*pivot; let nd=2*nodes.len()-1;
        let map=|i:usize| if i==fixed {None} else {Some(i-usize::from(i>fixed))};
        let radius=|x:f64| {
            let i=stations.partition_point(|p|p.x_m<=x).saturating_sub(1).min(stations.len()-2);
            let t=(x-stations[i].x_m)/(stations[i+1].x_m-stations[i].x_m);
            (1.-t)*stations[i].radius_m+t*stations[i+1].radius_m
        };
        let mut k=vec![0.;nd*nd]; let mut m=vec![0.;nd*nd];
        for (a,x) in nodes.windows(2).enumerate() {
            let (ke,me)=round_element(x[1]-x[0],[radius(x[0]),radius(x[1])],s.young_pa,s.density_kg_m3)?;
            // Scale slopes by total length to avoid mixing metre and radian
            // magnitudes in the generalized eigensolver. Ports undo this scale.
            for i in 0..4 { for j in 0..4 {
                if let (Some(row),Some(col))=(map(2*a+i),map(2*a+j)) {
                    let scale=(if i%2==0 {1.} else {1./length})*(if j%2==0 {1.} else {1./length});
                    k[row*nd+col]+=scale*ke[4*i+j]; m[row*nd+col]+=scale*me[4*i+j];
                }
            }}
        }
        for a in attached {
            let cell=nodes.partition_point(|x|*x<=a.x_m).saturating_sub(1).min(nodes.len()-2);
            let l=nodes[cell+1]-nodes[cell];
            let (n,d,_)=hermite((a.x_m-nodes[cell])/l,l);
            for i in 0..4 { for j in 0..4 {
                if let (Some(row),Some(col))=(map(2*cell+i),map(2*cell+j)) {
                    let scale=(if i%2==0 {1.} else {1./length})*(if j%2==0 {1.} else {1./length});
                    m[row*nd+col]+=scale*(a.mass_kg*n[i]*n[j]+a.rotary_kg_m2*d[i]*d[j]);
                }
            }}
        }
        if m.iter().any(|v|!v.is_finite()) {return Err(bad("attached beam mass matrix overflow"));}
        let mut rigid=vec![0.;nd];
        for (i,&x) in nodes.iter().enumerate() {
            if let Some(j)=map(2*i) {rigid[j]=x-s.pivot_m;}
            rigid[map(2*i+1).unwrap()]=length;
        }
        let mr:Vec<f64>=(0..nd).map(|i|(0..nd).map(|j|m[i*nd+j]*rigid[j]).sum()).collect();
        let inertia=mr.iter().zip(&rigid).map(|(a,b)|a*b).sum::<f64>();
        if !inertia.is_finite() || inertia<=0. {return Err(bad("beam pivot inertia is not representable"));}
        // The omitted pivot-slope entry has a nonzero rigid coefficient. The
        // remaining unit columns, projected M-orthogonally off rigid rotation,
        // span its entire complement. K*r=0 analytically; no tiny spring is added.
        let omit=map(fixed+1).unwrap(); let ids:Vec<_>=(0..nd).filter(|i|*i!=omit).collect(); let n=ids.len();
        let mut kr=vec![0.;n*n]; let mut mm=vec![0.;n*n];
        for (i,&a) in ids.iter().enumerate() { for (j,&b) in ids.iter().enumerate() {
            kr[i*n+j]=k[a*nd+b]; mm[i*n+j]=m[a*nd+b]-mr[a]*mr[b]/inertia;
        }}
        let modes=fs_modal::eigh_gen_dense(&kr,&mm,n).map_err(|_|bad("round beam eigensolve refused"))?;
        let root=inertia.sqrt();
        let mut shapes=vec![vec![0.;2*nodes.len()]];
        for (i,&x) in nodes.iter().enumerate() { shapes[0][2*i]=(x-s.pivot_m)/root; shapes[0][2*i+1]=1./root; }
        let mut omega=vec![0.];
        for mode in modes {
            if !mode.lambda.is_finite() || mode.lambda<=0. || !mode.residual.is_finite()
                || mode.phi.len()!=n || mode.phi.iter().any(|v|!v.is_finite()) {
                return Err(bad("round beam elastic spectrum is not finite positive"));
            }
            let w=mode.lambda.sqrt(); if w>std::f64::consts::TAU*s.maximum_hz {continue;}
            if omega.len()==s.maximum_modes {return Err(bad("round beam frequency slice exceeds mode budget; no truncation"));}
            let c=ids.iter().zip(&mode.phi).map(|(&i,p)|mr[i]*p).sum::<f64>()/inertia;
            let mut full:Vec<_>=rigid.iter().map(|r|-r*c).collect();
            for (&i,&p) in ids.iter().zip(&mode.phi) {full[i]+=p;}
            let mut physical=vec![0.;2*nodes.len()];
            for (i,out) in physical.iter_mut().enumerate() {if let Some(j)=map(i){*out=full[j]/if i%2==0{1.}else{length};}}
            let mut norm=0.; let mut error=0.;
            for i in 0..nd {
                let mv=(0..nd).map(|j|m[i*nd+j]*full[j]).sum::<f64>();
                let kv=(0..nd).map(|j|k[i*nd+j]*full[j]).sum::<f64>();
                norm+=full[i]*mv; error+= (kv-mode.lambda*mv).powi(2);
            }
            let scale=mode.lambda*full.iter().map(|v|v*v).sum::<f64>().sqrt()
                *m.iter().map(|v|v.abs()).fold(0.,f64::max)*nd as f64;
            if !norm.is_finite() || (norm-1.).abs()>1e-7 || error.sqrt()>1e-6*scale {
                return Err(bad("round beam eigenpair failed original-pencil validation"));
            }
            shapes.push(physical); omega.push(w);
        }
        if omega.len()==1 {return Err(bad("flexible shaft must retain at least one elastic mode"));}
        let mut result=Self {omega, tip:Vec::new(), hand:Vec::new(),pivot_inertia_kg_m2:inertia,nodes_m:nodes,shapes};
        result.tip=result.point(s.contact_m)?; result.hand=result.point(s.hand_m)?;
        // The exact positive rigid participation must not inherit an arbitrary
        // eigenvector sign. Initial launch and SI hand force use this coordinate.
        result.tip[0]=(s.contact_m-s.pivot_m)/root; result.hand[0]=(s.hand_m-s.pivot_m)/root;
        Ok(result)
    }
    /// Conjugate displacement/force row at a material station of the same beam.
    pub fn point(&self, x: f64) -> Result<Vec<f64>, PlateError> {
        self.project(x, false)
    }
    /// Section rotation / modal displacement. Its transpose applies a physical
    /// bending moment, or the lever-arm moment of a distributed face force.
    pub fn slope(&self, x: f64) -> Result<Vec<f64>, PlateError> {
        self.project(x, true)
    }
    fn project(&self, x: f64, rotation: bool) -> Result<Vec<f64>, PlateError> {
        if !x.is_finite() || x<self.nodes_m[0] || x>self.nodes_m[self.nodes_m.len()-1] {
            return Err(bad("beam point lies outside supplied geometry"));
        }
        let i=self.nodes_m.partition_point(|p|*p<=x).saturating_sub(1).min(self.nodes_m.len()-2);
        let l=self.nodes_m[i+1]-self.nodes_m[i]; let (n,d,_)=hermite((x-self.nodes_m[i])/l,l);
        let n=if rotation {d}else{n};
        let out:Vec<f64>=self.shapes.iter().map(|p|(0..4).map(|j|n[j]*p[2*i+j]).sum()).collect();
        if out.iter().any(|x|!x.is_finite()) {return Err(bad("beam point projection overflow"));} Ok(out)
    }
}

#[cfg(test)]
#[path="beam_tests.rs"]
mod tests;
