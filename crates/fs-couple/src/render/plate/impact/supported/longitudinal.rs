//! Axial waves and geometric endpoint extension in the existing support owner.
//! All strain polynomials and their derivatives come from fs-nlmodal. The
//! support supplies one rigid translation; its axial projection is geometry,
//! never an authored gain. The reference is already tensioned and equilibrated:
//! static prestress reactions are not applied a second time here.
use super::*;
use fs_nlmodal::{KcStringParams, longitudinal::LongitudinalString};
mod viscosity;
pub mod relaxation;

/// A physical span and the direction of the common translating termination.
/// The first transverse axis fixes the source's existing two-plane frame; the
/// second is string_axis × first_transverse_axis. Inputs are kept verbatim.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StringSupportGeometry {
    endpoints_m: [[f64; 3]; 2],
    first_transverse_axis: [f64; 3],
    support_direction: [f64; 3],
    components: [f64; 3],
    length_m: f64,
}
fn dot(a: &[f64], b: &[f64]) -> f64 { a.iter().zip(b).map(|(a,b)|a*b).sum() }
impl StringSupportGeometry {
    /// Resolve the two transverse and axial projections from actual geometry.
    pub fn try_new(endpoints_m: [[f64; 3]; 2], first_transverse_axis: [f64; 3],
        support_direction: [f64; 3]) -> Result<Self, ImpactError> {
        if endpoints_m.iter().flatten().chain(&first_transverse_axis).chain(&support_direction)
            .any(|v| !v.is_finite()) {
            return Err(invalid("support geometry requires finite positions and directions"));
        }
        let delta=std::array::from_fn::<_,3,_>(|i|endpoints_m[1][i]-endpoints_m[0][i]);
        let length_m=det::sqrt(dot(&delta,&delta));
        if !length_m.is_finite() || length_m<=0.0
            || (dot(&first_transverse_axis,&first_transverse_axis)-1.0).abs()>1e-12
            || (dot(&support_direction,&support_direction)-1.0).abs()>1e-12 {
            return Err(invalid("support geometry needs a nonzero span and unit physical axes"));
        }
        let axis=delta.map(|v|v/length_m);
        if dot(&axis,&first_transverse_axis).abs()>1e-12 {
            return Err(invalid("string transverse frame must be perpendicular to its actual span"));
        }
        let second=[axis[1]*first_transverse_axis[2]-axis[2]*first_transverse_axis[1],
            axis[2]*first_transverse_axis[0]-axis[0]*first_transverse_axis[2],
            axis[0]*first_transverse_axis[1]-axis[1]*first_transverse_axis[0]];
        let components=[dot(&support_direction,&first_transverse_axis),
            dot(&support_direction,&second),dot(&support_direction,&axis)];
        Ok(Self {endpoints_m,first_transverse_axis,support_direction,components,length_m})
    }
    pub const fn endpoints_m(&self) -> [[f64; 3]; 2] { self.endpoints_m }
    pub const fn first_transverse_axis(&self) -> [f64; 3] { self.first_transverse_axis }
    pub const fn support_direction(&self) -> [f64; 3] { self.support_direction }
    /// Components in source order: first transverse, second transverse, axial.
    pub const fn components(&self) -> [f64; 3] { self.components }
    pub const fn length_m(&self) -> f64 { self.length_m }
}

/// Explicit axial resolution and solid inertia. EA remains the original
/// StringStretching input; bending/gas damping is never reused as axial loss.
#[derive(Clone, Debug, PartialEq)]
pub struct LongitudinalSupportModes {
    pub solid_linear_density_kg_m: f64,
    pub modes: Vec<usize>,
    pub maximum_strain: f64,
}
impl LongitudinalSupportModes {
    pub fn validate(&self) -> Result<(), ImpactError> {
        if !self.solid_linear_density_kg_m.is_finite() || self.solid_linear_density_kg_m<=0.0
            || self.modes.is_empty() || self.modes.len()>fs_nlmodal::longitudinal::MAX_LONGITUDINAL_MODES
            || self.modes[0]==0 || self.modes.last().is_some_and(|&n|n>fs_nlmodal::longitudinal::MAX_LONGITUDINAL_INDEX)
            || self.modes.windows(2).any(|w|w[0]>=w[1])
            || !self.maximum_strain.is_finite() || !(0.0<self.maximum_strain && self.maximum_strain<=0.1) {
            return Err(invalid("axial support needs positive solid mass, explicit ordered modes and a finite strain limit"));
        }
        Ok(())
    }
}

/// Bounds on the complete retained axial strain, including endpoint extension.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LongitudinalSupportObservation {
    pub mean_dynamic_strain: f64,
    pub maximum_abs_dynamic_strain: f64,
    pub minimum_tension_n: f64,
    pub maximum_tension_n: f64,
    pub strain_energy_j: f64,
}

impl SupportedStrings {
    /// Bind explicit small-strain axial viscous rigidity C_a [N s] once, in
    /// original string order. The source owner must admit its material bands.
    ///
    /// For u=h*d_axial*X+sum(phi_n*q_n), cosine orthogonality gives
    /// integral(C_a*u_xt^2 dx) = sum(C_a/EA*omega_n^2*qdot_n^2)
    /// + C_a*L*(h_x*d_axial)^2*Xdot^2. The original kinetic pullback retains
    /// every qdot_n = zdot_n-c_n*Xdot cross term; the affine endpoint term is
    /// essential when one end moves. Equal endpoint translation has no strain
    /// and therefore contributes no extra endpoint resistance.
    ///
    /// This is the same linear axial Kelvin--Voigt approximation as the fixed
    /// endpoint bank. Geometric y_x*y_xt loss, heat and axial gas drag are not
    /// supplied. Zero is explicit lossless data; a missing string is an error.
    /// No state, mass, elastic force, transverse or torsional loss is replaced.
    pub fn with_axial_viscosity(self, viscous_rigidity_n_s: Vec<f64>)
        -> Result<Self, ImpactError>
    {
        self.bind_axial_viscosity(viscous_rigidity_n_s,false)
    }

    /// Bind Kelvin--Voigt loss to the SAME retained moderate-slope strain as
    /// storage: epsilon=u_x+(y_x^2+z_x^2)/2. Its rate includes both transverse
    /// planes, longitudinal waves and actual affine endpoint motion. The
    /// reciprocal force and analytic tangent enter the existing joint solve;
    /// their positive work enters its existing dissipation ledger.
    ///
    /// This replaces the linear axial approximation, never adds it twice.
    /// It preserves the original bending, torsional, gas and support losses.
    /// Source bands must cover the whole coupled image. No thermal feedback,
    /// slip, hereditary axial law or geometrically exact large strain is implied.
    pub fn with_finite_strain_axial_viscosity(self, viscous_rigidity_n_s: Vec<f64>)
        -> Result<Self, ImpactError>
    {
        self.bind_axial_viscosity(viscous_rigidity_n_s,true)
    }

    fn bind_axial_viscosity(mut self, viscous_rigidity_n_s: Vec<f64>, finite_strain: bool)
        -> Result<Self, ImpactError>
    {
        if self.0.axial_viscosity_n_s.is_some()
            || viscous_rigidity_n_s.len()!=self.wire_count()
            || viscous_rigidity_n_s.iter().any(|v|!v.is_finite() || *v<0.0) {
            return Err(invalid("axial viscosity needs one cold finite nonnegative rigidity per supported string"));
        }
        let mut damping=self.0.damping.clone();
        let mut endpoint_damping=0.0;
        for (index,&rigidity) in viscous_rigidity_n_s.iter().enumerate() {
            let (Some(range),SupportedPotential::Longitudinal(potential))=
                (&self.0.axial_ranges[index],&self.0.potentials[index]) else {
                    return Err(invalid("axial viscosity requires the original longitudinal support coordinates"));
                };
            let ratio=rigidity/potential.ea;
            for (index,&omega) in range.clone().zip(potential.axial_omegas()) {
                let rate=ratio*omega*omega;
                if !rate.is_finite() || (rigidity>0.0 && rate==0.0) {
                    return Err(invalid("supported axial viscous rate is unrepresentable"));
                }
                if !finite_strain {damping[index]=rate;}
            }
            let chord=potential.chord[2];
            let resistance=rigidity*self.0.wires[index].length_m()*chord*chord;
            if !resistance.is_finite() || (rigidity>0.0 && chord!=0.0 && resistance==0.0) {
                return Err(invalid("supported axial endpoint resistance is unrepresentable"));
            }
            endpoint_damping+=resistance;
        }
        let normalized=endpoint_damping*self.0.weight*self.0.weight;
        if !normalized.is_finite() || (endpoint_damping>0.0 && normalized==0.0) {
            return Err(invalid("supported axial resistance normalization is unrepresentable"));
        }
        let data=Rc::make_mut(&mut self.0);
        data.damping=damping;
        data.axial_endpoint_damping_n_s_m=if finite_strain {0.0} else {endpoint_damping};
        data.axial_viscosity_n_s=Some(viscous_rigidity_n_s);
        data.finite_strain_axial_viscosity=finite_strain;
        Ok(self)
    }

    /// Original per-string axial viscous rigidities [N s], without folding
    /// them into the independent support's ground resistance or gas drag.
    pub fn axial_viscosity_n_s(&self) -> Option<&[f64]> {
        self.0.axial_viscosity_n_s.as_deref()
    }
}

#[derive(Clone, Debug)]
pub(super) struct LongitudinalPotential {
    model: LongitudinalString,
    omegas: Vec<f64>,
    scales: Vec<f64>,
    chord: [f64; 3],
    lift_stiffness: f64,
    ea: f64,
    tension: f64,
    maximum_slope: f64,
    maximum_strain: f64,
}
struct Strain {
    slopes: Vec<f64>,
    mean: f64,
    harmonics: Vec<f64>,
    energy: f64,
    spread: f64,
}
impl LongitudinalPotential {
    pub(super) fn new(wire: &WireSpan, law: StringStretching, selection: &LongitudinalSupportModes,
        geometry: StringSupportGeometry, ends: StringEndMotion, omegas: Vec<f64>) -> Result<Self, ImpactError> {
        selection.validate()?;
        if (geometry.length_m-wire.length_m()).abs()>1e-12*wire.length_m()
            || omegas.len()!=2*wire.damping_per_s.len() || law.axial_rigidity_n<=0.0 {
            return Err(invalid("axial support must retain the original span, both planes and positive material EA"));
        }
        let length=wire.length_m();
        let sign=match ends {StringEndMotion::Both=>0.0,StringEndMotion::First=>-1.0,StringEndMotion::Second=>1.0};
        let chord=geometry.components.map(|d|sign*d/length);
        let model=LongitudinalString::new(KcStringParams {length,lin_density:selection.solid_linear_density_kg_m,
            tension:wire.tension_n,ea:law.axial_rigidity_n},wire.damping_per_s.len(),2,&selection.modes)
            .map_err(|e|ImpactError::Owner(e.to_string()))?;
        Ok(Self {model,omegas,scales:(1..=wire.damping_per_s.len()).map(|n|
            det::sqrt(2.0/(wire.linear_density_kg_m*length))*n as f64*core::f64::consts::PI/length).collect(),
            chord,lift_stiffness:wire.tension_n*length*(chord[0]*chord[0]+chord[1]*chord[1]),
            ea:law.axial_rigidity_n,tension:wire.tension_n,maximum_slope:law.maximum_slope,
            maximum_strain:selection.maximum_strain})
    }
    pub(super) fn axial_omegas(&self) -> &[f64] { self.model.axial_omegas() }
    pub(super) fn axial_shape(&self,station:f64,length:f64) -> Vec<f64> {
        self.model.axial_modes().iter().zip(self.model.axial_slope_scales()).map(|(&mode,&scale)| {
            let phase=mode as f64*core::f64::consts::PI;
            scale*length/phase*det::sin(phase*station)
        }).collect()
    }
    pub(super) fn endpoint_stiffness(&self) -> f64 {
        self.lift_stiffness+2.0*self.model.strain_energy_coefficient()*self.chord[2]*self.chord[2]
    }
    fn count(&self) -> usize { self.omegas.len()+self.model.axial_modes().len()+1 }
    fn slopes(&self,q:&[f64]) -> Vec<f64> {
        let n=self.scales.len();let mut slopes=vec![0.0;2*(n+1)];
        for plane in 0..2 {
            slopes[plane*(n+1)]=self.chord[plane]*q[self.count()-1];
            for mode in 0..n {slopes[plane*(n+1)+mode+1]=self.scales[mode]*q[plane*n+mode];}
        }
        slopes
    }
    fn strain(&self,q:&[f64]) -> Result<Strain,ImpactError> {
        if q.len()!=self.count() || q.iter().any(|v|!v.is_finite()) {return Err(invalid("complete finite supported axial coordinates"));}
        let slopes=self.slopes(q);let mut harmonics=vec![0.0;self.model.axial_modes().len()];
        let mean=self.model.geometric_strain(&slopes,&mut harmonics).map_err(|e|ImpactError::Owner(e.to_string()))?
            +self.chord[2]*q[q.len()-1];
        let mut strain=vec![0.0;harmonics.len()];
        let (energy,spread)=self.model.storage_from_strain(mean,&harmonics,&q[self.omegas.len()..q.len()-1],&mut strain)
            .map_err(|e|ImpactError::Owner(e.to_string()))?;
        Ok(Strain {slopes,mean,harmonics:strain,energy,spread})
    }
    pub(super) fn potential(&self,q:&[f64]) -> f64 {
        self.strain(q).map_or(f64::NAN,|s|s.energy+0.5*self.lift_stiffness*q[q.len()-1].powi(2)
            +self.omegas.iter().zip(q).map(|(w,q)|0.5*w*w*q*q).sum::<f64>())
    }
    fn pull_gradient(&self,slopes:&[f64],mean:f64,harmonics:&[f64],linear:&[f64],out:&mut[f64]) {
        let n=self.scales.len();let c=self.model.strain_energy_coefficient();
        for plane in 0..2 {for mode in 0..n {
            let i=plane*n+mode;out[i]=self.omegas[i].powi(2)*linear[i]+self.scales[mode]*slopes[plane*(n+1)+mode+1];
        }}
        for (i,(&scale,&strain)) in self.model.axial_slope_scales().iter().zip(harmonics).enumerate() {
            out[2*n+i]=c*scale*strain;
        }
        let last=out.len()-1;
        out[last]=self.lift_stiffness*linear[last]+self.chord[0]*slopes[0]
            +self.chord[1]*slopes[n+1]+2.0*c*mean*self.chord[2];
    }
    pub(super) fn gradient(&self,q:&[f64],out:&mut[f64]) {
        let Ok(s)=self.strain(q) else {out.fill(f64::NAN);return;};
        if out.len()!=self.count() {out.fill(f64::NAN);return;}
        let c=self.model.strain_energy_coefficient();let weights:Vec<_>=s.harmonics.iter().map(|v|c*v).collect();
        let mut gradient=vec![0.0;s.slopes.len()];
        if self.model.geometric_strain_gradient(&s.slopes,2.0*c*s.mean,&weights,&mut gradient).is_err() {
            out.fill(f64::NAN);return;
        }
        self.pull_gradient(&gradient,s.mean,&s.harmonics,q,out);
    }
    pub(super) fn hessian_vector(&self,q:&[f64],d:&[f64],out:&mut[f64]) {
        let Ok(s)=self.strain(q) else {out.fill(f64::NAN);return;};
        if d.len()!=self.count() || out.len()!=self.count() || d.iter().any(|v|!v.is_finite()) {out.fill(f64::NAN);return;}
        let ds=self.slopes(d);let count=s.harmonics.len();let c=self.model.strain_energy_coefficient();
        let mut weights=vec![0.0;count];let mut g=vec![0.0;ds.len()];
        // Differentiate the SAME owner's quadratic strain forms; no finite
        // differences or separately reimplemented Fourier products are used.
        if self.model.geometric_strain_gradient(&s.slopes,1.0,&weights,&mut g).is_err() {out.fill(f64::NAN);return;}
        let mean=dot(&g,&ds)+self.chord[2]*d[d.len()-1];
        let mut harmonics=vec![0.0;count];
        for mode in 0..count {
            weights[mode]=1.0;
            if self.model.geometric_strain_gradient(&s.slopes,0.0,&weights,&mut g).is_err() {out.fill(f64::NAN);return;}
            harmonics[mode]=dot(&g,&ds)+self.model.axial_slope_scales()[mode]*d[self.omegas.len()+mode];
            weights[mode]=0.0;
        }
        for (w,value) in weights.iter_mut().zip(&harmonics) {*w=c*value;}
        if self.model.geometric_strain_gradient(&s.slopes,2.0*c*mean,&weights,&mut g).is_err() {out.fill(f64::NAN);return;}
        for (w,value) in weights.iter_mut().zip(&s.harmonics) {*w=c*value;}
        let mut second=vec![0.0;ds.len()];
        if self.model.geometric_strain_gradient(&ds,2.0*c*s.mean,&weights,&mut second).is_err() {out.fill(f64::NAN);return;}
        for (g,h) in g.iter_mut().zip(second) {*g+=h;}
        self.pull_gradient(&g,mean,&harmonics,d,out);
    }
    pub(super) fn observe(&self,q:&[f64]) -> Result<(StringObservation,LongitudinalSupportObservation),ImpactError> {
        let s=self.strain(q)?;let n=self.scales.len();
        let mut slope=det::hypot(s.slopes[0],s.slopes[n+1]);
        for i in 1..=n {slope+=det::hypot(s.slopes[i],s.slopes[n+1+i]);}
        let maximum=s.mean.abs()+s.spread;let minimum_tension=self.tension+self.ea*(s.mean-s.spread);
        if !slope.is_finite() || slope>self.maximum_slope || maximum>self.maximum_strain || minimum_tension<=0.0 {
            return Err(invalid("supported axial strain exceeds its slope/strain envelope or loses local tension"));
        }
        Ok((StringObservation {slope_bound:slope,additional_strain:s.mean,
            tension_n:self.tension+self.ea*s.mean,stretching_energy_j:s.energy},
            LongitudinalSupportObservation {mean_dynamic_strain:s.mean,maximum_abs_dynamic_strain:maximum,
                minimum_tension_n:minimum_tension,maximum_tension_n:self.tension+self.ea*(s.mean+s.spread),strain_energy_j:s.energy}))
    }
}
