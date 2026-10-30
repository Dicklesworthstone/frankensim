//! A finite-mass support carrying selected ends of tensioned filaments.
//!
//! Relative string motion and support motion have kinetic cross terms. Retain
//! them by an exact triangular mass normalization, not by moving contact gaps:
//! w_i(s)=h_i(s)*X+sum(phi_in(s)*q_in), z_in=q_in+c_in*X, z_C=sqrt(M_r)*X,
//! c_in=integral(mu_i*phi_in*h_i ds), M_r=M_C+sum(integral(mu_i*h_i^2))-sum(c_in^2).
//! h=1 carries both ends; h=s/L or 1-s/L fixes the opposite end. Affine
//! lifts add their tension and geometric strain to the ORIGINAL string storage.
//! Bending remains simply supported (zero moment), not clamped at the bridge.
//! Then T=|zdot|^2/2. Existing string energy acts on q=z-c*X; its transpose
//! supplies support reactions. The same transformation applies to contact and
//! relative modal damping. No new constitutive law or time integrator.
use std::{ops::Range, rc::Rc};
use fs_dcontact::Obstacle;
use fs_math::det;
use crate::modal_acoustic_time::ModalAcousticState;
use super::{BodyPotential, ImpactBody, ImpactError, MAX_IMPACT_MODES, invalid};
use super::linear::wire::{LineContact, WireSpan};
use super::string::{StringObservation, StringPotential, StringStretching};

pub mod longitudinal;
use longitudinal::{LongitudinalPotential, LongitudinalSupportModes, LongitudinalSupportObservation,
    StringSupportGeometry};

/// One physical rail/bridge. Its selected wire endpoints share this transverse
/// displacement, mapped into each wire's declared transverse frame.
#[derive(Clone, Copy, Debug)]
pub struct TranslatingSupport {
    /// Rail/actuator effective mass only [kg]. Wire mass is added exactly once.
    pub mass_kg: f64,
    /// Optional ground stiffness [N/m], with rest at X=0.
    pub stiffness_n_m: f64,
    /// Optional ground resistance [N s/m].
    pub damping_n_s_m: f64,
    pub initial_position_m: f64,
    pub initial_velocity_m_s: f64,
    /// Symmetric |X| limit [m]; refusal, not a stop or clamp.
    pub maximum_travel_m: f64,
    /// Absolute wire-slope bound in (0,0.3], including the endpoint lift.
    pub maximum_slope: f64,
}
impl TranslatingSupport {
    pub fn validate(self) -> Result<(), ImpactError> {
        if [self.mass_kg,self.maximum_travel_m,self.maximum_slope].iter().any(|v| !v.is_finite() || *v<=0.0)
            || [self.stiffness_n_m,self.damping_n_s_m].iter().any(|v| !v.is_finite() || *v<0.0)
            || !self.initial_position_m.is_finite() || !self.initial_velocity_m_s.is_finite()
            || self.initial_position_m.abs()>self.maximum_travel_m || self.maximum_slope>0.3 {
            return Err(invalid("moving support requires explicit finite mass, motion, losses and travel/slope bounds"));
        }
        Ok(())
    }
}

/// Aggregate observations, with all energy already included in the host ledger.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SupportObservation {
    pub position_m: f64,
    pub velocity_m_s: f64,
    pub maximum_wire_slope: f64,
    pub maximum_wire_tension_n: f64,
    pub wire_stretching_energy_j: f64,
}

/// Endpoints that follow the common support. Unselected endpoints stay fixed.
/// The same lift applies to transverse and, when retained, axial displacement.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StringEndMotion {
    /// Both ends translate equally: the original rail image, h(s)=1.
    Both,
    /// Only endpoints_m[0] moves: h(s)=1-s/L.
    First,
    /// Only endpoints_m[1] moves: h(s)=s/L (fixed nut, moving bridge).
    Second,
}
impl StringEndMotion {
    fn lift(self, fraction: f64) -> f64 {
        match self { Self::Both => 1.0, Self::First => 1.0-fraction, Self::Second => fraction }
    }
    fn affine(self) -> bool { self != Self::Both }
    fn mass_fraction(self) -> f64 { if self.affine() { 1.0/3.0 } else { 1.0 } }
    fn participation(self, mode: usize, mass: f64) -> f64 {
        let scale = det::sqrt(2.0*mass)/(mode as f64*core::f64::consts::PI);
        match self {
            // Keep original arithmetic on the legacy specialization.
            Self::Both => if mode%2==0 {0.0} else {2.0*det::sqrt(2.0*mass)/(mode as f64*core::f64::consts::PI)},
            Self::First => scale,
            Self::Second => if mode%2==0 {-scale} else {scale},
        }
    }
}

/// Retained transverse displacement components and physical support direction.
/// A two-plane string still has one mass and one geometric extension. Its
/// bridge direction is independent of the direction of any bow or contact.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum StringTransverseMotion {
    /// Original scalar displacement, with positive support motion along it.
    OnePlane,
    /// Both components in one orthonormal transverse frame. The supplied unit
    /// vector maps the common support translation into that frame, verbatim.
    TwoPlanes { support_direction: [f64; 2] },
    /// Both transverse planes and the axial projection of one actual spatial
    /// support direction. Requires explicit longitudinal modes at construction.
    Spatial { geometry: StringSupportGeometry },
}
impl StringTransverseMotion {
    fn planes(self) -> usize {
        match self { Self::OnePlane => 1, Self::TwoPlanes { .. } | Self::Spatial { .. } => 2 }
    }
    fn direction(self) -> [f64; 2] {
        match self { Self::OnePlane => [1.0, 0.0], Self::TwoPlanes { support_direction } => support_direction,
            Self::Spatial { geometry } => [geometry.components()[0],geometry.components()[1]] }
    }
    fn axial_direction(self) -> f64 {
        match self { Self::Spatial { geometry } => geometry.components()[2], _ => 0.0 }
    }
    fn validate(self) -> Result<(), ImpactError> {
        match self { Self::Spatial { .. } => Ok(()), _ => validate_direction(self.direction()) }
    }
}
fn validate_direction(direction: [f64; 2]) -> Result<(), ImpactError> {
    if direction.iter().any(|v| !v.is_finite())
        || (direction[0]*direction[0] + direction[1]*direction[1] - 1.0).abs() > 1e-12 {
        return Err(invalid("supported transverse direction must be a finite unit vector"));
    }
    Ok(())
}

/// Passive drag on projections of ABSOLUTE string velocity, not curvature rate.
/// Original sine projections and their orthogonal endpoint-lift remainder are
/// separate channels. A uniform per-mass rate makes this exactly the spatial
/// local-drag integral; frequency-selected rates retain a modal approximation.
/// These rates must NOT also be included in WireSpan's relative damping.
#[derive(Clone, Debug, PartialEq)]
pub struct AbsoluteStringDrag {
    /// One nonnegative rate [1/s] per retained relative displacement mode,
    /// in plane-major order when both transverse components are retained.
    pub modal_per_s: Vec<f64>,
    /// Nonnegative rate [1/s] for the unresolved affine/constant lift field.
    pub remainder_per_s: f64,
}

#[derive(Clone, Debug)]
struct Data {
    support: TranslatingSupport,
    wires: Vec<WireSpan>,
    potentials: Vec<SupportedPotential>,
    ranges: Vec<Range<usize>>,
    storage_ranges: Vec<Range<usize>>,
    axial_ranges: Vec<Option<Range<usize>>>,
    translation: Vec<f64>,
    damping: Vec<f64>,
    axial_viscosity_n_s: Option<Vec<f64>>,
    finite_strain_axial_viscosity: bool,
    axial_endpoint_damping_n_s_m: f64,
    end_motion: Vec<StringEndMotion>,
    remainder_mass: Vec<f64>,
    absolute_drag: Option<Vec<AbsoluteStringDrag>>,
    transverse: Vec<StringTransverseMotion>,
    weight: f64,
    total_mass_kg: f64,
    // A Gershgorin upper bound for the coupled SMALL-SIGNAL pencil, not its
    // eigenfrequencies. Only the host's cold Nyquist admission consumes this.
    clock_bounds: Vec<f64>,
}

#[derive(Clone, Debug)]
enum SupportedPotential {
    Transverse(StringPotential),
    Longitudinal(LongitudinalPotential),
}
impl SupportedPotential {
    fn endpoint_argument(&self, ends: StringEndMotion) -> bool {
        matches!(self,Self::Longitudinal(_)) || ends.affine()
    }
    fn potential(&self,q:&[f64]) -> f64 {match self {
        Self::Transverse(p)=>p.potential(q),Self::Longitudinal(p)=>p.potential(q)}}
    fn gradient(&self,q:&[f64],out:&mut[f64]) {match self {
        Self::Transverse(p)=>p.gradient(q,out),Self::Longitudinal(p)=>p.gradient(q,out)}}
    fn hessian_vector(&self,q:&[f64],d:&[f64],out:&mut[f64]) {match self {
        Self::Transverse(p)=>p.hessian_vector(q,d,out),Self::Longitudinal(p)=>p.hessian_vector(q,d,out)}}
    fn observe(&self,q:&[f64]) -> Result<StringObservation,ImpactError> {match self {
        Self::Transverse(p)=>p.observe(q),Self::Longitudinal(p)=>p.observe(q).map(|v|v.0)}}
}

/// Immutable shared reduction. Relative string coordinates come first, followed
/// by exactly one support coordinate. Only the physical string laws are reused;
/// no additional string mass, tension spring, prescribed motion or reset.
#[derive(Clone, Debug)]
pub struct SupportedStrings(Rc<Data>);
impl SupportedStrings {
    pub fn new(wires: &[WireSpan], stretching: Option<StringStretching>, support: TranslatingSupport)
        -> Result<Self, ImpactError>
    {
        Self::prepare(wires, support, |_| stretching, |_| StringEndMotion::Both,
            |_| StringTransverseMotion::OnePlane, |_| None)
    }

    /// Use one independent stretching law per wire, in original wire order.
    /// Mass, bending, tension and damping still come from each WireSpan. None
    /// explicitly selects that wire's linear limit; no neighboring wire's EA
    /// or loss is borrowed. The support's slope ceiling also applies to each.
    /// This is the SAME mass transformation, storage and support reaction as new.
    pub fn with_individual_stretching(wires: &[WireSpan], stretching: &[Option<StringStretching>],
        support: TranslatingSupport) -> Result<Self, ImpactError>
    {
        if stretching.len() != wires.len() {
            return Err(invalid("each supported wire requires its own explicit stretching selection"));
        }
        Self::prepare(wires, support, |i| stretching[i], |_| StringEndMotion::Both,
            |_| StringTransverseMotion::OnePlane, |_| None)
    }

    /// Bind actual endpoint kinematics, including fixed-nut/moving-bridge spans.
    /// Each string keeps its independent material stretching. All motion enters
    /// one kinetic normalization and the existing stress-channel potential;
    /// merely moving contact gaps would omit inertia and tension reactions.
    /// WireSpan damping is RELATIVE modal loss; use with_absolute_drag for gas.
    pub fn with_endpoint_motion(wires: &[WireSpan], stretching: &[Option<StringStretching>],
        motion: &[StringEndMotion], support: TranslatingSupport) -> Result<Self, ImpactError>
    {
        if stretching.len()!=wires.len() || motion.len()!=wires.len() {
            return Err(invalid("each supported wire requires explicit stretching and endpoint kinematics"));
        }
        Self::prepare(wires, support, |i| stretching[i], |i| motion[i],
            |_| StringTransverseMotion::OnePlane, |_| None)
    }

    /// Retain one or both physical transverse components on the SAME support.
    /// For direction d, w_alpha=h*d_alpha*X+sum(phi*q_alpha), so the exact
    /// kinetic change is z_alpha=q_alpha+d_alpha*c*X. The Schur complement
    /// retains the unresolved lift mass once; two planes do not double mass.
    /// Both planes feed one EA stress channel, including affine bridge strain.
    /// The shared support has one translation; independent rotation or a second
    /// support translation is not inferred from this representation.
    pub fn with_transverse_motion(wires: &[WireSpan], stretching: &[Option<StringStretching>],
        motion: &[StringEndMotion], transverse: &[StringTransverseMotion],
        support: TranslatingSupport) -> Result<Self, ImpactError>
    {
        if stretching.len()!=wires.len() || motion.len()!=wires.len() || transverse.len()!=wires.len() {
            return Err(invalid("each supported wire requires complete stretching, endpoint and transverse data"));
        }
        Self::prepare(wires, support, |i| stretching[i], |i| motion[i], |i| transverse[i], |_| None)
    }

    /// Exact kinetic coupling of both transverse planes AND retained axial
    /// waves to one oblique translating termination. Geometry supplies every
    /// projection. The solid axial mass and the unresolved affine lift are
    /// counted once; positive finite-strain storage supplies its own reaction.
    /// The reference is already strung and statically balanced. Axial waves are
    /// conservative until explicit axial viscosity is bound with
    /// `with_axial_viscosity`; no axial gas damping is inferred.
    pub fn with_longitudinal_motion(wires:&[WireSpan], stretching:&[Option<StringStretching>],
        motion:&[StringEndMotion], geometry:&[StringSupportGeometry],
        axial:&[LongitudinalSupportModes], support:TranslatingSupport) -> Result<Self,ImpactError> {
        if [stretching.len(),motion.len(),geometry.len(),axial.len()].iter().any(|&n|n!=wires.len()) {
            return Err(invalid("each spatially supported wire requires its full geometry and axial selection"));
        }
        Self::prepare(wires,support,|i|stretching[i],|i|motion[i],
            |i|StringTransverseMotion::Spatial {geometry:geometry[i]},|i|Some(&axial[i]))
    }

    fn prepare<'a>(wires: &[WireSpan], support: TranslatingSupport,
        stretching: impl Fn(usize) -> Option<StringStretching>,
        motion: impl Fn(usize) -> StringEndMotion,
        transverse: impl Fn(usize) -> StringTransverseMotion,
        longitudinal: impl Fn(usize) -> Option<&'a LongitudinalSupportModes>) -> Result<Self, ImpactError>
    {
        support.validate()?;
        let count=wires.iter().enumerate().try_fold(0usize,|n,(i,w)|
            w.damping_per_s.len().checked_mul(transverse(i).planes())
                .and_then(|m|m.checked_add(longitudinal(i).map_or(0,|a|a.modes.len())))
                .and_then(|m| n.checked_add(m)))
            .ok_or_else(||invalid("supported wire mode count overflow"))?;
        if wires.is_empty() || wires.len()>32 || count==0 || count>=MAX_IMPACT_MODES {
            return Err(invalid("supported wires exceed the bounded complete mechanical basis"));
        }
        let mut potentials=Vec::with_capacity(wires.len());let mut ranges=Vec::with_capacity(wires.len());
        let mut storage_ranges=Vec::with_capacity(wires.len());let mut axial_ranges=Vec::with_capacity(wires.len());
        let mut translation=Vec::with_capacity(count);let mut damping=Vec::with_capacity(count);
        let mut end_motion=Vec::with_capacity(wires.len());let mut remainder_mass=Vec::with_capacity(wires.len());
        let mut transverse_motion=Vec::with_capacity(wires.len());
        let mut endpoint_stiffness=0.0;
        let mut total_mass=support.mass_kg;let mut residual_mass=support.mass_kg;let mut frequencies=Vec::with_capacity(count);
        for (index, wire) in wires.iter().enumerate() {
            let selected = stretching(index);
            if let Some(law) = selected { law.validate()?; }
            let law = StringStretching {
                axial_rigidity_n: selected.map_or(0.0, |s| s.axial_rigidity_n),
                maximum_slope: selected.map_or(support.maximum_slope,
                    |s| s.maximum_slope.min(support.maximum_slope)),
            };
            let transverse=transverse(index);
            transverse.validate()?;
            let axial=longitudinal(index);
            if axial.is_some()!=matches!(transverse,StringTransverseMotion::Spatial {..}) {
                return Err(invalid("spatial support projection and axial inertia must be selected together"));
            }
            if let Some(selection)=axial {selection.validate()?;}
            let per_plane=wire.damping_per_s.len();
            let n=per_plane*transverse.planes();
            let initial=vec![ModalAcousticState::default();n];
            let transverse_law=if axial.is_some() {StringStretching {axial_rigidity_n:0.0,..law}} else {law};
            let body=if transverse.planes()==2 {
                wire.stretching_body_two_planes(initial,transverse_law)?
            } else { wire.stretching_body(initial,transverse_law)? };
            let BodyPotential::String(potential)=body.potential else {return Err(invalid("wire preparation changed its storage image"));};
            let ends=motion(index);
            frequencies.extend_from_slice(potential.omegas());
            let potential=if let Some(selection)=axial {
                let StringTransverseMotion::Spatial {geometry}=transverse else {unreachable!()};
                let p=LongitudinalPotential::new(wire,law,selection,geometry,ends,potential.omegas().to_vec())?;
                endpoint_stiffness+=p.endpoint_stiffness();
                frequencies.extend_from_slice(p.axial_omegas());
                SupportedPotential::Longitudinal(p)
            } else {
                let p=if ends.affine() {potential.with_affine_endpoint(wire.length_m())?} else {potential};
                if ends.affine() {endpoint_stiffness+=p.omegas()[n].powi(2);}
                SupportedPotential::Transverse(p)
            };
            end_motion.push(ends);
            transverse_motion.push(transverse);
            let mass=wire.linear_density_kg_m*wire.length_m();
            let start=translation.len();let mut represented=0.0;
            for plane in 0..transverse.planes() {
                for k in 1..=per_plane {
                    let c=ends.participation(k,mass)*transverse.direction()[plane];
                    represented+=c*c;translation.push(c);
                }
            }
            let direction=transverse.direction();
            let transverse_fraction=if axial.is_some() {direction[0]*direction[0]+direction[1]*direction[1]} else {1.0};
            let remainder=mass*transverse_fraction*ends.mass_fraction()-represented;
            if !mass.is_finite() || mass<=0.0 || !remainder.is_finite() || remainder<0.0
                || (transverse_fraction>0.0 && remainder==0.0) {
                return Err(invalid("supported wire translation mass is unrepresentable"));
            }
            // Keep the full lift inertia, including its unresolved sine tail.
            // A single moving endpoint contributes mu*L/3, not mu*L or zero.
            remainder_mass.push(remainder);
            residual_mass+=remainder;total_mass+=mass*transverse_fraction;ranges.push(start..start+n);
            potentials.push(potential);
            for _ in 0..transverse.planes() { damping.extend_from_slice(&wire.damping_per_s); }
            if let Some(selection)=axial {
                let mass=selection.solid_linear_density_kg_m*wire.length_m();
                let axis=transverse.axial_direction();let mut represented=0.0;
                let first=translation.len();
                for &mode in &selection.modes {
                    let c=ends.participation(mode,mass)*axis;represented+=c*c;
                    translation.push(c);damping.push(0.0);
                }
                let remainder=mass*axis*axis*ends.mass_fraction()-represented;
                if !remainder.is_finite() || remainder<0.0 || (axis!=0.0 && remainder==0.0) {
                    return Err(invalid("supported axial lift mass is unrepresentable"));
                }
                residual_mass+=remainder;total_mass+=mass*axis*axis;
                axial_ranges.push(Some(first..translation.len()));
            } else {axial_ranges.push(None);}
            storage_ranges.push(start..translation.len());
        }
        let weight=1.0/det::sqrt(residual_mass);
        if !total_mass.is_finite() || !weight.is_finite() || weight<=0.0 {
            return Err(invalid("moving support mass normalization overflow"));
        }
        let mut support_row=support.stiffness_n_m*weight*weight;let mut bound=0.0_f64;
        if endpoint_stiffness!=0.0 {support_row+=endpoint_stiffness*weight*weight;}
        for (index,r) in ranges.iter().enumerate() {
            let planes=transverse_motion[index].planes();let per_plane=r.len()/planes;
            for mode in 0..per_plane {
                // Degenerate transverse partners can be rotated into the
                // support direction for this eigenvalue bound. Only that
                // combination couples to X; its orthogonal partner keeps k.
                // Summing both coordinate-wise row norms would make clock
                // admission depend on an arbitrary transverse frame.
                let c=if planes==2 {det::hypot(translation[r.start+mode],
                    translation[r.start+per_plane+mode])} else {translation[r.start+mode]};
                let a=c*weight;let omega=frequencies[r.start+mode];let k=omega*omega;
                bound=bound.max(k*(1.0+a.abs()));support_row+=k*(a*a+a.abs());
            }
            if let Some(axial)=&axial_ranges[index] {for mode in axial.clone() {
                let a=translation[mode]*weight;let k=frequencies[mode]*frequencies[mode];
                bound=bound.max(k*(1.0+a.abs()));support_row+=k*(a*a+a.abs());
            }}
        }
        bound=bound.max(support_row);
        if !bound.is_finite() {return Err(invalid("moving support linear clock bound overflow"));}
        let result=Self(Rc::new(Data {support,wires:wires.to_vec(),potentials,ranges,storage_ranges,axial_ranges,translation,damping,
            axial_viscosity_n_s:None,finite_strain_axial_viscosity:false,axial_endpoint_damping_n_s_m:0.0,
            end_motion,remainder_mass,absolute_drag:None,transverse:transverse_motion,weight,
            total_mass_kg:total_mass,clock_bounds:vec![det::sqrt(bound);count+1]}));
        // Check initial energy/motion representation before returning an adapter.
        result.observe_interleaved(&result.initial_state(),0)?;
        Ok(result)
    }
    pub fn mode_count(&self)->usize {self.0.translation.len()+1}
    pub fn wire_count(&self)->usize {self.0.wires.len()}
    pub fn support_coordinate(&self)->usize {self.mode_count()-1}
    /// Generalized force is F * this weight, applied ONLY at support_coordinate.
    pub fn force_weight(&self)->f64 {self.0.weight}
    /// Full effective translating mass in the support direction, not the
    /// one-end lift inertia. Spatial axial/transverse components do not
    /// duplicate the solid mass; transverse added mass stays directional.
    pub fn total_mass_kg(&self)->f64 {self.0.total_mass_kg}
    /// Original endpoint choice, in wire order.
    pub fn endpoint_motion(&self)->&[StringEndMotion] {&self.0.end_motion}
    /// Physical support directions and retained plane counts, in source order.
    pub fn transverse_motion(&self)->&[StringTransverseMotion] {&self.0.transverse}
    /// Retained axial coordinate ranges in original string order.
    pub fn axial_ranges(&self)->&[Option<Range<usize>>] {&self.0.axial_ranges}
    /// Material axial angular frequencies, before endpoint kinetic coupling.
    pub fn axial_angular_frequencies(&self,wire:usize)->Result<&[f64],ImpactError> {
        match self.0.potentials.get(wire) {
            Some(SupportedPotential::Longitudinal(p))=>Ok(p.axial_omegas()),
            _=>Err(invalid("wire has no retained supported axial modes")),
        }
    }
    /// Original absolute-velocity loss declaration; absent on the legacy path.
    pub fn absolute_drag(&self)->Option<&[AbsoluteStringDrag]> {self.0.absolute_drag.as_deref()}

    /// Add absolute-velocity drag once to this immutable cold image. Bending
    /// viscosity stays on the relative curvature velocity; gas drag must be
    /// removed from those rates by its source owner before selecting this path.
    /// z_i=q_i+c_i X are exactly the mass-normalized absolute sine projections.
    /// The remaining lift is orthogonal, so its resistance is rate*M_rem/M_r
    /// on the shared coordinate. No drag cross terms or unresolved mass vanish.
    pub fn with_absolute_drag(mut self, drag: Vec<AbsoluteStringDrag>) -> Result<Self, ImpactError> {
        if self.0.absolute_drag.is_some() || drag.len()!=self.wire_count() {
            return Err(invalid("absolute string drag requires one cold declaration per wire"));
        }
        let mut support_rate=0.0;
        for ((law, range), mass) in drag.iter().zip(&self.0.ranges).zip(&self.0.remainder_mass) {
            if law.modal_per_s.len()!=range.len() || !law.remainder_per_s.is_finite()
                || law.remainder_per_s<0.0 || law.modal_per_s.iter().any(|v| !v.is_finite() || *v<0.0) {
                return Err(invalid("absolute drag needs every finite nonnegative modal and remainder rate"));
            }
            support_rate+=law.remainder_per_s*mass*self.0.weight*self.0.weight;
        }
        if !support_rate.is_finite() {return Err(invalid("absolute support drag overflows"));}
        Rc::make_mut(&mut self.0).absolute_drag=Some(drag);
        Ok(self)
    }
    pub fn clock_bounds(&self)->&[f64] {&self.0.clock_bounds}
    fn initial_state(&self)->Vec<f64> {
        let mut state=vec![0.0;2*self.mode_count()];let s=self.0.support;
        for (i,&c) in self.0.translation.iter().enumerate() {
            state[2*i]=c*s.initial_position_m;state[2*i+1]=c*s.initial_velocity_m_s;
        }
        let i=self.support_coordinate();state[2*i]=s.initial_position_m/self.0.weight;
        state[2*i+1]=s.initial_velocity_m_s/self.0.weight;state
    }
    /// Initial relative modes are at rest; the support and its endpoint lifts
    /// have the declared motion. Relative internal damping is transformed
    /// by the host; adding diagonal modal drag here would double-count it.
    pub fn body(&self)->ImpactBody {
        let x=self.initial_state();
        ImpactBody {potential:BodyPotential::SupportedStrings(self.clone()),
            initial:x.chunks_exact(2).map(|s|ModalAcousticState {
                displacement_m_sqrt_kg:s[0],velocity_m_sqrt_kg_per_s:s[1]}).collect(),
            damping_per_s:vec![0.0;self.mode_count()]}
    }
    /// Initialize independent relative string motion without dropping the rail's
    /// declared translation or velocity. Input order is wire order, then each
    /// wire's transverse planes then its selected axial modes. This is a cold coordinate change,
    /// not a state reset, equilibrium solve or running-history transfer.
    /// Kinetic cross terms are kept by z_i=q_i+c_i*X and their velocity analogue.
    /// The host still owns time integration and the complete energy ceiling.
    pub fn body_with_relative_states(&self, relative: &[ModalAcousticState])
        -> Result<ImpactBody, ImpactError>
    {
        if relative.len() != self.support_coordinate() || relative.iter().any(|s|
            !s.displacement_m_sqrt_kg.is_finite() || !s.velocity_m_sqrt_kg_per_s.is_finite()) {
            return Err(invalid("supported initialization needs every finite relative wire coordinate"));
        }
        let mut body = self.body();
        for (absolute, local) in body.initial.iter_mut().zip(relative) {
            absolute.displacement_m_sqrt_kg += local.displacement_m_sqrt_kg;
            absolute.velocity_m_sqrt_kg_per_s += local.velocity_m_sqrt_kg_per_s;
        }
        let state: Vec<_> = body.initial.iter().flat_map(|s|
            [s.displacement_m_sqrt_kg, s.velocity_m_sqrt_kg_per_s]).collect();
        self.observe_interleaved(&state, 0)?;
        let q: Vec<_> = body.initial.iter().map(|s| s.displacement_m_sqrt_kg).collect();
        let energy = self.potential(&q) + body.initial.iter()
            .map(|s| 0.5 * s.velocity_m_sqrt_kg_per_s.powi(2)).sum::<f64>();
        if !energy.is_finite() {
            return Err(invalid("supported initial kinetic or potential energy is unrepresentable"));
        }
        Ok(body)
    }

    /// Work-conjugate ABSOLUTE point displacement/force row [1/sqrt(kg)].
    /// Its dot product with accepted normalized velocity is physical point
    /// velocity. Its transpose applies a signed force without omitting the
    /// support reaction. Only the declared endpoints follow the common support.
    /// Interior sine shapes come from the original string-collocation owner.
    /// This convenience method observes the first transverse plane; use
    /// `absolute_shape_along` for an independent physical contact direction.
    pub fn absolute_shape(&self, wire: usize, station_fraction: f64)
        -> Result<Vec<f64>, ImpactError>
    {
        self.absolute_shape_along(wire, station_fraction, [1.0, 0.0])
    }

    /// Absolute displacement and transposed force in a physical transverse
    /// direction, including the common support reaction. A one-plane source
    /// refuses a nonzero missing component instead of projecting it away.
    pub fn absolute_shape_along(&self, wire: usize, station_fraction: f64,
        direction: [f64; 2]) -> Result<Vec<f64>, ImpactError>
    {
        let r = self.0.ranges.get(wire).ok_or_else(|| invalid("unknown supported wire"))?;
        validate_direction(direction)?;
        let transverse=self.0.transverse[wire];
        if transverse.planes()==1 && direction[1]!=0.0 {
            return Err(invalid("a one-plane supported string cannot carry a second-plane force"));
        }
        if !station_fraction.is_finite() || !(0.0..=1.0).contains(&station_fraction) {
            return Err(invalid("supported point station must lie on its physical span"));
        }
        let support_direction=transverse.direction();
        let projection=direction[0]*support_direction[0]+direction[1]*support_direction[1];
        let mut row = vec![0.0; self.mode_count()];
        if station_fraction == 0.0 || station_fraction == 1.0 {
            row[self.support_coordinate()] = projection*self.0.end_motion[wire].lift(station_fraction)*self.force_weight();
            return Ok(row);
        }
        let source = &self.0.wires[wire];
        let per_plane=source.damping_per_s.len();
        let shape = fs_dcontact::string_collocation(source.length_m(), source.linear_density_kg_m,
            &[station_fraction * source.length_m()], per_plane)
            .map_err(|e| ImpactError::Owner(e.to_string()))?;
        for plane in 0..transverse.planes() {
            for (k,value) in shape.iter().enumerate() {
                row[r.start+plane*per_plane+k]=direction[plane]*value;
            }
        }
        let cross = row[r.clone()].iter().zip(&self.0.translation[r.clone()]).map(|(b,c)| b*c).sum::<f64>();
        row[self.support_coordinate()] = (projection*self.0.end_motion[wire].lift(station_fraction) - cross) * self.force_weight();
        if row.iter().any(|v| !v.is_finite()) {
            return Err(invalid("supported absolute point projection overflows"));
        }
        Ok(row)
    }

    /// Absolute axial point displacement/force row, including the SAME moving
    /// support reaction. Sparse axial selections retain their actual indices.
    pub fn absolute_axial_shape(&self,wire:usize,station_fraction:f64)->Result<Vec<f64>,ImpactError> {
        let Some(SupportedPotential::Longitudinal(p))=self.0.potentials.get(wire) else {
            return Err(invalid("axial point motion requires retained longitudinal inertia"));
        };
        if !station_fraction.is_finite() || !(0.0..=1.0).contains(&station_fraction) {
            return Err(invalid("supported axial station must lie on its physical span"));
        }
        let range=self.0.axial_ranges[wire].as_ref().ok_or_else(||invalid("missing supported axial basis"))?;
        let mut row=vec![0.0;self.mode_count()];
        if station_fraction!=0.0 && station_fraction!=1.0 {
            row[range.clone()].copy_from_slice(&p.axial_shape(station_fraction,self.0.wires[wire].length_m()));
        }
        let cross=row[range.clone()].iter().zip(&self.0.translation[range.clone()]).map(|(b,c)|b*c).sum::<f64>();
        row[self.support_coordinate()]=(self.0.transverse[wire].axial_direction()
            *self.0.end_motion[wire].lift(station_fraction)-cross)*self.force_weight();
        if row.iter().any(|v|!v.is_finite()) {return Err(invalid("supported axial point projection overflows"));}
        Ok(row)
    }

    fn relative(&self,z:&[f64],q:&mut[f64])->bool {
        if z.len()!=self.mode_count() || q.len()!=self.support_coordinate()
            || z.iter().any(|v|!v.is_finite()) {return false;}
        let x=z[self.support_coordinate()]*self.0.weight;
        for ((q,z),c) in q.iter_mut().zip(z).zip(&self.0.translation) {*q=*z-c*x;}
        q.iter().all(|v|v.is_finite())
    }
    // Extend each original relative modal slice by physical endpoint motion
    // only when needed. The extra coordinate is a potential argument, NOT a
    // duplicated body or an added kinetic-energy degree of freedom.
    fn local_coordinates(&self, wire: usize, q: &[f64], x: f64,
        out: &mut [f64; MAX_IMPACT_MODES]) -> usize {
        let range=&self.0.storage_ranges[wire];let n=range.len();
        out[..n].copy_from_slice(&q[range.clone()]);
        if self.0.potentials[wire].endpoint_argument(self.0.end_motion[wire]) {out[n]=x;n+1} else {n}
    }
    pub(super) fn potential(&self,z:&[f64])->f64 {
        let mut q=[0.0;MAX_IMPACT_MODES];let n=self.support_coordinate();
        if !self.relative(z,&mut q[..n]) {return f64::NAN;}
        let x=z[n]*self.0.weight;
        let mut energy=0.5*self.0.support.stiffness_n_m*x*x;
        let mut local=[0.0;MAX_IMPACT_MODES];
        for (i,p) in self.0.potentials.iter().enumerate() {
            let count=self.local_coordinates(i,&q,x,&mut local);
            energy+=p.potential(&local[..count]);
        }
        energy
    }
    pub(super) fn gradient(&self,z:&[f64],out:&mut[f64]) {
        let mut q=[0.0;MAX_IMPACT_MODES];let n=self.support_coordinate();
        if out.len()!=self.mode_count() || !self.relative(z,&mut q[..n]) {out.fill(f64::NAN);return;}
        let (mut local,mut g)=([0.0;MAX_IMPACT_MODES],[0.0;MAX_IMPACT_MODES]);
        let mut endpoint=0.0;
        for (i,(p,r)) in self.0.potentials.iter().zip(&self.0.storage_ranges).enumerate() {
            let count=self.local_coordinates(i,&q,z[n]*self.0.weight,&mut local);
            p.gradient(&local[..count],&mut g[..count]);
            out[r.clone()].copy_from_slice(&g[..r.len()]);
            if p.endpoint_argument(self.0.end_motion[i]) {endpoint+=g[r.len()];}
        }
        out[n]=self.0.support.stiffness_n_m*z[n]*self.0.weight*self.0.weight
            -self.0.weight*self.0.translation.iter().zip(&out[..n]).map(|(c,g)|c*g).sum::<f64>();
        if endpoint!=0.0 {out[n]+=self.0.weight*endpoint;}
    }
    pub(super) fn hessian_vector(&self,z:&[f64],d:&[f64],out:&mut[f64]) {
        let (mut q,mut v)=([0.0;MAX_IMPACT_MODES],[0.0;MAX_IMPACT_MODES]);let n=self.support_coordinate();
        if out.len()!=self.mode_count() || !self.relative(z,&mut q[..n]) || !self.relative(d,&mut v[..n]) {
            out.fill(f64::NAN);return;
        }
        let (mut local,mut direction,mut hv)=([0.0;MAX_IMPACT_MODES],[0.0;MAX_IMPACT_MODES],[0.0;MAX_IMPACT_MODES]);
        let mut endpoint=0.0;
        for (i,(p,r)) in self.0.potentials.iter().zip(&self.0.storage_ranges).enumerate() {
            let count=self.local_coordinates(i,&q,z[n]*self.0.weight,&mut local);
            self.local_coordinates(i,&v,d[n]*self.0.weight,&mut direction);
            p.hessian_vector(&local[..count],&direction[..count],&mut hv[..count]);
            out[r.clone()].copy_from_slice(&hv[..r.len()]);
            if p.endpoint_argument(self.0.end_motion[i]) {endpoint+=hv[r.len()];}
        }
        out[n]=self.0.support.stiffness_n_m*d[n]*self.0.weight*self.0.weight
            -self.0.weight*self.0.translation.iter().zip(&out[..n]).map(|(c,g)|c*g).sum::<f64>();
        if endpoint!=0.0 {out[n]+=self.0.weight*endpoint;}
    }
    /// Preserve the existing line-contact law and receiver row. Its relative
    /// wire row is transformed to ABSOLUTE motion before the implicit solve.
    #[allow(clippy::too_many_arguments)]
    pub fn contact(&self,wire:usize,line:&LineContact,receiver_shapes:&[Vec<f64>],
        receiver:Range<usize>,first:usize,total:usize)->Result<Obstacle,ImpactError>
    {
        self.contact_along(wire,line,receiver_shapes,receiver,first,total,[1.0,0.0])
    }
    /// The same distributed contact law along a declared physical direction.
    /// Receiver rows must observe that same direction in their own basis.
    #[allow(clippy::too_many_arguments)]
    pub fn contact_along(&self,wire:usize,line:&LineContact,receiver_shapes:&[Vec<f64>],
        receiver:Range<usize>,first:usize,total:usize,direction:[f64;2])->Result<Obstacle,ImpactError>
    {
        let r=self.0.ranges.get(wire).ok_or_else(||invalid("unknown supported wire"))?;
        validate_direction(direction)?;
        let transverse=self.0.transverse[wire];
        if transverse.planes()==1 && direction[1]!=0.0 {
            return Err(invalid("a one-plane supported string cannot carry a second-plane contact"));
        }
        let end=first.checked_add(self.mode_count()).ok_or_else(||invalid("supported wire address overflow"))?;
        if end>total || receiver.start<end && first<receiver.end {
            return Err(invalid("supported wire and receiver ranges must be disjoint and complete"));
        }
        // The original wire owner admits geometry, contact law, receiver and
        // quadrature. Replace its scalar point row by the SAME directed row
        // used by physical force and velocity observation in this reduction.
        let per_plane=self.0.wires[wire].damping_per_s.len();
        let ob=self.0.wires[wire].contact(line,receiver_shapes,receiver,first+r.start..first+r.start+per_plane,total)?;
        let mut rows=ob.collocation().to_vec();
        for (row,station) in rows.chunks_exact_mut(total).zip(&line.stations_m) {
            let original=row[first+r.start..first+r.start+per_plane].to_vec();
            for plane in 0..transverse.planes() {
                for (k,value) in original.iter().enumerate() {
                    row[first+r.start+plane*per_plane+k]=direction[plane]*value;
                }
            }
            let cross=row[first+r.start..first+r.end].iter().zip(&self.0.translation[r.clone()])
                .map(|(b,c)|b*c).sum::<f64>();
            let support_direction=transverse.direction();
            let projection=direction[0]*support_direction[0]+direction[1]*support_direction[1];
            let lift=self.0.end_motion[wire].lift(station/self.0.wires[wire].length_m());
            row[first+self.support_coordinate()]=(-projection*lift-cross)*self.0.weight;
        }
        Obstacle::new(rows,ob.n_points(),total,ob.gaps().to_vec(),ob.weights().to_vec(),
            ob.stiffness(),ob.alpha(),ob.provenance().into())
            .and_then(|o|o.with_internal_loss(ob.internal_loss())).map_err(|e|ImpactError::Owner(e.to_string()))
    }
    // Internal relative modal loss is part of the body, not a collection of
    // user-added localized mufflers. Do not consume/increase that separate cap.
    pub(super) fn add_resistance(&self,first:usize,dim:usize,r:&mut[f64])->Result<(),ImpactError> {
        let n=self.support_coordinate();let p=2*(first+n)+1;
        if p>=dim || r.len()!=dim*dim {return Err(invalid("supported damping dimensions"));}
        for i in 0..n {
            let j=2*(first+i)+1;let a=self.0.translation[i]*self.0.weight;let drag=self.0.damping[i];
            r[j*dim+j]+=drag;r[j*dim+p]-=drag*a;r[p*dim+j]-=drag*a;r[p*dim+p]+=drag*a*a;
        }
        r[p*dim+p]+=self.0.support.damping_n_s_m*self.0.weight*self.0.weight;
        if self.0.axial_endpoint_damping_n_s_m!=0.0 {
            r[p*dim+p]+=self.0.axial_endpoint_damping_n_s_m*self.0.weight*self.0.weight;
        }
        if let Some(drag)=&self.0.absolute_drag {
            for ((law,range),mass) in drag.iter().zip(&self.0.ranges).zip(&self.0.remainder_mass) {
                for (i,rate) in range.clone().zip(&law.modal_per_s) {
                    let j=2*(first+i)+1;r[j*dim+j]+=rate;
                }
                r[p*dim+p]+=law.remainder_per_s*mass*self.0.weight*self.0.weight;
            }
        }
        if r.iter().any(|v|!v.is_finite()) {return Err(invalid("supported damping overflow"));}Ok(())
    }
    pub fn string_observation(&self,wire:usize,z:&[f64])->Result<StringObservation,ImpactError> {
        self.0.ranges.get(wire).ok_or_else(||invalid("unknown supported wire"))?;
        let mut q=[0.0;MAX_IMPACT_MODES];
        if !self.relative(z,&mut q[..self.support_coordinate()]) {return Err(invalid("invalid supported wire state"));}
        let mut local=[0.0;MAX_IMPACT_MODES];
        let count=self.local_coordinates(wire,&q,z[self.support_coordinate()]*self.0.weight,&mut local);
        self.0.potentials[wire].observe(&local[..count])
    }
    /// Continuous-span strain and tension bounds, including actual endpoint
    /// extension and the retained propagating axial strains.
    pub fn longitudinal_observation(&self,wire:usize,z:&[f64])->Result<LongitudinalSupportObservation,ImpactError> {
        let Some(SupportedPotential::Longitudinal(p))=self.0.potentials.get(wire) else {
            return Err(invalid("wire has no supported longitudinal strain image"));
        };
        let mut q=[0.0;MAX_IMPACT_MODES];
        if !self.relative(z,&mut q[..self.support_coordinate()]) {return Err(invalid("invalid supported axial state"));}
        let mut local=[0.0;MAX_IMPACT_MODES];
        let count=self.local_coordinates(wire,&q,z[self.support_coordinate()]*self.0.weight,&mut local);
        p.observe(&local[..count]).map(|v|v.1)
    }
    pub fn observe_interleaved(&self,state:&[f64],first:usize)->Result<SupportObservation,ImpactError> {
        let end=first.checked_add(self.mode_count()).and_then(|n|n.checked_mul(2))
            .ok_or_else(||invalid("supported observation address overflow"))?;
        if state.len()<end || state[2*first..end].iter().any(|v|!v.is_finite()) {
            return Err(invalid("supported observation requires complete finite motion"));
        }
        let n=self.support_coordinate();let x=state[2*(first+n)]*self.0.weight;
        if x.abs()>self.0.support.maximum_travel_m {return Err(invalid("support travel exceeds its declared validity"));}
        let mut q=[0.0;MAX_IMPACT_MODES];
        for i in 0..n {q[i]=state[2*(first+i)]-self.0.translation[i]*x;}
        let mut result=SupportObservation {position_m:x,velocity_m_s:state[2*(first+n)+1]*self.0.weight,
            ..SupportObservation::default()};
        let mut local=[0.0;MAX_IMPACT_MODES];
        for (i,p) in self.0.potentials.iter().enumerate() {
            let count=self.local_coordinates(i,&q,x,&mut local);
            let (observation,maximum_tension)=match p {
                SupportedPotential::Transverse(p)=>{let o=p.observe(&local[..count])?;(o,o.tension_n)},
                SupportedPotential::Longitudinal(p)=>{let (o,b)=p.observe(&local[..count])?;(o,b.maximum_tension_n)},
            };
            result.maximum_wire_slope=result.maximum_wire_slope.max(observation.slope_bound);
            result.maximum_wire_tension_n=result.maximum_wire_tension_n.max(maximum_tension);
            result.wire_stretching_energy_j+=observation.stretching_energy_j;
        }
        if !result.wire_stretching_energy_j.is_finite() {return Err(invalid("supported wire energy overflow"));}
        Ok(result)
    }
}

#[cfg(test)]
#[path = "supported_tests.rs"]
mod tests;

#[cfg(test)]
#[path = "supported_heterogeneous_tests.rs"]
mod heterogeneous_tests;

#[cfg(test)]
#[path = "supported_bridge_tests.rs"]
mod bridge_tests;

#[cfg(test)]
#[path = "supported_longitudinal_tests.rs"]
mod longitudinal_tests;
