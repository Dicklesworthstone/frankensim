//! Actual finite microphone geometry and Green-field observation.
//! One prepared quadrature row per receiver/frequency serves every source mode.
//! The unchanged boundary solution also supplies optional radiation reaction.
use super::{Error, Receiver, Medium, RadiationSolution, SpherePanels, C64,
    receiver_response, shift_to_enclosing_sphere};
use fs_bem::near_field::{FirstOrder, Geometry, Options, Prepared};
use fs_bem::reciprocal_radiation as reciprocal;
use fs_exec::CancelGate;

#[path="microphone_spec.rs"]
pub mod input;

/// Explicit near-field accuracy/work and physical microphone selection. This
/// is an ideal pressure/particle-velocity observer, not a measured capsule.
#[derive(Clone, Copy, Debug)]
pub struct Microphone {
    position_m: [f64; 3],
    minimum_clearance_m: f64,
    quadrature: Options,
    pattern: FirstOrder,
}
impl Microphone {
    pub fn new(position_m: [f64; 3], minimum_clearance_m: f64, quadrature: Options,
        pattern: FirstOrder) -> Result<Self, Error>
    {
        // Match the numerical owner's cold envelope, before a mesh is available.
        // Geometry::new and prepare still perform the authoritative admission.
        if position_m.iter().any(|v| !v.is_finite() || v.abs()>1e9)
            || !minimum_clearance_m.is_finite() || minimum_clearance_m<=0.0
            || !quadrature.relative_tolerance.is_finite()
            || !(1e-12..=1e-3).contains(&quadrature.relative_tolerance)
            || quadrature.maximum_depth>16
            || !(80..=20_000_000).contains(&quadrature.maximum_kernel_evaluations) {
            return Err("near microphone requires finite SI position, positive clearance and bounded quadrature controls".into());
        }
        Ok(Self {position_m, minimum_clearance_m, quadrature, pattern})
    }
    pub fn position_m(self) -> [f64;3] {self.position_m}
    pub fn pressure_fraction(self) -> f64 {self.pattern.pressure_fraction()}
}

// Explicit delay lines require >=2 output samples. For a closer receiver peel
// NOTHING: its entire physical propagation remains in the fitted Green transfer.
// Zero here is not zero flight time, a relocated receiver or an added latency.
fn peeled_delay(clearance: f64, medium: Medium, dt: f64) -> Result<f64, Error> {
    if !dt.is_finite() || dt<=0.0 || !medium.sound_speed.is_finite() || medium.sound_speed<=0.0
        || !medium.density.is_finite() || medium.density<=0.0
        || !clearance.is_finite() || clearance<=0.0 {
        return Err("near microphone needs a finite positive medium, clearance and output clock".into());
    }
    let flight=clearance/medium.sound_speed;
    if !flight.is_finite() {return Err("near microphone flight time overflow".into());}
    Ok(if flight>=2.0*dt {flight} else {0.0})
}

/// Cold admitted observers; no physical state, source gain or runtime filters.
pub(super) struct Scene<'a> {
    surface: &'a SpherePanels,
    receivers: &'a [Receiver],
    geometry: Vec<Option<Geometry<'a>>>,
    delays: Vec<f64>,
    gains: Vec<f64>,
    medium: Medium,
    radius: f64,
}
impl<'a> Scene<'a> {
    pub(super) fn new(surface: &'a SpherePanels, receivers: &'a [Receiver], radius: f64,
        medium: Medium, dt: f64, gate: &CancelGate) -> Result<Self, Error>
    {
        if !(1..=2).contains(&receivers.len()) {return Err("one or two microphone receivers are required".into());}
        let mut geometry=Vec::with_capacity(receivers.len());
        let mut delays=Vec::with_capacity(receivers.len());
        let mut gains=Vec::with_capacity(receivers.len());
        for &receiver in receivers {
            if gate.is_requested() {return Err("microphone geometry preparation cancelled".into());}
            match receiver {
                Receiver::NearField(mic)=>{
                    let g=Geometry::new(surface, &[mic.position_m], mic.minimum_clearance_m)?;
                    let clearance=g.clearances_m()[0];
                    delays.push(peeled_delay(clearance,medium,dt)?);gains.push(1.0);
                    eprintln!("near microphone: position_m={:?}, clearance_m={clearance}, alpha={}, front={:?}, quadrature_tolerance={}, depth={}, maximum_kernel_evaluations_per_frequency={}; Pa-equivalent, no capsule/electronics model",
                        mic.position_m,mic.pattern.pressure_fraction(),mic.pattern.front_axis(),
                        mic.quadrature.relative_tolerance,mic.quadrature.maximum_depth,mic.quadrature.maximum_kernel_evaluations);
                    geometry.push(Some(g));
                }
                _=>{let (delay,gain)=receiver.propagation(radius,medium,dt)?;
                    delays.push(delay);gains.push(gain);geometry.push(None);}
            }
        }
        Ok(Self {surface, receivers, geometry, delays, gains, medium, radius})
    }
    pub(super) fn delay(&self, channel: usize) -> f64 {self.delays[channel]}
    pub(super) fn gain(&self, channel: usize) -> f64 {self.gains[channel]}
    pub(super) fn prepare(&self, k: f64, gate: &CancelGate) -> Result<Frequency<'_, 'a>, Error> {
        let mut rows=Vec::with_capacity(self.receivers.len());
        for (&receiver, geometry) in self.receivers.iter().zip(&self.geometry) {
            if gate.is_requested() {return Err("microphone quadrature preparation cancelled".into());}
            let prepared=match (receiver,geometry) {
                (Receiver::NearField(mic),Some(g))=>Some(if mic.pressure_fraction()==1.0 {
                    g.prepare(k,self.medium,mic.quadrature)?
                }else{g.prepare_velocity(k,self.medium,mic.quadrature)?}),
                (_,None)=>None,
                _=>return Err("microphone geometry and observation selection disagree".into()),
            };
            rows.push(prepared);
        }
        Ok(Frequency {scene:self, rows})
    }
    /// Bind adaptive receivers to the actual P1/P0 source once for the sweep.
    /// Existing gains and causal phase peeling remain those admitted above.
    pub(super) fn reciprocal<'s,'p>(&'s self, owner:&'p reciprocal::Prepared<'a>,
        gate:&CancelGate)->Result<ReciprocalScene<'s,'p,'a>,Error> {
        let mut geometry=Vec::with_capacity(self.receivers.len());
        for &receiver in self.receivers {
            if gate.is_requested() {return Err("reciprocal microphone geometry cancelled".into());}
            geometry.push(match receiver {
                Receiver::NearField(mic)=>Some(owner.receivers(&[mic.position_m],mic.minimum_clearance_m)?),
                _=>None,
            });
        }
        Ok(ReciprocalScene {scene:self,geometry})
    }
}

pub(super) struct Frequency<'s, 'a> {scene:&'s Scene<'a>, rows:Vec<Option<Prepared<'a>>>}
impl Frequency<'_, '_> {
    pub(super) fn response(&self, channel:usize, solution:&RadiationSolution) -> Result<C64,Error> {
        let s=self.scene;
        match (s.receivers[channel],&self.rows[channel]) {
            (Receiver::NearField(mic),Some(row))=>{
                let pressure=if mic.pressure_fraction()==1.0 {
                    row.evaluate(solution)?.pressure[0]
                }else{
                    let fields=row.evaluate_velocity(solution)?;
                    mic.pattern.observe(fields.scalar.pressure[0],fields.particle_velocity_m_s[0],s.medium)?
                };
                Ok(shift_to_enclosing_sphere(pressure,solution.k*s.medium.sound_speed,-s.delays[channel]))
            }
            (_,None)=>receiver_response(s.surface,solution,s.medium,s.receivers[channel],s.radius,s.delays[channel]),
            _=>Err("unprepared microphone response".into()),
        }
    }
}

pub(super) struct ReciprocalScene<'s,'p,'a> {
    scene:&'s Scene<'a>,geometry:Vec<Option<reciprocal::ReceiverGeometry<'p,'a>>>,
}
impl<'s,'p,'a> ReciprocalScene<'s,'p,'a> {
    /// Refine the near-field kernels once per frequency, not per source RHS.
    /// No P0 observer is prepared for this source representation.
    pub(super) fn prepare(&self,k:f64,gate:&CancelGate)->Result<ReciprocalFrequency<'s,'p,'a>,Error> {
        let mut rows=Vec::with_capacity(self.scene.receivers.len());
        for (&receiver,geometry) in self.scene.receivers.iter().zip(&self.geometry) {
            if gate.is_requested() {return Err("reciprocal microphone integration cancelled".into());}
            rows.push(match (receiver,geometry) {
                (Receiver::NearField(mic),Some(g))=>Some(g.prepare(k,self.scene.medium,mic.quadrature,
                    mic.pressure_fraction()!=1.,||gate.is_requested())?),
                (_,None)=>None,
                _=>return Err("reciprocal microphone geometry and source disagree".into()),
            });
        }
        Ok(ReciprocalFrequency {scene:self.scene,rows})
    }
}
pub(super) struct ReciprocalFrequency<'s,'p,'a> {
    scene:&'s Scene<'a>,rows:Vec<Option<reciprocal::PreparedReceivers<'p,'a>>>,
}
impl ReciprocalFrequency<'_,'_,'_> {
    pub(super) fn response(&self,channel:usize,solution:&reciprocal::Solution<'_,'_>)
        ->Result<C64,Error> {
        let scene=self.scene;let w=solution.k*scene.medium.sound_speed;
        match (scene.receivers[channel],&self.rows[channel]) {
            (Receiver::NearField(mic),Some(row))=>{
                let pressure=if mic.pressure_fraction()==1. {
                    row.evaluate(solution)?.pressure[0]
                }else{
                    let field=row.evaluate_velocity(solution)?;
                    mic.pattern.observe(field.scalar.pressure[0],field.particle_velocity_m_s[0],scene.medium)?
                };
                Ok(shift_to_enclosing_sphere(pressure,w,-scene.delays[channel]))
            }
            (Receiver::FarField(position),None)=>Ok(shift_to_enclosing_sphere(
                solution.far_field(position)?,w,scene.radius/scene.medium.sound_speed)),
            (Receiver::FinitePoint(position),None)=>Ok(shift_to_enclosing_sphere(
                solution.pressure_at(position)?,w,-scene.delays[channel])),
            _=>Err("unprepared reciprocal microphone response".into()),
        }
    }
}

#[cfg(test)]
#[path="receivers/tests.rs"]
mod tests;
