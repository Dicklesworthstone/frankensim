//! Material-derived nonlinear filaments on ONE inertially coupled support.
//! The same material section lowers EA, EI, mass and loss; SupportedStrings
//! owns the kinetic transformation and reactions; ImpactSystem owns every step.
//! No string is rebuilt from a pitch, averaged with a neighbor or muted.
//!
//! Endpoint motion and gas frame are explicit: the legacy co-moving rail carries
//! both ends, while stationary-gas construction permits a fixed nut and moving
//! bridge. Both transverse components use the source's physical bridge direction
//! and share one stretching channel. Material curvature loss and absolute-velocity gas loss remain distinct.
//! Gas loss is a passive spectral approximation, not unsteady fluid added mass
//! or full broadband hydrodynamics. Support rotation, live material changes,
//! winding slip, torsional gas drag and a bow friction law remain absent.
use std::ops::Range;
use super::{NonlinearMaterialString, BodyPotential, ImpactBody, ImpactError, MAX_IMPACT_MODES, invalid};
use crate::render::plate::impact::{ImpactConfig, ImpactFrame, ImpactSystem, PreparedImpactSystem};
use crate::render::plate::impact::string::StringStretching;
use crate::render::plate::impact::supported::{SupportedStrings, TranslatingSupport, StringEndMotion,
    StringTransverseMotion};
use crate::bowed_string::runtime::reciprocal::bank::{MaterialStringSource,
    polarization::TransverseDirection};
use fs_exec::CancelGate;

mod stationary;
pub mod longitudinal;
use longitudinal::MaterialLongitudinalSupport;

/// Explicit physical frame and approximation used for string/gas resistance.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SupportGasFrame {
    /// Original modal resistance relative to a gas translating with both ends.
    CoMoving,
    /// Stationary gas: original natural-frequency resistances act on absolute
    /// sine velocity projections. The orthogonal lift remainder uses the SAME
    /// cylinder law at this explicit frequency. It has no fictitious zero drag.
    /// This spectral approximation requires modal/time/frequency refinement;
    /// it is not an exact causal hydrodynamic impedance for arbitrary motion.
    StationaryModal {
        /// Evaluation frequency [Hz] for the unresolved endpoint-lift field.
        remainder_frequency_hz: f64,
    },
}

/// Immutable source geometry/materials and their complete coupled initial image.
/// Body zero contains each string's bending and selected axial coordinates,
/// followed by the rail coordinate.
/// Requested torsion bodies follow in original string order. Linear torsional
/// inertia is independent of a translating, nonrotating rail; it is not dropped.
#[derive(Clone, Debug)]
pub struct SupportedMaterialStrings {
    sources: Vec<NonlinearMaterialString>,
    support: TranslatingSupport,
    supported: SupportedStrings,
    bodies: Vec<ImpactBody>,
    torsion: Vec<Option<Range<usize>>>,
    mode_count: usize,
    sample_rate_hz: u32,
    gas_frame: SupportGasFrame,
    longitudinal: Option<Vec<MaterialLongitudinalSupport>>,
    axial_material_loss: Option<Vec<fs_material::section::axial::ConcentricAxialViscosity>>,
    axial_material_relaxation: Option<longitudinal::relaxation::SupportedMaterialRelaxation>,
}
fn poll(gate: &CancelGate) -> Result<(), ImpactError> {
    if gate.is_requested() { Err(ImpactError::Cancelled) } else { Ok(()) }
}
impl SupportedMaterialStrings {
    /// Join already source-bound filaments without losing independent EA, EI,
    /// density, prestress, modal material/gas loss, initial motion or torsion.
    /// Each existing source's slope limit and the support's limit both apply.
    /// This is cold construction, not transfer from running fixed-end histories.
    ///
    /// The gas is explicitly CO-MOVING with the translating support: inherited
    /// cylinder drag damps relative filament motion, not rigid translation.
    /// All sources must share that coherent gas and their actual physical clock.
    /// A laboratory-stationary atmosphere requires a different absolute-flow
    /// resistance model and cannot be silently represented by this constructor.
    /// The support's independent mass, spring and loss are explicit input data.
    ///
    /// # Errors
    /// Empty/oversized basis, mixed clocks/gases, source-shape mismatch, invalid
    /// initial motion or support, or cancellation. No physical state advances.
    pub fn in_comoving_gas(sources: Vec<NonlinearMaterialString>, support: TranslatingSupport,
        gate: &CancelGate) -> Result<Self, ImpactError>
    {
        Self::build(sources, support, None, SupportGasFrame::CoMoving, None, gate)
    }

    /// Source-defined strings with selected moving ends in a stationary gas.
    /// `Second` fixes the nut at station zero and moves the bridge at station
    /// one; `First` reverses that assignment. A shared bridge can thus join a
    /// speaking span and a distinct afterlength without duplicating its mass.
    /// All spans share one bridge translation, along each source's declared
    /// bridge direction when it retains two planes. This is not a full
    /// rotating bridge or a carved body. The support input is its effective
    /// mass/stiffness/loss, independent from all string material definitions.
    ///
    /// Material curvature damping retains the original relative sine basis.
    /// Gas drag is separated using the original cylinder-resistance owner and
    /// applied to ABSOLUTE velocity projections, including the unresolved lift.
    /// No total damping is added twice. A numerically unresolved separation
    /// refuses rather than clipping negative resistance or dropping source loss.
    /// `remainder_frequency_hz` explicitly samples the gas law for the orthogonal
    /// lift tail; resolved modes retain their original evaluation frequencies.
    /// This is frozen spectral gas resistance, not a broadband fluid solution.
    ///
    /// # Errors
    /// Original source/support admission, wrong endpoint count, invalid tail
    /// frequency, loss separation, mixed gases/clocks or cancellation. Original
    /// material records, initial modes and torsion remain attached throughout.
    pub fn in_stationary_gas(sources: Vec<NonlinearMaterialString>, support: TranslatingSupport,
        motion: Vec<StringEndMotion>, remainder_frequency_hz: f64, gate: &CancelGate)
        -> Result<Self, ImpactError>
    {
        Self::build(sources, support, Some(&motion),
            SupportGasFrame::StationaryModal { remainder_frequency_hz }, None, gate)
    }

    fn build(sources: Vec<NonlinearMaterialString>, support: TranslatingSupport,
        motion: Option<&[StringEndMotion]>, gas_frame: SupportGasFrame,
        longitudinal:Option<&[MaterialLongitudinalSupport]>, gate: &CancelGate)
        -> Result<Self, ImpactError>
    {
        poll(gate)?;
        support.validate()?;
        if sources.is_empty() || sources.len() > 32 {
            return Err(invalid("material support needs 1..=32 original string sources"));
        }
        let sample_rate_hz = sources[0].sample_rate_hz();
        let gas = *sources[0].gas();
        if motion.is_some_and(|m| m.len()!=sources.len())
            || longitudinal.is_some_and(|l|l.len()!=sources.len()) {
            return Err(invalid("material support requires an endpoint declaration for each string"));
        }
        if let SupportGasFrame::StationaryModal { remainder_frequency_hz } = gas_frame {
            if !remainder_frequency_hz.is_finite() || remainder_frequency_hz<=0.0
                || remainder_frequency_hz>=0.5*f64::from(sample_rate_hz) {
                return Err(invalid("stationary string gas requires an explicit positive sub-Nyquist remainder frequency"));
            }
        }
        let mode_count = sources.iter().enumerate().try_fold(1usize, |n,(i,s)|
            n.checked_add(s.mode_count()).and_then(|n|n.checked_add(longitudinal.map_or(0,|a|a[i].resolution.axial_modes.len()))))
            .filter(|&n| n <= MAX_IMPACT_MODES)
            .ok_or_else(|| invalid("material strings plus support exceed the complete mode budget"))?;
        let mut wires = Vec::with_capacity(sources.len());
        let mut laws = Vec::with_capacity(sources.len());
        let mut relative = Vec::with_capacity(mode_count - 1);
        let mut drag = Vec::with_capacity(sources.len());
        let mut transverse = Vec::with_capacity(sources.len());
        let mut axial=Vec::new();let mut bands=Vec::new();let mut geometry=Vec::new();
        for (index,source) in sources.iter().enumerate() {
            poll(gate)?;
            if source.sample_rate_hz() != sample_rate_hz || *source.gas() != gas {
                return Err(invalid("supported material strings require the same clock and complete common gas"));
            }
            let [bending, tail @ ..] = source.bodies() else {
                return Err(invalid("supported source is missing its bending image"));
            };
            let BodyPotential::String(potential) = &bending.potential else {
                return Err(invalid("supported source must retain its original nonlinear string law"));
            };
            let motion=match source.source() {
                MaterialStringSource::WithPolarizations(p) => StringTransverseMotion::TwoPlanes {
                    support_direction: p.directions().bridge(),
                },
                _ => StringTransverseMotion::OnePlane,
            };
            let planes=if matches!(motion,StringTransverseMotion::TwoPlanes {..}) {2} else {1};
            if tail.len() > 1 || potential.transverse_planes()!=planes
                || bending.initial.len() != planes*source.card().mode_count {
                return Err(invalid("supported source has an incompatible bending/torsion layout"));
            }
            transverse.push(motion);
            match gas_frame {
                SupportGasFrame::CoMoving => wires.push(source.wire_span().clone()),
                SupportGasFrame::StationaryModal { remainder_frequency_hz } => {
                    let (wire, law)=stationary::split_loss(source, remainder_frequency_hz, gate)?;
                    wires.push(wire);drag.push(law);
                }
            }
            laws.push(Some(StringStretching { axial_rigidity_n: source.axial_rigidity_n(),
                maximum_slope: longitudinal.map_or(potential.slope_limit(),|a|
                    potential.slope_limit().min(a[index].resolution.maximum_slope)) }));
            relative.extend_from_slice(&bending.initial);
            if let Some(attachments)=longitudinal {
                let (selection,band)=attachments[index].lower(source)?;
                axial.push(selection);bands.push(band);geometry.push(attachments[index].geometry);
                relative.extend_from_slice(&attachments[index].initial_axial);
            }
        }
        let default_ends=vec![StringEndMotion::Both;sources.len()];
        let supported=if longitudinal.is_some() {
            SupportedStrings::with_longitudinal_motion(&wires,&laws,motion.unwrap_or(&default_ends),
                &geometry,&axial,support)?
        } else {SupportedStrings::with_transverse_motion(&wires,&laws,
            motion.unwrap_or(&default_ends),&transverse,support)?};
        for (index,band) in bands.iter().enumerate() {
            if supported.axial_angular_frequencies(index)?.iter().any(|&omega|
                band.is_some_and(|(lo,hi)|omega<lo || omega>hi)) {
                return Err(invalid("supported axial mode exceeds original material elasticity/density applicability"));
            }
        }
        let supported=if matches!(gas_frame,SupportGasFrame::StationaryModal {..}) {
            supported.with_absolute_drag(drag)?
        } else { supported };
        let mut bodies = vec![supported.body_with_relative_states(&relative)?];
        let mut torsion = Vec::with_capacity(sources.len());
        let mut at = supported.mode_count();
        for source in &sources {
            poll(gate)?;
            if let Some(twist) = source.bodies().get(1) {
                let end = at + twist.initial.len();
                torsion.push(Some(at..end)); bodies.push(twist.clone()); at = end;
            } else { torsion.push(None); }
        }
        if at != mode_count { return Err(invalid("supported material basis lost source coordinates")); }
        poll(gate)?;
        Ok(Self { sources, support, supported, bodies, torsion, mode_count, sample_rate_hz, gas_frame,
            longitudinal:longitudinal.map(|a|a.to_vec()),axial_material_loss:None,axial_material_relaxation:None })
    }

    /// Physical gas frame and the explicit unresolved-field approximation.
    pub const fn gas_frame(&self) -> SupportGasFrame { self.gas_frame }
    /// Original independent cards, complete material receipts and relative states.
    pub fn sources(&self) -> &[NonlinearMaterialString] { &self.sources }
    /// Original rail input, excluding wire mass that the shared owner adds once.
    pub const fn support(&self) -> &TranslatingSupport { &self.support }
    /// Existing support image, for its exact mass, contact and motion observations.
    pub const fn supported(&self) -> &SupportedStrings { &self.supported }
    /// Complete bodies for composition with other objects in the existing owner.
    /// Do not additionally insert the original fixed-end bending bodies.
    pub fn bodies(&self) -> &[ImpactBody] { &self.bodies }
    /// Complete shared-bending/support/torsion coordinate count.
    pub const fn mode_count(&self) -> usize { self.mode_count }
    /// Actual shared mechanical clock, not a proposed audio resampling rate.
    pub const fn sample_rate_hz(&self) -> u32 { self.sample_rate_hz }
    /// Original string order; absent torsion remains absent, not fabricated.
    pub fn torsion_ranges(&self) -> &[Option<Range<usize>>] { &self.torsion }

    /// Absolute physical transverse point map along the source's bow direction
    /// in the complete normalized basis. Only selected endpoints move;
    /// transverse centreline force applies no torque.
    pub fn transverse_shape(&self, string: usize, station: f64) -> Result<Vec<f64>, ImpactError> {
        let source=self.sources.get(string).ok_or_else(|| invalid("unknown supported material string"))?;
        let direction=match source.source() {
            MaterialStringSource::WithPolarizations(p) => p.directions().bow(),
            _ => [1.0,0.0],
        };
        self.transverse_shape_along(string,station,
            TransverseDirection::try_new(direction).map_err(|e|ImpactError::Owner(e.to_string()))?)
    }
    /// Physical force/velocity row in an independently declared transverse
    /// direction. Both planes and their common support reaction are retained;
    /// a normal centreline load carries no invented torsional moment.
    pub fn transverse_shape_along(&self, string: usize, station: f64,
        direction: TransverseDirection) -> Result<Vec<f64>, ImpactError> {
        let mut row = self.supported.absolute_shape_along(string, station, direction.components())?;
        row.resize(self.mode_count, 0.0);
        Ok(row)
    }
    /// Tangential surface map, including the source-derived r*theta contribution.
    /// Uses the SAME row for physical velocity and transposed force/torque. This
    /// declares an actuator port, not a friction coefficient or contact force.
    pub fn surface_shape(&self, string: usize, station: f64) -> Result<Vec<f64>, ImpactError> {
        let mut row = self.transverse_shape(string, station)?;
        if let Some(range) = &self.torsion[string] {
            // Fixed twist at either endpoint is exactly zero, not sin(pi) roundoff.
            if station != 0.0 && station != 1.0 {
                let source = &self.sources[string];
                let original = source.surface_shape(station)?;
                row[range.clone()].copy_from_slice(&original[source.transverse_mode_count()..]);
            }
        }
        Ok(row)
    }
    /// Physical support-force row, with all filament and torsion entries zero.
    pub fn support_shape(&self) -> Vec<f64> {
        let mut row = vec![0.0; self.mode_count];
        row[self.supported.support_coordinate()] = self.supported.force_weight(); row
    }

    /// Prepare the original nonlinear time owner while retaining every source.
    /// External generalized forces use this image's work-conjugate point maps.
    /// No excitation is invented; initial stored vibration is preserved. To
    /// include other bodies or contacts, compose `bodies()` with ImpactSystem.
    pub fn prepare(self, config: ImpactConfig, gate: &CancelGate) -> Result<MaterialSupportRun, ImpactError> {
        poll(gate)?;
        if config.dt_s.to_bits() != (1.0/f64::from(self.sample_rate_hz)).to_bits() {
            return Err(invalid("material support runtime must use its actual source clock"));
        }
        let system = ImpactSystem::new(self.bodies.clone(), vec![], vec![], vec![], config)?.prepare_analytic()?;
        poll(gate)?;
        Ok(MaterialSupportRun { source: self, system })
    }
}

/// Source-retaining runtime. All nonlinear trials, storage, dissipation, work
/// checks and rollback remain in the original prepared ImpactSystem.
pub struct MaterialSupportRun {
    source: SupportedMaterialStrings,
    system: PreparedImpactSystem,
}
impl MaterialSupportRun {
    /// Original materials, gas, support prescription, coordinate and force maps.
    pub const fn source(&self) -> &SupportedMaterialStrings { &self.source }
    /// Read accepted physical state, time, storage and support diagnostics.
    pub const fn system(&self) -> &PreparedImpactSystem { &self.system }
    /// One original transactional mechanical step; no extra clock or integrator.
    pub fn step(&mut self, external: &[f64], gate: &CancelGate) -> Result<ImpactFrame, ImpactError> {
        self.system.step(external, gate)
    }
}
