//! Material-derived nonlinear filaments on ONE inertially coupled support.
//! The same material section lowers EA, EI, mass and loss; SupportedStrings
//! owns the kinetic transformation and reactions; ImpactSystem owns every step.
//! No string is rebuilt from a pitch, averaged with a neighbor or muted.
//!
//! This image carries both endpoints on a common translating rail, not one
//! moving violin bridge and one fixed nut. Gas drag is the original relative
//! modal image in a CO-MOVING gas frame, explicitly selected by the constructor.
//! It does not model translation through stationary gas, fluid added mass,
//! support rotation, live material changes, winding slip or a bow friction law.
use std::ops::Range;
use super::{NonlinearMaterialString, BodyPotential, ImpactBody, ImpactError, MAX_IMPACT_MODES, invalid};
use crate::render::plate::impact::{ImpactConfig, ImpactFrame, ImpactSystem, PreparedImpactSystem};
use crate::render::plate::impact::string::StringStretching;
use crate::render::plate::impact::supported::{SupportedStrings, TranslatingSupport};
use fs_exec::CancelGate;

/// Immutable source geometry/materials and their complete coupled initial image.
/// Body zero contains every bending coordinate followed by the rail coordinate.
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
        poll(gate)?;
        support.validate()?;
        if sources.is_empty() || sources.len() > 32 {
            return Err(invalid("material support needs 1..=32 original string sources"));
        }
        let sample_rate_hz = sources[0].sample_rate_hz();
        let gas = *sources[0].gas();
        let mode_count = sources.iter().try_fold(1usize, |n,s| n.checked_add(s.mode_count()))
            .filter(|&n| n <= MAX_IMPACT_MODES)
            .ok_or_else(|| invalid("material strings plus support exceed the complete mode budget"))?;
        let mut wires = Vec::with_capacity(sources.len());
        let mut laws = Vec::with_capacity(sources.len());
        let mut relative = Vec::with_capacity(mode_count - 1);
        for source in &sources {
            poll(gate)?;
            if source.sample_rate_hz() != sample_rate_hz || *source.gas() != gas {
                return Err(invalid("supported material strings require the same clock and co-moving gas"));
            }
            let [bending, tail @ ..] = source.bodies() else {
                return Err(invalid("supported source is missing its bending image"));
            };
            let BodyPotential::String(potential) = &bending.potential else {
                return Err(invalid("supported source must retain its original nonlinear string law"));
            };
            if tail.len() > 1 || bending.initial.len() != source.card().mode_count {
                return Err(invalid("supported source has an incompatible bending/torsion layout"));
            }
            wires.push(source.wire_span().clone());
            laws.push(Some(StringStretching { axial_rigidity_n: source.axial_rigidity_n(),
                maximum_slope: potential.slope_limit() }));
            relative.extend_from_slice(&bending.initial);
        }
        let supported = SupportedStrings::with_individual_stretching(&wires, &laws, support)?;
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
        Ok(Self { sources, support, supported, bodies, torsion, mode_count, sample_rate_hz })
    }

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

    /// Absolute physical transverse point map in the complete normalized basis.
    /// Both endpoints move with the support; transverse force applies no torque.
    pub fn transverse_shape(&self, string: usize, station: f64) -> Result<Vec<f64>, ImpactError> {
        let mut row = self.supported.absolute_shape(string, station)?;
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
                row[range.clone()].copy_from_slice(&original[source.card().mode_count..]);
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
