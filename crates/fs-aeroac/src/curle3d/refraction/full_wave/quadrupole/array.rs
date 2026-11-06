//! Original volume-stress cells through one refracting fluid interface.
//!
//! Every cell keeps its own six physical stress moments, position and phase.
//! The native single-source owner prepares each transfer; observations sum
//! complex fields before forming intensity or statistically specified moments.

use super::super::reflection::{PreparedReflectedFullWaveQuadrupole, field, line, singular_total};
use super::super::{PlanarFluidInterface, QuadratureLimits, poll, positive};
use super::{FullWaveQuadrupoleDiagnostics, FullWaveQuadrupoleLimits, PreparedFullWaveQuadrupole};
use crate::AeroacError;
use crate::curle3d::quadrupole::QuadrupolePhasor;
use crate::curle3d::statistics::{ComplexSum, StatisticsSum, finite, invalid, nonfinite};
use crate::curle3d::{CoherenceBudget, CoherentCompactField, CompactFieldStatistics};
use crate::observation::HarmonicFieldLine;
use fs_math::c64::C64;

mod recorded;

/// The physical source and propagation assumptions are shared by all views.
pub const FULL_WAVE_QUADRUPOLE_ARRAY_SCOPE: &str = "Prescribed stationary symmetric \
Lighthill stress moments [N m] in their original source positions and clock; coherent \
full-wave transmission across one infinite planar interface between stationary \
homogeneous lossless fluids. Full tensor pressure and vector velocity retain near, \
critical and evanescent fields. Source statistics are explicitly supplied; per-source \
quadrature estimates and spectral-tail bounds give no relative guarantee at array \
nulls. No source-size or SGS validation, absorption, mean flow, finite boundaries, \
source loading, interval certificate or measured-SPL qualification.";

/// Boundary returns retain every cell's physical tensor and mutual coherence.
pub const FULL_WAVE_REFLECTED_QUADRUPOLE_ARRAY_SCOPE: &str = "Prescribed stationary \
symmetric Lighthill stress moments [N m] at their original positions and clock; \
coherent reflected-only pressure and vector velocity in the original incident \
fluid at one infinite planar lossless fluid interface. Finite reflected self-fields \
are admitted. Source statistics are explicitly supplied; per-source quadrature \
estimates and analytic tails give no relative guarantee at array nulls. No direct \
self-field, source loading, source-size or SGS validation, absorption, mean flow, \
finite boundaries, interval certificate or measured-SPL qualification.";

/// Direct and returning fields interfere before any energy or statistics query.
pub const FULL_WAVE_TOTAL_QUADRUPOLE_ARRAY_SCOPE: &str = "Prescribed stationary \
symmetric Lighthill stress moments [N m] at their original positions and clock; \
coherent direct-plus-reflected pressure and vector velocity in the original incident \
fluid at one infinite planar lossless fluid interface. Source statistics are \
explicitly supplied. Numerical relative control applies to reflected transfers; \
their absolute estimates carry to total fields, without a relative guarantee at \
interference nulls. No singular total self-field, source loading, source-size or \
SGS validation, absorption, mean flow, finite boundaries, interval certificate or \
measured-SPL qualification.";

/// Physical field cached at the original receiver, including its direct-path policy.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FullWaveQuadrupoleField {
    Transmitted,
    Reflected,
    IncidentTotal,
}

/// One source-count and quadrature budget for the complete physical array.
#[derive(Clone, Copy, Debug)]
pub struct FullWaveQuadrupoleSourceLimits {
    pub max_sources: usize,
    /// Positive aggregate cap, including discarded adaptive quadrature panels.
    /// Exact equal-speed transfers require no spectral evaluations.
    pub maximum_total_evaluations: usize,
    /// Accuracy and work limits of each original unit stress transfer.
    pub quadrupole: FullWaveQuadrupoleLimits,
}

#[derive(Clone, Debug)]
struct Transfer {
    source_m: [f64; 3],
    field: [[C64; 6]; 4],
    diagnostics: FullWaveQuadrupoleDiagnostics,
}

/// Prepared coherent observation of original stress cells across a fluid plane.
/// Preparation owns bounded allocation and quadrature. Warm coherent and
/// low-rank statistical observations allocate nothing and use fixed-order DD
/// sums; neither a new propagation delay nor a plane-wave impedance is applied.
#[derive(Clone, Debug)]
pub struct PreparedFullWaveQuadrupoleObserver {
    interface: PlanarFluidInterface,
    frequency_hz: f64,
    observer_m: [f64; 3],
    field_kind: FullWaveQuadrupoleField,
    transfers: Vec<Transfer>,
    evaluations: usize,
}

impl PreparedFullWaveQuadrupoleObserver {
    /// Prepare each source strictly in the incident fluid and one receiver
    /// strictly in the receiving fluid. Empty arrays give silence only after
    /// admitting the original observer, frequency and numerical controls.
    /// Per-source errors multiply the physical stress Frobenius norm; they do
    /// not imply relative accuracy after destructive inter-cell interference.
    pub fn new(
        interface: PlanarFluidInterface,
        frequency_hz: f64,
        observer_m: [f64; 3],
        source_positions_m: &[[f64; 3]],
        limits: FullWaveQuadrupoleSourceLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        Self::prepare(
            interface,
            frequency_hz,
            observer_m,
            source_positions_m,
            limits,
            FullWaveQuadrupoleField::Transmitted,
            cancelled,
        )
    }

    /// Observe only the boundary returns in the original incident fluid.
    /// A cell may coincide with the receiver: its reflected self-field is finite.
    pub fn new_reflected(
        interface: PlanarFluidInterface,
        frequency_hz: f64,
        observer_m: [f64; 3],
        source_positions_m: &[[f64; 3]],
        limits: FullWaveQuadrupoleSourceLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        Self::prepare(
            interface,
            frequency_hz,
            observer_m,
            source_positions_m,
            limits,
            FullWaveQuadrupoleField::Reflected,
            cancelled,
        )
    }

    /// Observe coherent direct-plus-reflected fields in the original incident
    /// fluid. Every source must be distinct from the receiver. The original
    /// single-source owner adds both paths before this array's coherent or
    /// statistical reduction; their interference is never discarded.
    pub fn new_total(
        interface: PlanarFluidInterface,
        frequency_hz: f64,
        observer_m: [f64; 3],
        source_positions_m: &[[f64; 3]],
        limits: FullWaveQuadrupoleSourceLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        Self::prepare(
            interface,
            frequency_hz,
            observer_m,
            source_positions_m,
            limits,
            FullWaveQuadrupoleField::IncidentTotal,
            cancelled,
        )
    }

    fn prepare(
        interface: PlanarFluidInterface,
        frequency_hz: f64,
        observer_m: [f64; 3],
        source_positions_m: &[[f64; 3]],
        limits: FullWaveQuadrupoleSourceLimits,
        field_kind: FullWaveQuadrupoleField,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        poll(cancelled)?;
        QuadratureLimits::from(limits.quadrupole).admit()?;
        if source_positions_m.len() > limits.max_sources || limits.maximum_total_evaluations == 0 {
            return Err(invalid("full-wave stress array source or work budget"));
        }
        let distance = interface.distance(observer_m).to_f64();
        let clearance = match field_kind {
            FullWaveQuadrupoleField::Transmitted => distance,
            FullWaveQuadrupoleField::Reflected | FullWaveQuadrupoleField::IncidentTotal => {
                -distance
            }
        };
        if !positive(core::f64::consts::TAU * frequency_hz)
            || observer_m.iter().any(|v| !v.is_finite())
            || !positive(clearance)
        {
            return Err(invalid(
                "full-wave stress array frequency or observer fluid side",
            ));
        }
        let needs_quadrature =
            interface.incident.sound_speed_m_s() != interface.receiving.sound_speed_m_s();
        let mut transfers = Vec::new();
        transfers
            .try_reserve_exact(source_positions_m.len())
            .map_err(|_| invalid("full-wave stress array allocation"))?;
        let mut evaluations = 0usize;
        for &source_m in source_positions_m {
            poll(cancelled)?;
            if field_kind == FullWaveQuadrupoleField::IncidentTotal && source_m == observer_m {
                return Err(singular_total());
            }
            let remaining = limits.maximum_total_evaluations - evaluations;
            if needs_quadrature && remaining == 0 {
                return Err(invalid(
                    "full-wave stress array exhausted its aggregate work budget",
                ));
            }
            let mut source_limits = limits.quadrupole;
            if needs_quadrature {
                source_limits.maximum_evaluations =
                    source_limits.maximum_evaluations.min(remaining);
            }
            let (field, diagnostics) = match field_kind {
                FullWaveQuadrupoleField::Transmitted => {
                    let prepared = PreparedFullWaveQuadrupole::new(
                        interface,
                        frequency_hz,
                        source_m,
                        observer_m,
                        source_limits,
                        cancelled,
                    )?;
                    (prepared.transfer, *prepared.diagnostics())
                }
                FullWaveQuadrupoleField::Reflected | FullWaveQuadrupoleField::IncidentTotal => {
                    let prepared = PreparedReflectedFullWaveQuadrupole::new(
                        interface,
                        frequency_hz,
                        source_m,
                        observer_m,
                        source_limits,
                        cancelled,
                    )?;
                    // Extract the native cached map through six physical unit
                    // moments. This performs no extra quadrature or derivatives
                    // and keeps all direct/image tensor conventions with its owner.
                    let mut field = [[C64::ZERO; 6]; 4];
                    for column in 0..6 {
                        let mut unit = QuadrupolePhasor::default();
                        unit.stress_moment_n_m[column] = C64::ONE;
                        let line = if field_kind == FullWaveQuadrupoleField::Reflected {
                            prepared.evaluate_reflected_harmonic_line(unit)?
                        } else {
                            prepared.evaluate_total_harmonic_line(unit)?
                        };
                        field[0][column] = line.pressure_pa;
                        for axis in 0..3 {
                            field[axis + 1][column] = line.particle_velocity_m_s[axis];
                        }
                    }
                    (field, *prepared.diagnostics())
                }
            };
            evaluations = evaluations
                .checked_add(diagnostics.evaluations)
                .filter(|n| *n <= limits.maximum_total_evaluations)
                .ok_or_else(|| {
                    invalid("full-wave stress array exceeded its aggregate work budget")
                })?;
            transfers.push(Transfer {
                source_m,
                field,
                diagnostics,
            });
        }
        poll(cancelled)?;
        Ok(Self {
            interface,
            frequency_hz,
            observer_m,
            field_kind,
            transfers,
            evaluations,
        })
    }

    pub const fn interface(&self) -> PlanarFluidInterface {
        self.interface
    }
    pub const fn frequency_hz(&self) -> f64 {
        self.frequency_hz
    }
    pub const fn observer_position_m(&self) -> [f64; 3] {
        self.observer_m
    }
    pub const fn field_kind(&self) -> FullWaveQuadrupoleField {
        self.field_kind
    }
    pub fn source_count(&self) -> usize {
        self.transfers.len()
    }
    pub const fn evaluations(&self) -> usize {
        self.evaluations
    }
    pub const fn scope(&self) -> &'static str {
        match self.field_kind {
            FullWaveQuadrupoleField::Transmitted => FULL_WAVE_QUADRUPOLE_ARRAY_SCOPE,
            FullWaveQuadrupoleField::Reflected => FULL_WAVE_REFLECTED_QUADRUPOLE_ARRAY_SCOPE,
            FullWaveQuadrupoleField::IncidentTotal => FULL_WAVE_TOTAL_QUADRUPOLE_ARRAY_SCOPE,
        }
    }
    pub fn source_position_m(&self, index: usize) -> Option<[f64; 3]> {
        self.transfers.get(index).map(|source| source.source_m)
    }
    /// Transmission diagnostics for transmitted arrays; reflected-transfer
    /// diagnostics for both incident-side views. The native direct contribution
    /// adds no quadrature error, but a total-field null has no relative guarantee.
    pub fn source_diagnostics(&self, index: usize) -> Option<&FullWaveQuadrupoleDiagnostics> {
        self.transfers.get(index).map(|source| &source.diagnostics)
    }

    /// Sum all tensor components and cells coherently before computing flux.
    pub fn evaluate(
        &self,
        sources: &[QuadrupolePhasor],
    ) -> Result<CoherentCompactField, AeroacError> {
        self.evaluate_cancellable(sources, || false)
    }

    pub fn evaluate_cancellable(
        &self,
        sources: &[QuadrupolePhasor],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CoherentCompactField, AeroacError> {
        field(self.amplitudes(sources, &mut cancelled)?, self.scope())
    }

    /// The native line retains the original source clock and full propagation phase.
    pub fn evaluate_harmonic_line(
        &self,
        sources: &[QuadrupolePhasor],
    ) -> Result<HarmonicFieldLine<3>, AeroacError> {
        Ok(line(
            self.frequency_hz,
            self.amplitudes(sources, &mut || false)?,
        ))
    }

    /// Explicitly independent cells, with coherent tensor entries inside each.
    /// Independence is a caller assertion, never inferred from cell separation.
    pub fn evaluate_independent(
        &self,
        sources: &[QuadrupolePhasor],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CompactFieldStatistics, AeroacError> {
        poll(&mut cancelled)?;
        self.admit_count(sources.len())?;
        let mut statistics = StatisticsSum::default();
        for (i, (transfer, source)) in self.transfers.iter().zip(sources).enumerate() {
            if i.is_multiple_of(256) {
                poll(&mut cancelled)?;
            }
            let mut sums = [ComplexSum::default(); 4];
            add_source(transfer, source, &mut sums)?;
            statistics.add(sums.map(ComplexSum::value))?;
        }
        poll(&mut cancelled)?;
        statistics.finish(self.scope())
    }

    /// Low-rank source coherence factors expressed as peak stress moments.
    /// Each mode retains all cell/tensor interference before quadratic reduction.
    pub fn evaluate_coherence_modes(
        &self,
        modes: &[&[QuadrupolePhasor]],
        budget: CoherenceBudget,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CompactFieldStatistics, AeroacError> {
        poll(&mut cancelled)?;
        self.admit_modes(modes.len(), budget)?;
        if modes.iter().any(|mode| mode.len() != self.source_count()) {
            return Err(invalid("full-wave stress coherence factor source count"));
        }
        let mut statistics = StatisticsSum::default();
        for mode in modes {
            statistics.add(self.amplitudes(mode, &mut cancelled)?)?;
        }
        poll(&mut cancelled)?;
        statistics.finish(self.scope())
    }

    fn admit_count(&self, count: usize) -> Result<(), AeroacError> {
        if count != self.source_count() {
            return Err(invalid(
                "one physical stress moment is required per prepared cell",
            ));
        }
        Ok(())
    }

    fn admit_modes(&self, modes: usize, budget: CoherenceBudget) -> Result<(), AeroacError> {
        if modes > budget.max_modes
            || self
                .source_count()
                .checked_mul(modes)
                .is_none_or(|work| work > budget.max_source_mode_pairs)
        {
            return Err(invalid("full-wave stress coherence exceeds work budget"));
        }
        Ok(())
    }

    fn amplitudes(
        &self,
        sources: &[QuadrupolePhasor],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<[C64; 4], AeroacError> {
        poll(cancelled)?;
        self.admit_count(sources.len())?;
        let mut sums = [ComplexSum::default(); 4];
        for (i, (transfer, source)) in self.transfers.iter().zip(sources).enumerate() {
            if i.is_multiple_of(256) {
                poll(cancelled)?;
            }
            add_source(transfer, source, &mut sums)?;
        }
        poll(cancelled)?;
        let values = sums.map(ComplexSum::value);
        if values.iter().any(|v| !finite(*v)) {
            return Err(nonfinite("accumulated full-wave stress field"));
        }
        Ok(values)
    }
}

fn add_source(
    transfer: &Transfer,
    source: &QuadrupolePhasor,
    sums: &mut [ComplexSum; 4],
) -> Result<(), AeroacError> {
    for (entry, value) in source.stress_moment_n_m.into_iter().enumerate() {
        if !finite(value) {
            return Err(nonfinite("full-wave stress moment amplitude"));
        }
        if value != C64::ZERO {
            for (row, sum) in transfer.field.iter().zip(&mut *sums) {
                sum.add(row[entry] * value)?;
            }
        }
    }
    Ok(())
}
