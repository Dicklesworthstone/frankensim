//! Original volume-stress coherence at a passive impedance wall.
//!
//! Native per-cell 4-by-6 maps feed the shared stress-array reduction. Tensor
//! entries and direct/reflected paths remain coherent inside every physical
//! source, including when the caller explicitly declares independent cells.

use super::super::super::array::remaining_budget;
use super::super::super::quadrupole::array::{
    ArrayField, FullWaveQuadrupoleSourceLimits, Transfer,
};
use super::super::super::reflection::singular_total;
use super::super::super::{QuadratureLimits, poll, positive};
use super::super::{ImpedanceWallFieldKind, PlanarImpedanceWall};
use super::{ImpedanceWallQuadrupoleDiagnostics, PreparedImpedanceWallQuadrupole};
use crate::AeroacError;
use crate::curle3d::quadrupole::QuadrupolePhasor;
use crate::curle3d::statistics::invalid;
use crate::curle3d::{CoherenceBudget, CoherentCompactField, CompactFieldStatistics};
use crate::observation::HarmonicFieldLine;
use fs_math::c64::C64;

/// Every returning cell retains its original tensor, position and source clock.
pub const IMPEDANCE_WALL_REFLECTED_QUADRUPOLE_ARRAY_SCOPE: &str = "Prescribed stationary \
symmetric 3D Lighthill stress moments [N m] at their original positions and clock; \
reflected-only spherical pressure and vector velocity at an infinite passive locally \
reacting impedance wall in a homogeneous lossless fluid. Full trace, physical shear, \
near and evanescent fields; explicitly supplied source statistics. Per-cell quadrature \
estimates and analytic tails give no relative guarantee at array nulls. No direct \
field, lossless compliant pole, finite-source or SGS validation, finite edges, nonlocal \
wall response, mean flow, bulk absorption, source feedback, interval certificate or \
measured-SPL qualification.";

/// Direct and reflected paths interfere within each original cell before statistics.
pub const IMPEDANCE_WALL_TOTAL_QUADRUPOLE_ARRAY_SCOPE: &str = "Prescribed stationary \
symmetric 3D Lighthill stress moments [N m] at their original positions and clock; \
coherent direct-plus-reflected spherical pressure and vector velocity at an infinite \
passive locally reacting impedance wall in a homogeneous lossless fluid. Explicitly \
supplied source statistics retain all tensor and per-cell direct/reflected cross terms. \
Reflected-transfer absolute estimates give no relative guarantee at interference nulls. \
No singular total self-field, lossless compliant pole, finite-source or SGS validation, \
finite edges, nonlocal wall response, mean flow, bulk absorption, source feedback, \
interval certificate or measured-SPL qualification.";

/// Prepared observation of physical volume-stress cells beside a passive wall.
/// Preparation owns bounded allocation and quadrature. Warm application reuses
/// the shared compensated stress reduction without allocation, quadrature,
/// transcendental evaluation or another propagation delay.
#[derive(Clone, Debug)]
pub struct PreparedImpedanceWallQuadrupoleObserver {
    wall: PlanarImpedanceWall,
    field_kind: ImpedanceWallFieldKind,
    frequency_hz: f64,
    observer_m: [f64; 3],
    transfers: Vec<Transfer<ImpedanceWallQuadrupoleDiagnostics>>,
    evaluations: usize,
}

impl PreparedImpedanceWallQuadrupoleObserver {
    /// Sources lie strictly in the fluid; the receiver may lie on the wall.
    /// A coincident cell retains its finite reflected self-field, without
    /// claiming source loading. Empty arrays are validated silence. One shared
    /// work cap includes every cell; ideal rigid/release images consume none.
    pub fn new_reflected(
        wall: PlanarImpedanceWall,
        frequency_hz: f64,
        observer_m: [f64; 3],
        source_positions_m: &[[f64; 3]],
        limits: FullWaveQuadrupoleSourceLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        Self::prepare(
            wall,
            ImpedanceWallFieldKind::Reflected,
            frequency_hz,
            observer_m,
            source_positions_m,
            limits,
            cancelled,
        )
    }

    /// Each cell's direct and reflected fields interfere before any coherent
    /// or statistical reduction. A source at the receiver refuses its singular
    /// direct field regardless of later amplitude. Cancellation or exhaustion
    /// publishes no partially prepared array.
    pub fn new_total(
        wall: PlanarImpedanceWall,
        frequency_hz: f64,
        observer_m: [f64; 3],
        source_positions_m: &[[f64; 3]],
        limits: FullWaveQuadrupoleSourceLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        Self::prepare(
            wall,
            ImpedanceWallFieldKind::Total,
            frequency_hz,
            observer_m,
            source_positions_m,
            limits,
            cancelled,
        )
    }

    fn prepare(
        wall: PlanarImpedanceWall,
        field_kind: ImpedanceWallFieldKind,
        frequency_hz: f64,
        observer_m: [f64; 3],
        source_positions_m: &[[f64; 3]],
        limits: FullWaveQuadrupoleSourceLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        poll(cancelled)?;
        QuadratureLimits::from(limits.quadrupole).admit()?;
        if source_positions_m.len() > limits.max_sources || limits.maximum_total_evaluations == 0 {
            return Err(invalid("impedance wall stress array source or work budget"));
        }
        let distance = wall.geometry.distance(observer_m).to_f64();
        if !positive(core::f64::consts::TAU * frequency_hz)
            || observer_m.iter().any(|v| !v.is_finite())
            || !distance.is_finite()
            || distance > 0.0
        {
            return Err(invalid(
                "impedance wall stress array frequency or observer half-space",
            ));
        }
        let exact = wall.normalized.is_none() || wall.normalized == Some(C64::ZERO);
        let mut transfers = Vec::new();
        transfers
            .try_reserve_exact(source_positions_m.len())
            .map_err(|_| invalid("impedance wall stress array allocation"))?;
        let mut evaluations = 0usize;
        for &source_m in source_positions_m {
            poll(cancelled)?;
            if field_kind == ImpedanceWallFieldKind::Total && source_m == observer_m {
                return Err(singular_total());
            }
            let mut source_limits = limits.quadrupole;
            source_limits.maximum_evaluations = remaining_budget(
                source_limits.maximum_evaluations,
                limits.maximum_total_evaluations - evaluations,
                exact,
            )?;
            let prepared = PreparedImpedanceWallQuadrupole::new(
                wall,
                frequency_hz,
                source_m,
                observer_m,
                source_limits,
                cancelled,
            )?;
            // Reuse the native physical map, including both shear contractions
            // and coherent direct/returned fields; no extra source derivatives.
            let field = match field_kind {
                ImpedanceWallFieldKind::Reflected => prepared.reflected,
                ImpedanceWallFieldKind::Total => prepared.total.ok_or_else(singular_total)?,
            };
            let diagnostics = *prepared.diagnostics();
            evaluations = evaluations
                .checked_add(diagnostics.evaluations)
                .filter(|count| *count <= limits.maximum_total_evaluations)
                .ok_or_else(|| invalid("impedance wall stress array aggregate work budget"))?;
            transfers.push(Transfer {
                source_m,
                field,
                diagnostics,
            });
        }
        poll(cancelled)?;
        Ok(Self {
            wall,
            field_kind,
            frequency_hz,
            observer_m,
            transfers,
            evaluations,
        })
    }

    pub const fn wall(&self) -> PlanarImpedanceWall {
        self.wall
    }
    pub const fn field_kind(&self) -> ImpedanceWallFieldKind {
        self.field_kind
    }
    pub const fn frequency_hz(&self) -> f64 {
        self.frequency_hz
    }
    pub const fn observer_position_m(&self) -> [f64; 3] {
        self.observer_m
    }
    pub const fn evaluations(&self) -> usize {
        self.evaluations
    }
    pub fn source_count(&self) -> usize {
        self.transfers.len()
    }
    pub fn source_position_m(&self, index: usize) -> Option<[f64; 3]> {
        self.transfers.get(index).map(|source| source.source_m)
    }
    /// Native reflected-transfer estimates also apply to total fields. Multiply
    /// each by the physical cell's Frobenius stress norm, counting shear twice;
    /// array cancellation does not supply a relative error guarantee.
    pub fn source_diagnostics(&self, index: usize) -> Option<&ImpedanceWallQuadrupoleDiagnostics> {
        self.transfers.get(index).map(|source| &source.diagnostics)
    }
    pub const fn scope(&self) -> &'static str {
        match self.field_kind {
            ImpedanceWallFieldKind::Reflected => IMPEDANCE_WALL_REFLECTED_QUADRUPOLE_ARRAY_SCOPE,
            ImpedanceWallFieldKind::Total => IMPEDANCE_WALL_TOTAL_QUADRUPOLE_ARRAY_SCOPE,
        }
    }

    /// Sum all physical tensor entries and cells before computing signed flux.
    pub fn evaluate(
        &self,
        sources: &[QuadrupolePhasor],
    ) -> Result<CoherentCompactField, AeroacError> {
        self.evaluate_cancellable(sources, || false)
    }
    pub fn evaluate_cancellable(
        &self,
        sources: &[QuadrupolePhasor],
        cancelled: impl FnMut() -> bool,
    ) -> Result<CoherentCompactField, AeroacError> {
        self.view().evaluate_cancellable(sources, cancelled)
    }
    /// The native line retains the original source clock and propagation phase.
    pub fn evaluate_harmonic_line(
        &self,
        sources: &[QuadrupolePhasor],
    ) -> Result<HarmonicFieldLine<3>, AeroacError> {
        self.view().evaluate_harmonic_line(sources)
    }
    /// Explicitly independent cells retain coherent tensor entries and both
    /// paths inside each cell. Independence is never inferred from separation.
    pub fn evaluate_independent(
        &self,
        sources: &[QuadrupolePhasor],
        cancelled: impl FnMut() -> bool,
    ) -> Result<CompactFieldStatistics, AeroacError> {
        self.view().evaluate_independent(sources, cancelled)
    }
    /// Peak-unit low-rank factors preserve all cell/tensor phases within each
    /// factor. Quadratic moments add only between independent factors; the
    /// result is a spectral-line second moment, not a power density per hertz.
    pub fn evaluate_coherence_modes(
        &self,
        modes: &[&[QuadrupolePhasor]],
        budget: CoherenceBudget,
        cancelled: impl FnMut() -> bool,
    ) -> Result<CompactFieldStatistics, AeroacError> {
        self.view()
            .evaluate_coherence_modes(modes, budget, cancelled)
    }

    fn view(&self) -> ArrayField<'_, ImpedanceWallQuadrupoleDiagnostics> {
        ArrayField {
            transfers: &self.transfers,
            frequency_hz: self.frequency_hz,
            scope: self.scope(),
        }
    }
}
