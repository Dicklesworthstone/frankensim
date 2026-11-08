//! Original compact-source coherence at a locally reacting impedance wall.
//!
//! Preparation folds the native Q/F wall kernels into one 4-by-4 Cartesian
//! transfer per physical source. The shared array application forms pressure
//! and full vector velocity before any quadratic observations, preserving
//! direct/reflected, Q/F and cross-source interference.

use super::super::array::{
    ArrayField, FullWaveCompactSourceLimits, Transfer, remaining_budget, set_column,
};
use super::super::{admit, force, poll, positive};
use super::{
    ImpedanceWallDiagnostics, ImpedanceWallForceDiagnostics, PlanarImpedanceWall,
    PreparedImpedanceWallForce, PreparedImpedanceWallMonopole,
};
use crate::AeroacError;
use crate::curle3d::statistics::invalid;
use crate::curle3d::{
    CoherenceBudget, CoherentCompactField, CompactFieldStatistics, CompactSourcePhasor,
};
use crate::observation::HarmonicFieldLine;
use fs_math::c64::C64;

mod recorded;

/// Returning fields retain the original source identity and joint Q/F phase.
pub const IMPEDANCE_WALL_REFLECTED_ARRAY_SCOPE: &str = "Prescribed stationary compact \
3D volume flows [m³/s] and forces on fluid [N]; reflected-only spherical fields from \
one infinite passive locally reacting impedance plane in a homogeneous lossless fluid. \
Propagating and evanescent components, original source phases and explicitly supplied \
source statistics. Per-source quadrature estimates and tail bounds are not interval \
certificates or relative guarantees at array nulls. No direct field, lossless compliant \
pole, finite edges, nonlocal wall response, mean flow, bulk absorption, source feedback \
or measured-SPL qualification.";

/// Each original source retains its own direct/reflected interference in statistics.
pub const IMPEDANCE_WALL_TOTAL_ARRAY_SCOPE: &str = "Prescribed stationary compact 3D \
volume flows [m³/s] and forces on fluid [N]; coherent direct-plus-reflected spherical \
fields at one infinite passive locally reacting impedance plane in a homogeneous \
lossless fluid. Original Q/F and cross-source phases and explicitly supplied statistics. \
Reflected-transfer error estimates give no relative guarantee at interference nulls. \
No singular direct self-field, lossless compliant pole, finite edges, nonlocal wall \
response, mean flow, bulk absorption, source feedback, interval certificate or \
measured-SPL qualification.";

/// Which physical field is present in every cached source transfer.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ImpedanceWallFieldKind {
    Reflected,
    Total,
}

/// Native unit-Q and unit-force evidence, with their distinct physical units.
/// Multiply Q errors by |Q| and force norm errors by ||F||, then sum over
/// original sources for an absolute field estimate. These estimates concern
/// the reflected transfer even when the prepared field contains the direct
/// contribution; unknown physical amplitudes or source coherence are not inferred.
#[derive(Clone, Copy, Debug)]
pub struct ImpedanceWallCompactSourceDiagnostics {
    pub monopole: ImpedanceWallDiagnostics,
    pub force: ImpedanceWallForceDiagnostics,
}

/// Prepared observation of a compact Q/F source assembly beside a passive wall.
///
/// Every source contributes one cached 4-by-4 transfer. Application is O(sources)
/// or O(sources × coherence factors), with compensated reductions and no
/// allocation, quadrature, transcendental evaluation or extra propagation delay.
#[derive(Clone, Debug)]
pub struct PreparedImpedanceWallCompactSourceObserver {
    wall: PlanarImpedanceWall,
    field_kind: ImpedanceWallFieldKind,
    frequency_hz: f64,
    observer_m: [f64; 3],
    transfers: Vec<Transfer<ImpedanceWallCompactSourceDiagnostics>>,
    evaluations: usize,
}

impl PreparedImpedanceWallCompactSourceObserver {
    /// Prepare the returning field with sources strictly in the fluid and the
    /// receiver in the fluid or on the wall. Coincident source/receiver points
    /// retain their finite reflected response, without claiming source loading.
    /// The shared cap includes every Q/F integrand evaluation across sources;
    /// ideal rigid/release images consume none. Empty arrays are validated silence.
    pub fn new_reflected(
        wall: PlanarImpedanceWall,
        frequency_hz: f64,
        observer_m: [f64; 3],
        source_positions_m: &[[f64; 3]],
        limits: FullWaveCompactSourceLimits,
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

    /// Prepare coherent direct-plus-reflected source transfers. Statistical
    /// independence applies only between caller-declared physical sources: each
    /// source's Q/F and direct/reflected cross terms survive. A coincident
    /// source refuses its singular direct field regardless of later amplitude.
    /// Failure or cancellation publishes no partial array.
    pub fn new_total(
        wall: PlanarImpedanceWall,
        frequency_hz: f64,
        observer_m: [f64; 3],
        source_positions_m: &[[f64; 3]],
        limits: FullWaveCompactSourceLimits,
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
        limits: FullWaveCompactSourceLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        poll(cancelled)?;
        admit(limits.monopole)?;
        force::admit(limits.force)?;
        if source_positions_m.len() > limits.max_sources || limits.maximum_total_evaluations == 0 {
            return Err(invalid(
                "impedance wall array exceeds source or evaluation budget",
            ));
        }
        let distance = wall.geometry.distance(observer_m).to_f64();
        if !positive(core::f64::consts::TAU * frequency_hz)
            || observer_m.iter().any(|v| !v.is_finite())
            || !distance.is_finite()
            || distance > 0.0
        {
            return Err(invalid(
                "impedance wall array frequency or observer half-space",
            ));
        }
        let exact = wall.normalized.is_none() || wall.normalized == Some(C64::ZERO);
        let mut transfers = Vec::new();
        transfers
            .try_reserve_exact(source_positions_m.len())
            .map_err(|_| invalid("impedance wall array allocation exceeds available memory"))?;
        let mut evaluations = 0;
        for &source_m in source_positions_m {
            poll(cancelled)?;
            if field_kind == ImpedanceWallFieldKind::Total && source_m == observer_m {
                return Err(invalid(
                    "impedance wall array has a singular direct self-field",
                ));
            }
            let mut monopole_limits = limits.monopole;
            monopole_limits.maximum_evaluations = remaining_budget(
                monopole_limits.maximum_evaluations,
                limits.maximum_total_evaluations - evaluations,
                exact,
            )?;
            let monopole = PreparedImpedanceWallMonopole::new(
                wall,
                frequency_hz,
                source_m,
                observer_m,
                monopole_limits,
                cancelled,
            )?;
            let mut field = [[C64::ZERO; 4]; 4];
            let line = match field_kind {
                ImpedanceWallFieldKind::Reflected => {
                    monopole.evaluate_reflected_harmonic_line(C64::ONE)?
                }
                ImpedanceWallFieldKind::Total => monopole.evaluate_total_harmonic_line(C64::ONE)?,
            };
            set_column(&mut field, 0, line);
            let monopole = *monopole.diagnostics();
            evaluations += monopole.evaluations;
            let mut force_limits = limits.force;
            force_limits.maximum_evaluations = remaining_budget(
                force_limits.maximum_evaluations,
                limits.maximum_total_evaluations - evaluations,
                exact,
            )?;
            let force = PreparedImpedanceWallForce::new(
                wall,
                frequency_hz,
                source_m,
                observer_m,
                force_limits,
                cancelled,
            )?;
            for column in 0..3 {
                let unit = core::array::from_fn(|i| if i == column { C64::ONE } else { C64::ZERO });
                let line = match field_kind {
                    ImpedanceWallFieldKind::Reflected => {
                        force.evaluate_reflected_harmonic_line(unit)?
                    }
                    ImpedanceWallFieldKind::Total => force.evaluate_total_harmonic_line(unit)?,
                };
                set_column(&mut field, column + 1, line);
            }
            let force = *force.diagnostics();
            evaluations += force.evaluations;
            transfers.push(Transfer {
                source_m,
                field,
                diagnostics: ImpedanceWallCompactSourceDiagnostics { monopole, force },
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
    pub const fn scope(&self) -> &'static str {
        match self.field_kind {
            ImpedanceWallFieldKind::Reflected => IMPEDANCE_WALL_REFLECTED_ARRAY_SCOPE,
            ImpedanceWallFieldKind::Total => IMPEDANCE_WALL_TOTAL_ARRAY_SCOPE,
        }
    }
    pub fn source_count(&self) -> usize {
        self.transfers.len()
    }
    pub fn source_position_m(&self, index: usize) -> Option<[f64; 3]> {
        self.transfers.get(index).map(|t| t.source_m)
    }
    pub fn source_diagnostics(
        &self,
        index: usize,
    ) -> Option<&ImpedanceWallCompactSourceDiagnostics> {
        self.transfers.get(index).map(|t| &t.diagnostics)
    }

    /// Sum coherent peak Q/F phasors on the original source clock, then form
    /// signed active/reactive intensity from the complete pressure and velocity.
    pub fn evaluate(
        &self,
        sources: &[CompactSourcePhasor],
    ) -> Result<CoherentCompactField, AeroacError> {
        self.view().evaluate(sources)
    }

    /// All propagation and source phases are already present. Feed the native
    /// finite-window receiver without adding a ray delay or spectral scaling.
    pub fn evaluate_harmonic_line(
        &self,
        sources: &[CompactSourcePhasor],
    ) -> Result<HarmonicFieldLine<3>, AeroacError> {
        self.view().evaluate_harmonic_line(sources)
    }

    /// Propagate explicitly supplied S=Bξ, E[ξξ*]=I. All Q/F and cross-source
    /// phases survive within each factor; moments add only between independent
    /// factors. Peak-unit factors yield line powers, not PSD/Hz. A one-sided
    /// source PSD over df needs factors of 2 df S. Work is bounded before use.
    pub fn evaluate_coherence_modes(
        &self,
        modes: &[&[CompactSourcePhasor]],
        budget: CoherenceBudget,
    ) -> Result<CompactFieldStatistics, AeroacError> {
        self.view().evaluate_coherence_modes(modes, budget)
    }

    /// O(N) explicitly independent physical sources. Components and paths of
    /// each source stay coherent; an unknown relative phase is not independence.
    pub fn evaluate_independent_sources(
        &self,
        sources: &[CompactSourcePhasor],
    ) -> Result<CompactFieldStatistics, AeroacError> {
        self.view().evaluate_independent_sources(sources)
    }

    fn view(&self) -> ArrayField<'_, ImpedanceWallCompactSourceDiagnostics> {
        ArrayField {
            transfers: &self.transfers,
            frequency_hz: self.frequency_hz,
            scope: self.scope(),
        }
    }
}
