//! Coherent and statistically specified compact sources at the same fluid plane.
//!
//! Preparation folds the original volume and force Green functions into one
//! 4-by-4 transfer per physical source. Application preserves every Q/F and
//! cross-source phase before forming quadratic observations. A full-wave field
//! has no single ray delay; all propagation phase is already in these kernels.

use super::reflection::{PreparedReflectedFullWaveForce, PreparedReflectedFullWaveMonopole};
use super::{
    FullWaveForceDiagnostics, FullWaveForceLimits, FullWaveMonopoleDiagnostics,
    FullWaveMonopoleLimits, PlanarFluidInterface, PreparedFullWaveForce, PreparedFullWaveMonopole,
    admit, force, poll, positive,
};
use crate::AeroacError;
use crate::curle3d::statistics::{ComplexSum, StatisticsSum, finite, invalid, nonfinite};
use crate::curle3d::{
    CoherenceBudget, CoherentCompactField, CompactFieldStatistics, CompactSourcePhasor,
};
use crate::observation::HarmonicFieldLine;
use fs_math::c64::C64;

mod recorded;

/// Applicability of both the linear field and its explicitly supplied statistics.
pub const FULL_WAVE_COMPACT_ARRAY_SCOPE: &str = "Prescribed stationary compact 3D volume \
flows [m³/s] and forces on fluid [N]; coherent full-wave transmission across one infinite \
plane between two stationary homogeneous lossless fluids, including reactive and \
evanescent fields. Source statistics are caller-supplied. Numerical quadrature estimates \
and analytic spectral-tail bounds are per source, not interval certificates; no flow, \
absorption, finite boundaries, source feedback or measured-SPL qualification.";

/// Incident-side returning fields without the original direct contribution.
pub const FULL_WAVE_REFLECTED_ARRAY_SCOPE: &str = "Prescribed stationary compact 3D volume \
flows [m³/s] and forces on fluid [N]; coherent reflected-only full-wave field on the \
incident side of one infinite plane between original stationary homogeneous lossless \
fluids, including reactive and evanescent fields. Source statistics are caller-supplied. \
Per-source quadrature estimates and analytic spectral-tail bounds are not interval \
certificates; no direct field, source feedback, absorption, finite edges, moving fluids \
or measured-SPL qualification.";

/// Direct and reflected amplitudes remain coherent within every original source.
pub const FULL_WAVE_TOTAL_ARRAY_SCOPE: &str = "Prescribed stationary compact 3D volume \
flows [m³/s] and forces on fluid [N]; coherent direct-plus-reflected full-wave field on \
the incident side of one infinite plane between original stationary homogeneous lossless \
fluids. Source statistics are caller-supplied. Reflected-transfer error estimates give no \
relative guarantee at direct/reflected or cross-source interference nulls. No singular \
direct self-field, source feedback, absorption, finite edges, moving fluids, interval \
certificate or measured-SPL qualification.";

/// Which physical field the cached array and all of its statistical views contain.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FullWaveCompactFieldKind {
    /// Original sources in the incident fluid, observer in the receiving fluid.
    Transmitted,
    /// Original sources and observer in the incident fluid; returning field only.
    Reflected,
    /// Original sources and observer in the incident fluid; direct plus returning.
    IncidentTotal,
}

impl FullWaveCompactFieldKind {
    const fn scope(self) -> &'static str {
        match self {
            Self::Transmitted => FULL_WAVE_COMPACT_ARRAY_SCOPE,
            Self::Reflected => FULL_WAVE_REFLECTED_ARRAY_SCOPE,
            Self::IncidentTotal => FULL_WAVE_TOTAL_ARRAY_SCOPE,
        }
    }
}

/// Preparation limits for an entire array and for each physical source mechanism.
#[derive(Clone, Copy, Debug)]
pub struct FullWaveCompactSourceLimits {
    /// Maximum number of original source locations, including coincident sources.
    pub max_sources: usize,
    /// Positive total integrand-evaluation cap shared by every Q and F kernel.
    /// Each kernel also obeys its own cap. Exact equal-speed kernels cost zero.
    pub maximum_total_evaluations: usize,
    pub monopole: FullWaveMonopoleLimits,
    pub force: FullWaveForceLimits,
}

/// Unit-amplitude error information for one original source's two mechanisms.
/// The separate units must be retained: multiply Q errors by |Q| and force
/// vector/tensor norm errors by ||F||, then sum over sources for an absolute
/// field estimate. Relative per-source accuracy does not guarantee relative
/// accuracy at a coherently cancelling array null. Incident-side diagnostics
/// control only the reflected transfer, including when the prepared field is total.
#[derive(Clone, Copy, Debug)]
pub struct FullWaveCompactSourceDiagnostics {
    pub monopole: FullWaveMonopoleDiagnostics,
    pub force: FullWaveForceDiagnostics,
}

#[derive(Clone, Debug)]
pub(super) struct Transfer<D> {
    pub(super) source_m: [f64; 3],
    // Rows p,ux,uy,uz; columns Q,Fx,Fy,Fz in their original physical units.
    pub(super) field: [[C64; 4]; 4],
    pub(super) diagnostics: D,
}

/// Full-wave source-array observation with fixed, allocation-free application.
///
/// Each source may contain coherent volume and force amplitudes. Distinct
/// sources may be coherent, explicitly independent, or supplied as low-rank
/// coherence factors. No source is converted to a plane wave and no relative
/// phase is inferred. Preparation owns all quadrature and bounded allocation.
#[derive(Clone, Debug)]
pub struct PreparedFullWaveCompactSourceObserver {
    interface: PlanarFluidInterface,
    field_kind: FullWaveCompactFieldKind,
    frequency_hz: f64,
    observer_m: [f64; 3],
    transfers: Vec<Transfer<FullWaveCompactSourceDiagnostics>>,
    evaluations: usize,
}

impl PreparedFullWaveCompactSourceObserver {
    /// Prepare sources strictly in the incident half-space and one receiver
    /// strictly in the receiving half-space. The shared cap includes discarded
    /// adaptive panels. Exhaustion or cancellation publishes no partial array.
    ///
    /// The absolute tolerances have different Q/F units and apply per source;
    /// preparation cannot infer a physical source amplitude or a cancellation
    /// error budget. An empty array is admitted as silence after validating the
    /// observer, frequency and controls, and consumes no quadrature work.
    pub fn new(
        interface: PlanarFluidInterface,
        frequency_hz: f64,
        observer_m: [f64; 3],
        source_positions_m: &[[f64; 3]],
        limits: FullWaveCompactSourceLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        Self::prepare(
            interface,
            FullWaveCompactFieldKind::Transmitted,
            frequency_hz,
            observer_m,
            source_positions_m,
            limits,
            cancelled,
        )
    }

    /// Prepare only the returning field with every source and observer strictly
    /// in the original incident fluid. Coincident sources/observer are admitted:
    /// their reflected fields are finite, without claiming self-consistent loading.
    /// Limits and diagnostics apply to the reflected Q/F transfer of each source.
    pub fn new_reflected(
        interface: PlanarFluidInterface,
        frequency_hz: f64,
        observer_m: [f64; 3],
        source_positions_m: &[[f64; 3]],
        limits: FullWaveCompactSourceLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        Self::prepare(
            interface,
            FullWaveCompactFieldKind::Reflected,
            frequency_hz,
            observer_m,
            source_positions_m,
            limits,
            cancelled,
        )
    }

    /// Prepare the coherent direct-plus-reflected field in the incident fluid.
    /// Each original Q/F column already contains both propagation paths. Even
    /// explicitly independent sources retain their own direct/reflected cross
    /// terms. A coincident source refuses its singular direct field at admission,
    /// irrespective of later source amplitude. Numerical limits control only
    /// the reflected transfer; total fields inherit absolute, not relative,
    /// error estimates near an interference null.
    pub fn new_total(
        interface: PlanarFluidInterface,
        frequency_hz: f64,
        observer_m: [f64; 3],
        source_positions_m: &[[f64; 3]],
        limits: FullWaveCompactSourceLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        Self::prepare(
            interface,
            FullWaveCompactFieldKind::IncidentTotal,
            frequency_hz,
            observer_m,
            source_positions_m,
            limits,
            cancelled,
        )
    }

    fn prepare(
        interface: PlanarFluidInterface,
        field_kind: FullWaveCompactFieldKind,
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
                "full-wave compact array exceeds source or evaluation budget",
            ));
        }
        if !positive(core::f64::consts::TAU * frequency_hz)
            || observer_m.iter().any(|v| !v.is_finite())
        {
            return Err(invalid(
                "full-wave compact array frequency or observer position",
            ));
        }
        let observer_distance = interface.distance(observer_m).to_f64();
        let correct_side = match field_kind {
            FullWaveCompactFieldKind::Transmitted => positive(observer_distance),
            _ => positive(-observer_distance),
        };
        if !correct_side {
            return Err(invalid(
                "full-wave compact array frequency or observer half-space",
            ));
        }
        let exact = interface.incident_state().sound_speed_m_s()
            == interface.receiving_state().sound_speed_m_s();
        let mut transfers = Vec::new();
        transfers
            .try_reserve_exact(source_positions_m.len())
            .map_err(|_| invalid("full-wave compact array allocation exceeds available memory"))?;
        let mut evaluations = 0;
        for &source_m in source_positions_m {
            poll(cancelled)?;
            if field_kind == FullWaveCompactFieldKind::IncidentTotal && source_m == observer_m {
                return Err(invalid(
                    "total full-wave compact array has a singular direct self-field",
                ));
            }
            let mut monopole_limits = limits.monopole;
            monopole_limits.maximum_evaluations = remaining_budget(
                monopole_limits.maximum_evaluations,
                limits.maximum_total_evaluations - evaluations,
                exact,
            )?;
            let mut field = [[C64::ZERO; 4]; 4];
            let monopole_diagnostics = if field_kind == FullWaveCompactFieldKind::Transmitted {
                let monopole = PreparedFullWaveMonopole::new(
                    interface,
                    frequency_hz,
                    source_m,
                    observer_m,
                    monopole_limits,
                    cancelled,
                )?;
                for (row, coefficient) in field.iter_mut().zip(monopole.transfer) {
                    row[0] = coefficient;
                }
                *monopole.diagnostics()
            } else {
                let monopole = PreparedReflectedFullWaveMonopole::new(
                    interface,
                    frequency_hz,
                    source_m,
                    observer_m,
                    monopole_limits,
                    cancelled,
                )?;
                let line = if field_kind == FullWaveCompactFieldKind::Reflected {
                    monopole.evaluate_reflected_harmonic_line(C64::ONE)?
                } else {
                    monopole.evaluate_total_harmonic_line(C64::ONE)?
                };
                set_column(&mut field, 0, line);
                *monopole.diagnostics()
            };
            evaluations += monopole_diagnostics.evaluations;
            let mut force_limits = limits.force;
            force_limits.maximum_evaluations = remaining_budget(
                force_limits.maximum_evaluations,
                limits.maximum_total_evaluations - evaluations,
                exact,
            )?;
            let force_diagnostics = if field_kind == FullWaveCompactFieldKind::Transmitted {
                let force = PreparedFullWaveForce::new(
                    interface,
                    frequency_hz,
                    source_m,
                    observer_m,
                    force_limits,
                    cancelled,
                )?;
                for (row, coefficients) in field.iter_mut().zip(force.transfer) {
                    row[1..].copy_from_slice(&coefficients);
                }
                *force.diagnostics()
            } else {
                let force = PreparedReflectedFullWaveForce::new(
                    interface,
                    frequency_hz,
                    source_m,
                    observer_m,
                    force_limits,
                    cancelled,
                )?;
                for column in 0..3 {
                    let unit =
                        core::array::from_fn(|i| if i == column { C64::ONE } else { C64::ZERO });
                    let line = if field_kind == FullWaveCompactFieldKind::Reflected {
                        force.evaluate_reflected_harmonic_line(unit)?
                    } else {
                        force.evaluate_total_harmonic_line(unit)?
                    };
                    set_column(&mut field, column + 1, line);
                }
                *force.diagnostics()
            };
            evaluations += force_diagnostics.evaluations;
            transfers.push(Transfer {
                source_m,
                field,
                diagnostics: FullWaveCompactSourceDiagnostics {
                    monopole: monopole_diagnostics,
                    force: force_diagnostics,
                },
            });
        }
        poll(cancelled)?;
        Ok(Self {
            interface,
            field_kind,
            frequency_hz,
            observer_m,
            transfers,
            evaluations,
        })
    }

    pub const fn interface(&self) -> PlanarFluidInterface {
        self.interface
    }
    pub const fn field_kind(&self) -> FullWaveCompactFieldKind {
        self.field_kind
    }
    pub const fn scope(&self) -> &'static str {
        self.field_kind.scope()
    }
    pub const fn frequency_hz(&self) -> f64 {
        self.frequency_hz
    }
    pub const fn observer_position_m(&self) -> [f64; 3] {
        self.observer_m
    }
    pub fn source_count(&self) -> usize {
        self.transfers.len()
    }
    pub fn source_position_m(&self, index: usize) -> Option<[f64; 3]> {
        self.transfers.get(index).map(|t| t.source_m)
    }
    pub const fn evaluations(&self) -> usize {
        self.evaluations
    }
    pub fn source_diagnostics(&self, index: usize) -> Option<&FullWaveCompactSourceDiagnostics> {
        self.transfers.get(index).map(|t| &t.diagnostics)
    }

    /// Coherent peak Q/F phasors on the original source clock. All components
    /// and sources are accumulated before active/reactive intensity is formed.
    /// Double-double reductions retain residuals between strong contributions.
    pub fn evaluate(
        &self,
        sources: &[CompactSourcePhasor],
    ) -> Result<CoherentCompactField, AeroacError> {
        self.view().evaluate(sources)
    }

    /// Feed finite-window observation with all source and propagation phases
    /// already present; no additional ray delay or spectral normalization.
    pub fn evaluate_harmonic_line(
        &self,
        sources: &[CompactSourcePhasor],
    ) -> Result<HarmonicFieldLine<3>, AeroacError> {
        self.view().evaluate_harmonic_line(sources)
    }

    /// Propagate explicit coherence factors S=Bξ, E[ξξ*]=I. Each factor keeps
    /// all Q/F and cross-source phases; quadratic moments add only across the
    /// declared independent factors. Work is O(sources × factors), bounded
    /// before evaluation. Factors use peak amplitudes; these are line powers,
    /// not PSD/Hz. A one-sided source PSD over df requires factors of 2 df S.
    pub fn evaluate_coherence_modes(
        &self,
        modes: &[&[CompactSourcePhasor]],
        budget: CoherenceBudget,
    ) -> Result<CompactFieldStatistics, AeroacError> {
        self.view().evaluate_coherence_modes(modes, budget)
    }

    /// O(N) observation of explicitly independent physical sources. Q and F
    /// remain coherent within each source. Unknown phase alone does not imply
    /// independence; the caller must supply that statistical assumption.
    pub fn evaluate_independent_sources(
        &self,
        sources: &[CompactSourcePhasor],
    ) -> Result<CompactFieldStatistics, AeroacError> {
        self.view().evaluate_independent_sources(sources)
    }

    fn view(&self) -> ArrayField<'_, FullWaveCompactSourceDiagnostics> {
        ArrayField {
            transfers: &self.transfers,
            frequency_hz: self.frequency_hz,
            scope: self.scope(),
        }
    }
}

/// Borrow the same fixed Cartesian transfers for every propagation owner.
/// Diagnostics retain their native units/type but do not enter hot evaluation.
pub(super) struct ArrayField<'a, D> {
    pub(super) transfers: &'a [Transfer<D>],
    pub(super) frequency_hz: f64,
    pub(super) scope: &'static str,
}

impl<D> ArrayField<'_, D> {
    pub(super) fn evaluate(
        &self,
        sources: &[CompactSourcePhasor],
    ) -> Result<CoherentCompactField, AeroacError> {
        let [pressure, ux, uy, uz] = self.amplitudes(sources)?;
        let velocity = [ux, uy, uz];
        let intensity = velocity.map(|v| pressure.scale(0.5) * v.conj());
        if intensity.iter().any(|v| !finite(*v)) {
            return Err(nonfinite("full-wave compact array intensity"));
        }
        Ok(CoherentCompactField {
            pressure_pa: pressure,
            particle_velocity_m_s: velocity,
            active_intensity_w_m2: intensity.map(|v| v.re),
            reactive_intensity_w_m2: intensity.map(|v| v.im),
            scope: self.scope,
        })
    }

    pub(super) fn evaluate_harmonic_line(
        &self,
        sources: &[CompactSourcePhasor],
    ) -> Result<HarmonicFieldLine<3>, AeroacError> {
        let [pressure, ux, uy, uz] = self.amplitudes(sources)?;
        Ok(HarmonicFieldLine {
            frequency_hz: self.frequency_hz,
            pressure_pa: pressure,
            particle_velocity_m_s: [ux, uy, uz],
        })
    }

    pub(super) fn evaluate_coherence_modes(
        &self,
        modes: &[&[CompactSourcePhasor]],
        budget: CoherenceBudget,
    ) -> Result<CompactFieldStatistics, AeroacError> {
        if modes.len() > budget.max_modes
            || modes
                .len()
                .checked_mul(self.transfers.len())
                .is_none_or(|pairs| pairs > budget.max_source_mode_pairs)
        {
            return Err(invalid("full-wave compact coherence exceeds work budget"));
        }
        if modes.iter().any(|mode| mode.len() != self.transfers.len()) {
            return Err(invalid(
                "each full-wave coherence factor requires one vector per source",
            ));
        }
        let mut sum = StatisticsSum::default();
        for mode in modes {
            sum.add(self.amplitudes(mode)?)?;
        }
        sum.finish(self.scope)
    }

    pub(super) fn evaluate_independent_sources(
        &self,
        sources: &[CompactSourcePhasor],
    ) -> Result<CompactFieldStatistics, AeroacError> {
        self.check_count(sources)?;
        let mut sum = StatisticsSum::default();
        for (transfer, source) in self.transfers.iter().zip(sources) {
            let mut field = [ComplexSum::default(); 4];
            transfer.accumulate(*source, &mut field)?;
            sum.add(field.map(ComplexSum::value))?;
        }
        sum.finish(self.scope)
    }

    fn check_count(&self, sources: &[CompactSourcePhasor]) -> Result<(), AeroacError> {
        if sources.len() != self.transfers.len() {
            return Err(invalid(
                "one vector is required per full-wave compact source",
            ));
        }
        Ok(())
    }

    fn amplitudes(&self, sources: &[CompactSourcePhasor]) -> Result<[C64; 4], AeroacError> {
        self.check_count(sources)?;
        let mut sums = [ComplexSum::default(); 4];
        for (transfer, source) in self.transfers.iter().zip(sources) {
            transfer.accumulate(*source, &mut sums)?;
        }
        let field = sums.map(ComplexSum::value);
        if field.iter().any(|v| !finite(*v)) {
            return Err(nonfinite("full-wave compact array amplitude"));
        }
        Ok(field)
    }
}

pub(super) fn set_column(field: &mut [[C64; 4]; 4], column: usize, line: HarmonicFieldLine<3>) {
    field[0][column] = line.pressure_pa;
    for axis in 0..3 {
        field[axis + 1][column] = line.particle_velocity_m_s[axis];
    }
}

pub(super) fn remaining_budget(
    requested: usize,
    remaining: usize,
    exact: bool,
) -> Result<usize, AeroacError> {
    if remaining == 0 && !exact {
        return Err(invalid(
            "full-wave compact array exhausted its shared evaluation budget",
        ));
    }
    // Constructors require a positive cap; the exact branch uses no samples.
    Ok(requested.min(remaining.max(1)))
}

impl<D> Transfer<D> {
    fn accumulate(
        &self,
        source: CompactSourcePhasor,
        sums: &mut [ComplexSum; 4],
    ) -> Result<(), AeroacError> {
        let amplitudes = [
            source.volume_flow_m3_s,
            source.force_on_fluid_n[0],
            source.force_on_fluid_n[1],
            source.force_on_fluid_n[2],
        ];
        for (column, amplitude) in amplitudes.into_iter().enumerate() {
            if !finite(amplitude) {
                return Err(nonfinite("full-wave compact source amplitude"));
            }
            if amplitude != C64::ZERO {
                for (row, sum) in self.field.iter().zip(sums.iter_mut()) {
                    sum.add(row[column] * amplitude)?;
                }
            }
        }
        Ok(())
    }
}
