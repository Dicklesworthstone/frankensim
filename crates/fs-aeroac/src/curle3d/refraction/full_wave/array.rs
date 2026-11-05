//! Coherent and statistically specified compact sources across the same fluid plane.
//!
//! Preparation folds the original volume and force Green functions into one
//! 4-by-4 transfer per physical source. Application preserves every Q/F and
//! cross-source phase before forming quadratic observations. A full-wave field
//! has no single ray delay; all propagation phase is already in these kernels.

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

/// Applicability of both the linear field and its explicitly supplied statistics.
pub const FULL_WAVE_COMPACT_ARRAY_SCOPE: &str = "Prescribed stationary compact 3D volume \
flows [m³/s] and forces on fluid [N]; coherent full-wave transmission across one infinite \
plane between two stationary homogeneous lossless fluids, including reactive and \
evanescent fields. Source statistics are caller-supplied. Numerical quadrature estimates \
and analytic spectral-tail bounds are per source, not interval certificates; no flow, \
absorption, finite boundaries, source feedback or measured-SPL qualification.";

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
/// accuracy at a coherently cancelling array null.
#[derive(Clone, Copy, Debug)]
pub struct FullWaveCompactSourceDiagnostics {
    pub monopole: FullWaveMonopoleDiagnostics,
    pub force: FullWaveForceDiagnostics,
}

#[derive(Clone, Debug)]
struct Transfer {
    source_m: [f64; 3],
    // Rows p,ux,uy,uz; columns Q,Fx,Fy,Fz in their original physical units.
    field: [[C64; 4]; 4],
    diagnostics: FullWaveCompactSourceDiagnostics,
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
    frequency_hz: f64,
    observer_m: [f64; 3],
    transfers: Vec<Transfer>,
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
            || !positive(interface.distance(observer_m).to_f64())
        {
            return Err(invalid(
                "full-wave compact array frequency or receiving position",
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
            let mut monopole_limits = limits.monopole;
            monopole_limits.maximum_evaluations = remaining_budget(
                monopole_limits.maximum_evaluations,
                limits.maximum_total_evaluations - evaluations,
                exact,
            )?;
            let monopole = PreparedFullWaveMonopole::new(
                interface,
                frequency_hz,
                source_m,
                observer_m,
                monopole_limits,
                cancelled,
            )?;
            evaluations += monopole.diagnostics().evaluations;
            let mut force_limits = limits.force;
            force_limits.maximum_evaluations = remaining_budget(
                force_limits.maximum_evaluations,
                limits.maximum_total_evaluations - evaluations,
                exact,
            )?;
            let force = PreparedFullWaveForce::new(
                interface,
                frequency_hz,
                source_m,
                observer_m,
                force_limits,
                cancelled,
            )?;
            evaluations += force.diagnostics().evaluations;
            let field = core::array::from_fn(|row| {
                [
                    monopole.transfer[row],
                    force.transfer[row][0],
                    force.transfer[row][1],
                    force.transfer[row][2],
                ]
            });
            transfers.push(Transfer {
                source_m,
                field,
                diagnostics: FullWaveCompactSourceDiagnostics {
                    monopole: *monopole.diagnostics(),
                    force: *force.diagnostics(),
                },
            });
        }
        poll(cancelled)?;
        Ok(Self {
            interface,
            frequency_hz,
            observer_m,
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
            scope: FULL_WAVE_COMPACT_ARRAY_SCOPE,
        })
    }

    /// Feed finite-window observation with all source and propagation phases
    /// already present; no additional ray delay or spectral normalization.
    pub fn evaluate_harmonic_line(
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
        if modes.len() > budget.max_modes
            || modes
                .len()
                .checked_mul(self.source_count())
                .is_none_or(|pairs| pairs > budget.max_source_mode_pairs)
        {
            return Err(invalid("full-wave compact coherence exceeds work budget"));
        }
        if modes.iter().any(|mode| mode.len() != self.source_count()) {
            return Err(invalid(
                "each full-wave coherence factor requires one vector per source",
            ));
        }
        let mut sum = StatisticsSum::default();
        for mode in modes {
            sum.add(self.amplitudes(mode)?)?;
        }
        sum.finish(FULL_WAVE_COMPACT_ARRAY_SCOPE)
    }

    /// O(N) observation of explicitly independent physical sources. Q and F
    /// remain coherent within each source. Unknown phase alone does not imply
    /// independence; the caller must supply that statistical assumption.
    pub fn evaluate_independent_sources(
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
        sum.finish(FULL_WAVE_COMPACT_ARRAY_SCOPE)
    }

    fn check_count(&self, sources: &[CompactSourcePhasor]) -> Result<(), AeroacError> {
        if sources.len() != self.source_count() {
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

fn remaining_budget(requested: usize, remaining: usize, exact: bool) -> Result<usize, AeroacError> {
    if remaining == 0 && !exact {
        return Err(invalid(
            "full-wave compact array exhausted its shared evaluation budget",
        ));
    }
    // Constructors require a positive cap; the exact branch uses no samples.
    Ok(requested.min(remaining.max(1)))
}

impl Transfer {
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
