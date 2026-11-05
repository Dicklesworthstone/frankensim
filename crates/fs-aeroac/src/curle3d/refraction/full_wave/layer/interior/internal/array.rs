//! Coherent original Q/F sources and receivers inside the same fluid layer.
//!
//! Native internal monopole and force owners supply every reflected return
//! and, for total observations, the analytic direct field. Their Cartesian
//! columns are cached once per physical source and applied by the shared DD
//! array owner. Direct and returning paths never become independent sources.

use super::super::super::super::array::{
    ArrayField, FullWaveCompactSourceLimits, Transfer, remaining_budget, set_column,
};
use super::super::super::super::{admit, force, poll, positive};
use super::super::super::PlanarFluidLayer;
use super::super::super::array::FullWaveLayerCompactSourceDiagnostics;
use super::{PreparedLayerInternalFullWaveForce, PreparedLayerInternalFullWaveMonopole};
use crate::AeroacError;
use crate::curle3d::statistics::invalid;
use crate::curle3d::{
    CoherenceBudget, CoherentCompactField, CompactFieldStatistics, CompactSourcePhasor,
};
use crate::observation::HarmonicFieldLine;
use fs_math::{c64::C64, dd::Dd};

pub const FULL_WAVE_LAYER_INTERNAL_TOTAL_ARRAY_SCOPE: &str = "Prescribed stationary \
compact 3D volume flows [m³/s] and forces ON fluid [N], with sources and receiver \
strictly inside one infinite planar homogeneous lossless fluid layer. Native direct \
fields and all coherent reflected, critical, reactive and evanescent contributions \
retain each physical source's phase. Source statistics are caller-supplied. Per-source \
reflected-transfer estimates and analytic tails are not interval certificates or \
relative array accuracy guarantees at interference nulls. Singular total self-fields \
and slower-than-both guided layers refuse; no source loading, absorption, moving \
fluids, finite edges or measured-SPL qualification.";

pub const FULL_WAVE_LAYER_INTERNAL_REFLECTED_ARRAY_SCOPE: &str = "Prescribed stationary \
compact 3D volume flows [m³/s] and forces ON fluid [N], with sources and receiver \
strictly inside one infinite planar homogeneous lossless fluid layer. Complete coherent \
reflected pressure and actual particle velocity retain every internal return, critical, \
reactive and evanescent contribution; the direct field is excluded. Finite reflected \
self-fields do not imply source loading. Source statistics are caller-supplied; \
slower-than-both guided layers refuse. Per-source numerical estimates and analytic \
tails are not interval certificates or relative array accuracy guarantees at nulls; \
no absorption, moving fluids, finite edges or measured-SPL qualification.";

/// Original layer-interior Q/F arrays with one Cartesian 4-by-4 transfer per source.
/// Fixed-order DD application costs O(sources), or O(sources × coherence factors),
/// without allocation, additional quadrature or a second propagation delay.
#[derive(Clone, Debug)]
pub struct PreparedLayerInternalFullWaveCompactSourceObserver {
    layer: PlanarFluidLayer,
    frequency_hz: f64,
    observer_m: [f64; 3],
    include_direct: bool,
    transfers: Vec<Transfer<FullWaveLayerCompactSourceDiagnostics>>,
    evaluations: usize,
}

impl PreparedLayerInternalFullWaveCompactSourceObserver {
    /// Prepare the coherent direct-plus-returning field of original Q/F sources.
    /// Both sources and observer lie strictly between the original layer faces.
    /// A source at the receiver refuses before that source's spectral preparation.
    /// Empty arrays still validate the original DD receiver domain, frequency
    /// and numerical controls. No partially prepared array is published.
    pub fn new(
        layer: PlanarFluidLayer,
        frequency_hz: f64,
        observer_m: [f64; 3],
        source_positions_m: &[[f64; 3]],
        limits: FullWaveCompactSourceLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        Self::prepare(
            layer,
            frequency_hz,
            observer_m,
            source_positions_m,
            true,
            limits,
            cancelled,
        )
    }

    /// Prepare only the coherent returning field. Coincident sources retain
    /// their finite boundary response; this does not solve source loading.
    pub fn new_reflected(
        layer: PlanarFluidLayer,
        frequency_hz: f64,
        observer_m: [f64; 3],
        source_positions_m: &[[f64; 3]],
        limits: FullWaveCompactSourceLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        Self::prepare(
            layer,
            frequency_hz,
            observer_m,
            source_positions_m,
            false,
            limits,
            cancelled,
        )
    }

    fn prepare(
        layer: PlanarFluidLayer,
        frequency_hz: f64,
        observer_m: [f64; 3],
        source_positions_m: &[[f64; 3]],
        include_direct: bool,
        limits: FullWaveCompactSourceLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        poll(cancelled)?;
        admit(limits.monopole)?;
        force::admit(limits.force)?;
        if source_positions_m.len() > limits.max_sources || limits.maximum_total_evaluations == 0 {
            return Err(invalid(
                "internal compact array exceeds source or evaluation budget",
            ));
        }
        if !positive(core::f64::consts::TAU * frequency_hz)
            || observer_m.iter().any(|v| !v.is_finite())
        {
            return Err(invalid(
                "internal compact array frequency or observer position",
            ));
        }
        let depth = layer.exterior.distance(observer_m);
        let remaining = Dd::from_f64(layer.thickness_m()) - depth;
        if !positive(depth.to_f64()) || !positive(remaining.to_f64()) {
            return Err(invalid(
                "internal compact array observer must lie strictly between the layer faces",
            ));
        }
        // Match the native owners exactly. Equal sound speeds alone do not
        // remove two density interfaces or their coherent internal returns.
        let exact = layer.layer_state() == layer.receiving_state()
            && layer.incident_state().sound_speed_m_s() == layer.layer_state().sound_speed_m_s();
        let mut transfers = Vec::new();
        transfers
            .try_reserve_exact(source_positions_m.len())
            .map_err(|_| invalid("internal compact array allocation exceeds available memory"))?;
        let mut evaluations = 0;
        for &source_m in source_positions_m {
            poll(cancelled)?;
            if include_direct && source_m == observer_m {
                return Err(invalid(
                    "total internal compact array has a singular direct self-field",
                ));
            }
            let mut monopole_limits = limits.monopole;
            monopole_limits.maximum_evaluations = remaining_budget(
                monopole_limits.maximum_evaluations,
                limits.maximum_total_evaluations - evaluations,
                exact,
            )?;
            let monopole = PreparedLayerInternalFullWaveMonopole::new(
                layer,
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
            let force = PreparedLayerInternalFullWaveForce::new(
                layer,
                frequency_hz,
                source_m,
                observer_m,
                force_limits,
                cancelled,
            )?;
            evaluations += force.diagnostics().evaluations;
            let mut field = [[C64::ZERO; 4]; 4];
            set_column(
                &mut field,
                0,
                if include_direct {
                    monopole.evaluate_total_harmonic_line(C64::ONE)?
                } else {
                    monopole.evaluate_reflected_harmonic_line(C64::ONE)?
                },
            );
            for column in 0..3 {
                let unit = core::array::from_fn(|i| if i == column { C64::ONE } else { C64::ZERO });
                let line = if include_direct {
                    force.evaluate_total_harmonic_line(unit)?
                } else {
                    force.evaluate_reflected_harmonic_line(unit)?
                };
                set_column(&mut field, column + 1, line);
            }
            transfers.push(Transfer {
                source_m,
                field,
                diagnostics: FullWaveLayerCompactSourceDiagnostics {
                    monopole: *monopole.diagnostics(),
                    force: *force.diagnostics(),
                },
            });
        }
        poll(cancelled)?;
        Ok(Self {
            layer,
            frequency_hz,
            observer_m,
            include_direct,
            transfers,
            evaluations,
        })
    }

    pub const fn layer(&self) -> PlanarFluidLayer {
        self.layer
    }
    pub const fn frequency_hz(&self) -> f64 {
        self.frequency_hz
    }
    pub const fn observer_position_m(&self) -> [f64; 3] {
        self.observer_m
    }
    pub const fn includes_direct_field(&self) -> bool {
        self.include_direct
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
    /// Native reflected-transfer estimates per unit Q and per unit force retain
    /// their different physical units. Multiply by the physical amplitudes and
    /// sum absolute estimates across sources; no relative null guarantee follows.
    pub fn source_diagnostics(
        &self,
        index: usize,
    ) -> Option<&FullWaveLayerCompactSourceDiagnostics> {
        self.transfers.get(index).map(|t| &t.diagnostics)
    }
    pub const fn scope(&self) -> &'static str {
        if self.include_direct {
            FULL_WAVE_LAYER_INTERNAL_TOTAL_ARRAY_SCOPE
        } else {
            FULL_WAVE_LAYER_INTERNAL_REFLECTED_ARRAY_SCOPE
        }
    }

    /// Preserve Q/F, direct/return and cross-source phase before forming p/u flux.
    pub fn evaluate(
        &self,
        sources: &[CompactSourcePhasor],
    ) -> Result<CoherentCompactField, AeroacError> {
        self.view().evaluate(sources)
    }

    /// Original complex peak-source clock, with all propagation already included.
    pub fn evaluate_harmonic_line(
        &self,
        sources: &[CompactSourcePhasor],
    ) -> Result<HarmonicFieldLine<3>, AeroacError> {
        self.view().evaluate_harmonic_line(sources)
    }

    /// Caller-supplied factors S=Bξ with E[ξξ*]=I in physical peak-amplitude units.
    /// Moments add only across independent factors; Q/F and all path phases stay
    /// coherent inside each factor. These are line powers, not PSD/Hz; source
    /// PSD in width df requires factors of 2 df S.
    pub fn evaluate_coherence_modes(
        &self,
        modes: &[&[CompactSourcePhasor]],
        budget: CoherenceBudget,
    ) -> Result<CompactFieldStatistics, AeroacError> {
        self.view().evaluate_coherence_modes(modes, budget)
    }

    /// Explicitly independent physical sources retain coherent Q/F and all paths
    /// within each source. Unknown phase alone does not imply independence.
    pub fn evaluate_independent_sources(
        &self,
        sources: &[CompactSourcePhasor],
    ) -> Result<CompactFieldStatistics, AeroacError> {
        self.view().evaluate_independent_sources(sources)
    }

    fn view(&self) -> ArrayField<'_, FullWaveLayerCompactSourceDiagnostics> {
        ArrayField {
            transfers: &self.transfers,
            frequency_hz: self.frequency_hz,
            scope: self.scope(),
        }
    }
}
