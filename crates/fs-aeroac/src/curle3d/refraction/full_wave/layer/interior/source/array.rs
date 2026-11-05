//! Coherent sources inside the original layer, observed in its incident exterior.
//!
//! Native source-side monopole and force owners retain every internal return
//! and differentiate the actual exterior receiver. Their Cartesian columns are
//! cached once per physical source and applied by the common DD array owner;
//! no reciprocal velocity substitution, extra delay or per-evaluation spectral
//! integration is introduced.

use super::super::super::super::array::{
    ArrayField, FullWaveCompactSourceLimits, Transfer, remaining_budget, set_column,
};
use super::super::super::super::{admit, force, poll, positive};
use super::super::super::PlanarFluidLayer;
use super::super::super::array::FullWaveLayerCompactSourceDiagnostics;
use super::super::reduces_to_interface;
use super::{PreparedInteriorSourceFullWaveForce, PreparedInteriorSourceFullWaveMonopole};
use crate::AeroacError;
use crate::curle3d::statistics::invalid;
use crate::curle3d::{
    CoherenceBudget, CoherentCompactField, CompactFieldStatistics, CompactSourcePhasor,
};
use crate::observation::HarmonicFieldLine;
use fs_math::c64::C64;

/// Original layer-source radiation into the incident exterior, without source loading.
pub const FULL_WAVE_INTERIOR_SOURCE_ARRAY_SCOPE: &str = "Prescribed stationary compact 3D \
volume flows [m³/s] and forces on fluid [N] strictly inside one infinite planar homogeneous \
lossless fluid layer; coherent full-wave pressure and actual particle velocity strictly \
in the original incident exterior. All internal returns and critical, reactive and \
evanescent contributions retain each physical source's phase. Source statistics are \
caller-supplied; slower-than-both trapping layers refuse. Per-source quadrature estimates \
and analytic spectral-tail bounds are not interval certificates or relative array error \
guarantees at interference nulls. No source loading, zero-thickness interior, face trace, \
absorption, finite edges, moving fluids or measured-SPL qualification.";

/// Prepared layer Q/F source array measured in the original incident exterior.
/// Each physical location retains one 4-by-4 Cartesian transfer and native
/// diagnostics. Fixed-order DD application costs O(sources), or O(sources ×
/// factors) for explicit coherence, without allocation or further integration.
#[derive(Clone, Debug)]
pub struct PreparedInteriorSourceFullWaveCompactSourceObserver {
    layer: PlanarFluidLayer,
    frequency_hz: f64,
    observer_m: [f64; 3],
    transfers: Vec<Transfer<FullWaveLayerCompactSourceDiagnostics>>,
    evaluations: usize,
}

impl PreparedInteriorSourceFullWaveCompactSourceObserver {
    /// Prepare sources strictly between the original faces and a receiver
    /// strictly before the first face. Empty arrays still require a nonempty
    /// layer interior, a valid exterior receiver, positive frequency and valid
    /// numerical controls. Per-Q/F limits and a shared evaluation cap bound
    /// preparation; failure or cancellation publishes no partial array.
    ///
    /// A layer identical to the receiving fluid reduces to the first interface,
    /// with the actual layer fluid as source fluid. Its equal-speed case needs
    /// no quadrature, including different exterior density. Equal speeds alone
    /// do not remove both density interfaces. To observe the other exterior,
    /// explicitly reverse the original faces and exterior states.
    pub fn new(
        layer: PlanarFluidLayer,
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
                "interior-source compact array exceeds source or evaluation budget",
            ));
        }
        if !positive(core::f64::consts::TAU * frequency_hz)
            || observer_m.iter().any(|v| !v.is_finite())
        {
            return Err(invalid(
                "interior-source compact array frequency or observer position",
            ));
        }
        if !positive(layer.thickness_m())
            || !positive((-layer.exterior.distance(observer_m)).to_f64())
        {
            return Err(invalid(
                "interior-source compact array requires positive thickness and an observer strictly before the first face",
            ));
        }
        let exact = reduces_to_interface(layer)
            && layer.incident_state().sound_speed_m_s() == layer.layer_state().sound_speed_m_s();
        let mut transfers = Vec::new();
        transfers
            .try_reserve_exact(source_positions_m.len())
            .map_err(|_| {
                invalid("interior-source compact array allocation exceeds available memory")
            })?;
        let mut evaluations = 0;
        for &source_m in source_positions_m {
            poll(cancelled)?;
            let mut monopole_limits = limits.monopole;
            monopole_limits.maximum_evaluations = remaining_budget(
                monopole_limits.maximum_evaluations,
                limits.maximum_total_evaluations - evaluations,
                exact,
            )?;
            let monopole = PreparedInteriorSourceFullWaveMonopole::new(
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
            let force = PreparedInteriorSourceFullWaveForce::new(
                layer,
                frequency_hz,
                source_m,
                observer_m,
                force_limits,
                cancelled,
            )?;
            evaluations += force.diagnostics().evaluations;
            let mut field = [[C64::ZERO; 4]; 4];
            set_column(&mut field, 0, monopole.evaluate_harmonic_line(C64::ONE)?);
            for column in 0..3 {
                let unit = core::array::from_fn(|i| if i == column { C64::ONE } else { C64::ZERO });
                set_column(&mut field, column + 1, force.evaluate_harmonic_line(unit)?);
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
    pub fn source_count(&self) -> usize {
        self.transfers.len()
    }
    pub fn source_position_m(&self, index: usize) -> Option<[f64; 3]> {
        self.transfers.get(index).map(|t| t.source_m)
    }
    pub const fn evaluations(&self) -> usize {
        self.evaluations
    }
    /// Native unit-Q and unit-force errors retain their different physical units.
    /// Multiply Q errors by |Q| and force-vector/tensor errors by ||F||, then
    /// sum over physical sources for an absolute field estimate. Per-source
    /// relative tolerances do not imply relative accuracy at an array null.
    pub fn source_diagnostics(
        &self,
        index: usize,
    ) -> Option<&FullWaveLayerCompactSourceDiagnostics> {
        self.transfers.get(index).map(|t| &t.diagnostics)
    }
    pub const fn scope(&self) -> &'static str {
        FULL_WAVE_INTERIOR_SOURCE_ARRAY_SCOPE
    }

    /// Preserve physical-source and Q/F phase before computing actual p/u flux.
    pub fn evaluate(
        &self,
        sources: &[CompactSourcePhasor],
    ) -> Result<CoherentCompactField, AeroacError> {
        self.view().evaluate(sources)
    }

    /// Original complex peak-source clock; all propagation is already present.
    pub fn evaluate_harmonic_line(
        &self,
        sources: &[CompactSourcePhasor],
    ) -> Result<HarmonicFieldLine<3>, AeroacError> {
        self.view().evaluate_harmonic_line(sources)
    }

    /// Explicit source factors S=Bξ with E[ξξ*]=I and physical peak-amplitude units.
    /// Moments add only across independent factors; Q/F and internal-path phases
    /// remain coherent inside each factor. Results are line powers, not PSD/Hz;
    /// a source PSD in width df requires factors of 2 df S.
    pub fn evaluate_coherence_modes(
        &self,
        modes: &[&[CompactSourcePhasor]],
        budget: CoherenceBudget,
    ) -> Result<CompactFieldStatistics, AeroacError> {
        self.view().evaluate_coherence_modes(modes, budget)
    }

    /// Independent physical sources retain coherent Q/F and all internal returns
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
