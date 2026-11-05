//! Original compact-source coherence through the complete finite-layer spectrum.
//!
//! Preparation combines the native layered Q/F kernels into one Cartesian
//! 4-by-4 transfer per physical source. Internal returns stay coherent with
//! each other and with that source's complete pressure/velocity field. The
//! shared array application retains Q/F and cross-source phases without another
//! quadrature, propagation delay, allocation, or invented phase independence.

use super::super::array::{
    ArrayField, FullWaveCompactFieldKind, FullWaveCompactSourceLimits, Transfer, remaining_budget,
    set_column,
};
use super::super::{admit, force, poll, positive};
use super::reflection::{
    PreparedLayeredReflectedFullWaveForce, PreparedLayeredReflectedFullWaveMonopole,
    reduces_to_first_interface,
};
use super::{
    FullWaveLayerForceDiagnostics, FullWaveLayerMonopoleDiagnostics, PlanarFluidLayer,
    PreparedLayeredFullWaveForce, PreparedLayeredFullWaveMonopole,
};
use crate::AeroacError;
use crate::curle3d::statistics::invalid;
use crate::curle3d::{
    CoherenceBudget, CoherentCompactField, CompactFieldStatistics, CompactSourcePhasor,
};
use crate::observation::HarmonicFieldLine;
use fs_math::{c64::C64, dd::Dd};

/// All source statistics are explicit, and every original layer limit remains.
pub const FULL_WAVE_LAYER_ARRAY_SCOPE: &str = "Prescribed stationary compact 3D volume \
flows [m³/s] and forces on fluid [N]; coherent full-wave transmission through one \
infinite planar homogeneous lossless fluid layer between original stationary fluids. \
All internal returns, reactive and evanescent components remain coherent. Sources and \
observer are strictly outside opposite faces; slower-than-both trapping layers refuse. \
Source statistics are caller-supplied. Per-source quadrature estimates and analytic \
spectral-tail bounds are not interval certificates or relative array error guarantees \
at interference nulls; no absorption, finite edges, moving fluids, source feedback or \
measured-SPL qualification.";

/// Incident-side returning fields retain every internal layer history.
pub const FULL_WAVE_LAYER_REFLECTED_ARRAY_SCOPE: &str = "Prescribed stationary compact \
3D volume flows [m³/s] and forces on fluid [N]; reflected-only full-wave field on the \
incident side of one infinite planar homogeneous lossless fluid layer between original \
stationary fluids, including coherent internal returns and evanescent/reactive fields. \
Source statistics are caller-supplied; slower-than-both trapping layers refuse. \
Per-source quadrature estimates and analytic spectral-tail bounds are not interval \
certificates; no direct field, source feedback, absorption, finite edges, moving fluids \
or measured-SPL qualification.";

/// Original direct and complete layer-reflected amplitudes precede statistics.
pub const FULL_WAVE_LAYER_TOTAL_ARRAY_SCOPE: &str = "Prescribed stationary compact \
3D volume flows [m³/s] and forces on fluid [N]; coherent direct-plus-reflected full-wave \
field on the incident side of one infinite planar homogeneous lossless fluid layer \
between original stationary fluids. All internal returns retain their phases. Source \
statistics are caller-supplied; slower-than-both trapping layers refuse. Reflected-transfer \
error estimates give no relative guarantee at path or cross-source interference nulls. \
No singular direct self-field, source feedback, absorption, finite edges, moving fluids, \
interval certificate or measured-SPL qualification.";

/// Unit-Q and unit-force errors retain their original, distinct physical units.
/// Multiply Q errors by |Q| and force norm errors by ||F||, then sum across
/// original sources for an absolute field estimate. Unknown amplitudes or
/// coherence cannot be inferred from geometry or the fluid snapshots.
/// Incident-side diagnostics control the reflected transfer, also for total fields.
#[derive(Clone, Copy, Debug)]
pub struct FullWaveLayerCompactSourceDiagnostics {
    pub monopole: FullWaveLayerMonopoleDiagnostics,
    pub force: FullWaveLayerForceDiagnostics,
}

/// Prepared physical source array observed outside the original finite layer.
/// Cached application is O(sources), or O(sources × coherence factors), with
/// the same fixed-order DD accumulation as single-interface observations.
#[derive(Clone, Debug)]
pub struct PreparedLayeredFullWaveCompactSourceObserver {
    layer: PlanarFluidLayer,
    field_kind: FullWaveCompactFieldKind,
    frequency_hz: f64,
    observer_m: [f64; 3],
    transfers: Vec<Transfer<FullWaveLayerCompactSourceDiagnostics>>,
    evaluations: usize,
}

impl PreparedLayeredFullWaveCompactSourceObserver {
    /// Sources must lie strictly before the first face and the receiver strictly
    /// beyond the second face. Original per-Q/F limits and a shared total work
    /// cap apply to preparation. Failure or cancellation publishes no partial
    /// array. Empty input is validated silence, with no quadrature work.
    ///
    /// Equal sound speeds alone do not remove a finite density layer: its
    /// coherent internal returns still require the native spectral kernel.
    /// Zero thickness and three identical fluids inherit exact original limits.
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
            FullWaveCompactFieldKind::Transmitted,
            frequency_hz,
            observer_m,
            source_positions_m,
            limits,
            cancelled,
        )
    }

    /// Prepare the complete returning field with sources and observer strictly
    /// before the first face. Coincident points retain their finite reflected
    /// field without implying a source-loading or feedback model.
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
            FullWaveCompactFieldKind::Reflected,
            frequency_hz,
            observer_m,
            source_positions_m,
            limits,
            cancelled,
        )
    }

    /// Prepare coherent original direct plus layer-reflected transfer columns.
    /// Independence between physical sources never removes their own internal
    /// path interference. A coincident source refuses its singular direct field
    /// at admission, irrespective of later amplitudes. Numerical diagnostics
    /// control reflected transfer only; total fields inherit absolute estimates
    /// without a relative guarantee at an interference null.
    pub fn new_total(
        layer: PlanarFluidLayer,
        frequency_hz: f64,
        observer_m: [f64; 3],
        source_positions_m: &[[f64; 3]],
        limits: FullWaveCompactSourceLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, AeroacError> {
        Self::prepare(
            layer,
            FullWaveCompactFieldKind::IncidentTotal,
            frequency_hz,
            observer_m,
            source_positions_m,
            limits,
            cancelled,
        )
    }

    fn prepare(
        layer: PlanarFluidLayer,
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
                "layered compact array exceeds source or evaluation budget",
            ));
        }
        if !positive(core::f64::consts::TAU * frequency_hz)
            || observer_m.iter().any(|v| !v.is_finite())
        {
            return Err(invalid(
                "layered compact array frequency or observer position",
            ));
        }
        let distance = layer.exterior.distance(observer_m);
        let clearance = match field_kind {
            FullWaveCompactFieldKind::Transmitted => distance - Dd::from_f64(layer.thickness_m()),
            _ => -distance,
        };
        if !positive(clearance.to_f64()) {
            return Err(invalid(
                "layered compact array observer is outside its requested exterior half-space",
            ));
        }
        let reduced = match field_kind {
            FullWaveCompactFieldKind::Transmitted => layer.reduces_to_interface(),
            _ => reduces_to_first_interface(layer),
        };
        let exact = reduced
            && layer.incident_state().sound_speed_m_s()
                == layer.receiving_state().sound_speed_m_s();
        let mut transfers = Vec::new();
        transfers
            .try_reserve_exact(source_positions_m.len())
            .map_err(|_| invalid("layered compact array allocation exceeds available memory"))?;
        let mut evaluations = 0;
        for &source_m in source_positions_m {
            poll(cancelled)?;
            if field_kind == FullWaveCompactFieldKind::IncidentTotal && source_m == observer_m {
                return Err(invalid(
                    "total layered compact array has a singular direct self-field",
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
                let monopole = PreparedLayeredFullWaveMonopole::new(
                    layer,
                    frequency_hz,
                    source_m,
                    observer_m,
                    monopole_limits,
                    cancelled,
                )?;
                set_column(&mut field, 0, monopole.evaluate_harmonic_line(C64::ONE)?);
                *monopole.diagnostics()
            } else {
                let monopole = PreparedLayeredReflectedFullWaveMonopole::new(
                    layer,
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
                let force = PreparedLayeredFullWaveForce::new(
                    layer,
                    frequency_hz,
                    source_m,
                    observer_m,
                    force_limits,
                    cancelled,
                )?;
                for column in 0..3 {
                    let unit =
                        core::array::from_fn(|i| if i == column { C64::ONE } else { C64::ZERO });
                    set_column(&mut field, column + 1, force.evaluate_harmonic_line(unit)?);
                }
                *force.diagnostics()
            } else {
                let force = PreparedLayeredReflectedFullWaveForce::new(
                    layer,
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
                diagnostics: FullWaveLayerCompactSourceDiagnostics {
                    monopole: monopole_diagnostics,
                    force: force_diagnostics,
                },
            });
        }
        poll(cancelled)?;
        Ok(Self {
            layer,
            field_kind,
            frequency_hz,
            observer_m,
            transfers,
            evaluations,
        })
    }

    pub const fn layer(&self) -> PlanarFluidLayer {
        self.layer
    }
    pub const fn field_kind(&self) -> FullWaveCompactFieldKind {
        self.field_kind
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
    pub fn source_diagnostics(
        &self,
        index: usize,
    ) -> Option<&FullWaveLayerCompactSourceDiagnostics> {
        self.transfers.get(index).map(|t| &t.diagnostics)
    }
    pub const fn scope(&self) -> &'static str {
        match self.field_kind {
            FullWaveCompactFieldKind::Transmitted => FULL_WAVE_LAYER_ARRAY_SCOPE,
            FullWaveCompactFieldKind::Reflected => FULL_WAVE_LAYER_REFLECTED_ARRAY_SCOPE,
            FullWaveCompactFieldKind::IncidentTotal => FULL_WAVE_LAYER_TOTAL_ARRAY_SCOPE,
        }
    }

    /// All Q/F and cross-source contributions combine before energy flux.
    pub fn evaluate(
        &self,
        sources: &[CompactSourcePhasor],
    ) -> Result<CoherentCompactField, AeroacError> {
        self.view().evaluate(sources)
    }

    /// Original peak-source clock with all layer/exterior phase already present.
    pub fn evaluate_harmonic_line(
        &self,
        sources: &[CompactSourcePhasor],
    ) -> Result<HarmonicFieldLine<3>, AeroacError> {
        self.view().evaluate_harmonic_line(sources)
    }

    /// Explicit low-rank source factors with their original physical phases.
    /// Moments add only across independent factors, retaining all internal
    /// reflection interference. Peak-amplitude factors describe line powers,
    /// not PSD/Hz; a source PSD in a bin of width df requires factors of 2 df S.
    pub fn evaluate_coherence_modes(
        &self,
        modes: &[&[CompactSourcePhasor]],
        budget: CoherenceBudget,
    ) -> Result<CompactFieldStatistics, AeroacError> {
        self.view().evaluate_coherence_modes(modes, budget)
    }

    /// Independent original physical sources; Q/F and layer histories remain
    /// coherent within each source. Independence is an explicit caller model.
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
