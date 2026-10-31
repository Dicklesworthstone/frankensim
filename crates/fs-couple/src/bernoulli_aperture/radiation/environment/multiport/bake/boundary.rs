//! Material-body diffraction loading through the original shared-pole fitter.
//!
//! Original finite unit-flow traces retain their complete self/mutual reaction.
//! The boundary solve, rational identification and exact-ZOH terminal each keep
//! their own error and passivity limitations; no observation window is a load.

use super::{
    BaffledMimoBakeError, BaffledMimoBakeOptions, BaffledMimoBakeReport, MimoRadiationLimits,
    admit_common, fit_source, realization, scalar,
};
use fs_bem::helmholtz::{IterativeOptions, Medium};
use fs_bem::observation_transfers::{
    IterativeBoundaryObservation, IterativeObservationError, SourceId,
};
use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::CancelGate;
use fs_material::gas::GasState;
use fs_math::c64::C64;
use fs_vfit::mimo::load::{MimoLoadError, MimoLoadLimits, PreparedMimoLoad};
use std::cell::RefCell;

/// Complete source-sampling allowances, separate from fit and runtime memory.
#[derive(Clone, Copy, Debug)]
pub struct BoundaryMimoSamplingLimits {
    /// Both pressure and normal-velocity coefficients at one physical frequency.
    pub maximum_retained_boundary_entries: usize,
    /// Every original source pair at one frequency, without truncation.
    pub maximum_port_matrix_entries: usize,
    /// All source solves across training, withheld and causality grids.
    pub maximum_source_solves: usize,
    /// All source-pair/panel products across those same physical grids.
    pub maximum_projection_panel_products: u64,
    /// Original transpose defect allowance, at most one tenth of the fit error.
    /// The reciprocal fit still measures its error against the unmodified data.
    pub maximum_reciprocity_defect: f64,
}

/// Actual original-equation statistics, not an estimate of continuum/FMM error.
#[derive(Clone, Debug)]
pub struct BoundaryMimoSamplingReport {
    pub frequency_samples: usize,
    pub source_solves: usize,
    pub projection_panel_products: u64,
    pub maximum_gmres_relative_residual: f64,
    pub maximum_reciprocity_defect: f64,
    pub minimum_passivity_margin: f64,
    /// Hash of every sampled frequency and original complex physical matrix.
    pub sampled_operator_hash: ContentHash,
}

/// Source, units and approximation remain attached to the admitted runtime.
#[derive(Clone, Debug)]
pub struct BoundaryMimoRecord {
    sources: Vec<SourceId>,
    exterior: GasState,
    operator_hash: ContentHash,
    aperture_support_areas_m2: Vec<f64>,
    band_hz: [f64; 2],
    source_sample_band_hz: [f64; 2],
    time_step_s: f64,
    options: BaffledMimoBakeOptions,
    limits: MimoRadiationLimits,
    sampling: BoundaryMimoSamplingLimits,
    solver: IterativeOptions,
    report: BaffledMimoBakeReport,
    sampling_report: BoundaryMimoSamplingReport,
    declaration_identity: [u8; 32],
    realization_identity: [u8; 32],
}
impl BoundaryMimoRecord {
    pub fn sources(&self) -> &[SourceId] {
        &self.sources
    }
    pub const fn gas(&self) -> &GasState {
        &self.exterior
    }
    pub const fn operator_hash(&self) -> ContentHash {
        self.operator_hash
    }
    /// Original nonzero trace-support areas, without changing their profiles.
    pub fn aperture_support_areas_m2(&self) -> &[f64] {
        &self.aperture_support_areas_m2
    }
    /// Fitted/withheld band; sampled causality probes can extend below its floor.
    pub const fn band_hz(&self) -> [f64; 2] {
        self.band_hz
    }
    pub const fn source_sample_band_hz(&self) -> [f64; 2] {
        self.source_sample_band_hz
    }
    pub const fn time_step_s(&self) -> f64 {
        self.time_step_s
    }
    pub const fn options(&self) -> BaffledMimoBakeOptions {
        self.options
    }
    pub const fn limits(&self) -> MimoRadiationLimits {
        self.limits
    }
    pub const fn sampling_limits(&self) -> BoundaryMimoSamplingLimits {
        self.sampling
    }
    pub fn solver_options(&self) -> &IterativeOptions {
        &self.solver
    }
    pub const fn report(&self) -> &BaffledMimoBakeReport {
        &self.report
    }
    pub const fn sampling_report(&self) -> &BoundaryMimoSamplingReport {
        &self.sampling_report
    }
    pub const fn declaration_identity(&self) -> [u8; 32] {
        self.declaration_identity
    }
    pub const fn realization_identity(&self) -> [u8; 32] {
        self.realization_identity
    }
}

/// One quiescent original terminal load plus its immutable physical binding.
#[derive(Clone, Debug)]
pub struct BoundaryMimoLoad {
    load: PreparedMimoLoad,
    record: BoundaryMimoRecord,
}
impl BoundaryMimoLoad {
    pub const fn load(&self) -> &PreparedMimoLoad {
        &self.load
    }
    pub const fn record(&self) -> &BoundaryMimoRecord {
        &self.record
    }
    pub fn into_parts(self) -> (PreparedMimoLoad, BoundaryMimoRecord) {
        (self.load, self.record)
    }
}

/// A refused source, original solve, realization or explicit budget.
#[derive(Debug)]
pub enum BoundaryMimoBakeError {
    Invalid(&'static str),
    Source(IterativeObservationError),
    Sample {
        frequency_hz: f64,
        source: IterativeObservationError,
    },
    Work {
        kind: &'static str,
        required: u64,
        allowed: u64,
    },
    Fit(BaffledMimoBakeError),
    Runtime(MimoLoadError),
    Cancelled,
}
impl From<BaffledMimoBakeError> for BoundaryMimoBakeError {
    fn from(error: BaffledMimoBakeError) -> Self {
        match error {
            BaffledMimoBakeError::Cancelled => Self::Cancelled,
            other => Self::Fit(other),
        }
    }
}
impl core::fmt::Display for BoundaryMimoBakeError {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(formatter, "material boundary MIMO bake: {self:?}")
    }
}
impl std::error::Error for BoundaryMimoBakeError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Source(error) | Self::Sample { source: error, .. } => Some(error),
            Self::Fit(error) => Some(error),
            Self::Runtime(error) => Some(error),
            _ => None,
        }
    }
}
fn poll(gate: &CancelGate) -> Result<(), BoundaryMimoBakeError> {
    if gate.is_requested() {
        Err(BoundaryMimoBakeError::Cancelled)
    } else {
        Ok(())
    }
}
fn invalid(what: &'static str) -> BoundaryMimoBakeError {
    BoundaryMimoBakeError::Invalid(what)
}
fn work(kind: &'static str, required: u64, allowed: u64) -> Result<(), BoundaryMimoBakeError> {
    if required > allowed {
        Err(BoundaryMimoBakeError::Work {
            kind,
            required,
            allowed,
        })
    } else {
        Ok(())
    }
}

/// Build a proper rational terminal from the actual finite material-body solve.
///
/// Sources must be the original unit-volume-flow aperture coordinates. The
/// matrix already has pressure/flow units: only conjugation from exp(-iwt) to
/// s=+iw is performed. All obstacle and passive-patch scattering remains in it.
/// The full requested sampling grids and mechanical Nyquist band must fit the
/// original geometry/material domain before the first expensive source solve.
///
/// This reuses the existing three-order fit, source-data causality, withheld
/// error, bounded residue repair and continuous/mean-ZOH passivity owners. Raw
/// sample asymmetry remains in the fit error; raw sampled activity refuses.
/// Delay is approximated by shared poles, not separately extracted. GridOnly
/// stays GridOnly; this is neither a storage certificate nor the full W3b
/// closed-loop/continuum qualification. Below-band and DC source fidelity are
/// unresolved even though the runtime's sampled passivity grid includes DC.
#[allow(clippy::too_many_arguments)]
pub fn bake_boundary_mimo(
    observer: &IterativeBoundaryObservation<'_>,
    exterior: GasState,
    mechanics_rate_hz: u32,
    band_hz: [f64; 2],
    options: BaffledMimoBakeOptions,
    limits: MimoRadiationLimits,
    sampling: BoundaryMimoSamplingLimits,
    gate: &CancelGate,
) -> Result<BoundaryMimoLoad, BoundaryMimoBakeError> {
    poll(gate)?;
    let n = observer.sources().len();
    admit_common(&exterior, n, mechanics_rate_hz, band_hz, options, limits)?;
    if observer.medium() != Medium::from(&exterior) {
        return Err(invalid(
            "boundary and terminal exterior gas coefficients must match exactly",
        ));
    }
    if !sampling.maximum_reciprocity_defect.is_finite()
        || sampling.maximum_reciprocity_defect < 0.0
        || sampling.maximum_reciprocity_defect > 0.1 * options.maximum_relative_error
    {
        return Err(invalid(
            "raw reciprocity allowance must fit within a tenth of the total fit error",
        ));
    }
    let source_sample_band_hz = [
        band_hz[0].min(band_hz[1] / (2 * options.frequency_samples - 1) as f64),
        band_hz[1],
    ];
    observer
        .admit_loading_band_hz(source_sample_band_hz)
        .map_err(BoundaryMimoBakeError::Source)?;
    let panels = observer.surface_panel_count();
    let retained = n
        .checked_mul(panels)
        .and_then(|v| v.checked_mul(2))
        .ok_or_else(|| invalid("retained boundary shape overflow"))?;
    work(
        "retained boundary entries",
        retained as u64,
        sampling.maximum_retained_boundary_entries as u64,
    )?;
    work(
        "source matrix entries",
        (n * n) as u64,
        sampling.maximum_port_matrix_entries as u64,
    )?;
    let frequencies = 3 * options.frequency_samples - 1;
    work(
        "complete source solves",
        (n * frequencies) as u64,
        sampling.maximum_source_solves as u64,
    )?;
    let products_per_frequency = (n * n) as u64 * panels as u64;
    work(
        "complete source-pair panel products",
        products_per_frequency * frequencies as u64,
        sampling.maximum_projection_panel_products,
    )?;
    let aperture_support_areas_m2 = observer
        .aperture_support_areas_m2(gate)
        .map_err(BoundaryMimoBakeError::Source)?;
    let sampler = BoundarySampler {
        observer,
        sampling,
        products_per_frequency,
        statistics: RefCell::new(SamplingState {
            frequencies: 0,
            solves: 0,
            products: 0,
            residual: 0.0,
            reciprocity: 0.0,
            margin: f64::INFINITY,
            hash: DomainHasher::new("org.frankensim.fs-couple.material-boundary-samples.v1"),
        }),
    };
    let (model, report) = fit_source(
        n,
        mechanics_rate_hz,
        band_hz,
        options,
        limits,
        gate,
        |hz, gate| sampler.rows(hz, gate),
    )?;
    poll(gate)?;
    let time_step_s = 1.0 / f64::from(mechanics_rate_hz);
    let load = PreparedMimoLoad::new(
        &model,
        time_step_s,
        MimoLoadLimits {
            band_rad_s: [
                core::f64::consts::TAU * band_hz[0],
                core::f64::consts::TAU * band_hz[1],
            ],
            sampled_passivity_points: limits.sampled_passivity_points,
            max_states: limits.max_states,
            max_memory_bytes: limits.max_memory_bytes,
        },
    )
    .map_err(BoundaryMimoBakeError::Runtime)?;
    poll(gate)?;
    let sampling_report = sampler.into_report();
    let declaration_identity = declaration(
        observer,
        &exterior,
        mechanics_rate_hz,
        band_hz,
        limits,
        sampling,
        &sampling_report,
    );
    let realization_identity = realization(declaration_identity, &model, options);
    Ok(BoundaryMimoLoad {
        load,
        record: BoundaryMimoRecord {
            sources: observer.sources().to_vec(),
            exterior,
            operator_hash: observer.operator_hash(),
            aperture_support_areas_m2,
            band_hz,
            source_sample_band_hz,
            time_step_s,
            options,
            limits,
            sampling,
            solver: observer.solver_options().clone(),
            report,
            sampling_report,
            declaration_identity,
            realization_identity,
        },
    })
}

struct SamplingState {
    frequencies: usize,
    solves: usize,
    products: u64,
    residual: f64,
    reciprocity: f64,
    margin: f64,
    hash: DomainHasher,
}
struct BoundarySampler<'a, 's> {
    observer: &'a IterativeBoundaryObservation<'s>,
    sampling: BoundaryMimoSamplingLimits,
    products_per_frequency: u64,
    statistics: RefCell<SamplingState>,
}
impl BoundarySampler<'_, '_> {
    fn rows(
        &self,
        frequencies: &[f64],
        gate: &CancelGate,
    ) -> Result<Vec<Vec<C64>>, BoundaryMimoBakeError> {
        let mut rows = Vec::with_capacity(frequencies.len());
        for &frequency_hz in frequencies {
            poll(gate)?;
            let n = self.observer.sources().len();
            {
                let state = self.statistics.borrow();
                work(
                    "complete source solves",
                    (state.solves + n) as u64,
                    self.sampling.maximum_source_solves as u64,
                )?;
                work(
                    "complete source-pair panel products",
                    state.products + self.products_per_frequency,
                    self.sampling.maximum_projection_panel_products,
                )?;
            }
            let source_error = |source| BoundaryMimoBakeError::Sample {
                frequency_hz,
                source,
            };
            let response = self
                .observer
                .solve_frequency(
                    core::f64::consts::TAU * frequency_hz,
                    self.sampling.maximum_retained_boundary_entries,
                    gate,
                )
                .map_err(source_error)?;
            let matrix = response
                .velocity_port_impedance(
                    self.sampling.maximum_port_matrix_entries,
                    self.products_per_frequency,
                    gate,
                )
                .map_err(source_error)?;
            let original = matrix.acoustic_impedance().map_err(source_error)?;
            let margin = matrix
                .passivity_margin()
                .ok_or_else(|| invalid("original boundary passivity diagnostic is unresolved"))?;
            if matrix.reciprocity_defect() > self.sampling.maximum_reciprocity_defect
                || margin < -1e-12
            {
                return Err(invalid(
                    "original boundary matrix fails declared reciprocity or sampled passivity",
                ));
            }
            let mut state = self.statistics.borrow_mut();
            state.frequencies += 1;
            state.solves += n;
            state.products += self.products_per_frequency;
            state.reciprocity = state.reciprocity.max(matrix.reciprocity_defect());
            state.margin = state.margin.min(margin);
            for report in response.reports() {
                state.residual = state.residual.max(report.solve.rel_residual);
            }
            scalar(&mut state.hash, frequency_hz);
            for value in original {
                scalar(&mut state.hash, value.re);
                scalar(&mut state.hash, value.im);
            }
            // Unit-flow normalization belongs to the original boundary basis.
            // Preserve all matrix entries and conjugate exactly once.
            rows.push(original.iter().map(|value| value.conj()).collect());
        }
        poll(gate)?;
        Ok(rows)
    }
    fn into_report(self) -> BoundaryMimoSamplingReport {
        let state = self.statistics.into_inner();
        BoundaryMimoSamplingReport {
            frequency_samples: state.frequencies,
            source_solves: state.solves,
            projection_panel_products: state.products,
            maximum_gmres_relative_residual: state.residual,
            maximum_reciprocity_defect: state.reciprocity,
            minimum_passivity_margin: state.margin,
            sampled_operator_hash: state.hash.finalize(),
        }
    }
}

fn declaration(
    observer: &IterativeBoundaryObservation<'_>,
    gas: &GasState,
    rate: u32,
    band: [f64; 2],
    limits: MimoRadiationLimits,
    sampling: BoundaryMimoSamplingLimits,
    report: &BoundaryMimoSamplingReport,
) -> [u8; 32] {
    let mut hash = DomainHasher::new("org.frankensim.fs-couple.material-boundary-source.v1");
    hash.update(observer.operator_hash().as_bytes());
    hash.update(report.sampled_operator_hash.as_bytes());
    for value in [
        gas.temperature,
        gas.pressure,
        gas.density,
        gas.sound_speed,
        gas.dynamic_viscosity,
        gas.thermal_conductivity,
        gas.gamma,
        gas.specific_gas_constant,
        gas.specific_heat_cp,
        gas.prandtl,
        gas.characteristic_impedance,
        gas.water_mole_fraction,
        band[0],
        band[1],
        sampling.maximum_reciprocity_defect,
    ] {
        scalar(&mut hash, value);
    }
    for value in [
        u64::from(rate),
        limits.max_states as u64,
        limits.max_memory_bytes as u64,
        limits.sampled_passivity_points as u64,
        sampling.maximum_retained_boundary_entries as u64,
        sampling.maximum_port_matrix_entries as u64,
        sampling.maximum_source_solves as u64,
        sampling.maximum_projection_panel_products,
    ] {
        hash.update(&value.to_le_bytes());
    }
    hash.finalize().0
}

#[cfg(test)]
pub(crate) mod tests;
