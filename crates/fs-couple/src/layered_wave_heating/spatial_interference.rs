//! Coherent finite-footprint heat on the original layered thermal cells.
//!
//! The existing `LayeredWaveHeating` stack already divides the optical films
//! at the native thermal cell faces and admits their original material cards
//! and fixed query coordinates. This adapter integrates the actual internal
//! field in EACH of those cells, over one declared continuous rectangle. Its
//! area must equal the original thermal body's planar area. It is therefore a
//! finite-volume, through-thickness thermal column of that footprint; lateral
//! temperature gradients, patch-edge scattering and lateral heat transfer are
//! not resolved. Each cell retains its native one-temperature cross section;
//! area equality does not establish the validity of that thermal approximation
//! under lateral fringes. Different patch positions can have different heat.
//!
//! Original optical coefficients and geometry stay frozen. Positive authored
//! depth quadrature must resolve attenuation and standing waves; neither a
//! material refresh nor a depth-convergence certificate is inferred. The native
//! thermal owner retains nonlinear caloric/conductivity laws, contacts, baths,
//! authored sources, event refusal and energy acceptance. Its existing material
//! domain applies; this adapter adds no phase-change or thermo-optical law.
//!
//! Analytic source joules use the accepted DD source clock and the ACTUAL native
//! thermal interval. A separate positive sum of cell primitive-ripple bounds
//! controls replacing that input by its mean during one implicit step. No beat
//! or carrier is removed on the basis of the thermal step size. Statistical
//! sources drive the deterministic thermal problem with expected heat; the
//! nonlinear result is not an ensemble-mean temperature.
//! These optical loss joules must not also be authored as a thermal source.

use super::LayeredWaveHeating;
use fs_em::plane_wave::angular_spectrum::{
    AngularSpectrumAperture, AngularSpectrumError, AngularSpectrumExposureEnsemble,
    AngularSpectrumExposureHeatChannel, AngularSpectrumExposureLayerLimits,
    AngularSpectrumHeatAreaQuadrature, AngularSpectrumHeatDepthSample,
    AngularSpectrumHeatVolumeQuadrature, AngularSpectrumIntegratedHeatLimits,
    PreparedAngularSpectrumExposureLayer,
};
use fs_fft::harmonic_window::{PreparedHarmonicProductSpectrum, ProductSpectrumError};
use fs_material::thermal::{
    conduction::LayerConductionConfig,
    transport::stack::transient::implicit::{
        LayeredThermalSolveControls, LayeredThermalSolveError, LayeredThermalSolveReport,
        LayeredThermalStep, LayeredThermalWorkspace,
    },
};
use fs_math::dd::Dd;
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation};

/// Each native cell has its own original positive depth samples, in metres
/// relative to that cell's optical entrance face. The complete cell thickness
/// is integrated; weights are never normalized or copied from another cell.
#[derive(Clone, Copy, Debug)]
pub struct SpatialInterferenceCellQuadrature<'q> {
    pub depth_samples: &'q [AngularSpectrumHeatDepthSample],
}

#[derive(Clone, Copy, Debug)]
pub struct SpatialInterferenceHeatingLimits {
    pub maximum_cells: usize,
    /// Conservative sum of each cell's complete native heat-preparation
    /// workspace, plus retained source descriptors and interval vectors. This
    /// deliberately also counts scratch released between cell preparations.
    /// Borrowed inputs, source payloads and allocator overhead are excluded.
    pub maximum_preparation_workspace_bytes: usize,
    /// Original per-cell field/material cache admission.
    pub layer: AngularSpectrumExposureLayerLimits,
    /// Native work bounds apply independently to every original cell.
    pub cell_heat: AngularSpectrumIntegratedHeatLimits,
}

/// Temporal approximation controls, independent of the native nonlinear and
/// physical-energy residual tolerances. A refusal asks for a shorter interval.
#[derive(Clone, Copy, Debug)]
pub struct SpatialInterferenceHeatingControls {
    pub maximum_thermal_step_s: f64,
    pub maximum_source_ripple_j: f64,
}

#[derive(Debug)]
pub enum SpatialInterferenceHeatingError {
    Invalid(&'static str),
    Electromagnetic(AngularSpectrumError),
    Spectrum(ProductSpectrumError),
    Thermal(DepositionError),
    Solve(LayeredThermalSolveError),
    SourceRipple { bound_j: f64, tolerance_j: f64 },
    Cancelled,
}
impl core::fmt::Display for SpatialInterferenceHeatingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "spatial interference heating: {self:?}")
    }
}
impl std::error::Error for SpatialInterferenceHeatingError {}
impl From<AngularSpectrumError> for SpatialInterferenceHeatingError {
    fn from(value: AngularSpectrumError) -> Self {
        Self::Electromagnetic(value)
    }
}
impl From<ProductSpectrumError> for SpatialInterferenceHeatingError {
    fn from(value: ProductSpectrumError) -> Self {
        Self::Spectrum(value)
    }
}
impl From<DepositionError> for SpatialInterferenceHeatingError {
    fn from(value: DepositionError) -> Self {
        Self::Thermal(value)
    }
}
impl From<LayeredThermalSolveError> for SpatialInterferenceHeatingError {
    fn from(value: LayeredThermalSolveError) -> Self {
        Self::Solve(value)
    }
}

/// Source-qualified immutable cell watt spectra. Retaining original prepared
/// layers and quadrature keeps both the field statistics and physical support
/// inspectable after their scalar spectra have been prepared.
pub struct SpatialInterferenceHeating<'h, 'b, 'a, 'm> {
    owner: &'h LayeredWaveHeating<'b, 'a, 'm>,
    ensemble: &'h AngularSpectrumExposureEnsemble<'h, 'h, 'h, 'm>,
    layers: Vec<PreparedAngularSpectrumExposureLayer<'h, 'h, 'h, 'm>>,
    quadrature: &'h [SpatialInterferenceCellQuadrature<'h>],
    footprint: AngularSpectrumAperture,
    source_clock: &'h str,
    spectra: Vec<PreparedHarmonicProductSpectrum>,
    preparation_workspace_bytes: usize,
}

/// An accepted native timeline and its independent phase epoch. The clock
/// label binds all prescribed lines, factors and successive steps. Fields are
/// private so continuation cannot silently reset the phase after a heat step.
pub struct SpatialInterferenceState<'c, 't> {
    thermal: DepositionSimulation<'t>,
    source_clock: &'c str,
    source_time_s: Dd,
}
impl<'t> SpatialInterferenceState<'_, 't> {
    pub fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub const fn source_time_s(&self) -> Dd {
        self.source_time_s
    }
    pub const fn source_clock(&self) -> &str {
        self.source_clock
    }
    pub fn into_thermal(self) -> DepositionSimulation<'t> {
        self.thermal
    }
}

/// Published only after all cell sources and the native thermal solve succeed.
pub struct SpatialInterferenceStep<'c, 't> {
    pub state: SpatialInterferenceState<'c, 't>,
    pub heat_report: LayeredThermalSolveReport,
    pub source_start_s: Dd,
    pub duration_s: f64,
    /// Sum of the original analytic cell integrals before f64 thermal input.
    pub absorbed_energy_j: Dd,
    /// Joules supplied once to the original native cell-energy seam, in its
    /// CellId order. Native energy auditing remains authoritative for rounding
    /// in the subsequent energy-to-mean-power conversion and residual solve.
    pub cell_energy_j: Vec<(CellId, f64)>,
    pub source_ripple_bound_j: f64,
}

impl<'h, 'b, 'a, 'm> SpatialInterferenceHeating<'h, 'b, 'a, 'm> {
    /// Admit the exact original cell mapping before any field reconstruction.
    /// Prepare every cell from ONE ensemble anchored on owner.optical_layers(),
    /// including when its other lines have different frequencies. The original
    /// exposure owner admits those lines' common sources/geometry/statistics;
    /// adjacent thermal cells cannot accidentally receive incompatible drives.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        owner: &'h LayeredWaveHeating<'b, 'a, 'm>,
        ensemble: &'h AngularSpectrumExposureEnsemble<'h, 'h, 'h, 'm>,
        quadrature: &'h [SpatialInterferenceCellQuadrature<'h>],
        footprint: AngularSpectrumAperture,
        source_clock: &'h str,
        limits: SpatialInterferenceHeatingLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, SpatialInterferenceHeatingError> {
        poll(&mut cancelled)?;
        let n = owner.cell_ids().len();
        if n == 0 || n > limits.maximum_cells || quadrature.len() != n {
            return Err(bad(
                "one depth quadrature per native thermal cell required within budget",
            ));
        }
        if source_clock.trim().is_empty() {
            return Err(bad("all source harmonics require one named physical clock"));
        }
        let area = footprint
            .extent_m()
            .into_iter()
            .fold(Dd::ONE, |a, x| a * Dd::from_f64(x));
        match_measure(area, Dd::from_f64(owner.area_m2()))?;
        let metadata = n
            .checked_mul(
                core::mem::size_of::<PreparedHarmonicProductSpectrum>()
                    + core::mem::size_of::<PreparedAngularSpectrumExposureLayer<'_, '_, '_, '_>>()
                    + core::mem::size_of::<(CellId, f64)>()
                    + 2 * core::mem::size_of::<f64>(),
            )
            .ok_or(bad("cell source storage overflows"))?;
        if metadata > limits.maximum_preparation_workspace_bytes {
            return Err(bad("cell source workspace exceeds its budget"));
        }
        let mut spectra = reserve(n)?;
        let mut layers = reserve(n)?;
        let mut workspace = metadata;
        for (index, rule) in quadrature.iter().enumerate() {
            poll(&mut cancelled)?;
            let film = owner
                .optical_layers()
                .layer(index)
                .ok_or(bad("original cell film is absent"))?;
            let remaining = limits.maximum_preparation_workspace_bytes - workspace;
            let mut layer_limits = limits.layer;
            layer_limits.maximum_workspace_bytes =
                layer_limits.maximum_workspace_bytes.min(remaining);
            let layer = ensemble.prepare_layer(index, layer_limits, &mut cancelled)?;
            if !core::ptr::eq(layer.source(), owner.optical_layers()) {
                return Err(bad(
                    "the ensemble must retain the exact original cell-split optical anchor",
                ));
            }
            let mut cell_limits = limits.cell_heat;
            cell_limits.maximum_workspace_bytes =
                cell_limits.maximum_workspace_bytes.min(remaining);
            let load = layer.integrated_heat_spectrum(
                AngularSpectrumHeatVolumeQuadrature {
                    depth_span_m: [0.0, film.thickness_m],
                    depth_samples: rule.depth_samples,
                    transverse: AngularSpectrumHeatAreaQuadrature::ContinuousRectangle(footprint),
                },
                AngularSpectrumExposureHeatChannel::Total,
                cell_limits,
                &mut cancelled,
            )?;
            // The field owner retains actual weights. Refuse mismatched volumes
            // instead of renormalizing heat or spreading it over a larger body.
            match_measure(
                load.represented_volume_m3_dd(),
                Dd::from_f64(owner.area_m2()) * Dd::from_f64(film.thickness_m),
            )?;
            workspace = workspace
                .checked_add(load.preparation_workspace_bytes())
                .filter(|&v| v <= limits.maximum_preparation_workspace_bytes)
                .ok_or(bad("cell source workspace exceeds its budget"))?;
            spectra.push(load.into_spectrum());
            layers.push(layer);
        }
        poll(&mut cancelled)?;
        Ok(Self {
            owner,
            ensemble,
            layers,
            quadrature,
            footprint,
            source_clock,
            spectra,
            preparation_workspace_bytes: workspace,
        })
    }

    pub fn original_owner(&self) -> &'h LayeredWaveHeating<'b, 'a, 'm> {
        self.owner
    }
    pub fn original_ensemble(&self) -> &'h AngularSpectrumExposureEnsemble<'h, 'h, 'h, 'm> {
        self.ensemble
    }
    pub fn original_layers(&self) -> &[PreparedAngularSpectrumExposureLayer<'h, 'h, 'h, 'm>] {
        &self.layers
    }
    pub fn quadrature(&self) -> &'h [SpatialInterferenceCellQuadrature<'h>] {
        self.quadrature
    }
    pub const fn footprint(&self) -> AngularSpectrumAperture {
        self.footprint
    }
    pub const fn source_clock(&self) -> &'h str {
        self.source_clock
    }
    /// Already in W. Never multiply these spectra by cell volume or area again.
    pub fn cell_spectra(&self) -> &[PreparedHarmonicProductSpectrum] {
        &self.spectra
    }
    pub const fn preparation_workspace_bytes(&self) -> usize {
        self.preparation_workspace_bytes
    }

    /// Explicitly declare the source epoch for an existing accepted thermal
    /// state. The epoch may be large or offset from native elapsed time; future
    /// calls advance it by the actual accepted native interval without reset.
    pub fn initial_state<'t>(
        &self,
        thermal: DepositionSimulation<'t>,
        source_time_s: Dd,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<SpatialInterferenceState<'h, 't>, SpatialInterferenceHeatingError> {
        poll(&mut cancelled)?;
        let source_time_s = checked_clock(source_time_s)?;
        self.owner
            .thermal_body()
            .cell_temperatures(&thermal, || checkpoint(&mut cancelled))?;
        poll(&mut cancelled)?;
        Ok(SpatialInterferenceState {
            thermal,
            source_clock: self.source_clock,
            source_time_s,
        })
    }

    /// One explicit, frozen-optics interval. The native backward-Euler step
    /// owns physical duration, full candidate construction and acceptance.
    /// No source or partial temperature is published on ripple/domain/solver
    /// refusal or cancellation. Caller scratch may change on refusal.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_to<'t>(
        &self,
        initial: &SpatialInterferenceState<'_, 't>,
        target_s: f64,
        temporal: SpatialInterferenceHeatingControls,
        conduction: LayerConductionConfig,
        solve: LayeredThermalSolveControls,
        workspace: &mut LayeredThermalWorkspace,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<SpatialInterferenceStep<'h, 't>, SpatialInterferenceHeatingError> {
        poll(&mut cancelled)?;
        if initial.source_clock != self.source_clock {
            return Err(bad(
                "accepted state and optical sources must share their declared clock",
            ));
        }
        let thermal_step = LayeredThermalStep::new(
            self.owner.thermal_body(),
            &initial.thermal,
            target_s,
            conduction,
            || checkpoint(&mut cancelled),
        )?;
        let dt = thermal_step.duration_s();
        let source = integrate_cells(
            &self.spectra,
            initial.source_time_s,
            dt,
            temporal,
            &mut cancelled,
        )?;
        let thermal_step = thermal_step
            .with_cell_energy_j(&source.cell_energy_j, || checkpoint(&mut cancelled))?;
        let (thermal, _, heat_report) =
            thermal_step.solve(solve, workspace, || checkpoint(&mut cancelled))?;
        let mut cell_energy_j = reserve(self.spectra.len())?;
        for (&id, &energy) in self.owner.cell_ids().iter().zip(&source.cell_energy_j) {
            poll(&mut cancelled)?;
            cell_energy_j.push((id, energy));
        }
        poll(&mut cancelled)?;
        Ok(SpatialInterferenceStep {
            state: SpatialInterferenceState {
                thermal,
                source_clock: self.source_clock,
                source_time_s: source.end_s,
            },
            heat_report,
            source_start_s: initial.source_time_s,
            duration_s: dt,
            absorbed_energy_j: source.total_j,
            cell_energy_j,
            source_ripple_bound_j: source.ripple_j,
        })
    }
}

struct CellIntegrals {
    end_s: Dd,
    total_j: Dd,
    cell_energy_j: Vec<f64>,
    ripple_j: f64,
}

/// The actual native interval is authoritative. Every integral is irreversible
/// reservoir heat; an unrepresentable/negative result refuses, never clips.
fn integrate_cells(
    spectra: &[PreparedHarmonicProductSpectrum],
    start_s: Dd,
    duration_s: f64,
    controls: SpatialInterferenceHeatingControls,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<CellIntegrals, SpatialInterferenceHeatingError> {
    poll(cancelled)?;
    if !duration_s.is_finite()
        || duration_s <= 0.0
        || !controls.maximum_thermal_step_s.is_finite()
        || controls.maximum_thermal_step_s <= 0.0
        || duration_s > controls.maximum_thermal_step_s
        || !controls.maximum_source_ripple_j.is_finite()
        || controls.maximum_source_ripple_j < 0.0
    {
        return Err(bad(
            "a positive interval within the thermal step cap and a finite nonnegative ripple limit are required",
        ));
    }
    let start_s = checked_clock(start_s)?;
    let end_s = checked_clock(start_s + Dd::from_f64(duration_s))?;
    let elapsed = end_s - start_s;
    if elapsed.hi < 0.0 || (elapsed.hi == 0.0 && elapsed.lo <= 0.0) {
        return Err(bad("interval does not advance the DD source clock"));
    }
    let mut cell_energy_j = reserve(spectra.len())?;
    let mut ripple_j = 0.0;
    let mut total_j = Dd::ZERO;
    for spectrum in spectra {
        poll(cancelled)?;
        let bound = spectrum.input_ripple_bound_cancellable(duration_s, &mut *cancelled)?;
        ripple_j = outward_sum(ripple_j, bound);
        if !ripple_j.is_finite() || ripple_j > controls.maximum_source_ripple_j {
            return Err(SpatialInterferenceHeatingError::SourceRipple {
                bound_j: ripple_j,
                tolerance_j: controls.maximum_source_ripple_j,
            });
        }
        let energy = spectrum
            .causal_integral_from_clock_dd_cancellable(start_s, duration_s, 0.0, &mut *cancelled)?
            .integral_dd();
        let value = energy.to_f64();
        if !finite(energy)
            || !value.is_finite()
            || energy.hi < 0.0
            || (energy.hi == 0.0 && energy.lo < 0.0)
            || (value == 0.0 && (energy.hi != 0.0 || energy.lo != 0.0))
        {
            return Err(bad(
                "irreversible cell heat lost finite nonnegative representation",
            ));
        }
        total_j = total_j + energy;
        if !finite(total_j) {
            return Err(bad("total irreversible cell heat overflows"));
        }
        cell_energy_j.push(value);
    }
    poll(cancelled)?;
    Ok(CellIntegrals {
        end_s,
        total_j,
        cell_energy_j,
        ripple_j,
    })
}

fn bad(message: &'static str) -> SpatialInterferenceHeatingError {
    SpatialInterferenceHeatingError::Invalid(message)
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), SpatialInterferenceHeatingError> {
    if cancelled() {
        Err(SpatialInterferenceHeatingError::Cancelled)
    } else {
        Ok(())
    }
}
fn checkpoint(cancelled: &mut impl FnMut() -> bool) -> Result<(), DepositionError> {
    if cancelled() {
        Err(DepositionError::Cancelled)
    } else {
        Ok(())
    }
}
fn finite(value: Dd) -> bool {
    value.hi.is_finite() && value.lo.is_finite()
}
fn checked_clock(value: Dd) -> Result<Dd, SpatialInterferenceHeatingError> {
    if !finite(value) {
        return Err(bad("source clock must be finite"));
    }
    let value = Dd::from_pair(value.hi, value.lo);
    if !finite(value) {
        return Err(bad("source clock cannot be normalized"));
    }
    Ok(value)
}
fn reserve<T>(n: usize) -> Result<Vec<T>, SpatialInterferenceHeatingError> {
    let mut result = Vec::new();
    result
        .try_reserve_exact(n)
        .map_err(|_| bad("cell source allocation failed"))?;
    Ok(result)
}
fn match_measure(actual: Dd, expected: Dd) -> Result<(), SpatialInterferenceHeatingError> {
    if !finite(actual)
        || !finite(expected)
        || actual.to_f64() <= 0.0
        || expected.to_f64() <= 0.0
        || ((actual - expected) / expected).to_f64().abs() > 64.0 * f64::EPSILON
    {
        return Err(bad(
            "optical footprint and cell volume must match the original planar thermal geometry",
        ));
    }
    Ok(())
}
// Positive L1 sum: sources in different cells cannot hide local forcing by
// cancelling their phases. Outward rounding preserves zero only for zero.
fn outward_sum(a: f64, b: f64) -> f64 {
    if a == 0.0 {
        b
    } else if b == 0.0 {
        a
    } else {
        (a + b).next_up()
    }
}

#[cfg(test)]
#[path = "spatial_interference/tests.rs"]
mod tests;
