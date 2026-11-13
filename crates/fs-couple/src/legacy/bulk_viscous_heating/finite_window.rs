//! Finite-time heat from the original coherent Newtonian strain-rate field.
//!
//! For peak phasors on one authored clock, the physical source is
//! `eta [Re(sum_j du_j/dx exp(-i omega_j t))]^2`. Both difference-frequency
//! beats and carrier-sum terms enter the actual native thermal interval.
//! No cycle-mean factor is applied after squaring the real field, and no
//! finite-window acoustic boundary-flux difference is treated as heat.
//!
//! Positive composite Gauss integration resolves the original spatial fields
//! during preparation. Its explicit phase/attenuation resolution and analytic
//! single-line error check are separate from temporal integration. All spatial
//! products are contracted in double-double arithmetic into the existing
//! `fs_fft` harmonic spectrum. Thermal steps integrate that retained spectrum
//! analytically, without acoustic queries, spatial sampling, or time sampling.
//! Exact equal-frequency fields coalesce before products; distinct frequencies
//! never merge by tolerance. Negative unresolved mean heat refuses, without
//! clipping or an inferred incoherent approximation.
//!
//! The supplied responses must describe one original equilibrium and cylinder;
//! same-frequency responses must also have the same unit-drive load solution.
//! Frequency-dependent loads remain caller-authored. Every phasor uses the
//! SAME absolute seconds as the thermal timeline, `Re(A exp(-i omega t))`.
//! The nonempty clock identifier declares that correspondence, not a phase
//! alignment algorithm. These already-established harmonic fields do not
//! supply causal startup, material feedback, molecular heat, or mean flow.

use super::{
    BulkHeatCell, BulkTractionAbsorptionResponse, BulkViscousHeating, BulkViscousHeatingError,
    CellId, Dd, DepositionError, DepositionObserver, DepositionSimulation, EnergyAudit, Event,
    StepConfig,
};
use fs_fft::harmonic_window::{
    HarmonicWindow, PreparedHarmonicProductSpectrum, ProductSpectrumError, WindowError,
};
use std::collections::BTreeMap;

#[path = "finite_window/prepare.rs"]
mod prepare;

/// Preparation ceilings and independent spatial resolution controls.
#[derive(Clone, Copy, Debug)]
pub struct BulkViscousWindowLimits {
    pub maximum_lines: usize,
    pub maximum_cells: usize,
    /// Total positive Gauss nodes across the complete thermal partition.
    pub maximum_spatial_points: usize,
    /// Nodes times frequency-pair contractions, before spectral coalescing.
    pub maximum_pair_samples: usize,
    /// Logical retained/scratch vector payload and source-map entries,
    /// including the original cell binding; allocator/node overhead is excluded.
    pub maximum_workspace_bytes: usize,
    /// Each panel obeys max_j |k_j| * width <= this value, in (0,4].
    /// Refine this independently of the native thermal timestep.
    pub maximum_panel_exponent: f64,
    /// Per original line/cell, quadrature must match the original analytic
    /// cycle-mean viscous integral to this relative tolerance, in (0,1).
    /// This is a diagonal check, not a relative guarantee near coherent nulls.
    pub maximum_diagonal_relative_error: f64,
}
impl Default for BulkViscousWindowLimits {
    fn default() -> Self {
        Self {
            maximum_lines: 32,
            maximum_cells: 256,
            maximum_spatial_points: 65_536,
            maximum_pair_samples: 4_194_304,
            maximum_workspace_bytes: 16 * 1024 * 1024,
            maximum_panel_exponent: 1.0,
            maximum_diagonal_relative_error: 1e-8,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum BulkViscousWindowError {
    Heating(BulkViscousHeatingError),
    Spectrum(ProductSpectrumError),
    SpatialResolution {
        cell: CellId,
        line: usize,
        relative_error: f64,
        maximum: f64,
    },
    /// A signed roundoff residual is never silently deposited or clipped.
    NegativeHeat {
        cell: CellId,
        mean_power_w: f64,
    },
}
impl core::fmt::Display for BulkViscousWindowError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite-window bulk viscous heating refused: {self:?}")
    }
}
impl std::error::Error for BulkViscousWindowError {}
impl From<BulkViscousHeatingError> for BulkViscousWindowError {
    fn from(e: BulkViscousHeatingError) -> Self {
        Self::Heating(e)
    }
}
impl From<DepositionError> for BulkViscousWindowError {
    fn from(e: DepositionError) -> Self {
        Self::Heating(e.into())
    }
}
impl From<fs_duct::fluid::FluidDuctError> for BulkViscousWindowError {
    fn from(e: fs_duct::fluid::FluidDuctError) -> Self {
        Self::Heating(e.into())
    }
}
impl From<fs_material::MaterialError> for BulkViscousWindowError {
    fn from(e: fs_material::MaterialError) -> Self {
        Self::from(fs_duct::fluid::FluidDuctError::from(e))
    }
}
impl From<ProductSpectrumError> for BulkViscousWindowError {
    fn from(e: ProductSpectrumError) -> Self {
        Self::Spectrum(e)
    }
}
impl From<WindowError> for BulkViscousWindowError {
    fn from(e: WindowError) -> Self {
        Self::Spectrum(e.into())
    }
}
fn invalid(reason: &'static str) -> BulkViscousWindowError {
    BulkViscousHeatingError::Invalid(reason).into()
}
fn budget() -> BulkViscousWindowError {
    DepositionError::Budget("finite-window bulk viscous heat preparation").into()
}

/// Spatially integrated physical watts on the declared source clock.
#[derive(Clone, Debug)]
pub struct BulkViscousCellSpectrum {
    pub cell: CellId,
    pub start_distance_m: f64,
    pub end_distance_m: f64,
    pub spatial_points: usize,
    /// Worst analytic single-line comparison in this cell.
    pub maximum_diagonal_relative_error: f64,
    spectrum: PreparedHarmonicProductSpectrum,
}
impl BulkViscousCellSpectrum {
    /// Original native sum/difference spectrum, in watts; its DC alone is not
    /// the finite-window heat. Numerical null limitations remain explicit.
    pub fn spectrum(&self) -> &PreparedHarmonicProductSpectrum {
        &self.spectrum
    }

    /// Uniform physical-interval mean, already spatially integrated [W].
    /// The temporal owner retains both original start and duration separately.
    pub fn mean_heat_w(
        &self,
        start_s: f64,
        duration_s: f64,
    ) -> Result<f64, BulkViscousWindowError> {
        self.mean(HarmonicWindow::new(start_s, duration_s)?, &mut || Ok(()))
    }
    fn mean(
        &self,
        window: HarmonicWindow,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<f64, BulkViscousWindowError> {
        let value = spectral_call(checkpoint, |cancelled| {
            self.spectrum.mean_cancellable(window, cancelled)
        })?;
        if value < 0.0 {
            Err(BulkViscousWindowError::NegativeHeat {
                cell: self.cell,
                mean_power_w: value,
            })
        } else {
            Ok(value)
        }
    }
}

/// Immutable source spectra bound to original fixed-volume caloric cells.
pub struct BulkViscousWindowHeating<'a> {
    binding: BulkViscousHeating<'a>,
    responses: &'a [BulkTractionAbsorptionResponse],
    source_clock: &'a str,
    cells: Vec<BulkViscousCellSpectrum>,
}
impl<'a> BulkViscousWindowHeating<'a> {
    /// Prepare physical coherent products once. Source amplitudes are peak
    /// volume-flow phasors already retained by each original loaded response.
    /// All original thermal binding and frozen-temperature rules still apply.
    pub fn try_new(
        responses: &'a [BulkTractionAbsorptionResponse],
        source_clock: &'a str,
        thermal: &DepositionSimulation<'_>,
        partition: &[BulkHeatCell],
        limits: BulkViscousWindowLimits,
        maximum_temperature_change_k: f64,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, BulkViscousWindowError> {
        checkpoint()?;
        prepare::admit(responses, source_clock, partition.len(), limits)?;
        let binding = BulkViscousHeating::try_new(
            &responses[0],
            thermal,
            partition,
            limits.maximum_cells,
            maximum_temperature_change_k,
            &mut checkpoint,
        )?;
        let cells = prepare::cells(responses, &binding, limits, &mut checkpoint)?;
        checkpoint()?;
        Ok(Self {
            binding,
            responses,
            source_clock,
            cells,
        })
    }
    pub const fn responses(&self) -> &'a [BulkTractionAbsorptionResponse] {
        self.responses
    }
    pub const fn source_clock(&self) -> &'a str {
        self.source_clock
    }
    pub fn cells(&self) -> &[BulkViscousCellSpectrum] {
        &self.cells
    }

    /// Integrate this established coherent field on EACH native step's exact
    /// interval. Interval-average sources are passed once to the original
    /// thermal owner; source discontinuities require a newly prepared source
    /// at an authored advance boundary. There is no sampled acoustic timestep.
    pub fn advance_to<'t>(
        &self,
        initial: &DepositionSimulation<'t>,
        target_s: f64,
        config: StepConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<BulkViscousWindowStep<'t>, BulkViscousWindowError> {
        checkpoint()?;
        self.binding.validate(initial, &mut checkpoint)?;
        let observer = WindowHeat {
            owner: self,
            energy: vec![Dd::ZERO; self.cells.len()],
            pending: vec![Dd::ZERO; self.cells.len()],
        };
        let (thermal, observer, heat_report) =
            initial.advance_observed(&observer, target_s, config, &mut checkpoint)?;
        self.binding.validate(&thermal, &mut checkpoint)?;
        let mut total = Dd::ZERO;
        let mut cell_energy_j = Vec::with_capacity(self.cells.len());
        for (cell, energy) in self.cells.iter().zip(observer.energy) {
            checkpoint()?;
            total = total + energy;
            cell_energy_j.push((cell.cell, super::finite(energy.to_f64())?));
        }
        Ok(BulkViscousWindowStep {
            thermal,
            heat_report,
            cell_energy_j,
            viscous_heat_j: super::finite(total.to_f64())?,
        })
    }
}

/// Actual native interval heat, distinct from acoustic work or storage change.
pub struct BulkViscousWindowStep<'t> {
    pub thermal: DepositionSimulation<'t>,
    /// Includes these viscous joules once, plus original thermal exchanges.
    pub heat_report: EnergyAudit,
    pub viscous_heat_j: f64,
    pub cell_energy_j: Vec<(CellId, f64)>,
}

// Preserve a native checkpoint's actual refusal, not only a boolean label.
fn spectral_call<T>(
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    call: impl FnOnce(&mut dyn FnMut() -> bool) -> Result<T, ProductSpectrumError>,
) -> Result<T, BulkViscousWindowError> {
    let mut refusal = None;
    let result = call(&mut || match checkpoint() {
        Ok(()) => false,
        Err(error) => {
            refusal = Some(error);
            true
        }
    });
    if let Some(error) = refusal {
        Err(error.into())
    } else {
        Ok(result?)
    }
}

#[derive(Clone)]
struct WindowHeat<'h, 'a> {
    owner: &'h BulkViscousWindowHeating<'a>,
    energy: Vec<Dd>,
    pending: Vec<Dd>,
}
impl DepositionObserver for WindowHeat<'_, '_> {
    type Error = BulkViscousWindowError;
    fn requires_step_interval(&self) -> bool {
        true
    }
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        Ok(())
    }
    fn step_sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        start: f64,
        dt: f64,
        _: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.owner.binding.validate(thermal, checkpoint)?;
        let window = HarmonicWindow::new(start, dt)?;
        let mut sources = self.owner.binding.sources_w.clone();
        for (cell, energy) in self.owner.cells.iter().zip(&mut self.pending) {
            checkpoint()?;
            let power = cell.mean(window, checkpoint)?;
            *energy = super::energy(power, dt)?;
            *sources
                .get_mut(&cell.cell)
                .expect("original complete binding") = power;
        }
        Ok(sources)
    }
    fn step(
        &mut self,
        _: f64,
        _: f64,
        _: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.owner.binding.validate(thermal, checkpoint)?;
        for (sum, increment) in self.energy.iter_mut().zip(&self.pending) {
            *sum = *sum + *increment;
        }
        Ok(())
    }
}
