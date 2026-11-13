//! Finite-time coherent molecular heat on the original gas inventory.
//!
//! The original molecular law supplies its exact positive compression arms:
//! `p = K_eq epsilon + sum r_j`, `r'_j + r_j/tau_j = delta_K_j epsilon'`.
//! Their instantaneous irreversible heat is `sum r_j(t)^2/(delta_K_j tau_j)`.
//! Each harmonic strain is the ORIGINAL BEM pressure divided by its molecular
//! modulus. Arm stresses retain their own phases; taking a square root of the
//! scalar cycle-mean attenuation would discard this information.
//!
//! Positive original volume weights contract these real squares into the
//! existing exact harmonic spectrum. Both carrier sums and frequency beats
//! survive. Repeated interval queries perform no acoustic solve, field query,
//! spatial quadrature, or time sampling. Reservoir acoustic storage is not
//! deposited as heat, and finite-window boundary-work loss is not used.
//!
//! This is an already-established coherent harmonic field, with all source
//! phasors using `Re(A exp(-i omega t))` and one declared absolute-second clock.
//! It does not infer startup states, random-phase statistics or a new material
//! feedback solve. Original fixed-volume caloric storage, source-domain checks
//! and frozen-field temperature admission remain owned by molecular_heating.

use super::{
    MolecularHeatingError, finite_body_environment, heated_gas, poll,
    spectrum::MolecularGasSpectrumHeating,
};
use crate::{
    cavity_phs::piston::thermal::HeatExchangeSolve, finite_gas_thermal::FiniteGasThermalEnvironment,
};
use fs_bem::observation_transfers::{MolecularObservationError, TransferError};
use fs_conduction::lumped::LumpedThermalEnvironment;
use fs_fft::harmonic_window::{
    HarmonicWindow, PreparedHarmonicProductSpectrum, ProductSpectrumError, WindowError,
};
use fs_material::{
    MaterialError,
    relaxation::RelaxationState,
    thermochemical_gas::energy::{GasEnergySolve, GasHeatStep, IsochoricGasParcel},
};
use fs_math::c64::C64;

mod prepare;

/// Additional bounds for fixed coherent-field preparation.
#[derive(Clone, Copy, Debug)]
pub struct MolecularWindowLimits {
    pub maximum_bins: usize,
    pub maximum_receivers: usize,
    /// Conservative field and reservoir frequency-pair product count.
    pub maximum_products: usize,
    /// Logical vector/clock payload, including retained absorption workspaces.
    /// Borrowed response banks, allocator overhead and the original bounded
    /// <=32-reservoir cold material lowering are excluded.
    pub maximum_workspace_bytes: usize,
    /// Original caloric loss versus positive-arm loss, at each frequency.
    /// This is a constitutive check, not a spatial quadrature certificate.
    pub maximum_constitutive_relative_error: f64,
}
impl Default for MolecularWindowLimits {
    fn default() -> Self {
        Self {
            maximum_bins: 32,
            maximum_receivers: 4096,
            maximum_products: 4_194_304,
            maximum_workspace_bytes: 16 * 1024 * 1024,
            maximum_constitutive_relative_error: 1e-8,
        }
    }
}

#[derive(Debug)]
pub enum MolecularWindowError {
    Heating(MolecularHeatingError),
    Spectrum(ProductSpectrumError),
    Invalid(&'static str),
    ConstitutiveLoss {
        bin: usize,
        relative_error: f64,
        maximum: f64,
    },
    /// Unresolved signed roundoff is refused, never clipped or deposited.
    NegativeHeat {
        mean_power_w: f64,
    },
}
impl core::fmt::Display for MolecularWindowError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite-window molecular heating refused: {self:?}")
    }
}
impl std::error::Error for MolecularWindowError {}
impl From<MolecularHeatingError> for MolecularWindowError {
    fn from(e: MolecularHeatingError) -> Self {
        Self::Heating(e)
    }
}
impl From<MaterialError> for MolecularWindowError {
    fn from(e: MaterialError) -> Self {
        Self::Heating(e.into())
    }
}
impl From<TransferError> for MolecularWindowError {
    fn from(e: TransferError) -> Self {
        Self::Heating(MolecularObservationError::from(e).into())
    }
}
impl From<ProductSpectrumError> for MolecularWindowError {
    fn from(e: ProductSpectrumError) -> Self {
        Self::Spectrum(e)
    }
}
impl From<WindowError> for MolecularWindowError {
    fn from(e: WindowError) -> Self {
        Self::Spectrum(e.into())
    }
}

/// Original source-bound coherent molecular heating, in physical watts.
///
/// Distinct frequencies share exactly the same receiver identities/positions
/// and original volume quadrature. Sources combine coherently WITHIN each
/// bank before arm squares; duplicate frequency banks are already refused by
/// MolecularGasSpectrumHeating. All sources at one frequency belong in that
/// original bank. Volume equality alone cannot establish spatial overlap.
///
/// The exact spectral contraction uses DD sums. Extremely small coherent
/// nulls may still be unresolved after contraction; a negative interval mean
/// refuses instead of inventing heat or discarding cross terms.
pub struct MolecularGasWindowHeating<'a> {
    binding: MolecularGasSpectrumHeating<'a>,
    peak_rates: Vec<Vec<C64>>,
    source_clock: String,
    spectrum: PreparedHarmonicProductSpectrum,
    workspace_bytes: usize,
}

impl<'a> MolecularGasSpectrumHeating<'a> {
    /// Freeze original per-bank coherent peak rates on one physical clock.
    /// Preparation is atomic and cancellable; failure never changes the gas.
    /// The caller retains responsibility for BEM/spatial convergence and for
    /// using already-established fields throughout every requested interval.
    pub fn prepare_coherent_window(
        self,
        peak_rates: &[&[C64]],
        source_clock: &str,
        limits: MolecularWindowLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularGasWindowHeating<'a>, MolecularWindowError> {
        poll(&mut cancelled)?;
        let (spectrum, workspace_bytes) =
            prepare::spectrum(&self, peak_rates, source_clock, limits, &mut cancelled)?;
        let rates = peak_rates.iter().map(|row| row.to_vec()).collect();
        poll(&mut cancelled)?;
        Ok(MolecularGasWindowHeating {
            binding: self,
            peak_rates: rates,
            source_clock: source_clock.to_owned(),
            spectrum,
            workspace_bytes,
        })
    }
}

impl<'a> MolecularGasWindowHeating<'a> {
    pub const fn binding(&self) -> &MolecularGasSpectrumHeating<'a> {
        &self.binding
    }
    /// Frozen peak source rates, in each original bank's source order.
    pub fn peak_rates(&self) -> &[Vec<C64>] {
        &self.peak_rates
    }
    pub fn source_clock(&self) -> &str {
        &self.source_clock
    }
    /// Native carrier/beat spectrum in watts; DC alone is not interval heat.
    pub const fn spectrum(&self) -> &PreparedHarmonicProductSpectrum {
        &self.spectrum
    }
    pub const fn workspace_bytes(&self) -> usize {
        self.workspace_bytes
    }

    /// Exact temporal integral of the retained, spatially quadratured heat.
    /// Start and duration remain separate, including on a large source clock.
    pub fn heat_load(
        &self,
        start_s: f64,
        duration_s: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularGasWindowHeatLoad<'_, 'a>, MolecularWindowError> {
        poll(&mut cancelled)?;
        super::admit_duration(duration_s)?;
        let window = HarmonicWindow::new(start_s, duration_s)?;
        let mean_power_w = self.spectrum.mean_cancellable(window, &mut cancelled)?;
        if mean_power_w < 0.0 {
            return Err(MolecularWindowError::NegativeHeat { mean_power_w });
        }
        let heat_j = super::absorbed_energy(mean_power_w, duration_s)?;
        poll(&mut cancelled)?;
        Ok(MolecularGasWindowHeatLoad {
            owner: self,
            start_s,
            duration_s,
            mean_power_w,
            heat_j,
        })
    }
}

/// One immutable physical interval, on the original accepted gas inventory.
#[derive(Clone, Copy)]
pub struct MolecularGasWindowHeatLoad<'w, 'a> {
    owner: &'w MolecularGasWindowHeating<'a>,
    start_s: f64,
    duration_s: f64,
    mean_power_w: f64,
    heat_j: f64,
}
impl<'w, 'a> MolecularGasWindowHeatLoad<'w, 'a> {
    pub const fn owner(&self) -> &'w MolecularGasWindowHeating<'a> {
        self.owner
    }
    pub const fn start_s(&self) -> f64 {
        self.start_s
    }
    pub const fn duration_s(&self) -> f64 {
        self.duration_s
    }
    pub const fn mean_absorbed_power_w(&self) -> f64 {
        self.mean_power_w
    }
    pub const fn absorbed_energy_j(&self) -> f64 {
        self.heat_j
    }

    /// The original nonlinear u(T), fixed mass, EOS and source-law refresh.
    /// The candidate does not update this immutable acoustic solution bank.
    pub fn heat_gas(
        &self,
        controls: GasEnergySolve,
        maximum_temperature_increment_k: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MolecularGasWindowHeatingStep<'w, 'a>, MolecularWindowError> {
        let (gas, next_relaxation, energy) = heated_gas(
            self.owner.binding.gas(),
            self.owner.binding.relaxation_law(),
            self.heat_j,
            controls,
            maximum_temperature_increment_k,
            &mut cancelled,
        )?;
        Ok(MolecularGasWindowHeatingStep {
            load: *self,
            gas,
            next_relaxation,
            energy,
        })
    }

    /// Feed these joules once into the existing reciprocal gas/body solver.
    /// It begins from the ORIGINAL gas, never an already-heated candidate.
    pub fn finite_body_environment<'b>(
        &'b self,
        body: &'b LumpedThermalEnvironment,
        gas_solve: GasEnergySolve,
        exchange: HeatExchangeSolve,
        maximum_gas_temperature_increment_k: f64,
    ) -> Result<FiniteGasThermalEnvironment<'b>, MolecularWindowError> {
        Ok(finite_body_environment(
            self.owner.binding.gas(),
            self.heat_j,
            body,
            gas_solve,
            exchange,
            maximum_gas_temperature_increment_k,
        )?)
    }
}

/// Atomic candidate with the native gas-storage balance and next material state.
pub struct MolecularGasWindowHeatingStep<'w, 'a> {
    load: MolecularGasWindowHeatLoad<'w, 'a>,
    gas: IsochoricGasParcel,
    next_relaxation: RelaxationState,
    energy: GasHeatStep,
}
impl<'w, 'a> MolecularGasWindowHeatingStep<'w, 'a> {
    pub const fn load(&self) -> &MolecularGasWindowHeatLoad<'w, 'a> {
        &self.load
    }
    pub const fn gas(&self) -> &IsochoricGasParcel {
        &self.gas
    }
    pub const fn next_relaxation(&self) -> &RelaxationState {
        &self.next_relaxation
    }
    pub const fn energy_balance(&self) -> &GasHeatStep {
        &self.energy
    }
    pub fn into_gas(self) -> IsochoricGasParcel {
        self.gas
    }
}
