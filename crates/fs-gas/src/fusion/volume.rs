//! Reversible volume work on the original two-temperature fusion plasma.
//!
//! At fixed populations each admitted monatomic pool obeys `p = 2 U/(3 V)`.
//! A quasistatic isotropic volume change therefore gives `T V^(2/3) = constant`
//! separately for ions and electrons. Mechanical work is retained separately
//! from heat, fusion, Coulomb exchange and radiation; this endpoint operation
//! advances no clock. The caller owns the boundary motion and its scheduling
//! with other physics. Piston inertia, shocks and anisotropic magnetic
//! compression are outside this homogeneous thermal model.

use super::solve::Sum;
use super::{DtFusionError, poll, positive};
use fs_math::det::{exp, expm1, ln};
use fs_qty::{Energy, Temperature, Volume};
use fs_thermochem::equilibrium::composition::{
    FrozenIonizedGasStateErrorV1, FrozenIonizedGasStateV1,
};
use fs_thermochem::fusion::plasma::{
    DtFusionPlasmaErrorV1, DtFusionPlasmaLimitsV1, DtFusionPlasmaV1, LightIonFusionPlasmaV1,
};

mod history;
pub use history::{FusionVolumeHistory, FusionVolumeKnot};

/// Original isotope-to-material binding for a supported thermal fusion bank.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum FusionPlasmaRef<'a> {
    /// Four original D+, T+, He4++ and electron sources.
    Dt(&'a DtFusionPlasmaV1),
    /// Six original p+, D+, T+, He3++, He4++ and electron sources.
    LightIon(&'a LightIonFusionPlasmaV1),
}

impl<'a> From<&'a DtFusionPlasmaV1> for FusionPlasmaRef<'a> {
    fn from(plasma: &'a DtFusionPlasmaV1) -> Self {
        Self::Dt(plasma)
    }
}
impl<'a> From<&'a LightIonFusionPlasmaV1> for FusionPlasmaRef<'a> {
    fn from(plasma: &'a LightIonFusionPlasmaV1) -> Self {
        Self::LightIon(plasma)
    }
}

impl FusionPlasmaRef<'_> {
    fn admit(
        self,
        state: &FrozenIonizedGasStateV1,
        limits: DtFusionPlasmaLimitsV1,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), DtFusionError> {
        match self {
            Self::Dt(plasma) => plasma.admit_state_with_limits(state, limits, cancelled),
            Self::LightIon(plasma) => plasma.admit_state_with_limits(state, limits, cancelled),
        }
        .map(|_| ())
        .map_err(plasma_error)
    }

    fn admit_interval(
        self,
        initial: [Temperature; 2],
        final_temperatures: [Temperature; 2],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), DtFusionError> {
        let lower = core::array::from_fn::<_, 2, _>(|i| {
            Temperature::new(initial[i].value().min(final_temperatures[i].value()))
        });
        let upper = core::array::from_fn::<_, 2, _>(|i| {
            Temperature::new(initial[i].value().max(final_temperatures[i].value()))
        });
        match self {
            Self::Dt(plasma) => {
                plasma.admit_temperature_interval(lower[0], upper[0], lower[1], upper[1], cancelled)
            }
            Self::LightIon(plasma) => {
                plasma.admit_temperature_interval(lower[0], upper[0], lower[1], upper[1], cancelled)
            }
        }
        .map_err(plasma_error)
    }
}

/// Accepted fixed-inventory motion and its reciprocal mechanical-work receipt.
/// Both states retain the exact original material bank and canonical amounts.
#[derive(Clone, Debug, PartialEq)]
pub struct FusionVolumeResult {
    initial: FrozenIonizedGasStateV1,
    state: FrozenIonizedGasStateV1,
    work_on_pools: [Energy; 2],
    energy_residual_j: f64,
}

impl FusionVolumeResult {
    /// Accepted original-bank state at the prescribed new volume.
    pub const fn state(&self) -> &FrozenIonizedGasStateV1 {
        &self.state
    }
    /// Unmodified input state.
    pub const fn initial_state(&self) -> &FrozenIonizedGasStateV1 {
        &self.initial
    }
    /// Mechanical work on heavy particles [J]; positive for compression.
    pub const fn boundary_work_on_heavy(&self) -> Energy {
        self.work_on_pools[0]
    }
    /// Mechanical work on electrons [J]; positive for compression.
    pub const fn boundary_work_on_electrons(&self) -> Energy {
        self.work_on_pools[1]
    }
    /// Total work on the gas [J], opposite to the boundary's mechanical receipt.
    pub fn total_boundary_work(&self) -> Energy {
        Energy::new(self.work_on_pools[0].value() + self.work_on_pools[1].value())
    }
    /// Original-source energy change minus mechanical work [J]. This is a
    /// floating-point closure diagnostic, not a physical error certificate.
    pub const fn energy_residual_j(&self) -> f64 {
        self.energy_residual_j
    }
}

/// Change a fusion parcel's volume reversibly, with both thermal pools retained.
///
/// Every source, including absent possible products, must support the complete
/// monotone temperature path without a gap. Both endpoints must obey the
/// authored classical-plasma limits. These admit the complete ideal adiabatic
/// path: coupling scales as `V^(1/3)`, thermal phase-space density is constant,
/// and thermal/rest-energy ratio scales as `V^(-2/3)` at fixed populations.
///
/// No populations or original source cards are changed. An exact volume no-op
/// returns the original state and exactly zero work after source/domain
/// admission. Refusal or cancellation publishes no result.
pub fn advance_adiabatic_fusion_volume(
    initial: &FrozenIonizedGasStateV1,
    plasma: FusionPlasmaRef<'_>,
    final_volume: Volume,
    limits: DtFusionPlasmaLimitsV1,
    mut cancelled: impl FnMut() -> bool,
) -> Result<FusionVolumeResult, DtFusionError> {
    poll(&mut cancelled)?;
    let logarithm = log_volume_ratio(initial.volume(), final_volume)?;
    plasma.admit(initial, limits, &mut cancelled)?;
    if logarithm == 0.0 {
        poll(&mut cancelled)?;
        return Ok(FusionVolumeResult {
            initial: initial.clone(),
            state: initial.clone(),
            work_on_pools: [Energy::new(0.0); 2],
            energy_residual_j: 0.0,
        });
    }
    let before = [
        initial.heavy_internal_energy().value(),
        initial.electron_internal_energy().value(),
    ];
    let work = adiabatic_work_on_pools(before, initial.volume(), final_volume)?;
    let old_temperatures = [initial.heavy_temperature(), initial.electron_temperature()];
    // Evaluate the positive endpoint directly, rather than subtracting nearly
    // all of U after a deep expansion. expm1 retains small mechanical work.
    let factor = exp((-2.0 / 3.0) * logarithm);
    let temperatures = old_temperatures.map(|t| Temperature::new(t.value() * factor));
    if temperatures.iter().any(|t| !positive(t.value())) {
        return Err(DtFusionError::Arithmetic("adiabatic fusion temperature"));
    }
    plasma.admit_interval(old_temperatures, temperatures, &mut cancelled)?;
    let state = FrozenIonizedGasStateV1::new(
        initial.thermodynamics(),
        initial.amounts(),
        final_volume,
        temperatures[0],
        temperatures[1],
        64.0 * f64::EPSILON,
        &mut cancelled,
    )
    .map_err(|error| match error {
        FrozenIonizedGasStateErrorV1::Cancelled => DtFusionError::Cancelled,
        error => DtFusionError::State(error),
    })?;
    plasma.admit(&state, limits, &mut cancelled)?;

    // Independently evaluate changes through the original caloric sources.
    let mut changes = [Sum::default(); 2];
    for (i, (amount, source)) in initial
        .amounts()
        .iter()
        .zip(initial.thermodynamics().models())
        .enumerate()
    {
        poll(&mut cancelled)?;
        if amount.value() == 0.0 {
            continue;
        }
        let pool = usize::from(Some(i) == initial.electron_index());
        let change = source
            .molar_internal_energy_change(old_temperatures[pool], temperatures[pool])
            .map_err(|error| {
                DtFusionError::Plasma(DtFusionPlasmaErrorV1::Source {
                    species: source.species().clone(),
                    source: error,
                })
            })?;
        changes[pool].add(amount.value() * change.value());
    }
    let after = [
        state.heavy_internal_energy().value(),
        state.electron_internal_energy().value(),
    ];
    let mut residual = Sum::default();
    for i in 0..2 {
        let defect = changes[i].value - work[i];
        let storage_defect = (after[i] - before[i]) - changes[i].value;
        let tolerance = 512.0 * f64::EPSILON * before[i].abs().max(after[i].abs());
        if !defect.is_finite()
            || !storage_defect.is_finite()
            || defect.abs() > tolerance
            || storage_defect.abs() > tolerance
        {
            return Err(DtFusionError::Balance("adiabatic fusion pool energy"));
        }
        residual.add(defect);
    }
    if !(work[0] + work[1]).is_finite() || !residual.value.is_finite() {
        return Err(DtFusionError::Arithmetic("adiabatic fusion total work"));
    }
    poll(&mut cancelled)?;
    Ok(FusionVolumeResult {
        initial: initial.clone(),
        state,
        work_on_pools: work.map(Energy::new),
        energy_residual_j: residual.value,
    })
}

/// Signed work on fixed-stock nonrelativistic isotropic components. A kinetic
/// component may be empty; each thermal owner still admits positive storage.
pub(super) fn adiabatic_work_on_pools<const N: usize>(
    energies: [f64; N],
    from: Volume,
    to: Volume,
) -> Result<[f64; N], DtFusionError> {
    if energies
        .iter()
        .any(|&energy| !energy.is_finite() || energy < 0.0)
    {
        return Err(DtFusionError::Invalid(
            "nonnegative finite adiabatic component energies",
        ));
    }
    let logarithm = log_volume_ratio(from, to)?;
    if logarithm == 0.0 {
        return Ok([0.0; N]);
    }
    let fraction = expm1((-2.0 / 3.0) * logarithm);
    let work = energies.map(|energy| energy * fraction);
    if !fraction.is_finite() || work.iter().any(|value| !value.is_finite()) {
        return Err(DtFusionError::Arithmetic("adiabatic fusion work"));
    }
    Ok(work)
}

/// Stable logarithmic volume change, preserving adjacent floating-point volumes.
pub(super) fn log_volume_ratio(from: Volume, to: Volume) -> Result<f64, DtFusionError> {
    let (from, to) = (from.value(), to.value());
    if !positive(from) || !positive(to) {
        return Err(DtFusionError::Invalid("positive finite fusion volumes"));
    }
    if from == to {
        return Ok(0.0);
    }
    let relative = (to - from) / from;
    let result = if relative.abs() <= 1e-4 {
        // Fifth-order ln(1+x): remainder <= |x|^6/[6(1-|x|)], below
        // 1.667e-25 here, before f64 rounding. This bounds only polynomial
        // truncation, not source/model error or all floating-point arithmetic.
        relative
            * (1.0
                + relative * (-0.5 + relative * (1.0 / 3.0 + relative * (-0.25 + relative * 0.2))))
    } else {
        let ratio = to / from;
        if ratio.is_normal() {
            ln(ratio)
        } else {
            // A ratio can overflow or lose relative accuracy in subnormal
            // arithmetic even though both authored volumes are finite.
            ln(to) - ln(from)
        }
    };
    if !result.is_finite() {
        return Err(DtFusionError::Arithmetic("logarithmic fusion volume ratio"));
    }
    Ok(result)
}

fn plasma_error(error: DtFusionPlasmaErrorV1) -> DtFusionError {
    match error {
        DtFusionPlasmaErrorV1::Cancelled => DtFusionError::Cancelled,
        error => DtFusionError::Plasma(error),
    }
}

#[cfg(test)]
mod tests;
