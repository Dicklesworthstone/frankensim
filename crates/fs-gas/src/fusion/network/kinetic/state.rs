use super::*;
use crate::fusion::solve::Sum;
use fs_qty::{Pressure, Temperature, Volume};
use fs_thermochem::fusion::plasma::{DtFusionPlasmaErrorV1, DtFusionPlasmaLimitsV1};

/// Original-material plasma with a Maxwellian background and a separate
/// isotropic He4 energy distribution. Total canonical amounts include the
/// kinetic alphas, while background energy and heat capacity exclude them.
#[derive(Clone, Debug, PartialEq)]
pub struct KineticFusionState {
    pub(super) reference: FrozenIonizedGasStateV1,
    plasma: LightIonFusionPlasmaV1,
    grid: AlphaEnergyGrid,
    pub(super) alpha_amounts: [f64; MAX_ALPHA_BINS],
    pub(super) thermal_helium4_amount: f64,
    thermal_energy: [Energy; 2],
    alpha_amount: Amount,
    alpha_energy: Energy,
    pressure: Pressure,
    elapsed: Time,
}
impl KineticFusionState {
    /// Interpret a declared subset of the reference's He4 inventory as kinetic
    /// populations at the supplied energy centers. The reference retains total
    /// quasineutral amounts, background temperatures and original sources; its
    /// all-thermal pressure/energy are not the resulting hybrid observables.
    /// Initial kinetic energy is authored initial-condition energy, not heat
    /// silently deposited into either Maxwellian pool.
    pub fn new(
        reference: &FrozenIonizedGasStateV1,
        plasma: &LightIonFusionPlasmaV1,
        grid: &AlphaEnergyGrid,
        alpha_populations: &[Amount],
        plasma_limits: DtFusionPlasmaLimitsV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, FusionNetworkError> {
        poll(&mut cancelled)?;
        if alpha_populations.len() != grid.len() {
            return Err(FusionNetworkError::Invalid(
                "one alpha population per energy center required",
            ));
        }
        let mut alpha_amounts = [0.0; MAX_ALPHA_BINS];
        for (out, n) in alpha_amounts.iter_mut().zip(alpha_populations) {
            if !n.value().is_finite() || n.value() < 0.0 {
                return Err(FusionNetworkError::Invalid(
                    "finite nonnegative alpha populations required",
                ));
            }
            *out = n.value();
        }
        let amount = moments(grid, &alpha_amounts).0;
        let total = reference
            .amounts()
            .get(plasma.indices()[4])
            .ok_or(FusionNetworkError::Invalid(
                "original canonical fusion amounts required",
            ))?
            .value();
        if amount > total {
            return Err(FusionNetworkError::Invalid(
                "kinetic He4 must be a subset of total He4 inventory",
            ));
        }
        Self::from_parts(
            reference,
            plasma,
            grid,
            alpha_amounts,
            total - amount,
            amount,
            Time::new(0.0),
            plasma_limits,
            &mut cancelled,
        )
    }

    /// Start a kinetic burn from an ordinary admitted thermal plasma, with no
    /// artificial seed in any alpha energy bin.
    pub fn from_thermal(
        initial: &FrozenIonizedGasStateV1,
        plasma: &LightIonFusionPlasmaV1,
        grid: &AlphaEnergyGrid,
        plasma_limits: DtFusionPlasmaLimitsV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, FusionNetworkError> {
        let amount = initial
            .amounts()
            .get(plasma.indices()[4])
            .ok_or(FusionNetworkError::Invalid(
                "original canonical fusion amounts required",
            ))?
            .value();
        Self::from_parts(
            initial,
            plasma,
            grid,
            [0.0; MAX_ALPHA_BINS],
            amount,
            0.0,
            Time::new(0.0),
            plasma_limits,
            &mut cancelled,
        )
    }

    pub(super) fn from_parts(
        reference: &FrozenIonizedGasStateV1,
        plasma: &LightIonFusionPlasmaV1,
        grid: &AlphaEnergyGrid,
        alpha_amounts: [f64; MAX_ALPHA_BINS],
        thermal_helium4_amount: f64,
        alpha_inventory_amount: f64,
        elapsed: Time,
        plasma_limits: DtFusionPlasmaLimitsV1,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self, FusionNetworkError> {
        plasma
            .admit_state_with_limits(reference, plasma_limits, &mut *cancelled)
            .map_err(|e| match e {
                DtFusionPlasmaErrorV1::Cancelled => FusionNetworkError::Cancelled,
                e => FusionNetworkError::Plasma(e),
            })?;
        if !elapsed.value().is_finite()
            || elapsed.value() < 0.0
            || !thermal_helium4_amount.is_finite()
            || thermal_helium4_amount < 0.0
            || !alpha_inventory_amount.is_finite()
            || alpha_inventory_amount < 0.0
            || alpha_amounts.iter().any(|n| !n.is_finite() || *n < 0.0)
            || alpha_amounts[grid.len()..].iter().any(|n| *n != 0.0)
        {
            return Err(FusionNetworkError::Invalid(
                "finite hybrid time and canonical populations",
            ));
        }
        let (alpha_amount, alpha_energy) = moments(grid, &alpha_amounts);
        let mut energy = [Sum::default(); 2];
        for (i, (source, amount)) in plasma
            .thermodynamics()
            .models()
            .iter()
            .zip(reference.amounts())
            .enumerate()
        {
            poll(cancelled)?;
            let n = if i == plasma.indices()[4] {
                thermal_helium4_amount
            } else {
                amount.value()
            };
            if n == 0.0 {
                continue;
            }
            let electron = i == plasma.electron_index();
            let temperature = if electron {
                reference.electron_temperature()
            } else {
                reference.heavy_temperature()
            };
            let properties = source.evaluate_properties(temperature).map_err(|e| {
                FusionNetworkError::Plasma(DtFusionPlasmaErrorV1::Source {
                    species: source.species().clone(),
                    source: e,
                })
            })?;
            energy[usize::from(electron)].add(n * properties.u().value());
        }
        let thermal_energy = energy.map(|u| Energy::new(u.value));
        let pressure = (2.0 / 3.0)
            * ((energy[0].value + energy[1].value + alpha_energy) / reference.volume().value());
        if thermal_energy.iter().any(|u| !positive(u.value()))
            || !positive(pressure)
            || !alpha_amount.is_finite()
            || !alpha_energy.is_finite()
        {
            return Err(FusionNetworkError::Arithmetic(
                "hybrid source energy or pressure",
            ));
        }
        poll(cancelled)?;
        Ok(Self {
            reference: reference.clone(),
            plasma: plasma.clone(),
            grid: grid.clone(),
            alpha_amounts,
            thermal_helium4_amount,
            thermal_energy,
            alpha_amount: Amount::new(alpha_inventory_amount),
            alpha_energy: Energy::new(alpha_energy),
            pressure: Pressure::new(pressure),
            elapsed,
        })
    }

    pub const fn plasma(&self) -> &LightIonFusionPlasmaV1 {
        &self.plasma
    }
    pub const fn energy_grid(&self) -> &AlphaEnergyGrid {
        &self.grid
    }
    /// Canonical total amounts, including both thermal and kinetic He4.
    pub fn amounts(&self) -> &[Amount] {
        self.reference.amounts()
    }
    /// Amount [mol] at each original energy center; no bins are seeded or removed.
    pub fn alpha_populations_mol(&self) -> &[f64] {
        &self.alpha_amounts[..self.grid.len()]
    }
    /// Independently retained kinetic alpha inventory; distribution roundoff
    /// against this count is checked and reported by the evolution owner.
    pub const fn alpha_amount(&self) -> Amount {
        self.alpha_amount
    }
    pub const fn alpha_energy(&self) -> Energy {
        self.alpha_energy
    }
    pub const fn thermal_heavy_energy(&self) -> Energy {
        self.thermal_energy[0]
    }
    pub const fn thermal_electron_energy(&self) -> Energy {
        self.thermal_energy[1]
    }
    pub fn internal_energy(&self) -> Energy {
        Energy::new(
            self.thermal_energy[0].value()
                + self.thermal_energy[1].value()
                + self.alpha_energy.value(),
        )
    }
    /// Temperature of the Maxwellian heavy background, excluding kinetic He4.
    pub fn background_heavy_temperature(&self) -> Temperature {
        self.reference.heavy_temperature()
    }
    pub fn electron_temperature(&self) -> Temperature {
        self.reference.electron_temperature()
    }
    pub fn volume(&self) -> Volume {
        self.reference.volume()
    }
    /// Sum of background and isotropic nonrelativistic alpha kinetic pressures.
    pub const fn pressure(&self) -> Pressure {
        self.pressure
    }
    pub const fn elapsed_time(&self) -> Time {
        self.elapsed
    }
    pub fn thermal_amounts(&self) -> [Amount; 6] {
        let mut amounts = core::array::from_fn(|i| self.reference.amounts()[i]);
        amounts[self.plasma.indices()[4]] = Amount::new(self.thermal_helium4_amount);
        amounts
    }
}

pub(super) fn moments(grid: &AlphaEnergyGrid, amounts: &[f64; MAX_ALPHA_BINS]) -> (f64, f64) {
    let mut number = Sum::default();
    let mut energy = Sum::default();
    for (n, e) in amounts.iter().zip(grid.energies()) {
        number.add(*n);
        energy.add((*n * AVOGADRO) * e.value());
    }
    (number.value, energy.value)
}
