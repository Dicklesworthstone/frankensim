//! Reversible isotropic volume work with kinetic alpha storage retained.
//!
//! Momentum/energy advection follows Parker (1965), printed p. 15, Eq. (2):
//! https://polaris.nipr.ac.jp/~ryuho/coem/Parker1965PSS.pdf

use super::*;
use crate::fusion::solve::Sum;
use crate::fusion::volume::{
    adiabatic_work_on_pools, advance_adiabatic_fusion_volume, log_volume_ratio,
};
use fs_math::det::exp;
use fs_qty::Volume;
use fs_thermochem::coulomb::CoulombScreeningV1;
use fs_thermochem::fusion::plasma::{DtFusionPlasmaErrorV1, DtFusionPlasmaLimitsV1};

/// Accepted fixed-inventory motion with separate work on all three stores.
/// Positive work is supplied by the boundary during compression. Neither
/// collision heat nor nuclear energy is included in these receipts.
#[derive(Clone, Debug, PartialEq)]
pub struct KineticFusionVolumeResult {
    initial: KineticFusionState,
    state: KineticFusionState,
    work_on_pools: [Energy; 3],
    energy_residual_j: f64,
}

impl KineticFusionVolumeResult {
    pub const fn initial_state(&self) -> &KineticFusionState {
        &self.initial
    }
    pub const fn state(&self) -> &KineticFusionState {
        &self.state
    }
    /// Mechanical work on the Maxwellian heavy background [J].
    pub const fn boundary_work_on_heavy(&self) -> Energy {
        self.work_on_pools[0]
    }
    /// Mechanical work on the Maxwellian electrons [J].
    pub const fn boundary_work_on_electrons(&self) -> Energy {
        self.work_on_pools[1]
    }
    /// Mechanical work on the retained isotropic alpha distribution [J].
    pub const fn boundary_work_on_alphas(&self) -> Energy {
        self.work_on_pools[2]
    }
    /// Total work on the parcel, opposite to the boundary's mechanical receipt.
    pub fn total_boundary_work(&self) -> Energy {
        let mut total = Sum::default();
        for work in self.work_on_pools {
            total.add(work.value());
        }
        Energy::new(total.value)
    }
    /// Independently evaluated source/grid energy change minus boundary work.
    /// This is a floating-point closure diagnostic, not a physical error bound.
    pub const fn energy_residual_j(&self) -> f64 {
        self.energy_residual_j
    }
}

/// Compress or expand an explicitly confined, co-moving isotropic hybrid parcel.
///
/// Nonrelativistic particle momenta scale as `V^(-1/3)`, so every alpha energy
/// center and both thermal temperatures scale as `V^(-2/3)`. Bin amounts are
/// unchanged: geometry neither diffuses the distribution nor hands particles
/// to the thermal pool. All three pressures obey `p_s = 2 U_s/(3 V)` and their
/// separate work receipts use the same stable adiabatic map.
///
/// The original source reference admits the complete monotone temperature path,
/// including absent products and source gaps. Alpha transport supplies the
/// original projectile identity and explicit dilute/nonrelativistic grid bounds;
/// no collision, radiation, reaction or physical clock advances here. Endpoint
/// admission covers this fixed-inventory path: number and energy dilution ratios
/// and upper-bin occupancy are invariant, and maximum grid energy is monotone.
/// The caller owns boundary motion and scheduling with other physics. This
/// closure supplies no confinement law, piston dynamics, shocks or anisotropic
/// magnetized compression. Cancellation or refusal publishes no partial result.
pub fn advance_adiabatic_kinetic_fusion_volume(
    initial: &KineticFusionState,
    alpha_transport: &KineticAlphaTransport,
    final_volume: Volume,
    plasma_limits: DtFusionPlasmaLimitsV1,
    mut cancelled: impl FnMut() -> bool,
) -> Result<KineticFusionVolumeResult, FusionNetworkError> {
    poll(&mut cancelled)?;
    let logarithm = log_volume_ratio(initial.volume(), final_volume)?;
    admit_alpha_state(initial, alpha_transport, &mut cancelled)?;

    // This internal reference has total canonical He4 at the background
    // temperature. Reuse its source/path admission, never its all-thermal work
    // or energy as the physical hybrid observable.
    let reference_motion = advance_adiabatic_fusion_volume(
        &initial.reference,
        initial.plasma().into(),
        final_volume,
        plasma_limits,
        &mut cancelled,
    )?;
    if logarithm == 0.0 {
        poll(&mut cancelled)?;
        return Ok(KineticFusionVolumeResult {
            initial: initial.clone(),
            state: initial.clone(),
            work_on_pools: [Energy::new(0.0); 3],
            energy_residual_j: 0.0,
        });
    }

    let before = [
        initial.thermal_heavy_energy().value(),
        initial.thermal_electron_energy().value(),
        initial.alpha_energy().value(),
    ];
    let work = adiabatic_work_on_pools(before, initial.volume(), final_volume)?;
    let factor = exp((-2.0 / 3.0) * logarithm);
    let grid = initial.energy_grid().scaled(factor, &mut cancelled)?;
    let state = KineticFusionState::from_parts(
        reference_motion.state(),
        initial.plasma(),
        &grid,
        initial.alpha_amounts,
        initial.thermal_helium4_amount,
        initial.alpha_amount().value(),
        initial.elapsed_time(),
        plasma_limits,
        &mut cancelled,
    )?;
    admit_alpha_state(&state, alpha_transport, &mut cancelled)?;

    let mut changes = [Sum::default(); 3];
    let old_temperatures = [
        initial.background_heavy_temperature(),
        initial.electron_temperature(),
    ];
    let new_temperatures = [
        state.background_heavy_temperature(),
        state.electron_temperature(),
    ];
    for (i, (source, amount)) in initial
        .plasma()
        .thermodynamics()
        .models()
        .iter()
        .zip(initial.thermal_amounts())
        .enumerate()
    {
        poll(&mut cancelled)?;
        if amount.value() == 0.0 {
            continue;
        }
        let pool = usize::from(i == initial.plasma().electron_index());
        let change = source
            .molar_internal_energy_change(old_temperatures[pool], new_temperatures[pool])
            .map_err(|error| {
                FusionNetworkError::Plasma(DtFusionPlasmaErrorV1::Source {
                    species: source.species().clone(),
                    source: error,
                })
            })?;
        changes[pool].add(amount.value() * change.value());
    }
    for ((amount, before), after) in initial
        .alpha_populations_mol()
        .iter()
        .zip(initial.energy_grid().energies())
        .zip(state.energy_grid().energies())
    {
        poll(&mut cancelled)?;
        changes[2].add((*amount * AVOGADRO) * (after.value() - before.value()));
    }

    let after = [
        state.thermal_heavy_energy().value(),
        state.thermal_electron_energy().value(),
        state.alpha_energy().value(),
    ];
    let mut residual = Sum::default();
    let mut total_work = Sum::default();
    for i in 0..3 {
        let defect = changes[i].value - work[i];
        let storage_defect = (after[i] - before[i]) - changes[i].value;
        let tolerance = 512.0 * f64::EPSILON * before[i].abs().max(after[i].abs());
        if !defect.is_finite()
            || !storage_defect.is_finite()
            || defect.abs() > tolerance
            || storage_defect.abs() > tolerance
        {
            return Err(FusionNetworkError::Balance(
                "adiabatic hybrid source/grid energy and separate pressure work",
            ));
        }
        residual.add(defect);
        total_work.add(work[i]);
    }
    if !residual.value.is_finite() || !total_work.value.is_finite() {
        return Err(FusionNetworkError::Arithmetic(
            "finite adiabatic hybrid work and energy closure",
        ));
    }
    poll(&mut cancelled)?;
    Ok(KineticFusionVolumeResult {
        initial: initial.clone(),
        state,
        work_on_pools: work.map(Energy::new),
        energy_residual_j: residual.value,
    })
}

/// Geometry preserves these ratios, but its physical energy support changes.
/// Collision-time/source preparation remains the evolution owner's operation.
pub(super) fn admit_alpha_state(
    state: &KineticFusionState,
    alpha: &KineticAlphaTransport,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), FusionNetworkError> {
    poll(cancelled)?;
    let source = alpha.collisions();
    if source.background().initial_state().thermodynamics() != state.plasma().thermodynamics()
        || source.projectile_index() != state.plasma().indices()[4]
        || source.background().screening() != CoulombScreeningV1::ElectronDebye
    {
        return Err(FusionNetworkError::Invalid(
            "original He4 collision source with ElectronDebye screening required",
        ));
    }
    let grid = state.energy_grid();
    let (number, energy) = super::state::moments(grid, &state.alpha_amounts);
    let thermal = state.thermal_amounts();
    let mut thermal_ions = Sum::default();
    for i in state.plasma().indices() {
        thermal_ions.add(thermal[i].value());
    }
    let thermal_energy =
        state.thermal_heavy_energy().value() + state.thermal_electron_energy().value();
    if !number.is_finite()
        || !energy.is_finite()
        || grid.energies()[grid.len() - 1] > source.maximum_energy()
        || number > alpha.maximum_number_ratio() * thermal_ions.value
        || energy > alpha.maximum_energy_ratio() * thermal_energy
        || (number > 0.0
            && state.alpha_populations_mol()[grid.len() - 1]
                > alpha.maximum_upper_bin_fraction() * number)
    {
        return Err(FusionNetworkError::OutsideDomain);
    }
    poll(cancelled)?;
    Ok(())
}

#[cfg(test)]
mod tests;
