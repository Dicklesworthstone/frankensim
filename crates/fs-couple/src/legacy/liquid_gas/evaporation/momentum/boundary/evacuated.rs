//! Exact zero-inventory startup with retained molecular bulk momentum.
//!
//! The original evacuated-volume owner admits the captured species mass. Its
//! one energy callback removes bulk kinetic energy in the declared stream
//! frame before the ONE original caloric inversion. No absent-gas velocity,
//! temperature, seed mass or additional mixing heat is assigned to the vacuum.

use super::{MovingGasBoundaryError, MovingGasBoundarySolve, MovingGasParcel};
use fs_material::thermochemical_gas::{
    ThermochemicalGasModel,
    energy::{
        GasEnergyError,
        parcel::evacuated::{EvacuatedGasFillStep, EvacuatedGasVolume},
    },
};
use fs_math::dd::Dd;
use fs_thermochem::SpeciesId;
use std::collections::BTreeMap;

#[derive(Clone, Debug, PartialEq)]
pub struct MovingGasFillStep {
    pub caloric: EvacuatedGasFillStep,
    /// Momentum uses the original inertial axes; the initially absent gas has
    /// exactly zero momentum in every inertial frame.
    pub received_momentum_n_s: [f64; 3],
    pub energy_frame_velocity_m_s: [f64; 3],
    pub kinetic_energy_j: f64,
    pub laboratory_kinetic_energy_j: f64,
    /// Actual source-relative U plus actual bulk K minus the received energy.
    pub total_energy_residual_j: f64,
    /// Native caloric defect plus its independent total-energy reconstruction
    /// correction. Rounding the displayed stored energy cannot refund either.
    pub absolute_energy_defect_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct MovingGasFillCandidate {
    gas: MovingGasParcel,
    step: MovingGasFillStep,
}
impl MovingGasFillCandidate {
    pub const fn gas(&self) -> &MovingGasParcel {
        &self.gas
    }
    pub const fn step(&self) -> &MovingGasFillStep {
        &self.step
    }
    pub fn into_parts(self) -> (MovingGasParcel, MovingGasFillStep) {
        (self.gas, self.step)
    }
}

// Empty-to-filled energy, with the frame removed from momentum BEFORE the
// positive admitted mass is divided. The positive-mass changing-K owner and
// its contract are untouched; there is no division by an initial zero mass.
fn kinetic(mass: f64, momentum: [f64; 3], frame: [f64; 3]) -> Result<f64, GasEnergyError> {
    let mass = Dd::from_f64(mass);
    let mut energy = Dd::ZERO;
    for d in 0..3 {
        let relative = Dd::from_f64(momentum[d]) - mass * Dd::from_f64(frame[d]);
        energy = energy + (relative / mass) * relative * Dd::from_f64(0.5);
    }
    let rounded = energy.to_f64();
    if !energy.hi.is_finite()
        || !energy.lo.is_finite()
        || !rounded.is_finite()
        || (energy != Dd::ZERO && rounded == 0.0)
    {
        Err(GasEnergyError::Invalid(
            "captured bulk kinetic energy is unrepresentable",
        ))
    } else {
        Ok(rounded)
    }
}

impl MovingGasParcel {
    /// Construct a moving gas from conserved capture into exactly empty geometry.
    ///
    /// The energy coordinate is COMPLETE received energy in the declared
    /// inertial frame minus sum_i m_i*u_i(T_reference). Its thermal and kinetic
    /// terms must be formed in that frame directly, before large formation or
    /// common-boost totals are rounded. Only positive capture constructs a gas;
    /// an enclosing zero-transfer owner must retain the evacuated volume.
    #[allow(clippy::too_many_arguments)]
    pub fn fill_evacuated_candidate(
        empty: &EvacuatedGasVolume,
        catalog: &ThermochemicalGasModel,
        inventory_kg: &BTreeMap<SpeciesId, f64>,
        reference_temperature_k: f64,
        reference_relative_total_energy_j: f64,
        received_momentum_n_s: [f64; 3],
        energy_frame_velocity_m_s: [f64; 3],
        controls: MovingGasBoundarySolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MovingGasFillCandidate, MovingGasBoundaryError> {
        if cancelled() {
            return Err(MovingGasBoundaryError::Cancelled);
        }
        if received_momentum_n_s
            .into_iter()
            .chain(energy_frame_velocity_m_s)
            .any(|x| !x.is_finite())
            || !controls.momentum_absolute_n_s.is_finite()
            || controls.momentum_absolute_n_s <= 0.0
            || !controls.momentum_relative_tolerance.is_finite()
            || !(0.0..1.0).contains(&controls.momentum_relative_tolerance)
        {
            return Err(MovingGasBoundaryError::Invalid(
                "finite capture momentum/frame and positive bounds required",
            ));
        }
        let mut evaluated = None;
        let candidate = empty.fill_candidate_with_energy(
            catalog,
            inventory_kg,
            reference_temperature_k,
            reference_relative_total_energy_j,
            |mass| {
                let energy = kinetic(mass, received_momentum_n_s, energy_frame_velocity_m_s)?;
                evaluated = Some((mass, energy));
                Ok(-energy)
            },
            controls.gas.gas,
            controls.gas.mass_tolerance_kg,
            controls.gas.energy_tolerance_j,
            &mut cancelled,
        )?;
        let (parcel, caloric) = candidate.into_parts();
        let (mass, kinetic_energy_j) = evaluated.ok_or(MovingGasBoundaryError::Invalid(
            "original evacuated owner omitted the admitted-mass energy callback",
        ))?;
        if mass != parcel.mass_kg() {
            return Err(MovingGasBoundaryError::Invalid(
                "original fill changed its admitted mass",
            ));
        }
        let laboratory_kinetic_energy_j = kinetic(mass, received_momentum_n_s, [0.0; 3])?;
        let residual = Dd::from_f64(caloric.stored_relative_energy_j)
            + Dd::from_f64(kinetic_energy_j)
            - Dd::from_f64(reference_relative_total_energy_j);
        let total_energy_residual_j = residual.to_f64();
        let native = Dd::from_f64(caloric.energy_residual_j);
        let defect = native.abs() + (residual - native).abs();
        let absolute_energy_defect_j = defect.to_f64();
        if !total_energy_residual_j.is_finite()
            || Dd::from_f64(controls.gas.energy_tolerance_j).lt(residual.abs())
            || !absolute_energy_defect_j.is_finite()
            || Dd::from_f64(controls.gas.energy_tolerance_j).lt(defect)
        {
            return Err(MovingGasBoundaryError::Invalid(
                "captured total energy exceeds its fixed budget",
            ));
        }
        let gas = MovingGasParcel::from_momentum(parcel, received_momentum_n_s)
            .map_err(MovingGasBoundaryError::Motion)?;
        let step = MovingGasFillStep {
            caloric,
            received_momentum_n_s,
            energy_frame_velocity_m_s,
            kinetic_energy_j,
            laboratory_kinetic_energy_j,
            total_energy_residual_j,
            absolute_energy_defect_j,
        };
        if cancelled() {
            return Err(MovingGasBoundaryError::Cancelled);
        }
        Ok(MovingGasFillCandidate { gas, step })
    }
}
