//! Finite-reservoir stream exchange followed by native Gibbs membrane motion.
//!
//! The original nozzle transfers species and upstream total enthalpy at fixed
//! volumes. The original membrane step then owns heat and moving-boundary work.
//! This is a first-order split over one physical interval, not a resolved jet.

use super::*;
use fs_material::thermochemical_gas::energy::parcel::transfer::GasMassTransferStep;
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::{
    GasNozzle, GasNozzleFlow, GasNozzleStepSolve,
};
use fs_material::thermochemical_gas::nonideal::energy::parcel::transfer::nozzle::{
    RealGasNozzleFlow, RealGasNozzleStepSolve,
};

/// Original source controls; both connected parcels must use this source family.
/// The native donor/receiver brackets, donor-fraction bound and pressure-crossing
/// refusal apply without clipping a stream or changing the physical interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum GasMembraneNozzleSolve {
    Ideal(GasNozzleStepSolve),
    PengRobinson(RealGasNozzleStepSolve),
}

impl GasMembraneNozzleSolve {
    fn mass_tolerance_kg(self) -> f64 {
        match self {
            Self::Ideal(controls) => controls.transfer.mass_tolerance_kg,
            Self::PengRobinson(controls) => controls.transfer.mass_tolerance_kg,
        }
    }
}

/// Actual upstream side selected from the two accepted absolute pressures.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum GasMembraneNozzleDirection {
    IntoChamber,
    OutOfChamber,
}

/// Original source-resolved throat, choking and stream enthalpy receipt.
/// The PR branch keeps its nonideal acoustic and caloric state.
#[derive(Clone, Debug, PartialEq)]
pub enum GasMembraneNozzleFlow {
    Ideal(GasNozzleFlow),
    PengRobinson(RealGasNozzleFlow),
}

/// One accepted finite-reservoir transfer and Gibbs membrane interval.
///
/// Transported enthalpy is internal to the reported pair. The original transfer
/// residual cancels conserved species energy datums before rounding; summing the
/// separately rounded extensive storage fields need not recover that residual.
#[derive(Clone, Debug)]
pub struct GasMembraneNozzleStep {
    /// None when the original source transfers zero mass.
    pub direction: Option<GasMembraneNozzleDirection>,
    pub flow: GasMembraneNozzleFlow,
    /// Original per-species inventory, upstream enthalpy and transfer balances.
    /// Ideal and PR transfer owners use this same source receipt.
    pub transfer: GasMassTransferStep,
    /// Native Gibbs mechanics starting from the post-transfer chamber.
    /// Its prescribed heat contains no copy of transported enthalpy.
    pub motion: GasMembraneDiscreteGradientStep,
    /// Whole-interval first-chamber storage: stream plus heat/volume stages [J].
    pub chamber_energy_change_j: f64,
    /// Fixed-volume finite-reservoir storage over the interval [J].
    pub reservoir_energy_change_j: f64,
    /// Independently measured net mass change of the two accepted endpoints [kg].
    pub total_mass_residual_kg: f64,
    /// Delta U of both connected gases and any sealed opposing gas, plus Delta K
    /// and membrane storage, plus exterior work, minus dead work and prescribed
    /// heats. Assembled from the original reference-safe transfer residual and
    /// the independently gated native mechanical integration defect [J].
    pub total_energy_residual_j: f64,
}

struct Transfer {
    chamber: MembraneGas,
    reservoir: MembraneGas,
    inflow: bool,
    flow: GasMembraneNozzleFlow,
    report: GasMassTransferStep,
}

fn transfer(
    chamber: &MembraneGas,
    reservoir: &MembraneGas,
    nozzle: GasNozzle,
    duration_s: f64,
    controls: GasMembraneNozzleSolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Transfer, GasMembraneMotionError> {
    poll(cancelled)?;
    match (chamber, reservoir, controls) {
        (
            MembraneGas::Ideal(chamber),
            MembraneGas::Ideal(reservoir),
            GasMembraneNozzleSolve::Ideal(controls),
        ) => {
            let inflow = reservoir.state().state().pressure >= chamber.state().state().pressure;
            let (candidate, flow) = if inflow {
                reservoir.through_nozzle_to(
                    chamber,
                    nozzle,
                    duration_s,
                    controls,
                    &mut *cancelled,
                )?
            } else {
                chamber.through_nozzle_to(
                    reservoir,
                    nozzle,
                    duration_s,
                    controls,
                    &mut *cancelled,
                )?
            };
            let (donor, receiver, report) = candidate.into_parts();
            let (chamber, reservoir) = if inflow {
                (receiver, donor)
            } else {
                (donor, receiver)
            };
            Ok(Transfer {
                chamber: MembraneGas::Ideal(chamber),
                reservoir: MembraneGas::Ideal(reservoir),
                inflow,
                flow: GasMembraneNozzleFlow::Ideal(flow),
                report,
            })
        }
        (
            MembraneGas::PengRobinson(chamber),
            MembraneGas::PengRobinson(reservoir),
            GasMembraneNozzleSolve::PengRobinson(controls),
        ) => {
            // Move the unchanged source, composition, branch and auxiliary axes
            // through the original transfer owner; no gas is re-fitted.
            let chamber = chamber.clone().into_isochoric();
            let reservoir = reservoir.clone().into_isochoric();
            let inflow = reservoir.state().pressure_pa() >= chamber.state().pressure_pa();
            let (candidate, flow) = if inflow {
                reservoir.through_nozzle_to(
                    &chamber,
                    nozzle,
                    duration_s,
                    controls,
                    &mut *cancelled,
                )?
            } else {
                chamber.through_nozzle_to(
                    &reservoir,
                    nozzle,
                    duration_s,
                    controls,
                    &mut *cancelled,
                )?
            };
            let (donor, receiver, report) = candidate.into_parts();
            let (chamber, reservoir) = if inflow {
                (receiver, donor)
            } else {
                (donor, receiver)
            };
            Ok(Transfer {
                chamber: MembraneGas::PengRobinson(ClosedRealGasParcel::from_isochoric(chamber)),
                reservoir: MembraneGas::PengRobinson(ClosedRealGasParcel::from_isochoric(
                    reservoir,
                )),
                inflow,
                flow: GasMembraneNozzleFlow::PengRobinson(flow),
                report,
            })
        }
        _ => Err(GasMembraneMotionError::Invalid(
            "membrane nozzle needs two parcels and controls from the same original gas family",
        )),
    }
}

impl GasMembraneThermalMotion<'_, '_> {
    /// Exchange with a finite reservoir, then advance native pressure-work motion.
    ///
    /// The original source nozzle selects inflow/outflow from accepted pressure,
    /// retains choking and actual species, and transfers upstream total enthalpy
    /// once. The resulting chamber drives the original discrete Gibbs work
    /// pressure and exact swept-volume traction. A sealed opposing chamber keeps
    /// its independent source and heat control. The reservoir volume stays fixed.
    ///
    /// Both stages approximate their part of the SAME duration. This first-order
    /// transfer-then-motion split requires physical-time refinement; its energy
    /// gate is not a splitting-error estimate. No second heat copy of incoming
    /// enthalpy, idealized PR correction or jet thrust is added. The opening
    /// belongs to the rigid enclosure; jet/pipe momentum and its reaction-force
    /// receipt are unresolved. Membrane impingement needs a resolved momentum load.
    ///
    /// Reservoir, both chamber states, mechanics, clock and original history
    /// publish together only after source, native work/energy, whole-pair balance
    /// and final cancellation checks. A closed nozzle is the unchanged native
    /// Gibbs step. Subsequent closed, open or thermal-contact calls continue the
    /// accepted mixture and trajectory, including all original source receipts.
    #[allow(clippy::too_many_arguments)]
    pub fn step_discrete_gradient_with_nozzle(
        &mut self,
        reservoir: &mut MembraneGas,
        nozzle: GasNozzle,
        duration_s: f64,
        dead_load_n: &[[f64; 3]],
        heat_j: f64,
        exterior_heat_j: Option<f64>,
        exchange: GasMembraneNozzleSolve,
        controls: GasMembraneThermalControls,
        total_energy_tolerance_j: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<GasMembraneNozzleStep, GasMembraneMotionError> {
        poll(&mut cancelled)?;
        if !total_energy_tolerance_j.is_finite() || total_energy_tolerance_j < 0.0 {
            return Err(GasMembraneMotionError::Invalid(
                "finite nonnegative nozzle/membrane whole-pair energy allowance required",
            ));
        }
        let stream = transfer(
            &self.accepted.gas,
            reservoir,
            nozzle,
            duration_s,
            exchange,
            &mut cancelled,
        )?;
        // Keep accepted history in place until joint publication. The private
        // trajectory starts with the existing algorithmic acceleration, which
        // the rho_inf=1 native midpoint solve retains without a clock restart.
        let mut next = Self {
            model: self.model,
            state: SecondOrderState::new(self.state.t, &self.state.q, &self.state.v, &self.state.a),
            accepted: GasMembraneThermalEvaluation {
                pressure: stream
                    .chamber
                    .pressure_response()
                    .map_err(GasMembraneError::from)?,
                gas: stream.chamber,
                displacement_m: self.accepted.displacement_m.clone(),
                exterior_gas: self.accepted.exterior_gas.clone(),
                exterior_pressure_pa: self.accepted.exterior_pressure_pa,
                residual_n: self.accepted.residual_n.clone(),
                strain_energy_j: self.accepted.strain_energy_j,
            },
            masses_kg: self.masses_kg.clone(),
        };
        next.state.steps = self.state.steps;
        // The original operator uses this parcel and the unchanged geometry;
        // its final source trial rebuilds complete physical pressure/reactions.
        let motion = next.step_discrete_gradient(
            duration_s,
            dead_load_n,
            heat_j,
            exterior_heat_j,
            controls,
            &mut cancelled,
        )?;
        let (stream_chamber, stream_reservoir) = if stream.inflow {
            (
                stream.report.receiver_energy_change_j,
                stream.report.donor_energy_change_j,
            )
        } else {
            (
                stream.report.donor_energy_change_j,
                stream.report.receiver_energy_change_j,
            )
        };
        let chamber_energy_change_j = (Dd::from_f64(stream_chamber)
            + Dd::from_f64(motion.gas.stored_energy_change_j))
        .to_f64();
        let reservoir_energy_change_j = stream_reservoir;
        let total_mass_residual_kg = (Dd::from_f64(next.accepted.gas.mass_kg())
            - Dd::from_f64(self.accepted.gas.mass_kg())
            + Dd::from_f64(stream.reservoir.mass_kg())
            - Dd::from_f64(reservoir.mass_kg()))
        .to_f64();
        let total_energy_residual_j = (Dd::from_f64(stream.report.energy_residual_j)
            + Dd::from_f64(motion.integration_defect_j))
        .to_f64();
        if ![
            chamber_energy_change_j,
            reservoir_energy_change_j,
            total_mass_residual_kg,
            total_energy_residual_j,
        ]
        .iter()
        .all(|value| value.is_finite())
            || total_mass_residual_kg.abs() > exchange.mass_tolerance_kg()
        {
            return Err(GasMembraneMotionError::Invalid(
                "nozzle/membrane endpoint mass or energy account is unrepresentable or outside its mass allowance",
            ));
        }
        if total_energy_residual_j.abs() > total_energy_tolerance_j {
            return Err(GasMembraneMotionError::Energy {
                defect_j: total_energy_residual_j,
                allowance_j: total_energy_tolerance_j,
            });
        }
        let direction = (stream.report.transferred_mass_kg != 0.0).then_some(if stream.inflow {
            GasMembraneNozzleDirection::IntoChamber
        } else {
            GasMembraneNozzleDirection::OutOfChamber
        });
        let report = GasMembraneNozzleStep {
            direction,
            flow: stream.flow,
            transfer: stream.report,
            motion,
            chamber_energy_change_j,
            reservoir_energy_change_j,
            total_mass_residual_kg,
            total_energy_residual_j,
        };
        poll(&mut cancelled)?;
        let entry = next
            .state
            .history
            .pop()
            .expect("one native accepted interval");
        next.state.history = std::mem::take(&mut self.state.history);
        next.state.history.push(entry);
        *reservoir = stream.reservoir;
        *self = next;
        Ok(report)
    }
}

#[cfg(test)]
mod tests;
