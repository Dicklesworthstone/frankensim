//! A moving metering edge connects the two ORIGINAL finite gas inventories.
//!
//! The retained wall coordinate uncovers a rectangular window. Native nozzle
//! flow moves species and donor total enthalpy at the original opening; the
//! existing two-gas midpoint owner then solves both pressure works and the ONE
//! wall momentum. Neither gas is an infinite reservoir or an acoustic-flow
//! proxy. This first-order transfer-then-motion split needs time refinement.

use super::super::{GasEnergyError, PistonError, invalid, poll, positive};
use super::{DoubleActingGasPiston, DoublePistonInput, DoublePistonSolve, DoublePistonStep};
use fs_material::thermochemical_gas::energy::parcel::transfer::GasMassTransferStep;
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::{
    GasNozzle, GasNozzleFlow, GasNozzleStepSolve,
};
use fs_qty::Temperature;

/// One straight metering edge uncovering a fixed rectangular window.
///
/// This is overlap geometry, not a displacement clamp or a contact stop. The
/// reference volume merely locates the edge on the original wall coordinate;
/// it neither changes a gas volume nor creates another position state. The
/// pressure-exposed piston areas remain those of the original mechanical owner.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RectangularMeteringValve {
    /// Primary gas volume at the authored edge-position datum [m3].
    pub reference_primary_volume_m3: f64,
    /// Signed edge height above the window's lower end at that datum [m].
    /// Nonpositive covers the whole window; height above its upper end opens it.
    pub edge_height_at_reference_m: f64,
    /// Fixed window width perpendicular to wall motion [m].
    pub width_m: f64,
    /// Fixed window extent along wall motion [m].
    pub window_height_m: f64,
    /// True adds primary-expansion displacement to edge height; false subtracts.
    pub opens_on_primary_expansion: bool,
    /// Original nozzle coefficient in [0,1]; zero closes transport explicitly.
    pub discharge_coefficient: f64,
}

/// Actual opening evaluated from one accepted native wall position.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MeteringValveOpening {
    pub displacement_from_reference_m: f64,
    pub edge_height_m: f64,
    pub open_height_m: f64,
    pub area_m2: f64,
}

impl RectangularMeteringValve {
    /// Evaluate the exact rectangular overlap without moving either gas.
    pub fn opening(
        &self,
        piston: &DoubleActingGasPiston,
    ) -> Result<MeteringValveOpening, PistonError> {
        if !positive(self.reference_primary_volume_m3)
            || !self.edge_height_at_reference_m.is_finite()
            || !positive(self.width_m)
            || !positive(self.window_height_m)
            || !positive(self.width_m * self.window_height_m)
            || !self.discharge_coefficient.is_finite()
            || !(0.0..=1.0).contains(&self.discharge_coefficient)
        {
            return Err(invalid(
                "metering valve needs finite edge geometry, positive reference volume/window dimensions and a coefficient in [0,1]",
            ));
        }
        let displacement = (piston.primary.gas().volume_m3() - self.reference_primary_volume_m3)
            / piston.primary.area_m2();
        let edge = self.edge_height_at_reference_m
            + if self.opens_on_primary_expansion {
                displacement
            } else {
                -displacement
            };
        // Intersection of the uncovered half-plane with the fixed window.
        let height = edge.max(0.0).min(self.window_height_m);
        let area = self.width_m * height;
        if ![displacement, edge, height, area]
            .iter()
            .all(|v| v.is_finite())
            || (height > 0.0 && area == 0.0)
        {
            return Err(invalid(
                "metering edge position or exposed area is unrepresentable",
            ));
        }
        Ok(MeteringValveOpening {
            displacement_from_reference_m: displacement,
            edge_height_m: edge,
            open_height_m: height,
            area_m2: area,
        })
    }
}

/// Donor selected from the two original absolute pressures for this interval.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ValveDirection {
    PrimaryToOpposing,
    OpposingToPrimary,
}

/// Native transport and motion, with one enclosing energy and inventory audit.
#[derive(Clone, Debug, PartialEq)]
pub struct DoublePistonValveStep {
    /// Original opening used by the held nozzle flow over this interval.
    pub opening: MeteringValveOpening,
    /// Opening at the final wall position, used by the next accepted interval.
    pub endpoint_opening: MeteringValveOpening,
    /// None for closed or initially equilibrated transport.
    pub direction: Option<ValveDirection>,
    /// Native source-bearing throat state, mass rate and choking admission.
    pub flow: GasNozzleFlow,
    /// Original species and total-enthalpy transfer in its actual donor order.
    pub transfer: GasMassTransferStep,
    /// Both gas pressure works and the single shared wall momentum.
    pub piston: DoublePistonStep,
    /// Transfer plus subsequent native caloric storage, [primary, opposing] [J].
    /// These displayed values retain each original species' formation datum.
    pub gas_energy_change_j: [f64; 2],
    /// Actual final total mass minus original total mass [kg].
    pub total_mass_residual_kg: f64,
    /// dU_primary+dU_opposing+dK-Q_primary-Q_opposing-F*dx [J].
    /// Conserved formation datums cancel in the original transfer owner before
    /// arithmetic; do not reconstruct this small residual from displayed totals.
    pub total_energy_residual_j: f64,
}

/// Complete immutable proposal; neither accepted gas nor wall has advanced.
#[derive(Clone, Debug, PartialEq)]
pub struct DoublePistonValveCandidate {
    piston: DoubleActingGasPiston,
    step: DoublePistonValveStep,
}
impl DoublePistonValveCandidate {
    pub const fn piston(&self) -> &DoubleActingGasPiston {
        &self.piston
    }
    pub const fn step(&self) -> &DoublePistonValveStep {
        &self.step
    }
    pub fn into_parts(self) -> (DoubleActingGasPiston, DoublePistonValveStep) {
        (self.piston, self.step)
    }
}

// The finite-thermal adapter uses this same transport and final account. Only
// the original two-gas mechanical/thermal step between them varies.
pub(super) struct ValveTransfer {
    opening: MeteringValveOpening,
    direction: Option<ValveDirection>,
    flow: GasNozzleFlow,
    transfer: GasMassTransferStep,
    temperature_k: [f64; 2],
    mass_kg: [f64; 2],
}

impl DoubleActingGasPiston {
    pub(super) fn prepare_valve(
        &self,
        geometry: RectangularMeteringValve,
        duration_s: f64,
        exchange: GasNozzleStepSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Self, ValveTransfer), PistonError> {
        poll(cancelled)?;
        let opening = geometry.opening(self)?;
        let nozzle = GasNozzle {
            area_m2: opening.area_m2,
            discharge_coefficient: geometry.discharge_coefficient,
        };
        let primary = self.primary.gas();
        let forward = primary.state().state().pressure >= self.opposing.state().state().pressure;
        let (candidate, flow) = if forward {
            primary.through_nozzle_to(
                &self.opposing,
                nozzle,
                duration_s,
                exchange,
                &mut *cancelled,
            )?
        } else {
            self.opposing.through_nozzle_to(
                primary,
                nozzle,
                duration_s,
                exchange,
                &mut *cancelled,
            )?
        };
        let (donor, receiver, transfer) = candidate.into_parts();
        let (primary, opposing) = if forward {
            (donor, receiver)
        } else {
            (receiver, donor)
        };
        let temperature_k = [
            primary.state().state().temperature,
            opposing.state().state().temperature,
        ];
        let mass_kg = [primary.mass_kg(), opposing.mass_kg()];
        let direction = (transfer.transferred_mass_kg > 0.0).then_some(if forward {
            ValveDirection::PrimaryToOpposing
        } else {
            ValveDirection::OpposingToPrimary
        });
        let mut next = self.clone();
        next.primary.gas = primary;
        next.opposing = opposing;
        poll(cancelled)?;
        Ok((
            next,
            ValveTransfer {
                opening,
                direction,
                flow,
                transfer,
                temperature_k,
                mass_kg,
            },
        ))
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn finish_valve(
        &self,
        geometry: RectangularMeteringValve,
        input: DoublePistonInput,
        exchange: GasNozzleStepSolve,
        total_energy_tolerance_j: f64,
        transfer: ValveTransfer,
        mut next: Self,
        piston: DoublePistonStep,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DoublePistonValveCandidate, PistonError> {
        poll(cancelled)?;
        if !positive(total_energy_tolerance_j) {
            return Err(invalid(
                "positive finite metering-valve total energy tolerance required",
            ));
        }
        let endpoint_opening = geometry.opening(&next)?;
        let gases = [next.primary.gas(), &next.opposing];
        let pressures = gases.map(|g| g.state().state().pressure);
        if matches!(transfer.direction, Some(ValveDirection::PrimaryToOpposing))
            && pressures[0] < pressures[1]
            || matches!(transfer.direction, Some(ValveDirection::OpposingToPrimary))
                && pressures[1] < pressures[0]
        {
            return Err(invalid(
                "metered transfer crosses pressure equilibrium after motion/heat; reduce the physical step",
            ));
        }
        let mut motion_energy = [0.0; 2];
        for i in 0..2 {
            poll(cancelled)?;
            if gases[i].mass_kg().to_bits() != transfer.mass_kg[i].to_bits() {
                return Err(invalid(
                    "native closed wall motion changed the transferred gas inventory",
                ));
            }
            // The original mechanical owner preserves the post-transfer model
            // and mass. Re-evaluate its ACTUAL endpoint caloric increment, not
            // a copied first-law residual or a subtraction of absolute energies.
            motion_energy[i] = gases[i].mass_kg()
                * gases[i]
                    .model()
                    .thermodynamics()
                    .specific_internal_energy_change(
                        Temperature::new(transfer.temperature_k[i]),
                        Temperature::new(gases[i].state().state().temperature),
                    )
                    .map_err(|error| {
                        GasEnergyError::Material(fs_material::MaterialError::State {
                            what: format!("metering-valve endpoint caloric increment: {error}"),
                        })
                    })?
                    .value();
        }
        let transferred = if transfer.direction == Some(ValveDirection::OpposingToPrimary) {
            [
                transfer.transfer.receiver_energy_change_j,
                transfer.transfer.donor_energy_change_j,
            ]
        } else {
            [
                transfer.transfer.donor_energy_change_j,
                transfer.transfer.receiver_energy_change_j,
            ]
        };
        let gas_energy_change_j = std::array::from_fn(|i| transferred[i] + motion_energy[i]);
        let mass_residual = (gases[0].mass_kg() - self.primary.gas().mass_kg())
            + (gases[1].mass_kg() - self.opposing.mass_kg());
        let velocities = [self.primary.velocity_m_s(), next.primary.velocity_m_s()];
        let kinetic = (0.5 * self.primary.mass_kg() * (velocities[1] - velocities[0]))
            * (velocities[1] + velocities[0]);
        let displacement =
            (gases[0].volume_m3() - self.primary.gas().volume_m3()) / self.primary.area_m2();
        let force_work = input.external_force_n * displacement;
        // The paired transport is internal: its reference-invariant native
        // residual is its whole-pair storage change. Summing two displayed
        // formation-energy totals here would destroy a small thermal signal.
        let gas_change =
            (motion_energy[0] + motion_energy[1]) + transfer.transfer.energy_residual_j;
        let residual = ((gas_change + kinetic) - input.heat_j[0] - input.heat_j[1]) - force_work;
        if !motion_energy
            .iter()
            .chain(&gas_energy_change_j)
            .chain(&[
                mass_residual,
                kinetic,
                displacement,
                force_work,
                gas_change,
                residual,
            ])
            .all(|value| value.is_finite())
            || mass_residual.abs() > exchange.transfer.mass_tolerance_kg
            || residual.abs() > total_energy_tolerance_j
        {
            return Err(invalid(
                "metered two-gas endpoints do not close whole-system mass and energy",
            ));
        }
        let report = DoublePistonValveStep {
            opening: transfer.opening,
            endpoint_opening,
            direction: transfer.direction,
            flow: transfer.flow,
            transfer: transfer.transfer,
            piston,
            gas_energy_change_j,
            total_mass_residual_kg: mass_residual,
            total_energy_residual_j: residual,
        };
        if report.direction.is_some() {
            self.record_valve_account(
                &mut next,
                input,
                &report,
                gas_change + kinetic,
                kinetic,
                force_work,
                total_energy_tolerance_j,
            );
        }
        poll(cancelled)?;
        Ok(DoublePistonValveCandidate {
            piston: next,
            step: report,
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn record_valve_account(
        &self,
        next: &mut Self,
        input: DoublePistonInput,
        report: &DoublePistonValveStep,
        stored: f64,
        kinetic: f64,
        force_work: f64,
        tolerance_j: f64,
    ) {
        use crate::port_adoption::{PortTransfer, StepAccount};
        let work = report.piston.gas.each_ref().map(|gas| gas.work_by_gas_j);
        let temperature = report
            .piston
            .midpoint
            .each_ref()
            .map(|gas| gas.state().temperature);
        // Preserve actual mechanical/heat observations and the original tick.
        // The native mass-transfer report has no mixing-entropy account; the
        // closed typed audit must not silently omit that part of this interval.
        next.ports.accept(StepAccount { tick: self.ports.next_tick(), duration_s: input.duration_s,
            storage_change_j: Some(stored), loss_j: Some(0.0), included_work_j: force_work,
            entropy_change_j_k: 0.0, tolerance_j: Some(tolerance_j),
            unavailable: Some("metered gas transfer has no native mixing-entropy account; use the paired valve energy/species report"),
        }, [PortTransfer::work(work[0], false), PortTransfer::work(work[1], false),
            PortTransfer::work(-(work[0] + work[1]), false),
            PortTransfer::heat(input.heat_j[0], temperature[0]),
            PortTransfer::heat(input.heat_j[1], temperature[1])]);
        next.primary.ports.accept(StepAccount { tick: self.primary.ports.next_tick(),
            duration_s: input.duration_s,
            storage_change_j: Some(report.gas_energy_change_j[0] + kinetic),
            loss_j: Some(0.0), included_work_j: 0.0, entropy_change_j_k: 0.0,
            tolerance_j: Some(tolerance_j),
            unavailable: Some("metered primary gas has an open species/enthalpy stream; use the paired valve report"),
        }, [PortTransfer::work(force_work + work[1], true),
            PortTransfer::heat(input.heat_j[0], temperature[0])]);
    }

    /// Preview bidirectional compressible metering and the original shared wall.
    ///
    /// Geometry and both gases are sampled before transfer. The native nozzle
    /// selects unchoked/choked flow from source caloric/transport data and moves
    /// species with donor TOTAL enthalpy once. `exchange` retains its explicit
    /// donor/receiver solve roles; these roles follow the actual pressure order.
    /// Heat and applied force then enter the original two-gas midpoint solve.
    /// The final pressure-order gate includes that heating and motion.
    ///
    /// This is a first-order frozen-opening split. It supplies no jet recoil,
    /// viscous valve calibration, metering-edge friction, contact stop, resolved
    /// valve fluid storage or acoustic junction law. In particular, donor jet
    /// energy and flow work are already in enthalpy and are not credited again.
    /// Nonzero transfer retains actual port observations but the full typed
    /// stream/entropy audit remains unavailable; this report closes pair energy
    /// and native species balances. Closed transport retains the native account.
    ///
    /// All source, inventory, geometry, clock, pressure-order and conservation
    /// gates precede publication. A late refusal or cancellation returns neither
    /// gas nor wall endpoint. A caller enclosing finite counterparts must commit
    /// this candidate with those counterparts, not independently.
    #[allow(clippy::too_many_arguments)]
    pub fn valve_candidate(
        &self,
        geometry: RectangularMeteringValve,
        input: DoublePistonInput,
        exchange: GasNozzleStepSolve,
        controls: DoublePistonSolve,
        total_energy_tolerance_j: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DoublePistonValveCandidate, PistonError> {
        poll(&mut cancelled)?;
        if !positive(total_energy_tolerance_j) {
            return Err(invalid(
                "positive finite metering-valve total energy tolerance required",
            ));
        }
        let (mut next, transfer) =
            self.prepare_valve(geometry, input.duration_s, exchange, &mut cancelled)?;
        let piston = next.step(input, controls, &mut cancelled)?;
        self.finish_valve(
            geometry,
            input,
            exchange,
            total_energy_tolerance_j,
            transfer,
            next,
            piston,
            &mut cancelled,
        )
    }

    /// Atomically publish both inventories, their actual volume change and the
    /// one wall momentum after the complete candidate and final cancellation.
    #[allow(clippy::too_many_arguments)]
    pub fn step_with_valve(
        &mut self,
        geometry: RectangularMeteringValve,
        input: DoublePistonInput,
        exchange: GasNozzleStepSolve,
        controls: DoublePistonSolve,
        total_energy_tolerance_j: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<DoublePistonValveStep, PistonError> {
        let candidate = self.valve_candidate(
            geometry,
            input,
            exchange,
            controls,
            total_energy_tolerance_j,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        let (next, report) = candidate.into_parts();
        *self = next;
        Ok(report)
    }
}

#[cfg(test)]
mod tests;
