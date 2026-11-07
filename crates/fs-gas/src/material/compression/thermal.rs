//! Reacting compression with source-owned heat exchange and a finite-rate stroke.
//!
//! The same symmetric volume/rigid-stage composition as adiabatic compression
//! uses a prepared material reactor for simultaneous reaction and heat exchange.
//! Gas NASA calorics, formation energy, composition, transport and acoustics stay
//! in their original owners. A finite body retains its original enthalpy chart,
//! latent inventory, density and isobaric expansion work through every interval.
//!
//! The endpoint methods hold caller-supplied conductance constant over a stroke.
//! [`schedule`] additionally resolves prescribed time-varying bath temperature
//! and conductance. Conductance is total [W/K], not a conductivity or an inferred
//! moving-wall area law.
//! These homogeneous states have no piston inertia, spatial wall gradients or
//! gas phase change. Refine the prescribed substeps to assess splitting error;
//! native chemistry tolerances apply separately to each rigid-stage call.

use super::{
    ClosedGasParcel, CompressionChemicalStep, CompressionTotals, FixedCompressionPath,
    MaterialCompressionError, accumulate, advance_interval_with, poll,
};
use crate::chemistry::heat::{HeatBath, material::FiniteHeatBody};
use crate::material::isochoric::{MaterialIsochoricReactor, MaterialIsochoricStep};
use crate::material::{MaterialKineticsSolve, ThermochemicalGasModel};
use fs_material::thermochemical_gas::energy::entropy::GasEntropySolve;
use fs_qty::{Time, Volume};
use fs_thermochem::kinetics::IdealGasKineticsV1;

/// Adaptive gas, heat and finite-body enthalpy comparisons on the same interval kernel.
pub mod adaptive;
/// Physical-time bath and finite-body conductance programs on the same kernels.
pub mod schedule;

type Error = MaterialCompressionError;

/// Volume inversion and caller-bounded compression/heat/reaction splitting.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialThermalCompressionSolve {
    pub volume: GasEntropySolve,
    /// Equal-time intervals on the prescribed logarithmic volume path, 1..=4096.
    /// A fixed-volume request always uses one original reacting thermal call.
    pub substeps: usize,
}

/// Gas work and signed heat from the accepted original-source operations.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialThermalCompressionStep {
    pub elapsed_time: Time,
    pub initial_volume_m3: f64,
    pub final_volume_m3: f64,
    /// Work delivered to the gas; positive during compression [J].
    pub boundary_work_on_gas_j: f64,
    /// Heat delivered by the bath/body; formation heat is already in U [J].
    pub heat_input_j: f64,
    /// Sum of original-source gas stored-energy increments [J].
    pub internal_energy_change_j: f64,
    /// Delta Ugas - W_on_gas - Qgas [J], accumulated directly from the original
    /// rigid-stage material residuals rather than subtracting large totals.
    pub energy_residual_j: f64,
    /// Original volume-map plus native reactor entropy changes [J/K].
    /// Chemical entries precede material endpoint reconciliation, whose entropy
    /// correction is excluded. This is a source-stage diagnostic, not an
    /// independently bounded entropy difference of the published gas parcel.
    pub gas_entropy_change_j_k: f64,
    /// Actual bath entropy (-Q/Tbath) or original finite-body entropy change.
    pub thermal_partner_entropy_change_j_k: f64,
    /// Original native gas-plus-partner production plus volume-map entropy
    /// increments [J/K]. It has the same material-reconciliation limitation as
    /// gas_entropy_change_j_k; it is not a global accuracy/passivity certificate.
    pub entropy_production_j_k: f64,
    pub substeps: usize,
    pub chemical_accepted_steps: usize,
    pub chemical_rejected_steps: usize,
}

/// Finite source body and compressed gas accepted in one transaction.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialThermalCompressionBodyStep<'m> {
    pub gas: MaterialThermalCompressionStep,
    pub body: FiniteHeatBody<'m>,
    pub body_enthalpy_change_j: f64,
    pub body_internal_energy_change_j: f64,
    /// Original held-body-pressure volume work, positive outward [J].
    /// This is separate from prescribed gas compression work.
    pub body_boundary_work_j: f64,
    pub body_entropy_change_j_k: f64,
    /// Delta Ugas + Delta Hbody - W_on_gas [J].
    pub combined_heat_residual_j: f64,
    /// Delta Ugas + Delta Ubody + W_by_body - W_on_gas [J].
    /// Both combined residuals accumulate the native coupled residuals directly.
    pub combined_energy_residual_j: f64,
}

/// Prepared simultaneous kinetics/thermal stages reused throughout and between
/// strokes. Only numerical workspace changes on failure; physical state does not.
pub struct MaterialThermalCompressionReactor {
    reactor: MaterialIsochoricReactor,
}

impl MaterialThermalCompressionReactor {
    pub fn new(
        bank: ThermochemicalGasModel,
        mechanism: IdealGasKineticsV1,
        controls: MaterialKineticsSolve,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Self, Error> {
        Ok(Self {
            reactor: MaterialIsochoricReactor::new(bank, mechanism, controls, cancelled)?,
        })
    }

    pub const fn bank(&self) -> &ThermochemicalGasModel {
        self.reactor.bank()
    }
    pub const fn mechanism(&self) -> &IdealGasKineticsV1 {
        self.reactor.mechanism()
    }
    pub const fn kinetics_controls(&self) -> MaterialKineticsSolve {
        self.reactor.controls()
    }
    /// Reusable numerical workspace; source banks and endpoint allocations excluded.
    pub fn workspace_bytes(&self) -> usize {
        self.reactor.workspace_bytes()
    }

    /// A held reservoir exchanges G(Tbath-Tgas) while the original mechanism
    /// reacts. Infinite reservoir capacity does not imply infinite conductance.
    /// Any intermediate refusal or cancellation preserves the complete parcel.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_heat_bath(
        &mut self,
        parcel: &mut ClosedGasParcel,
        final_volume: Volume,
        duration: Time,
        bath: HeatBath,
        controls: MaterialThermalCompressionSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MaterialThermalCompressionStep, Error> {
        let (candidate, report) = self.propose(
            parcel,
            final_volume,
            duration,
            controls,
            &mut cancelled,
            |reactor, candidate, duration, cancelled| {
                let step = reactor.advance_with_heat_bath(candidate, duration, bath, cancelled)?;
                let partner_entropy = -step.chemistry.balance().entropy_transfer_j_k;
                Ok(ThermalInterval::from_step(&step, partner_entropy))
            },
        )?;
        poll(&mut cancelled)?;
        *parcel = candidate;
        Ok(report)
    }

    /// Exchange heat with an original finite caloric/phase body while compressing
    /// or expanding the gas. Body mass, source and held pressure remain original.
    /// Its density-derived volume never replaces the prescribed gas volume.
    /// Gas and returned body publish only after the complete stroke succeeds.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_finite_body<'m>(
        &mut self,
        parcel: &mut ClosedGasParcel,
        final_volume: Volume,
        duration: Time,
        mut body: FiniteHeatBody<'m>,
        controls: MaterialThermalCompressionSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MaterialThermalCompressionBodyStep<'m>, Error> {
        let mut body_totals = [(0.0, 0.0); 5];
        let (candidate, gas) = self.propose(
            parcel,
            final_volume,
            duration,
            controls,
            &mut cancelled,
            |reactor, candidate, duration, cancelled| {
                let step =
                    reactor.advance_with_finite_body(candidate, duration, body, cancelled)?;
                for (total, value) in body_totals.iter_mut().zip([
                    step.body_report.enthalpy_change_j,
                    step.body_report.internal_energy_change_j,
                    step.body_report.boundary_work_j,
                    step.combined_heat_residual_j,
                    step.combined_energy_residual_j,
                ]) {
                    accumulate(total, value);
                }
                body = step.body;
                Ok(ThermalInterval::from_step(
                    &step.gas,
                    step.body_entropy_change_j_k,
                ))
            },
        )?;
        if body_totals.iter().any(|value| !value.0.is_finite()) {
            return Err(Error::Invalid(
                "accumulated finite-body energy, work or balance is unrepresentable",
            ));
        }
        let report = MaterialThermalCompressionBodyStep {
            gas,
            body,
            body_enthalpy_change_j: body_totals[0].0,
            body_internal_energy_change_j: body_totals[1].0,
            body_boundary_work_j: body_totals[2].0,
            body_entropy_change_j_k: gas.thermal_partner_entropy_change_j_k,
            combined_heat_residual_j: body_totals[3].0,
            combined_energy_residual_j: body_totals[4].0,
        };
        poll(&mut cancelled)?;
        *parcel = candidate;
        Ok(report)
    }

    // The whole proposed stroke remains private until its public owner has
    // finished every gas/body account and the final publication checkpoint.
    #[allow(clippy::too_many_arguments)]
    fn propose<F: FnMut() -> bool>(
        &mut self,
        parcel: &ClosedGasParcel,
        final_volume: Volume,
        duration: Time,
        controls: MaterialThermalCompressionSolve,
        cancelled: &mut F,
        mut exchange: impl FnMut(
            &mut MaterialIsochoricReactor,
            &mut ClosedGasParcel,
            Time,
            &mut F,
        ) -> Result<ThermalInterval, Error>,
    ) -> Result<(ClosedGasParcel, MaterialThermalCompressionStep), Error> {
        poll(cancelled)?;
        let path = FixedCompressionPath::new(
            parcel.volume_m3(),
            final_volume,
            duration,
            controls.substeps,
        )?;
        let mut candidate = parcel.clone();
        candidate.set_volume_isentropic(
            path.initial_volume_m3,
            controls.volume,
            &mut *cancelled,
        )?;
        let mut gas = CompressionTotals::default();
        let mut heat = (0.0, 0.0);
        let mut partner_entropy = (0.0, 0.0);
        let mut production = (0.0, 0.0);
        for index in 0..path.substeps {
            poll(cancelled)?;
            let volume_entropy = advance_interval_with(
                &mut candidate,
                path.interval(index)?,
                controls.volume,
                &mut gas,
                cancelled,
                |candidate, duration, cancelled| {
                    let step = exchange(&mut self.reactor, candidate, duration, cancelled)?;
                    accumulate(&mut heat, step.heat_input_j);
                    accumulate(&mut partner_entropy, step.partner_entropy_change_j_k);
                    accumulate(&mut production, step.entropy_production_j_k);
                    Ok(step.chemical)
                },
            )?;
            for change in volume_entropy {
                accumulate(&mut production, change);
            }
        }
        let gas = gas.finish(path.initial_volume_m3, path.final_volume_m3, duration)?;
        if [heat.0, partner_entropy.0, production.0]
            .iter()
            .any(|value| !value.is_finite())
        {
            return Err(Error::Invalid(
                "accumulated heat or counterpart entropy is unrepresentable",
            ));
        }
        Ok((
            candidate,
            MaterialThermalCompressionStep {
                elapsed_time: duration,
                initial_volume_m3: path.initial_volume_m3,
                final_volume_m3: path.final_volume_m3,
                boundary_work_on_gas_j: gas.boundary_work_on_gas_j,
                heat_input_j: heat.0,
                internal_energy_change_j: gas.internal_energy_change_j,
                energy_residual_j: gas.energy_residual_j,
                gas_entropy_change_j_k: gas.entropy_change_j_k,
                thermal_partner_entropy_change_j_k: partner_entropy.0,
                entropy_production_j_k: production.0,
                substeps: gas.substeps,
                chemical_accepted_steps: gas.chemical_accepted_steps,
                chemical_rejected_steps: gas.chemical_rejected_steps,
            },
        ))
    }
}

struct ThermalInterval {
    chemical: CompressionChemicalStep,
    heat_input_j: f64,
    partner_entropy_change_j_k: f64,
    entropy_production_j_k: f64,
}

impl ThermalInterval {
    fn from_step(step: &MaterialIsochoricStep<'_>, partner_entropy_change_j_k: f64) -> Self {
        let balance = step.chemistry.balance();
        Self {
            chemical: CompressionChemicalStep {
                internal_energy_change_j: step.material.stored_energy_change_j,
                energy_residual_j: step.material.energy_residual_j,
                entropy_change_j_k: balance.entropy_change_j_k,
                accepted_steps: balance.accepted_steps,
                rejected_steps: balance.rejected_steps,
            },
            heat_input_j: step.material.supplied_heat_j,
            partner_entropy_change_j_k,
            entropy_production_j_k: balance.entropy_production_j_k,
        }
    }
}
