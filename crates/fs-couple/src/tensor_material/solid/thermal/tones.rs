//! Continuous nodal forces, retained tensor memory and finite-cell self-heating.
//!
//! The existing prepared exact flow owns all mechanical startup motion,
//! interference, inertia and positive Kelvin/Maxwell loss. Its one physical
//! interval is evaluated once. The original caloric timeline distributes that
//! interval's mean heat power over its own event/CFL/observer substeps.

use super::{
    ElementThermalCell, Receiver, Result, SolidThermalConfig, SolidThermalError,
    balance, check, number,
};
use super::spatial::{deposit_element_heat, same_inventory};
use crate::tensor_material::prepared::PreparedTensorError;
use crate::tensor_material::solid::prepared::{
    PreparedSolidError, PreparedTensorSolidBody, PreparedTensorSolidHistory,
    PreparedTensorSolidToneStep,
};
use fs_exec::CancelGate;
use fs_math::dd::Dd;
use fs_solid::linear3::{TetElasticError, visco::ViscoElementError};
use fs_thermal::enthalpy::{DepositionObserver, DepositionSimulation, EnergyAudit};
use std::collections::BTreeSet;

/// The original continuous-force report and its independently audited heat handoff.
#[derive(Clone, Debug)]
pub struct PreparedToneThermalReport {
    /// Full original nodal force moments, kinetic/elastic/memory storage, work
    /// and positive element loss. Fixed supports do no mechanical work.
    pub mechanical: PreparedTensorSolidToneStep,
    pub thermal: EnergyAudit,
    /// Physical heat in original element order [J].
    pub element_heat_j: Vec<f64>,
    pub material_heat_j: f64,
    pub deposited_heat_j: f64,
    /// Native external heat minus this INTERNAL mechanical heat transfer [J].
    /// Incoming material energy remains separately recorded in `thermal`.
    pub other_external_heat_j: f64,
    /// Delta(kinetic + elastic + Maxwell storage) + heat - force work [J].
    pub mechanical_energy_residual_j: f64,
    pub caloric_energy_residual_j: f64,
    pub heat_handoff_residual_j: f64,
    pub combined_energy_residual_j: f64,
}

/// A complete unpublished mechanical/thermal endpoint and observer history.
pub struct PreparedToneThermalCandidate<'a, 's, 't, O> {
    solid: PreparedTensorSolidHistory<'a, 's>,
    thermal: DepositionSimulation<'t>,
    observer: O,
    report: PreparedToneThermalReport,
}

impl<'a, 's, 't, O> PreparedToneThermalCandidate<'a, 's, 't, O> {
    pub fn solid(&self) -> &PreparedTensorSolidHistory<'a, 's> {
        &self.solid
    }
    pub fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub fn observer(&self) -> &O {
        &self.observer
    }
    pub fn report(&self) -> &PreparedToneThermalReport {
        &self.report
    }
    pub fn into_parts(self) -> (
        PreparedTensorSolidHistory<'a, 's>,
        DepositionSimulation<'t>,
        O,
        PreparedToneThermalReport,
    ) {
        (self.solid, self.thermal, self.observer, self.report)
    }
}

fn mechanical<E>(error: PreparedSolidError) -> SolidThermalError<E> {
    match error {
        PreparedSolidError::Cancelled
        | PreparedSolidError::Dynamics(PreparedTensorError::Cancelled)
        | PreparedSolidError::Source(ViscoElementError::Elastic(TetElasticError::Cancelled)) => {
            SolidThermalError::Cancelled
        }
        other => SolidThermalError::Prepared(other),
    }
}

impl<'s> PreparedTensorSolidBody<'s> {
    /// Advance all original continuous nodal-force tones and deposit their
    /// physical element heat on the original finite thermal timeline.
    ///
    /// This body must be the accepted history's exact original prepared owner.
    /// Each cosine/sine field follows original full node-major xyz coordinates;
    /// coefficients are relative to this interval's start. Use the existing
    /// `advance_force_phase` when continuing an unchanged physical drive.
    ///
    /// Motion, Maxwell memory and original work/loss caches advance only once,
    /// independently of caloric event or CFL cuts. No sum of separate-tone
    /// heats, settled response or prescribed-displacement substitute is used.
    /// The supplied basis retains its original full/reduced approximation scope.
    ///
    /// Every original element needs a unique whole-reference-mass/volume cell
    /// and its exact current caloric law. The explicit fully-solid temperature
    /// interval admits FROZEN mechanical coefficients; temperature feedback,
    /// thermal expansion and moving supports are outside this transaction.
    /// Original endpoint strain/frequency admission remains unchanged and does
    /// not certify the interior trajectory or its startup spectral bandwidth.
    ///
    /// The interval's mean physical heat enters the existing explicit enthalpy
    /// owner. Refine the prepared mechanical duration as well as thermal cuts
    /// for heating/transport accuracy. A small energy residual alone does not
    /// bound that temporal approximation. The original positive fixed clock,
    /// all prepared work/storage bounds and caloric substep limits remain.
    ///
    /// The observer's existing source, transport, interval, event and timestep
    /// hooks stay active; its Clone must isolate mutation as required by the
    /// native caloric API. A late refusal or cancellation publishes neither
    /// half-state and never consumes accepted material heat.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_tones_thermal_candidate<'a, 't, O: DepositionObserver>(
        &'a self,
        accepted: &PreparedTensorSolidHistory<'a, 's>,
        thermal: &DepositionSimulation<'t>,
        observer: &O,
        cosine_nodal_force_n: &[Vec<[f64; 3]>],
        sine_nodal_force_n: &[Vec<[f64; 3]>],
        cells: &[ElementThermalCell<'_>],
        config: SolidThermalConfig,
        gate: &CancelGate,
    ) -> Result<PreparedToneThermalCandidate<'a, 's, 't, O>, O::Error> {
        check(gate)?;
        let source = self.projection().source();
        if !core::ptr::eq(accepted.material_history().body(), self.material_body())
            || self.angular_frequencies_rad_s().is_none()
            || accepted.material_history().time_s().to_bits() != thermal.time_s().to_bits()
            || cells.len() != source.element_count()
            || !config.energy_tolerance_j.is_finite()
            || config.energy_tolerance_j <= 0.0
            || !config.thermal.energy_tolerance_j.is_finite()
            || config.thermal.energy_tolerance_j < 0.0
            || config.thermal.energy_tolerance_j > config.energy_tolerance_j
        {
            return Err(SolidThermalError::Invalid(
                "original tone owner, one shared clock, complete receivers and bounded energy controls required",
            ));
        }
        let mut used = BTreeSet::new();
        let mut receivers = Vec::with_capacity(cells.len());
        for (element, &binding) in cells.iter().enumerate() {
            check(gate)?;
            let [low, high] = binding.temperature_interval_k;
            if !used.insert(binding.cell)
                || !low.is_finite()
                || !high.is_finite()
                || low <= 0.0
                || high < low
            {
                return Err(SolidThermalError::Invalid(
                    "unique caloric receivers with explicit positive temperature intervals required",
                ));
            }
            let original = thermal.cells().get(&binding.cell)
                .ok_or(SolidThermalError::Invalid(
                    "inactive prepared-tone material-heat destination",
                ))?.clone();
            let mass = source.element_reference_mass_kg(element)
                .ok_or(SolidThermalError::Invalid("missing original element mass"))?;
            let volume = source.element_reference_volume_m3(element)
                .ok_or(SolidThermalError::Invalid("missing original element volume"))?;
            if !same_inventory(original.mass_kg, mass)
                || !same_inventory(original.geometry.volume_m3, volume)
            {
                return Err(SolidThermalError::Invalid(
                    "each caloric cell must retain its whole original element reference mass and volume",
                ));
            }
            let receiver = Receiver { binding, original };
            receiver.admit(thermal)?;
            receivers.push(receiver);
        }

        // A clone borrows the same immutable source/caches and preserves the
        // accepted physical history until every caloric hook and audit succeeds.
        let mut next = accepted.clone();
        let step = next.advance_tones(
            cosine_nodal_force_n, sine_nodal_force_n, gate,
        ).map_err(mechanical)?;
        let target = next.material_history().time_s();
        let duration = number(target - thermal.time_s())?;
        if duration <= 0.0 || step.motion.element_heat_j().len() != receivers.len() {
            return Err(SolidThermalError::Invalid(
                "prepared motion changed its element loss inventory or endpoint clock",
            ));
        }
        let element_heat_j = step.motion.element_heat_j().to_vec();
        let mut sum = Dd::ZERO;
        for &heat in &element_heat_j {
            check(gate)?;
            let heat = number(heat)?;
            if heat < 0.0 {
                return Err(SolidThermalError::Invalid("negative original material heat"));
            }
            sum = sum + Dd::from_f64(heat);
        }
        let energy = step.motion.material.energy;
        let material_heat_j = number(energy.dissipated_j)?;
        if material_heat_j < 0.0 {
            return Err(SolidThermalError::Invalid("negative total original material heat"));
        }
        balance((sum - Dd::from_f64(material_heat_j)).to_f64(), config.energy_tolerance_j)?;
        // The original owner computes this from stable coordinate increments,
        // not subtraction of rounded total storage or an inferred heat source.
        let mechanical_energy_residual_j = number(-energy.balance_error_j)?;
        balance(mechanical_energy_residual_j, config.energy_tolerance_j)?;

        let (next_thermal, next_observer, audit, deposited_heat_j) = deposit_element_heat(
            thermal, observer, receivers, &element_heat_j, target, config,
            || gate.is_requested(),
        )?;
        let deposited = Dd::from_f64(deposited_heat_j);
        let handoff = deposited - Dd::from_f64(material_heat_j);
        let caloric = Dd::from_f64(audit.stored_energy_change_j)
            - Dd::from_f64(audit.incoming_energy_j)
            - Dd::from_f64(audit.external_heat_j);
        let other_external_heat_j =
            number((Dd::from_f64(audit.external_heat_j) - deposited).to_f64())?;
        let heat_handoff_residual_j = number(handoff.to_f64())?;
        let caloric_energy_residual_j = number(caloric.to_f64())?;
        // Rm = Delta H + D - W and Rt = Delta U - incoming - external.
        // The material transfer is internal to the combined physical system.
        let combined_energy_residual_j = number(
            (caloric + handoff + Dd::from_f64(mechanical_energy_residual_j)).to_f64(),
        )?;
        for residual in [
            heat_handoff_residual_j,
            caloric_energy_residual_j,
            combined_energy_residual_j,
        ] {
            balance(residual, config.energy_tolerance_j)?;
        }
        check(gate)?;
        Ok(PreparedToneThermalCandidate {
            solid: next,
            thermal: next_thermal,
            observer: next_observer,
            report: PreparedToneThermalReport {
                mechanical: step,
                thermal: audit,
                element_heat_j,
                material_heat_j,
                deposited_heat_j,
                other_external_heat_j,
                mechanical_energy_residual_j,
                caloric_energy_residual_j,
                heat_handoff_residual_j,
                combined_energy_residual_j,
            },
        })
    }
}
