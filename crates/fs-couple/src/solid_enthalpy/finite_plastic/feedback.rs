//! Shared finite thermoplastic feedback on the original thermal timeline.
//!
//! Each mechanical callback proposes from its same accepted state and interval.
//! The native thermal owner retains storage, transport, events and observer
//! history; only the original material's signed caloric heat is added once.

use std::collections::BTreeMap;

use crate::iqn_ils::{IqnIls, driver::CouplingMethod};
use fs_exec::Cx;
use fs_material::phase::{EquilibriumPhaseState, PhaseStateError, ResolvedHeatCapacityCurve};
use fs_solid::hyper3::plastic::FinitePlasticLawField;
use fs_solid::hyper3::plastic::thermal::FiniteThermoplasticTetEnergy;
use fs_solid::plastic3::dynamics::InertialEnergy;
use fs_thermal::enthalpy::{
    CellId, DepositedCell, DepositionError, DepositionObserver, DepositionSimulation, Event,
    HeatTransfer, StepConfig, ThermalLaw,
};

use super::{
    FiniteThermoplasticConfig, FiniteThermoplasticError, FiniteThermoplasticMaterial,
    FiniteThermoplasticReport, finite, poll, positive,
};

/// Caloric ownership shared by native finite constitutive families. Mechanical
/// history remains entirely inside each original solid proposal.
pub(super) trait CaloricMaterial {
    fn caloric_source(&self) -> &ResolvedHeatCapacityCurve;
    fn state_at_specific_enthalpy(&self, h: f64) -> Result<EquilibriumPhaseState, PhaseStateError> {
        self.caloric_source().curve().state_at_specific_enthalpy(h)
    }
}

impl CaloricMaterial for FiniteThermoplasticMaterial<'_> {
    fn caloric_source(&self) -> &ResolvedHeatCapacityCurve {
        FiniteThermoplasticMaterial::caloric_source(self)
    }
}

impl CaloricMaterial for super::porous::material::FiniteThermalPorousMaterialBinding<'_> {
    fn caloric_source(&self) -> &ResolvedHeatCapacityCurve {
        super::porous::material::FiniteThermalPorousMaterialBinding::caloric_source(self)
    }
}

/// Exact original cell/material pairing, independent of mechanical dimension.
pub(super) struct ThermalBinding<'b, 'm, M> {
    pub(super) cells: &'b [DepositedCell],
    pub(super) materials: FinitePlasticLawField<'m, M>,
}

impl<M> Copy for ThermalBinding<'_, '_, M> {}
impl<M> Clone for ThermalBinding<'_, '_, M> {
    fn clone(&self) -> Self {
        *self
    }
}

impl<M: CaloricMaterial> ThermalBinding<'_, '_, M> {
    fn validate_count<E>(&self) -> Result<(), FiniteThermoplasticError<E>> {
        if let FinitePlasticLawField::PerElement(materials) = self.materials
            && materials.len() != self.cells.len()
        {
            return Err(FiniteThermoplasticError::Invalid(
                "exactly one original finite material per structural cell required",
            ));
        }
        Ok(())
    }

    fn material(&self, element: usize) -> &M {
        match self.materials {
            FinitePlasticLawField::Uniform(material) => material,
            FinitePlasticLawField::PerElement(materials) => &materials[element],
        }
    }

    pub(super) fn enthalpies<E>(
        &self,
        thermal: &DepositionSimulation<'_>,
        accepted: bool,
        cx: &Cx<'_>,
    ) -> Result<Vec<f64>, FiniteThermoplasticError<E>> {
        self.validate_count()?;
        let mut result = Vec::with_capacity(self.cells.len());
        for (i, old) in self.cells.iter().enumerate() {
            poll(cx)?;
            let material = self.material(i);
            let cell = thermal
                .cells()
                .get(&old.id)
                .ok_or(FiniteThermoplasticError::Invalid(
                    "finite structural cell became inactive",
                ))?;
            if cell.id != old.id
                || cell.segment_id != old.segment_id
                || cell.body_id != old.body_id
                || cell.material_id != old.material_id
                || cell.geometry != old.geometry
                || cell.mass_kg.to_bits() != old.mass_kg.to_bits()
                || accepted && cell.specific_h_j_kg.to_bits() != old.specific_h_j_kg.to_bits()
                || !thermal.uses_material_law(
                    cell.material_id,
                    material.caloric_source().curve() as &dyn ThermalLaw,
                )
            {
                return Err(FiniteThermoplasticError::Invalid(
                    "finite structural identity, reference geometry, mass or accepted caloric state changed",
                ));
            }
            material
                .state_at_specific_enthalpy(cell.specific_h_j_kg)
                .map_err(FiniteThermoplasticError::Caloric)?;
            result.push(cell.specific_h_j_kg);
        }
        Ok(result)
    }

    pub(super) fn temperatures<E>(
        &self,
        enthalpies: &[f64],
    ) -> Result<Vec<f64>, FiniteThermoplasticError<E>> {
        self.validate_count()?;
        if enthalpies.len() != self.cells.len() {
            return Err(FiniteThermoplasticError::Invalid(
                "exactly one finite caloric enthalpy per structural cell required",
            ));
        }
        enthalpies
            .iter()
            .enumerate()
            .map(|(i, h)| {
                self.material(i)
                    .state_at_specific_enthalpy(*h)
                    .map(|s| s.temperature_k())
                    .map_err(FiniteThermoplasticError::Caloric)
            })
            .collect()
    }
}

/// Original mechanics proposal plus already integrated physical cell heat.
pub(super) struct MechanicalProposal<S, E> {
    pub(super) state: S,
    pub(super) equilibrium: E,
    pub(super) energy: FiniteThermoplasticTetEnergy,
    pub(super) inertia: Option<InertialEnergy>,
    pub(super) heat_j: Vec<f64>,
}

/// Complete accepted thermal/mechanical/observer endpoint.
pub(super) struct FeedbackCandidate<'t, S, E, O> {
    pub(super) thermal: DepositionSimulation<'t>,
    pub(super) state: S,
    pub(super) observer: O,
    pub(super) cells: Vec<DepositedCell>,
    pub(super) report: FiniteThermoplasticReport<E>,
}

fn admit_heat<E>(heat_j: &[f64], cells: usize) -> Result<(), FiniteThermoplasticError<E>> {
    if heat_j.len() != cells || heat_j.iter().any(|q| !q.is_finite()) {
        return Err(FiniteThermoplasticError::Invalid(
            "exactly one finite integrated caloric heat per structural cell required",
        ));
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
pub(super) fn advance<'t, M: CaloricMaterial, S, E, O: DepositionObserver>(
    binding: ThermalBinding<'_, '_, M>,
    time_s: f64,
    thermal: &DepositionSimulation<'t>,
    observer: &O,
    target_s: f64,
    thermal_config: StepConfig,
    config: FiniteThermoplasticConfig,
    cx: &Cx<'_>,
    mechanics: impl FnMut(
        &[f64],
        f64,
    ) -> Result<MechanicalProposal<S, E>, FiniteThermoplasticError<O::Error>>,
) -> Result<FeedbackCandidate<'t, S, E, O>, FiniteThermoplasticError<O::Error>> {
    advance_with_coupling(
        binding,
        time_s,
        thermal,
        observer,
        target_s,
        thermal_config,
        config,
        CouplingMethod::RelaxedPicard,
        cx,
        mechanics,
    )
}

#[allow(clippy::too_many_arguments)]
pub(super) fn advance_with_coupling<'t, M: CaloricMaterial, S, E, O: DepositionObserver>(
    binding: ThermalBinding<'_, '_, M>,
    time_s: f64,
    thermal: &DepositionSimulation<'t>,
    observer: &O,
    target_s: f64,
    thermal_config: StepConfig,
    config: FiniteThermoplasticConfig,
    method: CouplingMethod,
    cx: &Cx<'_>,
    mut mechanics: impl FnMut(
        &[f64],
        f64,
    )
        -> Result<MechanicalProposal<S, E>, FiniteThermoplasticError<O::Error>>,
) -> Result<FeedbackCandidate<'t, S, E, O>, FiniteThermoplasticError<O::Error>> {
    poll(cx)?;
    let duration_s = target_s - time_s;
    if thermal.time_s() != time_s
        || !positive(duration_s)
        || !target_s.is_finite()
        || config.maximum_iterations == 0
        || config.maximum_iterations > 256
        || !positive(config.temperature_tolerance_k)
        || !positive(config.energy_tolerance_j)
        || !positive(config.enthalpy_relaxation)
        || config.enthalpy_relaxation > 1.0
        || !thermal_config.energy_tolerance_j.is_finite()
        || thermal_config.energy_tolerance_j < 0.0
        || thermal_config.energy_tolerance_j > config.energy_tolerance_j
    {
        return Err(FiniteThermoplasticError::Invalid(
            "matching clocks, positive interval and bounded independent thermal coupling controls required",
        ));
    }
    let mut enthalpies = binding.enthalpies(thermal, true, cx)?;
    // The thermal/mechanical map changes between physical intervals; never
    // carry its interface secants into the next accepted material history.
    let mut accelerator = match method {
        CouplingMethod::RelaxedPicard => None,
        CouplingMethod::IqnIls(policy) => Some(
            IqnIls::new(enthalpies.len(), policy).map_err(FiniteThermoplasticError::Accelerator)?,
        ),
    };
    let mut defects = [f64::INFINITY; 4];
    // Every proposal and cached endpoint starts from the SAME accepted owner.
    // Reuse an actual-temperature solve only at its exact evaluated enthalpies;
    // an accelerated iterate normally needs its own native equilibrium solve.
    let mut next_trial = None;
    for iteration in 1..=config.maximum_iterations {
        poll(cx)?;
        let trial_temperatures = binding.temperatures(&enthalpies)?;
        let trial = if let Some(cached) = next_trial.take() {
            cached
        } else {
            mechanics(&trial_temperatures, duration_s)?
        };
        admit_heat(&trial.heat_j, binding.cells.len())?;
        let applied_heat_j = finite(trial.heat_j.iter().sum())?;
        let heat = FiniteHeat {
            binding,
            inner: observer.clone(),
            cx,
            watts: binding
                .cells
                .iter()
                .zip(&trial.heat_j)
                .map(|(cell, q)| Ok((cell.id, finite(*q / duration_s)?)))
                .collect::<Result<_, FiniteThermoplasticError<O::Error>>>()?,
        };
        let (next_thermal, heat, thermal_audit) =
            thermal.advance_observed(&heat, target_s, thermal_config, || {
                cx.checkpoint().map_err(|_| DepositionError::Cancelled)
            })?;
        let final_h = binding.enthalpies(&next_thermal, false, cx)?;
        let final_temperatures = binding.temperatures(&final_h)?;
        let (actual, heat_defect_j) = if final_temperatures == trial_temperatures {
            (trial, 0.0)
        } else {
            let actual = mechanics(&final_temperatures, duration_s)?;
            admit_heat(&actual.heat_j, binding.cells.len())?;
            let defect = finite(
                trial
                    .heat_j
                    .iter()
                    .zip(&actual.heat_j)
                    .map(|(a, b)| (a - b).abs())
                    .sum(),
            )?;
            (actual, defect)
        };
        let caloric_heat_j = finite(actual.heat_j.iter().sum())?;
        let external_heat_j = finite(thermal_audit.external_heat_j - applied_heat_j)?;
        let thermal_residual_j = finite(
            thermal_audit.stored_energy_change_j
                - thermal_audit.incoming_energy_j
                - external_heat_j
                - caloric_heat_j,
        )?;
        let inertia_account = actual.inertia.as_ref().map_or(0.0, |inertia| {
            inertia.kinetic_change_j + inertia.numerical_loss_j
        });
        let combined_residual_j = finite(
            actual.energy.internal_energy_change_j
                + thermal_audit.stored_energy_change_j
                + actual.energy.endpoint_work_defect_j
                + inertia_account
                - actual.energy.dead_load_endpoint_work_j
                - actual.energy.support_endpoint_work_j
                - external_heat_j
                - thermal_audit.incoming_energy_j,
        )?;
        defects = [
            trial_temperatures
                .iter()
                .zip(&final_temperatures)
                .map(|(a, b)| (a - b).abs())
                .fold(0.0, f64::max),
            heat_defect_j,
            thermal_residual_j,
            combined_residual_j,
        ];
        if defects[0] <= config.temperature_tolerance_k
            && defects[1] <= config.energy_tolerance_j
            && defects[2].abs() <= config.energy_tolerance_j
            && defects[3].abs() <= config.energy_tolerance_j
        {
            let cells = binding
                .cells
                .iter()
                .map(|cell| next_thermal.cells()[&cell.id].clone())
                .collect();
            let report = FiniteThermoplasticReport {
                iterations: iteration,
                equilibrium: actual.equilibrium,
                material_energy: actual.energy,
                inertia: actual.inertia,
                thermal_audit,
                applied_heat_j,
                external_heat_j,
                temperature_defect_k: defects[0],
                heat_defect_j: defects[1],
                thermal_residual_j,
                combined_residual_j,
            };
            poll(cx)?;
            return Ok(FeedbackCandidate {
                thermal: next_thermal,
                state: actual.state,
                observer: heat.inner,
                cells,
                report,
            });
        }
        if iteration == config.maximum_iterations {
            break;
        }
        poll(cx)?;
        if let Some(accelerator) = &mut accelerator {
            let proposal = accelerator
                .step(&enthalpies, &final_h, config.enthalpy_relaxation)
                .map_err(FiniteThermoplasticError::Accelerator)?;
            poll(cx)?;
            if proposal
                .values
                .iter()
                .zip(&final_h)
                .all(|(proposed, evaluated)| proposed.to_bits() == evaluated.to_bits())
            {
                next_trial = Some(actual);
            }
            enthalpies = proposal.values;
        } else if config.enthalpy_relaxation == 1.0 {
            enthalpies = final_h;
            next_trial = Some(actual);
        } else {
            for (h, new_h) in enthalpies.iter_mut().zip(final_h) {
                *h = finite(
                    (1.0 - config.enthalpy_relaxation) * *h + config.enthalpy_relaxation * new_h,
                )?;
            }
        }
    }
    poll(cx)?;
    Err(FiniteThermoplasticError::NotConverged {
        iterations: config.maximum_iterations,
        temperature_defect_k: defects[0],
        heat_defect_j: defects[1],
        thermal_residual_j: defects[2],
        combined_residual_j: defects[3],
    })
}

struct FiniteHeat<'b, 'm, 'cx, M, O> {
    binding: ThermalBinding<'b, 'm, M>,
    inner: O,
    cx: &'b Cx<'cx>,
    watts: BTreeMap<CellId, f64>,
}
impl<M, O: Clone> Clone for FiniteHeat<'_, '_, '_, M, O> {
    fn clone(&self) -> Self {
        Self {
            binding: self.binding,
            inner: self.inner.clone(),
            cx: self.cx,
            watts: self.watts.clone(),
        }
    }
}
impl<M: CaloricMaterial, O: DepositionObserver> DepositionObserver
    for FiniteHeat<'_, '_, '_, M, O>
{
    type Error = FiniteThermoplasticError<O::Error>;
    fn begin(&mut self) -> Result<(), Self::Error> {
        self.inner
            .begin()
            .map_err(FiniteThermoplasticError::Observer)
    }
    fn requires_step_interval(&self) -> bool {
        self.inner.requires_step_interval()
    }
    fn next_transition_time_s(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Option<f64>, Self::Error> {
        self.inner
            .next_transition_time_s(thermal, checkpoint)
            .map_err(FiniteThermoplasticError::Observer)
    }
    fn sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.binding.enthalpies(thermal, false, self.cx)?;
        let mut rates = self
            .inner
            .sources(thermal, checkpoint)
            .map_err(FiniteThermoplasticError::Observer)?;
        for (&cell, &watts) in &self.watts {
            checkpoint()?;
            let rate = rates.entry(cell).or_insert(0.0);
            *rate = finite(*rate + watts)?;
        }
        Ok(rates)
    }
    fn step_sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        start: f64,
        dt: f64,
        end: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.inner
            .step_sources(thermal, start, dt, end, checkpoint)
            .map_err(FiniteThermoplasticError::Observer)
    }
    fn conductances(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<u64, f64>, Self::Error> {
        self.inner
            .conductances(thermal, checkpoint)
            .map_err(FiniteThermoplasticError::Observer)
    }
    fn reservoir_conductances(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.inner
            .reservoir_conductances(thermal, checkpoint)
            .map_err(FiniteThermoplasticError::Observer)
    }
    fn ambient_conductances(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.inner
            .ambient_conductances(thermal, checkpoint)
            .map_err(FiniteThermoplasticError::Observer)
    }
    fn heat_transfers(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<HeatTransfer>, Self::Error> {
        self.inner
            .heat_transfers(thermal, checkpoint)
            .map_err(FiniteThermoplasticError::Observer)
    }
    fn limit_timestep_s(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        dt: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<f64, Self::Error> {
        self.inner
            .limit_timestep_s(thermal, dt, checkpoint)
            .map_err(FiniteThermoplasticError::Observer)
    }
    fn event(
        &mut self,
        event: &Event,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<(), Self::Error> {
        self.inner
            .event(event, thermal)
            .map_err(FiniteThermoplasticError::Observer)?;
        self.binding.enthalpies(thermal, false, self.cx).map(|_| ())
    }
    fn step(
        &mut self,
        start: f64,
        duration: f64,
        old: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.inner
            .step(start, duration, old, thermal, checkpoint)
            .map_err(FiniteThermoplasticError::Observer)?;
        self.binding.enthalpies(thermal, false, self.cx).map(|_| ())
    }
}
