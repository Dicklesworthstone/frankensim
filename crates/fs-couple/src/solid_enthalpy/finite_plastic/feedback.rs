//! Shared finite thermoplastic feedback on the original thermal timeline.
//!
//! Each mechanical callback proposes from its same accepted state and interval.
//! The native thermal owner retains storage, transport, events and observer
//! history; only the original material's signed caloric heat is added once.

use std::collections::BTreeMap;

use crate::iqn_ils::{IqnIls, driver::CouplingMethod};
use fs_exec::Cx;
use fs_material::phase::{EquilibriumPhaseState, PhaseStateError, ResolvedHeatCapacityCurve};
use fs_math::dd::Dd;
use fs_solid::hyper3::plastic::FinitePlasticLawField;
use fs_solid::hyper3::plastic::thermal::FiniteThermoplasticTetEnergy;
use fs_solid::plastic3::dynamics::InertialEnergy;
use fs_thermal::enthalpy::{
    CellId, DepositedCell, DepositionError, DepositionObserver, DepositionSimulation, Event,
    HeatTransfer, StepConfig, ThermalLaw, collect_source_rates,
};

use super::material::ReferenceCaloricMaterial;
use super::{
    FiniteThermoplasticConfig, FiniteThermoplasticError, FiniteThermoplasticMaterialFamily,
    FiniteThermoplasticReport, finite, poll, positive,
};

/// Original caloric storage, independent of the mechanical constitutive family.
/// Every implementation retains its prepared source and exact chart object.
pub(in crate::solid_enthalpy) trait CaloricMaterial {
    fn caloric_source(&self) -> &ResolvedHeatCapacityCurve;

    fn state_at_specific_enthalpy(&self, h: f64) -> Result<EquilibriumPhaseState, PhaseStateError> {
        self.caloric_source().curve().state_at_specific_enthalpy(h)
    }
}

impl<M: FiniteThermoplasticMaterialFamily> CaloricMaterial for M {
    fn caloric_source(&self) -> &ResolvedHeatCapacityCurve {
        FiniteThermoplasticMaterialFamily::caloric_source(self)
    }
}
impl CaloricMaterial for ReferenceCaloricMaterial<'_> {
    fn caloric_source(&self) -> &ResolvedHeatCapacityCurve {
        ReferenceCaloricMaterial::caloric_source(self)
    }
}
impl CaloricMaterial for super::porous::material::FiniteThermalPorousMaterialBinding<'_> {
    fn caloric_source(&self) -> &ResolvedHeatCapacityCurve {
        super::porous::material::FiniteThermalPorousMaterialBinding::caloric_source(self)
    }
}

/// The original mechanical owner's solid and explicitly coupled storage account.
/// Nonthermal irreversible energy and independent gas storage remain explicit;
/// neither is silently deposited into the separate solid caloric chart.
pub(in crate::solid_enthalpy) trait MechanicalEnergy {
    fn internal_energy_change_j(&self) -> f64;
    fn endpoint_work_defect_j(&self) -> f64;
    fn dead_load_endpoint_work_j(&self) -> f64;
    fn support_endpoint_work_j(&self) -> f64;
    /// Path work from held follower pressure, separate from dead/support loads.
    fn follower_pressure_work_j(&self) -> f64 {
        0.0
    }
    /// Actual work of prescribed moving contact tools, separate from support work.
    fn prescribed_tool_work_j(&self) -> f64 {
        0.0
    }
    /// Direct stored-energy change of independently sealed source gases.
    fn cavity_gas_storage_change_j(&self) -> f64 {
        0.0
    }
    /// Work of the cavities' explicitly supplied exterior pressure fields.
    fn cavity_exterior_work_j(&self) -> f64 {
        0.0
    }
    fn nonthermal_loss_j(&self) -> f64;
}

impl MechanicalEnergy for FiniteThermoplasticTetEnergy {
    fn internal_energy_change_j(&self) -> f64 {
        self.internal_energy_change_j
    }
    fn endpoint_work_defect_j(&self) -> f64 {
        self.endpoint_work_defect_j
    }
    fn dead_load_endpoint_work_j(&self) -> f64 {
        self.dead_load_endpoint_work_j
    }
    fn support_endpoint_work_j(&self) -> f64 {
        self.support_endpoint_work_j
    }
    fn nonthermal_loss_j(&self) -> f64 {
        0.0
    }
}

/// Exact original cell/material pairing, independent of mechanical dimension.
pub(in crate::solid_enthalpy) struct ThermalBinding<'b, 'm, M> {
    pub(in crate::solid_enthalpy) cells: &'b [DepositedCell],
    pub(in crate::solid_enthalpy) materials: FinitePlasticLawField<'m, M>,
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

    pub(in crate::solid_enthalpy) fn enthalpies<E>(
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

    pub(in crate::solid_enthalpy) fn temperatures<E>(
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
pub(in crate::solid_enthalpy) struct MechanicalProposal<S, E, A = FiniteThermoplasticTetEnergy> {
    pub(in crate::solid_enthalpy) state: S,
    pub(in crate::solid_enthalpy) equilibrium: E,
    pub(in crate::solid_enthalpy) energy: A,
    pub(in crate::solid_enthalpy) inertia: Option<InertialEnergy>,
    pub(in crate::solid_enthalpy) heat_j: Vec<f64>,
}

/// Complete accepted thermal/mechanical/observer endpoint.
pub(in crate::solid_enthalpy) struct FeedbackCandidate<
    't,
    S,
    E,
    O,
    A = FiniteThermoplasticTetEnergy,
> {
    pub(in crate::solid_enthalpy) thermal: DepositionSimulation<'t>,
    pub(in crate::solid_enthalpy) state: S,
    pub(in crate::solid_enthalpy) observer: O,
    pub(in crate::solid_enthalpy) cells: Vec<DepositedCell>,
    pub(in crate::solid_enthalpy) report: FiniteThermoplasticReport<E, A>,
}

fn admit_heat<E>(heat_j: &[f64], cells: usize) -> Result<(), FiniteThermoplasticError<E>> {
    if heat_j.len() != cells || heat_j.iter().any(|q| !q.is_finite()) {
        return Err(FiniteThermoplasticError::Invalid(
            "exactly one finite integrated caloric heat per structural cell required",
        ));
    }
    Ok(())
}


/// Compose two physical source contributions without losing either silently.
/// The absolute rounding discrepancy remains separate from physical heat.
fn compose_source_power<E>(
    previous_w: f64,
    added_w: f64,
) -> Result<(f64, Dd), FiniteThermoplasticError<E>> {
    let next_w = finite(previous_w + added_w)?;
    if (previous_w != 0.0 && next_w == added_w)
        || (added_w != 0.0 && next_w == previous_w)
    {
        return Err(FiniteThermoplasticError::Invalid(
            "finite thermoplastic source addition loses a nonzero heat contribution",
        ));
    }
    let exact = Dd::from_f64(previous_w) + Dd::from_f64(added_w);
    let defect_w = finite((Dd::from_f64(next_w) - exact).to_f64())?.abs();
    Ok((next_w, Dd::from_f64(defect_w)))
}

fn check_source_composition<E>(
    defect_j: Dd,
    maximum_j: f64,
) -> Result<(), FiniteThermoplasticError<E>> {
    let defect_j = finite(defect_j.to_f64())?;
    if defect_j < 0.0 || defect_j > maximum_j {
        return Err(FiniteThermoplasticError::Invalid(
            "integrated finite thermoplastic source-composition defect exceeds the energy gate",
        ));
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
pub(in crate::solid_enthalpy) fn advance<
    't,
    S,
    E,
    O: DepositionObserver,
    M: CaloricMaterial,
    A: MechanicalEnergy,
>(
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
    )
        -> Result<MechanicalProposal<S, E, A>, FiniteThermoplasticError<O::Error>>,
) -> Result<FeedbackCandidate<'t, S, E, O, A>, FiniteThermoplasticError<O::Error>> {
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
pub(in crate::solid_enthalpy) fn advance_with_coupling<
    't,
    S,
    E,
    O: DepositionObserver,
    M: CaloricMaterial,
    A: MechanicalEnergy,
>(
    binding: ThermalBinding<'_, '_, M>,
    time_s: f64,
    thermal: &DepositionSimulation<'t>,
    observer: &O,
    target_s: f64,
    thermal_config: StepConfig,
    config: FiniteThermoplasticConfig,
    method: CouplingMethod,
    cx: &Cx<'_>,
    mechanics: impl FnMut(
        &[f64],
        f64,
    )
        -> Result<MechanicalProposal<S, E, A>, FiniteThermoplasticError<O::Error>>,
) -> Result<FeedbackCandidate<'t, S, E, O, A>, FiniteThermoplasticError<O::Error>> {
    advance_with_transport(
        binding,
        time_s,
        thermal,
        observer,
        target_s,
        thermal_config,
        config,
        method,
        None,
        cx,
        mechanics,
        |original, _| Ok(original.clone()),
        |_, _| Ok(0.0),
    )
    .map(|(candidate, _)| candidate)
}

/// Keep exact caloric coordinates for native source families whose histories
/// are expressed in specific enthalpy. The same feedback and transport loop is
/// reused; no inverse temperature-to-enthalpy chart is constructed.
#[allow(clippy::too_many_arguments)]
pub(in crate::solid_enthalpy) fn advance_with_enthalpy<
    't,
    S,
    E,
    O: DepositionObserver,
    M: CaloricMaterial,
    A: MechanicalEnergy,
>(
    binding: ThermalBinding<'_, '_, M>,
    time_s: f64,
    thermal: &DepositionSimulation<'t>,
    observer: &O,
    target_s: f64,
    thermal_config: StepConfig,
    config: FiniteThermoplasticConfig,
    method: CouplingMethod,
    cx: &Cx<'_>,
    mechanics: impl FnMut(
        &[f64],
        f64,
    )
        -> Result<MechanicalProposal<S, E, A>, FiniteThermoplasticError<O::Error>>,
) -> Result<FeedbackCandidate<'t, S, E, O, A>, FiniteThermoplasticError<O::Error>> {
    advance_with_input(
        binding,
        time_s,
        thermal,
        observer,
        target_s,
        thermal_config,
        config,
        method,
        MechanicalInput::Enthalpy,
        None,
        cx,
        mechanics,
        |original, _| Ok(original.clone()),
        |_, _| Ok(0.0),
    )
    .map(|(candidate, _)| candidate)
}

/// The same feedback owner with a transport observer prepared from each actual
/// native mechanical proposal. A separate geometric gate checks the transport
/// used for heat against the mechanical state at the returned caloric endpoint.
/// Every replay still clones the SAME original observer and accepted timeline.
#[allow(clippy::too_many_arguments)]
pub(in crate::solid_enthalpy) fn advance_with_transport<
    't,
    S,
    E,
    O: DepositionObserver,
    M: CaloricMaterial,
    A: MechanicalEnergy,
>(
    binding: ThermalBinding<'_, '_, M>,
    time_s: f64,
    thermal: &DepositionSimulation<'t>,
    observer: &O,
    target_s: f64,
    thermal_config: StepConfig,
    config: FiniteThermoplasticConfig,
    method: CouplingMethod,
    maximum_geometry_relative_defect: Option<f64>,
    cx: &Cx<'_>,
    mechanics: impl FnMut(
        &[f64],
        f64,
    )
        -> Result<MechanicalProposal<S, E, A>, FiniteThermoplasticError<O::Error>>,
    prepare_transport: impl FnMut(&O, &S) -> Result<O, FiniteThermoplasticError<O::Error>>,
    geometry_difference: impl FnMut(&O, &O) -> Result<f64, FiniteThermoplasticError<O::Error>>,
) -> Result<(FeedbackCandidate<'t, S, E, O, A>, f64), FiniteThermoplasticError<O::Error>> {
    advance_with_input(
        binding,
        time_s,
        thermal,
        observer,
        target_s,
        thermal_config,
        config,
        method,
        MechanicalInput::Temperature,
        maximum_geometry_relative_defect,
        cx,
        mechanics,
        prepare_transport,
        geometry_difference,
    )
}

#[derive(Clone, Copy)]
enum MechanicalInput {
    Temperature,
    Enthalpy,
}
impl MechanicalInput {
    fn values<'a>(self, enthalpies: &'a [f64], temperatures: &'a [f64]) -> &'a [f64] {
        match self {
            Self::Temperature => temperatures,
            Self::Enthalpy => enthalpies,
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn advance_with_input<
    't,
    S,
    E,
    O: DepositionObserver,
    M: CaloricMaterial,
    A: MechanicalEnergy,
>(
    binding: ThermalBinding<'_, '_, M>,
    time_s: f64,
    thermal: &DepositionSimulation<'t>,
    observer: &O,
    target_s: f64,
    thermal_config: StepConfig,
    config: FiniteThermoplasticConfig,
    method: CouplingMethod,
    input: MechanicalInput,
    maximum_geometry_relative_defect: Option<f64>,
    cx: &Cx<'_>,
    mut mechanics: impl FnMut(
        &[f64],
        f64,
    ) -> Result<MechanicalProposal<S, E, A>, FiniteThermoplasticError<O::Error>>,
    mut prepare_transport: impl FnMut(&O, &S) -> Result<O, FiniteThermoplasticError<O::Error>>,
    mut geometry_difference: impl FnMut(&O, &O) -> Result<f64, FiniteThermoplasticError<O::Error>>,
) -> Result<(FeedbackCandidate<'t, S, E, O, A>, f64), FiniteThermoplasticError<O::Error>> {
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
        || maximum_geometry_relative_defect
            .is_some_and(|limit| !limit.is_finite() || !(0.0..1.0).contains(&limit))
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
    let mut geometry_relative_defect = 0.0;
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
            mechanics(input.values(&enthalpies, &trial_temperatures), duration_s)?
        };
        admit_heat(&trial.heat_j, binding.cells.len())?;
        let applied_heat_j = finite(trial.heat_j.iter().sum())?;
        let heat = FiniteHeat {
            binding,
            inner: prepare_transport(observer, &trial.state)?,
            cx,
            watts: binding
                .cells
                .iter()
                .zip(&trial.heat_j)
                .map(|(cell, q)| Ok((cell.id, finite(*q / duration_s)?)))
                .collect::<Result<_, FiniteThermoplasticError<O::Error>>>()?,
            source_composition_defect_w: Dd::ZERO,
            source_composition_defect_j: Dd::ZERO,
        };
        let (next_thermal, heat, thermal_audit) =
            thermal.advance_observed(&heat, target_s, thermal_config, || {
                cx.checkpoint().map_err(|_| DepositionError::Cancelled)
            })?;
        // A rounded source loss can otherwise be absorbed into the derived
        // external-heat account while both aggregate first-law residuals pass.
        check_source_composition(
            heat.source_composition_defect_j,
            config.energy_tolerance_j,
        )?;
        let final_h = binding.enthalpies(&next_thermal, false, cx)?;
        let final_temperatures = binding.temperatures(&final_h)?;
        let same_input = match input {
            MechanicalInput::Temperature => final_temperatures == trial_temperatures,
            // The enthalpy owner returns exact caloric coordinates, including
            // the sign of zero; reuse only the identical evaluated endpoint.
            MechanicalInput::Enthalpy => final_h
                .iter()
                .zip(&enthalpies)
                .all(|(actual, trial)| actual.to_bits() == trial.to_bits()),
        };
        let (actual, heat_defect_j) = if same_input {
            (trial, 0.0)
        } else {
            let actual = mechanics(input.values(&final_h, &final_temperatures), duration_s)?;
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
        if maximum_geometry_relative_defect.is_some() {
            let actual_transport = prepare_transport(observer, &actual.state)?;
            geometry_relative_defect = geometry_difference(&heat.inner, &actual_transport)?;
            if !geometry_relative_defect.is_finite() || geometry_relative_defect < 0.0 {
                return Err(FiniteThermoplasticError::Invalid(
                    "finite nonnegative transport geometry discrepancy required",
                ));
            }
        }
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
        let pressure_work_j = finite(actual.energy.follower_pressure_work_j())?;
        let tool_work_j = finite(actual.energy.prescribed_tool_work_j())?;
        let gas_storage_j = finite(actual.energy.cavity_gas_storage_change_j())?;
        let cavity_exterior_work_j = finite(actual.energy.cavity_exterior_work_j())?;
        let combined_residual_j =
            if pressure_work_j == 0.0 && tool_work_j == 0.0
                && gas_storage_j == 0.0 && cavity_exterior_work_j == 0.0 {
                // Preserve the exact arithmetic of existing unpressured owners.
                finite(
                    actual.energy.internal_energy_change_j()
                        + thermal_audit.stored_energy_change_j
                        + actual.energy.endpoint_work_defect_j()
                        + inertia_account
                        - actual.energy.dead_load_endpoint_work_j()
                        - actual.energy.support_endpoint_work_j()
                        - external_heat_j
                        - thermal_audit.incoming_energy_j,
                )?
            } else {
                // Gas storage and external pressure work enter once. Internal gas
                // work cancels between solid and gas; its endpoint correction is
                // already in the complete external-work defect. Compensated sums
                // retain small heat balances against large opposing boundary work.
                let kinetic = actual.inertia.as_ref().map_or(Dd::ZERO, |inertia| {
                    Dd::from_f64(inertia.kinetic_change_j) + Dd::from_f64(inertia.numerical_loss_j)
                });
                finite(
                    (Dd::from_f64(actual.energy.internal_energy_change_j())
                        + Dd::from_f64(thermal_audit.stored_energy_change_j)
                        + Dd::from_f64(gas_storage_j)
                        + Dd::from_f64(actual.energy.endpoint_work_defect_j())
                        + kinetic
                        - Dd::from_f64(actual.energy.dead_load_endpoint_work_j())
                        - Dd::from_f64(actual.energy.support_endpoint_work_j())
                        - Dd::from_f64(pressure_work_j)
                        - Dd::from_f64(tool_work_j)
                        - Dd::from_f64(cavity_exterior_work_j)
                        - Dd::from_f64(external_heat_j)
                        - Dd::from_f64(thermal_audit.incoming_energy_j))
                    .to_f64(),
                )?
            };
        let nonthermal_loss_j = finite(actual.energy.nonthermal_loss_j())?;
        if nonthermal_loss_j < 0.0 {
            return Err(FiniteThermoplasticError::Invalid(
                "nonthermal irreversible energy must be finite and nonnegative",
            ));
        }
        // Preserve the original all-thermal account's exact arithmetic path.
        let combined_residual_j = if nonthermal_loss_j == 0.0 {
            combined_residual_j
        } else {
            finite(combined_residual_j + nonthermal_loss_j)?
        };
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
            && maximum_geometry_relative_defect
                .is_none_or(|limit| geometry_relative_defect <= limit)
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
            return Ok((FeedbackCandidate {
                thermal: next_thermal,
                state: actual.state,
                observer: heat.inner,
                cells,
                report,
            }, geometry_relative_defect));
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
    if let Some(limit) = maximum_geometry_relative_defect
        && geometry_relative_defect > limit
    {
        return Err(FiniteThermoplasticError::TransportNotConverged {
            iterations: config.maximum_iterations,
            geometry_relative_defect,
            geometry_tolerance: limit,
        });
    }
    Err(FiniteThermoplasticError::NotConverged {
        iterations: config.maximum_iterations,
        temperature_defect_k: defects[0],
        heat_defect_j: defects[1],
        thermal_residual_j: defects[2],
        combined_residual_j: defects[3],
    })
}

struct FiniteHeat<'b, 'm, 'cx, O, M> {
    binding: ThermalBinding<'b, 'm, M>,
    inner: O,
    cx: &'b Cx<'cx>,
    watts: BTreeMap<CellId, f64>,
    source_composition_defect_w: Dd,
    source_composition_defect_j: Dd,
}
impl<O: Clone, M> Clone for FiniteHeat<'_, '_, '_, O, M> {
    fn clone(&self) -> Self {
        Self {
            binding: self.binding,
            inner: self.inner.clone(),
            cx: self.cx,
            watts: self.watts.clone(),
            source_composition_defect_w: self.source_composition_defect_w,
            source_composition_defect_j: self.source_composition_defect_j,
        }
    }
}
impl<O: DepositionObserver, M: CaloricMaterial> DepositionObserver
    for FiniteHeat<'_, '_, '_, O, M>
{
    type Error = FiniteThermoplasticError<O::Error>;
    fn begin(&mut self) -> Result<(), Self::Error> {
        self.source_composition_defect_w = Dd::ZERO;
        self.source_composition_defect_j = Dd::ZERO;
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
        let mut rates = collect_source_rates(thermal, |emit| {
            self.inner.visit_sources(thermal, emit, checkpoint)
        }).map_err(FiniteThermoplasticError::Observer)?;
        let mut defect_w = Dd::ZERO;
        for (&cell, &watts) in &self.watts {
            checkpoint()?;
            let rate = rates.entry(cell).or_insert(0.0);
            let (next, material_defect) = compose_source_power(*rate, watts)?;
            // Retain the scalar map view's original two-addition defect
            // account. Native visit_sources emits the independent terms and
            // leaves their final combination to the thermal owner instead.
            let authored = thermal.cells().get(&cell)
                .ok_or(FiniteThermoplasticError::Invalid(
                    "finite thermoplastic source cell became inactive",
                ))?.boundary.source_w;
            let (_, boundary_defect) = compose_source_power(authored, next)?;
            defect_w = defect_w + material_defect + boundary_defect;
            *rate = next;
        }
        finite::<O::Error>(defect_w.to_f64())?;
        self.source_composition_defect_w = defect_w;
        Ok(rates)
    }
    fn visit_sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        emit: &mut dyn FnMut(CellId, f64) -> Result<(), DepositionError>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.binding.enthalpies(thermal, false, self.cx)?;
        // The shared thermal owner retains the original summands and accounts
        // for their final rounding. Do not charge a second, hypothetical merge.
        self.source_composition_defect_w = Dd::ZERO;
        self.inner
            .visit_sources(thermal, emit, checkpoint)
            .map_err(FiniteThermoplasticError::Observer)?;
        for (&cell, &watts) in &self.watts {
            checkpoint()?;
            emit(cell, watts)?;
        }
        Ok(())
    }
    fn visit_step_sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        start: f64,
        dt: f64,
        end: f64,
        emit: &mut dyn FnMut(CellId, f64) -> Result<(), DepositionError>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.inner
            .visit_step_sources(thermal, start, dt, end, emit, checkpoint)
            .map_err(FiniteThermoplasticError::Observer)
    }
    fn step_sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        start: f64,
        dt: f64,
        end: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        collect_source_rates(thermal, |emit| {
            self.visit_step_sources(thermal, start, dt, end, emit, checkpoint)
        })
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
        self.source_composition_defect_j = self.source_composition_defect_j
            + self.source_composition_defect_w * Dd::from_f64(duration);
        finite::<O::Error>(self.source_composition_defect_j.to_f64())?;
        self.binding.enthalpies(thermal, false, self.cx).map(|_| ())
    }
}

#[cfg(test)]
mod source_composition_tests {
    use super::*;
    use std::convert::Infallible;

    #[test]
    fn weak_interface_heat_cannot_disappear_beside_material_heat() {
        let material_w = 2.0_f64.powi(54);
        // A one-watt physical film is below half an ulp of this material rate.
        // An aggregate audit can misclassify that lost heat as an external sink.
        assert_eq!(material_w + 1.0, material_w);
        for (before, added) in [(1.0, material_w), (material_w, 1.0), (-1.0, material_w)] {
            assert!(compose_source_power::<Infallible>(before, added).is_err());
        }
        let (power, error) = compose_source_power::<Infallible>(7.0, -7.0).unwrap();
        assert_eq!(power, 0.0);
        assert_eq!(error.to_f64(), 0.0);
    }

    #[test]
    fn opposite_rounding_errors_spend_the_absolute_actual_interval_budget() {
        let large = 2.0_f64.powi(53);
        // Round-to-even makes +3 round upward by 1 W and +5 downward by 1 W.
        let (a, error_a) = compose_source_power::<Infallible>(large, 3.0).unwrap();
        let (b, error_b) = compose_source_power::<Infallible>(large, 5.0).unwrap();
        assert_eq!(a, large + 4.0);
        assert_eq!(b, large + 4.0);
        assert_eq!(error_a.to_f64(), 1.0);
        assert_eq!(error_b.to_f64(), 1.0);
        // These are the native accepted substep durations, not the proposed
        // maximum dt. Opposite cell/source errors must not cancel.
        let integrated = error_a * Dd::from_f64(0.125)
            + error_b * Dd::from_f64(0.375);
        assert_eq!(integrated.to_f64(), 0.5);
        assert!(check_source_composition::<Infallible>(integrated, 0.25).is_err());
        assert!(check_source_composition::<Infallible>(integrated, 0.5).is_ok());
    }
}

