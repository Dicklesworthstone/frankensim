//! Native tensor heat and DC current on the same original caloric endpoint.
//!
//! Every coupling trial restarts the same physical storage. Only the accepted
//! electrical source deposits heat; endpoint current is separate evidence.

use super::{invalid, validate_joule, validate_temperature};
use crate::conjugate::radiation::{RadiationConfig, Radiators};
use crate::conjugate::tensor_conduction::advance_radiating_caloric_prepared;
use crate::conjugate::transient::TransientRecord;
use crate::conjugate::{
    ChtError, ElectricSetup, ElectricSolution, ElectricalContact, SolidMaterial,
    TensorCaloricEndpoint, TensorCaloricRadiationReport, TensorCaloricStep, TensorCaloricStorage,
    TensorConductivityLaw, TensorNonlinearConfig, ThermalSetup, Voxel, VoxelDomain,
    advance_caloric_tensor_conduction, poll,
};
use fs_conduction::enthalpy::DepositionError;
use fs_exec::CancelGate;
use fs_math::dd::Dd;
use std::borrow::Cow;

mod adaptive;
mod sdirk;
pub use adaptive::{
    CaloricElectrothermalAdaptiveConfig, CaloricElectrothermalAdaptiveReport,
    CaloricElectrothermalAdaptiveSolution, CaloricElectrothermalForcing,
    advance_caloric_tensor_electrothermal_adaptive,
};
pub use sdirk::{
    CaloricElectrothermalSdirk2Step, advance_caloric_tensor_electrothermal_radiating_sdirk2,
    advance_caloric_tensor_electrothermal_sdirk2,
};

/// Native per-solve limits and bounded original-energy/current feedback.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CaloricElectrothermalConfig {
    pub thermal: TensorNonlinearConfig,
    /// At most this many thermal proposals. Current solves are stages times
    /// (proposals + 1); each refreshed stage field seeds the next proposal.
    /// Radiation retains its original bounded inner heat solves per proposal.
    /// Each current callback is evaluated at most once per non-electrode cell.
    pub maximum_passes: usize,
    /// Absolute tolerance on outer specific-energy changes, J/kg; positive.
    pub specific_energy_tolerance_j_kg: f64,
    /// Relative energy scaling uses increments from the physical opening state.
    /// The same relative allowance scales each local Joule-power comparison.
    pub relative_tolerance: f64,
    pub power_tolerance_w: f64,
    /// Sum of source mismatch, source-transfer rounding, applied/endpoint
    /// electrical defects and absolute native thermal closure, in joules.
    pub energy_tolerance_j: f64,
}

impl Default for CaloricElectrothermalConfig {
    fn default() -> Self {
        Self {
            thermal: TensorNonlinearConfig::default(),
            maximum_passes: 64,
            specific_energy_tolerance_j_kg: 1e-6,
            relative_tolerance: 1e-8,
            power_tolerance_w: 1e-9,
            energy_tolerance_j: 1e-6,
        }
    }
}

/// Physical source readback and numerical acceptance defects for one interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CaloricElectrothermalReport {
    pub passes: usize,
    pub current_solves: usize,
    pub maximum_energy_change_j_kg: f64,
    pub maximum_normalized_energy_change: f64,
    pub maximum_joule_source_change_w: f64,
    pub applied_electrode_energy_j: f64,
    /// All applied bulk/contact heat, including the thermostat channel below.
    pub applied_joule_energy_j: f64,
    /// Applied Joule heat on prescribed-energy cells, removed by their
    /// thermostats. The native record omits this equal source/removal pair.
    /// It is already part of applied_joule_energy_j, not additional heat.
    pub thermostat_joule_energy_j: f64,
    /// Positive tableau-weighted integral of absolute local applied/refreshed
    /// Joule-power differences; neither cells nor stages cancel each other.
    pub source_mismatch_j: f64,
    /// Combining authored/Joule powers and reducing the native source ledger.
    pub source_transfer_defect_j: f64,
    /// Electrode/total-Joule and total-Joule/deposited-field discrepancies.
    pub applied_electrical_defect_j: f64,
    /// The same interval-scaled diagnostic for the unapplied endpoint current.
    pub endpoint_electrical_defect_j: f64,
    pub thermal_closure_j: f64,
    pub absolute_energy_defect_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct CaloricElectrothermalStep {
    /// Original native endpoint and ledger; Joule heat is already included.
    pub thermal: TensorCaloricStep,
    /// The electrical solve whose cell powers actually supplied the heat solve.
    pub applied_electric: ElectricSolution,
    /// Current re-evaluated at the returned original energy and temperature.
    /// It supplied no additional heat; its difference is explicitly bounded.
    pub endpoint_electric: ElectricSolution,
    pub report: CaloricElectrothermalReport,
}

/// Advance native tensor conduction and original-energy-dependent Joule heat.
///
/// The resistivity callback receives (cell, occupancy, original specific energy,
/// temperature). None means explicit insulation, following the electrical
/// owner's semantics; unsupported material data must return an error. Thus a
/// latent fraction can change current and heating at exactly constant Kelvin.
/// Original caloric and tensor laws, electrode drives and contacts remain fixed.
///
/// Each outer proposal calls the SAME native backward-Euler owner with the
/// SAME opening masses and energies. It changes only the electrical source
/// estimate, never advances a discarded physical state or applies heat twice.
/// Acceptance checks original-energy increments and local source changes, plus
/// an absolute joule gate that cannot hide opposing cell defects. Source work,
/// native nonlinear/Krylov work and cancellation remain with their owners;
/// maximum_passes bounds their product. No new electrical/heat stencil is used.
///
/// The native heat ledger excludes authored powers on prescribed-energy cells.
/// Electrical heat deposited there is retained as thermostat_joule_energy_j:
/// an equal source and thermostat removal outside that ledger. Free-cell Joule
/// heat is already in thermal.record.source_j. Do not add either heat twice.
///
/// This is fixed-mass, fixed-geometry DC coupling and first-order time stepping.
/// The bounded source iteration need not converge for arbitrary feedback laws.
/// It adds no expansion work, chemistry, advection or electrical energy storage.
/// Any refusal returns no partial interval and leaves all supplied inventories
/// unchanged; mutable diagnostic state captured by callbacks remains theirs.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn advance_caloric_tensor_electrothermal(
    domain: &VoxelDomain,
    solids: &[SolidMaterial],
    electric: &ElectricSetup,
    contacts: &[ElectricalContact],
    thermal: &ThermalSetup,
    conductivity_laws: &[&dyn TensorConductivityLaw],
    storage: TensorCaloricStorage<'_>,
    prescribed_specific_energy_j_kg: &[(usize, f64)],
    start_time_s: f64,
    time_step_s: f64,
    config: &CaloricElectrothermalConfig,
    mut resistivity: impl FnMut(usize, Voxel, f64, f64) -> Result<Option<f64>, ChtError>,
    gate: &CancelGate,
) -> Result<CaloricElectrothermalStep, ChtError> {
    advance_with_heat(
        domain,
        solids,
        [electric],
        [contacts],
        [thermal],
        conductivity_laws,
        storage,
        start_time_s,
        time_step_s,
        [start_time_s + time_step_s],
        [1.0],
        config,
        |_, cell, occupancy, energy, temperature| resistivity(cell, occupancy, energy, temperature),
        |heated| {
            advance_caloric_tensor_conduction(
                domain,
                solids,
                heated[0],
                conductivity_laws,
                storage,
                prescribed_specific_energy_j_kg,
                start_time_s,
                time_step_s,
                &config.thermal,
                gate,
            )
        },
        gate,
    )
    .map(CoupledStep::backward_euler)
}

/// Advance spatial caloric/current feedback with original gray radiation.
///
/// Each electrical-source proposal converges the original radiating thermal
/// owner from the SAME physical opening inventory. Surface exchange and heat
/// escaping to the declared surroundings enter that owner's physical ledger;
/// neither is a second Joule deposit. Fixed-cell electrical heat retains the
/// explicit equal thermostat source/removal convention of the nonradiating API.
///
/// Immutable ray geometry is prepared once for the entire coupled interval.
/// Every proposal retains the native radiation, Newton, Krylov and backtracking
/// limits, while maximum_passes bounds the number of proposals. Original-law
/// query budgets stay with the supplied laws across all proposals. The returned
/// radiation report belongs to the accepted physical endpoint only.
///
/// This retains cell-centered gray/diffuse/opaque surface radiation, transparent
/// gas and fixed geometry/emissivities/surroundings. Authored sinks must not
/// already contain the same radiation. Time stepping remains backward Euler;
/// unsupported source data, nonconvergence and cancellation return no interval.
#[allow(clippy::too_many_arguments)]
pub fn advance_caloric_tensor_electrothermal_radiating(
    domain: &VoxelDomain,
    solids: &[SolidMaterial],
    electric: &ElectricSetup,
    contacts: &[ElectricalContact],
    thermal: &ThermalSetup,
    conductivity_laws: &[&dyn TensorConductivityLaw],
    storage: TensorCaloricStorage<'_>,
    prescribed_specific_energy_j_kg: &[(usize, f64)],
    start_time_s: f64,
    time_step_s: f64,
    config: &CaloricElectrothermalConfig,
    radiation: &RadiationConfig,
    mut resistivity: impl FnMut(usize, Voxel, f64, f64) -> Result<Option<f64>, ChtError>,
    gate: &CancelGate,
) -> Result<(CaloricElectrothermalStep, TensorCaloricRadiationReport), ChtError> {
    poll(gate)?;
    admit(config, start_time_s, time_step_s, domain.cell_count())?;
    domain.check_materials(solids.len())?;
    config
        .maximum_passes
        .checked_mul(radiation.max_iterations)
        .ok_or_else(|| refusal("radiating thermal work bound is unrepresentable"))?;
    let radiators = Radiators::build(domain, solids, radiation, gate)?;
    let mut accepted_radiation = None;
    let step = advance_with_heat(
        domain,
        solids,
        [electric],
        [contacts],
        [thermal],
        conductivity_laws,
        storage,
        start_time_s,
        time_step_s,
        [start_time_s + time_step_s],
        [1.0],
        config,
        |_, cell, occupancy, energy, temperature| resistivity(cell, occupancy, energy, temperature),
        |heated| {
            let (step, report) = advance_radiating_caloric_prepared(
                domain,
                solids,
                heated[0],
                conductivity_laws,
                storage,
                prescribed_specific_energy_j_kg,
                start_time_s,
                time_step_s,
                &config.thermal,
                radiation,
                Some(&radiators),
                gate,
            )?;
            accepted_radiation = Some(report);
            Ok(step)
        },
        gate,
    )?;
    let report = accepted_radiation
        .ok_or_else(|| refusal("accepted radiating thermal proposal is missing its report"))?;
    Ok((step.backward_euler(), report))
}

/// Borrow native evidence and move stage energies without copying whole solves.
trait HeatStep<const STAGES: usize> {
    fn endpoints(&self) -> [&TensorCaloricEndpoint; STAGES];
    fn record(&self) -> &TransientRecord;
    fn into_energies(self) -> [Vec<f64>; STAGES];
}

impl HeatStep<1> for TensorCaloricStep {
    fn endpoints(&self) -> [&TensorCaloricEndpoint; 1] {
        [&self.endpoint]
    }

    fn record(&self) -> &TransientRecord {
        &self.record
    }

    fn into_energies(self) -> [Vec<f64>; 1] {
        [self.endpoint.specific_energy_j_kg]
    }
}

struct CoupledStep<Step, const STAGES: usize> {
    thermal: Step,
    applied_electric: [ElectricSolution; STAGES],
    endpoint_electric: [ElectricSolution; STAGES],
    report: CaloricElectrothermalReport,
}

impl CoupledStep<TensorCaloricStep, 1> {
    fn backward_euler(self) -> CaloricElectrothermalStep {
        let [applied_electric] = self.applied_electric;
        let [endpoint_electric] = self.endpoint_electric;
        CaloricElectrothermalStep {
            thermal: self.thermal,
            applied_electric,
            endpoint_electric,
            report: self.report,
        }
    }
}

/// One stage-aware electrical-source iteration around the original heat owner.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn advance_with_heat<const STAGES: usize, Step: HeatStep<STAGES>>(
    domain: &VoxelDomain,
    solids: &[SolidMaterial],
    electric: [&ElectricSetup; STAGES],
    contacts: [&[ElectricalContact]; STAGES],
    thermal: [&ThermalSetup; STAGES],
    conductivity_laws: &[&dyn TensorConductivityLaw],
    storage: TensorCaloricStorage<'_>,
    start_time_s: f64,
    time_step_s: f64,
    stage_time_s: [f64; STAGES],
    stage_weights: [f64; STAGES],
    config: &CaloricElectrothermalConfig,
    mut resistivity: impl FnMut(f64, usize, Voxel, f64, f64) -> Result<Option<f64>, ChtError>,
    mut heat: impl FnMut([&ThermalSetup; STAGES]) -> Result<Step, ChtError>,
    gate: &CancelGate,
) -> Result<CoupledStep<Step, STAGES>, ChtError> {
    poll(gate)?;
    let cells = domain.cell_count();
    admit(config, start_time_s, time_step_s, cells)?;
    config
        .maximum_passes
        .checked_add(1)
        .and_then(|count| count.checked_mul(STAGES))
        .and_then(|count| count.checked_mul(cells))
        .ok_or_else(|| refusal("stage current callback work bound is unrepresentable"))?;
    domain.check_materials(solids.len())?;
    if storage.laws.len() != cells
        || storage.reference_mass_kg.len() != cells
        || storage.previous_specific_energy_j_kg.len() != cells
        || conductivity_laws.len() != cells
        || thermal
            .iter()
            .any(|setup| !setup.power_w.is_empty() && setup.power_w.len() != cells)
    {
        return Err(refusal(
            "one original law, mass, energy and tensor per cell is required",
        ));
    }
    let mut temperature = Vec::new();
    temperature
        .try_reserve_exact(cells)
        .map_err(|_| refusal("opening caloric temperature allocation refused"))?;
    let mut fixed: [Vec<bool>; STAGES] = std::array::from_fn(|_| Vec::new());
    for stage in 0..STAGES {
        fixed[stage]
            .try_reserve_exact(cells)
            .map_err(|_| refusal("prescribed-cell allocation refused"))?;
        fixed[stage].resize(cells, false);
        for (index, &(cell, _)) in thermal[stage].fixed_temperature.iter().enumerate() {
            checkpoint(gate, index)?;
            let entry = fixed[stage]
                .get_mut(cell)
                .ok_or_else(|| refusal("prescribed cell is outside the domain"))?;
            if *entry {
                return Err(refusal("prescribed cells must be unique"));
            }
            *entry = true;
        }
    }
    for cell in 0..cells {
        checkpoint(gate, cell)?;
        let mass = storage.reference_mass_kg[cell];
        let energy = storage.previous_specific_energy_j_kg[cell];
        let law = storage.laws[cell];
        let bounds = law.specific_energy_bounds_j_kg();
        if !mass.is_finite()
            || mass <= 0.0
            || !energy.is_finite()
            || bounds.iter().any(|e| !e.is_finite())
            || bounds[0] >= bounds[1]
            || !(bounds[0]..=bounds[1]).contains(&energy)
        {
            return Err(refusal(
                "opening inventory leaves its original caloric bounds",
            ));
        }
        law.validate_reference_cell(mass, domain.volume(cell))
            .map_err(caloric_error)?;
        temperature.push(law.temperature_k(energy).map_err(caloric_error)?);
        if thermal.iter().any(|setup| {
            setup
                .power_w
                .get(cell)
                .is_some_and(|power| !power.is_finite())
        }) {
            return Err(refusal("authored cell powers must be finite"));
        }
    }
    validate_temperature(&temperature, cells)?;
    let mut current = |stage: usize, energy: &[f64], temperature: &[f64]| {
        poll(gate)?;
        let field = electric[stage].solve_with_cell_contacts(
            domain,
            temperature,
            contacts[stage],
            |cell, occupancy, t| {
                poll(gate)?;
                let value = resistivity(stage_time_s[stage], cell, occupancy, energy[cell], t)?;
                poll(gate)?;
                Ok(value)
            },
            gate,
        )?;
        poll(gate)?;
        Ok::<_, ChtError>(field)
    };
    let mut trial_energy: [Cow<'_, [f64]>; STAGES] =
        std::array::from_fn(|_| Cow::Borrowed(storage.previous_specific_energy_j_kg));
    let mut applied = stage_array(|stage| current(stage, &trial_energy[stage], &temperature))?;
    let mut heated = thermal.map(Clone::clone);
    for setup in &mut heated {
        setup.power_w.resize(cells, 0.0);
    }
    let dt = Dd::from_f64(time_step_s);
    let weighted_time = stage_weights.map(|weight| dt * Dd::from_f64(weight));
    let mut residual = f64::INFINITY;
    for pass in 1..=config.maximum_passes {
        poll(gate)?;
        let mut applied_joule = Dd::ZERO;
        let mut applied_electrode = Dd::ZERO;
        let mut applied_electrical_defect = Dd::ZERO;
        let mut transfer = Dd::ZERO;
        let mut native_source = Dd::ZERO;
        let mut thermostat_joule = Dd::ZERO;
        for stage in 0..STAGES {
            let weight = weighted_time[stage];
            let (joule, defect) = electrical_account(&applied[stage], cells, weight, gate)?;
            applied_joule = applied_joule + joule;
            applied_electrical_defect = applied_electrical_defect + defect;
            applied_electrode =
                applied_electrode + weight * Dd::from_f64(applied[stage].delivered_w);
            let mut stage_transfer = Dd::ZERO;
            let mut stage_source = Dd::ZERO;
            let mut stage_thermostat = Dd::ZERO;
            for cell in 0..cells {
                checkpoint(gate, cell)?;
                let authored = thermal[stage].power_w.get(cell).copied().unwrap_or(0.0);
                let original = Dd::from_f64(authored) + Dd::from_f64(applied[stage].joule_w[cell]);
                heated[stage].power_w[cell] = number(original)?;
                if fixed[stage][cell] {
                    stage_thermostat =
                        stage_thermostat + Dd::from_f64(applied[stage].joule_w[cell]);
                } else {
                    let represented = Dd::from_f64(heated[stage].power_w[cell]);
                    stage_transfer = stage_transfer + magnitude(represented - original);
                    stage_source = stage_source + represented;
                }
            }
            transfer = transfer + weight * stage_transfer;
            native_source = native_source + weight * stage_source;
            thermostat_joule = thermostat_joule + weight * stage_thermostat;
        }
        let step = heat(heated.each_ref())?;
        poll(gate)?;
        let endpoints = step.endpoints();
        let next = stage_array(|stage| {
            current(
                stage,
                &endpoints[stage].specific_energy_j_kg,
                &endpoints[stage].solution.solution.temperature,
            )
        })?;
        let mut endpoint_electrical_defect = Dd::ZERO;
        let mut maximum_energy_change: f64 = 0.0;
        let mut normalized_energy: f64 = 0.0;
        let mut maximum_power_change: f64 = 0.0;
        let mut normalized_power: f64 = 0.0;
        let mut source_mismatch = Dd::ZERO;
        for stage in 0..STAGES {
            let (_, defect) = electrical_account(&next[stage], cells, weighted_time[stage], gate)?;
            endpoint_electrical_defect = endpoint_electrical_defect + defect;
            let mut stage_mismatch = Dd::ZERO;
            for cell in 0..cells {
                checkpoint(gate, cell)?;
                let opening = Dd::from_f64(storage.previous_specific_energy_j_kg[cell]);
                let previous = Dd::from_f64(trial_energy[stage][cell]);
                let endpoint = Dd::from_f64(endpoints[stage].specific_energy_j_kg[cell]);
                let change = number(magnitude(endpoint - previous))?;
                let scale = number(magnitude(endpoint - opening))?
                    .max(number(magnitude(previous - opening))?);
                let limit =
                    config.specific_energy_tolerance_j_kg + config.relative_tolerance * scale;
                maximum_energy_change = maximum_energy_change.max(change);
                normalized_energy = normalized_energy.max(ratio(change, limit)?);
                let power_change = magnitude(
                    Dd::from_f64(applied[stage].joule_w[cell])
                        - Dd::from_f64(next[stage].joule_w[cell]),
                );
                let power_change_w = number(power_change)?;
                let power_limit = config.power_tolerance_w
                    + config.relative_tolerance
                        * applied[stage].joule_w[cell].max(next[stage].joule_w[cell]);
                maximum_power_change = maximum_power_change.max(power_change_w);
                normalized_power = normalized_power.max(ratio(power_change_w, power_limit)?);
                stage_mismatch = stage_mismatch + power_change;
            }
            source_mismatch = source_mismatch + weighted_time[stage] * stage_mismatch;
        }
        let record = step.record();
        transfer = transfer + magnitude(Dd::from_f64(record.source_j) - native_source);
        let base_defect = transfer
            + applied_electrical_defect
            + endpoint_electrical_defect
            + magnitude(Dd::from_f64(record.closure_j));
        let total_defect = base_defect + source_mismatch;
        // These are original solver/representation defects, not source guesses.
        // More coupling passes cannot authorize an already excessive heat error.
        if exceeds(base_defect, config.energy_tolerance_j)? {
            return Err(refusal(
                "electrical, transfer or native thermal joule allowance exceeded",
            ));
        }
        residual = normalized_energy
            .max(normalized_power)
            .max(ratio(number(total_defect)?, config.energy_tolerance_j)?);
        if normalized_energy <= 1.0
            && normalized_power <= 1.0
            && !exceeds(total_defect, config.energy_tolerance_j)?
        {
            let report = CaloricElectrothermalReport {
                passes: pass,
                current_solves: STAGES * (pass + 1),
                maximum_energy_change_j_kg: maximum_energy_change,
                maximum_normalized_energy_change: normalized_energy,
                maximum_joule_source_change_w: maximum_power_change,
                applied_electrode_energy_j: number(applied_electrode)?,
                applied_joule_energy_j: number(applied_joule)?,
                thermostat_joule_energy_j: number(thermostat_joule)?,
                source_mismatch_j: number(source_mismatch)?,
                source_transfer_defect_j: number(transfer)?,
                applied_electrical_defect_j: number(applied_electrical_defect)?,
                endpoint_electrical_defect_j: number(endpoint_electrical_defect)?,
                thermal_closure_j: record.closure_j,
                absolute_energy_defect_j: number(total_defect)?,
            };
            poll(gate)?;
            return Ok(CoupledStep {
                thermal: step,
                applied_electric: applied,
                endpoint_electric: next,
                report,
            });
        }
        trial_energy = step.into_energies().map(Cow::Owned);
        applied = next;
    }
    Err(ChtError::SolverNotConverged {
        system: "caloric electrothermal coupling",
        iterations: config.maximum_passes,
        relative_residual: residual,
        tolerance: 1.0,
    })
}

fn stage_array<T, const STAGES: usize>(
    mut value: impl FnMut(usize) -> Result<T, ChtError>,
) -> Result<[T; STAGES], ChtError> {
    let mut stages = Vec::new();
    stages
        .try_reserve_exact(STAGES)
        .map_err(|_| refusal("stage field allocation refused"))?;
    for stage in 0..STAGES {
        stages.push(value(stage)?);
    }
    stages
        .try_into()
        .map_err(|_| refusal("native stage field shape is inconsistent"))
}

fn admit(
    config: &CaloricElectrothermalConfig,
    start: f64,
    dt: f64,
    cells: usize,
) -> Result<(), ChtError> {
    if !start.is_finite()
        || start < 0.0
        || !dt.is_finite()
        || dt <= 0.0
        || !(start + dt).is_finite()
        || start + dt <= start
        || config.maximum_passes == 0
        || !config.specific_energy_tolerance_j_kg.is_finite()
        || config.specific_energy_tolerance_j_kg <= 0.0
        || !config.relative_tolerance.is_finite()
        || !(0.0..1.0).contains(&config.relative_tolerance)
        || [config.power_tolerance_w, config.energy_tolerance_j]
            .iter()
            .any(|x| !x.is_finite() || *x < 0.0)
    {
        return Err(refusal(
            "invalid physical interval or caloric/current coupling controls",
        ));
    }
    config
        .maximum_passes
        .checked_add(1)
        .and_then(|n| n.checked_mul(cells))
        .ok_or_else(|| refusal("current callback work bound is unrepresentable"))?;
    for cap in [
        config.thermal.max_newton_iterations,
        config.thermal.max_backtracks,
        config.thermal.linear.max_iterations,
    ] {
        config
            .maximum_passes
            .checked_mul(cap)
            .ok_or_else(|| refusal("native thermal work bound is unrepresentable"))?;
    }
    Ok(())
}

fn electrical_account(
    field: &ElectricSolution,
    cells: usize,
    dt: Dd,
    gate: &CancelGate,
) -> Result<(Dd, Dd), ChtError> {
    validate_joule(field, cells)?;
    if !field.total_joule_w.is_finite()
        || field.total_joule_w < 0.0
        || !field.delivered_w.is_finite()
    {
        return Err(refusal(
            "electrical power accounting is nonfinite or nonphysical",
        ));
    }
    let mut deposited = Dd::ZERO;
    for (cell, &power) in field.joule_w.iter().enumerate() {
        checkpoint(gate, cell)?;
        deposited = deposited + Dd::from_f64(power);
    }
    let total = Dd::from_f64(field.total_joule_w);
    let defect = magnitude(Dd::from_f64(field.delivered_w) - total) + magnitude(total - deposited);
    let deposited = dt * deposited;
    let defect = dt * defect;
    number(deposited)?;
    number(defect)?;
    Ok((deposited, defect))
}

fn magnitude(value: Dd) -> Dd {
    if value.hi < 0.0 || (value.hi == 0.0 && value.lo < 0.0) {
        Dd::ZERO - value
    } else {
        value
    }
}

fn number(value: Dd) -> Result<f64, ChtError> {
    let scalar = value.hi + value.lo;
    if scalar.is_finite() {
        Ok(scalar)
    } else {
        Err(refusal("unrepresentable coupling accounting"))
    }
}

fn exceeds(value: Dd, allowance: f64) -> Result<bool, ChtError> {
    number(value)?;
    let excess = value - Dd::from_f64(allowance);
    Ok(excess.hi > 0.0 || (excess.hi == 0.0 && excess.lo > 0.0))
}

fn ratio(change: f64, allowance: f64) -> Result<f64, ChtError> {
    if !change.is_finite() || !allowance.is_finite() {
        return Err(refusal("unrepresentable coupling tolerance comparison"));
    }
    Ok(if change == 0.0 {
        0.0
    } else {
        change / allowance
    })
}

fn checkpoint(gate: &CancelGate, cell: usize) -> Result<(), ChtError> {
    if cell.is_multiple_of(256) {
        poll(gate)?;
    }
    Ok(())
}

fn caloric_error(error: DepositionError) -> ChtError {
    match error {
        DepositionError::Cancelled => ChtError::Cancelled,
        other => invalid(
            "material_coupling.caloric.original_source",
            &other.to_string(),
        ),
    }
}

fn refusal(reason: &str) -> ChtError {
    invalid("material_coupling.caloric", reason)
}
