//! Native tensor heat and DC current on the same original caloric endpoint.
//!
//! Every coupling trial restarts the same physical storage. Only the accepted
//! electrical source deposits heat; endpoint current is separate evidence.

use super::{invalid, validate_joule, validate_temperature};
use crate::conjugate::{
    ChtError, ElectricSetup, ElectricSolution, ElectricalContact, SolidMaterial, TensorCaloricStep,
    TensorCaloricStorage, TensorConductivityLaw, TensorNonlinearConfig, ThermalSetup, Voxel,
    VoxelDomain, advance_caloric_tensor_conduction, poll,
};
use fs_conduction::enthalpy::DepositionError;
use fs_exec::CancelGate;
use fs_math::dd::Dd;
use std::borrow::Cow;

/// Native per-solve limits and bounded original-energy/current feedback.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CaloricElectrothermalConfig {
    pub thermal: TensorNonlinearConfig,
    /// At most this many native heat solves and one additional current solve.
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
    /// dt * sum |applied local Joule power - refreshed endpoint local power|.
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
    poll(gate)?;
    let cells = domain.cell_count();
    admit(config, start_time_s, time_step_s, cells)?;
    domain.check_materials(solids.len())?;
    if storage.laws.len() != cells
        || storage.reference_mass_kg.len() != cells
        || storage.previous_specific_energy_j_kg.len() != cells
        || conductivity_laws.len() != cells
        || (!thermal.power_w.is_empty() && thermal.power_w.len() != cells)
    {
        return Err(refusal(
            "one original law, mass, energy and tensor per cell is required",
        ));
    }
    let mut temperature = Vec::new();
    temperature
        .try_reserve_exact(cells)
        .map_err(|_| refusal("opening caloric temperature allocation refused"))?;
    let mut fixed = Vec::new();
    fixed
        .try_reserve_exact(cells)
        .map_err(|_| refusal("prescribed-cell allocation refused"))?;
    fixed.resize(cells, false);
    for (index, &(cell, _)) in thermal.fixed_temperature.iter().enumerate() {
        checkpoint(gate, index)?;
        let entry = fixed
            .get_mut(cell)
            .ok_or_else(|| refusal("prescribed cell is outside the domain"))?;
        if *entry {
            return Err(refusal("prescribed cells must be unique"));
        }
        *entry = true;
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
        if thermal
            .power_w
            .get(cell)
            .is_some_and(|power| !power.is_finite())
        {
            return Err(refusal("authored cell powers must be finite"));
        }
    }
    validate_temperature(&temperature, cells)?;
    let mut current = |energy: &[f64], temperature: &[f64]| {
        poll(gate)?;
        let field = electric.solve_with_cell_contacts(
            domain,
            temperature,
            contacts,
            |cell, occupancy, t| {
                poll(gate)?;
                let value = resistivity(cell, occupancy, energy[cell], t)?;
                poll(gate)?;
                Ok(value)
            },
            gate,
        )?;
        poll(gate)?;
        Ok::<_, ChtError>(field)
    };
    let mut trial_energy = Cow::Borrowed(storage.previous_specific_energy_j_kg);
    let mut applied = current(&trial_energy, &temperature)?;
    let mut heated = thermal.clone();
    heated.power_w.resize(cells, 0.0);
    let dt = Dd::from_f64(time_step_s);
    let mut residual = f64::INFINITY;
    for pass in 1..=config.maximum_passes {
        poll(gate)?;
        let (applied_joule, applied_electrical_defect) =
            electrical_account(&applied, cells, dt, gate)?;
        let mut transfer = Dd::ZERO;
        let mut native_source = Dd::ZERO;
        let mut thermostat_joule = Dd::ZERO;
        for cell in 0..cells {
            checkpoint(gate, cell)?;
            let authored = thermal.power_w.get(cell).copied().unwrap_or(0.0);
            let original = Dd::from_f64(authored) + Dd::from_f64(applied.joule_w[cell]);
            heated.power_w[cell] = number(original)?;
            if fixed[cell] {
                thermostat_joule = thermostat_joule + Dd::from_f64(applied.joule_w[cell]);
            } else {
                let represented = Dd::from_f64(heated.power_w[cell]);
                transfer = transfer + magnitude(represented - original);
                native_source = native_source + represented;
            }
        }
        let step = advance_caloric_tensor_conduction(
            domain,
            solids,
            &heated,
            conductivity_laws,
            storage,
            prescribed_specific_energy_j_kg,
            start_time_s,
            time_step_s,
            &config.thermal,
            gate,
        )?;
        poll(gate)?;
        let endpoint_energy = &step.endpoint.specific_energy_j_kg;
        let next = current(
            endpoint_energy,
            &step.endpoint.solution.solution.temperature,
        )?;
        let (_, endpoint_electrical_defect) = electrical_account(&next, cells, dt, gate)?;
        let mut maximum_energy_change: f64 = 0.0;
        let mut normalized_energy: f64 = 0.0;
        let mut maximum_power_change: f64 = 0.0;
        let mut normalized_power: f64 = 0.0;
        let mut source_mismatch = Dd::ZERO;
        for cell in 0..cells {
            checkpoint(gate, cell)?;
            let opening = Dd::from_f64(storage.previous_specific_energy_j_kg[cell]);
            let previous = Dd::from_f64(trial_energy[cell]);
            let endpoint = Dd::from_f64(endpoint_energy[cell]);
            let change = number(magnitude(endpoint - previous))?;
            let scale =
                number(magnitude(endpoint - opening))?.max(number(magnitude(previous - opening))?);
            let limit = config.specific_energy_tolerance_j_kg + config.relative_tolerance * scale;
            maximum_energy_change = maximum_energy_change.max(change);
            normalized_energy = normalized_energy.max(ratio(change, limit)?);
            let power_change =
                magnitude(Dd::from_f64(applied.joule_w[cell]) - Dd::from_f64(next.joule_w[cell]));
            let power_change_w = number(power_change)?;
            let power_limit = config.power_tolerance_w
                + config.relative_tolerance * applied.joule_w[cell].max(next.joule_w[cell]);
            maximum_power_change = maximum_power_change.max(power_change_w);
            normalized_power = normalized_power.max(ratio(power_change_w, power_limit)?);
            source_mismatch = source_mismatch + power_change;
        }
        source_mismatch = dt * source_mismatch;
        transfer =
            dt * transfer + magnitude(Dd::from_f64(step.record.source_j) - dt * native_source);
        let base_defect = transfer
            + applied_electrical_defect
            + endpoint_electrical_defect
            + magnitude(Dd::from_f64(step.record.closure_j));
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
                current_solves: pass + 1,
                maximum_energy_change_j_kg: maximum_energy_change,
                maximum_normalized_energy_change: normalized_energy,
                maximum_joule_source_change_w: maximum_power_change,
                applied_electrode_energy_j: number(dt * Dd::from_f64(applied.delivered_w))?,
                applied_joule_energy_j: number(applied_joule)?,
                thermostat_joule_energy_j: number(dt * thermostat_joule)?,
                source_mismatch_j: number(source_mismatch)?,
                source_transfer_defect_j: number(transfer)?,
                applied_electrical_defect_j: number(applied_electrical_defect)?,
                endpoint_electrical_defect_j: number(endpoint_electrical_defect)?,
                thermal_closure_j: step.record.closure_j,
                absolute_energy_defect_j: number(total_defect)?,
            };
            poll(gate)?;
            return Ok(CaloricElectrothermalStep {
                thermal: step,
                applied_electric: applied,
                endpoint_electric: next,
                report,
            });
        }
        trial_energy = Cow::Owned(step.endpoint.specific_energy_j_kg);
        applied = next;
    }
    Err(ChtError::SolverNotConverged {
        system: "caloric electrothermal coupling",
        iterations: config.maximum_passes,
        relative_residual: residual,
        tolerance: 1.0,
    })
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
