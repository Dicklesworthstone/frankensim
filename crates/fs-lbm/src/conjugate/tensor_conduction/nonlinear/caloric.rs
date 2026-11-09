//! Original caloric storage, including sharp latent plateaus, on tensor faces.

use fs_conduction::enthalpy::{DepositionError, ThermalLaw};
use fs_exec::CancelGate;

use super::super::sdirk::{interval_record_with_storage, stage_carry_with_storage, stage_clock};
use super::sdirk::remaining_controls;
use super::{
    TensorConductivityLaw, TensorNonlinearConfig, TensorNonlinearSolution, solve_coordinate,
};
use crate::conjugate::tensor_conduction::invalid;
use crate::conjugate::transient::{TransientRecord, record_conduction_step};
use crate::conjugate::{
    ChtError, SolidMaterial, ThermalSetup, VoxelDomain, finite, finite_positive, poll,
};

mod adaptive;
mod material;
mod radiation;
mod storage;
pub(in crate::conjugate) use adaptive::{
    CaloricAdaptiveForcing, CaloricAdaptiveTrial, advance_caloric_adaptive_with, embedded_ratio,
};
pub(in crate::conjugate) use radiation::{
    advance_radiating_caloric_prepared, advance_radiating_sdirk2_prepared,
};
pub use radiation::{
    TensorCaloricRadiationReport, advance_caloric_tensor_conduction_radiating,
    advance_caloric_tensor_conduction_radiating_sdirk2,
};
pub use adaptive::{
    TensorCaloricAdaptiveConfig, TensorCaloricAdaptiveError, TensorCaloricAdaptiveReport,
    TensorCaloricAdaptiveSolution, TensorCaloricForcing, TensorCaloricForcingSample,
    advance_caloric_tensor_conduction_adaptive, TensorCaloricRadiatingAdaptiveSolution,
    advance_caloric_tensor_conduction_radiating_adaptive,
};
use storage::{CaloricStage, PreparedCaloric};

/// The original caloric law plus its own bounded, analytic temperature tangent.
///
/// The canonical ThermalLaw still owns the specific-energy basis, datum,
/// temperature inverse and source domain. This trait adds no conductivity:
/// the separately supplied TensorConductivityLaw remains its sole owner.
/// Temperature derivatives are nonnegative; zero represents an actual latent
/// plateau. At a genuine corner the source selects its one-sided derivative.
/// No numerical differentiation, apparent Cp, extrapolation or phase smearing
/// is supplied by the tensor solver.
pub trait TensorCaloricLaw: ThermalLaw {
    /// A continuous admitted interval in the original specific energy, J/kg.
    fn specific_energy_bounds_j_kg(&self) -> [f64; 2];

    /// Original dT/de in kg K/J, including exact zero on latent intervals.
    fn temperature_slope_kg_k_j(&self, specific_energy_j_kg: f64) -> Result<f64, DepositionError>;

    /// A fused query may reuse the original inverse once. It must return the
    /// SAME temperature and derivative as the original separate methods.
    fn temperature_and_slope_at_specific_energy(
        &self,
        specific_energy_j_kg: f64,
    ) -> Result<[f64; 2], DepositionError> {
        Ok([
            self.temperature_k(specific_energy_j_kg)?,
            self.temperature_slope_kg_k_j(specific_energy_j_kg)?,
        ])
    }

    /// Preserve any original mass/volume normalization before solving.
    ///
    /// The default is for geometry-independent specific-energy charts.
    /// Adapters for amount-bound or volume-bound storage MUST enforce their
    /// original reference constraints here; the tensor owner does not invent
    /// deposition IDs or silently bypass a source's inventory restrictions.
    fn validate_reference_cell(
        &self,
        _reference_mass_kg: f64,
        _reference_volume_m3: f64,
    ) -> Result<(), DepositionError> {
        Ok(())
    }
}

/// Borrowed original cell inventory and specific-energy state.
#[derive(Clone, Copy)]
pub struct TensorCaloricStorage<'a> {
    /// One original caloric law per cell, immutable throughout the interval.
    pub laws: &'a [&'a dyn TensorCaloricLaw],
    /// Conserved reference mass in kg per cell, never inferred from current
    /// density or multiplied by the voxel volume a second time.
    pub reference_mass_kg: &'a [f64],
    /// Physical start state in the law's declared energy basis, J/kg.
    /// The same state appears in both SDIRK2 stage equations.
    pub previous_specific_energy_j_kg: &'a [f64],
}

/// One accepted original caloric endpoint with actual tensor face heat.
#[derive(Debug, Clone, PartialEq)]
pub struct TensorCaloricEndpoint {
    pub solution: TensorNonlinearSolution,
    /// Original h, u or declared sensible-storage coordinate from ThermalLaw.
    pub specific_energy_j_kg: Vec<f64>,
}

/// One backward-Euler interval with exact original-law stored-energy change.
#[derive(Debug, Clone, PartialEq)]
pub struct TensorCaloricStep {
    pub endpoint: TensorCaloricEndpoint,
    pub record: TransientRecord,
}

/// Both physical SDIRK2 stages and their whole-interval heat ledger.
#[derive(Debug, Clone, PartialEq)]
pub struct TensorCaloricSdirk2Step {
    pub stage_time_s: [f64; 2],
    pub stages: [TensorCaloricEndpoint; 2],
    pub record: TransientRecord,
}

/// Advance original-law tensor conduction with caloric backward-Euler storage.
///
/// Cells solve in centered specific-energy increments; face traces remain
/// Kelvin. Storage is m*(e-e0)/dt and its derivative is m/dt, including on a
/// latent plateau where dT/de is zero. Only the cell column of the original
/// Kirchhoff tensor tangent is composed with dT/de; face derivatives remain
/// unchanged. Static elimination therefore remains regular on a plateau
/// without dividing by dT/de or replacing latent heat with apparent Cp.
///
/// Each setup.fixed_temperature cell MUST have exactly one matching entry in
/// prescribed_specific_energy_j_kg. Its original T(e) must agree with the
/// authored temperature to roundoff. Temperature alone cannot select a phase
/// fraction; no caloric inverse or arbitrary plateau endpoint is inferred.
/// Prescribed boundary-face temperatures need no cell-energy prescription.
///
/// Actual stored-energy differences feed the same independent source and
/// reservoir heat ledger as sensible conduction. Fixed-cell storage is counted
/// once. The mass and geometric inventory remain fixed. The separately supplied
/// tensor law retains its original K(T) integrals or explicitly selects cellwise
/// K(e) with analytic dK/de, including changing phase inventory at constant
/// temperature. This entry does not move geometric faces or add expansion work
/// or advective energy.
///
/// Backward Euler is first order. Source-support, nonlinear, Krylov and
/// backtracking limits remain active, and no partial step is returned.
///
/// # Errors
/// Refuses invalid inventory, unsupported caloric states or tangents, missing
/// or inconsistent fixed energies, original tensor/source refusals, exhausted
/// work budgets, nonconvergence, unresolved time/storage, and cancellation.
#[allow(clippy::too_many_arguments)]
pub fn advance_caloric_tensor_conduction(
    domain: &VoxelDomain,
    solids: &[SolidMaterial],
    setup: &ThermalSetup,
    conductivity_laws: &[&dyn TensorConductivityLaw],
    storage: TensorCaloricStorage<'_>,
    prescribed_specific_energy_j_kg: &[(usize, f64)],
    start_time_s: f64,
    time_step_s: f64,
    config: &TensorNonlinearConfig,
    gate: &CancelGate,
) -> Result<TensorCaloricStep, ChtError> {
    poll(gate)?;
    let time_s = end_time(start_time_s, time_step_s)?;
    let prepared =
        PreparedCaloric::new(domain, conductivity_laws.len(), storage, time_step_s, gate)?;
    let stage = prepared.stage(setup, prescribed_specific_energy_j_kg, gate)?;
    let endpoint = solve_stage(
        domain,
        solids,
        setup,
        conductivity_laws,
        &stage,
        &vec![0.0; domain.cell_count()],
        None,
        config,
        gate,
    )?;
    poll(gate)?;
    let record = stage_record(domain, setup, &prepared, &endpoint, time_s, time_step_s)?;
    poll(gate)?;
    Ok(TensorCaloricStep { endpoint, record })
}

/// Advance original caloric storage through two original-law SDIRK2 stages.
///
/// Uses the existing L-stable two-stage tableau with immutable physical e0
/// in BOTH storage equations. Both stages reevaluate the original T(e),
/// dT/de and tensor conductivity integrals. Stage two carries the actual
/// first-stage heat internally; it is not an external source or a fabricated
/// old energy. The shared interval ledger integrates actual stage heat and
/// corrects prescribed-cell reservoir storage once, even when masks change.
///
/// Supply the two physical source/boundary/contact/sink samples at the returned
/// stage times, and one explicit fixed-energy map for each stage. All Newton,
/// Krylov and backtracking work limits are shared across the interval.
/// The initial material inventory and original law constraints remain fixed.
///
/// Second order applies within smooth caloric/forcing segments with fixed
/// constraints. A sharp phase corner or forcing discontinuity requires time
/// refinement or interval splitting; no second-order claim is made across an
/// unresolved event. Latent plateaus retain zero dT/de without smoothing.
/// L-stability does not imply arbitrary-step positivity, unconditional nonlinear
/// convergence, or a maximum principle for general rotated tensor laws.
///
/// # Errors
/// The backward-Euler caloric and spatial admissions apply to both stages.
/// Either stage's refusal or cancellation returns no partial interval.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn advance_caloric_tensor_conduction_sdirk2(
    domain: &VoxelDomain,
    solids: &[SolidMaterial],
    setups: [&ThermalSetup; 2],
    conductivity_laws: &[&dyn TensorConductivityLaw],
    storage: TensorCaloricStorage<'_>,
    prescribed_specific_energy_j_kg: [&[(usize, f64)]; 2],
    start_time_s: f64,
    time_step_s: f64,
    config: &TensorNonlinearConfig,
    gate: &CancelGate,
) -> Result<TensorCaloricSdirk2Step, ChtError> {
    poll(gate)?;
    let (a, stage_time_s) = stage_clock(start_time_s, time_step_s).ok_or_else(|| {
        invalid("caloric SDIRK2 needs finite nonnegative and strictly ordered stage times")
    })?;
    let prepared = PreparedCaloric::new(domain, conductivity_laws.len(), storage, a, gate)?;
    let first = {
        let stage = prepared.stage(setups[0], prescribed_specific_energy_j_kg[0], gate)?;
        solve_stage(
            domain,
            solids,
            setups[0],
            conductivity_laws,
            &stage,
            &vec![0.0; domain.cell_count()],
            None,
            config,
            gate,
        )?
    };
    let first_record = stage_record(domain, setups[0], &prepared, &first, stage_time_s[0], a)?;
    let carried_power_w = stage_carry_with_storage(
        domain,
        setups[0],
        &first.solution.solution,
        |cell| prepared.storage_rate_w(cell, first.specific_energy_j_kg[cell]),
        gate,
    )?;
    let second_config = remaining_controls(config, &first.solution)?;
    let initial: Vec<_> = first
        .specific_energy_j_kg
        .iter()
        .enumerate()
        .map(|(cell, &energy)| energy - storage.previous_specific_energy_j_kg[cell])
        .collect();
    let second = {
        let stage = prepared.stage(setups[1], prescribed_specific_energy_j_kg[1], gate)?;
        solve_stage(
            domain,
            solids,
            setups[1],
            conductivity_laws,
            &stage,
            &initial,
            Some(&carried_power_w),
            &second_config,
            gate,
        )?
    };
    let iterations = first
        .solution
        .solution
        .report
        .iterations
        .checked_add(second.solution.solution.report.iterations)
        .filter(|&count| count <= config.linear.max_iterations)
        .ok_or_else(|| invalid("caloric SDIRK2 stages exceeded the shared Krylov budget"))?;
    poll(gate)?;
    let record = interval_record_with_storage(
        domain,
        setups[1],
        &first_record,
        &second.solution.solution,
        |cell| prepared.stored_j(cell, second.specific_energy_j_kg[cell]),
        &carried_power_w,
        stage_time_s[1],
        a,
        iterations,
        gate,
    )?;
    poll(gate)?;
    Ok(TensorCaloricSdirk2Step {
        stage_time_s,
        stages: [first, second],
        record,
    })
}

#[allow(clippy::too_many_arguments)]
fn solve_stage(
    domain: &VoxelDomain,
    solids: &[SolidMaterial],
    setup: &ThermalSetup,
    laws: &[&dyn TensorConductivityLaw],
    stage: &CaloricStage<'_, '_>,
    initial: &[f64],
    carried_power_w: Option<&[f64]>,
    config: &TensorNonlinearConfig,
    gate: &CancelGate,
) -> Result<TensorCaloricEndpoint, ChtError> {
    let (solution, coordinates) = solve_coordinate(
        domain,
        solids,
        setup,
        laws,
        initial,
        config,
        stage,
        carried_power_w,
        gate,
    )?;
    let coordinates = coordinates.ok_or_else(|| {
        invalid("caloric tensor solve did not retain its original energy coordinate")
    })?;
    let mut specific_energy_j_kg = Vec::with_capacity(coordinates.len());
    for (cell, value) in coordinates.into_iter().enumerate() {
        if cell.is_multiple_of(512) {
            poll(gate)?;
        }
        specific_energy_j_kg.push(stage.specific_energy(cell, value.value)?);
    }
    Ok(TensorCaloricEndpoint {
        solution,
        specific_energy_j_kg,
    })
}

fn stage_record(
    domain: &VoxelDomain,
    setup: &ThermalSetup,
    storage: &PreparedCaloric<'_>,
    endpoint: &TensorCaloricEndpoint,
    time_s: f64,
    interval_s: f64,
) -> Result<TransientRecord, ChtError> {
    record_conduction_step(
        domain,
        setup,
        &endpoint.solution.solution.temperature,
        endpoint.solution.solution.report.balance,
        endpoint.solution.solution.report.iterations,
        |cell| storage.stored_j(cell, endpoint.specific_energy_j_kg[cell]),
        (
            time_s,
            interval_s,
            (0..domain.cell_count()).any(|cell| !domain.is_fluid(cell)),
        ),
        0.0,
        0.0,
    )
}

fn end_time(start: f64, interval: f64) -> Result<f64, ChtError> {
    finite("tensor.caloric.start_time_s", start)?;
    finite_positive("tensor.caloric.time_step_s", interval)?;
    let end = start + interval;
    finite("tensor.caloric.end_time_s", end)?;
    if start < 0.0 || end <= start {
        return Err(invalid(
            "caloric tensor time requires a nonnegative start and representable forward interval",
        ));
    }
    Ok(end)
}

fn source_error(error: DepositionError) -> ChtError {
    match error {
        DepositionError::Cancelled => ChtError::Cancelled,
        other => ChtError::InvalidInput {
            // Only a trial outside the constitutive chart may be backtracked.
            // Work, contract, arithmetic and heat refusals are not new seeds.
            field: if matches!(&other, DepositionError::MaterialDomain) {
                "thermal.tensor_caloric.original_source"
            } else {
                "thermal.tensor_caloric.original_source_refusal"
            },
            reason: other.to_string(),
        },
    }
}

impl TensorCaloricLaw for fs_material::phase::EquilibriumEnthalpyPhaseCurve {
    fn specific_energy_bounds_j_kg(&self) -> [f64; 2] {
        let knots = self.knots();
        [
            knots[0].specific_enthalpy_j_kg,
            knots[knots.len() - 1].specific_enthalpy_j_kg,
        ]
    }

    fn temperature_slope_kg_k_j(&self, energy: f64) -> Result<f64, DepositionError> {
        self.temperature_derivative_at_specific_enthalpy(energy)
            .map_err(|_| DepositionError::MaterialDomain)
    }
}
