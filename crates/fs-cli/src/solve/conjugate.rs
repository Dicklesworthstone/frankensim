//! Conjugate airflow exchange for the conduction stage (bead
//! frankensim-s93ej.3): each declared vent branch feeds its own ordered air
//! path, and ALL paths exchange heat through the same conduction solve.
//! Coefficients come from the branch Reynolds number and a validity-gated
//! fs-convection card; Robin reference temperatures come from the air march.
//!
//! Branches have independent prescribed inlets and frozen mass flows. There
//! is no inferred air mixing, recirculation, buoyancy, or momentum feedback.
//! Dry-air transport properties remain frozen at 300 K, 1 atm; density is the
//! envelope estimate retained by the flow-network stage. Flow uncertainty is
//! disclosed, not propagated through the thermal fixed point.
//!
//! The parent stage still owns meshing, the solid solve, and ledger/checkpoint
//! discipline. This adapter preserves its flat target/reference seam, while
//! retaining branch identity in the solver and receipt instead of collapsing
//! separate inlets or flows into a fictional single path.
//!
//! Global IQN-ILS accelerates the complete branch-major interface. Its bounded
//! history, rank filter, and startup relaxation are recorded in the receipt;
//! the temperature, branch watt, and decomposition gates remain independent.

use std::collections::{BTreeMap, BTreeSet};

use fs_airflow::conjugate::{
    AirPath, AirSegment, ConjugateConfig, ConjugateSolution, IqnIlsConfig,
    Relaxation, SolidRegionState,
};
use fs_airflow::graph::thermal::{ConjugateBranchesSolution, solve_conjugate_branches_iqn};
use fs_airflow::{AirflowError, OperatingPoint};
use fs_convection::{CorrelationId, ThermalConductivity, evaluate};
use fs_exec::Cx;
use fs_project::{ConductionSetup, ThermalBoundaryCondition};
use fs_qty::{Area, Density, DynViscosity, Length};

use super::{SolveRefusal, canonical_f64, conduction_error};
use crate::import::json_string;

/// Dry-air dynamic viscosity at 300 K, 1 atm (Pa·s).
pub(super) const AIR_DYNAMIC_VISCOSITY_PA_S: f64 = 1.846e-5;
/// Dry-air thermal conductivity at 300 K, 1 atm (W/(m·K)).
pub(super) const AIR_THERMAL_CONDUCTIVITY_W_M_K: f64 = 26.3e-3;
/// Dry-air Prandtl number at 300 K.
pub(super) const AIR_PRANDTL: f64 = 0.707;
/// Dry-air specific heat at 300 K (J/(kg·K)).
pub(super) const AIR_SPECIFIC_HEAT_J_KG_K: f64 = 1007.0;
/// The provenance every receipt cites for the frozen transport properties.
pub(super) const AIR_PROPERTY_SOURCE: &str = "dry air at 300 K, 1 atm (Incropera & DeWitt, Fundamentals of Heat and Mass Transfer, Table A.4); frozen across the exchange, not re-evaluated at film or bulk temperature";

// One policy supplies BOTH the driver and receipt; no environment-dependent
// selection or hidden solver retry can change the iteration/memory budget.
const CONJUGATE_IQN_CONFIG: IqnIlsConfig = IqnIlsConfig {
    max_history: 8,
    relative_rank_tolerance: 1.0e-10,
};
const CONJUGATE_STARTUP_OMEGA: f64 = 1.0;

/// Receipt-level authority statement for the exchange.
pub(super) const CONJUGATE_AUTHORITY: &str = "partitioned solid/air fixed point over derived Robin rows: card-derived coefficient at the flow-network midpoint, exponential-law marched reference temperature, kelvin convergence plus an independent watt balance gate and an fs-conduction decomposition cross-check";

/// Retained single-branch receipt statement, unchanged for existing projects.
pub(super) const CONJUGATE_NO_CLAIM: &str = "one branch only; air properties and the coefficient are frozen across the exchange; the air is a 1-D stream-wise chain with no recirculation, buoyancy, redistribution, or momentum feedback; the flow bracket width is disclosed, not propagated; the card's model-form discrepancy band is an engineering allowance, not a validated interval; no experimental validation and no maturity claim";

const BRANCHES_NO_CLAIM: &str = "independently supplied 1-D air branches coupled through one shared solid; no air mixing, recirculation, buoyancy, flow redistribution, or momentum feedback; air properties and coefficients are frozen; flow and card discrepancy bands are disclosed, not propagated; branch watt gates are separate, while the final conduction decomposition cross-check is aggregate only; no experimental validation, interval certification, or maturity claim";

/// One lowered airflow-convection law.
#[derive(Debug, Clone, PartialEq)]
pub(super) struct AirflowLaw {
    pub target: String,
    pub branch: String,
    pub order: u32,
    pub inlet_temperature_k: f64,
    pub hydraulic_diameter_m: f64,
    pub flow_area_m2: f64,
    pub channel_length_m: f64,
    pub correlation: CorrelationId,
}

/// Lower the project laws into canonical branch-major, stream-wise order.
/// Inlet agreement and order uniqueness apply WITHIN a branch, not across
/// independent branches. A solid target still has exactly one Robin owner.
pub(super) fn airflow_laws(setup: &ConductionSetup) -> Result<Vec<AirflowLaw>, SolveRefusal> {
    let mut laws = Vec::new();
    for boundary in &setup.boundaries {
        let ThermalBoundaryCondition::AirflowConvection {
            branch,
            order,
            inlet_temperature,
            hydraulic_diameter,
            flow_area,
            channel_length,
            correlation,
        } = &boundary.condition
        else {
            continue;
        };
        let card = CorrelationId::ALL
            .iter()
            .copied()
            .find(|id| id.name() == correlation)
            .ok_or_else(|| {
                conduction_error(
                    "cli-solve-conduction-airflow-correlation",
                    format!(
                        "airflow convection on `{}` names `{correlation}`, which is not an fs-convection card",
                        boundary.target
                    ),
                    "name one card from the fs-convection catalog by its `convection.*` id",
                )
            })?;
        laws.push(AirflowLaw {
            target: boundary.target.clone(),
            branch: branch.clone(),
            order: *order,
            inlet_temperature_k: inlet_temperature.value,
            hydraulic_diameter_m: hydraulic_diameter.value,
            flow_area_m2: flow_area.value,
            channel_length_m: channel_length.value,
            correlation: card,
        });
    }
    order_laws(&mut laws)?;
    Ok(laws)
}

fn order_laws(laws: &mut [AirflowLaw]) -> Result<(), SolveRefusal> {
    laws.sort_by(|a, b| (&a.branch, a.order).cmp(&(&b.branch, b.order)));
    let mut targets = BTreeSet::new();
    for law in laws.iter() {
        if law.branch.trim().is_empty() || law.target.trim().is_empty() {
            return Err(conduction_error(
                "cli-solve-conduction-airflow-target",
                "an airflow law has an empty branch or target",
                "name the vent branch and a nonempty solid boundary target",
            ));
        }
        if !targets.insert(law.target.as_str()) {
            return Err(conduction_error(
                "cli-solve-conduction-airflow-target",
                format!("airflow target `{}` is owned by more than one law", law.target),
                "assign each solid Robin target to exactly one branch and segment",
            ));
        }
        for (what, value) in [
            ("inlet temperature", law.inlet_temperature_k),
            ("hydraulic diameter", law.hydraulic_diameter_m),
            ("flow area", law.flow_area_m2),
            ("channel length", law.channel_length_m),
        ] {
            if !(value.is_finite() && value > 0.0) {
                return Err(conduction_error(
                    "cli-solve-conduction-airflow-range",
                    format!("airflow convection {what} on `{}` is {value}", law.target),
                    "declare a finite positive quantity",
                ));
            }
        }
    }
    for pair in laws.windows(2) {
        if pair[0].branch != pair[1].branch {
            continue;
        }
        if pair[0].inlet_temperature_k.to_bits() != pair[1].inlet_temperature_k.to_bits() {
            return Err(conduction_error(
                "cli-solve-conduction-airflow-inlet",
                format!(
                    "airflow convection on `{}` declares inlet {} K while `{}` declares {} K on the same branch",
                    pair[0].target, pair[0].inlet_temperature_k,
                    pair[1].target, pair[1].inlet_temperature_k
                ),
                "one branch has one inlet temperature; declare it identically on every law of the branch",
            ));
        }
        if pair[0].order == pair[1].order {
            return Err(conduction_error(
                "cli-solve-conduction-airflow-order",
                format!(
                    "airflow convection on `{}` and `{}` share stream-wise order {} on branch `{}`",
                    pair[0].target, pair[1].target, pair[0].order, pair[0].branch
                ),
                "give every law on one branch a distinct stream-wise order",
            ));
        }
    }
    Ok(())
}

/// Per-segment derivation evidence retained in the receipt.
#[derive(Debug, Clone, PartialEq)]
pub(super) struct SegmentDerivation {
    pub target: String,
    pub order: u32,
    pub wetted_area_m2: f64,
    pub velocity_m_s: f64,
    pub reynolds: f64,
    pub length_over_hydraulic_diameter: f64,
    pub nusselt: f64,
    pub htc_w_m2_k: f64,
    pub card: CorrelationId,
    pub in_domain: bool,
}

#[derive(Debug, Clone)]
struct ConjugateBranch {
    branch: String,
    path_name: String,
    flow_mid_m3_s: f64,
    flow_lo_m3_s: f64,
    flow_hi_m3_s: f64,
    air_density_kg_m3: f64,
    mass_flow_kg_s: f64,
    inlet_temperature_k: f64,
    segments: Vec<SegmentDerivation>,
    air_path: AirPath,
}

/// Flat solid-side ordering paired with the independent air-side branches.
/// No aggregate inlet temperature or representative mass flow is invented.
#[derive(Debug, Clone)]
pub(super) struct ConjugatePath {
    pub segments: Vec<SegmentDerivation>,
    branches: Vec<ConjugateBranch>,
}

impl ConjugatePath {
    /// Preserve the independently supplied paths in the same branch-major
    /// order as the solid's Robin rows for the coupled goal linearization.
    pub(super) fn air_paths(&self) -> Vec<AirPath> {
        self.branches.iter().map(|branch| branch.air_path.clone()).collect()
    }
}

fn exchange_config(adaptive: bool) -> ConjugateConfig {
    let mut config = ConjugateConfig {
        relaxation: Relaxation::Fixed { omega: CONJUGATE_STARTUP_OMEGA },
        ..ConjugateConfig::default()
    };
    if adaptive {
        // The goal comparison independently checks the fully coupled solid
        // residual. Tighten the reference solve without weakening its watt gate.
        config.temperature_tolerance_k = 1.0e-11;
    }
    config
}

pub(super) fn goal_config(linear: fs_conduction::LinearConfig)
    -> fs_airflow::conjugate::goal::CoupledGoalConfig
{
    use fs_airflow::conjugate::goal::{CoupledGoalConfig, InterfaceSolveConfig};
    CoupledGoalConfig {
        solid: linear,
        primal: exchange_config(true),
        interface: InterfaceSolveConfig {
            max_iterations: 32,
            absolute_tolerance: 1.0e-12,
            relative_tolerance: 1.0e-12,
            relaxation: 1.0,
        },
        acceleration: CONJUGATE_IQN_CONFIG,
    }
}

/// Derive every branch at its own solved flow, then retain the same ordering
/// on the flat solid side and the branch-major air side.
///
/// `htc_scale` multiplies every card-derived coefficient before the air path
/// is built, so a model-form propagation perturbs the coefficient the
/// conjugate fixed point actually uses. The nominal derivation passes `1.0`.
pub(super) fn derive_air_path(
    laws: &[AirflowLaw],
    operating: &OperatingPoint,
    air_density_kg_m3: f64,
    wetted_area_m2: impl Fn(&str) -> Option<f64>,
    htc_scale: f64,
) -> Result<ConjugatePath, SolveRefusal> {
    if laws.is_empty() {
        return Err(conduction_error(
            "cli-solve-conduction-airflow-empty",
            "no airflow-convection law to derive",
            "report the driver defect; the exchange runs only when a law is declared",
        ));
    }
    let mut laws = laws.to_vec();
    order_laws(&mut laws)?;
    let mut branches = Vec::new();
    let mut segments = Vec::with_capacity(laws.len());
    let mut start = 0;
    while start < laws.len() {
        let mut end = start + 1;
        while end < laws.len() && laws[end].branch == laws[start].branch {
            end += 1;
        }
        let branch = derive_branch(
            &laws[start..end], operating, air_density_kg_m3, &wetted_area_m2, htc_scale,
        )?;
        segments.extend(branch.segments.iter().cloned());
        branches.push(branch);
        start = end;
    }
    Ok(ConjugatePath { segments, branches })
}

/// Each branch uses its own `vent:<name>` handoff for every channel segment.
fn derive_branch(
    laws: &[AirflowLaw],
    operating: &OperatingPoint,
    air_density_kg_m3: f64,
    wetted_area_m2: impl Fn(&str) -> Option<f64>,
    htc_scale: f64,
) -> Result<ConjugateBranch, SolveRefusal> {
    let first = &laws[0];
    if !(air_density_kg_m3.is_finite() && air_density_kg_m3 > 0.0) {
        return Err(conduction_error(
            "cli-solve-conduction-airflow-density",
            format!("the flow-network air density estimate is {air_density_kg_m3} kg/m^3"),
            "check the envelope ambient temperature and pressure",
        ));
    }
    let path_name = format!("vent:{}", first.branch);
    let mut segments = Vec::with_capacity(laws.len());
    let mut air_segments = Vec::with_capacity(laws.len());
    let mut mass_flow = None;
    let mut flow_bracket = None;
    for law in laws {
        let handoff = operating
            .correlation_handoff(
                &path_name,
                Area::new(law.flow_area_m2),
                Density::new(air_density_kg_m3),
                DynViscosity::new(AIR_DYNAMIC_VISCOSITY_PA_S),
                Length::new(law.hydraulic_diameter_m),
                AIR_PRANDTL,
            )
            .map_err(|error| {
                conduction_error(
                    "cli-solve-conduction-airflow-handoff",
                    format!(
                        "airflow convection on `{}` cannot take the branch `{}` handoff: {error}",
                        law.target, law.branch
                    ),
                    "declare the law on a vent the flow network actually solved",
                )
            })?;
        let ratio = law.channel_length_m / law.hydraulic_diameter_m;
        let inputs = handoff.correlation_inputs.with_length_ratio(ratio);
        let nusselt = evaluate(law.correlation, inputs).map_err(|error| {
            conduction_error(
                "cli-solve-conduction-airflow-card-domain",
                format!(
                    "card `{}` refuses the branch `{}` regime on `{}` (Re {:.4e}, L/Dh {:.4e}): {error}",
                    law.correlation.name(), law.branch, law.target, handoff.reynolds, ratio
                ),
                "pick a card whose declared domain covers this operating point, or change the channel declaration",
            )
        })?;
        let in_domain = nusselt.evidence().model.in_domain;
        if !in_domain {
            return Err(conduction_error(
                "cli-solve-conduction-airflow-card-domain",
                format!(
                    "card `{}` evaluated outside its declared domain on `{}`",
                    law.correlation.name(), law.target
                ),
                "the derivation never extrapolates a card; pick one whose domain covers the point",
            ));
        }
        let htc = nusselt
            .heat_transfer_coefficient(
                ThermalConductivity::new(AIR_THERMAL_CONDUCTIVITY_W_M_K),
                Length::new(law.hydraulic_diameter_m),
            )
            .map_err(|error| {
                conduction_error(
                    "cli-solve-conduction-airflow-card",
                    format!(
                        "card `{}` cannot lower Nu to a coefficient on `{}`: {error}",
                        law.correlation.name(), law.target
                    ),
                    "report the card defect",
                )
            })?
            .value
            .value()
            * htc_scale;
        let wetted = wetted_area_m2(&law.target).ok_or_else(|| {
            conduction_error(
                "cli-solve-conduction-airflow-area",
                format!("airflow convection target `{}` owns no exterior face area", law.target),
                "select a nonempty exterior face set for every airflow-convection target",
            )
        })?;
        if !(wetted.is_finite() && wetted > 0.0) {
            return Err(conduction_error(
                "cli-solve-conduction-airflow-area",
                format!("airflow convection target `{}` has wetted area {wetted} m^2", law.target),
                "the retained mesh must give every airflow target positive exterior area",
            ));
        }
        let branch_flow = handoff.branch_flow.value.value();
        let this_mass_flow = air_density_kg_m3 * branch_flow;
        match mass_flow {
            None => {
                mass_flow = Some(this_mass_flow);
                flow_bracket = Some((
                    branch_flow,
                    handoff.branch_flow.numerical.lo,
                    handoff.branch_flow.numerical.hi,
                ));
            }
            Some(previous) if previous.to_bits() != this_mass_flow.to_bits() => {
                return Err(conduction_error(
                    "cli-solve-conduction-airflow-handoff",
                    format!("branch `{}` handed off two different flows to laws on the same path", law.branch),
                    "report the driver defect; one branch has one solved flow",
                ));
            }
            Some(_) => {}
        }
        air_segments.push(AirSegment::new(&law.target, wetted, htc).map_err(|error| {
            conduction_error(
                "cli-solve-conduction-airflow-segment",
                format!("air segment `{}` refused (area {wetted} m^2, h {htc} W/m^2K): {error}", law.target),
                "check the derived coefficient and the target's exterior area",
            )
        })?);
        segments.push(SegmentDerivation {
            target: law.target.clone(),
            order: law.order,
            wetted_area_m2: wetted,
            velocity_m_s: handoff.velocity.value.value(),
            reynolds: handoff.reynolds,
            length_over_hydraulic_diameter: ratio,
            nusselt: nusselt.evidence().value,
            htc_w_m2_k: htc,
            card: law.correlation,
            in_domain,
        });
    }
    let mass_flow_kg_s = mass_flow.expect("at least one law");
    let (flow_mid, flow_lo, flow_hi) = flow_bracket.expect("at least one law");
    let air_path = AirPath::new(
        first.inlet_temperature_k,
        mass_flow_kg_s,
        AIR_SPECIFIC_HEAT_J_KG_K,
        air_segments,
    )
    .map_err(|error| {
        conduction_error(
            "cli-solve-conduction-airflow-path",
            format!("the branch `{}` air path refused: {error}", first.branch),
            "check the inlet temperature, branch flow, and segment declarations",
        )
    })?;
    Ok(ConjugateBranch {
        branch: first.branch.clone(),
        path_name,
        flow_mid_m3_s: flow_mid,
        flow_lo_m3_s: flow_lo,
        flow_hi_m3_s: flow_hi,
        air_density_kg_m3,
        mass_flow_kg_s,
        inlet_temperature_k: first.inlet_temperature_k,
        segments,
        air_path,
    })
}

/// The solid callback consumes one coefficient/reference pair per target.
pub(super) fn derived_coefficients(path: &ConjugatePath) -> BTreeMap<String, f64> {
    path.segments.iter()
        .map(|segment| (segment.target.clone(), segment.htc_w_m2_k))
        .collect()
}

#[derive(Debug, Clone)]
pub(super) struct ConjugateOutcome {
    pub solution: ConjugateBranchesSolution,
    pub balance_tolerance_w: f64,
    decomposition_residual_w: Option<f64>,
}

/// One common solid callback per iteration, with every branch's Robin rows.
/// The returned region sequence uses the same flat ordering as the parent
/// stage's boundary partition. Solid-side typed refusals propagate unchanged.
pub(super) fn run_exchange(
    cx: &Cx<'_>,
    path: &ConjugatePath,
    adaptive: bool,
    mut solid: impl FnMut(
        &Cx<'_>, &BTreeMap<String, f64>,
    ) -> Result<Vec<SolidRegionState>, SolveRefusal>,
) -> Result<ConjugateOutcome, SolveRefusal> {
    let config = exchange_config(adaptive);
    let mut stashed: Option<SolveRefusal> = None;
    let targets: Vec<String> = path.segments.iter().map(|s| s.target.clone()).collect();
    let air_paths = path.air_paths();
    let result = solve_conjugate_branches_iqn(cx, &air_paths, &config, CONJUGATE_IQN_CONFIG, |cx, references| {
        let by_target: BTreeMap<String, f64> = targets.iter().cloned()
            .zip(references.iter().copied()).collect();
        match solid(cx, &by_target) {
            Ok(states) => Ok(states),
            Err(refusal) => {
                stashed = Some(refusal);
                Err(AirflowError::Cancelled {
                    iteration: 0, references_k: references.to_vec(),
                })
            }
        }
    });
    if let Some(refusal) = stashed {
        return Err(refusal);
    }
    let solution = result.map_err(|error| match error {
        AirflowError::ConjugateBalanceUnclosed {
            iterations, max_region_imbalance_bits, tolerance_bits,
        } => conduction_error(
            "cli-solve-conduction-airflow-balance",
            format!(
                "a branch's solid/air exchange met its temperature criterion but heat rates disagree by {:.4e} W against its {:.4e} W gate (probe iterations {iterations}); another branch's power cannot relax this gate",
                f64::from_bits(max_region_imbalance_bits), f64::from_bits(tolerance_bits)
            ),
            "inspect the boundary partition, wetted areas, and each branch's mass flow",
        ),
        AirflowError::ConjugateNotConverged { .. } => conduction_error(
            "cli-solve-conduction-airflow-unconverged",
            format!("the solid/air fixed point did not converge: {error}"),
            "check the segment declarations; a runaway reference usually means an inconsistent wetted area or flow",
        ),
        AirflowError::Cancelled { .. } => conduction_error(
            "cli-solve-conduction-airflow-cancelled",
            "the solid/air exchange was cancelled".to_string(),
            "rerun or resume the stage",
        ),
        other => conduction_error(
            "cli-solve-conduction-airflow-exchange",
            format!("the solid/air exchange refused: {other}"),
            "inspect the derived air paths and the solid response",
        ),
    })?;
    // Sum branch-specific absolute thresholds for the aggregate cross-check.
    // Signed heat cancellation must not make this threshold spuriously tiny.
    // The solver already enforced each branch's OWN region-level watt gate.
    let balance_tolerance_w = finite_total(
        solution.branches.iter().map(branch_balance_tolerance),
        "aggregate balance tolerance",
    )?;
    Ok(ConjugateOutcome { solution, balance_tolerance_w, decomposition_residual_w: None })
}

fn branch_balance_tolerance(solution: &ConjugateSolution) -> f64 {
    let config = ConjugateConfig::default();
    let scale = solution.balance.solid_total_w.abs().max(solution.balance.air_total_w.abs());
    config.balance_tolerance_w.max(config.balance_relative_tolerance * scale)
}

/// Gate the final published solve's independently accumulated Robin total.
/// The parent reports a whole-domain total, so this is an AGGREGATE check;
/// multi-branch receipts deliberately leave branch decomposition checks null.
pub(super) fn cross_check_decomposition(
    mut outcome: ConjugateOutcome,
    final_robin_out_total_w: f64,
    final_robin_out_off_path_w: f64,
) -> Result<ConjugateOutcome, SolveRefusal> {
    let on_path = final_robin_out_total_w - final_robin_out_off_path_w;
    let solid_total = finite_total(
        outcome.solution.branches.iter().map(|b| b.balance.solid_total_w),
        "aggregate solid heat rate",
    )?;
    let residual = solid_total - on_path;
    if !(residual.is_finite() && residual.abs() <= outcome.balance_tolerance_w) {
        return Err(conduction_error(
            "cli-solve-conduction-airflow-decomposition",
            format!(
                "the paths' summed solid heat rate differs from fs-conduction's own Robin accumulation by {residual:.4e} W against a {:.4e} W gate",
                outcome.balance_tolerance_w
            ),
            "a Robin face may be unowned or counted twice; inspect the boundary partition",
        ));
    }
    outcome.decomposition_residual_w = Some(residual);
    if outcome.solution.branches.len() == 1 {
        outcome.solution.branches[0].balance.decomposition_residual_w = Some(residual);
    }
    Ok(outcome)
}

fn num(value: f64, what: &str) -> Result<String, SolveRefusal> {
    canonical_f64(value).ok_or_else(|| {
        conduction_error(
            "cli-solve-conduction-nonfinite",
            format!("conjugate exchange value `{what}` is {value}"),
            "report the solver defect; non-finite exchange values are never published",
        )
    })
}

fn finite_total(values: impl IntoIterator<Item = f64>, what: &str) -> Result<f64, SolveRefusal> {
    let mut total = 0.0;
    for value in values {
        total += value;
        num(total, what)?;
    }
    Ok(total)
}

fn acceleration_receipt_fragment() -> Result<String, SolveRefusal> {
    Ok(format!(
        "{{\"method\":\"iqn-ils\",\"scope\":\"all-branch-interfaces\",\"max_history\":{},\"relative_rank_tolerance\":{},\"fallback\":\"fixed\",\"fallback_omega\":{}}}",
        CONJUGATE_IQN_CONFIG.max_history,
        num(CONJUGATE_IQN_CONFIG.relative_rank_tolerance, "IQN rank tolerance")?,
        num(CONJUGATE_STARTUP_OMEGA, "IQN fallback omega")?,
    ))
}

/// Single-branch receipts retain their branch fields; both forms now record
/// the acceleration policy. Multiple branches use an explicit tagged object,
/// never a misleading scalar inlet, outlet, or mass flow. The aggregate must
/// have passed the publication gate.
pub(super) fn receipt_fragment(
    path: &ConjugatePath,
    outcome: &ConjugateOutcome,
) -> Result<String, SolveRefusal> {
    if path.branches.is_empty() || path.branches.len() != outcome.solution.branches.len() {
        return Err(conduction_error(
            "cli-solve-conduction-airflow-exchange",
            "the conjugate receipt does not have one solution per declared branch",
            "report the driver defect",
        ));
    }
    let residual = outcome.decomposition_residual_w.ok_or_else(|| conduction_error(
        "cli-solve-conduction-airflow-decomposition",
        "the final shared-solid decomposition has not been checked",
        "cross-check the final published conduction solve before rendering the receipt",
    ))?;
    if path.branches.len() == 1 {
        return branch_receipt_fragment(
            &path.branches[0], &outcome.solution.branches[0],
            outcome.balance_tolerance_w, Some(residual), CONJUGATE_NO_CLAIM,
        );
    }
    let mut fragments = Vec::with_capacity(path.branches.len());
    for (branch, solution) in path.branches.iter().zip(&outcome.solution.branches) {
        fragments.push(branch_receipt_fragment(
            branch, solution, branch_balance_tolerance(solution), None, BRANCHES_NO_CLAIM,
        )?);
    }
    let solid_total = finite_total(outcome.solution.branches.iter().map(|b| b.balance.solid_total_w), "solid_total_w")?;
    let air_total = finite_total(outcome.solution.branches.iter().map(|b| b.balance.air_total_w), "air_total_w")?;
    Ok(format!(
        "{{\"schema\":\"independent-branches-shared-solid-v1\",\"branch_count\":{},\"branches\":[{}],\"iterations\":{},\"solid_total_w\":{},\"air_total_w\":{},\"interface_imbalance_w\":{},\"decomposition_residual_w\":{},\"balance_tolerance_w\":{},\"acceleration\":{},\"authority\":{},\"no_claim\":{}}}",
        path.branches.len(), fragments.join(","), outcome.solution.iterations,
        num(solid_total, "solid_total_w")?, num(air_total, "air_total_w")?,
        num(solid_total - air_total, "interface_imbalance_w")?,
        num(residual, "decomposition_residual_w")?,
        num(outcome.balance_tolerance_w, "balance_tolerance_w")?,
        acceleration_receipt_fragment()?,
        json_string(CONJUGATE_AUTHORITY), json_string(BRANCHES_NO_CLAIM),
    ))
}

fn branch_receipt_fragment(
    path: &ConjugateBranch,
    solution: &ConjugateSolution,
    balance_tolerance_w: f64,
    decomposition_residual_w: Option<f64>,
    no_claim: &str,
) -> Result<String, SolveRefusal> {
    if path.segments.len() != solution.march.segments.len()
        || path.segments.len() != solution.balance.regions.len()
    {
        return Err(conduction_error(
            "cli-solve-conduction-airflow-exchange",
            "the branch receipt has mismatched derivation, march, or balance rows",
            "report the driver defect",
        ));
    }
    let mut segments = String::new();
    for (index, (derivation, state)) in path.segments.iter()
        .zip(solution.march.segments.iter()).enumerate()
    {
        if index > 0 { segments.push(','); }
        let balance = &solution.balance.regions[index];
        if state.region != derivation.target || balance.region != derivation.target {
            return Err(conduction_error(
                "cli-solve-conduction-airflow-exchange",
                "a branch receipt's segment identities do not agree",
                "report the driver defect",
            ));
        }
        segments.push_str(&format!(
            "{{\"target\":{},\"order\":{},\"card\":{},\"in_domain\":{},\"wetted_area_m2\":{},\"velocity_m_s\":{},\"reynolds\":{},\"length_over_hydraulic_diameter\":{},\"nusselt\":{},\"htc_w_m2_k\":{},\"air_in_k\":{},\"air_out_k\":{},\"reference_k\":{},\"ntu\":{},\"effectiveness\":{},\"solid_heat_rate_w\":{},\"air_heat_rate_w\":{},\"imbalance_w\":{}}}",
            json_string(&derivation.target), derivation.order,
            json_string(derivation.card.name()), derivation.in_domain,
            num(derivation.wetted_area_m2, "wetted_area_m2")?,
            num(derivation.velocity_m_s, "velocity_m_s")?,
            num(derivation.reynolds, "reynolds")?,
            num(derivation.length_over_hydraulic_diameter, "length_over_hydraulic_diameter")?,
            num(derivation.nusselt, "nusselt")?, num(derivation.htc_w_m2_k, "htc_w_m2_k")?,
            num(state.inlet_temperature_k, "air_in_k")?, num(state.outlet_temperature_k, "air_out_k")?,
            num(state.reference_temperature_k, "reference_k")?, num(state.ntu, "ntu")?,
            num(state.effectiveness, "effectiveness")?,
            num(balance.solid_heat_rate_w, "solid_heat_rate_w")?,
            num(balance.air_heat_rate_w, "air_heat_rate_w")?, num(balance.imbalance_w, "imbalance_w")?,
        ));
    }
    let audit = &solution.balance;
    let decomposition = match decomposition_residual_w {
        Some(value) => num(value, "decomposition_residual_w")?,
        None => "null".to_string(),
    };
    Ok(format!(
        "{{\"branch\":{},\"path\":{},\"flow_m3_s\":{{\"lo\":{},\"mid\":{},\"hi\":{}}},\"air_density_kg_m3\":{},\"mass_flow_kg_s\":{},\"inlet_k\":{},\"outlet_k\":{},\"air_properties\":{{\"dynamic_viscosity_pa_s\":{},\"thermal_conductivity_w_m_k\":{},\"prandtl\":{},\"specific_heat_j_kg_k\":{},\"source\":{}}},\"segments\":[{}],\"iterations\":{},\"solid_total_w\":{},\"air_total_w\":{},\"interface_imbalance_w\":{},\"max_region_imbalance_w\":{},\"worst_recorded_imbalance_w\":{},\"decomposition_residual_w\":{},\"balance_tolerance_w\":{},\"acceleration\":{},\"authority\":{},\"no_claim\":{}}}",
        json_string(&path.branch), json_string(&path.path_name),
        num(path.flow_lo_m3_s, "flow_lo")?, num(path.flow_mid_m3_s, "flow_mid")?,
        num(path.flow_hi_m3_s, "flow_hi")?, num(path.air_density_kg_m3, "air_density")?,
        num(path.mass_flow_kg_s, "mass_flow")?, num(path.inlet_temperature_k, "inlet_k")?,
        num(solution.march.outlet_temperature_k, "outlet_k")?,
        num(AIR_DYNAMIC_VISCOSITY_PA_S, "mu")?, num(AIR_THERMAL_CONDUCTIVITY_W_M_K, "k")?,
        num(AIR_PRANDTL, "Pr")?, num(AIR_SPECIFIC_HEAT_J_KG_K, "cp")?,
        json_string(AIR_PROPERTY_SOURCE), segments, solution.iterations,
        num(audit.solid_total_w, "solid_total_w")?, num(audit.air_total_w, "air_total_w")?,
        num(audit.interface_imbalance_w, "interface_imbalance_w")?,
        num(audit.max_region_imbalance_w, "max_region_imbalance_w")?,
        num(solution.worst_recorded_imbalance_w, "worst_recorded_imbalance_w")?,
        decomposition, num(balance_tolerance_w, "balance_tolerance_w")?,
        acceleration_receipt_fragment()?,
        json_string(CONJUGATE_AUTHORITY), json_string(no_claim),
    ))
}

#[cfg(test)]
mod tests;
