//! Conjugate airflow exchange for the conduction stage (bead
//! frankensim-s93ej.3): the solved flow-network operating point feeds a
//! declared `(airflow-convection ...)` boundary law through an
//! `fs-convection` card into a stream-wise `fs_airflow::conjugate` air path,
//! and the solid/air fixed point is closed against the real heterogeneous
//! conduction solve. The Robin rows the conduction operator finally sees are
//! DERIVED here (coefficient from the card at the branch Reynolds number,
//! reference temperature from the marched air), never declared.
//!
//! What this module owns: lowering the laws, deriving the air path with its
//! per-segment card evidence, running the partitioned exchange, mapping its
//! typed refusals, and rendering the receipt fragment. What it does not own:
//! the conduction solve itself (the caller supplies it as the solid response)
//! and the stage's ledger/checkpoint discipline.
//!
//! # No-claim boundaries (carried verbatim into the receipt)
//!
//! * ONE conjugate branch per project in this slice. Laws naming two or more
//!   branches refuse by name rather than being solved as an unmodeled
//!   multi-branch network.
//! * Air transport properties are dry air at 300 K, 1 atm, frozen across the
//!   exchange; density alone comes from the envelope ideal-gas estimate the
//!   flow-network stage already publishes. `h` is therefore frozen too
//!   (`fs_airflow::conjugate`'s own boundary).
//! * The air is a 1-D chain of well-mixed segments in declared stream-wise
//!   order over the declared channel geometry; wetted area is the retained
//!   mesh's exterior face area of the target, not a declared value.
//! * The flow used is the interval-certified nominal root's MIDPOINT; the
//!   bracket width is disclosed, not propagated through the card.

use std::collections::BTreeMap;

use fs_airflow::conjugate::{
    AirPath, AirSegment, ConjugateConfig, ConjugateSolution, SolidRegionState, solve_conjugate,
};
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

/// Receipt-level authority statement for the exchange.
pub(super) const CONJUGATE_AUTHORITY: &str = "partitioned solid/air fixed point over derived Robin rows: card-derived coefficient at the flow-network midpoint, exponential-law marched reference temperature, kelvin convergence plus an independent watt balance gate and an fs-conduction decomposition cross-check";

/// Receipt-level no-claim statement for the exchange.
pub(super) const CONJUGATE_NO_CLAIM: &str = "one branch only; air properties and the coefficient are frozen across the exchange; the air is a 1-D stream-wise chain with no recirculation, buoyancy, redistribution, or momentum feedback; the flow bracket width is disclosed, not propagated; the card's model-form discrepancy band is an engineering allowance, not a validated interval; no experimental validation and no maturity claim";

/// One lowered airflow-convection law.
#[derive(Debug, Clone, PartialEq)]
pub(super) struct AirflowLaw {
    /// Boundary target (assignment target; the Robin region name).
    pub target: String,
    /// Declared vent region whose branch carries the air.
    pub branch: String,
    /// Stream-wise order on the branch.
    pub order: u32,
    /// Declared inlet temperature, K.
    pub inlet_temperature_k: f64,
    /// Channel hydraulic diameter, m.
    pub hydraulic_diameter_m: f64,
    /// Channel free-flow area, m².
    pub flow_area_m2: f64,
    /// Channel stream-wise length, m.
    pub channel_length_m: f64,
    /// The named card.
    pub correlation: CorrelationId,
}

/// Lower every `AirflowConvection` law of the setup in stream-wise order,
/// enforcing this slice's single-branch boundary and a shared inlet.
///
/// Structural validation already checked branch existence, order
/// uniqueness, dimensions, and the card name; this re-derives the facts the
/// exchange needs and refuses (never repairs) anything it cannot honour.
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
        for (what, value) in [
            ("inlet temperature", inlet_temperature.value),
            ("hydraulic diameter", hydraulic_diameter.value),
            ("flow area", flow_area.value),
            ("channel length", channel_length.value),
        ] {
            if !(value.is_finite() && value > 0.0) {
                return Err(conduction_error(
                    "cli-solve-conduction-airflow-range",
                    format!(
                        "airflow convection {what} on `{}` is {value}",
                        boundary.target
                    ),
                    "declare a finite positive quantity",
                ));
            }
        }
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
    laws.sort_by(|a, b| (&a.branch, a.order).cmp(&(&b.branch, b.order)));
    if let Some(first) = laws.first() {
        if let Some(other) = laws.iter().find(|law| law.branch != first.branch) {
            return Err(conduction_error(
                "cli-solve-conduction-airflow-multi-branch",
                format!(
                    "airflow convection laws name branches `{}` and `{}`; this slice closes one conjugate branch per project",
                    first.branch, other.branch
                ),
                "declare every airflow-convection law on one vent branch, or split the study",
            ));
        }
        if let Some(other) = laws
            .iter()
            .find(|law| law.inlet_temperature_k.to_bits() != first.inlet_temperature_k.to_bits())
        {
            return Err(conduction_error(
                "cli-solve-conduction-airflow-inlet",
                format!(
                    "airflow convection on `{}` declares inlet {} K while `{}` declares {} K on the same branch",
                    first.target,
                    first.inlet_temperature_k,
                    other.target,
                    other.inlet_temperature_k
                ),
                "one branch has one inlet temperature; declare it identically on every law of the branch",
            ));
        }
        for pair in laws.windows(2) {
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
    }
    Ok(laws)
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

/// The lowered branch: everything the exchange consumes, plus what the
/// receipt discloses about how it was derived.
#[derive(Debug, Clone)]
pub(super) struct ConjugatePath {
    pub branch: String,
    pub path_name: String,
    pub flow_mid_m3_s: f64,
    pub flow_lo_m3_s: f64,
    pub flow_hi_m3_s: f64,
    pub air_density_kg_m3: f64,
    pub mass_flow_kg_s: f64,
    pub inlet_temperature_k: f64,
    pub segments: Vec<SegmentDerivation>,
    pub air_path: AirPath,
}

/// Derive the branch air path from the solved operating point: one card
/// evaluation per law at the branch Reynolds number, `h = Nu k / D_h`, and
/// the retained wetted area of each target.
///
/// `wetted_area_m2(target)` is the exterior face area the conduction
/// boundary partition owns for that target; `None` refuses.
pub(super) fn derive_air_path(
    laws: &[AirflowLaw],
    operating: &OperatingPoint,
    air_density_kg_m3: f64,
    wetted_area_m2: impl Fn(&str) -> Option<f64>,
) -> Result<ConjugatePath, SolveRefusal> {
    let first = laws.first().ok_or_else(|| {
        conduction_error(
            "cli-solve-conduction-airflow-empty",
            "no airflow-convection law to derive",
            "report the driver defect; the exchange runs only when a law is declared",
        )
    })?;
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
                    law.correlation.name(),
                    law.branch,
                    law.target,
                    handoff.reynolds,
                    ratio
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
                    law.correlation.name(),
                    law.target
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
                        law.correlation.name(),
                        law.target
                    ),
                    "report the card defect",
                )
            })?
            .value
            .value();
        let wetted = wetted_area_m2(&law.target).ok_or_else(|| {
            conduction_error(
                "cli-solve-conduction-airflow-area",
                format!(
                    "airflow convection target `{}` owns no exterior face area",
                    law.target
                ),
                "select a nonempty exterior face set for every airflow-convection target",
            )
        })?;
        if !(wetted.is_finite() && wetted > 0.0) {
            return Err(conduction_error(
                "cli-solve-conduction-airflow-area",
                format!(
                    "airflow convection target `{}` has wetted area {wetted} m^2",
                    law.target
                ),
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
                    format!(
                        "branch `{}` handed off two different flows to laws on the same path",
                        law.branch
                    ),
                    "report the driver defect; one branch has one solved flow",
                ));
            }
            Some(_) => {}
        }
        air_segments.push(AirSegment::new(&law.target, wetted, htc).map_err(|error| {
            conduction_error(
                "cli-solve-conduction-airflow-segment",
                format!(
                    "air segment `{}` refused (area {wetted} m^2, h {htc} W/m^2K): {error}",
                    law.target
                ),
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
    Ok(ConjugatePath {
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

/// Derived Robin coefficients per target, in path order, for the solid side
/// to lower into `fs_conduction::ThermalBc::robin(h, T_ref)`.
pub(super) fn derived_coefficients(path: &ConjugatePath) -> BTreeMap<String, f64> {
    path.segments
        .iter()
        .map(|segment| (segment.target.clone(), segment.htc_w_m2_k))
        .collect()
}

/// The closed exchange.
#[derive(Debug, Clone)]
pub(super) struct ConjugateOutcome {
    pub solution: ConjugateSolution,
    /// The effective watt tolerance the balance gate applied.
    pub balance_tolerance_w: f64,
}

/// Run the partitioned exchange against the caller's solid response, mapping
/// every `fs_airflow` refusal to a typed stage refusal. The solid response
/// receives the per-target reference temperatures and must return the
/// per-target Robin decomposition in PATH ORDER. A refusal raised by the
/// solid side itself (a conduction refusal, a cancellation) is returned
/// unchanged, never laundered through the airflow error.
pub(super) fn run_exchange(
    cx: &Cx<'_>,
    path: &ConjugatePath,
    mut solid: impl FnMut(
        &Cx<'_>,
        &BTreeMap<String, f64>,
    ) -> Result<Vec<SolidRegionState>, SolveRefusal>,
) -> Result<ConjugateOutcome, SolveRefusal> {
    let config = ConjugateConfig::default();
    let mut stashed: Option<SolveRefusal> = None;
    let targets: Vec<String> = path
        .segments
        .iter()
        .map(|segment| segment.target.clone())
        .collect();
    let result = solve_conjugate(cx, &path.air_path, &config, |cx, references| {
        let by_target: BTreeMap<String, f64> = targets
            .iter()
            .cloned()
            .zip(references.iter().copied())
            .collect();
        match solid(cx, &by_target) {
            Ok(states) => Ok(states),
            Err(refusal) => {
                stashed = Some(refusal);
                // The driver stops on this variant; the stashed refusal is
                // what the caller sees, so the payload here is never read.
                Err(AirflowError::Cancelled {
                    iteration: 0,
                    references_k: references.to_vec(),
                })
            }
        }
    });
    if let Some(refusal) = stashed {
        return Err(refusal);
    }
    let solution = result.map_err(|error| match error {
        AirflowError::ConjugateBalanceUnclosed {
            iterations,
            max_region_imbalance_bits,
            tolerance_bits,
        } => conduction_error(
            "cli-solve-conduction-airflow-balance",
            format!(
                "the solid/air exchange converged in temperature after {iterations} iteration(s) but a region's heat rates disagree by {:.4e} W against a {:.4e} W gate; a kelvin residual cannot bound a watt imbalance, so the result is refused",
                f64::from_bits(max_region_imbalance_bits),
                f64::from_bits(tolerance_bits)
            ),
            "inspect the boundary partition: a face owned by no path region, a wrong wetted area, or a mass-flow mismatch",
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
            "inspect the derived air path and the solid response",
        ),
    })?;
    // The same scale-free watt gate the driver applied per region also gates
    // the decomposition cross-check the stage runs on its published solve.
    let scale = solution
        .balance
        .solid_total_w
        .abs()
        .max(solution.balance.air_total_w.abs());
    let balance_tolerance_w = config
        .balance_tolerance_w
        .max(config.balance_relative_tolerance * scale);
    Ok(ConjugateOutcome {
        solution,
        balance_tolerance_w,
    })
}

/// Apply the decomposition cross-check against the FINAL published solve's
/// total Robin heat rate and gate it. Split from [`run_exchange`] because the
/// stage publishes one last solve at the converged references, and the check
/// must bind to the bytes it publishes.
pub(super) fn cross_check_decomposition(
    outcome: ConjugateOutcome,
    final_robin_out_total_w: f64,
    final_robin_out_off_path_w: f64,
) -> Result<ConjugateOutcome, SolveRefusal> {
    let on_path = final_robin_out_total_w - final_robin_out_off_path_w;
    let solution = outcome.solution.with_decomposition_cross_check(on_path);
    let residual = solution
        .balance
        .decomposition_residual_w
        .unwrap_or(f64::NAN);
    if !(residual.is_finite() && residual.abs() <= outcome.balance_tolerance_w) {
        return Err(conduction_error(
            "cli-solve-conduction-airflow-decomposition",
            format!(
                "the path's summed solid heat rate differs from fs-conduction's own Robin accumulation by {residual:.4e} W against a {:.4e} W gate",
                outcome.balance_tolerance_w
            ),
            "a Robin face owned by no path region or counted twice; inspect the boundary partition",
        ));
    }
    Ok(ConjugateOutcome {
        solution,
        balance_tolerance_w: outcome.balance_tolerance_w,
    })
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

/// Render the receipt's `conjugate` object: the derivation disclosures per
/// segment, the converged march, the balance audit, and the authority and
/// no-claim statements.
pub(super) fn receipt_fragment(
    path: &ConjugatePath,
    outcome: &ConjugateOutcome,
) -> Result<String, SolveRefusal> {
    let mut segments = String::new();
    for (index, (derivation, state)) in path
        .segments
        .iter()
        .zip(outcome.solution.march.segments.iter())
        .enumerate()
    {
        if index > 0 {
            segments.push(',');
        }
        let balance = outcome.solution.balance.regions.get(index).ok_or_else(|| {
            conduction_error(
                "cli-solve-conduction-airflow-exchange",
                format!("the balance audit has no row for segment {index}"),
                "report the driver defect",
            )
        })?;
        segments.push_str(&format!(
            "{{\"target\":{},\"order\":{},\"card\":{},\"in_domain\":{},\"wetted_area_m2\":{},\"velocity_m_s\":{},\"reynolds\":{},\"length_over_hydraulic_diameter\":{},\"nusselt\":{},\"htc_w_m2_k\":{},\"air_in_k\":{},\"air_out_k\":{},\"reference_k\":{},\"ntu\":{},\"effectiveness\":{},\"solid_heat_rate_w\":{},\"air_heat_rate_w\":{},\"imbalance_w\":{}}}",
            json_string(&derivation.target),
            derivation.order,
            json_string(derivation.card.name()),
            derivation.in_domain,
            num(derivation.wetted_area_m2, "wetted_area_m2")?,
            num(derivation.velocity_m_s, "velocity_m_s")?,
            num(derivation.reynolds, "reynolds")?,
            num(derivation.length_over_hydraulic_diameter, "length_over_hydraulic_diameter")?,
            num(derivation.nusselt, "nusselt")?,
            num(derivation.htc_w_m2_k, "htc_w_m2_k")?,
            num(state.inlet_temperature_k, "air_in_k")?,
            num(state.outlet_temperature_k, "air_out_k")?,
            num(state.reference_temperature_k, "reference_k")?,
            num(state.ntu, "ntu")?,
            num(state.effectiveness, "effectiveness")?,
            num(balance.solid_heat_rate_w, "solid_heat_rate_w")?,
            num(balance.air_heat_rate_w, "air_heat_rate_w")?,
            num(balance.imbalance_w, "imbalance_w")?,
        ));
    }
    let audit = &outcome.solution.balance;
    Ok(format!(
        "{{\"branch\":{},\"path\":{},\"flow_m3_s\":{{\"lo\":{},\"mid\":{},\"hi\":{}}},\"air_density_kg_m3\":{},\"mass_flow_kg_s\":{},\"inlet_k\":{},\"outlet_k\":{},\"air_properties\":{{\"dynamic_viscosity_pa_s\":{},\"thermal_conductivity_w_m_k\":{},\"prandtl\":{},\"specific_heat_j_kg_k\":{},\"source\":{}}},\"segments\":[{}],\"iterations\":{},\"solid_total_w\":{},\"air_total_w\":{},\"interface_imbalance_w\":{},\"max_region_imbalance_w\":{},\"worst_recorded_imbalance_w\":{},\"decomposition_residual_w\":{},\"balance_tolerance_w\":{},\"authority\":{},\"no_claim\":{}}}",
        json_string(&path.branch),
        json_string(&path.path_name),
        num(path.flow_lo_m3_s, "flow_lo")?,
        num(path.flow_mid_m3_s, "flow_mid")?,
        num(path.flow_hi_m3_s, "flow_hi")?,
        num(path.air_density_kg_m3, "air_density")?,
        num(path.mass_flow_kg_s, "mass_flow")?,
        num(path.inlet_temperature_k, "inlet_k")?,
        num(outcome.solution.march.outlet_temperature_k, "outlet_k")?,
        num(AIR_DYNAMIC_VISCOSITY_PA_S, "mu")?,
        num(AIR_THERMAL_CONDUCTIVITY_W_M_K, "k")?,
        num(AIR_PRANDTL, "Pr")?,
        num(AIR_SPECIFIC_HEAT_J_KG_K, "cp")?,
        json_string(AIR_PROPERTY_SOURCE),
        segments,
        outcome.solution.iterations,
        num(audit.solid_total_w, "solid_total_w")?,
        num(audit.air_total_w, "air_total_w")?,
        num(audit.interface_imbalance_w, "interface_imbalance_w")?,
        num(audit.max_region_imbalance_w, "max_region_imbalance_w")?,
        num(
            outcome.solution.worst_recorded_imbalance_w,
            "worst_recorded_imbalance_w"
        )?,
        num(
            audit.decomposition_residual_w.unwrap_or(f64::NAN),
            "decomposition_residual_w"
        )?,
        num(outcome.balance_tolerance_w, "balance_tolerance_w")?,
        json_string(CONJUGATE_AUTHORITY),
        json_string(CONJUGATE_NO_CLAIM),
    ))
}
