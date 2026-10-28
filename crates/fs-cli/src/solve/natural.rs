//! Natural (buoyancy-driven) convection boundary law (fsim v10).
//!
//! A `(natural-convection ...)` row lowers to a Robin row whose coefficient is
//! DERIVED from a named fs-convection natural-convection card at the solved
//! mean wall-to-ambient difference. The coefficient depends on that
//! difference, so the stage iterates: solve with `h_k`, read the target's
//! area-mean wall temperature from the Robin decomposition, re-evaluate the
//! card, and stop when `h` moves by at most `TOLERANCE_REL`. The card is never
//! extrapolated: a point outside its Rayleigh domain refuses. The vertical
//! Churchill-Chu law admits heated and cooled walls: Ra uses |Tw-Ta|, while
//! the film temperature, Robin flux and retained temperature difference keep
//! their physical sign. Exact thermal equilibrium remains explicitly refused.

use std::collections::BTreeMap;

use fs_convection::{CorrelationId, CorrelationInputs, ThermalConductivity, evaluate};
use fs_project::{ConductionSetup, ThermalBoundaryCondition};
use fs_qty::Length;

use super::conjugate::{
    AIR_DYNAMIC_VISCOSITY_PA_S, AIR_PRANDTL, AIR_PROPERTY_SOURCE, AIR_THERMAL_CONDUCTIVITY_W_M_K,
};
use super::{SolveRefusal, canonical_f64, conduction_error};
use crate::import::json_string;

pub(super) mod adaptive;

/// Relative change in every coefficient that ends the fixed point.
pub(super) const TOLERANCE_REL: f64 = 1e-10;
/// Fixed-point iteration budget (the map contracts by about 1/4 per step).
pub(super) const MAX_ITERATIONS: usize = 80;
/// The native natural-convection law with an admitted classical derivative.
pub(super) const DIFFERENTIABLE_CARD: &str = "convection.churchill-chu-vertical-plate";
/// Initial wall-to-ambient difference the first coefficient is evaluated at.
const INITIAL_DELTA_T_K: f64 = 10.0;
/// Standard gravity, m/s².
const GRAVITY_M_S2: f64 = 9.806_65;
/// Dry-air specific gas constant, J/(kg·K).
const AIR_GAS_CONSTANT_J_KG_K: f64 = 287.05;

const AUTHORITY: &str = "card-derived natural-convection Robin coefficient at the solved area-mean wall-to-ambient difference, by fixed-point iteration to a declared relative tolerance; the card's Rayleigh domain gates every evaluation";
const NO_CLAIM: &str = "air transport properties are frozen at 300 K (density is ideal-gas at the film temperature and envelope pressure, beta = 1/T_film); the card's isothermal-plate idealization is applied to a non-isothermal wall through its area-mean temperature; no radiation-convection interaction beyond the declared laws, no enclosure or neighbouring-body effects, and no experimental validation are claimed";

/// One declared natural-convection row.
#[derive(Debug, Clone)]
pub(super) struct NaturalLaw {
    pub(super) target: String,
    length_m: f64,
    pub(super) ambient_k: f64,
    card: CorrelationId,
    /// Explicit engineering perturbation of h, not of the card's Ra or Nu.
    htc_multiplier: f64,
}

/// The rows of `setup` that are natural-convection laws, with their cards
/// resolved.
pub(super) fn natural_laws(setup: &ConductionSetup) -> Result<Vec<NaturalLaw>, SolveRefusal> {
    natural_laws_scaled(setup, 1.0)
}

/// Lower one complete discrepancy vertex. Its absolute coefficient multiplier
/// travels with each law, so initialization AND every wall-temperature update
/// use h_s(T) = s h_card(T). Nominal adjoints use natural_laws and retain s = 1.
pub(super) fn natural_laws_scaled(
    setup: &ConductionSetup,
    htc_multiplier: f64,
) -> Result<Vec<NaturalLaw>, SolveRefusal> {
    if !(htc_multiplier.is_finite() && htc_multiplier > 0.0) {
        return Err(conduction_error(
            "cli-solve-conduction-natural-scale",
            "a natural-convection coefficient multiplier must be finite and positive",
            "use a card allowance that preserves a positive coefficient",
        ));
    }
    let mut laws = Vec::new();
    for row in &setup.boundaries {
        let ThermalBoundaryCondition::NaturalConvection {
            characteristic_length,
            ambient_temperature,
            correlation,
        } = &row.condition
        else {
            continue;
        };
        let card = CorrelationId::ALL
            .iter()
            .copied()
            .find(|id| id.name() == correlation)
            .ok_or_else(|| {
                conduction_error(
                    "cli-solve-conduction-natural-correlation",
                    format!(
                        "natural convection on `{}` names `{correlation}`, which is not an fs-convection card",
                        row.target
                    ),
                    "name a natural-convection card such as `convection.churchill-chu-vertical-plate`",
                )
            })?;
        // A forced card cannot be evaluated from (Ra, Pr): refuse it up front.
        if evaluate(card, CorrelationInputs::natural(1.0e6, AIR_PRANDTL)).is_err() {
            return Err(conduction_error(
                "cli-solve-conduction-natural-correlation",
                format!(
                    "card `{correlation}` on `{}` is not a natural-convection card",
                    row.target
                ),
                "name a card evaluated from Rayleigh and Prandtl numbers",
            ));
        }
        laws.push(NaturalLaw {
            target: row.target.clone(),
            length_m: characteristic_length.value,
            ambient_k: ambient_temperature.value,
            card,
            htc_multiplier,
        });
    }
    Ok(laws)
}

/// The card evaluated at one wall-to-ambient difference.
#[derive(Debug, Clone, Copy)]
pub(super) struct Coefficient {
    pub(super) htc_w_m2_k: f64,
    rayleigh: f64,
    nusselt: f64,
}

/// Evaluate `law`'s card at `delta_t_k` (wall minus ambient) and envelope
/// pressure `pressure_pa`.
pub(super) fn coefficient(
    law: &NaturalLaw,
    delta_t_k: f64,
    pressure_pa: f64,
) -> Result<Coefficient, SolveRefusal> {
    if !(delta_t_k.is_finite() && delta_t_k != 0.0) {
        return Err(conduction_error(
            "cli-solve-conduction-natural-unheated",
            format!(
                "natural convection on `{}` needs a finite nonzero wall-to-ambient difference; found {delta_t_k} K",
                law.target
            ),
            "declare a heated or cooled vertical wall; equilibrium requires a separately admitted limiting law",
        ));
    }
    // Churchill-Chu is a vertical-plate law for heating OR cooling. Only Ra
    // uses the magnitude: the film temperature and the Robin heat flux must
    // retain the signed wall-minus-ambient difference. Other orientations
    // cannot be promoted to the cooled branch without their own admission.
    if delta_t_k < 0.0 && law.card != CorrelationId::ChurchillChuVerticalPlate {
        return Err(conduction_error(
            "cli-solve-conduction-natural-cooled-card",
            "the cooled-wall branch is admitted only for the Churchill-Chu vertical-plate card",
            "declare a supported vertical-plate law",
        ));
    }
    let film_k = law.ambient_k + 0.5 * delta_t_k;
    if ![law.length_m, law.ambient_k, law.ambient_k + delta_t_k, film_k, pressure_pa]
        .iter().all(|v| v.is_finite() && *v > 0.0)
    {
        return Err(conduction_error(
            "cli-solve-conduction-natural-input",
            "natural convection requires positive finite absolute temperatures, pressure and length",
            "inspect the declared law and the actual wall temperature",
        ));
    }
    let density = pressure_pa / (AIR_GAS_CONSTANT_J_KG_K * film_k);
    let kinematic = AIR_DYNAMIC_VISCOSITY_PA_S / density;
    let diffusivity = kinematic / AIR_PRANDTL;
    let length3 = law.length_m * law.length_m * law.length_m;
    let rayleigh = GRAVITY_M_S2 * (1.0 / film_k) * delta_t_k.abs() * length3 / (kinematic * diffusivity);
    let nusselt = evaluate(law.card, CorrelationInputs::natural(rayleigh, AIR_PRANDTL)).map_err(|error| {
        conduction_error(
            "cli-solve-conduction-natural-card",
            format!("card `{}` refused Ra {rayleigh:e} on `{}`: {error}", law.card.name(), law.target),
            "report the card defect or pick a card whose inputs cover this point",
        )
    })?;
    if !nusselt.evidence().model.in_domain {
        return Err(conduction_error(
            "cli-solve-conduction-natural-card-domain",
            format!(
                "card `{}` evaluated outside its declared domain on `{}` (Ra {rayleigh:e})",
                law.card.name(),
                law.target
            ),
            "the derivation never extrapolates a card; pick one whose domain covers the point",
        ));
    }
    let nu_value = nusselt.evidence().value;
    let htc = nusselt
        .heat_transfer_coefficient(
            ThermalConductivity::new(AIR_THERMAL_CONDUCTIVITY_W_M_K),
            Length::new(law.length_m),
        )
        .map_err(|error| {
            conduction_error(
                "cli-solve-conduction-natural-card",
                format!("card `{}` cannot lower Nu to a coefficient on `{}`: {error}", law.card.name(), law.target),
                "report the card defect",
            )
        })?
        .value
        .value() * law.htc_multiplier;
    if !(htc.is_finite() && htc > 0.0) {
        return Err(conduction_error(
            "cli-solve-conduction-natural-card",
            format!("card `{}` with multiplier {} produced coefficient {htc} on `{}`", law.card.name(), law.htc_multiplier, law.target),
            "report the card defect",
        ));
    }
    Ok(Coefficient { htc_w_m2_k: htc, rayleigh, nusselt: nu_value })
}

/// Complete wall and reservoir partials of the declared vertical-plate law.
/// Ra is proportional to |Tw-Ta|/Tfilm^3 at fixed pressure and transport data.
/// Its signed difference and BOTH film-density derivatives remain active.
pub(super) fn coefficient_partials(
    length: f64, wall: f64, ambient: f64, h: f64,
) -> Result<[f64; 2], SolveRefusal> {
    let bad = |message| conduction_error(
        "cli-solve-conduction-natural-derivative", message,
        "use an interior, non-equilibrium Churchill-Chu state with positive physical law values",
    );
    let finite = |value: f64| if value.is_finite() { Ok(value) } else {
        Err(bad("natural-convection derivative arithmetic is nonfinite"))
    };
    if ![length, wall, ambient, h].iter().all(|v| v.is_finite() && *v > 0.0)
        || wall == ambient
    {
        return Err(bad("natural adjoint needs a finite non-equilibrium wall and positive retained law values"));
    }
    let delta = finite(wall - ambient)?;
    let film = finite(ambient + 0.5 * delta)?;
    let nu = finite(h * length / AIR_THERMAL_CONDUCTIVITY_W_M_K)?;
    let exponent = finite((1.0 - 0.825 / fs_math::det::sqrt(nu)) / 3.0)?;
    if !(exponent > 0.0 && exponent < 1.0 / 3.0) {
        return Err(bad("retained natural coefficient is outside the Churchill-Chu differential branch"));
    }
    Ok([
        finite(h * exponent * (1.0 / delta - 1.5 / film))?,
        finite(h * exponent * (-1.0 / delta - 1.5 / film))?,
    ])
}

/// Re-evaluate the physical card and its two partials at one accepted wall.
/// Nominal trajectories differentiate the unscaled, in-domain card only.
/// Both controls require admitted adjacent representable temperatures: five
/// bounded card queries on success, with the other temperature held fixed.
pub(super) fn differentiated_coefficient(
    law: &NaturalLaw, wall_k: f64, pressure_pa: f64,
) -> Result<(Coefficient, [f64; 2]), SolveRefusal> {
    if law.card != CorrelationId::ChurchillChuVerticalPlate || law.htc_multiplier != 1.0 {
        return Err(conduction_error(
            "cli-solve-conduction-natural-derivative",
            "native natural adjoints require the nominal Churchill-Chu vertical-plate law",
            "declare the differentiated natural-convection card without a discrepancy multiplier",
        ));
    }
    let value = coefficient(law, wall_k - law.ambient_k, pressure_pa)?;
    // A nominal point alone cannot admit a classical two-sided partial at a
    // positivity, equilibrium or correlation-domain boundary. Reuse the
    // actual card owner for both neighborhoods, including its film density.
    for nearby_wall in [wall_k.next_down(), wall_k.next_up()] {
        coefficient(law, nearby_wall - law.ambient_k, pressure_pa)?;
    }
    for nearby_ambient in [law.ambient_k.next_down(), law.ambient_k.next_up()] {
        let mut nearby_law = law.clone();
        nearby_law.ambient_k = nearby_ambient;
        coefficient(&nearby_law, wall_k - nearby_ambient, pressure_pa)?;
    }
    let partials = coefficient_partials(law.length_m, wall_k, law.ambient_k, value.htc_w_m2_k)?;
    Ok((value, partials))
}

/// The first coefficient of every law, at the initial difference.
pub(super) fn initial_coefficients(
    laws: &[NaturalLaw],
    pressure_pa: f64,
) -> Result<BTreeMap<String, f64>, SolveRefusal> {
    laws.iter()
        .map(|law| Ok((law.target.clone(), coefficient(law, INITIAL_DELTA_T_K, pressure_pa)?.htc_w_m2_k)))
        .collect()
}

/// One converged law, for the receipt.
pub(super) struct Converged {
    pub(super) law: NaturalLaw,
    pub(super) coefficient: Coefficient,
    pub(super) mean_wall_k: f64,
    pub(super) heat_rate_w: f64,
}

/// Receipt fragment of a converged natural-convection fixed point.
pub(super) fn receipt_fragment(converged: &[Converged], iterations: usize) -> Result<String, SolveRefusal> {
    let num = |name: &str, value: f64| {
        canonical_f64(value).ok_or_else(|| {
            conduction_error(
                "cli-solve-conduction-nonfinite",
                format!("natural-convection field `{name}` is non-finite ({value})"),
                "report the solver defect",
            )
        })
    };
    let mut rows = Vec::with_capacity(converged.len());
    for row in converged {
        let multiplier = if row.law.htc_multiplier == 1.0 {
            String::new()
        } else {
            format!(",\"htc_multiplier\":{}", num("htc_multiplier", row.law.htc_multiplier)?)
        };
        rows.push(format!(
            "{{\"target\":{},\"card\":{},\"characteristic_length_m\":{},\"ambient_k\":{},\"htc_w_m2_k\":{},\"mean_wall_k\":{},\"delta_t_k\":{},\"rayleigh\":{},\"nusselt\":{},\"heat_rate_w\":{},\"in_domain\":true{multiplier}}}",
            json_string(&row.law.target),
            json_string(row.law.card.name()),
            num("length", row.law.length_m)?,
            num("ambient", row.law.ambient_k)?,
            num("htc", row.coefficient.htc_w_m2_k)?,
            num("mean_wall", row.mean_wall_k)?,
            num("delta_t", row.mean_wall_k - row.law.ambient_k)?,
            num("rayleigh", row.coefficient.rayleigh)?,
            num("nusselt", row.coefficient.nusselt)?,
            num("heat_rate", row.heat_rate_w)?,
        ));
    }
    Ok(format!(
        "{{\"laws\":[{}],\"iterations\":{iterations},\"tolerance_rel\":{},\"air_properties\":{{\"dynamic_viscosity_pa_s\":{},\"thermal_conductivity_w_m_k\":{},\"prandtl\":{},\"source\":{},\"density_basis\":\"ideal gas at film temperature and envelope pressure\",\"beta_basis\":\"1/T_film\"}},\"authority\":{},\"no_claim\":{}}}",
        rows.join(","),
        num("tolerance", TOLERANCE_REL)?,
        num("mu", AIR_DYNAMIC_VISCOSITY_PA_S)?,
        num("k", AIR_THERMAL_CONDUCTIVITY_W_M_K)?,
        num("pr", AIR_PRANDTL)?,
        json_string(AIR_PROPERTY_SOURCE),
        json_string(AUTHORITY),
        json_string(NO_CLAIM),
    ))
}

pub(super) fn has_mixed_boundaries(setup: &ConductionSetup) -> bool {
    setup.boundaries.iter().any(|row|
        matches!(row.condition, ThermalBoundaryCondition::NaturalConvection { .. }))
        && setup.boundaries.iter().any(|row|
            matches!(row.condition, ThermalBoundaryCondition::AirflowConvection { .. }))
}

pub(super) struct MixedSteady {
    pub exchange: super::conjugate::SteadyExchange,
    pub natural_receipt: String,
}

/// Picard on each natural coefficient around a COMPLETE shared-solid/air
/// exchange. Natural reservoirs never enter the air march. Acceptance checks
/// the actual wall law on the same final field as the air decomposition.
pub(super) fn solve_mixed_steady(
    cx: &fs_exec::Cx<'_>,
    laws: &[NaturalLaw],
    pressure_pa: f64,
    path: &super::conjugate::ConjugatePath,
    adaptive_goal_tolerance: Option<f64>,
    mut solve: impl FnMut(
        &BTreeMap<String, (f64, f64)>,
    ) -> Result<super::radiation::SolidSolution, SolveRefusal>,
) -> Result<MixedSteady, SolveRefusal> {
    let poll = || cx.checkpoint().map_err(|_| conduction_error(
        "cli-solve-conduction-mixed-cancelled",
        "the mixed natural/airflow steady solve was cancelled",
        "resume the accepted pipeline prefix",
    ));
    let finite = |value: f64| -> Result<f64, SolveRefusal> {
        if value.is_finite() { Ok(value) } else { Err(conduction_error(
            "cli-solve-conduction-nonfinite",
            "mixed natural/airflow feedback produced nonfinite arithmetic",
            "inspect the declared wall laws and applied heat rates",
        )) }
    };
    poll()?;
    if adaptive_goal_tolerance.is_some_and(|tolerance|
        !tolerance.is_finite() || tolerance <= 0.0 || tolerance >= 1.0)
    {
        return Err(conduction_error(
            "cli-solve-conduction-mixed-adaptive",
            "the mixed adaptive goal requires a finite residual tolerance in (0, 1)",
            "declare an admitted positive solver tolerance",
        ));
    }
    if laws.is_empty() || path.segments.is_empty()
        || laws.iter().any(|law| path.segments.iter().any(|segment| segment.target == law.target))
    {
        return Err(conduction_error(
            "cli-solve-conduction-mixed-targets",
            "mixed convection requires natural and airflow laws on distinct targets",
            "partition the exterior into separate natural-reservoir and airflow targets",
        ));
    }
    let adaptive = adaptive_goal_tolerance.is_some();
    let config = super::conjugate::exchange_config(adaptive);
    // The comparison independently rebuilds the complete physical residual.
    // Leave room below its gate for both the solid solve and changed wall laws.
    let coefficient_tolerance = adaptive_goal_tolerance.map_or(TOLERANCE_REL,
        |tolerance| TOLERANCE_REL.min(tolerance * 0.1));
    let mut htc = initial_coefficients(laws, pressure_pa)?;
    let mut total_air_iterations = 0usize;
    let mut total_solid_responses = 0usize;
    for iteration in 1..=MAX_ITERATIONS {
        poll()?;
        let natural_rows: BTreeMap<String, (f64, f64)> = laws.iter()
            .map(|law| (law.target.clone(), (htc[&law.target], law.ambient_k)))
            .collect();
        let mut exchange = super::conjugate::solve_steady_exchange(cx, path, adaptive, |air_rows| {
            let mut derived = natural_rows.clone();
            derived.extend(air_rows.iter().map(|(target, values)| (target.clone(), *values)));
            solve(&derived)
        })?;
        total_air_iterations += exchange.outcome.solution.iterations;
        total_solid_responses += exchange.solid_solves;
        let mut next = BTreeMap::new();
        let mut converged = Vec::with_capacity(laws.len());
        let mut worst_relative = 0.0f64;
        let mut worst_heat_w = 0.0f64;
        let mut worst_heat_ratio = 0.0f64;
        for law in laws {
            poll()?;
            let flux = exchange.solid.convective_robin_fluxes().iter()
                .find(|flux| flux.region == law.target)
                .ok_or_else(|| conduction_error(
                    "cli-solve-conduction-natural-region",
                    format!("the conduction report has no Robin decomposition for `{}`", law.target),
                    "report the driver defect; every natural-convection target lowers to a Robin region",
                ))?;
            let point = coefficient(law, flux.mean_wall_temperature_k - law.ambient_k, pressure_pa)?;
            let applied = htc[&law.target];
            worst_relative = worst_relative.max(finite((point.htc_w_m2_k - applied).abs() / applied)?);
            let physical_heat_w = finite(
                point.htc_w_m2_k * flux.area_m2 * (flux.mean_wall_temperature_k - law.ambient_k),
            )?;
            let mismatch_w = finite((physical_heat_w - flux.heat_rate_w).abs())?;
            // Each natural reservoir has its own scale: a large forced-air
            // load cannot hide a failed passive-boundary heat gate.
            let tolerance_w = config.balance_tolerance_w.max(
                config.balance_relative_tolerance * physical_heat_w.abs().max(flux.heat_rate_w.abs()),
            );
            worst_heat_w = worst_heat_w.max(mismatch_w);
            worst_heat_ratio = worst_heat_ratio.max(finite(mismatch_w / tolerance_w)?);
            next.insert(law.target.clone(), point.htc_w_m2_k);
            converged.push(Converged {
                law: law.clone(), coefficient: point,
                mean_wall_k: flux.mean_wall_temperature_k, heat_rate_w: flux.heat_rate_w,
            });
        }
        if worst_relative <= coefficient_tolerance && worst_heat_ratio <= 1.0 {
            // Retain the rows actually applied to the published field, rather
            // than silently replacing them with the next coefficient iterate.
            exchange.derived.extend(natural_rows);
            let mut natural_receipt = receipt_fragment(&converged, iteration)?;
            natural_receipt.pop();
            natural_receipt.push_str(&format!(
                ",\"coupling\":{{\"method\":\"natural-picard-around-complete-air-exchange\",\"max_outer_iterations\":{MAX_ITERATIONS},\"total_air_iterations\":{total_air_iterations},\"total_solid_responses\":{total_solid_responses},\"coefficient_tolerance_rel\":{},\"max_relative_coefficient_change\":{},\"max_natural_heat_mismatch_w\":{},\"max_natural_heat_gate_ratio\":{},\"scope\":{}}}}}",
                canonical_f64(coefficient_tolerance).expect("finite positive coefficient tolerance"),
                canonical_f64(worst_relative).expect("finite relative change"),
                canonical_f64(worst_heat_w).expect("finite heat mismatch"),
                canonical_f64(worst_heat_ratio).expect("finite heat gate ratio"),
                json_string("Steady natural reservoirs and forced-air surfaces exchange through the complete shared solid. Each outer iteration closes every air branch; natural heat and radiation remain outside air transport. Counts cover solid responses; inner material/radiation work retains its existing separate caps. Adaptive goals, when requested, separately check the complete physical residual and transpose; this heat balance alone is not a spatial error estimate."),
            ));
            poll()?;
            return Ok(MixedSteady { exchange, natural_receipt });
        }
        if iteration == MAX_ITERATIONS {
            return Err(conduction_error(
                "cli-solve-conduction-mixed-unconverged",
                format!("mixed steady feedback exhausted {MAX_ITERATIONS} outer iterations; relative h change {worst_relative:e}, natural heat mismatch {worst_heat_w:e} W, heat gate ratio {worst_heat_ratio:e}"),
                "inspect the separate reservoir temperatures, convection cards and solver tolerance; no partial field was published",
            ));
        }
        htc = next;
    }
    unreachable!("the nonzero bounded loop returns a solution or refusal")
}

/// Admit only natural cards whose full coefficient tangent is implemented.
/// Base and ladder keep their original card surface. Mixed goals retain
/// natural coefficients and the full air-reference feedback on distinct traces.
pub(super) fn admit_fidelity(
    spec: &fs_project::ProjectSpec,
    setup: &ConductionSetup,
) -> Result<(), SolveRefusal> {
    if !spec.solver.as_ref().is_some_and(|solver|
        solver.fidelity == super::SOLVER_FIDELITY_ADAPTIVE)
    {
        return Ok(());
    }
    let laws = natural_laws(setup)?;
    if laws.len() > 64 || laws.iter().any(|law|
        law.card != CorrelationId::ChurchillChuVerticalPlate)
    {
        return Err(conduction_error(
            "cli-solve-conduction-natural-adaptive",
            "adaptive natural goals require at most 64 differentiated Churchill-Chu vertical-plate laws",
            "use supported vertical-plate laws or select base or ladder fidelity",
        ));
    }
    Ok(())
}

#[cfg(test)]
mod discrepancy_tests;

#[cfg(test)]
mod tests {
    use super::*;

    fn law(ambient_k: f64) -> NaturalLaw {
        NaturalLaw {
            target: "wall".into(), length_m: 0.06, ambient_k,
            card: CorrelationId::ChurchillChuVerticalPlate,
            htc_multiplier: 1.0,
        }
    }

    #[test]
    fn reversing_wall_and_ambient_preserves_the_card_but_reverses_heat_flow() {
        // Both states have exactly the same film temperature and |delta|.
        let hot = coefficient(&law(295.0), 10.0, 101325.0).unwrap();
        let cold = coefficient(&law(305.0), -10.0, 101325.0).unwrap();
        assert_eq!(hot.htc_w_m2_k.to_bits(), cold.htc_w_m2_k.to_bits());
        assert_eq!(hot.rayleigh.to_bits(), cold.rayleigh.to_bits());
        assert_eq!(hot.nusselt.to_bits(), cold.nusselt.to_bits());
        assert!(cold.htc_w_m2_k > 0.0);
        assert_eq!(hot.htc_w_m2_k * 10.0, -(cold.htc_w_m2_k * -10.0));
        // At a fixed ambient, density changes with the SIGNED film state.
        let same_ambient_hot = coefficient(&law(305.0), 10.0, 101325.0).unwrap();
        assert!(cold.htc_w_m2_k > same_ambient_hot.htc_w_m2_k);
    }

    #[test]
    fn differentiated_card_checks_both_temperature_controls_and_film_density() {
        for ambient in [290.0, 320.0] {
            for delta in [-80.0, -25.0, -1.0, 1.0, 25.0, 80.0] {
                for pressure in [80_000.0, 101_325.0] {
                    let input = law(ambient);
                    let wall = ambient + delta;
                    let (point, [dw, da]) =
                        differentiated_coefficient(&input, wall, pressure).unwrap();
                    let evaluate = |tw: f64, ta: f64| {
                        let mut nearby = input.clone();
                        nearby.ambient_k = ta;
                        coefficient(&nearby, tw - ta, pressure).unwrap().htc_w_m2_k
                    };
                    assert_eq!(point.htc_w_m2_k.to_bits(), evaluate(wall, ambient).to_bits());
                    let step = 1e-3;
                    let expected_w =
                        (evaluate(wall + step, ambient) - evaluate(wall - step, ambient))
                            / (2.0 * step);
                    let expected_a =
                        (evaluate(wall, ambient + step) - evaluate(wall, ambient - step))
                            / (2.0 * step);
                    for (actual, expected) in [(dw, expected_w), (da, expected_a)] {
                        assert!((actual - expected).abs() < 2e-6 * expected.abs().max(1e-3),
                            "{actual:e} vs {expected:e}; wall {wall}, ambient {ambient}, pressure {pressure}");
                    }
                    // Freezing film density would incorrectly make da = -dw.
                    assert!(dw + da < 0.0, "both partials must retain film-density feedback");
                }
            }
        }
        assert!(differentiated_coefficient(&law(300.0), 300.0, 101_325.0).is_err());
        let mut unsupported = law(300.0);
        unsupported.card = CorrelationId::Gnielinski;
        assert!(differentiated_coefficient(&unsupported, 310.0, 101_325.0).is_err());
        unsupported = law(300.0);
        unsupported.htc_multiplier = 1.1;
        assert!(differentiated_coefficient(&unsupported, 310.0, 101_325.0).is_err());
        let mut outside = law(300.0);
        outside.length_m = 1.0e9;
        assert!(differentiated_coefficient(&outside, 290.0, 101_325.0).is_err());

        // These nominal card points are valid, but one adjacent wall or
        // reservoir is not. Each claimed two-sided control must refuse.
        let mut wall_edge = law(300.0);
        wall_edge.length_m = 60.0;
        let wall = wall_edge.ambient_k.next_up();
        assert!(coefficient(&wall_edge, wall - wall_edge.ambient_k, 101_325.0).is_ok());
        assert!(differentiated_coefficient(&wall_edge, wall, 101_325.0).is_err());
        let ambient_edge = law(f64::from_bits(1));
        assert!(coefficient(&ambient_edge, 300.0 - ambient_edge.ambient_k, 101_325.0).is_ok());
        assert!(differentiated_coefficient(&ambient_edge, 300.0, 101_325.0).is_err());
    }

    #[test]
    fn cooling_does_not_admit_equilibrium_nonphysical_inputs_or_other_cards() {
        for delta in [0.0, -0.0, f64::NAN, f64::INFINITY, -305.0, -306.0] {
            assert!(coefficient(&law(305.0), delta, 101325.0).is_err());
        }
        for pressure in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(coefficient(&law(305.0), -10.0, pressure).is_err());
        }
        let mut other = law(305.0);
        other.card = CorrelationId::Gnielinski;
        assert!(coefficient(&other, -10.0, 101325.0).is_err());
        other = law(305.0);
        other.length_m = 1.0e9;
        assert!(coefficient(&other, -10.0, 101325.0).is_err(),
            "a cold wall cannot bypass the actual Rayleigh domain");
    }
}
