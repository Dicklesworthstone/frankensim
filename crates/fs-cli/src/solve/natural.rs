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

/// Admit only natural cards whose full coefficient tangent is implemented.
/// Base and ladder keep their original card surface. Radiation retains its
/// separate capability gate; no combined-law derivative is inferred here.
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
