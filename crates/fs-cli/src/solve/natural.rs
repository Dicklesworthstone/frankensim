//! Natural (buoyancy-driven) convection boundary law (fsim v10).
//!
//! A `(natural-convection ...)` row lowers to a Robin row whose coefficient is
//! DERIVED from a named fs-convection natural-convection card at the solved
//! mean wall-to-ambient difference. The coefficient depends on that
//! difference, so the stage iterates: solve with `h_k`, read the target's
//! area-mean wall temperature from the Robin decomposition, re-evaluate the
//! card, and stop when `h` moves by at most `TOLERANCE_REL`. The card is never
//! extrapolated: a point outside its Rayleigh domain refuses.

use std::collections::BTreeMap;

use fs_convection::{CorrelationId, CorrelationInputs, ThermalConductivity, evaluate};
use fs_project::{ConductionSetup, ThermalBoundaryCondition};
use fs_qty::Length;

use super::conjugate::{
    AIR_DYNAMIC_VISCOSITY_PA_S, AIR_PRANDTL, AIR_PROPERTY_SOURCE, AIR_THERMAL_CONDUCTIVITY_W_M_K,
};
use super::{SolveRefusal, canonical_f64, conduction_error};
use crate::import::json_string;

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
}

/// The rows of `setup` that are natural-convection laws, with their cards
/// resolved.
pub(super) fn natural_laws(setup: &ConductionSetup) -> Result<Vec<NaturalLaw>, SolveRefusal> {
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
    if !(delta_t_k.is_finite() && delta_t_k > 0.0) {
        return Err(conduction_error(
            "cli-solve-conduction-natural-unheated",
            format!(
                "natural convection on `{}` needs a wall warmer than the {} K ambient; the solved mean wall-to-ambient difference is {delta_t_k} K",
                law.target, law.ambient_k
            ),
            "natural convection is driven by a heated wall: declare power or a warmer boundary, or use a fixed-coefficient convection law",
        ));
    }
    let film_k = law.ambient_k + 0.5 * delta_t_k;
    let density = pressure_pa / (AIR_GAS_CONSTANT_J_KG_K * film_k);
    let kinematic = AIR_DYNAMIC_VISCOSITY_PA_S / density;
    let diffusivity = kinematic / AIR_PRANDTL;
    let length3 = law.length_m * law.length_m * law.length_m;
    let rayleigh = GRAVITY_M_S2 * (1.0 / film_k) * delta_t_k * length3 / (kinematic * diffusivity);
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
        .value();
    if !(htc.is_finite() && htc > 0.0) {
        return Err(conduction_error(
            "cli-solve-conduction-natural-card",
            format!("card `{}` produced coefficient {htc} on `{}`", law.card.name(), law.target),
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
        rows.push(format!(
            "{{\"target\":{},\"card\":{},\"characteristic_length_m\":{},\"ambient_k\":{},\"htc_w_m2_k\":{},\"mean_wall_k\":{},\"delta_t_k\":{},\"rayleigh\":{},\"nusselt\":{},\"heat_rate_w\":{},\"in_domain\":true}}",
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
