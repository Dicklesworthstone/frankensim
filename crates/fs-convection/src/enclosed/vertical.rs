//! Wright's sealed vertical cavity relation, ASHRAE Transactions 102(1),
//! 940-946 (1996), Eqs. (9a-c) and (7), combined as specified under
//! "Low Aspect Ratio" on p. 945. Ra and Nu use the full normal gas gap.

use super::{
    Ambition, CorrelationCard, CorrelationId, DiscrepancyBasis, ModelCard, SourceProvenance,
    ValidityDomain,
};
use fs_math::det;

const LOW_COEFFICIENT: f64 = 1.75967e-10;
const LOW_EXPONENT: f64 = 2.2984755;
const HEIGHT_EXPONENT: f64 = 0.272;

fn components(rayleigh: f64, height_over_gap: f64) -> (f64, f64) {
    let tall = if rayleigh <= 1.0e4 {
        1.0 + LOW_COEFFICIENT * det::pow(rayleigh, LOW_EXPONENT)
    } else if rayleigh <= 5.0e4 {
        0.028154 * det::pow(rayleigh, 0.4134)
    } else {
        0.0673838 * det::pow(rayleigh, 1.0 / 3.0)
    };
    let finite_height = 0.242 * det::pow(rayleigh / height_over_gap, HEIGHT_EXPONENT);
    (tall, finite_height)
}

pub(crate) fn nusselt(rayleigh: f64, height_over_gap: f64) -> f64 {
    let (tall, finite_height) = components(rayleigh, height_over_gap);
    tall.max(finite_height)
}

/// Absolute partials, without reconstructing small derivatives from Nu-1.
/// The parent admits domains and excludes derivatives across their endpoints.
pub(crate) fn partial(rayleigh: f64, height_over_gap: f64, axis: &str) -> Option<f64> {
    if axis == "Re" {
        return None;
    }
    if axis != "Ra" && axis != "aspect_ratio" {
        return Some(0.0);
    }
    let (tall, finite_height) = components(rayleigh, height_over_gap);
    if tall == finite_height {
        return None;
    }
    if axis == "aspect_ratio" {
        return Some(if finite_height > tall {
            -HEIGHT_EXPONENT * finite_height / height_over_gap
        } else {
            0.0
        });
    }
    if rayleigh == 1.0e4 || rayleigh == 5.0e4 {
        // Preserve source branch boundaries and their rounded coefficients;
        // do not invent a smooth join or an unsourced transition blend.
        return None;
    }
    Some(if finite_height > tall {
        HEIGHT_EXPONENT * finite_height / rayleigh
    } else if rayleigh < 1.0e4 {
        LOW_COEFFICIENT * LOW_EXPONENT * det::pow(rayleigh, LOW_EXPONENT - 1.0)
    } else if rayleigh < 5.0e4 {
        0.4134 * tall / rayleigh
    } else {
        tall / (3.0 * rayleigh)
    })
}

pub(crate) fn card() -> CorrelationCard {
    let id = CorrelationId::WrightVerticalAirLayer;
    CorrelationCard {
        id,
        source: SourceProvenance {
            citation: "J. L. Wright, A Correlation to Quantify Convective Heat Transfer Between Vertical Window Glazings, ASHRAE Transactions 102(1), 1996, 940-946, Eqs. (9a-c) and (7), Low Aspect Ratio section",
            identifier: "ASHRAE AT-96-16-2; https://dspacemainprd01.lib.uwaterloo.ca/server/api/core/bitstreams/2c7962f1-23f2-45ce-855c-5ca72c711eda/content",
        },
        model: ModelCard::new(
            id.name(), "1.0.0", Ambition::Solid,
            [
                "sealed vertical rectangular air cavity between parallel isothermal walls, without imposed throughflow or wall shear",
                "Ra is nonnegative, based on full normal gap and absolute wall temperature difference; properties use mean wall temperature",
                "aspect_ratio is whole vertical height/gap, not transverse width or the smaller lateral span",
                "Ra in [0,1e6] and height/gap in [5,110] conservatively retain source measurement limits; Ra=0 is total conduction Nu=1",
                "Pr in [0.70,0.72] is an implementation admission around source air near 0.71, not a published validation interval or proof of fluid identity",
                "quasi-steady Boussinesq gas transfer; finite wall storage belongs to the thermal owner",
                "Nu includes gas conduction once; independent radiation and solid sidewall thermal bridges are separate",
                "the finite-height term is the source's specified maximum, not an additive heat path",
                "the discrepancy is an engineering allowance, not a measured confidence interval",
            ].into_iter().map(str::to_string).collect(),
            ValidityDomain::unconstrained()
                .with("Ra", 0.0, 1.0e6)
                .with("Pr", 0.70, 0.72)
                .with("aspect_ratio", 5.0, 110.0),
            [
                "tilted gaps, open ventilation, phase change, curved or non-isothermal walls",
                "non-air media inferred from Pr alone, rarefaction, significant gas storage or non-Boussinesq property variation",
                "unresolved transverse confinement or edge thermal conditions inconsistent with the cavity model",
                "smooth derivatives at Ra=1e4, Ra=5e4 or equal competing Nusselt branches",
            ].into_iter().map(str::to_string).collect(),
            0.20,
        ),
        discrepancy_basis: DiscrepancyBasis::EngineeringAllowance,
    }
}
