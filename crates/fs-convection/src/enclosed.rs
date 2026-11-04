//! Horizontal limit of the original Hollands air-layer correlation.
//!
//! Source: Hollands, Raithby and Unny, report COO-2597-2 (1976), §4,
//! Eq. (4-1), printed p.15; doi:10.2172/7125623. The authors reproduce
//! Hollands et al., J. Heat Transfer 98(2), 189–193 (1976),
//! doi:10.1115/1.3450517. At zero inclination the two positive-part
//! expressions reduce to the branches below. The pure conduction branch
//! also applies to stable horizontal stratification; this is an analytic
//! limit, not an extrapolation of a buoyant power law to negative Ra.

use super::{
    Ambition, CorrelationCard, CorrelationId, DiscrepancyBasis, ModelCard, SourceProvenance,
    ValidityDomain, det,
};

pub(super) fn nusselt(signed_rayleigh: f64) -> f64 {
    if signed_rayleigh <= 1708.0 {
        return 1.0;
    }
    // Subtraction before division retains the small onset excess. This is
    // algebraically the source's 1 - 1708/Ra, not a smoothed transition.
    let onset = 1.44 * ((signed_rayleigh - 1708.0) / signed_rayleigh);
    let high_rayleigh = if signed_rayleigh > 5830.0 {
        det::pow(signed_rayleigh / 5830.0, 1.0 / 3.0) - 1.0
    } else {
        0.0
    };
    1.0 + onset + high_rayleigh
}

/// Absolute partials avoid an elasticity divided by signed/zero Rayleigh.
/// The source's positive parts are continuous but kink at 1708 and 5830.
pub(super) fn partial(signed_rayleigh: f64, axis: &str) -> Option<f64> {
    if axis == "Re" {
        return None;
    }
    if axis != "Ra" {
        return Some(0.0);
    }
    if signed_rayleigh == 1708.0 || signed_rayleigh == 5830.0 {
        return None;
    }
    if signed_rayleigh < 1708.0 {
        return Some(0.0);
    }
    let onset = (1.44 * 1708.0 / signed_rayleigh) / signed_rayleigh;
    let high_rayleigh = if signed_rayleigh > 5830.0 {
        det::pow(signed_rayleigh / 5830.0, 1.0 / 3.0) / (3.0 * signed_rayleigh)
    } else {
        0.0
    };
    Some(onset + high_rayleigh)
}

pub(super) fn card() -> CorrelationCard {
    let id = CorrelationId::HollandsHorizontalAirLayer;
    CorrelationCard {
        id,
        source: SourceProvenance {
            citation: "K. G. T. Hollands, G. D. Raithby and T. E. Unny, Studies on Methods of Reducing Heat Losses from Flat Plate Solar Collectors, report COO-2597-2, 1976, Sec. 4 Eq. (4-1), horizontal limit, pp. 15-17 and nomenclature p. 30; reproduces Hollands et al., J. Heat Transfer 98(2), 1976, 189-193",
            identifier: "doi:10.2172/7125623; doi:10.1115/1.3450517",
        },
        model: ModelCard::new(
            id.name(),
            "1.0.0",
            Ambition::Solid,
            [
                "parallel rigid horizontal isothermal walls containing a continuum air layer, with no imposed flow",
                "properties at the mean of the actual lower and upper wall temperatures; characteristic length is the full normal gas gap",
                "signed Ra = g beta (T_lower-T_upper) gap^3/(nu alpha); stable stratification and zero driving use the analytic total conduction limit Nu=1",
                "Nu includes gas conduction and buoyant convection together; radiation is separate and no second Fourier gas conductance may be added",
                "Pr interval [0.70,0.72] is a conservative implementation admission around source air Pr near 0.71, not a source-published validation range",
                "gap/minimum whole-plate span <0.1 applies the report's thin-layer assumption conservatively to both spans; mesh faces do not restart the layer geometry",
                "unstable Rayleigh numbers are bounded by the cited experiments at Ra<=1e5; finite negative Ra describes stable conduction rather than that experimental branch",
                "quasi-steady gas transfer with negligible layer heat storage and bounded property/density variation; wall energy storage remains with the thermal owner",
            ].into_iter().map(str::to_string).collect(),
            ValidityDomain::unconstrained()
                .with("Ra", -f64::MAX, 1.0e5)
                .with("Pr", 0.70, 0.72)
                .with("aspect_ratio", 0.0, 0.1_f64.next_down()),
            [
                "inclined or vertical cavities, curved or non-isothermal walls, small lateral spans or sidewall thermal bridges",
                "forced throughflow, moving-wall shear, evaporation, radiation included as gas heat, or rarefied gas",
                "non-air fluids inferred only from Pr, non-Boussinesq property variation, or transient gas storage comparable to wall storage",
                "smooth derivatives across the physical branch switches Ra=1708 and Ra=5830",
            ].into_iter().map(str::to_string).collect(),
            0.20,
        ),
        discrepancy_basis: DiscrepancyBasis::EngineeringAllowance,
    }
}
