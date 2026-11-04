//! Sealed air-layer correlations: horizontal/inclined Hollands and vertical Wright.
//!
//! Source: Hollands, Raithby and Unny, report COO-2597-2 (1976), §4,
//! Eq. (4-1), printed p.15; doi:10.2172/7125623. The authors reproduce
//! Hollands et al., J. Heat Transfer 98(2), 189–193 (1976),
//! doi:10.1115/1.3450517. At zero inclination the two positive-part
//! expressions reduce to the branches below. The pure conduction branch
//! also applies to stable horizontal stratification; this is an analytic
//! limit, not an extrapolation of a buoyant power law to negative Ra.
//! The inclined implementation cites the EnergyPlus Engineering Reference's
//! explicit reproduction of the equation; see that card's distinct source.

use super::{
    Ambition, CorrelationCard, CorrelationId, DiscrepancyBasis, ModelCard, SourceProvenance,
    ValidityDomain,
};

mod transfer;
pub(super) mod vertical;
pub(super) use transfer::{
    horizontal_nusselt as nusselt, horizontal_partial as partial, inclined_nusselt, inclined_partial,
};

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

pub(super) fn inclined_card() -> CorrelationCard {
    let id = CorrelationId::HollandsInclinedAirLayer;
    CorrelationCard {
        id,
        source: SourceProvenance {
            citation: "Hollands, Unny, Raithby and Konicek, Free Convective Heat Transfer Across Inclined Air Layers, J. Heat Transfer 98(2), 1976, 189-193; equation reproduced in EnergyPlus 8.3 Engineering Reference, Solar Collectors, Integral-collector-storage, Heat Transfer Coefficients",
            identifier: "doi:10.1115/1.3450517; https://bigladdersoftware.com/epx/docs/8-3/engineering-reference/solar-collectors.html#heat-transfer-coefficients",
        },
        model: ModelCard::new(
            id.name(), "1.0.0", Ambition::Solid,
            [
                "parallel rigid isothermal walls containing continuum air heated from below, without imposed flow",
                "full-gravity Ra = g beta (T_lower-T_upper) gap^3/(nu alpha); angle measured from horizontal and represented explicitly in degrees",
                "properties at the mean of the actual two wall temperatures; length is the full normal gas gap",
                "the source's sin(1.8 theta)^1.6 correction is retained; replacing only Ra with Ra cos(theta) is not this model",
                "Nu includes gas conduction and buoyant convection together; radiation is separate and no second gas Fourier conductance may be added",
                "0..=60 degrees is a conservative implementation subset of the original reported inclination measurements, not an additional experimental validation claim",
                "Pr interval [0.70,0.72] is a conservative implementation admission around air Pr near 0.71, not a source-published validation range; Pr alone cannot identify air",
                "gap/minimum whole-plate span <0.1 is an explicit thin-layer admission on both spans, not a mesh-cell dimension",
                "Ra<=1e5 retains the original published experimental ceiling; subcritical projected Ra and zero driving use the total conduction limit Nu=1",
                "quasi-steady transfer neglects gas thermal storage and requires bounded property/density variation; wall storage remains with the thermal owner",
            ].into_iter().map(str::to_string).collect(),
            ValidityDomain::unconstrained()
                .with("Ra", 0.0, 1.0e5)
                .with("Pr", 0.70, 0.72)
                .with("aspect_ratio", 0.0, 0.1_f64.next_down())
                .with("inclination_degrees", 0.0, 60.0),
            [
                "inclined heating from above, non-air media inferred from Pr alone, or inclination above 60 degrees",
                "non-isothermal/curved walls, small lateral spans, sidewall thermal bridges, throughflow, wall shear or phase change",
                "non-Boussinesq property variation, rarefied gas or resolved transient gas storage",
                "smooth derivatives at projected-Rayleigh branch switches Ra cos(theta)=1708 or 5830",
            ].into_iter().map(str::to_string).collect(),
            0.20,
        ),
        discrepancy_basis: DiscrepancyBasis::EngineeringAllowance,
    }
}
