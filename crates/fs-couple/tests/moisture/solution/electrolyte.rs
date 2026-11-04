//! Actual source-backed salt chemistry on the native finite-water endpoint.
use super::*;
use fs_material::osmosis::solubility::electrolyte::{
    ElectrolytePhaseEquilibrium, ElectrolyteSolubilitySolution, SolventActivityBasis,
};
use fs_thermochem::stoichiometry::thermodynamics::acid_base::{
    SpectatorIon,
    precipitation::{PrecipitationOptions, PrecipitationStatus},
};

// Share the original atom/charge/source-qualified manufactured chemistry;
// these values are analytical references, not a measured salt calibration.
#[path = "../../../../fs-material/src/osmosis/solubility/electrolyte/tests/fixture.rs"]
mod fixture;

fn basis() -> SolventActivityBasis {
    SolventActivityBasis::new(WATER_VOLUME, AmountConcentration::new(1000.0)).unwrap()
}
fn chemistry_options() -> PrecipitationOptions {
    let mut options = PrecipitationOptions::default();
    options.aqueous.ph_tolerance = fs_qty::Dimensionless::new(1e-13);
    options.aqueous.absolute_charge_tolerance = Amount::new(1e-18);
    options.aqueous.relative_charge_tolerance = fs_qty::Dimensionless::new(1e-12);
    options.aqueous.max_iterations = 256;
    options.log_activity_tolerance = fs_qty::Dimensionless::new(1e-11);
    options.log_saturation_tolerance = fs_qty::Dimensionless::new(1e-11);
    options
}
fn electrolyte(volume: f64, cation: f64, family: f64) -> ElectrolyteSolubilitySolution {
    ElectrolyteSolubilitySolution::from_salt(
        fixture::salt(293.15, 4e-6),
        basis(),
        Volume::new(volume),
        Amount::new(cation),
        vec![(fixture::FAMILY.into(), Amount::new(family))],
        vec![],
        chemistry_options(),
        || false,
    )
    .unwrap()
}

#[test]
fn native_water_exchange_creates_and_exhausts_non_unit_ionic_salt() {
    with_cx(|cx| {
        let final_volume = 0.001;
        let slope = 0.2;
        let diffusion = 0.01;
        let dt = 2.0;
        let alpha = 4.0 * dt * diffusion;
        // At pH7, half of the selected acid family is the terminal anion.
        // Ca(CN)2 saturation is aCa*aCN^2=4e-6. Each case has an independent
        // pH7 endpoint; the second is slightly undersaturated after dilution.
        for (initial_m, cation_total, family_total, final_cation, final_family, final_solid) in [
            (0.12, 0.01001, 0.04002, 0.01, 0.04, 1e-5),
            (0.25, 0.00999, 0.03996, 0.00999, 0.03996, 0.0),
        ] {
            let particles = final_cation + final_family + 2e-7;
            let activity = (-WATER_VOLUME * particles / final_volume).exp();
            let initial = slab(card([diffusion; 2], 0.3), 1.0, vec![initial_m], cx);
            let final_m = (initial_m + alpha * slope * activity) / (1.0 + alpha);
            let water = initial.dry_mass_kg() * (final_m - initial_m);
            let initial_volume = final_volume + water / WATER_DENSITY;
            let bath = FiniteMoistureSolution::with_electrolyte_solubility(
                electrolyte(initial_volume, cation_total, family_total),
                WATER_MOLAR_MASS,
            )
            .unwrap();
            let original = bath.electrolyte_solubility().unwrap();
            let start_solid = original.equilibrium().solid_amount().value();
            let faces = [
                Some(atmosphere::sorption(slope)),
                Some(atmosphere::sorption(slope)),
            ];
            let mut solve = controls(initial_volume);
            solve.mass_tolerance_kg = 1e-11;
            let accepted = bath
                .advance(&initial, dt, &faces, config(), solve, cx)
                .unwrap();
            let endpoint = accepted.reservoir().electrolyte_solubility().unwrap();
            let ElectrolytePhaseEquilibrium::Salt(phase) = endpoint.equilibrium() else {
                panic!("original ionic salt owner must survive moisture exchange")
            };
            near(endpoint.solution().solvent_volume_m3(), final_volume, 2e-14);
            near(accepted.moisture().moisture()[0], final_m, 2e-11);
            near(accepted.step().water_activity[1], activity, 2e-11);
            near(accepted.step().water_to_solid_kg, water, 1e-11);
            near(accepted.step().water_balance_residual_kg, 0.0, 1e-11);
            near(phase.aqueous.ph.value(), 7.0, 2e-8);
            near(phase.dissolved_cation.value(), final_cation, 2e-10);
            near(phase.dissolved_family.value(), final_family, 2e-10);
            near(phase.solid_amount.value(), final_solid, 2e-10);
            near(phase.cation_balance_residual.value(), 0.0, 1e-15);
            near(phase.family_balance_residual.value(), 0.0, 1e-15);
            assert_eq!(endpoint.component_total(), original.component_total());
            assert_eq!(endpoint.family_totals(), original.family_totals());
            assert_eq!(endpoint.spectators(), original.spectators());
            assert_eq!(endpoint.species(), original.species());
            assert!(accepted.reservoir().solubility().is_none());
            if water > 0.0 {
                assert_eq!(start_solid, 0.0);
                assert_eq!(phase.status, PrecipitationStatus::SolidPresent);
                assert!(phase.solid_amount.value() > 0.0);
            } else {
                assert!(start_solid > 0.0);
                assert_eq!(phase.status, PrecipitationStatus::SolidAbsent);
                assert_eq!(phase.solid_amount.value(), 0.0);
            }
            for source in accepted.step().transport.surface_sorption.iter().flatten() {
                near(source.query().axes()["relative_humidity"], activity, 2e-11);
            }
            assert_eq!(initial.moisture(), &[initial_m]);
            assert_eq!(original.solution().solvent_volume_m3(), initial_volume);
            assert_eq!(original.equilibrium().solid_amount().value(), start_solid);
        }
    });
}

#[test]
fn finite_water_owner_rejects_unaccounted_reaction_water() {
    let hydrolysis = ElectrolyteSolubilitySolution::from_hydrolysis(
        fixture::hydrolysis(293.15, 1e-9),
        basis(),
        Volume::new(0.001),
        Amount::new(1.1e-4),
        vec![],
        vec![SpectatorIon {
            species: fixture::id("Cl[-](aq)"),
            charge: -1,
            amount: Amount::new(2e-5),
        }],
        chemistry_options(),
        || false,
    )
    .unwrap();
    assert!(matches!(
        FiniteMoistureSolution::with_electrolyte_solubility(hydrolysis, WATER_MOLAR_MASS),
        Err(SolutionMoistureError::Invalid(_))
    ));
}
