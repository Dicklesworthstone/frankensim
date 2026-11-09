use super::*;
use crate::osmosis::solubility::electrolyte::ElectrolyteSolubilityError;
use fs_qty::{AmountConcentration, Dimensionless};
use fs_thermochem::{
    UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
    stoichiometry::thermodynamics::acid_base::{
        AcidBaseError,
        precipitation::kinetics::{MineralGrowthPolicyV1, MineralSurfaceRateV1},
    },
};

#[path = "../tests/fixture.rs"]
mod fixture;
use fixture::*;

const T: f64 = 298.15;
const VOLUME: f64 = 1e-3;
const VM: f64 = 18e-6;

fn basis() -> SolventActivityBasis {
    SolventActivityBasis::new(VM, AmountConcentration::new(1000.0)).unwrap()
}
fn aqueous_options() -> SpeciationOptions {
    let mut options = SpeciationOptions::default();
    options.absolute_charge_tolerance = Amount::new(1e-18);
    options.relative_charge_tolerance = Dimensionless::new(1e-12);
    options.ph_tolerance = Dimensionless::new(1e-13);
    options.max_iterations = 256;
    options
}
fn controls() -> MineralKineticControlsV1 {
    MineralKineticControlsV1 {
        aqueous: aqueous_options(),
        initial_step: Time::new(0.1),
        minimum_step: Time::new(1e-7),
        maximum_step: Time::new(0.25),
        relative_tolerance: Dimensionless::new(1e-5),
        amount_absolute_tolerance: Amount::new(1e-12),
        extent_absolute_tolerance: Amount::new(1e-12),
        maximum_stage_iterations: 128,
        maximum_attempts: 8192,
    }
}
fn salt_material(
    dissolved: f64,
    solid: f64,
    family: f64,
    policy: MineralGrowthPolicyV1,
) -> KineticElectrolyteSolubilitySolution {
    let rate = MineralSurfaceRateV1::new(1e-4, policy).unwrap();
    let model = MineralKineticsV1::from_salt(salt(T, 4e-6), rate).unwrap();
    KineticElectrolyteSolubilitySolution::new(
        model,
        basis(),
        Volume::new(VOLUME),
        Amount::new(dissolved),
        Amount::new(solid),
        vec![(FAMILY.into(), Amount::new(family))],
        vec![],
        aqueous_options(),
        || false,
    )
    .unwrap()
}
fn near(actual: f64, expected: f64, absolute: f64) {
    assert!((actual - expected).abs() <= absolute, "{actual:e} != {expected:e}");
}
fn snapshot(material: &KineticElectrolyteSolubilitySolution) -> Vec<u64> {
    let mut result = vec![
        material.solution().solvent_volume_m3().to_bits(),
        material.dissolved_cation().value().to_bits(),
        material.solid_amount().value().to_bits(),
        material.aqueous().ph.value().to_bits(),
        material.component_total().value().to_bits(),
        material.state().amount_at_unit_activity().value().to_bits(),
    ];
    result.extend(material.solution().solute_moles().iter().map(|x| x.to_bits()));
    result.extend(material.family_totals().iter().map(|(_, n)| n.value().to_bits()));
    result
}
fn check_salt_balances(material: &KineticElectrolyteSolubilitySolution) {
    near(
        material.dissolved_cation().value() + material.solid_amount().value(),
        material.component_total().value(),
        2e-14,
    );
    near(
        material.dissolved_family().unwrap().value() + 2.0 * material.solid_amount().value(),
        material.family_totals()[0].1.value(),
        3e-14,
    );
    let dissolved_particles = material.dissolved_cation().value()
        + material.dissolved_family().unwrap().value()
        + material.aqueous().proton_amount.value()
        + material.aqueous().hydroxide_amount.value();
    near(
        material.solution().osmotic_pressure_pa().unwrap(),
        R * T * dissolved_particles / material.solution().solvent_volume_m3(),
        1e-7,
    );
}

#[test]
fn actual_salt_phases_dissolve_at_finite_rate_and_publish_only_speciated_particles() {
    let original = salt_material(0.001, 0.109, 0.022, MineralGrowthPolicyV1::ExistingSolidOnly);
    let before = snapshot(&original);
    assert_eq!(original.dissolved_cation(), Amount::new(0.001));
    assert_eq!(original.solid_amount(), Amount::new(0.109));
    assert_eq!(original.dissolved_family(), Some(Amount::new(0.022)));
    near(original.component_total().value(), 0.11, 1e-16);
    near(original.family_totals()[0].1.value(), 0.24, 1e-16);
    assert!(original.aqueous().ph.value() < 7.0);
    assert!(original.state().log_saturation_ratio().unwrap().value() < 0.0);
    assert_eq!(original.species(), [id(H), id(OH), id(HC), id(CN), id(CA)]);
    assert_eq!(original.solvent(), &id(WATER));
    assert!(!original.species().contains(&id(SALT)));
    assert!(!original.species().contains(&id(WATER)));
    check_salt_balances(&original);

    let mut control = controls();
    control.initial_step = Time::new(0.05);
    control.maximum_step = Time::new(0.05);
    let step = original.advance(Time::new(50.0), control, || false).unwrap();
    let endpoint = step.solution();
    assert!(endpoint.dissolved_cation().value() > 0.001);
    assert!(endpoint.dissolved_cation().value() < 0.01); // independently manufactured equilibrium
    assert!(endpoint.solid_amount().value() < 0.109 && endpoint.solid_amount().value() > 0.1);
    assert!(endpoint.aqueous().ph.value() > original.aqueous().ph.value());
    // Independent charge-cubic + RK4 oracle for this authored source/rate,
    // evaluated at dt=.01 and .005 s (dissolved amounts agree within 2e-17).
    // This checks the finite trajectory, not just its eventual equilibrium.
    near(endpoint.dissolved_cation().value(), 0.00570340484820815, 3e-7);
    near(endpoint.solid_amount().value(), 0.10429659515179185, 3e-7);
    near(endpoint.aqueous().ph.value(), 6.75614127313096, 1e-4);
    assert_eq!(endpoint.component_total(), original.component_total());
    assert_eq!(endpoint.family_totals(), original.family_totals());
    assert!(std::ptr::eq(endpoint.species(), original.species()));
    assert!(std::ptr::eq(endpoint.model().source(), original.model().source()));
    assert!(step.kinetics().transferred_to_solution().value() > 0.0);
    assert!(step.kinetics().accepted_steps() > 0 && step.kinetics().aqueous_solves() > 0);
    assert_eq!(endpoint.dissolved_cation(), step.kinetics().final_state().dissolved_cation());
    check_salt_balances(endpoint);
    assert_eq!(snapshot(&original), before);
}

#[test]
fn volume_changes_freeze_crystals_then_native_kinetics_resolves_actual_saturation() {
    let original = salt_material(0.01, 0.1, 0.04, MineralGrowthPolicyV1::ExistingSolidOnly);
    let before = snapshot(&original);
    let concentrated = original.with_solvent_volume_m3(5e-4).unwrap();
    assert_eq!(concentrated.solid_amount(), original.solid_amount());
    assert_eq!(concentrated.dissolved_cation(), original.dissolved_cation());
    assert_eq!(concentrated.dissolved_family(), original.dissolved_family());
    assert_eq!(concentrated.component_total(), original.component_total());
    assert_eq!(concentrated.family_totals(), original.family_totals());
    assert_eq!(concentrated.state().amount_at_unit_activity(), Amount::new(0.5));
    assert!(concentrated.state().log_saturation_ratio().unwrap().value() > 0.0);
    assert!(std::ptr::eq(concentrated.model().source(), original.model().source()));
    assert!(std::ptr::eq(concentrated.species(), original.species()));
    check_salt_balances(&concentrated);
    let growth = concentrated.advance(Time::new(0.5), controls(), || false).unwrap();
    assert!(growth.solution().solid_amount().value() > 0.1);
    assert!(growth.kinetics().transferred_to_solution().value() < 0.0);
    check_salt_balances(growth.solution());

    let diluted = growth.solution().with_solvent_volume_m3(0.01).unwrap();
    assert_eq!(diluted.solid_amount(), growth.solution().solid_amount());
    assert_eq!(diluted.dissolved_cation(), growth.solution().dissolved_cation());
    assert!(diluted.state().log_saturation_ratio().unwrap().value() < 0.0);
    let dissolution = diluted.advance(Time::new(0.5), controls(), || false).unwrap();
    assert!(dissolution.solution().solid_amount().value() < diluted.solid_amount().value());
    assert!(dissolution.kinetics().transferred_to_solution().value() > 0.0);
    assert_eq!(dissolution.solution().family_totals(), original.family_totals());
    check_salt_balances(dissolution.solution());
    assert_eq!(snapshot(&original), before);
}

#[test]
fn native_hydrolysis_adoption_preserves_ph_counterions_and_fixed_solvent_convention() {
    let rate = MineralSurfaceRateV1::new(1e-6, MineralGrowthPolicyV1::ExistingSolidOnly).unwrap();
    let model = MineralKineticsV1::from_hydrolysis(hydrolysis(T, 1e-9), rate).unwrap();
    let chloride = SpectatorIon { species: id("Cl[-](aq)"), charge: -1, amount: Amount::new(2e-5) };
    let state = model.prepare(
        Amount::new(2e-6), Amount::new(1.08e-4), vec![], vec![chloride.clone()],
        Amount::new(1.0), aqueous_options(), || false,
    ).unwrap();
    assert!(KineticElectrolyteSolubilitySolution::from_state(
        basis(), Volume::new(2e-3), state.clone(), || false,
    ).is_err());
    let original = KineticElectrolyteSolubilitySolution::from_state(
        basis(), Volume::new(VOLUME), state, || false,
    ).unwrap();
    assert_eq!(original.species(), [id(H), id(OH), id(CA), id("Cl[-](aq)")]);
    assert_eq!(original.spectators(), &[chloride]);
    assert_eq!(original.dissolved_family(), None);
    assert!(original.family_totals().is_empty());
    let advanced = original.advance(Time::new(0.5), controls(), || false).unwrap();
    let endpoint = advanced.solution();
    assert!(endpoint.dissolved_cation().value() > 2e-6 && endpoint.dissolved_cation().value() < 1e-5);
    assert!(endpoint.solid_amount().value() < 1.08e-4);
    assert!(endpoint.aqueous().ph.value() > original.aqueous().ph.value());
    assert_eq!(endpoint.spectators(), original.spectators());
    assert_eq!(endpoint.solution().solvent_volume_m3(), VOLUME);
    assert_eq!(endpoint.solution().temperature_k(), T);
    assert_eq!(endpoint.state().amount_at_unit_activity(), Amount::new(1.0));
    assert!(std::ptr::eq(endpoint.model().source(), original.model().source()));
    near(endpoint.dissolved_cation().value() + endpoint.solid_amount().value(), 1.1e-4, 1e-18);
    let particles = endpoint.dissolved_cation().value() + 2e-5
        + endpoint.aqueous().proton_amount.value() + endpoint.aqueous().hydroxide_amount.value();
    near(endpoint.solution().osmotic_pressure_pa().unwrap(), R * T * particles / VOLUME, 1e-8);
}

#[test]
fn common_ion_dosing_preserves_solid_then_drives_finite_rate_precipitation() {
    let rate = MineralSurfaceRateV1::new(1e-4, MineralGrowthPolicyV1::ExistingSolidOnly).unwrap();
    let model = MineralKineticsV1::from_salt(salt(T, 4e-6), rate).unwrap();
    let sodium = SpectatorIon { species: id("Na+(aq)"), charge: 1, amount: Amount::new(0.0) };
    let original = KineticElectrolyteSolubilitySolution::new(
        model, basis(), Volume::new(VOLUME), Amount::new(0.01), Amount::new(0.1),
        vec![(FAMILY.into(), Amount::new(0.04))], vec![sodium], aqueous_options(), || false,
    ).unwrap();
    let before = snapshot(&original);
    // Add 0.02 mol of the already-declared Na+ and cyanide family. The actual
    // cation/crystals stay put while pH and free-anion activity respond.
    let dosed = original.with_dissolved_analytical_totals(
        Amount::new(0.01), &[Amount::new(0.06)], &[Amount::new(0.02)], || false,
    ).unwrap();
    assert_eq!(dosed.solid_amount(), original.solid_amount());
    assert_eq!(dosed.dissolved_cation(), original.dissolved_cation());
    assert_eq!(dosed.component_total(), original.component_total());
    assert_eq!(dosed.dissolved_family(), Some(Amount::new(0.06)));
    near(dosed.family_totals()[0].1.value(), 0.26, 1e-16);
    assert_eq!(dosed.spectators()[0].amount, Amount::new(0.02));
    assert_eq!(dosed.spectators()[0].charge, 1);
    assert_eq!(dosed.spectators()[0].species, original.spectators()[0].species);
    assert!(dosed.aqueous().ph.value() > original.aqueous().ph.value());
    assert!(dosed.state().log_saturation_ratio().unwrap().value() > 1.0);
    assert!(std::ptr::eq(dosed.species(), original.species()));
    assert!(std::ptr::eq(dosed.model().source(), original.model().source()));
    let advanced = dosed.advance(Time::new(0.5), controls(), || false).unwrap();
    assert!(advanced.solution().solid_amount().value() > 0.1);
    assert!(advanced.solution().dissolved_cation().value() < 0.01);
    assert_eq!(advanced.solution().family_totals(), dosed.family_totals());
    assert_eq!(advanced.solution().spectators(), dosed.spectators());
    near(advanced.solution().dissolved_family().unwrap().value()
        + 2.0 * advanced.solution().solid_amount().value(), 0.26, 3e-14);
    assert!(original.with_dissolved_analytical_totals(
        Amount::new(0.01), &[], &[Amount::new(0.02)], || false,
    ).is_err());
    assert!(original.with_dissolved_analytical_totals(
        Amount::new(0.01), &[Amount::new(0.06)], &[], || false,
    ).is_err());
    assert_eq!(snapshot(&original), before);
}

#[test]
fn nucleation_policy_and_late_cancellation_never_trigger_implicit_equilibrium() {
    let absent = salt_material(0.03, 0.0, 0.08, MineralGrowthPolicyV1::ExistingSolidOnly);
    let blocked = absent.advance(Time::new(0.1), controls(), || false).unwrap();
    assert_eq!(blocked.solution().solid_amount(), Amount::new(0.0));
    assert_eq!(blocked.solution().dissolved_cation(), Amount::new(0.03));
    assert!(blocked.kinetics().precipitation_blocked());
    let seeded = salt_material(0.03, 0.0, 0.08, MineralGrowthPolicyV1::SeededSurface);
    let checkpoint = snapshot(&seeded);
    let mut calls = 0;
    let accepted = seeded.advance(Time::new(0.1), controls(), || {
        calls += 1;
        false
    }).unwrap();
    assert!(accepted.solution().solid_amount().value() > 0.0);
    assert!(accepted.solution().dissolved_cation().value() < 0.03);
    check_salt_balances(accepted.solution());
    let mut observed = 0;
    assert!(matches!(seeded.advance(Time::new(0.1), controls(), || {
        observed += 1;
        observed == calls
    }), Err(ElectrolyteSolubilityError::Speciation(AcidBaseError::Cancelled))));
    assert_eq!(observed, calls);
    assert_eq!(snapshot(&seeded), checkpoint);
    let mut observed = 0;
    assert!(seeded.with_solvent_volume_m3_cancelled(2e-3, || {
        observed += 1;
        observed > 2
    }).is_err());
    assert_eq!(snapshot(&seeded), checkpoint);
    assert_eq!(snapshot(seeded.advance(Time::new(0.1), controls(), || false).unwrap().solution()),
        snapshot(accepted.solution()));
}
