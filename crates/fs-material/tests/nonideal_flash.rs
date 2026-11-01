//! Original-source multicomponent PT/PH equilibrium; no copied production EOS.
#[path = "common/nonideal_flash.rs"]
mod common;
#[path = "nonideal_flash/energy_volume.rs"]
mod energy_volume;
#[path = "nonideal_flash/parcel.rs"]
mod parcel;
#[path = "nonideal_flash/phase_equilibrium.rs"]
mod phase_equilibrium;

use common::{METHANE, MIXTURE_MASS, controls, model, source};
use fs_material::thermochemical_gas::nonideal::ThermochemicalPengRobinsonModel;
use fs_material::thermochemical_gas::nonideal::flash::{
    ResolvedThermochemicalFlashState, ThermochemicalFlashError as E,
    enthalpy::FlashEnthalpyControls,
};
use fs_qty::{Pressure, Temperature};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;
use fs_thermochem::eos::peng_robinson::PengRobinsonPhaseRootV1;

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}, difference {}",
        actual - expected
    );
}
fn pt(
    source: &ThermochemicalPengRobinsonModel,
    t: f64,
    p: f64,
) -> ResolvedThermochemicalFlashState {
    source
        .evaluate_two_phase_flash(
            Temperature::new(t),
            Pressure::new(p),
            controls().flash_options,
            || false,
        )
        .unwrap()
}

// Reconstruct each full phase from the original individual species, then mix
// those phase values. Production instead retains the exact feed reference once
// and adds the phase-redistribution correction before converting to mass basis.
fn independent_phase_calorics(
    source: &ThermochemicalPengRobinsonModel,
    root: &PengRobinsonPhaseRootV1,
    fractions: &[f64],
    t: f64,
    p: f64,
) -> (f64, f64, f64) {
    let mut enthalpy = 0.0;
    let mut entropy = 0.0;
    for (index, species) in source
        .source()
        .thermodynamics()
        .components()
        .iter()
        .enumerate()
    {
        let original = species.evaluate_properties(Temperature::new(t)).unwrap();
        let amount = fractions[index];
        enthalpy += amount * original.h().value();
        entropy += amount
            * (original.s().value()
                - R * (amount * p / species.convention().reference_pressure().value()).ln()
                - R * root.log_fugacity_coefficients()[index]);
    }
    enthalpy += root.residual_enthalpy_j_mol();
    entropy += root.residual_enthalpy_j_mol() / t;
    (
        enthalpy,
        enthalpy - p * root.molar_volume_m3_per_mol(),
        entropy,
    )
}

#[test]
fn g1_binary_pt_closes_independent_phase_energy_entropy_and_inventory() {
    let source = source();
    let state = pt(&source, 220.0, 3e6);
    let flash = state.flash();
    // Independent 50-digit simultaneous solve of both PR cubic equations and
    // both fugacity equalities, not production Wilson/Rachford-Rice iteration.
    close(flash.vapor_fraction_mol(), 0.302_585_915_611_944_35, 2e-10);
    close(
        flash.liquid_mole_fractions()[METHANE],
        0.364_782_477_343_904_27,
        2e-10,
    );
    close(
        flash.vapor_mole_fractions()[METHANE],
        0.811_655_631_973_835_5,
        2e-10,
    );
    close(
        flash.liquid().molar_volume_m3_per_mol(),
        0.000_056_108_753_799_670_64,
        2e-14,
    );
    close(
        flash.vapor().molar_volume_m3_per_mol(),
        0.000_453_008_645_446_556_15,
        2e-13,
    );
    close(
        state.molar_volume_m3_mol(),
        0.000_176_205_070_919_925,
        2e-13,
    );
    close(state.density_kg_m3(), 130.846_120_827_461_9, 2e-7);
    close(state.vapor_mass_fraction(), 0.245_214_412_632_239_24, 2e-10);
    close(state.vapor_volume_fraction(), 0.777_923_331_303_36, 2e-10);
    close(
        state.liquid_mass_fraction() + state.vapor_mass_fraction(),
        1.0,
        3e-14,
    );
    close(
        state.specific_enthalpy_j_kg(),
        -16_172.250_165_027_297,
        3e-5,
    );
    close(state.specific_entropy_j_kg_k(), 5507.749_287_855_677, 2e-7);
    close(
        (state.specific_entropy_j_kg_k()
            - state.reference().properties().s().value() / MIXTURE_MASS)
            * MIXTURE_MASS,
        -29.879_185_632_927_96,
        5e-9,
    );
    assert_eq!(state.constitutive_law(), source.mechanics());
    assert_eq!(
        state.reference(),
        &source
            .source()
            .thermodynamics()
            .evaluate(Temperature::new(220.0), Pressure::new(3e6))
            .unwrap()
    );
    let liquid = independent_phase_calorics(
        &source,
        flash.liquid(),
        flash.liquid_mole_fractions(),
        220.0,
        3e6,
    );
    let vapor = independent_phase_calorics(
        &source,
        flash.vapor(),
        flash.vapor_mole_fractions(),
        220.0,
        3e6,
    );
    let beta = flash.vapor_fraction_mol();
    close(
        state.specific_enthalpy_j_kg(),
        ((1.0 - beta) * liquid.0 + beta * vapor.0) / MIXTURE_MASS,
        1e-8,
    );
    close(
        state.specific_internal_energy_j_kg(),
        ((1.0 - beta) * liquid.1 + beta * vapor.1) / MIXTURE_MASS,
        1e-8,
    );
    close(
        state.specific_entropy_j_kg_k(),
        ((1.0 - beta) * liquid.2 + beta * vapor.2) / MIXTURE_MASS,
        1e-9,
    );
    close(
        state.specific_enthalpy_j_kg() - state.specific_internal_energy_j_kg(),
        state.pressure_pa() / state.density_kg_m3(),
        1e-8,
    );
    for (index, feed) in source.mechanics().mole_fractions().iter().enumerate() {
        close(
            (1.0 - beta) * flash.liquid_mole_fractions()[index]
                + beta * flash.vapor_mole_fractions()[index],
            *feed,
            1e-12,
        );
    }
}

#[test]
fn g1_isenthalpic_pressure_drop_matches_independent_simultaneous_equilibrium() {
    let source = source();
    let initial = pt(&source, 220.0, 3e6);
    let control = FlashEnthalpyControls::new([Temperature::new(200.0), Temperature::new(220.0)]);
    let solved = source
        .two_phase_flash_pressure_enthalpy_increment(
            Temperature::new(220.0),
            Pressure::new(3e6),
            Pressure::new(2e6),
            0.0,
            control,
            || false,
        )
        .unwrap();
    let state = solved.state();
    // Independent simultaneous five-variable cubic/fugacity/enthalpy Newton
    // solve, 50 decimal digits, with pressure continuation from 3 to 2 MPa.
    close(state.temperature_k(), 209.606_128_901_491_32, 2e-7);
    close(
        state.flash().vapor_fraction_mol(),
        0.392_287_810_793_803_86,
        2e-9,
    );
    close(
        state.flash().liquid_mole_fractions()[METHANE],
        0.290_967_906_203_411_7,
        2e-9,
    );
    close(
        state.flash().vapor_mole_fractions()[METHANE],
        0.823_821_816_126_350_2,
        2e-9,
    );
    close(state.density_kg_m3(), 74.238_856_531_517_66, 2e-6);
    close(
        solved.entropy_change_j_kg_k().unwrap(),
        47.320_830_152_714_31,
        2e-6,
    );
    close(
        solved.accepted_increment_j_kg().unwrap(),
        0.0,
        control.enthalpy_tolerance_j_kg,
    );
    close(
        state.specific_enthalpy_j_kg(),
        initial.specific_enthalpy_j_kg(),
        control.enthalpy_tolerance_j_kg * 1.01,
    );
    close(
        state.specific_entropy_j_kg_k() - initial.specific_entropy_j_kg_k(),
        solved.entropy_change_j_kg_k().unwrap(),
        3e-10,
    );
    assert!(solved.entropy_change_j_kg_k().unwrap() > 0.0);
    assert!(solved.flash_evaluations() <= control.maximum_flash_evaluations);
    assert!(solved.eos_evaluation_upper_bound() <= control.maximum_eos_evaluation_bound().unwrap());
    let absolute = source
        .two_phase_flash_pressure_enthalpy(
            Pressure::new(2e6),
            -16_172.250_165_027_297,
            control,
            || false,
        )
        .unwrap();
    close(
        absolute.state().temperature_k(),
        209.606_128_901_491_32,
        2e-7,
    );
    assert_eq!(absolute.accepted_increment_j_kg(), None);
    assert_eq!(absolute.entropy_change_j_kg_k(), None);
}

#[test]
fn g3_small_heat_and_entropy_changes_survive_independent_species_reference_shifts() {
    let ordinary = source();
    let shifted = model([1e20, -3e20], [2e20, -1e20], [100.0, 450.0]);
    for owner in [&ordinary, &shifted] {
        let initial = pt(owner, 220.0, 3e6);
        let mut one_call = controls();
        one_call.maximum_flash_evaluations = 1;
        let fixed = owner
            .two_phase_flash_pressure_enthalpy_increment(
                Temperature::new(220.0),
                Pressure::new(3e6),
                Pressure::new(3e6),
                0.0,
                one_call,
                || false,
            )
            .unwrap();
        assert_eq!(fixed.state(), &initial);
        assert_eq!(fixed.accepted_increment_j_kg(), Some(0.0));
        assert_eq!(fixed.entropy_change_j_kg_k(), Some(0.0));
        assert_eq!(fixed.residual_j_kg(), 0.0);
        assert_eq!(fixed.flash_evaluations(), 1);
    }
    let solve = |source: &ThermochemicalPengRobinsonModel| {
        source
            .two_phase_flash_pressure_enthalpy_increment(
                Temperature::new(220.0),
                Pressure::new(3e6),
                Pressure::new(3e6),
                10.0,
                controls(),
                || false,
            )
            .unwrap()
    };
    let normal = solve(&ordinary);
    let changed = solve(&shifted);
    assert!(normal.state().temperature_k() > 220.0);
    assert!(
        normal.state().flash().vapor_fraction_mol()
            > pt(&ordinary, 220.0, 3e6).flash().vapor_fraction_mol()
    );
    assert_eq!(
        changed.state().temperature_k(),
        normal.state().temperature_k()
    );
    assert_eq!(changed.state().flash(), normal.state().flash());
    assert_eq!(
        changed.accepted_increment_j_kg(),
        normal.accepted_increment_j_kg()
    );
    assert_eq!(
        changed.entropy_change_j_kg_k(),
        normal.entropy_change_j_kg_k()
    );
    close(
        changed.accepted_increment_j_kg().unwrap(),
        10.0,
        controls().enthalpy_tolerance_j_kg,
    );
    assert!(matches!(
        shifted.two_phase_flash_pressure_enthalpy(
            Pressure::new(3e6),
            pt(&shifted, 220.0, 3e6).specific_enthalpy_j_kg(),
            controls(),
            || false
        ),
        Err(E::Invalid(_))
    ));
    let restored = shifted
        .two_phase_flash_pressure_enthalpy_increment(
            Temperature::new(changed.state().temperature_k()),
            Pressure::new(3e6),
            Pressure::new(3e6),
            -10.0,
            controls(),
            || false,
        )
        .unwrap();
    close(restored.state().temperature_k(), 220.0, 2e-8);
}

#[test]
fn g0_species_reference_offsets_are_counted_once_at_the_exact_feed_composition() {
    let original = pt(&source(), 220.0, 3e6);
    let changed = pt(
        &model([500.0, -100.0], [2.0, -3.0], [100.0, 450.0]),
        220.0,
        3e6,
    );
    let mass_weighted = |methane: f64, ethane: f64| {
        (0.5 * 0.016_042_46 * methane + 0.5 * 0.030_069_04 * ethane) / MIXTURE_MASS
    };
    assert_eq!(changed.flash(), original.flash());
    close(
        changed.specific_enthalpy_j_kg() - original.specific_enthalpy_j_kg(),
        mass_weighted(500.0, -100.0),
        1e-10,
    );
    close(
        changed.specific_internal_energy_j_kg() - original.specific_internal_energy_j_kg(),
        mass_weighted(500.0, -100.0),
        1e-10,
    );
    close(
        changed.specific_entropy_j_kg_k() - original.specific_entropy_j_kg_k(),
        mass_weighted(2.0, -3.0),
        3e-12,
    );
}

#[test]
fn g4_phase_source_brackets_and_whole_call_work_limits_refuse_explicitly() {
    let source = source();
    assert!(matches!(
        source.evaluate_two_phase_flash(
            Temperature::new(400.0),
            Pressure::new(1e5),
            controls().flash_options,
            || false
        ),
        Err(E::Phase(_))
    ));
    assert!(matches!(pt_source_domain_refusal(), Err(E::Material(_))));
    let mut control = controls();
    control.maximum_flash_evaluations = 2;
    assert!(matches!(
        source.two_phase_flash_pressure_enthalpy(
            Pressure::new(3e6),
            -16_172.250_165_027_297,
            control,
            || false
        ),
        Err(E::WorkLimit {
            flash_evaluations: 2
        })
    ));
    // The reference PT solve consumes the same total budget as the endpoints.
    control.maximum_flash_evaluations = 3;
    assert!(matches!(
        source.two_phase_flash_pressure_enthalpy_increment(
            Temperature::new(220.0),
            Pressure::new(3e6),
            Pressure::new(3e6),
            10.0,
            control,
            || false
        ),
        Err(E::WorkLimit {
            flash_evaluations: 3
        })
    ));
    assert!(matches!(
        source.two_phase_flash_pressure_enthalpy(Pressure::new(3e6), 1e8, controls(), || false),
        Err(E::OutsideEnthalpyBracket)
    ));
    let invalid_phase =
        FlashEnthalpyControls::new([Temperature::new(210.0), Temperature::new(400.0)]);
    assert!(matches!(
        source.two_phase_flash_pressure_enthalpy(
            Pressure::new(3e6),
            -16_172.250_165_027_297,
            invalid_phase,
            || false
        ),
        Err(E::Phase(_))
    ));
    control.maximum_flash_evaluations = 0;
    assert!(control.validate().is_err());
}
fn pt_source_domain_refusal() -> Result<ResolvedThermochemicalFlashState, E> {
    model([0.0; 2], [0.0; 2], [225.0, 250.0]).evaluate_two_phase_flash(
        Temperature::new(220.0),
        Pressure::new(3e6),
        controls().flash_options,
        || false,
    )
}

#[test]
fn g4_final_cancellation_publishes_no_caloric_or_enthalpy_state() {
    let source = source();
    let mut polls = 0;
    source
        .evaluate_two_phase_flash(
            Temperature::new(220.0),
            Pressure::new(3e6),
            controls().flash_options,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut replay = 0;
    assert!(matches!(
        source.evaluate_two_phase_flash(
            Temperature::new(220.0),
            Pressure::new(3e6),
            controls().flash_options,
            || {
                replay += 1;
                replay == polls
            }
        ),
        Err(E::Cancelled)
    ));
    assert_eq!(replay, polls);
    let mut polls = 0;
    source
        .two_phase_flash_pressure_enthalpy(
            Pressure::new(3e6),
            -16_172.250_165_027_297,
            controls(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut replay = 0;
    assert!(matches!(
        source.two_phase_flash_pressure_enthalpy(
            Pressure::new(3e6),
            -16_172.250_165_027_297,
            controls(),
            || {
                replay += 1;
                replay == polls
            }
        ),
        Err(E::Cancelled)
    ));
    assert_eq!(replay, polls);
}
