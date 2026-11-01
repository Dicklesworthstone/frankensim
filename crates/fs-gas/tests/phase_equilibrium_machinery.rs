//! G0/G1/G3/G4: actual shaft work and entropy across original phase equilibrium.
#![cfg(feature = "material")]

// Reuse only source input data. Independent cubic/fugacity caloric reference
// values below do not call the production phase or temperature-root engine.
#[allow(dead_code)]
#[path = "../../fs-material/tests/common/nonideal_flash.rs"]
mod common;

use common::{METHANE, model, source};
use fs_gas::material::phase_machinery::{
    MAX_MACHINE_PT_EVALUATIONS, MixtureMachineControl, MixtureMachineError as E,
    MixtureMachineKind as Kind, MixtureMachineStage, MixtureMachineStep, adiabatic_mixture_machine,
};
use fs_material::thermochemical_gas::nonideal::ThermochemicalPengRobinsonModel;
use fs_material::thermochemical_gas::nonideal::flash::ThermochemicalFlashError;
use fs_qty::{MassFlowRate, Pressure, Temperature};

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}; difference {}",
        actual - expected,
    );
}
fn controls() -> MixtureMachineControl {
    MixtureMachineControl::new([Temperature::new(200.0), Temperature::new(235.0)])
}
fn expand(
    model: &ThermochemicalPengRobinsonModel,
    efficiency: f64,
    control: MixtureMachineControl,
    cancelled: impl FnMut() -> bool,
) -> Result<MixtureMachineStep, E> {
    adiabatic_mixture_machine(
        model,
        Kind::Expander,
        Temperature::new(220.0),
        Pressure::new(3e6),
        Pressure::new(2e6),
        efficiency,
        control,
        cancelled,
    )
}

#[test]
fn g1_reversible_two_phase_expansion_matches_independent_equilibrium_and_work() {
    let original = source();
    let result = expand(&original, 1.0, controls(), || false).unwrap();
    // Independently solve both PR cubics, equal component fugacities, material
    // balance and the phase-resolved entropy equation at the outlet pressure.
    close(
        result.outlet().temperature_k(),
        208.202_765_109_369_35,
        3e-7,
    );
    close(result.outlet().density_kg_m3(), 78.021_034_395_541_12, 2e-6);
    let phase = result.outlet().phase().two_phase().unwrap();
    close(phase.vapor_fraction_mol(), 0.372_659_298_790_543_3, 3e-9);
    close(
        phase.liquid_mole_fractions()[METHANE],
        0.301_426_580_035_008_8,
        3e-9,
    );
    close(
        phase.vapor_mole_fractions()[METHANE],
        0.834_281_712_348_777_5,
        3e-9,
    );
    close(result.shaft_work_to_fluid_j_kg, -9885.548_034_785_53, 1e-4);
    close(result.entropy_production_j_kg_k, 0.0, 1e-8);
    assert_eq!(result.outlet(), result.isentropic_outlet());
    assert_eq!(result.enthalpy_residual_j_kg, 0.0);
    assert_eq!(result.outlet().constitutive_law(), original.mechanics());
    let direct = original
        .phase_equilibrium_pressure_entropy_increment(
            Temperature::new(220.0),
            Pressure::new(3e6),
            Pressure::new(2e6),
            0.0,
            controls().isentropic,
            || false,
        )
        .unwrap();
    assert_eq!(result.pt_evaluations, direct.pt_evaluations());
    assert_eq!(
        result.eos_evaluation_upper_bound,
        direct.eos_evaluation_upper_bound(),
    );
    let power = result.shaft_power_to_fluid(MassFlowRate::new(2.0)).unwrap();
    close(power.value(), 2.0 * result.shaft_work_to_fluid_j_kg, 1e-10);
    assert!(power.value() < 0.0);
    // Shaft work uses stream enthalpy, not Delta u with another p*v debit.
    assert!((result.shaft_work_to_fluid_j_kg - -12591.965_499_954_124).abs() > 2000.0);
}

#[test]
fn g1_efficiency_changes_actual_work_and_generates_entropy_in_both_directions() {
    let original = source();
    let reversible = expand(&original, 1.0, controls(), || false).unwrap();
    let expanded = expand(&original, 0.8, controls(), || false).unwrap();
    close(
        expanded.shaft_work_to_fluid_j_kg,
        0.8 * reversible.shaft_work_to_fluid_j_kg,
        1e-5,
    );
    assert!(expanded.outlet().temperature_k() > reversible.outlet().temperature_k());
    assert!(expanded.entropy_production_j_kg_k > 0.0);
    assert!(expanded.enthalpy_residual_j_kg.abs() <= 1e-5);
    assert!(expanded.pt_evaluations > reversible.pt_evaluations);
    let compress = |efficiency| {
        adiabatic_mixture_machine(
            &original,
            Kind::Compressor,
            Temperature::new(reversible.outlet().temperature_k()),
            Pressure::new(2e6),
            Pressure::new(3e6),
            efficiency,
            controls(),
            || false,
        )
        .unwrap()
    };
    let reverse = compress(1.0);
    close(reverse.outlet().temperature_k(), 220.0, 3e-7);
    close(
        reverse.shaft_work_to_fluid_j_kg + reversible.shaft_work_to_fluid_j_kg,
        0.0,
        1e-4,
    );
    let compressed = compress(0.8);
    close(
        compressed.shaft_work_to_fluid_j_kg,
        reverse.shaft_work_to_fluid_j_kg / 0.8,
        1e-5,
    );
    assert!(compressed.outlet().temperature_k() > reverse.outlet().temperature_k());
    assert!(compressed.entropy_production_j_kg_k > 0.0);
    assert!(compressed.shaft_work_to_fluid_j_kg > 0.0);
    assert!(compressed.outlet().phase().two_phase().is_some());
}

#[test]
fn g3_independent_species_reference_datums_do_not_change_machine_work() {
    let plain = source();
    let shifted = model([1e22, -3e21], [2e18, -7e17], [100.0, 450.0]);
    let a = expand(&plain, 0.8, controls(), || false).unwrap();
    let b = expand(&shifted, 0.8, controls(), || false).unwrap();
    assert_eq!(a.outlet().temperature_k(), b.outlet().temperature_k());
    assert_eq!(a.outlet().phase(), b.outlet().phase());
    assert_eq!(a.shaft_work_to_fluid_j_kg, b.shaft_work_to_fluid_j_kg);
    assert_eq!(a.entropy_production_j_kg_k, b.entropy_production_j_kg_k);
    assert_eq!(a.enthalpy_residual_j_kg, b.enthalpy_residual_j_kg);
    assert_eq!(a.pt_evaluations, b.pt_evaluations);
    assert_eq!(a.eos_evaluation_upper_bound, b.eos_evaluation_upper_bound);
}

#[test]
fn g0_equal_pressure_is_exact_zero_work_with_one_original_pt_query() {
    let original = source();
    for kind in [Kind::Compressor, Kind::Expander] {
        let mut control = controls();
        control.maximum_pt_evaluations = 1;
        let result = adiabatic_mixture_machine(
            &original,
            kind,
            Temperature::new(220.0),
            Pressure::new(3e6),
            Pressure::new(3e6),
            0.65,
            control,
            || false,
        )
        .unwrap();
        assert_eq!(result.outlet().temperature_k(), 220.0);
        assert_eq!(result.shaft_work_to_fluid_j_kg, 0.0);
        assert_eq!(result.entropy_production_j_kg_k, 0.0);
        assert_eq!(result.enthalpy_residual_j_kg, 0.0);
        assert_eq!(result.pt_evaluations, 1);
    }
}

#[test]
fn g0_shared_work_budget_counts_both_reference_states_and_both_roots() {
    let original = source();
    let accepted = expand(&original, 0.8, controls(), || false).unwrap();
    let mut control = controls();
    control.maximum_pt_evaluations = accepted.pt_evaluations;
    assert_eq!(expand(&original, 0.8, control, || false).unwrap(), accepted);
    control.maximum_pt_evaluations -= 1;
    match expand(&original, 0.8, control, || false).unwrap_err() {
        E::Equilibrium {
            stage: MixtureMachineStage::ActualOutlet,
            completed_pt_evaluations,
            source: ThermochemicalFlashError::WorkLimit { flash_evaluations },
        } => assert_eq!(
            completed_pt_evaluations + flash_evaluations,
            control.maximum_pt_evaluations,
        ),
        other => panic!("unexpected whole-machine budget refusal: {other:?}"),
    }
}

#[test]
fn g0_invalid_machine_inputs_and_conflicting_phase_policy_refuse() {
    let original = source();
    for efficiency in [0.0, -0.1, 1.00001, f64::NAN] {
        assert!(matches!(
            expand(&original, efficiency, controls(), || false),
            Err(E::Invalid(_)),
        ));
    }
    assert!(matches!(
        adiabatic_mixture_machine(
            &original,
            Kind::Compressor,
            Temperature::new(220.0),
            Pressure::new(3e6),
            Pressure::new(2e6),
            0.8,
            controls(),
            || false,
        ),
        Err(E::Invalid(_)),
    ));
    let mut control = controls();
    control.outlet.stability.tangent_plane_tolerance *= 2.0;
    assert!(matches!(
        expand(&original, 0.8, control, || false),
        Err(E::Invalid(_))
    ));
    control = controls();
    control.maximum_pt_evaluations = MAX_MACHINE_PT_EVALUATIONS + 1;
    assert!(matches!(
        expand(&original, 0.8, control, || false),
        Err(E::Invalid(_))
    ));
    let result = expand(&original, 0.8, controls(), || false).unwrap();
    assert!(matches!(
        result.shaft_power_to_fluid(MassFlowRate::new(-1.0)),
        Err(E::Invalid(_)),
    ));
}

#[test]
fn g4_final_cancellation_never_publishes_an_outlet_and_replay_is_identical() {
    let original = source();
    let before = original.clone();
    let mut calls = 0;
    let expected = expand(&original, 0.8, controls(), || {
        calls += 1;
        false
    })
    .unwrap();
    for cutoff in [1, calls] {
        let mut current = 0;
        assert_eq!(
            expand(&original, 0.8, controls(), || {
                current += 1;
                current == cutoff
            }),
            Err(E::Cancelled),
        );
        assert_eq!(original, before);
    }
    assert_eq!(
        expand(&original, 0.8, controls(), || false).unwrap(),
        expected
    );
}
