//! G0/G1/G3/G4: original-source entropy/volume closure and pressure-work endpoints.
use super::{close, common};
use common::{model, source};
use fs_material::thermochemical_gas::nonideal::{
    ThermochemicalPengRobinsonModel,
    flash::{
        ThermochemicalFlashError as E,
        equilibrium::{
            ResolvedThermochemicalPhaseState,
            entropy::{
                PhaseEntropyControls,
                volume::{PhaseEntropyVolumeControls, PhaseEntropyVolumeSolution},
            },
        },
    },
};
use fs_qty::{Pressure, Temperature};

fn controls() -> PhaseEntropyVolumeControls {
    PhaseEntropyVolumeControls::new(
        [Pressure::new(1.75e6), Pressure::new(4e6)],
        PhaseEntropyControls::new([Temperature::new(170.0), Temperature::new(330.0)]),
    )
}
fn pt(
    source: &ThermochemicalPengRobinsonModel,
    t: f64,
    p: f64,
) -> ResolvedThermochemicalPhaseState {
    source
        .evaluate_phase_equilibrium(
            Temperature::new(t),
            Pressure::new(p),
            controls().phase.phase_options(),
            || false,
        )
        .unwrap()
}
fn volume(state: &ResolvedThermochemicalPhaseState) -> f64 {
    state.molar_volume_m3_mol() / state.constitutive_law().molar_mass().value()
}
fn change(
    source: &ThermochemicalPengRobinsonModel,
    initial: &ResolvedThermochemicalPhaseState,
    target_volume: f64,
    delta_s: f64,
    controls: PhaseEntropyVolumeControls,
) -> PhaseEntropyVolumeSolution {
    source
        .phase_equilibrium_entropy_volume_increment(
            Temperature::new(initial.temperature_k()),
            Pressure::new(initial.pressure_pa()),
            target_volume,
            delta_s,
            controls,
            || false,
        )
        .unwrap()
}

#[test]
fn g1_sv_recovers_independent_liquid_coexistence_and_vapor_oracles() {
    let source = source();
    // Independent cubic/fugacity thermodynamic oracles retained by g1/g7;
    // target pressure/temperature are withheld from the S,V inversion.
    for (t, s, rho, split) in [
        (180.0, 4589.992_780_575_417, 496.760_385_342_204_13, false),
        (220.0, 5507.749_287_855_677, 130.846_120_827_461_9, true),
        (280.0, 7011.990_134_658_653, 36.471_975_375_322_785, false),
    ] {
        let v = 1.0 / rho;
        let solved = source
            .phase_equilibrium_entropy_volume(v, s, controls(), || false)
            .unwrap();
        close(solved.state().temperature_k(), t, 3e-7);
        close(solved.state().pressure_pa(), 3e6, 0.02);
        close(volume(solved.state()), v, 1e-12);
        close(solved.state().specific_entropy_j_kg_k(), s, 1e-8);
        close(
            solved.entropy_residual_j_kg_k(),
            solved.state().specific_entropy_j_kg_k() - s,
            0.0,
        );
        close(
            solved.volume_residual_m3_kg(),
            volume(solved.state()) - v,
            0.0,
        );
        assert_eq!(solved.state().phase().two_phase().is_some(), split);
        assert_eq!(solved.state().constitutive_law(), source.mechanics());
        assert_eq!(solved.accepted_increment_j_kg_k(), None);
        assert_eq!(solved.enthalpy_change_j_kg(), None);
        assert_eq!(solved.internal_energy_change_j_kg(), None);
        assert!(solved.pressure_evaluations() >= 3);
        assert!(solved.pt_evaluations() > solved.pressure_evaluations());
        assert!(solved.pt_evaluations() <= controls().maximum_pt_evaluations);
        assert!(
            solved.eos_evaluation_upper_bound()
                <= controls().maximum_eos_evaluation_bound().unwrap()
        );
    }
}

#[test]
fn g1_sv_adiabatic_expansion_and_reverse_match_independent_pressure_work() {
    let source = source();
    let initial = pt(&source, 220.0, 3e6);
    // Independent coupled cubic/fugacity/entropy root used by the PS fixture.
    // Only its specific volume is prescribed here; pressure must be recovered.
    let target_v = 1.0 / 78.021_034_395_541_12;
    let expanded = change(&source, &initial, target_v, 0.0, controls());
    close(
        expanded.state().temperature_k(),
        208.202_765_109_369_35,
        3e-7,
    );
    close(expanded.state().pressure_pa(), 2e6, 0.02);
    close(
        expanded.enthalpy_change_j_kg().unwrap(),
        -9885.548_034_785_53,
        1e-4,
    );
    close(
        expanded.internal_energy_change_j_kg().unwrap(),
        -12591.965_499_954_124,
        1e-4,
    );
    close(expanded.accepted_increment_j_kg_k().unwrap(), 0.0, 1e-8);
    close(
        expanded.enthalpy_change_j_kg().unwrap() - expanded.internal_energy_change_j_kg().unwrap(),
        expanded.state().pressure_pa() * volume(expanded.state())
            - initial.pressure_pa() * volume(&initial),
        1e-7,
    );
    assert!(
        expanded.state().vapor_mass_fraction().unwrap() > initial.vapor_mass_fraction().unwrap()
    );
    let restored = change(&source, expanded.state(), volume(&initial), 0.0, controls());
    close(
        restored.state().temperature_k(),
        initial.temperature_k(),
        3e-7,
    );
    close(restored.state().pressure_pa(), initial.pressure_pa(), 0.03);
    close(
        restored.enthalpy_change_j_kg().unwrap() + expanded.enthalpy_change_j_kg().unwrap(),
        0.0,
        2e-4,
    );
    close(
        restored.internal_energy_change_j_kg().unwrap()
            + expanded.internal_energy_change_j_kg().unwrap(),
        0.0,
        2e-4,
    );
}

#[test]
fn g3_sv_phase_crossing_and_energy_changes_ignore_independent_species_datums() {
    let plain = source();
    let shifted = model([1e22, -3e21], [2e18, -7e17], [100.0, 450.0]);
    let start = pt(&plain, 180.0, 3e6);
    let shifted_start = pt(&shifted, 180.0, 3e6);
    let v = 1.0 / 36.471_975_375_322_785;
    let delta_s = 2421.997_354_083_236;
    let a = change(&plain, &start, v, delta_s, controls());
    let b = change(&shifted, &shifted_start, v, delta_s, controls());
    close(a.state().temperature_k(), 280.0, 3e-7);
    close(a.state().pressure_pa(), 3e6, 0.02);
    close(a.accepted_increment_j_kg_k().unwrap(), delta_s, 1e-8);
    close(
        a.enthalpy_change_j_kg().unwrap(),
        550185.794_810_218_8,
        1e-4,
    );
    close(
        a.internal_energy_change_j_kg().unwrap(),
        473969.987_719_206_8,
        1e-4,
    );
    assert_eq!(a.state().temperature_k(), b.state().temperature_k());
    assert_eq!(a.state().pressure_pa(), b.state().pressure_pa());
    assert_eq!(a.state().phase(), b.state().phase());
    assert_eq!(a.accepted_increment_j_kg_k(), b.accepted_increment_j_kg_k());
    assert_eq!(a.enthalpy_change_j_kg(), b.enthalpy_change_j_kg());
    assert_eq!(
        a.internal_energy_change_j_kg(),
        b.internal_energy_change_j_kg()
    );
    assert_eq!(a.pt_evaluations(), b.pt_evaluations());
    let reverse = change(&plain, a.state(), volume(&start), -delta_s, controls());
    close(reverse.state().temperature_k(), 180.0, 5e-7);
    close(reverse.state().pressure_pa(), 3e6, 0.04);
    assert!(matches!(
        shifted.phase_equilibrium_entropy_volume(
            v,
            b.state().specific_entropy_j_kg_k(),
            controls(),
            || false,
        ),
        Err(E::Invalid(_))
    ));
}

#[test]
fn g0_sv_exact_rest_and_whole_call_budgets_include_the_single_reference() {
    let source = source();
    for temperature in [180.0, 220.0, 280.0] {
        let state = pt(&source, temperature, 3e6);
        let mut control = controls();
        control.maximum_pressure_evaluations = 1;
        control.maximum_pt_evaluations = 1;
        control.phase.maximum_pt_evaluations = 1;
        let rest = change(&source, &state, volume(&state), 0.0, control);
        assert_eq!(rest.state(), &state);
        assert_eq!(rest.pressure_evaluations(), 0);
        assert_eq!(rest.pt_evaluations(), 1);
        assert_eq!(rest.accepted_increment_j_kg_k(), Some(0.0));
        assert_eq!(rest.enthalpy_change_j_kg(), Some(0.0));
        assert_eq!(rest.internal_energy_change_j_kg(), Some(0.0));
    }
    let state = pt(&source, 220.0, 3e6);
    let v = 1.0 / 78.021_034_395_541_12;
    let solved = change(&source, &state, v, 0.0, controls());
    let mut cap = controls();
    cap.maximum_pressure_evaluations = solved.pressure_evaluations();
    cap.maximum_pt_evaluations = solved.pt_evaluations();
    assert_eq!(change(&source, &state, v, 0.0, cap), solved);
    cap.maximum_pt_evaluations -= 1;
    assert!(matches!(source.phase_equilibrium_entropy_volume_increment(
        Temperature::new(220.0), Pressure::new(3e6), v, 0.0, cap, || false,
    ), Err(E::WorkLimit { flash_evaluations }) if flash_evaluations == cap.maximum_pt_evaluations));
    cap = controls();
    cap.maximum_pressure_evaluations = solved.pressure_evaluations() - 1;
    assert!(matches!(source.phase_equilibrium_entropy_volume_increment(
        Temperature::new(220.0), Pressure::new(3e6), v, 0.0, cap, || false,
    ), Err(E::PressureWorkLimit { pressure_evaluations }) if pressure_evaluations == cap.maximum_pressure_evaluations));
}

#[test]
fn g0_sv_refuses_unsupported_source_brackets_and_unrepresentable_constraints() {
    let source = source();
    let v = 1.0 / 130.846_120_827_461_9;
    let s = 5507.749_287_855_677;
    for (volume, entropy) in [(0.0, s), (f64::INFINITY, s), (v, f64::NAN), (v, 1e30)] {
        assert!(matches!(
            source.phase_equilibrium_entropy_volume(volume, entropy, controls(), || false),
            Err(E::Invalid(_))
        ));
    }
    assert!(matches!(
        source.phase_equilibrium_entropy_volume(10.0 * v, s, controls(), || false),
        Err(E::OutsideVolumeBracket)
    ));
    assert!(matches!(
        source.phase_equilibrium_entropy_volume(v, 1e6, controls(), || false),
        Err(E::OutsideEntropyBracket)
    ));
    let limited = model([0.0; 2], [0.0; 2], [210.0, 230.0]);
    assert!(matches!(
        limited.phase_equilibrium_entropy_volume(v, s, controls(), || false),
        Err(E::Material(_))
    ));
    let mut control = controls();
    control.maximum_pt_evaluations = 16385;
    assert!(matches!(control.validate(), Err(E::Invalid(_))));
    control = controls();
    control.maximum_pressure_evaluations = 129;
    assert!(matches!(control.validate(), Err(E::Invalid(_))));
}

#[test]
fn g4_sv_final_cancellation_publishes_no_endpoint_and_replays_deterministically() {
    let source = source();
    let original = source.clone();
    let v = 1.0 / 78.021_034_395_541_12;
    let mut polls = 0;
    let accepted = source
        .phase_equilibrium_entropy_volume_increment(
            Temperature::new(220.0),
            Pressure::new(3e6),
            v,
            0.0,
            controls(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    for cancel_at in [1, polls / 2, polls] {
        let mut count = 0;
        assert_eq!(
            source.phase_equilibrium_entropy_volume_increment(
                Temperature::new(220.0),
                Pressure::new(3e6),
                v,
                0.0,
                controls(),
                || {
                    count += 1;
                    count == cancel_at
                },
            ),
            Err(E::Cancelled)
        );
        assert_eq!(source, original);
    }
    let replay = source
        .phase_equilibrium_entropy_volume_increment(
            Temperature::new(220.0),
            Pressure::new(3e6),
            v,
            0.0,
            controls(),
            || false,
        )
        .unwrap();
    assert_eq!(replay, accepted);
}
