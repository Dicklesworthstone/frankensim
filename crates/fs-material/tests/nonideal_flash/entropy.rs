//! Original-source entropy constraints across pressure and phase boundaries.
use super::{close, common, independent_phase_calorics};
use common::{METHANE, MIXTURE_MASS, model, source};
use fs_material::thermochemical_gas::nonideal::{
    ThermochemicalPengRobinsonModel,
    flash::{
        ThermochemicalFlashError as E,
        equilibrium::{
            ResolvedThermochemicalPhaseState,
            entropy::{PhaseEntropyControls, PhaseEntropySolution},
        },
    },
};
use fs_qty::{Pressure, Temperature};

fn controls() -> PhaseEntropyControls {
    PhaseEntropyControls::new([Temperature::new(175.0), Temperature::new(285.0)])
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
            controls().phase_options(),
            || false,
        )
        .unwrap()
}

fn increment(
    source: &ThermochemicalPengRobinsonModel,
    t: f64,
    initial_p: f64,
    final_p: f64,
    change: f64,
    control: PhaseEntropyControls,
) -> PhaseEntropySolution {
    source
        .phase_equilibrium_pressure_entropy_increment(
            Temperature::new(t),
            Pressure::new(initial_p),
            Pressure::new(final_p),
            change,
            control,
            || false,
        )
        .unwrap()
}

// Independent assembly from individual original species and actual phase
// compositions. The production caloric owner instead evaluates the exact feed
// once and adds a redistribution correction; these two paths must agree.
fn reconstructed_calorics(
    source: &ThermochemicalPengRobinsonModel,
    state: &ResolvedThermochemicalPhaseState,
) -> (f64, f64, f64) {
    let t = state.temperature_k();
    let p = state.pressure_pa();
    let molar = if let Some(split) = state.phase().two_phase() {
        let liquid =
            independent_phase_calorics(source, split.liquid(), split.liquid_mole_fractions(), t, p);
        let vapor =
            independent_phase_calorics(source, split.vapor(), split.vapor_mole_fractions(), t, p);
        let beta = split.vapor_fraction_mol();
        (
            (1.0 - beta) * liquid.0 + beta * vapor.0,
            (1.0 - beta) * liquid.1 + beta * vapor.1,
            (1.0 - beta) * liquid.2 + beta * vapor.2,
        )
    } else {
        independent_phase_calorics(
            source,
            state.phase().homogeneous().unwrap(),
            source.mechanics().mole_fractions(),
            t,
            p,
        )
    };
    (
        molar.0 / MIXTURE_MASS,
        molar.1 / MIXTURE_MASS,
        molar.2 / MIXTURE_MASS,
    )
}

#[test]
fn g0_ps_absolute_recovers_independent_homogeneous_and_coexistence_states() {
    let source = source();
    // Independent cubic/fugacity caloric oracles already retained by g1/g7.
    for (t, s, h, split) in [
        (180.0, 4589.992_780_575_417, -203_472.419_617_408_38, false),
        (220.0, 5507.749_287_855_677, -16_172.250_165_027_297, true),
        (280.0, 7011.990_134_658_653, 346_713.375_192_810_4, false),
    ] {
        let solved = source
            .phase_equilibrium_pressure_entropy(Pressure::new(3e6), s, controls(), || false)
            .unwrap();
        let state = solved.state();
        close(state.temperature_k(), t, 3e-7);
        close(state.specific_entropy_j_kg_k(), s, 1e-8);
        close(state.specific_enthalpy_j_kg(), h, 1e-4);
        assert_eq!(state.phase().two_phase().is_some(), split);
        assert_eq!(state.constitutive_law(), source.mechanics());
        assert_eq!(state.phase().options(), controls().phase_options());
        let independent = reconstructed_calorics(&source, state);
        close(state.specific_enthalpy_j_kg(), independent.0, 1e-7);
        close(state.specific_internal_energy_j_kg(), independent.1, 1e-7);
        close(state.specific_entropy_j_kg_k(), independent.2, 1e-8);
        close(
            solved.residual_j_kg_k(),
            state.specific_entropy_j_kg_k() - s,
            0.0,
        );
        assert_eq!(solved.accepted_increment_j_kg_k(), None);
        assert_eq!(solved.enthalpy_change_j_kg(), None);
        assert_eq!(solved.internal_energy_change_j_kg(), None);
        assert!(solved.pt_evaluations() >= 2);
        assert!(solved.pt_evaluations() <= controls().maximum_pt_evaluations);
        assert!(
            solved.eos_evaluation_upper_bound()
                <= controls().maximum_eos_evaluation_bound().unwrap()
        );
    }
}

#[test]
fn g1_ps_reversible_pressure_drop_matches_independent_equilibrium_and_work() {
    let source = source();
    let control = PhaseEntropyControls::new([Temperature::new(200.0), Temperature::new(230.0)]);
    let initial = pt(&source, 220.0, 3e6);
    let expanded = increment(&source, 220.0, 3e6, 2e6, 0.0, control);
    let state = expanded.state();
    let split = state.phase().two_phase().unwrap();
    // Independent SciPy cubic roots and simultaneous fugacity equalities,
    // followed by a bracketed entropy solve, using the original synthetic
    // NASA Cp/R values 3.5 and 4.5. No production PT/PS iteration supplied these.
    close(state.temperature_k(), 208.202_765_109_369_35, 3e-7);
    close(split.vapor_fraction_mol(), 0.372_659_298_790_543_3, 2e-9);
    close(
        split.liquid_mole_fractions()[METHANE],
        0.301_426_580_035_008_8,
        2e-9,
    );
    close(
        split.vapor_mole_fractions()[METHANE],
        0.834_281_712_348_777_5,
        2e-9,
    );
    close(state.density_kg_m3(), 78.021_034_395_541_12, 2e-6);
    close(
        state.specific_enthalpy_j_kg(),
        -26_057.798_199_812_783,
        1e-4,
    );
    close(
        expanded.enthalpy_change_j_kg().unwrap(),
        -9885.548_034_785_53,
        1e-4,
    );
    close(
        expanded.internal_energy_change_j_kg().unwrap(),
        -12_591.965_499_954_124,
        1e-4,
    );
    close(expanded.accepted_increment_j_kg_k().unwrap(), 0.0, 1e-8);
    assert!(expanded.residual_j_kg_k().abs() <= 1e-8);
    assert!(expanded.pt_evaluations() > 1);
    let independent = reconstructed_calorics(&source, state);
    close(independent.2, initial.specific_entropy_j_kg_k(), 2e-8);
    close(
        expanded.enthalpy_change_j_kg().unwrap(),
        independent.0 - initial.specific_enthalpy_j_kg(),
        1e-7,
    );
    // Original state-function work identity, with actual equilibrium volume.
    // There is no extra latent source, ideal-gas gamma, or mixed heat capacity.
    close(
        expanded.enthalpy_change_j_kg().unwrap() - expanded.internal_energy_change_j_kg().unwrap(),
        state.pressure_pa() / state.density_kg_m3()
            - initial.pressure_pa() / initial.density_kg_m3(),
        1e-7,
    );
    let compressed = increment(&source, state.temperature_k(), 2e6, 3e6, 0.0, control);
    close(compressed.state().temperature_k(), 220.0, 3e-7);
    close(
        compressed.state().density_kg_m3(),
        initial.density_kg_m3(),
        2e-6,
    );
    close(compressed.accepted_increment_j_kg_k().unwrap(), 0.0, 1e-8);
    close(
        compressed.enthalpy_change_j_kg().unwrap() + expanded.enthalpy_change_j_kg().unwrap(),
        0.0,
        2e-4,
    );
    close(
        compressed.internal_energy_change_j_kg().unwrap()
            + expanded.internal_energy_change_j_kg().unwrap(),
        0.0,
        2e-4,
    );
}

#[test]
fn g1_ps_entropy_changes_cross_both_phase_boundaries_in_both_directions() {
    let source = source();
    // Independent bubble/dew temperatures at 3 MPa: 205.7464092585643 K,
    // 250.68921833557063 K. These targets enclose both changes in phase count.
    for (initial_t, delta_s, final_t, delta_h, split) in [
        (
            180.0,
            917.756_507_280_260,
            220.0,
            187_300.169_452_381_1,
            true,
        ),
        (
            220.0,
            1504.240_846_802_976,
            280.0,
            362_885.625_357_837_7,
            false,
        ),
        (
            180.0,
            2421.997_354_083_236,
            280.0,
            550_185.794_810_218_8,
            false,
        ),
        (
            280.0,
            -2421.997_354_083_236,
            180.0,
            -550_185.794_810_218_8,
            false,
        ),
    ] {
        let solved = increment(&source, initial_t, 3e6, 3e6, delta_s, controls());
        close(solved.state().temperature_k(), final_t, 3e-7);
        close(solved.accepted_increment_j_kg_k().unwrap(), delta_s, 1e-8);
        close(solved.enthalpy_change_j_kg().unwrap(), delta_h, 1e-4);
        assert_eq!(solved.state().phase().two_phase().is_some(), split);
        assert!(solved.residual_j_kg_k().abs() <= 1e-8);
        assert!(solved.pt_evaluations() >= 3);
        assert!(
            solved.eos_evaluation_upper_bound()
                <= controls().maximum_eos_evaluation_bound().unwrap()
        );
    }
}

#[test]
fn g3_ps_phase_crossing_and_pressure_work_ignore_independent_species_datums() {
    let ordinary = source();
    let shifted = model([1e22, -3e21], [2e18, -7e17], [100.0, 450.0]);
    for (t, initial_p, final_p, delta_s, control) in [
        (180.0, 3e6, 3e6, 917.756_507_280_260, controls()),
        (
            220.0,
            3e6,
            2e6,
            0.0,
            PhaseEntropyControls::new([Temperature::new(200.0), Temperature::new(230.0)]),
        ),
    ] {
        let a = increment(&ordinary, t, initial_p, final_p, delta_s, control);
        let b = increment(&shifted, t, initial_p, final_p, delta_s, control);
        assert_eq!(
            a.state().temperature_k().to_bits(),
            b.state().temperature_k().to_bits()
        );
        assert_eq!(
            a.state().density_kg_m3().to_bits(),
            b.state().density_kg_m3().to_bits()
        );
        assert_eq!(a.state().phase(), b.state().phase());
        assert_eq!(a.accepted_increment_j_kg_k(), b.accepted_increment_j_kg_k());
        assert_eq!(a.enthalpy_change_j_kg(), b.enthalpy_change_j_kg());
        assert_eq!(
            a.internal_energy_change_j_kg(),
            b.internal_energy_change_j_kg()
        );
        assert_eq!(a.pt_evaluations(), b.pt_evaluations());
        assert_eq!(
            a.eos_evaluation_upper_bound(),
            b.eos_evaluation_upper_bound()
        );
    }
    assert!(matches!(
        shifted.phase_equilibrium_pressure_entropy(
            Pressure::new(3e6),
            pt(&shifted, 220.0, 3e6).specific_entropy_j_kg_k(),
            controls(),
            || false,
        ),
        Err(E::Invalid(_))
    ));
}

#[test]
fn g3_ps_small_nonzero_entropy_survives_loose_tolerances_and_large_datums() {
    let ordinary = source();
    let shifted = model([1e22, -3e21], [2e18, -7e17], [100.0, 450.0]);
    for (requested, bracket) in [(1e-4, [220.0, 220.1]), (-1e-4, [219.9, 220.0])] {
        let mut control = PhaseEntropyControls::new(bracket.map(Temperature::new));
        control.entropy_tolerance_j_kg_k = 1.0;
        control.temperature_tolerance_k = 1.0;
        let a = increment(&ordinary, 220.0, 3e6, 3e6, requested, control);
        let b = increment(&shifted, 220.0, 3e6, 3e6, requested, control);
        let accepted = a.accepted_increment_j_kg_k().unwrap();
        close(accepted, requested, 0.5 * requested.abs());
        assert!(accepted * requested > 0.0);
        assert!((a.state().temperature_k() - 220.0) * requested > 0.0);
        assert!(a.enthalpy_change_j_kg().unwrap() * requested > 0.0);
        assert_eq!(
            a.state().temperature_k().to_bits(),
            b.state().temperature_k().to_bits()
        );
        assert_eq!(a.accepted_increment_j_kg_k(), b.accepted_increment_j_kg_k());
        assert_eq!(a.enthalpy_change_j_kg(), b.enthalpy_change_j_kg());
        assert_eq!(
            a.internal_energy_change_j_kg(),
            b.internal_energy_change_j_kg()
        );
    }
}

#[test]
fn g0_ps_rest_is_exact_after_one_original_pt_selection() {
    let source = model([1e22, -3e21], [2e18, -7e17], [100.0, 450.0]);
    let mut control = controls();
    control.maximum_pt_evaluations = 1;
    for t in [180.0, 220.0, 280.0] {
        let initial = pt(&source, t, 3e6);
        let solved = increment(&source, t, 3e6, 3e6, 0.0, control);
        assert_eq!(solved.state(), &initial);
        assert_eq!(solved.accepted_increment_j_kg_k(), Some(0.0));
        assert_eq!(solved.enthalpy_change_j_kg(), Some(0.0));
        assert_eq!(solved.internal_energy_change_j_kg(), Some(0.0));
        assert_eq!(solved.residual_j_kg_k(), 0.0);
        assert_eq!(solved.pt_evaluations(), 1);
        assert_eq!(
            solved.eos_evaluation_upper_bound(),
            initial.phase().eos_evaluation_upper_bound()
        );
    }
}

#[test]
fn g4_ps_retains_source_phase_bracket_and_whole_call_work_refusals() {
    let source = source();
    let mut capped = controls();
    capped.maximum_pt_evaluations = 2;
    assert!(matches!(
        source.phase_equilibrium_pressure_entropy_increment(
            Temperature::new(180.0),
            Pressure::new(3e6),
            Pressure::new(3e6),
            917.756_507_280_260,
            capped,
            || false,
        ),
        Err(E::WorkLimit {
            flash_evaluations: 2
        })
    ));
    let run = |control| {
        source.phase_equilibrium_pressure_entropy_increment(
            Temperature::new(220.0),
            Pressure::new(3e6),
            Pressure::new(2e6),
            0.0,
            control,
            || false,
        )
    };
    let mut unresolved = controls();
    unresolved.stability.maximum_iterations = 1;
    assert!(matches!(run(unresolved), Err(E::Phase(_))));
    unresolved = controls();
    unresolved.flash_options.maximum_iterations = 1;
    assert!(matches!(run(unresolved), Err(E::Phase(_))));
    assert!(matches!(
        source.phase_equilibrium_pressure_entropy(Pressure::new(3e6), 1e6, controls(), || false),
        Err(E::OutsideEntropyBracket)
    ));
    let limited = model([0.0; 2], [0.0; 2], [200.0, 260.0]);
    assert!(matches!(
        limited.phase_equilibrium_pressure_entropy(
            Pressure::new(3e6),
            5507.749_287_855_677,
            controls(),
            || false,
        ),
        Err(E::Material(_))
    ));
    for invalid in [f64::NAN, f64::INFINITY] {
        assert!(matches!(
            source.phase_equilibrium_pressure_entropy(
                Pressure::new(3e6),
                invalid,
                controls(),
                || false,
            ),
            Err(E::Invalid(_))
        ));
    }
    let mut invalid = controls();
    invalid.maximum_pt_evaluations = 129;
    assert!(matches!(run(invalid), Err(E::Invalid(_))));
}

#[test]
fn g4_ps_final_cancellation_publishes_no_state_and_replay_is_deterministic() {
    let source = source();
    let mut checkpoints = 0;
    let solved = source
        .phase_equilibrium_pressure_entropy(
            Pressure::new(3e6),
            5507.749_287_855_677,
            controls(),
            || {
                checkpoints += 1;
                false
            },
        )
        .unwrap();
    let mut seen = 0;
    assert_eq!(
        source
            .phase_equilibrium_pressure_entropy(
                Pressure::new(3e6),
                5507.749_287_855_677,
                controls(),
                || {
                    seen += 1;
                    seen == checkpoints
                },
            )
            .unwrap_err(),
        E::Cancelled
    );
    assert_eq!(seen, checkpoints);
    let replay = source
        .phase_equilibrium_pressure_entropy(
            Pressure::new(3e6),
            5507.749_287_855_677,
            controls(),
            || false,
        )
        .unwrap();
    assert_eq!(solved, replay);
}
