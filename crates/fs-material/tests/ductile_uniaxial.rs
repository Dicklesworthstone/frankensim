//! G0/G1/G3 native axial ductility, independent scalar work and tangent oracles.
//! Coefficients are synthetic verification data, not an alloy calibration.
#![deny(unsafe_code)]

use fs_material::ductile::evolution::{
    DuctileEvolutionError, DuctileJ2Conditions, DuctileJ2Settings, DuctileWorkEvolution,
};
use fs_material::ductile::uniaxial::{UniaxialDuctileJ2, UniaxialDuctileJ2State};
use fs_material::ductile::{
    DuctileDomain, JohnsonCookInitiation, JohnsonCookParameters, TabulatedDuctileLocus,
};
use fs_material::elastic::IsotropicElastic;
use fs_material::fiber::Uniaxial;
use fs_material::plastic::J2Plasticity;
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::reduction::UniaxialJ2;

fn native(hardening: f64, kinematic: f64) -> UniaxialJ2 {
    UniaxialJ2::new(
        J2Plasticity::new_mixed(
            IsotropicElastic::new(2000.0, 0.3, 0.05).unwrap(),
            2.0,
            hardening,
            kinematic,
        )
        .unwrap(),
    )
    .unwrap()
}
fn domain() -> DuctileDomain {
    DuctileDomain {
        triaxiality: [-1.0, 1.0],
        strain_rate_s_inv: [1e-10, 1e3],
        temperature_k: [290.0, 310.0],
    }
}
fn locus(strain: f64) -> TabulatedDuctileLocus {
    TabulatedDuctileLocus::new(vec![(-1.0, strain), (1.0, strain)], domain()).unwrap()
}
fn settings(length: f64, work: f64) -> DuctileJ2Settings {
    DuctileJ2Settings {
        characteristic_length_m: length,
        post_initiation_plastic_work_j_m2: work,
        maximum_damage_increment: 1.0,
        energy_tolerances: EnergyTolerances {
            energy_j_m3: 1e-10,
            yield_pa: 1e-9,
        },
    }
}
fn conditions() -> DuctileJ2Conditions {
    DuctileJ2Conditions {
        duration_s: 1.0,
        temperature_k: 300.0,
    }
}
fn strain_at_alpha(native: &UniaxialJ2, alpha: f64) -> f64 {
    let b = native.backbone();
    alpha + (b.yield_stress + (b.hardening + b.kinematic_hardening) * alpha) / b.elastic.youngs
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a} != {b}; difference {}",
        a - b
    );
}
fn unchanged_physical(a: &UniaxialDuctileJ2State, b: &UniaxialDuctileJ2State) {
    assert_eq!(a.strain(), b.strain());
    assert_eq!(a.stress_pa(), b.stress_pa());
    assert_eq!(a.native_state(), b.native_state());
    assert_eq!(a.initiation(), b.initiation());
    assert_eq!(a.integrity(), b.integrity());
    assert_eq!(a.stored_energy_j_m3(), b.stored_energy_j_m3());
    assert_eq!(a.dissipated_plastic_j_m3(), b.dissipated_plastic_j_m3());
    assert_eq!(a.dissipated_damage_j_m3(), b.dissipated_damage_j_m3());
    assert_eq!(
        a.numerical_dissipation_j_m3(),
        b.numerical_dissipation_j_m3()
    );
}

#[test]
fn intact_axial_response_preserves_the_original_native_return_and_energy() {
    let model =
        UniaxialDuctileJ2::new(native(50.0, 30.0), locus(1.0), settings(0.02, 0.0002)).unwrap();
    let mut state = model.initial_state();
    for strain in [0.0002, 0.006, 0.005, -0.004, 0.008] {
        let expected = model.native().response(strain, state.native_state());
        let step = model.trial(strain, &state, conditions()).unwrap();
        assert_eq!(step.stress_pa, expected.stress_pa);
        assert_eq!(step.tangent_pa, expected.tangent_pa);
        assert_eq!(step.state.native_state(), &expected.state);
        assert_eq!(step.state.integrity(), 1.0);
        assert_eq!(step.energy.damage_release_j_m3, 0.0);
        assert_eq!(step.state.last_step_energy(), &step.energy);
        close(step.energy.residual_j_m3, 0.0, 1e-12);
        state = step.state;
    }
}

#[test]
fn native_onset_integral_gives_softening_and_the_exact_axial_chain_rule() {
    let model =
        UniaxialDuctileJ2::new(native(50.0, 30.0), locus(0.002), settings(0.02, 0.0002)).unwrap();
    let old = model.initial_state();
    let alpha = 0.008;
    let strain = strain_at_alpha(model.native(), alpha);
    let step = model.trial(strain, &old, conditions()).unwrap();
    let work = 0.02 * (2.0 * (alpha - 0.002) + 25.0 * (alpha * alpha - 0.002_f64.powi(2)));
    let g = (-work / 0.0002).exp();
    let sigma = 2.0 + 80.0 * alpha;
    let da_de = 2000.0 / 2080.0;
    let tangent = g * (2000.0 * 80.0 / 2080.0 - sigma * 100.0 * (2.0 + 50.0 * alpha) * da_de);
    close(step.state.equivalent_plastic_strain(), alpha, 1e-16);
    close(step.state.equivalent_plastic_work_j_m2(), work, 1e-17);
    close(step.state.integrity(), g, 1e-13);
    close(step.stress_pa, g * sigma, 1e-12);
    close(step.tangent_pa, tangent, 1e-9);
    assert!(step.tangent_pa < 0.0);
    assert_eq!(old.initiation().utilization(), 0.0);
    let h = 1e-8;
    let finite_difference = (model
        .trial(strain + h, &old, conditions())
        .unwrap()
        .stress_pa
        - model
            .trial(strain - h, &old, conditions())
            .unwrap()
            .stress_pa)
        / (2.0 * h);
    close(step.tangent_pa, finite_difference, 1e-5);
}

#[test]
fn moving_rate_dependent_onset_uses_actual_interval_and_consistent_tangent() {
    let initiation = JohnsonCookInitiation::new(
        JohnsonCookParameters {
            d1: 0.0018,
            d2: 0.0004,
            d3: 0.5,
            d4: 0.04,
            d5: 0.0,
            reference_strain_rate_s_inv: 0.001,
            transition_temperature_k: 290.0,
            melting_temperature_k: 1200.0,
        },
        domain(),
    )
    .unwrap();
    let model =
        UniaxialDuctileJ2::new(native(50.0, 30.0), initiation, settings(0.02, 0.0002)).unwrap();
    let old = model.initial_state();
    let strain = strain_at_alpha(model.native(), 0.008);
    let prepared = model.prepared(conditions()).unwrap();
    let fast = prepared
        .with_conditions(DuctileJ2Conditions {
            duration_s: 0.1,
            ..conditions()
        })
        .unwrap();
    assert!(prepared.model().same_model(fast.model()));
    let slow = prepared.response(strain, &old).unwrap();
    let rapid = fast.response(strain, &old).unwrap();
    assert_eq!(slow.state.native_state(), rapid.state.native_state());
    assert!(rapid.stress_pa > slow.stress_pa);
    assert!(rapid.state.damage() < slow.state.damage());
    for law in [&prepared, &fast] {
        let step = law.response(strain, &old).unwrap();
        let h = 1e-8;
        let fd = (law.response(strain + h, &old).unwrap().stress_pa
            - law.response(strain - h, &old).unwrap().stress_pa)
            / (2.0 * h);
        close(step.tangent_pa, fd, 1e-5);
        assert_eq!(law.stress(strain, &old), step.stress_pa);
        assert_eq!(law.tangent(strain, &old), step.tangent_pa);
        assert_eq!(law.update_state(strain, &old), step.state);
    }
}

#[test]
fn reverse_flow_spends_more_ductility_and_all_work_accounts_close() {
    let model =
        UniaxialDuctileJ2::new(native(50.0, 30.0), locus(0.002), settings(0.01, 0.0005)).unwrap();
    let mut state = model.initial_state();
    let mut work = 0.0;
    let mut damages = Vec::new();
    for strain in [0.006, 0.0058, -0.007, 0.009] {
        let before = state.clone();
        let step = model.trial(strain, &state, conditions()).unwrap();
        let h = 1e-8;
        let fd = (model
            .trial(strain + h, &state, conditions())
            .unwrap()
            .stress_pa
            - model
                .trial(strain - h, &state, conditions())
                .unwrap()
                .stress_pa)
            / (2.0 * h);
        close(step.tangent_pa, fd, 1e-5);
        assert_eq!(state, before);
        assert!(step.state.equivalent_plastic_strain() >= state.equivalent_plastic_strain());
        assert!(step.state.damage() >= state.damage());
        assert!(step.energy.plastic_dissipation_j_m3 >= 0.0);
        assert!(step.energy.damage_release_j_m3 >= 0.0);
        assert!(step.energy.numerical_dissipation_j_m3 >= 0.0);
        work += step.energy.mechanical_work_j_m3;
        damages.push(step.state.damage());
        state = step.state;
        let replay = model.trial(strain, &state, conditions()).unwrap();
        unchanged_physical(&replay.state, &state);
        assert_eq!(replay.energy.mechanical_work_j_m3, 0.0);
        assert_eq!(replay.state.last_step_energy().mechanical_work_j_m3, 0.0);
        assert_eq!(replay.state.last_step_energy().damage_release_j_m3, 0.0);
    }
    assert_eq!(damages[0], damages[1]);
    assert!(damages[2] > damages[1]);
    close(
        work,
        state.stored_energy_j_m3()
            + state.dissipated_plastic_j_m3()
            + state.dissipated_damage_j_m3()
            + state.numerical_dissipation_j_m3(),
        1e-11,
    );
}

#[test]
fn regularization_preserves_work_capacity_and_finite_support_reaches_exact_failure() {
    let mut responses = Vec::new();
    for length in [0.01, 0.02] {
        let model =
            UniaxialDuctileJ2::new(native(0.0, 0.0), locus(0.002), settings(length, 0.0002))
                .unwrap();
        let alpha = 0.002 + 0.0001 / length;
        let step = model
            .trial(
                strain_at_alpha(model.native(), alpha),
                &model.initial_state(),
                conditions(),
            )
            .unwrap();
        close(step.state.equivalent_plastic_work_j_m2(), 0.0002, 1e-17);
        responses.push((
            step.stress_pa,
            step.state.degraded_equivalent_plastic_work_j_m2(),
        ));
    }
    close(responses[0].0, responses[1].0, 1e-12);
    close(responses[0].1, responses[1].1, 1e-17);

    let model = UniaxialDuctileJ2::new(native(0.0, 0.0), locus(0.002), settings(0.02, 0.0001))
        .unwrap()
        .with_evolution(DuctileWorkEvolution::FiniteSupportQuadratic);
    let first = model
        .trial(
            strain_at_alpha(model.native(), 0.004),
            &model.initial_state(),
            conditions(),
        )
        .unwrap();
    let failed = model
        .trial(
            strain_at_alpha(model.native(), 0.015),
            &first.state,
            conditions(),
        )
        .unwrap();
    assert!(failed.state.is_fully_failed());
    assert_eq!(failed.state.integrity(), 0.0);
    assert_eq!(failed.stress_pa, 0.0);
    assert_eq!(failed.tangent_pa, 0.0);
    assert_eq!(failed.state.degraded_equivalent_plastic_work_j_m2(), 0.0001);
    close(
        failed.energy.damage_release_j_m3,
        first.state.stored_energy_j_m3(),
        1e-12,
    );
    let compression = model.trial(-0.01, &failed.state, conditions()).unwrap();
    assert_eq!(compression.stress_pa, 0.0);
    assert_eq!(compression.tangent_pa, 0.0);
    assert!(
        compression.state.equivalent_plastic_strain() > failed.state.equivalent_plastic_strain()
    );
    assert!(compression.state.is_fully_failed());
    assert_eq!(compression.energy.plastic_dissipation_j_m3, 0.0);
    assert_eq!(
        compression.state.degraded_equivalent_plastic_work_j_m2(),
        0.0001
    );
    // A finite terminal coordinate can have an unrepresentable inactive
    // derivative. Exact failure must not evaluate zero times that derivative.
    let stiff = UniaxialJ2::new(
        J2Plasticity::new_mixed(
            IsotropicElastic::new(1e6, 0.3, 0.05).unwrap(),
            2000.0,
            0.0,
            0.0,
        )
        .unwrap(),
    )
    .unwrap();
    let terminal = UniaxialDuctileJ2::new(stiff, locus(0.002), settings(0.02, 2e-309))
        .unwrap()
        .with_evolution(DuctileWorkEvolution::FiniteSupportQuadratic);
    let end = terminal
        .trial(
            strain_at_alpha(terminal.native(), 0.008),
            &terminal.initial_state(),
            conditions(),
        )
        .unwrap();
    assert!(end.state.is_fully_failed());
    assert_eq!(end.tangent_pa, 0.0);
}

#[test]
fn owned_restart_preserves_model_and_refusals_never_clip_or_change_history() {
    let model =
        UniaxialDuctileJ2::new(native(50.0, 30.0), locus(0.002), settings(0.02, 0.0002)).unwrap();
    let accepted = model
        .trial(0.005, &model.initial_state(), conditions())
        .unwrap()
        .state;
    let restarted = model.clone();
    assert!(restarted.same_model(&model));
    assert_eq!(
        model.trial(0.007, &accepted, conditions()).unwrap(),
        restarted.trial(0.007, &accepted, conditions()).unwrap()
    );
    let foreign =
        UniaxialDuctileJ2::new(native(50.0, 30.0), locus(0.002), settings(0.02, 0.0002)).unwrap();
    assert!(matches!(
        foreign.trial(0.007, &accepted, conditions()),
        Err(DuctileEvolutionError::StateMismatch)
    ));
    let changed = model
        .clone()
        .with_evolution(DuctileWorkEvolution::FiniteSupportQuadratic);
    assert!(matches!(
        changed.trial(0.007, &accepted, conditions()),
        Err(DuctileEvolutionError::StateMismatch)
    ));
    let before = accepted.clone();
    for strain in [f64::NAN, 0.2] {
        assert!(model.trial(strain, &accepted, conditions()).is_err());
        assert_eq!(accepted, before);
    }
    assert!(
        model
            .prepared(DuctileJ2Conditions {
                duration_s: 0.0,
                ..conditions()
            })
            .is_err()
    );
    assert!(
        model
            .prepared(DuctileJ2Conditions {
                temperature_k: 400.0,
                ..conditions()
            })
            .is_err()
    );
    let limited = UniaxialDuctileJ2::new(
        native(50.0, 30.0),
        locus(0.002),
        DuctileJ2Settings {
            maximum_damage_increment: 0.01,
            ..settings(0.02, 0.0002)
        },
    )
    .unwrap();
    let zero = limited.initial_state();
    assert!(matches!(
        limited.trial(0.008, &zero, conditions()),
        Err(DuctileEvolutionError::DamageIncrement { .. })
    ));
    assert_eq!(zero.damage(), 0.0);
    assert_eq!(zero.equivalent_plastic_strain(), 0.0);
}
