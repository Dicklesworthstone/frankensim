use super::super::{DuctileDomain, JohnsonCookParameters};
use super::*;
use crate::elastic::IsotropicElastic;

fn native(hardening: f64, kinematic: f64) -> J2Plasticity {
    J2Plasticity::new_mixed(
        IsotropicElastic::new(2000.0, 0.3, 0.05).unwrap(),
        2.0,
        hardening,
        kinematic,
    )
    .unwrap()
}
fn domain() -> DuctileDomain {
    DuctileDomain {
        triaxiality: [-5.0, 5.0],
        strain_rate_s_inv: [1e-8, 1e3],
        temperature_k: [290.0, 310.0],
    }
}
fn locus(strain: f64) -> TabulatedDuctileLocus {
    TabulatedDuctileLocus::new(vec![(-5.0, strain), (5.0, strain)], domain()).unwrap()
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

// Independent monotonic J2 references, including Prager's contribution to the
// actual uniaxial/shear yield stress. alpha is accumulated equivalent strain.
fn shear_at_alpha(law: &J2Plasticity, alpha: f64) -> Voigt {
    let mu = law.elastic.youngs / (2.0 * (1.0 + law.elastic.poisson));
    let stress = law.yield_stress + (law.hardening + law.kinematic_hardening) * alpha;
    [
        0.0,
        0.0,
        0.0,
        stress / (2.0 * mu * 3.0_f64.sqrt()) + 0.5 * 3.0_f64.sqrt() * alpha,
        0.0,
        0.0,
    ]
}
fn axial_at_alpha(law: &J2Plasticity, alpha: f64) -> Voigt {
    let elastic =
        (law.yield_stress + (law.hardening + law.kinematic_hardening) * alpha) / law.elastic.youngs;
    [
        elastic + alpha,
        -law.elastic.poisson * elastic - 0.5 * alpha,
        -law.elastic.poisson * elastic - 0.5 * alpha,
        0.0,
        0.0,
        0.0,
    ]
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a} != {b}; delta={}", a - b);
}

#[test]
fn intact_response_reuses_native_j2_exactly() {
    let material = DuctileJ2::new(native(50.0, 25.0), locus(0.02), settings(0.02, 0.0002)).unwrap();
    let old = material.initial_state();
    for strain in [
        [0.0002, -0.0001, 0.0, 0.0001, 0.0, 0.0],
        shear_at_alpha(material.native(), 0.001),
    ] {
        let step = material.trial(&strain, &old, conditions()).unwrap();
        assert_eq!(
            step.stress_pa,
            material.native().stress(&strain, old.native_state())
        );
        assert_eq!(
            step.branch_tangent_pa,
            material.native().tangent(&strain, old.native_state())
        );
        assert_eq!(
            step.state.native_state(),
            &material.native().update_state(&strain, old.native_state())
        );
        assert_eq!(step.state.integrity(), 1.0);
        assert_eq!(step.state.damage(), 0.0);
        assert_eq!(step.energy.damage_release_j_m3, 0.0);
        assert_eq!(old.initiation().utilization(), 0.0);
    }
}

#[test]
fn onset_fraction_and_hardening_integral_control_actual_stress() {
    let material =
        DuctileJ2::new(native(50.0, 30.0), locus(0.002), settings(0.02, 0.0002)).unwrap();
    let alpha = 0.008;
    let strain = shear_at_alpha(material.native(), alpha);
    let step = material
        .trial(&strain, &material.initial_state(), conditions())
        .unwrap();
    let work = 0.02 * (2.0 * (alpha - 0.002) + 0.5 * 50.0 * (alpha * alpha - 0.002 * 0.002));
    let integrity = (-work / 0.0002).exp();
    close(step.state.native_state().alpha, alpha, 1e-16);
    close(step.state.initiation().utilization(), 4.0, 1e-13);
    close(step.state.equivalent_plastic_work_j_m2(), work, 1e-17);
    close(step.state.integrity(), integrity, 1e-13);
    close(
        step.stress_pa[3],
        integrity * (2.0 + 80.0 * alpha) / 3.0_f64.sqrt(),
        1e-12,
    );
    close(
        step.state.degraded_equivalent_plastic_work_j_m2(),
        0.0002 * (1.0 - integrity),
        1e-17,
    );
    // Prager storage is retained in the native energy, while Gp explicitly
    // follows the isotropic resistance integral, not a hidden total-work alias.
    assert!(step.state.stored_energy_j_m3() > 0.0);
    assert!(step.energy.numerical_dissipation_j_m3 > 0.0);
}

#[test]
fn characteristic_length_preserves_post_onset_plastic_displacement_work() {
    let work_scale = 0.00002;
    let plastic_displacement = 0.00005;
    let expected = (-2.0_f64 * plastic_displacement / work_scale).exp();
    let mut stresses = Vec::new();
    for length in [0.01, 0.04] {
        let material =
            DuctileJ2::new(native(0.0, 0.0), locus(0.002), settings(length, work_scale)).unwrap();
        let alpha = 0.002 + plastic_displacement / length;
        let step = material
            .trial(
                &shear_at_alpha(material.native(), alpha),
                &material.initial_state(),
                conditions(),
            )
            .unwrap();
        close(step.state.integrity(), expected, 1e-14);
        close(
            step.state.degraded_equivalent_plastic_work_j_m2(),
            work_scale * (1.0 - expected),
            1e-17,
        );
        assert!(
            step.state.integrity() > 0.0,
            "no hidden complete-failure threshold"
        );
        stresses.push(step.stress_pa[3]);
    }
    close(stresses[0], stresses[1], 1e-13);
    // Independent trapezoidal traction/plastic-displacement work verifies the
    // stated area measure; it does not claim total fracture-energy objectivity.
    let count = 10_000;
    let dx = plastic_displacement / count as f64;
    let mut integrated = 0.0;
    for i in 0..count {
        let a = i as f64 * dx;
        let b = (i + 1) as f64 * dx;
        integrated +=
            0.5 * dx * (2.0 * (-2.0 * a / work_scale).exp() + 2.0 * (-2.0 * b / work_scale).exp());
    }
    close(integrated, work_scale * (1.0 - expected), 1e-12);
}

fn check_tangent<L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq>(
    material: &DuctileJ2<L>,
    old: &DuctileJ2State<'_, L>,
    strain: Voigt,
    conditions: DuctileJ2Conditions,
) {
    let base = material.trial(&strain, old, conditions).unwrap();
    let h = 1e-8;
    for j in 0..6 {
        let mut plus = strain;
        let mut minus = strain;
        plus[j] += h;
        minus[j] -= h;
        let a = material.trial(&plus, old, conditions).unwrap().stress_pa;
        let b = material.trial(&minus, old, conditions).unwrap().stress_pa;
        for i in 0..6 {
            let difference = (a[i] - b[i]) / (2.0 * h);
            let expected = base.branch_tangent_pa[i][j];
            close(difference, expected, 3e-6 * (1.0 + expected.abs()));
        }
    }
}

#[test]
fn analytic_tangent_includes_moving_onset_rate_and_compressive_reclosure() {
    let initiation = JohnsonCookInitiation::new(
        JohnsonCookParameters {
            d1: 0.0008,
            d2: 0.0008,
            d3: 0.8,
            d4: 0.1,
            d5: 0.1,
            reference_strain_rate_s_inv: 0.004,
            transition_temperature_k: 250.0,
            melting_temperature_k: 1000.0,
        },
        DuctileDomain {
            strain_rate_s_inv: [1e-6, 1e3],
            ..domain()
        },
    )
    .unwrap();
    let material = DuctileJ2::new(native(50.0, 20.0), initiation, settings(0.02, 0.0005)).unwrap();
    let old = material.initial_state();
    let first = [0.006, -0.001, -0.001, 0.002, 0.0004, -0.0006];
    let physical = DuctileJ2Conditions {
        duration_s: 0.25,
        ..conditions()
    };
    let initiated = material.trial(&first, &old, physical).unwrap();
    assert!(initiated.state.initiation().initiated());
    assert!(initiated.state.damage() > 0.0);
    check_tangent(&material, &old, first, physical);
    let compressed = [-0.006, -0.003, -0.004, 0.008, 0.002, -0.001];
    let later = material
        .trial(&compressed, &initiated.state, physical)
        .unwrap();
    assert!(later.state.native_state().alpha > initiated.state.native_state().alpha);
    assert_eq!(later.state.initiation(), initiated.state.initiation());
    check_tangent(&material, &initiated.state, compressed, physical);
}

#[test]
fn compression_retains_bulk_stiffness_without_healing_or_requerying_onset() {
    let initiation = TabulatedDuctileLocus::new(
        vec![(0.2, 0.001), (0.5, 0.001)],
        DuctileDomain {
            triaxiality: [0.2, 0.5],
            ..domain()
        },
    )
    .unwrap();
    let material = DuctileJ2::new(native(50.0, 0.0), initiation, settings(0.02, 0.0002)).unwrap();
    let started = material
        .trial(
            &axial_at_alpha(material.native(), 0.004),
            &material.initial_state(),
            conditions(),
        )
        .unwrap();
    let mut compression = started.state.native_state().plastic_strain;
    for x in &mut compression[..3] {
        *x -= 0.001;
    }
    let closed = material
        .trial(&compression, &started.state, conditions())
        .unwrap();
    let bulk = 2000.0 / (3.0 * (1.0 - 2.0 * 0.3));
    for i in 0..3 {
        close(closed.stress_pa[i], -0.003 * bulk, 1e-12);
    }
    assert_eq!(closed.state.damage(), started.state.damage());
    assert_eq!(closed.state.initiation(), started.state.initiation());
    assert_eq!(
        closed.state.native_state().alpha,
        started.state.native_state().alpha
    );
    // The actual mean stress is now compression and outside the onset locus;
    // the completed initiation event is retained, not illegally reevaluated.
    let mut tension = started.state.native_state().plastic_strain;
    for x in &mut tension[..3] {
        *x += 0.001;
    }
    let reopened = material
        .trial(&tension, &closed.state, conditions())
        .unwrap();
    for i in 0..3 {
        close(
            reopened.stress_pa[i],
            started.state.integrity() * 0.003 * bulk,
            1e-12,
        );
    }
    assert_eq!(reopened.state.damage(), started.state.damage());
    assert_eq!(reopened.energy.damage_release_j_m3, 0.0);
    assert!(reopened.energy.numerical_dissipation_j_m3 >= 0.0);
}

#[test]
fn endpoint_energy_closes_over_plastic_reversal_and_damage() {
    let material =
        DuctileJ2::new(native(40.0, 30.0), locus(0.0005), settings(0.02, 0.0005)).unwrap();
    let mut state = material.initial_state();
    let mut actual_work = 0.0;
    for strain in [
        [0.004, -0.001, -0.001, 0.004, 0.0, 0.0],
        [0.0; 6],
        [-0.004, -0.002, -0.003, -0.004, 0.001, 0.0],
        [0.003, 0.002, 0.001, 0.006, -0.001, 0.002],
    ] {
        let step = material.trial(&strain, &state, conditions()).unwrap();
        assert_eq!(step, material.trial(&strain, &state, conditions()).unwrap());
        for i in 0..6 {
            actual_work +=
                step.stress_pa[i] * (strain[i] - state.strain()[i]) * if i < 3 { 1.0 } else { 2.0 };
        }
        assert!(step.energy.plastic_dissipation_j_m3 >= 0.0);
        assert!(step.energy.damage_release_j_m3 >= 0.0);
        assert!(step.energy.numerical_dissipation_j_m3 >= 0.0);
        assert!(step.state.damage() >= state.damage());
        state = step.state;
    }
    assert!(state.dissipated_damage_j_m3() > 0.0);
    assert!(state.dissipated_plastic_j_m3() > 0.0);
    close(
        actual_work,
        state.stored_energy_j_m3()
            + state.dissipated_plastic_j_m3()
            + state.dissipated_damage_j_m3()
            + state.numerical_dissipation_j_m3(),
        1e-10,
    );
}

#[test]
fn large_damage_steps_foreign_histories_and_invalid_conditions_refuse_immutably() {
    let mut bounded = settings(0.02, 0.0002);
    bounded.maximum_damage_increment = 0.01;
    let material = DuctileJ2::new(native(0.0, 0.0), locus(0.001), bounded).unwrap();
    let old = material.initial_state();
    let saved = old.clone();
    let strain = shear_at_alpha(material.native(), 0.005);
    assert!(matches!(
        material.trial(&strain, &old, conditions()),
        Err(DuctileEvolutionError::DamageIncrement { .. })
    ));
    assert_eq!(old, saved);
    let foreign = material.clone();
    assert!(matches!(
        foreign.trial(&[0.0; 6], &old, conditions()),
        Err(DuctileEvolutionError::StateMismatch)
    ));
    for physical in [
        DuctileJ2Conditions {
            duration_s: 0.0,
            ..conditions()
        },
        DuctileJ2Conditions {
            temperature_k: 311.0,
            ..conditions()
        },
    ] {
        assert!(material.small_strain_law(physical).is_err());
        assert!(material.trial(&[0.0; 6], &old, physical).is_err());
    }
    let adapter = material.small_strain_law(conditions()).unwrap();
    assert!(adapter.stress(&strain, &old).iter().all(|x| x.is_nan()));
    assert!(!adapter.update_state(&strain, &old).is_valid());
    assert_eq!(old, saved);
}

#[test]
fn constant_locus_refinement_preserves_hardening_work_and_damaged_response() {
    let material =
        DuctileJ2::new(native(50.0, 20.0), locus(0.002), settings(0.02, 0.0002)).unwrap();
    let target = shear_at_alpha(material.native(), 0.008);
    let direct = material
        .trial(&target, &material.initial_state(), conditions())
        .unwrap();
    let mut split = material.initial_state();
    for alpha in [0.001, 0.002, 0.003, 0.005, 0.008] {
        split = material
            .trial(
                &shear_at_alpha(material.native(), alpha),
                &split,
                DuctileJ2Conditions {
                    duration_s: 0.2,
                    ..conditions()
                },
            )
            .unwrap()
            .state;
    }
    close(
        split.native_state().alpha,
        direct.state.native_state().alpha,
        1e-16,
    );
    close(split.integrity(), direct.state.integrity(), 1e-13);
    close(
        split.equivalent_plastic_work_j_m2(),
        direct.state.equivalent_plastic_work_j_m2(),
        1e-16,
    );
    let held = material.trial(&target, &split, conditions()).unwrap();
    for i in 0..6 {
        close(held.stress_pa[i], direct.stress_pa[i], 1e-12);
    }
    // The damage response agrees; backward-Euler mechanical losses need not be
    // step independent or monotone under arbitrary coarse partitions, because
    // endpoint integrity weights them. On resolved paths they approach zero.
    let mut losses = Vec::new();
    for count in [32, 64] {
        let mut refined = material.initial_state();
        for i in 1..=count {
            let strain = target.map(|v| v * (i as f64 / count as f64));
            refined = material
                .trial(
                    &strain,
                    &refined,
                    DuctileJ2Conditions {
                        duration_s: 1.0 / count as f64,
                        ..conditions()
                    },
                )
                .unwrap()
                .state;
        }
        close(refined.integrity(), direct.state.integrity(), 1e-12);
        losses.push(refined.numerical_dissipation_j_m3());
    }
    assert!(losses[1] > 0.0 && losses[1] < 0.65 * losses[0]);
}

#[test]
fn exact_accepted_strain_replay_does_not_invent_pre_onset_plastic_rate() {
    let backbone =
        J2Plasticity::new(IsotropicElastic::new(1200.0, 0.2, 0.05).unwrap(), 3.0, 0.0).unwrap();
    let material = DuctileJ2::new(backbone, locus(0.002), settings(1.0 / 3.0, 0.004)).unwrap();
    // Includes the native tetrahedron fixture gamma=.003 and nearby elastic/
    // plastic accepted endpoints. An equilibrium solver evaluates this exact
    // accepted strain before applying its next prescribed-displacement step.
    for index in 1..=12 {
        let strain = [0.0, 0.0, 0.0, 0.003 * index as f64 / 12.0, 0.0, 0.0];
        let accepted = material
            .trial(&strain, &material.initial_state(), conditions())
            .unwrap();
        assert!(!accepted.state.initiation().initiated());
        let held = material
            .trial(
                &strain,
                &accepted.state,
                DuctileJ2Conditions {
                    duration_s: 0.1,
                    ..conditions()
                },
            )
            .unwrap();
        assert_eq!(held.state, accepted.state);
        for i in 0..6 {
            close(held.stress_pa[i], accepted.stress_pa[i], 1e-13);
        }
        assert_eq!(held.energy.mechanical_work_j_m3, 0.0);
        assert_eq!(held.energy.plastic_dissipation_j_m3, 0.0);
        assert_eq!(held.energy.damage_release_j_m3, 0.0);
        assert_eq!(held.energy.numerical_dissipation_j_m3, 0.0);
    }
    let strain = [0.0, 0.0, 0.0, 0.003, 0.0, 0.0];
    let accepted = material
        .trial(&strain, &material.initial_state(), conditions())
        .unwrap();
    let mut real_loading = strain;
    real_loading[3] += 1e-12;
    assert!(matches!(
        material.trial(
            &real_loading,
            &accepted.state,
            DuctileJ2Conditions {
                duration_s: 0.1,
                ..conditions()
            }
        ),
        Err(DuctileEvolutionError::Initiation(DuctileError::Domain {
            quantity: "plastic_strain_rate_s_inv",
            ..
        }))
    ));
}
