//! Same original bodies, fixed-source bond and clock through partial damage.
use super::*;
use fs_couple::cohesive_contact::CohesiveContactConfig;
use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_couple::render::schedule::force::coupled::contact::cohesive::{
    CohesiveModalConfig, CohesiveModalError, CohesiveModalSystem,
};
use fs_couple::render::schedule::force::coupled::{
    CoupledModalSystem, ModalAttachment, ModalConnection, ModalCouplingConfig, ModalCouplingError,
};
use fs_exec::CancelGate;
use fs_math::{c64::C64, det};

const MASSES: [f64; 2] = [2.0, 3.0];

fn native(source: &ResolvedBilinearCohesive, ceiling: f64) -> CohesiveModalSystem<'_> {
    native_with_connection(source, ceiling, true)
}

fn native_with_connection(
    source: &ResolvedBilinearCohesive,
    ceiling: f64,
    connected: bool,
) -> CohesiveModalSystem<'_> {
    native_at_rate(source, ceiling, connected, 5000)
}

fn native_at_rate(
    source: &ResolvedBilinearCohesive,
    ceiling: f64,
    connected: bool,
    sample_rate: u32,
) -> CohesiveModalSystem<'_> {
    let models = MASSES
        .map(|mass| {
            let mut model = ModalAcousticTimeModel::try_new_with_free_coordinates(
                sample_rate,
                vec![
                    ModalAcousticMode {
                        angular_frequency_rad_s: 0.0,
                        damping_ratio: 0.0,
                        pressure_per_modal_velocity: C64::ZERO,
                    };
                    3
                ],
                ModalAcousticTimeBudget::audible_reference(),
            )
            .unwrap();
            model
                .restore_states(&[0.012, -0.006, 0.003].map(|velocity| ModalAcousticState {
                    displacement_m_sqrt_kg: 0.0,
                    velocity_m_sqrt_kg_per_s: det::sqrt(mass) * velocity,
                }))
                .unwrap();
            model
        })
        .into_iter()
        .collect();
    let shape = [0.7, 0.3, -0.2];
    let network = CoupledModalSystem::new(
        models,
        if connected {
            vec![ModalConnection {
                left: ModalAttachment {
                    component: 0,
                    shapes: shape.map(|x| x / det::sqrt(MASSES[0])).to_vec(),
                },
                right: ModalAttachment {
                    component: 1,
                    shapes: shape.map(|x| x / det::sqrt(MASSES[1])).to_vec(),
                },
                stiffness_n_m: 1e8,
                damping_n_s_m: 0.0,
                rest_extension_m: 0.0,
            }]
        } else {
            vec![]
        },
        ModalCouplingConfig {
            max_modes: 6,
            max_connections: 2,
            max_setup_terms: 4096,
            nyquist_guard_fraction: 0.9,
            maximum_total_energy_j: ceiling,
            maximum_abs_pressure_pa: 1000.0,
            maximum_abs_connection_force_n: 1e5,
            solve_relative_tolerance: 1e-11,
            energy_absolute_tolerance_j: 1e-12,
            energy_relative_tolerance: 1e-10,
        },
        &CancelGate::new(),
    )
    .unwrap();
    let columns = core::array::from_fn(|i| {
        let mut column = vec![0.0; 6];
        column[i] = -1.0 / det::sqrt(MASSES[0]);
        column[3 + i] = 1.0 / det::sqrt(MASSES[1]);
        column
    });
    CohesiveModalSystem::from_source(
        network,
        source,
        receiving(source),
        AREA,
        columns,
        CohesiveModalConfig {
            contact: CohesiveContactConfig {
                max_iterations: 48,
                max_material_evaluations: 600,
                max_backtracks: 12,
                jump_absolute_tolerance_m: 1e-14,
                jump_relative_tolerance: 1e-11,
                maximum_jump_increment_m: 0.01,
                maximum_force_n: 100.0,
                material_energy_tolerance_j_m2: 1e-8,
            },
            force_absolute_tolerance_n: 1e-9,
            force_relative_tolerance: 1e-9,
        },
        &CancelGate::new(),
    )
    .unwrap()
}

fn states(system: &CohesiveModalSystem<'_>) -> Vec<ModalAcousticState> {
    system
        .components()
        .iter()
        .flat_map(|model| model.states().iter().copied())
        .collect()
}
fn momentum(system: &CohesiveModalSystem<'_>) -> [f64; 3] {
    core::array::from_fn(|i| {
        (0..2)
            .map(|side| {
                det::sqrt(MASSES[side])
                    * system.components()[side].states()[i].velocity_m_sqrt_kg_per_s
            })
            .sum()
    })
}

fn enabled<'a>(
    source: &'a ResolvedBilinearCohesive,
    friction: &'a ResidualCoulombSource,
    ceiling: f64,
) -> CohesiveModalSystem<'a> {
    let mut system = native(source, ceiling);
    let before = (
        states(&system),
        system.point(),
        system.samples_rendered(),
        system.last_frame().copied(),
    );
    system
        .enable_frictional_damage(
            friction,
            receiving(source),
            ReclosureActivation::AtRest,
            config(),
            &CancelGate::new(),
        )
        .unwrap();
    assert_eq!(
        (
            states(&system),
            system.point(),
            system.samples_rendered(),
            system.last_frame().copied()
        ),
        before
    );
    assert!(core::ptr::eq(system.source().unwrap(), source));
    assert!(core::ptr::eq(
        system.damage_friction_source().unwrap(),
        friction
    ));
    system
}

#[test]
fn g1_native_progressive_damage_sliding_opening_failure_and_reclosure_close_actual_work() {
    let (source, friction) = sources();
    let mut system = enabled(&source, &friction, 1.0);
    let initial = system.total_energy_j().unwrap();
    let initial_momentum = momentum(&system);
    let s = system.compliance_m_per_n();
    assert!(s[0][1].abs() > 1e-10 && s[0][2].abs() > 1e-10 && s[1][2].abs() > 1e-10);
    let mut work = 0.0;
    let mut loss = 0.0;
    let mut partial_sliding = false;
    let mut damage_growth_sliding = false;
    let mut opened = false;
    let mut reclosed = false;
    let mut failed_history = None;
    for drive in [
        [-1.0, 0.0, 0.0],
        [0.0, 1.6, 0.8],
        [0.0, -0.12, 0.03],
        [1.6, 0.04, 0.0],
        [0.8, 0.07, -0.04],
        [-2.2, 0.18, 0.12],
        [0.0, 5.0, 2.0],
        [3.0, 0.2, 0.1],
        [1.0, 0.2, 0.1],
        [-4.0, -0.2, 0.1],
    ] {
        let old = system.point();
        let frame = *system.step(&[0.0; 6], drive).unwrap();
        let extra = frame.damage_friction.unwrap();
        let (original, mean, weighted) = old
            .trial_over_time_with_damage_weighted_compression(
                frame.jump_m,
                system.sample_period_s(),
                1e-8,
            )
            .unwrap();
        assert_eq!(system.point(), original.equilibrium.point);
        assert!(frame.damage >= old.damage());
        close(extra.damaged_compression_force_n, AREA * weighted, 1e-12);
        let slip = (frame.jump_m[1] - old.jump_m()[1]).hypot(frame.jump_m[2] - old.jump_m()[2]);
        if extra.regime == DamageFrictionRegime::Sliding {
            let cap = 0.3 * AREA * weighted;
            close(extra.force_n[0].hypot(extra.force_n[1]), cap, 1e-11);
            close(extra.dissipation_j, cap * slip, 1e-13);
            for i in 1..3 {
                close(
                    frame.mean_force_n[i],
                    AREA * mean[i] + cap * (frame.jump_m[i] - old.jump_m()[i]) / slip,
                    1e-9,
                );
            }
            partial_sliding |=
                frame.damage > 0.0 && frame.damage < 1.0 && extra.dissipation_j > 0.0;
            damage_growth_sliding |=
                frame.damage > old.damage() && frame.damage < 1.0 && extra.dissipation_j > 0.0;
        }
        assert_eq!(frame.residual_friction_regime, None);
        assert_eq!(frame.residual_friction_dissipation_j, 0.0);
        assert!(extra.work_defect_j.abs() < 1e-12);
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        assert!(frame.material_evaluations <= config().contact.max_material_evaluations);
        work += frame.external_work_j + frame.prescribed_work_j;
        loss += frame.network_dissipation_j
            + frame.fracture_dissipation_j
            + frame.viscous_dissipation_j
            + extra.dissipation_j;
        close(
            system.total_energy_j().unwrap() + loss,
            initial + work,
            3e-12,
        );
        for i in 0..3 {
            close(momentum(&system)[i], initial_momentum[i], 5e-15);
        }
        if old.jump_m()[0] >= 0.0 && frame.jump_m[0] >= 0.0 {
            assert_eq!(extra.dissipation_j, 0.0);
            assert_eq!(extra.regime, DamageFrictionRegime::Inactive);
            opened = true;
        }
        reclosed |= opened && frame.jump_m[0] < 0.0 && extra.dissipation_j > 0.0;
        if let Some(history) = failed_history {
            assert!(system.point().failed());
            close(system.point().dissipated_j_m2(), history, 1e-10);
            assert_eq!(frame.fracture_dissipation_j, 0.0);
        } else if system.point().failed() {
            failed_history = Some(system.point().dissipated_j_m2());
        }
    }
    assert!(partial_sliding && damage_growth_sliding && opened && reclosed);
    assert!(failed_history.is_some());
    assert_eq!(system.samples_rendered(), 10);
}

#[test]
fn g4_native_source_cancellation_and_late_energy_refusal_leave_retry_identical() {
    let (source, friction) = sources();
    let mut system = native(&source, 0.001);
    let snapshot = (
        states(&system),
        system.point(),
        system.last_frame().copied(),
        system.samples_rendered(),
    );
    let materials = receiving(&source);
    assert!(
        system
            .enable_frictional_damage(
                &friction,
                [materials[1], materials[0]],
                ReclosureActivation::AtRest,
                config(),
                &CancelGate::new()
            )
            .is_err()
    );
    assert_eq!(
        (
            states(&system),
            system.point(),
            system.last_frame().copied(),
            system.samples_rendered()
        ),
        snapshot
    );
    let cancelled = CancelGate::new();
    cancelled.request();
    assert!(
        system
            .enable_frictional_damage(
                &friction,
                materials,
                ReclosureActivation::AtRest,
                config(),
                &cancelled
            )
            .is_err()
    );
    assert!(system.damage_friction_source().is_none());
    system
        .enable_frictional_damage(
            &friction,
            materials,
            ReclosureActivation::AtRest,
            config(),
            &CancelGate::new(),
        )
        .unwrap();
    let mut replay = enabled(&source, &friction, 0.001);
    for owner in [&mut system, &mut replay] {
        owner.step(&[0.0; 6], [-1.0, 0.0, 0.0]).unwrap();
        owner.step(&[0.0; 6], [0.0, 1.6, 0.8]).unwrap();
    }
    assert!(system.point().damage() > 0.0 && system.point().damage() < 1.0);
    let snapshot = (
        states(&system),
        system.point(),
        system.last_frame().copied(),
        system.samples_rendered(),
        system.damage_friction_occupancy(),
    );
    assert!(
        system
            .step_under_gate(&[0.0; 6], [0.0, 0.1, 0.0], &cancelled)
            .is_err()
    );
    let refused = system.step(&[0.0; 6], [-10.0, 0.4, 0.0]);
    assert!(
        matches!(
            refused,
            Err(CohesiveModalError::Network(ModalCouplingError::Budget {
                what: "total energy including damaged cohesive interface",
                ..
            }))
        ),
        "{refused:?}"
    );
    assert_eq!(
        (
            states(&system),
            system.point(),
            system.last_frame().copied(),
            system.samples_rendered(),
            system.damage_friction_occupancy()
        ),
        snapshot
    );
    assert_eq!(
        *system.step(&[0.0; 6], [0.0, -0.12, 0.03]).unwrap(),
        *replay.step(&[0.0; 6], [0.0, -0.12, 0.03]).unwrap()
    );
    assert_eq!(states(&system), states(&replay));
    assert_eq!(system.point(), replay.point());
    assert_eq!(
        system.damage_friction_occupancy(),
        replay.damage_friction_occupancy()
    );
}

#[test]
fn g1_partial_bond_zero_net_slip_retains_actual_moving_endpoint_occupancy() {
    let (source, friction) = sources();
    let mut system = native_with_connection(&source, 1.0, false);
    system.step(&[0.0; 6], [2.0, 0.0, 0.0]).unwrap();
    let dt = system.sample_period_s();
    let close_rate = (-1e-4 - system.point().jump_m()[0]) / dt;
    system.step(&[0.0; 6], [close_rate, 0.0, 0.0]).unwrap();
    assert!(system.point().damage() > 0.0 && system.point().damage() < 1.0);
    system
        .enable_frictional_damage(
            &friction,
            receiving(&source),
            ReclosureActivation::AtRest,
            config(),
            &CancelGate::new(),
        )
        .unwrap();
    let velocity = |s: &CohesiveModalSystem<'_>, i: usize| {
        s.components()[1].states()[i].velocity_m_sqrt_kg_per_s / det::sqrt(MASSES[1])
            - s.components()[0].states()[i].velocity_m_sqrt_kg_per_s / det::sqrt(MASSES[0])
    };
    let old = system.point().jump_m();
    let compliance = system.compliance_m_per_n();
    let pressure = -AREA * COMPRESSION * old[0];
    let damaged = system.point().damage() * pressure;
    let reaction = 0.2 * damaged;
    let drive = [
        -velocity(&system, 0) - compliance[0][0] * pressure / dt,
        -velocity(&system, 1) + compliance[1][1] * reaction / dt,
        -velocity(&system, 2),
    ];
    let frame = *system.step(&[0.0; 6], drive).unwrap();
    assert_eq!(
        frame.damage_friction.unwrap().regime,
        DamageFrictionRegime::Sticking
    );
    assert_eq!(frame.damage_friction.unwrap().dissipation_j, 0.0);
    close(frame.jump_m[1], old[1], 1e-14);
    close(frame.mean_force_n[1], reaction, 1e-10);
    assert_ne!(velocity(&system, 1) + drive[1], 0.0);
    assert_eq!(
        system.damage_friction_occupancy(),
        Some(DamageFrictionRegime::Sliding)
    );
    let drive = [
        -velocity(&system, 0) - compliance[0][0] * pressure / dt,
        -velocity(&system, 1) + compliance[1][1] * 0.45 * damaged / dt,
        -velocity(&system, 2),
    ];
    let frame = *system.step(&[0.0; 6], drive).unwrap();
    let friction = frame.damage_friction.unwrap();
    assert_eq!(friction.regime, DamageFrictionRegime::Sliding);
    close(friction.force_n[0], 0.3 * damaged, 2e-10);
    assert!(friction.dissipation_j > 0.0);
    assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
}

#[test]
fn g2_refining_the_actual_mechanical_timestep_converges_sliding_work_and_motion() {
    let (source, friction) = sources();
    let mut results = Vec::new();
    for rate in [5000, 10_000, 20_000, 40_000] {
        let mut system = native_at_rate(&source, 1.0, true, rate);
        system
            .enable_frictional_damage(
                &friction,
                receiving(&source),
                ReclosureActivation::AtRest,
                config(),
                &CancelGate::new(),
            )
            .unwrap();
        let mut friction_work = 0.0;
        for _ in 0..rate / 250 {
            let frame = *system.step(&[0.0; 6], [-0.08, 0.4, 0.18]).unwrap();
            friction_work += frame.damage_friction.unwrap().dissipation_j;
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        }
        assert!(system.point().failed());
        assert!(friction_work > 0.0);
        results.push((friction_work, system.point().jump_m(), states(&system)));
    }
    let reference = &results[3];
    let errors: Vec<f64> = results[..3]
        .iter()
        .map(|r| {
            let work = (r.0 - reference.0).abs() / reference.0;
            let motion =
                r.2.iter()
                    .zip(&reference.2)
                    .map(|(a, b)| (a.velocity_m_sqrt_kg_per_s - b.velocity_m_sqrt_kg_per_s).powi(2))
                    .sum::<f64>()
                    .sqrt();
            work + motion
        })
        .collect();
    assert!(
        errors[1] < 0.55 * errors[0] && errors[2] < 0.55 * errors[1],
        "{errors:?}"
    );
    assert!(errors[2] < 1e-4, "{errors:?}");
}
