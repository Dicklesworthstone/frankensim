//! In-place handoff on original mechanics, through actual rupture and closure.
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
    let models = MASSES
        .map(|mass| {
            let mut model = ModalAcousticTimeModel::try_new_with_free_coordinates(
                5000,
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
fn rupture(system: &mut CohesiveModalSystem<'_>) {
    let frame = *system.step(&[0.0; 6], [12.0, 0.0, 0.0]).unwrap();
    assert!(system.point().failed());
    assert!(frame.fracture_dissipation_j > 0.0);
    assert_eq!(frame.residual_friction_regime, None);
    assert_eq!(frame.residual_friction_dissipation_j, 0.0);
    assert_eq!(frame.residual_friction_work_defect_j, 0.0);
}

#[test]
fn g1_native_original_body_source_clock_and_fracture_survive_handoff_and_sliding() {
    let (source, friction) = sources();
    let mut system = native(&source, 1.0);
    let initial = system.total_energy_j().unwrap();
    let initial_momentum = momentum(&system);
    let s = system.compliance_m_per_n();
    assert!(s[0][1].abs() > 1e-10 && s[0][2].abs() > 1e-10 && s[1][2].abs() > 1e-10);
    rupture(&mut system);
    let last = system.last_frame().copied();
    let mut accumulated_work = last.unwrap().prescribed_work_j;
    let mut accumulated_loss =
        last.unwrap().fracture_dissipation_j + last.unwrap().network_dissipation_j;
    let before_states = states(&system);
    let before_point = system.point();
    let before_clock = system.samples_rendered();
    let receipts = source.resolved().properties().as_ptr();
    let residual_receipts = friction.resolved().properties().as_ptr();
    system
        .enable_frictional_reclosure(
            &friction,
            receiving(&source),
            ReclosureActivation::AtRest,
            config(),
            &CancelGate::new(),
        )
        .unwrap();
    assert_eq!(states(&system), before_states);
    assert_eq!(system.point(), before_point);
    assert_eq!(system.samples_rendered(), before_clock);
    assert_eq!(system.last_frame().copied(), last);
    assert!(core::ptr::eq(system.source().unwrap(), &source));
    assert!(core::ptr::eq(
        system.residual_friction_source().unwrap(),
        &friction
    ));
    assert_eq!(
        system.source().unwrap().resolved().properties().as_ptr(),
        receipts
    );
    assert_eq!(
        system
            .residual_friction_source()
            .unwrap()
            .resolved()
            .properties()
            .as_ptr(),
        residual_receipts
    );
    let mut friction_work = 0.0;
    for drive in [
        [-13.0, 0.2, -0.1],
        [0.0, 0.08, -0.06],
        [0.0, -0.05, 0.04],
        [4.0, 0.1, 0.04],
        [1.0, 0.05, -0.02],
    ] {
        let old = system.point();
        let frame = *system.step(&[0.0; 6], drive).unwrap();
        let (_, original_mean) = old
            .trial_over_time_with_mean_traction(frame.jump_m, system.sample_period_s(), 1e-8)
            .unwrap();
        let slip = (frame.jump_m[1] - old.jump_m()[1]).hypot(frame.jump_m[2] - old.jump_m()[2]);
        let pressure = -AREA * original_mean[0];
        if pressure > 0.0 {
            assert_eq!(
                frame.residual_friction_regime,
                Some(ReclosureRegime::Sliding)
            );
            close(
                frame.mean_force_n[1].hypot(frame.mean_force_n[2]),
                0.3 * pressure,
                2e-10,
            );
            close(
                frame.residual_friction_dissipation_j,
                0.3 * pressure * slip,
                1e-13,
            );
        }
        assert_eq!(frame.fracture_dissipation_j, 0.0);
        assert_eq!(frame.viscous_dissipation_j, 0.0);
        assert_eq!(
            system.point().dissipated_j_m2(),
            before_point.dissipated_j_m2()
        );
        assert_eq!(
            system.point().mixed_mode_history(),
            before_point.mixed_mode_history()
        );
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        assert!(frame.residual_friction_work_defect_j.abs() <= 1e-13);
        assert!(frame.material_evaluations <= config().max_evaluations);
        accumulated_work += frame.prescribed_work_j;
        accumulated_loss += frame.network_dissipation_j + frame.residual_friction_dissipation_j;
        friction_work += frame.residual_friction_dissipation_j;
        close(
            system.total_energy_j().unwrap() + accumulated_loss,
            initial + accumulated_work,
            2e-12,
        );
        for i in 0..3 {
            close(momentum(&system)[i], initial_momentum[i], 5e-15);
        }
    }
    assert!(friction_work > 0.0);
    assert_eq!(system.samples_rendered(), 6);
    assert_eq!(
        system.residual_contact_occupancy(),
        Some(ReclosureRegime::Open)
    );
    let open = system.last_frame().unwrap();
    assert_eq!(open.residual_friction_regime, Some(ReclosureRegime::Open));
    assert_eq!(open.mean_force_n, [0.0; 3]);
    assert_eq!(open.residual_friction_dissipation_j, 0.0);
}

#[test]
fn g4_native_enable_and_late_energy_refusals_keep_exact_history_frame_and_retry() {
    let (source, friction) = sources();
    let build = || native(&source, 6e-4);
    let mut system = build();
    let initial_states = states(&system);
    let initial_point = system.point();
    assert!(matches!(
        system.enable_frictional_reclosure(
            &friction,
            receiving(&source),
            ReclosureActivation::AtRest,
            config(),
            &CancelGate::new(),
        ),
        Err(CohesiveModalError::Reclosure(
            ReclosureError::NotFullyFailed
        ))
    ));
    assert_eq!(states(&system), initial_states);
    assert_eq!(system.point(), initial_point);
    assert_eq!(system.samples_rendered(), 0);
    assert_eq!(system.last_frame(), None);
    rupture(&mut system);
    let before = (
        states(&system),
        system.point(),
        system.samples_rendered(),
        system.last_frame().copied(),
    );
    let cancelled = CancelGate::new();
    cancelled.request();
    for (materials, gate, residual_config) in [
        (
            [receiving(&source)[1], receiving(&source)[0]],
            &CancelGate::new(),
            config(),
        ),
        (receiving(&source), &cancelled, config()),
        (
            receiving(&source),
            &CancelGate::new(),
            ReclosureConfig {
                max_evaluations: 2,
                ..config()
            },
        ),
    ] {
        assert!(
            system
                .enable_frictional_reclosure(
                    &friction,
                    materials,
                    ReclosureActivation::AtRest,
                    residual_config,
                    gate
                )
                .is_err()
        );
        assert_eq!(
            (
                states(&system),
                system.point(),
                system.samples_rendered(),
                system.last_frame().copied()
            ),
            before
        );
        assert!(system.residual_friction_source().is_none());
    }
    let mut replay = build();
    rupture(&mut replay);
    for owner in [&mut system, &mut replay] {
        owner
            .enable_frictional_reclosure(
                &friction,
                receiving(&source),
                ReclosureActivation::AtRest,
                config(),
                &CancelGate::new(),
            )
            .unwrap();
        owner.step(&[0.0; 6], [-13.0, 0.2, -0.1]).unwrap();
    }
    let before = (
        states(&system),
        system.point(),
        system.samples_rendered(),
        system.last_frame().copied(),
        system.residual_contact_occupancy(),
    );
    assert!(
        system
            .enable_frictional_reclosure(
                &friction,
                receiving(&source),
                ReclosureActivation::Sliding,
                config(),
                &CancelGate::new()
            )
            .is_err()
    );
    assert!(
        system
            .step_under_gate(&[0.0; 6], [0.0, 0.1, 0.0], &cancelled)
            .is_err()
    );
    let refused = system.step(&[0.0; 6], [-5.0, 0.2, 0.0]);
    assert!(matches!(
        refused,
        Err(CohesiveModalError::Network(ModalCouplingError::Budget {
            what: "total energy including failed interface",
            ..
        }))
    ));
    assert_eq!(
        (
            states(&system),
            system.point(),
            system.samples_rendered(),
            system.last_frame().copied(),
            system.residual_contact_occupancy()
        ),
        before
    );
    let drive = [0.0, -0.05, 0.04];
    assert_eq!(
        *system.step(&[0.0; 6], drive).unwrap(),
        *replay.step(&[0.0; 6], drive).unwrap()
    );
    assert_eq!(states(&system), states(&replay));
    assert_eq!(system.point(), replay.point());
    assert_eq!(
        system.residual_contact_occupancy(),
        replay.residual_contact_occupancy()
    );
}

#[test]
fn g1_zero_net_slip_does_not_relock_a_moving_native_interface() {
    let (source, friction) = sources();
    let mut system = native_with_connection(&source, 1.0, false);
    rupture(&mut system);
    let dt = system.sample_period_s();
    let close_rate = (-1e-4 - system.point().jump_m()[0]) / dt;
    system.step(&[0.0; 6], [close_rate, 0.0, 0.0]).unwrap();
    assert!(system.point().jump_m()[0] < 0.0);
    system
        .enable_frictional_reclosure(
            &friction,
            receiving(&source),
            ReclosureActivation::AtRest,
            config(),
            &CancelGate::new(),
        )
        .unwrap();
    let relative_velocity = |system: &CohesiveModalSystem<'_>, axis: usize| {
        system.components()[1].states()[axis].velocity_m_sqrt_kg_per_s / det::sqrt(MASSES[1])
            - system.components()[0].states()[axis].velocity_m_sqrt_kg_per_s / det::sqrt(MASSES[0])
    };
    let compliance = system.compliance_m_per_n();
    let old = system.point().jump_m();
    let pressure = -AREA * COMPRESSION * old[0];
    let force = 0.2 * pressure;
    let drive = [
        -relative_velocity(&system, 0) - compliance[0][0] * pressure / dt,
        -relative_velocity(&system, 1) + compliance[1][1] * force / dt,
        -relative_velocity(&system, 2),
    ];
    let frame = *system.step(&[0.0; 6], drive).unwrap();
    assert_eq!(
        frame.residual_friction_regime,
        Some(ReclosureRegime::Sticking)
    );
    close(frame.jump_m[1], old[1], 1e-14);
    close(frame.mean_force_n[1], force, 1e-11);
    close(
        relative_velocity(&system, 1) + drive[1],
        -compliance[1][1] * force / dt,
        1e-15,
    );
    assert_ne!(relative_velocity(&system, 1) + drive[1], 0.0);
    assert_eq!(
        system.residual_contact_occupancy(),
        Some(ReclosureRegime::Sliding)
    );
    assert_eq!(frame.residual_friction_dissipation_j, 0.0);
    assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);

    // This next stopping demand lies between kinetic and static capacity.
    // Motion must keep the kinetic law despite the previous zero net slip.
    let drive = [
        -relative_velocity(&system, 0) - compliance[0][0] * pressure / dt,
        -relative_velocity(&system, 1) + compliance[1][1] * 0.45 * pressure / dt,
        -relative_velocity(&system, 2),
    ];
    let frame = *system.step(&[0.0; 6], drive).unwrap();
    assert_eq!(
        frame.residual_friction_regime,
        Some(ReclosureRegime::Sliding)
    );
    close(frame.mean_force_n[1], 0.3 * pressure, 2e-10);
    assert!(frame.residual_friction_dissipation_j > 0.0);
    assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
}
