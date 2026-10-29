//! Original native motion, sourced adhesion/shear and independent work closure.
#[path = "../../fs-contact/tests/adhesive_dmt_coupled/common.rs"]
mod common;

use common::{Z, close, config, free, sources};
use fs_contact::normal_patch::dmt::coupled::{DmtCoupledLoadLine, DmtCoupledRegime};
use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_couple::render::schedule::force::coupled::contact::adhesive_dmt::{
    AdhesiveDmtModalError, AdhesiveDmtModalSystem,
};
use fs_couple::render::schedule::force::coupled::{
    CoupledModalSystem, ModalAttachment, ModalConnection, ModalCouplingConfig, ModalCouplingError,
};
use fs_exec::CancelGate;
use fs_math::{c64::C64, det};

const MA: f64 = 2e-9;
const MB: f64 = 3e-9;
const RATE: u32 = 5_000_000;

fn columns() -> [Vec<f64>; 3] {
    core::array::from_fn(|i| {
        let mut column = vec![0.0; 6];
        column[i] = -1.0 / det::sqrt(MA);
        column[3 + i] = 1.0 / det::sqrt(MB);
        column
    })
}
fn network(relative: [f64; 3], cross: bool, rate: u32) -> CoupledModalSystem {
    let models = [MA, MB]
        .into_iter()
        .enumerate()
        .map(|(side, mass)| {
            let mut model = ModalAcousticTimeModel::try_new_with_free_coordinates(
                rate,
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
            let states: Vec<_> = (0..3)
                .map(|i| ModalAcousticState {
                    displacement_m_sqrt_kg: 0.0,
                    velocity_m_sqrt_kg_per_s: det::sqrt(mass)
                        * if side == 0 {
                            -MB / (MA + MB)
                        } else {
                            MA / (MA + MB)
                        }
                        * relative[i],
                })
                .collect();
            model.restore_states(&states).unwrap();
            model
        })
        .collect();
    let links = if cross {
        vec![ModalConnection {
            left: ModalAttachment {
                component: 0,
                shapes: vec![
                    0.7 / det::sqrt(MA),
                    0.3 / det::sqrt(MA),
                    -0.2 / det::sqrt(MA),
                ],
            },
            right: ModalAttachment {
                component: 1,
                shapes: vec![
                    0.7 / det::sqrt(MB),
                    0.3 / det::sqrt(MB),
                    -0.2 / det::sqrt(MB),
                ],
            },
            stiffness_n_m: 3e4,
            damping_n_s_m: 0.0,
            rest_extension_m: 0.0,
        }]
    } else {
        vec![]
    };
    CoupledModalSystem::new(
        models,
        links,
        ModalCouplingConfig {
            max_modes: 6,
            max_connections: 2,
            max_setup_terms: 4096,
            nyquist_guard_fraction: 0.9,
            maximum_total_energy_j: 1e-10,
            maximum_abs_pressure_pa: 1000.0,
            maximum_abs_connection_force_n: 1.0,
            solve_relative_tolerance: 1e-11,
            energy_absolute_tolerance_j: 1e-26,
            energy_relative_tolerance: 1e-9,
        },
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}
fn make<'a>(
    n: &'a fs_contact::normal_patch::dmt::SourceDmtSphereContact,
    s: &'a fs_contact::interface_binding::adhesive_shear::SourceAdhesiveShear,
    velocity: [f64; 3],
    initial: [f64; 3],
    cross: bool,
    rate: u32,
) -> AdhesiveDmtModalSystem<'a> {
    let materials = n.normal_state().elastic().interface().surface_materials();
    let mut bounds = config();
    // Cross-axis stored energy can detach the junction. Continue the native
    // free flight inside an explicitly declared travel range after attraction
    // has become exactly zero; this is not source-law extrapolation.
    bounds.minimum_indentation_m = -100.0 * Z;
    AdhesiveDmtModalSystem::new(
        network(velocity, cross, rate),
        n,
        s,
        [&materials[0], &materials[1]],
        initial,
        columns(),
        bounds,
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}
fn momentum(system: &AdhesiveDmtModalSystem<'_>) -> [f64; 3] {
    core::array::from_fn(|i| {
        det::sqrt(MA) * system.components()[0].states()[i].velocity_m_sqrt_kg_per_s
            + det::sqrt(MB) * system.components()[1].states()[i].velocity_m_sqrt_kg_per_s
    })
}
fn states(system: &AdhesiveDmtModalSystem<'_>) -> Vec<ModalAcousticState> {
    system
        .components()
        .iter()
        .flat_map(|model| model.states().iter().copied())
        .collect()
}

#[test]
fn native_two_body_adhesion_and_slip_preserve_momenta_and_complete_work() {
    let (n, s) = sources(0.01, 1e7, 6e6);
    for cross in [false, true] {
        let mut system = make(
            &n,
            &s,
            [0.0, 1e-3, -0.5e-3],
            [0.1 * Z, 0.0, 0.0],
            cross,
            RATE,
        );
        assert_eq!(system.mode_count(), 6);
        assert_eq!(system.state().regime(), DmtCoupledRegime::Sliding);
        assert!(std::ptr::eq(system.normal(), &n) && std::ptr::eq(system.shear(), &s));
        let energy = system.total_energy_j().unwrap();
        let initial_momentum = momentum(&system);
        let mut friction = 0.0;
        let mut network_loss = 0.0;
        for sample in 0..150 {
            let old_position = system.state().position_m();
            let frame = *system.step(&[0.0; 6], [0.0; 3]).unwrap_or_else(|error| {
                panic!("cross={cross}, sample={sample}, position={old_position:?}: {error}")
            });
            friction += frame.friction_dissipation_j;
            network_loss += frame.network_dissipation_j;
            assert!(frame.friction_dissipation_j >= 0.0);
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
            assert!(frame.evaluations <= config().maximum_evaluations);
            assert_eq!(frame.position_m, system.state().position_m());
            for i in 0..3 {
                close(momentum(&system)[i], initial_momentum[i], 2e-24, 1e-10);
                assert!(frame.compatibility_residual_m[i].abs() <= frame.compatibility_tolerance_m);
            }
            close(
                system.total_energy_j().unwrap() + friction + network_loss,
                energy,
                1e-25,
                2e-8,
            );
        }
        assert!(friction > 0.0);
        close(network_loss, 0.0, 1e-26, 0.0);
        assert_eq!(system.samples_rendered(), 150);
    }
}

#[test]
fn native_cross_axis_compliance_prescribed_motion_and_actual_area_close_work() {
    let (n, s) = sources(0.01, 1e7, 6e6);
    let old = [0.2 * Z, 0.0, 0.0];
    let mut system = make(&n, &s, [0.0; 3], old, true, RATE);
    let c = system.compliance_m_per_n();
    assert!(c[0][1] != 0.0 && c[0][2] != 0.0 && c[1][2] != 0.0);
    let owner = DmtCoupledLoadLine::new(&n, &s, c, config()).unwrap();
    let target = [0.25 * Z, 0.2 * Z, -0.1 * Z];
    let path = owner.path(old[0], target[0]).unwrap();
    let distance = target[1].hypot(target[2]);
    let force = [
        path.mean_normal_force_n,
        path.kinetic_capacity_n * target[1] / distance,
        path.kinetic_capacity_n * target[2] / distance,
    ];
    let predicted = free(c, target, force);
    let rate = core::array::from_fn(|i| (predicted[i] - old[i]) / system.sample_period_s());
    let before = system.total_energy_j().unwrap();
    let frame = *system.step(&[0.0; 6], rate).unwrap();
    assert_eq!(frame.regime, DmtCoupledRegime::Sliding);
    for i in 0..3 {
        close(frame.position_m[i], target[i], 1e-21, 1e-9);
        close(frame.mean_force_n[i], force[i], 1e-16, 1e-8);
    }
    close(
        frame.path.mean_real_area_m2,
        path.mean_real_area_m2,
        1e-29,
        1e-9,
    );
    assert!(frame.path.mean_real_area_m2 < frame.path.contact.contact_area_m2);
    close(
        system.total_energy_j().unwrap() - before
            + frame.network_dissipation_j
            + frame.friction_dissipation_j,
        frame.prescribed_work_j,
        1e-26,
        2e-9,
    );
    assert_ne!(frame.network_energy_j, 0.0);
}

#[test]
fn native_gap_attraction_has_no_junction_friction_or_heat() {
    let (n, s) = sources(0.01, 1e7, 6e6);
    let mut system = make(
        &n,
        &s,
        [0.0, 1e-3, -0.5e-3],
        [-0.8 * Z, 0.0, 0.0],
        false,
        RATE,
    );
    let before = system.total_energy_j().unwrap();
    let frame = *system.step(&[0.0; 6], [0.0; 3]).unwrap();
    assert_eq!(frame.regime, DmtCoupledRegime::Released);
    assert!(frame.mean_force_n[0] < 0.0 && frame.path.normal_work_j < 0.0);
    assert_eq!([frame.mean_force_n[1], frame.mean_force_n[2]], [0.0; 2]);
    assert_eq!(frame.path.mean_real_area_m2, 0.0);
    assert_eq!(frame.friction_dissipation_j, 0.0);
    close(system.total_energy_j().unwrap(), before, 1e-26, 1e-9);
}

#[test]
fn native_initial_motion_retains_kinetic_branch_below_static_strength() {
    let (n, s) = sources(0.01, 1e7, 4e6);
    for velocity in [0.0, 1e-5] {
        let mut system = make(
            &n,
            &s,
            [0.0, velocity, 0.0],
            [0.2 * Z, 0.0, 0.0],
            false,
            RATE,
        );
        let c = system.compliance_m_per_n();
        let owner = DmtCoupledLoadLine::new(&n, &s, c, config()).unwrap();
        let path = owner.path(0.2 * Z, 0.2 * Z).unwrap();
        let effort = [
            path.mean_normal_force_n,
            0.7 * path.static_capacity_n - system.sample_period_s() * velocity / c[1][1],
            0.0,
        ];
        let columns = columns();
        let external: Vec<f64> = (0..6)
            .map(|i| (0..3).map(|j| columns[j][i] * effort[j]).sum())
            .collect();
        let frame = *system.step(&external, [0.0; 3]).unwrap();
        assert_eq!(
            frame.regime,
            if velocity == 0.0 {
                DmtCoupledRegime::Sticking
            } else {
                DmtCoupledRegime::Sliding
            }
        );
        assert_eq!(frame.friction_dissipation_j > 0.0, velocity != 0.0);
    }
}

#[test]
fn native_zero_mean_slip_does_not_relock_a_moving_endpoint() {
    let (n, s) = sources(0.01, 1e7, 6e6);
    let initial = [0.2 * Z, 0.0, 0.0];
    let mut system = make(&n, &s, [0.0, 1e-7, 0.0], initial, false, RATE);
    let c = system.compliance_m_per_n();
    let owner = DmtCoupledLoadLine::new(&n, &s, c, config()).unwrap();
    let path = owner.path(initial[0], initial[0]).unwrap();
    let columns = columns();
    let force = |effort: [f64; 3]| -> Vec<f64> {
        (0..6)
            .map(|i| (0..3).map(|j| columns[j][i] * effort[j]).sum())
            .collect()
    };
    let first = *system
        .step(&force([path.mean_normal_force_n, 0.0, 0.0]), [0.0; 3])
        .unwrap();
    assert_eq!(first.regime, DmtCoupledRegime::Sticking);
    close(first.position_m[1], initial[1], 1e-22, 0.0);
    let velocity = columns[1]
        .iter()
        .zip(states(&system))
        .map(|(weight, state)| weight * state.velocity_m_sqrt_kg_per_s)
        .sum::<f64>();
    close(velocity, -1e-7, 1e-20, 1e-10);
    assert_eq!(system.state().regime(), DmtCoupledRegime::Sliding);
    // The next stopping demand lies between kinetic and static strength.
    // Promoting the first zero-increment path to static occupancy would lock it.
    let effort = [
        path.mean_normal_force_n,
        0.8 * path.static_capacity_n - system.sample_period_s() * velocity / c[1][1],
        0.0,
    ];
    let second = *system.step(&force(effort), [0.0; 3]).unwrap();
    assert_eq!(second.regime, DmtCoupledRegime::Sliding);
    assert!(second.friction_dissipation_j > 0.0);
    assert!(second.energy_residual_j.abs() <= second.energy_tolerance_j);
}

#[test]
fn native_cancellation_source_mismatch_and_range_refusal_preserve_all_history() {
    let (n, s) = sources(0.01, 1e7, 6e6);
    let mut system = make(&n, &s, [0.0, 1e-3, 0.0], [0.1 * Z, 0.0, 0.0], false, RATE);
    system.step(&[0.0; 6], [0.0; 3]).unwrap();
    let old = (
        states(&system),
        system.state(),
        *system.last_frame().unwrap(),
        system.samples_rendered(),
    );
    let gate = CancelGate::new_clock_free();
    gate.request();
    assert!(matches!(
        system.step_under_gate(&[0.0; 6], [0.0; 3], &gate),
        Err(AdhesiveDmtModalError::Network(
            ModalCouplingError::Cancelled
        ))
    ));
    assert!(system.step(&[0.0; 6], [1e9, 0.0, 0.0]).is_err());
    assert_eq!(
        (
            states(&system),
            system.state(),
            *system.last_frame().unwrap(),
            system.samples_rendered()
        ),
        old
    );
    let materials = n.normal_state().elastic().interface().surface_materials();
    assert!(
        AdhesiveDmtModalSystem::new(
            network([0.0; 3], false, RATE),
            &n,
            &s,
            [&materials[1], &materials[0]],
            [0.1 * Z, 0.0, 0.0],
            columns(),
            config(),
            &CancelGate::new_clock_free()
        )
        .is_err()
    );
    let (foreign, _) = sources(0.02, 1e7, 6e6);
    assert!(
        AdhesiveDmtModalSystem::new(
            network([0.0; 3], false, RATE),
            &foreign,
            &s,
            [&materials[0], &materials[1]],
            [0.1 * Z, 0.0, 0.0],
            columns(),
            config(),
            &CancelGate::new_clock_free()
        )
        .is_err()
    );
    let mut replay = make(&n, &s, [0.0, 1e-3, 0.0], [0.1 * Z, 0.0, 0.0], false, RATE);
    replay.step(&[0.0; 6], [0.0; 3]).unwrap();
    assert_eq!(
        system.step(&[0.0; 6], [0.0; 3]).unwrap(),
        replay.step(&[0.0; 6], [0.0; 3]).unwrap()
    );
}

#[test]
fn native_physical_time_refinement_converges_without_losing_work_closure() {
    let (n, s) = sources(0.01, 1e7, 6e6);
    let results: Vec<f64> = [1, 2, 4, 8]
        .into_iter()
        .map(|refinement| {
            let mut system = make(
                &n,
                &s,
                [0.0, 1e-3, -0.5e-3],
                [0.1 * Z, 0.0, 0.0],
                false,
                RATE * refinement,
            );
            for _ in 0..36 * refinement {
                let frame = *system.step(&[0.0; 6], [0.0; 3]).unwrap();
                assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
            }
            system.state().position_m()[0]
        })
        .collect();
    let errors: Vec<_> = results[..3]
        .iter()
        .map(|v| (v - results[3]).abs())
        .collect();
    assert!(
        errors[1] < 0.4 * errors[0] && errors[2] < 0.4 * errors[1],
        "normal coordinate refinement errors: {errors:?}"
    );
}
