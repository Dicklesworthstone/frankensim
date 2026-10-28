//! Actual two-axis modal recoil and retained shear, with independent limits.
use fs_couple::modal_acoustic_time::{
    ModalAcousticMode, ModalAcousticState, ModalAcousticTimeBudget, ModalAcousticTimeModel,
};
use fs_couple::render::schedule::force::coupled::contact::friction::ModalFrictionConfig;
use fs_couple::render::schedule::force::coupled::contact::friction::tangential::planar::{
    PlanarTangentialModalDrive, PlanarTangentialModalSystem,
};
use fs_couple::render::schedule::force::coupled::{
    CoupledModalSystem, ModalAttachment, ModalConnection, ModalCouplingConfig,
};
use fs_exec::CancelGate;
use fs_math::c64::C64;
use fs_tribo::layer_contact::ThirdBodyLayerContact;
use fs_tribo::tangential::matrix::MatrixTangentialConfig;
use fs_tribo::{FrictionLaw, FrictionRegime, TangentialLaw};

const OMEGAS: [[f64; 2]; 2] = [[800.0, 1100.0], [1700.0, 2200.0]];
const SHAPES: [[[f64; 2]; 2]; 2] = [[[8.0, 2.0], [3.0, 6.0]], [[4.0, -1.0], [2.0, 5.0]]];

fn config() -> ModalFrictionConfig {
    ModalFrictionConfig {
        max_iterations: 96,
        maximum_force_n: 100.0,
        maximum_slip_speed_m_s: 10.0,
        force_absolute_tolerance_n: 1e-12,
        force_relative_tolerance: 1e-11,
        sticking_tolerance_m: 1e-12,
    }
}

fn matrix_config() -> MatrixTangentialConfig {
    MatrixTangentialConfig {
        maximum_evaluations: 512,
        displacement_tolerance_m: 1e-13,
        energy_tolerance_j: 1e-12,
        matrix_relative_tolerance: 1e-12,
    }
}

fn network(connected: bool) -> CoupledModalSystem {
    let models = OMEGAS
        .iter()
        .map(|omegas| {
            ModalAcousticTimeModel::try_new(
                48000,
                omegas
                    .iter()
                    .map(|&omega| ModalAcousticMode {
                        angular_frequency_rad_s: omega,
                        damping_ratio: 0.0,
                        pressure_per_modal_velocity: C64::new(1.0, 0.0),
                    })
                    .collect(),
                ModalAcousticTimeBudget::audible_reference(),
            )
            .unwrap()
        })
        .collect();
    let connections = if connected {
        vec![ModalConnection {
            left: ModalAttachment {
                component: 0,
                shapes: vec![1.0, 0.7],
            },
            right: ModalAttachment {
                component: 1,
                shapes: vec![1.0, -0.5],
            },
            stiffness_n_m: 3000.0,
            damping_n_s_m: 0.1,
            rest_extension_m: 0.0,
        }]
    } else {
        vec![]
    };
    CoupledModalSystem::new(
        models,
        connections,
        ModalCouplingConfig {
            max_modes: 8,
            max_connections: 4,
            max_setup_terms: 4096,
            nyquist_guard_fraction: 0.9,
            maximum_total_energy_j: 100.0,
            maximum_abs_pressure_pa: 1e6,
            maximum_abs_connection_force_n: 1e6,
            solve_relative_tolerance: 1e-11,
            energy_absolute_tolerance_j: 1e-12,
            energy_relative_tolerance: 1e-8,
        },
        &CancelGate::new(),
    )
    .unwrap()
}

fn attachments(component: usize, rotation: [[f64; 2]; 2]) -> [ModalAttachment; 2] {
    std::array::from_fn(|axis| ModalAttachment {
        component,
        shapes: (0..2)
            .map(|mode| {
                rotation[axis][0] * SHAPES[component][0][mode]
                    + rotation[axis][1] * SHAPES[component][1][mode]
            })
            .collect(),
    })
}

fn layer(yield_force: f64) -> ThirdBodyLayerContact {
    ThirdBodyLayerContact::new(2e6, yield_force, 0.5, 300.0, 0.01).unwrap()
}

fn system<L: TangentialLaw>(
    law: L,
    connected: bool,
    rotation: [[f64; 2]; 2],
) -> PlanarTangentialModalSystem<L> {
    PlanarTangentialModalSystem::new(
        network(connected),
        attachments(0, rotation),
        Some(attachments(1, rotation)),
        law,
        config(),
        matrix_config(),
        &CancelGate::new(),
    )
    .unwrap()
}

fn drive(speed: [f64; 2]) -> PlanarTangentialModalDrive {
    PlanarTangentialModalDrive {
        speed_m_s: speed,
        normal_force_n: 1.0,
        contact_temperature_k: Some(300.0),
        nominal_pressure_pa: None,
    }
}

fn states<L: TangentialLaw>(system: &PlanarTangentialModalSystem<L>) -> Vec<ModalAcousticState> {
    system
        .network()
        .components()
        .iter()
        .flat_map(|m| m.states().iter().copied())
        .collect()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, tolerance {tolerance:e}"
    );
}

const IDENTITY: [[f64; 2]; 2] = [[1.0, 0.0], [0.0, 1.0]];

#[test]
fn full_modal_compliance_matches_independent_coupled_elastic_solution() {
    let mut owner = system(layer(100.0), false, IDENTITY);
    let dt = owner.network().sample_period_s();
    let mut compliance = [[0.0; 2]; 2];
    for i in 0..2 {
        for j in 0..2 {
            for body in 0..2 {
                for mode in 0..2 {
                    let omega = OMEGAS[body][mode];
                    compliance[i][j] +=
                        SHAPES[body][i][mode] * SHAPES[body][j][mode] * (1.0 - (omega * dt).cos())
                            / (omega * omega);
                }
            }
            close(owner.compliance_m_per_n()[i][j], compliance[i][j], 1e-20);
        }
    }
    assert!(compliance[0][1].abs() > 0.3 * compliance[0][0]);
    let speed = [0.01, -0.007];
    let q = 1e6; // Midpoint effort of the original k=2e6 elastic layer.
    let matrix = [
        [1.0 + q * compliance[0][0], q * compliance[0][1]],
        [q * compliance[1][0], 1.0 + q * compliance[1][1]],
    ];
    let rhs = [q * speed[0] * dt, q * speed[1] * dt];
    let determinant = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
    let expected = [
        (rhs[0] * matrix[1][1] - rhs[1] * matrix[0][1]) / determinant,
        (matrix[0][0] * rhs[1] - matrix[1][0] * rhs[0]) / determinant,
    ];
    let frame = *owner.step(&[0.0; 4], drive(speed)).unwrap();
    for i in 0..2 {
        close(frame.traction_n[i], expected[i], 2e-10);
        close(
            owner.state().elastic_slip_m()[i],
            frame.slip_distance_m[i],
            1e-13,
        );
    }
    assert_eq!(frame.friction_dissipation_j, 0.0);
    assert!(frame.contact_energy_j > 0.0);
    assert!(frame.network_energy_j > 0.0);
    close(owner.total_energy_j().unwrap(), frame.drive_work_j, 1e-12);
}

#[test]
fn nonproportional_layer_history_drives_both_bodies_and_closes_actual_work() {
    let mut owner = system(layer(0.03), true, IDENTITY);
    let mut supplied = 0.0;
    let mut loss = 0.0;
    let mut released_storage = false;
    for k in 0..400 {
        let speed = match k / 100 {
            0 => [0.1, 0.0],
            1 => [0.0, 0.1],
            2 => [-0.1, 0.0],
            _ => [0.0, -0.1],
        };
        let frame = *owner
            .step(&[0.001, -0.002, 0.003, -0.001], drive(speed))
            .unwrap();
        supplied += frame.external_work_j + frame.drive_work_j;
        loss += frame.network_dissipation_j + frame.friction_dissipation_j;
        released_storage |= frame.contact_stored_change_j < -1e-12;
        assert!(frame.friction_dissipation_j >= 0.0);
        assert!(frame.contact_work_residual_j.abs() <= frame.energy_tolerance_j);
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
    }
    assert!(released_storage);
    assert!(loss > 1e-6);
    assert!(
        owner
            .state()
            .plastic_slip_m()
            .iter()
            .all(|s| s.abs() > 1e-7)
    );
    assert!(
        owner
            .network()
            .components()
            .iter()
            .all(|m| m.states().iter().any(|s| s.displacement_m_sqrt_kg != 0.0))
    );
    close(owner.total_energy_j().unwrap() + loss, supplied, 1e-10);
}

#[test]
fn tangent_basis_rotation_preserves_original_body_motion_and_work() {
    let angle: f64 = 0.63;
    let rotation = [[angle.cos(), -angle.sin()], [angle.sin(), angle.cos()]];
    let rotate =
        |x: [f64; 2]| std::array::from_fn(|i| rotation[i][0] * x[0] + rotation[i][1] * x[1]);
    let mut a = system(layer(0.03), true, IDENTITY);
    let mut b = system(layer(0.03), true, rotation);
    for k in 0..160 {
        let speed = if k < 80 { [0.08, -0.04] } else { [-0.03, 0.09] };
        let fa = *a.step(&[0.0; 4], drive(speed)).unwrap();
        let fb = *b.step(&[0.0; 4], drive(rotate(speed))).unwrap();
        let expected_force = rotate(fa.traction_n);
        for i in 0..2 {
            close(fb.traction_n[i], expected_force[i], 2e-8);
        }
        close(fb.friction_dissipation_j, fa.friction_dissipation_j, 1e-12);
        close(fb.contact_energy_j, fa.contact_energy_j, 1e-12);
    }
    for (sa, sb) in states(&a).iter().zip(states(&b)) {
        close(sa.displacement_m_sqrt_kg, sb.displacement_m_sqrt_kg, 1e-11);
        close(
            sa.velocity_m_sqrt_kg_per_s,
            sb.velocity_m_sqrt_kg_per_s,
            1e-8,
        );
    }
}

#[test]
fn anisotropic_coulomb_sticks_then_slides_against_actual_relative_motion() {
    let law = FrictionLaw::Coulomb {
        static_mu: 0.6,
        kinetic_mu: 0.4,
    };
    let mut owner = system(law, false, IDENTITY);
    let stuck = *owner.step(&[0.0; 4], drive([1e-5, 2e-5])).unwrap();
    assert_eq!(stuck.regime, FrictionRegime::Sticking);
    assert!(stuck.slip_distance_m.iter().all(|s| s.abs() < 1e-12));
    assert_eq!(stuck.friction_dissipation_j, 0.0);
    let frame = *owner.step(&[0.0; 4], drive([0.2, 0.15])).unwrap();
    assert_eq!(frame.regime, FrictionRegime::Sliding);
    close(frame.traction_n[0].hypot(frame.traction_n[1]), 0.4, 1e-10);
    close(
        frame.traction_n[0] * frame.slip_distance_m[1]
            - frame.traction_n[1] * frame.slip_distance_m[0],
        0.0,
        1e-12,
    );
    close(
        frame.friction_dissipation_j,
        frame.traction_n[0] * frame.slip_distance_m[0]
            + frame.traction_n[1] * frame.slip_distance_m[1],
        1e-12,
    );
}

#[test]
fn invalid_or_cancelled_trials_preserve_both_axis_histories_and_network_clock() {
    let mut owner = system(layer(0.03), true, IDENTITY);
    owner.step(&[0.0; 4], drive([0.05, -0.04])).unwrap();
    let state = owner.state().clone();
    let bodies = states(&owner);
    let frame = *owner.last_frame().unwrap();
    let sample = owner.network().samples_rendered();
    let mut bad = drive([0.02, 0.01]);
    bad.contact_temperature_k = Some(301.0);
    assert!(owner.step(&[0.0; 4], bad).is_err());
    assert!(owner.step(&[f64::NAN; 4], drive([0.0, 0.0])).is_err());
    let gate = CancelGate::new();
    gate.request();
    assert!(
        owner
            .step_under_gate(&[0.0; 4], drive([0.02, 0.01]), &gate)
            .is_err()
    );
    assert_eq!(owner.state(), &state);
    assert_eq!(states(&owner), bodies);
    assert_eq!(owner.last_frame(), Some(&frame));
    assert_eq!(owner.network().samples_rendered(), sample);
    let retried = *owner.step(&[0.0; 4], drive([0.02, 0.01])).unwrap();
    assert_eq!(retried.sample, sample + 1);
}

#[test]
fn a_fixed_transverse_direction_keeps_the_physical_singular_response() {
    let left = [
        ModalAttachment {
            component: 0,
            shapes: vec![8.0, 2.0],
        },
        ModalAttachment {
            component: 0,
            shapes: vec![0.0, 0.0],
        },
    ];
    let mut owner = PlanarTangentialModalSystem::new(
        network(false),
        left,
        None,
        layer(100.0),
        config(),
        matrix_config(),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(owner.compliance_m_per_n()[1], [0.0; 2]);
    let frame = *owner.step(&[0.0; 4], drive([0.01, 0.02])).unwrap();
    close(frame.traction_n[1], 1e6 * 0.02 / 48000.0, 2e-10);
    close(frame.slip_distance_m[1], 0.02 / 48000.0, 1e-15);
    assert!(frame.contact_energy_j > 0.0);
    assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
}
