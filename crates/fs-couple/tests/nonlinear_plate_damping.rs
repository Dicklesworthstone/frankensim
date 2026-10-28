//! G1/G3: physical force and surface ports retain complete nonlinear-body loss.
use fs_couple::thin_plate::VkBody;
use fs_nlmodal::{SosModalStorage, StressChannel};

fn storage(omegas: &[f64]) -> SosModalStorage {
    SosModalStorage {
        omegas: omegas.to_vec(),
        channels: vec![],
    }
}

#[test]
fn full_damping_body_retains_common_motion_and_differential_decay() {
    let omega: f64 = 8.0;
    let damping: f64 = 1.2;
    let force: f64 = 0.6;
    let dt: f64 = 0.2;
    let density: f64 = 1.2;
    let distance: f64 = 1.7;
    let area: f64 = 0.02;
    for sign in [1.0, -1.0] {
        // A physical damper on q1-q2 must leave equal motion undamped.
        let mut body = VkBody::from_modal_storage_with_damping(
            storage(&[omega, omega]),
            &[damping, -damping, -damping, damping],
            &[1.0, sign],
            vec![area, sign * area],
        )
        .unwrap();
        let expected_velocity = if sign == 1.0 {
            force / omega * (omega * dt).sin()
        } else {
            let damped = (omega * omega - damping * damping).sqrt();
            force / damped * (-damping * dt).exp() * (damped * dt).sin()
        };
        let expected_flow = 2.0 * area * expected_velocity;
        let (pressure, flow) = body
            .drive_and_radiate(force, dt, density, distance)
            .unwrap();
        assert!((flow - expected_flow).abs() < 2e-13);
        assert_eq!(flow, body.volume_velocity());
        // This compact observer returns step-mean acceleration; at rest its
        // pressure is rho/(2*pi*r) times the accepted volume-flow increment/dt.
        let expected_pressure =
            density * expected_flow / (2.0 * core::f64::consts::PI * distance * dt);
        assert!((pressure - expected_pressure).abs() < 2e-13);
    }
}

#[test]
fn full_damping_body_reuses_the_exact_linear_limit_under_subdivision() {
    let make = || {
        VkBody::from_modal_storage_with_damping(
            storage(&[2.0, 5.0]),
            &[0.6, 0.2, 0.2, 0.1],
            &[0.3, -0.2],
            vec![0.02, -0.01],
        )
        .unwrap()
    };
    let mut whole = make();
    let mut split = make();
    whole.drive_and_radiate(0.7, 0.2, 1.2, 1.0).unwrap();
    for _ in 0..8 {
        split.drive_and_radiate(0.7, 0.025, 1.2, 1.0).unwrap();
    }
    assert!((whole.volume_velocity() - split.volume_velocity()).abs() < 2e-13);
}

#[test]
fn public_nonlinear_body_admits_storage_and_all_physical_projections() {
    let make = |storage, force: &[f64], areas| {
        VkBody::from_modal_storage_with_damping(storage, &[0.0; 4], force, areas)
    };
    assert!(make(storage(&[2.0, 5.0]), &[1.0], vec![0.02, -0.01]).is_err());
    assert!(make(storage(&[2.0, 5.0]), &[1.0, 0.0], vec![0.02]).is_err());
    assert!(make(storage(&[2.0, 5.0]), &[f64::NAN, 0.0], vec![0.02, -0.01]).is_err());
    assert!(make(storage(&[2.0, 5.0]), &[1.0, 0.0], vec![0.02, f64::INFINITY]).is_err());
    for channel in [
        StressChannel {
            coefficient: -1.0,
            coupling: vec![1.0, 0.0, 0.0, 1.0],
        },
        StressChannel {
            coefficient: 1.0,
            coupling: vec![1.0],
        },
        StressChannel {
            coefficient: 1.0,
            coupling: vec![1.0, 1.0, 0.0, 1.0],
        },
    ] {
        let mut invalid = storage(&[2.0, 5.0]);
        invalid.channels.push(channel);
        assert!(make(invalid, &[1.0, 0.0], vec![0.02, -0.01]).is_err());
    }
}
