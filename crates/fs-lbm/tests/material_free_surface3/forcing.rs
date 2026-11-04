//! Physical acceleration histories through the original source and timestep.

use super::*;
use fs_lbm::d3q19::RigidFrame3;
use fs_lbm::d3q19::freesurface3::ConversionRedistribution3;
use fs_lbm::material::free_surface::MaterialRigidFrame3;

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn si_rigid_motion_matches_native_units_and_keeps_its_source_identity() {
    let source = liquid(0, 300.0, None);
    let gate = CancelGate::new_clock_free();
    for (dx, dt, conversion) in [
        (0.001, 0.0002, ConversionRedistribution3::Local),
        (0.002, 0.0004, ConversionRedistribution3::ConservativeLocal),
    ] {
        let setup = plan_free_surface_from_material(&source, dx, dt, 0.1).unwrap();
        let mut fluid = setup
            .free_surface3(domain(), |_, _, z| z < 6)
            .unwrap()
            .with_conversion_redistribution(conversion);
        let mut native = fluid.lattice().clone();
        let frame = MaterialRigidFrame3 {
            origin_m: [7.5 * dx, 7.5 * dx, 3.5 * dx],
            angular_velocity_rad_s: [2.0, -3.0, 1.0],
            angular_acceleration_rad_s2: [-1.0, 2.0, 3.0],
            origin_acceleration_m_s2: [0.4, -0.2, 0.1],
        };
        let body = [0.3, -0.1, -0.5];
        let lattice = RigidFrame3 {
            origin: frame.origin_m.map(|v| v / dx),
            angular_velocity: frame.angular_velocity_rad_s.map(|v| v * dt),
            angular_acceleration: frame.angular_acceleration_rad_s2.map(|v| (v * dt) * dt),
            origin_acceleration: frame.origin_acceleration_m_s2.map(|v| (v * dt) / (dx / dt)),
        };
        lattice
            .prepare_free_surface_step(&mut native, body.map(|v| (v * dt) / (dx / dt)), &gate)
            .unwrap();
        native.step_with_gate(&gate).unwrap();
        fluid.step_in_rigid_frame_si(frame, body, &gate).unwrap();
        assert_eq!(
            fluid.lattice().grid().state_bits(),
            native.grid().state_bits()
        );
        assert_eq!(fluid.lattice().ledger_mass(), native.ledger_mass());
        assert!(std::ptr::eq(fluid.scaling().material(), &source));
        let physical = fluid.sample_at_si(8, 8, 2).unwrap().unwrap();
        let numerical = native.sample_at(8, 8, 2).unwrap();
        for d in 0..3 {
            near(
                physical.velocity_m_s[d],
                numerical.velocity[d] * dx / dt,
                1e-15,
            );
        }
    }
}

#[test]
fn si_acceleration_switches_drive_the_native_trajectory_and_freefall_cancels_gravity() {
    let source = liquid(0, 300.0, None);
    let gate = CancelGate::new_clock_free();
    for dt in [0.0001, 0.0002] {
        let setup = plan_free_surface_from_material(&source, 0.001, dt, 0.1).unwrap();
        let mut input = domain();
        input.gravity_m_s2 = [0.3, -0.2, 0.1];
        let mut fluid = setup.free_surface3(input, |_, _, _| true).unwrap();
        // The historical binary constructor seeds raw rest populations, so
        // the accepted physical velocity includes its initial half source.
        // Free fall preserves that velocity; each later source supplies only
        // the independently known dt*(body-origin) increment.
        let mut velocity = fluid.sample_at_si(8, 8, 6).unwrap().unwrap().velocity_m_s;
        for (body, origin) in [
            (input.gravity_m_s2, input.gravity_m_s2),
            ([0.5, -0.4, 0.2], [0.1, -0.2, 0.05]),
            ([-0.3, 0.1, -0.2], [-0.1, 0.05, 0.1]),
        ] {
            fluid
                .step_in_rigid_frame_si(
                    MaterialRigidFrame3 {
                        origin_acceleration_m_s2: origin,
                        ..Default::default()
                    },
                    body,
                    &gate,
                )
                .unwrap();
            let actual = fluid.sample_at_si(8, 8, 6).unwrap().unwrap();
            for d in 0..3 {
                velocity[d] += dt * (body[d] - origin[d]);
                near(actual.velocity_m_s[d], velocity[d], 2e-14);
            }
        }
    }
}

#[test]
fn si_conversion_refusal_and_cancellation_preserve_the_complete_accepted_state() {
    let source = liquid(0, 300.0, None);
    let gate = CancelGate::new_clock_free();
    let mut fluid = plan(&source)
        .free_surface3(domain(), |_, _, z| z < 6)
        .unwrap();
    fluid
        .step_in_rigid_frame_si(MaterialRigidFrame3::default(), [0.2, -0.1, 0.1], &gate)
        .unwrap();
    let before = fluid.lattice().grid().state_bits();
    let mass = fluid.ledger_mass_kg().unwrap();
    let force = fluid.lattice().force_density_at(4, 4, 3);
    for value in [f64::NAN, f64::MAX, f64::from_bits(1)] {
        let frame = MaterialRigidFrame3 {
            angular_acceleration_rad_s2: [value, 0.0, 0.0],
            ..Default::default()
        };
        assert!(
            fluid
                .step_in_rigid_frame_si(frame, [0.0; 3], &gate)
                .is_err()
        );
        assert_eq!(fluid.lattice().grid().state_bits(), before);
        assert_eq!(fluid.ledger_mass_kg().unwrap(), mass);
        assert_eq!(fluid.lattice().force_density_at(4, 4, 3), force);
    }
    gate.request();
    assert!(matches!(
        fluid.step_in_rigid_frame_si(MaterialRigidFrame3::default(), [0.0; 3], &gate),
        Err(MaterialLbmError::FreeSurface(FreeSurfaceError3::Cancelled))
    ));
    assert_eq!(fluid.lattice().grid().state_bits(), before);
    assert_eq!(fluid.lattice().force_density_at(4, 4, 3), force);
}
