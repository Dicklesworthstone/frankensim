//! G1/G3/G4: acceleration and the original source viscosity share real motion.

use super::*;
use fs_lbm::material::free_surface::MaterialRigidFrame3;

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[test]
fn accelerating_vessel_closes_the_original_changing_viscosity_and_midpoint_force() {
    let gate = CancelGate::new_clock_free();
    let source = card(0.003, 1000.0, 0.001);
    let binding = law(&source);
    let body = [0.0, 0.0, -0.02];
    with_cx(&gate, |cx| {
        let mut flow = MaterialFreeSurfaceRheology3::from_cells(
            &binding,
            DX,
            DT,
            0.1,
            domain(body),
            &initial(true),
            solve(),
            cx,
            &gate,
        )
        .unwrap();
        let mut stationary = flow.clone();
        let mass = flow.ledger_mass_kg().unwrap();
        let initial_tau = flow.rheology_at(CELL).unwrap().collision_tau;
        let mut maximum_tau_change = 0.0_f64;
        for step in 0..6 {
            let midpoint_s = (f64::from(step) + 0.5) * DT;
            let frame = MaterialRigidFrame3 {
                origin_m: [2.5 * DX, 4.5 * DX, 1.5 * DX],
                angular_velocity_rad_s: [1.5, -2.0, 4.0 + 0.2 * midpoint_s],
                angular_acceleration_rad_s2: [0.0, 0.0, 0.2],
                origin_acceleration_m_s2: [0.1, -0.08, 0.03],
            };
            let before = flow.lattice().sample_at(CELL[0], CELL[1], CELL[2]).unwrap();
            flow.step_in_rigid_frame_si(frame, body, &gate).unwrap();
            stationary
                .step_in_rigid_frame_si(MaterialRigidFrame3::default(), body, &gate)
                .unwrap();

            let force = flow
                .lattice()
                .force_density_at(CELL[0], CELL[1], CELL[2])
                .unwrap();
            let acceleration = force.map(|f| f / before.density * DX / (DT * DT));
            let midpoint_velocity = core::array::from_fn(|axis| {
                before.velocity[axis] * DX / DT + 0.5 * DT * acceleration[axis]
            });
            let offset =
                core::array::from_fn(|axis| (CELL[axis] as f64 + 0.5) * DX - frame.origin_m[axis]);
            let coriolis = cross(frame.angular_velocity_rad_s, midpoint_velocity);
            let euler = cross(frame.angular_acceleration_rad_s2, offset);
            let centrifugal = cross(
                frame.angular_velocity_rad_s,
                cross(frame.angular_velocity_rad_s, offset),
            );
            // Independent SI momentum-source equation, without invoking the
            // production midpoint inverse or its coordinate conversion helper.
            for axis in 0..3 {
                near(
                    acceleration[axis]
                        + 2.0 * coriolis[axis]
                        + euler[axis]
                        + centrifugal[axis]
                        + frame.origin_acceleration_m_s2[axis],
                    body[axis],
                    2e-11,
                );
            }

            // The accepted force-corrected population stress and mu=.1+.003*g
            // give a quadratic constitutive root, independent of solve_local.
            let rhs = numerator_with_force(&flow, force, CELL) / DT;
            let rate = 2.0 * rhs / (0.8 + (0.64 + 0.036 * rhs).sqrt());
            let actual = flow.rheology_at(CELL).unwrap();
            near(actual.shear_rate_per_s, rate, 2e-9);
            near(actual.collision_tau, 0.8 + 0.009 * rate, 2e-11);
            maximum_tau_change = maximum_tau_change.max((actual.collision_tau - initial_tau).abs());
            let material = flow.material_at(CELL).unwrap();
            near(
                material.property("dynamic_viscosity").unwrap().value_si(),
                0.1 + 0.003 * rate,
                2e-11,
            );
            assert_eq!(material.card_identity(), source.content_hash());
            assert!(material.query_point().contains(&("T".into(), 400.0)));
            assert!(
                material
                    .query_point()
                    .contains(&("pressure".into(), 100000.0))
            );
            near(flow.ledger_mass_kg().unwrap(), mass, 1e-10 * mass);
        }
        assert!(
            maximum_tau_change > 1e-5,
            "viscosity must evolve with actual strain"
        );
        let moving = flow
            .sample_at_si(CELL[0], CELL[1], CELL[2])
            .unwrap()
            .unwrap();
        let fixed = stationary
            .sample_at_si(CELL[0], CELL[1], CELL[2])
            .unwrap()
            .unwrap();
        assert!(
            moving
                .velocity_m_s
                .iter()
                .zip(fixed.velocity_m_s)
                .any(|(a, b)| (a - b).abs() > 1e-6),
            "frame acceleration must alter the physical liquid trajectory"
        );
    });
}

#[test]
fn freefall_and_changing_uniform_acceleration_do_not_become_material_shear() {
    let gate = CancelGate::new_clock_free();
    let source = card(0.003, 1000.0, 0.001);
    let binding = law(&source);
    let gravity = [0.2, -0.1, 0.05];
    let mut cells = initial(true);
    for cell in &mut cells {
        cell.velocity_m_s = [0.0; 3];
    }
    with_cx(&gate, |cx| {
        let mut flow = MaterialFreeSurfaceRheology3::from_cells(
            &binding,
            DX,
            DT,
            0.1,
            domain(gravity),
            &cells,
            solve(),
            cx,
            &gate,
        )
        .unwrap();
        let mut expected = flow
            .sample_at_si(CELL[0], CELL[1], CELL[2])
            .unwrap()
            .unwrap()
            .velocity_m_s;
        for (body, origin) in [
            (gravity, gravity),
            ([0.35, -0.15, 0.08], [0.12, 0.2, -0.05]),
            ([-0.3, 0.1, -0.2], [-0.1, 0.05, 0.1]),
        ] {
            flow.step_in_rigid_frame_si(
                MaterialRigidFrame3 {
                    origin_acceleration_m_s2: origin,
                    ..Default::default()
                },
                body,
                &gate,
            )
            .unwrap();
            let actual = flow
                .sample_at_si(CELL[0], CELL[1], CELL[2])
                .unwrap()
                .unwrap();
            for axis in 0..3 {
                expected[axis] += DT * (body[axis] - origin[axis]);
                near(actual.velocity_m_s[axis], expected[axis], 3e-14);
            }
            let sample = flow.rheology_at(CELL).unwrap();
            near(sample.shear_rate_per_s, 0.0, 5e-9);
            near(sample.collision_tau, 0.8, 1e-10);
        }
    });
}

#[test]
fn frame_source_or_late_constitutive_refusal_preserves_all_accepted_state() {
    let gate = CancelGate::new_clock_free();
    with_cx(&gate, |cx| {
        for (max_rate, maximum_iterations) in [(1e-6, 128), (1000.0, 1)] {
            let source = card(0.003, max_rate, 0.001);
            let binding = law(&source);
            let mut flow = MaterialFreeSurfaceRheology3::from_cells(
                &binding,
                DX,
                DT,
                0.1,
                domain([0.0; 3]),
                &initial(true),
                RheologySolve {
                    maximum_iterations,
                    ..solve()
                },
                cx,
                &gate,
            )
            .unwrap();
            let before = flow.lattice().grid().state_bits();
            let force = flow.lattice().force_density_at(CELL[0], CELL[1], CELL[2]);
            let mass = flow.ledger_mass_kg().unwrap();
            let conversions = flow.lattice().conversions();
            let sample = *flow.rheology_at(CELL).unwrap();
            let receipt = flow.material_at(CELL).unwrap().clone();
            let frame = MaterialRigidFrame3 {
                angular_velocity_rad_s: [0.0, 0.0, 3.0],
                origin_acceleration_m_s2: [0.1, -0.1, 0.05],
                ..Default::default()
            };
            let refused = flow
                .step_in_rigid_frame_si(frame, [0.0; 3], &gate)
                .unwrap_err();
            assert!(matches!(
                refused,
                FreeSurfaceRheologyError::Material(_)
                    | FreeSurfaceRheologyError::Convergence { .. }
            ));
            assert_eq!(flow.lattice().grid().state_bits(), before);
            assert_eq!(
                flow.lattice().force_density_at(CELL[0], CELL[1], CELL[2]),
                force
            );
            assert_eq!(flow.ledger_mass_kg().unwrap(), mass);
            assert_eq!(flow.lattice().conversions(), conversions);
            assert_eq!(*flow.rheology_at(CELL).unwrap(), sample);
            assert_eq!(*flow.material_at(CELL).unwrap(), receipt);
            assert!(
                flow.step_in_rigid_frame_si(frame, [f64::NAN, 0.0, 0.0], &gate)
                    .is_err()
            );
            let cancelled = CancelGate::new_clock_free();
            cancelled.request();
            assert!(matches!(
                flow.step_in_rigid_frame_si(frame, [0.0; 3], &cancelled),
                Err(FreeSurfaceRheologyError::Flow(
                    MaterialLbmError::FreeSurface(FreeSurfaceError3::Cancelled)
                ))
            ));
            assert_eq!(flow.lattice().grid().state_bits(), before);
            assert_eq!(
                flow.lattice().force_density_at(CELL[0], CELL[1], CELL[2]),
                force
            );
        }
    });
}
