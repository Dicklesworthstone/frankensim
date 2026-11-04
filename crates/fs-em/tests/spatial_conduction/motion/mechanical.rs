//! G1/G3/G4: the accepted moving-current field reaches native mechanics.

use super::*;
use fs_em::spatial::conduction::motion::{SpatialLorentzLoadError, SpatialLorentzLoadLimits};
use fs_ga::{Mat3, Se3, Twist, Vec3, Wrench};

fn load_limits(problem: &SpatialConduction) -> SpatialLorentzLoadLimits {
    SpatialLorentzLoadLimits {
        max_vertices: problem.positions_m().len(),
        max_cells: problem.mesh().tets.len(),
    }
}

#[test]
fn bar_loads_preserve_analytic_force_moment_power_and_refinement() {
    for xs in [&[0.0, 1.0][..], &[0.0, 0.25, 1.0][..]] {
        let problem = scalar_problem(xs, 2.0);
        let cells = problem.mesh().tets.len();
        let selected: Vec<_> = (0..cells).collect();
        let reverse: Vec<_> = selected.iter().copied().rev().collect();
        let nodal_velocity = vec![[0.0, 0.8, 0.0]; problem.positions_m().len()];
        for voltage in [0.0, 0.4, 1.8] {
            let field = problem
                .solve_with_motion(
                    &ends(&problem, 1.0, Some(voltage), 0.0),
                    SpatialConductionMotion {
                        velocity_m_s: &vec![[0.0, 0.8, 0.0]; cells],
                        flux_density_t: &vec![[0.0, 0.0, 1.5]; cells],
                    },
                    controls(),
                    &mut || false,
                )
                .unwrap();
            let loads = field
                .lorentz_loads_for_velocity(
                    &selected,
                    &nodal_velocity,
                    0.0,
                    load_limits(&problem),
                    &mut || false,
                )
                .unwrap();
            let reordered = field
                .lorentz_loads(&reverse, load_limits(&problem), &mut || false)
                .unwrap();
            assert_eq!(loads.nodal_forces_n(), reordered.nodal_forces_n());
            // Unit-volume uniform bar: F_y=-B*sigma*(u*B-V), acting at
            // (1/2,1/2,1/2). These closed forms do not use the load assembly.
            let force_y = -3.0 * (1.2 - voltage);
            let expected = [-0.5 * force_y, 0.0, 0.5 * force_y, 0.0, force_y, 0.0];
            let wrench = loads.wrench(&Se3::identity(), [0.0; 3]).unwrap();
            for (a, b) in wrench.to_array().into_iter().zip(expected) {
                near(a, b);
            }
            near(
                loads.power_w(&nodal_velocity, &mut || false).unwrap(),
                0.8 * force_y,
            );
            near(field.mechanical_power_w(), 0.8 * force_y);
            if xs.len() == 3 {
                // The sum of the four P1 basis functions on the shared plane
                // is a one-dimensional hat with area (0.25+0.75)/2.
                let plane_force: f64 = problem
                    .positions_m()
                    .iter()
                    .zip(loads.nodal_forces_n())
                    .filter(|(p, _)| p[0] == 0.25)
                    .map(|(_, f)| f[1])
                    .sum();
                near(plane_force, 0.5 * force_y);
                let last_layer: Vec<_> = (6..12).collect();
                let partial = field
                    .lorentz_loads_for_velocity(
                        &last_layer,
                        &nodal_velocity,
                        0.0,
                        load_limits(&problem),
                        &mut || false,
                    )
                    .unwrap();
                let wrench = partial.wrench(&Se3::identity(), [0.0; 3]).unwrap();
                near(wrench.force.y, 0.75 * force_y);
                near(wrench.torque.x, -0.375 * force_y);
                near(wrench.torque.z, 0.75 * 0.625 * force_y);
                near(
                    partial.power_w(&nodal_velocity, &mut || false).unwrap(),
                    last_layer
                        .iter()
                        .map(|&cell| field.cell_mechanical_power_w()[cell])
                        .sum(),
                );
            }
        }
    }
}

fn affine_velocity(p: [f64; 3]) -> [f64; 3] {
    [
        0.2 + 0.3 * p[0] - 0.1 * p[1],
        0.8 - 0.2 * p[0] + 0.4 * p[2],
        -0.1 + 0.2 * p[1],
    ]
}

#[test]
fn p1_velocity_means_preserve_physical_power_and_arbitrary_rigid_frame_work() {
    let problem = scalar_problem(&[0.0, 0.35, 1.0], 3.0);
    let nodal: Vec<_> = problem
        .positions_m()
        .iter()
        .copied()
        .map(affine_velocity)
        .collect();
    let velocity: Vec<_> = problem
        .mesh()
        .tets
        .iter()
        .map(|tet| {
            let centroid = std::array::from_fn(|d| {
                tet.iter()
                    .map(|&v| problem.positions_m()[v as usize][d])
                    .sum::<f64>()
                    / 4.0
            });
            affine_velocity(centroid)
        })
        .collect();
    let field = problem
        .solve_with_motion(
            &ends(&problem, 1.0, Some(0.7), 0.0),
            SpatialConductionMotion {
                velocity_m_s: &velocity,
                flux_density_t: &[[0.4, -0.6, 0.8]; 12],
            },
            controls(),
            &mut || false,
        )
        .unwrap();
    let selected: Vec<_> = (0..12).collect();
    let loads = field
        .lorentz_loads_for_velocity(&selected, &nodal, 1e-14, load_limits(&problem), &mut || {
            false
        })
        .unwrap();
    near(
        loads.power_w(&nodal, &mut || false).unwrap(),
        field.mechanical_power_w(),
    );
    near(
        field.terminal_power_w(),
        field.joule_power_w() + field.mechanical_power_w(),
    );

    let placement = Se3::exp(Twist::new(
        Vec3::new(0.2, -0.4, 0.3),
        Vec3::new(2.0, -1.0, 3.0),
    ))
    .unwrap();
    let origin = Vec3::new(4.0, -2.0, 1.0);
    let twist = Twist::new(Vec3::new(0.7, -0.2, 0.4), Vec3::new(-0.5, 0.3, 1.1));
    let inverse = placement.inverse().unwrap();
    let virtual_velocity: Vec<_> = problem
        .positions_m()
        .iter()
        .map(|p| {
            let arm = placement
                .transform_point(Vec3::new(p[0], p[1], p[2]))
                .unwrap()
                - origin;
            let v = inverse
                .transform_vector(twist.linear + twist.angular.cross(arm))
                .unwrap();
            [v.x, v.y, v.z]
        })
        .collect();
    let wrench = loads
        .wrench(&placement, [origin.x, origin.y, origin.z])
        .unwrap();
    near(
        wrench.pairing(twist),
        loads.power_w(&virtual_velocity, &mut || false).unwrap(),
    );
}

#[test]
fn motoring_current_drives_native_body_momentum_and_exact_impulse_work() {
    use fs_mbd::articulated::{
        ArticulatedModel, FreeFloatingBaseState, JointModel, Link, SpatialInertia,
        apply_free_floating_impulses,
    };
    let problem = scalar_problem(&[0.0, 1.0], 2.0);
    let field = problem
        .solve_with_motion(
            &ends(&problem, 1.0, Some(1.8), 0.0),
            SpatialConductionMotion {
                velocity_m_s: &[[0.0, 0.8, 0.0]; 6],
                flux_density_t: &[[0.0, 0.0, 1.5]; 6],
            },
            controls(),
            &mut || false,
        )
        .unwrap();
    let loads = field
        .lorentz_loads_for_velocity(
            &[0, 1, 2, 3, 4, 5],
            &[[0.0, 0.8, 0.0]; 8],
            0.0,
            load_limits(&problem),
            &mut || false,
        )
        .unwrap();
    let wrench = loads.wrench(&Se3::identity(), [0.0; 3]).unwrap();
    let base = FreeFloatingBaseState::new(
        Se3::identity(),
        Twist::new(Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.8, 0.0)),
    );
    near(wrench.pairing(base.twist_body), field.mechanical_power_w());
    let inertia = SpatialInertia::new(
        2.0,
        Vec3::new(0.0, 0.0, 0.0),
        Mat3 {
            m: [3.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 5.0],
        },
    )
    .unwrap();
    let model = ArticulatedModel::new(vec![Link::new(
        "moving conductor",
        None,
        Se3::identity(),
        JointModel::FIXED,
        inertia,
    )])
    .unwrap();
    let dt = 0.07;
    let impulse = Wrench::new(wrench.torque.scale(dt), wrench.force.scale(dt));
    let response = apply_free_floating_impulses(&model, base, &[], &[], &[], &[impulse]).unwrap();
    let momentum = Wrench::new(
        response.angular_momentum_change_world,
        response.linear_momentum_change_world,
    );
    for (a, b) in momentum.to_array().into_iter().zip(impulse.to_array()) {
        near(a, b);
    }
    // The solved field gives F_y=1.8 N and torque=(-0.9,0,0.9) N m.
    near(
        response.base_after.twist_body.linear.y,
        0.8 + dt * 1.8 / 2.0,
    );
    near(response.base_after.twist_body.angular.x, -dt * 0.9 / 3.0);
    near(response.base_after.twist_body.angular.z, dt * 0.9 / 5.0);
    // A held-load impulse changes velocity. Its actual kinetic work includes
    // the quadratic term; the DC snapshot supplies only the initial rate.
    let quadratic =
        0.5 * dt * dt * (1.8_f64.powi(2) / 2.0 + 0.9_f64.powi(2) / 3.0 + 0.9_f64.powi(2) / 5.0);
    near(
        response.kinetic_energy_change,
        dt * field.mechanical_power_w() + quadratic,
    );
    near(response.kinetic_work_residual, 0.0);
    assert!(field.joule_power_w() > 0.0 && field.mechanical_power_w() > 0.0);
}

#[test]
fn incompatible_velocity_and_invalid_selections_refuse_without_mutating_field() {
    let problem = scalar_problem(&[0.0, 0.5, 1.0], 2.0);
    let field = problem
        .solve_with_motion(
            &ends(&problem, 1.0, Some(0.4), 0.0),
            SpatialConductionMotion {
                velocity_m_s: &[[0.0, 0.8, 0.0]; 12],
                flux_density_t: &[[0.0, 0.0, 1.5]; 12],
            },
            controls(),
            &mut || false,
        )
        .unwrap();
    let selected: Vec<_> = (0..12).collect();
    let limits = load_limits(&problem);
    let nodal = vec![[0.0, 0.8, 0.0]; problem.positions_m().len()];
    let before_current = field.current_density_a_m2().to_vec();
    let mut calls = 0;
    let before = field
        .lorentz_loads_for_velocity(&selected, &nodal, 0.0, limits, &mut || {
            calls += 1;
            false
        })
        .unwrap();
    for cutoff in [1, calls / 2, calls] {
        let mut count = 0;
        assert!(matches!(
            field.lorentz_loads_for_velocity(&selected, &nodal, 0.0, limits, &mut || {
                count += 1;
                count == cutoff
            },),
            Err(SpatialLorentzLoadError::Cancelled)
        ));
    }
    for invalid in [&[][..], &[0, 0][..], &[12][..]] {
        assert!(field.lorentz_loads(invalid, limits, &mut || false).is_err());
    }
    assert!(
        field
            .lorentz_loads(
                &selected,
                SpatialLorentzLoadLimits {
                    max_cells: 11,
                    ..limits
                },
                &mut || false
            )
            .is_err()
    );
    assert!(
        field
            .lorentz_loads(
                &selected,
                SpatialLorentzLoadLimits {
                    max_vertices: 11,
                    ..limits
                },
                &mut || false
            )
            .is_err()
    );
    let incompatible: Vec<_> = problem
        .positions_m()
        .iter()
        .map(|p| [0.0, 0.8 + 0.1 * (2.0 * p[0] - 1.0), 0.0])
        .collect();
    // The wrong trace has exactly the same total work in this symmetric bar.
    // Checking each selected cell still detects the physical mismatch.
    near(
        before.power_w(&incompatible, &mut || false).unwrap(),
        field.mechanical_power_w(),
    );
    assert!(matches!(
        field.lorentz_loads_for_velocity(&selected, &incompatible, 1e-12, limits, &mut || false,),
        Err(SpatialLorentzLoadError::Invalid(
            "mechanical velocity differs from prescribed cell motion"
        ))
    ));
    for tolerance in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(
            field
                .lorentz_loads_for_velocity(&selected, &nodal, tolerance, limits, &mut || false)
                .is_err()
        );
    }
    for invalid in [&nodal[..1], &vec![[f64::NAN, 0.0, 0.0]; nodal.len()][..]] {
        assert!(
            field
                .lorentz_loads_for_velocity(&selected, invalid, 0.0, limits, &mut || false)
                .is_err()
        );
    }
    assert_eq!(field.current_density_a_m2(), before_current);
    assert_eq!(
        field
            .lorentz_loads(&selected, limits, &mut || false)
            .unwrap()
            .nodal_forces_n(),
        before.nodal_forces_n()
    );
}
