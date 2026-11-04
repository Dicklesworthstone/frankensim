//! G0/G1/G3/G4: analytic three-dimensional motional-EMF power conversion.
use super::*;
use fs_em::spatial::conduction::motion::SpatialConductionMotion;

#[path = "motion/mechanical.rs"]
mod mechanical;

fn scalar_problem(xs: &[f64], sigma: f64) -> SpatialConduction {
    let source = ResolvedConductorState::try_from_material(material(
        ElectricalPropertyBasis::Conductivity,
        &[sigma],
    ))
    .unwrap();
    let (positions, tets) = slab(xs);
    let cells = tets.len();
    SpatialConduction::new(positions, tets, vec![source; cells], limits(), &mut || {
        false
    })
    .unwrap()
}

#[test]
fn a_three_dimensional_bar_shorts_generates_and_motors_with_physical_heat() {
    for sigma in [2.0, 5.0] {
        let problem = scalar_problem(&[0.0, 0.4, 1.0], sigma);
        let cells = problem.mesh().tets.len();
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
            // Uniform-bar closed form, independent of tetrahedral assembly.
            // u cross B = +1.2 x; the two end electrodes impose E=-V x.
            let j = sigma * (1.2 - voltage);
            let heat = j * (1.2 - voltage);
            let mechanical = -1.2 * j;
            let electrical = -voltage * j;
            for (v, p) in problem.positions_m().iter().enumerate() {
                near(field.voltage_v()[v], voltage * p[0]);
            }
            for cell in 0..cells {
                near(field.current_density_a_m2()[cell][0], j);
                near(field.current_density_a_m2()[cell][1], 0.0);
                near(field.current_density_a_m2()[cell][2], 0.0);
                near(field.lorentz_force_density_n_m3()[cell][1], -1.5 * j);
                let volume = problem.cell_volume_m3(cell).unwrap();
                near(field.cell_joule_power_w()[cell], volume * heat);
                near(field.cell_mechanical_power_w()[cell], volume * mechanical);
            }
            near(field.joule_power_w(), heat);
            near(field.mechanical_power_w(), mechanical);
            near(field.electrical_power_w(), electrical);
            near(field.terminal_power_w(), electrical);
            near(
                field.joule_power_w() + field.mechanical_power_w(),
                electrical,
            );
            near(field.power_balance_residual_w(), 0.0);
            near(field.work_pairing_residual_w(), 0.0);
            near(field.hodge_power_residual_w(), 0.0);
            assert!(field.joule_power_w() > 0.0);
            if voltage == 0.0 {
                near(field.terminal_power_w(), 0.0);
                assert!(field.mechanical_power_w() < 0.0);
            } else if voltage < 1.2 {
                assert!(field.terminal_power_w() < 0.0);
                assert!(field.mechanical_power_w() < 0.0);
            } else {
                assert!(field.terminal_power_w() > 0.0);
                assert!(field.mechanical_power_w() > 0.0);
            }
        }
    }
}

#[test]
fn a_floating_anisotropic_conductor_screens_a_conservative_motional_field() {
    let s = std::f64::consts::FRAC_1_SQRT_2;
    let basis = ElectricalPropertyBasis::Conductivity;
    let source = ResolvedOrthotropicConductor3::try_from_material(
        material(basis, &[2.0, 8.0, 3.0]),
        basis,
        std::array::from_fn(|i| PropertyKey::new(format!("principal_{i}"), basis.dimensions())),
        [[s, s, 0.0], [-s, s, 0.0], [0.0, 0.0, 1.0]],
    )
    .unwrap();
    let (positions, tets) = slab(&[0.0, 0.4, 1.0]);
    let problem =
        SpatialConduction::new(positions, tets, vec![source; 12], limits(), &mut || false).unwrap();
    let field = problem
        .solve_with_motion(
            &SpatialConductionBoundary::default(),
            SpatialConductionMotion {
                velocity_m_s: &[[2.0, -1.0, 0.5]; 12],
                flux_density_t: &[[0.3, 0.6, -0.2]; 12],
            },
            controls(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(field.gauges(), &[0]);
    // The open body develops V=(u cross B).x while J vanishes everywhere.
    for (&voltage, p) in field.voltage_v().iter().zip(problem.positions_m()) {
        near(voltage, -0.1 * p[0] + 0.55 * p[1] + 1.5 * p[2]);
    }
    for cell in 0..12 {
        for d in 0..3 {
            near(field.conductor_field_v_m()[cell][d], 0.0);
            near(field.current_density_a_m2()[cell][d], 0.0);
            near(field.lorentz_force_density_n_m3()[cell][d], 0.0);
        }
    }
    near(field.joule_power_w(), 0.0);
    near(field.mechanical_power_w(), 0.0);
    near(field.terminal_power_w(), 0.0);
    near(field.boundary_current_a(), 0.0);
    assert!(field.nodal_boundary_current_a().iter().all(Option::is_none));
}

#[test]
fn full_tensor_power_and_force_rotate_with_the_scene_and_scale_with_volume() {
    let s = std::f64::consts::FRAC_1_SQRT_2;
    for rotated in [false, true] {
        let rotate = |v: [f64; 3]| {
            if rotated {
                [s * (v[0] - v[2]), v[1], s * (v[0] + v[2])]
            } else {
                v
            }
        };
        for scale in [1.0, 0.25] {
            let basis = ElectricalPropertyBasis::Resistivity;
            let source = ResolvedOrthotropicConductor3::try_from_material(
                material(basis, &[0.5, 0.125, 1.0 / 3.0]),
                basis,
                std::array::from_fn(|i| {
                    PropertyKey::new(format!("principal_{i}"), basis.dimensions())
                }),
                [[s, s, 0.0], [-s, s, 0.0], [0.0, 0.0, 1.0]].map(rotate),
            )
            .unwrap();
            let (positions, tets) = slab(&[0.0, 1.0]);
            let positions: Vec<_> = positions
                .into_iter()
                .map(|p| rotate(p.map(|v| scale * v)))
                .collect();
            let problem =
                SpatialConduction::new(positions, tets, vec![source; 6], limits(), &mut || false)
                    .unwrap();
            let electric = rotate([-1.0, 2.0, -3.0]);
            let boundary = SpatialConductionBoundary {
                voltage: problem
                    .exterior_faces()
                    .map(|face| ConductionVoltageFace {
                        face,
                        voltage_v: problem.mesh().faces[face].map(|v| {
                            let p = problem.positions_m()[v as usize];
                            13.0 - (0..3).map(|d| electric[d] * p[d]).sum::<f64>()
                        }),
                    })
                    .collect(),
                current: vec![],
            };
            let field = problem
                .solve_with_motion(
                    &boundary,
                    SpatialConductionMotion {
                        velocity_m_s: &[rotate([1.0, -2.0, 0.5]); 6],
                        flux_density_t: &[rotate([0.4, -0.6, 0.8]); 6],
                    },
                    controls(),
                    &mut || false,
                )
                .unwrap();
            // Original-frame sigma=[[5,-3,0],[-3,5,0],[0,0,3]].
            // E'=(-2.3,1.4,-2.8); independently J=(-15.7,13.9,-8.4).
            let current = rotate([-15.7, 13.9, -8.4]);
            let force = rotate([6.08, 9.2, 3.86]);
            for cell in 0..6 {
                for d in 0..3 {
                    near(field.electric_field_v_m()[cell][d], electric[d]);
                    near(field.current_density_a_m2()[cell][d], current[d]);
                    near(field.lorentz_force_density_n_m3()[cell][d], force[d]);
                }
            }
            near(field.joule_power_w(), 79.09 * scale.powi(3));
            near(field.mechanical_power_w(), -10.39 * scale.powi(3));
            near(field.electrical_power_w(), 68.7 * scale.powi(3));
            near(field.terminal_power_w(), 68.7 * scale.powi(3));
        }
    }
}

#[test]
fn prescribed_current_and_layered_sources_keep_heat_separate_from_motion_work() {
    let (positions, tets) = slab(&[0.0, 1.0, 2.0]);
    let sources = [2.0, 5.0].map(contacts::conductor);
    let problem = SpatialConduction::new(
        positions,
        tets,
        (0..12).map(|cell| sources[cell / 6].clone()).collect(),
        limits(),
        &mut || false,
    )
    .unwrap();
    let motion: Vec<_> = (0..12)
        .map(|cell| [0.0, [0.3, 1.2][cell / 6] / 2.0, 0.0])
        .collect();
    for current in [1.25, -0.1] {
        for reference in [None, Some(0.0)] {
            let field = problem
                .solve_with_motion(
                    &ends(&problem, 2.0, reference, current),
                    SpatialConductionMotion {
                        velocity_m_s: &motion,
                        flux_density_t: &[[0.0, 0.0, 2.0]; 12],
                    },
                    controls(),
                    &mut || false,
                )
                .unwrap();
            for (&voltage, p) in field.voltage_v().iter().zip(problem.positions_m()) {
                near(
                    voltage,
                    (2.0 * current + 0.3) * p[0].min(1.0)
                        + (5.0 * current + 1.2) * (p[0] - 1.0).max(0.0),
                );
            }
            for cell in 0..12 {
                near(field.current_density_a_m2()[cell][0], -current);
                near(
                    field.cell_joule_power_w()[cell],
                    current * current * [2.0, 5.0][cell / 6] / 6.0,
                );
                near(
                    field.cell_mechanical_power_w()[cell],
                    current * [0.3, 1.2][cell / 6] / 6.0,
                );
            }
            near(field.joule_power_w(), 7.0 * current * current);
            near(field.mechanical_power_w(), 1.5 * current);
            near(
                field.terminal_power_w(),
                7.0 * current * current + 1.5 * current,
            );
            assert_eq!(field.gauges().is_empty(), reference.is_some());
            if current < 0.0 {
                assert!(field.terminal_power_w() < 0.0);
            }
        }
    }
}

#[test]
fn zero_motion_matches_the_native_field_and_refusals_preserve_the_owner() {
    let problem = scalar_problem(&[0.0, 0.5, 1.0], 3.0);
    let boundary = ends(&problem, 1.0, Some(2.0), 0.0);
    let stationary = problem.solve(&boundary, controls(), &mut || false).unwrap();
    let motion = SpatialConductionMotion {
        velocity_m_s: &[[0.0; 3]; 12],
        flux_density_t: &[[1.0, 2.0, 3.0]; 12],
    };
    let mut calls = 0;
    let field = problem
        .solve_with_motion(&boundary, motion, controls(), &mut || {
            calls += 1;
            false
        })
        .unwrap();
    assert_eq!(field.voltage_v(), stationary.voltage_v);
    assert_eq!(
        field.current_density_a_m2(),
        stationary.current_density_a_m2
    );
    assert_eq!(field.cell_joule_power_w(), stationary.cell_power_w);
    assert_eq!(
        field
            .electric_field_v_m()
            .iter()
            .map(|v| v.map(f64::to_bits))
            .collect::<Vec<_>>(),
        stationary
            .electric_field_v_m
            .iter()
            .map(|v| v.map(f64::to_bits))
            .collect::<Vec<_>>(),
    );
    assert_eq!(field.joule_power_w(), stationary.power_w);
    assert_eq!(field.terminal_power_w(), stationary.supplied_power_w);
    assert_eq!(field.mechanical_power_w(), 0.0);
    assert_eq!(field.problem().conductors(), problem.conductors());
    assert_eq!(field.flux_density_t(), motion.flux_density_t);
    for cutoff in [1, calls / 2, calls] {
        let mut count = 0;
        assert!(matches!(
            problem.solve_with_motion(&boundary, motion, controls(), &mut || {
                count += 1;
                count == cutoff
            }),
            Err(SpatialConductionError::Cancelled)
        ));
    }
    for velocity in [&[[0.0; 3]][..], &[[f64::NAN, 0.0, 0.0]; 12][..]] {
        assert!(
            problem
                .solve_with_motion(
                    &boundary,
                    SpatialConductionMotion {
                        velocity_m_s: velocity,
                        ..motion
                    },
                    controls(),
                    &mut || false,
                )
                .is_err()
        );
    }
    assert!(
        problem
            .solve_with_motion(
                &boundary,
                SpatialConductionMotion {
                    velocity_m_s: &[[0.0, f64::MAX, 0.0]; 12],
                    flux_density_t: &[[0.0, 0.0, 2.0]; 12],
                },
                controls(),
                &mut || false,
            )
            .is_err()
    );
    assert_eq!(
        problem
            .solve(&boundary, controls(), &mut || false)
            .unwrap()
            .voltage_v,
        stationary.voltage_v,
    );
}

#[test]
fn finite_contact_motion_requires_its_own_surface_work_law() {
    let sources = [2.0, 5.0].map(contacts::conductor);
    let source = contacts::interface(&sources[0], &sources[1], 1e5, 1.0);
    let (mut positions, mut tets) = slab(&[0.0, 1.0]);
    let (other, other_tets) = slab(&[1.0, 2.0]);
    positions.extend(other);
    tets.extend(other_tets.into_iter().map(|tet| tet.map(|v| v + 8)));
    let problem = SpatialConduction::new_with_contacts(
        positions,
        tets,
        (0..12).map(|cell| sources[cell / 6].clone()).collect(),
        contacts::contacts(source),
        limits(),
        &mut || false,
    )
    .unwrap();
    assert!(matches!(
        problem.solve_with_motion(
            &ends(&problem, 2.0, Some(1.0), 0.0),
            SpatialConductionMotion {
                velocity_m_s: &[[0.0, 1.0, 0.0]; 12],
                flux_density_t: &[[0.0, 0.0, 1.0]; 12],
            },
            controls(),
            &mut || false,
        ),
        Err(SpatialConductionError::Invalid(
            "moving electrical contacts require an owned surface work and heat law"
        ))
    ));
}

#[test]
fn owned_motion_snapshot_moves_native_buffers_and_keeps_power_semantics() {
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
    let buffers = (
        field.voltage_v().as_ptr(),
        field.current_density_a_m2().as_ptr(),
        field.conductor_field_v_m().as_ptr(),
        field.cell_joule_power_w().as_ptr(),
        field.cell_mechanical_power_w().as_ptr(),
        field.lorentz_force_density_n_m3().as_ptr(),
        field.velocity_m_s().as_ptr(),
        field.flux_density_t().as_ptr(),
    );
    let had_report = field.solve_report().is_some();
    let snapshot = field.into_snapshot();
    // A coupled candidate can now move this exact material/topology owner.
    let retained_owner = problem;
    assert_eq!(
        snapshot.current_density_a_m2().len(),
        retained_owner.mesh().tets.len()
    );
    assert_eq!(
        buffers,
        (
            snapshot.voltage_v().as_ptr(),
            snapshot.current_density_a_m2().as_ptr(),
            snapshot.conductor_field_v_m().as_ptr(),
            snapshot.cell_joule_power_w().as_ptr(),
            snapshot.cell_mechanical_power_w().as_ptr(),
            snapshot.lorentz_force_density_n_m3().as_ptr(),
            snapshot.velocity_m_s().as_ptr(),
            snapshot.flux_density_t().as_ptr(),
        )
    );
    near(snapshot.joule_power_w(), 1.28);
    near(snapshot.mechanical_power_w(), -1.92);
    near(snapshot.electrical_power_w(), -0.64);
    near(snapshot.terminal_power_w(), -0.64);
    near(snapshot.boundary_current_a(), 0.0);
    near(snapshot.maximum_free_current_residual_a(), 0.0);
    near(snapshot.power_balance_residual_w(), 0.0);
    near(snapshot.work_pairing_residual_w(), 0.0);
    near(snapshot.hodge_power_residual_w(), 0.0);
    assert_eq!(snapshot.solve_report().is_some(), had_report);
    assert!(snapshot.gauges().is_empty());
    assert_eq!(
        snapshot.nodal_applied_current_a().len(),
        retained_owner.positions_m().len()
    );
    assert!(
        snapshot
            .nodal_boundary_current_a()
            .iter()
            .any(Option::is_some)
    );
    for electric in snapshot.electric_field_v_m() {
        near(electric[0], -0.4);
    }
}
