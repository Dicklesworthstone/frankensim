//! G1/G3/G4: independent transport solutions, conservation and atomic refusal.
use super::*;

fn controls(dt: f64) -> MaterialAdvectionControls {
    MaterialAdvectionControls {
        projection: DensityProjectionControls {
            duration_s: dt,
            relative_tolerance: 1e-12,
            maximum_volume_defect_m3_s: 1e-10,
            maximum_energy_defect_j: 1e-8,
        },
        max_cells: 100,
        max_connections: 1000,
        max_iterations: 1000,
        relative_tolerance: 1e-12,
        maximum_mass_defect_kg: 1e-9,
        maximum_enthalpy_defect_j: 1e-7,
    }
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 2e-11 * expected.abs().max(1.0),
        "{actual:.17e} != {expected:.17e}"
    );
}

fn input<'a>(
    mass: &'a [f64],
    h: &'a [f64],
    source: &'a [MaterialVolumeSource],
) -> MaterialAdvectionInput<'a> {
    MaterialAdvectionInput {
        mass_kg: mass,
        specific_enthalpy_j_kg: h,
        volume_sources: source,
        boundary_inflow: [None; 6],
    }
}

#[test]
fn graded_two_cell_mass_and_enthalpy_match_independent_implicit_balances_on_every_axis() {
    // V=[6,12], dt*Q=3. Density: rho0=(12+30)/9=14/3,
    // rho1=(96+3*rho0)/15=22/3. Energy density solves the SAME
    // volume operator: E0=(240+3000)/9=360, E1=(7680+3*360)/15=584.
    for axis in 0..3 {
        let mut widths = [vec![2.0], vec![3.0], vec![2.0]];
        widths[axis] = vec![1.0, 2.0];
        widths[(axis + 1) % 3] = vec![2.0];
        widths[(axis + 2) % 3] = vec![3.0];
        let domain = VoxelDomain::graded(widths).unwrap();
        let mut faces = [FlowFace::Wall; 6];
        faces[2 * axis] = FlowFace::Fixed;
        faces[2 * axis + 1] = FlowFace::Fixed;
        let field =
            FlowField::from_face_velocity(&domain, faces, |a, _| if a == axis { 1.0 } else { 0.0 })
                .unwrap();
        let mass = [12.0, 96.0];
        let h = [20.0, 80.0];
        let source = [MaterialVolumeSource::None; 2];
        let mut data = input(&mass, &h, &source);
        data.boundary_inflow[2 * axis] = Some(AdvectedMaterial {
            density_kg_m3: 10.0,
            specific_enthalpy_j_kg: 100.0,
        });
        let next = field
            .project_and_advect_material(
                &domain,
                faces,
                data,
                controls(0.5),
                &CancelGate::new_clock_free(),
            )
            .unwrap();
        close(next.mass_kg[0], 28.0);
        close(next.mass_kg[1], 88.0);
        close(next.specific_enthalpy_j_kg[0], 540.0 / 7.0);
        close(next.specific_enthalpy_j_kg[1], 876.0 / 11.0);
        close(next.report.incoming_mass_kg, 30.0);
        close(next.report.outgoing_mass_kg, 22.0);
        close(next.report.mass_after_kg - next.report.mass_before_kg, 8.0);
        close(next.report.enthalpy_change_j, 1248.0);
        close(next.report.carried_enthalpy_j, 1248.0);
        assert_eq!(next.report.mass_routes, 3);
        assert_eq!(field, next.projection.field);
    }
}

#[test]
fn volume_injection_and_extraction_move_actual_material_through_a_sealed_component() {
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let faces = [FlowFace::Wall; 6];
    let field = FlowField::quiescent(&domain);
    let mass = [2.0, 4.0];
    let h = [10.0, 50.0];
    let source = [
        MaterialVolumeSource::Inflow {
            volume_rate_m3_s: 1.0,
            material: AdvectedMaterial {
                density_kg_m3: 8.0,
                specific_enthalpy_j_kg: 100.0,
            },
        },
        MaterialVolumeSource::Outflow {
            volume_rate_m3_s: 1.0,
        },
    ];
    let next = field
        .project_and_advect_material(
            &domain,
            faces,
            input(&mass, &h, &source),
            controls(0.5),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    close(next.projection.field.flux_x(1, 0, 0), 1.0);
    close(next.mass_kg[0], 4.0);
    close(next.mass_kg[1], 4.0);
    close(next.specific_enthalpy_j_kg[0], 70.0);
    close(next.specific_enthalpy_j_kg[1], 170.0 / 3.0);
    close(next.report.incoming_mass_kg, 4.0);
    close(next.report.outgoing_mass_kg, 2.0);
    close(next.report.enthalpy_change_j, 860.0 / 3.0);
    close(next.report.carried_enthalpy_j, 860.0 / 3.0);
    assert!(next.projection.report.projection_loss_j > 0.0);
    // The pressure projection's numerical loss has not entered material h.
    close(
        next.mass_kg
            .iter()
            .zip(&next.specific_enthalpy_j_kg)
            .map(|(m, h)| m * h)
            .sum::<f64>(),
        1520.0 / 3.0,
    );
    assert_eq!(field, FlowField::quiescent(&domain));
    assert_eq!(mass, [2.0, 4.0]);
    assert_eq!(h, [10.0, 50.0]);
}

#[test]
fn uniform_density_and_large_enthalpy_datums_preserve_the_same_physical_mixing() {
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let faces = [
        FlowFace::Fixed,
        FlowFace::Fixed,
        FlowFace::Wall,
        FlowFace::Wall,
        FlowFace::Wall,
        FlowFace::Wall,
    ];
    let field =
        FlowField::from_face_velocity(&domain, faces, |axis, _| if axis == 0 { 1.0 } else { 0.0 })
            .unwrap();
    let mass = [2.0; 2];
    let source = [MaterialVolumeSource::None; 2];
    for datum in [0.0, 2.0_f64.powi(50), -2.0_f64.powi(50)] {
        let h = [datum + 4096.0, datum + 12_288.0];
        let mut data = input(&mass, &h, &source);
        data.boundary_inflow[0] = Some(AdvectedMaterial {
            density_kg_m3: 2.0,
            specific_enthalpy_j_kg: datum,
        });
        let next = field
            .project_and_advect_material(
                &domain,
                faces,
                data,
                controls(1.0),
                &CancelGate::new_clock_free(),
            )
            .unwrap();
        assert_eq!(next.mass_kg, mass);
        assert_eq!(next.specific_enthalpy_j_kg[0] - datum, 2048.0);
        assert_eq!(next.specific_enthalpy_j_kg[1] - datum, 7168.0);
        assert_eq!(next.report.enthalpy_change_j, -14_336.0);
        assert_eq!(next.report.carried_enthalpy_j, -14_336.0);
    }
    let h = [12345.0; 2];
    let mut data = input(&mass, &h, &source);
    data.boundary_inflow[0] = Some(AdvectedMaterial {
        density_kg_m3: 2.0,
        specific_enthalpy_j_kg: h[0],
    });
    let next = field
        .project_and_advect_material(
            &domain,
            faces,
            data,
            controls(1.0),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert_eq!(next.mass_kg, mass);
    assert_eq!(next.specific_enthalpy_j_kg, h);
}

fn circulation(domain: &VoxelDomain) -> FlowField {
    FlowField::from_face_velocity(domain, [FlowFace::Wall; 6], |axis, p| match axis {
        0 => {
            if p[1] < 1.0 {
                1.0
            } else {
                -1.0
            }
        }
        1 => {
            if p[0] < 1.0 {
                -1.0
            } else {
                1.0
            }
        }
        _ => 0.0,
    })
    .unwrap()
}

#[test]
fn continuing_closed_circulation_preserves_original_mass_and_enthalpy() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let faces = [FlowFace::Wall; 6];
    let mut field = circulation(&domain);
    let mut mass = vec![1.0, 2.0, 3.0, 4.0];
    let mut h = vec![10.0, 30.0, 20.0, 40.0];
    let source = [MaterialVolumeSource::None; 4];
    for _ in 0..32 {
        let next = field
            .project_and_advect_material(
                &domain,
                faces,
                input(&mass, &h, &source),
                controls(0.25),
                &CancelGate::new_clock_free(),
            )
            .unwrap();
        close(next.report.mass_before_kg, 10.0);
        close(next.report.mass_after_kg, 10.0);
        close(
            next.mass_kg
                .iter()
                .zip(&next.specific_enthalpy_j_kg)
                .map(|(m, h)| m * h)
                .sum::<f64>(),
            290.0,
        );
        assert_eq!(next.report.incoming_mass_kg, 0.0);
        assert_eq!(next.report.outgoing_mass_kg, 0.0);
        assert!(next.mass_kg.iter().all(|rho| (1.0..=4.0).contains(rho)));
        assert!(
            next.specific_enthalpy_j_kg
                .iter()
                .all(|h| (10.0..=40.0).contains(h))
        );
        field = next.projection.field;
        mass = next.mass_kg;
        h = next.specific_enthalpy_j_kg;
    }
    assert!((mass[0] - 2.5).abs() < 0.01);
}

#[test]
fn independent_mass_readback_rejects_a_loose_transport_solve() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let field = circulation(&domain);
    let original = field.clone();
    let mass = [1.0, 2.0, 3.0, 4.0];
    let h = [10.0, 30.0, 20.0, 40.0];
    let source = [MaterialVolumeSource::None; 4];
    let error = field
        .project_and_advect_material(
            &domain,
            [FlowFace::Wall; 6],
            input(&mass, &h, &source),
            MaterialAdvectionControls {
                relative_tolerance: 0.9,
                maximum_mass_defect_kg: 1e-12,
                ..controls(0.1)
            },
            &CancelGate::new_clock_free(),
        )
        .unwrap_err();
    assert!(matches!(
        error,
        ChtError::InvalidInput {
            field: "flow.material_advection",
            ..
        }
    ));
    assert_eq!(field, original);
}

#[test]
fn invalid_inventories_missing_inflow_budgets_and_cancellation_return_no_candidate() {
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let faces = [
        FlowFace::Fixed,
        FlowFace::Fixed,
        FlowFace::Wall,
        FlowFace::Wall,
        FlowFace::Wall,
        FlowFace::Wall,
    ];
    let field =
        FlowField::from_face_velocity(&domain, faces, |axis, _| if axis == 0 { 1.0 } else { 0.0 })
            .unwrap();
    let original = field.clone();
    let mass = [2.0, 4.0];
    let h = [10.0, 20.0];
    let sources = [MaterialVolumeSource::None; 2];
    let data = input(&mass, &h, &sources);
    assert!(
        field
            .project_and_advect_material(
                &domain,
                faces,
                data,
                controls(0.5),
                &CancelGate::new_clock_free()
            )
            .is_err()
    );
    let mut data = data;
    data.boundary_inflow[0] = Some(AdvectedMaterial {
        density_kg_m3: 2.0,
        specific_enthalpy_j_kg: 10.0,
    });
    for bad in [
        [0.0, 4.0],
        [-1.0, 4.0],
        [f64::NAN, 4.0],
        [f64::INFINITY, 4.0],
    ] {
        assert!(
            field
                .project_and_advect_material(
                    &domain,
                    faces,
                    MaterialAdvectionInput {
                        mass_kg: &bad,
                        ..data
                    },
                    controls(0.5),
                    &CancelGate::new_clock_free()
                )
                .is_err()
        );
    }
    for budget in [
        MaterialAdvectionControls {
            max_cells: 1,
            ..controls(0.5)
        },
        MaterialAdvectionControls {
            max_connections: 2,
            ..controls(0.5)
        },
        MaterialAdvectionControls {
            max_iterations: 0,
            ..controls(0.5)
        },
        MaterialAdvectionControls {
            maximum_enthalpy_defect_j: f64::NAN,
            ..controls(0.5)
        },
    ] {
        assert!(
            field
                .project_and_advect_material(
                    &domain,
                    faces,
                    data,
                    budget,
                    &CancelGate::new_clock_free()
                )
                .is_err()
        );
    }
    let gate = CancelGate::new_clock_free();
    gate.request();
    assert!(matches!(
        field.project_and_advect_material(&domain, faces, data, controls(0.5), &gate),
        Err(ChtError::Cancelled)
    ));
    assert_eq!(field, original);
    assert_eq!(mass, [2.0, 4.0]);
    assert_eq!(h, [10.0, 20.0]);
}
