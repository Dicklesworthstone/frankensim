use super::*;
use crate::conjugate::domain::Voxel;
use crate::conjugate::flow::{
    MaterialAdvectionControls, MaterialAdvectionInput, MaterialVolumeSource,
};

fn force_controls() -> DensityForceControls {
    DensityForceControls {
        max_cells: 100,
        maximum_impulse_defect_kg_m_s: 1e-10,
        maximum_energy_defect_j: 1e-10,
    }
}

fn projection_controls() -> DensityProjectionControls {
    DensityProjectionControls {
        duration_s: 0.4,
        relative_tolerance: 1e-12,
        maximum_volume_defect_m3_s: 1e-11,
        maximum_energy_defect_j: 1e-10,
    }
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 3e-11 * expected.abs().max(1.0),
        "{actual:.17e} != {expected:.17e}"
    );
}

fn input<'a>(
    density: &'a [f64],
    force: &'a [[f64; 3]],
    acceleration: [f64; 3],
) -> DensityForceInput<'a> {
    DensityForceInput {
        density_kg_m3: density,
        acceleration_m_s2: acceleration,
        force_density_n_m3: force,
    }
}

fn profile(values: &[Vec<Option<f64>>; 3]) -> StaggeredBoundaryPressure<'_> {
    StaggeredBoundaryPressure {
        x_pa: &values[0],
        y_pa: &values[1],
        z_pa: &values[2],
    }
}

#[test]
fn reservoir_head_and_body_force_accelerate_actual_graded_inertia() {
    let gate = CancelGate::new_clock_free();
    for axis in 0..3 {
        let mut widths = std::array::from_fn(|_| vec![1.0]);
        widths[axis] = vec![1.0, 2.0];
        widths[(axis + 1) % 3] = vec![2.0];
        widths[(axis + 2) % 3] = vec![3.0];
        let domain = VoxelDomain::graded(widths).unwrap();
        let initial = FlowField::quiescent(&domain);
        let mut faces = [FlowFace::Wall; 6];
        faces[2 * axis] = FlowFace::Free;
        faces[2 * axis + 1] = FlowFace::Free;
        let mut acceleration = [0.0; 3];
        acceleration[axis] = 0.75;
        for difference in [-16.0, 16.0] {
            let mut baseline: Option<ForcedDensityProjection> = None;
            for datum in [0.0, 1e12] {
                let mut pressure = initial.face_arrays().map(|q| vec![None; q.len()]);
                pressure[axis][initial.slot(&domain, 0, 2 * axis).0] = Some(datum + difference);
                pressure[axis][initial.slot(&domain, 1, 2 * axis + 1).0] = Some(datum);
                let result = initial
                    .force_and_project_density_with_pressure_profile(
                        &domain,
                        faces,
                        profile(&pressure),
                        input(&[2.0, 5.0], &[[0.0; 3]; 2], acceleration),
                        &[0.0; 2],
                        projection_controls(),
                        force_controls(),
                        &gate,
                    )
                    .unwrap();
                // A=6; rho*h sums to 12 kg/m^2, so the entire column receives
                // du=dt*(a+delta_p/12). The three dual masses sum to 72 kg.
                let dt = projection_controls().duration_s;
                let speed = dt * (0.75 + difference / 12.0);
                for c in 0..2 {
                    close(
                        result.projection.field.outward(&domain, c, 2 * axis),
                        -6.0 * speed,
                    );
                    close(
                        result.projection.field.outward(&domain, c, 2 * axis + 1),
                        6.0 * speed,
                    );
                }
                close(result.report.momentum_change_kg_m_s[axis], 72.0 * speed);
                close(
                    result.report.pressure_impulse_kg_m_s[axis],
                    dt * 6.0 * difference,
                );
                close(result.report.kinetic_energy_change_j, 36.0 * speed * speed);
                close(
                    result.projection.report.open_boundary_pressure_work_j,
                    -dt * 6.0 * speed * difference,
                );
                close(result.report.energy_defect_j, 0.0);
                if let Some(original) = &baseline {
                    // The fractional cell head is not exactly representable
                    // after adding 1e12 Pa. Fresh impulses must use the original
                    // component-relative solution, without relaxed budgets.
                    assert_eq!(result.projection.field, original.projection.field);
                    assert_eq!(result.projection.report, original.projection.report);
                    assert_eq!(result.report, original.report);
                } else {
                    baseline = Some(result);
                }
            }
        }
    }
}

#[test]
fn zero_pressure_profile_preserves_force_bits_and_invalid_profiles_refuse() {
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let initial = FlowField::quiescent(&domain);
    let faces = [
        FlowFace::Free,
        FlowFace::Free,
        FlowFace::Wall,
        FlowFace::Wall,
        FlowFace::Wall,
        FlowFace::Wall,
    ];
    let gate = CancelGate::new_clock_free();
    let mut pressure = initial.face_arrays().map(|q| vec![None; q.len()]);
    let run = |pressure: StaggeredBoundaryPressure<'_>, gate: &CancelGate| {
        initial.force_and_project_density_with_pressure_profile(
            &domain,
            faces,
            pressure,
            input(&[2.0, 5.0], &[[0.0; 3]; 2], [0.25, 0.0, 0.0]),
            &[0.0; 2],
            projection_controls(),
            force_controls(),
            gate,
        )
    };
    let original = initial
        .force_and_project_density(
            &domain,
            faces,
            input(&[2.0, 5.0], &[[0.0; 3]; 2], [0.25, 0.0, 0.0]),
            &[0.0; 2],
            projection_controls(),
            force_controls(),
            &gate,
        )
        .unwrap();
    assert_eq!(run(profile(&pressure), &gate).unwrap(), original);
    pressure[0][1] = Some(0.0); // Interior faces cannot prescribe reservoir pressure.
    assert!(run(profile(&pressure), &gate).is_err());
    pressure[0][1] = None;
    pressure[0][0] = Some(f64::NAN);
    assert!(run(profile(&pressure), &gate).is_err());
    pressure[0][0] = Some(1.0);
    assert!(
        run(
            StaggeredBoundaryPressure {
                x_pa: &[],
                ..profile(&pressure)
            },
            &gate
        )
        .is_err()
    );
    gate.request();
    assert_eq!(run(profile(&pressure), &gate), Err(ChtError::Cancelled));
    assert_eq!(initial, FlowField::quiescent(&domain));
}

#[test]
fn graded_sharp_stratification_balances_weight_and_endpoint_work_on_every_axis() {
    for axis in 0..3 {
        let mut widths = [vec![2.0], vec![2.0], vec![2.0]];
        widths[axis] = vec![1.0, 2.0];
        widths[(axis + 1) % 3] = vec![3.0];
        let domain = VoxelDomain::graded(widths).unwrap();
        let initial = FlowField::quiescent(&domain);
        for free_side in [Some(0), Some(1), None] {
            let mut faces = [FlowFace::Wall; 6];
            if let Some(side) = free_side {
                faces[2 * axis + side] = FlowFace::Free;
            }
            let mut acceleration = [0.0; 3];
            acceleration[axis] = if free_side == Some(0) { 3.25 } else { -3.25 };
            let result = initial
                .force_and_project_density(
                    &domain,
                    faces,
                    input(&[2.0, 5.0], &[[0.0; 3]; 2], acceleration),
                    &[0.0; 2],
                    projection_controls(),
                    force_controls(),
                    &CancelGate::new_clock_free(),
                )
                .unwrap();
            // The pressure increment across the internal face is six times
            // g: rho0*h0/2 + rho1*h1/2 = 1+5. Free pressure anchors at its
            // boundary half-cell; sealed pressure pins cell zero instead.
            let expected = match free_side {
                Some(0) => [3.25, 7.0 * 3.25],
                Some(1) => [11.0 * 3.25, 5.0 * 3.25],
                _ => [0.0, -6.0 * 3.25],
            };
            for (actual, expected) in result
                .projection
                .pressure_correction_pa
                .iter()
                .zip(expected)
            {
                close(*actual, expected);
            }
            for fluxes in result.projection.field.face_arrays() {
                for &q in fluxes {
                    close(q, 0.0);
                }
            }
            let weight_impulse = 0.4 * acceleration[axis] * (6.0 * 2.0 + 12.0 * 5.0);
            close(result.forcing.applied_impulse_kg_m_s[axis], weight_impulse);
            close(
                result.report.pressure_and_constraint_impulse_kg_m_s[axis],
                -weight_impulse,
            );
            close(result.report.momentum_change_kg_m_s[axis], 0.0);
            assert!(result.forcing.applied_work_j > 0.0);
            assert!(result.projection.report.projection_loss_j > 0.0);
            close(result.report.applied_endpoint_work_j, 0.0);
            close(result.report.constraint_endpoint_work_j, 0.0);
            close(
                result.report.force_pressure_split_work_j,
                -result.projection.report.projection_loss_j,
            );
            close(result.report.kinetic_energy_change_j, 0.0);
            close(result.report.energy_defect_j, 0.0);
            assert_eq!(initial, FlowField::quiescent(&domain));
        }
    }
}

fn column_contrast(force_axis: usize) -> (VoxelDomain, Vec<f64>, [f64; 3]) {
    let across = (force_axis + 2) % 3;
    let mut dims = [1; 3];
    dims[force_axis] = 2;
    dims[across] = 2;
    let domain = VoxelDomain::new(dims[0], dims[1], dims[2], 1.0).unwrap();
    let density = (0..4)
        .map(|c| {
            if domain.coords(c)[across] == 0 {
                1.0
            } else {
                3.0
            }
        })
        .collect();
    let mut acceleration = [0.0; 3];
    acceleration[force_axis] = -2.0;
    (domain, density, acceleration)
}

#[test]
fn horizontal_density_contrast_drives_exact_baroclinic_circulation_in_all_planes() {
    for axis in 0..3 {
        let (domain, density, acceleration) = column_contrast(axis);
        let across = (axis + 2) % 3;
        let result = FlowField::quiescent(&domain)
            .force_and_project_density(
                &domain,
                [FlowFace::Wall; 6],
                input(&density, &[[0.0; 3]; 4], acceleration),
                &[0.0; 4],
                projection_controls(),
                force_controls(),
                &CancelGate::new_clock_free(),
            )
            .unwrap();
        // The single loop has masses 1+3+2+2=8. Its net gravitational
        // circulation impulse is dt*g*(3-1), hence speed dt*g/4 = 0.2.
        let speed = 0.4 * 2.0 / 4.0;
        for c in 0..4 {
            let position = domain.coords(c);
            if position[axis] == 0 {
                let expected = if position[across] == 0 { speed } else { -speed };
                close(
                    result.projection.field.outward(&domain, c, 2 * axis + 1),
                    expected,
                );
            }
            if position[across] == 0 {
                let expected = if position[axis] == 0 { -speed } else { speed };
                close(
                    result.projection.field.outward(&domain, c, 2 * across + 1),
                    expected,
                );
            }
            close(result.projection.field.net_outflow(&domain, c), 0.0);
        }
        close(
            result.report.kinetic_energy_change_j,
            0.5 * 8.0 * speed * speed,
        );
        close(result.report.applied_endpoint_work_j, 0.4 * 2.0 * speed);
        close(result.report.energy_defect_j, 0.0);
    }
}

#[test]
fn graded_cell_force_uses_physical_half_volumes_and_actual_endpoint_inertia() {
    let domain = VoxelDomain::graded([vec![1.0, 2.0], vec![2.0], vec![3.0]]).unwrap();
    let faces = [
        FlowFace::Free,
        FlowFace::Free,
        FlowFace::Wall,
        FlowFace::Wall,
        FlowFace::Wall,
        FlowFace::Wall,
    ];
    let initial =
        FlowField::from_face_velocity(&domain, faces, |axis, _| if axis == 0 { 0.25 } else { 0.0 })
            .unwrap();
    let result = initial
        .predict_density_body_force(
            &domain,
            faces,
            input(&[2.0, 5.0], &[[4.0, 0.0, 0.0], [5.0, 0.0, 0.0]], [0.0; 3]),
            0.5,
            force_controls(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    // Three independent dual masses: 6,36,30 kg; half-volume force
    // impulses: 6,21,15 kg m/s. Both boundary faces use world +x signs.
    close(result.field.flux_x(0, 0, 0) / 6.0, 1.25);
    close(result.field.flux_x(1, 0, 0) / 6.0, 0.25 + 7.0 / 12.0);
    close(result.field.flux_x(2, 0, 0) / 6.0, 0.75);
    close(result.report.applied_impulse_kg_m_s[0], 42.0);
    close(result.report.momentum_change_kg_m_s[0], 42.0);
    close(result.report.constraint_impulse_kg_m_s[0], 0.0);
    close(result.report.kinetic_energy_before_j, 2.25);
    close(result.report.kinetic_energy_change_j, 23.375);
    close(result.report.applied_work_j, 23.375);
    close(result.report.energy_defect_j, 0.0);
}

#[test]
fn fixed_flow_controller_reacts_with_equal_opposite_impulse_and_work() {
    let domain = VoxelDomain::graded([vec![2.0], vec![3.0], vec![4.0]]).unwrap();
    let faces = [
        FlowFace::Fixed,
        FlowFace::Fixed,
        FlowFace::Wall,
        FlowFace::Wall,
        FlowFace::Wall,
        FlowFace::Wall,
    ];
    let initial =
        FlowField::from_face_velocity(&domain, faces, |axis, _| if axis == 0 { 2.0 } else { 0.0 })
            .unwrap();
    let mut ctl = projection_controls();
    ctl.duration_s = 0.5;
    let result = initial
        .force_and_project_density(
            &domain,
            faces,
            input(&[5.0], &[[0.0; 3]], [3.0, -1.0, 2.0]),
            &[0.0],
            ctl,
            force_controls(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert_eq!(result.projection.field, initial);
    for (axis, expected) in [180.0, -60.0, 120.0].into_iter().enumerate() {
        close(result.forcing.applied_impulse_kg_m_s[axis], expected);
        close(result.forcing.constraint_impulse_kg_m_s[axis], -expected);
        close(result.report.momentum_change_kg_m_s[axis], 0.0);
    }
    close(result.forcing.applied_work_j, 360.0);
    close(result.forcing.constraint_work_j, -360.0);
    close(result.report.applied_endpoint_work_j, 360.0);
    close(result.report.constraint_endpoint_work_j, -360.0);
    close(result.report.force_pressure_split_work_j, 0.0);
    close(result.report.energy_defect_j, 0.0);
}

#[test]
fn walls_and_both_solid_interface_orientations_retain_every_fluid_half_dual() {
    let mut domain = VoxelDomain::new(3, 1, 1, 1.0).unwrap();
    domain.set(1, 0, 0, Voxel::Solid(0));
    let field = FlowField::quiescent(&domain);
    let result = field
        .predict_density_body_force(
            &domain,
            [FlowFace::Wall; 6],
            input(
                &[2.0, f64::NAN, 3.0],
                &[[1.0, 2.0, 3.0], [0.0; 3], [4.0, 5.0, 6.0]],
                [1.0, -2.0, 3.0],
            ),
            0.5,
            force_controls(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert_eq!(result.field, field);
    for (axis, expected) in [5.0, -1.5, 12.0].into_iter().enumerate() {
        close(result.report.applied_impulse_kg_m_s[axis], expected);
        close(result.report.constraint_impulse_kg_m_s[axis], -expected);
        close(result.report.momentum_change_kg_m_s[axis], 0.0);
    }
    assert_eq!(result.report.applied_work_j, 0.0);
    assert_eq!(result.report.constraint_work_j, 0.0);
}

#[test]
fn zero_force_preserves_the_native_density_projection_exactly() {
    let domain = VoxelDomain::graded([vec![1.0, 2.0], vec![2.0], vec![3.0]]).unwrap();
    let mut faces = [FlowFace::Wall; 6];
    faces[1] = FlowFace::Free;
    let field = FlowField::quiescent(&domain);
    let gate = CancelGate::new_clock_free();
    let data = input(&[2.0, 5.0], &[[0.0; 3]; 2], [0.0; 3]);
    let predicted = field
        .predict_density_body_force(&domain, faces, data, 0.4, force_controls(), &gate)
        .unwrap();
    assert_eq!(predicted.field, field);
    let native = field
        .project_density_volume_source(
            &domain,
            faces,
            data.density_kg_m3,
            &[0.125, 0.25],
            projection_controls(),
            &gate,
        )
        .unwrap();
    let composed = field
        .force_and_project_density(
            &domain,
            faces,
            data,
            &[0.125, 0.25],
            projection_controls(),
            force_controls(),
            &gate,
        )
        .unwrap();
    assert_eq!(composed.projection, native);
    assert_eq!(composed.report.applied_endpoint_work_j, 0.0);
    assert_eq!(composed.report.force_pressure_split_work_j, 0.0);
}

#[test]
fn force_predictor_feeds_the_actual_conservative_material_owner_once() {
    let (domain, density, acceleration) = column_contrast(1);
    let field = FlowField::quiescent(&domain);
    let faces = [FlowFace::Wall; 6];
    let gate = CancelGate::new_clock_free();
    let force = [[0.0; 3]; 4];
    let data = input(&density, &force, acceleration);
    let predicted = field
        .predict_density_body_force(&domain, faces, data, 0.4, force_controls(), &gate)
        .unwrap();
    let controls = MaterialAdvectionControls {
        projection: projection_controls(),
        max_cells: 100,
        max_connections: 100,
        max_iterations: 100,
        relative_tolerance: 1e-12,
        maximum_mass_defect_kg: 1e-10,
        maximum_enthalpy_defect_j: 1e-10,
    };
    let transported = predicted
        .field
        .project_and_advect_material(
            &domain,
            faces,
            MaterialAdvectionInput {
                mass_kg: &density,
                specific_enthalpy_j_kg: &[11.0; 4],
                volume_sources: &[MaterialVolumeSource::None; 4],
                boundary_inflow: [None; 6],
            },
            controls,
            &gate,
        )
        .unwrap();
    let direct = field
        .force_and_project_density(
            &domain,
            faces,
            data,
            &[0.0; 4],
            projection_controls(),
            force_controls(),
            &gate,
        )
        .unwrap();
    assert_eq!(transported.projection, direct.projection);
    close(transported.mass_kg.iter().sum(), 8.0);
    assert!(
        transported
            .mass_kg
            .iter()
            .zip(&density)
            .any(|(after, before)| (after - before).abs() > 1e-3)
    );
    for h in transported.specific_enthalpy_j_kg {
        close(h, 11.0);
    }
    close(transported.report.enthalpy_change_j, 0.0);
}

#[test]
fn refusal_cancellation_and_independent_readback_leave_the_original_unchanged() {
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let faces = [FlowFace::Free; 6];
    let field =
        FlowField::from_face_velocity(
            &domain,
            faces,
            |axis, _| if axis == 0 { 2f64.powi(54) } else { 0.0 },
        )
        .unwrap();
    let original = field.clone();
    let data = input(&[1.0], &[[0.0; 3]], [1.0, 0.0, 0.0]);
    let gate = CancelGate::new_clock_free();
    // A rounded-away physical impulse must not pass merely because the
    // returned enormous velocity is finite or because its net outflow is zero.
    assert!(matches!(
        field.predict_density_body_force(&domain, faces, data, 0.5, force_controls(), &gate),
        Err(ChtError::InvalidInput {
            field: "flow.force_impulse_budget",
            ..
        })
    ));
    let mut controls = force_controls();
    controls.max_cells = 0;
    assert!(matches!(
        field.predict_density_body_force(&domain, faces, data, 0.5, controls, &gate),
        Err(ChtError::InvalidInput {
            field: "flow.force_cell_budget",
            ..
        })
    ));
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            field
                .predict_density_body_force(
                    &domain,
                    faces,
                    input(&[bad], &[[0.0; 3]], [0.0; 3]),
                    0.5,
                    force_controls(),
                    &gate
                )
                .is_err()
        );
    }
    assert!(
        field
            .predict_density_body_force(
                &domain,
                faces,
                input(&[1.0], &[], [0.0; 3]),
                0.5,
                force_controls(),
                &gate
            )
            .is_err()
    );
    gate.request();
    assert_eq!(
        field.predict_density_body_force(&domain, faces, data, 0.5, force_controls(), &gate),
        Err(ChtError::Cancelled)
    );
    assert_eq!(field, original);

    // Pressure readback is independent of the pressure-corrected flux. An
    // edited pressure cannot inherit the already accepted native report.
    let quiet = FlowField::quiescent(&domain);
    let gate = CancelGate::new_clock_free();
    let prediction = quiet
        .predict_density_body_force(&domain, faces, data, 0.4, force_controls(), &gate)
        .unwrap();
    let mut projected = prediction
        .field
        .project_density_volume_source(&domain, faces, &[1.0], &[0.0], projection_controls(), &gate)
        .unwrap();
    projected.pressure_correction_pa[0] += 1.0;
    assert!(matches!(
        audit_combined(
            &quiet,
            &prediction,
            &projected,
            &domain,
            faces,
            None,
            None,
            data,
            0.4,
            force_controls(),
            &gate
        ),
        Err(ChtError::InvalidInput {
            field: "flow.combined_impulse_budget",
            ..
        })
    ));
}
