use super::*;
use crate::conjugate::domain::Voxel;

fn controls() -> DensityProjectionControls {
    DensityProjectionControls {
        duration_s: 0.4,
        relative_tolerance: 1e-12,
        maximum_volume_defect_m3_s: 1e-11,
        maximum_energy_defect_j: 1e-10,
    }
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 2e-11 * expected.abs().max(1.0),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn graded_density_jump_drives_actual_pressure_flux_and_work_on_all_axes() {
    // Independent two-cell solution: A=6, half-widths=.5 and 1,
    // M_internal=6*(2*.5+5*1)=36, M_out=6*5*1=30.
    // Conductances .4 and .48 require p1=1/.48, p0=p1+.8/.4.
    for axis in 0..3 {
        let mut widths = [vec![2.0], vec![3.0], vec![2.0]];
        widths[axis] = vec![1.0, 2.0];
        widths[(axis + 1) % 3] = vec![2.0];
        widths[(axis + 2) % 3] = vec![3.0];
        let domain = VoxelDomain::graded(widths).unwrap();
        let mut faces = [FlowFace::Wall; 6];
        faces[2 * axis + 1] = FlowFace::Free;
        let initial = FlowField::quiescent(&domain);
        let result = initial
            .project_density_volume_source(
                &domain,
                faces,
                &[2.0, 5.0],
                &[0.8, 0.2],
                controls(),
                &CancelGate::new(),
            )
            .unwrap();
        close(result.field.outward(&domain, 0, 2 * axis + 1), 0.8);
        close(result.field.outward(&domain, 1, 2 * axis + 1), 1.0);
        close(result.field.net_outflow(&domain, 0), 0.8);
        close(result.field.net_outflow(&domain, 1), 0.2);
        close(result.pressure_correction_pa[0], 49.0 / 12.0);
        close(result.pressure_correction_pa[1], 25.0 / 12.0);
        close(result.report.kinetic_energy_after_j, 221.0 / 300.0);
        close(result.report.kinetic_energy_change_j, 221.0 / 300.0);
        close(result.report.projection_loss_j, 221.0 / 300.0);
        close(result.report.pressure_source_work_j, 221.0 / 150.0);
        close(result.report.fixed_boundary_pressure_work_j, 0.0);
        assert!(result.report.energy_defect_j.abs() < 1e-12);
        assert_eq!(initial, FlowField::quiescent(&domain));
    }
}

#[test]
fn density_contrast_changes_actual_flow_partition_between_competing_outlets() {
    // Two neighbouring cells, both open at y-max, with one source in cell 0.
    // For rho=[1,R], dt=1, A=h=1, direct 2x2 elimination gives
    // p=[(1+2R)/(4+4R), R/(4+4R)] and outlet-1 Q=1/(2+2R).
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let mut faces = [FlowFace::Wall; 6];
    faces[3] = FlowFace::Free;
    for ratio in [1.0, 9.0, 999.0] {
        let result = FlowField::quiescent(&domain)
            .project_density_volume_source(
                &domain,
                faces,
                &[1.0, ratio],
                &[1.0, 0.0],
                DensityProjectionControls {
                    duration_s: 1.0,
                    ..controls()
                },
                &CancelGate::new(),
            )
            .unwrap();
        let remote_outflow = 1.0 / (2.0 + 2.0 * ratio);
        close(result.field.flux_x(1, 0, 0), remote_outflow);
        close(result.field.flux_y(0, 1, 0), 1.0 - remote_outflow);
        close(result.field.flux_y(1, 1, 0), remote_outflow);
        close(
            result.pressure_correction_pa[0],
            (1.0 + 2.0 * ratio) / (4.0 + 4.0 * ratio),
        );
        close(
            result.pressure_correction_pa[1],
            ratio / (4.0 + 4.0 * ratio),
        );
        close(
            result.report.kinetic_energy_after_j,
            0.5 * result.pressure_correction_pa[0],
        );
    }
}

#[test]
fn zero_source_projection_removes_only_the_density_weighted_gradient_energy() {
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let faces = [FlowFace::Wall; 6];
    let initial =
        FlowField::from_face_velocity(&domain, faces, |axis, _| if axis == 0 { 0.5 } else { 0.0 })
            .unwrap();
    let result = initial
        .project_density_volume_source(
            &domain,
            faces,
            &[2.0, 6.0],
            &[0.0, 0.0],
            controls(),
            &CancelGate::new(),
        )
        .unwrap();
    close(result.field.max_flux(), 0.0);
    close(result.report.kinetic_energy_before_j, 0.5);
    close(result.report.kinetic_energy_after_j, 0.0);
    close(result.report.kinetic_energy_change_j, -0.5);
    close(result.report.projection_loss_j, 0.5);
    assert_eq!(result.report.pressure_source_work_j, 0.0);
    assert_eq!(result.report.fixed_boundary_pressure_work_j, 0.0);
}

#[test]
fn homogeneous_zero_source_retains_original_projection_bits() {
    for domain in [
        VoxelDomain::new(4, 3, 2, 0.25).unwrap(),
        VoxelDomain::graded([
            vec![0.2, 0.7, 0.3, 0.5],
            vec![0.1, 0.4, 0.2],
            vec![0.3, 0.8],
        ])
        .unwrap(),
    ] {
        let faces = [
            FlowFace::Fixed,
            FlowFace::Free,
            FlowFace::Wall,
            FlowFace::Wall,
            FlowFace::Wall,
            FlowFace::Wall,
        ];
        let initial = FlowField::from_face_velocity(&domain, faces, |axis, p| {
            (axis as f64 + 0.4) * (0.1 + p[0] * p[1] + p[2])
        })
        .unwrap();
        let mut original = initial.clone();
        let native = original
            .project(&domain, faces, 1e-12, &CancelGate::new())
            .unwrap();
        let result = initial
            .project_density_volume_source(
                &domain,
                faces,
                &vec![997.0; domain.cell_count()],
                &vec![0.0; domain.cell_count()],
                DensityProjectionControls {
                    maximum_energy_defect_j: 1e-7,
                    ..controls()
                },
                &CancelGate::new(),
            )
            .unwrap();
        for axis in 0..3 {
            assert_eq!(
                original.face_arrays()[axis]
                    .iter()
                    .map(|q| q.to_bits())
                    .collect::<Vec<_>>(),
                result.field.face_arrays()[axis]
                    .iter()
                    .map(|q| q.to_bits())
                    .collect::<Vec<_>>()
            );
        }
        assert_eq!(native, result.report.projection);
    }
}

#[test]
fn held_inflow_pressure_work_closes_real_kinetic_change() {
    // Two unit cells, rho=[2,6], dt=.4. Initial fixed inlet Q=.3,
    // all other Q=0, sources=[.2,.5]: final Q=[.3,.5,1].
    // M_internal=4, M_outlet=3 => p=[12.5,7.5].
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let faces = [
        FlowFace::Fixed,
        FlowFace::Free,
        FlowFace::Wall,
        FlowFace::Wall,
        FlowFace::Wall,
        FlowFace::Wall,
    ];
    let initial = FlowField::from_face_velocity(&domain, faces, |axis, p| {
        if axis == 0 && p[0] == 0.0 { 0.3 } else { 0.0 }
    })
    .unwrap();
    let result = initial
        .project_density_volume_source(
            &domain,
            faces,
            &[2.0, 6.0],
            &[0.2, 0.5],
            controls(),
            &CancelGate::new(),
        )
        .unwrap();
    close(result.field.flux_x(0, 0, 0), 0.3);
    close(result.field.flux_x(1, 0, 0), 0.5);
    close(result.field.flux_x(2, 0, 0), 1.0);
    close(result.pressure_correction_pa[0], 12.5);
    close(result.pressure_correction_pa[1], 7.5);
    close(result.report.kinetic_energy_before_j, 0.045);
    close(result.report.kinetic_energy_after_j, 2.045);
    close(result.report.kinetic_energy_change_j, 2.0);
    close(result.report.projection_loss_j, 2.0);
    close(result.report.pressure_source_work_j, 2.5);
    close(result.report.fixed_boundary_pressure_work_j, -1.5);
    close(result.report.energy_defect_j, 0.0);
}

#[test]
fn balanced_sealed_source_and_sink_move_volume_without_losing_an_equation() {
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let initial = FlowField::quiescent(&domain);
    let result = initial
        .project_density_volume_source(
            &domain,
            [FlowFace::Wall; 6],
            &[2.0, 6.0],
            &[0.25, -0.25],
            controls(),
            &CancelGate::new(),
        )
        .unwrap();
    assert_eq!(result.report.projection.pinned_components, 1);
    close(result.field.flux_x(1, 0, 0), 0.25);
    close(result.pressure_correction_pa[0], 0.0);
    close(result.pressure_correction_pa[1], -2.5);
    close(result.report.kinetic_energy_after_j, 0.125);
    close(result.report.pressure_source_work_j, 0.25);
    for c in 0..2 {
        close(result.field.net_outflow(&domain, c), [0.25, -0.25][c]);
    }
}

#[test]
fn disconnected_components_cannot_hide_incompatible_volume_sources() {
    let domain = VoxelDomain::from_fn(3, 1, 1, 1.0, |p| {
        if p[0] > 1.0 && p[0] < 2.0 {
            Voxel::Solid(0)
        } else {
            Voxel::Fluid
        }
    })
    .unwrap();
    let initial = FlowField::quiescent(&domain);
    let error = initial
        .project_density_volume_source(
            &domain,
            [FlowFace::Wall; 6],
            &[2.0, 0.0, 5.0],
            &[0.25, 0.0, -0.25],
            controls(),
            &CancelGate::new(),
        )
        .unwrap_err();
    assert!(matches!(
        error,
        ChtError::InvalidInput {
            field: "flow.volume_source_compatibility",
            ..
        }
    ));
    assert_eq!(initial, FlowField::quiescent(&domain));
}

#[test]
fn source_can_balance_prescribed_outflow_without_a_free_boundary() {
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let faces = [
        FlowFace::Wall,
        FlowFace::Fixed,
        FlowFace::Wall,
        FlowFace::Wall,
        FlowFace::Wall,
        FlowFace::Wall,
    ];
    let initial =
        FlowField::from_face_velocity(&domain, faces, |axis, _| if axis == 0 { 0.5 } else { 0.0 })
            .unwrap();
    let result = initial
        .project_density_volume_source(
            &domain,
            faces,
            &[2.0],
            &[0.5],
            controls(),
            &CancelGate::new(),
        )
        .unwrap();
    assert_eq!(result.field, initial);
    assert_eq!(result.report.projection.pinned_components, 1);
    assert_eq!(result.report.projection.iterations, 0);
    assert_eq!(result.report.max_volume_defect_after_m3_s, 0.0);
}

#[test]
fn blocked_flux_invalid_density_solid_source_and_cancellation_refuse_unchanged() {
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let open = [
        FlowFace::Fixed,
        FlowFace::Free,
        FlowFace::Wall,
        FlowFace::Wall,
        FlowFace::Wall,
        FlowFace::Wall,
    ];
    let flowing =
        FlowField::from_face_velocity(&domain, open, |axis, _| if axis == 0 { 1.0 } else { 0.0 })
            .unwrap();
    let unchanged = flowing.clone();
    let error = flowing
        .project_density_volume_source(
            &domain,
            [FlowFace::Wall; 6],
            &[2.0, 6.0],
            &[0.0, 0.0],
            controls(),
            &CancelGate::new(),
        )
        .unwrap_err();
    assert!(matches!(
        error,
        ChtError::InvalidInput {
            field: "flow.blocked_face_flux",
            ..
        }
    ));
    for density in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            flowing
                .project_density_volume_source(
                    &domain,
                    open,
                    &[2.0, density],
                    &[0.0, 0.0],
                    controls(),
                    &CancelGate::new()
                )
                .is_err()
        );
    }
    let solid = VoxelDomain::from_fn(2, 1, 1, 1.0, |p| {
        if p[0] > 1.0 {
            Voxel::Solid(0)
        } else {
            Voxel::Fluid
        }
    })
    .unwrap();
    let solid_flow = FlowField::quiescent(&solid);
    assert!(matches!(
        solid_flow.project_density_volume_source(
            &solid,
            [FlowFace::Wall; 6],
            &[2.0, 0.0],
            &[0.0, 0.1],
            controls(),
            &CancelGate::new(),
        ),
        Err(ChtError::InvalidInput {
            field: "flow.solid_volume_source",
            ..
        })
    ));
    // Reusing a field against a newly solid cell must refuse its held flux.
    assert!(
        flowing
            .project_density_volume_source(
                &solid,
                open,
                &[2.0, 0.0],
                &[0.0, 0.0],
                controls(),
                &CancelGate::new()
            )
            .is_err()
    );
    let gate = CancelGate::new();
    gate.request();
    assert_eq!(
        flowing.project_density_volume_source(
            &domain,
            open,
            &[2.0, 6.0],
            &[0.0, 0.0],
            controls(),
            &gate
        ),
        Err(ChtError::Cancelled)
    );
    assert_eq!(flowing, unchanged);
}

#[test]
fn independent_volume_budget_rejects_a_loose_native_solve() {
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let initial = FlowField::quiescent(&domain);
    let faces = [
        FlowFace::Wall,
        FlowFace::Free,
        FlowFace::Wall,
        FlowFace::Wall,
        FlowFace::Wall,
        FlowFace::Wall,
    ];
    let error = initial
        .project_density_volume_source(
            &domain,
            faces,
            &[2.0, 6.0],
            &[0.2, 0.5],
            DensityProjectionControls {
                relative_tolerance: 2.0,
                ..controls()
            },
            &CancelGate::new(),
        )
        .unwrap_err();
    assert!(matches!(
        error,
        ChtError::InvalidInput {
            field: "flow.volume_defect_budget",
            ..
        }
    ));
    assert_eq!(initial, FlowField::quiescent(&domain));
}

#[test]
fn independent_energy_readback_refuses_inconsistent_pressure() {
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let initial = FlowField::quiescent(&domain);
    let faces = [
        FlowFace::Wall,
        FlowFace::Free,
        FlowFace::Wall,
        FlowFace::Wall,
        FlowFace::Wall,
        FlowFace::Wall,
    ];
    let result = initial
        .project_density_volume_source(
            &domain,
            faces,
            &[2.0, 6.0],
            &[0.2, 0.5],
            controls(),
            &CancelGate::new(),
        )
        .unwrap();
    let model = FrozenDensity {
        density_kg_m3: &[2.0, 6.0],
        target_outflow_m3_s: &[0.2, 0.5],
        duration_s: 0.4,
        uniform_density: None,
    };
    let pressure = result
        .pressure_correction_pa
        .iter()
        .map(|p| p + 1.0)
        .collect::<Vec<_>>();
    let error = audit(
        &initial,
        &result.field,
        &domain,
        faces,
        &model,
        &pressure,
        result.report.projection,
        controls(),
        &CancelGate::new(),
    )
    .unwrap_err();
    assert!(matches!(
        error,
        ChtError::InvalidInput {
            field: "flow.energy_defect_budget",
            ..
        }
    ));
}
