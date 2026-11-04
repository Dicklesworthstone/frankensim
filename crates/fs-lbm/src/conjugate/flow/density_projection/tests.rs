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
        boundary_pressure: BoundaryPressure::Planes([0.0; 6]),
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
        &ProjectionPotential {
            cells: vec![0, 1],
            values: Vec::new(),
            component: vec![0, 0],
            reference_pressure_pa: vec![0.0],
        },
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

#[test]
fn reservoir_pressure_difference_accelerates_graded_density_in_either_direction() {
    // A=6, h=[1,2], rho=[2,5]. The two boundary half-cells and the
    // shared dual give total series inertia sum(M)/A = 12 kg/m^2.
    // Therefore Q=dt*A*delta_p/12=.2*delta_p, K=loss=Q^2 and
    // reservoir work=-dt*Q*delta_p. These are pressure impulses, not
    // steady Poiseuille solutions (no viscous resistance is supplied).
    for axis in 0..3 {
        let mut widths = [vec![2.0], vec![3.0], vec![2.0]];
        widths[axis] = vec![1.0, 2.0];
        widths[(axis + 1) % 3] = vec![2.0];
        widths[(axis + 2) % 3] = vec![3.0];
        let domain = VoxelDomain::graded(widths).unwrap();
        let mut faces = [FlowFace::Wall; 6];
        faces[2 * axis] = FlowFace::Free;
        faces[2 * axis + 1] = FlowFace::Free;
        for difference in [-12.0, 12.0] {
            let mut pressures = [None; 6];
            pressures[2 * axis] = Some(80.0 + difference);
            pressures[2 * axis + 1] = Some(80.0);
            let initial = FlowField::quiescent(&domain);
            let result = initial
                .project_density_volume_source_with_boundary_pressure(
                    &domain,
                    faces,
                    pressures,
                    &[2.0, 5.0],
                    &[0.0, 0.0],
                    controls(),
                    &CancelGate::new(),
                )
                .unwrap();
            let flux = 0.2 * difference;
            close(result.field.outward(&domain, 0, 2 * axis), -flux);
            close(result.field.outward(&domain, 0, 2 * axis + 1), flux);
            close(result.field.outward(&domain, 1, 2 * axis + 1), flux);
            close(result.field.net_outflow(&domain, 0), 0.0);
            close(result.field.net_outflow(&domain, 1), 0.0);
            close(
                result.pressure_correction_pa[0],
                80.0 + 11.0 * difference / 12.0,
            );
            close(
                result.pressure_correction_pa[1],
                80.0 + 5.0 * difference / 12.0,
            );
            close(result.report.kinetic_energy_change_j, flux * flux);
            close(result.report.projection_loss_j, flux * flux);
            close(
                result.report.open_boundary_pressure_work_j,
                -0.4 * flux * difference,
            );
            close(result.report.pressure_source_work_j, 0.0);
            close(result.report.fixed_boundary_pressure_work_j, 0.0);
            close(result.report.energy_defect_j, 0.0);
            assert_eq!(result.report.projection.pinned_components, 0);
            assert_eq!(initial, FlowField::quiescent(&domain));
        }
    }
}

#[test]
fn large_common_reservoir_pressure_preserves_small_driving_flux_and_work() {
    let domain = VoxelDomain::graded([vec![1.0, 2.0], vec![2.0], vec![3.0]]).unwrap();
    let faces = [
        FlowFace::Free,
        FlowFace::Free,
        FlowFace::Wall,
        FlowFace::Wall,
        FlowFace::Wall,
        FlowFace::Wall,
    ];
    let initial = FlowField::quiescent(&domain);
    for density in [[2.0, 5.0], [5.0, 5.0]] {
        let project = |datum: f64, difference: f64| {
            initial
                .project_density_volume_source_with_boundary_pressure(
                    &domain,
                    faces,
                    [
                        Some(datum + difference),
                        Some(datum),
                        None,
                        None,
                        None,
                        None,
                    ],
                    &density,
                    &[0.0, 0.0],
                    controls(),
                    &CancelGate::new(),
                )
                .unwrap()
        };
        let baseline = project(0.0, 16.0);
        let shifted = project(1e12, 16.0);
        assert_eq!(baseline.field, shifted.field);
        assert_eq!(baseline.report, shifted.report);
        assert!(shifted.report.kinetic_energy_change_j > 0.0);
        let resting = project(1e12, 0.0);
        assert_eq!(resting.field, initial);
        assert_eq!(resting.pressure_correction_pa, vec![1e12, 1e12]);
        assert_eq!(resting.report.projection.iterations, 0);
        assert_eq!(resting.report.open_boundary_pressure_work_j, 0.0);
    }
}

#[test]
fn disconnected_reservoirs_choose_independent_pressure_references() {
    let domain = VoxelDomain::from_fn(5, 1, 1, 1.0, |p| {
        if p[0] > 2.0 && p[0] < 3.0 {
            Voxel::Solid(0)
        } else {
            Voxel::Fluid
        }
    })
    .unwrap();
    let faces = [
        FlowFace::Free,
        FlowFace::Free,
        FlowFace::Wall,
        FlowFace::Wall,
        FlowFace::Wall,
        FlowFace::Wall,
    ];
    let initial = FlowField::quiescent(&domain);
    let project = |pressures| {
        initial
            .project_density_volume_source_with_boundary_pressure(
                &domain,
                faces,
                pressures,
                &[2.0, 5.0, 0.0, 3.0, 7.0],
                &[0.1, -0.1, 0.0, -0.1, 0.1],
                controls(),
                &CancelGate::new(),
            )
            .unwrap()
    };
    let baseline = project([None; 6]);
    let shifted = project([Some(1e12), Some(-2e12), None, None, None, None]);
    assert_eq!(baseline.field, shifted.field);
    assert_eq!(baseline.report, shifted.report);
    assert_eq!(shifted.pressure_correction_pa[0], 1e12);
    assert_eq!(shifted.pressure_correction_pa[4], -2e12);
    assert_eq!(shifted.pressure_correction_pa[2], 0.0);
    close(shifted.field.flux_x(1, 0, 0), 0.1);
    close(shifted.field.flux_x(4, 0, 0), -0.1);
}

#[test]
fn reservoir_pressure_rejects_nonfinite_overconstrained_or_cancelled_inputs() {
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let initial = FlowField::quiescent(&domain);
    let open = [
        FlowFace::Free,
        FlowFace::Free,
        FlowFace::Wall,
        FlowFace::Wall,
        FlowFace::Wall,
        FlowFace::Wall,
    ];
    for blocked in [FlowFace::Wall, FlowFace::Fixed] {
        let mut faces = open;
        faces[0] = blocked;
        let error = initial
            .project_density_volume_source_with_boundary_pressure(
                &domain,
                faces,
                [Some(0.0), None, None, None, None, None],
                &[2.0, 5.0],
                &[0.0, 0.0],
                controls(),
                &CancelGate::new(),
            )
            .unwrap_err();
        assert!(matches!(
            error,
            ChtError::InvalidInput {
                field: "flow.boundary_pressure_pa",
                ..
            }
        ));
    }
    for pressures in [
        [Some(f64::NAN), None, None, None, None, None],
        [Some(f64::INFINITY), None, None, None, None, None],
        [Some(-f64::MAX), Some(f64::MAX), None, None, None, None],
    ] {
        assert!(
            initial
                .project_density_volume_source_with_boundary_pressure(
                    &domain,
                    open,
                    pressures,
                    &[2.0, 5.0],
                    &[0.0, 0.0],
                    controls(),
                    &CancelGate::new(),
                )
                .is_err()
        );
    }
    let gate = CancelGate::new();
    gate.request();
    assert_eq!(
        initial.project_density_volume_source_with_boundary_pressure(
            &domain,
            open,
            [None; 6],
            &[2.0, 5.0],
            &[0.0, 0.0],
            controls(),
            &gate,
        ),
        Err(ChtError::Cancelled)
    );
    assert_eq!(initial, FlowField::quiescent(&domain));
}

fn profile_view(buffers: &[Vec<Option<f64>>; 3]) -> StaggeredBoundaryPressure<'_> {
    StaggeredBoundaryPressure {
        x_pa: &buffers[0],
        y_pa: &buffers[1],
        z_pa: &buffers[2],
    }
}

fn sample_profile(
    field: &FlowField,
    domain: &VoxelDomain,
    faces: [FlowFace; 6],
    pressure: impl Fn(usize, usize, [f64; 3]) -> f64,
) -> [Vec<Option<f64>>; 3] {
    let mut buffers = field.face_arrays().map(|flux| vec![None; flux.len()]);
    for c in 0..domain.cell_count() {
        if !domain.is_fluid(c) {
            continue;
        }
        for f in 0..6 {
            if faces[f] == FlowFace::Free && domain.neighbor(c, f).is_none() {
                let [x, y, z] = domain.coords(c);
                let mut position = domain.center(x, y, z);
                position[f / 2] = domain.face_coord(f / 2, domain.coords(c)[f / 2] + f % 2);
                let (slot, _) = field.slot(domain, c, f);
                buffers[f / 2][slot] = Some(pressure(c, f, position));
            }
        }
    }
    buffers
}

#[test]
fn spatial_affine_reservoir_head_accelerates_uniform_flow_and_supplies_exact_work() {
    let domain =
        VoxelDomain::graded([vec![0.5, 1.0], vec![0.25, 0.75, 1.5], vec![0.75, 1.25]]).unwrap();
    let faces = [FlowFace::Free; 6];
    let initial = FlowField::quiescent(&domain);
    let gradient = [1.25, -2.0, 0.5];
    let head = |position: [f64; 3]| {
        32.0 + (0..3)
            .map(|axis| gradient[axis] * position[axis])
            .sum::<f64>()
    };
    let profile = sample_profile(&initial, &domain, faces, |_, _, p| head(p));
    let rho = 4.0;
    let result = initial
        .project_density_volume_source_with_pressure_profile(
            &domain,
            faces,
            profile_view(&profile),
            &vec![rho; domain.cell_count()],
            &vec![0.0; domain.cell_count()],
            controls(),
            &CancelGate::new(),
        )
        .unwrap();
    let velocity = gradient.map(|g| -controls().duration_s * g / rho);
    for c in 0..domain.cell_count() {
        let [x, y, z] = domain.coords(c);
        close(
            result.pressure_correction_pa[c],
            head(domain.center(x, y, z)),
        );
        for f in 0..6 {
            let sign = if f % 2 == 0 { -1.0 } else { 1.0 };
            close(
                result.field.outward(&domain, c, f),
                sign * domain.face_area(c, f / 2) * velocity[f / 2],
            );
        }
    }
    let volume: f64 = (0..domain.cell_count()).map(|c| domain.volume(c)).sum();
    let kinetic = 0.5 * rho * volume * velocity.iter().map(|u| u * u).sum::<f64>();
    close(result.report.kinetic_energy_after_j, kinetic);
    close(result.report.projection_loss_j, kinetic);
    close(result.report.open_boundary_pressure_work_j, -2.0 * kinetic);
    close(result.report.energy_defect_j, 0.0);
}

#[test]
fn spatial_hydrostatic_head_balances_actual_gravity_predictor_with_sharp_strata() {
    // The side opening spans both densities. A single plane pressure would
    // drive circulation; its resolved head integral rho*g*dh must retain rest.
    for axis in 0..3 {
        let opening_axis = (axis + 1) % 3;
        let mut widths = [vec![0.8], vec![0.8], vec![0.8]];
        widths[axis] = vec![0.5, 1.25];
        widths[opening_axis] = vec![0.75, 1.0];
        let domain = VoxelDomain::graded(widths).unwrap();
        let mut faces = [FlowFace::Wall; 6];
        faces[2 * opening_axis + 1] = FlowFace::Free;
        let initial = FlowField::quiescent(&domain);
        let density: Vec<f64> = (0..domain.cell_count())
            .map(|c| {
                if domain.coords(c)[axis] == 0 {
                    2.0
                } else {
                    5.0
                }
            })
            .collect();
        let mut acceleration = [0.0; 3];
        acceleration[axis] = -3.25;
        let head = |height: f64| {
            60.0 + acceleration[axis] * (2.0 * height.min(0.5) + 5.0 * (height - 0.5).max(0.0))
        };
        let profile = sample_profile(&initial, &domain, faces, |_, _, p| head(p[axis]));
        let predicted = initial
            .predict_density_body_force(
                &domain,
                faces,
                DensityForceInput {
                    density_kg_m3: &density,
                    acceleration_m_s2: acceleration,
                    force_density_n_m3: &vec![[0.0; 3]; domain.cell_count()],
                },
                controls().duration_s,
                DensityForceControls {
                    max_cells: 64,
                    maximum_impulse_defect_kg_m_s: 1e-10,
                    maximum_energy_defect_j: 1e-10,
                },
                &CancelGate::new(),
            )
            .unwrap();
        let result = predicted
            .field
            .project_density_volume_source_with_pressure_profile(
                &domain,
                faces,
                profile_view(&profile),
                &density,
                &vec![0.0; domain.cell_count()],
                controls(),
                &CancelGate::new(),
            )
            .unwrap();
        assert!(predicted.report.kinetic_energy_after_j > 0.0);
        assert!(result.field.max_flux() < 1e-11);
        for c in 0..domain.cell_count() {
            let [x, y, z] = domain.coords(c);
            close(
                result.pressure_correction_pa[c],
                head(domain.center(x, y, z)[axis]),
            );
        }
        close(result.report.kinetic_energy_after_j, 0.0);
        close(
            result.report.projection_loss_j,
            predicted.report.kinetic_energy_after_j,
        );
        close(result.report.open_boundary_pressure_work_j, 0.0);
        close(result.report.energy_defect_j, 0.0);
        assert_eq!(initial, FlowField::quiescent(&domain));
    }
}

#[test]
fn graded_nonuniform_reservoir_profiles_resolve_counterflow_and_face_work() {
    // Two coupled cells stacked in y; each has distinct x- and x+ reservoirs.
    // These are native x-face slots for nx=1, ny=2, nz=1, independently authored.
    let domain = VoxelDomain::graded([vec![2.0], vec![1.0, 3.0], vec![4.0]]).unwrap();
    let faces = [
        FlowFace::Free,
        FlowFace::Free,
        FlowFace::Wall,
        FlowFace::Wall,
        FlowFace::Wall,
        FlowFace::Wall,
    ];
    let initial = FlowField::quiescent(&domain);
    let profile = [
        vec![Some(6.0), Some(1.0), Some(-2.0), Some(3.0)],
        vec![None; 3],
        vec![None; 4],
    ];
    let result = initial
        .project_density_volume_source_with_pressure_profile(
            &domain,
            faces,
            profile_view(&profile),
            &[2.0, 5.0],
            &[0.0; 2],
            controls(),
            &CancelGate::new(),
        )
        .unwrap();
    // Direct two-node hydraulic elimination, using each half-cell inertia.
    let g0 = 0.8;
    let g1 = 0.96;
    let shared = 3.2 / 8.5;
    let a = 2.0 * g0 + shared;
    let d = 2.0 * g1 + shared;
    let determinant = a * d - shared * shared;
    let p0 = (g0 * 7.0 * d + shared * g1) / determinant;
    let p1 = (g1 * a + shared * g0 * 7.0) / determinant;
    let q = [
        g0 * (p0 - 6.0),
        g0 * (p0 - 1.0),
        g1 * (p1 + 2.0),
        g1 * (p1 - 3.0),
    ];
    let transfer = shared * (p0 - p1);
    close(result.pressure_correction_pa[0], p0);
    close(result.pressure_correction_pa[1], p1);
    for c in 0..2 {
        for side in 0..2 {
            close(result.field.outward(&domain, c, side), q[2 * c + side]);
        }
        close(result.field.net_outflow(&domain, c), 0.0);
    }
    close(result.field.flux_y(0, 1, 0), transfer);
    assert!(q[0] < 0.0 && q[2] > 0.0 && q[1] > 0.0 && q[3] < 0.0);
    let kinetic = 0.5
        * (8.0 * ((q[0] / 4.0).powi(2) + (q[1] / 4.0).powi(2))
            + 60.0 * ((q[2] / 12.0).powi(2) + (q[3] / 12.0).powi(2))
            + 68.0 * (transfer / 8.0).powi(2));
    close(result.report.kinetic_energy_after_j, kinetic);
    close(
        result.report.open_boundary_pressure_work_j,
        0.4 * (6.0 * q[0] + q[1] - 2.0 * q[2] + 3.0 * q[3]),
    );
    close(result.report.open_boundary_pressure_work_j, -2.0 * kinetic);
}

#[test]
fn staggered_constant_and_zero_profiles_retain_plane_api_bits() {
    let domain = VoxelDomain::graded([vec![0.5, 1.0], vec![0.25, 0.75], vec![2.0]]).unwrap();
    let faces = [
        FlowFace::Free,
        FlowFace::Free,
        FlowFace::Wall,
        FlowFace::Wall,
        FlowFace::Wall,
        FlowFace::Wall,
    ];
    let initial = FlowField::quiescent(&domain);
    for density in [[2.0; 4], [2.0, 5.0, 2.0, 5.0]] {
        for planes in [[Some(8.0), Some(-2.0), None, None, None, None], [None; 6]] {
            let mut profile =
                sample_profile(&initial, &domain, faces, |_, f, _| planes[f].unwrap_or(0.0));
            if planes == [None; 6] {
                for axis in &mut profile {
                    axis.fill(None);
                }
            }
            let expected = initial
                .project_density_volume_source_with_boundary_pressure(
                    &domain,
                    faces,
                    planes,
                    &density,
                    &[0.25, -0.1, 0.0, 0.05],
                    controls(),
                    &CancelGate::new(),
                )
                .unwrap();
            let actual = initial
                .project_density_volume_source_with_pressure_profile(
                    &domain,
                    faces,
                    profile_view(&profile),
                    &density,
                    &[0.25, -0.1, 0.0, 0.05],
                    controls(),
                    &CancelGate::new(),
                )
                .unwrap();
            assert_eq!(actual, expected);
        }
    }
}

#[test]
fn disconnected_profiles_on_the_same_plane_keep_independent_large_references() {
    let mut domain =
        VoxelDomain::graded([vec![1.0, 2.0], vec![0.5, 1.0, 0.75], vec![1.0]]).unwrap();
    for x in 0..2 {
        domain.set(x, 1, 0, Voxel::Solid(0));
    }
    let faces = [
        FlowFace::Free,
        FlowFace::Free,
        FlowFace::Wall,
        FlowFace::Wall,
        FlowFace::Wall,
        FlowFace::Wall,
    ];
    let initial = FlowField::quiescent(&domain);
    let project = |datum: [f64; 2]| {
        let profile = sample_profile(&initial, &domain, faces, |c, f, _| {
            datum[usize::from(domain.coords(c)[1] == 2)] + if f == 0 { 16.0 } else { 0.0 }
        });
        initial
            .project_density_volume_source_with_pressure_profile(
                &domain,
                faces,
                profile_view(&profile),
                &[2.0, 2.0, 0.0, 0.0, 5.0, 5.0],
                &[0.0; 6],
                controls(),
                &CancelGate::new(),
            )
            .unwrap()
    };
    let baseline = project([0.0; 2]);
    let shifted = project([1e12, -2e12]);
    assert_eq!(shifted.field, baseline.field);
    assert_eq!(shifted.report, baseline.report);
    close(
        shifted.field.flux_x(1, 0, 0),
        0.4 * 0.5 * 16.0 / (2.0 * 3.0),
    );
    close(
        shifted.field.flux_x(1, 2, 0),
        0.4 * 0.75 * 16.0 / (5.0 * 3.0),
    );
    assert_eq!(shifted.pressure_correction_pa[2..4], [0.0; 2]);
}

#[test]
fn staggered_profiles_refuse_invalid_shape_role_data_and_failed_solve_unchanged() {
    let domain = VoxelDomain::new(2, 2, 2, 1.0).unwrap();
    let initial = FlowField::quiescent(&domain);
    let empty = initial.face_arrays().map(|flux| vec![None; flux.len()]);
    let project = |profile: &[Vec<Option<f64>>; 3], faces, controls, gate: &CancelGate| {
        initial.project_density_volume_source_with_pressure_profile(
            &domain,
            faces,
            profile_view(profile),
            &[2.0; 8],
            &[0.0; 8],
            controls,
            gate,
        )
    };
    for axis in 0..3 {
        let mut wrong = empty.clone();
        wrong[axis].push(None);
        assert!(matches!(
            project(&wrong, [FlowFace::Free; 6], controls(), &CancelGate::new()),
            Err(ChtError::InvalidInput {
                field: "flow.boundary_pressure_shape",
                ..
            })
        ));
        let mut interior = empty.clone();
        let (slot, _) = initial.slot(&domain, 0, 2 * axis + 1);
        interior[axis][slot] = Some(0.0);
        assert!(matches!(
            project(
                &interior,
                [FlowFace::Free; 6],
                controls(),
                &CancelGate::new()
            ),
            Err(ChtError::InvalidInput {
                field: "flow.boundary_pressure_role",
                ..
            })
        ));
    }
    for blocked in [FlowFace::Wall, FlowFace::Fixed] {
        let mut profile = empty.clone();
        profile[0][0] = Some(0.0);
        let mut faces = [FlowFace::Free; 6];
        faces[0] = blocked;
        assert!(matches!(
            project(&profile, faces, controls(), &CancelGate::new()),
            Err(ChtError::InvalidInput {
                field: "flow.boundary_pressure_role",
                ..
            })
        ));
    }
    for value in [f64::NAN, f64::INFINITY, -f64::INFINITY] {
        let mut profile = empty.clone();
        profile[0][0] = Some(value);
        assert!(matches!(
            project(
                &profile,
                [FlowFace::Free; 6],
                controls(),
                &CancelGate::new()
            ),
            Err(ChtError::InvalidInput {
                field: "flow.boundary_pressure_pa",
                ..
            })
        ));
    }
    let mut solid_domain = domain.clone();
    solid_domain.set(0, 0, 0, Voxel::Solid(0));
    let solid_field = FlowField::quiescent(&solid_domain);
    let mut profile = empty.clone();
    profile[0][0] = Some(0.0);
    assert!(matches!(
        solid_field.project_density_volume_source_with_pressure_profile(
            &solid_domain,
            [FlowFace::Free; 6],
            profile_view(&profile),
            &[2.0; 8],
            &[0.0; 8],
            controls(),
            &CancelGate::new(),
        ),
        Err(ChtError::InvalidInput {
            field: "flow.boundary_pressure_role",
            ..
        })
    ));
    let profile = sample_profile(&initial, &domain, [FlowFace::Free; 6], |_, _, p| {
        p[0] - 2.0 * p[1]
    });
    let retained = profile.clone();
    assert!(
        project(
            &profile,
            [FlowFace::Free; 6],
            DensityProjectionControls {
                relative_tolerance: 2.0,
                ..controls()
            },
            &CancelGate::new()
        )
        .is_err()
    );
    let gate = CancelGate::new();
    gate.request();
    assert_eq!(
        project(&profile, [FlowFace::Free; 6], controls(), &gate),
        Err(ChtError::Cancelled)
    );
    assert_eq!(profile, retained);
    assert_eq!(initial, FlowField::quiescent(&domain));
}
