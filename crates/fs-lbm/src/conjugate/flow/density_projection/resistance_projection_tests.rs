use super::*;
use crate::conjugate::domain::Voxel;

fn controls() -> DensityProjectionControls {
    DensityProjectionControls {
        duration_s: 0.5,
        relative_tolerance: 1e-12,
        maximum_volume_defect_m3_s: 1e-12,
        maximum_energy_defect_j: 1e-12,
    }
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 2e-12 * expected.abs().max(1.0),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn resistance_mobility_integrates_both_graded_half_volumes_on_every_axis() {
    for axis in 0..3 {
        let mut widths = [vec![2.0], vec![2.0], vec![2.0]];
        widths[axis] = vec![1.0, 3.0];
        widths[(axis + 1) % 3] = vec![4.0];
        let domain = VoxelDomain::graded(widths).unwrap();
        let model = FrozenDensity {
            density_kg_m3: &[2.0, 6.0],
            resistance_pa_s_m2: Some(&[4.0, 10.0]),
            target_outflow_m3_s: &[0.0; 2],
            duration_s: 0.5,
            uniform_density: None,
            boundary_pressure: BoundaryPressure::Planes([0.0; 6]),
        };
        // A=8, half-widths=0.5 and 1.5: actual M=80 kg, B=136 kg/s.
        // dt A^2/(M+dt B) = 32/148. At the outlet M=72 and B=120.
        close(model.face_mass(&domain, 0, 2 * axis + 1, Some(1)), 80.0);
        close(
            face_resistance_kg_s(&domain, &[4.0, 10.0], 0, 2 * axis + 1, Some(1)),
            136.0,
        );
        close(
            model.conductance(&domain, 0, 2 * axis + 1, Some(1)),
            8.0 / 37.0,
        );
        close(
            model.conductance(&domain, 1, 2 * axis + 1, None),
            8.0 / 33.0,
        );
    }
}

#[test]
fn resistive_pressure_keeps_actual_mass_and_signed_correction_work() {
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let faces = [
        FlowFace::Free,
        FlowFace::Free,
        FlowFace::Wall,
        FlowFace::Wall,
        FlowFace::Wall,
        FlowFace::Wall,
    ];
    for (u0, head) in [(0.0, 1.0), (1.0, -1.0)] {
        let initial =
            FlowField::from_face_velocity(
                &domain,
                faces,
                |axis, _| {
                    if axis == 0 { u0 } else { 0.0 }
                },
            )
            .unwrap();
        let (result, _) = initial
            .project_density_volume_source_with_resistance_readback(
                &domain,
                faces,
                BoundaryPressure::Planes([head, 0.0, 0.0, 0.0, 0.0, 0.0]),
                &[2.0],
                &[0.0],
                &[4.0],
                controls(),
                &CancelGate::new(),
            )
            .unwrap();
        // Each boundary dual has physical M=1 and B=2. The mobility mass
        // is 2, so delta_u=head/8; it must not double the kinetic inertia.
        let du = head / 8.0;
        let u1 = u0 + du;
        close(result.field.flux_x(0, 0, 0), u1);
        close(result.field.flux_x(1, 0, 0), u1);
        close(result.pressure_correction_pa[0], head / 2.0);
        close(result.report.kinetic_energy_before_j, u0 * u0);
        close(result.report.kinetic_energy_after_j, u1 * u1);
        close(result.report.kinetic_energy_change_j, (u0 + u1) * du);
        close(result.report.projection_loss_j, du * du);
        close(result.report.resistance_correction_work_j, 2.0 * du * u1);
        close(
            result.report.open_boundary_pressure_work_j,
            -0.5 * head * u1,
        );
        close(result.report.energy_defect_j, 0.0);
        assert_eq!(result.report.resistance_correction_work_j < 0.0, head < 0.0);
    }
}

#[test]
fn zero_resistance_preserves_original_uniform_and_variable_density_bits() {
    let domain = VoxelDomain::graded([vec![0.5, 1.0], vec![0.25, 0.75], vec![2.0]]).unwrap();
    let initial = FlowField::quiescent(&domain);
    let faces = [FlowFace::Free; 6];
    let boundary = BoundaryPressure::Planes([8.0, -2.0, 0.0, 0.0, 0.0, 0.0]);
    let source = [0.25, -0.1, 0.0, 0.05];
    for density in [[2.0; 4], [2.0, 5.0, 2.0, 5.0]] {
        let (expected, old_pressure) = initial
            .project_density_volume_source_with_pressure_readback(
                &domain,
                faces,
                boundary,
                &density,
                &source,
                controls(),
                &CancelGate::new(),
            )
            .unwrap();
        let (actual, new_pressure) = initial
            .project_density_volume_source_with_resistance_readback(
                &domain,
                faces,
                boundary,
                &density,
                &source,
                &[0.0; 4],
                controls(),
                &CancelGate::new(),
            )
            .unwrap();
        assert_eq!(actual, expected);
        assert_eq!(new_pressure.potential.values, old_pressure.potential.values);
        assert_eq!(
            new_pressure.uniform_pressure_scale,
            old_pressure.uniform_pressure_scale
        );
        assert_eq!(
            actual.report.resistance_correction_work_j.to_bits(),
            0.0f64.to_bits()
        );
    }
}

#[test]
fn resistance_projection_refuses_invalid_or_unrepresentable_metrics_unchanged() {
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let initial = FlowField::quiescent(&domain);
    let project = |resistance: &[f64], duration_s| {
        initial.project_density_volume_source_with_resistance_readback(
            &domain,
            [FlowFace::Free; 6],
            BoundaryPressure::Planes([0.0; 6]),
            &[2.0],
            &[0.0],
            resistance,
            DensityProjectionControls {
                duration_s,
                ..controls()
            },
            &CancelGate::new(),
        )
    };
    for resistance in [&[][..], &[0.0, 0.0], &[-1.0], &[f64::NAN], &[f64::INFINITY]] {
        assert!(project(resistance, 0.5).is_err());
    }
    assert!(matches!(
        project(&[f64::MAX], 4.0),
        Err(ChtError::InvalidInput {
            field: "flow.pressure_mobility_mass_kg",
            ..
        })
    ));
    let mut solid = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    solid.set(1, 0, 0, Voxel::Solid(0));
    assert!(matches!(
        admit_resistance_field(&solid, &[0.0, 1.0], &CancelGate::new()),
        Err(ChtError::InvalidInput {
            field: "flow.resistance_solid",
            ..
        })
    ));
    assert_eq!(initial, FlowField::quiescent(&domain));
}
