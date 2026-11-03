//! G1/G3/G4: extensive momentum reaches actual native faces and pressure.
use super::*;
use crate::conjugate::domain::Voxel;
use crate::conjugate::flow::{
    AdvectedMaterial, DensityProjectionControls, MaterialAdvectionControls, MaterialAdvectionInput,
    MaterialMomentumControls, MaterialMomentumInput, MaterialVolumeSource,
};

fn controls() -> DensityMomentumReconstructionControls {
    DensityMomentumReconstructionControls {
        max_cells: 100,
        maximum_mass_defect_kg: 1e-10,
        maximum_momentum_defect_kg_m_s: 1e-9,
        maximum_energy_defect_j: 1e-8,
    }
}

fn projection_controls() -> DensityProjectionControls {
    DensityProjectionControls {
        duration_s: 0.5,
        relative_tolerance: 1e-12,
        maximum_volume_defect_m3_s: 1e-10,
        maximum_energy_defect_j: 1e-8,
    }
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 3e-11 * expected.abs().max(1.0),
        "{actual:.17e} != {expected:.17e}"
    );
}

fn along(axis: usize, minus: FlowFace, plus: FlowFace) -> [FlowFace; 6] {
    let mut faces = [FlowFace::Wall; 6];
    faces[2 * axis] = minus;
    faces[2 * axis + 1] = plus;
    faces
}

fn graded(axis: usize) -> VoxelDomain {
    let mut widths = [vec![2.0], vec![2.0], vec![2.0]];
    widths[axis] = vec![1.0, 2.0];
    widths[(axis + 1) % 3] = vec![3.0];
    VoxelDomain::graded(widths).unwrap()
}

#[test]
fn graded_unequal_density_momentum_has_exact_weighted_velocity_and_jensen_loss() {
    for axis in 0..3 {
        let domain = graded(axis);
        let field = FlowField::quiescent(&domain);
        let faces = along(axis, FlowFace::Free, FlowFace::Free);
        let mut momentum = [[0.0; 3]; 2];
        momentum[0][axis] = 24.0;
        momentum[1][axis] = -60.0;
        let result = field
            .reconstruct_material_momentum(
                &domain,
                faces,
                &[12.0, 60.0],
                &momentum,
                controls(),
                &CancelGate::new_clock_free(),
            )
            .unwrap();
        // Volumes 6,12; densities 2,5; cell velocities 2,-1.
        // Internal p=(24-60)/2=-18, M=(12+60)/2=36, so u=-1/2.
        // The three face momenta 12,-18,-30 sum to original P=-36.
        close(result.field.outward(&domain, 0, 2 * axis), -12.0);
        close(result.field.outward(&domain, 0, 2 * axis + 1), -3.0);
        close(result.field.outward(&domain, 1, 2 * axis + 1), -6.0);
        assert_eq!(result.density_kg_m3, [2.0, 5.0]);
        close(result.report.cell_momentum_kg_m_s[axis], -36.0);
        close(result.report.face_momentum_kg_m_s[axis], -36.0);
        close(result.report.kinetic_energy_before_j, 54.0);
        close(result.report.kinetic_energy_after_j, 31.5);
        close(result.report.kinetic_energy_change_j, -22.5);
        // Half masses 6,30: L=0.5*(6*30/36)*(2-(-1))^2=22.5.
        close(result.report.numerical_mixing_loss_j, 22.5);
        assert_eq!(result.report.constraint_impulse_kg_m_s, [0.0; 3]);
        assert_eq!(result.report.numerical_constraint_loss_j, 0.0);
        assert_eq!(result.report.constraint_reservoir_work_j, 0.0);
        assert_eq!(result.report.mass_representation_energy_j, 0.0);
        assert_eq!(result.report.absolute_mass_defect_kg, 0.0);
        close(result.report.energy_defect_j, 0.0);
        for mass in result.report.native_face_mass_kg {
            close(mass, 72.0);
        }
        assert_eq!(field, FlowField::quiescent(&domain));
    }
}

#[test]
fn uniform_velocity_survives_graded_density_changes_and_native_projection() {
    let domain = VoxelDomain::graded([vec![1.0, 2.0], vec![2.0, 4.0], vec![0.5]]).unwrap();
    let field = FlowField::quiescent(&domain);
    let mass = [1.0, 6.0, 4.0, 20.0];
    let velocity = [2.0, -3.0, 0.25];
    let momentum = mass.map(|m| velocity.map(|u| m * u));
    let faces = [FlowFace::Free; 6];
    let gate = CancelGate::new_clock_free();
    let result = field
        .reconstruct_material_momentum(&domain, faces, &mass, &momentum, controls(), &gate)
        .unwrap();
    for c in 0..4 {
        for f in 0..6 {
            let sign = if f % 2 == 0 { -1.0 } else { 1.0 };
            assert_eq!(
                result.field.outward(&domain, c, f),
                sign * domain.face_area(c, f / 2) * velocity[f / 2]
            );
        }
    }
    assert_eq!(result.report.numerical_mixing_loss_j, 0.0);
    assert_eq!(result.report.kinetic_energy_change_j, 0.0);
    assert_eq!(result.report.absolute_momentum_defect_kg_m_s, 0.0);
    let projected = result
        .field
        .project_density_volume_source(
            &domain,
            faces,
            &result.density_kg_m3,
            &[0.0; 4],
            projection_controls(),
            &gate,
        )
        .unwrap();
    assert_eq!(projected.field, result.field);
}

#[test]
fn held_flow_has_explicit_reservoir_work_and_stationary_wall_projection_loss() {
    let domain = VoxelDomain::graded([vec![2.0], vec![3.0], vec![4.0]]).unwrap();
    let faces = along(0, FlowFace::Fixed, FlowFace::Fixed);
    let field =
        FlowField::from_face_velocity(&domain, faces, |axis, _| if axis == 0 { 2.0 } else { 0.0 })
            .unwrap();
    let result = field
        .reconstruct_material_momentum(
            &domain,
            faces,
            &[120.0],
            &[[120.0, -240.0, 360.0]],
            controls(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert_eq!(result.field, field);
    assert_eq!(
        result.report.constraint_impulse_kg_m_s,
        [120.0, 240.0, -360.0]
    );
    assert_eq!(result.report.face_momentum_kg_m_s, [240.0, 0.0, 0.0]);
    close(result.report.kinetic_energy_before_j, 840.0);
    close(result.report.kinetic_energy_after_j, 240.0);
    close(result.report.kinetic_energy_change_j, -600.0);
    // Only the moving x reservoir supplies work: I_x * u_held=120*2.
    // Constraint losses are 60 (x) + 240 (y) + 540 (z), separately.
    close(result.report.constraint_reservoir_work_j, 240.0);
    close(result.report.numerical_constraint_loss_j, 840.0);
    close(result.report.energy_defect_j, 0.0);
}

#[test]
fn both_solid_interface_orientations_keep_all_mass_momentum_and_support_accounts() {
    let mut domain = VoxelDomain::new(3, 1, 1, 1.0).unwrap();
    domain.set(1, 0, 0, Voxel::Solid(0));
    let field = FlowField::quiescent(&domain);
    let result = field
        .reconstruct_material_momentum(
            &domain,
            [FlowFace::Wall; 6],
            &[2.0, 0.0, 3.0],
            &[[2.0, 4.0, 6.0], [0.0; 3], [-3.0, -6.0, -9.0]],
            controls(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert_eq!(result.field, field);
    assert_eq!(result.density_kg_m3, [2.0, 0.0, 3.0]);
    assert_eq!(result.report.native_face_mass_kg, [5.0; 3]);
    assert_eq!(result.report.cell_momentum_kg_m_s, [-1.0, -2.0, -3.0]);
    assert_eq!(result.report.constraint_impulse_kg_m_s, [1.0, 2.0, 3.0]);
    assert_eq!(result.report.face_momentum_kg_m_s, [0.0; 3]);
    assert_eq!(result.report.constraint_reservoir_work_j, 0.0);
    assert_eq!(result.report.numerical_mixing_loss_j, 0.0);
    close(result.report.kinetic_energy_change_j, -35.0);
    close(result.report.numerical_constraint_loss_j, 35.0);
    close(result.report.energy_defect_j, 0.0);
}

#[test]
fn large_velocity_datum_does_not_erase_the_actual_mixing_loss() {
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let faces = along(0, FlowFace::Free, FlowFace::Free);
    for datum in [0.0, 1e12] {
        let result = field
            .reconstruct_material_momentum(
                &domain,
                faces,
                &[2.0, 6.0],
                &[
                    [2.0 * (datum - 1.0), 0.0, 0.0],
                    [6.0 * (datum + 1.0), 0.0, 0.0],
                ],
                controls(),
                &CancelGate::new_clock_free(),
            )
            .unwrap();
        assert_eq!(result.field.flux_x(1, 0, 0), datum + 0.5);
        assert_eq!(result.report.numerical_mixing_loss_j, 1.5);
        assert_eq!(result.report.kinetic_energy_change_j, -1.5);
        assert_eq!(result.report.energy_defect_j, 0.0);
        assert_eq!(result.report.absolute_momentum_defect_kg_m_s, 0.0);
    }
}

#[test]
fn native_density_mass_roundoff_is_separate_and_cannot_evade_energy_or_mass_budgets() {
    let domain = VoxelDomain::graded([vec![0.1], vec![0.2], vec![0.3]]).unwrap();
    let field = FlowField::quiescent(&domain);
    let faces = [FlowFace::Free; 6];
    let mass = [1.0];
    let momentum = [[1_000_000.0, 0.0, 0.0]];
    let gate = CancelGate::new_clock_free();
    let relaxed = DensityMomentumReconstructionControls {
        maximum_energy_defect_j: 1.0,
        ..controls()
    };
    let result = field
        .reconstruct_material_momentum(&domain, faces, &mass, &momentum, relaxed, &gate)
        .unwrap();
    let rho = 1.0 / (0.1 * 0.2 * 0.3);
    let native_half = (0.2 * 0.3) * (rho * 0.05);
    // Two x half-duals; the other axes carry no momentum. No donor is
    // mixed and no boundary is held, so every nonzero energy account is
    // native-mass/flux representation, never physical or mixing dissipation.
    let expected = (0.5 - native_half) * (500_000.0 / native_half) * 1_000_000.0;
    close(result.report.mass_representation_energy_j, expected);
    assert!(result.report.mass_representation_energy_j.abs() > 1e-9);
    assert!(result.report.absolute_mass_defect_kg > 0.0);
    assert_eq!(result.report.numerical_mixing_loss_j, 0.0);
    assert_eq!(result.report.numerical_constraint_loss_j, 0.0);
    assert_eq!(result.report.constraint_reservoir_work_j, 0.0);
    assert!(matches!(
        field.reconstruct_material_momentum(&domain, faces, &mass, &momentum, controls(), &gate),
        Err(ChtError::InvalidInput {
            field: "flow.reconstruction_energy_budget",
            ..
        })
    ));
    assert!(matches!(
        field.reconstruct_material_momentum(
            &domain,
            faces,
            &mass,
            &momentum,
            DensityMomentumReconstructionControls {
                maximum_mass_defect_kg: 0.0,
                ..relaxed
            },
            &gate
        ),
        Err(ChtError::InvalidInput {
            field: "flow.reconstruction_mass_budget",
            ..
        })
    ));
}

#[test]
fn accepted_mass_and_momentum_transport_feed_real_reconstruction_and_pressure() {
    let domain = graded(0);
    let transport_faces = along(0, FlowFace::Fixed, FlowFace::Fixed);
    let field =
        FlowField::from_face_velocity(
            &domain,
            transport_faces,
            |axis, _| if axis == 0 { 1.0 } else { 0.0 },
        )
        .unwrap();
    let mut incoming = [None; 6];
    incoming[0] = Some(AdvectedMaterial {
        density_kg_m3: 10.0,
        specific_enthalpy_j_kg: 100.0,
    });
    let gate = CancelGate::new_clock_free();
    let material = field
        .project_and_advect_material(
            &domain,
            transport_faces,
            MaterialAdvectionInput {
                mass_kg: &[12.0, 96.0],
                specific_enthalpy_j_kg: &[20.0, 80.0],
                volume_sources: &[MaterialVolumeSource::None; 2],
                boundary_inflow: incoming,
            },
            MaterialAdvectionControls {
                projection: projection_controls(),
                max_cells: 100,
                max_connections: 1000,
                max_iterations: 1000,
                relative_tolerance: 1e-12,
                maximum_mass_defect_kg: 1e-9,
                maximum_enthalpy_defect_j: 1e-7,
            },
            &gate,
        )
        .unwrap();
    let saved = material.clone();
    let mut inlet_velocity = [None; 6];
    inlet_velocity[0] = Some([10.0, 0.0, 0.0]);
    let momentum = material
        .advect_momentum(
            MaterialMomentumInput {
                momentum_kg_m_s: &[[24.0, 0.0, 0.0], [768.0, 0.0, 0.0]],
                impulse_kg_m_s: &[[0.0; 3]; 2],
                source_velocity_m_s: &[None; 2],
                boundary_velocity_m_s: inlet_velocity,
            },
            MaterialMomentumControls {
                max_iterations: 1000,
                relative_tolerance: 1e-12,
                maximum_momentum_defect_kg_m_s: 1e-8,
                maximum_energy_defect_j: 1e-7,
            },
            &gate,
        )
        .unwrap();
    // Exact mass routes: incoming 30, internal 14, outgoing 22 kg.
    // End masses 28,88; (28+14)u0=24+30*10, (88+22)u1=768+14*u0.
    close(material.mass_kg[0], 28.0);
    close(material.mass_kg[1], 88.0);
    close(momentum.momentum_kg_m_s[0][0], 216.0);
    close(momentum.momentum_kg_m_s[1][0], 700.8);
    let faces = along(0, FlowFace::Fixed, FlowFace::Free);
    let reconstructed = material
        .projection
        .field
        .reconstruct_material_momentum(
            &domain,
            faces,
            &material.mass_kg,
            &momentum.momentum_kg_m_s,
            controls(),
            &gate,
        )
        .unwrap();
    close(reconstructed.field.flux_x(1, 0, 0) / 6.0, 458.4 / 58.0);
    close(reconstructed.field.flux_x(2, 0, 0) / 6.0, 350.4 / 44.0);
    assert_eq!(reconstructed.field.flux_x(0, 0, 0), 6.0);
    let projected = reconstructed
        .field
        .project_density_volume_source(
            &domain,
            faces,
            &reconstructed.density_kg_m3,
            &[0.0; 2],
            projection_controls(),
            &gate,
        )
        .unwrap();
    for i in 0..=2 {
        close(projected.field.flux_x(i, 0, 0), 6.0);
    }
    // J_internal=58-458.4; J_outlet=44-350.4. dt*A=3.
    close(projected.pressure_correction_pa[1], (44.0 - 350.4) / 3.0);
    close(projected.pressure_correction_pa[0], -235.6);
    assert_eq!(material, saved);
    assert_eq!(momentum.momentum_kg_m_s.len(), 2);
}

#[test]
fn reconstructed_material_momentum_drives_density_weighted_closed_circulation() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let faces = [FlowFace::Wall; 6];
    let gate = CancelGate::new_clock_free();
    let result = field
        .reconstruct_material_momentum(
            &domain,
            faces,
            &[1.0, 3.0, 1.0, 3.0],
            &[[0.0, 1.0, 0.0], [0.0; 3], [0.0, 1.0, 0.0], [0.0; 3]],
            controls(),
            &gate,
        )
        .unwrap();
    // Internal velocities before pressure: left rises at 1, the other three
    // faces are still. Loop masses 1+3+2+2=8 imply q=1/8, not q=1/4
    // from an unweighted pressure projection.
    close(result.field.flux_y(0, 1, 0), 1.0);
    close(result.report.numerical_constraint_loss_j, 0.5);
    let projected = result
        .field
        .project_density_volume_source(
            &domain,
            faces,
            &result.density_kg_m3,
            &[0.0; 4],
            projection_controls(),
            &gate,
        )
        .unwrap();
    close(projected.field.flux_y(0, 1, 0), 0.125);
    close(projected.field.flux_y(1, 1, 0), -0.125);
    close(projected.field.flux_x(1, 0, 0), -0.125);
    close(projected.field.flux_x(1, 1, 0), 0.125);
    close(projected.report.kinetic_energy_after_j, 0.0625);
    close(projected.report.projection_loss_j, 0.4375);
}

#[test]
fn malformed_inventory_readback_budget_and_cancellation_refuse_without_input_mutation() {
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let original = field.clone();
    let faces = [FlowFace::Free; 6];
    let gate = CancelGate::new_clock_free();
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            field
                .reconstruct_material_momentum(
                    &domain,
                    faces,
                    &[bad],
                    &[[0.0; 3]],
                    controls(),
                    &gate
                )
                .is_err()
        );
    }
    assert!(
        field
            .reconstruct_material_momentum(&domain, faces, &[1.0], &[], controls(), &gate)
            .is_err()
    );
    assert!(
        field
            .reconstruct_material_momentum(
                &domain,
                faces,
                &[1.0],
                &[[f64::NAN, 0.0, 0.0]],
                controls(),
                &gate
            )
            .is_err()
    );
    assert!(matches!(
        field.reconstruct_material_momentum(
            &domain,
            faces,
            &[1.0],
            &[[0.0; 3]],
            DensityMomentumReconstructionControls {
                max_cells: 0,
                ..controls()
            },
            &gate
        ),
        Err(ChtError::InvalidInput {
            field: "flow.reconstruction_cell_budget",
            ..
        })
    ));
    let result = field
        .reconstruct_material_momentum(&domain, faces, &[1.0], &[[0.0; 3]], controls(), &gate)
        .unwrap();
    let mut bad = result.field.clone();
    *bad.face_mut(0, 0) = 1.0;
    *bad.face_mut(0, 1) = -1.0;
    // The tampering cancels globally, but each returned face violates its
    // original half-momentum input. Independent local readback must refuse.
    assert!(matches!(
        audit(
            &field,
            &bad,
            &domain,
            faces,
            &[1.0],
            &[[0.0; 3]],
            &[1.0],
            controls(),
            &gate
        ),
        Err(ChtError::InvalidInput {
            field: "flow.reconstruction_momentum_budget",
            ..
        })
    ));
    assert!(matches!(
        audit(
            &field,
            &bad,
            &domain,
            [FlowFace::Fixed; 6],
            &[1.0],
            &[[0.0; 3]],
            &[1.0],
            controls(),
            &gate
        ),
        Err(ChtError::InvalidInput {
            field: "flow.reconstruction_held_flux",
            ..
        })
    ));
    gate.request();
    assert_eq!(
        field.reconstruct_material_momentum(&domain, faces, &[1.0], &[[0.0; 3]], controls(), &gate),
        Err(ChtError::Cancelled)
    );
    assert_eq!(field, original);

    let mut solid = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    solid.set(1, 0, 0, Voxel::Solid(0));
    assert!(matches!(
        FlowField::quiescent(&solid).reconstruct_material_momentum(
            &solid,
            [FlowFace::Wall; 6],
            &[1.0, 0.0],
            &[[0.0; 3], [1.0, 0.0, 0.0]],
            controls(),
            &CancelGate::new_clock_free()
        ),
        Err(ChtError::InvalidInput {
            field: "flow.reconstruction_solid_inventory",
            ..
        })
    ));
}
