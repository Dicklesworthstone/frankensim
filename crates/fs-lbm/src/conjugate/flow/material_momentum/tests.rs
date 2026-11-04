//! G1/G3/G4: composed native mechanics and accepted material momentum routes.

mod pressure;

use super::*;
use crate::conjugate::flow::{
    AdvectedMaterial, DensityForceControls, DensityMomentumReconstructionControls,
    DensityProjectionControls,
};

fn controls() -> MaterialMomentumAdvanceControls {
    let projection = DensityProjectionControls {
        duration_s: 0.5,
        relative_tolerance: 1e-12,
        maximum_volume_defect_m3_s: 1e-10,
        maximum_energy_defect_j: 1e-7,
    };
    MaterialMomentumAdvanceControls {
        mechanics: MaterialDensityMomentumControls {
            reconstruction: DensityMomentumReconstructionControls {
                max_cells: 100,
                maximum_mass_defect_kg: 1e-10,
                maximum_momentum_defect_kg_m_s: 1e-9,
                maximum_energy_defect_j: 1e-7,
            },
            forcing: DensityForceControls {
                max_cells: 100,
                maximum_impulse_defect_kg_m_s: 1e-9,
                maximum_energy_defect_j: 1e-7,
            },
            projection,
            maximum_feedback_momentum_defect_kg_m_s: 1e-8,
            maximum_feedback_angular_impulse_defect_kg_m2_s: 1e-8,
            maximum_feedback_energy_defect_j: 1e-7,
        },
        material: MaterialAdvectionControls {
            projection,
            max_cells: 100,
            max_connections: 1000,
            max_iterations: 1000,
            relative_tolerance: 1e-12,
            maximum_mass_defect_kg: 1e-9,
            maximum_enthalpy_defect_j: 1e-7,
        },
        momentum: MaterialMomentumControls {
            max_iterations: 1000,
            relative_tolerance: 1e-12,
            maximum_momentum_defect_kg_m_s: 1e-8,
            maximum_energy_defect_j: 1e-7,
        },
    }
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 4e-10 * expected.abs().max(1.0),
        "{actual:.17e} != {expected:.17e}"
    );
}

fn open_x() -> [FlowFace; 6] {
    [
        FlowFace::Free,
        FlowFace::Free,
        FlowFace::Wall,
        FlowFace::Wall,
        FlowFace::Wall,
        FlowFace::Wall,
    ]
}

#[test]
fn opposed_volume_sources_keep_gross_exchange_and_one_pressure_impulse() {
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let original_field = field.clone();
    let mass = [2.0, 4.0];
    let enthalpy = [10.0, 50.0];
    let zero = [[0.0; 3]; 2];
    let sources = [
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
    let source_velocity = [Some([6.0, 0.0, 0.0]), None];
    let input = MaterialMomentumAdvanceInput {
        material: MaterialAdvectionInput {
            mass_kg: &mass,
            specific_enthalpy_j_kg: &enthalpy,
            volume_sources: &sources,
            boundary_inflow: [None; 6],
        },
        momentum_kg_m_s: &zero,
        acceleration_m_s2: [0.0; 3],
        force_density_n_m3: &zero,
        source_velocity_m_s: &source_velocity,
        boundary_velocity_m_s: [None; 6],
    };
    let mut missing_velocity = input;
    missing_velocity.source_velocity_m_s = &[None; 2];
    assert!(
        field
            .advance_material_momentum(
                &domain,
                [FlowFace::Wall; 6],
                missing_velocity,
                controls(),
                &CancelGate::new_clock_free(),
            )
            .is_err()
    );
    let result = field
        .advance_material_momentum(
            &domain,
            [FlowFace::Wall; 6],
            input,
            controls(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();

    // Net source volume is zero, but per-cell continuity requires Q01=1.
    // Native dual mass=3 gives pressure impulse 3, shared as P*=[1,2].
    // Accepted route masses are injection=4, interior=2, extraction=2.
    close(result.transport.projection.field.flux_x(1, 0, 0), 1.0);
    close(result.force_projection.pressure_impulse_kg_m_s[0], 3.0);
    close(result.feedback.cell_momentum_change_kg_m_s[0], 3.0);
    close(result.transport.report.incoming_mass_kg, 4.0);
    close(result.transport.report.outgoing_mass_kg, 2.0);
    assert_eq!(result.transport.report.mass_routes, 3);
    for mass in &result.transport.mass_kg {
        close(*mass, 4.0);
    }
    // (4+2)h0=2*10+4*100; (4+2)h1=4*50+2*h0.
    close(result.transport.specific_enthalpy_j_kg[0], 70.0);
    close(result.transport.specific_enthalpy_j_kg[1], 170.0 / 3.0);
    close(result.transport.report.enthalpy_change_j, 860.0 / 3.0);
    close(result.transport.report.carried_enthalpy_j, 860.0 / 3.0);
    // (4+2)u0=1+4*6; (4+2)u1=2+2*u0.
    close(result.momentum.momentum_kg_m_s[0][0], 50.0 / 3.0);
    close(result.momentum.momentum_kg_m_s[1][0], 62.0 / 9.0);
    close(result.momentum.report.carried_momentum_kg_m_s[0], 185.0 / 9.0);
    close(result.momentum.report.momentum_after_kg_m_s[0], 212.0 / 9.0);
    assert_eq!(result.momentum.report.represented_impulse_kg_m_s, [0.0; 3]);
    close(
        result.momentum.report.momentum_after_kg_m_s[0],
        result.feedback.cell_momentum_change_kg_m_s[0]
            + result.momentum.report.carried_momentum_kg_m_s[0],
    );
    assert_eq!(field, original_field);
    assert_eq!(mass, [2.0, 4.0]);
    assert_eq!(enthalpy, [10.0, 50.0]);
    assert_eq!(zero, [[0.0; 3]; 2]);
}

#[test]
fn next_step_uses_endpoint_momentum_instead_of_stale_transport_faces() {
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let faces = open_x();
    let field = FlowField::quiescent(&domain);
    let mass = [2.0];
    let enthalpy = [10.0];
    let zero = [[0.0; 3]];
    let sources = [MaterialVolumeSource::None];
    let source_velocity = [None];
    let mut boundary_inflow = [None; 6];
    boundary_inflow[0] = Some(AdvectedMaterial {
        density_kg_m3: 2.0,
        specific_enthalpy_j_kg: 10.0,
    });
    let mut boundary_velocity = [None; 6];
    boundary_velocity[0] = Some([0.0; 3]);
    let input = MaterialMomentumAdvanceInput {
        material: MaterialAdvectionInput {
            mass_kg: &mass,
            specific_enthalpy_j_kg: &enthalpy,
            volume_sources: &sources,
            boundary_inflow,
        },
        momentum_kg_m_s: &zero,
        acceleration_m_s2: [2.0, 0.0, 0.0],
        force_density_n_m3: &zero,
        source_velocity_m_s: &source_velocity,
        boundary_velocity_m_s: boundary_velocity,
    };
    let first = field
        .advance_material_momentum(
            &domain,
            faces,
            input,
            controls(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    // Physical impulse=2. Transport speed=1; zero-speed inlet brings mass=1.
    // Endpoint (2+1)u=2, so the transport field is not endpoint P/M.
    close(first.forcing.applied_impulse_kg_m_s[0], 2.0);
    close(first.feedback.cell_momentum_change_kg_m_s[0], 2.0);
    close(first.transport.projection.field.flux_x(0, 0, 0), 1.0);
    close(first.transport.mass_kg[0], 2.0);
    close(first.momentum.momentum_kg_m_s[0][0], 4.0 / 3.0);
    close(first.momentum.velocity_m_s[0][0], 2.0 / 3.0);
    assert!(first.momentum.report.numerical_transport_loss_j > 0.0);
    close(first.transport.specific_enthalpy_j_kg[0], 10.0);

    let next_input = MaterialMomentumAdvanceInput {
        material: MaterialAdvectionInput {
            mass_kg: &first.transport.mass_kg,
            specific_enthalpy_j_kg: &first.transport.specific_enthalpy_j_kg,
            ..input.material
        },
        momentum_kg_m_s: &first.momentum.momentum_kg_m_s,
        acceleration_m_s2: [0.0; 3],
        ..input
    };
    let second = first
        .transport
        .projection
        .field
        .advance_material_momentum(
            &domain,
            faces,
            next_input,
            controls(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    // Reconstruct speed=2/3 from returned P/M. New route mass=2/3,
    // giving (2+2/3)u=4/3 and endpoint P=1.
    close(second.transport.projection.field.flux_x(0, 0, 0), 2.0 / 3.0);
    close(second.momentum.momentum_kg_m_s[0][0], 1.0);
    close(second.momentum.velocity_m_s[0][0], 0.5);
    close(second.transport.mass_kg[0], 2.0);
    close(second.transport.specific_enthalpy_j_kg[0], 10.0);
    close(second.feedback.cell_momentum_change_kg_m_s[0], 0.0);

    let unrelated_guess =
        FlowField::from_face_velocity(&domain, faces, |axis, _| if axis == 0 { 17.0 } else { 0.0 })
            .unwrap();
    let same = unrelated_guess
        .advance_material_momentum(
            &domain,
            faces,
            next_input,
            controls(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert_eq!(same, second);
}

#[test]
fn refusal_and_cancellation_leave_a_valid_retry_unchanged() {
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let faces = open_x();
    let field = FlowField::quiescent(&domain);
    let original = field.clone();
    let mass = [2.0];
    let enthalpy = [10.0];
    let zero = [[0.0; 3]];
    let sources = [MaterialVolumeSource::None];
    let source_velocity = [None];
    let mut boundary_inflow = [None; 6];
    boundary_inflow[0] = Some(AdvectedMaterial {
        density_kg_m3: 2.0,
        specific_enthalpy_j_kg: 10.0,
    });
    let mut boundary_velocity = [None; 6];
    boundary_velocity[0] = Some([0.0; 3]);
    let input = MaterialMomentumAdvanceInput {
        material: MaterialAdvectionInput {
            mass_kg: &mass,
            specific_enthalpy_j_kg: &enthalpy,
            volume_sources: &sources,
            boundary_inflow,
        },
        momentum_kg_m_s: &zero,
        acceleration_m_s2: [2.0, 0.0, 0.0],
        force_density_n_m3: &zero,
        source_velocity_m_s: &source_velocity,
        boundary_velocity_m_s: boundary_velocity,
    };
    let run = |data, controls, gate: &CancelGate| {
        field.advance_material_momentum(&domain, faces, data, controls, gate)
    };
    let expected = run(input, controls(), &CancelGate::new_clock_free()).unwrap();
    let mut mismatched_clock = controls();
    mismatched_clock.material.projection.duration_s = 0.25;
    let mut mismatched_gate = controls();
    mismatched_gate.material.projection.maximum_energy_defect_j *= 2.0;
    let mut cell_budget = controls();
    cell_budget.material.max_cells = 0;
    let mut late_solve_budget = controls();
    late_solve_budget.momentum.max_iterations = 0;
    for control in [
        mismatched_clock,
        mismatched_gate,
        cell_budget,
        late_solve_budget,
    ] {
        assert!(run(input, control, &CancelGate::new_clock_free()).is_err());
    }
    let mut missing_boundary = input;
    missing_boundary.boundary_velocity_m_s = [None; 6];
    assert!(run(missing_boundary, controls(), &CancelGate::new_clock_free()).is_err());
    let mut nonfinite_source = input;
    nonfinite_source.source_velocity_m_s = &[Some([f64::NAN, 0.0, 0.0])];
    assert!(run(nonfinite_source, controls(), &CancelGate::new_clock_free()).is_err());
    let gate = CancelGate::new_clock_free();
    gate.request();
    assert!(matches!(run(input, controls(), &gate), Err(ChtError::Cancelled)));
    assert_eq!(field, original);
    assert_eq!(mass, [2.0]);
    assert_eq!(enthalpy, [10.0]);
    assert_eq!(zero, [[0.0; 3]]);
    assert_eq!(
        run(input, controls(), &CancelGate::new_clock_free()).unwrap(),
        expected
    );
}
