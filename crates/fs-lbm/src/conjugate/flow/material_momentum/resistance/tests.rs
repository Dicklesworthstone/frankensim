//! G1/G3/G4: actual-density drag impulses, physical heat and shared mass routes.

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
fn graded_density_column_receives_implicit_drag_and_gauge_invariant_pressure_once() {
    let domain = VoxelDomain::graded([vec![1.0, 2.0], vec![2.0], vec![3.0]]).unwrap();
    let field = FlowField::quiescent(&domain);
    let gate = CancelGate::new_clock_free();
    for scale in [1.0, 4.0] {
        let mass = [12.0 * scale, 60.0 * scale];
        for difference in [-16.0, 16.0] {
            let mut baseline: Option<ResistiveMaterialDensityMomentumStep> = None;
            for datum in [0.0, 1e12] {
                let mut pressure = [None; 6];
                pressure[0] = Some(datum + difference);
                pressure[1] = Some(datum);
                let result = field
                    .force_project_material_momentum_with_resistance_and_reservoir_pressure(
                        &domain,
                        open_x(),
                        pressure,
                        MaterialDensityMomentumInput {
                            mass_kg: &mass,
                            momentum_kg_m_s: &[[0.0; 3]; 2],
                            acceleration_m_s2: [0.0; 3],
                            force_density_n_m3: &[[0.0; 3]; 2],
                            target_outflow_m3_s: &[0.0; 2],
                        },
                        &[4.0, 10.0],
                        controls().mechanics,
                        &gate,
                    )
                    .unwrap();
                // Every x half-cell contributes once: M=72*scale, B=144.
                // One common endpoint speed obeys (M+dt B)u=dt A delta_p.
                let dt = controls().mechanics.projection.duration_s;
                let speed = dt * 6.0 * difference / (72.0 * scale + dt * 144.0);
                for c in 0..2 {
                    close(result.momentum_kg_m_s[c][0], mass[c] * speed);
                    close(result.density_kg_m3[c], mass[c] / domain.volume(c));
                }
                close(
                    result.report.cell_momentum_change_kg_m_s[0],
                    72.0 * scale * speed,
                );
                close(
                    result.evolution.report.pressure_impulse_kg_m_s[0],
                    dt * 6.0 * difference,
                );
                close(
                    result.evolution.report.resistance_impulse_kg_m_s[0],
                    -dt * 144.0 * speed,
                );
                close(
                    result.evolution.report.resistance_heat_j,
                    dt * 144.0 * speed * speed,
                );
                close(
                    result.evolution.report.backward_euler_loss_j,
                    36.0 * scale * speed * speed,
                );
                close(
                    result.report.kinetic_energy_change_j,
                    36.0 * scale * speed * speed,
                );
                close(result.report.numerical_scatter_loss_j, 0.0);
                close(result.report.composed_energy_defect_j, 0.0);
                if let Some(original) = &baseline {
                    assert_eq!(result.momentum_kg_m_s, original.momentum_kg_m_s);
                    assert_eq!(result.report, original.report);
                    assert_eq!(result.evolution.report, original.evolution.report);
                } else {
                    baseline = Some(result);
                }
            }
        }
    }
}

#[test]
fn held_drag_has_controller_heat_but_only_reconstruction_changes_cell_momentum() {
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let mut faces = [FlowFace::Wall; 6];
    faces[0] = FlowFace::Fixed;
    faces[1] = FlowFace::Fixed;
    let field =
        FlowField::from_face_velocity(&domain, faces, |axis, _| if axis == 0 { 2.0 } else { 0.0 })
            .unwrap();
    let result = field
        .force_project_material_momentum_with_resistance(
            &domain,
            faces,
            MaterialDensityMomentumInput {
                mass_kg: &[2.0],
                momentum_kg_m_s: &[[6.0, 0.0, 0.0]],
                acceleration_m_s2: [1.0, 0.0, 0.0],
                force_density_n_m3: &[[0.0; 3]],
                target_outflow_m3_s: &[0.0],
            },
            &[4.0],
            controls().mechanics,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert_eq!(result.evolution.projection.field, field);
    assert_eq!(result.momentum_kg_m_s, [[4.0, 0.0, 0.0]]);
    close(result.report.cell_momentum_change_kg_m_s[0], -2.0);
    close(result.evolution.report.applied_impulse_kg_m_s[0], 1.0);
    close(result.evolution.report.resistance_impulse_kg_m_s[0], -4.0);
    close(result.evolution.report.constraint_impulse_kg_m_s[0], 3.0);
    close(result.evolution.report.resistance_heat_j, 8.0);
    close(result.evolution.report.body_endpoint_work_j, 2.0);
    close(result.evolution.report.constraint_endpoint_work_j, 6.0);
    close(result.evolution.report.backward_euler_loss_j, 0.0);
    close(result.reconstruction.numerical_constraint_loss_j, 1.0);
    close(result.reconstruction.constraint_reservoir_work_j, -4.0);
    close(result.report.composed_energy_defect_j, 0.0);
}

#[test]
fn opposed_sources_heat_original_cells_and_transport_the_same_kicked_momentum() {
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
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
    let result = field
        .advance_material_momentum_with_resistance(
            &domain,
            [FlowFace::Wall; 6],
            input,
            &[2.0, 6.0],
            controls(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    // The imposed source volume requires Q=1. Native M=3, B=4, so
    // pressure supplies 5 and drag removes 2, leaving the same impulse 3.
    close(result.transport.projection.field.flux_x(1, 0, 0), 1.0);
    close(result.resistance.pressure_impulse_kg_m_s[0], 5.0);
    close(result.resistance.resistance_impulse_kg_m_s[0], -2.0);
    close(result.feedback.cell_momentum_change_kg_m_s[0], 3.0);
    close(result.resistance.cell_resistance_heat_j[0], 0.5);
    close(result.resistance.cell_resistance_heat_j[1], 1.5);
    close(result.represented_resistance_heat_j, 2.0);
    close(result.absolute_heat_insertion_defect_j, 0.0);
    close(result.resistance.backward_euler_loss_j, 1.5);
    close(result.feedback.numerical_scatter_loss_j, 0.75);
    close(result.feedback.composed_energy_defect_j, 0.0);
    assert_eq!(result.transport.report.mass_routes, 3);
    for value in &result.transport.mass_kg {
        close(*value, 4.0);
    }
    // Only physical drag heat precedes the same accepted implicit donors:
    // (4+2)h0=2*10+0.5+4*100; (4+2)h1=4*50+1.5+2*h0.
    let h0 = 420.5 / 6.0;
    let h1 = (201.5 + 2.0 * h0) / 6.0;
    close(result.transport.specific_enthalpy_j_kg[0], h0);
    close(result.transport.specific_enthalpy_j_kg[1], h1);
    close(result.momentum.momentum_kg_m_s[0][0], 50.0 / 3.0);
    close(result.momentum.momentum_kg_m_s[1][0], 62.0 / 9.0);
    assert_eq!(result.momentum.report.represented_impulse_kg_m_s, [0.0; 3]);
    close(
        result.transport.report.enthalpy_change_j,
        result.transport.report.carried_enthalpy_j,
    );

    // With no drag the original material owner and the new route owner agree.
    let original = field
        .advance_material_momentum(
            &domain,
            [FlowFace::Wall; 6],
            input,
            controls(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    let zero_drag = field
        .advance_material_momentum_with_resistance(
            &domain,
            [FlowFace::Wall; 6],
            input,
            &[0.0; 2],
            controls(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert_eq!(zero_drag.transport.mass_kg, original.transport.mass_kg);
    assert_eq!(
        zero_drag.transport.specific_enthalpy_j_kg,
        original.transport.specific_enthalpy_j_kg
    );
    assert_eq!(zero_drag.momentum, original.momentum);
    assert_eq!(zero_drag.resistance.resistance_heat_j, 0.0);
    assert_eq!(zero_drag.represented_resistance_heat_j, 0.0);
}

#[test]
fn resistance_heat_representation_and_late_transport_refusals_are_atomic() {
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let original = field.clone();
    let mass = [2.0];
    let enthalpy = [10.0];
    let momentum = [[2.0, 0.0, 0.0]];
    let zero = [[0.0; 3]];
    let mut inflow = [None; 6];
    inflow[0] = Some(AdvectedMaterial {
        density_kg_m3: 2.0,
        specific_enthalpy_j_kg: 10.0,
    });
    let mut velocity = [None; 6];
    velocity[0] = Some([0.0; 3]);
    let input = MaterialMomentumAdvanceInput {
        material: MaterialAdvectionInput {
            mass_kg: &mass,
            specific_enthalpy_j_kg: &enthalpy,
            volume_sources: &[MaterialVolumeSource::None],
            boundary_inflow: inflow,
        },
        momentum_kg_m_s: &momentum,
        acceleration_m_s2: [0.0; 3],
        force_density_n_m3: &zero,
        source_velocity_m_s: &[None],
        boundary_velocity_m_s: velocity,
    };
    let run = |input, beta, controls, gate: &CancelGate| {
        field.advance_material_momentum_with_resistance(
            &domain,
            open_x(),
            input,
            beta,
            controls,
            gate,
        )
    };
    let expected = run(input, &[4.0], controls(), &CancelGate::new_clock_free()).unwrap();
    close(expected.resistance.resistance_heat_j, 0.5);
    close(expected.resistance.backward_euler_loss_j, 0.25);
    let mut mismatch = controls();
    mismatch.material.projection.duration_s = 0.25;
    let mut late_refusal = controls();
    late_refusal.momentum.max_iterations = 0;
    for control in [mismatch, late_refusal] {
        assert!(run(input, &[4.0], control, &CancelGate::new_clock_free()).is_err());
    }
    let invalid = [[-1.0], [f64::NAN], [f64::INFINITY]];
    for beta in &invalid {
        assert!(run(input, beta, controls(), &CancelGate::new_clock_free()).is_err());
    }
    let unrepresentable_heat = MaterialMomentumAdvanceInput {
        material: MaterialAdvectionInput {
            specific_enthalpy_j_kg: &[1e300],
            ..input.material
        },
        ..input
    };
    assert!(
        run(
            unrepresentable_heat,
            &[4.0],
            controls(),
            &CancelGate::new_clock_free()
        )
        .is_err()
    );
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        run(input, &[4.0], controls(), &cancelled),
        Err(ChtError::Cancelled)
    ));
    assert_eq!(field, original);
    assert_eq!(mass, [2.0]);
    assert_eq!(enthalpy, [10.0]);
    assert_eq!(momentum, [[2.0, 0.0, 0.0]]);
    assert_eq!(
        run(input, &[4.0], controls(), &CancelGate::new_clock_free()).unwrap(),
        expected
    );
}
