//! Source heat must drive the same accepted material impulse and mass routes.

use super::*;
use crate::conjugate::flow::{
    ConductiveThermalVolumeMomentumControls, DensityForceControls,
    DensityMomentumReconstructionControls, MaterialDensityMomentumControls,
    MaterialMomentumControls, StaggeredBoundaryPressure, ThermalVolumeMomentumControls,
    ThermalVolumeMomentumInput,
};

fn coupled_controls(dt: f64) -> ConductiveThermalVolumeMomentumControls {
    let heat = controls(dt);
    ConductiveThermalVolumeMomentumControls {
        evolution: ThermalVolumeMomentumControls {
            volume: heat.volume,
            mechanics: MaterialDensityMomentumControls {
                reconstruction: DensityMomentumReconstructionControls {
                    max_cells: 100,
                    maximum_mass_defect_kg: 1e-10,
                    maximum_momentum_defect_kg_m_s: 1e-9,
                    maximum_energy_defect_j: 1e-8,
                },
                forcing: DensityForceControls {
                    max_cells: 100,
                    maximum_impulse_defect_kg_m_s: 1e-9,
                    maximum_energy_defect_j: 1e-8,
                },
                projection: heat.volume.thermal.volume.transport.projection,
                maximum_feedback_momentum_defect_kg_m_s: 1e-9,
                maximum_feedback_angular_impulse_defect_kg_m2_s: 1e-9,
                maximum_feedback_energy_defect_j: 1e-8,
            },
            momentum: MaterialMomentumControls {
                max_iterations: 1000,
                relative_tolerance: 1e-12,
                maximum_momentum_defect_kg_m_s: 1e-9,
                maximum_energy_defect_j: 1e-8,
            },
        },
        max_conductive_faces: heat.max_conductive_faces,
    }
}

#[test]
fn original_source_heat_drives_density_volume_and_cell_momentum_together() {
    // T=300+h, k=T-299, v=1+h/8: at h=[1,3] the harmonic conductance
    // is 8/3. dt=3/32 transfers 1/2 J and drives Q=2/3 m3/s.
    // The actual face mass 80/99 gives impulse 160/297. Half-cell
    // scatter gives both cells velocity 1/3 before the SAME mass routes.
    let material = card("T", &[(300.0, 1.0), (304.0, 5.0)]);
    for axis in 0..3 {
        for datum in [0.0, 1024.0] {
            let caloric = curve(
                &material,
                &[
                    (datum, 300.0, 1.0, 1.0),
                    (datum + 4.0, 304.0, 1.0, 2.0 / 3.0),
                ],
            );
            let conductivity = source(&material, &caloric, [datum, datum + 4.0]);
            let mut dims = [1; 3];
            dims[axis] = 2;
            let domain = VoxelDomain::new(dims[0], dims[1], dims[2], 1.0).unwrap();
            let field = FlowField::quiescent(&domain);
            let mass = [8.0 / 9.0, 8.0 / 11.0];
            let states =
                [1.0, 3.0].map(|h| Some(caloric.state_at_specific_enthalpy(datum + h).unwrap()));
            let zero = [[0.0; 3]; 2];
            let data = ThermalVolumeMomentumInput {
                volume: input(&mass, &states, &[0.0; 2]),
                momentum_kg_m_s: &zero,
                acceleration_m_s2: [0.0; 3],
                force_density_n_m3: &zero,
                boundary_velocity_m_s: [None; 6],
            };
            let gate = CancelGate::new_clock_free();
            let result = field
                .advance_isotropic_conductive_thermal_volume_momentum(
                    &domain,
                    [FlowFace::Wall; 6],
                    &conductivity,
                    data,
                    coupled_controls(3.0 / 32.0),
                    &gate,
                )
                .unwrap();
            let thermal_only = field
                .advance_isotropic_conductive_thermal_volume(
                    &domain,
                    [FlowFace::Wall; 6],
                    &conductivity,
                    data.volume,
                    controls(3.0 / 32.0),
                    &gate,
                )
                .unwrap();
            assert_eq!(result.evolution.thermal, thermal_only.thermal);
            assert_eq!(result.report, thermal_only.report);
            assert_eq!(
                result.conductive_enthalpy_increment_j,
                thermal_only.conductive_enthalpy_increment_j
            );
            close(result.report.exchanged_heat_j, 0.5);
            close(result.report.maximum_enthalpy_mixing_weight, 11.0 / 32.0);
            let evolution = &result.evolution;
            let volume = &evolution.thermal.volume;
            let transport = &volume.transport;
            for (c, expected_mass) in [128.0 / 153.0, 1312.0 / 1683.0].into_iter().enumerate() {
                close(transport.mass_kg[c], expected_mass);
                close(
                    transport.mass_kg[c] / volume.phase_states[c].unwrap().bulk_density_kg_m3(),
                    1.0,
                );
                for component in 0..3 {
                    close(
                        evolution.momentum.momentum_kg_m_s[c][component],
                        if component == axis {
                            expected_mass / 3.0
                        } else {
                            0.0
                        },
                    );
                }
                close(evolution.momentum.velocity_m_s[c][axis], 1.0 / 3.0);
            }
            close(
                evolution.force_projection.pressure_impulse_kg_m_s[axis],
                160.0 / 297.0,
            );
            close(
                evolution.feedback.cell_momentum_change_kg_m_s[axis],
                160.0 / 297.0,
            );
            close(evolution.feedback.kinetic_energy_change_j, 80.0 / 891.0);
            close(evolution.feedback.numerical_scatter_loss_j, 80.0 / 891.0);
            close(transport.projection.report.projection_loss_j, 160.0 / 891.0);
            close(volume.report.required_source_pressure_work_j, 320.0 / 891.0);
            close(evolution.feedback.composed_energy_defect_j, 0.0);
            close(evolution.momentum.report.numerical_transport_loss_j, 0.0);
            assert_eq!(
                evolution.momentum.report.represented_impulse_kg_m_s,
                [0.0; 3]
            );
            close(volume.report.material_enthalpy_change_j, 0.0);
            close(transport.report.incoming_mass_kg, 0.0);
            close(transport.report.outgoing_mass_kg, 0.0);
            assert_eq!(field, FlowField::quiescent(&domain));
            assert_eq!(zero, [[0.0; 3]; 2]);
        }
    }
}

#[test]
fn latent_plateau_conducts_no_heat_while_authored_heat_drives_actual_momentum() {
    let material = card("T", &[(300.0, 3.0), (500.0, 3.0)]);
    let caloric = curve(
        &material,
        &[
            (-1.0, 300.0, 0.0, 2.0),
            (0.0, 400.0, 0.0, 2.0),
            (2.0, 400.0, 1.0, 1.0),
            (3.0, 500.0, 1.0, 1.0),
        ],
    );
    let conductivity = source(&material, &caloric, [-1.0, 3.0]);
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let states = [0.0, 2.0].map(|h| Some(caloric.state_at_specific_enthalpy(h).unwrap()));
    let zero = [[0.0; 3]; 2];
    let data = ThermalVolumeMomentumInput {
        volume: input(&[2.0, 1.0], &states, &[0.2, -0.2]),
        momentum_kg_m_s: &zero,
        acceleration_m_s2: [0.0; 3],
        force_density_n_m3: &zero,
        boundary_velocity_m_s: [None; 6],
    };
    let controls = coupled_controls(0.1);
    let gate = CancelGate::new_clock_free();
    let result = field
        .advance_isotropic_conductive_thermal_volume_momentum(
            &domain,
            [FlowFace::Wall; 6],
            &conductivity,
            data,
            controls,
            &gate,
        )
        .unwrap();
    let authored_only = field
        .advance_piecewise_thermal_volume_momentum(
            &domain,
            [FlowFace::Wall; 6],
            &caloric,
            data,
            controls.evolution,
            &gate,
        )
        .unwrap();
    assert_eq!(result.evolution, authored_only);
    assert_eq!(result.conductive_enthalpy_increment_j, [0.0; 2]);
    close(result.report.exchanged_heat_j, 0.0);
    close(result.report.maximum_enthalpy_mixing_weight, 0.0);
    close(result.evolution.momentum.momentum_kg_m_s[0][0], 10.0 / 21.0);
    close(result.evolution.momentum.momentum_kg_m_s[1][0], 23.0 / 84.0);
    close(
        result
            .evolution
            .thermal
            .volume
            .report
            .required_source_pressure_work_j,
        0.375,
    );
}

#[test]
fn sourced_isothermal_reservoir_flow_retains_gauge_and_refuses_missing_inflow_velocity() {
    let material = card("T", &[(300.0, 2.0), (304.0, 2.0)]);
    let caloric = curve(&material, &[(0.0, 300.0, 1.0, 2.0), (4.0, 304.0, 1.0, 2.0)]);
    let conductivity = source(&material, &caloric, [0.0, 4.0]);
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let states = [Some(caloric.state_at_specific_enthalpy(1.0).unwrap())];
    let mass = [2.0];
    let zero = [[0.0; 3]];
    let faces = [
        FlowFace::Free,
        FlowFace::Free,
        FlowFace::Wall,
        FlowFace::Wall,
        FlowFace::Wall,
        FlowFace::Wall,
    ];
    for direction in [-1.0, 1.0] {
        for reference in [0.0, 1e12] {
            let x = if direction > 0.0 {
                [Some(reference + 4.0), Some(reference)]
            } else {
                [Some(reference), Some(reference + 4.0)]
            };
            let transverse = [None; 2];
            let profile = StaggeredBoundaryPressure {
                x_pa: &x,
                y_pa: &transverse,
                z_pa: &transverse,
            };
            let data = ThermalVolumeMomentumInput {
                volume: CoexistenceVolumeInput {
                    boundary_inflow: [states[0]; 6],
                    ..input(&mass, &states, &[0.0])
                },
                momentum_kg_m_s: &zero,
                acceleration_m_s2: [0.0; 3],
                force_density_n_m3: &zero,
                boundary_velocity_m_s: [Some([3.0 * direction, 0.0, 0.0]); 6],
            };
            let gate = CancelGate::new_clock_free();
            let run = |data| {
                field.advance_isotropic_conductive_thermal_volume_momentum_with_pressure_profile(
                    &domain,
                    faces,
                    profile,
                    &conductivity,
                    data,
                    coupled_controls(0.5),
                    &gate,
                )
            };
            let accepted = run(data).unwrap();
            // Total half-cell inertia is M=2: u=dt*dp/M=+/-1.
            // Inlet and outlet each move 1 kg; (2+1)u_new=2*u+1*(3*u).
            close(
                accepted.evolution.momentum.momentum_kg_m_s[0][0],
                direction * 10.0 / 3.0,
            );
            close(
                accepted.evolution.feedback.cell_momentum_change_kg_m_s[0],
                2.0 * direction,
            );
            let transport = &accepted.evolution.thermal.volume.transport;
            close(transport.mass_kg[0], 2.0);
            close(transport.specific_enthalpy_j_kg[0], 1.0);
            close(transport.report.incoming_mass_kg, 1.0);
            close(transport.report.outgoing_mass_kg, 1.0);
            close(
                transport.projection.report.open_boundary_pressure_work_j,
                -2.0,
            );
            close(accepted.report.exchanged_heat_j, 0.0);
            assert!(
                run(ThermalVolumeMomentumInput {
                    boundary_velocity_m_s: [None; 6],
                    ..data
                })
                .is_err()
            );
            assert_eq!(run(data).unwrap(), accepted);
        }
    }
    assert_eq!(mass, [2.0]);
    assert_eq!(zero, [[0.0; 3]]);
    assert_eq!(field, FlowField::quiescent(&domain));
}
