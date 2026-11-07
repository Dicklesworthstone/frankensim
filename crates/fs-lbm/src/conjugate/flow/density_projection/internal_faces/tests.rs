use super::*;
use crate::conjugate::Voxel;

fn controls() -> DensityProjectionControls {
    DensityProjectionControls {
        duration_s: 0.4,
        relative_tolerance: 1e-12,
        maximum_volume_defect_m3_s: 1e-11,
        maximum_energy_defect_j: 1e-10,
    }
}

fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 2e-11 * b.abs().max(1.0), "{a:e} != {b:e}");
}

#[test]
fn g1_internal_gas_face_reuses_physical_density_flux_pressure_and_work() {
    // A=6; fluid half-widths .5 and 1. The internal dual mass is
    // 6*(2*.5+5*1)=36; outlet mass is 6*5*1=30, independent of gas.
    // dt=.4 gives G_internal=.4, G_out=.48. Sources [.8,.2] give
    // Q_internal=.8, Q_out=1, p1=25/12, p0=49/12 above gas pressure.
    for axis in 0..3 {
        let mut widths = [vec![2.0], vec![3.0], vec![2.0]];
        widths[axis] = vec![1.0, 2.0, 9.0];
        widths[(axis + 1) % 3] = vec![2.0];
        widths[(axis + 2) % 3] = vec![3.0];
        let mut domain = VoxelDomain::graded(widths).unwrap();
        let gas = domain.coords(2);
        domain.set(gas[0], gas[1], gas[2], Voxel::Solid(0));
        let initial = FlowField::quiescent(&domain);
        let mut reference = None;
        for offset in [0.0, 1.0e12] {
            let opening = InternalPressureFace {
                cell: 1,
                face: Face3::ALL[2 * axis + 1],
                pressure_pa: offset,
            };
            let result = initial
                .project_density_volume_source_with_internal_pressure_faces(
                    &domain,
                    [FlowFace::Wall; 6],
                    &[opening],
                    &[2.0, 5.0, f64::NAN],
                    &[0.8, 0.2, 0.0],
                    controls(),
                    &CancelGate::new_clock_free(),
                )
                .unwrap();
            close(result.field.outward(&domain, 0, 2 * axis + 1), 0.8);
            close(result.field.outward(&domain, 1, 2 * axis + 1), 1.0);
            close(result.field.net_outflow(&domain, 0), 0.8);
            close(result.field.net_outflow(&domain, 1), 0.2);
            assert_eq!(result.pressure_correction_pa[2], 0.0);
            assert_eq!(result.report.projection.pinned_components, 0);
            close(result.report.kinetic_energy_after_j, 221.0 / 300.0);
            close(result.report.projection_loss_j, 221.0 / 300.0);
            close(result.report.pressure_source_work_j, 221.0 / 150.0);
            close(result.report.open_boundary_pressure_work_j, 0.0);
            assert!(result.report.energy_defect_j.abs() < 1e-12);
            if offset == 0.0 {
                close(result.pressure_correction_pa[0], 49.0 / 12.0);
                close(result.pressure_correction_pa[1], 25.0 / 12.0);
                reference = Some(result.field.clone());
            } else {
                assert_eq!(result.field, reference.clone().unwrap());
            }
            // An accepted opening flux is legal on the same declared topology.
            let again = result
                .field
                .project_density_volume_source_with_internal_pressure_faces(
                    &domain,
                    [FlowFace::Wall; 6],
                    &[opening],
                    &[2.0, 5.0, f64::NAN],
                    &[0.8, 0.2, 0.0],
                    controls(),
                    &CancelGate::new_clock_free(),
                )
                .unwrap();
            close(again.field.outward(&domain, 1, 2 * axis + 1), 1.0);
        }
        assert_eq!(initial, FlowField::quiescent(&domain));
    }
}

#[test]
fn g1_opposite_gas_pressures_and_predictor_use_only_the_liquid_side() {
    let mut domain = VoxelDomain::new(3, 1, 1, 1.0).unwrap();
    domain.set(0, 0, 0, Voxel::Solid(0));
    domain.set(2, 0, 0, Voxel::Solid(0));
    for offset in [0.0, 1.0e12] {
        let openings = [
            InternalPressureFace {
                cell: 1,
                face: Face3::XMin,
                pressure_pa: offset,
            },
            InternalPressureFace {
                cell: 1,
                face: Face3::XMax,
                pressure_pa: offset + 2.0,
            },
        ];
        let gate = CancelGate::new_clock_free();
        let initial = FlowField::from_cell_velocity_predictor_with_internal_pressure_faces(
            &domain,
            &[[f64::NAN; 3], [0.0; 3], [f64::NAN; 3]],
            [FlowFace::Wall; 6],
            &openings,
            &gate,
        )
        .unwrap();
        let result = initial
            .project_density_volume_source_with_internal_pressure_faces(
                &domain,
                [FlowFace::Wall; 6],
                &openings,
                &[f64::NAN, 4.0, f64::NAN],
                &[0.0; 3],
                DensityProjectionControls {
                    duration_s: 0.5,
                    ..controls()
                },
                &gate,
            )
            .unwrap();
        close(result.field.outward(&domain, 1, 0), 0.25);
        close(result.field.outward(&domain, 1, 1), -0.25);
        close(result.pressure_correction_pa[1] - offset, 1.0);
        close(result.report.kinetic_energy_after_j, 0.125);
        close(result.report.open_boundary_pressure_work_j, -0.25);
        close(result.report.energy_defect_j, 0.0);
        let moving = FlowField::from_cell_velocity_predictor_with_internal_pressure_faces(
            &domain,
            &[[f64::NAN; 3], [0.3, 0.0, 0.0], [f64::NAN; 3]],
            [FlowFace::Wall; 6],
            &openings,
            &gate,
        )
        .unwrap();
        close(moving.outward(&domain, 1, 0), -0.3);
        close(moving.outward(&domain, 1, 1), 0.3);
    }
}

#[test]
fn g0_each_disconnected_region_keeps_its_own_volume_compatibility() {
    let mut domain = VoxelDomain::new(3, 1, 1, 1.0).unwrap();
    domain.set(1, 0, 0, Voxel::Solid(0));
    let initial = FlowField::quiescent(&domain);
    let openings = [
        InternalPressureFace {
            cell: 0,
            face: Face3::XMax,
            pressure_pa: 10.0,
        },
        InternalPressureFace {
            cell: 2,
            face: Face3::XMin,
            pressure_pa: 20.0,
        },
    ];
    let solve = |open: &[InternalPressureFace]| {
        initial.project_density_volume_source_with_internal_pressure_faces(
            &domain,
            [FlowFace::Wall; 6],
            open,
            &[2.0, f64::NAN, 5.0],
            &[1.0, 0.0, -1.0],
            controls(),
            &CancelGate::new_clock_free(),
        )
    };
    // Global source cancellation does not authorize flow through the gas.
    assert!(solve(&[]).is_err());
    assert!(solve(&openings[..1]).is_err());
    let result = solve(&openings).unwrap();
    close(result.field.outward(&domain, 0, 1), 1.0);
    close(result.field.outward(&domain, 2, 0), -1.0);
    close(result.report.energy_defect_j, 0.0);
    assert_eq!(initial, FlowField::quiescent(&domain));
}

#[test]
fn g0_pressure_roles_duplicates_and_cancellation_refuse_without_mutation() {
    let mut domain = VoxelDomain::new(3, 1, 1, 1.0).unwrap();
    domain.set(2, 0, 0, Voxel::Solid(0));
    let initial = FlowField::quiescent(&domain);
    let opening = InternalPressureFace {
        cell: 1,
        face: Face3::XMax,
        pressure_pa: 0.0,
    };
    let gate = CancelGate::new_clock_free();
    for invalid in [
        vec![opening, opening],
        vec![InternalPressureFace { cell: 0, ..opening }],
        vec![InternalPressureFace { cell: 3, ..opening }],
        vec![InternalPressureFace {
            pressure_pa: f64::NAN,
            ..opening
        }],
        vec![InternalPressureFace {
            face: Face3::YMax,
            ..opening
        }],
    ] {
        assert!(
            initial
                .project_density_volume_source_with_internal_pressure_faces(
                    &domain,
                    [FlowFace::Wall; 6],
                    &invalid,
                    &[1.0, 1.0, f64::NAN],
                    &[0.0; 3],
                    controls(),
                    &gate,
                )
                .is_err()
        );
    }
    gate.request();
    assert!(matches!(
        initial.project_density_volume_source_with_internal_pressure_faces(
            &domain,
            [FlowFace::Wall; 6],
            &[opening],
            &[1.0, 1.0, f64::NAN],
            &[0.0; 3],
            controls(),
            &gate,
        ),
        Err(ChtError::Cancelled)
    ));
    assert_eq!(initial, FlowField::quiescent(&domain));
}

#[test]
fn g0_accepted_pressure_opening_cannot_enter_ordinary_solid_heat_transport() {
    use crate::conjugate::{
        EnergyConfig, FluidProperties, SolidMaterial, ThermalFace, ThermalSetup, solve_energy,
    };
    let mut domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    domain.set(1, 0, 0, Voxel::Solid(0));
    let gate = CancelGate::new_clock_free();
    let fluid = FluidProperties {
        density_kg_m3: 1.0,
        specific_heat_j_kg_k: 1.0,
        conductivity_w_m_k: 1.0,
        kinematic_viscosity_m2_s: 1.0,
    };
    let solids = [SolidMaterial::new("inactive neighbour", 1.0)];
    let mut setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
    setup.fixed_temperature.push((0, 300.0));
    for source in [-0.1, 0.1] {
        let projected = FlowField::quiescent(&domain)
            .project_density_volume_source_with_internal_pressure_faces(
                &domain,
                [FlowFace::Wall; 6],
                &[InternalPressureFace {
                    cell: 0,
                    face: Face3::XMax,
                    pressure_pa: 0.0,
                }],
                &[1.0, f64::NAN],
                &[source, 0.0],
                controls(),
                &gate,
            )
            .unwrap();
        close(projected.field.outward(&domain, 0, 1), source);
        assert!(matches!(
            solve_energy(
                &domain,
                &fluid,
                &solids,
                &projected.field,
                &setup,
                &EnergyConfig::default(),
                &gate
            ),
            Err(ChtError::InvalidInput {
                field: "thermal.blocked_face_flux",
                ..
            })
        ));
    }
}
