//! G1/G3/G4: independent momentum/kinetic balances on actual accepted mass flux.
use super::super::{
    AdvectedMaterial, DensityProjectionControls, FlowFace, FlowField, MaterialAdvectionControls,
    MaterialAdvectionInput, MaterialVolumeSource, VoxelDomain,
};
use super::*;

fn material_controls(dt: f64) -> MaterialAdvectionControls {
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

fn controls() -> MaterialMomentumControls {
    MaterialMomentumControls {
        max_iterations: 1000,
        relative_tolerance: 1e-12,
        maximum_momentum_defect_kg_m_s: 1e-8,
        maximum_energy_defect_j: 1e-7,
    }
}

fn input<'a>(
    p: &'a [[f64; 3]],
    impulse: &'a [[f64; 3]],
    sources: &'a [Option<[f64; 3]>],
) -> MaterialMomentumInput<'a> {
    MaterialMomentumInput {
        momentum_kg_m_s: p,
        impulse_kg_m_s: impulse,
        source_velocity_m_s: sources,
        boundary_velocity_m_s: [None; 6],
    }
}

fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() < 3e-11 * b.abs().max(1.0),
        "{a:.17e} != {b:.17e}"
    );
}

fn throughflow(axis: usize) -> [FlowFace; 6] {
    let mut faces = [FlowFace::Wall; 6];
    faces[2 * axis] = FlowFace::Fixed;
    faces[2 * axis + 1] = FlowFace::Fixed;
    faces
}

#[test]
fn graded_unequal_density_transport_matches_independent_momentum_and_energy_on_all_axes() {
    for axis in 0..3 {
        let mut widths = [vec![2.0], vec![3.0], vec![2.0]];
        widths[axis] = vec![1.0, 2.0];
        widths[(axis + 1) % 3] = vec![2.0];
        widths[(axis + 2) % 3] = vec![3.0];
        let domain = VoxelDomain::graded(widths).unwrap();
        let faces = throughflow(axis);
        let field =
            FlowField::from_face_velocity(&domain, faces, |a, _| if a == axis { 1.0 } else { 0.0 })
                .unwrap();
        let mut inflow = [None; 6];
        inflow[2 * axis] = Some(AdvectedMaterial {
            density_kg_m3: 10.0,
            specific_enthalpy_j_kg: 100.0,
        });
        let step = field
            .project_and_advect_material(
                &domain,
                faces,
                MaterialAdvectionInput {
                    mass_kg: &[12.0, 96.0],
                    specific_enthalpy_j_kg: &[20.0, 80.0],
                    volume_sources: &[MaterialVolumeSource::None; 2],
                    boundary_inflow: inflow,
                },
                material_controls(0.5),
                &CancelGate::new_clock_free(),
            )
            .unwrap();
        let saved = step.clone();
        let p = [[24.0, -12.0, 36.0], [768.0, 384.0, -192.0]];
        let mut data = input(&p, &[[0.0; 3]; 2], &[None; 2]);
        data.boundary_velocity_m_s[2 * axis] = Some([10.0, 2.0, 1.0]);
        let next = step
            .advect_momentum(data, controls(), &CancelGate::new_clock_free())
            .unwrap();
        // Actual accepted route masses are inlet=30, internal=14, outlet=22.
        // (28+14)u0 = P0+30*u_in; (88+22)u1 = P1+14*u0.
        let expected = [
            [54.0 / 7.0, 8.0 / 7.0, 11.0 / 7.0],
            [438.0 / 55.0, 40.0 / 11.0, -17.0 / 11.0],
        ];
        for c in 0..2 {
            for k in 0..3 {
                close(next.velocity_m_s[c][k], expected[c][k]);
                close(next.momentum_kg_m_s[c][k], step.mass_kg[c] * expected[c][k]);
            }
        }
        for (k, value) in [124.8, -20.0, 64.0].into_iter().enumerate() {
            close(next.report.carried_momentum_kg_m_s[k], value);
            close(
                next.report.momentum_after_kg_m_s[k] - next.report.momentum_before_kg_m_s[k],
                value,
            );
        }
        let final_k = 14.0 * dot(expected[0], expected[0]) + 44.0 * dot(expected[1], expected[1]);
        let carried_k = 1575.0 - 11.0 * dot(expected[1], expected[1]);
        close(next.report.kinetic_energy_before_j, 4116.0);
        close(next.report.kinetic_energy_after_j, final_k);
        close(next.report.carried_kinetic_energy_j, carried_k);
        close(
            next.report.numerical_transport_loss_j,
            carried_k - (final_k - 4116.0),
        );
        assert!(next.report.numerical_transport_loss_j > 0.0);
        assert_eq!(step, saved);
    }
}

#[test]
fn local_material_injection_and_extraction_carry_their_authored_momentum() {
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let step = FlowField::quiescent(&domain)
        .project_and_advect_material(
            &domain,
            [FlowFace::Wall; 6],
            MaterialAdvectionInput {
                mass_kg: &[2.0, 4.0],
                specific_enthalpy_j_kg: &[10.0, 50.0],
                volume_sources: &[
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
                ],
                boundary_inflow: [None; 6],
            },
            material_controls(0.5),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    let next = step
        .advect_momentum(
            input(
                &[[2.0, 4.0, 6.0], [-8.0, 12.0, 4.0]],
                &[[0.0; 3]; 2],
                &[Some([6.0, -1.0, 2.0]), None],
            ),
            controls(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    let expected = [[13.0 / 3.0, 0.0, 7.0 / 3.0], [1.0 / 9.0, 2.0, 13.0 / 9.0]];
    for c in 0..2 {
        for k in 0..3 {
            close(next.velocity_m_s[c][k], expected[c][k]);
        }
    }
    for k in 0..3 {
        close(
            next.report.carried_momentum_kg_m_s[k],
            4.0 * [6.0, -1.0, 2.0][k] - 2.0 * expected[1][k],
        );
    }
    assert!(step.projection.report.projection_loss_j > 0.0);
    // The original pressure impulse is not silently applied to cell momentum.
    assert_eq!(next.report.represented_impulse_kg_m_s, [0.0; 3]);
    close(step.specific_enthalpy_j_kg[0], 70.0);
    close(step.specific_enthalpy_j_kg[1], 170.0 / 3.0);
}

#[test]
fn zero_flow_preserves_authoritative_momentum_and_impulses_have_exact_work_account() {
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let mass = [0.3, 1.7];
    let step = FlowField::quiescent(&domain)
        .project_and_advect_material(
            &domain,
            [FlowFace::Wall; 6],
            MaterialAdvectionInput {
                mass_kg: &mass,
                specific_enthalpy_j_kg: &[17.0, 29.0],
                volume_sources: &[MaterialVolumeSource::None; 2],
                boundary_inflow: [None; 6],
            },
            material_controls(0.5),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    let original = [[0.17, -0.31, 0.53], [2.17, -1.31, 3.53]];
    let mut p = original;
    for _ in 0..32 {
        let next = step
            .advect_momentum(
                input(&p, &[[0.0; 3]; 2], &[None; 2]),
                controls(),
                &CancelGate::new_clock_free(),
            )
            .unwrap();
        assert_eq!(next.momentum_kg_m_s.as_slice(), &original);
        assert_eq!(next.report.kinetic_energy_change_j, 0.0);
        assert_eq!(next.report.numerical_transport_loss_j, 0.0);
        p.copy_from_slice(&next.momentum_kg_m_s);
    }
    let impulse = [[0.06, -0.03, 0.09], [-0.3, 0.2, -0.1]];
    let next = step
        .advect_momentum(
            input(&p, &impulse, &[None; 2]),
            controls(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    let mut work = 0.0;
    for c in 0..2 {
        for k in 0..3 {
            let expected = p[c][k] + impulse[c][k];
            assert_eq!(next.momentum_kg_m_s[c][k], expected);
            work += 0.5 * (expected - p[c][k]) * (expected + p[c][k]) / mass[c];
        }
    }
    close(next.report.impulse_work_j, work);
    close(next.report.kinetic_energy_change_j, work);
    assert_eq!(next.report.numerical_transport_loss_j, 0.0);
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
fn continuing_closed_circulation_conserves_momentum_and_reports_numerical_mixing_loss() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let mut field = circulation(&domain);
    let mut mass = vec![1.0, 2.0, 3.0, 4.0];
    let mut h = vec![10.0, 30.0, 20.0, 40.0];
    let mut p = vec![
        [2.0, 3.0, -1.0],
        [7.0, 11.0, 5.0],
        [-1.0, 2.0, 3.0],
        [9.0, -1.0, 4.0],
    ];
    for _ in 0..32 {
        let step = field
            .project_and_advect_material(
                &domain,
                [FlowFace::Wall; 6],
                MaterialAdvectionInput {
                    mass_kg: &mass,
                    specific_enthalpy_j_kg: &h,
                    volume_sources: &[MaterialVolumeSource::None; 4],
                    boundary_inflow: [None; 6],
                },
                material_controls(0.25),
                &CancelGate::new_clock_free(),
            )
            .unwrap();
        let next = step
            .advect_momentum(
                input(&p, &[[0.0; 3]; 4], &[None; 4]),
                controls(),
                &CancelGate::new_clock_free(),
            )
            .unwrap();
        for (k, expected) in [17.0, 15.0, 11.0].into_iter().enumerate() {
            close(next.report.momentum_after_kg_m_s[k], expected);
        }
        close(
            next.report.kinetic_energy_change_j,
            -next.report.numerical_transport_loss_j,
        );
        assert!(next.report.numerical_transport_loss_j >= 0.0);
        assert_eq!(next.report.carried_momentum_kg_m_s, [0.0; 3]);
        field = step.projection.field;
        mass = step.mass_kg;
        h = step.specific_enthalpy_j_kg;
        p = next.momentum_kg_m_s;
    }
}

#[test]
fn velocity_datum_does_not_set_solver_tolerance_or_change_mixing_loss() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let mass = [1.0, 2.0, 3.0, 4.0];
    let step = circulation(&domain)
        .project_and_advect_material(
            &domain,
            [FlowFace::Wall; 6],
            MaterialAdvectionInput {
                mass_kg: &mass,
                specific_enthalpy_j_kg: &[0.0; 4],
                volume_sources: &[MaterialVolumeSource::None; 4],
                boundary_inflow: [None; 6],
            },
            material_controls(0.25),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    let p = [
        [2.0, 3.0, -1.0],
        [7.0, 11.0, 5.0],
        [-1.0, 2.0, 3.0],
        [9.0, -1.0, 4.0],
    ];
    let base = step
        .advect_momentum(
            input(&p, &[[0.0; 3]; 4], &[None; 4]),
            controls(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    let boost = [1024.0, -2048.0, 4096.0];
    let shifted: Vec<_> = (0..4)
        .map(|c| std::array::from_fn(|k| mass[c].mul_add(boost[k], p[c][k])))
        .collect();
    let moved = step
        .advect_momentum(
            input(&shifted, &[[0.0; 3]; 4], &[None; 4]),
            MaterialMomentumControls {
                maximum_energy_defect_j: 1e-5,
                ..controls()
            },
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    for c in 0..4 {
        for k in 0..3 {
            close(moved.velocity_m_s[c][k] - boost[k], base.velocity_m_s[c][k]);
        }
    }
    close(
        moved.report.numerical_transport_loss_j,
        base.report.numerical_transport_loss_j,
    );
    assert!(
        (moved.report.kinetic_energy_change_j - base.report.kinetic_energy_change_j).abs() < 1e-6
    );
}

#[test]
fn melting_phase_output_reuses_its_original_mass_flux_and_keeps_pressure_work_separate() {
    use super::super::super::{CoexistenceVolumeControls, CoexistenceVolumeInput};
    use fs_blake3::ContentHash;
    use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
    let chart = EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([91; 32]),
        [(0.0, 0.0, 2.0), (4.0, 1.0, 1.0)]
            .into_iter()
            .map(|(h, fraction, rho)| EnthalpyPhaseKnot {
                specific_enthalpy_j_kg: h,
                temperature_k: 400.0,
                liquid_mass_fraction: fraction,
                bulk_density_kg_m3: rho,
            })
            .collect(),
    )
    .unwrap();
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let mut faces = [FlowFace::Wall; 6];
    faces[1] = FlowFace::Free;
    let phase = FlowField::quiescent(&domain)
        .advance_coexistence_volume(
            &domain,
            faces,
            &chart,
            0,
            CoexistenceVolumeInput {
                mass_kg: &[4.0 / 3.0],
                phase_states: &[Some(chart.state_at_specific_enthalpy(2.0).unwrap())],
                net_enthalpy_increment_j: &[4.0 / 3.0],
                boundary_inflow: [None; 6],
            },
            CoexistenceVolumeControls {
                transport: material_controls(0.5),
                maximum_cell_volume_defect_m3: 1e-10,
            },
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    let next = phase
        .transport
        .advect_momentum(
            input(&[[8.0 / 3.0, -4.0, 16.0 / 3.0]], &[[0.0; 3]], &[None]),
            controls(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    for (k, speed) in [2.0, -3.0, 4.0].into_iter().enumerate() {
        close(next.velocity_m_s[0][k], speed);
        close(next.momentum_kg_m_s[0][k], 8.0 / 7.0 * speed);
        close(next.report.carried_momentum_kg_m_s[k], -4.0 / 21.0 * speed);
    }
    close(next.report.carried_kinetic_energy_j, -58.0 / 21.0);
    close(next.report.numerical_transport_loss_j, 0.0);
    assert_eq!(next.report.impulse_work_j, 0.0);
    assert!(phase.report.required_source_pressure_work_j > 0.0);
    close(phase.transport.specific_enthalpy_j_kg[0], 3.0);
}

#[test]
fn retained_mass_roundoff_is_an_energy_defect_and_never_physical_dissipation() {
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let faces = throughflow(0);
    let field =
        FlowField::from_face_velocity(&domain, faces, |axis, _| if axis == 0 { 0.1 } else { 0.0 })
            .unwrap();
    let mut incoming = [None; 6];
    incoming[0] = Some(AdvectedMaterial {
        density_kg_m3: 0.7,
        specific_enthalpy_j_kg: 0.0,
    });
    let step = field
        .project_and_advect_material(
            &domain,
            faces,
            MaterialAdvectionInput {
                mass_kg: &[0.3],
                specific_enthalpy_j_kg: &[0.0],
                volume_sources: &[MaterialVolumeSource::None],
                boundary_inflow: incoming,
            },
            material_controls(0.1),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    let mut data = input(&[[300_000.0, 0.0, 0.0]], &[[0.0; 3]], &[None]);
    data.boundary_velocity_m_s[0] = Some([1_000_000.0, 0.0, 0.0]);
    let next = step
        .advect_momentum(
            data,
            MaterialMomentumControls {
                maximum_energy_defect_j: 1.0,
                ..controls()
            },
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert!(next.report.mass_roundoff_energy_j.abs() > 1e-9);
    assert!(next.report.numerical_transport_loss_j < 1e-15);
    assert!(
        step.advect_momentum(
            data,
            MaterialMomentumControls {
                maximum_energy_defect_j: 0.5 * next.report.mass_roundoff_energy_j.abs(),
                ..controls()
            },
            &CancelGate::new_clock_free()
        )
        .is_err()
    );
}

#[test]
fn missing_inflow_changed_inventory_loose_solves_and_cancellation_refuse_atomically() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let step = circulation(&domain)
        .project_and_advect_material(
            &domain,
            [FlowFace::Wall; 6],
            MaterialAdvectionInput {
                mass_kg: &[1.0, 2.0, 3.0, 4.0],
                specific_enthalpy_j_kg: &[0.0; 4],
                volume_sources: &[MaterialVolumeSource::None; 4],
                boundary_inflow: [None; 6],
            },
            material_controls(0.25),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    let saved = step.clone();
    let p = [
        [2.0, 3.0, -1.0],
        [7.0, 11.0, 5.0],
        [-1.0, 2.0, 3.0],
        [9.0, -1.0, 4.0],
    ];
    let data = input(&p, &[[0.0; 3]; 4], &[None; 4]);
    for budget in [
        MaterialMomentumControls {
            max_iterations: 0,
            ..controls()
        },
        MaterialMomentumControls {
            relative_tolerance: 0.9,
            maximum_momentum_defect_kg_m_s: 1e-12,
            ..controls()
        },
        MaterialMomentumControls {
            maximum_energy_defect_j: f64::NAN,
            ..controls()
        },
    ] {
        assert!(
            step.advect_momentum(data, budget, &CancelGate::new_clock_free())
                .is_err()
        );
    }
    let mut changed = step.clone();
    changed.mass_kg[0] *= 2.0;
    assert!(
        changed
            .advect_momentum(data, controls(), &CancelGate::new_clock_free())
            .is_err()
    );
    let gate = CancelGate::new_clock_free();
    gate.request();
    assert!(matches!(
        step.advect_momentum(data, controls(), &gate),
        Err(ChtError::Cancelled)
    ));
    let mut bad = p;
    bad[0][0] = f64::NAN;
    assert!(
        step.advect_momentum(
            MaterialMomentumInput {
                momentum_kg_m_s: &bad,
                ..data
            },
            controls(),
            &CancelGate::new_clock_free()
        )
        .is_err()
    );
    assert_eq!(step, saved);

    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let faces = throughflow(0);
    let field =
        FlowField::from_face_velocity(&domain, faces, |axis, _| if axis == 0 { 1.0 } else { 0.0 })
            .unwrap();
    let mut incoming = [None; 6];
    incoming[0] = Some(AdvectedMaterial {
        density_kg_m3: 1.0,
        specific_enthalpy_j_kg: 0.0,
    });
    let step = field
        .project_and_advect_material(
            &domain,
            faces,
            MaterialAdvectionInput {
                mass_kg: &[1.0],
                specific_enthalpy_j_kg: &[0.0],
                volume_sources: &[MaterialVolumeSource::None],
                boundary_inflow: incoming,
            },
            material_controls(0.1),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert!(
        step.advect_momentum(
            input(&[[1.0; 3]], &[[0.0; 3]], &[None]),
            controls(),
            &CancelGate::new_clock_free()
        )
        .is_err()
    );
}
