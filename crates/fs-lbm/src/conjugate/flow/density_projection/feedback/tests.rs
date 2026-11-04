//! G1/G3/G4: native pressure/force impulses reach actual transported momentum.
use super::*;
use crate::conjugate::domain::Voxel;
use crate::conjugate::flow::{
    AdvectedMaterial, MaterialAdvectionControls, MaterialAdvectionInput, MaterialMomentumControls,
    MaterialMomentumInput, MaterialVolumeSource,
};

fn controls() -> MaterialDensityMomentumControls {
    MaterialDensityMomentumControls {
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
        projection: DensityProjectionControls {
            duration_s: 0.5,
            relative_tolerance: 1e-12,
            maximum_volume_defect_m3_s: 1e-10,
            maximum_energy_defect_j: 1e-7,
        },
        maximum_feedback_momentum_defect_kg_m_s: 1e-8,
        maximum_feedback_angular_impulse_defect_kg_m2_s: 1e-8,
        maximum_feedback_energy_defect_j: 1e-7,
    }
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 4e-10 * expected.abs().max(1.0),
        "{actual:.17e} != {expected:.17e}"
    );
}

fn along_x(minus: FlowFace, plus: FlowFace) -> [FlowFace; 6] {
    [
        minus,
        plus,
        FlowFace::Wall,
        FlowFace::Wall,
        FlowFace::Wall,
        FlowFace::Wall,
    ]
}

fn profile(values: &[Vec<Option<f64>>; 3]) -> StaggeredBoundaryPressure<'_> {
    StaggeredBoundaryPressure {
        x_pa: &values[0],
        y_pa: &values[1],
        z_pa: &values[2],
    }
}

#[test]
fn reservoir_impulse_reaches_unequal_cell_masses_with_gauge_invariant_work() {
    let domain = VoxelDomain::graded([vec![1.0, 2.0], vec![2.0], vec![3.0]]).unwrap();
    let field = FlowField::quiescent(&domain);
    let faces = along_x(FlowFace::Free, FlowFace::Free);
    let gate = CancelGate::new_clock_free();
    for scale in [1.0, 4.0] {
        let mass = [12.0 * scale, 60.0 * scale];
        for difference in [-16.0, 16.0] {
            let mut baseline: Option<MaterialDensityMomentumStep> = None;
            for datum in [0.0, 1e12] {
                let mut pressure = field.face_arrays().map(|q| vec![None; q.len()]);
                pressure[0][0] = Some(datum + difference);
                pressure[0][2] = Some(datum);
                let result = field
                    .force_project_material_momentum_with_pressure_profile(
                        &domain,
                        faces,
                        profile(&pressure),
                        MaterialDensityMomentumInput {
                            mass_kg: &mass,
                            momentum_kg_m_s: &[[0.0; 3]; 2],
                            acceleration_m_s2: [0.0; 3],
                            force_density_n_m3: &[[0.0; 3]; 2],
                            target_outflow_m3_s: &[0.0; 2],
                        },
                        controls(),
                        &gate,
                    )
                    .unwrap();
                // One common column velocity follows from series half-cell
                // inertia. At fixed pressure its speed falls as density grows,
                // while the actual total pressure impulse is independent of rho.
                let dt = controls().projection.duration_s;
                let speed = dt * difference / (12.0 * scale);
                for c in 0..2 {
                    close(result.momentum_kg_m_s[c][0], mass[c] * speed);
                    close(result.momentum_kg_m_s[c][1], 0.0);
                    close(result.momentum_kg_m_s[c][2], 0.0);
                }
                close(
                    result.report.cell_momentum_change_kg_m_s[0],
                    dt * 6.0 * difference,
                );
                close(
                    result.report.kinetic_energy_change_j,
                    36.0 * scale * speed * speed,
                );
                close(result.report.numerical_scatter_loss_j, 0.0);
                close(result.report.composed_energy_defect_j, 0.0);
                close(
                    result
                        .evolution
                        .projection
                        .report
                        .open_boundary_pressure_work_j,
                    -dt * 6.0 * speed * difference,
                );
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
fn stratified_reservoir_head_holds_a_vertical_opening_at_rest() {
    let domain = VoxelDomain::graded([vec![2.0], vec![1.0, 2.0], vec![3.0]]).unwrap();
    let field = FlowField::quiescent(&domain);
    let mut pressure = field.face_arrays().map(|q| vec![None; q.len()]);
    // rho=[5,2], vertical half-cell inertia per area=5/2+2=4.5.
    // g_y=-3 therefore requires a 13.5 Pa drop between the cell centers.
    // Both x-normal openings expose that same resolved hydrostatic head.
    pressure[0] = vec![
        Some(1e12 + 20.0),
        Some(1e12 + 20.0),
        Some(1e12 + 6.5),
        Some(1e12 + 6.5),
    ];
    let result = field
        .force_project_material_momentum_with_pressure_profile(
            &domain,
            along_x(FlowFace::Free, FlowFace::Free),
            profile(&pressure),
            MaterialDensityMomentumInput {
                mass_kg: &[30.0, 24.0],
                momentum_kg_m_s: &[[0.0; 3]; 2],
                acceleration_m_s2: [0.0, -3.0, 0.0],
                force_density_n_m3: &[[0.0; 3]; 2],
                target_outflow_m3_s: &[0.0; 2],
            },
            controls(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    for momentum in result.momentum_kg_m_s {
        for component in momentum {
            close(component, 0.0);
        }
    }
    for flux in result
        .evolution
        .projection
        .field
        .face_arrays()
        .into_iter()
        .flatten()
    {
        close(*flux, 0.0);
    }
    close(
        result
            .evolution
            .report
            .pressure_and_constraint_impulse_kg_m_s[1],
        81.0,
    );
    close(result.evolution.report.applied_endpoint_work_j, 0.0);
    close(result.report.kinetic_energy_change_j, 0.0);
    close(result.report.composed_energy_defect_j, 0.0);
}

#[test]
fn zero_impulse_retains_original_unequal_subface_momenta_without_jensen_debit() {
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let mass = [2.0, 6.0];
    let momentum = [[-2.0, 0.0, -0.0], [6.0, 0.0, -0.0]];
    let result = field
        .force_project_material_momentum(
            &domain,
            along_x(FlowFace::Free, FlowFace::Free),
            MaterialDensityMomentumInput {
                mass_kg: &mass,
                momentum_kg_m_s: &momentum,
                acceleration_m_s2: [0.0; 3],
                force_density_n_m3: &[[0.0; 3]; 2],
                // Reconstructed x velocities are -1, 1/2, 1: these targets
                // require no new pressure impulse. The enclosing source owner
                // would separately supply their material volume evolution.
                target_outflow_m3_s: &[1.5, 0.5],
            },
            controls(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    for (actual, expected) in result.momentum_kg_m_s.iter().zip(momentum) {
        assert_eq!(actual.map(f64::to_bits), expected.map(f64::to_bits));
    }
    assert_eq!(result.reconstruction.numerical_mixing_loss_j, 1.5);
    assert_eq!(result.report.numerical_scatter_loss_j, 0.0);
    assert_eq!(result.report.kinetic_energy_change_j, 0.0);
    assert_eq!(result.report.face_impulse_kg_m_s, [0.0; 3]);
    assert_eq!(result.report.feedback_energy_defect_j, 0.0);
    assert_eq!(result.report.composed_energy_defect_j, 0.0);
}

#[test]
fn graded_density_jump_receives_uniform_acceleration_and_independent_angular_impulse() {
    let domain = VoxelDomain::graded([vec![1.0, 2.0], vec![3.0, 1.0], vec![2.0, 4.0]]).unwrap();
    let field = FlowField::quiescent(&domain);
    let mass: Vec<_> = (0..8)
        .map(|c| domain.volume(c) * (1 + c % 3) as f64)
        .collect();
    let zeros = vec![[0.0; 3]; 8];
    let acceleration = [0.25, -0.5, 0.125];
    let result = field
        .force_project_material_momentum(
            &domain,
            [FlowFace::Free; 6],
            MaterialDensityMomentumInput {
                mass_kg: &mass,
                momentum_kg_m_s: &zeros,
                acceleration_m_s2: acceleration,
                force_density_n_m3: &zeros,
                target_outflow_m3_s: &[0.0; 8],
            },
            controls(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    let mut expected_torque = [0.0; 3];
    let mut returned_torque = [0.0; 3];
    let mut face_torque = [0.0; 3];
    let mut expected_energy = 0.0;
    for c in 0..8 {
        let expected_p = acceleration.map(|g| mass[c] * 0.5 * g);
        for axis in 0..3 {
            close(result.momentum_kg_m_s[c][axis], expected_p[axis]);
        }
        expected_energy +=
            0.5 * mass[c] * acceleration.iter().map(|g| (0.5 * g).powi(2)).sum::<f64>();
        let [x, y, z] = domain.coords(c);
        let r = domain.center(x, y, z);
        // Independent cell and face moments, not the product's cross helper.
        for (target, p) in [
            (&mut expected_torque, expected_p),
            (&mut returned_torque, result.momentum_kg_m_s[c]),
        ] {
            target[0] += r[1] * p[2] - r[2] * p[1];
            target[1] += r[2] * p[0] - r[0] * p[2];
            target[2] += r[0] * p[1] - r[1] * p[0];
        }
        for f in 0..6 {
            let n = domain.neighbor(c, f);
            if n.is_some() && f % 2 == 0 {
                continue;
            }
            let axis = f / 2;
            let sign = if f % 2 == 0 { -1.0 } else { 1.0 };
            let inertia = 0.5 * mass[c] + n.map_or(0.0, |n| 0.5 * mass[n]);
            let u = sign * result.evolution.projection.field.outward(&domain, c, f)
                / domain.face_area(c, axis);
            let mut face_position = r;
            face_position[axis] += sign * 0.5 * domain.widths(c)[axis];
            let mut impulse = [0.0; 3];
            impulse[axis] = inertia * u;
            face_torque[0] += face_position[1] * impulse[2] - face_position[2] * impulse[1];
            face_torque[1] += face_position[2] * impulse[0] - face_position[0] * impulse[2];
            face_torque[2] += face_position[0] * impulse[1] - face_position[1] * impulse[0];
        }
    }
    assert!(expected_torque.iter().any(|v| v.abs() > 1.0));
    for axis in 0..3 {
        close(returned_torque[axis], expected_torque[axis]);
        close(face_torque[axis], expected_torque[axis]);
        close(
            result.report.face_angular_impulse_kg_m2_s[axis],
            expected_torque[axis],
        );
        close(
            result.report.cell_angular_momentum_change_kg_m2_s[axis],
            expected_torque[axis],
        );
    }
    close(result.report.kinetic_energy_change_j, expected_energy);
    close(result.report.numerical_scatter_loss_j, 0.0);
    close(result.report.mass_representation_work_j, 0.0);
    close(result.report.composed_energy_defect_j, 0.0);
}

#[test]
fn wall_and_moving_fixed_reactions_reach_cells_and_keep_reservoir_work_explicit() {
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    for held_velocity in [0.0, 2.0] {
        let faces = if held_velocity == 0.0 {
            [FlowFace::Wall; 6]
        } else {
            along_x(FlowFace::Fixed, FlowFace::Fixed)
        };
        let field = FlowField::from_face_velocity(&domain, faces, |axis, _| {
            if axis == 0 { held_velocity } else { 0.0 }
        })
        .unwrap();
        let result = field
            .force_project_material_momentum(
                &domain,
                faces,
                MaterialDensityMomentumInput {
                    mass_kg: &[8.0],
                    momentum_kg_m_s: &[[8.0, -16.0, 24.0]],
                    acceleration_m_s2: [1.0, -2.0, 3.0],
                    force_density_n_m3: &[[4.0, -8.0, 12.0]],
                    target_outflow_m3_s: &[0.0],
                },
                controls(),
                &CancelGate::new_clock_free(),
            )
            .unwrap();
        assert_eq!(result.momentum_kg_m_s, [[8.0 * held_velocity, 0.0, 0.0]]);
        assert_eq!(
            result.report.face_impulse_kg_m_s,
            [8.0 * (held_velocity - 1.0), 16.0, -24.0]
        );
        assert_eq!(result.evolution.projection.field, field);
        close(result.report.kinetic_energy_before_j, 56.0);
        close(
            result.report.kinetic_energy_after_j,
            4.0 * held_velocity * held_velocity,
        );
        close(
            result.reconstruction.constraint_reservoir_work_j,
            8.0 * (held_velocity - 1.0) * held_velocity,
        );
        close(
            result.evolution.report.applied_endpoint_work_j
                + result.evolution.report.constraint_endpoint_work_j,
            0.0,
        );
        close(result.report.numerical_scatter_loss_j, 0.0);
        close(result.report.composed_energy_defect_j, 0.0);
    }
    let mut solids = VoxelDomain::new(3, 1, 1, 1.0).unwrap();
    solids.set(1, 0, 0, Voxel::Solid(0));
    let result = FlowField::quiescent(&solids)
        .force_project_material_momentum(
            &solids,
            [FlowFace::Wall; 6],
            MaterialDensityMomentumInput {
                mass_kg: &[2.0, 0.0, 3.0],
                momentum_kg_m_s: &[[2.0, 4.0, 6.0], [0.0; 3], [-3.0, -6.0, -9.0]],
                acceleration_m_s2: [0.0; 3],
                force_density_n_m3: &[[0.0; 3]; 3],
                target_outflow_m3_s: &[0.0; 3],
            },
            controls(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    assert_eq!(result.momentum_kg_m_s, [[0.0; 3]; 3]);
    assert_eq!(result.report.cell_momentum_change_kg_m_s, [1.0, 2.0, 3.0]);
    close(result.report.kinetic_energy_change_j, -35.0);
    close(result.report.composed_energy_defect_j, 0.0);
}

#[test]
fn density_baroclinic_pressure_returns_the_analytic_cell_impulses_and_scatter_loss() {
    let domain = VoxelDomain::new(2, 2, 1, 1.0).unwrap();
    let result = FlowField::quiescent(&domain)
        .force_project_material_momentum(
            &domain,
            [FlowFace::Wall; 6],
            MaterialDensityMomentumInput {
                mass_kg: &[1.0, 3.0, 1.0, 3.0],
                momentum_kg_m_s: &[[0.0; 3]; 4],
                acceleration_m_s2: [0.0, -2.0, 0.0],
                force_density_n_m3: &[[0.0; 3]; 4],
                target_outflow_m3_s: &[0.0; 4],
            },
            controls(),
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    // dt*g=1; loop inertia 1+3+2+2=8 gives q=(3-1)/8=1/4.
    let expected = [
        [-0.125, 0.125, 0.0],
        [-0.375, -0.375, 0.0],
        [0.125, 0.125, 0.0],
        [0.375, -0.375, 0.0],
    ];
    for (actual, expected) in result.momentum_kg_m_s.iter().zip(expected) {
        for axis in 0..3 {
            close(actual[axis], expected[axis]);
        }
    }
    close(result.evolution.projection.field.flux_y(0, 1, 0), 0.25);
    close(result.evolution.projection.field.flux_y(1, 1, 0), -0.25);
    close(result.report.cell_momentum_change_kg_m_s[1], -0.5);
    // Face-center orbital moment: 1/4-3/4+1/8-9/8=-3/2.
    close(result.report.face_angular_impulse_kg_m2_s[2], -1.5);
    close(result.report.cell_angular_momentum_change_kg_m2_s[2], -1.5);
    close(result.report.native_impulse_work_j, 0.25);
    close(result.report.kinetic_energy_after_j, 0.125);
    close(result.report.numerical_scatter_loss_j, 0.125);
    close(result.report.composed_energy_defect_j, 0.0);
}

#[test]
fn large_common_velocity_keeps_nonzero_feedback_loss_and_original_momentum_modes() {
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let mut options = controls();
    options.projection.duration_s = 1.0;
    for datum in [0.0, 1e12] {
        let result = FlowField::quiescent(&domain)
            .force_project_material_momentum(
                &domain,
                along_x(FlowFace::Free, FlowFace::Free),
                MaterialDensityMomentumInput {
                    mass_kg: &[2.0, 6.0],
                    momentum_kg_m_s: &[
                        [2.0 * (datum - 1.0), 0.0, 0.0],
                        [6.0 * (datum + 1.0), 0.0, 0.0],
                    ],
                    acceleration_m_s2: [0.0; 3],
                    force_density_n_m3: &[[2.0, 0.0, 0.0], [-6.0, 0.0, 0.0]],
                    target_outflow_m3_s: &[0.0; 2],
                },
                options,
                &CancelGate::new_clock_free(),
            )
            .unwrap();
        // Face velocity kicks 1,-1/2,-1 make every final face equal datum.
        // Cell kicks remain +1/4,-3/4, rather than resetting cells to datum.
        assert_eq!(result.momentum_kg_m_s[0][0], 2.0 * datum - 1.5);
        assert_eq!(result.momentum_kg_m_s[1][0], 6.0 * datum + 1.5);
        assert_eq!(result.report.numerical_scatter_loss_j, 0.75);
        assert_eq!(result.report.kinetic_energy_change_j, -4.0 * datum - 3.25);
        assert_eq!(result.report.native_impulse_work_j, -4.0 * datum - 2.5);
        assert_eq!(result.report.feedback_energy_defect_j, 0.0);
        assert_eq!(result.report.composed_energy_defect_j, 0.0);
    }
}

#[test]
fn native_mass_representation_work_is_budgeted_even_when_the_identity_closes() {
    let domain = VoxelDomain::graded([vec![0.1], vec![0.2], vec![0.3]]).unwrap();
    let field = FlowField::quiescent(&domain);
    let input = MaterialDensityMomentumInput {
        mass_kg: &[1.0],
        momentum_kg_m_s: &[[0.0; 3]],
        acceleration_m_s2: [1e6, 0.0, 0.0],
        force_density_n_m3: &[[0.0; 3]],
        target_outflow_m3_s: &[0.0],
    };
    let mut options = controls();
    options.projection.duration_s = 1.0;
    options.forcing.maximum_impulse_defect_kg_m_s = 1e-6;
    options.forcing.maximum_energy_defect_j = 1.0;
    options.maximum_feedback_momentum_defect_kg_m_s = 1e-6;
    options.maximum_feedback_angular_impulse_defect_kg_m2_s = 1e-6;
    options.maximum_feedback_energy_defect_j = 1.0;
    let result = field
        .force_project_material_momentum(
            &domain,
            [FlowFace::Free; 6],
            input,
            options,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    let rho = 1.0 / (0.1 * 0.2 * 0.3);
    let native_half = (0.2 * 0.3) * (rho * 0.05);
    let impulse = native_half * 1e6;
    let expected = 2.0 * (native_half - 0.5) * (impulse / 1.0) * 1e6;
    close(result.report.mass_representation_work_j, expected);
    assert!(result.report.mass_representation_work_j.abs() > 1e-8);
    assert_eq!(result.report.numerical_scatter_loss_j, 0.0);
    options.maximum_feedback_energy_defect_j = 1e-8;
    assert!(matches!(
        field.force_project_material_momentum(
            &domain,
            [FlowFace::Free; 6],
            input,
            options,
            &CancelGate::new_clock_free()
        ),
        Err(ChtError::InvalidInput {
            field: "flow.feedback_energy_budget",
            ..
        })
    ));
}

#[test]
fn repeated_accepted_transport_pressure_and_feedback_follow_independent_two_cell_oracle() {
    let domain = VoxelDomain::graded([vec![1.0, 2.0], vec![3.0], vec![2.0]]).unwrap();
    let transport_faces = along_x(FlowFace::Fixed, FlowFace::Fixed);
    let feedback_faces = along_x(FlowFace::Fixed, FlowFace::Free);
    let mut field =
        FlowField::from_face_velocity(
            &domain,
            transport_faces,
            |axis, _| if axis == 0 { 1.0 } else { 0.0 },
        )
        .unwrap();
    let mut mass = vec![12.0, 96.0];
    let mut enthalpy = vec![20.0, 80.0];
    let mut momentum = vec![[24.0, 0.0, 0.0], [768.0, 0.0, 0.0]];
    let mut incoming = [None; 6];
    incoming[0] = Some(AdvectedMaterial {
        density_kg_m3: 10.0,
        specific_enthalpy_j_kg: 100.0,
    });
    let mut velocity = [None; 6];
    velocity[0] = Some([10.0, 0.0, 0.0]);
    let gate = CancelGate::new_clock_free();
    for _ in 0..4 {
        // dt*Q=3; backward-Euler donor inventories solve two scalar rows.
        let rho0 = (mass[0] + 30.0) / 9.0;
        let rho1 = (mass[1] + 3.0 * rho0) / 15.0;
        let expected_mass = [6.0 * rho0, 12.0 * rho1];
        let u0 = (momentum[0][0] + 300.0) / (mass[0] + 30.0);
        let u1 = (momentum[1][0] + 3.0 * rho0 * u0) / (mass[1] + 3.0 * rho0);
        let transported = [expected_mass[0] * u0, expected_mass[1] * u1];
        let internal =
            0.5 * (expected_mass[0] + expected_mass[1] - transported[0] - transported[1]);
        let left = 0.5 * (expected_mass[0] - transported[0]);
        let right = 0.5 * (expected_mass[1] - transported[1]);
        let expected_p = [
            transported[0] + left + expected_mass[0] / expected_mass.iter().sum::<f64>() * internal,
            transported[1]
                + right
                + expected_mass[1] / expected_mass.iter().sum::<f64>() * internal,
        ];
        let material = field
            .project_and_advect_material(
                &domain,
                transport_faces,
                MaterialAdvectionInput {
                    mass_kg: &mass,
                    specific_enthalpy_j_kg: &enthalpy,
                    volume_sources: &[MaterialVolumeSource::None; 2],
                    boundary_inflow: incoming,
                },
                MaterialAdvectionControls {
                    projection: controls().projection,
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
        let carried = material
            .advect_momentum(
                MaterialMomentumInput {
                    momentum_kg_m_s: &momentum,
                    impulse_kg_m_s: &[[0.0; 3]; 2],
                    source_velocity_m_s: &[None; 2],
                    boundary_velocity_m_s: velocity,
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
        let saved = material.clone();
        let updated = material
            .projection
            .field
            .force_project_material_momentum(
                &domain,
                feedback_faces,
                MaterialDensityMomentumInput {
                    mass_kg: &material.mass_kg,
                    momentum_kg_m_s: &carried.momentum_kg_m_s,
                    acceleration_m_s2: [0.0; 3],
                    force_density_n_m3: &[[0.0; 3]; 2],
                    target_outflow_m3_s: &[0.0; 2],
                },
                controls(),
                &gate,
            )
            .unwrap();
        for c in 0..2 {
            close(material.mass_kg[c], expected_mass[c]);
            close(carried.momentum_kg_m_s[c][0], transported[c]);
            close(updated.momentum_kg_m_s[c][0], expected_p[c]);
        }
        for i in 0..=2 {
            close(updated.evolution.projection.field.flux_x(i, 0, 0), 6.0);
        }
        close(
            updated.evolution.projection.pressure_correction_pa[1],
            right / 3.0,
        );
        close(
            updated.evolution.projection.pressure_correction_pa[0],
            (right + internal) / 3.0,
        );
        close(updated.report.composed_energy_defect_j, 0.0);
        assert_eq!(material, saved); // No numerical loss is routed into h.
        mass = material.mass_kg;
        enthalpy = material.specific_enthalpy_j_kg;
        momentum = updated.momentum_kg_m_s;
        field = updated.evolution.projection.field;
    }
}

#[test]
fn independent_field_impulse_and_cell_readbacks_refuse_tampering_without_mutation() {
    let domain = VoxelDomain::new(2, 1, 1, 1.0).unwrap();
    let field = FlowField::quiescent(&domain);
    let original = field.clone();
    let faces = [FlowFace::Free; 6];
    let input = MaterialDensityMomentumInput {
        mass_kg: &[2.0, 6.0],
        momentum_kg_m_s: &[[0.0; 3]; 2],
        acceleration_m_s2: [0.0; 3],
        force_density_n_m3: &[[0.0; 3]; 2],
        target_outflow_m3_s: &[0.0; 2],
    };
    let gate = CancelGate::new_clock_free();
    let result = field
        .force_project_material_momentum(&domain, faces, input, controls(), &gate)
        .unwrap();
    let reconstructed = field
        .reconstruct_material_momentum(
            &domain,
            faces,
            input.mass_kg,
            input.momentum_kg_m_s,
            controls().reconstruction,
            &gate,
        )
        .unwrap();
    let exchange = Exchange {
        held: &field,
        domain: &domain,
        faces,
        boundary_pressure: None,
        pressure_readback: None,
        input,
        controls: controls(),
        reconstructed: &reconstructed,
        projection: &result.evolution.projection,
    };
    exchange.admit(&gate).unwrap();
    assert!(matches!(
        exchange.audit(&[[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]], &gate),
        Err(ChtError::InvalidInput {
            field: "flow.feedback_momentum_budget",
            ..
        })
    ));
    let mut unrelated = result.evolution.clone();
    unrelated.projection.pressure_correction_pa[0] = 1.0;
    let wrong_pressure = Exchange {
        projection: &unrelated.projection,
        ..exchange
    };
    let candidate = wrong_pressure.scatter(&gate).unwrap();
    assert!(matches!(
        wrong_pressure.audit(&candidate, &gate),
        Err(ChtError::InvalidInput {
            field: "flow.feedback_momentum_budget",
            ..
        })
    ));
    let mut wrong_start = reconstructed.clone();
    *wrong_start.field.face_mut(0, 1) = 1.0;
    let wrong_start = Exchange {
        reconstructed: &wrong_start,
        projection: &result.evolution.projection,
        ..wrong_pressure
    };
    assert!(matches!(
        wrong_start.scatter(&gate),
        Err(ChtError::InvalidInput {
            field: "flow.feedback_starting_field",
            ..
        })
    ));
    let mut wrong_density = reconstructed.clone();
    wrong_density.density_kg_m3[0] *= 2.0;
    let wrong_density = Exchange {
        reconstructed: &wrong_density,
        ..wrong_start
    };
    assert!(matches!(
        wrong_density.admit(&gate),
        Err(ChtError::InvalidInput {
            field: "flow.feedback_density",
            ..
        })
    ));
    let mut limited = controls();
    limited.maximum_feedback_energy_defect_j = -1.0;
    assert!(
        field
            .force_project_material_momentum(&domain, faces, input, limited, &gate)
            .is_err()
    );
    gate.request();
    assert_eq!(
        field.force_project_material_momentum(&domain, faces, input, controls(), &gate),
        Err(ChtError::Cancelled)
    );
    assert_eq!(field, original);
}
