//! Complete original solids driven by continuous affine nodal force records.
use super::*;

fn tet10(cx: &Cx<'_>) -> (Vec<[f64; 3]>, ViscoElements) {
    let (nodes, cells) = promote_p1_to_p2(&CORNERS, &[[0, 1, 2, 3]], None, cx).unwrap();
    let law = material(0.0001, 0.4);
    let source = ViscoElements::new(
        SolidMeshRef::Tet10 {
            nodes_m: &nodes,
            elements: &cells,
            quadrature: TetQuadrature {
                stiffness_points_per_dir: 3,
                mass_points_per_dir: 4,
            },
        },
        TetViscoField::Uniform(&law),
        &[],
        TetAssemblyBudget::standard(),
        ViscoStateBudget {
            max_internal_variables: 1000,
        },
        cx,
    )
    .unwrap();
    (nodes, source)
}

fn force_for_acceleration(source: &ViscoElements, acceleration: [f64; 3]) -> Vec<[f64; 3]> {
    let values: Vec<_> = (0..30).map(|i| acceleration[i % 3]).collect();
    let mut force = vec![0.0; 30];
    source.equilibrium().mass.spmv(&values, &mut force);
    force.chunks_exact(3).map(|f| [f[0], f[1], f[2]]).collect()
}

#[test]
fn g1_full_tet10_affine_body_force_matches_rigid_motion_work_and_moments() {
    with_cx(|cx, gate| {
        let (_, source) = tet10(cx);
        let dt = 1.0 / 4800.0;
        let body = PreparedTensorSolidBody::new_affine(
            source
                .project_basis(&identity(30), projection_budget(), cx)
                .unwrap(),
            dt,
            1.0,
            controls(),
            gate,
        )
        .unwrap();
        assert_eq!(body.material_body().state_dimension(), 90);
        assert_eq!(body.material_body().propagation_state_dimension(), 120);
        assert_eq!(body.projection().material_points().len(), 27);
        let initial_velocity = [0.03, -0.02, 0.01];
        let v0: Vec<_> = (0..30).map(|i| initial_velocity[i % 3]).collect();
        let mut history = body
            .history(&[0.0; 30], &v0, InitialMemory::Relaxed, gate)
            .unwrap();
        let a0 = [2.0, -1.0, 0.3];
        let a1 = [-0.5, 1.8, -0.2];
        let f0 = force_for_acceleration(&source, a0);
        let f1 = force_for_acceleration(&source, a1);
        let report = history.advance_affine(&f0, &f1, gate).unwrap();
        let u = history.nodal_displacement_m(gate).unwrap();
        let v = history.nodal_velocity_mps(gate).unwrap();
        let expected_v: [f64; 3] =
            core::array::from_fn(|j| initial_velocity[j] + 0.5 * dt * (a0[j] + a1[j]));
        let expected_u: [f64; 3] = core::array::from_fn(|j| {
            dt * initial_velocity[j] + dt * dt * (a0[j] / 3.0 + a1[j] / 6.0)
        });
        let first: [f64; 3] = core::array::from_fn(|j| {
            0.5 * initial_velocity[j] + dt * (a0[j] / 3.0 + (a1[j] - a0[j]) / 8.0)
        });
        for node in 0..10 {
            for axis in 0..3 {
                assert!((u[node][axis] - expected_u[axis]).abs() < 3e-13);
                assert!((v[node][axis] - expected_v[axis]).abs() < 2e-10);
                assert!(
                    (report.motion.mean_nodal_velocity_mps[node][axis] - expected_u[axis] / dt)
                        .abs()
                        < 2e-10
                );
                assert!(
                    (report.first_nodal_velocity_moment_mps[node][axis] - first[axis]).abs()
                        < 2e-10
                );
            }
        }
        let mass = source.equilibrium().total_mass_kg;
        let expected_work = 0.5
            * mass
            * (0..3)
                .map(|j| expected_v[j] * expected_v[j] - initial_velocity[j] * initial_velocity[j])
                .sum::<f64>();
        assert!((report.motion.material.energy.supplied_work_j - expected_work).abs() < 3e-13);
        assert!(report.motion.element_heat_j()[0].abs() < 2e-15);
        assert!(report.motion.material.energy.equilibrium_after_j.abs() < 2e-15);
        assert!(report.motion.material.energy.memory_after_j.abs() < 2e-15);
        assert!(report.motion.nodal_work_error_j.abs() < 3e-13);
        let flat_velocity: Vec<_> = v.iter().flatten().copied().collect();
        let mut momentum = vec![0.0; 30];
        source
            .equilibrium()
            .mass
            .spmv(&flat_velocity, &mut momentum);
        for axis in 0..3 {
            let actual: f64 = (0..10).map(|node| momentum[3 * node + axis]).sum();
            let impulse: f64 = f0
                .iter()
                .zip(&f1)
                .map(|(a, b)| 0.5 * dt * (a[axis] + b[axis]))
                .sum();
            assert!((actual - mass * initial_velocity[axis] - impulse).abs() < 3e-12);
        }
        println!(
            "affine TET10 rigid motion: physical_states=90, propagation_states=120, exact_work_j={expected_work:e}, nodal_work_defect_j={:e}",
            report.motion.nodal_work_error_j
        );
    });
}

fn triangular_force(time_in_intervals: f64) -> [[f64; 3]; 10] {
    let magnitude = if time_in_intervals <= 8.0 {
        0.2 * time_in_intervals / 8.0
    } else if time_in_intervals < 16.0 {
        0.2 * (16.0 - time_in_intervals) / 8.0
    } else {
        0.0
    };
    let mut force = [[0.0; 3]; 10];
    force[0][0] = -magnitude;
    force[1][0] = magnitude;
    force
}

#[test]
fn g3_full_tet10_affine_pulse_keeps_motion_work_and_heat_under_interval_splitting() {
    with_cx(|cx, gate| {
        let (nodes, source) = tet10(cx);
        let initial: Vec<_> = nodes.iter().flat_map(|p| [1e-5 * p[0], 0.0, 0.0]).collect();
        let initial_energy = 0.5 * (0.2 * 0.2 / 6.0) * 120e6 * 1.4 * 1e-10;
        let mut results = Vec::new();
        for subdivisions in [1, 2, 4] {
            let dt = 1.0 / (48000.0 * subdivisions as f64);
            let body = PreparedTensorSolidBody::new_affine(
                source
                    .project_basis(&identity(30), projection_budget(), cx)
                    .unwrap(),
                dt,
                0.2,
                controls(),
                gate,
            )
            .unwrap();
            let mut history = body
                .history(&initial, &[0.0; 30], InitialMemory::Instantaneous, gate)
                .unwrap();
            let mut independent_work = 0.0;
            let mut heat = 0.0;
            let mut maximum_balance = 0.0_f64;
            for step in 0..32 * subdivisions {
                let f0 = triangular_force(step as f64 / subdivisions as f64);
                let f1 = triangular_force((step + 1) as f64 / subdivisions as f64);
                let report = history.advance_affine(&f0, &f1, gate).unwrap();
                let mut nodal_work = 0.0;
                for node in 0..10 {
                    for axis in 0..3 {
                        nodal_work += dt
                            * (f0[node][axis] * report.motion.mean_nodal_velocity_mps[node][axis]
                                + (f1[node][axis] - f0[node][axis])
                                    * report.first_nodal_velocity_moment_mps[node][axis]);
                    }
                }
                independent_work += nodal_work;
                assert!(report.motion.element_heat_j()[0] >= 0.0);
                heat += report.motion.element_heat_j()[0];
                let physical = history.material_history();
                let equilibrium =
                    quadratic(&source.equilibrium().stiffness, physical.displacement_m());
                let kinetic = quadratic(&source.equilibrium().mass, physical.velocity_mps());
                let memory: f64 = body
                    .projection()
                    .material_points()
                    .iter()
                    .flat_map(|p| p.branches())
                    .flat_map(|b| b.rows())
                    .map(|row| {
                        let z: f64 = row
                            .iter()
                            .zip(physical.memory_displacement_m())
                            .map(|(a, y)| a * y)
                            .sum();
                        0.5 * z * z
                    })
                    .sum();
                let balance =
                    equilibrium + kinetic + memory + heat - initial_energy - independent_work;
                maximum_balance = maximum_balance.max(balance.abs());
                assert!(balance.abs() < 1e-10, "original solid balance {balance:e}");
            }
            let physical = history.material_history();
            assert!(heat > 0.0);
            assert!((physical.supplied_work_j() - independent_work).abs() < 2e-14);
            assert!((physical.dissipated_heat_j() - heat).abs() < 2e-14);
            println!(
                "affine TET10 pulse: subdivisions={subdivisions}, maximum_original_balance_j={maximum_balance:e}, heat_j={heat:e}, nodal_work_j={independent_work:e}"
            );
            results.push((physical.state().to_vec(), independent_work, heat));
        }
        let reference = &results[2];
        for result in &results[..2] {
            for (actual, expected) in result.0.iter().zip(&reference.0) {
                assert!(
                    (actual - expected).abs() < 2e-12,
                    "split trajectory {actual:e} != {expected:e}"
                );
            }
            assert!((result.1 - reference.1).abs() < 2e-14);
            assert!((result.2 - reference.2).abs() < 2e-14);
        }
    });
}

#[test]
fn g4_affine_solid_keeps_fixed_rows_powerless_and_rejects_partial_force_fields() {
    with_cx(|cx, gate| {
        let cells = [[0, 1, 2, 3]];
        let fixed: Vec<_> = (0..11).collect();
        let law = material(0.0002, 0.4);
        let source = ViscoElements::new(
            SolidMeshRef::Tet4 {
                nodes_m: &CORNERS,
                elements: &cells,
            },
            TetViscoField::Uniform(&law),
            &fixed,
            TetAssemblyBudget::standard(),
            ViscoStateBudget {
                max_internal_variables: 1000,
            },
            cx,
        )
        .unwrap();
        let body = PreparedTensorSolidBody::new_affine(
            source
                .project_basis(&identity(1), projection_budget(), cx)
                .unwrap(),
            1e-5,
            1.0,
            controls(),
            gate,
        )
        .unwrap();
        let mut a = body
            .history(&[1e-6], &[0.005], InitialMemory::Instantaneous, gate)
            .unwrap();
        let mut b = a.clone();
        let mut f0 = [[0.0; 3]; 4];
        let mut f1 = [[0.0; 3]; 4];
        f0[3][2] = 0.15;
        f1[3][2] = -0.2;
        let reference = b.advance_affine(&f0, &f1, gate).unwrap();
        f0[0] = [100.0; 3];
        f1[0] = [-200.0; 3];
        let report = a.advance_affine(&f0, &f1, gate).unwrap();
        assert_eq!(a.material_history().state(), b.material_history().state());
        assert_eq!(
            report.motion.material.energy,
            reference.motion.material.energy
        );
        assert_eq!(report.motion.mean_nodal_velocity_mps[0], [0.0; 3]);
        assert_eq!(report.first_nodal_velocity_moment_mps[0], [0.0; 3]);
        let snapshot =
            |h: &fs_couple::tensor_material::solid::prepared::PreparedTensorSolidHistory<
                '_,
                '_,
            >| {
                let m = h.material_history();
                let mut result: Vec<_> = m
                    .state()
                    .iter()
                    .chain(m.velocity_mps())
                    .map(|x| x.to_bits())
                    .collect();
                result.extend(
                    [
                        m.time_s(),
                        m.stored_energy_j(),
                        m.supplied_work_j(),
                        m.dissipated_heat_j(),
                    ]
                    .map(f64::to_bits),
                );
                result.push(m.steps());
                result
            };
        let before = snapshot(&a);
        assert!(a.advance_affine(&f0, &f1[..3], gate).is_err());
        assert_eq!(snapshot(&a), before);
        f1[0][0] = f64::NAN; // Even a fixed row must have valid supplied data.
        assert!(a.advance_affine(&f0, &f1, gate).is_err());
        assert_eq!(snapshot(&a), before);
    });
}
