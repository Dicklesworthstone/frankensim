use super::*;
use crate::tensor_material::InitialMemory;
use fs_material::{IsotropicElastic, visco_tensor::TensorPronyBranch};

fn projection() -> TensorProjectionBudget {
    TensorProjectionBudget {
        max_coordinates: 128,
        max_scalar_terms: 100_000_000,
    }
}
fn dynamics() -> TensorMaterialBudget {
    TensorMaterialBudget {
        max_state_coordinates: 128,
        max_material_points: 128,
    }
}
fn material(memory: bool) -> MembraneFacetMaterial {
    let c = IsotropicElastic::new(1000.0, 0.25, 0.2)
        .unwrap()
        .stiffness();
    let branches = if memory {
        vec![TensorPronyBranch::from_stiffness_mandel(c.map(|r| r.map(|v| 0.4 * v)), 0.02).unwrap()]
    } else {
        vec![]
    };
    MembraneFacetMaterial {
        law: fs_material::visco_tensor::TensorProny::from_tensor_shear(c, branches, 0.2).unwrap(),
        thickness_m: 0.1,
        density_kg_m3: 10.0,
    }
}
fn triangle() -> ShellMesh {
    ShellMesh::new(
        vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        vec![[0, 1, 2]],
    )
    .unwrap()
}
const FIXED: [usize; 6] = [0, 1, 2, 4, 5, 8];
fn near(actual: f64, expected: f64, absolute: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * expected.abs(),
        "{actual:e} != {expected:e}"
    );
}
fn prepare(density: [f64; 3], mass: MembraneMass, memory: bool) -> TensorMembrane {
    TensorMembrane::new_with_reference_areal_density(
        triangle(),
        vec![material(memory)],
        vec![density],
        &FIXED,
        mass,
        projection(),
        dynamics(),
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}

#[test]
fn constant_inventory_retains_scalar_mass_and_original_native_history() {
    let gate = CancelGate::new_clock_free();
    for selector in [MembraneMass::Lumped, MembraneMass::Consistent] {
        let original = match selector {
            MembraneMass::Lumped => TensorMembrane::new(
                triangle(),
                vec![material(true)],
                &FIXED,
                projection(),
                dynamics(),
                &gate,
            ),
            MembraneMass::Consistent => TensorMembrane::new_with_consistent_mass(
                triangle(),
                vec![material(true)],
                &FIXED,
                projection(),
                dynamics(),
                &gate,
            ),
        }
        .unwrap();
        let field = prepare([1.0; 3], selector, true);
        assert!(original.reference_corner_areal_density_kg_m2().is_none());
        assert_eq!(
            field.reference_corner_areal_density_kg_m2(),
            Some([[1.0; 3]].as_slice())
        );
        assert_eq!(original.dof_map(), field.dof_map());
        assert_eq!(
            original.body().state_dimension(),
            field.body().state_dimension()
        );
        let mass = original.body().mass_matrix_kg();
        for (a, b) in mass.iter().zip(field.body().mass_matrix_kg()) {
            near(*a, b, 0.0, 3e-15);
        }
        // The legacy preparation still executes its original scalar order.
        let share: f64 = 10.0 * (0.1 * 0.5) / 3.0;
        assert!(
            original
                .nodal_masses_kg()
                .iter()
                .all(|m| m.to_bits() == share.to_bits())
        );
        if selector == MembraneMass::Consistent {
            let weight: f64 = 10.0 * (0.1 * 0.5) / 12.0;
            assert_eq!(
                mass,
                [
                    2.0 * weight,
                    weight,
                    0.0,
                    weight,
                    2.0 * weight,
                    0.0,
                    0.0,
                    0.0,
                    2.0 * weight
                ]
            );
        }
        let zero = [[0.0; 3]; 3];
        let displacement = [[0.0; 3], [0.001, 0.0, 0.0], [0.0006, -0.0002, 0.0]];
        let mut a = original
            .history(&displacement, &zero, InitialMemory::Instantaneous, &gate)
            .unwrap();
        let mut b = field
            .history(&displacement, &zero, InitialMemory::Instantaneous, &gate)
            .unwrap();
        for _ in 0..3 {
            a.advance(&[0.01, -0.003, 0.002], 0.002, &gate).unwrap();
            b.advance(&[0.01, -0.003, 0.002], 0.002, &gate).unwrap();
            for (a, b) in a.state().iter().zip(b.state()) {
                near(*a, *b, 2e-13, 2e-12);
            }
            near(a.dissipated_heat_j(), b.dissipated_heat_j(), 2e-14, 2e-12);
        }
    }
}

// Independent degree-three triangle quadrature, including the centroid's
// negative weight. The production owner uses analytic barycentric moments.
fn kinetic_by_quadrature(density: [f64; 3], velocity: [[f64; 3]; 3]) -> f64 {
    let points = [
        ([1.0 / 3.0; 3], -27.0 / 48.0),
        ([0.6, 0.2, 0.2], 25.0 / 48.0),
        ([0.2, 0.6, 0.2], 25.0 / 48.0),
        ([0.2, 0.2, 0.6], 25.0 / 48.0),
    ];
    points
        .into_iter()
        .map(|(shape, weight)| {
            let rho: f64 = shape.iter().zip(density).map(|(n, rho)| n * rho).sum();
            let v: [f64; 3] = core::array::from_fn(|axis| {
                (0..3).map(|node| shape[node] * velocity[node][axis]).sum()
            });
            0.5 * 0.5 * weight * rho * v.iter().map(|v| v * v).sum::<f64>()
        })
        .sum()
}

#[test]
fn graded_consistent_inertia_matches_cubic_velocity_energy_and_corner_order() {
    let gate = CancelGate::new_clock_free();
    let density = [1.0, 4.0, 7.0];
    let original = prepare(density, MembraneMass::Consistent, true);
    let velocity = [[0.0; 3], [0.003, 0.0, 0.0], [-0.002, 0.004, 0.0]];
    let history = original
        .history(&[[0.0; 3]; 3], &velocity, InitialMemory::Relaxed, &gate)
        .unwrap();
    near(
        history.stored_energy_j(),
        kinetic_by_quadrature(density, velocity),
        1e-20,
        3e-14,
    );
    let mass = original.body().mass_matrix_kg();
    let v = original.project_nodal_forces(&velocity).unwrap();
    for i in 0..3 {
        let momentum: f64 = (0..3).map(|j| mass[3 * i + j] * v[j]).sum();
        near(history.state()[3 + i], momentum, 1e-19, 2e-14);
    }
    for (actual, want) in
        original
            .nodal_masses_kg()
            .iter()
            .zip([13.0 / 24.0, 16.0 / 24.0, 19.0 / 24.0])
    {
        near(*actual, want, 0.0, 2e-15);
    }
    assert!(mass[0] < original.nodal_masses_kg()[1]);
    assert!(mass[1] > 0.0);
    assert_eq!(history.dissipated_heat_j(), 0.0);
    assert_eq!(history.parameter_work_j(), 0.0);

    // Reversing winding with its density corner map preserves the same global
    // inertia; a reference area is positive, not an oriented mass measure.
    let mut reversed = triangle();
    reversed.tris[0] = [0, 2, 1];
    let reversed = TensorMembrane::new_with_reference_areal_density(
        reversed,
        vec![material(true)],
        vec![[density[0], density[2], density[1]]],
        &FIXED,
        MembraneMass::Consistent,
        projection(),
        dynamics(),
        &gate,
    )
    .unwrap();
    for (a, b) in mass.into_iter().zip(reversed.body().mass_matrix_kg()) {
        near(a, b, 0.0, 2e-15);
    }
    for (a, b) in original
        .nodal_masses_kg()
        .iter()
        .zip(reversed.nodal_masses_kg())
    {
        near(*a, *b, 0.0, 2e-15);
    }
}

#[test]
fn equal_total_inventory_has_distinct_native_acceleration_and_oscillation() {
    let gate = CancelGate::new_clock_free();
    // Only node 1's x displacement is free: the original plane-stress owner
    // reduces to one independently known spring, without modal projection.
    let fixed = [0, 1, 2, 4, 5, 6, 7, 8];
    let stiffness = 1000.0 / (1.0 - 0.25_f64.powi(2)) * 0.1 * 0.5;
    let force = 0.07;
    let dt = 0.03;
    for selector in [MembraneMass::Lumped, MembraneMass::Consistent] {
        let mut responses = Vec::new();
        for rho in [[1.0, 4.0, 1.0], [4.0, 1.0, 1.0]] {
            let mesh = TensorMembrane::new_with_reference_areal_density(
                triangle(),
                vec![material(false)],
                vec![rho],
                &fixed,
                selector,
                projection(),
                dynamics(),
                &gate,
            )
            .unwrap();
            near(mesh.nodal_masses_kg().iter().sum(), 1.0, 0.0, 2e-15);
            let inertia = match selector {
                MembraneMass::Lumped => 0.5 * (6.0 + rho[1]) / 12.0,
                MembraneMass::Consistent => 0.5 * (6.0 + 2.0 * rho[1]) / 30.0,
            };
            near(mesh.body().mass_matrix_kg()[0], inertia, 0.0, 2e-15);
            near(
                mesh.body().plane_stress_points()[0].energy_matrix()[0],
                stiffness,
                0.0,
                2e-14,
            );
            let mut history = mesh
                .history(
                    &[[0.0; 3]; 3],
                    &[[0.0; 3]; 3],
                    InitialMemory::Relaxed,
                    &gate,
                )
                .unwrap();
            history.advance(&[force], dt, &gate).unwrap();
            let omega = (stiffness / inertia).sqrt();
            let q = force / stiffness * (1.0 - (omega * dt).cos());
            let v = force / (inertia * omega) * (omega * dt).sin();
            near(history.displacement_m()[0], q, 2e-14, 3e-11);
            near(history.velocity_m_s()[0], v, 2e-13, 3e-11);
            near(history.state()[1], inertia * v, 2e-14, 3e-11);
            near(
                history.stored_energy_j(),
                0.5 * (inertia * v * v + stiffness * q * q),
                1e-16,
                5e-11,
            );
            near(history.supplied_work_j(), force * q, 1e-16, 5e-11);
            assert_eq!(history.dissipated_heat_j(), 0.0);
            assert_eq!(history.parameter_work_j(), 0.0);
            responses.push(q);
        }
        assert!((responses[0] - responses[1]).abs() > 1e-5);
    }
}

#[test]
fn graded_inventory_keeps_full_hereditary_work_and_material_heat() {
    let gate = CancelGate::new_clock_free();
    let rho = [1.0, 3.0, 8.0];
    let mut source = material(true);
    let mut gram = [[0.0; 6]; 6];
    gram[0] = [8.0, 2.0, 1.0, 3.0, 0.0, 0.0];
    source.law = fs_material::visco_tensor::TensorProny::new_mandel(
        *source.law.equilibrium_mandel(),
        vec![TensorPronyBranch::from_gram_rows(gram, 0.02).unwrap()],
        0.2,
    )
    .unwrap();
    let mesh = TensorMembrane::new_with_reference_areal_density(
        triangle(),
        vec![source],
        vec![rho],
        &FIXED,
        MembraneMass::Consistent,
        projection(),
        dynamics(),
        &gate,
    )
    .unwrap();
    let velocity = [[0.0; 3], [0.003, 0.0, 0.0], [-0.002, 0.004, 0.0]];
    let mut history = mesh
        .history(&[[0.0; 3]; 3], &velocity, InitialMemory::Relaxed, &gate)
        .unwrap();
    let initial = kinetic_by_quadrature(rho, velocity);
    near(history.stored_energy_j(), initial, 1e-20, 3e-14);
    let mut work = 0.0;
    for tick in 0..12 {
        let force = if tick < 4 {
            [0.01, -0.003, 0.002]
        } else {
            [0.0; 3]
        };
        let old = history.displacement_m().to_vec();
        history.advance(&force, 0.002, &gate).unwrap();
        work += force
            .iter()
            .zip(history.displacement_m())
            .zip(old)
            .map(|((f, new), old)| f * (new - old))
            .sum::<f64>();
        near(history.supplied_work_j(), work, 3e-14, 3e-11);
        near(
            history.stored_energy_j() + history.dissipated_heat_j(),
            initial + work,
            3e-14,
            3e-11,
        );
        for traction in history.plane_stress_normal_stress_pa(0).unwrap() {
            near(traction, 0.0, 2e-10, 0.0);
        }
        assert_eq!(history.parameter_work_j(), 0.0);
        assert_eq!(
            mesh.reference_corner_areal_density_kg_m2(),
            Some([rho].as_slice())
        );
    }
    assert!(history.dissipated_heat_j() > 0.0);
    let saved = history.state().to_vec();
    let time = history.time_s();
    let heat = history.dissipated_heat_j();
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(history.advance(&[0.01; 3], 0.002, &cancelled).is_err());
    assert_eq!(history.state(), saved);
    assert_eq!(history.time_s(), time);
    assert_eq!(history.dissipated_heat_j(), heat);
}

#[test]
fn incident_density_jumps_retain_native_rigid_translation_and_mass_moments() {
    let gate = CancelGate::new_clock_free();
    let geometry = ShellMesh::new(
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 0.0],
        ],
        vec![[0, 1, 2], [0, 2, 3]],
    )
    .unwrap();
    // Different incident corner values at nodes 0 and 2 remain distinct.
    let rho = vec![[1.0, 2.0, 3.0], [10.0, 11.0, 12.0]];
    let body = TensorMembrane::new_with_reference_areal_density(
        geometry,
        vec![material(true); 2],
        rho.clone(),
        &[2, 5, 8, 11],
        MembraneMass::Consistent,
        projection(),
        dynamics(),
        &gate,
    )
    .unwrap();
    let expected = [50.0 / 24.0, 8.0 / 24.0, 53.0 / 24.0, 45.0 / 24.0];
    for (actual, want) in body.nodal_masses_kg().iter().zip(expected) {
        near(*actual, want, 0.0, 2e-15);
    }
    let total: f64 = body.nodal_masses_kg().iter().sum();
    near(total, 6.5, 0.0, 2e-15);
    let acceleration = [0.2, -0.1, 0.0];
    let nodal = body
        .distributed_nodal_forces(&[[[0.0; 3]; 3]; 2], &[], acceleration, &gate)
        .unwrap();
    for (force, mass) in nodal.iter().zip(expected) {
        for axis in 0..3 {
            near(force[axis], mass * acceleration[axis], 0.0, 3e-15);
        }
    }
    let initial_velocity = [0.03, -0.02, 0.0];
    let mut history = body
        .history(
            &[[0.0; 3]; 4],
            &[initial_velocity; 4],
            InitialMemory::Relaxed,
            &gate,
        )
        .unwrap();
    let initial = 0.5 * total * initial_velocity.iter().map(|v| v * v).sum::<f64>();
    let dt = 0.01;
    body.advance_nodal(&mut history, &nodal, dt, &gate).unwrap();
    let displacement = body.expand_coordinates(history.displacement_m()).unwrap();
    let velocity = body.expand_coordinates(&history.velocity_m_s()).unwrap();
    for node in 0..4 {
        for axis in 0..3 {
            near(
                displacement[node][axis],
                initial_velocity[axis] * dt + 0.5 * acceleration[axis] * dt * dt,
                2e-12,
                2e-11,
            );
            near(
                velocity[node][axis],
                initial_velocity[axis] + acceleration[axis] * dt,
                2e-11,
                2e-10,
            );
        }
    }
    let final_kinetic = 0.5
        * total
        * (0..3)
            .map(|axis| (initial_velocity[axis] + acceleration[axis] * dt).powi(2))
            .sum::<f64>();
    near(history.stored_energy_j(), final_kinetic, 2e-12, 2e-10);
    near(
        history.supplied_work_j(),
        final_kinetic - initial,
        2e-12,
        2e-10,
    );
    near(history.dissipated_heat_j(), 0.0, 2e-12, 0.0);
    assert_eq!(
        body.reference_corner_areal_density_kg_m2(),
        Some(rho.as_slice())
    );
}

#[test]
fn shared_moments_admit_contrast_and_refuse_incomplete_or_unrepresentable_fields() {
    // Summing these densities first overflows, although every physical moment
    // is ordinary-sized. Scaling is shared with the original solid membrane.
    let rho = [1e308, 1e308, 1e-308];
    let area = 1e-308;
    let rows = p1_reference_row_masses(area, rho).unwrap();
    let matrix = p1_reference_consistent_mass(area, rho).unwrap();
    for (actual, want) in rows.into_iter().zip([0.25, 0.25, 1.0 / 6.0]) {
        near(actual, want, 0.0, 4e-15);
    }
    for i in 0..3 {
        near(matrix[i].iter().sum(), rows[i], 0.0, 4e-15);
        for j in 0..3 {
            assert_eq!(matrix[i][j], matrix[j][i]);
            assert!(matrix[i][j] > 0.0);
        }
    }
    for invalid in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(p1_reference_row_masses(invalid, [1.0; 3]).is_err());
        assert!(p1_reference_consistent_mass(invalid, [1.0; 3]).is_err());
        assert!(p1_reference_row_masses(0.5, [1.0, invalid, 2.0]).is_err());
        assert!(p1_reference_consistent_mass(0.5, [1.0, invalid, 2.0]).is_err());
    }
    assert!(p1_reference_row_masses(0.5, [f64::from_bits(1); 3]).is_err());
    assert!(p1_reference_consistent_mass(0.5, [f64::from_bits(1); 3]).is_err());
    assert!(p1_reference_consistent_mass(f64::MAX, [f64::MAX, 1.0, 2.0]).is_err());

    let gate = CancelGate::new_clock_free();
    let build = |density, gate: &CancelGate| {
        TensorMembrane::new_with_reference_areal_density(
            triangle(),
            vec![material(false)],
            density,
            &FIXED,
            MembraneMass::Consistent,
            projection(),
            dynamics(),
            gate,
        )
    };
    assert!(build(vec![], &gate).is_err());
    assert!(build(vec![[1.0; 3]; 2], &gate).is_err());
    for invalid in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::from_bits(1)] {
        assert!(build(vec![[invalid; 3]], &gate).is_err());
    }
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        build(vec![[1.0; 3]], &cancelled),
        Err(MembraneError::Dynamics(
            crate::tensor_material::TensorMaterialError::Cancelled
        ))
    ));
    let contrast = build(vec![[1e-12, 1.0, 1e12]], &gate).unwrap();
    assert!(
        contrast
            .body()
            .mass_matrix_kg()
            .iter()
            .all(|m| m.is_finite())
    );
    near(
        contrast.nodal_masses_kg().iter().sum(),
        0.5 * (1e12 + 1.0 + 1e-12) / 3.0,
        0.0,
        3e-15,
    );
}
