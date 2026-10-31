//! Closed distributed current, gauge/polarity and actual native 3D circuit work.
//! Manufactured linear cards test the equations, not measured coil fidelity.
mod common;
use common::{material, near};
use fs_em::spatial::{
    SpatialMagneticError, SpatialMagneticMaterial,
    vector_potential::{
        SpatialVectorMagnetostatics, VectorMagneticBoundary, VectorMagneticControls,
        VectorPotentialFace,
        winding::{SpatialWindingBasis, SpatialWindingConfig, SpatialWindingError},
    },
};
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_solver::CheckedCgConfig;

fn controls() -> VectorMagneticControls {
    let solve = CheckedCgConfig {
        tolerance: 1e-12,
        max_iterations: 2000,
        max_corrections: 2,
    };
    VectorMagneticControls {
        field: solve,
        gauge: solve,
        maximum_gauge_product_terms: 2_000_000,
        current_tolerance_a: 2e-8,
        flux_tolerance_wb: 2e-9,
        energy_tolerance_j: 3e-8,
    }
}
fn config() -> SpatialWindingConfig {
    SpatialWindingConfig {
        max_windings: 4,
        max_basis_values: 100_000,
        pairing_tolerance_h: 2e-9,
    }
}
fn cube(mu: f64) -> SpatialVectorMagnetostatics {
    let mut positions = Vec::new();
    for k in 0..=2 {
        for j in 0..=2 {
            for i in 0..=2 {
                positions.push([f64::from(i) / 2.0, f64::from(j) / 2.0, f64::from(k) / 2.0]);
            }
        }
    }
    let id = |i, j, k| ((k * 3 + j) * 3 + i) as u32;
    let mut tets = Vec::new();
    for k in 0..2 {
        for j in 0..2 {
            for i in 0..2 {
                let (a, b, c, d) = (
                    id(i, j, k),
                    id(i + 1, j, k),
                    id(i, j + 1, k),
                    id(i + 1, j + 1, k),
                );
                let (e, f, g, h) = (
                    id(i, j, k + 1),
                    id(i + 1, j, k + 1),
                    id(i, j + 1, k + 1),
                    id(i + 1, j + 1, k + 1),
                );
                tets.extend([
                    [a, b, d, h],
                    [a, d, c, h],
                    [a, c, g, h],
                    [a, g, e, h],
                    [a, e, f, h],
                    [a, f, b, h],
                ]);
            }
        }
    }
    let materials = vec![SpatialMagneticMaterial::Linear(material(mu)); tets.len()];
    SpatialVectorMagnetostatics::new(
        positions,
        tets,
        materials,
        WeightedAssemblyLimits {
            max_cells: 1000,
            max_dofs: 4000,
            max_triplets: 100_000,
        },
        &mut || false,
    )
    .unwrap()
}
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|k| a[k] - b[k])
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

// Independently derive the central P1 hat gradient from its opposite face.
// J_unit=2 curl(N e_z) carries exactly -1 A through x=1/2,y>1/2:
// its integral is -2 integral N(1/2,1/2,z) dz = -1. The second winding
// additionally carries one identically normalized x-axis loop in series.
fn profiles(problem: &SpatialVectorMagnetostatics) -> Vec<Vec<[f64; 3]>> {
    let origin = problem
        .positions_m()
        .iter()
        .position(|&p| p == [0.5; 3])
        .unwrap() as u32;
    let mut profiles = vec![Vec::new(), Vec::new()];
    for tet in &problem.mesh().tets {
        let gradient = if tet.contains(&origin) {
            let face: Vec<_> = tet.iter().copied().filter(|&v| v != origin).collect();
            let p = problem.positions_m()[face[0] as usize];
            let normal = cross(
                sub(problem.positions_m()[face[1] as usize], p),
                sub(problem.positions_m()[face[2] as usize], p),
            );
            let height = dot(normal, sub(problem.positions_m()[origin as usize], p));
            normal.map(|v| v / height)
        } else {
            [0.0; 3]
        };
        profiles[0].push(cross(gradient, [0.0, 0.0, 2.0]));
        profiles[1].push(cross(gradient, [2.0, 0.0, 2.0]));
    }
    profiles
}
fn boundary(problem: &SpatialVectorMagnetostatics) -> VectorMagneticBoundary {
    VectorMagneticBoundary {
        potential: problem
            .exterior_faces()
            .into_iter()
            .map(|face| VectorPotentialFace {
                face,
                edge_integrals_wb: [0.0; 3],
            })
            .collect(),
    }
}
fn basis<'a>(
    problem: &'a SpatialVectorMagnetostatics,
    p: Vec<Vec<[f64; 3]>>,
    ends: Vec<f64>,
) -> SpatialWindingBasis<'a> {
    SpatialWindingBasis::new(
        problem,
        p,
        boundary(problem),
        ends,
        config(),
        controls(),
        &mut || false,
    )
    .unwrap()
}
fn sum_current(profiles: &[Vec<[f64; 3]>], currents: &[f64]) -> Vec<[f64; 3]> {
    (0..profiles[0].len())
        .map(|cell| {
            std::array::from_fn(|k| {
                profiles
                    .iter()
                    .zip(currents)
                    .map(|(p, i)| p[cell][k] * i)
                    .sum()
            })
        })
        .collect()
}

#[test]
fn closed_current_pairing_gauge_and_polarity_reconstruct_the_original_native_field() {
    let problem = cube(2.0);
    let p = profiles(&problem);
    let mut cut_current = 0.0;
    for face in &problem.mesh().faces {
        let points = face.map(|v| problem.positions_m()[v as usize]);
        if points.iter().all(|v| v[0] == 0.5)
            && points.iter().map(|v| v[1]).sum::<f64>() / 3.0 > 0.5
        {
            let area = 0.5 * cross(sub(points[1], points[0]), sub(points[2], points[0]))[0].abs();
            let cell = problem
                .mesh()
                .tets
                .iter()
                .position(|tet| face.iter().all(|v| tet.contains(v)))
                .unwrap();
            cut_current += p[0][cell][0] * area;
        }
    }
    near(cut_current, -1.0);
    let ends = vec![0.02, 0.03];
    let winding = basis(&problem, p.clone(), ends.clone());
    let currents = [0.7, -0.2];
    let reconstructed = winding
        .field_for_currents(&currents, &mut || false)
        .unwrap();
    let direct = problem
        .solve(
            &sum_current(&p, &currents),
            &boundary(&problem),
            controls(),
            &mut || false,
        )
        .unwrap();
    for (a, b) in reconstructed
        .flux_density_t
        .iter()
        .flatten()
        .zip(direct.flux_density_t().iter().flatten())
    {
        near(*a, *b);
    }
    for (a, b) in reconstructed
        .field_strength_a_m
        .iter()
        .flatten()
        .zip(direct.field_strength_a_m().iter().flatten())
    {
        near(*a, *b);
    }
    near(reconstructed.field_energy_j, direct.stored_energy_j());
    near(
        reconstructed.external_energy_j,
        0.5 * (ends[0] * currents[0].powi(2) + ends[1] * currents[1].powi(2)),
    );
    near(reconstructed.energy_residual_j, 0.0);
    assert!(reconstructed.maximum_divergence_flux_wb <= controls().flux_tolerance_wb);
    assert!(reconstructed.maximum_recovered_flux_defect_wb <= controls().flux_tolerance_wb);
    assert!(winding.maximum_pairing_defect_h() <= config().pairing_tolerance_h);
    let l = winding.inductance_h();
    near(l[1], l[2]);
    assert!(l[1].abs() > 1e-5);
    let paired = winding
        .pair_potential_wb(&reconstructed.potential_wb, &mut || false)
        .unwrap();
    let psi = |v: u32| f64::from((v * 7) % 11) * 0.17;
    let shifted: Vec<_> = problem
        .mesh()
        .edges
        .iter()
        .zip(&reconstructed.potential_wb)
        .map(|(&[a, b], &value)| value + psi(b) - psi(a))
        .collect();
    let gauge = winding.pair_potential_wb(&shifted, &mut || false).unwrap();
    for i in 0..2 {
        near(paired[i], gauge[i]);
        near(
            reconstructed.flux_linkages_wb[i],
            paired[i] + ends[i] * currents[i],
        );
    }
    let mut reversed = p;
    for j in &mut reversed[1] {
        *j = j.map(|v| -v);
    }
    let reverse = basis(&problem, reversed, ends);
    for (index, sign) in [1.0, -1.0, -1.0, 1.0].into_iter().enumerate() {
        near(reverse.inductance_h()[index], sign * l[index]);
    }
    let reverse_field = reverse
        .field_for_currents(&[currents[0], -currents[1]], &mut || false)
        .unwrap();
    for (a, b) in reverse_field
        .flux_density_t
        .iter()
        .flatten()
        .zip(reconstructed.flux_density_t.iter().flatten())
    {
        near(*a, *b);
    }

    let changed = cube(4.0);
    let changed_basis = basis(&changed, profiles(&changed), vec![0.02, 0.03]);
    for (&a, &b) in changed_basis
        .active_inductance_h()
        .iter()
        .zip(winding.active_inductance_h())
    {
        near(a, 2.0 * b);
    }
}

fn graph() -> CircuitGraph {
    CircuitGraph {
        node_count: 6,
        branches: vec![
            (1, 0, Branch::VoltageSource { port: 0 }),
            (1, 2, Branch::Resistor { ohms: 0.7 }),
            (2, 0, Branch::Inductor { henries: f64::NAN }),
            (3, 0, Branch::VoltageSource { port: 2 }),
            (3, 4, Branch::Resistor { ohms: 1.1 }),
            (4, 0, Branch::Inductor { henries: f64::NAN }),
            (5, 0, Branch::Inductor { henries: 0.4 }),
            (1, 5, Branch::Resistor { ohms: 0.6 }),
            (5, 0, Branch::Capacitor { farads: 0.3 }),
        ],
        transformers: vec![],
    }
}

#[test]
fn full_3d_mutual_matrix_changes_actual_native_rlc_currents_and_work() {
    let problem = cube(2.0);
    let winding = basis(&problem, profiles(&problem), vec![0.02, 0.03]);
    let l = winding.inductance_h().to_vec();
    // Preserve nonascending original branch bindings and the gap at source 1.
    let circuit = winding
        .into_circuit(graph(), vec![5, 2], &mut || false)
        .unwrap();
    let inputs = [1.3, 17.0, -0.6];
    let mut currents = [0.2, -0.15];
    let mut state = circuit
        .consistent_initial_state(&currents, &[0.08], &[0.06], &inputs, &mut || false)
        .unwrap();
    assert_eq!(circuit.branch_indices(), &[5, 2]);
    let read = circuit.currents_a(&state).unwrap();
    near(read[0], currents[0]);
    near(read[1], currents[1]);
    near(state[circuit.dae().flux_index[2]], 0.08);
    near(state[circuit.dae().charge_index[0]], 0.06);
    let dt = 0.005;
    for _ in 0..6 {
        // Independent 2x2 trapezoidal RL equations, in declared winding order.
        let r = [1.1, 0.7];
        let voltage = [-0.6, 1.3];
        let a = l[0] + 0.5 * dt * r[0];
        let d = l[3] + 0.5 * dt * r[1];
        let b0 = (l[0] - 0.5 * dt * r[0]) * currents[0] + l[1] * currents[1] + dt * voltage[0];
        let b1 = l[2] * currents[0] + (l[3] - 0.5 * dt * r[1]) * currents[1] + dt * voltage[1];
        let determinant = a * d - l[1] * l[2];
        let expected = [
            (d * b0 - l[1] * b1) / determinant,
            (a * b1 - l[2] * b0) / determinant,
        ];
        let before = circuit.dae().system.hamiltonian(&state);
        let (step, _) = circuit
            .step_audited(&state, &inputs, dt, &mut || false)
            .unwrap();
        let actual = circuit.currents_a(&step.x).unwrap();
        for i in 0..2 {
            near(actual[i], expected[i]);
        }
        let middle: Vec<_> = state
            .iter()
            .zip(&step.x)
            .map(|(&a, &b)| f64::midpoint(a, b))
            .collect();
        let potential = |node: usize| {
            if node == 0 {
                0.0
            } else {
                middle[circuit.dae().node_potential_index[node - 1]]
            }
        };
        let heat: f64 = circuit
            .graph()
            .branches
            .iter()
            .map(|&(a, b, branch)| match branch {
                Branch::Resistor { ohms } => dt * (potential(a) - potential(b)).powi(2) / ohms,
                _ => 0.0,
            })
            .sum();
        let supply: f64 = dt * inputs.iter().zip(&step.y).map(|(u, y)| u * y).sum::<f64>();
        let after = circuit.dae().system.hamiltonian(&step.x);
        near(step.dissipated, heat);
        near(after - before + heat, supply);
        near(step.balance_residual(), 0.0);
        currents = [actual[0], actual[1]];
        state = step.x;
    }
    let field = circuit.field_for_state(&state, &mut || false).unwrap();
    let other_flux = state[circuit.dae().flux_index[2]];
    let charge = state[circuit.dae().charge_index[0]];
    near(
        circuit.dae().system.hamiltonian(&state),
        field.stored_energy_j + 0.5 * other_flux * other_flux / 0.4 + 0.5 * charge * charge / 0.3,
    );
    for (winding, &coordinate) in circuit.flux_coordinates().iter().enumerate() {
        near(state[coordinate], field.flux_linkages_wb[winding]);
    }
}

#[test]
fn current_return_boundary_budget_and_dependent_modes_refuse_without_regularization() {
    let problem = cube(2.0);
    let p = profiles(&problem);
    let open = vec![vec![[1.0, 0.0, 0.0]; problem.mesh().tets.len()]];
    assert!(
        SpatialWindingBasis::new(
            &problem,
            open,
            boundary(&problem),
            vec![0.0],
            config(),
            controls(),
            &mut || false
        )
        .is_err()
    );
    let mut biased = boundary(&problem);
    biased.potential[0].edge_integrals_wb[0] = 0.1;
    assert!(
        SpatialWindingBasis::new(
            &problem,
            p.clone(),
            biased,
            vec![0.0; 2],
            config(),
            controls(),
            &mut || false
        )
        .is_err()
    );
    let mut limited = config();
    limited.max_basis_values = 1;
    assert!(
        SpatialWindingBasis::new(
            &problem,
            p.clone(),
            boundary(&problem),
            vec![0.0; 2],
            limited,
            controls(),
            &mut || false
        )
        .is_err()
    );
    assert!(matches!(
        SpatialWindingBasis::new(
            &problem,
            p.clone(),
            boundary(&problem),
            vec![0.0; 2],
            config(),
            controls(),
            &mut || true
        ),
        Err(SpatialWindingError::Magnetic(
            SpatialMagneticError::Cancelled
        ))
    ));
    let duplicate = basis(&problem, vec![p[0].clone(), p[0].clone()], vec![0.0; 2]);
    assert!(
        duplicate
            .into_circuit(graph(), vec![2, 5], &mut || false)
            .is_err()
    );
    let explicit = basis(&problem, vec![p[0].clone(), p[0].clone()], vec![0.02, 0.03]);
    let circuit = explicit
        .into_circuit(graph(), vec![2, 5], &mut || false)
        .unwrap();
    let inputs = [0.0; 3];
    let state = circuit
        .consistent_initial_state(&[0.2, -0.1], &[0.0], &[0.0], &inputs, &mut || false)
        .unwrap();
    let saved = state.clone();
    let mut calls = 0;
    let first = circuit
        .step_audited(&state, &inputs, 0.01, &mut || {
            calls += 1;
            false
        })
        .unwrap();
    let mut late = 0;
    assert!(matches!(
        circuit.step_audited(&state, &inputs, 0.01, &mut || {
            late += 1;
            late == calls
        }),
        Err(SpatialWindingError::Magnetic(
            SpatialMagneticError::Cancelled
        ))
    ));
    let retry = circuit
        .step_audited(&state, &inputs, 0.01, &mut || false)
        .unwrap();
    assert_eq!(first.0.x, retry.0.x);
    assert_eq!(state, saved);
}
