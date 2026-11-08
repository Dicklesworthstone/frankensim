//! Physical modal motion and exact warped-reference inertia for native Hex8 plasticity.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::plastic::rate::LinearOverstressJ2;
use fs_material::plastic::{J2Plasticity, J2State};
use fs_material::{IsotropicElastic, SmallStrainLaw};
use fs_solid::plastic3::hex::{HexAssemblyBudget, SmallStrainHexProblem};
use fs_solid::{Formulation, IncrementSettings, SmallStrainLawField, SmallStrainTetError};

const NODES: [[f64; 3]; 8] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [1.0, 1.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
    [1.0, 0.0, 1.0],
    [1.0, 1.0, 1.0],
    [0.0, 1.0, 1.0],
];
const CELLS: [[usize; 8]; 1] = [[0, 1, 2, 3, 4, 5, 6, 7]];
const E: f64 = 1000.0;
const NU: f64 = 0.3;
const Y: f64 = 1.0;
const H: f64 = 100.0;
const DT: f64 = 0.1;

fn metal() -> J2Plasticity {
    J2Plasticity::new(IsotropicElastic::new(E, NU, 0.05).unwrap(), Y, H).unwrap()
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        actual.is_finite() && (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, tolerance {tolerance:e}"
    );
}
fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 19,
                kernel_id: 94,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ));
    });
}
fn problem<'a, L: SmallStrainLaw>(
    nodes: &'a [[f64; 3]],
    law: &'a L,
    formulation: Formulation,
    fixed: &'a [(usize, f64)],
    forces: &'a [(usize, f64)],
) -> SmallStrainHexProblem<'a, L> {
    SmallStrainHexProblem {
        nodes_m: nodes,
        hexahedra: &CELLS,
        laws: SmallStrainLawField::Uniform(law),
        formulation,
        strain_limit: 0.05,
        rotation_limit_rad: 0.05,
        prescribed_m: fixed,
        nodal_forces_n: forces,
        load_path: &[1.0],
        budget: HexAssemblyBudget::default(),
        settings: IncrementSettings::default(),
    }
}

fn shear_motion<L: SmallStrainLaw<State = J2State>>(law: &L, viscosity: f64, cx: &Cx<'_>) {
    let (rho, gamma) = (2.0, 0.008);
    let mu = E / (2.0 * (1.0 + NU));
    let alpha = (3.0_f64.sqrt() * mu * gamma - Y) / (3.0 * mu + H + viscosity / DT);
    let permanent = 3.0_f64.sqrt() * alpha;
    let tau = mu * (gamma - permanent);
    // u_x=gamma*z is one uniform shear mode. Its consistent generalized
    // mass is integral rho*z^2 = rho/3, not a lumped rho/2.
    let mass = rho / 3.0;
    let total_force = mass * gamma / (DT * DT) + tau;
    let fixed: Vec<_> = (0..24)
        .filter(|d| d % 3 != 0 || NODES[d / 3][2] == 0.0)
        .map(|d| (d, 0.0))
        .collect();
    let forces: Vec<_> = (4..8).map(|node| (3 * node, total_force / 4.0)).collect();
    for formulation in [Formulation::Standard, Formulation::BBar] {
        let states = vec![law.initial_state(); 8];
        let saved = states.clone();
        let mut native = problem(&NODES, law, formulation, &fixed, &forces);
        native.settings.max_newton_iterations = 1;
        assert!(matches!(
            native.advance_backward_euler(&[0.0; 24], &[[0.0; 3]; 8], &states, &[rho], DT, cx),
            Err(SmallStrainTetError::NewtonStalled { .. })
        ));
        assert_eq!(states, saved);
        native.settings = IncrementSettings::default();
        let mut step = native
            .advance_backward_euler(&[0.0; 24], &[[0.0; 3]; 8], &states, &[rho], DT, cx)
            .unwrap();
        assert_eq!(step.solution.quadrature_states.len(), 8);
        for ((state, stress), strain) in step
            .solution
            .quadrature_states
            .iter()
            .zip(&step.solution.quadrature_stress_pa)
            .zip(&step.solution.quadrature_strain)
        {
            near(state.alpha, alpha, 2e-11);
            near(state.plastic_strain[5], permanent / 2.0, 2e-11);
            near(strain[5], gamma / 2.0, 2e-11);
            near(stress[5], tau, 2e-8);
        }
        for (node, u) in step.solution.displacement_m.iter().enumerate() {
            near(u[0], gamma * NODES[node][2], 2e-11);
            near(step.velocity_mps[node][0], u[0] / DT, 2e-11);
        }
        near(step.energy.reference_mass_kg, rho, 2e-13);
        near(
            step.energy.kinetic_after_j,
            0.5 * mass * (gamma / DT).powi(2),
            2e-12,
        );
        near(
            step.energy.numerical_loss_j,
            step.energy.kinetic_after_j,
            2e-12,
        );
        near(
            total_force * gamma - tau * gamma,
            step.energy.inertial_work_j,
            2e-12,
        );
        near(
            step.energy.inertial_work_j,
            step.energy.kinetic_change_j + step.energy.numerical_loss_j,
            2e-12,
        );

        // The published rate-dependent stress must belong to this ONE
        // interval; recomputing with the newly committed history relaxes twice.
        if viscosity > 0.0 {
            near(
                3.0_f64.sqrt() * tau - Y - H * alpha,
                viscosity * alpha / DT,
                2e-8,
            );
        }

        // Select physical forces that unload to the permanent strain, then
        // stop there. The old point histories must survive both intervals.
        let (mut old_gamma, mut old_velocity) = (gamma, gamma / DT);
        for target in [permanent, permanent] {
            let velocity = (target - old_gamma) / DT;
            let force = mass * (velocity - old_velocity) / DT;
            let forces: Vec<_> = (4..8).map(|node| (3 * node, force / 4.0)).collect();
            let next = problem(&NODES, law, formulation, &fixed, &forces)
                .advance_backward_euler(
                    step.solution.displacement_m.as_flattened(),
                    &step.velocity_mps,
                    &step.solution.quadrature_states,
                    &[rho],
                    DT,
                    cx,
                )
                .unwrap();
            for (state, stress) in next
                .solution
                .quadrature_states
                .iter()
                .zip(&next.solution.quadrature_stress_pa)
            {
                near(state.alpha, alpha, 2e-11);
                near(state.plastic_strain[5], permanent / 2.0, 2e-11);
                for value in stress {
                    near(*value, 0.0, 2e-8);
                }
            }
            for (node, u) in next.solution.displacement_m.iter().enumerate() {
                near(u[0], target * NODES[node][2], 2e-11);
                near(next.velocity_mps[node][0], velocity * NODES[node][2], 3e-10);
            }
            near(
                force * (target - old_gamma),
                next.energy.inertial_work_j,
                2e-12,
            );
            near(
                next.energy.inertial_work_j,
                next.energy.kinetic_change_j + next.energy.numerical_loss_j,
                2e-12,
            );
            old_gamma = target;
            old_velocity = velocity;
            step = next;
        }
        near(step.energy.kinetic_after_j, 0.0, 2e-12);
    }
}

#[test]
fn standard_and_bbar_motion_retain_plastic_and_timed_overstress_history() {
    with_cx(|cx| {
        shear_motion(&metal(), 0.0, cx);
        let timed = LinearOverstressJ2::new(metal(), 20.0, DT).unwrap();
        shear_motion(&timed, 20.0, cx);
    });
}

// For x=X+a*Y*Z, y=Y+b*X*Z, z=Z+c*X*Y on [0,1]^3,
// det(J)=1-bc*X^2-ac*Y^2-ab*Z^2+2abc*X*Y*Z.
// Density is 2+3X. These exact polynomial moments independently integrate
// each nodal mass row. A two-point rule misses the quartic X terms.
const WARP: [f64; 3] = [0.12, 0.08, 0.10];
fn reference_nodal_mass(corner: [f64; 3]) -> f64 {
    let [a, b, c] = WARP;
    let moment = |axis: usize, power: usize| {
        if corner[axis] == 1.0 {
            1.0 / (power + 2) as f64
        } else {
            1.0 / ((power + 1) * (power + 2)) as f64
        }
    };
    [
        (1.0, [0, 0, 0]),
        (-b * c, [2, 0, 0]),
        (-a * c, [0, 2, 0]),
        (-a * b, [0, 0, 2]),
        (2.0 * a * b * c, [1, 1, 1]),
    ]
    .iter()
    .map(|(coefficient, p)| {
        coefficient
            * (2.0 * moment(0, p[0]) + 3.0 * moment(0, p[0] + 1))
            * moment(1, p[1])
            * moment(2, p[2])
    })
    .sum()
}

#[test]
fn warped_graded_inventory_drives_free_and_prescribed_rigid_motion() {
    with_cx(|cx| {
        let [a, b, c] = WARP;
        let nodes = NODES.map(|[x, y, z]| [x + a * y * z, y + b * x * z, z + c * x * y]);
        let profile = [NODES.map(|x| 2.0 + 3.0 * x[0])];
        let row_mass = NODES.map(reference_nodal_mass);
        let mass: f64 = row_mass.iter().sum();
        let law = metal();
        let states = vec![law.initial_state(); 8];
        let (dt, translation, velocity, acceleration) = (
            0.02,
            [0.02, -0.03, 0.04],
            [0.1, -0.2, 0.05],
            [0.2, -0.1, 0.3],
        );
        let final_velocity = std::array::from_fn::<_, 3, _>(|j| velocity[j] + dt * acceleration[j]);
        let target = std::array::from_fn::<_, 3, _>(|j| translation[j] + dt * final_velocity[j]);
        for formulation in [Formulation::Standard, Formulation::BBar] {
            for constrained in [false, true] {
                let fixed: Vec<_> = if constrained {
                    (0..24).map(|d| (d, target[d % 3])).collect()
                } else {
                    vec![]
                };
                let forces: Vec<_> = if constrained {
                    vec![]
                } else {
                    row_mass
                        .iter()
                        .enumerate()
                        .flat_map(|(node, m)| {
                            (0..3).map(move |j| (3 * node + j, m * acceleration[j]))
                        })
                        .collect()
                };
                let mut native = problem(&nodes, &law, formulation, &fixed, &forces);
                native.budget.maximum_tangent_entries = 768;
                let step = native
                    .advance_backward_euler_with_corner_density(
                        [translation; 8].as_flattened(),
                        &[velocity; 8],
                        &states,
                        &profile,
                        dt,
                        cx,
                    )
                    .unwrap();
                near(step.energy.reference_mass_kg, mass, 2e-13);
                assert_eq!(step.solution.quadrature_states, states);
                for (u, v) in step.solution.displacement_m.iter().zip(&step.velocity_mps) {
                    for j in 0..3 {
                        near(u[j], target[j], 3e-12);
                        near(v[j], final_velocity[j], 2e-10);
                    }
                }
                let work: f64 = if constrained {
                    fixed
                        .iter()
                        .zip(&step.solution.increments[0].reactions_n)
                        .map(|((dof, _), reaction)| {
                            near(*reaction, row_mass[dof / 3] * acceleration[dof % 3], 2e-9);
                            reaction * (target[dof % 3] - translation[dof % 3])
                        })
                        .sum()
                } else {
                    forces
                        .iter()
                        .map(|(dof, force)| {
                            force
                                * (step.solution.displacement_m[dof / 3][dof % 3]
                                    - translation[dof % 3])
                        })
                        .sum()
                };
                near(
                    step.energy.kinetic_before_j,
                    0.5 * mass * velocity.iter().map(|v| v * v).sum::<f64>(),
                    2e-12,
                );
                near(
                    step.energy.kinetic_after_j,
                    0.5 * mass * final_velocity.iter().map(|v| v * v).sum::<f64>(),
                    2e-12,
                );
                near(
                    step.energy.numerical_loss_j,
                    0.5 * mass * dt * dt * acceleration.iter().map(|v| v * v).sum::<f64>(),
                    2e-12,
                );
                near(work, step.energy.inertial_work_j, 2e-12);
                near(
                    work,
                    step.energy.kinetic_change_j + step.energy.numerical_loss_j,
                    2e-12,
                );
            }
        }
    });
}

#[test]
fn motion_requires_full_inventory_one_physical_endpoint_and_combined_matrix_budget() {
    with_cx(|cx| {
        let law = metal();
        let states = vec![law.initial_state(); 8];
        let mut native = problem(&NODES, &law, Formulation::BBar, &[], &[]);
        let advance = |p: &SmallStrainHexProblem<'_, J2Plasticity>, density: &[f64], dt| {
            p.advance_backward_euler(&[0.0; 24], &[[0.0; 3]; 8], &states, density, dt, cx)
        };
        assert!(advance(&native, &[], DT).is_err());
        assert!(advance(&native, &[0.0], DT).is_err());
        assert!(advance(&native, &[2.0], 0.0).is_err());
        native.load_path = &[0.5, 1.0];
        assert!(advance(&native, &[2.0], DT).is_err());
        native.load_path = &[1.0];
        native.budget.maximum_tangent_entries = 767;
        assert!(advance(&native, &[2.0], DT).is_err());
        native.budget.maximum_tangent_entries = 768;
        assert!(advance(&native, &[2.0], DT).is_ok());
        let mut invalid_profile = [[2.0; 8]];
        invalid_profile[0][7] = -1.0;
        assert!(
            native
                .advance_backward_euler_with_corner_density(
                    &[0.0; 24],
                    &[[0.0; 3]; 8],
                    &states,
                    &invalid_profile,
                    DT,
                    cx,
                )
                .is_err()
        );
        assert_eq!(states, vec![law.initial_state(); 8]);
    });
}
