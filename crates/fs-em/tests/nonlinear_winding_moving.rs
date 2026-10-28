//! Independent joint P1/circuit and work oracles for a loaded stationary coil
//! surrounding a moving conducting channel. Synthetic constitutive sources
//! verify the numerical coupling, not a measured engineering material.
#[path = "nonlinear_winding_circuit/common.rs"]
mod common;

use common::material;
use fs_em::{
    planar::{
        eddy::{
            moving::MovingEddyControls,
            nonlinear::{PlanarNonlinearEddy, moving::NonlinearMovingEddyControls},
        },
        nonlinear::PlanarNonlinearMagnetostatics,
    },
    winding::{
        nonlinear::{
            NonlinearWindingCircuit, NonlinearWindingControls, NonlinearWindingError,
            NonlinearWindingInitial, NonlinearWindings,
        },
        ohmic::WindingConductor,
    },
};
use fs_matdb::MaterialCard;
use fs_material::conductor::ElectricalPropertyBasis;
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_rep_mesh::{Metric2, TriComplex2};

const DEPTH: f64 = 0.3;
const SIGMA: f64 = 2.0;
const END_H: f64 = 0.08;
const LOAD_OHM: f64 = 0.4;

fn mesh() -> TriComplex2 {
    let vertices = (0..6)
        .flat_map(|y| (0..4).map(move |x| [f64::from(x), f64::from(y)]))
        .collect();
    let mut faces = Vec::new();
    for y in 0..5 {
        for x in 0..3 {
            let a = 4 * y + x;
            faces.extend([[a, a + 1, a + 5], [a, a + 5, a + 4]]);
        }
    }
    TriComplex2::from_indexed_triangles(
        "stationary-strands-around-moving-channel",
        vertices,
        faces,
        Metric2::planar(DEPTH).unwrap(),
    )
    .unwrap()
}
fn profile(cell: usize) -> f64 {
    match cell / 6 {
        0 => 1.0 / 3.0,
        4 => -1.0 / 3.0,
        _ => 0.0,
    }
}
fn boundary(mesh: &TriComplex2, background_b: f64) -> Vec<(usize, f64)> {
    mesh.vertices()
        .iter()
        .enumerate()
        .filter(|(_, p)| p[0] == 0.0 || p[0] == 3.0 || p[1] == 0.0 || p[1] == 5.0)
        .map(|(i, p)| (i, background_b * p[0]))
        .collect()
}
fn velocity(mesh: &TriComplex2) -> Vec<[f64; 2]> {
    // P1 shear is divergence free and preserves the fixed Eulerian material
    // bands. Unequal interior speeds give the top/bottom pickup nonzero flux.
    mesh.vertices()
        .iter()
        .map(|p| {
            [
                if p[1] == 2.0 {
                    0.55
                } else if p[1] == 3.0 {
                    0.15
                } else {
                    0.0
                },
                0.0,
            ]
        })
        .collect()
}
fn controls() -> NonlinearWindingControls<NonlinearMovingEddyControls> {
    let base = common::controls();
    NonlinearWindingControls {
        field: NonlinearMovingEddyControls {
            moving: MovingEddyControls {
                relative_tolerance: 1e-12,
                absolute_residual_a_m: 1e-13,
                restart: 8,
                max_cycles: 50,
                max_krylov_scalars: 10000,
                max_cell_magnetic_peclet: 2.0,
                energy_tolerance_j: 1e-9,
                current_balance_tolerance_a: 1e-9,
            },
            absolute_residual_tolerance_a_m: 1e-11,
            relative_residual_tolerance: 1e-11,
            max_newton_iterations: 30,
            max_backtracks: 30,
        },
        maximum_field_trials: base.maximum_field_trials,
        maximum_newton_iterations: base.maximum_newton_iterations,
        difference_step_a: base.difference_step_a,
        terminal_voltage_tolerance_v: base.terminal_voltage_tolerance_v,
        power_tolerance_w: base.power_tolerance_w,
        energy_tolerance_j: base.energy_tolerance_j,
        newton: base.newton,
    }
}
fn circuit<'m>(
    problem: &'m PlanarNonlinearMagnetostatics,
    source: &MaterialCard,
    background_b: f64,
    initial_current: f64,
) -> NonlinearWindingCircuit<'m> {
    let mesh = problem.mesh();
    let boundary = boundary(mesh, background_b);
    let profiles: Vec<_> = (0..mesh.faces().len()).map(profile).collect();
    let j: Vec<_> = profiles.iter().map(|p| p * initial_current).collect();
    let initial = problem
        .solve(
            &j,
            &boundary,
            &mesh
                .vertices()
                .iter()
                .map(|p| background_b * p[0])
                .collect::<Vec<_>>(),
            material::controls().nonlinear,
            &mut || false,
        )
        .unwrap();
    let bulk = material::conductors(
        source,
        mesh.faces().len(),
        300.0,
        ElectricalPropertyBasis::Conductivity,
    )
    .into_iter()
    .filter(|cell| profile(cell.triangle) == 0.0)
    .collect();
    let insulating: Vec<_> = (0..mesh.faces().len())
        .filter(|&cell| profile(cell) != 0.0)
        .collect();
    let field = PlanarNonlinearEddy::new_with_impressed_current(
        &initial,
        bulk,
        &insulating,
        &boundary,
        &j,
        0.0,
        material::controls(),
        &mut || false,
    )
    .unwrap();
    let conductor = material::conductor(source, 300.0, ElectricalPropertyBasis::Conductivity);
    let strands = insulating
        .into_iter()
        .map(|triangle| WindingConductor {
            triangle,
            conductor: conductor.clone(),
        })
        .collect();
    NonlinearWindingCircuit::new(
        field,
        CircuitGraph {
            node_count: 3,
            branches: vec![
                (1, 0, Branch::VoltageSource { port: 0 }),
                (1, 2, Branch::Inductor { henries: 123.0 }),
                (2, 0, Branch::Resistor { ohms: LOAD_OHM }),
            ],
            transformers: vec![],
        },
        NonlinearWindings {
            branch_indices: vec![1],
            profiles_per_m2: vec![profiles],
            end_inductance_h: vec![END_H],
        },
        strands,
        NonlinearWindingInitial {
            currents_a: vec![initial_current],
            remaining_fluxes_wb: vec![],
            charges_c: vec![],
            external_inputs: vec![0.0],
        },
        common::config(),
        &mut || false,
    )
    .unwrap()
}
fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 3e-9 * (1.0 + expected.abs()),
        "{actual:e} != {expected:e}"
    );
}

fn geometry(mesh: &TriComplex2, face: &[u32; 3]) -> ([usize; 3], f64, [[f64; 2]; 3]) {
    let ids = face.map(|v| v as usize);
    let p = ids.map(|v| mesh.vertices()[v]);
    let det = (p[1][0] - p[0][0]) * (p[2][1] - p[0][1]) - (p[2][0] - p[0][0]) * (p[1][1] - p[0][1]);
    let grad = [
        [(p[1][1] - p[2][1]) / det, (p[2][0] - p[1][0]) / det],
        [(p[2][1] - p[0][1]) / det, (p[0][0] - p[2][0]) / det],
        [(p[0][1] - p[1][1]) / det, (p[1][0] - p[0][0]) / det],
    ];
    (ids, 0.5 * det.abs() * DEPTH, grad)
}

fn transfers(
    mesh: &TriComplex2,
    velocity: &[[f64; 2]],
    before: &[f64],
    after: &[f64],
    dt: f64,
) -> (f64, f64) {
    let mut cells = Vec::new();
    let mut integral = 0.0;
    let mut total_volume = 0.0;
    for (cell, face) in mesh.faces().iter().enumerate() {
        if profile(cell) != 0.0 {
            continue;
        }
        let (ids, volume, grad) = geometry(mesh, face);
        let gradient: [f64; 2] =
            std::array::from_fn(|axis| (0..3).map(|a| after[ids[a]] * grad[a][axis]).sum());
        let advected = ids.map(|v| velocity[v][0] * gradient[0] + velocity[v][1] * gradient[1]);
        let q: [f64; 3] =
            std::array::from_fn(|a| (after[ids[a]] - before[ids[a]]) / dt + advected[a]);
        integral += volume * q.iter().sum::<f64>() / 3.0;
        total_volume += volume;
        cells.push((volume, q, advected));
    }
    let mean = integral / total_volume;
    let mut heat = 0.0;
    let mut work = 0.0;
    for (volume, q, advected) in cells {
        let current = q.map(|q| -SIGMA * (q - mean));
        let sum_j = current.iter().sum::<f64>();
        heat += dt * volume / (12.0 * SIGMA)
            * (sum_j * sum_j + current.iter().map(|j| j * j).sum::<f64>());
        // (J_z e_z cross B).v = J_z * v.grad(A_z), with the same J_z.
        work += dt * volume / 12.0
            * (sum_j * advected.iter().sum::<f64>()
                + current
                    .iter()
                    .zip(advected)
                    .map(|(j, advected)| j * advected)
                    .sum::<f64>());
    }
    (heat, work)
}

/// Elementary barycentric gradients and consistent P1 integrals, independent
/// of FEEC, conductor assembly, magnetic recovery and the production solvers.
struct Oracle {
    k: Vec<Vec<f64>>,
    c: Vec<Vec<f64>>,
    d: Vec<Vec<f64>>,
    s: Vec<f64>,
    winding_ohm: f64,
}
impl Oracle {
    fn new(mesh: &TriComplex2, velocity: &[[f64; 2]]) -> Self {
        let n = mesh.vertices().len();
        let mut k = vec![vec![0.0; n]; n];
        let mut c = k.clone();
        let mut d = k.clone();
        let mut s = vec![0.0; n];
        let mut mean = vec![0.0; n];
        let mut advected_mean = vec![0.0; n];
        let mut total_sigma_volume = 0.0;
        let mut winding_ohm = 0.0;
        for (cell, face) in mesh.faces().iter().enumerate() {
            let (ids, volume, grad) = geometry(mesh, face);
            let density = profile(cell);
            winding_ohm += volume * density * density / SIGMA;
            for a in 0..3 {
                s[ids[a]] += volume * density / 3.0;
                for b in 0..3 {
                    k[ids[a]][ids[b]] +=
                        volume * (grad[a][0] * grad[b][0] + grad[a][1] * grad[b][1]);
                }
            }
            if density != 0.0 {
                continue;
            }
            total_sigma_volume += SIGMA * volume;
            for a in 0..3 {
                mean[ids[a]] += SIGMA * volume / 3.0;
                for b in 0..3 {
                    c[ids[a]][ids[b]] +=
                        SIGMA * volume * if a == b { 1.0 / 6.0 } else { 1.0 / 12.0 };
                    for l in 0..3 {
                        let v_dot_grad =
                            velocity[ids[l]][0] * grad[b][0] + velocity[ids[l]][1] * grad[b][1];
                        d[ids[a]][ids[b]] += SIGMA
                            * volume
                            * (if a == l { 1.0 / 6.0 } else { 1.0 / 12.0 })
                            * v_dot_grad;
                    }
                }
                for l in 0..3 {
                    advected_mean[ids[a]] += SIGMA * volume / 3.0
                        * (velocity[ids[l]][0] * grad[a][0] + velocity[ids[l]][1] * grad[a][1]);
                }
            }
        }
        // One connected floating bulk conductor: remove its conductivity-
        // weighted mean current. The offset acts on BOTH Adot and v.grad A.
        for a in 0..n {
            for b in 0..n {
                c[a][b] -= mean[a] * mean[b] / total_sigma_volume;
                d[a][b] -= mean[a] * advected_mean[b] / total_sigma_volume;
            }
        }
        Self {
            k,
            c,
            d,
            s,
            winding_ohm,
        }
    }
    fn step(
        &self,
        initial: &NonlinearWindingCircuit<'_>,
        dt: f64,
        voltage: f64,
    ) -> (Vec<f64>, f64) {
        let mesh = initial.field().problem().mesh();
        let a0 = initial.field().potential_wb_m();
        let fixed: Vec<_> = mesh
            .vertices()
            .iter()
            .map(|p| p[0] == 0.0 || p[0] == 3.0 || p[1] == 0.0 || p[1] == 5.0)
            .collect();
        let free: Vec<_> = (0..a0.len()).filter(|&v| !fixed[v]).collect();
        let n = free.len();
        let mut matrix = vec![vec![0.0; n + 1]; n + 1];
        let mut rhs = vec![0.0; n + 1];
        for (r, &a) in free.iter().enumerate() {
            for b in 0..a0.len() {
                let value = self.k[a][b] + self.c[a][b] / dt + self.d[a][b];
                rhs[r] += self.c[a][b] * a0[b] / dt;
                if fixed[b] {
                    rhs[r] -= value * a0[b];
                }
            }
            for (col, &b) in free.iter().enumerate() {
                matrix[r][col] = self.k[a][b] + self.c[a][b] / dt + self.d[a][b];
            }
            matrix[r][n] = -self.s[a];
            matrix[n][r] = self.s[a] / dt;
        }
        matrix[n][n] = LOAD_OHM + self.winding_ohm + END_H / dt;
        rhs[n] = voltage + initial.flux_linkages_wb()[0] / dt;
        for b in 0..a0.len() {
            if fixed[b] {
                rhs[n] -= self.s[b] * a0[b] / dt;
            }
        }
        let result = solve(matrix, rhs);
        let mut potential = a0.to_vec();
        for (&v, &value) in free.iter().zip(&result) {
            potential[v] = value;
        }
        (potential, result[n])
    }
}
fn solve(mut matrix: Vec<Vec<f64>>, mut rhs: Vec<f64>) -> Vec<f64> {
    let n = rhs.len();
    for a in 0..n {
        let pivot = (a..n)
            .max_by(|&i, &j| matrix[i][a].abs().total_cmp(&matrix[j][a].abs()))
            .unwrap();
        matrix.swap(a, pivot);
        rhs.swap(a, pivot);
        assert!(matrix[a][a].abs() > 1e-12);
        for b in a + 1..n {
            let factor = matrix[b][a] / matrix[a][a];
            for c in a..n {
                matrix[b][c] -= factor * matrix[a][c];
            }
            rhs[b] -= factor * rhs[a];
        }
    }
    let mut result = vec![0.0; n];
    for a in (0..n).rev() {
        result[a] =
            (rhs[a] - (a + 1..n).map(|b| matrix[a][b] * result[b]).sum::<f64>()) / matrix[a][a];
    }
    result
}

#[test]
fn moving_channel_drives_real_load_and_matches_joint_field_circuit_oracle() {
    let source = common::source(false, SIGMA);
    let problem = material::problem(&source, mesh());
    let initial = circuit(&problem, &source, 0.65, 0.0);
    let velocity = velocity(problem.mesh());
    let oracle = Oracle::new(problem.mesh(), &velocity);
    let dt = 0.2;
    let (expected_a, expected_i) = oracle.step(&initial, dt, 0.0);
    let result = initial
        .step_with_eulerian_velocity_candidate(dt, &velocity, &[0.0], controls(), &mut || false)
        .unwrap();
    for (&actual, &expected) in result
        .circuit()
        .field()
        .potential_wb_m()
        .iter()
        .zip(&expected_a)
    {
        near(actual, expected);
    }
    let step = result.step();
    let (expected_heat, expected_work) = transfers(
        problem.mesh(),
        &velocity,
        initial.field().potential_wb_m(),
        &expected_a,
        dt,
    );
    near(step.currents_a[0], expected_i);
    near(step.field.heat_j, expected_heat);
    near(step.mechanical_work_j, expected_work);
    near(initial.resistance_ohm()[0], oracle.winding_ohm);
    assert!(
        expected_i.abs() > 1e-5,
        "actual pickup must be excited: {expected_i}"
    );
    near(step.terminal_voltage_v[0], -LOAD_OHM * expected_i);
    assert!(
        step.graph_resistors
            .iter()
            .map(|cell| cell.heat_j)
            .sum::<f64>()
            > 0.0
    );
    near(
        step.heat_j - step.field.heat_j,
        dt * (LOAD_OHM + oracle.winding_ohm) * expected_i * expected_i,
    );
    assert!(step.mechanical_work_j < 0.0);
    assert!(step.field.heat_j > 0.0);
    near(step.supplied_work_j, 0.0);
    near(
        step.stored_after_j - step.stored_before_j + step.heat_j + step.numerical_dissipation_j,
        -step.mechanical_work_j,
    );
    near(step.energy_residual_j, 0.0);
    near(step.field.region_currents_a[0].1, 0.0);
    near(initial.time_s(), 0.0);
    near(initial.currents_a()[0], 0.0);
}

#[test]
fn zero_velocity_preserves_stationary_terminal_solution_and_account() {
    let source = common::source(false, SIGMA);
    let problem = material::problem(&source, mesh());
    let initial = circuit(&problem, &source, 0.65, 0.1);
    let stationary = initial
        .step_to_candidate(
            0.2,
            initial.field().conductors().to_vec(),
            initial.conductors().to_vec(),
            &[0.3],
            common::controls(),
            &mut || false,
        )
        .unwrap();
    let moving = initial
        .step_with_eulerian_velocity_candidate(
            0.2,
            &vec![[0.0; 2]; problem.mesh().vertices().len()],
            &[0.3],
            controls(),
            &mut || false,
        )
        .unwrap();
    near(moving.step().currents_a[0], stationary.step().currents_a[0]);
    near(moving.step().heat_j, stationary.step().heat_j);
    near(
        moving.step().stored_after_j,
        stationary.step().stored_after_j,
    );
    near(
        moving.step().supplied_work_j,
        stationary.step().supplied_work_j,
    );
    near(
        moving.step().numerical_dissipation_j,
        stationary.step().numerical_dissipation_j,
    );
    near(moving.step().mechanical_work_j, 0.0);
    near(stationary.step().mechanical_work_j, 0.0);
}

#[test]
fn saturated_motional_feedback_retains_physical_heat_and_separate_mechanical_work() {
    let source = common::source(true, SIGMA);
    let problem = material::problem(&source, mesh());
    let initial = circuit(&problem, &source, 1.2, 0.0);
    let result = initial
        .step_with_eulerian_velocity_candidate(
            0.15,
            &velocity(problem.mesh()),
            &[0.4],
            controls(),
            &mut || false,
        )
        .unwrap();
    let step = result.step();
    assert!(step.field.report.newton_iterations > 0);
    assert!(step.heat_j > 0.0 && step.field.heat_j > 0.0);
    assert!(step.numerical_dissipation_j >= 0.0);
    near(step.mechanical_work_j, step.field.mechanical_work_j);
    near(
        step.stored_after_j - step.stored_before_j
            + step.heat_j
            + step.numerical_dissipation_j
            + step.mechanical_work_j,
        step.supplied_work_j,
    );
    near(step.energy_residual_j, 0.0);
    assert!(step.terminal_residual_v[0].abs() < controls().terminal_voltage_tolerance_v);
    // The candidate remains the original state type, so stopping the motion
    // continues the accepted field/circuit history without a reset or remap.
    let (accepted, _) = result.into_parts();
    let next = accepted
        .step_to_candidate(
            0.3,
            accepted.field().conductors().to_vec(),
            accepted.conductors().to_vec(),
            &[0.0],
            common::controls(),
            &mut || false,
        )
        .unwrap();
    near(
        next.step().stored_before_j,
        accepted.field().energy_j() + 0.5 * END_H * accepted.currents_a()[0].powi(2),
    );
    near(next.circuit().time_s(), 0.3);
}

#[test]
fn unsupported_winding_motion_budgets_and_cancellation_leave_history_untouched() {
    let source = common::source(false, SIGMA);
    let problem = material::problem(&source, mesh());
    let initial = circuit(&problem, &source, 0.65, 0.0);
    let v = velocity(problem.mesh());
    let mut moving_strand = v.clone();
    moving_strand[0] = [0.1, 0.0];
    assert!(matches!(
        initial.step_with_eulerian_velocity_candidate(
            0.2,
            &moving_strand,
            &[0.0],
            controls(),
            &mut || false
        ),
        Err(NonlinearWindingError::Invalid(
            "fixed stranded winding cells require zero material velocity"
        ))
    ));
    assert!(
        initial
            .step_with_eulerian_velocity_candidate(
                0.2,
                &v[..v.len() - 1],
                &[0.0],
                controls(),
                &mut || false
            )
            .is_err()
    );
    let mut limits = controls();
    limits.maximum_field_trials = 1;
    assert!(matches!(
        initial.step_with_eulerian_velocity_candidate(0.2, &v, &[0.0], limits, &mut || false),
        Err(NonlinearWindingError::TrialBudget)
    ));
    let mut checks = 0;
    assert!(matches!(
        initial.step_with_eulerian_velocity_candidate(0.2, &v, &[0.0], controls(), &mut || {
            checks += 1;
            checks > 300
        }),
        Err(NonlinearWindingError::Cancelled)
    ));
    near(initial.time_s(), 0.0);
    near(initial.currents_a()[0], 0.0);
    near(initial.flux_linkages_wb()[0], 0.0);
}
