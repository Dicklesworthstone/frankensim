//! Exact one-mode spatial winding balances; manufactured, not measured material.
#[allow(dead_code)]
#[path = "../spatial_induction/hysteretic/common.rs"]
pub mod material;

use fs_em::{
    spatial::vector_potential::induction::{
        SpatialInductionLimits,
        hysteretic::{
            SpatialHystereticInduction, SpatialHystereticProblem,
            winding::{SpatialHystereticWindingCircuit, SpatialHystereticWindingControls},
        },
        nonlinear::winding::{
            SpatialNonlinearWindingConfig, SpatialNonlinearWindings, SpatialWindingConductor,
        },
    },
    winding::{
        end::EndWindingMode,
        nonlinear::{NonlinearWindingControls, NonlinearWindingInitial, NonlinearWindingPort},
    },
};
use fs_matdb::MaterialCard;
use fs_material::conductor::ResolvedConductorState;
use fs_phs::circuit::{Branch, CircuitGraph};
use fs_solver::NewtonKrylovConfig;
use std::sync::Arc;

pub const END_L: f64 = 0.11;
pub const GRAPH_R: f64 = 0.7;
pub const CAPACITY: f64 = 10.0;
pub const DT: f64 = 0.1;

pub fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.16e} != {b:.16e}");
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.into_iter().zip(b).map(|(x, y)| x * y).sum()
}
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|k| a[k] - b[k])
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
pub fn cube(card: &Arc<MaterialCard>) -> SpatialHystereticProblem {
    let mut positions = Vec::new();
    for k in 0..=2 {
        for j in 0..=2 {
            for i in 0..=2 {
                positions.push([f64::from(i) / 2.0, f64::from(j) / 2.0, f64::from(k) / 2.0]);
            }
        }
    }
    let index = |i, j, k| ((k * 3 + j) * 3 + i) as u32;
    let mut tets = Vec::new();
    for k in 0..2 {
        for j in 0..2 {
            for i in 0..2 {
                let (a, b, c, d) = (
                    index(i, j, k),
                    index(i + 1, j, k),
                    index(i, j + 1, k),
                    index(i + 1, j + 1, k),
                );
                let (e, f, g, h) = (
                    index(i, j, k + 1),
                    index(i + 1, j, k + 1),
                    index(i, j + 1, k + 1),
                    index(i + 1, j + 1, k + 1),
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
    material::problem(card, positions, tets)
}
pub fn profile(problem: &SpatialHystereticProblem) -> Vec<[f64; 3]> {
    // curl(N_center*(-2 e_x)) has a complete closed return, no normal exterior
    // current, and one ampere through its transverse half-section.
    problem
        .mesh()
        .tets
        .iter()
        .map(|tet| {
            if !tet.contains(&13) {
                return [0.0; 3];
            }
            let face: Vec<_> = tet.iter().copied().filter(|&v| v != 13).collect();
            let origin = problem.positions_m()[face[0] as usize];
            let normal = cross(
                sub(problem.positions_m()[face[1] as usize], origin),
                sub(problem.positions_m()[face[2] as usize], origin),
            );
            let distance = dot(normal, sub(problem.positions_m()[13], origin));
            cross(normal.map(|v| v / distance), [-2.0, 0.0, 0.0])
        })
        .collect()
}
pub fn field<'m>(
    problem: &'m SpatialHystereticProblem,
    conductor: &ResolvedConductorState,
    eddy: bool,
) -> SpatialHystereticInduction<'m> {
    // This constrained FE manufactured problem leaves the [9,22] body
    // diagonal free. Its six-cell star has volume 1/8 and B=8a. Constraining
    // the other A cochains preserves a closed source and a nonzero free row.
    let diagonal = problem.mesh().edges.binary_search(&[9, 22]).unwrap();
    let prescribed = (0..problem.mesh().edges.len())
        .filter(|&edge| edge != diagonal)
        .collect::<Vec<_>>();
    let core = (0..48)
        .map(|cell| (eddy && (26..30).contains(&cell)).then(|| conductor.clone()))
        .collect();
    SpatialHystereticInduction::new(
        problem,
        core,
        &prescribed,
        SpatialInductionLimits {
            maximum_unknowns: 2000,
            maximum_conducting_cells: 100,
            maximum_assembly_terms: 100_000,
        },
        1e-10,
        &mut || false,
    )
    .unwrap()
}
pub fn strands(
    profiles: &[Vec<[f64; 3]>],
    source: &ResolvedConductorState,
    fill: f64,
) -> Vec<SpatialWindingConductor> {
    profiles
        .iter()
        .enumerate()
        .flat_map(|(winding, profile)| {
            profile
                .iter()
                .enumerate()
                .filter(|(_, j)| **j != [0.0; 3])
                .map(move |(cell, _)| SpatialWindingConductor {
                    winding,
                    cell,
                    fill_fraction: fill,
                    conductor: source.clone(),
                })
        })
        .collect()
}
pub fn controls() -> SpatialHystereticWindingControls {
    NonlinearWindingControls {
        field: material::field_controls(),
        maximum_field_trials: 200,
        maximum_newton_iterations: 25,
        difference_step_a: 1e-5,
        terminal_voltage_tolerance_v: 2e-7,
        power_tolerance_w: 2e-7,
        energy_tolerance_j: 4e-7,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-9,
            relative_tolerance: 1e-12,
            forcing_minimum: 1e-4,
            linear_restart: 8,
            ..NewtonKrylovConfig::default()
        },
    }
}
pub fn config() -> SpatialNonlinearWindingConfig {
    SpatialNonlinearWindingConfig {
        max_windings: 4,
        max_profile_values: 10_000,
        current_balance_tolerance_turns: 1e-10,
    }
}
pub fn initial(count: usize) -> NonlinearWindingInitial {
    NonlinearWindingInitial {
        currents_a: vec![0.0; count],
        remaining_fluxes_wb: vec![],
        charges_c: vec![],
        external_inputs: vec![0.0; count],
    }
}
pub fn circuit<'m>(
    problem: &'m SpatialHystereticProblem,
    source: &ResolvedConductorState,
) -> SpatialHystereticWindingCircuit<'m> {
    let profiles = vec![profile(problem)];
    SpatialHystereticWindingCircuit::new_with_end_modes(
        field(problem, source, true),
        CircuitGraph {
            node_count: 3,
            branches: vec![
                (1, 0, Branch::VoltageSource { port: 0 }),
                (1, 2, Branch::Resistor { ohms: GRAPH_R }),
                (2, 0, Branch::Inductor { henries: f64::NAN }),
            ],
            transformers: vec![],
        },
        SpatialNonlinearWindings {
            branch_indices: vec![2],
            profiles_per_m2: profiles.clone(),
            // The diagonal and shared mode both contribute: 0.03 + 0.02*2².
            end_inductance_h: vec![0.03],
        },
        vec![NonlinearWindingPort::CurrentInput],
        vec![EndWindingMode {
            inductance_h: 0.02,
            turns: vec![2.0],
        }],
        strands(&profiles, source, 1.0),
        initial(1),
        config(),
        &mut || false,
    )
    .unwrap()
}

#[derive(Clone, Copy, Debug, Default)]
pub struct Mode {
    pub h: f64,
    pub z: f64,
    pub b: f64,
    pub i: f64,
}
pub fn at_field(old: Mode, t: f64, h: f64) -> Mode {
    let z = material::play(old.z, h);
    let b = material::MU * h + z;
    // Exact integrals: int Junit.W=1/12; int |Junit|²=8. The four
    // conducting tetrahedra have M_aa=1/15 and scalar-potential Schur
    // correction 1/24, hence condensed eddy mass sigma/40. Their nodal
    // projection is 0 at vertex9, 1/2 at22, and 1/4 at all other vertices.
    let i = 12.0 * (h + (b - old.b) / (320.0 * DT * material::resistivity(t)));
    Mode { h, z, b, i }
}
pub fn voltage(old: Mode, next: Mode, t: f64, graph_r: f64) -> f64 {
    (next.b - old.b) / (96.0 * DT)
        + END_L * (next.i - old.i) / DT
        + (8.0 * material::resistivity(t) + graph_r) * next.i
}
pub fn voltage_endpoint(old: Mode, t: f64, input: f64, graph_r: f64) -> Mode {
    let h = material::bracket(-29.0, 29.0, |h| {
        voltage(old, at_field(old, t, h), t, graph_r) - input
    });
    at_field(old, t, h)
}
pub fn core_heat(old: Mode, next: Mode, t: f64) -> f64 {
    (next.b - old.b).powi(2) / (2560.0 * DT * material::resistivity(t))
}
pub fn pinning_heat(old: Mode, next: Mode) -> f64 {
    material::PIN * (next.z - old.z).abs() / 8.0
}
pub fn stored(mode: Mode) -> f64 {
    material::stored(mode.h, mode.z) / 8.0
}
pub fn assert_field(field: &SpatialHystereticInduction<'_>, mode: Mode) {
    let edge = field
        .problem()
        .mesh()
        .edges
        .binary_search(&[9, 22])
        .unwrap();
    near(field.potential_wb()[edge], mode.b / 8.0, 3e-9);
    for (cell, memory) in field.cell_states().iter().enumerate() {
        let direction = if (24..30).contains(&cell) {
            material::DIRECTIONS[cell - 24]
        } else {
            [0.0; 3]
        };
        for axis in 0..3 {
            near(memory.field_a_m()[axis], mode.h * direction[axis], 4e-7);
            near(
                memory.flux_density_t()[axis],
                mode.b * direction[axis],
                3e-8,
            );
            near(
                memory.polarizations_t()[0][axis],
                mode.z * direction[axis],
                3e-8,
            );
        }
    }
}
