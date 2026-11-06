//! G0/G1: directional yielding and permanent forming through native solid solvers.
//!
//! The independent free-sided uniaxial oracle uses the Hill flow normal and
//! isotropic elastic compliance. Coefficients are synthetic, not calibrated
//! rolling/formability data; all cases stay in the small-strain domain.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::plastic::hill::{Hill48, Hill48Coefficients};
use fs_material::plastic::reduction::{ReducedState, ReductionConfig};
use fs_material::{IsotropicElastic, SmallStrainLaw};
use fs_solid::{
    Formulation, IncrementSettings, Mesh2, PlaneAssemblyBudget, PlaneKind, SmallStrainLawField,
    SmallStrainPlaneProblem, SmallStrainTetProblem, TetAssemblyBudget,
};

const E: f64 = 1000.0;
const NU: f64 = 0.3;
const YIELD: f64 = 1.0;
const HARDENING: f64 = 60.0;
const THICKNESS: f64 = 0.2;
const F: f64 = 0.35;
const G: f64 = 0.8;
const H: f64 = 0.45;

fn metal() -> Hill48 {
    Hill48::new(
        IsotropicElastic::new(E, NU, 0.05).unwrap(),
        YIELD,
        HARDENING,
        Hill48Coefficients {
            f: F,
            g: G,
            h: H,
            l: 1.5,
            m: 1.1,
            n: 2.0,
        },
    )
    .unwrap()
}

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 29,
                kernel_id: 73,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        actual.is_finite() && (actual - expected).abs() <= tolerance,
        "actual {actual:.16e}, expected {expected:.16e}, tolerance {tolerance:.3e}"
    );
}

struct AxialOracle {
    stress: f64,
    strain: [f64; 6],
    plastic: [f64; 6],
    alpha: f64,
}

fn axial_oracle(axis: usize, strain: f64) -> AxialOracle {
    // q² = sigma^T Q sigma. For one nonzero normal stress, its flow normal
    // remains fixed: n_j = sign(sigma) Q_ja / sqrt(Q_aa). Hence the exact
    // free-sided return has denominator E Q_aa + H, not the J2 denominator.
    let q = [[G + H, -H, -G], [-H, F + H, -F], [-G, -F, F + G]];
    let a = q[axis][axis].sqrt();
    let alpha = (a * E * strain.abs() - YIELD).max(0.0) / (E * a * a + HARDENING);
    let mut plastic = [0.0; 6];
    for j in 0..3 {
        plastic[j] = strain.signum() * alpha * q[j][axis] / a;
    }
    let stress = E * (strain - plastic[axis]);
    let mut total = plastic;
    for (j, value) in total.iter_mut().take(3).enumerate() {
        *value += if j == axis {
            stress / E
        } else {
            -NU * stress / E
        };
    }
    AxialOracle {
        stress,
        strain: total,
        plastic,
        alpha,
    }
}

fn plane_supports(mesh: &Mesh2, axis: usize, axial_strain: Option<f64>) -> Vec<(usize, f64)> {
    let mut supports = Vec::new();
    for (node, x) in mesh.nodes.iter().enumerate() {
        for (d, &coordinate) in x.iter().enumerate() {
            if coordinate == 0.0 {
                supports.push((2 * node + d, 0.0));
            }
        }
        if x[axis] == 1.0
            && let Some(strain) = axial_strain
        {
            supports.push((2 * node + axis, strain));
        }
    }
    supports
}

fn plane_problem<'a>(
    mesh: &'a Mesh2,
    law: &'a Hill48,
    supports: &'a [(usize, f64)],
) -> SmallStrainPlaneProblem<'a, Hill48> {
    SmallStrainPlaneProblem {
        mesh,
        laws: SmallStrainLawField::Uniform(law),
        plane: PlaneKind::Stress,
        formulation: Formulation::Standard,
        thickness_m: THICKNESS,
        reduction: ReductionConfig {
            max_iterations: 40,
            stress_tolerance_pa: 1e-10,
            strain_limit: 0.02,
        },
        rotation_limit_rad: 0.02,
        prescribed_m: supports,
        nodal_forces_n: &[],
        load_path: &[1.0],
        budget: PlaneAssemblyBudget::default(),
        settings: IncrementSettings::default(),
    }
}

#[test]
fn plane_stress_reduction_preserves_distinct_directional_yield_thresholds() {
    let law = metal();
    // Y/(E sqrt(G+H)) < 0.001 < Y/(E sqrt(F+H)): x has yielded while
    // the same y-strain is still elastic, despite identical elastic moduli.
    for mesh in [
        Mesh2::quads(1.0, 1.0, 2, 2),
        Mesh2::triangles(1.0, 1.0, 2, 2),
    ] {
        for axis in 0..2 {
            let supports = plane_supports(&mesh, axis, Some(0.001));
            let solved = with_cx(|cx| plane_problem(&mesh, &law, &supports).solve(cx)).unwrap();
            let exact = axial_oracle(axis, 0.001);
            for (state, stress) in solved
                .quadrature_states
                .iter()
                .zip(&solved.quadrature_stress_pa)
            {
                assert_eq!(state.history.alpha > 0.0, axis == 0);
                near(state.history.alpha, exact.alpha, 2e-10);
                for j in 0..6 {
                    near(state.strain[j], exact.strain[j], 2e-10);
                    near(stress[j], if j == axis { exact.stress } else { 0.0 }, 2e-8);
                }
            }
        }
    }
}

#[test]
fn p1_q1_forming_and_unloading_retain_anisotropic_width_and_thickness() {
    let law = metal();
    for mesh in [
        Mesh2::quads(1.0, 1.0, 2, 2),
        Mesh2::triangles(1.0, 1.0, 2, 2),
    ] {
        for axis in 0..2 {
            for sign in [-1.0, 1.0] {
                let supports = plane_supports(&mesh, axis, Some(sign * 0.006));
                let loaded = with_cx(|cx| plane_problem(&mesh, &law, &supports).solve(cx)).unwrap();
                let exact = axial_oracle(axis, sign * 0.006);
                assert!(exact.alpha > 0.0);
                for (state, stress) in loaded
                    .quadrature_states
                    .iter()
                    .zip(&loaded.quadrature_stress_pa)
                {
                    near(state.history.alpha, exact.alpha, 2e-10);
                    for j in 0..6 {
                        near(state.strain[j], exact.strain[j], 2e-10);
                        near(state.history.plastic_strain[j], exact.plastic[j], 2e-10);
                        near(stress[j], if j == axis { exact.stress } else { 0.0 }, 2e-8);
                    }
                }
                let axial_reaction: f64 = supports
                    .iter()
                    .zip(&loaded.increments[0].reactions_n)
                    .filter(|((d, _), _)| d % 2 == axis && mesh.nodes[d / 2][axis] == 1.0)
                    .map(|(_, reaction)| reaction)
                    .sum();
                near(axial_reaction, THICKNESS * exact.stress, 4e-9);

                // Release the driven face, retaining its original accepted history.
                // Width and thickness spring back elastically but keep distinct
                // plastic contractions (or expansions under compression).
                let released = plane_supports(&mesh, axis, None);
                let unloaded = with_cx(|cx| {
                    plane_problem(&mesh, &law, &released).solve_from(
                        loaded.displacement_m.as_flattened(),
                        &loaded.quadrature_states,
                        cx,
                    )
                })
                .unwrap();
                for (x, u) in mesh.nodes.iter().zip(&unloaded.displacement_m) {
                    for j in 0..2 {
                        near(u[j], exact.plastic[j] * x[j], 2e-10);
                    }
                }
                for (state, stress) in unloaded
                    .quadrature_states
                    .iter()
                    .zip(&unloaded.quadrature_stress_pa)
                {
                    near(state.history.alpha, exact.alpha, 2e-10);
                    for j in 0..6 {
                        near(state.strain[j], exact.plastic[j], 2e-10);
                        near(state.history.plastic_strain[j], exact.plastic[j], 2e-10);
                        near(stress[j], 0.0, 2e-8);
                    }
                    near(state.strain[..3].iter().sum(), 0.0, 5e-10);
                    let width = 1 - axis;
                    near(
                        state.strain[width] / state.strain[2],
                        if axis == 0 { H / G } else { H / F },
                        2e-7,
                    );
                    assert!(state.strain[2] * sign < 0.0);
                }
            }
        }
    }
}

#[test]
fn native_tetrahedron_uses_the_same_directional_flow_and_permanent_set() {
    let law = metal();
    let nodes = [[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let tets = [[0, 1, 2, 3]];
    // Suppress rigid motion/shear; the three diagonal stretches remain available.
    let released = [
        (0, 0.0),
        (1, 0.0),
        (2, 0.0),
        (4, 0.0),
        (5, 0.0),
        (6, 0.0),
        (8, 0.0),
        (9, 0.0),
        (10, 0.0),
    ];
    for axis in 0..3 {
        for sign in [-1.0, 1.0] {
            let strain = sign * 0.006;
            let mut supports = released.to_vec();
            supports.push((3 * (axis + 1) + axis, strain));
            let problem = SmallStrainTetProblem {
                nodes_m: &nodes,
                tetrahedra: &tets,
                laws: SmallStrainLawField::Uniform(&law),
                prescribed_m: &supports,
                nodal_forces_n: &[],
                load_path: &[1.0],
                budget: TetAssemblyBudget::standard(),
                settings: IncrementSettings::default(),
            };
            let loaded = with_cx(|cx| problem.solve(cx)).unwrap();
            let exact = axial_oracle(axis, strain);
            near(loaded.element_states[0].alpha, exact.alpha, 2e-10);
            for j in 0..6 {
                near(loaded.element_strain[0][j], exact.strain[j], 2e-10);
                near(
                    loaded.element_stress[0][j],
                    if j == axis { exact.stress } else { 0.0 },
                    2e-8,
                );
            }
            near(
                *loaded.increments[0].reactions_n.last().unwrap(),
                exact.stress / 6.0,
                4e-9,
            );
            let unloading = SmallStrainTetProblem {
                prescribed_m: &released,
                ..problem
            };
            let unloaded = with_cx(|cx| {
                unloading.solve_from(
                    loaded.displacement_m.as_flattened(),
                    &loaded.element_states,
                    cx,
                )
            })
            .unwrap();
            near(unloaded.element_states[0].alpha, exact.alpha, 2e-10);
            for j in 0..6 {
                near(unloaded.element_strain[0][j], exact.plastic[j], 2e-10);
                near(
                    unloaded.element_states[0].plastic_strain[j],
                    exact.plastic[j],
                    2e-10,
                );
                near(unloaded.element_stress[0][j], 0.0, 2e-8);
            }
            for (x, u) in nodes.iter().zip(&unloaded.displacement_m) {
                for j in 0..3 {
                    near(u[j], exact.plastic[j] * x[j], 2e-10);
                }
            }
        }
    }
}

#[test]
fn condensed_hill_tangent_differentiates_native_skew_q1_forces() {
    let law = metal();
    let mesh = Mesh2::mapped_quads(2, 1, &|s, t| [2.0 * s + 0.3 * t, t + 0.2 * s]);
    let problem = plane_problem(&mesh, &law, &[]);
    let states: Vec<_> = (0..8)
        .map(|_| ReducedState {
            strain: [0.0; 6],
            history: law.initial_state(),
        })
        .collect();
    let u: Vec<_> = mesh
        .nodes
        .iter()
        .flat_map(|x| {
            [
                0.004 * x[0] + 0.001 * x[1] + 0.001 * x[0] * x[1],
                -0.0003 * x[1] + 0.0004 * x[0],
            ]
        })
        .collect();
    let direction: Vec<_> = (0..u.len())
        .map(|i| ((i * 7 % 11) as f64 - 5.0) / 7.0)
        .collect();
    let h = 1e-7;
    let plus: Vec<_> = u.iter().zip(&direction).map(|(u, d)| u + h * d).collect();
    let minus: Vec<_> = u.iter().zip(&direction).map(|(u, d)| u - h * d).collect();
    let (matrix, plus, minus) = with_cx(|cx| {
        (
            problem.tangent_matrix(&u, &states, cx).unwrap(),
            problem.internal_force(&plus, &states, cx).unwrap(),
            problem.internal_force(&minus, &states, cx).unwrap(),
        )
    });
    let mut analytic = vec![0.0; u.len()];
    matrix.spmv(&direction, &mut analytic);
    let mut error2 = 0.0;
    let mut scale2 = 0.0;
    for ((&value, plus), minus) in analytic.iter().zip(plus).zip(minus) {
        error2 += (value - (plus - minus) / (2.0 * h)).powi(2);
        scale2 += value * value;
    }
    assert!(scale2 > 0.0);
    assert!(
        (error2 / scale2).sqrt() < 2e-6,
        "relative force tangent error: {}",
        (error2 / scale2).sqrt()
    );
}
