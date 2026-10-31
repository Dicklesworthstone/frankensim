//! Independent mechanics checks for the shared-Newton hexahedral plastic owner.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::plastic::{J2Plasticity, J2State};
use fs_material::{IsotropicElastic, SmallStrainLaw};
use fs_solid::plastic3::hex::{HexAssemblyBudget, SmallStrainHexProblem, SmallStrainHexSolution};
use fs_solid::{Formulation, IncrementSettings, SmallStrainLawField, SmallStrainTetError};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 7,
                kernel_id: 31,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

const E: f64 = 1000.0;
const NU: f64 = 0.3;
const SY: f64 = 1.0;
const HI: f64 = 100.0;
const HK: f64 = 50.0;

fn metal() -> J2Plasticity {
    J2Plasticity::new_mixed(IsotropicElastic::new(E, NU, 0.05).unwrap(), SY, HI, HK).unwrap()
}

fn mesh(
    n: [usize; 3],
    map: impl Fn(f64, f64, f64) -> [f64; 3],
) -> (Vec<[f64; 3]>, Vec<[usize; 8]>) {
    let mut nodes = Vec::new();
    let mut cells = Vec::new();
    let id = |i, j, k| (k * (n[1] + 1) + j) * (n[0] + 1) + i;
    for k in 0..=n[2] {
        for j in 0..=n[1] {
            for i in 0..=n[0] {
                nodes.push(map(
                    i as f64 / n[0] as f64,
                    j as f64 / n[1] as f64,
                    k as f64 / n[2] as f64,
                ));
            }
        }
    }
    for k in 0..n[2] {
        for j in 0..n[1] {
            for i in 0..n[0] {
                cells.push([
                    id(i, j, k),
                    id(i + 1, j, k),
                    id(i + 1, j + 1, k),
                    id(i, j + 1, k),
                    id(i, j, k + 1),
                    id(i + 1, j, k + 1),
                    id(i + 1, j + 1, k + 1),
                    id(i, j + 1, k + 1),
                ]);
            }
        }
    }
    (nodes, cells)
}

fn problem<'a, L: SmallStrainLaw>(
    nodes: &'a [[f64; 3]],
    cells: &'a [[usize; 8]],
    law: &'a L,
    prescribed: &'a [(usize, f64)],
    path: &'a [f64],
) -> SmallStrainHexProblem<'a, L> {
    SmallStrainHexProblem {
        nodes_m: nodes,
        hexahedra: cells,
        laws: SmallStrainLawField::Uniform(law),
        formulation: Formulation::BBar,
        strain_limit: 0.02,
        rotation_limit_rad: 0.02,
        prescribed_m: prescribed,
        nodal_forces_n: &[],
        load_path: path,
        budget: HexAssemblyBudget::default(),
        settings: IncrementSettings::default(),
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual {actual:.15e}, expected {expected:.15e}, tolerance {tolerance:.3e}"
    );
}

fn supports(nodes: &[[f64; 3]], end: Option<f64>) -> Vec<(usize, f64)> {
    let mut fixed = Vec::new();
    for (i, x) in nodes.iter().enumerate() {
        if x[0] == 0.0 {
            fixed.push((3 * i, 0.0));
        }
        if x[0] == 1.0
            && let Some(strain) = end
        {
            fixed.push((3 * i, strain));
        }
        if x[1] == 0.0 {
            fixed.push((3 * i + 1, 0.0));
        }
        if x[2] == 0.0 {
            fixed.push((3 * i + 2, 0.0));
        }
    }
    fixed
}

fn axial_return(strain: f64, plastic: &mut f64, alpha: &mut f64, viscosity: f64) -> f64 {
    let trial = E * (strain - *plastic) - HK * *plastic;
    let increment = (trial.abs() - SY - HI * *alpha).max(0.0) / (E + HI + HK + viscosity);
    *plastic += trial.signum() * increment;
    *alpha += increment;
    E * (strain - *plastic)
}

#[test]
fn hex8_geometry_and_cyclic_plastic_patch_preserve_transverse_flow_and_permanent_set() {
    use fs_solid::linear3::hex::{HEX27_REFERENCE_NODES, hex8_shapes};
    for (a, reference) in HEX27_REFERENCE_NODES[..8].iter().enumerate() {
        let (shape, _) = hex8_shapes(reference.map(f64::from));
        for (b, value) in shape.iter().enumerate() {
            close(*value, if a == b { 1.0 } else { 0.0 }, 1e-15);
        }
    }
    for p in [[0.0; 3], [0.31, -0.68, 0.23]] {
        let (shape, gradient) = hex8_shapes(p);
        close(shape.iter().sum(), 1.0, 1e-15);
        for i in 0..3 {
            close(gradient.iter().map(|g| g[i]).sum(), 0.0, 1e-15);
        }
    }
    let (nodes, cells) = mesh([2, 1, 1], |s, t, v| [s, t, v]);
    let law = metal();
    let fixed = supports(&nodes, Some(1.0));
    for formulation in [Formulation::Standard, Formulation::BBar] {
        let mut accepted: Option<SmallStrainHexSolution<J2State>> = None;
        let (mut plastic, mut alpha) = (0.0, 0.0);
        for strain in [0.0005, 0.004, 0.001, -0.004] {
            let path = [strain];
            let mut p = problem(&nodes, &cells, &law, &fixed, &path);
            p.formulation = formulation;
            let next = with_cx(|cx| match &accepted {
                None => p.solve(cx),
                Some(old) => p.solve_from(
                    old.displacement_m.as_flattened(),
                    &old.quadrature_states,
                    cx,
                ),
            })
            .unwrap();
            let stress = axial_return(strain, &mut plastic, &mut alpha, 0.0);
            let transverse = -NU * stress / E - 0.5 * plastic;
            assert_eq!(next.element_offsets, [0, 8, 16]);
            close(next.quadrature_volume_m3.iter().sum(), 1.0, 2e-14);
            for ((eps, s), state) in next
                .quadrature_strain
                .iter()
                .zip(&next.quadrature_stress_pa)
                .zip(&next.quadrature_states)
            {
                close(eps[0], strain, 2e-11);
                close(eps[1], transverse, 2e-11);
                close(eps[2], transverse, 2e-11);
                close(state.plastic_strain[0], plastic, 2e-11);
                close(state.alpha, alpha, 2e-11);
                close(s[0], stress, 2e-8);
                for value in &s[1..] {
                    close(*value, 0.0, 2e-8);
                }
            }
            if strain == 0.004 {
                let released = supports(&nodes, None);
                let mut unload = problem(&nodes, &cells, &law, &released, &[0.0]);
                unload.formulation = formulation;
                let residual = with_cx(|cx| {
                    unload.solve_from(
                        next.displacement_m.as_flattened(),
                        &next.quadrature_states,
                        cx,
                    )
                })
                .unwrap();
                for (x, u) in nodes.iter().zip(&residual.displacement_m) {
                    close(u[0], plastic * x[0], 2e-11);
                    close(u[1], -0.5 * plastic * x[1], 2e-11);
                    close(u[2], -0.5 * plastic * x[2], 2e-11);
                }
            }
            accepted = Some(next);
        }
    }
}

#[test]
fn distorted_bbar_hex_keeps_local_histories_and_the_full_nonproportional_tangent() {
    use fs_material::plastic::chaboche::{BackstressBranch, ChabocheJ2};
    use fs_material::plastic::voce::VoceJ2;
    let law = ChabocheJ2::new(
        VoceJ2::new(metal(), 0.0, 1.0).unwrap(),
        [BackstressBranch::new(80.0, 500.0).unwrap()],
    )
    .unwrap();
    let (nodes, cells) = mesh([1, 1, 1], |s, t, v| {
        [s + 0.18 * s * t, t + 0.1 * t * v, v + 0.12 * s * v]
    });
    let preload: Vec<_> = nodes
        .iter()
        .flat_map(|x| [0.004 * x[0], -0.001 * x[1], -0.001 * x[2]])
        .enumerate()
        .collect();
    let old = with_cx(|cx| problem(&nodes, &cells, &law, &preload, &[1.0]).solve(cx)).unwrap();
    let u: Vec<_> = nodes
        .iter()
        .flat_map(|x| {
            [
                0.003 * x[0] + 0.003 * x[1] * x[2],
                -0.001 * x[1] + 0.001 * x[0] * x[2],
                0.002 * x[0] * x[1] - 0.001 * x[2],
            ]
        })
        .collect();
    let fixed: Vec<_> = nodes
        .iter()
        .enumerate()
        .filter(|(_, x)| x[0] == 0.0)
        .flat_map(|(n, _)| {
            (0..3).map({
                let u = &u;
                move |i| (3 * n + i, u[3 * n + i])
            })
        })
        .collect();
    let mut p = problem(&nodes, &cells, &law, &fixed, &[1.0]);
    let k = with_cx(|cx| p.tangent_matrix(&u, &old.quadrature_states, cx)).unwrap();
    let (mut error, mut scale, mut skew) = (0.0_f64, 0.0_f64, 0.0_f64);
    let h = 1e-8;
    for j in 0..u.len() {
        let (mut plus, mut minus) = (u.clone(), u.clone());
        plus[j] += h;
        minus[j] -= h;
        let fp = with_cx(|cx| p.internal_force(&plus, &old.quadrature_states, cx)).unwrap();
        let fm = with_cx(|cx| p.internal_force(&minus, &old.quadrature_states, cx)).unwrap();
        for i in 0..u.len() {
            let derivative = (fp[i] - fm[i]) / (2.0 * h);
            error = error.hypot(k.get(i, j) - derivative);
            scale = scale.hypot(derivative);
            skew = skew.hypot(k.get(i, j) - k.get(j, i));
        }
    }
    assert!(error / scale < 3e-6, "full tangent error {}", error / scale);
    assert!(
        skew / scale > 1e-5,
        "fixture must exercise a nonsymmetric tangent"
    );
    let forces: Vec<_> = with_cx(|cx| p.internal_force(&u, &old.quadrature_states, cx))
        .unwrap()
        .into_iter()
        .enumerate()
        .collect();
    p.nodal_forces_n = &forces;
    let result = with_cx(|cx| {
        p.solve_from(
            old.displacement_m.as_flattened(),
            &old.quadrature_states,
            cx,
        )
    })
    .unwrap();
    for (actual, expected) in result.displacement_m.as_flattened().iter().zip(&u) {
        close(*actual, *expected, 2e-10);
    }
    let alphas: Vec<_> = result
        .quadrature_states
        .iter()
        .map(|s| s.plastic.alpha)
        .collect();
    assert!(
        alphas.iter().any(|a| (a - alphas[0]).abs() > 1e-5),
        "local histories were averaged"
    );
    let trace: Vec<_> = result
        .quadrature_strain
        .iter()
        .map(|e| e[0] + e[1] + e[2])
        .collect();
    for value in &trace {
        close(*value, trace[0], 2e-14);
    }
    assert!(
        result
            .quadrature_strain
            .iter()
            .any(|e| (e[3] - result.quadrature_strain[0][3]).abs() > 1e-5)
    );
}

#[test]
fn timed_overstress_advances_once_per_increment_and_publishes_its_endpoint() {
    use fs_material::plastic::rate::LinearOverstressJ2;
    let (nodes, cells) = mesh([1, 1, 1], |s, t, v| [s, t, v]);
    let law = LinearOverstressJ2::new(metal(), 200.0, 0.4).unwrap();
    let fixed = supports(&nodes, Some(0.005));
    let result =
        with_cx(|cx| problem(&nodes, &cells, &law, &fixed, &[1.0, 1.0, 1.0]).solve(cx)).unwrap();
    let (mut plastic, mut alpha, mut stress) = (0.0, 0.0, 0.0);
    for _ in 0..3 {
        stress = axial_return(0.005, &mut plastic, &mut alpha, 500.0);
    }
    for (state, s) in result
        .quadrature_states
        .iter()
        .zip(&result.quadrature_stress_pa)
    {
        close(state.alpha, alpha, 2e-11);
        close(s[0], stress, 2e-8);
        close(E * (0.005 - state.plastic_strain[0]), s[0], 2e-8);
    }
    assert!(result.increments.iter().all(|r| r.evolving_elements == 1));
}

#[test]
fn heterogeneous_hexes_keep_distinct_yield_surfaces() {
    let (nodes, cells) = mesh([2, 1, 1], |s, t, v| [s, t, v]);
    let laws = [
        metal(),
        J2Plasticity::new_mixed(IsotropicElastic::new(E, NU, 0.05).unwrap(), 20.0, HI, HK).unwrap(),
    ];
    let gamma = 0.008;
    let fixed: Vec<_> = nodes
        .iter()
        .flat_map(|x| [0.5 * gamma * x[1], 0.5 * gamma * x[0], 0.0])
        .enumerate()
        .collect();
    let mut p = problem(&nodes, &cells, &laws[0], &fixed, &[1.0]);
    p.laws = SmallStrainLawField::PerElement(&laws);
    let result = with_cx(|cx| p.solve(cx)).unwrap();
    let shear = E / (2.0 * (1.0 + NU));
    let q = 3.0_f64.sqrt() * shear * gamma;
    for e in 0..2 {
        let increment = (q - if e == 0 { SY } else { 20.0 }).max(0.0) / (3.0 * shear + HI + HK);
        let tau = shear * gamma * (1.0 - 3.0 * shear * increment / q);
        for point in 8 * e..8 * (e + 1) {
            close(result.quadrature_states[point].alpha, increment, 2e-12);
            close(result.quadrature_stress_pa[point][3], tau, 2e-9);
        }
    }
    assert_eq!(result.increments[0].evolving_elements, 1);
}

fn cylinder(nr: usize, nt: usize, formulation: Formulation) -> f64 {
    let height = 0.2;
    let (nodes, cells) = mesh([nr, nt, 1], |s, t, v| {
        let angle = 0.5 * core::f64::consts::PI * t;
        [(1.0 + s) * angle.cos(), (1.0 + s) * angle.sin(), height * v]
    });
    let id = |i, j, k| (k * (nt + 1) + j) * (nr + 1) + i;
    let mut fixed = Vec::new();
    for k in 0..=1 {
        for j in 0..=nt {
            for i in 0..=nr {
                let node = id(i, j, k);
                fixed.push((3 * node + 2, 0.0)); // Axially constrained, plane-strain Lamé reference.
                if j == 0 {
                    fixed.push((3 * node + 1, 0.0));
                }
                if j == nt {
                    fixed.push((3 * node, 0.0));
                }
            }
        }
    }
    let mut forces = Vec::new();
    for j in 0..nt {
        let a = id(0, j, 0);
        let b = id(0, j + 1, 0);
        let dx = nodes[b][0] - nodes[a][0];
        let dy = nodes[b][1] - nodes[a][1];
        // Exact constant pressure work on the planar inner reference face;
        // its area vector is shared equally by the four bilinear shape functions.
        for node in [a, b, id(0, j, 1), id(0, j + 1, 1)] {
            forces.push((3 * node, 0.001 * height * dy / 4.0));
            forces.push((3 * node + 1, -0.001 * height * dx / 4.0));
        }
    }
    let law = IsotropicElastic::new(E, 0.4999, 0.05).unwrap();
    let mut p = problem(&nodes, &cells, &law, &fixed, &[1.0]);
    p.formulation = formulation;
    p.nodal_forces_n = &forces;
    p.settings.linear_tolerance = 1e-9;
    let result = with_cx(|cx| p.solve(cx))
        .unwrap_or_else(|error| panic!("cylinder {nr}x{nt} {formulation:?}: {error}"));
    let mut radial = 0.0;
    for k in 0..=1 {
        for j in 0..=nt {
            let n = id(0, j, k);
            let x = nodes[n];
            let u = result.displacement_m[n];
            radial += (x[0] * u[0] + x[1] * u[1]) / x[0].hypot(x[1]);
        }
    }
    radial / (2 * (nt + 1)) as f64
}

#[test]
fn three_dimensional_bbar_cylinder_reduces_volumetric_locking_and_refines() {
    let nu = 0.4999;
    let exact = (1.0 + nu) / E * ((1.0 - 2.0 * nu) * 0.001 / 3.0 + 0.004 / 3.0);
    let standard = (cylinder(2, 8, Formulation::Standard) / exact - 1.0).abs();
    let coarse = (cylinder(2, 8, Formulation::BBar) / exact - 1.0).abs();
    let fine = (cylinder(4, 16, Formulation::BBar) / exact - 1.0).abs();
    eprintln!(
        "3D Lamé cylinder relative displacement error: standard2x8 {standard:.9e}, BBar2x8 {coarse:.9e}, BBar4x16 {fine:.9e}"
    );
    assert!(
        standard > 0.2,
        "fixture did not exercise volumetric locking"
    );
    assert!(coarse < 0.04, "coarse B-bar error {coarse}");
    assert!(
        fine < 0.006 && fine < 0.3 * coarse,
        "refinement: coarse {coarse}, fine {fine}"
    );
}

#[test]
fn exhausted_hex_budget_refuses_without_changing_plastic_retry_history() {
    let (nodes, cells) = mesh([3, 2, 1], |s, t, v| [s + 0.1 * s * t, t, v]);
    let law = metal();
    let preload: Vec<_> = nodes
        .iter()
        .flat_map(|x| [0.004 * x[0], -0.001 * x[1], -0.001 * x[2]])
        .enumerate()
        .collect();
    let loaded = with_cx(|cx| problem(&nodes, &cells, &law, &preload, &[1.0]).solve(cx)).unwrap();
    assert!(loaded.quadrature_states.iter().all(|s| s.alpha > 0.0));
    let preserved = loaded.clone();
    let fixed: Vec<_> = nodes
        .iter()
        .enumerate()
        .filter(|(_, x)| x[0] == 0.0)
        .flat_map(|(n, x)| {
            [
                (3 * n, 0.0),
                (3 * n + 1, -0.001 * x[1]),
                (3 * n + 2, -0.001 * x[2]),
            ]
        })
        .collect();
    let tip = nodes.len() - 1;
    let forces = [
        (3 * tip, 0.003),
        (3 * tip + 1, -0.002),
        (3 * tip + 2, 0.001),
    ];
    let mut p = problem(&nodes, &cells, &law, &fixed, &[1.0]);
    p.nodal_forces_n = &forces;
    p.settings.max_linear_iterations = 1;
    let error = with_cx(|cx| {
        p.solve_from(
            loaded.displacement_m.as_flattened(),
            &loaded.quadrature_states,
            cx,
        )
    })
    .unwrap_err();
    assert!(
        matches!(error,SmallStrainTetError::LinearSolveFailed { relative_residual,.. }
        if relative_residual.is_finite() && relative_residual>1e3*p.settings.linear_tolerance)
    );
    assert_eq!(loaded, preserved);
    p.settings.max_linear_iterations = IncrementSettings::default().max_linear_iterations;
    let resumed = with_cx(|cx| {
        p.solve_from(
            loaded.displacement_m.as_flattened(),
            &loaded.quadrature_states,
            cx,
        )
    })
    .unwrap();
    assert_ne!(resumed.displacement_m, loaded.displacement_m);
    assert_eq!(loaded, preserved);
}

#[test]
fn hex_admission_bounds_memory_geometry_raw_strain_and_rotation() {
    let (nodes, cells) = mesh([1, 1, 1], |s, t, v| [s, t, v]);
    let law = metal();
    let virgin = vec![law.initial_state(); 8];
    let mut p = problem(&nodes, &cells, &law, &[], &[1.0]);
    p.budget.maximum_quadrature_points = 7;
    assert!(with_cx(|cx| p.solve(cx)).is_err());
    p.budget = HexAssemblyBudget::default();
    p.budget.maximum_tangent_entries = 575;
    assert!(with_cx(|cx| p.solve(cx)).is_err());
    p.budget = HexAssemblyBudget::default();
    let translated: Vec<_> = nodes.iter().flat_map(|_| [1000.0, -20.0, 3.0]).collect();
    let force = with_cx(|cx| p.internal_force(&translated, &virgin, cx)).unwrap();
    assert!(force.iter().all(|x| *x == 0.0));
    let spin: Vec<_> = nodes
        .iter()
        .flat_map(|x| [-0.5 * x[1], 0.5 * x[0], 0.0])
        .collect();
    assert!(with_cx(|cx| p.internal_force(&spin, &virgin, cx)).is_err());
    // Raw maximum exceeds0.02 although its volumetric projection does not.
    let high = 0.5 + 0.5 / 3.0_f64.sqrt();
    let raw = 0.0215 * high * 1.5_f64.sqrt();
    let projected = 0.0215 * ((7.0 / 6.0) * high * high + 1.0 / 12.0).sqrt();
    assert!(raw > p.strain_limit && projected < p.strain_limit);
    let hidden: Vec<_> = nodes
        .iter()
        .flat_map(|x| [0.0215 * x[0] * x[1], 0.0, 0.0])
        .collect();
    assert!(with_cx(|cx| p.internal_force(&hidden, &virgin, cx)).is_err());
    let mut folded = nodes.clone();
    folded[cells[0][6]] = [0.1, 0.1, 0.1];
    p.nodes_m = &folded;
    assert!(with_cx(|cx| p.internal_force(&vec![0.0; 24], &virgin, cx)).is_err());
}
