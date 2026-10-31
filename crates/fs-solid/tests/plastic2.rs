//! Independent mechanics oracles for the native P1/Q1 plane plastic solver.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::plastic::rate::LinearOverstressJ2;
use fs_material::plastic::reduction::{ReducedState, ReductionConfig};
use fs_material::plastic::{J2Plasticity, J2State};
use fs_material::{IsotropicElastic, SmallStrainLaw};
use fs_solid::{
    IncrementSettings, Mesh2, PlaneAssemblyBudget, PlaneKind, SmallStrainLawField,
    SmallStrainPlaneProblem, SmallStrainPlaneSolution,
};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 7,
                kernel_id: 29,
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
const YIELD: f64 = 1.0;
const HI: f64 = 100.0;
const HK: f64 = 50.0;

fn metal() -> J2Plasticity {
    J2Plasticity::new_mixed(IsotropicElastic::new(E, NU, 0.05).unwrap(), YIELD, HI, HK).unwrap()
}

fn problem<'a, L: SmallStrainLaw>(
    mesh: &'a Mesh2,
    law: &'a L,
    prescribed: &'a [(usize, f64)],
    path: &'a [f64],
) -> SmallStrainPlaneProblem<'a, L> {
    SmallStrainPlaneProblem {
        mesh,
        laws: SmallStrainLawField::Uniform(law),
        plane: PlaneKind::Stress,
        thickness_m: 0.2,
        reduction: ReductionConfig {
            max_iterations: 40,
            stress_tolerance_pa: 1e-12,
            strain_limit: 0.02,
        },
        rotation_limit_rad: 0.02,
        prescribed_m: prescribed,
        nodal_forces_n: &[],
        load_path: path,
        budget: PlaneAssemblyBudget::default(),
        settings: IncrementSettings::default(),
    }
}

fn axial_supports(mesh: &Mesh2, length: f64, end_strain: Option<f64>) -> Vec<(usize, f64)> {
    let mut fixed = Vec::new();
    for (node, p) in mesh.nodes.iter().enumerate() {
        if p[0] == 0.0 {
            fixed.push((2 * node, 0.0));
        }
        if p[0] == length
            && let Some(strain) = end_strain
        {
            fixed.push((2 * node, length * strain));
        }
        if p[1] == 0.0 {
            fixed.push((2 * node + 1, 0.0));
        }
    }
    fixed
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual {actual:.15e}, expected {expected:.15e}, tolerance {tolerance:.3e}"
    );
}

/// Independent signed uniaxial plastic increment, including the Prager shift.
fn axial_return(strain: f64, plastic: &mut f64, alpha: &mut f64, viscosity_modulus: f64) -> f64 {
    let relative_trial = E * (strain - *plastic) - HK * *plastic;
    let increment =
        (relative_trial.abs() - YIELD - HI * *alpha).max(0.0) / (E + HI + HK + viscosity_modulus);
    *plastic += relative_trial.signum() * increment;
    *alpha += increment;
    E * (strain - *plastic)
}

#[test]
fn plane_stress_p1_and_q1_follow_cyclic_yield_and_keep_permanent_deformation() {
    let law = metal();
    for mesh in [
        Mesh2::quads(2.0, 1.0, 3, 2),
        Mesh2::triangles(2.0, 1.0, 3, 2),
    ] {
        let supports = axial_supports(&mesh, 2.0, Some(1.0));
        let mut accepted: Option<SmallStrainPlaneSolution<J2State>> = None;
        let (mut plastic, mut alpha) = (0.0, 0.0);
        for strain in [0.0005, 0.004, 0.001, -0.004] {
            let path = [strain];
            let p = problem(&mesh, &law, &supports, &path);
            let next = with_cx(|cx| match &accepted {
                None => p.solve(cx),
                Some(old) => p.solve_from(
                    old.displacement_m.as_flattened(),
                    &old.quadrature_states,
                    cx,
                ),
            })
            .unwrap();
            let sigma = axial_return(strain, &mut plastic, &mut alpha, 0.0);
            let transverse = -NU * sigma / E - 0.5 * plastic;
            for (state, stress) in next
                .quadrature_states
                .iter()
                .zip(&next.quadrature_stress_pa)
            {
                close(state.strain[0], strain, 2e-11);
                close(state.strain[1], transverse, 2e-11);
                close(state.strain[2], transverse, 2e-11);
                close(state.history.plastic_strain[0], plastic, 2e-11);
                close(state.history.alpha, alpha, 2e-11);
                close(stress[0], sigma, 2e-8);
                for &s in &stress[1..] {
                    close(s, 0.0, 2e-8);
                }
            }
            let right_reaction: f64 = supports
                .iter()
                .zip(&next.increments[0].reactions_n)
                .filter(|((d, _), _)| d % 2 == 0 && mesh.nodes[d / 2][0] == 2.0)
                .map(|(_, reaction)| reaction)
                .sum();
            close(right_reaction, 0.2 * sigma, 2e-9);
            if strain == 0.004 {
                let released = axial_supports(&mesh, 2.0, None);
                let p = problem(&mesh, &law, &released, &[0.0]);
                let unloaded = with_cx(|cx| {
                    p.solve_from(
                        next.displacement_m.as_flattened(),
                        &next.quadrature_states,
                        cx,
                    )
                })
                .unwrap();
                for (node, u) in mesh.nodes.iter().zip(&unloaded.displacement_m) {
                    close(u[0], plastic * node[0], 2e-11);
                    close(u[1], -0.5 * plastic * node[1], 2e-11);
                }
                for stress in &unloaded.quadrature_stress_pa {
                    for &s in stress {
                        close(s, 0.0, 2e-8);
                    }
                }
            }
            accepted = Some(next);
        }
    }
}

#[test]
fn plane_strain_retains_confining_stress_and_out_of_plane_plastic_history() {
    let mesh = Mesh2::quads(1.0, 1.0, 2, 2);
    let law = metal();
    let strain = 0.005;
    let prescribed: Vec<_> = mesh
        .nodes
        .iter()
        .enumerate()
        .flat_map(|(n, x)| [(2 * n, strain * x[0]), (2 * n + 1, 0.0)])
        .collect();
    let mut p = problem(&mesh, &law, &prescribed, &[1.0]);
    p.plane = PlaneKind::Strain;
    let result = with_cx(|cx| p.solve(cx)).unwrap();
    let mu = E / (2.0 * (1.0 + NU));
    let bulk = E / (3.0 * (1.0 - 2.0 * NU));
    let plastic = (2.0 * mu * strain - YIELD) / (3.0 * mu + HI + HK);
    let sx = bulk * strain + (4.0 / 3.0) * mu * strain - 2.0 * mu * plastic;
    let sy = bulk * strain - (2.0 / 3.0) * mu * strain + mu * plastic;
    for (state, stress) in result
        .quadrature_states
        .iter()
        .zip(&result.quadrature_stress_pa)
    {
        close(state.strain[2], 0.0, 0.0);
        close(state.history.plastic_strain[2], -0.5 * plastic, 1e-14);
        close(stress[0], sx, 1e-12);
        close(stress[1], sy, 1e-12);
        close(stress[2], sy, 1e-12);
    }
}

#[test]
fn plastic_plane_stress_tangent_differentiates_global_force_on_skew_q1_elements() {
    let mesh = Mesh2::mapped_quads(2, 1, &|s, t| [2.0 * s + 0.3 * t, t + 0.2 * s]);
    let law = metal();
    let p = problem(&mesh, &law, &[], &[1.0]);
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
    let (matrix, f_plus, f_minus) = with_cx(|cx| {
        (
            p.tangent_matrix(&u, &states, cx).unwrap(),
            p.internal_force(&plus, &states, cx).unwrap(),
            p.internal_force(&minus, &states, cx).unwrap(),
        )
    });
    let mut derivative = vec![0.0; u.len()];
    matrix.spmv(&direction, &mut derivative);
    let mut error = 0.0_f64;
    let mut scale = 0.0_f64;
    for i in 0..u.len() {
        let difference = (f_plus[i] - f_minus[i]) / (2.0 * h);
        error = error.hypot(derivative[i] - difference);
        scale = scale.hypot(difference);
        for j in 0..u.len() {
            close(matrix.get(i, j), matrix.get(j, i), 2e-11);
        }
    }
    assert!(
        error / scale < 3e-7,
        "relative consistent-tangent error {}",
        error / scale
    );
    let prescribed: Vec<_> = u.into_iter().enumerate().collect();
    let loaded = with_cx(|cx| problem(&mesh, &law, &prescribed, &[1.0]).solve(cx)).unwrap();
    let alphas: Vec<_> = loaded
        .quadrature_states
        .iter()
        .map(|s| s.history.alpha)
        .collect();
    assert!(
        alphas[..4].iter().any(|a| (a - alphas[0]).abs() > 1e-5),
        "Q1 histories were incorrectly averaged over the element"
    );
}

#[test]
fn heterogeneous_layers_keep_their_own_yield_and_contraction() {
    let mesh = Mesh2::quads(1.0, 1.0, 1, 2);
    let laws = [
        metal(),
        J2Plasticity::new(IsotropicElastic::new(2.0 * E, NU, 0.05).unwrap(), 20.0, 0.0).unwrap(),
    ];
    let supports = axial_supports(&mesh, 1.0, Some(0.004));
    let mut p = problem(&mesh, &laws[0], &supports, &[1.0]);
    p.laws = SmallStrainLawField::PerElement(&laws);
    let result = with_cx(|cx| p.solve(cx)).unwrap();
    let plastic = (E * 0.004 - YIELD) / (E + HI + HK);
    let stress = [E * (0.004 - plastic), 2.0 * E * 0.004];
    for e in 0..2 {
        for q in result.element_offsets[e]..result.element_offsets[e + 1] {
            close(result.quadrature_stress_pa[q][0], stress[e], 2e-8);
            close(
                result.quadrature_states[q].history.alpha,
                if e == 0 { plastic } else { 0.0 },
                2e-11,
            );
        }
    }
    assert_eq!(result.increments[0].evolving_elements, 1);
}

#[test]
fn repeated_load_factors_advance_overstress_time_once_and_publish_the_same_endpoint() {
    let mesh = Mesh2::quads(1.0, 1.0, 2, 1);
    let viscosity = 200.0;
    let dt = 0.4;
    let law = LinearOverstressJ2::new(metal(), viscosity, dt).unwrap();
    let supports = axial_supports(&mesh, 1.0, Some(0.005));
    let p = problem(&mesh, &law, &supports, &[1.0, 1.0, 1.0]);
    let result = with_cx(|cx| p.solve(cx)).unwrap();
    let (mut plastic, mut alpha, mut sigma) = (0.0, 0.0, 0.0);
    for _ in 0..3 {
        sigma = axial_return(0.005, &mut plastic, &mut alpha, viscosity / dt);
    }
    for (state, stress) in result
        .quadrature_states
        .iter()
        .zip(&result.quadrature_stress_pa)
    {
        close(state.history.alpha, alpha, 2e-11);
        close(stress[0], sigma, 2e-8);
        let elastic_strain =
            std::array::from_fn(|i| state.strain[i] - state.history.plastic_strain[i]);
        let from_history = law.backbone().elastic.stress(&elastic_strain, &());
        for i in 0..6 {
            close(stress[i], from_history[i], 2e-10);
        }
    }
    assert_eq!(result.increments.len(), 3);
}

#[test]
fn continuation_matches_one_path_and_domain_refusal_preserves_retry_history() {
    let mesh = Mesh2::quads(1.0, 1.0, 2, 2);
    let law = metal();
    let supports = axial_supports(&mesh, 1.0, Some(1.0));
    let p = problem(&mesh, &law, &supports, &[0.004, 0.002, -0.003]);
    let entire = with_cx(|cx| p.solve(cx)).unwrap();
    let loaded = with_cx(|cx| problem(&mesh, &law, &supports, &[0.004]).solve(cx)).unwrap();
    let preserved = loaded.clone();
    let failed = problem(&mesh, &law, &supports, &[0.1]);
    assert!(
        with_cx(|cx| failed.solve_from(
            loaded.displacement_m.as_flattened(),
            &loaded.quadrature_states,
            cx
        ))
        .is_err()
    );
    assert_eq!(loaded, preserved);
    let resumed = with_cx(|cx| {
        problem(&mesh, &law, &supports, &[0.002, -0.003]).solve_from(
            loaded.displacement_m.as_flattened(),
            &loaded.quadrature_states,
            cx,
        )
    })
    .unwrap();
    assert_eq!(resumed.displacement_m, entire.displacement_m);
    assert_eq!(resumed.quadrature_states, entire.quadrature_states);
    assert_eq!(resumed.quadrature_stress_pa, entire.quadrature_stress_pa);
}

#[test]
fn folded_quads_and_quadrature_budget_refuse_before_evolution() {
    let mut mesh = Mesh2::quads(1.0, 1.0, 1, 1);
    let law = metal();
    let mut p = problem(&mesh, &law, &[], &[1.0]);
    p.budget.maximum_quadrature_points = 3;
    assert!(with_cx(|cx| p.solve(cx)).is_err());
    mesh.nodes[3] = [0.1, 0.1]; // Positive area at one corner, folded at another.
    assert!(with_cx(|cx| problem(&mesh, &law, &[], &[1.0]).solve(cx)).is_err());
}

#[test]
fn rigid_spin_has_an_independent_small_rotation_gate() {
    let mesh = Mesh2::quads(1.0, 1.0, 1, 1);
    let law = metal();
    for (rotation, admitted) in [(0.002, true), (0.5, false)] {
        let prescribed: Vec<_> = mesh
            .nodes
            .iter()
            .enumerate()
            .flat_map(|(n, x)| {
                [
                    (2 * n, 10.0 - rotation * x[1]),
                    (2 * n + 1, -7.0 + rotation * x[0]),
                ]
            })
            .collect();
        let result = with_cx(|cx| problem(&mesh, &law, &prescribed, &[1.0]).solve(cx));
        assert_eq!(result.is_ok(), admitted);
        if let Ok(result) = result {
            for stress in result.quadrature_stress_pa {
                for s in stress {
                    close(s, 0.0, 1e-10);
                }
            }
        }
    }
}

#[test]
fn saturating_voce_hardening_drives_the_plane_continuum() {
    use fs_material::plastic::voce::VoceJ2;
    let mesh = Mesh2::quads(1.0, 1.0, 2, 2);
    let saturation = 1.8;
    let rate = 500.0;
    let law = VoceJ2::new(metal(), saturation, rate).unwrap();
    let strain = 0.006;
    let supports = axial_supports(&mesh, 1.0, Some(strain));
    let result =
        with_cx(|cx| problem(&mesh, &law, &supports, &[0.25, 0.5, 0.75, 1.0]).solve(cx)).unwrap();
    // Solve the independent monotone axial balance by bisection, without
    // invoking either the tensor return map or its axial reduction.
    let (mut lo, mut hi) = (0.0, strain);
    for _ in 0..80 {
        let p = 0.5 * (lo + hi);
        let residual =
            E * (strain - p) - YIELD - (HI + HK) * p - saturation * (1.0 - (-rate * p).exp());
        if residual > 0.0 {
            lo = p;
        } else {
            hi = p;
        }
    }
    let plastic = 0.5 * (lo + hi);
    let sigma = E * (strain - plastic);
    for (state, stress) in result
        .quadrature_states
        .iter()
        .zip(&result.quadrature_stress_pa)
    {
        close(state.history.alpha, plastic, 2e-11);
        close(state.history.plastic_strain[0], plastic, 2e-11);
        close(stress[0], sigma, 2e-8);
        close(state.strain[2], -NU * sigma / E - 0.5 * plastic, 2e-11);
    }
}
