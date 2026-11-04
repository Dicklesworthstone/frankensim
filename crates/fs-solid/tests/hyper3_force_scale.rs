//! Free-equation force balance must not inherit unrelated support scales.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::hyper::{Hyperelastic, HyperelasticModel};
use fs_solid::hyper3::{
    HyperTetProblem, HyperTetSettings, HyperTetSolution, HyperelasticLawField,
};
use fs_solid::linear3::TetAssemblyBudget;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 3,
                kernel_id: 17,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn law(scale: f64) -> Hyperelastic {
    Hyperelastic::new(
        HyperelasticModel::NeoHookean {
            mu: scale,
            lambda: 2.0 * scale,
        },
        4.0,
    )
    .unwrap()
}

/// One axial free coordinate on a unit tetrahedron, optionally beside a
/// disconnected, fully prescribed body with a trillion times its modulus.
/// The added body's axial stretch is 1.2; it cannot change the first body's
/// mechanics. A rigid rotation about the loading axis keeps the scalar
/// reference equation independent of the production stress/tangent routines.
fn solve(
    companion: bool,
    force_at_support: f64,
    translation: [f64; 3],
    rotated: bool,
    axial_force: f64,
) -> HyperTetSolution {
    let mut nodes = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ];
    let mut tetrahedra = vec![[0, 1, 2, 3]];
    let mut materials = vec![law(1.0)];
    if companion {
        nodes.extend([
            [3.0, 0.0, 0.0],
            [4.0, 0.0, 0.0],
            [3.0, 1.0, 0.0],
            [3.0, 0.0, 1.0],
        ]);
        tetrahedra.push([4, 5, 6, 7]);
        materials.push(law(1e12));
    }
    let initial: Vec<f64> = nodes
        .iter()
        .flat_map(|&[_, y, z]| {
            let rotation = if rotated {
                [0.0, -z - y, y - z]
            } else {
                [0.0; 3]
            };
            core::array::from_fn::<_, 3, _>(|i| rotation[i] + translation[i])
        })
        .collect();
    let prescribed: Vec<(usize, f64)> = initial
        .iter()
        .enumerate()
        .filter(|&(dof, _)| dof != 3)
        .map(|(dof, &value)| {
            let stretch = if companion && dof == 15 { 0.2 } else { 0.0 };
            (dof, value + stretch)
        })
        .collect();
    let forces = [(3, axial_force), (0, force_at_support)];
    let problem = HyperTetProblem {
        nodes_m: &nodes,
        tetrahedra: &tetrahedra,
        material: HyperelasticLawField::PerElement(&materials),
        prescribed_m: &prescribed,
        nodal_forces_n: &forces,
        body_force_n_m3: [0.0; 3],
        budget: TetAssemblyBudget::standard(),
        settings: HyperTetSettings::default(),
    };
    with_cx(|cx| problem.solve_from(&initial, cx)).unwrap()
}

/// For F = R diag(s,1,1) and a rotation about x, native tetrahedron volume
/// 1/6 and unit reference gradient give f_x = P_xx/6. This explicit
/// Neo-Hookean equation uses neither the material API nor the FE assembly.
fn independent_force(extension: f64) -> f64 {
    let stretch = 1.0 + extension;
    (stretch - stretch.recip() + 2.0 * stretch.ln() / stretch) / 6.0
}

fn independent_extension(force: f64) -> f64 {
    let (mut low, mut high) = (0.0, 2.0);
    assert!(independent_force(high) > force);
    for _ in 0..80 {
        let middle = 0.5 * low + 0.5 * high;
        if independent_force(middle) < force {
            low = middle;
        } else {
            high = middle;
        }
    }
    0.5 * low + 0.5 * high
}

fn assert_equilibrium(solution: &HyperTetSolution, force: f64) -> f64 {
    let extension = solution.displacement_m[1][0] - solution.displacement_m[0][0];
    let relative_force_error = (independent_force(extension) - force).abs() / force;
    assert!(
        relative_force_error < 2e-10,
        "reported equilibrium has relative physical force error {relative_force_error:e}; extension {extension:e}"
    );
    assert!((extension - independent_extension(force)).abs() < 2e-10);
    let relative = *solution.steps.last().unwrap().residual_history.last().unwrap();
    assert!(relative.is_finite() && relative <= HyperTetSettings::default().relative_tolerance);
    extension
}

#[test]
fn disconnected_prescribed_stiffness_cannot_hide_compliant_force_imbalance() {
    let force = 0.2;
    let isolated = solve(false, 0.0, [0.0; 3], false, force);
    let beside_stiff_body = solve(true, 0.0, [0.0; 3], false, force);
    let expected = assert_equilibrium(&isolated, force);
    let actual = assert_equilibrium(&beside_stiff_body, force);
    assert!((actual - expected).abs() < 2e-10);
    // The old all-DOF normalization accepted its first linearized update
    // (extension 0.3), despite a roughly 22% free-force imbalance.
    assert!((actual - 0.3).abs() > 0.1);
}

#[test]
fn load_applied_only_to_a_prescribed_dof_changes_only_its_reaction() {
    let force = 0.2;
    let isolated = solve(false, 0.0, [0.0; 3], false, force);
    let support_loaded = solve(false, 1e12, [0.0; 3], false, force);
    let expected = assert_equilibrium(&isolated, force);
    let actual = assert_equilibrium(&support_loaded, force);
    assert!((actual - expected).abs() < 2e-10);
    let reaction = |solution: &HyperTetSolution| {
        solution.reactions_n.iter().find(|&&(dof, _)| dof == 0).unwrap().1
    };
    assert!((reaction(&support_loaded) - reaction(&isolated) + 1e12).abs() < 1e-3);
}

#[test]
fn connected_force_scale_preserves_translation_and_rigid_rotation() {
    let force = 0.2;
    let baseline = solve(true, 0.0, [0.0; 3], false, force);
    let expected = assert_equilibrium(&baseline, force);
    for rotated in [false, true] {
        let transformed = solve(true, 0.0, [8.0, -4.0, 2.0], rotated, force);
        let actual = assert_equilibrium(&transformed, force);
        assert!((actual - expected).abs() < 2e-10);
    }
}

#[test]
fn zero_force_reference_converges_without_an_artificial_force_floor() {
    let solution = solve(false, 0.0, [0.0; 3], false, 0.0);
    assert!(solution.displacement_m.iter().flatten().all(|&value| value == 0.0));
    assert_eq!(solution.steps[0].residual_history, vec![0.0]);
    assert_eq!(solution.steps[0].linear_iterations, 0);
}
