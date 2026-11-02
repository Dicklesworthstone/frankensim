//! Positive graded reference density in the native quadratic mass operator.
//! Analytic barycentric moments are independent of the assembly quadrature.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::ContentHash;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_modal::eigh_gen_dense;
use fs_solid::linear3::tet10::{QuadraticTetProblem, TET10_EDGES, TetQuadrature};
use fs_solid::linear3::{
    QuadraticThermalStrainField, TetAssemblyBudget, TetElasticAssembly, TetElasticError,
    TetElasticMaterial, TetMaterialField, TetThermalStrainState,
};

const NODES: [[f64; 3]; 10] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
    [0.5, 0.0, 0.0],
    [0.5, 0.5, 0.0],
    [0.0, 0.5, 0.0],
    [0.0, 0.0, 0.5],
    [0.5, 0.0, 0.5],
    [0.0, 0.5, 0.5],
];
const ELEMENTS: [[usize; 10]; 1] = [[0, 1, 2, 3, 4, 5, 6, 7, 8, 9]];
const GRADED: [[f64; 4]; 1] = [[1000.0, 1000.0, 1000.0, 5000.0]];

fn with_cx(run: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 59,
                kernel_id: 59,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        run(&cx, &gate);
    });
}

fn problem<'a>(
    nodes: &'a [[f64; 3]],
    material: &'a TetElasticMaterial,
    fixed: &'a [usize],
    mass_points: usize,
) -> QuadraticTetProblem<'a> {
    QuadraticTetProblem {
        nodes_m: nodes,
        elements: &ELEMENTS,
        materials: TetMaterialField::Uniform(material),
        fixed_dofs: fixed,
        budget: TetAssemblyBudget::standard(),
        quadrature: TetQuadrature {
            stiffness_points_per_dir: 3,
            mass_points_per_dir: mass_points,
        },
    }
}

fn material() -> TetElasticMaterial {
    TetElasticMaterial::try_new(2000.0, 12e9, 0.25, [7; 32]).unwrap()
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-12 * expected.abs().max(1.0),
        "{actual:.17e} != {expected:.17e}"
    );
}

type Polynomial = Vec<([usize; 4], f64)>;

fn multiply(a: &Polynomial, b: &Polynomial) -> Polynomial {
    a.iter()
        .flat_map(|(p, c)| {
            b.iter()
                .map(move |(q, d)| (core::array::from_fn(|i| p[i] + q[i]), c * d))
        })
        .collect()
}

fn integrate(poly: &Polynomial) -> f64 {
    let factorial = |n: usize| (1..=n).map(|i| i as f64).product::<f64>();
    poly.iter()
        .map(|(p, c)| {
            c * p.iter().map(|&i| factorial(i)).product::<f64>()
                / factorial(3 + p.iter().sum::<usize>())
        })
        .sum()
}

fn shapes() -> Vec<Polynomial> {
    (0..10)
        .map(|i| {
            let mut p = [0; 4];
            if i < 4 {
                p[i] = 1;
                let mut q = p;
                q[i] = 2;
                vec![(q, 2.0), (p, -1.0)]
            } else {
                let [a, b] = TET10_EDGES[i - 4];
                p[a] = 1;
                p[b] = 1;
                vec![(p, 4.0)]
            }
        })
        .collect()
}

fn analytic_measure(curvature: [f64; 3]) -> Polynomial {
    let density: Polynomial = GRADED[0]
        .iter()
        .enumerate()
        .map(|(i, &rho)| {
            let mut p = [0; 4];
            p[i] = 1;
            (p, rho)
        })
        .collect();
    // X_i = lambda_(i+1) + a_i lambda_(i+1)^2 has diagonal J,
    // det(J) = product(1 + 2 a_i lambda_(i+1)).
    curvature.iter().enumerate().fold(density, |poly, (i, a)| {
        let mut p = [0; 4];
        p[i + 1] = 1;
        multiply(&poly, &vec![([0; 4], 1.0), (p, 2.0 * a)])
    })
}

fn check_all_mass_entries(assembly: &TetElasticAssembly, measure: &Polynomial) {
    let mass = assembly.mass.to_dense();
    let shapes = shapes();
    for i in 0..10 {
        for j in 0..10 {
            let exact = integrate(&multiply(&multiply(&shapes[i], &shapes[j]), measure));
            for c in 0..3 {
                for d in 0..3 {
                    close(
                        mass[(3 * i + c) * 30 + 3 * j + d],
                        if c == d { exact } else { 0.0 },
                    );
                }
            }
        }
    }
    close(assembly.total_mass_kg, integrate(measure));
    // Density is positive throughout, so the FULL consistent matrix is SPD,
    // despite the negative row sums at the quadratic corner basis functions.
    fs_la::factor::cholesky(&mass, 30).unwrap();
    assert!(assembly.mass.row(0).1.iter().sum::<f64>() < 0.0);
}

#[test]
fn affine_mass_retains_degree_five_moments_and_rigid_inertia() {
    with_cx(|cx, _| {
        let material = material();
        let a = problem(&NODES, &material, &[], 4)
            .assemble_with_reference_density(&GRADED, cx)
            .unwrap();
        let measure = analytic_measure([0.0; 3]);
        check_all_mass_entries(&a, &measure);
        let translation: Vec<_> = (0..30)
            .map(|i| if i % 3 == 0 { 1.0 } else { 0.0 })
            .collect();
        let mut momentum = vec![0.0; 30];
        a.mass.spmv(&translation, &mut momentum);
        for axis in 0..3 {
            let first_moment: f64 = NODES
                .iter()
                .enumerate()
                .map(|(i, p)| p[axis] * momentum[3 * i])
                .sum();
            let mut power = [0; 4];
            power[axis + 1] = 1;
            close(
                first_moment,
                integrate(&multiply(&vec![(power, 1.0)], &measure)),
            );
        }
        // Rotation about x: |omega cross X|² = y² + z². The native mass
        // therefore must retain the exact graded rotational inertia too.
        let rotation: Vec<_> = NODES.iter().flat_map(|p| [0.0, -p[2], p[1]]).collect();
        a.mass.spmv(&rotation, &mut momentum);
        let twice_energy: f64 = rotation.iter().zip(momentum).map(|(u, p)| u * p).sum();
        close(
            twice_energy,
            integrate(&multiply(
                &vec![([0, 0, 2, 0], 1.0), ([0, 0, 0, 2], 1.0)],
                &measure,
            )),
        );
    });
}

#[test]
fn curved_mass_requires_and_uses_the_degree_eight_rule() {
    with_cx(|cx, _| {
        let material = material();
        let curvature = [0.3, 0.2, 0.4];
        let nodes = NODES.map(|p| core::array::from_fn(|i| p[i] + curvature[i] * p[i] * p[i]));
        for points in [4, 5] {
            assert!(matches!(
                problem(&nodes, &material, &[], points)
                    .assemble_with_reference_density(&GRADED, cx),
                Err(TetElasticError::QuadratureInadmissible { .. })
            ));
        }
        let a = problem(&nodes, &material, &[], 6)
            .assemble_with_reference_density(&GRADED, cx)
            .unwrap();
        check_all_mass_entries(&a, &analytic_measure(curvature));
        let refined = problem(&nodes, &material, &[], 7)
            .assemble_with_reference_density(&GRADED, cx)
            .unwrap();
        for (actual, expected) in a.mass.to_dense().iter().zip(refined.mass.to_dense()) {
            close(*actual, expected);
        }
    });
}

#[test]
fn density_location_changes_supported_modes_at_identical_total_mass() {
    with_cx(|cx, _| {
        let material = material();
        // Only the (0,3) midside node moves: each free component has shape
        // 4 lambda_0 lambda_3. Uniform and graded M are scalar multiples of I,
        // so the exact frequency shift is independent of stiffness orientation.
        let fixed: Vec<_> = (0..30).filter(|i| !(21..24).contains(i)).collect();
        let p = problem(&NODES, &material, &fixed, 4);
        let uniform = p.assemble(cx).unwrap();
        let graded = p.assemble_with_reference_density(&GRADED, cx).unwrap();
        close(graded.total_mass_kg, uniform.total_mass_kg);
        assert_eq!(graded.free_dofs, [21, 22, 23]);
        assert_eq!(graded.stiffness.to_dense(), uniform.stiffness.to_dense());
        let mass = graded.mass.to_dense();
        let original = uniform.mass.to_dense();
        for i in 0..3 {
            for j in 0..3 {
                close(mass[3 * i + j], if i == j { 2000.0 / 63.0 } else { 0.0 });
                close(mass[3 * i + j], 1.25 * original[3 * i + j]);
            }
        }
        let modes = |a: &TetElasticAssembly| {
            eigh_gen_dense(&a.stiffness.to_dense(), &a.mass.to_dense(), 3).unwrap()
        };
        for (after, before) in modes(&graded).iter().zip(modes(&uniform)) {
            close(after.lambda / before.lambda, 0.8);
            assert!(after.residual <= 1e-12 * after.lambda);
        }
    });
}

#[test]
fn constant_reference_profile_preserves_existing_scalar_bits() {
    with_cx(|cx, _| {
        let material = material();
        let mut curved = NODES;
        curved[4][2] += 0.03;
        for nodes in [&NODES, &curved] {
            let p = problem(nodes, &material, &[], 4);
            let old = p.assemble(cx).unwrap();
            let new = p
                .assemble_with_reference_density(&[[material.density_kg_m3; 4]], cx)
                .unwrap();
            assert_eq!(old.mass.to_dense(), new.mass.to_dense());
            assert_eq!(old.stiffness.to_dense(), new.stiffness.to_dense());
            assert_eq!(old.total_mass_kg.to_bits(), new.total_mass_kg.to_bits());
        }
    });
}

#[test]
fn thermoelastic_tensors_change_stiffness_without_refreshing_graded_mass() {
    with_cx(|cx, _| {
        let material = material();
        let hot = TetElasticMaterial::try_new(2000.0, 6e9, 0.25, [8; 32]).unwrap();
        let state = TetThermalStrainState::try_new(
            600.0,
            300.0,
            [0.001; 6],
            [8; 32],
            ContentHash([9; 32]),
            0.01,
        )
        .unwrap();
        let states = vec![state; 27];
        let materials = vec![hot; 27];
        let field = QuadraticThermalStrainField::ThermomechanicalAtQuadrature {
            states: &states,
            materials: &materials,
            points_per_element: 27,
        };
        let p = problem(&NODES, &material, &[], 4);
        let cold = p.assemble_with_reference_density(&GRADED, cx).unwrap();
        let warmed = p
            .assemble_thermoelastic_with_reference_density(field, &GRADED, cx)
            .unwrap();
        assert_eq!(cold.mass.to_dense(), warmed.mass.to_dense());
        assert_eq!(cold.total_mass_kg.to_bits(), warmed.total_mass_kg.to_bits());
        for (k, hot_k) in cold
            .stiffness
            .to_dense()
            .iter()
            .zip(warmed.stiffness.to_dense())
        {
            close(0.5 * k, hot_k);
        }
    });
}

#[test]
fn invalid_density_budget_and_cancellation_refuse_without_a_partial_assembly() {
    with_cx(|cx, gate| {
        let material = material();
        let mut p = problem(&NODES, &material, &[], 4);
        assert!(matches!(
            p.assemble_with_reference_density(&[], cx),
            Err(TetElasticError::MaterialCountMismatch { .. })
        ));
        for rho in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::from_bits(1)] {
            let density = [[rho; 4]];
            assert!(matches!(
                p.assemble_with_reference_density(&density, cx),
                Err(TetElasticError::InvalidMaterial { .. })
            ));
        }
        let large = NODES.map(|p| p.map(|x| 12.0 * x));
        p.nodes_m = &large;
        assert!(matches!(
            p.assemble_with_reference_density(&[[f64::MAX; 4]], cx),
            Err(TetElasticError::InvalidMaterial { .. })
        ));
        p.budget.maximum_nodes = 9;
        assert!(matches!(
            p.assemble_with_reference_density(&GRADED, cx),
            Err(TetElasticError::BudgetExceeded { what: "nodes", .. })
        ));
        gate.request();
        assert!(matches!(
            p.assemble_with_reference_density(&GRADED, cx),
            Err(TetElasticError::Cancelled)
        ));
    });
}
