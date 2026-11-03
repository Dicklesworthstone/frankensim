//! G0/G1/G3: exact graded reference inventory in native finite-strain inertia.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::hyper::{Hyperelastic, HyperelasticModel};
use fs_solid::{HyperTetError, HyperTetProblem, HyperTetSettings, HyperelasticLawField};
use fs_sparse::Csr;

const NODES: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];
// rho(X) = 3 + 12 X, with mass 1 kg on the unit tetrahedron.
const DENSITY: [f64; 4] = [3.0, 15.0, 3.0, 3.0];

fn with_cx(run: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        run(
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 31,
                    kernel_id: 17,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ),
            &gate,
        );
    });
}

fn material() -> Hyperelastic {
    Hyperelastic::new(
        HyperelasticModel::NeoHookean {
            mu: 2.0,
            lambda: 3.0,
        },
        5.0,
    )
    .unwrap()
}

fn problem<'a>(
    nodes: &'a [[f64; 3]],
    tets: &'a [[usize; 4]],
    material: &'a Hyperelastic,
) -> HyperTetProblem<'a> {
    HyperTetProblem {
        nodes_m: nodes,
        tetrahedra: tets,
        material: HyperelasticLawField::Uniform(material),
        prescribed_m: &[],
        nodal_forces_n: &[],
        body_force_n_m3: [0.0; 3],
        budget: Default::default(),
        settings: HyperTetSettings::default(),
    }
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-13 * expected.abs().max(1.0),
        "{actual:.17e} != {expected:.17e}"
    );
}

fn kinetic(mass: &Csr, velocity: &[f64]) -> f64 {
    let mut momentum = vec![0.0; velocity.len()];
    mass.spmv(velocity, &mut momentum);
    0.5 * velocity
        .iter()
        .zip(momentum)
        .map(|(v, p)| v * p)
        .sum::<f64>()
}

#[test]
fn exact_mass_centroid_and_cubic_moments_retain_supported_coordinates() {
    with_cx(|cx, _| {
        let law = material();
        let mut p = problem(&NODES, &TETS, &law);
        p.prescribed_m = &[(0, 0.0), (1, 0.0), (2, 0.0)];
        let mass = p
            .reference_mass_matrix_with_nodal_density(&DENSITY, cx)
            .unwrap();
        assert_eq!(mass.nrows(), 12);
        // Independent barycentric monomial integrals: product(alpha_i!)/6!
        // for a cubic monomial on this unit tetrahedron.
        let factorial = [1.0, 1.0, 2.0, 6.0];
        for a in 0..4 {
            for b in 0..4 {
                let exact: f64 = DENSITY
                    .iter()
                    .enumerate()
                    .map(|(i, rho)| {
                        let mut powers = [0; 4];
                        powers[i] += 1;
                        powers[a] += 1;
                        powers[b] += 1;
                        rho * powers
                            .into_iter()
                            .map(|power| factorial[power])
                            .product::<f64>()
                            / 720.0
                    })
                    .sum();
                for component in 0..3 {
                    close(mass.get(3 * a + component, 3 * b + component), exact);
                    close(
                        mass.get(3 * a + component, 3 * b + (component + 1) % 3),
                        0.0,
                    );
                }
            }
        }
        let nodal_mass: Vec<_> = (0..4)
            .map(|node| mass.row(3 * node).1.iter().sum::<f64>())
            .collect();
        close(nodal_mass.iter().sum(), 1.0);
        for (actual, expected) in nodal_mass.iter().zip([0.225, 0.325, 0.225, 0.225]) {
            close(*actual, expected);
        }
        for axis in 0..3 {
            close(
                (0..4)
                    .map(|node| nodal_mass[node] * NODES[node][axis])
                    .sum(),
                [0.325, 0.225, 0.225][axis],
            );
        }
        // Equal total mass does not imply the same mass distribution.
        let averaged = p.reference_mass_matrix(&[6.0], cx).unwrap();
        close(averaged.row(3).1.iter().sum(), 0.25);
        assert!((mass.get(3, 3) - averaged.get(3, 3)).abs() > 0.04);
    });
}

#[test]
fn positive_inertia_matches_independent_physical_kinetic_energy() {
    with_cx(|cx, _| {
        let law = material();
        let p = problem(&NODES, &TETS, &law);
        let mass = p
            .reference_mass_matrix_with_nodal_density(&DENSITY, cx)
            .unwrap();
        // Positive Cholesky pivots check the whole coupled shape mass, beyond
        // positive diagonal entries or one chosen velocity.
        let mut lower = [[0.0; 12]; 12];
        for i in 0..12 {
            for j in 0..=i {
                assert_eq!(mass.get(i, j).to_bits(), mass.get(j, i).to_bits());
                let remainder =
                    mass.get(i, j) - (0..j).map(|k| lower[i][k] * lower[j][k]).sum::<f64>();
                lower[i][j] = if i == j {
                    assert!(remainder > 0.0);
                    remainder.sqrt()
                } else {
                    remainder / lower[j][j]
                };
            }
        }
        let velocity = |x: [f64; 3]| {
            [
                1.0 + 2.0 * x[0] - 0.5 * x[1],
                -0.7 + 1.2 * x[2] + 0.3 * x[0],
                0.2 + x[1] - x[2],
            ]
        };
        let nodal_velocity: Vec<_> = NODES.into_iter().flat_map(velocity).collect();
        // Independent tensor Gauss integration after the Duffy map. Three
        // points integrate the transformed degree-five polynomial exactly.
        let offset = (0.6_f64).sqrt() / 2.0;
        let gauss = [
            (0.5 - offset, 5.0 / 18.0),
            (0.5, 4.0 / 9.0),
            (0.5 + offset, 5.0 / 18.0),
        ];
        let mut exact = 0.0;
        for (u, wu) in gauss {
            for (v, wv) in gauss {
                for (w, ww) in gauss {
                    let x = [u, (1.0 - u) * v, (1.0 - u) * (1.0 - v) * w];
                    let jacobian = (1.0 - u).powi(2) * (1.0 - v);
                    let speed2 = velocity(x).iter().map(|value| value * value).sum::<f64>();
                    exact += 0.5 * (3.0 + 12.0 * x[0]) * speed2 * jacobian * wu * wv * ww;
                }
            }
        }
        close(kinetic(&mass, &nodal_velocity), exact);
        assert!(exact > 0.0);
    });
}

#[test]
fn uniform_bits_and_shared_vertex_scatter_survive_element_permutation() {
    with_cx(|cx, _| {
        let law = material();
        let nodes = [NODES[0], NODES[1], NODES[2], NODES[3], [0.0, 0.0, -1.0]];
        let tets = [[0, 1, 2, 3], [0, 2, 1, 4]];
        let mut p = problem(&nodes, &tets, &law);
        let old = p.reference_mass_matrix(&[7.3; 2], cx).unwrap();
        let uniform = p
            .reference_mass_matrix_with_nodal_density(&[7.3; 5], cx)
            .unwrap();
        for (a, b) in old.to_dense().iter().zip(uniform.to_dense()) {
            assert_eq!(a.to_bits(), b.to_bits());
        }
        let rho = [3.0, 15.0, 3.0, 3.0, 9.0];
        let original = p
            .reference_mass_matrix_with_nodal_density(&rho, cx)
            .unwrap();
        let permuted = [[2, 1, 0, 4], [1, 2, 0, 3]];
        p.tetrahedra = &permuted;
        let reordered = p
            .reference_mass_matrix_with_nodal_density(&rho, cx)
            .unwrap();
        for (a, b) in original.to_dense().iter().zip(reordered.to_dense()) {
            close(*a, b);
        }
        close(
            (0..5)
                .map(|node| original.row(3 * node).1.iter().sum::<f64>())
                .sum(),
            2.25,
        );
        // The two tetrahedra contribute opposite z-centroid moments, with
        // greater mass below z=0; shared vertices do not duplicate inventory.
        close(
            (0..5)
                .map(|node| nodes[node][2] * original.row(3 * node).1.iter().sum::<f64>())
                .sum(),
            -0.1,
        );
    });
}

#[test]
fn invalid_inventory_geometry_budget_and_cancellation_return_no_matrix() {
    with_cx(|cx, gate| {
        let law = material();
        let mut p = problem(&NODES, &TETS, &law);
        for density in [
            vec![1.0; 3],
            vec![1.0, 2.0, 3.0, 0.0],
            vec![1.0, 2.0, 3.0, -1.0],
            vec![1.0, 2.0, 3.0, f64::NAN],
            vec![1.0, 2.0, 3.0, f64::INFINITY],
            vec![f64::from_bits(1); 4],
        ] {
            assert!(matches!(
                p.reference_mass_matrix_with_nodal_density(&density, cx),
                Err(HyperTetError::InvalidInput { .. })
            ));
        }
        let orphan = [NODES[0], NODES[1], NODES[2], NODES[3], [1.0; 3]];
        p.nodes_m = &orphan;
        assert!(matches!(
            p.reference_mass_matrix_with_nodal_density(&[1.0; 5], cx),
            Err(HyperTetError::InvalidInput { .. })
        ));
        let large = NODES.map(|point| point.map(|x| 12.0 * x));
        p.nodes_m = &large;
        assert!(matches!(
            p.reference_mass_matrix_with_nodal_density(&[f64::MAX; 4], cx),
            Err(HyperTetError::InvalidInput { .. })
        ));
        p.nodes_m = &NODES;
        p.budget.maximum_nodes = 3;
        assert!(matches!(
            p.reference_mass_matrix_with_nodal_density(&DENSITY, cx),
            Err(HyperTetError::InvalidInput { .. })
        ));
        gate.request();
        assert!(matches!(
            p.reference_mass_matrix_with_nodal_density(&DENSITY, cx),
            Err(HyperTetError::Cancelled)
        ));
    });
}
