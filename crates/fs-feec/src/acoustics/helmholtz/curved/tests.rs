use super::*;
use crate::acoustics::helmholtz::heterogeneous::HeterogeneousHelmholtzFe;
use crate::acoustics::helmholtz::heterogeneous::time::exact::ExactTimeLimits;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_math::det;
use fs_sparse::Csr;

const MEDIUM: Medium = Medium {
    density_kg_m3: 1.0,
    sound_speed_m_s: 1.0,
};
const CORNERS: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
];

fn run(f: impl FnOnce(&CancelGate, &Cx<'_>)) {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(
            &gate,
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 0x43555256,
                    kernel_id: 1,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ),
        );
    });
}

fn map(x: [f64; 3], a: f64) -> [f64; 3] {
    [x[0], x[1], x[2] + a * x[0] * x[1]]
}

fn midpoints(mesh: &TetComplex, a: f64) -> Vec<[f64; 3]> {
    mesh.edges
        .iter()
        .map(|edge| {
            map(
                std::array::from_fn(|axis| {
                    f64::midpoint(
                        CORNERS[edge[0] as usize][axis],
                        CORNERS[edge[1] as usize][axis],
                    )
                }),
                a,
            )
        })
        .collect()
}

// Independent nodal interpolation into the original hierarchical P2 basis.
// Its edge basis is lambda_i lambda_j, whose midpoint value is one quarter.
fn physical_linear(fe: &HelmholtzFe<'_>, mids: &[[f64; 3]], axis: usize) -> Vec<f64> {
    let mut p = vec![0.0; fe.ndof()];
    for i in 0..4 {
        p[i] = CORNERS[i][axis];
    }
    for (i, edge) in fe.space.complex.edges.iter().enumerate() {
        p[fe.space.edge_off + i] =
            4.0 * mids[i][axis] - 2.0 * (p[edge[0] as usize] + p[edge[1] as usize]);
    }
    p
}

fn form(matrix: &Csr, left: &[f64], right: &[f64]) -> f64 {
    let mut result = vec![0.0; right.len()];
    matrix.spmv(right, &mut result);
    left.iter().zip(result).map(|(a, b)| a * b).sum()
}

fn relative_matrices(a: &Csr, b: &Csr) -> f64 {
    let mut error = 0.0;
    let mut scale = 0.0;
    for j in 0..a.ncols() {
        let mut unit = vec![0.0; a.ncols()];
        unit[j] = 1.0;
        let mut av = vec![0.0; a.nrows()];
        let mut bv = av.clone();
        a.spmv(&unit, &mut av);
        b.spmv(&unit, &mut bv);
        for (x, y) in av.iter().zip(bv) {
            error += (x - y) * (x - y);
            scale += y * y;
        }
    }
    det::sqrt(error / scale)
}

#[test]
fn affine_quadratic_geometry_recovers_all_original_pressure_orders() {
    run(|_, cx| {
        let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
        let mids = midpoints(&mesh, 0.0);
        for order in 1..=4 {
            let original = HelmholtzFe::new(&mesh, &CORNERS, order, MEDIUM).unwrap();
            let curved = HelmholtzFe::new_quadratic(
                &mesh,
                &CORNERS,
                &mids,
                order,
                order + 2,
                MEDIUM,
                QuadraticGeometryControls::default(),
                cx,
            )
            .unwrap();
            assert_eq!(curved.ndof(), original.ndof());
            assert!(relative_matrices(curved.stiffness(), original.stiffness()) < 3e-12);
            assert!(relative_matrices(curved.mass(), original.mass()) < 3e-12);
            let faces = original.boundary_faces();
            assert!(
                relative_matrices(
                    &curved.boundary_mass(&faces).unwrap(),
                    &original.boundary_mass(&faces).unwrap()
                ) < 3e-12
            );
            assert!(
                relative_matrices(
                    &curved.boundary_tangential_stiffness(&faces).unwrap(),
                    &original.boundary_tangential_stiffness(&faces).unwrap()
                ) < 3e-12
            );
            for face in faces {
                assert_eq!(
                    curved.boundary_geometry(face),
                    original.boundary_geometry(face)
                );
            }
        }
    });
}

#[test]
fn curved_physical_fields_use_original_volume_gradients_and_material_points() {
    run(|_, cx| {
        let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
        let a = 1.25;
        let mids = midpoints(&mesh, a);
        let fe = HelmholtzFe::new_quadratic(
            &mesh,
            &CORNERS,
            &mids,
            2,
            6,
            MEDIUM,
            QuadraticGeometryControls::default(),
            cx,
        )
        .unwrap();
        let mut one = vec![0.0; fe.ndof()];
        one[..4].fill(1.0);
        let x = physical_linear(&fe, &mids, 0);
        let y = physical_linear(&fe, &mids, 1);
        let z = physical_linear(&fe, &mids, 2);
        // X=(xi,eta,zeta+a xi eta) has det J=1. These are independent
        // simplex moments, not measurements of the production quadrature.
        let volume = 1.0 / 6.0;
        for coordinate in [&x, &y, &z] {
            assert!((form(fe.stiffness(), coordinate, coordinate) - volume).abs() < 2e-12);
        }
        assert!(form(fe.stiffness(), &x, &z).abs() < 2e-12);
        assert!(form(fe.stiffness(), &y, &z).abs() < 2e-12);
        let z_squared_integral = 1.0 / 60.0 + a / 360.0 + a * a / 1260.0;
        assert!((form(fe.mass(), &z, &z) - z_squared_integral).abs() < 2e-13);
        assert!((form(fe.mass(), &one, &one) - volume).abs() < 2e-13);

        // rho^-1=1+X and beta=1+Z test material queries in PHYSICAL space.
        let heterogeneous = HeterogeneousHelmholtzFe::new(
            &fe,
            6,
            &mut |_, p| {
                Ok::<_, AcousticFeError>(Medium {
                    density_kg_m3: 1.0 / (1.0 + p[0]),
                    sound_speed_m_s: det::sqrt((1.0 + p[0]) / (1.0 + p[2])),
                })
            },
            &mut || false,
        )
        .unwrap();
        assert!((form(heterogeneous.stiffness(), &z, &z) - 5.0 / 24.0).abs() < 2e-12);
        let weighted_volume = 1.0 / 6.0 + 1.0 / 24.0 + a / 120.0;
        assert!((form(heterogeneous.mass(), &one, &one) - weighted_volume).abs() < 2e-12);

        // Liveness: retaining the same pressure coefficients but flattening
        // geometry gives the wrong physical kinetic energy by a^2/30.
        let flat = HelmholtzFe::new(&mesh, &CORNERS, 2, MEDIUM).unwrap();
        let flattened = form(flat.stiffness(), &z, &z);
        assert!((flattened - volume - a * a / 30.0).abs() < 2e-12);
        assert!((flattened - volume).abs() > 0.04);
    });
}

// Independent composite Simpson integration of the graph z=a*x*y. Both
// integrands follow directly from the graph metric and the tangential
// projection of the PHYSICAL unit x gradient, without FE geometry helpers.
fn graph_integrals(a: f64) -> [f64; 2] {
    let n = 512;
    let mut result = [0.0; 2];
    let weight = |i| {
        if i == 0 || i == n {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        }
    };
    for i in 0..=n {
        let x = i as f64 / n as f64;
        for j in 0..=n {
            let y = (1.0 - x) * (j as f64 / n as f64);
            let jac = det::sqrt(1.0 + a * a * (x * x + y * y));
            let w = weight(i) * weight(j) * (1.0 - x) / (9.0 * (n * n) as f64);
            result[0] += w * jac;
            result[1] += w * (1.0 + a * a * x * x) / jac;
        }
    }
    result
}

#[test]
fn curved_wall_area_and_tangential_viscous_form_have_live_quadrature_ladders() {
    run(|_, cx| {
        let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
        let mids = midpoints(&mesh, 2.0);
        let reference = graph_integrals(2.0);
        let face = [0, 1, 2];
        let mut previous = [f64::INFINITY; 2];
        for q in [4, 6, 8] {
            let fe = HelmholtzFe::new_quadratic(
                &mesh,
                &CORNERS,
                &mids,
                2,
                q,
                MEDIUM,
                QuadraticGeometryControls::default(),
                cx,
            )
            .unwrap();
            let x = physical_linear(&fe, &mids, 0);
            let mut one = vec![0.0; fe.ndof()];
            one[..4].fill(1.0);
            let surface = fe.boundary_mass(&[face]).unwrap();
            let tangent = fe.boundary_tangential_stiffness(&[face]).unwrap();
            let (load, area) = fe.boundary_load(&[face]).unwrap();
            assert!((area - form(&surface, &one, &one)).abs() < 2e-12);
            assert!((area - one.iter().zip(load).map(|(a, b)| a * b).sum::<f64>()).abs() < 2e-12);
            let energy = form(&tangent, &x, &x);
            let errors = [
                (area / reference[0] - 1.0).abs(),
                (energy / reference[1] - 1.0).abs(),
            ];
            for i in 0..2 {
                assert!(errors[i] < previous[i] * 0.1);
            }
            previous = errors;
            // The opposite vertex has identically zero trace; no normal
            // derivative may leak into its viscous surface energy.
            let mut normal = vec![0.0; fe.ndof()];
            normal[3] = 1.0;
            assert_eq!(form(&tangent, &normal, &normal), 0.0);
            assert!(fe.boundary_geometry(face).is_err());
            assert!(fe.boundary_geometry([0, 1, 3]).is_ok());
        }
        assert!(previous.iter().all(|e| *e < 5e-8));
    });
}

#[test]
fn curved_original_field_reaches_exact_transient_evolution_and_budget_refusals() {
    run(|gate, cx| {
        let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
        let mids = midpoints(&mesh, 1.0);
        let fe = HelmholtzFe::new_quadratic(
            &mesh,
            &CORNERS,
            &mids,
            2,
            5,
            MEDIUM,
            QuadraticGeometryControls::default(),
            cx,
        )
        .unwrap();
        let z = physical_linear(&fe, &mids, 2);
        let model = HeterogeneousHelmholtzFe::new(
            &fe,
            5,
            &mut |_, _| Ok::<_, AcousticFeError>(MEDIUM),
            &mut || false,
        )
        .unwrap();
        let prepared = model
            .prepare_exact_time(0.02, &[], &[], ExactTimeLimits::default(), &mut || false)
            .unwrap();
        let mut state = prepared.state(&vec![0.0; model.ndof()], &z).unwrap();
        let expected = (1.0 / 60.0 + 1.0 / 360.0 + 1.0 / 1260.0) / 2.0;
        assert!((prepared.energy(&state).unwrap().total_j() / expected - 1.0).abs() < 3e-12);
        for _ in 0..16 {
            let next = prepared.step(&state, &[], &mut || false).unwrap();
            assert!(next.work.balance_residual_j().abs() < expected * 2e-10);
            state = next.state;
        }
        assert!((prepared.energy(&state).unwrap().total_j() / expected - 1.0).abs() < 2e-9);
        for q in [3, 17] {
            assert!(
                HelmholtzFe::new_quadratic(
                    &mesh,
                    &CORNERS,
                    &mids,
                    2,
                    q,
                    MEDIUM,
                    QuadraticGeometryControls::default(),
                    cx
                )
                .is_err()
            );
        }
        assert!(
            HelmholtzFe::new_quadratic(
                &mesh,
                &CORNERS,
                &mids[..5],
                2,
                4,
                MEDIUM,
                QuadraticGeometryControls::default(),
                cx
            )
            .is_err()
        );
        gate.request();
        assert!(matches!(
            HelmholtzFe::new_quadratic(
                &mesh,
                &CORNERS,
                &mids,
                2,
                4,
                MEDIUM,
                QuadraticGeometryControls::default(),
                cx
            ),
            Err(AcousticFeError::Cancelled)
        ));
    });
}
