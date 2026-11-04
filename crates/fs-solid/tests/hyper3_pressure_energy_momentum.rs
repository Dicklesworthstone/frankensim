//! Native prepared forces: endpoint pressure work, spatial symmetry, tangents
//! and degenerate free-plane limits. These are authored analytical fixtures.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::hyper::{Hyperelastic, HyperelasticModel};
use fs_solid::hyper3::pressure::{ClosedTetPressureSurface, TetPressure};
use fs_solid::hyper3::{HyperTetError, HyperTetProblem, HyperTetSettings, HyperelasticLawField};
use fs_solid::linear3::TetAssemblyBudget;

const NODES: [[f64; 3]; 4] = [[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];

fn with_cx(f: impl FnOnce(&CancelGate, &Cx<'_>)) {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(
            &gate,
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 19,
                    kernel_id: 58,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ),
        );
    });
}
fn law() -> Hyperelastic {
    Hyperelastic::new(
        HyperelasticModel::NeoHookean {
            mu: 2.0,
            lambda: 3.0,
        },
        5.0,
    )
    .unwrap()
}
fn problem(law: &Hyperelastic) -> HyperTetProblem<'_> {
    HyperTetProblem {
        nodes_m: &NODES,
        tetrahedra: &TETS,
        material: HyperelasticLawField::Uniform(law),
        prescribed_m: &[],
        nodal_forces_n: &[],
        body_force_n_m3: [0.0; 3],
        budget: TetAssemblyBudget::standard(),
        settings: HyperTetSettings::default(),
    }
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}; tolerance {tolerance:e}"
    );
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    (0..3).map(|i| a[i] * b[i]).sum()
}
fn state(f: [f64; 9], translation: [f64; 3]) -> [f64; 12] {
    core::array::from_fn(|d| {
        let (a, i) = (d / 3, d % 3);
        translation[i] + (0..3).map(|j| f[3 * i + j] * NODES[a][j]).sum::<f64>() - NODES[a][i]
    })
}
fn endpoints() -> ([f64; 12], [f64; 12]) {
    (
        state(
            [1.08, 0.06, -0.04, 0.02, 0.94, 0.05, -0.03, 0.04, 1.11],
            [0.1, -0.2, 0.6],
        ),
        state(
            [0.91, 0.1, 0.02, -0.04, 1.06, -0.08, 0.02, -0.06, 0.98],
            [-0.05, 0.04, -0.4],
        ),
    )
}
fn field(clipped: bool) -> TetPressure<'static> {
    let (origin_m, pressure_pa, gradient_pa_m) = ([0.0; 3], 0.12, [0.23, -0.31, -0.89]);
    if clipped {
        TetPressure::FreeSurface {
            origin_m,
            pressure_pa,
            gradient_pa_m,
        }
    } else {
        TetPressure::Affine {
            origin_m,
            pressure_pa,
            gradient_pa_m,
        }
    }
}
// Native Neo-Hookean parameters inserted in their independent closed form;
// the pressure oracle is a truncated-power divided difference (separated nodes).
fn potential(u: &[f64; 12], clipped: bool) -> f64 {
    let x: [[f64; 3]; 4] =
        core::array::from_fn(|a| core::array::from_fn(|i| NODES[a][i] + u[3 * a + i]));
    let f: [f64; 9] = core::array::from_fn(|ij| x[ij % 3 + 1][ij / 3] - x[0][ij / 3]);
    let det = f[0] * (f[4] * f[8] - f[5] * f[7]) - f[1] * (f[3] * f[8] - f[5] * f[6])
        + f[2] * (f[3] * f[7] - f[4] * f[6]);
    let log = det.ln();
    let stored = (f.iter().map(|a| a * a).sum::<f64>() - 3.0 - 2.0 * log + 1.5 * log * log) / 6.0;
    let p = x.map(|v| 0.12 + dot(v, [0.23, -0.31, -0.89]));
    let mean = if clipped {
        (0..4)
            .map(|i| {
                p[i].max(0.0).powi(4)
                    / (4.0
                        * (0..4)
                            .filter(|&j| i != j)
                            .map(|j| p[i] - p[j])
                            .product::<f64>())
            })
            .sum()
    } else {
        p.iter().sum::<f64>() / 4.0
    };
    stored + det * mean / 6.0
}

#[test]
fn spatial_discrete_gradient_closes_native_endpoint_work_and_admissible_momenta() {
    with_cx(|_, cx| {
        let law = law();
        let p = problem(&law);
        let prepared = p.prepare(cx).unwrap();
        let surface =
            ClosedTetPressureSurface::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx)
                .unwrap();
        let (before, after) = endpoints();
        for clipped in [false, true] {
            let loading = field(clipped);
            let force = prepared
                .energy_momentum_force_with_pressure(&before, &after, &surface, loading, cx)
                .unwrap();
            let reverse = prepared
                .energy_momentum_force_with_pressure(&after, &before, &surface, loading, cx)
                .unwrap();
            let work = force
                .iter()
                .enumerate()
                .map(|(i, f)| f * (after[i] - before[i]))
                .sum::<f64>();
            near(
                work,
                potential(&after, clipped) - potential(&before, clipped),
                2e-13,
            );
            let native_work = surface.work_j(&before, &after, loading, cx).unwrap();
            near(
                work,
                p.stored_energy_j(&after, cx).unwrap()
                    - p.stored_energy_j(&before, cx).unwrap()
                    - native_work,
                2e-13,
            );
            let gradient = [0.23, -0.31, -0.89];
            let resultant = core::array::from_fn(|i| (0..4).map(|a| force[3 * a + i]).sum());
            for residual in cross(resultant, gradient) {
                near(residual, 0.0, 3e-15);
            }
            let mut torque = [0.0; 3];
            for a in 0..4 {
                let position = core::array::from_fn(|i| {
                    NODES[a][i] + 0.5 * before[3 * a + i] + 0.5 * after[3 * a + i]
                });
                let moment = cross(position, core::array::from_fn(|i| force[3 * a + i]));
                for i in 0..3 {
                    torque[i] += moment[i];
                }
            }
            near(dot(torque, gradient), 0.0, 3e-15);
            for i in 0..12 {
                near(force[i], reverse[i], 2e-14);
            }
            assert!(
                dot(resultant, gradient) > 0.0,
                "nonzero hydrostatic force must remain"
            );
        }
    });
}

#[test]
fn native_spatial_endpoint_tangent_and_coalescent_surface_limit_agree() {
    with_cx(|_, cx| {
        let law = law();
        let p = problem(&law);
        let prepared = p.prepare(cx).unwrap();
        let surface =
            ClosedTetPressureSurface::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx)
                .unwrap();
        let (before, after) = endpoints();
        for clipped in [false, true] {
            let loading = field(clipped);
            let tangent = prepared
                .energy_momentum_tangent_with_pressure(&before, &after, &surface, loading, cx)
                .unwrap();
            for j in 0..12 {
                let (mut plus, mut minus) = (after, after);
                let h = 1e-6;
                plus[j] += h;
                minus[j] -= h;
                let fp = prepared
                    .energy_momentum_force_with_pressure(&before, &plus, &surface, loading, cx)
                    .unwrap();
                let fm = prepared
                    .energy_momentum_force_with_pressure(&before, &minus, &surface, loading, cx)
                    .unwrap();
                for i in 0..12 {
                    near(tangent.get(i, j), (fp[i] - fm[i]) / (2.0 * h), 3e-9);
                }
            }
            let coalescent = prepared
                .energy_momentum_force_with_pressure(&after, &after, &surface, loading, cx)
                .unwrap();
            let original = p.internal_force(&after, cx).unwrap();
            let surface_force = surface.force_n(&after, loading, cx).unwrap();
            for i in 0..12 {
                near(coalescent[i], original[i] - surface_force[i], 3e-14);
            }
        }
    });
}

#[test]
fn coplanar_edge_arrival_keeps_native_pressure_work_and_a_finite_tangent() {
    with_cx(|_, cx| {
        let law = law();
        let p = problem(&law);
        let prepared = p.prepare(cx).unwrap();
        let surface =
            ClosedTetPressureSurface::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx)
                .unwrap();
        let loading = TetPressure::FreeSurface {
            origin_m: [0.0; 3],
            pressure_pa: 1.0,
            gradient_pa_m: [-1.0, -2.0, -2.0],
        };
        let before = [0.0; 12];
        let mut after = [0.0; 12];
        after[..3].copy_from_slice(&[0.25, 0.125, 0.25]); // node zero lands exactly on p=0, node one already lies there.
        let force = prepared
            .energy_momentum_force_with_pressure(&before, &after, &surface, loading, cx)
            .unwrap();
        let tangent = prepared
            .energy_momentum_tangent_with_pressure(&before, &after, &surface, loading, cx)
            .unwrap();
        let work = (0..12).map(|i| force[i] * after[i]).sum::<f64>();
        near(
            work,
            p.stored_energy_j(&after, cx).unwrap()
                - p.stored_energy_j(&before, cx).unwrap()
                - surface.work_j(&before, &after, loading, cx).unwrap(),
            2e-13,
        );
        for row in 0..12 {
            assert!(tangent.row(row).1.iter().all(|x| x.is_finite()));
        }
        // Across this two-vertex plane event the complete potential remains
        // twice differentiable; the native endpoint tangent retains that limit.
        for j in 0..12 {
            let (mut plus, mut minus) = (after, after);
            let h = 1e-6;
            plus[j] += h;
            minus[j] -= h;
            let fp = prepared
                .energy_momentum_force_with_pressure(&before, &plus, &surface, loading, cx)
                .unwrap();
            let fm = prepared
                .energy_momentum_force_with_pressure(&before, &minus, &surface, loading, cx)
                .unwrap();
            for i in 0..12 {
                near(tangent.get(i, j), (fp[i] - fm[i]) / (2.0 * h), 8e-8);
            }
        }
    });
}

#[test]
fn uniform_bits_and_conservative_field_admission_remain_native() {
    with_cx(|gate, cx| {
        let law = law();
        let p = problem(&law);
        let prepared = p.prepare(cx).unwrap();
        let surface =
            ClosedTetPressureSurface::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx)
                .unwrap();
        let (before, after) = endpoints();
        for pressure in [0.0, 0.7, -0.3] {
            let old = prepared
                .energy_momentum_force_with_uniform_pressure(
                    &before, &after, &surface, pressure, cx,
                )
                .unwrap();
            let old_tangent = prepared
                .energy_momentum_tangent_with_uniform_pressure(
                    &before, &after, &surface, pressure, cx,
                )
                .unwrap();
            for loading in [
                TetPressure::Uniform(pressure),
                TetPressure::Affine {
                    origin_m: [0.0; 3],
                    pressure_pa: pressure,
                    gradient_pa_m: [0.0; 3],
                },
            ] {
                let new = prepared
                    .energy_momentum_force_with_pressure(&before, &after, &surface, loading, cx)
                    .unwrap();
                assert_eq!(
                    old.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
                    new.iter().map(|v| v.to_bits()).collect::<Vec<_>>()
                );
                let tangent = prepared
                    .energy_momentum_tangent_with_pressure(&before, &after, &surface, loading, cx)
                    .unwrap();
                for i in 0..12 {
                    for j in 0..12 {
                        assert_eq!(old_tangent.get(i, j).to_bits(), tangent.get(i, j).to_bits());
                    }
                }
            }
        }
        assert!(
            prepared
                .energy_momentum_force_with_pressure(
                    &before,
                    &after,
                    &surface,
                    TetPressure::PerFace(&[0.0; 4]),
                    cx
                )
                .is_err()
        );
        gate.request();
        assert!(matches!(
            prepared.energy_momentum_force_with_pressure(
                &before,
                &after,
                &surface,
                field(true),
                cx
            ),
            Err(HyperTetError::Cancelled)
        ));
    });
}
