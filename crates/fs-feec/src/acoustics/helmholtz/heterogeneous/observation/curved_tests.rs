//! G1/G3: original physical positions, variable Jacobians and curved fluxes.
//! These manufactured maps have independent analytic inverses and integrals;
//! no curved-cylinder or full-instrument qualification is implied.

use super::*;
use crate::acoustics::helmholtz::heterogeneous::time::HeterogeneousTimeLimits;
use crate::acoustics::helmholtz::{HelmholtzFe, QuadraticGeometryControls};
use crate::highorder::simplex::LocalFn;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_rep_mesh::TetComplex;
use std::cell::Cell;

const MEDIUM: Medium = Medium {
    density_kg_m3: 2.0,
    sound_speed_m_s: 3.0,
};
const REFERENCE: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
];

fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig {
        chunk_bytes: 4096,
        max_chunk_bytes: 4096,
        limit_bytes: Some(1 << 20),
        ..fs_alloc::ArenaConfig::default()
    })
    .scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 93,
                kernel_id: 4,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ));
    });
}

fn mapped_nodes(
    mesh: &TetComplex,
    map: impl Fn([f64; 3]) -> [f64; 3],
) -> (Vec<[f64; 3]>, Vec<[f64; 3]>) {
    let positions = REFERENCE.iter().copied().map(&map).collect();
    let midpoints = mesh
        .edges
        .iter()
        .map(|edge| {
            map(core::array::from_fn(|axis| {
                f64::midpoint(
                    REFERENCE[edge[0] as usize][axis],
                    REFERENCE[edge[1] as usize][axis],
                )
            }))
        })
        .collect();
    (positions, midpoints)
}

fn uniform<'f, 'c>(fe: &'f HelmholtzFe<'c>) -> HeterogeneousHelmholtzFe<'f, 'c> {
    HeterogeneousHelmholtzFe::new(
        fe,
        fe.order() + 2,
        &mut |_, _| Ok::<_, AcousticFeError>(MEDIUM),
        &mut || false,
    )
    .unwrap()
}

// A physical affine field pulls back to a quadratic reference polynomial.
// Convert its physical nodal values to the existing hierarchical H1 basis;
// this uses no production derivative, inverse-map or quadrature routine.
fn physical_affine_field(
    fe: &HelmholtzFe<'_>,
    midpoints: &[[f64; 3]],
    field: impl Fn([f64; 3]) -> C64,
) -> Vec<C64> {
    let mut coefficients = vec![C64::ZERO; fe.ndof()];
    for (element, tet) in fe.tetrahedra().iter().enumerate() {
        for (dof, shape) in fe.space.element_dofs(element) {
            coefficients[dof] = match shape {
                LocalFn::Vertex(a) => field(fe.positions()[tet[a] as usize]),
                LocalFn::Edge(a, b, 0) => {
                    let mut edge = [tet[a], tet[b]];
                    edge.sort_unstable();
                    let index = fe.space.complex.edges.binary_search(&edge).unwrap();
                    field(midpoints[index]).scale(4.0)
                        - (field(fe.positions()[edge[0] as usize])
                            + field(fe.positions()[edge[1] as usize]))
                        .scale(2.0)
                }
                _ => C64::ZERO,
            };
        }
    }
    coefficients
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

#[test]
fn curved_probes_recover_physical_affine_fields_and_local_material() {
    with_cx(|cx| {
        let a = 1.6;
        let shear = |q: [f64; 3]| [q[0], q[1], q[2] + a * q[0] * q[1]];
        let linear = [C64::new(0.7, -0.2), C64::new(-0.4, 1.1), C64::new(1.3, 0.6)];
        let field = |x: [f64; 3]| {
            C64::new(1.2, -0.3)
                + linear[0].scale(x[0])
                + linear[1].scale(x[1])
                + linear[2].scale(x[2])
        };
        let lambda = [0.13, 0.21, 0.29, 0.37];
        let expected_position = shear([lambda[1], lambda[2], lambda[3]]);
        for tet in [[0, 1, 2, 3], [3, 1, 0, 2]] {
            let mesh = TetComplex::from_tets(4, vec![tet]);
            let (positions, midpoints) = mapped_nodes(&mesh, shear);
            for order in [2, 3, 4] {
                let fe = HelmholtzFe::new_quadratic(
                    &mesh,
                    &positions,
                    &midpoints,
                    order,
                    8,
                    MEDIUM,
                    QuadraticGeometryControls::default(),
                    cx,
                )
                .unwrap();
                let material = |x: [f64; 3]| Medium {
                    density_kg_m3: 2.0 + x[2],
                    sound_speed_m_s: 3.0 + x[0],
                };
                let model = HeterogeneousHelmholtzFe::new(
                    &fe,
                    order + 2,
                    &mut |_, x| Ok::<_, AcousticFeError>(material(x)),
                    &mut || false,
                )
                .unwrap();
                let coefficients = physical_affine_field(&fe, &midpoints, field);
                let calls = Cell::new(0);
                let probe = model
                    .probe(0, tet.map(|v| lambda[v as usize]), &mut |element, x| {
                        calls.set(calls.get() + 1);
                        assert_eq!(element, 0);
                        for axis in 0..3 {
                            close(x[axis], expected_position[axis], 3e-15);
                        }
                        Ok::<_, AcousticFeError>(material(x))
                    })
                    .unwrap();
                let sample = probe.harmonic(&coefficients, 0.7).unwrap();
                assert!((sample.pressure_pa - field(expected_position)).abs() < 3e-14);
                for axis in 0..3 {
                    assert!((sample.pressure_gradient_pa_m[axis] - linear[axis]).abs() < 3e-14);
                    let expected_velocity = C64::new(-linear[axis].im, linear[axis].re)
                        .scale(1.0 / (0.7 * core::f64::consts::TAU * probe.medium().density_kg_m3));
                    assert!((sample.particle_velocity_m_s[axis] - expected_velocity).abs() < 3e-14);
                }
                assert_eq!(probe.harmonic(&coefficients, 0.7).unwrap(), sample);
                assert_eq!(
                    calls.get(),
                    1,
                    "warm observations must not query the source"
                );
                let constant = physical_affine_field(&fe, &midpoints, |_| C64::ONE);
                for gradient in probe
                    .harmonic(&constant, 0.7)
                    .unwrap()
                    .pressure_gradient_pa_m
                {
                    assert!(gradient.abs() < 2e-14);
                }
            }
        }
    });
}

#[test]
fn physical_probe_location_is_bounded_and_samples_only_after_convergence() {
    with_cx(|cx| {
        let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
        let map = |q: [f64; 3]| [q[0], q[1], q[2] + 1.6 * q[0] * q[1]];
        let (positions, midpoints) = mapped_nodes(&mesh, map);
        let fe = HelmholtzFe::new_quadratic(
            &mesh,
            &positions,
            &midpoints,
            2,
            8,
            MEDIUM,
            QuadraticGeometryControls::default(),
            cx,
        )
        .unwrap();
        let model = uniform(&fe);
        let limits = ProbeLocationLimits {
            max_iterations: 64,
            max_backtracks: 32,
            absolute_tolerance_m: 1e-12,
        };
        let calls = Cell::new(0);
        let mut source = |_, _| {
            calls.set(calls.get() + 1);
            Ok::<_, AcousticFeError>(MEDIUM)
        };
        for lambda in [
            [0.13, 0.21, 0.29, 0.37],
            [0.1, 0.2, 0.7, 0.0],
            [1.0, 0.0, 0.0, 0.0],
        ] {
            let requested = map([lambda[1], lambda[2], lambda[3]]);
            let before = calls.get();
            let probe = model
                .probe_at_position(0, requested, limits, &mut source, &mut || false)
                .unwrap();
            assert_eq!(calls.get(), before + 1);
            for axis in 0..3 {
                close(
                    probe.position_m()[axis],
                    requested[axis],
                    limits.absolute_tolerance_m,
                );
            }
            for (actual, expected) in probe.reference_barycentric().into_iter().zip(lambda) {
                close(actual, expected, 4e-12);
                assert!((0.0..=1.0).contains(&actual));
            }
        }
        let before = calls.get();
        for (point, bad_limits) in [
            ([-0.1, 0.2, 0.3], limits),
            ([f64::NAN, 0.0, 0.0], limits),
            (
                map([0.21, 0.29, 0.37]),
                ProbeLocationLimits {
                    max_iterations: 1,
                    ..limits
                },
            ),
            (
                [0.0; 3],
                ProbeLocationLimits {
                    max_iterations: 0,
                    ..limits
                },
            ),
            (
                [0.0; 3],
                ProbeLocationLimits {
                    max_backtracks: 65,
                    ..limits
                },
            ),
            (
                [0.0; 3],
                ProbeLocationLimits {
                    absolute_tolerance_m: 0.0,
                    ..limits
                },
            ),
        ] {
            assert!(
                model
                    .probe_at_position(0, point, bad_limits, &mut source, &mut || false)
                    .is_err()
            );
            assert_eq!(calls.get(), before);
        }
        assert!(matches!(
            model.probe_at_position(0, map([0.2, 0.3, 0.1]), limits, &mut source, &mut || true),
            Err(AcousticFeError::Cancelled)
        ));
        assert_eq!(calls.get(), before);
    });
}

fn simplex_moment(powers: [u32; 3]) -> f64 {
    let factorial = |n: u32| (1..=n).map(f64::from).product::<f64>();
    powers.into_iter().map(factorial).product::<f64>()
        / factorial(powers.into_iter().sum::<u32>() + 3)
}

#[test]
fn curved_volume_source_uses_physical_density_and_variable_volume_jacobian() {
    with_cx(|cx| {
        // X=(xi,eta,zeta*(1+a*xi)), detJ=1+a*xi. The source's spatial
        // dependence and the variable volume measure are independent axes.
        let a = 0.8;
        let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
        let (positions, midpoints) = mapped_nodes(&mesh, |q| [q[0], q[1], q[2] * (1.0 + a * q[0])]);
        let fe = HelmholtzFe::new_quadratic(
            &mesh,
            &positions,
            &midpoints,
            2,
            8,
            MEDIUM,
            QuadraticGeometryControls::default(),
            cx,
        )
        .unwrap();
        let model = uniform(&fe);
        let z_field = physical_affine_field(&fe, &midpoints, |x| C64::from_re(x[2]));
        let z_terms = [(1.0, [0, 0, 1]), (a, [1, 0, 1])];
        let source_terms = [
            (1.0, [0, 0, 0]),
            (2.0, [1, 0, 0]),
            (-3.0, [0, 1, 0]),
            (4.0, [0, 0, 1]),
            (4.0 * a, [1, 0, 1]),
        ];
        let measure_terms = [(1.0, [0, 0, 0]), (a, [1, 0, 0])];
        let mut exact_work = 0.0;
        for (z, iz) in z_terms {
            for (s, is) in source_terms {
                for (j, ij) in measure_terms {
                    exact_work += z
                        * s
                        * j
                        * simplex_moment(core::array::from_fn(|axis| {
                            iz[axis] + is[axis] + ij[axis]
                        }));
                }
            }
        }
        let exact_integral = 7.0 / 24.0 + 7.0 * a / 60.0 + a * a / 90.0;
        for quadrature in [4, 6, 8] {
            let source = model
                .volume_velocity_source(
                    quadrature,
                    &mut |_, x| {
                        Ok::<_, HeterogeneousError>(1.0 + 2.0 * x[0] - 3.0 * x[1] + 4.0 * x[2])
                    },
                    &mut || false,
                )
                .unwrap();
            close(source.integrated_shape(), exact_integral, 3e-14);
            close(source.harmonic_pressure(&z_field).re, exact_work, 3e-14);
            let uniform = model
                .volume_velocity_source(
                    quadrature,
                    &mut |_, _| Ok::<_, HeterogeneousError>(1.0),
                    &mut || false,
                )
                .unwrap();
            close(uniform.integrated_shape(), 1.0 / 6.0 + a / 24.0, 2e-14);
        }
        assert!(
            (exact_integral - 7.0 / 24.0).abs() > 0.09,
            "a constant corner volume or affine position must change the result"
        );
    });
}

#[test]
fn curved_boundary_source_and_surface_power_resolve_area_and_normal_independently() {
    with_cx(|cx| {
        // The bottom face is Z=2 X Y, with outward n dS=(2Y,2X,-1)dX dY.
        // Independent one-dimensional integrals (analytic inner y integral,
        // composite Simpson 4096/16384) give the constants below.
        const AREA: f64 = 0.751_156_358_570_360_3;
        const SHAPE_AREA: f64 = 0.266_601_992_364_454_75;
        const SHAPE_X_AREA: f64 = 0.109_131_401_902_760_1;
        let a = 2.0;
        let mesh = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
        let (positions, midpoints) = mapped_nodes(&mesh, |q| [q[0], q[1], q[2] + a * q[0] * q[1]]);
        let fe = HelmholtzFe::new_quadratic(
            &mesh,
            &positions,
            &midpoints,
            2,
            8,
            MEDIUM,
            QuadraticGeometryControls::default(),
            cx,
        )
        .unwrap();
        let model = uniform(&fe);
        let face = [0, 1, 2];
        let trace = [(face, [0.0, 0.0, 0.0, 0.0, 1.0, 0.0])]; // g=4 X Y
        let x_field = physical_affine_field(&fe, &midpoints, |x| C64::from_re(x[0]));
        let mut errors = Vec::new();
        for order in [4, 8, 16] {
            let source = model
                .quadratic_velocity_source_with_quadrature(&trace, order, &mut || false)
                .unwrap();
            let surface = model
                .surface_probe(
                    &[ElementFace {
                        element: 0,
                        vertices: face,
                    }],
                    order,
                    SurfaceProbeLimits::default(),
                    &mut |_, _| Ok::<_, AcousticFeError>(MEDIUM),
                    &mut || false,
                )
                .unwrap();
            errors.push(
                (surface.area_m2() - AREA).abs()
                    + (source.integrated_shape_m2() - SHAPE_AREA).abs(),
            );
            if order == 16 {
                close(surface.area_m2(), AREA, 2e-12);
                close(source.integrated_shape_m2(), SHAPE_AREA, 2e-12);
                close(source.harmonic_pressure(&x_field).re, SHAPE_X_AREA, 2e-12);
            }
        }
        assert!(
            errors[1] < errors[0] && errors[2] < errors[1],
            "quadrature ladder={errors:?}"
        );
        assert!(
            AREA - 0.5 > 0.25,
            "the affine corner triangle must be a live contrast"
        );
        let default = model
            .quadratic_velocity_source(&trace, &mut || false)
            .unwrap();
        let same = model
            .quadratic_velocity_source_with_quadrature(&trace, 8, &mut || false)
            .unwrap();
        assert_eq!(
            default.load, same.load,
            "default must retain the declared geometry order"
        );
        assert_eq!(default.support_faces(), &[face]);
        assert!(default.has_nonzero_projection());
        // The antisymmetric shape X-Y has zero net integral on this symmetric
        // curved face, but remains a live pressure-space excitation. A zero
        // second face must also remain visible to boundary-law admission.
        let signed = model
            .quadratic_velocity_source_with_quadrature(
                &[
                    ([3, 0, 1], [0.0; 6]),
                    ([2, 1, 0], [-1.0, 1.0, 0.0, 0.0, 0.5, -0.5]),
                ],
                16,
                &mut || false,
            )
            .unwrap();
        assert_eq!(signed.support_faces(), &[[0, 1, 2], [0, 1, 3]]);
        close(signed.integrated_shape_m2(), 0.0, 2e-12);
        assert!(signed.has_nonzero_projection());
        let zero = model
            .quadratic_velocity_source(&[(face, [0.0; 6])], &mut || false)
            .unwrap();
        assert_eq!(zero.support_faces(), &[face]);
        assert!(!zero.has_nonzero_projection());
        let reversed = model
            .quadratic_velocity_source_with_quadrature(
                &[([2, 1, 0], [0.0, 0.0, 0.0, 1.0, 0.0, 0.0])],
                8,
                &mut || false,
            )
            .unwrap();
        for (&x, &y) in same.load.iter().zip(&reversed.load) {
            close(x, y, 2e-14);
        }
        assert!(
            model
                .quadratic_velocity_source_with_quadrature(&trace, 3, &mut || false)
                .is_err()
        );
        assert!(
            model
                .quadratic_velocity_source_with_quadrature(&trace, 17, &mut || false)
                .is_err()
        );
        let surface = model
            .surface_probe(
                &[ElementFace {
                    element: 0,
                    vertices: face,
                }],
                8,
                SurfaceProbeLimits::default(),
                &mut |_, _| Ok::<_, AcousticFeError>(MEDIUM),
                &mut || false,
            )
            .unwrap();
        let omega = 2.0;
        let velocity = 0.7;
        let pressure = physical_affine_field(&fe, &midpoints, |x| {
            C64::new(1.0, -omega * MEDIUM.density_kg_m3 * velocity * x[2])
        });
        let power = surface
            .harmonic(&pressure, omega / core::f64::consts::TAU, &mut || false)
            .unwrap();
        close(power.active_power_w, -velocity / 4.0, 2e-12);
        close(
            power.reactive_power_w,
            omega * MEDIUM.density_kg_m3 * velocity * velocity * a / 48.0,
            2e-12,
        );
        let faces = [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]].map(|vertices| ElementFace {
            element: 0,
            vertices,
        });
        let closed = model
            .surface_probe(
                &faces,
                8,
                SurfaceProbeLimits::default(),
                &mut |_, _| Ok::<_, AcousticFeError>(MEDIUM),
                &mut || false,
            )
            .unwrap();
        let closed_power = closed
            .harmonic(&pressure, omega / core::f64::consts::TAU, &mut || false)
            .unwrap();
        close(closed_power.active_power_w, 0.0, 2e-12);
        close(
            closed_power.reactive_power_w,
            -omega * MEDIUM.density_kg_m3 * velocity * velocity / 12.0,
            2e-12,
        );
        let owner = model
            .prepare_time(0.01, &[], HeterogeneousTimeLimits::default(), &mut || false)
            .unwrap();
        let potential = physical_affine_field(&fe, &midpoints, |x| {
            C64::from_re(17.0 - MEDIUM.density_kg_m3 * velocity * x[2])
        })
        .into_iter()
        .map(|v| v.re)
        .collect::<Vec<_>>();
        let pressure = physical_affine_field(&fe, &midpoints, |_| C64::ONE)
            .into_iter()
            .map(|v| v.re)
            .collect::<Vec<_>>();
        let state = owner.state(&potential, &pressure).unwrap();
        close(
            surface.transient(&state, &mut || false).unwrap(),
            -velocity / 2.0,
            2e-12,
        );
        close(closed.transient(&state, &mut || false).unwrap(), 0.0, 2e-12);
    });
}
