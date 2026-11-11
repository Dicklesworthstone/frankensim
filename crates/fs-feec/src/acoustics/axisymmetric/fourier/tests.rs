//! Independent cylindrical Bessel fields, full-turn energy and angular phase.
//! These are G0/G1/G3 discretization/implementation checks, not measured validation.
use super::*;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};

const RHO: f64 = 1.3;
const SPEED: f64 = 2.0;
const LIMITS: WeightedAssemblyLimits = WeightedAssemblyLimits {
    max_cells: 2048,
    max_dofs: 2048,
    max_triplets: 18432,
};
fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 91,
                kernel_id: 24,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}
fn controls() -> HarmonicResponseControls {
    HarmonicResponseControls {
        maximum_coordinates: 256,
        maximum_payload_bytes: 32 << 20,
        maximum_relative_residual: 5e-12,
    }
}
fn quadrature() -> FourierQuadrature {
    FourierQuadrature {
        points: 12,
        maximum_relative_difference: 1e-10,
    }
}
fn materials(mesh: &TriComplex2) -> Vec<AxisymmetricFluid> {
    vec![
        AxisymmetricFluid {
            density_kg_m3: RHO,
            compressibility_pa_inv: 1.0 / (RHO * SPEED * SPEED)
        };
        mesh.faces().len()
    ]
}
fn rectangle(a: f64, b: f64, length: f64, nr: usize, nz: usize, reverse: bool) -> TriComplex2 {
    let mut vertices = Vec::new();
    for z in 0..=nz {
        for r in 0..=nr {
            vertices.push([
                a + (b - a) * r as f64 / nr as f64,
                length * z as f64 / nz as f64,
            ]);
        }
    }
    let mut faces = Vec::new();
    for z in 0..nz {
        for r in 0..nr {
            let p = (z * (nr + 1) + r) as u32;
            for mut face in [
                [p, p + 1, p + nr as u32 + 2],
                [p, p + nr as u32 + 2, p + nr as u32 + 1],
            ] {
                if reverse {
                    face.swap(1, 2);
                }
                faces.push(face);
            }
        }
    }
    TriComplex2::from_indexed_triangles(
        "fourier-acoustic-fixture",
        vertices,
        faces,
        Metric2::axisymmetric(core::f64::consts::TAU).unwrap(),
    )
    .unwrap()
}
fn edges(model: &AxisymmetricPressureFe<'_>, axis: usize, coordinate: f64) -> Vec<[u32; 2]> {
    model
        .boundary_edges()
        .iter()
        .filter(|edge| {
            edge.vertices
                .iter()
                .all(|&v| model.mesh().vertices()[v as usize][axis] == coordinate)
        })
        .map(|e| e.vertices)
        .collect()
}
fn drive(edges: &[[u32; 2]], value: C64) -> Vec<AxisymmetricVelocityTrace> {
    edges.iter().map(|&e| (e, [value; 2])).collect()
}
fn near(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1.0),
        "{actual:?} != {expected:?}"
    );
}

// The convergent defining series, independent of the FE matrices, radial
// quadrature and native harmonic solve. Here |x|<=3.1 and m is 1 or 2.
fn bessel_j(order: usize, x: f64) -> f64 {
    let mut term = 1.0;
    for n in 1..=order {
        term *= 0.5 * x / n as f64;
    }
    let mut sum = term;
    for n in 1..48 {
        term *= -0.25 * x * x / (n * (n + order)) as f64;
        sum += term;
    }
    sum
}

#[test]
fn original_m_zero_is_unchanged_and_angular_forms_keep_exact_swept_moments() {
    with_cx(&CancelGate::new(), |cx| {
        for inner in [0.0, 0.2] {
            for reverse in [false, true] {
                let mesh = rectangle(inner, 1.0, 0.4, 3, 2, reverse);
                let data = materials(&mesh);
                let original = AxisymmetricPressureFe::new(&mesh, &data, LIMITS, cx).unwrap();
                let zero =
                    AxisymmetricPressureFe::new_fourier(&mesh, &data, 0, quadrature(), LIMITS, cx)
                        .unwrap();
                assert_eq!(original.mass(), zero.mass());
                assert_eq!(original.stiffness(), zero.stiffness());
                assert!(zero.axis_pressure_constraints().is_empty());
                assert!(zero.azimuthal_stiffness().is_none());
                let trace = drive(&edges(&original, 1, 0.0), C64::new(0.3, 0.2));
                let baseline = original
                    .prepare(0.3, &[], &[], controls(), cx)
                    .unwrap()
                    .solve(&trace, cx)
                    .unwrap();
                let result = zero
                    .prepare(0.3, &[], &[], controls(), cx)
                    .unwrap()
                    .solve(&trace, cx)
                    .unwrap();
                assert_eq!(baseline.pressure_pa(), result.pressure_pa());
                assert_eq!(baseline.kinetic_energy_j, result.kinetic_energy_j);
                assert_eq!(result.azimuthal_kinetic_energy_j, 0.0);
                for order in [1, -1, 3] {
                    let model = AxisymmetricPressureFe::new_fourier(
                        &mesh,
                        &data,
                        order,
                        quadrature(),
                        LIMITS,
                        cx,
                    )
                    .unwrap();
                    let radial = mesh
                        .vertices()
                        .iter()
                        .map(|p| C64::from_re(p[0]))
                        .collect::<Vec<_>>();
                    // |d_phi(r exp(i*m*phi))/r|²=m² at every physical point.
                    let expected = f64::from(order).powi(2) * model.volume_m3() / RHO;
                    let actual =
                        quadratic(model.azimuthal_stiffness().unwrap(), &radial, cx).unwrap();
                    assert!((actual / expected - 1.0).abs() < 5e-12);
                    let total = quadratic(model.stiffness(), &radial, cx).unwrap();
                    assert!((total / (expected + model.volume_m3() / RHO) - 1.0).abs() < 5e-12);
                    assert_eq!(
                        model.axis_pressure_constraints().len(),
                        if inner == 0.0 { 3 } else { 0 }
                    );
                    assert!(model.maximum_angular_quadrature_difference() < 1e-10);
                }
            }
        }
    });
}

#[test]
fn nonaxisymmetric_cylinder_fields_and_bessel_mode_frequencies_converge_in_h() {
    with_cx(&CancelGate::new(), |cx| {
        let k = 0.85;
        let omega = SPEED * k;
        // First positive zeros of J'_1 and J'_2; no FE fit determines them.
        for (m, root) in [(1_usize, 1.8411837813406593), (2, 3.0542369282271404)] {
            let derivative = 0.5 * k * (bessel_j(m - 1, k) - bessel_j(m + 1, k));
            let mut previous = [f64::INFINITY; 2];
            for nr in [4, 8, 16] {
                // h-refinement in BOTH directions. With nz fixed at 2 the axial cell size stays 0.2 while nr doubles, and
                // the first run of this test measured a rate of 1.22 at nr = 8 (band 1.5..2.5 below). nz = nr / 2 refines
                // the axial direction with the radial one and keeps the cell aspect h_z / h_r = 0.8 of the nr = 4 mesh at
                // every level.
                let mesh = rectangle(0.0, 1.0, 0.4, nr, nr / 2, false);
                let model = AxisymmetricPressureFe::new_fourier(
                    &mesh,
                    &materials(&mesh),
                    m as i32,
                    quadrature(),
                    LIMITS,
                    cx,
                )
                .unwrap();
                let trace = drive(&edges(&model, 0, 1.0), C64::from_re(derivative / RHO));
                let prepared = model
                    .prepare(omega / core::f64::consts::TAU, &[], &[], controls(), cx)
                    .unwrap();
                let field = prepared.solve(&trace, cx).unwrap();
                let reference = mesh
                    .vertices()
                    .iter()
                    .map(|p| C64::new(0.0, omega * bessel_j(m, k * p[0])))
                    .collect::<Vec<_>>();
                let error = field
                    .pressure_pa()
                    .iter()
                    .zip(&reference)
                    .map(|(&p, &q)| p - q)
                    .collect::<Vec<_>>();
                let field_error = (quadratic(model.mass(), &error, cx).unwrap()
                    / quadratic(model.mass(), &reference, cx).unwrap())
                .sqrt();
                let mode = mesh
                    .vertices()
                    .iter()
                    .map(|p| C64::from_re(bessel_j(m, root * p[0])))
                    .collect::<Vec<_>>();
                let rayleigh = quadratic(model.stiffness(), &mode, cx).unwrap()
                    / quadratic(model.mass(), &mode, cx).unwrap();
                let mode_error = (rayleigh.sqrt() / (SPEED * root) - 1.0).abs();
                for (axis, error) in [field_error, mode_error].into_iter().enumerate() {
                    if previous[axis].is_finite() {
                        let rate = (previous[axis] / error).log2();
                        assert!(
                            (1.5..2.5).contains(&rate),
                            "m={m}, nr={nr}, axis={axis}, rate={rate}"
                        );
                    }
                    previous[axis] = error;
                }
                for &axis in model.axis_pressure_constraints() {
                    assert!(!prepared.free_dofs().contains(&axis));
                    assert_eq!(field.pressure_pa()[axis], C64::ZERO);
                    assert_eq!(field.pressure_potential_pa_s()[axis], C64::ZERO);
                }
                assert_eq!(field.pressure_pa().len(), model.ndof());
                assert!(field.azimuthal_kinetic_energy_j > 0.0);
                assert!(field.kinetic_energy_j > field.azimuthal_kinetic_energy_j);
            }
            assert!(
                previous[0] < 0.005 && previous[1] < 0.005,
                "m={m}: {previous:?}"
            );
        }
    });
}

#[test]
fn native_boundary_work_signed_orders_and_rotated_observations_agree() {
    with_cx(&CancelGate::new(), |cx| {
        let mesh = rectangle(0.2, 1.0, 0.4, 4, 3, false);
        let data = materials(&mesh);
        let plus =
            AxisymmetricPressureFe::new_fourier(&mesh, &data, 2, quadrature(), LIMITS, cx).unwrap();
        let minus = AxisymmetricPressureFe::new_fourier(&mesh, &data, -2, quadrature(), LIMITS, cx)
            .unwrap();
        assert_eq!(plus.stiffness(), minus.stiffness());
        let walls = [AxisymmetricAdmittance {
            edges: plus.boundary_edges().iter().map(|e| e.vertices).collect(),
            admittance_m_pa_s: 0.2 / (RHO * SPEED),
        }];
        let trace = drive(&edges(&plus, 1, 0.0), C64::new(0.7, -0.3));
        let positive = plus
            .prepare(0.4, &walls, &[], controls(), cx)
            .unwrap()
            .solve(&trace, cx)
            .unwrap();
        let negative = minus
            .prepare(0.4, &walls, &[], controls(), cx)
            .unwrap()
            .solve(&trace, cx)
            .unwrap();
        assert_eq!(positive.pressure_pa(), negative.pressure_pa());
        assert_eq!(positive.kinetic_energy_j, negative.kinetic_energy_j);
        assert!(positive.boundary_dissipation_w > 0.0 && positive.azimuthal_kinetic_energy_j > 0.0);
        assert!(positive.power_balance_residual_w().abs() < 1e-10 * positive.input_power_w.abs());
        let p = positive.probe_fourier(7, [0.2, 0.3, 0.5], 0.0).unwrap();
        let n = negative.probe_fourier(7, [0.2, 0.3, 0.5], 0.0).unwrap();
        near(
            p.particle_velocity_r_phi_z_m_s[1],
            n.particle_velocity_r_phi_z_m_s[1].scale(-1.0),
            1e-13,
        );
        let shift = 0.63_f64;
        let phase = C64::new((-2.0 * shift).cos(), (-2.0 * shift).sin());
        let shifted = trace
            .iter()
            .map(|&(edge, v)| (edge, v.map(|x| x * phase)))
            .collect::<Vec<_>>();
        let rotated = plus
            .prepare(0.4, &walls, &[], controls(), cx)
            .unwrap()
            .solve(&shifted, cx)
            .unwrap();
        let q = rotated.probe_fourier(7, [0.2, 0.3, 0.5], shift).unwrap();
        near(q.pressure_pa, p.pressure_pa, 3e-12);
        let [vx, vy, vz] = p.particle_velocity_xyz_m_s;
        near(
            q.particle_velocity_xyz_m_s[0],
            vx.scale(shift.cos()) - vy.scale(shift.sin()),
            3e-12,
        );
        near(
            q.particle_velocity_xyz_m_s[1],
            vx.scale(shift.sin()) + vy.scale(shift.cos()),
            3e-12,
        );
        near(q.particle_velocity_xyz_m_s[2], vz, 3e-12);
        assert!((rotated.kinetic_energy_j / positive.kinetic_energy_j - 1.0).abs() < 1e-11);
        assert!(
            (rotated.boundary_dissipation_w / positive.boundary_dissipation_w - 1.0).abs() < 1e-11
        );
        for axis in 0..3 {
            let expected = 0.5 * (q.pressure_pa * q.particle_velocity_xyz_m_s[axis].conj()).re;
            assert_eq!(q.active_intensity_xyz_w_m2[axis], expected);
        }
    });
}

// This analytic nodal field tests reconstruction in isolation: psi=r gives
// psi*exp(±i phi)=x±iy, whose Cartesian gradient is exactly (1,±i,0).
fn linear_field<'f, 'm>(
    model: &'f AxisymmetricPressureFe<'m>,
    amplitude: f64,
) -> AxisymmetricHarmonicField<'f, 'm> {
    let potential = model
        .mesh()
        .vertices()
        .iter()
        .map(|p| C64::from_re(amplitude * p[0]))
        .collect::<Vec<_>>();
    let pressure = potential.iter().map(|&p| C64::new(0.0, 2.0) * p).collect();
    AxisymmetricHarmonicField {
        model,
        potential,
        pressure,
        input_power_w: 0.0,
        boundary_dissipation_w: 0.0,
        kinetic_energy_j: 0.0,
        azimuthal_kinetic_energy_j: 0.0,
        compressional_energy_j: 0.0,
        relative_residual: 0.0,
    }
}

#[test]
fn coherent_positive_and_negative_components_reconstruct_real_cosine_without_double_power() {
    with_cx(&CancelGate::new(), |cx| {
        let mesh = rectangle(0.0, 1.0, 0.4, 2, 2, false);
        let data = materials(&mesh);
        let plus =
            AxisymmetricPressureFe::new_fourier(&mesh, &data, 1, quadrature(), LIMITS, cx).unwrap();
        let minus = AxisymmetricPressureFe::new_fourier(&mesh, &data, -1, quadrature(), LIMITS, cx)
            .unwrap();
        let positive = linear_field(&plus, 0.5);
        let negative = linear_field(&minus, 0.5);
        for angle in [0.0, 0.7, 2.3, 5.2] {
            let a = positive.probe_fourier(3, [0.2, 0.3, 0.5], angle).unwrap();
            let b = negative.probe_fourier(3, [0.2, 0.3, 0.5], angle).unwrap();
            near(
                a.pressure_potential_pa_s + b.pressure_potential_pa_s,
                C64::from_re(a.position_xyz_m[0]),
                2e-13,
            );
            near(
                a.particle_velocity_xyz_m_s[0] + b.particle_velocity_xyz_m_s[0],
                C64::from_re(-1.0 / RHO),
                2e-13,
            );
            near(
                a.particle_velocity_xyz_m_s[1] + b.particle_velocity_xyz_m_s[1],
                C64::ZERO,
                2e-13,
            );
            near(
                a.particle_velocity_xyz_m_s[2] + b.particle_velocity_xyz_m_s[2],
                C64::ZERO,
                2e-13,
            );
        }
        // Integrating |grad x|² over the original physical volume gives this
        // energy. Each independently normalized signed component has half.
        let sum = 0.25
            * (quadratic(plus.stiffness(), &positive.potential, cx).unwrap()
                + quadratic(minus.stiffness(), &negative.potential, cx).unwrap());
        assert!((sum / (plus.volume_m3() / (4.0 * RHO)) - 1.0).abs() < 1e-12);
        assert!(positive.probe_fourier(0, [1.0, 0.0, 0.0], 0.0).is_err());
        assert_eq!(
            positive.probe(0, [1.0, 0.0, 0.0]).unwrap().pressure_pa,
            C64::ZERO
        );
    });
}

#[test]
fn native_admission_quadrature_refinement_and_axis_restrictions_remain_explicit() {
    with_cx(&CancelGate::new(), |cx| {
        let mesh = rectangle(0.0, 1.0, 0.4, 2, 2, false);
        let data = materials(&mesh);
        for order in [i32::MIN, MAX_AZIMUTHAL_ORDER as i32 + 1] {
            assert!(
                AxisymmetricPressureFe::new_fourier(&mesh, &data, order, quadrature(), LIMITS, cx)
                    .is_err()
            );
        }
        assert!(
            AxisymmetricPressureFe::new_fourier(
                &mesh,
                &data,
                MAX_AZIMUTHAL_ORDER as i32,
                quadrature(),
                LIMITS,
                cx
            )
            .is_ok()
        );
        for invalid in [
            FourierQuadrature {
                points: 1,
                ..quadrature()
            },
            FourierQuadrature {
                points: MAX_FOURIER_QUADRATURE_POINTS + 1,
                ..quadrature()
            },
            FourierQuadrature {
                maximum_relative_difference: f64::NAN,
                ..quadrature()
            },
        ] {
            assert!(
                AxisymmetricPressureFe::new_fourier(&mesh, &data, 1, invalid, LIMITS, cx).is_err()
            );
        }
        let model =
            AxisymmetricPressureFe::new_fourier(&mesh, &data, 2, quadrature(), LIMITS, cx).unwrap();
        let release = edges(&model, 1, 0.4);
        let prepared = model.prepare(0.2, &[], &release, controls(), cx).unwrap();
        let solved = prepared
            .solve(&drive(&edges(&model, 1, 0.0), C64::ONE), cx)
            .unwrap();
        for (i, p) in mesh.vertices().iter().enumerate() {
            if p[0] == 0.0 || p[1] == 0.4 {
                assert_eq!(solved.pressure_pa()[i], C64::ZERO);
            }
        }
        assert!(solved.probe_fourier(0, [0.2, 0.3, 0.5], f64::NAN).is_err());
        assert!(
            model
                .prepare(
                    0.2,
                    &[],
                    &[],
                    HarmonicResponseControls {
                        maximum_coordinates: 1,
                        ..controls()
                    },
                    cx
                )
                .is_err()
        );
        // A strongly varying original radius needs declared quadrature work;
        // lower q is not replaced with an averaged radius or a successful flag.
        let wide = rectangle(0.05, 1.0, 0.4, 1, 1, false);
        let coefficients = materials(&wide);
        let mut previous = f64::INFINITY;
        for points in [2, 4, 8] {
            let m = AxisymmetricPressureFe::new_fourier(
                &wide,
                &coefficients,
                1,
                FourierQuadrature {
                    points,
                    maximum_relative_difference: 0.9,
                },
                LIMITS,
                cx,
            )
            .unwrap();
            let ones = vec![C64::ONE; m.ndof()];
            let integral = quadratic(m.azimuthal_stiffness().unwrap(), &ones, cx).unwrap();
            let exact = core::f64::consts::TAU * 0.4 * (1.0 / 0.05_f64).ln() / RHO;
            let error = (integral / exact - 1.0).abs();
            assert!(error < previous);
            previous = error;
        }
        assert!(previous < 1e-5);
        assert!(
            AxisymmetricPressureFe::new_fourier(
                &wide,
                &coefficients,
                1,
                FourierQuadrature {
                    points: 2,
                    maximum_relative_difference: 1e-12
                },
                LIMITS,
                cx
            )
            .is_err()
        );
    });
    let gate = CancelGate::new();
    gate.request();
    with_cx(&gate, |cx| {
        let mesh = rectangle(0.0, 1.0, 0.4, 1, 1, false);
        assert!(matches!(
            AxisymmetricPressureFe::new_fourier(
                &mesh,
                &materials(&mesh),
                1,
                quadrature(),
                LIMITS,
                cx
            ),
            Err(AxisymmetricError::Cancelled)
        ));
    });
}
