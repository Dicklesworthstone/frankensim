use super::*;
use crate::acoustics::helmholtz::HelmholtzFe;
use crate::acoustics::helmholtz::heterogeneous::time::HeterogeneousTimeLimits;
use crate::highorder::simplex::LocalFn;
use fs_math::det;
use fs_rep_mesh::TetComplex;

const MEDIUM: Medium = Medium {
    density_kg_m3: 4.7,
    sound_speed_m_s: 1.4,
};
const GEOMETRY_MEDIUM: Medium = Medium {
    density_kg_m3: 0.7,
    sound_speed_m_s: 3.2,
};
const POSITIONS: [[f64; 3]; 4] = [
    [0.2, -0.3, 0.4],
    [2.2, 0.7, 0.4],
    [0.2, 2.7, 0.4],
    [1.2, -0.3, 4.4],
];
// Analytic gradients for this skew affine tetrahedron, in GLOBAL vertex order.
const GRAD_LAMBDA: [[f64; 3]; 4] = [
    [-1.0 / 3.0, -1.0 / 3.0, -1.0 / 6.0],
    [0.5, 0.0, -0.125],
    [-1.0 / 6.0, 1.0 / 3.0, 1.0 / 24.0],
    [0.0, 0.0, 0.25],
];

fn uniform<'f, 'c>(geometry: &'f HelmholtzFe<'c>) -> HeterogeneousHelmholtzFe<'f, 'c> {
    HeterogeneousHelmholtzFe::new(
        geometry,
        geometry.order() + 2,
        &mut |_, _| Ok::<_, AcousticFeError>(MEDIUM),
        &mut || false,
    )
    .unwrap()
}

fn close(actual: C64, expected: C64) {
    assert!((actual - expected).abs() < 2e-13, "{actual:?} != {expected:?}");
}

#[test]
fn high_order_probes_recover_physical_polynomial_gradients_on_reordered_skew_tets() {
    let lambda = [0.13, 0.21, 0.29, 0.37];
    let linear = [
        C64::new(1.2, -0.4),
        C64::new(-0.7, 0.6),
        C64::new(2.1, -0.2),
        C64::new(0.3, 1.1),
    ];
    let [a, b, c, d, e, f] = [
        C64::new(0.4, -0.8),
        C64::new(-0.6, 0.2),
        C64::new(0.7, 0.9),
        C64::new(-1.1, 0.3),
        C64::new(0.2, -0.5),
        C64::new(-0.9, -0.4),
    ];
    for tet in [[0, 1, 2, 3], [3, 1, 0, 2]] {
        for order in 1..=4 {
            let complex = TetComplex::from_tets(4, vec![tet]);
            let geometry = HelmholtzFe::new(&complex, &POSITIONS, order, GEOMETRY_MEDIUM).unwrap();
            let model = uniform(&geometry);
            let mut coefficients = vec![C64::ZERO; geometry.ndof()];
            for (dof, shape) in geometry.space.element_dofs(0) {
                coefficients[dof] = match shape {
                    LocalFn::Vertex(v) => linear[tet[v] as usize],
                    LocalFn::Edge(u, v, k) if [tet[u], tet[v]] == [0, 1] => [a, b, d][k],
                    LocalFn::Face(u, v, w, 0, 0) if [tet[u], tet[v], tet[w]] == [0, 1, 2] => c,
                    LocalFn::Face(u, v, w, 1, 0) if [tet[u], tet[v], tet[w]] == [0, 1, 2] => e,
                    LocalFn::Interior(0, 0, 0) => f,
                    _ => C64::ZERO,
                };
            }
            // Independent expanded monomials, including nonconstant edge/face
            // kernels and the P4 cell bubble. No FE derivative routine is used.
            let mut terms = vec![
                (linear[0], [1, 0, 0, 0]),
                (linear[1], [0, 1, 0, 0]),
                (linear[2], [0, 0, 1, 0]),
                (linear[3], [0, 0, 0, 1]),
            ];
            if order >= 2 {
                terms.push((a, [1, 1, 0, 0]));
            }
            if order >= 3 {
                terms.extend([(b, [1, 2, 0, 0]), (-b, [2, 1, 0, 0]), (c, [1, 1, 1, 0])]);
            }
            if order >= 4 {
                terms.extend([
                    (d.scale(1.5), [1, 3, 0, 0]),
                    (d.scale(-3.0), [2, 2, 0, 0]),
                    (d.scale(1.5), [3, 1, 0, 0]),
                    (d.scale(-0.5), [1, 1, 0, 0]),
                    (e, [1, 2, 1, 0]),
                    (-e, [2, 1, 1, 0]),
                    (f, [1, 1, 1, 1]),
                ]);
            }
            let mut exact_p = C64::ZERO;
            let mut exact_gradient = [C64::ZERO; 3];
            for (coefficient, powers) in terms {
                let monomial = |powers: [i32; 4]| -> f64 {
                    lambda.iter().zip(powers).map(|(&x, power)| f64::powi(x, power)).product()
                };
                exact_p = exact_p + coefficient.scale(monomial(powers));
                for v in 0..4 {
                    if powers[v] == 0 {
                        continue;
                    }
                    let mut differentiated = powers;
                    differentiated[v] -= 1;
                    let derivative = f64::from(powers[v]) * monomial(differentiated);
                    for axis in 0..3 {
                        exact_gradient[axis] = exact_gradient[axis]
                            + coefficient.scale(derivative * GRAD_LAMBDA[v][axis]);
                    }
                }
            }
            let barycentric = tet.map(|v| lambda[v as usize]);
            let mut calls = 0;
            let probe = model
                .probe(0, barycentric, &mut |element, position| {
                    calls += 1;
                    assert_eq!(element, 0);
                    for axis in 0..3 {
                        let expected: f64 = (0..4).map(|v| lambda[v] * POSITIONS[v][axis]).sum();
                        assert!((position[axis] - expected).abs() < 5e-16);
                    }
                    Ok::<_, AcousticFeError>(MEDIUM)
                })
                .unwrap();
            assert_eq!(probe.medium(), MEDIUM);
            assert_eq!(probe.element(), 0);
            let frequency = 1.25;
            let sample = probe.harmonic(&coefficients, frequency).unwrap();
            close(sample.pressure_pa, exact_p);
            close(sample.pressure_pa, model.pressure_at(&coefficients, 0, barycentric).unwrap());
            for axis in 0..3 {
                close(sample.pressure_gradient_pa_m[axis], exact_gradient[axis]);
                let expected_v = C64::new(-exact_gradient[axis].im, exact_gradient[axis].re)
                    .scale(1.0 / (core::f64::consts::TAU * frequency * MEDIUM.density_kg_m3));
                close(sample.particle_velocity_m_s[axis], expected_v);
                let expected_i = exact_p.scale(0.5) * expected_v.conj();
                close(
                    C64::new(sample.active_intensity_w_m2[axis], sample.reactive_intensity_w_m2[axis]),
                    expected_i,
                );
            }
            assert_eq!(probe.harmonic(&coefficients, frequency).unwrap(), sample);
            assert_eq!(calls, 1, "repeated fields must not resample material");
        }
    }
}

#[test]
fn wave_gradients_preserve_refracted_power_and_evanescent_reactive_flow() {
    let complex = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let geometry = HelmholtzFe::new(&complex, &POSITIONS, 1, GEOMETRY_MEDIUM).unwrap();
    let receiver = Medium { density_kg_m3: 2.0, sound_speed_m_s: 2.0 };
    let mut source = |_, _| Ok::<_, AcousticFeError>(receiver);
    let model = HeterogeneousHelmholtzFe::new(&geometry, 3, &mut source, &mut || false).unwrap();
    let probe = model.probe(0, [0.25; 4], &mut source).unwrap();
    let center = probe.position_m();
    let omega = 2.0;
    let pressure = C64::new(1.2, -0.7);
    for (normal_k, tangent_k) in [
        (C64::from_re(0.8), 0.6),
        (C64::new(0.0, -det::sqrt(0.96)), 1.4),
    ] {
        // Exact first-order Taylor jet of the transmitted Snell wave at the
        // observation point: interpolation/gradient are exact, with no mesh
        // dispersion allowance masking sign, direction or normalization errors.
        let gradient = [
            C64::new(normal_k.im, -normal_k.re) * pressure,
            C64::new(0.0, -tangent_k) * pressure,
            C64::ZERO,
        ];
        let coefficients = POSITIONS.map(|point| {
            (0..3).fold(pressure, |p, axis| p + gradient[axis].scale(point[axis] - center[axis]))
        });
        let sample = probe.harmonic(&coefficients, omega / core::f64::consts::TAU).unwrap();
        let base = pressure.norm_sq() / (2.0 * omega * receiver.density_kg_m3);
        let expected_active = [normal_k.re * base, tangent_k * base, 0.0];
        let expected_reactive = [-normal_k.im * base, 0.0, 0.0];
        for axis in 0..3 {
            close(C64::from_re(sample.active_intensity_w_m2[axis]), C64::from_re(expected_active[axis]));
            close(C64::from_re(sample.reactive_intensity_w_m2[axis]), C64::from_re(expected_reactive[axis]));
        }
    }
    // A real standing-wave pressure/gradient stores reactive energy, while
    // its time-averaged transported power vanishes in every direction.
    let coefficients = POSITIONS.map(|p| C64::from_re(2.0 + 0.3 * p[0]));
    let sample = probe.harmonic(&coefficients, omega / core::f64::consts::TAU).unwrap();
    assert!(sample.active_intensity_w_m2.iter().all(|i| i.abs() < 1e-14));
    let expected = -(2.0 + 0.3 * center[0]) * 0.3 / (2.0 * omega * receiver.density_kg_m3);
    assert!((sample.reactive_intensity_w_m2[0] - expected).abs() < 1e-14);
}

#[test]
fn interface_probes_keep_one_sided_density_and_continuous_normal_velocity() {
    let positions = [
        [0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0],
        [-1.0, 0.0, 0.0], [1.0, 0.0, 0.0],
    ];
    let complex = TetComplex::from_tets(5, vec![[0, 1, 2, 3], [0, 2, 1, 4]]);
    let geometry = HelmholtzFe::new(&complex, &positions, 1, GEOMETRY_MEDIUM).unwrap();
    let media = [Medium { density_kg_m3: 2.0, sound_speed_m_s: 1.5 },
        Medium { density_kg_m3: 8.0, sound_speed_m_s: 3.0 }];
    let mut source = |element: usize, _: [f64; 3]| Ok::<_, AcousticFeError>(media[element]);
    let model = HeterogeneousHelmholtzFe::new(&geometry, 3, &mut source, &mut || false).unwrap();
    let omega = 2.0;
    // Continuous pressure, discontinuous normal gradient. Momentum requires
    // the same v_x on both sides, not the same pressure gradient.
    let coefficients = [C64::ONE, C64::ONE, C64::ONE,
        C64::new(1.0, omega * media[0].density_kg_m3),
        C64::new(1.0, -omega * media[1].density_kg_m3)];
    for (element, medium) in media.into_iter().enumerate() {
        let probe = model.probe(element, [1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0, 0.0], &mut source).unwrap();
        assert_eq!(probe.medium(), medium);
        let sample = probe.harmonic(&coefficients, omega / core::f64::consts::TAU).unwrap();
        close(sample.pressure_pa, C64::ONE);
        close(sample.pressure_gradient_pa_m[0], C64::new(0.0, -omega * medium.density_kg_m3));
        close(sample.particle_velocity_m_s[0], C64::ONE);
        assert!((sample.active_intensity_w_m2[0] - 0.5).abs() < 1e-14);
    }
}

#[test]
fn transient_probes_use_potential_velocity_and_instantaneous_flux_without_half_factor() {
    let complex = TetComplex::from_tets(4, vec![[3, 1, 0, 2]]);
    let geometry = HelmholtzFe::new(&complex, &POSITIONS, 2, GEOMETRY_MEDIUM).unwrap();
    let model = uniform(&geometry);
    let prepared = model.prepare_time(0.01, &[], HeterogeneousTimeLimits::default(), &mut || false).unwrap();
    let mut potential = vec![0.0; model.ndof()];
    let mut pressure = vec![0.0; model.ndof()];
    for (v, point) in POSITIONS.iter().enumerate() {
        potential[v] = 0.5 * point[0] * point[0];
        pressure[v] = 2.0 + point[1];
    }
    for (dof, shape) in geometry.space.element_dofs(0) {
        if let LocalFn::Edge(a, b, 0) = shape {
            let tet = geometry.tetrahedra()[0];
            let dx = POSITIONS[tet[a] as usize][0] - POSITIONS[tet[b] as usize][0];
            potential[dof] = -0.5 * dx * dx;
        }
    }
    let probe = model.probe(0, [0.25; 4], &mut |_, _| Ok::<_, AcousticFeError>(MEDIUM)).unwrap();
    let [x, y, _] = probe.position_m();
    let expected_pressure = 2.0 + y;
    let expected_v = -x / MEDIUM.density_kg_m3;
    for gauge in [0.0, 17.0] {
        let mut shifted = potential.clone();
        for v in &mut shifted[..POSITIONS.len()] {
            *v += gauge;
        }
        let state = prepared.state(&shifted, &pressure).unwrap();
        let sample = probe.transient(&state).unwrap();
        assert!((sample.pressure_pa - expected_pressure).abs() < 1e-13);
        assert!((sample.pressure_potential_pa_s - (0.5 * x * x + gauge)).abs() < 1e-13);
        assert!((sample.potential_gradient_pa_s_m[0] - x).abs() < 1e-13);
        assert!((sample.particle_velocity_m_s[0] - expected_v).abs() < 1e-13);
        assert!((sample.instantaneous_intensity_w_m2[0] - expected_pressure * expected_v).abs() < 1e-13);
        assert!(sample.particle_velocity_m_s[1..].iter().all(|v| v.abs() < 1e-13));
    }
    let other_model = uniform(&geometry);
    let other_time = other_model.prepare_time(0.01, &[], HeterogeneousTimeLimits::default(), &mut || false).unwrap();
    assert!(probe.transient(&other_time.zero_state()).is_err());
}

#[test]
fn probes_preserve_source_errors_and_refuse_invalid_material_location_or_fields() {
    #[derive(Debug, PartialEq)]
    enum SourceError { Original, Acoustic(AcousticFeError) }
    impl From<AcousticFeError> for SourceError {
        fn from(value: AcousticFeError) -> Self { Self::Acoustic(value) }
    }
    let complex = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let geometry = HelmholtzFe::new(&complex, &POSITIONS, 1, GEOMETRY_MEDIUM).unwrap();
    let model = uniform(&geometry);
    assert!(matches!(model.probe(0, [0.25; 4], &mut |_, _| Err::<Medium, _>(SourceError::Original)), Err(SourceError::Original)));
    for (element, point) in [(1, [0.25; 4]), (0, [-0.1, 0.4, 0.4, 0.3]), (0, [0.2; 4]), (0, [f64::NAN; 4])] {
        let mut calls = 0;
        assert!(model.probe(element, point, &mut |_, _| { calls += 1; Ok::<_, AcousticFeError>(MEDIUM) }).is_err());
        assert_eq!(calls, 0);
    }
    for medium in [
        Medium { density_kg_m3: 0.0, ..MEDIUM },
        Medium { density_kg_m3: f64::NAN, ..MEDIUM },
        Medium { sound_speed_m_s: -1.0, ..MEDIUM },
        Medium { sound_speed_m_s: f64::INFINITY, ..MEDIUM },
        Medium { sound_speed_m_s: f64::MAX, ..MEDIUM },
    ] {
        assert!(model.probe(0, [0.25; 4], &mut |_, _| Ok::<_, AcousticFeError>(medium)).is_err());
    }
    let probe = model.probe(0, [0.25; 4], &mut |_, _| Ok::<_, AcousticFeError>(MEDIUM)).unwrap();
    for frequency in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::MAX] {
        assert!(probe.harmonic(&[C64::ONE; 4], frequency).is_err());
    }
    assert!(probe.harmonic(&[C64::ONE; 3], 1.0).is_err());
    assert!(probe.harmonic(&[C64::new(f64::NAN, 0.0); 4], 1.0).is_err());
    assert!(probe.harmonic(&[C64::new(f64::MAX, f64::MAX), C64::ZERO, C64::ZERO, C64::ZERO], 1.0).is_err());
}
