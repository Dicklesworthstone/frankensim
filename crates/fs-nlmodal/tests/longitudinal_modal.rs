//! Physical continuum and coupled-harmonic checks for the complete modal owner.
use fs_nlmodal::{
    KcStringParams, kirchhoff_carrier_string, longitudinal::LongitudinalModalStorage,
};
use fs_phs::{ConstitutiveDiscreteGradient, ContinuousStorageTangent, Storage};

fn parameters() -> KcStringParams {
    KcStringParams {
        length: 1.3,
        tension: 17.0,
        lin_density: 0.7,
        ea: 120.0,
    }
}

fn close(actual: f64, expected: f64, scale: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * scale,
        "actual={actual:.17e}, expected={expected:.17e}, scale={scale:.17e}"
    );
}

fn fixture(storage: &LongitudinalModalStorage, phase: f64) -> Vec<f64> {
    let transverse = storage.plane_count() * storage.transverse_modes_per_plane();
    (0..2 * storage.mode_count())
        .map(|i| {
            let magnitude = if i % 2 == 1 {
                0.03
            } else if i / 2 < transverse {
                0.004
            } else {
                0.0003
            };
            magnitude * ((i + 1) as f64 * 0.7 + phase).sin()
        })
        .collect()
}

// Independent physical-space quadrature. Omitted axial strain harmonics are
// eliminated by their spatial L2 minimizer, not by the production cosine formulas.
fn continuum(
    p: KcStringParams,
    n: usize,
    bending: &[f64],
    axial_modes: &[usize],
    state: &[f64],
    samples: usize,
) -> (f64, Vec<f64>) {
    let normalization = (2.0 / (p.lin_density * p.length)).sqrt();
    let count = n * bending.len() + axial_modes.len();
    let mut geometric = vec![0.0; samples];
    for (i, value) in geometric.iter_mut().enumerate() {
        let fraction = (i as f64 + 0.5) / samples as f64;
        for plane in 0..bending.len() {
            let slope: f64 = (1..=n)
                .map(|j| {
                    let k = j as f64 * std::f64::consts::PI / p.length;
                    normalization
                        * k
                        * state[2 * (plane * n + j - 1)]
                        * (k * p.length * fraction).cos()
                })
                .sum();
            *value += 0.5 * slope * slope;
        }
    }
    let omitted: Vec<(usize, f64)> = (1..=2 * n)
        .filter(|j| !axial_modes.contains(j))
        .map(|j| {
            (
                j,
                2.0 * geometric
                    .iter()
                    .enumerate()
                    .map(|(i, s)| {
                        s * (j as f64 * std::f64::consts::PI * (i as f64 + 0.5) / samples as f64)
                            .cos()
                    })
                    .sum::<f64>()
                    / samples as f64,
            )
        })
        .collect();
    let mut energy = 0.0;
    let mut gradient = vec![0.0; 2 * count];
    for mode in 0..count {
        energy += 0.5 * state[2 * mode + 1] * state[2 * mode + 1];
        gradient[2 * mode + 1] = state[2 * mode + 1];
    }
    for (i, &geometric) in geometric.iter().enumerate() {
        let fraction = (i as f64 + 0.5) / samples as f64;
        let mut strain = geometric;
        for (slot, &j) in axial_modes.iter().enumerate() {
            let k = j as f64 * std::f64::consts::PI / p.length;
            strain += normalization
                * k
                * state[2 * (n * bending.len() + slot)]
                * (k * p.length * fraction).cos();
        }
        for &(j, s) in &omitted {
            strain -= s * (j as f64 * std::f64::consts::PI * fraction).cos();
        }
        energy += 0.5 * p.ea * strain * strain * p.length / samples as f64;
        for (plane, &ei) in bending.iter().enumerate() {
            let mut slope = 0.0;
            let mut curvature = 0.0;
            for j in 1..=n {
                let k = j as f64 * std::f64::consts::PI / p.length;
                let q = state[2 * (plane * n + j - 1)];
                slope += normalization * k * q * (k * p.length * fraction).cos();
                curvature -= normalization * k * k * q * (k * p.length * fraction).sin();
            }
            energy += 0.5 * (p.tension * slope * slope + ei * curvature * curvature) * p.length
                / samples as f64;
            for j in 1..=n {
                let k = j as f64 * std::f64::consts::PI / p.length;
                let dphi = normalization * k * (k * p.length * fraction).cos();
                let ddphi = -normalization * k * k * (k * p.length * fraction).sin();
                gradient[2 * (plane * n + j - 1)] += ((p.tension + p.ea * strain) * slope * dphi
                    + ei * curvature * ddphi)
                    * p.length
                    / samples as f64;
            }
        }
        for (slot, &j) in axial_modes.iter().enumerate() {
            let k = j as f64 * std::f64::consts::PI / p.length;
            gradient[2 * (n * bending.len() + slot)] +=
                p.ea * strain * normalization * k * (k * p.length * fraction).cos() * p.length
                    / samples as f64;
        }
    }
    (energy, gradient)
}

#[test]
fn full_physical_energy_and_force_match_resolved_spatial_quadrature() {
    let p = parameters();
    let bending = [0.004, 0.012];
    for modes in [vec![], vec![2, 5], (1..=6).collect()] {
        let source = LongitudinalModalStorage::new(p, 3, &bending, &modes).unwrap();
        let state = fixture(&source, 0.2);
        let mut force = vec![0.0; state.len()];
        source.gradient(&state, &mut force);
        for samples in [64, 128, 256] {
            // The largest product order here is 12, so every rung resolves
            // spatial strain-energy products, including selected mode five/six.
            let (energy, expected) = continuum(p, 3, &bending, &modes, &state, samples);
            close(source.hamiltonian(&state), energy, energy, 3e-13);
            for (a, b) in force.iter().zip(&expected) {
                close(*a, *b, b.abs().max(1e-6), 4e-12);
            }
        }
        let (aliased, _) = continuum(p, 3, &bending, &modes, &state, 3);
        assert!((source.hamiltonian(&state) - aliased).abs() > 1e-5);
    }
}

#[test]
fn empty_axial_selection_is_the_original_kirchhoff_carrier_law() {
    let p = parameters();
    let old = kirchhoff_carrier_string(&p, 4).unwrap();
    let storage = LongitudinalModalStorage::new(p, 4, &[0.0], &[]).unwrap();
    let state = fixture(&storage, 0.3);
    close(
        storage.hamiltonian(&state),
        old.hamiltonian(&state),
        old.hamiltonian(&state),
        3e-14,
    );
    let mut a = vec![0.0; state.len()];
    let mut b = a.clone();
    storage.gradient(&state, &mut a);
    old.gradient(&state, &mut b);
    for (a, b) in a.into_iter().zip(b) {
        close(a, b, b.abs().max(1e-8), 3e-14);
    }
}

#[test]
fn continuous_hessian_and_finite_endpoint_jacobian_are_original_force_derivatives() {
    let storage =
        LongitudinalModalStorage::new(parameters(), 3, &[0.004, 0.012], &[1, 2, 5, 6]).unwrap();
    let a = fixture(&storage, 0.1);
    let b = fixture(&storage, 0.8);
    let direction = fixture(&storage, 1.9);
    let other_direction = fixture(&storage, 2.7);
    let mut action = vec![0.0; b.len()];
    let mut other_action = action.clone();
    let boxed: Box<dyn Storage> = Box::new(storage.clone());
    let tangent = boxed.continuous_tangent().unwrap();
    assert!(tangent.apply(&b, &direction, &mut action));
    assert!(tangent.apply(&b, &other_direction, &mut other_action));
    let work_a: f64 = action
        .iter()
        .zip(&other_direction)
        .map(|(a, d)| a * d)
        .sum();
    let work_b: f64 = other_action
        .iter()
        .zip(&direction)
        .map(|(a, d)| a * d)
        .sum();
    close(work_a, work_b, work_a.abs().max(work_b.abs()), 5e-14);
    let epsilon = 0.125;
    let perturbed = |factor: f64| {
        b.iter()
            .zip(&direction)
            .map(|(b, d)| b + factor * epsilon * d)
            .collect::<Vec<_>>()
    };
    let four: Vec<_> = [-2.0, -1.0, 1.0, 2.0]
        .iter()
        .map(|factor| {
            let mut out = vec![0.0; b.len()];
            boxed.gradient(&perturbed(*factor), &mut out);
            out
        })
        .collect();
    // This five-point derivative is exact for the coded cubic force polynomial;
    // epsilon need not be a cancellation-prone tiny perturbation.
    for j in 0..b.len() {
        let numerical =
            (four[0][j] - 8.0 * four[1][j] + 8.0 * four[2][j] - four[3][j]) / (12.0 * epsilon);
        close(action[j], numerical, numerical.abs().max(1e-7), 3e-12);
    }
    let rule = boxed.discrete_gradient_rule().unwrap();
    assert!(rule.end_jacobian_action(&a, &b, &direction, &mut action));
    let four: Vec<_> = [-2.0, -1.0, 1.0, 2.0]
        .iter()
        .map(|factor| {
            let end = perturbed(*factor);
            let mid: Vec<_> = a
                .iter()
                .zip(&end)
                .map(|(a, b)| f64::midpoint(*a, *b))
                .collect();
            let mut out = vec![0.0; b.len()];
            rule.evaluate(&a, &end, &mid, &mut out);
            out
        })
        .collect();
    for j in 0..b.len() {
        let numerical =
            (four[0][j] - 8.0 * four[1][j] + 8.0 * four[2][j] - four[3][j]) / (12.0 * epsilon);
        close(action[j], numerical, numerical.abs().max(1e-7), 4e-12);
    }
    assert!(rule.end_jacobian_action(&b, &b, &direction, &mut action));
    assert!(tangent.apply(&b, &direction, &mut other_action));
    for (a, b) in action.iter().zip(&other_action) {
        close(2.0 * a, *b, b.abs().max(1e-8), 3e-14);
    }
    let mut transverse_direction = vec![0.0; b.len()];
    transverse_direction[0] = 1.0;
    assert!(tangent.apply(&b, &transverse_direction, &mut action));
    assert!(
        action[2 * storage.transverse_modes_per_plane()..]
            .iter()
            .step_by(2)
            .any(|x| x.abs() > 1e-5)
    );
}

#[test]
fn source_discrete_gradient_preserves_complete_kinetic_and_elastic_work() {
    let storage =
        LongitudinalModalStorage::new(parameters(), 4, &[0.003, 0.009], &[1, 2, 4, 7, 8]).unwrap();
    let a = fixture(&storage, 0.4);
    let b = fixture(&storage, 1.2);
    let mid: Vec<_> = a
        .iter()
        .zip(&b)
        .map(|(a, b)| f64::midpoint(*a, *b))
        .collect();
    let mut effort = vec![0.0; a.len()];
    ConstitutiveDiscreteGradient::evaluate(&storage, &a, &b, &mid, &mut effort);
    let work: f64 = effort
        .iter()
        .zip(a.iter().zip(&b))
        .map(|(e, (a, b))| e * (b - a))
        .sum();
    let difference = storage.hamiltonian(&b) - storage.hamiltonian(&a);
    close(
        work,
        difference,
        storage.hamiltonian(&a) + storage.hamiltonian(&b),
        4e-15,
    );
    let mut ordinary = vec![0.0; a.len()];
    storage.gradient(&b, &mut ordinary);
    ConstitutiveDiscreteGradient::evaluate(&storage, &b, &b, &b, &mut effort);
    for (ordinary, effort) in ordinary.iter().zip(&effort) {
        close(*ordinary, *effort, effort.abs().max(1e-8), 3e-14);
    }
}

#[test]
fn circular_polarization_is_an_exact_coupled_orbit_with_axial_dc() {
    let storage = LongitudinalModalStorage::new(parameters(), 1, &[0.006, 0.006], &[2]).unwrap();
    let amplitude = 0.06;
    let d = storage.transverse_slope_scales()[0];
    let c = storage.source().strain_energy_coefficient();
    let mean_strain = 0.25 * d * d * amplitude * amplitude;
    let axial_dc = -mean_strain / storage.source().axial_slope_scales()[0];
    let omega2 = storage.omegas()[0] * storage.omegas()[0] + c * mean_strain * d * d;
    let omega = omega2.sqrt();
    let mut reference_energy = 0.0;
    for j in 0..32 {
        let phase = std::f64::consts::TAU * j as f64 / 32.0;
        let state = [
            amplitude * phase.cos(),
            -amplitude * omega * phase.sin(),
            amplitude * phase.sin(),
            amplitude * omega * phase.cos(),
            axial_dc,
            0.0,
        ];
        let mut force = [0.0; 6];
        storage.gradient(&state, &mut force);
        close(force[0], omega2 * state[0], omega2 * amplitude, 3e-14);
        close(force[2], omega2 * state[2], omega2 * amplitude, 3e-14);
        close(force[4], 0.0, omega2 * amplitude, 3e-14);
        assert_eq!(force[1], state[1]);
        assert_eq!(force[3], state[3]);
        if j == 0 {
            reference_energy = storage.hamiltonian(&state);
        }
        close(
            storage.hamiltonian(&state),
            reference_energy,
            reference_energy,
            3e-14,
        );
    }
    assert!(axial_dc < 0.0);
    assert!(omega > storage.omegas()[0]);
}

#[test]
fn sum_frequency_drives_the_true_axial_resonance_with_fourth_order_time_refinement() {
    let mut p = parameters();
    p.ea = 121.0 * p.tension;
    let storage = LongitudinalModalStorage::new(p, 6, &[0.0], &[1]).unwrap();
    let w = storage.omegas()[0];
    let axial_w = storage.omegas()[6];
    close(axial_w, 11.0 * w, axial_w, 3e-15);
    let (a5, a6) = (3e-4, 2e-4);
    let scales = storage.transverse_slope_scales();
    let forcing = -0.25
        * storage.source().strain_energy_coefficient()
        * storage.source().axial_slope_scales()[0]
        * scales[4]
        * scales[5]
        * a5
        * a6;
    let analytic = |t: f64| {
        let nonresonant = forcing / (axial_w * axial_w - w * w);
        let resonant = forcing / (2.0 * axial_w);
        [
            nonresonant * ((w * t).cos() - (axial_w * t).cos())
                + resonant * t * (axial_w * t).sin(),
            nonresonant * (-w * (w * t).sin() + axial_w * (axial_w * t).sin())
                + resonant * ((axial_w * t).sin() + axial_w * t * (axial_w * t).cos()),
        ]
    };
    // Two externally prescribed transverse modes isolate a physical axial
    // response oracle. This is not a one-way approximation of a free string:
    // production gradient retains the reciprocal transverse reaction.
    let rhs = |t: f64, u: [f64; 2]| {
        let mut state = [0.0; 14];
        state[8] = a5 * (5.0 * w * t).cos();
        state[10] = a6 * (6.0 * w * t).cos();
        state[12] = u[0];
        state[13] = u[1];
        let mut gradient = [0.0; 14];
        storage.gradient(&state, &mut gradient);
        [u[1], -gradient[12]]
    };
    let duration = std::f64::consts::TAU / w;
    let scale = (forcing * duration / (2.0 * axial_w)).abs();
    let mut errors = Vec::new();
    for steps in [1024, 2048, 4096] {
        let dt = duration / steps as f64;
        let mut u = [0.0; 2];
        let mut error = 0.0_f64;
        for step in 0..steps {
            let time = step as f64 * dt;
            let a = rhs(time, u);
            let b = rhs(
                time + 0.5 * dt,
                [u[0] + 0.5 * dt * a[0], u[1] + 0.5 * dt * a[1]],
            );
            let c = rhs(
                time + 0.5 * dt,
                [u[0] + 0.5 * dt * b[0], u[1] + 0.5 * dt * b[1]],
            );
            let d = rhs(time + dt, [u[0] + dt * c[0], u[1] + dt * c[1]]);
            for j in 0..2 {
                u[j] += dt / 6.0 * (a[j] + 2.0 * b[j] + 2.0 * c[j] + d[j]);
            }
            let exact = analytic(time + dt);
            error = error
                .max((u[0] - exact[0]).abs() / scale)
                .max((u[1] - exact[1]).abs() / (scale * axial_w));
        }
        errors.push(error);
        assert!(u[1].abs() > 0.99 * (forcing * duration / 2.0).abs());
    }
    assert!(
        errors[0] / errors[1] > 14.0 && errors[1] / errors[2] > 14.0,
        "{errors:?}"
    );
    assert!(errors[2] < 1e-7, "{errors:?}");
}

#[test]
fn original_linear_metadata_state_order_and_bounded_admission_survive_boxing() {
    let p = parameters();
    let source = LongitudinalModalStorage::new(p, 2, &[0.004, 0.012], &[1, 8]).unwrap();
    assert_eq!(source.mode_count(), 6);
    assert_eq!(source.source().axial_modes(), &[1, 8]);
    let storage: Box<dyn Storage> = Box::new(source.clone());
    let q = storage.quadratic_part().unwrap();
    assert!(!q.exact);
    let direction: Vec<_> = (0..12).map(|i| 0.01 * (i + 1) as f64).collect();
    let mut action = [0.0; 12];
    assert!(
        storage
            .continuous_tangent()
            .unwrap()
            .apply(&[0.0; 12], &direction, &mut action)
    );
    for i in 0..12 {
        let expected: f64 = q.matrix[i * 12..(i + 1) * 12]
            .iter()
            .zip(&direction)
            .map(|(a, b)| a * b)
            .sum();
        close(action[i], expected, expected.abs().max(1e-8), 4e-15);
    }
    let clone = storage.shared_clone().unwrap();
    assert_eq!(
        storage.hamiltonian(&direction),
        clone.hamiltonian(&direction)
    );
    let axial_only = LongitudinalModalStorage::new(p, 0, &[0.0], &[1, 4]).unwrap();
    assert!(axial_only.quadratic_part().unwrap().exact);
    assert!(LongitudinalModalStorage::new(p, 0, &[0.0], &[]).is_err());
    assert!(LongitudinalModalStorage::new(p, 2, &[-0.1], &[1]).is_err());
    assert!(LongitudinalModalStorage::new(p, 2, &[0.0, 0.0, 0.0], &[1]).is_err());
    assert!(LongitudinalModalStorage::new(p, 513, &[0.0], &[1]).is_err());
    assert!(LongitudinalModalStorage::new(p, 2, &[0.0], &[2, 1]).is_err());
    assert!(!ContinuousStorageTangent::apply(
        &source,
        &[0.0; 11],
        &direction,
        &mut action
    ));
    assert!(!ContinuousStorageTangent::apply(
        &source,
        &direction,
        &[f64::NAN; 12],
        &mut action
    ));
    source.gradient(&[0.0; 11], &mut action);
    assert!(action.iter().all(|x| x.is_nan()));
}
