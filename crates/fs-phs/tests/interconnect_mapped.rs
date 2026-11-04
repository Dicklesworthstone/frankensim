//! Physical references use elementary coupled modes and a separate nonlinear
//! pressure ODE, not the composed structure matrix or its prepared exponential.
use fs_phs::linear_exact::{
    DenseExponentialBudget, LinearSplitBudget, PreparedExactTones, QuadraticPart,
};
use fs_phs::{
    ConstitutiveDiscreteGradient, InterconnectionBudget, PhsError, PortHamiltonian,
    QuadraticStorage, ResistanceWitness, Storage, interconnect, interconnect_mapped, step,
    transformer,
};

fn budget() -> InterconnectionBudget {
    InterconnectionBudget {
        max_states: 32,
        max_total_ports: 128,
        max_dense_bytes: 16 * 1024 * 1024,
        max_coupling_multiply_terms: 1_000_000,
    }
}
fn split() -> LinearSplitBudget {
    LinearSplitBudget {
        max_modes: 32,
        max_omega_dt: 0.1,
    }
}
fn dense() -> DenseExponentialBudget {
    DenseExponentialBudget {
        max_dimension: 64,
        max_dense_bytes: 128 * 1024 * 1024,
        max_multiply_terms: 20_000_000_000,
    }
}
fn near(actual: f64, expected: f64, relative: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * expected.abs(),
        "{actual:.17e} != {expected:.17e}"
    );
}
fn inertias(masses: &[f64], decay: f64) -> PortHamiltonian {
    let n = masses.len();
    let mut q = vec![0.0; n * n];
    let mut r = vec![0.0; n * n];
    let mut roots = vec![0.0; n * n];
    let mut g = vec![0.0; n * 2 * n];
    for i in 0..n {
        q[i * n + i] = 1.0 / masses[i];
        r[i * n + i] = decay * masses[i];
        roots[i * n + i] = r[i * n + i].sqrt();
        g[i * 2 * n + i] = 1.0;
        g[i * 2 * n + n + i] = 1.0;
    }
    let witness = ResistanceWitness::try_from_source_rows(n, r.clone(), roots).unwrap();
    PortHamiltonian::new(
        n,
        2 * n,
        vec![0.0; n * n],
        r,
        g,
        Box::new(QuadraticStorage::new(q, n).unwrap()),
    )
    .unwrap()
    .with_resistance_witness(witness)
    .unwrap()
}
fn chamber(compliance: f64, decay: f64) -> PortHamiltonian {
    let r = decay * compliance;
    PortHamiltonian::new(
        1,
        2,
        vec![0.0],
        vec![r],
        vec![1.0, 1.0],
        Box::new(QuadraticStorage::new(vec![1.0 / compliance], 1).unwrap()),
    )
    .unwrap()
    .with_resistance_witness(
        ResistanceWitness::try_from_source_rows(1, vec![r], vec![r.sqrt()]).unwrap(),
    )
    .unwrap()
}
fn simpson(f: impl Fn(f64) -> f64, duration: f64) -> f64 {
    const N: usize = 4096;
    let h = duration / N as f64;
    let mut sum = f(0.0) + f(duration);
    for i in 1..N {
        sum += if i % 2 == 0 { 2.0 } else { 4.0 } * f(i as f64 * h);
    }
    h * sum / 3.0
}
fn signal(t: f64, frequencies: &[f64], cosine: &[f64], sine: &[f64]) -> f64 {
    frequencies
        .iter()
        .zip(cosine)
        .zip(sine)
        .map(|((&w, &c), &s)| c * (w * t).cos() + s * (w * t).sin())
        .sum()
}

const MASSES: [f64; 2] = [2.0, 3.0];
const AREAS: [f64; 2] = [0.1, -0.07];
const COMPLIANCE: f64 = 0.004;
const INITIAL: [f64; 3] = [0.8, -0.3, 0.003];

// p_i = A_i*z + n_i*exp(-decay*t), where sum A_i*n_i/m_i = 0.
// The physical pressure mode satisfies z'=-decay*z-P+f(t),
// P'=-decay*P+omega^2*z. Solve this two-variable oscillator analytically.
fn piston_motion(t: f64, decay: f64, frequencies: &[f64], c: &[f64], s: &[f64]) -> [f64; 3] {
    let stiffness: f64 = (0..2).map(|i| AREAS[i] * AREAS[i] / MASSES[i]).sum();
    let omega2 = stiffness / COMPLIANCE;
    let omega = omega2.sqrt();
    let z0: f64 = (0..2)
        .map(|i| AREAS[i] * INITIAL[i] / MASSES[i])
        .sum::<f64>()
        / stiffness;
    let particular = |t: f64| {
        let mut z = 0.0;
        let mut pressure = 0.0;
        for ((&w, &c), &s) in frequencies.iter().zip(c).zip(s) {
            let re = decay * decay - w * w + omega2;
            let im = 2.0 * decay * w;
            let denominator = re * re + im * im;
            let fr = (c * re - s * im) / denominator;
            let fi = (-s * re - c * im) / denominator;
            let zr = decay * fr - w * fi;
            let zi = decay * fi + w * fr;
            z += zr * (w * t).cos() - zi * (w * t).sin();
            pressure += omega2 * (fr * (w * t).cos() - fi * (w * t).sin());
        }
        [z, pressure]
    };
    let start = particular(0.0);
    let part = particular(t);
    let a = z0 - start[0];
    let b = INITIAL[2] / COMPLIANCE - start[1];
    let attenuation = (-decay * t).exp();
    let z = attenuation * (a * (omega * t).cos() - b / omega * (omega * t).sin()) + part[0];
    let pressure = attenuation * (b * (omega * t).cos() + omega * a * (omega * t).sin()) + part[1];
    [
        AREAS[0] * z + (INITIAL[0] - AREAS[0] * z0) * attenuation,
        AREAS[1] * z + (INITIAL[1] - AREAS[1] * z0) * attenuation,
        COMPLIANCE * pressure,
    ]
}

#[test]
fn continuous_forces_drive_the_complete_shared_pressure_mode_and_original_heat() {
    let omega = ((AREAS[0].powi(2) / MASSES[0] + AREAS[1].powi(2) / MASSES[1]) / COMPLIANCE).sqrt();
    let frequencies = [0.0, 0.61 * omega, 1.73 * omega];
    let c = [-0.4, 0.2, -0.05];
    let s = [0.0, -0.1, 0.08];
    let mut cosine = vec![0.0; 9];
    let mut sine = vec![0.0; 9];
    for i in 0..2 {
        for k in 0..3 {
            cosine[3 * i + k] = AREAS[i] * c[k];
            sine[3 * i + k] = AREAS[i] * s[k];
        }
    }
    let mut regenerative = false;
    for decay in [0.0, 0.35] {
        for duration in [0.17, 0.6, 1.3] {
            let system = interconnect_mapped(
                inertias(&MASSES, decay),
                chamber(COMPLIANCE, decay),
                &[0, 1],
                &[0],
                &AREAS,
                budget(),
            )
            .unwrap();
            assert_eq!(system.state_dim(), 3);
            assert_eq!(system.port_dim(), 3);
            assert!(system.resistance_witness().is_some());
            let mut prepared =
                PreparedExactTones::new(&system, &frequencies, duration, split(), dense()).unwrap();
            let mut endpoint = [0.0; 3];
            let mut output = [0.0; 3];
            let record = prepared
                .step_into(&INITIAL, &cosine, &sine, &mut endpoint, &mut output)
                .unwrap();
            let expected = piston_motion(duration, decay, &frequencies, &c, &s);
            for i in 0..3 {
                near(endpoint[i], expected[i], 3e-11, 2e-12);
            }
            near(
                AREAS[1] * endpoint[0] - AREAS[0] * endpoint[1],
                (-decay * duration).exp() * (AREAS[1] * INITIAL[0] - AREAS[0] * INITIAL[1]),
                2e-11,
                2e-13,
            );
            let motion = |t| piston_motion(t, decay, &frequencies, &c, &s);
            let work = simpson(
                |t| {
                    let x = motion(t);
                    signal(t, &frequencies, &c, &s)
                        * (AREAS[0] * x[0] / MASSES[0] + AREAS[1] * x[1] / MASSES[1])
                },
                duration,
            );
            let heat = simpson(
                |t| {
                    let x = motion(t);
                    decay
                        * (x[0] * x[0] / MASSES[0]
                            + x[1] * x[1] / MASSES[1]
                            + x[2] * x[2] / COMPLIANCE)
                },
                duration,
            );
            near(record.supplied, work, 5e-11, 5e-13);
            near(record.dissipated, heat, 5e-11, 5e-13);
            near(
                record.delta_h + record.dissipated,
                record.supplied,
                3e-11,
                2e-12,
            );
            for i in 0..3 {
                let scale = if i < 2 { MASSES[i] } else { COMPLIANCE };
                near(
                    output[i],
                    simpson(|t| motion(t)[i] / scale, duration) / duration,
                    4e-11,
                    2e-12,
                );
            }
            regenerative |= record.supplied < 0.0;
            assert!(record.dissipated >= 0.0);
        }
    }
    assert!(regenerative);
}

fn matrix_system(n: usize, m: usize, g: Vec<f64>) -> PortHamiltonian {
    let mut q = vec![0.0; n * n];
    for i in 0..n {
        q[i * n + i] = 1.0 + i as f64;
    }
    PortHamiltonian::new(
        n,
        m,
        vec![0.0; n * n],
        vec![0.0; n * n],
        g,
        Box::new(QuadraticStorage::new(q, n).unwrap()),
    )
    .unwrap()
}

#[test]
fn rectangular_contractions_preserve_port_order_and_reverse_with_minus_transpose() {
    let ga = vec![1.0, 0.3, -0.2, 0.5, -0.4, 1.1];
    let gb = vec![
        -0.2, 1.2, 0.4, 0.0, 1.0, -0.3, -0.6, 0.8, 0.7, 0.2, 1.1, -0.5,
    ];
    let left = [0, 2];
    let right = [1, 3, 0];
    let map = [0.5, -1.2, 0.0, 0.3, 2.0, 0.8];
    let compose = |left: &[usize], right: &[usize], map: &[f64]| {
        interconnect_mapped(
            matrix_system(2, 3, ga.clone()),
            matrix_system(3, 4, gb.clone()),
            left,
            right,
            map,
            budget(),
        )
        .unwrap()
    };
    let original = compose(&left, &right, &map);
    let permuted = compose(&[2, 0], &[0, 1, 3], &[0.8, 2.0, -1.2, 0.5, 0.3, 0.0]);
    let (j, _, g) = original.structure();
    assert_eq!(g, &[0.3, 0.0, -0.4, 0.0, 0.0, 0.4, 0.0, -0.6, 0.0, 1.1]);
    for i in 0..2 {
        for k in 0..3 {
            let mut physical = 0.0;
            for r in 0..3 {
                for l in 0..2 {
                    physical += ga[i * 3 + left[l]] * map[r * 2 + l] * gb[k * 4 + right[r]];
                }
            }
            near(j[i * 5 + 2 + k], -physical, 2e-14, 1e-15);
            assert_eq!(j[i * 5 + 2 + k], -j[(2 + k) * 5 + i]);
        }
    }
    for (a, b) in j.iter().zip(permuted.structure().0) {
        near(*a, *b, 2e-14, 1e-15);
    }
    assert_eq!(original.structure().2, permuted.structure().2);
    let reversed = interconnect_mapped(
        matrix_system(3, 4, gb),
        matrix_system(2, 3, ga),
        &right,
        &left,
        &[-0.5, 0.0, -2.0, 1.2, -0.3, -0.8],
        budget(),
    )
    .unwrap();
    let permutation = [2, 3, 4, 0, 1];
    for i in 0..5 {
        for k in 0..5 {
            near(
                reversed.structure().0[i * 5 + k],
                j[permutation[i] * 5 + permutation[k]],
                2e-14,
                1e-15,
            );
        }
        for port in 0..2 {
            assert_eq!(
                reversed.structure().2[i * 2 + port],
                g[permutation[i] * 2 + 1 - port]
            );
        }
    }
}

#[test]
fn canonical_and_transformer_entries_retain_their_physical_signs_and_external_ports() {
    let m = 2.0;
    let compliance = 0.004;
    for mapping in [1.0_f64, -0.125, 0.0] {
        let system = if mapping == 1.0 {
            interconnect(inertias(&[m], 0.0), chamber(compliance, 0.0), &[(0, 0)]).unwrap()
        } else {
            transformer(
                inertias(&[m], 0.0),
                chamber(compliance, 0.0),
                0,
                0,
                -mapping,
            )
            .unwrap()
        };
        assert_eq!(system.port_dim(), 2);
        assert_eq!(system.structure().2, &[1.0, 0.0, 0.0, 1.0]);
        let mapped = interconnect_mapped(
            inertias(&[m], 0.0),
            chamber(compliance, 0.0),
            &[0],
            &[0],
            &[mapping],
            budget(),
        )
        .unwrap();
        assert_eq!(system.structure(), mapped.structure());
        let duration = 0.13;
        let mut prepared =
            PreparedExactTones::new(&system, &[0.0], duration, split(), dense()).unwrap();
        let mut endpoint = [0.0; 2];
        let mut output = [0.0; 2];
        let record = prepared
            .step_into(
                &[1.0, 0.0],
                &[0.0; 2],
                &[0.0; 2],
                &mut endpoint,
                &mut output,
            )
            .unwrap();
        if mapping == 0.0 {
            assert_eq!(endpoint, [1.0, 0.0]);
        } else {
            let omega = mapping.abs() / (m * compliance).sqrt();
            near(endpoint[0], (omega * duration).cos(), 2e-12, 2e-13);
            near(
                endpoint[1],
                mapping * (omega * duration).sin() / (m * omega),
                2e-12,
                2e-13,
            );
        }
        assert_eq!(record.dissipated, 0.0);
    }
    let direct = interconnect(inertias(&[m], 0.0), chamber(compliance, 0.0), &[]).unwrap();
    assert_eq!(direct.port_dim(), 4);
    assert_eq!(
        direct.structure().2,
        &[1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0]
    );
    let mapped = interconnect_mapped(
        inertias(&[m], 0.0),
        chamber(compliance, 0.0),
        &[],
        &[],
        &[],
        budget(),
    )
    .unwrap();
    assert_eq!(direct.structure(), mapped.structure());
}

struct NonlinearPressure;
impl Storage for NonlinearPressure {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        x[0] * x[0] + 10.0 * x[0].powi(4)
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        out[0] = 2.0 * x[0] + 40.0 * x[0].powi(3);
    }
    fn quadratic_part(&self) -> Option<QuadraticPart> {
        Some(QuadraticPart {
            matrix: vec![2.0],
            exact: false,
        })
    }
    fn discrete_gradient_rule(&self) -> Option<&dyn ConstitutiveDiscreteGradient> {
        Some(self)
    }
}
impl ConstitutiveDiscreteGradient for NonlinearPressure {
    fn evaluate(&self, a: &[f64], b: &[f64], _midpoint: &[f64], out: &mut [f64]) {
        out[0] = a[0] + b[0] + 10.0 * (a[0] + b[0]) * (a[0] * a[0] + b[0] * b[0]);
    }
}
fn nonlinear_rhs(x: [f64; 5]) -> [f64; 5] {
    let pressure = 2.0 * x[2] + 40.0 * x[2].powi(3);
    let velocities = [x[0] / 2.0, x[1] / 3.0];
    [
        0.1 - 0.4 * pressure - 0.1 * x[0],
        -0.05 + 0.3 * pressure - 0.1 * x[1],
        0.02 + 0.4 * velocities[0] - 0.3 * velocities[1] - 0.02 * pressure,
        0.1 * velocities[0] - 0.05 * velocities[1] + 0.02 * pressure,
        0.1 * (x[0] * velocities[0] + x[1] * velocities[1]) + 0.02 * pressure * pressure,
    ]
}
fn nonlinear_reference(duration: f64) -> [f64; 5] {
    let mut x = [0.8, -0.3, 0.2, 0.0, 0.0];
    let h = duration / 4096.0;
    for _ in 0..4096 {
        let a = nonlinear_rhs(x);
        let b = nonlinear_rhs(std::array::from_fn(|i| x[i] + 0.5 * h * a[i]));
        let c = nonlinear_rhs(std::array::from_fn(|i| x[i] + 0.5 * h * b[i]));
        let d = nonlinear_rhs(std::array::from_fn(|i| x[i] + h * c[i]));
        for i in 0..5 {
            x[i] += h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
    x
}

#[test]
fn nonlinear_storage_and_original_loss_survive_complete_composition_and_converge() {
    let make = || {
        interconnect_mapped(
            inertias(&MASSES, 0.1),
            PortHamiltonian::new(
                1,
                2,
                vec![0.0],
                vec![0.02],
                vec![1.0, 1.0],
                Box::new(NonlinearPressure),
            )
            .unwrap(),
            &[0, 1],
            &[0],
            &[0.4, -0.3],
            budget(),
        )
        .unwrap()
    };
    let system = make();
    assert!(system.storage().discrete_gradient_rule().is_some());
    assert!(!system.storage().quadratic_part().unwrap().exact);
    assert!(PreparedExactTones::new(&system, &[1.0], 0.1, split(), dense()).is_err());
    let duration = 0.4;
    let truth = nonlinear_reference(duration);
    let mut errors = Vec::new();
    for steps in [20, 40, 80] {
        let mut x = vec![0.8, -0.3, 0.2];
        let mut work = 0.0;
        let mut heat = 0.0;
        for _ in 0..steps {
            let r = step(&system, &x, &[0.1, -0.05, 0.02], duration / steps as f64).unwrap();
            near(r.balance_residual(), 0.0, 0.0, 1e-12);
            assert!(r.dissipated >= 0.0);
            work += r.supplied;
            heat += r.dissipated;
            x = r.x;
        }
        x.extend([work, heat]);
        errors.push(
            x.iter()
                .zip(truth)
                .map(|(x, truth)| (x - truth).powi(2))
                .sum::<f64>()
                .sqrt(),
        );
    }
    assert!(
        errors[0] / errors[1] > 3.8 && errors[1] / errors[2] > 3.8,
        "{errors:?}"
    );
}

#[test]
fn sparse_wide_maps_keep_every_port_and_nonzero_coefficient_without_a_rank_cutoff() {
    let ports: Vec<_> = (0..40).collect();
    let mut map = vec![0.0; 40];
    map[39] = 1e-200;
    let make = || {
        (
            matrix_system(1, 41, vec![1.0; 41]),
            matrix_system(1, 2, vec![1.0; 2]),
        )
    };
    let bounded = InterconnectionBudget {
        max_coupling_multiply_terms: 6,
        ..budget()
    };
    let (a, b) = make();
    let system = interconnect_mapped(a, b, &ports, &[0], &map, bounded).unwrap();
    assert_eq!(system.state_dim(), 2);
    assert_eq!(system.port_dim(), 2);
    assert_eq!(system.structure().0, &[0.0, -1e-200, 1e-200, 0.0]);
    let (a, b) = make();
    assert!(
        interconnect_mapped(
            a,
            b,
            &ports,
            &[0],
            &map,
            InterconnectionBudget {
                max_coupling_multiply_terms: 5,
                ..bounded
            }
        )
        .is_err()
    );
    map.fill(0.0);
    let (a, b) = make();
    let zero = interconnect_mapped(
        a,
        b,
        &ports,
        &[0],
        &map,
        InterconnectionBudget {
            max_coupling_multiply_terms: 0,
            ..bounded
        },
    )
    .unwrap();
    assert_eq!(zero.structure().0, &[0.0; 4]);
    assert_eq!(zero.state_dim(), 2);
}

#[test]
fn intermediate_compensation_preserves_a_complete_unit_frequency_mode() {
    // Both contractions must retain the 1 after a 1e16 + 1 dot product.
    // Rounding intermediate dots to f64 would turn this physical mode into
    // zero coupling. The second case makes G_right*M the cheaper contraction.
    for left_states in [1, 2] {
        let ga = if left_states == 1 {
            vec![1.0, 1.0, 1.0]
        } else {
            vec![1.0, 1.0, 1.0, 0.0, 0.0, 1.0]
        };
        let map = if left_states == 1 {
            [1e16, 1.0, -1e16, 0.0]
        } else {
            [1e16, -1e16, 1.0, 0.0]
        };
        let system = interconnect_mapped(
            matrix_system(left_states, 3, ga),
            matrix_system(1, 3, vec![1.0; 3]),
            &[0, 1],
            &[0, 1],
            &map,
            budget(),
        )
        .unwrap();
        let n = left_states + 1;
        assert_eq!(system.structure().0[left_states], -1.0);
        assert_eq!(system.structure().0[left_states * n], 1.0);
        let duration: f64 = 0.7;
        let mut initial = vec![0.0; n];
        initial[0] = 1.0;
        if left_states == 2 {
            initial[1] = 0.3;
        }
        let mut endpoint = vec![0.0; n];
        let mut output = [0.0; 2];
        let mut prepared =
            PreparedExactTones::new(&system, &[0.0], duration, split(), dense()).unwrap();
        let record = prepared
            .step_into(&initial, &[0.0; 2], &[0.0; 2], &mut endpoint, &mut output)
            .unwrap();
        near(endpoint[0], duration.cos(), 2e-12, 2e-13);
        near(endpoint[left_states], duration.sin(), 2e-12, 2e-13);
        if left_states == 2 {
            near(endpoint[1], 0.3, 0.0, 1e-15);
        }
        assert_eq!(record.supplied, 0.0);
        assert_eq!(record.dissipated, 0.0);
    }
}

#[test]
fn invalid_maps_capacities_and_unrepresentable_products_refuse() {
    let make = || (inertias(&MASSES, 0.0), chamber(COMPLIANCE, 0.0));
    for (left, right, values) in [
        (vec![0, 0], vec![0], vec![0.1, -0.07]),
        (vec![4], vec![0], vec![0.1]),
        (vec![0], vec![2], vec![0.1]),
        (vec![0, 1], vec![0], vec![0.1]),
        (vec![0, 1], vec![0], vec![0.1, f64::NAN]),
        (vec![], vec![0], vec![]),
    ] {
        let (a, b) = make();
        assert!(interconnect_mapped(a, b, &left, &right, &values, budget()).is_err());
    }
    for controls in [
        InterconnectionBudget {
            max_states: 2,
            ..budget()
        },
        InterconnectionBudget {
            max_total_ports: 5,
            ..budget()
        },
        InterconnectionBudget {
            max_dense_bytes: 1,
            ..budget()
        },
        InterconnectionBudget {
            max_coupling_multiply_terms: 5,
            ..budget()
        },
    ] {
        let (a, b) = make();
        assert!(matches!(
            interconnect_mapped(a, b, &[0, 1], &[0], &AREAS, controls),
            Err(PhsError::Unadmitted { .. })
        ));
    }
    for (gain, coefficient) in [(1e-300, 1e-300), (1e300, 1e300), (f64::NAN, 1.0)] {
        // Deliberately exercise composition's original-matrix admission even
        // for a trusted/raw source that bypassed its own finite-data checks.
        let a = PortHamiltonian::from_raw_parts(
            1,
            1,
            vec![0.0],
            vec![0.0],
            vec![gain],
            Box::new(QuadraticStorage::new(vec![1.0], 1).unwrap()),
        );
        let b = matrix_system(1, 1, vec![1.0]);
        assert!(interconnect_mapped(a, b, &[0], &[0], &[coefficient], budget()).is_err());
    }
}
