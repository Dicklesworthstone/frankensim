//! G0/G1/G3/G4: source work uses the complete original physical trajectory.
use super::{close, distance, exact, free_mass, inertias, nonlinear, pair, piston_transfer};
use fs_phs::linear_exact::DenseExponentialBudget;
use fs_phs::multirate::{
    ContinuousSourcePlan, MultirateError, MultiratePair, SolveConfig, SourceInput,
};
use fs_phs::{duffing_oscillator, mass_spring_damper};

fn plan(frequencies: &[f64], affine: bool) -> ContinuousSourcePlan<'_> {
    ContinuousSourcePlan {
        angular_frequencies_rad_s: frequencies,
        affine_envelopes: affine,
        budget: DenseExponentialBudget {
            max_dimension: 64,
            max_dense_bytes: 128 * 1024 * 1024,
            max_multiply_terms: 20_000_000_000,
        },
    }
}

#[test]
fn zero_mean_continuous_force_recovers_real_regenerative_work_across_clocks() {
    // Independently integrate pL'=A*sin(w*t)+cL and pR'=cR. At a full
    // source period average(vL)=pL/mL+A/(mL*w)+h*cL/(2*mL).
    let (left, right) = (free_mass(2.0), free_mass(3.0));
    let (h, w, amplitude) = (1.0, std::f64::consts::TAU, 4.0);
    let c_right = (1.0 + amplitude / (2.0 * w) + 0.125) / (1.0 + 1.0 / 24.0);
    let c_left = 0.5 - c_right / 6.0;
    let source_work = -amplitude * c_left / (2.0 * w);
    for (nl, nr) in [(1, 1), (2, 4), (4, 16)] {
        let mut pair = MultiratePair::new(
            &left,
            &right,
            &pair(),
            h,
            exact(nl),
            exact(nr),
            SolveConfig::default(),
        )
        .unwrap()
        .with_continuous_sources(Some(plan(&[w], false)), None)
        .unwrap();
        let step = pair
            .step_sources(
                &[2.0],
                &[-1.5],
                SourceInput::Tones {
                    cosine: &[0.0],
                    sine: &[amplitude],
                },
                SourceInput::Held(&[0.0]),
            )
            .unwrap();
        close(step.left.state[0], 2.0 + c_left, 5e-13);
        close(step.right.state[0], -1.5 + c_right, 5e-13);
        close(step.left.average_output[0], c_right, 5e-13);
        close(step.right.average_output[0], -c_left, 5e-13);
        close(step.external_work_j[0], source_work, 5e-13);
        assert!(step.external_work_j[0] < -0.01);
        assert_eq!(step.external_work_j[1], 0.0);
        close(step.exchange_work_j[0], c_left * c_right, 5e-13);
        close(step.exchange_work_j.iter().sum(), 0.0, 5e-13);
        close(step.balance_residual_j, 0.0, 5e-13);
        assert_eq!(step.left.dissipated_j, 0.0);
        assert_eq!(step.right.dissipated_j, 0.0);
    }
}

#[test]
fn continuous_piston_forces_and_chamber_flow_keep_rectangular_transfer_work() {
    // Two independent sinusoidal piston forces and a sinusoidal volume source
    // share one compressible chamber. Integrate one full forcing period:
    // average(v_i)=p_i/m_i+F_i/(m_i*w)-h*A_i*average(P)/(2*m_i).
    // This independently determines both reciprocal held interface inputs and
    // each external source's signed work, without a sampled-force surrogate.
    let masses = [2.0, 3.0];
    let compliance = 0.004;
    let (left, right) = (inertias(&masses), free_mass(compliance));
    let (p, volume, amplitudes, flow_amplitude, h) = ([0.8, -0.3], 0.003, [0.3, -0.8], 0.01, 0.2);
    let w = std::f64::consts::TAU / h;
    let mut transfer = piston_transfer();
    transfer.left_input_scales = vec![0.3, 2.0];
    let areas = &transfer.right_from_left;
    let stiffness: f64 = (0..2).map(|i| areas[i] * areas[i] / masses[i]).sum();
    let free_flow: f64 = (0..2)
        .map(|i| areas[i] * (p[i] + amplitudes[i] / w) / masses[i])
        .sum();
    let pressure = (volume / compliance
        + flow_amplitude / (compliance * w)
        + h * free_flow / (2.0 * compliance))
        / (1.0 + h * h * stiffness / (4.0 * compliance));
    let velocity: [f64; 2] = std::array::from_fn(|i| {
        (p[i] + amplitudes[i] / w - h * areas[i] * pressure / 2.0) / masses[i]
    });
    let flow: f64 = (0..2).map(|i| areas[i] * velocity[i]).sum();
    let piston_work: f64 = (0..2)
        .map(|i| amplitudes[i] * areas[i] * pressure * h / (masses[i] * w))
        .sum();
    let chamber_work = -flow_amplitude * flow * h / (compliance * w);
    for (nl, nr) in [(1, 1), (2, 8), (16, 4)] {
        let mut coupling = MultiratePair::new_mapped(
            &left,
            &right,
            &transfer,
            h,
            exact(nl),
            exact(nr),
            SolveConfig::default(),
        )
        .unwrap()
        .with_continuous_sources(Some(plan(&[w], false)), Some(plan(&[w], false)))
        .unwrap();
        let r = coupling
            .step_sources(
                &p,
                &[volume],
                SourceInput::Tones {
                    cosine: &[0.0; 2],
                    sine: &amplitudes,
                },
                SourceInput::Tones {
                    cosine: &[0.0],
                    sine: &[flow_amplitude],
                },
            )
            .unwrap();
        for i in 0..2 {
            close(r.left.state[i], p[i] - h * areas[i] * pressure, 3e-13);
            close(r.left.average_output[i], velocity[i], 3e-13);
        }
        close(r.right.state[0], volume + h * flow, 3e-13);
        close(r.right.average_output[0], pressure, 3e-13);
        close(r.external_work_j[0], piston_work, 3e-13);
        close(r.external_work_j[1], chamber_work, 3e-13);
        close(r.exchange_work_j[0], -h * pressure * flow, 3e-13);
        close(r.exchange_work_j[1], h * pressure * flow, 3e-13);
        close(r.balance_residual_j, 0.0, 3e-13);
        assert!(r.external_work_j[1] < 0.0);
    }
}

fn rk4<const N: usize>(
    mut x: [f64; N],
    duration: f64,
    steps: usize,
    rhs: impl Fn(f64, &[f64; N]) -> [f64; N],
) -> [f64; N] {
    let h = duration / steps as f64;
    for step in 0..steps {
        let t = step as f64 * h;
        let a = rhs(t, &x);
        let b = rhs(t + h / 2.0, &std::array::from_fn(|i| x[i] + h * a[i] / 2.0));
        let c = rhs(t + h / 2.0, &std::array::from_fn(|i| x[i] + h * b[i] / 2.0));
        let d = rhs(t + h, &std::array::from_fn(|i| x[i] + h * c[i]));
        for i in 0..N {
            x[i] += h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
    x
}

fn force(t: f64, h: f64, w: &[f64], coefficients: &[&[f64]; 4]) -> f64 {
    (0..w.len())
        .map(|k| {
            let c = (1.0 - t / h) * coefficients[0][k] + (t / h) * coefficients[2][k];
            let s = (1.0 - t / h) * coefficients[1][k] + (t / h) * coefficients[3][k];
            c * (w[k] * t).cos() + s * (w[k] * t).sin()
        })
        .sum()
}

fn affine<'a>(c: &[&'a [f64]; 4]) -> SourceInput<'a> {
    SourceInput::AffineTones {
        cosine_start: c[0],
        sine_start: c[1],
        cosine_end: c[2],
        sine_end: c[3],
    }
}

#[test]
fn affine_beating_sources_match_independent_motion_heat_and_work_on_both_native_clocks() {
    let (left, right) = (
        mass_spring_damper(1.2, 3.7, 0.17).unwrap(),
        mass_spring_damper(0.8, 2.3, 0.31).unwrap(),
    );
    let h = 0.4;
    let (wl, wr) = ([0.0, 1.7, 3.6], [0.0, 2.3]);
    let l: [&[f64]; 4] = [
        &[0.2, 0.6, -0.3],
        &[0.0, -0.4, 0.2],
        &[-0.1, 0.1, 0.5],
        &[0.0, 0.3, -0.2],
    ];
    let r: [&[f64]; 4] = [&[-0.1, 0.4], &[0.0, 0.2], &[0.3, -0.2], &[0.0, -0.5]];
    let evolve = |side: usize, coupled: f64| {
        let (m, k, damping, q, p) = if side == 0 {
            (1.2, 3.7, 0.17, 0.3, -0.1)
        } else {
            (0.8, 2.3, 0.31, -0.2, 0.4)
        };
        rk4([q, p, 0.0, 0.0, 0.0, 0.0], h, 2048, |t, x| {
            let f = if side == 0 {
                force(t, h, &wl, &l)
            } else {
                force(t, h, &wr, &r)
            };
            let v = x[1] / m;
            [
                v,
                f + coupled - k * x[0] - damping * v,
                v,
                f * v,
                coupled * v,
                damping * v * v,
            ]
        })
    };
    // A separate physical ODE determines the affine interface equation.
    let (l0, r0) = (evolve(0, 0.0), evolve(1, 0.0));
    let (a, b) = (
        (evolve(0, 1.0)[2] - l0[2]) / h,
        (evolve(1, 1.0)[2] - r0[2]) / h,
    );
    let c_right = (l0[2] / h - a * r0[2] / h) / (1.0 + a * b);
    let c_left = -r0[2] / h - b * c_right;
    let truth = [evolve(0, c_left), evolve(1, c_right)];
    for (nl, nr) in [(1, 1), (2, 4), (4, 8)] {
        let mut coupling = MultiratePair::new(
            &left,
            &right,
            &pair(),
            h,
            exact(nl),
            exact(nr),
            SolveConfig::default(),
        )
        .unwrap()
        .with_continuous_sources(Some(plan(&wl, true)), Some(plan(&wr, true)))
        .unwrap();
        let actual = coupling
            .step_sources(&[0.3, -0.1], &[-0.2, 0.4], affine(&l), affine(&r))
            .unwrap();
        for (side, trajectory) in [&actual.left, &actual.right].into_iter().enumerate() {
            for i in 0..2 {
                close(trajectory.state[i], truth[side][i], 2e-12);
            }
            close(trajectory.average_output[0], truth[side][2] / h, 2e-12);
            close(actual.external_work_j[side], truth[side][3], 2e-12);
            close(actual.exchange_work_j[side], truth[side][4], 2e-12);
            close(trajectory.dissipated_j, truth[side][5], 2e-12);
        }
        close(actual.exchange_work_j.iter().sum(), 0.0, 2e-12);
        close(actual.balance_residual_j, 0.0, 2e-12);
    }
}

#[test]
fn tiny_external_work_is_not_the_subtraction_of_large_interface_work() {
    let (left, right) = (free_mass(2.0), free_mass(3.0));
    let mut coupling = MultiratePair::new(
        &left,
        &right,
        &pair(),
        1.0,
        exact(2),
        exact(4),
        SolveConfig::default(),
    )
    .unwrap()
    .with_continuous_sources(Some(plan(&[0.0], false)), None)
    .unwrap();
    let force = 1e-24;
    let r = coupling
        .step_sources(
            &[2.0],
            &[-1.5],
            SourceInput::Tones {
                cosine: &[force],
                sine: &[0.0],
            },
            SourceInput::Held(&[0.0]),
        )
        .unwrap();
    assert!(r.exchange_work_j[0].abs() > 0.1);
    assert!(r.external_work_j[0] > 0.0);
    close(
        r.external_work_j[0] / force,
        r.left.average_output[0],
        3e-13,
    );
    // The deliberately tempting remainder has already lost this physical term.
    assert!(
        (r.left.supplied_j - r.exchange_work_j[0]).abs() > 1e6 * r.external_work_j[0]
            || r.left.supplied_j == r.exchange_work_j[0]
    );
}

#[test]
fn continuous_source_can_drive_an_original_nonlinear_partner_without_changing_its_clock() {
    let (left, right) = (
        mass_spring_damper(1.0, 2.0, 0.07).unwrap(),
        duffing_oscillator(1.2, 1.3, 0.6, 0.02).unwrap(),
    );
    let mut coupling = MultiratePair::new(
        &left,
        &right,
        &pair(),
        0.02,
        exact(2),
        nonlinear(4),
        SolveConfig::default(),
    )
    .unwrap()
    .with_continuous_sources(Some(plan(&[2.7], false)), None)
    .unwrap();
    let (mut l, mut r) = (vec![0.2, -0.1], vec![-0.4, 0.3]);
    let initial = left.hamiltonian(&l) + right.hamiltonian(&r);
    let (mut c, mut s, mut work, mut heat) = (0.7, -0.2, 0.0, 0.0);
    for _ in 0..60 {
        let step = coupling
            .step_sources(
                &l,
                &r,
                SourceInput::Tones {
                    cosine: &[c],
                    sine: &[s],
                },
                SourceInput::Substeps(&[0.1, -0.05, 0.0, 0.2]),
            )
            .unwrap();
        assert_eq!((step.left.substeps, step.right.substeps), (2, 4));
        close(step.exchange_work_j.iter().sum(), 0.0, 2e-12);
        work += step.external_work_j.iter().sum::<f64>();
        heat += step.left.dissipated_j + step.right.dissipated_j;
        l = step.left.state;
        r = step.right.state;
        (c, s) = coupling.advance_left_force_phase(0, c, s).unwrap();
    }
    close(
        left.hamiltonian(&l) + right.hamiltonian(&r) - initial + heat - work,
        0.0,
        2e-11,
    );
    assert!(heat > 0.001);
}

#[test]
fn refining_the_continuously_forced_interface_converges_to_the_coupled_physical_ode() {
    let (left, right) = (
        mass_spring_damper(1.0, 2.0, 0.1).unwrap(),
        mass_spring_damper(1.3, 2.7, 0.2).unwrap(),
    );
    let duration = 1.2;
    let initial = [0.2, -0.1, -0.3, 0.4];
    let truth = rk4(initial, duration, 8192, |t, x| {
        let (vl, vr) = (x[1], x[3] / 1.3);
        [
            vl,
            0.7 * (2.9 * t).cos() - 0.2 * (2.9 * t).sin() - 2.0 * x[0] - 0.1 * vl - vr,
            vr,
            -0.4 * (1.7 * t).cos() + 0.3 * (1.7 * t).sin() - 2.7 * x[2] - 0.2 * vr + vl,
        ]
    });
    let mut errors = Vec::new();
    for windows in [12, 24, 48] {
        let mut coupling = MultiratePair::new(
            &left,
            &right,
            &pair(),
            duration / windows as f64,
            exact(2),
            exact(4),
            SolveConfig::default(),
        )
        .unwrap()
        .with_continuous_sources(Some(plan(&[2.9], false)), Some(plan(&[1.7], false)))
        .unwrap();
        let (mut a, mut b) = (initial[..2].to_vec(), initial[2..].to_vec());
        let (mut cl, mut sl, mut cr, mut sr) = (0.7, -0.2, -0.4, 0.3);
        for _ in 0..windows {
            let step = coupling
                .step_sources(
                    &a,
                    &b,
                    SourceInput::Tones {
                        cosine: &[cl],
                        sine: &[sl],
                    },
                    SourceInput::Tones {
                        cosine: &[cr],
                        sine: &[sr],
                    },
                )
                .unwrap();
            a = step.left.state;
            b = step.right.state;
            (cl, sl) = coupling.advance_left_force_phase(0, cl, sl).unwrap();
            (cr, sr) = coupling.advance_right_force_phase(0, cr, sr).unwrap();
        }
        a.extend(b);
        errors.push(distance(&a, &truth));
    }
    assert!(
        errors
            .windows(2)
            .all(|v| v[0] / v[1] > 3.8 && v[0] / v[1] < 4.2),
        "{errors:?}"
    );
}

#[test]
fn continuous_source_admission_and_cancellation_leave_the_previous_interface_reusable() {
    let (left, right) = (free_mass(2.0), free_mass(3.0));
    let base = || {
        MultiratePair::new(
            &left,
            &right,
            &pair(),
            0.1,
            exact(2),
            exact(4),
            SolveConfig::default(),
        )
        .unwrap()
    };
    let make = || {
        base()
            .with_continuous_sources(Some(plan(&[3.0], false)), None)
            .unwrap()
    };
    let (a, b) = ([2.0], [-1.5]);
    let input = SourceInput::Tones {
        cosine: &[0.2],
        sine: &[-0.3],
    };
    assert!(
        base()
            .step_sources(&a, &b, input, SourceInput::Held(&[0.0]))
            .is_err()
    );
    let mut coupling = make();
    let mut reference = make();
    let accepted = coupling
        .step_sources(&a, &b, input, SourceInput::Held(&[0.0]))
        .unwrap();
    assert_eq!(
        accepted,
        reference
            .step_sources(&a, &b, input, SourceInput::Held(&[0.0]))
            .unwrap()
    );
    for invalid in [
        SourceInput::Tones {
            cosine: &[],
            sine: &[0.0],
        },
        SourceInput::Tones {
            cosine: &[f64::NAN],
            sine: &[0.0],
        },
        SourceInput::AffineTones {
            cosine_start: &[0.0],
            sine_start: &[0.0],
            cosine_end: &[1.0],
            sine_end: &[0.0],
        },
    ] {
        assert!(matches!(
            coupling.step_sources(&a, &b, invalid, SourceInput::Held(&[0.0])),
            Err(MultirateError::Invalid(_))
        ));
    }
    let mut polls = 0;
    assert!(matches!(
        coupling.step_sources_controlled(&a, &b, input, SourceInput::Held(&[0.0]), || {
            polls += 1;
            polls >= 7
        }),
        Err(MultirateError::Cancelled)
    ));
    assert_eq!(
        coupling
            .step_sources(&a, &b, input, SourceInput::Held(&[0.0]))
            .unwrap(),
        reference
            .step_sources(&a, &b, input, SourceInput::Held(&[0.0]))
            .unwrap()
    );
    assert!(
        base()
            .with_continuous_sources(Some(plan(&[f64::NAN], false)), None)
            .is_err()
    );
    let insufficient = ContinuousSourcePlan {
        budget: DenseExponentialBudget {
            max_dimension: 2,
            ..plan(&[3.0], false).budget
        },
        ..plan(&[3.0], false)
    };
    assert!(
        base()
            .with_continuous_sources(Some(insufficient), None)
            .is_err()
    );
    assert!(matches!(
        base().with_continuous_sources_controlled(Some(plan(&[3.0], false)), None, || true),
        Err(MultirateError::Cancelled)
    ));
    let nonlinear_pair = MultiratePair::new(
        &left,
        &right,
        &pair(),
        0.1,
        nonlinear(2),
        exact(4),
        SolveConfig::default(),
    )
    .unwrap();
    assert!(matches!(
        nonlinear_pair.with_continuous_sources(Some(plan(&[3.0], false)), None),
        Err(MultirateError::Invalid(_))
    ));
    assert!(coupling.advance_left_force_phase(1, 0.0, 0.0).is_err());
    assert!(coupling.advance_right_force_phase(0, 0.0, 0.0).is_err());
}
