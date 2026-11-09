//! Analytic physical frequency/parameter jets, without frequency perturbations.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_math::c64::C64;
use fs_modal::harmonic::{
    HarmonicFrequencyParameterControls, HarmonicPortControls, HarmonicResponseControls,
    HarmonicResponseError, PreparedHarmonicResponse,
};
use fs_modal::quadratic_modes::QuadraticParameterDirection;

fn with_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            gate,
            arena,
            StreamKey {
                seed: 79,
                kernel_id: 251,
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
        maximum_coordinates: 8,
        maximum_payload_bytes: 1 << 20,
        maximum_relative_residual: 2e-11,
    }
}
fn limits() -> HarmonicFrequencyParameterControls {
    HarmonicFrequencyParameterControls {
        maximum_additional_payload_bytes: 1 << 20,
        maximum_work: 1_000_000,
    }
}
fn close(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-280),
        "{actual:?} != {expected:?}",
    );
}

/// Closed determinant/quotient-rule oracle; this never calls a response solver.
#[derive(Clone, Copy)]
struct Jet {
    v: C64,
    w: C64,
    ww: C64,
    t: C64,
    wt: C64,
}
impl Jet {
    fn sub(self, other: Self) -> Self {
        Self {
            v: self.v - other.v,
            w: self.w - other.w,
            ww: self.ww - other.ww,
            t: self.t - other.t,
            wt: self.wt - other.wt,
        }
    }
    fn mul(self, other: Self) -> Self {
        Self {
            v: self.v * other.v,
            w: self.w * other.v + self.v * other.w,
            ww: self.ww * other.v + (self.w * other.w).scale(2.0) + self.v * other.ww,
            t: self.t * other.v + self.v * other.t,
            wt: self.wt * other.v + self.w * other.t + self.t * other.w + self.v * other.wt,
        }
    }
    fn div(self, denominator: Self) -> Self {
        let d2 = denominator.v * denominator.v;
        let d3 = d2 * denominator.v;
        Self {
            v: self.v / denominator.v,
            w: self.w / denominator.v - self.v * denominator.w / d2,
            ww: self.ww / denominator.v
                - (self.w * denominator.w).scale(2.0) / d2
                - self.v * denominator.ww / d2
                + (self.v * denominator.w * denominator.w).scale(2.0) / d3,
            t: self.t / denominator.v - self.v * denominator.t / d2,
            wt: self.wt / denominator.v
                - (self.w * denominator.t + self.t * denominator.w + self.v * denominator.wt) / d2
                + (self.v * denominator.w * denominator.t).scale(2.0) / d3,
        }
    }
}
fn coefficient(m: f64, c: f64, k: f64, dm: f64, dc: f64, dk: f64, w: f64) -> Jet {
    Jet {
        v: C64::new(k - w * w * m, -w * c),
        w: C64::new(-2.0 * w * m, -c),
        ww: C64::from_re(-2.0 * m),
        t: C64::new(dk - w * w * dm, -w * dc),
        wt: C64::new(-2.0 * w * dm, -dc),
    }
}
fn source() -> Jet {
    Jet {
        v: C64::new(1.7, -0.6),
        w: C64::new(-0.2, 0.8),
        ww: C64::new(0.3, -0.4),
        t: C64::new(0.1, 0.7),
        wt: C64::new(-0.6, 0.2),
    }
}

#[test]
fn all_scalar_source_and_coefficient_jets_retain_high_q_phase_and_dc_velocity() {
    with_cx(&CancelGate::new(), |cx| {
        let (m, c, k) = (2.0, 0.00002, 8.0);
        let (dm, dc, dk) = (0.3, -0.004, 1.7);
        let f = source();
        for w in [0.0, 0.3, 1.99999, 2.0, 2.00001, 10.0] {
            let prepared =
                PreparedHarmonicResponse::try_new(&[m], &[c], &[k], 1, w, controls(), cx).unwrap();
            let result = prepared
                .solve_with_frequency_parameter_derivative(
                    &[f.v],
                    QuadraticParameterDirection {
                        mass: &[dm],
                        damping: &[dc],
                        stiffness: &[dk],
                    },
                    Some(&[f.w]),
                    Some(&[f.ww]),
                    Some(&[f.t]),
                    Some(&[f.wt]),
                    limits(),
                    cx,
                )
                .unwrap();
            let expected = f.div(coefficient(m, c, k, dm, dc, dk, w));
            close(
                result
                    .frequency_response
                    .frequency_response
                    .response
                    .displacement[0],
                expected.v,
                3e-11,
            );
            close(
                result
                    .frequency_response
                    .frequency_response
                    .displacement_derivative[0],
                expected.w,
                3e-11,
            );
            close(
                result.frequency_response.displacement_second_derivative[0],
                expected.ww,
                3e-11,
            );
            close(
                result.parameter_response.displacement_derivative[0],
                expected.t,
                3e-11,
            );
            close(result.displacement_mixed_derivative[0], expected.wt, 3e-11);
            close(
                result.velocity_mixed_derivative[0],
                C64::new(0.0, -1.0) * (expected.t + expected.wt.scale(w)),
                3e-11,
            );
            assert!(result.relative_mixed_derivative_residual < 2e-13);
        }
    });
}

#[test]
fn passive_two_mass_antiresonance_keeps_the_zero_and_physical_notch_shift() {
    with_cx(&CancelGate::new(), |cx| {
        // A grounded primary mass and lossless attached secondary mass.
        // The receiver q[0] vanishes exactly at w=1 while q[1]=-1.
        // Changing the attachment spring moves the notch to sqrt(1+theta).
        let m = [1.0, 0.0, 0.0, 1.0];
        let c = [0.04, 0.0, 0.0, 0.0];
        let k = [2.0, -1.0, -1.0, 1.0];
        let zero = [0.0; 4];
        let dk = [1.0, -1.0, -1.0, 1.0];
        let prepared =
            PreparedHarmonicResponse::try_new(&m, &c, &k, 2, 1.0, controls(), cx).unwrap();
        let response = prepared
            .solve_with_frequency_parameter_derivative(
                &[C64::ONE, C64::ZERO],
                QuadraticParameterDirection {
                    mass: &zero,
                    damping: &zero,
                    stiffness: &dk,
                },
                None,
                None,
                None,
                None,
                limits(),
                cx,
            )
            .unwrap();
        let q = response
            .frequency_response
            .frequency_response
            .response
            .displacement[0];
        let qw = response
            .frequency_response
            .frequency_response
            .displacement_derivative[0];
        let qt = response.parameter_response.displacement_derivative[0];
        assert!(q.abs() < 2e-15);
        close(qw, C64::from_re(2.0), 2e-13);
        close(qt, C64::from_re(-1.0), 2e-13);
        close(
            response.displacement_mixed_derivative[0],
            C64::new(0.0, -0.16),
            2e-12,
        );
        let notch_shift = -(qw.conj() * qt).re / qw.norm_sq();
        assert!((notch_shift - 0.5).abs() < 2e-13);
        assert!(response.relative_mixed_derivative_residual < 2e-13);
    });
}

#[test]
fn dense_nonproportional_gyroscopic_directions_keep_binary_coordinate_units() {
    with_cx(&CancelGate::new(), |cx| {
        let m = [2.0, 0.2, 0.2, 1.1];
        let c = [0.4, 0.7, -0.3, 0.6];
        let k = [7.0, -1.0, -1.0, 3.0];
        let dm = [0.1, 0.02, 0.02, -0.04];
        let dc = [0.03, 0.15, -0.11, -0.02];
        let dk = [0.6, -0.2, -0.2, 0.3];
        let w = 1.4;
        let f = [
            source(),
            Jet {
                v: C64::new(-0.2, 0.4),
                w: C64::new(0.6, -0.2),
                ww: C64::new(-0.1, 0.3),
                t: C64::new(0.4, 0.2),
                wt: C64::new(0.3, -0.1),
            },
        ];
        let a: [Jet; 4] =
            core::array::from_fn(|i| coefficient(m[i], c[i], k[i], dm[i], dc[i], dk[i], w));
        let det = a[0].mul(a[3]).sub(a[1].mul(a[2]));
        let expected = [
            a[3].mul(f[0]).sub(a[1].mul(f[1])).div(det),
            a[0].mul(f[1]).sub(a[2].mul(f[0])).div(det),
        ];
        for exponents in [[0, 0], [450, -450]] {
            let units = exponents.map(|exponent| 2.0_f64.powi(exponent));
            let transform = |values: [f64; 4]| {
                core::array::from_fn::<_, 4, _>(|i| values[i] * units[i / 2] * units[i % 2])
            };
            let (mass, damping, stiffness) = (transform(m), transform(c), transform(k));
            let (mass_direction, damping_direction, stiffness_direction) =
                (transform(dm), transform(dc), transform(dk));
            let force = core::array::from_fn::<_, 2, _>(|i| f[i].v.scale(units[i]));
            let fw = core::array::from_fn::<_, 2, _>(|i| f[i].w.scale(units[i]));
            let fww = core::array::from_fn::<_, 2, _>(|i| f[i].ww.scale(units[i]));
            let ft = core::array::from_fn::<_, 2, _>(|i| f[i].t.scale(units[i]));
            let fwt = core::array::from_fn::<_, 2, _>(|i| f[i].wt.scale(units[i]));
            let prepared = PreparedHarmonicResponse::try_new(
                &mass,
                &damping,
                &stiffness,
                2,
                w,
                controls(),
                cx,
            )
            .unwrap();
            let result = prepared
                .solve_with_frequency_parameter_derivative(
                    &force,
                    QuadraticParameterDirection {
                        mass: &mass_direction,
                        damping: &damping_direction,
                        stiffness: &stiffness_direction,
                    },
                    Some(&fw),
                    Some(&fww),
                    Some(&ft),
                    Some(&fwt),
                    limits(),
                    cx,
                )
                .unwrap();
            for i in 0..2 {
                close(
                    result
                        .frequency_response
                        .frequency_response
                        .response
                        .displacement[i]
                        .scale(units[i]),
                    expected[i].v,
                    3e-11,
                );
                close(
                    result
                        .frequency_response
                        .frequency_response
                        .displacement_derivative[i]
                        .scale(units[i]),
                    expected[i].w,
                    3e-11,
                );
                close(
                    result.frequency_response.displacement_second_derivative[i].scale(units[i]),
                    expected[i].ww,
                    3e-11,
                );
                close(
                    result.parameter_response.displacement_derivative[i].scale(units[i]),
                    expected[i].t,
                    3e-11,
                );
                close(
                    result.displacement_mixed_derivative[i].scale(units[i]),
                    expected[i].wt,
                    3e-11,
                );
                close(
                    result.velocity_mixed_derivative[i].scale(units[i]),
                    C64::new(0.0, -1.0) * (expected[i].t + expected[i].wt.scale(w)),
                    3e-11,
                );
            }
            assert!(result.relative_mixed_derivative_residual < 2e-13);
        }
    });
}

#[test]
fn mixed_jet_refusals_are_explicit_and_leave_the_factor_reusable() {
    with_cx(&CancelGate::new(), |cx| {
        let prepared =
            PreparedHarmonicResponse::try_new(&[1.0], &[0.1], &[1.0], 1, 0.8, controls(), cx)
                .unwrap();
        let zero = [0.0];
        let direction = QuadraticParameterDirection {
            mass: &zero,
            damping: &zero,
            stiffness: &zero,
        };
        let solve = |limits, mixed_force: Option<&[C64]>, context: &Cx<'_>| {
            prepared.solve_with_frequency_parameter_derivative(
                &[C64::ONE],
                direction,
                None,
                None,
                None,
                mixed_force,
                limits,
                context,
            )
        };
        assert!(matches!(
            solve(
                HarmonicFrequencyParameterControls {
                    maximum_additional_payload_bytes: 3071,
                    ..limits()
                },
                None,
                cx
            ),
            Err(HarmonicResponseError::Budget {
                requested: 3072,
                ..
            }),
        ));
        assert!(matches!(
            solve(
                HarmonicFrequencyParameterControls {
                    maximum_work: 3071,
                    ..limits()
                },
                None,
                cx
            ),
            Err(HarmonicResponseError::WorkBudget {
                requested: 3072,
                ..
            }),
        ));
        assert!(matches!(
            solve(limits(), Some(&[]), cx),
            Err(HarmonicResponseError::Invalid(_))
        ));
        assert!(matches!(
            solve(limits(), Some(&[C64::new(f64::NAN, 0.0)]), cx),
            Err(HarmonicResponseError::Invalid(_)),
        ));
        let gate = CancelGate::new();
        gate.request();
        with_cx(&gate, |cancelled| {
            assert!(matches!(
                solve(limits(), None, cancelled),
                Err(HarmonicResponseError::Cancelled)
            ));
        });
        let baseline = prepared.solve(&[C64::ONE], cx).unwrap();
        let after = solve(limits(), None, cx).unwrap();
        assert_eq!(
            after
                .frequency_response
                .frequency_response
                .response
                .displacement,
            baseline.displacement
        );
        assert_eq!(after.displacement_mixed_derivative, vec![C64::ZERO]);

        let loaded = PreparedHarmonicResponse::try_with_port_impedance(
            &[1.0],
            &[0.1],
            &[1.0],
            1,
            &[1.0],
            1,
            &[C64::from_re(0.2)],
            0.8,
            controls(),
            HarmonicPortControls {
                maximum_ports: 1,
                maximum_projection_terms: 100,
                maximum_payload_bytes: 1 << 20,
            },
            cx,
        )
        .unwrap();
        assert!(matches!(
            loaded.solve_with_frequency_parameter_derivative(
                &[C64::ONE],
                direction,
                None,
                None,
                None,
                None,
                limits(),
                cx,
            ),
            Err(HarmonicResponseError::Invalid(
                "mixed frequency/parameter derivatives require an unloaded mechanical pencil"
            )),
        ));
    });
}

fn port_limits() -> HarmonicPortControls {
    HarmonicPortControls {
        maximum_ports: 4,
        maximum_projection_terms: 100_000,
        maximum_payload_bytes: 1 << 20,
    }
}

fn check_physical_jet(
    actual: &fs_modal::harmonic::HarmonicFrequencyParameterResponse,
    expected: &[Jet],
    w: f64,
    tolerance: f64,
) {
    for (i, expected) in expected.iter().enumerate() {
        close(
            actual
                .frequency_response
                .frequency_response
                .response
                .displacement[i],
            expected.v,
            tolerance,
        );
        close(
            actual
                .frequency_response
                .frequency_response
                .displacement_derivative[i],
            expected.w,
            tolerance,
        );
        close(
            actual.frequency_response.displacement_second_derivative[i],
            expected.ww,
            tolerance,
        );
        close(
            actual.parameter_response.displacement_derivative[i],
            expected.t,
            tolerance,
        );
        close(
            actual.displacement_mixed_derivative[i],
            expected.wt,
            tolerance,
        );
        close(
            actual.velocity_mixed_derivative[i],
            C64::new(0.0, -1.0) * (expected.t + expected.wt.scale(w)),
            tolerance,
        );
    }
    assert!(actual.relative_mixed_derivative_residual < 2e-12);
}

#[test]
fn dispersive_fixed_load_mixed_jets_match_original_added_mass_stiffness_and_damping() {
    with_cx(&CancelGate::new(), |cx| {
        let (m, c, k) = (1.7, 0.07, 4.3);
        let (dm, dc, dk) = (0.12, -0.03, 0.4);
        let (lm, lc, lk, b) = (0.7, 0.4, 0.2, 1.8);
        let f = source();
        for w in [0.04, 0.5, 1.4, 2.0, 8.0] {
            let z = C64::new(lc, lk / w - w * lm);
            let zw = C64::new(0.0, -lk / (w * w) - lm);
            let zww = C64::new(0.0, 2.0 * lk / (w * w * w));
            let owner = PreparedHarmonicResponse::try_with_port_impedance(
                &[m],
                &[c],
                &[k],
                1,
                &[b],
                1,
                &[z],
                w,
                controls(),
                port_limits(),
                cx,
            )
            .unwrap();
            let result = owner
                .solve_with_loaded_frequency_parameter_derivative(
                    &[f.v],
                    QuadraticParameterDirection {
                        mass: &[dm],
                        damping: &[dc],
                        stiffness: &[dk],
                    },
                    Some(&[f.w]),
                    Some(&[f.ww]),
                    Some(&[f.t]),
                    Some(&[f.wt]),
                    &[zw],
                    &[zww],
                    limits(),
                    cx,
                )
                .unwrap();
            let expected = f.div(coefficient(
                m + b * b * lm,
                c + b * b * lc,
                k + b * b * lk,
                dm,
                dc,
                dk,
                w,
            ));
            check_physical_jet(&result, &[expected], w, 3e-11);
            assert!(result.frequency_response.port_load.is_some());
            assert!(result.parameter_response.port_load.is_some());
        }
    });
}

#[test]
fn maxwell_memory_mixed_jets_match_independent_original_cubic_quotient_at_dc_and_resonance() {
    with_cx(&CancelGate::new(), |cx| {
        let (m, c, k, kb, tau, b) = (1.7, 0.06, 4.3, 2.1, 0.8, 1.3);
        let (dm, dc, dk) = (0.2, -0.01, 0.3);
        let f = source();
        let zero = Jet {
            v: C64::ZERO,
            w: C64::ZERO,
            ww: C64::ZERO,
            t: C64::ZERO,
            wt: C64::ZERO,
        };
        for w in [0.0, 0.25, 1.5, 2.0, 10.0] {
            let n = Jet {
                v: C64::new(1.0, -w * tau),
                w: C64::new(0.0, -tau),
                ..zero
            };
            let z = C64::from_re(kb * tau) / n.v;
            let zw = C64::new(0.0, kb * tau * tau) / (n.v * n.v);
            let zww = C64::from_re(-2.0 * kb * tau * tau * tau) / (n.v * n.v * n.v);
            let owner = PreparedHarmonicResponse::try_with_port_impedance(
                &[m],
                &[c],
                &[k],
                1,
                &[b],
                1,
                &[z],
                w,
                controls(),
                port_limits(),
                cx,
            )
            .unwrap();
            let result = owner
                .solve_with_loaded_frequency_parameter_derivative(
                    &[f.v],
                    QuadraticParameterDirection {
                        mass: &[dm],
                        damping: &[dc],
                        stiffness: &[dk],
                    },
                    Some(&[f.w]),
                    Some(&[f.ww]),
                    Some(&[f.t]),
                    Some(&[f.wt]),
                    &[zw],
                    &[zww],
                    limits(),
                    cx,
                )
                .unwrap();
            // Multiplying out the retained memory state gives the actual
            // original cubic denominator, including the fixed port-map square.
            let mechanical = coefficient(m, c, k, dm, dc, dk, w);
            let memory = Jet {
                v: C64::new(0.0, w * b * b * kb * tau),
                w: C64::new(0.0, b * b * kb * tau),
                ..zero
            };
            let expected = f.mul(n).div(mechanical.mul(n).sub(memory));
            check_physical_jet(&result, &[expected], w, 3e-11);
        }
    });
}

#[test]
fn finite_fixed_damper_retains_loaded_mixed_kinematics_at_dc() {
    with_cx(&CancelGate::new(), |cx| {
        let (m, c, k, b, z) = (2.0, 0.1, 8.0, 3.0, 0.2);
        let (dm, dc, dk) = (0.3, 0.04, -0.2);
        let f = source();
        for w in [0.0, 0.1, 2.0] {
            let owner = PreparedHarmonicResponse::try_with_port_impedance(
                &[m],
                &[c],
                &[k],
                1,
                &[b],
                1,
                &[C64::from_re(z)],
                w,
                controls(),
                port_limits(),
                cx,
            )
            .unwrap();
            let result = owner
                .solve_with_loaded_frequency_parameter_derivative(
                    &[f.v],
                    QuadraticParameterDirection {
                        mass: &[dm],
                        damping: &[dc],
                        stiffness: &[dk],
                    },
                    Some(&[f.w]),
                    Some(&[f.ww]),
                    Some(&[f.t]),
                    Some(&[f.wt]),
                    &[C64::ZERO],
                    &[C64::ZERO],
                    limits(),
                    cx,
                )
                .unwrap();
            let expected = f.div(coefficient(m, c + b * b * z, k, dm, dc, dk, w));
            check_physical_jet(&result, &[expected], w, 3e-11);
            if w == 0.0 {
                close(
                    result.velocity_mixed_derivative[0],
                    C64::new(0.0, -1.0) * result.parameter_response.displacement_derivative[0],
                    2e-13,
                );
            }
        }
    });
}

#[test]
fn compact_mutual_nonreciprocal_loaded_mixed_jets_match_full_effective_mechanics() {
    with_cx(&CancelGate::new(), |cx| {
        let natural = [1.0, 2.0, 3.0];
        let damping = [0.0, 0.04, 0.08];
        let b = [1.0, 0.2, -0.4, 0.6, 0.7, -0.3];
        let lm = [0.2, 0.03, 0.03, 0.15];
        let lc = [0.4, 0.12, -0.04, 0.3];
        let lk = [0.1, -0.02, -0.02, 0.2];
        // The load's positive symmetric loss and full skew velocity coupling
        // both survive B*Z*B^T; parameter directions may also couple all modes.
        let dm = [0.1, 0.02, -0.01, 0.02, -0.04, 0.03, -0.01, 0.03, 0.06];
        let dc = [0.03, 0.15, 0.01, -0.11, -0.02, 0.04, -0.01, 0.02, 0.01];
        let dk = [0.6, -0.2, 0.1, -0.2, 0.3, -0.04, 0.1, -0.04, 0.2];
        let project = |a: &[f64; 4], symmetric: bool| {
            let mut values = [0.0; 9];
            for i in 0..3 {
                for j in 0..3 {
                    if symmetric && j < i {
                        values[3 * i + j] = values[3 * j + i];
                        continue;
                    }
                    for p in 0..2 {
                        for q in 0..2 {
                            values[3 * i + j] += b[2 * i + p] * a[2 * p + q] * b[2 * j + q];
                        }
                    }
                }
            }
            values
        };
        let mut mass = project(&lm, true);
        let mut velocity = project(&lc, false);
        let mut stiffness = project(&lk, true);
        for i in 0..3 {
            mass[3 * i + i] += 1.0;
            velocity[3 * i + i] += damping[i];
            stiffness[3 * i + i] += natural[i] * natural[i];
        }
        let force = [C64::new(1.0, 0.2), C64::new(-0.3, 0.8), C64::new(0.7, -0.5)];
        let fw = [C64::new(-0.2, 0.1), C64::new(0.6, -0.2), C64::new(0.1, 0.3)];
        let fww = [C64::new(0.1, 0.4), C64::new(-0.1, 0.3), C64::new(-0.2, 0.1)];
        let ft = [C64::new(0.4, 0.2), C64::new(0.2, -0.1), C64::new(-0.1, 0.3)];
        let fwt = [
            C64::new(0.3, -0.1),
            C64::new(-0.1, 0.2),
            C64::new(0.1, -0.4),
        ];
        let direction = QuadraticParameterDirection {
            mass: &dm,
            damping: &dc,
            stiffness: &dk,
        };
        // At w=1 the bare first structural diagonal is exactly zero.
        // The physical port load regularizes it; no artificial loss is added.
        for w in [0.7, 1.0, 1.000001, 1.4, 2.2] {
            let z = core::array::from_fn::<_, 4, _>(|i| C64::new(lc[i], lk[i] / w - w * lm[i]));
            let zw = core::array::from_fn::<_, 4, _>(|i| C64::new(0.0, -lk[i] / (w * w) - lm[i]));
            let zww = core::array::from_fn::<_, 4, _>(|i| C64::new(0.0, 2.0 * lk[i] / (w * w * w)));
            let compact = PreparedHarmonicResponse::try_diagonal_modal_with_port_impedance(
                &natural,
                &damping,
                &b,
                2,
                &z,
                w,
                controls(),
                port_limits(),
                cx,
            )
            .unwrap();
            let actual = compact
                .solve_with_loaded_frequency_parameter_derivative(
                    &force,
                    direction,
                    Some(&fw),
                    Some(&fww),
                    Some(&ft),
                    Some(&fwt),
                    &zw,
                    &zww,
                    limits(),
                    cx,
                )
                .unwrap();
            let dense = PreparedHarmonicResponse::try_new(
                &mass,
                &velocity,
                &stiffness,
                3,
                w,
                controls(),
                cx,
            )
            .unwrap()
            .solve_with_frequency_parameter_derivative(
                &force,
                direction,
                Some(&fw),
                Some(&fww),
                Some(&ft),
                Some(&fwt),
                limits(),
                cx,
            )
            .unwrap();
            let expected: Vec<_> = (0..3)
                .map(|i| Jet {
                    v: dense
                        .frequency_response
                        .frequency_response
                        .response
                        .displacement[i],
                    w: dense
                        .frequency_response
                        .frequency_response
                        .displacement_derivative[i],
                    ww: dense.frequency_response.displacement_second_derivative[i],
                    t: dense.parameter_response.displacement_derivative[i],
                    wt: dense.displacement_mixed_derivative[i],
                })
                .collect();
            check_physical_jet(&actual, &expected, w, 8e-11);
        }
    });
}

#[test]
fn loaded_mixed_admission_requires_complete_dispersion_and_preserves_the_owner() {
    with_cx(&CancelGate::new(), |cx| {
        let owner = PreparedHarmonicResponse::try_with_port_impedance(
            &[1.0],
            &[0.1],
            &[1.0],
            1,
            &[1.0],
            1,
            &[C64::from_re(0.2)],
            0.8,
            controls(),
            port_limits(),
            cx,
        )
        .unwrap();
        let direction = QuadraticParameterDirection {
            mass: &[0.1],
            damping: &[0.02],
            stiffness: &[0.3],
        };
        let zero = [C64::ZERO];
        let solve = |zw: &[C64], zww: &[C64], caps, context: &Cx<'_>| {
            owner.solve_with_loaded_frequency_parameter_derivative(
                &[C64::ONE],
                direction,
                None,
                None,
                None,
                None,
                zw,
                zww,
                caps,
                context,
            )
        };
        assert!(matches!(
            solve(&[], &zero, limits(), cx),
            Err(HarmonicResponseError::Invalid(_))
        ));
        assert!(matches!(
            solve(&zero, &[], limits(), cx),
            Err(HarmonicResponseError::Invalid(_))
        ));
        assert!(matches!(
            solve(&[C64::new(f64::NAN, 0.0)], &zero, limits(), cx),
            Err(HarmonicResponseError::Invalid(_)),
        ));
        assert!(matches!(
            solve(&zero, &[C64::new(0.0, f64::INFINITY)], limits(), cx),
            Err(HarmonicResponseError::Invalid(_)),
        ));
        assert!(matches!(
            solve(
                &zero,
                &zero,
                HarmonicFrequencyParameterControls {
                    maximum_additional_payload_bytes: 4615,
                    ..limits()
                },
                cx
            ),
            Err(HarmonicResponseError::Budget {
                requested: 4616,
                ..
            }),
        ));
        assert!(matches!(
            solve(
                &zero,
                &zero,
                HarmonicFrequencyParameterControls {
                    maximum_work: 6143,
                    ..limits()
                },
                cx
            ),
            Err(HarmonicResponseError::WorkBudget {
                requested: 6144,
                ..
            }),
        ));
        let gate = CancelGate::new();
        gate.request();
        with_cx(&gate, |cancelled| {
            assert!(matches!(
                solve(&zero, &zero, limits(), cancelled),
                Err(HarmonicResponseError::Cancelled),
            ));
        });
        let baseline = owner.solve(&[C64::ONE], cx).unwrap();
        let after = solve(&zero, &zero, limits(), cx).unwrap();
        assert_eq!(
            after
                .frequency_response
                .frequency_response
                .response
                .displacement,
            baseline.displacement,
        );
        let unloaded =
            PreparedHarmonicResponse::try_new(&[1.0], &[0.1], &[1.0], 1, 0.8, controls(), cx)
                .unwrap();
        assert!(matches!(
            unloaded.solve_with_loaded_frequency_parameter_derivative(
                &[C64::ONE],
                direction,
                None,
                None,
                None,
                None,
                &zero,
                &zero,
                limits(),
                cx,
            ),
            Err(HarmonicResponseError::Invalid(
                "loaded mixed frequency/parameter derivatives require an attached port load"
            )),
        ));
    });
}
