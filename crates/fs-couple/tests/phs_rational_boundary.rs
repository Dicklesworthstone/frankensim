//! Full original causal radiation-model states on the existing nonlinear owners.
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_couple::phs_orbit::{
    PeriodicPortInput, PhsOrbitError,
    rational::{PhsMimoError, PhsMimoLimits, PreparedPhsMimo},
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_math::c64::C64;
use fs_orbit::{
    HbBudget, HbScales, OrbitProblem, ShootBudget, VariationalShootControls, VariationalShootError,
};
use fs_phs::{PortHamiltonian, Storage};
use fs_vfit::{
    mimo::{
        MimoModel,
        load::{MimoLoadLimits, PreparedMimoLoad},
    },
    model::PoleTerm,
};
use std::{cell::Cell, rc::Rc};

const LIMITS: PhsMimoLimits = PhsMimoLimits {
    maximum_states: 32,
    maximum_ports: 8,
    maximum_harmonics: 31,
    maximum_hb_unknowns: 1024,
    maximum_products: 200_000_000,
    maximum_workspace_bytes: 4 << 20,
};
fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 419,
                    kernel_id: 97,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ),
            &gate,
        )
    });
}
fn source(model: &MimoModel) -> PreparedMimoLoad {
    PreparedMimoLoad::new(
        model,
        0.05,
        MimoLoadLimits {
            band_rad_s: [0.1, 100.0],
            sampled_passivity_points: 33,
            max_states: 32,
            max_memory_bytes: 1 << 20,
        },
    )
    .unwrap()
}
fn scalar(d: f64, dark: bool) -> MimoModel {
    MimoModel {
        ports: 1,
        poles: if dark {
            vec![
                PoleTerm::Real {
                    pole: -2.0,
                    residue: 0.0,
                },
                PoleTerm::Real {
                    pole: -3.0,
                    residue: 0.0,
                },
            ]
        } else {
            vec![PoleTerm::Real {
                pole: -3.0,
                residue: 0.0,
            }]
        },
        residues: if dark {
            vec![vec![C64::from_re(2.0)], vec![C64::ZERO]]
        } else {
            vec![vec![C64::from_re(2.0)]]
        },
        d: vec![d],
        e: vec![0.0],
        reciprocal: true,
    }
}
fn hb(harmonics: usize) -> HbBudget {
    HbBudget {
        harmonics,
        max_newton: 50,
        tolerance: 1e-11,
        deflate_equilibrium: false,
    }
}
fn scales(n: usize, w: f64) -> HbScales {
    HbScales {
        states: vec![1.0; n],
        equations: vec![1.0; n],
        angular_frequency: w,
    }
}
fn shoot_controls(n: usize, w: f64) -> VariationalShootControls {
    VariationalShootControls {
        state_scales: vec![1.0; n],
        period_scale: core::f64::consts::TAU / w,
        minimum_autonomous_excursion: 0.01,
        maximum_dimension: 32,
        maximum_payload_bytes: 8 << 20,
        maximum_work: 100_000_000_000,
        maximum_period_maps: 64,
    }
}
fn close(a: f64, b: f64, t: f64) {
    assert!((a - b).abs() < t * (1.0 + b.abs()), "{a:e} != {b:e}");
}
fn complex(a: C64, b: C64, t: f64) {
    assert!((a - b).abs() < t * (1.0 + b.abs()), "{a:?} != {b:?}");
}

#[test]
fn g1_full_mutual_original_real_pair_states_and_history_survive_loaded_hb() {
    with_cx(|cx, _| {
        let model = MimoModel {
            ports: 2,
            poles: vec![
                PoleTerm::Real {
                    pole: -2.0,
                    residue: 0.0,
                },
                PoleTerm::Pair {
                    pole: C64::new(-0.4, 1.8),
                    residue: C64::ZERO,
                },
            ],
            residues: vec![
                [0.8, -0.25, -0.25, 0.6].map(C64::from_re).to_vec(),
                [0.4, 0.12, 0.12, 0.3].map(C64::from_re).to_vec(),
            ],
            d: vec![0.9, 0.15, 0.15, 0.7],
            e: vec![0.0; 4],
            reciprocal: true,
        };
        let mut load = source(&model);
        load.step_mean_into(&[1.0, -0.4], &mut [0.0; 2]).unwrap();
        let history = load.pole_states().to_vec();
        let host = fs_phs::mass_spring_damper(2.0, 8.0, 0.3).unwrap();
        let w = 1.5;
        let f = C64::new(0.3, -0.1);
        let input = PeriodicPortInput::new(w, vec![vec![C64::from_re(0.4), f]]).unwrap();
        let prepared =
            PreparedPhsMimo::new(&host, &load, &input, &[0.0, 0.0, 1.0, -0.4], LIMITS, cx).unwrap();
        assert!(std::ptr::eq(prepared.source(), &load));
        assert!(std::ptr::eq(prepared.reports().0, load.reports().0));
        assert!(std::ptr::eq(prepared.reports().1, load.reports().1));
        assert_eq!(prepared.source_dimension(), 6);
        assert_eq!(prepared.dimension(), 8);
        let initial = prepared.initial_state(&[0.1, -0.2]).unwrap();
        assert_eq!(
            &initial[2..],
            load.continuous_operator().real_coordinates().as_slice()
        );
        assert!(prepared.time_domain_rhs().is_some());
        assert!(prepared.time_domain_state_jacobian().is_some());
        let result = prepared
            .solve_scaled_seeded(vec![0.0; 24], &hb(1), &scales(8, w))
            .unwrap();
        let s = C64::new(0.0, w);
        let map = [1.0, -0.4];
        let impedance = model.eval(s);
        let mut z = C64::ZERO;
        for a in 0..2 {
            for b in 0..2 {
                z = z + impedance[2 * a + b].scale(map[a] * map[b]);
            }
        }
        let q = f / (C64::from_re(8.0) + s.scale(0.3) + s * s.scale(2.0) + s * z);
        complex(result.orbit.coeffs[0][0], C64::from_re(0.05), 1e-10);
        complex(result.orbit.coeffs[0][1], q, 1e-10);
        complex(result.orbit.coeffs[1][1], s * q.scale(2.0), 1e-10);
        let velocity = s * q;
        let pair_den = (s + C64::from_re(0.4)) * (s + C64::from_re(0.4)) + C64::from_re(1.8 * 1.8);
        for p in 0..2 {
            let flow = velocity.scale(map[p]);
            complex(
                result.orbit.coeffs[2 + p][1],
                flow / (s + C64::from_re(2.0)),
                1e-10,
            );
            complex(
                result.orbit.coeffs[4 + p][1],
                flow * (s + C64::from_re(0.4)) / pair_den,
                1e-10,
            );
            complex(
                result.orbit.coeffs[6 + p][1],
                flow.scale(-1.8) / pair_den,
                1e-10,
            );
        }
        let audit = prepared.audit(&result.orbit, 128).unwrap();
        close(
            audit.mean_host_resistance_loss_w,
            0.6 * velocity.norm_sq(),
            1e-11,
        );
        close(
            audit.mean_source_exchange_power_w,
            2.0 * z.re * velocity.norm_sq(),
            1e-11,
        );
        assert!(audit.relative_power_residual < 1e-9);
        close(audit.mean_exchange_residual_w, 0.0, 1e-12);
        assert_eq!(load.pole_states(), history);
        let candidate = prepared.source_candidate(&initial).unwrap();
        assert_eq!(candidate.source_model(), load.source_model());
        assert_eq!(candidate.pole_states(), load.pole_states());
    });
}

fn amplification(pole: C64, period: f64, steps: usize) -> C64 {
    let z = pole.scale(period / steps as f64);
    let z2 = z * z;
    let step =
        C64::ONE + z + z2.scale(0.5) + (z2 * z).scale(1.0 / 6.0) + (z2 * z2).scale(1.0 / 24.0);
    (0..steps).fold(C64::ONE, |value, _| value * step)
}

#[test]
fn g1_shooting_retains_all_loaded_and_dark_floquet_modes() {
    with_cx(|cx, _| {
        let host = fs_phs::mass_spring_damper(1.0, 2.0, 0.25).unwrap();
        let model = scalar(0.75, true);
        let load = source(&model);
        let w = 1.2;
        let input =
            PeriodicPortInput::new(w, vec![vec![C64::from_re(0.2), C64::new(0.1, -0.04)]]).unwrap();
        let prepared = PreparedPhsMimo::new(&host, &load, &input, &[0.0, 1.0], LIMITS, cx).unwrap();
        let steps = 512;
        let result = prepared
            .solve_shooting_variational(
                &[0.0; 4],
                &ShootBudget {
                    steps_per_period: steps,
                    max_newton: 20,
                    tolerance: 1e-11,
                },
                &shoot_controls(4, w),
            )
            .unwrap();
        let s = C64::new(0.0, w);
        let q = C64::new(0.1, -0.04)
            / (C64::from_re(2.0) + s * s + s.scale(0.25) + s * model.eval(s)[0]);
        let p = s * q;
        let z = p / (s + C64::from_re(2.0));
        let dark = p / (s + C64::from_re(3.0));
        for (&actual, expected) in
            result
                .orbit
                .x0
                .iter()
                .zip([0.1 + 2.0 * q.re, 2.0 * p.re, 2.0 * z.re, 2.0 * dark.re])
        {
            close(actual, expected, 3e-8);
        }
        // The coupled characteristic is (s+1)(s^2+2s+4); the retained
        // zero-residue source pole contributes a separate multiplier at -3.
        assert_eq!(result.orbit.multipliers.len(), 4);
        for pole in [
            C64::from_re(-1.0),
            C64::new(-1.0, 3.0_f64.sqrt()),
            C64::new(-1.0, -3.0_f64.sqrt()),
            C64::from_re(-3.0),
        ] {
            let expected = amplification(pole, result.orbit.period, steps);
            let error = result
                .orbit
                .multipliers
                .iter()
                .map(|&(re, im)| (C64::new(re, im) - expected).abs())
                .fold(f64::INFINITY, f64::min);
            assert!(
                error < 2e-10,
                "missing {expected:?}: {:?}",
                result.orbit.multipliers
            );
        }
        assert_eq!(result.rhs_evaluations, 4 * steps * result.period_maps);
    });
}

struct Cubic {
    calls: Rc<Cell<usize>>,
}
impl Storage for Cubic {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        4.0 * x[0] * x[0] + 1.25 * x[0].powi(4) + 0.25 * x[1] * x[1]
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        self.calls.set(self.calls.get() + 1);
        out[0] = 8.0 * x[0] + 5.0 * x[0].powi(3);
        out[1] = x[1] / 2.0;
    }
}

#[test]
fn g1_nonlinear_storage_uses_full_hessian_and_original_boundary_memory() {
    with_cx(|cx, _| {
        let calls = Rc::new(Cell::new(0));
        let host = PortHamiltonian::new(
            2,
            1,
            vec![0.0, 1.0, -1.0, 0.0],
            vec![0.0, 0.0, 0.0, 0.6],
            vec![0.0, 1.0],
            Box::new(Cubic {
                calls: Rc::clone(&calls),
            }),
        )
        .unwrap();
        let model = scalar(0.4, false);
        let load = source(&model);
        let (w, a, map) = (2.0, 0.3, 0.7);
        let s = C64::new(0.0, w);
        let q = C64::from_re(a / 2.0);
        let f = q
            * (C64::from_re(8.0 - 2.0 * w * w + 0.75 * 5.0 * a * a)
                + s.scale(0.6)
                + s * model.eval(s)[0].scale(map * map));
        let input = PeriodicPortInput::new(
            w,
            vec![vec![C64::ZERO, f, C64::ZERO, (q * q * q).scale(5.0)]],
        )
        .unwrap();
        let prepared = PreparedPhsMimo::new(&host, &load, &input, &[0.0, map], LIMITS, cx).unwrap();
        let shooting = ShootBudget {
            steps_per_period: 512,
            max_newton: 20,
            tolerance: 1e-11,
        };
        assert!(matches!(
            prepared.solve_shooting_variational(&[0.0; 3], &shooting, &shoot_controls(3, w)),
            Err(PhsMimoError::Regime(PhsOrbitError::Variational(
                VariationalShootError::TangentUnavailable
            )))
        ));
        assert_eq!(calls.get(), 0);
        let hessian = |x: &[f64], direction: &[f64], out: &mut [f64]| {
            out[0] = (8.0 + 15.0 * x[0] * x[0]) * direction[0];
            out[1] = direction[1] / 2.0;
            true
        };
        let prepared = prepared.with_storage_hessian(&hessian);
        let mut jacobian = [0.0; 9];
        assert!(prepared.time_domain_state_jacobian().unwrap()(
            &prepared,
            0.0,
            &[0.2, 0.3, -0.1],
            &mut jacobian
        ));
        let expected = [
            0.0,
            0.5,
            0.0,
            -8.6,
            -(0.6 + map * map * 0.4) / 2.0,
            -2.0 * map,
            0.0,
            map / 2.0,
            -3.0,
        ];
        for (&actual, expected) in jacobian.iter().zip(expected) {
            close(actual, expected, 1e-13);
        }
        let response = prepared
            .solve_scaled_seeded(vec![0.0; 21], &hb(3), &scales(3, w))
            .unwrap();
        let z = s * q.scale(map) / (s + C64::from_re(3.0));
        complex(response.orbit.coeffs[0][1], q, 2e-10);
        complex(response.orbit.coeffs[1][1], s * q.scale(2.0), 2e-10);
        complex(response.orbit.coeffs[2][1], z, 2e-10);
        for row in &response.orbit.coeffs {
            assert!(row[2].abs() < 1e-10 && row[3].abs() < 1e-10);
        }
        let audit = prepared.audit(&response.orbit, 128).unwrap();
        close(
            audit.mean_host_stored_energy_j,
            0.25 * (8.0 + 2.0 * w * w) * a * a + 15.0 * a.powi(4) / 32.0,
            1e-10,
        );
        close(
            audit.mean_host_resistance_loss_w,
            0.3 * w * w * a * a,
            1e-10,
        );
        close(
            audit.mean_source_exchange_power_w,
            2.0 * model.eval(s)[0].re * (s * q.scale(map)).norm_sq(),
            1e-10,
        );
        assert!(audit.relative_power_residual < 1e-9);
        let shot = prepared
            .solve_shooting_variational(
                &[0.9 * a, 0.0, 2.0 * z.re],
                &shooting,
                &shoot_controls(3, w),
            )
            .unwrap();
        close(shot.orbit.x0[0], a, 3e-8);
        close(shot.orbit.x0[1], 0.0, 3e-8);
        close(shot.orbit.x0[2], 2.0 * z.re, 3e-8);
    });
}

#[test]
fn g3_dc_singular_direct_term_band_budget_and_cancellation_are_explicit() {
    with_cx(|cx, gate| {
        let host = PortHamiltonian::new(
            1,
            1,
            vec![0.0],
            vec![1.0],
            vec![1.0],
            Box::new(fs_phs::QuadraticStorage::new(vec![2.0], 1).unwrap()),
        )
        .unwrap();
        let model = scalar(0.0, false);
        let load = source(&model);
        let input = PeriodicPortInput::new(1.0, vec![vec![C64::from_re(1.0)]]).unwrap();
        let prepared = PreparedPhsMimo::new(&host, &load, &input, &[1.0], LIMITS, cx).unwrap();
        let response = prepared
            .solve_scaled_seeded(vec![0.0; 6], &hb(1), &scales(2, 1.0))
            .unwrap();
        // Actual original DC law: y=2x, z=y/3, 0=-y+1-2z.
        close(response.orbit.coeffs[0][0].re, 0.3, 1e-11);
        close(response.orbit.coeffs[1][0].re, 0.2, 1e-11);
        let audit = prepared.audit(&response.orbit, 32).unwrap();
        close(audit.mean_supplied_power_w, 0.6, 1e-11);
        close(audit.mean_host_resistance_loss_w, 0.36, 1e-11);
        close(audit.mean_source_exchange_power_w, 0.24, 1e-11);
        assert!(audit.relative_power_residual < 1e-10);
        let static_model = MimoModel {
            ports: 1,
            poles: vec![],
            residues: vec![],
            d: vec![0.5],
            e: vec![0.0],
            reciprocal: true,
        };
        let static_load = source(&static_model);
        let static_boundary =
            PreparedPhsMimo::new(&host, &static_load, &input, &[1.0], LIMITS, cx).unwrap();
        assert_eq!(static_boundary.dimension(), 1);
        let static_response = static_boundary
            .solve_scaled_seeded(vec![0.0; 3], &hb(1), &scales(1, 1.0))
            .unwrap();
        close(static_response.orbit.coeffs[0][0].re, 1.0 / 3.0, 1e-11);
        let high = PeriodicPortInput::new(60.0, vec![vec![C64::ZERO, C64::from_re(0.1)]]).unwrap();
        let high = PreparedPhsMimo::new(&host, &load, &high, &[1.0], LIMITS, cx).unwrap();
        assert!(matches!(
            high.solve_scaled_seeded(vec![0.0; 10], &hb(2), &scales(2, 60.0)),
            Err(PhsMimoError::Source(_))
        ));
        let tiny = PreparedPhsMimo::new(
            &host,
            &load,
            &input,
            &[1.0],
            PhsMimoLimits {
                maximum_products: 1,
                ..LIMITS
            },
            cx,
        )
        .unwrap();
        assert!(matches!(
            tiny.solve_scaled_seeded(vec![0.0; 6], &hb(1), &scales(2, 1.0)),
            Err(PhsMimoError::Budget(_))
        ));
        let mut improper = model.clone();
        improper.e[0] = 1.0;
        assert!(PreparedMimoLoad::new(&improper, 0.05, *load.limits()).is_err());
        gate.request();
        assert!(matches!(
            prepared.audit(&response.orbit, 32),
            Err(PhsMimoError::Cancelled)
        ));
        assert!(matches!(
            prepared.initial_state(&[0.0]),
            Err(PhsMimoError::Cancelled)
        ));
        assert_eq!(load.pole_states(), &[C64::ZERO]);
    });
}

#[test]
fn g1_original_host_and_causal_reaction_keep_one_joint_rounding() {
    with_cx(|cx, _| {
        let host = PortHamiltonian::new(
            1,
            1,
            vec![0.0],
            vec![1e20],
            vec![1.0],
            Box::new(fs_phs::QuadraticStorage::new(vec![1.0], 1).unwrap()),
        )
        .unwrap();
        let model = MimoModel {
            ports: 1,
            poles: vec![PoleTerm::Real {
                pole: -1.0,
                residue: 0.0,
            }],
            residues: vec![vec![C64::ONE]],
            d: vec![0.0],
            e: vec![0.0],
            reciprocal: true,
        };
        let load = source(&model);
        let input = PeriodicPortInput::new(1.0, vec![vec![C64::ONE]]).unwrap();
        let boundary = PreparedPhsMimo::new(&host, &load, &input, &[1.0], LIMITS, cx).unwrap();
        let mut rate = [0.0; 2];
        boundary.time_domain_rhs().unwrap()(&boundary, 0.0, &[1.0, -1e20], &mut rate);
        assert_eq!(rate[0], 1.0);
        assert!(rate[1].is_finite());

        // Complete native Hessian and direct-load derivative must retain the
        // same cancellation; applying a rounded host Jacobian first loses -1.
        let host = PortHamiltonian::new(
            2,
            1,
            vec![0.0, 1e20, -1e20, 0.0],
            vec![1e-20, 0.0, 0.0, 0.0],
            vec![0.0, 0.0],
            Box::new(fs_phs::QuadraticStorage::new(vec![1e20, 1.0, 1.0, 1.0], 2).unwrap()),
        )
        .unwrap();
        let model = MimoModel {
            ports: 1,
            poles: vec![],
            residues: vec![],
            d: vec![1.0],
            e: vec![0.0],
            reciprocal: true,
        };
        let load = source(&model);
        let input = PeriodicPortInput::new(1.0, vec![vec![C64::ZERO]]).unwrap();
        let boundary = PreparedPhsMimo::new(&host, &load, &input, &[1.0, 0.0], LIMITS, cx).unwrap();
        let mut jacobian = [0.0; 4];
        assert!(boundary.time_domain_state_jacobian().unwrap()(
            &boundary,
            0.0,
            &[1.0, 0.0],
            &mut jacobian
        ));
        assert_eq!(jacobian[0], -1.0);
        boundary.time_domain_rhs().unwrap()(&boundary, 0.0, &[1.0, 0.0], &mut rate);
        assert_eq!(rate[0], -1.0);
    });
}
