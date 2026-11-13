//! G1/G3/G4: complete forced branches, physical outputs and numerical stability.

use core::cell::Cell;
use std::rc::Rc;

use fs_couple::phs_orbit::{
    PeriodicPortInput, PhsFloquetControls, PhsForcedBranchStabilityControls,
    PhsForcedPortBranchControls, PhsForcedPortBranchError, PhsNumericalStability, PhsOrbitError,
    continue_forced_port_response_analytic_scaled_seeded,
    continue_forced_port_response_scaled_seeded, spectrum::PortSpectrumControls,
};
use fs_math::{c64::C64, det};
use fs_orbit::{
    ContinuationBudget, HbAnchor, HbBudget, HbOrbit, HbScales, OrbitError, ShootBudget,
    VariationalShootControls,
};
use fs_phs::{PortHamiltonian, RelaxationBranch, Storage};

fn controls(n: usize, m: usize, h: usize, output: usize) -> PhsForcedPortBranchControls {
    PhsForcedPortBranchControls {
        continuation: ContinuationBudget {
            max_steps: 500,
            initial_step: 0.02,
            min_step: 1e-6,
            hb: HbBudget {
                harmonics: h,
                max_newton: 40,
                tolerance: 1e-10,
                deflate_equilibrium: false,
            },
        },
        harmonic_balance_scales: HbScales {
            states: vec![1.0; n],
            equations: vec![1.0; n],
            angular_frequency: 1.0,
        },
        audit_samples: 256,
        spectrum: PortSpectrumControls {
            harmonics: output,
            initial_samples: 64,
            max_samples: 256,
            max_source_evaluations: 448,
            max_workspace_bytes: 4 << 20,
            rms_tolerances: vec![2e-8; m],
        },
        stability: None,
        maximum_retained_bytes: 16 << 20,
    }
}

fn stability(n: usize) -> PhsForcedBranchStabilityControls {
    PhsForcedBranchStabilityControls {
        floquet: PhsFloquetControls {
            shooting: ShootBudget {
                steps_per_period: 256,
                max_newton: 20,
                tolerance: 1e-11,
            },
            state_tolerances: vec![5e-6; n],
            spectral_radius_tolerance: 1e-6,
            unit_circle_guard: 1e-3,
        },
        variational: VariationalShootControls {
            state_scales: vec![1.0; n],
            period_scale: 1.0,
            minimum_autonomous_excursion: 0.01,
            maximum_dimension: n,
            maximum_payload_bytes: 2_000_000,
            maximum_work: 10_000_000_000,
            maximum_period_maps: 64,
        },
    }
}

fn near(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
        "actual={actual:?}, expected={expected:?}"
    );
}

#[test]
fn original_maxwell_branch_retains_phases_memories_work_and_all_decay_directions() {
    // The unforced characteristic polynomial is (s+1)(s+2)(s+3).
    let (k, arm, tau, damping) = (1.5, 1.5_f64, 0.25, 2.0);
    let system = fs_phs::mass_spring_damper(1.0, k, damping)
        .unwrap()
        .with_relaxation_branches(vec![RelaxationBranch {
            projection: vec![1.0, 0.0],
            stiffness: arm,
            relaxation_time_s: tau,
        }])
        .unwrap();
    let forcing = vec![
        C64::from_re(0.02),
        C64::new(0.2, -0.1),
        C64::new(0.02, 0.01),
        C64::new(-0.007, 0.003),
    ];
    let input = PeriodicPortInput::new(4.5, vec![forcing.clone()]).unwrap();
    let mut c = controls(3, 1, 3, 3);
    c.harmonic_balance_scales.angular_frequency = 5.0;
    c.stability = Some(stability(3));
    let branch = continue_forced_port_response_scaled_seeded(
        &system,
        &input,
        7.0,
        vec![0.0; 21],
        &c,
        || false,
    )
    .unwrap();
    assert!(branch.len() > 4);
    assert_eq!(branch[0].response.point.point.lambda, 4.5);
    assert!(branch.last().unwrap().response.point.point.lambda >= 7.0);
    for point in &branch {
        let orbit = &point.response.point.point.orbit;
        assert_eq!(orbit.coeffs.len(), 3);
        assert_eq!(point.ports.omega_rad_s.to_bits(), orbit.omega.to_bits());
        let mut heat = 0.0;
        for (h, &force) in forcing.iter().enumerate() {
            let w = h as f64 * orbit.omega;
            let s = C64::new(0.0, w);
            let relaxation = C64::ONE + s.scale(tau);
            let q =
                force / (C64::new(k - w * w, damping * w) + (s.scale(tau) / relaxation).scale(arm));
            let p = s * q;
            let z = q.scale(det::sqrt(arm)) / relaxation;
            near(orbit.coeffs[0][h], q, 2e-9);
            near(orbit.coeffs[1][h], p, 2e-9);
            near(orbit.coeffs[2][h], z, 2e-9);
            if h > 0 {
                near(
                    point.ports.peak_phasors[0][h - 1],
                    p.conj().scale(2.0),
                    5e-9,
                );
                heat += 2.0 * damping * p.norm_sq()
                    + 2.0 * (q.scale(det::sqrt(arm)) - z).norm_sq() / tau;
            } else {
                assert!(point.ports.dc[0].abs() < 1e-9);
            }
        }
        assert!((point.response.audit.mean_resistance_loss_w - heat).abs() < 2e-8);
        assert!(point.response.audit.relative_power_residual < 2e-8);
        let resolved = point.stability.as_ref().unwrap();
        assert_eq!(resolved.classification, PhsNumericalStability::Stable);
        assert_eq!(resolved.levels.len(), 3);
        assert!(resolved.dominant_growth_rate_per_s.unwrap() < 0.0);
        for level in &resolved.levels {
            assert_eq!(level.orbit.orbit.multipliers.len(), 3);
            assert_eq!(level.orbit.monodromy.len(), 9);
            assert_eq!(level.orbit.orbit.x0.len(), 3);
            let period = core::f64::consts::TAU / orbit.omega;
            assert_eq!(level.orbit.orbit.period.to_bits(), period.to_bits());
            for pole in [-1.0, -2.0, -3.0] {
                let exact = det::exp(pole * period);
                assert!(
                    level
                        .orbit
                        .orbit
                        .multipliers
                        .iter()
                        .any(|&(re, im)| { (C64::new(re, im) - C64::from_re(exact)).abs() < 2e-6 })
                );
            }
        }
    }
    assert_eq!(input.omega_rad_s(), 4.5);
}

struct CubicStorage {
    cubic: f64,
    calls: Rc<Cell<usize>>,
}
impl Storage for CubicStorage {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        0.5 * (x[0] * x[0] + x[1] * x[1]) + 0.25 * self.cubic * x[0].powi(4)
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        self.calls.set(self.calls.get() + 1);
        out[0] = x[0] + self.cubic * x[0].powi(3);
        out[1] = x[1];
    }
}

fn cubic_system(cubic: f64, damping: f64) -> (PortHamiltonian, Rc<Cell<usize>>) {
    let calls = Rc::new(Cell::new(0));
    let system = PortHamiltonian::new(
        2,
        2,
        vec![0.0, 1.0, -1.0, 0.0],
        vec![0.0, 0.0, 0.0, damping],
        vec![1.0, 0.0, 0.0, 1.0],
        Box::new(CubicStorage {
            cubic,
            calls: calls.clone(),
        }),
    )
    .unwrap();
    (system, calls)
}

fn signed(row: &[C64], harmonic: isize) -> C64 {
    if harmonic.unsigned_abs() >= row.len() {
        return C64::ZERO;
    }
    if harmonic < 0 {
        row[-harmonic as usize].conj()
    } else {
        row[harmonic as usize]
    }
}

fn cube_coefficient(row: &[C64], harmonic: isize) -> C64 {
    // Independent finite Fourier convolution of the ORIGINAL cubic effort.
    let h = row.len() as isize - 1;
    let mut result = C64::ZERO;
    for a in -h..=h {
        for b in -h..=h {
            result = result + signed(row, a) * signed(row, b) * signed(row, harmonic - a - b);
        }
    }
    result
}

#[test]
fn folded_branch_keeps_nonlinear_port_effort_and_the_complete_memory_hessian() {
    let (cubic, damping, arm, tau) = (0.15, 0.1, 0.005_f64, 0.2);
    let (base, _) = cubic_system(cubic, damping);
    let system = base
        .with_relaxation_branches(vec![RelaxationBranch {
            projection: vec![1.0, 0.0],
            stiffness: arm,
            relaxation_time_s: tau,
        }])
        .unwrap();
    let input = PeriodicPortInput::new(
        0.6,
        vec![
            vec![C64::ZERO, C64::ZERO],
            vec![C64::ZERO, C64::from_re(0.15)],
        ],
    )
    .unwrap();
    let c = controls(3, 2, 9, 27);
    let hessian_calls = Cell::new(0);
    let hessian = |x: &[f64], v: &[f64], out: &mut [f64]| {
        hessian_calls.set(hessian_calls.get() + 1);
        let coupling = det::sqrt(arm);
        out[0] = (1.0 + 3.0 * cubic * x[0] * x[0] + arm) * v[0] - coupling * v[2];
        out[1] = v[1];
        out[2] = v[2] - coupling * v[0];
        true
    };
    let branch = continue_forced_port_response_analytic_scaled_seeded(
        &system,
        &input,
        1.8,
        vec![0.0; 57],
        &c,
        &hessian,
        || false,
    )
    .unwrap();
    assert!(hessian_calls.get() > 0, "the complete Hessian was bypassed");
    let reversals = branch
        .windows(3)
        .filter(|p| {
            (p[1].ports.omega_rad_s - p[0].ports.omega_rad_s)
                * (p[2].ports.omega_rad_s - p[1].ports.omega_rad_s)
                < 0.0
        })
        .count();
    assert_eq!(
        reversals, 2,
        "the folded branch was sorted or its fold was lost"
    );
    assert!(branch.last().unwrap().ports.omega_rad_s >= 1.8);
    let mut nonlinear_difference = 0.0_f64;
    let mut third_harmonic = 0.0_f64;
    for point in &branch {
        let orbit = &point.response.point.point.orbit;
        assert_eq!(orbit.coeffs.len(), 3);
        assert!(point.stability.is_none());
        assert!(point.response.point.scaled_residual < c.continuation.hb.tolerance);
        assert!(point.response.audit.relative_power_residual < 3e-8);
        for h in [1, 3, 5] {
            let restoring = orbit.coeffs[0][h].scale(1.0 + arm)
                + cube_coefficient(&orbit.coeffs[0], h as isize).scale(cubic)
                - orbit.coeffs[2][h].scale(det::sqrt(arm));
            near(
                point.ports.peak_phasors[0][h - 1],
                restoring.conj().scale(2.0),
                5e-9,
            );
            near(
                point.ports.peak_phasors[1][h - 1],
                orbit.coeffs[1][h].conj().scale(2.0),
                5e-9,
            );
        }
        nonlinear_difference = nonlinear_difference
            .max((point.ports.peak_phasors[0][0] - orbit.coeffs[0][1].conj().scale(2.0)).abs());
        third_harmonic = third_harmonic.max(orbit.coeffs[0][3].abs());
        assert!(orbit.coeffs[2][1].abs() > 0.0);
    }
    assert!(nonlinear_difference > 0.2);
    assert!(third_harmonic > 1e-3);
}

#[test]
fn nonlinear_branch_classifies_the_original_parametric_instability_at_its_clock() {
    let (system, _) = cubic_system(1.0, 0.05);
    let hessian = |x: &[f64], v: &[f64], out: &mut [f64]| {
        out[0] = (1.0 + 3.0 * x[0] * x[0]) * v[0];
        out[1] = v[1];
        true
    };
    for (omega, expected, radius) in [
        (
            det::sqrt(2.5),
            PhsNumericalStability::Unstable,
            2.301_406_300_465_6,
        ),
        (3.0, PhsNumericalStability::Stable, 0.948_987_286_154),
    ] {
        // Exact original q=cos(theta), p=-omega*sin(theta). Its variational
        // equation is u''+.05u'+[2.5+1.5cos(2*omega*t)]u=0. The independent
        // midpoint matrix-exponential ladder gives the radii above; the first
        // frequency lies in its primary parametric-instability tongue.
        let mut rows = vec![vec![C64::ZERO; 8]; 2];
        rows[0][1] = C64::from_re(0.5);
        rows[1][1] = C64::new(0.0, 0.5 * omega);
        let seed = HbOrbit {
            omega,
            coeffs: rows,
            residual_trace: vec![],
            residual: 0.0,
            aft_samples: 0,
            sampling_defect: 0.0,
        }
        .packed_seed(&HbAnchor::Forced { omega });
        let input = PeriodicPortInput::new(
            omega,
            vec![
                vec![C64::ZERO; 4],
                vec![
                    C64::ZERO,
                    C64::new(0.5 * (1.75 - omega * omega), 0.025 * omega),
                    C64::ZERO,
                    C64::from_re(0.125),
                ],
            ],
        )
        .unwrap();
        let mut c = controls(2, 2, 7, 21);
        c.continuation.initial_step = 0.002;
        c.harmonic_balance_scales.angular_frequency = omega;
        c.stability = Some(stability(2));
        let branch = continue_forced_port_response_analytic_scaled_seeded(
            &system,
            &input,
            1.005 * omega,
            seed,
            &c,
            &hessian,
            || false,
        )
        .unwrap();
        assert!(branch.len() >= 2);
        let first = &branch[0];
        near(first.ports.peak_phasors[0][0], C64::from_re(1.75), 2e-9);
        near(first.ports.peak_phasors[0][2], C64::from_re(0.25), 2e-9);
        for point in &branch {
            let result = point.stability.as_ref().unwrap();
            assert_eq!(result.classification, expected);
            assert_eq!(
                result.dominant_growth_rate_per_s.unwrap() > 0.0,
                expected == PhsNumericalStability::Unstable
            );
        }
        for level in &first.stability.as_ref().unwrap().levels {
            assert!((level.spectral_radius - radius).abs() < 2e-6);
            let m = &level.orbit.monodromy;
            let determinant = m[0] * m[3] - m[1] * m[2];
            assert!((determinant - det::exp(-0.05 * level.orbit.orbit.period)).abs() < 2e-7);
        }
    }
}

#[test]
fn branch_refuses_missing_tangents_budget_and_cancellation_without_partial_results() {
    let (system, calls) = cubic_system(0.15, 0.1);
    let input = PeriodicPortInput::new(
        0.6,
        vec![
            vec![C64::ZERO, C64::ZERO],
            vec![C64::ZERO, C64::from_re(0.15)],
        ],
    )
    .unwrap();
    let mut c = controls(2, 2, 5, 15);
    let seed = vec![0.0; 22];
    let hessian = |x: &[f64], v: &[f64], out: &mut [f64]| {
        out[0] = (1.0 + 0.45 * x[0] * x[0]) * v[0];
        out[1] = v[1];
        true
    };
    calls.set(0);
    assert!(matches!(
        continue_forced_port_response_scaled_seeded(
            &system,
            &input,
            0.65,
            seed.clone(),
            &c,
            || false,
        ),
        Err(PhsForcedPortBranchError::Invalid(_))
    ));
    assert_eq!(calls.get(), 0);
    c.maximum_retained_bytes = 1;
    assert!(matches!(
        continue_forced_port_response_analytic_scaled_seeded(
            &system,
            &input,
            0.65,
            seed.clone(),
            &c,
            &hessian,
            || false,
        ),
        Err(PhsForcedPortBranchError::Budget { .. })
    ));
    assert_eq!(calls.get(), 0);
    c.maximum_retained_bytes = 16 << 20;
    assert_eq!(
        continue_forced_port_response_analytic_scaled_seeded(
            &system,
            &input,
            0.65,
            seed.clone(),
            &c,
            &hessian,
            || true,
        )
        .unwrap_err(),
        PhsForcedPortBranchError::Cancelled
    );
    assert_eq!(calls.get(), 0);
    let mut polls = 0;
    assert_eq!(
        continue_forced_port_response_analytic_scaled_seeded(
            &system,
            &input,
            0.65,
            seed.clone(),
            &c,
            &hessian,
            || {
                polls += 1;
                polls > 20
            },
        )
        .unwrap_err(),
        PhsForcedPortBranchError::Cancelled
    );
    assert!(
        calls.get() > 0,
        "cancellation was not exercised within continuation"
    );
    let refusing = |_: &[f64], _: &[f64], _: &mut [f64]| false;
    assert!(matches!(
        continue_forced_port_response_analytic_scaled_seeded(
            &system,
            &input,
            0.65,
            seed.clone(),
            &c,
            &refusing,
            || false,
        ),
        Err(PhsForcedPortBranchError::Continuation(
            PhsOrbitError::Orbit(OrbitError::TangentRefused { .. })
        ))
    ));
    let replay = || {
        continue_forced_port_response_analytic_scaled_seeded(
            &system,
            &input,
            0.65,
            seed.clone(),
            &c,
            &hessian,
            || false,
        )
        .unwrap()
    };
    let a = replay();
    let b = replay();
    assert_eq!(a.len(), b.len());
    for (a, b) in a.iter().zip(&b) {
        assert_eq!(
            a.response.point.point.orbit.coeffs,
            b.response.point.point.orbit.coeffs
        );
        assert_eq!(a.ports.peak_phasors, b.ports.peak_phasors);
    }
    assert_eq!(input.omega_rad_s(), 0.6);
}
