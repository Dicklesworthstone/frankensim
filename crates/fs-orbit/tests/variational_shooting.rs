//! Physical and analytic oracles for the complete discrete shooting tangent.
use core::f64::consts::TAU;
use fs_math::c64::C64;
use fs_orbit::{
    IslandStateJacobian, OrbitError, OrbitProblem, ShootBudget, TimeDomainRhs,
    TimeDomainStateJacobian, VariationalShootControls, VariationalShootError, solve_shooting,
    solve_shooting_variational,
};
use std::cell::Cell;

fn controls(scales: &[f64]) -> VariationalShootControls {
    VariationalShootControls {
        state_scales: scales.to_vec(),
        period_scale: TAU,
        minimum_autonomous_excursion: 0.1,
        maximum_dimension: 16,
        maximum_payload_bytes: 1 << 20,
        maximum_work: u64::MAX,
        maximum_period_maps: 96,
    }
}
fn budget(steps: usize) -> ShootBudget {
    ShootBudget {
        steps_per_period: steps,
        max_newton: 20,
        tolerance: 1.0e-11,
    }
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.17e}, expected={expected:.17e}, tolerance={tolerance:.3e}"
    );
}

// 2*x' = -x + scale*cos(t). The island Jacobian -1 is deliberately NOT
// the explicit flow Jacobian -1/2: an island tangent cannot replace the latter.
struct MassLoaded {
    scale: f64,
}
impl MassLoaded {
    fn flow(&self, t: f64, x: &[f64], out: &mut [f64]) {
        out[0] = 0.5 * (-x[0] + self.scale * t.cos());
    }
    fn flow_jacobian(&self, _: f64, _: &[f64], out: &mut [f64]) -> bool {
        out[0] = -0.5;
        true
    }
    fn island_jacobian(&self, _: f64, _: &[f64], out: &mut [f64]) -> bool {
        out[0] = -1.0;
        true
    }
}
impl OrbitProblem for MassLoaded {
    fn dim(&self) -> usize {
        1
    }
    fn island(&self, t: f64, x: &[f64], out: &mut [f64]) {
        out[0] = -x[0] + self.scale * t.cos();
    }
    fn port(&self, s: C64) -> Vec<C64> {
        vec![s * C64::from_re(2.0)]
    }
    fn island_state_jacobian(&self) -> Option<IslandStateJacobian<Self>> {
        Some(Self::island_jacobian)
    }
    fn time_domain_rhs(&self) -> Option<TimeDomainRhs<Self>> {
        Some(Self::flow)
    }
    fn time_domain_state_jacobian(&self) -> Option<TimeDomainStateJacobian<Self>> {
        Some(Self::flow_jacobian)
    }
    fn autonomous(&self) -> bool {
        false
    }
}

#[test]
fn mass_port_uses_original_flow_tangent_and_fourth_order_refinement() {
    let problem = MassLoaded { scale: 1.0 };
    let mut errors = Vec::new();
    for steps in [64, 128, 256] {
        let answer = solve_shooting_variational(
            &problem,
            &[0.0],
            TAU,
            &budget(steps),
            &controls(&[1.0]),
            || false,
        )
        .unwrap();
        errors.push((answer.orbit.x0[0] - 0.2).abs());
        let z = -0.5 * TAU / steps as f64;
        let amplification = 1.0 + z + z * z / 2.0 + z * z * z / 6.0 + z * z * z * z / 24.0;
        close(
            answer.monodromy[0],
            amplification.powi(steps as i32),
            2.0e-14,
        );
        close(answer.orbit.multipliers[0].0, answer.monodromy[0], 2.0e-14);
        close(answer.orbit.multipliers[0].1, 0.0, 0.0);
        assert!(answer.normalized_residual <= budget(steps).tolerance);
        assert_eq!(answer.rhs_evaluations, 4 * steps * answer.period_maps);
        assert!(answer.jacobian_evaluations < answer.rhs_evaluations);
        assert!(answer.period_derivative.is_none());
    }
    for pair in errors.windows(2) {
        assert!(
            pair[0] / pair[1] > 13.0 && pair[0] / pair[1] < 20.0,
            "{errors:?}"
        );
    }
    let old = solve_shooting(&problem, &[0.0], TAU, &budget(128)).unwrap();
    let exact = solve_shooting_variational(
        &problem,
        &[0.0],
        TAU,
        &budget(128),
        &controls(&[1.0]),
        || false,
    )
    .unwrap();
    close(old.x0[0], exact.orbit.x0[0], 1.0e-11);
    close(old.multipliers[0].0, (-core::f64::consts::PI).exp(), 1.0e-8);
}

#[test]
fn forced_response_and_floquet_are_covariant_in_original_state_units() {
    let reference = solve_shooting_variational(
        &MassLoaded { scale: 1.0 },
        &[0.0],
        TAU,
        &budget(128),
        &controls(&[1.0]),
        || false,
    )
    .unwrap();
    for scale in [1.0e-80, 1.0e80] {
        let answer = solve_shooting_variational(
            &MassLoaded { scale },
            &[0.0],
            TAU,
            &budget(128),
            &controls(&[scale]),
            || false,
        )
        .unwrap();
        close(answer.orbit.x0[0] / scale, reference.orbit.x0[0], 1.0e-14);
        close(answer.monodromy[0], reference.monodromy[0], 1.0e-14);
        assert!(answer.orbit.residual_trace.iter().all(|r| r.is_finite()));
    }
}

// Periodic, nonnormal linear flow with an analytic monodromy. Its exact zero
// forced response is valid and the off-diagonal entry retains the source phase.
struct PeriodicTriangular {
    scales: [f64; 2],
}
impl PeriodicTriangular {
    fn jacobian(&self, t: f64, _: &[f64], out: &mut [f64]) -> bool {
        out.copy_from_slice(&[
            -0.2,
            0.7 * t.cos() * (self.scales[0] / self.scales[1]),
            0.0,
            -0.45,
        ]);
        true
    }
}
impl OrbitProblem for PeriodicTriangular {
    fn dim(&self) -> usize {
        2
    }
    fn island(&self, t: f64, x: &[f64], out: &mut [f64]) {
        out[0] =
            self.scales[0] * (-0.2 * x[0] / self.scales[0] + 0.7 * t.cos() * x[1] / self.scales[1]);
        out[1] = -0.45 * x[1];
    }
    fn time_domain_rhs(&self) -> Option<TimeDomainRhs<Self>> {
        Some(Self::island)
    }
    fn time_domain_state_jacobian(&self) -> Option<TimeDomainStateJacobian<Self>> {
        Some(Self::jacobian)
    }
    fn autonomous(&self) -> bool {
        false
    }
}
#[test]
fn zero_forced_response_keeps_full_time_periodic_nonnormal_monodromy() {
    let a = (-0.2 * TAU).exp();
    let b = (-0.45 * TAU).exp();
    let off = 0.7 * (-0.25) * (b - a) / (1.0 + 0.25 * 0.25);
    for scales in [[1.0, 1.0], [1.0e-80, 1.0e80]] {
        let problem = PeriodicTriangular { scales };
        let answer = solve_shooting_variational(
            &problem,
            &[0.0, 0.0],
            TAU,
            &budget(512),
            &controls(&scales),
            || false,
        )
        .unwrap();
        assert_eq!(answer.orbit.x0, [0.0, 0.0]);
        assert_eq!(answer.period_maps, 1);
        assert_eq!(answer.maximum_scaled_excursion, 0.0);
        for (&x, y) in answer.monodromy.iter().zip([a, off, 0.0, b]) {
            close(x, y, 1.0e-10);
        }
    }
}

// A stable, nonisochronous Stuart-Landau cycle beside a static preload. The
// first state is constant on the cycle and cannot supply an autonomous phase.
struct Cycle {
    scales: [f64; 3],
    bias: f64,
}
impl Cycle {
    fn jacobian(&self, _: f64, x: &[f64], out: &mut [f64]) -> bool {
        let u = x[1] / self.scales[1];
        let v = x[2] / self.scales[2];
        let a = 0.15 * (1.0 - u * u - v * v);
        let w = 1.3 + 0.4 * (u * u + v * v - 1.0);
        let matrix = [
            -0.7,
            0.0,
            0.0,
            0.0,
            a - 0.3 * u * u - 0.8 * u * v,
            -0.3 * u * v - w - 0.8 * v * v,
            0.0,
            w + 0.8 * u * u - 0.3 * u * v,
            a - 0.3 * v * v + 0.8 * u * v,
        ];
        for i in 0..3 {
            for j in 0..3 {
                out[3 * i + j] = matrix[3 * i + j] * (self.scales[i] / self.scales[j]);
            }
        }
        true
    }
}
impl OrbitProblem for Cycle {
    fn dim(&self) -> usize {
        3
    }
    fn island(&self, _: f64, x: &[f64], out: &mut [f64]) {
        let u = x[1] / self.scales[1];
        let v = x[2] / self.scales[2];
        let a = 0.15 * (1.0 - u * u - v * v);
        let w = 1.3 + 0.4 * (u * u + v * v - 1.0);
        out[0] = -0.7 * (x[0] - self.bias);
        out[1] = self.scales[1] * (a * u - w * v);
        out[2] = self.scales[2] * (w * u + a * v);
    }
    fn time_domain_rhs(&self) -> Option<TimeDomainRhs<Self>> {
        Some(Self::island)
    }
    fn time_domain_state_jacobian(&self) -> Option<TimeDomainStateJacobian<Self>> {
        Some(Self::jacobian)
    }
    fn autonomous(&self) -> bool {
        true
    }
}

#[test]
fn autonomous_phase_uses_moving_full_state_despite_static_first_coordinate() {
    let period = TAU / 1.3;
    for (scales, angle) in [([1.0, 1.0, 1.0], 0.0_f64), ([1.0, 1.0e-12, 1.0e6], 0.7)] {
        let problem = Cycle {
            scales,
            bias: 1.0e12,
        };
        let seed = [
            problem.bias,
            1.08 * angle.cos() * scales[1],
            1.08 * angle.sin() * scales[2],
        ];
        let answer = solve_shooting_variational(
            &problem,
            &seed,
            period * 1.03,
            &budget(256),
            &controls(&scales),
            || false,
        )
        .unwrap();
        close(answer.orbit.period, period, 1.0e-7);
        assert_eq!(answer.orbit.x0[0], problem.bias);
        let u = answer.orbit.x0[1] / scales[1];
        let v = answer.orbit.x0[2] / scales[2];
        close(u.hypot(v), 1.0, 5.0e-8);
        assert!(answer.normalized_residual <= 1.0e-11);
        assert!(answer.phase_residual <= 1.0e-11);
        assert!(answer.maximum_scaled_excursion > 1.0);
        let mut multipliers = answer.orbit.multipliers.clone();
        multipliers.sort_by(|a, b| a.0.total_cmp(&b.0));
        for (actual, expected) in
            multipliers
                .iter()
                .zip([(-0.7 * period).exp(), (-0.3 * period).exp(), 1.0])
        {
            close(actual.0, expected, 2.0e-7);
            close(actual.1, 0.0, 1.0e-12);
        }
        let mut flow = [0.0; 3];
        problem.island(0.0, &answer.orbit.x0, &mut flow);
        for (i, &actual) in answer
            .period_derivative
            .as_ref()
            .unwrap()
            .iter()
            .enumerate()
        {
            close(actual / scales[i], flow[i] / scales[i], 3.0e-7);
        }
        assert_eq!(answer.rhs_evaluations, 1 + 4 * 256 * answer.period_maps);
    }
}

// Manufactured nonlinear forced oscillator q=cos(t), v=-sin(t), with the
// exact cubic stiffness retained in both original flow and state derivative.
struct Duffing;
impl Duffing {
    fn jacobian(&self, _: f64, x: &[f64], out: &mut [f64]) -> bool {
        out.copy_from_slice(&[0.0, 1.0, -1.7 - 1.2 * x[0] * x[0], -0.3]);
        true
    }
}
impl OrbitProblem for Duffing {
    fn dim(&self) -> usize {
        2
    }
    fn island(&self, t: f64, x: &[f64], out: &mut [f64]) {
        let c = t.cos();
        let force = 0.7 * c - 0.3 * t.sin() + 0.4 * c * c * c;
        out[0] = x[1];
        out[1] = -1.7 * x[0] - 0.3 * x[1] - 0.4 * x[0] * x[0] * x[0] + force;
    }
    fn time_domain_rhs(&self) -> Option<TimeDomainRhs<Self>> {
        Some(Self::island)
    }
    fn time_domain_state_jacobian(&self) -> Option<TimeDomainStateJacobian<Self>> {
        Some(Self::jacobian)
    }
    fn autonomous(&self) -> bool {
        false
    }
}
#[test]
fn nonlinear_forced_turning_point_and_floquet_volume_match_original_law() {
    let mut errors = Vec::new();
    let mut volume_errors = Vec::new();
    for steps in [128, 256, 512] {
        let answer = solve_shooting_variational(
            &Duffing,
            &[0.8, 0.2],
            TAU,
            &budget(steps),
            &controls(&[1.0, 1.0]),
            || false,
        )
        .unwrap();
        errors.push((answer.orbit.x0[0] - 1.0).hypot(answer.orbit.x0[1]));
        let m = &answer.monodromy;
        // Liouville's formula: exact continuous determinant = exp(integral tr J).
        volume_errors.push((m[0] * m[3] - m[1] * m[2] - (-0.3 * TAU).exp()).abs());
        assert!(
            answer
                .orbit
                .multipliers
                .iter()
                .all(|&(r, i)| r.hypot(i) < 1.0)
        );
    }
    for measured in [&errors, &volume_errors] {
        for pair in measured.windows(2) {
            assert!(
                pair[0] / pair[1] > 13.0 && pair[0] / pair[1] < 20.0,
                "{measured:?}"
            );
        }
        assert!(measured[2] < 2.0e-9, "{measured:?}");
    }
}

#[derive(Clone, Copy)]
enum Fault {
    None,
    NoRhs,
    NoJacobian,
    RefusedJacobian,
    UnwrittenJacobian,
    UnwrittenRhs,
}
struct Probe {
    fault: Fault,
    calls: Cell<usize>,
    auto: bool,
    drift: bool,
}
impl Probe {
    fn new(fault: Fault) -> Self {
        Self {
            fault,
            calls: Cell::new(0),
            auto: false,
            drift: false,
        }
    }
    fn jacobian(&self, _: f64, _: &[f64], out: &mut [f64]) -> bool {
        self.calls.set(self.calls.get() + 1);
        if matches!(self.fault, Fault::RefusedJacobian) {
            return false;
        }
        if !matches!(self.fault, Fault::UnwrittenJacobian) {
            out[0] = if self.drift { 0.0 } else { -1.0 };
        }
        true
    }
}
impl OrbitProblem for Probe {
    fn dim(&self) -> usize {
        1
    }
    fn island(&self, _: f64, x: &[f64], out: &mut [f64]) {
        self.calls.set(self.calls.get() + 1);
        if !matches!(self.fault, Fault::UnwrittenRhs) {
            out[0] = if self.drift { 1.0 } else { -x[0] };
        }
    }
    fn time_domain_rhs(&self) -> Option<TimeDomainRhs<Self>> {
        if matches!(self.fault, Fault::NoRhs) {
            None
        } else {
            Some(Self::island)
        }
    }
    fn time_domain_state_jacobian(&self) -> Option<TimeDomainStateJacobian<Self>> {
        if matches!(self.fault, Fault::NoJacobian) {
            None
        } else {
            Some(Self::jacobian)
        }
    }
    fn autonomous(&self) -> bool {
        self.auto
    }
}

#[test]
fn unavailable_or_refused_exact_sources_never_fall_back_to_trajectory_probes() {
    for (fault, expected) in [
        (Fault::NoRhs, VariationalShootError::TimeDomainUnavailable),
        (Fault::NoJacobian, VariationalShootError::TangentUnavailable),
    ] {
        let problem = Probe::new(fault);
        assert_eq!(
            solve_shooting_variational(
                &problem,
                &[1.0],
                TAU,
                &budget(16),
                &controls(&[1.0]),
                || false
            )
            .unwrap_err(),
            expected
        );
        assert_eq!(problem.calls.get(), 0);
    }
    for (fault, tangent, calls) in [
        (Fault::RefusedJacobian, true, 2),
        (Fault::UnwrittenJacobian, true, 2),
        (Fault::UnwrittenRhs, false, 1),
    ] {
        let problem = Probe::new(fault);
        assert_eq!(
            solve_shooting_variational(
                &problem,
                &[1.0],
                TAU,
                &budget(16),
                &controls(&[1.0]),
                || false
            )
            .unwrap_err(),
            VariationalShootError::SourceRefused { time: 0.0, tangent }
        );
        assert_eq!(problem.calls.get(), calls);
    }
}

#[test]
fn resource_admission_and_autonomous_excursion_precede_false_cycle_acceptance() {
    let problem = Probe::new(Fault::None);
    let mut config = controls(&[1.0]);
    config.maximum_payload_bytes = 1;
    assert!(matches!(
        solve_shooting_variational(&problem, &[1.0], TAU, &budget(16), &config, || false),
        Err(VariationalShootError::Budget { .. })
    ));
    config = controls(&[1.0]);
    config.maximum_work = 1;
    assert!(matches!(
        solve_shooting_variational(&problem, &[1.0], TAU, &budget(16), &config, || false),
        Err(VariationalShootError::Budget { .. })
    ));
    config = controls(&[0.0]);
    assert!(matches!(
        solve_shooting_variational(&problem, &[1.0], TAU, &budget(16), &config, || false),
        Err(VariationalShootError::Invalid(_))
    ));
    config = controls(&[1.0]);
    config.maximum_period_maps = usize::MAX;
    assert!(matches!(
        solve_shooting_variational(&problem, &[1.0], TAU, &budget(usize::MAX), &config, || {
            false
        }),
        Err(VariationalShootError::Invalid(_))
    ));
    assert_eq!(problem.calls.get(), 0);

    let drift = Probe {
        auto: true,
        drift: true,
        ..Probe::new(Fault::None)
    };
    let mut config = controls(&[1.0]);
    config.minimum_autonomous_excursion = 1.0e-12;
    assert!(matches!(
        solve_shooting_variational(
            &drift,
            &[0.0],
            1.0e-10,
            &ShootBudget {
                tolerance: 1.0e-8,
                ..budget(16)
            },
            &config,
            || false
        ),
        Err(VariationalShootError::Invalid(_))
    ));
    assert_eq!(drift.calls.get(), 0);
    config.minimum_autonomous_excursion = 0.1;
    assert_eq!(
        solve_shooting_variational(
            &drift,
            &[0.0],
            1.0e-10,
            &ShootBudget {
                tolerance: 1.0e-8,
                ..budget(16)
            },
            &config,
            || false
        )
        .unwrap_err(),
        VariationalShootError::Numerical(OrbitError::TrivialCollapse)
    );
    let equilibrium = Probe {
        auto: true,
        ..Probe::new(Fault::None)
    };
    assert_eq!(
        solve_shooting_variational(&equilibrium, &[0.0], TAU, &budget(16), &config, || false)
            .unwrap_err(),
        VariationalShootError::Numerical(OrbitError::TrivialCollapse)
    );
}

#[test]
fn cancellation_and_map_budget_are_atomic_and_replay_is_deterministic() {
    let problem = Probe::new(Fault::None);
    let mut polls = 0;
    assert_eq!(
        solve_shooting_variational(
            &problem,
            &[1.0],
            TAU,
            &budget(16),
            &controls(&[1.0]),
            || {
                polls += 1;
                polls >= 25
            }
        )
        .unwrap_err(),
        VariationalShootError::Cancelled
    );
    assert!(problem.calls.get() > 0);
    let mut config = controls(&[1.0]);
    config.maximum_period_maps = 1;
    assert!(matches!(
        solve_shooting_variational(&problem, &[1.0], TAU, &budget(16), &config, || false),
        Err(VariationalShootError::Budget {
            quantity: "period maps",
            requested: 2,
            cap: 1
        })
    ));
    let one = solve_shooting_variational(
        &problem,
        &[1.0],
        TAU,
        &budget(16),
        &controls(&[1.0]),
        || false,
    )
    .unwrap();
    let two = solve_shooting_variational(
        &problem,
        &[1.0],
        TAU,
        &budget(16),
        &controls(&[1.0]),
        || false,
    )
    .unwrap();
    assert_eq!(one.orbit.x0, two.orbit.x0);
    assert_eq!(one.orbit.residual_trace, two.orbit.residual_trace);
    assert_eq!(one.orbit.multipliers, two.orbit.multipliers);
    assert_eq!(one.monodromy, two.monodromy);
    assert_eq!(one.period_maps, two.period_maps);
}
