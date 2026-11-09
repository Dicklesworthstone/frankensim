//! G1/G3/G4: continuously evolving friction heat, finite cooling, and rollback.
use super::*;
use fs_couple::friction_heat::hertz::march::*;
use fs_mbd::tangential_friction::hertz::transient::HertzTransientControls;
use fs_tribo::sliding_patch::{HertzSlidingPatch, HertzSlidingPatchParameters};

fn march_controls(scale: f64) -> HertzHeatMarchControls {
    HertzHeatMarchControls {
        initial_step_s: 0.3,
        minimum_step_s: 1e-7,
        maximum_step_s: 0.3,
        maximum_attempts: 4_000,
        maximum_accepted_phases: 8_000,
        maximum_thermal_trials: 400_000,
        tolerance: HertzHeatMarchTolerance {
            linear_momentum_ns: 2e-5 * scale,
            angular_momentum_nms: 2e-6 * scale,
            cell_energy_j: 2e-6 * scale,
            scaled_travel_m: 2e-6 * scale,
            dissipated_j: 2e-6 * scale,
            relative: 0.0,
        },
        energy_tolerance_j: 5e-8,
    }
}
fn path_controls() -> HertzTransientControls {
    HertzTransientControls {
        maximum_step_s: 0.02,
        minimum_step_s: 1e-15,
        maximum_attempts: 100_000,
        maximum_evaluations: 400_001,
        scaled_rate_tolerance_m_s: 1e-11,
        scaled_impulse_tolerance_ns: 1e-10,
        scaled_travel_tolerance_m: 1e-10,
        heat_tolerance_j: 1e-12,
        stopping_time_tolerance_s: 1e-7,
    }
}
fn run_march<'m>(
    pair: &Pair,
    owner: &SourceThermalCoulomb,
    thermal: &DepositionSimulation<'m>,
    route: &FrictionHeatRoute,
    duration: f64,
    model: HertzHeatMarchModel,
    march: HertzHeatMarchControls,
    resolve: impl Fn(&DepositionSimulation<'m>, f64) -> Result<RigidHertzHeatSource<f64>, &'static str>,
    gate: &CancelGate,
) -> Result<HeatedRigidHertzMarchCandidate<'m, f64>, HertzHeatMarchError<&'static str>> {
    propose_rigid_hertz_heat_march(
        &pair.bodies,
        &pair.properties,
        pair.contact(owner),
        thermal,
        route,
        duration,
        controls(),
        model,
        march,
        300.0,
        resolve,
        gate,
    )
}

#[test]
fn refined_original_phases_resolve_continuous_fade_arrest_and_rest_heat() {
    let pair = Pair::new(0.5, 0.0);
    let owner = law(0);
    let thermal = thermal([300.0; 3], So3::identity());
    let route = routes(&thermal)[0].clone();
    let reference = continuous_reference(0.3);
    let mut errors = Vec::new();
    for scale in [1.0, 0.0625] {
        let candidate = run_march(
            &pair,
            &owner,
            &thermal,
            &route,
            0.3,
            HertzHeatMarchModel::InvariantDirection,
            march_controls(scale),
            source,
            &CancelGate::new(),
        )
        .unwrap();
        let temperature_error = (0..3)
            .map(|i| {
                (candidate.thermal.temperature_k(CellId(i as u64)).unwrap()
                    - 300.0
                    - reference[i + 1])
                    .abs()
            })
            .fold(0.0_f64, f64::max);
        let travel_error = (candidate.slip_increment_m[0] - reference[4]).abs();
        errors.push([travel_error, temperature_error]);
        eprintln!(
            "Hertz thermal march: scale={scale}, phases={}, rejected={}, thermal_trials={}, travel_error={travel_error}, temperature_error={temperature_error}",
            candidate.phases.len(),
            candidate.rejected_attempts,
            candidate.thermal_trials
        );
        assert!(candidate.rejected_attempts > 0);
        assert_eq!(
            candidate.phases.len(),
            2 * (candidate.attempts - candidate.rejected_attempts)
        );
        assert!(candidate.maximum_accepted_normalized_defect <= 1.0);
        assert_eq!(candidate.thermal.time_s(), 0.3);
        close(candidate.dissipated_j, reference[5]);
        close(candidate.thermal_energy_change_j, reference[5]);
        close(candidate.kinetic_energy_change_j, -reference[5]);
        for residual in candidate.energy_residual_j {
            assert!(residual.abs() < 5e-8);
        }
        let mut clock = 0.0;
        let mut stop = None;
        let mut heat = 0.0;
        let mut trials = 0usize;
        let mut temperatures = Vec::new();
        for phase in &candidate.phases {
            assert_eq!(phase.start_s, clock);
            clock = phase.end_s;
            trials += phase.trials;
            heat += phase.heat.dissipated_j;
            temperatures.push(phase.source);
            assert_eq!(phase.friction_source.owner_identity(), owner.identity());
            assert_eq!(
                phase
                    .friction_source
                    .friction_source()
                    .state_coordinate("T"),
                Some(phase.source)
            );
            close(
                phase.source,
                phase.heat.interfaces[0].interface_temperature_k,
            );
            close(
                phase.friction_source.law().parameters().kinetic_mu,
                MU[0] - SLOPE[0] * (phase.source - 300.0),
            );
            let HertzHeatMarchMotion::InvariantDirection(m) = &phase.mechanical else {
                panic!("wrong native owner");
            };
            if stop.is_some() {
                assert!(phase.heat.dissipated_j < 1e-18);
            }
            if stop.is_none()
                && let Some(t) = m.relaxation.kinetic_phase_end_s
            {
                stop = Some(phase.start_s + t);
            }
        }
        assert!(stop.unwrap() > 0.05 && stop.unwrap() < 0.2);
        assert!(temperatures.windows(2).any(|w| (w[1] - w[0]).abs() > 1e-3));
        assert!(
            candidate.thermal_trials > trials,
            "coarse/refused trials are charged too"
        );
        close(heat, candidate.dissipated_j);
        assert_eq!(clock, 0.3);
    }
    for i in 0..2 {
        assert!(errors[1][i] < 0.5 * errors[0][i], "{errors:?}");
    }
    assert!(errors[1][0] < 1e-5 && errors[1][1] < 1e-4, "{errors:?}");
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(pair.bodies, Pair::new(0.5, 0.0).bodies);
}

#[test]
fn resting_cells_refine_their_own_cooling_even_with_zero_mechanical_error() {
    let pair = Pair::new(0.0, 0.0);
    let owner = law(0);
    let thermal = thermal([310.0, 290.0, 300.0], So3::identity());
    let route = routes(&thermal)[0].clone();
    let total = 0.03;
    let mut reference = [0.0, 10.0, -10.0, 0.0, 0.0, 0.0];
    for _ in 0..20_000 {
        reference = rk4(reference, total / 20_000.0, false);
    }
    let candidate = run_march(
        &pair,
        &owner,
        &thermal,
        &route,
        total,
        HertzHeatMarchModel::InvariantDirection,
        march_controls(1.0),
        source,
        &CancelGate::new(),
    )
    .unwrap();
    assert!(candidate.rejected_attempts > 0);
    assert_eq!(candidate.bodies_after, pair.bodies);
    assert_eq!(candidate.dissipated_j, 0.0);
    assert_eq!(candidate.slip_increment_m, [0.0; 2]);
    close(candidate.thermal_energy_change_j, 0.0);
    let mut maximum_error = 0.0_f64;
    for i in 0..3 {
        maximum_error = maximum_error.max(
            (candidate.thermal.temperature_k(CellId(i as u64)).unwrap() - 300.0 - reference[i + 1])
                .abs(),
        );
    }
    assert!(maximum_error < 0.02, "{maximum_error}");
    assert!(candidate.body_heat_j[0] < 0.0 && candidate.body_heat_j[1] > 0.0);
    assert!(candidate.phases.iter().all(|p| p.heat.dissipated_j == 0.0));
    assert_eq!(candidate.thermal.time_s(), total);
}

// Independent continuous two-rate/three-cell equations. Only the original
// pressure integral supplies instantaneous traction; no native transient,
// thermal Newton, or adaptive march is used by this small-step RK4 reference.
fn curved_reference(total: f64) -> [f64; 8] {
    let unit = HertzSlidingPatch::new(HertzSlidingPatchParameters {
        normal_force_n: 10.0,
        contact_radius_m: radius(),
        friction_coefficient: 1.0,
    })
    .unwrap();
    let derivative = |y: [f64; 8]| {
        let wrench = unit.evaluate([y[0], 0.0], y[1]).unwrap();
        let g0 = K[0] * AREA / A;
        let g1 = K[1] * AREA / A;
        let g02 = AREA / (A / K[0] + A / K[2]);
        let p = wrench.dissipated_power_w;
        let theta = (g0 * y[2] + g1 * y[3] + MU[0] * p) / (g0 + g1 + SLOPE[0] * p);
        let mu = MU[0] - SLOPE[0] * theta;
        [
            mobility(false) * mu * wrench.force_n[0],
            mobility(true) * mu * wrench.spin_moment_nm,
            (g0 * (theta - y[2]) + g02 * (y[4] - y[2])) / CAP[0],
            g1 * (theta - y[3]) / CAP[1],
            g02 * (y[2] - y[4]) / CAP[2],
            y[0],
            y[1],
            mu * p,
        ]
    };
    let mut y = [0.1, 3.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0];
    let dt = total / 20_000.0;
    for _ in 0..20_000 {
        let a = derivative(y);
        let b = derivative(std::array::from_fn(|i| y[i] + 0.5 * dt * a[i]));
        let c = derivative(std::array::from_fn(|i| y[i] + 0.5 * dt * b[i]));
        let d = derivative(std::array::from_fn(|i| y[i] + dt * c[i]));
        y = std::array::from_fn(|i| y[i] + dt * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0);
    }
    y
}

#[test]
fn curved_native_momenta_and_heat_follow_continuously_changing_temperature() {
    let pair = Pair::new(0.1, 3.0);
    let owner = law(0);
    let thermal = thermal([300.0; 3], So3::identity());
    let route = routes(&thermal)[0].clone();
    let total = 0.08;
    let reference = curved_reference(total);
    let candidate = run_march(
        &pair,
        &owner,
        &thermal,
        &route,
        total,
        HertzHeatMarchModel::Curved(path_controls()),
        march_controls(0.25),
        source,
        &CancelGate::new(),
    )
    .unwrap();
    let HertzHeatMarchMotion::Curved(m) = &candidate.phases.last().unwrap().mechanical else {
        panic!("wrong native owner");
    };
    assert!(!m.transient.exact_ray);
    assert!((m.slip_velocity_m_s[1][0] - reference[0]).abs() < 1e-5);
    assert!((m.spin_rad_s[1] - reference[1]).abs() < 1e-4);
    assert!((candidate.slip_increment_m[0] - reference[5]).abs() < 2e-6);
    assert!((candidate.spin_increment_rad - reference[6]).abs() < 1e-5);
    for i in 0..3 {
        let error =
            (candidate.thermal.temperature_k(CellId(i as u64)).unwrap() - 300.0 - reference[i + 2])
                .abs();
        assert!(error < 2e-4, "cell {i}: {error}");
    }
    let initial = 0.5 * (0.1 * 0.1 / mobility(false) + 3.0 * 3.0 / mobility(true));
    let remaining = 0.5
        * (m.slip_velocity_m_s[1][0].powi(2) / mobility(false)
            + m.spin_rad_s[1].powi(2) / mobility(true));
    close(candidate.dissipated_j, initial - remaining);
    close(candidate.thermal_energy_change_j, candidate.dissipated_j);
    assert!((candidate.dissipated_j - reference[7]).abs() < 1e-5);
    assert!(candidate.rejected_attempts > 0);
    eprintln!(
        "curved Hertz thermal march: phases={}, rejected={}, thermal_trials={}",
        candidate.phases.len(),
        candidate.rejected_attempts,
        candidate.thermal_trials
    );
}

#[test]
fn late_source_cancellation_and_work_limits_leave_the_whole_span_retryable() {
    let pair = Pair::new(0.5, 0.0);
    let owner = law(0);
    let thermal = thermal([300.0; 3], So3::identity());
    let before = thermal.cells().clone();
    let route = routes(&thermal)[0].clone();
    let mut limit = march_controls(1e6);
    limit.initial_step_s = 0.02;
    limit.maximum_step_s = 0.02;
    let model = HertzHeatMarchModel::InvariantDirection;
    let first = run_march(
        &pair,
        &owner,
        &thermal,
        &route,
        0.02,
        model,
        limit,
        source,
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(first.attempts, 1);
    assert_eq!(first.phases.len(), 2);
    let reference = run_march(
        &pair,
        &owner,
        &thermal,
        &route,
        0.06,
        model,
        limit,
        source,
        &CancelGate::new(),
    )
    .unwrap();
    let calls = Cell::new(0usize);
    let failed = run_march(
        &pair,
        &owner,
        &thermal,
        &route,
        0.06,
        model,
        limit,
        |state, t| {
            calls.set(calls.get() + 1);
            if calls.get() > first.thermal_trials {
                Err("late source refusal")
            } else {
                source(state, t)
            }
        },
        &CancelGate::new(),
    );
    assert!(matches!(
        failed,
        Err(HertzHeatMarchError::InvariantDirection(
            RigidHertzHeatError::Source("late source refusal")
        ))
    ));
    assert_eq!(calls.get(), first.thermal_trials + 1);
    let cancel = CancelGate::new();
    calls.set(0);
    assert!(
        run_march(
            &pair,
            &owner,
            &thermal,
            &route,
            0.06,
            model,
            limit,
            |state, t| {
                calls.set(calls.get() + 1);
                if calls.get() > first.thermal_trials {
                    cancel.request();
                }
                source(state, t)
            },
            &cancel
        )
        .is_err()
    );
    for budget in 0..3 {
        let mut bounded = limit;
        match budget {
            0 => bounded.maximum_attempts = 1,
            1 => bounded.maximum_accepted_phases = 2,
            _ => bounded.maximum_thermal_trials = first.thermal_trials,
        }
        assert!(matches!(
            run_march(
                &pair,
                &owner,
                &thermal,
                &route,
                0.06,
                model,
                bounded,
                source,
                &CancelGate::new()
            ),
            Err(HertzHeatMarchError::Budget(_))
        ));
    }
    let mut minimum = march_controls(1e-4);
    minimum.minimum_step_s = minimum.initial_step_s;
    assert!(matches!(
        run_march(
            &pair,
            &owner,
            &thermal,
            &route,
            0.3,
            model,
            minimum,
            source,
            &CancelGate::new()
        ),
        Err(HertzHeatMarchError::Unresolved { .. })
    ));
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(thermal.cells(), &before);
    assert_eq!(pair.bodies, Pair::new(0.5, 0.0).bodies);
    let retry = run_march(
        &pair,
        &owner,
        &thermal,
        &route,
        0.06,
        model,
        limit,
        source,
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(retry.bodies_after, reference.bodies_after);
    assert_eq!(retry.thermal.cells(), reference.thermal.cells());
    assert_eq!(retry.thermal_trials, reference.thermal_trials);
    assert_eq!(retry.energy_residual_j, reference.energy_residual_j);
}
