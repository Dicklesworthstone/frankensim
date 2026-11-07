//! G1/G3/G4: one recorded interpolation, original clocks, complete retarded fields.

use fs_aeroac::AeroacError;
use fs_aeroac::curle3d::array::{IdealAcousticPlane, IdealPlaneBoundary};
use fs_aeroac::curle3d::quadrupole::recorded::QuadrupoleSampleClock;
use fs_aeroac::curle3d::quadrupole::transient::history::{
    PreparedQuadrupoleHistory, QuadrupoleHistoryBudget, QuadrupoleHistoryEnds,
};
use fs_aeroac::curle3d::quadrupole::transient::{
    PreparedTransientQuadrupoleObserver, RetardedQuadrupoleField, TransientQuadrupoleSource,
};
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::dd::Dd;

const DURATION: f64 = 0.04;
const POSITIONS: [[f64; 3]; 2] = [[0.1, 0.2, 0.15], [-0.2, 0.1, 0.3]];

fn medium() -> FluidAcousticState {
    FluidAcousticState::try_new(1.23, 320.0).unwrap()
}

fn budget() -> QuadrupoleHistoryBudget {
    QuadrupoleHistoryBudget {
        maximum_sources: 8,
        maximum_samples_per_source: 4096,
        maximum_workspace_bytes: 8 << 20,
    }
}

fn close(actual: f64, expected: f64, scale: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * scale.abs().max(1e-25),
        "actual={actual:e}, expected={expected:e}, scale={scale:e}"
    );
}

// Six independent analytic cubics, with a nonzero physical initial integral.
// No spline, finite difference or propagation formula supplies this oracle.
fn cubic(source: usize, elapsed: f64) -> TransientQuadrupoleSource {
    let mut result = TransientQuadrupoleSource::default();
    let z = elapsed / DURATION;
    for j in 0..6 {
        let sign = if j % 2 == 0 { 1.0 } else { -1.0 };
        let a = sign * (j + 1) as f64 * 0.11 + source as f64 * 0.07;
        let b = sign * (j + source + 2) as f64 * 0.31;
        let c = -sign * (j + 3) as f64 * 0.19;
        let d = sign * (source + 1) as f64 * 0.13;
        result.stress_moment_n_m[j] = a + b * z + c * z.powi(2) + d * z.powi(3);
        result.stress_rate_n_m_s[j] = (b + 2.0 * c * z + 3.0 * d * z.powi(2)) / DURATION;
        result.stress_acceleration_n_m_s2[j] = (2.0 * c + 6.0 * d * z) / DURATION.powi(2);
        result.stress_integral_n_m_s[j] = sign * (j + source + 1) as f64 * 0.003
            + DURATION * (a * z + b * z.powi(2) / 2.0 + c * z.powi(3) / 3.0 + d * z.powi(4) / 4.0);
    }
    result
}

#[derive(Clone)]
struct Record {
    positions: Vec<[f64; 3]>,
    clock: QuadrupoleSampleClock,
    times: Vec<f64>,
    moments: Vec<[f64; 6]>,
    ends: Vec<QuadrupoleHistoryEnds>,
}

impl Record {
    fn cubic(epoch: f64, samples: usize) -> Self {
        let h = DURATION / (samples - 1) as f64;
        let clock = QuadrupoleSampleClock::new(epoch, h, 0.0).unwrap();
        Self {
            positions: POSITIONS.to_vec(),
            clock,
            times: (0..samples).map(|i| clock.time_at(i).unwrap()).collect(),
            moments: (0..samples)
                .flat_map(|i| {
                    (0..POSITIONS.len()).map(move |s| cubic(s, i as f64 * h).stress_moment_n_m)
                })
                .collect(),
            ends: (0..POSITIONS.len())
                .map(|s| QuadrupoleHistoryEnds {
                    initial_rate_n_m_s: cubic(s, 0.0).stress_rate_n_m_s,
                    final_rate_n_m_s: cubic(s, (samples - 1) as f64 * h).stress_rate_n_m_s,
                    initial_integral_n_m_s: cubic(s, 0.0).stress_integral_n_m_s,
                })
                .collect(),
        }
    }

    fn prepare(&self) -> Result<PreparedQuadrupoleHistory, AeroacError> {
        PreparedQuadrupoleHistory::new(
            &self.positions,
            self.clock,
            &self.times,
            &self.moments,
            &self.ends,
            budget(),
            || false,
        )
    }
}

fn check_field(actual: RetardedQuadrupoleField, expected: RetardedQuadrupoleField) {
    close(
        actual.pressure_pa,
        expected.pressure_pa,
        expected.pressure_pa.abs().max(1e-10),
        2e-10,
    );
    let velocity_scale = expected
        .particle_velocity_m_s
        .iter()
        .map(|v| v.abs())
        .fold(0.0, f64::max);
    for axis in 0..3 {
        close(
            actual.particle_velocity_m_s[axis],
            expected.particle_velocity_m_s[axis],
            velocity_scale,
            2e-10,
        );
        assert_eq!(
            actual.instantaneous_intensity_w_m2[axis],
            actual.pressure_pa * actual.particle_velocity_m_s[axis]
        );
        close(
            actual.instantaneous_intensity_w_m2[axis],
            expected.instantaneous_intensity_w_m2[axis],
            expected.pressure_pa.abs() * velocity_scale,
            4e-10,
        );
    }
}

#[test]
fn physical_cubic_record_recovers_all_derivatives_and_the_original_integral() {
    for samples in [2, 9, 23] {
        let record = Record::cubic(0.0, samples);
        let history = record.prepare().unwrap();
        assert_eq!(history.clock(), record.clock);
        assert_eq!(history.sample_count(), samples);
        assert_eq!(history.source_positions_m(), POSITIONS);
        let span = Dd::from_f64(record.clock.timestep_s()) * Dd::from_f64((samples - 1) as f64);
        for source in 0..POSITIONS.len() {
            for j in 0..=72 {
                let time = span * Dd::from_f64(j as f64 / 72.0);
                let elapsed = time.to_f64();
                let actual = history.source_at_clock_dd(source, time).unwrap();
                let expected = cubic(source, elapsed);
                for component in 0..6 {
                    close(
                        actual.stress_moment_n_m[component],
                        expected.stress_moment_n_m[component],
                        4.0,
                        3e-13,
                    );
                    close(
                        actual.stress_rate_n_m_s[component],
                        expected.stress_rate_n_m_s[component],
                        128.0,
                        2e-12,
                    );
                    close(
                        actual.stress_acceleration_n_m_s2[component],
                        expected.stress_acceleration_n_m_s2[component],
                        8192.0,
                        3e-12,
                    );
                    close(
                        actual.stress_integral_n_m_s[component],
                        expected.stress_integral_n_m_s[component],
                        0.1,
                        3e-13,
                    );
                }
            }
        }
    }
}

#[test]
fn recorded_sources_drive_coherent_pressure_velocity_and_every_echo() {
    let record = Record::cubic(0.0, 9);
    let history = record.prepare().unwrap();
    let receiver = [0.8, -0.4, 0.7];
    for boundary in [
        IdealPlaneBoundary::Rigid,
        IdealPlaneBoundary::PressureRelease,
    ] {
        let plane = IdealAcousticPlane::new([0.0; 3], [0.0, 0.0, 1.0], boundary).unwrap();
        let raw = PreparedTransientQuadrupoleObserver::new(medium(), receiver, &POSITIONS, 2)
            .unwrap()
            .with_orthogonal_planes(&[plane], 4)
            .unwrap();
        let observer = history.bind_observer(raw.clone(), || false).unwrap();
        assert_eq!(observer.propagation().path_count(), 4);
        assert_eq!(observer.history().source_positions_m(), POSITIONS);
        let [first, last] = observer.receiver_time_bounds_s();
        for fraction in [0.1, 0.4, 0.75] {
            let time = first + (last - first) * Dd::from_f64(fraction);
            let actual = observer.evaluate_at_clock_dd(time, || false).unwrap();
            let expected = raw
                .evaluate_at_clock_dd(time, |s, t| Ok(cubic(s, t.to_f64())))
                .unwrap();
            check_field(actual, expected);
        }
        for time in [first - Dd::from_f64(1e-8), last + Dd::from_f64(1e-8)] {
            assert!(observer.evaluate_at_clock_dd(time, || false).is_err());
        }
        let middle = (first + last) * Dd::from_f64(0.5);
        assert!(matches!(
            observer.evaluate_at_clock_dd(middle, || true),
            Err(AeroacError::Cancelled { .. })
        ));
        let after = observer.evaluate_at_clock_dd(middle, || false).unwrap();
        let expected = raw
            .evaluate_at_clock_dd(middle, |s, t| Ok(cubic(s, t.to_f64())))
            .unwrap();
        check_field(after, expected);
    }
}

#[test]
fn original_dd_clock_keeps_large_epoch_fields_and_query_order_consistent() {
    let base = Record::cubic(0.0, 9).prepare().unwrap();
    let epoch = 2.0_f64.powi(30);
    let shifted = Record::cubic(epoch, 9).prepare().unwrap();
    let a = base
        .prepare_observer(medium(), [0.8, -0.4, 0.7], || false)
        .unwrap();
    let b = shifted
        .prepare_observer(medium(), [0.8, -0.4, 0.7], || false)
        .unwrap();
    // Deliberately nonmonotone queries: there is no moving integral cursor.
    for t in [0.03, 0.0123456789, 0.0223456789, 0.03] {
        let original_time = Dd::from_f64(t);
        let absolute_time = Dd::from_f64(epoch) + original_time;
        let actual = b.evaluate_at_clock_dd(absolute_time, || false).unwrap();
        let expected = a.evaluate_at_clock_dd(original_time, || false).unwrap();
        check_field(actual, expected);
    }
    for time in shifted.source_time_bounds_s() {
        for source in 0..POSITIONS.len() {
            shifted.source_at_clock_dd(source, time).unwrap();
        }
    }
    // Both endpoints of a single-path covered receiver interval remain usable.
    let history = PreparedQuadrupoleHistory::new(
        &[[0.0; 3]],
        QuadrupoleSampleClock::new(epoch, 0.125, 0.0).unwrap(),
        &[epoch, epoch + 0.125],
        &[[0.0; 6], [1.0; 6]],
        &[QuadrupoleHistoryEnds {
            initial_rate_n_m_s: [0.0; 6],
            final_rate_n_m_s: [0.0; 6],
            initial_integral_n_m_s: [0.0; 6],
        }],
        budget(),
        || false,
    )
    .unwrap();
    let observer = history
        .prepare_observer(
            FluidAcousticState::try_new(1.0, 8.0).unwrap(),
            [1.0, 0.0, 0.0],
            || false,
        )
        .unwrap();
    for time in observer.receiver_time_bounds_s() {
        observer.evaluate_at_clock_dd(time, || false).unwrap();
    }
}

#[test]
fn radiation_driving_second_derivative_converges_under_original_clock_refinement() {
    let omega = core::f64::consts::TAU * 40.0;
    let duration = 0.1;
    let mut errors = Vec::new();
    for samples in [17, 33, 65] {
        let h = duration / (samples - 1) as f64;
        let clock = QuadrupoleSampleClock::new(0.0, h, 0.0).unwrap();
        let times: Vec<_> = (0..samples).map(|i| clock.time_at(i).unwrap()).collect();
        let moments: Vec<_> = times.iter().map(|&t| [(omega * t).sin(); 6]).collect();
        let end_time = clock.time_at(samples - 1).unwrap();
        let history = PreparedQuadrupoleHistory::new(
            &[[0.0; 3]],
            clock,
            &times,
            &moments,
            &[QuadrupoleHistoryEnds {
                initial_rate_n_m_s: [omega; 6],
                final_rate_n_m_s: [omega * (omega * end_time).cos(); 6],
                initial_integral_n_m_s: [0.0; 6],
            }],
            budget(),
            || false,
        )
        .unwrap();
        let mut error = 0.0_f64;
        for j in 1..1000 {
            let t = duration * j as f64 / 1000.0;
            let actual = history.source_at_clock_dd(0, Dd::from_f64(t)).unwrap();
            let expected = -omega * omega * (omega * t).sin();
            error = error
                .max((actual.stress_acceleration_n_m_s2[0] - expected).abs() / (omega * omega));
        }
        errors.push(error);
    }
    assert!(
        errors[0] > 1e-2,
        "the refinement fixture must resolve a live source error"
    );
    for pair in errors.windows(2) {
        assert!(
            pair[1] < 0.34 * pair[0],
            "second-derivative refinement: {errors:?}"
        );
    }
}

#[test]
fn reported_non_dyadic_arrivals_are_usable_at_large_epochs() {
    let epoch = 1e12;
    let clock = QuadrupoleSampleClock::new(epoch, 0.003, 0.0).unwrap();
    let times: Vec<_> = (0..3).map(|i| clock.time_at(i).unwrap()).collect();
    let history = PreparedQuadrupoleHistory::new(
        &[[0.0; 3]],
        clock,
        &times,
        &[[0.0; 6], [0.4; 6], [0.9; 6]],
        &[QuadrupoleHistoryEnds {
            initial_rate_n_m_s: [0.0; 6],
            final_rate_n_m_s: [1.0; 6],
            initial_integral_n_m_s: [0.04; 6],
        }],
        budget(),
        || false,
    )
    .unwrap();
    let observer = history
        .prepare_observer(
            FluidAcousticState::try_new(1.2, 343.0).unwrap(),
            [0.007, 0.0, 0.0],
            || false,
        )
        .unwrap();
    let arrivals = observer.receiver_time_bounds_s();
    for (arrival, emission) in arrivals.into_iter().zip(history.source_time_bounds_s()) {
        let expected = observer
            .propagation()
            .evaluate_at_retarded(&[history.source_at_clock_dd(0, emission).unwrap()])
            .unwrap();
        check_field(
            observer.evaluate_at_clock_dd(arrival, || false).unwrap(),
            expected,
        );
    }
    // These offsets are below an absolute f64 ULP but resolve on the DD clock.
    assert!(
        observer
            .evaluate_at_clock_dd(arrivals[0] - Dd::from_f64(1e-10), || false)
            .is_err()
    );
    assert!(
        observer
            .evaluate_at_clock_dd(arrivals[1] + Dd::from_f64(1e-10), || false)
            .is_err()
    );
    // A tiny delay may resolve at both dyadic endpoints while disappearing
    // at an interior DD timestamp with a nonzero low part. Preserve the
    // original propagator's per-query strict-retardation admission.
    let tiny_clock = QuadrupoleSampleClock::new(epoch, 0.125, 0.0).unwrap();
    let tiny = PreparedQuadrupoleHistory::new(
        &[[0.0; 3]],
        tiny_clock,
        &[epoch, epoch + 0.125],
        &[[0.0; 6]; 2],
        &[QuadrupoleHistoryEnds {
            initial_rate_n_m_s: [0.0; 6],
            final_rate_n_m_s: [0.0; 6],
            initial_integral_n_m_s: [0.0; 6],
        }],
        budget(),
        || false,
    )
    .unwrap();
    let tiny_observer = tiny
        .prepare_observer(
            FluidAcousticState::try_new(1.2, 343.0).unwrap(),
            [3.43e-22, 0.0, 0.0],
            || false,
        )
        .unwrap();
    assert!(
        tiny_observer
            .evaluate_at_clock_dd(Dd::from_f64(epoch) + Dd::from_f64(0.0635), || false,)
            .is_err()
    );
}

#[test]
fn sampled_history_preserves_all_three_momentum_equations() {
    let history = Record::cubic(0.0, 9).prepare().unwrap();
    let position = [0.8, -0.4, 0.7];
    let time = 0.02;
    let dt = 1e-7;
    let dx = 2e-5;
    let owner = history
        .prepare_observer(medium(), position, || false)
        .unwrap();
    let before = owner
        .evaluate_at_clock_dd(Dd::from_f64(time - dt), || false)
        .unwrap();
    let after = owner
        .evaluate_at_clock_dd(Dd::from_f64(time + dt), || false)
        .unwrap();
    for axis in 0..3 {
        let mut left = position;
        let mut right = position;
        left[axis] -= dx;
        right[axis] += dx;
        let p0 = history
            .prepare_observer(medium(), left, || false)
            .unwrap()
            .evaluate_at_clock_dd(Dd::from_f64(time), || false)
            .unwrap()
            .pressure_pa;
        let p1 = history
            .prepare_observer(medium(), right, || false)
            .unwrap()
            .evaluate_at_clock_dd(Dd::from_f64(time), || false)
            .unwrap()
            .pressure_pa;
        let gradient = (p1 - p0) / (2.0 * dx);
        let acceleration = medium().density_kg_m3()
            * (after.particle_velocity_m_s[axis] - before.particle_velocity_m_s[axis])
            / (2.0 * dt);
        close(
            acceleration + gradient,
            0.0,
            acceleration.abs() + gradient.abs(),
            3e-7,
        );
    }
}

#[test]
fn source_clock_geometry_coverage_budgets_and_cancellation_refuse_without_fabrication() {
    let good = Record::cubic(0.0, 9);
    for case in 0..6 {
        let mut bad = good.clone();
        match case {
            0 => bad.times[4] += 0.25 * bad.clock.timestep_s(),
            1 => bad.moments[7][4] = f64::NAN,
            2 => bad.ends[1].initial_rate_n_m_s[0] = f64::INFINITY,
            3 => bad.ends[1].initial_integral_n_m_s[3] = f64::NAN,
            4 => bad.positions[0][2] = f64::NAN,
            _ => {
                bad.moments.pop();
            }
        }
        assert!(bad.prepare().is_err());
    }
    let history = good.prepare().unwrap();
    let mut small = budget();
    small.maximum_workspace_bytes = history.workspace_bytes() - 1;
    assert!(
        PreparedQuadrupoleHistory::new(
            &good.positions,
            good.clock,
            &good.times,
            &good.moments,
            &good.ends,
            small,
            || false
        )
        .is_err()
    );
    let mut polls = 0;
    assert!(matches!(
        PreparedQuadrupoleHistory::new(
            &good.positions,
            good.clock,
            &good.times,
            &good.moments,
            &good.ends,
            budget(),
            || {
                polls += 1;
                polls == 7
            }
        ),
        Err(AeroacError::Cancelled { .. })
    ));
    assert!(history.source_at_clock_dd(0, Dd::from_f64(-1e-9)).is_err());
    assert!(
        history
            .source_at_clock_dd(0, Dd::from_f64(DURATION + 1e-9))
            .is_err()
    );
    assert!(history.source_at_clock_dd(2, Dd::from_f64(0.02)).is_err());
    assert!(
        history
            .source_at_clock_dd(0, Dd::from_f64(f64::NAN))
            .is_err()
    );
    let wrong = PreparedTransientQuadrupoleObserver::new(
        medium(),
        [0.8, -0.4, 0.7],
        &[POSITIONS[1], POSITIONS[0]],
        2,
    )
    .unwrap();
    assert!(history.bind_observer(wrong, || false).is_err());
    // An echo with no common recorded receiver interval must not become silence.
    let far_plane = IdealAcousticPlane::new(
        [0.0, 0.0, -100.0],
        [0.0, 0.0, 1.0],
        IdealPlaneBoundary::Rigid,
    )
    .unwrap();
    let echoes =
        PreparedTransientQuadrupoleObserver::new(medium(), [0.8, -0.4, 0.7], &POSITIONS, 2)
            .unwrap()
            .with_orthogonal_planes(&[far_plane], 4)
            .unwrap();
    assert!(history.bind_observer(echoes, || false).is_err());
}
