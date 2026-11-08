use super::super::*;
use super::*;
use scalar::LagrangianScalarConfig;
use scalar::components::{LagrangianComponentsConfig, LagrangianDynamicComponents};

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-12 * expected.abs().max(1e-280),
        "{actual:e} != {expected:e}"
    );
}

fn adaptive() -> AdaptiveLagrangianMemory {
    AdaptiveLagrangianMemory {
        coefficient: 0.02,
        minimum_time: 0.005,
    }
}

fn config() -> LagrangianSgsConfig {
    LagrangianSgsConfig {
        memory_time: 2.0,
        maximum_departure_cells: 1.0,
        maximum_coefficient_squared: 100.0,
        maximum_eddy_to_molecular: 1000.0,
        usage: SgsUse::Exploration,
    }
}

fn update(
    previous: Option<f64>,
    current: Option<f64>,
    elapsed: f64,
    maximum: f64,
) -> HistoryUpdate {
    let flow = PeriodicDynamicSgs::new([3; 3], 1.0, usize::MAX).unwrap();
    HistoryUpdate::new(&flow, true, previous, current, elapsed, maximum, 1.0).unwrap()
}

#[test]
fn invariant_time_preserves_velocity_length_and_density_similarity() {
    let controls = AdaptiveLagrangianMemory {
        coefficient: 1.5,
        minimum_time: 1e-200,
    };
    let pair = [0.25, 1.0];
    // Physical contractions are [16,64] and rho=2. No production logarithmic
    // arithmetic is used in this independent moderate-scale dimensional oracle.
    let expected = 1.5 * 0.3 * (16.0_f64 * 64.0 / 2.0_f64.powi(4)).powf(-0.125);
    for velocity in [1e-40_f64, 0.01, 1.0, 100.0, 1e40] {
        for length in [1e-20_f64, 1.0, 1e20] {
            for density in [1e-100_f64, 1.0, 1e100] {
                let log = 64.0_f64.ln() + 4.0 * velocity.ln() + 2.0 * density.ln();
                let base = update(Some(log), Some(-700.0), 0.1, 1e200);
                let rule = AdaptiveUpdate::new(
                    controls,
                    1e200,
                    0.3 * length,
                    Some(log),
                    Some(-700.0),
                    &base,
                );
                close(
                    rule.memory_time(pair, 2.0 * density),
                    expected * length / velocity,
                );
            }
        }
    }
}

#[test]
fn advected_history_controls_local_time_and_exact_frozen_source_relaxation() {
    let mut history = vec![[0.0, 1.0]; 27];
    history[12] = [0.16, 3.2];
    history[13] = [0.04, 0.8];
    let mut current = vec![[0.0; 2]; 27];
    current[13] = [0.2, 0.7];
    let old_log = 9.0_f64.ln();
    let new_log = 4.0_f64.ln();
    let base = update(Some(old_log), Some(new_log), 1.0, 10.0);
    let rule = AdaptiveUpdate::new(
        AdaptiveLagrangianMemory {
            coefficient: 1.5,
            minimum_time: 0.01,
        },
        10.0,
        1.0,
        Some(old_log),
        Some(new_log),
        &base,
    );
    // Half-cell backtracking samples BOTH upstream centers. Using the current
    // center, or the current amplitude's logarithmic scale, gives a wrong T.
    let advected = [0.1, 2.0];
    let expected_time = 1.5 * (9.0_f64.powi(2) * advected[0] * advected[1]).powf(-0.125);
    let (pair, time) = rule
        .pair(&base, 13, [0.5, 0.0, 0.0], 1.0, &history, &current)
        .unwrap();
    close(time, expected_time);
    let keep = (-1.0 / expected_time).exp();
    for k in 0..2 {
        let expected = keep * 9.0 * advected[k] + (1.0 - keep) * 4.0 * current[13][k];
        close(pair[k] * rule.common_log.exp(), expected);
    }
}

#[test]
fn authored_bounds_do_not_erase_signed_history_or_require_representable_products() {
    let base = update(Some(1000.0), Some(0.0), 1000.0, 1.0);
    let rule = AdaptiveUpdate::new(
        AdaptiveLagrangianMemory {
            coefficient: 1.0,
            minimum_time: 1.0,
        },
        1.0,
        1.0,
        Some(1000.0),
        Some(0.0),
        &base,
    );
    let history = vec![[-0.5, 1.0]; 27];
    let current = vec![[0.25, 0.5]; 27];
    let (pair, time) = rule
        .pair(&base, 13, [0.0; 3], 1.0, &history, &current)
        .unwrap();
    // exp(1000)*exp(-1000)=1 even though neither separate exponential is
    // representable. Signed I_LM is retained through the bounded fallback.
    assert_eq!(time, 1.0);
    assert_eq!(rule.common_log, 0.0);
    assert_eq!(pair, [-0.25, 1.5]);
    assert_eq!(rule.memory_time([1.0, 1.0], 1.0), 1.0);

    let base = update(Some(0.0), Some(0.0), 0.1, 2.0);
    let rule = AdaptiveUpdate::new(adaptive(), 2.0, 1.0, Some(0.0), Some(0.0), &base);
    assert_eq!(rule.memory_time([0.0, 1.0], 1.0), 2.0);
    assert_eq!(rule.memory_time([-1.0, 1.0], 1.0), 2.0);
    assert_eq!(
        rule.memory_time([1e100, 1e100], 1.0),
        adaptive().minimum_time
    );
    assert_eq!(rule.memory_time([1e-100, 1e-100], 1.0), 2.0);
}

fn field(amplitude: f64) -> (Vec<f64>, Vec<[f64; 3]>, Vec<f64>) {
    let mut rho = Vec::new();
    let mut u = Vec::new();
    let mut temperature = Vec::new();
    for z in 0..5 {
        for y in 0..5 {
            for x in 0..5 {
                let [x, y, z] = [x, y, z].map(|i| core::f64::consts::TAU * (i as f64 + 0.5) / 5.0);
                rho.push(1.2 + 0.2 * x.sin());
                u.push([
                    amplitude * (0.3 * x.sin() + 0.2 * (y + z).sin()),
                    amplitude * (0.2 * y.cos() + 0.3 * (z + x).sin()),
                    amplitude * (0.2 * z.sin() + 0.3 * (x + y).cos()),
                ]);
                temperature.push(2.0 + 0.2 * (x + 2.0 * y).sin() + 0.1 * z.cos());
            }
        }
    }
    (rho, u, temperature)
}

fn prepare_owner<'a>(
    owner: &'a mut LagrangianDynamicSgs,
    config: LagrangianSgsConfig,
    time: f64,
    rho: &[f64],
    velocity: &[[f64; 3]],
    cancelled: impl FnMut() -> bool,
) -> Result<LagrangianSgsTrial<'a>, DynamicSgsError> {
    if let Some(normal) = owner.normal_axis() {
        let mut walls = [[0.03; 3], [-0.01; 3]];
        walls[0][normal] = 0.0;
        walls[1][normal] = 0.0;
        owner.prepare_with_walls(config, time, rho, velocity, walls, cancelled)
    } else {
        owner.prepare(config, time, rho, velocity, cancelled)
    }
}

#[test]
fn adaptive_periodic_and_channel_histories_publish_only_with_their_consumer() {
    for normal in [None, Some(0), Some(1), Some(2)] {
        let make = || {
            let owner = match normal {
                None => LagrangianDynamicSgs::new([5; 3], [1.0; 3], usize::MAX),
                Some(axis) => LagrangianDynamicSgs::new_channel([5; 3], [1.0; 3], axis, usize::MAX),
            }
            .unwrap();
            let bytes = owner.workspace_bytes();
            let owner = owner.with_adaptive_memory(adaptive()).unwrap();
            assert_eq!(owner.workspace_bytes(), bytes);
            owner
        };
        let mut tested = make();
        let mut clean = make();
        let mut buffers = [tested.history.as_ptr(), tested.candidate.as_ptr()];
        buffers.sort();
        let (rho, velocity, _) = field(1.0);
        for owner in [&mut tested, &mut clean] {
            let trial = prepare_owner(owner, config(), 0.0, &rho, &velocity, || false).unwrap();
            assert_eq!(trial.memory_time_range(), [config().memory_time; 2]);
            trial.commit();
        }
        let accepted = (
            tested.history.clone(),
            tested.history_log_scale,
            tested.accepted_time,
            tested.accepted_memory_time_range(),
        );
        let (rho, velocity, _) = field(0.7);
        let trial = prepare_owner(&mut tested, config(), 0.2, &rho, &velocity, || false).unwrap();
        assert!(trial.memory_time_range()[0] < config().memory_time);
        assert!(trial.memory_time_range()[0] >= adaptive().minimum_time);
        drop(trial); // a consuming flow operation refused
        let denied = LagrangianSgsConfig {
            maximum_departure_cells: 1e-6,
            ..config()
        };
        assert!(prepare_owner(&mut tested, denied, 0.2, &rho, &velocity, || false).is_err());
        let mut polls = 0;
        assert!(
            prepare_owner(&mut tested, config(), 0.2, &rho, &velocity, || {
                polls += 1;
                polls == 5
            })
            .is_err()
        );
        assert_eq!(
            (
                tested.history.clone(),
                tested.history_log_scale,
                tested.accepted_time,
                tested.accepted_memory_time_range()
            ),
            accepted
        );
        for owner in [&mut tested, &mut clean] {
            prepare_owner(owner, config(), 0.2, &rho, &velocity, || false)
                .unwrap()
                .commit();
        }
        assert_eq!(tested.history, clean.history);
        assert_eq!(tested.history_log_scale, clean.history_log_scale);
        assert_eq!(
            tested.accepted_memory_time_range(),
            clean.accepted_memory_time_range()
        );
        let mut final_buffers = [tested.history.as_ptr(), tested.candidate.as_ptr()];
        final_buffers.sort();
        assert_eq!(buffers, final_buffers);
        assert!(tested.with_adaptive_memory(adaptive()).is_err());
    }
}

#[test]
fn adaptive_momentum_does_not_replace_independent_heat_or_species_memory() {
    let heat = LagrangianScalarConfig {
        memory_time: 0.4,
        maximum_departure_cells: 1.0,
        maximum_coefficient: 100.0,
        maximum_diffusivity_to_molecular: 1000.0,
        usage: SgsUse::Exploration,
    };
    let species = LagrangianComponentsConfig {
        memory_time: 0.7,
        maximum_departure_cells: 1.0,
        maximum_coefficient: 100.0,
        usage: SgsUse::Exploration,
    };
    let mut fixed = LagrangianDynamicComponents::new([5; 3], [1.0; 3], 2, usize::MAX).unwrap();
    let mut varied = LagrangianDynamicComponents::new([5; 3], [1.0; 3], 2, usize::MAX)
        .unwrap()
        .with_adaptive_momentum_memory(adaptive())
        .unwrap();
    assert_eq!(fixed.workspace_bytes(), varied.workspace_bytes());
    for (step, amplitude) in [1.0, 0.7, 0.3].into_iter().enumerate() {
        let (rho, u, t) = field(amplitude);
        let component = |cell: usize, k: usize| {
            let first = 0.5 + 0.2 * (cell as f64 * 0.13).sin();
            if k == 0 { first } else { 1.0 - first }
        };
        let a = fixed
            .prepare(
                config(),
                heat,
                species,
                step as f64 * 0.2,
                &rho,
                &u,
                &t,
                component,
                || false,
            )
            .unwrap();
        let b = varied
            .prepare(
                config(),
                heat,
                species,
                step as f64 * 0.2,
                &rho,
                &u,
                &t,
                component,
                || false,
            )
            .unwrap();
        assert_eq!(a.scalar_report(), b.scalar_report());
        assert_eq!(a.components_report(), b.components_report());
        assert_eq!(a.scalar_diffusivities(), b.scalar_diffusivities());
        assert_eq!(a.component_diffusivities(), b.component_diffusivities());
        if step > 0 {
            assert!(b.momentum_memory_time_range()[0] < a.momentum_memory_time_range()[0]);
        }
        a.commit();
        b.commit();
    }
}

#[test]
fn adaptive_controls_refuse_invalid_or_inconsistent_times_before_preparation() {
    for coefficient in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            AdaptiveLagrangianMemory {
                coefficient,
                ..adaptive()
            }
            .validate()
            .is_err()
        );
    }
    for minimum_time in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            AdaptiveLagrangianMemory {
                minimum_time,
                ..adaptive()
            }
            .validate()
            .is_err()
        );
    }
    let (rho, u, _) = field(1.0);
    let mut owner = LagrangianDynamicSgs::new([5; 3], [1.0; 3], usize::MAX)
        .unwrap()
        .with_adaptive_memory(AdaptiveLagrangianMemory {
            minimum_time: 3.0,
            ..adaptive()
        })
        .unwrap();
    let before = owner.flow.contraction.clone();
    assert!(owner.prepare(config(), 0.0, &rho, &u, || false).is_err());
    assert_eq!(owner.flow.contraction, before);
    assert_eq!(owner.accepted_time(), None);
    assert_eq!(owner.accepted_memory_time_range(), None);
}
