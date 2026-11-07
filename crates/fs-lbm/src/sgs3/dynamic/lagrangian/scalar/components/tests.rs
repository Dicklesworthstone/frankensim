use super::*;

const SHAPE: [usize; 3] = [5; 3];
const SPACING: [f64; 3] = [0.8, 1.0, 1.2];
const CELLS: usize = 125;

fn controls() -> (
    LagrangianSgsConfig,
    LagrangianScalarConfig,
    LagrangianComponentsConfig,
) {
    (
        LagrangianSgsConfig {
            memory_time: 1.5,
            maximum_departure_cells: 1.0,
            maximum_coefficient_squared: 100.0,
            maximum_eddy_to_molecular: 1000.0,
            usage: SgsUse::Exploration,
        },
        LagrangianScalarConfig {
            memory_time: 0.4,
            maximum_departure_cells: 1.0,
            maximum_coefficient: 100.0,
            maximum_diffusivity_to_molecular: 1000.0,
            usage: SgsUse::Exploration,
        },
        LagrangianComponentsConfig {
            memory_time: 0.7,
            maximum_departure_cells: 1.0,
            maximum_coefficient: 100.0,
            usage: SgsUse::Exploration,
        },
    )
}

fn field(speed: f64, amplitude: f64) -> (Vec<f64>, Vec<[f64; 3]>, Vec<f64>) {
    let mut rho = Vec::new();
    let mut velocity = Vec::new();
    let mut heat = Vec::new();
    for cell in 0..CELLS {
        let at = coordinates(SHAPE, cell);
        let [x, y, z] = core::array::from_fn(|axis| {
            core::f64::consts::TAU * (at[axis] as f64 + 0.5) / SHAPE[axis] as f64
        });
        rho.push(1.2 + 0.2 * x.sin() + 0.1 * y.cos());
        velocity.push([
            speed * (0.3 * x.sin() + 0.2 * (y + z).sin() + 0.1 * (2.0 * x - y).cos()),
            speed * (0.2 * y.cos() + 0.3 * (z + x).sin() + 0.1 * (x - 2.0 * z).sin()),
            speed * (0.2 * z.sin() + 0.1 * (x + y).cos() + 0.3 * (2.0 * y - z).sin()),
        ]);
        heat.push(
            amplitude
                * (0.7 * (x + y).sin() + 0.9 * (2.0 * y - z).cos() + 0.4 * (x - 2.0 * z).sin()),
        );
    }
    (rho, velocity, heat)
}

fn fractions(heat: &[f64]) -> [Vec<f64>; 4] {
    let a: Vec<_> = heat.iter().map(|v| 0.25 + 0.03 * v).collect();
    let b: Vec<_> = (0..CELLS)
        .map(|i| 0.35 + 0.02 * heat[(i + 17) % CELLS])
        .collect();
    let c = a.iter().zip(&b).map(|(a, b)| 1.0 - a - b).collect();
    [a, b, c, vec![0.0; CELLS]]
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

fn owner(count: usize) -> LagrangianDynamicComponents {
    LagrangianDynamicComponents::new(SHAPE, SPACING, count, usize::MAX).unwrap()
}

#[test]
fn g1_one_component_reproduces_scalar_history_and_preserves_heat_and_momentum_bits() {
    let (momentum, heat_config, mut component_config) = controls();
    component_config.memory_time = heat_config.memory_time;
    let mut all = owner(1);
    let mut original = LagrangianDynamicScalar::new(SHAPE, SPACING, usize::MAX).unwrap();
    for (time, speed, amplitude) in [(0.0, 1.0, 2.0), (0.1, 0.4, 0.05), (0.25, 0.7, 3.0)] {
        let (rho, velocity, heat) = field(speed, amplitude);
        let trial = all
            .prepare(
                momentum,
                heat_config,
                component_config,
                time,
                &rho,
                &velocity,
                &heat,
                |_, cell| heat[cell],
                || false,
            )
            .unwrap();
        let baseline = original
            .prepare(momentum, heat_config, time, &rho, &velocity, &heat, || {
                false
            })
            .unwrap();
        assert_eq!(trial.momentum_report(), baseline.momentum_report());
        assert_eq!(trial.scalar_report(), baseline.scalar_report());
        assert_eq!(trial.components_report(), baseline.scalar_report());
        assert_eq!(trial.viscosities(), baseline.viscosities());
        assert_eq!(trial.scalar_diffusivities(), baseline.diffusivities());
        assert_eq!(trial.component_diffusivities(), baseline.diffusivities());
        assert_eq!(trial.component_coefficients().len(), CELLS);
        assert_eq!(trial.component_coefficient(CELLS), None);
        for (cell, coefficient) in trial.component_coefficients().enumerate() {
            assert_eq!(trial.component_coefficient(cell), Some(coefficient));
            assert_eq!(Some(coefficient), baseline.scalar_coefficient(cell));
            assert_eq!(
                trial.scalar_coefficient(cell),
                baseline.scalar_coefficient(cell)
            );
            assert_eq!(
                trial.coefficient_squared(cell),
                baseline.coefficient_squared(cell)
            );
        }
        trial.commit();
        baseline.commit();
        assert_eq!(all.accepted_time(), Some(time));
    }
}

#[test]
fn g3_component_order_units_and_trace_weights_preserve_the_common_objective() {
    let (momentum, heat_config, component_config) = controls();
    let mut base = owner(4);
    let mut reordered = owner(4);
    let mut dominant = owner(1);
    let mut trace = owner(3);
    let mut equal = owner(2);
    for (time, speed, amplitude) in [(0.0, 1.0, 1.0), (0.15, 0.5, 0.2)] {
        let (rho, velocity, heat) = field(speed, amplitude);
        let components = fractions(&heat);
        let a = base
            .prepare(
                momentum,
                heat_config,
                component_config,
                time,
                &rho,
                &velocity,
                &heat,
                |k, cell| components[k][cell],
                || false,
            )
            .unwrap();
        // One permanent change of units applies to the whole bank, while
        // current amplitudes change between samples. The absent species moves.
        let b = reordered
            .prepare(
                momentum,
                heat_config,
                component_config,
                time,
                &rho,
                &velocity,
                &heat,
                |k, cell| 1e100 * components[[3, 2, 0, 1][k]][cell],
                || false,
            )
            .unwrap();
        assert_eq!(a.momentum_report(), b.momentum_report());
        assert_eq!(a.scalar_report(), b.scalar_report());
        for (a, b) in a.component_coefficients().zip(b.component_coefficients()) {
            close(a, b, 3e-12);
        }
        a.commit();
        b.commit();

        let a = dominant
            .prepare(
                momentum,
                heat_config,
                component_config,
                time,
                &rho,
                &velocity,
                &heat,
                |_, cell| heat[cell],
                || false,
            )
            .unwrap();
        let b = trace
            .prepare(
                momentum,
                heat_config,
                component_config,
                time,
                &rho,
                &velocity,
                &heat,
                |k, cell| match k {
                    0 => heat[cell],
                    1 => 1e-6 * heat[(cell + 17) % CELLS],
                    _ => 0.0,
                },
                || false,
            )
            .unwrap();
        let c = equal
            .prepare(
                momentum,
                heat_config,
                component_config,
                time,
                &rho,
                &velocity,
                &heat,
                |k, cell| heat[(cell + 17 * k) % CELLS],
                || false,
            )
            .unwrap();
        for (a, b) in a.component_coefficients().zip(b.component_coefficients()) {
            close(a, b, 2e-11);
        }
        assert!(
            a.component_coefficients()
                .zip(c.component_coefficients())
                .any(|(a, c)| (a - c).abs() > 1e-5)
        );
        a.commit();
        b.commit();
        c.commit();
    }
}

#[test]
fn g1_component_memory_relaxes_the_sum_of_physical_sources_independently() {
    let (momentum, heat_config, component_config) = controls();
    let (rho, velocity, heat) = field(1.0, 1.0);
    let mut all = owner(4);
    let original = fractions(&heat);
    all.prepare(
        momentum,
        heat_config,
        component_config,
        0.0,
        &rho,
        &velocity,
        &heat,
        |k, cell| original[k][cell],
        || false,
    )
    .unwrap()
    .commit();
    let old = [0.03, 2.0];
    all.components.history.fill(old);
    all.components.history_log_scale = Some(0.0);

    let (rho, velocity, heat) = field(0.4, 3.0);
    let rho: Vec<_> = rho.into_iter().map(|v| 1.7 * v).collect();
    let components = fractions(&heat);
    let mut flow = PeriodicDynamicSgs::new_cartesian(SHAPE, SPACING, usize::MAX).unwrap();
    let (scales, _) = flow
        .prepare_contractions(&rho, &velocity, [[0.0; 3]; 2], || false)
        .unwrap();
    let mut scalar = PeriodicDynamicScalar::new(SHAPE, usize::MAX).unwrap();
    let mut sources = vec![[0.0; 2]; CELLS];
    // Independent physical reconstruction: each scalar has its own reference
    // normalization here; multiply back to physical contractions before summing.
    for values in &components {
        let amplitude = scalar.normalize(values, &mut || false).unwrap();
        if amplitude == 0.0 {
            continue;
        }
        scalar
            .prepare_contractions(&flow, &rho, scales, &mut || false)
            .unwrap();
        let physical = scales.density.powi(2) * scales.velocity.powi(2) * amplitude.powi(2);
        for (total, pair) in sources.iter_mut().zip(&scalar.contraction) {
            for k in 0..2 {
                total[k] += physical * pair[k];
            }
        }
    }
    let dt = 0.2;
    let keep = (-dt / component_config.memory_time).exp();
    let fresh = -(-dt / component_config.memory_time).exp_m1();
    let trial = all
        .prepare(
            momentum,
            heat_config,
            component_config,
            dt,
            &rho,
            &velocity,
            &heat,
            |k, cell| components[k][cell],
            || false,
        )
        .unwrap();
    let physical = det::exp(trial.log_scale.unwrap());
    for (cell, source) in sources.iter().enumerate() {
        for k in 0..2 {
            close(
                trial.components.candidate[cell][k] * physical,
                keep * old[k] + fresh * source[k],
                5e-13,
            );
        }
    }
}

fn snapshot(
    all: &LagrangianDynamicComponents,
) -> (Option<f64>, [Option<f64>; 3], [Vec<[f64; 2]>; 3]) {
    (
        all.accepted_time(),
        [
            all.flow.flow.history_log_scale,
            all.flow.scalar.history_log_scale,
            all.components.history_log_scale,
        ],
        [
            all.flow.flow.history.clone(),
            all.flow.scalar.history.clone(),
            all.components.history.clone(),
        ],
    )
}

#[test]
fn g4_dropped_cancelled_and_refused_trials_leave_all_three_histories_unchanged() {
    let (momentum, heat_config, component_config) = controls();
    let (rho, velocity, heat) = field(1.0, 1.0);
    let mut all = owner(4);
    let initial = fractions(&heat);
    all.prepare(
        momentum,
        heat_config,
        component_config,
        0.0,
        &rho,
        &velocity,
        &heat,
        |k, cell| initial[k][cell],
        || false,
    )
    .unwrap()
    .commit();
    let accepted = snapshot(&all);
    let bytes = all.workspace_bytes();
    let (rho, velocity, heat) = field(1.4, 0.2);
    let components = fractions(&heat);
    let mut calls = 0usize;
    let trial = all
        .prepare(
            momentum,
            heat_config,
            component_config,
            0.1,
            &rho,
            &velocity,
            &heat,
            |k, cell| components[k][cell],
            || {
                calls += 1;
                false
            },
        )
        .unwrap();
    let expected: Vec<_> = trial.component_coefficients().collect();
    assert!(trial.components_report().maximum_coefficient > 0.0);
    drop(trial);
    assert_eq!(snapshot(&all), accepted);

    let refused = LagrangianComponentsConfig {
        maximum_coefficient: 0.0,
        ..component_config
    };
    assert!(matches!(
        all.prepare(
            momentum,
            heat_config,
            refused,
            0.1,
            &rho,
            &velocity,
            &heat,
            |k, cell| components[k][cell],
            || false
        ),
        Err(DynamicSgsError::Coefficient { .. })
    ));
    assert_eq!(snapshot(&all), accepted);
    let departure = LagrangianComponentsConfig {
        maximum_departure_cells: 1e-8,
        ..component_config
    };
    assert!(matches!(
        all.prepare(
            momentum,
            heat_config,
            departure,
            0.1,
            &rho,
            &velocity,
            &heat,
            |k, cell| components[k][cell],
            || false
        ),
        Err(DynamicSgsError::Cell { .. })
    ));
    assert_eq!(snapshot(&all), accepted);
    let mut polls = 0usize;
    assert!(calls > 1);
    assert!(matches!(
        all.prepare(
            momentum,
            heat_config,
            component_config,
            0.1,
            &rho,
            &velocity,
            &heat,
            |k, cell| components[k][cell],
            || {
                polls += 1;
                polls == calls - 1
            }
        ),
        Err(DynamicSgsError::Cancelled)
    ));
    assert_eq!(snapshot(&all), accepted);

    let retry = all
        .prepare(
            momentum,
            heat_config,
            component_config,
            0.1,
            &rho,
            &velocity,
            &heat,
            |k, cell| components[k][cell],
            || false,
        )
        .unwrap();
    assert_eq!(retry.component_coefficients().collect::<Vec<_>>(), expected);
    retry.commit();
    assert_eq!(all.accepted_time(), Some(0.1));
    let updated = snapshot(&all);
    for k in 0..3 {
        assert_ne!(updated.2[k], accepted.2[k]);
    }
    assert_eq!(all.workspace_bytes(), bytes);
}

#[test]
fn g0_component_count_does_not_grow_storage_and_admission_is_bounded() {
    let bytes = LagrangianDynamicComponents::required_bytes(CELLS).unwrap();
    assert_eq!(bytes, 376 * CELLS);
    for count in [1, 4, 1024] {
        let all = LagrangianDynamicComponents::new(SHAPE, SPACING, count, bytes).unwrap();
        assert_eq!(all.workspace_bytes(), bytes);
        assert_eq!(all.component_count(), count);
        assert_eq!(all.shape(), SHAPE);
        assert_eq!(all.accepted_time(), None);
        assert!(all.filter_width() > 0.0);
    }
    assert!(matches!(
        LagrangianDynamicComponents::new(SHAPE, SPACING, 4, bytes - 1),
        Err(DynamicSgsError::WorkspaceBytes { .. })
    ));
    assert!(LagrangianDynamicComponents::required_bytes(usize::MAX).is_err());
    assert!(LagrangianDynamicComponents::new(SHAPE, SPACING, 0, usize::MAX).is_err());
    assert!(LagrangianDynamicComponents::new(SHAPE, SPACING, usize::MAX, usize::MAX).is_err());
    let (momentum, heat_config, config) = controls();
    for bad in [
        LagrangianComponentsConfig {
            memory_time: 0.0,
            ..config
        },
        LagrangianComponentsConfig {
            maximum_departure_cells: 0.0,
            ..config
        },
        LagrangianComponentsConfig {
            maximum_departure_cells: 1.1,
            ..config
        },
        LagrangianComponentsConfig {
            maximum_coefficient: -1.0,
            ..config
        },
        LagrangianComponentsConfig {
            usage: SgsUse::AcceptanceRender,
            ..config
        },
    ] {
        assert!(bad.validate().is_err());
    }
    let mut all = owner(4);
    let (rho, velocity, heat) = field(1.0, 1.0);
    assert!(matches!(
        all.prepare(
            momentum,
            heat_config,
            config,
            0.0,
            &rho,
            &velocity,
            &heat,
            |k, cell| if k == 3 && cell == 100 { f64::NAN } else { 0.0 },
            || false
        ),
        Err(DynamicSgsError::Cell { .. })
    ));
    assert_eq!(all.accepted_time(), None);
}

#[test]
fn g1_component_coefficients_survive_centered_nulls_and_uniform_source_advection() {
    let (momentum, heat_config, config) = controls();
    let (rho, mut velocity, heat) = field(1.0, 1.0);
    let components = fractions(&heat);
    let mut all = owner(4);
    let trial = all
        .prepare(
            momentum,
            heat_config,
            config,
            0.0,
            &rho,
            &velocity,
            &heat,
            |k, cell| components[k][cell],
            || false,
        )
        .unwrap();
    let null_cell = trial
        .component_coefficients()
        .position(|c| c > 1e-6)
        .unwrap();
    trial.commit();
    let at = coordinates(SHAPE, null_cell);
    for axis in 0..3 {
        for direction in [1, SHAPE[axis] - 1] {
            let mut neighbor = at;
            neighbor[axis] = (neighbor[axis] + direction) % SHAPE[axis];
            velocity[index(SHAPE, neighbor)] = [0.0; 3];
        }
    }
    assert!(velocity[null_cell].iter().any(|v| *v != 0.0));
    let trial = all
        .prepare(
            momentum,
            heat_config,
            config,
            1e-6,
            &rho,
            &velocity,
            &heat,
            |k, cell| components[k][cell],
            || false,
        )
        .unwrap();
    assert_eq!(trial.component_diffusivities()[null_cell], 0.0);
    assert!(trial.component_coefficient(null_cell).unwrap() > 0.0);
    trial.commit();
    assert!(
        all.components
            .workspace
            .scalar
            .contraction
            .iter()
            .any(|p| p[1] > 0.0)
    );
    let accepted = all.components.history.clone();
    let old_log = all.components.history_log_scale.unwrap();
    let uniform = vec![[2.0, 0.0, 0.0]; CELLS];
    let trial = all
        .prepare(
            momentum,
            heat_config,
            config,
            0.400001,
            &rho,
            &uniform,
            &vec![7.0; CELLS],
            |k, _| k as f64,
            || false,
        )
        .unwrap();
    assert!(trial.components_report().maximum_coefficient > 0.0);
    assert!(trial.component_diffusivities().iter().all(|v| *v == 0.0));
    assert!(
        trial
            .components
            .workspace
            .scalar
            .contraction
            .iter()
            .all(|p| *p == [0.0; 2])
    );
    close(
        trial.log_scale.unwrap(),
        old_log - 0.4 / config.memory_time,
        3e-14,
    );
    for cell in 0..CELLS {
        let [x, y, z] = coordinates(SHAPE, cell);
        let pair = accepted[index(SHAPE, [(x + SHAPE[0] - 1) % SHAPE[0], y, z])];
        let expected = if pair[1] == 0.0 {
            0.0
        } else {
            (pair[0] / pair[1]).max(0.0)
        };
        close(trial.component_coefficient(cell).unwrap(), expected, 3e-14);
    }
}
