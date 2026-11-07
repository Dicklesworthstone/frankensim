use super::*;

const SHAPE: [usize; 3] = [5, 4, 6];
const SPACING: [f64; 3] = [0.7, 1.1, 0.9];

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

fn noise(cell: usize, salt: u64) -> f64 {
    let mut h = ((cell as u64 + 1).wrapping_mul(0x9e37_79b9_7f4a_7c15)) ^ salt;
    h = (h ^ (h >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    h = (h ^ (h >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    h ^= h >> 31;
    (h >> 11) as f64 / 9_007_199_254_740_992.0 * 2.0 - 1.0
}

fn fields(normal: usize) -> (Vec<f64>, Vec<[f64; 3]>, Vec<f64>) {
    let mut rho = Vec::new();
    let mut velocity = Vec::new();
    let mut scalar = Vec::new();
    for cell in 0..SHAPE.iter().product() {
        let p = coordinates(SHAPE, cell);
        let eta = (p[normal] as f64 + 0.5) / SHAPE[normal] as f64;
        let envelope = 4.0 * eta * (1.0 - eta);
        rho.push(1.0 + 0.15 * noise(cell, 911));
        velocity.push(core::array::from_fn(|axis| {
            0.02 * envelope * noise(cell, 31 * axis as u64 + 7)
        }));
        scalar.push(1.0 - eta + 0.05 * envelope * noise(cell, 101));
    }
    (rho, velocity, scalar)
}

fn walls(normal: usize) -> [[f64; 3]; 2] {
    let mut walls = [[0.0; 3]; 2];
    walls[0][(normal + 1) % 3] = 0.01;
    walls[1][(normal + 1) % 3] = -0.02;
    walls
}

fn species(heat: &[f64], component: usize, cell: usize) -> f64 {
    let first = 0.4 + 0.1 * heat[cell] + 0.02 * noise(cell, 291);
    match component {
        0 => first,
        1 => 1.0 - first,
        _ => 0.0,
    }
}

fn owner(normal: usize) -> LagrangianDynamicComponents {
    LagrangianDynamicComponents::new_channel(SHAPE, SPACING, normal, 3, usize::MAX).unwrap()
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-13 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

#[test]
fn g1_channel_interpolation_preserves_wall_separation_and_tangential_periodicity() {
    for normal in 0..3 {
        let tangent = (normal + 1) % 3;
        let history: Vec<_> = (0..SHAPE.iter().product())
            .map(|cell| {
                let p = coordinates(SHAPE, cell);
                [10.0 * p[normal] as f64 + p[tangent] as f64, 1.0]
            })
            .collect();
        for (layer, outward) in [(0, 0.5), (SHAPE[normal] - 1, -0.5)] {
            let mut p = [1; 3];
            p[normal] = layer;
            p[tangent] = 0;
            let cell = index(SHAPE, p);
            let mut departure = [0.0; 3];
            departure[normal] = outward;
            departure[tangent] = 1.0;
            let value = sample_history(
                FilterGeometry::Channel(normal),
                SHAPE,
                cell,
                departure,
                &history,
            )
            .unwrap();
            assert_eq!(
                value,
                [10.0 * layer as f64 + (SHAPE[tangent] - 1) as f64, 1.0]
            );
            departure[normal] = outward * (1.0 + 1e-12);
            assert!(matches!(
                sample_history(FilterGeometry::Channel(normal), SHAPE, cell, departure, &history),
                Err(DynamicSgsError::Cell { cell: at, reason:
                    "Lagrangian departure crosses impermeable channel wall" }) if at == cell
            ));
        }
        let cell = index(SHAPE, [2; 3]);
        let mut departure = [0.0; 3];
        departure[normal] = 0.25;
        departure[tangent] = -0.5;
        assert_eq!(
            sample_history(
                FilterGeometry::Channel(normal),
                SHAPE,
                cell,
                departure,
                &history
            )
            .unwrap(),
            [20.0, 1.0]
        );
        // A periodic coordinate can round to exactly the extent after modulo.
        let mut departure = [0.0; 3];
        departure[tangent] = 1e-300;
        assert_eq!(
            sample_history(
                FilterGeometry::Channel(normal),
                SHAPE,
                0,
                departure,
                &history
            )
            .unwrap(),
            history[0]
        );
    }
}

// Explicit two-axis sum, independent of FilterGeometry's sample traversal.
// This intentionally differs from the instantaneous channel's whole-plane fit.
fn local_pair(normal: usize, cell: usize, sources: &[[f64; 2]]) -> [f64; 2] {
    let tangent: Vec<_> = (0..3).filter(|axis| *axis != normal).collect();
    let at = coordinates(SHAPE, cell);
    let mut pair = [0.0; 2];
    for a in 0..3 {
        for b in 0..3 {
            let mut p = at;
            p[tangent[0]] = (p[tangent[0]] + SHAPE[tangent[0]] + a - 1) % SHAPE[tangent[0]];
            p[tangent[1]] = (p[tangent[1]] + SHAPE[tangent[1]] + b - 1) % SHAPE[tangent[1]];
            let source = sources[index(SHAPE, p)];
            pair[0] += source[0] / 9.0;
            pair[1] += source[1] / 9.0;
        }
    }
    pair
}

#[test]
fn g1_channel_initial_joint_fit_uses_local_unaveraged_physical_sources() {
    let (momentum, heat_config, component_config) = controls();
    for normal in 0..3 {
        let (rho, velocity, heat) = fields(normal);
        let wall_velocity = walls(normal);
        let heat_walls = [ScalarWall::Value(1.2), ScalarWall::Value(-0.1)];
        let geometry = FilterGeometry::Channel(normal);
        let mut flow =
            PeriodicDynamicSgs::new_cartesian_geometry(SHAPE, SPACING, geometry, usize::MAX)
                .unwrap();
        let (scales, _) = flow
            .prepare_contractions(&rho, &velocity, wall_velocity, || false)
            .unwrap();
        let mut scalar = PeriodicDynamicScalar::new_geometry(SHAPE, geometry, usize::MAX).unwrap();
        scalar
            .normalize_with_walls(&heat, Some(heat_walls), &mut || false)
            .unwrap();
        scalar
            .prepare_contractions(&flow, &rho, scales, &mut || false)
            .unwrap();
        let mut components =
            PeriodicDynamicComponents::new_geometry(SHAPE, geometry, 3, usize::MAX).unwrap();
        let values = |k, cell| species(&heat, k, cell);
        let scale = components.measure_scale(&values, &mut || false).unwrap();
        components
            .prepare_contractions(&flow, &rho, scales, scale, &values, &mut || false)
            .unwrap();
        let mut history = owner(normal);
        let trial = history
            .prepare_with_walls(
                momentum,
                heat_config,
                component_config,
                0.0,
                &rho,
                &velocity,
                wall_velocity,
                &heat,
                heat_walls,
                values,
                || false,
            )
            .unwrap();
        assert!(trial.components_report().maximum_coefficient > 0.0);
        let sources = [
            &flow.contraction,
            &scalar.contraction,
            &components.scalar.contraction,
        ];
        for cell in 0..rho.len() {
            let actual = [
                trial.coefficient_squared(cell).unwrap(),
                trial.scalar_coefficient(cell).unwrap(),
                trial.component_coefficient(cell).unwrap(),
            ];
            for (fit, source) in actual.into_iter().zip(sources) {
                let pair = local_pair(normal, cell, source);
                close(
                    fit,
                    if pair[1] == 0.0 {
                        0.0
                    } else {
                        (pair[0] / pair[1]).max(0.0)
                    },
                );
            }
        }
        // Distinct same-plane coefficients prove the component objective was
        // not prematurely collapsed by the instantaneous plane-average helper.
        let first = trial.component_coefficient(0).unwrap();
        assert!((0..rho.len()).any(|cell| {
            coordinates(SHAPE, cell)[normal] == 0
                && (trial.component_coefficient(cell).unwrap() - first).abs() > 1e-8
        }));
        trial.commit();
        assert_eq!(history.accepted_time(), Some(0.0));
    }
}

#[derive(Debug, PartialEq)]
struct Accepted {
    momentum: Vec<[f64; 2]>,
    heat: Vec<[f64; 2]>,
    components: Vec<[f64; 2]>,
    scales: [Option<f64>; 3],
    time: Option<f64>,
}

fn accepted(owner: &LagrangianDynamicComponents) -> Accepted {
    Accepted {
        momentum: owner.flow.flow.history.clone(),
        heat: owner.flow.scalar.history.clone(),
        components: owner.components.history.clone(),
        scales: [
            owner.flow.flow.history_log_scale,
            owner.flow.scalar.history_log_scale,
            owner.components.history_log_scale,
        ],
        time: owner.accepted_time(),
    }
}

#[test]
fn g4_channel_refused_cancelled_and_dropped_joint_trials_preserve_retry() {
    let (momentum, heat_config, component_config) = controls();
    let normal = 1;
    let (rho, velocity, heat) = fields(normal);
    let heat_walls = [ScalarWall::Value(1.0), ScalarWall::Interior];
    let mut history = owner(normal);
    let mut reference = owner(normal);
    for filter in [&mut history, &mut reference] {
        filter
            .prepare_with_walls(
                momentum,
                heat_config,
                component_config,
                0.0,
                &rho,
                &velocity,
                walls(normal),
                &heat,
                heat_walls,
                |k, cell| species(&heat, k, cell),
                || false,
            )
            .unwrap()
            .commit();
    }
    let saved = accepted(&history);
    let mut crossing = vec![[0.0; 3]; rho.len()];
    for u in &mut crossing {
        u[normal] = 0.6 * SPACING[normal];
    }
    assert!(matches!(
        history.prepare_with_walls(
            momentum,
            heat_config,
            component_config,
            1.0,
            &rho,
            &crossing,
            walls(normal),
            &heat,
            heat_walls,
            |k, cell| species(&heat, k, cell),
            || false,
        ),
        Err(DynamicSgsError::Cell {
            reason: "Lagrangian departure crosses impermeable channel wall",
            ..
        })
    ));
    assert_eq!(accepted(&history), saved);
    let mut polls = 0;
    drop(
        history
            .prepare_with_walls(
                momentum,
                heat_config,
                component_config,
                0.1,
                &rho,
                &velocity,
                walls(normal),
                &heat,
                heat_walls,
                |k, cell| species(&heat, k, cell),
                || {
                    polls += 1;
                    false
                },
            )
            .unwrap(),
    );
    assert_eq!(accepted(&history), saved);
    let mut seen = 0;
    assert!(matches!(
        history.prepare_with_walls(
            momentum,
            heat_config,
            component_config,
            0.1,
            &rho,
            &velocity,
            walls(normal),
            &heat,
            heat_walls,
            |k, cell| species(&heat, k, cell),
            || {
                seen += 1;
                seen >= polls - 1
            },
        ),
        Err(DynamicSgsError::Cancelled)
    ));
    assert_eq!(accepted(&history), saved);
    assert!(matches!(
        history.prepare_with_walls(
            momentum,
            heat_config,
            LagrangianComponentsConfig {
                maximum_coefficient: 0.0,
                ..component_config
            },
            0.1,
            &rho,
            &velocity,
            walls(normal),
            &heat,
            heat_walls,
            |k, cell| species(&heat, k, cell),
            || false,
        ),
        Err(DynamicSgsError::Coefficient { .. })
    ));
    assert_eq!(accepted(&history), saved);
    let retry = history
        .prepare_with_walls(
            momentum,
            heat_config,
            component_config,
            0.1,
            &rho,
            &velocity,
            walls(normal),
            &heat,
            heat_walls,
            |k, cell| species(&heat, k, cell),
            || false,
        )
        .unwrap();
    let expected = reference
        .prepare_with_walls(
            momentum,
            heat_config,
            component_config,
            0.1,
            &rho,
            &velocity,
            walls(normal),
            &heat,
            heat_walls,
            |k, cell| species(&heat, k, cell),
            || false,
        )
        .unwrap();
    assert_eq!(retry.momentum_report(), expected.momentum_report());
    assert_eq!(retry.scalar_report(), expected.scalar_report());
    assert_eq!(retry.components_report(), expected.components_report());
    assert_eq!(retry.viscosities(), expected.viscosities());
    assert_eq!(
        retry.scalar_diffusivities(),
        expected.scalar_diffusivities()
    );
    assert_eq!(
        retry.component_diffusivities(),
        expected.component_diffusivities()
    );
    retry.commit();
    expected.commit();
    assert_eq!(accepted(&history), accepted(&reference));
}

#[test]
fn g1_channel_zero_sources_advect_signed_memory_in_the_moving_wall_frame() {
    let (momentum, heat_config, component_config) = controls();
    for normal in 0..3 {
        let (rho, velocity, heat) = fields(normal);
        let mut history =
            LagrangianDynamicComponents::new_channel(SHAPE, [1.0; 3], normal, 3, usize::MAX)
                .unwrap();
        history
            .prepare_with_walls(
                momentum,
                heat_config,
                component_config,
                0.0,
                &rho,
                &velocity,
                walls(normal),
                &heat,
                [ScalarWall::Interior; 2],
                |k, cell| species(&heat, k, cell),
                || false,
            )
            .unwrap()
            .commit();
        let tangent = (normal + 1) % 3;
        for cell in 0..rho.len() {
            let p = coordinates(SHAPE, cell);
            let value = 0.01 * (p[tangent] as f64 - 2.0 + p[normal] as f64);
            history.flow.flow.history[cell] = [value, 1.0];
            history.flow.scalar.history[cell] = [2.0 * value, 2.0];
            history.components.history[cell] = [3.0 * value, 3.0];
        }
        history.flow.flow.history_log_scale = Some(0.0);
        history.flow.scalar.history_log_scale = Some(0.0);
        history.components.history_log_scale = Some(0.0);
        let old = accepted(&history);
        let mut translation = [0.0; 3];
        translation[tangent] = 1.0;
        let trial = history
            .prepare_with_walls(
                momentum,
                heat_config,
                component_config,
                1.0,
                &vec![1.0; rho.len()],
                &vec![translation; rho.len()],
                [translation; 2],
                &vec![7.0; rho.len()],
                [ScalarWall::Value(7.0); 2],
                |k, _| [0.25, 0.75, 0.0][k],
                || false,
            )
            .unwrap();
        assert!(trial.viscosities().iter().all(|value| *value == 0.0));
        assert!(
            trial
                .scalar_diffusivities()
                .iter()
                .all(|value| *value == 0.0)
        );
        assert!(
            trial
                .component_diffusivities()
                .iter()
                .all(|value| *value == 0.0)
        );
        assert!(trial.momentum_report().minimum_raw_coefficient_squared < 0.0);
        assert!(trial.scalar_report().minimum_raw_coefficient < 0.0);
        assert!(trial.components_report().minimum_raw_coefficient < 0.0);
        for (source, candidate, log, memory) in [
            (
                &old.momentum,
                &trial.joint.momentum.owner.candidate,
                trial.joint.momentum.log_scale,
                momentum.memory_time,
            ),
            (
                &old.heat,
                &trial.joint.scalar.candidate,
                trial.joint.log_scale,
                heat_config.memory_time,
            ),
            (
                &old.components,
                &trial.components.candidate,
                trial.log_scale,
                component_config.memory_time,
            ),
        ] {
            let scale = det::exp(log.unwrap());
            let decay = (-1.0 / memory).exp();
            for cell in 0..rho.len() {
                let mut p = coordinates(SHAPE, cell);
                p[tangent] = (p[tangent] + SHAPE[tangent] - 1) % SHAPE[tangent];
                let upstream = index(SHAPE, p);
                for k in 0..2 {
                    close(candidate[cell][k] * scale, decay * source[upstream][k]);
                }
            }
        }
        trial.commit();
        assert_eq!(history.accepted_time(), Some(1.0));
    }
}

#[test]
fn g0_channel_wall_geometry_and_storage_admission_remain_explicit() {
    let cells = SHAPE.iter().product();
    let momentum_bytes = LagrangianDynamicSgs::required_bytes(cells).unwrap();
    let heat_bytes = LagrangianDynamicScalar::required_bytes(cells).unwrap();
    let joint_bytes = LagrangianDynamicComponents::required_bytes(cells).unwrap();
    assert_eq!(
        (momentum_bytes, heat_bytes, joint_bytes),
        (168 * cells, 264 * cells, 376 * cells)
    );
    for normal in 0..3 {
        let flow =
            LagrangianDynamicSgs::new_channel(SHAPE, SPACING, normal, momentum_bytes).unwrap();
        let scalar =
            LagrangianDynamicScalar::new_channel(SHAPE, SPACING, normal, heat_bytes).unwrap();
        let joint =
            LagrangianDynamicComponents::new_channel(SHAPE, SPACING, normal, 1000, joint_bytes)
                .unwrap();
        assert_eq!(
            (
                flow.normal_axis(),
                scalar.normal_axis(),
                joint.normal_axis()
            ),
            (Some(normal), Some(normal), Some(normal))
        );
        assert_eq!(
            (
                flow.workspace_bytes(),
                scalar.workspace_bytes(),
                joint.workspace_bytes()
            ),
            (momentum_bytes, heat_bytes, joint_bytes)
        );
        assert!(
            LagrangianDynamicComponents::new_channel(SHAPE, SPACING, normal, 3, joint_bytes - 1)
                .is_err()
        );
    }
    assert!(LagrangianDynamicSgs::new_channel(SHAPE, SPACING, 3, usize::MAX).is_err());
    assert!(
        LagrangianDynamicComponents::new_channel([3, 2, 3], SPACING, 1, 3, usize::MAX).is_err()
    );
    let (momentum, scalar_config, _) = controls();
    let mut two = LagrangianDynamicScalar::new_channel([3, 2, 3], SPACING, 1, usize::MAX).unwrap();
    let rho = vec![1.0; 18];
    let u = vec![[0.0; 3]; 18];
    let heat = vec![7.0; 18];
    assert!(
        two.prepare(momentum, scalar_config, 0.0, &rho, &u, &heat, || false)
            .is_err()
    );
    for policy in [ScalarWall::ZeroGradient, ScalarWall::Interior] {
        assert!(
            two.prepare_with_walls(
                momentum,
                scalar_config,
                0.0,
                &rho,
                &u,
                [[0.0; 3]; 2],
                &heat,
                [policy; 2],
                || false,
            )
            .is_err()
        );
    }
    two.prepare_with_walls(
        momentum,
        scalar_config,
        0.0,
        &rho,
        &u,
        [[0.0; 3]; 2],
        &heat,
        [ScalarWall::Value(7.0); 2],
        || false,
    )
    .unwrap()
    .commit();
    let mut periodic = LagrangianDynamicSgs::new(SHAPE, SPACING, usize::MAX).unwrap();
    assert_eq!(periodic.normal_axis(), None);
    assert!(
        periodic
            .prepare_with_walls(
                momentum,
                0.0,
                &vec![1.0; cells],
                &vec![[0.0; 3]; cells],
                [[0.0; 3]; 2],
                || false,
            )
            .is_err()
    );
}

#[test]
fn g3_channel_scalar_walls_and_all_history_units_follow_original_sources() {
    let (momentum, heat_config, component_config) = controls();
    let normal = 1;
    let (rho, velocity, heat) = fields(normal);
    let wall_velocity = walls(normal);
    let heat_walls = [ScalarWall::Value(1.2), ScalarWall::ZeroGradient];
    let mut reference = owner(normal);
    let length = 1e30;
    let thermal = 1e50;
    let mut rescaled = LagrangianDynamicComponents::new_channel(
        SHAPE,
        SPACING.map(|v| v * length),
        normal,
        3,
        usize::MAX,
    )
    .unwrap();
    let density: Vec<_> = rho.iter().map(|v| v * 1e40).collect();
    let speed: Vec<_> = velocity.iter().map(|u| u.map(|v| v * length)).collect();
    let temperature: Vec<_> = heat.iter().map(|v| v * thermal).collect();
    for time in [0.0, 0.1] {
        let a = reference
            .prepare_with_walls(
                momentum,
                heat_config,
                component_config,
                time,
                &rho,
                &velocity,
                wall_velocity,
                &heat,
                heat_walls,
                |k, cell| species(&heat, k, cell),
                || false,
            )
            .unwrap();
        let b = rescaled
            .prepare_with_walls(
                momentum,
                heat_config,
                component_config,
                time,
                &density,
                &speed,
                wall_velocity.map(|u| u.map(|v| v * length)),
                &temperature,
                [ScalarWall::Value(1.2 * thermal), ScalarWall::ZeroGradient],
                |k, cell| species(&heat, k, cell) * 1e20,
                || false,
            )
            .unwrap();
        for (a, b) in a.coefficients_squared().zip(b.coefficients_squared()) {
            close(a, b);
        }
        for (a, b) in a.scalar_coefficients().zip(b.scalar_coefficients()) {
            close(a, b);
        }
        for (a, b) in a.component_coefficients().zip(b.component_coefficients()) {
            close(a, b);
        }
        for (a, b) in a
            .component_diffusivities()
            .iter()
            .zip(b.component_diffusivities())
        {
            close(*a, b / (length * length));
        }
        a.commit();
        b.commit();
    }
    // Uniform fluid values still have a physical wall-driven gradient. Its
    // nonzero denominator history must survive even with zero scalar covariance.
    let mut scalar =
        LagrangianDynamicScalar::new_channel(SHAPE, SPACING, normal, usize::MAX).unwrap();
    scalar
        .prepare_with_walls(
            momentum,
            heat_config,
            0.0,
            &rho,
            &velocity,
            wall_velocity,
            &vec![3.0; rho.len()],
            [ScalarWall::Value(1.0), ScalarWall::Value(5.0)],
            || false,
        )
        .unwrap()
        .commit();
    assert!(scalar.scalar.history_log_scale.is_some());
    assert!(scalar.scalar.history.iter().any(|pair| pair[1] > 0.0));
    assert!(scalar.scalar.history.iter().all(|pair| pair[0] == 0.0));
}
