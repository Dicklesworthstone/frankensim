use super::*;

fn momentum_config() -> LagrangianSgsConfig {
    LagrangianSgsConfig {
        memory_time: 1.5,
        maximum_departure_cells: 1.0,
        maximum_coefficient_squared: 100.0,
        maximum_eddy_to_molecular: 1000.0,
        usage: SgsUse::Exploration,
    }
}

fn scalar_config() -> LagrangianScalarConfig {
    LagrangianScalarConfig {
        memory_time: 0.4,
        maximum_departure_cells: 1.0,
        maximum_coefficient: 100.0,
        maximum_diffusivity_to_molecular: 1000.0,
        usage: SgsUse::Exploration,
    }
}

fn field(shape: [usize; 3], speed: f64, amplitude: f64) -> (Vec<f64>, Vec<[f64; 3]>, Vec<f64>) {
    let mut rho = Vec::new();
    let mut u = Vec::new();
    let mut phi = Vec::new();
    for z in 0..shape[2] {
        for y in 0..shape[1] {
            for x in 0..shape[0] {
                let [x, y, z] = core::array::from_fn::<_, 3, _>(|axis| {
                    core::f64::consts::TAU * ([x, y, z][axis] as f64 + 0.5) / shape[axis] as f64
                });
                rho.push(1.2 + 0.2 * x.sin() + 0.1 * y.cos());
                u.push([
                    speed * (0.3 * x.sin() + 0.2 * (y + z).sin() + 0.1 * (2.0 * x - y).cos()),
                    speed * (0.2 * y.cos() + 0.3 * (z + x).sin() + 0.1 * (x - 2.0 * z).sin()),
                    speed * (0.2 * z.sin() + 0.1 * (x + y).cos() + 0.3 * (2.0 * y - z).sin()),
                ]);
                phi.push(
                    amplitude
                        * (0.7 * (x + y).sin()
                            + 0.9 * (2.0 * y - z).cos()
                            + 0.4 * (x - 2.0 * z).sin()),
                );
            }
        }
    }
    (rho, u, phi)
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

#[test]
fn g3_native_refresh_preserves_history_and_public_time_admission() {
    let shape = [5, 6, 7];
    let (rho, velocity, scalar) = field(shape, 0.05, 2.0);
    let mut owner = LagrangianDynamicScalar::new(shape, [1.0; 3], usize::MAX).unwrap();
    owner.prepare(momentum_config(), scalar_config(), 0.0,
        &rho, &velocity, &scalar, || false).unwrap().commit();
    let momentum_history = owner.flow.history.clone();
    let scalar_history = owner.scalar.history.clone();
    let scales = (owner.flow.history_log_scale, owner.scalar.history_log_scale);
    let zero_velocity = vec![[0.0; 3]; rho.len()];
    let zero_scalar = vec![0.0; rho.len()];
    let refreshed = owner.prepare_or_refresh(
        momentum_config(), scalar_config(), 0.0, &rho, &zero_velocity, None,
        &zero_scalar, None, || false,
    ).unwrap();
    assert!(refreshed.viscosities().iter().all(|&value| value == 0.0));
    assert!(refreshed.diffusivities().iter().all(|&value| value == 0.0));
    drop(refreshed);
    assert_eq!(owner.flow.history, momentum_history);
    assert_eq!(owner.scalar.history, scalar_history);
    assert_eq!((owner.flow.history_log_scale, owner.scalar.history_log_scale), scales);
    assert_eq!(owner.accepted_time(), Some(0.0));
    assert!(owner.prepare(momentum_config(), scalar_config(), 0.0,
        &rho, &velocity, &scalar, || false).is_err());
}

#[test]
fn g3_native_history_fork_excludes_rejected_scratch_and_evolves_independently() {
    let shape = [5, 6, 7];
    let (rho, velocity, scalar) = field(shape, 0.05, 2.0);
    let mut owner = LagrangianDynamicScalar::new(shape, [0.7, 1.1, 0.9], usize::MAX)
        .unwrap().with_adaptive_momentum_memory(AdaptiveLagrangianMemory {
            coefficient: 0.15, minimum_time: 0.05,
        }).unwrap();
    owner.prepare(momentum_config(), scalar_config(), 0.0,
        &rho, &velocity, &scalar, || false).unwrap().commit();
    let (_, other_velocity, other_scalar) = field(shape, 0.1, 0.5);
    drop(owner.prepare(momentum_config(), scalar_config(), 0.2,
        &rho, &other_velocity, &other_scalar, || false).unwrap());
    assert!(owner.fork_for_trial(owner.workspace_bytes() - 1).is_err());
    let mut fork = owner.fork_for_trial(owner.workspace_bytes()).unwrap();
    assert_eq!(fork.flow.history, owner.flow.history);
    assert_eq!(fork.scalar.history, owner.scalar.history);
    assert_eq!(fork.flow.flow.spacing, [0.7, 1.1, 0.9]);
    let original = owner.prepare(momentum_config(), scalar_config(), 0.4,
        &rho, &velocity, &scalar, || false).unwrap();
    let copied = fork.prepare(momentum_config(), scalar_config(), 0.4,
        &rho, &velocity, &scalar, || false).unwrap();
    assert_eq!(original.viscosities(), copied.viscosities());
    assert_eq!(original.diffusivities(), copied.diffusivities());
    assert_eq!(original.momentum_memory_time_range(), copied.momentum_memory_time_range());
    copied.commit();
    drop(original);
    assert_eq!(owner.accepted_time(), Some(0.0));
    assert_eq!(fork.accepted_time(), Some(0.4));
}

#[test]
fn g0_initial_joint_fit_matches_original_local_box_momentum_and_scalar() {
    let shape = [5, 6, 7];
    let spacing = [0.7, 1.1, 0.9];
    let (rho, u, phi) = field(shape, 1.0, 2.0);
    let mut instantaneous = PeriodicDynamicSgs::new_cartesian(shape, spacing, usize::MAX).unwrap();
    let mut scalar = PeriodicDynamicScalar::new(shape, usize::MAX).unwrap();
    let expected = instantaneous
        .prepare_with_scalar(
            DynamicSgsConfig {
                averaging: DynamicAveraging::LocalBox,
                maximum_coefficient_squared: momentum_config().maximum_coefficient_squared,
                maximum_eddy_to_molecular: momentum_config().maximum_eddy_to_molecular,
                usage: SgsUse::Exploration,
            },
            DynamicScalarConfig {
                averaging: DynamicAveraging::LocalBox,
                maximum_coefficient: scalar_config().maximum_coefficient,
                maximum_diffusivity_to_molecular: scalar_config().maximum_diffusivity_to_molecular,
                usage: SgsUse::Exploration,
            },
            &rho,
            &u,
            &phi,
            &mut scalar,
            || false,
        )
        .unwrap();
    let mut history = LagrangianDynamicScalar::new(shape, spacing, usize::MAX).unwrap();
    assert_eq!(history.shape(), shape);
    assert_eq!(history.filter_width(), instantaneous.filter_width());
    assert_eq!(history.accepted_time(), None);
    let trial = history
        .prepare(
            momentum_config(),
            scalar_config(),
            0.0,
            &rho,
            &u,
            &phi,
            || false,
        )
        .unwrap();
    assert_eq!(trial.momentum_report(), expected.0);
    assert_eq!(trial.scalar_report(), expected.1);
    assert_eq!(trial.viscosities(), instantaneous.viscosities());
    assert_eq!(trial.diffusivities(), scalar.diffusivities());
    assert!(expected.1.maximum_diffusivity > 0.0);
    assert_eq!(trial.coefficients_squared().len(), rho.len());
    assert_eq!(trial.scalar_coefficients().len(), rho.len());
    for (cell, (stress, scalar)) in trial
        .coefficients_squared()
        .zip(trial.scalar_coefficients())
        .enumerate()
    {
        assert_eq!(trial.coefficient_squared(cell), Some(stress));
        assert_eq!(trial.scalar_coefficient(cell), Some(scalar));
        assert_eq!(Some(stress), instantaneous.coefficient_squared(cell));
    }
    assert_eq!(trial.coefficient_squared(rho.len()), None);
    assert_eq!(trial.scalar_coefficient(rho.len()), None);
    trial.commit();
    assert_eq!(history.accepted_time(), Some(0.0));
}

#[test]
fn g0_scalar_relaxation_uses_its_own_physical_contractions_and_memory_time() {
    let shape = [5; 3];
    let spacing = [0.8, 1.0, 1.2];
    let (rho, u, phi) = field(shape, 0.7, 3.0);
    // The independently owned instantaneous filter supplies the original
    // Germano source; the time-integration oracle uses direct physical powers,
    // not logarithmic normalization or the production history helper.
    let mut flow = PeriodicDynamicSgs::new_cartesian(shape, spacing, usize::MAX).unwrap();
    let (scales, _) = flow
        .prepare_contractions(&rho, &u, [[0.0; 3]; 2], || false)
        .unwrap();
    let mut scalar = PeriodicDynamicScalar::new(shape, usize::MAX).unwrap();
    let phi_scale = scalar.normalize(&phi, &mut || false).unwrap();
    scalar
        .prepare_contractions(&flow, &rho, scales, &mut || false)
        .unwrap();
    let physical_scale = scales.density.powi(2) * scales.velocity.powi(2) * phi_scale.powi(2);
    let sources: Vec<_> = scalar
        .contraction
        .iter()
        .map(|p| p.map(|v| v * physical_scale))
        .collect();

    let mut history = LagrangianDynamicScalar::new(shape, spacing, usize::MAX).unwrap();
    history
        .prepare(
            momentum_config(),
            scalar_config(),
            0.0,
            &rho,
            &u,
            &phi,
            || false,
        )
        .unwrap()
        .commit();
    let old = [0.03, 2.0];
    history.scalar.history.fill(old);
    history.scalar.history_log_scale = Some(0.0);
    let dt = 0.25;
    let keep = (-dt / scalar_config().memory_time).exp();
    let fresh = -(-dt / scalar_config().memory_time).exp_m1();
    let trial = history
        .prepare(
            momentum_config(),
            scalar_config(),
            dt,
            &rho,
            &u,
            &phi,
            || false,
        )
        .unwrap();
    let scale = det::exp(trial.log_scale.unwrap());
    for (cell, source) in sources.iter().enumerate() {
        for component in 0..2 {
            let expected = keep * old[component] + fresh * source[component];
            let actual = trial.scalar.candidate[cell][component] * scale;
            close(actual, expected, 3e-14);
        }
        let pair = trial.scalar.candidate[cell];
        close(
            trial.scalar_coefficient(cell).unwrap(),
            (pair[0] / pair[1]).max(0.0),
            1e-15,
        );
    }
    trial.commit();
    assert_eq!(history.accepted_time(), Some(dt));
}

#[test]
fn g0_uniform_sources_advect_signed_history_and_retain_long_decay() {
    let shape = [5, 4, 3];
    let cells = 60;
    let mut history = LagrangianDynamicScalar::new(shape, [0.5, 1.0, 2.0], usize::MAX).unwrap();
    history.flow.accepted_time = Some(0.0);
    history.flow.history_log_scale = Some(0.0);
    history.flow.history.fill([0.02, 1.0]);
    history.scalar.history_log_scale = Some(0.0);
    for (cell, pair) in history.scalar.history.iter_mut().enumerate() {
        let [x, y, z] = coordinates(shape, cell);
        *pair = [0.01 * (x as f64 - 2.0 + y as f64 - z as f64), 1.0];
    }
    let original = history.scalar.history.clone();
    let rho = vec![1.0; cells];
    // One exact periodic x-cell departure. The uniform Galilean translation
    // has no new stress/scalar contractions but still transports old memory.
    let u = vec![[1.0, 0.0, 0.0]; cells];
    let phi = vec![12.0; cells];
    let config = LagrangianScalarConfig {
        memory_time: 0.0125,
        ..scalar_config()
    };
    let trial = history
        .prepare(momentum_config(), config, 0.5, &rho, &u, &phi, || false)
        .unwrap();
    assert!(trial.viscosities().iter().all(|value| *value == 0.0));
    assert!(trial.diffusivities().iter().all(|value| *value == 0.0));
    assert!(trial.scalar_report().backscatter_limited_cells > 0);
    assert!(trial.scalar_report().maximum_coefficient > 0.0);
    close(trial.log_scale.unwrap(), -40.0, 2e-15);
    for cell in 0..cells {
        let [x, y, z] = coordinates(shape, cell);
        let upstream = index(shape, [(x + 4) % 5, y, z]);
        let pair = trial.scalar.candidate[cell];
        close(pair[0] / pair[1], original[upstream][0], 1e-15);
        close(
            trial.scalar_coefficient(cell).unwrap(),
            original[upstream][0].max(0.0),
            1e-15,
        );
        assert!(pair[1] > 0.0);
    }
    trial.commit();
    let (rho, u, phi) = field(shape, 1e-5, 1e-4);
    let continued = history
        .prepare(momentum_config(), config, 0.501, &rho, &u, &phi, || false)
        .unwrap();
    assert!(continued.scalar_report().maximum_diffusivity > 0.0);
    assert!(continued.scalar.candidate.iter().all(|pair| pair[1] > 0.0));
}

#[test]
fn g3_joint_memory_preserves_units_and_original_momentum_trajectory() {
    let shape = [5; 3];
    let sequence = [(0.0, 1.0, 1.0), (0.1, 0.4, 0.05), (0.25, 0.7, 2.0)];
    let mut base = LagrangianDynamicScalar::new(shape, [1.0; 3], usize::MAX).unwrap();
    let mut momentum = LagrangianDynamicSgs::new(shape, [1.0; 3], usize::MAX).unwrap();
    let mut expected = Vec::new();
    for (time, speed, amplitude) in sequence {
        let (rho, u, phi) = field(shape, speed, amplitude);
        let a = base
            .prepare(
                momentum_config(),
                scalar_config(),
                time,
                &rho,
                &u,
                &phi,
                || false,
            )
            .unwrap();
        let b = momentum
            .prepare(momentum_config(), time, &rho, &u, || false)
            .unwrap();
        assert_eq!(a.momentum_report(), b.report());
        assert_eq!(a.viscosities(), b.viscosities());
        assert!(a.coefficients_squared().eq(b.coefficients_squared()));
        expected.push((
            a.viscosities().to_vec(),
            a.diffusivities().to_vec(),
            a.coefficients_squared().collect::<Vec<_>>(),
            a.scalar_coefficients().collect::<Vec<_>>(),
        ));
        a.commit();
        b.commit();
    }
    for (length, density_scale, scalar_scale) in [(1e-40, 1e40, -1e-60), (1e40, 1e-40, 1e60)] {
        let mut scaled = LagrangianDynamicScalar::new(shape, [length; 3], usize::MAX).unwrap();
        for (step, (time, speed, amplitude)) in sequence.into_iter().enumerate() {
            let (rho, u, phi) = field(shape, speed, amplitude);
            let rho: Vec<_> = rho.into_iter().map(|v| v * density_scale).collect();
            let u: Vec<_> = u.into_iter().map(|v| v.map(|c| c * length)).collect();
            let phi: Vec<_> = phi.into_iter().map(|v| (25.0 + v) * scalar_scale).collect();
            let trial = scaled
                .prepare(
                    momentum_config(),
                    scalar_config(),
                    time,
                    &rho,
                    &u,
                    &phi,
                    || false,
                )
                .unwrap();
            for (actual, expected) in trial.viscosities().iter().zip(&expected[step].0) {
                close(actual / (length * length), *expected, 3e-10);
            }
            for (actual, expected) in trial.diffusivities().iter().zip(&expected[step].1) {
                close(actual / (length * length), *expected, 3e-10);
            }
            for (actual, expected) in trial.coefficients_squared().zip(&expected[step].2) {
                close(actual, *expected, 3e-10);
            }
            for (actual, expected) in trial.scalar_coefficients().zip(&expected[step].3) {
                close(actual, *expected, 3e-10);
            }
            trial.commit();
        }
    }
}

#[test]
fn g4_cancelled_refused_and_dropped_joint_trials_preserve_both_histories() {
    let shape = [5; 3];
    let (rho, u, phi) = field(shape, 1.0, 1.0);
    let mut history = LagrangianDynamicScalar::new(shape, [1.0; 3], usize::MAX).unwrap();
    history
        .prepare(
            momentum_config(),
            scalar_config(),
            0.0,
            &rho,
            &u,
            &phi,
            || false,
        )
        .unwrap()
        .commit();
    let accepted_momentum = history.flow.history.clone();
    let accepted_scalar = history.scalar.history.clone();
    let old_scales = (
        history.flow.history_log_scale,
        history.scalar.history_log_scale,
    );
    let addresses = (
        history.flow.history.as_ptr(),
        history.scalar.history.as_ptr(),
        history.scalar.workspace.basis.as_ptr(),
        history.scalar.workspace.diffusivity.as_ptr(),
    );
    let mut polls = 0;
    let trial = history
        .prepare(
            momentum_config(),
            scalar_config(),
            0.2,
            &rho,
            &u,
            &phi,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let reports = (trial.momentum_report(), trial.scalar_report());
    let viscosity = trial.viscosities().to_vec();
    let diffusivity = trial.diffusivities().to_vec();
    let scalar_coefficients = trial.scalar_coefficients().collect::<Vec<_>>();
    assert!(reports.1.maximum_coefficient > 0.0);
    drop(trial);
    for cutoff in [1, polls / 2, polls] {
        let mut calls = 0;
        assert!(matches!(
            history.prepare(
                momentum_config(),
                scalar_config(),
                0.2,
                &rho,
                &u,
                &phi,
                || {
                    calls += 1;
                    calls == cutoff
                }
            ),
            Err(DynamicSgsError::Cancelled)
        ));
        assert_eq!(history.flow.history, accepted_momentum);
        assert_eq!(history.scalar.history, accepted_scalar);
        assert_eq!(
            (
                history.flow.history_log_scale,
                history.scalar.history_log_scale
            ),
            old_scales
        );
        assert_eq!(history.accepted_time(), Some(0.0));
    }
    let denied = LagrangianScalarConfig {
        maximum_coefficient: 0.0,
        ..scalar_config()
    };
    assert!(matches!(
        history.prepare(momentum_config(), denied, 0.2, &rho, &u, &phi, || false),
        Err(DynamicSgsError::Coefficient { .. })
    ));
    assert_eq!(history.flow.history, accepted_momentum);
    assert_eq!(history.scalar.history, accepted_scalar);
    assert_eq!(history.accepted_time(), Some(0.0));
    assert!(
        history
            .prepare(
                momentum_config(),
                scalar_config(),
                0.0,
                &rho,
                &u,
                &phi,
                || false
            )
            .is_err()
    );
    assert!(
        history
            .prepare(
                momentum_config(),
                scalar_config(),
                4.0,
                &rho,
                &u,
                &phi,
                || false
            )
            .is_err()
    );
    assert_eq!(
        (
            history.flow.history.as_ptr(),
            history.scalar.history.as_ptr(),
            history.scalar.workspace.basis.as_ptr(),
            history.scalar.workspace.diffusivity.as_ptr()
        ),
        addresses
    );
    let retried = history
        .prepare(
            momentum_config(),
            scalar_config(),
            0.2,
            &rho,
            &u,
            &phi,
            || false,
        )
        .unwrap();
    assert_eq!(
        (retried.momentum_report(), retried.scalar_report()),
        reports
    );
    assert_eq!(retried.viscosities(), viscosity);
    assert_eq!(retried.diffusivities(), diffusivity);
    assert!(retried.scalar_coefficients().eq(scalar_coefficients));
    retried.commit();
    assert_eq!(history.accepted_time(), Some(0.2));
}

#[test]
fn g0_joint_byte_budget_and_scalar_control_admission() {
    let shape = [5; 3];
    let bytes = LagrangianDynamicScalar::required_bytes(125).unwrap();
    assert_eq!(bytes, 264 * 125);
    assert!(matches!(
        LagrangianDynamicScalar::new(shape, [1.0; 3], bytes - 1),
        Err(DynamicSgsError::WorkspaceBytes { required, maximum })
            if required == bytes && maximum == bytes - 1
    ));
    let mut history = LagrangianDynamicScalar::new(shape, [1.0; 3], bytes).unwrap();
    assert_eq!(history.workspace_bytes(), bytes);
    assert!(LagrangianDynamicScalar::required_bytes(usize::MAX).is_err());
    assert!(LagrangianDynamicScalar::new([2, 3, 3], [1.0; 3], usize::MAX).is_err());
    assert!(LagrangianDynamicScalar::new(shape, [0.0, 1.0, 1.0], usize::MAX).is_err());
    for value in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(
            LagrangianScalarConfig {
                memory_time: value,
                ..scalar_config()
            }
            .validate()
            .is_err()
        );
    }
    for value in [0.0, -1.0, 1.01, f64::INFINITY, f64::NAN] {
        assert!(
            LagrangianScalarConfig {
                maximum_departure_cells: value,
                ..scalar_config()
            }
            .validate()
            .is_err()
        );
    }
    assert!(
        LagrangianScalarConfig {
            usage: SgsUse::AcceptanceRender,
            ..scalar_config()
        }
        .validate()
        .is_err()
    );
    let (rho, u, mut phi) = field(shape, 1.0, 1.0);
    assert!(
        history
            .prepare(
                momentum_config(),
                scalar_config(),
                0.0,
                &rho,
                &u,
                &phi[..124],
                || false
            )
            .is_err()
    );
    phi[7] = f64::NAN;
    assert!(
        history
            .prepare(
                momentum_config(),
                scalar_config(),
                0.0,
                &rho,
                &u,
                &phi,
                || false
            )
            .is_err()
    );
    assert_eq!(history.accepted_time(), None);
}
