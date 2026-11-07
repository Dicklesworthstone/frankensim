use super::*;

fn config() -> LagrangianSgsConfig {
    LagrangianSgsConfig {
        memory_time: 2.0,
        maximum_departure_cells: 1.0,
        maximum_coefficient_squared: 100.0,
        maximum_eddy_to_molecular: 1000.0,
        usage: SgsUse::Exploration,
    }
}

fn field(n: usize, amplitude: f64) -> (Vec<f64>, Vec<[f64; 3]>) {
    let mut density = Vec::new();
    let mut velocity = Vec::new();
    for z in 0..n {
        for y in 0..n {
            for x in 0..n {
                let [x, y, z] =
                    [x, y, z].map(|i| core::f64::consts::TAU * (i as f64 + 0.5) / n as f64);
                density.push(1.2 + 0.2 * x.sin());
                velocity.push([
                    amplitude * (0.3 * x.sin() + 0.2 * (y + z).sin()),
                    amplitude * (0.2 * y.cos() + 0.3 * (z + x).sin()),
                    amplitude * (0.2 * z.sin() + 0.3 * (x + y).cos()),
                ]);
            }
        }
    }
    (density, velocity)
}

#[test]
fn initial_fit_is_the_existing_local_box_model() {
    let (rho, u) = field(5, 1.0);
    let mut local = PeriodicDynamicSgs::new([5; 3], 1.0, usize::MAX).unwrap();
    let expected = local
        .prepare(
            DynamicSgsConfig {
                averaging: DynamicAveraging::LocalBox,
                maximum_coefficient_squared: config().maximum_coefficient_squared,
                maximum_eddy_to_molecular: config().maximum_eddy_to_molecular,
                usage: SgsUse::Exploration,
            },
            &rho,
            &u,
            || false,
        )
        .unwrap();
    let mut flow = LagrangianDynamicSgs::new([5; 3], [1.0; 3], usize::MAX).unwrap();
    let trial = flow.prepare(config(), 0.0, &rho, &u, || false).unwrap();
    assert_eq!(trial.report(), expected);
    assert_eq!(trial.viscosities(), local.viscosities());
    assert!(expected.maximum_eddy_viscosity > 0.0);
    trial.commit();
    assert_eq!(flow.accepted_time(), Some(0.0));
}

#[test]
fn physical_advection_translates_history_including_uniform_flow() {
    let n = 5;
    let mut flow = LagrangianDynamicSgs::new([n; 3], [1.0; 3], usize::MAX).unwrap();
    // Synthetic accepted Germano history supplies a direction-sensitive oracle.
    // Every denominator is positive and a rigid translation has no new strain.
    for cell in 0..n.pow(3) {
        let [x, y, z] = coordinates([n; 3], cell);
        flow.history[cell] = [0.01 * (1 + x + 2 * y + 3 * z) as f64, 1.0];
    }
    let original = flow.history.clone();
    flow.accepted_time = Some(0.0);
    flow.history_log_scale = Some(0.0);
    let rho = vec![1.0; n.pow(3)];
    let u = vec![[1.0, -1.0, 0.0]; n.pow(3)];
    let trial = flow.prepare(config(), 1.0, &rho, &u, || false).unwrap();
    assert!(trial.viscosities().iter().all(|v| *v == 0.0));
    for cell in 0..n.pow(3) {
        let [x, y, z] = coordinates([n; 3], cell);
        let upstream = index([n; 3], [(x + n - 1) % n, (y + 1) % n, z]);
        let pair = trial.owner.candidate[cell];
        assert!((pair[0] / pair[1] - original[upstream][0]).abs() < 1e-15);
    }
    let physical_amplitude = det::exp(trial.log_scale.unwrap());
    assert!((physical_amplitude - (-0.5_f64).exp()).abs() < 1e-14);
    trial.commit();
    assert_eq!(flow.accepted_time(), Some(1.0));
}

#[test]
fn trilinear_departures_reproduce_affine_history_and_wrap() {
    let shape = [5; 3];
    let values: Vec<_> = (0..125)
        .map(|cell| {
            let [x, y, z] = coordinates(shape, cell).map(|i| i as f64);
            [x + 2.0 * y - 3.0 * z, 1.0 + x + y + z]
        })
        .collect();
    let value = sample_history(
        FilterGeometry::Periodic,
        shape,
        index(shape, [2, 2, 2]),
        [0.25, -0.5, 0.75],
        &values,
    )
    .unwrap();
    assert_eq!(
        value,
        [1.75 + 2.0 * 2.5 - 3.0 * 1.25, 1.0 + 1.75 + 2.5 + 1.25]
    );
    let wrapped = sample_history(
        FilterGeometry::Periodic,
        shape,
        index(shape, [0, 0, 0]),
        [0.5, 0.0, 0.0],
        &values,
    )
    .unwrap();
    assert_eq!(wrapped, [2.0, 3.0]);
}

#[test]
fn long_relaxation_retains_representable_history_in_quiescent_flow() {
    let mut flow = LagrangianDynamicSgs::new([3; 3], [1.0; 3], usize::MAX).unwrap();
    flow.history.fill([0.02, 1.0]);
    flow.accepted_time = Some(0.0);
    flow.history_log_scale = Some(0.0);
    let rho = vec![1.0; 27];
    let u = vec![[0.0; 3]; 27];
    let trial = flow.prepare(config(), 80.0, &rho, &u, || false).unwrap();
    assert!(trial.viscosities().iter().all(|v| *v == 0.0));
    assert!((trial.log_scale.unwrap() + 40.0).abs() < 1e-13);
    for pair in &trial.owner.candidate {
        assert!((pair[0] / pair[1] - 0.02).abs() < 1e-16);
        assert!(pair[1] > 0.0);
    }
    trial.commit();
    let (density, mut velocity) = field(3, 1e-5);
    for vector in &mut velocity {
        *vector = vector.map(|v| v * 0.1);
    }
    let continued = flow
        .prepare(config(), 80.1, &density, &velocity, || false)
        .unwrap();
    assert!(continued.report().maximum_eddy_viscosity > 0.0);
    assert!(continued.owner.candidate.iter().all(|p| p[1] > 0.0));
}

#[test]
fn frozen_source_relaxation_uses_physical_contractions_not_ratios() {
    let (rho, u) = field(5, 0.5);
    let mut flow = LagrangianDynamicSgs::new([5; 3], [1.0; 3], usize::MAX).unwrap();
    let (scales, _) = flow
        .flow
        .prepare_contractions(&rho, &u, [[0.0; 3]; 2], || false)
        .unwrap();
    let physical_scale = scales.density.powi(2) * scales.velocity.powi(4);
    let sources: Vec<_> = flow
        .flow
        .contraction
        .iter()
        .map(|p| p.map(|v| v * physical_scale))
        .collect();
    let old = [0.02, 2.0];
    flow.history.fill(old);
    flow.history_log_scale = Some(0.0);
    flow.accepted_time = Some(0.0);
    let dt = 0.4;
    let keep = (-dt / config().memory_time).exp();
    let trial = flow.prepare(config(), dt, &rho, &u, || false).unwrap();
    let scale = det::exp(trial.log_scale.unwrap());
    for (cell, source) in sources.iter().enumerate() {
        for k in 0..2 {
            let expected = keep * old[k] + (1.0 - keep) * source[k];
            let actual = trial.owner.candidate[cell][k] * scale;
            assert!(
                (actual - expected).abs() < 2e-14,
                "cell {cell} component {k}: {actual} != {expected}"
            );
        }
    }
}

#[test]
fn changing_density_and_velocity_normalizations_preserve_units() {
    let mut reference = LagrangianDynamicSgs::new([5; 3], [1.0; 3], usize::MAX).unwrap();
    let fields = [field(5, 1.0), field(5, 0.4), field(5, 0.7)];
    let mut outputs = Vec::new();
    for (step, (rho, u)) in fields.iter().enumerate() {
        let trial = reference
            .prepare(config(), step as f64 * 0.25, rho, u, || false)
            .unwrap();
        outputs.push(trial.viscosities().to_vec());
        trial.commit();
    }
    for length in [1e-50, 1e50] {
        let mut scaled = LagrangianDynamicSgs::new([5; 3], [length; 3], usize::MAX).unwrap();
        for (step, (rho, u)) in fields.iter().enumerate() {
            let density: Vec<_> = rho.iter().map(|v| v * 1e50).collect();
            let velocity: Vec<_> = u.iter().map(|v| v.map(|c| c * length)).collect();
            let trial = scaled
                .prepare(config(), step as f64 * 0.25, &density, &velocity, || false)
                .unwrap();
            for (&actual, &expected) in trial.viscosities().iter().zip(&outputs[step]) {
                let actual = actual / (length * length);
                assert!(
                    (actual - expected).abs() < 2e-11 * expected.abs().max(1e-8),
                    "step {step}: {actual} != {expected}"
                );
            }
            trial.commit();
        }
    }
}

#[test]
fn refused_cancelled_and_dropped_trials_do_not_advance_memory() {
    let (rho, u) = field(5, 1.0);
    let mut flow = LagrangianDynamicSgs::new([5; 3], [1.0; 3], usize::MAX).unwrap();
    flow.prepare(config(), 0.0, &rho, &u, || false)
        .unwrap()
        .commit();
    let saved = flow.history.clone();
    let saved_scale = flow.history_log_scale;
    let mut poll_count = 0;
    assert!(matches!(
        flow.prepare(config(), 0.5, &rho, &u, || {
            poll_count += 1;
            poll_count == 8
        }),
        Err(DynamicSgsError::Cancelled)
    ));
    assert_eq!(flow.history, saved);
    assert_eq!(flow.history_log_scale, saved_scale);
    assert_eq!(flow.accepted_time(), Some(0.0));
    assert!(flow.prepare(config(), 0.0, &rho, &u, || false).is_err());
    assert!(flow.prepare(config(), 10.0, &rho, &u, || false).is_err());
    let report;
    let viscosity;
    {
        let trial = flow.prepare(config(), 0.5, &rho, &u, || false).unwrap();
        report = trial.report();
        viscosity = trial.viscosities().to_vec();
    }
    assert_eq!(flow.history, saved);
    assert_eq!(flow.accepted_time(), Some(0.0));
    let retry = flow.prepare(config(), 0.5, &rho, &u, || false).unwrap();
    assert_eq!(retry.report(), report);
    assert_eq!(retry.viscosities(), viscosity);
    retry.commit();
    assert_eq!(flow.accepted_time(), Some(0.5));
}

#[test]
fn byte_model_and_admission_match_real_payload() {
    let bytes = LagrangianDynamicSgs::required_bytes(125).unwrap();
    assert!(matches!(
        LagrangianDynamicSgs::new([5; 3], [1.0; 3], bytes - 1),
        Err(DynamicSgsError::WorkspaceBytes { .. })
    ));
    let mut flow = LagrangianDynamicSgs::new([5; 3], [1.0; 3], bytes).unwrap();
    assert_eq!(flow.workspace_bytes(), bytes);
    for memory_time in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            LagrangianSgsConfig {
                memory_time,
                ..config()
            }
            .validate()
            .is_err()
        );
    }
    assert!(
        LagrangianSgsConfig {
            usage: SgsUse::AcceptanceRender,
            ..config()
        }
        .validate()
        .is_err()
    );
    let (rho, u) = field(5, 1.0);
    assert!(
        flow.prepare(
            LagrangianSgsConfig {
                maximum_coefficient_squared: 0.0,
                ..config()
            },
            0.0,
            &rho,
            &u,
            || false
        )
        .is_err()
    );
    assert_eq!(flow.accepted_time(), None);
    flow.prepare(config(), 0.0, &rho, &u, || false)
        .unwrap()
        .commit();
}
