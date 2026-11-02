use super::*;

fn provenance() -> Provenance {
    Provenance {
        source: "authored ideal-tank fixture; no measured plant calibration".into(),
        license: "project fixture".into(),
        artifact: None,
    }
}
fn train(n: usize, volume: f64, concentration: f64) -> StirredTankCascade {
    StirredTankCascade::uniform(
        n,
        Volume::new(volume),
        Density::new(concentration),
        provenance(),
    )
    .unwrap()
}
fn advance(train: &mut StirredTankCascade, dt: f64, q: f64, feed: f64) -> TankInterval {
    train
        .advance(
            Time::new(dt),
            VolumetricFlowRate::new(q),
            Density::new(feed),
        )
        .unwrap()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a} != {b}");
}

#[test]
fn single_tank_endpoint_and_independent_mean_known_answers() {
    let mut tanks = train(1, 2.0, 1.0);
    let interval = advance(&mut tanks, 4.0, 0.5, 5.0);
    close(
        tanks.tank_concentrations()[0].value(),
        3.5284822353142307,
        2e-15,
    );
    close(
        interval.outlet_mean.unwrap().value(),
        2.4715177646857693,
        2e-15,
    );
    close(interval.outlet_mass.value(), 4.943035529371539, 4e-15);
    assert_eq!(interval.inlet_mass.value(), 10.0);
    close(interval.final_inventory.value(), 7.056964470628461, 4e-15);
    close(interval.conservation_residual.value(), 0.0, 8e-15);
    assert_eq!(tanks.inventory(), interval.final_inventory);
    assert_eq!(tanks.provenance(), &provenance());
}

#[test]
fn small_exposure_keeps_positive_endpoint_and_integrated_outlet() {
    let mut tanks = train(1, 1.0, 0.0);
    let x = 1e-8;
    let interval = advance(&mut tanks, x, 1.0, 1.0);
    close(
        tanks.tank_concentrations()[0].value(),
        -(-x).exp_m1(),
        2e-23,
    );
    let independent_mean = x / 2.0 - x * x / 6.0 + x * x * x / 24.0;
    close(
        interval.outlet_mean.unwrap().value(),
        independent_mean,
        2e-23,
    );
    close(interval.outlet_mass.value(), x * independent_mean, 2e-31);
    assert!(interval.outlet_mass.value() > 0.0);
}

#[test]
fn source_scaling_preserves_responses_when_bare_probabilities_underflow() {
    let mut single = train(1, 1.0, 0.0);
    let result = advance(&mut single, 1e-200, 1.0, 1e200);
    close(single.tank_concentrations()[0].value(), 1.0, 2e-16);
    close(result.outlet_mean.unwrap().value(), 0.5, 1e-16);
    close(result.outlet_mass.value(), 5e-201, 2e-216);

    // Same tiny exposure through two tanks, but enough active volume for the
    // integrated outlet mass itself to be representable as well as C2.
    let mut two = train(2, 2e200, 0.0);
    let result = advance(&mut two, 1.0, 1.0, 1e200);
    close(two.tank_concentrations()[1].value(), 5e-201, 2e-216);
    close(result.outlet_mean.unwrap().value(), 1e-200 / 6.0, 1e-216);
    assert!(result.outlet_mass.value() > 0.0);
}

#[test]
fn truly_unrepresentable_positive_paths_refuse_without_publishing_zero() {
    let mut two = train(2, 2.0, 0.0);
    let before = two.clone();
    for feed in [1.0, 1e200] {
        assert!(matches!(
            two.advance(
                Time::new(1e-200),
                VolumetricFlowRate::new(1.0),
                Density::new(feed)
            ),
            Err(ResidenceError::Arithmetic(_))
        ));
        assert_eq!(two, before);
    }
}

#[test]
fn stopped_flow_preserves_inventory_and_variable_flow_uses_actual_volume() {
    let mut varied = StirredTankCascade::new(
        Volume::new(50.0),
        [1.0, 0.0, 3.0, 0.5, 4.0].map(Density::new).to_vec(),
        provenance(),
    )
    .unwrap();
    let mut reference = varied.clone();
    let first = advance(&mut varied, 10.0, 0.5, 2.0);
    let stopped_state = varied.clone();
    let stop = advance(&mut varied, 100_000.0, 0.0, 99.0);
    assert_eq!(varied, stopped_state);
    assert_eq!(stop.outlet_mean, None);
    assert_eq!(stop.inlet_mass.value(), 0.0);
    assert_eq!(stop.outlet_mass.value(), 0.0);
    assert_eq!(stop.initial_inventory, stop.final_inventory);
    let zero_time = advance(&mut varied, 0.0, 100.0, 99.0);
    assert_eq!(zero_time.outlet_mean, None);
    assert_eq!(varied, stopped_state);
    let second = advance(&mut varied, 5.0, 2.0, 2.0);
    let whole = advance(&mut reference, 15.0, 1.0, 2.0);
    for (a, b) in varied
        .tank_concentrations()
        .iter()
        .zip(reference.tank_concentrations())
    {
        close(a.value(), b.value(), 3e-15);
    }
    close(
        first.outlet_mass.value() + second.outlet_mass.value(),
        whole.outlet_mass.value(),
        2e-14,
    );
}

// Independent Erlang CDF from its finite survival polynomial, with the standard
// library exp. Production uses positive infinite Poisson tails and det::exp.
fn erlang_step(hours: f64) -> f64 {
    if hours <= 0.0 {
        return 0.0;
    }
    let x = hours / 3.0;
    let (mut term, mut sum) = (1.0, 1.0);
    for j in 1..8 {
        term *= x / f64::from(j);
        sum += term;
    }
    1.0 - (-x).exp() * sum
}

#[test]
fn eight_tank_pulse_exposes_shift_recovery_bias_without_fitting() {
    // Synthetic GB41-shaped illustration, not observations from a named mine.
    // 1000 kg dry solids/m³ gives density[kg/m³] = grade[g/t] * 0.001.
    // A fixed explicit 88% extraction leaves the passive tail constituent.
    let mut tanks = train(8, 86_400.0, 0.36e-3);
    let references = [
        0.36605285770815293,
        0.5488988905925538,
        0.7139539382608006,
        0.5040462968338287,
    ];
    let mut total_in = 0.0;
    let mut total_out = 0.0;
    let initial_inventory = tanks.inventory().value();
    for (shift, expected_grade) in references.into_iter().enumerate() {
        let head_grade = if shift == 0 { 9.0 } else { 3.0 };
        let interval = advance(&mut tanks, 43_200.0, 1.0, head_grade * 0.12e-3);
        let grade = interval.outlet_mean.unwrap().value() * 1000.0;
        close(grade, expected_grade, 3e-14);
        // Integrate the analytic pulse response independently, without using
        // initial/final inventories to construct an outlet balance plug.
        let n = 4096;
        let mut sum = 0.0;
        for i in 0..=n {
            let t = 12.0 * shift as f64 + 12.0 * f64::from(i) / f64::from(n);
            let outlet_grade = 0.36 + 0.72 * (erlang_step(t) - erlang_step(t - 12.0));
            let weight = if i == 0 || i == n {
                1.0
            } else if i % 2 == 0 {
                2.0
            } else {
                4.0
            };
            sum += weight * outlet_grade;
        }
        close(grade, sum / (3.0 * f64::from(n)), 2e-13);
        let apparent = 1.0 - grade / head_grade;
        assert!((apparent - 0.88).abs() > 0.04);
        total_in += interval.inlet_mass.value();
        total_out += interval.outlet_mass.value();
    }
    close(
        tanks.inventory().value() - initial_inventory + total_out - total_in,
        0.0,
        3e-13,
    );
    close(erlang_step(12.0), 0.051133615792847475, 2e-15);
    close(
        erlang_step(36.0) - erlang_step(12.0),
        0.8593618873669767,
        2e-15,
    );
}

#[test]
fn maximum_train_retains_positive_mass_and_long_interval_matches_split() {
    let mut whole = train(MAX_TANKS, 64.0, 0.25);
    let mut split = whole.clone();
    let result = advance(&mut whole, 80.0, 1.0, 1.0);
    let mut split_out = 0.0;
    for _ in 0..20 {
        split_out += advance(&mut split, 4.0, 1.0, 1.0).outlet_mass.value();
    }
    assert_eq!(result.substeps, 10);
    for (a, b) in whole
        .tank_concentrations()
        .iter()
        .zip(split.tank_concentrations())
    {
        assert!((0.25..=1.0 + 2e-15).contains(&a.value()));
        close(a.value(), b.value(), 2e-14);
    }
    close(result.outlet_mass.value(), split_out, 5e-13);
    assert!(result.conservation_residual.value().abs() < 2e-12);
}

#[test]
fn refusals_leave_existing_physical_state_unchanged() {
    let mut tanks = train(8, 24.0, 0.3);
    let before = tanks.clone();
    for (dt, q, feed) in [
        (1.0, -1.0, 1.0),
        (-1.0, 1.0, 1.0),
        (1.0, 1.0, f64::NAN),
        (1.0, f64::INFINITY, 1.0),
        (1e8, 1.0, 1.0),
        (f64::MAX, 2.0, 1.0),
        (1.0, 2.0, f64::MAX),
        (f64::MIN_POSITIVE, f64::MIN_POSITIVE, 1.0),
    ] {
        assert!(
            tanks
                .advance(
                    Time::new(dt),
                    VolumetricFlowRate::new(q),
                    Density::new(feed)
                )
                .is_err()
        );
        assert_eq!(tanks, before);
    }
    for (n, v, c) in [
        (0, 1.0, 1.0),
        (65, 1.0, 1.0),
        (1, 0.0, 1.0),
        (1, f64::NAN, 1.0),
        (1, 1.0, -1.0),
        (1, 1.0, f64::INFINITY),
        (2, f64::MAX, 2.0),
    ] {
        assert!(
            StirredTankCascade::uniform(n, Volume::new(v), Density::new(c), provenance()).is_err()
        );
    }
    let mut missing = provenance();
    missing.source.clear();
    assert!(StirredTankCascade::uniform(1, Volume::new(1.0), Density::new(0.0), missing).is_err());
}
