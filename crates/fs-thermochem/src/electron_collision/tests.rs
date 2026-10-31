use super::*;

fn id(s: &str) -> SpeciesId {
    SpeciesId::new(s).unwrap()
}
fn controls() -> MaxwellianRateControlsV1 {
    MaxwellianRateControlsV1 {
        relative_tolerance: 1e-9,
        maximum_subdivisions: 4096,
    }
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs(),
        "{actual:.17e} != {expected:.17e}"
    );
}
fn momentum_table(upper_ev: f64, first: f64, last: f64) -> ElectronCrossSectionV1 {
    ElectronCrossSectionV1::tabulated(
        id("neutral"),
        ElectronCollisionKindV1::ElasticMomentumTransfer,
        Provenance {
            source: "manufactured analytic collision fixture, not measured data".into(),
            license: "fixture".into(),
            artifact: None,
        },
        vec![
            ElectronCrossSectionPointV1 {
                energy: Energy::new(0.0),
                cross_section: Area::new(first),
            },
            ElectronCrossSectionPointV1 {
                energy: Energy::new(upper_ev * ELECTRONVOLT_J),
                cross_section: Area::new(last),
            },
        ],
        Area::new(last),
    )
    .unwrap()
}

#[test]
fn generalized_beb_retains_published_parameters_and_thresholds() {
    // Independent 50-decimal evaluation of the printed NIST Eq.1 with its
    // published B,U,N,Q and constants. These are theoretical-formula checks,
    // not experimental validation or a claim to reproduce the BED model.
    for (source, energy, expected) in [
        (nist::hydrogen(id("H")), 20.0, 2.974_727_565_291_051_6e-21),
        (nist::hydrogen(id("H")), 50.0, 6.112_468_194_002_238e-21),
        (nist::hydrogen(id("H")), 1000.0, 1.187_921_933_687_247_7e-21),
        (nist::helium(id("He")), 50.0, 2.526_268_151_265_614_7e-21),
        (nist::helium(id("He")), 100.0, 3.601_288_815_626_636e-21),
    ] {
        near(
            source
                .cross_section(Energy::new(energy * ELECTRONVOLT_J))
                .unwrap()
                .value(),
            expected,
            8e-14,
        );
        assert_eq!(
            source.cross_section(source.threshold()).unwrap().value(),
            0.0
        );
        assert_eq!(source.cross_section(Energy::new(0.0)).unwrap().value(), 0.0);
        assert_eq!(
            source.cross_section(Energy::new(10_001.0 * ELECTRONVOLT_J)),
            Err(ElectronCollisionErrorV1::OutsideEnergyDomain)
        );
    }
    let h = nist::hydrogen(id("original-H"));
    assert_eq!(h.threshold().value(), 13.6057 * ELECTRONVOLT_J);
    assert!(h.provenance().source.contains("SRD 107"));
    let just_above = Energy::new(f64::from_bits(h.threshold().value().to_bits() + 1));
    assert!(h.cross_section(just_above).unwrap().value() > 0.0);
}

#[test]
fn nist_maxwellian_rate_and_derivative_match_independent_integral() {
    // Independent composite Simpson integration in incident energy, 131072
    // uniform panels over 128 thermal energies above threshold. The values
    // below integrate the generalized theoretical BEB formula, not a fit.
    for (source, temperature, expected, derivative) in [
        (
            nist::hydrogen(id("H")),
            10_000.0,
            7.574_964_457_176_377e-22,
            1.236_213_866_094_943e-24,
        ),
        (
            nist::hydrogen(id("H")),
            100_000.0,
            3.676_538_059_223_906e-15,
            7.440_481_923_236_229e-20,
        ),
        (
            nist::helium(id("He")),
            10_000.0,
            9.087_366_353_639_62e-28,
            2.644_172_189_570_746_2e-30,
        ),
        (
            nist::helium(id("He")),
            100_000.0,
            5.190_629_180_569_401e-16,
            1.791_743_736_775_225_8e-20,
        ),
    ] {
        let evaluation = source
            .evaluate(
                Temperature::new(temperature),
                MaxwellianRateWeightingV1::Event,
                controls(),
                || false,
            )
            .unwrap();
        near(evaluation.coefficient_m3_s().unwrap(), expected, 2e-11);
        near(evaluation.derivative_m3_s_k().unwrap(), derivative, 3e-11);
        assert_eq!(evaluation.source(), &source);
        assert_eq!(
            evaluation.electron_temperature(),
            Temperature::new(temperature)
        );
        assert!(evaluation.properties().relative_tail_bound() < controls().relative_tolerance);
        assert!(
            evaluation
                .properties()
                .estimated_relative_quadrature_error()
                < controls().relative_tolerance
        );
    }
}

#[test]
fn drag_has_correct_maxwellian_weight_not_event_count() {
    let t = 12_000.0;
    let cross = 3e-20;
    let accurate = MaxwellianRateControlsV1 {
        relative_tolerance: 1e-12,
        ..controls()
    };
    let source = momentum_table(200.0, cross, cross);
    let p = source
        .evaluate_properties(
            Temperature::new(t),
            MaxwellianRateWeightingV1::WeakDriftMomentum,
            accurate,
            || false,
        )
        .unwrap();
    let mean_speed =
        (8.0 * BOLTZMANN_J_PER_K * t / (core::f64::consts::PI * ELECTRON_MASS_KG)).sqrt();
    let expected = (4.0 / 3.0) * cross * mean_speed;
    near(p.coefficient_m3_s().unwrap(), expected, 2e-12);
    near(p.derivative_m3_s_k().unwrap(), expected / (2.0 * t), 2e-12);
    // A linear energy dependence changes both the weighted average and slope.
    let linear = momentum_table(200.0, 0.0, 200.0 * cross);
    let p = linear
        .evaluate_properties(
            Temperature::new(t),
            MaxwellianRateWeightingV1::WeakDriftMomentum,
            accurate,
            || false,
        )
        .unwrap();
    let expected_linear = 4.0 * cross * (BOLTZMANN_J_PER_K * t / ELECTRONVOLT_J) * mean_speed;
    near(p.coefficient_m3_s().unwrap(), expected_linear, 2e-12);
    near(
        p.derivative_m3_s_k().unwrap(),
        1.5 * expected_linear / t,
        2e-12,
    );
}

#[test]
fn source_coverage_weighting_work_and_cancellation_fail_closed() {
    let source = momentum_table(1.0, 1e-20, 1e-20);
    assert!(matches!(
        source.evaluate_properties(
            Temperature::new(100_000.0),
            MaxwellianRateWeightingV1::WeakDriftMomentum,
            controls(),
            || false
        ),
        Err(ElectronCollisionErrorV1::UnresolvedTail { .. })
    ));
    assert!(
        source
            .evaluate_properties(
                Temperature::new(50.0),
                MaxwellianRateWeightingV1::WeakDriftMomentum,
                controls(),
                || false
            )
            .is_ok()
    );
    assert!(matches!(
        source.evaluate_properties(
            Temperature::new(300.0),
            MaxwellianRateWeightingV1::Event,
            controls(),
            || false
        ),
        Err(ElectronCollisionErrorV1::Invalid(_))
    ));
    let h = nist::hydrogen(id("H"));
    assert!(matches!(
        h.evaluate_properties(
            Temperature::new(300.0),
            MaxwellianRateWeightingV1::WeakDriftMomentum,
            controls(),
            || false
        ),
        Err(ElectronCollisionErrorV1::Invalid(_))
    ));
    assert_eq!(
        h.evaluate_properties(
            Temperature::new(10000.0),
            MaxwellianRateWeightingV1::Event,
            controls(),
            || true
        ),
        Err(ElectronCollisionErrorV1::Cancelled)
    );
    let tiny = MaxwellianRateControlsV1 {
        maximum_subdivisions: 1,
        ..controls()
    };
    assert_eq!(
        h.evaluate_properties(
            Temperature::new(10000.0),
            MaxwellianRateWeightingV1::Event,
            tiny,
            || false
        ),
        Err(ElectronCollisionErrorV1::QuadratureBudget)
    );
    assert!(matches!(
        h.evaluate_properties(
            Temperature::new(1e9),
            MaxwellianRateWeightingV1::Event,
            controls(),
            || false
        ),
        Err(ElectronCollisionErrorV1::UnresolvedTail { .. })
    ));
}

#[test]
fn rare_ionization_keeps_finite_log_rate_when_direct_value_underflows() {
    let p = nist::hydrogen(id("H"))
        .evaluate_properties(
            Temperature::new(100.0),
            MaxwellianRateWeightingV1::Event,
            controls(),
            || false,
        )
        .unwrap();
    assert_eq!(p.coefficient_m3_s(), None);
    assert!(p.log_coefficient_m3_s().is_finite());
    assert!(p.logarithmic_temperature_derivative_per_k() > 0.0);
}

#[test]
fn every_table_knot_is_integrated_even_for_a_narrow_collision_band() {
    let source = ElectronCrossSectionV1::tabulated(
        id("H"),
        ElectronCollisionKindV1::Ionization,
        Provenance {
            source: "synthetic triangular event cross section".into(),
            license: "fixture".into(),
            artifact: None,
        },
        vec![(1.0, 0.0), (1.0001, 1e-20), (1.0002, 0.0)]
            .into_iter()
            .map(|(e, s)| ElectronCrossSectionPointV1 {
                energy: Energy::new(e * ELECTRONVOLT_J),
                cross_section: Area::new(s),
            })
            .collect(),
        Area::new(0.0),
    )
    .unwrap();
    let t = 1.0 / (BOLTZMANN_J_PER_K / ELECTRONVOLT_J);
    let p = source
        .evaluate_properties(
            Temperature::new(t),
            MaxwellianRateWeightingV1::Event,
            controls(),
            || false,
        )
        .unwrap();
    assert!(p.coefficient_m3_s().unwrap() > 0.0);
    assert_eq!(p.relative_tail_bound(), 0.0);
    let step = 0.1;
    let plus = source
        .evaluate_properties(
            Temperature::new(t + step),
            MaxwellianRateWeightingV1::Event,
            controls(),
            || false,
        )
        .unwrap()
        .coefficient_m3_s()
        .unwrap();
    let minus = source
        .evaluate_properties(
            Temperature::new(t - step),
            MaxwellianRateWeightingV1::Event,
            controls(),
            || false,
        )
        .unwrap()
        .coefficient_m3_s()
        .unwrap();
    near(
        p.derivative_m3_s_k().unwrap(),
        (plus - minus) / (2.0 * step),
        2e-6,
    );
}
