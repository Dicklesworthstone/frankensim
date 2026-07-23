//! G0/G3 conformance for entity-bound sensor declarations and compiled
//! observation rows.

use std::collections::BTreeSet;

use fs_assimilate::{Belief, Observation, assimilate};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey, VirtualClock};
use fs_scenario::{
    EntityDeclaration, EntityKind, EntityRef, KindExpectation, ObservationSupport, ObservationTerm,
    PlacementUncertainty, ScenarioSensor, SensorCalibration, SensorDynamics, SensorError,
    SensorKind, SensorLocation, SensorMount, SensorQuantity,
};

const TEST_STREAM: StreamKey = StreamKey {
    seed: 0x0053_454e_534f_5201,
    kernel_id: 0x5345_4e53,
    tile: 0,
    iteration: 0,
};

#[derive(Clone, Copy)]
struct FixtureEntities {
    region: fs_scenario::EntityId,
    surface: fs_scenario::EntityId,
}

fn entities() -> FixtureEntities {
    let assembly = EntityDeclaration::assembly("instrumented-rig").identity();
    let part = EntityDeclaration::part(assembly, "heated-block").identity();
    let region = EntityDeclaration::region(part, "bulk").identity();
    let surface = EntityDeclaration::surface(part, "top-face").identity();
    FixtureEntities { region, surface }
}

fn physical_sensor(
    name: &str,
    kind: SensorKind,
    target: fs_scenario::EntityId,
    expectation: KindExpectation,
    support: ObservationSupport,
    mount: SensorMount,
    placement: PlacementUncertainty,
    candidate: bool,
) -> ScenarioSensor {
    ScenarioSensor::new(
        name,
        kind,
        SensorLocation::new(
            EntityRef::new(target, expectation),
            [0.01, 0.02, 0.0],
            placement,
        )
        .expect("location"),
        support,
        mount,
        SensorDynamics::first_order(0.25, "instrument-data-sheet").expect("dynamics"),
        SensorCalibration::physical("cal-cert-2026-0042", "2026-06-30", "accredited-lab", 0.04)
            .expect("calibration"),
        candidate,
    )
    .expect("sensor")
}

#[test]
fn every_sensor_family_has_a_closed_quantity_and_entity_contract() {
    let e = entities();
    let exact = PlacementUncertainty::declared_exact("fixture-location").expect("placement");
    let point = || ObservationSupport::point(4, 2).expect("point");
    let mount = || SensorMount::declared_ideal("declared-perfect-mount").expect("mount");

    let cases = [
        (
            SensorKind::Thermocouple,
            e.region,
            KindExpectation::Domain,
            SensorQuantity::Temperature,
        ),
        (
            SensorKind::Rtd,
            e.surface,
            KindExpectation::Boundary,
            SensorQuantity::Temperature,
        ),
        (
            SensorKind::FlowMeter,
            e.surface,
            KindExpectation::Boundary,
            SensorQuantity::VolumetricFlow,
        ),
        (
            SensorKind::PressureTap,
            e.region,
            KindExpectation::Domain,
            SensorQuantity::Pressure,
        ),
        (
            SensorKind::IrCameraRegion,
            e.surface,
            KindExpectation::Exact(EntityKind::Surface),
            SensorQuantity::Temperature,
        ),
    ];

    for (kind, target, expectation, quantity) in cases {
        let sensor = physical_sensor(
            kind.label(),
            kind,
            target,
            expectation,
            point(),
            mount(),
            exact.clone(),
            false,
        );
        assert_eq!(sensor.quantity(), quantity);
        assert_eq!(
            quantity.dims(),
            match quantity {
                SensorQuantity::Temperature => [0, 0, 0, 1, 0, 0],
                SensorQuantity::VolumetricFlow => [3, 0, -1, 0, 0, 0],
                SensorQuantity::Pressure => [-1, 1, -2, 0, 0, 0],
            }
        );
    }

    let wrong = ScenarioSensor::new(
        "camera-in-volume",
        SensorKind::IrCameraRegion,
        SensorLocation::new(
            EntityRef::new(e.region, KindExpectation::Domain),
            [0.0; 3],
            exact,
        )
        .expect("location"),
        point(),
        mount(),
        SensorDynamics::declared_instantaneous("steady-only").expect("dynamics"),
        SensorCalibration::virtual_sensor("virtual-camera").expect("virtual"),
        false,
    );
    assert!(matches!(
        wrong,
        Err(SensorError::UnsupportedEntityKind {
            kind: SensorKind::IrCameraRegion,
            entity_kind: EntityKind::Region
        })
    ));
}

#[test]
fn mount_and_placement_uncertainty_compile_into_one_assimilation_row() {
    let e = entities();
    let sensor = physical_sensor(
        "tc-top",
        SensorKind::Thermocouple,
        e.surface,
        KindExpectation::Boundary,
        ObservationSupport::point(3, 1).expect("point"),
        SensorMount::affine(0.98, 1.5, "bead-contact-model").expect("mount"),
        PlacementUncertainty::axis_aligned(
            [0.001, 0.002, 0.0],
            [10.0, 5.0, 0.0],
            "placement-survey",
        )
        .expect("placement"),
        true,
    );

    let compiled = sensor.compile().expect("compile");
    assert_eq!(compiled.operator(), &[0.0, 0.98, 0.0]);
    assert_eq!(compiled.offset(), 1.5);
    assert!((compiled.placement_variance() - 0.0002).abs() < 1.0e-15);
    assert_eq!(compiled.instrument_variance(), Some(0.04));
    assert!(compiled.is_placement_candidate());
    assert!(!compiled.is_virtual());

    let state = [300.0, 310.0, 320.0];
    let predicted = compiled.predict(&state).expect("prediction");
    assert!((predicted - 305.3).abs() < 1.0e-12);

    let comparison = compiled.compare(304.5, &state).expect("comparison");
    assert_eq!(comparison.measured(), 304.5);
    assert_eq!(comparison.predicted(), predicted);
    assert!((comparison.residual() + 0.8).abs() < 1.0e-12);

    let parts = compiled.observation_parts(304.5).expect("handoff");
    assert_eq!(parts.operator(), compiled.operator());
    assert_eq!(parts.adjusted_value(), 303.0);
    assert!((parts.noise_variance() - 0.0402).abs() < 1.0e-15);
    assert_eq!(
        parts.instrument_identity(),
        compiled.sensor_identity().to_hex()
    );
    assert_eq!(parts.instrument_identity().len(), 64);

    let observation = Observation::new(
        parts.operator().to_vec(),
        parts.adjusted_value(),
        parts.noise_variance(),
        parts.instrument_identity(),
    )
    .expect("fs-assimilate observation");
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let gate = CancelGate::new();
    let clock = VirtualClock::new();
    let (prior, posterior) = pool
        .scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                TEST_STREAM,
                Budget::INFINITE,
                ExecMode::Deterministic,
            )
            .with_time_source(&clock);
            let prior = Belief::diagonal(vec![300.0; 3], &[25.0; 3], &cx)?;
            let posterior = assimilate(&prior, &observation, &cx)?;
            Ok::<_, fs_assimilate::AssimError>((prior, posterior))
        })
        .expect("assimilation");
    assert!(posterior.component_mean(1).expect("mean") > prior.component_mean(1).expect("mean"));
    assert!(posterior.variance(1).expect("variance") < prior.variance(1).expect("variance"));
    assert!(pool.stats().quiescent());

    println!(
        "{{\"suite\":\"fs-scenario/sensors\",\"case\":\"instrumented-temperature\",\"status\":\"pass\",\"raw_field\":310.0,\"predicted\":{predicted},\"measured\":{},\"residual\":{},\"noise_variance\":{},\"posterior_mean\":{},\"posterior_variance\":{},\"placement_candidate\":true,\"claim\":\"compiled-linear-observation-and-assimilation-not-calibration-authentication\"}}",
        comparison.measured(),
        comparison.residual(),
        parts.noise_variance(),
        posterior.component_mean(1).expect("mean"),
        posterior.variance(1).expect("variance"),
    );
}

#[test]
fn patch_virtual_sensor_unifies_qoi_and_comparison_but_refuses_fake_noise() {
    let e = entities();
    let support = ObservationSupport::patch_average(
        4,
        vec![ObservationTerm::new(2, 0.75), ObservationTerm::new(0, 0.25)],
    )
    .expect("patch");
    assert_eq!(
        support,
        ObservationSupport::patch_average(
            4,
            vec![ObservationTerm::new(0, 0.25), ObservationTerm::new(2, 0.75)],
        )
        .expect("canonical patch")
    );
    let sensor = ScenarioSensor::new(
        "ir-patch-qoi",
        SensorKind::IrCameraRegion,
        SensorLocation::new(
            EntityRef::new(e.surface, KindExpectation::Boundary),
            [0.03, 0.04, 0.0],
            PlacementUncertainty::declared_exact("virtual-mesh-support").expect("placement"),
        )
        .expect("location"),
        support,
        SensorMount::declared_ideal("virtual-ideal-mount").expect("mount"),
        SensorDynamics::declared_instantaneous("steady-qoi").expect("dynamics"),
        SensorCalibration::virtual_sensor("qoi-definition:ir-patch-v1").expect("virtual"),
        true,
    )
    .expect("sensor");

    let compiled = sensor.compile().expect("compile");
    assert_eq!(compiled.operator(), &[0.25, 0.0, 0.75, 0.0]);
    assert_eq!(compiled.predict(&[300.0, 0.0, 320.0, 0.0]), Ok(315.0));
    let comparison = compiled
        .compare(314.0, &[300.0, 0.0, 320.0, 0.0])
        .expect("comparison");
    assert_eq!(comparison.predicted(), 315.0);
    assert_eq!(comparison.measured(), 314.0);
    assert_eq!(comparison.residual(), -1.0);
    assert!(compiled.is_virtual());
    assert!(compiled.is_placement_candidate());
    assert_eq!(
        compiled.observation_parts(314.0),
        Err(SensorError::VirtualSensorHasNoNoiseAuthority)
    );
}

#[test]
fn malformed_support_uncertainty_mount_dynamics_and_calibration_refuse() {
    let e = entities();
    assert!(matches!(
        ObservationSupport::point(0, 0),
        Err(SensorError::InvalidStateDimension { .. })
    ));
    assert!(matches!(
        ObservationSupport::point(2, 2),
        Err(SensorError::ComponentOutOfRange { .. })
    ));
    assert!(matches!(
        ObservationSupport::patch_average(
            2,
            vec![ObservationTerm::new(0, 0.5), ObservationTerm::new(0, 0.5)]
        ),
        Err(SensorError::DuplicateComponent { component: 0 })
    ));
    assert!(matches!(
        ObservationSupport::patch_average(
            2,
            vec![ObservationTerm::new(0, 0.4), ObservationTerm::new(1, 0.4)]
        ),
        Err(SensorError::PatchWeightsDoNotSumToOne { .. })
    ));
    assert!(matches!(
        PlacementUncertainty::axis_aligned([0.0; 3], [1.0; 3], "source"),
        Err(SensorError::InvalidPlacementUncertainty { .. })
    ));
    assert!(matches!(
        SensorMount::affine(0.0, 0.0, "source"),
        Err(SensorError::InvalidMount { .. })
    ));
    assert!(matches!(
        SensorDynamics::first_order(-1.0, "source"),
        Err(SensorError::InvalidDynamics { .. })
    ));
    assert_eq!(
        SensorCalibration::physical("cert", "2025-02-29", "lab", 1.0),
        Err(SensorError::InvalidCalibrationDate)
    );
    assert_eq!(
        SensorCalibration::physical("cert", "2024-02-29", "lab", 0.0),
        Err(SensorError::InvalidInstrumentVariance)
    );
    assert!(matches!(
        SensorLocation::new(
            EntityRef::new(e.surface, KindExpectation::Domain),
            [0.0; 3],
            PlacementUncertainty::declared_exact("source").expect("placement")
        ),
        Err(SensorError::EntityExpectationMismatch { .. })
    ));
}

#[test]
fn compiled_operator_refuses_wrong_shape_and_non_finite_state() {
    let e = entities();
    let sensor = physical_sensor(
        "pressure",
        SensorKind::PressureTap,
        e.region,
        KindExpectation::Domain,
        ObservationSupport::point(2, 0).expect("point"),
        SensorMount::declared_ideal("flush-tap").expect("mount"),
        PlacementUncertainty::declared_exact("machined-port").expect("placement"),
        false,
    );
    let compiled = sensor.compile().expect("compile");
    assert!(matches!(
        compiled.predict(&[1.0]),
        Err(SensorError::StateDimensionMismatch {
            expected: 2,
            actual: 1
        })
    ));
    assert_eq!(
        compiled.predict(&[f64::NAN, 0.0]),
        Err(SensorError::NonFiniteField {
            field: "sensor state"
        })
    );
    assert_eq!(
        compiled.compare(f64::INFINITY, &[1.0, 0.0]),
        Err(SensorError::NonFiniteField {
            field: "measured reading"
        })
    );
}

#[test]
fn sensor_identity_binds_every_semantic_family() {
    let e = entities();
    let build = |name: &str,
                 mount: SensorMount,
                 dynamics: SensorDynamics,
                 calibration: SensorCalibration,
                 placement: PlacementUncertainty,
                 candidate: bool| {
        ScenarioSensor::new(
            name,
            SensorKind::Thermocouple,
            SensorLocation::new(
                EntityRef::new(e.surface, KindExpectation::Boundary),
                [0.01, 0.02, 0.0],
                placement,
            )
            .expect("location"),
            ObservationSupport::point(3, 1).expect("point"),
            mount,
            dynamics,
            calibration,
            candidate,
        )
        .expect("sensor")
        .identity()
    };
    let mount = || SensorMount::declared_ideal("mount-a").expect("mount");
    let dynamics = || SensorDynamics::declared_instantaneous("dynamics-a").expect("dynamics");
    let calibration =
        || SensorCalibration::physical("cert-a", "2026-01-02", "lab-a", 0.25).expect("calibration");
    let placement = || PlacementUncertainty::declared_exact("placement-a").expect("placement");

    let identities = [
        build(
            "sensor-a",
            mount(),
            dynamics(),
            calibration(),
            placement(),
            false,
        ),
        build(
            "sensor-b",
            mount(),
            dynamics(),
            calibration(),
            placement(),
            false,
        ),
        build(
            "sensor-a",
            SensorMount::affine(1.0, 0.1, "mount-a").expect("mount"),
            dynamics(),
            calibration(),
            placement(),
            false,
        ),
        build(
            "sensor-a",
            mount(),
            SensorDynamics::first_order(1.0, "dynamics-a").expect("dynamics"),
            calibration(),
            placement(),
            false,
        ),
        build(
            "sensor-a",
            mount(),
            dynamics(),
            SensorCalibration::physical("cert-b", "2026-01-02", "lab-a", 0.25)
                .expect("calibration"),
            placement(),
            false,
        ),
        build(
            "sensor-a",
            mount(),
            dynamics(),
            calibration(),
            PlacementUncertainty::axis_aligned([0.001, 0.0, 0.0], [1.0, 0.0, 0.0], "placement-a")
                .expect("placement"),
            false,
        ),
        build(
            "sensor-a",
            mount(),
            dynamics(),
            calibration(),
            placement(),
            true,
        ),
    ];
    assert_eq!(
        identities.into_iter().collect::<BTreeSet<_>>().len(),
        identities.len()
    );
}

#[test]
fn retained_first_order_dynamics_are_explicit_but_not_silently_applied() {
    let e = entities();
    let sensor = physical_sensor(
        "rtd-dynamic",
        SensorKind::Rtd,
        e.region,
        KindExpectation::Domain,
        ObservationSupport::point(1, 0).expect("point"),
        SensorMount::declared_ideal("embedded-rtd").expect("mount"),
        PlacementUncertainty::declared_exact("molded-location").expect("placement"),
        false,
    );
    assert_eq!(sensor.dynamics().time_constant_s(), Some(0.25));
    assert_eq!(
        sensor.compile().expect("compile").predict(&[350.0]),
        Ok(350.0)
    );
}
