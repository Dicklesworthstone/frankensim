//! G1/G3 joint normal/sliding/thermal references; sources are synthetic fixtures.
use super::*;
use fs_contact::interface_binding::bind_dry_interface_system_card;
use fs_contact::interface_binding::thermal_friction::elliptic::SourceThermalEllipticCoulomb;
use fs_contact::interface_binding::thermal_friction::{
    FrictionQueryAxis, FrictionQueryAxisKind, SourceThermalCoulomb, ThermalFrictionAxes,
};
use fs_couple::tribothermal::normal::directional::*;
use fs_evidence::ValidityDomain;
use fs_matdb::*;
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
};
use fs_material::thermal::interface::ThermalInterfacePoint;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, Pressure, QuantitySpec, Temperature};
use fs_tribo::TangentialSolveConfig;
use std::cell::Cell;

const NAMES: [&str; 4] = ["static-x", "static-y", "kinetic-x", "kinetic-y"];
const SPEED: Dims = Dims([1, 0, -1, 0, 0, 0]);
const LENGTH: Dims = Dims([1, 0, 0, 0, 0, 0]);
fn quantity(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn context() -> QueryPoint {
    QueryPoint::new()
        .with_quantity("Ta", quantity(QuantityKind::AbsoluteTemperature), 300.0)
        .unwrap()
        .with_quantity("Tb", quantity(QuantityKind::AbsoluteTemperature), 300.0)
        .unwrap()
        .with_quantity("p", quantity(QuantityKind::Pressure), 10.0)
        .unwrap()
        .with_quantity("speed", QuantitySpec::dimensional(SPEED), 0.0)
        .unwrap()
        .with_quantity("path", QuantitySpec::dimensional(LENGTH), 0.0)
        .unwrap()
}
fn fixed(mu: f64) -> PropertyValue {
    PropertyValue::Scalar {
        value: mu,
        dims: Dims::NONE,
    }
}
fn curve(axis: &str, dims: Dims, knots: [(f64, f64); 2]) -> PropertyValue {
    PropertyValue::Curve {
        abscissa: axis.into(),
        abscissa_dims: dims,
        knots: knots.to_vec(),
        dims: Dims::NONE,
    }
}
fn interface(
    cards: &[MaterialCard; 2],
    x: PropertyValue,
    pressure_limit: f64,
) -> InterfaceSystemCard {
    let mut claims = ClaimSet::new();
    for (name, value) in NAMES
        .into_iter()
        .zip([fixed(2.0), fixed(2.0), x, fixed(0.4)])
    {
        let interpolation = if matches!(value, PropertyValue::Scalar { .. }) {
            InterpolationPolicy::ConstantWithinValidity
        } else {
            InterpolationPolicy::LinearInside
        };
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, Dims::NONE),
                value,
                interpolation,
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "Ta",
                        quantity(QuantityKind::AbsoluteTemperature),
                        200.0,
                        600.0,
                    )
                    .with_quantity(
                        "Tb",
                        quantity(QuantityKind::AbsoluteTemperature),
                        200.0,
                        600.0,
                    )
                    .with_quantity("p", quantity(QuantityKind::Pressure), 0.0, pressure_limit)
                    .with_quantity("speed", QuantitySpec::dimensional(SPEED), 0.0, 10.0)
                    .with_quantity("path", QuantitySpec::dimensional(LENGTH), 0.0, 10.0),
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "authored normal/directional coupling fixture".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    InterfaceSystemCard::assemble(
        SurfaceSpec {
            material: cards[0].id().clone(),
            texture_frame: "A direction".into(),
        },
        SurfaceSpec {
            material: cards[1].id().clone(),
            texture_frame: "B direction".into(),
        },
        SystemContext {
            medium: "dry".into(),
            third_body: None,
            environment: "fixture".into(),
            history: "virgin".into(),
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn tangent() -> TangentialSolveConfig {
    TangentialSolveConfig {
        max_iterations: 128,
        maximum_force_n: 1e8,
        force_absolute_tolerance_n: 1e-12,
        force_relative_tolerance: 1e-13,
    }
}
#[allow(clippy::too_many_arguments)]
fn joint<'a>(
    thermal: &DepositionSimulation<'_>,
    cards: &'a [MaterialCard; 2],
    curves: &'a [ResolvedHeatCapacityCurve; 2],
    card: &'a InterfaceSystemCard,
    point: &'a QueryPoint,
    indent: f64,
    damping: f64,
    slopes: [f64; 2],
    partition: HeatPartition,
    axis: [f64; 2],
) -> NormalDirectionalContact<'a> {
    let bound = bind_dry_interface_system_card(card, InputAuthority::SyntheticFixture).unwrap();
    let mut normal_request = request(indent, damping);
    normal_request.interface = bound.interface().clone();
    let normal_model =
        AffineNormalModulus::new(normal_request, [300.0; 2], slopes, [[200.0, 600.0]; 2]).unwrap();
    let normal = bind(thermal, curves, normal_model);
    let sources = [0, 1].map(|i| {
        SourceThermalCoulomb::resolve_with_query_axes(
            bound.clone(),
            point.clone(),
            ScalarPropertyRequirement::try_new(
                NAMES[i],
                Dims::NONE,
                ScalarAdmissibility::NonNegative,
            )
            .unwrap(),
            ScalarPropertyRequirement::try_new(
                NAMES[i + 2],
                Dims::NONE,
                ScalarAdmissibility::NonNegative,
            )
            .unwrap(),
            MaterialPropertySelection::SingleClaimOnly,
            ThermalFrictionAxes {
                temperature: "Ta".into(),
                nominal_pressure: Some("p".into()),
            },
            vec![
                FrictionQueryAxis {
                    name: "Tb".into(),
                    kind: FrictionQueryAxisKind::Temperature,
                },
                FrictionQueryAxis {
                    name: "speed".into(),
                    kind: FrictionQueryAxisKind::SlidingSpeed,
                },
                FrictionQueryAxis {
                    name: "path".into(),
                    kind: FrictionQueryAxisKind::AccumulatedDistance,
                },
            ],
        )
        .unwrap()
    });
    NormalDirectionalContact::new(
        thermal,
        normal,
        NormalDirectionalInterface {
            point: ThermalInterfacePoint {
                interface: card,
                materials: [&cards[0], &cards[1]],
                context: point,
                temperature_axes: ["Ta", "Tb"],
                liquid_fraction_axes: [None, None],
                mean_temperature_axis: None,
            },
            tangent_frame_id: "original tangent plane",
            partition,
            outside_destination: Some("named exterior reservoir"),
        },
        SourceThermalEllipticCoulomb::from_principal_sources(sources, axis).unwrap(),
        0.0,
        tangent(),
    )
    .unwrap()
}
fn temperatures(
    thermal: &DepositionSimulation<'_>,
    curves: &[ResolvedHeatCapacityCurve; 2],
) -> [f64; 2] {
    [0, 1].map(|i| {
        curves[i]
            .curve()
            .state_at_specific_enthalpy(thermal.cells()[&CELLS[i]].specific_h_j_kg)
            .unwrap()
            .temperature_k()
    })
}

#[test]
fn held_patch_recovers_directional_work_and_exports_only_the_declared_share() {
    let (a, ca) = common::material(1.0, 2.0, 0.0);
    let (b, cb) = common::material(1.0, 3.0, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let point = context();
    let card = interface(&cards, fixed(0.2), 1e6);
    let thermal = common::initial(&curves, 0.0, None);
    let contact = joint(
        &thermal,
        &cards,
        &curves,
        &card,
        &point,
        0.002,
        0.0,
        [0.0; 2],
        HeatPartition::new(0.3, 0.5, 0.2).unwrap(),
        [1.0, 0.0],
    );
    let slip = [0.003, 0.004];
    let next = contact
        .step_candidate(
            &thermal,
            0.1,
            0.002,
            slip,
            controls(),
            tangent(),
            &guess(&thermal),
            other,
            || false,
        )
        .unwrap();
    let normal = next.report().normal.endpoint.normal_force_n;
    let work = normal * (0.2_f64 * slip[0]).hypot(0.4 * slip[1]);
    let endpoint = &next.report().normal.source;
    near(endpoint.friction.dissipated_j, work);
    near(next.report().tangential_work_j, work);
    near(next.report().normal.work.work_into_contact_j, 0.0);
    near(
        next.report().normal.thermal.stored_energy_change_j,
        0.8 * work,
    );
    near(endpoint.partition.other_work_j(), 0.2 * work);
    near(next.report().external_heat_j, 0.0);
    near(next.report().total_energy_residual_j, 0.0);
    assert_eq!(endpoint.friction.slip_increment_m, slip);
    let t = temperatures(next.thermal(), &curves);
    near(t[0], 300.0 + 0.3 * work / 2.0);
    near(t[1], 300.0 + 0.5 * work / 3.0);
    for source in next.contact().state().friction_sources() {
        near(source.state_coordinate("path").unwrap(), 0.005);
        near(source.state_coordinate("speed").unwrap(), 0.05);
    }
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(contact.accumulated_distance_m(), 0.0);
}

#[test]
fn evolving_hertz_area_and_pressure_drive_the_original_friction_claim() {
    let (a, ca) = common::material(1.0, 10.0, 0.0);
    let (b, cb) = common::material(1.0, 20.0, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let point = context();
    let card = interface(
        &cards,
        curve("p", Pressure::DIMS, [(0.0, 0.1), (1e6, 1.1)]),
        1e6,
    );
    let thermal = common::initial(&curves, 0.0, None);
    let contact = joint(
        &thermal,
        &cards,
        &curves,
        &card,
        &point,
        0.002,
        2.0,
        [0.0; 2],
        HeatPartition::new(0.4, 0.6, 0.0).unwrap(),
        [1.0, 0.0],
    );
    let next = contact
        .step_candidate(
            &thermal,
            1.0,
            0.008,
            [0.01, 0.0],
            controls(),
            tangent(),
            &guess(&thermal),
            other,
            || false,
        )
        .unwrap();
    let force = (4.0 / 3.0) * 1e6 * 0.008_f64.powf(1.5) * (1.0 + 2.0 * 0.006);
    let area = std::f64::consts::PI * 0.008;
    let pressure = force / area;
    let heat = (0.1 + pressure / 1e6) * force * 0.01;
    let report = &next.report().normal;
    let endpoint = &report.source;
    near(report.endpoint.normal_force_n, force);
    near(endpoint.contact_area_m2, area);
    near(endpoint.pressure_pa, pressure);
    near(endpoint.friction.dissipated_j, heat);
    assert!(report.work.physical_dissipation_j > 0.0);
    near(
        report.thermal.stored_energy_change_j,
        heat + report.work.physical_dissipation_j,
    );
    near(
        report.work.elastic_work_j,
        1e6 * (basis(0.008) - basis(0.002)),
    );
    near(next.report().total_energy_residual_j, 0.0);
    assert_ne!(
        contact.normal().endpoint().patch_radius_m,
        report.endpoint.patch_radius_m
    );
}

#[test]
fn friction_heat_refreshes_source_temperature_and_reciprocal_normal_stiffness() {
    let (a, ca) = common::material(1.0, 2.0, 0.0);
    let (b, cb) = common::material(1.0, 3.0, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let point = context();
    let card = interface(
        &cards,
        curve("Ta", Temperature::DIMS, [(200.0, 0.1), (600.0, 0.9)]),
        1e6,
    );
    let thermal = common::initial(&curves, 0.0, None);
    let constant = joint(
        &thermal,
        &cards,
        &curves,
        &card,
        &point,
        0.002,
        0.0,
        [0.0; 2],
        HeatPartition::new(0.4, 0.6, 0.0).unwrap(),
        [1.0, 0.0],
    );
    let next = constant
        .step_candidate(
            &thermal,
            1.0,
            0.002,
            [0.02, 0.0],
            controls(),
            tangent(),
            &guess(&thermal),
            other,
            || false,
        )
        .unwrap();
    let ns = constant.normal().endpoint().normal_force_n * 0.02;
    let expected = 0.3 * ns / (1.0 - 0.002 * 0.4 * ns / 2.0);
    near(next.report().tangential_work_j, expected);
    let soft = joint(
        &thermal,
        &cards,
        &curves,
        &card,
        &point,
        0.002,
        0.0,
        [-1000.0, -500.0],
        HeatPartition::new(0.4, 0.6, 0.0).unwrap(),
        [1.0, 0.0],
    );
    let warm = soft
        .step_candidate(
            &thermal,
            1.0,
            0.008,
            [0.02, 0.0],
            controls(),
            tangent(),
            &guess(&thermal),
            other,
            || false,
        )
        .unwrap();
    let cold = soft
        .step_candidate(
            &thermal,
            1.0,
            0.008,
            [0.0; 2],
            controls(),
            tangent(),
            &guess(&thermal),
            other,
            || false,
        )
        .unwrap();
    assert!(
        warm.contact().normal().temperature_k()[0] > cold.contact().normal().temperature_k()[0]
    );
    assert!(
        warm.contact().normal().endpoint().normal_force_n
            < cold.contact().normal().endpoint().normal_force_n
    );
    let t = warm.contact().normal().temperature_k();
    near(
        warm.report().normal.modulus.modulus_pa,
        1e6 - 1000.0 * (t[0] - 300.0) - 500.0 * (t[1] - 300.0),
    );
    near(
        warm.contact()
            .state()
            .constitutive_law()
            .kinetic_coefficients()[0],
        0.1 + 0.002 * (t[0] - 200.0),
    );
    near(warm.report().total_energy_residual_j, 0.0);
}

#[test]
fn changing_load_refines_to_continuous_work_and_preserves_rotated_path_history() {
    let (a, ca) = common::material(1.0, 10.0, 0.0);
    let (b, cb) = common::material(1.0, 20.0, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let point = context();
    let card = interface(&cards, fixed(0.2), 1e6);
    let mut errors = Vec::new();
    let expected = 0.2
        * 0.01
        * (4.0 / 3.0)
        * 1e6
        * (2.0 / (5.0 * 0.006))
        * (0.008_f64.powf(2.5) - 0.002_f64.powf(2.5));
    for steps in [4, 8, 16] {
        let mut thermal = common::initial(&curves, 0.0, None);
        let mut contact = joint(
            &thermal,
            &cards,
            &curves,
            &card,
            &point,
            0.002,
            0.0,
            [0.0; 2],
            HeatPartition::new(0.4, 0.6, 0.0).unwrap(),
            [0.6, 0.8],
        );
        let mut total = 0.0;
        for i in 1..=steps {
            let t = i as f64 / steps as f64;
            let next = contact
                .step_candidate(
                    &thermal,
                    t,
                    0.002 + 0.006 * t,
                    [0.006 / steps as f64, 0.008 / steps as f64],
                    controls(),
                    tangent(),
                    &guess(&thermal),
                    other,
                    || false,
                )
                .unwrap();
            total += next.report().tangential_work_j;
            for source in next.contact().state().friction_sources() {
                assert!(
                    source.state_coordinate("path").unwrap()
                        <= next.contact().accumulated_distance_m()
                );
            }
            (thermal, contact, _) = next.into_parts();
        }
        near(contact.accumulated_distance_m(), 0.01);
        errors.push((total - expected).abs());
    }
    assert!(
        errors[1] < 0.6 * errors[0] && errors[2] < 0.6 * errors[1],
        "{errors:?}"
    );
}

#[test]
fn failed_source_cancel_and_budget_preserve_normal_friction_and_caloric_owners() {
    let (a, ca) = common::material(1.0, 10.0, 0.0);
    let (b, cb) = common::material(1.0, 20.0, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let point = context();
    let card = interface(&cards, fixed(0.2), 1e6);
    let thermal = common::initial(&curves, 0.0, None);
    let contact = joint(
        &thermal,
        &cards,
        &curves,
        &card,
        &point,
        0.001,
        0.0,
        [0.0; 2],
        HeatPartition::new(0.4, 0.6, 0.0).unwrap(),
        [1.0, 0.0],
    );
    let state = contact.state().clone();
    let patch = contact.normal().endpoint().clone();
    let cells = thermal.cells().clone();
    let mut limited = controls();
    limited.maximum_contact_trials = 1;
    assert!(matches!(
        contact.step_candidate(
            &thermal,
            1.0,
            0.008,
            [0.01, 0.0],
            limited,
            tangent(),
            &guess(&thermal),
            other,
            || false
        ),
        Err(NormalThermalError::TrialBudget)
    ));
    assert!(matches!(
        contact.step_candidate(
            &thermal,
            1.0,
            0.008,
            [0.01, 0.0],
            controls(),
            tangent(),
            &guess(&thermal),
            |_| Err::<StepCoefficients<()>, _>("other heat unavailable"),
            || false
        ),
        Err(NormalThermalError::Source(
            DirectionalEndpointError::Source("other heat unavailable")
        ))
    ));
    let calls = Cell::new(0usize);
    let good = contact
        .step_candidate(
            &thermal,
            1.0,
            0.008,
            [0.01, 0.0],
            controls(),
            tangent(),
            &guess(&thermal),
            other,
            || {
                calls.set(calls.get() + 1);
                false
            },
        )
        .unwrap();
    let last = calls.get();
    calls.set(0);
    assert!(matches!(
        contact.step_candidate(
            &thermal,
            1.0,
            0.008,
            [0.01, 0.0],
            controls(),
            tangent(),
            &guess(&thermal),
            other,
            || {
                calls.set(calls.get() + 1);
                calls.get() == last
            }
        ),
        Err(NormalThermalError::Cancelled)
    ));
    let bounded = interface(&cards, fixed(0.2), 15000.0);
    let limited = joint(
        &thermal,
        &cards,
        &curves,
        &bounded,
        &point,
        0.001,
        0.0,
        [0.0; 2],
        HeatPartition::new(0.4, 0.6, 0.0).unwrap(),
        [1.0, 0.0],
    );
    assert!(matches!(
        limited.step_candidate(
            &thermal,
            1.0,
            0.008,
            [0.01, 0.0],
            controls(),
            tangent(),
            &guess(&thermal),
            other,
            || false
        ),
        Err(NormalThermalError::Source(
            DirectionalEndpointError::Friction(_)
        ))
    ));
    let retry = contact
        .step_candidate(
            &thermal,
            1.0,
            0.008,
            [0.01, 0.0],
            controls(),
            tangent(),
            &guess(&thermal),
            other,
            || false,
        )
        .unwrap();
    assert_eq!(retry.contact().state(), good.contact().state());
    assert_eq!(contact.state(), &state);
    assert_eq!(contact.normal().endpoint(), &patch);
    assert_eq!(thermal.cells(), &cells);
}
