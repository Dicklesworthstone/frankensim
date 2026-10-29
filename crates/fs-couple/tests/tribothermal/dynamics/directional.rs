//! Original directional source claims, finite body inertia and caloric heat.
use super::*;
use fs_contact::interface_binding::bind_dry_interface_system_card;
use fs_contact::interface_binding::thermal_friction::elliptic::SourceThermalEllipticCoulomb;
use fs_contact::interface_binding::thermal_friction::{
    FrictionQueryAxis, FrictionQueryAxisKind, SourceThermalCoulomb, ThermalFrictionAxes,
};
use fs_couple::tribothermal::directional::dynamics::DirectionalPair;
use fs_couple::tribothermal::directional::{
    DirectionalSlidingContact, DirectionalSlidingInterface,
};
use fs_couple::tribothermal::dynamics::planar::PlanarMotionSolve;
use fs_material::state_point::{ScalarAdmissibility, ScalarPropertyRequirement};

const STATIC: [&str; 2] = ["principal-static-x", "principal-static-y"];
const KINETIC: [&str; 2] = ["principal-kinetic-x", "principal-kinetic-y"];
fn fixed(mu: f64) -> PropertyValue {
    PropertyValue::Scalar {
        value: mu,
        dims: Dims::NONE,
    }
}
fn absolute_temperature() -> QuantitySpec {
    use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}
fn static_pressure() -> QuantitySpec {
    use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
    QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static))
}
fn curve(axis: &str, dims: Dims, knots: [(f64, f64); 2]) -> PropertyValue {
    PropertyValue::Curve {
        abscissa: axis.into(),
        abscissa_dims: dims,
        knots: knots.to_vec(),
        dims: Dims::NONE,
    }
}
fn query(dynamic: bool) -> QueryPoint {
    let mut q = QueryPoint::new()
        .with_quantity("Ta", absolute_temperature(), 300.0)
        .unwrap()
        .with_quantity("Tb", absolute_temperature(), 300.0)
        .unwrap()
        .with_quantity("pressure", static_pressure(), 10.0)
        .unwrap();
    if dynamic {
        q = q
            .with_quantity("speed", QuantitySpec::dimensional(SPEED), 1.0)
            .unwrap()
            .with_quantity("distance", QuantitySpec::dimensional(LENGTH), 0.0)
            .unwrap();
    }
    q
}
fn directional_card(
    cards: &[MaterialCard; 2],
    static_mu: [f64; 2],
    kinetic: [PropertyValue; 2],
    dynamic: bool,
    upper_temperature: f64,
) -> InterfaceSystemCard {
    let mut claims = ClaimSet::new();
    for (i, name) in STATIC.into_iter().chain(KINETIC).enumerate() {
        let value = if i < 2 {
            fixed(static_mu[i])
        } else {
            kinetic[i - 2].clone()
        };
        let interpolation = if matches!(value, PropertyValue::Scalar { .. }) {
            InterpolationPolicy::ConstantWithinValidity
        } else {
            InterpolationPolicy::LinearInside
        };
        let mut validity = ValidityDomain::unconstrained()
            .with_quantity("Ta", absolute_temperature(), 200.0, upper_temperature)
            .with_quantity("Tb", absolute_temperature(), 200.0, upper_temperature)
            .with_quantity("pressure", static_pressure(), 0.0, 1e5);
        if dynamic {
            validity = validity
                .with_quantity("speed", QuantitySpec::dimensional(SPEED), 0.0, 10.0)
                .with_quantity("distance", QuantitySpec::dimensional(LENGTH), 0.0, 100.0);
        }
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, Dims::NONE),
                value,
                validity,
                interpolation,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic coupled directional fixture, not calibrated data".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    InterfaceSystemCard::assemble(
        SurfaceSpec {
            material: cards[0].id().clone(),
            texture_frame: "directional-A".into(),
        },
        SurfaceSpec {
            material: cards[1].id().clone(),
            texture_frame: "directional-B".into(),
        },
        SystemContext {
            medium: "dry".into(),
            third_body: None,
            environment: "synthetic".into(),
            history: "virgin".into(),
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn solve_controls(step: f64) -> PlanarMotionSolve {
    let mut sliding = controls();
    sliding.tangential.max_iterations = 128;
    sliding.tangential.force_absolute_tolerance_n = 1e-12;
    sliding.tangential.force_relative_tolerance = 1e-13;
    PlanarMotionSolve {
        sliding,
        maximum_step_s: step,
        maximum_segments: 512,
        maximum_motion_trials: 20_000,
        maximum_halvings: 512,
        maximum_bisections: 96,
        impulse_tolerance_n_s: 1e-8,
        event_time_tolerance_s: 1e-10,
        total_energy_tolerance_j: 1e-6,
    }
}
#[allow(clippy::too_many_arguments)]
fn contact<'a>(
    thermal: &DepositionSimulation<'_>,
    cards: &'a [MaterialCard; 2],
    curves: &'a [ResolvedHeatCapacityCurve; 2],
    card: &'a InterfaceSystemCard,
    point: &'a QueryPoint,
    axis: [f64; 2],
    dynamic: bool,
    initial_speed: f64,
    partition: HeatPartition,
) -> DirectionalSlidingContact<'a> {
    let mut owned = vec![FrictionQueryAxis {
        name: "Tb".into(),
        kind: FrictionQueryAxisKind::Temperature,
    }];
    if dynamic {
        owned.push(FrictionQueryAxis {
            name: "speed".into(),
            kind: FrictionQueryAxisKind::SlidingSpeed,
        });
        owned.push(FrictionQueryAxis {
            name: "distance".into(),
            kind: FrictionQueryAxisKind::AccumulatedDistance,
        });
    }
    let sources = [0, 1].map(|i| {
        SourceThermalCoulomb::resolve_with_query_axes(
            bind_dry_interface_system_card(card, InputAuthority::SyntheticFixture).unwrap(),
            point.clone(),
            ScalarPropertyRequirement::try_new(
                STATIC[i],
                Dims::NONE,
                ScalarAdmissibility::NonNegative,
            )
            .unwrap(),
            ScalarPropertyRequirement::try_new(
                KINETIC[i],
                Dims::NONE,
                ScalarAdmissibility::NonNegative,
            )
            .unwrap(),
            MaterialPropertySelection::SingleClaimOnly,
            ThermalFrictionAxes {
                temperature: "Ta".into(),
                nominal_pressure: Some("pressure".into()),
            },
            owned.clone(),
        )
        .unwrap()
    });
    DirectionalSlidingContact::new(
        thermal,
        DirectionalSlidingInterface {
            point: ThermalInterfacePoint {
                interface: card,
                materials: [&cards[0], &cards[1]],
                context: point,
                temperature_axes: ["Ta", "Tb"],
                liquid_fraction_axes: [None, None],
                mean_temperature_axis: None,
            },
            caloric: [curves[0].curve(), curves[1].curve()],
            cells: CELLS,
            nominal_area_m2: 1.0,
            tangent_frame_id: "declared test tangent plane",
            partition,
            outside_destination: Some("declared external thermal body"),
        },
        SourceThermalEllipticCoulomb::from_principal_sources(sources, axis).unwrap(),
        10.0,
        initial_speed,
        solve_controls(1.0).sliding.tangential,
    )
    .unwrap()
}
fn free_drive() -> MotionDrive {
    MotionDrive {
        normal_load_n: 10.0,
        body_forces_n: [[0.0; 2]; 2],
    }
}
fn half_heat() -> HeatPartition {
    HeatPartition::new(0.5, 0.5, 0.0).unwrap()
}

#[test]
fn an_unmapped_primary_pressure_is_not_invented_for_the_source_owner() {
    let (a, ca) = common::material(1.0, 2.0, 0.0);
    let (b, cb) = common::material(1.0, 3.0, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let point = query(false);
    let thermal = inventory(&curves, [1.0, 1.0]);
    let card = directional_card(&cards, [0.6, 0.6], [fixed(0.2), fixed(0.4)], false, 2000.0);
    let sources = [0, 1].map(|i| {
        SourceThermalCoulomb::resolve_with_query_axes(
            bind_dry_interface_system_card(&card, InputAuthority::SyntheticFixture).unwrap(),
            point.clone(),
            ScalarPropertyRequirement::try_new(
                STATIC[i],
                Dims::NONE,
                ScalarAdmissibility::NonNegative,
            )
            .unwrap(),
            ScalarPropertyRequirement::try_new(
                KINETIC[i],
                Dims::NONE,
                ScalarAdmissibility::NonNegative,
            )
            .unwrap(),
            MaterialPropertySelection::SingleClaimOnly,
            ThermalFrictionAxes {
                temperature: "Ta".into(),
                nominal_pressure: None,
            },
            vec![FrictionQueryAxis {
                name: "Tb".into(),
                kind: FrictionQueryAxisKind::Temperature,
            }],
        )
        .unwrap()
    });
    let owner = DirectionalSlidingContact::new(
        &thermal,
        DirectionalSlidingInterface {
            point: ThermalInterfacePoint {
                interface: &card,
                materials: [&cards[0], &cards[1]],
                context: &point,
                temperature_axes: ["Ta", "Tb"],
                liquid_fraction_axes: [None, None],
                mean_temperature_axis: None,
            },
            caloric: [curves[0].curve(), curves[1].curve()],
            cells: CELLS,
            nominal_area_m2: 1.0,
            tangent_frame_id: "explicit plane",
            partition: half_heat(),
            outside_destination: None,
        },
        SourceThermalEllipticCoulomb::from_principal_sources(sources, [1.0, 0.0]).unwrap(),
        10.0,
        1.0,
        solve_controls(1.0).sliding.tangential,
    )
    .unwrap();
    let next = owner
        .step_candidate(
            &thermal,
            0.1,
            SlidingDrive {
                normal_load_n: 10.0,
                relative_increment_m: [0.03, 0.04],
            },
            solve_controls(0.1).sliding,
            &[0.0; 2],
            other,
            || false,
        )
        .unwrap();
    close(next.report().deposited_heat_j, 0.06_f64.hypot(0.16));
    close(next.report().total_energy_residual_j, 0.0);
    assert_eq!(
        next.contact().source_law().principal_sources()[0]
            .axes()
            .nominal_pressure,
        None
    );
}

#[test]
fn principal_arrest_retains_sourced_masses_receipts_and_full_physical_heat() {
    for axis in [[1.0, 0.0], [0.6, 0.8]] {
        let (a, ca, body_a) = sourced_body(2.0, 2.0);
        let (b, cb, body_b) = sourced_body(3.0, 3.0);
        let cards = [a, b];
        let curves = [ca, cb];
        let point = query(false);
        let thermal = inventory_with_geometry(&curves, [2.0, 3.0], [1.0, 1.0]);
        let card = directional_card(&cards, [0.6, 0.8], [fixed(0.2), fixed(0.5)], false, 2000.0);
        let contact = contact(
            &thermal,
            &cards,
            &curves,
            &card,
            &point,
            axis,
            false,
            1.0,
            HeatPartition::new(0.4, 0.6, 0.0).unwrap(),
        );
        let pair = DirectionalPair::from_material_bodies(
            &thermal,
            contact,
            [body_a, body_b],
            [axis, [0.0; 2]],
        )
        .unwrap();
        let next = pair
            .step_planar_candidate(
                &thermal,
                1.0,
                free_drive(),
                solve_controls(1.0),
                &[0.0; 2],
                other,
                || false,
            )
            .unwrap();
        assert!(next.pair().material_bodies().is_some());
        assert_eq!(next.report().segments.len(), 2);
        assert!(next.report().segments[0].reached_stop);
        close(next.report().segments[0].end_time_s, 0.6);
        close(next.report().kinetic_energy_change_j, -0.6);
        close(next.report().deposited_heat_j, 0.6);
        close(next.report().thermal_energy_change_j, 0.6);
        close(next.pair().contact().accumulated_distance_m(), 0.3);
        close(next.thermal().temperature_k(CELLS[0]).unwrap(), 300.06);
        close(next.thermal().temperature_k(CELLS[1]).unwrap(), 300.04);
        for velocity in next.pair().velocities_mps() {
            for i in 0..2 {
                close(velocity[i], 0.4 * axis[i]);
            }
        }
        for source in next.pair().contact().state().friction_sources() {
            close(source.state_coordinate("Ta").unwrap(), 300.06);
            close(source.state_coordinate("Tb").unwrap(), 300.04);
        }
        assert_eq!(
            pair.contact().state().constitutive_state().slip_m(),
            [0.0; 2]
        );
    }
}

fn turning_reference() -> [f64; 5] {
    fn rhs(x: [f64; 5]) -> [f64; 5] {
        let support = (2.0 * x[0]).hypot(4.0 * x[1]);
        [
            -8.0 * x[0] / support,
            -32.0 * x[1] / support,
            x[0],
            x[1],
            x[0].hypot(x[1]),
        ]
    }
    let mut x = [0.8, 0.6, 0.0, 0.0, 0.0];
    let h = 0.08 / 16_384.0;
    for _ in 0..16_384 {
        let a = rhs(x);
        let b = rhs(std::array::from_fn(|i| x[i] + 0.5 * h * a[i]));
        let c = rhs(std::array::from_fn(|i| x[i] + 0.5 * h * b[i]));
        let d = rhs(std::array::from_fn(|i| x[i] + h * c[i]));
        for i in 0..5 {
            x[i] += h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0;
        }
    }
    x
}
#[test]
fn nonparallel_traction_turns_the_body_and_refines_to_independent_continuous_motion() {
    let (a, ca) = common::material(1.0, 2.0, 0.0);
    let (b, cb) = common::material(1.0, 3.0, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let point = query(false);
    let thermal = inventory(&curves, [1.0, 1.0]);
    let card = directional_card(&cards, [0.7, 0.8], [fixed(0.2), fixed(0.4)], false, 2000.0);
    let expected = turning_reference();
    let rotate = |v: [f64; 2]| [0.6 * v[0] - 0.8 * v[1], 0.8 * v[0] + 0.6 * v[1]];
    let mut errors = Vec::new();
    for h in [0.01, 0.005, 0.0025] {
        let mut original = None;
        for rotated in [false, true] {
            let transform = |v| if rotated { rotate(v) } else { v };
            let contact = contact(
                &thermal,
                &cards,
                &curves,
                &card,
                &point,
                transform([1.0, 0.0]),
                false,
                1.0,
                half_heat(),
            );
            let pair = DirectionalPair::new(
                &thermal,
                contact,
                [transform([0.4, 0.3]), transform([-0.4, -0.3])],
            )
            .unwrap();
            let next = pair
                .step_planar_candidate(
                    &thermal,
                    0.08,
                    free_drive(),
                    solve_controls(h),
                    &[0.0; 2],
                    other,
                    || false,
                )
                .unwrap();
            let report = next.report();
            let v = next.pair().velocities_mps();
            let x = [
                v[0][0] - v[1][0],
                v[0][1] - v[1][1],
                report.body_displacements_m[0][0] - report.body_displacements_m[1][0],
                report.body_displacements_m[0][1] - report.body_displacements_m[1][1],
                next.pair().contact().accumulated_distance_m(),
            ];
            close(
                report.deposited_heat_j,
                0.25 * (1.0 - x[0] * x[0] - x[1] * x[1]),
            );
            close(report.total_energy_residual_j, 0.0);
            let first = &report.segments[0].sliding.friction;
            assert!(
                (first.traction_n[0] * first.slip_increment_m[1]
                    - first.traction_n[1] * first.slip_increment_m[0])
                    .abs()
                    > 1e-5
            );
            if rotated {
                let old: [f64; 5] = original.unwrap();
                let v = rotate([old[0], old[1]]);
                let p = rotate([old[2], old[3]]);
                for (a, b) in x.into_iter().zip([v[0], v[1], p[0], p[1], old[4]]) {
                    close(a, b);
                }
            } else {
                errors.push(
                    x.into_iter()
                        .zip(expected)
                        .map(|(a, b)| (a - b).abs())
                        .fold(0.0, f64::max),
                );
                original = Some(x);
            }
        }
    }
    assert!(
        errors[1] < 0.35 * errors[0] && errors[2] < 0.35 * errors[1],
        "{errors:?}"
    );
    assert!(errors[2] < 1e-4, "{errors:?}");
}

#[test]
fn vector_static_ellipse_and_moving_kinetic_branch_remain_distinct() {
    let (a, ca) = common::material(1.0, 2.0, 0.0);
    let (b, cb) = common::material(1.0, 3.0, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let point = query(false);
    let thermal = inventory(&curves, [1.0, 1.0]);
    let card = directional_card(&cards, [0.6, 0.3], [fixed(0.2), fixed(0.1)], false, 2000.0);
    let owner = contact(
        &thermal,
        &cards,
        &curves,
        &card,
        &point,
        [1.0, 0.0],
        false,
        0.0,
        half_heat(),
    );
    let rest = DirectionalPair::new(&thermal, owner.clone(), [[0.0; 2]; 2]).unwrap();
    for (force, sticking) in [([4.0, 1.0], true), ([6.0, 1.0], false)] {
        let next = rest
            .step_planar_candidate(
                &thermal,
                0.1,
                MotionDrive {
                    normal_load_n: 10.0,
                    body_forces_n: [force, force.map(|x| -x)],
                },
                solve_controls(0.1),
                &[0.0; 2],
                other,
                || false,
            )
            .unwrap();
        let response = &next.report().segments[0].sliding.friction;
        assert_eq!(
            response.regime == fs_tribo::FrictionRegime::Sticking,
            sticking
        );
        if sticking {
            assert_eq!(next.pair().velocities_mps(), [[0.0; 2]; 2]);
            close(next.report().deposited_heat_j, 0.0);
        } else {
            close(
                (response.traction_n[0] / 2.0).powi(2) + response.traction_n[1].powi(2),
                1.0,
            );
        }
        close(
            next.report().kinetic_energy_change_j + next.report().deposited_heat_j,
            next.report().external_work_j,
        );
    }
    // The midpoint released traction is 5 N: below static 6, above kinetic 2.
    // Existing motion must keep sliding rather than be captured by static-first.
    let moving = DirectionalPair::new(&thermal, owner, [[0.005, 0.0], [-0.005, 0.0]]).unwrap();
    let next = moving
        .step_planar_candidate(
            &thermal,
            0.002,
            free_drive(),
            solve_controls(0.002),
            &[0.0; 2],
            other,
            || false,
        )
        .unwrap();
    assert_eq!(next.report().segments.len(), 1);
    close(
        next.pair().velocities_mps()[0][0] - next.pair().velocities_mps()[1][0],
        0.002,
    );
    close(next.report().deposited_heat_j, 0.000024);
    let unloaded = rest
        .step_planar_candidate(
            &thermal,
            0.1,
            MotionDrive {
                normal_load_n: 0.0,
                body_forces_n: [[2.0, 0.0], [2.0, 0.0]],
            },
            solve_controls(0.1),
            &[0.0; 2],
            other,
            || false,
        )
        .unwrap();
    for v in unloaded.pair().velocities_mps() {
        close(v[0], 0.2);
        close(v[1], 0.0);
    }
    close(unloaded.report().deposited_heat_j, 0.0);
    close(unloaded.report().kinetic_energy_change_j, 0.04);
    close(unloaded.report().external_work_j, 0.04);
}

#[test]
fn both_original_surface_temperatures_refresh_directional_sources_and_change_motion() {
    let mut speeds = Vec::new();
    for capacity in [0.1, 1.0] {
        let (a, ca) = common::material(1.0, capacity, 0.0);
        let (b, cb) = common::material(1.0, 3.0, 0.0);
        let cards = [a, b];
        let curves = [ca, cb];
        let point = query(false);
        let thermal = inventory(&curves, [1.0, 1.0]);
        let card = directional_card(
            &cards,
            [2.0, 2.0],
            [
                curve("Ta", Temperature::DIMS, [(299.0, 0.1), (310.0, 1.2)]),
                curve("Tb", Temperature::DIMS, [(299.0, 0.05), (310.0, 0.6)]),
            ],
            false,
            2000.0,
        );
        let owner = contact(
            &thermal,
            &cards,
            &curves,
            &card,
            &point,
            [1.0, 0.0],
            false,
            1.0,
            half_heat(),
        );
        let pins = owner.state().friction_sources().each_ref().map(|s| {
            s.properties()
                .iter()
                .map(|p| p.answer().receipt.selected)
                .collect::<Vec<_>>()
        });
        let pair = DirectionalPair::new(&thermal, owner, [[0.5, 0.0], [-0.5, 0.0]]).unwrap();
        let next = pair
            .step_planar_candidate(
                &thermal,
                0.1,
                free_drive(),
                solve_controls(0.1),
                &[0.0; 2],
                other,
                || false,
            )
            .unwrap();
        let k = 0.05 / capacity;
        let mu = 0.4 / ((1.0 - k) + ((1.0 - k) * (1.0 - k) + 0.8 * k).sqrt());
        let heat = mu * (1.0 - mu);
        let velocity = next.pair().velocities_mps();
        let speed = velocity[0][0] - velocity[1][0];
        speeds.push(speed);
        close(speed, 1.0 - 2.0 * mu);
        close(next.report().deposited_heat_j, heat);
        let ta = next.thermal().temperature_k(CELLS[0]).unwrap();
        let tb = next.thermal().temperature_k(CELLS[1]).unwrap();
        close(ta, 300.0 + heat / (2.0 * capacity));
        close(tb, 300.0 + heat / 6.0);
        for (i, source) in next
            .pair()
            .contact()
            .state()
            .friction_sources()
            .iter()
            .enumerate()
        {
            close(source.state_coordinate("Ta").unwrap(), ta);
            close(source.state_coordinate("Tb").unwrap(), tb);
            assert_eq!(
                source
                    .properties()
                    .iter()
                    .map(|p| p.answer().receipt.selected)
                    .collect::<Vec<_>>(),
                pins[i]
            );
        }
        close(
            next.pair()
                .contact()
                .state()
                .constitutive_law()
                .kinetic_coefficients()[1],
            0.1 + 0.05 * (tb - 300.0),
        );
    }
    assert!(speeds[1] > speeds[0] + 0.1);
}

#[test]
fn actual_speed_and_unsigned_path_close_the_same_native_motion_and_heat() {
    let (a, ca) = common::material(1.0, 2.0, 0.0);
    let (b, cb) = common::material(1.0, 3.0, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let point = query(true);
    let thermal = inventory(&curves, [1.0, 1.0]);
    let card = directional_card(
        &cards,
        [2.0, 2.0],
        [
            curve("speed", SPEED, [(0.0, 0.1), (4.0, 0.9)]),
            curve("distance", LENGTH, [(0.0, 0.3), (1.0, 0.7)]),
        ],
        true,
        2000.0,
    );
    let owner = contact(
        &thermal,
        &cards,
        &curves,
        &card,
        &point,
        [1.0, 0.0],
        true,
        1.0,
        half_heat(),
    );
    let pair = DirectionalPair::new(&thermal, owner, [[0.5, 0.0], [-0.5, 0.0]]).unwrap();
    let next = pair
        .step_planar_candidate(
            &thermal,
            0.1,
            free_drive(),
            solve_controls(0.1),
            &[0.0; 2],
            other,
            || false,
        )
        .unwrap();
    close(
        next.pair().velocities_mps()[0][0] - next.pair().velocities_mps()[1][0],
        0.5,
    );
    close(next.pair().contact().accumulated_distance_m(), 0.075);
    close(next.report().deposited_heat_j, 0.1875);
    let mu = next
        .pair()
        .contact()
        .state()
        .constitutive_law()
        .kinetic_coefficients();
    close(mu[0], 0.25);
    close(mu[1], 0.33);
    for source in next.pair().contact().state().friction_sources() {
        close(source.state_coordinate("speed").unwrap(), 0.75);
        close(source.state_coordinate("distance").unwrap(), 0.075);
    }
    assert!(next.report().segments[0].source_motion_residual_n_s.abs() < 1e-8);
    // Repeated source refresh, then arrest and static continuation, must not
    // rewind the original unsigned history through scalar-root roundoff.
    let follow = next
        .pair()
        .step_planar_candidate(
            next.thermal(),
            0.5,
            free_drive(),
            solve_controls(0.02),
            &[0.0; 2],
            other,
            || false,
        )
        .unwrap();
    let v = follow.pair().velocities_mps();
    close(v[0][0] - v[1][0], 0.0);
    close(
        next.report().deposited_heat_j + follow.report().deposited_heat_j,
        0.25,
    );
    for source in follow.pair().contact().state().friction_sources() {
        close(
            source.state_coordinate("distance").unwrap(),
            follow.pair().contact().accumulated_distance_m(),
        );
    }
}

#[test]
fn frictionless_axis_survives_other_axis_arrest_without_false_static_capture() {
    let (a, ca) = common::material(1.0, 2.0, 0.0);
    let (b, cb) = common::material(1.0, 3.0, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let point = query(false);
    let thermal = inventory(&curves, [1.0, 1.0]);
    let card = directional_card(&cards, [0.6, 0.6], [fixed(0.0), fixed(0.2)], false, 2000.0);
    for free_speed in [1e-10_f64, 1e15] {
        let owner = contact(
            &thermal,
            &cards,
            &curves,
            &card,
            &point,
            [1.0, 0.0],
            false,
            free_speed.hypot(1.0),
            half_heat(),
        );
        let pair = DirectionalPair::new(
            &thermal,
            owner,
            [[0.5 * free_speed, 0.5], [-0.5 * free_speed, -0.5]],
        )
        .unwrap();
        let next = pair
            .step_planar_candidate(
                &thermal,
                1.0,
                free_drive(),
                solve_controls(1.0),
                &[0.0; 2],
                other,
                || false,
            )
            .unwrap();
        let v = next.pair().velocities_mps();
        assert_eq!(v[0][0] - v[1][0], free_speed);
        close(v[0][1] - v[1][1], 0.0);
        close(next.report().deposited_heat_j, 0.25);
        assert!(!next.report().segments.iter().any(|s| s.reached_stop));
        assert!(
            (next.pair().contact().state().constitutive_state().slip_m()[0] - free_speed).abs()
                < 1e-12 * free_speed
        );
    }
}

#[test]
fn source_refusal_cancellation_and_motion_budget_never_publish_partial_owners() {
    let (a, ca) = common::material(1.0, 0.1, 0.0);
    let (b, cb) = common::material(1.0, 3.0, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let point = query(false);
    let thermal = inventory(&curves, [1.0, 1.0]);
    let card = directional_card(&cards, [0.6, 0.6], [fixed(0.2), fixed(0.4)], false, 2000.0);
    let owner = contact(
        &thermal,
        &cards,
        &curves,
        &card,
        &point,
        [1.0, 0.0],
        false,
        1.0,
        half_heat(),
    );
    let pair = DirectionalPair::new(&thermal, owner, [[0.5, 0.0], [-0.5, 0.0]]).unwrap();
    let old = thermal.cells().clone();
    let state = pair.contact().state().clone();
    let mut limited = solve_controls(1.0);
    limited.maximum_motion_trials = 1;
    assert!(matches!(
        pair.step_planar_candidate(
            &thermal,
            1.0,
            free_drive(),
            limited,
            &[0.0; 2],
            other,
            || false
        ),
        Err(SlidingError::MotionBudget)
    ));
    assert!(matches!(
        pair.step_planar_candidate(
            &thermal,
            1.0,
            free_drive(),
            solve_controls(1.0),
            &[0.0; 2],
            |_| Err::<StepCoefficients<()>, _>("unavailable source"),
            || false
        ),
        Err(SlidingError::Source("unavailable source"))
    ));
    let polls = Cell::new(0usize);
    let good = pair
        .step_planar_candidate(
            &thermal,
            1.0,
            free_drive(),
            solve_controls(1.0),
            &[0.0; 2],
            other,
            || {
                polls.set(polls.get() + 1);
                false
            },
        )
        .unwrap();
    let count = polls.get();
    polls.set(0);
    assert!(matches!(
        pair.step_planar_candidate(
            &thermal,
            1.0,
            free_drive(),
            solve_controls(1.0),
            &[0.0; 2],
            other,
            || {
                polls.set(polls.get() + 1);
                polls.get() == count
            }
        ),
        Err(SlidingError::Cancelled)
    ));
    let retry = pair
        .step_planar_candidate(
            &thermal,
            1.0,
            free_drive(),
            solve_controls(1.0),
            &[0.0; 2],
            other,
            || false,
        )
        .unwrap();
    assert_eq!(retry.pair().velocities_mps(), good.pair().velocities_mps());
    assert_eq!(thermal.cells(), &old);
    assert_eq!(pair.contact().state(), &state);
    let bounded = directional_card(&cards, [0.6, 0.6], [fixed(0.2), fixed(0.4)], false, 300.01);
    let contact = contact(
        &thermal,
        &cards,
        &curves,
        &bounded,
        &point,
        [1.0, 0.0],
        false,
        1.0,
        half_heat(),
    );
    let bounded_pair = DirectionalPair::new(&thermal, contact, [[0.5, 0.0], [-0.5, 0.0]]).unwrap();
    assert!(matches!(
        bounded_pair.step_planar_candidate(
            &thermal,
            0.1,
            free_drive(),
            solve_controls(0.1),
            &[0.0; 2],
            other,
            || false
        ),
        Err(SlidingError::DirectionalSource(_))
    ));
    assert_eq!(
        bounded_pair.contact().state().constitutive_state().slip_m(),
        [0.0; 2]
    );
}

#[test]
fn oblique_source_paths_never_advance_beyond_accepted_physical_motion() {
    let (a, ca) = common::material(1.0, 2.0, 0.0);
    let (b, cb) = common::material(1.0, 3.0, 0.0);
    let cards = [a, b];
    let curves = [ca, cb];
    let point = query(true);
    let mut thermal = inventory(&curves, [1.0, 1.0]);
    let card = directional_card(
        &cards,
        [2.0, 2.0],
        [
            curve("speed", SPEED, [(0.0, 0.1), (4.0, 0.9)]),
            curve("distance", LENGTH, [(0.0, 0.3), (1.0, 0.7)]),
        ],
        true,
        2000.0,
    );
    let owner = contact(
        &thermal,
        &cards,
        &curves,
        &card,
        &point,
        [0.6, 0.8],
        true,
        1.2_f64.hypot(0.7),
        half_heat(),
    );
    let mut pair = DirectionalPair::new(&thermal, owner, [[0.6, 0.35], [-0.6, -0.35]]).unwrap();
    for i in 1..=8 {
        let next = pair
            .step_planar_candidate(
                &thermal,
                0.017 * i as f64,
                free_drive(),
                solve_controls(0.017),
                &[0.0; 2],
                other,
                || false,
            )
            .unwrap();
        let mut path = pair.contact().accumulated_distance_m();
        for segment in &next.report().segments {
            let slip = segment.sliding.friction.slip_increment_m;
            path += fs_math::det::hypot(slip[0], slip[1]);
            for source in segment.sliding.friction.state.friction_sources() {
                assert!(source.state_coordinate("distance").unwrap() <= path);
            }
        }
        assert_eq!(path, next.pair().contact().accumulated_distance_m());
        (thermal, pair, _) = next.into_parts();
    }
}
