//! Authored multiaxis sources and independent curved-storage temperature limits.

use super::*;
use fs_couple::moisture::{
    atmosphere::{DryBasisSorption, FiniteMoistureAtmosphere},
    radial::{MoistureRadial, MoistureRadialGeometry, RadialMoistureSymmetry},
};

fn curved(
    material: MaterialCard,
    symmetry: RadialMoistureSymmetry,
    faces_m: Vec<f64>,
    initial: Vec<f64>,
    cx: &Cx<'_>,
) -> MoistureRadial {
    MoistureRadial::new(
        material,
        point(0.1),
        "moisture_diffusivity",
        MaterialPropertySelection::SingleClaimOnly,
        MoistureRadialGeometry {
            symmetry,
            faces_m,
            maximum_moisture: 0.3,
            maximum_cells: 256,
        },
        initial,
        cx,
    )
    .unwrap()
}

fn sorption_at(temperature: f64) -> DryBasisSorption {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::with_quantity(
                "sorption",
                qty(QuantityKind::Moisture(MoistureBasis::DryMass)),
            ),
            value: PropertyValue::Curve {
                abscissa: "relative_humidity".into(),
                abscissa_dims: Dims::NONE,
                knots: vec![(0.0, 0.0), (1.0, 0.25)],
                dims: Dims::NONE,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity(
                    "temperature",
                    qty(QuantityKind::AbsoluteTemperature),
                    temperature,
                    temperature,
                )
                .with_quantity(
                    "relative_humidity",
                    qty(QuantityKind::Moisture(MoistureBasis::RelativeHumidity)),
                    0.0,
                    1.0,
                ),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "authored temperature-qualified M_eq=0.25*RH boundary".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    let material = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "authored radial temperature boundary".into(),
            phase: "solid".into(),
            process: "verification".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let point = QueryPoint::new()
        .with_quantity(
            "temperature",
            qty(QuantityKind::AbsoluteTemperature),
            temperature,
        )
        .unwrap()
        .with_quantity(
            "relative_humidity",
            qty(QuantityKind::Moisture(MoistureBasis::RelativeHumidity)),
            0.5,
        )
        .unwrap();
    DryBasisSorption::new(
        material,
        point,
        "sorption",
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

#[test]
fn radial_temperature_events_preserve_unequal_shell_inventory_and_original_sources() {
    with_cx(|cx| {
        for symmetry in [
            RadialMoistureSymmetry::Sphere,
            RadialMoistureSymmetry::Cylinder { length_m: 0.3 },
        ] {
            let start = curved(
                thermal_card(false, 0.3),
                symmetry,
                vec![0.02, 0.035, 0.1],
                vec![0.05, 0.2],
                cx,
            );
            let sealed = [MoistureFace::Sealed, MoistureFace::Sealed];
            let accepted = start
                .advance(0.001, &sealed, config(), nonlinear(), cx)
                .unwrap()
                .0;
            let warm = accepted.at_prescribed_temperature(HOT, cx).unwrap();
            assert_eq!(warm.moisture(), accepted.moisture());
            assert_eq!(warm.time_s(), accepted.time_s());
            assert_eq!(warm.dry_masses_kg(), accepted.dry_masses_kg());
            assert_eq!(warm.water_mass_kg(), accepted.water_mass_kg());
            assert_eq!(warm.wet_mass_kg(), accepted.wet_mass_kg());
            assert_eq!(warm.identity(), accepted.identity());
            assert_eq!(warm.reference().identity(), accepted.reference().identity());
            assert_eq!(warm.geometry().faces_m, accepted.geometry().faces_m);
            assert_eq!(warm.temperature_kelvin(), HOT);
            assert_eq!(warm.reference().state_coordinate("temperature"), Some(COLD));
            let pin = start
                .reference()
                .property("moisture_diffusivity")
                .unwrap()
                .answer()
                .receipt
                .selected;
            assert_eq!(
                warm.transport_context()
                    .property("moisture_diffusivity")
                    .unwrap()
                    .answer()
                    .receipt
                    .selected,
                pin
            );
            near(
                warm.transport_context()
                    .property("moisture_diffusivity")
                    .unwrap()
                    .value_si(),
                0.03 + 0.1 * warm.moisture()[0],
                3e-15,
            );

            let (next, report) = warm
                .advance(0.001, &sealed, config(), nonlinear(), cx)
                .unwrap();
            let cold = accepted
                .advance(0.001, &sealed, config(), nonlinear(), cx)
                .unwrap()
                .0;
            assert!(next.moisture()[0] > cold.moisture()[0]);
            near(next.water_mass_kg(), accepted.water_mass_kg(), 2e-12);
            for state in report.transport_start.iter().chain(&report.transport_end) {
                assert_eq!(state.state_coordinate("temperature"), Some(HOT));
                assert_eq!(
                    state
                        .property("moisture_diffusivity")
                        .unwrap()
                        .answer()
                        .receipt
                        .selected,
                    pin
                );
            }
            let returned = next.at_prescribed_temperature(COLD, cx).unwrap();
            assert_eq!(returned.moisture(), next.moisture());
            assert_eq!(returned.time_s(), next.time_s());
            assert_eq!(returned.dry_masses_kg(), start.dry_masses_kg());
            assert_eq!(
                returned.reference().identity(),
                start.reference().identity()
            );
            near(
                returned
                    .transport_context()
                    .property("moisture_diffusivity")
                    .unwrap()
                    .value_si(),
                0.01 + 0.1 * next.moisture()[0],
                3e-15,
            );
        }
    });
}

fn positive_quadratic(coefficient: f64, rhs: f64) -> f64 {
    2.0 * rhs / (1.0 + (1.0 + 4.0 * coefficient * rhs).sqrt())
}

#[test]
fn radial_hot_dry_front_and_both_time_methods_use_the_current_boundary_temperature() {
    with_cx(|cx| {
        let radius: f64 = 0.1;
        let original = curved(
            thermal_card(true, 0.3),
            RadialMoistureSymmetry::Sphere,
            vec![0.0, radius],
            vec![0.0],
            cx,
        );
        let warm = original.at_prescribed_temperature(HOT, cx).unwrap();
        let boundary = [
            MoistureFace::Sealed,
            MoistureFace::SourceSorption {
                law: sorption_at(HOT),
                relative_humidity: 0.8,
            },
        ];
        let dt = 0.01;
        let equilibrium: f64 = 0.2;
        // Original exact sphere resistance and volume give dM/dt =
        // 3/R^2 * (K(M_eq)-K(M)), with hot D=.3M and K=.15M^2.
        let a = 0.45 * dt / radius.powi(2);
        let (be, be_report) = warm
            .advance(dt, &boundary, config(), nonlinear(), cx)
            .unwrap();
        near(
            be.moisture()[0],
            positive_quadratic(a, a * equilibrium.powi(2)),
            3e-12,
        );
        let gamma = 1.0 - core::f64::consts::FRAC_1_SQRT_2;
        let first = positive_quadratic(gamma * a, gamma * a * equilibrium.powi(2));
        let rhs = (1.0 - gamma) / gamma * first + gamma * a * equilibrium.powi(2);
        let (rk, report) = warm
            .advance_sdirk2(dt, &boundary, config(), nonlinear(), cx)
            .unwrap();
        near(rk.moisture()[0], positive_quadratic(gamma * a, rhs), 3e-12);
        for (field, step) in [(&be, &be_report), (&rk, &report.endpoint)] {
            assert!(field.moisture()[0] > 0.0);
            near(step.boundary_water_kg[1], field.water_mass_kg(), 3e-12);
            near(step.water_balance_residual_kg, 0.0, 3e-12);
            assert_eq!(
                step.transport_end[0].state_coordinate("temperature"),
                Some(HOT)
            );
            assert_eq!(
                step.surface_sorption[1]
                    .as_ref()
                    .unwrap()
                    .state_coordinate("temperature"),
                Some(HOT)
            );
        }
        let stale = [
            MoistureFace::Sealed,
            MoistureFace::SourceSorption {
                law: sorption_at(COLD),
                relative_humidity: 0.8,
            },
        ];
        assert!(warm.advance(dt, &stale, config(), nonlinear(), cx).is_err());
        assert!(
            warm.advance_sdirk2(dt, &stale, config(), nonlinear(), cx)
                .is_err()
        );
        assert_eq!(warm.moisture(), &[0.0]);
        assert_eq!(warm.time_s(), 0.0);
        assert_eq!(original.temperature_kelvin(), COLD);
    });
}

#[test]
fn radial_finite_gas_uses_current_temperature_and_retains_original_water_account() {
    with_cx(|cx| {
        let original = curved(
            thermal_card(true, 0.3),
            RadialMoistureSymmetry::Sphere,
            vec![0.0, 0.1],
            vec![0.0],
            cx,
        );
        let warm = original.at_prescribed_temperature(HOT, cx).unwrap();
        let cold_air = atmosphere::atmosphere(0.8, 1.0);
        let (heated, _) = cold_air.gas().at_temperature(HOT).unwrap();
        let hot_air = FiniteMoistureAtmosphere::new(heated).unwrap();
        let boundaries = [None, Some(sorption_at(HOT))];
        let total_water = hot_air.vapor_mass_kg().unwrap() + warm.water_mass_kg();
        let total_mass = hot_air.gas().mass_kg() + warm.wet_mass_kg();
        for second_order in [false, true] {
            let (air, field, transfer) = if second_order {
                let candidate = hot_air
                    .advance_radial_sdirk2(
                        &warm,
                        0.01,
                        &boundaries,
                        config(),
                        nonlinear(),
                        atmosphere::controls(),
                        cx,
                    )
                    .unwrap();
                let (air, field, report) = candidate.into_parts();
                (air, field, report.water_to_solid_kg)
            } else {
                let candidate = hot_air
                    .advance_radial(
                        &warm,
                        0.01,
                        &boundaries,
                        config(),
                        nonlinear(),
                        atmosphere::controls(),
                        cx,
                    )
                    .unwrap();
                let (air, field, report) = candidate.into_parts();
                (air, field, report.water_to_solid_kg)
            };
            assert!(transfer > 0.0);
            assert!(air.relative_humidity().unwrap() < hot_air.relative_humidity().unwrap());
            assert_eq!(field.temperature_kelvin(), HOT);
            assert_eq!(
                field.reference().state_coordinate("temperature"),
                Some(COLD)
            );
            assert_eq!(field.dry_masses_kg(), original.dry_masses_kg());
            near(
                air.vapor_mass_kg().unwrap() + field.water_mass_kg(),
                total_water,
                3e-12,
            );
            near(air.gas().mass_kg() + field.wet_mass_kg(), total_mass, 3e-12);
            near(field.water_mass_kg(), transfer, 3e-12);
        }
        assert!(
            cold_air
                .advance_radial(
                    &warm,
                    0.01,
                    &[None, None],
                    config(),
                    nonlinear(),
                    atmosphere::controls(),
                    cx,
                )
                .is_err()
        );
        assert!(
            cold_air
                .advance_radial_sdirk2(
                    &warm,
                    0.01,
                    &[None, None],
                    config(),
                    nonlinear(),
                    atmosphere::controls(),
                    cx,
                )
                .is_err()
        );
        assert_eq!(warm.moisture(), &[0.0]);
        assert_eq!(warm.time_s(), 0.0);
        assert_eq!(hot_air.gas().state().state().temperature, HOT);
    });
}

#[test]
fn radial_temperature_event_refuses_unsupported_cells_and_cancellation_atomically() {
    with_cx(|cx| {
        let original = curved(
            thermal_card(false, 0.15),
            RadialMoistureSymmetry::Sphere,
            vec![0.0, 0.02, 0.1],
            vec![0.1, 0.25],
            cx,
        );
        // The first occupied point is supported at HOT. The second is outside
        // the real source triangulation despite passing its enclosing box.
        assert!(original.at_prescribed_temperature(HOT, cx).is_err());
        for invalid in [0.0, -1.0, f64::NAN, f64::INFINITY, HOT + 1.0] {
            assert!(original.at_prescribed_temperature(invalid, cx).is_err());
        }
        assert_eq!(original.moisture(), &[0.1, 0.25]);
        assert_eq!(original.temperature_kelvin(), COLD);
        assert_eq!(original.time_s(), 0.0);
    });
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 3,
                kernel_id: 17,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let original = curved(
            thermal_card(false, 0.3),
            RadialMoistureSymmetry::Sphere,
            vec![0.0, 0.1],
            vec![0.1],
            &cx,
        );
        gate.request();
        assert!(matches!(
            original.at_prescribed_temperature(HOT, &cx),
            Err(MoistureError::Cancelled)
        ));
        assert_eq!(original.moisture(), &[0.1]);
        assert_eq!(original.temperature_kelvin(), COLD);
        assert_eq!(original.time_s(), 0.0);
    });
}
