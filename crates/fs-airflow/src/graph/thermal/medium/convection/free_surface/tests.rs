use super::*;
use fs_lbm::material::free_surface::enthalpy::conduction::{
    FreeSurfaceGasFace, FreeSurfaceInterfaceFace,
};
use fs_material::thermochemical_gas::ResolvedThermochemicalGasState;

// Existing canonical source-pack fixture; no new constitutive model or gas defaults.
#[path = "../../../../../../tests/finite_gas_convection/source.rs"]
mod source;

#[path = "mixed_tests.rs"]
mod mixed;
#[path = "orientation_tests.rs"]
mod orientation;

fn query<'a>(
    wall: f64,
    bulk: &'a ResolvedThermochemicalGasState,
    film: &'a ResolvedThermochemicalGasState,
) -> FreeSurfaceConvectionQuery<'a> {
    FreeSurfaceConvectionQuery {
        patch_index: 0,
        surface: FreeSurfaceInterfaceFace {
            face: FreeSurfaceGasFace {
                liquid: [3, 3, 3],
                gas: [4, 3, 3],
                area_m2: 1e-6,
                distances_m: [0.0005, 0.05],
            },
            emissivity: 0.0.into(),
            surroundings_k: 300.0,
        },
        outward_normal: [1.0, 0.0, 0.0],
        surface_temperature_k: wall,
        bulk,
        film,
    }
}

fn forced() -> PlanarFreeSurfaceConvection {
    PlanarFreeSurfaceConvection {
        outward_normal: [1.0, 0.0, 0.0],
        characteristic_length: Length::new(0.1),
        law: PlanarFreeSurfaceHeatLaw::Forced {
            regime: PlanarForcedConvectionRegime::Laminar,
            relative_velocity_m_s: [0.0, 12.0, 0.0],
        },
    }
}

fn natural(gravity_m_s2: [f64; 3]) -> PlanarFreeSurfaceConvection {
    PlanarFreeSurfaceConvection {
        law: PlanarFreeSurfaceHeatLaw::Natural { gravity_m_s2 },
        ..forced()
    }
}

#[test]
fn forced_plate_retains_the_original_film_correlation_and_pressure_response() {
    let model = source::model_transport(0.75, true);
    let mut coefficients = Vec::new();
    for pressure in [100000.0, 400000.0] {
        let bulk = model.evaluate(&source::point(320.0, pressure)).unwrap();
        let film = model.evaluate(&source::point(400.0, pressure)).unwrap();
        let actual = forced()
            .coefficient(query(480.0, &bulk, &film))
            .unwrap()
            .unwrap();
        let direct = evaluate_forced_from_material(
            CorrelationId::FlatPlateLaminarAverage,
            &film,
            12.0,
            Length::new(0.1),
            ForcedConvectionSettings::default(),
        )
        .unwrap();
        assert_eq!(&actual, direct.coefficient());
        let gas = film.state();
        let re = gas.density * 12.0 * 0.1 / gas.dynamic_viscosity;
        let pr = gas.specific_heat_cp * gas.dynamic_viscosity / gas.thermal_conductivity;
        let analytic = 0.664 * re.sqrt() * pr.cbrt() * gas.thermal_conductivity / 0.1;
        assert!((actual.value.value() / analytic - 1.0).abs() < 1e-12);
        // The gas film distance is geometry input to the lower owner, not a
        // second resistance or a new plate entrance in this correlation law.
        let mut different_film_distance = query(480.0, &bulk, &film);
        different_film_distance.surface.face.distances_m[1] = 5.0;
        assert_eq!(
            forced()
                .coefficient(different_film_distance)
                .unwrap()
                .unwrap(),
            actual
        );
        coefficients.push(actual.value.value());
    }
    assert!((coefficients[1] / coefficients[0] - 2.0).abs() < 1e-12);
}

#[test]
fn source_expansion_and_actual_gravity_orientation_select_natural_plate_heat() {
    let model = source::model(0.75);
    let film = model.evaluate(&source::point(350.0, 100000.0)).unwrap();
    let beta = film
        .thermal_transport_state()
        .unwrap()
        .volumetric_expansion_per_k()
        .unwrap();
    for (wall, ambient) in [(400.0, 300.0), (300.0, 400.0)] {
        let bulk = model.evaluate(&source::point(ambient, 100000.0)).unwrap();
        for (gravity, geometry) in [
            ([0.0, -9.81, 0.0], NaturalConvectionGeometry::VerticalPlate),
            (
                [-9.81, 0.0, 0.0],
                NaturalConvectionGeometry::HorizontalPlate(HorizontalPlateFacing::Up),
            ),
            (
                [9.81, 0.0, 0.0],
                NaturalConvectionGeometry::HorizontalPlate(HorizontalPlateFacing::Down),
            ),
        ] {
            let actual = natural(gravity)
                .coefficient(query(wall, &bulk, &film))
                .unwrap()
                .unwrap();
            let direct = evaluate_natural_from_material(
                geometry,
                &film,
                9.81,
                beta,
                Temperature::new(wall),
                Temperature::new(ambient),
                Length::new(0.1),
            )
            .unwrap();
            assert_eq!(&actual, direct.coefficient());
            assert_eq!(
                (actual.value.value() * (ambient - wall)).is_sign_positive(),
                ambient > wall
            );
        }
    }
}

#[test]
fn exact_equilibrium_admits_physical_inputs_before_returning_zero_heat() {
    let model = source::model(0.75);
    let gas = model.evaluate(&source::point(350.0, 100000.0)).unwrap();
    let q = query(350.0, &gas, &gas);
    assert_eq!(natural([0.0, -9.81, 0.0]).coefficient(q).unwrap(), None);
    for law in [
        natural([0.0; 3]),
        natural([1.0, -9.81, 0.0]),
        natural([f64::NAN, 0.0, 0.0]),
        PlanarFreeSurfaceConvection {
            characteristic_length: Length::new(0.0),
            ..natural([0.0, -9.81, 0.0])
        },
        PlanarFreeSurfaceConvection {
            outward_normal: [-1.0, 0.0, 0.0],
            ..natural([0.0, -9.81, 0.0])
        },
    ] {
        assert!(
            law.coefficient(q).is_err(),
            "equilibrium must not bypass {law:?}"
        );
    }
}

#[test]
fn relative_flow_tangency_regime_and_actual_mean_film_are_enforced() {
    let model = source::model(0.75);
    let bulk = model.evaluate(&source::point(300.0, 100000.0)).unwrap();
    let film = model.evaluate(&source::point(350.0, 100000.0)).unwrap();
    let q = query(400.0, &bulk, &film);
    for (regime, relative_velocity_m_s) in [
        (PlanarForcedConvectionRegime::Laminar, [1.0, 12.0, 0.0]),
        (PlanarForcedConvectionRegime::Laminar, [0.0; 3]),
        (PlanarForcedConvectionRegime::Turbulent, [0.0, 0.01, 0.0]),
        (PlanarForcedConvectionRegime::Laminar, [0.0, 10000.0, 0.0]),
    ] {
        let law = PlanarFreeSurfaceConvection {
            law: PlanarFreeSurfaceHeatLaw::Forced {
                regime,
                relative_velocity_m_s,
            },
            ..forced()
        };
        assert!(law.coefficient(q).is_err());
    }
    assert!(
        forced().coefficient(query(400.0, &bulk, &bulk)).is_err(),
        "bulk properties cannot silently replace actual mean-film properties"
    );
}
