use super::*;

fn law(
    plate_up: [f64; 3],
    relative_velocity_m_s: [f64; 3],
    gravity_m_s2: [f64; 3],
    maximum_relative_density_departure: f64,
) -> PlanarFreeSurfaceConvection {
    PlanarFreeSurfaceConvection {
        law: PlanarFreeSurfaceHeatLaw::AidingVerticalPlate {
            plate_up,
            relative_velocity_m_s,
            gravity_m_s2,
            maximum_relative_density_departure,
        },
        ..forced()
    }
}

fn mixed(speed_up_m_s: f64) -> PlanarFreeSurfaceConvection {
    law(
        [0.0, 3.0, 0.0],
        [0.0, speed_up_m_s, 0.0],
        [0.0, -9.81, 0.0],
        0.1,
    )
}

// Independent dimensional algebra; the original correlation still owns Nu.
fn expected(
    film: &ResolvedThermochemicalGasState,
    wall: f64,
    bulk: f64,
    speed: f64,
    gravity: f64,
) -> f64 {
    let gas = film.state();
    let beta = film
        .thermal_transport_state()
        .unwrap()
        .volumetric_expansion_per_k()
        .unwrap();
    let length = 0.1_f64;
    let re = gas.density * speed.abs() * length / gas.dynamic_viscosity;
    let ra = gravity
        * beta
        * (wall - bulk).abs()
        * length.powi(3)
        * gas.density.powi(2)
        * gas.specific_heat_cp
        / (gas.dynamic_viscosity * gas.thermal_conductivity);
    let pr = gas.dynamic_viscosity * gas.specific_heat_cp / gas.thermal_conductivity;
    fs_convection::evaluate(
        CorrelationId::ManningQureshiAidingVerticalPlate,
        CorrelationInputs::mixed(re, ra, pr),
    )
    .unwrap()
    .evidence()
    .value
        * gas.thermal_conductivity
        / length
}

#[test]
fn mixed_plate_uses_actual_film_full_contrast_and_signed_aiding_flow() {
    let model = source::model_transport(0.75, true);
    let mut pressure_coefficients = Vec::new();
    for pressure in [100000.0, 400000.0] {
        for (wall, ambient, speed) in [(360.0, 340.0, 0.2), (340.0, 360.0, -0.2)] {
            let bulk = model.evaluate(&source::point(ambient, pressure)).unwrap();
            let film = model.evaluate(&source::point(350.0, pressure)).unwrap();
            let actual = mixed(speed)
                .coefficient(query(wall, &bulk, &film))
                .unwrap()
                .unwrap();
            let direct = expected(&film, wall, ambient, speed, 9.81);
            assert!((actual.value.value() / direct - 1.0).abs() < 1e-12);
            assert_eq!(
                (actual.value.value() * (ambient - wall)).is_sign_positive(),
                ambient > wall
            );
            assert!(
                mixed(-speed)
                    .coefficient(query(wall, &bulk, &film))
                    .is_err()
            );
            if wall > ambient {
                pressure_coefficients.push(actual.value.value());
            }
        }
    }
    assert!(pressure_coefficients[1] > pressure_coefficients[0]);
}

#[test]
fn mixed_plate_retains_original_pure_limits_and_exact_resting_equilibrium() {
    let model = source::model(0.75);
    let bulk = model.evaluate(&source::point(340.0, 100000.0)).unwrap();
    let film = model.evaluate(&source::point(350.0, 100000.0)).unwrap();
    let q = query(360.0, &bulk, &film);
    for (provider, speed, gravity) in [
        (mixed(0.0), 0.0, 9.81),
        (
            law([0.0, 1.0, 0.0], [0.0, -0.2, 0.0], [0.0; 3], 0.1),
            -0.2,
            0.0,
        ),
    ] {
        let actual = provider.coefficient(q).unwrap().unwrap();
        assert!(
            (actual.value.value() / expected(&film, 360.0, 340.0, speed, gravity) - 1.0).abs()
                < 1e-12
        );
    }
    let equilibrium = query(350.0, &film, &film);
    assert_eq!(mixed(0.0).coefficient(equilibrium).unwrap(), None);
    let no_driving = law([0.0, 1.0, 0.0], [0.0; 3], [0.0; 3], 0.1);
    assert_eq!(no_driving.coefficient(equilibrium).unwrap(), None);
    assert!(
        no_driving.coefficient(q).is_err(),
        "no stagnant-gas coefficient is invented"
    );
    let flowing_equilibrium = mixed(0.2).coefficient(equilibrium).unwrap().unwrap();
    assert!(
        (flowing_equilibrium.value.value() / expected(&film, 350.0, 350.0, 0.2, 9.81) - 1.0).abs()
            < 1e-12
    );
}

#[test]
fn mixed_plate_frame_is_tangent_to_the_actual_face_and_preserves_world_orientation() {
    let model = source::model(0.75);
    let bulk = model.evaluate(&source::point(340.0, 100000.0)).unwrap();
    let film = model.evaluate(&source::point(350.0, 100000.0)).unwrap();
    let q = query(360.0, &bulk, &film);
    let expected_h = mixed(0.2).coefficient(q).unwrap().unwrap().value.value();
    // A non-Cartesian vertical direction is valid within a Cartesian face.
    let sqrt_two = 2.0_f64.sqrt();
    for provider in [
        law([0.0, -3.0, 0.0], [0.0, -0.2, 0.0], [0.0, 9.81, 0.0], 0.1),
        law(
            [0.0, 1.0, 1.0],
            [0.0, 0.2 / sqrt_two, 0.2 / sqrt_two],
            [0.0, -9.81 / sqrt_two, -9.81 / sqrt_two],
            0.1,
        ),
    ] {
        let actual = provider.coefficient(q).unwrap().unwrap();
        assert!((actual.value.value() / expected_h - 1.0).abs() < 1e-12);
    }
    for provider in [
        // Aligned gravity/flow cannot make an axis normal to the face vertical.
        law([1.0, 0.0, 0.0], [0.2, 0.0, 0.0], [-9.81, 0.0, 0.0], 0.1),
        law([0.0, 1.0, 0.0], [0.01, 0.2, 0.0], [0.0, -9.81, 0.0], 0.1),
        law([0.0, 1.0, 0.0], [0.0, 0.2, 0.01], [0.0, -9.81, 0.0], 0.1),
        law([0.0, 1.0, 0.0], [0.0, 0.2, 0.0], [0.0, -9.81, 0.01], 0.1),
    ] {
        assert!(
            provider.coefficient(q).is_err(),
            "invalid frame admitted: {provider:?}"
        );
    }
    let mut other_face = q;
    other_face.outward_normal = [0.0, 0.0, 1.0];
    other_face.surface.face.gas = [3, 3, 4];
    assert!(mixed(0.2).coefficient(other_face).is_err());
    let rotated_patch = PlanarFreeSurfaceConvection {
        outward_normal: other_face.outward_normal,
        ..mixed(0.2)
    };
    assert_eq!(
        rotated_patch
            .coefficient(other_face)
            .unwrap()
            .unwrap()
            .value
            .value(),
        expected_h
    );
}

#[test]
fn mixed_plate_never_bypasses_density_frame_source_or_original_regime_admission() {
    let model = source::model_transport(0.75, true);
    let bulk = model.evaluate(&source::point(340.0, 100000.0)).unwrap();
    let film = model.evaluate(&source::point(350.0, 100000.0)).unwrap();
    let q = query(360.0, &bulk, &film);
    let beta = film
        .thermal_transport_state()
        .unwrap()
        .volumetric_expansion_per_k()
        .unwrap();
    let departure = beta * 20.0;
    assert!(
        law(
            [0.0, 1.0, 0.0],
            [0.0, 0.2, 0.0],
            [0.0, -9.81, 0.0],
            departure
        )
        .coefficient(q)
        .is_ok()
    );
    for (velocity, gravity) in [
        ([0.0, 0.2, 0.0], [0.0, -9.81, 0.0]),
        ([0.0; 3], [0.0, -9.81, 0.0]),
        ([0.0, 0.2, 0.0], [0.0; 3]),
    ] {
        assert!(
            law([0.0, 1.0, 0.0], velocity, gravity, departure.next_down())
                .coefficient(q)
                .is_err()
        );
    }
    assert!(
        mixed(10000.0).coefficient(q).is_err(),
        "original laminar Re bound remains active"
    );
    assert!(
        PlanarFreeSurfaceConvection {
            characteristic_length: Length::new(100.0),
            ..mixed(0.0)
        }
        .coefficient(q)
        .is_err(),
        "original laminar Ra bound remains active"
    );
    assert!(mixed(0.2).coefficient(query(360.0, &bulk, &bulk)).is_err());

    let equilibrium = query(350.0, &film, &film);
    for provider in [
        law([0.0; 3], [0.0; 3], [0.0; 3], 0.1),
        law([1.0, 0.0, 0.0], [0.0; 3], [0.0; 3], 0.1),
        law([0.0, 1.0, 0.0], [0.0; 3], [0.0, -9.81, 0.1], 0.1),
        law([0.0, 1.0, 0.0], [0.0; 3], [0.0, -9.81, 0.0], f64::NAN),
        law([0.0, 1.0, 0.0], [0.0; 3], [0.0, -9.81, 0.0], 1.0),
        PlanarFreeSurfaceConvection {
            characteristic_length: Length::new(0.0),
            ..mixed(0.0)
        },
    ] {
        assert!(
            provider.coefficient(equilibrium).is_err(),
            "equilibrium bypassed {provider:?}"
        );
    }
}
