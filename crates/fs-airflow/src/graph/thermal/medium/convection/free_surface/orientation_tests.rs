use super::*;

// A proper orthogonal rotation with rational entries and no Cartesian column.
// Applying it twice exercises a second unrelated physical orientation.
fn rotate(v: [f64; 3]) -> [f64; 3] {
    [
        (-10.0 * v[0] + 2.0 * v[1] + 11.0 * v[2]) / 15.0,
        (2.0 * v[0] - v[1] + 2.0 * v[2]) / 3.0,
        (5.0 * v[0] + 14.0 * v[1] + 2.0 * v[2]) / 15.0,
    ]
}

fn rotated(mut plate: PlanarFreeSurfaceConvection) -> PlanarFreeSurfaceConvection {
    plate.outward_normal = rotate(plate.outward_normal);
    plate.law = match plate.law {
        PlanarFreeSurfaceHeatLaw::Forced {
            regime,
            relative_velocity_m_s,
        } => PlanarFreeSurfaceHeatLaw::Forced {
            regime,
            relative_velocity_m_s: rotate(relative_velocity_m_s),
        },
        PlanarFreeSurfaceHeatLaw::Natural { gravity_m_s2 } => PlanarFreeSurfaceHeatLaw::Natural {
            gravity_m_s2: rotate(gravity_m_s2),
        },
        PlanarFreeSurfaceHeatLaw::AidingVerticalPlate {
            plate_up,
            relative_velocity_m_s,
            gravity_m_s2,
            maximum_relative_density_departure,
        } => PlanarFreeSurfaceHeatLaw::AidingVerticalPlate {
            plate_up: rotate(plate_up),
            relative_velocity_m_s: rotate(relative_velocity_m_s),
            gravity_m_s2: rotate(gravity_m_s2),
            maximum_relative_density_departure,
        },
    };
    plate
}

fn aiding(speed: f64, gravity: f64) -> PlanarFreeSurfaceConvection {
    PlanarFreeSurfaceConvection {
        law: PlanarFreeSurfaceHeatLaw::AidingVerticalPlate {
            plate_up: [0.0, 3.0, 0.0],
            relative_velocity_m_s: [0.0, speed, 0.0],
            gravity_m_s2: [0.0, -gravity, 0.0],
            maximum_relative_density_departure: 0.1,
        },
        ..forced()
    }
}

fn check_rotations(plate: PlanarFreeSurfaceConvection, q: FreeSurfaceConvectionQuery<'_>) {
    let expected = plate.coefficient(q).unwrap();
    let mut plate = plate;
    let mut q = q;
    for _ in 0..2 {
        plate = rotated(plate);
        // The interface owner supplies the physical normal, independently of
        // the Cartesian indices of cells containing a reconstructed interface.
        q.outward_normal = rotate(q.outward_normal);
        let actual = plate.coefficient(q).unwrap();
        match (&expected, actual) {
            (Some(expected), Some(actual)) => {
                assert!((actual.value.value() / expected.value.value() - 1.0).abs() < 1e-12);
            }
            (None, None) => {}
            _ => panic!("rotation changed equilibrium or correlation availability"),
        }
    }
}

#[test]
fn forced_plate_is_world_rotation_invariant_at_both_admitted_regimes() {
    let model = source::model_transport(1.0, true);
    for (regime, speed, pressure) in [
        (PlanarForcedConvectionRegime::Laminar, 12.0, 100000.0),
        (PlanarForcedConvectionRegime::Turbulent, 30.0, 1000000.0),
    ] {
        let bulk = model.evaluate(&source::point(300.0, pressure)).unwrap();
        let film = model.evaluate(&source::point(350.0, pressure)).unwrap();
        let plate = PlanarFreeSurfaceConvection {
            law: PlanarFreeSurfaceHeatLaw::Forced {
                regime,
                relative_velocity_m_s: [0.0, speed, 0.0],
            },
            ..forced()
        };
        check_rotations(plate, query(400.0, &bulk, &film));
    }
}

#[test]
fn natural_vertical_and_horizontal_branches_follow_rotated_gravity_in_both_heat_directions() {
    let model = source::model(0.75);
    let film = model.evaluate(&source::point(350.0, 100000.0)).unwrap();
    for (wall, ambient) in [(400.0, 300.0), (300.0, 400.0), (350.0, 350.0)] {
        let bulk = model.evaluate(&source::point(ambient, 100000.0)).unwrap();
        for gravity in [[0.0, -9.81, 0.0], [-9.81, 0.0, 0.0], [9.81, 0.0, 0.0]] {
            check_rotations(natural(gravity), query(wall, &bulk, &film));
        }
    }
}

#[test]
fn aiding_plate_rotates_with_actual_relative_flow_including_pure_limits_and_equilibrium() {
    let model = source::model_transport(0.75, true);
    let film = model.evaluate(&source::point(350.0, 100000.0)).unwrap();
    for (wall, ambient, speed, gravity) in [
        (360.0, 340.0, 0.2, 9.81),
        (340.0, 360.0, -0.2, 9.81),
        (360.0, 340.0, 0.0, 9.81),
        (360.0, 340.0, -0.2, 0.0),
        (350.0, 350.0, 0.0, 9.81),
        (350.0, 350.0, 0.0, 0.0),
        (350.0, 350.0, 0.2, 9.81),
    ] {
        let bulk = model.evaluate(&source::point(ambient, 100000.0)).unwrap();
        check_rotations(aiding(speed, gravity), query(wall, &bulk, &film));
    }
}

#[test]
fn rotated_normal_alignment_tangency_and_inclination_admit_only_roundoff() {
    let model = source::model(0.75);
    let bulk = model.evaluate(&source::point(340.0, 100000.0)).unwrap();
    let film = model.evaluate(&source::point(350.0, 100000.0)).unwrap();
    let plate = rotated(forced());
    let mut q = query(360.0, &bulk, &film);
    q.outward_normal = plate.outward_normal;
    let expected = plate.coefficient(q).unwrap().unwrap();
    // Small normalization roundoff is admitted; a material angular tolerance
    // is not. These perturbations remain far inside/outside the declared bound.
    let rounded_plate = PlanarFreeSurfaceConvection {
        outward_normal: plate.outward_normal.map(|x| x * (1.0 + 8.0 * f64::EPSILON)),
        ..plate
    };
    let mut rounded_query = q;
    rounded_query.outward_normal = q.outward_normal.map(|x| x * (1.0 - 8.0 * f64::EPSILON));
    assert_eq!(
        rounded_plate.coefficient(rounded_query).unwrap().unwrap(),
        expected
    );
    for normal in [
        [0.0; 3],
        [f64::NAN, 0.0, 0.0],
        [f64::INFINITY, 0.0, 0.0],
        q.outward_normal.map(|x| 2.0 * x),
        q.outward_normal.map(|x| x * (1.0 + 128.0 * f64::EPSILON)),
    ] {
        let invalid_plate = PlanarFreeSurfaceConvection {
            outward_normal: normal,
            ..plate
        };
        assert!(invalid_plate.coefficient(q).is_err());
        assert!(
            plate
                .coefficient(FreeSurfaceConvectionQuery {
                    outward_normal: normal,
                    ..q
                })
                .is_err()
        );
    }
    // A 1e-8-radian mismatch can have a dot product rounded to one. It must
    // still refuse, as must the opposite liquid-to-gas orientation.
    for normal in [rotate([1.0, 1e-8, 0.0]), q.outward_normal.map(|x| -x)] {
        assert!(
            plate
                .coefficient(FreeSurfaceConvectionQuery {
                    outward_normal: normal,
                    ..q
                })
                .is_err()
        );
    }
    for provider in [
        PlanarFreeSurfaceConvection {
            law: PlanarFreeSurfaceHeatLaw::Forced {
                regime: PlanarForcedConvectionRegime::Laminar,
                relative_velocity_m_s: [0.001, 12.0, 0.0],
            },
            ..forced()
        },
        natural([1e-8, -9.81, 0.0]),
        natural([-9.81, 1e-8, 0.0]),
        aiding(-0.2, 9.81),
        PlanarFreeSurfaceConvection {
            law: PlanarFreeSurfaceHeatLaw::AidingVerticalPlate {
                plate_up: [0.0, 3.0, 0.0],
                relative_velocity_m_s: [0.0, 0.2, 0.001],
                gravity_m_s2: [0.0, -9.81, 0.0],
                maximum_relative_density_departure: 0.1,
            },
            ..forced()
        },
    ] {
        assert!(
            rotated(provider).coefficient(q).is_err(),
            "unsupported orientation admitted: {provider:?}"
        );
    }
    let mut equilibrium = query(350.0, &film, &film);
    equilibrium.outward_normal = q.outward_normal;
    assert!(
        rotated(natural([-9.81, 1e-8, 0.0]))
            .coefficient(equilibrium)
            .is_err()
    );
    let tilted_rest = PlanarFreeSurfaceConvection {
        law: PlanarFreeSurfaceHeatLaw::AidingVerticalPlate {
            plate_up: [1e-8, 3.0, 0.0],
            relative_velocity_m_s: [0.0; 3],
            gravity_m_s2: [0.0; 3],
            maximum_relative_density_departure: 0.1,
        },
        ..forced()
    };
    assert!(rotated(tilted_rest).coefficient(equilibrium).is_err());
}
