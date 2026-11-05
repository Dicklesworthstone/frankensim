//! Real plate laws on native finite-solid faces, including signed outer walls.
use super::*;
use fs_lbm::material::free_surface::enthalpy::conduction::{
    FreeSurfaceBodyGasConvection, FreeSurfaceBodyGasFace,
};

fn body_query<'a>(
    surface: f64,
    bulk: &'a ResolvedThermochemicalGasState,
    film: &'a ResolvedThermochemicalGasState,
) -> FreeSurfaceBodyConvectionQuery<'a> {
    FreeSurfaceBodyConvectionQuery {
        body_id: 17,
        face: FreeSurfaceBodyGasFace {
            wall: [-1, 2, 2],
            gas: [0, 2, 2],
            area_m2: 1e-6,
            distances_m: [0.002, 0.0],
        },
        outward_normal: [1.0, 0.0, 0.0],
        surface_temperature_k: surface,
        bulk,
        film,
    }
}

#[test]
fn solid_plate_uses_actual_surface_film_and_original_pressure_in_both_heat_directions() {
    let model = source::model_transport(0.75, true);
    let mut coefficients = Vec::new();
    for pressure in [100000.0, 400000.0] {
        for (surface, ambient) in [(480.0, 320.0), (320.0, 480.0)] {
            let bulk = model.evaluate(&source::point(ambient, pressure)).unwrap();
            let film = model.evaluate(&source::point(400.0, pressure)).unwrap();
            let sample = body_query(surface, &bulk, &film);
            let actual = FreeSurfaceBodyGasConvection::coefficient(&forced(), sample)
                .unwrap().unwrap();
            let gas = film.state();
            let re = gas.density * 12.0 * 0.1 / gas.dynamic_viscosity;
            let pr = gas.specific_heat_cp * gas.dynamic_viscosity / gas.thermal_conductivity;
            let expected = 0.664 * re.sqrt() * pr.cbrt() * gas.thermal_conductivity / 0.1;
            assert!((actual.value.value() / expected - 1.0).abs() < 1e-12);
            // Material kind and signed native coordinates do not change the
            // original plate correlation or create a fresh entrance per face.
            let liquid = FreeSurfaceGasConvection::coefficient(
                &forced(), query(surface, &bulk, &film),
            ).unwrap().unwrap();
            assert_eq!(actual, liquid);
            if surface > ambient { coefficients.push(actual.value.value()); }
            let mut cell_temperature_film = sample;
            cell_temperature_film.film = &bulk;
            assert!(FreeSurfaceBodyGasConvection::coefficient(
                &forced(), cell_temperature_film,
            ).is_err());
        }
    }
    assert!((coefficients[1] / coefficients[0] - 2.0).abs() < 1e-12);
}

#[test]
fn solid_natural_and_aiding_plates_keep_buoyancy_sign_orientation_and_equilibrium() {
    let model = source::model_transport(0.75, true);
    for (surface, ambient, speed) in [(360.0, 340.0, 0.2), (340.0, 360.0, -0.2)] {
        let bulk = model.evaluate(&source::point(ambient, 100000.0)).unwrap();
        let film = model.evaluate(&source::point(350.0, 100000.0)).unwrap();
        let q = body_query(surface, &bulk, &film);
        let mixed = PlanarFreeSurfaceConvection {
            law: PlanarFreeSurfaceHeatLaw::AidingVerticalPlate {
                plate_up: [0.0, 1.0, 0.0],
                relative_velocity_m_s: [0.0, speed, 0.0],
                gravity_m_s2: [0.0, -9.81, 0.0],
                maximum_relative_density_departure: 0.1,
            },
            ..forced()
        };
        for law in [natural([0.0, -9.81, 0.0]), natural([-9.81, 0.0, 0.0]), mixed] {
            let actual = FreeSurfaceBodyGasConvection::coefficient(&law, q).unwrap().unwrap();
            let liquid = FreeSurfaceGasConvection::coefficient(
                &law, query(surface, &bulk, &film),
            ).unwrap().unwrap();
            assert_eq!(actual, liquid);
            assert!(actual.value.value() > 0.0);
        }
        let opposing = PlanarFreeSurfaceConvection {
            law: PlanarFreeSurfaceHeatLaw::AidingVerticalPlate {
                plate_up: [0.0, 1.0, 0.0],
                relative_velocity_m_s: [0.0, -speed, 0.0],
                gravity_m_s2: [0.0, -9.81, 0.0],
                maximum_relative_density_departure: 0.1,
            },
            ..forced()
        };
        assert!(FreeSurfaceBodyGasConvection::coefficient(&opposing, q).is_err());
        let mut reversed = q;
        reversed.outward_normal = [-1.0, 0.0, 0.0];
        assert!(FreeSurfaceBodyGasConvection::coefficient(&mixed, reversed).is_err());
        let eq = body_query(350.0, &film, &film);
        assert_eq!(FreeSurfaceBodyGasConvection::coefficient(
            &natural([0.0, -9.81, 0.0]), eq,
        ).unwrap(), None);
        assert!(FreeSurfaceBodyGasConvection::coefficient(&natural([0.0; 3]), eq).is_err());
    }
}
