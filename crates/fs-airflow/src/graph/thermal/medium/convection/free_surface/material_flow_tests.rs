//! Original plate correlations on the stationary material/phase-flow owner.
use super::*;
use fs_conduction::transient::deposition::CellId;
use fs_lbm::conjugate::enthalpy::advection::feedback::environment::{
    FiniteGasContact, FiniteGasConvection, FiniteGasConvectionContact, FiniteGasConvectionQuery,
};

fn phase_query<'a>(
    surface: f64,
    bulk: &'a ResolvedThermochemicalGasState,
    film: &'a ResolvedThermochemicalGasState,
) -> FiniteGasConvectionQuery<'a> {
    FiniteGasConvectionQuery {
        contact: FiniteGasConvectionContact {
            contact: FiniteGasContact {
                cell: CellId(7), area_m2: 1e-6, distances_m: [0.002, 0.0],
            },
            outward_normal: [1.0, 0.0, 0.0],
        },
        condensed_temperature_k: surface + 20.0 * (surface - bulk.state().temperature).signum(),
        surface_temperature_k: surface,
        bulk,
        film,
    }
}

#[test]
fn phase_flow_plate_uses_the_solved_surface_and_original_pressure_not_cell_temperature() {
    let model = source::model_transport(0.75, true);
    let mut coefficients = Vec::new();
    for pressure in [100000.0, 400000.0] {
        for (surface, ambient) in [(480.0, 320.0), (320.0, 480.0)] {
            let bulk = model.evaluate(&source::point(ambient, pressure)).unwrap();
            let film = model.evaluate(&source::point(400.0, pressure)).unwrap();
            let q = phase_query(surface, &bulk, &film);
            let actual = FiniteGasConvection::coefficient(&forced(), q).unwrap().unwrap();
            let gas = film.state();
            let re = gas.density * 12.0 * 0.1 / gas.dynamic_viscosity;
            let pr = gas.specific_heat_cp * gas.dynamic_viscosity / gas.thermal_conductivity;
            let independent = 0.664 * re.sqrt() * pr.cbrt() * gas.thermal_conductivity / 0.1;
            assert!((actual.value.value() / independent - 1.0).abs() < 1e-12);
            assert_eq!(actual, FreeSurfaceGasConvection::coefficient(
                &forced(), query(surface, &bulk, &film),
            ).unwrap().unwrap());
            if surface > ambient { coefficients.push(actual.value.value()); }
            let cell_film = model.evaluate(&source::point(
                0.5 * q.condensed_temperature_k + 0.5 * ambient, pressure,
            )).unwrap();
            assert!(FiniteGasConvection::coefficient(&forced(), FiniteGasConvectionQuery {
                film: &cell_film, ..q
            }).is_err());
        }
    }
    assert!((coefficients[1] / coefficients[0] - 2.0).abs() < 1e-12);
}

#[test]
fn authored_phase_contact_normals_keep_natural_and_aiding_orientation_and_equilibrium() {
    let model = source::model_transport(0.75, true);
    for (surface, ambient, speed) in [(360.0, 340.0, 0.2), (340.0, 360.0, -0.2)] {
        let bulk = model.evaluate(&source::point(ambient, 100000.0)).unwrap();
        let film = model.evaluate(&source::point(350.0, 100000.0)).unwrap();
        let q = phase_query(surface, &bulk, &film);
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
            let actual = FiniteGasConvection::coefficient(&law, q).unwrap().unwrap();
            assert_eq!(actual, FreeSurfaceGasConvection::coefficient(
                &law, query(surface, &bulk, &film),
            ).unwrap().unwrap());
        }
        let mut reversed = q;
        reversed.contact.outward_normal = [-1.0, 0.0, 0.0];
        assert!(FiniteGasConvection::coefficient(&mixed, reversed).is_err());
        let opposing = PlanarFreeSurfaceConvection {
            law: PlanarFreeSurfaceHeatLaw::AidingVerticalPlate {
                plate_up: [0.0, 1.0, 0.0],
                relative_velocity_m_s: [0.0, -speed, 0.0],
                gravity_m_s2: [0.0, -9.81, 0.0],
                maximum_relative_density_departure: 0.1,
            },
            ..mixed
        };
        assert!(FiniteGasConvection::coefficient(&opposing, q).is_err());
        let equilibrium = phase_query(350.0, &film, &film);
        assert_eq!(FiniteGasConvection::coefficient(
            &natural([0.0, -9.81, 0.0]), equilibrium,
        ).unwrap(), None);
        assert!(FiniteGasConvection::coefficient(&natural([0.0; 3]), equilibrium).is_err());
    }
}
