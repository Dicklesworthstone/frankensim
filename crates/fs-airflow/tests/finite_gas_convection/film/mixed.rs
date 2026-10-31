//! Mixed-convection source states feed the original gas/body/motion residual.
use super::*;
use fs_airflow::graph::thermal::medium::convection::AidingVerticalPlateSettings;
use fs_airflow::graph::thermal::medium::convection::finite_gas::film::AidingVerticalPlateGasSurface;

fn patch(speed: f64) -> AidingVerticalPlateGasSurface {
    AidingVerticalPlateGasSurface {
        area_m2: 0.02,
        characteristic_length: Length::new(0.1),
        settings: AidingVerticalPlateSettings {
            speed_up_m_s: speed,
            gravity_m_s2: 9.80665,
            maximum_relative_density_departure: 0.1,
        },
    }
}

#[test]
fn aiding_film_fan_stop_and_running_flow_preserve_gas_body_motion_energy() {
    let curve = source::body_curve();
    let old = piston(0.75, 400.0, 100000.0);
    let saved = old.clone();
    let thermal = body(&curve, 120000.0); // 420 K, below the phase transition.
    let mut last_conductance = 0.0;
    for speed in [0.0, 0.25, 1.0] {
        let film = FilmGasSurface {
            model: old.gas().model(),
            law: FilmGasHeatLaw::AidingVerticalPlate(patch(speed)),
        };
        let fixed = film.resolve(old.gas().state(), thermal.initial).unwrap();
        assert!(fixed.conductance_w_k > last_conductance);
        last_conductance = fixed.conductance_w_k;
        let sample = fixed.source.interface.evaluation.as_ref().unwrap();
        assert_eq!(
            sample.nusselt.card().id,
            CorrelationId::ManningQureshiAidingVerticalPlate
        );
        assert_eq!(sample.nusselt.groups()["Re"] == 0.0, speed == 0.0);
        let next = film
            .finite_body_candidate(&old, input(), thermal, controls(), exchange(), || false)
            .unwrap();
        let source = &next.thermal().report;
        let property = source.property_gas.as_ref().unwrap();
        let gas = &source.interface.gas;
        near(
            property.state().temperature,
            0.5 * (gas.state().temperature + source.interface.body.temperature_k()),
            1e-12,
        );
        assert_eq!(property.state().pressure, gas.state().pressure);
        assert_eq!(
            property.thermodynamics().receipt().declared_fraction_bits(),
            gas.thermodynamics().receipt().declared_fraction_bits()
        );
        let du = old.gas().mass_kg()
            * (next.piston().gas().state().specific_internal_energy_j_kg()
                - old.gas().state().specific_internal_energy_j_kg());
        let dh = thermal.mass_kg
            * (next.thermal().state.specific_enthalpy_j_kg()
                - thermal.initial.specific_enthalpy_j_kg());
        let kinetic = 0.5 * old.mass_kg() * next.piston().velocity_m_s().powi(2);
        let work = input().exterior_pressure_pa
            * (next.piston().gas().volume_m3() - old.gas().volume_m3());
        near(du + dh + kinetic + work, 0.0, 1e-7);
        near(dh + next.step().exchanged_heat_j, 0.0, 1e-7);
        assert!(next.step().exchanged_heat_j > 0.0);
    }
    assert_eq!(old, saved);
}

#[test]
fn opposing_or_excessive_density_departure_cannot_publish_a_gas_candidate() {
    let curve = source::body_curve();
    let old = piston(0.75, 400.0, 100000.0);
    let saved = old.clone();
    let thermal = body(&curve, 120000.0);
    let mut excessive = patch(0.25);
    excessive.settings.maximum_relative_density_departure = 0.01;
    for surface in [patch(-0.25), excessive] {
        let film = FilmGasSurface {
            model: old.gas().model(),
            law: FilmGasHeatLaw::AidingVerticalPlate(surface),
        };
        assert!(
            film.finite_body_candidate(&old, input(), thermal, controls(), exchange(), || false)
                .is_err()
        );
        assert_eq!(old, saved);
    }
}
