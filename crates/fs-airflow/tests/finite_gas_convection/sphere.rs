//! G0/G1/G3/G4: source-backed sphere heat reaches the original finite solver.
use super::*;
use fs_airflow::graph::thermal::medium::convection::finite_gas::sphere::ForcedSphereGasSurface;
use fs_material::fluid_transport::FluidTransport;
use fs_qty::Temperature;
use fs_qty::semantic::{Composition, CompositionBasis};

fn hot_piston(fraction: f64) -> GasPiston {
    GasPiston::new(
        ClosedGasParcel::new(
            source::model_transport(fraction, true),
            source::point(600.0, 100000.0),
            0.001,
        )
        .unwrap(),
        1.0,
        0.01,
        0.0,
    )
    .unwrap()
}

fn sphere(old: &GasPiston) -> ForcedSphereGasSurface<'_> {
    ForcedSphereGasSurface {
        model: old.gas().model(),
        area_m2: 0.02,
        diameter: Length::new(0.1),
        speed_m_s: 1.0,
    }
}

#[test]
fn g1_bulk_transport_and_surface_viscosity_keep_the_actual_sources() {
    let bank = source::model_transport(0.25, true);
    let actual = bank
        .with_composition(Composition::new(CompositionBasis::MoleFraction, vec![1.0, 0.0]).unwrap())
        .unwrap();
    let gas = actual.evaluate(&source::point(600.0, 200000.0)).unwrap();
    let curve = source::body_curve();
    let wall = body(&curve, 180000.0).initial;
    let surface = ForcedSphereGasSurface {
        model: &bank,
        area_m2: 0.02,
        diameter: Length::new(0.1),
        speed_m_s: 1.0,
    };
    let result = surface.resolve(&gas, wall).unwrap();
    let at_wall = result.source.surface_gas.as_ref().unwrap();
    assert_eq!(at_wall.state().temperature, 450.0);
    assert_eq!(at_wall.state().pressure, gas.state().pressure);
    assert_eq!(at_wall.source_packs(), gas.source_packs());
    assert_eq!(at_wall.components().len(), 1);
    assert_eq!(
        at_wall.thermodynamics().receipt().declared_fraction_bits(),
        gas.thermodynamics().receipt().declared_fraction_bits()
    );
    assert_eq!(
        at_wall.components()[0].card_identity(),
        gas.components()[0].card_identity()
    );
    for (property, original) in at_wall.components()[0]
        .properties()
        .iter()
        .zip(gas.components()[0].properties())
    {
        assert_eq!(
            property.answer().receipt.selected,
            original.answer().receipt.selected
        );
    }
    let bulk = gas.transport_state().unwrap();
    let re = bulk.density_kg_m3() * surface.speed_m_s * surface.diameter.value()
        / bulk.dynamic_viscosity_pa_s();
    let ratio =
        bulk.dynamic_viscosity_pa_s() / at_wall.transport_state().unwrap().dynamic_viscosity_pa_s();
    let nu = 2.0
        + (0.4 * re.sqrt() + 0.06 * re.powf(2.0 / 3.0))
            * bulk.prandtl().powf(0.4)
            * ratio.powf(0.25);
    let evaluation = result.source.interface.evaluation.as_ref().unwrap();
    assert_eq!(evaluation.transport, bulk);
    near(evaluation.nusselt.groups()["Re"], re, 1e-10);
    near(
        evaluation.nusselt.groups()["mu_bulk_over_surface"],
        4.0 / 3.0,
        1e-14,
    );
    near(
        result.conductance_w_k,
        surface.area_m2 * nu * bulk.thermal_conductivity_w_m_k() / surface.diameter.value(),
        1e-11,
    );
    assert_eq!(result.source.interface.gas, gas);
    assert_eq!(result.source.interface.body, wall);

    let colder = surface
        .resolve(&gas, body(&curve, 50000.0).initial)
        .unwrap();
    assert!(colder.conductance_w_k > result.conductance_w_k);
    let lower_pressure = actual.evaluate(&source::point(600.0, 100000.0)).unwrap();
    assert!(
        surface
            .resolve(&lower_pressure, wall)
            .unwrap()
            .conductance_w_k
            < result.conductance_w_k
    );
}

#[test]
fn g1_sphere_heat_melts_the_body_and_closes_gas_motion_energy() {
    let curve = source::body_curve();
    for fraction in [0.0, 1.0] {
        let old = hot_piston(fraction);
        let saved = old.clone();
        let thermal = body(&curve, 180000.0);
        let surface = sphere(&old);
        let next = surface
            .finite_body_candidate(&old, input(), thermal, controls(), exchange(), || false)
            .unwrap();
        let sample = &next.thermal().report;
        let report = next.step();
        assert!(report.exchanged_heat_j < 0.0);
        assert!(next.thermal().state.liquid_mass_fraction() > 0.5);
        near(next.thermal().state.temperature_k(), 450.0, 1e-12);
        let bulk = &sample.interface.gas;
        let at_wall = sample.surface_gas.as_ref().unwrap();
        assert_eq!(bulk, &report.piston.midpoint);
        near(
            at_wall.state().temperature,
            report.thermal_midpoint.temperature_k(),
            1e-12,
        );
        assert_eq!(at_wall.state().pressure, bulk.state().pressure);
        let evaluation = sample.interface.evaluation.as_ref().unwrap();
        assert_eq!(evaluation.transport, bulk.transport_state().unwrap());
        near(
            report.exchanged_heat_j,
            DT * report.conductance_w_k
                * (sample.interface.body.temperature_k() - bulk.state().temperature),
            exchange().source_tolerance_j,
        );
        let gas_energy = old.gas().mass_kg()
            * (next.piston().gas().state().specific_internal_energy_j_kg()
                - old.gas().state().specific_internal_energy_j_kg());
        let body_energy = thermal.mass_kg
            * (next.thermal().state.specific_enthalpy_j_kg()
                - thermal.initial.specific_enthalpy_j_kg());
        let kinetic = 0.5 * old.mass_kg() * next.piston().velocity_m_s().powi(2);
        let exterior = input().exterior_pressure_pa
            * (next.piston().gas().volume_m3() - old.gas().volume_m3());
        near(gas_energy + body_energy + kinetic + exterior, 0.0, 1e-7);
        near(body_energy + report.exchanged_heat_j, 0.0, 1e-7);
        let replay = surface
            .finite_body_candidate(&old, input(), thermal, controls(), exchange(), || false)
            .unwrap();
        assert_eq!(replay.step(), report);
        assert_eq!(replay.thermal().state, next.thermal().state);
        assert_eq!(old, saved);
    }
}

#[test]
fn g0_isolation_and_source_refusals_do_not_invent_a_second_heat_model() {
    let old = hot_piston(1.0);
    let curve = source::body_curve();
    let thermal = body(&curve, 180000.0);
    let surface = sphere(&old);
    let isolated = ForcedSphereGasSurface {
        area_m2: 0.0,
        speed_m_s: 0.0,
        ..surface
    };
    let next = isolated
        .finite_body_candidate(&old, input(), thermal, controls(), exchange(), || false)
        .unwrap();
    let closed = old
        .thermal_candidate(
            input(),
            thermal,
            controls(),
            exchange(),
            |_, _| {
                Ok::<_, std::convert::Infallible>(ConductanceSample {
                    conductance_w_k: 0.0,
                    source: (),
                })
            },
            || false,
        )
        .unwrap();
    assert_eq!(next.piston(), closed.piston());
    assert_eq!(next.thermal().state, closed.thermal().state);
    assert_eq!(next.step(), closed.step());
    assert!(next.thermal().report.surface_gas.is_none());
    assert!(next.thermal().report.interface.evaluation.is_none());

    let cold = old
        .gas()
        .model()
        .evaluate_temperature_from_snapshot(old.gas().state(), Temperature::new(400.0))
        .unwrap();
    // A colder ordinary gas has mu_bulk/mu_surface < 1; no unit-ratio clamp.
    assert!(surface.resolve(&cold, thermal.initial).is_err());
    assert!(
        ForcedSphereGasSurface {
            speed_m_s: 0.0,
            ..surface
        }
        .resolve(old.gas().state(), thermal.initial)
        .is_err()
    );
    assert!(
        ForcedSphereGasSurface {
            area_m2: 0.04,
            ..surface
        }
        .resolve(old.gas().state(), thermal.initial)
        .is_err()
    );
    assert!(
        ForcedSphereGasSurface {
            area_m2: 0.025,
            ..surface
        }
        .finite_body_candidate(&old, input(), thermal, controls(), exchange(), || false)
        .is_err()
    );
}

#[test]
fn g4_cancelled_sphere_heat_leaves_both_original_states_unchanged() {
    let old = hot_piston(1.0);
    let saved = old.clone();
    let curve = source::body_curve();
    let thermal = body(&curve, 180000.0);
    let surface = sphere(&old);
    let mut calls = 0;
    surface
        .finite_body_candidate(&old, input(), thermal, controls(), exchange(), || {
            calls += 1;
            false
        })
        .unwrap();
    for stop in [1, calls / 2, calls] {
        let mut count = 0;
        assert!(matches!(
            surface.finite_body_candidate(&old, input(), thermal, controls(), exchange(), || {
                count += 1;
                count == stop
            }),
            Err(ThermalPistonError::Cancelled)
        ));
        assert_eq!(old, saved);
        assert_eq!(
            thermal.initial,
            curve.state_at_specific_enthalpy(180000.0).unwrap()
        );
    }
}
