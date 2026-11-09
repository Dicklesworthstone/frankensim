//! Source/coupling invariants for an engineering closure, not gas validation.
use super::*;
use fs_airflow::graph::thermal::medium::convection::finite_gas::tube::ForcedTubeGasSurface;
use fs_material::fluid_transport::FluidTransport;

fn gas(temperature_k: f64) -> GasPiston {
    GasPiston::new(
        ClosedGasParcel::new(
            source::model_transport(1.0, true),
            source::point(temperature_k, 100000.0),
            0.001,
        )
        .unwrap(),
        1.0,
        0.01,
        0.0,
    )
    .unwrap()
}

fn tube(old: &GasPiston) -> ForcedTubeGasSurface<'_> {
    let bulk = old.gas().state().transport_state().unwrap();
    ForcedTubeGasSurface {
        model: old.gas().model(),
        area_m2: 0.02,
        diameter: Length::new(0.01),
        heated_length: Length::new(1.0),
        speed_m_s: 20000.0 * bulk.dynamic_viscosity_pa_s() / (bulk.density_kg_m3() * 0.01),
    }
}

#[test]
fn wall_viscosity_changes_transfer_in_both_directions_without_replacing_bulk_properties() {
    let old = gas(400.0);
    let bulk = old.gas().state();
    // An unrelated initial composition cannot replace the actual pure species.
    let bank = source::model_transport(0.25, true);
    let surface = ForcedTubeGasSurface {
        model: &bank,
        ..tube(&old)
    };
    let curve = source::body_curve();
    let transport = bulk.transport_state().unwrap();
    let mut previous = f64::INFINITY;
    for enthalpy in [80000.0, 100000.0, 120000.0] {
        let wall = body(&curve, enthalpy).initial;
        let result = surface.resolve(bulk, wall).unwrap();
        let at_wall = result.source.surface_gas.as_ref().unwrap();
        let evaluated = result.source.interface.evaluation.as_ref().unwrap();
        assert_eq!(result.source.interface.gas, *bulk);
        assert_eq!(result.source.interface.body, wall);
        assert_eq!(evaluated.transport, transport);
        assert_eq!(
            evaluated.nusselt.card().id,
            CorrelationId::SiederTateTurbulent
        );
        assert_eq!(at_wall.state().temperature, wall.temperature_k());
        assert_eq!(at_wall.state().pressure, bulk.state().pressure);
        assert_eq!(at_wall.source_packs(), bulk.source_packs());
        assert_eq!(at_wall.components().len(), 1);
        assert_eq!(
            at_wall.thermodynamics().receipt().mole_fraction_bits(),
            bulk.thermodynamics().receipt().mole_fraction_bits()
        );
        assert_eq!(
            at_wall.components()[0].card_identity(),
            bulk.components()[0].card_identity()
        );
        for (query, original) in at_wall.components()[0]
            .properties()
            .iter()
            .zip(bulk.components()[0].properties())
        {
            assert_eq!(
                query.answer().receipt.selected,
                original.answer().receipt.selected
            );
        }
        let ratio = 400.0 / wall.temperature_k();
        let re = transport.density_kg_m3() * surface.speed_m_s * surface.diameter.value()
            / transport.dynamic_viscosity_pa_s();
        let nu = 0.027 * re.powf(0.8) * transport.prandtl().powf(1.0 / 3.0) * ratio.powf(0.14);
        near(evaluated.nusselt.groups()["Re"], 20000.0, 1e-10);
        assert_eq!(evaluated.nusselt.groups()["L_over_Dh"], 100.0);
        near(
            evaluated.nusselt.groups()["mu_bulk_over_surface"],
            ratio,
            1e-13,
        );
        near(evaluated.nusselt.evidence().value, nu, 1e-11);
        near(
            result.conductance_w_k,
            surface.area_m2 * nu * transport.thermal_conductivity_w_m_k()
                / surface.diameter.value(),
            1e-12,
        );
        assert!(result.conductance_w_k < previous);
        previous = result.conductance_w_k;
        if enthalpy == 100000.0 {
            // Equal temperatures remove the driving force, not the conductance.
            assert_eq!(at_wall.transport_state().unwrap(), transport);
            assert!(result.conductance_w_k > 0.0);
        } else {
            assert_ne!(
                at_wall
                    .transport_state()
                    .unwrap()
                    .thermal_conductivity_w_m_k(),
                transport.thermal_conductivity_w_m_k()
            );
        }
    }
}

#[test]
fn tube_heat_closes_piston_and_body_energy_including_a_latent_plateau() {
    let curve = source::body_curve();
    for (temperature, enthalpy) in [(400.0, 80000.0), (400.0, 120000.0), (470.0, 180000.0)] {
        let old = gas(temperature);
        let saved = old.clone();
        let thermal = body(&curve, enthalpy);
        let surface = tube(&old);
        let next = surface
            .finite_body_candidate(&old, input(), thermal, controls(), exchange(), || false)
            .unwrap();
        let sample = &next.thermal().report;
        let report = next.step();
        assert_eq!(sample.interface.gas, report.piston.midpoint);
        assert_eq!(sample.interface.body, report.thermal_midpoint);
        near(
            sample.surface_gas.as_ref().unwrap().state().temperature,
            report.thermal_midpoint.temperature_k(),
            1e-12,
        );
        assert_eq!(
            report.exchanged_heat_j > 0.0,
            thermal.initial.temperature_k() > temperature
        );
        near(
            report.exchanged_heat_j,
            DT * report.conductance_w_k
                * (sample.interface.body.temperature_k()
                    - sample.interface.gas.state().temperature),
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
        if enthalpy == 180000.0 {
            assert_eq!(next.thermal().state.temperature_k(), 450.0);
            assert!(next.thermal().state.liquid_mass_fraction() > 0.5);
        }
        assert_eq!(old, saved);
    }
}

#[test]
fn source_and_tube_admission_refuse_while_explicit_insulation_stays_exact() {
    let old = gas(400.0);
    let saved = old.clone();
    let curve = source::body_curve();
    let thermal = body(&curve, 100000.0);
    let surface = tube(&old);
    // The wrong bank is numerically identical at 400 K but has different claims.
    let wrong_bank = source::model_transport(1.0, false);
    assert!(
        ForcedTubeGasSurface {
            model: &wrong_bank,
            ..surface
        }
        .resolve(old.gas().state(), thermal.initial)
        .is_err()
    );
    for invalid in [
        ForcedTubeGasSurface {
            speed_m_s: 0.0,
            ..surface
        },
        ForcedTubeGasSurface {
            speed_m_s: surface.speed_m_s / 10.0,
            ..surface
        },
        ForcedTubeGasSurface {
            heated_length: Length::new(0.5),
            area_m2: 0.01,
            ..surface
        },
        ForcedTubeGasSurface {
            area_m2: 0.04,
            ..surface
        },
    ] {
        assert!(invalid.resolve(old.gas().state(), thermal.initial).is_err());
    }
    assert!(
        ForcedTubeGasSurface {
            area_m2: 0.025,
            ..surface
        }
        .finite_body_candidate(&old, input(), thermal, controls(), exchange(), || false)
        .is_err()
    );
    let insulated = ForcedTubeGasSurface {
        area_m2: 0.0,
        speed_m_s: 0.0,
        ..surface
    }
    .finite_body_candidate(&old, input(), thermal, controls(), exchange(), || false)
    .unwrap();
    assert_eq!(insulated.step().exchanged_heat_j, 0.0);
    assert_eq!(insulated.thermal().state, thermal.initial);
    assert!(insulated.thermal().report.surface_gas.is_none());
    assert!(insulated.thermal().report.interface.evaluation.is_none());
    assert!(matches!(
        surface.finite_body_candidate(&old, input(), thermal, controls(), exchange(), || true),
        Err(ThermalPistonError::Cancelled)
    ));
    assert_eq!(old, saved);
}

#[path = "tube/laminar.rs"]
mod laminar;
