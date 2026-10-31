//! Material-driven vertical convection reaches the original finite heat/motion owner.
use super::{DT, body, controls, exchange, input, near, piston, source};
use fs_airflow::graph::thermal::medium::convection::finite_gas::NaturalGasSurface;
use fs_couple::cavity_phs::piston::thermal::{ConductanceSample, ThermalPistonError};
use fs_material::fluid_transport::FluidTransport;
use fs_qty::Length;

fn surface() -> NaturalGasSurface {
    NaturalGasSurface {
        geometry: fs_convection::NaturalConvectionGeometry::VerticalPlate,
        area_m2: 0.02,
        characteristic_length: Length::new(0.1),
        gravity_m_s2: 9.80665,
    }
}

#[test]
fn actual_species_pressure_and_transport_supply_buoyancy_and_heat_conductance() {
    let curve = source::body_curve();
    let hot = body(&curve, 260000.0).initial; // 500 K
    let cold = body(&curve, 0.0).initial; // 300 K, same contrast magnitude
    let old = piston(0.75, 400.0, 100000.0);
    let gas = old.gas().state();
    let patch = surface();
    let resolved = patch.resolve(gas, hot).unwrap();
    let evaluation = resolved.source.evaluation.as_ref().unwrap();
    let transport = gas.transport_state().unwrap();
    // Independently reconstruct Ra from rho, cp, mu and k. The material
    // evaluator supplies both dimensionless groups and dimensional h.
    let rho = transport.density_kg_m3();
    let cp = transport.specific_heat_j_kg_k();
    let mu = transport.dynamic_viscosity_pa_s();
    let k = transport.thermal_conductivity_w_m_k();
    let pr = mu * cp / k;
    let ra = patch.gravity_m_s2 * (hot.temperature_k() - 400.0) / 400.0
        * patch.characteristic_length.value().powi(3)
        * rho
        * rho
        * cp
        / (mu * k);
    let nu = (0.825
        + 0.387 * ra.powf(1.0 / 6.0) / (1.0 + (0.492 / pr).powf(9.0 / 16.0)).powf(8.0 / 27.0))
    .powi(2);
    near(evaluation.nusselt.groups()["Ra"], ra, ra * 1e-12);
    near(evaluation.nusselt.groups()["Pr"], pr, pr * 1e-12);
    near(
        resolved.conductance_w_k,
        patch.area_m2 * k * nu / patch.characteristic_length.value(),
        1e-11,
    );
    assert_eq!(evaluation.transport, transport);
    assert_eq!(resolved.source.gas, *gas);
    assert_eq!(
        resolved.source.gas.source_packs(),
        old.gas().model().source_packs()
    );
    assert_eq!(resolved.source.body, hot);
    // The vertical-plate law has the same positive h for equal hot/cold
    // contrasts. The reciprocal solver, not the coefficient, owns heat sign.
    near(
        patch.resolve(gas, cold).unwrap().conductance_w_k,
        resolved.conductance_w_k,
        1e-13,
    );
    let high = piston(0.75, 400.0, 200000.0);
    let dense = patch.resolve(high.gas().state(), hot).unwrap();
    near(
        dense.source.evaluation.as_ref().unwrap().nusselt.groups()["Ra"] / ra,
        4.0,
        1e-11,
    );
    assert!(dense.conductance_w_k > resolved.conductance_w_k);
    let other = piston(0.25, 400.0, 100000.0);
    let substituted = patch.resolve(other.gas().state(), hot).unwrap();
    assert!((substituted.conductance_w_k - resolved.conductance_w_k).abs() > 1e-3);
}

#[test]
fn natural_heat_uses_actual_midpoints_and_closes_body_gas_motion_and_latent_storage() {
    let curve = source::body_curve();
    for (gas_k, enthalpy, heating_gas) in [(400.0, 260000.0, true), (500.0, 180000.0, false)] {
        let old = piston(0.75, gas_k, 100000.0);
        let saved = old.clone();
        let thermal = body(&curve, enthalpy);
        let next = surface()
            .finite_body_candidate(&old, input(), thermal, controls(), exchange(), || false)
            .unwrap();
        let report = next.step();
        let sample = &next.thermal().report;
        assert_eq!(report.exchanged_heat_j > 0.0, heating_gas);
        assert_eq!(sample.gas, report.piston.midpoint);
        assert_eq!(sample.body, report.thermal_midpoint);
        let independent = surface()
            .resolve(&report.piston.midpoint, report.thermal_midpoint)
            .unwrap();
        near(report.conductance_w_k, independent.conductance_w_k, 1e-14);
        near(
            report.exchanged_heat_j,
            DT * report.conductance_w_k
                * (sample.body.temperature_k() - sample.gas.state().temperature),
            exchange().source_tolerance_j,
        );
        let du = old.gas().mass_kg()
            * (next.piston().gas().state().specific_internal_energy_j_kg()
                - old.gas().state().specific_internal_energy_j_kg());
        let dh = thermal.mass_kg
            * (next.thermal().state.specific_enthalpy_j_kg()
                - thermal.initial.specific_enthalpy_j_kg());
        let dk = 0.5 * old.mass_kg() * next.piston().velocity_m_s().powi(2);
        let exterior_work = input().exterior_pressure_pa
            * (next.piston().gas().volume_m3() - old.gas().volume_m3());
        near(du + dh + dk + exterior_work, 0.0, 1e-7);
        near(dh + report.exchanged_heat_j, 0.0, 1e-7);
        if !heating_gas {
            near(next.thermal().state.temperature_k(), 450.0, 1e-12);
            assert!(next.thermal().state.liquid_mass_fraction() > 0.5);
            near(
                thermal.mass_kg * 60000.0 * (next.thermal().state.liquid_mass_fraction() - 0.5),
                -report.exchanged_heat_j,
                1e-7,
            );
        }
        assert_eq!(old, saved);
    }
}

#[test]
fn exact_equilibrium_and_zero_area_are_explicit_zero_transfer_without_a_card_claim() {
    let curve = source::body_curve();
    let old = piston(0.75, 400.0, 100000.0);
    let same_temperature = body(&curve, 100000.0).initial;
    let equilibrium = surface()
        .resolve(old.gas().state(), same_temperature)
        .unwrap();
    assert_eq!(equilibrium.conductance_w_k, 0.0);
    assert_eq!(equilibrium.source.area_m2, surface().area_m2);
    assert_eq!(equilibrium.source.gas, *old.gas().state());
    assert!(equilibrium.source.evaluation.is_none());
    // Equality is exact; a nonzero contrast below the native card's Ra
    // window must refuse instead of being rounded into equilibrium.
    let almost = curve.state_at_specific_enthalpy(100000.0 + 1e-7).unwrap();
    assert_ne!(
        almost.temperature_k(),
        old.gas().state().state().temperature
    );
    assert!(surface().resolve(old.gas().state(), almost).is_err());
    let thermal = body(&curve, 260000.0);
    let isolated = NaturalGasSurface {
        area_m2: 0.0,
        ..surface()
    }
    .finite_body_candidate(&old, input(), thermal, controls(), exchange(), || false)
    .unwrap();
    let explicit = old
        .thermal_candidate(
            input(),
            thermal,
            controls(),
            exchange(),
            |_, _| {
                Ok::<_, core::convert::Infallible>(ConductanceSample {
                    conductance_w_k: 0.0,
                    source: (),
                })
            },
            || false,
        )
        .unwrap();
    assert_eq!(isolated.piston(), explicit.piston());
    assert_eq!(isolated.thermal().state, explicit.thermal().state);
    assert_eq!(isolated.step(), explicit.step());
    assert!(isolated.thermal().report.evaluation.is_none());
}

#[test]
fn geometry_gravity_domain_and_late_cancellation_refuse_without_publishing() {
    let curve = source::body_curve();
    let old = piston(0.75, 400.0, 100000.0);
    let saved = old.clone();
    let thermal = body(&curve, 260000.0);
    for patch in [
        NaturalGasSurface {
            gravity_m_s2: 0.0,
            ..surface()
        },
        NaturalGasSurface {
            gravity_m_s2: -9.8,
            ..surface()
        },
        NaturalGasSurface {
            gravity_m_s2: f64::NAN,
            ..surface()
        },
        NaturalGasSurface {
            characteristic_length: Length::new(0.0),
            ..surface()
        },
        NaturalGasSurface {
            characteristic_length: Length::new(1e8),
            ..surface()
        },
        NaturalGasSurface {
            area_m2: 0.03,
            ..surface()
        },
    ] {
        assert!(
            patch
                .finite_body_candidate(&old, input(), thermal, controls(), exchange(), || false)
                .is_err()
        );
        assert_eq!(old, saved);
    }
    let mut count = 0;
    surface()
        .finite_body_candidate(&old, input(), thermal, controls(), exchange(), || {
            count += 1;
            false
        })
        .unwrap();
    for stop in [1, count / 2, count] {
        let mut calls = 0;
        let next =
            surface().finite_body_candidate(&old, input(), thermal, controls(), exchange(), || {
                calls += 1;
                calls == stop
            });
        assert!(matches!(next, Err(ThermalPistonError::Cancelled)));
        assert_eq!(old, saved);
        assert_eq!(
            thermal.initial,
            curve.state_at_specific_enthalpy(260000.0).unwrap()
        );
    }
}
