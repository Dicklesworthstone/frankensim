//! Source-owned real-gas expansion reaches the existing buoyant heat/motion residual.
use super::*;
use fs_airflow::graph::thermal::medium::convection::{
    MaterialConvectionError, NaturalConvectionGeometry,
    evaluate_vertical_plate_from_thermal_material, finite_gas::NaturalGasSurface,
};
use fs_material::fluid_transport::ThermalFluidTransport;
use fs_qty::Length;

fn natural_surface() -> NaturalGasSurface {
    NaturalGasSurface {
        geometry: NaturalConvectionGeometry::VerticalPlate,
        area_m2: 0.02,
        characteristic_length: Length::new(0.1),
        gravity_m_s2: 9.80665,
    }
}

#[test]
fn dense_gas_expansion_and_density_jointly_determine_buoyancy_and_conductance() {
    let curve = source::body_curve();
    let old = gas(350.0, 0.0);
    let state = old.state();
    let patch = natural_surface();
    let wall = body(&curve, 55000.0).initial; // 355 K: a small bulk-property contrast.
    let resolved = patch.resolve(state, wall).unwrap();
    let evaluated = resolved.source.evaluation.as_ref().unwrap();
    let thermal = state.thermal_transport_state().unwrap();
    let transport = thermal.transport();
    let beta = thermal.volumetric_expansion_per_k().unwrap();
    assert_eq!(thermal.temperature_k(), state.temperature_k());
    assert_eq!(beta, state.volumetric_expansion_per_k());
    assert!((beta * state.temperature_k() - 1.0).abs() > 0.01);
    assert_eq!(transport, state.transport_state().unwrap());
    assert_eq!(transport.density_kg_m3(), state.density_kg_m3());
    assert_ne!(transport.density_kg_m3(), state.source().state().density);
    assert_eq!(
        transport.specific_heat_j_kg_k(),
        state.specific_heat_cp_j_kg_k()
    );

    // Reconstruct the groups and Churchill-Chu response independently from the
    // accepted source values, including the actual EOS derivative in Ra.
    let rho = transport.density_kg_m3();
    let cp = transport.specific_heat_j_kg_k();
    let mu = transport.dynamic_viscosity_pa_s();
    let k = transport.thermal_conductivity_w_m_k();
    let length = patch.characteristic_length.value();
    let pr = mu * cp / k;
    let ra = patch.gravity_m_s2
        * beta
        * (wall.temperature_k() - state.temperature_k())
        * length.powi(3)
        * rho.powi(2)
        * cp
        / (mu * k);
    let nusselt = |ra: f64| {
        (0.825
            + 0.387 * ra.powf(1.0 / 6.0) / (1.0 + (0.492 / pr).powf(9.0 / 16.0)).powf(8.0 / 27.0))
        .powi(2)
    };
    near(evaluated.nusselt.groups()["Ra"], ra, ra * 1e-12);
    near(evaluated.nusselt.groups()["Pr"], pr, pr * 1e-12);
    near(
        resolved.conductance_w_k,
        patch.area_m2 * k * nusselt(ra) / length,
        1e-11,
    );
    let ideal_expansion_ra = ra / (beta * state.temperature_k());
    let ideal_expansion_g = patch.area_m2 * k * nusselt(ideal_expansion_ra) / length;
    assert!((resolved.conductance_w_k / ideal_expansion_g - 1.0).abs() > 1e-3);
    assert_eq!(evaluated.transport, transport);
    assert_eq!(resolved.source.gas, *state);
    assert_eq!(resolved.source.body, wall);
    assert_eq!(
        resolved.source.gas.source().source_packs(),
        old.model().source().source_packs()
    );
    near(
        patch
            .resolve(state, body(&curve, 45000.0).initial)
            .unwrap()
            .conductance_w_k,
        resolved.conductance_w_k,
        1e-13,
    );
}

#[test]
fn buoyant_heat_closes_real_gas_body_and_piston_energy_for_both_heat_directions() {
    let curve = source::body_curve();
    let patch = natural_surface();
    for enthalpy in [45000.0, 55000.0] {
        // 345/355 K around a 350 K gas.
        let old = GasPiston::new(gas(350.0, 0.0), 2.0, 0.01, 0.15).unwrap();
        let saved = old.clone();
        let thermal = UniformEnthalpyStepInput {
            duration_s: moving_input().duration_s,
            ..body(&curve, enthalpy)
        };
        let candidate = patch
            .finite_body_candidate(
                &old,
                moving_input(),
                thermal,
                moving_controls(),
                super::super::exchange(),
                || false,
            )
            .unwrap();
        let next = candidate.piston();
        let report = candidate.step();
        let sample = &candidate.thermal().report;
        let midpoint = &report.piston.midpoint;
        let evaluated = sample.evaluation.as_ref().unwrap();
        assert_eq!(sample.gas, *midpoint);
        assert_eq!(sample.body, report.thermal_midpoint);
        assert_eq!(evaluated.transport, midpoint.transport_state().unwrap());
        assert_eq!(
            sample.gas.source().source_packs(),
            old.gas().model().source().source_packs()
        );
        let checked = evaluate_vertical_plate_from_thermal_material(
            midpoint,
            patch.gravity_m_s2,
            Temperature::new(sample.body.temperature_k()),
            patch.characteristic_length,
        )
        .unwrap();
        assert_eq!(evaluated.coefficient.value, checked.coefficient().value);
        let heat = report.exchanged_heat_j;
        assert_eq!(heat > 0.0, enthalpy > 50000.0);
        near(
            heat,
            moving_input().duration_s
                * report.conductance_w_k
                * (sample.body.temperature_k() - midpoint.temperature_k()),
            super::super::exchange().source_tolerance_j,
        );
        assert_eq!(candidate.thermal().external_heat_j, -heat);
        let gas_energy = old.gas().mass_kg()
            * old
                .gas()
                .model()
                .specific_internal_energy_change_between_densities(
                    old.gas().density_kg_m3(),
                    old.gas().state().temperature_k(),
                    next.gas().density_kg_m3(),
                    next.gas().state().temperature_k(),
                )
                .unwrap();
        let body_energy = thermal.mass_kg
            * (candidate.thermal().state.specific_enthalpy_j_kg()
                - thermal.initial.specific_enthalpy_j_kg());
        let kinetic = 0.5
            * old.mass_kg()
            * (next.velocity_m_s() - old.velocity_m_s())
            * (next.velocity_m_s() + old.velocity_m_s());
        let dv = next.gas().volume_m3() - old.gas().volume_m3();
        near(
            gas_energy + body_energy + kinetic + moving_input().exterior_pressure_pa * dv,
            0.0,
            super::super::exchange().energy_tolerance_j,
        );
        near(
            body_energy + heat,
            0.0,
            super::super::exchange().energy_tolerance_j,
        );
        assert!(
            report.total_energy_residual_j.abs() <= super::super::exchange().energy_tolerance_j
        );
        assert_eq!(next.gas().mass_kg(), old.gas().mass_kg());
        assert_eq!(
            next.gas().density_kg_m3(),
            next.gas().mass_kg() / next.gas().volume_m3()
        );
        assert!(next.gas().volume_m3() > old.gas().volume_m3());
        assert_eq!(old, saved);
        assert_eq!(
            thermal.initial,
            curve.state_at_specific_enthalpy(enthalpy).unwrap()
        );
    }
}

#[test]
fn real_gas_limits_domain_refusals_and_late_cancellation_preserve_inputs() {
    let curve = source::body_curve();
    let old = GasPiston::new(gas(350.0, 0.0), 2.0, 0.01, 0.15).unwrap();
    let saved = old.clone();
    let thermal = UniformEnthalpyStepInput {
        duration_s: moving_input().duration_s,
        ..body(&curve, 55000.0)
    };
    let patch = natural_surface();
    let equilibrium = patch
        .resolve(old.gas().state(), body(&curve, 50000.0).initial)
        .unwrap();
    let isolated = NaturalGasSurface {
        area_m2: 0.0,
        ..patch
    }
    .resolve(old.gas().state(), thermal.initial)
    .unwrap();
    for sample in [equilibrium, isolated] {
        assert_eq!(sample.conductance_w_k, 0.0);
        assert!(sample.source.evaluation.is_none());
        assert_eq!(sample.source.gas, *old.gas().state());
    }
    let outside_domain = NaturalGasSurface {
        characteristic_length: Length::new(1e-8),
        ..patch
    };
    assert!(matches!(
        outside_domain.resolve(old.gas().state(), thermal.initial),
        Err(MaterialConvectionError::Correlation(_))
    ));
    for refused in [
        outside_domain,
        NaturalGasSurface {
            geometry: NaturalConvectionGeometry::HorizontalCylinder,
            ..patch
        },
    ] {
        assert!(matches!(
            refused.finite_body_candidate(
                &old,
                moving_input(),
                thermal,
                moving_controls(),
                super::super::exchange(),
                || false,
            ),
            Err(ThermalPistonError::Interface(_))
        ));
    }
    let mut total = 0;
    patch
        .finite_body_candidate(
            &old,
            moving_input(),
            thermal,
            moving_controls(),
            super::super::exchange(),
            || {
                total += 1;
                false
            },
        )
        .unwrap();
    let mut polls = 0;
    assert!(matches!(
        patch.finite_body_candidate(
            &old,
            moving_input(),
            thermal,
            moving_controls(),
            super::super::exchange(),
            || {
                polls += 1;
                polls == total
            },
        ),
        Err(ThermalPistonError::Cancelled)
    ));
    assert_eq!(old, saved);
    assert_eq!(
        thermal.initial,
        curve.state_at_specific_enthalpy(55000.0).unwrap()
    );
}
