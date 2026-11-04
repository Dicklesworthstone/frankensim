//! Physical leading-edge history on the existing reciprocal finite-gas solves.
use super::*;
use fs_airflow::graph::thermal::medium::convection::{
    evaluate_flat_plate_section_from_material,
    finite_gas::film::{FilmGasHeatLaw, FilmGasSurface, ForcedPlateGasSurface},
};
use fs_couple::finite_gas_thermal::WallConvectionSolve;
use fs_material::fluid_transport::FluidTransport;
use fs_qty::semantic::{Composition, CompositionBasis};

fn section(start: f64, end: f64, area_m2: f64) -> ForcedPlateGasSurface {
    ForcedPlateGasSurface {
        area_m2,
        start: Length::new(start),
        end: Length::new(end),
        speed_m_s: 1.0,
        correlation: CorrelationId::FlatPlateLaminarAverage,
    }
}

fn wall_solve() -> WallConvectionSolve {
    WallConvectionSolve {
        temperature_increment_bracket_k: [-15.0, 15.0],
        temperature_tolerance_k: 1e-8,
        maximum_iterations: 128,
    }
}

#[test]
fn physical_half_plates_recover_whole_heat_without_restarting_the_boundary_layer() {
    let bank = source::model_transport(0.75, true);
    let actual = bank
        .with_composition(
            Composition::new(CompositionBasis::MoleFraction, vec![0.25, 0.75]).unwrap(),
        )
        .unwrap();
    let bulk = actual.evaluate(&source::point(400.0, 200000.0)).unwrap();
    let curve = source::body_curve();
    let wall = body(&curve, 260000.0).initial;
    let resolve = |patch| {
        FilmGasSurface {
            model: &bank,
            law: FilmGasHeatLaw::ForcedPlateSection(patch),
        }
        .resolve(&bulk, wall)
        .unwrap()
    };
    let whole = resolve(section(0.0, 0.2, 0.02));
    let first = resolve(section(0.0, 0.1, 0.01));
    let second = resolve(section(0.1, 0.2, 0.01));
    let original = FilmGasSurface {
        model: &bank,
        law: FilmGasHeatLaw::Forced(ForcedGasSurface {
            characteristic_length: Length::new(0.2),
            ..surface()
        }),
    }
    .resolve(&bulk, wall)
    .unwrap();
    assert_eq!(whole.conductance_w_k, original.conductance_w_k);

    // Independent integral of h(x) proportional to x^(-1/2).
    let root_two = 2.0_f64.sqrt();
    near(first.conductance_w_k / whole.conductance_w_k, root_two / 2.0, 2e-13);
    near(second.conductance_w_k / whole.conductance_w_k, 1.0 - root_two / 2.0, 2e-13);
    near(
        (first.conductance_w_k + second.conductance_w_k) / whole.conductance_w_k,
        1.0,
        2e-13,
    );
    // Restarting at both cell boundaries would overpredict by sqrt(2).
    near(2.0 * first.conductance_w_k / whole.conductance_w_k, root_two, 2e-13);

    for result in [&whole, &first, &second] {
        let sample = &result.source;
        let film = sample.property_gas.as_ref().unwrap();
        let evaluation = sample.interface.evaluation.as_ref().unwrap();
        let strip = evaluation.section.as_ref().unwrap();
        assert_eq!(sample.interface.gas, bulk);
        assert_eq!(film.state().temperature, 450.0);
        assert_eq!(film.state().pressure, bulk.state().pressure);
        assert_eq!(film.source_packs(), bulk.source_packs());
        assert_eq!(
            film.thermodynamics().receipt().declared_fraction_bits(),
            bulk.thermodynamics().receipt().declared_fraction_bits()
        );
        let transport = film.transport_state().unwrap();
        assert_eq!(evaluation.transport, transport);
        near(
            evaluation.effective_nusselt(),
            evaluation.coefficient.value.value() * strip.end.value()
                / transport.thermal_conductivity_w_m_k(),
            1e-11,
        );
        let re_end = transport.density_kg_m3() * strip.end.value()
            / transport.dynamic_viscosity_pa_s();
        near(evaluation.nusselt.groups()["Re"], re_end, re_end * 2e-13);
        if strip.start.value() == 0.0 {
            assert!(strip.upstream.is_none());
            assert_eq!(evaluation.effective_nusselt(), evaluation.nusselt.evidence().value);
            assert!(evaluation.coefficient.model.discrepancy_rel.is_finite());
        } else {
            let upstream = strip.upstream.as_ref().unwrap();
            near(
                upstream.groups()["Re"],
                re_end * strip.start.value() / strip.end.value(),
                re_end * 2e-13,
            );
            assert!(evaluation.effective_nusselt() < evaluation.nusselt.evidence().value);
            assert!(evaluation.coefficient.model.discrepancy_rel.is_infinite());
        }
    }
}

#[test]
fn section_heat_closes_body_gas_and_motion_in_both_directions_and_through_latent_heat() {
    let curve = source::body_curve();
    let old = piston(0.75, 400.0, 100000.0);
    let saved = old.clone();
    let patch = section(0.1, 0.2, 0.02);
    let law = FilmGasSurface {
        model: old.gas().model(),
        law: FilmGasHeatLaw::ForcedPlateSection(patch),
    };
    for enthalpy in [50000.0, 180000.0, 260000.0] {
        let thermal = body(&curve, enthalpy);
        let next = law
            .finite_body_candidate(&old, input(), thermal, controls(), exchange(), || false)
            .unwrap();
        let report = next.step();
        let sample = &next.thermal().report;
        let film = sample.property_gas.as_ref().unwrap();
        let evaluation = sample.interface.evaluation.as_ref().unwrap();
        assert!(evaluation.section.as_ref().unwrap().upstream.is_some());
        assert_eq!(sample.interface.gas, report.piston.midpoint);
        assert_eq!(sample.interface.body, report.thermal_midpoint);
        near(
            film.state().temperature,
            0.5 * (report.thermal_midpoint.temperature_k() + report.piston.midpoint.state().temperature),
            1e-12,
        );
        assert_eq!(film.state().pressure, report.piston.midpoint.state().pressure);
        assert_eq!(film.source_packs(), old.gas().state().source_packs());
        near(
            report.exchanged_heat_j,
            DT * patch.area_m2 * evaluation.coefficient.value.value()
                * (report.thermal_midpoint.temperature_k()
                    - report.piston.midpoint.state().temperature),
            exchange().source_tolerance_j,
        );
        assert!(report.exchanged_heat_j * (thermal.initial.temperature_k() - 400.0) > 0.0);
        let gas_heat = old.gas().mass_kg()
            * (next.piston().gas().state().specific_internal_energy_j_kg()
                - old.gas().state().specific_internal_energy_j_kg());
        let body_heat = thermal.mass_kg
            * (next.thermal().state.specific_enthalpy_j_kg()
                - thermal.initial.specific_enthalpy_j_kg());
        let kinetic = 0.5 * old.mass_kg() * next.piston().velocity_m_s().powi(2);
        let outside_work = input().exterior_pressure_pa
            * (next.piston().gas().volume_m3() - old.gas().volume_m3());
        near(gas_heat + body_heat + kinetic + outside_work, 0.0, 1e-7);
        near(body_heat + report.exchanged_heat_j, 0.0, 1e-7);
        if enthalpy == 180000.0 {
            assert_eq!(next.thermal().state.temperature_k(), 450.0);
            assert!(next.thermal().state.liquid_mass_fraction() < 0.5);
        }
        assert_eq!(old, saved);
    }
}

#[test]
fn fixed_volume_section_refreshes_wall_film_and_preserves_both_energy_stores() {
    scope(|cx, _| {
        let old = IsochoricGasParcel::new(
            source::model_transport(0.75, true),
            source::point(350.0, 100000.0),
            0.001,
        )
        .unwrap();
        let saved = old.clone();
        let curve = source::body_curve();
        let boundary = uniform_boundary();
        let mut environment = uniform_environment(&old, &boundary);
        environment.gas_heat_j = 0.2;
        let thermal = UniformEnthalpyStepInput {
            duration_s: 0.1,
            internal_heat_j: 0.1,
            ..body(&curve, 260000.0)
        };
        let patch = section(0.1, 0.2, thermal.surface_area_m2);
        let law = FilmGasSurface {
            model: old.model(),
            law: FilmGasHeatLaw::ForcedPlateSection(patch),
        };
        let initial = law.resolve(old.state(), thermal.initial).unwrap();
        let next = law
            .isochoric_body_candidate(cx, &environment, thermal, wall_solve())
            .unwrap();
        let report = &next.report;
        let endpoint = report.next_gas.state();
        let wall = &report.convection.source;
        let sample = &wall.source;
        let film = sample.property_gas.as_ref().unwrap();
        let evaluation = sample.interface.evaluation.as_ref().unwrap();
        assert_eq!(&sample.interface.gas, endpoint);
        assert_eq!(sample.interface.body.value(), wall.temperature_k);
        assert_eq!(wall.endpoint_source.interface.body.value(), next.state.temperature_k());
        near(
            film.state().temperature,
            0.5 * (wall.temperature_k + endpoint.state().temperature),
            1e-12,
        );
        assert_eq!(film.state().pressure, endpoint.state().pressure);
        assert_eq!(film.source_packs(), old.state().source_packs());
        let fresh = evaluate_flat_plate_section_from_material(
            patch.correlation,
            film,
            patch.speed_m_s,
            patch.start,
            patch.end,
        )
        .unwrap();
        assert_eq!(report.convection.coefficient_w_per_m2_k, fresh.coefficient().value.value());
        assert_eq!(evaluation.coefficient.value.value(), fresh.coefficient().value.value());
        assert!(
            (report.convection.coefficient_w_per_m2_k
                - initial.conductance_w_k / patch.area_m2)
                .abs()
                > 1e-4
        );
        near(
            report.solid_to_gas_heat_j,
            thermal.duration_s * thermal.surface_area_m2
                * report.convection.coefficient_w_per_m2_k
                * (next.state.temperature_k() - endpoint.state().temperature),
            environment.exchange.source_tolerance_j,
        );
        let gas_heat = old.mass_kg()
            * (endpoint.specific_internal_energy_j_kg()
                - old.state().specific_internal_energy_j_kg());
        let body_heat = thermal.mass_kg
            * (next.state.specific_enthalpy_j_kg() - thermal.initial.specific_enthalpy_j_kg());
        near(gas_heat, report.solid_to_gas_heat_j + environment.gas_heat_j, 1e-7);
        near(body_heat, thermal.internal_heat_j - report.solid_to_gas_heat_j, 1e-7);
        near(gas_heat + body_heat, environment.gas_heat_j + thermal.internal_heat_j, 1e-7);
        assert!(wall.temperature_residual_k.abs() <= wall_solve().temperature_tolerance_k);
        assert!(wall.boundary_residual_j.abs() <= environment.exchange.source_tolerance_j);
        assert_eq!(old, saved);
    });
}

#[test]
fn unsupported_section_geometry_and_domains_refuse_without_publishing_state() {
    let old = piston(0.75, 400.0, 100000.0);
    let saved = old.clone();
    let curve = source::body_curve();
    let thermal = body(&curve, 260000.0);
    let valid = section(0.1, 0.2, 0.02);
    for patch in [
        ForcedPlateGasSurface { start: Length::new(-0.1), ..valid },
        ForcedPlateGasSurface { start: valid.end, ..valid },
        ForcedPlateGasSurface { end: Length::new(f64::INFINITY), ..valid },
        ForcedPlateGasSurface { start: Length::new(1e-12), ..valid },
        ForcedPlateGasSurface { speed_m_s: 1e20, ..valid },
        ForcedPlateGasSurface { area_m2: 0.03, ..valid },
        ForcedPlateGasSurface { correlation: CorrelationId::ChurchillBernsteinCylinder, ..valid },
    ] {
        let law = FilmGasSurface {
            model: old.gas().model(),
            law: FilmGasHeatLaw::ForcedPlateSection(patch),
        };
        assert!(law
            .finite_body_candidate(&old, input(), thermal, controls(), exchange(), || false)
            .is_err());
        assert_eq!(old, saved);
    }
    let law = FilmGasSurface {
        model: old.gas().model(),
        law: FilmGasHeatLaw::ForcedPlateSection(valid),
    };
    assert!(matches!(
        law.finite_body_candidate(&old, input(), thermal, controls(), exchange(), || true),
        Err(ThermalPistonError::Cancelled)
    ));
    assert_eq!(old, saved);
}

#[test]
fn zero_area_section_keeps_the_original_closed_solver_and_skips_property_queries() {
    let old = piston(0.75, 400.0, 100000.0);
    let curve = source::body_curve();
    let thermal = body(&curve, 260000.0);
    let law = FilmGasSurface {
        model: old.gas().model(),
        law: FilmGasHeatLaw::ForcedPlateSection(ForcedPlateGasSurface {
            speed_m_s: 0.0,
            ..section(0.1, 0.2, 0.0)
        }),
    };
    let isolated = law
        .finite_body_candidate(&old, input(), thermal, controls(), exchange(), || false)
        .unwrap();
    let original = old
        .thermal_candidate(
            input(), thermal, controls(), exchange(),
            |_, _| Ok::<_, core::convert::Infallible>(ConductanceSample {
                conductance_w_k: 0.0,
                source: (),
            }),
            || false,
        )
        .unwrap();
    assert_eq!(isolated.piston(), original.piston());
    assert_eq!(isolated.thermal().state, original.thermal().state);
    assert_eq!(isolated.step(), original.step());
    assert!(isolated.thermal().report.property_gas.is_none());
    assert!(isolated.thermal().report.interface.evaluation.is_none());
}
