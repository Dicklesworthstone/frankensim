//! Source-resolved wall laws on the original fixed-volume reciprocal heat solve.
use super::*;
use fs_airflow::graph::thermal::medium::convection::{
    AidingVerticalPlateSettings,
    finite_gas::{
        film::{AidingVerticalPlateGasSurface, FilmGasHeatLaw, FilmGasSurface},
        sphere::ForcedSphereGasSurface,
    },
};
use fs_convection::{HorizontalPlateFacing, NaturalConvectionGeometry};
use fs_couple::finite_gas_thermal::{
    FiniteGasThermalReport, WallConvectionError, WallConvectionSolve,
};
use fs_material::fluid_transport::FluidTransport;
use fs_material::phase::UniformEnthalpyStep;

fn wall_solve() -> WallConvectionSolve {
    WallConvectionSolve {
        temperature_increment_bracket_k: [-15.0, 15.0],
        temperature_tolerance_k: 1e-8,
        maximum_iterations: 128,
    }
}

fn gas(temperature_k: f64) -> IsochoricGasParcel {
    IsochoricGasParcel::new(
        source::model_transport(1.0, true),
        source::point(temperature_k, 100000.0),
        0.003,
    )
    .unwrap()
}

fn natural(geometry: NaturalConvectionGeometry) -> FilmGasHeatLaw {
    FilmGasHeatLaw::Natural(NaturalGasSurface {
        geometry,
        area_m2: 0.02,
        characteristic_length: Length::new(0.1),
        gravity_m_s2: 9.80665,
    })
}

fn closed<S>(
    old: &IsochoricGasParcel,
    input: UniformEnthalpyStepInput<'_>,
    next: &UniformEnthalpyStep<FiniteGasThermalReport<S>>,
) {
    let report = &next.report;
    let gas_heat = old.mass_kg()
        * (report.next_gas.state().specific_internal_energy_j_kg()
            - old.state().specific_internal_energy_j_kg());
    let body_heat = input.mass_kg
        * (next.state.specific_enthalpy_j_kg() - input.initial.specific_enthalpy_j_kg());
    near(
        gas_heat,
        report.solid_to_gas_heat_j + report.external_gas_heat_j,
        1e-7,
    );
    near(
        body_heat,
        input.internal_heat_j + report.radiation_heat_j - report.solid_to_gas_heat_j,
        1e-7,
    );
    assert!(report.combined_energy_residual_j.abs() <= 1e-7);
    assert!(report.entropy.convection_production_j_k >= 0.0);
    near(report.entropy.balance_residual_j_k, 0.0, 1e-12);
}

#[test]
fn isochoric_natural_films_use_current_wall_and_full_rayleigh_contrast() {
    scope(|cx, _| {
        let old = gas(400.0);
        let saved = old.clone();
        let curve = source::body_curve();
        let thermal = UniformEnthalpyStepInput {
            duration_s: 0.2,
            ..body(&curve, 260000.0)
        };
        let boundary = uniform_boundary();
        let environment = uniform_environment(&old, &boundary);
        for geometry in [
            NaturalConvectionGeometry::VerticalPlate,
            NaturalConvectionGeometry::HorizontalPlate(HorizontalPlateFacing::Up),
            NaturalConvectionGeometry::HorizontalPlate(HorizontalPlateFacing::Down),
            NaturalConvectionGeometry::HorizontalCylinder,
            NaturalConvectionGeometry::Sphere,
        ] {
            let next = FilmGasSurface {
                model: old.model(),
                law: natural(geometry),
            }
            .isochoric_body_candidate(cx, &environment, thermal, wall_solve())
            .unwrap();
            let report = &next.report;
            let wall = &report.convection.source;
            let sample = &wall.source;
            let bulk = report.next_gas.state();
            let film = sample.property_gas.as_ref().unwrap();
            let evaluation = sample.interface.evaluation.as_ref().unwrap();
            assert_eq!(&sample.interface.gas, bulk);
            assert_eq!(sample.interface.body.value(), wall.temperature_k);
            assert_eq!(
                wall.endpoint_source.interface.body.value(),
                next.state.temperature_k()
            );
            near(
                film.state().temperature,
                0.5 * (wall.temperature_k + bulk.state().temperature),
                1e-12,
            );
            assert_eq!(film.state().pressure, bulk.state().pressure);
            assert_eq!(film.source_packs(), bulk.source_packs());
            assert_eq!(
                film.thermodynamics().receipt().declared_fraction_bits(),
                bulk.thermodynamics().receipt().declared_fraction_bits()
            );
            let transport = film.transport_state().unwrap();
            assert_eq!(evaluation.transport, transport);
            let ra = 9.80665 / film.state().temperature
                * (wall.temperature_k - bulk.state().temperature).abs()
                * 0.1_f64.powi(3)
                * transport.density_kg_m3().powi(2)
                * transport.specific_heat_j_kg_k()
                / (transport.dynamic_viscosity_pa_s() * transport.thermal_conductivity_w_m_k());
            near(evaluation.nusselt.groups()["Ra"], ra, ra * 1e-12);
            near(
                report.solid_to_gas_heat_j,
                thermal.duration_s
                    * thermal.surface_area_m2
                    * report.convection.coefficient_w_per_m2_k
                    * (next.state.temperature_k() - bulk.state().temperature),
                environment.exchange.source_tolerance_j,
            );
            assert!(wall.temperature_residual_k.abs() <= wall_solve().temperature_tolerance_k);
            assert!(wall.boundary_residual_j.abs() <= environment.exchange.source_tolerance_j);
            assert!(report.solid_to_gas_heat_j > 0.0);
            assert!(next.state.temperature_k() < thermal.initial.temperature_k());
            closed(&old, thermal, &next);
        }
        assert_eq!(old, saved);
    });
}

#[test]
fn isochoric_forced_and_aiding_films_share_the_existing_heat_solver() {
    scope(|cx, _| {
        let old = gas(400.0);
        let curve = source::body_curve();
        let thermal = UniformEnthalpyStepInput {
            duration_s: 0.2,
            ..body(&curve, 120000.0)
        };
        let boundary = uniform_boundary();
        let mut environment = uniform_environment(&old, &boundary);
        environment.exchange.heat_bracket_j = [-2.0, 2.0];
        let wall = WallConvectionSolve {
            temperature_increment_bracket_k: [-5.0, 5.0],
            ..wall_solve()
        };
        let forced = FilmGasSurface {
            model: old.model(),
            law: FilmGasHeatLaw::Forced(surface()),
        }
        .isochoric_body_candidate(cx, &environment, thermal, wall)
        .unwrap();
        closed(&old, thermal, &forced);
        let mut previous_heat = 0.0;
        for speed_up_m_s in [0.0, 0.25, 1.0] {
            let next = FilmGasSurface {
                model: old.model(),
                law: FilmGasHeatLaw::AidingVerticalPlate(AidingVerticalPlateGasSurface {
                    area_m2: thermal.surface_area_m2,
                    characteristic_length: Length::new(0.1),
                    settings: AidingVerticalPlateSettings {
                        speed_up_m_s,
                        gravity_m_s2: 9.80665,
                        maximum_relative_density_departure: 0.1,
                    },
                }),
            }
            .isochoric_body_candidate(cx, &environment, thermal, wall)
            .unwrap();
            let evaluated = next
                .report
                .convection
                .source
                .source
                .interface
                .evaluation
                .as_ref()
                .unwrap();
            assert_eq!(
                evaluated.nusselt.card().id,
                CorrelationId::ManningQureshiAidingVerticalPlate
            );
            assert_eq!(evaluated.nusselt.groups()["Re"] == 0.0, speed_up_m_s == 0.0);
            assert!(next.report.solid_to_gas_heat_j > previous_heat);
            previous_heat = next.report.solid_to_gas_heat_j;
            closed(&old, thermal, &next);
        }
    });
}

#[test]
fn isochoric_sphere_keeps_bulk_transport_surface_viscosity_and_actual_latent_heat() {
    scope(|cx, _| {
        let old = gas(600.0);
        let saved = old.clone();
        let curve = source::body_curve();
        let boundary = uniform_boundary();
        let environment = uniform_environment(&old, &boundary);
        let sphere = ForcedSphereGasSurface {
            model: old.model(),
            area_m2: 0.02,
            diameter: Length::new(0.1),
            speed_m_s: 1.0,
        };
        for enthalpy in [180000.0, 260000.0] {
            let thermal = UniformEnthalpyStepInput {
                duration_s: 0.2,
                ..body(&curve, enthalpy)
            };
            let next = sphere
                .isochoric_body_candidate(cx, &environment, thermal, wall_solve())
                .unwrap();
            let report = &next.report;
            let wall = &report.convection.source;
            let sample = &wall.source;
            let bulk = report.next_gas.state();
            let surface = sample.surface_gas.as_ref().unwrap();
            let evaluated = sample.interface.evaluation.as_ref().unwrap();
            let transport = bulk.transport_state().unwrap();
            assert_eq!(evaluated.transport, transport);
            assert_eq!(&sample.interface.gas, bulk);
            assert_eq!(surface.state().temperature, wall.temperature_k);
            assert_eq!(surface.state().pressure, bulk.state().pressure);
            assert_eq!(surface.source_packs(), bulk.source_packs());
            near(
                wall.endpoint_source
                    .surface_gas
                    .as_ref()
                    .unwrap()
                    .state()
                    .temperature,
                next.state.temperature_k(),
                1e-12,
            );
            let re = transport.density_kg_m3() * sphere.speed_m_s * sphere.diameter.value()
                / transport.dynamic_viscosity_pa_s();
            let ratio = transport.dynamic_viscosity_pa_s()
                / surface.transport_state().unwrap().dynamic_viscosity_pa_s();
            let nu = 2.0
                + (0.4 * re.sqrt() + 0.06 * re.powf(2.0 / 3.0))
                    * transport.prandtl().powf(0.4)
                    * ratio.powf(0.25);
            near(
                evaluated.nusselt.groups()["mu_bulk_over_surface"],
                ratio,
                1e-13,
            );
            near(
                report.convection.coefficient_w_per_m2_k,
                nu * transport.thermal_conductivity_w_m_k() / sphere.diameter.value(),
                1e-12,
            );
            near(
                report.solid_to_gas_heat_j,
                thermal.duration_s
                    * thermal.surface_area_m2
                    * report.convection.coefficient_w_per_m2_k
                    * (next.state.temperature_k() - bulk.state().temperature),
                environment.exchange.source_tolerance_j,
            );
            assert!(wall.temperature_residual_k.abs() <= wall_solve().temperature_tolerance_k);
            assert!(wall.boundary_residual_j.abs() <= environment.exchange.source_tolerance_j);
            assert!(report.solid_to_gas_heat_j < 0.0);
            if enthalpy == 180000.0 {
                assert_eq!(next.state.temperature_k(), 450.0);
                assert!(next.state.liquid_mass_fraction() > 0.5);
            } else {
                assert!(next.state.temperature_k() > thermal.initial.temperature_k());
            }
            closed(&old, thermal, &next);
        }
        assert_eq!(old, saved);
    });
}

#[test]
fn isochoric_wall_adapters_preserve_inputs_on_source_area_biot_and_cancellation_refusal() {
    scope(|cx, gate| {
        let old = gas(400.0);
        let saved = old.clone();
        let curve = source::body_curve();
        let thermal = UniformEnthalpyStepInput {
            duration_s: 0.2,
            ..body(&curve, 260000.0)
        };
        let boundary = uniform_boundary();
        let environment = uniform_environment(&old, &boundary);
        let film = FilmGasSurface {
            model: old.model(),
            law: FilmGasHeatLaw::Forced(surface()),
        };
        let sphere = ForcedSphereGasSurface {
            model: old.model(),
            area_m2: 0.02,
            diameter: Length::new(0.1),
            speed_m_s: 1.0,
        };
        let partial = FilmGasSurface {
            law: FilmGasHeatLaw::Forced(ForcedGasSurface {
                area_m2: 0.01,
                ..surface()
            }),
            ..film
        };
        assert!(matches!(
            partial.isochoric_body_candidate(cx, &environment, thermal, wall_solve()),
            Err(ConvectionError::Thermal(ThermalPistonError::Invalid(_)))
        ));
        assert!(matches!(
            (ForcedSphereGasSurface {
                area_m2: 0.01,
                ..sphere
            })
            .isochoric_body_candidate(cx, &environment, thermal, wall_solve()),
            Err(ConvectionError::Thermal(ThermalPistonError::Invalid(_)))
        ));
        // The actual cold-gas/hot-wall viscosity ratio is outside Whitaker's card.
        assert!(matches!(
            sphere.isochoric_body_candidate(cx, &environment, thermal, wall_solve()),
            Err(ConvectionError::Convection(WallConvectionError::Source(_)))
        ));
        let heat_flux = FilmGasSurface {
            law: FilmGasHeatLaw::Forced(ForcedGasSurface {
                correlation: CorrelationId::CircularDuctLaminarChf,
                ..surface()
            }),
            ..film
        };
        assert!(
            heat_flux
                .isochoric_body_candidate(cx, &environment, thermal, wall_solve())
                .is_err()
        );
        let resistive = LumpedThermalEnvironment {
            transport: LumpedThermalTransport::try_declared(1e-6, 0.0).unwrap(),
            ..boundary.clone()
        };
        assert!(
            film.isochoric_body_candidate(
                cx,
                &uniform_environment(&old, &resistive),
                thermal,
                wall_solve()
            )
            .is_err()
        );
        let insulated = ForcedSphereGasSurface {
            area_m2: 0.0,
            speed_m_s: 0.0,
            ..sphere
        }
        .isochoric_body_candidate(cx, &environment, thermal, wall_solve())
        .unwrap();
        assert_eq!(insulated.state, thermal.initial);
        assert_eq!(insulated.report.next_gas, old);
        assert!(
            insulated
                .report
                .convection
                .source
                .source
                .surface_gas
                .is_none()
        );
        gate.request();
        assert!(matches!(
            film.isochoric_body_candidate(cx, &environment, thermal, wall_solve()),
            Err(ConvectionError::Thermal(ThermalPistonError::Cancelled))
        ));
        assert!(matches!(
            sphere.isochoric_body_candidate(cx, &environment, thermal, wall_solve()),
            Err(ConvectionError::Thermal(ThermalPistonError::Cancelled))
        ));
        assert_eq!(old, saved);
    });
}
