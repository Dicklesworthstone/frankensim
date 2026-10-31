//! Real film properties share the original source, correlation and heat/motion owners.
use super::*;
use fs_airflow::graph::thermal::medium::convection::{
    AidingVerticalPlateSettings,
    finite_gas::{
        NaturalGasSurface,
        film::{
            AidingVerticalPlateGasSurface, FilmGasHeatLaw, FilmGasSurface, PengRobinsonFilmSource,
        },
    },
};
use fs_convection::{CorrelationId, HorizontalPlateFacing, NaturalConvectionGeometry};
use fs_material::fluid_transport::ThermalFluidTransport;
use fs_qty::Length;

fn forced() -> FilmGasHeatLaw {
    FilmGasHeatLaw::Forced(ForcedGasSurface {
        characteristic_length: Length::new(0.01),
        speed_m_s: 0.5,
        ..surface()
    })
}

fn natural(geometry: NaturalConvectionGeometry) -> FilmGasHeatLaw {
    FilmGasHeatLaw::Natural(NaturalGasSurface {
        geometry,
        area_m2: 0.02,
        characteristic_length: Length::new(0.01),
        gravity_m_s2: 9.80665,
    })
}

fn mixed(speed: f64, departure: f64) -> FilmGasHeatLaw {
    FilmGasHeatLaw::AidingVerticalPlate(AidingVerticalPlateGasSurface {
        area_m2: 0.02,
        characteristic_length: Length::new(0.01),
        settings: AidingVerticalPlateSettings {
            speed_up_m_s: speed,
            gravity_m_s2: 9.80665,
            maximum_relative_density_departure: departure,
        },
    })
}

#[test]
fn forced_oriented_natural_and_aiding_films_use_real_density_cp_and_expansion() {
    let old = gas(350.0, 0.0);
    let saved = old.clone();
    let curve = source::body_curve();
    let wall = body(&curve, 55000.0).initial; // 355 K; the full contrast is 5 K.
    let source = PengRobinsonFilmSource {
        model: old.model(),
        phase: PengRobinsonPhaseV1::Vapor,
    };
    for law in [
        forced(),
        natural(NaturalConvectionGeometry::VerticalPlate),
        natural(NaturalConvectionGeometry::HorizontalPlate(
            HorizontalPlateFacing::Up,
        )),
        natural(NaturalConvectionGeometry::HorizontalPlate(
            HorizontalPlateFacing::Down,
        )),
        natural(NaturalConvectionGeometry::HorizontalCylinder),
        mixed(0.0, 0.1),
        mixed(0.25, 0.1),
    ] {
        let resolved = FilmGasSurface {
            model: &source,
            law,
        }
        .resolve(old.state(), wall)
        .unwrap();
        let film = resolved.source.property_gas.as_ref().unwrap();
        let evaluated = resolved.source.interface.evaluation.as_ref().unwrap();
        let thermal = film.thermal_transport_state().unwrap();
        let transport = thermal.transport();
        let beta = thermal.volumetric_expansion_per_k().unwrap();
        assert_eq!(film.temperature_k(), 352.5);
        assert_eq!(film.pressure_pa(), old.state().pressure_pa());
        assert!(std::ptr::eq(
            film.constitutive_law(),
            old.state().constitutive_law()
        ));
        assert_eq!(
            film.source().source_packs(),
            old.state().source().source_packs()
        );
        assert_eq!(
            film.source()
                .thermodynamics()
                .receipt()
                .declared_fraction_bits(),
            old.state()
                .source()
                .thermodynamics()
                .receipt()
                .declared_fraction_bits()
        );
        assert_eq!(resolved.source.interface.gas, *old.state());
        assert_eq!(resolved.source.interface.body, wall);
        assert_eq!(evaluated.transport, transport);
        assert_eq!(transport.density_kg_m3(), film.density_kg_m3());
        assert_eq!(
            transport.specific_heat_j_kg_k(),
            film.specific_heat_cp_j_kg_k()
        );
        assert_ne!(transport.density_kg_m3(), film.source().state().density);
        assert!((beta * film.temperature_k() - 1.0).abs() > 0.01);
        let rho = transport.density_kg_m3();
        let cp = transport.specific_heat_j_kg_k();
        let mu = transport.dynamic_viscosity_pa_s();
        let k = transport.thermal_conductivity_w_m_k();
        let pr = mu * cp / k;
        let ra = 9.80665 * beta * 5.0 * 0.01_f64.powi(3) * rho.powi(2) * cp / (mu * k);
        let groups = evaluated.nusselt.groups();
        let nu = match law {
            FilmGasHeatLaw::Forced(surface) => {
                let re = rho * surface.speed_m_s * 0.01 / mu;
                near(groups["Re"], re, re * 1e-12);
                near(groups["Pr"], pr, pr * 1e-12);
                0.664 * re.sqrt() * pr.powf(1.0 / 3.0)
            }
            FilmGasHeatLaw::Natural(surface) => {
                near(groups["Ra"], ra, ra * 1e-12);
                match surface.geometry {
                    NaturalConvectionGeometry::VerticalPlate => (0.825
                        + 0.387 * ra.powf(1.0 / 6.0)
                            / (1.0 + (0.492 / pr).powf(9.0 / 16.0)).powf(8.0 / 27.0))
                    .powi(2),
                    NaturalConvectionGeometry::HorizontalCylinder => (0.6
                        + 0.387 * ra.powf(1.0 / 6.0)
                            / (1.0 + (0.559 / pr).powf(9.0 / 16.0)).powf(8.0 / 27.0))
                    .powi(2),
                    NaturalConvectionGeometry::HorizontalPlate(HorizontalPlateFacing::Up) => {
                        assert_eq!(
                            evaluated.nusselt.card().id,
                            CorrelationId::HorizontalPlateUnstableLaminar
                        );
                        0.54 * ra.powf(0.25)
                    }
                    NaturalConvectionGeometry::HorizontalPlate(HorizontalPlateFacing::Down) => {
                        assert_eq!(
                            evaluated.nusselt.card().id,
                            CorrelationId::HorizontalPlateStable
                        );
                        0.27 * ra.powf(0.25)
                    }
                    _ => {
                        unreachable!("only the declared plate/cylinder fixtures enter this oracle")
                    }
                }
            }
            FilmGasHeatLaw::AidingVerticalPlate(surface) => {
                let re = rho * surface.settings.speed_up_m_s * 0.01 / mu;
                near(groups["Re"], re, re.max(1.0) * 1e-12);
                near(groups["Ra"], ra, ra * 1e-12);
                near(groups["Pr"], pr, pr * 1e-12);
                let forced = 2.0 * 0.339 * re.sqrt() * pr.powf(1.0 / 3.0)
                    / (1.0 + 0.100 * pr.powf(-0.75)).powf(2.0 / 9.0);
                let natural = (4.0 * 0.503 / 3.0) * ra.powf(0.25)
                    / (1.0 + 0.670 * pr.powf(-9.0 / 16.0)).powf(4.0 / 9.0);
                let exponent = 3.5 * pr.powf(0.075);
                (forced.powf(exponent) + natural.powf(exponent)).powf(1.0 / exponent)
            }
        };
        near(resolved.conductance_w_k, 0.02 * k * nu / 0.01, 1e-10);
    }
    assert_eq!(old, saved);
}

#[test]
fn film_heat_closes_real_gas_body_and_piston_storage_in_both_directions() {
    let curve = source::body_curve();
    for enthalpy in [45000.0, 55000.0] {
        let old = GasPiston::new(gas(350.0, 0.0), 2.0, 0.01, 0.15).unwrap();
        let saved = old.clone();
        let thermal = UniformEnthalpyStepInput {
            duration_s: moving_input().duration_s,
            ..body(&curve, enthalpy)
        };
        let source = PengRobinsonFilmSource {
            model: old.gas().model(),
            phase: PengRobinsonPhaseV1::Vapor,
        };
        let speed = if enthalpy > 50000.0 { 0.25 } else { -0.25 };
        for law in [
            forced(),
            natural(NaturalConvectionGeometry::HorizontalCylinder),
            mixed(speed, 0.1),
        ] {
            let candidate = FilmGasSurface {
                model: &source,
                law,
            }
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
            let bulk = &sample.interface.gas;
            let film = sample.property_gas.as_ref().unwrap();
            assert_eq!(*bulk, report.piston.midpoint);
            assert_eq!(sample.interface.body, report.thermal_midpoint);
            near(
                film.temperature_k(),
                0.5 * bulk.temperature_k() + 0.5 * sample.interface.body.temperature_k(),
                1e-12,
            );
            assert_eq!(film.pressure_pa(), bulk.pressure_pa());
            assert_ne!(film.density_kg_m3(), bulk.density_kg_m3());
            let heat = report.exchanged_heat_j;
            assert_eq!(heat > 0.0, enthalpy > 50000.0);
            near(
                heat,
                moving_input().duration_s
                    * report.conductance_w_k
                    * (sample.interface.body.temperature_k() - bulk.temperature_k()),
                super::super::exchange().source_tolerance_j,
            );
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
            let work = moving_input().exterior_pressure_pa
                * (next.gas().volume_m3() - old.gas().volume_m3());
            near(
                gas_energy + body_energy + kinetic + work,
                0.0,
                super::super::exchange().energy_tolerance_j,
            );
            near(
                body_energy + heat,
                0.0,
                super::super::exchange().energy_tolerance_j,
            );
            assert_eq!(candidate.thermal().external_heat_j, -heat);
            assert_eq!(next.gas().mass_kg(), old.gas().mass_kg());
            assert_eq!(
                next.gas().density_kg_m3(),
                next.gas().mass_kg() / next.gas().volume_m3()
            );
            assert_eq!(
                next.gas().state().acoustic_state().unwrap().density_kg_m3(),
                next.gas().density_kg_m3()
            );
            assert_eq!(old, saved);
        }
    }
}

#[test]
fn ambiguous_phase_density_policy_and_late_cancellation_cannot_publish_storage() {
    let curve = source::body_curve();
    let old = GasPiston::new(gas(350.0, 0.0), 2.0, 0.01, 0.15).unwrap();
    let saved = old.clone();
    let thermal = UniformEnthalpyStepInput {
        duration_s: moving_input().duration_s,
        ..body(&curve, 55000.0)
    };
    for (phase, law) in [
        (PengRobinsonPhaseV1::Stable, forced()),
        (PengRobinsonPhaseV1::Vapor, mixed(-0.25, 0.1)),
        (PengRobinsonPhaseV1::Vapor, mixed(0.25, 0.001)),
    ] {
        let source = PengRobinsonFilmSource {
            model: old.gas().model(),
            phase,
        };
        assert!(matches!(
            FilmGasSurface {
                model: &source,
                law
            }
            .finite_body_candidate(
                &old,
                moving_input(),
                thermal,
                moving_controls(),
                super::super::exchange(),
                || false,
            ),
            Err(ThermalPistonError::Interface(_))
        ));
        assert_eq!(old, saved);
    }
    let source = PengRobinsonFilmSource {
        model: old.gas().model(),
        phase: PengRobinsonPhaseV1::Vapor,
    };
    let film = FilmGasSurface {
        model: &source,
        law: forced(),
    };
    let mut total = 0;
    film.finite_body_candidate(
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
        film.finite_body_candidate(
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
