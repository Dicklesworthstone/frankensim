//! G0 identity/admission and G1 independent interface/film and energy balances.
use super::*;
use crate::evaporation::{LiquidSurfaceState, SolverConfig, advance_surface, surface_flux};

#[derive(Clone, Copy)]
struct Liquid {
    pressure: f64,
    slope: f64,
}
impl EvaporationLaw for Liquid {
    fn temperature_range_k(&self) -> [f64; 2] {
        [280.0, 330.0]
    }
    fn properties(&self, t: f64) -> Result<EvaporationProperties, EvaporationError> {
        Ok(EvaporationProperties {
            saturation_pressure_pa: self.pressure + self.slope * (t - 300.0),
            saturation_pressure_derivative_pa_k: self.slope,
            latent_heat_j_kg: 2.4e6,
            liquid_density_kg_m3: 1000.0,
            molar_mass_kg_mol: 0.018,
        })
    }
}
impl ThermalLaw for Liquid {
    fn storage_basis(&self) -> ThermalStorageBasis {
        ThermalStorageBasis::InternalEnergy
    }
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        Some([23; 32])
    }
    fn validate_cell(&self, _: &DepositedCell) -> Result<(), DepositionError> {
        Err(DepositionError::MaterialDomain)
    }
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        Ok(300.0 + h / 4000.0)
    }
    fn bulk_density_kg_m3(&self, _: f64) -> Result<Option<f64>, DepositionError> {
        Ok(Some(1000.0))
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / 4000.0
    }
    fn specific_entropy_change_j_kg_k(
        &self,
        a: f64,
        b: f64,
    ) -> Result<Option<f64>, DepositionError> {
        Ok(Some(
            4000.0 * ((300.0 + b / 4000.0) / (300.0 + a / 4000.0)).ln(),
        ))
    }
}
fn film(vapor: f64, transfer: f64) -> GasFilm {
    GasFilm {
        total_pressure_pa: 100_000.0,
        vapor_partial_pressure_pa: vapor,
        temperature_k: 300.0,
        mass_transfer_coefficient_m_s: transfer,
    }
}
fn hk(alpha: f64) -> InterfaceKinetics {
    InterfaceKinetics::HertzKnudsen {
        accommodation_coefficient: alpha,
    }
}
fn k(alpha: f64, t: f64) -> f64 {
    alpha * (0.018 / (2.0 * core::f64::consts::PI * MOLAR_GAS_CONSTANT * t)).sqrt()
}
fn near(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * actual.abs().max(expected.abs()),
        "{actual:e} != {expected:e}, relative tolerance {relative:e}"
    );
}

#[test]
fn g1_manufactured_interface_pressure_closes_both_nonlinear_fluxes() {
    // Choose p_i first, then derive h_m from independent kinetic and logarithmic
    // film equations. No production nonlinear solver supplies the oracle.
    for (ps, pv, pi) in [(50_000.0, 1000.0, 20_000.0), (3000.0, 50_000.0, 20_000.0)] {
        for schrage in [false, true] {
            let alpha = 0.003;
            let (model, factor) = if schrage {
                (
                    InterfaceKinetics::Schrage {
                        accommodation_coefficient: alpha,
                    },
                    2.0 * alpha / (2.0 - alpha),
                )
            } else {
                (hk(alpha), alpha)
            };
            let source = Liquid {
                pressure: ps,
                slope: 100.0,
            };
            let law = KineticEvaporation::try_new(&source, model).unwrap();
            let expected = k(factor, 300.0) * (ps - pi);
            let logarithm = ((100_000.0_f64 - pv) / (100_000.0 - pi)).ln();
            let transfer = expected * MOLAR_GAS_CONSTANT * 300.0 / (0.018 * 100_000.0 * logarithm);
            let actual = surface_flux(&law, 300.0, film(pv, transfer)).unwrap();
            near(actual.mass_flux_kg_m2_s, expected, 3e-14);
            near(actual.latent_heat_flux_w_m2, expected * 2.4e6, 3e-14);
            near(
                actual.outward_normal_velocity_m_s,
                -expected / 1000.0,
                3e-14,
            );
            let recovered_pi = ps - actual.mass_flux_kg_m2_s / k(factor, 300.0);
            near(recovered_pi, pi, 1e-13);
            let independent_film = transfer * 0.018 * 100_000.0 / (MOLAR_GAS_CONSTANT * 300.0)
                * ((100_000.0 - pv) / (100_000.0 - recovered_pi)).ln();
            near(independent_film, actual.mass_flux_kg_m2_s, 1e-13);
        }
    }
}

#[test]
fn g1_limits_preserve_film_and_kinetic_resistance() {
    let source = Liquid {
        pressure: 4000.0,
        slope: 100.0,
    };
    let fast = KineticEvaporation::try_new(&source, hk(1.0)).unwrap();
    let slow = KineticEvaporation::try_new(&source, hk(1e-12)).unwrap();
    for pv in [1000.0, 6000.0] {
        let diffusion = film(pv, 1e-10);
        let equilibrium = surface_flux(&source, 300.0, diffusion).unwrap();
        let limited = surface_flux(&fast, 300.0, diffusion).unwrap();
        near(
            limited.mass_flux_kg_m2_s,
            equilibrium.mass_flux_kg_m2_s,
            1e-11,
        );
        let molecular = surface_flux(&slow, 300.0, film(pv, 1.0)).unwrap();
        near(
            molecular.mass_flux_kg_m2_s,
            k(1e-12, 300.0) * (4000.0 - pv),
            3e-10,
        );
        assert!(
            molecular.mass_flux_kg_m2_s.abs()
                < surface_flux(&source, 300.0, film(pv, 1.0))
                    .unwrap()
                    .mass_flux_kg_m2_s
                    .abs()
        );
    }
}

#[test]
fn g1_temperature_tangent_includes_the_kinetic_prefactor_derivative() {
    for slope in [0.0, 500.0] {
        let source = Liquid {
            pressure: 50_000.0,
            slope,
        };
        for alpha in [1e-6, 1e-4, 1.0] {
            for model in [
                hk(alpha),
                InterfaceKinetics::Schrage {
                    accommodation_coefficient: alpha,
                },
            ] {
                let law = KineticEvaporation::try_new(&source, model).unwrap();
                for pv in [1000.0, 70_000.0] {
                    let bath = film(pv, 0.01);
                    let dt = 0.002;
                    let before = surface_flux(&law, 300.0 - dt, bath).unwrap();
                    let after = surface_flux(&law, 300.0 + dt, bath).unwrap();
                    let reference =
                        (after.mass_flux_kg_m2_s - before.mass_flux_kg_m2_s) / (2.0 * dt);
                    let exact = surface_flux(&law, 300.0, bath).unwrap();
                    near(exact.mass_flux_derivative_kg_m2_s_k, reference, 2e-5);
                    if slope == 0.0 {
                        assert_eq!(
                            exact.mass_flux_derivative_kg_m2_s_k.is_sign_negative(),
                            pv < 50_000.0
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn g0_near_equilibrium_retains_pressure_gap_below_total_pressure_ulp() {
    let source = Liquid {
        pressure: 0.001,
        slope: 0.0,
    };
    let alpha = 1e-5;
    let law = KineticEvaporation::try_new(&source, hk(alpha)).unwrap();
    for pv in [0.001_f64.next_down(), 0.001_f64.next_up()] {
        let bath = GasFilm {
            total_pressure_pa: 1e12,
            ..film(pv, 0.01)
        };
        assert_eq!(
            bath.total_pressure_pa - pv,
            bath.total_pressure_pa - source.pressure
        );
        let actual = surface_flux(&law, 300.0, bath).unwrap();
        let g = 0.01 * 0.018 / (MOLAR_GAS_CONSTANT * 300.0);
        let kinetic = k(alpha, 300.0);
        let expected = kinetic * g / (kinetic + g) * (source.pressure - pv);
        near(actual.mass_flux_kg_m2_s, expected, 2e-14);
        assert_ne!(actual.mass_flux_kg_m2_s, 0.0);
    }
}

#[test]
fn g0_equilibrium_and_seals_have_distinct_tangents() {
    let source = Liquid {
        pressure: 4000.0,
        slope: 100.0,
    };
    let alpha = 0.0001;
    let law = KineticEvaporation::try_new(&source, hk(alpha)).unwrap();
    let bath = film(4000.0, 0.01);
    let equilibrium = surface_flux(&law, 300.0, bath).unwrap();
    assert_eq!(equilibrium.mass_flux_kg_m2_s, 0.0);
    let g = 0.01 * 0.018 * 100_000.0 / (MOLAR_GAS_CONSTANT * 300.0 * 96_000.0);
    near(
        equilibrium.mass_flux_derivative_kg_m2_s_k,
        k(alpha, 300.0) * g / (k(alpha, 300.0) + g) * 100.0,
        2e-14,
    );
    for model in [
        hk(0.0),
        InterfaceKinetics::Schrage {
            accommodation_coefficient: 0.0,
        },
    ] {
        let sealed_law = KineticEvaporation::try_new(&source, model).unwrap();
        assert_eq!(
            surface_flux(&sealed_law, 300.0, film(0.0, 0.01)).unwrap(),
            sealed()
        );
    }
    assert_eq!(surface_flux(&law, 300.0, film(0.0, 0.0)).unwrap(), sealed());
}

#[test]
fn g0_default_interface_preserves_source_and_flux_bitwise() {
    let source = Liquid {
        pressure: 4000.0,
        slope: 100.0,
    };
    let wrapped = KineticEvaporation::try_new(&source, InterfaceKinetics::Equilibrium).unwrap();
    assert!(core::ptr::eq(wrapped.source(), &source));
    assert_eq!(wrapped.temperature_range_k(), source.temperature_range_k());
    for t in [280.0, 300.0, 330.0] {
        assert_eq!(wrapped.properties(t), source.properties(t));
        for pv in [0.0, 4000.0, 6000.0] {
            assert_eq!(
                surface_flux(&wrapped, t, film(pv, 0.01)),
                surface_flux(&source, t, film(pv, 0.01))
            );
        }
    }
}

#[test]
fn g0_invalid_kinetics_and_unsupported_domains_refuse() {
    let source = Liquid {
        pressure: 4000.0,
        slope: 100.0,
    };
    for alpha in [-1.0, 1.0_f64.next_up(), f64::NAN, f64::INFINITY] {
        assert!(KineticEvaporation::try_new(&source, hk(alpha)).is_err());
        assert!(
            KineticEvaporation::try_new(
                &source,
                InterfaceKinetics::Schrage {
                    accommodation_coefficient: alpha,
                }
            )
            .is_err()
        );
    }
    let law = KineticEvaporation::try_new(&source, hk(0.1)).unwrap();
    assert!(matches!(
        surface_flux(&law, 279.0, film(0.0, 0.01)),
        Err(EvaporationError::TemperatureOutOfRange { .. })
    ));
    assert!(matches!(
        surface_flux(
            &law,
            300.0,
            GasFilm {
                total_pressure_pa: 4000.0,
                ..film(0.0, 0.01)
            }
        ),
        Err(EvaporationError::Boiling { .. })
    ));
    assert!(surface_flux(&law, 300.0, film(100_000.0, 0.01)).is_err());
    struct Invalid;
    impl EvaporationLaw for Invalid {
        fn temperature_range_k(&self) -> [f64; 2] {
            [280.0, 330.0]
        }
        fn properties(&self, t: f64) -> Result<EvaporationProperties, EvaporationError> {
            Liquid {
                pressure: 4000.0,
                slope: 100.0,
            }
            .properties(t)
        }
        fn interface_kinetics(&self) -> InterfaceKinetics {
            hk(f64::NAN)
        }
    }
    assert!(surface_flux(&Invalid, 300.0, film(0.0, 0.0)).is_err());
}

#[test]
fn g1_implicit_latent_feedback_uses_the_same_finite_interface_rate() {
    let source = Liquid {
        pressure: 4000.0,
        slope: 100.0,
    };
    let law = KineticEvaporation::try_new(&source, hk(1e-6)).unwrap();
    for pv in [0.0, 6000.0] {
        let initial = LiquidSurfaceState {
            mass_kg: 0.1,
            specific_enthalpy_j_kg: 0.0,
            area_m2: 0.01,
        };
        let dt = 100.0;
        let next = advance_surface(
            &law,
            film(pv, 0.01),
            initial,
            dt,
            SolverConfig {
                specific_enthalpy_bracket_j_kg: [-80_000.0, 120_000.0],
                relative_tolerance: 1e-11,
                absolute_energy_tolerance_j: 1e-9,
                max_iterations: 128,
            },
            |h| Ok(300.0 + h / 4000.0),
            || Ok(()),
        )
        .unwrap();
        let rate = surface_flux(&law, next.surface_temperature_k, film(pv, 0.01)).unwrap();
        near(
            next.transferred_mass_kg,
            rate.mass_flux_kg_m2_s * initial.area_m2 * dt,
            1e-10,
        );
        assert_eq!(
            next.transferred_mass_kg,
            initial.mass_kg - next.liquid.mass_kg
        );
        assert!(next.energy_residual_j.abs() < 1e-8);
        assert!(
            (next.liquid.mass_kg * next.liquid.specific_enthalpy_j_kg + next.vapor_energy_j).abs()
                < 1e-8
        );
        assert_eq!(next.transferred_mass_kg.is_sign_positive(), pv < 4000.0);
    }
}

#[test]
fn g0_borrowed_wrapper_preserves_all_original_thermal_capabilities() {
    let source = Liquid {
        pressure: 4000.0,
        slope: 100.0,
    };
    let wrapper = KineticEvaporation::try_new(&source, hk(0.01)).unwrap();
    assert_eq!(wrapper.storage_basis(), source.storage_basis());
    assert_eq!(wrapper.caloric_identity(), source.caloric_identity());
    assert_eq!(wrapper.temperature_k(4000.0), source.temperature_k(4000.0));
    assert_eq!(
        wrapper.maximum_temperature_slope(),
        source.maximum_temperature_slope()
    );
    assert_eq!(
        wrapper.bulk_density_kg_m3(4000.0),
        source.bulk_density_kg_m3(4000.0)
    );
    assert_eq!(
        wrapper.specific_entropy_change_j_kg_k(0.0, 4000.0),
        source.specific_entropy_change_j_kg_k(0.0, 4000.0)
    );
    assert!(wrapper.advected_composition().is_none());
    assert!(wrapper.advected_entropy_anchor().is_none());
}
