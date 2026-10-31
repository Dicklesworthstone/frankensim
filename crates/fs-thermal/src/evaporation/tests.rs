use super::*;

#[derive(Clone, Copy)]
struct Liquid {
    pressure_pa: f64,
    slope_pa_k: f64,
}
impl EvaporationLaw for Liquid {
    fn temperature_range_k(&self) -> [f64; 2] {
        [280.0, 330.0]
    }
    fn properties(&self, t: f64) -> Result<EvaporationProperties, EvaporationError> {
        Ok(EvaporationProperties {
            saturation_pressure_pa: self.pressure_pa + self.slope_pa_k * (t - 300.0),
            saturation_pressure_derivative_pa_k: self.slope_pa_k,
            latent_heat_j_kg: 2.4e6,
            liquid_density_kg_m3: 1000.0,
            molar_mass_kg_mol: 0.018,
        })
    }
}
fn law() -> Liquid {
    Liquid {
        pressure_pa: 4000.0,
        slope_pa_k: 100.0,
    }
}
fn film(pv: f64) -> GasFilm {
    GasFilm {
        total_pressure_pa: 100_000.0,
        vapor_partial_pressure_pa: pv,
        temperature_k: 300.0,
        mass_transfer_coefficient_m_s: 0.01,
    }
}
fn liquid() -> LiquidSurfaceState {
    LiquidSurfaceState {
        mass_kg: 0.1,
        specific_enthalpy_j_kg: 0.0,
        area_m2: 0.01,
    }
}
fn controls(offset: f64) -> SolverConfig {
    SolverConfig {
        specific_enthalpy_bracket_j_kg: [-80_000.0 + offset, 120_000.0 + offset],
        relative_tolerance: 1e-11,
        absolute_energy_tolerance_j: 1e-9,
        max_iterations: 128,
    }
}
fn near(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol, "{a} != {b}, tolerance {tol}");
}

#[test]
fn g1_stefan_flux_matches_integrated_binary_diffusion_and_dilute_limit() {
    let x = surface_flux(&law(), 300.0, film(1000.0)).unwrap();
    let independent =
        0.01 * 100_000.0 * 0.018 / (MOLAR_GAS_CONSTANT * 300.0) * (99_000.0_f64 / 96_000.0).ln();
    near(x.mass_flux_kg_m2_s, independent, 2e-18);
    near(x.latent_heat_flux_w_m2, independent * 2.4e6, 1e-11);
    near(x.outward_normal_velocity_m_s, -independent / 1000.0, 2e-21);
    let dilute = Liquid {
        pressure_pa: 0.004,
        slope_pa_k: 0.0,
    };
    let flux = surface_flux(&dilute, 300.0, film(0.001)).unwrap();
    let fick = 0.01 * 0.018 / (MOLAR_GAS_CONSTANT * 300.0) * 0.003;
    near(flux.mass_flux_kg_m2_s / fick, 1.0, 3e-8);
    // The nondilute correction must survive; dilute Fick would miss almost 39%.
    let concentrated = Liquid {
        pressure_pa: 50_000.0,
        slope_pa_k: 0.0,
    };
    let flux = surface_flux(&concentrated, 300.0, film(0.0)).unwrap();
    let fick = 0.01 * 0.018 / (MOLAR_GAS_CONSTANT * 300.0) * 50_000.0;
    near(flux.mass_flux_kg_m2_s / fick, 2.0 * 2.0_f64.ln(), 2e-15);
}

#[test]
fn g0_equilibrium_condensation_and_zero_transport_keep_physical_signs() {
    let equilibrium = surface_flux(&law(), 300.0, film(4000.0)).unwrap();
    assert_eq!(equilibrium.mass_flux_kg_m2_s, 0.0);
    assert_eq!(equilibrium.latent_heat_flux_w_m2, 0.0);
    let condensation = surface_flux(&law(), 300.0, film(6000.0)).unwrap();
    assert!(condensation.mass_flux_kg_m2_s < 0.0);
    assert!(condensation.latent_heat_flux_w_m2 < 0.0);
    assert!(condensation.outward_normal_velocity_m_s > 0.0);
    let sealed = GasFilm {
        mass_transfer_coefficient_m_s: 0.0,
        ..film(0.0)
    };
    assert_eq!(
        surface_flux(&law(), 300.0, sealed)
            .unwrap()
            .mass_flux_kg_m2_s,
        0.0
    );
    assert_eq!(
        surface_flux(&law(), 300.0, sealed)
            .unwrap()
            .mass_flux_derivative_kg_m2_s_k,
        0.0
    );
}

#[test]
fn g1_surface_temperature_tangent_retains_stefan_denominator() {
    for t in [281.0, 300.0, 329.0] {
        let delta = 1e-3;
        let a = surface_flux(&law(), t - delta, film(1000.0)).unwrap();
        let b = surface_flux(&law(), t + delta, film(1000.0)).unwrap();
        let derivative = (b.mass_flux_kg_m2_s - a.mass_flux_kg_m2_s) / (2.0 * delta);
        near(
            surface_flux(&law(), t, film(1000.0))
                .unwrap()
                .mass_flux_derivative_kg_m2_s_k,
            derivative,
            1e-14,
        );
    }
}

#[test]
fn g0_near_equilibrium_retains_gap_smaller_than_total_pressure_ulp() {
    let liquid = Liquid {
        pressure_pa: 0.001,
        slope_pa_k: 0.0,
    };
    let bath = GasFilm {
        total_pressure_pa: 1e12,
        vapor_partial_pressure_pa: f64::from_bits(0.001_f64.to_bits() - 1),
        ..film(0.0)
    };
    assert_eq!(
        bath.total_pressure_pa - liquid.pressure_pa,
        bath.total_pressure_pa - bath.vapor_partial_pressure_pa
    );
    let rate = surface_flux(&liquid, 300.0, bath)
        .unwrap()
        .mass_flux_kg_m2_s;
    assert!(rate > 0.0);
    let expected = 0.01 * 0.018 / (MOLAR_GAS_CONSTANT * 300.0)
        * (liquid.pressure_pa - bath.vapor_partial_pressure_pa);
    near(rate, expected, expected * 3e-15);
}

#[test]
fn g3_diffusivity_length_and_pressure_scaling_follow_one_flux_law() {
    let diffusion = GasFilm::from_diffusion(100_000.0, 1000.0, 300.0, 2.5e-5, 0.0025).unwrap();
    let base = surface_flux(&law(), 300.0, diffusion)
        .unwrap()
        .mass_flux_kg_m2_s;
    assert_eq!(diffusion, film(1000.0));
    let twice_pressure = GasFilm {
        total_pressure_pa: 200_000.0,
        vapor_partial_pressure_pa: 2000.0,
        mass_transfer_coefficient_m_s: 0.005,
        ..diffusion
    };
    let scaled = Liquid {
        pressure_pa: 8000.0,
        slope_pa_k: 200.0,
    };
    assert_eq!(
        base,
        surface_flux(&scaled, 300.0, twice_pressure)
            .unwrap()
            .mass_flux_kg_m2_s
    );
    assert!(GasFilm::from_diffusion(100_000.0, 0.0, 300.0, 2.5e-5, 0.0).is_err());
    assert!(GasFilm::from_diffusion(100_000.0, 0.0, 300.0, -1.0, 0.1).is_err());
}

#[test]
fn g1_implicit_latent_feedback_conserves_mass_and_full_stream_enthalpy() {
    for pv in [0.0, 6000.0] {
        let initial = liquid();
        let next = advance_surface(
            &law(),
            film(pv),
            initial,
            20.0,
            controls(0.0),
            |h| Ok(300.0 + h / 4000.0),
            || Ok(()),
        )
        .unwrap();
        assert_eq!(
            initial.mass_kg - next.liquid.mass_kg,
            next.transferred_mass_kg
        );
        assert!(next.energy_residual_j.abs() < 2e-8);
        near(
            next.liquid.mass_kg * next.liquid.specific_enthalpy_j_kg + next.vapor_energy_j,
            initial.mass_kg * initial.specific_enthalpy_j_kg,
            2e-8,
        );
        if pv == 0.0 {
            assert!(next.surface_temperature_k < 300.0);
            assert!(next.transferred_mass_kg > 0.0);
        } else {
            assert!(next.surface_temperature_k > 300.0);
            assert!(next.transferred_mass_kg < 0.0);
        }
        let rate = surface_flux(&law(), next.surface_temperature_k, film(pv))
            .unwrap()
            .mass_flux_kg_m2_s;
        near(
            next.transferred_mass_kg,
            20.0 * initial.area_m2 * rate,
            1e-17,
        );
        near(
            next.recession_m,
            next.transferred_mass_kg / initial.area_m2 / 1000.0,
            1e-20,
        );
    }
}

#[test]
fn g3_caloric_reference_shift_preserves_cooling_and_evaporation() {
    let run = |offset| {
        advance_surface(
            &law(),
            film(0.0),
            LiquidSurfaceState {
                specific_enthalpy_j_kg: offset,
                ..liquid()
            },
            20.0,
            controls(offset),
            |h| Ok(300.0 + (h - offset) / 4000.0),
            || Ok(()),
        )
        .unwrap()
    };
    let original = run(0.0);
    let shifted = run(1e7);
    near(
        original.surface_temperature_k,
        shifted.surface_temperature_k,
        5e-10,
    );
    near(
        original.transferred_mass_kg,
        shifted.transferred_mass_kg,
        2e-16,
    );
    near(
        shifted.vapor_energy_j - original.vapor_energy_j,
        1e7 * shifted.transferred_mass_kg,
        5e-8,
    );
}

#[test]
fn g1_backward_euler_converges_to_independently_integrated_cooling_trajectory() {
    // Independent explicit midpoint reference of the same continuous balance:
    // dm_l/dt=-A*j(T), m_l*cp*dT/dt=-A*j(T)*L.
    let mut reference_m = liquid().mass_kg;
    let mut reference_t = 300.0;
    let step = 0.001;
    let rate = |t: f64| {
        0.01 * 100_000.0 * 0.018 / (MOLAR_GAS_CONSTANT * 300.0)
            * (100_000.0 / (100_000.0 - (4000.0 + 100.0 * (t - 300.0)))).ln()
    };
    for _ in 0..20_000 {
        let r = rate(reference_t) * liquid().area_m2;
        let mid_m = reference_m - 0.5 * step * r;
        let mid_t = reference_t - 0.5 * step * r * 2.4e6 / (reference_m * 4000.0);
        let mid_r = rate(mid_t) * liquid().area_m2;
        reference_m -= step * mid_r;
        reference_t -= step * mid_r * 2.4e6 / (mid_m * 4000.0);
    }
    let run = |count: usize| {
        let mut state = liquid();
        let mut t = 300.0;
        for _ in 0..count {
            let next = advance_surface(
                &law(),
                film(0.0),
                state,
                20.0 / count as f64,
                controls(0.0),
                |h| Ok(300.0 + h / 4000.0),
                || Ok(()),
            )
            .unwrap();
            state = next.liquid;
            t = next.surface_temperature_k;
        }
        (state.mass_kg, t)
    };
    let coarse = run(10);
    let fine = run(20);
    assert!((fine.1 - reference_t).abs() < 0.55 * (coarse.1 - reference_t).abs());
    assert!((fine.0 - reference_m).abs() < 0.55 * (coarse.0 - reference_m).abs());
}

#[test]
fn g4_boiling_domains_inventory_and_cancellation_refuse_transactionally() {
    let initial = liquid();
    for bath in [
        GasFilm {
            total_pressure_pa: 4000.0,
            ..film(0.0)
        },
        GasFilm {
            total_pressure_pa: 3000.0,
            ..film(0.0)
        },
    ] {
        assert!(matches!(
            surface_flux(&law(), 300.0, bath),
            Err(EvaporationError::Boiling { .. })
        ));
    }
    assert!(surface_flux(&law(), 279.0, film(0.0)).is_err());
    assert!(surface_flux(&law(), f64::NAN, film(0.0)).is_err());
    assert!(surface_flux(&law(), 300.0, film(100_000.0)).is_err());
    let bounded = SolverConfig {
        specific_enthalpy_bracket_j_kg: [0.0, 120_000.0],
        ..controls(0.0)
    };
    assert!(matches!(
        advance_surface(
            &law(),
            film(0.0),
            initial,
            20.0,
            bounded,
            |h| Ok(300.0 + h / 4000.0),
            || Ok(())
        ),
        Err(EvaporationError::Unbracketed { .. })
    ));
    assert_eq!(
        advance_surface(
            &law(),
            film(0.0),
            initial,
            20.0,
            controls(0.0),
            |h| Ok(300.0 + h / 4000.0),
            || Err(EvaporationError::Cancelled)
        ),
        Err(EvaporationError::Cancelled)
    );
    let large = SolverConfig {
        specific_enthalpy_bracket_j_kg: [-1e9, 1e9],
        ..controls(0.0)
    };
    let Err(EvaporationError::InventoryLimit { maximum_duration_s }) = advance_surface(
        &law(),
        film(0.0),
        initial,
        1e6,
        large,
        |h| Ok(300.0 + h / 1e8),
        || Ok(()),
    ) else {
        panic!("exhausted liquid must refuse with an inventory retry estimate");
    };
    // The rejected endpoint is cooler/slower; its rate would overstate a retry
    // interval. The estimate must use the original, hotter evaporation rate.
    let initial_rate = surface_flux(&law(), 300.0, film(0.0))
        .unwrap()
        .mass_flux_kg_m2_s
        * initial.area_m2;
    near(maximum_duration_s * initial_rate, initial.mass_kg, 1e-15);
    assert_eq!(initial, liquid());
}
