//! G0/G1/G3: conservation, exact reservoir relaxation and surface tangents.
use super::*;

#[derive(Clone, Copy)]
struct Law;
impl SublimationLaw for Law {
    fn temperature_range_k(&self) -> [f64; 2] {
        [50.0, 273.16]
    }
    fn properties(&self, t: f64) -> Result<SublimationProperties, SublimationError> {
        let p = 10.0 * (0.06 * (t - 230.0)).exp();
        Ok(SublimationProperties {
            saturation_pressure_pa: p,
            saturation_pressure_derivative_pa_k: 0.06 * p,
            latent_heat_j_kg: 2.8e6,
            latent_heat_derivative_j_kg_k: 0.0,
            solid_density_kg_m3: 900.0,
            molar_mass_kg_mol: 0.018,
        })
    }
}
fn kinetics() -> SurfaceKinetics {
    SurfaceKinetics {
        accommodation_coefficient: 1.0,
        mass_transfer_coefficient_m_s: None,
    }
}
fn bath(p: f64) -> VaporReservoir {
    VaporReservoir::Prescribed {
        partial_pressure_pa: p,
        temperature_k: 230.0,
    }
}
fn solid() -> SolidSurfaceState {
    SolidSurfaceState {
        mass_kg: 0.01,
        specific_enthalpy_j_kg: 460_000.0,
        area_m2: 0.01,
    }
}
fn near(a: f64, b: f64, atol: f64, rtol: f64) {
    assert!(
        (a - b).abs() <= atol + rtol * a.abs().max(b.abs()),
        "{a:e} != {b:e}"
    );
}
fn run(s: SolidSurfaceState, r: VaporReservoir, dt: f64) -> SurfaceStep {
    advance_surface(
        &Law,
        kinetics(),
        s,
        r,
        dt,
        SolverConfig::default(),
        |h| Ok(h / 2000.0),
        || Ok(()),
    )
    .unwrap()
}
#[test]
fn equilibrium_and_signed_deposition_share_all_flux_accounts() {
    assert_eq!(
        surface_flux(&Law, 230.0, kinetics(), bath(10.0))
            .unwrap()
            .mass_flux_kg_m2_s,
        0.0
    );
    let loss = surface_flux(&Law, 230.0, kinetics(), bath(0.0)).unwrap();
    let growth = surface_flux(&Law, 230.0, kinetics(), bath(20.0)).unwrap();
    assert!(loss.mass_flux_kg_m2_s > 0.0 && loss.outward_normal_velocity_m_s < 0.0);
    assert_eq!(loss.mass_flux_kg_m2_s, -growth.mass_flux_kg_m2_s);
    assert_eq!(loss.latent_heat_flux_w_m2, -growth.latent_heat_flux_w_m2);
    assert_eq!(
        loss.outward_normal_velocity_m_s,
        -growth.outward_normal_velocity_m_s
    );
}
#[test]
fn unequal_gas_temperature_uses_incoming_molecular_speed() {
    let gas = VaporReservoir::Prescribed {
        partial_pressure_pa: 10.0 * (260.0_f64 / 230.0).sqrt(),
        temperature_k: 260.0,
    };
    near(
        surface_flux(&Law, 230.0, kinetics(), gas)
            .unwrap()
            .mass_flux_kg_m2_s,
        0.0,
        1e-17,
        0.0,
    );
}
#[test]
fn analytic_surface_tangent_includes_latent_flux_and_kinetic_temperature() {
    let k = SurfaceKinetics {
        mass_transfer_coefficient_m_s: Some(0.02),
        ..kinetics()
    };
    let gas = VaporReservoir::Prescribed {
        partial_pressure_pa: 4.0,
        temperature_k: 240.0,
    };
    let f = surface_flux(&Law, 225.0, k, gas).unwrap();
    let plus = surface_flux(&Law, 225.0001, k, gas).unwrap();
    let minus = surface_flux(&Law, 224.9999, k, gas).unwrap();
    near(
        f.mass_flux_derivative_kg_m2_s_k,
        (plus.mass_flux_kg_m2_s - minus.mass_flux_kg_m2_s) / 0.0002,
        1e-15,
        1e-8,
    );
    near(
        f.latent_heat_flux_derivative_w_m2_k,
        (plus.latent_heat_flux_w_m2 - minus.latent_heat_flux_w_m2) / 0.0002,
        1e-10,
        1e-8,
    );
}
#[test]
fn film_resistance_recovers_diffusion_limit_and_sealed_boundaries() {
    let vacuum = surface_flux(&Law, 230.0, kinetics(), bath(0.0)).unwrap();
    let h = 1e-5;
    let film = SurfaceKinetics {
        mass_transfer_coefficient_m_s: Some(h),
        ..kinetics()
    };
    let actual = surface_flux(&Law, 230.0, film, bath(0.0))
        .unwrap()
        .mass_flux_kg_m2_s;
    let diff = h * 0.018 / (MOLAR_GAS_CONSTANT * 230.0) * 10.0;
    assert!(actual < diff && diff < vacuum.mass_flux_kg_m2_s);
    near(actual, diff, 0.0, 1e-6);
    for k in [
        SurfaceKinetics {
            accommodation_coefficient: 0.0,
            ..kinetics()
        },
        SurfaceKinetics {
            mass_transfer_coefficient_m_s: Some(0.0),
            ..kinetics()
        },
    ] {
        assert_eq!(
            surface_flux(&Law, 230.0, k, bath(0.0))
                .unwrap()
                .mass_flux_kg_m2_s,
            0.0
        );
    }
}
#[test]
fn finite_bath_exact_relaxation_never_crosses_equilibrium() {
    let volume = 0.001;
    let eq = 10.0 * volume * 0.018 / (MOLAR_GAS_CONSTANT * 230.0);
    for initial in [0.0, 0.5 * eq, 2.0 * eq] {
        let r = VaporReservoir::Finite {
            vapor_mass_kg: initial,
            volume_m3: volume,
            temperature_k: 230.0,
        };
        let dm = isothermal_mass_transfer(&Law, 230.0, kinetics(), solid(), r, 1e8).unwrap();
        near(initial + dm, eq, 1e-20, 1e-13);
        assert!(dm >= -initial && dm <= solid().mass_kg);
    }
}
#[test]
fn finite_bath_relaxation_has_semigroup_and_small_step_limit() {
    let r = VaporReservoir::Finite {
        vapor_mass_kg: 0.0,
        volume_m3: 1.0,
        temperature_k: 230.0,
    };
    let s = solid();
    let dt = 0.005;
    let first = isothermal_mass_transfer(&Law, 230.0, kinetics(), s, r, dt).unwrap();
    let next = VaporReservoir::Finite {
        vapor_mass_kg: first,
        volume_m3: 1.0,
        temperature_k: 230.0,
    };
    let second = isothermal_mass_transfer(&Law, 230.0, kinetics(), s, next, dt).unwrap();
    let all = isothermal_mass_transfer(&Law, 230.0, kinetics(), s, r, 2.0 * dt).unwrap();
    near(first + second, all, 1e-20, 1e-13);
    let tiny = isothermal_mass_transfer(&Law, 230.0, kinetics(), s, r, 1e-12).unwrap();
    let rate = surface_flux(&Law, 230.0, kinetics(), r)
        .unwrap()
        .mass_flux_kg_m2_s;
    near(tiny, rate * s.area_m2 * 1e-12, 1e-28, 1e-10);
}
#[test]
fn solid_inventory_limits_the_same_mass_that_drives_recession_and_energy() {
    let s = SolidSurfaceState {
        mass_kg: 1e-12,
        ..solid()
    };
    assert_eq!(
        isothermal_mass_transfer(&Law, 230.0, kinetics(), s, bath(0.0), 1e8).unwrap(),
        s.mass_kg
    );
}
#[test]
fn implicit_cooling_and_deposition_close_mass_energy_and_recession() {
    for pressure in [0.0, 20.0] {
        let s = solid();
        let step = run(s, bath(pressure), 2.0);
        let dm = step.transferred_mass_kg;
        near(step.solid.mass_kg + dm, s.mass_kg, 1e-17, 0.0);
        near(
            step.solid.mass_kg * step.solid.specific_enthalpy_j_kg
                - s.mass_kg * s.specific_enthalpy_j_kg
                + step.vapor_energy_j,
            0.0,
            1e-7,
            0.0,
        );
        near(step.latent_energy_j, dm * 2.8e6, 1e-12, 1e-14);
        near(step.recession_m, dm / s.area_m2 / 900.0, 1e-20, 1e-14);
        assert_eq!(step.surface_temperature_k < 230.0, pressure == 0.0);
        assert_eq!(dm > 0.0, pressure == 0.0);
        assert!(step.energy_residual_j.abs() < 1e-7);
    }
}
#[test]
fn coupled_finite_reservoir_conserves_total_mass_with_cooling() {
    let s = solid();
    let r = VaporReservoir::Finite {
        vapor_mass_kg: 0.0,
        volume_m3: 0.1,
        temperature_k: 230.0,
    };
    let step = run(s, r, 1e5);
    let VaporReservoir::Finite { vapor_mass_kg, .. } = step.reservoir else {
        panic!()
    };
    near(step.solid.mass_kg + vapor_mass_kg, s.mass_kg, 1e-17, 1e-14);
    assert!(step.surface_temperature_k < 230.0);
    let equilibrium = Law
        .properties(step.surface_temperature_k)
        .unwrap()
        .saturation_pressure_pa
        * (230.0 / step.surface_temperature_k).sqrt();
    near(
        step.reservoir.partial_pressure_pa(0.018).unwrap(),
        equilibrium,
        1e-10,
        1e-9,
    );
}
#[test]
fn zero_duration_area_and_equilibrium_are_noops() {
    for (s, r, dt) in [
        (solid(), bath(0.0), 0.0),
        (
            SolidSurfaceState {
                area_m2: 0.0,
                ..solid()
            },
            bath(0.0),
            10.0,
        ),
        (solid(), bath(10.0), 10.0),
    ] {
        let step = run(s, r, dt);
        assert_eq!(step.solid, s);
        assert_eq!(step.reservoir, r);
        assert_eq!(step.transferred_mass_kg, 0.0);
        assert_eq!(step.vapor_energy_j, 0.0);
    }
}
#[test]
fn caloric_datum_does_not_change_the_physical_endpoint() {
    let s = solid();
    let offset = 1e6;
    let baseline = run(s, bath(0.0), 1.0);
    let shifted = advance_surface(
        &Law,
        kinetics(),
        SolidSurfaceState {
            specific_enthalpy_j_kg: s.specific_enthalpy_j_kg + offset,
            ..s
        },
        bath(0.0),
        1.0,
        SolverConfig::default(),
        |h| Ok((h - offset) / 2000.0),
        || Ok(()),
    )
    .unwrap();
    near(
        baseline.surface_temperature_k,
        shifted.surface_temperature_k,
        1e-8,
        0.0,
    );
    near(
        baseline.transferred_mass_kg,
        shifted.transferred_mass_kg,
        1e-13,
        1e-8,
    );
    near(
        shifted.vapor_energy_j - baseline.vapor_energy_j,
        offset * baseline.transferred_mass_kg,
        1e-7,
        1e-8,
    );
}
#[test]
fn refinement_reduces_first_order_coupled_error() {
    fn evolve(n: usize) -> SolidSurfaceState {
        let mut s = solid();
        for _ in 0..n {
            s = run(s, bath(0.0), 1.0 / n as f64).solid;
        }
        s
    }
    let reference = evolve(512);
    let coarse = (evolve(8).specific_enthalpy_j_kg - reference.specific_enthalpy_j_kg).abs();
    let fine = (evolve(16).specific_enthalpy_j_kg - reference.specific_enthalpy_j_kg).abs();
    assert!(fine < 0.6 * coarse, "coarse {coarse}, fine {fine}");
}
#[test]
fn bad_inputs_cancellation_and_nonlinear_budget_refuse() {
    for bad in [f64::NAN, -1.0, f64::INFINITY] {
        assert!(surface_flux(&Law, 230.0, kinetics(), bath(bad)).is_err());
        assert!(
            isothermal_mass_transfer(&Law, 230.0, kinetics(), solid(), bath(0.0), bad).is_err()
        );
    }
    assert!(surface_flux(&Law, 280.0, kinetics(), bath(0.0)).is_err());
    assert_eq!(
        advance_surface(
            &Law,
            kinetics(),
            solid(),
            bath(0.0),
            10.0,
            SolverConfig::default(),
            |h| Ok(h / 2000.0),
            || Err(SublimationError::Cancelled)
        ),
        Err(SublimationError::Cancelled)
    );
    let config = SolverConfig {
        max_iterations: 1,
        ..SolverConfig::default()
    };
    assert_eq!(
        advance_surface(
            &Law,
            kinetics(),
            solid(),
            bath(0.0),
            10.0,
            config,
            |h| Ok(h / 2000.0),
            || Ok(())
        ),
        Err(SublimationError::NonConverged)
    );
}
