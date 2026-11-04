use super::*;
use crate::water::if97::{If97Region, If97State};
use fs_qty::Pressure;

fn response(t: f64, rho: f64, cp: f64, rho_t: f64, rho_p: f64) -> DensityThermodynamicResponseV1 {
    DensityThermodynamicResponseV1::try_new(
        Temperature::new(t),
        Density::new(rho),
        MassSpecificThermalQuantityV1::new(cp),
        DensityTemperatureDerivativeV1::new(rho_t),
        DensityPressureDerivativeV1::new(rho_p),
    )
    .unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs(),
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

#[test]
fn g1_ideal_gas_density_derivatives_recover_calorics_and_adiabatic_sound() {
    for r in [2077.0, 287.05, 188.9] {
        for t in [250.0, 700.0, 1500.0] {
            for p in [1000.0, 101325.0, 5.0e6] {
                let cp = 3.5 * r;
                let rho = p / (r * t);
                let s = response(t, rho, cp, -rho / t, rho / p);
                close(s.volumetric_expansion_per_k(), 1.0 / t, 3e-15);
                close(s.isothermal_compressibility_per_pa(), 1.0 / p, 3e-15);
                close(s.heat_capacity_difference_j_kg_k(), r, 3e-15);
                close(s.specific_heat_cv_j_kg_k(), 2.5 * r, 3e-15);
                close(s.heat_capacity_ratio(), 1.4, 3e-15);
                close(s.adiabatic_bulk_modulus_pa(), 1.4 * p, 3e-15);
                close(s.sound_speed_m_s(), (1.4 * r * t).sqrt(), 3e-15);
            }
        }
    }
}

#[test]
fn g1_quadratic_gibbs_potential_matches_its_independent_hessian() {
    // g(T,p) = -C[T ln(T/T0)-T] + v0 dp + A dT dp - B dp²/2.
    // Its original Hessian gives Cp=C, v=v0+A dT-B dp,
    // Cv=C-T A²/B, and c²=v²/(B-T A²/C), independently of rho algebra.
    let (c, a, b) = (4000.0, 2e-7, 5e-13);
    for t in [280.0, 320.0, 400.0] {
        for p in [1e5, 1e7, 5e7] {
            let v = 0.001 + a * (t - 300.0) - b * (p - 1e5);
            let rho = 1.0 / v;
            let s = response(t, rho, c, -a / (v * v), b / (v * v));
            close(s.specific_heat_cv_j_kg_k(), c - t * a * a / b, 3e-15);
            close(s.sound_speed_m_s(), v / (b - t * a * a / c).sqrt(), 3e-15);
            close(s.volumetric_expansion_per_k(), a / v, 3e-15);
            close(
                s.adiabatic_compressibility_per_pa(),
                (b - t * a * a / c) / v,
                3e-15,
            );
        }
    }
}

#[test]
fn g1_original_if97_gibbs_derivatives_recover_liquid_and_steam_response() {
    for (p, t, region) in [
        (1e5, 274.0, If97Region::Region1), // negative thermal expansion
        (3e6, 300.0, If97Region::Region1),
        (80e6, 300.0, If97Region::Region1),
        (3e6, 500.0, If97Region::Region1),
        (3500.0, 300.0, If97Region::Region2),
        (30e6, 700.0, If97Region::Region2),
    ] {
        let original =
            If97State::at_pressure_temperature(Pressure::new(p), Temperature::new(t), region)
                .unwrap();
        let rho = original.density().value();
        let s = response(
            t,
            rho,
            original.cp().value(),
            -rho * rho * original.dv_dt_at_pressure().value(),
            -rho * rho * original.dv_dp_at_temperature().value(),
        );
        close(s.specific_heat_cv_j_kg_k(), original.cv().value(), 5e-14);
        close(s.sound_speed_m_s(), original.sound_speed().value(), 5e-14);
        close(
            s.volumetric_expansion_per_k(),
            original.thermal_expansion().value(),
            5e-14,
        );
        close(
            s.isothermal_compressibility_per_pa(),
            original.isothermal_compressibility().value(),
            5e-14,
        );
        if t == 274.0 {
            assert!(s.volumetric_expansion_per_k() < 0.0);
        }
        if t == 300.0 && p == 3e6 {
            // Original IAPWS R7-97(2012) Table 5 sound speed, published precision.
            assert!((s.sound_speed_m_s() - 1507.73921).abs() < 5e-6);
        }
    }
}

#[test]
fn g0_density_maximum_and_negative_expansion_are_not_gas_approximations() {
    let at_maximum = response(277.0, 1000.0, 4200.0, 0.0, 4.5e-7);
    assert_eq!(at_maximum.volumetric_expansion_per_k(), 0.0);
    assert_eq!(at_maximum.specific_heat_cv_j_kg_k(), 4200.0);
    assert_eq!(at_maximum.heat_capacity_difference_j_kg_k(), 0.0);
    assert_eq!(at_maximum.heat_capacity_ratio(), 1.0);
    assert_eq!(
        at_maximum.isothermal_compressibility_per_pa(),
        at_maximum.adiabatic_compressibility_per_pa()
    );
    close(at_maximum.sound_speed_m_s(), 1.0 / 4.5e-7_f64.sqrt(), 3e-15);
    let expanding = response(280.0, 1000.0, 4200.0, -0.1, 4.5e-7);
    let contracting = response(280.0, 1000.0, 4200.0, 0.1, 4.5e-7);
    assert_eq!(
        expanding.volumetric_expansion_per_k(),
        -contracting.volumetric_expansion_per_k()
    );
    assert_eq!(expanding.sound_speed_m_s(), contracting.sound_speed_m_s());
    assert_eq!(
        expanding.specific_heat_cv_j_kg_k(),
        contracting.specific_heat_cv_j_kg_k()
    );
    assert!(contracting.heat_capacity_difference_j_kg_k() > 0.0);
}

#[test]
fn g3_density_scaling_at_fixed_volumetric_heat_capacity_preserves_compressibility() {
    let reference = response(300.0, 1000.0, 4000.0, -0.2, 5e-7);
    // Keep rho*Cp fixed while changing density: the thermal pressure correction
    // and both compressibilities stay fixed, and sound scales as 1/sqrt(rho).
    for scale in [1e-6, 0.1, 1000.0, 1e6] {
        let s = response(
            300.0,
            1000.0 * scale,
            4000.0 / scale,
            -0.2 * scale,
            5e-7 * scale,
        );
        close(
            s.sound_speed_m_s(),
            reference.sound_speed_m_s() / scale.sqrt(),
            4e-15,
        );
        close(
            s.isothermal_compressibility_per_pa(),
            reference.isothermal_compressibility_per_pa(),
            3e-15,
        );
        close(
            s.specific_heat_cv_j_kg_k() * scale,
            reference.specific_heat_cv_j_kg_k(),
            4e-15,
        );
    }
}

#[test]
fn g0_local_thermodynamic_instability_and_invalid_inputs_refuse() {
    let run = |values: [f64; 5]| {
        DensityThermodynamicResponseV1::try_new(
            Temperature::new(values[0]),
            Density::new(values[1]),
            MassSpecificThermalQuantityV1::new(values[2]),
            DensityTemperatureDerivativeV1::new(values[3]),
            DensityPressureDerivativeV1::new(values[4]),
        )
    };
    let good = [300.0, 1000.0, 4000.0, -0.2, 5e-7];
    for axis in 0..5 {
        for invalid in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let mut values = good;
            values[axis] = invalid;
            assert!(run(values).is_err());
        }
    }
    for axis in [0, 1, 2, 4] {
        for invalid in [0.0, -1.0] {
            let mut values = good;
            values[axis] = invalid;
            assert!(run(values).is_err());
        }
    }
    // Positive rho, Cp and rho_p alone do not establish adiabatic stability.
    assert!(matches!(
        run([300.0, 1000.0, 10.0, -0.2, 5e-7]),
        Err(DensityResponseErrorV1::Unstable(_))
    ));
}

#[test]
fn g0_small_density_slopes_retain_their_caloric_correction() {
    // Squaring beta first would underflow and falsely return Cv=Cp.
    let s = response(1.0, 1.0, 2.0, -1e-170, 1e-300);
    close(s.heat_capacity_difference_j_kg_k(), 1e-40, 3e-15);
    assert_eq!(s.specific_heat_cv_j_kg_k(), 2.0);
    // bulk/rho overflows, but the sound speed itself is representable.
    let s = response(1.0, 1e-200, 1.0, 0.0, 1e-310);
    close(s.sound_speed_m_s(), 1e155, 1e-13);
}

#[test]
fn g0_roundoff_at_the_caloric_stability_boundary_cannot_create_sound() {
    // T beta²/rho_p = 24 J/(kg K). Square-root rounding may leave a tiny
    // positive Cp-delta; it is not evidence of a stable finite acoustic branch.
    for cp in [24.0, 24.0_f64.next_up(), 24.0_f64.next_down()] {
        assert!(
            DensityThermodynamicResponseV1::try_new(
                Temperature::new(300.0),
                Density::new(1000.0),
                MassSpecificThermalQuantityV1::new(cp),
                DensityTemperatureDerivativeV1::new(-0.2),
                DensityPressureDerivativeV1::new(5e-7),
            )
            .is_err()
        );
    }
    assert!(
        response(300.0, 1000.0, 24.1, -0.2, 5e-7)
            .sound_speed_m_s()
            .is_finite()
    );
}
