use super::*;
use crate::fluid_acoustics::{FluidAcousticState, thermoviscous::ThermoviscousFluidState};
use crate::fluid_transport::{FluidTransportState, ThermalFluidTransportState};

pub(super) fn medium(expansion: f64, viscosity: f64, bulk: f64) -> BulkThermoviscousEquilibrium {
    // rho=2, c=10, Cp=8, T0=12.5. Expansion is independently specified;
    // Cv follows the thermodynamic identity, not an ideal-gas assumption.
    let transport = FluidTransportState::try_new(2.0, 8.0, 0.6, viscosity).unwrap();
    let beta_s = 1.0 / 200.0;
    let beta_t = beta_s + expansion * expansion * 12.5 / 16.0;
    BulkThermoviscousEquilibrium::from_states(
        ThermoviscousFluidState::try_new(
            FluidAcousticState::try_new(2.0, 10.0).unwrap(),
            transport,
            8.0 * beta_s / beta_t,
        )
        .unwrap(),
        ThermalFluidTransportState::try_new(transport, 12.5, Some(expansion)).unwrap(),
        bulk,
    )
    .unwrap()
}

pub(super) fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-300),
        "actual={actual:e}, expected={expected:e}"
    );
}
pub(super) fn near_complex(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-300),
        "actual={actual:?}, expected={expected:?}"
    );
}

#[test]
fn g1_bulk_fields_and_absorption_close_the_original_nsf_equations() {
    for expansion in [-0.08, -0.04, 0.04, 0.08] {
        let equilibrium = medium(expansion, 0.003, 0.002);
        let cv = equilibrium.absorption().state().specific_heat_cv_j_kg_k();
        let gamma = 8.0 / cv;
        let beta_t = equilibrium.isothermal_compressibility_pa_inv();
        for omega in [0.01, 100.0, 10_000.0] {
            let prepared = equilibrium.prepare_wave_absorption(omega).unwrap();
            assert_eq!(prepared.equilibrium(), &equilibrium);
            assert_eq!(prepared.angular_frequency_rad_s(), omega);
            assert_eq!(
                prepared.wave(),
                equilibrium.absorption().thermoviscous_wave(omega).unwrap()
            );
            // Independent ORIGINAL dispersion quadratic, then the original
            // continuity, momentum, EOS and heat equations. This does not use
            // the production root's retained caloric/thermal expressions.
            let (x, y) = (omega * 0.003 / 100.0, omega * 0.0375 / 100.0);
            let a = C64::new(gamma * x * y, y);
            let b = C64::new(1.0, -(x + gamma * y));
            let z = C64::from_re(2.0) / (b + (b * b + a.scale(4.0)).sqrt());
            let k = z.sqrt().scale(omega / 10.0);
            let ik = C64::new(-k.im, k.re);
            let origin = C64::new(0.2, -0.3);
            let distance = 0.7 / k.abs();
            let (sin, cos) = det::sin_cos(k.re * distance);
            let u = origin * C64::new(cos, sin).scale(det::exp(-k.im * distance));
            let traction = C64::from_re(2.0 * omega) / k * u;
            let pressure = traction + ik.scale(0.006) * u;
            let temperature = -ik.scale(expansion * 12.5 / beta_t) * u
                / ((k * k).scale(0.6) + C64::new(0.0, -omega * 2.0 * cv));
            let heat_flux = -(ik * temperature).scale(0.6);
            let field = prepared.at_distance(origin, distance).unwrap();
            for (actual, expected) in [
                (field.velocity_m_s, u),
                (field.pressure_pa, pressure),
                (field.normal_traction_pa, traction),
                (field.temperature_perturbation_k, temperature),
                (field.conductive_heat_flux_w_m2, heat_flux),
            ] {
                near_complex(actual, expected, 2e-12);
            }
            near_complex(
                pressure.scale(beta_t) - temperature.scale(expansion),
                k * u / C64::from_re(omega),
                2e-12,
            );
            let viscous = 0.003 * (k * u).norm_sq();
            let entropy = 0.6 / (2.0 * 12.5 * 12.5) * (k * temperature).norm_sq();
            let mechanical_flux = 0.5 * (traction * u.conj()).re;
            let thermal_flux = 0.5 / 12.5 * (temperature * heat_flux.conj()).re;
            near(field.viscous_heating_w_m3, viscous, 3e-12);
            near(field.conductive_entropy_w_m3_k, entropy, 3e-12);
            near(field.thermal_availability_loss_w_m3, 12.5 * entropy, 3e-12);
            near(field.mechanical_power_flux_w_m2, mechanical_flux, 3e-12);
            // At weak loss the independent rounded phasors are almost in
            // quadrature: Re(T q*) is ill-conditioned even though both fields
            // and the physical flux are resolved. Bound that reference's
            // absolute roundoff; constitutive and energy checks stay tight.
            let flux_roundoff =
                16.0 * f64::EPSILON * temperature.abs() * heat_flux.abs() / (2.0 * 12.5);
            assert!(
                (field.thermal_availability_flux_w_m2 - thermal_flux).abs()
                    <= 3e-12 * thermal_flux.abs() + flux_roundoff
            );
            near(
                2.0 * k.im * field.acoustic_power_flux_w_m2(),
                viscous + 12.5 * entropy,
                4e-12,
            );
        }
    }
}

#[test]
fn g3_signed_expansion_and_zero_thermal_coupling_remain_physical() {
    let positive = medium(0.04, 0.003, 0.002)
        .prepare_wave_absorption(400.0)
        .unwrap();
    let negative = medium(-0.04, 0.003, 0.002)
        .prepare_wave_absorption(400.0)
        .unwrap();
    let a = positive.at_distance(C64::ONE, 0.02).unwrap();
    let b = negative.at_distance(C64::ONE, 0.02).unwrap();
    assert_eq!(positive.wave(), negative.wave());
    assert_eq!(a.temperature_perturbation_k, -b.temperature_perturbation_k);
    assert_eq!(a.conductive_heat_flux_w_m2, -b.conductive_heat_flux_w_m2);
    assert_eq!(a.conductive_entropy_w_m3_k, b.conductive_entropy_w_m3_k);
    assert_eq!(
        a.thermal_availability_flux_w_m2,
        b.thermal_availability_flux_w_m2
    );
    let zero = medium(0.0, 0.003, 0.002)
        .prepare_wave_absorption(400.0)
        .unwrap();
    let field = zero.at_distance(C64::ONE, 0.02).unwrap();
    assert_eq!(field.temperature_perturbation_k, C64::ZERO);
    assert_eq!(field.conductive_heat_flux_w_m2, C64::ZERO);
    assert_eq!(field.conductive_entropy_w_m3_k, 0.0);
    assert_eq!(field.thermal_availability_flux_w_m2, 0.0);
    near(
        2.0 * zero.wave().wavenumber_rad_m().im * field.mechanical_power_flux_w_m2,
        field.viscous_heating_w_m3,
        2e-14,
    );
}

#[test]
fn g1_thermal_field_survives_coincident_strong_diffusion_scales() {
    // nu_L=gamma*chi=0.075, gamma=2: at large omega the retained
    // thermal factor tends sqrt(2), even though q+i*gamma*y cancels if
    // reconstructed from rounded k. Pressure/velocity tends rho*c/sqrt(2).
    let prepared = medium(0.08, 0.075, 0.05)
        .prepare_wave_absorption(1e22 / 0.075)
        .unwrap();
    let field = prepared.at_distance(C64::from_re(1e-8), 0.0).unwrap();
    near_complex(
        field.temperature_perturbation_k,
        C64::from_re(1.25 * core::f64::consts::SQRT_2 * 1e-8),
        3e-10,
    );
    near_complex(
        field.pressure_pa,
        C64::from_re(20.0 / core::f64::consts::SQRT_2 * 1e-8),
        3e-10,
    );
    near(
        2.0 * prepared.wave().wavenumber_rad_m().im * field.acoustic_power_flux_w_m2(),
        field.viscous_heating_w_m3 + field.thermal_availability_loss_w_m3,
        3e-14,
    );
}

#[test]
fn g1_cell_integrals_retain_weak_loss_opaque_limits_and_partition_additivity() {
    let prepared = medium(0.08, 0.003, 0.002)
        .prepare_wave_absorption(100.0)
        .unwrap();
    let alpha = prepared.wave().wavenumber_rad_m().im;
    let origin = C64::new(0.3, -0.4);
    let field = prepared.at_distance(origin, 0.0).unwrap();
    let area = 0.7;
    let end = 2.0 / alpha;
    let full = prepared.absorption_between(origin, 0.0, end, area).unwrap();
    let mut sums = [0.0; 3];
    for i in 0..31 {
        let cell = prepared
            .absorption_between(
                origin,
                end * f64::from(i) / 31.0,
                end * f64::from(i + 1) / 31.0,
                area,
            )
            .unwrap();
        sums[0] += cell.viscous_heat_w;
        sums[1] += cell.conductive_entropy_w_k;
        sums[2] += cell.acoustic_power_loss_w;
    }
    for (actual, expected) in sums.into_iter().zip([
        full.viscous_heat_w,
        full.conductive_entropy_w_k,
        full.acoustic_power_loss_w,
    ]) {
        near(actual, expected, 2e-14);
    }
    let reference_length = (1.0 - (-4.0_f64).exp()) / (2.0 * alpha);
    near(
        full.viscous_heat_w,
        area * reference_length * field.viscous_heating_w_m3,
        3e-15,
    );
    assert!(full.balance_residual_w.abs() < 3e-14 * full.incoming_acoustic_power_w);
    near(
        full.thermal_availability_loss_w,
        12.5 * full.conductive_entropy_w_k,
        3e-15,
    );
    let thin = prepared
        .absorption_between(origin, 0.0, 1e-19 / alpha, area)
        .unwrap();
    assert_eq!(
        thin.incoming_acoustic_power_w,
        thin.outgoing_acoustic_power_w
    );
    assert!(thin.viscous_heat_w > 0.0 && thin.conductive_entropy_w_k > 0.0);
    near(
        thin.acoustic_power_loss_w,
        2e-19 * thin.incoming_acoustic_power_w,
        3e-14,
    );
    assert_eq!(thin.balance_residual_w, -thin.acoustic_power_loss_w);
    let opaque = prepared
        .absorption_between(origin, 0.0, 1000.0 / alpha, area)
        .unwrap();
    assert_eq!(opaque.outgoing_acoustic_power_w, 0.0);
    near(
        opaque.acoustic_power_loss_w,
        opaque.incoming_acoustic_power_w,
        3e-14,
    );
}

#[test]
fn g3_thermal_acoustic_flux_is_not_pressure_work_or_an_extra_heat_source() {
    let prepared = medium(0.08, 0.003, 0.002)
        .prepare_wave_absorption(10_000.0)
        .unwrap();
    let field = prepared.at_distance(C64::ONE, 0.0).unwrap();
    let pressure_work = 0.5 * (field.pressure_pa * field.velocity_m_s.conj()).re;
    assert!(field.acoustic_power_flux_w_m2() > 10.0 * pressure_work);
    assert!(field.thermal_availability_loss_w_m3 > 100.0 * field.viscous_heating_w_m3);
    assert_ne!(pressure_work, field.mechanical_power_flux_w_m2);
    let cell = prepared
        .absorption_between(C64::ONE, 0.0, 0.01, 0.7)
        .unwrap();
    near(
        cell.acoustic_power_loss_w,
        cell.viscous_heat_w + 12.5 * cell.conductive_entropy_w_k,
        3e-15,
    );
}

#[test]
fn g0_original_domain_refusals_and_zero_drives_do_not_manufacture_loss() {
    let equilibrium = medium(0.04, 0.003, 0.002);
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::MAX] {
        assert!(equilibrium.prepare_wave_absorption(bad).is_err());
    }
    let prepared = equilibrium.prepare_wave_absorption(100.0).unwrap();
    for bad in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(prepared.at_distance(C64::ONE, bad).is_err());
    }
    for bad in [C64::new(f64::NAN, 0.0), C64::new(1.0, f64::INFINITY)] {
        assert!(prepared.at_distance(bad, 0.0).is_err());
        assert!(prepared.absorption_between(bad, 0.0, 1.0, 1.0).is_err());
    }
    for (start, end, area) in [
        (1.0, 0.0, 1.0),
        (0.0, f64::NAN, 1.0),
        (0.0, 1.0, 0.0),
        (0.0, 1.0, f64::INFINITY),
    ] {
        assert!(
            prepared
                .absorption_between(C64::ONE, start, end, area)
                .is_err()
        );
    }
    let zero = prepared
        .absorption_between(C64::ZERO, 0.0, f64::MAX, 1.0)
        .unwrap();
    assert_eq!(zero.acoustic_power_loss_w, 0.0);
    assert_eq!(zero.balance_residual_w, 0.0);
    let empty = prepared
        .absorption_between(C64::ONE, 0.1, 0.1, 1.0)
        .unwrap();
    assert_eq!(empty.acoustic_power_loss_w, 0.0);
    assert_eq!(
        empty.incoming_acoustic_power_w,
        empty.outgoing_acoustic_power_w
    );
    assert_eq!(
        prepared
            .at_distance(C64::ONE, f64::MAX)
            .unwrap()
            .velocity_m_s,
        C64::ZERO
    );
}
