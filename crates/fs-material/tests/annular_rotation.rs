//! G0/G1/G3 cylindrical momentum, work and independent constitutive controls.
use fs_material::MaterialError;
use fs_material::fluid_transport::rheology::annular::{
    CarreauAnnularControls, LaminarCarreauAnnulus,
};
use fs_material::fluid_transport::rheology::{
    CarreauYasuda, ThermalCarreauYasuda, ViscosityShiftMode,
};
use fs_material::visco::temperature::{TemperatureShift, TemperatureShiftModel};

fn state_error(what: &'static str) -> MaterialError {
    MaterialError::State { what: what.into() }
}

fn law(constant: bool, mode: ViscosityShiftMode) -> ThermalCarreauYasuda {
    let curve = if constant {
        [3.0, 3.0, 0.6, 0.5, 1.0]
    } else {
        [12.0, 0.0, 0.6, 0.5, 1.0]
    };
    ThermalCarreauYasuda::new(
        CarreauYasuda::new(curve[0], curve[1], curve[2], curve[3], curve[4], [0.0, 1e8]).unwrap(),
        TemperatureShift::new(
            TemperatureShiftModel::Arrhenius {
                activation_energy_j_mol: 12_000.0,
            },
            300.0,
            [280.0, 360.0],
        )
        .unwrap(),
        mode,
    )
}

fn owner(constant: bool) -> LaminarCarreauAnnulus {
    LaminarCarreauAnnulus::new(
        law(constant, ViscosityShiftMode::VerticalAndRate),
        300.0,
        0.01,
        0.04,
        0.5,
        1e5,
    )
    .unwrap()
}

fn tight() -> CarreauAnnularControls {
    let mut controls = CarreauAnnularControls::default();
    controls.local.relative_constitutive_tolerance = 1e-13;
    controls.local.relative_quadrature_tolerance = 2e-11;
    controls.relative_wall_velocity_tolerance = 1e-12;
    controls
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(actual.is_finite() && expected.is_finite());
    let scale = actual.abs().max(expected.abs()).max(1e-20);
    assert!(
        (actual - expected).abs() <= tolerance * scale,
        "actual {actual:e}, expected {expected:e}, tolerance {tolerance:e}"
    );
}

#[test]
fn g1_newtonian_rotation_recovers_cylindrical_torque_velocity_and_heat() {
    let annulus = owner(true);
    let (a, b, length, mu) = (0.01_f64, 0.04_f64, 0.5, 3.0);
    let resistance = 4.0 * core::f64::consts::PI * mu * length * a * a * b * b / (b * b - a * a);
    for walls in [[10.0, -2.0], [-7.0, 12.0], [5.0, 5.0], [0.0, 0.0]] {
        let profile = annulus.rotation_at_wall_speeds(walls, tight()).unwrap();
        let rotation = profile.rotation();
        let delta = walls[1] - walls[0];
        close(rotation.wall_torques_n_m[1], resistance * delta, 2e-14);
        assert_eq!(rotation.wall_torques_n_m[0], -rotation.wall_torques_n_m[1]);
        close(rotation.differential_resistance_n_m_s, resistance, 2e-14);
        close(
            rotation.dissipated_power_w,
            resistance * delta * delta,
            2e-14,
        );
        assert_eq!(rotation.wall_stress_iterations, 0);
        for fraction in [0.0, 0.1, 0.4, 0.8, 1.0] {
            let r = a + fraction * (b - a);
            let sample = profile.sample(r).unwrap();
            let angular = walls[0] + delta * (1.0 - a * a / (r * r)) / (1.0 - a * a / (b * b));
            let shear = 2.0 * delta * a * a * b * b / ((b * b - a * a) * r * r);
            assert!((sample.angular_velocity_rad_s - angular).abs() <= 1e-12);
            close(sample.azimuthal_velocity_m_s, r * angular, 1e-12);
            close(sample.shear_rate_per_s, shear, 2e-12);
            close(sample.radial_shear_stress_pa, mu * shear, 2e-12);
            close(sample.dissipation_w_m3, mu * shear * shear, 4e-12);
        }
    }
}

#[test]
fn g0_constant_rotation_admits_the_exact_source_endpoint_without_extrapolation() {
    let (a, b, maximum_rate) = (0.001_f64, 0.013_f64, 3.0_f64);
    let annulus = LaminarCarreauAnnulus::new(
        law(true, ViscosityShiftMode::VerticalAndRate),
        300.0,
        a,
        b,
        0.5,
        maximum_rate,
    )
    .unwrap();
    let span = ((b - a) / b) * (1.0 + a / b);
    let maximum_speed = 0.5 * span * maximum_rate;
    // Direct inversion overshoots the authored ceiling by one representable rate.
    assert!(2.0 * maximum_speed / span > maximum_rate);
    let endpoint = annulus
        .rotation_at_wall_speeds([0.0, maximum_speed], tight())
        .unwrap();
    assert_eq!(endpoint.rotation().maximum_shear_rate_per_s, maximum_rate);
    assert_eq!(endpoint.rotation().angular_velocity_residual_rad_s, 0.0);
    let above = f64::from_bits(maximum_speed.to_bits() + 1);
    assert!(
        annulus
            .rotation_at_wall_speeds([0.0, above], tight())
            .is_err()
    );
}

#[test]
fn g0_common_rigid_rotation_has_no_viscous_stress_or_heat() {
    // du_theta/dr is nonzero here. Using it as shear would falsely damp spin.
    let profile = owner(false)
        .rotation_at_wall_speeds([123.0, 123.0], tight())
        .unwrap();
    assert_eq!(profile.rotation().wall_power_w, 0.0);
    assert_eq!(profile.rotation().dissipated_power_w, 0.0);
    assert_eq!(profile.rotation().wall_torques_n_m, [-0.0, 0.0]);
    assert!(profile.rotation().differential_resistance_n_m_s > 0.0);
    for fraction in [0.0, 0.2, 0.5, 0.9, 1.0] {
        let sample = profile.sample_gap_fraction(fraction).unwrap();
        assert_eq!(sample.angular_velocity_rad_s, 123.0);
        assert_eq!(sample.shear_rate_per_s, 0.0);
        assert_eq!(sample.radial_shear_stress_pa, 0.0);
        assert_eq!(sample.dissipation_w_m3, 0.0);
    }
}

// Independent exact stress inverse for mu=12/sqrt(1+0.6*gamma).
fn inverse(stress: f64) -> f64 {
    let z = stress * 0.6 / 12.0;
    (z * z + z * (z * z + 4.0).sqrt()) / 1.2
}

fn simpson(f: impl Fn(f64) -> f64, a: f64, b: f64) -> f64 {
    let n = 8192;
    let mut sum = 0.0;
    for i in 0..=n {
        let weight = if i == 0 || i == n {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        sum += weight * f(a + (b - a) * f64::from(i) / f64::from(n));
    }
    (b - a) * sum / (3.0 * f64::from(n))
}

#[test]
fn g1_nonlinear_rotation_matches_independent_curved_stress_inverse_and_volume_heat() {
    let annulus = owner(false);
    let (a, b, length, inner_stress) = (0.01_f64, 0.04, 0.5, 60.0);
    let stress = |r: f64| inner_stress * (a / r).powi(2);
    let delta = simpson(|r| inverse(stress(r)) / r, a, b);
    let profile = annulus
        .rotation_at_wall_speeds([2.0, 2.0 + delta], tight())
        .unwrap();
    let rotation = profile.rotation();
    close(rotation.wall_tractions_pa[0], -inner_stress, 5e-10);
    close(rotation.relative_angular_velocity_rad_s, delta, 3e-10);
    let heat =
        core::f64::consts::TAU * length * simpson(|r| r * stress(r) * inverse(stress(r)), a, b);
    close(rotation.dissipated_power_w, heat, 6e-10);
    assert!(rotation.work_residual_w.abs() <= 3e-10 * heat);
    for r in [a, 0.013, 0.025, 0.037, b] {
        let sample = profile.sample(r).unwrap();
        close(sample.radial_shear_stress_pa, stress(r), 5e-10);
        close(sample.shear_rate_per_s, inverse(stress(r)), 1e-9);
        close(
            sample.angular_velocity_rad_s,
            2.0 + simpson(|s| inverse(stress(s)) / s, a, r),
            1e-9,
        );
        close(
            sample.dissipation_w_m3,
            stress(r) * inverse(stress(r)),
            2e-9,
        );
    }
    // Independent finite-difference torque derivative at a physical speed.
    let step = delta * 2e-5;
    let plus = annulus
        .rotation_at_wall_speeds([2.0, 2.0 + delta + step], tight())
        .unwrap();
    let minus = annulus
        .rotation_at_wall_speeds([2.0, 2.0 + delta - step], tight())
        .unwrap();
    let derivative =
        (plus.rotation().wall_torques_n_m[1] - minus.rotation().wall_torques_n_m[1]) / (2.0 * step);
    close(rotation.differential_resistance_n_m_s, derivative, 2e-7);
}

#[test]
fn g3_reversal_and_common_spin_preserve_the_same_irreversible_power() {
    let annulus = owner(false);
    let direct = annulus
        .rotation_at_wall_speeds([0.0, 13.0], tight())
        .unwrap();
    let shifted = annulus
        .rotation_at_wall_speeds([1e9, 1e9 + 13.0], tight())
        .unwrap();
    let reversed = annulus
        .rotation_at_wall_speeds([0.0, -13.0], tight())
        .unwrap();
    assert_eq!(
        direct.rotation().wall_torques_n_m,
        shifted.rotation().wall_torques_n_m
    );
    assert_eq!(
        direct.rotation().dissipated_power_w,
        shifted.rotation().dissipated_power_w
    );
    assert_eq!(
        direct.rotation().wall_power_w,
        shifted.rotation().wall_power_w
    );
    assert_eq!(
        direct.rotation().work_residual_w,
        shifted.rotation().work_residual_w
    );
    assert_eq!(
        direct.rotation().wall_torques_n_m.map(|v| -v),
        reversed.rotation().wall_torques_n_m
    );
    assert_eq!(
        direct.rotation().dissipated_power_w,
        reversed.rotation().dissipated_power_w
    );
    assert_eq!(
        direct.rotation().differential_resistance_n_m_s,
        reversed.rotation().differential_resistance_n_m_s
    );
}

#[test]
fn g1_narrow_gap_retains_nonlinear_shear_and_differential_resistance() {
    let a = 1.0;
    let b = a + 1e-9;
    let gap = b - a;
    let length = 0.5;
    let annulus = LaminarCarreauAnnulus::new(
        law(false, ViscosityShiftMode::VerticalAndRate),
        300.0,
        a,
        b,
        length,
        1e5,
    )
    .unwrap();
    let profile = annulus
        .rotation_at_wall_speeds([0.0, 2.0 * gap], tight())
        .unwrap();
    let local = annulus.law().evaluate(300.0, 2.0).unwrap();
    let sample = profile.sample_gap_fraction(0.5).unwrap();
    close(sample.shear_rate_per_s, 2.0, 1e-8);
    close(
        sample.dynamic_viscosity_pa_s,
        local.dynamic_viscosity_pa_s(),
        1e-8,
    );
    let resistance = core::f64::consts::TAU * length * local.stress_tangent_pa_s() / gap;
    close(
        profile.rotation().differential_resistance_n_m_s,
        resistance,
        1e-8,
    );
}

#[test]
fn g3_temperature_shift_is_retained_by_the_actual_rotational_consumer() {
    let law = law(false, ViscosityShiftMode::Vertical);
    let cold = LaminarCarreauAnnulus::new(law, 300.0, 0.01, 0.04, 0.5, 1e5).unwrap();
    let warm = LaminarCarreauAnnulus::new(law, 330.0, 0.01, 0.04, 0.5, 1e5).unwrap();
    let factor = law.temperature_shift().evaluate(330.0).unwrap().factor();
    let cold = cold.rotation_at_wall_speeds([1.0, 8.0], tight()).unwrap();
    let warm = warm.rotation_at_wall_speeds([1.0, 8.0], tight()).unwrap();
    close(
        warm.rotation().wall_torques_n_m[1],
        factor * cold.rotation().wall_torques_n_m[1],
        2e-9,
    );
    close(
        warm.rotation().dissipated_power_w,
        factor * cold.rotation().dissipated_power_w,
        2e-9,
    );
    close(
        warm.sample_gap_fraction(0.4)
            .unwrap()
            .angular_velocity_rad_s,
        cold.sample_gap_fraction(0.4)
            .unwrap()
            .angular_velocity_rad_s,
        2e-9,
    );
}

#[test]
fn g3_temperature_rate_shift_preserves_torque_tangent_and_local_heating() {
    let law = law(false, ViscosityShiftMode::VerticalAndRate);
    let reference = LaminarCarreauAnnulus::new(law, 300.0, 0.01, 0.04, 0.5, 1e5).unwrap();
    let warm = LaminarCarreauAnnulus::new(law, 330.0, 0.01, 0.04, 0.5, 1e5).unwrap();
    let factor = law.temperature_shift().evaluate(330.0).unwrap().factor();
    let reference = reference
        .rotation_at_wall_speeds([0.0, factor * 7.0], tight())
        .unwrap();
    let warm = warm.rotation_at_wall_speeds([0.0, 7.0], tight()).unwrap();
    close(
        warm.rotation().wall_torques_n_m[1],
        reference.rotation().wall_torques_n_m[1],
        2e-9,
    );
    close(
        warm.rotation().differential_resistance_n_m_s,
        factor * reference.rotation().differential_resistance_n_m_s,
        2e-9,
    );
    close(
        warm.rotation().dissipated_power_w,
        reference.rotation().dissipated_power_w / factor,
        2e-9,
    );
    for fraction in [0.0, 0.3, 0.8, 1.0] {
        let hot = warm.sample_gap_fraction(fraction).unwrap();
        let cold = reference.sample_gap_fraction(fraction).unwrap();
        close(
            hot.radial_shear_stress_pa,
            cold.radial_shear_stress_pa,
            2e-9,
        );
        close(hot.shear_rate_per_s, cold.shear_rate_per_s / factor, 2e-9);
        close(
            hot.angular_velocity_rad_s,
            cold.angular_velocity_rad_s / factor,
            2e-9,
        );
        close(
            hot.dynamic_viscosity_pa_s,
            factor * cold.dynamic_viscosity_pa_s,
            2e-9,
        );
        close(hot.dissipation_w_m3, cold.dissipation_w_m3 / factor, 2e-9);
    }
}

#[test]
fn g0_rotation_preserves_source_support_budgets_and_transactional_cancellation() {
    let annulus = owner(false);
    assert!(
        annulus
            .rotation_at_wall_speeds([0.0, 1e9], tight())
            .is_err()
    );
    assert!(
        annulus
            .rotation_at_wall_speeds([f64::NAN, 0.0], tight())
            .is_err()
    );
    let mut controls = tight();
    controls.local.maximum_quadrature_refinements = 2;
    assert!(
        annulus
            .rotation_at_wall_speeds([0.0, 13.0], controls)
            .is_err()
    );
    let mut calls = 0;
    let rejected = annulus.rotation_with_check([0.0, 13.0], tight(), || {
        calls += 1;
        if calls == 7 {
            Err(state_error("rotation interrupted"))
        } else {
            Ok(())
        }
    });
    assert!(rejected.is_err());
    assert_eq!(calls, 7);
    let retry = annulus
        .rotation_at_wall_speeds([0.0, 13.0], tight())
        .unwrap();
    assert_eq!(
        retry,
        annulus
            .rotation_at_wall_speeds([0.0, 13.0], tight())
            .unwrap()
    );
    assert!(retry.sample(f64::NAN).is_err());
    assert!(retry.sample_gap_fraction(-0.1).is_err());
    assert!(
        retry
            .sample_gap_fraction_with_check(0.5, || Err(state_error("sampling interrupted")))
            .is_err()
    );
}
