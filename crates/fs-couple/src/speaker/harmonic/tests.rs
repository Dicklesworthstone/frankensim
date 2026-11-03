use super::*;
use crate::speaker::declared_air;

fn budget() -> HarmonicBudget {
    HarmonicBudget {
        max_dimension: 5,
        max_dense_bytes: 32_768,
        max_dense_work: 10_000,
        relative_residual_tolerance: 2e-12,
    }
}
fn driver(panel: Option<PanelMode>) -> BaffledTsDriver {
    BaffledTsDriver::new_with_media(
        TsCard::datasheet_class_6p5(),
        Some(0.02),
        panel,
        declared_air(),
        Some(declared_air()),
    )
    .unwrap()
}
fn close(a: C64, b: C64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * a.abs().max(b.abs()).max(1e-20),
        "{a:?} != {b:?}"
    );
}

#[test]
fn full_frequency_load_matches_independent_electromechanical_circuit() {
    let card = TsCard::datasheet_class_6p5();
    let air = declared_air();
    let model = driver(None);
    let voltage = C64::new(0.2, -0.1);
    let stiffness =
        1.0 / card.cms_m_per_n + air.density * air.sound_speed.powi(2) * card.sd_m2.powi(2) / 0.02;
    for f in [5.0, 30.0, 80.0, 200.0, 1000.0, 6000.0, 18000.0] {
        let w = core::f64::consts::TAU * f;
        let response = model.response(w, voltage, budget(), || false).unwrap();
        let zm = C64::new(card.rms_n_s_m, stiffness / w - w * card.mms_kg)
            + response.radiation_impedance_n_s_m;
        let ze = C64::new(card.re_ohm, -w * card.le_h) + C64::from_re(card.bl.powi(2)) / zm;
        let current = voltage / ze;
        let velocity = current.scale(card.bl) / zm;
        close(response.coil_current_a, current, 3e-12);
        close(response.cone_velocity_m_s, velocity, 3e-12);
        close(
            response.cone_displacement_m,
            velocity / C64::new(0.0, -w),
            3e-12,
        );
        assert!(response.power_residual_w.abs() < 3e-12 * response.supplied_power_w);
        assert!(response.radiated_power_w > 0.0);
    }
}

#[test]
fn full_piston_retains_low_frequency_mass_and_frequency_squared_radiation() {
    let model = driver(None);
    let air = declared_air();
    let card = TsCard::datasheet_class_6p5();
    let radius = det::sqrt(card.sd_m2 / core::f64::consts::PI);
    let z_scale = air.density * air.sound_speed * card.sd_m2;
    let mut resistances = Vec::new();
    for ka in [1e-4, 2e-4] {
        let w = ka * air.sound_speed / radius;
        let response = model
            .response(w, C64::from_re(0.01), budget(), || false)
            .unwrap();
        let z = response.radiation_impedance_n_s_m;
        assert!((z.re / (0.5 * z_scale * ka * ka) - 1.0).abs() < 1e-7);
        let air_mass = 8.0 * air.density * radius.powi(3) / 3.0;
        assert!((-z.im / (w * air_mass) - 1.0).abs() < 1e-7);
        resistances.push(z.re);
    }
    assert!((resistances[1] / resistances[0] - 4.0).abs() < 1e-6);
    let high = model
        .response(
            50.0 * air.sound_speed / radius,
            C64::from_re(0.01),
            budget(),
            || false,
        )
        .unwrap();
    assert!((high.radiation_impedance_n_s_m.re / z_scale - 1.0).abs() < 0.002);
    assert!((high.radiation_impedance_n_s_m.im / z_scale).abs() < 0.02);
}

#[test]
fn cavity_panel_resonance_uses_original_coupling_and_closes_power() {
    let panel = PanelMode {
        frequency_hz: 180.0,
        mass_kg: 0.35,
        area_m2: 0.06,
        damping_ratio: 0.02,
    };
    let model = driver(Some(panel));
    let card = TsCard::datasheet_class_6p5();
    let air = declared_air();
    let bulk = air.density * air.sound_speed.powi(2) / 0.02;
    for f in [20.0, 150.0, 180.0, 250.0, 1000.0, 10_000.0] {
        let w = core::f64::consts::TAU * f;
        let response = model
            .response(w, C64::from_re(0.1), budget(), || false)
            .unwrap();
        let xp = response.panel_displacement_m.unwrap();
        let vp = response.panel_velocity_m_s.unwrap();
        let wp = core::f64::consts::TAU * panel.frequency_hz;
        let panel_dynamic = C64::new(
            panel.mass_kg * (wp * wp - w * w) + bulk * panel.area_m2.powi(2),
            -2.0 * w * panel.damping_ratio * panel.mass_kg * wp,
        );
        close(
            panel_dynamic * xp,
            response
                .cone_displacement_m
                .scale(bulk * panel.area_m2 * card.sd_m2),
            2e-11,
        );
        close(vp, C64::new(0.0, -w) * xp, 3e-12);
        assert!(response.panel_heat_w > 0.0);
        assert!(response.power_residual_w.abs() < 3e-12 * response.supplied_power_w);
    }
}

#[test]
fn undamped_panel_antiresonance_does_not_invert_a_singular_component() {
    let panel = PanelMode {
        frequency_hz: 180.0,
        mass_kg: 0.35,
        area_m2: 0.06,
        damping_ratio: 0.0,
    };
    let model = driver(Some(panel));
    let air = declared_air();
    let card = TsCard::datasheet_class_6p5();
    let wp = core::f64::consts::TAU * panel.frequency_hz;
    let bulk = air.density * air.sound_speed.powi(2) / 0.02;
    let w = det::sqrt(wp * wp + bulk * panel.area_m2.powi(2) / panel.mass_kg);
    let response = model
        .response(w, C64::from_re(0.1), budget(), || false)
        .unwrap();
    assert!(response.cone_displacement_m.abs() < 1e-17);
    assert!(response.panel_displacement_m.unwrap().abs() > 1e-7);
    close(
        response.coil_current_a,
        C64::from_re(0.1) / C64::new(card.re_ohm, -w * card.le_h),
        2e-12,
    );
    assert_eq!(response.panel_heat_w, 0.0);
}

#[test]
fn far_field_directivity_and_half_space_power_match_the_actual_radiation_load() {
    let model = driver(None);
    let radius = model.radius_m;
    // ka=6 includes the first directivity null and a signed sidelobe.
    let w = 6.0 * declared_air().sound_speed / radius;
    let response = model
        .response(w, C64::new(0.1, 0.03), budget(), || false)
        .unwrap();
    let range = 100.0;
    let on_axis = response
        .far_field_pressure(range, 0.0, 0.01, 0.001)
        .unwrap();
    let angle = (5.0_f64 / 6.0).asin();
    let sidelobe = response
        .far_field_pressure(range, angle, 0.01, 0.001)
        .unwrap();
    let ratio = sidelobe.pressure_pa / on_axis.pressure_pa;
    assert!(ratio.re < 0.0 && ratio.im.abs() < 1e-12);
    close(
        ratio,
        C64::from_re(2.0 * fs_math::special::bessel_j1(5.0) / 5.0),
        2e-12,
    );
    let two_range = response
        .far_field_pressure(2.0 * range, 0.0, 0.01, 0.001)
        .unwrap();
    let (s, c) = det::sin_cos(w * range / declared_air().sound_speed);
    close(
        two_range.pressure_pa / on_axis.pressure_pa,
        C64::new(c, s).scale(0.5),
        3e-12,
    );
    // Independent angular quadrature of intensity over the front hemisphere.
    let n = 2048;
    let dtheta = core::f64::consts::FRAC_PI_2 / n as f64;
    let mut power = 0.0;
    for i in 0..n {
        let theta = (i as f64 + 0.5) * dtheta;
        let pressure = response
            .far_field_pressure(range, theta, 0.01, 0.001)
            .unwrap()
            .pressure_pa;
        power += pressure.norm_sq() * det::sin(theta) * dtheta;
    }
    let air = declared_air();
    power *= core::f64::consts::PI * range * range / (air.density * air.sound_speed);
    assert!((power / response.radiated_power_w - 1.0).abs() < 2e-6);
}

#[test]
fn invalid_frequency_travel_budget_and_receiver_domains_refuse() {
    let model = driver(None);
    for w in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(model.response(w, C64::ONE, budget(), || false).is_err());
    }
    assert!(matches!(
        model.response(100.0, C64::ONE, budget(), || true),
        Err(BaffledDriverError::Harmonic(HarmonicError::Cancelled))
    ));
    let mut tiny = budget();
    tiny.max_dimension = 2;
    assert!(model.response(100.0, C64::ONE, tiny, || false).is_err());
    assert!(matches!(
        model.response(100.0, C64::from_re(1e6), budget(), || false),
        Err(BaffledDriverError::Source(
            SpeakerError::ExcursionExceeded { .. }
        ))
    ));
    let response = model
        .response(1000.0, C64::from_re(0.1), budget(), || false)
        .unwrap();
    assert!(
        response
            .far_field_pressure(model.radius_m * 0.5, 0.0, 0.1, 0.1)
            .is_err()
    );
    assert!(response.far_field_pressure(0.1, 0.0, 1e-6, 1e-6).is_err());
    assert!(
        response
            .far_field_pressure(100.0, core::f64::consts::PI, 0.1, 0.1)
            .is_err()
    );
    let changed = GasState::try_new(
        &fs_material::gas::GasSpec::dry_air_ussa1976(),
        293.15,
        202_650.0,
    )
    .unwrap();
    assert!(
        BaffledTsDriver::new_with_media(
            TsCard::datasheet_class_6p5(),
            Some(0.02),
            None,
            declared_air(),
            Some(changed)
        )
        .is_err()
    );
}
