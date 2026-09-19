//! Analytical checks of the physical R-L-C boundary, not a fitted-filter golden.
use fs_vfit::impedance::{SeriesImpedance, SeriesImpedanceSpec};

fn spec(r: f64, l: f64, c: Option<f64>) -> SeriesImpedanceSpec {
    SeriesImpedanceSpec { resistance_pa_s_m3: r, inertance_pa_s2_m3: l, compliance_m3_pa: c }
}

#[test]
fn resistive_limit_matches_characteristic_reflection_without_storage() {
    for z in [0.5, 1.0, 1e6] {
        for r in [0.0, 0.5 * z, z, 5.0 * z] {
            let mut load = SeriesImpedance::new(spec(r, 0.0, None), z, 1e-3).unwrap();
            for a in [-20.0, 0.0, 7.0, 20.0] {
                let f = load.step(a).unwrap();
                let expected = a * (r - z) / (r + z);
                assert!((f.reflected_pressure_pa - expected).abs() < 1e-13 * (1.0 + a.abs()));
                assert_eq!(f.stored_energy_j, 0.0);
                assert_eq!(f.state.inertive_flow_m3_s, 0.0);
                assert_eq!(f.state.compliance_pressure_pa, 0.0);
                assert!(f.balance_residual_j().abs() < 1e-12 * (1.0 + f.supplied_work_j.abs()));
            }
        }
    }
}

#[test]
fn independent_storage_and_constitutive_equations_close_each_step() {
    for (r, l, c) in [(0.0, 0.3, None), (0.0, 0.0, Some(0.7)),
        (0.4, 0.3, Some(0.7)), (0.0, 0.3, Some(0.7))] {
        let mut load = SeriesImpedance::new(spec(r, l, c), 1.0, 0.02).unwrap();
        let mut released = false;
        for n in 0..2000 {
            let old = load.state();
            let before = load.stored_energy_j();
            let a = if n < 100 { f64::from(n % 13) - 6.0 } else { 0.0 };
            let f = load.step(a).unwrap();
            let u = f.flow_m3_s;
            let stored = 0.5 * l * f.state.inertive_flow_m3_s.powi(2)
                + 0.5 * c.unwrap_or(0.0) * f.state.compliance_pressure_pa.powi(2);
            let scale = (before + stored + f.supplied_work_j.abs() + f.dissipated_energy_j)
                .max(f64::MIN_POSITIVE);
            assert!((stored - f.stored_energy_j).abs() < 1e-12 * scale);
            assert!(f.balance_residual_j().abs() < 1e-12 * scale);
            assert!(f.dissipated_energy_j >= 0.0);
            let pressure = r * u + l * (f.state.inertive_flow_m3_s - old.inertive_flow_m3_s) / 0.02
                + 0.5 * (old.compliance_pressure_pa + f.state.compliance_pressure_pa);
            assert!((pressure - f.pressure_pa).abs() <= 1e-11 * (1.0 + pressure.abs()));
            if l > 0.0 {
                assert!((u - 0.5 * (old.inertive_flow_m3_s + f.state.inertive_flow_m3_s)).abs() < 1e-13);
            }
            if let Some(c) = c {
                assert!((c * (f.state.compliance_pressure_pa - old.compliance_pressure_pa) / 0.02 - u).abs() < 1e-13);
            }
            if n >= 100 {
                assert!(stored <= before + 1e-12 * scale);
                released |= f.supplied_work_j < 0.0;
            }
        }
        assert!(released, "reactive load must return stored energy, not classify it as active gain");
    }
}

#[test]
fn frequency_response_matches_the_independent_bilinear_impedance() {
    let (r, l, c, z, dt) = (0.4, 0.3, 0.7, 1.0, 0.02);
    for bin in [3, 19, 79, 200] {
        let theta = 2.0 * core::f64::consts::PI * f64::from(bin) / 2048.0;
        let mut load = SeriesImpedance::new(spec(r, l, Some(c)), z, dt).unwrap();
        let (mut re, mut im) = (0.0, 0.0);
        for n in 0..12288 {
            let phase = theta * f64::from(n);
            let f = load.step(phase.cos()).unwrap();
            if n >= 10240 {
                re += f.reflected_pressure_pa * phase.cos() / 1024.0;
                im -= f.reflected_pressure_pa * phase.sin() / 1024.0;
            }
        }
        let omega = 2.0 / dt * (0.5 * theta).tan();
        let x = omega * l - 1.0 / (omega * c);
        let den = (r + z).powi(2) + x * x;
        let expected_re = (r * r - z * z + x * x) / den;
        let expected_im = 2.0 * z * x / den;
        assert!((re - expected_re).hypot(im - expected_im) < 1e-10,
            "bin={bin}: ({re},{im}) != ({expected_re},{expected_im})");
    }
}

#[test]
fn preview_and_nonfinite_failures_preserve_reactive_state() {
    let mut a = SeriesImpedance::new(spec(0.4, 0.3, Some(0.7)), 1.0, 0.02).unwrap();
    a.step(20.0).unwrap();
    let mut b = a;
    let old = a.state();
    let before = a.stored_energy_j();
    let preview = a.preview_step(7.0).unwrap();
    for bad in [f64::NAN, f64::INFINITY, f64::MAX] { assert!(a.step(bad).is_err()); }
    assert_eq!(a.state(), old);
    assert_eq!(a.stored_energy_j(), before);
    assert_eq!(a.step(7.0).unwrap(), preview);
    assert_eq!(a.step(-2.0).unwrap(), { b.step(7.0).unwrap(); b.step(-2.0).unwrap() });
}

#[test]
fn admission_rejects_active_invalid_and_unrepresentable_elements() {
    for s in [spec(-1.0, 1.0, None), spec(1.0, -1.0, None), spec(f64::NAN, 1.0, None),
        spec(1.0, f64::INFINITY, None), spec(0.0, 1.0, Some(0.0)),
        spec(0.0, 1.0, Some(-1.0)), spec(0.0, 1.0, Some(f64::INFINITY))] {
        assert!(SeriesImpedance::new(s, 1.0, 0.02).is_err());
    }
    for dt in [0.0, -1.0, f64::NAN, f64::MIN_POSITIVE] {
        assert!(SeriesImpedance::new(spec(0.0, 10.0, None), 1.0, dt).is_err());
    }
}
