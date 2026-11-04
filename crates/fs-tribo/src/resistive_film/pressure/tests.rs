use super::super::{FilmCell, FilmChannel, FilmLimits, GapPort};
use super::*;

fn gap(h: f64) -> GapPort {
    GapPort {
        reference_m: h,
        closure: vec![1.0],
    }
}

fn native() -> ResistiveFilm {
    ResistiveFilm::new(
        vec![FilmCell {
            area_m2: 0.002,
            gap: gap(0.001),
        }],
        vec![
            FilmChannel {
                from: 0,
                to: None,
                width_m: 0.02,
                length_m: 0.01,
                gap: gap(0.001),
            },
            FilmChannel {
                from: 0,
                to: None,
                width_m: 0.03,
                length_m: 0.02,
                gap: gap(0.001),
            },
        ],
        1,
        0.02,
        FilmLimits {
            minimum_cell_gap_m: 1e-6,
            maximum_gap_m: 0.01,
            maximum_pressure_pa: 1e7,
        },
    )
    .unwrap()
}

fn close(a: f64, b: f64, relative: f64, absolute: f64) {
    assert!(
        (a - b).abs() <= absolute + relative * a.abs().max(b.abs()),
        "{a:e} != {b:e}"
    );
}

fn evaluate(
    model: &PressureDrivenFilm,
    q: f64,
    v: f64,
) -> (PressureDrivenReport, f64, f64, [FilmChannelReport; 2]) {
    let (mut force, mut pressure, mut channels) = ([0.0], [0.0], [FilmChannelReport::default(); 2]);
    let report = model
        .evaluate_into(&[q], &[v], &mut force, &mut pressure, &mut channels)
        .unwrap();
    (report, force[0], pressure[0], channels)
}

#[test]
fn g1_external_pressure_drives_flow_and_supports_a_stationary_load() {
    let film = native();
    let g0 = 0.02 * 0.001_f64.powi(3) / (12.0 * 0.02 * 0.01);
    let g1 = 0.03 * 0.001_f64.powi(3) / (12.0 * 0.02 * 0.02);
    let model = PressureDrivenFilm::new(film, &[Some(2000.0), Some(-500.0)]).unwrap();
    let (report, force, p, channels) = evaluate(&model, 0.0, 0.0);
    let expected_p = (g0 * 2000.0 - g1 * 500.0) / (g0 + g1);
    let q = (2000.0 + 500.0) / (1.0 / g0 + 1.0 / g1);
    close(p, expected_p, 2e-14, 0.0);
    close(force, 0.002 * expected_p, 2e-14, 0.0);
    close(channels[0].volume_flow_m3_s, -q, 3e-14, 0.0);
    close(channels[1].volume_flow_m3_s, q, 3e-14, 0.0);
    close(report.reservoir_power_w, 2500.0 * q, 3e-14, 0.0);
    close(
        report.film.dissipated_power_w,
        report.reservoir_power_w,
        3e-14,
        0.0,
    );
    assert_eq!(report.mechanical_power_w, 0.0);
    assert!(force > 0.0);
}

#[test]
fn g1_pressure_can_deliver_mechanical_work_while_dissipation_stays_positive() {
    let model = PressureDrivenFilm::new(native(), &[Some(2000.0), Some(1000.0)]).unwrap();
    let (report, force, _, channels) = evaluate(&model, 0.0, -0.001);
    assert!(force > 0.0);
    assert!(report.mechanical_power_w < 0.0);
    assert!(report.reservoir_power_w > report.film.dissipated_power_w);
    close(
        report.mechanical_power_w + report.reservoir_power_w,
        report.film.dissipated_power_w,
        3e-14,
        0.0,
    );
    close(
        channels.iter().map(|c| c.volume_flow_m3_s).sum(),
        -0.002 * 0.001,
        3e-14,
        0.0,
    );
}

#[test]
fn g3_common_pressure_is_exact_at_rest_and_zero_supplies_recover_native_bits() {
    for alpha in [0.0, 2e-4] {
        let film = native().with_barus_viscosity(alpha).unwrap();
        let model = PressureDrivenFilm::new(film.clone(), &[Some(1234.5); 2]).unwrap();
        let (report, force, pressure, channels) = evaluate(&model, 0.0, 0.0);
        assert_eq!(pressure, 1234.5);
        assert_eq!(force, 0.002 * 1234.5);
        assert_eq!(channels, [FilmChannelReport::default(); 2]);
        assert_eq!(report.film.dissipated_power_w, 0.0);
        let zero = model.with_reservoir_pressures(&[Some(0.0); 2]).unwrap();
        let (a, force, pressure, channels) = evaluate(&zero, 0.0001, 0.003);
        let (mut f, mut p, mut c) = ([0.0], [0.0], [FilmChannelReport::default(); 2]);
        let b = film
            .evaluate_with_channels_into(&[0.0001], &[0.003], &mut f, &mut p, &mut c)
            .unwrap();
        assert_eq!(force.to_bits(), f[0].to_bits());
        assert_eq!(pressure.to_bits(), p[0].to_bits());
        assert_eq!(channels, c);
        close(a.film.dissipated_power_w, b.dissipated_power_w, 2e-15, 0.0);
    }
}

#[test]
fn g1_barus_reservoir_flux_uses_original_zero_gauge_viscosity() {
    let alpha = 2e-4;
    let model = PressureDrivenFilm::new(
        native().with_barus_viscosity(alpha).unwrap(),
        &[Some(2000.0), Some(-500.0)],
    )
    .unwrap();
    let (report, force, p, channels) = evaluate(&model, 0.0, 0.0);
    let g0 = 0.02 * 0.001_f64.powi(3) / (12.0 * 0.02 * 0.01);
    let g1 = 0.03 * 0.001_f64.powi(3) / (12.0 * 0.02 * 0.02);
    let psi = |p: f64| -(-alpha * p).exp_m1() / alpha;
    let target = (g0 * psi(2000.0) + g1 * psi(-500.0)) / (g0 + g1);
    let expected_p = -(1.0 - alpha * target).ln() / alpha;
    let q = (psi(2000.0) - psi(-500.0)) / (1.0 / g0 + 1.0 / g1);
    close(p, expected_p, 3e-14, 0.0);
    close(force, 0.002 * p, 2e-15, 0.0);
    close(channels[0].volume_flow_m3_s, -q, 3e-14, 0.0);
    close(channels[1].volume_flow_m3_s, q, 3e-14, 0.0);
    close(report.film.dissipated_power_w, q * 2500.0, 3e-14, 0.0);
}

#[test]
fn g1_full_tangent_includes_aperture_motion_and_supply_pressure() {
    for alpha in [0.0, 2e-4] {
        let model = PressureDrivenFilm::new(
            native().with_barus_viscosity(alpha).unwrap(),
            &[Some(2000.0), Some(-500.0)],
        )
        .unwrap();
        let (q, v, dq, dv) = (0.0001, 0.002, 0.0003, -0.001);
        let dp = [700.0, -350.0];
        let mut tangent = [0.0];
        model
            .tangent_into(&[q], &[v], &[dq], &[dv], &dp, &mut tangent)
            .unwrap();
        let eps = 1e-5;
        let plus = model
            .with_reservoir_pressures(&[Some(2000.0 + eps * dp[0]), Some(-500.0 + eps * dp[1])])
            .unwrap();
        let minus = model
            .with_reservoir_pressures(&[Some(2000.0 - eps * dp[0]), Some(-500.0 - eps * dp[1])])
            .unwrap();
        let fp = evaluate(&plus, q + eps * dq, v + eps * dv).1;
        let fm = evaluate(&minus, q - eps * dq, v - eps * dv).1;
        close(tangent[0], (fp - fm) / (2.0 * eps), 2e-8, 0.0);
    }
}

#[test]
fn g1_annular_hydrostatic_bearing_matches_radial_flow_and_load_refinement() {
    let (ri, ro, h, mu, ps): (f64, f64, f64, f64, f64) = (0.01, 0.03, 0.0001, 0.02, 100000.0);
    let log = (ro / ri).ln();
    let pi = std::f64::consts::PI;
    let exact_flow = 2.0 * pi * h.powi(3) * ps / (12.0 * mu * log);
    let exact_force = pi * ps * ((ro * ro - ri * ri) / (2.0 * log) - ri * ri);
    let mut errors = Vec::new();
    for n in [8, 16, 32] {
        let faces: Vec<f64> = (0..=n)
            .map(|i| ri * (log * i as f64 / n as f64).exp())
            .collect();
        let centers: Vec<f64> = (0..n).map(|i| (faces[i] * faces[i + 1]).sqrt()).collect();
        let cells = (0..n)
            .map(|i| FilmCell {
                area_m2: pi * (faces[i + 1].powi(2) - faces[i].powi(2)),
                gap: gap(h),
            })
            .collect();
        // Exact integrated radial Poiseuille resistance. width/length is
        // 2*pi/log(r2/r1), independent of the pressure solver under test.
        let channel = |from, to, a: f64, b: f64| FilmChannel {
            from,
            to,
            width_m: 2.0 * pi * (b - a) / (b / a).ln(),
            length_m: b - a,
            gap: gap(h),
        };
        let mut edges = vec![channel(0, None, ri, centers[0])];
        edges.extend((0..n - 1).map(|i| channel(i, Some(i + 1), centers[i], centers[i + 1])));
        edges.push(channel(n - 1, None, centers[n - 1], ro));
        let film = ResistiveFilm::new(
            cells,
            edges,
            1,
            mu,
            FilmLimits {
                minimum_cell_gap_m: 1e-6,
                maximum_gap_m: 0.001,
                maximum_pressure_pa: 2.0 * ps,
            },
        )
        .unwrap();
        let mut boundaries = vec![None; n + 1];
        boundaries[0] = Some(ps);
        boundaries[n] = Some(0.0);
        let model = PressureDrivenFilm::new(film, &boundaries).unwrap();
        let (mut force, mut pressure, mut channels) = (
            [0.0],
            vec![0.0; n],
            vec![FilmChannelReport::default(); n + 1],
        );
        let report = model
            .evaluate_into(&[0.0], &[0.0], &mut force, &mut pressure, &mut channels)
            .unwrap();
        close(channels[0].volume_flow_m3_s, -exact_flow, 3e-12, 0.0);
        close(channels[n].volume_flow_m3_s, exact_flow, 3e-12, 0.0);
        close(report.film.dissipated_power_w, ps * exact_flow, 3e-12, 0.0);
        for (i, p) in pressure.iter().enumerate() {
            close(*p, ps * (ro / centers[i]).ln() / log, 3e-12, 1e-8);
        }
        errors.push((force[0] - exact_force).abs());
    }
    assert!(
        errors[0] / errors[1] > 3.9 && errors[1] / errors[2] > 3.9,
        "{errors:?}"
    );
    assert!(errors[2] / exact_force < 1e-3);
}

#[test]
fn g0_closed_supply_and_failed_admission_publish_no_partial_outputs() {
    let mut film = native();
    std::sync::Arc::make_mut(&mut film.channels)[0].gap = gap(-0.001);
    let model = PressureDrivenFilm::new(film, &[Some(2000.0), Some(500.0)]).unwrap();
    let (report, _, p, channels) = evaluate(&model, 0.0, 0.0);
    close(p, 500.0, 2e-15, 0.0);
    assert_eq!(channels[0], FilmChannelReport::default());
    close(report.film.dissipated_power_w, 0.0, 0.0, 1e-25);
    for boundaries in [
        vec![Some(0.0)],
        vec![None, Some(0.0)],
        vec![Some(f64::NAN), Some(0.0)],
        vec![Some(2e7), Some(0.0)],
    ] {
        assert!(PressureDrivenFilm::new(native(), &boundaries).is_err());
    }
    let (mut force, mut pressure, mut channels) = (
        [9.0],
        [8.0],
        [FilmChannelReport {
            volume_flow_m3_s: 7.0,
            dissipated_power_w: 6.0,
        }; 2],
    );
    let saved = channels;
    assert!(
        model
            .evaluate_into(&[0.002], &[0.0], &mut force, &mut pressure, &mut channels)
            .is_err()
    );
    assert_eq!(force, [9.0]);
    assert_eq!(pressure, [8.0]);
    assert_eq!(channels, saved);
    assert!(
        model
            .tangent_into(&[0.0], &[0.0], &[0.0], &[0.0], &[f64::NAN, 0.0], &mut force)
            .is_err()
    );
    assert_eq!(force, [9.0]);
}
