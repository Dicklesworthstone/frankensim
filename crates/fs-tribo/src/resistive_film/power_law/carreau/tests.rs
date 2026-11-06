//! Physical controls and mathematical oracles, not a fluid calibration.
use super::*;
use crate::resistive_film::shear::{
    ShearChannelReport, ShearedFilm, pressure::PressureDrivenShearedFilm,
};
use crate::resistive_film::{FilmCell, FilmChannel, FilmLimits, GapPort};
use fs_material::fluid_transport::rheology::{
    CarreauYasuda, ThermalCarreauYasuda, ViscosityShiftMode,
};
use fs_material::visco::temperature::{TemperatureShift, TemperatureShiftModel};

const H: f64 = 2e-4;
const W: f64 = 0.01;
const L: f64 = 0.02;
const AREA: f64 = 1e-4;

fn near(a: f64, b: f64, absolute: f64, relative: f64) {
    assert!(a.is_finite() && b.is_finite());
    assert!(
        (a - b).abs() <= absolute + relative * a.abs().max(b.abs()),
        "{a:.16e} != {b:.16e}"
    );
}

fn controls() -> CarreauFilmControls {
    let mut c = CarreauFilmControls::default();
    c.pressure.absolute_flow_tolerance_m3_s = 1e-22;
    c.pressure.relative_flow_tolerance = 1e-10;
    c.pressure.absolute_power_tolerance_w = 1e-20;
    c.pressure.relative_power_tolerance = 2e-9;
    c.passage.local.relative_constitutive_tolerance = 1e-13;
    c.passage.local.relative_quadrature_tolerance = 1e-10;
    c.passage.relative_wall_velocity_tolerance = 1e-11;
    c
}

fn passage(mu: f64, thinning: bool, temperature: f64, width: f64) -> LaminarCarreauSlit {
    let law = ThermalCarreauYasuda::new(
        CarreauYasuda::new(
            mu,
            if thinning { 0.0 } else { mu },
            0.8,
            0.5,
            1.0,
            [0.0, 1e8],
        )
        .unwrap(),
        TemperatureShift::new(
            TemperatureShiftModel::Arrhenius {
                activation_energy_j_mol: 12_000.0,
            },
            300.0,
            [280.0, 360.0],
        )
        .unwrap(),
        ViscosityShiftMode::Vertical,
    );
    LaminarCarreauSlit::new(law, temperature, width, H, L, 1e7).unwrap()
}

fn base() -> ResistiveFilm {
    let gap = GapPort {
        reference_m: H,
        closure: vec![1.0, 0.0, 0.0],
    };
    ResistiveFilm::new(
        vec![FilmCell {
            area_m2: AREA,
            gap: gap.clone(),
        }],
        vec![
            FilmChannel {
                from: 0,
                to: None,
                width_m: W,
                length_m: L,
                gap
            };
            2
        ],
        3,
        0.02,
        FilmLimits {
            minimum_cell_gap_m: 1e-6,
            maximum_gap_m: 1e-3,
            maximum_pressure_pa: 1e7,
        },
    )
    .unwrap()
}

fn walls() -> Vec<WallVelocityPorts> {
    vec![
        WallVelocityPorts {
            lower: vec![0.0, -1.0, 0.0],
            upper: vec![0.0, 0.0, -1.0],
        },
        WallVelocityPorts {
            lower: vec![0.0, 1.0, 0.0],
            upper: vec![0.0, 0.0, 1.0],
        },
    ]
}

fn evaluate(
    film: &CarreauFilm,
    q: &[f64],
    velocity: &[f64],
) -> (
    Vec<f64>,
    Vec<f64>,
    Vec<CarreauFilmChannelReport>,
    CarreauFilmReport,
) {
    let mut f = vec![0.0; film.film.port_count()];
    let mut p = vec![0.0; film.film.cell_count()];
    let mut channels = vec![CarreauFilmChannelReport::default(); film.film.channel_count()];
    let report = film
        .evaluate_into(q, velocity, &mut f, &mut p, &mut channels)
        .unwrap();
    (f, p, channels, report)
}

#[test]
fn g1_newtonian_source_matches_complete_supplied_sliding_reynolds_owner() {
    let geometry = base();
    let supplies = [Some(1200.0), Some(-800.0)];
    let film = CarreauFilm::new(
        &geometry,
        vec![passage(0.02, false, 300.0, W); 2],
        &supplies,
        walls(),
        controls(),
    )
    .unwrap();
    let native = PressureDrivenShearedFilm::new(
        ShearedFilm::new(geometry.clone(), walls()).unwrap(),
        &supplies,
    )
    .unwrap();
    let q = [2e-5, 0.0, 0.0];
    let velocity = [2e-6, 0.013, 0.021];
    let (f, p, ch, report) = evaluate(&film, &q, &velocity);
    let mut nf = [0.0; 3];
    let mut np = [0.0];
    let mut nc = [ShearChannelReport::default(); 2];
    let nr = native
        .evaluate_into(&q, &velocity, &mut nf, &mut np, &mut nc)
        .unwrap();
    for i in 0..3 {
        near(f[i], nf[i], 1e-16, 2e-13);
    }
    near(p[0], np[0], 1e-12, 2e-13);
    for i in 0..2 {
        near(ch[i].volume_flow_m3_s, nc[i].volume_flow_m3_s, 1e-22, 2e-13);
        near(
            ch[i].dissipated_power_w,
            nc[i].dissipated_power_w(),
            1e-18,
            2e-13,
        );
    }
    near(
        report.film.dissipated_power_w,
        nr.film.dissipated_power_w,
        1e-18,
        2e-13,
    );
    near(
        report.mechanical_power_w + report.reservoir_power_w,
        report.film.dissipated_power_w,
        1e-17,
        3e-13,
    );
    assert_eq!(report.iterations, 0);
    assert_eq!(
        film.passages()[0].law(),
        passage(0.02, false, 300.0, W).law()
    );
}

fn simpson(f: impl Fn(f64) -> f64) -> f64 {
    let panels = 2048;
    let dx = 1.0 / panels as f64;
    let mut sum = f(0.0) + f(1.0);
    for i in 1..panels {
        sum += if i % 2 == 0 { 2.0 } else { 4.0 } * f(i as f64 * dx);
    }
    sum * dx / 3.0
}

// Independent closed-form inverse of tau=mu0*gamma/sqrt(1+0.8*|gamma|).
// The oracle integrates in physical gap position, unlike the production
// inverse-tangent moments, and solves only its own scalar no-slip condition.
fn inverse(tau: f64, mu: f64) -> f64 {
    let s = 0.8 * tau.abs() / mu;
    (0.5 * (s * s + s * (s * s + 4.0).sqrt()) / 0.8).copysign(tau)
}

fn oracle(dp: f64, mu: f64, h: f64, lower: f64, upper: f64) -> (f64, f64, [f64; 2]) {
    let half = dp * h / (2.0 * L);
    let relative = upper - lower;
    let (mut low, mut high) = (-1000.0, 1000.0);
    for _ in 0..90 {
        let c = 0.5 * (low + high);
        let r = h * simpson(|s| inverse(c + half * (1.0 - 2.0 * s), mu));
        if r < relative {
            low = c;
        } else {
            high = c;
        }
    }
    let c = 0.5 * (low + high);
    let q = W * h * (0.5 * lower + 0.5 * upper)
        + W * h * h * simpson(|s| (0.5 - s) * inverse(c + half * (1.0 - 2.0 * s), mu));
    let heat = W
        * L
        * h
        * simpson(|s| {
            let tau = c + half * (1.0 - 2.0 * s);
            tau * inverse(tau, mu)
        });
    (q, heat, [-c - half, c - half])
}

#[test]
fn g1_nonlinear_squeeze_supply_and_sliding_match_independent_gap_integrals() {
    let geometry = base();
    let supplies = [Some(100.0), Some(-200.0)];
    let film = CarreauFilm::new(
        &geometry,
        vec![passage(0.4, true, 300.0, W), passage(0.22, true, 300.0, W)],
        &supplies,
        walls(),
        controls(),
    )
    .unwrap();
    let expected_pressure = 500.0;
    let lower = 0.012;
    let upper = 0.028;
    let a = oracle(expected_pressure - 100.0, 0.4, H, -lower, -upper);
    let b = oracle(expected_pressure + 200.0, 0.22, H, lower, upper);
    let velocity = [(a.0 + b.0) / AREA, lower, upper];
    let (f, p, ch, r) = evaluate(&film, &[0.0; 3], &velocity);
    near(p[0], expected_pressure, 1e-6, 2e-8);
    near(f[0], AREA * expected_pressure, 1e-12, 2e-8);
    let mut expected_wall_force = [0.0; 2];
    for (i, original) in [a, b].into_iter().enumerate() {
        near(ch[i].volume_flow_m3_s, original.0, 1e-16, 2e-8);
        near(ch[i].dissipated_power_w, original.1, 1e-12, 3e-8);
        for t in 0..2 {
            near(ch[i].wall_tractions_pa[t], original.2[t], 1e-8, 3e-8);
            expected_wall_force[t] += if i == 0 { -1.0 } else { 1.0 } * W * L * original.2[t];
        }
    }
    for t in 0..2 {
        near(f[t + 1], expected_wall_force[t], 1e-11, 3e-8);
    }
    near(
        ch[0].volume_flow_m3_s + ch[1].volume_flow_m3_s,
        AREA * velocity[0],
        1e-18,
        1e-10,
    );
    near(r.film.dissipated_power_w, a.1 + b.1, 1e-12, 3e-8);
    near(
        r.mechanical_power_w + r.reservoir_power_w,
        r.film.dissipated_power_w,
        1e-15,
        3e-9,
    );
    assert!(r.iterations > 0);
    assert!(r.constitutive_work_defect_w.abs() < 1e-12);
}

#[test]
fn g1_complete_geometry_velocity_and_supply_tangent_matches_independent_differences() {
    let supplies = [Some(130.0), Some(-70.0)];
    let film = CarreauFilm::new(
        &base(),
        vec![passage(0.4, true, 300.0, W); 2],
        &supplies,
        walls(),
        controls(),
    )
    .unwrap();
    let q = [2e-5, 0.0, 0.0];
    let v = [2e-6, 0.009, 0.02];
    let dq = [2e-6, 0.0, 0.0];
    let dv = [2e-7, -0.001, 0.002];
    let dp = [30.0, -20.0];
    let mut tangent = [0.0; 3];
    film.tangent_into(&q, &v, &dq, &dv, &dp, &mut tangent)
        .unwrap();
    let epsilon = 2e-4;
    let shifted = |sign: f64| {
        let qq = core::array::from_fn::<_, 3, _>(|i| q[i] + sign * epsilon * dq[i]);
        let vv = core::array::from_fn::<_, 3, _>(|i| v[i] + sign * epsilon * dv[i]);
        let pp = core::array::from_fn::<_, 2, _>(|i| {
            Some(supplies[i].unwrap() + sign * epsilon * dp[i])
        });
        evaluate(&film.with_reservoir_pressures(&pp).unwrap(), &qq, &vv).0
    };
    let minus = shifted(-1.0);
    let plus = shifted(1.0);
    for i in 0..3 {
        near(
            tangent[i],
            (plus[i] - minus[i]) / (2.0 * epsilon),
            1e-9,
            2e-5,
        );
    }
    // The velocity-only Schur complement is symmetric because pressure and
    // wall traction use the same original inverse-stress tangent.
    let mut matrix = [[0.0; 3]; 3];
    for j in 0..3 {
        let mut direction = [0.0; 3];
        direction[j] = 1.0;
        film.tangent_into(&q, &v, &[0.0; 3], &direction, &[0.0; 2], &mut matrix[j])
            .unwrap();
    }
    for i in 0..3 {
        for j in 0..3 {
            near(matrix[i][j], matrix[j][i], 1e-10, 2e-8);
        }
    }
}

#[test]
fn g3_multicell_continuity_full_work_and_reversal_use_one_material_profile() {
    let gap = |row: Vec<f64>| GapPort {
        reference_m: H,
        closure: row,
    };
    let geometry = ResistiveFilm::new(
        vec![
            FilmCell {
                area_m2: AREA,
                gap: gap(vec![1.0, 0.0, 0.0, 0.0]),
            },
            FilmCell {
                area_m2: 0.8 * AREA,
                gap: gap(vec![0.0, 1.0, 0.0, 0.0]),
            },
        ],
        vec![
            FilmChannel {
                from: 0,
                to: None,
                width_m: W,
                length_m: L,
                gap: gap(vec![1.0, 0.0, 0.0, 0.0]),
            },
            FilmChannel {
                from: 0,
                to: Some(1),
                width_m: 0.8 * W,
                length_m: L,
                gap: gap(vec![0.5, 0.5, 0.0, 0.0]),
            },
            FilmChannel {
                from: 1,
                to: None,
                width_m: W,
                length_m: L,
                gap: gap(vec![0.0, 1.0, 0.0, 0.0]),
            },
        ],
        4,
        0.1,
        FilmLimits {
            minimum_cell_gap_m: 1e-6,
            maximum_gap_m: 1e-3,
            maximum_pressure_pa: 1e7,
        },
    )
    .unwrap();
    let ws = [-1.0, 1.0, 1.0]
        .into_iter()
        .map(|s| WallVelocityPorts {
            lower: vec![0.0, 0.0, s, 0.0],
            upper: vec![0.0, 0.0, 0.0, s],
        })
        .collect();
    let film = CarreauFilm::new(
        &geometry,
        vec![
            passage(0.4, true, 300.0, W),
            passage(0.3, true, 300.0, 0.8 * W),
            passage(0.2, true, 300.0, W),
        ],
        &[Some(0.0), None, Some(0.0)],
        ws,
        controls(),
    )
    .unwrap();
    let q = [1e-5, -2e-5, 0.0, 0.0];
    let v = [1e-6, -3e-7, 0.01, -0.003];
    let (f, p, ch, r) = evaluate(&film, &q, &v);
    near(
        ch[0].volume_flow_m3_s + ch[1].volume_flow_m3_s,
        AREA * v[0],
        1e-17,
        2e-8,
    );
    near(
        ch[2].volume_flow_m3_s - ch[1].volume_flow_m3_s,
        0.8 * AREA * v[1],
        1e-17,
        2e-8,
    );
    near(
        r.mechanical_power_w + r.reservoir_power_w,
        r.film.dissipated_power_w,
        1e-15,
        3e-9,
    );
    let reversed = core::array::from_fn::<_, 4, _>(|i| -v[i]);
    let (rf, rp, rc, rr) = evaluate(&film, &q, &reversed);
    for i in 0..4 {
        near(rf[i], -f[i], 1e-12, 3e-8);
    }
    for i in 0..2 {
        near(rp[i], -p[i], 1e-7, 3e-8);
    }
    for i in 0..3 {
        near(rc[i].volume_flow_m3_s, -ch[i].volume_flow_m3_s, 1e-18, 3e-8);
    }
    near(
        rr.film.dissipated_power_w,
        r.film.dissipated_power_w,
        1e-14,
        3e-8,
    );
}

#[test]
fn g3_thermal_source_refresh_changes_the_actual_graph_response() {
    let cold = CarreauFilm::new(
        &base(),
        vec![passage(0.4, true, 300.0, W); 2],
        &[Some(0.0); 2],
        walls(),
        controls(),
    )
    .unwrap();
    let hot = cold
        .with_passages(vec![passage(0.4, true, 330.0, W); 2])
        .unwrap();
    let v = [2e-5, 0.011, 0.023];
    let (_, cp, _, cr) = evaluate(&cold, &[0.0; 3], &v);
    let (_, hp, _, hr) = evaluate(&hot, &[0.0; 3], &v);
    assert!(cp[0] > hp[0] && hp[0] > 0.0);
    assert!(cr.film.dissipated_power_w > hr.film.dissipated_power_w);
    for (owner, t) in [(&cold, 300.0), (&hot, 330.0)] {
        for (p, mu) in owner
            .passages()
            .iter()
            .zip(owner.pressure_operator().channel_viscosities_pa_s())
        {
            assert_eq!(p.temperature_k(), t);
            assert_eq!(*mu, p.zero_shear_viscosity_pa_s());
            assert_eq!(p.law(), cold.passages()[0].law());
        }
    }
    let (_, cp2, _, cr2) = evaluate(&cold, &[0.0; 3], &v);
    assert_eq!(cp2, cp);
    assert_eq!(cr2.film.dissipated_power_w, cr.film.dissipated_power_w);
}

#[test]
fn g3_common_translation_does_not_erase_small_sliding_or_create_heat() {
    let ws = [-1.0, 1.0]
        .into_iter()
        .map(|s| WallVelocityPorts {
            lower: vec![0.0, s, -0.5 * s],
            upper: vec![0.0, s, 0.5 * s],
        })
        .collect();
    let film = CarreauFilm::new(
        &base(),
        vec![passage(0.4, true, 300.0, W); 2],
        &[Some(0.0); 2],
        ws,
        controls(),
    )
    .unwrap();
    let (f, p, ch, r) = evaluate(&film, &[0.0; 3], &[0.0, 0.0, 1e-7]);
    let (ff, pp, cc, rr) = evaluate(&film, &[0.0; 3], &[0.0, 1e12, 1e-7]);
    assert_eq!(ff, f);
    assert_eq!(pp, p);
    assert_eq!(rr.film.dissipated_power_w, r.film.dissipated_power_w);
    assert!(r.film.dissipated_power_w > 0.0);
    for i in 0..2 {
        assert_eq!(cc[i].pressure_flow_m3_s, 0.0);
        assert_eq!(cc[i].wall_tractions_pa, ch[i].wall_tractions_pa);
        assert_eq!(cc[i].dissipated_power_w, ch[i].dissipated_power_w);
        assert!(cc[i].volume_flow_m3_s.abs() > 1e6);
    }
}

#[test]
fn g4_source_limits_and_late_cancellation_are_atomic_and_exactly_retryable() {
    let film = CarreauFilm::new(
        &base(),
        vec![passage(0.4, true, 300.0, W); 2],
        &[Some(10.0); 2],
        walls(),
        controls(),
    )
    .unwrap();
    let q = [0.0; 3];
    let v = [2e-5, 0.01, 0.03];
    let (expected_f, expected_p, expected_ch, expected_report) = evaluate(&film, &q, &v);
    let mut calls = 0;
    let mut f = [0.0; 3];
    let mut p = [0.0];
    let mut ch = [CarreauFilmChannelReport::default(); 2];
    film.evaluate_with_check_into(&q, &v, &mut f, &mut p, &mut ch, || {
        calls += 1;
        Ok(())
    })
    .unwrap();
    assert!(calls > 20);
    let sentinel = CarreauFilmChannelReport {
        volume_flow_m3_s: 123.0,
        ..Default::default()
    };
    for stop in [calls / 2, calls] {
        let mut f = [77.0; 3];
        let mut p = [88.0];
        let mut ch = [sentinel; 2];
        let mut n = 0;
        let result = film.evaluate_with_check_into(&q, &v, &mut f, &mut p, &mut ch, || {
            n += 1;
            if n == stop {
                Err(FilmError("cancel Carreau trial"))
            } else {
                Ok(())
            }
        });
        assert_eq!(result.unwrap_err(), FilmError("cancel Carreau trial"));
        assert_eq!(f, [77.0; 3]);
        assert_eq!(p, [88.0]);
        assert_eq!(ch, [sentinel; 2]);
    }
    let limited = LaminarCarreauSlit::new(*film.passages()[0].law(), 300.0, W, H, L, 1.0).unwrap();
    let bounded = film.with_passages(vec![limited; 2]).unwrap();
    let mut f = [77.0; 3];
    let mut p = [88.0];
    let mut ch = [sentinel; 2];
    assert!(
        bounded
            .evaluate_into(&q, &v, &mut f, &mut p, &mut ch)
            .is_err()
    );
    assert_eq!(f, [77.0; 3]);
    assert_eq!(p, [88.0]);
    assert_eq!(ch, [sentinel; 2]);
    let (f, p, ch, r) = evaluate(&film, &q, &v);
    assert_eq!(f, expected_f);
    assert_eq!(p, expected_p);
    assert_eq!(ch, expected_ch);
    assert_eq!(
        r.film.dissipated_power_w,
        expected_report.film.dissipated_power_w
    );
    let geometry = base().with_barus_viscosity(1e-6).unwrap();
    assert!(
        CarreauFilm::new(
            &geometry,
            film.passages().to_vec(),
            &[Some(0.0); 2],
            walls(),
            controls()
        )
        .is_err()
    );
    assert!(
        CarreauFilm::new(
            &base(),
            vec![passage(0.4, true, 300.0, W * 2.0); 2],
            &[Some(0.0); 2],
            walls(),
            controls()
        )
        .is_err()
    );
}

#[test]
fn g4_closed_fixed_passages_keep_zero_flow_and_moving_closure_refuses() {
    let mut geometry = base();
    let mut channels = geometry.channels.to_vec();
    channels.push(FilmChannel {
        from: 0,
        to: None,
        width_m: W,
        length_m: L,
        gap: GapPort {
            reference_m: 0.0,
            closure: vec![0.0; 3],
        },
    });
    geometry =
        ResistiveFilm::new(geometry.cells.to_vec(), channels, 3, 0.02, geometry.limits).unwrap();
    let mut ws = walls();
    ws.push(WallVelocityPorts {
        lower: vec![0.0; 3],
        upper: vec![0.0; 3],
    });
    let film = CarreauFilm::new(
        &geometry,
        vec![passage(0.4, true, 300.0, W); 3],
        &[Some(0.0); 3],
        ws.clone(),
        controls(),
    )
    .unwrap();
    let (_, _, ch, _) = evaluate(&film, &[0.0; 3], &[1e-6, 0.01, 0.02]);
    assert_eq!(ch[2], CarreauFilmChannelReport::default());
    ws[2].upper[2] = 1.0;
    let moving = CarreauFilm::new(
        &geometry,
        film.passages().to_vec(),
        &[Some(0.0); 3],
        ws,
        controls(),
    )
    .unwrap();
    let mut f = [77.0; 3];
    let mut p = [88.0];
    let mut c = [CarreauFilmChannelReport::default(); 3];
    assert!(
        moving
            .evaluate_into(&[0.0; 3], &[1e-6, 0.0, 0.0], &mut f, &mut p, &mut c)
            .is_err()
    );
    assert_eq!(f, [77.0; 3]);
    assert_eq!(p, [88.0]);
}
