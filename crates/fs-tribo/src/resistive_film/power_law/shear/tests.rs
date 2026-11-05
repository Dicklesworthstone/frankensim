use super::super::{PowerLawControls, PowerLawFilm};
use super::*;
use crate::resistive_film::shear::{
    ShearChannelReport, ShearedFilm, pressure::PressureDrivenShearedFilm,
};
use crate::resistive_film::{
    FilmCell, FilmChannel, FilmChannelReport, FilmLimits, GapPort, ResistiveFilm,
};

const H: f64 = 2e-4;
const W: f64 = 0.01;
const L: f64 = 0.02;
const AREA: f64 = 1e-4;

fn controls() -> PowerLawControls {
    PowerLawControls {
        maximum_iterations: 80,
        maximum_line_search_steps: 48,
        absolute_flow_tolerance_m3_s: 1e-23,
        relative_flow_tolerance: 3e-12,
        absolute_power_tolerance_w: 1e-21,
        relative_power_tolerance: 5e-12,
    }
}

fn law(n: f64) -> RegularizedPowerLaw {
    RegularizedPowerLaw {
        consistency_pa_s_n: 0.1 / 5.0_f64.powf(n),
        flow_index: n,
        minimum_shear_rate_per_s: 5.0,
    }
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

fn binding(
    base: &ResistiveFilm,
    laws: Vec<RegularizedPowerLaw>,
    supplies: &[Option<f64>],
) -> ShearedPowerLawFilm {
    ShearedPowerLawFilm::new(
        PressureDrivenPowerLawFilm::new(
            PowerLawFilm::new(base, laws, controls()).unwrap(),
            supplies,
        )
        .unwrap(),
        walls(),
    )
    .unwrap()
}

#[track_caller]
fn near(a: f64, b: f64, abs: f64, rel: f64) {
    assert!(
        (a - b).abs() <= abs + rel * a.abs().max(b.abs()),
        "{a:.16e} != {b:.16e}"
    );
}

fn run(
    film: &ShearedPowerLawFilm,
    q: &[f64],
    v: &[f64],
) -> (
    [f64; 3],
    [f64; 1],
    [ShearedPowerLawChannelReport; 2],
    PressureDrivenPowerLawReport,
) {
    let (mut f, mut p, mut c) = (
        [0.0; 3],
        [0.0],
        [ShearedPowerLawChannelReport::default(); 2],
    );
    let report = film.evaluate_into(q, v, &mut f, &mut p, &mut c).unwrap();
    (f, p, c, report)
}

// Independent reconstruction of steady momentum, the wall constraint and
// dissipation. It numerically integrates g(tau) and tau*g(tau), not the
// production antiderivatives, stress moments or pressure-flow formulas.
fn reference(
    material: RegularizedPowerLaw,
    h: f64,
    dp: f64,
    lower: f64,
    upper: f64,
) -> (f64, f64, f64, f64) {
    let mu = material.consistency_pa_s_n
        * material
            .minimum_shear_rate_per_s
            .powf(material.flow_index - 1.0);
    let core = mu * material.minimum_shear_rate_per_s;
    let inverse = |stress: f64| {
        if stress.abs() <= core {
            stress / mu
        } else {
            (stress.abs() / material.consistency_pa_s_n)
                .powf(1.0 / material.flow_index)
                .copysign(stress)
        }
    };
    let s = dp * h / (2.0 * L);
    let integral = |c: f64, f: &dyn Fn(f64, f64) -> f64| {
        let mut cuts = vec![-1.0, 1.0];
        if s != 0.0 {
            for stress in [-core, core] {
                let x = (c - stress) / s;
                if x > -1.0 && x < 1.0 {
                    cuts.push(x);
                }
            }
        }
        cuts.sort_by(f64::total_cmp);
        let mut result = 0.0;
        for pair in cuts.windows(2) {
            let count = 2048;
            let step = (pair[1] - pair[0]) / f64::from(count);
            let at = |x: f64| f(x, c - s * x);
            let mut sum = at(pair[0]) + at(pair[1]);
            for i in 1..count {
                sum += (if i % 2 == 0 { 2.0 } else { 4.0 }) * at(pair[0] + f64::from(i) * step);
            }
            result += sum * step / 3.0;
        }
        0.5 * result
    };
    let rate = (upper - lower) / h;
    let pure = if rate.abs() <= material.minimum_shear_rate_per_s {
        mu * rate.abs()
    } else {
        material.consistency_pa_s_n * rate.abs().powf(material.flow_index)
    };
    let (mut a, mut b) = (-pure - s.abs(), pure + s.abs());
    for _ in 0..80 {
        let c = 0.5 * (a + b);
        if integral(c, &|_, stress| inverse(stress)) < rate {
            a = c;
        } else {
            b = c;
        }
    }
    let c = 0.5 * (a + b);
    let flux = W * h * 0.5 * (lower + upper)
        - 0.5 * W * h * h * integral(c, &|x, stress| x * inverse(stress));
    let heat = W * L * h * integral(c, &|_, stress| stress * inverse(stress));
    (flux, -(c + s), c - s, heat)
}

#[test]
fn g3_newtonian_moving_limit_recovers_native_supplied_film_and_tangent() {
    let base = base();
    let materials = vec![law(1.0); 2];
    let supplies = [Some(230.0), Some(-70.0)];
    let film = binding(&base, materials.clone(), &supplies);
    let native = PressureDrivenShearedFilm::new(
        ShearedFilm::new(
            base.with_channel_viscosities(&[materials[0].consistency_pa_s_n; 2])
                .unwrap(),
            walls(),
        )
        .unwrap(),
        &supplies,
    )
    .unwrap();
    let q = [1e-5, 0.0, 0.0];
    for v in [[2e-5, 0.01, -0.03], [-2e-5, -0.01, 0.03], [0.0; 3]] {
        let (f, p, c, report) = run(&film, &q, &v);
        let (mut ef, mut ep, mut ec) = ([0.0; 3], [0.0], [ShearChannelReport::default(); 2]);
        let expected = native
            .evaluate_into(&q, &v, &mut ef, &mut ep, &mut ec)
            .unwrap();
        for (a, b) in f.into_iter().zip(ef) {
            near(a, b, 1e-15, 2e-13);
        }
        assert_eq!(p, ep);
        assert_eq!(report.power_law.iterations, 0);
        near(
            report.mechanical_power_w,
            expected.mechanical_power_w,
            1e-20,
            2e-13,
        );
        near(
            report.reservoir_power_w,
            expected.reservoir_power_w,
            1e-20,
            2e-13,
        );
        for (a, b) in c.into_iter().zip(ec) {
            near(a.volume_flow_m3_s, b.volume_flow_m3_s, 1e-23, 2e-13);
            near(a.dissipated_power_w(), b.dissipated_power_w(), 1e-20, 2e-13);
        }
        let (dq, dv, db) = ([2e-5, 0.0, 0.0], [1e-5, -0.002, 0.003], [12.0, -9.0]);
        let (mut action, mut expected) = ([0.0; 3], [0.0; 3]);
        film.tangent_into(&q, &v, &dq, &dv, &db, &mut action)
            .unwrap();
        native
            .tangent_into(&q, &v, &dq, &dv, &db, &mut expected)
            .unwrap();
        for (a, b) in action.into_iter().zip(expected) {
            near(a, b, 1e-15, 3e-12);
        }
    }
}

#[test]
fn g1_pure_couette_has_exact_power_law_traction_and_differential_pressure_viscosity() {
    for n in [0.4, 0.7, 1.6, 3.0] {
        let material = law(n);
        let film = binding(&base(), vec![material; 2], &[Some(0.0); 2]);
        for relative in [0.0002, 0.02, -0.02] {
            let v = [0.0, 0.017 - 0.5 * relative, 0.017 + 0.5 * relative];
            let (f, p, c, report) = run(&film, &[0.0; 3], &v);
            let rate = relative / H;
            let viscosity = material.consistency_pa_s_n * rate.abs().max(5.0).powf(n - 1.0);
            let stress = viscosity * rate;
            assert_eq!(p, [0.0]);
            near(f[1], -2.0 * W * L * stress, 1e-16, 2e-13);
            near(f[2], 2.0 * W * L * stress, 1e-16, 2e-13);
            near(c[1].volume_flow_m3_s, W * H * 0.017, 1e-22, 2e-13);
            near(
                report.power_law.film.dissipated_power_w,
                2.0 * W * L * stress * relative,
                1e-20,
                2e-13,
            );
            let original = film.pressure_film().power_law_film();
            let local = channel::evaluate(
                original,
                1,
                H,
                original.film.factors[1] * H * H * H,
                0.0,
                Motion {
                    mean: 0.017,
                    relative,
                    entrainment: W * H * 0.017,
                },
                &mut || Ok(()),
            )
            .unwrap();
            let differential_mu = if rate.abs() <= 5.0 {
                viscosity
            } else {
                n * viscosity
            };
            near(
                local.differential,
                W * H.powi(3) / (12.0 * L * differential_mu),
                1e-25,
                2e-13,
            );
        }
    }
}

#[test]
fn g1_combined_profiles_match_independent_momentum_and_heat_with_active_wall_work() {
    let mut motor = false;
    let mut pump = false;
    for n in [0.4, 0.7, 1.6, 3.0] {
        let material = law(n);
        let relative = H * 5.0 * 10.0_f64.powf(1.0 / n);
        let film = binding(&base(), vec![material; 2], &[Some(600.0), Some(0.0)]);
        let mut previous_heat = None;
        for mean in [-0.3, 0.0, 0.3] {
            let lower = mean - 0.5 * relative;
            let upper = mean + 0.5 * relative;
            let (f, p, c, report) = run(&film, &[0.0; 3], &[0.0, lower, upper]);
            let expected = reference(material, H, 300.0, lower, upper);
            near(p[0], 300.0, 2e-8, 2e-10);
            near(f[0], AREA * 300.0, 1e-12, 2e-10);
            near(c[1].volume_flow_m3_s, expected.0, 2e-20, 2e-9);
            near(c[1].lower_wall_traction_pa, expected.1, 1e-10, 2e-9);
            near(c[1].upper_wall_traction_pa, expected.2, 1e-10, 2e-9);
            near(c[1].dissipated_power_w(), expected.3, 1e-16, 2e-9);
            near(c[0].volume_flow_m3_s, -c[1].volume_flow_m3_s, 1e-22, 5e-12);
            near(
                report.mechanical_power_w + report.reservoir_power_w,
                report.power_law.film.dissipated_power_w,
                1e-19 + 5e-12 * (report.mechanical_power_w.abs() + report.reservoir_power_w.abs()),
                0.0,
            );
            near(
                c[0].wall_power_w + c[1].wall_power_w,
                report.mechanical_power_w,
                1e-19,
                3e-12,
            );
            if let Some(heat) = previous_heat {
                near(report.power_law.film.dissipated_power_w, heat, 1e-18, 5e-12);
            }
            previous_heat = Some(report.power_law.film.dissipated_power_w);
            motor |= report.mechanical_power_w < 0.0;
            pump |= report.reservoir_power_w < 0.0;
            if mean == 0.0 {
                let reversed = film
                    .with_reservoir_pressures(&[Some(0.0), Some(600.0)])
                    .unwrap();
                let (_, _, reversed, rr) = run(&reversed, &[0.0; 3], &[0.0, -lower, -upper]);
                near(
                    reversed[1].volume_flow_m3_s,
                    -c[1].volume_flow_m3_s,
                    1e-20,
                    3e-10,
                );
                near(
                    reversed[1].lower_wall_traction_pa,
                    -c[1].lower_wall_traction_pa,
                    1e-10,
                    3e-10,
                );
                near(
                    reversed[1].upper_wall_traction_pa,
                    -c[1].upper_wall_traction_pa,
                    1e-10,
                    3e-10,
                );
                near(
                    rr.power_law.film.dissipated_power_w,
                    report.power_law.film.dissipated_power_w,
                    1e-18,
                    3e-10,
                );
            }
        }
    }
    assert!(
        motor && pump,
        "both signs of actual external power must be exercised"
    );
}

#[test]
fn g3_fixed_wall_binding_keeps_original_supplied_arithmetic() {
    let original = PressureDrivenPowerLawFilm::new(
        PowerLawFilm::new(&base(), vec![law(0.5); 2], controls()).unwrap(),
        &[Some(230.0), Some(-70.0)],
    )
    .unwrap();
    let fixed = ShearedPowerLawFilm::new(
        original.clone(),
        vec![
            WallVelocityPorts {
                lower: vec![0.0; 3],
                upper: vec![0.0; 3],
            },
            WallVelocityPorts {
                lower: vec![0.0; 3],
                upper: vec![0.0; 3],
            },
        ],
    )
    .unwrap();
    let (q, v) = ([1e-5, 0.0, 0.0], [2e-5, 0.01, -0.03]);
    let (f, p, c, report) = run(&fixed, &q, &v);
    let (mut ef, mut ep, mut ec) = ([0.0; 3], [0.0], [FilmChannelReport::default(); 2]);
    let expected = original
        .evaluate_into(&q, &v, &mut ef, &mut ep, &mut ec)
        .unwrap();
    assert_eq!(f, ef);
    assert_eq!(p, ep);
    assert_eq!(report.mechanical_power_w, expected.mechanical_power_w);
    assert_eq!(report.reservoir_power_w, expected.reservoir_power_w);
    assert_eq!(
        report.power_law.film.dissipated_power_w,
        expected.power_law.film.dissipated_power_w
    );
    for (a, b) in c.into_iter().zip(ec) {
        assert_eq!(a.volume_flow_m3_s, b.volume_flow_m3_s);
        assert_eq!(a.dissipated_power_w(), b.dissipated_power_w);
        assert_eq!(a.wall_power_w, 0.0);
    }
}

#[test]
fn g1_combined_geometry_velocity_and_supply_tangent_is_consistent_and_reciprocal() {
    let original = base();
    let mut edges = original.channels.to_vec();
    edges[0].gap.closure = vec![0.4, 0.2, -0.1];
    edges[1].gap.closure = vec![1.0, 0.1, 0.1];
    let base =
        ResistiveFilm::new(original.cells.to_vec(), edges, 3, 0.02, original.limits).unwrap();
    let film = binding(&base, vec![law(0.5), law(1.7)], &[Some(200.0), Some(-50.0)]);
    let q = [1e-5, -2e-5, 1e-5];
    let v = [2e-5, 0.003, -0.006];
    let dq = [2e-5, -1e-5, 0.5e-5];
    let dv = [1e-5, 0.0008, -0.0004];
    let db = [12.0, -9.0];
    let mut action = [0.0; 3];
    film.tangent_into(&q, &v, &dq, &dv, &db, &mut action)
        .unwrap();
    let epsilon = 2e-5;
    let perturb = |sign: f64| {
        let supplies = [
            Some(200.0 + sign * epsilon * db[0]),
            Some(-50.0 + sign * epsilon * db[1]),
        ];
        let shifted = film.with_reservoir_pressures(&supplies).unwrap();
        let qq = std::array::from_fn::<_, 3, _>(|j| q[j] + sign * epsilon * dq[j]);
        let vv = std::array::from_fn::<_, 3, _>(|j| v[j] + sign * epsilon * dv[j]);
        run(&shifted, &qq, &vv).0
    };
    let plus = perturb(1.0);
    let minus = perturb(-1.0);
    for j in 0..3 {
        near(
            action[j],
            (plus[j] - minus[j]) / (2.0 * epsilon),
            2e-10,
            3e-7,
        );
    }
    let mut jacobian = [[0.0; 3]; 3];
    for (j, column) in jacobian.iter_mut().enumerate() {
        let mut direction = [0.0; 3];
        direction[j] = 1.0;
        film.tangent_into(&q, &v, &[0.0; 3], &direction, &[0.0; 2], column)
            .unwrap();
    }
    for i in 0..3 {
        for j in 0..3 {
            near(jacobian[i][j], jacobian[j][i], 1e-10, 5e-12);
        }
    }
    let quadratic: f64 = (0..3)
        .map(|i| (0..3).map(|j| dv[i] * jacobian[i][j] * dv[j]).sum::<f64>())
        .sum();
    assert!(
        quadratic > 0.0,
        "incremental mechanical resistance must be positive"
    );
    let (_, p, c, report) = run(&film, &q, &v);
    near(
        c[0].volume_flow_m3_s + c[1].volume_flow_m3_s,
        AREA * v[0],
        1e-22,
        5e-12,
    );
    let normal = AREA * p[0] * v[0];
    near(
        normal + c[0].wall_power_w + c[1].wall_power_w,
        report.mechanical_power_w,
        1e-20,
        5e-12,
    );
}

#[test]
fn g1_extreme_pressure_to_sliding_ratios_preserve_channel_derivatives_and_symmetries() {
    for n in [0.35, 0.8, 1.5, 3.0] {
        let original = PowerLawFilm::new(&base(), vec![law(n); 2], controls()).unwrap();
        for (dp, relative) in [
            (300.0, 1e-12),
            (1e-12, 0.01),
            (300.0, 0.01),
            (-300.0, -0.01),
        ] {
            let motion = Motion {
                mean: 0.0,
                relative,
                entrainment: 0.0,
            };
            let response = channel::evaluate(
                &original,
                0,
                H,
                original.film.factors[0] * H * H * H,
                dp,
                motion,
                &mut || Ok(()),
            )
            .unwrap();
            let reversed = channel::evaluate(
                &original,
                0,
                H,
                original.film.factors[0] * H * H * H,
                -dp,
                Motion {
                    relative: -relative,
                    ..motion
                },
                &mut || Ok(()),
            )
            .unwrap();
            near(response.flow, -reversed.flow, 1e-25, 2e-13);
            near(
                response.wall.relative_effort,
                -reversed.wall.relative_effort,
                1e-20,
                2e-13,
            );
            assert!(response.differential > 0.0 && response.wall.relative_effort_derivative > 0.0);
            // A small pressure gradient on a strongly sliding film uses the
            // differential (n*mu) viscosity, even when dp is below stress ULPs.
            if dp == 1e-12 {
                let material = law(n);
                let mu = material.consistency_pa_s_n * (relative / H).powf(n - 1.0);
                near(
                    response.flow / dp,
                    W * H.powi(3) / (12.0 * L * n * mu),
                    1e-24,
                    3e-12,
                );
            }
        }
    }
    // Relative rows must be compiled before evaluating velocities: separate
    // upper/lower dot products would erase this finite slide in common motion.
    let fixed = PressureDrivenPowerLawFilm::new(
        PowerLawFilm::new(&base(), vec![law(0.5); 2], controls()).unwrap(),
        &[Some(0.0); 2],
    )
    .unwrap();
    let film = ShearedPowerLawFilm::new(
        fixed,
        vec![
            WallVelocityPorts {
                lower: vec![0.0, -1.0, 0.0],
                upper: vec![0.0, -1.0, -1.0],
            },
            WallVelocityPorts {
                lower: vec![0.0, 1.0, 0.0],
                upper: vec![0.0, 1.0, 1.0],
            },
        ],
    )
    .unwrap();
    let (force, pressure, _, report) = run(&film, &[0.0; 3], &[0.0, 1e16, 0.01]);
    let stress = law(0.5).consistency_pa_s_n * (0.01 / H).sqrt();
    assert_eq!(pressure, [0.0]);
    assert_eq!(force[1], 0.0);
    near(force[2], 2.0 * W * L * stress, 1e-18, 3e-13);
    near(
        report.power_law.film.dissipated_power_w,
        2.0 * W * L * stress * 0.01,
        1e-20,
        3e-13,
    );
    // Supplied pressure and a nonzero swept source must survive in both the
    // nonlinear owner and its native n=1 fast path. Total transported Q is too
    // large to represent its tiny pressure excess, so compare that excess too.
    for n in [0.5, 1.0] {
        let supplied = film
            .with_laws(vec![law(n); 2])
            .unwrap()
            .with_reservoir_pressures(&[Some(600.0), Some(0.0)])
            .unwrap();
        let still = [2e-5, 0.0, 0.01];
        let moving = [2e-5, 1e16, 0.01];
        let (expected_force, expected_pressure, expected_channels, expected_report) =
            run(&supplied, &[0.0; 3], &still);
        let (force, pressure, channels, report) = run(&supplied, &[0.0; 3], &moving);
        near(pressure[0], expected_pressure[0], 1e-9, 5e-12);
        assert!(
            pressure[0] > 300.0,
            "the small normal swept source must raise film pressure"
        );
        for (a, b) in force.into_iter().zip(expected_force) {
            near(a, b, 1e-13, 5e-12);
        }
        for (a, b) in channels.into_iter().zip(expected_channels) {
            near(a.pressure_flow_m3_s, b.pressure_flow_m3_s, 1e-22, 5e-12);
            near(a.dissipated_power_w(), b.dissipated_power_w(), 1e-19, 5e-12);
        }
        near(
            report.power_law.film.dissipated_power_w,
            expected_report.power_law.film.dissipated_power_w,
            1e-19,
            5e-12,
        );
        let (dq, dv, db) = ([1e-5, 0.0, 0.0], [1e-5, 0.002, -0.003], [11.0, -7.0]);
        let (mut tangent, mut expected) = ([0.0; 3], [0.0; 3]);
        supplied
            .tangent_into(&[0.0; 3], &moving, &dq, &dv, &db, &mut tangent)
            .unwrap();
        supplied
            .tangent_into(&[0.0; 3], &still, &dq, &dv, &db, &mut expected)
            .unwrap();
        for (a, b) in tangent.into_iter().zip(expected) {
            near(a, b, 1e-12, 2e-10);
        }
    }
}

#[test]
fn g4_wall_roots_domains_cancellation_and_core_corner_refuse_atomically() {
    let base = base();
    let film = binding(&base, vec![law(0.5); 2], &[Some(600.0), Some(0.0)]);
    let (mut f, mut p, mut c) = (
        [71.0; 3],
        [72.0],
        [ShearedPowerLawChannelReport {
            volume_flow_m3_s: 73.0,
            ..ShearedPowerLawChannelReport::default()
        }; 2],
    );
    let before = (f, p, c);
    let mut calls = 0;
    let error = film
        .evaluate_with_check_into(
            &[0.0; 3],
            &[0.0, 0.001, 0.02],
            &mut f,
            &mut p,
            &mut c,
            || {
                calls += 1;
                if calls == 6 {
                    Err(FilmError("caller wall cancellation"))
                } else {
                    Ok(())
                }
            },
        )
        .unwrap_err();
    assert_eq!(error, FilmError("caller wall cancellation"));
    assert_eq!((f, p, c), before);
    let mut short = controls();
    short.maximum_iterations = 1;
    let limited = ShearedPowerLawFilm::new(
        PressureDrivenPowerLawFilm::new(
            PowerLawFilm::new(&base, vec![law(0.5); 2], short).unwrap(),
            &[Some(600.0), Some(0.0)],
        )
        .unwrap(),
        walls(),
    )
    .unwrap();
    assert!(
        limited
            .evaluate_into(&[0.0; 3], &[0.0, 0.001, 0.02], &mut f, &mut p, &mut c)
            .is_err()
    );
    assert_eq!((f, p, c), before);
    let mut closed_edges = base.channels.to_vec();
    closed_edges[0].gap = GapPort {
        reference_m: 0.0,
        closure: vec![0.0; 3],
    };
    let closed_base =
        ResistiveFilm::new(base.cells.to_vec(), closed_edges, 3, 0.02, base.limits).unwrap();
    let closed = binding(&closed_base, vec![law(0.5); 2], &[Some(0.0); 2]);
    assert!(
        closed
            .evaluate_into(&[0.0; 3], &[0.0; 3], &mut f, &mut p, &mut c)
            .is_err()
    );
    assert_eq!((f, p, c), before);
    let ambient = film.with_reservoir_pressures(&[Some(0.0); 2]).unwrap();
    let corner = [0.0, 0.0, H * 5.0];
    run(&ambient, &[0.0; 3], &corner);
    let mut out = [74.0; 3];
    assert_eq!(
        ambient.tangent_into(
            &[0.0; 3],
            &corner,
            &[0.0; 3],
            &[0.0, 0.0, 1.0],
            &[0.0; 2],
            &mut out
        ),
        Err(FilmError(
            "power-law wall tangent is undefined at the pure-Couette core corner"
        ))
    );
    assert_eq!(out, [74.0; 3]);
    let mut invalid_walls = walls();
    invalid_walls[0].lower[1] = f64::NAN;
    assert!(ShearedPowerLawFilm::new(film.pressure_film().clone(), invalid_walls).is_err());
}
