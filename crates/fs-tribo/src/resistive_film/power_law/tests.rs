use super::super::{FilmCell, FilmChannel, FilmLimits, GapPort};
use super::*;

const H: f64 = 2e-4;
const W: f64 = 0.01;
const L: f64 = 0.02;
const AREA: f64 = 1e-4;

fn controls() -> PowerLawControls {
    PowerLawControls {
        maximum_iterations: 80,
        maximum_line_search_steps: 48,
        absolute_flow_tolerance_m3_s: 1e-23,
        relative_flow_tolerance: 2e-12,
        absolute_power_tolerance_w: 1e-21,
        relative_power_tolerance: 3e-12,
    }
}

fn law(k: f64, n: f64, core: f64) -> RegularizedPowerLaw {
    RegularizedPowerLaw {
        consistency_pa_s_n: k,
        flow_index: n,
        minimum_shear_rate_per_s: core,
    }
}

fn limits() -> FilmLimits {
    FilmLimits {
        minimum_cell_gap_m: 1e-6,
        maximum_gap_m: 1e-3,
        maximum_pressure_pa: 1e7,
    }
}

fn one_channel() -> ResistiveFilm {
    let gap = GapPort {
        reference_m: H,
        closure: vec![1.0],
    };
    ResistiveFilm::new(
        vec![FilmCell {
            area_m2: AREA,
            gap: gap.clone(),
        }],
        vec![FilmChannel {
            from: 0,
            to: None,
            width_m: W,
            length_m: L,
            gap,
        }],
        1,
        0.01,
        limits(),
    )
    .unwrap()
}

fn network() -> ResistiveFilm {
    ResistiveFilm::new(
        vec![
            FilmCell {
                area_m2: 2.0 * AREA,
                gap: GapPort {
                    reference_m: H,
                    closure: vec![1.0, 0.0, 0.0],
                },
            },
            FilmCell {
                area_m2: 1.5 * AREA,
                gap: GapPort {
                    reference_m: H,
                    closure: vec![0.0, 1.0, 0.0],
                },
            },
        ],
        vec![
            FilmChannel {
                from: 0,
                to: Some(1),
                width_m: W,
                length_m: L,
                gap: GapPort {
                    reference_m: H,
                    closure: vec![0.3, 0.2, 1.0],
                },
            },
            FilmChannel {
                from: 1,
                to: None,
                width_m: 2.0 * W,
                length_m: 1.2 * L,
                gap: GapPort {
                    reference_m: H,
                    closure: vec![0.0, 0.4, 0.2],
                },
            },
            FilmChannel {
                from: 0,
                to: None,
                width_m: 0.7 * W,
                length_m: 2.0 * L,
                gap: GapPort {
                    reference_m: H,
                    closure: vec![0.6, 0.0, -0.3],
                },
            },
            FilmChannel {
                from: 0,
                to: None,
                width_m: W,
                length_m: L,
                gap: GapPort {
                    reference_m: 0.0,
                    closure: vec![0.0; 3],
                },
            },
        ],
        3,
        0.01,
        limits(),
    )
    .unwrap()
}

fn near(a: f64, b: f64, absolute: f64, relative: f64) {
    let allowance = absolute + relative * a.abs().max(b.abs());
    assert!(
        (a - b).abs() <= allowance,
        "{a} != {b}, allowance {allowance}"
    );
}

fn simpson(a: f64, b: f64, f: impl Fn(f64) -> f64) -> f64 {
    if a == b {
        return 0.0;
    }
    let count = 4096;
    let step = (b - a) / f64::from(count);
    let mut sum = f(a) + f(b);
    for i in 1..count {
        sum += f(a + f64::from(i) * step) * if i % 2 == 0 { 2.0 } else { 4.0 };
    }
    sum * step / 3.0
}

/// Independent no-slip momentum quadrature. This reconstructs the local
/// shear rate from physical stress and integrates y*gamma(y), without using
/// the implementation's compiled flux, blended fractions or derivatives.
fn reference_flow(law: RegularizedPowerLaw, width: f64, length: f64, h: f64, dp: f64) -> f64 {
    if h <= 0.0 || dp == 0.0 {
        return 0.0;
    }
    let a = 0.5 * h;
    let gradient = dp.abs() / length;
    let mu0 = law.consistency_pa_s_n * law.minimum_shear_rate_per_s.powf(law.flow_index - 1.0);
    let tau0 = mu0 * law.minimum_shear_rate_per_s;
    let interface = a.min(tau0 / gradient);
    let core = simpson(0.0, interface, |y| y * gradient * y / mu0);
    let outer = simpson(interface, a, |y| {
        y * (gradient * y / law.consistency_pa_s_n).powf(1.0 / law.flow_index)
    });
    (2.0 * width * (core + outer)).copysign(dp)
}

fn mixed_laws() -> Vec<RegularizedPowerLaw> {
    vec![
        law(0.02, 0.5, 1.0),
        law(0.04, 1.7, 1.0),
        law(0.03, 0.8, 1.0),
        law(0.01, 1.2, 1.0),
    ]
}

fn manufactured_rates(
    film: &ResistiveFilm,
    laws: &[RegularizedPowerLaw],
    q: &[f64],
    p: [f64; 2],
) -> ([f64; 3], [f64; 4]) {
    let mut swept = [0.0; 2];
    let mut flows = [0.0; 4];
    for (i, edge) in film.channels.iter().enumerate() {
        let h = edge.gap.reference_m - dot(&edge.gap.closure, q);
        let dp = p[edge.from] - edge.to.map_or(0.0, |j| p[j]);
        let flow = reference_flow(laws[i], edge.width_m, edge.length_m, h, dp);
        flows[i] = flow;
        swept[edge.from] += flow;
        if let Some(j) = edge.to {
            swept[j] -= flow;
        }
    }
    (
        [swept[0] / (2.0 * AREA), swept[1] / (1.5 * AREA), 0.0],
        flows,
    )
}

#[test]
fn g3_newtonian_laws_recover_original_force_channel_work_and_tangent_bits() {
    let base = network();
    let viscosities = [0.003, 0.01, 0.02, 0.5];
    let laws: Vec<_> = viscosities
        .into_iter()
        .enumerate()
        .map(|(i, k)| law(k, 1.0, 10.0_f64.powi(i as i32 - 1)))
        .collect();
    let nonlinear = PowerLawFilm::new(&base, laws.clone(), controls()).unwrap();
    let native = base.with_channel_viscosities(&viscosities).unwrap();
    assert_eq!(nonlinear.laws(), laws);
    assert_eq!(
        nonlinear.pressure_operator().channel_viscosities_pa_s(),
        viscosities
    );
    let q = [1e-5, -2e-5, 5e-6];
    for v in [[2e-4, -1e-4, 0.0], [-2e-4, 1e-4, 0.0], [0.0; 3]] {
        let (mut force, mut pressure, mut channels) =
            ([0.0; 3], [0.0; 2], [FilmChannelReport::default(); 4]);
        let (mut expected_force, mut expected_pressure, mut expected_channels) =
            ([0.0; 3], [0.0; 2], [FilmChannelReport::default(); 4]);
        let report = nonlinear
            .evaluate_with_channels_into(&q, &v, &mut force, &mut pressure, &mut channels)
            .unwrap();
        let expected = native
            .evaluate_with_channels_into(
                &q,
                &v,
                &mut expected_force,
                &mut expected_pressure,
                &mut expected_channels,
            )
            .unwrap();
        assert_eq!(force, expected_force);
        assert_eq!(pressure, expected_pressure);
        assert_eq!(channels, expected_channels);
        assert_eq!(report.film.dissipated_power_w, expected.dissipated_power_w);
        assert_eq!(report.film.minimum_gap_m, expected.minimum_gap_m);
        assert_eq!(report.film.peak_pressure_pa, expected.peak_pressure_pa);
        assert_eq!(report.iterations, 0);
        let dq = [0.2, -0.3, 0.5];
        let dv = [-0.01, 0.03, 0.02];
        let (mut tangent, mut expected_tangent) = ([0.0; 3], [0.0; 3]);
        nonlinear
            .tangent_into(&q, &v, &dq, &dv, &mut tangent)
            .unwrap();
        native
            .tangent_into(&q, &v, &dq, &dv, &mut expected_tangent)
            .unwrap();
        assert_eq!(tangent, expected_tangent);
    }
}

#[test]
fn g1_thinning_and_thickening_pressures_match_independent_no_slip_profile_integrals() {
    let film = one_channel();
    for n in [0.5, 1.7, 3.0] {
        let material = law(0.04, n, 3.0);
        let transition_pressure =
            material.consistency_pa_s_n * material.minimum_shear_rate_per_s.powf(n) * L / (0.5 * H);
        let model = PowerLawFilm::new(&film, vec![material], controls()).unwrap();
        for multiple in [0.05, 0.999, 1.0, 1.001, 4.0, -4.0] {
            let expected_pressure = transition_pressure * multiple;
            let flux = reference_flow(material, W, L, H, expected_pressure);
            let (mut force, mut pressure, mut channels) =
                ([0.0], [0.0], [FilmChannelReport::default()]);
            let report = model
                .evaluate_with_channels_into(
                    &[0.0],
                    &[flux / AREA],
                    &mut force,
                    &mut pressure,
                    &mut channels,
                )
                .unwrap();
            near(pressure[0], expected_pressure, 1e-9, 3e-10);
            near(force[0], AREA * expected_pressure, 1e-13, 3e-10);
            near(channels[0].volume_flow_m3_s, flux, 2e-23, 4e-12);
            near(
                report.film.dissipated_power_w,
                flux * expected_pressure,
                1e-20,
                3e-10,
            );
            assert!(report.film.dissipated_power_w > 0.0);
            assert!(report.iterations > 0);
            assert!(
                report.maximum_flow_defect_m3_s
                    <= controls().absolute_flow_tolerance_m3_s
                        + controls().relative_flow_tolerance * flux.abs()
            );
        }
    }
}

#[test]
fn g1_mixed_series_and_parallel_passages_recover_manufactured_pressures_and_local_continuity() {
    let base = network();
    let laws = mixed_laws();
    let model = PowerLawFilm::new(&base, laws.clone(), controls()).unwrap();
    let q = [1e-5, -2e-5, 5e-6];
    for target in [[400.0, -100.0], [-400.0, 100.0], [40.0, 18.0]] {
        let (rates, flows) = manufactured_rates(&base, &laws, &q, target);
        let (mut force, mut pressure, mut channels) =
            ([0.0; 3], [0.0; 2], [FilmChannelReport::default(); 4]);
        let report = model
            .evaluate_with_channels_into(&q, &rates, &mut force, &mut pressure, &mut channels)
            .unwrap();
        for (actual, expected) in pressure.into_iter().zip(target) {
            near(actual, expected, 1e-8, 3e-8);
        }
        for (channel, expected) in channels.iter().zip(flows) {
            near(channel.volume_flow_m3_s, expected, 1e-22, 1e-8);
            assert!(channel.dissipated_power_w >= 0.0);
        }
        assert_eq!(channels[3], FilmChannelReport::default());
        let residual = [
            channels[0].volume_flow_m3_s + channels[2].volume_flow_m3_s - 2.0 * AREA * rates[0],
            channels[1].volume_flow_m3_s - channels[0].volume_flow_m3_s - 1.5 * AREA * rates[1],
        ];
        for r in residual {
            assert!(
                r.abs()
                    <= 3e-23
                        + 3e-12
                            * (2.0 * AREA * rates[0])
                                .abs()
                                .max((1.5 * AREA * rates[1]).abs())
            );
        }
        near(
            dot(&force, &rates),
            report.film.dissipated_power_w,
            1e-20,
            4e-12,
        );
        near(
            channels.iter().map(|r| r.dissipated_power_w).sum(),
            report.film.dissipated_power_w,
            0.0,
            1e-15,
        );
    }
}

#[test]
fn g3_full_directional_tangent_is_reciprocal_in_velocity_and_matches_geometry_differences() {
    let base = network();
    let laws = mixed_laws();
    let model = PowerLawFilm::new(&base, laws.clone(), controls()).unwrap();
    let q = [1e-5, -2e-5, 5e-6];
    let (rates, _) = manufactured_rates(&base, &laws, &q, [40.0, -15.0]);
    let a = [0.02, -0.01, 0.005];
    let b = [-0.003, 0.007, -0.002];
    let (mut ta, mut tb) = ([0.0; 3], [0.0; 3]);
    model
        .tangent_into(&q, &rates, &[0.0; 3], &a, &mut ta)
        .unwrap();
    model
        .tangent_into(&q, &rates, &[0.0; 3], &b, &mut tb)
        .unwrap();
    near(dot(&a, &tb), dot(&b, &ta), 1e-12, 3e-13);
    assert!(dot(&a, &ta) > 0.0);
    let dq = [1e-4, -5e-5, -2e-5];
    let dv = [2e-4, -3e-4, 1e-4];
    let mut tangent = [0.0; 3];
    model
        .tangent_into(&q, &rates, &dq, &dv, &mut tangent)
        .unwrap();
    let epsilon = 1e-4;
    let mut probes = [[0.0; 3]; 2];
    for (slot, sign) in [-1.0, 1.0].into_iter().enumerate() {
        let state = std::array::from_fn::<_, 3, _>(|i| q[i] + sign * epsilon * dq[i]);
        let velocity = std::array::from_fn::<_, 3, _>(|i| rates[i] + sign * epsilon * dv[i]);
        model
            .evaluate_into(&state, &velocity, &mut probes[slot], &mut [0.0; 2])
            .unwrap();
    }
    for i in 0..3 {
        let independent = (probes[1][i] - probes[0][i]) / (2.0 * epsilon);
        near(tangent[i], independent, 1e-9, 3e-6);
    }
}

#[test]
fn g3_consistency_scaling_changes_pressure_work_and_force_without_changing_volume_flow() {
    let base = network();
    let laws = mixed_laws();
    let scale = 37.0;
    let scaled_laws = laws
        .iter()
        .map(|law| RegularizedPowerLaw {
            consistency_pa_s_n: scale * law.consistency_pa_s_n,
            ..*law
        })
        .collect();
    let first = PowerLawFilm::new(&base, laws.clone(), controls()).unwrap();
    let second = PowerLawFilm::new(&base, scaled_laws, controls()).unwrap();
    let q = [1e-5, -2e-5, 5e-6];
    let (rates, _) = manufactured_rates(&base, &laws, &q, [40.0, -15.0]);
    let mut results = Vec::new();
    for model in [&first, &second] {
        let (mut force, mut pressure, mut channels) =
            ([0.0; 3], [0.0; 2], [FilmChannelReport::default(); 4]);
        let report = model
            .evaluate_with_channels_into(&q, &rates, &mut force, &mut pressure, &mut channels)
            .unwrap();
        results.push((force, pressure, channels, report.film.dissipated_power_w));
    }
    for i in 0..3 {
        near(results[1].0[i], scale * results[0].0[i], 1e-12, 3e-10);
    }
    for i in 0..2 {
        near(results[1].1[i], scale * results[0].1[i], 1e-8, 3e-10);
    }
    for i in 0..4 {
        near(
            results[1].2[i].volume_flow_m3_s,
            results[0].2[i].volume_flow_m3_s,
            1e-22,
            3e-11,
        );
    }
    near(results[1].3, scale * results[0].3, 1e-18, 3e-10);
}

#[test]
fn work_domain_law_and_late_callback_refusals_leave_all_outputs_unchanged() {
    let base = one_channel();
    let valid = law(0.04, 0.5, 3.0);
    for bad in [
        law(0.0, 0.5, 3.0),
        law(0.04, 0.0, 3.0),
        law(0.04, 0.5, 0.0),
        law(f64::NAN, 0.5, 3.0),
        law(0.04, f64::INFINITY, 3.0),
    ] {
        assert!(PowerLawFilm::new(&base, vec![bad], controls()).is_err());
    }
    assert!(PowerLawFilm::new(&base, vec![], controls()).is_err());
    assert!(
        PowerLawFilm::new(
            &base.with_barus_viscosity(1e-4).unwrap(),
            vec![valid],
            controls()
        )
        .is_err()
    );
    for bad in [
        PowerLawControls {
            maximum_iterations: 0,
            ..controls()
        },
        PowerLawControls {
            maximum_iterations: 129,
            ..controls()
        },
        PowerLawControls {
            maximum_line_search_steps: 65,
            ..controls()
        },
        PowerLawControls {
            absolute_flow_tolerance_m3_s: -1.0,
            ..controls()
        },
        PowerLawControls {
            relative_power_tolerance: 1.0,
            ..controls()
        },
    ] {
        assert!(PowerLawFilm::new(&base, vec![valid], bad).is_err());
    }
    let flux = reference_flow(valid, W, L, H, 100.0);
    let model = PowerLawFilm::new(&base, vec![valid], controls()).unwrap();
    let mut checks = 0;
    let (mut expected_force, mut expected_pressure, mut expected_channels) =
        ([0.0], [0.0], [FilmChannelReport::default()]);
    model
        .evaluate_with_channels_and_check_into(
            &[0.0],
            &[flux / AREA],
            &mut expected_force,
            &mut expected_pressure,
            &mut expected_channels,
            || {
                checks += 1;
                Ok(())
            },
        )
        .unwrap();
    assert!(checks > 3);
    let sentinel = FilmChannelReport {
        volume_flow_m3_s: 42.0,
        dissipated_power_w: 43.0,
    };
    let (mut force, mut pressure, mut channels) = ([42.0], [42.0], [sentinel]);
    let mut calls = 0;
    let cancellation = FilmError("test cancellation");
    assert_eq!(
        model
            .evaluate_with_channels_and_check_into(
                &[0.0],
                &[flux / AREA],
                &mut force,
                &mut pressure,
                &mut channels,
                || {
                    calls += 1;
                    if calls == checks {
                        Err(cancellation)
                    } else {
                        Ok(())
                    }
                },
            )
            .unwrap_err(),
        cancellation
    );
    assert_eq!((force, pressure, channels), ([42.0], [42.0], [sentinel]));
    let limited = PowerLawFilm::new(
        &base,
        vec![valid],
        PowerLawControls {
            maximum_iterations: 1,
            ..controls()
        },
    )
    .unwrap();
    assert!(
        limited
            .evaluate_with_channels_into(
                &[0.0],
                &[flux / AREA],
                &mut force,
                &mut pressure,
                &mut channels,
            )
            .is_err()
    );
    let mut small_domain = base.clone();
    small_domain.limits.maximum_pressure_pa = 1.0;
    let small_domain = PowerLawFilm::new(&small_domain, vec![valid], controls()).unwrap();
    assert!(
        small_domain
            .evaluate_with_channels_into(
                &[0.0],
                &[flux / AREA],
                &mut force,
                &mut pressure,
                &mut channels,
            )
            .is_err()
    );
    assert!(
        model
            .evaluate_with_channels_into(
                &[H],
                &[flux / AREA],
                &mut force,
                &mut pressure,
                &mut channels,
            )
            .is_err()
    );
    assert_eq!((force, pressure, channels), ([42.0], [42.0], [sentinel]));
    let mut sealed = base.clone();
    std::sync::Arc::make_mut(&mut sealed.channels)[0]
        .gap
        .reference_m = 0.0;
    let sealed = PowerLawFilm::new(&sealed, vec![valid], controls()).unwrap();
    assert!(
        sealed
            .evaluate_with_channels_into(&[0.0], &[0.0], &mut force, &mut pressure, &mut channels,)
            .is_err()
    );
    assert_eq!((force, pressure, channels), ([42.0], [42.0], [sentinel]));
    let mut tangent = [42.0];
    assert_eq!(
        model.tangent_with_check_into(
            &[0.0],
            &[flux / AREA],
            &[0.0],
            &[1.0],
            &mut tangent,
            || Err(cancellation),
        ),
        Err(cancellation)
    );
    assert_eq!(tangent, [42.0]);
    model
        .evaluate_with_channels_into(
            &[0.0],
            &[flux / AREA],
            &mut force,
            &mut pressure,
            &mut channels,
        )
        .unwrap();
    assert_eq!(
        (force, pressure, channels),
        (expected_force, expected_pressure, expected_channels)
    );
}

#[test]
fn force_and_tangent_share_strict_projected_work_admission() {
    // G0: unit channel conductance isolates the one-ULP difference between
    // p*(A*B*v) and v*(B^T*A*p). The latter is the actual mechanical port work,
    // and force/tangent must share the explicitly requested zero-defect gate.
    let base = ResistiveFilm::new(
        vec![FilmCell {
            area_m2: 0.1,
            gap: GapPort {
                reference_m: 1.0,
                closure: vec![0.3],
            },
        }],
        vec![FilmChannel {
            from: 0,
            to: None,
            width_m: 12.0,
            length_m: 1.0,
            gap: GapPort {
                reference_m: 1.0,
                closure: vec![0.0],
            },
        }],
        1,
        1.0,
        FilmLimits {
            minimum_cell_gap_m: 0.1,
            maximum_gap_m: 2.0,
            maximum_pressure_pa: 1.0,
        },
    )
    .unwrap();
    let q = [0.0];
    let velocity = [0.7];
    let (mut native_force, mut native_pressure) = ([0.0], [0.0]);
    let native_report = base
        .evaluate_into(&q, &velocity, &mut native_force, &mut native_pressure)
        .unwrap();
    assert_ne!(
        native_force[0] * velocity[0],
        native_report.dissipated_power_w
    );
    let mut native_tangent = [0.0];
    base.tangent_into(&q, &velocity, &[0.0], &[1.0], &mut native_tangent)
        .unwrap();
    for n in [1.0, 2.0] {
        // n=2 stays inside its declared Newtonian core for this state. Both
        // native and nonlinear pressure paths therefore encounter this gate.
        let strict = PowerLawControls {
            absolute_flow_tolerance_m3_s: 0.0,
            relative_flow_tolerance: 0.0,
            absolute_power_tolerance_w: 0.0,
            relative_power_tolerance: 0.0,
            ..controls()
        };
        let strict = PowerLawFilm::new(&base, vec![law(1.0, n, 1.0)], strict).unwrap();
        let (mut force, mut pressure, mut tangent) = ([42.0], [43.0], [44.0]);
        let force_error = strict
            .evaluate_into(&q, &velocity, &mut force, &mut pressure)
            .unwrap_err();
        let tangent_error = strict
            .tangent_into(&q, &velocity, &[0.0], &[1.0], &mut tangent)
            .unwrap_err();
        assert_eq!(force_error, tangent_error);
        assert_eq!((force, pressure, tangent), ([42.0], [43.0], [44.0]));
        let admitted = PowerLawFilm::new(
            &base,
            vec![law(1.0, n, 1.0)],
            PowerLawControls {
                absolute_power_tolerance_w: 1e-18,
                ..strict.controls()
            },
        )
        .unwrap();
        admitted
            .evaluate_into(&q, &velocity, &mut force, &mut pressure)
            .unwrap();
        admitted
            .tangent_into(&q, &velocity, &[0.0], &[1.0], &mut tangent)
            .unwrap();
        assert_eq!(
            (force, pressure, tangent),
            (native_force, native_pressure, native_tangent)
        );
    }
}
