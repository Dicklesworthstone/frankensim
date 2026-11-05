use super::*;
use crate::resistive_film::pressure::PressureDrivenFilm;

fn supplied_rates(
    film: &ResistiveFilm,
    laws: &[RegularizedPowerLaw],
    q: &[f64],
    pressure: [f64; 2],
    supplies: &[Option<f64>],
) -> ([f64; 3], [f64; 4]) {
    let mut swept = [0.0; 2];
    let mut flows = [0.0; 4];
    for (k, edge) in film.channels.iter().enumerate() {
        let h = edge.gap.reference_m - dot(&edge.gap.closure, q);
        let dp = pressure[edge.from] - edge.to.map_or(supplies[k].unwrap_or(0.0), |j| pressure[j]);
        let flow = reference_flow(laws[k], edge.width_m, edge.length_m, h, dp);
        flows[k] = flow;
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
fn g3_supplied_newtonian_limit_recovers_native_pressure_flow_work_and_tangent() {
    let base = network();
    let viscosities = [0.003, 0.01, 0.02, 0.5];
    let supplies = [None, Some(230.0), Some(700.0), Some(-20_000.0)];
    let nonlinear = PressureDrivenPowerLawFilm::new(
        PowerLawFilm::new(
            &base,
            viscosities.map(|mu| law(mu, 1.0, 1.0)).to_vec(),
            controls(),
        )
        .unwrap(),
        &supplies,
    )
    .unwrap();
    let native = PressureDrivenFilm::new(
        base.with_channel_viscosities(&viscosities).unwrap(),
        &supplies,
    )
    .unwrap();
    let q = [1e-5, -2e-5, 5e-6];
    for v in [[2e-4, -1e-4, 0.0], [-2e-4, 1e-4, 0.0], [0.0; 3]] {
        let (mut force, mut pressure, mut channels) =
            ([0.0; 3], [0.0; 2], [FilmChannelReport::default(); 4]);
        let (mut ef, mut ep, mut ec) = ([0.0; 3], [0.0; 2], [FilmChannelReport::default(); 4]);
        let actual = nonlinear
            .evaluate_into(&q, &v, &mut force, &mut pressure, &mut channels)
            .unwrap();
        let expected = native
            .evaluate_into(&q, &v, &mut ef, &mut ep, &mut ec)
            .unwrap();
        assert_eq!(force, ef);
        assert_eq!(pressure, ep);
        assert_eq!(channels, ec);
        assert_eq!(actual.power_law.iterations, 0);
        assert_eq!(actual.mechanical_power_w, expected.mechanical_power_w);
        assert_eq!(actual.reservoir_power_w, expected.reservoir_power_w);
        assert_eq!(
            actual.power_law.film.dissipated_power_w,
            expected.film.dissipated_power_w
        );
        let (dq, dv, db) = (
            [2e-5, -1e-5, 3e-5],
            [1e-4, -2e-4, 0.0],
            [0.0, 13.0, -17.0, 1.0],
        );
        let (mut tangent, mut et) = ([0.0; 3], [0.0; 3]);
        nonlinear
            .tangent_into(&q, &v, &dq, &dv, &db, &mut tangent)
            .unwrap();
        native.tangent_into(&q, &v, &dq, &dv, &db, &mut et).unwrap();
        for (a, b) in tangent.into_iter().zip(et) {
            near(a, b, 1e-14, 3e-14);
        }
    }
}

#[test]
fn g1_pressure_feed_drives_thinning_and_thickening_flow_at_stationary_walls() {
    let mut base = one_channel();
    let mut passages = base.channels.to_vec();
    passages.push(passages[0].clone());
    base = ResistiveFilm::new(base.cells.to_vec(), passages, 1, 0.01, limits()).unwrap();
    for n in [0.5, 1.7, 3.0] {
        let material = law(0.04, n, 3.0);
        let original = PowerLawFilm::new(&base, vec![material; 2], controls()).unwrap();
        let expected_flow = reference_flow(material, W, L, H, 100.0);
        for (supplies, sign) in [
            ([Some(200.0), Some(0.0)], 1.0),
            ([Some(0.0), Some(200.0)], -1.0),
        ] {
            let supplied = PressureDrivenPowerLawFilm::new(original.clone(), &supplies).unwrap();
            let (mut f, mut p, mut channels) = ([0.0], [0.0], [FilmChannelReport::default(); 2]);
            let r = supplied
                .evaluate_into(&[0.0], &[0.0], &mut f, &mut p, &mut channels)
                .unwrap();
            near(p[0], 100.0, 1e-9, 2e-10);
            near(f[0], AREA * 100.0, 1e-13, 2e-10);
            near(
                channels[0].volume_flow_m3_s,
                -sign * expected_flow,
                1e-22,
                2e-10,
            );
            near(
                channels[1].volume_flow_m3_s,
                sign * expected_flow,
                1e-22,
                2e-10,
            );
            near(
                channels[0].volume_flow_m3_s + channels[1].volume_flow_m3_s,
                0.0,
                2e-23 + 4e-12 * expected_flow,
                0.0,
            );
            assert_eq!(r.mechanical_power_w, 0.0);
            near(r.reservoir_power_w, 200.0 * expected_flow, 1e-19, 3e-10);
            near(
                r.power_law.film.dissipated_power_w,
                r.reservoir_power_w,
                1e-20,
                3e-12,
            );
            assert!(r.power_law.iterations > 0);
        }
        let uniform = PressureDrivenPowerLawFilm::new(original, &[Some(-50.0); 2]).unwrap();
        let (mut f, mut p, mut channels) = ([0.0], [0.0], [FilmChannelReport::default(); 2]);
        let r = uniform
            .evaluate_into(&[0.0], &[0.0], &mut f, &mut p, &mut channels)
            .unwrap();
        assert_eq!(p, [-50.0]);
        assert_eq!(f, [-50.0 * AREA]);
        assert_eq!(r.power_law.iterations, 0);
        assert_eq!(r.power_law.film.dissipated_power_w, 0.0);
        assert_eq!(channels, [FilmChannelReport::default(); 2]);
    }
}

#[test]
fn g1_supplied_mixed_network_closes_each_cell_and_the_complete_active_work() {
    let base = network();
    let laws = mixed_laws();
    let supplies = [None, Some(180.0), Some(-90.0), Some(1e5)];
    let supplied = PressureDrivenPowerLawFilm::new(
        PowerLawFilm::new(&base, laws.clone(), controls()).unwrap(),
        &supplies,
    )
    .unwrap();
    let q = [1e-5, -2e-5, 5e-6];
    for target in [[400.0, -100.0], [-400.0, 100.0]] {
        let (v, expected_flows) = supplied_rates(&base, &laws, &q, target, &supplies);
        let (mut f, mut p, mut channels) = ([0.0; 3], [0.0; 2], [FilmChannelReport::default(); 4]);
        let r = supplied
            .evaluate_into(&q, &v, &mut f, &mut p, &mut channels)
            .unwrap();
        for (actual, expected) in p.into_iter().zip(target) {
            near(actual, expected, 2e-8, 3e-8);
        }
        for (actual, expected) in channels.iter().zip(expected_flows) {
            near(actual.volume_flow_m3_s, expected, 1e-22, 1e-8);
            assert!(actual.dissipated_power_w >= 0.0);
        }
        assert_eq!(channels[3], FilmChannelReport::default());
        near(
            channels[0].volume_flow_m3_s + channels[2].volume_flow_m3_s,
            2.0 * AREA * v[0],
            2e-22,
            5e-12,
        );
        near(
            channels[1].volume_flow_m3_s - channels[0].volume_flow_m3_s,
            1.5 * AREA * v[1],
            2e-22,
            5e-12,
        );
        let supply = -180.0 * channels[1].volume_flow_m3_s + 90.0 * channels[2].volume_flow_m3_s;
        near(r.reservoir_power_w, supply, 0.0, 1e-15);
        near(r.mechanical_power_w, dot(&v, &f), 0.0, 1e-15);
        near(
            r.mechanical_power_w + supply,
            r.power_law.film.dissipated_power_w,
            1e-20,
            5e-12,
        );
    }
    // A supply can do mechanical work, or receive work from a closing wall.
    let base = one_channel();
    let material = law(0.02, 0.5, 1.0);
    let supplied = PressureDrivenPowerLawFilm::new(
        PowerLawFilm::new(&base, vec![material], controls()).unwrap(),
        &[Some(1000.0)],
    )
    .unwrap();
    for dp in [-100.0, 100.0] {
        let flow = reference_flow(material, W, L, H, dp);
        let (mut f, mut p, mut channels) = ([0.0], [0.0], [FilmChannelReport::default()]);
        let r = supplied
            .evaluate_into(&[0.0], &[flow / AREA], &mut f, &mut p, &mut channels)
            .unwrap();
        near(p[0], 1000.0 + dp, 1e-8, 1e-10);
        assert_eq!(
            r.mechanical_power_w.is_sign_positive(),
            dp.is_sign_positive()
        );
        assert_eq!(
            r.reservoir_power_w.is_sign_positive(),
            dp.is_sign_negative()
        );
        near(
            r.mechanical_power_w + r.reservoir_power_w,
            flow * dp,
            1e-18,
            2e-10,
        );
        assert!(r.power_law.film.dissipated_power_w > 0.0);
    }
}

#[test]
fn g3_supplied_tangent_matches_geometry_rate_and_pressure_differences() {
    let base = network();
    let laws = mixed_laws();
    let supplies = [None, Some(180.0), Some(-90.0), Some(1e5)];
    let supplied = PressureDrivenPowerLawFilm::new(
        PowerLawFilm::new(&base, laws.clone(), controls()).unwrap(),
        &supplies,
    )
    .unwrap();
    let q = [1e-5, -2e-5, 5e-6];
    let (v, _) = supplied_rates(&base, &laws, &q, [40.0, -15.0], &supplies);
    let (dq, dv, db) = (
        [2e-5, -1e-5, 3e-5],
        [1e-4, -2e-4, 0.0],
        [0.0, 13.0, -17.0, 1.0],
    );
    let mut tangent = [0.0; 3];
    supplied
        .tangent_into(&q, &v, &dq, &dv, &db, &mut tangent)
        .unwrap();
    let epsilon = 1e-4;
    let mut forces = [[0.0; 3]; 2];
    for (index, sign) in [-1.0, 1.0].into_iter().enumerate() {
        let boundary: Vec<_> = supplies
            .iter()
            .zip(db)
            .map(|(p, d)| p.map(|p| p + sign * epsilon * d))
            .collect();
        let model = supplied.with_reservoir_pressures(&boundary).unwrap();
        let shifted_q = std::array::from_fn::<_, 3, _>(|i| q[i] + sign * epsilon * dq[i]);
        let shifted_v = std::array::from_fn::<_, 3, _>(|i| v[i] + sign * epsilon * dv[i]);
        model
            .evaluate_force_into(&shifted_q, &shifted_v, &mut forces[index], &mut [0.0; 2])
            .unwrap();
    }
    for i in 0..3 {
        near(
            tangent[i],
            (forces[1][i] - forces[0][i]) / (2.0 * epsilon),
            2e-9,
            2e-5,
        );
    }
    let mut bias = [0.0; 3];
    supplied
        .tangent_into(
            &q,
            &v,
            &[0.0; 3],
            &[0.0; 3],
            &[0.0, 1.0, 1.0, 1.0],
            &mut bias,
        )
        .unwrap();
    for (a, b) in bias.into_iter().zip([2.0 * AREA, 1.5 * AREA, 0.0]) {
        near(a, b, 1e-17, 2e-13);
    }
    let (a, b) = ([0.01, -0.003, 0.0], [-0.005, 0.002, 0.0]);
    let (mut ta, mut tb) = ([0.0; 3], [0.0; 3]);
    supplied
        .tangent_into(&q, &v, &[0.0; 3], &a, &[0.0; 4], &mut ta)
        .unwrap();
    supplied
        .tangent_into(&q, &v, &[0.0; 3], &b, &[0.0; 4], &mut tb)
        .unwrap();
    near(dot(&a, &tb), dot(&b, &ta), 1e-12, 4e-13);
    assert!(dot(&a, &ta) > 0.0);
}

#[test]
fn g3_ambient_recovery_and_rheology_refresh_keep_the_original_owner() {
    let base = network();
    let native = PowerLawFilm::new(&base, mixed_laws(), controls()).unwrap();
    let ambient = [None, Some(0.0), Some(0.0), Some(0.0)];
    let supplied = PressureDrivenPowerLawFilm::new(native.clone(), &ambient).unwrap();
    let (q, v) = ([1e-5, -2e-5, 5e-6], [1e-4, -2e-4, 0.0]);
    let (mut f, mut p, mut c) = ([0.0; 3], [0.0; 2], [FilmChannelReport::default(); 4]);
    let (mut ef, mut ep, mut ec) = ([0.0; 3], [0.0; 2], [FilmChannelReport::default(); 4]);
    let r = supplied
        .evaluate_into(&q, &v, &mut f, &mut p, &mut c)
        .unwrap();
    let expected = native
        .evaluate_with_channels_into(&q, &v, &mut ef, &mut ep, &mut ec)
        .unwrap();
    assert_eq!(f, ef);
    assert_eq!(p, ep);
    assert_eq!(c, ec);
    assert_eq!(r.power_law.iterations, expected.iterations);
    assert_eq!(r.power_law.backtracks, expected.backtracks);
    assert_eq!(
        r.power_law.maximum_flow_defect_m3_s,
        expected.maximum_flow_defect_m3_s
    );
    assert_eq!(
        r.power_law.film.dissipated_power_w,
        expected.film.dissipated_power_w
    );
    let (mut tangent, mut et) = ([0.0; 3], [0.0; 3]);
    supplied
        .tangent_into(&q, &v, &[0.2; 3], &[0.1; 3], &[0.0; 4], &mut tangent)
        .unwrap();
    native
        .tangent_into(&q, &v, &[0.2; 3], &[0.1; 3], &mut et)
        .unwrap();
    assert_eq!(tangent, et);
    let active = supplied
        .with_reservoir_pressures(&[None, Some(180.0), Some(-90.0), Some(0.0)])
        .unwrap();
    let laws = vec![law(0.04, 1.7, 1.0); 4];
    let refreshed = active.with_laws(laws.clone()).unwrap();
    let independent = PressureDrivenPowerLawFilm::new(
        PowerLawFilm::new(&base, laws.clone(), controls()).unwrap(),
        active.reservoir_pressures_pa(),
    )
    .unwrap();
    assert_eq!(refreshed.power_law_film().laws(), laws);
    assert_eq!(
        refreshed.reservoir_pressures_pa(),
        active.reservoir_pressures_pa()
    );
    refreshed
        .evaluate_into(&q, &v, &mut f, &mut p, &mut c)
        .unwrap();
    independent
        .evaluate_into(&q, &v, &mut ef, &mut ep, &mut ec)
        .unwrap();
    assert_eq!(f, ef);
    assert_eq!(p, ep);
    assert_eq!(c, ec);
}

#[test]
fn g4_supplied_pressure_domains_budgets_and_checks_refuse_atomically() {
    let base = network();
    let original = PowerLawFilm::new(&base, mixed_laws(), controls()).unwrap();
    for supplies in [
        vec![None; 4],
        vec![Some(0.0); 4],
        vec![None, Some(f64::NAN), Some(0.0), Some(0.0)],
        vec![
            None,
            Some(limits().maximum_pressure_pa * 2.0),
            Some(0.0),
            Some(0.0),
        ],
    ] {
        assert!(PressureDrivenPowerLawFilm::new(original.clone(), &supplies).is_err());
    }
    let model =
        PressureDrivenPowerLawFilm::new(original, &[None, Some(180.0), Some(-90.0), Some(0.0)])
            .unwrap();
    let (q, v) = ([1e-5, -2e-5, 5e-6], [1e-4, -2e-4, 0.0]);
    let mut checks = 0;
    model
        .evaluate_with_check_into(
            &q,
            &v,
            &mut [0.0; 3],
            &mut [0.0; 2],
            &mut [FilmChannelReport::default(); 4],
            || {
                checks += 1;
                Ok(())
            },
        )
        .unwrap();
    for stop in [1, 2, checks / 2, checks] {
        let (mut f, mut p, mut c) = (
            [71.0; 3],
            [72.0; 2],
            [FilmChannelReport {
                volume_flow_m3_s: 73.0,
                dissipated_power_w: 74.0,
            }; 4],
        );
        let before = (f, p, c);
        let mut calls = 0;
        let error = model
            .evaluate_with_check_into(&q, &v, &mut f, &mut p, &mut c, || {
                calls += 1;
                if calls == stop {
                    Err(FilmError("cancelled supplied film"))
                } else {
                    Ok(())
                }
            })
            .unwrap_err();
        assert_eq!(error, FilmError("cancelled supplied film"));
        assert_eq!((f, p, c), before);
    }
    let mut tangent = [75.0; 3];
    for db in [[1.0, 0.0, 0.0, 0.0], [0.0, f64::NAN, 0.0, 0.0]] {
        assert!(
            model
                .tangent_into(&q, &v, &[0.0; 3], &[0.0; 3], &db, &mut tangent)
                .is_err()
        );
        assert_eq!(tangent, [75.0; 3]);
    }
    let base = one_channel();
    let domain = PressureDrivenPowerLawFilm::new(
        PowerLawFilm::new(&base, vec![law(0.02, 0.5, 1.0)], controls()).unwrap(),
        &[Some(limits().maximum_pressure_pa)],
    )
    .unwrap();
    let (mut f, mut p) = ([81.0], [82.0]);
    assert!(
        domain
            .evaluate_force_into(&[0.0], &[1e-4], &mut f, &mut p)
            .is_err()
    );
    assert_eq!((f, p), ([81.0], [82.0]));
    let mut small_budget = controls();
    small_budget.maximum_iterations = 1;
    let limited = PressureDrivenPowerLawFilm::new(
        PowerLawFilm::new(&base, vec![law(0.02, 0.5, 1.0)], small_budget).unwrap(),
        &[Some(100.0)],
    )
    .unwrap();
    assert!(
        limited
            .evaluate_force_into(&[0.0], &[1e-3], &mut f, &mut p)
            .is_err()
    );
    assert_eq!((f, p), ([81.0], [82.0]));
}
