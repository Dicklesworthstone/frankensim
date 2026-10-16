use super::*;

fn gap(h: f64, b: &[f64]) -> GapPort {
    GapPort {
        reference_m: h,
        closure: b.to_vec(),
    }
}
fn limits() -> FilmLimits {
    FilmLimits {
        minimum_cell_gap_m: 1e-6,
        maximum_gap_m: 0.01,
        maximum_pressure_pa: 1e6,
    }
}
fn single_parts() -> (Vec<FilmCell>, Vec<FilmChannel>) {
    let g = gap(1e-3, &[1.0, -1.0]);
    (
        vec![FilmCell {
            area_m2: 0.002,
            gap: g.clone(),
        }],
        vec![FilmChannel {
            from: 0,
            to: None,
            width_m: 0.05,
            length_m: 0.01,
            gap: g,
        }],
    )
}
fn single() -> ResistiveFilm {
    let (cells, channels) = single_parts();
    ResistiveFilm::new(cells, channels, 2, 1.8e-5, limits()).unwrap()
}
fn network_parts() -> (Vec<FilmCell>, Vec<FilmChannel>) {
    let a = gap(0.001, &[1.0, -1.0, 0.3]);
    let b = gap(0.0014, &[0.7, -0.7, -0.2]);
    (
        vec![
            FilmCell {
                area_m2: 0.002,
                gap: a.clone(),
            },
            FilmCell {
                area_m2: 0.003,
                gap: b.clone(),
            },
        ],
        vec![
            FilmChannel {
                from: 0,
                to: None,
                width_m: 0.04,
                length_m: 0.01,
                gap: a,
            },
            FilmChannel {
                from: 1,
                to: None,
                width_m: 0.05,
                length_m: 0.02,
                gap: b,
            },
            FilmChannel {
                from: 0,
                to: Some(1),
                width_m: 0.03,
                length_m: 0.012,
                gap: gap(0.0012, &[0.85, -0.85, 0.05]),
            },
        ],
    )
}
fn network() -> ResistiveFilm {
    let (cells, channels) = network_parts();
    ResistiveFilm::new(cells, channels, 3, 1.8e-5, limits()).unwrap()
}
fn near(a: f64, b: f64, tol: f64) {
    assert!(
        (a - b).abs() <= tol * (a.abs() + b.abs()).max(1e-12),
        "{a} != {b}"
    );
}

#[test]
fn exact_poiseuille_pressure_equal_opposite_reaction_and_power() {
    let film = single();
    let mut f = [0.0; 2];
    let mut p = [0.0; 1];
    let r = film
        .evaluate_into(&[0.0; 2], &[0.02, 0.0], &mut f, &mut p)
        .unwrap();
    let conductance = 0.05 * 1e-9 / (12.0 * 1.8e-5 * 0.01);
    near(p[0], 0.002 * 0.02 / conductance, 1e-13);
    near(f[0], 0.002 * p[0], 1e-13);
    assert_eq!(f[0], -f[1]);
    near(r.dissipated_power_w, 0.02 * f[0], 1e-13);
}

#[test]
fn common_translation_has_no_air_force_and_opening_preserves_suction() {
    let film = single();
    let mut f = [99.0; 2];
    let mut p = [99.0; 1];
    let r = film
        .evaluate_into(&[0.0; 2], &[0.3, 0.3], &mut f, &mut p)
        .unwrap();
    assert_eq!(f, [0.0; 2]);
    assert_eq!(p, [0.0; 1]);
    assert_eq!(r.dissipated_power_w, 0.0);
    film.evaluate_into(&[0.0; 2], &[-0.02, 0.0], &mut f, &mut p)
        .unwrap();
    assert!(p[0] < 0.0);
    assert!(f[0] < 0.0);
    assert_eq!(f[0], -f[1]);
}

#[test]
fn halving_hydraulic_gap_multiplies_resistance_by_eight() {
    let film = single();
    let mut f = [0.0; 2];
    let mut p = [0.0; 1];
    film.evaluate_into(&[0.0; 2], &[0.01, 0.0], &mut f, &mut p)
        .unwrap();
    let initial = f[0];
    film.evaluate_into(&[0.0005, 0.0], &[0.01, 0.0], &mut f, &mut p)
        .unwrap();
    near(f[0], 8.0 * initial, 1e-13);
}

#[test]
fn exact_tangent_includes_gap_and_independent_effort_directions() {
    let film = network();
    let q = [0.0001, -0.00005, 0.00002];
    let v = [0.01, -0.03, 0.02];
    let dq = [0.0002, -0.0001, 0.00003];
    let dv = [0.04, -0.01, 0.03];
    let mut analytic = [0.0; 3];
    film.tangent_into(&q, &v, &dq, &dv, &mut analytic).unwrap();
    let epsilon = 1e-5;
    let mut plus = [0.0; 3];
    let mut minus = [0.0; 3];
    let mut pressure = [0.0; 2];
    let qp = std::array::from_fn::<_, 3, _>(|i| q[i] + epsilon * dq[i]);
    let qm = std::array::from_fn::<_, 3, _>(|i| q[i] - epsilon * dq[i]);
    let vp = std::array::from_fn::<_, 3, _>(|i| v[i] + epsilon * dv[i]);
    let vm = std::array::from_fn::<_, 3, _>(|i| v[i] - epsilon * dv[i]);
    film.evaluate_into(&qp, &vp, &mut plus, &mut pressure)
        .unwrap();
    film.evaluate_into(&qm, &vm, &mut minus, &mut pressure)
        .unwrap();
    for i in 0..3 {
        near(analytic[i], (plus[i] - minus[i]) / (2.0 * epsilon), 2e-8);
    }
}

#[test]
fn distributed_pressure_is_reciprocal_and_passive_not_diagonal_modal_loss() {
    let film = network();
    let q = [0.0001, 0.0, 0.0];
    let mut matrix = [[0.0; 3]; 3];
    for j in 0..3 {
        let mut v = [0.0; 3];
        v[j] = 1.0;
        let mut p = [0.0; 2];
        film.evaluate_into(&q, &v, &mut matrix[j], &mut p).unwrap();
    }
    for i in 0..3 {
        for j in 0..3 {
            near(matrix[i][j], matrix[j][i], 1e-12);
        }
    }
    assert!(matrix[0][1] < 0.0);
    assert!(matrix[0][2].abs() > 0.0);
    for v in [[0.02, -0.01, 0.03], [-0.02, 0.03, -0.01], [0.1, 0.1, 0.0]] {
        let mut f = [0.0; 3];
        let mut p = [0.0; 2];
        let r = film.evaluate_into(&q, &v, &mut f, &mut p).unwrap();
        let work = v.iter().zip(f).map(|(v, f)| v * f).sum();
        assert!(r.dissipated_power_w >= 0.0);
        near(work, r.dissipated_power_w, 1e-12);
    }
}

#[test]
fn contact_shuts_a_passage_without_inventing_a_gap_floor() {
    let (cells, mut channels) = network_parts();
    channels[2].gap = gap(-1e-4, &[0.0; 3]);
    let edge = &channels[0];
    let conductance = edge.width_m * edge.gap.reference_m.powi(3) / (12.0 * 1.8e-5 * edge.length_m);
    let expected_pressure = cells[0].area_m2 * 0.01 / conductance;
    let film = ResistiveFilm::new(cells, channels, 3, 1.8e-5, limits()).unwrap();
    let mut f = [0.0; 3];
    let mut p = [0.0; 2];
    film.evaluate_into(&[0.0; 3], &[0.01, 0.0, 0.0], &mut f, &mut p)
        .unwrap();
    near(p[0], expected_pressure, 1e-12);
    let mut tangent = [0.0; 3];
    film.tangent_into(
        &[0.0; 3],
        &[0.01, 0.0, 0.0],
        &[0.0; 3],
        &[0.01, 0.0, 0.0],
        &mut tangent,
    )
    .unwrap();
    for i in 0..3 {
        near(tangent[i], f[i], 1e-12);
    }
}

#[test]
fn trapped_gas_refuses_instead_of_creating_compressibility_or_leakage() {
    let (cells, mut channels) = network_parts();
    channels[0].gap.reference_m = 0.0;
    channels[1].gap.reference_m = 0.0;
    let film = ResistiveFilm::new(cells, channels, 3, 1.8e-5, limits()).unwrap();
    let mut f = [91.0; 3];
    let mut p = [92.0; 2];
    assert!(
        film.evaluate_into(&[0.0; 3], &[0.01, 0.0, 0.0], &mut f, &mut p)
            .is_err()
    );
    assert_eq!(f, [91.0; 3]);
    assert_eq!(p, [92.0; 2]);
}

#[test]
fn pressure_gap_nonfinite_and_dimension_refusals_preserve_outputs() {
    let film = single();
    let mut f = [31.0; 2];
    let mut p = [41.0; 1];
    for (q, v) in [
        ([0.002, 0.0], [0.0; 2]),
        ([f64::NAN, 0.0], [0.0; 2]),
        ([0.0; 2], [f64::INFINITY, 0.0]),
        ([0.0; 2], [1e9, 0.0]),
    ] {
        assert!(film.evaluate_into(&q, &v, &mut f, &mut p).is_err());
        assert_eq!(f, [31.0; 2]);
        assert_eq!(p, [41.0; 1]);
    }
    assert!(film.evaluate_into(&[], &[0.0; 2], &mut f, &mut p).is_err());
    assert!(
        film.tangent_into(&[0.0; 2], &[0.0; 2], &[f64::NAN, 0.0], &[0.0; 2], &mut f)
            .is_err()
    );
    assert_eq!(f, [31.0; 2]);
    assert_eq!(p, [41.0; 1]);
}

#[test]
fn construction_rejects_bad_geometry_without_normalizing_it() {
    let (cells, channels) = single_parts();
    let make = |cells, edges, mu| ResistiveFilm::new(cells, edges, 2, mu, limits());
    assert!(make(cells.clone(), channels.clone(), 0.0).is_err());
    let mut bad = cells.clone();
    bad[0].area_m2 = -1.0;
    assert!(make(bad, channels.clone(), 1.8e-5).is_err());
    let mut bad = channels.clone();
    bad[0].to = Some(99);
    assert!(make(cells.clone(), bad, 1.8e-5).is_err());
    let mut bad = channels.clone();
    bad[0].gap.closure.pop();
    assert!(make(cells.clone(), bad, 1.8e-5).is_err());
}
