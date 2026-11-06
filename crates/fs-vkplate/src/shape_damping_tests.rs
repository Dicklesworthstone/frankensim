//! Measured shape damping on a head (`Body::attach_shape_damping`, MB 7247).
use super::*;

const RADIUS: f64 = 0.15;
const TENSION: f64 = 2000.0;
const THICKNESS: f64 = 2.54e-4;
const DENSITY: f64 = 1390.0;

fn head() -> Body {
    Body::plate_without_surface_loss("shape head", &PlateSpec {
        radius: RADIUS, hole: 0.0, grid: 0.01, edge: Edge::Clamped, thickness: &|_| THICKNESS,
        rest_height: &|_| 0.0, young: 4e9, poisson: 0.38, density: DENSITY, tension: TENSION,
        bending: true, nonlinear: Nonlinear::None, t60: [(100.0, 1.0), (1000.0, 1.0)],
        support: None, muffle: None,
    }).unwrap()
}

/// The (0,1) shape on the head's nodes and its in-vacuo angular frequency (Rayleigh quotient).
fn first_shape(body: &Body) -> (Vec<f64>, f64) {
    let shape = radiation::ModeShape::new(0, 1);
    let phi: Vec<f64> = body.xy.iter().map(|p| shape.value(p[0] / RADIUS, p[1] / RADIUS, false)).collect();
    let mut k_phi = vec![0.0; body.n];
    body.stiff.mul_add(&phi, 1.0, &mut k_phi);
    body.lap.mul_add(&phi, body.tension, &mut k_phi);
    let modal_mass: f64 = phi.iter().zip(&body.mass).map(|(a, w)| a * a * w).sum();
    (phi.clone(), (dot(&phi, &k_phi) / modal_mass).sqrt())
}

#[test]
fn shape_damping_decays_its_shape_at_the_added_rate_and_closes_energy() {
    let (dt, zeta) = (1.0 / 96_000.0, 0.02);
    let mut body = head();
    let mut bare = head();
    let receipt = body.attach_shape_damping(RADIUS, &[(0, 1, zeta), (1, 1, 0.01)], "test rows").unwrap();
    assert!(receipt[0].contains("(0,1) ") && receipt[0].contains("(1,1)c") && receipt[0].contains("(1,1)s"), "{receipt:?}");
    assert!(receipt[1].contains("largest off-diagonal Gram term") && receipt[1].contains("3x3"), "{receipt:?}");
    let (phi, omega) = first_shape(&body);
    let membrane = 2.404826 * (TENSION / (DENSITY * THICKNESS)).sqrt() / RADIUS;
    assert!((omega / membrane - 1.0).abs() < 0.03, "in-vacuo (0,1) {omega:.1} vs membrane {membrane:.1} rad/s");
    let damping = body.shape_damping.as_ref().unwrap();
    let largest = damping.two_sigma.iter().cloned().fold(0.0_f64, f64::max);
    assert!((damping.two_sigma[0] - 2.0 * zeta * omega).abs() <= 1e-12 * omega);
    assert!(damping.rate_bound >= largest && damping.rate_bound <= 1.01 * largest, "{} vs {largest}", damping.rate_bound);
    // from rest in the (0,1) shape
    for b in [&mut body, &mut bare] {
        for i in 0..b.n {
            b.w[i] = 1e-4 * phi[i];
            b.wp[i] = b.w[i];
        }
        b.sync();
    }
    let initial = body.discrete_energy(dt);
    assert!(initial > 0.0 && (bare.discrete_energy(dt) - initial).abs() <= 1e-15 * initial);
    // the shape's own coordinate, q = phi' M w
    let coordinate_of = |b: &Body| -> f64 { b.w.iter().zip(&b.mass).zip(&phi).map(|((w, m), p)| w * m * p).sum() };
    let mut coordinate = vec![coordinate_of(&body)];
    let mut previous = initial;
    let steps = 4800;
    for _ in 0..steps {
        for b in [&mut body, &mut bare] {
            b.linear_forces();
            b.predict(dt);
            b.finish(dt);
            b.rotate();
        }
        coordinate.push(coordinate_of(&body));
        let energy = body.discrete_energy(dt);
        assert!(energy >= 0.0 && energy <= previous + 2e-12 * initial);
        assert!((energy + body.shape_damping_dissipated_j() - initial).abs() < 2e-10 * initial,
            "E {energy:e} + D {:e} vs E0 {initial:e}", body.shape_damping_dissipated_j());
        previous = energy;
    }
    // without the row nothing is lost
    assert!((bare.discrete_energy(dt) - initial).abs() < 1e-9 * initial);
    assert_eq!(bare.shape_damping_dissipated_j(), 0.0);
    // With it the shape's coordinate falls at sigma = zeta omega, read as the least-squares slope of ln|q| through its
    // half-cycle peaks. The rate is not read on the energy: the sampled Bessel shape is not exactly a grid mode, and its
    // remainder in higher grid modes, which the row barely touches, keeps about 2 % of the start energy on this 1 cm
    // grid, so the energy's mean rate over 0.05 s reads low (44.3 against 47.9 1/s when first run).
    let quarter = (std::f64::consts::FRAC_PI_2 / omega / dt).round() as usize;
    let peaks: Vec<(f64, f64)> = (0..coordinate.len() - quarter)
        .filter(|&i| coordinate[i.saturating_sub(quarter)..=i + quarter].iter().all(|x| x.abs() <= coordinate[i].abs()))
        .map(|i| (i as f64 * dt, coordinate[i].abs().ln()))
        .collect();
    assert!(peaks.len() >= 15, "{} half-cycle peaks", peaks.len());
    let (mean_t, mean_y) = (
        peaks.iter().map(|p| p.0).sum::<f64>() / peaks.len() as f64,
        peaks.iter().map(|p| p.1).sum::<f64>() / peaks.len() as f64,
    );
    let slope = peaks.iter().map(|p| (p.0 - mean_t) * (p.1 - mean_y)).sum::<f64>()
        / peaks.iter().map(|p| (p.0 - mean_t).powi(2)).sum::<f64>();
    assert!((-slope / (zeta * omega) - 1.0).abs() < 0.03, "the shape falls at {:.2} 1/s, expected {:.2}", -slope, zeta * omega);
}

#[test]
fn shape_damping_refuses_bad_rows_a_second_attach_and_untensioned_plates() {
    let mut body = head();
    let cases: [&[(usize, usize, f64)]; 8] = [
        &[],
        &[(0, 1, 0.0)],
        &[(0, 1, 1.0)],
        &[(0, 1, -0.01)],
        &[(0, 1, f64::NAN)],
        &[(0, 0, 0.01)],
        &[(0, 1, 0.01), (0, 1, 0.02)],
        // j_06 = 18.07 exceeds this grid's 2 pi a / (6 h) = 15.7
        &[(0, 6, 0.01)],
    ];
    for rows in cases {
        assert!(body.attach_shape_damping(RADIUS, rows, "bad").is_err(), "{rows:?}");
        assert!(body.shape_damping.is_none());
    }
    assert!(body.attach_shape_damping(-RADIUS, &[(0, 1, 0.01)], "bad").is_err());
    body.attach_shape_damping(RADIUS, &[(0, 1, 0.01)], "once").unwrap();
    assert!(body.attach_shape_damping(RADIUS, &[(1, 1, 0.01)], "twice").is_err());
    let mut plate = Body::plate_without_surface_loss("plate", &PlateSpec {
        radius: 0.2, hole: 0.0, grid: 0.02, edge: Edge::Free, thickness: &|_| 1e-3, rest_height: &|_| 0.0,
        young: 1e11, poisson: 0.33, density: 8600.0, tension: 0.0, bending: true, nonlinear: Nonlinear::None,
        t60: [(100.0, 1.0), (1000.0, 1.0)], support: None, muffle: None,
    }).unwrap();
    assert!(plate.attach_shape_damping(0.2, &[(0, 1, 0.01)], "plate").is_err());
}

#[test]
#[should_panic(expected = "measured shape damping")]
fn exported_stepper_refuses_shape_damping() {
    let mut body = head();
    body.attach_shape_damping(RADIUS, &[(0, 1, 0.01)], "export").unwrap();
    let mut instrument = Instrument::new(96_000.0);
    instrument.add(body);
    let _ = instrument.image();
}
