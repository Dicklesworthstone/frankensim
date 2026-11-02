//! Known-answer tests of the closed-body bake (pexfp.378), with planted
//! negatives asserted to FAIL the same gates:
//! - T1: a rigid sphere with two polar caps, in phase and in anti-phase,
//!   against the exact series; the tolerance is the route's own mesh ladder.
//!   N1: today's baffled microphone on the same body must fail it.
//! - T1b: the closed cylinder's monopole limit and its dipole slope in ka.
//! - T1c: the card drums against the independent numpy route
//!   (`closedbody.py`, hard-body CBIE by reciprocity) per line within 1 dB.
//!   N2: a sign-flipped resonant head, or the image body dropped (the floor
//!   kept only as the microphone's image), must miss it.
//!
//! Gas for every test: rho 1.2041 kg/m^3, c 343.2 m/s (`closedbody.py`'s).
use super::*;

const GAS: (f64, f64) = (1.2041, 343.2);

fn c(re: f64, im: f64) -> C64 {
    C64::new(re, im)
}

fn rel(a: C64, b: C64) -> f64 {
    (a - b).abs() / b.abs()
}

// ------------------------------------------------------------------ T1: sphere

/// A lat-long sphere whose rings fall exactly on both cap edges. Tag +1 top
/// cap (theta < cap), -1 bottom cap (theta > pi - cap), 0 between.
fn sphere(radius: f64, cap: f64, cap_bands: usize, mid_bands: usize, sectors: usize) -> (Vec<[[f64; 3]; 3]>, Vec<i8>) {
    let mut theta = vec![0.0];
    for k in 1..=cap_bands {
        theta.push(cap * k as f64 / cap_bands as f64);
    }
    for k in 1..=mid_bands {
        theta.push(cap + (PI - 2.0 * cap) * k as f64 / mid_bands as f64);
    }
    for k in 1..=cap_bands {
        theta.push(if k == cap_bands { PI } else { PI - cap + cap * k as f64 / cap_bands as f64 });
    }
    let ring = |t: f64| -> Vec<[f64; 3]> {
        let (st, ct) = det::sin_cos(t);
        circle(sectors).into_iter().map(|(c, s)| [radius * st * c, radius * st * s, radius * ct]).collect()
    };
    let bands = theta.len() - 1;
    let tag = |band: usize| -> i8 {
        if band < cap_bands {
            1
        } else if band >= cap_bands + mid_bands {
            -1
        } else {
            0
        }
    };
    let mut tris = Vec::new();
    let mut tags = Vec::new();
    let mut push = |mut t: [[f64; 3]; 3], g: i8| {
        let out = centroid(&t);
        if dot(cross(sub(t[1], t[0]), sub(t[2], t[0])), out) < 0.0 {
            t.swap(1, 2);
        }
        tris.push(t);
        tags.push(g);
    };
    let north = [0.0, 0.0, radius];
    let south = [0.0, 0.0, -radius];
    let mut previous = ring(theta[1]);
    for j in 0..sectors {
        push([north, previous[j], previous[(j + 1) % sectors]], tag(0));
    }
    for band in 1..bands - 1 {
        let next = ring(theta[band + 1]);
        for j in 0..sectors {
            let j1 = (j + 1) % sectors;
            push([previous[j], previous[j1], next[j1]], tag(band));
            push([previous[j], next[j1], next[j]], tag(band));
        }
        previous = next;
    }
    for j in 0..sectors {
        push([south, previous[(j + 1) % sectors], previous[j]], tag(bands - 1));
    }
    (tris, tags)
}

fn legendre(nmax: usize, x: f64) -> Vec<f64> {
    let mut p = vec![1.0, x];
    for n in 1..nmax {
        p.push(((2 * n + 1) as f64 * x * p[n] - n as f64 * p[n - 1]) / (n + 1) as f64);
    }
    p
}

/// Outgoing spherical Hankel functions h_n(x) and their derivatives, n <= nmax.
fn hankel(nmax: usize, x: f64) -> (Vec<C64>, Vec<C64>) {
    let (s, co) = det::sin_cos(x);
    let e = c(co, s);
    // h0 = -i e^{ix} / x, h1 = -e^{ix} (x + i) / x^2
    let mut h = vec![c(0.0, -1.0 / x) * e, (e * c(x, 1.0)).scale(-1.0 / (x * x))];
    for n in 1..nmax {
        let next = h[n].scale((2 * n + 1) as f64 / x) - h[n - 1];
        h.push(next);
    }
    let mut d = vec![h[1].scale(-1.0)];
    for n in 1..=nmax {
        d.push(h[n - 1] - h[n].scale((n + 1) as f64 / x));
    }
    (h, d)
}

/// Rigid sphere of radius `r0`: outward velocity 1 on theta < cap and `s` on
/// theta > pi - cap; pressure at (r, theta), exp(-i omega t).
fn sphere_exact(r0: f64, cap: f64, s: f64, k: f64, r: f64, theta: f64) -> C64 {
    let nmax = 60;
    let pc = legendre(nmax + 1, det::cos(cap));
    let (_, dh) = hankel(nmax, k * r0);
    let (h, _) = hankel(nmax, k * r);
    let pt = legendre(nmax, det::cos(theta));
    let omega = k * GAS.1;
    let mut total = C64::ZERO;
    for n in 0..=nmax {
        let top = if n == 0 { 1.0 - det::cos(cap) } else { (pc[n - 1] - pc[n + 1]) / (2 * n + 1) as f64 };
        let parity = if n % 2 == 0 { 1.0 } else { -1.0 };
        let vn = (2 * n + 1) as f64 / 2.0 * top * (1.0 + s * parity);
        if vn == 0.0 {
            continue;
        }
        // dp/dr = i omega rho v_r on r = r0
        let a = c(0.0, omega * GAS.0 * vn) * (c(1.0, 0.0) * dh[n].conj()).scale(1.0 / (k * dh[n].abs() * dh[n].abs()));
        total = total + a * h[n].scale(pt[n]);
    }
    total
}

#[test]
fn t1_sphere_caps_match_the_exact_series_within_their_own_ladder_and_n1_fails() {
    let (r0, cap, r_obs) = (0.1, 40.0_f64.to_radians(), 0.6);
    let rungs = [(4, 10, 24), (6, 15, 36)];
    let ids = vec![
        SourceId::StructuralMode("caps in phase".into()),
        SourceId::StructuralMode("caps in anti-phase".into()),
    ];
    let angles = [0.0, 150.0_f64.to_radians()];
    let receivers: Vec<[f64; 3]> = angles.iter().map(|&t| [r_obs * det::sin(t), 0.0, r_obs * det::cos(t)]).collect();
    let mut worst = (0.0_f64, String::new());
    let mut n1 = None;
    for ka in [0.1, 0.5, 1.0, 2.0] {
        let k = ka / r0;
        let omega = k * GAS.1;
        let mut results = Vec::new();
        let mut fine_mesh = None;
        for &(cb, mb, sectors) in &rungs {
            let (tris, tags) = sphere(r0, cap, cb, mb, sectors);
            let traces: Vec<Vec<f64>> = [1.0, -1.0]
                .iter()
                .map(|&s| tags.iter().map(|&g| if g == 1 { 1.0 } else if g == -1 { s } else { 0.0 }).collect())
                .collect();
            let surface = SpherePanels::from_triangles(tris.clone()).unwrap();
            results.push(forward(&surface, GAS, &ids, &traces, &receivers, omega).unwrap());
            fine_mesh = Some((tris, traces));
        }
        let (tris, traces) = fine_mesh.unwrap();
        let surface = SpherePanels::from_triangles(tris.clone()).unwrap();
        let mesh = Mesh { part: vec![Part::Batter; tris.len()], body: vec![0; tris.len()], tris };
        for (si, s) in [1.0, -1.0].into_iter().enumerate() {
            for (ri, &theta) in angles.iter().enumerate() {
                let exact = sphere_exact(r0, cap, s, k, r_obs, theta);
                let (coarse, fine) = (results[0][si][ri], results[1][si][ri]);
                let ladder = rel(fine, coarse);
                let error = rel(fine, exact);
                let what = format!(
                    "ka {ka} {} {:.0} deg: error {error:.2e} ({:+.3} dB), ladder {ladder:.2e}",
                    if s > 0.0 { "in" } else { "anti" },
                    theta.to_degrees(),
                    db(fine) - db(exact)
                );
                eprintln!("T1 {what}");
                assert!(error <= 2.0 * ladder, "T1 error exceeds its own ladder: {what}");
                assert!(error <= 0.0116, "T1 error above 0.1 dB / 0.66 deg: {what}");
                if error > worst.0 {
                    worst = (error, what);
                }
                // N1: today's baffled microphone on the same panels and tolerance.
                let today = baffled(&mesh, &surface, GAS, &traces[si], &[receivers[ri]], omega);
                if ka == 0.1 && s < 0.0 && theta == 0.0 {
                    n1 = Some((rel(today, exact), db(today) - db(exact), 2.0 * ladder));
                }
            }
        }
    }
    eprintln!("T1 worst: {}", worst.1);
    let (error, level, tolerance) = n1.unwrap();
    eprintln!("N1 baffled anti-phase on axis at ka 0.1: {level:+.1} dB, error {error:.2e} against tolerance {tolerance:.2e}");
    assert!(error > tolerance && error > 0.0116, "N1 must fail T1's gate");
    assert!(level > 6.0, "N1: the baffled microphone must read the anti-phase pair far too loud");
}

// ------------------------------------------------------------------ T1b: cylinder limits

#[test]
fn t1b_closed_cylinder_monopole_limit_and_dipole_slope() {
    let (radius, depth) = (0.1703, 0.1524);
    let mesh = drum_body(radius, depth, SECTORS, CAP_RINGS);
    let surface = SpherePanels::from_triangles(mesh.tris.clone()).unwrap();
    let pistons = |s: f64| trace(&mesh, radius, &|part, _, _| if part == Part::Batter { 1.0 } else { s });
    let traces = vec![pistons(1.0), pistons(-1.0)];
    let ids = vec![
        SourceId::StructuralMode("both heads out".into()),
        SourceId::StructuralMode("heads together".into()),
    ];
    let q: f64 = traces[0].iter().zip(surface.areas()).map(|(v, a)| v * a).sum();
    let centre = [0.0, 0.0, -0.5 * depth];
    let direction = [0.48, 0.6, 0.64];
    let r = 100.0;
    let x = [centre[0] + r * direction[0], centre[1] + r * direction[1], centre[2] + r * direction[2]];
    let mut ratio = Vec::new();
    for ka in [0.05, 0.1] {
        let k = ka / radius;
        let omega = k * GAS.1;
        let p = forward(&surface, GAS, &ids, &traces, &[x], omega).unwrap();
        let (s, co) = det::sin_cos(k * r);
        let monopole = c(0.0, -omega * GAS.0 * q) * c(co, s).scale(1.0 / (4.0 * PI * r));
        let error = rel(p[0][0], monopole);
        eprintln!("T1b ka {ka}: in phase / monopole - 1 = {error:.2e}; anti / in {:.4e}", p[1][0].abs() / p[0][0].abs());
        if ka == 0.05 {
            assert!(error <= 0.01, "T1b: in phase is not the monopole within 1 %: {error:.3e}");
        }
        ratio.push(p[1][0].abs() / p[0][0].abs());
    }
    let slope = (ratio[1] / ratio[0]).ln() / 2.0_f64.ln();
    eprintln!("T1b anti/in slope in ka {slope:.4}");
    assert!((slope - 1.0).abs() <= 0.02, "T1b dipole slope {slope}");
}

// ------------------------------------------------------------------ T1c: card drums

/// closedbody.py's J0 zeros (rounded as there).
const J0Z: [f64; 3] = [2.404826, 5.520078, 8.653728];

/// One lumped-model line: frequency, batter and resonant J0 amplitudes, the
/// numpy change (closed body with image body minus baffled, dB, 17 mm
/// panels), and whether it is gated (a primary prediction line, or above
/// -20 dB re the strongest closed-body line in numpy).
struct Line(f64, [f64; 3], [f64; 3], f64, bool);

fn t1c_lines(piece: &str) -> Vec<Line> {
    match piece {
        "kick" => vec![
            Line(40.6244, [1.971017711e+01, -7.925604209e-01, 2.355633316e-01], [-1.988943353e+01, 7.807258612e-01, -2.318225766e-01], -21.1, true),
            Line(98.5820, [1.881598129e+01, 4.093933572e+01, -3.767698366e+00], [2.547125490e+01, 5.851319122e+01, -5.188722221e+00], 0.1, true),
            Line(106.7175, [6.306774016e+00, 6.331128022e+01, -1.689093062e+00], [-5.157461503e+00, -6.397871966e+01, 1.406337808e+00], -15.4, false),
            Line(149.7653, [-1.212028355e+01, 1.555519539e+01, 1.843470044e+01], [-1.528400380e+01, 1.942451935e+01, 2.423530327e+01], 0.2, true),
            Line(169.6730, [-3.526941873e+00, 3.700308480e+00, 7.277023144e+01], [2.520241826e+00, -2.628450561e+00, -8.821204113e+01], -8.2, false),
            Line(179.0738, [1.009709504e+01, -9.959435421e+00, 6.159338948e+01], [1.233940211e+01, -1.211126621e+01, 6.799094485e+01], 0.1, true),
        ],
        "snare-off" => vec![
            Line(135.4236, [5.482940355e+01, 1.752039365e+00, -4.854883428e-01], [-4.193237628e+01, 4.261909605e+00, -1.350160177e+00], -7.6, true),
            Line(272.8239, [1.254191484e+01, 7.729203422e+01, -3.172643999e+00], [2.047788026e+01, -2.165003656e-01, 6.248457404e-02], -7.3, true),
            Line(374.6525, [1.339665222e+01, -1.869291096e+01, -1.528786857e+01], [3.865090175e+01, 5.525741449e+00, -1.386773578e+00], -0.8, true),
            Line(446.8013, [-1.326583234e+00, 1.367580701e+00, -1.251190403e+02], [-1.661307658e+01, -5.963314817e+00, 1.264554587e+00], -9.1, true),
            Line(626.9486, [5.308933907e-02, -4.268594109e-02, 5.480347518e-02], [-5.004118961e-01, -3.695256005e+00, 1.296685176e-01], -2.3, false),
            Line(1001.7392, [2.629594090e-03, -1.871855417e-03, 1.715739576e-03], [-8.291725811e-02, 8.742013405e-02, 1.280226076e+00], -2.5, false),
        ],
        "tom-12" => vec![
            Line(135.5512, [5.369821423e+01, -1.998414579e+00, 5.921805975e-01], [-5.476034522e+01, 2.037942486e+00, -6.038937125e-01], -7.4, true),
            Line(277.3614, [5.405754796e+01, 3.567158297e+01, -6.220611007e+00], [6.655767609e+01, 4.392018792e+01, -7.659049071e+00], -3.2, true),
            Line(352.1059, [5.135862896e+00, 5.996927114e+01, -1.393757023e+00], [-5.865642461e+00, -6.849059453e+01, 1.591802690e+00], -1.1, true),
            Line(377.2355, [-2.773611816e+00, 1.712478511e+01, 9.842055842e-01], [-3.183177026e+00, 1.965351540e+01, 1.129538238e+00], -4.0, false),
            Line(559.4727, [2.415974974e+00, -2.525629664e+00, -6.629729923e+01], [-3.377648073e+00, 3.530950553e+00, 9.268678170e+01], 3.5, true),
            Line(563.4463, [6.610489016e-01, -6.850738565e-01, -9.156551662e+01], [6.291366831e-01, -6.520018303e-01, -8.714517984e+01], -5.7, true),
        ],
        "tom-14" => vec![
            Line(60.6833, [3.260430725e+01, -1.177425879e+00, 3.485153308e-01], [-3.327365213e+01, 1.179122334e+00, -3.487814191e-01], -10.6, true),
            Line(141.4693, [1.767121962e+01, 3.121847250e+01, -3.292785626e+00], [2.989201709e+01, 5.448599575e+01, -5.634023725e+00], -0.8, true),
            Line(156.6982, [3.304837136e+00, 3.746995209e+01, -8.946838088e-01], [-3.172514461e+00, -4.199817928e+01, 8.692913438e-01], -3.9, true),
            Line(199.4116, [7.549922300e+00, -1.245386592e+01, -5.978236750e+00], [1.143859366e+01, -1.866741131e+01, -9.224259807e+00], 0.0, true),
            Line(248.7513, [1.606861697e+00, -1.683094967e+00, -3.735904996e+01], [-1.742126265e+00, 1.817531985e+00, 5.867823297e+01], -1.1, true),
            Line(255.0379, [-2.224698741e+00, 2.261683958e+00, -3.815987950e+01], [-3.142317893e+00, 3.182940520e+00, -4.434045936e+01], -13.9, false),
        ],
        _ => unreachable!(),
    }
}

fn j0_sum(coef: &[f64; 3], x: f64, y: f64) -> f64 {
    let r = (x * x + y * y).sqrt();
    if r >= 1.0 {
        return 0.0;
    }
    coef.iter().zip(J0Z).map(|(c, z)| c * bessel_j(0, z * r)).sum()
}

#[test]
fn t1c_card_drums_match_the_numpy_route_per_line_and_n2_misses() {
    let mut flipped_miss = 0.0_f64;
    let mut dropped_miss = 0.0_f64;
    let mut failures = Vec::new();
    for piece in ["kick", "snare-off", "tom-12", "tom-14"] {
        let card = card(piece).unwrap();
        let two = scene(&card, true);
        let one = scene(&card, false);
        let s2 = SpherePanels::from_triangles(two.tris.clone()).unwrap();
        let s1 = SpherePanels::from_triangles(one.tris.clone()).unwrap();
        let image_mic = mirror_point(card.mic, card.floor.0, card.floor.1);
        let ids = vec![
            SourceId::StructuralMode("line head motion".into()),
            SourceId::StructuralMode("line head motion, resonant head sign-flipped (N2)".into()),
        ];
        for Line(f, ab, ar, expected, gated) in t1c_lines(piece) {
            let omega = TAU * f;
            let motion = |mesh: &Mesh, flip: f64| {
                trace(mesh, card.radius, &|part, x, y| match part {
                    Part::Batter => j0_sum(&ab, x, y),
                    Part::Resonant => flip * j0_sum(&ar, x, y),
                    Part::Shell => 0.0,
                })
            };
            let closed = forward(&s2, GAS, &ids, &[motion(&two, 1.0), motion(&two, -1.0)], &[card.mic], omega).unwrap();
            let single = forward(&s1, GAS, &ids[..1], &[motion(&one, 1.0)], &[card.mic, image_mic], omega).unwrap();
            let now = baffled(&one, &s1, GAS, &motion(&one, 1.0), &[card.mic, image_mic], omega);
            let now_flipped = baffled(&one, &s1, GAS, &motion(&one, -1.0), &[card.mic, image_mic], omega);
            let change = db(closed[0][0]) - db(now);
            let flipped = db(closed[1][0]) - db(now_flipped);
            let dropped = db(single[0][0] + single[0][1]) - db(now);
            eprintln!(
                "T1c {piece} {f:.1} Hz: change {change:+.2} dB (numpy {expected:+.1}){}; N2 flipped {flipped:+.2}, image body dropped {dropped:+.2}",
                if gated { "" } else { ", not gated" }
            );
            if gated {
                if (change - expected).abs() > 1.0 {
                    failures.push(format!("{piece} {f:.1} Hz: {change:+.2} against {expected:+.1}"));
                }
                flipped_miss = flipped_miss.max((flipped - expected).abs());
                dropped_miss = dropped_miss.max((dropped - expected).abs());
            }
        }
    }
    assert!(failures.is_empty(), "T1c lines beyond 1 dB: {failures:?}");
    eprintln!("N2 largest miss on gated lines: flipped {flipped_miss:.2} dB, image body dropped {dropped_miss:.2} dB");
    assert!(flipped_miss > 1.0, "N2: a sign-flipped resonant head must miss T1c");
    assert!(dropped_miss > 1.0, "N2: dropping the image body must miss T1c");
}

// ------------------------------------------------------------------ mesh admission

#[test]
fn drum_scenes_are_closed_disjoint_bodies_with_the_kit_band_tops() {
    for (piece, f_expected) in [("kick", 532.0), ("snare", 1067.0), ("tom-12", 1018.0), ("tom-14", 687.0)] {
        let card = card(piece).unwrap();
        let mesh = scene(&card, true);
        let surface = SpherePanels::from_triangles(mesh.tris.clone()).unwrap();
        let policy = GeometryPolicy::new(&surface).unwrap();
        assert_eq!(policy.first_dirichlet_bounds().len(), 2, "{piece}: the drum and its image are two bodies");
        let f_x = policy.plain_cbie_limit() * GAS.1 / TAU;
        eprintln!("{piece}: {} panels, f_x {f_x:.1} Hz", mesh.tris.len());
        assert!((f_x - f_expected).abs() < 2.0, "{piece}: f_x {f_x}");
        // Outward closure: the divergence theorem gives the enclosed volume.
        let volume: f64 = surface
            .centroids()
            .iter()
            .zip(surface.normals())
            .zip(surface.areas())
            .map(|((c, n), a)| dot(*c, *n) * a / 3.0)
            .sum();
        let exact = 2.0 * PI * card.radius * card.radius * card.depth;
        assert!((volume / exact - 1.0).abs() < 0.01, "{piece}: volume {volume} against {exact}");
    }
}

/// FIR sizing. The first bakes sized the FIR on the complement alone, and all
/// five refused NotConverged at their last rung, tail ~ wrap: pre-ring past
/// the latency. At the snare-off's band (96 kHz, pass 0.6 f_x) an analytic
/// (i omega)^6 row reproduces it there (N), and the steep sizing holds rows
/// up to (i omega)^8 at delays it was not sized on.
#[test]
fn fir_sizing_holds_steep_rows_where_the_complement_alone_refuses() {
    let (fs, stop) = (96_000.0, 1066.0);
    let pass = 0.6 * stop;
    let spread = 700.0;
    let steep = |order: usize, delay: f64| {
        move |omega: f64| -> Result<Vec<C64>, HelmholtzError> {
            let (s, c) = det::sin_cos(omega * delay / fs);
            let mut z = C64::new(1.0, 0.0);
            for _ in 0..order {
                z = z * C64::new(0.0, omega / (TAU * stop));
            }
            Ok(vec![z * C64::new(c, s), C64::ZERO, C64::ZERO, C64::ZERO])
        }
    };
    let (n0, l0) = size_fir(fs, pass, stop, spread, 0).expect("complement-only sizing");
    let refused = bake_observation_transfers(steep(6, 0.0), &probe_declaration(fs, pass, stop, n0, l0, 1), || false);
    eprintln!("complement-only sizing: FFT {n0}, latency {l0}; (i omega)^6 at once: {:?}", refused.as_ref().err());
    assert!(
        matches!(refused, Err(TransferError::NotConverged { .. })),
        "N: the complement-only sizing must refuse a steep row"
    );
    let (n8, l8) = size_fir(fs, pass, stop, spread, STEEP_PROBE_ORDER).expect("steep sizing");
    eprintln!("steep sizing (p <= {STEEP_PROBE_ORDER}): FFT {n8}, latency {l8}");
    for (order, delay) in [(6, 0.0), (7, 0.5 * spread), (8, 0.25 * spread)] {
        let held = bake_observation_transfers(steep(order, delay), &probe_declaration(fs, pass, stop, n8, l8, 1), || false);
        let last = held
            .as_ref()
            .map(|b| b.ladder.last().map(|r| (r.tail_energy_fraction, r.wrap_energy_fraction, r.dtft_relative_change)));
        eprintln!("(i omega)^{order} behind {delay:.0} samples: {last:?}");
        assert!(held.is_ok(), "(i omega)^{order} behind {delay:.0} samples must hold at the steep sizing");
    }
}
