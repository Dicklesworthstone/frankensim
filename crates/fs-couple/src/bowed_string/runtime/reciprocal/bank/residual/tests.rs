//! Residual flexibility against closed forms and against the exact
//! point-supported string (MagentaBadger 7592, gate 1). The exact roots come
//! from shooting the piecewise-sine solution here, with std trigonometry; they
//! share nothing with the residual owner. Predictions: sharp/pred-r3.md.

use super::*;

const SPAN_M: f64 = 0.383;
const A_TENSION_N: f64 = 55.486_491_050_429_706;
const DENSITY_KG_M: f64 = 1e-3;
const BRIDGE: f64 = 0.856_396_866_840_731_1;
const BRIDGE_N_M: f64 = 2e5;

fn card(bending: f64, modes: usize) -> BowedStringCard {
    BowedStringCard { length_m: SPAN_M, tension_n: A_TENSION_N, linear_density_kg_m: DENSITY_KG_M,
        bending_stiffness_n_m2: bending, viscous_bending_n_m2_s: 0.0, mode_count: modes,
        zetas: vec![0.0; modes], sample_rate_hz: 48_000 }
}

#[test]
fn static_green_is_the_complete_modal_sum() {
    // With EI the omitted tail falls as 1/N^3: 40 000 sines leave about 2e-10 of it.
    let beam = card(2e-5, 40_000);
    for (x, y) in [(0.1, 0.1), (0.1, BRIDGE), (0.5, 0.52), (BRIDGE, BRIDGE), (0.93, 0.07)] {
        let exact = static_green_m_n(&beam, x, y);
        let sum = retained_compliance_m_n(&beam, x, y);
        assert!((exact-sum).abs() <= 1e-8*exact.abs(), "x {x} y {y}: {exact:e} vs {sum:e}");
        assert_eq!(exact.to_bits(), static_green_m_n(&beam, y, x).to_bits());
        assert_eq!(sum.to_bits(), retained_compliance_m_n(&beam, y, x).to_bits());
    }
    // The taut string's own diagonal tail is 4/(pi^2 N) of it.
    let string = card(0.0, 1_000_000);
    let exact = static_green_m_n(&string, 0.3, 0.3);
    let textbook = 0.3*0.7*SPAN_M/A_TENSION_N;
    assert!((exact-textbook).abs() <= 4e-15*textbook, "{exact:e} vs {textbook:e}");
    let sum = retained_compliance_m_n(&string, 0.3, 0.3);
    assert!(exact-sum > 0.0 && exact-sum <= 1e-6*exact, "{exact:e} vs {sum:e}");
}

// Independent sine eigenfunction sum, using std trigonometry and summing the
// smallest terms first. Starting above the retained modes checks their omitted tail
// without subtracting the retained compliance from the production Green function.
fn modal_compliance_reference(card: &BowedStringCard, x: f64, y: f64, first: usize) -> f64 {
    (first..=40_000).rev().map(|n| {
        let order = n as f64*core::f64::consts::PI;
        let wave = order/card.length_m;
        2.0*(order*x).sin()*(order*y).sin()
            /(card.length_m*(card.tension_n*wave.powi(2)+card.bending_stiffness_n_m2*wave.powi(4)))
    }).sum()
}

#[test]
fn bending_dominated_green_matches_modes_and_is_continuous_at_the_series_switch() {
    let mut beam = card(1e-3, 8);
    beam.length_m = 0.5;
    for (x, y) in [(0.1, 0.1), (0.1, BRIDGE), (0.5, 0.52), (0.8, 0.8), (0.93, 0.07)] {
        let mut below_switch = 0.0;
        for q in [1e-20, 1e-12, 1e-8, 1e-3*(1.0-1e-8), 1e-3*(1.0+1e-8)] {
            beam.tension_n = q*beam.bending_stiffness_n_m2/beam.length_m.powi(2);
            let actual = static_green_m_n(&beam, x, y);
            let reference = modal_compliance_reference(&beam, x, y, 1);
            assert!(actual > 0.0 && (actual-reference).abs() <= 5e-11*reference,
                "q {q:e}, ({x}, {y}): Green {actual:e}, independent modes {reference:e}");
            assert_eq!(actual.to_bits(), static_green_m_n(&beam, y, x).to_bits());
            if q < 1e-3 { below_switch = actual; }
            else {
                assert!((actual-below_switch).abs() <= 5e-11*below_switch,
                    "({x}, {y}): switch jump from {below_switch:e} to {actual:e}");
            }
        }
    }
}

#[test]
fn bending_dominated_green_has_the_pinned_beam_limit_even_near_the_ends() {
    let mut beam = card(1e-3, 8);
    beam.length_m = 0.5;
    beam.tension_n = 1e-18;
    let (l, ei) = (beam.length_m, beam.bending_stiffness_n_m2);
    for (x, y) in [(1e-10_f64, 1e-10), (1e-10, 0.8), (0.1, 0.1), (0.8, 0.8),
        (1.0-1e-10, 1.0-1e-10)] {
        let (a, b) = (x.min(y)*l, x.max(y)*l);
        // The piecewise cubic solution of EI*w''''=delta, w=w''=0 at
        // both ends: a*(L-b)*(L^2-a^2-(L-b)^2)/(6*EI*L).
        let reference = a*(l-b)*(2.0*b*(l-b)+(b-a)*(b+a))/(6.0*ei*l);
        let actual = static_green_m_n(&beam, x, y);
        assert!(actual > 0.0 && (actual-reference).abs() <= 2e-14*reference,
            "({x}, {y}): Green {actual:e}, pinned beam {reference:e}");
    }
}

#[test]
fn bending_dominated_polarized_bridge_is_passive_and_rotation_invariant() {
    let mut beam = card(1e-3, 4);
    beam.length_m = 0.5;
    beam.tension_n = 1e-18;
    let (bridge_station, point_station) = (0.8, 0.3);
    let condensed = condense(residual_compliance_m_n(&beam, bridge_station, bridge_station).unwrap(),
        BRIDGE_N_M, 3.0).unwrap();
    assert!(condensed.stiffness_n_m > 0.0 && condensed.stiffness_n_m < BRIDGE_N_M);
    assert!(condensed.damping_n_s_m > 0.0 && condensed.damping_n_s_m < 3.0);
    let tail = |x, y| modal_compliance_reference(&beam, x, y, beam.mode_count+1);
    let (bb, sb, ss) = (tail(bridge_station, bridge_station), tail(point_station, bridge_station),
        tail(point_station, point_station));
    let reference_k = BRIDGE_N_M/(1.0+BRIDGE_N_M*bb);
    let reference_axis = [0.6, 0.8];
    let reference_cross = reference_axis.map(|axis| reference_k*axis*sb);
    let reference_compliance = [ss-reference_k*(0.6*sb).powi(2), -reference_k*0.6*0.8*sb*sb,
        -reference_k*0.6*0.8*sb*sb, ss-reference_k*(0.8*sb).powi(2)];
    let mut baseline: Option<(Vec<f64>, Vec<f64>)> = None;
    for angle in [0.0_f64, 0.51, 1.7] {
        let (s, c) = angle.sin_cos();
        let rotate = |v: [f64; 2]| [c*v[0]-s*v[1], s*v[0]+c*v[1]];
        let bridge = ResidualBridge { card: &beam, station_fraction: bridge_station,
            axis: rotate(reference_axis), condensed };
        let points = [[1.0, 0.0], [0.0, 1.0]].map(|direction| ResidualStop {
            string: 0, station_fraction: point_station, direction: rotate(direction) });
        let (cross, g) = stop_residual(&[bridge], &points).unwrap();
        for (actual, reference) in cross.iter().zip(reference_cross) {
            assert!((actual-reference).abs() <= 2e-9*reference.abs());
        }
        for (actual, reference) in g.iter().zip(reference_compliance) {
            assert!((actual-reference).abs() <= 2e-9*ss);
        }
        assert_eq!(g[1].to_bits(), g[2].to_bits());
        assert!(g[0] > 0.0 && g[3] > 0.0 && g[0]*g[3]-g[1]*g[2] > 0.0);
        if let Some((base_cross, base_g)) = &baseline {
            for (actual, original) in cross.iter().zip(base_cross) {
                assert!((actual-original).abs() <= 2e-14*reference_cross[1].abs());
            }
            for (actual, original) in g.iter().zip(base_g) {
                assert!((actual-original).abs() <= 2e-14*ss);
            }
        } else { baseline = Some((cross, g)); }
    }
}

#[test]
fn residual_is_positive_symmetric_and_refuses_supports_off_the_span() {
    let a = card(1.6e-5, 25);
    let stations = [0.09, 0.2, BRIDGE];
    for &x in &stations {
        assert!(residual_compliance_m_n(&a, x, x).unwrap() > 0.0);
        for &y in &stations {
            assert_eq!(residual_compliance_m_n(&a, x, y).unwrap().to_bits(),
                residual_compliance_m_n(&a, y, x).unwrap().to_bits());
        }
    }
    for (x, y) in [(0.0, 0.5), (0.5, 1.0), (-0.1, 0.5), (0.5, f64::NAN)] {
        assert!(residual_compliance_m_n(&a, x, y).is_err());
    }
    let mut slack = a.clone();
    slack.tension_n = 0.0;
    assert!(residual_compliance_m_n(&slack, 0.5, 0.5).is_err());
    assert!(condense(-1e-9, BRIDGE_N_M, 0.0).is_err());
    // k G = 5: the spring and dashpot see 1/6 and 1/36 of themselves.
    let condensed = condense(2.5e-5, BRIDGE_N_M, 3.0).unwrap();
    assert!((condensed.stiffness_n_m*6.0/BRIDGE_N_M-1.0).abs() <= 1e-15);
    assert!((condensed.damping_n_s_m*36.0/3.0-1.0).abs() <= 1e-15);
}

// y(L) for wavenumber kw: y = sin(kw x) from the nut, each spring K to ground
// adding the slope jump K y / T. Springs are in ascending station order.
fn shoot(kw: f64, springs: &[(f64, f64)]) -> f64 {
    let (mut y, mut dy, mut from) = (0.0_f64, kw, 0.0);
    for &(station, k) in springs {
        let x = station*SPAN_M;
        let (s, c) = (kw*(x-from)).sin_cos();
        (y, dy) = (y*c + dy/kw*s, -y*kw*s + dy*c);
        dy += k*y/A_TENSION_N;
        from = x;
    }
    let (s, c) = (kw*(SPAN_M-from)).sin_cos();
    y*c + dy/kw*s
}

// The lowest `count` exact frequencies [Hz], by a fine wavenumber scan and bisection.
fn exact(springs: &[(f64, f64)], count: usize) -> Vec<f64> {
    let step = 1e-3;
    let mut roots = Vec::new();
    let (mut lo, mut f_lo) = (step, shoot(step, springs));
    while roots.len() < count {
        let hi = lo+step;
        let f_hi = shoot(hi, springs);
        if f_lo == 0.0 || f_lo.signum() != f_hi.signum() {
            let (mut a, mut b, f_a) = (lo, hi, f_lo);
            for _ in 0..200 {
                let mid = 0.5*(a+b);
                if mid <= a || mid >= b { break; }
                if shoot(mid, springs).signum() == f_a.signum() { a = mid; } else { b = mid; }
            }
            roots.push(0.5*(a+b));
        }
        (lo, f_lo) = (hi, f_hi);
        assert!(lo < 5e3, "scan found only {} roots", roots.len());
    }
    roots.iter().map(|kw| kw*(A_TENSION_N/DENSITY_KG_M).sqrt()/(2.0*core::f64::consts::PI)).collect()
}

// Mass-normalized sine shape of `card` at `station`.
fn shape(card: &BowedStringCard, station: f64) -> Vec<f64> {
    let scale = (2.0/(card.linear_density_kg_m*card.length_m)).sqrt();
    (1..=card.mode_count).map(|n| scale*(n as f64*core::f64::consts::PI*station).sin()).collect()
}

// Ascending frequencies [Hz] of the retained sines plus the supports, each a
// rank-one stiffness k c c^T.
fn retained(card: &BowedStringCard, supports: &[(f64, Vec<f64>)]) -> Vec<f64> {
    let n = card.mode_count;
    let mut stiffness = vec![0.0; n*n];
    for m in 0..n {
        let k = (m+1) as f64*core::f64::consts::PI/card.length_m;
        stiffness[m*n+m] = (card.tension_n*k*k + card.bending_stiffness_n_m2*k*k*k*k)/card.linear_density_kg_m;
    }
    for (k, column) in supports {
        for i in 0..n {
            for j in 0..n { stiffness[i*n+j] += k*column[i]*column[j]; }
        }
    }
    let (squared, _) = fs_la::eigen::jacobi_eigh(&stiffness, n);
    squared.iter().map(|w2| w2.sqrt()/(2.0*core::f64::consts::PI)).collect()
}

fn cents(f: f64, reference: f64) -> f64 { 1200.0*(f/reference).log2() }

#[test]
fn residual_flexibility_restores_the_exact_point_supported_string() {
    // The taut A string at N 25 with B4's or D5's stop (its summed K at that
    // station, to ground), alone and beside the 2e5 N/m bridge. Gate: f1..f4
    // within 1 c of exact and at most a quarter of the uncorrected error.
    let a = card(0.0, 25);
    for (note, semitones) in [("B4", 2.0_f64), ("D5", 5.0)] {
        let stop = BRIDGE*(1.0-(-semitones/12.0).exp2());
        for bridge_n_m in [0.0, BRIDGE_N_M] {
            for stop_n_m in [1e5, 1e6, 1e7] {
                let mut springs = vec![(stop, stop_n_m)];
                if bridge_n_m > 0.0 { springs.push((BRIDGE, bridge_n_m)); }
                let truth = exact(&springs, 10);
                let (phi_s, phi_b) = (shape(&a, stop), shape(&a, BRIDGE));
                let plain = retained(&a, &[(stop_n_m, phi_s.clone()), (bridge_n_m, phi_b.clone())]);

                let bridge = ResidualBridge { card: &a, station_fraction: BRIDGE, axis: [1.0, 0.0],
                    condensed: condense(residual_compliance_m_n(&a, BRIDGE, BRIDGE).unwrap(), bridge_n_m, 0.0).unwrap() };
                let (cross, compliance) = stop_residual(&[bridge],
                    &[ResidualStop { string: 0, station_fraction: stop, direction: [1.0, 0.0] }]).unwrap();
                let column: Vec<f64> = phi_s.iter().zip(&phi_b).map(|(s, b)| s-cross[0]*b).collect();
                let corrected = retained(&a, &[(stop_n_m/(1.0+stop_n_m*compliance[0]), column),
                    (bridge.condensed.stiffness_n_m, phi_b)]);

                let before: Vec<f64> = plain.iter().zip(&truth).map(|(f, e)| cents(*f, *e)).collect();
                let after: Vec<f64> = corrected.iter().zip(&truth).map(|(f, e)| cents(*f, *e)).collect();
                eprintln!("{note} bridge {bridge_n_m:e} K {stop_n_m:e}: uncorrected f1..4 {:+.3?}; corrected f1..4 \
                    {:+.3?}; corrected n5..10 (ungated) {:+.2?}", &before[..4], &after[..4], &after[4..10]);
                for (n, (after, before)) in after.iter().zip(&before).take(4).enumerate() {
                    assert!(after.abs() <= 1.0 && after.abs() <= 0.25*before.abs(),
                        "{note} bridge {bridge_n_m:e} K {stop_n_m:e} f{}: {after:+.3} c corrected, {before:+.3} c uncorrected",
                        n+1);
                }
            }
        }
    }
}

#[test]
fn bow_residual_enrichment_converges_to_the_exact_driven_string() {
    // Independent dynamic Green function of a taut pinned string, below its
    // first resonance, with the same spring bridge to ground. This exercises
    // the bow/bridge cross-map as well as the bow's own omitted compliance.
    let bow = 0.73;
    let omega = 0.4*core::f64::consts::PI*(A_TENSION_N/DENSITY_KG_M).sqrt()/SPAN_M;
    let wave = omega*(DENSITY_KG_M/A_TENSION_N).sqrt();
    let exact = |x:f64,y:f64| {
        let (a,b)=(x.min(y)*SPAN_M,x.max(y)*SPAN_M);
        (wave*a).sin()*(wave*(SPAN_M-b)).sin()/(A_TENSION_N*wave*(wave*SPAN_M).sin())
    };
    let truth = exact(bow,bow)-BRIDGE_N_M*exact(bow,BRIDGE).powi(2)
        /(1.0+BRIDGE_N_M*exact(BRIDGE,BRIDGE));
    let mut previous = f64::INFINITY;
    let mut previous_tail = f64::INFINITY;
    for count in [4,8,16] {
        let card = card(0.0,count);
        let bridge = ResidualBridge { card:&card,station_fraction:BRIDGE,axis:[1.0,0.0],
            condensed:condense(residual_compliance_m_n(&card,BRIDGE,BRIDGE).unwrap(),BRIDGE_N_M,0.0).unwrap() };
        let (cross,g)=stop_residual(&[bridge],&[ResidualStop {string:0,station_fraction:bow,direction:[1.0,0.0]}]).unwrap();
        let response=|x:f64,y:f64| (1..=count).map(|n| {
            let order=n as f64*core::f64::consts::PI;
            2.0*(order*x).sin()*(order*y).sin()
                /(SPAN_M*(A_TENSION_N*(order/SPAN_M).powi(2)-DENSITY_KG_M*omega*omega))
        }).sum::<f64>();
        let (aa,ab,bb)=(response(bow,bow),response(bow,BRIDGE),response(BRIDGE,BRIDGE));
        let k=bridge.condensed.stiffness_n_m;
        let corrected=aa-2.0*cross[0]*ab+cross[0]*cross[0]*bb
            -k*(ab-cross[0]*bb).powi(2)/(1.0+k*bb)+g[0];
        let plain=aa-BRIDGE_N_M*ab*ab/(1.0+BRIDGE_N_M*bb);
        let error=(corrected-truth).abs();
        assert!(error<0.01*(plain-truth).abs(),"N={count}: corrected {error:e}, plain {:e}",(plain-truth).abs());
        assert!(error<previous && g[0]>0.0 && g[0]<previous_tail);
        previous=error;previous_tail=g[0];
    }
}

#[test]
fn a_stop_and_a_bow_share_one_exact_static_green_function() {
    // MagentaBadger 12160 (a): one point list [stop, bow]. The retained sines through the condensed
    // bridge plus the joint omitted compliance, stop/bow block included, are the exact static Green
    // function with the bridge spring, at each point and between them. A rigid stop then leaves the
    // bow on the exact stopped span; for the taut string that span's own closed form is checked too.
    let (stop, bow) = (0.25, 0.73);
    let x = [stop, bow];
    for (bending, count) in [(0.0, 64), (0.0, 128), (2e-5, 64), (2e-5, 128)] {
        let card = card(bending, count);
        let bridge = ResidualBridge { card:&card, station_fraction:BRIDGE, axis:[1.0,0.0],
            condensed:condense(residual_compliance_m_n(&card,BRIDGE,BRIDGE).unwrap(),BRIDGE_N_M,0.0).unwrap() };
        let points = x.map(|station| ResidualStop { string:0, station_fraction:station, direction:[1.0,0.0] });
        let (cross, g) = stop_residual(&[bridge], &points).unwrap();
        assert_eq!(g[1].to_bits(), g[2].to_bits());
        let k = bridge.condensed.stiffness_n_m;
        let r = |a: f64, b: f64| retained_compliance_m_n(&card, a, b);
        // Each point's retained column is its sine row minus cross times the bridge row; the
        // condensed spring k joins the bridge row to ground (Sherman-Morrison).
        let model = |i: usize, j: usize| {
            let (a, b) = (x[i], x[j]);
            let own = r(a,b) - cross[i]*r(BRIDGE,b) - cross[j]*r(a,BRIDGE) + cross[i]*cross[j]*r(BRIDGE,BRIDGE);
            let (u, v) = (r(a,BRIDGE)-cross[i]*r(BRIDGE,BRIDGE), r(b,BRIDGE)-cross[j]*r(BRIDGE,BRIDGE));
            own - k*u*v/(1.0+k*r(BRIDGE,BRIDGE)) + g[i*2+j]
        };
        let green = |a: f64, b: f64| static_green_m_n(&card, a, b);
        let exact = |i: usize, j: usize| green(x[i],x[j])
            - BRIDGE_N_M*green(x[i],BRIDGE)*green(BRIDGE,x[j])/(1.0+BRIDGE_N_M*green(BRIDGE,BRIDGE));
        for i in 0..2 {
            for j in 0..2 {
                let (m, e) = (model(i, j), exact(i, j));
                assert!((m-e).abs() <= 1e-9*e.abs(), "EI {bending} N {count} ({i},{j}): {m:e} vs {e:e}");
            }
        }
        let stopped = |c: &dyn Fn(usize, usize) -> f64| c(1,1) - c(1,0)*c(0,1)/c(0,0);
        let (m, e) = (stopped(&model), stopped(&exact));
        assert!((m-e).abs() <= 1e-9*e.abs(), "EI {bending} N {count} stopped: {m:e} vs {e:e}");
        // Dropping only the stop/bow block leaves the bow far from the stopped span.
        let open = model(1,1) - (model(1,0)-g[2])*(model(0,1)-g[1])/model(0,0);
        assert!((open-e).abs() > 1e3*(m-e).abs(), "EI {bending} N {count}: without the block {open:e}, exact {e:e}");
        if bending == 0.0 {
            // The taut stopped span [stop, 1] with the same bridge spring, in closed form.
            let span = (1.0-stop)*SPAN_M;
            let taut = |a: f64, b: f64| {
                let (a, b) = ((a.min(b)-stop)*SPAN_M, (a.max(b)-stop)*SPAN_M);
                a*(span-b)/(A_TENSION_N*span)
            };
            let closed = taut(bow,bow) - BRIDGE_N_M*taut(bow,BRIDGE).powi(2)/(1.0+BRIDGE_N_M*taut(BRIDGE,BRIDGE));
            assert!((m-closed).abs() <= 1e-9*closed, "N {count}: stopped {m:e} vs closed form {closed:e}");
        }
    }
}
