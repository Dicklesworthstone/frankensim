use super::*;

fn crossings_hz(x: &[f64], rate: f64) -> f64 {
    let mean = x.iter().sum::<f64>() / x.len() as f64;
    let mut idx = vec![];
    for i in 1..x.len() {
        let (a, b) = (x[i - 1] - mean, x[i] - mean);
        if a < 0.0 && b >= 0.0 {
            idx.push(i as f64 - 1.0 + a / (a - b));
        }
    }
    assert!(idx.len() >= 3, "too few cycles");
    (idx.len() - 1) as f64 / ((idx[idx.len() - 1] - idx[0]) / rate)
}

const NO_LOSS: [(f64, f64); 2] = [(100.0, 1e12), (1000.0, 1e12)];

#[test]
fn membrane_fundamental_matches_bessel_zero() {
    let (r, t, th) = (0.15, 2000.0, 2.54e-4);
    let thick = |_r: f64| th;
    let z = |_r: f64| 0.0;
    let mut b = Body::plate("m", &PlateSpec {
        radius: r, hole: 0.0, grid: 2.0e-3, edge: Edge::Clamped, thickness: &thick, rest_height: &z,
        young: 4e9, poisson: 0.38, density: 1390.0, tension: t, bending: false, nonlinear: Nonlinear::None,
        t60: NO_LOSS, support: None, muffle: None,
    }).unwrap();
    for i in 0..b.n {
        let p = b.xy[i];
        let rr = (p[0] * p[0] + p[1] * p[1]).sqrt() / r;
        // close to J0(2.405 r/R): mostly the (0,1) mode
        let v = 1e-4 * (1.0 - rr * rr) * (1.0 - 0.3 * rr * rr);
        b.w[i] = v;
        b.wp[i] = v;
    }
    b.sync();
    let mut ins = Instrument::new(96_000.0);
    let c = ins.add(b);
    let centre = ins.bodies[c].node_at(0.0, 0.0).unwrap();
    let mut out = [];
    let mut trace = vec![];
    for _ in 0..9600 {
        ins.step(&mut out);
        trace.push(ins.bodies[c].w[centre]);
    }
    let f = crossings_hz(&trace, 96_000.0);
    let exact = 2.404826 * (t / (1390.0 * th)).sqrt() / (2.0 * PI * r);
    assert!((f / exact - 1.0).abs() < 0.015, "f01 {f:.2} Hz vs {exact:.2} Hz");
}

#[test]
fn free_plate_lowest_elastic_mode_matches_disc_theory() {
    // free circular plate, nodal diameters 2, no nodal circle: lambda^2 = 5.253 (nu = 0.33)
    let (r, th, e, nu, rho) = (0.2, 1.0e-3, 1.0e11, 0.33, 8600.0);
    let thick = |_r: f64| th;
    let z = |_r: f64| 0.0;
    let mut b = Body::plate("p", &PlateSpec {
        radius: r, hole: 0.0, grid: 8.0e-3, edge: Edge::Free, thickness: &thick, rest_height: &z,
        young: e, poisson: nu, density: rho, tension: 0.0, bending: true, nonlinear: Nonlinear::None,
        t60: NO_LOSS, support: None, muffle: None,
    }).unwrap();
    for i in 0..b.n {
        let p = b.xy[i];
        let v = 1e-5 * (p[0] * p[0] - p[1] * p[1]) / (r * r); // r^2 cos(2 theta)
        b.w[i] = v;
        b.wp[i] = v;
    }
    let probe = b.node_at(0.15, 0.0).unwrap();
    b.sync();
    let mut ins = Instrument::new(96_000.0);
    let c = ins.add(b);
    let (ratio, _) = ins.stability();
    assert!(ratio < 1.0, "unstable: {ratio}");
    let mut out = [];
    let mut trace = vec![];
    for _ in 0..19200 {
        ins.step(&mut out);
        trace.push(ins.bodies[c].w[probe]);
    }
    let f = crossings_hz(&trace, 96_000.0);
    let d = e * th.powi(3) / (12.0 * (1.0 - nu * nu));
    let exact = 5.253 / (2.0 * PI * r * r) * (d / (rho * th)).sqrt();
    assert!((f / exact - 1.0).abs() < 0.04, "(2,0) {f:.2} Hz vs {exact:.2} Hz");
}

#[test]
fn shallow_shell_struck_hard_stays_bounded_and_spreads_energy_upward() {
    let spec = kit::crash();
    let mut p = kit::cymbal(&spec, &[(0.0, 8.0, 0.95, 0.008)]).unwrap();
    let out = kit::render(&mut p, 0.12, |_| {});
    let x = &out[0];
    assert!(x.iter().all(|v| v.is_finite()));
    let peak = x.iter().fold(0.0f64, |a, &b| a.max(b.abs()));
    assert!(peak > 1e-3 && peak < 1e3, "peak {peak}");
}

#[test]
fn a_new_stroke_retires_the_sticks_earlier_contact() {
    // a hi-hat adds one stick contact per stroke; the stale one must stop pushing on the shell
    let mut p = kit::hihat(0.002, 0.0, &[(0.002, 2.0), (0.012, 2.0)]).unwrap();
    let first: Vec<usize> = (0..p.instrument.contacts.len()).filter(|&i| p.instrument.contacts[i].b.iter().any(|&(b, _, _)| p.instrument.bodies[b].name == "stick")).collect();
    assert_eq!(first.len(), 2, "one stick contact per stroke");
    let _ = kit::render(&mut p, 0.02, |_| {});
    let c = &p.instrument.contacts;
    assert!(!c[first[0]].active, "the first stroke's contact is still active after the second launch");
    assert!(c[first[1]].active && c[first[1]].peak_force > 0.0, "the second stroke never landed");
}

/// The auxiliary variable must keep tracking sqrt(2V): a shallow shell carries its curvature
/// stiffness in V, so psi / sqrt(2V) scales that stiffness. Without the relaxation this ratio
/// reached 1.48 and then 0.52 within 2 s on this live-grid ride (and the 8-16 kHz bands grew back).
#[test]
fn shell_auxiliary_variable_tracks_its_potential() {
    let mut p = live::Live::new("ride", 2.5).unwrap();
    p.strike(0.6 * p.radius, 0.0, 2.5);
    let mut blk = vec![0.0f32; 4800];
    let mut worst = 0.0f64;
    for _ in 0..20 {
        p.render(&mut blk, 1.0);
        let (psi, r) = p.ins.bodies[0].sav_state();
        worst = worst.max((psi / r - 1.0).abs());
    }
    assert!(worst < 0.08, "psi / sqrt(2V) drifted {:.1} % from 1 within 2 s", 100.0 * worst);
}
