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

fn point_contact_plate() -> Body {
    Body::plate("contact", &PlateSpec {
        radius: 0.1, hole: 0.0, grid: 0.02, edge: Edge::Free,
        thickness: &|_| 1e-3, rest_height: &|_| 0.0, young: 1e11,
        poisson: 0.33, density: 8600.0, tension: 0.0, bending: true,
        nonlinear: Nonlinear::None, t60: NO_LOSS, support: None, muffle: None,
    }).unwrap()
}

#[test]
fn point_contact_retains_its_position_and_reconstructs_an_affine_surface() {
    let b = point_contact_plate();
    let (x, y) = (1.35 * b.h, -0.4 * b.h);
    let footprint = b.footprint(x, y, 0.0);
    assert_eq!(footprint.len(), 4);
    assert!(footprint.iter().all(|&(_, w)| w > 0.0 && w < 1.0));
    let sum: f64 = footprint.iter().map(|&(_, w)| w).sum();
    let px: f64 = footprint.iter().map(|&(i, w)| w * b.xy[i][0]).sum();
    let py: f64 = footprint.iter().map(|&(i, w)| w * b.xy[i][1]).sum();
    let surface: f64 = footprint.iter().map(|&(i, w)| w * (2.0 + 3.0 * b.xy[i][0] - 4.0 * b.xy[i][1])).sum();
    assert!((sum - 1.0).abs() < 1e-14);
    assert!((px - x).abs() < 1e-14 && (py - y).abs() < 1e-14);
    assert!((surface - (2.0 + 3.0 * x - 4.0 * y)).abs() < 1e-14);
}

#[test]
fn point_contact_moves_continuously_across_the_nearest_node_boundary() {
    let b = point_contact_plate();
    let centres: Vec<f64> = [0.49, 0.51].iter().map(|x| {
        b.footprint(x * b.h, 0.0, 0.0).iter().map(|&(i, w)| w * b.xy[i][0]).sum()
    }).collect();
    assert!((centres[1] - centres[0] - 0.02 * b.h).abs() < 1e-14);
    let node = b.node_at(b.h, -b.h).unwrap();
    assert_eq!(b.footprint(b.h, -b.h, 0.0), vec![(node, 1.0)]);
}

#[test]
fn interpolated_point_force_and_surface_velocity_do_the_same_work() {
    let b = point_contact_plate();
    let (x, y, force) = (-0.3 * b.h, 1.7 * b.h, 23.0);
    let footprint = b.footprint(x, y, 0.0);
    let nodal_work: f64 = footprint.iter().map(|&(i, w)| {
        force * w * (0.3 + 2.0 * b.xy[i][0] - b.xy[i][1])
    }).sum();
    let point_work = force * (0.3 + 2.0 * x - y);
    assert!((nodal_work - point_work).abs() < 1e-13);
}

#[test]
fn finite_contact_faces_and_incomplete_grid_cells_keep_their_existing_rules() {
    let b = point_contact_plate();
    let face = b.footprint(0.0, 0.0, 1.1 * b.h);
    assert_eq!(face.len(), 5);
    assert!(face.iter().all(|&(_, w)| w == 0.2));
    // The free-edge plate retains the union of complete curvature stencils;
    // its axial grid boundary is one node inside this geometric radius.
    let nearest = b.node_at(4.0 * b.h, 0.0).unwrap();
    assert_eq!(b.footprint(6.0 * b.h, 0.0, 0.0), vec![(nearest, 1.0)]);
}

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

/// `image` (the export for the GPU renderer) only reads: exporting before the first step leaves every
/// output sample unchanged, and the exported Airy solve is symmetric positive definite, as A^-1 must be.
#[test]
fn exporting_an_instrument_changes_nothing_and_exports_a_spd_airy_inverse() {
    let run = |export: bool| {
        let mut p = live::Live::new("ride", 2.5).unwrap();
        p.strike(0.6 * p.radius, 0.0, 2.5);
        if export {
            let img = p.ins.image();
            assert_eq!(img.bodies.len(), p.ins.bodies.len());
            assert!(img.bodies.iter().all(|b| b.cf.len() == b.n && b.lap.ptr.len() == b.n + 1));
            assert!(img.mics.iter().all(|m| m.ring_len.is_power_of_two() && !m.taps.is_empty()));
        }
        let mut blk = vec![0.0f32; 960];
        p.render(&mut blk, 1.0);
        (blk.iter().map(|x| x.to_bits()).collect::<Vec<_>>(), p)
    };
    let (plain, _) = run(false);
    let (exported, mut p) = run(true);
    assert_eq!(plain, exported, "exporting changed the output");
    let img = p.ins.image();
    let (bi, nc) = img.bodies.iter().enumerate().find_map(|(i, b)| b.vk.as_ref().map(|v| (i, v.nb.len()))).unwrap();
    assert!(nc > 10, "no von Karman core found");
    let shell = &p.ins.bodies[bi];
    assert!(!shell.airy_solve(&mut vec![0.0; nc + 1]), "a wrong-length solve must refuse");
    let col = |j: usize| {
        let mut x = vec![0.0; nc];
        x[j] = 1.0;
        assert!(shell.airy_solve(&mut x));
        x
    };
    let (a, b) = (col(3), col(nc / 2));
    assert!(a[3] > 0.0 && b[nc / 2] > 0.0, "diagonal of A^-1 not positive");
    let rel = (a[nc / 2] - b[3]).abs() / a[nc / 2].abs().max(b[3].abs());
    assert!(rel < 1e-9, "A^-1 not symmetric: {rel:e}");
    // positive definite along a few rough vectors too (x' A^-1 x > 0)
    for seed in 1..4u64 {
        let x: Vec<f64> = (0..nc as u64).map(|i| (((i + 1) * 2_654_435_761 * seed) % 1000) as f64 / 500.0 - 1.0).collect();
        let mut y = x.clone();
        assert!(shell.airy_solve(&mut y));
        let q: f64 = x.iter().zip(&y).map(|(a, b)| a * b).sum();
        assert!(q > 0.0, "x' A^-1 x = {q:e} for seed {seed}");
    }
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

#[test]
fn stencil_runs_sum_every_row_exactly_as_its_csr_row() {
    // the grid interior is computed as stencils; every row must come out bit for bit as the CSR row
    // (same products, same ascending-column order), signed zeros included
    let mut covered = 0;
    for (name, scale) in [("kick", 2.5), ("snare", 7.0), ("hihat-closed", 5.0), ("crash", 4.0), ("ride", 3.0)] {
        let p = live::Live::new(name, scale).unwrap();
        for b in &p.ins.bodies {
            for m in [&b.lap, &b.stiff] {
                covered += m.runs.iter().map(|r| r.len).sum::<usize>();
                let val = |i: usize, salt: usize| match (i * 7919 + salt) % 13 {
                    0 => 0.0,
                    1 => -0.0,
                    k => (k as f64 - 6.5) * 1.37e-3 / (1.0 + (i % 5) as f64),
                };
                let x: Vec<f64> = (0..b.n).map(|i| val(i, 3)).collect();
                let y0: Vec<f64> = (0..b.n).map(|i| val(i, 11)).collect();
                for s in [1.0, -1.0, 0.37] {
                    let (mut ya, mut yb) = (y0.clone(), y0.clone());
                    m.mul_add(&x, s, &mut ya);
                    m.rows_add(0, m.n, &x, s, &mut yb);
                    for i in 0..b.n {
                        assert_eq!(ya[i].to_bits(), yb[i].to_bits(), "{name}@{scale} {} row {i}", b.name);
                    }
                }
            }
        }
    }
    assert!(covered > 1000, "the stencil runs cover only {covered} rows");
}

#[test]
fn queued_contact_solves_match_solving_one_at_a_time() {
    // the reference: every active contact in list order, each seeing the earlier ones' updates
    fn one_at_a_time(bodies: &mut [Body], contacts: &mut [Contact], k: f64) {
        for c in contacts.iter_mut().filter(|c| c.active) {
            let (mut d_pred, mut comp) = (-c.gap, 0.0);
            for &(b, i, w) in &c.a {
                d_pred += w * bodies[b].wn[i];
                comp += w * w * bodies[b].compliance(i, k);
            }
            for &(b, i, w) in &c.b {
                d_pred -= w * bodies[b].wn[i];
                comp += w * w * bodies[b].compliance(i, k);
            }
            let f = contact_force(d_pred, comp, c.force, c.d_now, c.law, c.fixed_point_exit, k);
            apply_force(bodies, c, f, d_pred, comp, k);
        }
    }
    // a snare whose wire list also holds contacts that share nodes with others (sqrt-law copies, which
    // must close the batch, and power-law copies, which are solved alone), struck and ringing; the
    // copies leave the fixed-point exit off, so a batch mixes lanes with and without it
    let build = || {
        let mut p = live::Live::new("snare", 7.0).unwrap();
        let n = p.ins.contacts.len();
        for (at, from, alpha) in [(3, 9, 1.5), (11, 4, 2.0), (n - 1, 20, 1.5), (7, 7, 2.5)] {
            let c = &p.ins.contacts[from];
            let Law::Power { k, mu, .. } = c.law;
            let copy = Contact::new(c.a.clone(), c.b.clone(), c.gap + 1e-7, Law::Power { k, alpha, mu });
            p.ins.contacts.insert(at, copy);
        }
        p.strike(0.06, 0.02, 4.5);
        p
    };
    let (mut a, mut b) = (build(), build());
    let k = a.ins.k();
    let mut out = [0.0f32; 64];
    let (mut solved, mut pressed) = (0, 0);
    for round in 0..120 {
        a.render(&mut out, 1.0);
        b.render(&mut out, 1.0);
        // the same predicted positions in both: free flight plus a small push, some apart, some pressed
        for (ins, solve) in [(&mut a.ins, true), (&mut b.ins, false)] {
            for (bi, body) in ins.bodies.iter_mut().enumerate() {
                for i in 0..body.n {
                    let salt = ((i * 2654435761 + bi * 97 + round * 31) % 1009) as f64 / 1009.0 - 0.45;
                    body.wn[i] = 2.0 * body.w[i] - body.wp[i] + salt * 4e-7;
                }
            }
            if solve {
                solve_contacts(&mut ins.bodies, &mut ins.contacts, &mut ins.cq, k);
            } else {
                one_at_a_time(&mut ins.bodies, &mut ins.contacts, k);
            }
        }
        for (ca, cb) in a.ins.contacts.iter().zip(&b.ins.contacts) {
            assert_eq!(ca.force.to_bits(), cb.force.to_bits(), "round {round}: force");
            assert_eq!(ca.d_now.to_bits(), cb.d_now.to_bits(), "round {round}: d_now");
            assert_eq!(ca.peak_force.to_bits(), cb.peak_force.to_bits(), "round {round}: peak");
            solved += 1;
            pressed += usize::from(ca.force > 0.0);
        }
        for (ba, bb) in a.ins.bodies.iter().zip(&b.ins.bodies) {
            assert!(ba.wn.iter().zip(&bb.wn).all(|(x, y)| x.to_bits() == y.to_bits()), "round {round}: {} moved differently", ba.name);
        }
    }
    assert!(pressed > solved / 4 && pressed < solved, "only {pressed} of {solved} solves pressed");
}

#[test]
fn empty_matrix_and_nan_penetration_do_nothing() {
    // an empty matrix adds nothing (a default Csr has no ptr[0] for the CSR rows to slice)
    let mut y = [0.25];
    Csr::default().mul_add(&[], 1.0, &mut []);
    Csr::default().mul_add(&[1.0], 2.0, &mut y);
    assert_eq!(y, [0.25]);
    // a NaN penetration exerts no force, solved alone or queued
    let law = Law::Power { k: 3.0e7, alpha: 1.5, mu: 0.35 };
    assert_eq!(contact_force(f64::NAN, 1e-6, 0.0, 0.0, law, true, 1e-5).to_bits(), 0.0f64.to_bits());
    let mut bodies = vec![Body::point("a", 0.01), Body::point("b", 0.01)];
    let mut contacts = vec![Contact::new(vec![(0, 0, 1.0)], vec![(1, 0, 1.0)], 0.0, law)];
    contacts[0].force = 1.0;
    bodies[0].wn[0] = f64::NAN;
    solve_contacts(&mut bodies, &mut contacts, &mut ContactQueue::default(), 1e-5);
    assert_eq!(contacts[0].force.to_bits(), 0.0f64.to_bits());
}

#[test]
fn rim_friction_sticks_inside_its_cone_slides_on_it_and_only_removes_energy() {
    // the closed hat's rim frictions on a coarse grid, each solved alone on a random step
    let rim = Law::Power { k: 2.0e9, alpha: 1.5, mu: 0.08 };
    let mut ins = kit::hihat_grid_mu(0.0, 80.0, &[], rim, 0.8, 6.0, 0.4).unwrap().instrument;
    let k = ins.k();
    for b in ins.bodies.iter_mut() {
        b.coefs(k);
    }
    assert!(ins.frictions.len() >= 8, "only {} rim frictions", ins.frictions.len());
    let mut seed = 7u64;
    let mut rnd = move || {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (seed >> 11) as f64 / (1u64 << 53) as f64 - 0.5
    };
    let slip = |bodies: &[Body], f: &Friction| {
        let s = |t: &[(usize, usize, f64)]| t.iter().map(|&(b, i, w)| w * (bodies[b].wn[i] - bodies[b].w[i])).sum::<f64>();
        (s(&f.sx), s(&f.sy))
    };
    let (mut stuck, mut slid) = (0, 0);
    for (normal, expect_stick) in [(1.0e7, true), (1.0e-6, false)] {
        for j in 0..ins.frictions.len() {
            for b in ins.bodies.iter_mut() {
                for i in 0..b.n {
                    b.w[i] = 0.0;
                    b.wn[i] = 1.0e-6 * rnd();
                }
            }
            let c = ins.frictions[j].contact;
            ins.contacts[c].active = true;
            ins.contacts[c].force = normal;
            let (bx, by) = slip(&ins.bodies, &ins.frictions[j]);
            solve_frictions(&mut ins.bodies, &ins.contacts, &mut ins.frictions[j..j + 1], k);
            let f = &ins.frictions[j];
            let (ax, ay) = slip(&ins.bodies, f);
            let before = (bx * bx + by * by).sqrt();
            assert!(before > 0.0);
            if expect_stick {
                assert!(f.stuck, "friction {j} slid under {normal} N");
                assert!((ax * ax + ay * ay).sqrt() <= 1e-9 * before, "stuck friction {j} left slip {ax:e}, {ay:e} of {before:e}");
                stuck += 1;
            } else {
                assert!(!f.stuck, "friction {j} stuck under {normal} N");
                assert!((f.force - 0.4 * normal).abs() <= 1e-12 * normal, "sliding force {} != mu N", f.force);
                // the slip it leaves points the same way as before (it is only slowed), so the friction
                // force, which opposes the slip, does no positive work on it
                assert!(ax * bx + ay * by >= 0.0 && (ax * ax + ay * ay).sqrt() <= before, "friction {j} reversed or grew the slip");
                slid += 1;
            }
        }
    }
    assert!(stuck > 0 && slid > 0);
}
