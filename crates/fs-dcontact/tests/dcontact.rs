//! fs-dcontact conformance battery: polyline geometry, the bouncing
//! analytic fixture (restitution + closed-form max penetration),
//! string-fret rattle energy exactness with an explicit-integrator
//! mutation, one-sidedness mutation, iteration budgets across a
//! velocity sweep, collocation refinement convergence, the jawari
//! spectral-enrichment casebook, determinism, and typed refusals.

use fs_dcontact::{ContactStorage, DContactError, Obstacle, polyline_heights, string_collocation};
use fs_phs::{PortHamiltonian, QuadraticStorage, Storage, step};

/// Free 1-DOF mass (ball): H = p^2/(2 m) — no spring.
struct FreeMass {
    m: f64,
}

impl Storage for FreeMass {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        x[1] * x[1] / (2.0 * self.m)
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        out[0] = 0.0;
        out[1] = x[1] / self.m;
    }
}

fn symplectic_j(n_modes: usize) -> Vec<f64> {
    let dim = 2 * n_modes;
    let mut j = vec![0.0; dim * dim];
    for k in 0..n_modes {
        j[(2 * k) * dim + 2 * k + 1] = 1.0;
        j[(2 * k + 1) * dim + 2 * k] = -1.0;
    }
    j
}

/// Modal string as a QuadraticStorage in the interleaved layout
/// (mass-normalized: H = 1/2 sum p^2 + w^2 q^2).
fn string_storage(
    length: f64,
    tension: f64,
    mu: f64,
    n_modes: usize,
) -> (QuadraticStorage, Vec<f64>) {
    let c = fs_math::det::sqrt(tension / mu);
    let pi = core::f64::consts::PI;
    let omegas: Vec<f64> = (1..=n_modes).map(|k| k as f64 * pi * c / length).collect();
    let dim = 2 * n_modes;
    let mut q = vec![0.0; dim * dim];
    for k in 0..n_modes {
        q[(2 * k) * dim + 2 * k] = omegas[k] * omegas[k];
        q[(2 * k + 1) * dim + 2 * k + 1] = 1.0;
    }
    (QuadraticStorage::new(q, dim).expect("q"), omegas)
}

#[test]
fn polyline_heights_match_analytic_geometry() {
    // A tent profile: heights interpolate linearly on each segment.
    let verts = [(0.0, 0.0), (1.0, 2.0), (3.0, 0.0)];
    let h = polyline_heights(&verts, &[0.0, 0.5, 1.0, 2.0, 3.0]).expect("heights");
    let expect = [0.0, 1.0, 2.0, 1.0, 0.0];
    for (a, b) in h.iter().zip(&expect) {
        assert!((a - b).abs() < 1.0e-14);
    }
    // Typed refusals: unsorted, short, out of span.
    assert!(matches!(
        polyline_heights(&[(1.0, 0.0), (0.5, 0.0)], &[0.7]),
        Err(DContactError::Parameter { .. })
    ));
    assert!(polyline_heights(&[(0.0, 0.0)], &[0.0]).is_err());
    assert!(polyline_heights(&verts, &[3.5]).is_err());
}

#[test]
fn bouncing_mass_restitution_and_max_penetration() {
    // Ball (1 kg) falling onto a power-law floor: pure potential
    // contact means restitution EXACTLY 1 (to solver tolerance), and
    // the max penetration follows the closed form
    //   E_impact = w K/(alpha+1) p_max^(alpha+1)
    // up to the gravity work during penetration (bounded in-test).
    let m = 1.0;
    let g = 9.81;
    let (k_c, alpha) = (1.0e6, 1.5);
    let floor = Obstacle::new(
        vec![1.0], // collocation: q IS the displacement toward the floor
        1,
        1,
        vec![0.5], // gap: floor 0.5 m below the start
        vec![1.0],
        k_c,
        alpha,
        "test-fixture: authored (K, alpha), no material claim".to_string(),
    )
    .expect("floor");
    let storage = ContactStorage::new(Box::new(FreeMass { m }), 1, vec![floor]).expect("storage");
    let sys = PortHamiltonian::new(
        2,
        1,
        symplectic_j(1),
        vec![0.0; 4],
        vec![0.0, 1.0],
        Box::new(storage),
    )
    .expect("sys");
    let dt = 1.0e-5;
    let mut x = vec![0.0, 0.0];
    let mut v_impact = 0.0f64;
    let mut v_out = 0.0f64;
    let mut max_pen = 0.0f64;
    let mut in_contact_prev = false;
    for _ in 0..40_000 {
        let rec = step(&sys, &x, &[m * g], dt).expect("step");
        x = rec.x;
        let pen = (x[0] - 0.5).max(0.0);
        max_pen = max_pen.max(pen);
        let in_contact = pen > 0.0;
        if in_contact && !in_contact_prev {
            v_impact = x[1] / m; // velocity entering contact
        }
        if !in_contact && in_contact_prev {
            v_out = x[1] / m; // velocity leaving contact
            break;
        }
        in_contact_prev = in_contact;
    }
    assert!(v_impact > 0.0, "ball never hit the floor");
    let restitution = (v_out / v_impact).abs(); // v_out is negative (rebound)
    assert!(
        (restitution - 1.0).abs() < 1.0e-3,
        "elastic potential restitution {restitution:.6} != 1"
    );
    // Closed-form max penetration from the impact kinetic energy.
    let e_impact = 0.5 * m * v_impact * v_impact;
    let p_closed = fs_math::det::pow(e_impact * (alpha + 1.0) / k_c, 1.0 / (alpha + 1.0));
    // Gravity adds m g p of work during penetration; bound the
    // correction and require agreement within it plus 2%.
    let gravity_correction = m * g * p_closed / e_impact;
    assert!(
        (max_pen - p_closed).abs() <= (0.02 + gravity_correction) * p_closed,
        "max penetration {max_pen:.5e} vs closed form {p_closed:.5e}"
    );
    println!(
        "{{\"suite\":\"fs-dcontact\",\"case\":\"bouncing\",\"restitution\":{restitution:.6},\"p_max\":{max_pen:.4e},\"p_closed\":{p_closed:.4e},\"verdict\":\"pass\"}}"
    );
}

/// The string + flat-fret fixture used by several tests.
fn fret_system(n_modes: usize, gap: f64, k_c: f64) -> (PortHamiltonian, ContactStorage, Vec<f64>) {
    let (length, tension, mu) = (0.65, 70.0, 5.0e-3);
    let (storage, omegas) = string_storage(length, tension, mu, n_modes);
    // 24 fret-line collocation points along the first two thirds.
    let points: Vec<f64> = (1i32..=24).map(|i| length * f64::from(i) / 36.0).collect();
    let coll = string_collocation(length, mu, &points, n_modes).expect("collocation");
    let seg = length / 36.0;
    let fret = Obstacle::new(
        coll.clone(),
        points.len(),
        n_modes,
        vec![gap; points.len()],
        vec![seg; points.len()],
        k_c,
        1.5,
        "test-fixture: authored fret line".to_string(),
    )
    .expect("fret");
    let contact = ContactStorage::new(
        Box::new(string_storage(0.65, 70.0, 5.0e-3, n_modes).0),
        n_modes,
        vec![fret.clone()],
    )
    .expect("contact");
    let sys = PortHamiltonian::new(
        2 * n_modes,
        0,
        symplectic_j(n_modes),
        vec![0.0; 4 * n_modes * n_modes],
        vec![],
        Box::new(ContactStorage::new(Box::new(storage), n_modes, vec![fret]).expect("cs")),
    )
    .expect("sys");
    (sys, contact, omegas)
}

/// Pluck-like initial state: displacement toward the fret in the low
/// modes.
fn pluck(n_modes: usize, amp: f64) -> Vec<f64> {
    let mut x = vec![0.0; 2 * n_modes];
    for k in 0..4.min(n_modes) {
        x[2 * k] = amp / (k + 1) as f64;
    }
    x
}

#[test]
fn fret_rattle_conserves_energy_and_explicit_mutation_grows_it() {
    let n_modes = 8;
    let (sys, contact, _) = fret_system(n_modes, 2.0e-4, 1.0e7);
    let x0 = pluck(n_modes, 6.0e-3);
    let h0 = sys.hamiltonian(&x0);
    let dt = 2.0e-6;
    let steps = 20_000;
    let mut x = x0.clone();
    let mut worst = 0.0f64;
    let mut contact_events = 0usize;
    let mut was_active = false;
    let mut max_pen = 0.0f64;
    for _ in 0..steps {
        let rec = step(&sys, &x, &[], dt).expect("step");
        x = rec.x;
        let probe = contact.probe(&x);
        max_pen = max_pen.max(probe.max_penetration);
        let active = probe.active_points > 0;
        if active && !was_active {
            contact_events += 1;
        }
        was_active = active;
        worst = worst.max((sys.hamiltonian(&x) - h0).abs());
    }
    assert!(
        contact_events >= 3,
        "rattle must actually rattle: {contact_events}"
    );
    assert!(
        worst <= 1.0e-8 * h0,
        "discrete-gradient contact drifted {worst:.3e} of {h0:.3e}"
    );
    // MUTATION: explicit symplectic-Euler on the same vector field at
    // the same dt visibly pumps energy through the stiff contact.
    let (_, contact2, _) = fret_system(n_modes, 2.0e-4, 1.0e7);
    let mut xe = x0.clone();
    let mut grad = vec![0.0; 2 * n_modes];
    let mut h_explicit_max = h0;
    for _ in 0..steps {
        contact2.gradient(&xe, &mut grad);
        for k in 0..n_modes {
            xe[2 * k + 1] -= dt * grad[2 * k];
        }
        contact2.gradient(&xe, &mut grad);
        for k in 0..n_modes {
            xe[2 * k] += dt * grad[2 * k + 1];
        }
        h_explicit_max = h_explicit_max.max(contact2.hamiltonian(&xe));
        if !h_explicit_max.is_finite() {
            break;
        }
    }
    let explicit_growth = (h_explicit_max - h0) / h0;
    assert!(
        !h_explicit_max.is_finite() || explicit_growth > 1.0e3 * (worst / h0).max(1.0e-12),
        "explicit integration must visibly grow energy: growth {explicit_growth:.3e} vs DG drift {:.3e}",
        worst / h0
    );
    println!(
        "{{\"suite\":\"fs-dcontact\",\"case\":\"fret-rattle\",\"events\":{contact_events},\"max_penetration\":{max_pen:.3e},\"dg_drift\":{:.3e},\"explicit_growth\":{explicit_growth:.3e},\"verdict\":\"pass\"}}",
        worst / h0
    );
}

/// Mutant potential WITHOUT the one-sided clamp: separation becomes
/// attraction — the exact bug class the [.]_+ exists to prevent.
struct TwoSidedMutant {
    inner: ContactStorage,
}

impl Storage for TwoSidedMutant {
    fn hamiltonian(&self, x: &[f64]) -> f64 {
        // Same structure energy; contact term evaluated WITHOUT clamp
        // (the |p|-weighted odd form the bug class produces).
        let ob = &self.inner.obstacles()[0];
        let nm = self.inner.n_modes();
        let mut h = self.inner.inner_storage().hamiltonian(x);
        for i in 0..ob.n_points() {
            let mut disp = 0.0;
            for k in 0..nm {
                disp += ob.collocation()[i * nm + k] * x[2 * k];
            }
            let p = disp - ob.gaps()[i]; // NO clamp
            h += ob.weights()[i] * ob.stiffness() / 3.0 * p * p * p.abs();
        }
        h
    }
    fn gradient(&self, x: &[f64], out: &mut [f64]) {
        self.inner.inner_storage().gradient(x, out);
        let ob = &self.inner.obstacles()[0];
        let nm = self.inner.n_modes();
        for i in 0..ob.n_points() {
            let mut disp = 0.0;
            for k in 0..nm {
                disp += ob.collocation()[i * nm + k] * x[2 * k];
            }
            let p = disp - ob.gaps()[i];
            let f = ob.weights()[i] * ob.stiffness() * p * p.abs();
            for k in 0..nm {
                out[2 * k] += f * ob.collocation()[i * nm + k];
            }
        }
    }
}

#[test]
fn dropped_one_sidedness_detected_as_attraction() {
    // Probe: a state SEPARATED from the obstacle must feel zero
    // contact force from the true storage, nonzero (attractive) from
    // the mutant.
    let floor = Obstacle::new(
        vec![1.0],
        1,
        1,
        vec![0.5],
        vec![1.0],
        1.0e6,
        2.0,
        "test".to_string(),
    )
    .expect("floor");
    let real =
        ContactStorage::new(Box::new(FreeMass { m: 1.0 }), 1, vec![floor.clone()]).expect("real");
    let mutant = TwoSidedMutant {
        inner: ContactStorage::new(Box::new(FreeMass { m: 1.0 }), 1, vec![floor]).expect("m"),
    };
    // Sweep separated depths: true force identically zero everywhere
    // on the free side, mutant force present and pulling TOWARD the
    // obstacle (positive q-gradient = restoring toward +q contact).
    for sep in [0.1, 0.3, 0.6] {
        let x_separated = vec![0.5 - sep - 0.5, 0.0]; // sep below the gap
        let mut g_real = vec![0.0; 2];
        let mut g_mut = vec![0.0; 2];
        real.gradient(&x_separated, &mut g_real);
        mutant.gradient(&x_separated, &mut g_mut);
        assert!(
            g_real[0].abs() < f64::MIN_POSITIVE,
            "separated true contact force must vanish at sep {sep}"
        );
        assert!(
            g_mut[0].abs() > 1.0e2,
            "mutant must exert spurious force at sep {sep}: {:.3e}",
            g_mut[0]
        );
    }
}

#[test]
fn iteration_budget_held_across_velocity_sweep_and_stall_is_typed() {
    let n_modes = 8;
    let (sys, _, _) = fret_system(n_modes, 2.0e-4, 1.0e7);
    let dt = 2.0e-6;
    let mut histogram = [0usize; 8]; // buckets of 4 iterations
    let mut max_iters = 0usize;
    for amp in [2.0e-3, 6.0e-3, 1.5e-2] {
        let mut x = pluck(n_modes, amp);
        for _ in 0..4000 {
            let rec = step(&sys, &x, &[], dt).expect("step within budget");
            max_iters = max_iters.max(rec.newton_iters);
            histogram[(rec.newton_iters / 4).min(7)] += 1;
            x = rec.x;
        }
    }
    // Interior bound (review finding: <= 50 is vacuous — a successful
    // step CANNOT report more than the fs-phs cap; a stall is already
    // a panic above). Histogram-measured headroom at authoring: all
    // steps <= 8 iterations.
    assert!(max_iters <= 24, "iteration budget regressed: {max_iters}");
    println!(
        "{{\"suite\":\"fs-dcontact\",\"case\":\"iteration-budget\",\"max_iters\":{max_iters},\"histogram\":{histogram:?},\"verdict\":\"pass\"}}"
    );
    // The stall path is TYPED: absurd stiffness with a huge step
    // refuses by name instead of looping or lying.
    let (sys_stiff, _, _) = fret_system(4, 1.0e-6, 1.0e18);
    let x = pluck(4, 5.0e-2);
    let out = step(&sys_stiff, &x, &[], 1.0e-2);
    assert!(
        matches!(out, Err(fs_phs::PhsError::NewtonStalled { .. })),
        "expected NewtonStalled, got {out:?}"
    );
}

#[test]
fn collocation_refinement_converges() {
    // Contact energy of a FIXED penetrating string shape under a flat
    // obstacle: refining the collocation quadrature converges (the
    // 8->32 difference dominates the 32->128 difference).
    let (length, mu) = (0.65, 5.0e-3);
    let n_modes = 4;
    let shape = |x: &mut Vec<f64>| {
        x[0] = 8.0e-3; // fundamental pushed into the obstacle
        x[4] = 2.0e-3;
    };
    let energy_at = |n_pts: usize| -> f64 {
        let pts: Vec<f64> = (1..=n_pts)
            .map(|i| length * i as f64 / (n_pts as f64 + 1.0))
            .collect();
        let coll = string_collocation(length, mu, &pts, n_modes).expect("collocation");
        let seg = length / (n_pts as f64 + 1.0);
        let ob = Obstacle::new(
            coll,
            n_pts,
            n_modes,
            vec![1.0e-4; n_pts],
            vec![seg; n_pts],
            1.0e7,
            1.5,
            "test".to_string(),
        )
        .expect("ob");
        let cs = ContactStorage::new(
            Box::new(string_storage(length, 70.0, mu, n_modes).0),
            n_modes,
            vec![ob],
        )
        .expect("cs");
        let mut x = vec![0.0; 2 * n_modes];
        shape(&mut x);
        cs.probe(&x).contact_energy
    };
    let (e8, e32, e128) = (energy_at(8), energy_at(32), energy_at(128));
    let d_coarse = (e32 - e8).abs();
    let d_fine = (e128 - e32).abs();
    assert!(e128 > 0.0, "fixture must penetrate");
    // Measured ratio 0.0099 at authoring (near O(h^2)); 0.1 keeps 10x
    // headroom while catching an O(h)-degradation (ratio ~0.25).
    assert!(
        d_fine < 0.1 * d_coarse.max(1.0e-30),
        "refinement not converging: {e8:.6e} {e32:.6e} {e128:.6e}"
    );
}

#[test]
fn jawari_bridge_enriches_high_band_vs_clean_termination() {
    // The sitar/tanpura jawari: a GRADED bridge surface near the
    // termination repeatedly grazes the string and pumps energy
    // upward. The CONTROL is the clean knife-edge termination — the
    // plain string with NO collision surface (a first attempt used a
    // tiny-gap rattling point as "hard bridge" and measured MORE high
    // band than the jawari, executed: a rattle point is itself a
    // collision exciter, not a termination). Quantified as the
    // high-band energy fraction after the transient.
    let (length, tension, mu) = (0.65, 70.0, 5.0e-3);
    // n = 12 modes / 16k steps at 2.5 us: ~4 fundamental periods of
    // grazing — enough for the enrichment signature while keeping the
    // debug-build implicit solve inside the suite budget (n = 20 with
    // 60k steps measured ~25 minutes).
    let n_modes = 12;
    let dt = 2.5e-6;
    let steps = 16_000;
    let run = |obstacle: Obstacle| -> Vec<f64> {
        let (storage, omegas) = string_storage(length, tension, mu, n_modes);
        let cs = ContactStorage::new(Box::new(storage), n_modes, vec![obstacle]).expect("cs");
        let sys = PortHamiltonian::new(
            2 * n_modes,
            0,
            symplectic_j(n_modes),
            vec![0.0; 4 * n_modes * n_modes],
            vec![],
            Box::new(cs),
        )
        .expect("sys");
        // Pluck AWAY from the obstacle side: the string swings back
        // through neutral and GRAZES the bridge each period (a
        // positive pluck started 20 gap-depths INSIDE the profile —
        // executed unphysical preload, energies 100x the pluck's).
        let mut x = pluck(n_modes, -4.0e-3);
        for _ in 0..steps {
            x = step(&sys, &x, &[], dt).expect("step").x;
        }
        (0..n_modes)
            .map(|k| {
                f64::midpoint(
                    x[2 * k + 1] * x[2 * k + 1],
                    omegas[k] * omegas[k] * x[2 * k] * x[2 * k],
                )
            })
            .collect()
    };
    // Jawari: parabolic graded profile over x in [0.02, 0.08] m with
    // gaps from 0.15 mm (near the end) to 1.5 mm, via the polyline
    // helper.
    let jaw_pts: Vec<f64> = (0i32..12)
        .map(|i| 0.02 + 0.06 * f64::from(i) / 11.0)
        .collect();
    let profile: Vec<(f64, f64)> = (0i32..12)
        .map(|i| {
            let t = f64::from(i) / 11.0;
            (0.02 + 0.06 * t, 1.5e-4 + 1.35e-3 * t * t)
        })
        .collect();
    let gaps = polyline_heights(&profile, &jaw_pts).expect("gaps");
    let seg = 0.06 / 12.0;
    let jawari = Obstacle::new(
        string_collocation(length, mu, &jaw_pts, n_modes).expect("collocation"),
        jaw_pts.len(),
        n_modes,
        gaps,
        vec![seg; jaw_pts.len()],
        1.0e8,
        // alpha = 2 (C^2 force): alpha = 1.5 has an UNBOUNDED contact
        // Hessian at the boundary (d2f ~ p^-0.5) and the tight-graze
        // jawari regime stalled the FD-Jacobian Newton on it
        // (executed); exponent 2 is inside the published musical-
        // collision range and removes the singularity.
        2.0,
        "test-fixture: jawari-class graded profile (authored)".to_string(),
    )
    .expect("jawari");
    // Clean termination: zero contact stiffness = the plain linear
    // string (the model's fixed ends ARE the knife-edge bridge).
    let clean = Obstacle::new(
        string_collocation(length, mu, &[0.02], n_modes).expect("collocation"),
        1,
        n_modes,
        vec![1.5e-4],
        vec![1.0],
        0.0,
        2.0,
        "test-fixture: zero-stiffness control (clean termination)".to_string(),
    )
    .expect("clean");
    let e_jaw = run(jawari);
    let e_hard = run(clean);
    let frac = |e: &[f64]| -> f64 {
        let hi: f64 = e[6..].iter().sum();
        let total: f64 = e.iter().sum();
        hi / total.max(f64::MIN_POSITIVE)
    };
    let (f_jaw, f_clean) = (frac(&e_jaw), frac(&e_hard));
    // The linear control PRESERVES the pluck's band split exactly
    // (nothing couples modes), so its high fraction is the initial
    // one; the jawari must sit far above it AND carry an absolute
    // high-band share.
    assert!(
        f_jaw > 10.0 * f_clean.max(1.0e-9) && f_jaw > 0.05,
        "jawari must enrich the high band: {f_jaw:.4} vs clean {f_clean:.4}"
    );
    println!(
        "{{\"suite\":\"fs-dcontact\",\"case\":\"jawari\",\"high_band_fraction_jawari\":{f_jaw:.4},\"high_band_fraction_clean\":{f_clean:.4},\"band_energies_jawari\":{e_jaw:?},\"band_energies_clean\":{e_hard:?},\"verdict\":\"pass\"}}"
    );
}

#[test]
fn contact_gradient_matches_finite_difference_of_h() {
    // THE oracle the conservation test cannot be (review finding: the
    // Gonzalez correction forces dg.dx = dH for ANY gradient, so
    // energy conservation is structurally blind to contact-gradient
    // errors): central finite differences of the coded H must match
    // gradient() on a MULTI-POINT, NON-UNIFORM-WEIGHT obstacle at a
    // state with a mixed active/inactive contact set.
    let (length, mu) = (0.65, 5.0e-3);
    let n_modes = 5;
    let pts: Vec<f64> = (1i32..=9).map(|i| length * f64::from(i) / 10.0).collect();
    let coll = string_collocation(length, mu, &pts, n_modes).expect("collocation");
    let weights: Vec<f64> = (1i32..=9).map(|i| 0.01 * f64::from(i)).collect(); // non-uniform
    let gaps: Vec<f64> = (0..9)
        .map(|i| if i % 2 == 0 { 1.0e-4 } else { 5.0e-2 }) // half active, half far
        .collect();
    let ob = Obstacle::new(
        coll,
        9,
        n_modes,
        gaps,
        weights,
        3.0e6,
        1.7,
        "test".to_string(),
    )
    .expect("ob");
    let cs = ContactStorage::new(
        Box::new(string_storage(length, 70.0, mu, n_modes).0),
        n_modes,
        vec![ob],
    )
    .expect("cs");
    let mut x = vec![0.0; 2 * n_modes];
    for k in 0..n_modes {
        x[2 * k] = 3.0e-3 / (k + 1) as f64;
        x[2 * k + 1] = 0.1 * (k as f64 - 2.0);
    }
    let probe = cs.probe(&x);
    assert!(
        probe.active_points >= 2 && probe.active_points <= 8,
        "fixture needs a MIXED active set, got {}",
        probe.active_points
    );
    let mut g = vec![0.0; 2 * n_modes];
    cs.gradient(&x, &mut g);
    let scale = g.iter().fold(1.0e-30f64, |a, &v| a.max(v.abs()));
    for i in 0..2 * n_modes {
        let h_step = 1.0e-7 * (1.0 + x[i].abs());
        let mut xp = x.clone();
        let mut xm = x.clone();
        xp[i] += h_step;
        xm[i] -= h_step;
        let fd = (cs.hamiltonian(&xp) - cs.hamiltonian(&xm)) / (2.0 * h_step);
        assert!(
            (fd - g[i]).abs() <= 1.0e-5 * scale,
            "gradient[{i}] {:.6e} vs FD {fd:.6e}",
            g[i]
        );
    }
    // Probe-vs-Hamiltonian consistency (review finding: duplicated
    // potential formula): contact energy equals the storage split.
    let inner_h = string_storage(length, 70.0, mu, n_modes).0.hamiltonian(&x);
    assert!(
        (probe.contact_energy - (cs.hamiltonian(&x) - inner_h)).abs()
            <= 1.0e-12 * probe.contact_energy.max(1.0),
        "probe energy diverged from the Hamiltonian split"
    );
}

#[test]
fn determinism_bitwise() {
    let n_modes = 6;
    let (sys, _, _) = fret_system(n_modes, 2.0e-4, 1.0e7);
    let x0 = pluck(n_modes, 6.0e-3);
    let runs: Vec<Vec<f64>> = (0..2)
        .map(|_| {
            let mut x = x0.clone();
            for _ in 0..2000 {
                x = step(&sys, &x, &[], 2.0e-6).expect("step").x;
            }
            x
        })
        .collect();
    for (a, b) in runs[0].iter().zip(&runs[1]) {
        assert_eq!(a.to_bits(), b.to_bits(), "nondeterministic step");
    }
}

#[test]
fn refusals_are_typed() {
    assert!(matches!(
        Obstacle::new(
            vec![1.0],
            1,
            1,
            vec![0.1],
            vec![1.0],
            -1.0,
            1.5,
            String::new()
        ),
        Err(DContactError::Parameter { .. })
    ));
    assert!(matches!(
        Obstacle::new(
            vec![1.0],
            1,
            1,
            vec![0.1],
            vec![1.0],
            1.0,
            0.5,
            String::new()
        ),
        Err(DContactError::Parameter { .. })
    ));
    assert!(matches!(
        Obstacle::new(
            vec![1.0],
            1,
            1,
            vec![0.1],
            vec![1.0],
            f64::NAN,
            1.5,
            String::new()
        ),
        Err(DContactError::Parameter { .. })
    ));
    assert!(
        Obstacle::new(
            vec![1.0, 2.0],
            1,
            1,
            vec![0.1],
            vec![1.0],
            1.0,
            1.5,
            String::new()
        )
        .is_err()
    );
}
