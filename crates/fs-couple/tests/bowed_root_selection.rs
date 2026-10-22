//! M4.3-W2 (IRP pexfp.152): the set-valued bow contact discloses its roots.
//!
//! `StribeckLoadLine` (e9133f454) applies McIntyre, Schumacher & Woodhouse
//! (1983) branch hysteresis where the discrete load line meets the friction
//! graph more than once. Every step now reports how many intersections there
//! are and which one it took, and the bow runtime keeps that receipt for every
//! substep. The fixture string is an exact-ZOH modal string bowed at its
//! Schelleng mid-range force (Euphonics 9.3.1 eqs. 4 and 9).

use fs_couple::bowed_string::runtime::{BowContactReceipt, BowedStringState};
use fs_couple::bowed_string::{
    BowGesture, BowedRunConfig, BowedStringCard, FrictionIsland, Termination,
};
use fs_couple::stribeck_friction::{ContactBranch, StribeckFriction, StribeckLoadLine};
use fs_math::det;
use std::fs;
use std::path::{Path, PathBuf};

// String fixture values follow the cello D string of the Galluzzo bowing rig
// (Euphonics 9.6.3 Table 2: 146.8 Hz, 111 N, 2.7 g/m). They are used here as
// fixture values for an ideal string.
const F0_HZ: f64 = 146.8;
const TENSION_N: f64 = 111.0;
const DENSITY_KG_M: f64 = 2.7e-3;
/// Bridge resistance over string impedance; Schelleng's analysis needs
/// `R >> Z0`, and beta >= Z0/(4R) for a Helmholtz region at beta = 1/20.
const R_OVER_Z0: f64 = 20.0;
const BOW_SPEED_M_S: f64 = 0.1;
/// Static and asymptotic kinetic limits of Smith & Woodhouse (2000) Eq. (2).
const MU_S: f64 = 1.2;
const MU_D: f64 = 0.35;
/// Stribeck decay width, a fixture value well below every Helmholtz slip
/// speed `v_b (1/beta - 1)`, so the slip coefficient is mu_d.
const DECAY_M_S: f64 = 0.05;
fn z0() -> f64 {
    det::sqrt(TENSION_N * DENSITY_KG_M)
}
fn wave_speed() -> f64 {
    det::sqrt(TENSION_N / DENSITY_KG_M)
}
fn length_m() -> f64 {
    wave_speed() / (2.0 * F0_HZ)
}
fn resistance() -> f64 {
    R_OVER_Z0 * z0()
}
fn schelleng_fmax(beta: f64) -> f64 {
    2.0 * z0() * BOW_SPEED_M_S / (beta * (MU_S - MU_D))
}
fn schelleng_fmin(beta: f64) -> f64 {
    z0() * z0() * BOW_SPEED_M_S / (2.0 * resistance() * beta * beta * (MU_S - MU_D))
}
fn nominal_force(beta: f64) -> f64 {
    det::sqrt(schelleng_fmin(beta) * schelleng_fmax(beta))
}

fn card(sample_rate_hz: u32, modes: usize) -> BowedStringCard {
    let (c, l) = (wave_speed(), length_m());
    let sigma = c / (2.0 * l) * det::ln((R_OVER_Z0 + 1.0) / (R_OVER_Z0 - 1.0));
    BowedStringCard {
        length_m: l,
        tension_n: TENSION_N,
        linear_density_kg_m: DENSITY_KG_M,
        bending_stiffness_n_m2: 0.0,
        viscous_bending_n_m2_s: 0.0,
        mode_count: modes,
        zetas: (1..=modes)
            .map(|k| {
                #[allow(clippy::cast_precision_loss)]
                let omega = k as f64 * core::f64::consts::PI * c / l;
                sigma / omega
            })
            .collect(),
        sample_rate_hz,
    }
}

fn config(
    sample_rate_hz: u32,
    modes: usize,
    substeps: usize,
    beta: f64,
    force: f64,
) -> BowedRunConfig {
    BowedRunConfig {
        card: card(sample_rate_hz, modes),
        island: FrictionIsland::Stribeck(
            StribeckFriction::try_new(MU_S, MU_D, DECAY_M_S, DECAY_M_S).unwrap(),
        ),
        gesture: BowGesture::admit(BOW_SPEED_M_S, force, beta).unwrap(),
        steps: 0,
        subsamples: substeps,
        termination: Termination::Rigid,
        listener_m: 1.0,
    }
}

fn stribeck() -> fs_tribo::FrictionLaw {
    fs_tribo::FrictionLaw::Stribeck {
        static_mu: MU_S,
        kinetic_mu: MU_D,
        characteristic_speed: DECAY_M_S,
        viscous_per_speed: 0.0,
    }
}

#[test]
fn g1_three_roots_hysteresis_follows_mcintyre_rule() {
    // Load line v = s - m F(v) at N = 1 N and m = 1/(2 Z0), with three
    // intersections over a window of free slip s. Independent oracle: every
    // sign change of g(v) = v - s + m N mu(v) on a dense grid, each refined
    // by bisection, plus the cone when s <= m N mu_s. MSW hysteresis: the
    // contact enters the window stuck on the way up and leaves it on the
    // high branch on the way down.
    let (normal, mobility) = (1.0, 1.0 / (2.0 * z0()));
    let law = stribeck();
    let line = StribeckLoadLine::new(&law, normal, mobility).unwrap();
    let mu = |v: f64| {
        if v == 0.0 {
            MU_S
        } else {
            law.kinetic_coefficient(v).unwrap()
        }
    };
    let oracle = |s: f64| -> (bool, Vec<f64>) {
        let g = |v: f64| v - s + mobility * normal * mu(v);
        let cells = 100_000_u32;
        let mut roots = Vec::new();
        let mut prev = (0.0, g(0.0));
        for i in 1..=cells {
            let v = s * f64::from(i) / f64::from(cells);
            let cur = g(v);
            if (prev.1 < 0.0) != (cur < 0.0) {
                let (mut lo, mut hi) = (prev.0, v);
                let below = prev.1 < 0.0;
                loop {
                    let mid = lo + 0.5 * (hi - lo);
                    if mid <= lo || mid >= hi {
                        break;
                    }
                    if (g(mid) < 0.0) == below {
                        lo = mid;
                    } else {
                        hi = mid;
                    }
                }
                roots.push(hi);
            }
            prev = (v, cur);
        }
        (s <= mobility * normal * MU_S, roots)
    };
    let speeds: Vec<f64> = (1..=300).map(|i| 0.01 * f64::from(i)).collect();
    let mut branch = ContactBranch::Fresh;
    let (mut up_window, mut down_window) = (Vec::new(), Vec::new());
    let mut worst = 0.0_f64;
    for (pass, sweep) in [speeds.clone(), speeds.iter().rev().copied().collect()]
        .into_iter()
        .enumerate()
    {
        for s in sweep {
            let step = line.step(s, branch).unwrap();
            let (cone, roots) = oracle(s);
            let expected = usize::from(cone) + roots.len();
            assert_eq!(
                usize::from(step.roots),
                expected,
                "intersections at s = {s}"
            );
            assert!(
                step.selected < step.roots,
                "index {} of {}",
                step.selected,
                step.roots
            );
            if step.branch == ContactBranch::Stuck {
                assert!(
                    cone && step.selected == 0,
                    "stuck outside the cone at s = {s}"
                );
            } else {
                let slip = step.relative_velocity_m_s;
                let index = usize::from(step.selected) - usize::from(cone);
                let e = (slip - roots[index]).abs() / roots[index];
                worst = worst.max(e);
            }
            if expected == 3 {
                if pass == 0 {
                    up_window.push(step.selected)
                } else {
                    down_window.push(step.selected)
                }
            }
            branch = step.branch;
        }
    }
    println!(
        "{{\"test\":\"three_roots_msw\",\"normal_n\":{normal},\"mobility\":{mobility},\"critical_slips\":{:?},\"three_root_steps_up\":{},\"three_root_steps_down\":{},\"up_selected\":{:?},\"down_selected\":{:?},\"root_rel_err\":{worst:e},\"tol\":1e-12}}",
        line.critical_slips_m_s(),
        up_window.len(),
        down_window.len(),
        up_window.iter().collect::<std::collections::BTreeSet<_>>(),
        down_window
            .iter()
            .collect::<std::collections::BTreeSet<_>>()
    );
    assert!(
        !up_window.is_empty() && !down_window.is_empty(),
        "the fixture must have a three-root window"
    );
    assert!(
        up_window.iter().all(|&i| i == 0),
        "rising: stuck through the window {up_window:?}"
    );
    assert!(
        down_window.iter().all(|&i| i == 2),
        "falling: high branch through the window {down_window:?}"
    );
    assert!(worst <= 1e-12, "root values {worst:e}");
}

#[test]
fn g0_root_disclosure_in_every_receipt() {
    let beta = 0.1;
    let mut state =
        BowedStringState::new(&config(48_000, 64, 3, beta, nominal_force(beta)), 1).unwrap();
    let mut counts = [0_usize; 4];
    for _ in 0..4_000 {
        let sample = state.step().unwrap();
        let receipts = state.last_contact_receipts();
        assert_eq!(receipts.len(), 3, "one receipt per substep");
        assert_eq!(sample.contact, receipts[2]);
        let stuck = receipts
            .iter()
            .filter(|r| r.branch == ContactBranch::Stuck)
            .count();
        assert_eq!(sample.stuck_substeps as usize, stuck);
        for r in receipts {
            assert!((1..=3).contains(&r.roots) && r.selected < r.roots, "{r:?}");
            match r.branch {
                ContactBranch::Stuck => assert_eq!(r.selected, 0),
                ContactBranch::SlidingHigh { .. } => assert_eq!(r.selected, r.roots - 1),
                ContactBranch::SlidingLow { .. } => {}
                ContactBranch::Fresh => panic!("a loaded contact is never fresh"),
            }
            counts[usize::from(r.roots)] += 1;
        }
    }
    state.set_bow(BOW_SPEED_M_S, 0.0, beta).unwrap();
    state.step().unwrap();
    assert!(
        state
            .last_contact_receipts()
            .iter()
            .all(|r| *r == BowContactReceipt::default())
    );
    let mut viscous = config(48_000, 64, 2, beta, 1.0);
    viscous.island = FrictionIsland::ViscousOnly {
        viscous_n_s_per_m: 0.5,
    };
    let mut falsifier = BowedStringState::new(&viscous, 1).unwrap();
    falsifier.step().unwrap();
    assert!(
        falsifier
            .last_contact_receipts()
            .iter()
            .all(|r| r.roots == 1 && r.selected == 0)
    );
    println!(
        "{{\"test\":\"receipt_disclosure\",\"substeps\":3,\"root_count_histogram\":{counts:?}}}"
    );
}

fn rust_files(dir: &Path, out: &mut Vec<PathBuf>) {
    let Ok(entries) = fs::read_dir(dir) else {
        return;
    };
    for entry in entries.flatten() {
        let path = entry.path();
        if path.is_dir() {
            rust_files(&path, out);
        } else if path.extension().is_some_and(|e| e == "rs") {
            out.push(path);
        }
    }
}

#[test]
fn g0_no_hand_typed_capture_constants() {
    // BR-04 #10: the bow path once decided capture with
    // max(stiction, 0.02 m/s) and a 0.15 m/s basin.
    let src = Path::new(env!("CARGO_MANIFEST_DIR")).join("src");
    let mut files = vec![
        src.join("bowed_string.rs"),
        src.join("stribeck_friction.rs"),
    ];
    rust_files(&src.join("bowed_string"), &mut files);
    for file in &files {
        let text = fs::read_to_string(file).unwrap();
        for needle in ["CAPTURE_BASIN", "capture_tol", "max(0.02)", ".max(0.02"] {
            assert!(
                !text.contains(needle),
                "{} contains {needle}",
                file.display()
            );
        }
    }
    println!(
        "{{\"test\":\"no_capture_constants\",\"files_scanned\":{}}}",
        files.len()
    );
}
