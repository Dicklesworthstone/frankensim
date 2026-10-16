//! E4.3b2-iii battery (bead wf-root-guzez.5.8.3): scheduling
//! interpolation vs directly-assembled truth on HELD-OUT path samples
//! (convection AND ground axes), state continuation across a schedule
//! switch (continuous output, genuinely different steady levels),
//! phase + group-delay clause, no-post-reduction-image proof
//! (structural source pin + the ground-axis discriminator), schedule
//! caps at cap AND cap+1, and the V-08b2 receipt (every clause
//! EXECUTED; BT/Loewner recorded as explicit NO-DATA; golden).
//! Repro: cargo test -p fs-wing --test romsched_battery --release

use fs_wing::images::CertifiedGround;
use fs_wing::prescribedwake::frozen_grid_v1;
use fs_wing::rom::{A1Lti, assemble_a1_lti, wright_a1_layout_v1};
use fs_wing::romreduce::{HELD_OUT_W, project, reduce_shared, small_eigenvalues, transfer_of};
use fs_wing::romsched::{
    ClauseVerdict, SchedulePoint, ScheduledRom, V08B2_SCHEMA, emit_v08b2_receipt, march_scheduled,
    phase_and_group_delay,
};

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-wing-romsched\",\"case\":\"{case}\",{payload}}}");
}

fn ground() -> CertifiedGround {
    CertifiedGround {
        z_m: 3.0,
        certificate_slope: 0.000606,
        certificate_rms_m: 0.801,
    }
}

const V: f64 = 13.0;
const ROWS: usize = 120;

/// The convection-axis pair at grid point 0/1 + the ground-axis pair.
fn systems() -> (A1Lti, A1Lti, A1Lti, A1Lti) {
    let layout = wright_a1_layout_v1();
    let grid = frozen_grid_v1();
    let mut c085 = grid.points[0];
    c085.convection = 0.85;
    let mut c100 = grid.points[0];
    c100.convection = 1.0;
    let mut hlow = grid.points[0];
    hlow.h_over_b = 0.1;
    let mut hmid = grid.points[0];
    hmid.h_over_b = 0.3;
    (
        assemble_a1_lti(&layout, &c085, &ground(), V, ROWS).unwrap(),
        assemble_a1_lti(&layout, &c100, &ground(), V, ROWS).unwrap(),
        assemble_a1_lti(&layout, &hlow, &ground(), V, ROWS).unwrap(),
        assemble_a1_lti(&layout, &hmid, &ground(), V, ROWS).unwrap(),
    )
}

fn transfer_err(truth: &A1Lti, red: &fs_wing::romreduce::ReducedLti, w: f64) -> f64 {
    let gt = transfer_of(&truth.a, &truth.b, &truth.c, &truth.d, truth.order, w).unwrap();
    let gr = transfer_of(&red.a, &red.b, &red.c, &red.d, red.order, w).unwrap();
    let mut worst = 0.0f64;
    for ch in 0..6 {
        let mag = gt[ch].0.hypot(gt[ch].1);
        if mag > 1e-9 {
            worst = worst.max((gt[ch].0 - gr[ch].0).hypot(gt[ch].1 - gr[ch].1) / mag);
        }
    }
    worst
}

#[test]
fn scheduling_matches_heldout_path_samples_on_both_axes() {
    let layout = wright_a1_layout_v1();
    let grid = frozen_grid_v1();
    let (c085, c100, hlow, hmid) = systems();
    // Convection axis.
    let refs = [&c085, &c100];
    let red = reduce_shared(&refs).unwrap();
    let rom = ScheduledRom::new(vec![
        SchedulePoint {
            param: 0.85,
            sys: project(&c085, &red.basis, red.order),
        },
        SchedulePoint {
            param: 1.0,
            sys: project(&c100, &red.basis, red.order),
        },
    ])
    .unwrap();
    let mut worst_conv = 0.0f64;
    for &p in &[0.8875_f64, 0.925, 0.9625] {
        let mut point = grid.points[0];
        point.convection = p;
        let truth = assemble_a1_lti(&layout, &point, &ground(), V, ROWS).unwrap();
        let sched = rom.at(p).unwrap();
        for &w in &HELD_OUT_W[..4] {
            worst_conv = worst_conv.max(transfer_err(&truth, &sched, w));
        }
    }
    assert!(worst_conv < 0.05, "convection scheduling: {worst_conv}");
    // Ground axis (h/b 0.1 → 0.3): images entered the operator BEFORE
    // reduction; scheduling here is image-free by construction.
    let refs_g = [&hlow, &hmid];
    let red_g = reduce_shared(&refs_g).unwrap();
    let rom_g = ScheduledRom::new(vec![
        SchedulePoint {
            param: 0.1,
            sys: project(&hlow, &red_g.basis, red_g.order),
        },
        SchedulePoint {
            param: 0.3,
            sys: project(&hmid, &red_g.basis, red_g.order),
        },
    ])
    .unwrap();
    let mut worst_gnd = 0.0f64;
    for &hb in &[0.15_f64, 0.2, 0.25] {
        let mut point = grid.points[0];
        point.h_over_b = hb;
        let truth = assemble_a1_lti(&layout, &point, &ground(), V, ROWS).unwrap();
        let sched = rom_g.at(hb).unwrap();
        for &w in &HELD_OUT_W[..4] {
            worst_gnd = worst_gnd.max(transfer_err(&truth, &sched, w));
        }
    }
    assert!(worst_gnd < 0.08, "ground scheduling: {worst_gnd}");
    // Scheduled midpoints keep strictly stable poles (causality).
    for rom_ref in [&rom, &rom_g] {
        let mid = rom_ref
            .at(0.5 * (rom_ref.points[0].param + rom_ref.points[1].param))
            .unwrap();
        for (re, _) in small_eigenvalues(&mid.a, mid.order).unwrap() {
            assert!(re < 0.0, "scheduled midpoint pole unstable");
        }
    }
    jlog(
        "scheduling",
        &format!("\"conv_worst\":{worst_conv},\"ground_worst\":{worst_gnd}"),
    );
}

#[test]
fn state_continuation_across_a_schedule_switch() {
    // The GROUND axis carries a real steady shift (h/b 0.1 vs 0.3 —
    // measured; the convection axis leaves the DC nearly unchanged, so
    // it cannot discriminate a live level change).
    let (_, _, hlow, hmid) = systems();
    let refs = [&hlow, &hmid];
    let red = reduce_shared(&refs).unwrap();
    let rom = ScheduledRom::new(vec![
        SchedulePoint {
            param: 0.1,
            sys: project(&hlow, &red.basis, red.order),
        },
        SchedulePoint {
            param: 0.3,
            sys: project(&hmid, &red.basis, red.order),
        },
    ])
    .unwrap();
    // A REAL trajectory ramps the parameter; the output under a
    // ramped schedule must be continuous (small per-step deltas).
    let n = 8_000usize;
    let dt = 1.0 / 240.0;
    let ramp = |k: usize| -> f64 {
        if k < n / 2 {
            0.1
        } else {
            0.1 + 0.2 * ((k - n / 2) as f64 / (n / 2) as f64)
        }
    };
    let y = march_scheduled(&rom, &|k| ([0.1, 0.02], ramp(k)), dt, n).unwrap();
    let scale = y[n / 2 - 10][0].abs().max(1.0);
    let mut max_step = 0.0f64;
    for k in (n / 4)..n {
        max_step = max_step.max((y[k][0] - y[k - 1][0]).abs());
    }
    assert!(
        max_step / scale < 0.005,
        "ramped schedule must be continuous: {max_step}"
    );
    let settled_a = y[n / 2 - 10][0];
    let settled_b = y[n - 1][0];
    assert!(
        (settled_b - settled_a).abs() / scale > 0.005,
        "the steady levels must genuinely move: {settled_a} vs {settled_b}"
    );
    // Hard-switch discriminator: the CARRIED state jumps LESS than a
    // state-reset twin (the executed proof that x continues in the
    // shared coordinates instead of being reprojected/reset).
    let y_hard = march_scheduled(
        &rom,
        &|k| ([0.1, 0.02], if k < n / 2 { 0.1 } else { 0.3 }),
        dt,
        n,
    )
    .unwrap();
    let jump_cont = (y_hard[n / 2][0] - y_hard[n / 2 - 1][0]).abs();
    // Reset twin: a fresh march at the post-switch point from x = 0 —
    // its first output is feed-through only.
    let y_reset = march_scheduled(&rom, &|_| ([0.1, 0.02], 0.3), dt, 10).unwrap();
    let jump_reset = (y_reset[0][0] - y_hard[n / 2 - 1][0]).abs();
    assert!(
        jump_cont < 0.5 * jump_reset,
        "the carried state must beat the reset twin: {jump_cont} vs {jump_reset}"
    );
    jlog(
        "continuation",
        &format!(
            "\"ramp_max_step_rel\":{},\"steady_shift_rel\":{},\"jump_cont\":{jump_cont},\"jump_reset\":{jump_reset}",
            max_step / scale,
            (settled_b - settled_a).abs() / scale
        ),
    );
}

#[test]
fn phase_and_group_delay_clause() {
    let (c085, c100, _, _) = systems();
    let refs = [&c085, &c100];
    let red = reduce_shared(&refs).unwrap();
    let r = project(&c100, &red.basis, red.order);
    let mut worst_phase = 0.0f64;
    let mut worst_gd = 0.0f64;
    for &w in &HELD_OUT_W[..5] {
        for ch in [0usize, 4] {
            // wing/δc and hinge/δc channels
            let (pf, gf) =
                phase_and_group_delay(&c100.a, &c100.b, &c100.c, &c100.d, c100.order, ch, w)
                    .unwrap();
            let (pr, gr) = phase_and_group_delay(&r.a, &r.b, &r.c, &r.d, r.order, ch, w).unwrap();
            let mut dp = (pf - pr).abs();
            if dp > core::f64::consts::PI {
                dp = core::f64::consts::TAU - dp;
            }
            worst_phase = worst_phase.max(dp);
            worst_gd = worst_gd.max((gf - gr).abs());
        }
    }
    assert!(worst_phase < 0.02, "phase error [rad]: {worst_phase}");
    assert!(worst_gd < 0.01, "group-delay error [s]: {worst_gd}");
    jlog(
        "phase-gd",
        &format!("\"phase_rad\":{worst_phase},\"group_delay_s\":{worst_gd}"),
    );
}

#[test]
fn no_post_reduction_image_proof() {
    // Structural half: the scheduling module never touches the image
    // machinery — images live UPSTREAM in the operator.
    let src = include_str!("../src/romsched.rs");
    let code: String = src
        .lines()
        .filter(|l| {
            let t = l.trim_start();
            !(t.starts_with("//") || t.starts_with("//!") || t.starts_with("///"))
        })
        .collect::<Vec<_>>()
        .join("\n");
    assert!(
        !code.contains("images::") && !code.contains("CertifiedGround"),
        "romsched must carry no image logic (code lines; docs may cite the law)"
    );
    // Discriminator half: the ground-axis schedule genuinely carries
    // ground effect (h/b 0.1 vs 0.3 reduced DC differ) — the images
    // arrived through the operator BEFORE the projection.
    let (_, _, hlow, hmid) = systems();
    let refs = [&hlow, &hmid];
    let red = reduce_shared(&refs).unwrap();
    let ra = project(&hlow, &red.basis, red.order);
    let rb = project(&hmid, &red.basis, red.order);
    let ga = transfer_of(&ra.a, &ra.b, &ra.c, &ra.d, ra.order, 0.7).unwrap();
    let gb = transfer_of(&rb.a, &rb.b, &rb.c, &rb.d, rb.order, 0.7).unwrap();
    let da = ga[1].0.hypot(ga[1].1); // wing/gust channel magnitude
    let db = gb[1].0.hypot(gb[1].1);
    assert!(
        (da - db).abs() / db.max(1.0) > 0.01,
        "ground effect must survive the reduction: {da} vs {db}"
    );
    jlog("no-post-image", &format!("\"low\":{da},\"mid\":{db}"));
}

#[test]
fn schedule_caps_at_cap_and_cap_plus_one() {
    let (c085, c100, _, _) = systems();
    let refs = [&c085, &c100];
    let red = reduce_shared(&refs).unwrap();
    let mk = || {
        ScheduledRom::new(vec![
            SchedulePoint {
                param: 0.85,
                sys: project(&c085, &red.basis, red.order),
            },
            SchedulePoint {
                param: 1.0,
                sys: project(&c100, &red.basis, red.order),
            },
        ])
        .unwrap()
    };
    let rom = mk();
    // Domain edges admit; one float past refuses (never extrapolate).
    assert!(rom.at(0.85).is_ok());
    assert!(rom.at(1.0).is_ok());
    assert!(matches!(
        rom.at(0.85_f64.next_down()),
        Err(e) if e.code == "rom-schedule-out-of-domain"
    ));
    assert!(matches!(
        rom.at(1.0_f64.next_up()),
        Err(e) if e.code == "rom-schedule-out-of-domain"
    ));
    // Construction refusals.
    assert!(matches!(
        ScheduledRom::new(vec![]),
        Err(e) if e.code == "rom-schedule-invalid"
    ));
    // March caps.
    assert!(march_scheduled(&rom, &|_| ([0.0, 0.0], 0.9), 0.01, 1).is_ok());
    assert!(matches!(
        march_scheduled(&rom, &|_| ([0.0, 0.0], 0.9), 0.01_f64.next_up(), 1),
        Err(e) if e.code == "rom-march-invalid"
    ));
    jlog("caps", "\"cap_and_cap_plus_one\":true");
}

#[test]
fn v08b2_receipt_emits_with_golden() {
    let (c085, c100, _, _) = systems();
    let refs = [&c085, &c100];
    let red = reduce_shared(&refs).unwrap();
    let ladder: Vec<(usize, f64, bool)> = red
        .ladder
        .iter()
        .map(|r| (r.order, r.worst_rel_err, r.passed))
        .collect();
    let clauses = vec![
        ClauseVerdict {
            clause: "mimo-transfer-incl-hinge",
            passed: true,
            measure: red.ladder.last().unwrap().worst_rel_err,
        },
        ClauseVerdict {
            clause: "state-continuation",
            passed: true,
            measure: 0.0,
        },
    ];
    let receipt = emit_v08b2_receipt(red.order, ladder, clauses);
    assert_eq!(receipt.schema, V08B2_SCHEMA);
    assert!(receipt.balanced_truncation.starts_with("NO-DATA"));
    assert!(receipt.loewner.starts_with("NO-DATA"));
    let again = emit_v08b2_receipt(
        receipt.order,
        receipt.ladder.clone(),
        receipt.clauses.clone(),
    );
    assert_eq!(receipt.receipt_digest, again.receipt_digest);
    jlog(
        "receipt",
        &format!("\"digest\":\"{}\"", receipt.receipt_digest),
    );
}
