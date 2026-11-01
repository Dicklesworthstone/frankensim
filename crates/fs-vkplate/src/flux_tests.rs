//! The uniform-mass skip of the mass-flux correction is output-identical: `Body::coefs` marks the
//! rows of the Laplacian of positive mass whose correction is a sum of signed zeros, and a marked
//! row of a bounded state takes the +0.0 that sum would (predict's s1 term runs a vector pass with
//! +0.0 in place of c/m, then sums only the unmarked rows). A per-thread hook rebuilds the marks as
//! no row (the full sum everywhere, the code before the skip) or as every row (a skip forced where
//! it is wrong), so the tests compare the skip against the full sum bit for bit and show a wrong
//! skip is seen.

use std::cell::Cell;
use std::collections::BTreeMap;

use crate::live::{Live, PIECES};
use crate::{kit, Body, Csr};

const MARKED: u8 = 0;
const FULL_SUM: u8 = 1;
const FORCED: u8 = 2;

thread_local! {
    static HOOK: Cell<u8> = const { Cell::new(MARKED) };
}

/// The mark `coefs` keeps for one row under this thread's hook.
pub(crate) fn hooked(marked: bool) -> bool {
    match HOOK.with(Cell::get) {
        FULL_SUM => false,
        FORCED => true,
        _ => marked,
    }
}

fn with_hook<T>(mode: u8, f: impl FnOnce() -> T) -> T {
    HOOK.with(|h| h.set(mode));
    let out = f();
    HOOK.with(|h| h.set(MARKED));
    out
}

/// The kit page's atmosphere record (kitStage builds every lane with `withGas`).
const GAS: &str = "293.15 101325 0";
const HAT_SCALES: [f64; 8] = [2.5, 3.0, 3.5, 4.0, 4.5, 6.0, 7.0, 8.0];
/// The worker's default list, which the crash and ride use.
const CYMBAL_SCALES: [f64; 8] = [2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 8.0];

/// kitStage's grids and stroke speeds for every live piece. The two pieces kitStage does not list
/// (wires off, pedal hi-hat) take the grids and speed of their kitStage sibling.
fn stage(piece: &str) -> (&'static [f64], f64) {
    match piece {
        "kick" => (&[6.0, 7.0, 8.0], 3.0),
        "snare" | "snare-off" => (&[6.0, 7.0, 8.0], 4.5),
        "tom-12" => (&[5.0, 6.0], 3.5),
        "tom-14" => (&[4.0, 5.0, 6.0], 3.5),
        "hihat" | "hihat-closed" | "hihat-open" => (&HAT_SCALES, 2.5),
        "crash" => (&CYMBAL_SCALES, 6.0),
        "ride" => (&CYMBAL_SCALES, 2.5),
        other => panic!("{other}: no kitStage grids listed for this piece"),
    }
}

#[derive(PartialEq, Debug)]
struct Trace {
    samples: Vec<u32>,
    block_force_kinetic: Vec<u64>,
    state: Vec<u64>,
    /// Rows marked for the skip, over every body, after the run.
    marked: usize,
}

/// Build the piece (which steps its settle), strike it at kitStage's speed after 2 blocks, and
/// render 40 blocks of 128 frames (about 0.1 s).
fn run(piece: &str, scale: f64, speed: f64, mode: u8) -> Result<Trace, String> {
    with_hook(mode, || {
        let mut live = Live::with_gas_record(piece, scale, GAS)?;
        let mut block = [0.0f32; 128];
        let mut samples = Vec::with_capacity(40 * 128);
        let mut block_force_kinetic = Vec::with_capacity(80);
        for k in 0..40 {
            if k == 2 {
                live.strike(0.35 * live.radius, 0.0, speed);
            }
            live.render(&mut block, 1.0);
            samples.extend(block.iter().map(|x| x.to_bits()));
            block_force_kinetic.push(live.force.to_bits());
            block_force_kinetic.push(live.kinetic.to_bits());
        }
        let state = live
            .ins
            .bodies
            .iter()
            .flat_map(|b| b.w.iter().chain(&b.wp).chain(&b.wn).chain(&b.f).map(|x| x.to_bits()))
            .collect();
        let marked = live.ins.bodies.iter().map(|b| b.flux_zero.iter().filter(|&&z| z).count()).sum();
        Ok(Trace { samples, block_force_kinetic, state, marked })
    })
}

/// A route's result, or the route's own panic.
fn joined<T>(h: std::thread::ScopedJoinHandle<'_, T>) -> T {
    h.join().unwrap_or_else(|e| std::panic::resume_unwind(e))
}

fn same_as_the_full_sum(piece: &str) {
    assert!(PIECES.contains(&piece), "{piece} is not a live piece");
    let (scales, speed) = stage(piece);
    // Each grid's two routes run on threads of their own (the hook is per thread), so the slowest
    // single build sets this test's time, not the sum of sixteen hi-hat builds.
    let routes: Vec<_> = std::thread::scope(|s| {
        let spawned: Vec<_> = scales
            .iter()
            .map(|&scale| {
                let full = s.spawn(move || run(piece, scale, speed, FULL_SUM));
                let skip = s.spawn(move || run(piece, scale, speed, MARKED));
                (scale, full, skip)
            })
            .collect();
        spawned.into_iter().map(|(scale, full, skip)| (scale, joined(full), joined(skip))).collect()
    });
    for (scale, full, skip) in routes {
        match (&full, &skip) {
            (Ok(a), Ok(b)) => {
                assert_eq!(a.marked, 0, "{piece} {scale}x: the full-sum hook left rows marked");
                eprintln!("{piece} {scale}x: {} rows marked for the skip", b.marked);
                if matches!(piece, "kick" | "snare" | "snare-off" | "tom-12" | "tom-14") {
                    // uniform heads: without marked rows this comparison would test nothing
                    assert!(b.marked > 0, "{piece} {scale}x: no row marked on a uniform-mass head");
                }
                assert!(a.samples.iter().any(|&s| f32::from_bits(s) != 0.0), "{piece} {scale}x: silent");
                assert!(a.samples == b.samples, "{piece} {scale}x: output samples differ");
                assert!(a.block_force_kinetic == b.block_force_kinetic, "{piece} {scale}x: force or kinetic energy differs");
                assert!(a.state == b.state, "{piece} {scale}x: body state differs");
            }
            (Err(a), Err(b)) => {
                assert_eq!(a, b, "{piece} {scale}x: refusals differ");
                eprintln!("{piece} {scale}x: refused alike: {a}");
            }
            _ => panic!("{piece} {scale}x: one route refused: full sum {:?}, skip {:?}", full.as_ref().err(), skip.as_ref().err()),
        }
    }
}

#[test]
fn kick_skip_matches_the_full_sum() { same_as_the_full_sum("kick") }
#[test]
fn snare_skip_matches_the_full_sum() { same_as_the_full_sum("snare") }
#[test]
fn snare_off_skip_matches_the_full_sum() { same_as_the_full_sum("snare-off") }
#[test]
fn tom_12_skip_matches_the_full_sum() { same_as_the_full_sum("tom-12") }
#[test]
fn tom_14_skip_matches_the_full_sum() { same_as_the_full_sum("tom-14") }
#[test]
fn hihat_skip_matches_the_full_sum() { same_as_the_full_sum("hihat") }
#[test]
fn hihat_closed_skip_matches_the_full_sum() { same_as_the_full_sum("hihat-closed") }
#[test]
fn hihat_open_skip_matches_the_full_sum() { same_as_the_full_sum("hihat-open") }
#[test]
fn crash_skip_matches_the_full_sum() { same_as_the_full_sum("crash") }
#[test]
fn ride_skip_matches_the_full_sum() { same_as_the_full_sum("ride") }

#[test]
fn every_live_piece_has_a_test_above() {
    let tested = ["kick", "snare", "snare-off", "tom-12", "tom-14", "hihat", "hihat-closed", "hihat-open", "crash", "ride"];
    assert_eq!(PIECES.len(), tested.len());
    for p in PIECES {
        assert!(tested.contains(&p), "{p} has no skip test");
    }
}

/// A smooth state with distinct velocity at every node, so every edge's (x_j - x_i) is nonzero.
fn set_state(b: &mut Body) {
    for i in 0..b.n {
        b.w[i] = 1e-6 * (0.37 * i as f64).sin();
        b.wp[i] = 0.8e-6 * (0.37 * i as f64 + 0.9).sin();
    }
}

/// One predict (and finish) of a fresh clone under `mode`: the predicted displacement and the
/// loss budget, as bits.
fn step_bits(b: &Body, k: f64, mode: u8) -> (Vec<u64>, u64, Vec<bool>) {
    with_hook(mode, || {
        let mut b = b.clone();
        b.kc = 0.0;
        b.cf.clear();
        b.sync();
        b.linear_forces();
        b.predict(k);
        b.finish(k);
        // predict's s1 term sums exactly the unmarked rows
        let unmarked: Vec<usize> = (0..b.n).filter(|&i| !b.flux_zero[i]).collect();
        assert_eq!(b.flux_rows, unmarked, "{}: the summed rows are not the unmarked rows", b.name);
        (b.wn.iter().map(|x| x.to_bits()).collect(), b.budget.to_bits(), b.flux_zero.clone())
    })
}

/// The mark takes equal, positive, finite masses only. A positive mass makes the skipped division
/// exact: +0.0 / m is +0.0 for m > 0, but -0.0 for m < 0, and NaN for m = 0 or NaN.
#[test]
fn the_mark_needs_equal_positive_finite_masses() {
    let lap = Csr::from_rows(&[
        BTreeMap::from([(0, 1.0), (1, -1.0)]),
        BTreeMap::from([(0, -1.0), (1, 1.0)]),
    ]);
    for m in [1.0, 1e-12, f64::MIN_POSITIVE, f64::MAX] {
        assert!(lap.flux_zero_row(0, &[m, m]) && lap.flux_zero_row(1, &[m, m]), "{m} not marked");
    }
    for m in [0.0, -0.0, -1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(!lap.flux_zero_row(0, &[m, m]) && !lap.flux_zero_row(1, &[m, m]), "{m} marked");
    }
    assert!(!lap.flux_zero_row(0, &[1.0, 1.0 + f64::EPSILON]));
    assert!(!lap.flux_zero_row(1, &[1.0 + f64::EPSILON, 1.0]));
    assert!((0.0f64 / -1.0).is_sign_negative());
}

/// The planted negative at predict: a snare batter head with one node's mass raised. The rows that
/// touch that node lose their mark and are exactly the rows predict sums (`step_bits` checks the
/// list), so the skip equals the full sum; and forcing the skip onto every row changes the
/// predicted head, so a wrong mark would not pass unseen.
#[test]
fn a_planted_mass_step_unmarks_its_rows_and_a_forced_skip_changes_the_head() {
    let mut spec = kit::snare();
    spec.grid *= 8.0;
    let piece = kit::drum(&spec, &[]).unwrap();
    let k = 1.0 / (48_000.0 * piece.oversample as f64);
    let mut head = piece
        .instrument
        .bodies
        .into_iter()
        .find(|b| b.s1 > 0.0 && b.n > 50)
        .expect("a head with s1 loss");
    set_state(&mut head);
    // the preset head is uniform: every row marked, and forcing the skip changes nothing
    let (full, _, marks) = step_bits(&head, k, FULL_SUM);
    let (skip, _, marks_skip) = step_bits(&head, k, MARKED);
    let (forced, _, _) = step_bits(&head, k, FORCED);
    assert!(marks.iter().all(|&z| !z));
    assert!(marks_skip.iter().all(|&z| z), "{}: the preset head has unmarked rows", head.name);
    assert!(skip == full && forced == full);

    let q = (0..head.n).max_by_key(|&i| head.lap.ptr[i + 1] - head.lap.ptr[i]).unwrap();
    head.mass[q] *= 1.5;
    let touches_q: Vec<bool> = (0..head.n)
        .map(|i| (head.lap.ptr[i]..head.lap.ptr[i + 1]).any(|p| head.lap.col[p] == q && i != q) || i == q)
        .collect();
    let (full, _, _) = step_bits(&head, k, FULL_SUM);
    let (skip, _, marks) = step_bits(&head, k, MARKED);
    let (forced, _, _) = step_bits(&head, k, FORCED);
    assert!(touches_q.iter().filter(|&&t| t).count() > 1);
    for i in 0..head.n {
        assert_eq!(marks[i], !touches_q[i], "{}: row {i} (planted node {q})", head.name);
    }
    assert!(skip == full, "{}: the skip changed a head with a planted mass step", head.name);
    // the air mass and the nonlinearity couple every node, so only the planted node is checked
    assert!(forced[q] != full[q], "{}: a forced skip left the planted node {q} unchanged", head.name);
}

/// The planted negative at finish: the von Karman loss budget on a tapered cymbal, whose rows of
/// unequal mass stay unmarked and still sum, while a skip forced onto every row changes it.
#[test]
fn a_tapered_shell_keeps_its_flux_in_the_loss_budget() {
    let mut spec = kit::splash();
    spec.grid = 0.025;
    let mut b = kit::cymbal(&spec, &[]).unwrap().instrument.bodies.remove(0);
    assert!(b.vk.is_some(), "{}: no von Karman budget to test", b.name);
    let k = 1.0 / 48_000.0;
    b.s1 = 0.02;
    set_state(&mut b);
    let (full_wn, full, _) = step_bits(&b, k, FULL_SUM);
    let (skip_wn, skip, marks) = step_bits(&b, k, MARKED);
    let (forced_wn, forced, _) = step_bits(&b, k, FORCED);
    let unmarked = marks.iter().filter(|&&z| !z).count();
    eprintln!("{}: {unmarked} of {} rows unmarked", b.name, b.n);
    assert!(unmarked > 0, "{}: a tapered shell with every row marked", b.name);
    assert!(skip == full && skip_wn == full_wn, "{}: the skip changed the shell", b.name);
    assert!(f64::from_bits(full).is_finite() && f64::from_bits(full) > 0.0, "{}: budget {}", b.name, f64::from_bits(full));
    assert!(forced != full, "{}: a forced skip left the loss budget unchanged", b.name);
    assert!(forced_wn != full_wn, "{}: a forced skip left the predicted shell unchanged", b.name);
}

/// Outside the state bound every row sums, so NaN reaches the head through the correction exactly
/// as it did before the skip.
#[test]
fn a_non_finite_state_takes_the_full_sum() {
    let mut spec = kit::snare();
    spec.grid *= 8.0;
    let piece = kit::drum(&spec, &[]).unwrap();
    let k = 1.0 / (48_000.0 * piece.oversample as f64);
    let mut head = piece.instrument.bodies.into_iter().find(|b| b.s1 > 0.0 && b.n > 50).unwrap();
    set_state(&mut head);
    for bad in [f64::NAN, f64::INFINITY, f64::MAX] {
        let mut b = head.clone();
        b.wp[b.n / 2] = bad;
        assert!(!b.flux_state_bounded());
        let (full, full_budget, _) = step_bits(&b, k, FULL_SUM);
        let (skip, skip_budget, _) = step_bits(&b, k, MARKED);
        assert!(skip == full && skip_budget == full_budget, "{bad}: the skip changed a non-finite state");
    }
    assert!(head.flux_state_bounded());
}
