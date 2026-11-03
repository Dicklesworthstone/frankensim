//! The flute D4 card of the posted flute D4,
//! `data/instruments/descriptions/flute/flute-d4-onset-noise-contraction.performance`, is bound to
//! what it is made from: the committed `flute-d4.performance` with two opt-in rows added, `mouth_loss
//! contraction` after its windway row and `jet_noise` after its `jet_limits` row, and its authored
//! 50 ms onset replaced by the player's measured D4 attack, one linear segment per pixel column of
//! de la Cuadra's thesis Fig. 3.4 (`players/flute-d4-onset-delacuadra-2005-fig3-4-reads.tsv`), before
//! the D4 card's own release. The card must parse, and its gesture must be those segments exactly.
//!
//! A data-form check and a parse, never a physical claim. An edit to the card, the reads table or the
//! D4 card that is not carried to the others fails here, so the posted render stays reproducible from
//! the committed files.

use fs_couple::flue_loop::network::performance::file::FlueFilePerformance;
use fs_exec::CancelGate;
use fs_scenario::gesture::GestureValue;

const D4: &str = include_str!("../../../data/instruments/descriptions/flute/flute-d4.performance");
const CARD: &str = include_str!(
    "../../../data/instruments/descriptions/flute/flute-d4-onset-noise-contraction.performance"
);
const READS: &str = include_str!(
    "../../../data/instruments/descriptions/players/flute-d4-onset-delacuadra-2005-fig3-4-reads.tsv"
);
const SCHEMA: &str = "frankensim.player-flute-onset-figure-reads.v1";
/// The D4 card's onset, replaced, and its event count.
const D4_ONSET: &str = "event\t0e0\t5e-2\tpressure-pa\t3.32e2";
const D4_EVENTS: &str = "events\t2";
/// The row after which `mouth_loss` is added.
const WINDWAY: &str = "windway none";
const MOUTH_LOSS: &str = "mouth_loss contraction 0.6 0.3 \
    Auvray2014-TableI-alpha_vc(Zen-On-recorder-mouth;flute-embouchure-not-measured)";
const JET_NOISE: &str = "jet_noise 14 1.3e-3 0.045 0.115 1700";
/// D4's held mouth pressure [Pa], the median over its sustained segment in the per-note reads.
const HOLD_PA: f64 = 332.0;
/// One pixel of the P_m panel [Pa].
const PIXEL_PA: f64 = 8.1;
/// The first column with ink above 0 Pa, at card time 0.
const FIRST_COLUMN: u32 = 102;
/// One pixel column of the figure [s].
const COLUMN_S: f64 = 1.0 / 43.0;

struct OnsetRead {
    column: u32,
    median_pa: f64,
    time_s: f64,
    ramp_s: f64,
    target_pa: f64,
}

fn reads() -> Vec<OnsetRead> {
    let mut lines = READS.lines();
    assert_eq!(lines.next(), Some(SCHEMA), "reads table schema");
    let head: Vec<&str> = lines
        .next()
        .expect("reads table header")
        .split('\t')
        .collect();
    let col = |name: &str| {
        head.iter()
            .position(|h| *h == name)
            .unwrap_or_else(|| panic!("reads table has no column {name}"))
    };
    let (x, med, t, ramp, target) = (
        col("column_x"),
        col("pm_med_pa"),
        col("card_t_s"),
        col("card_ramp_s"),
        col("card_target_pa"),
    );
    lines
        .map(|line| {
            let f: Vec<&str> = line.split('\t').collect();
            assert_eq!(f.len(), head.len(), "reads row: {line}");
            let num = |i: usize| {
                f[i].parse::<f64>()
                    .unwrap_or_else(|e| panic!("reads field {}: {e}", f[i]))
            };
            OnsetRead {
                column: f[x].parse().expect("column_x"),
                median_pa: num(med),
                time_s: num(t),
                ramp_s: num(ramp),
                target_pa: num(target),
            }
        })
        .collect()
}

/// Rounded to 1 ns, as the card's times are.
fn ns(s: f64) -> f64 {
    (s * 1e9).round() / 1e9
}

fn keyword(line: &str) -> Option<&str> {
    line.split_whitespace().next()
}

#[test]
fn onset_card_is_the_d4_card_with_two_opt_in_rows_and_the_fig_3_4_onset() {
    let reads = reads();
    assert_eq!(reads.len(), 12, "D4's attack spans columns 102 to 113");
    for (k, r) in reads.iter().enumerate() {
        assert_eq!(r.column, FIRST_COLUMN + k as u32, "columns are consecutive");
        assert_eq!(
            r.time_s,
            ns(k as f64 * COLUMN_S),
            "column {} starts its segment",
            r.column
        );
        assert_eq!(
            r.ramp_s,
            ns(COLUMN_S),
            "column {} ramps over one column",
            r.column
        );
    }
    let (last, onset) = reads.split_last().expect("reads");
    for r in onset {
        assert_eq!(
            r.target_pa, r.median_pa,
            "column {} ends on its median",
            r.column
        );
    }
    assert_eq!(
        last.target_pa, HOLD_PA,
        "the last segment ends on D4's hold"
    );
    assert!(
        (last.median_pa - HOLD_PA).abs() <= PIXEL_PA,
        "column {} reads {} Pa, more than one pixel from the hold",
        last.column,
        last.median_pa
    );

    for anchor in [D4_ONSET, D4_EVENTS, WINDWAY] {
        let n = D4.lines().filter(|l| *l == anchor).count();
        assert_eq!(
            n, 1,
            "the D4 card must hold `{anchor}` once; update this test"
        );
    }
    assert!(
        D4_ONSET.ends_with(&format!("\t{HOLD_PA:e}")),
        "the D4 card's hold moved"
    );
    for key in ["jet_limits", "mouth_loss", "jet_noise"] {
        let n = D4.lines().filter(|l| keyword(l) == Some(key)).count();
        assert_eq!(
            n,
            usize::from(key == "jet_limits"),
            "D4 card rows named `{key}`"
        );
    }

    let events: Vec<String> = reads
        .iter()
        .map(|r| {
            format!(
                "event\t{:e}\t{:e}\tpressure-pa\t{:e}",
                r.time_s, r.ramp_s, r.target_pa
            )
        })
        .collect();
    let mut expected: Vec<String> = Vec::new();
    for line in D4.lines() {
        match line {
            D4_EVENTS => expected.push(format!("events\t{}", events.len() + 1)),
            D4_ONSET => expected.extend(events.iter().cloned()),
            _ => {
                expected.push(line.to_owned());
                if line == WINDWAY {
                    expected.push(MOUTH_LOSS.to_owned());
                } else if keyword(line) == Some("jet_limits") {
                    expected.push(JET_NOISE.to_owned());
                }
            }
        }
    }
    let got: Vec<&str> = CARD.lines().collect();
    assert_eq!(
        got, expected,
        "the card must be flute-d4.performance with the two opt-in rows and the reads' onset"
    );
    assert!(CARD.ends_with('\n'), "the card must end with a newline");
}

#[test]
fn onset_card_parses_and_its_gesture_is_the_onset_then_the_d4_release() {
    let card = FlueFilePerformance::from_bytes(CARD.as_bytes(), 256, &CancelGate::new_clock_free())
        .unwrap_or_else(|e| panic!("the onset card must parse: {e:?}"));
    let tracks = card.renderer().schedule().tracks();
    assert_eq!(tracks.len(), 1, "one gesture track");
    let blow = &tracks[0];
    assert_eq!(blow.id, "blow");
    assert_eq!(blow.initial, GestureValue::PressurePa(0.0));
    let reads = reads();
    assert_eq!(
        blow.events.len(),
        reads.len() + 1,
        "the onset segments and the release"
    );
    for (e, r) in blow.events.iter().zip(&reads) {
        assert_eq!(
            e.time_s.to_bits(),
            r.time_s.to_bits(),
            "column {} time",
            r.column
        );
        assert_eq!(
            e.transition_s.to_bits(),
            r.ramp_s.to_bits(),
            "column {} ramp",
            r.column
        );
        assert_eq!(
            e.value,
            GestureValue::PressurePa(r.target_pa),
            "column {} target",
            r.column
        );
    }
    // The reads table rounds each time and ramp to 1 ns, so a segment's end and the next
    // start can differ by up to three half-nanoseconds.
    for pair in blow.events[..reads.len()].windows(2) {
        let gap = pair[0].time_s + pair[0].transition_s - pair[1].time_s;
        assert!(
            gap.abs() <= 1.5e-9,
            "segments at {} s leave a gap of {gap} s",
            pair[0].time_s
        );
    }
    let release = blow.events.last().expect("release");
    assert_eq!(
        (release.time_s, release.transition_s),
        (2.2, 1.05e-2),
        "the D4 release"
    );
    assert_eq!(release.value, GestureValue::PressurePa(0.0));
    let end = &blow.events[reads.len() - 1];
    println!(
        "onset card: {} segments, 0 Pa at {} s to {} Pa at {} s, released at {} s",
        reads.len(),
        blow.events[0].time_s,
        HOLD_PA,
        end.time_s + end.transition_s,
        release.time_s
    );
}
