//! The reed-card row `open_hole_orifice N`: absent or `1` keeps the first-hole shunt bit for bit, `2` and `3`
//! solve the first holes together, and anything else refuses.
use super::*;

const SOURCE: &str = include_str!("../../../../examples/reed-regional-mixture.performance");

/// A 7.45 mm bore with three open holes (two sizes) behind a flanged mouth, a massless reed and a 70 ms blow.
fn card(holes: usize, row: &str) -> String {
    let schedule = SOURCE.split_once("\nschedule\n").unwrap().1;
    let mut segments = String::from("cylinder 0.00745 0.12\n");
    let hole_rows = [
        "hole 0.00135 0.0127 0.00745 1\n",
        "hole 0.00235 0.0069 0.00745 1\n",
        "hole 0.0020 0.0050 0.00745 1\n",
    ];
    for hole in hole_rows.iter().take(holes) {
        segments += hole;
        segments += "cylinder 0.00745 0.02\n";
    }
    segments += "cylinder 0.00745 0.25\n";
    let count = 2 * holes + 2;
    format!(
        "frankensim-reed-performance-v1\naudio 48000 4801 20000\nambient 293.15 100000 0\n\
         reed 0.0004 0.013 6000 0 0 0.35\nlistener_m 1\ntermination flanged\nsegments {count}\n{segments}\
         {row}compile_limits 100000 1024\nschedule\n{schedule}"
    )
}

fn render(text: &str) -> Result<(Vec<f64>, [u64; 5]), ReedPerformanceError> {
    let block = 64;
    let mut p = ReedPerformance::from_bytes(text.as_bytes(), block)?;
    let mut out = vec![0.0; p.info().samples as usize];
    for part in out.chunks_mut(block) {
        p.block(part).unwrap();
    }
    let RenderVoice::ReedBore(voice) = &p.renderer.context().voices[0] else { panic!("expected a reed voice") };
    let motion = [voice.reed_y, voice.reed_v, voice.p_plus_prev, voice.f_jet_prev, voice.line.incoming()]
        .map(f64::to_bits);
    assert!(out.iter().all(|x| x.is_finite()));
    Ok((out, motion))
}

fn bits(samples: &[f64]) -> Vec<u64> {
    samples.iter().map(|x| x.to_bits()).collect()
}

#[test]
fn an_absent_row_and_row_one_are_bit_identical() {
    let none = render(&card(3, "")).unwrap();
    let one = render(&card(3, "open_hole_orifice 1\n")).unwrap();
    assert_eq!(bits(&none.0), bits(&one.0));
    assert_eq!(none.1, one.1);
    assert!(none.0.iter().any(|x| x.abs() > 1.0e-6), "the card must sound");
}

#[test]
fn a_single_open_hole_ignores_the_count() {
    let one = render(&card(1, "open_hole_orifice 1\n")).unwrap();
    let three = render(&card(1, "open_hole_orifice 3\n")).unwrap();
    assert_eq!(bits(&one.0), bits(&three.0));
    assert_eq!(one.1, three.1);
}

#[test]
fn two_and_three_coupled_holes_change_the_sound_and_differ_from_each_other() {
    let one = render(&card(3, "")).unwrap();
    let two = render(&card(3, "open_hole_orifice 2\n")).unwrap();
    let three = render(&card(3, "open_hole_orifice 3\n")).unwrap();
    assert_ne!(bits(&one.0), bits(&two.0));
    assert_ne!(bits(&two.0), bits(&three.0));
    // The extra holes only add loss: the level may fall but never rise past the first-hole render by much.
    let peak = |x: &[f64]| x.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
    assert!(peak(&three.0) <= 1.05 * peak(&one.0), "{} vs {}", peak(&three.0), peak(&one.0));
}

#[test]
fn counts_outside_one_to_three_refuse() {
    for row in ["open_hole_orifice 0\n", "open_hole_orifice 4\n", "open_hole_orifice\n", "open_hole_orifice 2 3\n"] {
        assert!(render(&card(3, row)).is_err(), "{row:?} must refuse");
    }
}
