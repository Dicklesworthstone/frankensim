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
fn counts_outside_one_to_the_cap_refuse() {
    // The cap is 12 (`fs_duct::exterior::MAX_ORIFICE_PORTS`); a count above the open holes is clamped, not refused.
    for row in ["open_hole_orifice 0\n", "open_hole_orifice 13\n", "open_hole_orifice\n", "open_hole_orifice 2 3\n"] {
        assert!(render(&card(3, row)).is_err(), "{row:?} must refuse");
    }
    let clamped = render(&card(3, "open_hole_orifice 12\n")).unwrap();
    let three = render(&card(3, "open_hole_orifice 3\n")).unwrap();
    assert_eq!(bits(&clamped.0), bits(&three.0), "12 on a three-hole bore is the three-port render");
}

// ---- N above 1 on a `gas_field continuous` binding: each orifice reads its own local gas ----

mod on_a_field {
    use super::super::*;
    use crate::driving_point::regional::field::tests::field;

    const FIXTURE: &str = include_str!("../../../../examples/reed-flanged-baffled.performance");

    /// The field tests' bore (one 0.57 m cylinder, a 280 to 300 K field) with three open holes; the lengths sum to 0.57 m.
    fn holed() -> String {
        FIXTURE
            .replace("observer_law stokes-kirchhoff\n", "")
            .replace("audio 48000 12000 2", "audio 48000 4801 2")
            .replace("ambient mixture 293.15 100000 mole 2 N2 0.7 O2 0.3", "ambient 300 101325 0")
            .replace("reed 0.0004 0.013 6000 0 0 0.35", "reed 0.0004 0.013 6000 0.00001 100 0.35")
            .replace("segments 2\ncylinder 0.0015 0.2\ncylinder 0.0015 0.3 gas mixture 293.15 100000 mole 3 N2 0.4 O2 0.2 CO2 0.4",
                "segments 7\ncylinder 0.0015 0.285\nhole 0.0013 0.0005 0.0015 1\ncylinder 0.0015 0.0475\n\
                 hole 0.0013 0.0005 0.0015 1\ncylinder 0.0015 0.0475\nhole 0.0012 0.0005 0.0015 1\ncylinder 0.0015 0.19")
            .replace("termination flanged", "exterior 300 101325 0\ntermination flanged")
    }
    fn with_row(text: &str, row: &str) -> String {
        text.replace("compile_limits", &format!("{row}\ncompile_limits"))
    }
    fn continuous(text: &str, outlet_temperature: f64) -> Result<ReedPerformance, ReedPerformanceError> {
        ReedPerformance::from_bytes_with_field_resolution(
            text.as_bytes(), 113, ObservationPolicy::Diagnostic, None, field(outlet_temperature),
            GasFieldResolution::Continuous(fs_duct::RiccatiSteps(16)),
        )
    }
    fn pressure(text: &str) -> Vec<f64> {
        let mut p = continuous(text, 280.0).unwrap();
        let mut out = vec![0.0; p.info.samples as usize];
        for chunk in out.chunks_mut(113) {
            p.pressure_block(chunk).unwrap();
        }
        assert!(out.iter().all(|x| x.is_finite()));
        out
    }
    fn distance(a: &[f64], b: &[f64]) -> f64 {
        a.iter().zip(b).map(|(x, y)| (x - y).powi(2)).sum::<f64>().sqrt()
    }
    fn norm(a: &[f64]) -> f64 {
        a.iter().map(|x| x * x).sum::<f64>().sqrt()
    }

    #[test]
    fn an_absent_row_and_row_one_are_bit_identical_on_a_continuous_field() {
        let base = holed();
        let plain = pressure(&base);
        assert!(norm(&plain) > 0.0);
        assert_eq!(plain, pressure(&with_row(&base, "open_hole_orifice 1")));
        assert_eq!(plain, pressure(&base), "deterministic");
    }

    #[test]
    fn two_and_three_coupled_holes_act_on_a_continuous_field_and_differ_from_each_other() {
        let base = holed();
        let (one, two, three) = (
            pressure(&base),
            pressure(&with_row(&base, "open_hole_orifice 2")),
            pressure(&with_row(&base, "open_hole_orifice 3")),
        );
        let scale = norm(&one);
        assert!(distance(&one, &two) > 1.0e-9 * scale, "N = 2 must act");
        assert!(distance(&two, &three) > 1.0e-9 * scale, "the third hole must act");
        assert_eq!(three, pressure(&with_row(&base, "open_hole_orifice 3")), "deterministic");
    }

    #[test]
    fn the_orifice_ports_compose_with_the_open_hole_interaction_and_the_printed_lengths() {
        let coupled = with_row(&holed(), "open_hole_interaction tmmi 2");
        let one = pressure(&coupled);
        let three = pressure(&with_row(&coupled, "open_hole_orifice 3"));
        assert!(distance(&one, &three) > 1.0e-9 * norm(&one), "N = 3 acts under the interaction");
        let lengths = with_row(&coupled, "tonehole_lengths lefebvre-scavone-2012");
        let with_lengths = pressure(&with_row(&lengths, "open_hole_orifice 3"));
        assert!(distance(&pressure(&lengths), &with_lengths) > 1.0e-9 * norm(&one), "N = 3 acts with the printed lengths");
    }

    #[test]
    fn a_midpoint_gas_field_still_refuses_above_one_by_name() {
        let text = with_row(&holed(), "open_hole_orifice 2");
        let error = ReedPerformance::from_bytes_with_gas_field(
            text.as_bytes(), 113, ObservationPolicy::Diagnostic, None, field(280.0), 8,
        )
        .err()
        .expect("a midpoint binding with N = 2 refuses");
        assert!(format!("{error:?}").contains("midpoint gas_field"), "{error:?}");
    }
}
