//! The card row `residual_flow AREA_M2` (the effective area of the channels between the reed blades that stay open on the lay,
//! Almeida, Vergez and Causse 2007, Sec. III.1 and Fig. 9), on a card whose reed is blown past its closing pressure so that it
//! is pressed shut. G0: the row reaches the voice and a card without it has none; the row refuses a bad area, a massless reed,
//! `slit_length`, `staple_recovery` and an area at or above the rest jet area. G1: a vanishing area changes the bytes by no more
//! than rounding. G2: a reed that lies shut and still across a whole step (opening and velocity 0 at both ends, after the blow
//! has become constant) passes the leak `S sqrt(2 dp / rho)` and nothing else, with and without a contraction, while the same card
//! without the row passes nothing. Synthetic numbers: a test of the row's arithmetic, not an instrument.
use super::*;

const HEAD: &str = "frankensim-reed-performance-v1\naudio 48000 4800 20000\nambient 293.15 101325 0\n";
// Blown at 9 kPa, 1.5 times the reed's 6 kPa closing pressure: the reed is pressed on its lay.
const TAIL: &str = "listener_m 1\ntermination flanged\nsegments 1\ncylinder 0.0073 0.5\n\
    compile_limits 100000 1024\nschedule\nfrankensim-gesture-schedule-v1\n\
    control_rate_hz\t700\ntracks\t1\ntrack\tblow\tblowing-pressure\n\
    initial\tpressure-pa\t0e0\nevents\t1\nevent\t0e0\t1e-2\tpressure-pa\t9e3\n";
// A massive reed with a very soft lay (stiffness 1e3, exponent 2, onset 0.2 mm, no loss).
const SOFT: &str = "reed 0.0004 0.013 6000 3.8e-6 1400 0.077 1e3 2 2e-4 0";
const MASSLESS: &str = "reed 0.0004 0.013 6000 0 0 0.35";
// The effective area of a 5 mm^2-class reed's residual channels, of the order of the printed one (0.097 mm^2 = 9.7e-8 m^2 for the
// reed of Almeida et al. 2007, Fig. 9): 5e-8 m^2 here, well under the rest jet area 5.2e-6 m^2.
const AREA: f64 = 5e-8;

fn card(reed_row: &str, rows: &str) -> String {
    format!("{HEAD}{reed_row}\n{rows}{TAIL}")
}

fn refusal(text: &str) -> &'static str {
    match ReedPerformance::from_bytes(text.as_bytes(), 64) {
        Err(ReedPerformanceError::Input { what, .. }) => what,
        other => panic!("expected an input refusal, got {:?}", other.err()),
    }
}

fn model_refusal(text: &str) -> &'static str {
    match ReedPerformance::from_bytes(text.as_bytes(), 64) {
        Err(ReedPerformanceError::Model(AcousticRealizeError::InvalidDescription { what })) => what,
        other => panic!("expected a model refusal, got {:?}", other.err()),
    }
}

/// The pressure bits of the whole window and the residual flow area the voice carries.
fn run(text: &str) -> (Vec<u64>, Option<f64>) {
    let mut performance = ReedPerformance::from_bytes(text.as_bytes(), 64).expect("the card must load");
    let mut bits = Vec::with_capacity(4800);
    let mut out = vec![0.0; 64];
    let mut area = None;
    for _ in 0..75 {
        performance.block(&mut out).expect("the short window must render");
        bits.extend(out.iter().map(|x| x.to_bits()));
        let RenderVoice::ReedBore(voice) = &performance.renderer.context().voices[0] else {
            panic!("expected the original reed voice");
        };
        area = voice.residual_flow_area_m2();
    }
    (bits, area)
}

#[test]
fn g0_the_row_reaches_the_voice_and_a_card_without_it_has_none() {
    let (_, absent) = run(&card(SOFT, ""));
    let (_, present) = run(&card(SOFT, &format!("residual_flow {AREA:e}\n")));
    assert_eq!(absent, None, "no row, no residual flow");
    assert_eq!(present, Some(AREA), "the row must reach the voice with its area");
    // With a contraction in front of it and the closure stop row between them, the order of the grammar holds.
    let (_, both) = run(&card(SOFT, &format!("contraction 0.611\nclosure_stop on\nresidual_flow {AREA:e}\n")));
    assert_eq!(both, Some(AREA));
}

#[test]
fn g0_residual_flow_refuses_a_bad_area_a_massless_reed_a_slit_a_staple_and_an_area_past_the_rest_jet() {
    assert_eq!(refusal(&card(SOFT, "residual_flow 0\n")), "residual_flow needs a finite positive area in m^2");
    assert_eq!(refusal(&card(SOFT, "residual_flow -1e-8\n")), "residual_flow needs a finite positive area in m^2");
    assert_eq!(refusal(&card(MASSLESS, "residual_flow 5e-8\n")), "residual_flow requires a massive reed");
    assert_eq!(
        refusal(&card(SOFT, "slit_length 0.0014\nresidual_flow 5e-8\n")),
        "residual_flow is defined on the inviscid reed law only, not with slit_length"
    );
    assert_eq!(
        refusal(&card(SOFT, "staple_recovery 0.8 4e-6\nresidual_flow 5e-8\n")),
        "residual_flow is not defined with staple_recovery"
    );
    // The rest jet area is C_c w h0 = 1 x 0.013 x 0.0004 = 5.2e-6 m^2: an area at or above it is refused by the voice.
    assert_eq!(
        model_refusal(&card(SOFT, "residual_flow 6e-6\n")),
        "the residual flow area must be smaller than the reed's rest jet area"
    );
    // A trailing field is refused by the row reader (its own message).
    let _ = refusal(&card(SOFT, "residual_flow 5e-8 extra\n"));
}

#[test]
fn g1_a_vanishing_area_changes_the_bytes_by_no_more_than_rounding() {
    let (plain, _) = run(&card(SOFT, ""));
    let (tiny, area) = run(&card(SOFT, "residual_flow 1e-20\n"));
    assert_eq!(area, Some(1e-20));
    let peak = plain.iter().fold(0.0_f64, |m, bits| m.max(f64::from_bits(*bits).abs()));
    assert!(peak > 0.0 && peak.is_finite(), "the plain window must sound: {peak:e}");
    let worst = plain.iter().zip(&tiny).fold(0.0_f64, |m, (a, b)| m.max((f64::from_bits(*a) - f64::from_bits(*b)).abs()));
    assert!(worst <= 1e-9 * peak, "a vanishing leak must not change the sound: {worst:e} against peak {peak:e}");
}

/// Step the card one sample at a time and collect, at every sample from `FROM` on at which the reed lay shut and still across the whole
/// step (its opening and velocity 0 at both ends, so the mid-step opening and swept flow are exactly 0 whatever iterate the implicit
/// solve took), the flow the junction returned, the pressure drop it saw and the density: `(flow, dp, rho)`. The blow is constant by
/// `FROM` (the attack ends at 10 ms, sample 480), so the blowing pressure read after the step is the one the step used.
#[allow(clippy::float_cmp)]
fn shut_samples(rows: &str) -> Vec<(f64, f64, f64)> {
    const FROM: usize = 700;
    let text = card(SOFT, rows);
    let mut performance = ReedPerformance::from_bytes(text.as_bytes(), 1).expect("the card must load");
    let mut out = vec![0.0; 1];
    let mut shut = Vec::new();
    for sample in 0..1200 {
        let (old_y, old_v) = {
            let RenderVoice::ReedBore(v) = &performance.renderer.context().voices[0] else { panic!("expected the reed voice") };
            (v.reed_y, v.reed_v)
        };
        performance.block(&mut out).expect("one sample must render");
        let RenderVoice::ReedBore(v) = &performance.renderer.context().voices[0] else { panic!("expected the reed voice") };
        if sample < FROM || old_y != 0.0 || old_v != 0.0 || v.reed_y != 0.0 || v.reed_v != 0.0 {
            continue;
        }
        let incoming = v.line.incoming();
        let flow = (v.p_plus_prev - incoming) / v.zc;
        let dp = v.reed.blowing_pressure_pa - (v.p_plus_prev + incoming);
        shut.push((flow, dp, v.rho));
    }
    shut
}

#[test]
fn g2_a_reed_lying_shut_passes_the_leak_and_nothing_else_with_or_without_a_contraction() {
    for contraction in ["", "contraction 0.611\n"] {
        let with = shut_samples(&format!("{contraction}residual_flow {AREA:e}\n"));
        let without = shut_samples(contraction);
        assert!(with.len() >= 100, "{contraction:?}: the reed must lie shut and still at many samples, found {}", with.len());
        assert!(without.len() >= 100, "{contraction:?}: the card without the row must lie shut too, found {}", without.len());
        // The effective leak area is S whatever the contraction: C_c w (S / (C_c w)) = S.
        let mut worst = 0.0_f64;
        for &(flow, dp, rho) in &with {
            let leak = AREA * (2.0 * dp.abs() / rho).sqrt() * dp.signum();
            assert!(leak.abs() > 0.0, "{contraction:?}: the shut reed must see a pressure drop");
            worst = worst.max((flow - leak).abs() / leak.abs());
        }
        assert!(worst < 1e-8, "{contraction:?}: a shut reed must pass exactly the leak S sqrt(2 dp / rho): worst relative gap {worst:e}");
        // The same card without the row passes nothing at those samples, against a leak of this size.
        let scale = AREA * (2.0 * 9.0e3 / with[0].2).sqrt();
        let leaked = without.iter().fold(0.0_f64, |m, &(flow, _, _)| m.max(flow.abs()));
        assert!(leaked < 1e-9 * scale, "{contraction:?}: without the row nothing may pass a shut reed: {leaked:e} against a leak of {scale:e}");
    }
}
