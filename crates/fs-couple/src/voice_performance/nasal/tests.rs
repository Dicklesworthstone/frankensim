//! N3 gates: (2) a shut port is the series tract, line and lip flow; (3) a nasal card's island line is passive, the
//! lips' and both nostrils' lines realize, and the render reaches the three-opening receiver and differs from the oral
//! card's; the nasal rows' refusals by name; the ostium's split of its passage section. The dimensions are test
//! geometry, not anatomy. Gate (1), a row-less card rendering md5-equal before and after, is a byte check on built
//! binaries. The tolerance of (2), 1e-9 of the largest tap, is N2's one-edge one: the two sides solve the same chain by
//! different arithmetic, so equal bits are not expected.

use super::*;
use crate::driving_point::{PASSIVITY_TOL, analyze_reflectance_ir};

const TAP_TOL: f64 = 1.0e-9;

/// The parent's test card (a short uniform-ish tract; fold values scaled from the two-mass pack), with `extra` rows.
fn fixture(extra: &str) -> String {
    let mut text = String::from(
        "frankensim-voice-performance-v1\n\
         audio 48000 2400 2\n\
         full_scale_pa 2.0\n\
         ambient moist_air 310.15 101325 1.0\n\
         fold two_mass 1.25e-4 2.5e-5 128 12.8 40 0.1 3.57e-4 1.4e-2 2.5e-3\n\
         fold_damping_upper 0.6\n\
         termination flanged\n\
         tract_source test-fixture\n\
         tract_licence test fixture, not data\n\
         microphone 0 0 0.3 4 16 4000\n\
         subglottal 0 0\n\
         subglottal 0.02 900\n",
    );
    for area in [0.5e-4, 0.6e-4, 1.5e-4, 2.5e-4, 3.0e-4, 3.0e-4, 2.5e-4, 2.0e-4] {
        text.push_str(&format!("section {area:e} 0.02\n"));
    }
    text.push_str(extra);
    text
}

/// Test geometry: the port after section 3, a three-section common passage, three pair sections and a sinus each side.
fn nasal_rows() -> String {
    let mut text = String::from(
        "nasal_source test-fixture\n\
         nasal_licence test fixture, not data\n\
         velum 3 2.0e-5 test-fixture\n",
    );
    for area in [1.0e-4, 2.0e-4, 2.5e-4] {
        text.push_str(&format!("nasal_section {area:e} 0.02\n"));
    }
    for _ in 0..3 {
        text.push_str("nasal_pair 1.2e-4 1.0e-4 0.02\n");
    }
    text.push_str("sinus left 0.03 1.0e-5 0.005 1.6e-5\n");
    text.push_str("sinus right 0.025 1.2e-5 0.006 1.7e-5\n");
    text.push_str("nostrils 0.008 0.028 test-fixture\n");
    text
}

const WALL: &str = "tract_wall 15 3e5 8000\n";

fn parse(text: &str) -> VoiceDescription {
    VoiceDescription::from_bytes(text.as_bytes()).expect("test card parses")
}

fn impulse(line: &DelayedFilter, n: usize) -> Vec<f64> {
    let mut probe = line.clone();
    (0..n).map(|k| probe.push(if k == 0 { 1.0 } else { 0.0 }).expect("push")).collect()
}

fn relative_difference(a: &[f64], b: &[f64]) -> f64 {
    assert_eq!(a.len(), b.len());
    let scale = a.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
    a.iter().zip(b).fold(0.0_f64, |m, (x, y)| m.max((x - y).abs())) / scale
}

fn render(description: &VoiceDescription) -> (Vec<f64>, VoiceInfo) {
    let mut performance = description.build().expect("the card builds");
    let mut pcm = vec![0.0; 1200];
    for block in pcm.chunks_mut(37) {
        performance.block(block).expect("physical pressure");
    }
    assert!(pcm.iter().all(|v| v.is_finite()));
    assert!(pcm.iter().any(|v| v.abs() > 1e-12));
    (pcm, performance.info().clone())
}

/// (2) With the port shut the graph is the oral tract split at the velum: the island's line and the lip flow are the
/// series duct's, on the same ladder, rigid and walled, against the same `rho c / S`.
#[test]
fn a_shut_port_is_the_series_tract() {
    let description = parse(&fixture(&nasal_rows()));
    let card = description.nasal.clone().expect("nasal rows");
    let gas = *description.tract_gas();
    let rate = description.output_rate_hz * description.mechanics_per_output;
    let duct = TractChart::try_new(description.sections.clone(), "test-fixture", "test fixture, not data")
        .expect("chart")
        .to_duct();
    // The series build's mouth line takes the default ladder (voice_performance.rs), so the reference flow does too.
    let ladder = LineLadder::default();
    let tissue = WallPin { surface_density: 15.0, stiffness_per_area: 3.0e5, resistance: 8000.0 };
    for (case, wall) in [("rigid", None), ("walled", Some(tissue))] {
        let lines = nasal_lines(
            &description.sections, &card, &gas, wall.as_ref(), Termination::FlangedOpen, rate, Velum::Shut,
        )
        .expect("shut-port lines");
        assert!(lines.nostrils.is_none());
        let (series_taps, zc) = glottal_tract_line_taps(
            &duct, &gas, None, MouthLoad::Analytic(Termination::FlangedOpen), rate, wall.as_ref(),
        )
        .expect("series line");
        assert_eq!(glottal_port_impedance(&gas, lines.inlet_radius_m).0.to_bits(), zc.to_bits());
        let (_, network_taps) = lines.line.quiescent_fir().expect("a plain FIR");
        assert_eq!(network_taps.len(), series_taps.len(), "{case}: the ladders stopped on different rungs");
        let line_error = relative_difference(&series_taps, network_taps);
        let (series_flow, series_receipt) =
            mouth_flow_line(&duct, &gas, Termination::FlangedOpen, rate, zc, wall.as_ref(), &ladder).expect("series flow");
        let (network_flow, network_receipt) = &lines.lips;
        assert_eq!(network_receipt.n_fft, series_receipt.n_fft, "{case}: the flow ladders stopped on different rungs");
        let n = series_receipt.n_fft;
        let flow_error = relative_difference(&impulse(&series_flow, n), &impulse(network_flow, n));
        eprintln!(
            "{{\"suite\":\"fs-couple\",\"case\":\"voice-nasal-shut-port-{case}\",\"line_taps\":{},\"flow_n_fft\":{n},\"line_max_relative_tap_difference\":{line_error:.3e},\"flow_max_relative_tap_difference\":{flow_error:.3e},\"tol\":{TAP_TOL:e}}}",
            series_taps.len()
        );
        assert!(line_error <= TAP_TOL, "{case}: the shut-port line differs from the series line by {line_error:.3e}");
        assert!(flow_error <= TAP_TOL, "{case}: the shut-port lip flow differs from the series flow by {flow_error:.3e}");
    }
}

/// (3) The walled nasal card: the island's line is passive on the 4x grid, the lips' and both nostrils' flow lines
/// realize, and the render reaches the three-opening receiver and differs from the same card without the nasal rows.
#[test]
fn a_nasal_card_is_passive_and_reaches_the_receiver() {
    let description = parse(&fixture(&format!("{WALL}{}", nasal_rows())));
    let card = description.nasal.clone().expect("nasal rows");
    let gas = *description.tract_gas();
    let rate = description.output_rate_hz * description.mechanics_per_output;
    let lines = nasal_lines(
        &description.sections, &card, &gas, description.wall.as_ref(), description.termination, rate, Velum::Open,
    )
    .expect("nasal lines");
    let (_, taps) = lines.line.quiescent_fir().expect("a plain FIR");
    let (peak, _) = analyze_reflectance_ir(&impulse(&lines.line, taps.len()), 1.0 / f64::from(rate), lines.round_trip_samples);
    let nostrils = lines.nostrils.as_ref().expect("an open port's nostrils");
    for (label, (line, receipt)) in [("lips", &lines.lips), ("nostril-left", &nostrils[0]), ("nostril-right", &nostrils[1])] {
        let ir = impulse(line, receipt.n_fft);
        assert!(ir.iter().all(|h| h.is_finite()) && ir.iter().any(|h| *h != 0.0), "{label}: empty flow line");
        let last = receipt.rungs.last().expect("a converged ladder has a rung");
        eprintln!(
            "{{\"suite\":\"fs-couple\",\"case\":\"voice-nasal-flow\",\"port\":\"{label}\",\"n_fft\":{},\"tail_energy\":{:.4e},\"tap_change\":{:.4e}}}",
            receipt.n_fft, last.tail_energy, last.tap_change
        );
    }
    eprintln!(
        "{{\"suite\":\"fs-couple\",\"case\":\"voice-nasal-line\",\"taps\":{},\"peak_reflection\":{peak:.12}}}",
        taps.len()
    );
    assert!(peak <= 1.0 + PASSIVITY_TOL, "the nasal card's line is active on the 4x grid: {peak}");
    let (nasal_pcm, info) = render(&description);
    assert_eq!(info.nostril_flow_taps, Some([nostrils[0].1.n_fft, nostrils[1].1.n_fft]));
    assert!(info.passive_formants_hz.is_err());
    let (oral_pcm, oral_info) = render(&parse(&fixture(WALL)));
    assert_eq!(oral_info.nostril_flow_taps, None);
    let change = relative_difference(&oral_pcm, &nasal_pcm);
    eprintln!("{{\"suite\":\"fs-couple\",\"case\":\"voice-nasal-render\",\"max_relative_change_from_oral\":{change:.3e}}}");
    assert!(change > 1.0e-3, "the nasal tract did not reach the radiated pressure: {change:.3e}");
}

/// The rows' admission and refusals, each by name.
#[test]
fn nasal_rows_are_refused_by_name() {
    let rows = nasal_rows();
    let without = |key: &str| rows.lines().filter(|l| !l.starts_with(key)).map(|l| format!("{l}\n")).collect::<String>();
    let cases = [
        (without("nostrils"), "a nasal card needs `nostrils`"),
        (without("velum"), "a nasal card needs `velum`"),
        (without("nasal_pair"), "a nasal card needs `nasal_pair`"),
        (without("nasal_licence"), "a nasal card needs `nasal_licence`"),
        (format!("{}velum 7 2.0e-5 test-fixture\n", without("velum")), "needs an oral section after it"),
        (format!("{}velum 3 0 test-fixture\n", without("velum")), "velum port area must be positive"),
        (format!("{rows}velum 3 2.0e-5 test-fixture\n"), "`velum` may only appear once"),
        (format!("{}sinus left 0.07 1.0e-5 0.005 1.6e-5\n", without("sinus left")), "sinus left opens"),
        (format!("{}sinus right 0.06 1.2e-5 0.006 1.7e-5\n", without("sinus right")), "sinus right opens"),
        (format!("{rows}sinus middle 0.01 1.0e-5 0.005 1.6e-5\n"), "sinus side must be `left` or `right`"),
        (format!("{rows}tract_move 0 0.01 card 1\n"), "refused with `tract_move`"),
        (format!("{rows}head_sphere 0.09 test-fixture\n"), "refused with `head_sphere`"),
        (format!("{rows}fricative_noise 1 2 5\n"), "refused with `fricative_noise`"),
        (format!("{rows}receiver free_field_monopole\n"), "refused with `receiver free_field_monopole`"),
    ];
    for (extra, needle) in cases {
        let error = VoiceDescription::from_bytes(fixture(&extra).as_bytes()).expect_err(needle);
        assert!(error.contains(needle), "expected {needle:?}, got {error:?}");
    }
    assert!(parse(&fixture("")).nasal.is_none());
    let card = parse(&fixture(&nasal_rows())).nasal.expect("nasal rows");
    assert_eq!((card.after_section, card.sections.len(), card.pairs.len()), (3, 3, 3));
    assert!(card.sinus_left.is_some() && card.sinus_right.is_some());
}

/// An ostium inside a passage section splits it into two of the same area; on a boundary nothing is split. The graph
/// with both sinuses: two oral edges, the passage, and three edges a side.
#[test]
fn an_ostium_splits_its_passage_section_exactly() {
    let chain = [(1.0e-4, 0.02), (2.0e-4, 0.02), (3.0e-4, 0.02)];
    let (back, front) = split_chain(&chain, 0.03);
    let areas = |c: &[(f64, f64)]| c.iter().map(|s| s.0).collect::<Vec<_>>();
    assert_eq!(areas(&back), vec![1.0e-4, 2.0e-4]);
    assert_eq!(areas(&front), vec![2.0e-4, 3.0e-4]);
    let total: f64 = back.iter().chain(&front).map(|s| s.1).sum();
    assert!((total - 0.06).abs() < 1e-15, "the split moved length: {total}");
    let (back, front) = split_chain(&chain, 0.02);
    assert_eq!((back.len(), front.len()), (1, 2));
    for at in [0.02 + 1.0e-13, 0.02 - 1.0e-13] {
        let (back, front) = split_chain(&chain, at);
        assert_eq!((back.len(), front.len()), (1, 2), "a cut {at} within rounding of a boundary left a sliver");
    }
    let description = parse(&fixture(&nasal_rows()));
    let card = description.nasal.clone().expect("nasal rows");
    let (network, terminations, radius) =
        nasal_graph(&description.sections, &card, Termination::FlangedOpen, Velum::Open);
    assert_eq!((network.nodes.len(), network.edges.len(), terminations.len()), (10, 9, 3));
    assert_eq!(radius[0].to_bits(), (1.2e-4 / core::f64::consts::PI).sqrt().to_bits());
    assert_eq!(radius[1].to_bits(), (1.0e-4 / core::f64::consts::PI).sqrt().to_bits());
}
