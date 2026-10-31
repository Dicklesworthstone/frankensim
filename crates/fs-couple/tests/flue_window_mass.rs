//! `window wall-cut` on the Clare D whistle's D5 card: the window derived from
//! the card's own channel, pocket and bore, its computed mass replacing the
//! legacy `0.61 a_eq` inertance, parallel cheeks and a labium ramp from the
//! row, the pocket radius tolerance, the refusals, and the window's radiation
//! from the head's outer cylinder in place of the authored loss; and
//! `window embouchure` on the G4 flute card: the hole under the lip, centred
//! or from the player's side, with the player's face behind it. Then the
//! recorder cards' window at Ernoult and Fabre 2017's mass, moved from their
//! edge plane to the card's junction, and their side holes' inner length:
//! a measured tonehole run end to end through the card's own hole code.

use fs_bem::wall_window::{
    CoverSide, WallWindow, WallWindowOptions, WindowCheeks, WindowCover, WindowRamp,
    WindowUpstream, wall_window_ladder,
};
use fs_couple::flue_loop::network::performance::file::{
    ErnoultFabreWindow, FlueFilePerformance, FlueInputError, HoleInnerLength,
};
use fs_couple::pcm_wav::observation::PressureRenderer;
use fs_exec::CancelGate;

const CARD: &str =
    include_str!("../../../data/instruments/descriptions/tin-whistle/tin-whistle-d5.performance");

/// The card's flue row: the window's radiation replaces the authored loss.
const FLUE: &str = "flue 0.00800 0.00148 0.00530 0 0 0\n";

/// The card's forced jet. These tests are about the window, so they keep the
/// relay jet; the forced jet's kernel also takes minutes to integrate per parse.
const FORCED: &str = "jet_response forced-bickley 0.1 90\n";

/// The card without its window radiation and window rows (lossless mouth) and
/// without its forced jet, cut to a 1024-sample window whose only gesture event
/// is the opening pressure ramp.
fn legacy() -> String {
    let rows: Vec<&str> = CARD
        .lines()
        .filter(|line| line.starts_with("window "))
        .collect();
    assert_eq!(rows.len(), 1, "the D5 card carries one window row");
    let radiating = format!("mouth_reaction window-cylinder\n{}\n", rows[0]);
    assert!(
        CARD.contains(&radiating) && CARD.contains(FLUE),
        "the D5 window no longer radiates through its row; update this test"
    );
    assert!(CARD.contains(FORCED), "the D5 jet changed; update this test");
    let text = CARD
        .replace(&radiating, "")
        .replace(FORCED, "")
        .replace("output 120000 5", "output 1024 1");
    let events = "events\t2\nevent\t0e0\t1.25e-2\tpressure-pa\t9.3447e1\nevent\t2.2e0\t1.05e-2\tpressure-pa\t0e0";
    assert!(
        text.contains(events),
        "the D5 gesture changed; update this test"
    );
    text.replace(events, "events\t1\nevent\t0e0\t5e-2\tpressure-pa\t9.3447e1")
}

fn with_window(text: &str, row: &str) -> String {
    assert!(text.contains("\ncontrols 4096"));
    text.replace("\ncontrols 4096", &format!("\n{row}\ncontrols 4096"))
}

fn parse(text: &str) -> Result<FlueFilePerformance, FlueInputError> {
    FlueFilePerformance::from_bytes(text.as_bytes(), 256, &CancelGate::new_clock_free())
}

fn render(file: FlueFilePerformance) -> Vec<f64> {
    let mut pressure = vec![0.0; file.info().samples as usize];
    let mut renderer = file.into_renderer();
    for part in pressure.chunks_mut(256) {
        renderer.block(part).unwrap();
    }
    assert!(pressure.iter().all(|p| p.is_finite()));
    pressure
}

#[test]
fn wall_cut_window_replaces_only_the_legacy_mouth_inertance() {
    let base = legacy();
    let old = parse(&base).unwrap();
    let new = parse(&with_window(&base, "window wall-cut 16")).unwrap();
    let record = new.window().expect("window record");
    let w = record.window;
    assert_eq!(
        (w.bore_radius, w.wall_thickness, w.arc_width, w.axial_length),
        (0.0061, 0.0024, 0.008, 0.0053),
        "geometry comes from the card's bore, channel and flue rows"
    );
    assert_eq!(
        w.upstream,
        fs_bem::wall_window::WindowUpstream::ClosedEnd { distance: 0.01065 }
    );
    assert_eq!(record.cells, [16, 24, 16]);
    assert_eq!((w.cheeks, w.ramp), (WindowCheeks::Radial, None));
    let bore_area = core::f64::consts::PI * 0.0061 * 0.0061;
    let length = record.ladder.mass_per_density * bore_area;
    // fs-bem's ignored three-rung ladder extrapolates to 15.71259 mm.
    assert!((length / 0.01571259 - 1.0).abs() < 1e-4, "{length}");
    // The bits the landed prism solver wrote into the whistle renders'
    // provenance (r02, rungs and bore length): the prism path is unchanged.
    let bits: Vec<u64> = record
        .ladder
        .rungs
        .iter()
        .map(|r| r.bore_length.to_bits())
        .collect();
    assert_eq!(
        bits,
        [
            1.5970602151823502e-2f64.to_bits(),
            1.581485190768185e-2f64.to_bits()
        ]
    );
    assert_eq!(length.to_bits(), 1.571237399621776e-2f64.to_bits());
    // The graph keeps the channel only as a one-sample inlet line; the card's
    // 2.4 mm channel stays the window's wall (asserted above).
    let c = new.renderer().voice().gas().sound_speed;
    assert_eq!(
        record.channel_line_m.to_bits(),
        (c * (1.0 / 288_000.0)).to_bits()
    );
    let riser_area = core::f64::consts::PI * 0.0036737 * 0.0036737;
    assert_eq!(
        record.riser_mass_per_density.to_bits(),
        (record.channel_line_m / riser_area).to_bits()
    );
    let inlet = |file: &FlueFilePerformance| {
        file.renderer()
            .voice()
            .network()
            .expect("physical graph")
            .source()
            .graph
            .sections[0]
            .length_m
    };
    assert_eq!(inlet(&new).to_bits(), record.channel_line_m.to_bits());
    assert_eq!(
        inlet(&old),
        0.0024,
        "without the row the channel stays a line"
    );

    let rho = new.renderer().voice().gas().density;
    let s_m = 0.008 * 0.0053;
    let legacy_inertance = rho * 0.61 * (s_m / core::f64::consts::PI).sqrt() / s_m;
    let old_inertance = old.renderer().voice().mouth_inertance_pa_s2_m3();
    assert!((old_inertance / legacy_inertance - 1.0).abs() < 1e-12);
    let inertance = new.renderer().voice().mouth_inertance_pa_s2_m3();
    assert_eq!(
        inertance.to_bits(),
        (rho * (record.ladder.mass_per_density - record.riser_mass_per_density)).to_bits()
    );
    eprintln!(
        "window branch {:.4} mm of bore; mouth inertance {inertance:.1} vs legacy {old_inertance:.1} Pa s^2/m^3",
        length * 1e3
    );
    assert!(inertance > old_inertance);
    assert_eq!(old.geometry(), new.geometry(), "the flue row is unchanged");
    assert!(old.window().is_none());
    assert_ne!(
        render(old),
        render(new),
        "the mouth mass enters the dynamics"
    );
    let none = parse(&with_window(&base, "window none")).unwrap();
    assert!(none.window().is_none());
    assert_eq!(
        none.renderer().voice().mouth_inertance_pa_s2_m3().to_bits(),
        old_inertance.to_bits()
    );
}

#[test]
fn inconsistent_windows_refuse_before_output() {
    let base = legacy();
    let refused = |text: &str| parse(text).is_err();
    assert!(refused(&with_window(&base, "window wall-cut 1")));
    assert!(refused(&with_window(&base, "window wall-cut 64")));
    assert!(refused(&with_window(&base, "window slot 16")));
    // The inlet section must be the channel through the wall, of the
    // window's area.
    assert!(refused(&with_window(
        &base.replace("flue 0.00800", "flue 0.00700"),
        "window wall-cut 16"
    )));
    // A channel shorter than one transit sample (1.19 mm at 288 kHz) has no
    // line to keep.
    assert!(refused(&with_window(
        &base.replace("duct_section 0 1 0.0024000", "duct_section 0 1 0.0011000"),
        "window wall-cut 16"
    )));
    // A stepped pocket is not represented.
    assert!(refused(&with_window(
        &base.replace(
            "duct_section 1 2 0.0106500 0.0061000",
            "duct_section 1 2 0.0106500 0.0063000"
        ),
        "window wall-cut 16"
    )));
    // A compact mouth reaction brings its own outer mass.
    assert!(refused(&with_window(
        &base.replace(
            "\ncontrols 4096",
            "\nmouth_reaction compact-baffled\ncontrols 4096"
        ),
        "window wall-cut 16"
    )));
}

#[test]
fn window_cylinder_radiates_the_solved_outer_face_in_place_of_the_authored_loss() {
    let lossless = legacy();
    let reaction = |text: &str, kind: &str| {
        with_window(text, &format!("mouth_reaction {kind}\nwindow wall-cut 8 cheeks parallel"))
    };
    // The card's own rows.
    let file = parse(&reaction(&lossless, "window-cylinder")).unwrap();
    let record = file.window().expect("window record");
    let voice = file.renderer().voice();
    let load = voice.window_radiation_load().expect("window radiation");
    let face = &record.ladder.rungs[record.ladder.rungs.len() - 1].outer_cells;
    assert_eq!(load.cells(), face.len());
    assert_eq!(load.outer_radius_m(), 0.0061 + 0.0024);
    assert!(load.band_error() <= 0.05);
    let rho = voice.gas().density;
    assert_eq!(
        voice.mouth_inertance_pa_s2_m3().to_bits(),
        (rho * (record.ladder.mass_per_density - record.riser_mass_per_density)
            - load.load().static_inertance_pa_s2_m3())
        .to_bits(),
        "the radiation's static mass comes out of the computed window mass"
    );
    let exact = load.exact_resistance(face, 1000.0).unwrap();
    assert!((load.realized_resistance(1000.0, false) / exact - 1.0).abs() <= 0.05);
    let bracket = parse(&reaction(&lossless, "window-cylinder-piston")).unwrap();
    let piston = bracket.renderer().voice().window_radiation_load().expect("bracket");
    assert_eq!(piston.cells(), 1);
    for (name, l) in [("outer face", load), ("cut-up piston", piston)] {
        eprintln!(
            "{name}: {} cells, {} R-L-C cells, M0 {:.3} Pa s^2/m^3 ({:.4} rho/sqrt(S), tail share {:.4}), realized {:.3}; band {:.4}, discrete {:.4}, tail {:.4} to {:.0} Hz; kA {:?}",
            l.cells(),
            l.load().cells().len(),
            l.static_mass_pa_s2_m3(),
            l.static_mass_pa_s2_m3() * l.open_area_m2().sqrt() / rho,
            l.static_mass_tail_share(),
            l.load().static_inertance_pa_s2_m3(),
            l.band_error(),
            l.discrete_band_error(),
            l.tail_error(),
            l.tail_top_hz(),
            l.ka_band()
        );
    }
    let authored = lossless.replace(FLUE, "flue 0.00800 0.00148 0.00530 0 0 0.02\n");
    assert!(
        parse(&reaction(&authored, "window-cylinder")).is_err(),
        "an authored loss must be declared zero"
    );
    assert!(
        parse(&with_window(&lossless, "mouth_reaction window-cylinder")).is_err(),
        "the radiation needs the window row"
    );
    render(file);
}

#[test]
fn parallel_cheeks_and_a_ramp_come_from_the_row() {
    let base = legacy();
    let file = parse(&with_window(
        &base,
        "window wall-cut 8 cheeks parallel ramp 0.0005 25",
    ))
    .unwrap();
    let record = file.window().expect("window record");
    let w = record.window;
    assert_eq!(w.cheeks, WindowCheeks::Parallel);
    // The cheeks stand the flue width apart and meet the bore on its arc.
    let a = 0.0061;
    assert_eq!(
        w.arc_width.to_bits(),
        (2.0 * a * fs_math::det::asin(0.5 * 0.008 / a)).to_bits()
    );
    assert_eq!(
        w.ramp,
        Some(WindowRamp {
            edge_height: 0.0005,
            angle: 25f64.to_radians(),
        })
    );
    let direct = wall_window_ladder(&w, &WallWindowOptions::rung(&w, 8)).expect("direct ladder");
    assert_eq!(record.ladder.rungs.len(), 3);
    assert_eq!(
        record.ladder.mass_per_density.to_bits(),
        direct.mass_per_density.to_bits()
    );
    let bore_area = core::f64::consts::PI * a * a;
    eprintln!(
        "parallel cheeks, 0.5 mm edge, 25 degree ramp: rungs {:?} mm, order {:.3}, {:.4} mm of bore",
        record
            .ladder
            .rungs
            .iter()
            .map(|r| r.bore_length * 1e3)
            .collect::<Vec<_>>(),
        record.ladder.order,
        record.ladder.mass_per_density * bore_area * 1e3
    );
    let rho = file.renderer().voice().gas().density;
    assert_eq!(
        file.renderer().voice().mouth_inertance_pa_s2_m3().to_bits(),
        (rho * (record.ladder.mass_per_density - record.riser_mass_per_density)).to_bits()
    );
}

#[test]
fn a_pocket_within_two_and_a_half_percent_of_the_bore_is_one_bore() {
    // 2 %, about the flute's cork cavity against its first head-joint section.
    let base = legacy();
    let file = parse(&with_window(
        &base.replace(
            "duct_section 1 2 0.0106500 0.0061000",
            "duct_section 1 2 0.0106500 0.0062220",
        ),
        "window wall-cut 16",
    ))
    .unwrap();
    let record = file.window().expect("window record");
    assert_eq!(record.window.bore_radius.to_bits(), 0.0061f64.to_bits());
    let bits: Vec<u64> = record
        .ladder
        .rungs
        .iter()
        .map(|r| r.bore_length.to_bits())
        .collect();
    assert_eq!(
        bits,
        [
            1.5970602151823502e-2f64.to_bits(),
            1.581485190768185e-2f64.to_bits()
        ],
        "the window solve sees the bore section's radius"
    );
}

#[test]
fn malformed_cut_rows_refuse() {
    let base = legacy();
    let refused = |row: &str| parse(&with_window(&base, row)).is_err();
    assert!(refused("window wall-cut 8 cheeks slanted"));
    assert!(refused("window wall-cut 8 cheeks"));
    assert!(refused("window wall-cut 8 ramp 0.0005"));
    assert!(refused("window wall-cut 8 ramp 0.0005 90"));
    assert!(refused("window wall-cut 8 ramp -0.0001 25"));
    // The edge must lie inside the 2.4 mm wall.
    assert!(refused("window wall-cut 8 ramp 0.003 25"));
    // Cheeks come before the ramp.
    assert!(refused("window wall-cut 8 ramp 0.0005 25 cheeks parallel"));
}

const FLUTE: &str =
    include_str!("../../../data/instruments/descriptions/flute/flute-g4.performance");

/// The G4 flute card cut to a 1024-sample window without its release, and
/// its `window` row's hole (across, along) and open strip [m].
fn flute() -> (String, [f64; 3], String) {
    let rows: Vec<&str> = FLUTE
        .lines()
        .filter(|line| line.starts_with("window embouchure "))
        .collect();
    assert_eq!(rows.len(), 1, "the G4 card carries one embouchure row");
    let fields: Vec<f64> = rows[0]
        .split_whitespace()
        .skip(3)
        .take(3)
        .map(|f| f.parse().unwrap())
        .collect();
    assert_eq!(fields.len(), 3);
    let events = "events\t2\nevent\t0e0\t5e-2\tpressure-pa\t3.79e2\nevent\t2.2e0\t1.05e-2\tpressure-pa\t0e0";
    assert!(FLUTE.contains(events), "the G4 gesture changed; update this test");
    (
        FLUTE
            .replace("output 120000 5", "output 1024 1")
            .replace(events, "events\t1\nevent\t0e0\t5e-2\tpressure-pa\t3.79e2"),
        [fields[0], fields[1], fields[2]],
        rows[0].to_owned(),
    )
}

#[test]
fn the_flute_window_is_its_embouchure_hole_under_the_lip() {
    let (text, [across, along, open], _) = flute();
    let file = parse(&text).unwrap();
    let record = file.window().expect("window record");
    let w = record.window;
    // The bore beside the hole, the 5 mm riser as the wall, the 17 mm cork
    // cavity as the pocket, parallel walls `across` apart, and the lip's
    // strip.
    let a = 0.0087744;
    assert_eq!(
        (w.bore_radius, w.wall_thickness, w.axial_length),
        (a, 0.005, along)
    );
    assert!((2.0 * a * (0.5 * w.arc_width / a).sin() / across - 1.0).abs() < 1e-12);
    assert_eq!((w.cheeks, w.ramp), (WindowCheeks::Parallel, None));
    assert_eq!(w.cover.map(|cover| cover.open_width), Some(open));
    assert_eq!(w.upstream, WindowUpstream::ClosedEnd { distance: 0.017 });
    let riser_area = core::f64::consts::PI * 0.0061025 * 0.0061025;
    assert!((across * along / riser_area - 1.0).abs() <= 1e-3);
    // A covered hole is a staircase: three rungs at their observed order.
    assert_eq!(record.ladder.rungs.len(), 3);
    assert!((0.5..=2.0).contains(&record.ladder.order));
    // The lip adds mass to the same hole.
    let uncovered = wall_window_ladder(
        &WallWindow { cover: None, ..w },
        &WallWindowOptions {
            workers: 8,
            ..WallWindowOptions::rung(&w, record.cells[0])
        },
    )
    .unwrap();
    let bore_area = core::f64::consts::PI * a * a;
    eprintln!(
        "flute G4 embouchure: {:.4} x {:.4} mm, {:.4} mm open: {:.4} mm of bore covered (order {:.3}), {:.4} mm uncovered (order {:.3})",
        across * 1e3,
        along * 1e3,
        open * 1e3,
        record.ladder.mass_per_density * bore_area * 1e3,
        record.ladder.order,
        uncovered.mass_per_density * bore_area * 1e3,
        uncovered.order
    );
    assert!(record.ladder.mass_per_density > uncovered.mass_per_density);
    // It radiates from the outer cylinder through the riser's top, with no
    // authored loss.
    let voice = file.renderer().voice();
    let load = voice.window_radiation_load().expect("window radiation");
    assert_eq!(load.outer_radius_m(), a + 0.005);
    assert_eq!(
        load.cells(),
        record.ladder.rungs[record.ladder.rungs.len() - 1].outer_cells.len()
    );
    assert!(load.band_error() <= 0.05);
    assert!(text.contains("\nflue 0.00594 0.00090 0.00708 0.00029 0 0\n"));
    render(file);
}

#[test]
fn malformed_embouchure_rows_refuse() {
    let (text, [across, along, open], row) = flute();
    let refused = |new: &str| parse(&text.replace(&row, new)).is_err();
    assert!(refused("window embouchure 1 0.0098740 0.0118488 0.0054355"));
    assert!(refused("window embouchure 4 0.0098740 0.0118488"));
    assert!(refused(&format!("window embouchure 4 {across} {along} {across}")));
    assert!(refused(&format!("window embouchure 4 {across} {along} 0")));
    // The hole must be the riser's area.
    assert!(refused(&format!(
        "window embouchure 4 {} {along} {open}",
        1.01 * across
    )));
    let hole = format!("window embouchure 4 {across} {along} {open}");
    for lip in [
        " lip",
        " lip sideways",
        " lip one-sided face",
        " lip one-sided face -0.001",
        " lip centred face 0.001",
        " lip one-sided face 0.001 0.002",
    ] {
        assert!(refused(&format!("{hole}{lip}")), "{lip}");
    }
}

#[test]
fn a_one_sided_lip_and_the_players_face_come_from_the_row() {
    // The lip covers the hole from the player's side, leaving its strip
    // against the far wall, and the face stands 1 mm behind the lip's edge.
    let (text, [across, along, open], row) = flute();
    let hole = format!("window embouchure 8 {across} {along} {open}");
    let solve = |lip: &str| {
        let file = parse(&text.replace(&row, &format!("{hole}{lip}"))).unwrap();
        let record = file.window().expect("window record").clone();
        let voice = file.renderer().voice();
        let load = voice.window_radiation_load().expect("window radiation");
        assert!(load.band_error() <= 0.05, "{lip}: band error {}", load.band_error());
        (record, voice.mouth_inertance_pa_s2_m3())
    };
    let (centred, _) = solve("");
    let (lip, lip_inertance) = solve(" lip one-sided");
    let (face, face_inertance) = solve(" lip one-sided face 0.001");
    assert_eq!(
        centred.window.cover,
        Some(WindowCover {
            open_width: open,
            side: CoverSide::Centred,
        })
    );
    assert_eq!(
        lip.window.cover,
        Some(WindowCover {
            open_width: open,
            side: CoverSide::OneSided { face: None },
        })
    );
    assert_eq!(
        face.window.cover.map(|cover| cover.side),
        Some(CoverSide::OneSided { face: Some(0.001) })
    );
    // The open strip runs OPEN wide from the far wall on the outer surface,
    // off the hole's centre.
    let cells = &lip.ladder.rungs[lip.ladder.rungs.len() - 1].outer_cells;
    let outer = 0.0087744 + 0.005;
    let far = cells.iter().map(|c| c.theta + 0.5 * c.arc).fold(f64::MIN, f64::max);
    let near = cells.iter().map(|c| c.theta - 0.5 * c.arc).fold(f64::MAX, f64::min);
    let edge = ((outer * far.sin() - open) / outer).asin();
    assert!((near - edge).abs() < 1e-9, "strip from {near} rad, edge at {edge} rad");
    let centroid = cells.iter().map(|c| c.theta * c.flux).sum::<f64>()
        / cells.iter().map(|c| c.flux).sum::<f64>();
    assert!(centroid > 0.0, "centroid {centroid}");
    // The face adds mass, and the voice's inertance carries it.
    assert!(face.ladder.mass_per_density > lip.ladder.mass_per_density);
    assert!(face_inertance > lip_inertance);
    let bore_area = core::f64::consts::PI * 0.0087744 * 0.0087744;
    eprintln!(
        "flute G4 at 8 axial cells, {:.4} mm open: centred {:.4} mm, one-sided {:.4} mm, face 1 mm behind {:.4} mm of bore",
        open * 1e3,
        centred.ladder.mass_per_density * bore_area * 1e3,
        lip.ladder.mass_per_density * bore_area * 1e3,
        face.ladder.mass_per_density * bore_area * 1e3
    );
}

const RECORDER_A1: &str =
    include_str!("../../../data/instruments/descriptions/recorder/recorder-a1.performance");

/// The alto recorder's a1 card, cut to a 1024-sample window with its onset
/// only (the release at 2.2 s lies outside that window).
fn recorder_a1() -> String {
    assert!(RECORDER_A1.contains("\noutput 120000 5 "));
    let events = "events\t2\nevent\t0e0\t5e-2\tpressure-pa\t2.1e2\nevent\t2.2e0\t1.05e-2\tpressure-pa\t0e0";
    assert!(RECORDER_A1.contains(events), "the a1 gesture changed; update this test");
    RECORDER_A1
        .replace("\noutput 120000 5 ", "\noutput 1024 1 ")
        .replace(events, "events\t1\nevent\t0e0\t5e-2\tpressure-pa\t2.1e2")
}

#[test]
fn ernoult_fabre_window_is_their_eq_18_moved_from_the_edge_to_the_junction() {
    use core::f64::consts::{FRAC_PI_2, PI};
    let (w, h) = (0.0046, 0.012);
    // At the right angle (no ramp) the ears are eq. 18's straight chimney.
    let (le, rp) = (0.0066503, 0.0097564);
    let right = ErnoultFabreWindow::new(w, h, le, None, rp);
    assert_eq!(right.ears_per_density.to_bits(), (le / (w * h)).to_bits());
    // Eq. 23 tends to it there, plus its own small-angle blend 32^-alpha 0.85 sqrt(l_e/H).
    let near = ErnoultFabreWindow::new(w, h, le, Some(FRAC_PI_2 - 1e-7), rp);
    let blend = fs_math::det::exp(-FRAC_PI_2 * fs_math::det::ln(32.0)) * 0.85
        * fs_math::det::sqrt(le / h) / h;
    assert!(((near.ears_per_density - right.ears_per_density) / blend - 1.0).abs() < 1e-3);
    // The edge plane is eq. 18's sum; the junction lies W/2 of bore upstream.
    let ramp = ErnoultFabreWindow::new(w, h, le, Some(18.6210f64.to_radians()), rp);
    let window_area = w * h;
    assert_eq!(
        ramp.edge_mass_per_density.to_bits(),
        (ramp.right_angle_per_density + ramp.ears_per_density
            + ramp.radiation_length_m / window_area)
            .to_bits()
    );
    let bore_area = PI * rp * rp;
    assert!(
        ((ramp.edge_mass_per_density - ramp.junction_mass_per_density) * bore_area - 0.5 * w)
            .abs()
            < 1e-15
    );
    // An independent evaluation (scratchpad passive_r2.py, posted in comment
    // 8343): 56.38 mm at the edge, 54.08 at the junction on the origin card's
    // bore, 53.08 at the junction on the card's area-equivalent head.
    let origin = ErnoultFabreWindow::new(w, h, 0.006575, Some(18.6210f64.to_radians()), 0.0098317);
    let origin_area = PI * 0.0098317 * 0.0098317;
    for (value, expected) in [
        (origin.edge_mass_per_density * origin_area, 0.056_383_309),
        (origin.junction_mass_per_density * origin_area, 0.054_083_309),
        (ramp.junction_mass_per_density * bore_area, 0.053_080_068),
    ] {
        assert!((value / expected - 1.0).abs() < 1e-7, "{value} vs {expected}");
    }
    assert_eq!(ramp.outside_spans(), None, "the card lies inside their FE set");
    assert!(ErnoultFabreWindow::new(w, h, le, Some(3f64.to_radians()), rp).outside_spans().is_some());
    assert!(ErnoultFabreWindow::new(w, 0.005, le, Some(0.3), rp).outside_spans().is_some());
}

#[test]
fn recorder_window_takes_ernoult_fabre_and_keeps_the_computed_mass_for_comparison() {
    let file = parse(&recorder_a1()).unwrap();
    let record = file.window().expect("window record");
    let ef = record.ernoult_fabre.expect("the card selects mass ernoult-fabre-2017");
    assert_eq!(
        (ef.width_m, ef.height_m, ef.ears_m, ef.bore_radius_m),
        (0.0046, 0.012, 0.0066503, 0.0097564),
        "W, H, l_e and r_p come from the card's flue row, channel and bore"
    );
    assert_eq!(ef.angle_rad.to_bits(), 18.6210f64.to_radians().to_bits());
    assert_eq!(record.mass_per_density().to_bits(), ef.junction_mass_per_density.to_bits());
    let voice = file.renderer().voice();
    let load = voice.window_radiation_load().expect("window radiation");
    let rho = voice.gas().density;
    assert_eq!(
        voice.mouth_inertance_pa_s2_m3().to_bits(),
        (rho * (ef.junction_mass_per_density - record.riser_mass_per_density)
            - load.load().static_inertance_pa_s2_m3())
        .to_bits(),
        "the network carries Ernoult and Fabre's mass at the junction"
    );
    let bore_area = core::f64::consts::PI * 0.0097564 * 0.0097564;
    let computed = record.ladder.mass_per_density * bore_area;
    assert!(computed > 0.0);
    eprintln!(
        "a1 window: Ernoult-Fabre {:.3} mm at the edge, {:.3} mm at the junction; computed (comparison) {:.3} mm",
        ef.edge_mass_per_density * bore_area * 1e3,
        ef.junction_mass_per_density * bore_area * 1e3,
        computed * 1e3
    );
}

#[test]
fn ernoult_fabre_refuses_outside_its_geometry() {
    let base = recorder_a1();
    let refused = |text: &str| parse(text).is_err();
    let row = "window wall-cut 4 cheeks parallel ramp 0.00051 18.6210 mass ernoult-fabre-2017";
    assert!(base.contains(row));
    assert!(refused(&base.replace(row, "window wall-cut 4 cheeks radial ramp 0.00051 18.6210 mass ernoult-fabre-2017")));
    assert!(refused(&base.replace(row, "window wall-cut 4 cheeks parallel ramp 0.00051 3.0 mass ernoult-fabre-2017")));
    assert!(refused(&base.replace(row, "window wall-cut 4 cheeks parallel ramp 0.00051 18.6210 mass ernoult-fabre")));
    // Their closed end is flush with the window's upstream edge.
    assert!(refused(&base.replace(
        "duct_section 1 2 0.0023000 0.0098200",
        "duct_section 1 2 0.0030000 0.0098200"
    )));
}

#[test]
fn inner_length_is_the_dubos_dalmont_law_on_open_side_holes_only() {
    // Garcia Mayen et al. (HAL hal-03230052, 2021): b = 4 mm on a = 7.3 mm, "t_i = 2.1 mm".
    let t = HoleInnerLength::law(0.004, 0.0073).unwrap();
    assert!((t - 0.0021).abs() < 0.05e-3, "{t}");
    // A small hole: the isolated-orifice 0.82 b.
    assert!((HoleInnerLength::law(1e-6, 1.0).unwrap() / 1e-6 - 0.82).abs() < 1e-5);
    assert_eq!(HoleInnerLength::law(0.0, 0.0073), None);
    assert_eq!(HoleInnerLength::law(0.008, 0.0073), None);

    let base = recorder_a1();
    let file = parse(&base).unwrap();
    let holes = file.hole_inner_lengths();
    assert_eq!(holes.len(), 2, "a1 opens holes 6 and 7");
    let graph = &file.renderer().voice().network().expect("graph").source().graph;
    for (hole, (section, b, a, chimney)) in holes.iter().zip([
        (18, 0.002675, 0.0075740, 0.0049260),
        (20, 0.002625, 0.0070111, 0.0089889),
    ]) {
        assert_eq!((hole.section, hole.hole_radius_m, hole.bore_radius_m, hole.chimney_m),
            (section, b, a, chimney));
        assert_eq!(hole.inner_length_m.to_bits(), HoleInnerLength::law(b, a).unwrap().to_bits());
        assert_eq!(graph.sections[section].length_m.to_bits(), (chimney + hole.inner_length_m).to_bits(),
            "t_i is added once, to the chimney");
    }
    let refused = |text: &str| parse(text).is_err();
    // A closed hole (rigid end), a bore section, and a hole wider than its bore.
    assert!(refused(&base.replace(
        "duct_section 16 17 0.0046850 0.0032000 1e-10\n",
        "duct_section 16 17 0.0046850 0.0032000 1e-10 inner 0.0080650\n"
    )));
    assert!(refused(&base.replace(
        "duct_section 18 20 0.0336000 0.0072673 1e-10\n",
        "duct_section 18 20 0.0336000 0.0072673 1e-10 inner 0.0080000\n"
    )));
    assert!(refused(&base.replace("inner 0.0075740", "inner 0.0020000")));
}

/// Garcia Mayen, Kergomard, Vergez, Guillemain, Jousserand, Pachebat and
/// Sanchez, "Characterization of open woodwind toneholes by the tube reversed
/// method" (HAL hal-03230052, 2021): one cylindrical tonehole, b = 4 mm,
/// t = 8.5 mm (on the hole's axis from the main tube's wall, their Fig. 1),
/// on a main tube of a = 7.3 mm, L1 = 44 mm and L2 = 74 mm. The card's own
/// hole code (the duct parser with `inner`, the compact baffled outlet
/// binding, the fractional waveguide) is driven at the inlet by an impulse;
/// Z_in = P/U from the network's own inlet pressure and flow; L1 and L2 (rigid
/// far end) are removed as lossless lines; `t_s = Im(Z_s) / (k Z_ch)` (their
/// eq. 13) is averaged over 400-600 Hz (their conclusion's band). Pass rule
/// and prediction, posted before the run (bead pexfp.138): inside their
/// measured 12.9-14.4 mm (Figs. 10-12, read by eye); with minus without
/// `inner` 2.158 mm within 0.05 mm; both runs within 0.05 mm of the predicted
/// 14.098 and 11.940 mm.
#[test]
fn a_measured_tonehole_through_the_cards_hole_code() {
    let ts = |inner: bool| -> f64 {
        let mut text = String::new();
        for line in recorder_a1().lines() {
            if line.starts_with("duct_node") || line.starts_with("duct_section")
                || line.starts_with("loss_section") || line.starts_with("outlet ")
                || line.starts_with("mouth_radiation") || line.starts_with("mouth_reaction")
                || line.starts_with("window ")
            {
                continue;
            }
            let hole = if inner { " inner 0.0073" } else { "" };
            let line = match line {
                "network 25 24" => format!("network 4 3\nduct_node inlet\nduct_node junction\nduct_node reflection -1\nduct_node reflection 1\nduct_section 0 1 0.044 0.0073 1e-10\nduct_section 1 2 0.0085 0.004 1e-10{hole}\nduct_section 1 3 0.074 0.0073 1e-10"),
                "viscothermal sections 13" => "viscothermal none".to_string(),
                "outlets 24 2" => "outlets 2 0".to_string(),
                other => other.to_string(),
            };
            text.push_str(&line);
            text.push('\n');
        }
        let file = parse(&text).unwrap();
        let network = file.renderer().voice().network().expect("graph");
        let gases = network.section_gases().to_vec();
        let (policy, _) = network.fractional_sections().expect("fractional travel");
        let dt = 1.0 / f64::from(file.info().mechanics_rate_hz);
        let (mut wave, _, _) = network
            .graph()
            .clone()
            .realize_fractional_waveguide(&gases, dt, policy, &CancelGate::new_clock_free())
            .unwrap();
        let n = 1 << 16;
        let (mut p, mut u) = (Vec::with_capacity(n), Vec::with_capacity(n));
        for i in 0..n {
            let frame = wave.step(if i == 0 { 1.0 } else { 0.0 }).unwrap();
            p.push(frame.inlet_pressure_pa);
            u.push(frame.inlet_flow_m3_s);
        }
        let peak = p.iter().fold(0.0f64, |m, v| m.max(v.abs()));
        let tail = p[n - 4096..].iter().fold(0.0f64, |m, v| m.max(v.abs()));
        assert!(tail < 1e-9 * peak, "the response must have died away: {tail} of {peak}");
        type C = (f64, f64);
        let div = |a: C, b: C| {
            let d = b.0 * b.0 + b.1 * b.1;
            ((a.0 * b.0 + a.1 * b.1) / d, (a.1 * b.0 - a.0 * b.1) / d)
        };
        let dtft = |x: &[f64], omega: f64| {
            x.iter().enumerate().fold((0.0, 0.0), |(re, im), (i, v)| {
                let phase = omega * i as f64;
                (re + v * fs_math::det::cos(phase), im - v * fs_math::det::sin(phase))
            })
        };
        let (main, hole) = (gases[0], gases[1]);
        let zc = main.density * main.sound_speed / (core::f64::consts::PI * 0.0073 * 0.0073);
        let zch = hole.density * hole.sound_speed / (core::f64::consts::PI * 0.004 * 0.004);
        let frequencies: Vec<f64> = (0..=20).map(|i| 400.0 + 10.0 * f64::from(i)).collect();
        let mean = frequencies.iter().map(|&f| {
            let omega = core::f64::consts::TAU * f;
            let z_in = div(dtft(&p, omega * dt), dtft(&u, omega * dt));
            let k = omega / main.sound_speed;
            let t1 = fs_math::det::tan(k * 0.044);
            let t2 = fs_math::det::tan(k * 0.074);
            // Remove L1: Z_j = Zc (Z_in - j Zc t1) / (Zc - j Z_in t1).
            let num = (z_in.0 * zc, (z_in.1 - zc * t1) * zc);
            let den = (zc + z_in.1 * t1, -z_in.0 * t1);
            let z_j = div(num, den);
            // Remove the rigid-ended L2 in parallel: Y_s = 1/Z_j - j t2 / Zc.
            let y_j = div((1.0, 0.0), z_j);
            let z_s = div((1.0, 0.0), (y_j.0, y_j.1 - t2 / zc));
            z_s.1 / (omega / hole.sound_speed * zch)
        }).sum::<f64>() / frequencies.len() as f64;
        mean
    };
    let with = ts(true);
    let without = ts(false);
    eprintln!(
        "Garcia Mayen tonehole through the card's code, t_s 400-600 Hz: with inner {:.3} mm, without {:.3} mm, difference {:.3} mm",
        with * 1e3, without * 1e3, (with - without) * 1e3
    );
    assert!((0.0129..=0.0144).contains(&with), "inside the measured envelope: {with}");
    assert!(((with - without) - 0.002158).abs() <= 0.05e-3, "t_i counted once: {}", with - without);
    assert!((with - 0.014098).abs() <= 0.05e-3, "prediction with inner: {with}");
    assert!((without - 0.011940).abs() <= 0.05e-3, "prediction without inner: {without}");
    assert!(without < 0.0129, "the planted negative: the origin hole code falls below the envelope");
}
