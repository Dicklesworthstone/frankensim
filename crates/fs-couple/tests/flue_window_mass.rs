//! `window wall-cut` on the Clare D whistle's D5 card: the window derived from
//! the card's own channel, pocket and bore, its computed mass replacing the
//! legacy `0.61 a_eq` inertance, parallel cheeks and a labium ramp from the
//! row, the pocket radius tolerance, the refusals, and the window's radiation
//! from the head's outer cylinder in place of the authored loss; and
//! `window embouchure` on the G4 flute card: the hole under the lip.

use fs_bem::wall_window::{
    WallWindow, WallWindowOptions, WindowCheeks, WindowCover, WindowRamp, WindowUpstream,
    wall_window_ladder,
};
use fs_couple::flue_loop::network::performance::file::{FlueFilePerformance, FlueInputError};
use fs_couple::pcm_wav::observation::PressureRenderer;
use fs_exec::CancelGate;

const CARD: &str =
    include_str!("../../../data/instruments/descriptions/tin-whistle/tin-whistle-d5.performance");

/// The card's flue row: the window's radiation replaces the authored loss.
const FLUE: &str = "flue 0.00800 0.00148 0.00530 0 0 0\n";

/// The card without its window radiation and window rows (lossless mouth),
/// cut to a 1024-sample window whose only gesture event is the opening
/// pressure ramp.
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
    let text = CARD
        .replace(&radiating, "")
        .replace("output 120000 5", "output 1024 1");
    let events = "events\t2\nevent\t0e0\t5e-2\tpressure-pa\t6.715e1\nevent\t2.2e0\t1.05e-2\tpressure-pa\t0e0";
    assert!(
        text.contains(events),
        "the D5 gesture changed; update this test"
    );
    text.replace(events, "events\t1\nevent\t0e0\t5e-2\tpressure-pa\t6.715e1")
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
    assert_eq!(w.cover, Some(WindowCover { open_width: open }));
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
}
