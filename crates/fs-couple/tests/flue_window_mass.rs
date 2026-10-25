//! `window wall-cut` on the Clare D whistle's D5 card: the window derived from
//! the card's own channel, pocket and bore, its computed mass replacing the
//! legacy `0.61 a_eq` inertance, and the refusals.

use fs_couple::flue_loop::network::performance::file::{FlueFilePerformance, FlueInputError};
use fs_couple::pcm_wav::observation::PressureRenderer;
use fs_exec::CancelGate;

const CARD: &str =
    include_str!("../../../data/instruments/descriptions/tin-whistle/tin-whistle-d5.performance");

/// The card without any window row, cut to a 1024-sample window whose only
/// gesture event is the opening pressure ramp.
fn legacy() -> String {
    let text = CARD
        .replace("window wall-cut 16\n", "")
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
    let bore_area = core::f64::consts::PI * 0.0061 * 0.0061;
    let length = record.ladder.mass_per_density * bore_area;
    // fs-bem's ignored three-rung ladder extrapolates to 15.71259 mm.
    assert!((length / 0.01571259 - 1.0).abs() < 1e-4, "{length}");
    let riser_area = core::f64::consts::PI * 0.0036737 * 0.0036737;
    assert_eq!(
        record.riser_mass_per_density.to_bits(),
        (0.0024 / riser_area).to_bits()
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
    // A stepped pocket is not represented.
    assert!(refused(&with_window(
        &base.replace(
            "duct_section 1 2 0.0106500 0.0061000",
            "duct_section 1 2 0.0106500 0.0062000"
        ),
        "window wall-cut 16"
    )));
    // A compact mouth reaction brings its own outer mass.
    assert!(refused(&with_window(
        &base
            .replace(
                "flue 0.00800 0.00148 0.00530 0 0 0.02",
                "flue 0.00800 0.00148 0.00530 0 0 0"
            )
            .replace(
                "\ncontrols 4096",
                "\nmouth_reaction compact-baffled\ncontrols 4096"
            ),
        "window wall-cut 16"
    )));
}
