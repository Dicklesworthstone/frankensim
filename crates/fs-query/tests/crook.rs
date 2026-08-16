//! Crook delta-L battery (music bead `frankensim-music-v8-root-3ez8g.3.2`).
//!
//! - gk-001 the 3-valve fixture: three bent-torus crook lumens with
//!   ANALYTIC centerline arc lengths (angle x major radius) extract to
//!   certified `ValveChartDelta`s inside the authored band, plus the
//!   slide-range certificate from two of them.
//! - gk-002 band refusals by name: a lying CAD length refuses
//!   `DeltaLengthOutOfBand`; an undeclared port-area step refuses
//!   `JunctionAreaMismatch` while the SAME geometry with the step
//!   DECLARED passes (the acoustic-feature-vs-error distinction).
//! - gk-003 THE COMPOSED ORACLE: inserting the extracted delta-L into a
//!   cylinder duct moves the TMM first impedance peak by the
//!   theoretical open-pipe ratio `f ~ 1/(L + dL + end corr)` within an
//!   authored cents band — the emergent-pitch precondition.

use asupersync::types::Budget;
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};
use fs_query::{
    BoreConfig, CrookCadIntent, CrookConfig, CrookError, certify_slide_range, extract_crook_delta,
};
use fs_rep_frep::{BoolOp, BoolStyle, FrepBuilder};
use fs_rep_mesh::{DcOptions, dual_contour};

const SUITE: &str = "fs-query/crook";

fn verdict(case: &str, pass: bool, detail: &str) {
    println!(
        "{{\"suite\":\"{SUITE}\",\"case\":\"{case}\",\"verdict\":\"{}\",\"detail\":\"{detail}\"}}",
        if pass { "pass" } else { "fail" }
    );
    assert!(pass, "case {case}: {detail}");
}

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0xC400,
                kernel_id: 41,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

/// A crook lumen: a torus segment of `angle` radians at major radius
/// `major`, tube radius `minor` — analytic centerline arc = angle*major.
fn crook_lumen(major: f64, minor: f64, angle: f64) -> fs_rep_frep::Frep {
    let mut b = FrepBuilder::new();
    let torus = b
        .torus(Point3::new(0.0, 0.0, 0.0), major, minor)
        .expect("torus");
    // Keep the wedge 0..angle: half-space x >= 0 rotated cuts. For angle
    // <= pi/2 the quadrant pair (x >= 0 rotated by angle) works: keep
    // y >= 0 and the rotated half-space.
    let hy = b.half_space(Vec3::new(0.0, -1.0, 0.0), 0.0).expect("hy");
    // Rotated plane with outward normal (sin a, -cos a, 0): keeps points
    // with azimuth <= angle.
    // Keep azimuth <= angle: n . p = r sin(theta - angle) <= 0.
    let hcut = b
        .half_space(Vec3::new(-angle.sin(), angle.cos(), 0.0), 0.0)
        .expect("cut");
    let wedge = b
        .boolean(BoolOp::Intersect, BoolStyle::Hard, hy, hcut)
        .expect("wedge");
    let root = b
        .boolean(BoolOp::Intersect, BoolStyle::Hard, torus, wedge)
        .expect("root");
    b.finish(root).expect("frep")
}

fn config() -> CrookConfig {
    CrookConfig {
        bore: BoreConfig::default(),
        length_band: 0.10,
        junction_band: 0.15,
    }
}

#[test]
fn gk_001_three_valve_fixture_and_slide_range() {
    with_cx(|cx| {
        // Trumpet-like valve set: valve 2 (short), valve 1, valve 3 —
        // arc lengths angle*major with a shared tube radius.
        let minor = 0.02f64;
        let main_area = core::f64::consts::PI * minor * minor;
        let valves = [
            ("valve-2", 0.30f64, 0.25f64 * core::f64::consts::PI), // ~0.236 m
            ("valve-1", 0.35, 0.40 * core::f64::consts::PI),       // ~0.440 m
            ("valve-3", 0.40, 0.45 * core::f64::consts::PI),       // ~0.565 m
        ];
        let cfg = config();
        let mut deltas = Vec::new();
        let mut worst_dev = 0.0f64;
        for (label, major, angle) in valves {
            let lumen = crook_lumen(major, minor, angle);
            let (soup, _) = dual_contour(&lumen, DcOptions::sharp(0.008), cx).expect("dc");
            let cad = CrookCadIntent {
                centerline_length_m: major * angle,
                port_radius_m: minor,
                declared_area_step: 1.0,
            };
            let delta = extract_crook_delta(&lumen, &soup, &cfg, &cad, main_area, label, cx)
                .expect("crook extracts");
            println!("{}", delta.debug_line());
            worst_dev = worst_dev.max(delta.length_deviation);
            deltas.push(delta);
        }
        // Deltas must ORDER like the valves (2 < 1 < 3), the ordering the
        // harmonic series depends on.
        let ordered =
            deltas[0].delta_l_m < deltas[1].delta_l_m && deltas[1].delta_l_m < deltas[2].delta_l_m;
        let slide = certify_slide_range("tuning-slide", &deltas[0], &deltas[2]).expect("range");
        let inverted = certify_slide_range("bad", &deltas[2], &deltas[0]).is_err();
        let pass =
            ordered && worst_dev < 0.10 && slide.max_delta_l_m > slide.min_delta_l_m && inverted;
        verdict(
            "gk-001-three-valves",
            pass,
            &format!(
                "deltas {:.4}/{:.4}/{:.4} m ordered {ordered}, worst CAD deviation \
                 {worst_dev:.4}, slide range [{:.4}, {:.4}] certified, inverted refused \
                 {inverted}",
                deltas[0].delta_l_m,
                deltas[1].delta_l_m,
                deltas[2].delta_l_m,
                slide.min_delta_l_m,
                slide.max_delta_l_m
            ),
        );
    });
}

#[test]
fn gk_002_band_refusals_by_name() {
    with_cx(|cx| {
        let minor = 0.02f64;
        let main_area = core::f64::consts::PI * minor * minor;
        let major = 0.30f64;
        let angle = 0.25 * core::f64::consts::PI;
        let lumen = crook_lumen(major, minor, angle);
        let (soup, _) = dual_contour(&lumen, DcOptions::sharp(0.008), cx).expect("dc");
        let cfg = config();
        // Lying CAD length (20% long) refuses by name.
        let lying_cad = CrookCadIntent {
            centerline_length_m: 1.2 * major * angle,
            port_radius_m: minor,
            declared_area_step: 1.0,
        };
        let length_refused = matches!(
            extract_crook_delta(&lumen, &soup, &cfg, &lying_cad, main_area, "lying", cx),
            Err(CrookError::DeltaLengthOutOfBand { .. })
        );
        // Undeclared port step: main bore area 1.5x the crook's refuses...
        let honest_cad = CrookCadIntent {
            centerline_length_m: major * angle,
            port_radius_m: minor,
            declared_area_step: 1.0,
        };
        let undeclared = matches!(
            extract_crook_delta(
                &lumen,
                &soup,
                &cfg,
                &honest_cad,
                1.5 * main_area,
                "step",
                cx
            ),
            Err(CrookError::JunctionAreaMismatch { .. })
        );
        // ...while the SAME geometry with the step DECLARED passes: the
        // acoustic-feature-vs-extraction-error distinction.
        let declared_cad = CrookCadIntent {
            centerline_length_m: major * angle,
            port_radius_m: minor,
            declared_area_step: 1.0 / 1.5,
        };
        let declared_ok = extract_crook_delta(
            &lumen,
            &soup,
            &cfg,
            &declared_cad,
            1.5 * main_area,
            "declared-step",
            cx,
        )
        .is_ok();
        let pass = length_refused && undeclared && declared_ok;
        verdict(
            "gk-002-band-refusals",
            pass,
            &format!(
                "lying CAD length refused {length_refused}; undeclared step refused \
                 {undeclared}; declared step passes {declared_ok}"
            ),
        );
    });
}

#[test]
fn gk_003_composed_peak_shift_oracle() {
    with_cx(|cx| {
        use fs_duct::{Duct, LossModel, Segment, Termination, impedance_peaks, impedance_sweep};
        use fs_material::gas::{GasSpec, GasState};
        // Extract a real crook delta, insert it into a cylinder duct, and
        // check the TMM first-peak shift against open-pipe theory
        // f1 ~ c / 2(L + dL + end corrections): the emergent-pitch
        // precondition.
        let minor = 0.02f64;
        let main_area = core::f64::consts::PI * minor * minor;
        let major = 0.35f64;
        let angle = 0.40 * core::f64::consts::PI;
        let lumen = crook_lumen(major, minor, angle);
        let (soup, _) = dual_contour(&lumen, DcOptions::sharp(0.008), cx).expect("dc");
        let cad = CrookCadIntent {
            centerline_length_m: major * angle,
            port_radius_m: minor,
            declared_area_step: 1.0,
        };
        let delta = extract_crook_delta(&lumen, &soup, &config(), &cad, main_area, "v1", cx)
            .expect("crook extracts");
        let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).expect("air");
        let base_l = 0.60f64;
        let first_peak = |length: f64| -> f64 {
            let duct = Duct {
                segments: vec![Segment::Cylinder {
                    radius: minor,
                    length,
                }],
            };
            let sweep = impedance_sweep(
                &duct,
                &gas,
                core::f64::consts::TAU * 60.0,
                core::f64::consts::TAU * 400.0,
                4000,
                LossModel::WideTube,
                Termination::UnflangedOpen,
            )
            .expect("sweep");
            let peaks = impedance_peaks(&sweep);
            sweep[peaks[0]].omega / core::f64::consts::TAU
        };
        let f_open = first_peak(base_l);
        let f_valved = first_peak(base_l + delta.delta_l_m);
        let measured_cents = 1200.0 * (f_valved / f_open).log2();
        // Theory: closed-open quarter-wave with the unflanged end
        // correction 0.6133 a on the open end.
        let corr = 0.6133 * minor;
        let predicted_cents = 1200.0 * ((base_l + corr) / (base_l + delta.delta_l_m + corr)).log2();
        let dev = (measured_cents - predicted_cents).abs();
        let pass = measured_cents < -400.0 && dev < 12.0;
        verdict(
            "gk-003-composed-peak-shift",
            pass,
            &format!(
                "inserting dL = {:.4} m: first peak {f_open:.2} -> {f_valved:.2} Hz \
                 ({measured_cents:.1} cents; theory {predicted_cents:.1}; dev {dev:.1} \
                 cents inside the authored 12-cent band)",
                delta.delta_l_m
            ),
        );
    });
}
