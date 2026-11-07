//! Per-opening quadrature: an opening the declared grid resolves keeps it bit for bit, an opening it does not resolve
//! (the 19 mm bell of the oboe's 20 kHz DG card) takes the smallest grid that does, and every other refusal stays.
use super::*;
use core::f64::consts::TAU;

fn medium() -> RayleighMedium {
    RayleighMedium { density: 1.2, sound_speed: 343.2 }
}

fn receiver(rings: usize, angles: usize, band_hz: f64) -> CircularOutletReceiver {
    CircularOutletReceiver {
        position_m: [0.0, 0.0, 0.5],
        radial_rings: rings,
        angular_points: angles,
        maximum_frequency_hz: band_hz,
    }
}

/// The openings of the oboe DG card: the 19 mm bell at the origin and the three OPEN chimneys
/// (the closed ones are not openings of the baffle).
fn card_openings() -> Vec<CircularOpening> {
    vec![
        CircularOpening { center_m: [0.0, 0.0], radius_m: 0.019 },
        CircularOpening { center_m: [0.171, 0.0], radius_m: 0.0038 },
        CircularOpening { center_m: [0.135, 0.0], radius_m: 0.0049 },
        CircularOpening { center_m: [0.095, 0.0], radius_m: 0.0054 },
    ]
}

/// What the build used for every opening before: the declared grid, translated to the opening's frame.
fn declared_locals(openings: &[CircularOpening], receiver: CircularOutletReceiver) -> Vec<CircularOutletReceiver> {
    openings
        .iter()
        .map(|opening| CircularOutletReceiver {
            position_m: [
                receiver.position_m[0] - opening.center_m[0],
                receiver.position_m[1] - opening.center_m[1],
                receiver.position_m[2],
            ],
            ..receiver
        })
        .collect()
}

fn brute_force_minimum(radius_m: f64, rings: usize, angles: usize, band_hz: f64) -> Option<(usize, usize)> {
    let mut best: Option<(usize, usize, usize)> = None;
    for r in rings..=128 {
        for a in angles..=512 {
            let phase = TAU * band_hz * (radius_m / (r as f64).sqrt() + TAU * radius_m / a as f64) / medium().sound_speed;
            if phase <= core::f64::consts::FRAC_PI_4 && best.is_none_or(|(count, _, _)| r * a < count) {
                best = Some((r * a, r, a));
            }
        }
    }
    best.map(|(_, r, a)| (r, a))
}

#[test]
fn openings_the_declared_grid_resolves_get_exactly_the_old_local_receivers() {
    let openings = card_openings();
    for (rings, angles, band) in [(32, 128, 9_600.0), (32, 128, 4_800.0), (128, 224, 19_000.0), (16, 64, 1_000.0)] {
        let declared = receiver(rings, angles, band);
        let locals = BaffledPressure::refined_local_receivers(&openings, declared, medium().sound_speed).unwrap();
        assert_eq!(locals, declared_locals(&openings, declared), "{rings} x {angles} at {band} Hz");
        assert_eq!(
            BaffledPressure::resolved_circular_grids(&openings, declared, medium()).unwrap(),
            vec![(rings, angles); 4]
        );
    }
}

#[test]
fn the_oboe_bell_at_20_khz_is_refused_by_the_shared_grid_and_refined_alone() {
    let openings = card_openings();
    let declared = receiver(32, 128, 20_000.0);
    // Today's refusal, from the untouched surface rule.
    let local = declared_locals(&openings, declared)[0];
    assert!(BaffledPressure::circular_surface(0.019, local, medium()).unwrap_err().contains("underresolves"));
    for (index, opening) in openings.iter().enumerate().skip(1) {
        let local = declared_locals(&openings, declared)[index];
        assert!(BaffledPressure::circular_surface(opening.radius_m, local, medium()).is_ok(), "chimney {index}");
    }
    let grids = BaffledPressure::resolved_circular_grids(&openings, declared, medium()).unwrap();
    let want = brute_force_minimum(0.019, 32, 128, 20_000.0).unwrap();
    assert_eq!(grids[0], want, "the bell takes the minimum grid by exhaustive search");
    assert!(want.0 * want.1 > 32 * 128 && want.0 <= 128 && want.1 <= 512);
    assert_eq!(grids[1..], [(32, 128); 3], "the chimneys keep the declared grid");
    let points: usize = grids.iter().map(|(r, a)| r * a).sum();
    assert!(points <= 120_000, "{points}");
    // The receiver builds on the refined grids (an unrefined shared grid cannot).
    let built = BaffledPressure::circular_openings_multirate_differentiating(&openings, 96_000, 48_000, declared, medium());
    assert!(built.is_ok(), "{:?}", built.err());
    println!(
        "{{\"suite\":\"fs-couple\",\"case\":\"opening-grids-oboe-20khz\",\"grids\":{grids:?},\"points\":{points}}}"
    );
}

#[test]
fn a_refined_sum_over_the_point_budget_and_an_unresolvable_band_keep_their_refusals() {
    // Six bells: the declared 32 x 128 sums to 24 576, each refines to the same 128 x 257-class grid and the sum passes 120 000.
    let bells: Vec<_> = (0..6)
        .map(|k| CircularOpening { center_m: [0.1 * k as f64, 0.0], radius_m: 0.019 })
        .collect();
    let declared = receiver(32, 128, 20_000.0);
    let error = BaffledPressure::resolved_circular_grids(&bells, declared, medium()).unwrap_err();
    assert!(error.contains("budget"), "{error}");
    // The array rule on the DECLARED grid is unchanged: 4 x 128 x 512 is over 120 000.
    let error = BaffledPressure::resolved_circular_grids(&card_openings(), receiver(128, 512, 9_600.0), medium()).unwrap_err();
    assert!(error.contains("budget"), "{error}");
    // No grid up to 128 x 512 resolves a 19 mm disk at 40 kHz: the declared grid is kept and the surface rule refuses as before.
    let band = receiver(32, 128, 40_000.0);
    let grids = BaffledPressure::resolved_circular_grids(&card_openings()[..1], band, medium()).unwrap();
    assert_eq!(grids, vec![(32, 128)]);
    let error = BaffledPressure::circular_openings_multirate_differentiating(
        &card_openings()[..1], 192_000, 96_000, band, medium(),
    )
    .err()
    .expect("refused");
    assert!(error.contains("underresolves"), "{error}");
}

#[test]
fn the_grid_test_agrees_with_the_surface_rule_on_a_sweep() {
    for radius in [1.0e-3, 5.4e-3, 19.0e-3, 40.0e-3] {
        for (rings, angles) in [(1, 8), (4, 32), (32, 128), (128, 512), (16, 512), (128, 8)] {
            for band in [500.0, 2_000.0, 9_600.0, 20_000.0] {
                let local = receiver(rings, angles, band);
                let by_rule = BaffledPressure::circular_surface(radius, local, medium()).is_ok();
                let by_helper = BaffledPressure::circular_grid_resolves(radius, rings, angles, band, medium().sound_speed);
                assert_eq!(by_helper, by_rule, "R {radius}, {rings} x {angles}, {band} Hz");
            }
        }
    }
}

#[test]
fn the_refined_bell_pressure_converges_to_a_finer_grid_off_axis() {
    // One 19 mm disk driven by a 15 kHz sine through the differentiating image, receiver off axis so rings AND angles matter:
    // the grid the rule picks for a 20 kHz band against an independent 128 x 512 grid. The bound is fixed before the run.
    let local = CircularOutletReceiver { position_m: [0.1, 0.0, 0.4], ..receiver(32, 128, 20_000.0) };
    let (rings, angles) = BaffledPressure::resolved_circular_grids(&card_openings()[..1], local, medium()).unwrap()[0];
    assert_eq!((rings, angles), brute_force_minimum(0.019, 32, 128, 20_000.0).unwrap());
    let fine = CircularOutletReceiver { radial_rings: 128, angular_points: 512, ..local };
    let refined = CircularOutletReceiver { radial_rings: rings, angular_points: angles, ..local };
    let mut a = BaffledPressure::circular_outlet_multirate_differentiating(0.019, 96_000, 48_000, refined, medium()).unwrap();
    let mut b = BaffledPressure::circular_outlet_multirate_differentiating(0.019, 96_000, 48_000, fine, medium()).unwrap();
    let (mut error, mut signal) = (0.0, 0.0);
    for frame in 0..1200 {
        let trace: Vec<f64> = (0..2)
            .map(|sub| 1e-5 * (TAU * 15_000.0 * (frame * 2 + sub) as f64 / 96_000.0).sin())
            .collect();
        let (pa, pb) = (a.step_trace(&trace).unwrap(), b.step_trace(&trace).unwrap());
        if frame > 400 {
            error += (pa - pb).powi(2);
            signal += pb.powi(2);
        }
    }
    let relative = (error / signal).sqrt();
    println!("{{\"suite\":\"fs-couple\",\"case\":\"opening-grid-convergence\",\"refined\":[{rings},{angles}],\"relative_rms_vs_128x512\":{relative:e}}}");
    assert!(relative < 0.02, "relative RMS {relative}");
}

#[test]
fn a_single_opening_is_still_the_disk_image() {
    let bell = [CircularOpening { center_m: [0.0, 0.0], radius_m: 0.019 }];
    let declared = receiver(32, 128, 9_600.0);
    let mut one = BaffledPressure::circular_openings_multirate_differentiating(&bell, 96_000, 48_000, declared, medium()).unwrap();
    let mut disk = BaffledPressure::circular_outlet_multirate_differentiating(0.019, 96_000, 48_000, declared, medium()).unwrap();
    for frame in 0..300 {
        let trace = [1e-5 * (frame as f64).sin(), 2e-5 * (frame as f64 * 0.7).cos()];
        assert_eq!(one.step_trace(&trace).unwrap().to_bits(), disk.step_trace(&trace).unwrap().to_bits());
    }
}
