//! Wall-window acoustic mass: the ladder, the inner end correction against
//! Dalmont's law, far-boundary and port invariance, an independent
//! implementation, worker determinism and the refusals.

use fs_bem::wall_window::{
    WallWindow, WallWindowError, WallWindowMass, WallWindowOptions, WindowOutside, WindowUpstream,
    wall_window_ladder, wall_window_mass,
};

/// The Clare D tin whistle's window (Davies 2006 rows of the instrument
/// ledger): bore 12.2 mm, wall 2.4 mm at the window, window 8 mm across and
/// 5.3 mm long, flat pocket end 10.65 mm from the window centre.
fn clare_d() -> WallWindow {
    WallWindow {
        bore_radius: 0.0061,
        wall_thickness: 0.0024,
        arc_width: 0.008,
        axial_length: 0.0053,
        upstream: WindowUpstream::ClosedEnd { distance: 0.01065 },
        outside: WindowOutside::Free,
    }
}

fn rung(axial: usize, arc: usize, layers: usize) -> WallWindowOptions {
    WallWindowOptions {
        axial_cells: axial,
        arc_cells: arc,
        wall_layers: layers,
        ..WallWindowOptions::default()
    }
}

fn solve(window: &WallWindow, options: &WallWindowOptions) -> WallWindowMass {
    let mass = wall_window_mass(window, options).expect("wall window solve");
    eprintln!(
        "a {:.2} mm t {:.2} mm window {:.3} x {:.2} mm, cells {}x{}x{} (N {}, bore {:?}, exterior {:?}, {} unknowns): bore length {:.6} mm, window length {:.6} mm, flux {:.6e} m, outer {:.6e} m, residual {:.1e}",
        window.bore_radius * 1e3,
        window.wall_thickness * 1e3,
        mass.resolved_arc_width * 1e3,
        window.axial_length * 1e3,
        options.axial_cells,
        options.arc_cells,
        options.wall_layers,
        mass.theta_cells,
        mass.bore_cells,
        mass.exterior_cells,
        mass.unknowns,
        mass.bore_length * 1e3,
        mass.window_length * 1e3,
        mass.flux,
        mass.outer_flux,
        mass.residual,
    );
    assert!(
        mass.residual < 1e-10,
        "face system residual {}",
        mass.residual
    );
    assert!(
        (mass.outer_flux - mass.flux).abs() <= 1e-9 * mass.flux,
        "the channel must conserve flux: {} in, {} out",
        mass.flux,
        mass.outer_flux
    );
    mass
}

/// Dalmont et al. (2002) inner length of a circular side hole of radius `b`
/// in a bore of radius `a`.
fn dalmont_inner(b: f64, a: f64) -> f64 {
    let d = b / a;
    ((((-0.71 * d + 1.27) * d - 1.09) * d - 0.193) * d + 0.82) * b
}

fn ladder(window: &WallWindow, rungs: &[WallWindowOptions]) -> (Vec<f64>, f64, f64) {
    let lengths: Vec<f64> = rungs
        .iter()
        .map(|options| solve(window, options).bore_length)
        .collect();
    let n = lengths.len();
    let (d1, d2) = (
        lengths[n - 3] - lengths[n - 2],
        lengths[n - 2] - lengths[n - 1],
    );
    assert!(d1 > 0.0 && d2 > 0.0, "the ladder must fall monotonically");
    let order = -(d2 / d1).log2();
    // Richardson at the theoretical order 4/3: 2^(4/3) - 1 = 1.5198...
    let extrapolated = lengths[n - 1] - d2 / 1.519_842_099_789_746_3;
    eprintln!(
        "ladder {:?} mm: observed order {order:.3}, extrapolated at 4/3 {:.5} mm",
        lengths.iter().map(|l| l * 1e3).collect::<Vec<_>>(),
        extrapolated * 1e3
    );
    (lengths, order, extrapolated)
}

#[test]
fn ladder_converges_at_the_re_entrant_edge_order() {
    // An arc of exactly a twelfth of the circumference, so that every rung
    // resolves the same window.
    let window = WallWindow {
        arc_width: 2.0 * core::f64::consts::PI * 0.0061 / 12.0,
        axial_length: 0.0032,
        upstream: WindowUpstream::ClosedEnd { distance: 0.008 },
        ..clare_d()
    };
    let (lengths, order, _) = ladder(&window, &[rung(4, 6, 4), rung(8, 12, 8), rung(16, 24, 16)]);
    // The velocity is singular like d^(-1/3) at the window's right-angle
    // edges, so the mass converges at order 4/3.
    assert!(
        (1.15..=1.5).contains(&order),
        "observed order {order}, expected about 4/3"
    );
    // The numpy prototype of the same discretization (see below).
    for (got, expected_mm) in lengths.iter().zip([52.3320511, 50.6949145, 50.0302013]) {
        assert!(
            (got * 1e3 - expected_mm).abs() <= 1e-8 * expected_mm,
            "Rust {:.10} mm vs prototype {expected_mm:.7} mm",
            got * 1e3
        );
    }
}

#[test]
#[ignore = "the Clare D window's production ladder to 64 x 96 cells: minutes, 1 GB"]
fn clare_d_window_ladder() {
    let fine = WallWindowOptions {
        max_unknowns: 8192,
        workers: 8,
        ..rung(64, 96, 64)
    };
    let (_, order, extrapolated) = ladder(&clare_d(), &[rung(16, 24, 16), rung(32, 48, 32), fine]);
    assert!((1.15..=1.5).contains(&order), "observed order {order}");
    // Longer than the 1-D card's riser plus 0.61 a_eq (12.80 mm of bore).
    assert!(extrapolated > 0.0128);
}

#[test]
fn matches_the_independent_dense_prototype() {
    // A numpy implementation of the same discretization, written first,
    // gives these values. It works in the window cells' basis with explicit
    // block inverses and a dense solve over all arc cells, with the port at
    // 40 mm and the exterior at 60 mm.
    let window = clare_d();
    for (options, expected_mm) in [
        (rung(8, 12, 8), 16.2288478951),
        (rung(16, 24, 16), 15.9707601308),
    ] {
        let options = WallWindowOptions {
            port_radii: 0.04 / 0.0061,
            exterior_radii: 0.06 / 0.0061,
            ..options
        };
        let got = solve(&window, &options).bore_length * 1e3;
        assert!(
            (got - expected_mm).abs() <= 1e-8 * expected_mm,
            "Rust {got:.10} mm vs prototype {expected_mm:.10} mm"
        );
    }
}

#[test]
fn inner_end_correction_matches_dalmont_within_shape_tolerance() {
    // A square window through a thin wall of a through bore (a symmetric T),
    // with a pressure release on its outer face. What is left after the
    // channel's own length is the inner correction. Dalmont's law is for a
    // circular hole; the equal-area square differs by the shape alone, about
    // 1 % for an isolated hole.
    let (a, t) = (0.0061, 0.0006);
    for side in [0.002, 0.004] {
        let window = WallWindow {
            bore_radius: a,
            wall_thickness: t,
            arc_width: side,
            axial_length: side,
            upstream: WindowUpstream::Port { distance: 0.04 },
            outside: WindowOutside::Release,
        };
        let options = WallWindowOptions {
            port_radii: 0.04 / a,
            ..rung(16, 16, 8)
        };
        let mass = solve(&window, &options);
        let b = (mass.window_area / core::f64::consts::PI).sqrt();
        let channel = a * (1.0 + t / a).ln();
        let inner = mass.window_length - channel;
        let law = dalmont_inner(b, a);
        eprintln!(
            "square {:.1} mm: b/a {:.3}, inner {:.4} b, Dalmont {:.4} b",
            side * 1e3,
            b / a,
            inner / b,
            law / b
        );
        assert!(
            (inner / law - 1.0).abs() < 0.03,
            "inner {inner} vs Dalmont {law}"
        );
    }
}

#[test]
fn far_boundaries_do_not_move_the_mass() {
    let window = clare_d();
    let base = solve(&window, &rung(8, 12, 8)).bore_length;
    let far_exterior = solve(
        &window,
        &WallWindowOptions {
            exterior_radii: 20.0,
            ..rung(8, 12, 8)
        },
    )
    .bore_length;
    let far_port = solve(
        &window,
        &WallWindowOptions {
            port_radii: 16.0,
            ..rung(8, 12, 8)
        },
    )
    .bore_length;
    eprintln!(
        "base {:.6} mm, exterior x2 {:.6} mm, port x2 {:.6} mm",
        base * 1e3,
        far_exterior * 1e3,
        far_port * 1e3
    );
    assert!((far_exterior / base - 1.0).abs() < 1e-3);
    // Plane flow is exact on the graded grid, so only the evanescent field
    // reaching the port can move the extrapolated junction.
    assert!((far_port / base - 1.0).abs() < 1e-4);
}

#[test]
fn a_flush_closed_end_lengthens_the_window() {
    let pocket = solve(&clare_d(), &rung(8, 12, 8)).bore_length;
    let flush = solve(
        &WallWindow {
            upstream: WindowUpstream::ClosedEnd { distance: 0.00265 },
            ..clare_d()
        },
        &rung(8, 12, 8),
    )
    .bore_length;
    eprintln!(
        "pocket {:.5} mm, flush end {:.5} mm",
        pocket * 1e3,
        flush * 1e3
    );
    assert!(flush > pocket);
}

#[test]
fn workers_do_not_change_a_bit() {
    let window = clare_d();
    let one = wall_window_mass(&window, &rung(8, 12, 8)).expect("one worker");
    let three = wall_window_mass(
        &window,
        &WallWindowOptions {
            workers: 3,
            ..rung(8, 12, 8)
        },
    )
    .expect("three workers");
    assert_eq!(one, three);
}

#[test]
fn two_rung_ladder_reaches_the_three_rung_limit() {
    let window = clare_d();
    let options = WallWindowOptions::rung(&window, 16);
    // 24 arc cells: N = 115 here and 230 at the doubling, both within 2e-4 of
    // exact, so both rungs resolve the same arc.
    assert_eq!(
        (options.axial_cells, options.arc_cells, options.wall_layers),
        (16, 24, 16)
    );
    let ladder = wall_window_ladder(&window, &options).expect("ladder");
    let area = core::f64::consts::PI * window.bore_radius * window.bore_radius;
    let length = ladder.mass_per_density * area;
    eprintln!(
        "rungs {:.6} {:.6} mm, extrapolated {:.6} mm",
        ladder.coarse.bore_length * 1e3,
        ladder.fine.bore_length * 1e3,
        length * 1e3
    );
    // The ignored clare_d_window_ladder extrapolates 16x24, 32x48 and 64x96
    // cells to 15.71259 mm.
    assert!((length / 0.01571259 - 1.0).abs() < 1e-4, "{length}");
}

#[test]
fn refusals() {
    let ok = clare_d();
    let options = rung(4, 6, 4);
    let invalid = |window: &WallWindow, options: &WallWindowOptions| {
        matches!(
            wall_window_mass(window, options),
            Err(WallWindowError::Invalid(_))
        )
    };
    assert!(invalid(
        &WallWindow {
            bore_radius: -1.0,
            ..ok
        },
        &options
    ));
    assert!(invalid(
        &WallWindow {
            wall_thickness: f64::NAN,
            ..ok
        },
        &options
    ));
    assert!(invalid(
        &WallWindow {
            arc_width: 2.0 * core::f64::consts::PI * ok.bore_radius,
            ..ok
        },
        &options
    ));
    assert!(invalid(
        &WallWindow {
            upstream: WindowUpstream::ClosedEnd { distance: 0.002 },
            ..ok
        },
        &options
    ));
    assert!(invalid(
        &WallWindow {
            upstream: WindowUpstream::Port { distance: 0.00265 },
            ..ok
        },
        &options
    ));
    assert!(invalid(
        &ok,
        &WallWindowOptions {
            axial_cells: 0,
            ..options
        }
    ));
    assert!(invalid(
        &ok,
        &WallWindowOptions {
            grading: 0.9,
            ..options
        }
    ));
    assert!(invalid(
        &ok,
        &WallWindowOptions {
            port_radii: 0.1,
            ..options
        }
    ));
    assert!(invalid(
        &ok,
        &WallWindowOptions {
            exterior_radii: 1.0,
            ..options
        }
    ));
    assert!(matches!(
        wall_window_mass(
            &ok,
            &WallWindowOptions {
                max_unknowns: 10,
                ..options
            }
        ),
        Err(WallWindowError::Budget {
            what: "coupled face unknowns",
            ..
        })
    ));
    assert!(matches!(
        wall_window_mass(
            &ok,
            &WallWindowOptions {
                grading: 1.0,
                max_line_cells: 32,
                ..options
            }
        ),
        Err(WallWindowError::Budget { .. })
    ));
}
