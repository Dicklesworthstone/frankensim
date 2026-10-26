//! Wall-window acoustic mass: the ladder, the inner end correction against
//! Dalmont's law, far-boundary and port invariance, independent
//! implementations of the prism and of the staircase cuts, the bounds that
//! radial cuts put on parallel cheeks, worker determinism and the refusals.

use fs_bem::wall_window::{
    WallWindow, WallWindowError, WallWindowMass, WallWindowOptions, WindowCheeks, WindowOutside,
    WindowRamp, WindowUpstream, wall_window_ladder, wall_window_mass,
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
        cheeks: WindowCheeks::Radial,
        ramp: None,
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
    let outer_r = window.bore_radius + window.wall_thickness;
    let cell_flux: f64 = mass.outer_cells.iter().map(|c| c.flux).sum();
    let cell_area: f64 = mass
        .outer_cells
        .iter()
        .map(|c| outer_r * c.arc * c.length)
        .sum();
    assert!(
        (cell_flux / mass.outer_flux - 1.0).abs() < 1e-12
            && (cell_area / mass.outer_area - 1.0).abs() < 1e-12,
        "outer cells carry {cell_flux} over {cell_area} m^2, the face {} over {} m^2",
        mass.outer_flux,
        mass.outer_area
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
            cheeks: WindowCheeks::Radial,
            ramp: None,
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
    assert_eq!(ladder.rungs.len(), 2);
    assert_eq!(ladder.order.to_bits(), (4.0f64 / 3.0).to_bits());
    let area = core::f64::consts::PI * window.bore_radius * window.bore_radius;
    let length = ladder.mass_per_density * area;
    eprintln!(
        "rungs {:.6} {:.6} mm, extrapolated {:.6} mm",
        ladder.rungs[0].bore_length * 1e3,
        ladder.rungs[1].bore_length * 1e3,
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

/// The Clare D window with parallel cheeks 8 mm apart: they meet the bore on
/// the arc `2 a asin(4 mm / a)`.
fn clare_d_parallel() -> WallWindow {
    let a = 0.0061;
    WallWindow {
        arc_width: 2.0 * a * (0.004f64 / a).asin(),
        cheeks: WindowCheeks::Parallel,
        ..clare_d()
    }
}

/// The Clare D labium as a 25 degree ramp rising from `edge_height` above
/// the bore.
fn wedge(edge_height: f64) -> Option<WindowRamp> {
    Some(WindowRamp {
        edge_height,
        angle: 25f64.to_radians(),
    })
}

#[test]
fn a_ramp_whose_edge_is_the_outer_surface_is_the_prism() {
    // The staircase path on the cut the prism path solves: its conjugate
    // gradients must land on the prism's LU solution, and its outer face's
    // cell fluxes (from the exterior's face map, or a release face's ground
    // links) on the prism's (from its diagonal channel maps).
    let release = WallWindow {
        outside: WindowOutside::Release,
        ..clare_d()
    };
    for (window, options) in [
        (clare_d(), rung(8, 12, 8)),
        (clare_d(), rung(16, 24, 16)),
        (release, rung(8, 12, 8)),
    ] {
        let prism = solve(&window, &options);
        let stair = solve(
            &WallWindow {
                ramp: wedge(0.0024),
                ..window
            },
            &options,
        );
        let largest = prism
            .outer_cells
            .iter()
            .map(|c| c.flux.abs())
            .fold(0.0, f64::max);
        assert_eq!(prism.outer_cells.len(), stair.outer_cells.len());
        for (p, s) in prism.outer_cells.iter().zip(&stair.outer_cells) {
            assert!(
                (p.theta, p.arc, p.axial, p.length) == (s.theta, s.arc, s.axial, s.length)
                    && (p.flux - s.flux).abs() < 1e-9 * largest,
                "prism cell {p:?} vs staircase {s:?}"
            );
        }
        // Against free space the flux crowds towards the outer face's
        // re-entrant edges, where the velocity is singular.
        let least = prism
            .outer_cells
            .iter()
            .map(|c| c.flux)
            .fold(f64::INFINITY, f64::min);
        eprintln!("outer cell flux: largest / least {:.4}", largest / least);
        if window.outside == WindowOutside::Free {
            assert!(largest > 1.5 * least, "edge {largest} vs centre {least}");
        }
        eprintln!(
            "prism {:.12} mm, staircase {:.12} mm after {} iterations",
            prism.bore_length * 1e3,
            stair.bore_length * 1e3,
            stair.iterations
        );
        assert_eq!(prism.iterations, 0);
        assert!(stair.iterations > 0);
        assert!(
            (stair.bore_length / prism.bore_length - 1.0).abs() < 1e-10,
            "staircase {} vs prism {}",
            stair.bore_length,
            prism.bore_length
        );
    }
}

#[test]
fn staircases_match_the_independent_dense_prototype() {
    // A numpy implementation written first, port at 40 mm and exterior at
    // 60 mm, solved by its own conjugate gradients to 1e-13: the cut cells,
    // the face potentials, and an exterior on its own grid that ends every
    // open part of the outer face on a cell edge. Spreading the outer face
    // over whole cells instead, it reproduces its direct solve of the cut
    // cells to ten digits. The release face, which has no exterior, comes
    // from that direct solve.
    let ramp = WallWindow {
        ramp: wedge(0.0005),
        ..clare_d()
    };
    let parallel = clare_d_parallel();
    let both = WallWindow {
        ramp: wedge(0.0005),
        ..clare_d_parallel()
    };
    let release = WallWindow {
        outside: WindowOutside::Release,
        ..clare_d_parallel()
    };
    for (window, options, expected_mm) in [
        (&ramp, rung(8, 12, 8), 15.1023968854),
        (&ramp, rung(16, 24, 16), 14.7738166561),
        (&parallel, rung(8, 13, 8), 18.6302048855),
        (&parallel, rung(16, 26, 16), 18.1511039750),
        (&both, rung(8, 13, 8), 17.0859619497),
        (&both, rung(16, 26, 16), 16.5541238992),
        (&release, rung(8, 13, 8), 11.7654317712),
    ] {
        let options = WallWindowOptions {
            port_radii: 0.04 / 0.0061,
            exterior_radii: 0.06 / 0.0061,
            ..options
        };
        let got = solve(window, &options).bore_length * 1e3;
        assert!(
            (got - expected_mm).abs() <= 1e-8 * expected_mm,
            "Rust {got:.10} mm vs prototype {expected_mm:.10} mm"
        );
    }
}

#[test]
fn a_grid_aligned_ramp_is_the_whole_cell_staircase() {
    // Edge at half a layer and a slope of two axial cells per layer put the
    // ramp on cell boundaries at every layer's mid-radius, so no cell is
    // cut. The numpy staircase of whole cells (open where the cell centre
    // lies inside the cut), checked against its dense twin, gives this.
    let (layer, cell) = (0.0024 / 8.0, 0.0053 / 8.0);
    let window = WallWindow {
        ramp: Some(WindowRamp {
            edge_height: 0.5 * layer,
            angle: (layer / (2.0 * cell)).atan(),
        }),
        ..clare_d()
    };
    let options = WallWindowOptions {
        port_radii: 0.04 / 0.0061,
        exterior_radii: 0.06 / 0.0061,
        ..rung(8, 12, 8)
    };
    let got = solve(&window, &options).bore_length * 1e3;
    let expected_mm = 14.6543202168;
    assert!(
        (got - expected_mm).abs() <= 1e-8 * expected_mm,
        "Rust {got:.10} mm vs whole-cell staircase {expected_mm:.10} mm"
    );
}

#[test]
fn parallel_cheeks_lie_between_the_radial_cuts_that_bound_them() {
    // On one grid, a radial cut over the window's arc contains the
    // parallel-cheeked cut, and one over the whole cells inside the cheeks'
    // angle at the outer surface, where it is narrowest, is contained in it,
    // so by Rayleigh's monotonicity their masses bracket it.
    let window = clare_d_parallel();
    let options = rung(8, 13, 8);
    let parallel = solve(&window, &options);
    let dth = 2.0 * core::f64::consts::PI / parallel.theta_cells as f64;
    let (a, t) = (window.bore_radius, window.wall_thickness);
    let half_chord = a * (6.5 * dth).sin();
    let narrowest = (2.0 * (half_chord / (a + t)).asin() / dth).floor() as usize;
    let radial = |cells: usize| {
        solve(
            &WallWindow {
                arc_width: cells as f64 * dth * a,
                ..clare_d()
            },
            &WallWindowOptions {
                arc_cells: cells,
                ..options
            },
        )
    };
    let containing = radial(13);
    let contained = radial(narrowest);
    eprintln!(
        "radial over 13 cells {:.6} mm < parallel {:.6} mm < radial over {narrowest} cells {:.6} mm",
        containing.bore_length * 1e3,
        parallel.bore_length * 1e3,
        contained.bore_length * 1e3
    );
    assert!(narrowest < 13);
    assert_eq!(containing.theta_cells, parallel.theta_cells);
    assert_eq!(contained.theta_cells, parallel.theta_cells);
    assert!(containing.bore_length < parallel.bore_length);
    assert!(parallel.bore_length < contained.bore_length);
}

#[test]
fn a_staircase_ladder_reports_the_order_it_observes() {
    let window = WallWindow {
        ramp: wedge(0.0005),
        ..clare_d_parallel()
    };
    let options = WallWindowOptions {
        workers: 4,
        ..WallWindowOptions::rung(&window, 8)
    };
    let ladder = wall_window_ladder(&window, &options).expect("staircase ladder");
    let area = core::f64::consts::PI * window.bore_radius * window.bore_radius;
    eprintln!(
        "rungs {:?} mm, observed order {:.3}, extrapolated {:.6} mm",
        ladder
            .rungs
            .iter()
            .map(|r| r.bore_length * 1e3)
            .collect::<Vec<_>>(),
        ladder.order,
        ladder.mass_per_density * area * 1e3
    );
    assert_eq!(ladder.rungs.len(), 3);
    assert!((0.5..=2.0).contains(&ladder.order));
    assert!(ladder.mass_per_density < ladder.rungs[2].mass_per_density);
}

#[test]
fn staircase_workers_do_not_change_a_bit() {
    let window = WallWindow {
        ramp: wedge(0.0005),
        ..clare_d_parallel()
    };
    let one = wall_window_mass(&window, &rung(8, 13, 8)).expect("one worker");
    let three = wall_window_mass(
        &window,
        &WallWindowOptions {
            workers: 3,
            ..rung(8, 13, 8)
        },
    )
    .expect("three workers");
    assert_eq!(one, three);
}

#[test]
fn staircase_refusals() {
    let ok = WallWindow {
        ramp: wedge(0.0005),
        ..clare_d_parallel()
    };
    let options = rung(4, 7, 4);
    assert!(wall_window_mass(&ok, &options).is_ok());
    let invalid = |window: &WallWindow, options: &WallWindowOptions| {
        matches!(
            wall_window_mass(window, options),
            Err(WallWindowError::Invalid(_))
        )
    };
    for ramp in [
        wedge(-1e-4),
        wedge(0.0025),
        wedge(f64::NAN),
        Some(WindowRamp {
            edge_height: 0.0005,
            angle: 0.0,
        }),
        Some(WindowRamp {
            edge_height: 0.0005,
            angle: 0.5 * core::f64::consts::PI,
        }),
    ] {
        assert!(invalid(&WallWindow { ramp, ..ok }, &options), "{ramp:?}");
    }
    assert!(invalid(
        &WallWindow {
            arc_width: core::f64::consts::PI * ok.bore_radius,
            ..ok
        },
        &options
    ));
    assert!(invalid(
        &ok,
        &WallWindowOptions {
            tolerance: 0.0,
            ..options
        }
    ));
    assert!(matches!(
        wall_window_mass(
            &ok,
            &WallWindowOptions {
                max_iterations: 3,
                ..options
            }
        ),
        Err(WallWindowError::Budget {
            what: "conjugate-gradient iterations",
            ..
        })
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
            what: "dense face unknowns",
            ..
        })
    ));
    assert!(matches!(
        wall_window_mass(
            &WallWindow {
                ramp: Some(WindowRamp {
                    edge_height: 0.0,
                    angle: 1e-12,
                }),
                ..ok
            },
            &options
        ),
        Err(WallWindowError::Budget {
            what: "axial cells over the outer face",
            ..
        })
    ));
    // A ramp so shallow that it runs past the exterior.
    assert!(invalid(
        &WallWindow {
            ramp: Some(WindowRamp {
                edge_height: 0.0,
                angle: 1e-3,
            }),
            ..ok
        },
        &WallWindowOptions {
            max_line_cells: 1 << 20,
            ..options
        }
    ));
}

/// The alto recorder's window (Bressan alto d-, Bouterse 2010 rows of the
/// instrument ledger): bore radius 9.8317 mm and wall 6.575 mm at the window,
/// parallel cheeks 12 mm apart, window 4.6 mm long, labium edge 0.51 mm above
/// the bore on a ramp rising 6.065 mm over 18.0 mm, flush closed end.
fn bressan_alto() -> WallWindow {
    let a = 0.0098317;
    WallWindow {
        bore_radius: a,
        wall_thickness: 0.006575,
        arc_width: 2.0 * a * (0.006f64 / a).asin(),
        axial_length: 0.0046,
        cheeks: WindowCheeks::Parallel,
        ramp: Some(WindowRamp {
            edge_height: 0.00051,
            angle: (0.006065f64 / 0.018).atan(),
        }),
        upstream: WindowUpstream::ClosedEnd { distance: 0.0023 },
        outside: WindowOutside::Free,
    }
}

/// The alto recorder's window with its labium edge at another height: the ramp
/// still meets the outer surface 18.0 mm beyond the window, as measured.
fn recorder_edge(edge_height: f64) -> WallWindow {
    let window = bressan_alto();
    WallWindow {
        ramp: Some(WindowRamp {
            edge_height,
            angle: ((window.wall_thickness - edge_height) / 0.018).atan(),
        }),
        ..window
    }
}

#[test]
#[ignore = "staircase ladders of the whistle and recorder windows to their finest rungs: a few minutes"]
fn staircase_ladders_by_rung() {
    let cases = [
        ("whistle parallel cheeks", clare_d_parallel(), 8, 4),
        (
            "whistle parallel cheeks, 0.5 mm wedge",
            WallWindow {
                ramp: wedge(0.0005),
                ..clare_d_parallel()
            },
            8,
            4,
        ),
        ("recorder", bressan_alto(), 4, 4),
        // The ends of the edge height's bracket (u-windway-height).
        ("recorder, edge at the bore", recorder_edge(0.0), 4, 4),
        ("recorder, edge 1.25 mm", recorder_edge(0.00125), 4, 4),
    ];
    for (name, window, coarse, count) in cases {
        let base = WallWindowOptions::rung(&window, coarse);
        let mut lengths = Vec::new();
        for k in 0..count {
            let options = WallWindowOptions {
                axial_cells: base.axial_cells << k,
                arc_cells: base.arc_cells << k,
                wall_layers: base.wall_layers << k,
                workers: 8,
                max_unknowns: 8192,
                ..base
            };
            let mass = solve(&window, &options);
            lengths.push(mass.bore_length);
            eprintln!(
                "{name}: rung {}/{}/{}: {:.6} mm, {} unknowns, {} iterations, outer {:.4} mm long",
                options.axial_cells,
                options.arc_cells,
                options.wall_layers,
                mass.bore_length * 1e3,
                mass.unknowns,
                mass.iterations,
                mass.outer_axial_length * 1e3
            );
            if let [.., l0, l1, l2] = lengths[..] {
                let (d1, d2) = (l0 - l1, l1 - l2);
                eprintln!(
                    "{name}: order {:.3}, extrapolated {:.6} mm",
                    (d1 / d2).log2(),
                    (l2 - d2 / (d1 / d2 - 1.0)) * 1e3
                );
            }
        }
        assert!(lengths.windows(2).all(|pair| pair[1] < pair[0]), "{name}");
    }
}
