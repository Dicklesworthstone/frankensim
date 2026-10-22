//! G1 of M2.5-W3a: an aperture matrix on an unbaffled closed tube body (two
//! wall holes plus the mouth) converges under a conforming panel ladder, and
//! the collocation and Galerkin routes converge together on a partitioned
//! sphere. G1 of M2.5-W3b realizes the same matrix as a delay-extracted,
//! shared-pole rational model. Large: run release, --test-threads=1, on
//! ts1/trj.
use fs_bem::helmholtz::aperture::{Aperture, ApertureRoute, aperture_impedance};
use fs_bem::helmholtz::{IterativeOptions, Medium};
use fs_bem::panel3d::SpherePanels;
use fs_bem::radiation_bake::lathe_profile;
use fs_bem::reciprocal_radiation::Budget;
use fs_exec::CancelGate;
use fs_math::{c64::C64, det};

/// Two disjoint caps and the rest of a unit sphere.
fn partition(surface: &SpherePanels) -> Vec<Aperture> {
    let cap = |axis: [f64; 3], cos_half_angle: f64| -> Vec<f64> {
        surface
            .centroids()
            .iter()
            .map(|c| {
                let r = det::sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]);
                let along = (c[0] * axis[0] + c[1] * axis[1] + c[2] * axis[2]) / r;
                f64::from(u8::from(along > cos_half_angle))
            })
            .collect()
    };
    let (a, b) = (cap([0.0, 0.0, 1.0], 0.85), cap([-1.0, 0.0, 0.0], 0.9));
    let rest = a.iter().zip(&b).map(|(x, y)| 1.0 - x - y).collect();
    [
        (a, [0.0, 0.0, 1.0]),
        (b, [-1.0, 0.0, 0.0]),
        (rest, [0.0, 1.0, 0.0]),
    ]
    .into_iter()
    .map(|(shape, axis)| Aperture {
        shape,
        reference_point: axis,
        normal: axis,
    })
    .collect()
}

#[test]
fn multi_aperture_routes_converge_together() {
    // rho = c = 1 on a unit sphere. The sum of every entry is the pulsating
    // sphere (the apertures partition it), exact for the discrete shapes.
    let medium = Medium {
        density: 1.0,
        sound_speed: 1.0,
    };
    let omegas = [0.3, 0.6, 0.9];
    let gate = CancelGate::new_clock_free();
    let budget = Budget {
        maximum_unknowns: 4096,
        maximum_cache_bytes: 1 << 30,
        maximum_work: 200_000_000_000,
        angular_refinement: 0,
    };
    let mut levels = Vec::new();
    for level in [2, 3] {
        let surface = SpherePanels::icosphere(1.0, level).unwrap();
        let apertures = partition(&surface);
        let area: f64 = surface.areas().iter().sum();
        let route = |route: &ApertureRoute| {
            aperture_impedance(&surface, &apertures, &omegas, medium, route, &gate).unwrap()
        };
        let galerkin = route(&ApertureRoute::Galerkin { budget });
        let collocation = route(&ApertureRoute::Collocation { iterative: None });
        let mut rows = Vec::new();
        for (g, c) in galerkin.rows.iter().zip(&collocation.rows) {
            assert!(g.reciprocity_defect <= 1e-12 && g.passivity_margin >= -1e-12);
            let ka = g.omega;
            let oracle = C64::new(ka * ka, -ka).scale(area / (1.0 + ka * ka));
            let total = |z: &[C64]| z.iter().fold(C64::ZERO, |s, &v| s + v);
            let eg = (total(&g.impedance) - oracle).abs() / oracle.abs();
            let ec = (total(&c.impedance) - oracle).abs() / oracle.abs();
            let gap: Vec<C64> = g
                .impedance
                .iter()
                .zip(&c.impedance)
                .map(|(&a, &b)| a - b)
                .collect();
            let gap = frobenius(&gap) / frobenius(&g.impedance);
            println!(
                "{{\"suite\":\"fs-bem-aperture\",\"case\":\"route-convergence\",\"level\":{level},\"panels\":{},\"ka\":{ka},\"galerkin_pulsating_rel\":{eg:.3e},\"collocation_pulsating_rel\":{ec:.3e},\"route_gap\":{gap:.3e},\"collocation_defect\":{:.3e},\"collocation_margin\":{:.3e}}}",
                surface.areas().len(),
                c.reciprocity_defect,
                c.passivity_margin
            );
            rows.push((eg, ec, gap));
        }
        levels.push(rows);
    }
    for (coarse, fine) in levels[0].iter().zip(&levels[1]) {
        assert!(
            fine.0 < coarse.0,
            "galerkin pulsating error must fall: {coarse:?} -> {fine:?}"
        );
        assert!(
            fine.1 < coarse.1,
            "collocation pulsating error must fall: {coarse:?} -> {fine:?}"
        );
        assert!(
            fine.2 < coarse.2,
            "route gap must fall: {coarse:?} -> {fine:?}"
        );
    }
}

const RADIUS: f64 = 0.015;
const LENGTH: f64 = 0.12;
/// Hole extents on the wall grid: axial [x0, x1] and angular [t0, t1] are
/// grid lines of every rung, so each hole is an exact union of panels.
const HOLES: [((f64, f64), (f64, f64)); 2] = [
    ((0.03, 0.04), (0.0, core::f64::consts::TAU / 12.0)),
    (
        (0.07, 0.08),
        (
            core::f64::consts::PI,
            core::f64::consts::PI + core::f64::consts::TAU / 12.0,
        ),
    ),
];

fn tube(rung: usize) -> (SpherePanels, Vec<Aperture>) {
    let stations = 12 << rung;
    let around = 12 << rung;
    // Equal-area faceting: each ring polygon encloses pi R^2. Inscribed rings
    // (radius R) were measured first: the mouth area ran 6.75e-4 / 6.99e-4 /
    // 7.05e-4 m^2 against pi R^2 = 7.07e-4 over 12 / 24 / 48 sides, and that
    // O(h^2) area deficit set the 500 Hz ladder steps (5.3e-2, 1.3e-2).
    let angle = core::f64::consts::TAU / around as f64;
    let radius = RADIUS * det::sqrt(angle / det::sin(angle));
    let profile: Vec<(f64, f64)> = (0..=stations)
        .map(|i| (LENGTH * i as f64 / stations as f64, radius))
        .collect();
    let body = lathe_profile(&profile, around).unwrap();
    let mouth: Vec<f64> = body
        .driven_mask()
        .iter()
        .map(|&d| f64::from(u8::from(d)))
        .collect();
    let surface = SpherePanels::from_triangles(body.triangles).unwrap();
    let mut apertures = Vec::new();
    for &((x0, x1), (t0, t1)) in &HOLES {
        let shape = surface
            .centroids()
            .iter()
            .map(|c| {
                let mut angle = det::atan2(c[2], c[1]);
                if angle < 0.0 {
                    angle += core::f64::consts::TAU;
                }
                let wall = (det::hypot(c[1], c[2]) - RADIUS).abs() < 0.05 * RADIUS;
                f64::from(u8::from(
                    wall && c[0] > x0 && c[0] < x1 && angle > t0 && angle < t1,
                ))
            })
            .collect();
        let middle = 0.5 * (t0 + t1);
        let (s, c) = det::sin_cos(middle);
        apertures.push(Aperture {
            shape,
            reference_point: [0.5 * (x0 + x1), RADIUS * c, RADIUS * s],
            normal: [0.0, c, s],
        });
    }
    apertures.push(Aperture {
        shape: mouth,
        reference_point: [LENGTH, 0.0, 0.0],
        normal: [1.0, 0.0, 0.0],
    });
    (surface, apertures)
}

fn frobenius(z: &[C64]) -> f64 {
    det::sqrt(z.iter().map(|v| v.norm_sq()).sum())
}

/// The W3a body ladder: monotone convergence at a consistent observed order,
/// with a Richardson cross-check. This is NOT the bead's 1e-3 G1 gate, which
/// the constant-panel owner does not reach at affordable rungs. trj b187
/// measured successive steps 2.07e-2 / 8.6e-3 / 3.6e-3 at 500 Hz, ratios
/// 2.39 and 2.37, about O(h^1.25): the piston velocity jumps at every
/// aperture edge.
#[test]
fn multi_aperture_body_panel_ladder_converges() {
    let medium = Medium::air();
    let omegas: Vec<f64> = [500.0, 2000.0]
        .iter()
        .map(|f| core::f64::consts::TAU * f)
        .collect();
    let gate = CancelGate::new_clock_free();
    // Workers and the M2L row cache are scheduling budgets only: every dense
    // and iterative solve is bit-identical at any worker count.
    let workers = std::thread::available_parallelism().map_or(1, |n| n.get().min(16));
    // 360 / 1488 / 5856 panels solve dense; the 23232-panel rung is above the
    // dense cap and routes to the FMM/GMRES owner. Its order-5 interpolation
    // error is far below the finest ladder step.
    let route = ApertureRoute::Collocation {
        iterative: Some(IterativeOptions {
            tol: 1e-8,
            fmm_order: 5,
            // The near field dominates this rung's matvecs (b180 profile):
            // cache it once per frequency; values are unchanged.
            near_field_cache_bytes: 8 << 30,
            workers,
            m2l_cache_bytes: 16 << 30,
            ..IterativeOptions::default()
        }),
    };
    let mut ladder = Vec::new();
    for rung in 0..4 {
        let (surface, apertures) = tube(rung);
        // Exact unions of panels; their areas converge with the faceting.
        let areas: Vec<f64> = apertures
            .iter()
            .map(|a| {
                a.shape
                    .iter()
                    .zip(surface.areas())
                    .map(|(s, a)| s * a)
                    .sum()
            })
            .collect();
        let started = std::time::Instant::now();
        let matrix = fs_bem::helmholtz::with_dense_workers(workers, || {
            aperture_impedance(&surface, &apertures, &omegas, medium, &route, &gate)
        })
        .unwrap();
        let expected = if rung < 3 { "dense" } else { "iterative" };
        assert!(matrix.rows.iter().all(|row| row.route == expected));
        for row in &matrix.rows {
            println!(
                "{{\"suite\":\"fs-bem-aperture\",\"case\":\"tube-ladder-row\",\"rung\":{rung},\"panels\":{},\"hz\":{:.1},\"route\":\"{}\",\"reciprocity_defect\":{:.3e},\"passivity_margin\":{:.3e},\"aperture_areas_m2\":{areas:?},\"wall_s\":{:.1}}}",
                surface.areas().len(),
                row.omega / core::f64::consts::TAU,
                row.route,
                row.reciprocity_defect,
                row.passivity_margin,
                started.elapsed().as_secs_f64()
            );
        }
        ladder.push(matrix);
    }
    // All rows are printed before any gate fires.
    let ladders: Vec<(Vec<f64>, Vec<f64>)> = (0..omegas.len())
        .map(|f| {
            let steps: Vec<f64> = ladder
                .windows(2)
                .map(|w| {
                    let (coarse, fine) = (&w[0].rows[f].impedance, &w[1].rows[f].impedance);
                    let gap: Vec<C64> = coarse.iter().zip(fine).map(|(&a, &b)| a - b).collect();
                    frobenius(&gap) / frobenius(fine)
                })
                .collect();
            let ratios: Vec<f64> = steps.windows(2).map(|s| s[0] / s[1]).collect();
            // Richardson with the last observed ratio from rungs (1, 2) and
            // (2, 3). Their disagreement is a model-free check of how far the
            // extrapolation itself can be trusted; nothing gates on it.
            let r = ratios[ratios.len() - 1];
            let extrapolate = |coarse: usize, fine: usize| -> Vec<C64> {
                let (c, x) = (&ladder[coarse].rows[f].impedance, &ladder[fine].rows[f].impedance);
                x.iter()
                    .zip(c)
                    .map(|(&x, &c)| x + (x - c).scale(1.0 / (r - 1.0)))
                    .collect()
            };
            let (r12, r23) = (extrapolate(1, 2), extrapolate(2, 3));
            let gap: Vec<C64> = r23.iter().zip(&r12).map(|(&a, &b)| a - b).collect();
            println!(
                "{{\"suite\":\"fs-bem-aperture\",\"case\":\"tube-ladder\",\"hz\":{:.1},\"successive_frobenius\":{steps:?},\"ratios\":{ratios:?},\"finest_error_estimate\":{:.3e},\"richardson_cross_check\":{:.3e}}}",
                omegas[f] / core::f64::consts::TAU,
                steps[steps.len() - 1] / (r - 1.0),
                frobenius(&gap) / frobenius(&r23)
            );
            (steps, ratios)
        })
        .collect();
    for ((steps, ratios), &omega) in ladders.iter().zip(&omegas) {
        assert!(
            steps.windows(2).all(|s| s[1] < s[0]),
            "ladder must converge at {omega}: {steps:?}"
        );
        // A consistent observed order: successive ratios within 15% and
        // above 2, i.e. better than first order in h.
        assert!(
            ratios.iter().all(|&q| q > 2.0) && (ratios[0] / ratios[1] - 1.0).abs() < 0.15,
            "inconsistent observed order at {omega}: {ratios:?}"
        );
    }
}

/// W3b G1: the coarse tube's aperture matrix (200 Hz to 4 kHz) as a
/// delay-extracted, shared-pole rational model. Delays are the reference
/// point distances over c. The order ladder must fall monotonically to the
/// declared 1e-3 matrix-relative operator tolerance; grid passivity of the
/// restored model is reported.
#[test]
fn mimo_rational_order_ladder() {
    use fs_vfit::FitOptions;
    use fs_vfit::mimo::{
        check_delayed_mimo_passivity, delayed_mimo_vector_fit, repair_delayed_mimo_passivity,
    };
    let (surface, apertures) = tube(0);
    let medium = Medium::air();
    let omegas: Vec<f64> = (0..120)
        .map(|i| core::f64::consts::TAU * 200.0 * det::exp(det::ln(20.0) * f64::from(i) / 119.0))
        .collect();
    let gate = CancelGate::new_clock_free();
    let matrix = aperture_impedance(
        &surface,
        &apertures,
        &omegas,
        medium,
        &ApertureRoute::Collocation { iterative: None },
        &gate,
    )
    .unwrap();
    // e^{-iwt} BEM data onto fs-vfit's s = i omega axis.
    let data: Vec<Vec<C64>> = matrix
        .rows
        .iter()
        .map(|r| r.impedance.iter().map(|z| z.conj()).collect())
        .collect();
    let n = apertures.len();
    let delays: Vec<f64> = (0..n * n)
        .map(|i| {
            let (a, b) = (
                apertures[i / n].reference_point,
                apertures[i % n].reference_point,
            );
            let d = [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
            det::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]) / medium.sound_speed
        })
        .collect();
    let mut errors = Vec::new();
    let mut finest = None;
    for order in [2, 4, 6, 8] {
        let opts = FitOptions {
            fit_d: false,
            ..FitOptions::new(order)
        };
        let (model, report) =
            delayed_mimo_vector_fit(&omegas, &data, n, &delays, &opts, true).unwrap();
        println!(
            "{{\"suite\":\"fs-bem-aperture\",\"case\":\"mimo-order-ladder\",\"order\":{order},\"max_rel\":{:.3e},\"rms_rel\":{:.3e},\"data_reciprocity\":{:.3e},\"stable\":{}}}",
            report.max_relative_error,
            report.rms_relative_error,
            report.data_reciprocity_defect,
            model.model.is_stable()
        );
        errors.push(report.max_relative_error);
        finest = Some(model);
    }
    let model = finest.unwrap();
    let band = (omegas[0], omegas[omegas.len() - 1]);
    let passivity = check_delayed_mimo_passivity(&model, band).unwrap();
    println!(
        "{{\"suite\":\"fs-bem-aperture\",\"case\":\"mimo-passivity\",\"passive\":{},\"worst\":{:.3e},\"worst_rad_s\":{:.1},\"violations\":{}}}",
        passivity.passive,
        passivity.worst.0,
        passivity.worst.1,
        passivity.violation_freqs.len()
    );
    assert!(
        errors.windows(2).all(|e| e[1] < e[0]),
        "order ladder must fall: {errors:?}"
    );
    assert!(errors[errors.len() - 1] <= 1e-3, "finest order: {errors:?}");
    // Repair moves residues only, and must stay within 0.1 x the 1e-3
    // operator tolerance over the band.
    let (repaired, repair) = repair_delayed_mimo_passivity(&model, band).unwrap();
    println!(
        "{{\"suite\":\"fs-bem-aperture\",\"case\":\"mimo-repair\",\"rounds\":{},\"relative_perturbation\":{:.3e},\"max_relative_change\":{:.3e},\"passive\":{},\"worst\":{:.3e}}}",
        repair.rounds,
        repair.relative_perturbation,
        repair.max_relative_change,
        repair.certificate.passive,
        repair.certificate.worst.0
    );
    assert!(repair.certificate.passive);
    assert_eq!(repaired.model.poles, model.model.poles);
    assert!(
        repair.max_relative_change <= 1e-4,
        "repair perturbation {:.3e} above 0.1 x tol",
        repair.max_relative_change
    );
}

/// The W3a matrix as a duct-network boundary: the coarse tube's two holes and
/// mouth terminate the plane-wave bore of the same geometry through one coupled
/// load (`fs_duct::network::CoupledTermination`). Peaks of the input impedance
/// with and without the mutual terms show the external interaction is live;
/// the even-mode identity in fs-duct checks the coupled solve itself.
#[test]
fn aperture_matrix_terminates_a_duct_network() {
    use fs_duct::network::{
        CoupledTermination, DuctNetwork, JunctionCorrections, NetworkEdge, NetworkNode,
        NodeBoundary, network_port_impedance,
    };
    use fs_duct::{Duct, LossModel, Segment};
    use fs_material::gas::{GasSpec, GasState};
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap();
    let medium = Medium::from(&gas);
    let (surface, apertures) = tube(0);
    let gate = CancelGate::new_clock_free();
    let node = |id: &str, port: bool| NetworkNode {
        id: id.into(),
        boundary: if port {
            NodeBoundary::Port { label: id.into() }
        } else {
            NodeBoundary::Junction
        },
    };
    let edge = |id: &str, from: usize, to: usize, radius: f64, length: f64| NetworkEdge {
        id: id.into(),
        from,
        to,
        duct: Duct {
            segments: vec![Segment::Cylinder { radius, length }],
        },
        corrections: JunctionCorrections::Disabled,
    };
    // Hole centres at 35 and 75 mm on a 120 mm bore; 3 mm chimneys of the
    // holes' equal-area radius.
    let hole_radius =
        det::sqrt(0.01 * RADIUS * core::f64::consts::TAU / 12.0 / core::f64::consts::PI);
    let net = DuctNetwork {
        nodes: vec![
            node("in", true),
            node("j1", false),
            node("j2", false),
            node("mouth", true),
            node("h1", true),
            node("h2", true),
        ],
        edges: vec![
            edge("bore0", 0, 1, RADIUS, 0.035),
            edge("bore1", 1, 2, RADIUS, 0.04),
            edge("bore2", 2, 3, RADIUS, 0.045),
            edge("chimney1", 1, 4, hole_radius, 0.003),
            edge("chimney2", 2, 5, hole_radius, 0.003),
        ],
    };
    // Input impedance at frequency f, with the full matrix or its diagonal.
    let z_in = |f: f64, mutual: bool| -> C64 {
        let omega = core::f64::consts::TAU * f;
        let matrix = aperture_impedance(
            &surface,
            &apertures,
            &[omega],
            medium,
            &ApertureRoute::Collocation { iterative: None },
            &gate,
        )
        .unwrap();
        let mut z = matrix.acoustic(0).unwrap();
        if !mutual {
            for (i, v) in z.iter_mut().enumerate() {
                if i / 3 != i % 3 {
                    *v = C64::ZERO;
                }
            }
        }
        let coupled = CoupledTermination {
            omega,
            labels: vec!["h1".into(), "h2".into(), "mouth".into()],
            z,
            source_id: format!("{:016x}", matrix.declaration),
        };
        network_port_impedance(&net, &gas, omega, LossModel::Bessel)
            .unwrap()
            .drive_coupled("in", &[], &coupled, &gas)
            .unwrap()
            .impedance
    };
    // Coarse log scan, then golden-section refinement of each |Z| maximum.
    let scan: Vec<f64> = (0..121)
        .map(|i| 300.0 * det::exp(det::ln(10.0) * f64::from(i) / 120.0))
        .collect();
    let values: Vec<C64> = scan.iter().map(|&f| z_in(f, true)).collect();
    for (&f, z) in scan.iter().zip(&values) {
        assert!(z.re >= 0.0, "input impedance must stay passive at {f} Hz");
    }
    let magnitudes: Vec<f64> = values.iter().map(|z| z.abs()).collect();
    // 28 golden steps shrink a ~40-cent bracket below 1e-4 cent.
    let refine = |lo: f64, hi: f64, mutual: bool| -> f64 {
        let ratio = 0.5 * (det::sqrt(5.0) - 1.0);
        let (mut a, mut b) = (lo, hi);
        for _ in 0..28 {
            let c = b - ratio * (b - a);
            let d = a + ratio * (b - a);
            if z_in(c, mutual).abs() > z_in(d, mutual).abs() {
                b = d;
            } else {
                a = c;
            }
        }
        0.5 * (a + b)
    };
    let mut shifts = Vec::new();
    for i in 1..scan.len() - 1 {
        if magnitudes[i] > magnitudes[i - 1] && magnitudes[i] > magnitudes[i + 1] {
            let coupled = refine(scan[i - 1], scan[i + 1], true);
            let diagonal = refine(scan[i - 1], scan[i + 1], false);
            let cents = 1200.0 * det::ln(coupled / diagonal) / core::f64::consts::LN_2;
            println!(
                "{{\"suite\":\"fs-bem-aperture\",\"case\":\"duct-network-peak\",\"coupled_hz\":{coupled:.3},\"self_only_hz\":{diagonal:.3},\"mutual_shift_cents\":{cents:.4}}}"
            );
            shifts.push(cents);
        }
    }
    assert!(
        shifts.len() >= 2,
        "input impedance peaks in 300-3000 Hz: {shifts:?}"
    );
    assert!(
        shifts.iter().any(|c| c.abs() > 0.01),
        "mutual radiation terms must move a resonance: {shifts:?}"
    );
}
