//! G0/G1/G2/G3/G5: open hypersingular screens and Robin Burton-Miller.
//! Run release, --test-threads=1, on ts1/trj. Strict 1e-3 gates are intentional;
//! a failed gate is a refusal of qualification, never a wider band. The disc
//! ladder gates the convergence of its sequence and its extrapolated limit at
//! 1e-3 (the 2026-10-16 .47 ruling); its finest rung keeps a declared 1% band.
use fs_bem::helmholtz::bc::{self, Admittance, BoundaryOperator, PlaneWave, SurfaceBc};
use fs_bem::helmholtz::iterative::{IterativeError, IterativeOptions, IterativeRadiation};
use fs_bem::helmholtz::{self, Formulation, HelmholtzError, Medium};
use fs_bem::panel3d::SpherePanels;
use fs_bem::reciprocal_radiation::open_screen::{OpenOptions, OpenPanels, PreparedScreen};
use fs_exec::CancelGate;
use fs_math::{c64::C64, det};
use std::f64::consts::{PI, TAU};
use std::time::Instant;

const TOL: f64 = 1e-3;
const FINEST_BAND: f64 = 1e-2;
// Independent 120/140-digit evaluation of Mellow & Karkkainen (2005),
// doi:10.1121/1.2000828, eqs.47-54, b=a. See disc_reference.py; final M24/M30
// discrepancy at ka=5 is 1.06e-30, and much smaller at the lower frequencies.
// Z/(rho*c*a^2), conjugated from the paper's exp(+iwt) convention.
const DISC: &[(f64, f64, f64)] = &[
    (0.01, 1.88634116895864947e-9, -2.66672000050792585e-2),
    (0.5, 1.27697364824787378e-2, -1.40149319530120287),
    (1.0, 2.57995671847969732e-1, -3.23312820270569601),
    (2.0, 5.79164884094304266, -6.28964477509052866),
    (5.0, 6.01693845181283926, -1.95142151494804250),
];
fn record(
    case: &str,
    provenance: &str,
    axis: &str,
    rung: usize,
    observed: C64,
    oracle: C64,
    error: f64,
    tol: f64,
    start: Instant,
    extra: &str,
) {
    println!(
        "{{\"suite\":\"fs-bem\",\"case\":\"{case}\",\"inputs\":{{\"fixture\":\"{case}\",\"provenance\":\"{provenance}\"}},\"axis\":\"{axis}\",\"rung\":{rung},\"observed\":[{:.17e},{:.17e}],\"oracle\":[{:.17e},{:.17e}],\"err\":{error:.17e},\"tol\":{tol:.17e},\"verdict\":\"{}\",\"elapsed_s\":{:.6}{extra}}}",
        observed.re,
        observed.im,
        oracle.re,
        oracle.im,
        if error <= tol { "pass" } else { "fail" },
        start.elapsed().as_secs_f64()
    );
}
fn options(level: u32) -> OpenOptions {
    OpenOptions {
        edge_level: level,
        maximum_edge_level: 4,
        triangle_order: 4,
        angular_refinement: 2,
        maximum_unknowns: 8192,
        maximum_cache_bytes: 1_000_000_000,
        maximum_work: u64::MAX,
        workers: 1,
    }
}
/// The disc ladders' options: 16 declared workers, whose results are the
/// serial owner's bits (`open_screen_workers_are_bitwise`); serial, the
/// finest rung costs about 500 s per solve.
fn ladder_options(level: u32) -> OpenOptions {
    OpenOptions {
        workers: 16,
        ..options(level)
    }
}
/// Actual mesh data, not a production disc-shaped solver. r=sin(pi*t/2)
/// gives quadratic grading of distance to the circular free edge.
fn disc(radial: usize, azimuth: usize) -> OpenPanels {
    let mut vertices = vec![[0.0; 3]];
    for j in 1..=radial {
        let r = det::sin(PI * j as f64 / (2.0 * radial as f64));
        for m in 0..azimuth {
            let phi = TAU * m as f64 / azimuth as f64;
            vertices.push([r * det::cos(phi), r * det::sin(phi), 0.0]);
        }
    }
    let at = |j: usize, m: usize| 1 + (j - 1) * azimuth + m % azimuth;
    let mut triangles = Vec::new();
    for m in 0..azimuth {
        triangles.push([0, at(1, m), at(1, m + 1)]);
    }
    for j in 1..radial {
        for m in 0..azimuth {
            triangles.push([at(j, m), at(j + 1, m), at(j + 1, m + 1)]);
            triangles.push([at(j, m), at(j + 1, m + 1), at(j, m + 1)]);
        }
    }
    OpenPanels::new(vertices, triangles).unwrap()
}
fn relative(a: C64, b: C64) -> f64 {
    (a - b).abs() / b.abs()
}
/// J1 from its power series, which needs no cancellation guard below x = 10.
fn bessel_j1(x: f64) -> f64 {
    let mut term = x / 2.0;
    let mut sum = term;
    for m in 1..40 {
        term *= -(x * x / 4.0) / (m * (m + 1)) as f64;
        sum += term;
    }
    sum
}
fn mean(surface: &SpherePanels, p: &[C64]) -> C64 {
    p.iter()
        .zip(surface.areas())
        .fold(C64::ZERO, |sum, (&p, &a)| sum + p.scale(a))
        .scale(1.0 / surface.areas().iter().sum::<f64>())
}
/// Signed Aitken limit of three rungs, exact for v_n = v + c q^n with complex
/// q, so monotone and oscillating convergence are both admitted.
fn aitken(values: &[C64]) -> (C64, C64) {
    let [a, b, c] = [values[0], values[1], values[2]];
    let q = (c - b) / (b - a);
    (c + (c - b) * q / (C64::ONE - q), q)
}
/// The ladder must contract (strictly decreasing changes, observed order in
/// [1, 3] for halved h) before its limit is compared with the oracle at TOL;
/// the finest rung is held to FINEST_BAND. Returns whether every gate held.
fn ladder_gate(case: &str, provenance: &str, ka: f64, values: &[C64], oracle: C64) -> bool {
    let start = Instant::now();
    let (limit, q) = aitken(values);
    let changes = [(values[1] - values[0]).abs(), (values[2] - values[1]).abs()];
    let order = -det::ln(q.abs()) / det::ln(2.0);
    let finest = relative(values[2], oracle);
    let error = relative(limit, oracle);
    let contracting = changes[1] < changes[0] && (1.0..=3.0).contains(&order);
    record(
        case,
        provenance,
        "edge-mesh-limit",
        2,
        limit,
        oracle,
        if contracting { error } else { f64::INFINITY },
        TOL,
        start,
        &format!(
            ",\"ka\":{ka},\"changes\":[{:.17e},{:.17e}],\"q\":[{:.17e},{:.17e}],\"observed_order\":{order},\"order_band\":[1,3],\"finest_err\":{finest:.17e},\"finest_band\":{FINEST_BAND},\"limit\":\"signed Aitken of rungs 0-2\"",
            changes[0], changes[1], q.re, q.im
        ),
    );
    contracting && error <= TOL && finest <= FINEST_BAND
}
fn bits(a: &[C64], b: &[C64]) {
    assert_eq!(a.len(), b.len());
    for (&a, &b) in a.iter().zip(b) {
        assert_eq!(
            (a.re.to_bits(), a.im.to_bits()),
            (b.re.to_bits(), b.im.to_bits())
        );
    }
}

#[test]
fn open_disc_radiation_and_scattering() {
    let medium = Medium::air();
    let rho_c = medium.density * medium.sound_speed;
    // Planted falsifier: a 2e-3 density bias must fail the limit gate.
    let biased = Medium {
        density: medium.density * 1.002,
        ..medium
    };
    let mut last = vec![f64::INFINITY; DISC.len()];
    let mut previous = last.clone();
    let mut impedances = vec![Vec::new(); DISC.len()];
    let mut sigmas = vec![Vec::new(); DISC.len()];
    let (mut masses, mut biased_masses) = (Vec::new(), Vec::new());
    for (rung, (nr, nt)) in [(4, 32), (8, 64), (16, 128)].into_iter().enumerate() {
        let surface = disc(nr, nt);
        let prepared =
            PreparedScreen::new(&surface, 5.0, ladder_options(rung as u32), || false).unwrap();
        let v = vec![C64::ONE; surface.vertices().len()];
        let scatter = vec![C64::from_re(-1.0 / rho_c); v.len()];
        for (index, &(ka, re, im)) in DISC.iter().enumerate() {
            let start = Instant::now();
            let solution = prepared
                .solve(
                    ka,
                    medium,
                    SurfaceBc {
                        drive_velocity: &v,
                        admittance: Admittance::Rigid,
                    },
                    || false,
                )
                .unwrap();
            let z = solution.force_n().scale(1.0 / rho_c);
            let oracle = C64::new(re, im);
            let error = relative(z, oracle);
            record(
                "open-disc-radiation",
                "doi:10.1121/1.2000828 eq47-54; independent Decimal series",
                "edge-mesh",
                rung,
                z,
                oracle,
                error,
                TOL,
                start,
                &format!(
                    ",\"ka\":{ka},\"panels\":{},\"unknowns\":{},\"radial\":{nr},\"azimuth\":{nt},\"grading\":2,\"quadrature\":4",
                    surface.triangles().len(),
                    prepared.unknowns()
                ),
            );
            previous[index] = last[index];
            last[index] = error;
            impedances[index].push(z);
            if ka == 0.01 {
                let mass = -solution.force_n().im / (ka * medium.sound_speed);
                let expected = (8.0 / 3.0) * medium.density;
                let err = (mass - expected).abs() / expected;
                masses.push(C64::from_re(mass));
                let force = prepared
                    .solve(
                        ka,
                        biased,
                        SurfaceBc {
                            drive_velocity: &v,
                            admittance: Admittance::Rigid,
                        },
                        || false,
                    )
                    .unwrap()
                    .force_n();
                biased_masses.push(C64::from_re(-force.im / (ka * medium.sound_speed)));
                record(
                    "open-disc-low-ka-mass",
                    "analytic: (8/3)*rho*a^3, a=1 m",
                    "edge-mesh",
                    rung,
                    C64::from_re(mass),
                    C64::from_re(expected),
                    err,
                    TOL,
                    start,
                    "",
                );
            } else {
                // Normal unit-amplitude incident wave has v_inc=1/(rho*c).
                // Solve the scattering RHS independently. Optical theorem:
                // sigma_scat/(a^2)=Re[Z/(rho*c*a^2)] for a lossless thin disc.
                let scattered = prepared
                    .solve(
                        ka,
                        medium,
                        SurfaceBc {
                            drive_velocity: &scatter,
                            admittance: Admittance::Rigid,
                        },
                        || false,
                    )
                    .unwrap();
                let sigma = 2.0 * rho_c * scattered.radiated_power_w().unwrap();
                let err = (sigma - re).abs() / re;
                sigmas[index].push(C64::from_re(sigma));
                record(
                    "open-disc-scattering",
                    "doi:10.1121/1.2000828 + optical theorem, normal incidence",
                    "edge-mesh",
                    rung,
                    C64::from_re(sigma),
                    C64::from_re(re),
                    err,
                    TOL,
                    start,
                    &format!(",\"ka\":{ka}"),
                );
                let radiation = solution.radiated_power_w().unwrap();
                let balance =
                    (radiation - solution.driving_power_w).abs() / solution.driving_power_w.abs();
                assert!(balance <= 1e-10, "open flux identity: {balance}");
            }
        }
    }
    for (&a, &b) in previous.iter().zip(&last) {
        assert!(
            b < a,
            "disc mesh ladder did not converge: {previous:?} -> {last:?}"
        );
    }
    let series = "doi:10.1121/1.2000828 eq47-54; independent Decimal series";
    let optical = "doi:10.1121/1.2000828 + optical theorem, normal incidence";
    let analytic = "analytic: (8/3)*rho*a^3, a=1 m";
    let mass = C64::from_re((8.0 / 3.0) * medium.density);
    let mut ladders = vec![("open-disc-low-ka-mass", analytic, 0.01, &masses, mass)];
    for (index, &(ka, re, im)) in DISC.iter().enumerate() {
        let z = C64::new(re, im);
        ladders.push(("open-disc-radiation", series, ka, &impedances[index], z));
        if ka != 0.01 {
            let sigma = C64::from_re(re);
            ladders.push(("open-disc-scattering", optical, ka, &sigmas[index], sigma));
        }
    }
    let refused: Vec<_> = ladders
        .into_iter()
        .filter(|&(case, source, ka, values, oracle)| {
            !ladder_gate(case, source, ka, values, oracle)
        })
        .map(|(case, _, ka, _, _)| format!("{case} ka={ka}"))
        .collect();
    assert!(refused.is_empty(), "disc ladder refused: {refused:?}");
    assert!(
        !ladder_gate("open-disc-falsifier", analytic, 0.01, &biased_masses, mass),
        "a 2e-3 density bias passed the limit gate: {biased_masses:?}"
    );
}

/// The angular chart (ANGULAR_REACH) against the rigid-disc control. kD 45 is
/// refused before any assembly, and so is kD 35 at r = 1: the chart reads the
/// requested refinement. Unguarded, kD 45 solved 24 to 29 dB high at a 1e-15
/// backward residual. At r = 2 a translating unit disc is admitted at kD 35
/// and radiates within 1.5 dB of the baffled piston on axis and 19.8 degrees
/// off it. The free disc is not the baffled piston: its edge diffraction
/// moves the ratio by about 1.1 dB even well inside the chart (kD 25 off axis
/// 1.12 dB, kD 35 on axis 1.10 dB), converged over the 16x128, 20x160, 24x160
/// and 28x144 discs to 0.04 dB, so 1.5 dB admits that physics and nothing
/// like the unresolved band.
#[test]
fn open_screen_refuses_a_band_its_angular_rule_does_not_resolve() {
    let medium = Medium::air();
    let surface = disc(16, 128);
    let prepare = |k: f64, angular_refinement: u32| {
        let options = OpenOptions {
            angular_refinement,
            ..ladder_options(0)
        };
        PreparedScreen::new(&surface, k, options, || false)
    };
    for (k, refinement) in [(22.5, 2), (17.5, 1)] {
        assert!(
            matches!(
                prepare(k, refinement),
                Err(HelmholtzError::BadParameter {
                    what: "open screen band exceeds its angular-resolution chart"
                })
            ),
            "kD {} at r = {refinement} was admitted",
            2.0 * k
        );
    }
    let start = Instant::now();
    let prepared = prepare(17.5, 2).expect("kD 35 at r = 2 is inside the chart");
    let v = vec![C64::ONE; surface.vertices().len()];
    let rigid = SurfaceBc {
        drive_velocity: &v,
        admittance: Admittance::Rigid,
    };
    let solved = prepared.solve(17.5, medium, rigid, || false).unwrap();
    let omega_rho = 17.5 * medium.sound_speed * medium.density;
    let off_axis = [0.3, 0.2, 1.0].map(|c: f64| c / 1.13_f64.sqrt());
    for (name, direction) in [("axis", [0.0, 0.0, 1.0]), ("off-axis", off_axis)] {
        // the baffled piston of unit radius: r|p| = (omega rho a^2 / 2)|2 J1(u)/u|
        let u = 17.5 * direction[0].hypot(direction[1]);
        let pattern = if u == 0.0 {
            1.0
        } else {
            2.0 * bessel_j1(u) / u
        };
        let baffled = 0.5 * omega_rho * pattern.abs();
        let free = solved.far_field(direction).unwrap().abs();
        let db = 20.0 * (free / baffled).log10();
        println!(
            "{{\"case\":\"open-angular-chart\",\"kD\":35,\"direction\":\"{name}\",\"free_over_baffled_db\":{db:.3},\"backward_residual\":{:.1e},\"seconds\":{:.1}}}",
            solved.backward_residual,
            start.elapsed().as_secs_f64()
        );
        assert!(db.abs() < 1.5, "kD 35 {name}: free/baffled {db:.2} dB");
    }
}

#[test]
fn open_disc_edge_grading_ladder() {
    let medium = Medium::air();
    let ka = 0.01;
    let mut hs = Vec::new();
    let mut errors = Vec::new();
    for (rung, (nr, nt)) in [(4, 32), (8, 64), (16, 128)].into_iter().enumerate() {
        let start = Instant::now();
        let surface = disc(nr, nt);
        let prepared =
            PreparedScreen::new(&surface, ka, ladder_options(rung as u32), || false).unwrap();
        let v = vec![C64::ONE; surface.vertices().len()];
        let solution = prepared
            .solve(
                ka,
                medium,
                SurfaceBc {
                    drive_velocity: &v,
                    admittance: Admittance::Rigid,
                },
                || false,
            )
            .unwrap();
        let mut sum = 0.0;
        let mut norm = 0.0;
        // The complete leading low-ka profile is -(4 i rho omega/pi)*sqrt(1-r^2).
        // This checks the free-edge singularity, beyond a force-only check.
        for (point, &p) in surface.vertices().iter().zip(solution.pressure_jump_pa()) {
            let r2 = point[0] * point[0] + point[1] * point[1];
            if r2 < 1.0 - 1e-12 {
                let expected = C64::new(
                    0.0,
                    -4.0 * medium.density * medium.sound_speed * ka / PI * (1.0 - r2).sqrt(),
                );
                sum += (p - expected).norm_sq();
                norm += expected.norm_sq();
            }
        }
        let error = (sum / norm).sqrt();
        hs.push(1.0 / nr as f64);
        errors.push(error);
        record(
            "open-disc-edge-profile",
            "analytic low-ka sqrt(distance) jump; ka=0.01",
            "edge-grading",
            rung,
            C64::from_re(error),
            C64::ZERO,
            error,
            TOL,
            start,
            &format!(",\"panels\":{},\"grading\":2", surface.triangles().len()),
        );
        for i in 0..surface.vertices().len() {
            if surface.is_free_edge_vertex(i) {
                assert_eq!(solution.pressure_jump_pa()[i], C64::ZERO);
            }
        }
    }
    let fit = fs_mms::fit_order(&fs_mms::RefinementLadder::new(hs, errors.clone()).unwrap());
    println!(
        "{{\"suite\":\"fs-bem\",\"case\":\"open-disc-edge-order\",\"inputs\":{{\"provenance\":\"analytic sqrt-edge jump\"}},\"axis\":\"edge-grading\",\"rung\":2,\"observed\":{},\"rms_residual\":{},\"err\":{},\"tol\":null,\"verdict\":\"{}\",\"elapsed_s\":0}}",
        fit.observed,
        fit.rms_residual,
        errors[2],
        if errors.windows(2).all(|w| w[1] < w[0]) {
            "pass"
        } else {
            "fail"
        }
    );
    assert!(
        errors.windows(2).all(|w| w[1] < w[0]),
        "edge ladder: {errors:?}"
    );
}

/// Independent spherical Bessel recurrence, not the BEM's disc/self kernels.
fn spherical(x: f64, maximum: usize) -> (Vec<f64>, Vec<f64>) {
    // Downward Miller j recurrence avoids catastrophic cancellation at ka=.1.
    let top = maximum + 48;
    let mut raw = vec![0.0; top + 2];
    raw[top] = 1.0;
    for n in (1..=top).rev() {
        raw[n - 1] = (2 * n + 1) as f64 / x * raw[n] - raw[n + 1];
    }
    let scale = det::sin(x) / x / raw[0];
    let j: Vec<_> = raw[..=maximum].iter().map(|v| v * scale).collect();
    let mut y = vec![-det::cos(x) / x, -det::cos(x) / (x * x) - det::sin(x) / x];
    for n in 1..maximum {
        y.push((2 * n + 1) as f64 / x * y[n] - y[n - 1]);
    }
    (j, y)
}
fn sphere_scattering(ka: f64, beta: C64, mu: f64) -> C64 {
    let maximum = 32;
    let (j, y) = spherical(ka, maximum + 1);
    let mut out = C64::ZERO;
    let mut previous = 1.0;
    let mut p = 1.0;
    for n in 0..=maximum {
        if n == 1 {
            p = mu;
        } else if n > 1 {
            let next = ((2 * n - 1) as f64 * mu * p - (n - 1) as f64 * previous) / n as f64;
            previous = p;
            p = next;
        }
        let derivative = if n == 0 {
            -j[1]
        } else {
            j[n - 1] - (n + 1) as f64 / ka * j[n]
        };
        let hy = if n == 0 {
            -y[1]
        } else {
            y[n - 1] - (n + 1) as f64 / ka * y[n]
        };
        let h = C64::new(j[n], y[n]);
        let hp = C64::new(derivative, hy);
        let ib = C64::new(0.0, 1.0) * beta;
        let a = (C64::from_re(derivative) + ib.scale(j[n])).scale(-1.0) / (hp + ib * h);
        out = out + a.scale((2 * n + 1) as f64 * p);
    }
    out * C64::new(0.0, -1.0 / ka)
}
fn iterative_options() -> IterativeOptions {
    IterativeOptions {
        tol: 1e-8,
        fmm_order: 6,
        max_iters: 900,
        // The 82k-panel rung's double and mixed near fields need about 2 GiB;
        // the cache is bitwise-neutral (bem_near_field_cache_is_bitwise).
        near_field_cache_bytes: 1 << 31,
        // At 81920 panels and order 6 the mixed operator counts about
        // 68720 pairs * 6^6 * 9 = 2.9e10 M2L multiply-adds per pass, past the
        // v1 8e9. The budget changes no bit; it admits that rung.
        max_translation_work: 40_000_000_000,
        ..IterativeOptions::default()
    }
}

/// The series' options: the shared iterative envelope, the second-order
/// triangle near rule (impedance_sphere_triangle_rule_ladder) and
/// bit-identical scheduling. 81920 panels hold about 3.4e6 near pairs, stored
/// as three operators each beside the block-Jacobi blocks.
fn series_options() -> IterativeOptions {
    IterativeOptions {
        panel_rule: helmholtz::PanelRule::TriangleNear { max_pairs: 1 << 23 },
        max_near_entries: 1 << 25,
        workers: 16,
        m2l_cache_bytes: 1 << 30,
        ..iterative_options()
    }
}

/// The series is ka in {0.1, 1, 5} by beta in {0, 0.5, 1, 2 + 1i}, one test
/// per case, and one per field for the ka = 5 Robin cases. Each case solves
/// alone, and in sequence they take more than 12 h, so the harness may run
/// them in parallel.
///
/// Choose `--test-threads` for memory. At 81920 panels a low-ka case may hold
/// up to `near_field_cache_bytes` (2 GiB) of cached near field beyond its FMM.
/// A ka = 5 case holds its declared directional envelopes. Each row records
/// the process's peak RSS (VmHWM), which cases run in one process share.
mod impedance_sphere_series {
    use super::*;

    /// The fields a case solves. A ka = 5 Robin case solving both takes
    /// longer than one 4 h lane (ka5_beta2_1i: 14400 s reached during the
    /// rung-6 scattering solve), so those run one test per field. Each field
    /// keeps the same gates and stops at its own first passing rung >= 5.
    #[derive(Clone, Copy, PartialEq)]
    enum Fields {
        Both,
        Pulsating,
        Scattering,
    }

    fn case(ka: f64, beta: C64) {
        case_fields(ka, beta, Fields::Both);
    }

    fn case_fields(ka: f64, beta: C64, fields: Fields) {
        let medium = Medium::air();
        let rho_c = medium.density * medium.sound_speed;
        let gate = CancelGate::new_clock_free();
        let directions: Vec<_> = (0..33)
            .map(|i| {
                let mu = -1.0 + 2.0 * i as f64 / 32.0;
                [(1.0 - mu * mu).sqrt(), 0.0, mu]
            })
            .collect();
        // Exact ka=5: use .108's explicit high-kD owner, including root-box
        // padding. Do not perturb the benchmark input to evade the W5a cap.
        // That plan stops deepening at 5120 panels. There the single layer
        // counts 2.45e7 C64 slots and 4.4e9 setup products, and the mixed
        // layer about 10x both (its ratio at 1280 panels), past the v1
        // ceilings. The declared envelopes leave about 2x over that.
        let options = if ka == 5.0 {
            let mut directional = fs_fmm::DirectionalOptions::for_tolerance(1e-4).unwrap();
            directional.max_coefficient_slots = 1 << 29;
            directional.max_setup_work = 100_000_000_000;
            // 81920 panels keep about 4.0e8 direct near pairs per operator
            // pass (the volume depth rule), past the v1 2.7e8.
            directional.max_near_pairs = 600_000_000;
            directional.setup_workers = 16;
            // At 81920 panels one directional pass counts 8.0e9 compressed
            // M2L multiply-adds (ts1, a2e2fb443, ka5_beta0 rung 6: 8000017095),
            // just past the v1 8e9.
            directional.max_translation_work = 40_000_000_000;
            // The static k = 0 rows delegate whole to W5a at this order 7.
            // At 5120 panels their mixed pass counts past 8e9 M2L
            // multiply-adds (ts1, 610ada520, ka5_beta0). Each static pass
            // runs once per solve.
            IterativeOptions {
                directional: Some(directional),
                max_translation_work: 200_000_000_000,
                ..series_options()
            }
        } else {
            series_options()
        };
        let rigid = C64::new(ka * ka, -ka).scale(1.0 / (1.0 + ka * ka));
        let oracle = (rigid / (C64::ONE + beta * rigid)).scale(rho_c);
        let mut errors = Vec::new();
        let mut scattering_errors = Vec::new();
        for rung in 3..=6 {
            let start = Instant::now();
            let surface = SpherePanels::icosphere(1.0, rung).unwrap();
            let n = surface.areas().len();
            let coefficients = vec![beta; n];
            let law = Admittance::Local {
                beta: &coefficients,
                declared_active: false,
            };
            let drive = vec![C64::ONE; n];
            let mut passed = true;
            if fields != Fields::Scattering {
                let (actual, report, receipt) = bc::solve_radiation_bc_iterative(
                    &surface,
                    ka,
                    medium,
                    SurfaceBc {
                        drive_velocity: &drive,
                        admittance: law,
                    },
                    Formulation::BurtonMiller,
                    &options,
                    &gate,
                )
                .unwrap();
                assert!(report.converged_euclidean());
                let value = mean(&surface, &actual.radiation.pressure);
                let error = relative(value, oracle);
                errors.push(error);
                record(
                    "impedance-sphere-pulsating",
                    "analytic l=0 Robin sphere: z=z0/(1+beta*z0)",
                    "sphere-mesh",
                    rung as usize,
                    value,
                    oracle,
                    error,
                    TOL,
                    start,
                    &format!(
                        ",\"rule\":\"triangle-near\",\"ka\":{ka},\"beta\":[{},{}],\"panels\":{n},\"iterations\":{},\"rel_residual\":{},\"near_field_bytes\":{},\"peak_rss_bytes\":{},\"fmm_stats\":[{},{}]",
                        beta.re,
                        beta.im,
                        report.iters,
                        report.rel_residual,
                        receipt.near_field_bytes,
                        receipt
                            .peak_rss_bytes
                            .map_or_else(|| "null".to_string(), |b| b.to_string()),
                        receipt.fmm_stats[0],
                        receipt.fmm_stats[1]
                    ),
                );
                assert!(actual.absorbed_power_w >= -1e-10);
                passed &= error <= TOL;
            }
            if fields != Fields::Pulsating {
                let (scatter, report, _) = bc::solve_scattering_bc_iterative(
                    &surface,
                    ka,
                    medium,
                    law,
                    PlaneWave {
                        direction: [0.0, 0.0, 1.0],
                        amplitude_pa: C64::ONE,
                    },
                    Formulation::BurtonMiller,
                    &options,
                    &gate,
                )
                .unwrap();
                assert!(report.converged_euclidean());
                let values = helmholtz::far_field_integrated(
                    &surface,
                    &scatter.scattered,
                    medium,
                    &directions,
                );
                let mut numerator = 0.0;
                let mut denominator = 0.0;
                for (&d, &got) in directions.iter().zip(&values) {
                    let want = sphere_scattering(ka, beta, d[2]);
                    numerator += (got - want).norm_sq();
                    denominator += want.norm_sq();
                }
                let scatter_error = (numerator / denominator).sqrt();
                scattering_errors.push(scatter_error);
                record(
                    "impedance-sphere-scattering",
                    "analytic spherical Hankel Robin partial-wave series, l<=32",
                    "sphere-mesh",
                    rung as usize,
                    C64::from_re(scatter_error),
                    C64::ZERO,
                    scatter_error,
                    TOL,
                    start,
                    &format!(
                        ",\"rule\":\"triangle-near\",\"ka\":{ka},\"beta\":[{},{}],\"panels\":{n}",
                        beta.re, beta.im
                    ),
                );
                assert!(scatter.absorbed_power_w >= -1e-10);
                passed &= scatter_error <= TOL;
            }
            // At least three spatial rungs; stop only after the strict gate
            // of every solved field.
            if rung >= 5 && passed {
                break;
            }
        }
        if fields != Fields::Scattering {
            assert!(
                *errors.last().unwrap() <= TOL,
                "ka={ka} beta={beta:?}: {errors:?}"
            );
            assert!(errors.last().unwrap() < &errors[0]);
        }
        if fields != Fields::Pulsating {
            assert!(
                *scattering_errors.last().unwrap() <= TOL,
                "ka={ka} beta={beta:?}: {scattering_errors:?}"
            );
            assert!(scattering_errors.last().unwrap() < &scattering_errors[0]);
        }
    }

    #[test]
    fn ka0_1_beta0() {
        case(0.1, C64::ZERO);
    }

    #[test]
    fn ka0_1_beta0_5() {
        case(0.1, C64::from_re(0.5));
    }

    #[test]
    fn ka0_1_beta1() {
        case(0.1, C64::ONE);
    }

    #[test]
    fn ka0_1_beta2_1i() {
        case(0.1, C64::new(2.0, 1.0));
    }

    #[test]
    fn ka1_beta0() {
        case(1.0, C64::ZERO);
    }

    #[test]
    fn ka1_beta0_5() {
        case(1.0, C64::from_re(0.5));
    }

    #[test]
    fn ka1_beta1() {
        case(1.0, C64::ONE);
    }

    #[test]
    fn ka1_beta2_1i() {
        case(1.0, C64::new(2.0, 1.0));
    }

    #[test]
    fn ka5_beta0() {
        case(5.0, C64::ZERO);
    }

    #[test]
    fn ka5_beta0_5_pulsating() {
        case_fields(5.0, C64::from_re(0.5), Fields::Pulsating);
    }

    #[test]
    fn ka5_beta0_5_scattering() {
        case_fields(5.0, C64::from_re(0.5), Fields::Scattering);
    }

    #[test]
    fn ka5_beta1_pulsating() {
        case_fields(5.0, C64::ONE, Fields::Pulsating);
    }

    #[test]
    fn ka5_beta1_scattering() {
        case_fields(5.0, C64::ONE, Fields::Scattering);
    }

    #[test]
    fn ka5_beta2_1i_pulsating() {
        case_fields(5.0, C64::new(2.0, 1.0), Fields::Pulsating);
    }

    #[test]
    fn ka5_beta2_1i_scattering() {
        case_fields(5.0, C64::new(2.0, 1.0), Fields::Scattering);
    }
}

/// G1: the declared triangle near rule raises the closed Burton-Miller rows
/// to second order. Dense icosphere rungs 2..4 at ka 1 and 3, rigid sphere;
/// ka 5 refuses rung 2 as too coarse, so its ladder is the series above. The
/// scattering far field must converge at observed order >= 1.8 under the
/// rule. The pulsating mean impedance keeps the one-point far weights' O(h)
/// remnant, so its order is recorded and it must only decrease. The centroid
/// rows run on the same meshes, and the rule must beat them at the finest
/// rung. PanelRule::Triangle (the same near rule plus a Radon seven-point far rule) is
/// held to the same claims.
#[test]
fn impedance_sphere_triangle_rule_ladder() {
    let medium = Medium::air();
    let rho_c = medium.density * medium.sound_speed;
    let directions: Vec<_> = (0..33)
        .map(|i| {
            let mu = -1.0 + 2.0 * i as f64 / 32.0;
            [(1.0 - mu * mu).sqrt(), 0.0, mu]
        })
        .collect();
    let rules = [
        ("centroid", helmholtz::PanelRule::Centroid),
        (
            "triangle-near",
            helmholtz::PanelRule::TriangleNear { max_pairs: 1 << 20 },
        ),
        (
            "triangle",
            helmholtz::PanelRule::Triangle { max_pairs: 1 << 20 },
        ),
    ];
    for ka in [1.0, 3.0] {
        let rigid = C64::new(ka * ka, -ka).scale(rho_c / (1.0 + ka * ka));
        let want: Vec<_> = directions
            .iter()
            .map(|d| sphere_scattering(ka, C64::ZERO, d[2]))
            .collect();
        let mut finest = Vec::new();
        for (label, rule) in rules {
            let (mut hs, mut pulsating, mut scattering) = (Vec::new(), Vec::new(), Vec::new());
            for rung in 2..=4 {
                let start = Instant::now();
                let surface = SpherePanels::icosphere(1.0, rung).unwrap();
                let drive = vec![C64::ONE; surface.areas().len()];
                // Rigid plane-wave scattering along +z: v = -(d.n) p_inc/(rho c).
                let incident: Vec<_> = surface
                    .centroids()
                    .iter()
                    .zip(surface.normals())
                    .map(|(x, normal)| {
                        C64::new(det::cos(ka * x[2]), det::sin(ka * x[2])).scale(-normal[2] / rho_c)
                    })
                    .collect();
                let solutions = helmholtz::with_dense_workers(16, || {
                    helmholtz::solve_radiation_batch_with_rule(
                        &surface,
                        ka,
                        medium,
                        &[&drive, &incident],
                        Formulation::BurtonMiller,
                        rule,
                    )
                })
                .unwrap();
                let value = mean(&surface, &solutions[0].pressure);
                let far =
                    helmholtz::far_field_integrated(&surface, &solutions[1], medium, &directions);
                let (numerator, denominator) = far
                    .iter()
                    .zip(&want)
                    .fold((0.0, 0.0), |(a, b), (&got, &want)| {
                        (a + (got - want).norm_sq(), b + want.norm_sq())
                    });
                let h = surface
                    .areas()
                    .iter()
                    .fold(0.0_f64, |m, &a| m.max(a))
                    .sqrt();
                hs.push(h);
                pulsating.push(relative(value, rigid));
                scattering.push((numerator / denominator).sqrt());
                println!(
                    "{{\"suite\":\"fs-bem\",\"case\":\"triangle-rule-sphere\",\"inputs\":{{\"fixture\":\"rigid icosphere\",\"provenance\":\"analytic l=0 impedance; spherical Hankel partial waves, l<=32\"}},\"axis\":\"sphere-mesh\",\"rung\":{rung},\"rule\":\"{label}\",\"ka\":{ka},\"panels\":{},\"h\":{h:.17e},\"observed\":[{:.17e},{:.17e}],\"oracle\":[{:.17e},{:.17e}],\"err\":{{\"pulsating\":{:.17e},\"scattering\":{:.17e}}},\"tol\":null,\"verdict\":\"measured\",\"elapsed_s\":{:.6}}}",
                    surface.areas().len(),
                    value.re,
                    value.im,
                    rigid.re,
                    rigid.im,
                    pulsating[pulsating.len() - 1],
                    scattering[scattering.len() - 1],
                    start.elapsed().as_secs_f64()
                );
            }
            let fit = |errors: &[f64]| {
                fs_mms::fit_order(
                    &fs_mms::RefinementLadder::new(hs.clone(), errors.to_vec()).unwrap(),
                )
            };
            let (p_fit, s_fit) = (fit(&pulsating), fit(&scattering));
            let triangle = label != "centroid";
            let held =
                !triangle || (s_fit.observed >= 1.8 && pulsating.windows(2).all(|w| w[1] < w[0]));
            println!(
                "{{\"suite\":\"fs-bem\",\"case\":\"triangle-rule-order\",\"inputs\":{{\"provenance\":\"fs_mms fit over rungs 2-4\"}},\"axis\":\"sphere-mesh\",\"rung\":4,\"rule\":\"{label}\",\"ka\":{ka},\"observed\":{{\"pulsating\":{},\"scattering\":{}}},\"rms_residual\":{{\"pulsating\":{},\"scattering\":{}}},\"err\":{:.17e},\"tol\":{},\"verdict\":\"{}\",\"elapsed_s\":0}}",
                p_fit.observed,
                s_fit.observed,
                p_fit.rms_residual,
                s_fit.rms_residual,
                scattering[2],
                if triangle {
                    "\"scattering order >= 1.8\""
                } else {
                    "null"
                },
                if !triangle {
                    "measured"
                } else if held {
                    "pass"
                } else {
                    "fail"
                }
            );
            assert!(held, "ka={ka} {label}: {pulsating:?} {scattering:?}");
            finest.push(scattering[2]);
        }
        assert!(
            finest[1..].iter().all(|&e| e < finest[0]),
            "ka={ka}: rules {finest:?}"
        );
    }
}

/// G1: the low-ka dipole radiation resistance, the damping of a rigid-body or
/// bending mode at small ka. At ka = 0.05 it is 6e-5 of |z|, so a complex
/// error that looks small can hide it. Under alpha = i/k, one-point far
/// weights leave an O(h) remnant that the resistance exposes: the centroid
/// rows overstate it about 11x, and the near rule alone stalls.
/// PanelRule::Triangle's seven-point far rule must make both its complex and its
/// resistance errors fall strictly over dense rungs 2-4, and beat both other
/// rules' resistance at the finest rung. The others are recorded.
#[test]
fn impedance_sphere_dipole_resistance_far_rule_ladder() {
    let medium = Medium::air();
    let rho_c = medium.density * medium.sound_speed;
    let ka = 0.05;
    let (j, y) = spherical(ka, 2);
    // z_1 = i rho c h_1 / h_1', with h_1' = h_0 - 2 h_1 / ka.
    let h = C64::new(j[1], y[1]);
    let hp = C64::new(j[0] - 2.0 / ka * j[1], y[0] - 2.0 / ka * y[1]);
    let exact = (C64::new(0.0, 1.0) * h / hp).scale(rho_c);
    let rules = [
        ("centroid", helmholtz::PanelRule::Centroid),
        (
            "triangle-near",
            helmholtz::PanelRule::TriangleNear { max_pairs: 1 << 20 },
        ),
        (
            "triangle",
            helmholtz::PanelRule::Triangle { max_pairs: 1 << 20 },
        ),
    ];
    let mut finest = Vec::new();
    for (label, rule) in rules {
        let (mut complex, mut resistance) = (Vec::new(), Vec::new());
        for rung in 2..=4 {
            let start = Instant::now();
            let surface = SpherePanels::icosphere(1.0, rung).unwrap();
            // The dipole (rigid translation along z): v = cos(theta).
            let v: Vec<_> = surface
                .centroids()
                .iter()
                .map(|c| C64::from_re(c[2] / (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt()))
                .collect();
            let solution = helmholtz::with_dense_workers(16, || {
                helmholtz::solve_radiation_with_rule(
                    &surface,
                    ka,
                    medium,
                    &v,
                    Formulation::BurtonMiller,
                    rule,
                )
            })
            .unwrap();
            let (pv, vv) = solution
                .pressure
                .iter()
                .zip(&v)
                .zip(surface.areas())
                .fold((C64::ZERO, 0.0), |(pv, vv), ((&p, &v), &a)| {
                    (pv + (p * v).scale(a), vv + v.re * v.re * a)
                });
            let z = pv.scale(1.0 / vv);
            complex.push(relative(z, exact));
            resistance.push((z.re - exact.re).abs() / exact.re);
            println!(
                "{{\"suite\":\"fs-bem\",\"case\":\"dipole-resistance\",\"inputs\":{{\"fixture\":\"rigid icosphere, v = cos(theta)\",\"provenance\":\"analytic z_1 = i rho c h_1/h_1'\"}},\"axis\":\"sphere-mesh\",\"rung\":{rung},\"rule\":\"{label}\",\"ka\":{ka},\"observed\":[{:.17e},{:.17e}],\"oracle\":[{:.17e},{:.17e}],\"err\":{{\"complex\":{:.17e},\"resistance\":{:.17e}}},\"tol\":null,\"verdict\":\"measured\",\"elapsed_s\":{:.6}}}",
                z.re,
                z.im,
                exact.re,
                exact.im,
                complex[complex.len() - 1],
                resistance[resistance.len() - 1],
                start.elapsed().as_secs_f64()
            );
        }
        if label == "triangle" {
            assert!(
                complex.windows(2).all(|w| w[1] < w[0]),
                "complex: {complex:?}"
            );
            assert!(
                resistance.windows(2).all(|w| w[1] < w[0]),
                "resistance: {resistance:?}"
            );
        }
        finest.push(resistance[2]);
    }
    assert!(
        finest[2] < finest[0] && finest[2] < finest[1],
        "finest resistance errors {finest:?}"
    );
}

/// G5: the iterative owner applies the dense owner's triangle near rule as
/// its sparse near rows. Sphere rungs 2-4 at k = 1 and pi agree within ten
/// times the GMRES tolerance at FMM order 7, the bound bem_iterative_matches
/// _dense holds the centroid rows to. The rule's refusals reach this owner,
/// and plain CBIE ignores the rule bit for bit.
#[test]
fn impedance_sphere_triangle_rule_iterative_matches_dense() {
    let gate = CancelGate::new_clock_free();
    let rule = helmholtz::PanelRule::TriangleNear { max_pairs: 1 << 20 };
    let options = IterativeOptions {
        tol: 1e-6,
        fmm_order: 7,
        leaf_capacity: 128,
        workers: 16,
        panel_rule: rule,
        ..IterativeOptions::default()
    };
    for rung in 2..=4 {
        let surface = SpherePanels::icosphere(1.0, rung).unwrap();
        let velocity: Vec<_> = surface
            .normals()
            .iter()
            .enumerate()
            .map(|(i, n)| C64::new(1.0 + 0.2 * n[0] + 0.03 * det::sin(i as f64), 0.1 * n[1]))
            .collect();
        for k in [1.0, PI] {
            let start = Instant::now();
            let dense = helmholtz::with_dense_workers(16, || {
                helmholtz::solve_radiation_with_rule(
                    &surface,
                    k,
                    Medium::air(),
                    &velocity,
                    Formulation::BurtonMiller,
                    rule,
                )
            })
            .unwrap();
            let (iterative, report, _) = helmholtz::solve_radiation_iterative(
                &surface,
                k,
                Medium::air(),
                &velocity,
                Formulation::BurtonMiller,
                &options,
                &gate,
            )
            .unwrap();
            assert!(report.converged_euclidean());
            let (numerator, denominator) = iterative
                .pressure
                .iter()
                .zip(&dense.pressure)
                .fold((0.0, 0.0), |(a, b), (&x, &y)| {
                    (a + (x - y).norm_sq(), b + y.norm_sq())
                });
            let error = (numerator / denominator).sqrt();
            let tol = 10.0 * options.tol;
            println!(
                "{{\"suite\":\"fs-bem\",\"case\":\"triangle-rule-iterative-dense\",\"inputs\":{{\"fixture\":\"icosphere, smooth complex velocity\",\"provenance\":\"dense LU under the same rule\"}},\"axis\":\"sphere-mesh\",\"rung\":{rung},\"k\":{k},\"panels\":{},\"iterations\":{},\"err\":{error:.17e},\"tol\":{tol:.17e},\"verdict\":\"{}\",\"elapsed_s\":{:.6}}}",
                surface.areas().len(),
                report.iters,
                if error <= tol { "pass" } else { "fail" },
                start.elapsed().as_secs_f64()
            );
            assert!(error <= tol, "rung {rung} k={k}: {error}");
        }
    }
    let surface = SpherePanels::icosphere(1.0, 2).unwrap();
    let velocity = vec![C64::ONE; surface.areas().len()];
    let solve = |formulation, o: &IterativeOptions| {
        helmholtz::solve_radiation_iterative(
            &surface,
            1.0,
            Medium::air(),
            &velocity,
            formulation,
            o,
            &gate,
        )
        .map(|(solution, _, _)| {
            solution
                .pressure
                .iter()
                .map(|p| (p.re.to_bits(), p.im.to_bits()))
                .collect::<Vec<_>>()
        })
    };
    let tight = IterativeOptions {
        panel_rule: helmholtz::PanelRule::TriangleNear { max_pairs: 1 },
        ..options.clone()
    };
    assert!(matches!(
        solve(Formulation::BurtonMiller, &tight),
        Err(IterativeError::Helmholtz(HelmholtzError::NearRuleBudget {
            maximum_pairs: 1,
            ..
        }))
    ));
    // Plain CBIE ignores the rule: one set of options can route rows of
    // either formulation, and the CBIE rows keep their own triangle table.
    let centroid = IterativeOptions {
        panel_rule: helmholtz::PanelRule::Centroid,
        ..options.clone()
    };
    assert_eq!(
        solve(Formulation::PlainCbie, &options).unwrap(),
        solve(Formulation::PlainCbie, &centroid).unwrap()
    );
    // The seven-point far cloud is W5a only.
    let directional = IterativeOptions {
        panel_rule: helmholtz::PanelRule::Triangle { max_pairs: 1 << 20 },
        directional: Some(fs_fmm::DirectionalOptions::for_tolerance(1e-4).unwrap()),
        ..options.clone()
    };
    assert!(matches!(
        solve(Formulation::BurtonMiller, &directional),
        Err(IterativeError::Invalid { .. })
    ));
}

/// G5: PanelRule::Triangle on the iterative owner, its far field on the Radon
/// seven-point cloud. It is the dense Triangle operator to FMM accuracy, at
/// rungs 2-4, ka 0.05 (where the far rule decides a dipole's resistance) and
/// ka pi (an interior resonance), within ten times the GMRES tolerance at FMM
/// order 7. The rule's exact-pair budget reaches this owner.
#[test]
fn impedance_sphere_triangle_far_cloud_matches_dense() {
    let gate = CancelGate::new_clock_free();
    let rule = helmholtz::PanelRule::Triangle { max_pairs: 1 << 20 };
    let options = IterativeOptions {
        tol: 1e-6,
        fmm_order: 7,
        leaf_capacity: 128,
        workers: 16,
        panel_rule: rule,
        // At 5120 panels and leaf capacity 128 (896 cloud points a leaf) the
        // three-operator near rows plus blocks need 18986576 entries (ts1,
        // ddcf0ceae), past 2^24; the budget changes no bit.
        max_near_entries: 1 << 25,
        ..IterativeOptions::default()
    };
    for rung in 2..=4 {
        let surface = SpherePanels::icosphere(1.0, rung).unwrap();
        // A pulsating part plus the dipole whose resistance the far rule fixes.
        let velocity: Vec<_> = surface
            .centroids()
            .iter()
            .map(|c| C64::new(1.0 + 0.5 * c[2], 0.2 * c[0]))
            .collect();
        for ka in [0.05, PI] {
            let start = Instant::now();
            let dense = helmholtz::with_dense_workers(16, || {
                helmholtz::solve_radiation_with_rule(
                    &surface,
                    ka,
                    Medium::air(),
                    &velocity,
                    Formulation::BurtonMiller,
                    rule,
                )
            })
            .unwrap();
            let (iterative, report, receipt) = helmholtz::solve_radiation_iterative(
                &surface,
                ka,
                Medium::air(),
                &velocity,
                Formulation::BurtonMiller,
                &options,
                &gate,
            )
            .unwrap();
            assert!(report.converged_euclidean());
            let (numerator, denominator) = iterative
                .pressure
                .iter()
                .zip(&dense.pressure)
                .fold((0.0, 0.0), |(a, b), (&x, &y)| {
                    (a + (x - y).norm_sq(), b + y.norm_sq())
                });
            let error = (numerator / denominator).sqrt();
            let tol = 10.0 * options.tol;
            println!(
                "{{\"suite\":\"fs-bem\",\"case\":\"triangle-far-cloud-dense\",\"inputs\":{{\"fixture\":\"icosphere, pulsating plus dipole velocity\",\"provenance\":\"dense LU under PanelRule::Triangle\"}},\"axis\":\"sphere-mesh\",\"rung\":{rung},\"ka\":{ka},\"panels\":{},\"iterations\":{},\"err\":{error:.17e},\"tol\":{tol:.17e},\"verdict\":\"{}\",\"elapsed_s\":{:.6},\"fmm_stats\":[{},{}]}}",
                surface.areas().len(),
                report.iters,
                if error <= tol { "pass" } else { "fail" },
                start.elapsed().as_secs_f64(),
                receipt.fmm_stats[0],
                receipt.fmm_stats[1]
            );
            assert!(error <= tol, "rung {rung} ka={ka}: {error}");
        }
    }
    let surface = SpherePanels::icosphere(1.0, 2).unwrap();
    let tight = IterativeOptions {
        panel_rule: helmholtz::PanelRule::Triangle { max_pairs: 1 },
        ..options
    };
    assert!(matches!(
        helmholtz::solve_radiation_iterative(
            &surface,
            1.0,
            Medium::air(),
            &vec![C64::ONE; surface.areas().len()],
            Formulation::BurtonMiller,
            &tight,
            &gate,
        ),
        Err(IterativeError::Helmholtz(HelmholtzError::NearRuleBudget {
            maximum_pairs: 1,
            ..
        }))
    ));
}

#[test]
fn impedance_bc_beta_zero_bitwise() {
    let start = Instant::now();
    let surface = SpherePanels::icosphere(1.0, 1).unwrap();
    let n = surface.areas().len();
    let drive: Vec<_> = (0..n)
        .map(|i| C64::new(1.0 + 0.2 * det::sin(i as f64), 0.1 * det::cos(i as f64)))
        .collect();
    let beta = vec![C64::ZERO; n];
    for formulation in [Formulation::PlainCbie, Formulation::BurtonMiller] {
        let legacy =
            helmholtz::solve_radiation(&surface, 1.0, Medium::air(), &drive, formulation).unwrap();
        for law in [
            Admittance::Rigid,
            Admittance::Local {
                beta: &beta,
                declared_active: false,
            },
        ] {
            let actual = bc::solve_radiation_bc(
                &surface,
                1.0,
                Medium::air(),
                SurfaceBc {
                    drive_velocity: &drive,
                    admittance: law,
                },
                formulation,
            )
            .unwrap();
            bits(&legacy.pressure, &actual.radiation.pressure);
            bits(&legacy.velocity, &actual.radiation.velocity);
            assert_eq!(
                legacy.radiated_power.to_bits(),
                actual.radiation.radiated_power.to_bits()
            );
        }
    }
    let gate = CancelGate::new_clock_free();
    let opts = iterative_options();
    let legacy = helmholtz::iterative::solve_radiation_iterative(
        &surface,
        1.0,
        Medium::air(),
        &drive,
        Formulation::BurtonMiller,
        &opts,
        &gate,
    )
    .unwrap()
    .0;
    let actual = bc::solve_radiation_bc_iterative(
        &surface,
        1.0,
        Medium::air(),
        SurfaceBc {
            drive_velocity: &drive,
            admittance: Admittance::Local {
                beta: &beta,
                declared_active: false,
            },
        },
        Formulation::BurtonMiller,
        &opts,
        &gate,
    )
    .unwrap()
    .0;
    bits(&legacy.pressure, &actual.radiation.pressure);
    bits(&legacy.velocity, &actual.radiation.velocity);
    record(
        "impedance-beta-zero",
        "legacy rigid owner",
        "none",
        0,
        C64::ZERO,
        C64::ZERO,
        0.0,
        0.0,
        start,
        "",
    );
}

/// Weighted graph Laplacian, a small explicit nonlocal test operator. It is
/// self-adjoint in the panel-area work metric: A_i Y_ij=A_j Y_ji.
struct GraphLaw {
    areas: Vec<f64>,
}
impl BoundaryOperator for GraphLaw {
    fn len(&self) -> usize {
        self.areas.len()
    }
    fn fingerprint(&self) -> u64 {
        self.areas.iter().fold(0xcbf29ce484222325, |h, v| {
            (h ^ v.to_bits()).wrapping_mul(0x100000001b3)
        })
    }
    fn apply(&self, p: &[C64], out: &mut [C64]) -> Result<(), HelmholtzError> {
        if p.len() != self.len() || out.len() != p.len() {
            return Err(HelmholtzError::ShapeMismatch {
                what: "graph admittance shape",
            });
        }
        for i in 0..p.len() {
            let previous = (i + p.len() - 1) % p.len();
            let next = (i + 1) % p.len();
            out[i] = (p[i].scale(2.0) - p[previous] - p[next]).scale(0.01 / self.areas[i]);
        }
        Ok(())
    }
}
struct MeanLaw {
    weights: Vec<f64>,
}
impl BoundaryOperator for MeanLaw {
    fn len(&self) -> usize {
        self.weights.len()
    }
    fn fingerprint(&self) -> u64 {
        self.weights.iter().fold(0x9483da11227677b9, |h, v| {
            (h ^ v.to_bits()).wrapping_mul(0x100000001b3)
        })
    }
    fn apply(&self, p: &[C64], out: &mut [C64]) -> Result<(), HelmholtzError> {
        if p.len() != self.len() || out.len() != p.len() {
            return Err(HelmholtzError::ShapeMismatch {
                what: "mean admittance shape",
            });
        }
        let mean = p
            .iter()
            .zip(&self.weights)
            .fold(C64::ZERO, |s, (&p, &w)| s + p.scale(w))
            .scale(0.5 / self.weights.iter().sum::<f64>());
        out.fill(mean);
        Ok(())
    }
}
#[test]
fn impedance_bc_passivity() {
    let start = Instant::now();
    let medium = Medium::air();
    let surface = SpherePanels::icosphere(1.0, 2).unwrap();
    let n = surface.areas().len();
    let drive: Vec<_> = surface
        .normals()
        .iter()
        .map(|n| C64::new(1.0 + 0.4 * n[2], 0.2 * n[0]))
        .collect();
    for b in [
        C64::ZERO,
        C64::from_re(0.5),
        C64::ONE,
        C64::new(2.0, 1.0),
        C64::new(0.0, 1.0),
    ] {
        let beta = vec![b; n];
        let actual = bc::solve_radiation_bc(
            &surface,
            1.0,
            medium,
            SurfaceBc {
                drive_velocity: &drive,
                admittance: Admittance::Local {
                    beta: &beta,
                    declared_active: false,
                },
            },
            Formulation::BurtonMiller,
        )
        .unwrap();
        let total = actual.radiation.radiated_power + actual.absorbed_power_w;
        assert!(total >= -1e-10 && actual.absorbed_power_w >= -1e-10);
        assert!((total - actual.driving_power_w).abs() <= 1e-10 * total.abs().max(1.0));
    }
    let graph = GraphLaw {
        areas: surface.areas().to_vec(),
    };
    let nonlocal = bc::solve_radiation_bc(
        &surface,
        1.0,
        medium,
        SurfaceBc {
            drive_velocity: &drive,
            admittance: Admittance::NonLocal {
                operator: &graph,
                declared_active: false,
            },
        },
        Formulation::BurtonMiller,
    )
    .unwrap();
    assert!(nonlocal.absorbed_power_w > 0.0 && nonlocal.radiation.radiated_power >= -1e-10);
    let gate = CancelGate::new_clock_free();
    let opts = iterative_options();
    let (iterative, report, _) = bc::solve_radiation_bc_iterative(
        &surface,
        1.0,
        medium,
        SurfaceBc {
            drive_velocity: &drive,
            admittance: Admittance::NonLocal {
                operator: &graph,
                declared_active: false,
            },
        },
        Formulation::BurtonMiller,
        &opts,
        &gate,
    )
    .unwrap();
    assert!(report.converged_euclidean());
    let relative = (iterative
        .radiation
        .pressure
        .iter()
        .zip(&nonlocal.radiation.pressure)
        .map(|(&a, &b)| (a - b).norm_sq())
        .sum::<f64>()
        / nonlocal
            .radiation
            .pressure
            .iter()
            .map(|p| p.norm_sq())
            .sum::<f64>())
    .sqrt();
    assert!(
        relative <= 1e-6,
        "dense/iterative nonlocal Robin mismatch: {relative}"
    );
    let negative = vec![C64::from_re(-1e-12); n];
    assert!(
        bc::solve_radiation_bc(
            &surface,
            1.0,
            medium,
            SurfaceBc {
                drive_velocity: &drive,
                admittance: Admittance::Local {
                    beta: &negative,
                    declared_active: false
                }
            },
            Formulation::BurtonMiller
        )
        .is_err()
    );
    let active = bc::solve_radiation_bc(
        &surface,
        1.0,
        medium,
        SurfaceBc {
            drive_velocity: &drive,
            admittance: Admittance::Local {
                beta: &negative,
                declared_active: true,
            },
        },
        Formulation::BurtonMiller,
    )
    .unwrap();
    assert!(active.declared_active);
    // Open-sheet local and nonlocal hooks, including nonuniform passive beta.
    let open = disc(3, 16);
    let p = PreparedScreen::new(&open, 1.0, options(0), || false).unwrap();
    let v = vec![C64::ONE; open.vertices().len()];
    let b: Vec<_> = open
        .vertices()
        .iter()
        .map(|x| C64::new(0.5 + 0.2 * x[0], 0.1 * x[1]))
        .collect();
    let result = p
        .solve(
            1.0,
            medium,
            SurfaceBc {
                drive_velocity: &v,
                admittance: Admittance::Local {
                    beta: &b,
                    declared_active: false,
                },
            },
            || false,
        )
        .unwrap();
    let radiated = result.radiated_power_w().unwrap();
    assert!(result.absorbed_power_w >= -1e-10);
    assert!(
        (radiated + result.absorbed_power_w - result.driving_power_w).abs()
            <= 1e-10 * result.driving_power_w.abs()
    );
    let mut weights = vec![0.0; open.vertices().len()];
    for t in open.triangles() {
        let a = open.vertices()[t[0]];
        let b = open.vertices()[t[1]];
        let c = open.vertices()[t[2]];
        let area = ((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])).abs() / 2.0;
        for &i in t {
            weights[i] += area / 3.0;
        }
    }
    let mean_law = MeanLaw { weights };
    let nonlocal = p
        .solve(
            1.0,
            medium,
            SurfaceBc {
                drive_velocity: &v,
                admittance: Admittance::NonLocal {
                    operator: &mean_law,
                    declared_active: false,
                },
            },
            || false,
        )
        .unwrap();
    assert!(nonlocal.absorbed_power_w > 0.0);
    assert!(
        (nonlocal.radiated_power_w().unwrap() + nonlocal.absorbed_power_w
            - nonlocal.driving_power_w)
            .abs()
            <= 1e-10 * nonlocal.driving_power_w.abs()
    );
    record(
        "impedance-passivity",
        "analytic weighted Robin work identity",
        "none",
        0,
        C64::ZERO,
        C64::ZERO,
        0.0,
        1e-10,
        start,
        "",
    );
}

#[test]
fn impedance_bc_fictitious_frequency_robust() {
    let start = Instant::now();
    // Level 2 (320 panels, ~10 per wavelength at ka = pi): level 1 sits at
    // 4.95, under the solver's MIN_PANELS_PER_WAVELENGTH floor of 6.
    let surface = SpherePanels::icosphere(1.0, 2).unwrap();
    let n = surface.areas().len();
    let drive = vec![C64::ONE; n];
    let gate = CancelGate::new_clock_free();
    let opts = IterativeOptions {
        leaf_capacity: n,
        ..iterative_options()
    };
    for beta in [C64::ZERO, C64::from_re(0.5), C64::new(2.0, 1.0)] {
        let b = vec![beta; n];
        let mut condition = Vec::new();
        for ka in [PI - 0.02, PI, PI + 0.02] {
            let solution = bc::solve_radiation_bc(
                &surface,
                ka,
                Medium::air(),
                SurfaceBc {
                    drive_velocity: &drive,
                    admittance: Admittance::Local {
                        beta: &b,
                        declared_active: false,
                    },
                },
                Formulation::BurtonMiller,
            )
            .unwrap();
            assert!(solution.radiation.condition_lower_bound.is_finite());
            // A lower-bound probe cannot prove good conditioning. Form the
            // actual original-equation operator by unit images, then compute
            // ||A||_1 ||A^-1||_1 with the EXISTING complex LU owner. At this
            // rung the one-leaf FMM is direct, so no far approximation enters.
            let op = IterativeRadiation::prepare_bc(
                &surface,
                ka,
                Medium::air(),
                Formulation::BurtonMiller,
                &opts,
                Admittance::Local {
                    beta: &b,
                    declared_active: false,
                },
                &gate,
            )
            .unwrap();
            let mut matrix = vec![C64::ZERO; n * n];
            let mut norm = 0.0_f64;
            for j in 0..n {
                let mut unit = vec![C64::ZERO; n];
                unit[j] = C64::ONE;
                let col = op.apply(&unit, &gate).unwrap();
                norm = norm.max(col.iter().map(|v| v.abs()).sum());
                for i in 0..n {
                    matrix[i * n + j] = col[i];
                }
            }
            let lu = fs_la::eigen_complex::lu_complex(&matrix, n).unwrap();
            let mut inverse_norm = 0.0_f64;
            for j in 0..n {
                let mut unit = vec![C64::ZERO; n];
                unit[j] = C64::ONE;
                lu.solve(&mut unit);
                inverse_norm = inverse_norm.max(unit.iter().map(|v| v.abs()).sum());
            }
            let exact_condition = norm * inverse_norm;
            assert!(
                exact_condition < 1e4,
                "Robin BM condition at ka={ka}: {exact_condition}"
            );
            condition.push(exact_condition);
        }
        assert!(condition.iter().all(|&v| v > 0.0));
        let spike = condition[1] / condition[0].max(condition[2]);
        assert!(
            spike < 4.0,
            "fictitious-frequency condition spike {condition:?}"
        );
        record(
            "impedance-fictitious-frequency",
            "analytic sphere Dirichlet ka=pi",
            "ka",
            1,
            C64::from_re(spike),
            C64::ONE,
            (spike - 1.0).abs(),
            3.0,
            start,
            "",
        );
    }
}

#[test]
fn open_and_impedance_refusals() {
    let start = Instant::now();
    assert_eq!(bc::closed_route(8192), bc::ClosedRoute::Dense);
    assert_eq!(bc::closed_route(8193), bc::ClosedRoute::Iterative);
    // Exercise the auto route with an actual over-cap closed mesh. Zero drive
    // has the exact zero solution, so this isolates routing without disguising
    // it as the nonzero continuum-accuracy sphere test above.
    let over = SpherePanels::icosphere(1.0, 5).unwrap();
    let zero = vec![C64::ZERO; over.areas().len()];
    let gate = CancelGate::new_clock_free();
    let opts = IterativeOptions {
        fmm_order: 2,
        ..IterativeOptions::default()
    };
    let run = bc::solve_radiation_bc_auto(
        &over,
        0.1,
        Medium::air(),
        SurfaceBc {
            drive_velocity: &zero,
            admittance: Admittance::Rigid,
        },
        Formulation::BurtonMiller,
        &opts,
        &gate,
    )
    .unwrap();
    assert!(run.iterative.is_some());
    assert!(
        run.solution
            .radiation
            .pressure
            .iter()
            .all(|&p| p == C64::ZERO)
    );
    let open = disc(2, 12);
    let triangles: Vec<_> = open
        .triangles()
        .iter()
        .map(|t| t.map(|i| open.vertices()[i]))
        .collect();
    let as_closed = SpherePanels::from_triangles(triangles).unwrap();
    let v = vec![C64::ONE; as_closed.areas().len()];
    assert!(
        bc::solve_radiation_bc(
            &as_closed,
            1.0,
            Medium::air(),
            SurfaceBc {
                drive_velocity: &v,
                admittance: Admittance::Rigid
            },
            Formulation::BurtonMiller
        )
        .is_err()
    );
    let sphere = SpherePanels::icosphere(1.0, 0).unwrap();
    let mut vertices = Vec::new();
    let mut faces = Vec::new();
    for t in sphere.triangles().unwrap() {
        let mut ids = [0; 3];
        for a in 0..3 {
            ids[a] = vertices.iter().position(|p| p == &t[a]).unwrap_or_else(|| {
                vertices.push(t[a]);
                vertices.len() - 1
            });
        }
        faces.push(ids);
    }
    assert!(OpenPanels::new(vertices, faces).is_err());
    let mut faces = open.triangles().to_vec();
    faces.push(faces[0]);
    assert!(OpenPanels::new(open.vertices().to_vec(), faces).is_err());
    assert!(
        OpenPanels::new(
            vec![
                [0., 0., 0.],
                [1., 0., 0.],
                [0., 1., 0.],
                [0., -1., 0.],
                [0., 0., 1.]
            ],
            vec![[0, 1, 2], [1, 0, 3], [0, 1, 4]]
        )
        .is_err()
    );
    let admitted = options(4);
    assert!(PreparedScreen::new(&open, 1.0, admitted, || false).is_ok());
    let refused = OpenOptions {
        edge_level: 5,
        ..admitted
    };
    assert!(matches!(
        PreparedScreen::new(&open, 1.0, refused, || false),
        Err(HelmholtzError::BadParameter {
            what: "open edge-grading work cap exceeded"
        })
    ));
    let area = as_closed.areas().iter().copied().fold(0.0_f64, f64::max);
    let k6 = TAU / (6.0 * area.sqrt());
    let p = PreparedScreen::new(&open, k6, options(0), || false);
    assert!(p.is_ok());
    assert!(matches!(
        PreparedScreen::new(&open, k6 * (1.0 + 1e-9), options(0), || false),
        Err(HelmholtzError::TooCoarse { .. })
    ));
    // 8192 triangles are admitted through geometry/budget checks; cancel before
    // the quadratic integration. 8193 is a named dense-route work refusal.
    for count in [8192, 8193] {
        let fan = disc(1, count);
        let mut polls = 0;
        let result = PreparedScreen::new(&fan, 1.0, options(0), || {
            polls += 1;
            polls > count + 1
        });
        if count == 8192 {
            assert!(matches!(result, Err(HelmholtzError::Cancelled)));
        } else {
            assert!(matches!(
                result,
                Err(HelmholtzError::WorkCap { panels: 8193 })
            ));
        }
    }
    let p = PreparedScreen::new(&open, 1.0, options(0), || false).unwrap();
    let byte_cap = p.cache_bytes();
    let work_cap = u64::try_from(p.frequency_work()).unwrap();
    assert!(
        PreparedScreen::new(
            &open,
            1.0,
            OpenOptions {
                maximum_cache_bytes: byte_cap,
                maximum_work: work_cap,
                ..options(0)
            },
            || false
        )
        .is_ok()
    );
    assert!(matches!(
        PreparedScreen::new(
            &open,
            1.0,
            OpenOptions {
                maximum_cache_bytes: byte_cap - 1,
                ..options(0)
            },
            || false
        ),
        Err(HelmholtzError::CacheBudget { .. })
    ));
    assert!(matches!(
        PreparedScreen::new(
            &open,
            1.0,
            OpenOptions {
                maximum_work: work_cap - 1,
                ..options(0)
            },
            || false
        ),
        Err(HelmholtzError::BadParameter {
            what: "open hypersingular frequency work cap exceeded"
        })
    ));
    let v = vec![C64::ONE; open.vertices().len()];
    let negative = vec![C64::from_re(-1e-12); v.len()];
    assert!(
        p.solve(
            1.0,
            Medium::air(),
            SurfaceBc {
                drive_velocity: &v,
                admittance: Admittance::Local {
                    beta: &negative,
                    declared_active: false
                }
            },
            || false
        )
        .is_err()
    );
    assert!(
        p.solve(
            1.0,
            Medium::air(),
            SurfaceBc {
                drive_velocity: &v,
                admittance: Admittance::Rigid
            },
            || true
        )
        .is_err()
    );
    assert!(
        p.solve(
            1.0 + 1e-9,
            Medium::air(),
            SurfaceBc {
                drive_velocity: &v,
                admittance: Admittance::Rigid
            },
            || false
        )
        .is_err()
    );
    let low_budget = OpenOptions {
        maximum_cache_bytes: 1,
        ..options(0)
    };
    assert!(matches!(
        PreparedScreen::new(&open, 1.0, low_budget, || false),
        Err(HelmholtzError::CacheBudget { .. })
    ));
    record(
        "open-impedance-refusals",
        "explicit topology, sign, work and ppw boundaries",
        "cap",
        0,
        C64::ZERO,
        C64::ZERO,
        0.0,
        0.0,
        start,
        "",
    );
}

#[test]
fn open_and_impedance_bitwise_repeat() {
    let start = Instant::now();
    let surface = disc(3, 16);
    let medium = Medium::air();
    let v = vec![C64::ONE; surface.vertices().len()];
    // The damped query below has |k| = sqrt(1 + (0.01/c)^2) > 1, so the
    // prepared band must reach past the real-frequency solves at k = 1.
    let p = PreparedScreen::new(&surface, 1.001, options(0), || false).unwrap();
    let a = p
        .solve(
            1.0,
            medium,
            SurfaceBc {
                drive_velocity: &v,
                admittance: Admittance::Rigid,
            },
            || false,
        )
        .unwrap();
    let b = p
        .solve(
            1.0,
            medium,
            SurfaceBc {
                drive_velocity: &v,
                admittance: Admittance::Rigid,
            },
            || false,
        )
        .unwrap();
    bits(a.pressure_jump_pa(), b.pressure_jump_pa());
    let fine = PreparedScreen::new(
        &surface,
        1.0,
        OpenOptions {
            triangle_order: 8,
            ..options(0)
        },
        || false,
    )
    .unwrap();
    let q = fine
        .solve(
            1.0,
            medium,
            SurfaceBc {
                drive_velocity: &v,
                admittance: Admittance::Rigid,
            },
            || false,
        )
        .unwrap();
    assert!(
        relative(q.force_n(), a.force_n()) <= TOL,
        "4/8 triangle quadrature is unresolved"
    );
    assert_eq!(
        a.radiated_power_w().unwrap().to_bits(),
        b.radiated_power_w().unwrap().to_bits()
    );
    let near = a
        .pressure_at(
            &[[0.1, 0.2, 2.0], [0.1, 0.2, -2.0]],
            1e-8,
            fs_bem::near_field::Options::default(),
            || false,
        )
        .unwrap();
    assert!(relative(near[0], near[1].scale(-1.0)) < 1e-10);
    assert!(
        a.pressure_at(
            &[[0.0; 3]],
            1e-8,
            fs_bem::near_field::Options::default(),
            || false
        )
        .is_err()
    );
    let damped = p
        .solve_laplace(
            C64::new(0.01, -medium.sound_speed),
            medium,
            SurfaceBc {
                drive_velocity: &v,
                admittance: Admittance::Rigid,
            },
            || false,
        )
        .unwrap();
    assert!(relative(damped.force_n(), a.force_n()) < 1e-3);
    let closed = SpherePanels::icosphere(1.0, 1).unwrap();
    let v = vec![C64::ONE; closed.areas().len()];
    let beta = vec![C64::new(0.5, 0.2); v.len()];
    let solve = || {
        bc::solve_radiation_bc(
            &closed,
            1.0,
            medium,
            SurfaceBc {
                drive_velocity: &v,
                admittance: Admittance::Local {
                    beta: &beta,
                    declared_active: false,
                },
            },
            Formulation::BurtonMiller,
        )
        .unwrap()
    };
    let a = solve();
    let b = solve();
    bits(&a.radiation.pressure, &b.radiation.pressure);
    bits(&a.radiation.velocity, &b.radiation.velocity);
    record(
        "open-impedance-repeat",
        "same-ISA deterministic arithmetic",
        "repeat",
        1,
        C64::ZERO,
        C64::ZERO,
        0.0,
        0.0,
        start,
        "",
    );
}

#[test]
fn open_screen_workers_are_bitwise() {
    // Declared workers assemble the pair moments, the outgoing Gram form and
    // the LU on threads; every bit must match the serial owner. Both real
    // and damped frequencies, and a row count that is not a block multiple.
    let start = Instant::now();
    let surface = disc(5, 24);
    let medium = Medium::air();
    let v: Vec<C64> = surface
        .vertices()
        .iter()
        .map(|x| C64::new(1.0 + x[0], 0.5 * x[1]))
        .collect();
    let bc = SurfaceBc {
        drive_velocity: &v,
        admittance: Admittance::Rigid,
    };
    let serial = PreparedScreen::new(&surface, 2.0, options(0), || false).unwrap();
    let a = serial.solve(2.0, medium, bc, || false).unwrap();
    let ad = serial
        .solve_laplace(
            C64::new(30.0, -1.5 * medium.sound_speed),
            medium,
            bc,
            || false,
        )
        .unwrap();
    for workers in [2, 3, 7] {
        let parallel = PreparedScreen::new(
            &surface,
            2.0,
            OpenOptions {
                workers,
                ..options(0)
            },
            || false,
        )
        .unwrap();
        let b = parallel.solve(2.0, medium, bc, || false).unwrap();
        bits(a.pressure_jump_pa(), b.pressure_jump_pa());
        assert_eq!(
            a.radiated_power_w().unwrap().to_bits(),
            b.radiated_power_w().unwrap().to_bits()
        );
        let bd = parallel
            .solve_laplace(
                C64::new(30.0, -1.5 * medium.sound_speed),
                medium,
                bc,
                || false,
            )
            .unwrap();
        bits(ad.pressure_jump_pa(), bd.pressure_jump_pa());
    }
    assert!(matches!(
        PreparedScreen::new(
            &surface,
            2.0,
            OpenOptions {
                workers: 0,
                ..options(0)
            },
            || false
        ),
        Err(HelmholtzError::BadParameter { .. })
    ));
    record(
        "open-screen-workers",
        "serial and declared-worker assembly, Gram and LU; bitwise",
        "workers",
        3,
        C64::ZERO,
        C64::ZERO,
        0.0,
        0.0,
        start,
        "",
    );
}
