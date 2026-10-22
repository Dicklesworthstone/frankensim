//! G0/G1/G2/G3/G5 directional Helmholtz tests. Native release on a quiet host.
use fs_fmm::{
    ComplexBudget, DirectionalOptions, FmmC64, FmmError, HelmholtzAdjointDouble, HelmholtzDouble,
    HelmholtzHypersingular, HelmholtzSingle, MAX_COMPLEX_COEFFICIENT_SLOTS,
    MAX_DECLARED_COEFFICIENT_SLOTS, MAX_DIRECTIONAL_CONES, MAX_INTERPOLATION_ORDER,
};
use fs_math::{c64::C64, det};
use std::{sync::Once, time::Instant};

const SUITE: &str = "fs-fmm/directional";
const SEED: u64 = 0x1080_2026_1016_0001;
fn emit(case: &str, json: String) {
    static RUN: Once = Once::new();
    RUN.call_once(|| {
        let sha=std::process::Command::new("git").args(["rev-parse","HEAD"]).output().unwrap();
        let dirty=std::process::Command::new("git").args(["status","--porcelain","--untracked-files=no"]).output().unwrap();
        assert!(sha.status.success()&&dirty.status.success());
        println!("{{\"suite\":\"{SUITE}\",\"case\":\"run\",\"git_sha\":\"{}\",\"dirty\":{},\"profile\":\"{}\",\"host\":{:?},\"cpu\":{:?},\"command\":{:?}}}",String::from_utf8(sha.stdout).unwrap().trim(),!dirty.stdout.is_empty(),if cfg!(debug_assertions){"debug"}else{"release"},std::fs::read_to_string("/proc/sys/kernel/hostname").unwrap_or_default().trim(), std::fs::read_to_string("/proc/cpuinfo").unwrap_or_default().lines().find(|s|s.starts_with("model name")).unwrap_or("unavailable"),std::env::args().collect::<Vec<_>>().join(" "));
    });
    let mut emitter = fs_obs::Emitter::new(SUITE, case);
    let event = emitter.emit(
        fs_obs::Severity::Info,
        fs_obs::EventKind::Custom {
            name: case.into(),
            json,
        },
        None,
    );
    fs_obs::lint_failure_record(&event).unwrap();
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).unwrap();
    println!("{line}");
}
fn rel(a: &[C64], b: &[C64]) -> f64 {
    let num = a
        .iter()
        .zip(b)
        .map(|(&a, &b)| (a - b).norm_sq())
        .sum::<f64>();
    let den = b.iter().map(|&v| v.norm_sq()).sum::<f64>();
    assert!(den > 0.0);
    det::sqrt(num / den)
}
fn bits(a: &[C64]) -> Vec<(u64, u64)> {
    a.iter().map(|v| (v.re.to_bits(), v.im.to_bits())).collect()
}
fn side(points: &[[f64; 3]]) -> f64 {
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for p in points {
        for a in 0..3 {
            lo[a] = lo[a].min(p[a]);
            hi[a] = hi[a].max(p[a]);
        }
    }
    (hi[0] - lo[0])
        .max(hi[1] - lo[1])
        .max(hi[2] - lo[2])
        .max(1e-12)
        * (1.0 + 1e-9)
}
fn charges(n: usize) -> Vec<C64> {
    let mut rng = SEED;
    (0..n)
        .map(|_| {
            rng = rng
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let x = (rng >> 11) as f64 / (1u64 << 53) as f64;
            rng = rng
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            C64::new(x - 0.5, (rng >> 11) as f64 / (1u64 << 53) as f64 - 0.5)
        })
        .collect()
}
fn fixture(n: usize) -> (Vec<[f64; 3]>, Vec<[f64; 3]>) {
    // Two compact 3D clusters. High-frequency M2L remains active at every kD;
    // far-only checks below prevent exact intra-cluster P2P hiding its error.
    let mut points = Vec::new();
    let mut normals = Vec::new();
    for i in 0..n {
        let j = i % (n / 2);
        points.push([
            if i < n / 2 { 0.0 } else { 1.0 } + 0.007 * j as f64 / (n / 2) as f64,
            0.007 * ((j * 3) % 11) as f64 / 11.0,
            0.007 * ((j * 5) % 13) as f64 / 13.0,
        ]);
        normals.push(match i % 4 {
            0 => [0.6, 0.0, 0.8],
            1 => [0.0, -0.8, 0.6],
            2 => [1.0, 0.0, 0.0],
            _ => [0.0, 1.0, 0.0],
        });
    }
    (points, normals)
}
fn options(eps: f64, rung: usize) -> DirectionalOptions {
    let mut o = DirectionalOptions::for_tolerance(eps).unwrap();
    o.low_freq_threshold = 1.0;
    o.cone_resolution = rung + 1;
    o.leaf_capacity = 16;
    o
}
fn operator<'a>(
    kind: usize,
    kernel: &'a HelmholtzSingle,
    points: &[[f64; 3]],
    normals: &[[f64; 3]],
    o: DirectionalOptions,
) -> FmmC64<'a> {
    match kind {
        0 => FmmC64::new_directional(kernel, points, o).unwrap(),
        1 => FmmC64::new_double_directional(
            HelmholtzDouble { k: kernel.k },
            points.to_vec(),
            normals.to_vec(),
            o,
        )
        .unwrap(),
        2 => FmmC64::new_adjoint_double_directional(
            HelmholtzAdjointDouble { k: kernel.k },
            points.to_vec(),
            normals.to_vec(),
            o,
        )
        .unwrap(),
        _ => FmmC64::new_hypersingular_directional(
            HelmholtzHypersingular { k: kernel.k },
            points.to_vec(),
            normals.to_vec(),
            o,
        )
        .unwrap(),
    }
}

#[test]
fn dfmm_high_kd_tolerance_ladder() {
    let (points, normals) = fixture(16);
    let mut q = charges(points.len());
    for v in &mut q[points.len() / 2..] {
        *v = C64::ZERO;
    }
    for kd in [20.0, 50.0, 130.0, 370.0] {
        for kind in 0..4 {
            let kernel = HelmholtzSingle {
                k: kd / side(&points),
            };
            let mut previous = f64::INFINITY;
            for (rung, eps) in [1e-3, 1e-5, 1e-7].into_iter().enumerate() {
                let started = Instant::now();
                let o = options(eps, rung);
                let op = operator(kind, &kernel, &points, &normals, o.clone());
                let oracle = op.direct(&q).unwrap();
                let got = op.potentials(&q).unwrap();
                let full = rel(&got, &oracle);
                let far = rel(&got[8..], &oracle[8..]);
                let error = full.max(far);
                let stats = op.directional_stats().unwrap();
                emit(
                    "dfmm-high-kd",
                    format!(
                        "{{\"case\":\"dfmm-high-kd\",\"inputs\":\"analytic outgoing Helmholtz, two 3D clusters, source charges on first cluster\",\"seed\":{SEED},\"axis\":\"cone-and-order\",\"rung\":{rung},\"kind\":{kind},\"n\":16,\"kd\":{kd},\"eps\":{eps},\"order\":{},\"observed\":{error},\"err\":{error},\"rel_l2\":{full},\"far_rel_l2\":{far},\"oracle\":\"direct point sum\",\"tol\":{eps},\"verdict\":\"{}\",\"elapsed_s\":{},\"stats\":{}}}",
                        o.order,
                        if error <= eps && error < previous {
                            "pass"
                        } else {
                            "fail"
                        },
                        started.elapsed().as_secs_f64(),
                        op.stats()
                    ),
                );
                assert!(
                    stats.high_frequency_pairs > 0 && stats.m2l_pairs > 0,
                    "directional factor must be live"
                );
                assert!(error <= eps, "kD={kd}, kind={kind}, eps={eps}: {error}");
                assert!(
                    error < previous,
                    "cone/order error not decreasing at kD={kd}, kind={kind}: {previous} -> {error}"
                );
                previous = error;
            }
        }
    }
    // A distributed 3D geometry also exercises the scalar interpolation. The
    // small N keeps the p12 unfused translations within the unchanged cap.
    let points = sphere(12);
    let q = charges(points.len());
    for kd in [20.0, 50.0, 130.0, 370.0] {
        let kernel = HelmholtzSingle {
            k: kd / side(&points),
        };
        let mut previous = f64::INFINITY;
        for (rung, eps) in [1e-3, 1e-5, 1e-7].into_iter().enumerate() {
            let started = Instant::now();
            let op = FmmC64::new_directional(&kernel, &points, options(eps, rung)).unwrap();
            let oracle = op.direct(&q).unwrap();
            let got = op.potentials(&q).unwrap();
            let error = rel(&got, &oracle);
            emit(
                "dfmm-sphere-high-kd",
                format!(
                    "{{\"case\":\"dfmm-sphere-high-kd\",\"inputs\":\"analytic Fibonacci sphere, complex seeded charges, direct point oracle\",\"seed\":{SEED},\"axis\":\"cone-and-order\",\"rung\":{rung},\"n\":12,\"kd\":{kd},\"observed\":{error},\"err\":{error},\"tol\":{eps},\"verdict\":\"{}\",\"elapsed_s\":{},\"stats\":{}}}",
                    if error <= eps && error < previous {
                        "pass"
                    } else {
                        "fail"
                    },
                    started.elapsed().as_secs_f64(),
                    op.stats()
                ),
            );
            assert!(
                error <= eps && error < previous,
                "distributed sphere kD={kd}, eps={eps}: {previous} -> {error}"
            );
            previous = error;
        }
    }
    // Exercise the >MAX_DIRECT_POINTS branch with the prescribed 1000-target
    // independent direct oracle. Uniform cloud and localized-source checks
    // above are distinct; this sparse two-cluster large case is not a timing
    // or generic geometry-complexity qualification.
    let (points, _) = fixture(16_400);
    let q = charges(points.len());
    let kd = 130.0;
    let kernel = HelmholtzSingle {
        k: kd / side(&points),
    };
    let o = options(1e-5, 1);
    let started = Instant::now();
    let op = FmmC64::new_directional(&kernel, &points, o).unwrap();
    assert!(matches!(
        op.direct(&q),
        Err(FmmError::WorkEnvelopeExceeded { .. })
    ));
    let targets: Vec<_> = (0..1000).map(|j| j * points.len() / 1000).collect();
    let oracle = op.direct_at_targets(&q, &targets).unwrap();
    let got = op.potentials(&q).unwrap();
    let sampled: Vec<_> = targets.iter().map(|&j| got[j]).collect();
    let error = rel(&sampled, &oracle);
    emit(
        "dfmm-sampled-direct",
        format!(
            "{{\"case\":\"dfmm-sampled-direct\",\"inputs\":\"analytic outgoing kernel, 16400 clustered points, 1000 fixed-stride targets\",\"axis\":\"target-count\",\"rung\":1000,\"n\":16400,\"kd\":{kd},\"observed\":{error},\"err\":{error},\"tol\":0.00001,\"verdict\":\"{}\",\"elapsed_s\":{},\"stats\":{}}}",
            if error <= 1e-5 { "pass" } else { "fail" },
            started.elapsed().as_secs_f64(),
            op.stats()
        ),
    );
    assert!(error <= 1e-5);
}

fn sphere(n: usize) -> Vec<[f64; 3]> {
    (0..n)
        .map(|i| {
            let z = 1.0 - 2.0 * (i as f64 + 0.5) / n as f64;
            let r = det::sqrt(1.0 - z * z);
            let (s, c) = det::sin_cos(i as f64 * core::f64::consts::PI * (3.0 - det::sqrt(5.0)));
            [r * c, r * s, z]
        })
        .collect()
}
fn exponent(x: &[f64], y: &[f64]) -> f64 {
    let lx: Vec<_> = x.iter().map(|&v| det::ln(v)).collect();
    let ly: Vec<_> = y.iter().map(|&v| det::ln(v)).collect();
    let mx = lx.iter().sum::<f64>() / lx.len() as f64;
    let my = ly.iter().sum::<f64>() / ly.len() as f64;
    lx.iter()
        .zip(&ly)
        .map(|(&a, &b)| (a - mx) * (b - my))
        .sum::<f64>()
        / lx.iter().map(|&a| (a - mx) * (a - mx)).sum::<f64>()
}
#[test]
fn dfmm_work_scaling() {
    let (mut sizes, mut times, mut works) = (Vec::new(), Vec::new(), Vec::new());
    let (mut total_times, mut total_works) = (Vec::new(), Vec::new());
    for (rung, n) in [4096, 8192, 16384, 32768].into_iter().enumerate() {
        // Radius, ring resolution, axial spacing and physical k stay fixed;
        // length grows. Both surface directions keep fixed points/wavelength.
        // This finite family is the declared scaling domain, not an isotropic
        // large-sphere or arbitrary aspect-ratio performance claim.
        let points: Vec<_> = (0..n)
            .map(|j| {
                let (s, c) = det::sin_cos(core::f64::consts::TAU * (j % 8) as f64 / 8.0);
                [0.02 * c, 0.02 * s, 0.125 * (j / 8) as f64]
            })
            .collect();
        let q = charges(n);
        let kd = 5.0 * side(&points);
        let kernel = HelmholtzSingle { k: 5.0 };
        let mut o = DirectionalOptions::for_tolerance(1e-3).unwrap();
        o.leaf_capacity = 16;
        let setup = Instant::now();
        let op = FmmC64::new_directional(&kernel, &points, o).unwrap();
        let setup_s = setup.elapsed().as_secs_f64();
        // Three native matvec timings, median, with no timing acceptance fit.
        let mut samples = Vec::new();
        let mut got = Vec::new();
        for _ in 0..3 {
            let t = Instant::now();
            got = op.potentials(&q).unwrap();
            samples.push(t.elapsed().as_secs_f64());
        }
        samples.sort_by(f64::total_cmp);
        let elapsed = samples[1];
        let error = if n <= fs_fmm::MAX_DIRECT_POINTS {
            rel(&got, &op.direct(&q).unwrap())
        } else {
            let targets: Vec<_> = (0..1000).map(|i| i * n / 1000).collect();
            let sampled: Vec<_> = targets.iter().map(|&i| got[i]).collect();
            rel(&sampled, &op.direct_at_targets(&q, &targets).unwrap())
        };
        let stats = op.directional_stats().unwrap();
        assert!(
            stats.maximum_rank > 0,
            "scaling must exercise compressed coupling blocks"
        );
        emit(
            "dfmm-work-scaling",
            format!(
                "{{\"case\":\"dfmm-work-scaling\",\"inputs\":\"analytic cylindrical surface, radius .02m, 8 points/ring, .125m axial spacing, k=5/m, complex seeded charges\",\"seed\":{SEED},\"axis\":\"N-fixed-points-per-wavelength\",\"rung\":{rung},\"n\":{n},\"kd\":{kd},\"observed\":{error},\"err\":{error},\"tol\":0.001,\"verdict\":\"{}\",\"elapsed_s\":{elapsed},\"setup_s\":{setup_s},\"stats\":{}}}",
                if error <= 1e-3 { "pass" } else { "fail" },
                op.stats()
            ),
        );
        assert!(error <= 1e-3);
        sizes.push(n as f64);
        times.push(elapsed);
        works.push(stats.work() as f64);
        total_times.push(setup_s + elapsed);
        total_works.push((stats.setup_work + stats.work()) as f64);
    }
    let time_p = exponent(&sizes, &times);
    let work_p = exponent(&sizes, &works);
    let total_time_p = exponent(&sizes, &total_times);
    let total_work_p = exponent(&sizes, &total_works);
    emit(
        "dfmm-scaling-exponent",
        format!(
            "{{\"case\":\"dfmm-scaling-exponent\",\"inputs\":\"four native quiet-host rungs, fixed surface sampling per wavelength\",\"axis\":\"least-squares-log-log\",\"rung\":4,\"observed\":{time_p},\"err\":{work_p},\"total_time_exponent\":{total_time_p},\"total_work_exponent\":{total_work_p},\"oracle\":\"matvec and preparation+matvec wall/work exponents <=1.3\",\"tol\":1.3,\"verdict\":\"{}\",\"elapsed_s\":{}}}",
            if time_p <= 1.3 && work_p <= 1.3 && total_time_p <= 1.3 && total_work_p <= 1.3 {
                "pass"
            } else {
                "fail"
            },
            times.iter().sum::<f64>()
        ),
    );
    assert!(
        time_p <= 1.3 && work_p <= 1.3 && total_time_p <= 1.3 && total_work_p <= 1.3,
        "apply wall/work {time_p}/{work_p}; setup+apply wall/work {total_time_p}/{total_work_p}; times={times:?}, work={works:?}"
    );
}

#[test]
fn dfmm_matches_bbfmm_in_overlap() {
    let (points, normals) = fixture(16);
    let q = charges(points.len());
    for kd in [1.0, 5.0, 10.0] {
        for eps in [1e-3, 1e-5, 1e-7] {
            for kind in 0..4 {
                let kernel = HelmholtzSingle {
                    k: kd / side(&points),
                };
                let o = options(eps, 0);
                let directional = operator(kind, &kernel, &points, &normals, o.clone());
                let bb = match kind {
                    0 => FmmC64::new(&kernel, points.clone(), 12, 1).unwrap(),
                    1 => FmmC64::new_double(
                        HelmholtzDouble { k: kernel.k },
                        points.clone(),
                        normals.clone(),
                        12,
                        1,
                    )
                    .unwrap(),
                    2 => FmmC64::new_adjoint_double(
                        HelmholtzAdjointDouble { k: kernel.k },
                        points.clone(),
                        normals.clone(),
                        12,
                        1,
                    )
                    .unwrap(),
                    _ => FmmC64::new_hypersingular(
                        HelmholtzHypersingular { k: kernel.k },
                        points.clone(),
                        normals.clone(),
                        12,
                        1,
                    )
                    .unwrap(),
                };
                let started = Instant::now();
                let oracle = bb.direct(&q).unwrap();
                let a = bb.potentials(&q).unwrap();
                let b = directional.potentials(&q).unwrap();
                let ea = rel(&a, &oracle);
                let eb = rel(&b, &oracle);
                let agreement = rel(&a, &b);
                emit(
                    "dfmm-overlap",
                    format!(
                        "{{\"case\":\"dfmm-overlap\",\"inputs\":\"analytic two-cluster Helmholtz, W5a p12 versus directional\",\"axis\":\"kD\",\"rung\":{kd},\"kind\":{kind},\"eps\":{eps},\"observed\":{agreement},\"err\":{},\"tol\":{},\"verdict\":\"{}\",\"elapsed_s\":{}}}",
                        ea.max(eb),
                        2.0 * eps,
                        if ea <= eps && eb <= eps && agreement <= 2.0 * eps {
                            "pass"
                        } else {
                            "fail"
                        },
                        started.elapsed().as_secs_f64()
                    ),
                );
                assert!(
                    ea <= eps && eb <= eps && agreement <= 2.0 * eps,
                    "overlap kD={kd} kind={kind} eps={eps}: {ea}, {eb}, {agreement}"
                );
            }
        }
    }
}

#[test]
fn dfmm_gradient_tolerance_ladder() {
    let (points, _) = fixture(16);
    let mut q = charges(points.len());
    q[8..].fill(C64::ZERO);
    for kd in [20.0, 50.0, 130.0, 370.0] {
        for axis in 0..3 {
            let kernel = fs_fmm::HelmholtzGradient::new(kd / side(&points), axis).unwrap();
            let mut previous = f64::INFINITY;
            for (rung, eps) in [1e-3, 1e-5, 1e-7].into_iter().enumerate() {
                let started = Instant::now();
                let op = FmmC64::new_directional(&kernel, &points, options(eps, rung)).unwrap();
                let oracle = op.direct(&q).unwrap();
                let actual = op.potentials(&q).unwrap();
                assert!(
                    actual
                        .iter()
                        .chain(&oracle)
                        .all(|z| z.re.is_finite() && z.im.is_finite()),
                    "nonfinite gradient kD={kd}, axis={axis}, eps={eps}"
                );
                let error = rel(&actual, &oracle).max(rel(&actual[8..], &oracle[8..]));
                emit(
                    "dfmm-gradient",
                    format!(
                        "{{\"case\":\"dfmm-gradient\",\"inputs\":\"analytic Cartesian outgoing Green derivative, far-only and full point oracle\",\"axis\":\"cone-and-order\",\"rung\":{rung},\"component\":{axis},\"kd\":{kd},\"observed\":{error},\"err\":{error},\"tol\":{eps},\"verdict\":\"{}\",\"elapsed_s\":{}}}",
                        if error <= eps && error < previous {
                            "pass"
                        } else {
                            "fail"
                        },
                        started.elapsed().as_secs_f64()
                    ),
                );
                assert!(
                    error <= eps && error < previous,
                    "gradient kD={kd}, axis={axis}, eps={eps}: {previous} -> {error}"
                );
                previous = error;
            }
        }
    }
}

#[test]
fn dfmm_transverse_gradient_retains_coordinate_factor() {
    // G1/G3: the original p5/kD50 failure and an axis permutation plus a
    // translation. The untagged, identical point kernel plants the old
    // interpolate-the-entire-component operator as an independent control.
    struct Unfactored(fs_fmm::HelmholtzGradient);
    impl fs_fmm::ComplexKernel for Unfactored {
        fn eval(&self, x: [f64; 3], y: [f64; 3]) -> C64 {
            fs_fmm::ComplexKernel::eval(&self.0, x, y)
        }
        fn wavenumber(&self) -> Option<f64> {
            fs_fmm::ComplexKernel::wavenumber(&self.0)
        }
        fn translation_invariant(&self) -> bool {
            true
        }
    }
    let (points, _) = fixture(16);
    let mut q = charges(points.len());
    q[8..].fill(C64::ZERO);
    for moved in [false, true] {
        let started = Instant::now();
        let (points, axis) = if moved {
            (
                points
                    .iter()
                    .map(|p| [p[2] + 8.0, p[0] - 4.0, p[1] + 2.0])
                    .collect::<Vec<_>>(),
                0,
            )
        } else {
            (points.clone(), 2)
        };
        let kernel = fs_fmm::HelmholtzGradient::new(50.0 / side(&points), axis).unwrap();
        let control = Unfactored(kernel);
        let o = options(1e-3, 0);
        let factored = FmmC64::new_directional(&kernel, &points, o.clone()).unwrap();
        let unfactored = FmmC64::new_directional(&control, &points, o.clone()).unwrap();
        let oracle = unfactored.direct(&q).unwrap();
        assert_eq!(bits(&factored.direct(&q).unwrap()), bits(&oracle));
        let actual = factored.potentials(&q).unwrap();
        let old = unfactored.potentials(&q).unwrap();
        let error = rel(&actual[8..], &oracle[8..]);
        let old_error = rel(&old[8..], &oracle[8..]);
        emit(
            "dfmm-gradient-coordinate-factor",
            format!(
                "{{\"case\":\"dfmm-gradient-coordinate-factor\",\"inputs\":\"analytic outgoing gradient; identical unfactored point kernel control\",\"seed\":{SEED},\"axis\":\"coordinate-factor\",\"rung\":{},\"moved\":{moved},\"kd\":50,\"component\":{axis},\"order\":{},\"err\":{error},\"unfactored_err\":{old_error},\"tol\":{},\"verdict\":\"{}\",\"elapsed_s\":{}}}",
                usize::from(moved),
                o.order,
                o.tolerance,
                if error <= o.tolerance && old_error > o.tolerance {
                    "pass"
                } else {
                    "fail"
                },
                started.elapsed().as_secs_f64()
            ),
        );
        assert!(
            actual
                .iter()
                .chain(&old)
                .all(|z| z.re.is_finite() && z.im.is_finite())
        );
        assert!(
            error <= o.tolerance,
            "factored gradient moved={moved}: {error}"
        );
        assert!(
            old_error > o.tolerance,
            "unfactored control must fail: {old_error}"
        );
        let stats = factored.directional_stats().unwrap();
        assert!(stats.high_frequency_pairs > 0);
        assert_eq!(stats.order, o.order);
        assert_eq!(
            stats.cones_per_level,
            unfactored.directional_stats().unwrap().cones_per_level
        );
    }
}

#[test]
fn dfmm_coplanar_transverse_gradient_is_exactly_zero() {
    // G0/G3: the affine factor vanishes identically. Neither interpolation
    // nor a translation of the plane may manufacture a transverse field.
    let (mut points, _) = fixture(16);
    for p in &mut points {
        p[2] = 2.5;
    }
    let kernel = fs_fmm::HelmholtzGradient::new(50.0 / side(&points), 2).unwrap();
    let op = FmmC64::new_directional(&kernel, &points, options(1e-3, 0)).unwrap();
    assert!(op.directional_stats().unwrap().high_frequency_pairs > 0);
    let q = charges(points.len());
    for value in op
        .potentials(&q)
        .unwrap()
        .into_iter()
        .chain(op.direct(&q).unwrap())
    {
        assert_eq!(value, C64::ZERO);
    }
}

#[test]
fn dfmm_translation_rotation_invariance() {
    let (points, _) = fixture(16);
    let q = charges(points.len());
    let kernel = HelmholtzSingle {
        k: 130.0 / side(&points),
    };
    let eps = 1e-5;
    let mut o = options(eps, 1);
    o.low_freq_threshold = 4.0;
    let op = FmmC64::new_directional(&kernel, &points, o.clone()).unwrap();
    let reference = op.potentials(&q).unwrap();
    let (s, c) = det::sin_cos(0.37);
    for (rung, rotated) in [false, true].into_iter().enumerate() {
        let moved: Vec<_> = points
            .iter()
            .map(|p| {
                if rotated {
                    [c * p[0] - s * p[1], s * p[0] + c * p[1], p[2]]
                } else {
                    [p[0] + 8.0, p[1] - 4.0, p[2] + 2.0]
                }
            })
            .collect();
        let started = Instant::now();
        let actual = FmmC64::new_directional(&kernel, &moved, o.clone())
            .unwrap()
            .potentials(&q)
            .unwrap();
        let error = rel(&actual, &reference);
        emit(
            "dfmm-rigid-transform",
            format!(
                "{{\"case\":\"dfmm-rigid-transform\",\"inputs\":\"analytic two-cluster cloud; rigid coordinate transform\",\"axis\":\"translation-rotation\",\"rung\":{rung},\"observed\":{error},\"err\":{error},\"tol\":{},\"verdict\":\"{}\",\"elapsed_s\":{}}}",
                10.0 * eps,
                if error < 10.0 * eps { "pass" } else { "fail" },
                started.elapsed().as_secs_f64()
            ),
        );
        assert!(error < 10.0 * eps);
    }
}

#[test]
fn dfmm_caps_refuse() {
    let started = Instant::now();
    let (points, normals) = fixture(16);
    let kernel = HelmholtzSingle {
        k: 370.0 / side(&points),
    };
    let o = options(1e-7, 2);
    let op = operator(3, &kernel, &points, &normals, o.clone());
    let needed = *op
        .directional_stats()
        .unwrap()
        .cones_per_level
        .iter()
        .max()
        .unwrap();
    let mut cap = o.clone();
    cap.max_cones_per_level = needed;
    assert!(FmmC64::new_directional(&kernel, &points, cap.clone()).is_ok());
    cap.max_cones_per_level = needed - 1;
    assert!(matches!(
        FmmC64::new_directional(&kernel, &points, cap),
        Err(FmmError::WorkEnvelopeExceeded {
            operation: "directional cones per level",
            ..
        })
    ));
    let mut at = o.clone();
    at.max_cones_per_level = MAX_DIRECTIONAL_CONES;
    at.order = MAX_INTERPOLATION_ORDER;
    assert!(at.validate().is_ok());
    at.order += 1;
    assert!(matches!(at.validate(), Err(FmmError::InvalidOrder { .. })));
    at = o.clone();
    at.max_cones_per_level = MAX_DIRECTIONAL_CONES + 1;
    assert!(matches!(
        at.validate(),
        Err(FmmError::WorkEnvelopeExceeded {
            operation: "directional cone grid",
            ..
        })
    ));
    assert!(FmmC64::validate_coefficient_slots(MAX_COMPLEX_COEFFICIENT_SLOTS).is_ok());
    assert!(matches!(
        FmmC64::validate_coefficient_slots(MAX_COMPLEX_COEFFICIENT_SLOTS + 1),
        Err(FmmError::WorkEnvelopeExceeded { .. })
    ));
    let mut too_low = o.clone();
    too_low.order = 5;
    assert!(matches!(
        too_low.validate(),
        Err(FmmError::DirectionalAccuracy { .. })
    ));
    let mut rank = o.clone();
    rank.max_rank = fs_fmm::MAX_DIRECTIONAL_RANK;
    assert!(rank.validate().is_ok());
    rank.max_rank += 1;
    assert!(matches!(
        rank.validate(),
        Err(FmmError::WorkEnvelopeExceeded {
            operation: "directional coupling rank",
            ..
        })
    ));
    assert!(matches!(
        DirectionalOptions::for_tolerance(1e-12),
        Err(FmmError::DirectionalAccuracy { .. })
    ));
    let q = charges(points.len());
    assert!(matches!(
        op.direct_at_targets(&q, &[points.len()]),
        Err(FmmError::ComplexRequest { .. })
    ));
    assert!(matches!(
        op.direct_at_targets(&q, &vec![0; 1001]),
        Err(FmmError::WorkEnvelopeExceeded {
            operation: "sampled-direct target",
            ..
        })
    ));
    assert!(op.direct_at_targets(&q, &vec![0; 1000]).is_ok());
    assert!(matches!(
        op.potentials(&q[..15]),
        Err(FmmError::InvalidCharges { .. })
    ));
    let mut nonfinite = q;
    nonfinite[0] = C64::new(f64::NAN, 0.0);
    assert!(matches!(
        op.potentials(&nonfinite),
        Err(FmmError::InvalidCharges { .. })
    ));
    assert!(FmmC64::new_directional(&kernel, &[], o).is_err());
    let mut depth_options = DirectionalOptions::for_tolerance(1e-3).unwrap();
    depth_options.low_freq_threshold = 4.0;
    let at_depth = HelmholtzSingle {
        k: 4.0 * f64::from(1u32 << fs_fmm::MAX_DIRECTIONAL_DEPTH) / side(&points) * (1.0 - 1e-12),
    };
    let at = FmmC64::new_directional(&at_depth, &points, depth_options.clone()).unwrap();
    assert_eq!(
        at.directional_stats().unwrap().cones_per_level.len(),
        fs_fmm::MAX_DIRECTIONAL_DEPTH as usize + 1
    );
    let past = HelmholtzSingle {
        k: 2.0 * at_depth.k,
    };
    assert!(matches!(
        FmmC64::new_directional(&past, &points, depth_options),
        Err(FmmError::WorkEnvelopeExceeded {
            operation: "directional depth",
            ..
        })
    ));
    emit(
        "dfmm-caps",
        format!(
            "{{\"case\":\"dfmm-caps\",\"inputs\":\"analytic two-cluster cloud; cone/order/rank/depth/coefficient/oracle boundaries\",\"axis\":\"cap-plus-one\",\"rung\":1,\"observed\":{needed},\"oracle\":\"cap accepted; one past refused by name\",\"err\":0,\"tol\":0,\"verdict\":\"refused\",\"elapsed_s\":{}}}",
            started.elapsed().as_secs_f64()
        ),
    );
}

#[test]
fn dfmm_declared_envelopes_are_neutral_and_enforced() {
    // The ka = 5 unit-sphere regime of the impedance sphere series (kD about
    // 10, shared low-frequency couplings compressed), at a size a test affords.
    let points = sphere(1024);
    let normals = points.clone();
    let q = charges(points.len());
    let kernel = HelmholtzSingle { k: 5.0 };
    let mut o = DirectionalOptions::for_tolerance(1e-3).unwrap();
    o.leaf_capacity = 16;
    let build = |kind: usize, o: DirectionalOptions| match kind {
        0 => FmmC64::new_directional(&kernel, &points, o),
        _ => FmmC64::new_double_directional(
            HelmholtzDouble { k: kernel.k },
            points.clone(),
            normals.clone(),
            o,
        ),
    };
    for kind in [0, 1] {
        let started = Instant::now();
        let base = build(kind, o.clone()).unwrap();
        let stats = base.directional_stats().unwrap().clone();
        assert!(
            stats.setup_work > 0 && stats.maximum_rank > 0,
            "the fixture must compress a shared coupling: {stats:?}"
        );
        let reference = bits(&base.potentials(&q).unwrap());
        // The widest declarable budgets admit the same plan, bit for bit.
        let mut wide = o.clone();
        wide.max_coefficient_slots = MAX_DECLARED_COEFFICIENT_SLOTS;
        wide.max_setup_work = usize::MAX;
        let wide_op = build(kind, wide).unwrap();
        assert_eq!(reference, bits(&wide_op.potentials(&q).unwrap()));
        assert_eq!(wide_op.stats(), base.stats());
        // Setup is admitted at exactly its counted products, refused one below.
        let mut setup = o.clone();
        setup.max_setup_work = stats.setup_work;
        assert_eq!(
            reference,
            bits(&build(kind, setup.clone()).unwrap().potentials(&q).unwrap())
        );
        setup.max_setup_work = stats.setup_work - 1;
        assert!(matches!(
            build(kind, setup),
            Err(FmmError::WorkEnvelopeExceeded {
                operation: "directional coupling setup",
                requested,
                max,
            }) if requested == stats.setup_work && max == stats.setup_work - 1
        ));
        // A slot budget below the retained factors refuses before allocation.
        let mut slots = o.clone();
        slots.max_coefficient_slots = stats.coefficient_slots - 1;
        assert!(matches!(
            build(kind, slots),
            Err(FmmError::WorkEnvelopeExceeded {
                operation: "complex coefficient-slot (two f64 per C64)"
                    | "directional basis coefficient-slot",
                ..
            })
        ));
        emit(
            "dfmm-declared-envelopes",
            format!(
                "{{\"case\":\"dfmm-declared-envelopes\",\"inputs\":\"Fibonacci unit sphere, 1024 points, k=5\",\"axis\":\"kernel-kind\",\"rung\":{kind},\"observed\":{{\"setup_work\":{},\"coefficient_slots\":{},\"maximum_rank\":{},\"coupling_blocks\":{}}},\"verdict\":\"pass\",\"elapsed_s\":{}}}",
                stats.setup_work,
                stats.coefficient_slots,
                stats.maximum_rank,
                stats.coupling_blocks,
                started.elapsed().as_secs_f64()
            ),
        );
    }
    let mut bad = o.clone();
    bad.max_coefficient_slots = 0;
    assert!(matches!(
        bad.validate(),
        Err(FmmError::ComplexRequest { .. })
    ));
    bad = o.clone();
    bad.max_setup_work = 0;
    assert!(matches!(
        bad.validate(),
        Err(FmmError::ComplexRequest { .. })
    ));
    bad = o.clone();
    bad.max_coefficient_slots = MAX_DECLARED_COEFFICIENT_SLOTS;
    assert!(bad.validate().is_ok());
    bad.max_coefficient_slots += 1;
    assert!(matches!(
        bad.validate(),
        Err(FmmError::WorkEnvelopeExceeded {
            operation: "directional declared coefficient-slot",
            ..
        })
    ));
    // k = 0 is below every low-frequency threshold, so the request delegates
    // whole to W5a, under the declared w5a_budget.
    let statics = |budget| {
        let mut o = o.clone();
        o.w5a_budget = budget;
        FmmC64::new_double_directional(
            HelmholtzDouble { k: 0.0 },
            points.clone(),
            normals.clone(),
            o,
        )
    };
    let delegated = statics(ComplexBudget::V1).unwrap();
    assert!(delegated.directional_stats().is_none());
    assert!(matches!(
        statics(ComplexBudget {
            max_translation_work: 1
        }),
        Err(FmmError::WorkEnvelopeExceeded {
            operation: "complex M2L channel multiply-add",
            max: 1,
            ..
        })
    ));
}

#[test]
fn dfmm_declared_near_pair_envelope_is_neutral_and_enforced() {
    // The near field is admitted at exactly its counted point pairs and
    // refused one below, during the pair traversal and so before any coupling
    // setup. Any admitting budget leaves every potential bit unchanged.
    let started = Instant::now();
    let points = sphere(1024);
    let q = charges(points.len());
    let kernel = HelmholtzSingle { k: 5.0 };
    let mut o = DirectionalOptions::for_tolerance(1e-3).unwrap();
    o.leaf_capacity = 16;
    let build = |o: DirectionalOptions| FmmC64::new_directional(&kernel, &points, o);
    let base = build(o.clone()).unwrap();
    let pairs = base.directional_stats().unwrap().near_pairs;
    assert!(pairs > 0, "the fixture must have a direct near field");
    let reference = bits(&base.potentials(&q).unwrap());
    for declared in [pairs, usize::MAX] {
        let mut admitted = o.clone();
        admitted.max_near_pairs = declared;
        assert_eq!(
            reference,
            bits(&build(admitted).unwrap().potentials(&q).unwrap())
        );
    }
    let mut refused = o.clone();
    refused.max_near_pairs = pairs - 1;
    assert!(matches!(
        build(refused),
        Err(FmmError::WorkEnvelopeExceeded {
            operation: "directional near-field point pair",
            requested,
            max,
        }) if requested == pairs && max == pairs - 1
    ));
    let mut zero = o.clone();
    zero.max_near_pairs = 0;
    assert!(matches!(
        zero.validate(),
        Err(FmmError::ComplexRequest { .. })
    ));
    emit(
        "dfmm-declared-near-pairs",
        format!(
            "{{\"case\":\"dfmm-declared-near-pairs\",\"inputs\":\"Fibonacci unit sphere, 1024 points, k=5\",\"axis\":\"near-pair budget\",\"rung\":0,\"observed\":{{\"near_pairs\":{pairs}}},\"verdict\":\"pass\",\"elapsed_s\":{}}}",
            started.elapsed().as_secs_f64()
        ),
    );
}

#[test]
fn dfmm_bitwise_repeat() {
    let (points, normals) = fixture(16);
    let q = charges(16);
    let kernel = HelmholtzSingle {
        k: 130.0 / side(&points),
    };
    for kind in 0..4 {
        let started = Instant::now();
        let o = options(1e-5, 1);
        let a = operator(kind, &kernel, &points, &normals, o.clone());
        let b = operator(kind, &kernel, &points, &normals, o);
        let reference = bits(&a.potentials(&q).unwrap());
        assert_eq!(reference, bits(&a.potentials(&q).unwrap()));
        assert_eq!(reference, bits(&b.potentials(&q).unwrap()));
        assert_eq!(a.stats(), b.stats());
        emit(
            "dfmm-bitwise",
            format!(
                "{{\"case\":\"dfmm-bitwise\",\"inputs\":\"analytic two-cluster cloud, repeated apply and independent reconstruction\",\"seed\":{SEED},\"axis\":\"kernel-kind\",\"rung\":{kind},\"observed\":0,\"err\":0,\"tol\":0,\"verdict\":\"pass\",\"elapsed_s\":{}}}",
                started.elapsed().as_secs_f64()
            ),
        );
    }
}
