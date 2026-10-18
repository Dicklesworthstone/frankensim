//! G0/G1/G2/G3/G5 realization checks. Analytic sphere realization error is
//! separated from BEM spatial error. CQ comparisons use the SAME actual BEM
//! panels/receiver rows on a sphere and an authored lathed finite-wall bell.
use fs_bem::{
    helmholtz::{self, Medium, laplace},
    observation_transfers::*,
    panel3d::SpherePanels,
};
use fs_blake3::hash_domain;
use fs_math::{c64::C64, det};
use fs_mms::{LadderSide, OrderGate, RefinementLadder};
use std::time::Instant;

const RHO: f64 = 1.2;
const C: f64 = 343.0;
const A: f64 = 0.05;
const R: f64 = 3.48;
const FS: f64 = 8000.0;
const PULSE_S: f64 = 0.04;
const ANALYTIC_INPUTS: &str = "analytic:sphere:a=.05:r=3.48:rho=1.2:c=343:pass=1000:stop=3000:latency=64;sin^8 pulse T=.04 SI";
fn medium() -> Medium {
    Medium {
        density: RHO,
        sound_speed: C,
    }
}
fn phase(x: f64) -> C64 {
    C64::new(det::cos(x), det::sin(x))
}
fn declaration() -> TransferDeclaration {
    TransferDeclaration {
        producer_hash: hash_domain(
            "analytic:pulsating-sphere:SI",
            b"a=.05,r=3.48,rho=1.2,c=343",
        ),
        sources: vec![SourceId::StructuralMode(
            "unit outward radial boundary velocity [m/s]".into(),
        )],
        receivers: vec![ReceiverSpec {
            id: "pressure+Cartesian-gradient:SI".into(),
            position_m: [0.0, 0.0, R],
        }],
        fs_hz: FS,
        geometric_delays_s: vec![(R - A) / C],
        pass_band_upper_hz: 1000.0,
        taper_stop_hz: 3000.0,
        filter_latency_samples: 64,
        fft_lengths: vec![128, 256, 512, 1024],
        tail_energy_tolerance: 1e-10,
        dtft_relative_tolerance: 1e-4,
        producer_uncertainty: ProducerUncertainty::ExactAnalytic,
        maximum_frequency_entries: MAX_OBSERVATION_FFT * 4,
    }
}
fn analytic_sphere(s: C64, a: f64, r: f64, delay_s: f64) -> Vec<C64> {
    let st = s.scale(a / C);
    let h = st
        * (C64::ONE + st).recip()
        * phase(-s.im * delay_s).scale(det::exp(-s.re * delay_s) * RHO * C * a / r);
    vec![
        h,
        C64::ZERO,
        C64::ZERO,
        h * (s.scale(-1.0 / C) - C64::from_re(1.0 / r)),
    ]
}
fn analytic_bake(d: &TransferDeclaration) -> ObservationTransferBake {
    bake_observation_transfers(
        |w| Ok(analytic_sphere(C64::new(0.0, -w), A, R, (R - A) / C)),
        d,
        || false,
    )
    .unwrap()
}
fn pulse(t: f64) -> f64 {
    if t <= 0.0 || t >= PULSE_S {
        0.0
    } else {
        det::sin(core::f64::consts::PI * t / PULSE_S).powi(8)
    }
}
// sin^8(theta)=(35-56 cos(2theta)+28 cos(4theta)-8 cos(6theta)+cos(8theta))/128.
// Integrating each harmonic against exp(-t/tau)/tau gives the exact near-field
// memory. This is independent of DFT, vector fit, BEM and the CQ recurrence.
fn sphere_exact(t: f64, a: f64, r: f64) -> [f64; 2] {
    if t <= 0.0 {
        return [0.0, 0.0];
    }
    let tau = a / C;
    let x = t.min(PULSE_S);
    let mut memory = 0.0;
    let mut u = 0.0;
    let mut du = 0.0;
    for (j, coefficient) in [35.0, -56.0, 28.0, -8.0, 1.0].into_iter().enumerate() {
        let w = core::f64::consts::TAU * j as f64 / PULSE_S;
        let b = w * tau;
        let coefficient = coefficient / 128.0;
        let angle = w * x;
        memory += coefficient * (det::cos(angle) + b * det::sin(angle) - det::exp(-x / tau))
            / (1.0 + b * b);
        u += coefficient * det::cos(angle);
        du -= coefficient * w * det::sin(angle);
    }
    if t >= PULSE_S {
        u = 0.0;
        du = 0.0;
        memory *= det::exp(-(t - PULSE_S) / tau);
    }
    let scale = RHO * C * a / r;
    let p = scale * (u - memory);
    let derivative = scale * (du - (u - memory) / tau);
    [p, -p / r - derivative / C]
}
fn relative(a: &[f64], b: &[f64]) -> f64 {
    assert_eq!(a.len(), b.len());
    let norm = b.iter().map(|v| v * v).sum::<f64>();
    assert!(norm > 0.0);
    det::sqrt(a.iter().zip(b).map(|(a, b)| (a - b) * (a - b)).sum::<f64>() / norm)
}
fn record(
    case: &str,
    inputs: &str,
    axis: &str,
    rung: usize,
    observed: f64,
    tol: f64,
    elapsed: f64,
) {
    println!(
        "{{\"suite\":\"fs-bem\",\"case\":\"{case}\",\"inputs\":\"{inputs}\",\"axis\":\"{axis}\",\"rung\":{rung},\"observed\":{observed},\"err\":{observed},\"tol\":{tol},\"verdict\":\"{}\",\"elapsed_s\":{elapsed}}}",
        if observed <= tol { "pass" } else { "fail" }
    );
}

#[test]
fn fir_transfer_causal() {
    let start = Instant::now();
    let d = declaration();
    let bake = analytic_bake(&d);
    let mut observer = bake.observer().unwrap();
    let samples = 2048;
    let mut pre = 0.0;
    let mut total = 0.0;
    let first = (d.geometric_delays_s[0] * FS - d.filter_latency_samples as f64)
        .ceil()
        .max(0.0) as usize;
    for n in 0..samples {
        let y = observer.step(&[if n == 0 { 1.0 } else { 0.0 }]).unwrap()[0];
        total += y * y;
        if n < first {
            pre += y * y;
        }
    }
    assert!(total > 0.0);
    let fraction = pre / total;
    record(
        "obs-fir-causal",
        ANALYTIC_INPUTS,
        "pre-arrival-energy",
        1024,
        fraction,
        1e-10,
        start.elapsed().as_secs_f64(),
    );
    assert!(fraction < 1e-10);
    assert!(bake.ladder.last().unwrap().wrap_energy_fraction < 1e-10);
    // Wrong Fourier sign must fail the causal/wrap gates; it cannot be fixed
    // by making a real kernel, delaying it circularly, or clipping its tail.
    let wrong = bake_observation_transfers(
        |w| Ok(analytic_sphere(C64::new(0.0, w), A, R, (R - A) / C)),
        &d,
        || false,
    );
    assert!(
        matches!(wrong, Err(TransferError::NotConverged { .. })),
        "wrong sign admitted"
    );
}

#[test]
fn fir_pulsating_sphere_exact() {
    let start = Instant::now();
    let d = declaration();
    let bake = analytic_bake(&d);
    let mut observer = bake.observer().unwrap();
    let mut actual = [Vec::new(), Vec::new()];
    let mut expected = [Vec::new(), Vec::new()];
    for n in 0..1024 {
        let t = n as f64 / FS;
        let y = observer.step(&[pulse(t)]).unwrap();
        let oracle = sphere_exact(
            t - d.geometric_delays_s[0] - d.filter_latency_samples as f64 / FS,
            A,
            R,
        );
        actual[0].push(y[0]);
        actual[1].push(y[3]);
        expected[0].push(oracle[0]);
        expected[1].push(oracle[1]);
    }
    for c in 0..2 {
        let error = relative(&actual[c], &expected[c]);
        record(
            "obs-fir-sphere-exact",
            ANALYTIC_INPUTS,
            if c == 0 { "p" } else { "gradient" },
            1024,
            error,
            1e-6,
            start.elapsed().as_secs_f64(),
        );
        assert!(error <= 1e-6, "component {c}: {error:e}");
    }
}

#[test]
fn fir_length_ladder_converges() {
    let start = Instant::now();
    let mut d = declaration();
    d.fft_lengths = vec![128, 256, 512, 1024];
    let bake = analytic_bake(&d);
    assert!(
        bake.ladder[0].tail_energy_fraction > d.tail_energy_tolerance,
        "FFT-length axis is dead"
    );
    for row in &bake.ladder {
        println!(
            "{{\"suite\":\"fs-bem\",\"case\":\"obs-fir-length\",\"inputs\":\"analytic sphere a=.05,r=3.48,rho=1.2,c=343,Fs=8000,pass=1000,stop=3000,latency=64\",\"axis\":\"n_fft\",\"rung\":{},\"n_fft\":{},\"tail_energy\":{},\"wrap_energy\":{},\"dtft_rel_change\":{},\"observed\":{},\"tol\":1e-10,\"verdict\":\"{}\",\"elapsed_s\":{}}}",
            row.n_fft,
            row.n_fft,
            row.tail_energy_fraction,
            row.wrap_energy_fraction,
            row.dtft_relative_change
                .map_or("null".to_string(), |v| v.to_string()),
            row.tail_energy_fraction,
            if row.tail_energy_fraction <= 1e-10 {
                "pass"
            } else {
                "fail"
            },
            start.elapsed().as_secs_f64()
        );
    }
    let last = bake.ladder.last().unwrap();
    assert!(last.tail_energy_fraction <= d.tail_energy_tolerance);
    assert!(last.wrap_energy_fraction < 1e-10);
    assert!(last.dtft_relative_change.unwrap() < 1e-4);
    assert!(last.dtft_relative_change.unwrap() < bake.ladder[1].dtft_relative_change.unwrap());
}

#[test]
fn cq_order_ladder() {
    let start = Instant::now();
    let mut hs = Vec::new();
    let mut errors = Vec::new();
    for rate in [2000usize, 4000, 8000, 16000] {
        let steps = rate * 64 / 1000;
        let dt = 1.0 / rate as f64;
        let input: Vec<_> = (0..steps).map(|n| pulse(n as f64 * dt)).collect();
        let nfft = (steps * 2).next_power_of_two();
        let output = cq_time_response(
            |s| Ok(vec![analytic_sphere(s, A, R, 0.0)[0]]),
            &[input],
            1,
            dt,
            nfft,
            20.0,
            || false,
        )
        .unwrap();
        let exact: Vec<_> = (0..steps)
            .map(|n| sphere_exact(n as f64 * dt, A, R)[0])
            .collect();
        let error = relative(&output[0], &exact);
        hs.push(dt);
        errors.push(error);
        record(
            "obs-cq-order",
            "analytic:delay-extracted sphere a=.05,r=3.48,rho=1.2,c=343;sin^8 pulse T=.04 SI;Lubich 1988 BDF2",
            "dt",
            rate,
            error,
            1e-3,
            start.elapsed().as_secs_f64(),
        );
    }
    assert!(errors.windows(2).all(|e| e[1] < e[0]));
    assert!(errors[3] < 1e-3);
    let verdict = OrderGate { theoretical: 2.0 }
        .check(
            "obs-cq-bdf2",
            LadderSide::Primal,
            &RefinementLadder::new(hs, errors).unwrap(),
        )
        .unwrap();
    println!("{}", verdict.json_line(true));
}

// Authored lathed solid bell, finite wall: throat cap, outer flare, annular lip,
// inner flare and mouth-facing throat disc. Geometry alone; no fitted acoustics.
fn bell() -> SpherePanels {
    let around = 8usize;
    let ring = |z: f64, r: f64| {
        (0..around)
            .map(|j| {
                let angle = core::f64::consts::TAU * j as f64 / around as f64;
                [r * det::cos(angle), r * det::sin(angle), z]
            })
            .collect::<Vec<_>>()
    };
    let rings = [
        ring(-0.00004, 0.00006),
        ring(0.00016, 0.00012),
        ring(0.00016, 0.00010),
        ring(0.0, 0.00004),
    ];
    let mut triangles = Vec::new();
    for pair in rings.windows(2) {
        for j in 0..around {
            let k = (j + 1) % around;
            triangles.extend([
                [pair[0][j], pair[0][k], pair[1][j]],
                [pair[0][k], pair[1][k], pair[1][j]],
            ]);
        }
    }
    for j in 0..around {
        let k = (j + 1) % around;
        triangles.push([[0.0, 0.0, -0.00004], rings[0][k], rings[0][j]]);
        triangles.push([[0.0, 0.0, 0.0], rings[3][j], rings[3][k]]);
    }
    SpherePanels::from_triangles(triangles).unwrap()
}

#[test]
fn fir_matches_convolution_quadrature() {
    let start = Instant::now();
    for (name, surface) in [
        ("sphere", SpherePanels::icosphere(0.0001, 1).unwrap()),
        ("lathed-bell", bell()),
    ] {
        let sources = vec![
            SourceId::StructuralMode(format!("{name}:authored unit outward velocity")),
            SourceId::ApertureVolumeFlow(format!("{name}:unit SI volume flow")),
        ];
        // A distinct driven cap, normalized by its physical SI area, exercises
        // aperture flow independently of the body's uniform mode coordinate.
        let cap_area = surface
            .areas()
            .iter()
            .zip(surface.normals())
            .filter(|(_, n)| n[2] > 0.5)
            .map(|(a, _)| a)
            .sum::<f64>();
        assert!(cap_area > 0.0);
        let aperture: Vec<_> = surface
            .normals()
            .iter()
            .map(|n| if n[2] > 0.5 { 1.0 / cap_area } else { 0.0 })
            .collect();
        let traces = vec![vec![1.0; surface.areas().len()], aperture];
        let receivers = vec![ReceiverSpec {
            id: "on-axis p+grad".into(),
            position_m: [0.0, 0.0, 0.01],
        }];
        let producer =
            BemObservation::new(&surface, medium(), &sources, &traces, &receivers, 0.001).unwrap();
        let mut d = declaration();
        d.sources = sources;
        d.receivers = receivers;
        d.producer_hash = producer.operator_hash();
        d.geometric_delays_s = producer.geometric_delays_s();
        d.producer_uncertainty = ProducerUncertainty::Unresolved {
            source_hash: producer.operator_hash(),
        };
        let bake = producer.bake(&d, || false).unwrap();
        let inputs = format!(
            "authored:{name}:panels={}:rho=1.2:c=343:receiver=(0,0,.01):unit mode+area-normalized cap;declaration-blake3={};producer-error=unresolved",
            surface.areas().len(),
            bake.declaration_hash.to_hex()
        );
        let mut observer = bake.observer().unwrap();
        let steps = 1024;
        let mut fir = [Vec::new(), Vec::new()];
        for n in 0..steps {
            let t = n as f64 / FS;
            let y = observer
                .step(&[pulse(t), cap_area * pulse(t - 0.003)])
                .unwrap();
            fir[0].push(y[0]);
            fir[1].push(y[3]);
        }
        let mut hs = Vec::new();
        let mut errors = Vec::new();
        for multiple in [1usize, 2, 4, 8] {
            let dt = 1.0 / (FS * multiple as f64);
            let count = steps * multiple;
            let input: Vec<_> = (0..count).map(|n| pulse(n as f64 * dt)).collect();
            let aperture_input: Vec<_> = (0..count)
                .map(|n| cap_area * pulse(n as f64 * dt - 0.003))
                .collect();
            // Both source coordinates drive both p and its radial gradient.
            // The BEM producer retains all four receiver components; no source
            // column is discarded at this reference comparison seam.
            let output = cq_time_response(
                |s| {
                    let row = producer.laplace_sample(s, &[0.0])?;
                    Ok(vec![row[0], row[3], row[4], row[7]])
                },
                &[input, aperture_input],
                2,
                dt,
                (count * 2).next_power_of_two(),
                20.0,
                || false,
            )
            .unwrap();
            let latency = d.filter_latency_samples;
            let mut error = 0.0_f64;
            for c in 0..2 {
                let aligned: Vec<_> = (0..steps)
                    .map(|n| {
                        if n < latency {
                            0.0
                        } else {
                            output[c][(n - latency) * multiple]
                        }
                    })
                    .collect();
                let component_error = relative(&aligned, &fir[c]);
                error = error.max(component_error);
                record(
                    &format!(
                        "obs-fir-cq-{name}-{}",
                        if c == 0 { "p" } else { "gradient" }
                    ),
                    &inputs,
                    "cq-dt",
                    multiple,
                    component_error,
                    1e-3,
                    start.elapsed().as_secs_f64(),
                );
            }
            hs.push(dt);
            errors.push(error);
            record(
                &format!("obs-fir-cq-{name}"),
                &inputs,
                "cq-dt",
                multiple,
                error,
                1e-3,
                start.elapsed().as_secs_f64(),
            );
        }
        assert!(errors.windows(2).all(|e| e[1] < e[0]), "{name}: {errors:?}");
        assert!(errors[3] < 1e-3, "{name}: {errors:?}");
        let verdict = OrderGate { theoretical: 2.0 }
            .check(
                &format!("obs-cq-bem-{name}"),
                LadderSide::Primal,
                &RefinementLadder::new(hs, errors).unwrap(),
            )
            .unwrap();
        println!("{}", verdict.json_line(true));
    }
}

#[test]
fn fir_transfer_refusals() {
    assert!(admit_observation_fft(MAX_OBSERVATION_FFT).is_ok());
    assert!(matches!(
        admit_observation_fft(MAX_OBSERVATION_FFT + 1),
        Err(TransferError::Cap { .. })
    ));
    let mut cap_decl = declaration();
    cap_decl.fft_lengths = vec![32, 64, 128, 256];
    cap_decl.maximum_frequency_entries = 4_194_304;
    cap_decl.sources = (0..4096)
        .map(|i| SourceId::StructuralMode(format!("SI-mode-{i}")))
        .collect();
    cap_decl.geometric_delays_s = vec![0.01; 4096];
    assert_eq!(cap_decl.admit_work().unwrap(), 4096 * 4);
    cap_decl
        .sources
        .push(SourceId::StructuralMode("cap+1".into()));
    assert!(matches!(
        cap_decl.admit_work(),
        Err(TransferError::Cap { .. })
    ));
    cap_decl = declaration();
    cap_decl.receivers = (0..64)
        .map(|i| ReceiverSpec {
            id: format!("SI-receiver-{i}"),
            position_m: [0.0, 0.0, R],
        })
        .collect();
    cap_decl.geometric_delays_s = vec![0.01; 64];
    assert_eq!(cap_decl.admit_work().unwrap(), 256);
    cap_decl.receivers.push(ReceiverSpec {
        id: "cap+1".into(),
        position_m: [0.0, 0.0, R],
    });
    assert!(matches!(
        cap_decl.admit_work(),
        Err(TransferError::Cap { .. })
    ));
    cap_decl = declaration();
    cap_decl.fft_lengths = (0..8).map(|i| 32 << i).collect();
    assert!(cap_decl.admit_work().is_ok());
    cap_decl.fft_lengths.push(8192);
    assert!(matches!(
        cap_decl.admit_work(),
        Err(TransferError::Cap { .. })
    ));
    let mut d = declaration();
    d.fft_lengths = vec![32, 64, 128, 256];
    d.filter_latency_samples = 100;
    assert!(matches!(
        bake_observation_transfers(
            |w| Ok(analytic_sphere(C64::new(0.0, -w), A, R, (R - A) / C)),
            &d,
            || false
        ),
        Err(TransferError::NotConverged { .. })
    ));
    d = declaration();
    d.fft_lengths[3] = MAX_OBSERVATION_FFT + 1;
    assert!(
        bake_observation_transfers(|_| panic!("inadmissible FFT called producer"), &d, || false)
            .is_err()
    );
    assert!(matches!(
        bake_observation_transfers(|_| panic!("cancelled producer"), &declaration(), || true),
        Err(TransferError::Cancelled)
    ));
    assert!(bake_observation_transfers(|_| Ok(vec![C64::ZERO]), &declaration(), || false).is_err());
    let surface = SpherePanels::icosphere(0.0001, 0).unwrap();
    let source = [SourceId::StructuralMode("unit SI".into())];
    let trace = vec![vec![1.0; surface.areas().len()]];
    for point in [[0.0; 3], [0.0, 0.0, 0.0001]] {
        assert!(
            BemObservation::new(
                &surface,
                medium(),
                &source,
                &trace,
                &[ReceiverSpec {
                    id: "refused".into(),
                    position_m: point
                }],
                0.0001
            )
            .is_err()
        );
    }
    let receiver = [ReceiverSpec {
        id: "explicit".into(),
        position_m: [0.0, 0.0, 0.01],
    }];
    assert!(
        BemObservation::new(
            &surface,
            medium(),
            &[SourceId::ApertureVolumeFlow("not normalized".into())],
            &trace,
            &receiver,
            0.001
        )
        .is_err()
    );
    let owner = BemObservation::new(&surface, medium(), &source, &trace, &receiver, 0.001).unwrap();
    assert!(owner.bake(&declaration(), || false).is_err());
    assert!(
        cq_time_response(
            |_| panic!("bad CQ called producer"),
            &[vec![0.0; 33]],
            1,
            1e-4,
            64,
            20.0,
            || false
        )
        .is_err()
    );
    let zero = vec![C64::ZERO; surface.areas().len()];
    let at_cap = vec![zero.as_slice(); 256];
    assert_eq!(
        laplace::solve_radiation_laplace_batch(&surface, C64::new(100.0, 200.0), medium(), &at_cap)
            .unwrap()
            .len(),
        256
    );
    let over_cap = vec![zero.as_slice(); 257];
    assert!(
        laplace::solve_radiation_laplace_batch(
            &surface,
            C64::new(100.0, 200.0),
            medium(),
            &over_cap
        )
        .is_err()
    );
    let solved = laplace::solve_radiation_laplace_batch(
        &surface,
        C64::new(100.0, 200.0),
        medium(),
        &[zero.as_slice()],
    )
    .unwrap();
    assert!(
        laplace::observe_laplace(&surface, &solved[0], &[[0.0, 0.0, 0.01]], 0.001, &[0.01])
            .is_err()
    );
}

#[test]
fn fir_transfer_bitwise_repeat() {
    let d = declaration();
    let a = analytic_bake(&d);
    let b = analytic_bake(&d);
    assert_eq!(a.declaration_hash, b.declaration_hash);
    assert_eq!(a.ladder, b.ladder);
    for (a, b) in a.kernels.iter().zip(&b.kernels) {
        assert_eq!(a.integer_delay_samples, b.integer_delay_samples);
        for (a, b) in a.taps.iter().zip(&b.taps) {
            assert_eq!(a.to_bits(), b.to_bits());
        }
    }
    let mut different = d.clone();
    different.filter_latency_samples += 1;
    assert_ne!(d.declaration_hash(), different.declaration_hash());
}

#[test]
fn fir_ladder_is_scale_invariant_and_rejects_unresolved_tails() {
    // G3: squared physical units may overflow/underflow while the response is
    // finite. Neither event may turn an unresolved tail into a passing gate.
    let d = declaration();
    let reference = analytic_bake(&d);
    for scale in [1e-180, 1e180] {
        let sample = |w: f64| {
            Ok(analytic_sphere(C64::new(0.0, -w), A, R, (R - A) / C)
                .into_iter()
                .map(|z| z.scale(scale))
                .collect())
        };
        let measured = bake_observation_transfers(sample, &d, || false).unwrap();
        let actual = measured.ladder.last().unwrap();
        let expected = reference.ladder.last().unwrap();
        assert!(actual.tail_energy_fraction.is_finite());
        assert!(actual.wrap_energy_fraction.is_finite());
        assert!(
            (actual.tail_energy_fraction - expected.tail_energy_fraction).abs()
                < 1e-3 * expected.tail_energy_fraction + 1e-20
        );
        let mut coarse = d.clone();
        coarse.fft_lengths = vec![32, 64, 128, 256];
        coarse.filter_latency_samples = 100;
        assert!(matches!(
            bake_observation_transfers(sample, &coarse, || false),
            Err(TransferError::NotConverged { .. })
        ));
    }
}

#[test]
fn zero_by_symmetry_components_are_judged_with_their_gradient_vector() {
    // An on-axis receiver's transverse gradient vanishes by symmetry; a real BEM
    // producer leaves roundoff there (~1e-17 of |grad p|), whose self-relative
    // tail is ~3/4 at every DFT size. Against the whole gradient vector that
    // roundoff passes, while 1e-3 of an unresolved Q=1000 resonance in the same
    // component is still refused.
    let d = declaration();
    let reference = analytic_bake(&d);
    let expected = reference.ladder.last().unwrap();
    let transverse = |extra: fn(f64, C64) -> C64| {
        move |w: f64| {
            let mut row = analytic_sphere(C64::new(0.0, -w), A, R, (R - A) / C);
            row[1] = extra(w, row[3]);
            Ok::<_, helmholtz::HelmholtzError>(row)
        }
    };
    let roundoff: fn(f64, C64) -> C64 = |w, g| phase(1.0e4 * w).scale(1e-17 * g.abs());
    let baked = bake_observation_transfers(transverse(roundoff), &d, || false).unwrap();
    let actual = baked.ladder.last().unwrap();
    assert!(
        (actual.tail_energy_fraction - expected.tail_energy_fraction).abs()
            < 1e-3 * expected.tail_energy_fraction + 1e-20
    );
    let resonance: fn(f64, C64) -> C64 = |w, g| {
        let w0 = core::f64::consts::TAU * 500.0;
        g.scale(1e-3 * w0 * w0) * C64::new(w0 * w0 - w * w, -w * w0 / 1000.0).recip()
    };
    assert!(matches!(
        bake_observation_transfers(transverse(resonance), &d, || false),
        Err(TransferError::NotConverged { .. })
    ));
}

#[test]
fn damped_boundary_and_receiver_match_the_original_real_frequency_owner() {
    let surface = SpherePanels::icosphere(0.001, 1).unwrap();
    let velocity = vec![C64::ONE; surface.areas().len()];
    let w = 4000.0;
    let point = [[0.0, 0.0, 0.1]];
    let real = helmholtz::solve_radiation(
        &surface,
        w / C,
        medium(),
        &velocity,
        helmholtz::Formulation::BurtonMiller,
    )
    .unwrap();
    let damped =
        laplace::solve_radiation_laplace_batch(&surface, C64::new(0.0, -w), medium(), &[&velocity])
            .unwrap();
    let reference =
        helmholtz::exterior_pressure_at_points(&surface, &real, medium(), &point).unwrap()[0];
    let observed =
        laplace::observe_laplace(&surface, &damped[0], &point, 0.001, &[0.0]).unwrap()[0];
    assert!((observed.pressure_pa - reference).abs() / reference.abs() < 1e-10);
    let eps = 1e-5;
    let plus = [[0.0, 0.0, 0.1 + eps]];
    let minus = [[0.0, 0.0, 0.1 - eps]];
    let derivative = (helmholtz::exterior_pressure_at_points(&surface, &real, medium(), &plus)
        .unwrap()[0]
        - helmholtz::exterior_pressure_at_points(&surface, &real, medium(), &minus).unwrap()[0])
        .scale(0.5 / eps);
    assert!((observed.gradient_pa_per_m[2] - derivative).abs() / derivative.abs() < 1e-6);
    assert!(
        laplace::solve_radiation_laplace_batch(
            &surface,
            C64::new(-1.0, 2.0),
            medium(),
            &[&velocity]
        )
        .is_err()
    );
}
