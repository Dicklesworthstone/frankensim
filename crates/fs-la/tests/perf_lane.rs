//! The fs-la GEMM PERF LANE (bead xdgf): packed BLIS-style gemm_f64
//! throughput against the MEASURED machine peak (fs-roofline axes) at
//! large square sizes, with the ≥75%-of-peak attainment gate. Run
//! explicitly in release:
//! `cargo test -p fs-la --release --test perf_lane -- --ignored --nocapture`
//!
//! The microkernel is the fs-simd capsule (NEON on aarch64), bitwise-
//! identical to the scalar twin — the golden 0x1d7a_a3c6_b631_7ef0 is
//! tier-invariant, verified by the gemm test suite, not here.

use fs_la::{gemm_f64, gemm_f64_parallel};
use fs_roofline::regress::{Cusum, GateSpec, GateVerdict, Night, gate, standardize};
use fs_roofline::{KernelSpec, MachineAxes, TargetAxis, Threading, attainment_for};

/// Best-of-3 measured GFLOP/s (2·m·n·k flops per GEMM).
fn measure(n: usize, reps: usize) -> f64 {
    let a: Vec<f64> = (0..n * n).map(|i| ((i as f64) * 0.13).sin()).collect();
    let b: Vec<f64> = (0..n * n).map(|i| ((i as f64) * 0.31).cos()).collect();
    let mut c = vec![0.0f64; n * n];
    gemm_f64(n, n, n, 1.0, &a, &b, 0.0, &mut c); // warm
    let mut best = f64::INFINITY;
    for _ in 0..3 {
        let t0 = std::time::Instant::now();
        for _ in 0..reps {
            gemm_f64(n, n, n, 1.0, &a, &b, 0.0, &mut c);
        }
        best = best.min(t0.elapsed().as_secs_f64() / reps as f64);
    }
    2.0 * (n * n * n) as f64 / best / 1e9
}

#[test]
#[ignore = "perf lane: run explicitly in release with --ignored"]
fn gemm_attainment() {
    let axes = MachineAxes::probe();
    println!(
        "{{\"metric\":\"axes\",\"cpu\":\"{}\",\"peak_single_gflops\":{:.1}}}",
        axes.cpu_brand, axes.peak_single_gflops
    );
    // Size ladder, ledgered; the gate reads the large square sizes.
    let mut large_best = 0.0f64;
    for &n in &[128usize, 256, 512, 1024] {
        let reps = (256 / (n / 128)).max(1) / 8 + 1;
        let gflops = measure(n, reps);
        let att = gflops / axes.peak_single_gflops;
        println!(
            "{{\"metric\":\"gemm-f64\",\"n\":{n},\"gflops\":{gflops:.2},\
             \"attainment_single\":{att:.3}}}"
        );
        if n >= 512 {
            large_best = large_best.max(att);
        }
    }
    // THE GATE: >= 75% of measured single-thread peak at large square
    // sizes on THIS machine (blocking constants MR=8 NR=4 KC=256
    // MC=128 NC=512). The second-ISA (x86-64/AVX) row is ARMED
    // PENDING hardware, per the recorded fleet census.
    println!(
        "{{\"metric\":\"gemm-gate\",\"attainment\":{large_best:.3},\"floor\":0.75,\
         \"machine\":\"{}-{}\"}}",
        std::env::consts::OS,
        std::env::consts::ARCH
    );
    assert!(
        large_best >= 0.75,
        "large-square gemm_f64 clears 75% of measured peak: {large_best:.3}"
    );
}

/// Best-of-3 all-core GFLOP/s via row-band parallel GEMM.
fn measure_parallel(n: usize, reps: usize, threads: usize) -> f64 {
    let a: Vec<f64> = (0..n * n).map(|i| ((i as f64) * 0.13).sin()).collect();
    let b: Vec<f64> = (0..n * n).map(|i| ((i as f64) * 0.31).cos()).collect();
    let mut c = vec![0.0f64; n * n];
    gemm_f64_parallel(n, n, n, 1.0, &a, &b, 0.0, &mut c, threads); // warm
    let mut best = f64::INFINITY;
    for _ in 0..3 {
        let t0 = std::time::Instant::now();
        for _ in 0..reps {
            gemm_f64_parallel(n, n, n, 1.0, &a, &b, 0.0, &mut c, threads);
        }
        best = best.min(t0.elapsed().as_secs_f64() / reps as f64);
    }
    2.0 * (n * n * n) as f64 / best / 1e9
}

/// The ALL-CORE attainment row (bead xlvx item 3): row-band parallel
/// GEMM against the measured all-core FMA axis. REPORT row by default;
/// FS_LA_ROOFLINE_GATE=1 asserts >= 0.5 (parallel GEMM leaves more on
/// the table than single-thread — memory bandwidth and band tails —
/// so the all-core floor is honest, not aspirational).
#[test]
#[ignore = "perf lane: run explicitly in release with --ignored"]
fn gemm_attainment_all_core() {
    let threads = std::env::var("FS_LA_THREADS")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .unwrap_or_else(|| std::thread::available_parallelism().map_or(8, std::num::NonZero::get));
    let axes = MachineAxes::probe();
    println!(
        "{{\"metric\":\"axes-all-core\",\"cpu\":\"{}\",\"threads\":{threads},\"peak_all_core_gflops\":{:.1}}}",
        axes.cpu_brand, axes.peak_all_core_gflops
    );
    for n in [512usize, 1024, 2048] {
        let reps = if n >= 2048 { 1 } else { 3 };
        let g = measure_parallel(n, reps, threads);
        // MIN-ROOF attainment (fs-roofline's two-axis model): per C
        // element, 2k flops; traffic model at BLIS blocking = A re-read
        // per NC column chunk + B once + C read+write per KC chunk:
        // bytes/elem = 8·(k·ceil(n/NC)/m_norm + k/m + 2·ceil(k/KC)).
        // On a bandwidth-starved box the MEMORY roof binds and the
        // compute axis is the wrong denominator (measured on ts1:
        // 219 GFLOP/s read 0.14 vs compute but the memory roof binds).
        let nf = n as f64;
        // NC_PAR_CAP=2048 (measured s5 defaults), KC=256 (bit contract).
        let (ncb, kcb) = (nf.min(2048.0), 256.0f64);
        let bytes_per_elem = 8.0 * (nf * (nf / ncb).ceil() / nf + 1.0 + 2.0 * (nf / kcb).ceil());
        let spec = KernelSpec {
            name: "gemm-f64-parallel",
            version: "v3-worksteal",
            bytes_per_elem,
            flops_per_elem: 2.0 * nf,
            threading: Threading::AllCore,
            target_axis: TargetAxis::BindingRoof,
            target_fraction: None,
        };
        let elems_per_sec = g * 1e9 / (2.0 * nf);
        let att = attainment_for(&spec, elems_per_sec, &axes);
        println!(
            "{{\"metric\":\"gemm-f64-parallel\",\"n\":{n},\"gflops\":{g:.2},\"roof\":\"{:?}\",\"attainment_minroof\":{:.3}}}",
            att.roof, att.attainment
        );
        if n == 2048 && std::env::var("FS_LA_ROOFLINE_GATE").as_deref() == Ok("1") {
            assert!(
                att.attainment >= 0.5,
                "all-core GEMM min-roof attainment {:.3} below the 50% floor",
                att.attainment
            );
        }
    }
}

/// Compute tonight's min-roof attainment at one n (the ledgered
/// scalar), returning (attainment, wall seconds).
fn nightly_attainment(n: usize, threads: usize, axes: &MachineAxes) -> (f64, f64) {
    let t0 = std::time::Instant::now();
    let g = measure_parallel(n, 1, threads);
    let wall = t0.elapsed().as_secs_f64();
    let nf = n as f64;
    let (ncb, kcb) = (nf.min(2048.0), 256.0f64);
    let bytes_per_elem = 8.0 * (nf * (nf / ncb).ceil() / nf + 1.0 + 2.0 * (nf / kcb).ceil());
    let spec = KernelSpec {
        name: "gemm-f64-parallel",
        version: "v3-worksteal",
        bytes_per_elem,
        flops_per_elem: 2.0 * nf,
        threading: Threading::AllCore,
        target_axis: TargetAxis::BindingRoof,
        target_fraction: None,
    };
    (
        attainment_for(&spec, g * 1e9 / (2.0 * nf), axes).attainment,
        wall,
    )
}

/// Pull `"key":<f64>` out of one ledger line (own fixed format).
fn ledger_f64(line: &str, key: &str) -> f64 {
    let pat = format!("\"{key}\":");
    let start = line
        .find(&pat)
        .unwrap_or_else(|| panic!("ledger line missing {key}: {line}"))
        + pat.len();
    let rest = &line[start..];
    let end = rest
        .find([',', '}'])
        .unwrap_or_else(|| panic!("unterminated {key} in ledger line: {line}"));
    rest[..end]
        .parse::<f64>()
        .unwrap_or_else(|e| panic!("bad {key} in ledger line ({e}): {line}"))
}

/// NIGHTLY PERF-REGRESSION wiring (xlvx item 6 / fz2.4 machinery):
/// each run measures the n=2048 all-core min-roof attainment, appends
/// it to the JSONL ledger at FS_LA_REGRESS_LEDGER (with per-n wall
/// times as the phase stream for attribution), then runs the
/// dispersion-aware GateSpec band over the ledgered nights and the
/// CUSUM drift detector over the expanding-baseline z-scores. A Red
/// verdict or a CUSUM alarm IS a test failure (perf regressions are
/// test failures — plan §14.4); Invalid evidence fails loudest.
/// Without the env var the lane reports and skips: the ledger's HOME
/// (per-machine persistent path) is the nightly CI's to choose.
#[test]
#[ignore = "perf lane: run explicitly in release with --ignored"]
fn gemm_regress_nightly() {
    let Ok(ledger) = std::env::var("FS_LA_REGRESS_LEDGER") else {
        println!(
            "{{\"metric\":\"gemm-regress\",\"verdict\":\"skip\",\"detail\":\"FS_LA_REGRESS_LEDGER unset; nightly CI owns the ledger path\"}}"
        );
        return;
    };
    let threads = std::env::var("FS_LA_THREADS")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .unwrap_or_else(|| std::thread::available_parallelism().map_or(8, std::num::NonZero::get));
    let axes = MachineAxes::probe();
    // Phases: per-n wall shares. Coarse but honest attribution — a Red
    // whose n512 share grew is per-chunk overhead; n2048 is the kernel.
    let (_, w512) = nightly_attainment(512, threads, &axes);
    let (_, w1024) = nightly_attainment(1024, threads, &axes);
    let (att, w2048) = nightly_attainment(2048, threads, &axes);
    let prior = std::fs::read_to_string(&ledger).unwrap_or_default();
    let night = prior.lines().filter(|l| !l.trim().is_empty()).count() as u64;
    let row = format!(
        "{{\"night\":{night},\"attainment\":{att},\"phases\":{{\"n512\":{w512},\"n1024\":{w1024},\"n2048\":{w2048}}}}}\n"
    );
    std::fs::write(&ledger, prior.clone() + &row).expect("append to FS_LA_REGRESS_LEDGER");
    // Reload the full history into regress::Night rows.
    let history: Vec<Night> = (prior + &row)
        .lines()
        .filter(|l| !l.trim().is_empty())
        .map(|line| Night {
            night: ledger_f64(line, "night") as u64,
            attainment: ledger_f64(line, "attainment"),
            phases: [
                ("n512".to_string(), ledger_f64(line, "n512")),
                ("n1024".to_string(), ledger_f64(line, "n1024")),
                ("n2048".to_string(), ledger_f64(line, "n2048")),
            ]
            .into(),
        })
        .collect();
    let verdict = gate(&history, GateSpec::default());
    let attainments: Vec<f64> = history.iter().map(|n| n.attainment).collect();
    let standardized = standardize(&attainments, GateSpec::default().min_baseline)
        .expect("the retained regression ledger is bounded");
    let alarm = Cusum::default().first_alarm(&standardized);
    println!(
        "{{\"metric\":\"gemm-regress\",\"nights\":{},\"tonight\":{att:.3},\"gate\":\"{}\",\"cusum_alarm\":{}}}",
        history.len(),
        match &verdict {
            GateVerdict::Green { z } => format!("green z={z:.2}"),
            GateVerdict::Red { z, .. } => format!("RED z={z:.2}"),
            GateVerdict::Invalid { reason } => format!("INVALID {reason}"),
        },
        alarm.map_or_else(|| "null".to_string(), |i| i.to_string()),
    );
    match verdict {
        GateVerdict::Green { .. } => {}
        GateVerdict::Red { z, attribution } => panic!(
            "nightly GEMM attainment regressed: z={z:.2}, top offender phase: {:?}",
            attribution.first()
        ),
        GateVerdict::Invalid { reason } => panic!("malformed regress evidence: {reason}"),
    }
    assert!(
        alarm.is_none(),
        "CUSUM drift alarm at night index {} — slow regression the per-night gate missed",
        alarm.unwrap_or_default()
    );
}

/// The MC/NC AUTOTUNE SWEEP (xlvx segment 5): report-only rows over the
/// bit-neutral blocking grid at n = 2048, all cores. "shipping" is what
/// gemm_f64_parallel actually ships (MC_PAR/NC_PAR_CAP); the grid brackets it.
/// Feeds the tuned-defaults decision — KC is NOT swept here (bit
/// contract; retuning it is a golden bump with justification).
#[test]
#[ignore = "perf lane: run explicitly in release with --ignored"]
fn gemm_tune_sweep() {
    let threads = std::env::var("FS_LA_THREADS")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .unwrap_or_else(|| std::thread::available_parallelism().map_or(8, std::num::NonZero::get));
    let n = 2048usize;
    let a: Vec<f64> = (0..n * n).map(|i| ((i as f64) * 0.13).sin()).collect();
    let b: Vec<f64> = (0..n * n).map(|i| ((i as f64) * 0.31).cos()).collect();
    let mut c = vec![0.0f64; n * n];
    let mut measure_with = |mc: usize, nc: usize| -> f64 {
        fs_la::gemm_f64_parallel_with(n, n, n, 1.0, &a, &b, 0.0, &mut c, threads, mc, nc); // warm
        let mut best = f64::INFINITY;
        for _ in 0..3 {
            let t0 = std::time::Instant::now();
            fs_la::gemm_f64_parallel_with(n, n, n, 1.0, &a, &b, 0.0, &mut c, threads, mc, nc);
            best = best.min(t0.elapsed().as_secs_f64());
        }
        2.0 * (n * n * n) as f64 / best / 1e9
    };
    for mc in [16usize, 32, 64, 128] {
        for nc in [256usize, 512, 1024, 2048] {
            let g = measure_with(mc, nc);
            println!(
                "{{\"metric\":\"gemm-tune\",\"threads\":{threads},\"mc\":{mc},\"nc\":{nc},\"gflops\":{g:.2}}}"
            );
        }
    }
    // The shipping-defaults row, for comparison against the grid.
    let g = {
        gemm_f64_parallel(n, n, n, 1.0, &a, &b, 0.0, &mut c, threads);
        let mut best = f64::INFINITY;
        for _ in 0..3 {
            let t0 = std::time::Instant::now();
            gemm_f64_parallel(n, n, n, 1.0, &a, &b, 0.0, &mut c, threads);
            best = best.min(t0.elapsed().as_secs_f64());
        }
        2.0 * (n * n * n) as f64 / best / 1e9
    };
    println!(
        "{{\"metric\":\"gemm-tune\",\"threads\":{threads},\"mc\":\"shipping\",\"nc\":\"shipping\",\"gflops\":{g:.2}}}"
    );
}
