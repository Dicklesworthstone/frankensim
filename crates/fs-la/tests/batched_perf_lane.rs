//! The fs-la BATCHED PERF LANE (bead 9ekv): batch_gemm attainment per
//! size class {4, 6, 8, 12, 16, 24, 32, 48} against the machine
//! ROOFLINE (fs-roofline conventions: limit = min(bandwidth·intensity,
//! compute). Each row reports both binding-roof `attainment` and
//! compute-peak `target_attainment`: the plan's 60% target uses the latter,
//! while the executable anti-collapse floor uses the former. See the note at
//! the bottom and bead 9ekv for the measured achieved-vs-target gap. Run
//! explicitly in release:
//! `cargo test -p fs-la --release --test batched_perf_lane -- --ignored --nocapture`
//!
//! Batch sizes put the working set at ~50 MB (memory-resident, the
//! FEM-assembly regime the layout doctrine targets). LU is reported
//! (flop model documented), gated only against pathological collapse.

use fs_la::batched::{BatchMat, batch_gemm, batch_lu};
use fs_roofline::{KernelSpec, MachineAxes, RooflineKernel, TargetAxis, Threading, measure};

fn measurement_json(metric: &str, k: usize, n: usize, receipt: &str) -> String {
    let receipt_fields = receipt
        .strip_prefix('{')
        .and_then(|fields| fields.strip_suffix('}'))
        .expect("roofline attainment receipt must be a JSON object");
    format!("{{\"metric\":\"{metric}\",\"k\":{k},\"n\":{n},{receipt_fields}}}")
}

#[test]
fn measurement_receipt_is_one_json_object() {
    assert_eq!(
        measurement_json(
            "batch-gemm",
            4,
            16,
            "{\"schema\":\"attainment-v1\",\"attainment\":0.5}"
        ),
        "{\"metric\":\"batch-gemm\",\"k\":4,\"n\":16,\"schema\":\"attainment-v1\",\"attainment\":0.5}"
    );
}

struct BatchGemmKernel {
    k: usize,
    a: BatchMat,
    b: BatchMat,
    c: BatchMat,
}

impl BatchGemmKernel {
    fn new(k: usize) -> BatchGemmKernel {
        let n = ((2usize << 20) / (k * k)).max(256);
        let f = |m: usize, i: usize, j: usize| ((m * 31 + i * 7 + j) % 17) as f64 * 0.125 - 1.0;
        BatchGemmKernel {
            k,
            a: BatchMat::from_fn(k, n, f),
            b: BatchMat::from_fn(k, n, |m, i, j| f(m + 3, j, i)),
            c: BatchMat::zeros(k, n),
        }
    }
}

impl RooflineKernel for BatchGemmKernel {
    fn spec(&self) -> KernelSpec {
        let kf = self.k as f64;
        KernelSpec {
            name: "batch-gemm",
            version: "9ekv",
            // Compulsory traffic per matrix: read A and B, write C
            // once (chunk-resident accumulator/planes, MBLK doctrine).
            bytes_per_elem: 3.0 * kf * kf * 8.0,
            flops_per_elem: 2.0 * kf * kf * kf,
            threading: Threading::SingleThread,
            target_axis: TargetAxis::ComputePeak,
            target_fraction: Some(0.60),
        }
    }
    fn elements(&self) -> usize {
        self.a.batch_len()
    }
    fn run_once(&mut self) {
        batch_gemm(1.0, &self.a, &self.b, 0.0, &mut self.c);
    }
}

struct BatchLuKernel {
    a: BatchMat,
}

impl RooflineKernel for BatchLuKernel {
    fn spec(&self) -> KernelSpec {
        let kf = self.a.k() as f64;
        KernelSpec {
            name: "batch-lu",
            version: "9ekv",
            // clone(A) + factor in place: ~3k² compulsory + k²/2 pivot
            // rescans; modeled as 4k² (documented approximation).
            bytes_per_elem: 4.0 * kf * kf * 8.0,
            // ~(2/3)k³ multiply-adds = (4/3)k³ flops + k² divides.
            flops_per_elem: 4.0 / 3.0 * kf * kf * kf,
            threading: Threading::SingleThread,
            target_axis: TargetAxis::BindingRoof,
            target_fraction: None, // reported; collapse-gated below
        }
    }
    fn elements(&self) -> usize {
        self.a.batch_len()
    }
    fn run_once(&mut self) {
        let out = batch_lu(&self.a);
        assert!(out.flags.is_empty(), "perf fixture must be nonsingular");
        std::hint::black_box(out.lu.get(0, 0, 0));
    }
}

#[test]
#[ignore = "perf lane: run explicitly in release with --ignored"]
fn batched_attainment() {
    let axes = MachineAxes::probe();
    println!("{}", axes.to_jsonl());
    let mut all_within = true;
    let mut floor_ok = true;
    let mut environment_valid = true;
    for &k in &[4usize, 6, 8, 12, 16, 24, 32, 48] {
        let mut kern = BatchGemmKernel::new(k);
        let att = measure(&mut kern, 1, 5, &axes);
        let receipt = measurement_json("batch-gemm", k, kern.elements(), &att.to_jsonl());
        println!("{receipt}");
        if att.verdict == fs_roofline::Verdict::EnvironmentInvalid {
            environment_valid = false;
            continue;
        }
        all_within &= att.target_attainment >= 0.60;
        floor_ok &= att.attainment >= 0.08;
    }
    // LU report rows (diagonally-dominant fixture, flag-free).
    for &k in &[4usize, 8, 16, 32] {
        let n = ((1usize << 20) / (k * k)).max(256);
        let a = BatchMat::from_fn(k, n, |m, i, j| {
            let base = ((m * 13 + i * 3 + j * 11) % 23) as f64 * 0.0625 - 0.7;
            if i == j { base + 3.0 * k as f64 } else { base }
        });
        let mut kern = BatchLuKernel { a };
        let att = measure(&mut kern, 1, 5, &axes);
        let receipt = measurement_json("batch-lu", k, n, &att.to_jsonl());
        println!("{receipt}");
        if att.verdict == fs_roofline::Verdict::EnvironmentInvalid {
            environment_valid = false;
            continue;
        }
        assert!(
            att.attainment >= 0.05,
            "batch-lu k={k} collapsed: attainment {:.3}",
            att.attainment
        );
    }
    // The 60% target is REPORTED per row (verdict field) but not yet
    // met on this machine: the plane-SoA lane walk is load-port/TLB
    // bound near 10-26 GFLOP/s depending on k (measured; the 4×4-tile
    // capsule already removed the accumulator round-trips). The
    // achieved-vs-target gap and the successor design notes live in
    // bead 9ekv. The ASSERTED gate here is the anti-collapse floor —
    // a regression, not an aspiration.
    let target_met = environment_valid && all_within;
    let floor_met = environment_valid && floor_ok;
    println!(
        "{{\"metric\":\"batched-gate\",\"target_axis\":\"compute_peak\",\
         \"target\":0.60,\"target_met\":{target_met},\
         \"floor_axis\":\"binding_roof\",\"floor\":0.08,\"floor_met\":{floor_met},\
         \"environment_valid\":{environment_valid},\"machine\":\"{}-{}\"}}",
        std::env::consts::OS,
        std::env::consts::ARCH
    );
    assert!(
        environment_valid,
        "batched roofline evidence rejected: contaminated environment"
    );
    assert!(
        floor_ok,
        "batch_gemm attainment collapsed below the 8% floor"
    );
}
