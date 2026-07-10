//! The fs-fft PERF LANE (bead 27d3): mixed radix-4/2 Stockham
//! throughput against the MEMORY-BOUND roofline (fs-substrate STREAM
//! triad via fs-roofline axes — the plan's denominator for this
//! kernel), ≥40% attainment gate at memory-resident sizes. Run
//! explicitly in release:
//! `cargo test -p fs-fft --release --test perf_lane -- --ignored --nocapture`
//!
//! One `run_once` is a forward+inverse ROUND TRIP (keeps values
//! bounded across repetitions); the byte model counts every Stockham
//! pass (32 B/element each), ping-pong copy-back passes, and the
//! inverse's 1/n scale pass — the honest traffic of THIS algorithm,
//! not a compulsory-miss fantasy.

use fs_fft::{C64, Fft};
use fs_roofline::{KernelSpec, MachineAxes, RooflineKernel, Threading, measure};

/// Stockham stage count for the mixed radix-4/2 formulation.
fn stages(n: usize) -> usize {
    let mut c = 0;
    let mut m = n;
    while m >= 4 {
        m /= 4;
        c += 1;
    }
    if m == 2 {
        c += 1;
    }
    c
}

struct FftRoundTrip {
    n: usize,
    plan: Fft,
    data: Vec<C64>,
    scratch: Vec<C64>,
}

impl FftRoundTrip {
    fn new(n: usize) -> FftRoundTrip {
        FftRoundTrip {
            n,
            plan: Fft::new(n),
            data: (0..n)
                .map(|i| C64::new(((i * 37) % 101) as f64 * 0.02 - 1.0, ((i * 53) % 97) as f64 * 0.02))
                .collect(),
            scratch: vec![C64::new(0.0, 0.0); n],
        }
    }
}

impl RooflineKernel for FftRoundTrip {
    fn spec(&self) -> KernelSpec {
        let st = stages(self.n) as f64;
        let copy = if stages(self.n) % 2 == 1 { 32.0 } else { 0.0 };
        KernelSpec {
            name: "fft-roundtrip",
            version: "27d3-r4",
            // Two transforms of `st` passes (32 B/elem each: read one
            // C64, write one C64) + copy-back per transform when the
            // stage count is odd + the inverse's scale pass.
            bytes_per_elem: 2.0 * (32.0 * st + copy) + 32.0,
            // Radix-4 butterfly ≈ 34 flops / 4 outputs = 8.5 per
            // element-stage; + 2 for the scale. Approximate — the roof
            // is bandwidth at this intensity either way.
            flops_per_elem: 2.0 * 8.5 * st + 2.0,
            threading: Threading::SingleThread,
            target_fraction: Some(0.40),
        }
    }
    fn elements(&self) -> usize {
        self.n
    }
    fn run_once(&mut self) {
        self.plan.forward(&mut self.data, &mut self.scratch);
        self.plan.inverse(&mut self.data, &mut self.scratch);
    }
}

#[test]
#[ignore = "perf lane: run explicitly in release with --ignored"]
fn fft_attainment() {
    let axes = MachineAxes::probe();
    println!("{}", axes.to_jsonl());
    // Size ladder: L2-resident (2^16) reported for context; the GATE
    // reads the memory-resident sizes (2^20, 2^22 — 32/128 MB working
    // sets against the DRAM STREAM axis).
    let mut gate_ok = true;
    for &(n, gated) in &[(1usize << 16, false), (1 << 20, true), (1 << 22, true)] {
        let mut kern = FftRoundTrip::new(n);
        let att = measure(&mut kern, 1, 5, &axes);
        println!(
            "{{\"metric\":\"fft-roundtrip\",\"n\":{n},\"gated\":{gated},{}}}",
            att.to_jsonl().trim_start_matches('{')
        );
        if gated {
            gate_ok &= att.attainment >= 0.40;
        }
    }
    println!(
        "{{\"metric\":\"fft-gate\",\"floor\":0.40,\"machine\":\"{}-{}\",\"pass\":{gate_ok}}}",
        std::env::consts::OS,
        std::env::consts::ARCH
    );
    assert!(
        gate_ok,
        "memory-resident FFT round trips clear 40% of the STREAM-bound roofline"
    );
}
