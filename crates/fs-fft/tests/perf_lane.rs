//! The fs-fft PERF LANE (bead 27d3): mixed radix-8/4/2 Stockham
//! throughput against the MEMORY-BOUND roofline (fs-substrate STREAM
//! triad via fs-roofline axes — the plan's denominator for this
//! kernel). The ≥40% plan target is reported honestly; until it lands,
//! the executable regression gate is a 15% anti-collapse floor at
//! memory-resident sizes. Run explicitly in release:
//! `cargo test -p fs-fft --release --test perf_lane -- --ignored --nocapture`
//!
//! One `run_once` is a forward+inverse ROUND TRIP (keeps values
//! bounded across repetitions); the byte model counts every Stockham
//! pass (32 B/element each), ping-pong copy-back passes, and the
//! inverse's 1/n scale pass — the honest traffic of THIS algorithm,
//! not a compulsory-miss fantasy.
//!
//! CITABLE GATES REQUIRE A TRUSTED BASELINE (beads dfh3/fz2.7, same
//! doctrine as the fs-feec lane): run with
//! `FRANKENSIM_BASELINE_STORE=<jsonl> FRANKENSIM_FIRMWARE_ID=<id>` —
//! the pre/post axes must be admitted against the promoted baseline
//! for this machine fingerprint or the lane FAILS with the receipt.

use fs_fft::{C64, Fft};
use fs_roofline::{
    AxisBaselinePolicy, BaselineIdentity, BaselineStore, KernelSpec, MachineAxes, RooflineKernel,
    TargetAxis, Threading, days_since_epoch_now, measure,
};

/// Stockham stage count for the mixed radix-8/4/2 formulation — MUST
/// mirror the transform's decomposition exactly or the traffic model
/// (and hence attainment) lies.
fn stages(n: usize) -> usize {
    let mut c = 0;
    let mut m = n;
    while m >= 8 {
        m /= 8;
        c += 1;
    }
    if m >= 2 {
        c += 1; // one radix-4 or radix-2 residue stage
    }
    c
}

/// Does `n` take the six-step path? MUST mirror `Fft::takes_sixstep`
/// (feature-gated; n ≥ 2^16 with even log₂). The default lane models
/// and measures the stage walk; enabling `frontier-sixstep` flips both
/// the kernel and this model together.
fn takes_sixstep(n: usize) -> bool {
    cfg!(feature = "frontier-sixstep") && n >= (1 << 16) && n.trailing_zeros().is_multiple_of(2)
}

/// Full-array DRAM passes per single transform: the six-step does six
/// (three out-of-place transposes + two row sweeps whose sub-transforms
/// are cache-resident + the final copy-back); the stage walk does one
/// per stage plus the odd-parity copy-back.
fn dram_passes(n: usize) -> f64 {
    if takes_sixstep(n) {
        6.0
    } else {
        let st = stages(n);
        st as f64 + f64::from(u8::from(st % 2 == 1))
    }
}

/// Butterfly element-stages actually executed (for the flop model):
/// the six-step runs the sub-plan (√n) twice per transform.
fn butterfly_stages(n: usize) -> f64 {
    if takes_sixstep(n) {
        let n1 = 1usize << (n.trailing_zeros() / 2);
        2.0 * butterfly_stages(n1)
    } else {
        stages(n) as f64
    }
}

fn measurement_json(n: usize, gated: bool, receipt: &str) -> String {
    let receipt_fields = receipt
        .strip_prefix('{')
        .and_then(|fields| fields.strip_suffix('}'))
        .expect("roofline attainment receipt must be a JSON object");
    format!("{{\"metric\":\"fft-roundtrip\",\"n\":{n},\"gated\":{gated},{receipt_fields}}}")
}

#[test]
fn measurement_receipt_is_one_json_object() {
    assert_eq!(
        measurement_json(16, true, "{\"schema\":\"attainment-v1\",\"value\":1}"),
        "{\"metric\":\"fft-roundtrip\",\"n\":16,\"gated\":true,\"schema\":\"attainment-v1\",\"value\":1}"
    );
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
                .map(|i| {
                    C64::new(
                        ((i * 37) % 101) as f64 * 0.02 - 1.0,
                        ((i * 53) % 97) as f64 * 0.02,
                    )
                })
                .collect(),
            scratch: vec![C64::new(0.0, 0.0); n],
        }
    }
}

impl RooflineKernel for FftRoundTrip {
    fn spec(&self) -> KernelSpec {
        let passes = dram_passes(self.n);
        let bf = butterfly_stages(self.n);
        // Six-step adds one fused complex twiddle multiply per element
        // per transform (6 flops).
        let twiddle = if takes_sixstep(self.n) { 6.0 } else { 0.0 };
        KernelSpec {
            name: "fft-roundtrip",
            version: if takes_sixstep(self.n) {
                "27d3-6s"
            } else {
                "27d3-r8"
            },
            // Two transforms of `passes` full-array DRAM passes
            // (32 B/elem each: read one C64, write one C64) + the
            // inverse's scale pass.
            bytes_per_elem: 2.0 * 32.0 * passes + 32.0,
            // Radix-8 butterfly ≈ 100 flops / 8 outputs = 12.5 per
            // element-stage; + 2 for the scale. Approximate — the roof
            // is bandwidth at this intensity either way.
            flops_per_elem: 2.0 * (12.5 * bf + twiddle) + 2.0,
            threading: Threading::SingleThread,
            target_axis: TargetAxis::BindingRoof,
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
    // Environment validity (bead 1n61): implausible axes poison both
    // the numerator and denominator of attainment — refuse up front.
    if let Some(reason) = axes.plausibility_error() {
        println!(
            "{{\"metric\":\"fft-gate\",\"verdict\":\"environment_invalid\",\
             \"reason\":\"{reason}\",\"machine\":\"{}-{}\"}}",
            std::env::consts::OS,
            std::env::consts::ARCH
        );
        panic!("FFT roofline evidence rejected: {reason}");
    }
    // Trusted-baseline admission (beads dfh3/fz2.7): a citable row needs
    // the promoted quiet-machine baseline for THIS fingerprint — static
    // floors and pre/post self-agreement cannot detect a host that is
    // consistently degraded through the whole run.
    let baseline_path = std::env::var("FRANKENSIM_BASELINE_STORE")
        .unwrap_or_else(|_| panic!("FRANKENSIM_BASELINE_STORE is required for a citable gate"));
    let firmware = std::env::var("FRANKENSIM_FIRMWARE_ID")
        .unwrap_or_else(|_| panic!("FRANKENSIM_FIRMWARE_ID is required for a citable gate"));
    let baseline_text = std::fs::read_to_string(&baseline_path)
        .unwrap_or_else(|error| panic!("cannot read baseline store {baseline_path:?}: {error}"));
    let baseline_store = BaselineStore::from_jsonl(&baseline_text)
        .unwrap_or_else(|error| panic!("invalid baseline store: {error}"));
    let identity = BaselineIdentity::current(&axes, firmware)
        .unwrap_or_else(|error| panic!("invalid baseline identity: {error}"));
    let now_day = days_since_epoch_now()
        .unwrap_or_else(|error| panic!("cannot establish baseline age: {error}"));
    let baseline_policy = AxisBaselinePolicy::new(
        baseline_store.for_fingerprint(axes.fingerprint),
        &identity,
        now_day,
    );
    // Size ladder: L2-resident (2^16) reported for context; the gate
    // rows are the memory-resident sizes (2^20, 2^22 — 32/128 MB
    // working sets against the DRAM STREAM axis).
    let mut gate_ok = true;
    let mut floor_ok = true;
    let mut env_ok = true;
    for &(n, gated) in &[(1usize << 16, false), (1 << 20, true), (1 << 22, true)] {
        let mut kern = FftRoundTrip::new(n);
        let att = measure(&mut kern, 1, 5, &axes);
        let receipt = att.to_jsonl();
        let measurement = measurement_json(n, gated, &receipt);
        println!("{measurement}");
        // An environment-invalid row contributes neither a target pass nor a
        // numerical regression failure. It poisons the evidence lane as a
        // whole below, so cargo cannot report a green citable run.
        if att.verdict == fs_roofline::Verdict::EnvironmentInvalid {
            env_ok = false;
            continue;
        }
        if gated {
            gate_ok &= att.attainment >= 0.40;
            floor_ok &= att.attainment >= 0.15;
        }
    }
    if !env_ok {
        println!(
            "{{\"metric\":\"fft-gate\",\"verdict\":\"environment_invalid\",             \"machine\":\"{}-{}\"}}",
            std::env::consts::OS,
            std::env::consts::ARCH
        );
        panic!("FFT roofline evidence rejected: contaminated environment");
    }
    // Post-run reprobe + baseline verdict: the run is citable only if
    // BOTH probes are admitted against the trusted baseline (drift
    // during the run, or a consistently-degraded host, both refuse).
    let post_axes = MachineAxes::probe();
    println!("{{\"metric\":\"axes-post\",\"axes\":{}}}", post_axes.to_jsonl());
    let baseline_verdict = baseline_policy.verdict(&axes, &post_axes);
    println!("{}", baseline_policy.receipt_json(&axes, &post_axes));
    if !baseline_verdict.trusted() {
        println!(
            "{{\"metric\":\"fft-gate\",\"verdict\":\"environment_invalid\",\
             \"reason\":\"historical baseline admission rejected\",\"machine\":\"{}-{}\"}}",
            std::env::consts::OS,
            std::env::consts::ARCH
        );
        panic!("FFT roofline evidence rejected: historical baseline admission rejected");
    }
    // The 0.40 target is REPORTED per row; measured 0.26–0.43 across
    // runs on this machine, dominated by axis and load noise from
    // concurrent agent builds (bead 27d3 records the numbers). The
    // ASSERTED gate is the anti-collapse floor.
    println!(
        "{{\"metric\":\"fft-gate\",\"target\":0.40,\"target_met\":{gate_ok},\
         \"floor\":0.15,\"floor_met\":{floor_ok},\"machine\":\"{}-{}\"}}",
        std::env::consts::OS,
        std::env::consts::ARCH
    );
    assert!(
        floor_ok,
        "memory-resident FFT round trips collapsed below the 15% floor"
    );
}
