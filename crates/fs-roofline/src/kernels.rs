//! Built-in registry kernels: fs-simd primitives plus the production
//! session-autotuned f64 GEMM route, and a deliberately de-optimized kernel
//! for harness meta-tests.
//!
//! Intensity models (the hand-calculation basis for `rf_002`):
//! - axpy `y = a·x + y`: reads x and y, writes y → 24 B/elem, 2 flop/elem.
//! - dot `Σ x·y`: reads x and y → 16 B/elem, 2 flop/elem.
//! - sum `Σ x`: reads x → 8 B/elem, 1 flop/elem.
//! - GEMM `C = A*B`: minimum traffic A+B+C → 24 B/output element for square
//!   matrices, 2k flop/output element. The actual timed path is
//!   `fs_session::gemm_f64_session`, so warmup closes measure → cache → model
//!   → dispatch and repetitions reuse the same validated tune row.
//!
//! Targets here are report-only in v0 (`target_fraction: None` except where
//! a band is deliberately claimed for meta-testing): CI runners are shared
//! machines and §14 bands belong to fingerprinted reference machines.

use crate::{KernelSpec, RooflineKernel, Threading};

/// Roofline wrapper version. The lower-layer implementation/tier/placement
/// identities are independently bound inside fs-session's tune key.
pub const GEMM_ROOFLINE_VERSION: &str = "1";

/// Production f64 GEMM benchmark routed through the session autotuner.
///
/// The kernel owns its [`fs_exec::Tuner`] and [`fs_exec::CancelGate`]. The
/// roofline registry invokes kernels sequentially through exclusive `&mut`
/// borrows, so tune/cache/dispatch state needs no wrapper lock. A failure is
/// fail-closed: [`RooflineKernel::run_once`] cannot return `Result`, therefore
/// an unexpected session diagnostic panics instead of emitting a normal-
/// looking timing row.
pub struct GemmKernel {
    m: usize,
    n: usize,
    k: usize,
    threads: usize,
    a: Vec<f64>,
    b: Vec<f64>,
    c: Vec<f64>,
    tuner: fs_exec::Tuner,
    tune_ledger: Option<fs_ledger::Ledger>,
    gate: fs_exec::CancelGate,
    dispatches: usize,
    sweeps: usize,
}

impl GemmKernel {
    /// A square production GEMM sized by one matrix edge.
    ///
    /// # Panics
    /// If the requested matrix extents overflow `usize`.
    #[must_use]
    pub fn square(side: usize, threads: usize, machine_fingerprint: u64) -> Self {
        Self::new(side, side, side, threads, machine_fingerprint, None)
    }

    #[track_caller]
    fn new(
        m: usize,
        n: usize,
        k: usize,
        threads: usize,
        machine_fingerprint: u64,
        tune_ledger: Option<fs_ledger::Ledger>,
    ) -> Self {
        assert!(
            m > 0 && n > 0 && k > 0,
            "roofline GEMM extents must be positive"
        );
        let a_len = m.checked_mul(k).expect("GEMM A extent overflow");
        let b_len = k.checked_mul(n).expect("GEMM B extent overflow");
        let c_len = m.checked_mul(n).expect("GEMM C extent overflow");
        let a = (0..a_len)
            .map(|i| ((i % 31) as f64 - 15.0) / 31.0)
            .collect();
        let b = (0..b_len)
            .map(|i| ((i % 29) as f64 - 14.0) / 29.0)
            .collect();
        Self {
            m,
            n,
            k,
            threads: threads.max(1),
            a,
            b,
            c: vec![0.0; c_len],
            tuner: fs_exec::Tuner::cold(machine_fingerprint),
            tune_ledger,
            gate: fs_exec::CancelGate::new(),
            dispatches: 0,
            sweeps: 0,
        }
    }

    /// Number of completed session dispatches (warmups included).
    #[must_use]
    pub fn dispatches(&self) -> usize {
        self.dispatches
    }

    /// Number of calls that performed the bounded measurement sweep. A stable
    /// kernel instance should report one after its cold first invocation.
    #[must_use]
    pub fn sweeps(&self) -> usize {
        self.sweeps
    }
}

impl RooflineKernel for GemmKernel {
    fn spec(&self) -> KernelSpec {
        KernelSpec {
            name: "gemm-f64",
            version: GEMM_ROOFLINE_VERSION,
            // Square production instances have one A, B, and C matrix. The
            // rectangular constructor exists only for the bounded regression.
            bytes_per_elem: 8.0 * (self.a.len() as f64 + self.b.len() as f64 + self.c.len() as f64)
                / self.c.len() as f64,
            flops_per_elem: 2.0 * self.k as f64,
            threading: Threading::AllCore,
            target_fraction: Some(0.75),
        }
    }

    fn elements(&self) -> usize {
        self.m * self.n
    }

    fn run_once(&mut self) {
        let dispatch = fs_session::gemm_f64_session(
            &mut self.tuner,
            self.tune_ledger.as_ref(),
            &self.gate,
            self.threads,
            self.m,
            self.n,
            self.k,
            1.0,
            &self.a,
            &self.b,
            0.0,
            &mut self.c,
        )
        .unwrap_or_else(|error| panic!("production roofline GEMM dispatch failed: {error}"));
        self.dispatches += 1;
        self.sweeps += usize::from(dispatch.swept);
        std::hint::black_box(self.c[self.c.len() / 2]);
    }
}

/// `y = a·x + y` through the dispatched fs-simd table.
pub struct AxpyKernel {
    x: Vec<f64>,
    y: Vec<f64>,
}

impl AxpyKernel {
    /// Buffers of `n` elements each (pick `n` large enough to stream past
    /// the last-level cache when measuring the bandwidth roof).
    #[must_use]
    pub fn new(n: usize) -> AxpyKernel {
        AxpyKernel {
            x: vec![1.5; n],
            y: vec![0.5; n],
        }
    }
}

impl RooflineKernel for AxpyKernel {
    fn spec(&self) -> KernelSpec {
        KernelSpec {
            name: "simd-axpy-f64",
            version: "1",
            bytes_per_elem: 24.0,
            flops_per_elem: 2.0,
            threading: Threading::SingleThread,
            target_fraction: None,
        }
    }

    fn elements(&self) -> usize {
        self.x.len()
    }

    fn run_once(&mut self) {
        (fs_simd::ops().axpy)(1.000_000_1, &self.x, &mut self.y);
        std::hint::black_box(self.y[self.y.len() / 2]);
    }
}

/// `Σ x·y` through the dispatched fs-simd table.
pub struct DotKernel {
    x: Vec<f64>,
    y: Vec<f64>,
    out: f64,
}

impl DotKernel {
    /// Buffers of `n` elements each.
    #[must_use]
    pub fn new(n: usize) -> DotKernel {
        DotKernel {
            x: vec![1.5; n],
            y: vec![0.5; n],
            out: 0.0,
        }
    }
}

impl RooflineKernel for DotKernel {
    fn spec(&self) -> KernelSpec {
        KernelSpec {
            name: "simd-dot-f64",
            version: "1",
            bytes_per_elem: 16.0,
            flops_per_elem: 2.0,
            threading: Threading::SingleThread,
            target_fraction: None,
        }
    }

    fn elements(&self) -> usize {
        self.x.len()
    }

    fn run_once(&mut self) {
        self.out = (fs_simd::ops().dot)(&self.x, &self.y);
        std::hint::black_box(self.out);
    }
}

/// `Σ x` through the dispatched fs-simd table.
pub struct SumKernel {
    x: Vec<f64>,
    out: f64,
}

impl SumKernel {
    /// A buffer of `n` elements.
    #[must_use]
    pub fn new(n: usize) -> SumKernel {
        SumKernel {
            x: vec![0.25; n],
            out: 0.0,
        }
    }
}

impl RooflineKernel for SumKernel {
    fn spec(&self) -> KernelSpec {
        KernelSpec {
            name: "simd-sum-f64",
            version: "1",
            bytes_per_elem: 8.0,
            flops_per_elem: 1.0,
            threading: Threading::SingleThread,
            target_fraction: None,
        }
    }

    fn elements(&self) -> usize {
        self.x.len()
    }

    fn run_once(&mut self) {
        self.out = (fs_simd::ops().sum)(&self.x);
        std::hint::black_box(self.out);
    }
}

/// Deliberately de-optimized kernel with a band it cannot meet: a serial
/// dependency chain strided across a buffer, claiming 90% of the bandwidth
/// roof. The harness meta-test asserts it reports `BelowBand` — proof that
/// a slow kernel is caught, not absorbed (bead acceptance criterion).
pub struct SeededSlowKernel {
    x: Vec<f64>,
    out: f64,
}

impl SeededSlowKernel {
    /// A buffer of `n` elements.
    #[must_use]
    pub fn new(n: usize) -> SeededSlowKernel {
        SeededSlowKernel {
            x: (0..n).map(|i| (i % 7) as f64).collect(),
            out: 0.0,
        }
    }
}

impl RooflineKernel for SeededSlowKernel {
    fn spec(&self) -> KernelSpec {
        KernelSpec {
            name: "seeded-slow",
            version: "1",
            bytes_per_elem: 8.0,
            flops_per_elem: 1.0,
            threading: Threading::SingleThread,
            target_fraction: Some(0.9),
        }
    }

    fn elements(&self) -> usize {
        self.x.len()
    }

    fn run_once(&mut self) {
        // Serial chain + division: nowhere near any roof, by construction.
        let mut acc = 1.0f64;
        for &v in &self.x {
            acc = (acc / 1.000_000_01) + v.sqrt().sin();
        }
        self.out = acc;
        std::hint::black_box(self.out);
    }
}

/// The default registry: everything that exists today.
#[must_use]
pub fn default_registry(n: usize) -> Vec<Box<dyn RooflineKernel>> {
    vec![
        Box::new(AxpyKernel::new(n)),
        Box::new(DotKernel::new(n)),
        Box::new(SumKernel::new(n)),
    ]
}

/// The registry used by the shipped `roofline` command. Vector kernels use
/// `n` elements; GEMM uses approximately the same per-matrix element count,
/// with a 256 edge floor so its production parallel/tuning route is real.
#[must_use]
pub fn production_registry(n: usize, axes: &crate::MachineAxes) -> Vec<Box<dyn RooflineKernel>> {
    production_registry_with_ledger(n, axes, None)
}

/// The shipped registry with an optional persistent tune-ledger connection.
/// Supplying one lets a cold process adopt the previous run's validated GEMM
/// row before timing. Ownership keeps fsqlite's deliberately `!Send`
/// connection on this synchronous registry thread.
#[must_use]
pub fn production_registry_with_ledger(
    n: usize,
    axes: &crate::MachineAxes,
    tune_ledger: Option<fs_ledger::Ledger>,
) -> Vec<Box<dyn RooflineKernel>> {
    let mut registry = default_registry(n);
    let side = n.isqrt().max(256);
    registry.push(Box::new(GemmKernel::new(
        side,
        side,
        side,
        axes.logical_cpus as usize,
        axes.fingerprint,
        tune_ledger,
    )));
    registry
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The shipped kernel wrapper must exercise a cold measurement exactly
    /// once, then reuse its in-memory row while continuing to dispatch through
    /// the session API. The narrow N/K shape keeps this control-plane proof
    /// fast while M=256 forces the same tuning route as production.
    #[test]
    fn gemm_kernel_closes_and_reuses_the_tune_loop() {
        let mut kernel = GemmKernel::new(256, 1, 1, 2, 0xC105_ED10, None);
        assert_eq!((kernel.dispatches(), kernel.sweeps()), (0, 0));
        kernel.run_once();
        assert_eq!((kernel.dispatches(), kernel.sweeps()), (1, 1));
        let first_decisions = kernel.tuner.decisions().len();
        assert_eq!(
            first_decisions, 1,
            "first dispatch must record its decision"
        );

        kernel.run_once();
        assert_eq!((kernel.dispatches(), kernel.sweeps()), (2, 1));
        assert_eq!(
            kernel.tuner.decisions().len(),
            first_decisions + 1,
            "warm row must dispatch again without another measurement sweep"
        );
    }

    #[test]
    fn production_registry_contains_session_gemm() {
        let axes = crate::MachineAxes {
            fingerprint: 0xA11_C0DE,
            cpu_brand: "fixture".to_string(),
            logical_cpus: 2,
            bandwidth_single_gbs: 10.0,
            bandwidth_all_core_gbs: 20.0,
            peak_single_gflops: 10.0,
            peak_all_core_gflops: 20.0,
        };
        let registry = production_registry(1, &axes);
        let specs: Vec<_> = registry.iter().map(|kernel| kernel.spec().name).collect();
        assert_eq!(
            specs,
            ["simd-axpy-f64", "simd-dot-f64", "simd-sum-f64", "gemm-f64"]
        );
    }

    /// A new process-local tuner must adopt the validated ledger row instead
    /// of re-measuring. Moving one in-memory ledger connection between kernel
    /// instances isolates the persistent-cache behavior without filesystem
    /// timing or cleanup.
    #[test]
    fn gemm_kernel_adopts_persisted_row_without_resweep() {
        let fingerprint = 0x1ED6_E2ED;
        let ledger = fs_ledger::Ledger::open(":memory:").expect("in-memory tune ledger");
        let mut first = GemmKernel::new(256, 1, 1, 2, fingerprint, Some(ledger));
        first.run_once();
        assert_eq!((first.dispatches(), first.sweeps()), (1, 1));
        let ledger = first.tune_ledger.take().expect("owned ledger");

        let mut replay = GemmKernel::new(256, 1, 1, 2, fingerprint, Some(ledger));
        replay.run_once();
        assert_eq!(
            (replay.dispatches(), replay.sweeps()),
            (1, 0),
            "cold tuner must adopt the persisted row before dispatch"
        );
        assert_eq!(replay.tuner.decisions().len(), 1);
    }
}
