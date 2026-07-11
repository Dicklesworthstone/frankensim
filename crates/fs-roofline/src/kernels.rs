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
//!   `fs_session::gemm_f64_session_with_pool`, so warmup closes measure → cache
//!   → model → dispatch and repetitions reuse the same validated tune row and
//!   TilePool.
//!
//! Targets here are report-only in v0 (`target_fraction: None` except where
//! a band is deliberately claimed for meta-testing): CI runners are shared
//! machines and §14 bands belong to fingerprinted reference machines.

use crate::{KernelExecutionBinding, KernelSpec, RooflineKernel, TargetAxis, Threading};

/// Roofline wrapper version. The lower-layer implementation/tier/placement
/// identities are independently bound inside fs-session's tune key.
pub const GEMM_ROOFLINE_VERSION: &str = "2";

const GEMM_ROOFLINE_POOL_SEED: u64 = 0x524F_4F46_4C49_4E45;

/// Production f64 GEMM benchmark routed through the session autotuner.
///
/// The kernel owns its [`fs_exec::Tuner`], reusable [`fs_exec::TilePool`], and
/// [`fs_exec::CancelGate`]. The
/// roofline registry invokes kernels sequentially through exclusive `&mut`
/// borrows, so tune/cache/dispatch state needs no wrapper lock. A failure is
/// fail-closed: [`RooflineKernel::run_once`] cannot return `Result`, therefore
/// an unexpected session diagnostic panics instead of emitting a normal-
/// looking timing row.
pub struct GemmKernel {
    m: usize,
    n: usize,
    k: usize,
    a: Vec<f64>,
    b: Vec<f64>,
    c: Vec<f64>,
    tuner: fs_exec::Tuner,
    pool: fs_exec::TilePool,
    tune_ledger: Option<fs_ledger::Ledger>,
    pending_tune_row: Option<fs_session::ValidatedGemmTuneRow>,
    active_tune_row: Option<fs_session::ValidatedGemmTuneRow>,
    last_binding: Option<KernelExecutionBinding>,
    gate: fs_exec::CancelGate,
    dispatches: usize,
    sweeps: usize,
    lifecycle_pending: bool,
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
            a,
            b,
            c: vec![0.0; c_len],
            tuner: fs_exec::Tuner::cold(machine_fingerprint),
            pool: fs_exec::TilePool::for_host(threads.max(1), GEMM_ROOFLINE_POOL_SEED),
            tune_ledger,
            pending_tune_row: None,
            active_tune_row: None,
            last_binding: None,
            gate: fs_exec::CancelGate::new(),
            dispatches: 0,
            sweeps: 0,
            lifecycle_pending: false,
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
            target_axis: TargetAxis::ComputePeak,
            target_fraction: Some(0.75),
        }
    }

    fn elements(&self) -> usize {
        self.m * self.n
    }

    fn run_once(&mut self) {
        let cache = self
            .tune_ledger
            .as_ref()
            .map_or(fs_session::GemmTuneCache::Disabled, |ledger| {
                fs_session::GemmTuneCache::ReadOnly(ledger)
            });
        let declared_run = fs_exec::RunId(
            u64::try_from(self.dispatches).expect("roofline GEMM dispatch count exceeds u64"),
        );
        let dispatch = fs_session::gemm_f64_session_with_pool_declared(
            &mut self.tuner,
            cache,
            &self.pool,
            &self.gate,
            declared_run,
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
        self.lifecycle_pending = true;
        self.dispatches += 1;
        self.sweeps += usize::from(dispatch.swept);
        if let Some(row) = dispatch.validated_tune_row.as_ref() {
            self.active_tune_row = Some(row.clone());
        }
        if let Some(row) = dispatch.new_tune_row.as_ref() {
            self.pending_tune_row = Some(row.clone());
        }
        let Some(active_row) = self.active_tune_row.as_ref() else {
            // Serial/small/no-product dispatches have no meaningful MC/NC row.
            // They may still be reported, but GEMM admission requires a sealed
            // decision binding and therefore refuses them as citable evidence.
            self.last_binding = None;
            std::hint::black_box(self.c[self.c.len() / 2]);
            return;
        };
        let expected_key =
            fs_session::gemm_tune::gemm_tune_key_with_pool(&self.pool, self.m, self.n, self.k)
                .unwrap_or_else(|error| panic!("roofline GEMM key construction failed: {error}"));
        assert_eq!(
            dispatch.kernel,
            expected_key.kernel(),
            "session dispatch returned a different scoped key"
        );
        assert_eq!(
            dispatch.shape_class,
            expected_key.shape_class(),
            "session dispatch returned a different shape class"
        );
        let source = match dispatch.source {
            fs_exec::TuneSource::Tuned => "tuned",
            fs_exec::TuneSource::Pinned => "pinned",
            fs_exec::TuneSource::ColdStart => "cold-start",
        };
        self.last_binding = Some(
            KernelExecutionBinding::gemm(
                dispatch.kernel.clone(),
                dispatch.shape_class.clone(),
                dispatch.plan.canonical(),
                source,
                expected_key.execution().isa_tier().to_string(),
                expected_key.execution().build().to_string(),
                active_row.clone(),
                self.tuner.machine(),
                dispatch.execution_receipt(),
            )
            .unwrap_or_else(|error| panic!("roofline GEMM receipt refused: {error}")),
        );
        std::hint::black_box(self.c[self.c.len() / 2]);
    }

    fn execution_binding(&self) -> Option<KernelExecutionBinding> {
        self.last_binding.clone()
    }

    fn pending_tune_publication(&self) -> Option<fs_session::ValidatedGemmTuneRow> {
        self.pending_tune_row.clone()
    }

    fn finalize_tuning(&mut self, admitted: bool) -> Result<(), String> {
        if !self.lifecycle_pending {
            return Err("GEMM registry run has no unfinalized execution state".to_string());
        }
        self.lifecycle_pending = false;
        self.pending_tune_row = None;
        if !admitted {
            self.active_tune_row = None;
            self.last_binding = None;
            self.tuner = fs_exec::Tuner::cold(self.tuner.machine());
        }
        Ok(())
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
            target_axis: TargetAxis::BindingRoof,
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
            target_axis: TargetAxis::BindingRoof,
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
            target_axis: TargetAxis::BindingRoof,
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
            target_axis: TargetAxis::MemoryBandwidth,
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

    fn receipt_axes(fingerprint: u64) -> crate::MachineAxes {
        crate::MachineAxes {
            fingerprint,
            cpu_brand: "receipt-fixture".to_string(),
            logical_cpus: 2,
            bandwidth_single_gbs: 1.0e9,
            bandwidth_all_core_gbs: 2.0e9,
            peak_single_gflops: 1.0e9,
            peak_all_core_gflops: 2.0e9,
        }
    }

    fn trusted_baseline(
        axes: &crate::MachineAxes,
    ) -> (crate::BaselineAxes, crate::BaselineIdentity) {
        let identity = crate::BaselineIdentity::current(axes, "test-firmware")
            .expect("valid synthetic identity");
        let candidates: Vec<_> = (0_u64..3)
            .map(|ordinal| {
                crate::BaselineCandidate::from_receipt(
                    axes.clone(),
                    identity.clone(),
                    fs_blake3::hash_domain(
                        "fs-roofline.test-baseline-source.v1",
                        &ordinal.to_le_bytes(),
                    ),
                )
                .expect("valid synthetic candidate")
            })
            .collect();
        let baseline = crate::promote_baseline(
            &candidates,
            "test-operator",
            "deterministic roofline receipt fixture",
            20_000,
            90,
        )
        .expect("valid synthetic baseline");
        (baseline, identity)
    }

    fn tamper_all_bindings(
        measured: &crate::Attainment,
        mut tamper: impl FnMut(&mut KernelExecutionBinding),
    ) -> crate::Attainment {
        let mut result = measured.clone();
        let crate::MeasurementOrigin::Timed {
            decision_bindings, ..
        } = &mut result.measurement_origin
        else {
            panic!("expected timed receipt");
        };
        for binding in decision_bindings {
            tamper(binding.as_mut().expect("GEMM binding"));
        }
        result
    }

    fn sealed_results_fixture(
        axes: &crate::MachineAxes,
        baseline: crate::AxisBaselinePolicy<'_>,
        results: &[crate::Attainment],
    ) -> crate::FinalizedRegistryRun {
        crate::FinalizedRegistryRun {
            receipt: crate::finalized_run_receipt(axes, axes, baseline, results),
            admitted: crate::run_admission_error(axes, axes, baseline, results).is_none(),
            consumed: false,
        }
    }

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
        kernel
            .finalize_tuning(true)
            .expect("the measured lifecycle finalizes once");
        assert!(
            kernel.finalize_tuning(true).is_err(),
            "a finalized execution state cannot be relabeled by a second registry run"
        );
    }

    #[test]
    fn rejected_gemm_lifecycle_invalidates_the_local_decision() {
        let mut kernel = GemmKernel::new(256, 1, 1, 2, 0xBAD0_CAFE, None);
        kernel.run_once();
        assert!(kernel.active_tune_row.is_some());
        assert!(kernel.last_binding.is_some());
        let sweeps_before_rejection = kernel.sweeps();

        kernel
            .finalize_tuning(false)
            .expect("rejected lifecycle drains once");
        assert!(kernel.active_tune_row.is_none());
        assert!(kernel.pending_tune_row.is_none());
        assert!(kernel.last_binding.is_none());

        kernel.run_once();
        assert_eq!(
            kernel.sweeps(),
            sweeps_before_rejection + 1,
            "reuse after rejection must revalidate through a fresh sweep"
        );
    }

    #[test]
    fn stale_results_cannot_finalize_a_newer_gemm_execution_state() {
        let axes = receipt_axes(0x51A1_E001);
        let (baseline, identity) = trusted_baseline(&axes);
        let baseline_policy = crate::AxisBaselinePolicy::new(Some(&baseline), &identity, 20_010);
        let mut registry: Vec<Box<dyn RooflineKernel>> = vec![Box::new(GemmKernel::new(
            256,
            1,
            1,
            2,
            axes.fingerprint,
            None,
        ))];
        let first = crate::run_registry(&mut registry, 1, 1, &axes);
        let second = crate::run_registry(&mut registry, 1, 1, &axes);
        assert_ne!(
            first[0].to_jsonl(),
            second[0].to_jsonl(),
            "declared execution runs distinguish the two result sets"
        );

        let error =
            crate::finalize_registry_tuning(&mut registry, &axes, &axes, baseline_policy, &first)
                .expect_err("an old result must not govern newer pending kernel state");
        assert!(
            error.contains("execution state changed after this result was measured"),
            "{error}"
        );
        assert!(
            registry[0].pending_tune_publication().is_none(),
            "the newer run's pending tune publication must be drained on mismatch"
        );
        assert!(
            registry[0].execution_binding().is_none(),
            "the newer local decision must be invalidated rather than surviving the refusal"
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

    /// An admitted cold run publishes its buffered row in the same transaction
    /// as its citable evidence; a new process-local tuner must then adopt that
    /// validated ledger row instead of re-measuring.
    /// Moving one in-memory ledger connection between kernel instances isolates
    /// the persistent-cache behavior without filesystem timing or cleanup.
    #[test]
    #[allow(clippy::too_many_lines)] // one end-to-end two-ledger provenance scenario
    fn gemm_kernel_adopts_persisted_row_without_resweep() {
        let fingerprint = 0x1ED6_E2ED;
        let axes = receipt_axes(fingerprint);
        let (baseline, identity) = trusted_baseline(&axes);
        let baseline_policy = crate::AxisBaselinePolicy::new(Some(&baseline), &identity, 20_010);
        let ledger = fs_ledger::Ledger::open(":memory:").expect("in-memory tune ledger");
        let mut first = GemmKernel::new(256, 1, 1, 2, fingerprint, Some(ledger));
        let mut measured = crate::measure(&mut first, 1, 2, &axes);
        let mut stale_fresh_clone = measured.clone();
        assert!(
            stale_fresh_clone.pending_tune_publication.is_some(),
            "regression fixture must duplicate the former overwrite marker"
        );
        assert_eq!((first.dispatches(), first.sweeps()), (3, 1));
        assert!(crate::run_is_citable(
            &axes,
            &axes,
            baseline_policy,
            std::slice::from_ref(&measured)
        ));
        assert_eq!(
            first
                .tune_ledger
                .as_ref()
                .expect("ledger")
                .table_count("tune")
                .expect("count"),
            0,
            "measurement remains process-local before admission"
        );
        first
            .finalize_tuning(true)
            .expect("admitted lifecycle finalizes");
        let mut finalized_first =
            sealed_results_fixture(&axes, baseline_policy, std::slice::from_ref(&measured));
        assert_eq!(
            first
                .tune_ledger
                .as_ref()
                .expect("ledger")
                .table_count("tune")
                .expect("count"),
            0,
            "finalization cannot publish outside the evidence transaction"
        );
        crate::record_run(
            first.tune_ledger.as_ref().expect("ledger"),
            &axes,
            &axes,
            baseline_policy,
            &mut finalized_first,
            std::slice::from_mut(&mut measured),
        )
        .expect("atomic evidence transaction");
        assert_eq!(
            first
                .tune_ledger
                .as_ref()
                .expect("ledger")
                .table_count("tune")
                .expect("count"),
            2,
            "one session tune row and one citable roofline row commit together"
        );
        assert!(first.pending_tune_row.is_none());
        assert!(
            measured.pending_tune_publication.is_none(),
            "successful commit consumes the one-shot fresh-row marker"
        );
        let ledger = first.tune_ledger.take().expect("owned ledger");

        let mut replay = GemmKernel::new(256, 1, 1, 2, fingerprint, Some(ledger));
        let mut adopted = crate::measure(&mut replay, 1, 2, &axes);
        assert_eq!(
            (replay.dispatches(), replay.sweeps()),
            (3, 0),
            "cold tuner must adopt the persisted row before dispatch"
        );
        assert_eq!(replay.tuner.decisions().len(), 3);
        assert!(replay.active_tune_row.is_some());
        assert!(replay.pending_tune_row.is_none());
        assert!(crate::run_is_citable(
            &axes,
            &axes,
            baseline_policy,
            std::slice::from_ref(&adopted)
        ));
        let original_binding = crate::stable_decision_binding(match &measured.measurement_origin {
            crate::MeasurementOrigin::Timed {
                decision_bindings, ..
            } => decision_bindings,
            crate::MeasurementOrigin::Analytic => panic!("expected timed receipt"),
        })
        .expect("original binding")
        .clone();
        let adopted_binding = crate::stable_decision_binding(match &adopted.measurement_origin {
            crate::MeasurementOrigin::Timed {
                decision_bindings, ..
            } => decision_bindings,
            crate::MeasurementOrigin::Analytic => panic!("expected timed receipt"),
        })
        .expect("adopted binding")
        .clone();
        assert_eq!(
            adopted_binding.gemm.tune_row_identity, original_binding.gemm.tune_row_identity,
            "adoption must bind the exact measured row identity"
        );
        let sealed_row_json = adopted_binding.validated_row().receipt_json();
        assert_eq!(
            fs_blake3::hash_domain(
                fs_session::GEMM_TUNE_ROW_RECEIPT_DOMAIN,
                sealed_row_json.as_bytes(),
            ),
            adopted_binding.gemm.tune_row_identity,
            "the embedded canonical row must be the exact identity preimage"
        );
        let historical_payload = adopted.to_jsonl();
        assert!(
            historical_payload.contains(&format!("\"tune_row\":{sealed_row_json}")),
            "the benchmark payload must embed its immutable tune-row preimage"
        );

        let evidence_ledger = fs_ledger::Ledger::open(":memory:").expect("evidence ledger B");
        replay
            .finalize_tuning(true)
            .expect("adopted lifecycle finalizes");
        let mut finalized_adopted =
            sealed_results_fixture(&axes, baseline_policy, std::slice::from_ref(&adopted));
        crate::record_run(
            &evidence_ledger,
            &axes,
            &axes,
            baseline_policy,
            &mut finalized_adopted,
            std::slice::from_mut(&mut adopted),
        )
        .expect("record adopted receipt into independent ledger");
        assert_eq!(
            evidence_ledger.table_count("tune").expect("count"),
            2,
            "independent evidence ledger receives the bound session and roofline rows"
        );
        let stored = evidence_ledger
            .tune_get(
                &adopted_binding.gemm.scoped_tune_key,
                &adopted_binding.gemm.shape_class,
                &fingerprint.to_le_bytes(),
            )
            .expect("query exact session row")
            .expect("session row committed with receipt");
        assert!(
            adopted_binding.validated_row().matches_ledger_row(&stored),
            "ledger B must contain the exact sealed row, not only its hash"
        );

        let mut newer = GemmKernel::new(256, 1, 1, 2, fingerprint, None);
        newer.run_once();
        let newer_row = newer
            .pending_tune_row
            .as_ref()
            .expect("independently validated newer row");
        assert_ne!(
            newer_row.receipt_identity(),
            adopted_binding.gemm.tune_row_identity,
            "independent measurement must carry distinct timing evidence"
        );
        let conflict_ledger =
            fs_ledger::Ledger::open(":memory:").expect("conflict evidence ledger");
        newer_row
            .publish_to_ledger(&conflict_ledger)
            .expect("seed newer destination row");
        let stale_fresh_conflict = crate::record_run(
            &conflict_ledger,
            &axes,
            &axes,
            baseline_policy,
            &mut sealed_results_fixture(
                &axes,
                baseline_policy,
                std::slice::from_ref(&stale_fresh_clone),
            ),
            std::slice::from_mut(&mut stale_fresh_clone),
        )
        .expect_err("a delayed clone of fresh evidence must not replace a newer row");
        assert!(
            stale_fresh_conflict
                .to_string()
                .contains("conflicting tune row")
        );
        assert!(
            stale_fresh_clone.pending_tune_publication.is_some(),
            "transaction rollback retains a retry marker without granting overwrite authority"
        );
        let conflict = crate::record_run(
            &conflict_ledger,
            &axes,
            &axes,
            baseline_policy,
            &mut sealed_results_fixture(&axes, baseline_policy, std::slice::from_ref(&adopted)),
            std::slice::from_mut(&mut adopted),
        )
        .expect_err("adopted evidence must not replace a conflicting destination row");
        assert!(conflict.to_string().contains("conflicting tune row"));
        let retained = conflict_ledger
            .tune_get(
                &adopted_binding.gemm.scoped_tune_key,
                &adopted_binding.gemm.shape_class,
                &fingerprint.to_le_bytes(),
            )
            .expect("query retained newer row")
            .expect("newer row survives conflict");
        assert!(
            newer_row.matches_ledger_row(&retained),
            "failed adopted receipt must leave the destination row unchanged"
        );
        assert!(
            historical_payload.contains(&sealed_row_json)
                && !historical_payload.contains(&newer_row.receipt_json()),
            "later cache replacement must not change the historical receipt preimage"
        );
    }

    #[test]
    fn rejected_gemm_run_discards_buffer_without_persistent_row() {
        let ledger = fs_ledger::Ledger::open(":memory:").expect("in-memory tune ledger");
        let mut kernel = GemmKernel::new(256, 1, 1, 2, 0xBAD_A11E5, Some(ledger));
        kernel.run_once();
        assert_eq!((kernel.dispatches(), kernel.sweeps()), (1, 1));
        assert!(kernel.pending_tune_row.is_some());

        kernel
            .finalize_tuning(false)
            .expect("rejection only discards local state");
        assert!(kernel.pending_tune_row.is_none());
        assert!(
            !kernel.tuner.has_gemm_row(
                &fs_session::gemm_tune::gemm_tune_key(2, 256, 1, 1).expect("tune key"),
            ),
            "a reusable registry must not retain rejected process-local tuning"
        );
        assert_eq!(
            kernel
                .tune_ledger
                .as_ref()
                .expect("ledger")
                .table_count("tune")
                .expect("count"),
            0,
            "rejected work must never contaminate the durable cache"
        );
    }

    #[test]
    fn admitted_finalize_consumes_pending_marker_before_registry_reuse() {
        let fingerprint = 0xA11D_771D;
        let mut kernel = GemmKernel::new(256, 1, 1, 2, fingerprint, None);
        kernel.run_once();
        assert!(kernel.pending_tune_row.is_some());
        kernel
            .finalize_tuning(true)
            .expect("admitted lifecycle finalizes");
        assert!(kernel.pending_tune_row.is_none());
        kernel.run_once();
        assert!(kernel.pending_tune_row.is_none());
        assert!(kernel.active_tune_row.is_some());
    }

    #[test]
    fn untuned_serial_gemm_reports_but_cannot_be_cited() {
        let axes = receipt_axes(0x51A1);
        let (baseline, identity) = trusted_baseline(&axes);
        let baseline_policy = crate::AxisBaselinePolicy::new(Some(&baseline), &identity, 20_010);
        let mut kernel = GemmKernel::square(8, 1, axes.fingerprint);
        let measured = crate::measure(&mut kernel, 1, 1, &axes);
        assert_eq!((kernel.dispatches(), kernel.sweeps()), (2, 0));
        assert!(kernel.last_binding.is_none());
        assert!(!crate::run_is_citable(
            &axes,
            &axes,
            baseline_policy,
            &[measured]
        ));
    }

    #[test]
    fn citable_gemm_receipt_rejects_every_bound_field_tamper() {
        let axes = receipt_axes(0x55_00_11);
        let (baseline, identity) = trusted_baseline(&axes);
        let baseline_policy = crate::AxisBaselinePolicy::new(Some(&baseline), &identity, 20_010);
        let mut kernel = GemmKernel::new(256, 1, 1, 2, axes.fingerprint, None);
        let measured = crate::measure(&mut kernel, 1, 3, &axes);
        assert!(crate::run_is_citable(
            &axes,
            &axes,
            baseline_policy,
            std::slice::from_ref(&measured)
        ));
        let json = measured.to_jsonl();
        for field in [
            "scoped_tune_key",
            "shape_class",
            "plan",
            "source",
            "operation_tier",
            "build_identity",
            "tune_row_identity",
            "execution_path_identity",
            "execution_path",
            "declared_run",
            "decision_binding_hashes",
            "warmup_runs",
        ] {
            assert!(json.contains(field), "receipt omitted {field}");
        }

        let rejected = [
            tamper_all_bindings(&measured, |binding| binding.gemm.scoped_tune_key.push('x')),
            tamper_all_bindings(&measured, |binding| binding.gemm.shape_class.push('x')),
            tamper_all_bindings(&measured, |binding| binding.gemm.canonical_plan.push('x')),
            tamper_all_bindings(&measured, |binding| binding.gemm.source = "pinned"),
            tamper_all_bindings(&measured, |binding| binding.gemm.operation_tier.push('x')),
            tamper_all_bindings(&measured, |binding| binding.gemm.build_identity.push('x')),
            tamper_all_bindings(&measured, |binding| {
                binding.gemm.tune_row_identity = fs_ledger::hash_bytes(b"tampered");
            }),
            tamper_all_bindings(&measured, |binding| {
                binding.gemm.execution_path.completed_tiles = 0;
            }),
            tamper_all_bindings(&measured, |binding| {
                binding.gemm.execution_path.panels[0].mode.push('x');
            }),
            tamper_all_bindings(&measured, |binding| {
                binding.gemm.execution_path.panels[0].declared_run = 1;
            }),
            tamper_all_bindings(&measured, |binding| {
                binding.gemm.execution_path_identity = fs_ledger::hash_bytes(b"tampered");
            }),
        ];
        for (index, result) in rejected.into_iter().enumerate() {
            assert!(
                !crate::run_is_citable(&axes, &axes, baseline_policy, &[result]),
                "tamper case {index} was admitted"
            );
        }

        let mut unstable = measured.clone();
        let crate::MeasurementOrigin::Timed {
            decision_bindings, ..
        } = &mut unstable.measurement_origin
        else {
            panic!("expected timed receipt");
        };
        decision_bindings[0]
            .as_mut()
            .expect("binding")
            .gemm
            .canonical_plan
            .push('x');
        assert!(!crate::run_is_citable(
            &axes,
            &axes,
            baseline_policy,
            &[unstable]
        ));

        let mut no_warmup = measured;
        let crate::MeasurementOrigin::Timed { warmup_runs, .. } = &mut no_warmup.measurement_origin
        else {
            panic!("expected timed receipt");
        };
        *warmup_runs = 0;
        assert!(!crate::run_is_citable(
            &axes,
            &axes,
            baseline_policy,
            &[no_warmup]
        ));
    }
}
