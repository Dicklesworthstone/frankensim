//! Sealed production-run protocol (bead fz2.5): the only path to citable
//! roofline evidence.
//!
//! The public `RooflineKernel`/`run_registry`/`record_run` surface treats
//! caller-supplied kernel implementations and `MachineAxes` as a trust root:
//! a caller can clone a valid GEMM `KernelExecutionBinding` into a fake
//! custom kernel named `gemm-f64`, or discard a drifted post-probe and pass
//! the pre-probe twice. That surface stays available for harness tests and
//! exploration, but everything it records is stamped
//! `"protocol":"custom-registry"` and is explicitly NON-CITABLE.
//!
//! The protocol is two opaque stages:
//!
//! 1. [`ProductionProbe::observe`] performs the pre-run axis probe and mints
//!    the per-run nonce. The caller may READ the observed axes (baseline
//!    selection needs them) but can never supply its own.
//! 2. [`ProductionProbe::run`] owns production registry selection, timed
//!    warmup/repetitions, the post-run axis probe (observed strictly after
//!    the timed loop), aggregate admission, and tune finalization, yielding
//!    a [`ProductionRun`]. [`ProductionRun::record`] commits atomically and
//!    consumes the run; the operation `ir` carries
//!    `"protocol":"production-v2"`, the nonce, content hashes of both
//!    observed axis receipts, and the retained dependency-receipt binding.
//!
//! Trust model: the nonce is a process-unique challenge, not cryptographic
//! proof. Type opacity prevents ordinary API consumers from constructing a
//! `ProductionRun`, but `fs-ledger` intentionally exposes general mutation
//! APIs. A trusted ledger writer can therefore mint or replace internally
//! consistent rows. External authentication of the ledger/package is a
//! separate proof obligation; this crate detects corruption inside that
//! trusted-writer boundary and makes no cryptographic-authority claim.

use fs_ledger::{Ledger, LedgerError};

use crate::kernels::production_registry_with_ledger;
use crate::{
    Attainment, AxisBaselinePolicy, DependencyReceiptBinding, FinalizedRegistryRun, MachineAxes,
    PRODUCTION_AXES_RECEIPT_DOMAIN, PRODUCTION_PROTOCOL_FIELD, RooflineKernel,
    finalize_registry_tuning, json_escape, record_run_with_protocol, run_admission_error,
    run_registry,
};

const RUN_NONCE_DOMAIN: &str = "org.frankensim.fs-roofline.production-run-nonce.v1";
#[cfg(test)]
const DEVELOPMENT_SALT_REFUSAL: &str = "dependency graph uses the development equivalence salt; production citation requires an exact operator-observed normal/build receipt";

static NONCE_COUNTER: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CitationAuthority {
    Receipt(DependencyReceiptBinding),
    Refused(&'static str),
}

impl CitationAuthority {
    fn from_build() -> Self {
        DependencyReceiptBinding::current().map_or_else(Self::Refused, Self::Receipt)
    }

    const fn refusal(self) -> Option<&'static str> {
        match self {
            Self::Receipt(_) => None,
            Self::Refused(reason) => Some(reason),
        }
    }

    const fn receipt(self) -> Option<DependencyReceiptBinding> {
        match self {
            Self::Receipt(binding) => Some(binding),
            Self::Refused(_) => None,
        }
    }
}

/// Sizing and repetition parameters for one production run.
#[derive(Debug, Clone, Copy)]
pub struct ProductionRunConfig {
    /// Vector-kernel element count (GEMM derives its edge from this).
    pub n: usize,
    /// Untimed warmup repetitions per kernel.
    pub warmup: usize,
    /// Timed repetitions per kernel.
    pub reps: usize,
}

/// Largest vector length accepted by the sealed production runner.
///
/// The shipped registry holds several `f64` vectors plus three approximately
/// `n`-element GEMM matrices, so this bounds aggregate allocation as well as
/// each individual kernel input.
pub const MAX_PRODUCTION_ELEMENTS: usize = crate::kernels::MAX_VECTOR_KERNEL_ELEMENTS;
/// Largest untimed warmup count accepted by the sealed production runner.
pub const MAX_PRODUCTION_WARMUP: usize = crate::MAX_MEASUREMENT_WARMUP;
/// Largest timed repetition count accepted by the sealed production runner.
pub const MAX_PRODUCTION_REPS: usize = crate::MAX_MEASUREMENT_REPS;
/// Largest `n * (warmup + reps)` admitted by one production invocation.
pub const MAX_PRODUCTION_ELEMENT_PASSES: usize = 1 << 28;

impl ProductionRunConfig {
    /// Validate resource-driving inputs before registry allocation or timing.
    ///
    /// # Errors
    /// Returns a stable diagnostic for zero or out-of-envelope inputs.
    pub fn validate(self) -> Result<(), String> {
        if self.n == 0 || self.n > MAX_PRODUCTION_ELEMENTS {
            return Err(format!(
                "production n must be in 1..={MAX_PRODUCTION_ELEMENTS}, got {}",
                self.n
            ));
        }
        if self.warmup > MAX_PRODUCTION_WARMUP {
            return Err(format!(
                "production warmup must be in 0..={MAX_PRODUCTION_WARMUP}, got {}",
                self.warmup
            ));
        }
        if self.reps == 0 || self.reps > MAX_PRODUCTION_REPS {
            return Err(format!(
                "production reps must be in 1..={MAX_PRODUCTION_REPS}, got {}",
                self.reps
            ));
        }
        let passes = self
            .warmup
            .checked_add(self.reps)
            .and_then(|count| self.n.checked_mul(count))
            .ok_or_else(|| "production element-pass budget overflowed usize".to_string())?;
        if passes > MAX_PRODUCTION_ELEMENT_PASSES {
            return Err(format!(
                "production n * (warmup + reps) must be at most {MAX_PRODUCTION_ELEMENT_PASSES}, got {passes}"
            ));
        }
        Ok(())
    }
}

/// Stage one of the sealed protocol: a pre-run axis probe this crate
/// performed itself, plus the minted per-run nonce.
///
/// No public constructor accepts axes; the only public way in is
/// [`ProductionProbe::observe`], which probes the actual machine.
pub struct ProductionProbe {
    axes: MachineAxes,
    nonce: fs_blake3::ContentHash,
}

impl ProductionProbe {
    /// Probe the machine and mint this run's nonce.
    #[must_use]
    pub fn observe() -> Self {
        Self::from_observed(MachineAxes::probe())
    }

    /// Test seam (`pub(crate)`): inject a synthetic pre-probe. Unreachable
    /// by API consumers, so a forged probe cannot enter the protocol.
    pub(crate) fn from_observed(axes: MachineAxes) -> Self {
        let nonce = mint_nonce(&axes);
        Self { axes, nonce }
    }

    /// The observed pre-run axes (read-only; baseline selection needs them).
    #[must_use]
    pub fn axes(&self) -> &MachineAxes {
        &self.axes
    }

    /// Run the production registry and finalize, consuming the probe.
    ///
    /// The tune ledger (optional) lets the GEMM kernel adopt a previously
    /// validated row; the registry (and with it fsqlite's `!Send`
    /// connection) is dropped before this returns, so the caller may reopen
    /// the same database for [`ProductionRun::record`].
    ///
    /// # Errors
    /// Structured diagnostics from tuning finalization; admission refusal is
    /// NOT an error — the run comes back with `citation_eligible() == false`
    /// and can be recorded as an explicit rejection.
    pub fn run(
        self,
        config: ProductionRunConfig,
        baseline: AxisBaselinePolicy<'_>,
        tune_ledger: Option<Ledger>,
    ) -> Result<ProductionRun<'_>, String> {
        config.validate()?;
        let registry = production_registry_with_ledger(config.n, &self.axes, tune_ledger)?;
        self.run_with_parts(config, baseline, registry, MachineAxes::probe)
    }

    /// Protocol core with injected registry and post-probe (`pub(crate)`
    /// test seam: drifted-post and finalizer-failure paths need determinism;
    /// API consumers cannot reach this to forge a run).
    pub(crate) fn run_with_parts(
        self,
        config: ProductionRunConfig,
        baseline: AxisBaselinePolicy<'_>,
        registry: Vec<Box<dyn RooflineKernel>>,
        post_probe: impl FnOnce() -> MachineAxes,
    ) -> Result<ProductionRun<'_>, String> {
        self.run_with_parts_and_authority(
            config,
            baseline,
            registry,
            post_probe,
            CitationAuthority::from_build(),
        )
    }

    fn run_with_parts_and_authority(
        self,
        config: ProductionRunConfig,
        baseline: AxisBaselinePolicy<'_>,
        mut registry: Vec<Box<dyn RooflineKernel>>,
        post_probe: impl FnOnce() -> MachineAxes,
        citation_authority: CitationAuthority,
    ) -> Result<ProductionRun<'_>, String> {
        config.validate()?;
        let build_identity = crate::read_executable_build_identity().map_err(|error| {
            format!("cannot capture pre-measurement executable identity: {error}")
        })?;
        let results = run_registry(&mut registry, config.warmup, config.reps, &self.axes)?;
        let post_axes = post_probe();
        let finalized =
            finalize_registry_tuning(&mut registry, &self.axes, &post_axes, baseline, &results)?;
        drop(registry);
        Ok(ProductionRun {
            axes: self.axes,
            post_axes,
            baseline,
            nonce: self.nonce,
            results,
            finalized,
            citation_authority,
            build_identity,
        })
    }

    #[cfg(test)]
    pub(crate) fn run_with_test_receipt<'a>(
        self,
        config: ProductionRunConfig,
        baseline: AxisBaselinePolicy<'a>,
        registry: Vec<Box<dyn RooflineKernel>>,
        post_probe: impl FnOnce() -> MachineAxes,
        receipt: &'static str,
    ) -> Result<ProductionRun<'a>, String> {
        let digest =
            fs_blake3::hash_domain(fs_session::GEMM_DEPGRAPH_RECEIPT_DOMAIN, receipt.as_bytes());
        let binding = DependencyReceiptBinding::from_parts(receipt, digest)
            .expect("test receipt digest was computed from the same bytes");
        self.run_with_parts_and_authority(
            config,
            baseline,
            registry,
            post_probe,
            CitationAuthority::Receipt(binding),
        )
    }
}

/// One complete, sealed production registry run.
///
/// Fields are private, there is no public constructor, and the value is
/// neither `Clone` nor `Copy`: the only way to obtain one is
/// [`ProductionProbe::run`], which performed both probes and timed the
/// production registry itself. [`ProductionRun::record`] consumes the run.
pub struct ProductionRun<'a> {
    axes: MachineAxes,
    post_axes: MachineAxes,
    baseline: AxisBaselinePolicy<'a>,
    nonce: fs_blake3::ContentHash,
    results: Vec<Attainment>,
    finalized: FinalizedRegistryRun,
    citation_authority: CitationAuthority,
    build_identity: fs_blake3::ContentHash,
}

impl std::fmt::Debug for ProductionRun<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("ProductionRun")
            .field(
                "fingerprint",
                &format_args!("{:016x}", self.axes.fingerprint),
            )
            .field("kernels", &self.results.len())
            .field("nonce", &self.nonce)
            .field("citation_eligible", &self.citation_eligible())
            .finish_non_exhaustive()
    }
}

impl ProductionRun<'_> {
    /// The pre-run axis probe observed by the protocol.
    #[must_use]
    pub fn axes(&self) -> &MachineAxes {
        &self.axes
    }

    /// The post-run axis probe, observed strictly after the timed loop.
    #[must_use]
    pub fn post_axes(&self) -> &MachineAxes {
        &self.post_axes
    }

    /// The measured result set in registry order.
    #[must_use]
    pub fn results(&self) -> &[Attainment] {
        &self.results
    }

    /// The per-run nonce bound into the recorded operation.
    #[must_use]
    pub fn nonce(&self) -> fs_blake3::ContentHash {
        self.nonce
    }

    /// Whether this run passed aggregate admission and has production receipt
    /// provenance. This is a pre-commit eligibility predicate, not a claim
    /// that evidence was durably recorded or remains fresh.
    #[must_use]
    pub fn citation_eligible(&self) -> bool {
        self.finalized.admitted() && self.citation_authority.refusal().is_none()
    }

    /// Why admission refused this run, if it did.
    #[must_use]
    pub fn admission_error(&self) -> Option<String> {
        run_admission_error(&self.axes, &self.post_axes, self.baseline, &self.results)
            .or_else(|| self.citation_authority.refusal().map(str::to_string))
    }

    /// The baseline-admission receipt for this run's exact probe pair.
    #[must_use]
    pub fn receipt_json(&self) -> String {
        self.baseline.receipt_json(&self.axes, &self.post_axes)
    }

    /// Record the run atomically, consuming it. The operation `ir` carries
    /// `"protocol":"production-v2"`, the per-run nonce, content hashes of
    /// both observed axis receipts, and dependency-receipt provenance.
    ///
    /// # Errors
    /// Ledger errors propagate and roll back the whole write set; the run is
    /// consumed either way (a failed transaction cannot be replayed into a
    /// different ledger with edited results).
    pub fn record(mut self, ledger: &Ledger) -> Result<i64, LedgerError> {
        let protocol_fields = match self.citation_authority {
            CitationAuthority::Receipt(binding) => format!(
                "{PRODUCTION_PROTOCOL_FIELD},\"run_nonce\":\"{}\",\"pre_axes_receipt\":\"{}\",\"post_axes_receipt\":\"{}\",\"dependency_graph_evidence\":\"operator-observed-receipt\",\"dependency_receipt_digest\":\"{}\",\"dependency_receipt_artifact\":\"{}\"",
                self.nonce,
                axes_receipt(&self.axes),
                axes_receipt(&self.post_axes),
                binding.domain_digest,
                binding.artifact_hash,
            ),
            CitationAuthority::Refused(reason) => format!(
                "{PRODUCTION_PROTOCOL_FIELD},\"run_nonce\":\"{}\",\"pre_axes_receipt\":\"{}\",\"post_axes_receipt\":\"{}\",\"dependency_graph_evidence\":\"development-equivalence-salt\",\"citation_refusal\":\"{}\"",
                self.nonce,
                axes_receipt(&self.axes),
                axes_receipt(&self.post_axes),
                json_escape(reason),
            ),
        };
        record_run_with_protocol(
            ledger,
            &self.axes,
            &self.post_axes,
            self.baseline,
            &mut self.finalized,
            &mut self.results,
            &protocol_fields,
            self.citation_authority.receipt(),
            self.citation_authority.refusal(),
            crate::EvidenceNamespace::Production,
            Some(self.build_identity),
        )
    }
}

/// Content hash of one observed probe's canonical JSONL receipt.
fn axes_receipt(axes: &MachineAxes) -> fs_blake3::ContentHash {
    fs_blake3::hash_domain(
        PRODUCTION_AXES_RECEIPT_DOMAIN,
        crate::machine_axes_receipt_json(axes).as_bytes(),
    )
}

/// Process-unique per-run challenge: wall clock, pid, a monotone counter,
/// and the pre-probe receipt. Uniqueness, not secrecy — see the module docs.
fn mint_nonce(axes: &MachineAxes) -> fs_blake3::ContentHash {
    let count = NONCE_COUNTER.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
    let mut material = Vec::new();
    material.extend_from_slice(&fs_ledger::now_wall_ns().to_le_bytes());
    material.extend_from_slice(&u64::from(std::process::id()).to_le_bytes());
    material.extend_from_slice(&count.to_le_bytes());
    material.extend_from_slice(axes.to_jsonl().as_bytes());
    fs_blake3::hash_domain(RUN_NONCE_DOMAIN, &material)
}

#[cfg(test)]
mod tests {
    use std::cell::Cell;
    use std::rc::Rc;

    use super::*;
    use crate::kernels::default_registry;
    use crate::{
        BaselineAxes, BaselineCandidate, BaselineIdentity, KernelSpec, STALENESS_MAX_AGE_NS,
        Staleness, TargetAxis, Threading, promote_baseline, roofline_machine_key, staleness,
        staleness_at,
    };

    fn synthetic_axes(fingerprint: u64) -> MachineAxes {
        // Roofs far above any real machine (bead xjhz): cache-resident test
        // kernels must never outrun the fixture roof.
        MachineAxes {
            fingerprint,
            cpu_brand: "synthetic".to_string(),
            logical_cpus: 8,
            bandwidth_single_gbs: 100_000.0,
            bandwidth_all_core_gbs: 400_000.0,
            peak_single_gflops: 50_000.0,
            peak_all_core_gflops: 300_000.0,
        }
    }

    fn trusted_baseline(axes: &MachineAxes) -> (BaselineAxes, BaselineIdentity) {
        let identity =
            BaselineIdentity::current(axes, "test-firmware").expect("valid synthetic identity");
        let candidates: Vec<_> = (0_u64..3)
            .map(|ordinal| {
                BaselineCandidate::from_receipt(
                    axes.clone(),
                    identity.clone(),
                    fs_blake3::hash_domain(
                        "fs-roofline.production-baseline-source.v1",
                        &ordinal.to_le_bytes(),
                    ),
                )
                .expect("valid synthetic candidate")
            })
            .collect();
        let baseline = promote_baseline(
            &candidates,
            "test-operator",
            "deterministic production-protocol fixture",
            20_000,
            90,
        )
        .expect("valid synthetic baseline");
        (baseline, identity)
    }

    fn temp_db(tag: &str) -> String {
        static NEXT: std::sync::atomic::AtomicU32 = std::sync::atomic::AtomicU32::new(0);
        let n = NEXT.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        std::env::temp_dir()
            .join(format!(
                "fs-roofline-prod-{tag}-{}-{n}.db",
                std::process::id()
            ))
            .display()
            .to_string()
    }

    const CONFIG: ProductionRunConfig = ProductionRunConfig {
        n: 1 << 10,
        warmup: 0,
        reps: 1,
    };

    const TEST_DEPGRAPH_RECEIPT: &str = "{\"schema\":\"fs-roofline-synthetic-dependency-receipt-v1\",\"purpose\":\"unit-test-only\"}";
    const SUBSTITUTE_DEPGRAPH_RECEIPT: &str = "{\"schema\":\"fs-roofline-synthetic-dependency-receipt-v1\",\"purpose\":\"substitution-attack\"}";

    fn test_receipt_binding() -> DependencyReceiptBinding {
        let digest = fs_blake3::hash_domain(
            fs_session::GEMM_DEPGRAPH_RECEIPT_DOMAIN,
            TEST_DEPGRAPH_RECEIPT.as_bytes(),
        );
        DependencyReceiptBinding::from_parts(TEST_DEPGRAPH_RECEIPT, digest)
            .expect("test receipt digest agrees")
    }

    fn test_receipt_authority() -> CitationAuthority {
        CitationAuthority::Receipt(test_receipt_binding())
    }

    #[allow(clippy::too_many_arguments)]
    fn receipt_staleness_at(
        ledger: &Ledger,
        kernel: &str,
        version: &str,
        fingerprint: u64,
        baseline: Option<fs_blake3::ContentHash>,
        observed_wall_ns: i64,
        dependency: DependencyReceiptBinding,
    ) -> Result<Staleness, fs_ledger::LedgerError> {
        crate::staleness_at_with_build_and_dependency(
            ledger,
            kernel,
            version,
            fingerprint,
            baseline,
            observed_wall_ns,
            crate::executable_build_identity()?,
            Some(dependency),
        )
    }

    struct CountingKernel {
        runs: Rc<Cell<usize>>,
        value: u64,
    }

    impl crate::RooflineKernel for CountingKernel {
        fn spec(&self) -> KernelSpec {
            KernelSpec {
                name: "counting-kernel",
                version: "1",
                bytes_per_elem: 8.0,
                flops_per_elem: 1.0,
                threading: Threading::SingleThread,
                target_axis: TargetAxis::BindingRoof,
                target_fraction: None,
            }
        }

        fn elements(&self) -> usize {
            64
        }

        fn run_once(&mut self) -> Result<(), String> {
            self.runs.set(self.runs.get() + 1);
            for _ in 0..64 {
                self.value = std::hint::black_box(
                    self.value
                        .wrapping_mul(6_364_136_223_846_793_005)
                        .wrapping_add(1),
                );
            }
            Ok(())
        }
    }

    struct FailingFinalizeKernel {
        inner: CountingKernel,
    }

    impl crate::RooflineKernel for FailingFinalizeKernel {
        fn spec(&self) -> KernelSpec {
            self.inner.spec()
        }
        fn elements(&self) -> usize {
            self.inner.elements()
        }
        fn run_once(&mut self) -> Result<(), String> {
            self.inner.run_once()
        }
        fn finalize_tuning(&mut self, _admitted: bool) -> Result<(), String> {
            Err("tune ledger unavailable mid-finalize".to_string())
        }
    }

    #[test]
    fn nonces_are_unique_per_probe() {
        let a = ProductionProbe::from_observed(synthetic_axes(0xA));
        let b = ProductionProbe::from_observed(synthetic_axes(0xA));
        assert_ne!(
            a.nonce, b.nonce,
            "identical axes must still mint distinct nonces"
        );
    }

    #[test]
    fn production_config_rejects_zero_and_unbounded_work_before_running() {
        for (config, expected) in [
            (
                ProductionRunConfig {
                    n: 0,
                    warmup: 0,
                    reps: 1,
                },
                "production n",
            ),
            (
                ProductionRunConfig {
                    n: MAX_PRODUCTION_ELEMENTS + 1,
                    warmup: 0,
                    reps: 1,
                },
                "production n",
            ),
            (
                ProductionRunConfig {
                    n: 1,
                    warmup: MAX_PRODUCTION_WARMUP + 1,
                    reps: 1,
                },
                "production warmup",
            ),
            (
                ProductionRunConfig {
                    n: 1,
                    warmup: 0,
                    reps: 0,
                },
                "production reps",
            ),
            (
                ProductionRunConfig {
                    n: 1,
                    warmup: 0,
                    reps: MAX_PRODUCTION_REPS + 1,
                },
                "production reps",
            ),
            (
                ProductionRunConfig {
                    n: MAX_PRODUCTION_ELEMENTS,
                    warmup: MAX_PRODUCTION_WARMUP,
                    reps: 1,
                },
                "n * (warmup + reps)",
            ),
        ] {
            let error = config.validate().expect_err("invalid config must fail");
            assert!(error.contains(expected), "unexpected diagnostic: {error}");
        }
        ProductionRunConfig {
            n: MAX_PRODUCTION_ELEMENTS,
            warmup: 0,
            reps: 1,
        }
        .validate()
        .expect("maximum allocation with one pass is admitted without allocating");
        ProductionRunConfig {
            n: 1,
            warmup: MAX_PRODUCTION_WARMUP,
            reps: MAX_PRODUCTION_REPS,
        }
        .validate()
        .expect("maximum repetition counts with a minimal input are admitted");
    }

    #[test]
    #[allow(clippy::too_many_lines)] // one complete retained-receipt lineage attack matrix
    fn dependency_receipt_verifier_requires_input_edge_and_intact_bytes() {
        let CitationAuthority::Receipt(binding) = test_receipt_authority() else {
            unreachable!("test helper is receipt-backed")
        };
        let db = temp_db("dependency-verifier");
        let ledger = Ledger::open(&db).expect("open ledger");
        let stored = ledger
            .put_artifact(
                crate::DEPGRAPH_RECEIPT_ARTIFACT_KIND,
                binding.bytes.as_bytes(),
                Some(crate::DEPGRAPH_RECEIPT_ARTIFACT_META),
            )
            .expect("store receipt");
        assert_eq!(stored.hash, binding.artifact_hash);
        let explicits = fs_ledger::FiveExplicits {
            seed: b"dependency-verifier",
            versions: "{}",
            budget: "{}",
            capability: "{}",
        };
        let op = ledger
            .begin_op(None, "{}", &explicits, 1)
            .expect("begin verifier fixture");
        ledger
            .finish_op(op, fs_ledger::OpOutcome::Ok, None, 2)
            .expect("finish verifier fixture");
        let placeholder = fs_blake3::hash_domain(
            "fs-roofline.dependency-verifier-placeholder.v1",
            b"placeholder",
        );
        let protocol = crate::CanonicalProductionOp {
            kernel_count: 1,
            fingerprint: 1,
            post_fingerprint: 1,
            run_nonce: placeholder,
            pre_axes_receipt: placeholder,
            post_axes_receipt: placeholder,
            dependency_receipt_digest: binding.domain_digest,
            dependency_receipt_artifact: binding.artifact_hash,
            finalized_run_receipt: placeholder,
            result_manifest: "{\"schema\":\"fs-roofline-run-manifest-v1\",\"entries\":[]}"
                .to_string(),
            baseline_admission: "{}".to_string(),
        };
        let mut missing_protocol = protocol.clone();
        missing_protocol.dependency_receipt_artifact =
            fs_blake3::hash_domain("fs-roofline.missing-dependency-receipt.v1", b"missing");
        assert!(
            !crate::dependency_receipt_is_structurally_valid(&ledger, op, &missing_protocol)
                .expect("missing-artifact verdict")
        );
        assert!(
            !crate::dependency_receipt_is_structurally_valid(&ledger, op, &protocol)
                .expect("missing-edge verdict"),
            "retained bytes without an op input edge are not lineage evidence"
        );
        ledger
            .link(op, &binding.artifact_hash, fs_ledger::EdgeRole::In)
            .expect("link receipt input");
        assert!(
            crate::dependency_receipt_is_structurally_valid(&ledger, op, &protocol)
                .expect("linked receipt verdict")
        );
        let substitute = DependencyReceiptBinding::from_parts(
            SUBSTITUTE_DEPGRAPH_RECEIPT,
            fs_blake3::hash_domain(
                fs_session::GEMM_DEPGRAPH_RECEIPT_DOMAIN,
                SUBSTITUTE_DEPGRAPH_RECEIPT.as_bytes(),
            ),
        )
        .expect("substitute digest agrees with its own bytes");
        ledger
            .put_artifact(
                crate::DEPGRAPH_RECEIPT_ARTIFACT_KIND,
                substitute.bytes.as_bytes(),
                Some(crate::DEPGRAPH_RECEIPT_ARTIFACT_META),
            )
            .expect("store internally consistent substitute");
        ledger
            .link(op, &substitute.artifact_hash, fs_ledger::EdgeRole::In)
            .expect("link substitute input");
        let mut substituted_protocol = protocol.clone();
        substituted_protocol.dependency_receipt_digest = substitute.domain_digest;
        substituted_protocol.dependency_receipt_artifact = substitute.artifact_hash;
        assert!(
            crate::dependency_receipt_is_structurally_valid(&ledger, op, &substituted_protocol)
                .expect("substitution structure verdict"),
            "a historical receipt is validated against its own retained bytes"
        );
        assert!(
            !crate::dependency_receipt_matches_binding(&substituted_protocol, Some(binding)),
            "a structurally sound historical receipt must not impersonate today's build receipt"
        );
        ledger
            .corrupt_artifact_for_test(&binding.artifact_hash)
            .expect("tamper receipt bytes");
        assert!(
            !crate::dependency_receipt_is_structurally_valid(&ledger, op, &protocol)
                .expect("tampered receipt verdict"),
            "content corruption must classify as invalid evidence, not a valid lineage edge"
        );
        cleanup_db(&db);
    }

    #[test]
    fn payload_artifact_envelope_requires_exact_kind_and_metadata() {
        let ledger = Ledger::open(":memory:").expect("in-memory ledger");
        let wrong_kind = ledger
            .put_artifact(
                "untrusted-result",
                b"wrong-kind",
                Some(crate::ROOFLINE_PAYLOAD_ARTIFACT_META),
            )
            .expect("store wrong-kind fixture");
        assert!(
            !crate::artifact_envelope_is_valid(
                &ledger,
                &wrong_kind.hash,
                crate::ROOFLINE_PAYLOAD_ARTIFACT_KIND,
                crate::ROOFLINE_PAYLOAD_ARTIFACT_META,
            )
            .expect("wrong-kind verdict")
        );

        let wrong_meta = ledger
            .put_artifact(
                crate::ROOFLINE_PAYLOAD_ARTIFACT_KIND,
                b"wrong-meta",
                Some("{\"schema\":\"attacker-controlled\"}"),
            )
            .expect("store wrong-metadata fixture");
        assert!(
            !crate::artifact_envelope_is_valid(
                &ledger,
                &wrong_meta.hash,
                crate::ROOFLINE_PAYLOAD_ARTIFACT_KIND,
                crate::ROOFLINE_PAYLOAD_ARTIFACT_META,
            )
            .expect("wrong-metadata verdict")
        );
    }

    #[test]
    fn post_probe_is_observed_strictly_after_every_timed_repetition() {
        let axes = synthetic_axes(0xB);
        let (baseline, identity) = trusted_baseline(&axes);
        let policy = AxisBaselinePolicy::new(Some(&baseline), &identity, 20_010);
        let runs = Rc::new(Cell::new(0_usize));
        let registry: Vec<Box<dyn crate::RooflineKernel>> = vec![Box::new(CountingKernel {
            runs: Rc::clone(&runs),
            value: 1,
        })];
        let probe = ProductionProbe::from_observed(axes.clone());
        let runs_at_post = Rc::new(Cell::new(usize::MAX));
        let observed = Rc::clone(&runs_at_post);
        let counter = Rc::clone(&runs);
        let config = ProductionRunConfig {
            n: 64,
            warmup: 2,
            reps: 3,
        };
        let run = probe
            .run_with_parts(config, policy, registry, move || {
                observed.set(counter.get());
                axes.clone()
            })
            .expect("protocol run");
        // warmup(2) + reps(3): the post-probe fired only after all five.
        assert_eq!(runs_at_post.get(), 5);
        assert_eq!(run.results().len(), 1);
    }

    #[test]
    fn drifted_post_probe_refuses_citation_and_records_a_rejection() {
        let axes = synthetic_axes(0xC);
        let (baseline, identity) = trusted_baseline(&axes);
        let policy = AxisBaselinePolicy::new(Some(&baseline), &identity, 20_010);
        let mut drifted = axes.clone();
        drifted.bandwidth_single_gbs *= 0.3;
        drifted.bandwidth_all_core_gbs *= 0.3;
        let probe = ProductionProbe::from_observed(axes);
        let run = probe
            .run_with_parts(
                CONFIG,
                policy,
                default_registry(1 << 10).expect("bounded registry fixture"),
                move || drifted,
            )
            .expect("protocol run");
        assert!(
            !run.citation_eligible(),
            "drifted post-probe must refuse citation eligibility"
        );
        let reason = run.admission_error().expect("admission diagnostic");
        assert!(
            reason.contains("baseline admission refused"),
            "unexpected diagnostic: {reason}"
        );

        let db = temp_db("drift");
        let ledger = Ledger::open(&db).expect("open ledger");
        let kernel = run.results()[0].kernel.clone();
        let version = run.results()[0].version.clone();
        let fingerprint = run.axes().fingerprint;
        let baseline_hash = policy.baseline_hash();
        let op = run.record(&ledger).expect("record rejection");
        let ir = ledger.op(op).unwrap().expect("op row").ir;
        assert!(ir.contains("\"protocol\":\"production-v2\""));
        assert!(ir.contains("\"admitted\":false"));
        // A rejected run publishes no tune evidence.
        assert_eq!(
            staleness_at(&ledger, &kernel, &version, fingerprint, baseline_hash, 1)
                .expect("staleness"),
            Staleness::NeverMeasured
        );
        cleanup_db(&db);
    }

    #[test]
    fn partial_finalizer_failure_yields_no_recordable_run() {
        let axes = synthetic_axes(0xD);
        let (baseline, identity) = trusted_baseline(&axes);
        let policy = AxisBaselinePolicy::new(Some(&baseline), &identity, 20_010);
        let registry: Vec<Box<dyn crate::RooflineKernel>> = vec![Box::new(FailingFinalizeKernel {
            inner: CountingKernel {
                runs: Rc::new(Cell::new(0)),
                value: 1,
            },
        })];
        let probe = ProductionProbe::from_observed(axes);
        let error = probe
            .run_with_parts(CONFIG, policy, registry, || synthetic_axes(0xD))
            .expect_err("finalizer failure must poison the whole run");
        assert!(
            error.contains("tune ledger unavailable mid-finalize"),
            "diagnostic must name the failing kernel's reason: {error}"
        );
        // No ProductionRun exists, so nothing can reach a ledger at all.
    }

    #[test]
    fn development_salt_is_report_only_even_when_measurements_admit() {
        let axes = synthetic_axes(0xD1);
        let (baseline, identity) = trusted_baseline(&axes);
        let policy = AxisBaselinePolicy::new(Some(&baseline), &identity, 20_010);
        let probe = ProductionProbe::from_observed(axes.clone());
        let post = axes.clone();
        let run = probe
            .run_with_parts_and_authority(
                CONFIG,
                policy,
                default_registry(1 << 10).expect("bounded registry fixture"),
                move || post,
                CitationAuthority::Refused(DEVELOPMENT_SALT_REFUSAL),
            )
            .expect("report-only run");
        assert!(run.finalized.admitted(), "numeric admission is independent");
        assert!(
            !run.citation_eligible(),
            "a development salt is never citation evidence"
        );
        assert_eq!(
            run.admission_error().as_deref(),
            Some(DEVELOPMENT_SALT_REFUSAL)
        );
        let kernel = run.results()[0].kernel.clone();
        let version = run.results()[0].version.clone();

        let db = temp_db("salt-refusal");
        let ledger = Ledger::open(&db).expect("open ledger");
        let op = run.record(&ledger).expect("record structured refusal");
        let row = ledger.op(op).unwrap().expect("refusal op");
        assert!(row.ir.contains("\"measurement_admitted\":true"));
        assert!(row.ir.contains("\"admitted\":false"));
        assert!(row.ir.contains("\"citation_refusal\":"));
        assert_eq!(row.outcome.as_deref(), Some("error"));
        assert_eq!(
            staleness_at(
                &ledger,
                &kernel,
                &version,
                axes.fingerprint,
                policy.baseline_hash(),
                row.t_end.expect("finished refusal"),
            )
            .expect("staleness"),
            Staleness::NeverMeasured,
        );
        cleanup_db(&db);
    }

    #[test]
    #[allow(clippy::too_many_lines)] // one end-to-end protocol and staleness state matrix
    fn successful_production_run_records_nonce_and_both_axis_receipts() {
        let axes = synthetic_axes(0xE);
        let (baseline, identity) = trusted_baseline(&axes);
        let policy = AxisBaselinePolicy::new(Some(&baseline), &identity, 20_010);
        let probe = ProductionProbe::from_observed(axes.clone());
        let nonce = probe.nonce;
        let post = axes.clone();
        let authority = test_receipt_authority();
        let CitationAuthority::Receipt(dependency_receipt) = authority else {
            unreachable!("test helper is receipt-backed")
        };
        let run = probe
            .run_with_parts_and_authority(
                CONFIG,
                policy,
                default_registry(1 << 10).expect("bounded registry fixture"),
                move || post,
                authority,
            )
            .expect("protocol run");
        assert!(
            run.citation_eligible(),
            "stable synthetic probes must pass the synthetic eligibility fixture"
        );
        assert_eq!(run.nonce(), nonce);

        let db = temp_db("ok");
        let ledger = Ledger::open(&db).expect("open ledger");
        let kernel = run.results()[0].kernel.clone();
        let version = run.results()[0].version.clone();
        let baseline_hash = policy.baseline_hash();
        let op = run.record(&ledger).expect("record production run");
        let row = ledger.op(op).unwrap().expect("op row");
        let recorded_at = row.t_end.expect("finished op");
        assert!(row.ir.contains("\"protocol\":\"production-v2\""));
        assert!(
            row.ir
                .contains("\"dependency_graph_evidence\":\"operator-observed-receipt\"")
        );
        assert!(row.ir.contains(&format!("\"run_nonce\":\"{nonce}\"")));
        assert!(
            row.ir
                .contains(&format!("\"pre_axes_receipt\":\"{}\"", axes_receipt(&axes)))
        );
        assert!(row.ir.contains(&format!(
            "\"post_axes_receipt\":\"{}\"",
            axes_receipt(&axes)
        )));
        assert!(row.ir.contains("\"admitted\":true"));
        assert!(
            ledger
                .edge_exists(
                    op,
                    &dependency_receipt.artifact_hash,
                    fs_ledger::EdgeRole::In
                )
                .expect("dependency receipt edge")
        );
        let dependency_info = ledger
            .artifact_info(&dependency_receipt.artifact_hash)
            .expect("dependency receipt metadata")
            .expect("retained dependency receipt");
        assert_eq!(dependency_info.kind, crate::DEPGRAPH_RECEIPT_ARTIFACT_KIND);
        assert_eq!(
            dependency_info.meta.as_deref(),
            Some(crate::DEPGRAPH_RECEIPT_ARTIFACT_META)
        );
        assert_eq!(
            ledger
                .get_artifact(&dependency_receipt.artifact_hash)
                .expect("dependency receipt bytes")
                .as_deref(),
            Some(TEST_DEPGRAPH_RECEIPT.as_bytes())
        );
        assert_eq!(
            receipt_staleness_at(
                &ledger,
                &kernel,
                &version,
                axes.fingerprint,
                baseline_hash,
                recorded_at + STALENESS_MAX_AGE_NS,
                dependency_receipt,
            )
            .expect("staleness"),
            Staleness::Fresh
        );
        assert_eq!(
            receipt_staleness_at(
                &ledger,
                &kernel,
                &version,
                axes.fingerprint,
                baseline_hash,
                recorded_at + STALENESS_MAX_AGE_NS + 1,
                dependency_receipt,
            )
            .expect("expired staleness"),
            Staleness::Expired
        );
        assert_eq!(
            receipt_staleness_at(
                &ledger,
                &kernel,
                &version,
                axes.fingerprint,
                baseline_hash,
                recorded_at - 1,
                dependency_receipt,
            )
            .expect("clock rollback staleness"),
            Staleness::ClockRollback
        );
        assert_eq!(
            staleness(&ledger, &kernel, &version, axes.fingerprint, None)
                .expect("missing baseline staleness"),
            Staleness::BaselineUnavailable
        );
        let foreign_baseline = fs_blake3::hash_domain("fs-roofline.foreign-baseline.v1", b"other");
        assert_eq!(
            staleness(
                &ledger,
                &kernel,
                &version,
                axes.fingerprint,
                Some(foreign_baseline),
            )
            .expect("baseline drift staleness"),
            Staleness::BaselineDrift
        );
        assert_eq!(
            staleness(
                &ledger,
                &kernel,
                "different-version",
                axes.fingerprint,
                baseline_hash,
            )
            .expect("version drift staleness"),
            Staleness::NeverMeasured
        );
        assert_eq!(
            staleness(&ledger, &kernel, &version, 0xFFFF, baseline_hash)
                .expect("fingerprint drift staleness"),
            Staleness::FingerprintDrift
        );

        let current_row = ledger
            .tune_rows(&kernel)
            .expect("production tune rows")
            .into_iter()
            .find(|row| {
                row.machine == roofline_machine_key(axes.fingerprint, baseline.content_hash())
            })
            .expect("current production row");
        ledger
            .tune_put(
                &current_row.kernel,
                &current_row.shape_class,
                &current_row.machine,
                &current_row.params,
                "{}",
            )
            .expect("inject valid-JSON payload corruption");
        assert_eq!(
            receipt_staleness_at(
                &ledger,
                &kernel,
                &version,
                axes.fingerprint,
                baseline_hash,
                recorded_at + 1,
                dependency_receipt,
            )
            .expect("corrupt production staleness"),
            Staleness::CorruptEvidence
        );
        cleanup_db(&db);
    }

    struct RecordedManifestRun {
        ledger: Ledger,
        baseline: BaselineAxes,
        kernels: Vec<(String, String)>,
        recorded_at: i64,
        dependency: DependencyReceiptBinding,
    }

    fn recorded_manifest_run(db: &str) -> RecordedManifestRun {
        let ledger = Ledger::open(db).expect("open ledger");
        let axes = synthetic_axes(0xBEEF);
        let (baseline, identity) = trusted_baseline(&axes);
        let policy = AxisBaselinePolicy::new(Some(&baseline), &identity, 20_010);
        let probe = ProductionProbe::from_observed(axes.clone());
        let post = axes.clone();
        let run = probe
            .run_with_parts_and_authority(
                CONFIG,
                policy,
                default_registry(1 << 10).expect("bounded registry fixture"),
                move || post,
                test_receipt_authority(),
            )
            .expect("sealed manifest fixture");
        assert!(run.citation_eligible());
        let kernels = run
            .results()
            .iter()
            .map(|result| (result.kernel.clone(), result.version.clone()))
            .collect();
        let op = run.record(&ledger).expect("record sealed manifest fixture");
        let recorded_at = ledger
            .op(op)
            .unwrap()
            .expect("recorded op")
            .t_end
            .expect("finished op");
        RecordedManifestRun {
            ledger,
            baseline,
            kernels,
            recorded_at,
            dependency: test_receipt_binding(),
        }
    }

    #[test]
    fn custom_registry_history_cannot_poison_a_fresh_production_row() {
        let db = temp_db("custom-history-isolation");
        let production = recorded_manifest_run(&db);
        let axes = synthetic_axes(0xBEEF);
        let identity = BaselineIdentity::current(&axes, "test-firmware")
            .expect("synthetic identity agrees with the retained baseline");
        let policy = AxisBaselinePolicy::new(Some(&production.baseline), &identity, 20_010);
        let mut registry = default_registry(1 << 10).expect("bounded registry fixture");
        let mut results =
            run_registry(&mut registry, 0, 1, &axes).expect("bounded exploratory registry run");
        let mut finalized = finalize_registry_tuning(&mut registry, &axes, &axes, policy, &results)
            .expect("finalize exploratory run");
        crate::record_run(
            &production.ledger,
            &axes,
            &axes,
            policy,
            &mut finalized,
            &mut results,
        )
        .expect("record exploratory row in its candidate namespace");

        let (kernel, version) = &production.kernels[0];
        let rows = production
            .ledger
            .tune_rows(kernel)
            .expect("query production and candidate rows");
        assert_eq!(
            rows.iter()
                .filter(|row| row.shape_class.starts_with(crate::TUNE_SHAPE_CLASS))
                .count(),
            1
        );
        assert_eq!(
            rows.iter()
                .filter(|row| row.shape_class.starts_with(crate::CUSTOM_TUNE_SHAPE_CLASS))
                .count(),
            1
        );
        assert_eq!(
            manifest_probe(&production, kernel, version),
            Staleness::Fresh,
            "candidate history must not enter the production staleness scan"
        );
        cleanup_db(&db);
    }

    fn manifest_probe(run: &RecordedManifestRun, kernel: &str, version: &str) -> Staleness {
        receipt_staleness_at(
            &run.ledger,
            kernel,
            version,
            0xBEEF,
            Some(run.baseline.content_hash()),
            run.recorded_at + 1,
            run.dependency,
        )
        .expect("manifest staleness probe")
    }

    fn roofline_row(ledger: &Ledger, kernel: &str) -> fs_ledger::TuneRow {
        let mut rows: Vec<_> = ledger
            .tune_rows(kernel)
            .expect("tune rows")
            .into_iter()
            .filter(|row| row.shape_class.contains(":run="))
            .collect();
        assert_eq!(rows.len(), 1, "expected one roofline row for {kernel}");
        rows.pop().expect("row")
    }

    fn splice_payload(ledger: &Ledger, row: &fs_ledger::TuneRow, new_measured: &str) {
        let old_hash = fs_ledger::hash_bytes(row.measured.as_bytes()).to_string();
        let new_hash = fs_ledger::hash_bytes(new_measured.as_bytes());
        let artifact = ledger
            .put_artifact(
                crate::ROOFLINE_PAYLOAD_ARTIFACT_KIND,
                new_measured.as_bytes(),
                Some("{\"schema\":\"fs-roofline-benchmark-result-v1\"}"),
            )
            .expect("store forged artifact");
        assert_eq!(artifact.hash, new_hash);
        let op: i64 = row
            .params
            .split_once("\"op\":")
            .and_then(|(_, rest)| rest.split_once(','))
            .and_then(|(digits, _)| digits.parse().ok())
            .expect("op id in params");
        ledger
            .link(op, &new_hash, fs_ledger::EdgeRole::Out)
            .expect("forged edge");
        let forged_params = row.params.replace(&old_hash, &new_hash.to_string());
        assert_ne!(forged_params, row.params);
        ledger
            .tune_put(
                &row.kernel,
                &row.shape_class,
                &row.machine,
                &forged_params,
                new_measured,
            )
            .expect("overwrite row");
    }

    fn altered_measured(measured: &str) -> String {
        let (before, after) = measured
            .split_once("\"dispersion\":")
            .expect("dispersion field");
        let end = after.find([',', '}']).expect("field end");
        let forged = format!("{before}\"dispersion\":9.5e-1{}", &after[end..]);
        assert_ne!(forged, measured);
        forged
    }

    #[test]
    fn manifest_replacement_poisons_the_row_and_its_siblings() {
        let db = temp_db("manifest-splice");
        let run = recorded_manifest_run(&db);
        let (kernel_a, version_a) = run.kernels[0].clone();
        let (kernel_b, version_b) = run.kernels[1].clone();
        assert_eq!(
            manifest_probe(&run, &kernel_a, &version_a),
            Staleness::Fresh
        );
        assert_eq!(
            manifest_probe(&run, &kernel_b, &version_b),
            Staleness::Fresh
        );

        let row = roofline_row(&run.ledger, &kernel_a);
        splice_payload(&run.ledger, &row, &altered_measured(&row.measured));
        assert_eq!(
            manifest_probe(&run, &kernel_a, &version_a),
            Staleness::CorruptEvidence
        );
        assert_eq!(
            manifest_probe(&run, &kernel_b, &version_b),
            Staleness::CorruptEvidence
        );
        cleanup_db(&db);
    }

    #[test]
    fn sibling_parameter_tamper_poisons_every_manifest_member() {
        let db = temp_db("manifest-sibling-params");
        let run = recorded_manifest_run(&db);
        let (kernel_a, _) = run.kernels[0].clone();
        let (kernel_b, version_b) = run.kernels[1].clone();
        assert_eq!(
            manifest_probe(&run, &kernel_b, &version_b),
            Staleness::Fresh
        );

        let row = roofline_row(&run.ledger, &kernel_a);
        let tampered_params = row.params.replace("\"reps\":1,", "\"reps\":2,");
        assert_ne!(
            tampered_params, row.params,
            "fixture must alter sibling params"
        );
        run.ledger
            .tune_put(
                &row.kernel,
                &row.shape_class,
                &row.machine,
                &tampered_params,
                &row.measured,
            )
            .expect("overwrite sibling params");

        assert_eq!(
            manifest_probe(&run, &kernel_b, &version_b),
            Staleness::CorruptEvidence,
            "querying an untouched row must still validate every sibling's canonical params"
        );
        cleanup_db(&db);
    }

    #[test]
    fn sibling_artifact_corruption_poisons_every_manifest_member() {
        let db = temp_db("manifest-sibling-artifact");
        let run = recorded_manifest_run(&db);
        let (kernel_a, _) = run.kernels[0].clone();
        let (kernel_b, version_b) = run.kernels[1].clone();
        let row = roofline_row(&run.ledger, &kernel_a);
        let params = crate::parse_roofline_row_params(&row.params).expect("canonical row params");
        run.ledger
            .corrupt_artifact_for_test(&params.payload_artifact)
            .expect("corrupt sibling artifact bytes");

        assert_eq!(
            manifest_probe(&run, &kernel_b, &version_b),
            Staleness::CorruptEvidence,
            "an untouched row cannot stay Fresh when a sibling artifact is corrupt"
        );
        cleanup_db(&db);
    }

    #[test]
    fn production_operation_parser_rejects_noncanonical_and_ambiguous_ir() {
        let db = temp_db("canonical-operation-parser");
        let run = recorded_manifest_run(&db);
        let (kernel, _) = run.kernels[0].clone();
        let row = roofline_row(&run.ledger, &kernel);
        let params = crate::parse_roofline_row_params(&row.params).expect("canonical row params");
        let ir = run
            .ledger
            .op(params.op)
            .expect("query op")
            .expect("recorded op")
            .ir;
        let parsed = crate::parse_canonical_production_op(&ir).expect("canonical production IR");
        assert_eq!(parsed.to_json(), ir);
        assert!(
            crate::validate_protocol_axes(&parsed, 0xBEEF, run.baseline.content_hash(),).is_some(),
            "canonical pre/post receipts must bind the recorded fingerprints, axes, and baseline"
        );

        let mut substituted_axes_receipt = parsed.clone();
        substituted_axes_receipt.pre_axes_receipt =
            fs_blake3::hash_domain("fs-roofline.substituted-axes-receipt.v1", b"substitute");
        assert!(
            crate::validate_protocol_axes(
                &substituted_axes_receipt,
                0xBEEF,
                run.baseline.content_hash(),
            )
            .is_none(),
            "a canonical operation cannot substitute a different pre-probe receipt"
        );

        let duplicate = ir.replacen(
            "\"measurement_admitted\":true,",
            "\"measurement_admitted\":true,\"measurement_admitted\":true,",
            1,
        );
        assert!(crate::parse_canonical_production_op(&duplicate).is_none());
        assert!(crate::parse_canonical_production_op(&format!("{ir} ")).is_none());
        let reordered = ir.replacen(
            "\"measurement_admitted\":true,\"admitted\":true",
            "\"admitted\":true,\"measurement_admitted\":true",
            1,
        );
        assert!(crate::parse_canonical_production_op(&reordered).is_none());
        cleanup_db(&db);
    }

    #[test]
    fn manifest_rejects_rows_added_after_finalization() {
        let db = temp_db("manifest-added");
        let run = recorded_manifest_run(&db);
        let (kernel, version) = run.kernels[0].clone();
        let row = roofline_row(&run.ledger, &kernel);
        let ghost = "ghost-kernel";
        let ghost_measured = row.measured.replace(
            &format!("\"kernel\":\"{kernel}\""),
            &format!("\"kernel\":\"{ghost}\""),
        );
        let ghost_hash = fs_ledger::hash_bytes(ghost_measured.as_bytes());
        run.ledger
            .put_artifact(
                crate::ROOFLINE_PAYLOAD_ARTIFACT_KIND,
                ghost_measured.as_bytes(),
                Some("{\"schema\":\"fs-roofline-benchmark-result-v1\"}"),
            )
            .expect("store ghost artifact");
        let op: i64 = row
            .params
            .split_once("\"op\":")
            .and_then(|(_, rest)| rest.split_once(','))
            .and_then(|(digits, _)| digits.parse().ok())
            .expect("op id in params");
        run.ledger
            .link(op, &ghost_hash, fs_ledger::EdgeRole::Out)
            .expect("ghost edge");
        let ghost_params = row.params.replace(
            &fs_ledger::hash_bytes(row.measured.as_bytes()).to_string(),
            &ghost_hash.to_string(),
        );
        run.ledger
            .tune_put(
                ghost,
                &row.shape_class,
                &row.machine,
                &ghost_params,
                &ghost_measured,
            )
            .expect("insert ghost row");
        assert_eq!(
            manifest_probe(&run, ghost, &version),
            Staleness::CorruptEvidence
        );
        assert_eq!(manifest_probe(&run, &kernel, &version), Staleness::Fresh);
        cleanup_db(&db);
    }

    #[test]
    fn identical_receipt_backed_rerun_history_stays_fresh() {
        let db = temp_db("manifest-rerun");
        let first = recorded_manifest_run(&db);
        let (first_kernel, first_version) = first.kernels[0].clone();
        assert_eq!(
            manifest_probe(&first, &first_kernel, &first_version),
            Staleness::Fresh
        );

        let axes = synthetic_axes(0xBEEF);
        let (baseline, identity) = trusted_baseline(&axes);
        let policy = AxisBaselinePolicy::new(Some(&baseline), &identity, 20_010);
        let probe = ProductionProbe::from_observed(axes.clone());
        let post = axes.clone();
        let run = probe
            .run_with_parts_and_authority(
                CONFIG,
                policy,
                default_registry(1 << 10).expect("bounded registry fixture"),
                move || post,
                test_receipt_authority(),
            )
            .expect("second sealed run");
        let second_op = run.record(&first.ledger).expect("record second run");
        let rerecorded_at = first
            .ledger
            .op(second_op)
            .unwrap()
            .expect("second op")
            .t_end
            .expect("finished second op");
        for (kernel, version) in &first.kernels {
            assert_eq!(
                receipt_staleness_at(
                    &first.ledger,
                    kernel,
                    version,
                    0xBEEF,
                    Some(first.baseline.content_hash()),
                    rerecorded_at + 1,
                    first.dependency,
                )
                .expect("staleness probe"),
                Staleness::Fresh
            );
        }
        cleanup_db(&db);
    }

    #[test]
    fn pre_dependency_receipt_rows_are_retired_as_corrupt() {
        let db = temp_db("manifest-v3-row");
        let run = recorded_manifest_run(&db);
        let (kernel, version) = run.kernels[0].clone();
        assert_eq!(manifest_probe(&run, &kernel, &version), Staleness::Fresh);
        let row = roofline_row(&run.ledger, &kernel);
        let old_params = row.params.replace(
            "\"schema\":\"fs-roofline-ledger-row-v4\"",
            "\"schema\":\"fs-roofline-ledger-row-v3\"",
        );
        assert_ne!(old_params, row.params);
        run.ledger
            .tune_put(
                &row.kernel,
                &row.shape_class,
                &row.machine,
                &old_params,
                &row.measured,
            )
            .expect("downgrade row schema");
        assert_eq!(
            manifest_probe(&run, &kernel, &version),
            Staleness::CorruptEvidence
        );
        cleanup_db(&db);
    }

    fn cleanup_db(path: &str) {
        for suffix in ["", "-wal", "-shm", ".fsqlite-wal", ".fsqlite-shm"] {
            let _ = std::fs::remove_file(format!("{path}{suffix}"));
        }
    }
}
