//! fs-roofline: the roofline harness (plan §14; Decalogue P6).
//!
//! Performance claims as FALSIFIABLE targets: every registered kernel is
//! benchmarked against its arithmetic-intensity-derived limit on the actual
//! machine — measured axes, never spec-sheet numbers — with dispersion
//! reported and results ledgered under the machine fingerprint. "A target
//! that was never re-measured is a lie waiting to happen."
//!
//! Layer: L6 (consumes fs-substrate probes, fs-simd primitives, and writes
//! fs-ledger records). Reporting-only in v0: attainment verdicts inform;
//! gating bands belong to nightly runs on ledgered reference machines.

pub mod authority;
pub mod axes;
pub mod baseline;
pub mod kernels;
pub mod stats;

pub use authority::{
    KeyVerdict, NoPromotionAuthority, PromotionAttestation, PromotionAuthorityVerifier,
    StaticKeyRegistry,
};
pub use axes::MachineAxes;
pub use baseline::{
    AttestedBaselineStore, BASELINE_SCHEMA_VERSION, BaselineAxes, BaselineCandidate,
    BaselineClockError, BaselineIdentity, BaselineProvenance, BaselineStore, BaselineVerdict,
    PromotionError, citable_axis_admission, citable_axis_admission_authorized,
    days_since_epoch_now, promote_baseline,
};

use fs_ledger::{EdgeRole, EventRow, FiveExplicits, Ledger, LedgerError, OpOutcome, now_wall_ns};

pub mod regress;

/// Crate version (compile-time stamp).
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

const FINALIZED_RUN_DOMAIN: &str = "org.frankensim.fs-roofline.finalized-run.v2";
const RESULT_MANIFEST_DOMAIN: &str = "org.frankensim.fs-roofline.run-result-manifest.v1";
const RESULT_MANIFEST_SCHEMA: &str = "fs-roofline-run-manifest-v1";

/// One-shot capability proving that every registry kernel observed the exact
/// run's aggregate admission decision and completed its tuning lifecycle.
/// Fields are private and the value is deliberately neither `Clone` nor `Copy`.
#[derive(Debug)]
pub struct FinalizedRegistryRun {
    receipt: fs_blake3::ContentHash,
    admitted: bool,
    consumed: bool,
}

impl FinalizedRegistryRun {
    /// Whether the finalized run passed citable admission.
    #[must_use]
    pub fn admitted(&self) -> bool {
        self.admitted
    }

    /// Content identity of the exact axes, baseline decision, and result set
    /// seen during finalization.
    #[must_use]
    pub fn receipt_identity(&self) -> fs_blake3::ContentHash {
        self.receipt
    }
}

/// Explicit historical-baseline policy for one roofline run.
///
/// `None` is not a permissive default: it represents a first/unbaselined run,
/// which may be reported as candidate evidence but cannot publish citable
/// metrics or tune rows. Without [`AxisBaselinePolicy::with_authority`], the
/// referenced store is an OPERATOR-TRUSTED root (tamper-evident, not
/// independently verified); binding gates bind an attestation and an
/// injected [`PromotionAuthorityVerifier`] (bead fz2.7), which re-verifies
/// the promotion signature before every citable decision.
#[derive(Clone, Copy)]
pub struct AxisBaselinePolicy<'a> {
    baseline: Option<&'a BaselineAxes>,
    identity: &'a BaselineIdentity,
    now_day: u64,
    authority: Option<(
        Option<&'a PromotionAttestation>,
        &'a dyn PromotionAuthorityVerifier,
    )>,
}

impl<'a> AxisBaselinePolicy<'a> {
    /// Bind the selected baseline (if any), declared current environment, and
    /// observed epoch day used for age admission.
    #[must_use]
    pub const fn new(
        baseline: Option<&'a BaselineAxes>,
        identity: &'a BaselineIdentity,
        now_day: u64,
    ) -> Self {
        Self {
            baseline,
            identity,
            now_day,
            authority: None,
        }
    }

    /// Bind the promotion attestation and an injected authority (bead
    /// fz2.7): the verdict then requires an AUTHORIZED signature over
    /// the baseline's content hash before any band math, and the
    /// receipt binds the verifying key identity.
    #[must_use]
    pub fn with_authority(
        mut self,
        attestation: Option<&'a PromotionAttestation>,
        authority: &'a dyn PromotionAuthorityVerifier,
    ) -> Self {
        self.authority = Some((attestation, authority));
        self
    }

    /// Evaluate the complete pre/probe/post baseline policy.
    #[must_use]
    pub fn verdict(&self, pre: &MachineAxes, post: &MachineAxes) -> BaselineVerdict {
        match self.authority {
            Some((attestation, authority)) => citable_axis_admission_authorized(
                pre,
                post,
                self.baseline,
                attestation,
                self.identity,
                self.now_day,
                authority,
            ),
            None => citable_axis_admission(pre, post, self.baseline, self.identity, self.now_day),
        }
    }

    /// Domain-separated identity of the selected baseline, if one exists.
    #[must_use]
    pub fn baseline_hash(&self) -> Option<fs_blake3::ContentHash> {
        self.baseline.map(BaselineAxes::content_hash)
    }

    /// Canonical, self-contained receipt for the baseline admission decision.
    #[must_use]
    pub fn receipt_json(&self, pre: &MachineAxes, post: &MachineAxes) -> String {
        let baseline = self
            .baseline
            .map_or_else(|| "null".to_string(), BaselineAxes::canonical_json);
        let baseline_hash = self.baseline_hash().map_or_else(
            || "null".to_string(),
            |hash| format!("\"{}\"", hash.to_hex()),
        );
        // Authority binding (bead fz2.7): the verifying key identity (or
        // the explicit operator-trusted tier) travels in the receipt.
        let authority = match self.authority {
            None => "\"operator-trusted\"".to_string(),
            Some((attestation, _)) => attestation.map_or_else(
                || "{\"key_id\":null}".to_string(),
                |a| format!("{{\"key_id\":\"{}\"}}", json_escape(a.key_id())),
            ),
        };
        format!(
            "{{\"schema\":\"fs-roofline-axis-admission-v1\",\"now_day\":{},\"identity\":{},\"pre\":{},\"post\":{},\"baseline_hash\":{},\"baseline\":{},\"authority\":{},\"verdict\":{}}}",
            self.now_day,
            baseline_identity_json(self.identity),
            machine_axes_receipt_json(pre),
            machine_axes_receipt_json(post),
            baseline_hash,
            baseline,
            authority,
            self.verdict(pre, post).to_jsonl(),
        )
    }
}

fn baseline_identity_json(identity: &BaselineIdentity) -> String {
    format!(
        "{{\"fingerprint\":\"{:016x}\",\"cpu_brand\":\"{}\",\"logical_cpus\":{},\"os\":\"{}\",\"arch\":\"{}\",\"firmware\":\"{}\"}}",
        identity.fingerprint(),
        json_escape(identity.cpu_brand()),
        identity.logical_cpus(),
        json_escape(identity.os()),
        json_escape(identity.arch()),
        json_escape(identity.firmware()),
    )
}

fn machine_axes_receipt_json(axes: &MachineAxes) -> String {
    format!(
        "{{\"fingerprint\":\"{:016x}\",\"cpu_brand\":\"{}\",\"logical_cpus\":{},\"bandwidth_single_bits\":\"{:016x}\",\"bandwidth_all_core_bits\":\"{:016x}\",\"peak_single_bits\":\"{:016x}\",\"peak_all_core_bits\":\"{:016x}\"}}",
        axes.fingerprint,
        json_escape(&axes.cpu_brand),
        axes.logical_cpus,
        axes.bandwidth_single_gbs.to_bits(),
        axes.bandwidth_all_core_gbs.to_bits(),
        axes.peak_single_gflops.to_bits(),
        axes.peak_all_core_gflops.to_bits(),
    )
}

/// Shape-class prefix under which versioned roofline rows land in the ledger
/// `tune` table.
pub const TUNE_SHAPE_CLASS: &str = "roofline-v6";

/// Versioned ledger shape-class key for a kernel implementation.
#[must_use]
pub fn tune_shape_class(version: &str) -> String {
    format!("{TUNE_SHAPE_CLASS}:{version}")
}

/// Append-only shape key for one finalized measurement operation.
#[must_use]
pub fn tune_measurement_shape_class(
    version: &str,
    run_receipt: fs_blake3::ContentHash,
    op: i64,
) -> String {
    format!("{}:run={run_receipt}:op={op}", tune_shape_class(version))
}

/// Composite machine key for baseline-bound roofline rows.
#[must_use]
pub fn roofline_machine_key(fingerprint: u64, baseline: fs_blake3::ContentHash) -> [u8; 40] {
    let mut key = [0_u8; 40];
    key[..8].copy_from_slice(&fingerprint.to_le_bytes());
    key[8..].copy_from_slice(baseline.as_bytes());
    key
}

/// Which machine axis a kernel is measured against.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Threading {
    /// One thread: per-core bandwidth/compute axes.
    SingleThread,
    /// All logical cores: aggregate axes.
    AllCore,
}

/// Machine quantity against which a declared target fraction is evaluated.
///
/// This is independent of [`RoofSide`]: the binding roof remains part of every
/// report, while a plan target such as "75% of measured peak FLOPs" must stay
/// compute-relative even on a machine where memory bandwidth binds first.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TargetAxis {
    /// Fraction of `min(bandwidth limit, compute limit)`.
    BindingRoof,
    /// Fraction of the measured floating-point compute peak.
    ComputePeak,
    /// Fraction of the measured memory-bandwidth peak.
    MemoryBandwidth,
}

impl TargetAxis {
    /// Stable receipt spelling.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::BindingRoof => "binding_roof",
            Self::ComputePeak => "compute_peak",
            Self::MemoryBandwidth => "memory_bandwidth",
        }
    }
}

/// Static description of a benchmarkable kernel: identity plus the
/// arithmetic-intensity model that derives its machine-specific limit.
#[derive(Debug, Clone, Copy)]
pub struct KernelSpec {
    /// Registry name (ledger key; kebab-case).
    pub name: &'static str,
    /// Kernel version (bumped when the implementation changes — attainment
    /// history is only comparable within one version).
    pub version: &'static str,
    /// Bytes moved to/from memory per element processed.
    pub bytes_per_elem: f64,
    /// Floating-point operations per element processed.
    pub flops_per_elem: f64,
    /// Measurement threading model.
    pub threading: Threading,
    /// Axis whose measured value defines `target_fraction`.
    pub target_axis: TargetAxis,
    /// Target as a fraction of [`KernelSpec::target_axis`]. `None` =
    /// report-only, no band claimed.
    pub target_fraction: Option<f64>,
}

/// One benchmarkable kernel: owns its buffers; `run_once` is the timed unit.
pub trait RooflineKernel {
    /// The kernel's spec (identity + intensity model + target).
    fn spec(&self) -> KernelSpec;
    /// Elements processed per `run_once` call.
    fn elements(&self) -> usize;
    /// Execute one timed repetition.
    fn run_once(&mut self);
    /// Sealed execution decision made by the most recent `run_once`, when this
    /// kernel has decision provenance beyond its static [`KernelSpec`].
    ///
    /// The default is `None`; the production GEMM route returns a binding that
    /// includes its exact tune key, plan, source, implementation tier, build,
    /// and validated tune-row identity.
    fn execution_binding(&self) -> Option<KernelExecutionBinding> {
        None
    }
    /// Newly measured, validated tune state awaiting the enclosing evidence
    /// transaction. Returning it does not publish it.
    fn pending_tune_publication(&self) -> Option<fs_session::ValidatedGemmTuneRow> {
        None
    }
    /// Complete process-local tuning state after the enclosing registry run's
    /// aggregate admission decision. The measured [`Attainment`] already owns
    /// any pending publication clone; kernels discard their pending marker on
    /// either outcome and additionally invalidate rejected local decisions.
    ///
    /// # Errors
    /// Returns a structured diagnostic if lifecycle cleanup fails. The default
    /// is a no-op for kernels without tune state.
    fn finalize_tuning(&mut self, _admitted: bool) -> Result<(), String> {
        Ok(())
    }
}

/// A sealed execution-decision receipt sampled after one timed repetition.
///
/// Fields are deliberately private. Built-in kernels construct this only from
/// producer-validated values; consumers can inspect the canonical JSON emitted
/// by [`Attainment::to_jsonl`] but cannot mint citable bindings directly.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct KernelExecutionBinding {
    gemm: GemmDecisionBinding,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct GemmDecisionBinding {
    scoped_tune_key: String,
    shape_class: String,
    canonical_plan: String,
    source: &'static str,
    operation_tier: String,
    build_identity: String,
    tune_row_identity: fs_ledger::ContentHash,
    validated_row: fs_session::ValidatedGemmTuneRow,
    execution_path: fs_session::GemmExecutionReceipt,
    execution_path_identity: fs_ledger::ContentHash,
}

const GEMM_EXECUTION_PATH_DOMAIN: &str = "org.frankensim.fs-roofline.gemm-execution-path.v3";
const EXECUTION_BINDING_DOMAIN: &str = "org.frankensim.fs-roofline.execution-binding.v3";

fn execution_path_is_complete(path: &fs_session::GemmExecutionReceipt) -> bool {
    path.is_complete()
}

fn execution_path_json(path: &fs_session::GemmExecutionReceipt) -> String {
    let panels = path
        .panels
        .iter()
        .map(|panel| {
            format!(
                "{{\"kernel\":\"{}\",\"mode\":\"{}\",\"declared_run\":{},\"completed\":{},\"total\":{}}}",
                json_escape(&panel.kernel),
                json_escape(&panel.mode),
                panel.declared_run,
                panel.completed,
                panel.total,
            )
        })
        .collect::<Vec<_>>()
        .join(",");
    format!(
        "{{\"declared_run\":{},\"completed_tiles\":{},\"total_tiles\":{},\"panel_count\":{},\"panels\":[{}]}}",
        path.declared_run,
        path.completed_tiles,
        path.total_tiles,
        path.panels.len(),
        panels,
    )
}

impl KernelExecutionBinding {
    #[allow(clippy::too_many_arguments)] // every independent provenance field stays explicit
    pub(crate) fn gemm(
        scoped_tune_key: String,
        shape_class: String,
        canonical_plan: String,
        source: &'static str,
        operation_tier: String,
        build_identity: String,
        validated_row: fs_session::ValidatedGemmTuneRow,
        machine: u64,
        execution_path: fs_session::GemmExecutionReceipt,
    ) -> Result<Self, String> {
        if !validated_row.matches_decision(&scoped_tune_key, &shape_class, machine, &canonical_plan)
        {
            return Err(
                "validated GEMM tune row does not bind the dispatched decision".to_string(),
            );
        }
        let tune_row_identity = validated_row.receipt_identity();
        let execution_path_identity = fs_blake3::hash_domain(
            GEMM_EXECUTION_PATH_DOMAIN,
            execution_path_json(&execution_path).as_bytes(),
        );
        let binding = Self {
            gemm: GemmDecisionBinding {
                scoped_tune_key,
                shape_class,
                canonical_plan,
                source,
                operation_tier,
                build_identity,
                tune_row_identity,
                validated_row,
                execution_path,
                execution_path_identity,
            },
        };
        if !binding.is_valid_for(machine) {
            return Err("GEMM execution binding is internally inconsistent".to_string());
        }
        Ok(binding)
    }

    fn is_valid_for(&self, machine: u64) -> bool {
        let gemm = &self.gemm;
        !gemm.scoped_tune_key.is_empty()
            && !gemm.shape_class.is_empty()
            && !gemm.canonical_plan.is_empty()
            && gemm.source == "tuned"
            && !gemm.operation_tier.is_empty()
            && !gemm.build_identity.is_empty()
            && gemm
                .scoped_tune_key
                .contains(&format!("/shape={}/", gemm.shape_class))
            && gemm
                .scoped_tune_key
                .contains(&format!("/tier={}/", gemm.operation_tier))
            && gemm
                .scoped_tune_key
                .ends_with(&format!("/build={}", gemm.build_identity))
            && gemm.tune_row_identity == gemm.validated_row.receipt_identity()
            && execution_path_is_complete(&gemm.execution_path)
            && gemm.execution_path_identity
                == fs_blake3::hash_domain(
                    GEMM_EXECUTION_PATH_DOMAIN,
                    execution_path_json(&gemm.execution_path).as_bytes(),
                )
            && gemm.validated_row.matches_decision(
                &gemm.scoped_tune_key,
                &gemm.shape_class,
                machine,
                &gemm.canonical_plan,
            )
    }

    fn declared_run(&self) -> u64 {
        self.gemm.execution_path.declared_run
    }

    fn stable_equivalent(&self, other: &Self) -> bool {
        let left = &self.gemm;
        let right = &other.gemm;
        left.scoped_tune_key == right.scoped_tune_key
            && left.shape_class == right.shape_class
            && left.canonical_plan == right.canonical_plan
            && left.source == right.source
            && left.operation_tier == right.operation_tier
            && left.build_identity == right.build_identity
            && left.tune_row_identity == right.tune_row_identity
            && left.validated_row == right.validated_row
            && execution_path_shape_eq(&left.execution_path, &right.execution_path)
    }

    fn canonical_json(&self) -> String {
        let gemm = &self.gemm;
        format!(
            "{{\"kind\":\"gemm-v3\",\"scoped_tune_key\":\"{}\",\"shape_class\":\"{}\",\"plan\":\"{}\",\"source\":\"{}\",\"operation_tier\":\"{}\",\"build_identity\":\"{}\",\"tune_row_identity\":\"{}\",\"tune_row\":{},\"execution_path_identity\":\"{}\",\"execution_path\":{}}}",
            json_escape(&gemm.scoped_tune_key),
            json_escape(&gemm.shape_class),
            json_escape(&gemm.canonical_plan),
            gemm.source,
            json_escape(&gemm.operation_tier),
            json_escape(&gemm.build_identity),
            gemm.tune_row_identity,
            gemm.validated_row.receipt_json(),
            gemm.execution_path_identity,
            execution_path_json(&gemm.execution_path),
        )
    }

    fn receipt_identity(&self) -> fs_ledger::ContentHash {
        fs_blake3::hash_domain(EXECUTION_BINDING_DOMAIN, self.canonical_json().as_bytes())
    }

    fn validated_row(&self) -> &fs_session::ValidatedGemmTuneRow {
        &self.gemm.validated_row
    }
}

fn execution_path_shape_eq(
    left: &fs_session::GemmExecutionReceipt,
    right: &fs_session::GemmExecutionReceipt,
) -> bool {
    left.is_complete()
        && right.is_complete()
        && left.completed_tiles == right.completed_tiles
        && left.total_tiles == right.total_tiles
        && left.panels.len() == right.panels.len()
        && left.panels.iter().zip(&right.panels).all(|(left, right)| {
            left.kernel == right.kernel
                && left.mode == right.mode
                && left.completed == right.completed
                && left.total == right.total
        })
}

/// Which side of the roofline binds the limit.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RoofSide {
    /// Memory-bandwidth-bound at this intensity.
    Bandwidth,
    /// Compute-bound at this intensity.
    Compute,
}

/// Verdict against the kernel's declared band.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Verdict {
    /// Attainment ≥ target fraction.
    WithinBand,
    /// Attainment < target fraction.
    BelowBand,
    /// No target declared: report-only.
    NoTarget,
    /// The MEASUREMENT ENVIRONMENT is invalid: the probed axes fail
    /// absolute plausibility floors, or the kernel "beat" its roofline
    /// by an impossible margin (stale/contention-crushed axes). A gate
    /// must never pass — or fail — on a machine that was useless while
    /// it was measured (bead 1n61: a load-68 window collapsed both the
    /// STREAM probe and the kernel ~1000× together, and the RATIO
    /// self-normalized to a vacuous within_band).
    EnvironmentInvalid,
}

impl Verdict {
    /// Stable lowercase name for logs and ledger rows.
    #[must_use]
    pub fn name(self) -> &'static str {
        match self {
            Verdict::WithinBand => "within_band",
            Verdict::BelowBand => "below_band",
            Verdict::NoTarget => "no_target",
            Verdict::EnvironmentInvalid => "environment_invalid",
        }
    }
}

/// One kernel's measured attainment against the machine roofline.
#[derive(Debug, Clone)]
pub struct Attainment {
    /// Kernel name.
    pub kernel: String,
    /// Kernel version.
    pub version: String,
    /// Median elements/second across repetitions.
    pub elems_per_sec: f64,
    /// Achieved memory traffic, GB/s.
    pub achieved_gbs: f64,
    /// Achieved compute, GFLOP/s.
    pub achieved_gflops: f64,
    /// Roofline limit in elements/second for this machine + intensity.
    pub limit_elems_per_sec: f64,
    /// Which axis binds.
    pub roof: RoofSide,
    /// `elems_per_sec / limit_elems_per_sec` (1.0 = at the roof).
    pub attainment: f64,
    /// Fraction of the spec's declared target axis. This equals `attainment`
    /// only when `target_axis == binding_roof` (or the declared physical axis
    /// happens to bind).
    pub target_attainment: f64,
    /// Relative interquartile dispersion of the repetition times
    /// ((p75 − p25) / median): a benchmark without variance bars is
    /// folklore.
    pub dispersion: f64,
    /// Repetitions measured (after warmup).
    pub reps: usize,
    /// Verdict against the declared band.
    pub verdict: Verdict,
    /// Why this row cannot support a verdict. Present exactly when
    /// `verdict == EnvironmentInvalid`.
    pub invalid_reason: Option<String>,
    axis_binding: AxisBinding,
    spec_binding: SpecBinding,
    measurement_origin: MeasurementOrigin,
    pending_tune_publication: Option<fs_session::ValidatedGemmTuneRow>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum MeasurementOrigin {
    Analytic,
    Timed {
        elements: usize,
        warmup_runs: usize,
        sample_seconds_bits: Vec<u64>,
        decision_bindings: Vec<Option<KernelExecutionBinding>>,
    },
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct AxisBinding {
    fingerprint: u64,
    logical_cpus: u32,
    bandwidth_single_bits: u64,
    bandwidth_all_core_bits: u64,
    peak_single_bits: u64,
    peak_all_core_bits: u64,
}

impl AxisBinding {
    fn new(axes: &MachineAxes) -> Self {
        Self {
            fingerprint: axes.fingerprint,
            logical_cpus: axes.logical_cpus,
            bandwidth_single_bits: axes.bandwidth_single_gbs.to_bits(),
            bandwidth_all_core_bits: axes.bandwidth_all_core_gbs.to_bits(),
            peak_single_bits: axes.peak_single_gflops.to_bits(),
            peak_all_core_bits: axes.peak_all_core_gflops.to_bits(),
        }
    }

    fn matches(&self, axes: &MachineAxes) -> bool {
        self == &Self::new(axes)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct SpecBinding {
    kernel: String,
    version: String,
    bytes_per_elem_bits: u64,
    flops_per_elem_bits: u64,
    threading: Threading,
    target_axis: TargetAxis,
    target_fraction_bits: Option<u64>,
}

impl SpecBinding {
    fn new(spec: &KernelSpec) -> Self {
        Self {
            kernel: spec.name.to_string(),
            version: spec.version.to_string(),
            bytes_per_elem_bits: spec.bytes_per_elem.to_bits(),
            flops_per_elem_bits: spec.flops_per_elem.to_bits(),
            threading: spec.threading,
            target_axis: spec.target_axis,
            target_fraction_bits: spec.target_fraction.map(f64::to_bits),
        }
    }
}

fn stable_decision_binding(
    bindings: &[Option<KernelExecutionBinding>],
) -> Option<&KernelExecutionBinding> {
    let first = bindings.first()?.as_ref()?;
    bindings
        .iter()
        .enumerate()
        .all(|(ordinal, candidate)| {
            let Some(candidate) = candidate.as_ref() else {
                return false;
            };
            let Ok(ordinal) = u64::try_from(ordinal) else {
                return false;
            };
            let Some(expected_run) = first.declared_run().checked_add(ordinal) else {
                return false;
            };
            candidate.declared_run() == expected_run && candidate.stable_equivalent(first)
        })
        .then_some(first)
}

impl Attainment {
    /// One JSON line for logs/agents (stable field order).
    #[must_use]
    #[allow(clippy::too_many_lines)] // one canonical receipt schema, kept visibly ordered
    pub fn to_jsonl(&self) -> String {
        let invalid_reason = self.invalid_reason.as_ref().map_or_else(
            || "null".to_string(),
            |reason| format!("\"{}\"", json_escape(reason)),
        );
        let target_bits = self
            .spec_binding
            .target_fraction_bits
            .map_or_else(|| "null".to_string(), |bits| format!("\"{bits:016x}\""));
        let measurement = match &self.measurement_origin {
            MeasurementOrigin::Analytic => "{\"origin\":\"analytic\"}".to_string(),
            MeasurementOrigin::Timed {
                elements,
                warmup_runs,
                sample_seconds_bits,
                decision_bindings,
            } => {
                let sample = stats::sample_from_times(
                    sample_seconds_bits
                        .iter()
                        .copied()
                        .map(f64::from_bits)
                        .collect(),
                );
                let sample_bits = sample_seconds_bits
                    .iter()
                    .map(|bits| format!("\"{bits:016x}\""))
                    .collect::<Vec<_>>()
                    .join(",");
                let decision_bits = decision_bindings
                    .iter()
                    .map(|binding| {
                        binding.as_ref().map_or_else(
                            || "null".to_string(),
                            |binding| format!("\"{}\"", binding.receipt_identity()),
                        )
                    })
                    .collect::<Vec<_>>()
                    .join(",");
                format!(
                    "{{\"origin\":\"timed\",\"elements\":{elements},\"warmup_runs\":{warmup_runs},\"sample_seconds_bits\":[{sample_bits}],\"decision_binding_hashes\":[{decision_bits}],\"median_seconds_bits\":\"{:016x}\",\"p25_seconds_bits\":\"{:016x}\",\"p75_seconds_bits\":\"{:016x}\",\"dispersion_bits\":\"{:016x}\"}}",
                    sample.median.to_bits(),
                    sample.p25.to_bits(),
                    sample.p75.to_bits(),
                    sample.dispersion.to_bits(),
                )
            }
        };
        let execution = match &self.measurement_origin {
            MeasurementOrigin::Timed {
                decision_bindings, ..
            } => stable_decision_binding(decision_bindings).map_or_else(
                || "null".to_string(),
                KernelExecutionBinding::canonical_json,
            ),
            MeasurementOrigin::Analytic => "null".to_string(),
        };
        format!(
            "{{\"receipt_version\":3,\"kernel\":\"{}\",\"version\":\"{}\",\"machine\":\"{:016x}\",\
             \"axes\":{{\"logical_cpus\":{},\"bandwidth_single_bits\":\"{:016x}\",\"bandwidth_all_core_bits\":\"{:016x}\",\"peak_single_bits\":\"{:016x}\",\"peak_all_core_bits\":\"{:016x}\"}},\
             \"spec\":{{\"bytes_per_elem_bits\":\"{:016x}\",\"flops_per_elem_bits\":\"{:016x}\",\"threading\":\"{}\",\"target_axis\":\"{}\",\"target_fraction_bits\":{}}},\"measurement\":{},\"execution\":{},\"elems_per_sec_bits\":\"{:016x}\",\"gbs_bits\":\"{:016x}\",\"gflops_bits\":\"{:016x}\",\"limit_elems_per_sec_bits\":\"{:016x}\",\"attainment_bits\":\"{:016x}\",\"target_attainment_bits\":\"{:016x}\",\"dispersion_bits\":\"{:016x}\",\"elems_per_sec_display\":{:.3e},\
             \"gbs\":{:.3},\"gflops\":{:.3},\"limit_elems_per_sec\":{:.3e},\
             \"roof\":\"{}\",\"attainment\":{:.4},\"target_attainment\":{:.4},\"dispersion\":{:.4},\
             \"reps\":{},\"verdict\":\"{}\",\"invalid_reason\":{}}}",
            json_escape(&self.kernel),
            json_escape(&self.version),
            self.axis_binding.fingerprint,
            self.axis_binding.logical_cpus,
            self.axis_binding.bandwidth_single_bits,
            self.axis_binding.bandwidth_all_core_bits,
            self.axis_binding.peak_single_bits,
            self.axis_binding.peak_all_core_bits,
            self.spec_binding.bytes_per_elem_bits,
            self.spec_binding.flops_per_elem_bits,
            match self.spec_binding.threading {
                Threading::SingleThread => "single_thread",
                Threading::AllCore => "all_core",
            },
            self.spec_binding.target_axis.name(),
            target_bits,
            measurement,
            execution,
            self.elems_per_sec.to_bits(),
            self.achieved_gbs.to_bits(),
            self.achieved_gflops.to_bits(),
            self.limit_elems_per_sec.to_bits(),
            self.attainment.to_bits(),
            self.target_attainment.to_bits(),
            self.dispersion.to_bits(),
            self.elems_per_sec,
            self.achieved_gbs,
            self.achieved_gflops,
            self.limit_elems_per_sec,
            match self.roof {
                RoofSide::Bandwidth => "bandwidth",
                RoofSide::Compute => "compute",
            },
            self.attainment,
            self.target_attainment,
            self.dispersion,
            self.reps,
            self.verdict.name(),
            invalid_reason,
        )
    }

    #[allow(clippy::too_many_lines)] // fail-closed rederivation is easier to audit as one flow
    fn is_citable_against(&self, axes: &MachineAxes) -> bool {
        if !self.axis_binding.matches(axes)
            || self.kernel != self.spec_binding.kernel
            || self.version != self.spec_binding.version
            || self.kernel.trim().is_empty()
            || self.version.trim().is_empty()
            || self.reps == 0
            || self.invalid_reason.is_some()
            || self.verdict == Verdict::EnvironmentInvalid
        {
            return false;
        }
        let MeasurementOrigin::Timed {
            elements,
            warmup_runs,
            sample_seconds_bits,
            decision_bindings,
        } = &self.measurement_origin
        else {
            return false;
        };
        let sample_times: Vec<_> = sample_seconds_bits
            .iter()
            .copied()
            .map(f64::from_bits)
            .collect();
        if *elements == 0
            || sample_times.is_empty()
            || sample_times.len() != self.reps
            || decision_bindings.len() != self.reps
            || sample_times
                .iter()
                .any(|seconds| !seconds.is_finite() || *seconds <= 0.0)
        {
            return false;
        }
        let sample = stats::sample_from_times(sample_times);
        if ((*elements as f64) / sample.median).to_bits() != self.elems_per_sec.to_bits()
            || sample.dispersion.to_bits() != self.dispersion.to_bits()
        {
            return false;
        }
        let bytes_per_elem = f64::from_bits(self.spec_binding.bytes_per_elem_bits);
        let flops_per_elem = f64::from_bits(self.spec_binding.flops_per_elem_bits);
        let target = self.spec_binding.target_fraction_bits.map(f64::from_bits);
        if !bytes_per_elem.is_finite()
            || !flops_per_elem.is_finite()
            || bytes_per_elem < 0.0
            || flops_per_elem < 0.0
            || (bytes_per_elem == 0.0 && flops_per_elem == 0.0)
            || target.is_some_and(|value| !value.is_finite() || value <= 0.0 || value > 1.0)
            || !self.elems_per_sec.is_finite()
            || self.elems_per_sec < 0.0
            || !self.dispersion.is_finite()
            || self.dispersion < 0.0
        {
            return false;
        }
        let (bandwidth_gbs, peak_gflops) = match self.spec_binding.threading {
            Threading::SingleThread => (axes.bandwidth_single_gbs, axes.peak_single_gflops),
            Threading::AllCore => (axes.bandwidth_all_core_gbs, axes.peak_all_core_gflops),
        };
        let bandwidth_limit = if bytes_per_elem > 0.0 {
            bandwidth_gbs * 1e9 / bytes_per_elem
        } else {
            f64::INFINITY
        };
        let compute_limit = if flops_per_elem > 0.0 {
            peak_gflops * 1e9 / flops_per_elem
        } else {
            f64::INFINITY
        };
        let (limit, roof) = if bandwidth_limit <= compute_limit {
            (bandwidth_limit, RoofSide::Bandwidth)
        } else {
            (compute_limit, RoofSide::Compute)
        };
        if !limit.is_finite() || limit <= 0.0 {
            return false;
        }
        let attainment = self.elems_per_sec / limit;
        let achieved_gbs = self.elems_per_sec * bytes_per_elem / 1e9;
        let achieved_gflops = self.elems_per_sec * flops_per_elem / 1e9;
        let target_attainment = match self.spec_binding.target_axis {
            TargetAxis::BindingRoof => attainment,
            TargetAxis::ComputePeak => achieved_gflops / peak_gflops,
            TargetAxis::MemoryBandwidth => achieved_gbs / bandwidth_gbs,
        };
        if !attainment.is_finite()
            || attainment > 1.5
            || !target_attainment.is_finite()
            || !achieved_gbs.is_finite()
            || !achieved_gflops.is_finite()
        {
            return false;
        }
        let expected_verdict = match target {
            None => Verdict::NoTarget,
            Some(fraction) if target_attainment >= fraction => Verdict::WithinBand,
            Some(_) => Verdict::BelowBand,
        };
        let decision_valid = if self.kernel == "gemm-f64" {
            *warmup_runs > 0
                && stable_decision_binding(decision_bindings)
                    .is_some_and(|binding| binding.is_valid_for(axes.fingerprint))
        } else {
            decision_bindings.iter().all(Option::is_none)
        };
        let pending_valid = self
            .pending_tune_publication
            .as_ref()
            .is_none_or(|pending| {
                self.kernel == "gemm-f64"
                    && stable_decision_binding(decision_bindings).is_some_and(|binding| {
                        pending.receipt_identity() == binding.validated_row().receipt_identity()
                    })
            });
        self.roof == roof
            && self.limit_elems_per_sec.to_bits() == limit.to_bits()
            && self.attainment.to_bits() == attainment.to_bits()
            && self.target_attainment.to_bits() == target_attainment.to_bits()
            && self.achieved_gbs.to_bits() == achieved_gbs.to_bits()
            && self.achieved_gflops.to_bits() == achieved_gflops.to_bits()
            && self.verdict == expected_verdict
            && decision_valid
            && pending_valid
    }
}

fn json_escape(value: &str) -> String {
    let mut out = String::with_capacity(value.len());
    for ch in value.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if c.is_control() => {
                use std::fmt::Write as _;
                let _ = write!(out, "\\u{:04x}", u32::from(c));
            }
            c => out.push(c),
        }
    }
    out
}

fn spec_error(spec: &KernelSpec) -> Option<&'static str> {
    if spec.name.trim().is_empty() || spec.version.trim().is_empty() {
        return Some("kernel name and version must be non-empty");
    }
    if !spec.bytes_per_elem.is_finite()
        || !spec.flops_per_elem.is_finite()
        || spec.bytes_per_elem < 0.0
        || spec.flops_per_elem < 0.0
        || (spec.bytes_per_elem == 0.0 && spec.flops_per_elem == 0.0)
    {
        return Some("kernel intensity must be finite, non-negative, and exercise an axis");
    }
    if let Some(target) = spec.target_fraction
        && (!target.is_finite() || target <= 0.0 || target > 1.0)
    {
        return Some("target fraction must be finite and in (0, 1]");
    }
    if spec.target_fraction.is_some()
        && match spec.target_axis {
            TargetAxis::BindingRoof => false,
            TargetAxis::ComputePeak => spec.flops_per_elem == 0.0,
            TargetAxis::MemoryBandwidth => spec.bytes_per_elem == 0.0,
        }
    {
        return Some("declared target axis is not exercised by this kernel");
    }
    None
}

/// Compute attainment from a measured rate and the machine axes. Pure
/// arithmetic — meta-tested against hand calculations (`rf_002`).
#[must_use]
pub fn attainment_for(spec: &KernelSpec, elems_per_sec: f64, axes: &MachineAxes) -> Attainment {
    attainment_with_dispersion(spec, elems_per_sec, 0.0, 0, axes)
}

/// [`attainment_for`] with measured dispersion and repetition count.
#[must_use]
pub fn attainment_with_dispersion(
    spec: &KernelSpec,
    elems_per_sec: f64,
    dispersion: f64,
    reps: usize,
    axes: &MachineAxes,
) -> Attainment {
    attainment_with_origin(
        spec,
        elems_per_sec,
        dispersion,
        reps,
        axes,
        MeasurementOrigin::Analytic,
    )
}

#[allow(clippy::too_many_lines)] // one roof decision per branch; splitting obscures the min-roof flow
fn attainment_with_origin(
    spec: &KernelSpec,
    elems_per_sec: f64,
    dispersion: f64,
    reps: usize,
    axes: &MachineAxes,
    measurement_origin: MeasurementOrigin,
) -> Attainment {
    let spec_error = spec_error(spec);
    let measurement_error = if !elems_per_sec.is_finite() || elems_per_sec < 0.0 {
        Some("measured element rate is non-finite or negative")
    } else if !dispersion.is_finite() || dispersion < 0.0 {
        Some("measured dispersion is non-finite or negative")
    } else {
        None
    };
    let safe_rate = if measurement_error.is_none() {
        elems_per_sec
    } else {
        0.0
    };
    let safe_bytes = if spec.bytes_per_elem.is_finite() && spec.bytes_per_elem >= 0.0 {
        spec.bytes_per_elem
    } else {
        0.0
    };
    let safe_flops = if spec.flops_per_elem.is_finite() && spec.flops_per_elem >= 0.0 {
        spec.flops_per_elem
    } else {
        0.0
    };
    let (bandwidth_gbs, peak_gflops) = match spec.threading {
        Threading::SingleThread => (axes.bandwidth_single_gbs, axes.peak_single_gflops),
        Threading::AllCore => (axes.bandwidth_all_core_gbs, axes.peak_all_core_gflops),
    };
    // Limits in elements/second on each axis; +inf when the kernel does not
    // exercise an axis (zero bytes or zero flops per element).
    let bw_limit = if safe_bytes > 0.0 && bandwidth_gbs.is_finite() && bandwidth_gbs > 0.0 {
        bandwidth_gbs * 1e9 / safe_bytes
    } else {
        f64::INFINITY
    };
    let comp_limit = if safe_flops > 0.0 && peak_gflops.is_finite() && peak_gflops > 0.0 {
        peak_gflops * 1e9 / safe_flops
    } else {
        f64::INFINITY
    };
    let (limit, roof) = if bw_limit <= comp_limit {
        (bw_limit, RoofSide::Bandwidth)
    } else {
        (comp_limit, RoofSide::Compute)
    };
    let raw_attainment = if limit.is_finite() && limit > 0.0 {
        safe_rate / limit
    } else {
        0.0
    };
    let raw_achieved_gbs = safe_rate * safe_bytes / 1e9;
    let raw_achieved_gflops = safe_rate * safe_flops / 1e9;
    let raw_target_attainment = match spec.target_axis {
        TargetAxis::BindingRoof => raw_attainment,
        TargetAxis::ComputePeak if peak_gflops.is_finite() && peak_gflops > 0.0 => {
            raw_achieved_gflops / peak_gflops
        }
        TargetAxis::MemoryBandwidth if bandwidth_gbs.is_finite() && bandwidth_gbs > 0.0 => {
            raw_achieved_gbs / bandwidth_gbs
        }
        TargetAxis::ComputePeak | TargetAxis::MemoryBandwidth => f64::NAN,
    };
    let derived_error = if !raw_attainment.is_finite()
        || !raw_target_attainment.is_finite()
        || !raw_achieved_gbs.is_finite()
        || !raw_achieved_gflops.is_finite()
    {
        Some("derived roofline quantities overflowed or became non-finite")
    } else {
        None
    };
    let attainment = if raw_attainment.is_finite() {
        raw_attainment
    } else {
        0.0
    };
    let target_attainment = if raw_target_attainment.is_finite() {
        raw_target_attainment
    } else {
        0.0
    };
    // Environment validity BEFORE band comparison (bead 1n61): the
    // ratio is meaningless when the axes are implausible, and an
    // attainment materially above 1 means the kernel outran its own
    // roofline — the axes were probed under different (crushed)
    // conditions than the kernel run. Refuse to gate either way.
    let invalid_reason = axes
        .plausibility_error()
        .or(spec_error)
        .or(measurement_error)
        .or(derived_error)
        .or_else(|| (attainment > 1.5).then_some("attainment exceeds the credible roofline band"));
    let verdict = if invalid_reason.is_some() {
        Verdict::EnvironmentInvalid
    } else {
        match spec.target_fraction {
            None => Verdict::NoTarget,
            Some(t) if target_attainment >= t => Verdict::WithinBand,
            Some(_) => Verdict::BelowBand,
        }
    };
    Attainment {
        kernel: spec.name.to_string(),
        version: spec.version.to_string(),
        elems_per_sec: safe_rate,
        achieved_gbs: if raw_achieved_gbs.is_finite() {
            raw_achieved_gbs
        } else {
            0.0
        },
        achieved_gflops: if raw_achieved_gflops.is_finite() {
            raw_achieved_gflops
        } else {
            0.0
        },
        limit_elems_per_sec: if limit.is_finite() { limit } else { 0.0 },
        roof,
        attainment,
        target_attainment,
        dispersion: if dispersion.is_finite() && dispersion >= 0.0 {
            dispersion
        } else {
            0.0
        },
        reps,
        verdict,
        invalid_reason: invalid_reason.map(str::to_string),
        axis_binding: AxisBinding::new(axes),
        spec_binding: SpecBinding::new(spec),
        measurement_origin,
        pending_tune_publication: None,
    }
}

/// Measure one kernel (warmup + repetitions) and compute its attainment.
pub fn measure(
    kernel: &mut dyn RooflineKernel,
    warmup: usize,
    reps: usize,
    axes: &MachineAxes,
) -> Attainment {
    let spec = kernel.spec();
    let elements = kernel.elements();
    let elems = elements as f64;
    for _ in 0..warmup {
        kernel.run_once();
    }
    let measured_reps = reps.max(1);
    let mut times = Vec::with_capacity(measured_reps);
    let mut decision_bindings = Vec::with_capacity(measured_reps);
    for _ in 0..measured_reps {
        let start = std::time::Instant::now();
        kernel.run_once();
        times.push(start.elapsed().as_secs_f64());
        decision_bindings.push(kernel.execution_binding());
    }
    let sample = stats::sample_from_times(times);
    let elems_per_sec = if sample.median > 0.0 {
        elems / sample.median
    } else {
        0.0
    };
    let mut result = attainment_with_origin(
        &spec,
        elems_per_sec,
        sample.dispersion,
        measured_reps,
        axes,
        MeasurementOrigin::Timed {
            elements,
            warmup_runs: warmup,
            sample_seconds_bits: sample.times.iter().map(|value| value.to_bits()).collect(),
            decision_bindings,
        },
    );
    result.pending_tune_publication = kernel.pending_tune_publication();
    result
}

/// Run every kernel in the registry.
pub fn run_registry(
    registry: &mut [Box<dyn RooflineKernel>],
    warmup: usize,
    reps: usize,
    axes: &MachineAxes,
) -> Vec<Attainment> {
    let mut results: Vec<_> = registry
        .iter_mut()
        .map(|k| measure(k.as_mut(), warmup, reps, axes))
        .collect();
    poison_invalid_run(&mut results);
    results
}

/// Whether a registry result set is admissible as citable performance
/// evidence for these exact measured axes.
#[must_use]
pub fn run_is_citable(
    axes: &MachineAxes,
    post_axes: &MachineAxes,
    baseline: AxisBaselinePolicy<'_>,
    results: &[Attainment],
) -> bool {
    run_admission_error(axes, post_axes, baseline, results).is_none()
}

fn finalized_run_receipt(
    axes: &MachineAxes,
    post_axes: &MachineAxes,
    baseline: AxisBaselinePolicy<'_>,
    results: &[Attainment],
) -> fs_blake3::ContentHash {
    fn push_field(payload: &mut Vec<u8>, bytes: &[u8]) {
        let len = u64::try_from(bytes.len()).expect("receipt field length fits u64");
        payload.extend_from_slice(&len.to_le_bytes());
        payload.extend_from_slice(bytes);
    }

    let mut payload = Vec::new();
    let baseline_receipt = baseline.receipt_json(axes, post_axes);
    push_field(&mut payload, baseline_receipt.as_bytes());
    let result_count = u64::try_from(results.len()).expect("result count fits u64");
    payload.extend_from_slice(&result_count.to_le_bytes());
    for result in results {
        let receipt = result.to_jsonl();
        push_field(&mut payload, receipt.as_bytes());
    }
    fs_blake3::hash_domain(FINALIZED_RUN_DOMAIN, &payload)
}

/// Complete the registry's tuning lifecycle at the same admission boundary
/// used for citable performance evidence.
///
/// Every kernel sees the single aggregate decision. Publication belongs to
/// [`record_run`]'s evidence transaction; this hook clears kernel-owned pending
/// markers and invalidates rejected local decisions.
///
/// # Errors
/// Calls every kernel, then returns all lifecycle diagnostics in deterministic
/// registry order.
pub fn finalize_registry_tuning(
    registry: &mut [Box<dyn RooflineKernel>],
    axes: &MachineAxes,
    post_axes: &MachineAxes,
    baseline: AxisBaselinePolicy<'_>,
    results: &[Attainment],
) -> Result<FinalizedRegistryRun, String> {
    let mut diagnostics = Vec::new();
    let mut registry_matches = registry.len() == results.len();
    if !registry_matches {
        diagnostics.push(format!(
            "registry/result length mismatch: {} kernels for {} result rows",
            registry.len(),
            results.len()
        ));
    }
    for (index, (kernel, result)) in registry.iter().zip(results).enumerate() {
        let spec = kernel.spec();
        if spec.name != result.kernel || spec.version != result.version {
            registry_matches = false;
            diagnostics.push(format!(
                "kernel[{index}] identity mismatch: registry {}/{} vs result {}/{}",
                spec.name, spec.version, result.kernel, result.version
            ));
        }
        let expected_binding = match &result.measurement_origin {
            MeasurementOrigin::Timed {
                decision_bindings, ..
            } => decision_bindings.last().cloned().flatten(),
            MeasurementOrigin::Analytic => None,
        };
        if kernel.execution_binding() != expected_binding {
            registry_matches = false;
            diagnostics.push(format!(
                "kernel[{index}] execution state changed after this result was measured"
            ));
        }
        if kernel.pending_tune_publication() != result.pending_tune_publication {
            registry_matches = false;
            diagnostics.push(format!(
                "kernel[{index}] pending tune publication changed after this result was measured"
            ));
        }
    }
    let admitted =
        registry_matches && run_admission_error(axes, post_axes, baseline, results).is_none();
    for (index, kernel) in registry.iter_mut().enumerate() {
        if let Err(error) = kernel.finalize_tuning(admitted) {
            diagnostics.push(format!("kernel[{index}]: {error}"));
        }
    }
    if !diagnostics.is_empty() {
        return Err(format!(
            "tuning lifecycle finalization failed for {} kernel(s): {}",
            diagnostics.len(),
            diagnostics.join("; ")
        ));
    }
    Ok(FinalizedRegistryRun {
        receipt: finalized_run_receipt(axes, post_axes, baseline, results),
        admitted,
        consumed: false,
    })
}

/// Explain why a result set cannot be admitted as citable performance
/// evidence. `None` means every row is a bound timed measurement.
#[must_use]
pub fn run_admission_error(
    axes: &MachineAxes,
    post_axes: &MachineAxes,
    baseline: AxisBaselinePolicy<'_>,
    results: &[Attainment],
) -> Option<String> {
    if results.is_empty() {
        return Some("registry produced no measured kernels".to_string());
    }
    let baseline_verdict = baseline.verdict(axes, post_axes);
    if !baseline_verdict.trusted() {
        return Some(format!(
            "historical baseline admission refused: {}",
            baseline_verdict.to_jsonl()
        ));
    }
    let mut identities = std::collections::BTreeSet::new();
    for (index, result) in results.iter().enumerate() {
        if !identities.insert(result.kernel.as_str()) {
            return Some(format!(
                "duplicate kernel identity at row {index}: {}/{}",
                result.kernel, result.version
            ));
        }
        if !result.is_citable_against(axes) {
            return Some(format!(
                "row {index} ({}/{}) is not a bound timed measurement for these axes",
                result.kernel, result.version
            ));
        }
    }
    None
}

fn poison_invalid_run(results: &mut [Attainment]) {
    let Some(origin) = results
        .iter()
        .find(|result| result.verdict == Verdict::EnvironmentInvalid)
    else {
        return;
    };
    let reason = format!(
        "registry invalidated by {}: {}",
        origin.kernel,
        origin
            .invalid_reason
            .as_deref()
            .unwrap_or("invalid evidence")
    );
    for result in results {
        result.verdict = Verdict::EnvironmentInvalid;
        result.invalid_reason = Some(reason.clone());
    }
}

// ---------------------------------------------------------------------------
// §14.1 target table as data
// ---------------------------------------------------------------------------

/// One row of the plan §14.1 target table. `landed = false` rows are
/// visible from day one so nothing is silently uncovered — they flip as the
/// owning kernels register.
#[derive(Debug, Clone, Copy)]
pub struct TargetRow {
    /// Kernel family name.
    pub kernel: &'static str,
    /// What the target means (unit and roof context).
    pub statement: &'static str,
    /// Whether an implementation is registered in this harness yet.
    pub landed: bool,
}

/// The §14.1 table. Statements, not claims: every value must be re-measured
/// on a fingerprinted machine before anyone may cite it.
pub const SECTION_14_1_TARGETS: &[TargetRow] = &[
    TargetRow {
        kernel: "lbm-d3q19-stream-collide",
        statement: "≥1.0 GLUP/s (M-class) / ≥0.6 GLUP/s (TR-class), bandwidth-bound",
        landed: false,
    },
    TargetRow {
        kernel: "gemm-f64",
        statement: "≥75% of measured peak FLOPs for the selected SIMD tier",
        landed: true,
    },
    TargetRow {
        kernel: "spmv-sell-c-sigma",
        statement: "≥85% of measured STREAM-class bandwidth",
        landed: false,
    },
    TargetRow {
        kernel: "feec-apply-p4",
        statement: "≥30% of peak FLOPs, sum-factorized",
        landed: false,
    },
    TargetRow {
        kernel: "batched-small-dense",
        statement: "≥60% of peak FLOPs, SIMD-across-elements",
        landed: false,
    },
    TargetRow {
        kernel: "fft-3d-pencil",
        statement: "≥40% of the memory-bound limit",
        landed: false,
    },
    TargetRow {
        kernel: "sdf-primary-rays",
        statement: "≥80 Mray/s (M-class) / ≥120 Mray/s (TR-class)",
        landed: false,
    },
];

const ROOFLINE_ROW_SCHEMA: &str = "fs-roofline-ledger-row-v2";
const ROOFLINE_PAYLOAD_ARTIFACT_KIND: &str = "roofline-benchmark-result";
const ROOFLINE_EXECUTABLE_DOMAIN: &str = "org.frankensim.fs-roofline.executable.v1";

/// Maximum age of a roofline measurement that can be reported as fresh.
pub const STALENESS_MAX_AGE_NS: i64 = 30 * 24 * 60 * 60 * 1_000_000_000;

fn versions_json(build_identity: fs_blake3::ContentHash) -> String {
    format!("{{\"frankensim_executable\":\"{build_identity}\",\"fs-roofline\":\"{VERSION}\"}}")
}

fn executable_build_identity() -> Result<fs_blake3::ContentHash, LedgerError> {
    static IDENTITY: std::sync::OnceLock<Result<fs_blake3::ContentHash, LedgerError>> =
        std::sync::OnceLock::new();
    IDENTITY.get_or_init(read_executable_build_identity).clone()
}

fn read_executable_build_identity() -> Result<fs_blake3::ContentHash, LedgerError> {
    use std::io::Read as _;

    let path = std::env::current_exe().map_err(|error| LedgerError::Invalid {
        field: "executable_identity".to_string(),
        problem: format!("cannot resolve current executable: {error}"),
    })?;
    let mut file = std::fs::File::open(&path).map_err(|error| LedgerError::Invalid {
        field: "executable_identity".to_string(),
        problem: format!("cannot open current executable {}: {error}", path.display()),
    })?;
    let mut hasher = fs_blake3::Blake3::new();
    let mut total = 0_u64;
    let mut chunk = [0_u8; 64 * 1024];
    loop {
        let read = file
            .read(&mut chunk)
            .map_err(|error| LedgerError::Invalid {
                field: "executable_identity".to_string(),
                problem: format!("cannot read current executable {}: {error}", path.display()),
            })?;
        if read == 0 {
            break;
        }
        total = total
            .checked_add(u64::try_from(read).expect("read chunk length fits u64"))
            .ok_or_else(|| LedgerError::Invalid {
                field: "executable_identity".to_string(),
                problem: "current executable length exceeds u64".to_string(),
            })?;
        hasher.update(&chunk[..read]);
    }
    let raw = hasher.finalize();
    let mut preimage = [0_u8; 40];
    preimage[..8].copy_from_slice(&total.to_le_bytes());
    preimage[8..].copy_from_slice(raw.as_bytes());
    Ok(fs_blake3::hash_domain(
        ROOFLINE_EXECUTABLE_DOMAIN,
        &preimage,
    ))
}

#[derive(Debug)]
struct RooflineRowParams {
    op: i64,
    run_receipt: fs_blake3::ContentHash,
    payload_artifact: fs_blake3::ContentHash,
    baseline_hash: fs_blake3::ContentHash,
    build_identity: fs_blake3::ContentHash,
    reps: u64,
    post_axis_bits: [u64; 4],
}

impl RooflineRowParams {
    fn to_json(&self) -> String {
        format!(
            "{{\"schema\":\"{ROOFLINE_ROW_SCHEMA}\",\"op\":{},\"run_receipt\":\"{}\",\"payload_artifact\":\"{}\",\"baseline_hash\":\"{}\",\"build_identity\":\"{}\",\"reps\":{},\"post_bandwidth_single_bits\":\"{:016x}\",\"post_bandwidth_all_core_bits\":\"{:016x}\",\"post_peak_single_bits\":\"{:016x}\",\"post_peak_all_core_bits\":\"{:016x}\"}}",
            self.op,
            self.run_receipt,
            self.payload_artifact,
            self.baseline_hash,
            self.build_identity,
            self.reps,
            self.post_axis_bits[0],
            self.post_axis_bits[1],
            self.post_axis_bits[2],
            self.post_axis_bits[3],
        )
    }
}

fn parse_roofline_row_params(text: &str) -> Option<RooflineRowParams> {
    fn take<'a>(rest: &mut &'a str, prefix: &str) -> Option<()> {
        *rest = rest.strip_prefix(prefix)?;
        Some(())
    }
    fn decimal(rest: &mut &str) -> Option<u64> {
        let end = rest
            .find(|ch: char| !ch.is_ascii_digit())
            .unwrap_or(rest.len());
        if end == 0 || (end > 1 && rest.starts_with('0')) {
            return None;
        }
        let (digits, tail) = rest.split_at(end);
        *rest = tail;
        digits.parse().ok()
    }
    fn hash(rest: &mut &str) -> Option<fs_blake3::ContentHash> {
        let (raw, tail) = rest.split_once('"')?;
        if raw.len() != 64 {
            return None;
        }
        *rest = tail;
        fs_blake3::ContentHash::from_hex(raw)
    }
    fn bits(rest: &mut &str) -> Option<u64> {
        let (raw, tail) = rest.split_once('"')?;
        if raw.len() != 16
            || !raw
                .bytes()
                .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
        {
            return None;
        }
        *rest = tail;
        u64::from_str_radix(raw, 16).ok()
    }

    let mut rest = text;
    take(
        &mut rest,
        &format!("{{\"schema\":\"{ROOFLINE_ROW_SCHEMA}\",\"op\":"),
    )?;
    let op = i64::try_from(decimal(&mut rest)?).ok()?;
    if op <= 0 {
        return None;
    }
    take(&mut rest, ",\"run_receipt\":\"")?;
    let run_receipt = hash(&mut rest)?;
    take(&mut rest, ",\"payload_artifact\":\"")?;
    let payload_artifact = hash(&mut rest)?;
    take(&mut rest, ",\"baseline_hash\":\"")?;
    let baseline_hash = hash(&mut rest)?;
    take(&mut rest, ",\"build_identity\":\"")?;
    let build_identity = hash(&mut rest)?;
    take(&mut rest, ",\"reps\":")?;
    let reps = decimal(&mut rest)?;
    if reps == 0 {
        return None;
    }
    take(&mut rest, ",\"post_bandwidth_single_bits\":\"")?;
    let bandwidth_single = bits(&mut rest)?;
    take(&mut rest, ",\"post_bandwidth_all_core_bits\":\"")?;
    let bandwidth_all_core = bits(&mut rest)?;
    take(&mut rest, ",\"post_peak_single_bits\":\"")?;
    let peak_single = bits(&mut rest)?;
    take(&mut rest, ",\"post_peak_all_core_bits\":\"")?;
    let peak_all_core = bits(&mut rest)?;
    take(&mut rest, "}")?;
    if !rest.is_empty() {
        return None;
    }
    let params = RooflineRowParams {
        op,
        run_receipt,
        payload_artifact,
        baseline_hash,
        build_identity,
        reps,
        post_axis_bits: [
            bandwidth_single,
            bandwidth_all_core,
            peak_single,
            peak_all_core,
        ],
    };
    (params.to_json() == text).then_some(params)
}

// ---------------------------------------------------------------------------
// Ledger integration and staleness
// ---------------------------------------------------------------------------

/// Record a harness run atomically in the ledger. Admitted timed rows receive
/// metrics, `benchmark_result` events, content-addressed output artifacts, and
/// versioned tune rows. Rejected input receives an exact baseline-admission
/// receipt plus one explicit rejection event and no normal-looking
/// measurements.
/// A successful commit consumes each result's one-shot fresh-row marker;
/// rollback retains it so the same transaction can be retried.
///
/// # Errors
/// Ledger errors propagate and roll back the whole write set.
#[allow(clippy::too_many_lines)] // one auditable all-or-nothing evidence transaction
pub fn record_run(
    ledger: &Ledger,
    axes: &MachineAxes,
    post_axes: &MachineAxes,
    baseline: AxisBaselinePolicy<'_>,
    finalized: &mut FinalizedRegistryRun,
    results: &mut [Attainment],
) -> Result<i64, LedgerError> {
    let admission_error = run_admission_error(axes, post_axes, baseline, results);
    let run_valid = admission_error.is_none();
    if finalized.consumed {
        return Err(LedgerError::Invalid {
            field: "finalized_run".to_string(),
            problem: "the finalized roofline run was already recorded".to_string(),
        });
    }
    let expected_receipt = finalized_run_receipt(axes, post_axes, baseline, results);
    if finalized.receipt != expected_receipt || finalized.admitted != run_valid {
        return Err(LedgerError::Invalid {
            field: "finalized_run".to_string(),
            problem:
                "axes, baseline decision, results, or admission changed after registry finalization"
                    .to_string(),
        });
    }
    let baseline_receipt = baseline.receipt_json(axes, post_axes);
    let build_identity = executable_build_identity()?;
    let versions = versions_json(build_identity);
    let explicits = FiveExplicits {
        seed: b"roofline",
        versions: &versions,
        budget: "{\"wall_s\":600}",
        capability: "{\"ops\":[\"perf.roofline\"]}",
    };
    let ir = format!(
        "{{\"op\":\"perf.roofline\",\"kernels\":{},\"fingerprint\":\"{:016x}\",\"post_fingerprint\":\"{:016x}\",\"admitted\":{run_valid},\"finalized_run_receipt\":\"{}\",\"baseline_admission\":{baseline_receipt}}}",
        results.len(),
        axes.fingerprint,
        post_axes.fingerprint,
        finalized.receipt,
    );
    ledger.begin()?;
    let write_result: Result<i64, LedgerError> = (|| {
        let op = ledger.begin_op(Some(b"roofline"), &ir, &explicits, now_wall_ns())?;
        ledger.append_event(&EventRow {
            session: Some(b"roofline"),
            t: 0,
            kind: "axis_baseline_admission",
            payload: Some(&baseline_receipt),
        })?;
        if run_valid {
            let baseline_hash = baseline
                .baseline_hash()
                .ok_or_else(|| LedgerError::Invalid {
                    field: "baseline".to_string(),
                    problem: "trusted roofline admission has no selected baseline".to_string(),
                })?;
            let machine_key = roofline_machine_key(axes.fingerprint, baseline_hash);
            for r in results.iter() {
                if let MeasurementOrigin::Timed {
                    decision_bindings, ..
                } = &r.measurement_origin
                    && let Some(binding) = stable_decision_binding(decision_bindings)
                {
                    // Publish the exact row from the sealed binding even when
                    // this measurement adopted it from another ledger. This
                    // evidence transaction is deliberately insert-only: a
                    // clone or delayed receipt must never replace a newer
                    // cache row. Refresh/replacement is a separate explicit
                    // protocol, never authority carried by Attainment.
                    binding
                        .validated_row()
                        .publish_if_absent_or_identical(ledger)?;
                }
                ledger.record_metric(
                    op,
                    0,
                    &format!("{}.elems_per_sec", r.kernel),
                    r.elems_per_sec,
                )?;
                ledger.record_metric(op, 0, &format!("{}.attainment", r.kernel), r.attainment)?;
                ledger.record_metric(op, 0, &format!("{}.dispersion", r.kernel), r.dispersion)?;
                let payload = r.to_jsonl();
                ledger.append_event(&EventRow {
                    session: Some(b"roofline"),
                    t: 0,
                    kind: "benchmark_result",
                    payload: Some(&payload),
                })?;
                let payload_artifact = ledger.put_artifact(
                    ROOFLINE_PAYLOAD_ARTIFACT_KIND,
                    payload.as_bytes(),
                    Some("{\"schema\":\"fs-roofline-benchmark-result-v1\"}"),
                )?;
                ledger.link(op, &payload_artifact.hash, EdgeRole::Out)?;
                let params = RooflineRowParams {
                    op,
                    run_receipt: finalized.receipt,
                    payload_artifact: payload_artifact.hash,
                    baseline_hash,
                    build_identity,
                    reps: u64::try_from(r.reps).map_err(|_| LedgerError::Invalid {
                        field: "reps".to_string(),
                        problem: "roofline repetition count exceeds u64".to_string(),
                    })?,
                    post_axis_bits: [
                        post_axes.bandwidth_single_gbs.to_bits(),
                        post_axes.bandwidth_all_core_gbs.to_bits(),
                        post_axes.peak_single_gflops.to_bits(),
                        post_axes.peak_all_core_gflops.to_bits(),
                    ],
                }
                .to_json();
                let shape_class = tune_measurement_shape_class(&r.version, finalized.receipt, op);
                ledger.tune_put_if_absent(
                    &r.kernel,
                    &shape_class,
                    &machine_key,
                    &params,
                    &payload,
                )?;
                let stored = ledger
                    .tune_get(&r.kernel, &shape_class, &machine_key)?
                    .ok_or_else(|| LedgerError::Invalid {
                        field: "tune".to_string(),
                        problem: "roofline insert-if-absent returned without a stored row"
                            .to_string(),
                    })?;
                if stored.params != params || stored.measured != payload {
                    return Err(LedgerError::Invalid {
                        field: "tune".to_string(),
                        problem: format!(
                            "refusing conflicting roofline evidence for kernel {:?}, shape {:?}",
                            r.kernel, shape_class
                        ),
                    });
                }
            }
            ledger.finish_op(op, OpOutcome::Ok, None, now_wall_ns())?;
        } else {
            let reason = admission_error
                .as_deref()
                .unwrap_or("unknown admission failure");
            let payload = format!(
                "{{\"code\":\"roofline_evidence_rejected\",\"reason\":\"{}\",\"effect\":\"no_measurements_or_tune_rows_published\"}}",
                json_escape(reason)
            );
            ledger.append_event(&EventRow {
                session: Some(b"roofline"),
                t: 0,
                kind: "roofline_run_rejected",
                payload: Some(&payload),
            })?;
            let diagnostic = format!(
                "{{\"code\":\"roofline_evidence_rejected\",\"reason\":\"{}\"}}",
                json_escape(reason)
            );
            ledger.finish_op(op, OpOutcome::Error, Some(&diagnostic), now_wall_ns())?;
        }
        Ok(op)
    })();
    match write_result {
        Ok(op) => match ledger.commit() {
            Ok(()) => {
                finalized.consumed = true;
                for result in results.iter_mut() {
                    result.pending_tune_publication = None;
                }
                Ok(op)
            }
            Err(error) => {
                let _ = ledger.rollback();
                Err(error)
            }
        },
        Err(error) => {
            let _ = ledger.rollback();
            Err(error)
        }
    }
}

/// Staleness state of one kernel's ledgered attainment on this machine.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Staleness {
    /// At least one semantically valid row from this exact executable is no
    /// more than [`STALENESS_MAX_AGE_NS`] old.
    Fresh,
    /// Matching current-build evidence exists, but its newest successful
    /// operation is older than [`STALENESS_MAX_AGE_NS`].
    Expired,
    /// The observation clock predates the newest matching operation receipt.
    ClockRollback,
    /// Rows exist, but none for the current fingerprint — the machine
    /// drifted and every cited number is stale until re-measured.
    FingerprintDrift,
    /// The current machine fingerprint matches historical rows, but no trusted
    /// baseline was selected for comparison.
    BaselineUnavailable,
    /// The machine matches but the selected historical baseline changed.
    BaselineDrift,
    /// Semantically valid rows exist for this machine and baseline, but only
    /// for a different executable-content identity.
    BuildDrift,
    /// A row exists under the exact current version, machine, and baseline key,
    /// but its canonical parameters, payload artifact/edge, operation receipt,
    /// or executable identity no longer agree. Corrupt evidence is never
    /// treated as fresh.
    CorruptEvidence,
    /// No roofline rows at all: never measured.
    NeverMeasured,
}

/// Check one kernel's staleness against the current fingerprint and admitted
/// historical baseline identity.
///
/// # Errors
/// Ledger errors propagate.
pub fn staleness(
    ledger: &Ledger,
    kernel: &str,
    version: &str,
    current_fingerprint: u64,
    current_baseline: Option<fs_blake3::ContentHash>,
) -> Result<Staleness, LedgerError> {
    staleness_at(
        ledger,
        kernel,
        version,
        current_fingerprint,
        current_baseline,
        now_wall_ns(),
    )
}

/// Deterministic form of [`staleness`] evaluated at an explicit wall-clock
/// nanosecond. Supplying time makes expiry and rollback tests replayable.
///
/// # Errors
/// Ledger and executable-identity errors propagate.
pub fn staleness_at(
    ledger: &Ledger,
    kernel: &str,
    version: &str,
    current_fingerprint: u64,
    current_baseline: Option<fs_blake3::ContentHash>,
    observed_wall_ns: i64,
) -> Result<Staleness, LedgerError> {
    let current_build = executable_build_identity()?;
    staleness_at_with_build(
        ledger,
        kernel,
        version,
        current_fingerprint,
        current_baseline,
        observed_wall_ns,
        current_build,
    )
}

fn staleness_at_with_build(
    ledger: &Ledger,
    kernel: &str,
    version: &str,
    current_fingerprint: u64,
    current_baseline: Option<fs_blake3::ContentHash>,
    observed_wall_ns: i64,
    current_build: fs_blake3::ContentHash,
) -> Result<Staleness, LedgerError> {
    let rows = ledger.tune_rows(kernel)?;
    let shape_prefix = format!("{}:run=", tune_shape_class(version));
    let roofline_rows: Vec<_> = rows
        .iter()
        .filter(|r| r.shape_class.starts_with(&shape_prefix))
        .collect();
    if roofline_rows.is_empty() {
        return Ok(Staleness::NeverMeasured);
    }
    let fp = current_fingerprint.to_le_bytes();
    let same_machine = roofline_rows
        .iter()
        .filter(|row| row.machine.get(..8) == Some(fp.as_slice()))
        .collect::<Vec<_>>();
    if same_machine.is_empty() {
        return Ok(Staleness::FingerprintDrift);
    }
    let Some(current_baseline) = current_baseline else {
        return Ok(Staleness::BaselineUnavailable);
    };
    let key = roofline_machine_key(current_fingerprint, current_baseline);
    let matching_rows = same_machine
        .into_iter()
        .filter(|row| row.machine == key)
        .collect::<Vec<_>>();
    if matching_rows.is_empty() {
        return Ok(Staleness::BaselineDrift);
    }
    let mut newest_current_build = None;
    let mut saw_foreign_build = false;
    for row in matching_rows {
        let Some(validated) = validate_roofline_row(
            ledger,
            row,
            kernel,
            version,
            current_fingerprint,
            current_baseline,
        )?
        else {
            return Ok(Staleness::CorruptEvidence);
        };
        if validated.build_identity == current_build {
            newest_current_build = Some(
                newest_current_build.map_or(validated.recorded_at_ns, |newest: i64| {
                    newest.max(validated.recorded_at_ns)
                }),
            );
        } else {
            saw_foreign_build = true;
        }
    }
    let Some(recorded_at_ns) = newest_current_build else {
        return Ok(if saw_foreign_build {
            Staleness::BuildDrift
        } else {
            Staleness::CorruptEvidence
        });
    };
    if observed_wall_ns < recorded_at_ns {
        return Ok(Staleness::ClockRollback);
    }
    if observed_wall_ns.saturating_sub(recorded_at_ns) > STALENESS_MAX_AGE_NS {
        return Ok(Staleness::Expired);
    }
    Ok(Staleness::Fresh)
}

struct ValidatedRooflineRow {
    build_identity: fs_blake3::ContentHash,
    recorded_at_ns: i64,
}

fn validate_roofline_row(
    ledger: &Ledger,
    row: &fs_ledger::TuneRow,
    kernel: &str,
    version: &str,
    current_fingerprint: u64,
    current_baseline: fs_blake3::ContentHash,
) -> Result<Option<ValidatedRooflineRow>, LedgerError> {
    let Some(params) = parse_roofline_row_params(&row.params) else {
        return Ok(None);
    };
    if params.baseline_hash != current_baseline
        || row.machine != roofline_machine_key(current_fingerprint, current_baseline)
        || row.shape_class != tune_measurement_shape_class(version, params.run_receipt, params.op)
        || params.payload_artifact != fs_ledger::hash_bytes(row.measured.as_bytes())
    {
        return Ok(None);
    }
    let Some(artifact_bytes) = ledger.get_artifact(&params.payload_artifact)? else {
        return Ok(None);
    };
    if artifact_bytes.as_slice() != row.measured.as_bytes()
        || !ledger.edge_exists(params.op, &params.payload_artifact, EdgeRole::Out)?
    {
        return Ok(None);
    }

    let measured_prefix = format!(
        "{{\"receipt_version\":3,\"kernel\":\"{}\",\"version\":\"{}\",\"machine\":\"{current_fingerprint:016x}\",",
        json_escape(kernel),
        json_escape(version),
    );
    if !row.measured.starts_with(&measured_prefix)
        || !row
            .measured
            .contains(&format!("\"reps\":{},\"verdict\":", params.reps))
    {
        return Ok(None);
    }

    let Some(op) = ledger.op(params.op)? else {
        return Ok(None);
    };
    let Some(recorded_at_ns) = op.t_end else {
        return Ok(None);
    };
    if op.id != params.op
        || op.session.as_deref() != Some(b"roofline".as_slice())
        || op.outcome.as_deref() != Some("ok")
        || recorded_at_ns < op.t_start
        || op.versions != versions_json(params.build_identity)
        || !op.ir.contains("\"op\":\"perf.roofline\"")
        || !op.ir.contains("\"admitted\":true")
        || !op
            .ir
            .contains(&format!("\"fingerprint\":\"{current_fingerprint:016x}\""))
        || !op.ir.contains(&format!(
            "\"finalized_run_receipt\":\"{}\"",
            params.run_receipt
        ))
        || !op
            .ir
            .contains(&format!("\"baseline_hash\":\"{}\"", params.baseline_hash))
    {
        return Ok(None);
    }

    let Some((_, post_and_later)) = op.ir.split_once("\"post\":") else {
        return Ok(None);
    };
    let Some((post, _)) = post_and_later.split_once(",\"baseline_hash\"") else {
        return Ok(None);
    };
    let post_names = [
        "bandwidth_single_bits",
        "bandwidth_all_core_bits",
        "peak_single_bits",
        "peak_all_core_bits",
    ];
    let post_matches = post_names
        .iter()
        .zip(params.post_axis_bits)
        .all(|(name, bits)| post.contains(&format!("\"{name}\":\"{bits:016x}\"")));
    Ok(post_matches.then_some(ValidatedRooflineRow {
        build_identity: params.build_identity,
        recorded_at_ns,
    }))
}

#[cfg(test)]
mod tests {
    use super::*;

    struct ReceiptKernel {
        elements: usize,
        value: u64,
    }

    impl RooflineKernel for ReceiptKernel {
        fn spec(&self) -> KernelSpec {
            KernelSpec {
                name: "receipt-kernel",
                version: "1",
                bytes_per_elem: 8.0,
                flops_per_elem: 1.0,
                threading: Threading::SingleThread,
                target_axis: TargetAxis::BindingRoof,
                target_fraction: None,
            }
        }

        fn elements(&self) -> usize {
            self.elements
        }

        fn run_once(&mut self) {
            for _ in 0..1024 {
                self.value = std::hint::black_box(
                    self.value
                        .wrapping_mul(6_364_136_223_846_793_005)
                        .wrapping_add(1),
                );
            }
        }
    }

    struct AdmissionProbeKernel {
        observed: std::rc::Rc<std::cell::Cell<Option<bool>>>,
    }

    impl RooflineKernel for AdmissionProbeKernel {
        fn spec(&self) -> KernelSpec {
            ReceiptKernel {
                elements: 1,
                value: 0,
            }
            .spec()
        }

        fn elements(&self) -> usize {
            1
        }

        fn run_once(&mut self) {}

        fn finalize_tuning(&mut self, admitted: bool) -> Result<(), String> {
            self.observed.set(Some(admitted));
            Ok(())
        }
    }

    struct FallibleAdmissionProbeKernel {
        id: usize,
        observed: std::rc::Rc<std::cell::RefCell<Vec<(usize, bool)>>>,
        failure: Option<&'static str>,
    }

    impl RooflineKernel for FallibleAdmissionProbeKernel {
        fn spec(&self) -> KernelSpec {
            let mut spec = ReceiptKernel {
                elements: 1,
                value: 0,
            }
            .spec();
            spec.name = ["receipt-kernel-0", "receipt-kernel-1", "receipt-kernel-2"]
                .get(self.id)
                .copied()
                .unwrap_or("receipt-kernel-out-of-range");
            spec
        }

        fn elements(&self) -> usize {
            1
        }

        fn run_once(&mut self) {
            let mut value = self.id as u64;
            for _ in 0..1024 {
                value = std::hint::black_box(
                    value
                        .wrapping_mul(6_364_136_223_846_793_005)
                        .wrapping_add(1),
                );
            }
            std::hint::black_box(value);
        }

        fn finalize_tuning(&mut self, admitted: bool) -> Result<(), String> {
            self.observed.borrow_mut().push((self.id, admitted));
            self.failure.map_or(Ok(()), |error| Err(error.to_string()))
        }
    }

    fn synthetic_axes() -> MachineAxes {
        MachineAxes {
            fingerprint: 0xABCD,
            cpu_brand: "synthetic".to_string(),
            logical_cpus: 8,
            bandwidth_single_gbs: 100.0,
            bandwidth_all_core_gbs: 400.0,
            peak_single_gflops: 50.0,
            peak_all_core_gflops: 300.0,
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
                        "fs-roofline.lib-test-baseline-source.v1",
                        &ordinal.to_le_bytes(),
                    ),
                )
                .expect("valid synthetic candidate")
            })
            .collect();
        let baseline = promote_baseline(
            &candidates,
            "test-operator",
            "deterministic lib receipt fixture",
            20_000,
            90,
        )
        .expect("valid synthetic baseline");
        (baseline, identity)
    }

    #[test]
    fn crushed_axes_cannot_gate_vacuously() {
        // Bead 1n61 counterexample replay: on a load-68 host both the
        // STREAM probe and the FFT kernel collapsed ~1000× together
        // (axes 0.2 GB/s, kernel 0.156 GB/s on a ~200 GB/s machine) and
        // the RATIO self-normalized to 0.89 = a vacuous within_band.
        let crushed = MachineAxes {
            fingerprint: 0xDEAD,
            cpu_brand: "crushed".to_string(),
            logical_cpus: 128,
            bandwidth_single_gbs: 0.2,
            bandwidth_all_core_gbs: 0.4,
            peak_single_gflops: 0.05,
            peak_all_core_gflops: 0.4,
        };
        assert!(!crushed.plausible());
        let spec = KernelSpec {
            name: "fft-roundtrip",
            version: "1n61",
            bytes_per_elem: 672.0,
            flops_per_elem: 172.0,
            threading: Threading::SingleThread,
            target_axis: TargetAxis::BindingRoof,
            target_fraction: Some(0.40),
        };
        // 0.156 GB/s effective — 89% of the crushed axis.
        let a = attainment_for(&spec, 0.156e9 / 672.0, &crushed);
        assert_eq!(
            a.verdict,
            Verdict::EnvironmentInvalid,
            "a crushed environment must refuse to gate (got attainment {:.3})",
            a.attainment
        );
    }

    #[test]
    fn over_roof_attainment_poisons_the_gate() {
        // Healthy axes but the kernel 'beats' its roofline by 2× — the
        // axes are stale relative to the run; refuse.
        let spec = KernelSpec {
            name: "axpy",
            version: "1n61",
            bytes_per_elem: 24.0,
            flops_per_elem: 2.0,
            threading: Threading::SingleThread,
            target_axis: TargetAxis::BindingRoof,
            target_fraction: Some(0.6),
        };
        let a = attainment_for(&spec, 2.0 * 100.0e9 / 24.0, &synthetic_axes());
        assert_eq!(a.verdict, Verdict::EnvironmentInvalid);
        // Slightly over 1 (measurement jitter) still gates normally.
        let b = attainment_for(&spec, 1.2 * 100.0e9 / 24.0, &synthetic_axes());
        assert_eq!(b.verdict, Verdict::WithinBand);
    }

    #[test]
    fn invalid_numeric_inputs_fail_closed_and_remain_json() {
        let base = KernelSpec {
            name: "probe\"escaped",
            version: "1",
            bytes_per_elem: 8.0,
            flops_per_elem: 1.0,
            threading: Threading::SingleThread,
            target_axis: TargetAxis::BindingRoof,
            target_fraction: Some(0.5),
        };
        let bad_rate = attainment_with_dispersion(&base, f64::NAN, 0.0, 3, &synthetic_axes());
        assert_eq!(bad_rate.verdict, Verdict::EnvironmentInvalid);
        assert!(bad_rate.elems_per_sec.is_finite());
        assert!(bad_rate.to_jsonl().contains("probe\\\"escaped"));
        assert!(!bad_rate.to_jsonl().contains("NaN"));

        let bad_dispersion =
            attainment_with_dispersion(&base, 1.0, f64::INFINITY, 3, &synthetic_axes());
        assert_eq!(bad_dispersion.verdict, Verdict::EnvironmentInvalid);
        let bad_target = KernelSpec {
            target_fraction: Some(f64::NAN),
            ..base
        };
        assert_eq!(
            attainment_for(&bad_target, 1.0, &synthetic_axes()).verdict,
            Verdict::EnvironmentInvalid
        );
    }

    #[test]
    fn one_invalid_row_poisons_every_registry_verdict() {
        let normal = KernelSpec {
            name: "normal",
            version: "1",
            bytes_per_elem: 24.0,
            flops_per_elem: 2.0,
            threading: Threading::SingleThread,
            target_axis: TargetAxis::BindingRoof,
            target_fraction: Some(0.5),
        };
        let impossible = KernelSpec {
            name: "impossible",
            ..normal
        };
        let mut results = vec![
            attainment_for(&normal, 100.0e9 / 24.0 * 0.8, &synthetic_axes()),
            attainment_for(&impossible, 100.0e9 / 24.0 * 2.0, &synthetic_axes()),
        ];
        assert_eq!(results[0].verdict, Verdict::WithinBand);
        assert_eq!(results[1].verdict, Verdict::EnvironmentInvalid);
        poison_invalid_run(&mut results);
        assert!(
            results
                .iter()
                .all(|row| row.verdict == Verdict::EnvironmentInvalid)
        );
        assert!(results.iter().all(|row| {
            row.invalid_reason
                .as_deref()
                .is_some_and(|r| r.contains("impossible"))
        }));
    }

    #[test]
    fn only_bound_timed_receipts_are_citable() {
        let axes = synthetic_axes();
        let (baseline, identity) = trusted_baseline(&axes);
        let baseline_policy = AxisBaselinePolicy::new(Some(&baseline), &identity, 20_010);
        let spec = ReceiptKernel {
            elements: 1,
            value: 0,
        }
        .spec();
        let analytic = attainment_with_dispersion(&spec, 1.0, 0.0, 1, &axes);
        assert!(
            !run_is_citable(&axes, &axes, baseline_policy, &[analytic]),
            "an analytic helper result is not measurement evidence"
        );

        let mut kernel = ReceiptKernel {
            elements: 1,
            value: 0,
        };
        let timed = measure(&mut kernel, 0, 3, &axes);
        assert!(run_is_citable(
            &axes,
            &axes,
            baseline_policy,
            std::slice::from_ref(&timed)
        ));
        assert!(timed.to_jsonl().contains("\"sample_seconds_bits\""));
        let mut drifted_post = axes.clone();
        drifted_post.bandwidth_single_gbs = 60.0;
        assert!(!run_is_citable(
            &axes,
            &drifted_post,
            baseline_policy,
            std::slice::from_ref(&timed)
        ));

        let mut tampered = timed.clone();
        tampered.dispersion += 0.01;
        assert!(!run_is_citable(&axes, &axes, baseline_policy, &[tampered]));
        let mut tampered_target = timed.clone();
        tampered_target.target_attainment += 0.01;
        assert!(!run_is_citable(
            &axes,
            &axes,
            baseline_policy,
            &[tampered_target]
        ));
        let mut tampered_axis = timed.clone();
        tampered_axis.spec_binding.target_axis = TargetAxis::ComputePeak;
        assert!(!run_is_citable(
            &axes,
            &axes,
            baseline_policy,
            &[tampered_axis]
        ));
        assert!(
            run_admission_error(
                &axes,
                &axes,
                baseline_policy,
                &[timed.clone(), timed.clone()],
            )
            .is_some()
        );

        let mut empty = ReceiptKernel {
            elements: 0,
            value: 0,
        };
        let empty_row = measure(&mut empty, 0, 1, &axes);
        assert!(!run_is_citable(&axes, &axes, baseline_policy, &[empty_row]));

        let unbaselined = AxisBaselinePolicy::new(None, &identity, 20_010);
        assert!(
            !run_is_citable(&axes, &axes, unbaselined, std::slice::from_ref(&timed)),
            "first-run candidate evidence must never authorize itself"
        );
    }

    #[test]
    fn stable_sustained_contention_cannot_cross_the_production_admission_boundary() {
        let quiet = synthetic_axes();
        let (baseline, identity) = trusted_baseline(&quiet);
        let baseline_policy = AxisBaselinePolicy::new(Some(&baseline), &identity, 20_010);
        let mut crushed = quiet.clone();
        crushed.bandwidth_single_gbs = 10.0;
        crushed.bandwidth_all_core_gbs = 40.0;
        crushed.peak_single_gflops = 10.0;
        crushed.peak_all_core_gflops = 60.0;
        assert!(crushed.plausible());
        assert!(crushed.reprobe_error(&crushed).is_none());

        let mut kernel = ReceiptKernel {
            elements: 1,
            value: 0,
        };
        let timed = measure(&mut kernel, 1, 3, &crushed);
        assert!(!run_is_citable(
            &crushed,
            &crushed,
            baseline_policy,
            std::slice::from_ref(&timed)
        ));
        let receipt = baseline_policy.receipt_json(&crushed, &crushed);
        assert!(receipt.contains("\"baseline\":\"degraded\""));
        assert!(receipt.contains(&baseline.content_hash().to_hex()));
        assert!(!receipt.contains("NaN"));
    }

    #[test]
    fn registry_tuning_hook_uses_the_complete_admission_decision() {
        let axes = synthetic_axes();
        let (baseline, identity) = trusted_baseline(&axes);
        let baseline_policy = AxisBaselinePolicy::new(Some(&baseline), &identity, 20_010);
        let mut measured = ReceiptKernel {
            elements: 1,
            value: 0,
        };
        let result = measure(&mut measured, 0, 3, &axes);
        let observed = std::rc::Rc::new(std::cell::Cell::new(None));
        let mut registry: Vec<Box<dyn RooflineKernel>> = vec![Box::new(AdmissionProbeKernel {
            observed: std::rc::Rc::clone(&observed),
        })];

        assert!(
            finalize_registry_tuning(
                &mut registry,
                &axes,
                &axes,
                baseline_policy,
                std::slice::from_ref(&result),
            )
            .expect("admitted hook")
            .admitted()
        );
        assert_eq!(observed.get(), Some(true));

        observed.set(None);
        let mut drifted_post = axes.clone();
        drifted_post.bandwidth_single_gbs *= 0.5;
        assert!(
            !finalize_registry_tuning(
                &mut registry,
                &axes,
                &drifted_post,
                baseline_policy,
                &[result],
            )
            .expect("rejected hook")
            .admitted()
        );
        assert_eq!(observed.get(), Some(false));
    }

    #[test]
    fn registry_tuning_hook_drains_every_kernel_after_a_middle_failure() {
        let axes = synthetic_axes();
        let (baseline, identity) = trusted_baseline(&axes);
        let baseline_policy = AxisBaselinePolicy::new(Some(&baseline), &identity, 20_010);
        let observed = std::rc::Rc::new(std::cell::RefCell::new(Vec::new()));
        let mut registry: Vec<Box<dyn RooflineKernel>> = (0..3)
            .map(|id| {
                Box::new(FallibleAdmissionProbeKernel {
                    id,
                    observed: std::rc::Rc::clone(&observed),
                    failure: (id == 1).then_some("middle cleanup failed"),
                }) as Box<dyn RooflineKernel>
            })
            .collect();
        let results = run_registry(&mut registry, 0, 3, &axes);

        let error =
            finalize_registry_tuning(&mut registry, &axes, &axes, baseline_policy, &results)
                .expect_err("middle failure must be reported after every hook drains");
        assert_eq!(
            observed.borrow().as_slice(),
            &[(0, true), (1, true), (2, true)],
            "first, failing middle, and last kernel must see the same admission decision"
        );
        assert_eq!(
            error,
            "tuning lifecycle finalization failed for 1 kernel(s): kernel[1]: middle cleanup failed"
        );
    }

    #[test]
    fn registry_finalization_refuses_missing_or_mismatched_kernels() {
        let axes = synthetic_axes();
        let (baseline, identity) = trusted_baseline(&axes);
        let baseline_policy = AxisBaselinePolicy::new(Some(&baseline), &identity, 20_010);
        let mut measured = ReceiptKernel {
            elements: 1,
            value: 0,
        };
        let result = measure(&mut measured, 0, 1, &axes);

        let missing = finalize_registry_tuning(
            &mut [],
            &axes,
            &axes,
            baseline_policy,
            std::slice::from_ref(&result),
        )
        .expect_err("an empty registry cannot finalize a nonempty result set");
        assert!(missing.contains("length mismatch"));

        let observed = std::rc::Rc::new(std::cell::Cell::new(None));
        let mut mismatched: Vec<Box<dyn RooflineKernel>> = vec![Box::new(AdmissionProbeKernel {
            observed: std::rc::Rc::clone(&observed),
        })];
        let mut mismatched_result = result;
        mismatched_result.version = "different-version".to_string();
        let mismatch = finalize_registry_tuning(
            &mut mismatched,
            &axes,
            &axes,
            baseline_policy,
            &[mismatched_result],
        )
        .expect_err("a different kernel cannot finalize the measured row");
        assert!(mismatch.contains("identity mismatch"));
        assert_eq!(
            observed.get(),
            Some(false),
            "identity mismatch must be part of the aggregate admission decision"
        );
    }

    #[test]
    fn attainment_matches_hand_calculation_bandwidth_bound() {
        // axpy: 24 B/elem, 2 flop/elem on axes (100 GB/s, 50 GFLOP/s):
        // bw limit = 100e9/24 = 4.1667e9 elem/s; compute = 50e9/2 = 25e9.
        // Bandwidth binds. At 2.0833e9 elem/s attainment = 0.5 exactly.
        let spec = KernelSpec {
            name: "axpy",
            version: "1",
            bytes_per_elem: 24.0,
            flops_per_elem: 2.0,
            threading: Threading::SingleThread,
            target_axis: TargetAxis::BindingRoof,
            target_fraction: Some(0.6),
        };
        let a = attainment_for(&spec, 100.0e9 / 24.0 / 2.0, &synthetic_axes());
        assert_eq!(a.roof, RoofSide::Bandwidth);
        assert!((a.attainment - 0.5).abs() < 1e-12, "got {}", a.attainment);
        assert!((a.achieved_gbs - 50.0).abs() < 1e-9);
        assert_eq!(a.verdict, Verdict::BelowBand);
    }

    #[test]
    fn attainment_matches_hand_calculation_compute_bound() {
        // High-intensity kernel: 1 B/elem, 100 flop/elem.
        // bw limit = 100e9 elem/s; compute = 50e9/100 = 0.5e9 → compute binds.
        let spec = KernelSpec {
            name: "dense",
            version: "1",
            bytes_per_elem: 1.0,
            flops_per_elem: 100.0,
            threading: Threading::SingleThread,
            target_axis: TargetAxis::BindingRoof,
            target_fraction: Some(0.5),
        };
        let a = attainment_for(&spec, 0.4e9, &synthetic_axes());
        assert_eq!(a.roof, RoofSide::Compute);
        assert!((a.attainment - 0.8).abs() < 1e-12);
        assert_eq!(a.verdict, Verdict::WithinBand);
        // All-core axes flip the limit.
        let all = KernelSpec {
            threading: Threading::AllCore,
            ..spec
        };
        let b = attainment_for(&all, 0.4e9, &synthetic_axes());
        assert!((b.limit_elems_per_sec - 3.0e9).abs() < 1.0);
    }

    #[test]
    fn compute_target_is_not_laundered_through_the_binding_bandwidth_roof() {
        let bandwidth_starved = MachineAxes {
            fingerprint: 0xC0DE,
            cpu_brand: "high-compute-low-bandwidth".to_string(),
            logical_cpus: 8,
            bandwidth_single_gbs: 10.0,
            bandwidth_all_core_gbs: 20.0,
            peak_single_gflops: 1_000.0,
            peak_all_core_gflops: 2_000.0,
        };
        let spec = KernelSpec {
            name: "gemm-target-model",
            version: "ss0n",
            bytes_per_elem: 8.0,
            flops_per_elem: 1.0,
            threading: Threading::SingleThread,
            target_axis: TargetAxis::ComputePeak,
            target_fraction: Some(0.75),
        };
        let above_bandwidth_target = attainment_for(&spec, 1.0e9, &bandwidth_starved);
        assert_eq!(above_bandwidth_target.roof, RoofSide::Bandwidth);
        assert!((above_bandwidth_target.attainment - 0.8).abs() < 1e-12);
        assert!((above_bandwidth_target.target_attainment - 0.001).abs() < 1e-12);
        assert_eq!(above_bandwidth_target.verdict, Verdict::BelowBand);
        assert!(
            above_bandwidth_target
                .to_jsonl()
                .contains("\"target_axis\":\"compute_peak\"")
        );

        let compute_bound = KernelSpec {
            bytes_per_elem: 1.0,
            flops_per_elem: 100.0,
            ..spec
        };
        let exact_boundary = attainment_for(&compute_bound, 0.375e9, &synthetic_axes());
        assert_eq!(exact_boundary.roof, RoofSide::Compute);
        assert_eq!(
            exact_boundary.target_attainment.to_bits(),
            0.75f64.to_bits()
        );
        assert_eq!(exact_boundary.verdict, Verdict::WithinBand);
    }

    #[test]
    fn no_target_reports_without_verdict() {
        let spec = KernelSpec {
            name: "probe",
            version: "1",
            bytes_per_elem: 8.0,
            flops_per_elem: 1.0,
            threading: Threading::SingleThread,
            target_axis: TargetAxis::BindingRoof,
            target_fraction: None,
        };
        let a = attainment_for(&spec, 1.0e9, &synthetic_axes());
        assert_eq!(a.verdict, Verdict::NoTarget);
        assert!(a.to_jsonl().contains("\"verdict\":\"no_target\""));
    }

    #[test]
    fn section_14_1_table_is_complete_and_honest() {
        assert_eq!(SECTION_14_1_TARGETS.len(), 7, "all §14.1 families present");
        let registry = kernels::production_registry(1, &synthetic_axes());
        let registered: std::collections::BTreeSet<_> =
            registry.iter().map(|kernel| kernel.spec().name).collect();
        for row in SECTION_14_1_TARGETS {
            assert_eq!(
                row.landed,
                registered.contains(row.kernel),
                "{} target registration and landed status drifted",
                row.kernel,
            );
        }
    }

    #[test]
    fn versions_json_binds_the_exact_executable_identity() {
        let identity = fs_blake3::hash_domain("fs-roofline.test-executable.v1", b"binary");
        assert_eq!(
            versions_json(identity),
            format!("{{\"frankensim_executable\":\"{identity}\",\"fs-roofline\":\"{VERSION}\"}}")
        );
    }

    #[test]
    fn staleness_refuses_a_different_current_executable() {
        let axes = synthetic_axes();
        let (baseline, identity) = trusted_baseline(&axes);
        let baseline_policy = AxisBaselinePolicy::new(Some(&baseline), &identity, 20_010);
        let mut registry: Vec<Box<dyn RooflineKernel>> = vec![Box::new(ReceiptKernel {
            elements: 1,
            value: 0,
        })];
        let mut results = run_registry(&mut registry, 0, 1, &axes);
        let mut finalized =
            finalize_registry_tuning(&mut registry, &axes, &axes, baseline_policy, &results)
                .expect("finalize fixture");
        let ledger = Ledger::open(":memory:").expect("in-memory ledger");
        let op = record_run(
            &ledger,
            &axes,
            &axes,
            baseline_policy,
            &mut finalized,
            &mut results,
        )
        .expect("record fixture");
        let recorded_at = ledger
            .op(op)
            .expect("query op")
            .expect("stored op")
            .t_end
            .expect("finished op");
        let foreign_build = fs_blake3::hash_domain(
            "org.frankensim.fs-roofline.foreign-executable.test.v1",
            b"rebuilt binary",
        );
        assert_eq!(
            staleness_at_with_build(
                &ledger,
                &results[0].kernel,
                &results[0].version,
                axes.fingerprint,
                Some(baseline.content_hash()),
                recorded_at,
                foreign_build,
            )
            .expect("classify"),
            Staleness::BuildDrift
        );
    }

    #[test]
    fn bit_identical_reruns_keep_distinct_operation_bound_rows() {
        let axes = synthetic_axes();
        let (baseline, identity) = trusted_baseline(&axes);
        let baseline_policy = AxisBaselinePolicy::new(Some(&baseline), &identity, 20_010);
        let mut first_registry: Vec<Box<dyn RooflineKernel>> = vec![Box::new(ReceiptKernel {
            elements: 1,
            value: 0,
        })];
        let mut first_results = run_registry(&mut first_registry, 0, 1, &axes);
        let mut second_results = first_results.clone();
        let mut first_finalized = finalize_registry_tuning(
            &mut first_registry,
            &axes,
            &axes,
            baseline_policy,
            &first_results,
        )
        .expect("finalize first identical run");
        let mut second_registry: Vec<Box<dyn RooflineKernel>> = vec![Box::new(ReceiptKernel {
            elements: 1,
            value: 0,
        })];
        let mut second_finalized = finalize_registry_tuning(
            &mut second_registry,
            &axes,
            &axes,
            baseline_policy,
            &second_results,
        )
        .expect("finalize second identical run");
        assert_eq!(
            first_finalized.receipt_identity(),
            second_finalized.receipt_identity(),
            "fixture must exercise an exact receipt collision"
        );

        let ledger = Ledger::open(":memory:").expect("in-memory ledger");
        let first_op = record_run(
            &ledger,
            &axes,
            &axes,
            baseline_policy,
            &mut first_finalized,
            &mut first_results,
        )
        .expect("record first identical run");
        let second_op = record_run(
            &ledger,
            &axes,
            &axes,
            baseline_policy,
            &mut second_finalized,
            &mut second_results,
        )
        .expect("record second identical run");
        assert_ne!(first_op, second_op);

        let rows = ledger
            .tune_rows("receipt-kernel")
            .expect("query retained run history");
        assert_eq!(rows.len(), 2);
        assert!(rows.iter().any(|row| {
            row.shape_class
                == tune_measurement_shape_class("1", first_finalized.receipt_identity(), first_op)
        }));
        assert!(rows.iter().any(|row| {
            row.shape_class
                == tune_measurement_shape_class("1", second_finalized.receipt_identity(), second_op)
        }));
    }
}
