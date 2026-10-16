//! Ignored D3Q19 sparse-sweep throughput lane (bead 712t).
//!
//! Run explicitly in release mode; ordinary test suites only exercise the
//! bounded setup/classification helpers:
//!
//! ```text
//! FRANKENSIM_BASELINE_STORE=<jsonl> \
//! FRANKENSIM_FIRMWARE_ID=<id> \
//! FRANKENSIM_PROMOTION_AUTHORITY_POLICY=<tsv> \
//! FRANKENSIM_RETAINED_SOURCE_RECEIPTS=<txt> \
//! cargo test -p fs-lbm --release --test perf_lane -- --ignored --nocapture
//! ```
//!
//! This first measurement driver intentionally emits report-only rows even
//! when the axis baseline is authority-admitted: no anti-collapse floor has
//! yet been calibrated and fs-roofline's external-gate recorder does not yet
//! own an LBM lane. Raw plan targets remain informational. A later retained
//! both-ISA calibration must authorize the floor and ledger schema before any
//! positive or negative performance verdict becomes citable.

use std::collections::BTreeSet;
use std::io::Read as _;

use fs_exec::{CancelGate, TilePool};
use fs_lbm::d3q19::sparse::{SparseGrid3, SparseSweepObservation, morton3};
use fs_lbm::perf::{
    D3Q19_PERF_MODEL_VERSION, D3q19PerfRow, D3q19TrafficModel, EvidenceClass, LaneShape,
    OccupancyClass, RATIO_PPM, ReferenceIsa, ThreadingClass, attribute_critical_path,
    sparse_sweep_task_samples,
};
use fs_roofline::authority::{ConfiguredPromotionAuthority, MAX_PROMOTION_AUTHORITY_POLICY_BYTES};
use fs_roofline::{
    Attainment, AttestedAxisBaselinePolicy, AttestedBaselineStore, AxisAdmissionSnapshot,
    AxisBaselinePolicy, BaselineAxes, BaselineIdentity, BaselineStore, ContentHash, KernelSpec,
    MachineAxes, RooflineKernel, TargetAxis, Threading, Verdict, days_since_epoch_now, measure,
};

const OBS_SUITE: &str = "fs-lbm/d3q19-perf-lane";
const WARMUP_REPETITIONS: usize = 1;
const TIMED_REPETITIONS: usize = 5;
const MEASUREMENT_INVOCATIONS: usize = WARMUP_REPETITIONS + TIMED_REPETITIONS;
const MAX_RETAINED_RECEIPT_INPUT_BYTES: usize = fs_roofline::baseline::MAX_BASELINE_STORE_BYTES;

fn emit_observation(identity: &str, name: &str, severity: fs_obs::Severity, json: String) {
    let mut emitter = fs_obs::Emitter::new(OBS_SUITE, identity);
    let event = emitter.emit(
        severity,
        fs_obs::EventKind::Custom {
            name: name.to_owned(),
            json,
        },
        None,
    );
    fs_obs::lint_failure_record(&event).expect("performance diagnostic must be replayable");
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("performance diagnostic must use the fs-obs wire schema");
    println!("{line}");
}

fn json_escape(value: &str) -> String {
    use core::fmt::Write as _;

    let mut escaped = String::with_capacity(value.len());
    for ch in value.chars() {
        match ch {
            '"' => escaped.push_str("\\\""),
            '\\' => escaped.push_str("\\\\"),
            '\n' => escaped.push_str("\\n"),
            '\r' => escaped.push_str("\\r"),
            '\t' => escaped.push_str("\\t"),
            c if c.is_control() => {
                let _ = write!(escaped, "\\u{:04x}", u32::from(c));
            }
            c => escaped.push(c),
        }
    }
    escaped
}

fn read_bounded_text(path: &str, kind: &str, limit: usize) -> Result<String, String> {
    let file = std::fs::File::open(path)
        .map_err(|error| format!("cannot read {kind} {path:?}: {error}"))?;
    let bounded_bytes = limit
        .checked_add(1)
        .ok_or_else(|| format!("{kind} read bound overflows usize"))?;
    let read_limit =
        u64::try_from(bounded_bytes).map_err(|_| format!("{kind} read bound does not fit u64"))?;
    let mut bytes = Vec::with_capacity(bounded_bytes);
    file.take(read_limit)
        .read_to_end(&mut bytes)
        .map_err(|error| format!("cannot read {kind} {path:?}: {error}"))?;
    if bytes.len() > limit {
        return Err(format!("{kind} {path:?} exceeds the {limit}-byte bound"));
    }
    String::from_utf8(bytes).map_err(|_| format!("{kind} {path:?} is not UTF-8"))
}

fn parse_retained_receipts(text: &str) -> Result<BTreeSet<ContentHash>, String> {
    let body = text.strip_suffix('\n').ok_or_else(|| {
        "retained source receipts must be canonical newline-terminated lowercase hex".to_owned()
    })?;
    if body.is_empty() {
        return Err("retained source receipts must contain at least one receipt".to_owned());
    }
    let mut receipts = BTreeSet::new();
    let mut previous = None;
    for (index, line) in body.split('\n').enumerate() {
        if line.len() != 64
            || !line
                .bytes()
                .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
        {
            return Err(format!(
                "retained source receipt line {} must be exactly 64 lowercase hexadecimal bytes",
                index + 1
            ));
        }
        let receipt = ContentHash::from_hex(line).ok_or_else(|| {
            format!(
                "retained source receipt line {} is not a content hash",
                index + 1
            )
        })?;
        if previous.is_some_and(|prior| receipt <= prior) {
            return Err(format!(
                "retained source receipt line {} is not in strict ascending order",
                index + 1
            ));
        }
        previous = Some(receipt);
        let inserted = receipts.insert(receipt);
        debug_assert!(inserted);
    }
    Ok(receipts)
}

enum PreparedAdmission {
    Attested(AttestedAxisBaselinePolicy),
    ReportOnly {
        baseline: Option<BaselineAxes>,
        identity: BaselineIdentity,
        now_day: u64,
        refusal: String,
    },
}

impl PreparedAdmission {
    fn report_only(
        baseline: Option<BaselineAxes>,
        identity: BaselineIdentity,
        now_day: u64,
        refusal: impl Into<String>,
    ) -> Self {
        Self::ReportOnly {
            baseline,
            identity,
            now_day,
            refusal: refusal.into(),
        }
    }

    fn snapshot(
        self,
        pre: &MachineAxes,
        post: &MachineAxes,
    ) -> (AxisAdmissionSnapshot, Option<String>) {
        match self {
            Self::Attested(policy) => (policy.decide(pre, post), None),
            Self::ReportOnly {
                baseline,
                identity,
                now_day,
                refusal,
            } => (
                AxisBaselinePolicy::new(baseline.as_ref(), &identity, now_day).snapshot(pre, post),
                Some(refusal),
            ),
        }
    }
}

#[derive(Debug, PartialEq, Eq)]
enum GateAdmission {
    Citable,
    ReportOnly(String),
    EnvironmentInvalid(String),
}

fn classify_gate_admission(
    snapshot: &AxisAdmissionSnapshot,
    configuration_refusal: Option<String>,
    pre: &MachineAxes,
    post: &MachineAxes,
) -> GateAdmission {
    if let Some(reason) = pre.reprobe_error(post) {
        return GateAdmission::EnvironmentInvalid(reason.to_owned());
    }
    if let Some(reason) = configuration_refusal {
        return GateAdmission::ReportOnly(reason);
    }
    match snapshot.baseline_citation_error() {
        Some(reason) => GateAdmission::EnvironmentInvalid(reason),
        None => GateAdmission::Citable,
    }
}

fn prepare_configured_attested_admission(
    store: &AttestedBaselineStore,
    identity: BaselineIdentity,
    now_day: u64,
    authority_text: &str,
    receipts_text: &str,
) -> PreparedAdmission {
    let candidate = store.for_fingerprint(identity.fingerprint()).cloned();
    let authority = match ConfiguredPromotionAuthority::from_text(authority_text) {
        Ok(authority) => authority,
        Err(error) => {
            return PreparedAdmission::report_only(
                candidate,
                identity,
                now_day,
                format!("invalid promotion-authority policy: {error}"),
            );
        }
    };
    let receipts = match parse_retained_receipts(receipts_text) {
        Ok(receipts) => receipts,
        Err(error) => {
            return PreparedAdmission::report_only(candidate, identity, now_day, error);
        }
    };
    match store.policy_for_run(&identity, &authority, &receipts) {
        Ok(policy) => PreparedAdmission::Attested(policy),
        Err(error) => PreparedAdmission::report_only(
            candidate,
            identity,
            now_day,
            format!("attested baseline authority refused: {error}"),
        ),
    }
}

fn report_only_day(refusal: impl Into<String>) -> (u64, String) {
    let refusal = refusal.into();
    match days_since_epoch_now() {
        Ok(day) => (day, refusal),
        Err(error) => (
            0,
            format!("{refusal}; cannot establish baseline age: {error}"),
        ),
    }
}

#[allow(clippy::too_many_lines)]
fn prepare_admission(axes: &MachineAxes) -> PreparedAdmission {
    let firmware = match std::env::var("FRANKENSIM_FIRMWARE_ID") {
        Ok(value) if !value.is_empty() => value,
        _ => {
            let identity = BaselineIdentity::current(axes, "unbaselined-candidate")
                .expect("plausible probed axes form a candidate identity");
            let (now_day, refusal) = report_only_day("FRANKENSIM_FIRMWARE_ID is missing or empty");
            return PreparedAdmission::report_only(None, identity, now_day, refusal);
        }
    };
    let identity = match BaselineIdentity::current(axes, firmware) {
        Ok(identity) => identity,
        Err(error) => {
            let identity = BaselineIdentity::current(axes, "unbaselined-candidate")
                .expect("plausible probed axes form a candidate identity");
            let (now_day, refusal) = report_only_day(format!("invalid baseline identity: {error}"));
            return PreparedAdmission::report_only(None, identity, now_day, refusal);
        }
    };
    let now_day = match days_since_epoch_now() {
        Ok(day) => day,
        Err(error) => {
            return PreparedAdmission::report_only(
                None,
                identity,
                0,
                format!("cannot establish baseline age: {error}"),
            );
        }
    };
    let baseline_path = match std::env::var("FRANKENSIM_BASELINE_STORE") {
        Ok(path) if !path.is_empty() => path,
        _ => {
            return PreparedAdmission::report_only(
                None,
                identity,
                now_day,
                "FRANKENSIM_BASELINE_STORE is missing or empty",
            );
        }
    };
    let baseline_text = match read_bounded_text(
        &baseline_path,
        "baseline store",
        fs_roofline::baseline::MAX_BASELINE_STORE_BYTES,
    ) {
        Ok(text) => text,
        Err(error) => {
            return PreparedAdmission::report_only(None, identity, now_day, error);
        }
    };
    if !baseline_text.starts_with("{\"record\":") {
        return match BaselineStore::from_jsonl(&baseline_text) {
            Ok(store) => PreparedAdmission::report_only(
                store.for_fingerprint(axes.fingerprint).cloned(),
                identity,
                now_day,
                "plain baseline stores are candidate/report-only inputs",
            ),
            Err(error) => PreparedAdmission::report_only(
                None,
                identity,
                now_day,
                format!("invalid plain baseline store: {error}"),
            ),
        };
    }
    let store = match AttestedBaselineStore::from_jsonl(&baseline_text) {
        Ok(store) => store,
        Err(error) => {
            return PreparedAdmission::report_only(
                None,
                identity,
                now_day,
                format!("invalid attested baseline store: {error}"),
            );
        }
    };
    let candidate = store.for_fingerprint(axes.fingerprint).cloned();
    let authority_path = match std::env::var("FRANKENSIM_PROMOTION_AUTHORITY_POLICY") {
        Ok(path) if !path.is_empty() => path,
        _ => {
            return PreparedAdmission::report_only(
                candidate,
                identity,
                now_day,
                "FRANKENSIM_PROMOTION_AUTHORITY_POLICY is missing or empty",
            );
        }
    };
    let authority_text = match read_bounded_text(
        &authority_path,
        "promotion-authority policy",
        MAX_PROMOTION_AUTHORITY_POLICY_BYTES,
    ) {
        Ok(text) => text,
        Err(error) => {
            return PreparedAdmission::report_only(candidate, identity, now_day, error);
        }
    };
    let receipts_path = match std::env::var("FRANKENSIM_RETAINED_SOURCE_RECEIPTS") {
        Ok(path) if !path.is_empty() => path,
        _ => {
            return PreparedAdmission::report_only(
                candidate,
                identity,
                now_day,
                "FRANKENSIM_RETAINED_SOURCE_RECEIPTS is missing or empty",
            );
        }
    };
    let receipts_text = match read_bounded_text(
        &receipts_path,
        "retained source receipts",
        MAX_RETAINED_RECEIPT_INPUT_BYTES,
    ) {
        Ok(text) => text,
        Err(error) => {
            return PreparedAdmission::report_only(candidate, identity, now_day, error);
        }
    };
    prepare_configured_attested_admission(
        &store,
        identity,
        now_day,
        &authority_text,
        &receipts_text,
    )
}

fn reference_isa_from(cpu_brand: &str, os: &str, arch: &str) -> ReferenceIsa {
    let brand = cpu_brand.to_ascii_lowercase();
    if os == "macos" && arch == "aarch64" && brand.contains("apple") {
        ReferenceIsa::AppleMClass
    } else if arch == "x86_64" && (brand.contains("threadripper") || brand.contains("epyc")) {
        ReferenceIsa::ThreadripperClass
    } else {
        ReferenceIsa::Other
    }
}

fn reference_isa(axes: &MachineAxes) -> ReferenceIsa {
    reference_isa_from(
        &axes.cpu_brand,
        std::env::consts::OS,
        std::env::consts::ARCH,
    )
}

fn active_coordinates(shape: LaneShape) -> Result<Vec<(u32, u32, u32)>, String> {
    let total = shape.total_tiles().map_err(|error| error.to_string())?;
    let active = shape.active_tiles().map_err(|error| error.to_string())?;
    let [nx, ny, nz] = shape.dims.map(|dim| dim / fs_lbm::d3q19::TILE);
    let ntx = u32::try_from(nx).map_err(|_| "x tile dimension exceeds u32".to_owned())?;
    let nty = u32::try_from(ny).map_err(|_| "y tile dimension exceeds u32".to_owned())?;
    let ntz = u32::try_from(nz).map_err(|_| "z tile dimension exceeds u32".to_owned())?;

    let mut coordinates = Vec::new();
    coordinates
        .try_reserve_exact(total)
        .map_err(|_| format!("cannot reserve {total} active-tile coordinates"))?;
    for tz in 0..ntz {
        for ty in 0..nty {
            for tx in 0..ntx {
                coordinates.push((tx, ty, tz));
            }
        }
    }
    if coordinates.len() != total {
        return Err("tile coordinate cardinality disagrees with lane shape".to_owned());
    }
    coordinates.sort_unstable_by_key(|&(tx, ty, tz)| morton3(tx, ty, tz));
    coordinates.truncate(active);
    Ok(coordinates)
}

fn kernel_name(shape: LaneShape) -> &'static str {
    match (shape.occupancy, shape.threading) {
        (OccupancyClass::DenseActive, ThreadingClass::SingleThread) => "d3q19-sweep-dense-single",
        (OccupancyClass::DenseActive, ThreadingClass::AllCore) => "d3q19-sweep-dense-all-core",
        (OccupancyClass::SparseTenPercent, ThreadingClass::SingleThread) => {
            "d3q19-sweep-sparse-single"
        }
        (OccupancyClass::SparseTenPercent, ThreadingClass::AllCore) => {
            "d3q19-sweep-sparse-all-core"
        }
    }
}

struct SparseSweepKernel {
    grid: SparseGrid3,
    pool: TilePool,
    gate: CancelGate,
    shape: LaneShape,
    placement_identity: String,
    observations: Vec<SparseSweepObservation>,
}

impl SparseSweepKernel {
    fn new(shape: LaneShape) -> Result<Self, String> {
        shape.validate().map_err(|error| error.to_string())?;
        let coordinates = active_coordinates(shape)?;
        let mut grid = SparseGrid3::new(shape.dims[0], shape.dims[1], shape.dims[2], 0.8, [0.0; 3])
            .map_err(|error| error.to_string())?;
        grid.activate_tiles(&coordinates)
            .map_err(|error| error.to_string())?;
        let active_tiles = shape.active_tiles().map_err(|error| error.to_string())?;
        let seed = 0x7120_d3a1_9000_0000_u64
            ^ u64::try_from(active_tiles).map_err(|_| "active tile count exceeds u64")?
            ^ u64::try_from(shape.workers).map_err(|_| "worker count exceeds u64")?;
        let pool = TilePool::for_host(shape.workers, seed);
        let placement_identity = format!(
            "{};tile-edge={};active-tiles={active_tiles};runner=spawn-per-pass",
            pool.placement_identity(),
            fs_lbm::d3q19::TILE,
        );
        let mut observations = Vec::new();
        observations
            .try_reserve_exact(MEASUREMENT_INVOCATIONS)
            .map_err(|_| "cannot reserve bounded sparse-sweep observations".to_owned())?;
        Ok(Self {
            grid,
            pool,
            gate: CancelGate::new(),
            shape,
            placement_identity,
            observations,
        })
    }

    fn take_timed_observations(&mut self) -> Result<Vec<SparseSweepObservation>, String> {
        if self.observations.len() != MEASUREMENT_INVOCATIONS {
            return Err(format!(
                "measurement retained {}/{} warmup+timed observations",
                self.observations.len(),
                MEASUREMENT_INVOCATIONS
            ));
        }
        self.observations.drain(..WARMUP_REPETITIONS);
        Ok(core::mem::take(&mut self.observations))
    }
}

impl RooflineKernel for SparseSweepKernel {
    fn spec(&self) -> KernelSpec {
        let model = D3q19TrafficModel::default();
        KernelSpec {
            name: kernel_name(self.shape),
            version: D3Q19_PERF_MODEL_VERSION,
            bytes_per_elem: model.bytes_per_cell(),
            flops_per_elem: f64::from(model.flops_per_cell),
            threading: match self.shape.threading {
                ThreadingClass::SingleThread => Threading::SingleThread,
                ThreadingClass::AllCore => Threading::AllCore,
            },
            target_axis: TargetAxis::MemoryBandwidth,
            target_fraction: None,
        }
    }

    fn elements(&self) -> usize {
        self.shape
            .active_cells()
            .expect("constructor validated the bounded lane shape")
    }

    fn run_once(&mut self) -> Result<(), String> {
        if self.observations.len() >= MEASUREMENT_INVOCATIONS {
            return Err("roofline requested more than the bounded sweep repetitions".to_owned());
        }
        let observation = self
            .grid
            .step_pooled_observed(&self.pool, &self.gate)
            .map_err(|error| error.to_string())?;
        self.observations.push(observation);
        Ok(())
    }
}

struct MeasuredLane {
    attainment: Attainment,
    shape: LaneShape,
    placement_identity: String,
    critical_paths: Vec<fs_lbm::perf::CriticalPathAttribution>,
}

fn measure_lane(shape: LaneShape, axes: &MachineAxes) -> Result<MeasuredLane, String> {
    let mut kernel = SparseSweepKernel::new(shape)?;
    let attainment = measure(&mut kernel, WARMUP_REPETITIONS, TIMED_REPETITIONS, axes)?;
    let observations = kernel.take_timed_observations()?;
    if observations.len() != attainment.reps {
        return Err("roofline repetition count disagrees with retained sweep telemetry".to_owned());
    }
    let mut critical_paths = Vec::new();
    critical_paths
        .try_reserve_exact(observations.len())
        .map_err(|_| "cannot reserve bounded critical-path rows".to_owned())?;
    for observation in &observations {
        let tasks =
            sparse_sweep_task_samples(None, observation).map_err(|error| error.to_string())?;
        critical_paths.push(attribute_critical_path(&tasks).map_err(|error| error.to_string())?);
    }
    Ok(MeasuredLane {
        attainment,
        shape,
        placement_identity: kernel.placement_identity,
        critical_paths,
    })
}

fn measurement_invalid_reason(attainment: &Attainment) -> Option<String> {
    if attainment.verdict == Verdict::EnvironmentInvalid {
        return Some(
            attainment
                .invalid_reason
                .clone()
                .unwrap_or_else(|| "roofline marked the timed row environment-invalid".to_owned()),
        );
    }
    if !attainment.elems_per_sec.is_finite() || attainment.elems_per_sec <= 0.0 {
        return Some("D3Q19 throughput is non-finite or non-positive".to_owned());
    }
    if !attainment.dispersion.is_finite()
        || attainment.dispersion < 0.0
        || attainment.dispersion > 1.0
    {
        return Some(format!(
            "D3Q19 interquartile dispersion is outside [0,1]: {}",
            attainment.dispersion
        ));
    }
    None
}

fn row_evidence(admission: &GateAdmission, attainment: &Attainment) -> EvidenceClass {
    if let Some(reason) = measurement_invalid_reason(attainment) {
        return EvidenceClass::EnvironmentInvalid { reason };
    }
    match admission {
        GateAdmission::Citable => EvidenceClass::ReportOnly {
            reason: "axis baseline admitted, but no LBM anti-collapse floor or durable external-gate lane is authorized"
                .to_owned(),
        },
        GateAdmission::ReportOnly(reason) => EvidenceClass::ReportOnly {
            reason: reason.clone(),
        },
        GateAdmission::EnvironmentInvalid(reason) => EvidenceClass::EnvironmentInvalid {
            reason: reason.clone(),
        },
    }
}

#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
fn dispersion_ppm(attainment: &Attainment) -> u32 {
    let bounded = if attainment.dispersion.is_finite() {
        attainment.dispersion.clamp(0.0, 1.0)
    } else {
        1.0
    };
    (bounded * f64::from(RATIO_PPM)).round() as u32
}

fn lane_identity(shape: LaneShape) -> String {
    format!("{}/{}", shape.occupancy.as_str(), shape.threading.as_str())
}

fn fail_invalid_environment(identity: &str, reason: &str) -> ! {
    emit_observation(
        identity,
        "environment-invalid",
        fs_obs::Severity::Error,
        format!(
            "{{\"metric\":\"lbm-d3q19-gate\",\"verdict\":\"environment_invalid\",\"reason\":\"{}\"}}",
            json_escape(reason)
        ),
    );
    panic!("D3Q19 performance evidence rejected: {reason}");
}

#[test]
fn reference_isa_classifier_is_conservative() {
    assert_eq!(
        reference_isa_from("Apple M4 Max", "macos", "aarch64"),
        ReferenceIsa::AppleMClass
    );
    assert_eq!(
        reference_isa_from("AMD Ryzen Threadripper PRO 7995WX", "linux", "x86_64"),
        ReferenceIsa::ThreadripperClass
    );
    assert_eq!(
        reference_isa_from("AMD EPYC 9654", "linux", "x86_64"),
        ReferenceIsa::ThreadripperClass
    );
    assert_eq!(
        reference_isa_from("Apple M4 Max", "linux", "aarch64"),
        ReferenceIsa::Other
    );
    assert_eq!(
        reference_isa_from("AMD Ryzen 9", "linux", "x86_64"),
        ReferenceIsa::Other
    );
}

#[test]
fn sparse_lane_selects_exact_morton_prefix() {
    let shape = LaneShape::memory_resident(
        OccupancyClass::SparseTenPercent,
        ThreadingClass::SingleThread,
        1,
    )
    .expect("shape admits");
    let coordinates = active_coordinates(shape).expect("coordinate plan admits");
    assert_eq!(
        coordinates.len(),
        shape.active_tiles().expect("active tiles")
    );
    assert!(coordinates.windows(2).all(|pair| {
        morton3(pair[0].0, pair[0].1, pair[0].2) < morton3(pair[1].0, pair[1].1, pair[1].2)
    }));
}

#[test]
#[ignore = "hardware perf lane: run explicitly in release with --ignored"]
fn d3q19_sparse_sweep_glups() {
    let pre_axes = MachineAxes::probe();
    emit_observation(
        "axes/pre/measurement",
        "axes-pre",
        fs_obs::Severity::Info,
        format!(
            "{{\"metric\":\"axes-pre\",\"axes\":{}}}",
            pre_axes.to_jsonl()
        ),
    );
    if let Some(reason) = pre_axes.plausibility_error() {
        fail_invalid_environment("terminal/pre-axes/environment-invalid", reason);
    }
    let admission = prepare_admission(&pre_axes);
    emit_observation(
        "header/measurement",
        "d3q19-perf-header",
        fs_obs::Severity::Info,
        format!(
            "{{\"metric\":\"lbm-d3q19-header\",\"model\":{},\"warmup\":{WARMUP_REPETITIONS},\"repetitions\":{TIMED_REPETITIONS},\"floor_glups\":null}}",
            D3q19TrafficModel::default().receipt_json()
        ),
    );

    let workers = usize::try_from(pre_axes.logical_cpus)
        .expect("plausible logical CPU count must fit usize")
        .max(1);
    let lane_axes = [
        (OccupancyClass::DenseActive, ThreadingClass::SingleThread, 1),
        (
            OccupancyClass::DenseActive,
            ThreadingClass::AllCore,
            workers,
        ),
        (
            OccupancyClass::SparseTenPercent,
            ThreadingClass::SingleThread,
            1,
        ),
        (
            OccupancyClass::SparseTenPercent,
            ThreadingClass::AllCore,
            workers,
        ),
    ];
    let mut measured = Vec::new();
    measured
        .try_reserve_exact(lane_axes.len())
        .expect("bounded lane table allocation");
    for (occupancy, threading, lane_workers) in lane_axes {
        let shape = LaneShape::memory_resident(occupancy, threading, lane_workers)
            .expect("fixed lane shape admits");
        measured.push(measure_lane(shape, &pre_axes).unwrap_or_else(|reason| {
            fail_invalid_environment(
                &format!("terminal/{}/measurement-invalid", lane_identity(shape)),
                &reason,
            )
        }));
    }

    let post_axes = MachineAxes::probe();
    emit_observation(
        "axes/post/measurement",
        "axes-post",
        fs_obs::Severity::Info,
        format!(
            "{{\"metric\":\"axes-post\",\"axes\":{}}}",
            post_axes.to_jsonl()
        ),
    );
    let (snapshot, configuration_refusal) = admission.snapshot(&pre_axes, &post_axes);
    emit_observation(
        "admission/decision",
        "axis-admission",
        fs_obs::Severity::Info,
        snapshot.receipt_json().to_owned(),
    );
    let gate_admission =
        classify_gate_admission(&snapshot, configuration_refusal, &pre_axes, &post_axes);
    let family = reference_isa(&pre_axes);
    let mut invalid_reasons = Vec::new();

    for lane in measured {
        let identity = lane_identity(lane.shape);
        emit_observation(
            &format!("{identity}/roofline"),
            "roofline-attainment",
            fs_obs::Severity::Info,
            lane.attainment.to_jsonl(),
        );
        for (repetition, critical_path) in lane.critical_paths.iter().enumerate() {
            emit_observation(
                &format!("{identity}/critical-path/{repetition}"),
                "critical-path",
                fs_obs::Severity::Info,
                critical_path.receipt_json(),
            );
        }
        let evidence = row_evidence(&gate_admission, &lane.attainment);
        if let EvidenceClass::EnvironmentInvalid { reason } = &evidence {
            invalid_reasons.push(format!("{identity}: {reason}"));
        }
        let row = D3q19PerfRow {
            reference_isa: family,
            shape: lane.shape,
            glups: lane.attainment.elems_per_sec / 1.0e9,
            dispersion_ppm: dispersion_ppm(&lane.attainment),
            floor_glups: None,
            evidence,
            placement_identity: lane.placement_identity,
            critical_paths: lane.critical_paths,
        };
        let receipt = row.receipt_json().unwrap_or_else(|error| {
            fail_invalid_environment(
                &format!("terminal/{identity}/receipt-invalid"),
                &error.to_string(),
            )
        });
        emit_observation(
            &format!("{identity}/measurement"),
            "lbm-d3q19-sweep",
            fs_obs::Severity::Info,
            receipt,
        );
    }
    if !invalid_reasons.is_empty() {
        fail_invalid_environment(
            "terminal/post-axes/environment-invalid",
            &invalid_reasons.join("; "),
        );
    }
}
