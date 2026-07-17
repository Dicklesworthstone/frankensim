//! Ignored production-scale battery scaffold for `frankensim-ei3b`.
//!
//! The scalar-field tranche exercises the arena/operation-lease path at exact
//! 128^3 and 256^3 shapes. The sparse tranche compares one serial and one real
//! TilePool sweep over exactly 1,000,000 active D3Q19 cells. Both preserve
//! explicit no-claim boundaries for NUMA, CCD, quiet-host, cross-ISA, and
//! macOS peak-RSS acceptance that still need admitted host evidence.
//!
//! Run explicitly in release mode with one exact host profile:
//! `FRANKENSIM_PRODUCTION_SCALE_PROFILE=m4-128 cargo test --locked
//! -p fs-flagship-e2e --release --test production_scale -- --ignored
//! --nocapture`

use core::mem::size_of;
use fs_alloc::{
    AllocError, ArenaConfig, ArenaPool, ArenaStats, HUGEPAGE_BYTES, HugepagePolicy, LeaseReceipt,
    OperationMemoryLease, PoolStats, Site,
};
use fs_exec::{CancelGate, KernelRunner, RunError, RunReport, TileKernel, TilePool};
use fs_lbm::d3q19::{
    D3Q19_BIT_SEMANTICS_VERSION, Q3, TILE,
    sparse::{SparseGrid3, morton3, state_bytes_per_tile},
};
use fs_obs::{CapabilityDecision, EventKind, Severity};
use std::cell::{Cell, RefCell};
use std::time::{Duration, Instant};

const SESSION: &str = "fs-flagship-e2e/production-scale";
const PROFILE_ENV: &str = "FRANKENSIM_PRODUCTION_SCALE_PROFILE";
const FIELD_SITE: Site = Site::named("production-scale/scalar-field");
const REFUSAL_SITE: Site = Site::named("production-scale/refusal-probe");
const SPARSE_DOMAIN_EDGE: usize = 200;
const SPARSE_ACTIVE_TILE_ORIGIN: u32 = 12;
const SPARSE_ACTIVE_TILE_EDGE: u32 = 25;
const SPARSE_ACTIVE_TILES: usize = 15_625;
const SPARSE_ACTIVE_CELLS: usize = 1_000_000;
const SPARSE_TAU: f64 = 0.8;
const SPARSE_FORCE: [f64; 3] = [0.0; 3];
const SPARSE_PERTURB_SEED: u64 = 0x5ca1_e001_0000_0001;
const SPARSE_PERTURB_AMPLITUDE: f64 = 0.01;
const SPARSE_POOL_SEED: u64 = 0x71e0_0001_0000_0001;
const SPARSE_SWEEP_GROUP_TILES: usize = 8;
const SPARSE_STATE_LEASE_SITE: &str = "production-scale/sparse-d3q19-retained-state";

#[derive(Clone, Copy)]
struct Profile {
    name: &'static str,
    edge: usize,
    os: &'static str,
    arch: &'static str,
    model_needle: &'static str,
}

impl Profile {
    const M4_128: Self = Self {
        name: "m4-128",
        edge: 128,
        os: "macos",
        arch: "aarch64",
        model_needle: "Apple M4",
    };

    const THREADRIPPER_256: Self = Self {
        name: "threadripper-256",
        edge: 256,
        os: "linux",
        arch: "x86_64",
        model_needle: "Threadripper",
    };

    fn parse(value: &str) -> Result<Self, String> {
        match value {
            "m4-128" => Ok(Self::M4_128),
            "threadripper-256" => Ok(Self::THREADRIPPER_256),
            other => Err(format!(
                "unknown {PROFILE_ENV} value {other:?}; expected m4-128 or threadripper-256"
            )),
        }
    }

    fn cells(self) -> Result<usize, String> {
        self.edge
            .checked_mul(self.edge)
            .and_then(|area| area.checked_mul(self.edge))
            .ok_or_else(|| format!("{}^3 cell count overflows usize", self.edge))
    }
}

struct Evidence {
    emitter: fs_obs::Emitter,
}

struct RecordingRunner<'a> {
    pool: &'a TilePool,
    reports: RefCell<Vec<RunReport>>,
}

impl<'a> RecordingRunner<'a> {
    fn new(pool: &'a TilePool) -> Self {
        Self {
            pool,
            reports: RefCell::new(Vec::new()),
        }
    }

    fn reports(&self) -> Vec<RunReport> {
        self.reports.borrow().clone()
    }
}

impl KernelRunner for RecordingRunner<'_> {
    fn workers(&self) -> usize {
        self.pool.workers()
    }

    fn run_with_gate<K: TileKernel>(
        &self,
        kernel: &K,
        gate: &CancelGate,
    ) -> (Result<K::Out, RunError>, RunReport) {
        let (outcome, report) = self.pool.run_with_gate(kernel, gate);
        self.reports.borrow_mut().push(report.clone());
        (outcome, report)
    }
}

impl Evidence {
    fn new(scope: &str) -> Self {
        Self {
            emitter: fs_obs::Emitter::new(SESSION, scope),
        }
    }

    fn emit(&mut self, severity: Severity, kind: EventKind, wall_ns: Option<u64>) {
        let event = self.emitter.emit(severity, kind, wall_ns);
        fs_obs::lint_failure_record(&event)
            .expect("production-scale evidence must reproduce failures from the log");
        let line = event.to_jsonl();
        fs_obs::validate_line(&line)
            .expect("production-scale evidence must use the fs-obs wire schema");
        println!("{line}");
    }

    fn decision(
        &mut self,
        capability: &str,
        domain: &str,
        decision: CapabilityDecision,
        detail: impl Into<String>,
    ) {
        let severity = match decision {
            CapabilityDecision::Admitted => Severity::Info,
            CapabilityDecision::Refused | CapabilityDecision::Restricted => Severity::Warn,
        };
        self.emit(
            severity,
            EventKind::CapabilityDomainDecision {
                capability: capability.to_string(),
                domain: domain.to_string(),
                decision,
                detail: detail.into(),
            },
            None,
        );
    }

    fn verdict(&mut self, case: &str, pass: bool, detail: impl Into<String>) {
        let detail = detail.into();
        self.emit(
            if pass {
                Severity::Info
            } else {
                Severity::Error
            },
            EventKind::ConformanceCase {
                suite: SESSION.to_string(),
                case: case.to_string(),
                pass,
                detail: detail.clone(),
                seed: 0,
            },
            None,
        );
        assert!(pass, "case {case}: {detail}");
    }
}

struct RunEvidence {
    allocation: Duration,
    sweep: Duration,
    reclaim: Duration,
    checksum: u64,
    expected_checksum: u64,
    content_verified: bool,
    arena: ArenaStats,
    lease: LeaseReceipt,
    pool: PoolStats,
}

fn arena_config(reservation_bytes: usize) -> ArenaConfig {
    ArenaConfig {
        chunk_bytes: HUGEPAGE_BYTES,
        max_chunk_bytes: reservation_bytes,
        limit_bytes: Some(reservation_bytes),
        free_list_max_bytes: 0,
        hugepage: HugepagePolicy::Never,
    }
}

fn duration_ns(duration: Duration) -> u64 {
    u64::try_from(duration.as_nanos()).unwrap_or(u64::MAX)
}

fn field_value(index: usize) -> f64 {
    let lane = u32::try_from(index & 1023).expect("masked lane fits u32");
    f64::from(lane) * 0.125
}

fn checksum_step(checksum: u64, index: usize, value: f64) -> u64 {
    checksum.rotate_left(7)
        ^ value
            .to_bits()
            .wrapping_add(u64::try_from(index).expect("field index fits u64"))
}

fn run_field(cells: usize, reservation_bytes: usize) -> Result<RunEvidence, AllocError> {
    let pool = ArenaPool::new(arena_config(reservation_bytes));
    let lease = OperationMemoryLease::bounded(
        u64::try_from(reservation_bytes).expect("reservation bytes fit u64"),
    );
    let arena = pool.arena_leased(&lease);

    let allocation_start = Instant::now();
    let field = arena.alloc_slice_fill(FIELD_SITE, cells, 0.0f64)?;
    std::hint::black_box(&*field);
    let allocation = allocation_start.elapsed();

    let sweep_start = Instant::now();
    for (index, value) in field.iter_mut().enumerate() {
        *value = field_value(index);
    }
    std::hint::black_box(&*field);
    let sweep = sweep_start.elapsed();

    let mut checksum = 0xcbf2_9ce4_8422_2325u64;
    let mut expected_checksum = 0xcbf2_9ce4_8422_2325u64;
    let mut content_verified = true;
    for (index, value) in field.iter().copied().enumerate() {
        let expected = field_value(index);
        content_verified &= value.to_bits() == expected.to_bits();
        checksum = checksum_step(checksum, index, value);
        expected_checksum = checksum_step(expected_checksum, index, expected);
    }
    std::hint::black_box((checksum, expected_checksum, content_verified));
    let arena_stats = arena.stats();

    let reclaim_start = Instant::now();
    drop(arena);
    let reclaim = reclaim_start.elapsed();

    Ok(RunEvidence {
        allocation,
        sweep,
        reclaim,
        checksum,
        expected_checksum,
        content_verified,
        arena: arena_stats,
        lease: lease.receipt(),
        pool: pool.stats(),
    })
}

fn refusal_is_pre_mutation(cells: usize, reservation_bytes: usize) -> (bool, String) {
    let pool = ArenaPool::new(arena_config(reservation_bytes));
    let before = pool.stats();
    let reservation_u64 = u64::try_from(reservation_bytes).expect("reservation bytes fit u64");
    let refusal_limit = reservation_u64
        .checked_sub(1)
        .expect("non-empty field reservation exceeds zero");
    let lease = OperationMemoryLease::bounded(refusal_limit);
    let mutated = Cell::new(false);

    let (result, arena_stats) = pool.scope_leased(&lease, |arena| {
        match arena.alloc_slice_fill(REFUSAL_SITE, cells, 0.0f64) {
            Ok(field) => {
                field[0] = 1.0;
                mutated.set(true);
                let stats = arena.stats();
                (
                    Ok(()),
                    (
                        stats.allocated_bytes,
                        stats.allocation_count,
                        stats.chunk_count,
                    ),
                )
            }
            Err(error) => {
                let stats = arena.stats();
                (
                    Err(error),
                    (
                        stats.allocated_bytes,
                        stats.allocation_count,
                        stats.chunk_count,
                    ),
                )
            }
        }
    });
    let mutated = mutated.get();

    let after = pool.stats();
    let receipt = lease.receipt();
    let exact_error = matches!(
        &result,
        Err(AllocError::LeaseExhausted {
            site: "production-scale/refusal-probe",
            requested_bytes,
            used_bytes: 0,
            limit_bytes,
        }) if *requested_bytes == reservation_u64 && *limit_bytes == refusal_limit
    );
    let exact_refusal = receipt.first_refusal.as_ref().is_some_and(|refusal| {
        refusal.what == "arena-chunk"
            && refusal.requested_bytes == reservation_u64
            && refusal.used_bytes == 0
            && refusal.limit_bytes == refusal_limit
    });
    let pass = exact_error
        && exact_refusal
        && !mutated
        && arena_stats == (0, 0, 0)
        && before == after
        && receipt.limit_bytes == Some(refusal_limit)
        && receipt.requested_bytes == 0
        && receipt.peak_bytes == 0
        && receipt.used_bytes == 0
        && receipt.refusals == 1
        && receipt.release_invariant_violations == 0;
    let error = match result {
        Ok(()) => "allocation unexpectedly succeeded".to_string(),
        Err(error) => error.to_string(),
    };
    (
        pass,
        format!(
            "pre-mutation refusal: error={error}; mutated={mutated}; arena={arena_stats:?}; \
             before={}; after={}; lease={}",
            before.to_json(),
            after.to_json(),
            receipt.to_json()
        ),
    )
}

fn configured_profile() -> Result<Option<Profile>, String> {
    let Some(raw) = std::env::var_os(PROFILE_ENV) else {
        return Ok(None);
    };
    let value = raw
        .into_string()
        .map_err(|_| format!("{PROFILE_ENV} is not valid UTF-8"))?;
    Profile::parse(&value).map(Some)
}

fn command_value(program: &str, args: &[&str]) -> Option<String> {
    let output = std::process::Command::new(program)
        .args(args)
        .output()
        .ok()?;
    if !output.status.success() {
        return None;
    }
    let value = String::from_utf8(output.stdout).ok()?.trim().to_string();
    (!value.is_empty()).then_some(value)
}

fn machine_model() -> Result<String, String> {
    let model = if cfg!(target_os = "macos") {
        command_value("/usr/sbin/sysctl", &["-n", "machdep.cpu.brand_string"])
            .or_else(|| command_value("/usr/sbin/sysctl", &["-n", "hw.model"]))
    } else if cfg!(target_os = "linux") {
        std::fs::read_to_string("/proc/cpuinfo")
            .ok()
            .and_then(|cpuinfo| {
                cpuinfo.lines().find_map(|line| {
                    let (key, value) = line.split_once(':')?;
                    (key.trim() == "model name")
                        .then(|| value.trim().to_string())
                        .filter(|value| !value.is_empty())
                })
            })
    } else {
        None
    }
    .ok_or_else(|| "could not determine a bounded machine model string".to_string())?;
    if model.len() > 256 || model.chars().any(char::is_control) {
        return Err("machine model is empty, oversized, or contains control bytes".to_string());
    }
    Ok(model)
}

fn admitted_profile(evidence: &mut Evidence) -> Option<(Profile, String)> {
    let profile = match configured_profile() {
        Ok(Some(profile)) => profile,
        Ok(None) => {
            evidence.decision(
                "production-scale-profile",
                "ignored scale battery",
                CapabilityDecision::Refused,
                format!("named skip: set {PROFILE_ENV}=m4-128 or threadripper-256 explicitly"),
            );
            return None;
        }
        Err(error) => {
            evidence.decision(
                "production-scale-profile",
                "ignored scale battery",
                CapabilityDecision::Refused,
                error.clone(),
            );
            panic!("{error}");
        }
    };

    let model = match machine_model() {
        Ok(model) => model,
        Err(error) => {
            evidence.decision(
                "production-scale-host",
                profile.name,
                CapabilityDecision::Refused,
                error.clone(),
            );
            panic!("{error}");
        }
    };
    let host_matches = std::env::consts::OS == profile.os
        && std::env::consts::ARCH == profile.arch
        && model.contains(profile.model_needle);
    let assertions_off = !cfg!(debug_assertions);
    if !host_matches || !assertions_off {
        let detail = format!(
            "requested profile {} requires debug assertions off on {}+{} with model containing \
             {:?}; detected {}+{} model {:?}, debug_assertions={}",
            profile.name,
            profile.os,
            profile.arch,
            profile.model_needle,
            std::env::consts::OS,
            std::env::consts::ARCH,
            model,
            cfg!(debug_assertions)
        );
        evidence.decision(
            "production-scale-host",
            profile.name,
            CapabilityDecision::Refused,
            detail.clone(),
        );
        panic!("{detail}");
    }
    evidence.decision(
        "production-scale-host",
        profile.name,
        CapabilityDecision::Restricted,
        "OS, architecture, model family, and debug-assertions-off mode match; Cargo profile name, optimization level, quiet-host state, and topology authority remain unproven",
    );
    Some((profile, model))
}

fn linux_peak_rss_bytes() -> Option<u64> {
    if !cfg!(target_os = "linux") {
        return None;
    }
    let status = std::fs::read_to_string("/proc/self/status").ok()?;
    let kib = status.lines().find_map(|line| {
        let rest = line.strip_prefix("VmHWM:")?;
        rest.split_whitespace().next()?.parse::<u64>().ok()
    })?;
    kib.checked_mul(1024)
}

fn json_string(value: &str) -> String {
    use core::fmt::Write as _;

    let mut out = String::with_capacity(value.len() + 2);
    out.push('"');
    for character in value.chars() {
        match character {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            control if control.is_control() => {
                let _ = write!(out, "\\u{:04x}", u32::from(control));
            }
            other => out.push(other),
        }
    }
    out.push('"');
    out
}

fn json_f64(value: f64) -> String {
    if value.is_finite() {
        value.to_string()
    } else {
        "null".to_string()
    }
}

fn emit_process_peak_rss(
    evidence: &mut Evidence,
    profile: Profile,
    kernel: &str,
    configuration: u64,
    before: Option<u64>,
    after: Option<u64>,
) {
    match (before, after) {
        (Some(before), Some(after)) => {
            evidence.emit(
                Severity::Info,
                EventKind::Custom {
                    name: "production-scale-process-memory".to_string(),
                    json: format!(
                        "{{\"kernel\":{},\
                         \"metric\":\"linux_process_lifetime_peak_rss_bytes\",\
                         \"value_bytes\":{after},\"unit\":\"bytes\",\
                         \"configuration_identity_root\":{configuration},\
                         \"source\":\"/proc/self/status:VmHWM\",\
                         \"process_lifetime\":true,\"report_only\":true}}",
                        json_string(kernel)
                    ),
                },
                None,
            );
            evidence.decision(
                "peak-rss-budget",
                profile.name,
                CapabilityDecision::Restricted,
                format!(
                    "Linux VmHWM is an authoritative process-lifetime high-water, but includes \
                     harness startup and cannot be reset; before={before}, after={after}, \
                     delta={}",
                    after.saturating_sub(before)
                ),
            );
        }
        _ => evidence.decision(
            "peak-rss-budget",
            profile.name,
            CapabilityDecision::Refused,
            "named skip: no authoritative in-process peak-RSS source is available on this platform; current RSS is not relabeled as peak",
        ),
    }
}

fn emit_phase(
    evidence: &mut Evidence,
    kernel: &str,
    metric: &str,
    duration: Duration,
    configuration: u64,
) {
    let ns = duration_ns(duration);
    evidence.emit(
        Severity::Info,
        EventKind::Custom {
            name: "production-scale-phase".to_string(),
            json: format!(
                "{{\"kernel\":{},\"metric\":{},\"value_ns\":{ns},\
                 \"unit\":\"ns\",\"configuration_identity_root\":{configuration},\
                 \"report_only\":true}}",
                json_string(kernel),
                json_string(metric)
            ),
        },
        None,
    );
}

fn named_arena_open_claims(evidence: &mut Evidence, profile: Profile) {
    let domain = format!("{} scalar-field production-scale tranche", profile.name);
    for (capability, detail) in [
        (
            "quiet-host-performance-admission",
            "this tranche records report-only phase times; no quiet-host attestation or performance envelope is available",
        ),
        (
            "numa-first-touch-placement",
            "arena initialization touches the field serially; worker-owned NUMA first-touch is not implemented in this tranche",
        ),
        (
            "tile-pool-ownership",
            "the scalar-field tranche runs directly through one leased arena and does not claim TilePool scheduling or worker ownership",
        ),
        (
            "per-ccd-bandwidth-uniformity",
            "no measured per-CCD bandwidth surface is available through the current direct dependencies",
        ),
        (
            "ccd-shaped-reduction-tree",
            "flat and CCD-shaped million-cell reductions remain a separate measured tranche",
        ),
        (
            "sparse-lbm-million-active-cells",
            "this tranche proves scalar-field arena and lease scale only; sparse D3Q19 execution remains open",
        ),
        (
            "cross-isa-bit-identity",
            "no paired M4 and Threadripper retained run exists for this scaffold",
        ),
        (
            "machine-fingerprint-authority",
            "phase rows bind a workload configuration identity, not an fs-substrate machine fingerprint; cross-run performance comparison remains refused",
        ),
        (
            "full-production-scale-acceptance",
            "the arena/lease tranche cannot close the production-scale battery while RSS gating, NUMA, CCD, and retained dual-host evidence remain open; sparse execution belongs to its companion tranche",
        ),
    ] {
        evidence.decision(capability, &domain, CapabilityDecision::Refused, detail);
    }
}

#[test]
#[ignore = "production-scale lane: requires debug assertions off and an exact admitted host profile"]
#[allow(clippy::too_many_lines)]
fn production_scale_arena_budget() {
    let mut evidence = Evidence::new("arena-budget");
    let Some((profile, model)) = admitted_profile(&mut evidence) else {
        return;
    };

    let cells = profile.cells().expect("profile cell count must fit");
    let payload_bytes = cells
        .checked_mul(size_of::<f64>())
        .expect("profile payload bytes must fit");
    let preflight_pool = ArenaPool::new(arena_config(payload_bytes));
    let reservation_bytes = preflight_pool
        .reservation_bytes_for_slice::<f64>(FIELD_SITE, cells)
        .expect("profile reservation must be representable");
    let reservation_u64 = u64::try_from(reservation_bytes).expect("reservation fits u64");
    let logical_cpus = std::thread::available_parallelism().map_or(1, usize::from);
    let identity = fs_obs::ident::IdentityBuilder::new("production-scale-arena-config-v1")
        .str("profile", profile.name)
        .str("build-mode", "debug-assertions-off-profile-unattested")
        .u64("edge", u64::try_from(profile.edge).expect("edge fits u64"))
        .u64("cells", u64::try_from(cells).expect("cells fit u64"))
        .u64(
            "element-bytes",
            u64::try_from(size_of::<f64>()).expect("element size fits u64"),
        )
        .u64(
            "payload-bytes",
            u64::try_from(payload_bytes).expect("payload fits u64"),
        )
        .u64("reservation-bytes", reservation_u64)
        .u64("pool-limit-bytes", reservation_u64)
        .u64("lease-limit-bytes", reservation_u64)
        .u64(
            "chunk-bytes",
            u64::try_from(HUGEPAGE_BYTES).expect("chunk fits u64"),
        )
        .u64("max-chunk-bytes", reservation_u64)
        .u64("free-list-max-bytes", 0)
        .str("hugepage-policy", "never")
        .str("initialization", "serial-arena-fill-f64-zero")
        .str("sweep", "serial-lane-mod-1024-times-one-eighth-v1")
        .str("os", std::env::consts::OS)
        .str("arch", std::env::consts::ARCH)
        .str("model", &model)
        .u64(
            "logical-cpus",
            u64::try_from(logical_cpus).expect("logical CPU count fits u64"),
        )
        .str("fs-alloc-version", fs_alloc::VERSION)
        .str("fs-obs-version", fs_obs::VERSION)
        .str("harness-version", env!("CARGO_PKG_VERSION"))
        .exclude(
            "phase-wall-ns",
            "report-only measurement, not configuration identity",
        )
        .exclude("process-rss", "run-varying process observation")
        .finish();
    evidence.emit(
        Severity::Info,
        EventKind::Custom {
            name: "production-scale-configuration".to_string(),
            json: format!(
                "{{\"profile\":{},\"build_mode\":\"debug-assertions-off-profile-unattested\",\
                 \"edge\":{},\"cells\":{cells},\"payload_bytes\":{payload_bytes},\
                 \"reservation_bytes\":{reservation_bytes},\"machine_model\":{},\
                 \"logical_cpus\":{logical_cpus},\"replay_identity\":{},\
                 \"replay_identity_root\":{},\"identity_exclusions\":2}}",
                json_string(profile.name),
                profile.edge,
                json_string(&model),
                json_string(&identity.hex()),
                identity.root()
            ),
        },
        None,
    );

    let (refusal_pass, refusal_detail) = refusal_is_pre_mutation(cells, reservation_bytes);
    evidence.verdict(
        "scale-arena-refusal-pre-mutation",
        refusal_pass,
        refusal_detail,
    );

    let rss_before = linux_peak_rss_bytes();
    let run = match run_field(cells, reservation_bytes) {
        Ok(run) => run,
        Err(error) => {
            evidence.verdict(
                "scale-arena-allocate-touch-reclaim",
                false,
                format!("admitted scale allocation refused unexpectedly: {error}"),
            );
            unreachable!("failing verdict asserts")
        }
    };
    let rss_after = linux_peak_rss_bytes();

    emit_phase(
        &mut evidence,
        "production-scale-arena",
        "allocate_initialize_first_touch_ns",
        run.allocation,
        identity.root(),
    );
    emit_phase(
        &mut evidence,
        "production-scale-arena",
        "serial_sweep_ns",
        run.sweep,
        identity.root(),
    );
    emit_phase(
        &mut evidence,
        "production-scale-arena",
        "arena_drop_reclaim_ns",
        run.reclaim,
        identity.root(),
    );

    emit_process_peak_rss(
        &mut evidence,
        profile,
        "production-scale-arena",
        identity.root(),
        rss_before,
        rss_after,
    );

    let payload_u64 = u64::try_from(payload_bytes).expect("payload bytes fit u64");
    let scale_pass = run.arena.allocated_bytes == payload_u64
        && run.content_verified
        && run.checksum == run.expected_checksum
        && run.arena.allocation_count == 1
        && run.arena.chunk_count == 1
        && run.lease.limit_bytes == Some(reservation_u64)
        && run.lease.requested_bytes == reservation_u64
        && run.lease.peak_bytes == reservation_u64
        && run.lease.used_bytes == 0
        && run.lease.refusals == 0
        && run.lease.release_invariant_violations == 0
        && run.pool.quiescent()
        && run.pool.arenas_live == 0
        && run.pool.reserved_bytes == 0
        && run.pool.free_bytes == 0
        && run.pool.chunks_created == 1
        && run.pool.chunks_recycled == 0;
    evidence.emit(
        if scale_pass {
            Severity::Info
        } else {
            Severity::Error
        },
        EventKind::Custom {
            name: "production-scale-arena-receipt".to_string(),
            json: format!(
                "{{\"profile\":{},\"cells\":{cells},\"payload_bytes\":{payload_bytes},\
                 \"reservation_bytes\":{reservation_bytes},\"checksum\":{},\
                 \"expected_checksum\":{},\"content_verified\":{},\"arena\":{},\
                 \"lease\":{},\"pool\":{},\"logical_budget_pass\":{scale_pass}}}",
                json_string(profile.name),
                run.checksum,
                run.expected_checksum,
                run.content_verified,
                run.arena.to_json(),
                run.lease.to_json(),
                run.pool.to_json()
            ),
        },
        None,
    );
    named_arena_open_claims(&mut evidence, profile);
    evidence.verdict(
        "scale-arena-allocate-touch-reclaim",
        scale_pass,
        format!(
            "{}^3 field: payload={payload_bytes} B, reservation={reservation_bytes} B, \
             checksum={:#018x}, expected_checksum={:#018x}, content_verified={}; \
             arena={}; lease={}; pool={}",
            profile.edge,
            run.checksum,
            run.expected_checksum,
            run.content_verified,
            run.arena.to_json(),
            run.lease.to_json(),
            run.pool.to_json()
        ),
    );
}

fn sparse_active_tiles() -> Vec<(u32, u32, u32)> {
    let end = SPARSE_ACTIVE_TILE_ORIGIN
        .checked_add(SPARSE_ACTIVE_TILE_EDGE)
        .expect("sparse active-tile range fits u32");
    (SPARSE_ACTIVE_TILE_ORIGIN..end)
        .flat_map(|tx| {
            (SPARSE_ACTIVE_TILE_ORIGIN..end)
                .flat_map(move |ty| (SPARSE_ACTIVE_TILE_ORIGIN..end).map(move |tz| (tx, ty, tz)))
        })
        .collect()
}

fn sparse_mass_roundoff_bound(population_values: usize) -> f64 {
    let population_values =
        u32::try_from(population_values).expect("production-scale population count fits u32");
    let scaled_epsilon = f64::from(population_values) * f64::EPSILON;
    assert!(
        scaled_epsilon < 1.0,
        "roundoff-bound denominator must stay positive"
    );
    8.0 * scaled_epsilon / (1.0 - scaled_epsilon)
}

fn named_sparse_open_claims(evidence: &mut Evidence, profile: Profile) {
    let domain = format!("{} sparse D3Q19 million-cell acceptance", profile.name);
    for (capability, detail) in [
        (
            "quiet-host-performance-admission",
            "serial and pooled phase times are report-only; no quiet-host attestation, GLUP/s floor, or roofline authority is present",
        ),
        (
            "numa-first-touch-placement",
            "SparseGrid3 activation and perturbation first-touch ordinary heap vectors serially before the pooled sweep; worker-owned placement is not implemented",
        ),
        (
            "per-ccd-bandwidth-uniformity",
            "the TilePool placement identity is retained, but no per-CCD byte counter, pin-success receipt, or bandwidth surface is measured",
        ),
        (
            "ccd-shaped-reduction-tree",
            "the million-cell mass oracle uses the grid's canonical flat reduction; a separately identified CCD-shaped reduction remains open",
        ),
        (
            "million-cell-mid-stream-cancellation",
            "this tranche runs one open-gate sweep; the scale cancellation-latency bead owns timed mid-stream request, drain, and finalize evidence",
        ),
        (
            "million-cell-worker-count-bit-identity",
            "pooled output is compared bit-for-bit with one serial reference, but no second pooled worker-count run is retained at this scale",
        ),
        (
            "cross-isa-bit-identity",
            "no paired retained M4 and Threadripper output exists for this production-scale protocol",
        ),
        (
            "machine-fingerprint-authority",
            "phase rows bind the workload and TilePool placement identities, not an admitted fs-substrate machine fingerprint; cross-run performance comparison remains refused",
        ),
        (
            "full-production-scale-acceptance",
            "million-cell sparse execution alone cannot close RSS-budget, NUMA, CCD, cancellation-latency, and retained dual-host evidence",
        ),
    ] {
        evidence.decision(capability, &domain, CapabilityDecision::Refused, detail);
    }
}

#[test]
#[ignore = "production-scale lane: requires debug assertions off and an exact admitted host profile"]
#[allow(clippy::too_many_lines)]
fn production_scale_sparse_lbm_million_cells() {
    let mut evidence = Evidence::new("sparse-lbm-million");
    let Some((profile, model)) = admitted_profile(&mut evidence) else {
        return;
    };

    let tiles = sparse_active_tiles();
    let mut reverse_tiles = tiles.clone();
    reverse_tiles.reverse();
    let tile_cells = TILE
        .checked_mul(TILE)
        .and_then(|area| area.checked_mul(TILE))
        .expect("sparse tile cell count fits usize");
    let active_cells = tiles
        .len()
        .checked_mul(tile_cells)
        .expect("sparse active-cell count fits usize");
    let bytes_per_tile = state_bytes_per_tile();
    let retained_state_bytes = tiles
        .len()
        .checked_mul(bytes_per_tile)
        .expect("sparse retained-state bytes fit usize");
    let published_values = active_cells
        .checked_mul(Q3)
        .expect("sparse published-population count fits usize");
    let published_bytes = published_values
        .checked_mul(size_of::<f64>())
        .expect("sparse published-population bytes fit usize");
    let harness_retained_state_bytes = retained_state_bytes
        .checked_mul(2)
        .expect("serial plus pooled retained-state bytes fit usize");
    let oracle_copy_bytes = published_bytes
        .checked_mul(2)
        .expect("serial plus pooled oracle-copy bytes fit usize");
    let expected_groups = tiles.len().div_ceil(SPARSE_SWEEP_GROUP_TILES);
    let mass_bound = sparse_mass_roundoff_bound(published_values);
    let tile_domain_edge = SPARSE_DOMAIN_EDGE / TILE;
    let active_end = SPARSE_ACTIVE_TILE_ORIGIN + SPARSE_ACTIVE_TILE_EDGE;

    let mut active_keys: Vec<u64> = tiles
        .iter()
        .map(|&(tx, ty, tz)| morton3(tx, ty, tz))
        .collect();
    active_keys.sort_unstable();
    let active_set_identity = active_keys
        .into_iter()
        .fold(
            fs_obs::ident::IdentityBuilder::new("sparse-d3q19-active-key-set-v1"),
            |builder, key| builder.u64("morton-key", key),
        )
        .finish();

    let logical_cpus = std::thread::available_parallelism().map_or(1, usize::from);
    let pool = TilePool::for_host(logical_cpus, SPARSE_POOL_SEED);
    let placement_identity = pool.placement_identity();
    let placement_round_trips = pool
        .admit_retained_placement_identity(
            TilePool::PLACEMENT_IDENTITY_VERSION,
            &placement_identity,
        )
        .is_ok();
    let identity = fs_obs::ident::IdentityBuilder::new("production-scale-sparse-d3q19-v1")
        .str("profile", profile.name)
        .str("build-mode", "debug-assertions-off-profile-unattested")
        .str("os", std::env::consts::OS)
        .str("arch", std::env::consts::ARCH)
        .str("model", &model)
        .u64(
            "logical-cpus",
            u64::try_from(logical_cpus).expect("logical CPU count fits u64"),
        )
        .u64(
            "domain-edge-cells",
            u64::try_from(SPARSE_DOMAIN_EDGE).expect("domain edge fits u64"),
        )
        .u64(
            "tile-edge-cells",
            u64::try_from(TILE).expect("tile edge fits u64"),
        )
        .u64(
            "tile-cells",
            u64::try_from(tile_cells).expect("tile cells fit u64"),
        )
        .u64("active-tile-origin", u64::from(SPARSE_ACTIVE_TILE_ORIGIN))
        .u64("active-tile-edge", u64::from(SPARSE_ACTIVE_TILE_EDGE))
        .u64(
            "active-tiles",
            u64::try_from(tiles.len()).expect("active tiles fit u64"),
        )
        .u64(
            "active-cells",
            u64::try_from(active_cells).expect("active cells fit u64"),
        )
        .u64(
            "populations-per-cell",
            u64::try_from(Q3).expect("Q fits u64"),
        )
        .u64(
            "population-element-bytes",
            u64::try_from(size_of::<f64>()).expect("f64 size fits u64"),
        )
        .u64("retained-population-buffers", 3)
        .u64(
            "retained-state-bytes-per-tile",
            u64::try_from(bytes_per_tile).expect("per-tile bytes fit u64"),
        )
        .u64(
            "workload-retained-state-bytes",
            u64::try_from(retained_state_bytes).expect("retained bytes fit u64"),
        )
        .u64(
            "serial-plus-pooled-retained-state-bytes",
            u64::try_from(harness_retained_state_bytes).expect("harness bytes fit u64"),
        )
        .u64(
            "oracle-copy-bytes",
            u64::try_from(oracle_copy_bytes).expect("oracle bytes fit u64"),
        )
        .str("collision-model", "D3Q19-BGK")
        .f64_bits("tau", SPARSE_TAU)
        .f64_bits("force-x", SPARSE_FORCE[0])
        .f64_bits("force-y", SPARSE_FORCE[1])
        .f64_bits("force-z", SPARSE_FORCE[2])
        .u64("perturb-seed", SPARSE_PERTURB_SEED)
        .f64_bits("perturb-amplitude", SPARSE_PERTURB_AMPLITUDE)
        .u64("pooled-run-seed", SPARSE_POOL_SEED)
        .u64("declared-run-id", 0)
        .u64("pooled-steps", 1)
        .str("serial-activation-order", "ascending-coordinate-input")
        .str("pooled-activation-order", "descending-coordinate-input")
        .str("canonical-active-order", "ascending-morton")
        .child("active-key-set", &active_set_identity)
        .str("serial-step", "SparseGrid3::step_serial-one-step")
        .str("pooled-step", "SparseGrid3::step_pooled-two-pass-v1")
        .str("state-oracle", "full-u64-bit-vector-equality")
        .str("mass-bound-policy", "sparse-mass-roundoff-envelope-v1")
        .str("mass-bound-formula", "8*n*epsilon/(1-n*epsilon)")
        .u64("mass-bound-multiplier", 8)
        .u64(
            "mass-reduction-population-values",
            u64::try_from(published_values).expect("population values fit u64"),
        )
        .f64_bits("mass-bound", mass_bound)
        .u64(
            "sweep-group-tiles",
            u64::try_from(SPARSE_SWEEP_GROUP_TILES).expect("group tiles fit u64"),
        )
        .u64(
            "expected-kernel-groups",
            u64::try_from(expected_groups).expect("group count fits u64"),
        )
        .u64(
            "pool-workers",
            u64::try_from(pool.workers()).expect("pool workers fit u64"),
        )
        .u64(
            "pool-placement-identity-version",
            u64::from(TilePool::PLACEMENT_IDENTITY_VERSION),
        )
        .str("pool-placement-identity", &placement_identity)
        .u64(
            "d3q19-bit-semantics-version",
            u64::from(D3Q19_BIT_SEMANTICS_VERSION),
        )
        .str("fs-exec-version", fs_exec::VERSION)
        .str("harness-version", env!("CARGO_PKG_VERSION"))
        .str(
            "memory-admission",
            "coarse-caller-held-two-grid-retained-population-state-only",
        )
        .exclude(
            "phase-wall-ns",
            "report-only measurement, not configuration identity",
        )
        .exclude("process-rss", "run-varying process observation")
        .exclude(
            "allocator-metadata",
            "ordinary Vec and BTreeMap implementation overhead is not surfaced",
        )
        .exclude(
            "activation-temporaries",
            "ordinary heap construction overhead is not surfaced",
        )
        .exclude(
            "pin-success",
            "observed timing fact, not requested placement identity",
        )
        .finish();

    evidence.emit(
        Severity::Info,
        EventKind::Custom {
            name: "production-scale-sparse-configuration".to_string(),
            json: format!(
                "{{\"profile\":{profile_json},\
                 \"build_mode\":\"debug-assertions-off-profile-unattested\",\
                 \"os\":{os_json},\
                 \"arch\":{arch_json},\"machine_model\":{model_json},\
                 \"domain_cells\":[{domain_edge},{domain_edge},{domain_edge}],\
                 \"tile_edge\":{tile_edge},\"active_cube_origin\":{active_origin},\
                 \"active_cube_edge_tiles\":{active_edge},\"active_tiles\":{active_tiles},\
                 \"active_cells\":{active_cells},\"bytes_per_tile\":{bytes_per_tile},\
                 \"workload_retained_state_bytes\":{retained_state_bytes},\
                 \"harness_retained_state_bytes\":{harness_retained_state_bytes},\
                 \"oracle_copy_bytes\":{oracle_copy_bytes},\
                 \"collision_model\":\"D3Q19-BGK\",\"tau_bits\":{tau_bits_json},\
                 \"force_bits\":[{force_x_json},{force_y_json},{force_z_json}],\
                 \"perturb_seed\":{perturb_seed_json},\
                 \"perturb_amplitude_bits\":{perturb_amplitude_json},\
                 \"pool_seed\":{pool_seed_json},\"steps\":1,\
                 \"expected_kernel_groups\":{expected_groups},\
                 \"mass_bound_policy\":\"sparse-mass-roundoff-envelope-v1\",\
                 \"mass_bound_formula\":\"8*n*epsilon/(1-n*epsilon)\",\
                 \"mass_bound_multiplier\":8,\"mass_bound\":{mass_bound_json},\
                 \"mass_bound_bits\":{mass_bound_bits_json},\
                 \"workers\":{logical_cpus},\
                 \"placement_identity\":{placement_json},\
                 \"active_key_set_identity\":{active_set_json},\
                 \"active_key_set_identity_root\":{active_set_root},\
                 \"replay_identity\":{identity_json},\
                 \"replay_identity_root\":{identity_root},\"identity_exclusions\":5}}",
                profile_json = json_string(profile.name),
                os_json = json_string(std::env::consts::OS),
                arch_json = json_string(std::env::consts::ARCH),
                model_json = json_string(&model),
                domain_edge = SPARSE_DOMAIN_EDGE,
                tile_edge = TILE,
                active_origin = SPARSE_ACTIVE_TILE_ORIGIN,
                active_edge = SPARSE_ACTIVE_TILE_EDGE,
                active_tiles = tiles.len(),
                tau_bits_json = json_string(&format!("0x{:016x}", SPARSE_TAU.to_bits())),
                force_x_json = json_string(&format!("0x{:016x}", SPARSE_FORCE[0].to_bits())),
                force_y_json = json_string(&format!("0x{:016x}", SPARSE_FORCE[1].to_bits())),
                force_z_json = json_string(&format!("0x{:016x}", SPARSE_FORCE[2].to_bits())),
                perturb_seed_json = json_string(&format!("0x{SPARSE_PERTURB_SEED:016x}")),
                perturb_amplitude_json =
                    json_string(&format!("0x{:016x}", SPARSE_PERTURB_AMPLITUDE.to_bits())),
                pool_seed_json = json_string(&format!("0x{SPARSE_POOL_SEED:016x}")),
                mass_bound_json = json_f64(mass_bound),
                mass_bound_bits_json = json_string(&format!("0x{:016x}", mass_bound.to_bits())),
                placement_json = json_string(&placement_identity),
                active_set_json = json_string(&active_set_identity.hex()),
                active_set_root = active_set_identity.root(),
                identity_json = json_string(&identity.hex()),
                identity_root = identity.root()
            ),
        },
        None,
    );

    let sparse_domain = format!("{} sparse D3Q19 million-cell acceptance", profile.name);
    evidence.decision(
        "sparse-state-memory-lease-authority",
        &sparse_domain,
        CapabilityDecision::Refused,
        "SparseGrid3 and step_pooled expose no lease-aware allocation path; the shadow receipt has zero allocation authority",
    );
    evidence.decision(
        "structured-sparse-heap-oom-refusal",
        &sparse_domain,
        CapabilityDecision::Refused,
        "SparseGrid3 retains ordinary infallible Vec/BTreeMap storage, so allocator pressure may abort before a typed refusal; this explicit ignored host lane does not claim OOM containment",
    );
    evidence.decision(
        "shadow-memory-preflight-ledger",
        &sparse_domain,
        CapabilityDecision::Restricted,
        format!(
            "the planned {harness_retained_state_bytes}-byte charge is detached accounting for two exact {retained_state_bytes}-byte retained population payloads only; actual grid, oracle, and TilePool allocations are not charged"
        ),
    );

    let rss_before = linux_peak_rss_bytes();
    let retained_lease_bytes =
        u64::try_from(harness_retained_state_bytes).expect("harness bytes fit u64");
    let retained_lease = OperationMemoryLease::bounded(retained_lease_bytes);
    let retained_charge =
        match retained_lease.reserve(SPARSE_STATE_LEASE_SITE, retained_lease_bytes) {
            Ok(charge) => charge,
            Err(error) => {
                evidence.verdict(
                    "scale-sparse-lbm-million-active",
                    false,
                    format!("coarse retained-state preflight refused unexpectedly: {error:?}"),
                );
                unreachable!("failing verdict asserts")
            }
        };

    let construct_start = Instant::now();
    let mut serial = match SparseGrid3::new(
        SPARSE_DOMAIN_EDGE,
        SPARSE_DOMAIN_EDGE,
        SPARSE_DOMAIN_EDGE,
        SPARSE_TAU,
        SPARSE_FORCE,
    ) {
        Ok(grid) => grid,
        Err(error) => {
            evidence.verdict(
                "scale-sparse-lbm-million-active",
                false,
                format!("serial sparse-grid construction refused: {error}"),
            );
            unreachable!("failing verdict asserts")
        }
    };
    let mut pooled = match SparseGrid3::new(
        SPARSE_DOMAIN_EDGE,
        SPARSE_DOMAIN_EDGE,
        SPARSE_DOMAIN_EDGE,
        SPARSE_TAU,
        SPARSE_FORCE,
    ) {
        Ok(grid) => grid,
        Err(error) => {
            evidence.verdict(
                "scale-sparse-lbm-million-active",
                false,
                format!("pooled sparse-grid construction refused: {error}"),
            );
            unreachable!("failing verdict asserts")
        }
    };
    let construct = construct_start.elapsed();
    emit_phase(
        &mut evidence,
        "sparse-d3q19-million-active",
        "construct_two_empty_grids_ns",
        construct,
        identity.root(),
    );

    let serial_activation_start = Instant::now();
    if let Err(error) = serial.activate_tiles(&tiles) {
        evidence.verdict(
            "scale-sparse-lbm-million-active",
            false,
            format!("serial million-cell activation refused: {error}"),
        );
    }
    let serial_activation = serial_activation_start.elapsed();
    emit_phase(
        &mut evidence,
        "sparse-d3q19-million-active",
        "serial_activate_equilibrium_ns",
        serial_activation,
        identity.root(),
    );

    let pooled_activation_start = Instant::now();
    if let Err(error) = pooled.activate_tiles(&reverse_tiles) {
        evidence.verdict(
            "scale-sparse-lbm-million-active",
            false,
            format!("pooled million-cell activation refused: {error}"),
        );
    }
    let pooled_activation = pooled_activation_start.elapsed();
    emit_phase(
        &mut evidence,
        "sparse-d3q19-million-active",
        "pooled_activate_equilibrium_ns",
        pooled_activation,
        identity.root(),
    );

    let perturb_start = Instant::now();
    serial.perturb(SPARSE_PERTURB_SEED, SPARSE_PERTURB_AMPLITUDE);
    pooled.perturb(SPARSE_PERTURB_SEED, SPARSE_PERTURB_AMPLITUDE);
    let perturb = perturb_start.elapsed();
    emit_phase(
        &mut evidence,
        "sparse-d3q19-million-active",
        "perturb_two_grids_ns",
        perturb,
        identity.root(),
    );

    let serial_initial_mass = serial.total_mass();
    let pooled_initial_mass = pooled.total_mass();
    let serial_state_bytes = serial.allocated_state_bytes();
    let pooled_state_bytes = pooled.allocated_state_bytes();
    let layout_pass = tiles.len() == SPARSE_ACTIVE_TILES
        && tile_cells == 64
        && active_cells == SPARSE_ACTIVE_CELLS
        && bytes_per_tile == 3 * Q3 * tile_cells * size_of::<f64>()
        && retained_state_bytes == 456_000_000
        && serial.active_tiles() == SPARSE_ACTIVE_TILES
        && pooled.active_tiles() == SPARSE_ACTIVE_TILES
        && serial_state_bytes == retained_state_bytes
        && pooled_state_bytes == retained_state_bytes
        && tiles
            .iter()
            .all(|&(tx, ty, tz)| serial.is_active(tx, ty, tz) && pooled.is_active(tx, ty, tz))
        && !serial.is_active(
            SPARSE_ACTIVE_TILE_ORIGIN - 1,
            SPARSE_ACTIVE_TILE_ORIGIN,
            SPARSE_ACTIVE_TILE_ORIGIN,
        )
        && !pooled.is_active(
            active_end,
            SPARSE_ACTIVE_TILE_ORIGIN,
            SPARSE_ACTIVE_TILE_ORIGIN,
        )
        && usize::try_from(active_end).is_ok_and(|end| end < tile_domain_edge)
        && serial_initial_mass.to_bits() == pooled_initial_mass.to_bits();
    evidence.emit(
        if layout_pass {
            Severity::Info
        } else {
            Severity::Error
        },
        EventKind::Custom {
            name: "production-scale-sparse-layout".to_string(),
            json: format!(
                "{{\"requested_tiles\":{},\"serial_active_tiles\":{},\
                 \"pooled_active_tiles\":{},\"active_cells\":{active_cells},\
                 \"bytes_per_tile\":{bytes_per_tile},\
                 \"serial_retained_state_bytes\":{serial_state_bytes},\
                 \"pooled_retained_state_bytes\":{pooled_state_bytes},\
                 \"serial_initial_mass\":{},\"pooled_initial_mass\":{},\
                 \"layout_pass\":{layout_pass}}}",
                tiles.len(),
                serial.active_tiles(),
                pooled.active_tiles(),
                json_f64(serial_initial_mass),
                json_f64(pooled_initial_mass)
            ),
        },
        None,
    );
    if !layout_pass {
        evidence.verdict(
            "scale-sparse-lbm-million-active",
            false,
            format!(
                "layout mismatch: requested={}, serial_active={}, pooled_active={}, active_cells={active_cells}, bytes_per_tile={bytes_per_tile}, serial_state={serial_state_bytes}, pooled_state={pooled_state_bytes}, serial_mass={serial_initial_mass:.17e}, pooled_mass={pooled_initial_mass:.17e}",
                tiles.len(),
                serial.active_tiles(),
                pooled.active_tiles()
            ),
        );
    }

    let serial_step_start = Instant::now();
    let serial_step = serial.step_serial();
    let serial_step_time = serial_step_start.elapsed();
    emit_phase(
        &mut evidence,
        "sparse-d3q19-million-active",
        "serial_reference_step_ns",
        serial_step_time,
        identity.root(),
    );
    if let Err(error) = serial_step {
        evidence.verdict(
            "scale-sparse-lbm-million-active",
            false,
            format!("serial million-cell reference sweep refused: {error}"),
        );
    }

    let gate = CancelGate::new();
    let runner = RecordingRunner::new(&pool);
    let pooled_step_start = Instant::now();
    let pooled_step = pooled.step_pooled(&runner, &gate);
    let pooled_step_time = pooled_step_start.elapsed();
    emit_phase(
        &mut evidence,
        "sparse-d3q19-million-active",
        "tilepool_step_ns",
        pooled_step_time,
        identity.root(),
    );
    if let Err(error) = pooled_step {
        evidence.verdict(
            "scale-sparse-lbm-million-active",
            false,
            format!("pooled million-cell sweep refused: {error}"),
        );
    }
    let reports = runner.reports();

    let oracle_start = Instant::now();
    let serial_bits = serial.state_bits();
    let pooled_bits = pooled.state_bits();
    let first_divergence = serial_bits
        .iter()
        .zip(&pooled_bits)
        .position(|(serial, pooled)| serial != pooled)
        .or_else(|| {
            (serial_bits.len() != pooled_bits.len())
                .then_some(serial_bits.len().min(pooled_bits.len()))
        });
    let state_bits_equal = first_divergence.is_none();
    let population_values_exact = pooled_bits.len() == published_values;
    let populations_finite = pooled_bits
        .iter()
        .all(|bits| f64::from_bits(*bits).is_finite());
    drop(serial_bits);
    drop(pooled_bits);
    let oracle_time = oracle_start.elapsed();
    emit_phase(
        &mut evidence,
        "sparse-d3q19-million-active",
        "exact_full_state_oracle_ns",
        oracle_time,
        identity.root(),
    );

    let serial_final_mass = serial.total_mass();
    let pooled_final_mass = pooled.total_mass();
    let mass_scale = serial_initial_mass.abs().max(1.0);
    let serial_mass_residual = (serial_final_mass - serial_initial_mass).abs() / mass_scale;
    let pooled_mass_residual = (pooled_final_mass - pooled_initial_mass).abs() / mass_scale;
    let mass_pass = serial_initial_mass.is_finite()
        && serial_final_mass.is_finite()
        && pooled_initial_mass.is_finite()
        && pooled_final_mass.is_finite()
        && serial_mass_residual <= mass_bound
        && pooled_mass_residual <= mass_bound
        && serial_final_mass.to_bits() == pooled_final_mass.to_bits();
    let report_pass = reports.len() == 2
        && reports[0].kernel == "fs-lbm/d3q19-sparse-collide"
        && reports[1].kernel == "fs-lbm/d3q19-sparse-stream"
        && reports.iter().all(|report| {
            report.mode == "deterministic"
                && report.declared_run.0 == 0
                && report.completed == u64::try_from(expected_groups).expect("groups fit u64")
                && report.total == u64::try_from(expected_groups).expect("groups fit u64")
                && report.cancel_latencies_ns.is_empty()
                && report.tiles_by_worker.len() == pool.workers()
                && report.tiles_by_worker.iter().sum::<u64>()
                    == u64::try_from(expected_groups).expect("groups fit u64")
        });
    let pool_stats = pool.arena_pool().stats();
    let lease_while_live = retained_lease.receipt();
    let lease_live_pass = retained_charge.bytes() == retained_lease_bytes
        && lease_while_live.limit_bytes == Some(retained_lease_bytes)
        && lease_while_live.requested_bytes == retained_lease_bytes
        && lease_while_live.peak_bytes == retained_lease_bytes
        && lease_while_live.used_bytes == retained_lease_bytes
        && lease_while_live.refusals == 0
        && lease_while_live.first_refusal.is_none()
        && lease_while_live.release_invariant_violations == 0;
    let step_pass = layout_pass
        && serial.steps() == 1
        && pooled.steps() == 1
        && serial.allocated_state_bytes() == retained_state_bytes
        && pooled.allocated_state_bytes() == retained_state_bytes
        && state_bits_equal
        && population_values_exact
        && populations_finite
        && mass_pass
        && report_pass
        && !gate.is_requested()
        && placement_round_trips
        && pool_stats.quiescent()
        && pool_stats.arenas_live == 0
        && lease_live_pass;

    let reclaim_start = Instant::now();
    drop(serial);
    drop(pooled);
    drop(retained_charge);
    let reclaim = reclaim_start.elapsed();
    let lease_after_drop = retained_lease.receipt();
    let release_pass = lease_after_drop.limit_bytes == Some(retained_lease_bytes)
        && lease_after_drop.used_bytes == 0
        && lease_after_drop.requested_bytes == retained_lease_bytes
        && lease_after_drop.peak_bytes == retained_lease_bytes
        && lease_after_drop.refusals == 0
        && lease_after_drop.first_refusal.is_none()
        && lease_after_drop.release_invariant_violations == 0;
    let scale_pass = step_pass && release_pass;
    emit_phase(
        &mut evidence,
        "sparse-d3q19-million-active",
        "drop_two_grids_and_shadow_lease_ns",
        reclaim,
        identity.root(),
    );

    let rss_after = linux_peak_rss_bytes();
    emit_process_peak_rss(
        &mut evidence,
        profile,
        "sparse-d3q19-million-active",
        identity.root(),
        rss_before,
        rss_after,
    );
    evidence.emit(
        if scale_pass {
            Severity::Info
        } else {
            Severity::Error
        },
        EventKind::Custom {
            name: "production-scale-sparse-receipt".to_string(),
            json: format!(
                "{{\"profile\":{},\"active_tiles\":{},\"active_cells\":{active_cells},\
                 \"published_population_values\":{published_values},\
                 \"retained_state_bytes_each\":{retained_state_bytes},\
                 \"serial_initial_mass\":{},\"serial_final_mass\":{},\
                 \"pooled_initial_mass\":{},\"pooled_final_mass\":{},\
                 \"serial_relative_mass_residual\":{},\
                 \"pooled_relative_mass_residual\":{},\"mass_roundoff_bound\":{},\
                 \"state_bits_equal\":{state_bits_equal},\"first_divergence\":{},\
                 \"populations_finite\":{populations_finite},\
                 \"placement_identity_round_trips\":{placement_round_trips},\
                 \"gate_requested\":{},\"reports\":[{},{}],\"pool\":{},\
                 \"shadow_lease_while_live\":{},\"shadow_lease_after_drop\":{},\
                 \"logical_scale_pass\":{scale_pass}}}",
                json_string(profile.name),
                SPARSE_ACTIVE_TILES,
                json_f64(serial_initial_mass),
                json_f64(serial_final_mass),
                json_f64(pooled_initial_mass),
                json_f64(pooled_final_mass),
                json_f64(serial_mass_residual),
                json_f64(pooled_mass_residual),
                json_f64(mass_bound),
                first_divergence.map_or_else(|| "null".to_string(), |index| index.to_string()),
                gate.is_requested(),
                reports
                    .first()
                    .map_or_else(|| "null".to_string(), RunReport::to_json),
                reports
                    .get(1)
                    .map_or_else(|| "null".to_string(), RunReport::to_json),
                pool_stats.to_json(),
                lease_while_live.to_json(),
                lease_after_drop.to_json()
            ),
        },
        None,
    );
    evidence.decision(
        "sparse-lbm-million-active-cells",
        &sparse_domain,
        if scale_pass {
            CapabilityDecision::Admitted
        } else {
            CapabilityDecision::Refused
        },
        format!(
            "exactly {SPARSE_ACTIVE_TILES} whole tiles / {active_cells} cells; serial/pool exact_state={state_bits_equal}; mass residuals={serial_mass_residual:e}/{pooled_mass_residual:e} <= {mass_bound:e}"
        ),
    );
    evidence.decision(
        "tile-pool-ownership",
        &sparse_domain,
        if scale_pass && report_pass {
            CapabilityDecision::Admitted
        } else {
            CapabilityDecision::Refused
        },
        format!(
            "one collide plus one stream pass retained {} reports over {expected_groups} canonical kernel groups with placement identity {placement_identity}",
            reports.len()
        ),
    );
    named_sparse_open_claims(&mut evidence, profile);
    evidence.verdict(
        "scale-sparse-lbm-million-active",
        scale_pass,
        format!(
            "active_tiles={}/{SPARSE_ACTIVE_TILES}; active_cells={active_cells}; state_bytes_each={retained_state_bytes}; exact_state={state_bits_equal}; first_divergence={first_divergence:?}; finite={populations_finite}; serial_mass={serial_initial_mass:.17e}->{serial_final_mass:.17e} ({serial_mass_residual:e}); pooled_mass={pooled_initial_mass:.17e}->{pooled_final_mass:.17e} ({pooled_mass_residual:e}); bound={mass_bound:e}; reports={}; pool={}; shadow_lease_live={}; shadow_lease_after={}",
            tiles.len(),
            reports.len(),
            pool_stats.to_json(),
            lease_while_live.to_json(),
            lease_after_drop.to_json()
        ),
    );
}
