//! Ignored production-scale battery scaffold for `frankensim-ei3b`.
//!
//! This first tranche proves the arena/operation-lease path at the exact
//! 128^3 and 256^3 scalar-field shapes. It deliberately does not claim the
//! sparse-LBM, NUMA, CCD, quiet-host, cross-ISA, or macOS peak-RSS acceptance
//! rows that still need their owning implementations and admitted host runs.
//!
//! Run explicitly in release mode with one exact host profile:
//! `FRANKENSIM_PRODUCTION_SCALE_PROFILE=m4-128 cargo test -p fs-flagship-e2e
//! --release --test production_scale -- --ignored --exact
//! production_scale_arena_budget --nocapture`

use core::mem::size_of;
use fs_alloc::{
    AllocError, ArenaConfig, ArenaPool, ArenaStats, HUGEPAGE_BYTES, HugepagePolicy, LeaseReceipt,
    OperationMemoryLease, PoolStats, Site,
};
use fs_obs::{CapabilityDecision, EventKind, Severity};
use std::time::{Duration, Instant};

const SESSION: &str = "fs-flagship-e2e/production-scale";
const PROFILE_ENV: &str = "FRANKENSIM_PRODUCTION_SCALE_PROFILE";
const FIELD_SITE: Site = Site::named("production-scale/scalar-field");
const REFUSAL_SITE: Site = Site::named("production-scale/refusal-probe");

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

impl Evidence {
    fn new() -> Self {
        Self {
            emitter: fs_obs::Emitter::new(SESSION, "scale-battery"),
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
    let mut mutated = false;
    let mut arena_stats = (u64::MAX, u64::MAX, usize::MAX);

    let result = pool.scope_leased(&lease, |arena| {
        match arena.alloc_slice_fill(REFUSAL_SITE, cells, 0.0f64) {
            Ok(field) => {
                field[0] = 1.0;
                mutated = true;
                let stats = arena.stats();
                arena_stats = (
                    stats.allocated_bytes,
                    stats.allocation_count,
                    stats.chunk_count,
                );
                Ok(())
            }
            Err(error) => {
                let stats = arena.stats();
                arena_stats = (
                    stats.allocated_bytes,
                    stats.allocation_count,
                    stats.chunk_count,
                );
                Err(error)
            }
        }
    });

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

fn emit_phase(evidence: &mut Evidence, metric: &str, duration: Duration, configuration: u64) {
    let ns = duration_ns(duration);
    evidence.emit(
        Severity::Info,
        EventKind::Custom {
            name: "production-scale-phase".to_string(),
            json: format!(
                "{{\"kernel\":\"production-scale-arena\",\"metric\":{},\"value_ns\":{ns},\
                 \"unit\":\"ns\",\"configuration_identity_root\":{configuration},\
                 \"report_only\":true}}",
                json_string(metric)
            ),
        },
        Some(ns),
    );
}

fn named_open_claims(evidence: &mut Evidence, profile: Profile) {
    let domain = format!("{} production-scale acceptance", profile.name);
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
            "the arena/lease tranche cannot close the production-scale battery while sparse LBM, RSS gating, NUMA, CCD, and retained dual-host evidence remain open",
        ),
    ] {
        evidence.decision(capability, &domain, CapabilityDecision::Refused, detail);
    }
}

#[test]
#[ignore = "production-scale lane: requires release mode and an exact admitted host profile"]
fn production_scale_arena_budget() {
    let mut evidence = Evidence::new();
    let profile = match configured_profile() {
        Ok(Some(profile)) => profile,
        Ok(None) => {
            evidence.decision(
                "production-scale-profile",
                "ignored scale battery",
                CapabilityDecision::Refused,
                format!("named skip: set {PROFILE_ENV}=m4-128 or threadripper-256 explicitly"),
            );
            return;
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
    let release = !cfg!(debug_assertions);
    if !host_matches || !release {
        let detail = format!(
            "requested profile {} requires release {}+{} with model containing {:?}; \
             detected {}+{} model {:?}, release={release}",
            profile.name,
            profile.os,
            profile.arch,
            profile.model_needle,
            std::env::consts::OS,
            std::env::consts::ARCH,
            model
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
        "OS, architecture, model family, and release mode match; quiet-host and topology authority remain unproven",
    );

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
                "{{\"profile\":{},\"edge\":{},\"cells\":{cells},\"payload_bytes\":{payload_bytes},\
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
        "allocate_initialize_first_touch_ns",
        run.allocation,
        identity.root(),
    );
    emit_phase(&mut evidence, "serial_sweep_ns", run.sweep, identity.root());
    emit_phase(
        &mut evidence,
        "arena_drop_reclaim_ns",
        run.reclaim,
        identity.root(),
    );

    match (rss_before, rss_after) {
        (Some(before), Some(after)) => {
            evidence.emit(
                Severity::Info,
                EventKind::Custom {
                    name: "production-scale-process-memory".to_string(),
                    json: format!(
                        "{{\"metric\":\"linux_process_lifetime_peak_rss_bytes\",\
                         \"value_bytes\":{after},\"unit\":\"bytes\",\
                         \"configuration_identity_root\":{},\
                         \"source\":\"/proc/self/status:VmHWM\",\
                         \"process_lifetime\":true,\"report_only\":true}}",
                        identity.root()
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
    named_open_claims(&mut evidence, profile);
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
