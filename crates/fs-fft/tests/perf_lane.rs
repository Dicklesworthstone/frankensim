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
//! A citable gate requires an attested store plus configured authority:
//! `FRANKENSIM_BASELINE_STORE=<jsonl> FRANKENSIM_FIRMWARE_ID=<id>
//! FRANKENSIM_PROMOTION_AUTHORITY_POLICY=<tsv>
//! FRANKENSIM_RETAINED_SOURCE_RECEIPTS=<txt>
//! FRANKENSIM_ROOFLINE_LEDGER=<db>`. Plain, refused, or unledgered inputs
//! still measure, but emit only a frozen report-only snapshot.

use fs_fft::{C64, Fft, FftNd, SIXSTEP_FULL_ARRAY_PASSES, SIXSTEP_PERFORMANCE_MODEL_VERSION};
use fs_roofline::authority::{ConfiguredPromotionAuthority, MAX_PROMOTION_AUTHORITY_POLICY_BYTES};
use fs_roofline::{
    AttestedAxisBaselinePolicy, AttestedBaselineStore, AxisAdmissionSnapshot, AxisBaselinePolicy,
    BaselineAxes, BaselineCandidate, BaselineIdentity, BaselineStore, ContentHash,
    EXTERNAL_PERF_GATE_LEDGER_ENV, ExternalPerfGateLane, KernelSpec, MachineAxes,
    PromotionAttestation, RooflineKernel, TargetAxis, Threading, days_since_epoch_now, measure,
    promote_baseline, record_external_perf_gate_at_path,
};
use std::collections::BTreeSet;
use std::io::Read as _;

const MAX_RETAINED_RECEIPT_INPUT_BYTES: usize = fs_roofline::baseline::MAX_BASELINE_STORE_BYTES;

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
        "retained source receipts must be canonical newline-terminated lowercase hex".to_string()
    })?;
    if body.is_empty() {
        return Err("retained source receipts must contain at least one receipt".to_string());
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
) -> GateAdmission {
    if let Some(reason) = configuration_refusal {
        return GateAdmission::ReportOnly(reason);
    }
    match snapshot.baseline_citation_error() {
        Some(reason) => GateAdmission::EnvironmentInvalid(reason),
        None => GateAdmission::Citable,
    }
}

fn attach_citable_ledger(
    admission: GateAdmission,
    ledger_path: Result<String, String>,
) -> (GateAdmission, Option<String>) {
    if admission != GateAdmission::Citable {
        return (admission, None);
    }
    match ledger_path {
        Ok(path) if !path.is_empty() && path != ":memory:" => (GateAdmission::Citable, Some(path)),
        Ok(path) if path == ":memory:" => (
            GateAdmission::ReportOnly(format!(
                "{EXTERNAL_PERF_GATE_LEDGER_ENV} must name a durable ledger, not :memory:"
            )),
            None,
        ),
        Ok(_) => (
            GateAdmission::ReportOnly(format!(
                "{EXTERNAL_PERF_GATE_LEDGER_ENV} is empty; authority-admitted evidence was not ledgered"
            )),
            None,
        ),
        Err(reason) => (GateAdmission::ReportOnly(reason), None),
    }
}

fn configured_citable_ledger(admission: GateAdmission) -> (GateAdmission, Option<String>) {
    let ledger_path = std::env::var(EXTERNAL_PERF_GATE_LEDGER_ENV).map_err(|_| {
        format!(
            "{EXTERNAL_PERF_GATE_LEDGER_ENV} is missing; authority-admitted evidence was not ledgered"
        )
    });
    attach_citable_ledger(admission, ledger_path)
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

fn fail_invalid_environment(reason: &str) -> ! {
    let escaped_reason = json_escape(reason);
    println!(
        "{{\"metric\":\"fft-gate\",\"verdict\":\"environment_invalid\",\
         \"reason\":\"{escaped_reason}\",\"machine\":\"{}-{}\"}}",
        std::env::consts::OS,
        std::env::consts::ARCH
    );
    panic!("FFT roofline evidence rejected: {reason}");
}

#[allow(clippy::too_many_lines)] // Keep the fail-closed input/refusal order auditable in one place.
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

/// Use the production dispatch predicate directly so feature, shape, and
/// power-of-two admission cannot drift from the measured implementation.
fn takes_sixstep(n: usize) -> bool {
    Fft::takes_sixstep(n)
}

/// Full-array DRAM passes per single transform. The fused six-step does
/// exactly the implementation-declared two passes; the stage walk does one
/// per stage plus the odd-parity copy-back.
fn dram_passes(n: usize) -> f64 {
    if takes_sixstep(n) {
        SIXSTEP_FULL_ARRAY_PASSES as f64
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

#[test]
fn retained_receipt_config_is_strict_and_canonical() {
    let first = "00".repeat(32);
    let second = "ab".repeat(32);
    assert_eq!(
        parse_retained_receipts(&format!("{first}\n{second}\n"))
            .expect("canonical retained receipts")
            .len(),
        2
    );
    for malformed in [
        first.clone(),
        format!("{}\n", "AB".repeat(32)),
        format!("{first}\n{first}\n"),
        format!("{first}\n\n"),
        format!("{second}\n{first}\n"),
    ] {
        assert!(parse_retained_receipts(&malformed).is_err());
    }
}

struct ConfiguredAdmissionFixture {
    axes: MachineAxes,
    identity: BaselineIdentity,
    now_day: u64,
    store_text: String,
    authority_text: String,
    receipts_text: String,
    policy_receipt: ContentHash,
}

const MAX_STABLE_DAY_ATTEMPTS: usize = 3;

fn configured_admission_fixture_for_day(now_day: u64) -> ConfiguredAdmissionFixture {
    let axes = MachineAxes {
        fingerprint: 0xFF70_C001,
        cpu_brand: "synthetic-perf-host".to_string(),
        logical_cpus: 8,
        bandwidth_single_gbs: 40.0,
        bandwidth_all_core_gbs: 120.0,
        peak_single_gflops: 60.0,
        peak_all_core_gflops: 240.0,
    };
    let identity = BaselineIdentity::current(&axes, "fixture-firmware").expect("fixture identity");
    let candidates: Vec<_> = (1u8..=3)
        .map(|byte| {
            BaselineCandidate::from_receipt(axes.clone(), identity.clone(), ContentHash([byte; 32]))
                .expect("fixture candidate")
        })
        .collect();
    let baseline = promote_baseline(
        &candidates,
        "fixture-operator",
        "synthetic configured-admission proof",
        now_day,
        90,
    )
    .expect("fixture baseline");
    let key_id = "ops/perf-fixture";
    let signature = "fixture-signature";
    let authority_text = format!(
        "authorize\t{key_id}\t{}\t{signature}\n",
        ContentHash(baseline.promotion_message())
    );
    let authority =
        ConfiguredPromotionAuthority::from_text(&authority_text).expect("fixture authority");
    let retained: BTreeSet<_> = baseline
        .provenance()
        .source_receipts()
        .iter()
        .copied()
        .collect();
    let mut receipts_text = String::new();
    for receipt in &retained {
        receipts_text.push_str(&receipt.to_string());
        receipts_text.push('\n');
    }
    let mut store = AttestedBaselineStore::new();
    store
        .admit_verified(
            baseline,
            PromotionAttestation::new(key_id, signature),
            &authority,
            &retained,
        )
        .expect("fixture attested store");
    ConfiguredAdmissionFixture {
        axes,
        identity,
        now_day,
        store_text: store.to_jsonl(),
        authority_text,
        receipts_text,
        policy_receipt: authority.policy_receipt(),
    }
}

fn configured_admission_fixture() -> ConfiguredAdmissionFixture {
    let now_day = days_since_epoch_now().expect("fixture clock");
    configured_admission_fixture_for_day(now_day)
}

fn with_stable_configured_admission<T>(
    mut operation: impl FnMut(ConfiguredAdmissionFixture) -> T,
) -> T {
    for attempt in 1..=MAX_STABLE_DAY_ATTEMPTS {
        let before = days_since_epoch_now().expect("configured-attestation precondition clock");
        // Only the test baseline's promotion day is injected. `operation`
        // still mints its opaque production policy through the live-clock
        // `policy_for_run` boundary, and a rollover remints the whole fixture.
        let result = operation(configured_admission_fixture_for_day(before));
        let after = days_since_epoch_now().expect("configured-attestation postcondition clock");
        if before == after {
            return result;
        }
        assert!(
            attempt < MAX_STABLE_DAY_ATTEMPTS,
            "UTC epoch day changed from {before} to {after} during all \
             {MAX_STABLE_DAY_ATTEMPTS} configured-attestation attempts"
        );
    }
    unreachable!("bounded stable-day loop always returns or panics")
}

#[test]
fn configured_attestation_mints_exact_snapshot_and_refusals_stay_report_only() {
    let (fixture, snapshot, admission) = with_stable_configured_admission(|fixture| {
        let happy_store =
            AttestedBaselineStore::from_jsonl(&fixture.store_text).expect("fixture store parses");
        let (snapshot, refusal) = prepare_configured_attested_admission(
            &happy_store,
            fixture.identity.clone(),
            fixture.now_day,
            &fixture.authority_text,
            &fixture.receipts_text,
        )
        .snapshot(&fixture.axes, &fixture.axes);
        let admission = classify_gate_admission(&snapshot, refusal);
        (fixture, snapshot, admission)
    });
    let ConfiguredAdmissionFixture {
        axes,
        identity,
        now_day,
        store_text,
        authority_text,
        receipts_text,
        policy_receipt,
    } = fixture;
    assert_eq!(admission, GateAdmission::Citable);
    assert!(snapshot.authority_admitted());
    assert!(snapshot.verdict().trusted());
    assert!(snapshot.receipt_json().contains("\"tier\":\"attested\""));
    assert!(
        snapshot
            .receipt_json()
            .contains("\"verdict\":\"authorized\"")
    );
    assert!(
        snapshot
            .receipt_json()
            .contains(&format!("\"policy_receipt\":\"{policy_receipt}\""))
    );

    let assert_report_only = |candidate_store: &str,
                              candidate_identity: BaselineIdentity,
                              policy: &str,
                              expected: &str| {
        let store =
            AttestedBaselineStore::from_jsonl(candidate_store).expect("structural attested store");
        let (snapshot, refusal) = prepare_configured_attested_admission(
            &store,
            candidate_identity,
            now_day,
            policy,
            &receipts_text,
        )
        .snapshot(&axes, &axes);
        let refusal = refusal.expect("refused configuration stays report-only");
        assert!(refusal.contains(expected), "unexpected refusal: {refusal}");
        assert!(!snapshot.authority_admitted());
        assert!(!snapshot.baseline_citation_eligible());
        assert!(snapshot.receipt_json().contains("\"tier\":\"candidate\""));
    };
    assert_report_only(&store_text, identity.clone(), "", "unknown-key");
    assert_report_only(
        &store_text,
        identity.clone(),
        "revoke\tops/perf-fixture\n",
        "revoked-key",
    );
    let tampered = store_text.replace("fixture-signature", "tampered-signature");
    assert_report_only(&tampered, identity, &authority_text, "wrong-signature");
    let cross_machine =
        BaselineIdentity::current(&axes, "other-firmware").expect("cross-machine identity");
    assert_report_only(
        &store_text,
        cross_machine,
        &authority_text,
        "does not match",
    );
}

#[test]
fn config_refusal_is_report_only_but_attested_drift_is_environment_invalid() {
    let ConfiguredAdmissionFixture {
        axes,
        identity,
        now_day,
        store_text,
        receipts_text,
        ..
    } = configured_admission_fixture();
    let store = AttestedBaselineStore::from_jsonl(&store_text).expect("fixture store parses");
    let (candidate_snapshot, configuration_refusal) =
        prepare_configured_attested_admission(&store, identity, now_day, "", &receipts_text)
            .snapshot(&axes, &axes);
    assert!(matches!(
        classify_gate_admission(&candidate_snapshot, configuration_refusal),
        GateAdmission::ReportOnly(reason) if reason.contains("unknown-key")
    ));

    let (without_ledger, path) = attach_citable_ledger(
        GateAdmission::Citable,
        Err(format!("{EXTERNAL_PERF_GATE_LEDGER_ENV} is missing")),
    );
    assert!(matches!(
        without_ledger,
        GateAdmission::ReportOnly(reason) if reason.contains(EXTERNAL_PERF_GATE_LEDGER_ENV)
    ));
    assert!(path.is_none());
    let (ephemeral, path) =
        attach_citable_ledger(GateAdmission::Citable, Ok(":memory:".to_string()));
    assert!(matches!(
        ephemeral,
        GateAdmission::ReportOnly(reason) if reason.contains("durable ledger")
    ));
    assert!(path.is_none());
    let (with_ledger, path) = attach_citable_ledger(
        GateAdmission::Citable,
        Ok("retained-fft-evidence.ledger".to_string()),
    );
    assert_eq!(with_ledger, GateAdmission::Citable);
    assert_eq!(path.as_deref(), Some("retained-fft-evidence.ledger"));

    let (attested_snapshot, configuration_refusal) = with_stable_configured_admission(|fixture| {
        let store =
            AttestedBaselineStore::from_jsonl(&fixture.store_text).expect("fixture store parses");
        let mut degraded = fixture.axes.clone();
        degraded.bandwidth_single_gbs *= 0.5;
        prepare_configured_attested_admission(
            &store,
            fixture.identity,
            fixture.now_day,
            &fixture.authority_text,
            &fixture.receipts_text,
        )
        .snapshot(&fixture.axes, &degraded)
    });
    assert!(configuration_refusal.is_none());
    assert!(matches!(
        classify_gate_admission(&attested_snapshot, configuration_refusal),
        GateAdmission::EnvironmentInvalid(reason)
            if reason.contains("historical baseline admission refused")
    ));
}

#[test]
fn missing_retained_source_refuses_authority_and_preserves_candidate_identity() {
    let ConfiguredAdmissionFixture {
        axes,
        identity,
        now_day,
        store_text,
        authority_text,
        receipts_text,
        ..
    } = configured_admission_fixture();
    let store = AttestedBaselineStore::from_jsonl(&store_text).expect("fixture store parses");
    let baseline_hash = store
        .for_fingerprint(identity.fingerprint())
        .expect("fixture baseline exists")
        .content_hash();
    let mut incomplete = parse_retained_receipts(&receipts_text).expect("fixture receipts parse");
    let missing = *incomplete
        .iter()
        .next()
        .expect("fixture promotion has source receipts");
    assert!(incomplete.remove(&missing));
    let mut incomplete_text = String::new();
    for receipt in incomplete {
        incomplete_text.push_str(&receipt.to_string());
        incomplete_text.push('\n');
    }

    let (snapshot, refusal) = prepare_configured_attested_admission(
        &store,
        identity,
        now_day,
        &authority_text,
        &incomplete_text,
    )
    .snapshot(&axes, &axes);
    let refusal = refusal.expect("missing retained source must stay report-only");
    assert!(
        refusal.contains(&format!(
            "source receipt {missing} named by the promotion is not available"
        )),
        "unexpected refusal: {refusal}"
    );
    assert!(!snapshot.authority_admitted());
    assert!(!snapshot.baseline_citation_eligible());
    assert_eq!(snapshot.baseline_hash(), Some(baseline_hash));
    assert!(snapshot.receipt_json().contains("\"tier\":\"candidate\""));
}

#[test]
#[allow(clippy::too_many_lines)] // Keep the rotation proof and its stable-day retry atomic.
fn key_rotation_requires_exact_same_baseline_reendorsement() {
    let (
        revoked_snapshot,
        revoked_refusal,
        snapshot,
        admission,
        baseline_hash,
        policy_receipt,
        rotated_policy_receipt,
    ) = with_stable_configured_admission(|fixture| {
        let ConfiguredAdmissionFixture {
            axes,
            identity,
            now_day,
            store_text,
            receipts_text,
            policy_receipt,
            ..
        } = fixture;
        let mut store =
            AttestedBaselineStore::from_jsonl(&store_text).expect("fixture store parses");
        let baseline = store
            .for_fingerprint(identity.fingerprint())
            .expect("fixture baseline exists")
            .clone();
        let baseline_hash = baseline.content_hash();
        let retained = parse_retained_receipts(&receipts_text).expect("fixture receipts parse");
        let rotated_policy_text = format!(
            "authorize\tops/perf-fixture-rotated\t{}\tfixture-rotated-signature\n\
             revoke\tops/perf-fixture\n",
            ContentHash(baseline.promotion_message())
        );
        let rotated_authority = ConfiguredPromotionAuthority::from_text(&rotated_policy_text)
            .expect("canonical rotated authority policy");
        let rotated_policy_receipt = rotated_authority.policy_receipt();
        let (revoked_snapshot, revoked_refusal) = prepare_configured_attested_admission(
            &store,
            identity.clone(),
            now_day,
            &rotated_policy_text,
            &receipts_text,
        )
        .snapshot(&axes, &axes);
        store
            .admit_verified(
                baseline,
                PromotionAttestation::new("ops/perf-fixture-rotated", "fixture-rotated-signature"),
                &rotated_authority,
                &retained,
            )
            .expect("same immutable baseline is explicitly re-endorsed");
        assert_eq!(
            store
                .for_fingerprint(identity.fingerprint())
                .expect("re-endorsed baseline remains present")
                .content_hash(),
            baseline_hash
        );
        let (snapshot, refusal) = prepare_configured_attested_admission(
            &store,
            identity,
            now_day,
            &rotated_policy_text,
            &receipts_text,
        )
        .snapshot(&axes, &axes);
        let admission = classify_gate_admission(&snapshot, refusal);
        (
            revoked_snapshot,
            revoked_refusal,
            snapshot,
            admission,
            baseline_hash,
            policy_receipt,
            rotated_policy_receipt,
        )
    });
    assert_ne!(rotated_policy_receipt, policy_receipt);
    let revoked_refusal = revoked_refusal.expect("old key must refuse after rotation");
    assert!(
        revoked_refusal.contains("revoked-key"),
        "unexpected refusal: {revoked_refusal}"
    );
    assert!(!revoked_snapshot.authority_admitted());
    assert!(!revoked_snapshot.baseline_citation_eligible());
    assert_eq!(revoked_snapshot.baseline_hash(), Some(baseline_hash));
    assert_eq!(admission, GateAdmission::Citable);
    assert!(snapshot.authority_admitted());
    assert!(snapshot.verdict().trusted());
    assert_eq!(snapshot.baseline_hash(), Some(baseline_hash));
    assert!(
        snapshot
            .receipt_json()
            .contains("\"key_id\":\"ops/perf-fixture-rotated\"")
    );
    assert!(
        snapshot
            .receipt_json()
            .contains("\"signature\":\"fixture-rotated-signature\"")
    );
    assert!(
        snapshot
            .receipt_json()
            .contains(&format!("\"policy_receipt\":\"{rotated_policy_receipt}\""))
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
                SIXSTEP_PERFORMANCE_MODEL_VERSION
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
    fn run_once(&mut self) -> Result<(), String> {
        self.plan.forward(&mut self.data, &mut self.scratch);
        self.plan.inverse(&mut self.data, &mut self.scratch);
        Ok(())
    }
}

/// N-D pooled roundtrip (bead 27d3): the executor-tiled pencil path,
/// all axes parallel — measured against the ALL-CORE axes since the
/// TilePool owns placement. Generic over the pool lane so the parked
/// crew (bead tkr7) serves every axis pass and every row without
/// respawning — the per-run spawn/join overhead that made the first
/// N-D rows report-only is out of the measured path.
struct FftNdRoundTrip<'p, P> {
    dims: Vec<usize>,
    plan: FftNd,
    data: Vec<C64>,
    pool: &'p P,
    gate: fs_exec::CancelGate,
}

impl<'p, P: fs_exec::KernelRunner> FftNdRoundTrip<'p, P> {
    fn new(dims: &[usize], pool: &'p P) -> FftNdRoundTrip<'p, P> {
        let plan = FftNd::new(dims);
        let total = plan.total();
        FftNdRoundTrip {
            dims: dims.to_vec(),
            plan,
            data: (0..total)
                .map(|i| {
                    C64::new(
                        ((i * 37) % 101) as f64 * 0.02 - 1.0,
                        ((i * 53) % 97) as f64 * 0.02,
                    )
                })
                .collect(),
            pool,
            gate: fs_exec::CancelGate::new(),
        }
    }
}

impl<P: fs_exec::KernelRunner> RooflineKernel for FftNdRoundTrip<'_, P> {
    fn spec(&self) -> KernelSpec {
        // Per axis pass: gather one C64 + scatter one C64 per element
        // (32 B); a roundtrip runs every axis twice. Line/scratch
        // traffic is cache-resident and deliberately uncounted — the
        // model stays a lower bound on traffic, which keeps attainment
        // honest (never inflated).
        let axes_count = self.dims.len() as f64;
        let bf: f64 = self.dims.iter().map(|&n| butterfly_stages(n)).sum();
        KernelSpec {
            name: "fftnd-roundtrip",
            version: "27d3-nd1",
            bytes_per_elem: 2.0 * 32.0 * axes_count,
            flops_per_elem: 2.0 * 12.5 * bf + 2.0,
            threading: Threading::AllCore,
            target_axis: TargetAxis::BindingRoof,
            target_fraction: Some(0.40),
        }
    }
    fn elements(&self) -> usize {
        self.plan.total()
    }
    fn run_once(&mut self) -> Result<(), String> {
        self.plan
            .forward_pooled(&mut self.data, self.pool, &self.gate)
            .map_err(|error| format!("pooled forward failed: {error}"))?;
        self.plan
            .inverse_pooled(&mut self.data, self.pool, &self.gate)
            .map_err(|error| format!("pooled inverse failed: {error}"))?;
        Ok(())
    }
}

#[test]
fn fused_sixstep_traffic_and_evidence_version_are_bound() {
    assert_eq!(SIXSTEP_FULL_ARRAY_PASSES, 2);
    assert_eq!(SIXSTEP_PERFORMANCE_MODEL_VERSION, "27d3-6s-fused2");

    if cfg!(feature = "frontier-sixstep") {
        let n = 1usize << 16;
        assert!(takes_sixstep(n));
        assert_eq!(
            dram_passes(n).to_bits(),
            (SIXSTEP_FULL_ARRAY_PASSES as f64).to_bits()
        );
        let spec = FftRoundTrip::new(n).spec();
        assert_eq!(spec.version, SIXSTEP_PERFORMANCE_MODEL_VERSION);
        assert_eq!(spec.bytes_per_elem.to_bits(), 160.0f64.to_bits());
    }
}

/// N-D pooled rows (bead 27d3): measured on the PARKED-CREW lane (bead
/// tkr7) — one crew parked for the whole sweep, so per-axis-pass worker
/// spawn/join (the overhead that made the first rows report-only:
/// 0.011 attainment at 256×256 on a 5995WX) is out of the measured
/// path. Rows stay REPORT-ONLY until both baseline machines clear the
/// 0.40 floor with band margin — floors assert settled claims.
/// Returns false when any row is environment-invalid.
fn fftnd_report_rows(axes: &MachineAxes) -> bool {
    // Diagnostic override (bead 3f6c): sweep worker counts on the
    // report-only rows to locate the small-kernel granularity peak.
    // Timing-only by construction (P2) — results are bitwise-identical
    // at every worker count.
    let workers = std::env::var("FRANKENSIM_FFTND_WORKERS")
        .ok()
        .and_then(|raw| raw.parse::<usize>().ok())
        .filter(|&w| w >= 1)
        .unwrap_or_else(|| std::thread::available_parallelism().map_or(8, std::num::NonZero::get));
    let pool = fs_exec::TilePool::new(fs_exec::PoolConfig::for_host(workers, 0xFD1D));
    pool.with_parked_crew_local(|parked| {
        let mut env_ok = true;
        for dims in [vec![256usize, 256], vec![1024, 1024], vec![128, 128, 64]] {
            let mut kern = FftNdRoundTrip::new(&dims, parked);
            let att = measure(&mut kern, 1, 5, axes).expect("bounded FFT-ND measurement");
            let receipt = att.to_jsonl();
            println!(
                "{{\"metric\":\"fftnd-roundtrip\",\"dims\":{dims:?},\"workers\":{workers},\
                 \"lane\":\"parked-crew\",\"gated\":false,\"receipt\":{receipt}}}"
            );
            env_ok &= att.verdict != fs_roofline::Verdict::EnvironmentInvalid;
        }
        env_ok
    })
}

#[test]
#[ignore = "perf lane: run explicitly in release with --ignored"]
fn fft_attainment() {
    let axes = MachineAxes::probe();
    println!("{}", axes.to_jsonl());
    // Environment validity (bead 1n61): implausible axes poison both
    // the numerator and denominator of attainment — refuse up front.
    if let Some(reason) = axes.plausibility_error() {
        fail_invalid_environment(reason);
    }
    // Static floors and pre/post self-agreement cannot detect a host that is
    // consistently degraded through the whole run. Only the attested path can
    // authorize a positive gate; every configuration refusal remains measured.
    let admission = prepare_admission(&axes);
    // Size ladder: L2-resident (2^16) reported for context; the gate
    // rows are the memory-resident sizes (2^20, 2^22 — 32/128 MB
    // working sets against the DRAM STREAM axis).
    let mut gate_ok = true;
    let mut floor_ok = true;
    let mut env_ok = true;
    for &(n, gated) in &[(1usize << 16, false), (1 << 20, true), (1 << 22, true)] {
        let mut kern = FftRoundTrip::new(n);
        let att = measure(&mut kern, 1, 5, &axes).expect("bounded FFT measurement");
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
    env_ok &= fftnd_report_rows(&axes);
    if !env_ok {
        fail_invalid_environment("contaminated environment");
    }
    // Post-run reprobe + one frozen authority/baseline decision.
    let post_axes = MachineAxes::probe();
    println!(
        "{{\"metric\":\"axes-post\",\"axes\":{}}}",
        post_axes.to_jsonl()
    );
    let (snapshot, configuration_refusal) = admission.snapshot(&axes, &post_axes);
    println!("{}", snapshot.receipt_json());
    let admission_json = snapshot.receipt_json();
    // The 0.40 target is REPORTED per row; measured 0.26–0.43 across
    // runs on this machine, dominated by axis and load noise from
    // concurrent agent builds (bead 27d3 records the numbers). The
    // ASSERTED gate is the anti-collapse floor, and only the authority-admitted
    // path can assert it. Report-only rows retain the observed comparison.
    let (gate_admission, ledger_path) =
        configured_citable_ledger(classify_gate_admission(&snapshot, configuration_refusal));
    let (citation_eligible, refusal) = match gate_admission {
        GateAdmission::Citable => (true, None),
        GateAdmission::ReportOnly(reason) => (false, Some(reason)),
        GateAdmission::EnvironmentInvalid(reason) => fail_invalid_environment(&reason),
    };
    let target_met = citation_eligible && gate_ok;
    let floor_met = citation_eligible && floor_ok;
    let refusal_json = refusal.as_deref().map_or_else(
        || "null".to_string(),
        |reason| format!("\"{}\"", json_escape(reason)),
    );
    let gate_json = format!(
        "{{\"metric\":\"fft-gate\",\"target\":0.40,\"target_met\":{target_met},\
         \"floor\":0.15,\"floor_met\":{floor_met},\"observed_floor_met\":{floor_ok},\
         \"citation_eligible\":{citation_eligible},\"recorded\":{citation_eligible},\
         \"report_only\":{},\
         \"reason\":{refusal_json},\"admission\":{admission_json},\
         \"machine\":\"{}-{}\"}}",
        !citation_eligible,
        std::env::consts::OS,
        std::env::consts::ARCH
    );
    if let Some(path) = ledger_path {
        let receipt = record_external_perf_gate_at_path(
            &path,
            ExternalPerfGateLane::Fft,
            &snapshot,
            &gate_json,
        )
        .unwrap_or_else(|error| {
            fail_invalid_environment(&format!(
                "cannot ledger authority-admitted FFT gate: {error}"
            ))
        });
        println!("{}", receipt.receipt_json());
    }
    println!("{gate_json}");
    if citation_eligible {
        assert!(
            floor_ok,
            "authority-admitted memory-resident FFT round trips collapsed below the 15% floor"
        );
    }
}
