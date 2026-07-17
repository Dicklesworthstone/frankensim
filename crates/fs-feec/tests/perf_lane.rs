//! The fs-feec HIGH-ORDER PERF LANE (bead cwjn): sum-factorized apply
//! throughput at p = 4 against the MEASURED machine peak (fs-roofline
//! axes), plus the apply-throughput-vs-p sweep. Run explicitly in
//! release:
//! `FRANKENSIM_BASELINE_STORE=<jsonl> FRANKENSIM_FIRMWARE_ID=<id>
//! FRANKENSIM_PROMOTION_AUTHORITY_POLICY=<tsv>
//! FRANKENSIM_RETAINED_SOURCE_RECEIPTS=<txt>
//! FRANKENSIM_ROOFLINE_LEDGER=<db> cargo test -p fs-feec
//! --release --test perf_lane -- --ignored --nocapture`
//!
//! GOLDEN CONSTRAINT: this lane only MEASURES the existing apply — the
//! 0xaaf1_076a_196c_6902 output golden is untouched by construction.

use fs_feec::highorder::hex::TensorSpace;
use fs_math::det;
use fs_roofline::authority::{ConfiguredPromotionAuthority, MAX_PROMOTION_AUTHORITY_POLICY_BYTES};
use fs_roofline::{
    AttestedAxisBaselinePolicy, AttestedBaselineStore, AxisAdmissionSnapshot, AxisBaselinePolicy,
    BaselineAxes, BaselineCandidate, BaselineIdentity, BaselineStore, ContentHash,
    EXTERNAL_PERF_GATE_LEDGER_ENV, ExternalPerfGateLane, MachineAxes, PromotionAttestation,
    days_since_epoch_now, promote_baseline, record_external_perf_gate_at_path,
};
use std::collections::BTreeSet;
use std::io::Read as _;

const MAX_RETAINED_RECEIPT_INPUT_BYTES: usize = fs_roofline::baseline::MAX_BASELINE_STORE_BYTES;
const OBS_SUITE: &str = "fs-feec/perf-lane";

fn emit_observation(identity: &str, name: &str, severity: fs_obs::Severity, json: String) {
    let mut emitter = fs_obs::Emitter::new(OBS_SUITE, identity);
    let event = emitter.emit(
        severity,
        fs_obs::EventKind::Custom {
            name: name.to_string(),
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

fn finite_json_2(value: f64) -> String {
    if value.is_finite() {
        format!("{value:.2}")
    } else {
        "null".to_string()
    }
}

fn finite_json_3(value: f64) -> String {
    if value.is_finite() {
        format!("{value:.3}")
    } else {
        "null".to_string()
    }
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

fn fail_invalid_environment(identity: &str, reason: &str, attainment: Option<f64>) -> ! {
    let escaped_reason = json_escape(reason);
    match attainment {
        Some(value) => emit_observation(
            identity,
            "environment-invalid",
            fs_obs::Severity::Error,
            format!(
                "{{\"metric\":\"feec-gate\",\"verdict\":\"environment_invalid\",\
                 \"reason\":\"{escaped_reason}\",\"attainment\":{},\
                 \"machine\":\"{}-{}\"}}",
                finite_json_3(value),
                std::env::consts::OS,
                std::env::consts::ARCH
            ),
        ),
        None => emit_observation(
            identity,
            "environment-invalid",
            fs_obs::Severity::Error,
            format!(
                "{{\"metric\":\"feec-gate\",\"verdict\":\"environment_invalid\",\
                 \"reason\":\"{escaped_reason}\",\"attainment\":null,\
                 \"machine\":\"{}-{}\"}}",
                std::env::consts::OS,
                std::env::consts::ARCH
            ),
        ),
    }
    panic!("FEEC roofline evidence rejected: {reason}");
}

fn fail_invalid_numerics(identity: &str, reason: &str) -> ! {
    emit_observation(
        identity,
        "numerical-invalid",
        fs_obs::Severity::Error,
        format!(
            "{{\"metric\":\"feec-gate\",\"verdict\":\"numerical_invalid\",\
             \"reason\":\"{}\",\"machine\":\"{}-{}\"}}",
            json_escape(reason),
            std::env::consts::OS,
            std::env::consts::ARCH
        ),
    );
    panic!("FEEC numerical measurement rejected: {reason}");
}

/// FLOPs per element per apply for degree r (p = r + 1): 9 axis
/// contractions of 2·p⁴ each, plus 3·p³ accumulate adds.
fn flops_per_element(r: usize) -> f64 {
    let p = (r + 1) as f64;
    18.0 * det::powi(p, 4) + 3.0 * det::powi(p, 3)
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
        fingerprint: 0xFEE0_C001,
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
        Ok("retained-feec-evidence.ledger".to_string()),
    );
    assert_eq!(with_ledger, GateAdmission::Citable);
    assert_eq!(path.as_deref(), Some("retained-feec-evidence.ledger"));

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

fn measure_apply(m: usize, r: usize, reps: usize) -> (f64, f64) {
    let space = TensorSpace::new(m, r);
    let n = space.ndof();
    let u: Vec<f64> = (0..n).map(|i| (i as f64 * 0.37).sin()).collect();
    // Warm.
    let mut sink = space.apply_stiffness(&u)[0];
    // Best of 3 trials: the attainment claim is about machine
    // capability, so scheduler/thermal noise must not deflate it.
    let mut best = f64::INFINITY;
    for _ in 0..3 {
        let t0 = std::time::Instant::now();
        for _ in 0..reps {
            sink += space.apply_stiffness(&u)[n / 2];
        }
        best = best.min(t0.elapsed().as_secs_f64());
    }
    let elements = (m * m * m * reps) as f64;
    let gflops = elements * flops_per_element(r) / best / 1e9;
    (gflops, sink)
}

#[test]
#[ignore = "perf lane: run explicitly in release with --ignored"]
fn sum_factorized_attainment() {
    let axes = MachineAxes::probe();
    emit_observation(
        "axes/pre/measurement",
        "axes-pre",
        fs_obs::Severity::Info,
        format!("{{\"metric\":\"axes-pre\",\"axes\":{}}}", axes.to_jsonl()),
    );
    // Environment validity (bead 1n61): implausible axes mean the probe
    // itself was contaminated (contended/throttled machine), so BOTH the
    // numerator and denominator of attainment are garbage — refuse to
    // gate rather than emit a vacuous pass or a false failure.
    if let Some(reason) = axes.plausibility_error() {
        fail_invalid_environment("terminal/pre-axes/environment-invalid", reason, None);
    }
    let admission = prepare_admission(&axes);
    // The p-sweep table (r = 1..6), retained as measurement JSON.
    for r in 1..=6usize {
        let m = (48 / (r + 1)).max(6);
        let (gflops, sink) = measure_apply(m, r, 3);
        if !gflops.is_finite() || gflops < 0.0 {
            fail_invalid_environment(
                &format!("terminal/sweep/r-{r}/environment-invalid"),
                "non-finite or negative FEEC throughput",
                None,
            );
        }
        if !sink.is_finite() {
            fail_invalid_numerics(
                &format!("terminal/sweep/r-{r}/numerical-invalid"),
                "non-finite FEEC apply sink",
            );
        }
        emit_observation(
            &format!("sweep/r-{r}/measurement"),
            "feec-apply",
            fs_obs::Severity::Info,
            format!(
                "{{\"metric\":\"feec-apply\",\"r\":{r},\"m\":{m},\"gflops\":{},\
                 \"attainment_single\":{},\"sink\":{}}}",
                finite_json_2(gflops),
                finite_json_3(gflops / axes.peak_single_gflops),
                finite_json_3(sink)
            ),
        );
    }
    // THE GATE at p = 4 (r = 3, per the bead's p-convention: degree-4
    // tensor basis = 4 points per axis): >= 30% of measured
    // single-thread peak on THIS machine. Bead cwjn requires separately
    // admitted rows on both reference ISAs before the cross-ISA claim lands.
    let (gflops, sink) = measure_apply(12, 3, 6);
    let attainment = gflops / axes.peak_single_gflops;
    if !attainment.is_finite() || attainment < 0.0 {
        fail_invalid_environment(
            "terminal/gate-measurement/environment-invalid",
            "non-finite or negative FEEC attainment",
            None,
        );
    }
    if !sink.is_finite() {
        fail_invalid_numerics(
            "terminal/gate-measurement/numerical-invalid",
            "non-finite FEEC gate sink",
        );
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
    let (snapshot, configuration_refusal) = admission.snapshot(&axes, &post_axes);
    println!("{}", snapshot.receipt_json());
    let admission_json = snapshot.receipt_json();
    // Over-roof poisoning (bead 1n61): a kernel "beating" the probed
    // peak by >1.5x means the PEAK probe was contaminated, not that the
    // kernel is fast — the whole measurement is invalid, both directions.
    if attainment > 1.5 {
        fail_invalid_environment(
            "terminal/over-roof/environment-invalid",
            "attainment exceeds the credible roofline band",
            Some(attainment),
        );
    }
    let (gate_admission, ledger_path) =
        configured_citable_ledger(classify_gate_admission(&snapshot, configuration_refusal));
    let (citation_eligible, refusal) = match gate_admission {
        GateAdmission::Citable => (true, None),
        GateAdmission::ReportOnly(reason) => (false, Some(reason)),
        GateAdmission::EnvironmentInvalid(reason) => {
            fail_invalid_environment(
                "terminal/admission/environment-invalid",
                &reason,
                Some(attainment),
            );
        }
    };
    let target_met = citation_eligible && attainment >= 0.30;
    let refusal_json = refusal.as_deref().map_or_else(
        || "null".to_string(),
        |reason| format!("\"{}\"", json_escape(reason)),
    );
    let gate_json = format!(
        "{{\"metric\":\"feec-gate\",\"r\":3,\"gflops\":{gflops:.2},\
         \"attainment\":{attainment:.3},\"target\":0.30,\"target_met\":{target_met},\
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
            ExternalPerfGateLane::Feec,
            &snapshot,
            &gate_json,
        )
        .unwrap_or_else(|error| {
            fail_invalid_environment(
                "terminal/ledger/environment-invalid",
                &format!("cannot ledger authority-admitted FEEC gate: {error}"),
                Some(attainment),
            )
        });
        println!("{}", receipt.receipt_json());
    }
    println!("{gate_json}");
    if citation_eligible {
        assert!(
            target_met,
            "the authority-admitted p=4 sum-factorized apply clears 30% of measured peak: {attainment:.3}"
        );
    }
}
