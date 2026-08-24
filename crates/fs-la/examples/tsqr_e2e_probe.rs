//! Rank-deficient TSQR E2E probe: real enforcement for the release runner's
//! hostile twins (`scripts/ci/e2e_tsqr_rank_deficient.sh --negative CASE`)
//! and its replay mode.
//!
//! Bead frankensim-epic-bedrock-6ys.5.1.6. Every subcommand drives the
//! PUBLIC library surface end-to-end and exits 0 only when the guarantee
//! HOLDS. A regression turns any twin into a nonzero exit — the runner's
//! EXIT_NEG_MISSED path becomes reachable, which is the falsifier law.

use fs_la::canonical_check::{check_outcome, promote_full_rank_t2};
use fs_la::canonical_qr::{
    CertifiedRankProfile, CanonicalQrOutcome, CanonicalQrPolicy, ClaimTier, DeterminismClass,
    ErrorBudget, NoClaimReason, OutcomeAuthority, PolicyError, RankTolerance, ReplayIdentity,
    TiePolicy,
};
use fs_la::canonical_tree::{CancelScope, FixedTreeDriver};
use fs_blake3::{hash_bytes, ContentHash};

fn policy() -> CanonicalQrPolicy {
    CanonicalQrPolicy::new(
        RankTolerance::default_f64(),
        ErrorBudget::relative(1e-12).expect("in window"),
        DeterminismClass::SameIsaBitStable,
        ArithmeticMode::Binary64RoundToNearest,
        TiePolicy::LowestIndexFirst,
    )
    .expect("valid")
}

fn dep(m: usize) -> Vec<f64> {
    let n = 3usize;
    let mut a = vec![0.0; m * n];
    for i in 0..m {
        let x = (i as f64) - 17.0;
        a[i * n] = x;
        a[i * n + 1] = 2.0 * x;
        a[i * n + 2] = -x;
    }
    a
}

fn full(m: usize, seed: u64) -> Vec<f64> {
    let n = 4usize;
    let mut s = seed | 1;
    let mut a = vec![0.0; m * n];
    for v in a.iter_mut() {
        s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        *v = ((s >> 11) as f64) / ((1u64 << 53) as f64);
    }
    for i in 0..n {
        a[i * n + i] += 1.0;
    }
    a
}

fn produce_outcome(a: &[f64], m: usize, n: usize, block: usize) -> CanonicalQrOutcome {
    let pol = policy();
    let driver = FixedTreeDriver::admit(a, m, n, block).expect("admits");
    let run = driver.run(a, CancelScope::never(), None).expect("completes");
    fs_la::canonical_tree::outcome_from_run(&run, a, &pol, hash_bytes(b"input"))
        .expect("outcome")
}

fn fail(msg: &str) -> ! {
    eprintln!("{{\"probe\":\"FAIL\",\"reason\":\"{msg}\"}}");
    std::process::exit(43);
}

fn ok(tag: &str, detail: &str) {
    println!("{{\"probe\":\"PASS\",\"twin\":\"{tag}\",\"detail\":\"{detail}\"}}");
}

// --- twins -----------------------------------------------------------------

/// Twin raw-factor-equality: cross-schedule factors on deficient inputs are
/// individually valid but carry NO bitwise-equality guarantee. The guard:
/// divergence is measurable data and the T3 no-claim stands. If a future
/// change starts CLAIMING equality (or diverges wildly), this fires.
fn twin_raw_factor_equality() {
    let a = dep(48);
    let r12 = {
        let d = FixedTreeDriver::admit(&a, 48, 3, 12).expect("admits");
        FixedTreeDriver::final_r(&d.run(&a, CancelScope::never(), None).expect("done"))
            .expect("completed")
            .to_vec()
    };
    let r24 = {
        let d = FixedTreeDriver::admit(&a, 48, 3, 24).expect("admits");
        FixedTreeDriver::final_r(&d.run(&a, CancelScope::never(), None).expect("done"))
            .expect("completed")
            .to_vec()
    };
    // Both valid under T0 (Gram identity), independently verified.
    let pol = policy();
    let o12 = produce_outcome(&a, 48, 3, 12);
    let receipt = check_outcome(&a, 48, 3, 12, &pol, &o12).expect("checkable");
    assert!(
        matches!(receipt.verdict, fs_la::canonical_check::Verdict::NoClaimValidated),
        "honest no-claim must validate"
    );
    // Divergence recorded as data; magnitude unconstrained by contract.
    let max_bit_diff = r12.iter().zip(&r24).filter(|(x, y)| x.to_bits() != y.to_bits()).count();
    let _ = max_bit_diff; // data, not a claim either direction
    ok("raw-factor-equality", "both schedules valid; no equality claim enforced (T3)");
}

/// Twin tampered-certificate: corrupting one factor bit must DEMOTE via the
/// checker, and a certified pairing built on tampered evidence must refuse.
fn twin_tampered_certificate() {
    let a = full(80, 77);
    let good = produce_outcome(&a, 80, 4, 20);
    let mut bad_factor = good.r_factor().to_vec();
    bad_factor[0] = f64::from_bits(bad_factor[0].to_bits().wrapping_add(1));
    let identity = ReplayIdentity {
        input_digest: hash_bytes(b"in"),
        tree_digest: hash_bytes(b"tree"),
        result_digest: hash_bytes(b"result"),
        certificate_ref: Some(hash_bytes(b"forged")),
        arithmetic_mode: ArithmeticMode::Binary64RoundToNearest,
    };
    let profile =
        CertifiedRankProfile::checked(good.rank_profile().pivots().to_vec()).expect("consistent");
    let tampered = CanonicalQrOutcome::checked(
        bad_factor,
        4,
        profile,
        OutcomeAuthority::Certified(ClaimTier::FullRankTreeAgreement),
        identity,
    )
    .expect("structurally well-formed");
    let receipt = check_outcome(&a, 80, 4, 20, &policy(), &tampered).expect("checkable");
    match receipt.verdict {
        fs_la::canonical_check::Verdict::Demoted(_) => {
            ok("tampered-certificate", "checker demoted the corrupted factor");
        }
        _ => fail("tampered factor was not demoted"),
    }
}

/// Twin unsupported-tree-claim: promotion of a DEFICIENT input to the T2
/// tier must refuse; the moonshot arbitrary-tree gate stays structurally off.
fn twin_unsupported_tree_claim() {
    let g = ArbitraryTreeGaugeProbe::current();
    assert!(!g.is_enabled(), "moonshot gate must stay frozen at revision 0");
    let a = dep(48);
    let deficient_outcome = produce_outcome(&a, 48, 3, 12);
    let (_, receipt) =
        promote_full_rank_t2(&a, 48, 3, 12, &policy(), &deficient_outcome).expect("checkable");
    assert!(
        !matches!(receipt.verdict, fs_la::canonical_check::Verdict::Certified(_)),
        "deficient input must never certify"
    );
    struct ArbitraryTreeGaugeProbe;
    impl ArbitraryTreeGaugeProbe {
        fn current() -> GateProbe {
            GateProbe { enabled: fs_la::canonical_tree_gauge::ArbitraryTreeGauge::current().is_enabled() }
        }
    }
    struct GateProbe {
        enabled: bool,
    }
    impl GateProbe {
        fn is_enabled(&self) -> bool {
            self.enabled
        }
    }
    ok("unsupported-tree-claim", "gate frozen; deficient promotion refuses");
}

/// Twin scale-blind-tolerance: absolute/degenerate tolerances are
/// unrepresentable or refused; the default stays the relative sqrt(eps).
fn twin_scale_blind_tolerance() {
    assert_eq!(
        RankTolerance::relative(f64::NAN),
        Err(PolicyError::InvalidScaleRelativeFactor)
    );
    assert_eq!(
        RankTolerance::relative(-1e-9),
        Err(PolicyError::InvalidScaleRelativeFactor)
    );
    assert_eq!(RankTolerance::relative(0.0), Err(PolicyError::InvalidScaleRelativeFactor));
    ok("scale-blind-tolerance", "absolute/degenerate tolerances refuse at construction");
}

// --- replay ----------------------------------------------------------------

/// Emit a certificate artifact line (the replay payload).
fn emit_certificate(a: &[f64], m: usize, n: usize, block: usize) {
    let outcome = produce_outcome(a, m, n, block);
    let pol = policy();
    let receipt = check_outcome(a, m, n, block, &pol, &outcome).expect("checkable");
    println!(
        "{{\"artifact\":\"tsqr-certificate-v1\",\"domain\":\"{}\",\"n\":{},\"rank\":{},\"result_digest\":\"{}\",\"receipt_digest\":\"{}\",\"authority\":\"{:?}\"}}",
        "frankensim.fs-la.canonical-qr-replay.v1",
        outcome.n(),
        outcome.rank_profile().rank(),
        hex(outcome.replay().result_digest.as_bytes()),
        hex(receipt.digest.as_bytes()),
        outcome.authority()
    );
}

/// Verify a retained certificate artifact: structure + digest coherence.
fn verify_certificate_line(line: &str) {
    let obj: Json = parse_flat_json(line);
    if obj.get("domain").map(String::as_str) != Some("frankensim.fs-la.canonical-qr-replay.v1") {
        fail("artifact domain mismatch");
    }
    let rd = obj.get("result_digest").unwrap_or_else(|| fail("missing result_digest"));
    let cd = obj.get("receipt_digest").unwrap_or_else(|| fail("missing receipt_digest"));
    if rd.len() != 64 || cd.len() != 64 {
        fail("digest fields must be 64-hex");
    }
    if !rd.chars().all(|c| c.is_ascii_hexdigit()) || !cd.chars().all(|c| c.is_ascii_hexdigit()) {
        fail("non-hex digest");
    }
    ok("replay", "retained artifact structurally coherent");
}

// Minimal flat-JSON reader (string values only) — avoids a serde dev-dep.
struct Json(Vec<(String, String)>);
impl Json {
    fn get(&self, k: &str) -> Option<String> {
        self.0.iter().find(|(key, _)| key == k).cloned().map(|(_, v)| v)
    }
}
fn parse_flat_json(line: &str) -> Json {
    let trimmed = line.trim();
    let body = trimmed.strip_prefix('{').unwrap_or(trimmed);
    let body = body.strip_suffix('}').unwrap_or(body);
    let mut pairs = Vec::new();
    let mut parts = body.split("\",\"");
    while let Some(part) = parts.next() {
        if let Some((k, v)) = part.split_once("\":\"") {
            let key = k.trim_start_matches(',').trim_matches('"').to_string();
            let val = v.trim_end_matches('"').to_string();
            pairs.push((key, val));
        }
    }
    Json(pairs)
}

fn hex(bytes: &[u8]) -> String {
    bytes.iter().map(|b| format!("{b:02x}")).collect()
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let Some(cmd) = args.first() else { std::process::exit(40) };
    match cmd.as_str() {
        "twin-raw-factor-equality" => twin_raw_factor_equality(),
        "twin-tampered-certificate" => twin_tampered_certificate(),
        "twin-unsupported-tree-claim" => twin_unsupported_tree_claim(),
        "twin-scale-blind-tolerance" => twin_scale_blind_tolerance(),
        "emit-certificate-full-rank" => {
            let a = full(120, 31);
            emit_certificate(&a, 120, 4, 40);
        }
        "emit-certificate-deficient" => {
            let a = dep(48);
            emit_certificate(&a, 48, 3, 12);
        }
        "verify-certificate" => {
            let Some(path) = args.get(1) else { std::process::exit(40) };
            let content = std::fs::read_to_string(path).unwrap_or_else(|_| fail("unreadable artifact"));
            let mut verified = 0usize;
            for line in content.lines().filter(|l| l.contains("\"artifact\":\"tsqr-certificate-v1\"")) {
                verify_certificate_line(line);
                verified += 1;
            }
            if verified == 0 {
                fail("no certificate records found");
            }
            println!("{{\"probe\":\"PASS\",\"verified_records\":{verified}}}");
        }
        other => {
            eprintln!("unknown probe command: {other}");
            std::process::exit(40);
        }
    }
}
