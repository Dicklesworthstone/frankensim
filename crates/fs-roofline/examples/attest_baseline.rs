//! Operator ceremony helper: turn one `roofline promote` plain-store record
//! into the attested-run input set — attested baseline store, promotion
//! authority policy, and retained-receipt set — for the production CLI's
//! `--baseline/--authority-policy/--retained-receipts` flags.
//!
//! Usage:
//!   attest_baseline --store <plain.jsonl> --fingerprint <hex16>
//!       --key-id <id> --signature <sig>
//!       --out-attested <path> --out-authority <path> --out-receipts <path>
//!
//! The attestation is an operator statement, not cryptography: the emitted
//! authority policy authorizes exactly this key/signature over exactly this
//! baseline's promotion message, so any edit to the record, key, or signature
//! is a named refusal at admission (see fs-roofline::authority).

use fs_roofline::authority::{ConfiguredPromotionAuthority, PromotionAttestation};
use fs_roofline::{AttestedBaselineStore, BaselineStore};
use std::collections::BTreeMap;
use std::collections::BTreeSet;

fn main() -> Result<(), String> {
    let raw: Vec<String> = std::env::args().collect();
    let mut values = BTreeMap::new();
    let mut args = raw.iter().skip(1);
    while let Some(flag) = args.next() {
        if !matches!(
            flag.as_str(),
            "--store"
                | "--fingerprint"
                | "--key-id"
                | "--signature"
                | "--out-attested"
                | "--out-authority"
                | "--out-receipts"
        ) {
            return Err(format!("unknown attest_baseline argument {flag:?}"));
        }
        let value = args
            .next()
            .filter(|value| !value.is_empty())
            .ok_or_else(|| format!("{flag} requires a non-empty value"))?;
        if values.insert(flag.clone(), value.clone()).is_some() {
            return Err(format!("duplicate attest_baseline argument {flag:?}"));
        }
    }
    let get = |flag: &str| {
        values
            .get(flag)
            .cloned()
            .ok_or_else(|| format!("attest_baseline requires {flag}"))
    };
    let store_path = get("--store")?;
    let fingerprint = u64::from_str_radix(&get("--fingerprint")?, 16)
        .map_err(|error| format!("--fingerprint must be 16 hex digits: {error}"))?;
    let key_id = get("--key-id")?;
    let signature = get("--signature")?;

    let store_text = std::fs::read_to_string(&store_path)
        .map_err(|error| format!("cannot read plain store {store_path:?}: {error}"))?;
    let store =
        BaselineStore::from_jsonl(&store_text).map_err(|error| format!("plain store: {error}"))?;
    let baseline = store
        .for_fingerprint(fingerprint)
        .ok_or_else(|| format!("no baseline for machine fingerprint {fingerprint:016x}"))?
        .clone();

    // The authority policy authorizes exactly this key over exactly this
    // baseline's promotion message.
    let message_hex = fs_roofline::ContentHash(baseline.promotion_message()).to_hex();
    let authority_text = format!("authorize\t{key_id}\t{message_hex}\t{signature}\n");
    let authority = ConfiguredPromotionAuthority::from_text(&authority_text)
        .map_err(|error| format!("emitted authority policy is not canonical: {error}"))?;

    // Evidence must outlive the promotion it backs: retain every source
    // receipt named by the baseline's provenance.
    let receipts: BTreeSet<_> = baseline
        .provenance()
        .source_receipts()
        .iter()
        .copied()
        .collect();
    if receipts.is_empty() {
        return Err("baseline provenance names no source receipts".to_string());
    }
    let mut receipts_text = String::new();
    for receipt in &receipts {
        receipts_text.push_str(&receipt.to_hex());
        receipts_text.push('\n');
    }

    let mut attested = AttestedBaselineStore::new();
    attested
        .admit_verified(
            baseline,
            PromotionAttestation::new(key_id.clone(), signature),
            &authority,
            &receipts,
        )
        .map_err(|error| format!("attested admission refused: {error}"))?;

    let write = |flag: &str, contents: &str| -> Result<(), String> {
        let path = get(flag)?;
        std::fs::write(&path, contents).map_err(|error| format!("cannot write {path:?}: {error}"))
    };
    write("--out-attested", &attested.to_jsonl())?;
    write("--out-authority", &authority_text)?;
    write("--out-receipts", &receipts_text)?;
    println!(
        "{{\"attest\":\"ok\",\"fingerprint\":\"{fingerprint:016x}\",\"key_id\":\"{key_id}\",\"message\":\"{message_hex}\",\"source_receipts\":{}}}",
        receipts.len()
    );
    Ok(())
}
