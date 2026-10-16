//! E4.3-ii battery (bead wf-root-guzez.5.6.2): the Wright v1 record
//! admits; the double-ownership HOSTILE TWIN refuses; missing effect
//! refuses; the noncirculatory class rule fires; the digest is
//! order-canonical AND content-sensitive; trim wiring through the
//! admitted record has no startup transient.
//! Repro: cargo test -p fs-flyer --test effectowners_battery

use fs_airfoil::unsteady::{UNSTEADY_SECTION_V1, UnsteadySectionState};
use fs_flyer::effectowners::{
    AeroEffect, AeroEffectOwners, OwnerAssignment, OwnerClass, wright_owners_v1,
};

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-flyer-e43ii\",\"case\":\"{case}\",{payload}}}");
}

#[test]
fn wright_v1_admits_and_digest_is_canonical() {
    let rec = wright_owners_v1();
    let admitted = rec.admit().unwrap();
    let d1 = admitted.digest();
    // Order canonicalization: reverse the assignment list — same digest.
    let mut rev = rec.clone();
    rev.assignments.reverse();
    let d2 = rev.admit().unwrap().digest();
    assert_eq!(d1, d2, "assignment-list order must not move the digest");
    // Content sensitivity: a different far-wake owner moves it.
    let mut other = rec.clone();
    other.assignments[5].owner_id = "fs-vpm.hybrid-wake-v1";
    other.assignments[5].class = OwnerClass::WakeModel;
    let d3 = other.admit().unwrap().digest();
    assert_ne!(d1, d3, "owner change must move the digest");
    jlog("digest", &format!("\"d1\":\"{}\"", &d1[..16]));
}

#[test]
fn double_ownership_hostile_twin_refuses() {
    // The classic double-count: an indicial component ALSO claiming the
    // noncirculatory (apparent-mass) effect alongside the added-mass
    // owner. This is exactly the bug class the record exists to kill.
    let mut rec = wright_owners_v1();
    rec.assignments.push(OwnerAssignment {
        effect: AeroEffect::Noncirculatory,
        owner_id: "fs-airfoil.unsteady.wagner-jones-2pole-v1",
        class: OwnerClass::IndicialKernel,
    });
    let err = rec.admit().unwrap_err();
    assert_eq!(err.code, "effect-owner-duplicate");
    assert!(err.message.contains("Noncirculatory"));
    jlog("hostile-twin", "\"duplicate_refused\":true");
}

#[test]
fn missing_effect_refuses() {
    let mut rec = wright_owners_v1();
    rec.assignments.retain(|a| a.effect != AeroEffect::FarWake);
    let err = rec.admit().unwrap_err();
    assert_eq!(err.code, "effect-owner-missing");
    assert!(err.message.contains("FarWake"));
    // And the empty record names the FIRST canonical effect.
    let empty = AeroEffectOwners {
        assignments: vec![],
    };
    assert_eq!(empty.admit().unwrap_err().code, "effect-owner-missing");
    jlog("missing", "\"orphan_refused\":true");
}

#[test]
fn noncirculatory_class_rule_fires() {
    // Same id count, same uniqueness — only the CLASS is wrong. The
    // admission must catch it (id checks alone cannot).
    let mut rec = wright_owners_v1();
    rec.assignments[2] = OwnerAssignment {
        effect: AeroEffect::Noncirculatory,
        owner_id: "fs-airfoil.unsteady.wagner-jones-2pole-v1",
        class: OwnerClass::IndicialKernel,
    };
    let err = rec.admit().unwrap_err();
    assert_eq!(err.code, "noncirculatory-owner-not-added-mass");
    // Empty owner id refuses regardless of structure.
    let mut rec2 = wright_owners_v1();
    rec2.assignments[0].owner_id = "";
    assert_eq!(rec2.admit().unwrap_err().code, "owner-id-empty");
    jlog("class-rule", "\"added_mass_only_enforced\":true");
}

#[test]
fn trim_wiring_has_no_startup_transient() {
    // The admitted record's indicial/separation owners are the E4.3-i
    // constructors; a section initialized at trim under the admitted
    // record must hold the static answer on the first tick (plan 5.1.5
    // trim-state initialization).
    let admitted = wright_owners_v1().admit().unwrap();
    assert_eq!(admitted.assignments().len(), 6);
    let alpha_trim = 0.06;
    let spec = UNSTEADY_SECTION_V1;
    let mut st = UnsteadySectionState::trim(&spec, alpha_trim);
    let first = st.advance(&spec, 0.1, alpha_trim, 0.0).unwrap();
    let static_cl = spec.cl_alpha * alpha_trim;
    assert!(
        (first.cl - static_cl).abs() < 1e-14,
        "startup transient under the admitted record: {} vs {static_cl}",
        first.cl
    );
    jlog(
        "trim-wiring",
        &format!("\"first_cl\":{},\"static\":{static_cl}", first.cl),
    );
}

#[test]
fn owners_golden_digest() {
    let digest = wright_owners_v1().admit().unwrap().digest();
    jlog("golden", &format!("\"digest\":\"{digest}\""));
    assert_eq!(
        digest, "bf68e945bbfd910e398254b094f584f0d9f76a0d1708d13a49009e90d699a3e9",
        "AeroEffectOwners digest moved — an ownership change is a MODEL \
         change and must ride the golden-bump protocol"
    );
}
