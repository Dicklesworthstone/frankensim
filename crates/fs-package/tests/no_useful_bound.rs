//! Package projection refuses to color or certify `NoUsefulBound`.

use fs_evidence::{
    BoundInterval, BoundOutcome, ClaimClass, NoUsefulBoundCause, UsefulnessCriterion,
};
use fs_package::NoUsefulBoundRecord;

#[test]
fn package_manifest_keeps_the_refusal_visibly_outside_certificate_claims() {
    let criterion = UsefulnessCriterion::try_new("thermal trip", "kelvin", 2.0).expect("criterion");
    let outcome = BoundOutcome::refuse(
        BoundInterval::try_new(80.0, 120.0).expect("interval"),
        criterion,
        NoUsefulBoundCause::HorizonTooLong,
        ClaimClass::LongHorizonMeanLoad,
    );
    let record = NoUsefulBoundRecord::try_new(
        "thermal/long-horizon-trajectory",
        outcome.no_useful_bound().expect("typed refusal").clone(),
    )
    .expect("package projection");

    assert!(!record.has_certificate_claim());
    let manifest = record.render_manifest();
    assert!(manifest.contains("package-record=no-useful-bound"));
    assert!(manifest.contains("cause=horizon-too-long"));
    assert!(manifest.contains("suggested-reformulation=long-horizon-mean-load"));
    assert!(manifest.contains("certificate-claim=none"));
    assert!(manifest.contains("engineering compliance and scientific color not established"));
}
