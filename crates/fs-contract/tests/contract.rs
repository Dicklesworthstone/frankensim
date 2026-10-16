//! Battery for assume-guarantee contracts (addendum Proposal E). Covers
//! interval validation + inclusive-boundary containment, envelope-containment
//! composition, the weakest-member soundness invariant, color discipline
//! (nonlinear cannot be verified), missing/outside conditions, cycle
//! detection, the diamond (shared sub-contract, not a cycle), and determinism.

use fs_contract::{
    Contract, ContractError, ContractLibrary, Envelope, Interval, OperatingConditions, compose,
};
use fs_evidence::{Color, ColorRank};
use fs_iface::SpaceType;

fn verified() -> Color {
    Color::Verified { lo: -1.0, hi: 1.0 }
}
fn estimated() -> Color {
    Color::Estimated {
        estimator: "surrogate".into(),
        dispersion: 3.0,
    }
}

fn traction_envelope(lo: f64, hi: f64) -> Envelope {
    Envelope::new().with("traction", SpaceType::HDiv, Interval::new(lo, hi).unwrap())
}

fn contract(name: &str, linear: bool, cert: Color, env: Envelope, requires: &[&str]) -> Contract {
    Contract {
        name: name.into(),
        interface: SpaceType::HDiv,
        linear,
        envelope: env,
        guarantee: format!("{name} holds"),
        certificate: cert,
        requires: requires.iter().map(|s| (*s).to_string()).collect(),
    }
}

#[test]
fn intervals_validate_and_contain_inclusively() {
    assert!(Interval::new(1.0, 2.0).is_ok());
    assert!(matches!(
        Interval::new(2.0, 1.0),
        Err(ContractError::BadInterval { .. })
    ));
    assert!(matches!(
        Interval::new(f64::NAN, 1.0),
        Err(ContractError::BadInterval { .. })
    ));
    let outer = Interval::new(0.0, 100.0).unwrap();
    assert!(outer.contains(&Interval::new(10.0, 50.0).unwrap()));
    // boundary is inclusive.
    assert!(outer.contains(&Interval::new(0.0, 100.0).unwrap()));
    // just outside is not contained.
    assert!(!outer.contains(&Interval::new(0.0, 100.0001).unwrap()));
    assert!(!outer.contains(&Interval::new(-0.0001, 50.0).unwrap()));
}

#[test]
fn contracts_compose_when_conditions_land_inside_envelopes() {
    let mut lib = ContractLibrary::new();
    lib.insert(contract(
        "beam",
        true,
        verified(),
        traction_envelope(0.0, 100.0),
        &[],
    ));
    lib.insert(contract(
        "joint",
        true,
        verified(),
        traction_envelope(0.0, 100.0),
        &[],
    ));
    lib.insert(contract(
        "assembly",
        true,
        verified(),
        Envelope::new(),
        &["beam", "joint"],
    ));
    let ops = OperatingConditions::new().with("traction", Interval::new(10.0, 50.0).unwrap());
    let claim = compose(&lib, "assembly", &ops).unwrap();
    assert_eq!(claim.members, vec!["assembly", "beam", "joint"]);
    // all verified -> composed verified.
    assert_eq!(claim.certificate.rank(), ColorRank::Verified);
}

#[test]
fn composed_claim_is_never_tighter_than_the_weakest_member() {
    // one verified, one estimated -> the SYSTEM claim is only estimated.
    let mut lib = ContractLibrary::new();
    lib.insert(contract(
        "beam",
        true,
        verified(),
        traction_envelope(0.0, 100.0),
        &[],
    ));
    lib.insert(contract(
        "joint",
        true,
        estimated(),
        traction_envelope(0.0, 100.0),
        &[],
    ));
    lib.insert(contract(
        "assembly",
        true,
        verified(),
        Envelope::new(),
        &["beam", "joint"],
    ));
    let ops = OperatingConditions::new().with("traction", Interval::new(10.0, 50.0).unwrap());
    let claim = compose(&lib, "assembly", &ops).unwrap();
    assert_eq!(
        claim.certificate.rank(),
        ColorRank::Estimated,
        "the weakest member (estimated) caps the system claim"
    );
}

#[test]
fn conditions_outside_the_envelope_are_rejected() {
    let mut lib = ContractLibrary::new();
    lib.insert(contract(
        "beam",
        true,
        verified(),
        traction_envelope(0.0, 100.0),
        &[],
    ));
    let ops = OperatingConditions::new().with("traction", Interval::new(10.0, 150.0).unwrap());
    match compose(&lib, "beam", &ops) {
        Err(ContractError::OutsideEnvelope { contract, quantity }) => {
            assert_eq!(contract, "beam");
            assert_eq!(quantity, "traction");
        }
        other => panic!("expected OutsideEnvelope, got {other:?}"),
    }
}

#[test]
fn a_missing_operating_condition_is_rejected() {
    let env = traction_envelope(0.0, 100.0).with(
        "temperature",
        SpaceType::L2,
        Interval::new(280.0, 350.0).unwrap(),
    );
    let mut lib = ContractLibrary::new();
    lib.insert(contract("beam", true, verified(), env, &[]));
    // ops provides traction but NOT temperature.
    let ops = OperatingConditions::new().with("traction", Interval::new(10.0, 50.0).unwrap());
    assert!(matches!(
        compose(&lib, "beam", &ops),
        Err(ContractError::MissingCondition { quantity, .. }) if quantity == "temperature"
    ));
}

#[test]
fn a_nonlinear_contract_cannot_be_verified_color() {
    let mut lib = ContractLibrary::new();
    // nonlinear + verified -> illegal (contracts are not exempt from the type system).
    lib.insert(contract(
        "plastic-joint",
        false,
        verified(),
        traction_envelope(0.0, 100.0),
        &[],
    ));
    let ops = OperatingConditions::new().with("traction", Interval::new(10.0, 50.0).unwrap());
    assert!(matches!(
        compose(&lib, "plastic-joint", &ops),
        Err(ContractError::ColorDiscipline { contract }) if contract == "plastic-joint"
    ));
    // the same nonlinear contract with estimated color is fine.
    let mut lib2 = ContractLibrary::new();
    lib2.insert(contract(
        "plastic-joint",
        false,
        estimated(),
        traction_envelope(0.0, 100.0),
        &[],
    ));
    assert!(compose(&lib2, "plastic-joint", &ops).is_ok());
}

#[test]
fn circular_contract_dependencies_are_rejected() {
    let mut lib = ContractLibrary::new();
    lib.insert(contract("a", true, verified(), Envelope::new(), &["b"]));
    lib.insert(contract("b", true, verified(), Envelope::new(), &["a"]));
    assert!(matches!(
        compose(&lib, "a", &OperatingConditions::new()),
        Err(ContractError::CircularDependency { .. })
    ));
}

#[test]
fn a_shared_sub_contract_is_not_a_cycle() {
    // diamond: root -> {left, right} -> shared. `shared` resolves once.
    let mut lib = ContractLibrary::new();
    lib.insert(contract("shared", true, verified(), Envelope::new(), &[]));
    lib.insert(contract(
        "left",
        true,
        verified(),
        Envelope::new(),
        &["shared"],
    ));
    lib.insert(contract(
        "right",
        true,
        verified(),
        Envelope::new(),
        &["shared"],
    ));
    lib.insert(contract(
        "root",
        true,
        verified(),
        Envelope::new(),
        &["left", "right"],
    ));
    let claim = compose(&lib, "root", &OperatingConditions::new()).unwrap();
    assert_eq!(claim.members, vec!["left", "right", "root", "shared"]);
    // shared appears exactly once.
    assert_eq!(claim.members.iter().filter(|m| *m == "shared").count(), 1);
}

#[test]
fn an_unknown_contract_is_rejected() {
    let lib = ContractLibrary::new();
    assert!(matches!(
        compose(&lib, "ghost", &OperatingConditions::new()),
        Err(ContractError::UnknownContract { name }) if name == "ghost"
    ));
}

#[test]
fn composition_is_deterministic() {
    let mut lib = ContractLibrary::new();
    lib.insert(contract(
        "beam",
        true,
        verified(),
        traction_envelope(0.0, 100.0),
        &[],
    ));
    lib.insert(contract(
        "joint",
        true,
        estimated(),
        traction_envelope(0.0, 100.0),
        &[],
    ));
    lib.insert(contract(
        "assembly",
        true,
        verified(),
        Envelope::new(),
        &["beam", "joint"],
    ));
    let ops = OperatingConditions::new().with("traction", Interval::new(10.0, 50.0).unwrap());
    assert_eq!(
        compose(&lib, "assembly", &ops),
        compose(&lib, "assembly", &ops)
    );
}
