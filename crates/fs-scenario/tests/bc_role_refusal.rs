//! Admission regression: machine-graph relations are never boundary kinds.

use fs_qty::{Dims, QtyAny};
use fs_scenario::ir::{parse_ir, write_ir};
use fs_scenario::{
    BcKind, BcValue, BoundaryCondition, Environment, Physics, Scenario, ScenarioError,
};

fn scenario_with_one_bc() -> Scenario {
    let mut scenario = Scenario::new("machine-role-refusal", 17, Environment::earth_lab());
    scenario.base_bcs.push(BoundaryCondition {
        region: "boundary/main".to_string(),
        physics: Physics::IncompressibleFlow,
        kind: BcKind::Dirichlet,
        value: Some(BcValue::Uniform(QtyAny::new(
            1.0,
            Dims([1, 0, -1, 0, 0, 0]),
        ))),
        compatibility: None,
        frame: 0,
    });
    scenario
}

#[test]
fn machine_graph_roles_are_structured_bc_parse_refusals() {
    let canonical = write_ir(&scenario_with_one_bc());
    assert!(canonical.contains(" incompressible-flow dirichlet "));

    for role in ["joint", "terminal", "controller", "reset"] {
        let smuggled = canonical.replacen(
            " incompressible-flow dirichlet ",
            &format!(" incompressible-flow {role} "),
            1,
        );
        let refusal = parse_ir(&smuggled).expect_err("machine role must not decode as a BC kind");
        assert_eq!(
            refusal,
            ScenarioError::ReservedBoundaryRole { role },
            "reserved role must retain its exact structured identity"
        );
        let diagnostic = refusal.to_string();
        assert!(
            diagnostic.contains("machine-graph role")
                && diagnostic.contains("boundary-condition kind")
                && diagnostic.contains("fs-ir machine relation"),
            "diagnostic did not identify the correct ownership boundary: {diagnostic}"
        );
    }
}

#[test]
fn ordinary_unknown_bc_kind_remains_distinct_from_reserved_role_refusal() {
    let canonical = write_ir(&scenario_with_one_bc());
    let unknown = canonical.replacen(
        " incompressible-flow dirichlet ",
        " incompressible-flow mystery-kind ",
        1,
    );
    let refusal = parse_ir(&unknown).expect_err("unknown BC kind must refuse");
    let ScenarioError::Parse { what, .. } = refusal else {
        panic!("expected parse refusal");
    };
    assert!(what.contains("unknown bc kind \"mystery-kind\""));
    assert!(!what.contains("reserved machine-graph role"));
}
