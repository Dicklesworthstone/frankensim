//! G0/G3: independently measured mass closure, including energy-neutral leaks.

#![allow(clippy::float_cmp)]

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cancelled, Cx, ExecMode, StreamKey};
use fs_ports::StableId;
use fs_ports::component::conformance::{Case, Failure, Observed, check_component};
use fs_ports::component::{
    Class, Component, ContractViolation, PortLayout, PortRole, PortSpec, Proposal, Receipt,
    audit_mass_step, audit_step,
};

fn stream(name: &str, role: PortRole) -> PortSpec {
    PortSpec::stream(StableId::new(name).expect("id"), role, 0).expect("stream")
}

#[test]
fn stream_roles_and_stateless_transport_use_mass_positive_into_the_owner() {
    for role in [PortRole::FlowIn, PortRole::EffortIn] {
        let ports = [stream("port", role)];
        let layout = PortLayout::new(&ports).expect("layout");
        let (input, output) = if role == PortRole::FlowIn {
            ([0.25, 0.0], [101_325.0, 0.0])
        } else {
            ([101_325.0, 0.0], [0.25, 0.0])
        };
        let audit = audit_mass_step(&ports, &layout, 2.0, &input, &output, Some(1.0), Some(1.5))
            .expect("mass closes")
            .expect("measured inventory");
        assert_eq!(audit.port_mass_in_kg, 0.5);
        assert_eq!(audit.stored_change_kg, 0.5);
        assert!(audit.balance_defect_kg.abs() <= audit.tolerance_kg);
    }
    let ports = [
        stream("left", PortRole::EffortIn),
        stream("right", PortRole::EffortIn),
    ];
    let layout = PortLayout::new(&ports).expect("layout");
    let input = [101_325.0, 0.0, 100_000.0, 0.0];
    for rate in [-0.4, 0.0, 0.4] {
        let audit = audit_mass_step(
            &ports,
            &layout,
            0.25,
            &input,
            &[rate, 0.0, -rate, 0.0],
            Some(0.0),
            Some(0.0),
        )
        .expect("transport closes")
        .expect("zero inventory is available");
        assert_eq!(audit.port_mass_in_kg, 0.0);
        assert_eq!(audit.balance_defect_kg, 0.0);
    }
}

struct Inventory {
    ports: [PortSpec; 1],
    leak_kg: f64,
}

impl Component for Inventory {
    type State = f64;
    type Error = Cancelled;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }
    fn class(&self) -> Class {
        Class::ExactLinear
    }
    fn stored_energy_j(&self, _: &f64) -> Result<f64, Cancelled> {
        Ok(0.0)
    }
    fn stored_mass_kg(&self, state: &f64) -> Result<Option<f64>, Cancelled> {
        Ok(Some(*state))
    }

    fn propose(
        &self,
        old: &f64,
        dt_s: f64,
        inputs: &[f64],
        new: &mut f64,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, Cancelled> {
        cx.checkpoint()?;
        *new = old + dt_s * inputs[0] + self.leak_kg;
        outputs.copy_from_slice(&[101_325.0, 0.0]);
        Ok(Proposal::Step(Receipt::lossless(0.0, 0.0)))
    }
}

#[test]
fn conformance_refuses_an_energy_balanced_mass_leak_and_preserves_cancellation() {
    let mut component = Inventory {
        ports: [stream("inlet", PortRole::FlowIn)],
        leak_kg: 0.0,
    };
    let old = 2.0;
    let case = Case {
        state: &old,
        stale: &7.0,
        dt_s: 0.5,
        inputs: &[0.25, 0.0],
    };
    let gate = CancelGate::new_clock_free();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let key = StreamKey {
            seed: 1,
            kernel_id: 2,
            tile: 0,
            iteration: 0,
        };
        let cx = Cx::new(&gate, arena, key, Budget::INFINITE, ExecMode::Deterministic);
        assert!(matches!(
            check_component(&component, &case, &cx),
            Ok(Observed::Stepped { .. })
        ));
        component.leak_kg = 0.05;
        let mut new = old;
        let mut output = [0.0; 2];
        let Proposal::Step(receipt) = component
            .propose(&old, case.dt_s, case.inputs, &mut new, &mut output, &cx)
            .expect("proposal")
        else {
            panic!("unexpected shortening")
        };
        let layout = PortLayout::new(component.ports()).expect("layout");
        assert!(
            audit_step(
                component.ports(),
                &layout,
                case.dt_s,
                case.inputs,
                &output,
                &receipt,
                0.0,
                0.0,
            )
            .is_ok(),
            "the leak is invisible to the energy audit"
        );
        assert!(matches!(
            check_component(&component, &case, &cx),
            Err(Failure::Contract(ContractViolation::MassBalance { .. }))
        ));
        assert_eq!(old, 2.0);
        gate.request();
        assert!(matches!(
            check_component(&component, &case, &cx),
            Err(Failure::Owner(Cancelled))
        ));
    });
}

#[test]
fn compensated_scaled_sums_retain_small_terms_and_balance_extreme_transfers() {
    let ports = [
        stream("first", PortRole::FlowIn),
        stream("small", PortRole::FlowIn),
        stream("last", PortRole::FlowIn),
    ];
    let layout = PortLayout::new(&ports).expect("layout");
    let output = [101_325.0, 0.0, 101_325.0, 0.0, 101_325.0, 0.0];
    let audit = audit_mass_step(
        &ports,
        &layout,
        1.0,
        &[1.0e16, 0.0, 1.0, 0.0, -1.0e16, 0.0],
        &output,
        Some(0.0),
        Some(1.0),
    )
    .expect("small transfer retained")
    .expect("inventory");
    assert_eq!(audit.port_mass_in_kg, 1.0);
    assert_eq!(audit.balance_defect_kg, 0.0);
    for rate in [1.0e308, 1.0e-300] {
        let audit = audit_mass_step(
            &ports,
            &layout,
            1.0,
            &[rate, 0.0, 0.0, 0.0, -rate, 0.0],
            &output,
            Some(0.0),
            Some(0.0),
        )
        .expect("balanced extreme transfers")
        .expect("inventory");
        assert_eq!(audit.port_mass_in_kg, 0.0);
        assert_eq!(audit.balance_defect_kg, 0.0);
        assert!(audit.tolerance_kg.is_finite());
    }
    let audit = audit_mass_step(
        &ports,
        &layout,
        1.0,
        &[-1.0e308, 0.0, 0.0, 0.0, 0.0, 0.0],
        &output,
        Some(1.0e308),
        Some(0.0),
    )
    .expect("finite inventory emptied")
    .expect("inventory");
    assert_eq!(audit.stored_change_kg, -1.0e308);
    assert_eq!(audit.balance_defect_kg, 0.0);
    let audit = audit_mass_step(
        &ports,
        &layout,
        1.0,
        &[0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
        &output,
        Some(1.0e16),
        Some(1.0e16),
    )
    .expect("roundoff-scale defect accepted")
    .expect("inventory");
    assert_eq!(
        audit.balance_defect_kg, 1.0,
        "accepted defects remain visible"
    );
}

#[test]
fn unavailable_inventories_still_require_valid_wires_and_finite_transfers() {
    let ports = [stream("port", PortRole::FlowIn)];
    let layout = PortLayout::new(&ports).expect("layout");
    let output = [101_325.0, 0.0];
    assert_eq!(
        audit_mass_step(&ports, &layout, 1.0, &[2.0, 0.0], &output, None, None),
        Ok(None)
    );
    assert!(matches!(
        audit_mass_step(&ports, &layout, 0.0, &[2.0, 0.0], &output, None, None),
        Err(ContractViolation::Step { .. })
    ));
    assert!(matches!(
        audit_mass_step(&ports, &layout, 1.0, &[2.0], &output, None, None),
        Err(ContractViolation::VectorLength {
            vector: "inputs",
            ..
        })
    ));
    assert!(matches!(
        audit_mass_step(&[], &layout, 1.0, &[2.0, 0.0], &output, None, None),
        Err(ContractViolation::VectorLength {
            vector: "ports",
            ..
        })
    ));
    assert!(matches!(
        audit_mass_step(
            &ports,
            &layout,
            1.0,
            &[2.0, 0.0],
            &[f64::NAN, 0.0],
            None,
            None
        ),
        Err(ContractViolation::NonFinite { .. })
    ));
    assert!(matches!(
        audit_mass_step(&ports, &layout, f64::MAX, &[2.0, 0.0], &output, None, None),
        Err(ContractViolation::NonFinite {
            quantity: "port mass transfer"
        })
    ));
}

#[test]
fn missing_or_invalid_inventories_never_masquerade_as_zero_storage() {
    let ports = [stream("port", PortRole::FlowIn)];
    let layout = PortLayout::new(&ports).expect("layout");
    let input = [0.0, 0.0];
    let output = [101_325.0, 0.0];
    for (old, new) in [(Some(0.0), None), (None, Some(0.0))] {
        assert!(matches!(
            audit_mass_step(&ports, &layout, 1.0, &input, &output, old, new),
            Err(ContractViolation::MassAvailability { old_available, new_available })
                if old_available == old.is_some() && new_available == new.is_some()
        ));
    }
    for bad in [-1.0, f64::NAN, f64::INFINITY] {
        for (old, new, expected) in [(bad, 0.0, "accepted"), (0.0, bad, "candidate")] {
            assert!(matches!(
                audit_mass_step(&ports, &layout, 1.0, &input, &output, Some(old), Some(new)),
                Err(ContractViolation::StoredMass { endpoint, .. }) if endpoint == expected
            ));
        }
    }
    assert!(matches!(
        audit_mass_step(&ports, &layout, 1.0, &input, &output, Some(1.0), Some(1.01)),
        Err(ContractViolation::MassBalance { .. })
    ));
}
