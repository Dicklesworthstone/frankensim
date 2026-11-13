//! G0/G3: transported coordinates are audited separately, with no assumed
//! chemistry, positivity or shared numerical scale.

#![allow(clippy::float_cmp)]

use fs_ports::StableId;
use fs_ports::component::{
    ContractViolation, PortLayout, PortRole, PortSpec, StreamBasis, StreamQuantity,
    audit_carried_step,
};
use fs_qty::chemistry::SpeciesId;
use fs_qty::{Dims, Mass};

fn port(name: &str, role: PortRole) -> PortSpec {
    PortSpec::stream(StableId::new(name).unwrap(), role, 2).unwrap()
}

fn basis(quantities: Vec<StreamQuantity>) -> StreamBasis {
    StreamBasis::new(
        StableId::new("gas/source-v1").unwrap(),
        StableId::new("enthalpy/reference-v1").unwrap(),
        quantities,
    )
    .unwrap()
}

fn species(name: &str) -> StreamQuantity {
    StreamQuantity::SpeciesMass(SpeciesId::new(name).unwrap())
}

#[test]
fn identified_streams_preserve_order_dimensions_and_wire_arithmetic() {
    let quantities = vec![
        species("O2"),
        species("N2"),
        StreamQuantity::Named {
            id: StableId::new("momentum/world-x").unwrap(),
            dimensions: Dims([1, 1, -1, 0, 0, 0]),
        },
    ];
    let declared = basis(quantities.clone());
    assert_eq!(declared.quantities(), quantities);
    assert_eq!(declared.quantities()[0].dimensions(), Mass::DIMS);
    for role in [PortRole::FlowIn, PortRole::EffortIn] {
        let plain = PortSpec::stream(StableId::new("p").unwrap(), role, 3).unwrap();
        let identified =
            PortSpec::identified_stream(StableId::new("p").unwrap(), role, declared.clone())
                .unwrap();
        assert_eq!(identified.stream_basis(), Some(&declared));
        assert_eq!(identified.width(), 5);
        assert_eq!(identified.input_width(), plain.input_width());
        assert_eq!(identified.output_width(), plain.output_width());
        let (input, output) = if role == PortRole::FlowIn {
            ([0.25, 30.0, 0.05, 0.20, -0.4], [1e5, 120.0, 0.2, 0.8, -1.6])
        } else {
            ([1e5, 120.0, 0.2, 0.8, -1.6], [0.25, 30.0, 0.05, 0.20, -0.4])
        };
        assert_eq!(
            identified.energy_in_j(2.0, &input, &output),
            plain.energy_in_j(2.0, &input, &output)
        );
        assert_eq!(identified.mass_in_kg(2.0, &input, &output), 0.5);
        assert_eq!(
            identified.carried_quantity_in(2.0, &input, &output, 2),
            -0.8
        );
    }
    let empty =
        PortSpec::identified_stream(StableId::new("p").unwrap(), PortRole::FlowIn, basis(vec![]))
            .unwrap();
    assert_eq!(empty.width(), 2);
    assert!(empty.stream_basis().is_some());
}

#[test]
fn ambiguous_stream_coordinates_are_refused_at_declaration() {
    assert!(matches!(
        StreamBasis::new(
            StableId::new("gas").unwrap(),
            StableId::new("h").unwrap(),
            vec![species("N2"), species("N2")]
        ),
        Err(ContractViolation::DuplicateStreamQuantity { first: 0, index: 1 })
    ));
    let id = StableId::new("tracer").unwrap();
    assert!(matches!(
        StreamBasis::new(
            StableId::new("gas").unwrap(),
            StableId::new("h").unwrap(),
            vec![
                StreamQuantity::Named {
                    id: id.clone(),
                    dimensions: Mass::DIMS
                },
                StreamQuantity::Named {
                    id,
                    dimensions: Dims::NONE
                },
            ]
        ),
        Err(ContractViolation::DuplicateStreamQuantity { first: 0, index: 1 })
    ));
    for role in [PortRole::SignalIn, PortRole::SignalOut] {
        assert!(matches!(
            PortSpec::identified_stream(StableId::new("p").unwrap(), role, basis(vec![])),
            Err(ContractViolation::PortRole { .. })
        ));
    }
    assert!(matches!(
        PortSpec::stream(StableId::new("p").unwrap(), PortRole::FlowIn, usize::MAX),
        Err(ContractViolation::StreamWidth { .. })
    ));
    let ports = [
        PortSpec::identified_stream(
            StableId::new("a").unwrap(),
            PortRole::EffortIn,
            basis(vec![species("N2"), species("O2")]),
        )
        .unwrap(),
        PortSpec::identified_stream(
            StableId::new("b").unwrap(),
            PortRole::EffortIn,
            basis(vec![species("O2"), species("N2")]),
        )
        .unwrap(),
    ];
    assert!(matches!(
        PortLayout::new(&ports),
        Err(ContractViolation::StreamCoordinate { carried: 0, .. })
    ));
}

#[test]
fn both_roles_and_signed_inventories_use_transport_into_the_owner() {
    for role in [PortRole::FlowIn, PortRole::EffortIn] {
        let ports = [port("stream", role)];
        let layout = PortLayout::new(&ports).unwrap();
        let (input, output) = if role == PortRole::FlowIn {
            ([0.0, 0.0, -0.25, 1e-20], [1e5, 0.0, 0.0, 0.0])
        } else {
            ([1e5, 0.0, 0.0, 0.0], [0.0, 0.0, -0.25, 1e-20])
        };
        let a = audit_carried_step(
            &ports,
            &layout,
            2.0,
            &input,
            &output,
            0,
            Some(0.2),
            Some(-0.3),
        )
        .unwrap()
        .unwrap();
        assert_eq!(a.port_in, -0.5);
        assert!(a.balance_defect.abs() <= a.tolerance);
        let b = audit_carried_step(
            &ports,
            &layout,
            2.0,
            &input,
            &output,
            1,
            Some(0.0),
            Some(2e-20),
        )
        .unwrap()
        .unwrap();
        assert_eq!(b.port_in, 2e-20);
        assert_eq!(b.balance_defect, 0.0);
        assert!(matches!(
            audit_carried_step(
                &ports,
                &layout,
                2.0,
                &input,
                &output,
                1,
                Some(0.0),
                Some(3e-20)
            ),
            Err(ContractViolation::CarriedBalance { carried: 1, .. })
        ));
    }
}

#[test]
fn missing_nonfinite_and_unrepresentable_measurements_are_not_closure() {
    let ports = [port("stream", PortRole::FlowIn)];
    let layout = PortLayout::new(&ports).unwrap();
    let values = [0.0; 4];
    assert_eq!(
        audit_carried_step(&ports, &layout, 1.0, &values, &values, 0, None, None),
        Ok(None)
    );
    for (old, new) in [(Some(0.0), None), (None, Some(0.0))] {
        assert!(matches!(
            audit_carried_step(&ports, &layout, 1.0, &values, &values, 0, old, new),
            Err(ContractViolation::CarriedAvailability { carried: 0, .. })
        ));
    }
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(matches!(
            audit_carried_step(
                &ports,
                &layout,
                1.0,
                &values,
                &values,
                1,
                Some(bad),
                Some(0.0)
            ),
            Err(ContractViolation::StoredCarried { carried: 1, .. })
        ));
    }
    assert!(matches!(
        audit_carried_step(&ports, &layout, 0.0, &values, &values, 0, None, None),
        Err(ContractViolation::Step { .. })
    ));
    assert!(matches!(
        audit_carried_step(&ports, &layout, 1.0, &values[..3], &values, 0, None, None),
        Err(ContractViolation::VectorLength { .. })
    ));
    assert!(matches!(
        audit_carried_step(
            &ports,
            &layout,
            f64::MAX,
            &[0.0, 0.0, 2.0, 0.0],
            &values,
            0,
            None,
            None
        ),
        Err(ContractViolation::NonFinite {
            quantity: "carried port transfer"
        })
    ));
}

#[test]
fn scaled_compensation_retains_small_transfers_and_closes_large_cancellation() {
    let ports = [
        port("a", PortRole::FlowIn),
        port("b", PortRole::FlowIn),
        port("c", PortRole::FlowIn),
    ];
    let layout = PortLayout::new(&ports).unwrap();
    let mut input = [0.0; 12];
    input[2] = 1e16;
    input[6] = 1.0;
    input[10] = -1e16;
    let audit = audit_carried_step(
        &ports,
        &layout,
        1.0,
        &input,
        &[0.0; 12],
        0,
        Some(0.0),
        Some(1.0),
    )
    .unwrap()
    .unwrap();
    assert_eq!(audit.port_in, 1.0);
    assert_eq!(audit.balance_defect, 0.0);
    for rate in [1e308, 1e-300] {
        input[2] = rate;
        input[6] = 0.0;
        input[10] = -rate;
        let audit = audit_carried_step(
            &ports,
            &layout,
            1.0,
            &input,
            &[0.0; 12],
            0,
            Some(0.0),
            Some(0.0),
        )
        .unwrap()
        .unwrap();
        assert_eq!(audit.port_in, 0.0);
        assert_eq!(audit.balance_defect, 0.0);
        assert!(audit.tolerance.is_finite());
    }
}
