//! G0/G3 fixtures for feature-gated I01.3 port equation lowering.

use core::num::NonZeroUsize;

use fs_couple::{
    AccountingBoundary, AmountFlowRate, BoundaryTreatment, ChemicalEnergyAccounting,
    ChemicalEnergyInput, ConservationRole, ConservativeJunction, CoordinateBinding,
    DeviatoricStressWork, DissipationEvidence, DissipationLaw, DissipativeRelation,
    EntropyFlowRate, FieldMeasureSide, MovingStreamEnthalpyChart, PortKind, PortOrientation,
    PortSchema, PortTimestamp, PortValueShape, PowerPairing, SourceClass, SourceOrReservoir,
    SpecificEnergy, StableId, StorageElement, StoragePotential, StreamChartBinding,
    StreamConstituentFlow, StreamConstituentId, StreamEnergyChart, StreamKinematics, StreamPort,
    StreamStressWorkConvention,
};
use fs_iface::SpaceType;
use fs_opdsl::{
    AccountingTermKind, CompiledInterfaceEquation, ConservativeJunctionSide, InterfaceEquationSpec,
    LossOwnershipId, OwnershipDisposition, PortDiscretization, PortEquationError,
    PortEquationSense, PortEquationSpec, PortPrimitiveKind, ReversibleSkewCoupling,
    ReversibleSkewSide, SpatialSupport, StreamEnergyChartKind, StreamEnergyOwnership,
    StreamEquationSpec, SystemExpr, compile_interface_equations, compile_port_equation,
    compile_port_equations, compile_stream_equation,
};
use fs_qty::chemistry::SpeciesId;
use fs_qty::{Dims, Force, MassFlowRate, Power, Velocity};

const POWER: Dims = Dims([2, 1, -3, 0, 0, 0]);

fn stable(value: &str) -> StableId {
    StableId::new(value).expect("fixture stable id")
}

fn coordinates(orientation: PortOrientation) -> CoordinateBinding {
    CoordinateBinding::new(stable("basis-main"), stable("frame-lab"), orientation)
}

fn scalar_schema(port: &str, orientation: PortOrientation) -> PortSchema {
    PortKind::MechanicalForceVelocity
        .scalar_seed_schema(
            stable(port),
            coordinates(orientation),
            PortTimestamp::new(stable("clock-main"), 17),
        )
        .expect("admitted scalar port")
}

fn spec(
    port: &str,
    sense: PortEquationSense,
    term_kind: AccountingTermKind,
    ownership: OwnershipDisposition,
) -> PortEquationSpec {
    PortEquationSpec::new(
        scalar_schema(port, PortOrientation::OutwardFromOwner),
        PortDiscretization::lumped(),
        sense,
        term_kind,
        ownership,
    )
}

fn admitted_stream(port: &str) -> StreamPort {
    let species = StreamConstituentId::Species(SpeciesId::new("H2").expect("species id"));
    let binding = StreamChartBinding::try_new(
        stable(port),
        stable("state-stream-v1"),
        stable("basis-stream-species-v1"),
        [species.clone()],
        stable("reference-stream-chemical-v1"),
        coordinates(PortOrientation::OutwardFromOwner),
        PortTimestamp::new(stable("clock-main"), 17),
        stable("gravity-stream-datum-v1"),
        StreamStressWorkConvention::CauchyTensionPositiveOutwardPower,
    )
    .expect("stream binding");
    let chemical_energy = ChemicalEnergyAccounting::try_new(
        &binding,
        ChemicalEnergyInput::IncludedInStatePotential {
            reference_state: binding.chemical_reference_state().clone(),
        },
    )
    .expect("exclusive chemical energy");
    let kinematics = StreamKinematics::try_new([Velocity::new(0.0); 3], SpecificEnergy::new(0.0))
        .expect("finite stream kinematics");
    let deviatoric_work = DeviatoricStressWork::try_new(
        &binding,
        Power::new(0.0),
        stable("operator-stream-deviatoric-v1"),
        stable("evidence-stream-deviatoric-v1"),
    )
    .expect("stream stress work");
    let chart = StreamEnergyChart::MovingStreamEnthalpy(Box::new(
        MovingStreamEnthalpyChart::try_new(
            binding.clone(),
            SpecificEnergy::new(10.0),
            kinematics,
            deviatoric_work,
            chemical_energy,
        )
        .expect("moving enthalpy chart"),
    ));
    StreamPort::try_new(
        binding,
        MassFlowRate::new(2.0),
        [
            StreamConstituentFlow::try_new(species, AmountFlowRate::new(0.5))
                .expect("constituent flow"),
        ],
        [Force::new(1.0), Force::new(2.0), Force::new(3.0)],
        Power::new(20.0),
        EntropyFlowRate::new(0.25),
        chart,
    )
    .expect("admitted stream")
}

fn stream_spec(port: &str, owner: &str) -> StreamEquationSpec {
    StreamEquationSpec::new(
        admitted_stream(port),
        stable(&format!("receipt-{port}")),
        StreamEnergyOwnership::Owned(stable(owner)),
    )
}

fn storage_spec(gradient: &str) -> PortEquationSpec {
    let primitive = StorageElement::new(
        stable("storage-relation"),
        scalar_schema("storage-port", PortOrientation::OutwardFromOwner),
        StoragePotential::Hamiltonian,
        stable("storage-state-v1"),
        NonZeroUsize::new(3).expect("nonzero storage state"),
        stable(gradient),
    )
    .expect("admitted storage primitive");
    PortEquationSpec::from_storage(
        primitive,
        PortDiscretization::lumped(),
        PortEquationSense::AsDeclared,
    )
}

fn dissipation_spec(evidence: DissipationEvidence) -> PortEquationSpec {
    let primitive = DissipativeRelation::new(
        stable("dissipation-relation"),
        scalar_schema("dissipation-port", PortOrientation::OutwardFromOwner),
        DissipationLaw::Viscous,
        stable("operator-viscous-v1"),
        evidence,
    )
    .expect("admitted dissipation primitive");
    PortEquationSpec::from_dissipation(
        primitive,
        PortDiscretization::lumped(),
        PortEquationSense::AsDeclared,
    )
}

fn source_spec(treatment: BoundaryTreatment) -> PortEquationSpec {
    let schema = scalar_schema("source-port", PortOrientation::OutwardFromOwner);
    let boundary = AccountingBoundary::new(
        stable("source-boundary"),
        schema.coordinates().clone(),
        treatment,
    );
    let primitive = SourceOrReservoir::new(
        stable("source-relation"),
        schema,
        SourceClass::Reservoir,
        boundary,
    )
    .expect("admitted source primitive");
    PortEquationSpec::from_source_or_reservoir(
        primitive,
        PortDiscretization::lumped(),
        PortEquationSense::AsDeclared,
    )
}

fn conservative_junction_specs(port_a: &str, port_b: &str) -> [PortEquationSpec; 2] {
    let primitive = ConservativeJunction::new(
        stable("conservative-junction"),
        scalar_schema(port_a, PortOrientation::OutwardFromOwner),
        scalar_schema(port_b, PortOrientation::OutwardFromOwner),
    )
    .expect("admitted conservative junction");
    PortEquationSpec::from_conservative_junction(primitive, PortDiscretization::lumped())
}

fn reversible_skew_specs(
    port_a: &str,
    port_b: &str,
    forward_operator: &str,
    adjoint_operator: &str,
    evidence: &str,
) -> [PortEquationSpec; 2] {
    let coupling = ReversibleSkewCoupling::new(
        stable("reversible-skew-coupling"),
        scalar_schema(port_a, PortOrientation::OutwardFromOwner),
        scalar_schema(port_b, PortOrientation::OutwardFromOwner),
        stable(forward_operator),
        stable(adjoint_operator),
        stable(evidence),
    )
    .expect("reversible skew descriptor");
    PortEquationSpec::from_reversible_skew_coupling(coupling, PortDiscretization::lumped())
}

#[test]
fn g0_scalar_port_lowers_to_typed_power_equation_and_receipt() {
    let compiled = compile_port_equation(spec(
        "port-mechanical-a",
        PortEquationSense::AsDeclared,
        AccountingTermKind::Reversible,
        OwnershipDisposition::NotApplicable,
    ))
    .expect("neutral scalar schema lowers");

    assert_eq!(compiled.system().fields().len(), 3);
    assert_eq!(compiled.system().equations().len(), 1);
    assert!(matches!(
        &compiled.system().equations()[0].rhs,
        SystemExpr::PortPair {
            kind: PortKind::MechanicalForceVelocity,
            measure_dims,
            ..
        } if *measure_dims == Dims::NONE
    ));
    assert_eq!(compiled.receipt().port_id(), "port-mechanical-a");
    assert_eq!(compiled.receipt().product_dims(), POWER);
    assert_eq!(compiled.receipt().sign(), 1);
    assert_eq!(
        compiled.receipt().term_kind(),
        AccountingTermKind::Reversible
    );
    assert!(!compiled.system().extension().is_empty());

    let json = compiled.receipt().to_json();
    assert!(json.contains("\"schema\":\"fs-opdsl-port-equation-receipt-v1\""));
    assert!(json.contains("\"authority\":\"structural-generated\""));
    assert!(json.contains("\"product_dims\":[2,1,-3,0,0,0]"));
    assert!(json.contains("closed-window conservation remain external"));
}

#[test]
fn g3_orientation_reversal_is_an_explicit_identity_bearing_negative_sign() {
    let declared = compile_port_equation(spec(
        "port-orientation",
        PortEquationSense::AsDeclared,
        AccountingTermKind::Reversible,
        OwnershipDisposition::NotApplicable,
    ))
    .expect("declared orientation");
    let reversed = compile_port_equation(spec(
        "port-orientation",
        PortEquationSense::Reversed,
        AccountingTermKind::Reversible,
        OwnershipDisposition::NotApplicable,
    ))
    .expect("reversed orientation");

    assert_ne!(
        declared.receipt().system_identity(),
        reversed.receipt().system_identity(),
        "orientation sense is semantic, not display metadata"
    );
    assert_eq!(reversed.receipt().sign(), -1);
    assert!(matches!(
        &reversed.system().equations()[0].rhs,
        SystemExpr::Scale(value, inner)
            if value.to_bits() == (-1.0f64).to_bits()
                && matches!(inner.as_ref(), SystemExpr::PortPair { .. })
    ));
}

#[test]
fn g0_field_duality_retains_space_roles_measure_and_component_shape() {
    let area = Dims([2, 0, 0, 0, 0, 0]);
    let kind = PortKind::MechanicalForceVelocity;
    let pointwise_effort = kind
        .canonical_effort_dimensions()
        .checked_minus(area)
        .expect("traction dimensions");
    let shape =
        PortValueShape::field(3, SpaceType::HGrad, SpaceType::HDiv).expect("nonempty field shape");
    let schema = PortSchema::try_new(
        stable("port-field"),
        kind,
        pointwise_effort,
        kind.canonical_flow_dimensions(),
        shape,
        coordinates(PortOrientation::OutwardFromOwner),
        PowerPairing::FieldDuality {
            measure_dimensions: area,
            measure_side: FieldMeasureSide::Effort,
        },
        PortTimestamp::new(stable("clock-main"), 23),
        [ConservationRole::Energy],
    )
    .expect("field duality schema");
    let compiled = compile_port_equation(PortEquationSpec::new(
        schema.clone(),
        PortDiscretization::field(12, 18).expect("nonempty dofs"),
        PortEquationSense::AsDeclared,
        AccountingTermKind::Storage,
        OwnershipDisposition::Owned(stable("storage-owner")),
    ))
    .expect("field schema lowers");

    assert_eq!(compiled.system().fields()[0].space.degree, 0);
    assert_eq!(compiled.system().fields()[0].space.n, 12);
    assert_eq!(compiled.system().fields()[1].space.degree, 2);
    assert_eq!(compiled.system().fields()[1].space.n, 18);
    assert!(
        compiled
            .system()
            .fields()
            .iter()
            .all(|field| field.support == SpatialSupport::BoundaryTrace)
    );
    assert!(matches!(
        &compiled.system().equations()[0].rhs,
        SystemExpr::PortPair {
            measure_dims,
            ..
        } if *measure_dims == area
    ));

    assert!(matches!(
        compile_port_equation(PortEquationSpec::new(
            schema.clone(),
            PortDiscretization::lumped(),
            PortEquationSense::AsDeclared,
            AccountingTermKind::Storage,
            OwnershipDisposition::Owned(stable("owner-two")),
        )),
        Err(PortEquationError::DiscretizationMismatch {
            expected: "field",
            actual: "lumped",
        })
    ));
    assert!(matches!(
        compile_port_equation(PortEquationSpec::new(
            schema,
            PortDiscretization::field(10, 18).expect("nonempty dofs"),
            PortEquationSense::AsDeclared,
            AccountingTermKind::Storage,
            OwnershipDisposition::Owned(stable("owner-three")),
        )),
        Err(PortEquationError::FieldComponentMismatch {
            variable: "effort",
            dofs: 10,
            components: 3,
        })
    ));
}

#[test]
fn g0_ownership_is_role_checked_and_unique_across_a_batch() {
    let duplicate_port = spec(
        "port-duplicate",
        PortEquationSense::AsDeclared,
        AccountingTermKind::Reversible,
        OwnershipDisposition::NotApplicable,
    );
    assert_eq!(
        compile_port_equations(vec![duplicate_port.clone(), duplicate_port])
            .expect_err("duplicate source identity must refuse"),
        PortEquationError::DuplicatePortId {
            port: "port-duplicate".to_string(),
        }
    );

    assert!(matches!(
        compile_port_equation(spec(
            "port-bad-reversible",
            PortEquationSense::AsDeclared,
            AccountingTermKind::Reversible,
            OwnershipDisposition::Owned(stable("impossible-owner")),
        )),
        Err(PortEquationError::OwnershipMismatch {
            term_kind: AccountingTermKind::Reversible,
            ..
        })
    ));
    assert!(matches!(
        compile_port_equation(spec(
            "port-bad-storage",
            PortEquationSense::AsDeclared,
            AccountingTermKind::Storage,
            OwnershipDisposition::ExplicitlyUnowned {
                rationale: stable("missing-storage-model"),
            },
        )),
        Err(PortEquationError::OwnershipMismatch {
            term_kind: AccountingTermKind::Storage,
            ..
        })
    ));

    let owned_loss = compile_port_equation(spec(
        "port-owned-loss",
        PortEquationSense::AsDeclared,
        AccountingTermKind::Dissipation,
        OwnershipDisposition::Owned(stable("loss-owner")),
    ))
    .expect("concrete loss ownership");
    let reversed_loss = compile_port_equation(spec(
        "port-owned-loss",
        PortEquationSense::Reversed,
        AccountingTermKind::Dissipation,
        OwnershipDisposition::Owned(stable("loss-owner")),
    ))
    .expect("orientation does not change physical loss ownership");
    let loss_id = owned_loss
        .receipt()
        .loss_ownership_id()
        .expect("owned loss mints a nominal identity");
    assert_eq!(Some(loss_id), reversed_loss.receipt().loss_ownership_id());
    assert_eq!(loss_id.to_hex().len(), 64);
    assert_eq!(LossOwnershipId::parse_hex(&loss_id.to_hex()), Some(loss_id));
    assert!(owned_loss.receipt().to_json().contains(&loss_id.to_hex()));

    let duplicate_owner = stable("shared-loss-owner");
    let error = compile_port_equations(vec![
        spec(
            "port-source",
            PortEquationSense::AsDeclared,
            AccountingTermKind::Source,
            OwnershipDisposition::Owned(duplicate_owner.clone()),
        ),
        spec(
            "port-dissipation",
            PortEquationSense::AsDeclared,
            AccountingTermKind::Dissipation,
            OwnershipDisposition::Owned(duplicate_owner),
        ),
    ])
    .expect_err("one owner cannot own two generated terms");
    assert!(matches!(
        error,
        PortEquationError::DuplicateOwnership {
            ref owner,
            ref first_port,
            ref second_port,
        } if owner == "shared-loss-owner"
            && first_port == "port-dissipation"
            && second_port == "port-source"
    ));
}

#[test]
fn g3_batch_order_is_canonical_and_explicit_unowned_loss_is_retained() {
    let source = spec(
        "port-z-source",
        PortEquationSense::AsDeclared,
        AccountingTermKind::Source,
        OwnershipDisposition::Owned(stable("source-owner")),
    );
    let loss = spec(
        "port-a-loss",
        PortEquationSense::Reversed,
        AccountingTermKind::Dissipation,
        OwnershipDisposition::ExplicitlyUnowned {
            rationale: stable("outside-model-scope"),
        },
    );
    let forward = compile_port_equations(vec![source.clone(), loss.clone()])
        .expect("forward declaration order");
    let reverse = compile_port_equations(vec![loss, source]).expect("reverse declaration order");

    let forward_ids: Vec<_> = forward
        .equations()
        .iter()
        .map(|equation| equation.receipt().system_identity())
        .collect();
    let reverse_ids: Vec<_> = reverse
        .equations()
        .iter()
        .map(|equation| equation.receipt().system_identity())
        .collect();
    assert_eq!(forward_ids, reverse_ids);
    assert_eq!(forward.equations()[0].receipt().port_id(), "port-a-loss");
    assert!(
        forward.equations()[0]
            .receipt()
            .to_json()
            .contains("explicitly-unowned:outside-model-scope")
    );
    assert_eq!(forward.equations()[0].receipt().loss_ownership_id(), None);
}

#[test]
fn g0_empty_shape_and_metadata_resource_bombs_refuse_before_generation() {
    assert_eq!(
        PortDiscretization::field(0, 1),
        Err(PortEquationError::ZeroFieldDofs { variable: "effort" })
    );
    assert!(matches!(
        compile_port_equations(Vec::new()),
        Err(PortEquationError::EmptyBatch)
    ));

    let oversized_id = "a".repeat(5_000);
    let oversized = PortKind::MechanicalForceVelocity
        .scalar_seed_schema(
            StableId::new(oversized_id).expect("upstream stable ids are not length capped"),
            coordinates(PortOrientation::OutwardFromOwner),
            PortTimestamp::new(stable("clock-main"), 29),
        )
        .expect("upstream schema admits the identifier");
    assert!(matches!(
        compile_port_equation(PortEquationSpec::new(
            oversized,
            PortDiscretization::lumped(),
            PortEquationSense::AsDeclared,
            AccountingTermKind::Reversible,
            OwnershipDisposition::NotApplicable,
        )),
        Err(PortEquationError::CompilerMetadataTooLarge { .. })
    ));

    let oversized_stream_receipt = StableId::new("r".repeat(5_000))
        .expect("upstream stable receipt ids are not length capped");
    assert!(matches!(
        compile_stream_equation(StreamEquationSpec::new(
            admitted_stream("stream-resource-bomb"),
            oversized_stream_receipt,
            StreamEnergyOwnership::Owned(stable("stream-resource-owner")),
        )),
        Err(PortEquationError::CompilerMetadataTooLarge { .. })
    ));
    assert!(matches!(
        compile_interface_equations(Vec::new()),
        Err(PortEquationError::EmptyBatch)
    ));
}

#[test]
fn g0_stream_bundle_lowers_all_conserved_rate_blocks_with_bound_provenance() {
    let compiled = compile_stream_equation(stream_spec("stream-bundle-a", "stream-owner-a"))
        .expect("admitted stream lowers");

    assert_eq!(compiled.system().fields().len(), 5);
    assert!(compiled.system().equations().is_empty());
    let spaces: Vec<_> = compiled
        .system()
        .fields()
        .iter()
        .map(|field| (field.space.n, field.space.dims))
        .collect();
    assert!(spaces.contains(&(1, Dims([0, 1, -1, 0, 0, 0]))));
    assert!(spaces.contains(&(1, Dims([0, 0, -1, 0, 0, 1]))));
    assert!(spaces.contains(&(3, Dims([1, 1, -2, 0, 0, 0]))));
    assert!(spaces.contains(&(1, POWER)));
    assert!(spaces.contains(&(1, Dims([2, 1, -3, -1, 0, 0]))));
    assert!(
        compiled
            .system()
            .fields()
            .iter()
            .all(|field| field.support == SpatialSupport::BoundaryTrace)
    );

    let receipt = compiled.receipt();
    assert_eq!(receipt.port_id(), "stream-bundle-a");
    assert_eq!(receipt.source_receipt(), "receipt-stream-bundle-a");
    assert_eq!(receipt.constituent_count(), 1);
    assert_eq!(
        receipt.chart_kind(),
        StreamEnergyChartKind::MovingStreamEnthalpy
    );
    assert_eq!(receipt.energy_flow_bits(), 20.0_f64.to_bits());
    assert!(matches!(
        receipt.energy_ownership(),
        StreamEnergyOwnership::Owned(owner) if owner.as_str() == "stream-owner-a"
    ));
    let json = receipt.to_json();
    assert!(json.contains("\"schema\":\"fs-opdsl-stream-equation-receipt-v1\""));
    assert!(json.contains("\"authority\":\"structural-generated\""));
    assert!(json.contains("upstream chart evidence"));
}

#[test]
fn g3_mixed_batch_is_permutation_stable_and_refuses_stream_energy_double_count() {
    let duplicated_port = spec(
        "shared-boundary",
        PortEquationSense::AsDeclared,
        AccountingTermKind::Reversible,
        OwnershipDisposition::NotApplicable,
    );
    let duplicated_stream = stream_spec("shared-boundary", "shared-stream-owner");
    assert_eq!(
        compile_interface_equations(vec![
            InterfaceEquationSpec::Port(duplicated_port),
            InterfaceEquationSpec::Stream(duplicated_stream),
        ])
        .expect_err("one stable boundary cannot carry energy twice"),
        PortEquationError::DuplicateEnergyCarrier {
            port: "shared-boundary".to_string(),
        }
    );

    let port = spec(
        "port-a",
        PortEquationSense::AsDeclared,
        AccountingTermKind::Reversible,
        OwnershipDisposition::NotApplicable,
    );
    let stream = stream_spec("stream-z", "stream-owner-z");
    let forward = compile_interface_equations(vec![
        InterfaceEquationSpec::Port(port.clone()),
        InterfaceEquationSpec::Stream(stream.clone()),
    ])
    .expect("forward mixed batch");
    let reverse = compile_interface_equations(vec![
        InterfaceEquationSpec::Stream(stream),
        InterfaceEquationSpec::Port(port),
    ])
    .expect("reverse mixed batch");
    let forward_ids: Vec<_> = forward
        .equations()
        .iter()
        .map(CompiledInterfaceEquation::port_id)
        .collect();
    let reverse_ids: Vec<_> = reverse
        .equations()
        .iter()
        .map(CompiledInterfaceEquation::port_id)
        .collect();
    assert_eq!(forward_ids, vec!["port-a", "stream-z"]);
    assert_eq!(forward_ids, reverse_ids);
    assert!(matches!(
        forward.equations(),
        [
            CompiledInterfaceEquation::Port(_),
            CompiledInterfaceEquation::Stream(_)
        ]
    ));

    let shared_owner = stable("cross-kind-owner");
    let owned_port = spec(
        "port-owned",
        PortEquationSense::AsDeclared,
        AccountingTermKind::Source,
        OwnershipDisposition::Owned(shared_owner.clone()),
    );
    let owned_stream = StreamEquationSpec::new(
        admitted_stream("stream-owned"),
        stable("receipt-stream-owned"),
        StreamEnergyOwnership::Owned(shared_owner),
    );
    assert!(matches!(
        compile_interface_equations(vec![
            InterfaceEquationSpec::Stream(owned_stream),
            InterfaceEquationSpec::Port(owned_port),
        ]),
        Err(PortEquationError::DuplicateOwnership {
            ref owner,
            ref first_port,
            ref second_port,
        }) if owner == "cross-kind-owner"
            && first_port == "port-owned"
            && second_port == "stream-owned"
    ));
}

#[test]
fn g3_closed_primitive_metadata_is_identity_bearing_and_owner_derived() {
    let storage_a = compile_port_equation(storage_spec("gradient-storage-a"))
        .expect("storage primitive lowers");
    let storage_b = compile_port_equation(storage_spec("gradient-storage-b"))
        .expect("alternate storage operator lowers");
    assert_ne!(
        storage_a.receipt().system_identity(),
        storage_b.receipt().system_identity(),
        "constitutive-gradient identity is semantic"
    );
    assert_eq!(
        storage_a.receipt().primitive_kind(),
        Some(PortPrimitiveKind::Storage)
    );
    assert_eq!(storage_a.receipt().primitive_id(), Some("storage-relation"));
    assert!(matches!(
        storage_a.receipt().ownership(),
        OwnershipDisposition::Owned(owner) if owner.as_str() == "storage-relation"
    ));

    let monotone = compile_port_equation(dissipation_spec(DissipationEvidence::Monotonicity(
        stable("evidence-monotone-v1"),
    )))
    .expect("monotone dissipation lowers");
    let production = compile_port_equation(dissipation_spec(
        DissipationEvidence::NonnegativeProduction(stable("evidence-production-v1")),
    ))
    .expect("production-evidenced dissipation lowers");
    assert_ne!(
        monotone.receipt().system_identity(),
        production.receipt().system_identity(),
        "evidence kind and receipt are semantic"
    );
    assert_ne!(
        monotone.receipt().loss_ownership_id(),
        production.receipt().loss_ownership_id(),
        "nominal loss ownership binds its constitutive evidence"
    );
    assert_eq!(
        monotone.receipt().primitive_kind(),
        Some(PortPrimitiveKind::Dissipation)
    );

    let included = compile_port_equation(source_spec(BoundaryTreatment::IncludedSourceTerm))
        .expect("included source lowers");
    let reservoir =
        compile_port_equation(source_spec(BoundaryTreatment::ExternalReservoirExchange))
            .expect("external reservoir lowers");
    assert_ne!(
        included.receipt().system_identity(),
        reservoir.receipt().system_identity(),
        "accounting-boundary treatment is semantic"
    );
    assert_eq!(
        reservoir.receipt().primitive_kind(),
        Some(PortPrimitiveKind::SourceOrReservoir)
    );
    assert_eq!(reservoir.receipt().primitive_id(), Some("source-relation"));
    let json = reservoir.receipt().to_json();
    assert!(json.contains("\"primitive_kind\":\"source-or-reservoir\""));
    assert!(json.contains("\"primitive_id\":\"source-relation\""));
}

#[test]
fn g3_conservative_junction_lowers_signed_sides_and_preserves_permutation() {
    let ab = compile_port_equations(conservative_junction_specs("junction-a", "junction-b").into())
        .expect("A/B junction sides lower");
    let ab_a = ab
        .equations()
        .iter()
        .find(|equation| equation.receipt().port_id() == "junction-a")
        .expect("A-side equation");
    let ab_b = ab
        .equations()
        .iter()
        .find(|equation| equation.receipt().port_id() == "junction-b")
        .expect("B-side equation");

    assert_eq!(ab_a.receipt().sign(), 1);
    assert_eq!(ab_b.receipt().sign(), -1);
    assert_eq!(
        ab_a.receipt().primitive_kind(),
        Some(PortPrimitiveKind::ConservativeJunction)
    );
    assert_eq!(
        ab_a.receipt().conservative_junction_side(),
        Some(ConservativeJunctionSide::A)
    );
    assert_eq!(
        ab_b.receipt().conservative_junction_side(),
        Some(ConservativeJunctionSide::B)
    );
    assert_eq!(ab_a.receipt().primitive_id(), Some("conservative-junction"));
    assert!(matches!(
        ab_a.receipt().ownership(),
        OwnershipDisposition::NotApplicable
    ));
    assert!(
        ab_a.receipt()
            .to_json()
            .contains("\"conservative_junction_side\":\"a\"")
    );

    let ba = compile_port_equations(conservative_junction_specs("junction-b", "junction-a").into())
        .expect("B/A junction sides lower");
    let ba_a = ba
        .equations()
        .iter()
        .find(|equation| equation.receipt().port_id() == "junction-a")
        .expect("permuted A port");
    let ba_b = ba
        .equations()
        .iter()
        .find(|equation| equation.receipt().port_id() == "junction-b")
        .expect("permuted B port");

    assert_eq!(ba_a.receipt().sign(), -1);
    assert_eq!(ba_b.receipt().sign(), 1);
    assert_eq!(
        ba_a.receipt().conservative_junction_side(),
        Some(ConservativeJunctionSide::B)
    );
    assert_eq!(
        ba_b.receipt().conservative_junction_side(),
        Some(ConservativeJunctionSide::A)
    );
    assert_ne!(
        ab_a.receipt().system_identity(),
        ba_a.receipt().system_identity(),
        "permuting a junction changes the exact side role and sign"
    );
    assert_eq!(
        i16::from(ab_a.receipt().sign()) + i16::from(ab_b.receipt().sign()),
        0,
        "the two structural power contributions retain balancing signs"
    );
}

#[test]
fn g3_reversible_skew_coupling_retains_actions_evidence_and_signed_roles() {
    let ab = compile_port_equations(
        reversible_skew_specs(
            "skew-a",
            "skew-b",
            "operator-a-effort-to-b-flow-v1",
            "operator-b-effort-to-a-flow-adjoint-v1",
            "evidence-skew-adjoint-v1",
        )
        .into(),
    )
    .expect("A/B reversible skew sides lower");
    let ab_a = ab
        .equations()
        .iter()
        .find(|equation| equation.receipt().port_id() == "skew-a")
        .expect("A-side skew equation");
    let ab_b = ab
        .equations()
        .iter()
        .find(|equation| equation.receipt().port_id() == "skew-b")
        .expect("B-side skew equation");

    assert_eq!(ab_a.receipt().sign(), -1);
    assert_eq!(ab_b.receipt().sign(), 1);
    assert_eq!(
        ab_a.receipt().primitive_kind(),
        Some(PortPrimitiveKind::ReversibleSkewCoupling)
    );
    assert_eq!(
        ab_a.receipt().reversible_skew_side(),
        Some(ReversibleSkewSide::A)
    );
    assert_eq!(
        ab_b.receipt().reversible_skew_side(),
        Some(ReversibleSkewSide::B)
    );
    assert_eq!(
        ab_a.receipt().reversible_skew_forward_operator(),
        Some("operator-a-effort-to-b-flow-v1")
    );
    assert_eq!(
        ab_a.receipt().reversible_skew_adjoint_operator(),
        Some("operator-b-effort-to-a-flow-adjoint-v1")
    );
    assert_eq!(
        ab_a.receipt().reversible_skew_evidence(),
        Some("evidence-skew-adjoint-v1")
    );
    assert!(matches!(
        ab_a.receipt().ownership(),
        OwnershipDisposition::NotApplicable
    ));
    assert_eq!(ab_a.receipt().loss_ownership_id(), None);
    assert_eq!(ab_b.receipt().loss_ownership_id(), None);
    assert_eq!(
        i16::from(ab_a.receipt().sign()) + i16::from(ab_b.receipt().sign()),
        0,
        "the structural skew contributions retain opposite polarity"
    );

    let json = ab_a.receipt().to_json();
    assert!(json.contains("\"primitive_kind\":\"reversible-skew-coupling\""));
    assert!(json.contains("\"reversible_skew_side\":\"a-adjoint-negative\""));
    assert!(json.contains("\"reversible_skew_evidence\":\"evidence-skew-adjoint-v1\""));

    let changed_operator = compile_port_equations(
        reversible_skew_specs(
            "skew-a",
            "skew-b",
            "operator-a-effort-to-b-flow-v2",
            "operator-b-effort-to-a-flow-adjoint-v1",
            "evidence-skew-adjoint-v1",
        )
        .into(),
    )
    .expect("changed operator remains structurally admissible");
    let changed_a = changed_operator
        .equations()
        .iter()
        .find(|equation| equation.receipt().port_id() == "skew-a")
        .expect("changed A-side equation");
    assert_ne!(
        ab_a.receipt().system_identity(),
        changed_a.receipt().system_identity(),
        "forward action identity is semantic"
    );

    let changed_evidence = compile_port_equations(
        reversible_skew_specs(
            "skew-a",
            "skew-b",
            "operator-a-effort-to-b-flow-v1",
            "operator-b-effort-to-a-flow-adjoint-v1",
            "evidence-skew-adjoint-v2",
        )
        .into(),
    )
    .expect("changed evidence remains structurally admissible");
    let changed_evidence_a = changed_evidence
        .equations()
        .iter()
        .find(|equation| equation.receipt().port_id() == "skew-a")
        .expect("evidence-mutated A-side equation");
    assert_ne!(
        ab_a.receipt().system_identity(),
        changed_evidence_a.receipt().system_identity(),
        "skew-adjoint evidence identity is semantic"
    );

    let ba = compile_port_equations(
        reversible_skew_specs(
            "skew-b",
            "skew-a",
            "operator-a-effort-to-b-flow-v1",
            "operator-b-effort-to-a-flow-adjoint-v1",
            "evidence-skew-adjoint-v1",
        )
        .into(),
    )
    .expect("permuted reversible skew sides lower");
    let ba_a = ba
        .equations()
        .iter()
        .find(|equation| equation.receipt().port_id() == "skew-a")
        .expect("permuted skew-a equation");
    assert_eq!(ba_a.receipt().sign(), 1);
    assert_eq!(
        ba_a.receipt().reversible_skew_side(),
        Some(ReversibleSkewSide::B)
    );
    assert_ne!(
        ab_a.receipt().system_identity(),
        ba_a.receipt().system_identity(),
        "permuting the ordered cross action changes side role and identity"
    );
}

#[test]
fn g3_reversible_skew_coupling_refuses_aliased_relation_or_port_ids() {
    let port_a = scalar_schema("skew-aliased", PortOrientation::OutwardFromOwner);
    let port_b = scalar_schema("skew-distinct", PortOrientation::OutwardFromOwner);
    assert!(matches!(
        ReversibleSkewCoupling::new(
            stable("skew-aliased"),
            port_a,
            port_b,
            stable("operator-forward"),
            stable("operator-adjoint"),
            stable("evidence-skew"),
        ),
        Err(PortEquationError::PrimitivePortIdentityAlias {
            ref primitive,
            ref port,
        }) if primitive == "skew-aliased" && port == "skew-aliased"
    ));

    let duplicate_port = scalar_schema("skew-duplicate", PortOrientation::OutwardFromOwner);
    assert!(matches!(
        ReversibleSkewCoupling::new(
            stable("skew-distinct-relation"),
            duplicate_port.clone(),
            duplicate_port,
            stable("operator-forward"),
            stable("operator-adjoint"),
            stable("evidence-skew"),
        ),
        Err(PortEquationError::DuplicatePortId { ref port })
            if port == "skew-duplicate"
    ));
}
