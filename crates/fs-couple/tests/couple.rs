//! Battery for port-Hamiltonian coupling (fs-couple). Covers power-conjugate
//! ports, the Dirac interconnection's exact power conservation, the energy
//! audit (passivity measured, not assumed), the Aitken relaxation factor, and
//! the load-bearing added-mass comparison: naive staggering diverges where
//! Aitken-relaxed coupling stays stable. The PortSchema v2 PR-1/PR-2 battery
//! pins scalar migration, four thermodynamic primitives, and extended kinds.

use core::num::NonZeroUsize;

use fs_couple::{
    AccountingBoundary, AitkenRelaxation, BoundaryTreatment, ConservationRole,
    ConservativeJunction, CoordinateBinding, CoupleError, DissipationEvidence, DissipationLaw,
    DissipativeRelation, EnergyAudit, FieldMeasureSide, FsiResult, PORT_SCHEMA_VERSION, Port,
    PortKind, PortOrientation, PortPrimitive, PortSchema, PortTimestamp, PortValueShape,
    PortVariable, PowerPairing, SourceClass, SourceOrReservoir, StableId, StorageElement,
    StoragePotential, interconnect, interface_power, iterate_aitken, iterate_fixed_relaxation,
};
use fs_iface::SpaceType;
use fs_qty::{Area, Dims, Force, Power, Pressure, Temperature, Velocity, VolumetricFlowRate};

fn stable(value: &str) -> StableId {
    StableId::new(value).unwrap_or_else(|error| panic!("invalid test ID {value:?}: {error:?}"))
}

fn scalar_schema(
    kind: PortKind,
    id: &str,
    frame: &str,
    orientation: PortOrientation,
    tick: u64,
) -> PortSchema {
    kind.scalar_seed_schema(
        stable(id),
        CoordinateBinding::new(stable("basis/si-scalar"), stable(frame), orientation),
        PortTimestamp::new(stable("clock/coupling-window"), tick),
    )
    .unwrap_or_else(|error| panic!("scalar seed migration failed for {id}: {error:?}"))
}

fn assert_fsi_result_bits_eq(actual: FsiResult, expected: FsiResult) {
    assert_eq!(actual.converged, expected.converged);
    assert_eq!(actual.steps, expected.steps);
    assert_eq!(actual.solution.to_bits(), expected.solution.to_bits());
    assert_eq!(
        actual.final_residual.to_bits(),
        expected.final_residual.to_bits()
    );
}

#[test]
fn port_schema_v2_scalar_seed_migration_goldens() {
    let cases = [
        (
            PortKind::MechanicalForceVelocity,
            Force::DIMS,
            Velocity::DIMS,
            &[ConservationRole::Energy][..],
        ),
        (
            PortKind::FluidPressureFlux,
            Pressure::DIMS,
            VolumetricFlowRate::DIMS,
            &[ConservationRole::Energy][..],
        ),
        (
            PortKind::ThermalTemperatureEntropy,
            Temperature::DIMS,
            Dims([2, 1, -3, -1, 0, 0]),
            &[ConservationRole::Energy][..],
        ),
    ];

    for (index, (kind, effort_dimensions, flow_dimensions, expected_roles)) in
        cases.into_iter().enumerate()
    {
        let a = scalar_schema(
            kind,
            &format!("port/scalar-{index}-a"),
            "frame/interface",
            PortOrientation::OutwardFromOwner,
            17,
        );
        let b = scalar_schema(
            kind,
            &format!("port/scalar-{index}-b"),
            "frame/interface",
            PortOrientation::OutwardFromOwner,
            17,
        );

        assert_eq!(a.version(), PORT_SCHEMA_VERSION, "schema={a:#?}");
        assert_eq!(a.kind(), kind, "schema={a:#?}");
        assert_eq!(a.effort_dimensions(), effort_dimensions, "schema={a:#?}");
        assert_eq!(a.flow_dimensions(), flow_dimensions, "schema={a:#?}");
        assert_eq!(
            a.effort_dimensions().checked_plus(a.flow_dimensions()),
            Some(Power::DIMS),
            "schema={a:#?}"
        );
        assert_eq!(a.shape(), PortValueShape::Scalar, "schema={a:#?}");
        assert_eq!(a.power_pairing(), PowerPairing::ScalarProduct);
        assert_eq!(a.coordinates().basis().as_str(), "basis/si-scalar");
        assert_eq!(a.coordinates().frame().as_str(), "frame/interface");
        assert_eq!(
            a.coordinates().orientation(),
            PortOrientation::OutwardFromOwner
        );
        assert_eq!(a.timestamp().clock().as_str(), "clock/coupling-window");
        assert_eq!(a.timestamp().tick(), 17);
        assert_eq!(a.conservation_roles(), expected_roles);

        let junction = ConservativeJunction::new(stable(&format!("junction/scalar-{index}")), a, b)
            .unwrap_or_else(|error| panic!("typed junction refused: {error:?}"));
        let typed = junction
            .interconnect_scalar(7.0, 3.0)
            .unwrap_or_else(|error| panic!("typed scalar seed refused: {error:?}"));
        let legacy = interconnect(kind, kind, 7.0, 3.0).unwrap();

        assert_eq!(
            typed.port_a.effort().to_bits(),
            legacy.port_a.effort().to_bits()
        );
        assert_eq!(
            typed.port_a.flow().to_bits(),
            legacy.port_a.flow().to_bits()
        );
        assert_eq!(
            typed.port_b.effort().to_bits(),
            legacy.port_b.effort().to_bits()
        );
        assert_eq!(
            typed.port_b.flow().to_bits(),
            legacy.port_b.flow().to_bits()
        );
        assert_eq!(
            typed.interface_power.to_bits(),
            legacy.interface_power.to_bits()
        );
        assert_eq!(
            typed.interface_power.to_bits(),
            0.0_f64.to_bits(),
            "junction={junction:#?}"
        );
    }
}

#[test]
fn extended_port_kinds_have_watt_pairings_and_kind_identity() {
    let cases = [
        (
            PortKind::RotationalTorqueAngularVelocity,
            Dims([2, 1, -2, 0, 0, 0]),
            Dims([0, 0, -1, 0, 0, 0]),
            &[ConservationRole::Energy, ConservationRole::AngularMomentum][..],
        ),
        (
            PortKind::ElectricalVoltageCurrent,
            Dims([2, 1, -3, 0, -1, 0]),
            Dims([0, 0, 0, 0, 1, 0]),
            &[ConservationRole::Energy, ConservationRole::ElectricCharge][..],
        ),
        (
            PortKind::MagneticMmfFluxRate,
            Dims([0, 0, 0, 0, 1, 0]),
            Dims([2, 1, -3, 0, -1, 0]),
            &[ConservationRole::Energy][..],
        ),
        (
            PortKind::ChemicalPotentialAmountFlow,
            Dims([2, 1, -2, 0, 0, -1]),
            Dims([0, 0, -1, 0, 0, 1]),
            &[ConservationRole::Energy, ConservationRole::Amount][..],
        ),
    ];

    for (index, (kind, effort, flow, expected_roles)) in cases.into_iter().enumerate() {
        let schema = scalar_schema(
            kind,
            &format!("port/extended-{index}"),
            "frame/extended",
            PortOrientation::OutwardFromOwner,
            52,
        );
        assert_eq!(schema.effort_dimensions(), effort, "schema={schema:#?}");
        assert_eq!(schema.flow_dimensions(), flow, "schema={schema:#?}");
        assert_eq!(
            effort.checked_plus(flow),
            Some(Power::DIMS),
            "schema={schema:#?}"
        );
        assert_eq!(schema.conservation_roles(), expected_roles);
    }

    let chemical_vector = PortSchema::try_new(
        stable("port/chemical-vector"),
        PortKind::ChemicalPotentialAmountFlow,
        PortKind::ChemicalPotentialAmountFlow.canonical_effort_dimensions(),
        PortKind::ChemicalPotentialAmountFlow.canonical_flow_dimensions(),
        PortValueShape::vector(4).unwrap(),
        CoordinateBinding::new(
            stable("basis/species-order-v1"),
            stable("frame/reactor"),
            PortOrientation::OutwardFromOwner,
        ),
        PowerPairing::EuclideanDot,
        PortTimestamp::new(stable("clock/reaction-window"), 52),
        [ConservationRole::Energy, ConservationRole::Amount],
    )
    .unwrap();
    assert!(matches!(
        chemical_vector.shape(),
        PortValueShape::Vector(components) if components.get() == 4
    ));
}

#[test]
fn extended_port_kinds_fail_closed_without_schema_admission() {
    let extended_kinds = [
        PortKind::RotationalTorqueAngularVelocity,
        PortKind::ElectricalVoltageCurrent,
        PortKind::MagneticMmfFluxRate,
        PortKind::ChemicalPotentialAmountFlow,
    ];
    for kind in extended_kinds {
        let raw_a = Port::new(1.0, 2.0, kind);
        let raw_b = Port::new(1.0, -2.0, kind);
        assert!(!raw_a.conjugate_to(&raw_b));
        assert_eq!(
            interconnect(kind, kind, 1.0, 2.0),
            Err(CoupleError::LegacyInterconnectRequiresSeedKind { kind })
        );
    }

    let wrong_kind = PortSchema::try_new(
        stable("port/electrical-with-mechanical-dimensions"),
        PortKind::ElectricalVoltageCurrent,
        Force::DIMS,
        Velocity::DIMS,
        PortValueShape::Scalar,
        CoordinateBinding::new(
            stable("basis/si-scalar"),
            stable("frame/electrical"),
            PortOrientation::OutwardFromOwner,
        ),
        PowerPairing::ScalarProduct,
        PortTimestamp::new(stable("clock/electrical"), 52),
        [ConservationRole::Energy],
    );
    assert!(matches!(
        wrong_kind,
        Err(CoupleError::PortKindDimensionMismatch {
            kind: PortKind::ElectricalVoltageCurrent,
            side: PortVariable::Effort,
            expected: Dims([2, 1, -3, 0, -1, 0]),
            actual: Dims([1, 1, -2, 0, 0, 0]),
        })
    ));

    let missing_charge = PortSchema::try_new(
        stable("port/electrical-missing-charge-role"),
        PortKind::ElectricalVoltageCurrent,
        PortKind::ElectricalVoltageCurrent.canonical_effort_dimensions(),
        PortKind::ElectricalVoltageCurrent.canonical_flow_dimensions(),
        PortValueShape::Scalar,
        CoordinateBinding::new(
            stable("basis/si-scalar"),
            stable("frame/electrical"),
            PortOrientation::OutwardFromOwner,
        ),
        PowerPairing::ScalarProduct,
        PortTimestamp::new(stable("clock/electrical"), 52),
        [ConservationRole::Energy],
    );
    assert_eq!(
        missing_charge,
        Err(CoupleError::MissingPortKindConservationRole {
            kind: PortKind::ElectricalVoltageCurrent,
            role: ConservationRole::ElectricCharge,
        })
    );
}

#[test]
fn field_measure_side_preserves_generalized_kind_dimensions() {
    let entropy_flux_density = Dims([0, 1, -3, -1, 0, 0]);
    let coordinates = CoordinateBinding::new(
        stable("basis/thermal-trace"),
        stable("frame/thermal-surface"),
        PortOrientation::OutwardFromOwner,
    );
    let timestamp = PortTimestamp::new(stable("clock/thermal-surface"), 61);
    let thermal_field = PortSchema::try_new(
        stable("port/thermal-surface"),
        PortKind::ThermalTemperatureEntropy,
        Temperature::DIMS,
        entropy_flux_density,
        PortValueShape::field(1, SpaceType::HGrad, SpaceType::HDiv).unwrap(),
        coordinates.clone(),
        PowerPairing::FieldDuality {
            measure_dimensions: Area::DIMS,
            measure_side: FieldMeasureSide::Flow,
        },
        timestamp.clone(),
        [ConservationRole::Energy, ConservationRole::Entropy],
    )
    .unwrap();
    assert_eq!(
        thermal_field.flow_dimensions(),
        entropy_flux_density,
        "schema={thermal_field:#?}"
    );

    let wrong_side = PortSchema::try_new(
        stable("port/thermal-surface-wrong-side"),
        PortKind::ThermalTemperatureEntropy,
        Temperature::DIMS,
        entropy_flux_density,
        PortValueShape::field(1, SpaceType::HGrad, SpaceType::HDiv).unwrap(),
        coordinates,
        PowerPairing::FieldDuality {
            measure_dimensions: Area::DIMS,
            measure_side: FieldMeasureSide::Effort,
        },
        timestamp,
        [ConservationRole::Energy, ConservationRole::Entropy],
    );
    assert!(matches!(
        wrong_side,
        Err(CoupleError::PortKindDimensionMismatch {
            kind: PortKind::ThermalTemperatureEntropy,
            side: PortVariable::Effort,
            ..
        })
    ));

    let measure_application_overflow = PortSchema::try_new(
        stable("port/measure-application-overflow"),
        PortKind::MechanicalForceVelocity,
        Dims([127, 0, 0, 0, 0, 0]),
        Dims([-127, 0, 0, 0, 0, 0]),
        PortValueShape::field(1, SpaceType::HDiv, SpaceType::HGrad).unwrap(),
        CoordinateBinding::new(
            stable("basis/overflow"),
            stable("frame/overflow"),
            PortOrientation::OutwardFromOwner,
        ),
        PowerPairing::FieldDuality {
            measure_dimensions: Power::DIMS,
            measure_side: FieldMeasureSide::Effort,
        },
        PortTimestamp::new(stable("clock/overflow"), 61),
        [ConservationRole::Energy],
    );
    assert!(matches!(
        measure_application_overflow,
        Err(CoupleError::PortMeasureApplicationOverflow {
            side: PortVariable::Effort,
            ..
        })
    ));
}

#[test]
fn added_mass_fixture_is_a_bitwise_v2_schema_migration_golden() {
    let junction = ConservativeJunction::new(
        stable("junction/added-mass-v2"),
        scalar_schema(
            PortKind::MechanicalForceVelocity,
            "port/structure-v2",
            "frame/fsi-interface",
            PortOrientation::OutwardFromOwner,
            31,
        ),
        scalar_schema(
            PortKind::MechanicalForceVelocity,
            "port/fluid-v2",
            "frame/fsi-interface",
            PortOrientation::OutwardFromOwner,
            31,
        ),
    )
    .unwrap();

    let fixed_legacy = iterate_fixed_relaxation(2.0, 3.0, 0.0, 1.0, 100, 1e-9);
    let fixed_v2 = junction
        .iterate_added_mass_fixed(2.0, 3.0, 0.0, 1.0, 100, 1e-9)
        .unwrap();
    assert_fsi_result_bits_eq(fixed_v2, fixed_legacy);

    let aitken_legacy = iterate_aitken(2.0, 3.0, 0.0, 0.5, 2.0, 100, 1e-9);
    let aitken_v2 = junction
        .iterate_added_mass_aitken(2.0, 3.0, 0.0, 0.5, 2.0, 100, 1e-9)
        .unwrap();
    assert_fsi_result_bits_eq(aitken_v2, aitken_legacy);

    let fluid = ConservativeJunction::new(
        stable("junction/not-added-mass"),
        scalar_schema(
            PortKind::FluidPressureFlux,
            "port/fluid-a",
            "frame/fluid-interface",
            PortOrientation::OutwardFromOwner,
            31,
        ),
        scalar_schema(
            PortKind::FluidPressureFlux,
            "port/fluid-b",
            "frame/fluid-interface",
            PortOrientation::OutwardFromOwner,
            31,
        ),
    )
    .unwrap();
    assert!(matches!(
        fluid.iterate_added_mass_fixed(2.0, 3.0, 0.0, 1.0, 100, 1e-9),
        Err(CoupleError::AddedMassFixtureRequiresMechanicalPort {
            kind: PortKind::FluidPressureFlux
        })
    ));
}

#[test]
#[allow(
    clippy::too_many_lines,
    reason = "one fail-closed matrix keeps shared metadata and error precedence visible"
)]
fn port_schema_v2_fails_closed_on_malformed_metadata() {
    assert!(matches!(
        StableId::new("contains whitespace"),
        Err(CoupleError::InvalidStableId { .. })
    ));
    assert!(matches!(
        PortValueShape::vector(0),
        Err(CoupleError::EmptyPortShape)
    ));

    let coordinates = CoordinateBinding::new(
        stable("basis/si"),
        stable("frame/common"),
        PortOrientation::OutwardFromOwner,
    );
    let timestamp = PortTimestamp::new(stable("clock/window"), 0);
    let wrong_power = PortSchema::try_new(
        stable("port/wrong-power"),
        PortKind::MechanicalForceVelocity,
        Dims::NONE,
        Dims::NONE,
        PortValueShape::Scalar,
        coordinates.clone(),
        PowerPairing::ScalarProduct,
        timestamp.clone(),
        [ConservationRole::Energy],
    );
    assert!(matches!(
        wrong_power,
        Err(CoupleError::PortPowerDimensionMismatch {
            product: Dims::NONE,
            ..
        })
    ));

    let overflow = PortSchema::try_new(
        stable("port/overflow"),
        PortKind::MechanicalForceVelocity,
        Dims([127, 0, 0, 0, 0, 0]),
        Dims([1, 0, 0, 0, 0, 0]),
        PortValueShape::Scalar,
        coordinates.clone(),
        PowerPairing::ScalarProduct,
        timestamp.clone(),
        [ConservationRole::Energy],
    );
    assert!(matches!(
        overflow,
        Err(CoupleError::PortDimensionOverflow { .. })
    ));

    let wrong_pairing = PortSchema::try_new(
        stable("port/wrong-pairing"),
        PortKind::MechanicalForceVelocity,
        Force::DIMS,
        Velocity::DIMS,
        PortValueShape::vector(3).unwrap(),
        coordinates.clone(),
        PowerPairing::ScalarProduct,
        timestamp.clone(),
        [ConservationRole::Energy],
    );
    assert!(matches!(
        wrong_pairing,
        Err(CoupleError::PairingShapeMismatch { .. })
    ));

    let missing_energy = PortSchema::try_new(
        stable("port/missing-energy-role"),
        PortKind::MechanicalForceVelocity,
        Force::DIMS,
        Velocity::DIMS,
        PortValueShape::Scalar,
        coordinates,
        PowerPairing::ScalarProduct,
        timestamp,
        [ConservationRole::Mass],
    );
    assert_eq!(
        missing_energy,
        Err(CoupleError::MissingEnergyConservationRole)
    );

    let canonical_roles = PortSchema::try_new(
        stable("port/canonical-roles"),
        PortKind::MechanicalForceVelocity,
        Force::DIMS,
        Velocity::DIMS,
        PortValueShape::Scalar,
        CoordinateBinding::new(
            stable("basis/si"),
            stable("frame/common"),
            PortOrientation::OutwardFromOwner,
        ),
        PowerPairing::ScalarProduct,
        PortTimestamp::new(stable("clock/window"), 0),
        [
            ConservationRole::Mass,
            ConservationRole::Energy,
            ConservationRole::Mass,
        ],
    )
    .unwrap();
    assert_eq!(
        canonical_roles.conservation_roles(),
        &[ConservationRole::Energy, ConservationRole::Mass]
    );
}

#[test]
fn conservative_junction_localizes_schema_mismatches() {
    let a = scalar_schema(
        PortKind::MechanicalForceVelocity,
        "port/a",
        "frame/shared",
        PortOrientation::AlongFrame,
        4,
    );
    let bad_orientation = scalar_schema(
        PortKind::MechanicalForceVelocity,
        "port/b",
        "frame/shared",
        PortOrientation::AlongFrame,
        4,
    );
    assert!(matches!(
        ConservativeJunction::new(
            stable("junction/bad-orientation"),
            a.clone(),
            bad_orientation
        ),
        Err(CoupleError::IncompatiblePortSchemas {
            field: "orientation",
            ..
        })
    ));

    let opposite_common_frame = scalar_schema(
        PortKind::MechanicalForceVelocity,
        "port/b",
        "frame/shared",
        PortOrientation::AgainstFrame,
        4,
    );
    assert!(matches!(
        ConservativeJunction::new(
            stable("junction/no-unproved-frame-pullback"),
            a,
            opposite_common_frame
        ),
        Err(CoupleError::IncompatiblePortSchemas {
            field: "orientation",
            ..
        })
    ));

    let outward_a = scalar_schema(
        PortKind::MechanicalForceVelocity,
        "port/outward-a",
        "frame/shared",
        PortOrientation::OutwardFromOwner,
        4,
    );
    let outward_b = scalar_schema(
        PortKind::MechanicalForceVelocity,
        "port/outward-b",
        "frame/shared",
        PortOrientation::OutwardFromOwner,
        4,
    );
    let junction =
        ConservativeJunction::new(stable("junction/good"), outward_a, outward_b).unwrap();
    assert!(matches!(
        junction.interconnect_scalar(f64::NAN, 1.0),
        Err(CoupleError::NonFinitePortValue { field: "effort" })
    ));
    assert!(matches!(
        junction.interconnect_scalar(f64::MAX, 2.0),
        Err(CoupleError::NonFinitePortValue { field: "power" })
    ));
}

#[test]
fn non_scalar_schema_is_typed_but_not_silently_run_by_scalar_seed() {
    let coordinates = CoordinateBinding::new(
        stable("basis/cartesian"),
        stable("frame/vector"),
        PortOrientation::OutwardFromOwner,
    );
    let timestamp = PortTimestamp::new(stable("clock/vector"), 9);
    let make = |id: &str| {
        PortSchema::try_new(
            stable(id),
            PortKind::MechanicalForceVelocity,
            Force::DIMS,
            Velocity::DIMS,
            PortValueShape::vector(3).unwrap(),
            coordinates.clone(),
            PowerPairing::EuclideanDot,
            timestamp.clone(),
            [ConservationRole::Energy, ConservationRole::LinearMomentum],
        )
        .unwrap()
    };
    let junction = ConservativeJunction::new(
        stable("junction/vector"),
        make("port/vector-a"),
        make("port/vector-b"),
    )
    .unwrap();
    assert!(matches!(
        junction.interconnect_scalar(1.0, 2.0),
        Err(CoupleError::ScalarOperationRequiresScalarPort {
            shape: PortValueShape::Vector(_),
            ..
        })
    ));

    let field = PortValueShape::field(3, SpaceType::HDiv, SpaceType::HGrad).unwrap();
    assert!(matches!(
        field,
        PortValueShape::Field {
            components,
            effort_space: SpaceType::HDiv,
            flow_space: SpaceType::HGrad,
        } if components.get() == 3
    ));

    let field_schema = PortSchema::try_new(
        stable("port/surface-traction"),
        PortKind::MechanicalForceVelocity,
        Pressure::DIMS,
        Velocity::DIMS,
        field,
        coordinates.clone(),
        PowerPairing::FieldDuality {
            measure_dimensions: Area::DIMS,
            measure_side: FieldMeasureSide::Effort,
        },
        timestamp.clone(),
        [ConservationRole::Energy, ConservationRole::LinearMomentum],
    )
    .unwrap();
    assert_eq!(
        field_schema.power_pairing(),
        PowerPairing::FieldDuality {
            measure_dimensions: Area::DIMS,
            measure_side: FieldMeasureSide::Effort,
        }
    );

    let missing_measure = PortSchema::try_new(
        stable("port/surface-traction-missing-measure"),
        PortKind::MechanicalForceVelocity,
        Pressure::DIMS,
        Velocity::DIMS,
        field,
        coordinates,
        PowerPairing::FieldDuality {
            measure_dimensions: Dims::NONE,
            measure_side: FieldMeasureSide::Effort,
        },
        timestamp,
        [ConservationRole::Energy, ConservationRole::LinearMomentum],
    );
    assert!(matches!(
        missing_measure,
        Err(CoupleError::PortPowerDimensionMismatch { .. })
    ));
}

#[test]
#[allow(
    clippy::too_many_lines,
    reason = "one primitive matrix keeps the four thermodynamic claims directly comparable"
)]
fn four_port_thermodynamic_primitives_keep_claims_separate() {
    let port = scalar_schema(
        PortKind::ThermalTemperatureEntropy,
        "port/thermal",
        "frame/thermal-boundary",
        PortOrientation::OutwardFromOwner,
        23,
    );
    let peer = scalar_schema(
        PortKind::ThermalTemperatureEntropy,
        "port/thermal-peer",
        "frame/thermal-boundary",
        PortOrientation::OutwardFromOwner,
        23,
    );
    let conservative =
        ConservativeJunction::new(stable("junction/thermal"), port.clone(), peer).unwrap();
    let storage = StorageElement::new(
        stable("storage/body"),
        port.clone(),
        StoragePotential::FreeEnergy,
        stable("state/body-temperature"),
        NonZeroUsize::new(4).unwrap(),
        stable("operator/free-energy-gradient-v1"),
    )
    .unwrap();
    let dissipation = DissipativeRelation::new(
        stable("dissipation/conduction"),
        port.clone(),
        DissipationLaw::Conductive,
        stable("operator/fourier-v1"),
        DissipationEvidence::NonnegativeProduction(stable("receipt/fourier-production")),
    )
    .unwrap();
    let mismatched_boundary = AccountingBoundary::new(
        stable("boundary/mismatched"),
        CoordinateBinding::new(
            stable("basis/other"),
            stable("frame/thermal-boundary"),
            PortOrientation::OutwardFromOwner,
        ),
        BoundaryTreatment::ExternalReservoirExchange,
    );
    assert!(matches!(
        SourceOrReservoir::new(
            stable("reservoir/mismatched"),
            port.clone(),
            SourceClass::Reservoir,
            mismatched_boundary,
        ),
        Err(CoupleError::AccountingBoundaryMismatch { .. })
    ));

    let boundary = AccountingBoundary::new(
        stable("boundary/ambient"),
        port.coordinates().clone(),
        BoundaryTreatment::ExternalReservoirExchange,
    );
    let source = SourceOrReservoir::new(
        stable("reservoir/ambient"),
        port,
        SourceClass::Reservoir,
        boundary,
    )
    .unwrap();

    assert_eq!(storage.potential(), StoragePotential::FreeEnergy);
    assert_eq!(storage.state_dimension().get(), 4);
    assert_eq!(
        storage.constitutive_gradient().as_str(),
        "operator/free-energy-gradient-v1"
    );
    assert!(matches!(
        dissipation.evidence(),
        DissipationEvidence::NonnegativeProduction(receipt)
            if receipt.as_str() == "receipt/fourier-production"
    ));
    assert_eq!(
        source.boundary().treatment(),
        BoundaryTreatment::ExternalReservoirExchange
    );
    assert_eq!(source.boundary().coordinates(), source.port().coordinates());

    let primitives = [
        PortPrimitive::ConservativeJunction(conservative),
        PortPrimitive::StorageElement(storage),
        PortPrimitive::DissipativeRelation(dissipation),
        PortPrimitive::SourceOrReservoir(source),
    ];
    assert!(matches!(
        &primitives[0],
        PortPrimitive::ConservativeJunction(_)
    ));
    assert!(matches!(&primitives[1], PortPrimitive::StorageElement(_)));
    assert!(matches!(
        &primitives[2],
        PortPrimitive::DissipativeRelation(_)
    ));
    assert!(matches!(
        &primitives[3],
        PortPrimitive::SourceOrReservoir(_)
    ));
}

#[test]
fn primitive_identity_aliases_fail_closed() {
    let port_a = scalar_schema(
        PortKind::MechanicalForceVelocity,
        "port/identity-a",
        "frame/identity",
        PortOrientation::OutwardFromOwner,
        41,
    );
    let port_b = scalar_schema(
        PortKind::MechanicalForceVelocity,
        "port/identity-b",
        "frame/identity",
        PortOrientation::OutwardFromOwner,
        41,
    );

    assert!(matches!(
        ConservativeJunction::new(port_a.id().clone(), port_a.clone(), port_b.clone()),
        Err(CoupleError::DuplicateIdentity { .. })
    ));
    assert!(matches!(
        ConservativeJunction::new(
            stable("junction/duplicate-port"),
            port_a.clone(),
            port_a.clone()
        ),
        Err(CoupleError::IncompatiblePortSchemas {
            field: "stable_id",
            ..
        })
    ));
    assert!(matches!(
        StorageElement::new(
            port_a.id().clone(),
            port_a.clone(),
            StoragePotential::Hamiltonian,
            stable("state/identity"),
            NonZeroUsize::new(1).unwrap(),
            stable("operator/identity-gradient"),
        ),
        Err(CoupleError::DuplicateIdentity { .. })
    ));
    assert!(matches!(
        DissipativeRelation::new(
            port_a.id().clone(),
            port_a.clone(),
            DissipationLaw::Frictional,
            stable("operator/identity-friction"),
            DissipationEvidence::Monotonicity(stable("receipt/identity-monotonicity")),
        ),
        Err(CoupleError::DuplicateIdentity { .. })
    ));

    let distinct_boundary = AccountingBoundary::new(
        stable("boundary/identity-distinct"),
        port_a.coordinates().clone(),
        BoundaryTreatment::IncludedSourceTerm,
    );
    assert!(matches!(
        SourceOrReservoir::new(
            port_a.id().clone(),
            port_a.clone(),
            SourceClass::PrescribedFlow,
            distinct_boundary,
        ),
        Err(CoupleError::DuplicateIdentity { .. })
    ));

    let aliased_boundary = AccountingBoundary::new(
        port_a.id().clone(),
        port_a.coordinates().clone(),
        BoundaryTreatment::IncludedSourceTerm,
    );
    assert!(matches!(
        SourceOrReservoir::new(
            stable("source/identity-distinct"),
            port_a,
            SourceClass::PrescribedFlow,
            aliased_boundary,
        ),
        Err(CoupleError::DuplicateIdentity { .. })
    ));

    let relation_boundary_id = stable("boundary/relation-alias");
    let relation_aliased_boundary = AccountingBoundary::new(
        relation_boundary_id.clone(),
        port_b.coordinates().clone(),
        BoundaryTreatment::IncludedSourceTerm,
    );
    assert!(matches!(
        SourceOrReservoir::new(
            relation_boundary_id,
            port_b,
            SourceClass::PrescribedEffort,
            relation_aliased_boundary,
        ),
        Err(CoupleError::DuplicateIdentity { .. })
    ));
}

#[test]
fn ports_are_power_conjugate_by_physical_type() {
    let force = Port::new(10.0, 2.0, PortKind::MechanicalForceVelocity);
    let force2 = Port::new(5.0, 1.0, PortKind::MechanicalForceVelocity);
    let pressure = Port::new(3.0, 4.0, PortKind::FluidPressureFlux);
    assert!((force.power() - 20.0).abs() < 1e-12); // effort × flow
    assert!(force.conjugate_to(&force2)); // same physical type
    assert!(!force.conjugate_to(&pressure)); // force can't couple to pressure
}

#[test]
fn the_dirac_interconnection_conserves_interface_power_exactly() {
    let c = interconnect(
        PortKind::MechanicalForceVelocity,
        PortKind::MechanicalForceVelocity,
        7.0,
        3.0,
    )
    .unwrap();
    // shared effort, opposite flow -> net interface power is exactly zero (G0).
    assert!(c.interface_power.abs() < 1e-15);
    assert!((c.port_a.effort() - c.port_b.effort()).abs() < 1e-15);
    assert!((c.port_a.flow() + c.port_b.flow()).abs() < 1e-15);
    // incompatible ports are refused at composition time.
    assert!(matches!(
        interconnect(
            PortKind::MechanicalForceVelocity,
            PortKind::FluidPressureFlux,
            1.0,
            1.0
        ),
        Err(CoupleError::IncompatiblePorts { .. })
    ));
}

#[test]
fn the_energy_audit_measures_passivity_and_alarms_on_generation() {
    let mut audit = EnergyAudit::new();
    // a correct interconnection conserves power.
    let good = interconnect(
        PortKind::FluidPressureFlux,
        PortKind::FluidPressureFlux,
        4.0,
        2.0,
    )
    .unwrap();
    audit.record(good.interface_power);
    assert!(audit.is_passive(1e-12));
    // a BROKEN coupling (both ports inject power) generates energy -> alarm.
    let broken = interface_power(&[
        Port::new(2.0, 1.0, PortKind::MechanicalForceVelocity),
        Port::new(2.0, 1.0, PortKind::MechanicalForceVelocity),
    ]);
    audit.record(broken);
    assert!(!audit.is_passive(1e-12));
    assert!((audit.max_generation() - 4.0).abs() < 1e-12);
}

#[test]
fn the_energy_audit_fails_closed_on_a_nan_interface_power() {
    // Regression: a NaN interface power is a hard numerical breakdown — the
    // worst thing the passivity audit exists to flag. `f64::max` drops NaN
    // (`f64::max(0.0, NaN) == 0.0`), so the old fold reported ZERO generation
    // and certified the blown-up coupling as passive — a false certificate.
    let mut audit = EnergyAudit::new();
    audit.record(0.0); // a clean, conserved exchange first
    assert!(audit.is_passive(1e-12), "a conserved exchange is passive");
    audit.record(f64::NAN); // then a diverged exchange
    assert!(
        audit.max_generation().is_nan(),
        "a NaN balance must poison the metric, not vanish"
    );
    assert!(
        !audit.is_passive(1e-12),
        "a NaN interface power must never certify as passive"
    );
    // An arbitrarily large tolerance cannot rescue a NaN, either.
    assert!(!audit.is_passive(f64::INFINITY));
}

#[test]
fn the_aitken_factor_follows_the_delta_squared_formula() {
    let mut a = AitkenRelaxation::new(0.5, 2.0);
    // first call returns the initial ω.
    assert!((a.next_omega(3.0) - 0.5).abs() < 1e-12);
    // ω₁ = −ω₀·r₀/(r₁−r₀) = −0.5·3/(−1.5−3) = 1/3.
    assert!((a.next_omega(-1.5) - 1.0 / 3.0).abs() < 1e-9);
}

#[test]
fn naive_staggering_diverges_where_aitken_stays_stable() {
    // dense fluid on a light structure: added-mass ratio μ = 2 (> 1).
    let (mu, c, x0) = (2.0, 3.0, 0.0);
    // naive Gauss-Seidel staggering (ω = 1) DIVERGES.
    let naive = iterate_fixed_relaxation(mu, c, x0, 1.0, 100, 1e-9);
    assert!(!naive.converged, "naive should diverge, got {naive:?}");
    // Aitken-relaxed coupling CONVERGES to the fixed point x* = c/(1+μ) = 1.
    let aitken = iterate_aitken(mu, c, x0, 0.5, 2.0, 100, 1e-9);
    assert!(aitken.converged);
    assert!((aitken.solution - 1.0).abs() < 1e-6);
    assert!(
        aitken.steps <= 5,
        "Aitken should converge fast, took {}",
        aitken.steps
    );
}

#[test]
fn aitken_accelerates_over_a_stable_fixed_relaxation() {
    let (mu, c, x0) = (2.0, 3.0, 0.0);
    // a stable but slower under-relaxation.
    let fixed = iterate_fixed_relaxation(mu, c, x0, 0.3, 200, 1e-12);
    let aitken = iterate_aitken(mu, c, x0, 0.5, 2.0, 200, 1e-12);
    assert!(fixed.converged && aitken.converged);
    assert!(
        aitken.steps <= fixed.steps,
        "Aitken {} !<= fixed {}",
        aitken.steps,
        fixed.steps
    );
}

#[test]
fn light_added_mass_converges_even_naively() {
    // μ < 1 (heavy structure): naive staggering is already stable.
    let r = iterate_fixed_relaxation(0.5, 3.0, 0.0, 1.0, 100, 1e-9);
    assert!(r.converged);
    assert!((r.solution - 2.0).abs() < 1e-6); // x* = 3/(1+0.5) = 2
}

#[test]
fn coupling_is_deterministic() {
    let a = iterate_aitken(2.0, 3.0, 0.0, 0.5, 2.0, 100, 1e-9);
    let b = iterate_aitken(2.0, 3.0, 0.0, 0.5, 2.0, 100, 1e-9);
    assert_eq!(a, b);
}
