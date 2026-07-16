//! Battery for port-Hamiltonian coupling (fs-couple). Covers power-conjugate
//! ports, the Dirac interconnection's exact power conservation, the energy
//! audit (passivity measured, not assumed), the Aitken relaxation factor, and
//! the load-bearing added-mass comparison: naive staggering diverges where
//! Aitken-relaxed coupling stays stable. The PortSchema v2 PR-1/PR-2,
//! StreamPort PR-3, and closed accounting-window PR-4 batteries pin typed
//! relations, chart admission, balance completeness, and retained red reports.

use core::num::NonZeroUsize;

use fs_couple::{
    ACCOUNTING_WINDOW_VERSION, AccountingBoundary, AccountingWindowInterval, AccountingWindowSpec,
    AitkenRelaxation, AmountFlowRate, BoundaryEntropyBreakdown, BoundaryTemperatureReference,
    BoundaryTreatment, ChemicalEnergyAccounting, ChemicalEnergyInput, ChemicalEnergyMode,
    ConjugatePotentialChart, ConjugatePotentialKind, ConservationRole, ConservativeJunction,
    CoordinateBinding, CoupleError, DeviatoricStressWork, DissipationEvidence, DissipationLaw,
    DissipativeRelation, EnergyAudit, Entropy, EntropyFlowRate, EulerLegendreCrosswalk,
    ExactIdentityProofRef, ExergyEnvironment, FieldMeasureSide, FsiResult,
    IntegratedWindowTransfer, InternalEnergyCauchyWorkChart, MovingStreamEnthalpyChart,
    PORT_SCHEMA_VERSION, Port, PortKind, PortOrientation, PortPrimitive, PortSchema, PortTimestamp,
    PortValueShape, PortVariable, PowerPairing, PressureWorkCrosswalk, STREAM_PORT_VERSION,
    SourceClass, SourceOrReservoir, SpecificEnergy, StableId, StorageElement, StoragePotential,
    StreamChartBinding, StreamConstituentFlow, StreamConstituentId, StreamEnergyChart,
    StreamIdentity, StreamKinematics, StreamPort, StreamStressWorkConvention, WindowAuditReport,
    WindowAuditTolerances, WindowBalance, WindowBalanceLaw, WindowBoundaryContribution,
    WindowChargeEvidenceSubject, WindowChargeSchema, WindowElementAmount, WindowElementSchema,
    WindowEndpoint, WindowEvidenceRef, WindowEvidenceRole, WindowInventorySnapshot,
    WindowManifestEntry, interconnect, interface_power, iterate_aitken, iterate_fixed_relaxation,
};
use fs_iface::SpaceType;
use fs_qty::chemistry::{ElementId, SpeciesId};
use fs_qty::{
    Amount, Area, Density, Dims, ElectricCharge, Energy, Force, Mass, MassFlowRate, Power,
    Pressure, Temperature, Velocity, VolumetricFlowRate,
};

fn stable(value: &str) -> StableId {
    StableId::new(value).unwrap_or_else(|error| panic!("invalid test ID {value:?}: {error:?}"))
}

fn species_id(id: &str) -> StreamConstituentId {
    StreamConstituentId::Species(SpeciesId::new(id).unwrap())
}

fn element_id(id: &str) -> StreamConstituentId {
    StreamConstituentId::Element(ElementId::new(id).unwrap())
}

fn try_stream_binding(
    port_id: &str,
    orientation: PortOrientation,
    constituent_axis: impl IntoIterator<Item = StreamConstituentId>,
) -> Result<StreamChartBinding, CoupleError> {
    StreamChartBinding::try_new(
        stable(port_id),
        stable("state/stream-v1"),
        stable(&format!("basis/{port_id}")),
        constituent_axis,
        stable("reference/chemical-v1"),
        CoordinateBinding::new(
            stable("basis/cartesian"),
            stable("frame/stream-boundary"),
            orientation,
        ),
        PortTimestamp::new(stable("clock/stream-window"), 71),
        stable("gravity/datum-v1"),
        StreamStressWorkConvention::CauchyTensionPositiveOutwardPower,
    )
}

fn stream_binding(
    port_id: &str,
    orientation: PortOrientation,
    constituent_axis: impl IntoIterator<Item = StreamConstituentId>,
) -> StreamChartBinding {
    try_stream_binding(port_id, orientation, constituent_axis).unwrap()
}

fn identity_proof(
    binding: &StreamChartBinding,
    identity: StreamIdentity,
    label: &str,
) -> ExactIdentityProofRef {
    ExactIdentityProofRef::new(
        identity,
        stable(&format!("receipt/{label}")),
        stable("verifier/stream-identities-v1"),
        stable(&format!("digest/{label}")),
        binding.clone(),
    )
}

fn embedded_chemical(binding: &StreamChartBinding) -> ChemicalEnergyAccounting {
    ChemicalEnergyAccounting::try_new(
        binding,
        ChemicalEnergyInput::IncludedInStatePotential {
            reference_state: binding.chemical_reference_state().clone(),
        },
    )
    .unwrap()
}

fn stream_kinematics(velocity: [f64; 3], gz: f64) -> StreamKinematics {
    StreamKinematics::try_new(velocity.map(Velocity::new), SpecificEnergy::new(gz)).unwrap()
}

fn deviatoric_work(binding: &StreamChartBinding, power_rate: f64) -> DeviatoricStressWork {
    DeviatoricStressWork::try_new(
        binding,
        Power::new(power_rate),
        stable("operator/deviatoric-work-v1"),
        stable("evidence/deviatoric-work-v1"),
    )
    .unwrap()
}

fn species_flow(id: &str, amount_flow: f64) -> StreamConstituentFlow {
    StreamConstituentFlow::try_new(species_id(id), AmountFlowRate::new(amount_flow)).unwrap()
}

fn canonical_moving_chart(binding: &StreamChartBinding) -> StreamEnergyChart {
    StreamEnergyChart::MovingStreamEnthalpy(Box::new(
        MovingStreamEnthalpyChart::try_new(
            binding.clone(),
            SpecificEnergy::new(10.0),
            stream_kinematics([3.0, 4.0, 0.0], 2.5),
            deviatoric_work(binding, 3.0),
            embedded_chemical(binding),
        )
        .unwrap(),
    ))
}

fn complete_stream(
    binding: StreamChartBinding,
    constituent_flows: impl IntoIterator<Item = StreamConstituentFlow>,
    energy_flow: f64,
    chart: StreamEnergyChart,
) -> Result<StreamPort, CoupleError> {
    StreamPort::try_new(
        binding,
        MassFlowRate::new(2.0),
        constituent_flows,
        [Force::new(1.0), Force::new(2.0), Force::new(3.0)],
        Power::new(energy_flow),
        EntropyFlowRate::new(0.5),
        chart,
    )
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

fn assert_f64_bits_eq(actual: f64, expected: f64) {
    assert_eq!(actual.to_bits(), expected.to_bits());
}

fn window_interval() -> AccountingWindowInterval {
    AccountingWindowInterval::try_new(
        stable("window/pr4"),
        PortTimestamp::new(stable("clock/coupling-window"), 10),
        PortTimestamp::new(stable("clock/coupling-window"), 12),
    )
    .unwrap()
}

fn window_evidence(
    label: &str,
    role: WindowEvidenceRole,
    interval: &AccountingWindowInterval,
) -> WindowEvidenceRef {
    WindowEvidenceRef::new(
        stable(&format!("receipt/window-{label}")),
        stable("verifier/window-v1"),
        stable(&format!("digest/window-{label}")),
        role,
        interval.clone(),
    )
}

fn window_closure_evidence(
    label: &str,
    interval: &AccountingWindowInterval,
    manifest: &[WindowManifestEntry],
) -> WindowEvidenceRef {
    window_evidence(
        label,
        WindowEvidenceRole::manifest_closure(manifest),
        interval,
    )
}

fn hydrogen_amount(value: f64) -> WindowElementAmount {
    WindowElementAmount::try_new(ElementId::new("H").unwrap(), Amount::new(value)).unwrap()
}

fn window_balance(
    energy: f64,
    mass: f64,
    hydrogen: f64,
    charge: f64,
    entropy: f64,
) -> WindowBalance {
    WindowBalance::try_new(
        Energy::new(energy),
        Mass::new(mass),
        [hydrogen_amount(hydrogen)],
        ElectricCharge::new(charge),
        Entropy::new(entropy),
    )
    .unwrap()
}

fn boundary_entropy(
    label: &str,
    interval: &AccountingWindowInterval,
    values: [f64; 4],
) -> BoundaryEntropyBreakdown {
    let temperature_reference = if values[2].to_bits() << 1 == 0 && values[3].to_bits() << 1 == 0 {
        BoundaryTemperatureReference::not_applicable(stable("rationale/no-heat-radiation"))
    } else {
        BoundaryTemperatureReference::try_constant(
            Temperature::new(400.0),
            window_evidence(
                &format!("temperature-{label}"),
                WindowEvidenceRole::BoundaryTemperature {
                    contribution_id: stable(&format!("contribution/{label}")),
                },
                interval,
            ),
        )
        .unwrap()
    };
    BoundaryEntropyBreakdown::try_new(
        Entropy::new(values[0]),
        Entropy::new(values[1]),
        Entropy::new(values[2]),
        Entropy::new(values[3]),
        temperature_reference,
    )
    .unwrap()
}

fn window_stream(label: &str) -> StreamPort {
    let binding = StreamChartBinding::try_new(
        stable(&format!("port/window-{label}")),
        stable("state/window-stream-v1"),
        stable("basis/window-species-v1"),
        [species_id("H")],
        stable("reference/window-chemical-v1"),
        CoordinateBinding::new(
            stable("basis/window-cartesian"),
            stable("frame/window"),
            PortOrientation::OutwardFromOwner,
        ),
        PortTimestamp::new(stable("clock/coupling-window"), 11),
        stable("gravity/window-datum-v1"),
        StreamStressWorkConvention::CauchyTensionPositiveOutwardPower,
    )
    .unwrap();
    let chart = canonical_moving_chart(&binding);
    complete_stream(binding, [species_flow("H", 1.0)], 53.0, chart).unwrap()
}

fn external_window_entry(label: &str, interval: &AccountingWindowInterval) -> WindowManifestEntry {
    let port = window_stream(label);
    let contribution_id = stable(&format!("contribution/{label}"));
    let boundary = AccountingBoundary::new(
        stable(&format!("boundary/window-{label}")),
        port.binding().coordinates().clone(),
        BoundaryTreatment::ExternalReservoirExchange,
    );
    WindowManifestEntry::external_stream_port(
        contribution_id.clone(),
        stable(&format!("exchange/{label}")),
        stable(&format!("reservoir/{label}")),
        &port,
        boundary,
        window_evidence(
            &format!("binding-{label}"),
            WindowEvidenceRole::ManifestBinding {
                contribution_id,
                local_port_id: port.binding().port_id().clone(),
            },
            interval,
        ),
    )
    .unwrap()
}

fn audited_window_elements(interval: &AccountingWindowInterval) -> WindowElementSchema {
    let elements = vec![ElementId::new("H").unwrap()];
    WindowElementSchema::try_audited(
        elements.clone(),
        window_evidence(
            "element-projection",
            WindowEvidenceRole::ElementProjection { elements },
            interval,
        ),
    )
    .unwrap()
}

fn direct_window_charge(interval: &AccountingWindowInterval) -> WindowChargeSchema {
    let basis = stable("basis/direct-coulomb-v1");
    WindowChargeSchema::DirectCoulomb {
        basis: basis.clone(),
        projection_evidence: window_evidence(
            "charge-projection",
            WindowEvidenceRole::ChargeProjection {
                subject: WindowChargeEvidenceSubject::DirectCoulomb { basis },
            },
            interval,
        ),
    }
}

fn zero_window_tolerances() -> WindowAuditTolerances {
    WindowAuditTolerances::try_new(
        Energy::new(0.0),
        Mass::new(0.0),
        Amount::new(0.0),
        ElectricCharge::new(0.0),
        Entropy::new(0.0),
    )
    .unwrap()
}

fn window_energy_reference() -> StableId {
    stable("reference/window-energy-v1")
}

fn window_spec(labels: &[&str], environment: Option<ExergyEnvironment>) -> AccountingWindowSpec {
    let interval = window_interval();
    let manifest: Vec<_> = labels
        .iter()
        .map(|label| external_window_entry(label, &interval))
        .collect();
    let closure_evidence = window_closure_evidence("closure", &interval, &manifest);
    AccountingWindowSpec::try_new(
        interval.clone(),
        window_energy_reference(),
        audited_window_elements(&interval),
        direct_window_charge(&interval),
        stable("convention/entropy-pr4"),
        manifest,
        zero_window_tolerances(),
        environment,
        closure_evidence,
    )
    .unwrap()
}

fn window_snapshot(
    spec: &AccountingWindowSpec,
    endpoint: WindowEndpoint,
    balance: WindowBalance,
) -> WindowInventorySnapshot {
    let label = match endpoint {
        WindowEndpoint::Initial => "initial-state",
        WindowEndpoint::Final => "final-state",
    };
    WindowInventorySnapshot::try_new(
        spec,
        endpoint,
        balance,
        window_evidence(
            label,
            WindowEvidenceRole::Inventory {
                endpoint,
                energy_reference: spec.energy_reference().clone(),
            },
            spec.interval(),
        ),
    )
    .unwrap()
}

fn window_contribution(
    spec: &AccountingWindowSpec,
    label: &str,
    values: [f64; 4],
    entropy: [f64; 4],
) -> WindowBoundaryContribution {
    let transfer = IntegratedWindowTransfer::try_new(
        Energy::new(values[0]),
        Mass::new(values[1]),
        [hydrogen_amount(values[2])],
        ElectricCharge::new(values[3]),
        boundary_entropy(label, spec.interval(), entropy),
        stable("convention/entropy-pr4"),
    )
    .unwrap();
    let contribution_id = stable(&format!("contribution/{label}"));
    let local_port_id = spec
        .manifest()
        .iter()
        .find(|entry| entry.contribution_id() == &contribution_id)
        .map(|entry| entry.local_port_id().clone())
        .unwrap();
    WindowBoundaryContribution::try_new(
        spec,
        &contribution_id,
        transfer,
        window_evidence(
            &format!("integration-{label}"),
            WindowEvidenceRole::IntegratedTransfer {
                contribution_id: contribution_id.clone(),
                local_port_id,
                energy_reference: spec.energy_reference().clone(),
            },
            spec.interval(),
        ),
    )
    .unwrap()
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
fn stream_port_moving_enthalpy_carries_the_complete_bundle() {
    let binding = stream_binding(
        "stream/canonical",
        PortOrientation::OutwardFromOwner,
        [species_id("O2"), element_id("O")],
    );
    let element_flow = StreamConstituentFlow::try_new(
        StreamConstituentId::Element(ElementId::new("O").unwrap()),
        AmountFlowRate::new(0.25),
    )
    .unwrap();
    let chart = canonical_moving_chart(&binding);
    let stream = complete_stream(
        binding,
        [element_flow, species_flow("O2", 0.5)],
        53.0,
        chart,
    )
    .unwrap();

    assert_eq!(stream.version(), STREAM_PORT_VERSION);
    assert_eq!(stream.mass_flow().value().to_bits(), 2.0_f64.to_bits());
    assert_eq!(stream.energy_flow().value().to_bits(), 53.0_f64.to_bits());
    assert_eq!(stream.constituent_flows().len(), 2);
    assert!(matches!(
        stream.constituent_flows()[0].id(),
        StreamConstituentId::Species(id) if id.as_str() == "O2"
    ));
    assert_eq!(
        stream.conservation_roles(),
        &[
            ConservationRole::Energy,
            ConservationRole::Mass,
            ConservationRole::Amount,
            ConservationRole::LinearMomentum,
            ConservationRole::Entropy,
        ]
    );
    assert!(matches!(
        stream.energy_chart(),
        StreamEnergyChart::MovingStreamEnthalpy(_)
    ));
}

#[test]
fn moving_enthalpy_admits_one_proved_explicit_chemical_term() {
    let binding = stream_binding(
        "stream/explicit-chemistry",
        PortOrientation::OutwardFromOwner,
        [species_id("H2")],
    );
    assert!(matches!(
        ChemicalEnergyAccounting::try_new(
            &binding,
            ChemicalEnergyInput::ExplicitSpeciesPotentials {
                reference_state: stable("reference/foreign-chemistry"),
                power_rate: Power::new(2.0),
                partition_proof: identity_proof(
                    &binding,
                    StreamIdentity::ChemicalEnergyPartition,
                    "foreign-reference",
                ),
            },
        ),
        Err(CoupleError::StreamChartBindingMismatch {
            field: "chemical_reference_state",
        })
    ));
    assert!(matches!(
        ChemicalEnergyAccounting::try_new(
            &binding,
            ChemicalEnergyInput::ExplicitSpeciesPotentials {
                reference_state: binding.chemical_reference_state().clone(),
                power_rate: Power::new(2.0),
                partition_proof: identity_proof(
                    &binding,
                    StreamIdentity::PressureWork,
                    "wrong-chemical-proof-kind",
                ),
            },
        ),
        Err(CoupleError::WrongStreamIdentityProof {
            expected: StreamIdentity::ChemicalEnergyPartition,
            actual: StreamIdentity::PressureWork,
        })
    ));
    let foreign_binding = stream_binding(
        "stream/foreign-explicit-chemistry",
        PortOrientation::OutwardFromOwner,
        [species_id("H2")],
    );
    assert!(matches!(
        ChemicalEnergyAccounting::try_new(
            &binding,
            ChemicalEnergyInput::ExplicitSpeciesPotentials {
                reference_state: binding.chemical_reference_state().clone(),
                power_rate: Power::new(2.0),
                partition_proof: identity_proof(
                    &foreign_binding,
                    StreamIdentity::ChemicalEnergyPartition,
                    "foreign-chemical-binding",
                ),
            },
        ),
        Err(CoupleError::StreamChartBindingMismatch { field: "port_id" })
    ));

    let chemical = ChemicalEnergyAccounting::try_new(
        &binding,
        ChemicalEnergyInput::ExplicitSpeciesPotentials {
            reference_state: binding.chemical_reference_state().clone(),
            power_rate: Power::new(2.0),
            partition_proof: identity_proof(
                &binding,
                StreamIdentity::ChemicalEnergyPartition,
                "explicit-chemical",
            ),
        },
    )
    .unwrap();
    assert_eq!(
        chemical.mode(),
        ChemicalEnergyMode::ExplicitSpeciesPotentials
    );
    let chart = StreamEnergyChart::MovingStreamEnthalpy(Box::new(
        MovingStreamEnthalpyChart::try_new(
            binding.clone(),
            SpecificEnergy::new(10.0),
            stream_kinematics([0.0; 3], 0.0),
            deviatoric_work(&binding, 0.0),
            chemical,
        )
        .unwrap(),
    ));
    let stream = complete_stream(binding, [species_flow("H2", 0.5)], 22.0, chart).unwrap();
    assert_eq!(stream.energy_flow().value().to_bits(), 22.0_f64.to_bits());
}

#[test]
fn internal_energy_chart_exactly_crosswalks_pressure_work() {
    let binding = stream_binding(
        "stream/internal",
        PortOrientation::OutwardFromOwner,
        [species_id("N2")],
    );
    let crosswalk = PressureWorkCrosswalk::try_new(
        identity_proof(&binding, StreamIdentity::PressureWork, "pressure-work"),
        MassFlowRate::new(2.0),
        Density::new(2.0),
        SpecificEnergy::new(8.0),
        SpecificEnergy::new(4.0),
        Pressure::new(8.0),
        VolumetricFlowRate::new(1.0),
    )
    .unwrap();
    let internal_chart = StreamEnergyChart::InternalEnergyCauchyWork(Box::new(
        InternalEnergyCauchyWorkChart::try_new(
            binding.clone(),
            crosswalk,
            stream_kinematics([0.0; 3], 0.0),
            deviatoric_work(&binding, 0.0),
            embedded_chemical(&binding),
        )
        .unwrap(),
    ));
    let internal =
        complete_stream(binding, [species_flow("N2", 1.0)], 16.0, internal_chart).unwrap();

    let enthalpy_binding = stream_binding(
        "stream/enthalpy-equivalent",
        PortOrientation::OutwardFromOwner,
        [species_id("N2")],
    );
    let enthalpy_chart = StreamEnergyChart::MovingStreamEnthalpy(Box::new(
        MovingStreamEnthalpyChart::try_new(
            enthalpy_binding.clone(),
            SpecificEnergy::new(8.0),
            stream_kinematics([0.0; 3], 0.0),
            deviatoric_work(&enthalpy_binding, 0.0),
            embedded_chemical(&enthalpy_binding),
        )
        .unwrap(),
    ));
    let enthalpy = complete_stream(
        enthalpy_binding,
        [species_flow("N2", 1.0)],
        16.0,
        enthalpy_chart,
    )
    .unwrap();
    assert_eq!(internal.energy_flow(), enthalpy.energy_flow());
}

#[test]
fn pressure_work_crosswalk_refuses_inexact_or_foreign_evidence() {
    let binding = stream_binding(
        "stream/pressure-refusal",
        PortOrientation::OutwardFromOwner,
        [species_id("N2")],
    );
    let inexact_enthalpy = SpecificEnergy::new(f64::from_bits(8.0_f64.to_bits() + 1));
    assert!(matches!(
        PressureWorkCrosswalk::try_new(
            identity_proof(&binding, StreamIdentity::PressureWork, "pressure-inexact"),
            MassFlowRate::new(2.0),
            Density::new(2.0),
            inexact_enthalpy,
            SpecificEnergy::new(4.0),
            Pressure::new(8.0),
            VolumetricFlowRate::new(1.0),
        ),
        Err(CoupleError::StreamIdentityMismatch {
            identity: StreamIdentity::PressureWork,
            check: "specific_enthalpy",
        })
    ));

    let foreign_binding = stream_binding(
        "stream/foreign",
        PortOrientation::OutwardFromOwner,
        [species_id("N2")],
    );
    let foreign_crosswalk = PressureWorkCrosswalk::try_new(
        identity_proof(
            &foreign_binding,
            StreamIdentity::PressureWork,
            "pressure-foreign",
        ),
        MassFlowRate::new(2.0),
        Density::new(2.0),
        SpecificEnergy::new(8.0),
        SpecificEnergy::new(4.0),
        Pressure::new(8.0),
        VolumetricFlowRate::new(1.0),
    )
    .unwrap();
    assert!(matches!(
        InternalEnergyCauchyWorkChart::try_new(
            binding.clone(),
            foreign_crosswalk,
            stream_kinematics([0.0; 3], 0.0),
            deviatoric_work(&binding, 0.0),
            embedded_chemical(&binding),
        ),
        Err(CoupleError::StreamChartBindingMismatch { field: "port_id" })
    ));

    let local_crosswalk = PressureWorkCrosswalk::try_new(
        identity_proof(&binding, StreamIdentity::PressureWork, "pressure-local"),
        MassFlowRate::new(2.0),
        Density::new(2.0),
        SpecificEnergy::new(8.0),
        SpecificEnergy::new(4.0),
        Pressure::new(8.0),
        VolumetricFlowRate::new(1.0),
    )
    .unwrap();
    assert!(matches!(
        InternalEnergyCauchyWorkChart::try_new(
            binding.clone(),
            local_crosswalk,
            stream_kinematics([0.0; 3], 0.0),
            deviatoric_work(&foreign_binding, 0.0),
            embedded_chemical(&binding),
        ),
        Err(CoupleError::StreamChartBindingMismatch { field: "port_id" })
    ));
}

#[test]
fn pressure_work_crosswalk_refuses_wrong_identity_and_nonpositive_density() {
    let binding = stream_binding(
        "stream/pressure-proof-refusal",
        PortOrientation::OutwardFromOwner,
        [species_id("N2")],
    );
    assert_eq!(
        PressureWorkCrosswalk::try_new(
            identity_proof(
                &binding,
                StreamIdentity::ChemicalEnergyPartition,
                "pressure-wrong-kind",
            ),
            MassFlowRate::new(2.0),
            Density::new(2.0),
            SpecificEnergy::new(8.0),
            SpecificEnergy::new(4.0),
            Pressure::new(8.0),
            VolumetricFlowRate::new(1.0),
        ),
        Err(CoupleError::WrongStreamIdentityProof {
            expected: StreamIdentity::PressureWork,
            actual: StreamIdentity::ChemicalEnergyPartition,
        })
    );
    assert_eq!(
        PressureWorkCrosswalk::try_new(
            identity_proof(
                &binding,
                StreamIdentity::PressureWork,
                "pressure-zero-density",
            ),
            MassFlowRate::new(2.0),
            Density::new(0.0),
            SpecificEnergy::new(8.0),
            SpecificEnergy::new(4.0),
            Pressure::new(8.0),
            VolumetricFlowRate::new(1.0),
        ),
        Err(CoupleError::NonPositiveStreamDensity)
    );
    assert!(matches!(
        PressureWorkCrosswalk::try_new(
            identity_proof(
                &binding,
                StreamIdentity::PressureWork,
                "pressure-rate-overflow",
            ),
            MassFlowRate::new(2.0),
            Density::new(1.0),
            SpecificEnergy::new(1.0e308),
            SpecificEnergy::new(0.0),
            Pressure::new(1.0e308),
            VolumetricFlowRate::new(2.0),
        ),
        Err(CoupleError::NonFiniteStreamValue {
            field: "pressure_work_rate",
            index: None,
        })
    ));
}

#[test]
fn conjugate_chart_requires_exact_euler_legendre_and_single_chemical_owner() {
    let binding = stream_binding(
        "stream/conjugate",
        PortOrientation::OutwardFromOwner,
        [species_id("CO2")],
    );
    let crosswalk = EulerLegendreCrosswalk::try_new(
        identity_proof(&binding, StreamIdentity::EulerLegendre, "euler-legendre"),
        SpecificEnergy::new(17.0),
        SpecificEnergy::new(8.0),
        SpecificEnergy::new(2.0),
        SpecificEnergy::new(11.0),
        SpecificEnergy::new(19.0),
        SpecificEnergy::new(9.0),
        SpecificEnergy::new(11.0),
    )
    .unwrap();
    for selected_potential in [
        ConjugatePotentialKind::InternalEnergy,
        ConjugatePotentialKind::Enthalpy,
        ConjugatePotentialKind::Helmholtz,
        ConjugatePotentialKind::Gibbs,
    ] {
        let chart = StreamEnergyChart::ConjugatePotential(Box::new(
            ConjugatePotentialChart::try_new(
                binding.clone(),
                crosswalk.clone(),
                selected_potential,
                stream_kinematics([0.0; 3], 0.0),
                deviatoric_work(&binding, 0.0),
                embedded_chemical(&binding),
            )
            .unwrap(),
        ));
        let stream =
            complete_stream(binding.clone(), [species_flow("CO2", 0.125)], 38.0, chart).unwrap();
        assert_eq!(stream.energy_flow().value().to_bits(), 38.0_f64.to_bits());
    }

    let partition = identity_proof(
        &binding,
        StreamIdentity::ChemicalEnergyPartition,
        "chemical-partition",
    );
    let explicit = ChemicalEnergyAccounting::try_new(
        &binding,
        ChemicalEnergyInput::ExplicitSpeciesPotentials {
            reference_state: binding.chemical_reference_state().clone(),
            power_rate: Power::new(2.0),
            partition_proof: partition.clone(),
        },
    )
    .unwrap();
    assert_eq!(
        ConjugatePotentialChart::try_new(
            binding.clone(),
            crosswalk,
            ConjugatePotentialKind::Gibbs,
            stream_kinematics([0.0; 3], 0.0),
            deviatoric_work(&binding, 0.0),
            explicit,
        ),
        Err(CoupleError::DoubleCountedChemicalEnergy)
    );
    assert!(matches!(
        ChemicalEnergyAccounting::try_new(
            &binding,
            ChemicalEnergyInput::IncludedAndExplicitSpeciesPotentials {
                reference_state: binding.chemical_reference_state().clone(),
                power_rate: Power::new(2.0),
                partition_proof: partition,
            },
        ),
        Err(CoupleError::DoubleCountedChemicalEnergy)
    ));

    let wrong_gibbs = SpecificEnergy::new(f64::from_bits(11.0_f64.to_bits() + 1));
    assert!(matches!(
        EulerLegendreCrosswalk::try_new(
            identity_proof(&binding, StreamIdentity::EulerLegendre, "euler-wrong-gibbs"),
            SpecificEnergy::new(17.0),
            SpecificEnergy::new(8.0),
            SpecificEnergy::new(2.0),
            SpecificEnergy::new(11.0),
            SpecificEnergy::new(19.0),
            SpecificEnergy::new(9.0),
            wrong_gibbs,
        ),
        Err(CoupleError::StreamIdentityMismatch {
            identity: StreamIdentity::EulerLegendre,
            check: "gibbs_legendre",
        })
    ));
}

#[test]
fn stream_constituent_and_chemical_basis_binding_fails_closed() {
    assert_eq!(
        try_stream_binding("stream/empty-axis", PortOrientation::OutwardFromOwner, [],),
        Err(CoupleError::EmptyStreamConstituents)
    );
    assert!(matches!(
        try_stream_binding(
            "stream/duplicate-axis",
            PortOrientation::OutwardFromOwner,
            [species_id("H2"), species_id("H2")],
        ),
        Err(CoupleError::DuplicateStreamConstituent { .. })
    ));

    let binding = stream_binding(
        "stream/axis-refusal",
        PortOrientation::OutwardFromOwner,
        [species_id("H2")],
    );
    assert_eq!(
        complete_stream(
            binding.clone(),
            [species_flow("O2", 0.5)],
            53.0,
            canonical_moving_chart(&binding),
        ),
        Err(CoupleError::StreamConstituentAxisMismatch {
            expected: vec!["species:H2".to_string()],
            actual: vec!["species:O2".to_string()],
        })
    );

    let element_binding = stream_binding(
        "stream/element-chemistry-refusal",
        PortOrientation::OutwardFromOwner,
        [element_id("H")],
    );
    assert_eq!(
        ChemicalEnergyAccounting::try_new(
            &element_binding,
            ChemicalEnergyInput::ExplicitSpeciesPotentials {
                reference_state: element_binding.chemical_reference_state().clone(),
                power_rate: Power::new(1.0),
                partition_proof: identity_proof(
                    &element_binding,
                    StreamIdentity::ChemicalEnergyPartition,
                    "element-explicit-chemical",
                ),
            },
        ),
        Err(CoupleError::StreamChemicalEnergyRequiresSpeciesAxis {
            constituent: "element:H".to_string(),
        })
    );

    let element_euler = EulerLegendreCrosswalk::try_new(
        identity_proof(
            &element_binding,
            StreamIdentity::EulerLegendre,
            "element-euler",
        ),
        SpecificEnergy::new(17.0),
        SpecificEnergy::new(8.0),
        SpecificEnergy::new(2.0),
        SpecificEnergy::new(11.0),
        SpecificEnergy::new(19.0),
        SpecificEnergy::new(9.0),
        SpecificEnergy::new(11.0),
    )
    .unwrap();
    assert_eq!(
        ConjugatePotentialChart::try_new(
            element_binding.clone(),
            element_euler,
            ConjugatePotentialKind::Enthalpy,
            stream_kinematics([0.0; 3], 0.0),
            deviatoric_work(&element_binding, 0.0),
            embedded_chemical(&element_binding),
        ),
        Err(CoupleError::StreamChemicalEnergyRequiresSpeciesAxis {
            constituent: "element:H".to_string(),
        })
    );
}

#[test]
fn stream_chart_terms_reject_nonfinite_values() {
    let binding = stream_binding(
        "stream/nonfinite-refusal",
        PortOrientation::OutwardFromOwner,
        [species_id("H2")],
    );
    assert!(matches!(
        StreamKinematics::try_new(
            [
                Velocity::new(0.0),
                Velocity::new(f64::NAN),
                Velocity::new(0.0)
            ],
            SpecificEnergy::new(0.0),
        ),
        Err(CoupleError::NonFiniteStreamValue {
            field: "velocity",
            index: Some(1),
        })
    ));
    assert!(matches!(
        StreamKinematics::try_new(
            [
                Velocity::new(1.0e308),
                Velocity::new(0.0),
                Velocity::new(0.0)
            ],
            SpecificEnergy::new(0.0),
        ),
        Err(CoupleError::NonFiniteStreamValue {
            field: "specific_transport_energy",
            index: None,
        })
    ));
    assert!(matches!(
        DeviatoricStressWork::try_new(
            &binding,
            Power::new(f64::INFINITY),
            stable("operator/nonfinite"),
            stable("evidence/nonfinite"),
        ),
        Err(CoupleError::NonFiniteStreamValue {
            field: "deviatoric_stress_work",
            index: None,
        })
    ));
    assert!(matches!(
        MovingStreamEnthalpyChart::try_new(
            binding.clone(),
            SpecificEnergy::new(f64::NAN),
            stream_kinematics([0.0; 3], 0.0),
            deviatoric_work(&binding, 0.0),
            embedded_chemical(&binding),
        ),
        Err(CoupleError::NonFiniteStreamValue {
            field: "specific_enthalpy",
            index: None,
        })
    ));
}

#[test]
fn stream_bundle_rates_reject_nonfinite_values() {
    let binding = stream_binding(
        "stream/nonfinite-bundle-refusal",
        PortOrientation::OutwardFromOwner,
        [species_id("H2")],
    );
    assert!(matches!(
        StreamPort::try_new(
            binding.clone(),
            MassFlowRate::new(2.0),
            [species_flow("H2", 0.5)],
            [Force::new(0.0), Force::new(f64::NAN), Force::new(0.0)],
            Power::new(53.0),
            EntropyFlowRate::new(0.5),
            canonical_moving_chart(&binding),
        ),
        Err(CoupleError::NonFiniteStreamValue {
            field: "momentum_flow",
            index: Some(1),
        })
    ));
    assert!(matches!(
        StreamPort::try_new(
            binding.clone(),
            MassFlowRate::new(2.0),
            [species_flow("H2", 0.5)],
            [Force::new(0.0); 3],
            Power::new(f64::NAN),
            EntropyFlowRate::new(0.5),
            canonical_moving_chart(&binding),
        ),
        Err(CoupleError::NonFiniteStreamValue {
            field: "energy_flow",
            index: None,
        })
    ));
    assert!(matches!(
        StreamPort::try_new(
            binding.clone(),
            MassFlowRate::new(2.0),
            [species_flow("H2", 0.5)],
            [Force::new(0.0); 3],
            Power::new(53.0),
            EntropyFlowRate::new(f64::NEG_INFINITY),
            canonical_moving_chart(&binding),
        ),
        Err(CoupleError::NonFiniteStreamValue {
            field: "entropy_flow",
            index: None,
        })
    ));
    assert!(matches!(
        StreamConstituentFlow::try_new(species_id("H2"), AmountFlowRate::new(f64::NAN)),
        Err(CoupleError::NonFiniteStreamValue {
            field: "constituent_amount_flow",
            index: None,
        })
    ));
}

#[test]
fn stream_bundle_fails_closed_on_missing_duplicate_or_nonfinite_rates() {
    let binding = stream_binding(
        "stream/refusal",
        PortOrientation::OutwardFromOwner,
        [species_id("H2")],
    );
    assert_eq!(
        complete_stream(binding.clone(), [], 53.0, canonical_moving_chart(&binding),),
        Err(CoupleError::EmptyStreamConstituents)
    );
    assert!(matches!(
        complete_stream(
            binding.clone(),
            [species_flow("H2", 0.5), species_flow("H2", 0.25)],
            53.0,
            canonical_moving_chart(&binding),
        ),
        Err(CoupleError::DuplicateStreamConstituent { .. })
    ));
    assert!(matches!(
        StreamPort::try_new(
            binding.clone(),
            MassFlowRate::new(f64::NAN),
            [species_flow("H2", 0.5)],
            [Force::new(0.0); 3],
            Power::new(53.0),
            EntropyFlowRate::new(0.5),
            canonical_moving_chart(&binding),
        ),
        Err(CoupleError::NonFiniteStreamValue {
            field: "mass_flow",
            index: None,
        })
    ));
    assert!(matches!(
        complete_stream(
            binding.clone(),
            [species_flow("H2", 0.5)],
            54.0,
            canonical_moving_chart(&binding),
        ),
        Err(CoupleError::StreamEnergyFlowMismatch { .. })
    ));

    let along = stream_binding(
        "stream/along",
        PortOrientation::AlongFrame,
        [species_id("H2")],
    );
    assert_eq!(
        complete_stream(
            along.clone(),
            [species_flow("H2", 0.5)],
            53.0,
            canonical_moving_chart(&along),
        ),
        Err(CoupleError::StreamPortRequiresOutwardOrientation {
            actual: PortOrientation::AlongFrame,
        })
    );
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
fn accounting_window_closes_all_pr4_balances_and_retains_exergy() {
    let environment =
        ExergyEnvironment::try_new(stable("environment/reference"), Temperature::new(300.0))
            .unwrap();
    let spec = window_spec(&["inlet", "outlet", "source"], Some(environment));
    let initial = window_snapshot(
        &spec,
        WindowEndpoint::Initial,
        window_balance(100.0, 10.0, 10.0, 0.0, 20.0),
    );
    let final_snapshot = window_snapshot(
        &spec,
        WindowEndpoint::Final,
        window_balance(115.0, 13.0, 13.0, 1.0, 24.0),
    );
    let inlet = window_contribution(
        &spec,
        "inlet",
        [20.0, 3.0, 3.0, 1.0],
        [1.5, -0.25, 0.75, -0.25],
    );
    let outlet = window_contribution(
        &spec,
        "outlet",
        [-10.0, -2.0, -2.0, -0.5],
        [-0.75, 0.0, 0.0, 0.0],
    );
    let source = window_contribution(
        &spec,
        "source",
        [5.0, 2.0, 2.0, 0.5],
        [1.0, 0.5, 0.25, 0.25],
    );

    let report =
        WindowAuditReport::audit(spec, initial, final_snapshot, [source, outlet, inlet]).unwrap();

    assert_eq!(report.spec().version(), ACCOUNTING_WINDOW_VERSION);
    assert!(report.is_green(), "report={report:#?}");
    assert_eq!(
        report.residuals().first_law().value().to_bits(),
        0.0_f64.to_bits()
    );
    assert_eq!(
        report.residuals().mass().value().to_bits(),
        0.0_f64.to_bits()
    );
    assert_eq!(
        report.residuals().charge().value().to_bits(),
        0.0_f64.to_bits()
    );
    assert_eq!(
        report.residuals().elements()[0].amount().value().to_bits(),
        0.0_f64.to_bits()
    );
    assert_f64_bits_eq(report.residuals().entropy_generation().value(), 1.0);
    assert_f64_bits_eq(report.exergy().unwrap().destruction().value(), 300.0);

    let contribution_ids: Vec<_> = report
        .contributions()
        .iter()
        .map(|entry| entry.manifest_entry().contribution_id().as_str())
        .collect();
    assert_eq!(
        contribution_ids,
        [
            "contribution/inlet",
            "contribution/outlet",
            "contribution/source"
        ]
    );
    assert_f64_bits_eq(
        report.contributions()[0].transfer().energy_into().value(),
        20.0,
    );
    assert_f64_bits_eq(
        report.contributions()[0]
            .transfer()
            .entropy_into()
            .total()
            .value(),
        1.75,
    );
    assert!(matches!(
        report.contributions()[0]
            .transfer()
            .entropy_into()
            .temperature_reference(),
        BoundaryTemperatureReference::Constant { .. }
    ));
}

#[test]
fn physical_window_failures_return_a_retained_red_report() {
    let spec = window_spec(&["inlet", "outlet", "source"], None);
    let initial = window_snapshot(
        &spec,
        WindowEndpoint::Initial,
        window_balance(100.0, 10.0, 10.0, 0.0, 20.0),
    );
    let final_snapshot = window_snapshot(
        &spec,
        WindowEndpoint::Final,
        window_balance(116.0, 13.0, 13.0, 1.0, 22.0),
    );
    let contributions = [
        window_contribution(
            &spec,
            "inlet",
            [20.0, 3.0, 3.0, 1.0],
            [1.5, -0.25, 0.75, -0.25],
        ),
        window_contribution(
            &spec,
            "outlet",
            [-10.0, -2.0, -2.0, -0.5],
            [-0.75, 0.0, 0.0, 0.0],
        ),
        window_contribution(
            &spec,
            "source",
            [5.0, 2.0, 2.0, 0.5],
            [1.0, 0.5, 0.25, 0.25],
        ),
    ];

    let report = WindowAuditReport::audit(spec, initial, final_snapshot, contributions).unwrap();
    assert!(!report.is_green());
    assert_f64_bits_eq(report.residuals().first_law().value(), 1.0);
    assert_f64_bits_eq(report.residuals().entropy_generation().value(), -1.0);
    assert!(
        report
            .violations()
            .iter()
            .any(|violation| matches!(violation.law(), WindowBalanceLaw::FirstLaw))
    );
    let entropy_violation = report
        .violations()
        .iter()
        .find(|violation| matches!(violation.law(), WindowBalanceLaw::EntropyGeneration))
        .unwrap();
    assert_f64_bits_eq(entropy_violation.residual_si(), -1.0);
    assert_eq!(entropy_violation.contribution_ids().len(), 3);
    assert_eq!(entropy_violation.port_ids().len(), 3);
    assert_eq!(entropy_violation.exchange_ids().len(), 3);
}

#[test]
fn accounting_window_admission_refuses_missing_or_foreign_evidence() {
    let spec = window_spec(&["inlet", "source"], None);
    let initial = window_snapshot(
        &spec,
        WindowEndpoint::Initial,
        window_balance(0.0, 0.0, 0.0, 0.0, 0.0),
    );
    let final_snapshot = window_snapshot(
        &spec,
        WindowEndpoint::Final,
        window_balance(1.0, 0.0, 0.0, 0.0, 0.0),
    );
    let inlet = window_contribution(&spec, "inlet", [1.0, 0.0, 0.0, 0.0], [0.0; 4]);
    let missing =
        WindowAuditReport::audit(spec.clone(), initial, final_snapshot, [inlet]).unwrap_err();
    assert!(matches!(
        missing,
        CoupleError::AccountingWindowContributionSetMismatch {
            missing,
            unexpected
        } if missing == vec!["contribution/source".to_string()] && unexpected.is_empty()
    ));

    let transfer = IntegratedWindowTransfer::try_new(
        Energy::new(1.0),
        Mass::new(0.0),
        [hydrogen_amount(0.0)],
        ElectricCharge::new(0.0),
        boundary_entropy("inlet", spec.interval(), [0.0; 4]),
        stable("convention/entropy-pr4"),
    )
    .unwrap();
    let wrong_role = WindowBoundaryContribution::try_new(
        &spec,
        &stable("contribution/inlet"),
        transfer,
        window_evidence(
            "wrong-integration-role",
            WindowEvidenceRole::manifest_closure(spec.manifest()),
            spec.interval(),
        ),
    )
    .unwrap_err();
    assert!(matches!(
        wrong_role,
        CoupleError::AccountingWindowEvidenceRoleMismatch {
            expected: "integrated_transfer",
            actual: "manifest_closure"
        }
    ));

    let wrong_closure = window_evidence(
        "wrong-closure-subject",
        WindowEvidenceRole::manifest_closure(&spec.manifest()[..1]),
        spec.interval(),
    );
    let closure_error = AccountingWindowSpec::try_new(
        spec.interval().clone(),
        spec.energy_reference().clone(),
        spec.element_schema().clone(),
        spec.charge_schema().clone(),
        spec.entropy_convention().clone(),
        spec.manifest().to_vec(),
        spec.tolerances(),
        None,
        wrong_closure,
    )
    .unwrap_err();
    assert!(matches!(
        closure_error,
        CoupleError::AccountingWindowContextMismatch {
            field: "evidence_subject"
        }
    ));
}

#[test]
fn accounting_context_cannot_be_rebound_across_window_specs() {
    let interval = window_interval();
    let make_charge_schema = |label: &str| {
        let basis = stable(&format!("basis/{label}"));
        WindowChargeSchema::DirectCoulomb {
            basis: basis.clone(),
            projection_evidence: window_evidence(
                &format!("charge-{label}"),
                WindowEvidenceRole::ChargeProjection {
                    subject: WindowChargeEvidenceSubject::DirectCoulomb { basis },
                },
                &interval,
            ),
        }
    };
    let spec_a = AccountingWindowSpec::try_new(
        interval.clone(),
        stable("reference/energy-a"),
        WindowElementSchema::not_applicable(stable("rationale/context-a")),
        make_charge_schema("charge-a"),
        stable("convention/entropy-a"),
        [],
        zero_window_tolerances(),
        None,
        window_closure_evidence("context-a-closure", &interval, &[]),
    )
    .unwrap();
    let spec_b = AccountingWindowSpec::try_new(
        interval.clone(),
        stable("reference/energy-b"),
        WindowElementSchema::not_applicable(stable("rationale/context-a")),
        make_charge_schema("charge-b"),
        stable("convention/entropy-b"),
        [],
        zero_window_tolerances(),
        None,
        window_closure_evidence("context-b-closure", spec_a.interval(), &[]),
    )
    .unwrap();
    let balance = WindowBalance::try_new(
        Energy::new(0.0),
        Mass::new(0.0),
        [],
        ElectricCharge::new(0.0),
        Entropy::new(0.0),
    )
    .unwrap();
    let initial = window_snapshot(&spec_a, WindowEndpoint::Initial, balance.clone());
    let final_snapshot = window_snapshot(&spec_a, WindowEndpoint::Final, balance);

    let error = WindowAuditReport::audit(spec_b, initial, final_snapshot, []).unwrap_err();
    assert!(matches!(
        error,
        CoupleError::AccountingWindowContextMismatch {
            field: "snapshot_energy_reference"
        }
    ));
}

#[test]
#[allow(clippy::too_many_lines)]
fn window_port_kinds_gate_their_admissible_balance_axes() {
    let interval = window_interval();
    let stream_entry = external_window_entry("schema-required", &interval);
    let stream_manifest = [stream_entry];
    let stream_closure =
        window_closure_evidence("missing-elements-closure", &interval, &stream_manifest);
    let missing_elements = AccountingWindowSpec::try_new(
        interval.clone(),
        window_energy_reference(),
        WindowElementSchema::not_applicable(stable("rationale/incorrectly-nonchemical")),
        direct_window_charge(&interval),
        stable("convention/entropy-pr4"),
        stream_manifest,
        zero_window_tolerances(),
        None,
        stream_closure,
    )
    .unwrap_err();
    assert!(matches!(
        missing_elements,
        CoupleError::AccountingWindowElementSchemaRequired
    ));

    let port = scalar_schema(
        PortKind::ThermalTemperatureEntropy,
        "port/window-thermal",
        "frame/window-thermal",
        PortOrientation::OutwardFromOwner,
        11,
    );
    let contribution_id = stable("contribution/thermal");
    let entry = WindowManifestEntry::external_power_port(
        contribution_id.clone(),
        stable("exchange/thermal"),
        stable("reservoir/thermal"),
        &port,
        AccountingBoundary::new(
            stable("boundary/window-thermal"),
            port.coordinates().clone(),
            BoundaryTreatment::ExternalReservoirExchange,
        ),
        window_evidence(
            "binding-thermal",
            WindowEvidenceRole::ManifestBinding {
                contribution_id: contribution_id.clone(),
                local_port_id: port.id().clone(),
            },
            &interval,
        ),
    )
    .unwrap();
    let thermal_manifest = [entry];
    let thermal_closure = window_closure_evidence("thermal-closure", &interval, &thermal_manifest);
    let spec = AccountingWindowSpec::try_new(
        interval.clone(),
        window_energy_reference(),
        WindowElementSchema::not_applicable(stable("rationale/nonchemical-thermal")),
        WindowChargeSchema::ProvenNeutral {
            projection_evidence: window_evidence(
                "thermal-neutrality",
                WindowEvidenceRole::ChargeProjection {
                    subject: WindowChargeEvidenceSubject::ProvenNeutral,
                },
                &interval,
            ),
        },
        stable("convention/entropy-pr4"),
        thermal_manifest,
        zero_window_tolerances(),
        None,
        thermal_closure,
    )
    .unwrap();
    let transfer = IntegratedWindowTransfer::try_new(
        Energy::new(0.0),
        Mass::new(1.0),
        [],
        ElectricCharge::new(0.0),
        boundary_entropy("thermal", &interval, [0.0; 4]),
        stable("convention/entropy-pr4"),
    )
    .unwrap();
    let error = WindowBoundaryContribution::try_new(
        &spec,
        &contribution_id,
        transfer,
        window_evidence(
            "integration-thermal",
            WindowEvidenceRole::IntegratedTransfer {
                contribution_id: contribution_id.clone(),
                local_port_id: port.id().clone(),
                energy_reference: spec.energy_reference().clone(),
            },
            &interval,
        ),
    )
    .unwrap_err();
    assert!(matches!(
        error,
        CoupleError::AccountingWindowPortTransferMismatch {
            field: "mass_on_power_port",
            ..
        }
    ));

    let thermal_identity_error = WindowBoundaryContribution::try_new(
        &spec,
        &contribution_id,
        IntegratedWindowTransfer::try_new(
            Energy::new(401.0),
            Mass::new(0.0),
            [],
            ElectricCharge::new(0.0),
            boundary_entropy("thermal", &interval, [0.0, 0.0, 1.0, 0.0]),
            stable("convention/entropy-pr4"),
        )
        .unwrap(),
        window_evidence(
            "integration-thermal-identity",
            WindowEvidenceRole::IntegratedTransfer {
                contribution_id: contribution_id.clone(),
                local_port_id: port.id().clone(),
                energy_reference: spec.energy_reference().clone(),
            },
            &interval,
        ),
    )
    .unwrap_err();
    assert!(matches!(
        thermal_identity_error,
        CoupleError::AccountingWindowPortTransferMismatch {
            field: "constant_temperature_energy_entropy_identity",
            ..
        }
    ));
    let constant_radiation_error = WindowBoundaryContribution::try_new(
        &spec,
        &contribution_id,
        IntegratedWindowTransfer::try_new(
            Energy::new(0.0),
            Mass::new(0.0),
            [],
            ElectricCharge::new(0.0),
            boundary_entropy("thermal", &interval, [0.0, 0.0, 0.0, 1.0]),
            stable("convention/entropy-pr4"),
        )
        .unwrap(),
        window_evidence(
            "integration-thermal-radiation",
            WindowEvidenceRole::IntegratedTransfer {
                contribution_id: contribution_id.clone(),
                local_port_id: port.id().clone(),
                energy_reference: spec.energy_reference().clone(),
            },
            &interval,
        ),
    )
    .unwrap_err();
    assert!(matches!(
        constant_radiation_error,
        CoupleError::AccountingWindowPortTransferMismatch {
            field: "constant_temperature_radiation_requires_profile",
            ..
        }
    ));

    let mechanical_port = scalar_schema(
        PortKind::MechanicalForceVelocity,
        "port/window-mechanical",
        "frame/window-mechanical",
        PortOrientation::OutwardFromOwner,
        11,
    );
    let mechanical_contribution_id = stable("contribution/mechanical");
    let mechanical_entry = WindowManifestEntry::external_power_port(
        mechanical_contribution_id.clone(),
        stable("exchange/mechanical"),
        stable("reservoir/mechanical"),
        &mechanical_port,
        AccountingBoundary::new(
            stable("boundary/window-mechanical"),
            mechanical_port.coordinates().clone(),
            BoundaryTreatment::ExternalReservoirExchange,
        ),
        window_evidence(
            "binding-mechanical",
            WindowEvidenceRole::ManifestBinding {
                contribution_id: mechanical_contribution_id.clone(),
                local_port_id: mechanical_port.id().clone(),
            },
            &interval,
        ),
    )
    .unwrap();
    let mechanical_manifest = [mechanical_entry];
    let mechanical_spec = AccountingWindowSpec::try_new(
        interval.clone(),
        window_energy_reference(),
        WindowElementSchema::not_applicable(stable("rationale/nonchemical-mechanical")),
        WindowChargeSchema::ProvenNeutral {
            projection_evidence: window_evidence(
                "mechanical-neutrality",
                WindowEvidenceRole::ChargeProjection {
                    subject: WindowChargeEvidenceSubject::ProvenNeutral,
                },
                &interval,
            ),
        },
        stable("convention/entropy-pr4"),
        mechanical_manifest.clone(),
        zero_window_tolerances(),
        None,
        window_closure_evidence("mechanical-closure", &interval, &mechanical_manifest),
    )
    .unwrap();
    let cancelling_entropy_error = WindowBoundaryContribution::try_new(
        &mechanical_spec,
        &mechanical_contribution_id,
        IntegratedWindowTransfer::try_new(
            Energy::new(0.0),
            Mass::new(0.0),
            [],
            ElectricCharge::new(0.0),
            boundary_entropy("mechanical", &interval, [1.0, -1.0, 0.0, 0.0]),
            stable("convention/entropy-pr4"),
        )
        .unwrap(),
        window_evidence(
            "integration-mechanical",
            WindowEvidenceRole::IntegratedTransfer {
                contribution_id: mechanical_contribution_id.clone(),
                local_port_id: mechanical_port.id().clone(),
                energy_reference: mechanical_spec.energy_reference().clone(),
            },
            &interval,
        ),
    )
    .unwrap_err();
    assert!(matches!(
        cancelling_entropy_error,
        CoupleError::AccountingWindowPortTransferMismatch {
            field: "entropy_on_nonthermal_power_port",
            ..
        }
    ));
}

#[test]
fn boundary_entropy_and_nonfinite_inputs_fail_closed() {
    let interval = window_interval();
    let wrong_element_subject = WindowElementSchema::try_audited(
        [ElementId::new("O").unwrap()],
        window_evidence(
            "wrong-element-subject",
            WindowEvidenceRole::ElementProjection {
                elements: vec![ElementId::new("H").unwrap()],
            },
            &interval,
        ),
    )
    .unwrap_err();
    assert!(matches!(
        wrong_element_subject,
        CoupleError::AccountingWindowContextMismatch {
            field: "evidence_subject"
        }
    ));

    let wrong_charge_subject = AccountingWindowSpec::try_new(
        interval.clone(),
        window_energy_reference(),
        WindowElementSchema::not_applicable(stable("rationale/nonchemical")),
        WindowChargeSchema::ProvenNeutral {
            projection_evidence: window_evidence(
                "wrong-charge-subject",
                WindowEvidenceRole::ChargeProjection {
                    subject: WindowChargeEvidenceSubject::DirectCoulomb {
                        basis: stable("basis/foreign-charge"),
                    },
                },
                &interval,
            ),
        },
        stable("convention/entropy-pr4"),
        [],
        zero_window_tolerances(),
        None,
        window_closure_evidence("wrong-charge-closure", &interval, &[]),
    )
    .unwrap_err();
    assert!(matches!(
        wrong_charge_subject,
        CoupleError::AccountingWindowContextMismatch {
            field: "evidence_subject"
        }
    ));

    assert!(matches!(
        BoundaryTemperatureReference::try_constant(
            Temperature::new(0.0),
            window_evidence(
                "zero-temperature",
                WindowEvidenceRole::BoundaryTemperature {
                    contribution_id: stable("contribution/temperature")
                },
                &interval,
            )
        ),
        Err(CoupleError::NonPositiveBoundaryTemperature)
    ));
    assert!(matches!(
        BoundaryEntropyBreakdown::try_new(
            Entropy::new(0.0),
            Entropy::new(0.0),
            Entropy::new(1.0),
            Entropy::new(0.0),
            BoundaryTemperatureReference::not_applicable(stable("rationale/missing-temperature"))
        ),
        Err(CoupleError::AccountingWindowTemperatureReferenceRequired)
    ));
    assert!(matches!(
        WindowBalance::try_new(
            Energy::new(f64::NAN),
            Mass::new(0.0),
            [hydrogen_amount(0.0)],
            ElectricCharge::new(0.0),
            Entropy::new(0.0)
        ),
        Err(CoupleError::NonFiniteAccountingWindowValue {
            field: "energy",
            ..
        })
    ));
}

#[test]
fn canonical_compensated_window_sum_is_permutation_stable() {
    let spec = window_spec(&["a-large", "b-cancel", "c-unit"], None);
    let initial = window_snapshot(
        &spec,
        WindowEndpoint::Initial,
        window_balance(0.0, 0.0, 0.0, 0.0, 0.0),
    );
    let final_snapshot = window_snapshot(
        &spec,
        WindowEndpoint::Final,
        window_balance(1.0, 0.0, 0.0, 0.0, 0.0),
    );
    let contributions = [
        window_contribution(&spec, "a-large", [1.0e16, 0.0, 0.0, 0.0], [0.0; 4]),
        window_contribution(&spec, "b-cancel", [-1.0e16, 0.0, 0.0, 0.0], [0.0; 4]),
        window_contribution(&spec, "c-unit", [1.0, 0.0, 0.0, 0.0], [0.0; 4]),
    ];
    let permutations = [
        [0, 1, 2],
        [0, 2, 1],
        [1, 0, 2],
        [1, 2, 0],
        [2, 0, 1],
        [2, 1, 0],
    ];
    let mut residual_bits = None;
    for permutation in permutations {
        let report = WindowAuditReport::audit(
            spec.clone(),
            initial.clone(),
            final_snapshot.clone(),
            permutation.map(|index| contributions[index].clone()),
        )
        .unwrap();
        assert!(
            report.is_green(),
            "permutation={permutation:?}, report={report:#?}"
        );
        let bits = report.residuals().first_law().value().to_bits();
        assert_eq!(bits, 0.0_f64.to_bits());
        assert_eq!(*residual_bits.get_or_insert(bits), bits);
    }
}

#[test]
fn explicitly_closed_isolated_window_audits_without_boundary_rows() {
    let interval = window_interval();
    let spec = AccountingWindowSpec::try_new(
        interval.clone(),
        window_energy_reference(),
        WindowElementSchema::not_applicable(stable("rationale/nonchemical-isolated")),
        WindowChargeSchema::ProvenNeutral {
            projection_evidence: window_evidence(
                "isolated-neutrality",
                WindowEvidenceRole::ChargeProjection {
                    subject: WindowChargeEvidenceSubject::ProvenNeutral,
                },
                &interval,
            ),
        },
        stable("convention/entropy-pr4"),
        [],
        zero_window_tolerances(),
        None,
        window_closure_evidence("isolated-closure", &interval, &[]),
    )
    .unwrap();
    let initial_balance = WindowBalance::try_new(
        Energy::new(10.0),
        Mass::new(1.0),
        [],
        ElectricCharge::new(0.0),
        Entropy::new(2.0),
    )
    .unwrap();
    let final_balance = WindowBalance::try_new(
        Energy::new(10.0),
        Mass::new(1.0),
        [],
        ElectricCharge::new(0.0),
        Entropy::new(2.5),
    )
    .unwrap();
    let initial = window_snapshot(&spec, WindowEndpoint::Initial, initial_balance);
    let final_snapshot = window_snapshot(&spec, WindowEndpoint::Final, final_balance);

    let report = WindowAuditReport::audit(spec, initial, final_snapshot, []).unwrap();
    assert!(report.is_green(), "report={report:#?}");
    assert!(report.contributions().is_empty());
    assert!(report.exchanges().is_empty());
    assert_f64_bits_eq(report.residuals().entropy_generation().value(), 0.5);
}

#[test]
#[allow(clippy::too_many_lines)]
fn hidden_internal_entropy_source_is_localized_to_its_exact_pair() {
    let interval = window_interval();
    let port_a = window_stream("pair-a");
    let port_b = window_stream("pair-b");
    let contribution_a = stable("contribution/pair-a");
    let contribution_b = stable("contribution/pair-b");
    let pair = WindowManifestEntry::internal_stream_pair(
        stable("exchange/internal-pair"),
        contribution_a.clone(),
        &port_a,
        contribution_b.clone(),
        &port_b,
        window_evidence(
            "binding-pair-a",
            WindowEvidenceRole::ManifestBinding {
                contribution_id: contribution_a,
                local_port_id: port_a.binding().port_id().clone(),
            },
            &interval,
        ),
        window_evidence(
            "binding-pair-b",
            WindowEvidenceRole::ManifestBinding {
                contribution_id: contribution_b,
                local_port_id: port_b.binding().port_id().clone(),
            },
            &interval,
        ),
    )
    .unwrap();

    let one_sided_manifest = [pair[0].clone()];
    let one_sided_closure =
        window_closure_evidence("one-sided-closure", &interval, &one_sided_manifest);
    let one_sided = AccountingWindowSpec::try_new(
        interval.clone(),
        window_energy_reference(),
        audited_window_elements(&interval),
        direct_window_charge(&interval),
        stable("convention/entropy-pr4"),
        one_sided_manifest,
        zero_window_tolerances(),
        None,
        one_sided_closure,
    )
    .unwrap_err();
    assert!(matches!(
        one_sided,
        CoupleError::MalformedAccountingWindowExchange {
            field: "internal_cardinality",
            ..
        }
    ));

    let mut full_manifest = pair.to_vec();
    full_manifest.push(external_window_entry("pair-mask", &interval));
    let pair_closure = window_closure_evidence("pair-closure", &interval, &full_manifest);
    let spec = AccountingWindowSpec::try_new(
        interval.clone(),
        window_energy_reference(),
        audited_window_elements(&interval),
        direct_window_charge(&interval),
        stable("convention/entropy-pr4"),
        full_manifest,
        zero_window_tolerances(),
        None,
        pair_closure,
    )
    .unwrap();
    let initial = window_snapshot(
        &spec,
        WindowEndpoint::Initial,
        window_balance(0.0, 0.0, 0.0, 0.0, 0.0),
    );
    let final_snapshot = window_snapshot(
        &spec,
        WindowEndpoint::Final,
        window_balance(0.0, 0.0, 0.0, 0.0, 0.0),
    );
    let pair_a = window_contribution(&spec, "pair-a", [5.0, 0.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0]);
    let pair_b = window_contribution(
        &spec,
        "pair-b",
        [-4.0, 0.0, 0.0, 0.0],
        [-0.75, 0.0, 0.0, 0.0],
    );
    let mask = window_contribution(&spec, "pair-mask", [-1.0, 0.0, 0.0, 0.0], [0.0; 4]);
    let report =
        WindowAuditReport::audit(spec, initial, final_snapshot, [mask, pair_b, pair_a]).unwrap();
    assert!(!report.is_green());
    assert_f64_bits_eq(report.residuals().first_law().value(), 0.0);
    assert_f64_bits_eq(report.residuals().entropy_generation().value(), -0.25);
    assert!(
        report
            .violations()
            .iter()
            .all(|violation| !matches!(violation.law(), WindowBalanceLaw::FirstLaw))
    );
    let internal_exchange = report
        .exchanges()
        .iter()
        .find(|exchange| exchange.exchange_id().as_str() == "exchange/internal-pair")
        .unwrap();
    assert_f64_bits_eq(internal_exchange.energy_into().value(), 1.0);
    assert_f64_bits_eq(internal_exchange.entropy_into().value(), 0.25);
    let energy_pair_violation = report
        .violations()
        .iter()
        .find(|violation| matches!(violation.law(), WindowBalanceLaw::InternalFirstLaw))
        .unwrap();
    assert_f64_bits_eq(energy_pair_violation.residual_si(), 1.0);
    assert_eq!(energy_pair_violation.exchange_ids().len(), 1);
    let pair_violation = report
        .violations()
        .iter()
        .find(|violation| matches!(violation.law(), WindowBalanceLaw::InternalEntropyTransfer))
        .unwrap();
    assert_f64_bits_eq(pair_violation.residual_si(), 0.25);
    assert_eq!(
        pair_violation
            .contribution_ids()
            .iter()
            .map(StableId::as_str)
            .collect::<Vec<_>>(),
        ["contribution/pair-a", "contribution/pair-b"]
    );
    assert_eq!(pair_violation.port_ids().len(), 2);
    assert_eq!(
        pair_violation.exchange_ids()[0].as_str(),
        "exchange/internal-pair"
    );
}

#[test]
fn coupling_is_deterministic() {
    let a = iterate_aitken(2.0, 3.0, 0.0, 0.5, 2.0, 100, 1e-9);
    let b = iterate_aitken(2.0, 3.0, 0.0, 0.5, 2.0, 100, 1e-9);
    assert_eq!(a, b);
}
