//! G0/G3 battery for I13.2a terminal-relative physical schema identities.
//!
//! The battery covers exact incidence, relative subcomplex admission,
//! canonical replay, port/Machine binding, nominal coefficient sectors, and
//! declared conversion boundaries.  It does not claim homology computation,
//! field transfer, or coil manufacturability.

#![cfg(feature = "terminal-relative")]

use fs_couple::{
    CoordinateBinding, PortKind, PortOrientation, PortSchema, PortTimestamp, StableId,
};
use fs_feec::terminal_relative::{
    BoundaryIncidence, CellRef, CellularSubcomplex, ConductorComponent, ConductorComponentId,
    ConversionMapId, DeclaredPhysicalMap, DeclaredPhysicalMapKind, DistributedCurrent,
    FiniteCellComplex, GeometricCoil, IncidenceSign, IntegralRelativeChain,
    IntegralRelativeCochain, IntegralWindingRepresentative, MachineBindingStatus,
    OrientationMapSign, PhaseId, PhysicalObjectId, PhysicalTerminal, PhysicalTerminalId,
    PresentedMachinePortRef, RealCurrentAmplitude, RealRelativeCochain, TerminalOrientation,
    TerminalPortCoordinate, TerminalPortTrivialization, TerminalRelativeCoefficientDomain,
    TerminalRelativeError, TerminalRelativePair, TerminalRole, TrivializationId,
};
use fs_qty::{Current, Dims};

fn stable(value: &str) -> StableId {
    StableId::new(value).expect("fixture stable id")
}

fn interval_complex() -> FiniteCellComplex {
    FiniteCellComplex::try_new(
        1,
        vec![2, 1],
        vec![
            BoundaryIncidence::new(
                CellRef::new(0, 0),
                CellRef::new(1, 0),
                IncidenceSign::Negative,
            ),
            BoundaryIncidence::new(
                CellRef::new(0, 1),
                CellRef::new(1, 0),
                IncidenceSign::Positive,
            ),
        ],
    )
    .expect("oriented interval")
}

fn subcomplex(
    ambient: &FiniteCellComplex,
    id: &str,
    cells: impl IntoIterator<Item = CellRef>,
) -> CellularSubcomplex {
    CellularSubcomplex::try_new(stable(id), cells, ambient).expect("fixture subcomplex")
}

fn electrical_port(id: &str, tick: u64) -> PortSchema {
    PortKind::ElectricalVoltageCurrent
        .scalar_seed_schema(
            stable(id),
            CoordinateBinding::new(
                stable("basis/winding-terminal"),
                stable("frame/winding-terminal"),
                PortOrientation::OutwardFromOwner,
            ),
            PortTimestamp::new(stable("clock/electrical"), tick),
        )
        .expect("electrical port")
}

fn terminal(
    ambient: &FiniteCellComplex,
    ordinal: u32,
    id: &str,
    role: TerminalRole,
    orientation: TerminalOrientation,
    sign: OrientationMapSign,
    tick: u64,
) -> PhysicalTerminal {
    terminal_for(
        ambient,
        ordinal,
        id,
        "component/winding",
        "phase/a",
        role,
        orientation,
        sign,
        tick,
    )
}

#[allow(clippy::too_many_arguments)]
fn terminal_for(
    ambient: &FiniteCellComplex,
    ordinal: u32,
    id: &str,
    component: &str,
    phase: &str,
    role: TerminalRole,
    orientation: TerminalOrientation,
    sign: OrientationMapSign,
    tick: u64,
) -> PhysicalTerminal {
    let port = electrical_port(&format!("port/{id}"), tick);
    PhysicalTerminal::new(
        PhysicalTerminalId::new(format!("terminal/{id}")).expect("terminal id"),
        subcomplex(
            ambient,
            &format!("support/{id}"),
            [CellRef::new(0, ordinal)],
        ),
        ConductorComponentId::new(component).expect("component id"),
        PhaseId::new(phase).expect("phase id"),
        role,
        orientation,
        TerminalPortCoordinate::Flow,
        port.clone(),
        PresentedMachinePortRef::try_new(
            stable("org.frankensim.fs-ir.machine.graph.v1"),
            1,
            [0x42; 32],
            stable("machine-owner/stator-winding"),
            stable(&format!("port/{id}")),
            stable(&format!("machine-terminal/{id}-voltage")),
            stable(&format!("machine-terminal/{id}-current")),
        )
        .expect("presented Machine-IR port"),
        TerminalPortTrivialization::new(
            TrivializationId::new(format!("trivialization/{id}")).expect("trivialization id"),
            port.id().clone(),
            sign,
            stable("voltage-reference/dc-link-negative"),
            stable(&format!("current-reference/{id}")),
        ),
    )
    .expect("physical terminal")
}

fn pair(tick: u64, reverse_declarations: bool) -> TerminalRelativePair {
    let complex = interval_complex();
    let conductor = subcomplex(
        &complex,
        "support/conductor",
        [CellRef::new(0, 0), CellRef::new(0, 1), CellRef::new(1, 0)],
    );
    let insulation = subcomplex(&complex, "support/insulation-empty", []);
    let relative = subcomplex(
        &complex,
        "support/terminal-relative",
        [CellRef::new(0, 0), CellRef::new(0, 1)],
    );
    let component = ConductorComponent::new(
        ConductorComponentId::new("component/winding").expect("component id"),
        subcomplex(
            &complex,
            "support/component-winding",
            [CellRef::new(0, 0), CellRef::new(0, 1), CellRef::new(1, 0)],
        ),
    )
    .expect("component");
    let driven = terminal(
        &complex,
        0,
        "a-positive",
        TerminalRole::Driven,
        TerminalOrientation::OutOfConductor,
        OrientationMapSign::Preserve,
        tick,
    );
    let return_reference = terminal(
        &complex,
        1,
        "a-return",
        TerminalRole::ReturnReference,
        TerminalOrientation::IntoConductor,
        OrientationMapSign::Reverse,
        tick,
    );
    let terminals = if reverse_declarations {
        vec![return_reference, driven]
    } else {
        vec![driven, return_reference]
    };
    TerminalRelativePair::try_new(
        complex,
        conductor,
        relative,
        insulation,
        vec![component],
        terminals,
    )
    .expect("terminal-relative pair")
}

fn terminal_cut_loop_pair() -> TerminalRelativePair {
    let complex = FiniteCellComplex::try_new(
        1,
        vec![4, 4],
        vec![
            BoundaryIncidence::new(
                CellRef::new(0, 0),
                CellRef::new(1, 0),
                IncidenceSign::Negative,
            ),
            BoundaryIncidence::new(
                CellRef::new(0, 1),
                CellRef::new(1, 0),
                IncidenceSign::Positive,
            ),
            BoundaryIncidence::new(
                CellRef::new(0, 1),
                CellRef::new(1, 1),
                IncidenceSign::Negative,
            ),
            BoundaryIncidence::new(
                CellRef::new(0, 2),
                CellRef::new(1, 1),
                IncidenceSign::Positive,
            ),
            BoundaryIncidence::new(
                CellRef::new(0, 1),
                CellRef::new(1, 2),
                IncidenceSign::Negative,
            ),
            BoundaryIncidence::new(
                CellRef::new(0, 2),
                CellRef::new(1, 2),
                IncidenceSign::Positive,
            ),
            BoundaryIncidence::new(
                CellRef::new(0, 2),
                CellRef::new(1, 3),
                IncidenceSign::Negative,
            ),
            BoundaryIncidence::new(
                CellRef::new(0, 3),
                CellRef::new(1, 3),
                IncidenceSign::Positive,
            ),
        ],
    )
    .expect("terminal-cut loop graph");
    let conductor = subcomplex(
        &complex,
        "support/conductor-loop",
        [
            CellRef::new(0, 0),
            CellRef::new(0, 1),
            CellRef::new(0, 2),
            CellRef::new(0, 3),
            CellRef::new(1, 0),
            CellRef::new(1, 1),
            CellRef::new(1, 2),
            CellRef::new(1, 3),
        ],
    );
    let component = ConductorComponent::new(
        ConductorComponentId::new("component/winding").unwrap(),
        conductor.clone(),
    )
    .unwrap();
    TerminalRelativePair::try_new(
        complex.clone(),
        conductor,
        subcomplex(
            &complex,
            "support/terminal-relative-loop",
            [CellRef::new(0, 0), CellRef::new(0, 3)],
        ),
        subcomplex(&complex, "support/insulation-empty-loop", []),
        vec![component],
        vec![
            terminal(
                &complex,
                0,
                "loop-positive",
                TerminalRole::Driven,
                TerminalOrientation::OutOfConductor,
                OrientationMapSign::Preserve,
                31,
            ),
            terminal(
                &complex,
                3,
                "loop-return",
                TerminalRole::ReturnReference,
                TerminalOrientation::IntoConductor,
                OrientationMapSign::Reverse,
                31,
            ),
        ],
    )
    .expect("terminal-cut loop pair")
}

#[test]
fn i13_2a_001_exact_incidence_accepts_a_triangle_and_rejects_d_squared_defect() {
    let vertices_and_edges = vec![
        BoundaryIncidence::new(
            CellRef::new(0, 0),
            CellRef::new(1, 0),
            IncidenceSign::Negative,
        ),
        BoundaryIncidence::new(
            CellRef::new(0, 1),
            CellRef::new(1, 0),
            IncidenceSign::Positive,
        ),
        BoundaryIncidence::new(
            CellRef::new(0, 1),
            CellRef::new(1, 1),
            IncidenceSign::Negative,
        ),
        BoundaryIncidence::new(
            CellRef::new(0, 2),
            CellRef::new(1, 1),
            IncidenceSign::Positive,
        ),
        BoundaryIncidence::new(
            CellRef::new(0, 0),
            CellRef::new(1, 2),
            IncidenceSign::Negative,
        ),
        BoundaryIncidence::new(
            CellRef::new(0, 2),
            CellRef::new(1, 2),
            IncidenceSign::Positive,
        ),
    ];
    let mut valid = vertices_and_edges.clone();
    valid.extend([
        BoundaryIncidence::new(
            CellRef::new(1, 0),
            CellRef::new(2, 0),
            IncidenceSign::Positive,
        ),
        BoundaryIncidence::new(
            CellRef::new(1, 1),
            CellRef::new(2, 0),
            IncidenceSign::Positive,
        ),
        BoundaryIncidence::new(
            CellRef::new(1, 2),
            CellRef::new(2, 0),
            IncidenceSign::Negative,
        ),
    ]);
    FiniteCellComplex::try_new(2, vec![3, 3, 1], valid).expect("triangle d squared is zero");

    let mut corrupt = vertices_and_edges;
    corrupt.extend([
        BoundaryIncidence::new(
            CellRef::new(1, 0),
            CellRef::new(2, 0),
            IncidenceSign::Positive,
        ),
        BoundaryIncidence::new(
            CellRef::new(1, 1),
            CellRef::new(2, 0),
            IncidenceSign::Positive,
        ),
        BoundaryIncidence::new(
            CellRef::new(1, 2),
            CellRef::new(2, 0),
            IncidenceSign::Positive,
        ),
    ]);
    assert!(matches!(
        FiniteCellComplex::try_new(2, vec![3, 3, 1], corrupt),
        Err(TerminalRelativeError::BoundarySquaredNonzero { .. })
    ));
}

#[test]
fn i13_2a_002_subcomplex_and_terminal_insulation_defects_fail_closed() {
    let complex = interval_complex();
    assert!(matches!(
        CellularSubcomplex::try_new(stable("support/not-closed"), [CellRef::new(1, 0)], &complex,),
        Err(TerminalRelativeError::NotASubcomplex { .. })
    ));

    let conductor = subcomplex(
        &complex,
        "support/conductor",
        [CellRef::new(0, 0), CellRef::new(0, 1), CellRef::new(1, 0)],
    );
    let insulation = subcomplex(&complex, "support/insulation-left", [CellRef::new(0, 0)]);
    let component = ConductorComponent::new(
        ConductorComponentId::new("component/winding").unwrap(),
        conductor.clone(),
    )
    .unwrap();
    let terminals = vec![
        terminal(
            &complex,
            0,
            "a-positive",
            TerminalRole::Driven,
            TerminalOrientation::OutOfConductor,
            OrientationMapSign::Preserve,
            7,
        ),
        terminal(
            &complex,
            1,
            "a-return",
            TerminalRole::ReturnReference,
            TerminalOrientation::IntoConductor,
            OrientationMapSign::Reverse,
            7,
        ),
    ];
    assert!(matches!(
        TerminalRelativePair::try_new(
            complex,
            conductor,
            subcomplex(
                &interval_complex(),
                "support/terminal-relative-overlap",
                [CellRef::new(0, 0), CellRef::new(0, 1)],
            ),
            insulation,
            vec![component],
            terminals,
        ),
        Err(TerminalRelativeError::TerminalInsulationOverlap { .. })
    ));
}

#[test]
fn i13_2a_003_phase_reference_and_orientation_semantics_are_mandatory() {
    let complex = interval_complex();
    let conductor = subcomplex(
        &complex,
        "support/conductor",
        [CellRef::new(0, 0), CellRef::new(0, 1), CellRef::new(1, 0)],
    );
    let component = ConductorComponent::new(
        ConductorComponentId::new("component/winding").unwrap(),
        conductor.clone(),
    )
    .unwrap();
    let both_driven = vec![
        terminal(
            &complex,
            0,
            "a-positive",
            TerminalRole::Driven,
            TerminalOrientation::OutOfConductor,
            OrientationMapSign::Preserve,
            7,
        ),
        terminal(
            &complex,
            1,
            "a-return",
            TerminalRole::Driven,
            TerminalOrientation::IntoConductor,
            OrientationMapSign::Reverse,
            7,
        ),
    ];
    assert!(matches!(
        TerminalRelativePair::try_new(
            complex.clone(),
            conductor.clone(),
            subcomplex(
                &complex,
                "support/terminal-relative-missing-role",
                [CellRef::new(0, 0), CellRef::new(0, 1)],
            ),
            subcomplex(&complex, "support/insulation-empty", []),
            vec![component.clone()],
            both_driven,
        ),
        Err(TerminalRelativeError::MissingPhaseRole {
            role: TerminalRole::ReturnReference,
            ..
        })
    ));

    let same_orientation = vec![
        terminal(
            &complex,
            0,
            "a-positive",
            TerminalRole::Driven,
            TerminalOrientation::OutOfConductor,
            OrientationMapSign::Preserve,
            7,
        ),
        terminal(
            &complex,
            1,
            "a-return",
            TerminalRole::ReturnReference,
            TerminalOrientation::OutOfConductor,
            OrientationMapSign::Preserve,
            7,
        ),
    ];
    assert!(matches!(
        TerminalRelativePair::try_new(
            complex.clone(),
            conductor,
            subcomplex(
                &complex,
                "support/terminal-relative-orientation",
                [CellRef::new(0, 0), CellRef::new(0, 1)],
            ),
            subcomplex(&complex, "support/insulation-empty-2", []),
            vec![component],
            same_orientation,
        ),
        Err(TerminalRelativeError::PhaseOrientationDoesNotClose { .. })
    ));
}

#[test]
fn i13_2a_004_port_and_machine_bindings_are_presented_and_identity_bearing() {
    let canonical = pair(17, false);
    let permuted = pair(17, true);
    let retimed = pair(18, false);
    assert_eq!(canonical.identity(), permuted.identity());
    assert_ne!(canonical.identity(), retimed.identity());
    assert!(canonical.canonical_bytes() > 0);
    let receipt = canonical.complex_receipt();
    assert_eq!(receipt.identity_receipt().id(), canonical.identity());
    assert_eq!(
        receipt.coefficient_domains(),
        &[
            TerminalRelativeCoefficientDomain::Integers,
            TerminalRelativeCoefficientDomain::FiniteReal,
        ]
    );
    assert_eq!(receipt.current_dimensions(), Current::DIMS);
    assert_eq!(receipt.terminal_bindings().len(), 2);
    assert_eq!(
        receipt.machine_binding(),
        MachineBindingStatus::PresentedOnly
    );

    let positive = canonical
        .terminals()
        .iter()
        .find(|terminal| terminal.role() == TerminalRole::Driven)
        .expect("driven terminal");
    let positive_receipt = receipt
        .terminal_bindings()
        .iter()
        .find(|binding| binding.terminal() == positive.id())
        .expect("driven terminal receipt");
    assert_eq!(positive_receipt.port_schema(), positive.port());
    assert_eq!(positive_receipt.machine(), positive.machine());
    assert_eq!(positive_receipt.trivialization(), positive.trivialization());
    assert_eq!(positive.phase().as_str(), "phase/a");
    assert_eq!(positive.port().kind(), PortKind::ElectricalVoltageCurrent);
    assert_eq!(positive.coordinate(), TerminalPortCoordinate::Flow);
    assert_eq!(positive.port().timestamp().tick(), 17);
    assert_eq!(
        positive.port().coordinates().orientation(),
        PortOrientation::OutwardFromOwner
    );
    assert_eq!(
        positive.machine().authority_domain().as_str(),
        "org.frankensim.fs-ir.machine.graph.v1"
    );
    assert_eq!(
        positive.machine().flow_terminal().as_str(),
        "machine-terminal/a-positive-current"
    );
    assert_eq!(
        positive.trivialization().sign(),
        OrientationMapSign::Preserve
    );
}

#[test]
fn i13_2a_005_integral_and_real_objects_remain_nominally_distinct() {
    let pair = pair(23, false);
    let phase = PhaseId::new("phase/a").unwrap();
    let chain = IntegralRelativeChain::try_new(&pair, phase.clone(), 1, vec![3]).unwrap();
    let boundary = pair.boundary(&chain).expect("relative boundary");
    assert_eq!(boundary.degree(), 0);
    assert!(boundary.coefficients().is_empty());

    let representative =
        IntegralWindingRepresentative::try_new(&pair, phase.clone(), vec![3]).unwrap();
    assert_eq!(representative.chain().coefficients(), &[3]);
    let scaled_representative =
        IntegralWindingRepresentative::try_new(&pair, phase.clone(), vec![4]).unwrap();
    assert_ne!(representative.identity(), scaled_representative.identity());
    assert_ne!(
        representative.object_ref().identity(),
        scaled_representative.object_ref().identity()
    );

    let amplitude = RealCurrentAmplitude::try_new(
        PhysicalObjectId::new("object/current-amplitude").unwrap(),
        &pair,
        phase.clone(),
        Current::new(2.5),
    )
    .unwrap();
    assert_eq!(amplitude.value().value().to_bits(), 2.5_f64.to_bits());

    let current_cochain =
        RealRelativeCochain::try_new(&pair, phase.clone(), 1, Current::DIMS, vec![2.5]).unwrap();
    let distributed = DistributedCurrent::new(
        PhysicalObjectId::new("object/distributed-current").unwrap(),
        current_cochain,
        stable("receipt/divergence-v1"),
        stable("receipt/terminal-constraint-v1"),
    )
    .unwrap();
    assert_eq!(distributed.cochain().values(), &[2.5]);

    let coil = GeometricCoil::try_new(
        PhysicalObjectId::new("object/geometric-coil").unwrap(),
        &pair,
        phase,
        ConductorComponentId::new("component/winding").unwrap(),
        stable("artifact/connectivity-v1"),
        stable("artifact/manufacturing-v1"),
    )
    .unwrap();

    let winding_realization = DeclaredPhysicalMap::try_new(
        ConversionMapId::new("map/winding-realization").unwrap(),
        DeclaredPhysicalMapKind::WindingRealization,
        representative.object_ref(),
        coil.object_ref(),
        stable("artifact/winding-realization-v1"),
    )
    .unwrap();
    assert_ne!(
        winding_realization.source().kind(),
        winding_realization.target().kind()
    );

    let current_realization = DeclaredPhysicalMap::try_new(
        ConversionMapId::new("map/current-realization").unwrap(),
        DeclaredPhysicalMapKind::CurrentRealization,
        amplitude.object_ref(),
        distributed.object_ref(),
        stable("artifact/current-realization-v1"),
    )
    .unwrap();
    assert_ne!(
        current_realization.source().kind(),
        current_realization.target().kind()
    );
}

#[test]
fn i13_2a_006_real_coboundary_is_typed_and_nonfinite_values_refuse() {
    let pair = pair(29, false);
    let phase = PhaseId::new("phase/a").unwrap();
    let zero_form =
        RealRelativeCochain::try_new(&pair, phase.clone(), 0, Dims::NONE, Vec::new()).unwrap();
    let one_form = pair.coboundary(&zero_form).unwrap();
    assert_eq!(one_form.degree(), 1);
    assert_eq!(one_form.values(), &[0.0]);

    assert_eq!(
        RealRelativeCochain::try_new(&pair, phase, 1, Current::DIMS, vec![f64::NAN]),
        Err(TerminalRelativeError::NonFiniteRealCoefficient { index: 0 })
    );
}

#[test]
fn i13_2a_007_duplicate_incidence_and_wrong_trivialization_refuse() {
    let duplicate = BoundaryIncidence::new(
        CellRef::new(0, 0),
        CellRef::new(1, 0),
        IncidenceSign::Negative,
    );
    assert!(matches!(
        FiniteCellComplex::try_new(1, vec![2, 1], vec![duplicate, duplicate]),
        Err(TerminalRelativeError::DuplicateIncidence { .. })
    ));

    let complex = interval_complex();
    let port = electrical_port("port/mismatch", 1);
    assert!(matches!(
        PhysicalTerminal::new(
            PhysicalTerminalId::new("terminal/mismatch").unwrap(),
            subcomplex(&complex, "support/mismatch", [CellRef::new(0, 0)]),
            ConductorComponentId::new("component/winding").unwrap(),
            PhaseId::new("phase/a").unwrap(),
            TerminalRole::Driven,
            TerminalOrientation::OutOfConductor,
            TerminalPortCoordinate::Flow,
            port,
            PresentedMachinePortRef::try_new(
                stable("org.frankensim.fs-ir.machine.graph.v1"),
                1,
                [0x24; 32],
                stable("machine-owner/stator-winding"),
                stable("port/mismatch"),
                stable("machine-terminal/mismatch-voltage"),
                stable("machine-terminal/mismatch-current"),
            )
            .unwrap(),
            TerminalPortTrivialization::new(
                TrivializationId::new("trivialization/mismatch").unwrap(),
                stable("port/foreign"),
                OrientationMapSign::Preserve,
                stable("voltage-reference/zero"),
                stable("current-reference/positive"),
            ),
        ),
        Err(TerminalRelativeError::TrivializationPortMismatch { .. })
    ));

    let port = electrical_port("port/effort-selected", 1);
    let port_id = port.id().clone();
    assert!(matches!(
        PhysicalTerminal::new(
            PhysicalTerminalId::new("terminal/effort-selected").unwrap(),
            subcomplex(&complex, "support/effort-selected", [CellRef::new(0, 0)],),
            ConductorComponentId::new("component/winding").unwrap(),
            PhaseId::new("phase/a").unwrap(),
            TerminalRole::Driven,
            TerminalOrientation::OutOfConductor,
            TerminalPortCoordinate::Effort,
            port,
            PresentedMachinePortRef::try_new(
                stable("org.frankensim.fs-ir.machine.graph.v1"),
                1,
                [0x25; 32],
                stable("machine-owner/stator-winding"),
                stable("port/effort-selected"),
                stable("machine-terminal/effort-selected-voltage"),
                stable("machine-terminal/effort-selected-current"),
            )
            .unwrap(),
            TerminalPortTrivialization::new(
                TrivializationId::new("trivialization/effort-selected").unwrap(),
                port_id,
                OrientationMapSign::Preserve,
                stable("voltage-reference/zero"),
                stable("current-reference/positive"),
            ),
        ),
        Err(TerminalRelativeError::TerminalRequiresFlowCoordinate { .. })
    ));

    let port = electrical_port("port/orientation-mismatch", 1);
    let port_id = port.id().clone();
    assert!(matches!(
        PhysicalTerminal::new(
            PhysicalTerminalId::new("terminal/orientation-mismatch").unwrap(),
            subcomplex(
                &complex,
                "support/orientation-mismatch",
                [CellRef::new(0, 0)],
            ),
            ConductorComponentId::new("component/winding").unwrap(),
            PhaseId::new("phase/a").unwrap(),
            TerminalRole::Driven,
            TerminalOrientation::IntoConductor,
            TerminalPortCoordinate::Flow,
            port,
            PresentedMachinePortRef::try_new(
                stable("org.frankensim.fs-ir.machine.graph.v1"),
                1,
                [0x27; 32],
                stable("machine-owner/stator-winding"),
                stable("port/orientation-mismatch"),
                stable("machine-terminal/orientation-mismatch-voltage"),
                stable("machine-terminal/orientation-mismatch-current"),
            )
            .unwrap(),
            TerminalPortTrivialization::new(
                TrivializationId::new("trivialization/orientation-mismatch").unwrap(),
                port_id,
                OrientationMapSign::Preserve,
                stable("voltage-reference/zero"),
                stable("current-reference/positive"),
            ),
        ),
        Err(TerminalRelativeError::TerminalOrientationMismatch { .. })
    ));

    assert!(matches!(
        PresentedMachinePortRef::try_new(
            stable("org.frankensim.fs-ir.machine-graph.v1"),
            1,
            [0x26; 32],
            stable("machine-owner/stator-winding"),
            stable("port/wrong-domain"),
            stable("machine-terminal/wrong-domain-voltage"),
            stable("machine-terminal/wrong-domain-current"),
        ),
        Err(TerminalRelativeError::MachineGraphSchemaMismatch { .. })
    ));
    assert!(matches!(
        PresentedMachinePortRef::try_new(
            stable("org.frankensim.fs-ir.machine.graph.v1"),
            1,
            [0; 32],
            stable("machine-owner/stator-winding"),
            stable("port/zero-graph"),
            stable("machine-terminal/zero-graph-voltage"),
            stable("machine-terminal/zero-graph-current"),
        ),
        Err(TerminalRelativeError::ZeroMachineGraphDigest)
    ));
}

#[test]
fn i13_2a_008_relative_subcomplex_is_explicit_and_contained() {
    let complex = interval_complex();
    let conductor = subcomplex(&complex, "support/conductor-left", [CellRef::new(0, 0)]);
    let outside = subcomplex(&complex, "support/relative-outside", [CellRef::new(0, 1)]);
    assert!(matches!(
        TerminalRelativePair::try_new(
            complex.clone(),
            conductor,
            outside,
            subcomplex(&complex, "support/insulation-empty", []),
            Vec::new(),
            Vec::new(),
        ),
        Err(TerminalRelativeError::RelativeOutsideConductor { .. })
    ));

    let conductor = subcomplex(
        &complex,
        "support/conductor",
        [CellRef::new(0, 0), CellRef::new(0, 1), CellRef::new(1, 0)],
    );
    let component = ConductorComponent::new(
        ConductorComponentId::new("component/winding").unwrap(),
        conductor.clone(),
    )
    .unwrap();
    let terminals = vec![
        terminal(
            &complex,
            0,
            "a-positive",
            TerminalRole::Driven,
            TerminalOrientation::OutOfConductor,
            OrientationMapSign::Preserve,
            7,
        ),
        terminal(
            &complex,
            1,
            "a-return",
            TerminalRole::ReturnReference,
            TerminalOrientation::IntoConductor,
            OrientationMapSign::Reverse,
            7,
        ),
    ];
    assert!(matches!(
        TerminalRelativePair::try_new(
            complex.clone(),
            conductor,
            subcomplex(&complex, "support/relative-only-left", [CellRef::new(0, 0)],),
            subcomplex(&complex, "support/insulation-empty-2", []),
            vec![component],
            terminals,
        ),
        Err(TerminalRelativeError::TerminalOutsideRelativeSubcomplex { .. })
    ));
}

#[test]
fn i13_2a_009_components_must_be_full_dimensional_closures() {
    let complex = interval_complex();
    let conductor = subcomplex(
        &complex,
        "support/conductor",
        [CellRef::new(0, 0), CellRef::new(0, 1), CellRef::new(1, 0)],
    );
    let winding = ConductorComponent::new(
        ConductorComponentId::new("component/winding").unwrap(),
        conductor.clone(),
    )
    .unwrap();
    let ghost = ConductorComponent::new(
        ConductorComponentId::new("component/ghost").unwrap(),
        subcomplex(&complex, "support/ghost", [CellRef::new(0, 0)]),
    )
    .unwrap();
    let terminals = vec![
        terminal(
            &complex,
            0,
            "a-positive",
            TerminalRole::Driven,
            TerminalOrientation::OutOfConductor,
            OrientationMapSign::Preserve,
            7,
        ),
        terminal(
            &complex,
            1,
            "a-return",
            TerminalRole::ReturnReference,
            TerminalOrientation::IntoConductor,
            OrientationMapSign::Reverse,
            7,
        ),
    ];
    assert!(matches!(
        TerminalRelativePair::try_new(
            complex.clone(),
            conductor,
            subcomplex(
                &complex,
                "support/terminal-relative",
                [CellRef::new(0, 0), CellRef::new(0, 1)],
            ),
            subcomplex(&complex, "support/insulation-empty", []),
            vec![winding, ghost],
            terminals,
        ),
        Err(TerminalRelativeError::ComponentHasNoTopCell { .. })
    ));
}

#[test]
fn i13_2a_010_integral_coboundary_satisfies_exact_stokes_pairing() {
    let pair = terminal_cut_loop_pair();
    let phase = PhaseId::new("phase/a").unwrap();
    let alpha = IntegralRelativeCochain::try_new(&pair, phase.clone(), 0, vec![2, 5])
        .expect("integral zero-cochain");
    let delta_alpha = pair
        .integral_coboundary(&alpha)
        .expect("exact integral coboundary");
    assert_eq!(delta_alpha.coefficients(), &[2, 3, 3, -5]);

    let arbitrary = IntegralRelativeChain::try_new(&pair, phase.clone(), 1, vec![2, -1, 3, 4])
        .expect("arbitrary relative chain");
    let boundary = pair.boundary(&arbitrary).expect("relative boundary");
    assert_eq!(boundary.coefficients(), &[0, -2]);
    assert_eq!(
        pair.integral_pairing(&delta_alpha, &arbitrary).unwrap(),
        pair.integral_pairing(&alpha, &boundary).unwrap()
    );

    for coefficients in [vec![1, 1, 0, 1], vec![1, 0, 1, 1], vec![0, 1, -1, 0]] {
        let cycle = IntegralRelativeChain::try_new(&pair, phase.clone(), 1, coefficients)
            .expect("relative cycle candidate");
        assert_eq!(pair.boundary(&cycle).unwrap().coefficients(), &[0, 0]);
        assert_eq!(pair.integral_pairing(&delta_alpha, &cycle).unwrap(), 0);
    }
}

fn disconnected_two_phase_pair() -> TerminalRelativePair {
    let complex = FiniteCellComplex::try_new(
        1,
        vec![4, 2],
        vec![
            BoundaryIncidence::new(
                CellRef::new(0, 0),
                CellRef::new(1, 0),
                IncidenceSign::Negative,
            ),
            BoundaryIncidence::new(
                CellRef::new(0, 1),
                CellRef::new(1, 0),
                IncidenceSign::Positive,
            ),
            BoundaryIncidence::new(
                CellRef::new(0, 2),
                CellRef::new(1, 1),
                IncidenceSign::Negative,
            ),
            BoundaryIncidence::new(
                CellRef::new(0, 3),
                CellRef::new(1, 1),
                IncidenceSign::Positive,
            ),
        ],
    )
    .unwrap();
    let conductor = subcomplex(
        &complex,
        "support/two-phase-conductor",
        [
            CellRef::new(0, 0),
            CellRef::new(0, 1),
            CellRef::new(0, 2),
            CellRef::new(0, 3),
            CellRef::new(1, 0),
            CellRef::new(1, 1),
        ],
    );
    let component_a = ConductorComponent::new(
        ConductorComponentId::new("component/a").unwrap(),
        subcomplex(
            &complex,
            "support/component-a",
            [CellRef::new(0, 0), CellRef::new(0, 1), CellRef::new(1, 0)],
        ),
    )
    .unwrap();
    let component_b = ConductorComponent::new(
        ConductorComponentId::new("component/b").unwrap(),
        subcomplex(
            &complex,
            "support/component-b",
            [CellRef::new(0, 2), CellRef::new(0, 3), CellRef::new(1, 1)],
        ),
    )
    .unwrap();
    TerminalRelativePair::try_new(
        complex.clone(),
        conductor,
        subcomplex(
            &complex,
            "support/two-phase-relative",
            [
                CellRef::new(0, 0),
                CellRef::new(0, 1),
                CellRef::new(0, 2),
                CellRef::new(0, 3),
            ],
        ),
        subcomplex(&complex, "support/two-phase-insulation-empty", []),
        vec![component_b, component_a],
        vec![
            terminal_for(
                &complex,
                0,
                "a-positive",
                "component/a",
                "phase/a",
                TerminalRole::Driven,
                TerminalOrientation::OutOfConductor,
                OrientationMapSign::Preserve,
                37,
            ),
            terminal_for(
                &complex,
                1,
                "a-return",
                "component/a",
                "phase/a",
                TerminalRole::ReturnReference,
                TerminalOrientation::IntoConductor,
                OrientationMapSign::Reverse,
                37,
            ),
            terminal_for(
                &complex,
                2,
                "b-positive",
                "component/b",
                "phase/b",
                TerminalRole::Driven,
                TerminalOrientation::OutOfConductor,
                OrientationMapSign::Preserve,
                37,
            ),
            terminal_for(
                &complex,
                3,
                "b-return",
                "component/b",
                "phase/b",
                TerminalRole::ReturnReference,
                TerminalOrientation::IntoConductor,
                OrientationMapSign::Reverse,
                37,
            ),
        ],
    )
    .expect("disconnected two-phase pair")
}

#[test]
fn i13_2a_011_phase_bases_restrict_top_cells_to_owned_components() {
    let pair = disconnected_two_phase_pair();
    let phase_a = PhaseId::new("phase/a").unwrap();
    let phase_b = PhaseId::new("phase/b").unwrap();
    assert_eq!(
        pair.phase_relative_basis(&phase_a, 1),
        Ok(vec![CellRef::new(1, 0)])
    );
    assert_eq!(
        pair.phase_relative_basis(&phase_b, 1),
        Ok(vec![CellRef::new(1, 1)])
    );
    assert!(matches!(
        IntegralRelativeChain::try_new(&pair, phase_a, 1, vec![1, 0]),
        Err(TerminalRelativeError::CoefficientArity {
            expected: 1,
            actual: 2
        })
    ));
}

#[test]
fn i13_2a_012_phase_component_admission_refuses_ambiguous_bindings() {
    let admitted = disconnected_two_phase_pair();
    let phase_a_terminals = admitted
        .terminals()
        .iter()
        .filter(|terminal| terminal.phase().as_str() == "phase/a")
        .cloned()
        .collect();
    assert_eq!(
        TerminalRelativePair::try_new(
            admitted.complex().clone(),
            admitted.conductor().clone(),
            admitted.relative().clone(),
            admitted.insulation().clone(),
            admitted.components().to_vec(),
            phase_a_terminals,
        ),
        Err(TerminalRelativeError::UnboundConductorComponent {
            component: "component/b".to_owned(),
        })
    );

    assert_eq!(
        TerminalRelativePair::try_new(
            admitted.complex().clone(),
            admitted.conductor().clone(),
            admitted.relative().clone(),
            admitted.insulation().clone(),
            admitted.components().to_vec(),
            vec![
                terminal_for(
                    admitted.complex(),
                    0,
                    "mixed-positive",
                    "component/a",
                    "phase/a",
                    TerminalRole::Driven,
                    TerminalOrientation::OutOfConductor,
                    OrientationMapSign::Preserve,
                    41,
                ),
                terminal_for(
                    admitted.complex(),
                    3,
                    "mixed-return",
                    "component/b",
                    "phase/a",
                    TerminalRole::ReturnReference,
                    TerminalOrientation::IntoConductor,
                    OrientationMapSign::Reverse,
                    41,
                ),
            ],
        ),
        Err(TerminalRelativeError::PhaseComponentMismatch {
            phase: "phase/a".to_owned(),
            driven_component: "component/a".to_owned(),
            return_component: "component/b".to_owned(),
        })
    );

    let shared_component = ConductorComponent::new(
        ConductorComponentId::new("component/shared").unwrap(),
        admitted.conductor().clone(),
    )
    .unwrap();
    assert_eq!(
        TerminalRelativePair::try_new(
            admitted.complex().clone(),
            admitted.conductor().clone(),
            admitted.relative().clone(),
            admitted.insulation().clone(),
            vec![shared_component],
            vec![
                terminal_for(
                    admitted.complex(),
                    0,
                    "shared-a-positive",
                    "component/shared",
                    "phase/a",
                    TerminalRole::Driven,
                    TerminalOrientation::OutOfConductor,
                    OrientationMapSign::Preserve,
                    43,
                ),
                terminal_for(
                    admitted.complex(),
                    1,
                    "shared-a-return",
                    "component/shared",
                    "phase/a",
                    TerminalRole::ReturnReference,
                    TerminalOrientation::IntoConductor,
                    OrientationMapSign::Reverse,
                    43,
                ),
                terminal_for(
                    admitted.complex(),
                    2,
                    "shared-b-positive",
                    "component/shared",
                    "phase/b",
                    TerminalRole::Driven,
                    TerminalOrientation::OutOfConductor,
                    OrientationMapSign::Preserve,
                    43,
                ),
                terminal_for(
                    admitted.complex(),
                    3,
                    "shared-b-return",
                    "component/shared",
                    "phase/b",
                    TerminalRole::ReturnReference,
                    TerminalOrientation::IntoConductor,
                    OrientationMapSign::Reverse,
                    43,
                ),
            ],
        ),
        Err(TerminalRelativeError::ComponentPhaseConflict {
            component: "component/shared".to_owned(),
            first_phase: "phase/a".to_owned(),
            second_phase: "phase/b".to_owned(),
        })
    );
}

#[test]
fn i13_2a_013_pairing_and_geometry_preserve_phase_component_bindings() {
    let pair = disconnected_two_phase_pair();
    let phase_a = PhaseId::new("phase/a").unwrap();
    let phase_b = PhaseId::new("phase/b").unwrap();
    let alpha_a = IntegralRelativeCochain::try_new(&pair, phase_a.clone(), 0, Vec::new()).unwrap();
    let zero_b = IntegralRelativeChain::try_new(&pair, phase_b, 0, Vec::new()).unwrap();
    assert_eq!(
        pair.integral_pairing(&alpha_a, &zero_b),
        Err(TerminalRelativeError::PairingPhaseMismatch {
            cochain: "phase/a".to_owned(),
            chain: "phase/b".to_owned(),
        })
    );

    let edge_a = IntegralRelativeChain::try_new(&pair, phase_a.clone(), 1, vec![1]).unwrap();
    assert_eq!(
        pair.integral_pairing(&alpha_a, &edge_a),
        Err(TerminalRelativeError::PairingDegreeMismatch {
            cochain: 0,
            chain: 1,
        })
    );

    assert_eq!(
        GeometricCoil::try_new(
            PhysicalObjectId::new("object/cross-phase-coil").unwrap(),
            &pair,
            phase_a,
            ConductorComponentId::new("component/b").unwrap(),
            stable("artifact/cross-phase-connectivity"),
            stable("artifact/cross-phase-manufacturing"),
        ),
        Err(TerminalRelativeError::CoilPhaseComponentMismatch {
            phase: "phase/a".to_owned(),
            expected_component: "component/a".to_owned(),
            actual_component: "component/b".to_owned(),
        })
    );
}

#[test]
fn i13_2a_014_shared_lower_closure_cells_remain_phase_tagged() {
    let complex = FiniteCellComplex::try_new(
        1,
        vec![5, 4],
        vec![
            BoundaryIncidence::new(
                CellRef::new(0, 0),
                CellRef::new(1, 0),
                IncidenceSign::Negative,
            ),
            BoundaryIncidence::new(
                CellRef::new(0, 4),
                CellRef::new(1, 0),
                IncidenceSign::Positive,
            ),
            BoundaryIncidence::new(
                CellRef::new(0, 1),
                CellRef::new(1, 1),
                IncidenceSign::Negative,
            ),
            BoundaryIncidence::new(
                CellRef::new(0, 4),
                CellRef::new(1, 1),
                IncidenceSign::Positive,
            ),
            BoundaryIncidence::new(
                CellRef::new(0, 2),
                CellRef::new(1, 2),
                IncidenceSign::Negative,
            ),
            BoundaryIncidence::new(
                CellRef::new(0, 4),
                CellRef::new(1, 2),
                IncidenceSign::Positive,
            ),
            BoundaryIncidence::new(
                CellRef::new(0, 3),
                CellRef::new(1, 3),
                IncidenceSign::Negative,
            ),
            BoundaryIncidence::new(
                CellRef::new(0, 4),
                CellRef::new(1, 3),
                IncidenceSign::Positive,
            ),
        ],
    )
    .unwrap();
    let conductor = subcomplex(
        &complex,
        "support/shared-closure-conductor",
        [
            CellRef::new(0, 0),
            CellRef::new(0, 1),
            CellRef::new(0, 2),
            CellRef::new(0, 3),
            CellRef::new(0, 4),
            CellRef::new(1, 0),
            CellRef::new(1, 1),
            CellRef::new(1, 2),
            CellRef::new(1, 3),
        ],
    );
    let component_a = ConductorComponent::new(
        ConductorComponentId::new("component/a").unwrap(),
        subcomplex(
            &complex,
            "support/shared-closure-a",
            [
                CellRef::new(0, 0),
                CellRef::new(0, 1),
                CellRef::new(0, 4),
                CellRef::new(1, 0),
                CellRef::new(1, 1),
            ],
        ),
    )
    .unwrap();
    let component_b = ConductorComponent::new(
        ConductorComponentId::new("component/b").unwrap(),
        subcomplex(
            &complex,
            "support/shared-closure-b",
            [
                CellRef::new(0, 2),
                CellRef::new(0, 3),
                CellRef::new(0, 4),
                CellRef::new(1, 2),
                CellRef::new(1, 3),
            ],
        ),
    )
    .unwrap();
    let pair = TerminalRelativePair::try_new(
        complex.clone(),
        conductor,
        subcomplex(
            &complex,
            "support/shared-closure-relative",
            [
                CellRef::new(0, 0),
                CellRef::new(0, 1),
                CellRef::new(0, 2),
                CellRef::new(0, 3),
            ],
        ),
        subcomplex(&complex, "support/shared-closure-insulation-empty", []),
        vec![component_a, component_b],
        vec![
            terminal_for(
                &complex,
                0,
                "shared-closure-a-positive",
                "component/a",
                "phase/a",
                TerminalRole::Driven,
                TerminalOrientation::OutOfConductor,
                OrientationMapSign::Preserve,
                47,
            ),
            terminal_for(
                &complex,
                1,
                "shared-closure-a-return",
                "component/a",
                "phase/a",
                TerminalRole::ReturnReference,
                TerminalOrientation::IntoConductor,
                OrientationMapSign::Reverse,
                47,
            ),
            terminal_for(
                &complex,
                2,
                "shared-closure-b-positive",
                "component/b",
                "phase/b",
                TerminalRole::Driven,
                TerminalOrientation::OutOfConductor,
                OrientationMapSign::Preserve,
                47,
            ),
            terminal_for(
                &complex,
                3,
                "shared-closure-b-return",
                "component/b",
                "phase/b",
                TerminalRole::ReturnReference,
                TerminalOrientation::IntoConductor,
                OrientationMapSign::Reverse,
                47,
            ),
        ],
    )
    .unwrap();

    let phase_a = PhaseId::new("phase/a").unwrap();
    let phase_b = PhaseId::new("phase/b").unwrap();
    let shared_vertex = vec![CellRef::new(0, 4)];
    assert_eq!(
        pair.phase_relative_basis(&phase_a, 0),
        Ok(shared_vertex.clone())
    );
    assert_eq!(pair.phase_relative_basis(&phase_b, 0), Ok(shared_vertex));
    let chain_a = IntegralRelativeChain::try_new(&pair, phase_a, 0, vec![1]).unwrap();
    let chain_b = IntegralRelativeChain::try_new(&pair, phase_b, 0, vec![1]).unwrap();
    assert_ne!(chain_a, chain_b);
}
