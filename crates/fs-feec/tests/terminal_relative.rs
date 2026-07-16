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
    BoundaryIncidence, CellRef, CellularSubcomplex, ComponentRelabelEntry, ConductorComponent,
    ConductorComponentId, ConversionMapId, DeclaredPhysicalMap, DeclaredPhysicalMapKind,
    DistributedCurrent, FiniteCellComplex, GeometricCoil, IncidenceSign, IntegralRelativeChain,
    IntegralRelativeCochain, IntegralWindingRepresentative, MachineBindingStatus,
    OrientationMapSign, PhaseCurrentSign, PhaseId, PhaseRelabelEntry, PhysicalObjectId,
    PhysicalTerminal, PhysicalTerminalId, PresentedMachinePortRef, RealCurrentAmplitude,
    RealRelativeCochain, SignedCellRelabelEntry, TerminalOrientation, TerminalPortCoordinate,
    TerminalPortTrivialization, TerminalRelabelEntry, TerminalRelativeCoefficientDomain,
    TerminalRelativeError, TerminalRelativePair, TerminalRelativePhysicalRelabel,
    TerminalRelativePhysicalRelabelError, TerminalRelativeSemanticPermutation,
    TerminalRelativeSignedRelabel, TerminalRelativeSignedRelabelError, TerminalRole,
    TrivializationId,
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
    terminal_for_with_current_reference(
        ambient,
        ordinal,
        id,
        component,
        phase,
        role,
        orientation,
        sign,
        tick,
        stable(&format!("current-reference/{id}")),
    )
}

#[allow(clippy::too_many_arguments)]
fn terminal_for_with_current_reference(
    ambient: &FiniteCellComplex,
    ordinal: u32,
    id: &str,
    component: &str,
    phase: &str,
    role: TerminalRole,
    orientation: TerminalOrientation,
    sign: OrientationMapSign,
    tick: u64,
    current_reference: StableId,
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
            current_reference,
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
    terminal_cut_loop_pair_with_terminals(0, 3)
}

fn terminal_cut_loop_pair_with_terminals(
    driven_ordinal: u32,
    return_ordinal: u32,
) -> TerminalRelativePair {
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
                driven_ordinal,
                "loop-positive",
                TerminalRole::Driven,
                TerminalOrientation::OutOfConductor,
                OrientationMapSign::Preserve,
                31,
            ),
            terminal(
                &complex,
                return_ordinal,
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

fn parallel_edge_relabel_entries() -> Vec<SignedCellRelabelEntry> {
    vec![
        SignedCellRelabelEntry::new(
            CellRef::new(0, 0),
            CellRef::new(0, 0),
            IncidenceSign::Positive,
        ),
        SignedCellRelabelEntry::new(
            CellRef::new(0, 1),
            CellRef::new(0, 1),
            IncidenceSign::Positive,
        ),
        SignedCellRelabelEntry::new(
            CellRef::new(0, 2),
            CellRef::new(0, 2),
            IncidenceSign::Positive,
        ),
        SignedCellRelabelEntry::new(
            CellRef::new(0, 3),
            CellRef::new(0, 3),
            IncidenceSign::Positive,
        ),
        SignedCellRelabelEntry::new(
            CellRef::new(1, 0),
            CellRef::new(1, 0),
            IncidenceSign::Positive,
        ),
        SignedCellRelabelEntry::new(
            CellRef::new(1, 1),
            CellRef::new(1, 2),
            IncidenceSign::Positive,
        ),
        SignedCellRelabelEntry::new(
            CellRef::new(1, 2),
            CellRef::new(1, 1),
            IncidenceSign::Positive,
        ),
        SignedCellRelabelEntry::new(
            CellRef::new(1, 3),
            CellRef::new(1, 3),
            IncidenceSign::Positive,
        ),
    ]
}

fn reflected_cut_loop_entries() -> Vec<SignedCellRelabelEntry> {
    (0_u32..4)
        .map(|ordinal| {
            SignedCellRelabelEntry::new(
                CellRef::new(0, ordinal),
                CellRef::new(0, 3 - ordinal),
                IncidenceSign::Positive,
            )
        })
        .chain((0_u32..4).map(|ordinal| {
            SignedCellRelabelEntry::new(
                CellRef::new(1, ordinal),
                CellRef::new(1, 3 - ordinal),
                IncidenceSign::Negative,
            )
        }))
        .collect()
}

fn two_phase_preserve_swap_cell_entries() -> Vec<SignedCellRelabelEntry> {
    vec![
        SignedCellRelabelEntry::new(
            CellRef::new(0, 0),
            CellRef::new(0, 2),
            IncidenceSign::Positive,
        ),
        SignedCellRelabelEntry::new(
            CellRef::new(0, 1),
            CellRef::new(0, 3),
            IncidenceSign::Positive,
        ),
        SignedCellRelabelEntry::new(
            CellRef::new(0, 2),
            CellRef::new(0, 0),
            IncidenceSign::Positive,
        ),
        SignedCellRelabelEntry::new(
            CellRef::new(0, 3),
            CellRef::new(0, 1),
            IncidenceSign::Positive,
        ),
        SignedCellRelabelEntry::new(
            CellRef::new(1, 0),
            CellRef::new(1, 1),
            IncidenceSign::Positive,
        ),
        SignedCellRelabelEntry::new(
            CellRef::new(1, 1),
            CellRef::new(1, 0),
            IncidenceSign::Positive,
        ),
    ]
}

fn two_phase_terminal_reverse_cell_entries() -> Vec<SignedCellRelabelEntry> {
    vec![
        SignedCellRelabelEntry::new(
            CellRef::new(0, 0),
            CellRef::new(0, 1),
            IncidenceSign::Positive,
        ),
        SignedCellRelabelEntry::new(
            CellRef::new(0, 1),
            CellRef::new(0, 0),
            IncidenceSign::Positive,
        ),
        SignedCellRelabelEntry::new(
            CellRef::new(0, 2),
            CellRef::new(0, 3),
            IncidenceSign::Positive,
        ),
        SignedCellRelabelEntry::new(
            CellRef::new(0, 3),
            CellRef::new(0, 2),
            IncidenceSign::Positive,
        ),
        SignedCellRelabelEntry::new(
            CellRef::new(1, 0),
            CellRef::new(1, 0),
            IncidenceSign::Negative,
        ),
        SignedCellRelabelEntry::new(
            CellRef::new(1, 1),
            CellRef::new(1, 1),
            IncidenceSign::Negative,
        ),
    ]
}

fn two_phase_preserve_swap_semantics() -> TerminalRelativeSemanticPermutation {
    TerminalRelativeSemanticPermutation::new(
        vec![
            ComponentRelabelEntry::new(
                ConductorComponentId::new("component/a").unwrap(),
                ConductorComponentId::new("component/b").unwrap(),
            ),
            ComponentRelabelEntry::new(
                ConductorComponentId::new("component/b").unwrap(),
                ConductorComponentId::new("component/a").unwrap(),
            ),
        ],
        vec![
            PhaseRelabelEntry::new(
                PhaseId::new("phase/a").unwrap(),
                PhaseId::new("phase/b").unwrap(),
                PhaseCurrentSign::Preserve,
            ),
            PhaseRelabelEntry::new(
                PhaseId::new("phase/b").unwrap(),
                PhaseId::new("phase/a").unwrap(),
                PhaseCurrentSign::Preserve,
            ),
        ],
        vec![
            TerminalRelabelEntry::new(
                PhysicalTerminalId::new("terminal/a-positive").unwrap(),
                PhysicalTerminalId::new("terminal/b-positive").unwrap(),
            ),
            TerminalRelabelEntry::new(
                PhysicalTerminalId::new("terminal/a-return").unwrap(),
                PhysicalTerminalId::new("terminal/b-return").unwrap(),
            ),
            TerminalRelabelEntry::new(
                PhysicalTerminalId::new("terminal/b-positive").unwrap(),
                PhysicalTerminalId::new("terminal/a-positive").unwrap(),
            ),
            TerminalRelabelEntry::new(
                PhysicalTerminalId::new("terminal/b-return").unwrap(),
                PhysicalTerminalId::new("terminal/a-return").unwrap(),
            ),
        ],
    )
}

fn two_phase_terminal_reverse_semantics() -> TerminalRelativeSemanticPermutation {
    TerminalRelativeSemanticPermutation::new(
        vec![
            ComponentRelabelEntry::new(
                ConductorComponentId::new("component/a").unwrap(),
                ConductorComponentId::new("component/a").unwrap(),
            ),
            ComponentRelabelEntry::new(
                ConductorComponentId::new("component/b").unwrap(),
                ConductorComponentId::new("component/b").unwrap(),
            ),
        ],
        vec![
            PhaseRelabelEntry::new(
                PhaseId::new("phase/a").unwrap(),
                PhaseId::new("phase/a").unwrap(),
                PhaseCurrentSign::Reverse,
            ),
            PhaseRelabelEntry::new(
                PhaseId::new("phase/b").unwrap(),
                PhaseId::new("phase/b").unwrap(),
                PhaseCurrentSign::Reverse,
            ),
        ],
        vec![
            TerminalRelabelEntry::new(
                PhysicalTerminalId::new("terminal/a-positive").unwrap(),
                PhysicalTerminalId::new("terminal/a-return").unwrap(),
            ),
            TerminalRelabelEntry::new(
                PhysicalTerminalId::new("terminal/a-return").unwrap(),
                PhysicalTerminalId::new("terminal/a-positive").unwrap(),
            ),
            TerminalRelabelEntry::new(
                PhysicalTerminalId::new("terminal/b-positive").unwrap(),
                PhysicalTerminalId::new("terminal/b-return").unwrap(),
            ),
            TerminalRelabelEntry::new(
                PhysicalTerminalId::new("terminal/b-return").unwrap(),
                PhysicalTerminalId::new("terminal/b-positive").unwrap(),
            ),
        ],
    )
}

fn reversed_semantic_rows(
    semantics: &TerminalRelativeSemanticPermutation,
) -> TerminalRelativeSemanticPermutation {
    let mut components = semantics.components().to_vec();
    let mut phases = semantics.phases().to_vec();
    let mut terminals = semantics.terminals().to_vec();
    components.reverse();
    phases.reverse();
    terminals.reverse();
    TerminalRelativeSemanticPermutation::new(components, phases, terminals)
}

fn two_phase_composed_cell_entries() -> Vec<SignedCellRelabelEntry> {
    vec![
        SignedCellRelabelEntry::new(
            CellRef::new(0, 0),
            CellRef::new(0, 3),
            IncidenceSign::Positive,
        ),
        SignedCellRelabelEntry::new(
            CellRef::new(0, 1),
            CellRef::new(0, 2),
            IncidenceSign::Positive,
        ),
        SignedCellRelabelEntry::new(
            CellRef::new(0, 2),
            CellRef::new(0, 1),
            IncidenceSign::Positive,
        ),
        SignedCellRelabelEntry::new(
            CellRef::new(0, 3),
            CellRef::new(0, 0),
            IncidenceSign::Positive,
        ),
        SignedCellRelabelEntry::new(
            CellRef::new(1, 0),
            CellRef::new(1, 1),
            IncidenceSign::Negative,
        ),
        SignedCellRelabelEntry::new(
            CellRef::new(1, 1),
            CellRef::new(1, 0),
            IncidenceSign::Negative,
        ),
    ]
}

fn two_phase_composed_semantics() -> TerminalRelativeSemanticPermutation {
    TerminalRelativeSemanticPermutation::new(
        two_phase_preserve_swap_semantics().components().to_vec(),
        vec![
            PhaseRelabelEntry::new(
                PhaseId::new("phase/a").unwrap(),
                PhaseId::new("phase/b").unwrap(),
                PhaseCurrentSign::Reverse,
            ),
            PhaseRelabelEntry::new(
                PhaseId::new("phase/b").unwrap(),
                PhaseId::new("phase/a").unwrap(),
                PhaseCurrentSign::Reverse,
            ),
        ],
        vec![
            TerminalRelabelEntry::new(
                PhysicalTerminalId::new("terminal/a-positive").unwrap(),
                PhysicalTerminalId::new("terminal/b-return").unwrap(),
            ),
            TerminalRelabelEntry::new(
                PhysicalTerminalId::new("terminal/a-return").unwrap(),
                PhysicalTerminalId::new("terminal/b-positive").unwrap(),
            ),
            TerminalRelabelEntry::new(
                PhysicalTerminalId::new("terminal/b-positive").unwrap(),
                PhysicalTerminalId::new("terminal/a-return").unwrap(),
            ),
            TerminalRelabelEntry::new(
                PhysicalTerminalId::new("terminal/b-return").unwrap(),
                PhysicalTerminalId::new("terminal/a-positive").unwrap(),
            ),
        ],
    )
}

fn two_phase_current_winding_product(
    winding: &IntegralWindingRepresentative,
    amplitude: &RealCurrentAmplitude,
) -> [f64; 2] {
    assert_eq!(winding.chain().phase(), amplitude.phase());
    let coefficient = winding.chain().coefficients()[0] as f64 * amplitude.value().value();
    match winding.chain().phase().as_str() {
        "phase/a" => [coefficient, 0.0],
        "phase/b" => [0.0, coefficient],
        phase => panic!("unexpected two-phase fixture phase {phase}"),
    }
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
    disconnected_two_phase_pair_with_current_references([
        "current-reference/two-phase",
        "current-reference/two-phase",
        "current-reference/two-phase",
        "current-reference/two-phase",
    ])
}

fn disconnected_two_phase_pair_with_current_references(
    current_references: [&str; 4],
) -> TerminalRelativePair {
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
            terminal_for_with_current_reference(
                &complex,
                0,
                "a-positive",
                "component/a",
                "phase/a",
                TerminalRole::Driven,
                TerminalOrientation::OutOfConductor,
                OrientationMapSign::Preserve,
                37,
                stable(current_references[0]),
            ),
            terminal_for_with_current_reference(
                &complex,
                1,
                "a-return",
                "component/a",
                "phase/a",
                TerminalRole::ReturnReference,
                TerminalOrientation::IntoConductor,
                OrientationMapSign::Reverse,
                37,
                stable(current_references[1]),
            ),
            terminal_for_with_current_reference(
                &complex,
                2,
                "b-positive",
                "component/b",
                "phase/b",
                TerminalRole::Driven,
                TerminalOrientation::OutOfConductor,
                OrientationMapSign::Preserve,
                37,
                stable(current_references[2]),
            ),
            terminal_for_with_current_reference(
                &complex,
                3,
                "b-return",
                "component/b",
                "phase/b",
                TerminalRole::ReturnReference,
                TerminalOrientation::IntoConductor,
                OrientationMapSign::Reverse,
                37,
                stable(current_references[3]),
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

#[test]
fn i13_2a_015_parallel_edge_relabel_is_canonical_invertible_and_composable() {
    let pair = terminal_cut_loop_pair();
    let canonical_entries = parallel_edge_relabel_entries();
    let mut reversed_entries = canonical_entries.clone();
    reversed_entries.reverse();

    let relabel = TerminalRelativeSignedRelabel::try_new(&pair, &pair, reversed_entries)
        .expect("parallel-edge permutation is an exact relabeling");
    let canonical_replay =
        TerminalRelativeSignedRelabel::try_new(&pair, &pair, canonical_entries.clone())
            .expect("canonical declaration replays");
    assert_eq!(relabel.entries(), canonical_entries.as_slice());
    assert_eq!(relabel.identity(), canonical_replay.identity());

    let phase = PhaseId::new("phase/a").unwrap();
    let chain = IntegralRelativeChain::try_new(&pair, phase.clone(), 1, vec![2, -1, 3, 4])
        .expect("fixture chain");
    let transported_chain = relabel
        .transport_integral_chain(&pair, &pair, &chain)
        .expect("transport chain");
    assert_eq!(transported_chain.coefficients(), &[2, 3, -1, 4]);

    let cochain = IntegralRelativeCochain::try_new(&pair, phase.clone(), 1, vec![11, 13, 17, 19])
        .expect("fixture cochain");
    let transported_cochain = relabel
        .transport_integral_cochain(&pair, &pair, &cochain)
        .expect("transport cochain");
    assert_eq!(transported_cochain.coefficients(), &[11, 17, 13, 19]);

    let winding = IntegralWindingRepresentative::try_new(&pair, phase, vec![1, 1, 0, 1])
        .expect("fixture winding cycle");
    let transported_winding = relabel
        .transport_winding_representative(&pair, &pair, &winding)
        .expect("transport winding representative");
    assert_eq!(transported_winding.chain().coefficients(), &[1, 0, 1, 1]);

    let inverse = relabel.inverse(&pair, &pair).expect("inverse relabeling");
    assert_eq!(inverse.identity(), relabel.identity());
    assert_eq!(
        inverse
            .transport_integral_chain(&pair, &pair, &transported_chain)
            .expect("inverse chain transport"),
        chain
    );
    assert_eq!(
        inverse
            .transport_winding_representative(&pair, &pair, &transported_winding)
            .expect("inverse winding transport")
            .identity(),
        winding.identity()
    );

    let composed = relabel
        .compose(&relabel, &pair, &pair, &pair)
        .expect("self-composition is the identity permutation");
    let identity_entries = canonical_entries
        .iter()
        .map(|entry| {
            SignedCellRelabelEntry::new(entry.source(), entry.source(), IncidenceSign::Positive)
        })
        .collect();
    let identity = TerminalRelativeSignedRelabel::try_new(&pair, &pair, identity_entries)
        .expect("explicit identity relabeling");
    assert_eq!(composed.identity(), identity.identity());
    assert_eq!(
        composed
            .transport_integral_chain(&pair, &pair, &chain)
            .expect("composed chain transport"),
        chain
    );
}

#[test]
fn i13_2a_016_orientation_reflection_preserves_relative_naturality() {
    let source = terminal_cut_loop_pair();
    let target = terminal_cut_loop_pair_with_terminals(3, 0);
    let reflection =
        TerminalRelativeSignedRelabel::try_new(&source, &target, reflected_cut_loop_entries())
            .expect("orientation reflection preserves terminal-relative semantics");
    let phase = PhaseId::new("phase/a").unwrap();

    let chain = IntegralRelativeChain::try_new(&source, phase.clone(), 1, vec![2, -1, 3, 4])
        .expect("fixture chain");
    let source_boundary = source.boundary(&chain).expect("source boundary");
    assert_eq!(source_boundary.coefficients(), &[0, -2]);
    let transported_chain = reflection
        .transport_integral_chain(&source, &target, &chain)
        .expect("reflect chain");
    assert_eq!(transported_chain.coefficients(), &[-4, -3, 1, -2]);
    let transported_boundary = reflection
        .transport_integral_chain(&source, &target, &source_boundary)
        .expect("reflect source boundary");
    assert_eq!(transported_boundary.coefficients(), &[-2, 0]);
    assert_eq!(
        target
            .boundary(&transported_chain)
            .expect("target boundary"),
        transported_boundary
    );

    let cochain = IntegralRelativeCochain::try_new(&source, phase.clone(), 0, vec![2, 5])
        .expect("fixture cochain");
    let source_coboundary = source
        .integral_coboundary(&cochain)
        .expect("source coboundary");
    assert_eq!(source_coboundary.coefficients(), &[2, 3, 3, -5]);
    let transported_cochain = reflection
        .transport_integral_cochain(&source, &target, &cochain)
        .expect("reflect cochain");
    assert_eq!(transported_cochain.coefficients(), &[5, 2]);
    let transported_coboundary = reflection
        .transport_integral_cochain(&source, &target, &source_coboundary)
        .expect("reflect source coboundary");
    assert_eq!(transported_coboundary.coefficients(), &[5, -3, -3, -2]);
    assert_eq!(
        target
            .integral_coboundary(&transported_cochain)
            .expect("target coboundary"),
        transported_coboundary
    );

    assert_eq!(source.integral_pairing(&source_coboundary, &chain), Ok(-10));
    assert_eq!(source.integral_pairing(&cochain, &source_boundary), Ok(-10));
    assert_eq!(
        target.integral_pairing(&transported_coboundary, &transported_chain),
        Ok(-10)
    );
    assert_eq!(
        target.integral_pairing(&transported_cochain, &transported_boundary),
        Ok(-10)
    );

    let winding = IntegralWindingRepresentative::try_new(&source, phase, vec![1, 1, 0, 1])
        .expect("fixture winding cycle");
    let transported_winding = reflection
        .transport_winding_representative(&source, &target, &winding)
        .expect("reflect winding representative");
    assert_eq!(transported_winding.chain().coefficients(), &[-1, 0, -1, -1]);

    let inverse = reflection
        .inverse(&source, &target)
        .expect("inverse reflection");
    assert_eq!(
        inverse
            .transport_integral_chain(&target, &source, &transported_chain)
            .expect("inverse chain transport"),
        chain
    );
    assert_eq!(
        inverse
            .transport_winding_representative(&target, &source, &transported_winding)
            .expect("inverse winding transport")
            .identity(),
        winding.identity()
    );
}

#[test]
fn i13_2a_017_signed_relabel_admission_refuses_partial_duplicate_and_non_chain_maps() {
    let pair = terminal_cut_loop_pair();

    let mut missing = parallel_edge_relabel_entries();
    missing.pop();
    assert_eq!(
        TerminalRelativeSignedRelabel::try_new(&pair, &pair, missing),
        Err(TerminalRelativeSignedRelabelError::EntryCountMismatch {
            expected: 8,
            actual: 7,
        })
    );

    let mut duplicate_source = parallel_edge_relabel_entries();
    duplicate_source[7] = SignedCellRelabelEntry::new(
        CellRef::new(1, 2),
        CellRef::new(1, 3),
        IncidenceSign::Positive,
    );
    assert_eq!(
        TerminalRelativeSignedRelabel::try_new(&pair, &pair, duplicate_source),
        Err(TerminalRelativeSignedRelabelError::DuplicateSourceCell {
            cell: CellRef::new(1, 2),
        })
    );

    let mut duplicate_target = parallel_edge_relabel_entries();
    duplicate_target[7] = SignedCellRelabelEntry::new(
        CellRef::new(1, 3),
        CellRef::new(1, 2),
        IncidenceSign::Positive,
    );
    assert_eq!(
        TerminalRelativeSignedRelabel::try_new(&pair, &pair, duplicate_target),
        Err(TerminalRelativeSignedRelabelError::DuplicateTargetCell {
            cell: CellRef::new(1, 2),
        })
    );

    let reflected_target = terminal_cut_loop_pair_with_terminals(3, 0);
    let mut wrong_sign = reflected_cut_loop_entries();
    wrong_sign[4] = SignedCellRelabelEntry::new(
        CellRef::new(1, 0),
        CellRef::new(1, 3),
        IncidenceSign::Positive,
    );
    assert!(matches!(
        TerminalRelativeSignedRelabel::try_new(&pair, &reflected_target, wrong_sign),
        Err(TerminalRelativeSignedRelabelError::MappedIncidenceMismatch { .. })
    ));
}

#[test]
fn i13_2a_018_reflection_to_same_pair_refuses_terminal_support_mismatch() {
    let pair = terminal_cut_loop_pair();
    assert_eq!(
        TerminalRelativeSignedRelabel::try_new(&pair, &pair, reflected_cut_loop_entries()),
        Err(TerminalRelativeSignedRelabelError::MappedSupportMismatch {
            role: "physical terminal support",
            owner: Some("terminal/loop-positive".to_owned()),
            cell: CellRef::new(0, 0),
            expected_mapped: false,
            actual_target: true,
        })
    );
}

#[test]
fn i13_2a_019_signed_transport_refuses_exact_coefficient_overflow() {
    let source = terminal_cut_loop_pair();
    let target = terminal_cut_loop_pair_with_terminals(3, 0);
    let reflection =
        TerminalRelativeSignedRelabel::try_new(&source, &target, reflected_cut_loop_entries())
            .expect("orientation reflection");
    let chain = IntegralRelativeChain::try_new(
        &source,
        PhaseId::new("phase/a").unwrap(),
        1,
        vec![i64::MIN, 0, 0, 0],
    )
    .expect("minimum exact coefficient remains an admitted source value");
    assert_eq!(
        reflection.transport_integral_chain(&source, &target, &chain),
        Err(TerminalRelativeSignedRelabelError::CoefficientOverflow {
            cell: CellRef::new(1, 0),
        })
    );
}

#[test]
fn i13_2a_020_physical_phase_swap_is_canonical_and_transports_each_coefficient_sector() {
    let pair = disconnected_two_phase_pair();
    let canonical_cells = two_phase_preserve_swap_cell_entries();
    let mut reversed_cells = canonical_cells.clone();
    reversed_cells.reverse();
    let canonical_semantics = two_phase_preserve_swap_semantics();
    let relabel = TerminalRelativePhysicalRelabel::try_new(
        &pair,
        &pair,
        reversed_cells,
        reversed_semantic_rows(&canonical_semantics),
    )
    .expect("physical phase swap");
    let canonical_replay = TerminalRelativePhysicalRelabel::try_new(
        &pair,
        &pair,
        canonical_cells.clone(),
        canonical_semantics.clone(),
    )
    .expect("canonical phase-swap replay");

    assert_eq!(relabel.cell_entries(), canonical_cells.as_slice());
    assert_eq!(
        relabel.semantic_permutation().components(),
        canonical_semantics.components()
    );
    assert_eq!(
        relabel.semantic_permutation().phases(),
        canonical_semantics.phases()
    );
    assert_eq!(
        relabel.semantic_permutation().terminals(),
        canonical_semantics.terminals()
    );
    assert_eq!(relabel.identity(), canonical_replay.identity());

    let phase_a = PhaseId::new("phase/a").unwrap();
    let phase_b = PhaseId::new("phase/b").unwrap();
    let chain =
        IntegralRelativeChain::try_new(&pair, phase_a.clone(), 1, vec![7]).expect("phase-a chain");
    let transported_chain = relabel
        .transport_integral_chain(&pair, &pair, &chain)
        .expect("swap chain phase");
    assert_eq!(transported_chain.phase(), &phase_b);
    assert_eq!(transported_chain.coefficients(), &[7]);

    let cochain = IntegralRelativeCochain::try_new(&pair, phase_a.clone(), 1, vec![11])
        .expect("phase-a cochain");
    let transported_cochain = relabel
        .transport_integral_cochain(&pair, &pair, &cochain)
        .expect("swap cochain phase");
    assert_eq!(transported_cochain.phase(), &phase_b);
    assert_eq!(transported_cochain.coefficients(), &[11]);

    let winding = IntegralWindingRepresentative::try_new(&pair, phase_a.clone(), vec![3])
        .expect("phase-a winding");
    let amplitude = RealCurrentAmplitude::try_new(
        PhysicalObjectId::new("object/current-a-before-phase-swap").unwrap(),
        &pair,
        phase_a,
        Current::new(2.5),
    )
    .expect("phase-a current amplitude");
    let transported_winding = relabel
        .transport_winding_representative(&pair, &pair, &winding)
        .expect("swap winding phase");
    let target_amplitude_id = PhysicalObjectId::new("object/current-b-after-phase-swap").unwrap();
    let transported_amplitude = relabel
        .transport_current_amplitude(&pair, &pair, &amplitude, target_amplitude_id.clone())
        .expect("swap current-amplitude phase");
    assert_eq!(transported_winding.chain().phase(), &phase_b);
    assert_eq!(transported_winding.chain().coefficients(), &[3]);
    assert_eq!(transported_amplitude.phase(), &phase_b);
    assert_eq!(
        transported_amplitude.value().value().to_bits(),
        2.5_f64.to_bits()
    );
    assert_eq!(transported_amplitude.id(), &target_amplitude_id);
    assert_eq!(
        two_phase_current_winding_product(&winding, &amplitude),
        [7.5, 0.0]
    );
    assert_eq!(
        two_phase_current_winding_product(&transported_winding, &transported_amplitude),
        [0.0, 7.5]
    );

    let winding_b = IntegralWindingRepresentative::try_new(&pair, phase_b.clone(), vec![-5])
        .expect("phase-b winding");
    let amplitude_b = RealCurrentAmplitude::try_new(
        PhysicalObjectId::new("object/current-b-before-phase-swap").unwrap(),
        &pair,
        phase_b,
        Current::new(1.5),
    )
    .expect("phase-b current amplitude");
    let transported_winding_b = relabel
        .transport_winding_representative(&pair, &pair, &winding_b)
        .expect("swap phase-b winding");
    let transported_amplitude_b = relabel
        .transport_current_amplitude(
            &pair,
            &pair,
            &amplitude_b,
            PhysicalObjectId::new("object/current-a-after-phase-swap").unwrap(),
        )
        .expect("swap phase-b current amplitude");
    assert_eq!(
        two_phase_current_winding_product(&winding_b, &amplitude_b),
        [0.0, -7.5]
    );
    assert_eq!(
        two_phase_current_winding_product(&transported_winding_b, &transported_amplitude_b),
        [-7.5, 0.0]
    );
}

#[test]
fn i13_2a_021_terminal_current_reversal_separates_cell_and_physical_signs() {
    let pair = disconnected_two_phase_pair();
    let reversal = TerminalRelativePhysicalRelabel::try_new(
        &pair,
        &pair,
        two_phase_terminal_reverse_cell_entries(),
        two_phase_terminal_reverse_semantics(),
    )
    .expect("terminal/current reversal");
    let phase_a = PhaseId::new("phase/a").unwrap();

    let chain =
        IntegralRelativeChain::try_new(&pair, phase_a.clone(), 1, vec![7]).expect("phase-a chain");
    assert_eq!(
        reversal
            .transport_integral_chain(&pair, &pair, &chain)
            .expect("reverse raw chain")
            .coefficients(),
        &[-7]
    );
    let cochain = IntegralRelativeCochain::try_new(&pair, phase_a.clone(), 1, vec![11])
        .expect("phase-a cochain");
    assert_eq!(
        reversal
            .transport_integral_cochain(&pair, &pair, &cochain)
            .expect("reverse raw cochain")
            .coefficients(),
        &[-11]
    );

    let winding = IntegralWindingRepresentative::try_new(&pair, phase_a.clone(), vec![3])
        .expect("phase-a winding");
    let amplitude = RealCurrentAmplitude::try_new(
        PhysicalObjectId::new("object/current-a-before-terminal-reversal").unwrap(),
        &pair,
        phase_a.clone(),
        Current::new(2.5),
    )
    .expect("phase-a amplitude");
    let transported_winding = reversal
        .transport_winding_representative(&pair, &pair, &winding)
        .expect("combine cell and current reversals");
    let target_amplitude_id =
        PhysicalObjectId::new("object/current-a-after-terminal-reversal").unwrap();
    let transported_amplitude = reversal
        .transport_current_amplitude(&pair, &pair, &amplitude, target_amplitude_id.clone())
        .expect("reverse physical current coordinate");
    assert_eq!(transported_winding.chain().coefficients(), &[3]);
    assert_eq!(transported_amplitude.phase(), &phase_a);
    assert_eq!(
        transported_amplitude.value().value().to_bits(),
        (-2.5_f64).to_bits()
    );
    assert_eq!(transported_amplitude.id(), &target_amplitude_id);
    assert_eq!(
        two_phase_current_winding_product(&winding, &amplitude),
        [7.5, 0.0]
    );
    assert_eq!(
        two_phase_current_winding_product(&transported_winding, &transported_amplitude),
        [-7.5, 0.0]
    );

    let minimum_winding =
        IntegralWindingRepresentative::try_new(&pair, phase_a.clone(), vec![i64::MIN])
            .expect("minimum exact winding coefficient");
    assert_eq!(
        reversal
            .transport_winding_representative(&pair, &pair, &minimum_winding)
            .expect("two reversals avoid an intermediate negation")
            .chain()
            .coefficients(),
        &[i64::MIN]
    );
    let minimum_chain = IntegralRelativeChain::try_new(&pair, phase_a, 1, vec![i64::MIN])
        .expect("minimum exact raw-chain coefficient");
    assert_eq!(
        reversal.transport_integral_chain(&pair, &pair, &minimum_chain),
        Err(TerminalRelativePhysicalRelabelError::CoefficientOverflow {
            cell: CellRef::new(1, 0),
        })
    );
}

#[test]
fn i13_2a_022_physical_relabel_inverse_and_composition_obey_the_same_exact_action() {
    let pair = disconnected_two_phase_pair();
    let phase_swap = TerminalRelativePhysicalRelabel::try_new(
        &pair,
        &pair,
        two_phase_preserve_swap_cell_entries(),
        two_phase_preserve_swap_semantics(),
    )
    .expect("phase swap");
    let terminal_reversal = TerminalRelativePhysicalRelabel::try_new(
        &pair,
        &pair,
        two_phase_terminal_reverse_cell_entries(),
        two_phase_terminal_reverse_semantics(),
    )
    .expect("terminal/current reversal");
    let direct_composition = TerminalRelativePhysicalRelabel::try_new(
        &pair,
        &pair,
        two_phase_composed_cell_entries(),
        two_phase_composed_semantics(),
    )
    .expect("direct composed action");
    let phase_then_terminal = phase_swap
        .compose(&terminal_reversal, &pair, &pair, &pair)
        .expect("terminal reversal after phase swap");
    let terminal_then_phase = terminal_reversal
        .compose(&phase_swap, &pair, &pair, &pair)
        .expect("phase swap after terminal reversal");
    assert_eq!(
        phase_then_terminal.identity(),
        direct_composition.identity()
    );
    assert_eq!(
        terminal_then_phase.identity(),
        direct_composition.identity()
    );

    assert_eq!(
        phase_swap
            .inverse(&pair, &pair)
            .expect("phase-swap inverse")
            .identity(),
        phase_swap.identity()
    );
    assert_eq!(
        terminal_reversal
            .inverse(&pair, &pair)
            .expect("terminal-reversal inverse")
            .identity(),
        terminal_reversal.identity()
    );
    assert_eq!(
        direct_composition
            .inverse(&pair, &pair)
            .expect("composed inverse")
            .identity(),
        direct_composition.identity()
    );
    let identity = TerminalRelativePhysicalRelabel::identity_on(&pair).expect("physical identity");
    assert_eq!(
        phase_swap
            .compose(&phase_swap, &pair, &pair, &pair)
            .expect("phase swap squared")
            .identity(),
        identity.identity()
    );
    assert_eq!(
        terminal_reversal
            .compose(&terminal_reversal, &pair, &pair, &pair)
            .expect("terminal reversal squared")
            .identity(),
        identity.identity()
    );
    assert_eq!(
        direct_composition
            .compose(&direct_composition, &pair, &pair, &pair)
            .expect("composed involution squared")
            .identity(),
        identity.identity()
    );

    let phase_a = PhaseId::new("phase/a").unwrap();
    let winding = IntegralWindingRepresentative::try_new(&pair, phase_a.clone(), vec![3])
        .expect("phase-a winding");
    let amplitude = RealCurrentAmplitude::try_new(
        PhysicalObjectId::new("object/current-a-before-composition").unwrap(),
        &pair,
        phase_a,
        Current::new(2.5),
    )
    .expect("phase-a amplitude");
    let intermediate_winding = phase_swap
        .transport_winding_representative(&pair, &pair, &winding)
        .expect("first winding transport");
    let sequential_winding = terminal_reversal
        .transport_winding_representative(&pair, &pair, &intermediate_winding)
        .expect("second winding transport");
    let direct_winding = direct_composition
        .transport_winding_representative(&pair, &pair, &winding)
        .expect("direct winding transport");
    assert_eq!(sequential_winding.identity(), direct_winding.identity());

    let intermediate_amplitude = phase_swap
        .transport_current_amplitude(
            &pair,
            &pair,
            &amplitude,
            PhysicalObjectId::new("object/current-b-between-composed-actions").unwrap(),
        )
        .expect("first amplitude transport");
    let final_amplitude_id =
        PhysicalObjectId::new("object/current-b-after-composed-actions").unwrap();
    let sequential_amplitude = terminal_reversal
        .transport_current_amplitude(
            &pair,
            &pair,
            &intermediate_amplitude,
            final_amplitude_id.clone(),
        )
        .expect("second amplitude transport");
    let direct_amplitude = direct_composition
        .transport_current_amplitude(&pair, &pair, &amplitude, final_amplitude_id.clone())
        .expect("direct amplitude transport");
    assert_eq!(sequential_amplitude, direct_amplitude);
    assert_eq!(direct_amplitude.id(), &final_amplitude_id);
    assert_eq!(
        two_phase_current_winding_product(&direct_winding, &direct_amplitude),
        [0.0, -7.5]
    );
}

#[test]
fn i13_2a_023_physical_relabel_admission_refuses_incomplete_or_noncommuting_data() {
    let pair = disconnected_two_phase_pair();
    let preserve_swap = two_phase_preserve_swap_semantics();

    let missing_component = TerminalRelativeSemanticPermutation::new(
        preserve_swap.components()[..1].to_vec(),
        preserve_swap.phases().to_vec(),
        preserve_swap.terminals().to_vec(),
    );
    assert_eq!(
        TerminalRelativePhysicalRelabel::try_new(
            &pair,
            &pair,
            two_phase_preserve_swap_cell_entries(),
            missing_component,
        ),
        Err(
            TerminalRelativePhysicalRelabelError::SemanticEntryCountMismatch {
                role: "conductor component",
                expected: 2,
                actual: 1,
            }
        )
    );

    let duplicate_phase_source = TerminalRelativeSemanticPermutation::new(
        preserve_swap.components().to_vec(),
        vec![
            PhaseRelabelEntry::new(
                PhaseId::new("phase/a").unwrap(),
                PhaseId::new("phase/b").unwrap(),
                PhaseCurrentSign::Preserve,
            ),
            PhaseRelabelEntry::new(
                PhaseId::new("phase/a").unwrap(),
                PhaseId::new("phase/a").unwrap(),
                PhaseCurrentSign::Preserve,
            ),
        ],
        preserve_swap.terminals().to_vec(),
    );
    assert_eq!(
        TerminalRelativePhysicalRelabel::try_new(
            &pair,
            &pair,
            two_phase_preserve_swap_cell_entries(),
            duplicate_phase_source,
        ),
        Err(
            TerminalRelativePhysicalRelabelError::DuplicateSemanticSource {
                role: "phase",
                id: "phase/a".to_owned(),
            }
        )
    );

    let duplicate_phase_target = TerminalRelativeSemanticPermutation::new(
        preserve_swap.components().to_vec(),
        vec![
            PhaseRelabelEntry::new(
                PhaseId::new("phase/a").unwrap(),
                PhaseId::new("phase/b").unwrap(),
                PhaseCurrentSign::Preserve,
            ),
            PhaseRelabelEntry::new(
                PhaseId::new("phase/b").unwrap(),
                PhaseId::new("phase/b").unwrap(),
                PhaseCurrentSign::Preserve,
            ),
        ],
        preserve_swap.terminals().to_vec(),
    );
    assert_eq!(
        TerminalRelativePhysicalRelabel::try_new(
            &pair,
            &pair,
            two_phase_preserve_swap_cell_entries(),
            duplicate_phase_target,
        ),
        Err(
            TerminalRelativePhysicalRelabelError::DuplicateSemanticTarget {
                role: "phase",
                id: "phase/b".to_owned(),
            }
        )
    );

    let wrong_current_parity = TerminalRelativeSemanticPermutation::new(
        preserve_swap.components().to_vec(),
        vec![
            PhaseRelabelEntry::new(
                PhaseId::new("phase/a").unwrap(),
                PhaseId::new("phase/b").unwrap(),
                PhaseCurrentSign::Reverse,
            ),
            PhaseRelabelEntry::new(
                PhaseId::new("phase/b").unwrap(),
                PhaseId::new("phase/a").unwrap(),
                PhaseCurrentSign::Reverse,
            ),
        ],
        preserve_swap.terminals().to_vec(),
    );
    assert_eq!(
        TerminalRelativePhysicalRelabel::try_new(
            &pair,
            &pair,
            two_phase_preserve_swap_cell_entries(),
            wrong_current_parity,
        ),
        Err(
            TerminalRelativePhysicalRelabelError::TerminalConventionMismatch {
                source_terminal: "terminal/a-positive".to_owned(),
                target_terminal: "terminal/b-positive".to_owned(),
                field: "terminal role/current-sign parity",
            }
        )
    );

    let changed_current_reference_target = disconnected_two_phase_pair_with_current_references([
        "current-reference/two-phase",
        "current-reference/two-phase",
        "current-reference/changed-only-for-target-b-positive",
        "current-reference/two-phase",
    ]);
    assert_eq!(
        TerminalRelativePhysicalRelabel::try_new(
            &pair,
            &changed_current_reference_target,
            two_phase_preserve_swap_cell_entries(),
            preserve_swap.clone(),
        ),
        Err(
            TerminalRelativePhysicalRelabelError::TerminalConventionMismatch {
                source_terminal: "terminal/a-positive".to_owned(),
                target_terminal: "terminal/b-positive".to_owned(),
                field: "current reference",
            }
        )
    );

    let identity_cells = TerminalRelativePhysicalRelabel::identity_on(&pair)
        .expect("physical identity")
        .cell_entries()
        .to_vec();
    assert_eq!(
        TerminalRelativePhysicalRelabel::try_new(
            &pair,
            &pair,
            identity_cells,
            preserve_swap.clone(),
        ),
        Err(
            TerminalRelativePhysicalRelabelError::MappedSemanticSupportMismatch {
                role: "conductor component support",
                source_owner: "component/a".to_owned(),
                target_owner: "component/b".to_owned(),
                cell: CellRef::new(0, 0),
                expected_mapped: true,
                actual_target: false,
            }
        )
    );

    let wrong_phase_square = TerminalRelativeSemanticPermutation::new(
        preserve_swap.components().to_vec(),
        vec![
            PhaseRelabelEntry::new(
                PhaseId::new("phase/a").unwrap(),
                PhaseId::new("phase/a").unwrap(),
                PhaseCurrentSign::Preserve,
            ),
            PhaseRelabelEntry::new(
                PhaseId::new("phase/b").unwrap(),
                PhaseId::new("phase/b").unwrap(),
                PhaseCurrentSign::Preserve,
            ),
        ],
        preserve_swap.terminals().to_vec(),
    );
    assert_eq!(
        TerminalRelativePhysicalRelabel::try_new(
            &pair,
            &pair,
            two_phase_preserve_swap_cell_entries(),
            wrong_phase_square,
        ),
        Err(
            TerminalRelativePhysicalRelabelError::PhaseComponentSquareMismatch {
                source_phase: "phase/a".to_owned(),
                target_phase: "phase/a".to_owned(),
                expected_target_component: "component/b".to_owned(),
                actual_target_component: "component/a".to_owned(),
            }
        )
    );
}
