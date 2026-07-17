//! G0/G4 first-tranche battery for I13.2b exact integral topology algebra.

#![cfg(feature = "moonshot-integral-topology")]

use fs_couple::{CoordinateBinding, PortKind, PortOrientation, PortTimestamp, StableId};
use fs_feec::integral_topology::{
    ExactAlgebraBudget, ExactIntegerMatrix, IntegralTopologyError, IntegralTopologyFailureClass,
    KernelCoordinateBudget, MatrixRole, SmithConstructionBudget, SmithNormalFormWitness,
    SmithWitnessStage, TerminalRelativeBoundaryBudget, TerminalRelativeBoundaryMatrix,
    TopologyApplicability, VerifiedSmithNormalForm, construct_smith_normal_form,
    construct_smith_normal_form_with_checkpoint, extract_terminal_relative_boundary_matrix,
    extract_terminal_relative_boundary_matrix_with_checkpoint, verify_smith_normal_form,
    verify_smith_normal_form_with_checkpoint, verify_terminal_relative_kernel_transport,
    verify_terminal_relative_kernel_transport_with_checkpoint,
};
use fs_feec::terminal_relative::{
    BoundaryIncidence, CellRef, CellularSubcomplex, ConductorComponent, ConductorComponentId,
    FiniteCellComplex, IncidenceSign, IntegralRelativeChain, IntegralRelativeCochain,
    OrientationMapSign, PhaseId, PhysicalTerminal, PhysicalTerminalId, PresentedMachinePortRef,
    TerminalOrientation, TerminalPortCoordinate, TerminalPortTrivialization, TerminalRelativePair,
    TerminalRole, TrivializationId,
};

fn budget(max_extent: usize, max_work: u128) -> ExactAlgebraBudget {
    let matrix_entries = max_extent * max_extent;
    ExactAlgebraBudget::new(
        max_extent,
        max_extent,
        matrix_entries,
        6 * matrix_entries,
        matrix_entries,
        max_work,
    )
}

fn matrix(rows: usize, cols: usize, entries: &[i128]) -> ExactIntegerMatrix {
    ExactIntegerMatrix::try_new(rows, cols, entries.to_vec(), ExactAlgebraBudget::default())
        .expect("fixture matrix")
}

fn identity(size: usize) -> ExactIntegerMatrix {
    let mut entries = vec![0_i128; size * size];
    for index in 0..size {
        entries[index * size + index] = 1;
    }
    matrix(size, size, &entries)
}

fn rank_one_witness() -> (ExactIntegerMatrix, SmithNormalFormWitness) {
    // Row 1 -= 2*row 0, then col 1 -= 2*col 0.
    let source = matrix(2, 2, &[2, 4, 4, 8]);
    let diagonal = matrix(2, 2, &[2, 0, 0, 0]);
    let left = matrix(2, 2, &[1, 0, -2, 1]);
    let left_inverse = matrix(2, 2, &[1, 0, 2, 1]);
    let right = matrix(2, 2, &[1, -2, 0, 1]);
    let right_inverse = matrix(2, 2, &[1, 2, 0, 1]);
    (
        source,
        SmithNormalFormWitness::new(diagonal, left, left_inverse, right, right_inverse),
    )
}

fn stable(value: &str) -> StableId {
    StableId::new(value).expect("fixture stable id")
}

fn subcomplex(
    ambient: &FiniteCellComplex,
    id: &str,
    cells: impl IntoIterator<Item = CellRef>,
) -> CellularSubcomplex {
    CellularSubcomplex::try_new(stable(id), cells, ambient).expect("fixture subcomplex")
}

fn loop_terminal(
    ambient: &FiniteCellComplex,
    ordinal: u32,
    id: &str,
    role: TerminalRole,
    orientation: TerminalOrientation,
    sign: OrientationMapSign,
) -> PhysicalTerminal {
    let port = PortKind::ElectricalVoltageCurrent
        .scalar_seed_schema(
            stable(&format!("port/{id}")),
            CoordinateBinding::new(
                stable("basis/winding-terminal"),
                stable("frame/winding-terminal"),
                PortOrientation::OutwardFromOwner,
            ),
            PortTimestamp::new(stable("clock/electrical"), 31),
        )
        .expect("electrical port");
    PhysicalTerminal::new(
        PhysicalTerminalId::new(format!("terminal/{id}")).expect("terminal id"),
        subcomplex(
            ambient,
            &format!("support/{id}"),
            [CellRef::new(0, ordinal)],
        ),
        ConductorComponentId::new("component/winding").expect("component id"),
        PhaseId::new("phase/a").expect("phase id"),
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

#[allow(clippy::too_many_lines)]
fn terminal_cut_loop_pair(reverse_declarations: bool) -> TerminalRelativePair {
    let mut incidences = vec![
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
    ];
    if reverse_declarations {
        incidences.reverse();
    }
    let complex =
        FiniteCellComplex::try_new(1, vec![4, 4], incidences).expect("terminal-cut loop graph");
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
        ConductorComponentId::new("component/winding").expect("component id"),
        conductor.clone(),
    )
    .expect("component");
    let mut terminals = vec![
        loop_terminal(
            &complex,
            0,
            "loop-positive",
            TerminalRole::Driven,
            TerminalOrientation::OutOfConductor,
            OrientationMapSign::Preserve,
        ),
        loop_terminal(
            &complex,
            3,
            "loop-return",
            TerminalRole::ReturnReference,
            TerminalOrientation::IntoConductor,
            OrientationMapSign::Reverse,
        ),
    ];
    if reverse_declarations {
        terminals.reverse();
    }
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
        terminals,
    )
    .expect("terminal-cut loop pair")
}

fn surface_terminal(
    ambient: &FiniteCellComplex,
    edge: u32,
    vertices: [u32; 2],
    id: &str,
    role: TerminalRole,
    orientation: TerminalOrientation,
    sign: OrientationMapSign,
) -> PhysicalTerminal {
    let port = PortKind::ElectricalVoltageCurrent
        .scalar_seed_schema(
            stable(&format!("port/{id}")),
            CoordinateBinding::new(
                stable("basis/surface-terminal"),
                stable("frame/surface-terminal"),
                PortOrientation::OutwardFromOwner,
            ),
            PortTimestamp::new(stable("clock/electrical"), 47),
        )
        .expect("surface electrical port");
    PhysicalTerminal::new(
        PhysicalTerminalId::new(format!("terminal/{id}")).expect("terminal id"),
        subcomplex(
            ambient,
            &format!("support/{id}"),
            [
                CellRef::new(0, vertices[0]),
                CellRef::new(0, vertices[1]),
                CellRef::new(1, edge),
            ],
        ),
        ConductorComponentId::new("component/surface").expect("component id"),
        PhaseId::new("phase/surface").expect("phase id"),
        role,
        orientation,
        TerminalPortCoordinate::Flow,
        port.clone(),
        PresentedMachinePortRef::try_new(
            stable("org.frankensim.fs-ir.machine.graph.v1"),
            1,
            [0x57; 32],
            stable("machine-owner/surface"),
            stable(&format!("port/{id}")),
            stable(&format!("machine-terminal/{id}-voltage")),
            stable(&format!("machine-terminal/{id}-current")),
        )
        .expect("presented Machine-IR port"),
        TerminalPortTrivialization::new(
            TrivializationId::new(format!("trivialization/{id}")).expect("trivialization id"),
            port.id().clone(),
            sign,
            stable("voltage-reference/surface-common"),
            stable("current-reference/surface-common"),
        ),
    )
    .expect("surface terminal")
}

#[allow(clippy::too_many_lines)]
fn centered_surface_pair() -> TerminalRelativePair {
    let complex = FiniteCellComplex::try_new(
        2,
        vec![5, 8, 4],
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
                CellRef::new(0, 2),
                CellRef::new(1, 2),
                IncidenceSign::Negative,
            ),
            BoundaryIncidence::new(
                CellRef::new(0, 3),
                CellRef::new(1, 2),
                IncidenceSign::Positive,
            ),
            BoundaryIncidence::new(
                CellRef::new(0, 3),
                CellRef::new(1, 3),
                IncidenceSign::Negative,
            ),
            BoundaryIncidence::new(
                CellRef::new(0, 0),
                CellRef::new(1, 3),
                IncidenceSign::Positive,
            ),
            BoundaryIncidence::new(
                CellRef::new(0, 4),
                CellRef::new(1, 4),
                IncidenceSign::Negative,
            ),
            BoundaryIncidence::new(
                CellRef::new(0, 0),
                CellRef::new(1, 4),
                IncidenceSign::Positive,
            ),
            BoundaryIncidence::new(
                CellRef::new(0, 4),
                CellRef::new(1, 5),
                IncidenceSign::Negative,
            ),
            BoundaryIncidence::new(
                CellRef::new(0, 1),
                CellRef::new(1, 5),
                IncidenceSign::Positive,
            ),
            BoundaryIncidence::new(
                CellRef::new(0, 4),
                CellRef::new(1, 6),
                IncidenceSign::Negative,
            ),
            BoundaryIncidence::new(
                CellRef::new(0, 2),
                CellRef::new(1, 6),
                IncidenceSign::Positive,
            ),
            BoundaryIncidence::new(
                CellRef::new(0, 4),
                CellRef::new(1, 7),
                IncidenceSign::Negative,
            ),
            BoundaryIncidence::new(
                CellRef::new(0, 3),
                CellRef::new(1, 7),
                IncidenceSign::Positive,
            ),
            BoundaryIncidence::new(
                CellRef::new(1, 0),
                CellRef::new(2, 0),
                IncidenceSign::Positive,
            ),
            BoundaryIncidence::new(
                CellRef::new(1, 4),
                CellRef::new(2, 0),
                IncidenceSign::Positive,
            ),
            BoundaryIncidence::new(
                CellRef::new(1, 5),
                CellRef::new(2, 0),
                IncidenceSign::Negative,
            ),
            BoundaryIncidence::new(
                CellRef::new(1, 1),
                CellRef::new(2, 1),
                IncidenceSign::Positive,
            ),
            BoundaryIncidence::new(
                CellRef::new(1, 5),
                CellRef::new(2, 1),
                IncidenceSign::Positive,
            ),
            BoundaryIncidence::new(
                CellRef::new(1, 6),
                CellRef::new(2, 1),
                IncidenceSign::Negative,
            ),
            BoundaryIncidence::new(
                CellRef::new(1, 2),
                CellRef::new(2, 2),
                IncidenceSign::Positive,
            ),
            BoundaryIncidence::new(
                CellRef::new(1, 6),
                CellRef::new(2, 2),
                IncidenceSign::Positive,
            ),
            BoundaryIncidence::new(
                CellRef::new(1, 7),
                CellRef::new(2, 2),
                IncidenceSign::Negative,
            ),
            BoundaryIncidence::new(
                CellRef::new(1, 3),
                CellRef::new(2, 3),
                IncidenceSign::Positive,
            ),
            BoundaryIncidence::new(
                CellRef::new(1, 7),
                CellRef::new(2, 3),
                IncidenceSign::Positive,
            ),
            BoundaryIncidence::new(
                CellRef::new(1, 4),
                CellRef::new(2, 3),
                IncidenceSign::Negative,
            ),
        ],
    )
    .expect("centered cellular surface");
    let all_cells = (0..5)
        .map(|ordinal| CellRef::new(0, ordinal))
        .chain((0..8).map(|ordinal| CellRef::new(1, ordinal)))
        .chain((0..4).map(|ordinal| CellRef::new(2, ordinal)))
        .collect::<Vec<_>>();
    let conductor = subcomplex(&complex, "support/surface-conductor", all_cells.clone());
    let component = ConductorComponent::new(
        ConductorComponentId::new("component/surface").expect("component id"),
        subcomplex(&complex, "support/component-surface", all_cells),
    )
    .expect("surface component");
    TerminalRelativePair::try_new(
        complex.clone(),
        conductor,
        subcomplex(
            &complex,
            "support/surface-relative",
            [
                CellRef::new(0, 0),
                CellRef::new(0, 1),
                CellRef::new(0, 2),
                CellRef::new(0, 3),
                CellRef::new(1, 0),
                CellRef::new(1, 2),
            ],
        ),
        subcomplex(&complex, "support/surface-insulation-empty", []),
        vec![component],
        vec![
            surface_terminal(
                &complex,
                0,
                [0, 1],
                "surface-driven",
                TerminalRole::Driven,
                TerminalOrientation::OutOfConductor,
                OrientationMapSign::Preserve,
            ),
            surface_terminal(
                &complex,
                2,
                [2, 3],
                "surface-return",
                TerminalRole::ReturnReference,
                TerminalOrientation::IntoConductor,
                OrientationMapSign::Reverse,
            ),
        ],
    )
    .expect("centered terminal-relative surface")
}

fn boundary_budget(
    max_rows: usize,
    max_cols: usize,
    max_entries: usize,
    max_retained: usize,
    max_component_visits: usize,
    max_incidence_visits: usize,
) -> TerminalRelativeBoundaryBudget {
    TerminalRelativeBoundaryBudget::new(
        max_rows,
        max_cols,
        max_entries,
        max_retained,
        max_component_visits,
        max_incidence_visits,
    )
}

fn kernel_budget(
    max_extent: usize,
    max_output_entries: usize,
    max_retained_entries: usize,
    max_binding_items: usize,
    max_scalar_operations: u128,
) -> KernelCoordinateBudget {
    KernelCoordinateBudget::new(
        max_extent,
        max_output_entries,
        max_retained_entries,
        max_binding_items,
        max_scalar_operations,
    )
}

fn construction_budget(
    max_live_entries: usize,
    max_elementary_operations: u128,
    max_entry_steps: u128,
    max_abs_coefficient: u128,
) -> SmithConstructionBudget {
    SmithConstructionBudget::new(
        max_live_entries,
        max_elementary_operations,
        max_entry_steps,
        max_abs_coefficient,
    )
}

fn centered_kernel_inputs() -> (
    TerminalRelativeBoundaryMatrix,
    TerminalRelativeBoundaryMatrix,
    VerifiedSmithNormalForm,
) {
    let pair = centered_surface_pair();
    let phase = PhaseId::new("phase/surface").expect("phase id");
    let outgoing = extract_terminal_relative_boundary_matrix(
        &pair,
        &phase,
        1,
        TerminalRelativeBoundaryBudget::default(),
    )
    .expect("surface outgoing boundary");
    let incoming = extract_terminal_relative_boundary_matrix(
        &pair,
        &phase,
        2,
        TerminalRelativeBoundaryBudget::default(),
    )
    .expect("surface incoming boundary");
    assert_eq!(outgoing.matrix().entries(), &[0, 0, -1, -1, -1, -1]);

    let diagonal = matrix(1, 6, &[1, 0, 0, 0, 0, 0]);
    let left = identity(1);
    let right = matrix(
        6,
        6,
        &[
            0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -1, 0, 0, -1, -1, -1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
            1, 0, 0, 0, 0, 0, 0, 1,
        ],
    );
    let right_inverse = matrix(
        6,
        6,
        &[
            0, 0, -1, -1, -1, -1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
            1, 0, 0, 0, 0, 0, 0, 1,
        ],
    );
    let smith = verify_smith_normal_form(
        outgoing.matrix().clone(),
        SmithNormalFormWitness::new(diagonal, left.clone(), left, right, right_inverse),
        budget(6, 476),
    )
    .expect("surface outgoing Smith witness");
    (outgoing, incoming, smith)
}

#[test]
fn it_001_complete_integer_witness_verifies_exactly() {
    let (source, witness) = rank_one_witness();
    let verified = verify_smith_normal_form(source.clone(), witness, budget(2, 48))
        .expect("complete exact witness");

    assert_eq!(verified.source(), &source);
    assert_eq!(verified.diagonal().entries(), &[2, 0, 0, 0]);
    assert_eq!(verified.invariant_factors(), &[2]);
    assert_eq!(verified.rank(), 1);
    assert_eq!(verified.scalar_operations(), 48);
    assert_eq!(
        verified.applicability(),
        TopologyApplicability::AbstractAlgebraOnly
    );
    assert_eq!(verified.left_transform().entries(), &[1, 0, -2, 1]);
    assert_eq!(verified.left_inverse().entries(), &[1, 0, 2, 1]);
    assert_eq!(verified.right_transform().entries(), &[1, -2, 0, 1]);
    assert_eq!(verified.right_inverse().entries(), &[1, 2, 0, 1]);
}

#[test]
fn it_002_canonical_diagonal_refuses_every_false_claim() {
    let identity3 = identity(3);
    for (diagonal, expected) in [
        (matrix(3, 3, &[2, 1, 0, 0, 4, 0, 0, 0, 0]), "off-diagonal"),
        (matrix(3, 3, &[-2, 0, 0, 0, 4, 0, 0, 0, 0]), "negative"),
        (matrix(3, 3, &[2, 0, 0, 0, 0, 0, 0, 0, 4]), "zero ordering"),
        (matrix(3, 3, &[4, 0, 0, 0, 6, 0, 0, 0, 0]), "divisibility"),
    ] {
        let witness = SmithNormalFormWitness::new(
            diagonal.clone(),
            identity3.clone(),
            identity3.clone(),
            identity3.clone(),
            identity3.clone(),
        );
        let error = verify_smith_normal_form(diagonal, witness, budget(3, 162))
            .expect_err("noncanonical diagonal must refuse");
        match expected {
            "off-diagonal" => assert!(matches!(
                error,
                IntegralTopologyError::OffDiagonalEntry { .. }
            )),
            "negative" => assert!(matches!(
                error,
                IntegralTopologyError::NegativeInvariantFactor { .. }
            )),
            "zero ordering" => assert!(matches!(
                error,
                IntegralTopologyError::NonzeroAfterZero { .. }
            )),
            "divisibility" => assert!(matches!(
                error,
                IntegralTopologyError::InvariantFactorDivisibility { .. }
            )),
            _ => unreachable!(),
        }
    }
}

#[test]
fn it_003_integer_inverse_is_required_even_when_zero_source_hides_transform() {
    let source = matrix(1, 1, &[0]);
    let witness = SmithNormalFormWitness::new(
        matrix(1, 1, &[0]),
        matrix(1, 1, &[2]),
        matrix(1, 1, &[1]),
        matrix(1, 1, &[1]),
        matrix(1, 1, &[1]),
    );
    let error = verify_smith_normal_form(source, witness, budget(1, 6))
        .expect_err("non-unimodular transform must be refuted");
    assert_eq!(error.failure_class(), IntegralTopologyFailureClass::Refuted);
    assert!(matches!(
        error,
        IntegralTopologyError::WitnessProductMismatch {
            stage: SmithWitnessStage::LeftTimesInverse,
            expected: 1,
            actual: 2,
            ..
        }
    ));
}

#[test]
fn it_004_transform_mutation_cannot_rebind_the_source() {
    let (source, witness) = rank_one_witness();
    let corrupted = SmithNormalFormWitness::new(
        matrix(2, 2, &[2, 0, 0, 2]),
        witness.left().clone(),
        witness.left_inverse().clone(),
        witness.right().clone(),
        witness.right_inverse().clone(),
    );
    assert!(matches!(
        verify_smith_normal_form(source, corrupted, budget(2, 48)),
        Err(IntegralTopologyError::WitnessProductMismatch {
            stage: SmithWitnessStage::DiagonalTransform,
            row: 1,
            col: 1,
            expected: 2,
            actual: 0,
        })
    ));
}

#[test]
fn it_005_shape_and_retained_storage_are_preflighted() {
    let (source, witness) = rank_one_witness();
    let wrong_left = SmithNormalFormWitness::new(
        witness.diagonal().clone(),
        matrix(1, 1, &[1]),
        witness.left_inverse().clone(),
        witness.right().clone(),
        witness.right_inverse().clone(),
    );
    assert!(matches!(
        verify_smith_normal_form(source.clone(), wrong_left, budget(2, 48)),
        Err(IntegralTopologyError::WitnessShape {
            role: MatrixRole::LeftTransform,
            expected_rows: 2,
            expected_cols: 2,
            actual_rows: 1,
            actual_cols: 1,
        })
    ));

    let tight_retained = ExactAlgebraBudget::new(2, 2, 4, 23, 4, 48);
    assert!(matches!(
        verify_smith_normal_form(source, witness, tight_retained),
        Err(IntegralTopologyError::RetainedEntryBudgetExceeded {
            requested: 24,
            max: 23,
        })
    ));
}

#[test]
fn it_006_exact_work_cap_and_limit_plus_one_refuse_before_execution() {
    let (source, witness) = rank_one_witness();
    let verified = verify_smith_normal_form(source.clone(), witness.clone(), budget(2, 48))
        .expect("exact scalar-work cap admits");
    assert_eq!(verified.scalar_operations(), 48);
    assert!(matches!(
        verify_smith_normal_form(source, witness, budget(2, 47)),
        Err(IntegralTopologyError::ScalarWorkBudgetExceeded {
            requested: 48,
            max: 47,
        })
    ));
}

#[test]
fn it_007_every_cancellation_poll_is_transactional() {
    let (source, witness) = rank_one_witness();
    let mut poll_count = 0_usize;
    let verified = verify_smith_normal_form_with_checkpoint(
        source.clone(),
        witness.clone(),
        budget(2, 48),
        &mut |_| {
            poll_count += 1;
            true
        },
    )
    .expect("uninterrupted witness");
    assert_eq!(verified.scalar_operations(), 48);

    for stop_at in 0..poll_count {
        let mut observed = 0_usize;
        let result = verify_smith_normal_form_with_checkpoint(
            source.clone(),
            witness.clone(),
            budget(2, 48),
            &mut |_| {
                let keep_running = observed != stop_at;
                observed += 1;
                keep_running
            },
        );
        let error = result.expect_err("every injected cancellation must refuse");
        assert!(
            matches!(&error, IntegralTopologyError::Cancelled { .. }),
            "poll {stop_at} must publish only cancellation: {error:?}"
        );
        if stop_at + 1 == poll_count {
            assert!(matches!(
                &error,
                IntegralTopologyError::Cancelled {
                    phase: "smith witness finalize",
                    completed_scalar_operations: 48,
                    planned_scalar_operations: 48,
                }
            ));
        }
    }
}

#[test]
fn it_008_checked_arithmetic_refuses_overflow() {
    let shear = i128::MAX;
    let source = matrix(2, 2, &[0, 0, 2, 0]);
    let witness = SmithNormalFormWitness::new(
        matrix(2, 2, &[0, 0, 0, 0]),
        matrix(2, 2, &[1, shear, 0, 1]),
        matrix(2, 2, &[1, -shear, 0, 1]),
        identity(2),
        identity(2),
    );
    let error = verify_smith_normal_form(source, witness, budget(2, 48))
        .expect_err("coefficient explosion must remain unknown");
    assert_eq!(error.failure_class(), IntegralTopologyFailureClass::Unknown);
    assert!(matches!(
        error,
        IntegralTopologyError::ArithmeticOverflow {
            stage: SmithWitnessStage::LeftTimesSource,
            row: 0,
            col: 0,
            term: 1,
        }
    ));
}

#[test]
fn it_009_empty_exact_complex_is_not_forced_nonempty() {
    let empty = matrix(0, 0, &[]);
    let witness = SmithNormalFormWitness::new(
        empty.clone(),
        empty.clone(),
        empty.clone(),
        empty.clone(),
        empty.clone(),
    );
    let verified =
        verify_smith_normal_form(empty, witness, ExactAlgebraBudget::new(0, 0, 0, 0, 0, 0))
            .expect("zero-dimensional exact algebra is valid");
    assert_eq!(verified.rank(), 0);
    assert!(verified.invariant_factors().is_empty());
    assert_eq!(verified.scalar_operations(), 0);
}

#[test]
fn it_010_matrix_extent_entry_count_and_entry_budget_refuse_exactly() {
    let small = ExactAlgebraBudget::new(2, 2, 4, 24, 4, 48);
    assert!(matches!(
        ExactIntegerMatrix::try_new(3, 1, vec![0; 3], small),
        Err(IntegralTopologyError::MatrixExtentExceeded { .. })
    ));
    assert!(matches!(
        ExactIntegerMatrix::try_new(2, 2, vec![0; 3], small),
        Err(IntegralTopologyError::MatrixEntryCount {
            expected: 4,
            actual: 3,
            ..
        })
    ));
    let entry_tight = ExactAlgebraBudget::new(2, 2, 3, 24, 4, 48);
    assert!(matches!(
        ExactIntegerMatrix::try_new(2, 2, vec![0; 4], entry_tight),
        Err(IntegralTopologyError::MatrixEntryBudgetExceeded {
            requested: 4,
            max: 3,
        })
    ));
}

#[test]
fn it_011_rectangular_and_empty_shapes_remain_exactly_distinct() {
    let rectangular = matrix(2, 3, &[2, 0, 0, 0, 6, 0]);
    let witness = SmithNormalFormWitness::new(
        rectangular.clone(),
        identity(2),
        identity(2),
        identity(3),
        identity(3),
    );
    let verified = verify_smith_normal_form(rectangular, witness, budget(3, 100))
        .expect("rectangular canonical Smith matrix");
    assert_eq!(verified.source().rows(), 2);
    assert_eq!(verified.source().cols(), 3);
    assert_eq!(verified.invariant_factors(), &[2, 6]);
    assert_eq!(verified.scalar_operations(), 100);

    let zero_by_three = matrix(0, 3, &[]);
    let witness = SmithNormalFormWitness::new(
        zero_by_three.clone(),
        matrix(0, 0, &[]),
        matrix(0, 0, &[]),
        identity(3),
        identity(3),
    );
    let verified = verify_smith_normal_form(zero_by_three, witness, budget(3, 54))
        .expect("zero-by-three matrix keeps its shape");
    assert_eq!((verified.source().rows(), verified.source().cols()), (0, 3));
    assert_eq!(verified.scalar_operations(), 54);

    let three_by_zero = matrix(3, 0, &[]);
    let witness = SmithNormalFormWitness::new(
        three_by_zero.clone(),
        identity(3),
        identity(3),
        matrix(0, 0, &[]),
        matrix(0, 0, &[]),
    );
    let verified = verify_smith_normal_form(three_by_zero, witness, budget(3, 54))
        .expect("three-by-zero matrix keeps its shape");
    assert_eq!((verified.source().rows(), verified.source().cols()), (3, 0));
    assert_eq!(verified.scalar_operations(), 54);
}

#[test]
fn it_012_terminal_relative_boundary_binds_pair_phase_component_and_bases() {
    let pair = terminal_cut_loop_pair(false);
    let phase = PhaseId::new("phase/a").expect("phase id");
    let boundary = extract_terminal_relative_boundary_matrix(
        &pair,
        &phase,
        1,
        boundary_budget(2, 4, 8, 14, 16, 8),
    )
    .expect("canonical terminal-relative boundary");

    assert_eq!(boundary.pair_id(), pair.identity());
    assert_eq!(boundary.phase(), &phase);
    assert_eq!(boundary.component().as_str(), "component/winding");
    assert_eq!(boundary.source_degree(), 1);
    assert_eq!(boundary.target_degree(), Some(0));
    assert_eq!(
        boundary.source_basis(),
        &[
            CellRef::new(1, 0),
            CellRef::new(1, 1),
            CellRef::new(1, 2),
            CellRef::new(1, 3),
        ]
    );
    assert_eq!(
        boundary.target_basis(),
        &[CellRef::new(0, 1), CellRef::new(0, 2)]
    );
    assert_eq!((boundary.matrix().rows(), boundary.matrix().cols()), (2, 4));
    assert_eq!(boundary.matrix().entries(), &[1, -1, -1, 0, 0, 1, 1, -1]);
    assert_eq!(boundary.work_items(), 24);
    assert_eq!(
        boundary.applicability(),
        TopologyApplicability::TerminalRelativeIncidenceOnly
    );

    let chain = IntegralRelativeChain::try_new(&pair, phase.clone(), 1, vec![2, -1, 3, 4])
        .expect("relative chain");
    assert_eq!(
        pair.boundary(&chain)
            .expect("exact boundary")
            .coefficients(),
        &[0, -2]
    );
    let cochain =
        IntegralRelativeCochain::try_new(&pair, phase, 0, vec![2, 5]).expect("relative cochain");
    assert_eq!(
        pair.integral_coboundary(&cochain)
            .expect("exact coboundary")
            .coefficients(),
        &[2, 3, 3, -5]
    );
}

#[test]
fn it_013_pair_boundary_replays_across_declaration_order() {
    let forward = terminal_cut_loop_pair(false);
    let reverse = terminal_cut_loop_pair(true);
    let phase = PhaseId::new("phase/a").expect("phase id");
    assert_eq!(forward.identity(), reverse.identity());

    let forward = extract_terminal_relative_boundary_matrix(
        &forward,
        &phase,
        1,
        TerminalRelativeBoundaryBudget::default(),
    )
    .expect("forward boundary");
    let reverse = extract_terminal_relative_boundary_matrix(
        &reverse,
        &phase,
        1,
        TerminalRelativeBoundaryBudget::default(),
    )
    .expect("reverse boundary");
    assert_eq!(forward, reverse);
}

#[test]
fn it_014_unaugmented_edge_zero_matrices_preserve_rectangular_shapes() {
    let pair = terminal_cut_loop_pair(false);
    let phase = PhaseId::new("phase/a").expect("phase id");

    let bottom = extract_terminal_relative_boundary_matrix(
        &pair,
        &phase,
        0,
        TerminalRelativeBoundaryBudget::default(),
    )
    .expect("bottom zero map");
    assert_eq!(bottom.source_degree(), 0);
    assert_eq!(bottom.target_degree(), None);
    assert_eq!((bottom.matrix().rows(), bottom.matrix().cols()), (0, 2));
    assert!(bottom.target_basis().is_empty());
    assert_eq!(
        bottom.source_basis(),
        &[CellRef::new(0, 1), CellRef::new(0, 2)]
    );

    let top = extract_terminal_relative_boundary_matrix(
        &pair,
        &phase,
        2,
        TerminalRelativeBoundaryBudget::default(),
    )
    .expect("top zero map");
    assert_eq!(top.source_degree(), 2);
    assert_eq!(top.target_degree(), Some(1));
    assert_eq!((top.matrix().rows(), top.matrix().cols()), (4, 0));
    assert!(top.source_basis().is_empty());
    assert_eq!(
        top.target_basis(),
        &[
            CellRef::new(1, 0),
            CellRef::new(1, 1),
            CellRef::new(1, 2),
            CellRef::new(1, 3),
        ]
    );
    assert!(bottom.matrix().entries().is_empty());
    assert!(top.matrix().entries().is_empty());
}

#[test]
fn it_015_pair_boundary_preflights_every_independent_budget() {
    let pair = terminal_cut_loop_pair(false);
    let phase = PhaseId::new("phase/a").expect("phase id");
    for (budget, expected) in [
        (boundary_budget(1, 4, 8, 14, 16, 8), "rows"),
        (boundary_budget(2, 3, 8, 14, 16, 8), "columns"),
        (boundary_budget(2, 4, 7, 14, 16, 8), "entries"),
        (boundary_budget(2, 4, 8, 13, 16, 8), "retained"),
        (boundary_budget(2, 4, 8, 14, 15, 8), "component visits"),
        (boundary_budget(2, 4, 8, 14, 16, 7), "incidence visits"),
    ] {
        let error = extract_terminal_relative_boundary_matrix(&pair, &phase, 1, budget)
            .expect_err("limit minus one must refuse");
        match expected {
            "rows" | "columns" => assert!(matches!(
                &error,
                IntegralTopologyError::MatrixExtentExceeded {
                    rows: 2,
                    cols: 4,
                    ..
                }
            )),
            "entries" => assert!(matches!(
                &error,
                IntegralTopologyError::MatrixEntryBudgetExceeded {
                    requested: 8,
                    max: 7,
                }
            )),
            "retained" => assert!(matches!(
                &error,
                IntegralTopologyError::RetainedEntryBudgetExceeded {
                    requested: 14,
                    max: 13,
                }
            )),
            "component visits" => assert!(matches!(
                &error,
                IntegralTopologyError::ComponentVisitBudgetExceeded {
                    requested: 16,
                    max: 15,
                }
            )),
            "incidence visits" => assert!(matches!(
                &error,
                IntegralTopologyError::IncidenceVisitBudgetExceeded {
                    requested: 8,
                    max: 7,
                }
            )),
            _ => unreachable!(),
        }
        assert_eq!(error.failure_class(), IntegralTopologyFailureClass::Unknown);
    }
}

#[test]
fn it_016_pair_boundary_cancellation_is_transactional_through_publication() {
    let pair = terminal_cut_loop_pair(false);
    let phase = PhaseId::new("phase/a").expect("phase id");
    let mut poll_count = 0_usize;
    let boundary = extract_terminal_relative_boundary_matrix_with_checkpoint(
        &pair,
        &phase,
        1,
        TerminalRelativeBoundaryBudget::default(),
        &mut |_| {
            poll_count += 1;
            true
        },
    )
    .expect("uninterrupted extraction");
    assert_eq!(boundary.work_items(), 24);

    for stop_at in 0..poll_count {
        let mut observed = 0_usize;
        let error = extract_terminal_relative_boundary_matrix_with_checkpoint(
            &pair,
            &phase,
            1,
            TerminalRelativeBoundaryBudget::default(),
            &mut |_| {
                let keep_running = observed != stop_at;
                observed += 1;
                keep_running
            },
        )
        .expect_err("injected cancellation must refuse");
        assert_eq!(error.failure_class(), IntegralTopologyFailureClass::Unknown);
        assert!(matches!(
            &error,
            IntegralTopologyError::PairBoundaryCancelled { .. }
        ));
        if stop_at + 1 == poll_count {
            assert!(matches!(
                &error,
                IntegralTopologyError::PairBoundaryCancelled {
                    phase: "terminal-relative boundary finalize",
                    completed_work_items: 24,
                    planned_work_items: 24,
                }
            ));
        }
    }
}

#[test]
fn it_017_pair_boundary_refuses_unknown_phase_and_excess_degree() {
    let pair = terminal_cut_loop_pair(false);
    let unknown = PhaseId::new("phase/unknown").expect("phase id");
    let error = extract_terminal_relative_boundary_matrix(
        &pair,
        &unknown,
        1,
        TerminalRelativeBoundaryBudget::default(),
    )
    .expect_err("unknown phase must refuse");
    assert_eq!(error.failure_class(), IntegralTopologyFailureClass::Refuted);
    assert!(matches!(
        error,
        IntegralTopologyError::UnknownTerminalRelativePhase { .. }
    ));

    let phase = PhaseId::new("phase/a").expect("phase id");
    let error = extract_terminal_relative_boundary_matrix(
        &pair,
        &phase,
        3,
        TerminalRelativeBoundaryBudget::default(),
    )
    .expect_err("degree above dimension plus one must refuse");
    assert_eq!(error.failure_class(), IntegralTopologyFailureClass::Refuted);
    assert!(matches!(
        error,
        IntegralTopologyError::BoundaryDegreeOutOfRange { degree: 3, max: 2 }
    ));
}

#[test]
fn it_018_adjacent_boundary_image_is_verified_in_smith_kernel_coordinates() {
    let (outgoing, incoming, smith) = centered_kernel_inputs();
    let transport = verify_terminal_relative_kernel_transport(
        outgoing,
        incoming,
        smith,
        kernel_budget(6, 20, 154, 12, 144),
    )
    .expect("exact adjacent-boundary transport");

    assert_eq!(transport.degree(), 1);
    assert_eq!(transport.phase().as_str(), "phase/surface");
    assert_eq!(transport.component().as_str(), "component/surface");
    assert_eq!(transport.chain_extent(), 6);
    assert_eq!(transport.outgoing_rank(), 1);
    assert_eq!(transport.kernel_dimension(), 5);
    assert_eq!(transport.binding_items(), 12);
    assert_eq!(transport.scalar_operations(), 144);
    assert_eq!(
        (
            transport.kernel_image().rows(),
            transport.kernel_image().cols()
        ),
        (5, 4)
    );
    assert_eq!(
        transport.kernel_image().entries(),
        &[
            0, 1, 0, 0, 0, 0, 0, 1, -1, 1, 0, 0, 0, -1, 1, 0, 0, 0, -1, 1
        ]
    );
    assert_eq!(
        transport.applicability(),
        TopologyApplicability::TerminalRelativeKernelCoordinatesOnly
    );
}

#[test]
fn it_019_bottom_edge_uses_inverse_smith_transform_not_cell_labels() {
    let pair = terminal_cut_loop_pair(false);
    let phase = PhaseId::new("phase/a").expect("phase id");
    let outgoing = extract_terminal_relative_boundary_matrix(
        &pair,
        &phase,
        0,
        TerminalRelativeBoundaryBudget::default(),
    )
    .expect("bottom outgoing boundary");
    let incoming = extract_terminal_relative_boundary_matrix(
        &pair,
        &phase,
        1,
        TerminalRelativeBoundaryBudget::default(),
    )
    .expect("incoming boundary");
    let empty = matrix(0, 0, &[]);
    let smith = verify_smith_normal_form(
        outgoing.matrix().clone(),
        SmithNormalFormWitness::new(
            matrix(0, 2, &[]),
            empty.clone(),
            empty,
            matrix(2, 2, &[1, 1, 0, 1]),
            matrix(2, 2, &[1, -1, 0, 1]),
        ),
        budget(2, 16),
    )
    .expect("bottom zero-map Smith witness");
    let transport = verify_terminal_relative_kernel_transport(
        outgoing,
        incoming,
        smith,
        kernel_budget(4, 8, 32, 2, 16),
    )
    .expect("bottom incoming transport");

    assert_eq!(transport.outgoing_rank(), 0);
    assert_eq!(transport.kernel_dimension(), 2);
    assert_eq!(
        transport.kernel_image().entries(),
        &[1, -2, -2, 1, 0, 1, 1, -1]
    );
    assert_ne!(
        transport.kernel_image(),
        transport.incoming_boundary().matrix(),
        "kernel rows are witness coordinates, not retained CellRef labels"
    );
}

#[test]
fn it_020_kernel_transport_refuses_source_and_every_budget_mismatch() {
    let (outgoing, incoming, _) = centered_kernel_inputs();
    let zero = matrix(1, 6, &[0, 0, 0, 0, 0, 0]);
    let identity6 = identity(6);
    let wrong_smith = verify_smith_normal_form(
        zero.clone(),
        SmithNormalFormWitness::new(zero, identity(1), identity(1), identity6.clone(), identity6),
        budget(6, 476),
    )
    .expect("valid Smith witness for the wrong source");
    let error = verify_terminal_relative_kernel_transport(
        outgoing.clone(),
        incoming.clone(),
        wrong_smith,
        kernel_budget(6, 20, 154, 12, 144),
    )
    .expect_err("same-shape source mutation must refuse");
    assert_eq!(error.failure_class(), IntegralTopologyFailureClass::Refuted);
    assert!(matches!(
        error,
        IntegralTopologyError::OutgoingSmithSourceEntryMismatch {
            row: 0,
            col: 2,
            expected: -1,
            actual: 0,
        }
    ));

    let (_, _, smith) = centered_kernel_inputs();
    for (budget, expected) in [
        (kernel_budget(5, 20, 154, 12, 144), "extent"),
        (kernel_budget(6, 19, 154, 12, 144), "output"),
        (kernel_budget(6, 20, 153, 12, 144), "retained"),
        (kernel_budget(6, 20, 154, 11, 144), "binding"),
        (kernel_budget(6, 20, 154, 12, 143), "scalar"),
    ] {
        let error = verify_terminal_relative_kernel_transport(
            outgoing.clone(),
            incoming.clone(),
            smith.clone(),
            budget,
        )
        .expect_err("limit minus one must refuse");
        match expected {
            "extent" => assert!(matches!(
                &error,
                IntegralTopologyError::KernelCoordinateExtentExceeded {
                    chain_extent: 6,
                    max: 5,
                    ..
                }
            )),
            "output" => assert!(matches!(
                &error,
                IntegralTopologyError::MatrixEntryBudgetExceeded {
                    requested: 20,
                    max: 19,
                }
            )),
            "retained" => assert!(matches!(
                &error,
                IntegralTopologyError::RetainedEntryBudgetExceeded {
                    requested: 154,
                    max: 153,
                }
            )),
            "binding" => assert!(matches!(
                &error,
                IntegralTopologyError::KernelBindingBudgetExceeded {
                    requested: 12,
                    max: 11,
                }
            )),
            "scalar" => assert!(matches!(
                &error,
                IntegralTopologyError::ScalarWorkBudgetExceeded {
                    requested: 144,
                    max: 143,
                }
            )),
            _ => unreachable!(),
        }
        assert_eq!(error.failure_class(), IntegralTopologyFailureClass::Unknown);
    }
}

#[test]
fn it_021_kernel_transport_cancellation_is_transactional_through_publication() {
    let (outgoing, incoming, smith) = centered_kernel_inputs();
    let mut poll_count = 0_usize;
    let transport = verify_terminal_relative_kernel_transport_with_checkpoint(
        outgoing.clone(),
        incoming.clone(),
        smith.clone(),
        kernel_budget(6, 20, 154, 12, 144),
        &mut |_| {
            poll_count += 1;
            true
        },
    )
    .expect("uninterrupted kernel transport");
    assert_eq!(transport.binding_items(), 12);
    assert_eq!(transport.scalar_operations(), 144);

    for stop_at in 0..poll_count {
        let mut observed = 0_usize;
        let error = verify_terminal_relative_kernel_transport_with_checkpoint(
            outgoing.clone(),
            incoming.clone(),
            smith.clone(),
            kernel_budget(6, 20, 154, 12, 144),
            &mut |_| {
                let keep_running = observed != stop_at;
                observed += 1;
                keep_running
            },
        )
        .expect_err("injected cancellation must refuse");
        assert_eq!(error.failure_class(), IntegralTopologyFailureClass::Unknown);
        assert!(matches!(
            &error,
            IntegralTopologyError::KernelCoordinateCancelled { .. }
        ));
        if stop_at + 1 == poll_count {
            assert!(matches!(
                &error,
                IntegralTopologyError::KernelCoordinateCancelled {
                    phase: "kernel-coordinate finalize",
                    completed_binding_items: 12,
                    planned_binding_items: 12,
                    completed_scalar_operations: 144,
                    planned_scalar_operations: 144,
                }
            ));
        }
    }
}

#[test]
fn it_022_constructive_smith_repairs_divisibility_and_replays_byte_exactly() {
    let cases = [
        (matrix(2, 2, &[2, 4, 4, 8]), vec![2, 0, 0, 0], vec![2]),
        (matrix(2, 2, &[4, 0, 0, 6]), vec![2, 0, 0, 12], vec![2, 12]),
        (
            matrix(3, 3, &[6, 0, 0, 0, 10, 0, 0, 0, 15]),
            vec![1, 0, 0, 0, 30, 0, 0, 0, 30],
            vec![1, 30, 30],
        ),
        (
            matrix(2, 3, &[2, 4, 0, 0, 6, 3]),
            vec![1, 0, 0, 0, 6, 0],
            vec![1, 6],
        ),
    ];
    for (source, expected_diagonal, expected_invariants) in cases {
        let constructed = construct_smith_normal_form(
            source.clone(),
            SmithConstructionBudget::default(),
            ExactAlgebraBudget::default(),
        )
        .expect("deterministic constructor must publish only a verified witness");
        assert_eq!(constructed.verified().source(), &source);
        assert_eq!(
            constructed.verified().diagonal().entries(),
            expected_diagonal
        );
        assert_eq!(
            constructed.verified().invariant_factors(),
            expected_invariants
        );
        assert!(
            constructed.entry_steps()
                >= u128::try_from(source.entries().len()).expect("fixture entry count")
        );
        assert_eq!(
            constructed.applicability(),
            TopologyApplicability::AbstractAlgebraOnly
        );
    }

    let source = matrix(2, 2, &[4, 0, 0, 6]);
    let first = construct_smith_normal_form(
        source.clone(),
        SmithConstructionBudget::default(),
        ExactAlgebraBudget::default(),
    )
    .expect("first deterministic replay");
    let second = construct_smith_normal_form(
        source,
        SmithConstructionBudget::default(),
        ExactAlgebraBudget::default(),
    )
    .expect("second deterministic replay");
    assert_eq!(first, second, "all transforms, inverses, and counts replay");
}

#[test]
fn it_023_constructive_smith_preserves_rectangles_signs_and_empty_shapes() {
    for (source, expected_diagonal, expected_invariants) in [
        (matrix(1, 1, &[-7]), vec![7], vec![7]),
        (matrix(1, 3, &[6, -10, 15]), vec![1, 0, 0], vec![1]),
        (matrix(3, 1, &[6, -10, 15]), vec![1, 0, 0], vec![1]),
        (
            matrix(2, 3, &[0, 0, 0, 0, 0, 0]),
            vec![0, 0, 0, 0, 0, 0],
            vec![],
        ),
        (matrix(0, 0, &[]), vec![], vec![]),
        (matrix(0, 3, &[]), vec![], vec![]),
        (matrix(3, 0, &[]), vec![], vec![]),
    ] {
        let shape = (source.rows(), source.cols());
        let constructed = construct_smith_normal_form(
            source,
            SmithConstructionBudget::default(),
            ExactAlgebraBudget::default(),
        )
        .expect("rectangular or empty exact source");
        assert_eq!(
            (
                constructed.verified().diagonal().rows(),
                constructed.verified().diagonal().cols()
            ),
            shape
        );
        assert_eq!(
            constructed.verified().diagonal().entries(),
            expected_diagonal
        );
        assert_eq!(
            constructed.verified().invariant_factors(),
            expected_invariants
        );
    }
}

#[test]
fn it_024_constructive_smith_enforces_every_data_dependent_budget() {
    let source = matrix(2, 2, &[4, 0, 0, 6]);
    let baseline = construct_smith_normal_form(
        source.clone(),
        SmithConstructionBudget::default(),
        budget(2, 48),
    )
    .expect("baseline construction");
    let operations = baseline.elementary_operations();
    let steps = baseline.entry_steps();
    assert_eq!(operations, 6);
    assert_eq!(steps, 71);

    construct_smith_normal_form(
        source.clone(),
        construction_budget(24, operations, steps, i128::MAX.unsigned_abs()),
        budget(2, 48),
    )
    .expect("every exact construction limit admits");

    let error = construct_smith_normal_form(
        source.clone(),
        construction_budget(24, operations, 3, i128::MAX.unsigned_abs()),
        budget(2, 48),
    )
    .expect_err("source coefficient scan itself is entry-budgeted");
    assert!(matches!(
        error,
        IntegralTopologyError::SmithConstructionEntryStepBudgetExceeded {
            requested: 4,
            max: 3,
        }
    ));

    for (construction_limits, expected) in [
        (
            construction_budget(23, operations, steps, i128::MAX.unsigned_abs()),
            "live",
        ),
        (
            construction_budget(24, operations - 1, steps, i128::MAX.unsigned_abs()),
            "operations",
        ),
        (
            construction_budget(24, operations, steps - 1, i128::MAX.unsigned_abs()),
            "steps",
        ),
        (construction_budget(24, operations, steps, 5), "coefficient"),
    ] {
        let error = construct_smith_normal_form(source.clone(), construction_limits, budget(2, 48))
            .expect_err("limit minus one must refuse construction");
        assert_eq!(error.failure_class(), IntegralTopologyFailureClass::Unknown);
        match expected {
            "live" => assert!(matches!(
                error,
                IntegralTopologyError::SmithConstructionLiveEntryBudgetExceeded {
                    requested: 24,
                    max: 23,
                }
            )),
            "operations" => assert!(matches!(
                error,
                IntegralTopologyError::SmithConstructionOperationBudgetExceeded { .. }
            )),
            "steps" => assert!(matches!(
                error,
                IntegralTopologyError::SmithConstructionEntryStepBudgetExceeded { .. }
            )),
            "coefficient" => assert!(matches!(
                error,
                IntegralTopologyError::SmithConstructionCoefficientMagnitudeExceeded {
                    role: MatrixRole::Source,
                    magnitude: 6,
                    max: 5,
                    ..
                }
            )),
            _ => unreachable!(),
        }
    }

    let error = construct_smith_normal_form(
        matrix(1, 1, &[i128::MIN]),
        SmithConstructionBudget::default(),
        ExactAlgebraBudget::default(),
    )
    .expect_err("unrepresentable positive magnitude must remain unknown");
    assert_eq!(error.failure_class(), IntegralTopologyFailureClass::Unknown);
    assert!(matches!(
        error,
        IntegralTopologyError::SmithConstructionCoefficientMagnitudeExceeded {
            role: MatrixRole::Source,
            row: 0,
            col: 0,
            magnitude,
            max,
        } if magnitude == i128::MIN.unsigned_abs() && max == i128::MAX.unsigned_abs()
    ));

    let error = construct_smith_normal_form(
        matrix(3, 3, &[6, 0, 0, 0, 10, 0, 0, 0, 15]),
        construction_budget(54, 1_000, 10_000, 29),
        budget(3, 162),
    )
    .expect_err("generated coefficient growth must obey its cap");
    assert_eq!(error.failure_class(), IntegralTopologyFailureClass::Unknown);
    assert!(matches!(
        error,
        IntegralTopologyError::SmithConstructionCoefficientMagnitudeExceeded {
            role,
            magnitude,
            max: 29,
            ..
        } if role != MatrixRole::Source && magnitude > 29
    ));
}

#[test]
fn it_025_constructive_smith_cancellation_covers_reduction_and_verification() {
    let source = matrix(2, 2, &[4, 0, 0, 6]);
    let mut poll_count = 0_usize;
    let complete = construct_smith_normal_form_with_checkpoint(
        source.clone(),
        SmithConstructionBudget::default(),
        budget(2, 48),
        &mut |_| {
            poll_count += 1;
            true
        },
    )
    .expect("uninterrupted construction and verification");
    assert_eq!(complete.verified().invariant_factors(), &[2, 12]);

    for stop_at in 0..poll_count {
        let mut observed = 0_usize;
        let error = construct_smith_normal_form_with_checkpoint(
            source.clone(),
            SmithConstructionBudget::default(),
            budget(2, 48),
            &mut |_| {
                let keep_running = observed != stop_at;
                observed += 1;
                keep_running
            },
        )
        .expect_err("every cancellation poll must refuse publication");
        assert_eq!(error.failure_class(), IntegralTopologyFailureClass::Unknown);
        assert!(matches!(
            &error,
            IntegralTopologyError::SmithConstructionCancelled { .. }
                | IntegralTopologyError::Cancelled { .. }
        ));
        if stop_at + 1 == poll_count {
            assert!(matches!(
                error,
                IntegralTopologyError::Cancelled {
                    phase: "smith witness finalize",
                    completed_scalar_operations: 48,
                    planned_scalar_operations: 48,
                }
            ));
        }
    }
}

fn gcd_u128(mut left: u128, mut right: u128) -> u128 {
    while right != 0 {
        let remainder = left % right;
        left = right;
        right = remainder;
    }
    left
}

#[test]
fn it_026_constructive_smith_exhausts_small_two_by_two_oracle() {
    for a in -2_i128..=2 {
        for b in -2_i128..=2 {
            for c in -2_i128..=2 {
                for d in -2_i128..=2 {
                    let source = matrix(2, 2, &[a, b, c, d]);
                    let constructed = construct_smith_normal_form(
                        source,
                        SmithConstructionBudget::default(),
                        budget(2, 48),
                    )
                    .expect("small exact matrix must remain inside i128");
                    let gcd = [a, b, c, d]
                        .into_iter()
                        .map(i128::unsigned_abs)
                        .fold(0_u128, gcd_u128);
                    let determinant = a * d - b * c;
                    let expected = if gcd == 0 {
                        Vec::new()
                    } else if determinant == 0 {
                        vec![i128::try_from(gcd).expect("small gcd")]
                    } else {
                        vec![
                            i128::try_from(gcd).expect("small gcd"),
                            i128::try_from(determinant.unsigned_abs() / gcd)
                                .expect("small determinant quotient"),
                        ]
                    };
                    assert_eq!(
                        constructed.verified().invariant_factors(),
                        expected,
                        "source [{a}, {b}; {c}, {d}]"
                    );
                }
            }
        }
    }
}
