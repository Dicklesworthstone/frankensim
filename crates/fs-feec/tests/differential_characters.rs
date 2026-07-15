//! G0/G3 schema battery for the RA.2a finite relative differential-character
//! conventions.  These tests exercise typed domains/codomains, coefficient
//! and degree refusals, relative/terminal semantics, and deterministic replay.
//! They do not claim constructive exactness or smooth convergence.

use core::num::NonZeroU16;
use fs_feec::differential_characters::{
    AlgebraBudget, AlgebraId, BilinearMapKind, BoundaryComponent, BoundaryOrientation,
    BoundaryRole, CharacterError, CoefficientLattice, CoefficientProductSchema, CoefficientSector,
    CoefficientSystem, CohomologicalDegree, ComplexLane, ExactSequenceKind, ExactnessStatus,
    FiniteComplexSchema, FiniteExactnessAssumption, MapKind, NilpotenceLaw, ObjectKind,
    ObjectSupport, Orientation, RelativeDifferentialCharacter, RelativeModel, RelativePairSchema,
    RelativeTrivialization, RepresentativeKind, TerminalTrivialization,
};
use fs_qty::Dims;

fn budget() -> AlgebraBudget {
    AlgebraBudget::new(100_000, 16, 64 * 1024).expect("fixture budget")
}

fn flux_lattice(scale: f64) -> CoefficientLattice {
    CoefficientLattice::new(
        "magnetic-flux-lattice",
        1,
        Dims([2, 1, -2, 0, -1, 0]),
        &[scale],
    )
    .expect("fixture lattice")
}

fn flux_product() -> CoefficientProductSchema {
    let coefficients = CoefficientSystem::RealModuloLattice(flux_lattice(1.0));
    CoefficientProductSchema::new(
        "flux-cup-product",
        1,
        coefficients.clone(),
        coefficients.clone(),
        coefficients,
        AlgebraId::from_bytes([0x42; 32]),
        budget(),
    )
    .expect("fixture coefficient product")
}

fn empty_relative(id: &str, dimension: u8, lane: ComplexLane) -> FiniteComplexSchema {
    FiniteComplexSchema::new(
        id,
        1,
        dimension,
        lane,
        Orientation::Positive,
        vec![0; usize::from(dimension) + 1],
    )
    .expect("empty relative complex")
}

fn closed_pair(id: &str, dimension: u8, cells: Vec<u64>, lane: ComplexLane) -> RelativePairSchema {
    let ambient = FiniteComplexSchema::new(
        format!("{id}-ambient"),
        1,
        dimension,
        lane,
        Orientation::Positive,
        cells,
    )
    .expect("ambient complex");
    RelativePairSchema::new(
        id,
        1,
        ambient,
        empty_relative(&format!("{id}-empty"), dimension, lane),
        RelativeModel::MappingCone,
        Vec::new(),
        budget(),
    )
    .expect("closed pair")
}

fn hopkins_singer(
    pair: RelativePairSchema,
    degree: u8,
    terminals: Vec<TerminalTrivialization>,
) -> RelativeDifferentialCharacter {
    let coefficients = CoefficientSystem::RealModuloLattice(flux_lattice(1.0));
    let relative_trivialization = whole_relative(&pair, degree, &coefficients);
    RelativeDifferentialCharacter::new(
        pair,
        CohomologicalDegree::new(degree),
        coefficients,
        RepresentativeKind::HopkinsSingerTriple,
        relative_trivialization,
        terminals,
    )
    .expect("valid Hopkins-Singer schema")
}

fn whole_relative(
    pair: &RelativePairSchema,
    degree: u8,
    coefficients: &CoefficientSystem,
) -> Option<RelativeTrivialization> {
    if pair.relative().is_empty() {
        None
    } else {
        let previous = degree
            .checked_sub(1)
            .expect("nonempty relative fixtures use positive degree");
        Some(
            RelativeTrivialization::new(
                CohomologicalDegree::new(previous),
                pair.ambient().lane(),
                coefficients.clone(),
                "whole-relative-zero",
            )
            .expect("fixture relative trivialization"),
        )
    }
}

fn interval_pair(lane: ComplexLane, boundaries: Vec<BoundaryComponent>) -> RelativePairSchema {
    let interval =
        FiniteComplexSchema::new("interval", 1, 1, lane, Orientation::Positive, vec![2, 1])
            .expect("interval");
    let endpoints = FiniteComplexSchema::new(
        "interval-endpoints",
        1,
        1,
        lane,
        Orientation::Positive,
        vec![2, 0],
    )
    .expect("endpoints");
    RelativePairSchema::new(
        "relative-interval",
        1,
        interval,
        endpoints,
        RelativeModel::MappingCone,
        boundaries,
        budget(),
    )
    .expect("relative interval")
}

fn terminal(id: &str, lane: ComplexLane) -> BoundaryComponent {
    BoundaryComponent::new(
        id,
        0,
        lane,
        BoundaryOrientation::Induced,
        BoundaryRole::Terminal,
    )
    .expect("terminal")
}

fn trivialization(id: &str, lane: ComplexLane, representative: &str) -> TerminalTrivialization {
    TerminalTrivialization::new(
        id,
        CohomologicalDegree::new(0),
        lane,
        CoefficientSystem::RealModuloLattice(flux_lattice(1.0)),
        representative,
    )
    .expect("terminal trivialization")
}

#[test]
fn ra2a_001_circle_and_torus_maps_are_typed_by_degree_and_coefficients() {
    let circle = hopkins_singer(
        closed_pair("circle", 1, vec![1, 1], ComplexLane::Primal),
        1,
        Vec::new(),
    );
    let curvature = circle.curvature_map().expect("curvature map");
    assert_eq!(curvature.kind(), MapKind::Curvature);
    assert_eq!(curvature.degree_shift(), 0);
    assert_eq!(curvature.domain().degree().get(), 1);
    assert_eq!(
        curvature.domain().kind(),
        ObjectKind::DifferentialCharacters
    );
    assert_eq!(
        curvature.codomain().kind(),
        ObjectKind::ClosedIntegralCurvatures
    );
    assert_eq!(
        curvature.codomain().coefficients().sector(),
        CoefficientSector::RealWithLatticePeriods
    );
    let CoefficientSystem::RealWithLatticePeriods(curvature_lattice) =
        curvature.codomain().coefficients()
    else {
        panic!("curvature must retain its period lattice");
    };
    assert_eq!(curvature_lattice.generator_scale(0), Some(1.0));
    let curvature_sequence = circle
        .curvature_exact_sequence()
        .expect("curvature sequence");
    assert_eq!(curvature.codomain(), &curvature_sequence.objects()[3]);

    let torus = hopkins_singer(
        closed_pair("torus", 2, vec![1, 2, 1], ComplexLane::Primal),
        2,
        Vec::new(),
    );
    let characteristic = torus
        .characteristic_class_map()
        .expect("characteristic map");
    assert_eq!(characteristic.kind(), MapKind::CharacteristicClass);
    assert_eq!(characteristic.degree_shift(), 0);
    assert_eq!(
        characteristic.codomain().coefficients().sector(),
        CoefficientSector::Integral
    );
    assert_eq!(characteristic.codomain().degree().get(), 2);
}

#[test]
fn ra2a_002_short_sequences_expose_image_kernel_obligations_without_false_proof() {
    let character = hopkins_singer(
        closed_pair("torus", 2, vec![1, 2, 1], ComplexLane::Primal),
        1,
        Vec::new(),
    );
    for (sequence, kind, middle_map) in [
        (
            character
                .curvature_exact_sequence()
                .expect("curvature sequence"),
            ExactSequenceKind::Curvature,
            MapKind::Curvature,
        ),
        (
            character
                .characteristic_exact_sequence()
                .expect("characteristic sequence"),
            ExactSequenceKind::CharacteristicClass,
            MapKind::CharacteristicClass,
        ),
    ] {
        assert_eq!(sequence.kind(), kind);
        assert_eq!(sequence.objects().len(), 5);
        assert_eq!(sequence.maps().len(), 4);
        assert_eq!(sequence.maps()[2].kind(), middle_map);
        assert_eq!(sequence.exactness_claims().len(), 3);
        sequence.validate_composable().expect("typed composition");
        for (offset, claim) in sequence.exactness_claims().iter().enumerate() {
            assert_eq!(claim.at_object, offset + 1);
            assert_eq!(claim.image_of_map + 1, claim.at_object);
            assert_eq!(claim.kernel_of_map, claim.at_object);
            assert_eq!(
                claim.status,
                ExactnessStatus::RequiresConstructiveWitness {
                    assumption: FiniteExactnessAssumption::ExactFiniteChainComplex,
                }
            );
        }
    }
    let characteristic = character
        .characteristic_exact_sequence()
        .expect("characteristic sequence");
    assert_eq!(
        characteristic.objects()[1].kind(),
        ObjectKind::RealCochainsModuloIntegralCocycles
    );
}

#[test]
fn ra2a_003_relative_boundary_window_has_the_correct_degree_shifts() {
    let ambient = FiniteComplexSchema::new(
        "solid-torus",
        7,
        3,
        ComplexLane::Primal,
        Orientation::Positive,
        vec![8, 16, 12, 4],
    )
    .expect("ambient");
    let relative = FiniteComplexSchema::new(
        "solid-torus-boundary",
        3,
        3,
        ComplexLane::Primal,
        Orientation::Positive,
        vec![4, 8, 4, 0],
    )
    .expect("relative");
    let boundary = BoundaryComponent::new(
        "outer-wall",
        2,
        ComplexLane::Primal,
        BoundaryOrientation::Induced,
        BoundaryRole::Relative,
    )
    .expect("boundary");
    let pair = RelativePairSchema::new(
        "solid-torus-pair",
        2,
        ambient,
        relative,
        RelativeModel::MappingCone,
        vec![boundary],
        budget(),
    )
    .expect("pair");
    let character = hopkins_singer(
        pair,
        1,
        vec![trivialization(
            "outer-wall",
            ComplexLane::Primal,
            "outer-relative-zero",
        )],
    );
    let sequence = character
        .boundary_exact_sequence()
        .expect("relative long-exact window");
    assert_eq!(sequence.kind(), ExactSequenceKind::RelativeBoundary);
    assert_eq!(
        sequence
            .maps()
            .iter()
            .map(|map| map.kind())
            .collect::<Vec<_>>(),
        vec![
            MapKind::Connecting,
            MapKind::ForgetRelative,
            MapKind::BoundaryRestriction,
            MapKind::Connecting,
        ]
    );
    assert_eq!(
        sequence
            .maps()
            .iter()
            .map(|map| map.degree_shift())
            .collect::<Vec<_>>(),
        vec![1, 0, 0, 1]
    );
    assert_eq!(
        sequence
            .objects()
            .iter()
            .map(|object| object.support())
            .collect::<Vec<_>>(),
        vec![
            ObjectSupport::RelativeSubcomplex,
            ObjectSupport::RelativePair,
            ObjectSupport::AmbientComplex,
            ObjectSupport::RelativeSubcomplex,
            ObjectSupport::RelativePair,
        ]
    );
    assert!(sequence.exactness_claims().iter().all(|claim| {
        claim.status
            == ExactnessStatus::RequiresConstructiveWitness {
                assumption: FiniteExactnessAssumption::AdmittedRelativeSubcomplex,
            }
    }));
}

#[test]
fn ra2a_004_terminal_trivializations_are_complete_typed_and_canonical() {
    let left = terminal("left", ComplexLane::Primal);
    let right = terminal("right", ComplexLane::Primal);
    let pair = interval_pair(ComplexLane::Primal, vec![right.clone(), left.clone()]);
    let coefficients = CoefficientSystem::RealModuloLattice(flux_lattice(1.0));

    let missing = RelativeDifferentialCharacter::new(
        pair.clone(),
        CohomologicalDegree::new(1),
        coefficients.clone(),
        RepresentativeKind::HopkinsSingerTriple,
        whole_relative(&pair, 1, &coefficients),
        vec![trivialization("left", ComplexLane::Primal, "left-zero")],
    );
    assert_eq!(
        missing,
        Err(CharacterError::MissingTerminalTrivialization {
            id: "right".to_owned(),
        })
    );

    let wrong_degree = TerminalTrivialization::new(
        "left",
        CohomologicalDegree::new(1),
        ComplexLane::Primal,
        coefficients.clone(),
        "left-wrong-degree",
    )
    .expect("metadata construction is context free");
    assert_eq!(
        RelativeDifferentialCharacter::new(
            pair.clone(),
            CohomologicalDegree::new(1),
            coefficients.clone(),
            RepresentativeKind::HopkinsSingerTriple,
            whole_relative(&pair, 1, &coefficients),
            vec![
                wrong_degree,
                trivialization("right", ComplexLane::Primal, "right-zero"),
            ],
        ),
        Err(CharacterError::BoundaryTrivializationDegreeMismatch {
            expected: 0,
            actual: 1,
        })
    );

    let wrong_lane = TerminalTrivialization::new(
        "left",
        CohomologicalDegree::new(0),
        ComplexLane::Dual,
        coefficients.clone(),
        "left-dual",
    )
    .expect("metadata construction is context free");
    assert!(matches!(
        RelativeDifferentialCharacter::new(
            pair.clone(),
            CohomologicalDegree::new(1),
            coefficients.clone(),
            RepresentativeKind::HopkinsSingerTriple,
            whole_relative(&pair, 1, &coefficients),
            vec![
                wrong_lane,
                trivialization("right", ComplexLane::Primal, "right-zero"),
            ],
        ),
        Err(CharacterError::LaneMismatch {
            object: "boundary trivialization",
            ..
        })
    ));

    let a = RelativeDifferentialCharacter::new(
        pair.clone(),
        CohomologicalDegree::new(1),
        coefficients.clone(),
        RepresentativeKind::HopkinsSingerTriple,
        whole_relative(&pair, 1, &coefficients),
        vec![
            trivialization("right", ComplexLane::Primal, "right-zero"),
            trivialization("left", ComplexLane::Primal, "left-zero"),
        ],
    )
    .expect("valid terminals");
    let pair_b = interval_pair(ComplexLane::Primal, vec![left, right]);
    let b = RelativeDifferentialCharacter::new(
        pair_b.clone(),
        CohomologicalDegree::new(1),
        coefficients.clone(),
        RepresentativeKind::HopkinsSingerTriple,
        whole_relative(&pair_b, 1, &coefficients),
        vec![
            trivialization("left", ComplexLane::Primal, "left-zero"),
            trivialization("right", ComplexLane::Primal, "right-zero"),
        ],
    )
    .expect("valid terminals");
    assert_eq!(a.pair().algebra_id(), b.pair().algebra_id());
    assert_eq!(a.algebra_id(), b.algebra_id());
    assert_eq!(a.boundary_trivializations()[0].boundary_id(), "left");
}

#[test]
fn ra2a_005_primal_dual_and_relative_pair_confusion_refuses() {
    let ambient = FiniteComplexSchema::new(
        "primal-circle",
        1,
        1,
        ComplexLane::Primal,
        Orientation::Positive,
        vec![1, 1],
    )
    .expect("ambient");
    let dual_relative = empty_relative("dual-empty", 1, ComplexLane::Dual);
    assert!(matches!(
        RelativePairSchema::new(
            "mixed-pair",
            1,
            ambient,
            dual_relative,
            RelativeModel::MappingCone,
            Vec::new(),
            budget(),
        ),
        Err(CharacterError::LaneMismatch {
            object: "relative subcomplex",
            ..
        })
    ));

    let primal = hopkins_singer(
        closed_pair("same-shape-primal", 2, vec![1, 2, 1], ComplexLane::Primal),
        1,
        Vec::new(),
    );
    let dual = hopkins_singer(
        closed_pair("same-shape-dual", 2, vec![1, 2, 1], ComplexLane::Dual),
        1,
        Vec::new(),
    );
    assert_eq!(
        primal.cup_product(&dual, &flux_product()),
        Err(CharacterError::PairMismatch)
    );
}

#[test]
fn ra2a_006_integral_torsion_real_and_quotient_sectors_never_coerce() {
    assert_eq!(
        fs_feec::differential_characters::CyclicModulus::new(1),
        Err(CharacterError::InvalidCyclicModulus { modulus: 1 })
    );
    let torsion = CoefficientSystem::Torsion(
        fs_feec::differential_characters::CyclicModulus::new(5).expect("Z/5"),
    );
    let moore_pair = closed_pair("moore-z5", 2, vec![1, 1, 1], ComplexLane::Primal);
    let torsion_object = RelativeDifferentialCharacter::new(
        moore_pair.clone(),
        CohomologicalDegree::new(1),
        torsion.clone(),
        RepresentativeKind::FlatTorsionCocycle,
        None,
        Vec::new(),
    )
    .expect("torsion cocycle object");
    assert_eq!(
        torsion_object.object_space().kind(),
        ObjectKind::TorsionCocycles
    );
    assert_eq!(
        torsion_object.curvature_map(),
        Err(CharacterError::OperationRequiresDifferentialCharacter)
    );
    assert_eq!(
        torsion_object.characteristic_class_map(),
        Err(CharacterError::OperationRequiresDifferentialCharacter)
    );
    assert_eq!(
        torsion
            .characteristic_coefficients()
            .expect("torsion class coefficients")
            .sector(),
        CoefficientSector::Torsion
    );

    let real = CoefficientSystem::Real {
        rank: NonZeroU16::new(1).expect("nonzero"),
        units: Dims::NONE,
    };
    let real_object = RelativeDifferentialCharacter::new(
        moore_pair.clone(),
        CohomologicalDegree::new(1),
        real.clone(),
        RepresentativeKind::RealCochain,
        None,
        Vec::new(),
    )
    .expect("real cochain schema");
    assert_eq!(real_object.object_space().kind(), ObjectKind::RealCochains);
    assert_eq!(
        real_object.characteristic_class_map(),
        Err(CharacterError::OperationRequiresDifferentialCharacter)
    );
    assert_eq!(
        real.characteristic_coefficients(),
        Err(CharacterError::RealHasNoCharacteristicClass)
    );

    assert!(matches!(
        RelativeDifferentialCharacter::new(
            moore_pair,
            CohomologicalDegree::new(1),
            torsion,
            RepresentativeKind::HopkinsSingerTriple,
            None,
            Vec::new(),
        ),
        Err(CharacterError::RepresentativeCoefficientMismatch {
            expected: CoefficientSector::RealModuloLattice,
            actual: CoefficientSector::Torsion,
            ..
        })
    ));
}

#[test]
fn ra2a_007_gauge_product_and_holonomy_degrees_are_explicit() {
    let pair = closed_pair("torus", 2, vec![1, 2, 1], ComplexLane::Primal);
    let first = hopkins_singer(pair.clone(), 1, Vec::new());
    let second = hopkins_singer(pair, 1, Vec::new());

    let gauge = first.gauge_equivalence().expect("gauge schema");
    assert_eq!(gauge.gauge_parameters().degree().get(), 0);
    assert_eq!(gauge.representatives().degree().get(), 1);
    assert_eq!(gauge.quotient().degree().get(), 1);
    assert_eq!(
        gauge.cellular_nilpotence(),
        NilpotenceLaw::CellularCoboundarySquaredZero
    );
    assert_eq!(
        gauge.de_rham_nilpotence(),
        NilpotenceLaw::DeRhamDifferentialSquaredZero
    );

    let wrong_product = CoefficientProductSchema::new(
        "wrong-left-coefficient",
        1,
        CoefficientSystem::Integral(flux_lattice(1.0)),
        second.coefficients().clone(),
        second.coefficients().clone(),
        AlgebraId::from_bytes([0x24; 32]),
        budget(),
    )
    .expect("well-formed but inapplicable product");
    assert_eq!(
        first.cup_product(&second, &wrong_product),
        Err(CharacterError::CoefficientProductMismatch { input: "left" })
    );

    let product_rule = flux_product();
    let product = first
        .cup_product(&second, &product_rule)
        .expect("cup product");
    assert_eq!(product.kind(), BilinearMapKind::CupProduct);
    assert_eq!(product.left().degree().get(), 1);
    assert_eq!(product.right().degree().get(), 1);
    assert_eq!(product.output().degree().get(), 2);
    assert_eq!(product.output().coefficients(), product_rule.output());
    assert_eq!(
        product_rule.map_artifact(),
        AlgebraId::from_bytes([0x42; 32])
    );

    let holonomy = first.holonomy_pairing().expect("holonomy");
    assert_eq!(holonomy.kind(), BilinearMapKind::HolonomyPairing);
    assert_eq!(holonomy.right().kind(), ObjectKind::RelativeCycles);
    assert_eq!(holonomy.right().degree().get(), 0);
    let CoefficientSystem::Integral(cycle_lattice) = holonomy.right().coefficients() else {
        panic!("holonomy must consume integral cycles");
    };
    assert_eq!(cycle_lattice.units(), Dims::NONE);
    assert_eq!(cycle_lattice.generator_scale(0), Some(1.0));
    assert_eq!(holonomy.output().support(), ObjectSupport::Point);

    let circle_pair = closed_pair("circle", 1, vec![1, 1], ComplexLane::Primal);
    let circle_degree_one = hopkins_singer(circle_pair.clone(), 1, Vec::new());
    let circle_degree_two = hopkins_singer(circle_pair, 2, Vec::new());
    assert_eq!(
        circle_degree_one
            .cup_product(&circle_degree_one, &flux_product())
            .expect("degree two is valid on a circle")
            .output()
            .degree()
            .get(),
        2
    );
    assert_eq!(
        circle_degree_two.cup_product(&circle_degree_one, &flux_product()),
        Err(CharacterError::ProductDegreeOutOfRange)
    );
}

#[test]
fn ra2a_008_identity_replay_is_stable_and_semantic_mutations_move_it() {
    let baseline = hopkins_singer(
        closed_pair("torus", 2, vec![1, 2, 1], ComplexLane::Primal),
        1,
        Vec::new(),
    );
    let replay = hopkins_singer(
        closed_pair("torus", 2, vec![1, 2, 1], ComplexLane::Primal),
        1,
        Vec::new(),
    );
    assert_eq!(baseline.canonical_bytes(), replay.canonical_bytes());
    assert_eq!(baseline.algebra_id(), replay.algebra_id());
    assert_eq!(baseline.algebra_id().to_hex().len(), 64);

    let changed_degree = hopkins_singer(
        closed_pair("torus", 2, vec![1, 2, 1], ComplexLane::Primal),
        2,
        Vec::new(),
    );
    assert_ne!(baseline.algebra_id(), changed_degree.algebra_id());

    let changed_scale = RelativeDifferentialCharacter::new(
        baseline.pair().clone(),
        CohomologicalDegree::new(1),
        CoefficientSystem::RealModuloLattice(flux_lattice(2.0)),
        RepresentativeKind::HopkinsSingerTriple,
        None,
        Vec::new(),
    )
    .expect("changed normalization");
    assert_ne!(baseline.algebra_id(), changed_scale.algebra_id());

    let negative_ambient = FiniteComplexSchema::new(
        "torus-ambient",
        1,
        2,
        ComplexLane::Primal,
        Orientation::Negative,
        vec![1, 2, 1],
    )
    .expect("negative torus");
    let changed_orientation = RelativePairSchema::new(
        "torus",
        1,
        negative_ambient,
        empty_relative("torus-empty", 2, ComplexLane::Primal),
        RelativeModel::MappingCone,
        Vec::new(),
        budget(),
    )
    .expect("orientation mutation");
    assert_ne!(
        baseline.pair().algebra_id(),
        changed_orientation.algebra_id()
    );

    let product = flux_product();
    let product_replay = flux_product();
    assert_eq!(product.canonical_bytes(), product_replay.canonical_bytes());
    assert_eq!(product.algebra_id(), product_replay.algebra_id());
    let product_coefficients = CoefficientSystem::RealModuloLattice(flux_lattice(1.0));
    let changed_map_artifact = CoefficientProductSchema::new(
        "flux-cup-product",
        1,
        product_coefficients.clone(),
        product_coefficients.clone(),
        product_coefficients,
        AlgebraId::from_bytes([0x43; 32]),
        budget(),
    )
    .expect("changed coefficient-map artifact");
    assert_ne!(product.algebra_id(), changed_map_artifact.algebra_id());
}

#[test]
fn ra2a_009_degree_coefficient_boundary_and_budget_mutations_refuse() {
    let pair = closed_pair("circle", 1, vec![1, 1], ComplexLane::Primal);
    let top_flat = RelativeDifferentialCharacter::new(
        pair.clone(),
        CohomologicalDegree::new(2),
        CoefficientSystem::RealModuloLattice(flux_lattice(1.0)),
        RepresentativeKind::HopkinsSingerTriple,
        None,
        Vec::new(),
    )
    .expect("degree dim+1 flat characters are representable");
    assert_eq!(top_flat.degree().get(), 2);
    assert_eq!(
        RelativeDifferentialCharacter::new(
            pair.clone(),
            CohomologicalDegree::new(3),
            CoefficientSystem::RealModuloLattice(flux_lattice(1.0)),
            RepresentativeKind::HopkinsSingerTriple,
            None,
            Vec::new(),
        ),
        Err(CharacterError::DegreeOutOfRange {
            degree: 3,
            maximum: 2,
        })
    );

    assert!(matches!(
        CoefficientLattice::new("bad", 1, Dims::NONE, &[f64::NAN]),
        Err(CharacterError::InvalidLatticeScale { index: 0 })
    ));
    assert_eq!(
        RelativePairSchema::new(
            "degree-representation-overflow",
            1,
            FiniteComplexSchema::new(
                "dimension-255-ambient",
                1,
                u8::MAX,
                ComplexLane::Primal,
                Orientation::Positive,
                vec![0; usize::from(u8::MAX) + 1],
            )
            .expect("finite dimension-255 schema"),
            empty_relative("dimension-255-empty", u8::MAX, ComplexLane::Primal),
            RelativeModel::MappingCone,
            Vec::new(),
            budget(),
        ),
        Err(CharacterError::CharacterDegreeRepresentationOverflow { dimension: u8::MAX })
    );
    assert_eq!(
        FiniteComplexSchema::new(
            "bad-arity",
            1,
            2,
            ComplexLane::Primal,
            Orientation::Positive,
            vec![1, 1],
        ),
        Err(CharacterError::CellCountArity {
            dimension: 2,
            actual: 2,
        })
    );

    let too_small = AlgebraBudget::new(1, 1, 128).expect("small budget");
    assert!(matches!(
        RelativePairSchema::new(
            "over-budget",
            1,
            FiniteComplexSchema::new(
                "ambient",
                1,
                1,
                ComplexLane::Primal,
                Orientation::Positive,
                vec![2, 1],
            )
            .expect("ambient"),
            empty_relative("empty", 1, ComplexLane::Primal),
            RelativeModel::MappingCone,
            Vec::new(),
            too_small,
        ),
        Err(CharacterError::CellBudgetExceeded { .. })
    ));

    let tiny_canonical = AlgebraBudget::new(100, 1, 16).expect("tiny canonical budget");
    assert!(matches!(
        RelativePairSchema::new(
            "canonical-over-budget",
            1,
            FiniteComplexSchema::new(
                "point",
                1,
                0,
                ComplexLane::Primal,
                Orientation::Positive,
                vec![1],
            )
            .expect("point"),
            empty_relative("empty-point", 0, ComplexLane::Primal),
            RelativeModel::MappingCone,
            Vec::new(),
            tiny_canonical,
        ),
        Err(CharacterError::CanonicalBudgetExceeded { limit: 16, .. })
    ));

    let product_coefficients = CoefficientSystem::RealModuloLattice(flux_lattice(1.0));
    assert!(matches!(
        CoefficientProductSchema::new(
            "canonical-product-over-budget",
            1,
            product_coefficients.clone(),
            product_coefficients.clone(),
            product_coefficients,
            AlgebraId::from_bytes([0x42; 32]),
            tiny_canonical,
        ),
        Err(CharacterError::CanonicalBudgetExceeded { limit: 16, .. })
    ));

    let relative_boundary = BoundaryComponent::new(
        "wall",
        0,
        ComplexLane::Primal,
        BoundaryOrientation::Induced,
        BoundaryRole::Relative,
    )
    .expect("wall");
    let relative_pair = interval_pair(ComplexLane::Primal, vec![relative_boundary]);
    assert_eq!(
        RelativeDifferentialCharacter::new(
            relative_pair,
            CohomologicalDegree::new(1),
            CoefficientSystem::RealModuloLattice(flux_lattice(1.0)),
            RepresentativeKind::HopkinsSingerTriple,
            None,
            Vec::new(),
        ),
        Err(CharacterError::MissingRelativeSubcomplexTrivialization)
    );

    let interface = BoundaryComponent::new(
        "interface",
        0,
        ComplexLane::Primal,
        BoundaryOrientation::Induced,
        BoundaryRole::Interface,
    )
    .expect("interface");
    let interface_pair = interval_pair(ComplexLane::Primal, vec![interface]);
    let interface_coefficients = CoefficientSystem::RealModuloLattice(flux_lattice(1.0));
    assert_eq!(
        RelativeDifferentialCharacter::new(
            interface_pair.clone(),
            CohomologicalDegree::new(1),
            interface_coefficients.clone(),
            RepresentativeKind::HopkinsSingerTriple,
            whole_relative(&interface_pair, 1, &interface_coefficients),
            vec![trivialization(
                "interface",
                ComplexLane::Primal,
                "bad-interface-trivialization",
            )],
        ),
        Err(CharacterError::BoundaryDoesNotAdmitTrivialization {
            id: "interface".to_owned(),
        })
    );
}

#[test]
fn ra2a_010_degree_zero_without_terminals_is_valid_but_negative_degree_maps_refuse() {
    let object = hopkins_singer(
        closed_pair("point", 0, vec![1], ComplexLane::Primal),
        0,
        Vec::new(),
    );
    assert_eq!(object.degree().get(), 0);
    assert_eq!(
        object.gauge_equivalence(),
        Err(CharacterError::DegreeZeroHasNoGaugeParameter)
    );
    assert_eq!(
        object.curvature_exact_sequence(),
        Err(CharacterError::DegreeZeroHasNoExactSequence)
    );
    assert_eq!(
        object.holonomy_pairing(),
        Err(CharacterError::DegreeZeroHasNoHolonomyCycle)
    );
}
