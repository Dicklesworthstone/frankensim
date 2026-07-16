//! G0 battery for the multi-field system IR (bead i94v.1.1.1):
//! type/algebra laws, identity metamorphics (rename/permutation
//! invariance, convention sensitivity), affine-quantity misuse,
//! clock/frame/pairing refusals, and depth/resource adversaries.

use fs_couple::{PortKind, PortOrientation};
use fs_opdsl::{
    AtomSignature, BlockEquation, ConventionRef, CoordinateConvention, FieldDecl, FieldId,
    FieldQuantity, ScalarConvention, Space, SpatialSupport, StateOwnership, SystemDef, SystemExpr,
    SystemTypeError,
};
use fs_qty::Dims;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

const KELVIN: Dims = Dims([0, 0, 0, 1, 0, 0]);
const VELOCITY: Dims = Dims([1, 0, -1, 0, 0, 0]);
const FORCE: Dims = Dims([1, 1, -2, 0, 0, 0]);

fn refs(basis: &str, frame: &str, clock: &str) -> (CoordinateConvention, ConventionRef) {
    (
        CoordinateConvention {
            basis: ConventionRef::new(basis).expect("basis ref"),
            frame: ConventionRef::new(frame).expect("frame ref"),
            orientation: PortOrientation::AlongFrame,
        },
        ConventionRef::new(clock).expect("clock ref"),
    )
}

fn field(
    name: &str,
    degree: u8,
    n: usize,
    quantity: FieldQuantity,
    conventions: (&str, &str, &str),
    slot: u32,
) -> FieldDecl {
    let (coordinates, clock) = refs(conventions.0, conventions.1, conventions.2);
    FieldDecl {
        name: name.to_string(),
        space: Space {
            degree,
            n,
            dims: quantity.dims(),
        },
        quantity,
        coordinates,
        clock,
        support: SpatialSupport::Interior,
        state: StateOwnership::Owned { slot },
    }
}

fn two_field_system(names: (&str, &str), declaration_order_swapped: bool) -> fs_opdsl::SystemId {
    let mut system = SystemDef::new();
    let velocity = field(
        names.0,
        1,
        64,
        FieldQuantity::Dimensional(VELOCITY),
        ("chart-a", "lab", "clk-main"),
        0,
    );
    let force = field(
        names.1,
        1,
        64,
        FieldQuantity::Dimensional(FORCE),
        ("chart-a", "lab", "clk-main"),
        1,
    );
    let (first, second) = if declaration_order_swapped {
        (force.clone(), velocity.clone())
    } else {
        (velocity.clone(), force.clone())
    };
    let a = system.declare_field(first).expect("first field");
    let b = system.declare_field(second).expect("second field");
    // Map the handles back to (velocity, force) regardless of order.
    let (v, f) = if declaration_order_swapped {
        (b, a)
    } else {
        (a, b)
    };
    system
        .add_equation(BlockEquation {
            name: "coupling-power".to_string(),
            target: v,
            rhs: SystemExpr::Scale(0.5, Box::new(SystemExpr::FieldRef(v))),
        })
        .expect("velocity equation");
    system
        .add_equation(BlockEquation {
            name: "force-balance".to_string(),
            target: f,
            rhs: SystemExpr::FieldRef(f),
        })
        .expect("force equation");
    system.admit().expect("system admits").identity()
}

#[test]
fn sys_001_identity_ignores_names_and_declaration_order() {
    let baseline = two_field_system(("velocity", "force"), false);
    let renamed = two_field_system(("u", "sigma"), false);
    let reordered = two_field_system(("velocity", "force"), true);
    assert_eq!(baseline, renamed, "display renaming must preserve identity");
    assert_eq!(
        baseline, reordered,
        "declaration/serialization order must preserve identity"
    );
}

#[test]
fn sys_002_convention_changes_move_identity() {
    let baseline = two_field_system(("velocity", "force"), false);

    // A different frame reference on one field is a meaningful change.
    let mut system = SystemDef::new();
    let v = system
        .declare_field(field(
            "velocity",
            1,
            64,
            FieldQuantity::Dimensional(VELOCITY),
            ("chart-a", "rotor", "clk-main"),
            0,
        ))
        .expect("field");
    let f = system
        .declare_field(field(
            "force",
            1,
            64,
            FieldQuantity::Dimensional(FORCE),
            ("chart-a", "rotor", "clk-main"),
            1,
        ))
        .expect("field");
    system
        .add_equation(BlockEquation {
            name: "coupling-power".to_string(),
            target: v,
            rhs: SystemExpr::Scale(0.5, Box::new(SystemExpr::FieldRef(v))),
        })
        .expect("eq");
    system
        .add_equation(BlockEquation {
            name: "force-balance".to_string(),
            target: f,
            rhs: SystemExpr::FieldRef(f),
        })
        .expect("eq");
    let rotor_frame = system.admit().expect("admits").identity();
    assert_ne!(
        baseline, rotor_frame,
        "a frame convention change must move identity"
    );

    // A unit/dims rescale (velocity -> force dims on the first field)
    // is likewise identity-bearing (covered by distinct dims already);
    // scalar-convention change moves identity too.
    let complex = {
        let mut system = SystemDef::new().scalar_convention(ScalarConvention::ComplexHermitian);
        let v = system
            .declare_field(field(
                "velocity",
                1,
                64,
                FieldQuantity::Dimensional(VELOCITY),
                ("chart-a", "lab", "clk-main"),
                0,
            ))
            .expect("field");
        let f = system
            .declare_field(field(
                "force",
                1,
                64,
                FieldQuantity::Dimensional(FORCE),
                ("chart-a", "lab", "clk-main"),
                1,
            ))
            .expect("field");
        system
            .add_equation(BlockEquation {
                name: "coupling-power".to_string(),
                target: v,
                rhs: SystemExpr::Scale(0.5, Box::new(SystemExpr::FieldRef(v))),
            })
            .expect("eq");
        system
            .add_equation(BlockEquation {
                name: "force-balance".to_string(),
                target: f,
                rhs: SystemExpr::FieldRef(f),
            })
            .expect("eq");
        system.admit().expect("admits").identity()
    };
    assert_ne!(
        baseline, complex,
        "the scalar convention is identity-bearing"
    );
}

#[test]
fn sys_003_affine_temperature_misuse_refuses_with_minimal_diagnostics() {
    let mut system = SystemDef::new();
    let temperature = system
        .declare_field(field(
            "temperature",
            0,
            32,
            FieldQuantity::Semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            ("chart-a", "lab", "clk-main"),
            0,
        ))
        .expect("temperature field");

    let scaled = system.add_equation(BlockEquation {
        name: "bad-scale".to_string(),
        target: temperature,
        rhs: SystemExpr::Scale(2.0, Box::new(SystemExpr::FieldRef(temperature))),
    });
    assert!(
        matches!(
            scaled,
            Err(SystemTypeError::AffineQuantityMisuse { operation: "scaling", ref field })
                if field == "temperature"
        ),
        "scaling an absolute temperature must refuse by name, got {scaled:?}"
    );

    let summed = system.add_equation(BlockEquation {
        name: "bad-sum".to_string(),
        target: temperature,
        rhs: SystemExpr::Add(
            Box::new(SystemExpr::FieldRef(temperature)),
            Box::new(SystemExpr::FieldRef(temperature)),
        ),
    });
    assert!(
        matches!(
            summed,
            Err(SystemTypeError::AffineQuantityMisuse {
                operation: "summation",
                ..
            })
        ),
        "summing absolute temperatures must refuse, got {summed:?}"
    );

    // A temperature DIFFERENCE field composes freely.
    let mut system = SystemDef::new();
    let delta = system
        .declare_field(field(
            "delta-t",
            0,
            32,
            FieldQuantity::Semantic(SemanticType::new(
                QuantityKind::TemperatureDifference,
                ValueForm::Static,
            )),
            ("chart-a", "lab", "clk-main"),
            0,
        ))
        .expect("difference field");
    system
        .add_equation(BlockEquation {
            name: "relax".to_string(),
            target: delta,
            rhs: SystemExpr::Scale(0.25, Box::new(SystemExpr::FieldRef(delta))),
        })
        .expect("temperature differences scale freely");
    system.admit().expect("difference system admits");
}

#[test]
fn sys_004_semantic_kind_dims_must_match_space() {
    let mut system = SystemDef::new();
    let (coordinates, clock) = refs("chart-a", "lab", "clk-main");
    let wrong = system.declare_field(FieldDecl {
        name: "temperature".to_string(),
        space: Space {
            degree: 0,
            n: 8,
            dims: VELOCITY, // wrong: kind expects kelvin
        },
        quantity: FieldQuantity::Semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature,
            ValueForm::Static,
        )),
        coordinates,
        clock,
        support: SpatialSupport::Interior,
        state: StateOwnership::Owned { slot: 0 },
    });
    assert!(
        matches!(
            wrong,
            Err(SystemTypeError::QuantityDimsMismatch { ref space_dims, ref kind_dims, .. })
                if *space_dims == VELOCITY && *kind_dims == KELVIN
        ),
        "kind/space dims disagreement must refuse, got {wrong:?}"
    );
}

#[test]
fn sys_005_clock_and_frame_mismatches_refuse_before_lowering() {
    let mut system = SystemDef::new();
    let fast = system
        .declare_field(field(
            "fast",
            0,
            16,
            FieldQuantity::Dimensional(VELOCITY),
            ("chart-a", "lab", "clk-fast"),
            0,
        ))
        .expect("fast field");
    let slow = system
        .declare_field(field(
            "slow",
            0,
            16,
            FieldQuantity::Dimensional(VELOCITY),
            ("chart-a", "lab", "clk-slow"),
            1,
        ))
        .expect("slow field");
    let mixed_clock = system.add_equation(BlockEquation {
        name: "mixed".to_string(),
        target: fast,
        rhs: SystemExpr::Add(
            Box::new(SystemExpr::FieldRef(fast)),
            Box::new(SystemExpr::FieldRef(slow)),
        ),
    });
    assert!(
        matches!(mixed_clock, Err(SystemTypeError::ClockMismatch { .. })),
        "cross-clock summation must refuse, got {mixed_clock:?}"
    );

    let mut system = SystemDef::new();
    let lab = system
        .declare_field(field(
            "lab-velocity",
            0,
            16,
            FieldQuantity::Dimensional(VELOCITY),
            ("chart-a", "lab", "clk-main"),
            0,
        ))
        .expect("lab field");
    let rotor = system
        .declare_field(field(
            "rotor-velocity",
            0,
            16,
            FieldQuantity::Dimensional(VELOCITY),
            ("chart-a", "rotor", "clk-main"),
            1,
        ))
        .expect("rotor field");
    let mixed_frame = system.add_equation(BlockEquation {
        name: "mixed-frame".to_string(),
        target: lab,
        rhs: SystemExpr::Add(
            Box::new(SystemExpr::FieldRef(lab)),
            Box::new(SystemExpr::FieldRef(rotor)),
        ),
    });
    assert!(
        matches!(mixed_frame, Err(SystemTypeError::ConventionMismatch { .. })),
        "cross-frame summation must refuse, got {mixed_frame:?}"
    );
}

#[test]
fn sys_006_port_pairings_require_power_conjugacy() {
    let mut system = SystemDef::new();
    let velocity = system
        .declare_field(field(
            "velocity",
            0,
            8,
            FieldQuantity::Dimensional(VELOCITY),
            ("chart-a", "lab", "clk-main"),
            0,
        ))
        .expect("velocity");
    let force = system
        .declare_field(field(
            "force",
            0,
            8,
            FieldQuantity::Dimensional(FORCE),
            ("chart-a", "lab", "clk-main"),
            1,
        ))
        .expect("force");
    let power = system
        .declare_field(field(
            "power",
            0,
            1,
            FieldQuantity::Dimensional(Dims([2, 1, -3, 0, 0, 0])),
            ("chart-a", "lab", "clk-main"),
            2,
        ))
        .expect("power");

    // force x velocity = power: conjugate, admitted.
    system
        .add_equation(BlockEquation {
            name: "interface-power".to_string(),
            target: power,
            rhs: SystemExpr::PortPair {
                kind: PortKind::MechanicalForceVelocity,
                effort: Box::new(SystemExpr::FieldRef(force)),
                flow: Box::new(SystemExpr::FieldRef(velocity)),
                measure_dims: Dims::NONE,
            },
        })
        .expect("conjugate pairing admits");

    // velocity x velocity is NOT power: refused with both sides named.
    let bad = system.add_equation(BlockEquation {
        name: "bad-pair".to_string(),
        target: power,
        rhs: SystemExpr::PortPair {
            kind: PortKind::MechanicalForceVelocity,
            effort: Box::new(SystemExpr::FieldRef(velocity)),
            flow: Box::new(SystemExpr::FieldRef(velocity)),
            measure_dims: Dims::NONE,
        },
    });
    assert!(
        matches!(
            bad,
            Err(SystemTypeError::NonConjugatePairing { ref effort_dims, .. })
                if *effort_dims == VELOCITY
        ),
        "non-conjugate pairing must refuse, got {bad:?}"
    );
}

#[test]
fn sys_007_indistinguishable_fields_and_duplicate_slots_refuse() {
    let mut system = SystemDef::new();
    system
        .declare_field(field(
            "a",
            0,
            16,
            FieldQuantity::Dimensional(VELOCITY),
            ("chart-a", "lab", "clk-main"),
            0,
        ))
        .expect("first");
    let duplicate_slot = system.declare_field(field(
        "b",
        0,
        16,
        FieldQuantity::Dimensional(FORCE),
        ("chart-a", "lab", "clk-main"),
        0,
    ));
    assert!(
        matches!(
            duplicate_slot,
            Err(SystemTypeError::DuplicateStateSlot { slot: 0, .. })
        ),
        "duplicate owned slot must refuse, got {duplicate_slot:?}"
    );

    // Byte-identical payloads (External state, same everything) are
    // ambiguous for canonical ordering and refuse at admit().
    let mut system = SystemDef::new();
    for name in ["left", "right"] {
        let (coordinates, clock) = refs("chart-a", "lab", "clk-main");
        system
            .declare_field(FieldDecl {
                name: name.to_string(),
                space: Space {
                    degree: 0,
                    n: 16,
                    dims: VELOCITY,
                },
                quantity: FieldQuantity::Dimensional(VELOCITY),
                coordinates,
                clock,
                support: SpatialSupport::Interior,
                state: StateOwnership::External,
            })
            .expect("declare");
    }
    let ambiguous = system.admit();
    assert!(
        matches!(
            ambiguous,
            Err(SystemTypeError::IndistinguishableFields { ref first, ref second })
                if first == "left" && second == "right"
        ),
        "identical payloads must refuse as ambiguous, got {:?}",
        ambiguous.map(|admitted| admitted.identity())
    );
}

#[test]
fn sys_008_depth_and_reference_adversaries_refuse_structurally() {
    let mut system = SystemDef::new();
    let base = system
        .declare_field(field(
            "base",
            0,
            4,
            FieldQuantity::Dimensional(VELOCITY),
            ("chart-a", "lab", "clk-main"),
            0,
        ))
        .expect("base field");

    // Nesting beyond the cap refuses without exhausting the stack.
    let mut deep = SystemExpr::FieldRef(base);
    for _ in 0..(fs_opdsl::system::MAX_SYSTEM_EXPR_DEPTH + 8) {
        deep = SystemExpr::Scale(0.5, Box::new(deep));
    }
    let too_deep = system.add_equation(BlockEquation {
        name: "deep".to_string(),
        target: base,
        rhs: deep,
    });
    assert!(
        matches!(too_deep, Err(SystemTypeError::DepthExceeded { .. })),
        "adversarial nesting must refuse, got a different outcome"
    );

    // Dangling references refuse by table name.
    let dangling = system.add_equation(BlockEquation {
        name: "dangling".to_string(),
        target: base,
        rhs: SystemExpr::FieldRef(FieldId(999)),
    });
    assert!(
        matches!(
            dangling,
            Err(SystemTypeError::UnknownId {
                what: "field",
                id: 999
            })
        ),
        "dangling field ref must refuse, got {dangling:?}"
    );

    // Non-finite scale constants refuse with the exact bits.
    let non_finite = system.add_equation(BlockEquation {
        name: "nan".to_string(),
        target: base,
        rhs: SystemExpr::Scale(f64::NAN, Box::new(SystemExpr::FieldRef(base))),
    });
    assert!(
        matches!(non_finite, Err(SystemTypeError::NonFiniteScale { .. })),
        "non-finite scale must refuse, got {non_finite:?}"
    );
}

#[test]
fn sys_009_atom_application_and_extension_semantics() {
    let mut system = SystemDef::new();
    let velocity = system
        .declare_field(field(
            "velocity",
            1,
            64,
            FieldQuantity::Dimensional(VELOCITY),
            ("chart-a", "lab", "clk-main"),
            0,
        ))
        .expect("velocity");
    let grad = system.register_atom(AtomSignature {
        name: "d0".to_string(),
        content: ConventionRef::new("feec:d0:v1").expect("content ref"),
        in_space: Space {
            degree: 1,
            n: 64,
            dims: VELOCITY,
        },
        out_space: Space {
            degree: 2,
            n: 96,
            dims: VELOCITY,
        },
    });
    // Applying to the wrong space refuses with both spaces.
    let wrong = system.add_equation(BlockEquation {
        name: "wrong-space".to_string(),
        target: velocity,
        rhs: SystemExpr::Apply {
            atom: grad,
            arg: Box::new(SystemExpr::Apply {
                atom: grad,
                arg: Box::new(SystemExpr::FieldRef(velocity)),
            }),
        },
    });
    assert!(
        matches!(
            wrong,
            Err(SystemTypeError::SpaceMismatch {
                context: "atom application",
                ..
            })
        ),
        "space-mismatched application must refuse, got {wrong:?}"
    );

    // The opaque extension is identity-bearing.
    let plain = {
        let mut system = SystemDef::new();
        let v = system
            .declare_field(field(
                "velocity",
                1,
                64,
                FieldQuantity::Dimensional(VELOCITY),
                ("chart-a", "lab", "clk-main"),
                0,
            ))
            .expect("velocity");
        system
            .add_equation(BlockEquation {
                name: "id".to_string(),
                target: v,
                rhs: SystemExpr::FieldRef(v),
            })
            .expect("eq");
        system.admit().expect("admits").identity()
    };
    let extended = {
        let mut system = SystemDef::new();
        let v = system
            .declare_field(field(
                "velocity",
                1,
                64,
                FieldQuantity::Dimensional(VELOCITY),
                ("chart-a", "lab", "clk-main"),
                0,
            ))
            .expect("velocity");
        system
            .add_equation(BlockEquation {
                name: "id".to_string(),
                target: v,
                rhs: SystemExpr::FieldRef(v),
            })
            .expect("eq");
        system
            .with_extension(b"dialect-x-v0".to_vec())
            .expect("bounded extension")
            .admit()
            .expect("admits")
            .identity()
    };
    assert_ne!(plain, extended, "the opaque extension is identity-bearing");
}

#[test]
fn sys_010_admission_is_deterministic() {
    let first = two_field_system(("velocity", "force"), false);
    let second = two_field_system(("velocity", "force"), false);
    assert_eq!(first, second, "identical systems mint identical identities");
}

#[test]
fn sys_011_transport_round_trips_and_preserves_identity() {
    let mut system = SystemDef::new();
    let velocity = system
        .declare_field(field(
            "velocity",
            1,
            64,
            FieldQuantity::Dimensional(VELOCITY),
            ("chart-a", "lab", "clk-main"),
            0,
        ))
        .expect("velocity");
    let temperature = system
        .declare_field(field(
            "delta-t",
            0,
            32,
            FieldQuantity::Semantic(SemanticType::new(
                QuantityKind::TemperatureDifference,
                ValueForm::Static,
            )),
            ("chart-a", "lab", "clk-main"),
            1,
        ))
        .expect("temperature difference");
    let grad = system.register_atom(AtomSignature {
        name: "d0".to_string(),
        content: ConventionRef::new("feec:d0:v1").expect("content ref"),
        in_space: Space {
            degree: 1,
            n: 64,
            dims: VELOCITY,
        },
        out_space: Space {
            degree: 1,
            n: 64,
            dims: VELOCITY,
        },
    });
    system
        .add_equation(BlockEquation {
            name: "momentum".to_string(),
            target: velocity,
            rhs: SystemExpr::Add(
                Box::new(SystemExpr::Apply {
                    atom: grad,
                    arg: Box::new(SystemExpr::FieldRef(velocity)),
                }),
                Box::new(SystemExpr::Scale(
                    0.25,
                    Box::new(SystemExpr::FieldRef(velocity)),
                )),
            ),
        })
        .expect("momentum equation");
    system
        .add_equation(BlockEquation {
            name: "heat".to_string(),
            target: temperature,
            rhs: SystemExpr::Scale(0.5, Box::new(SystemExpr::FieldRef(temperature))),
        })
        .expect("heat equation");
    let admitted = system
        .with_extension(b"golden-v1".to_vec())
        .expect("extension")
        .admit()
        .expect("admits");

    let text = fs_opdsl::system::transport::to_text(&admitted).expect("serializes");
    let reparsed = fs_opdsl::system::transport::from_text(&text)
        .expect("parses")
        .admit()
        .expect("re-admits");
    assert_eq!(
        admitted.identity(),
        reparsed.identity(),
        "transport round trip must preserve semantic identity"
    );
    let text_again = fs_opdsl::system::transport::to_text(&reparsed).expect("serializes again");
    assert_eq!(text, text_again, "the transport text is canonical");
}

#[test]
fn sys_012_transport_refuses_other_versions_and_malformed_lines() {
    let versioned = "fs-opdsl-system-transport-v1\nversion\t2\nconvention\treal-only\n";
    let refusal = fs_opdsl::system::transport::from_text(versioned);
    assert!(
        matches!(
            refusal,
            Err(SystemTypeError::VersionMismatch {
                found: 2,
                supported: 1
            })
        ),
        "a future IR version must refuse pending audited migration, got {refusal:?}"
    );

    let garbled =
        "fs-opdsl-system-transport-v1\nversion\t1\nconvention\treal-only\nfield\tonly-a-name\n";
    let refusal = fs_opdsl::system::transport::from_text(garbled);
    assert!(
        matches!(
            refusal,
            Err(SystemTypeError::TransportMalformed { line: 4, .. })
        ),
        "a malformed field record must refuse with its line, got {refusal:?}"
    );

    let bad_magic = "some-other-transport\nversion\t1\n";
    assert!(matches!(
        fs_opdsl::system::transport::from_text(bad_magic),
        Err(SystemTypeError::TransportMalformed { line: 1, .. })
    ));
}

#[test]
fn sys_013_tampered_transport_cannot_alias_the_original_identity() {
    let baseline = two_field_system(("velocity", "force"), false);
    let mut system = SystemDef::new();
    let v = system
        .declare_field(field(
            "velocity",
            1,
            64,
            FieldQuantity::Dimensional(VELOCITY),
            ("chart-a", "lab", "clk-main"),
            0,
        ))
        .expect("velocity");
    let f = system
        .declare_field(field(
            "force",
            1,
            64,
            FieldQuantity::Dimensional(FORCE),
            ("chart-a", "lab", "clk-main"),
            1,
        ))
        .expect("force");
    system
        .add_equation(BlockEquation {
            name: "coupling-power".to_string(),
            target: v,
            rhs: SystemExpr::Scale(0.5, Box::new(SystemExpr::FieldRef(v))),
        })
        .expect("eq");
    system
        .add_equation(BlockEquation {
            name: "force-balance".to_string(),
            target: f,
            rhs: SystemExpr::FieldRef(f),
        })
        .expect("eq");
    let admitted = system.admit().expect("admits");
    let text = fs_opdsl::system::transport::to_text(&admitted).expect("serializes");

    // Tamper: move one field to a different frame reference.
    let tampered = text.replace("\tlab\t", "\trotor\t");
    assert_ne!(text, tampered, "fixture must contain the frame reference");
    let outcome = fs_opdsl::system::transport::from_text(&tampered)
        .expect("still structurally valid")
        .admit();
    // Refusing outright would be equally fail-closed; when it admits,
    // the identity must have moved.
    if let Ok(reparsed) = outcome {
        assert_ne!(
            reparsed.identity(),
            baseline,
            "a tampered convention must never alias the original identity"
        );
    }
}

#[test]
fn sys_014_migration_golden_identity_is_pinned() {
    // The golden system: any change to canonical encoding, field payload
    // layout, remapping, or the identity schema moves this hex and must
    // arrive as a DELIBERATE version bump with a migration note
    // (docs/GOLDEN_POLICY.md discipline).
    let identity = two_field_system(("velocity", "force"), false);
    let hex = identity.to_string();
    assert_eq!(
        hex, "2db84c5e4cab57a2b1bac5c0ebe109062ddf21464ddd6365d4c71353c4865bc8",
        "system-identity golden moved: bump SYSTEM_IR_VERSION deliberately and record the cause"
    );
}

#[test]
fn sys_015_atom_identity_is_order_and_name_free_but_content_bearing() {
    let build = |swap_order: bool, grad_content: &str, names: (&str, &str)| {
        let mut system = SystemDef::new();
        let velocity = system
            .declare_field(field(
                "velocity",
                1,
                64,
                FieldQuantity::Dimensional(VELOCITY),
                ("chart-a", "lab", "clk-main"),
                0,
            ))
            .expect("velocity");
        let space = Space {
            degree: 1,
            n: 64,
            dims: VELOCITY,
        };
        let grad_sig = AtomSignature {
            name: names.0.to_string(),
            content: ConventionRef::new(grad_content).expect("content"),
            in_space: space,
            out_space: space,
        };
        let mass_sig = AtomSignature {
            name: names.1.to_string(),
            content: ConventionRef::new("feec:m1:v1").expect("content"),
            in_space: space,
            out_space: space,
        };
        let (grad, _mass) = if swap_order {
            let m = system.register_atom(mass_sig);
            let g = system.register_atom(grad_sig);
            (g, m)
        } else {
            let g = system.register_atom(grad_sig);
            let m = system.register_atom(mass_sig);
            (g, m)
        };
        system
            .add_equation(BlockEquation {
                name: "apply-grad".to_string(),
                target: velocity,
                rhs: SystemExpr::Apply {
                    atom: grad,
                    arg: Box::new(SystemExpr::FieldRef(velocity)),
                },
            })
            .expect("equation");
        system.admit().expect("admits").identity()
    };

    let baseline = build(false, "feec:d0:v1", ("d0", "m1"));
    assert_eq!(
        baseline,
        build(true, "feec:d0:v1", ("d0", "m1")),
        "atom declaration order must not move identity"
    );
    assert_eq!(
        baseline,
        build(false, "feec:d0:v1", ("gradient", "mass")),
        "atom display names must not move identity"
    );
    assert_ne!(
        baseline,
        build(false, "feec:d0:v2", ("d0", "m1")),
        "the atom content reference is identity-bearing"
    );

    // Byte-identical atom payloads refuse as ambiguous.
    let mut system = SystemDef::new();
    let space = Space {
        degree: 1,
        n: 64,
        dims: VELOCITY,
    };
    for name in ["first", "second"] {
        system.register_atom(AtomSignature {
            name: name.to_string(),
            content: ConventionRef::new("feec:d0:v1").expect("content"),
            in_space: space,
            out_space: space,
        });
    }
    let ambiguous = system.admit();
    assert!(
        matches!(
            ambiguous,
            Err(SystemTypeError::IndistinguishableAtoms { ref first, ref second })
                if first == "first" && second == "second"
        ),
        "duplicate atom payloads must refuse, got {:?}",
        ambiguous.map(|admitted| admitted.identity())
    );
}
