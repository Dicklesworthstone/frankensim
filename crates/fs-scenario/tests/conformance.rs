//! fs-scenario conformance suite (the tfz.1 bead). Acceptance: scenario
//! values round-trip IR ↔ memory ↔ ledger; compatibility checks catch
//! seeded violations with structured fixes; ensembles reproduce bitwise
//! from seed; Kanai–Tajimi realizations match the target spectrum
//! statistically; frame transforms obey G0 composition laws; unit
//! coherence holds through non-SI spellings (G3).

use fs_blake3::ContentHash;
use fs_ga::{Motor, Point, Quat, Vec3};
use fs_qty::{Dims, QtyAny};
use fs_scenario::{
    BcKind, BcValue, BoundaryCondition, ChebProfile, Combination, Compat, ContactLaw, ContactModel,
    Environment, Frame, FrameId, FrameMotion, Interp, LoadCase, Physics, RealizationBudget,
    Scenario, SpectrumModel, StochasticEnsemble, TimeSignal,
    ir::{
        FiveToSixRule, IrParseBudget, LEGACY_SCENARIO_IR_VERSION, SCENARIO_IR_VERSION, parse_ir,
        parse_ir_with_budget, write_ir,
    },
};

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-scenario/conformance\",\"case\":\"{case}\",\"verdict\":\"pass\",\
         \"detail\":\"{detail}\"}}"
    );
}

const MASS_FLOW: Dims = Dims([0, 1, -1, 0, 0, 0]);
const PRESSURE: Dims = Dims([-1, 1, -2, 0, 0, 0]);
const RATE: Dims = Dims([0, 0, -1, 0, 0, 0]);
const TIME: Dims = Dims([0, 0, 1, 0, 0, 0]);
const LEGACY_MINIMAL_IR: &str = concat!(
    "(scenario \"legacy-v1\" 7 (environment ",
    "(qty 0 1 0 -2 0 0) (qty 0 1 0 -2 0 0) (qty -9.80665 1 0 -2 0 0) ",
    "(qty 293.15 0 0 0 1 0) (qty 101325 -1 1 -2 0 0)) ",
    "(frames) (bcs) (cases) (combos) (ensembles) (contacts))"
);

/// The fixture's ensemble battery (all three model families).
fn fixture_ensembles() -> Vec<StochasticEnsemble> {
    vec![
        StochasticEnsemble {
            name: "gm".to_string(),
            seed: 42,
            members: 8,
            duration: QtyAny::new(20.48, TIME),
            dt: QtyAny::new(0.04, TIME),
            model: SpectrumModel::KanaiTajimi {
                s0: 0.011,
                omega_g: QtyAny::new(12.5, RATE),
                zeta_g: 0.6,
            },
        },
        StochasticEnsemble {
            name: "gusts".to_string(),
            seed: 42,
            members: 4,
            duration: QtyAny::new(10.0, TIME),
            dt: QtyAny::new(0.05, TIME),
            model: SpectrumModel::Dryden {
                sigma: QtyAny::new(1.8, Dims([1, 0, -1, 0, 0, 0])),
                length_scale: QtyAny::new(200.0, Dims([1, 0, 0, 0, 0, 0])),
                mean_speed: QtyAny::new(12.0, Dims([1, 0, -1, 0, 0, 0])),
            },
        },
        StochasticEnsemble {
            name: "rheology".to_string(),
            seed: 42,
            members: 16,
            duration: QtyAny::new(1.0, TIME),
            dt: QtyAny::new(1.0, TIME),
            model: SpectrumModel::CarreauBand {
                eta_zero: [
                    QtyAny::new(8.0, Dims([-1, 1, -1, 0, 0, 0])),
                    QtyAny::new(15.0, Dims([-1, 1, -1, 0, 0, 0])),
                ],
                eta_inf: [
                    QtyAny::new(0.08, Dims([-1, 1, -1, 0, 0, 0])),
                    QtyAny::new(0.15, Dims([-1, 1, -1, 0, 0, 0])),
                ],
                lambda: [QtyAny::new(0.5, TIME), QtyAny::new(2.0, TIME)],
                n: [0.3, 0.6],
            },
        },
    ]
}

/// The fixture's base boundary conditions (all value kinds).
fn fixture_base_bcs() -> Vec<BoundaryCondition> {
    let mut s = Scenario::new("tmp", 0, Environment::earth_lab());
    s.base_bcs.push(BoundaryCondition {
        region: "inlet".to_string(),
        physics: Physics::IncompressibleFlow,
        kind: BcKind::MassFlowInlet,
        value: Some(BcValue::Signal(TimeSignal::Table {
            times: vec![0.0, 1.0, 2.5],
            values: vec![0.0, 0.75, 0.75],
            dims: MASS_FLOW,
            interp: Interp::Linear,
        })),
        compatibility: Some(Compat::Incompressible),
        frame: 1,
    });
    s.base_bcs.push(BoundaryCondition {
        region: "spout".to_string(),
        physics: Physics::IncompressibleFlow,
        kind: BcKind::PressureOutlet,
        value: Some(BcValue::Uniform(QtyAny::new(101_325.0, PRESSURE))),
        compatibility: None,
        frame: 0,
    });
    s.base_bcs.push(BoundaryCondition {
        region: "walls".to_string(),
        physics: Physics::IncompressibleFlow,
        kind: BcKind::WallNoSlip,
        value: None,
        compatibility: None,
        frame: 0,
    });
    s.base_bcs.push(BoundaryCondition {
        region: "inlet".to_string(),
        physics: Physics::Thermal,
        kind: BcKind::Dirichlet,
        value: Some(BcValue::Profile(ChebProfile {
            cheb: fs_cheb::Cheb1::build(&|x: f64| 350.0 + 5.0 * x * x, 0.0, 1.0, 8),
            dims: Dims([0, 0, 0, 1, 0, 0]),
        })),
        compatibility: None,
        frame: 0,
    });
    s.base_bcs
}

/// The vessel-pour flavored fixture: every construct exercised at once.
fn rich_scenario() -> Scenario {
    let mut s = Scenario::new("spout-pour", 42, Environment::earth_lab());
    let tilt_end = fs_qty::parse::parse_qty("65deg").expect("fs-qty parses degrees");
    s.frames.add(Frame {
        id: FrameId(1),
        name: "vessel".to_string(),
        parent: FrameId(0),
        motion: FrameMotion::Tilt {
            axis: [0.0, 1.0, 0.0],
            center: Vec3::new(0.2, 0.0, 0.35),
            angle: TimeSignal::Ramp {
                t_start: 0.0,
                t_end: 3.0,
                from: QtyAny::dimensionless(0.0),
                to: tilt_end,
            },
        },
    });
    s.frames.add(Frame {
        id: FrameId(2),
        name: "stirrer".to_string(),
        parent: FrameId(1),
        motion: FrameMotion::Rotating {
            axis: [0.0, 0.0, 1.0],
            center: Vec3::new(0.0, 0.0, 0.0),
            rate: QtyAny::new(3.5, RATE),
        },
    });
    s.frames.add(Frame {
        id: FrameId(3),
        name: "gauge".to_string(),
        parent: FrameId(2),
        motion: FrameMotion::Fixed {
            orientation: Quat::from_axis_angle(Vec3::new(1.0, 0.0, 0.0), 0.25),
            translation: Vec3::new(0.0, 0.05, -0.125),
        },
    });
    s.base_bcs = fixture_base_bcs();
    s.cases.push(LoadCase {
        name: "dead".to_string(),
        bcs: vec![BoundaryCondition {
            region: "base".to_string(),
            physics: Physics::Elasticity,
            kind: BcKind::Traction,
            value: Some(BcValue::Uniform(QtyAny::new(-2.4e4, PRESSURE))),
            compatibility: None,
            frame: 0,
        }],
    });
    s.cases.push(LoadCase {
        name: "live".to_string(),
        bcs: vec![BoundaryCondition {
            region: "base".to_string(),
            physics: Physics::Elasticity,
            kind: BcKind::Traction,
            value: Some(BcValue::Uniform(QtyAny::new(-1.05e4, PRESSURE))),
            compatibility: None,
            frame: 0,
        }],
    });
    s.combinations.push(Combination {
        name: "1.2D+1.6L".to_string(),
        terms: vec![("dead".to_string(), 1.2), ("live".to_string(), 1.6)],
    });
    s.ensembles = fixture_ensembles();
    s.contacts.push(ContactLaw {
        region_a: "lid".to_string(),
        region_b: "rim".to_string(),
        model: ContactModel::Coulomb {
            mu_s: 0.5,
            mu_k: 0.38,
        },
    });
    s
}

#[test]
fn sc_001_ir_round_trip_memory_and_ledger() {
    let s = rich_scenario();
    assert_eq!(s.validate(), Vec::new(), "the fixture must be admissible");
    let text = write_ir(&s);
    assert!(
        text.starts_with("(scenario :version 2 "),
        "canonical scenario IR must carry its semantic version"
    );
    assert_eq!(text, write_ir(&s), "canonical IR must be byte-stable");
    let back = parse_ir(&text).expect("canonical IR parses");
    assert_eq!(back.source_version(), SCENARIO_IR_VERSION);
    assert_eq!(back.migration(), None);
    assert_eq!(back.scenario(), &s, "IR round trip must be lossless");
    // Ledger round trip (dev-dependency: the L6 integration exercised
    // from tests, keeping the L3 runtime dependency graph clean).
    let dir = std::env::temp_dir().join(format!("fs-scenario-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("mkdir");
    let db = dir.join("scenario.led");
    let ledger = fs_ledger::Ledger::open(db.to_str().expect("utf8 path")).expect("open ledger");
    let receipt = ledger
        .put_artifact("scenario-ir", text.as_bytes(), None)
        .expect("store scenario IR");
    let bytes = ledger
        .get_artifact(&receipt.hash)
        .expect("fetch")
        .expect("present");
    let from_ledger = parse_ir(std::str::from_utf8(&bytes).expect("utf8")).expect("parses");
    assert_eq!(
        from_ledger.scenario(),
        &s,
        "ledger round trip must be lossless"
    );
    let _ = std::fs::remove_dir_all(&dir);
    verdict(
        "sc-001",
        "scenario == parse(write(scenario)) in memory and through a ledger artifact",
    );
}

#[test]
fn sc_001a_legacy_five_base_ir_decodes_with_receipt() {
    let implicit = parse_ir(LEGACY_MINIMAL_IR).expect("implicit historical v1 parses");
    assert_eq!(implicit.source_version(), LEGACY_SCENARIO_IR_VERSION);
    let implicit_receipt = implicit
        .migration()
        .expect("legacy migration receipt is mandatory");
    assert_eq!(
        implicit_receipt.source_version(),
        LEGACY_SCENARIO_IR_VERSION
    );
    assert_eq!(implicit_receipt.target_version(), SCENARIO_IR_VERSION);
    assert_eq!(implicit_receipt.source_width(), 5);
    assert_eq!(implicit_receipt.target_width(), 6);
    assert_eq!(implicit_receipt.rule(), FiveToSixRule::AppendMoleZero);
    assert_eq!(
        implicit_receipt.old_hash(),
        ContentHash::from_hex("a9d2537cc4717dae81a760d18f11f5f4876b5f8fc77535156f29f44004f37a22")
            .expect("pinned implicit-v1 hash")
    );
    assert_eq!(
        implicit.scenario().environment.ambient_pressure.dims,
        Dims([-1, 1, -2, 0, 0, 0]),
        "the immutable crosswalk appends mol=0"
    );

    let explicit_text = LEGACY_MINIMAL_IR.replacen("(scenario ", "(scenario :version 1 ", 1);
    let explicit = parse_ir(&explicit_text).expect("explicit historical v1 parses");
    assert_eq!(explicit.source_version(), LEGACY_SCENARIO_IR_VERSION);
    let explicit_receipt = explicit
        .migration()
        .expect("explicit v1 also requires a receipt");
    assert_eq!(
        explicit_receipt.old_hash(),
        ContentHash::from_hex("1ca63b192cbac63086a4447cf9adfd5f87b081505aafb17107a0538e42428415")
            .expect("pinned explicit-v1 hash")
    );
    assert_ne!(explicit_receipt.old_hash(), implicit_receipt.old_hash());
    assert_eq!(explicit.scenario(), implicit.scenario());

    let canonical = write_ir(implicit.scenario());
    assert!(canonical.starts_with("(scenario :version 2 "));
    assert_ne!(
        canonical, LEGACY_MINIMAL_IR,
        "legacy bytes are never rewritten in place"
    );
    let pinned_canonical =
        ContentHash::from_hex("f74494d0152c0587c611cdab5dcc2f96baca05bfe50fa13a18f4a64ada15bffc")
            .expect("pinned canonical-v2 hash");
    assert_eq!(implicit_receipt.new_hash(), pinned_canonical);
    assert_eq!(explicit_receipt.new_hash(), pinned_canonical);
    assert!(implicit_receipt.verifies(LEGACY_MINIMAL_IR.as_bytes(), canonical.as_bytes()));
    assert!(explicit_receipt.verifies(explicit_text.as_bytes(), canonical.as_bytes()));
    assert!(!implicit_receipt.verifies(b"tampered", canonical.as_bytes()));
    let reparsed = parse_ir(&canonical).expect("migrated canonical v2 parses");
    assert_eq!(reparsed.source_version(), SCENARIO_IR_VERSION);
    assert_eq!(reparsed.migration(), None);
    assert_eq!(reparsed.scenario(), implicit.scenario());
}

#[test]
fn sc_001c_version_and_dimension_arity_fail_closed() {
    let explicit_v1 = LEGACY_MINIMAL_IR.replacen("(scenario ", "(scenario :version 1 ", 1);
    let v1_six = explicit_v1.replacen("(qty 0 1 0 -2 0 0)", "(qty 0 1 0 -2 0 0 0)", 1);
    assert_ne!(v1_six, explicit_v1, "the v1 fixture mutation must apply");
    let error = parse_ir(&v1_six).expect_err("v1 must reject six dimension exponents");
    assert!(
        error.to_string().contains("5 exponents"),
        "unexpected v1 arity diagnosis: {error}"
    );

    let legacy = parse_ir(LEGACY_MINIMAL_IR).expect("legacy fixture parses");
    let canonical = write_ir(legacy.scenario());
    let v2_five = canonical.replacen("(qty 0 1 0 -2 0 0 0)", "(qty 0 1 0 -2 0 0)", 1);
    assert_ne!(v2_five, canonical, "the v2 fixture mutation must apply");
    let error = parse_ir(&v2_five).expect_err("v2 must reject five dimension exponents");
    assert!(
        error.to_string().contains("6 exponents"),
        "unexpected v2 arity diagnosis: {error}"
    );

    let current_tag = format!(":version {SCENARIO_IR_VERSION}");
    for unsupported in [0, 3, u32::MAX] {
        let bad = canonical.replacen(&current_tag, &format!(":version {unsupported}"), 1);
        assert_ne!(bad, canonical, "the version mutation must apply");
        let error = parse_ir(&bad).expect_err("unsupported versions must fail closed");
        assert!(
            error
                .to_string()
                .contains("unsupported scenario IR version"),
            "unexpected version diagnosis for {unsupported}: {error}"
        );
    }
}

#[test]
fn sc_001d_legacy_crosswalk_reaches_nested_quantities_and_dims() {
    const LEGACY_NESTED: &str = concat!(
        "(scenario :version 1 \"nested-v1\" 7 ",
        "(environment ",
        "(qty 0 1 0 -2 0 0) (qty 0 1 0 -2 0 0) (qty -9.8 1 0 -2 0 0) ",
        "(qty 293 0 0 0 1 0) (qty 101325 -1 1 -2 0 0)) ",
        "(frames (frame 1 \"rotor\" 0 ",
        "(rotating (vec 0 0 1) (vec 0 0 0) (qty 2 0 0 -1 0 0)))) ",
        "(bcs (bc \"feed\" incompressible-flow mass-flow-inlet 1 ",
        "(signal (table linear (dims 0 1 -1 0 0) (times 0 1) (values 0 1))) ",
        "incompressible)) ",
        "(cases) (combos) ",
        "(ensembles (ensemble \"gust\" 9 2 ",
        "(qty 1 0 0 1 0 0) (qty 0.1 0 0 1 0 0) ",
        "(dryden (qty 1 1 0 -1 0 0) (qty 2 1 0 0 0 0) ",
        "(qty 3 1 0 -1 0 0)))) ",
        "(contacts))"
    );

    let decoded = parse_ir(LEGACY_NESTED).expect("nested explicit v1 parses");
    assert_eq!(decoded.source_version(), LEGACY_SCENARIO_IR_VERSION);
    let receipt = decoded
        .migration()
        .expect("nested legacy input requires a receipt");
    let scenario = decoded.scenario();

    let FrameMotion::Rotating { rate, .. } = &scenario.frames.frames[0].motion else {
        panic!("the nested fixture must retain its rotating frame");
    };
    assert_eq!(rate.dims.0[5], 0);

    let Some(BcValue::Signal(TimeSignal::Table { dims, .. })) = scenario.base_bcs[0].value.as_ref()
    else {
        panic!("the nested fixture must retain its table signal");
    };
    assert_eq!(dims.0[5], 0);

    let ensemble = &scenario.ensembles[0];
    let SpectrumModel::Dryden {
        sigma,
        length_scale,
        mean_speed,
    } = &ensemble.model
    else {
        panic!("the nested fixture must retain its Dryden model");
    };
    for quantity in [
        &ensemble.duration,
        &ensemble.dt,
        sigma,
        length_scale,
        mean_speed,
    ] {
        assert_eq!(quantity.dims.0[5], 0);
    }

    let canonical = write_ir(scenario);
    assert!(receipt.verifies(LEGACY_NESTED.as_bytes(), canonical.as_bytes()));
}

#[test]
fn sc_001e_v2_preserves_nonzero_mole_exponents_in_nested_paths() {
    let mut scenario = rich_scenario();

    let FrameMotion::Rotating { rate, .. } = &mut scenario.frames.frames[1].motion else {
        panic!("rich fixture must contain its rotating frame");
    };
    rate.dims.0[5] = 1;

    let Some(BcValue::Signal(TimeSignal::Table { dims, .. })) = scenario.base_bcs[0].value.as_mut()
    else {
        panic!("rich fixture must contain its table signal");
    };
    dims.0[5] = -1;

    let SpectrumModel::Dryden { sigma, .. } = &mut scenario.ensembles[1].model else {
        panic!("rich fixture must contain its Dryden ensemble");
    };
    sigma.dims.0[5] = 2;

    let text = write_ir(&scenario);
    let decoded = parse_ir(&text).expect("canonical v2 with mole exponents parses");
    assert_eq!(decoded.source_version(), SCENARIO_IR_VERSION);
    assert_eq!(decoded.migration(), None);
    assert_eq!(decoded.scenario(), &scenario);

    let FrameMotion::Rotating { rate, .. } = &decoded.scenario().frames.frames[1].motion else {
        panic!("decoded fixture must retain its rotating frame");
    };
    assert_eq!(rate.dims.0[5], 1);
    let Some(BcValue::Signal(TimeSignal::Table { dims, .. })) =
        decoded.scenario().base_bcs[0].value.as_ref()
    else {
        panic!("decoded fixture must retain its table signal");
    };
    assert_eq!(dims.0[5], -1);
    let SpectrumModel::Dryden { sigma, .. } = &decoded.scenario().ensembles[1].model else {
        panic!("decoded fixture must retain its Dryden ensemble");
    };
    assert_eq!(sigma.dims.0[5], 2);
}

#[test]
fn sc_002_compatibility_checks_catch_seeded_violations() {
    // (a) Incompressible inlets with no outlet and unbalanced flux.
    let mut s = Scenario::new("bad-flux", 1, Environment::earth_lab());
    s.base_bcs.push(BoundaryCondition {
        region: "in".to_string(),
        physics: Physics::IncompressibleFlow,
        kind: BcKind::MassFlowInlet,
        value: Some(BcValue::Uniform(QtyAny::new(0.9, MASS_FLOW))),
        compatibility: Some(Compat::Incompressible),
        frame: 0,
    });
    let violations = s.validate();
    let flux = violations
        .iter()
        .find(|v| v.code == "flux-imbalance")
        .expect("net-flux violation must be caught");
    assert!(
        flux.fix.contains("pressure outlet"),
        "fix must be actionable"
    );
    // Adding the outlet repairs it.
    s.base_bcs.push(BoundaryCondition {
        region: "out".to_string(),
        physics: Physics::IncompressibleFlow,
        kind: BcKind::PressureOutlet,
        value: Some(BcValue::Uniform(QtyAny::new(0.0, PRESSURE))),
        compatibility: None,
        frame: 0,
    });
    assert!(
        s.validate().iter().all(|v| v.code != "flux-imbalance"),
        "outlet must absorb the imbalance"
    );

    // (b) Inconsistent frame chains: a cycle and a dangling parent.
    let mut s2 = Scenario::new("bad-frames", 1, Environment::earth_lab());
    let spin = |id: u32, name: &str, parent: u32| Frame {
        id: FrameId(id),
        name: name.to_string(),
        parent: FrameId(parent),
        motion: FrameMotion::Rotating {
            axis: [0.0, 0.0, 1.0],
            center: Vec3::new(0.0, 0.0, 0.0),
            rate: QtyAny::new(1.0, RATE),
        },
    };
    s2.frames.add(spin(1, "a", 2));
    s2.frames.add(spin(2, "b", 1));
    s2.frames.add(spin(3, "c", 99));
    let v2 = s2.validate();
    assert!(v2.iter().any(|v| v.code == "frame-chain-cyclic"));
    assert!(v2.iter().any(|v| v.code == "frame-parent-missing"));

    // (c) Dimension violations + missing compatibility + bad combos.
    let mut s3 = Scenario::new("bad-dims", 1, Environment::earth_lab());
    s3.base_bcs.push(BoundaryCondition {
        region: "in".to_string(),
        physics: Physics::IncompressibleFlow,
        kind: BcKind::Dirichlet,
        value: Some(BcValue::Uniform(QtyAny::new(
            300.0,
            Dims([0, 0, 0, 1, 0, 0]),
        ))),
        compatibility: None,
        frame: 7, // undefined frame
    });
    s3.base_bcs.push(BoundaryCondition {
        region: "in2".to_string(),
        physics: Physics::IncompressibleFlow,
        kind: BcKind::MassFlowInlet,
        value: Some(BcValue::Uniform(QtyAny::new(0.1, MASS_FLOW))),
        compatibility: None, // flux inlet must declare a regime
        frame: 0,
    });
    s3.combinations.push(Combination {
        name: "c".to_string(),
        terms: vec![("ghost".to_string(), 1.4)],
    });
    s3.contacts.push(ContactLaw {
        region_a: "a".to_string(),
        region_b: "b".to_string(),
        model: ContactModel::Coulomb {
            mu_s: 0.2,
            mu_k: 0.5, // kinetic above static
        },
    });
    let v3 = s3.validate();
    for code in [
        "bc-dims",
        "bc-frame-missing",
        "bc-compat-missing",
        "combo-case-missing",
        "contact-coulomb-range",
    ] {
        let hit = v3.iter().find(|v| v.code == code);
        assert!(hit.is_some(), "expected violation {code}");
        assert!(
            !hit.expect("present").fix.is_empty(),
            "{code} must carry a fix"
        );
    }
    verdict(
        "sc-002",
        "flux imbalance, cyclic/dangling frames, dims, compat, combos, contacts all caught \
         with structured fixes",
    );
}

#[test]
fn sc_002a_ramp_and_table_from_zero_are_checked_after_t_zero() {
    let signals = [
        (
            "ramp",
            TimeSignal::Ramp {
                t_start: 0.0,
                t_end: 2.0,
                from: QtyAny::new(0.0, MASS_FLOW),
                to: QtyAny::new(0.75, MASS_FLOW),
            },
        ),
        (
            "table",
            TimeSignal::Table {
                times: vec![0.0, 1.0, 2.0],
                values: vec![0.0, 0.75, 0.75],
                dims: MASS_FLOW,
                interp: Interp::Linear,
            },
        ),
    ];
    for (name, signal) in signals {
        let mut scenario = Scenario::new(name, 1, Environment::earth_lab());
        scenario.base_bcs.push(BoundaryCondition {
            region: format!("{name}-inlet"),
            physics: Physics::IncompressibleFlow,
            kind: BcKind::MassFlowInlet,
            value: Some(BcValue::Signal(signal)),
            compatibility: Some(Compat::Incompressible),
            frame: 0,
        });
        let violations = scenario.validate();
        let flux = violations
            .iter()
            .find(|violation| violation.code == "flux-imbalance")
            .unwrap_or_else(|| panic!("{name} from zero must fail after t=0: {violations:#?}"));
        assert!(
            flux.what.contains("at t=") && !flux.what.contains("t=0.000000e0"),
            "the imbalance must be tied to a nonzero exact breakpoint: {flux:?}"
        );
    }
    verdict(
        "sc-002a",
        "ramp/table mass flows that vanish at t=0 are rejected at their deterministic nonzero breakpoints without a pressure outlet",
    );
}

#[test]
fn sc_002b_mixed_signal_flux_grid_is_deterministic() {
    let mut scenario = Scenario::new("mixed-flux-grid", 1, Environment::earth_lab());
    let signals = [
        (
            "ramp",
            TimeSignal::Ramp {
                t_start: 2.0,
                t_end: 6.0,
                from: QtyAny::new(0.0, MASS_FLOW),
                to: QtyAny::new(1.0, MASS_FLOW),
            },
        ),
        (
            "table",
            TimeSignal::Table {
                times: vec![2.0, 6.0],
                values: vec![0.0, -1.0],
                dims: MASS_FLOW,
                interp: Interp::Linear,
            },
        ),
        (
            "smooth",
            TimeSignal::Chebfun(ChebProfile {
                // On [2, 6], this is 1 - x_ref^2: zero at both exact
                // breakpoints and positive only inside the declared domain.
                cheb: fs_cheb::Cheb1::from_coeffs(2.0, 6.0, vec![0.5, 0.0, -0.5]),
                dims: MASS_FLOW,
            }),
        ),
    ];
    for (name, signal) in signals {
        scenario.base_bcs.push(BoundaryCondition {
            region: format!("{name}-inlet"),
            physics: Physics::IncompressibleFlow,
            kind: BcKind::MassFlowInlet,
            value: Some(BcValue::Signal(signal)),
            compatibility: Some(Compat::Incompressible),
            frame: 0,
        });
    }

    let first = scenario.validate();
    let replay = scenario.validate();
    assert_eq!(
        first, replay,
        "mixed-signal checkpoint traversal must replay"
    );
    let flux: Vec<_> = first
        .iter()
        .filter(|violation| violation.code == "flux-imbalance")
        .collect();
    assert_eq!(
        flux.len(),
        1,
        "unexpected mixed-signal findings: {first:#?}"
    );
    assert!(
        flux[0].what.contains("at t=") && !flux[0].what.contains("t=0.000000e0"),
        "the smooth interior imbalance must be found on its declared-domain grid: {:?}",
        flux[0]
    );
    verdict(
        "sc-002b",
        "ramp/table breakpoints plus the bounded Chebfun domain grid replay deterministically and expose an interior-only mixed-signal imbalance",
    );
}

#[test]
fn sc_001b_non_ascii_names_round_trip() {
    // Regression: the IR string parser decoded bytes as Latin-1 (`push(byte as
    // char)`), splitting every multi-byte UTF-8 code point, so any non-ASCII
    // name (scenario/region/frame) silently corrupted on parse — violating the
    // round-trip losslessness invariant. The writer emits proper UTF-8.
    let mut s = Scenario::new("Kármán–pour café ✓", 7, Environment::earth_lab());
    s.frames.add(Frame {
        id: FrameId(1),
        name: "vase — 花瓶".to_string(),
        parent: FrameId(0),
        motion: FrameMotion::Fixed {
            orientation: Quat::from_axis_angle(Vec3::new(1.0, 0.0, 0.0), 0.0),
            translation: Vec3::new(0.0, 0.0, 0.0),
        },
    });
    s.base_bcs.push(BoundaryCondition {
        region: "über-inlet — 入口".to_string(),
        physics: Physics::IncompressibleFlow,
        kind: BcKind::WallNoSlip,
        value: None,
        compatibility: None,
        frame: 0,
    });
    assert_eq!(
        s.validate(),
        Vec::new(),
        "non-ASCII exact identities are admissible when nonempty and unique"
    );
    let text = write_ir(&s);
    let back = parse_ir(&text).expect("non-ASCII IR parses");
    assert_eq!(
        back.scenario(),
        &s,
        "non-ASCII names must round-trip losslessly"
    );
    verdict(
        "sc-001b",
        "scenario/frame/region names with non-ASCII (é, ü, –, ✓, CJK) round-trip losslessly",
    );
}

#[test]
fn sc_001c_string_escapes_are_canonical_and_unambiguous() {
    let mut scenario = Scenario::new(r#"quoted "name" and \ path"#, 7, Environment::earth_lab());
    scenario.frames.add(Frame {
        id: FrameId(1),
        name: r#"frame \ "one""#.to_string(),
        parent: FrameId(0),
        motion: FrameMotion::Fixed {
            orientation: Quat::identity(),
            translation: Vec3::new(0.0, 0.0, 0.0),
        },
    });
    let canonical = write_ir(&scenario);
    let decoded = parse_ir(&canonical).expect("writer escapes must parse");
    assert_eq!(decoded.scenario(), &scenario);
    assert_eq!(write_ir(decoded.scenario()), canonical);

    let unknown_escape = canonical.replacen("quoted", r"quo\qted", 1);
    assert_ne!(unknown_escape, canonical, "escape mutation must apply");
    let error = parse_ir(&unknown_escape).expect_err("unknown escapes are non-canonical");
    assert!(
        error.to_string().contains("unsupported string escape"),
        "unexpected diagnostic: {error}"
    );
    verdict(
        "sc-001c",
        "IR accepts exactly the writer's quote/backslash escapes and rejects alias encodings",
    );
}

#[test]
fn sc_001d_canonical_ir_normalizes_physically_irrelevant_signed_zero() {
    let positive_zero = Scenario::new("zero", 11, Environment::earth_lab());
    let mut negative_zero = positive_zero.clone();
    negative_zero.environment.gravity[0].value = -0.0;

    assert_eq!(
        positive_zero, negative_zero,
        "derived semantic equality aliases signed zero"
    );
    assert_eq!(
        write_ir(&positive_zero),
        write_ir(&negative_zero),
        "equal scenarios must have one canonical byte identity"
    );
    let canonical = write_ir(&negative_zero);
    assert_eq!(
        write_ir(
            parse_ir(&canonical)
                .expect("canonical signed-zero form parses")
                .scenario()
        ),
        canonical,
        "canonical writer must be idempotent"
    );
    verdict(
        "sc-001d",
        "semantic float equality and canonical scenario identity agree for signed zero",
    );
}

/// Replace the balanced sub-form beginning at `head` (e.g. `"(fixed"`) with
/// just its head token (`"(fixed)"`) — a truncated, under-arity form.
fn truncate_form(ir: &str, head: &str) -> String {
    let start = ir.find(head).expect("truncation target form present in IR");
    let mut depth = 0i32;
    let mut end = ir.len();
    for (i, ch) in ir[start..].char_indices() {
        match ch {
            '(' => depth += 1,
            ')' => {
                depth -= 1;
                if depth == 0 {
                    end = start + i + 1; // ')' is one byte
                    break;
                }
            }
            _ => {}
        }
    }
    let mut out = String::with_capacity(ir.len());
    out.push_str(&ir[..start]);
    out.push_str(head); // "(fixed"
    out.push(')'); // -> "(fixed)"
    out.push_str(&ir[end..]);
    out
}

#[test]
fn sc_002b_truncated_ir_forms_error_not_panic() {
    // Regression: several parser arms indexed `rest[0]`/`rest[0..2]` without an
    // arity check, so a truncated form panicked (index out of bounds) instead
    // of returning `ScenarioError::Parse`. `parse_ir` is public and documented
    // to return `Result` on malformed input. Truncate each sub-form in a valid
    // scenario's IR to just its head and confirm a clean error, not a panic.
    let ir = write_ir(&rich_scenario());
    for head in [
        "(fixed",
        "(rotating",
        "(tilt",
        "(dryden",
        "(kanai-tajimi",
        "(uniform",
    ] {
        let broken = truncate_form(&ir, head);
        assert!(
            parse_ir(&broken).is_err(),
            "truncated {head}) form must return an error, not panic"
        );
    }
    verdict(
        "sc-002b",
        "truncated frame/model/bc IR forms return Parse errors, never an index-OOB panic",
    );
}

#[test]
fn sc_002c_ensemble_rejects_nan_producing_spectra() {
    // Regression: KanaiTajimi `check` validated S0/zeta_g but NOT omega_g (the
    // psd divisor `r = ω/ω_g`), and Dryden had no positivity check at all
    // (`mean_speed` is the divisor `x = l·ω/v`). A zero divisor makes every
    // realization NaN/inf, yet `validate()` admitted the ensemble. Both now
    // fail closed.
    let mut s = Scenario::new("nan-ensembles", 1, Environment::earth_lab());
    s.ensembles.push(StochasticEnsemble {
        name: "kt-zero-omega".to_string(),
        seed: 1,
        members: 2,
        duration: QtyAny::new(1.0, TIME),
        dt: QtyAny::new(0.1, TIME),
        model: SpectrumModel::KanaiTajimi {
            s0: 0.02,
            omega_g: QtyAny::new(0.0, RATE),
            zeta_g: 0.5,
        },
    });
    s.ensembles.push(StochasticEnsemble {
        name: "dryden-zero-speed".to_string(),
        seed: 2,
        members: 2,
        duration: QtyAny::new(1.0, TIME),
        dt: QtyAny::new(0.1, TIME),
        model: SpectrumModel::Dryden {
            sigma: QtyAny::new(1.0, Dims([1, 0, -1, 0, 0, 0])),
            length_scale: QtyAny::new(10.0, Dims([1, 0, 0, 0, 0, 0])),
            mean_speed: QtyAny::new(0.0, Dims([1, 0, -1, 0, 0, 0])),
        },
    });
    let v = s.validate();
    assert!(
        v.iter().any(|x| x.code == "ensemble-kt-params"),
        "zero omega_g (KT psd divisor) must be rejected"
    );
    assert!(
        v.iter().any(|x| x.code == "ensemble-dryden-params"),
        "zero mean_speed (Dryden psd divisor) must be rejected"
    );
    verdict(
        "sc-002c",
        "ensembles whose PSD divisor (omega_g / mean_speed) is zero fail closed, not admitted as NaN",
    );
}

#[test]
fn sc_002d_signal_failure_is_not_misreported_as_a_frame_cycle() {
    // Regression: the cycle check was `world_pose(...).is_err()`, which returns
    // Err both for a real cycle AND for a `local_pose` evaluation failure. An
    // ACYCLIC frame whose Tilt angle is an (independently flagged) empty table
    // therefore got a SPURIOUS `frame-chain-cyclic` violation on top of the
    // real `signal-table-shape` one, misdirecting repair. The check is now
    // structural (no pose evaluation).
    let mut s = Scenario::new("acyclic-with-bad-signal", 1, Environment::earth_lab());
    s.frames.add(Frame {
        id: FrameId(1),
        name: "base".to_string(),
        parent: FrameId(0),
        motion: FrameMotion::Fixed {
            orientation: Quat::from_axis_angle(Vec3::new(1.0, 0.0, 0.0), 0.0),
            translation: Vec3::new(0.0, 0.0, 0.0),
        },
    });
    // Acyclic child of frame 1 whose Tilt angle is an empty table.
    s.frames.add(Frame {
        id: FrameId(2),
        name: "tilter".to_string(),
        parent: FrameId(1),
        motion: FrameMotion::Tilt {
            axis: [0.0, 1.0, 0.0],
            center: Vec3::new(0.0, 0.0, 0.0),
            angle: TimeSignal::Table {
                interp: Interp::Linear,
                dims: Dims([0, 0, 0, 0, 0, 0]),
                times: vec![],
                values: vec![],
            },
        },
    });
    let v = s.validate();
    assert!(
        v.iter().any(|x| x.code == "signal-table-shape"),
        "the genuine empty-table defect must still be reported"
    );
    assert!(
        !v.iter().any(|x| x.code == "frame-chain-cyclic"),
        "an acyclic chain must NOT be flagged cyclic just because a signal failed to evaluate"
    );
    verdict(
        "sc-002d",
        "a local-pose/signal evaluation failure is reported as itself, not misattributed to a cycle",
    );
}

#[test]
fn sc_003_ensembles_reproduce_bitwise_from_seed() {
    let s = rich_scenario();
    for e in &s.ensembles {
        for member in 0..e.members.min(4) {
            let a = e.realize(member).expect("realize");
            let b = e.realize(member).expect("realize again");
            assert_eq!(a.values.len(), b.values.len());
            for (x, y) in a.values.iter().zip(b.values.iter()) {
                assert_eq!(
                    x.to_bits(),
                    y.to_bits(),
                    "{}: member must be bitwise",
                    e.name
                );
            }
            if member > 0 {
                let prev = e.realize(member - 1).expect("prev member");
                assert_ne!(
                    a.values, prev.values,
                    "{}: members must differ from each other",
                    e.name
                );
            }
        }
        // A different seed produces a different realization.
        let mut reseeded = e.clone();
        reseeded.seed ^= 0xDEAD_BEEF;
        assert_ne!(
            e.realize(0).expect("a").values,
            reseeded.realize(0).expect("b").values,
            "{}: seed must matter",
            e.name
        );
    }
    // Carreau draws stay inside their declared bands.
    let rheo = &s.ensembles[2];
    for member in 0..rheo.members {
        let r = rheo.realize(member).expect("carreau");
        assert!(r.values[0] >= 8.0 && r.values[0] <= 15.0, "eta0 in band");
        assert!(r.values[1] >= 0.08 && r.values[1] <= 0.15, "etainf in band");
        assert!(r.values[2] >= 0.5 && r.values[2] <= 2.0, "lambda in band");
        assert!(r.values[3] >= 0.3 && r.values[3] <= 0.6, "n in band");
    }
    verdict(
        "sc-003",
        "KT/Dryden/Carreau members replay bitwise from the complete retained ensemble spec; bands respected",
    );
}

#[test]
fn sc_004_kanai_tajimi_matches_target_spectrum() {
    let members = 48u32;
    let n = 512usize;
    let dt = 0.05f64;
    let ens = StochasticEnsemble {
        name: "kt".to_string(),
        seed: 2026,
        members,
        duration: QtyAny::new(n as f64 * dt, TIME),
        dt: QtyAny::new(dt, TIME),
        model: SpectrumModel::KanaiTajimi {
            s0: 0.02,
            omega_g: QtyAny::new(12.0, RATE),
            zeta_g: 0.55,
        },
    };
    let d_omega = 2.0 * std::f64::consts::PI / (n as f64 * dt);
    // Ensemble-averaged periodogram at a band of harmonics.
    let band: Vec<usize> = (4..40).collect();
    let mut est = vec![0.0f64; band.len()];
    for member in 0..members {
        let r = ens.realize(member).expect("realize");
        for (slot, &k) in est.iter_mut().zip(band.iter()) {
            let (mut re, mut im) = (0.0f64, 0.0f64);
            for (j, &x) in r.values.iter().enumerate() {
                let phase = 2.0 * std::f64::consts::PI * (k * j % n) as f64 / n as f64;
                re += x * phase.cos();
                im -= x * phase.sin();
            }
            *slot += 2.0 * (re * re + im * im) / (n as f64 * n as f64 * d_omega);
        }
    }
    let mut rel_sum = 0.0f64;
    let mut rel_max = 0.0f64;
    for (avg, &k) in est.iter().map(|e| e / f64::from(members)).zip(band.iter()) {
        let target = ens
            .model
            .try_psd(k as f64 * d_omega)
            .expect("validated fixture has a finite PSD");
        let rel = (avg - target).abs() / target;
        rel_sum += rel;
        rel_max = rel_max.max(rel);
        assert!(
            rel < 0.75,
            "periodogram bin k={k} off by {rel:.3} (est {avg:.4e} vs target {target:.4e})"
        );
    }
    let rel_mean = rel_sum / band.len() as f64;
    assert!(
        rel_mean < 0.15,
        "band-averaged spectral mismatch {rel_mean:.3} exceeds statistical tolerance"
    );
    println!(
        "{{\"suite\":\"fs-scenario/conformance\",\"metric\":\"kt_spectrum_match\",\
         \"members\":{members},\"rel_mean\":{rel_mean:.4},\"rel_max\":{rel_max:.4}}}"
    );
    verdict(
        "sc-004",
        "48-member KT ensemble periodogram matches the target PSD across the band",
    );
}

#[test]
fn sc_005_frame_g0_laws_and_moving_frames() {
    let s = rich_scenario();
    let t = 1.5f64;
    // G0 composition: world pose of a chain equals the manual composition
    // of the locals, top-down.
    let f1 = &s.frames.frames[0];
    let f2 = &s.frames.frames[1];
    let f3 = &s.frames.frames[2];
    let manual = fs_scenario::FrameTree::local_pose(f1, t)
        .expect("l1")
        .compose(&fs_scenario::FrameTree::local_pose(f2, t).expect("l2"))
        .compose(&fs_scenario::FrameTree::local_pose(f3, t).expect("l3"));
    let chained = s.frames.world_pose(FrameId(3), t).expect("chain");
    let p = Point::new(0.3, -0.2, 0.7);
    let a = manual.transform_point(p).expect("manual");
    let b = chained.transform_point(p).expect("chained");
    assert!(
        (Vec3::new(a.x - b.x, a.y - b.y, a.z - b.z)).norm() < 1e-12,
        "composition law"
    );
    // The tilt schedule clamps and matches a directly-built motor.
    let tilt_at = |time: f64| -> Motor {
        let angle = match &f1.motion {
            FrameMotion::Tilt { angle, .. } => angle.eval(time).expect("angle").value,
            _ => unreachable!(),
        };
        let c = Vec3::new(0.2, 0.0, 0.35);
        Motor::translator(c.x, c.y, c.z)
            .compose(&Motor::rotor([0.0, 1.0, 0.0], angle))
            .compose(&Motor::translator(-c.x, -c.y, -c.z))
    };
    for time in [0.0, 0.7, 1.5, 3.0, 5.0] {
        let via_tree = s.frames.world_pose(FrameId(1), time).expect("pose");
        let direct = tilt_at(time);
        let pa = via_tree.transform_point(p).expect("a");
        let pb = direct.transform_point(p).expect("b");
        assert!(
            (Vec3::new(pa.x - pb.x, pa.y - pb.y, pa.z - pb.z)).norm() < 1e-12,
            "tilt schedule pose at t={time}"
        );
    }
    // Clamp: past the ramp end the pose freezes at 65°.
    let end = s.frames.world_pose(FrameId(1), 3.0).expect("end");
    let later = s.frames.world_pose(FrameId(1), 30.0).expect("later");
    let pe = end.transform_point(p).expect("e");
    let pl = later.transform_point(p).expect("l");
    assert!(
        (Vec3::new(pe.x - pl.x, pe.y - pl.y, pe.z - pl.z)).norm() < 1e-12,
        "ramp must clamp"
    );
    // A rotating frame leaves its own axis fixed.
    let center_point = Point::new(0.0, 0.0, 0.9);
    let spin = s.frames.frames[1].clone();
    let local = fs_scenario::FrameTree::local_pose(&spin, 2.0).expect("spin");
    let moved = local.transform_point(center_point).expect("fixed axis");
    assert!(
        (Vec3::new(moved.x, moved.y, moved.z - 0.9)).norm() < 1e-12,
        "points on the rotation axis must not move"
    );
    verdict(
        "sc-005",
        "chain composition == manual motor product; tilt ramp matches, clamps; axis fixed",
    );
}

#[test]
fn sc_006_unit_rescaling_coherence() {
    // G3: the same physics expressed through different unit spellings
    // must land on the same SI values and identical canonical IR.
    let deg = fs_qty::parse::parse_qty("65deg").expect("deg");
    let rad = fs_qty::parse::parse_qty("1.1344640137963142rad").expect("rad");
    assert_eq!(deg.dims, rad.dims, "angles agree on dimensions");
    assert!((deg.value - rad.value).abs() < 1e-15, "same SI value");
    let mm = fs_qty::parse::parse_qty("350mm").expect("mm");
    let m = fs_qty::parse::parse_qty("0.35m").expect("m");
    assert_eq!(mm.dims, m.dims);
    assert!((mm.value - m.value).abs() < 1e-15);
    // A scenario whose inputs came through either spelling validates
    // identically (dimension checks see only SI exponents).
    let mut a = Scenario::new("g3", 9, Environment::earth_lab());
    a.base_bcs.push(BoundaryCondition {
        region: "in".to_string(),
        physics: Physics::Thermal,
        kind: BcKind::Dirichlet,
        value: Some(BcValue::Uniform(
            fs_qty::parse::parse_qty("350K").expect("kelvin"),
        )),
        compatibility: None,
        frame: 0,
    });
    assert_eq!(a.validate(), Vec::new());
    let mut b = a.clone();
    b.base_bcs[0].value = Some(BcValue::Uniform(QtyAny::new(
        350.0,
        Dims([0, 0, 0, 1, 0, 0]),
    )));
    assert_eq!(
        write_ir(&a),
        write_ir(&b),
        "spellings converge to one canonical IR"
    );
    verdict(
        "sc-006",
        "deg/rad and mm/m spellings converge; validation + canonical IR are spelling-blind",
    );
}

#[test]
fn sc_007_ir_budgets_and_chebyshev_inputs_fail_closed() {
    let minimal = Scenario::new("budget", 1, Environment::earth_lab());
    let canonical = write_ir(&minimal);
    let defaults = IrParseBudget::default();

    // Every configurable dimension has a demonstrated exact-pass / one-less
    // refusal boundary on this retained canonical fixture.
    let exact_bytes = IrParseBudget {
        max_bytes: canonical.len(),
        ..defaults
    };
    parse_ir_with_budget(&canonical, exact_bytes).expect("exact byte budget admits fixture");
    let error = parse_ir_with_budget(
        &canonical,
        IrParseBudget {
            max_bytes: canonical.len() - 1,
            ..defaults
        },
    )
    .expect_err("byte budget minus one must refuse");
    assert!(error.to_string().contains("byte budget"));

    for (dimension, exact, below) in [
        ("depth", 4usize, 3usize),
        ("atom", 11usize, 10usize),
        ("list", 12usize, 11usize),
    ] {
        let exact_budget = match dimension {
            "depth" => IrParseBudget {
                max_depth: exact,
                ..defaults
            },
            "atom" => IrParseBudget {
                max_atom_bytes: exact,
                ..defaults
            },
            "list" => IrParseBudget {
                max_list_items: exact,
                ..defaults
            },
            _ => unreachable!(),
        };
        parse_ir_with_budget(&canonical, exact_budget)
            .unwrap_or_else(|error| panic!("exact {dimension} budget must pass: {error}"));
        let below_budget = match dimension {
            "depth" => IrParseBudget {
                max_depth: below,
                ..defaults
            },
            "atom" => IrParseBudget {
                max_atom_bytes: below,
                ..defaults
            },
            "list" => IrParseBudget {
                max_list_items: below,
                ..defaults
            },
            _ => unreachable!(),
        };
        let error = parse_ir_with_budget(&canonical, below_budget)
            .expect_err("budget minus one must refuse");
        assert!(
            error.to_string().contains("budget"),
            "{dimension} refusal must identify its budget: {error}"
        );
    }

    let exact_nodes = (1usize..=256)
        .find(|&max_nodes| {
            parse_ir_with_budget(
                &canonical,
                IrParseBudget {
                    max_nodes,
                    ..defaults
                },
            )
            .is_ok()
        })
        .expect("minimal fixture has a small finite node count");
    assert_eq!(
        exact_nodes, 65,
        "wire-shape node count is retained evidence"
    );
    let error = parse_ir_with_budget(
        &canonical,
        IrParseBudget {
            max_nodes: exact_nodes - 1,
            ..defaults
        },
    )
    .expect_err("node budget minus one must refuse");
    assert!(error.to_string().contains("node budget"));
    let error = parse_ir_with_budget(
        &canonical,
        IrParseBudget {
            max_depth: fs_scenario::ir::MAX_IR_PARSE_DEPTH + 1,
            ..defaults
        },
    )
    .expect_err("callers cannot raise the recursive hard-safety ceiling");
    assert!(error.to_string().contains("hard safety limit"));

    // The public parser must never unwind through Cheb1's asserting
    // constructor. Exercise every constructor precondition using malformed
    // authority bytes derived from a valid retained fixture.
    let profile_ir = write_ir(&rich_scenario());
    let empty_coefficients = truncate_form(&profile_ir, "(coeffs");
    let reversed_domain = profile_ir.replacen(") 0 1 (coeffs", ") 1 0 (coeffs", 1);
    assert_ne!(reversed_domain, profile_ir, "domain mutation must apply");
    let nonfinite_coefficient = profile_ir.replacen("(coeffs ", "(coeffs NaN ", 1);
    assert_ne!(
        nonfinite_coefficient, profile_ir,
        "coefficient mutation must apply"
    );
    for (name, malformed) in [
        ("empty coefficients", empty_coefficients),
        ("reversed domain", reversed_domain),
        ("non-finite coefficient", nonfinite_coefficient),
    ] {
        let error = parse_ir(&malformed).expect_err("malformed profile must return Parse");
        assert!(
            matches!(&error, fs_scenario::ScenarioError::Parse { .. }),
            "{name} returned the wrong error family: {error}"
        );
    }
    verdict(
        "sc-007",
        "byte/depth/node/atom/list budgets enforce exact boundaries; malformed Chebyshev IR returns Parse without unwind",
    );
}

#[test]
fn sc_008_result_apis_and_validation_refuse_nonfinite_public_values() {
    let bad_signals = [
        TimeSignal::Constant(QtyAny::dimensionless(f64::NAN)),
        TimeSignal::Ramp {
            t_start: 0.0,
            t_end: 1.0,
            from: QtyAny::dimensionless(0.0),
            to: QtyAny::dimensionless(f64::INFINITY),
        },
        TimeSignal::Table {
            times: vec![0.0, f64::NAN],
            values: vec![0.0, 1.0],
            dims: Dims::NONE,
            interp: Interp::Linear,
        },
        TimeSignal::Table {
            times: vec![0.0, 1.0],
            values: vec![0.0, f64::NEG_INFINITY],
            dims: Dims::NONE,
            interp: Interp::Hold,
        },
    ];
    for signal in &bad_signals {
        assert!(
            signal.eval(0.5).is_err(),
            "direct signal evaluation must independently refuse {signal:?}"
        );
    }
    let overflowing_ramp = TimeSignal::Ramp {
        t_start: 0.0,
        t_end: 1.0,
        from: QtyAny::dimensionless(f64::MAX),
        to: QtyAny::dimensionless(-f64::MAX),
    };
    assert_eq!(
        overflowing_ramp
            .eval(-1.0)
            .expect("left clamp")
            .value
            .to_bits(),
        f64::MAX.to_bits(),
        "clamped endpoints must not evaluate the overflowing delta"
    );
    assert_eq!(
        overflowing_ramp
            .eval(2.0)
            .expect("right clamp")
            .value
            .to_bits(),
        (-f64::MAX).to_bits(),
        "clamped endpoints must be returned exactly"
    );
    assert_eq!(
        overflowing_ramp.eval(0.5).expect("stable midpoint").value,
        0.0,
        "interior convex interpolation must not overflow opposite-sign endpoints"
    );
    let huge_time_table = TimeSignal::Table {
        times: vec![-f64::MAX, f64::MAX],
        values: vec![f64::MAX, -f64::MAX],
        dims: Dims::NONE,
        interp: Interp::Linear,
    };
    assert_eq!(
        huge_time_table
            .eval(0.0)
            .expect("scaled affine ratio")
            .value,
        0.0,
        "table interpolation must remain finite when direct time/value deltas overflow"
    );
    let degenerate_ramp = TimeSignal::Ramp {
        t_start: 1.0,
        t_end: 1.0,
        from: QtyAny::dimensionless(0.0),
        to: QtyAny::dimensionless(1.0),
    };
    assert!(
        degenerate_ramp.eval(0.0).is_err(),
        "a linear ramp requires a strict nonempty time interval"
    );

    let fixed_bad_quaternion = Frame {
        id: FrameId(1),
        name: "bad-quaternion".to_string(),
        parent: FrameId(0),
        motion: FrameMotion::Fixed {
            orientation: Quat {
                w: 2.0,
                x: 0.0,
                y: 0.0,
                z: 0.0,
            },
            translation: Vec3::new(0.0, 0.0, 0.0),
        },
    };
    assert!(
        fs_scenario::FrameTree::local_pose(&fixed_bad_quaternion, 0.0).is_err(),
        "non-unit fixed quaternions must not reach Motor::from_parts"
    );
    let rotating_bad_axis = Frame {
        id: FrameId(2),
        name: "bad-axis".to_string(),
        parent: FrameId(0),
        motion: FrameMotion::Rotating {
            axis: [f64::NAN, 0.0, 1.0],
            center: Vec3::new(0.0, f64::INFINITY, 0.0),
            rate: QtyAny::new(f64::INFINITY, RATE),
        },
    };
    assert!(
        fs_scenario::FrameTree::local_pose(&rotating_bad_axis, 1.0).is_err(),
        "invalid rotating geometry/rate must fail before motor construction"
    );
    let derived_overflow = Frame {
        id: FrameId(3),
        name: "derived-overflow".to_string(),
        parent: FrameId(0),
        motion: FrameMotion::Rotating {
            axis: [
                0.577_350_269_189_625_8,
                0.577_350_269_189_625_8,
                0.577_350_269_189_625_8,
            ],
            center: Vec3::new(f64::MAX, -f64::MAX, 0.0),
            rate: QtyAny::new(std::f64::consts::PI, RATE),
        },
    };
    assert!(
        fs_scenario::FrameTree::local_pose(&derived_overflow, 1.0).is_err(),
        "finite center/rate inputs that overflow motor coefficients must be refused"
    );
    let mut overflowing_chain = fs_scenario::FrameTree::new();
    for (id, parent) in [(1, 0), (2, 1), (3, 2)] {
        overflowing_chain.add(Frame {
            id: FrameId(id),
            name: format!("huge-translation-{id}"),
            parent: FrameId(parent),
            motion: FrameMotion::Fixed {
                orientation: Quat::identity(),
                translation: Vec3::new(f64::MAX, 0.0, 0.0),
            },
        });
    }
    assert!(
        overflowing_chain.world_pose(FrameId(3), 0.0).is_err(),
        "finite local motors whose parent composition overflows must be refused"
    );
    assert!(
        fs_scenario::FrameTree::new()
            .world_pose(FrameId(0), f64::NAN)
            .is_err(),
        "even the world-frame fast path must reject non-finite time"
    );
    let duplicate_frame = Frame {
        id: FrameId(9),
        name: "duplicate-a".to_string(),
        parent: FrameId(0),
        motion: FrameMotion::Fixed {
            orientation: Quat::identity(),
            translation: Vec3::new(0.0, 0.0, 0.0),
        },
    };
    let mut ambiguous_tree = fs_scenario::FrameTree::new();
    ambiguous_tree.add(duplicate_frame.clone());
    ambiguous_tree.add(Frame {
        name: "duplicate-b".to_string(),
        ..duplicate_frame
    });
    assert!(
        ambiguous_tree.world_pose(FrameId(9), 0.0).is_err(),
        "the direct pose API must not resolve duplicate ids by first match"
    );

    let mut scenario = Scenario::new("nonfinite", 9, Environment::earth_lab());
    scenario.environment.gravity[0].value = f64::INFINITY;
    scenario.environment.ambient_temperature.value = -1.0;
    scenario.environment.ambient_pressure.value = f64::NAN;
    scenario.frames.add(fixed_bad_quaternion);
    scenario.frames.add(rotating_bad_axis);
    scenario.base_bcs.push(BoundaryCondition {
        region: "hot".to_string(),
        physics: Physics::Thermal,
        kind: BcKind::Dirichlet,
        value: Some(BcValue::Uniform(QtyAny::new(
            f64::INFINITY,
            Dims([0, 0, 0, 1, 0, 0]),
        ))),
        compatibility: Some(Compat::Incompressible),
        frame: 0,
    });
    scenario.base_bcs.push(BoundaryCondition {
        region: "trace".to_string(),
        physics: Physics::Thermal,
        kind: BcKind::Dirichlet,
        value: Some(BcValue::Signal(TimeSignal::Table {
            times: vec![0.0, 1.0],
            values: vec![293.0, f64::NAN],
            dims: Dims([0, 0, 0, 1, 0, 0]),
            interp: Interp::Linear,
        })),
        compatibility: None,
        frame: 0,
    });
    scenario.contacts.push(ContactLaw {
        region_a: "a".to_string(),
        region_b: "b".to_string(),
        model: ContactModel::Coulomb {
            mu_s: 0.5,
            mu_k: f64::NAN,
        },
    });
    let violations = scenario.validate();
    for code in [
        "env-gravity-nonfinite",
        "env-temperature-range",
        "env-pressure-range",
        "frame-orientation-invalid",
        "frame-axis-not-unit",
        "frame-center-nonfinite",
        "frame-rate-nonfinite",
        "bc-value-nonfinite",
        "bc-compat-forbidden",
        "signal-value-nonfinite",
        "contact-coulomb-range",
    ] {
        assert!(
            violations.iter().any(|violation| violation.code == code),
            "missing fail-closed validation code {code}; got {violations:#?}"
        );
    }

    let mut overflowing_flux = Scenario::new("overflowing-flux", 1, Environment::earth_lab());
    for region in ["a", "b"] {
        overflowing_flux.base_bcs.push(BoundaryCondition {
            region: region.to_string(),
            physics: Physics::IncompressibleFlow,
            kind: BcKind::MassFlowInlet,
            value: Some(BcValue::Uniform(QtyAny::new(f64::MAX, MASS_FLOW))),
            compatibility: Some(Compat::Incompressible),
            frame: 0,
        });
    }
    assert!(
        overflowing_flux
            .validate()
            .iter()
            .any(|violation| violation.code == "flux-aggregation-nonfinite"),
        "finite mass-flow values whose sum overflows must not bypass compatibility admission"
    );

    let profile_flow = BoundaryCondition {
        region: "profile-inlet".to_string(),
        physics: Physics::IncompressibleFlow,
        kind: BcKind::MassFlowInlet,
        value: Some(BcValue::Profile(ChebProfile {
            cheb: fs_cheb::Cheb1::from_coeffs(0.0, 1.0, vec![2.0]),
            dims: MASS_FLOW,
        })),
        compatibility: Some(Compat::Incompressible),
        frame: 0,
    };
    let mut profile_violations = Vec::new();
    profile_flow.check(&mut profile_violations);
    assert!(
        profile_violations
            .iter()
            .any(|violation| violation.code == "bc-mass-flow-profile"),
        "a profile is not a certified total kg/s contribution: {profile_violations:#?}"
    );
    assert!(
        profile_flow.mass_flow_at(0.0).is_err(),
        "the direct compatibility API must not silently reinterpret a profile as zero flow"
    );
    let missing_flow = BoundaryCondition {
        region: "missing-total".to_string(),
        physics: Physics::IncompressibleFlow,
        kind: BcKind::MassFlowInlet,
        value: None,
        compatibility: Some(Compat::Incompressible),
        frame: 0,
    };
    assert!(
        missing_flow.mass_flow_at(0.0).is_err(),
        "a direct total-flow query must not conflate a missing value with an irrelevant BC"
    );
    let wrong_physics_flow = BoundaryCondition {
        physics: Physics::Thermal,
        value: Some(BcValue::Uniform(QtyAny::new(1.0, MASS_FLOW))),
        ..missing_flow
    };
    assert!(
        wrong_physics_flow.mass_flow_at(0.0).is_err(),
        "a MassFlowInlet attached to unsupported physics must not disappear from a direct query"
    );

    let unevaluable_flow = BoundaryCondition {
        region: "overflowing-signal".to_string(),
        physics: Physics::IncompressibleFlow,
        kind: BcKind::MassFlowInlet,
        value: Some(BcValue::Signal(TimeSignal::Chebfun(ChebProfile {
            cheb: fs_cheb::Cheb1::from_coeffs(-1.0, 1.0, vec![f64::MAX, 0.0, -f64::MAX]),
            dims: MASS_FLOW,
        }))),
        compatibility: Some(Compat::Incompressible),
        frame: 0,
    };
    assert!(
        unevaluable_flow.mass_flow_at(0.0).is_err(),
        "finite signal coefficients can still overflow during evaluation"
    );
    let mut unevaluable_flux = Scenario::new("unevaluable-flux", 1, Environment::earth_lab());
    unevaluable_flux.base_bcs.push(unevaluable_flow.clone());
    let unevaluable_violations = unevaluable_flux.validate();
    assert!(
        unevaluable_violations
            .iter()
            .any(|violation| violation.code == "flux-evaluation"),
        "signal evaluation failure must block compatibility admission: {unevaluable_violations:#?}"
    );
    unevaluable_flux.base_bcs.push(BoundaryCondition {
        region: "pressure-relief".to_string(),
        physics: Physics::IncompressibleFlow,
        kind: BcKind::PressureOutlet,
        value: Some(BcValue::Uniform(QtyAny::new(101_325.0, PRESSURE))),
        compatibility: None,
        frame: 0,
    });
    let outlet_violations = unevaluable_flux.validate();
    assert!(
        outlet_violations
            .iter()
            .any(|violation| violation.code == "flux-evaluation"),
        "a pressure outlet may absorb finite imbalance but must not hide an unevaluable inlet: {outlet_violations:#?}"
    );
    verdict(
        "sc-008",
        "direct signal/frame/flux APIs and whole-scenario validation reject non-finite/domain-invalid public values",
    );
}

#[test]
fn sc_008a_realization_budgets_and_carreau_domains_fail_closed() {
    let spectral = fixture_ensembles()[0].clone();
    let requested_samples = 512usize;
    let requested_work =
        requested_samples * (requested_samples / 2) + requested_samples + requested_samples / 2;
    let error = spectral
        .realize_with_budget(
            0,
            RealizationBudget {
                max_samples: requested_samples - 1,
                max_work: usize::MAX,
            },
        )
        .expect_err("sample cap minus one must refuse before allocation");
    assert!(error.to_string().contains("samples exceeds budget"));
    let error = spectral
        .realize_with_budget(
            0,
            RealizationBudget {
                max_samples: requested_samples,
                max_work: requested_work - 1,
            },
        )
        .expect_err("work cap minus one must refuse before allocation");
    assert!(error.to_string().contains("work"));
    let admitted = spectral
        .realize_with_budget(
            0,
            RealizationBudget {
                max_samples: requested_samples,
                max_work: requested_work,
            },
        )
        .expect("exact sample/work budgets admit the realization");
    assert_eq!(admitted.times.len(), requested_samples);
    assert!(admitted.values.iter().all(|value| value.is_finite()));

    let mut nonfinite_duration = spectral.clone();
    nonfinite_duration.duration.value = f64::NAN;
    assert!(nonfinite_duration.realize(0).is_err());
    let mut one_sample_spectrum = spectral.clone();
    one_sample_spectrum.duration.value = one_sample_spectrum.dt.value;
    let one_sample_error = one_sample_spectrum
        .realize(0)
        .expect_err("a one-sample spectral grid has no random harmonic and must refuse");
    assert!(one_sample_error.to_string().contains("spectral"));
    let mut one_sample_violations = Vec::new();
    one_sample_spectrum.check(&mut one_sample_violations);
    assert!(
        one_sample_violations
            .iter()
            .any(|violation| violation.code == "ensemble-spectral-grid"),
        "whole-scenario admission must catch the same degenerate grid"
    );
    let mut unrepresentable_ratio = spectral.clone();
    unrepresentable_ratio.duration.value = f64::MAX;
    unrepresentable_ratio.dt.value = f64::MIN_POSITIVE;
    assert!(
        unrepresentable_ratio.realize(0).is_err(),
        "infinite duration/dt ratio must refuse without a capacity attempt"
    );
    let mut unrepresentable_violations = Vec::new();
    unrepresentable_ratio.check(&mut unrepresentable_violations);
    assert!(
        unrepresentable_violations
            .iter()
            .any(|violation| violation.code == "ensemble-spectral-grid"),
        "structural admission must reject the same unrepresentable spectral grid"
    );
    let mut invalid_model = spectral.clone();
    let SpectrumModel::KanaiTajimi { omega_g, .. } = &mut invalid_model.model else {
        unreachable!();
    };
    omega_g.value = 0.0;
    assert!(
        invalid_model.realize(0).is_err(),
        "direct realize must not rely on Scenario::validate"
    );
    assert!(
        invalid_model.model.try_psd(1.0).is_err(),
        "the direct PSD API must independently refuse invalid public model values"
    );
    assert!(
        spectral.model.try_psd(f64::NAN).is_err(),
        "the direct PSD API must refuse non-finite frequencies"
    );
    assert!(
        spectral.model.try_psd(-1.0).is_err(),
        "a one-sided PSD must refuse negative angular frequencies"
    );

    let valid_carreau = fixture_ensembles()[2].clone();
    assert!(
        valid_carreau.model.try_psd(1.0).is_err(),
        "a parameter-band model must not masquerade as a zero spectral density"
    );
    for member in 0..valid_carreau.members {
        let realization = valid_carreau.realize(member).expect("valid Carreau member");
        assert!(realization.values.iter().all(|value| value.is_finite()));
        assert!(
            realization.values[0] >= realization.values[1],
            "every independently sampled member must preserve eta_zero >= eta_inf"
        );
    }
    let mut nonfinite_band_grid = valid_carreau.clone();
    nonfinite_band_grid.duration.value = f64::NAN;
    let mut grid_violations = Vec::new();
    nonfinite_band_grid.check(&mut grid_violations);
    assert!(
        grid_violations
            .iter()
            .any(|violation| violation.code == "ensemble-time-range"),
        "ignored band-grid placeholders still travel in canonical IR and must be finite"
    );
    assert!(nonfinite_band_grid.realize(0).is_err());
    assert!(
        parse_ir(&write_ir(&Scenario {
            ensembles: vec![nonfinite_band_grid],
            ..Scenario::new("bad-band-grid", 1, Environment::earth_lab())
        }))
        .is_err(),
        "the retained malformed fixture demonstrates why admission must reject non-finite band-grid fields"
    );

    let mut nonpositive = valid_carreau.clone();
    let SpectrumModel::CarreauBand { eta_zero, .. } = &mut nonpositive.model else {
        unreachable!();
    };
    eta_zero[0].value = 0.0;
    let mut invalid_n = valid_carreau.clone();
    let SpectrumModel::CarreauBand { n, .. } = &mut invalid_n.model else {
        unreachable!();
    };
    *n = [0.0, 1.2];
    let mut crossing_viscosities = valid_carreau.clone();
    let SpectrumModel::CarreauBand {
        eta_zero, eta_inf, ..
    } = &mut crossing_viscosities.model
    else {
        unreachable!();
    };
    eta_zero[0].value = 0.1;
    eta_zero[1].value = 0.2;
    eta_inf[0].value = 0.15;
    eta_inf[1].value = 0.25;

    for (ensemble, expected_code) in [
        (nonpositive, "ensemble-carreau-positive-finite"),
        (invalid_n, "ensemble-carreau-power-index"),
        (crossing_viscosities, "ensemble-carreau-viscosity-order"),
    ] {
        let mut violations = Vec::new();
        ensemble.check(&mut violations);
        assert!(
            violations
                .iter()
                .any(|violation| violation.code == expected_code),
            "missing {expected_code}: {violations:#?}"
        );
        assert!(
            ensemble.realize(0).is_err(),
            "direct realization must reject {expected_code}"
        );
    }
    verdict(
        "sc-008a",
        "sample/work exact boundaries, ratio overflow, direct model validation, and physical Carreau domains fail closed",
    );
}

#[test]
fn sc_009_identity_and_unordered_contact_integrity_fail_closed() {
    let fixed = || FrameMotion::Fixed {
        orientation: Quat::identity(),
        translation: Vec3::new(0.0, 0.0, 0.0),
    };
    let mut scenario = Scenario::new("", 17, Environment::earth_lab());
    scenario.frames.add(Frame {
        id: FrameId(1),
        name: String::new(),
        parent: FrameId(0),
        motion: fixed(),
    });
    scenario.frames.add(Frame {
        id: FrameId(2),
        name: "fixture-frame".to_string(),
        parent: FrameId(0),
        motion: fixed(),
    });
    scenario.frames.add(Frame {
        id: FrameId(3),
        name: "fixture-frame".to_string(),
        parent: FrameId(0),
        motion: fixed(),
    });
    scenario.base_bcs.push(BoundaryCondition {
        region: String::new(),
        physics: Physics::IncompressibleFlow,
        kind: BcKind::WallNoSlip,
        value: None,
        compatibility: None,
        frame: 0,
    });
    scenario.cases.extend([
        LoadCase {
            name: String::new(),
            bcs: Vec::new(),
        },
        LoadCase {
            name: "load".to_string(),
            bcs: Vec::new(),
        },
        LoadCase {
            name: "load".to_string(),
            bcs: Vec::new(),
        },
    ]);
    scenario.combinations.extend([
        Combination {
            name: String::new(),
            terms: vec![
                ("load".to_string(), 1.0),
                ("load".to_string(), 2.0),
                (String::new(), 1.0),
            ],
        },
        Combination {
            name: "combo".to_string(),
            terms: Vec::new(),
        },
        Combination {
            name: "combo".to_string(),
            terms: Vec::new(),
        },
    ]);
    let mut ensembles = fixture_ensembles();
    ensembles[0].name.clear();
    ensembles[1].name = "ensemble".to_string();
    ensembles[2].name = "ensemble".to_string();
    scenario.ensembles = ensembles;
    scenario.contacts.extend([
        ContactLaw {
            region_a: "a".to_string(),
            region_b: "b".to_string(),
            model: ContactModel::Frictionless,
        },
        ContactLaw {
            region_a: "b".to_string(),
            region_b: "a".to_string(),
            model: ContactModel::Tied,
        },
        ContactLaw {
            region_a: "a".to_string(),
            region_b: "b".to_string(),
            model: ContactModel::Frictionless,
        },
        ContactLaw {
            region_a: String::new(),
            region_b: "c".to_string(),
            model: ContactModel::Frictionless,
        },
        ContactLaw {
            region_a: "self".to_string(),
            region_b: "self".to_string(),
            model: ContactModel::Tied,
        },
    ]);

    let violations = scenario.validate();
    for code in [
        "scenario-name-empty",
        "frame-name-empty",
        "frame-name-duplicate",
        "bc-region-empty",
        "case-name-empty",
        "case-name-duplicate",
        "combo-name-empty",
        "combo-name-duplicate",
        "combo-case-empty",
        "combo-term-duplicate",
        "ensemble-name-empty",
        "ensemble-name-duplicate",
        "contact-region-empty",
        "contact-self-pair",
        "contact-pair-conflict",
        "contact-pair-duplicate",
    ] {
        assert!(
            violations.iter().any(|violation| violation.code == code),
            "missing {code}: {violations:#?}"
        );
    }
    let repeated_term = violations
        .iter()
        .find(|violation| violation.code == "combo-term-duplicate")
        .expect("repeated term diagnosis");
    assert!(repeated_term.what.contains("terms 0 and 1"));
    let conflicting_pair = violations
        .iter()
        .find(|violation| violation.code == "contact-pair-conflict")
        .expect("unordered conflict diagnosis");
    assert!(conflicting_pair.what.contains("row 0") && conflicting_pair.what.contains("row 1"));
    verdict(
        "sc-009",
        "nonempty unique exact identities, duplicate combination terms, and unordered contact-pair semantics fail closed with declaration provenance",
    );
}
