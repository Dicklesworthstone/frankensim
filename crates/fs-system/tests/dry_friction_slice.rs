//! Slice C (`frankensim-rc-1109-root-3gc12.6.3`): dry-friction heating with
//! stick, slip, arrest, reversal, and heat partitioning on the non-smooth
//! `LoadLine` tier via `LoweringRegistry::with_kernel_defaults()`.

use std::f64::consts::PI;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::{ContentHash, hash_bytes};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb as m;
use fs_ports::component::conformance::{Case, check_component, check_load_line};
use fs_ports::component::{Component, LoadLine};
use fs_qty::{Dims, QtyAny};
use fs_scenario::assembly::draft as d;
use fs_scenario::assembly::{MaterialRef, ParameterKey, ReferenceState};
use fs_scenario::{EntityCatalog, EntityDeclaration, EntityId, FrameId, InterfacePair, Scenario};
use fs_system::elements::{Mass, Motion, Slider, SliderEnd, Slip, SourceAccount};
use fs_system::lowering::InterfaceKindTag;
use fs_system::{CardStore, LoweredSystem, LoweringError, LoweringRegistry, System, SystemBuilder};

const L: Dims = Dims([1, 0, 0, 0, 0, 0]);
const T_DIM: Dims = Dims([0, 0, 1, 0, 0, 0]);
const V_DIM: Dims = Dims([1, 0, -1, 0, 0, 0]);
const F_DIM: Dims = Dims([1, 1, -2, 0, 0, 0]);
const TEMP: Dims = Dims([0, 0, 0, 1, 0, 0]);

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let k = StreamKey {
        seed: 42,
        kernel_id: 1,
        tile: 0,
        iteration: 0,
    };
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|a| {
        f(&Cx::new(
            &gate,
            a,
            k,
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

fn kv(name: &str, val: f64, dims: Dims) -> (ParameterKey, QtyAny) {
    (ParameterKey::new(name, dims), QtyAny::new(val, dims))
}

fn claim(name: &str, value: f64, dims: Dims, validity: ValidityDomain) -> m::PropertyClaim {
    m::PropertyClaim {
        key: m::PropertyKey::new(name, dims),
        value: m::PropertyValue::Scalar { value, dims },
        validity,
        uncertainty: m::UncertaintyModel::Unstated,
        interpolation: m::InterpolationPolicy::ConstantWithinValidity,
        observations: Vec::new(),
        provenance: m::Provenance {
            source: name.to_owned(),
            license: "CC0-1.0".to_owned(),
            artifact: Some(hash_bytes(name.as_bytes())),
        },
    }
}

fn mat_id(c: &str) -> m::MaterialStateId {
    m::MaterialStateId {
        chemistry: c.into(),
        phase: "s".into(),
        process: "dry".into(),
        revision: 0,
    }
}

fn unit_mat(rho: f64, cp: f64) -> m::MaterialCard {
    let mut cs = m::ClaimSet::new();
    let v = ValidityDomain::unconstrained().with("T", 200.0, 800.0);
    for (k, x, dim) in [
        ("density", rho, Dims([-3, 1, 0, 0, 0, 0])),
        ("specific_heat_capacity", cp, Dims([2, 0, -2, -1, 0, 0])),
    ] {
        cs.insert_claim(claim(k, x, dim, v.clone())).unwrap();
    }
    m::MaterialCard::assemble(mat_id("m"), cs, Vec::new()).unwrap()
}

fn iface_card(cs: m::ClaimSet) -> m::InterfaceSystemCard {
    let s = m::SurfaceSpec {
        material: mat_id("52100"),
        texture_frame: "g".into(),
    };
    let ctx = m::SystemContext {
        medium: "dry".into(),
        third_body: None,
        environment: "air".into(),
        history: "1h".into(),
    };
    m::InterfaceSystemCard::assemble(s.clone(), s, ctx, cs, Vec::new()).unwrap()
}

fn friction_card(mu_s: f64, mu_k: f64, max_n: f64, max_v: f64) -> m::InterfaceSystemCard {
    let mut cs = m::ClaimSet::new();
    let val = ValidityDomain::unconstrained()
        .with("temperature", 250.0, 500.0)
        .with("normal_force", 1.0, max_n)
        .with("sliding_speed", 0.0, max_v);
    for (k, v) in [
        ("static-friction-coefficient", mu_s),
        ("kinetic-friction-coefficient", mu_k),
    ] {
        cs.insert_claim(claim(k, v, Dims::NONE, val.clone()))
            .unwrap();
    }
    iface_card(cs)
}

fn init_st(vx: f64, tk: Option<f64>) -> d::DraftInitialState {
    let q = QtyAny::new;
    d::DraftInitialState {
        position: [q(0.0, L); 3],
        orientation: [1.0, 0.0, 0.0, 0.0],
        linear_velocity: [q(vx, V_DIM), q(0.0, V_DIM), q(0.0, V_DIM)],
        angular_velocity: [q(0.0, Dims([0, 0, -1, 0, 0, 0])); 3],
        temperature: tk.map(|t| q(t, TEMP)),
    }
}

type Side = (
    &'static str,
    d::DraftPartKind,
    ContentHash,
    bool,
    Option<d::DraftInitialState>,
);

fn make_part(s: &Side, entity: EntityId) -> d::DraftPartDecl {
    use d::PhysicsDomain as P;
    let dim = QtyAny::new(1.0, L);
    let eb = d::FidelityErrorBudget {
        relative_tol: 1e-5,
        absolute_tol: QtyAny::dimensionless(1e-12),
        max_work_units: Some(1_000_000),
    };
    let mr = MaterialRef {
        card: s.2,
        query_point: vec![kv("T", 300.0, TEMP)],
        frame: FrameId(0),
        layer: None,
        unknown: Vec::new(),
    };
    d::DraftPartDecl {
        name: s.0.into(),
        entity,
        kind: s.1,
        geometry: d::DraftGeometrySource::Primitive(d::PrimitiveShape::Box {
            lx: dim,
            ly: dim,
            lz: dim,
        }),
        materials: vec![mr],
        physics: if s.3 {
            [P::Mechanics, P::Thermal].into()
        } else {
            [P::Mechanics].into()
        },
        fidelity: d::FidelityPolicy {
            rungs: vec![d::FidelityRung::Lumped],
            error_budget: eb,
        },
        uncertainty: vec![],
        initial: s.4.clone(),
    }
}

type Built = (
    Vec<d::JointDecl>,
    Vec<d::DraftInterfaceDecl>,
    Vec<d::SourceDecl>,
);

fn admit_and_lower(
    (root, catalog): (EntityId, EntityCatalog),
    parts: Vec<d::DraftPartDecl>,
    (joints, interfaces, sources): Built,
    (dt, t_end, t_amb): (f64, f64, f64),
    cards: &CardStore,
    reg: &LoweringRegistry,
) -> Result<LoweredSystem, LoweringError> {
    let mut env = d::DraftEnvironment::earth_lab_air();
    env.base.ambient_temperature = QtyAny::new(t_amb, TEMP);
    let scen = Scenario::new("p", 42, env.base.clone());
    let inp = d::DraftAssemblyInput {
        name: "p".into(),
        root,
        catalog,
        reference_state: Some(ReferenceState::Assembled),
        environment: env,
        time_plan: d::TimePlan {
            horizon: QtyAny::new(t_end, T_DIM),
            step_policy: d::StepPolicy::Fixed {
                dt: QtyAny::new(dt, T_DIM),
            },
            default_sample_rate: QtyAny::new(1.0 / dt, Dims([0, 0, -1, 0, 0, 0])),
            scheme: d::TimeScheme::Explicit,
        },
        parts,
        joints,
        gear_meshes: vec![],
        interfaces,
        sources,
        controllers: vec![],
        probes: vec![],
    };
    reg.lower(&d::DraftAssembly::admit(inp, &scen).unwrap(), cards)
}

fn lower_draft_ifaces(
    (side_a, side_b): (&Side, &Side),
    (dt, t_end): (f64, f64),
    cards: &CardStore,
    build: impl FnOnce(&mut EntityCatalog, EntityId, [EntityId; 2]) -> Built,
) -> Result<LoweredSystem, LoweringError> {
    let mut cat = EntityCatalog::new();
    let root = cat.declare(EntityDeclaration::assembly("p")).unwrap();
    let pa = cat
        .declare(EntityDeclaration::part(root, side_a.0))
        .unwrap();
    let pb = cat
        .declare(EntityDeclaration::part(root, side_b.0))
        .unwrap();
    let items = build(&mut cat, root, [pa, pb]);
    let parts = vec![make_part(side_a, pa), make_part(side_b, pb)];
    let reg = LoweringRegistry::with_kernel_defaults();
    admit_and_lower((root, cat), parts, items, (dt, t_end, 300.0), cards, &reg)
}

fn lower_pair(
    sides: (&Side, &Side),
    (card, third_body, state): (
        ContentHash,
        Option<d::ThirdBodyDecl>,
        Vec<(ParameterKey, QtyAny)>,
    ),
    (spring_k, belt): (Option<f64>, Option<d::TimeHistory>),
    plan: (f64, f64),
    cards: &CardStore,
) -> Result<LoweredSystem, LoweringError> {
    lower_draft_ifaces(sides, plan, cards, |cat, root, [pa, pb]| {
        let sa = cat.declare(EntityDeclaration::surface(pa, "sa")).unwrap();
        let sb = cat.declare(EntityDeclaration::surface(pb, "sb")).unwrap();
        let cid = cat
            .declare(EntityDeclaration::interface(
                root,
                "pad",
                InterfacePair::ordered(sa, sb),
            ))
            .unwrap();
        let joints = spring_k
            .map(|k| d::JointDecl {
                name: "spring".into(),
                parent: root,
                child: pa,
                frame: FrameId(0),
                kind: d::JointKind::Prismatic,
                axis: [1.0, 0.0, 0.0],
                pitch: None,
                limits: None,
                damping: None,
                stiffness: Some(QtyAny::new(k, Dims([0, 1, -2, 0, 0, 0]))),
            })
            .into_iter()
            .collect();
        let sources = belt
            .map(|history| d::SourceDecl {
                name: "belt_v".into(),
                target: d::SourceTarget::Port {
                    entity: pb,
                    port: "mechanics".into(),
                },
                history,
            })
            .into_iter()
            .collect();
        let kind = d::DraftInterfaceKind::MechanicalContact {
            law_card: card,
            third_body,
            state,
        };
        let iface = d::DraftInterfaceDecl {
            name: "pad".into(),
            entity: cid,
            feature_a: sa,
            feature_b: sb,
            kind,
        };
        (joints, vec![iface], sources)
    })
}

fn run_stick_slip(
    cards: &CardStore,
    mat: ContentHash,
    fcard: ContentHash,
    (dt, steps, check_rcpt): (f64, u32, bool),
) -> [f64; 5] {
    let (v_drive, stiff, n_f, share_a, x_k) = (0.1_f64, 100.0_f64, 10.0_f64, 0.6_f64, 0.03_f64);
    let blk = (
        "block",
        d::DraftPartKind::RigidBody,
        mat,
        true,
        Some(init_st(v_drive, Some(300.0))),
    );
    let drv = ("belt", d::DraftPartKind::DrivenBody, mat, false, None);
    let st = vec![
        kv("normal_force", n_f, F_DIM),
        kv("heat_partition_a", share_a, Dims::NONE),
    ];
    let drive = Some(d::TimeHistory::Constant(QtyAny::new(v_drive, V_DIM)));
    let h = f64::from(steps) * dt;
    let mut lo = lower_pair(
        (&blk, &drv),
        (fcard, None, st),
        (Some(stiff), drive),
        (dt, h),
        cards,
    )
    .unwrap();
    if check_rcpt {
        let r = lo.receipt.canonical_json();
        assert!(
            r.contains("static-friction-coefficient.held")
                && r.contains("kinetic-friction-coefficient.held")
                && r.contains("held: mu does not follow the interface's live state (one-way)")
                && r.contains("heat_partition_a")
                && r.contains("heat_partition_ambient")
        );
    }
    let (sl, bt) = (lo.node_map["interface.pad"], lo.node_map["block.thermal"]);
    let (bm, blt) = (lo.node_map["block.mechanics"], lo.node_map["belt"]);
    let (mut crosses, mut caught, mut was_stuck, mut x_prev) = (Vec::new(), Vec::new(), true, 0.0);
    let (mut e0, mut e1) = (0.0, 0.0);
    with_cx(|cx| {
        for step in 1..=steps {
            lo.system.step(dt, cx).unwrap();
            let pos = lo.system.state::<Motion>(bm).unwrap().position_m;
            let now = lo.system.state::<Slip>(sl).unwrap().is_stuck();
            if now && x_prev < x_k && pos >= x_k {
                crosses.push((f64::from(step - 1) + (x_k - x_prev) / (pos - x_prev)) * dt);
            }
            if now && !was_stuck {
                let ein = lo.system.state::<SourceAccount>(blt).unwrap().energy_in_j;
                if caught.is_empty() {
                    e0 = ein;
                }
                e1 = ein;
                caught.push(f64::from(step) * dt);
            }
            was_stuck = now;
            x_prev = pos;
        }
    });
    assert!(crosses.len() >= 3 && caught.len() >= 2);
    let n_cr = f64::from(u32::try_from(crosses.len() - 2).unwrap());
    let cyc = f64::from(u32::try_from(caught.len() - 1).unwrap());
    let t_interp = (crosses[crosses.len() - 1] - crosses[1]) / n_cr;
    let t_grid = (caught[caught.len() - 1] - caught[0]) / cyc;
    let q_tot = 100.0 * (*lo.system.state::<f64>(bt).unwrap() - 300.0) / share_a;
    let mot = *lo.system.state::<Motion>(bm).unwrap();
    let w_diss = -lo.system.state::<SourceAccount>(blt).unwrap().energy_in_j
        - (0.5 * (mot.velocity_m_s.powi(2) - v_drive * v_drive)
            + 0.5 * stiff * mot.position_m.powi(2));
    let rel_wh = (q_tot - w_diss).abs() / w_diss;
    let led = lo.system.ledger().energy_residual_j().abs() / w_diss;
    [t_interp, t_grid, -(e1 - e0) / cyc, rel_wh, led]
}

#[test]
fn block_on_belt_stick_slip_period_work_and_heat_partition() {
    let mut cards = CardStore::new();
    let mat = cards.insert_material(unit_mat(1.0, 100.0));
    let fcard = cards
        .insert_interface("dry-friction", friction_card(0.5, 0.3, 50.0, 5.0))
        .unwrap();
    let delta: f64 = (0.5 - 0.3) * 10.0 / 100.0;
    let exact_t = 2.0 * delta / 0.1 + 2.0 * (PI - (delta * 10.0 / 0.1).atan()) / 10.0;
    let exact_w = 0.3 * 10.0 * 0.1 * exact_t;
    for (idx, (dt, steps)) in [
        (2.0e-3_f64, 3_000_u32),
        (1.0e-3, 6_000),
        (5.0e-4, 12_000),
        (2.5e-4, 24_000),
        (1.25e-4, 48_000),
        (1.5625e-5, 166_400),
    ]
    .into_iter()
    .enumerate()
    {
        let [t_i, t_g, w_c, rel_wh, rel_l] =
            run_stick_slip(&cards, mat, fcard, (dt, steps, idx == 0));
        let (err_i, err_g) = ((t_i - exact_t).abs(), (t_g - exact_t).abs());
        let w_rel = (w_c - exact_w).abs() / exact_w;
        assert!(
            err_i < dt && err_g < dt,
            "dt={dt} T_i={t_i} T_g={t_g} exact={exact_t}"
        );
        assert!(w_rel <= dt / exact_t + 1e-9, "dt={dt} w_rel={w_rel}");
        assert!(rel_wh <= 1e-9 && rel_l <= 1e-12);
        eprintln!(
            "dt={dt:.7} T_i={t_i:.9} err_i={err_i:.3e} T_g={t_g:.9} err_g={err_g:.3e} W_rel={w_rel:.3e} rel_wh={rel_wh:.3e}"
        );
    }
}

#[test]
fn nasa_52100_kinetic_only_pack_refuses() {
    const NASA: &str = concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../data/matdb/seed-v1/nasa-52100-dry-air-interface/interface.tsv"
    );
    let (mut prop, mut num, mut val) = (String::new(), 0.0, ValidityDomain::unconstrained());
    for line in std::fs::read_to_string(NASA).unwrap().lines() {
        let c: Vec<&str> = line.split('\t').collect();
        if c.first() == Some(&"scalar") {
            c[3].clone_into(&mut prop);
            num = c[4].parse().unwrap();
        } else if c.first() == Some(&"validity") {
            val = val.with(c[2], c[3].parse().unwrap(), c[4].parse().unwrap());
        }
    }
    let mut cs = m::ClaimSet::new();
    cs.insert_claim(claim(&prop, num, Dims::NONE, val)).unwrap();
    let mut cards = CardStore::new();
    let mat = cards.insert_material(unit_mat(1.0, 100.0));
    let nasa = cards
        .insert_interface("dry-friction", iface_card(cs))
        .unwrap();
    let blk = ("a", d::DraftPartKind::RigidBody, mat, true, None);
    let drv = ("b", d::DraftPartKind::DrivenBody, mat, false, None);
    let st = vec![kv("normal_force", 10.0, F_DIM)];
    let msg = lower_pair(
        (&blk, &drv),
        (nasa, None, st),
        (None, None),
        (0.01, 0.1),
        &cards,
    )
    .err()
    .unwrap()
    .to_string();
    assert!(msg.contains("static-friction-coefficient"), "{msg}");
}

#[test]
fn two_free_thermal_bodies_arrest_partition_and_harmonic_entropy() {
    let mut cards = CardStore::new();
    let (mat_a, mat_b) = (
        cards.insert_material(unit_mat(1.0, 100.0)),
        cards.insert_material(unit_mat(3.0, 100.0)),
    );
    let fcard = cards
        .insert_interface("dry-friction", friction_card(0.5, 0.3, 50.0, 5.0))
        .unwrap();
    let (share_a, ta0, tb0) = (0.25_f64, 300.0_f64, 320.0_f64);
    let body_a = (
        "a",
        d::DraftPartKind::RigidBody,
        mat_a,
        true,
        Some(init_st(2.0, Some(ta0))),
    );
    let body_b = (
        "b",
        d::DraftPartKind::RigidBody,
        mat_b,
        true,
        Some(init_st(0.0, Some(tb0))),
    );
    let st = vec![
        kv("normal_force", 10.0, F_DIM),
        kv("heat_partition_a", share_a, Dims::NONE),
    ];
    let mut lo = lower_pair(
        (&body_a, &body_b),
        (fcard, None, st),
        (None, None),
        (1e-3, 0.8),
        &cards,
    )
    .unwrap();
    with_cx(|cx| {
        for _ in 0..800 {
            lo.system.step(1e-3, cx).unwrap();
        }
    });
    let va = lo
        .system
        .state::<Motion>(lo.node_map["a.mechanics"])
        .unwrap()
        .velocity_m_s;
    let vb = lo
        .system
        .state::<Motion>(lo.node_map["b.mechanics"])
        .unwrap()
        .velocity_m_s;
    assert!((va - 0.5).abs() <= 1e-12 && (vb - 0.5).abs() <= 1e-12);
    assert!(
        lo.system
            .state::<Slip>(lo.node_map["interface.pad"])
            .unwrap()
            .is_stuck()
    );
    let ta1 = *lo.system.state::<f64>(lo.node_map["a.thermal"]).unwrap();
    let tb1 = *lo.system.state::<f64>(lo.node_map["b.thermal"]).unwrap();
    let (qa, qb, q_ex) = (100.0 * (ta1 - ta0), 300.0 * (tb1 - tb0), 1.5_f64);
    let rel_wh = ((qa + qb) - q_ex).abs() / q_ex;
    let ds_ex = 100.0 * (ta1 / ta0).ln() + 300.0 * (tb1 / tb0).ln();
    let ds_led = lo.system.ledger().entropy_generated_j_per_k();
    assert!(
        (qa - share_a * q_ex).abs() / q_ex <= 1e-9
            && (qb - (1.0 - share_a) * q_ex).abs() / q_ex <= 1e-9
    );
    assert!(rel_wh <= 1e-9 && (ds_led - ds_ex).abs() / ds_ex <= 1e-9);
    assert!(
        lo.receipt
            .canonical_json()
            .contains("declared: side_a thermal")
    );
    eprintln!("two-body arrest: rel_wh={rel_wh:.3e} dS={ds_led:.9} exact={ds_ex:.9}");
}

#[test]
fn reversing_piecewise_linear_belt_tracks_breaks_reverses_and_arrests() {
    let mut cards = CardStore::new();
    let mat = cards.insert_material(unit_mat(1.0, 100.0));
    let fcard = cards
        .insert_interface("dry-friction", friction_card(0.5, 0.3, 50.0, 5.0))
        .unwrap();
    let q = QtyAny::new;
    let hist = d::TimeHistory::PiecewiseLinear {
        knots: vec![
            (q(0.0, T_DIM), q(0.0, V_DIM)),
            (q(0.2, T_DIM), q(0.6, V_DIM)),
            (q(0.4, T_DIM), q(-0.6, V_DIM)),
            (q(0.8, T_DIM), q(-0.6, V_DIM)),
        ],
    };
    let blk = (
        "block",
        d::DraftPartKind::RigidBody,
        mat,
        true,
        Some(init_st(0.0, Some(300.0))),
    );
    let drv = ("belt", d::DraftPartKind::DrivenBody, mat, false, None);
    let st = vec![kv("normal_force", 10.0, F_DIM)];
    let mut lo = lower_pair(
        (&blk, &drv),
        (fcard, None, st),
        (None, Some(hist)),
        (1e-3, 0.8),
        &cards,
    )
    .unwrap();
    let (bm, bt) = (lo.node_map["block.mechanics"], lo.node_map["block.thermal"]);
    let (sl, blt) = (lo.node_map["interface.pad"], lo.node_map["belt"]);
    with_cx(|cx| {
        for step in 1..=800 {
            lo.system.step(1e-3, cx).unwrap();
            if step == 200 {
                let vel = lo.system.state::<Motion>(bm).unwrap().velocity_m_s;
                assert!(
                    (vel - 0.6).abs() <= 1e-12 && lo.system.state::<Slip>(sl).unwrap().is_stuck()
                );
                assert!((*lo.system.state::<f64>(bt).unwrap() - 300.0).abs() <= 1e-12);
            }
        }
    });
    let v_end = lo.system.state::<Motion>(bm).unwrap().velocity_m_s;
    let q_blk = 100.0 * (*lo.system.state::<f64>(bt).unwrap() - 300.0);
    let w_diss = -lo.system.state::<SourceAccount>(blt).unwrap().energy_in_j - 0.5 * v_end * v_end;
    let rel_wh = (q_blk - w_diss).abs() / w_diss;
    assert!((v_end + 0.6).abs() <= 1e-12 && lo.system.state::<Slip>(sl).unwrap().is_stuck());
    assert!(q_blk > 0.1 && rel_wh <= 1e-9);
    let rcpt = lo.receipt.canonical_json();
    assert!(rcpt.contains("default: the only thermal side takes all the heat"));
    assert!(!rcpt.contains("declared: single thermal"));
    eprintln!("reversing belt: q={q_blk:.6} J rel_wh={rel_wh:.3e}");
}

fn refusal_msg(
    cards: &CardStore,
    mat: ContentHash,
    card: ContentHash,
    third: Option<d::ThirdBodyDecl>,
    state: Vec<(ParameterKey, QtyAny)>,
    two_thermal: bool,
) -> String {
    let b_kind = if two_thermal {
        d::DraftPartKind::RigidBody
    } else {
        d::DraftPartKind::DrivenBody
    };
    let body_a = ("a", d::DraftPartKind::RigidBody, mat, true, None);
    let body_b = ("b", b_kind, mat, two_thermal, None);
    lower_pair(
        (&body_a, &body_b),
        (card, third, state),
        (None, None),
        (0.01, 0.1),
        cards,
    )
    .err()
    .unwrap()
    .to_string()
}

#[test]
fn dry_friction_refusals_cover_validity_order_third_body_partition_and_slip() {
    let mut cards = CardStore::new();
    let mat = cards.insert_material(unit_mat(1.0, 100.0));
    let good = cards
        .insert_interface("dry-friction", friction_card(0.5, 0.3, 50.0, 5.0))
        .unwrap();
    let bad_mu = cards
        .insert_interface("dry-friction", friction_card(0.2, 0.4, 50.0, 5.0))
        .unwrap();
    let n10 = kv("normal_force", 10.0, F_DIM);
    let tb = d::ThirdBodyDecl {
        material_card: mat,
        film_thickness: QtyAny::new(1e-4, L),
        query_point: vec![],
    };
    let rows = [
        (
            good,
            None,
            vec![kv("normal_force", 90.0, F_DIM)],
            false,
            "outside the validity",
        ),
        (
            good,
            None,
            vec![n10.clone(), kv("sliding_speed", 9.0, V_DIM)],
            false,
            "outside the validity",
        ),
        (
            good,
            None,
            vec![n10.clone(), kv("sliding_speed", -0.5, V_DIM)],
            false,
            "sliding_speed",
        ),
        (
            bad_mu,
            None,
            vec![n10.clone()],
            false,
            "static-friction-coefficient",
        ),
        (good, None, vec![], false, "normal_force"),
        (
            good,
            None,
            vec![kv("normal_force", -1.0, F_DIM)],
            false,
            "normal_force",
        ),
        (good, Some(tb), vec![n10.clone()], false, "third-body"),
        (good, None, vec![n10.clone()], true, "heat_partition_a"),
        (
            good,
            None,
            vec![n10.clone(), kv("heat_partition_a", 1.2, Dims::NONE)],
            true,
            "heat_partition_a",
        ),
        (
            good,
            None,
            vec![n10, kv("initial_slip", 1.0, V_DIM)],
            false,
            "initial_slip",
        ),
    ];
    for (card, third, state, two_t, want) in rows {
        let msg = refusal_msg(&cards, mat, card, third, state, two_t);
        assert!(msg.contains(want), "expected `{want}` in: {msg}");
    }
    let ma = ("a", d::DraftPartKind::RigidBody, mat, false, None);
    let mb = ("b", d::DraftPartKind::DrivenBody, mat, false, None);
    let amb = lower_pair(
        (&ma, &mb),
        (good, None, vec![kv("normal_force", 10.0, F_DIM)]),
        (None, None),
        (0.01, 0.1),
        &cards,
    )
    .unwrap();
    let arcpt = amb.receipt.canonical_json();
    assert!(arcpt.contains("default: no thermal side; all heat to ambient"));
    assert!(!arcpt.contains("declared: ambient sink"));
}

#[test]
fn slider_identity_conformance_and_liftoff() {
    let id = |s: &Slider| {
        let mut w = Vec::new();
        s.identity_words(&mut w);
        w
    };
    let fixed = Slider::new(10.0, 0.5, 0.3).unwrap();
    assert_eq!(id(&fixed), [5.0_f64.to_bits(), 3.0_f64.to_bits(), 0, 0]);
    let sig = Slider::loaded_by_signal(0.5, 0.3).unwrap();
    assert!(Slider::loaded_by_signal(0.2, 0.3).is_err());
    let alt = Slider::loaded_by_signal(0.5, 0.25).unwrap();
    assert!(id(&fixed) != id(&sig) && id(&sig) != id(&alt));
    let told = sig.clone().told(SliderEnd::A);
    let names: Vec<_> = told.ports().iter().map(|p| p.id().as_str()).collect();
    assert_eq!(names, ["a", "b", "heat", "normal_n", "span_a"]);
    with_cx(|cx| {
        let (old, stale, line) = (Slip::sliding(-1.0), Slip::stuck(), [0.05, 0.0, 0.0, 0.0]);
        for n in [10.0, 0.0] {
            let inp = [-1.0, 0.0, 300.0, n];
            let case = Case {
                state: &old,
                stale: &stale,
                dt_s: 1.0e-3,
                inputs: &inp,
            };
            check_component(&sig, &case, cx).unwrap();
            check_load_line(&sig, &case, &line, cx).unwrap();
        }
        let run = |s: &Slider, inp: &[f64]| {
            let (mut st, mut out) = (stale, [1.0; 3]);
            s.propose_on_line(&old, 1.0e-3, inp, &line, &mut st, &mut out, cx)
                .map(|r| (st, out, r))
        };
        let res_fix = run(&fixed, &[-1.0, 0.0, 300.0]).unwrap();
        let res_sig = run(&sig, &[-1.0, 0.0, 300.0, 10.0]).unwrap();
        assert_eq!(res_fix, res_sig);
        let (s0, o0, _) = run(&sig, &[-1.0, 0.0, 300.0, 0.0]).unwrap();
        assert!(o0.iter().all(|&v| v == 0.0) && !s0.is_stuck());
        for bad in [-1.0e-12, f64::NAN, f64::INFINITY] {
            assert!(run(&sig, &[-1.0, 0.0, 300.0, bad]).is_err());
        }
        let build = |use_sig: bool| {
            let mut b = SystemBuilder::new();
            let (m0, mass) = (
                Motion {
                    position_m: 0.0,
                    velocity_m_s: 2.0,
                },
                Mass::new(1.0).unwrap(),
            );
            let m = b.add_with_tangent("m", mass, m0).unwrap();
            let sc = if use_sig {
                Slider::loaded_by_signal(0.4, 0.3)
            } else {
                Slider::new(10.0, 0.4, 0.3)
            };
            let sl = b
                .add_on_line("s", sc.unwrap(), Slip::sliding(-2.0))
                .unwrap();
            let p = |b: &SystemBuilder, n, k| b.port(n, k).unwrap();
            b.bond(p(&b, sl, "a"), p(&b, m, "p")).unwrap();
            b.boundary(p(&b, sl, "b"), Box::new(|_, _, v| v[0] = 0.0))
                .unwrap();
            b.boundary(p(&b, sl, "heat"), Box::new(|_, _, v| v[0] = 300.0))
                .unwrap();
            if use_sig {
                let f = Box::new(|t, _, v: &mut [f64]| v[0] = if t < 0.05 { 10.0 } else { 0.0 });
                b.scheduled_boundary(p(&b, sl, "normal_n"), f, vec![0.05])
                    .unwrap();
            }
            (b.seal().unwrap(), m)
        };
        let ((mut ss, ms), (mut sc, mc)) = (build(true), build(false));
        let diss = |s: &System| s.ledger().nodes[1].dissipated_j.value();
        let vel = |s: &System, m| s.state::<Motion>(m).unwrap().velocity_m_s;
        for _ in 0..50 {
            ss.step(1.0e-3, cx).unwrap();
            sc.step(1.0e-3, cx).unwrap();
            assert!((vel(&ss, ms) - vel(&sc, mc)).abs() <= 1.0e-13);
        }
        let (v0, q0) = (vel(&ss, ms), diss(&ss));
        for _ in 0..50 {
            ss.step(1.0e-3, cx).unwrap();
            assert!((vel(&ss, ms) - v0).abs() <= 1.0e-14 && diss(&ss).to_bits() == q0.to_bits());
        }
    });
}

fn add_iface(
    cat: &mut EntityCatalog,
    root: EntityId,
    name: &str,
    (fa, fb): (EntityId, EntityId),
    kind: d::DraftInterfaceKind,
) -> d::DraftInterfaceDecl {
    let pair = InterfacePair::ordered(fa, fb);
    let entity = cat
        .declare(EntityDeclaration::interface(root, name, pair))
        .unwrap();
    d::DraftInterfaceDecl {
        name: name.into(),
        entity,
        feature_a: fa,
        feature_b: fb,
        kind,
    }
}

#[test]
fn normal_load_and_hertz_lowering_refusals() {
    let mut cards = CardStore::new();
    let mat = cards.insert_material(unit_mat(1.0, 100.0));
    let fc = cards
        .insert_interface("dry-friction", friction_card(0.5, 0.3, 50.0, 5.0))
        .unwrap();
    let sa = ("a", d::DraftPartKind::RigidBody, mat, false, None);
    let sb = ("b", d::DraftPartKind::DrivenBody, mat, false, None);
    let hk = || d::DraftInterfaceKind::HertzContact {
        state: vec![
            kv("maximum_force", 1e4, F_DIM),
            kv("maximum_penetration", 1e-3, L),
        ],
    };
    let mc = |state| d::DraftInterfaceKind::MechanicalContact {
        law_card: fc,
        third_body: None,
        state,
    };
    let s = |c: &mut EntityCatalog, p, n| c.declare(EntityDeclaration::surface(p, n)).unwrap();
    let run = |mk: &dyn Fn(&mut EntityCatalog, EntityId, [EntityId; 2]) -> Vec<_>| {
        lower_draft_ifaces((&sa, &sb), (1e-3, 1e-2), &cards, |c, r, ps| {
            (vec![], mk(c, r, ps), vec![])
        })
        .err()
        .expect("refusal")
    };
    let err = run(&|c, r, [pa, pb]| {
        let pair = (s(c, pa, "sa"), s(c, pb, "sb"));
        vec![add_iface(c, r, "hertz", pair, hk())]
    });
    assert!(matches!(
        err,
        LoweringError::MissingInterfaceConstructor {
            ref interface,
            kind: InterfaceKindTag::HertzContact,
            ref card_family,
        } if interface == "hertz" && card_family == "hertz-contact"
    ));
    let n10 = || kv("normal_force", 10.0, F_DIM);
    let err = run(&|c, r, [pa, pb]| {
        let p1 = (s(c, pa, "a1"), s(c, pb, "b1"));
        let p2 = (s(c, pa, "a2"), s(c, pb, "b2"));
        let (so, sh, th) = (s(c, pa, "ao"), s(c, pa, "ah"), s(c, pb, "bh"));
        vec![
            add_iface(c, r, "pad1", p1, mc(vec![n10()])),
            add_iface(c, r, "pad2", p2, mc(vec![n10()])),
            add_iface(c, r, "nload", (so, sh), d::DraftInterfaceKind::NormalLoad),
            add_iface(c, r, "hertz", (sh, th), hk()),
        ]
    });
    assert!(
        matches!(err, LoweringError::Invalid { ref item, .. } if item == "interface:nload"),
        "{err}"
    );
    let r10 = || kv("reference_normal_force", 10.0, F_DIM);
    for (p_st, with_nl, needle) in [
        (vec![n10()], true, "both"),
        (vec![], false, "requires"),
        (vec![], true, "reference_normal_force"),
        (vec![n10(), r10()], false, "reference_normal_force"),
    ] {
        let err = run(&|c, r, [pa, pb]| {
            let (fa, fb) = (s(c, pa, "sa"), s(c, pb, "sb"));
            let (sh, th) = (s(c, pa, "sh"), s(c, pb, "th"));
            let mut v = vec![add_iface(c, r, "pad", (fa, fb), mc(p_st.clone()))];
            if with_nl {
                let nl = d::DraftInterfaceKind::NormalLoad;
                v.push(add_iface(c, r, "nload", (fa, sh), nl));
                v.push(add_iface(c, r, "hertz", (sh, th), hk()));
            }
            v
        });
        assert!(
            matches!(err, LoweringError::Invalid { ref item, ref detail }
                if item == "interface:pad" && detail.contains(needle)),
            "{err}"
        );
    }
}

#[cfg(feature = "slices")]
mod hertz {
    use std::collections::BTreeMap;

    use fs_dcontact::Obstacle;
    use fs_dcontact::component::{ContactInterface, ContactSolve, ContactState};
    use fs_ports::PortKind;
    use fs_system::elements::{Force, HeatCapacity, Waveform};
    use fs_system::lowering::owners::register_owner_defaults;
    use fs_system::{NodeId, RouteTemperature};

    use super::*;

    const T0: f64 = 293.15;
    const R_SPHERE: f64 = 0.05;
    const F_PRE: f64 = 5.0;
    const MU_S: f64 = 0.40;
    const MU_K: f64 = 0.25;
    const V_T0: f64 = 2.5e-3;
    const BETA_316: f64 = 0.363_297;
    const SEED: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/../../data/matdb/seed-v1/");

    type Kv = (ParameterKey, QtyAny);
    type MatSpec = (ContentHash, Vec<Kv>);

    fn load_seed_mat(dir: &str) -> (m::MaterialCard, BTreeMap<String, f64>) {
        let txt = std::fs::read_to_string(format!("{SEED}{dir}/properties.tsv")).unwrap();
        let (mut rows, mut vals, u) = (Vec::new(), BTreeMap::new(), ValidityDomain::unconstrained);
        for line in txt.lines().skip(1) {
            let c: Vec<&str> = line.split('\t').collect();
            if c[0] == "scalar" {
                let dim = match c[3] {
                    "density" => Dims([-3, 1, 0, 0, 0, 0]),
                    "young_modulus" => Dims([-1, 1, -2, 0, 0, 0]),
                    "poisson_ratio" => Dims::NONE,
                    "specific_heat_capacity" => Dims([2, 0, -2, -1, 0, 0]),
                    "thermal_conductivity" => Dims([1, 1, -3, -1, 0, 0]),
                    _ => continue,
                };
                let sc = if c[5] == "GPa" {
                    1e9
                } else if c[5] == "MPa" {
                    1e6
                } else {
                    1.0
                };
                rows.push((c[1], c[3], c[4].parse::<f64>().unwrap() * sc, dim));
            } else if c[0] == "validity" {
                let (lo, hi) = (c[3].parse().unwrap(), c[4].parse().unwrap());
                let e = vals.entry(c[1]).or_insert_with(u);
                *e = e.clone().with(c[2], lo, hi);
            }
        }
        let (mut cs, mut map) = (m::ClaimSet::new(), BTreeMap::new());
        for (id, prop, val, dim) in rows {
            cs.insert_claim(claim(prop, val, dim, vals.remove(id).unwrap_or_else(u)))
                .unwrap();
            map.insert(prop.to_owned(), val);
        }
        let sid = m::MaterialStateId {
            chemistry: dir.into(),
            phase: "s".into(),
            process: "r".into(),
            revision: 0,
        };
        (m::MaterialCard::assemble(sid, cs, Vec::new()).unwrap(), map)
    }

    fn qp_316() -> Vec<Kv> {
        vec![
            kv("T", T0, TEMP),
            kv("source-pressure-known", 0.0, Dims::NONE),
        ]
    }

    fn qp_a36() -> Vec<Kv> {
        const KEYS: [(&str, f64); 5] = [
            ("engineering_reference_model", 1.0),
            ("source_specimens_matched", 0.0),
            ("source_test_temperature_known", 0.0),
            ("source_test_method_known", 0.0),
            ("source_hot_rolled", 1.0),
        ];
        KEYS.map(|(k, v)| kv(k, v, Dims::NONE)).into()
    }

    fn shaped_part(s: &Side, e: EntityId, g: d::PrimitiveShape, qp: Vec<Kv>) -> d::DraftPartDecl {
        let mut p = make_part(s, e);
        p.geometry = d::DraftGeometrySource::Primitive(g);
        p.materials[0].query_point = qp;
        p
    }

    fn box_g() -> d::PrimitiveShape {
        d::PrimitiveShape::Box {
            lx: QtyAny::new(0.1, L),
            ly: QtyAny::new(0.1, L),
            lz: QtyAny::new(0.02, L),
        }
    }

    fn lower_two_body(
        cards: &CardStore,
        (fcard, a36): (ContentHash, ContentHash),
        (ga, (ma, qa)): (d::PrimitiveShape, MatSpec),
        h_st: Vec<Kv>,
        v0: f64,
    ) -> Result<LoweredSystem, LoweringError> {
        let mut cat = EntityCatalog::new();
        let root = cat.declare(EntityDeclaration::assembly("r")).unwrap();
        let pa = cat.declare(EntityDeclaration::part(root, "a")).unwrap();
        let pb = cat.declare(EntityDeclaration::part(root, "b")).unwrap();
        let mut sf = |p, n| cat.declare(EntityDeclaration::surface(p, n)).unwrap();
        let (fah, fbh, fap, fbp) = (sf(pa, "ah"), sf(pb, "bh"), sf(pa, "ap"), sf(pb, "bp"));
        let st = vec![
            kv("reference_normal_force", F_PRE, F_DIM),
            kv("heat_partition_a", 0.5, Dims::NONE),
        ];
        let mc = d::DraftInterfaceKind::MechanicalContact {
            law_card: fcard,
            third_body: None,
            state: st,
        };
        let hk = d::DraftInterfaceKind::HertzContact { state: h_st };
        let nl = d::DraftInterfaceKind::NormalLoad;
        let ifaces = vec![
            add_iface(&mut cat, root, "h", (fah, fbh), hk),
            add_iface(&mut cat, root, "p", (fap, fbp), mc),
            add_iface(&mut cat, root, "nl", (fap, fah), nl),
        ];
        let rb = d::DraftPartKind::RigidBody;
        let sa = ("a", rb, ma, true, Some(init_st(v0, Some(T0))));
        let sb = ("b", rb, a36, true, Some(init_st(0.0, Some(T0))));
        let parts = vec![
            shaped_part(&sa, pa, ga, qa),
            shaped_part(&sb, pb, box_g(), qp_a36()),
        ];
        let mut reg = LoweringRegistry::with_kernel_defaults();
        register_owner_defaults(&mut reg);
        admit_and_lower(
            (root, cat),
            parts,
            (vec![], ifaces, vec![]),
            (1e-5, 1e-4, T0),
            cards,
            &reg,
        )
    }

    #[test]
    fn hertz_lowering_refusals_and_lossy_shared_thermal_closure() {
        let mut cards = CardStore::new();
        let mut ins = |d| cards.insert_material(load_seed_mat(d).0);
        let s316 = ins("stainless-316-20c-engineering-reference");
        let a36 = ins("steel-a36-onlinemetals-reference");
        let wire = ins("steel-music-wire-string-heptonstall-ligo-p1400229");
        // Synthetic in-test friction card (not Katsivardi).
        let fcard = cards
            .insert_interface("dry-friction", friction_card(MU_S, MU_K, 20.0, 1.0))
            .unwrap();
        let ok_h = || {
            vec![
                kv("maximum_force", 1e6, F_DIM),
                kv("maximum_penetration", 1e-2, L),
            ]
        };
        let sph = d::PrimitiveShape::Sphere {
            radius: QtyAny::new(R_SPHERE, L),
        };
        let run = |g, m, h, v| lower_two_body(&cards, (fcard, a36), (g, m), h, v);
        let w_err = run(sph.clone(), (wire, vec![kv("T", T0, TEMP)]), ok_h(), 0.0)
            .err()
            .unwrap();
        assert!(w_err.to_string().contains("poisson_ratio"), "{w_err}");
        let flat = run(box_g(), (s316, qp_316()), ok_h(), 0.0).err().unwrap();
        assert!(flat.to_string().contains("flat-on-flat"), "{flat}");
        let mut lossy = ok_h();
        lossy.extend([
            kv("restitution", 0.7, Dims::NONE),
            kv("approach_speed", 0.2, V_DIM),
        ]);
        let mut lo = run(sph, (s316, qp_316()), lossy, -0.2).unwrap();
        let json = lo.receipt.canonical_json();
        assert!(
            json.contains("reference_normal_force")
                && json.contains("can leave card normal_force validity")
                && json.contains("state.maximum_penetration")
                && json.contains("state.approach_speed")
                && json.contains("\"gap\"")
                && json.contains("not declared: 0")
                && !json.contains("state.gap")
                && !json.contains("state.initial_closure")
        );
        with_cx(|cx| {
            for _ in 0..10 {
                lo.system.step(1e-5, cx).unwrap();
            }
        });
        let nodes = &lo.system.ledger().nodes;
        let by = |name| nodes.iter().find(|n| n.name == name).unwrap();
        let (qh, qp) = (
            by("interface.h").dissipated_j.value(),
            by("interface.p").dissipated_j.value(),
        );
        let du = by("a.thermal").stored_change_j.value() + by("b.thermal").stored_change_j.value();
        assert!(qh > 0.0 && qp > 0.0 && (du - (qh + qp)).abs() <= 1e-12);
        assert!(lo.system.ledger().energy_residual_j().abs() <= 1e-12);
    }

    fn hertz_params(
        stainless: &BTreeMap<String, f64>,
        carbon: &BTreeMap<String, f64>,
    ) -> (f64, f64, f64) {
        let (ea, nua) = (stainless["young_modulus"], stainless["poisson_ratio"]);
        let (eb, nub) = (carbon["young_modulus"], carbon["poisson_ratio"]);
        let e_star = 1.0 / ((1.0 - nua * nua) / ea + (1.0 - nub * nub) / eb);
        let k_h = (4.0 / 3.0) * e_star * R_SPHERE.sqrt();
        (e_star, k_h, (F_PRE / k_h).powf(2.0 / 3.0))
    }

    fn oracle1_draft(
        cards: &CardStore,
        [s316, a36, fcard]: [ContentHash; 3],
        delta_eq: f64,
        nl_first: bool,
    ) -> LoweredSystem {
        let mut cat = EntityCatalog::new();
        let root = cat.declare(EntityDeclaration::assembly("o1")).unwrap();
        let mut mk_p = |n| cat.declare(EntityDeclaration::part(root, n)).unwrap();
        let (nb, nf, ta, tb) = (mk_p("nb"), mk_p("nf"), mk_p("ta"), mk_p("tb"));
        let mut sf = |p, n| cat.declare(EntityDeclaration::surface(p, n)).unwrap();
        let (ball_face, flat_face) = (sf(nb, "h"), sf(nf, "h"));
        let (sled_face, base_face) = (sf(ta, "p"), sf(tb, "p"));
        let hk = d::DraftInterfaceKind::HertzContact {
            state: vec![
                kv("maximum_force", 1e6, F_DIM),
                kv("maximum_penetration", 1e-2, L),
                kv("gap", 0.0, L),
                kv("initial_closure", delta_eq, L),
            ],
        };
        let mc = d::DraftInterfaceKind::MechanicalContact {
            law_card: fcard,
            third_body: None,
            state: vec![
                kv("reference_normal_force", F_PRE, F_DIM),
                kv("heat_partition_a", BETA_316, Dims::NONE),
            ],
        };
        let mut ifaces = Vec::new();
        let nl = || d::DraftInterfaceKind::NormalLoad;
        if nl_first {
            ifaces.push(add_iface(
                &mut cat,
                root,
                "nl",
                (sled_face, ball_face),
                nl(),
            ));
        }
        ifaces.push(add_iface(&mut cat, root, "h", (ball_face, flat_face), hk));
        ifaces.push(add_iface(&mut cat, root, "p", (sled_face, base_face), mc));
        if !nl_first {
            ifaces.push(add_iface(
                &mut cat,
                root,
                "nl",
                (sled_face, ball_face),
                nl(),
            ));
        }
        let sph = d::PrimitiveShape::Sphere {
            radius: QtyAny::new(R_SPHERE, L),
        };
        let (rb, drv) = (d::DraftPartKind::RigidBody, d::DraftPartKind::DrivenBody);
        let ball_side = ("nb", rb, s316, false, Some(init_st(0.0, None)));
        let flat_side = ("nf", drv, a36, false, None);
        let sled_side = ("ta", rb, s316, true, Some(init_st(V_T0, Some(T0))));
        let base_side = ("tb", rb, a36, true, Some(init_st(0.0, Some(T0))));
        let parts = vec![
            shaped_part(&ball_side, nb, sph.clone(), qp_316()),
            shaped_part(&flat_side, nf, box_g(), qp_a36()),
            shaped_part(&sled_side, ta, sph, qp_316()),
            shaped_part(&base_side, tb, box_g(), qp_a36()),
        ];
        let src = d::SourceDecl {
            name: "preload".into(),
            target: d::SourceTarget::Port {
                entity: nb,
                port: "mechanics".into(),
            },
            history: d::TimeHistory::Constant(QtyAny::new(-F_PRE, F_DIM)),
        };
        let mut reg = LoweringRegistry::with_kernel_defaults();
        register_owner_defaults(&mut reg);
        admit_and_lower(
            (root, cat),
            parts,
            (vec![], ifaces, vec![src]),
            (1e-4, 4e-3, T0),
            cards,
            &reg,
        )
        .unwrap()
    }

    fn hertz_comp(e_star: f64, k_h: f64) -> ContactInterface {
        let prov = format!("hertz-contact:h:E*={e_star:.6e},R*={R_SPHERE:.6e}");
        let law = Obstacle::new(vec![-1.0], 1, 1, vec![0.0], vec![1.0], k_h, 1.5, prov).unwrap();
        let solve = ContactSolve {
            max_iterations: 64,
            maximum_force_n: 1e6,
            maximum_penetration_m: 1e-2,
            force_absolute_tolerance_n: 1e-12,
            force_relative_tolerance: 1e-12,
        };
        ContactInterface::new(law, solve)
            .unwrap()
            .publishing_reaction()
            .unwrap()
    }

    fn oracle1_kernel(
        stainless: &BTreeMap<String, f64>,
        carbon: &BTreeMap<String, f64>,
        delta_eq: f64,
        live: bool,
    ) -> (System, [Option<NodeId>; 6]) {
        let (e_star, k_h, _) = hertz_params(stainless, carbon);
        let (v_s, v_b) = ((4.0 / 3.0) * PI * R_SPHERE.powi(3), 0.1 * 0.1 * 0.02);
        let (m_a, m_b) = (stainless["density"] * v_s, carbon["density"] * v_b);
        let (c_a, c_b) = (
            m_a * stainless["specific_heat_capacity"],
            m_b * carbon["specific_heat_capacity"],
        );
        let mut b = SystemBuilder::new();
        let m0 = |v| Motion {
            position_m: 0.0,
            velocity_m_s: v,
        };
        let nb = live.then(|| {
            b.add_with_tangent("nb.mechanics", Mass::new(m_a).unwrap(), m0(0.0))
                .unwrap()
        });
        let sled_m = b
            .add_with_tangent("ta.mechanics", Mass::new(m_a).unwrap(), m0(V_T0))
            .unwrap();
        let sled_th = b
            .add_with_tangent("ta.thermal", HeatCapacity::new(c_a).unwrap(), T0)
            .unwrap();
        let base_m = b
            .add_with_tangent("tb.mechanics", Mass::new(m_b).unwrap(), m0(0.0))
            .unwrap();
        let base_th = b
            .add_with_tangent("tb.thermal", HeatCapacity::new(c_b).unwrap(), T0)
            .unwrap();
        let h_node = live.then(|| {
            b.add_on_line(
                "interface.h",
                hertz_comp(e_star, k_h),
                ContactState::at(delta_eq),
            )
            .unwrap()
        });
        let sl = if live {
            Slider::loaded_by_signal(MU_S, MU_K).unwrap()
        } else {
            Slider::new(F_PRE, MU_S, MU_K).unwrap()
        };
        let p_node = b
            .add_on_line("interface.p", sl, Slip::sliding(-V_T0))
            .unwrap();
        let pt = |b: &SystemBuilder, n, k| b.port(n, k).unwrap();
        if let (Some(nb), Some(hn)) = (nb, h_node) {
            let wf = Waveform::constant(-F_PRE).unwrap();
            let f = Force::in_domain(PortKind::MechanicalForceVelocity, wf);
            let fn_id = b
                .add_with_tangent_and_stage_residual("source.preload", f, SourceAccount::default())
                .unwrap();
            b.common_flow(pt(&b, nb, "p"), &[pt(&b, hn, "a"), pt(&b, fn_id, "p")])
                .unwrap();
            b.boundary(pt(&b, hn, "b"), Box::new(|_, _, o| o.fill(0.0)))
                .unwrap();
            b.signal(pt(&b, hn, "reaction_n"), pt(&b, p_node, "normal_n"))
                .unwrap();
        }
        b.common_flow(pt(&b, sled_m, "p"), &[pt(&b, p_node, "a")])
            .unwrap();
        b.common_flow(pt(&b, base_m, "p"), &[pt(&b, p_node, "b")])
            .unwrap();
        let sinks = [
            (pt(&b, sled_th, "heat"), BETA_316),
            (pt(&b, base_th, "heat"), 1.0 - BETA_316),
        ];
        b.heat_route(
            pt(&b, p_node, "heat"),
            &sinks,
            None,
            RouteTemperature::SinkEquivalent,
        )
        .unwrap();
        let p_heat = pt(&b, p_node, "heat");
        let mut sys = b.seal().unwrap();
        let _ = sys.guess(p_heat, &[T0]);
        let ids = [
            h_node,
            Some(p_node),
            Some(sled_m),
            Some(base_m),
            Some(sled_th),
            Some(base_th),
        ];
        (sys, ids)
    }

    fn snap_step(s: &System, [h, p, tam, tbm, tat, tbt]: [NodeId; 6]) -> ([u64; 8], Slip) {
        let cs = *s.state::<ContactState>(h).unwrap();
        let ma = *s.state::<Motion>(tam).unwrap();
        let mb = *s.state::<Motion>(tbm).unwrap();
        let ta = *s.state::<f64>(tat).unwrap();
        let tb = *s.state::<f64>(tbt).unwrap();
        let sl = *s.state::<Slip>(p).unwrap();
        let bits = [
            cs.closure_m,
            cs.reaction_n,
            ma.position_m,
            ma.velocity_m_s,
            mb.position_m,
            mb.velocity_m_s,
            ta,
            tb,
        ]
        .map(f64::to_bits);
        (bits, sl)
    }

    #[test]
    fn hertz_preload_bit_for_bit_parity_and_order_independence() {
        let mut cards = CardStore::new();
        let (ss_card, stainless) = load_seed_mat("stainless-316-20c-engineering-reference");
        let (carbon_card, carbon) = load_seed_mat("steel-a36-onlinemetals-reference");
        let (s316, a36) = (
            cards.insert_material(ss_card),
            cards.insert_material(carbon_card),
        );
        let fcard = cards
            .insert_interface("dry-friction", friction_card(MU_S, MU_K, 20.0, 1.0))
            .unwrap();
        let (_, _, delta_eq) = hertz_params(&stainless, &carbon);
        let mut lo_after = oracle1_draft(&cards, [s316, a36, fcard], delta_eq, false);
        let mut lo_before = oracle1_draft(&cards, [s316, a36, fcard], delta_eq, true);
        let (mut kern, k_ids) = oracle1_kernel(&stainless, &carbon, delta_eq, true);
        let (mut const_n, c_ids) = oracle1_kernel(&stainless, &carbon, delta_eq, false);
        let rcpt = lo_after.receipt.canonical_json();
        for needle in [
            "effective_modulus",
            "effective_radius",
            "state.gap",
            "state.initial_closure",
            "reference_normal_force",
            "can leave card normal_force validity",
            "heat_partition_a",
            "heat_partition_b",
        ] {
            assert!(rcpt.contains(needle), "missing {needle} in receipt");
        }
        let lo_ids = |lo: &LoweredSystem| {
            [
                "interface.h",
                "interface.p",
                "ta.mechanics",
                "tb.mechanics",
                "ta.thermal",
                "tb.thermal",
            ]
            .map(|n| lo.node_map[n])
        };
        let (a_ids, b_ids) = (lo_ids(&lo_after), lo_ids(&lo_before));
        let k6 = k_ids.map(Option::unwrap);
        with_cx(|cx| {
            for step in 1..=40 {
                lo_after.system.step(1e-4, cx).unwrap();
                lo_before.system.step(1e-4, cx).unwrap();
                kern.step(1e-4, cx).unwrap();
                const_n.step(1e-4, cx).unwrap();
                let sa = snap_step(&lo_after.system, a_ids);
                assert_eq!(
                    sa,
                    snap_step(&lo_before.system, b_ids),
                    "order diff at step {step}"
                );
                assert_eq!(sa, snap_step(&kern, k6), "kernel diff at step {step}");
                let cs = *lo_after.system.state::<ContactState>(a_ids[0]).unwrap();
                assert!(
                    (cs.reaction_n - F_PRE).abs() <= 1e-11,
                    "step {step} R={}",
                    cs.reaction_n
                );
                assert!(
                    (cs.closure_m - delta_eq).abs() <= 1e-15,
                    "step {step} x={}",
                    cs.closure_m
                );
                let va = lo_after
                    .system
                    .state::<Motion>(a_ids[2])
                    .unwrap()
                    .velocity_m_s;
                let vc = const_n
                    .state::<Motion>(c_ids[2].unwrap())
                    .unwrap()
                    .velocity_m_s;
                assert!((va - vc).abs() <= 1e-12, "step {step} va={va} vc={vc}");
            }
        });
        assert!(lo_after.system.state::<Slip>(a_ids[1]).unwrap().is_stuck());
        let nodes = &lo_after.system.ledger().nodes;
        let by = |name| nodes.iter().find(|n| n.name == name).unwrap();
        let (dua, dub) = (
            by("ta.thermal").stored_change_j.value(),
            by("tb.thermal").stored_change_j.value(),
        );
        let q = dua + dub;
        assert!(q > 0.0 && ((dua / q) - BETA_316).abs() <= 1e-12);
        assert!(lo_after.system.ledger().energy_residual_j().abs() <= 1e-12);
    }
}
