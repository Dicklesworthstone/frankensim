//! Slice C (`frankensim-rc-1109-root-3gc12.6.3`): dry-friction heating with
//! stick, slip, arrest, reversal, and heat partitioning on the non-smooth
//! `LoadLine` tier via `LoweringRegistry::with_kernel_defaults()`.

use std::f64::consts::PI;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::{ContentHash, hash_bytes};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb as m;
use fs_qty::{Dims, QtyAny};
use fs_scenario::assembly::draft as d;
use fs_scenario::assembly::{MaterialRef, ParameterKey, ReferenceState};
use fs_scenario::{
    EntityCatalog, EntityDeclaration, EntityId, Environment, FrameId, InterfacePair, Scenario,
};
use fs_system::elements::{Motion, Slip, SourceAccount};
use fs_system::{CardStore, LoweredSystem, LoweringError, LoweringRegistry};

const L: Dims = Dims([1, 0, 0, 0, 0, 0]);
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

fn lower_pair(
    (side_a, side_b): (&Side, &Side),
    (card, state): (ContentHash, Vec<(ParameterKey, QtyAny)>),
    (spring_k, belt): (Option<f64>, Option<d::TimeHistory>),
    (dt, t_end): (f64, f64),
    cards: &CardStore,
) -> Result<LoweredSystem, LoweringError> {
    let mut cat = EntityCatalog::new();
    let mut dec = |e| cat.declare(e).unwrap();
    let root = dec(EntityDeclaration::assembly("p"));
    let pa = dec(EntityDeclaration::part(root, side_a.0));
    let pb = dec(EntityDeclaration::part(root, side_b.0));
    let sa = dec(EntityDeclaration::surface(pa, "sa"));
    let sb = dec(EntityDeclaration::surface(pb, "sb"));
    let cid = dec(EntityDeclaration::interface(
        root,
        "pad",
        InterfacePair::ordered(sa, sb),
    ));
    let mut joints = Vec::new();
    if let Some(k) = spring_k {
        joints.push(d::JointDecl {
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
        });
    }
    let mut sources = Vec::new();
    if let Some(history) = belt {
        sources.push(d::SourceDecl {
            name: "belt_v".into(),
            target: d::SourceTarget::Port {
                entity: pb,
                port: "mechanics".into(),
            },
            history,
        });
    }
    let kind = d::DraftInterfaceKind::MechanicalContact {
        law_card: card,
        third_body: None,
        state,
    };
    let iface = d::DraftInterfaceDecl {
        name: "pad".into(),
        entity: cid,
        feature_a: sa,
        feature_b: sb,
        kind,
    };
    let mut env = d::DraftEnvironment::earth_lab_air();
    env.base.ambient_temperature = QtyAny::new(300.0, TEMP);
    let td = Dims([0, 0, 1, 0, 0, 0]);
    let inp = d::DraftAssemblyInput {
        name: "p".into(),
        root,
        catalog: cat,
        reference_state: Some(ReferenceState::Assembled),
        environment: env,
        time_plan: d::TimePlan {
            horizon: QtyAny::new(t_end, td),
            step_policy: d::StepPolicy::Fixed {
                dt: QtyAny::new(dt, td),
            },
            default_sample_rate: QtyAny::new(1.0 / dt, Dims([0, 0, -1, 0, 0, 0])),
            scheme: d::TimeScheme::Explicit,
        },
        parts: vec![make_part(side_a, pa), make_part(side_b, pb)],
        joints,
        gear_meshes: vec![],
        interfaces: vec![iface],
        sources,
        controllers: vec![],
        probes: vec![],
    };
    let scen = Scenario::new("p", 42, Environment::earth_lab());
    LoweringRegistry::with_kernel_defaults()
        .lower(&d::DraftAssembly::admit(inp, &scen).unwrap(), cards)
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
        (fcard, st),
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
    let msg = lower_pair((&blk, &drv), (nasa, st), (None, None), (0.01, 0.1), &cards)
        .err()
        .unwrap()
        .to_string();
    assert!(msg.contains("static-friction-coefficient"), "{msg}");
}
