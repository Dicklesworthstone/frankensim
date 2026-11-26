//! Adiabatic Rüchardt and finite-`G` thermal-wall modal oracles (`3gc12.6.2`).
//!
//! Uses `examples/assemblies/packs/nitrogen-n2-mcbride2002-lemmon2004.fsmcdpk`
//! and `steel-music-wire-heptonstall-ligo-p1400229.fsmcdpk`.

use core::f64::consts::PI;
use std::{collections::BTreeSet, path::Path};

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_gas::component::ClosedGasParcel;
use fs_matdb::{NormalizedInterfacePack, NormalizedMaterialCardPack};
use fs_qty::{Dims, QtyAny};
use fs_scenario::assembly::draft::{
    DraftPartKind as Pk, PhysicsDomain as Pd, PrimitiveShape as Ps, *,
};
use fs_scenario::assembly::{MaterialRef, ParameterKey, ReferenceState};
use fs_scenario::{
    EntityCatalog, EntityDeclaration as Ed, EntityId, FrameId, InterfacePair, Scenario,
};
use fs_system::elements::Motion;
use fs_system::lowering::owners::register_owner_defaults;
use fs_system::lowering::{CardStore, LoweredSystem, LoweringError, LoweringRegistry};
use fs_system::{
    AnalysisError, AnalysisOutcome, AnalysisRequest, ModalConfig, NodeId, SolveConfig, StateMap,
    System,
};

const LEN: Dims = Dims([1, 0, 0, 0, 0, 0]);
const AREA: Dims = Dims([2, 0, 0, 0, 0, 0]);
const SEC: Dims = Dims([0, 0, 1, 0, 0, 0]);
const TEMP: Dims = Dims([0, 0, 0, 1, 0, 0]);
const FORCE: Dims = Dims([1, 1, -2, 0, 0, 0]);
const PRESS: Dims = Dims([-1, 1, -2, 0, 0, 0]);

struct GasPistonMap(NodeId, NodeId);

impl StateMap for GasPistonMap {
    fn dim(&self) -> usize {
        3
    }
    fn read(&self, sys: &System, out: &mut [f64]) -> Result<(), AnalysisError> {
        let inv = |what| AnalysisError::Invalid { what };
        let m = sys.state::<Motion>(self.0).ok_or(inv("piston"))?;
        let g = sys.state::<ClosedGasParcel>(self.1).ok_or(inv("gas"))?;
        out.copy_from_slice(&[m.velocity_m_s, g.volume_m3(), g.state().state().temperature]);
        Ok(())
    }
    fn write(&self, sys: &mut System, s: &[f64]) -> Result<(), AnalysisError> {
        let inv = |what| AnalysisError::Invalid { what };
        let mut m = *sys.state::<Motion>(self.0).ok_or(inv("piston"))?;
        m.velocity_m_s = s[0];
        sys.set_accepted_state(self.0, m)?;
        let g = sys.state::<ClosedGasParcel>(self.1).ok_or(inv("gas"))?;
        let next = g
            .at_volume_temperature(s[1], s[2])
            .map_err(|_| inv("gas"))?;
        sys.set_accepted_state(self.1, next)?;
        Ok(())
    }
}

fn load_pack(name: &str) -> fs_matdb::MaterialCard {
    let dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/assemblies/packs");
    let pack = NormalizedMaterialCardPack::from_bytes(&std::fs::read(dir.join(name)).unwrap());
    pack.unwrap().card().clone()
}

fn load_atkins_pack() -> fs_matdb::InterfaceSystemCard {
    let p = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../examples/system/atkins1964/al6061t6-lathe-0p30um-vacuum.fsintpk");
    NormalizedInterfacePack::from_bytes(&std::fs::read(p).unwrap())
        .unwrap()
        .card()
        .clone()
}

fn with_cx<T>(seed: u64, f: impl FnOnce(&Cx<'_>) -> T) -> T {
    let gate = CancelGate::new_clock_free();
    let k = StreamKey {
        seed,
        kernel_id: 1,
        tile: 0,
        iteration: 0,
    };
    ArenaPool::new(ArenaConfig::default()).scope(|a| {
        f(&Cx::new(
            &gate,
            a,
            k,
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

fn mk_part(
    name: &str,
    entity: EntityId,
    shape: Ps,
    spec: (Pk, fs_blake3::ContentHash, (&str, QtyAny), BTreeSet<Pd>),
) -> DraftPartDecl {
    let (kind, card, (qk, qv), physics) = spec;
    let query_point = vec![(ParameterKey::new(qk, qv.dims), qv)];
    let error_budget = FidelityErrorBudget {
        relative_tol: 1e-6,
        absolute_tol: QtyAny::dimensionless(1e-12),
        max_work_units: Some(128),
    };
    DraftPartDecl {
        name: name.into(),
        entity,
        kind,
        geometry: DraftGeometrySource::Primitive(shape),
        materials: vec![MaterialRef {
            card,
            query_point,
            frame: FrameId(0),
            layer: None,
            unknown: vec![],
        }],
        physics,
        fidelity: FidelityPolicy {
            rungs: vec![FidelityRung::Lumped],
            error_budget,
        },
        uncertainty: vec![],
        initial: None,
    }
}

fn add_thermal_wall(
    cat: &mut EntityCatalog,
    cards: &mut CardStore,
    ids: (EntityId, EntityId),
    st: fs_blake3::ContentHash,
    bx: Ps,
    wall: (f64, bool),
) -> (DraftPartDecl, DraftInterfaceDecl) {
    let ((root, e_g), (wall_area, both_ground), q) = (ids, wall, QtyAny::new);
    let e_fr = cat.declare(Ed::part(root, "frame")).unwrap();
    let dom = BTreeSet::from([Pd::Mechanics]);
    let frame = mk_part(
        "frame",
        e_fr,
        bx,
        (Pk::DrivenBody, st, ("T", q(293.15, TEMP)), dom),
    );
    let wc = cards
        .insert_interface("thermal-contact-resistance", load_atkins_pack())
        .unwrap();
    let (s_gw, s_a) = (
        cat.declare(Ed::surface(e_g, "gw")).unwrap(),
        cat.declare(Ed::surface(e_fr, "amb")).unwrap(),
    );
    let fa = if both_ground {
        cat.declare(Ed::surface(e_fr, "amb2")).unwrap()
    } else {
        s_gw
    };
    let e_w = cat
        .declare(Ed::interface(
            root,
            "wall",
            InterfacePair::unordered(fa, s_a),
        ))
        .unwrap();
    let state = vec![
        (ParameterKey::new("area", AREA), q(wall_area, AREA)),
        (
            ParameterKey::new("normal_pressure", PRESS),
            q(313_000.0, PRESS),
        ),
    ];
    let kind = DraftInterfaceKind::ThermalContact {
        law_card: wc,
        third_body: None,
        state,
    };
    (
        frame,
        DraftInterfaceDecl {
            name: "wall".into(),
            entity: e_w,
            feature_a: fa,
            feature_b: s_a,
            kind,
        },
    )
}

fn lower_piston(
    with_ambient: bool,
    face_state: Vec<(ParameterKey, QtyAny)>,
    wall: Option<(f64, bool)>,
) -> Result<(LoweredSystem, f64, f64, f64, f64), LoweringError> {
    let (r, h, t0, p0) = (0.15_f64, 0.20_f64, 300.0_f64, 90_000.0_f64);
    let (area, v0, q) = (PI * r * r, PI * r * r * h, QtyAny::new);
    let mut cards = CardStore::new();
    let n2 = cards.insert_material(load_pack("nitrogen-n2-mcbride2002-lemmon2004.fsmcdpk"));
    let st = cards.insert_material(load_pack(
        "steel-music-wire-heptonstall-ligo-p1400229.fsmcdpk",
    ));
    let mut cat = EntityCatalog::new();
    let (root, e_p, e_g, s_p, s_g, e_f) = {
        let mut d = |e| cat.declare(e).unwrap();
        let root = d(Ed::assembly("gp"));
        let (e_p, e_g) = (d(Ed::part(root, "piston")), d(Ed::part(root, "gas")));
        let (s_p, s_g) = (d(Ed::surface(e_p, "fp")), d(Ed::surface(e_g, "fg")));
        let e_f = d(Ed::interface(
            root,
            "pf",
            InterfacePair::unordered(s_p, s_g),
        ));
        (root, e_p, e_g, s_p, s_g, e_f)
    };
    let d = q(0.1, LEN);
    let bx = Ps::Box {
        lx: d,
        ly: d,
        lz: d,
    };
    let cyl = Ps::Cylinder {
        radius: q(r, LEN),
        height: q(h, LEN),
    };
    let md = BTreeSet::from([Pd::Mechanics]);
    let mut parts = vec![mk_part(
        "piston",
        e_p,
        bx.clone(),
        (Pk::RigidBody, st, ("T", q(293.15, TEMP)), md),
    )];
    let mut g_dom = BTreeSet::from([Pd::Fluid]);
    let mut interfaces = vec![DraftInterfaceDecl {
        name: "pf".into(),
        entity: e_f,
        feature_a: s_p,
        feature_b: s_g,
        kind: DraftInterfaceKind::PistonFace { state: face_state },
    }];
    if let Some(w) = wall {
        g_dom.insert(Pd::Thermal);
        let (frame, iface) = add_thermal_wall(&mut cat, &mut cards, (root, e_g), st, bx, w);
        parts.push(frame);
        interfaces.push(iface);
    }
    parts.push(mk_part(
        "gas",
        e_g,
        cyl,
        (Pk::Cavity, n2, ("pressure", q(p0, PRESS)), g_dom),
    ));
    let target = SourceTarget::Port {
        entity: e_p,
        port: "mechanics".into(),
    };
    let history = TimeHistory::Constant(q(-p0 * area, FORCE));
    let sources = with_ambient
        .then(|| SourceDecl {
            name: "ambient_pressure".into(),
            target,
            history,
        })
        .into_iter()
        .collect();
    let mut env = DraftEnvironment::earth_lab_air();
    (env.base.ambient_temperature, env.base.ambient_pressure) = (q(t0, TEMP), q(p0, PRESS));
    let sc = Scenario::new("gp", 1, env.base.clone());
    let time_plan = TimePlan {
        horizon: q(0.05, SEC),
        step_policy: StepPolicy::Fixed { dt: q(1e-4, SEC) },
        default_sample_rate: q(1000.0, Dims([0, 0, -1, 0, 0, 0])),
        scheme: TimeScheme::Explicit,
    };
    let inp = DraftAssemblyInput {
        name: "gp".into(),
        root,
        catalog: cat,
        reference_state: Some(ReferenceState::Assembled),
        environment: env,
        time_plan,
        parts,
        joints: vec![],
        gear_meshes: vec![],
        interfaces,
        sources,
        controllers: vec![],
        probes: vec![],
    };
    let mut reg = LoweringRegistry::with_kernel_defaults();
    register_owner_defaults(&mut reg);
    let lowered = reg.lower(&DraftAssembly::admit(inp, &sc).unwrap(), &cards)?;
    let ps = &lowered.receipt.components[0].parameter_sources;
    let m_p = ps.iter().find(|p| p.parameter == "mass").unwrap().value;
    Ok((lowered, m_p, area, v0, p0))
}

/// Unique real root in `(-alpha, 0)` and complex `(sigma, omega_d, omega_n)`
/// of `lambda^3 + alpha*lambda^2 + omega_ad^2*lambda + alpha*omega_iso^2 = 0`.
fn solve_three_state_cubic(alpha: f64, w_iso2: f64, w_ad2: f64) -> (f64, f64, f64, f64) {
    let f = |x: f64| ((x + alpha) * x + w_ad2) * x + alpha * w_iso2;
    let df = |x: f64| (3.0 * x + 2.0 * alpha) * x + w_ad2;
    let (mut lo, mut hi) = (-alpha, 0.0);
    for _ in 0..80 {
        let m = 0.5 * (lo + hi);
        if f(m) > 0.0 {
            hi = m;
        } else {
            lo = m;
        }
    }
    let mut r = 0.5 * (lo + hi);
    for _ in 0..8 {
        r -= f(r) / df(r);
    }
    let sigma = 0.5 * (alpha + r);
    let wn2 = w_ad2 + r * (alpha + r);
    (-r, sigma, (wn2 - sigma * sigma).sqrt(), wn2.sqrt())
}

fn run_modal(lowered: &mut LoweredSystem, fd_eps: f64) -> Vec<fs_system::ModalMode> {
    let (piston, gas) = (lowered.node_map["piston"], lowered.node_map["gas"]);
    let cfg = SolveConfig {
        relative_tolerance: 1e-13,
        ..SolveConfig::default()
    };
    lowered.system.set_config(cfg);
    let (map, mut mcfg) = (GasPistonMap(piston, gas), ModalConfig::default());
    mcfg.linearize.fd_eps = fd_eps;
    with_cx(1, |cx| {
        let req = AnalysisRequest::Modal {
            state_map: &map,
            config: mcfg,
        };
        match lowered.system.run_analysis(req, cx).unwrap() {
            AnalysisOutcome::Modal(m) => m.modes,
            _ => unreachable!(),
        }
    })
}

/// Adiabatic Rüchardt small-amplitude eigenfrequency oracle: exact for the model's ideal gas
/// (`omega_ad = sqrt(gamma(T0) * p0 * A^2 / (m * V0))`, with `gamma(T0)` from the same NASA-9
/// `cp`). Real-gas departure is a model-form question, not part of this test.
#[test]
fn adiabatic_gas_piston_matches_exact_ideal_gas_ruchardt_frequency() {
    let (mut lowered, m_p, area, v0, p0) = lower_piston(true, vec![], None).unwrap();
    lowered.receipt.verify_completeness().unwrap();
    let gas = lowered.node_map["gas"];
    let gs = lowered
        .system
        .state::<ClosedGasParcel>(gas)
        .unwrap()
        .state()
        .state();
    let gamma = gs.specific_heat_cp / (gs.specific_heat_cp - gs.specific_gas_constant);
    let omega_ad = (gamma * p0 * area * area / (m_p * v0)).sqrt();
    let modes = run_modal(&mut lowered, ModalConfig::default().linearize.fd_eps);
    assert_eq!(modes.len(), 3);
    assert!(modes[0].natural_frequency_rad_s < 1e-5);
    let ac = &modes[2];
    let rel = (ac.angular_frequency_rad_s - omega_ad).abs() / omega_ad;
    assert!(
        rel <= 1e-6 && ac.damping_ratio.abs() <= 1e-6,
        "rel={rel:.3e}"
    );
}

#[test]
fn piston_face_has_no_hidden_ambient_force_honors_area_sign_and_refuses_bad_state() {
    let step_v = |state| {
        let (mut lowered, m_p, area, _, p0) = lower_piston(false, state, None).unwrap();
        with_cx(2, |cx| lowered.system.step(1e-4, cx).unwrap());
        let p = lowered.node_map["piston"];
        let v = lowered.system.state::<Motion>(p).unwrap().velocity_m_s;
        (v, m_p, area, p0)
    };
    let (v, m_p, area, p0) = step_v(vec![]);
    let v_exp = (p0 * area / m_p) * 1e-4;
    assert!((v - v_exp).abs() / v_exp < 1e-2, "thrust {v} vs {v_exp}");

    let (a_k, s_k) = (
        ParameterKey::new("area", AREA),
        ParameterKey::new("sign", Dims::NONE),
    );
    let rev = vec![
        (a_k.clone(), QtyAny::new(0.5 * area, AREA)),
        (s_k.clone(), QtyAny::dimensionless(-1.0)),
    ];
    let v_rev_exp = -(p0 * 0.5 * area / m_p) * 1e-4;
    assert!((step_v(rev).0 - v_rev_exp).abs() / v_rev_exp.abs() < 1e-2);

    let bad = |s| {
        matches!(
            lower_piston(true, s, None),
            Err(LoweringError::Invalid { .. })
        )
    };
    let bogus = (
        ParameterKey::new("bogus", Dims::NONE),
        QtyAny::dimensionless(1.0),
    );
    let bad_dim = (ParameterKey::new("area", LEN), QtyAny::new(1.0, LEN));
    assert!(bad(vec![bogus]) && bad(vec![(s_k, QtyAny::dimensionless(2.0))]));
    assert!(bad(vec![(a_k, QtyAny::new(-area, AREA))]) && bad(vec![bad_dim]));
}

#[test]
fn finite_conductance_gas_piston_matches_exact_three_state_cubic_across_g() {
    assert!(matches!(
        lower_piston(true, vec![], Some((0.015, true))),
        Err(LoweringError::Invalid { .. })
    ));
    let mut prev_gap = f64::INFINITY;
    // Atkins 1964 R'' = 3.1606e-4 m^2*K/W at 313 kPa -> G ~ 5, 50, 500 W/K.
    for wall_area in [0.001_580_3, 0.015_803, 0.158_03] {
        let (mut lowered, m_p, area, v0, p0) =
            lower_piston(true, vec![], Some((wall_area, false))).unwrap();
        lowered.receipt.verify_completeness().unwrap();
        let wps = &lowered.receipt.interfaces[1].parameter_sources;
        let g_w_k = wps
            .iter()
            .find(|p| p.parameter == "thermal_conductance")
            .unwrap()
            .value;
        let gas = lowered.node_map["gas"];
        let parcel = lowered.system.state::<ClosedGasParcel>(gas).unwrap();
        let gs = parcel.state().state();
        let cv = gs.specific_heat_cp - gs.specific_gas_constant;
        let alpha = g_w_k / (parcel.mass_kg() * cv);
        let w_iso2 = p0 * area * area / (m_p * v0);
        let (r_ex, sig_ex, wd_ex, wn_ex) =
            solve_three_state_cubic(alpha, w_iso2, (gs.specific_heat_cp / cv) * w_iso2);
        let modes = run_modal(&mut lowered, 1e-6);
        let real = modes
            .iter()
            .find(|m| m.angular_frequency_rad_s < 1e-3)
            .unwrap();
        let osc = modes
            .iter()
            .find(|m| m.angular_frequency_rad_s > 1.0)
            .unwrap();
        let re = (real.natural_frequency_rad_s - r_ex).abs() / r_ex;
        let we = (osc.angular_frequency_rad_s - wd_ex).abs() / wd_ex;
        let ze = (osc.damping_ratio - sig_ex / wn_ex).abs() / (sig_ex / wn_ex);
        assert!(
            re <= 1e-6 && we <= 1e-6 && ze <= 1e-6,
            "G={g_w_k}: ({re:.2e},{we:.2e},{ze:.2e})"
        );
        let gap = (osc.angular_frequency_rad_s - w_iso2.sqrt()).abs();
        assert!(gap < prev_gap, "G={g_w_k}: gap {gap} >= {prev_gap}");
        prev_gap = gap;
    }
}
