//! Modal-attachment Hertz contact with strike-station `phi(x, y)` participation (`.11.9.4`):
//! kernel parity + physical signs and shared-feature refusals (`.fsim` drum in `fs-cli`).

use core::num::NonZeroUsize;
use std::collections::BTreeSet;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_dcontact::{Obstacle, component::*};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{MaterialCard, NormalizedMaterialCardPack as Pack};
use fs_modal::{SliceOptions, slice_window, time::ModalAcousticState};
use fs_nlmodal::component::{ModalAcousticTimeBudget, ModalAttachment, ModalStructure};
use fs_plate::shell::{head::*, profile::ProfileBudget};
use fs_ports::{PortKind, StableId};
use fs_qty::{Dims, QtyAny};
use fs_scenario::{EntityDeclaration as Ed, assembly::draft::*, assembly::*, *};
use fs_system::{elements::*, lowering::owners::register_owner_defaults, *};

const MYLAR: &[u8] =
    include_bytes!("../../fs-vkplate/packs/head-film-pet-mylar-dupont-matsuzawa1959.fsmatpk");
const HICKORY: &[u8] =
    include_bytes!("../../fs-vkplate/packs/stick-hickory-shagbark-fpl-gtr282.fsmatpk");
const P: [f64; 5] = [0.1778, 0.00025, 2500.0, 293.15, 0.015];
const DL: [i8; 6] = [1, 0, 0, 0, 0, 0];
const DF: [i8; 6] = [1, 1, -2, 0, 0, 0];
const DT: [i8; 6] = [0, 0, 0, 1, 0, 0];

fn q(v: f64, d: [i8; 6]) -> QtyAny {
    QtyAny::new(v, Dims(d))
}
fn kv(n: &str, d: [i8; 6], v: f64) -> (ParameterKey, QtyAny) {
    (ParameterKey::new(n, Dims(d)), q(v, d))
}
fn mq(s: &System, n: NodeId) -> &[ModalAcousticState] {
    s.state::<Vec<ModalAcousticState>>(n).unwrap()
}
fn pv(s: &[lowering::ParameterProvenance], n: &str) -> f64 {
    s.iter().find(|p| p.parameter == n).unwrap().value
}

fn mk_parts(ids: [EntityId; 2], [mh, hh]: [fs_blake3::ContentHash; 2]) -> Vec<DraftPartDecl> {
    use {DraftGeometrySource as Gs, PrimitiveShape as Ps};
    let mr = |card, qp: &[(&str, [i8; 6], f64)]| {
        let query_point = qp.iter().map(|&(k, d, v)| kv(k, d, v)).collect();
        vec![MaterialRef {
            card,
            query_point,
            frame: FrameId(0),
            layer: None,
            unknown: vec![],
        }]
    };
    let (vl, om) = ([1, 0, -1, 0, 0, 0], [0, 0, -1, 0, 0, 0]);
    let eb = FidelityErrorBudget {
        relative_tol: 1e-6,
        absolute_tol: q(1e-12, [0; 6]),
        max_work_units: Some(128),
    };
    let init = DraftInitialState {
        position: [q(0.0, DL); 3],
        orientation: [1.0, 0.0, 0.0, 0.0],
        linear_velocity: [q(-2.0, vl), q(0.0, vl), q(0.0, vl)],
        angular_velocity: [q(0.0, om); 3],
        temperature: None,
    };
    let ps = DraftPartDecl {
        name: "stick".into(),
        entity: ids[0],
        kind: DraftPartKind::RigidBody,
        geometry: Gs::Primitive(Ps::Sphere {
            radius: q(P[4], DL),
        }),
        materials: mr(hh, &[("T", DT, P[3]), ("moisture_content", [0; 6], 0.12)]),
        physics: BTreeSet::from([PhysicsDomain::Mechanics]),
        fidelity: FidelityPolicy {
            rungs: vec![FidelityRung::Lumped],
            error_budget: eb,
        },
        uncertainty: vec![],
        initial: Some(init),
    };
    let mut pb = ps.clone();
    (pb.name, pb.entity, pb.initial, pb.kind) =
        ("batter".into(), ids[1], None, DraftPartKind::Shell);
    pb.materials = mr(mh, &[("T", DT, P[3])]);
    pb.fidelity.rungs = vec![FidelityRung::Modal { mode_count: 3 }];
    let disk = TensionedDiskGenerator {
        radius: q(P[0], DL),
        thickness: q(P[1], DL),
        tension: q(P[2], [0, 1, -2, 0, 0, 0]),
        radial_intervals: 6,
        azimuths: 24,
    };
    pb.geometry = Gs::TensionedDisk(disk);
    vec![ps, pb]
}

fn build(cards: (&MaterialCard, &MaterialCard), xy: [f64; 2]) -> (DraftAssemblyInput, CardStore) {
    let (mut store, mut cat) = (CardStore::new(), EntityCatalog::new());
    let hashes = [cards.0, cards.1].map(|c| store.insert_material(c.clone()));
    let root = cat.declare(Ed::assembly("ds")).unwrap();
    let [st, hd] = ["stick", "batter"].map(|n| cat.declare(Ed::part(root, n)).unwrap());
    let feature_a = cat.declare(Ed::surface(st, "tip")).unwrap();
    let feature_b = cat.declare(Ed::surface(hd, "hit")).unwrap();
    let entity = cat
        .declare(Ed::interface(
            root,
            "hit",
            InterfacePair::unordered(feature_a, feature_b),
        ))
        .unwrap();
    let state = vec![
        kv("maximum_force", DF, 1e5),
        kv("maximum_penetration", DL, 0.01),
        kv("gap", DL, 2e-4),
        kv("station_x", DL, xy[0]),
        kv("station_y", DL, xy[1]),
    ];
    let kind = DraftInterfaceKind::HertzContact { state };
    let hit = DraftInterfaceDecl {
        name: "hit".into(),
        entity,
        feature_a,
        feature_b,
        kind,
    };
    let mut env = DraftEnvironment::earth_lab_air();
    env.base.ambient_temperature = q(P[3], DT);
    let (tm, dt) = ([0, 0, 1, 0, 0, 0], q(1.0 / 48_000.0, [0, 0, 1, 0, 0, 0]));
    let sr = q(48_000.0, [0, 0, -1, 0, 0, 0]);
    let time_plan = TimePlan {
        horizon: q(0.01, tm),
        step_policy: StepPolicy::Fixed { dt },
        default_sample_rate: sr,
        scheme: TimeScheme::Explicit,
    };
    let inp = DraftAssemblyInput {
        name: "ds".into(),
        root,
        catalog: cat,
        reference_state: Some(ReferenceState::Assembled),
        environment: env,
        time_plan,
        parts: mk_parts([st, hd], hashes),
        joints: vec![],
        gear_meshes: vec![],
        interfaces: vec![hit],
        sources: vec![],
        controllers: vec![],
        probes: vec![],
    };
    (inp, store)
}

fn lower(inp: DraftAssemblyInput, store: &CardStore) -> Result<LoweredSystem, String> {
    let sc = Scenario::new("ds", 1, inp.environment.base.clone());
    let mut reg = LoweringRegistry::with_kernel_defaults();
    register_owner_defaults(&mut reg);
    let d = DraftAssembly::admit(inp, &sc).map_err(|e| e.to_string())?;
    reg.lower(&d, store).map_err(|e| e.to_string())
}

fn ref_sys(low: &LoweredSystem, xy: [f64; 2]) -> (System, [NodeId; 2], Vec<f64>) {
    let bps = &low.receipt.components[1].parameter_sources;
    let sps = &low.receipt.components[0].parameter_sources;
    let bdg = ProfileBudget {
        max_nodes: 10_000,
        max_triangles: 20_000,
        max_feature_evaluations: 100_000,
    };
    let spec = TensionedDiskSpec {
        radius_m: P[0],
        thickness_m: P[1],
        young_pa: pv(bps, "young_modulus"),
        poisson: pv(bps, "poisson_ratio"),
        density_kg_m3: pv(bps, "density"),
        tension_n_m: P[2],
        radial_intervals: 6,
        azimuths: 24,
    };
    let film = TensionedDisk::new(spec, bdg).unwrap();
    let rep = slice_window(
        &film.model.k,
        &film.model.m,
        (0.0, 1e7),
        &SliceOptions::default(),
    )
    .unwrap();
    let ms = &rep.modes[..3];
    let sm = |x, y| {
        ms.iter()
            .map(|m| film.sample_mode_at(&m.phi, x, y).unwrap())
            .collect()
    };
    let (phi_c, weights): (Vec<_>, Vec<_>) = (sm(0.0, 0.0), sm(xy[0], xy[1]));
    let z = |k| pv(bps, &format!("modal_damping_ratio[{k}]"));
    let pairs: Vec<_> = (0..3).map(|k| (ms[k].lambda.sqrt(), z(k))).collect();
    let ips = &low.receipt.interfaces[0].parameter_sources;
    let e_star = pv(ips, "effective_modulus");
    let mut b = SystemBuilder::new();
    let m = Mass::new(pv(sps, "mass")).unwrap();
    let m0 = Motion {
        position_m: 0.0,
        velocity_m_s: -2.0,
    };
    let sn = b.add_with_tangent_and_stage_residual("s", m, m0).unwrap();
    let id = StableId::new("k").unwrap();
    let att = ModalAttachment {
        id,
        kind: PortKind::MechanicalForceVelocity,
        width: NonZeroUsize::MIN,
        weights,
    };
    let mb = ModalAcousticTimeBudget::audible_reference();
    let st = ModalStructure::from_modes(48_000, &pairs, mb, vec![att]).unwrap();
    let hn = b.add_with_tangent("h", st.clone(), st.at_rest()).unwrap();
    let kh = (4.0 / 3.0) * e_star * P[4].sqrt();
    let law = Obstacle::new(vec![-1.0], 1, 1, vec![2e-4], vec![1.0], kh, 1.5, "r".into()).unwrap();
    let slv = ContactSolve {
        max_iterations: 64,
        maximum_force_n: 1e5,
        maximum_penetration_m: 0.01,
        force_absolute_tolerance_n: 1e-12,
        force_relative_tolerance: 1e-12,
    };
    let ci = ContactInterface::new(law, slv).unwrap();
    let ci = ci.publishing_reaction().unwrap();
    let cn = b.add_on_line("c", ci, ContactState::at(0.0)).unwrap();
    let p = |b: &SystemBuilder, n, k| b.port(n, k).unwrap();
    b.common_flow(p(&b, sn, "p"), &[p(&b, cn, "a")]).unwrap();
    b.common_flow(p(&b, hn, "k"), &[p(&b, cn, "b")]).unwrap();
    let (sk, rt) = (
        Some(heat_sink("a", P[3]).unwrap()),
        RouteTemperature::Fixed(P[3]),
    );
    b.heat_route(p(&b, hn, "heat"), &[], sk, rt).unwrap();
    (b.seal().unwrap(), [hn, cn], phi_c)
}

#[test]
fn modal_strike_parity_signs_and_shared_feature_refusals() {
    let (m, h) = (
        Pack::from_bytes(MYLAR).unwrap(),
        Pack::from_bytes(HICKORY).unwrap(),
    );
    let (cards, dt, xy) = ((m.card(), h.card()), 1.0 / 48_000.0, [0.06, 0.02]);
    let (inp, store) = build(cards, xy);
    let mut low = lower(inp, &store).unwrap();
    low.receipt.verify_completeness().unwrap();
    let (mut rsys, [hn, cn], phi_c) = ref_sys(&low, xy);
    let [ls, lh, lc] = ["stick", "batter", "interface.hit"].map(|k| low.node_map[k]);
    let (mut r_max, mut w_min, mut max_dr, mut max_dw) = (0.0_f64, 0.0_f64, 0.0_f64, 0.0_f64);
    let wc = |s: &System, n| {
        mq(s, n)
            .iter()
            .zip(&phi_c)
            .map(|(x, &c)| c * x.displacement_m_sqrt_kg)
            .sum::<f64>()
    };
    let (g, bdg, mode) = (
        CancelGate::new_clock_free(),
        Budget::INFINITE,
        ExecMode::Deterministic,
    );
    let key = StreamKey {
        seed: 7,
        kernel_id: 4,
        tile: 0,
        iteration: 0,
    };
    ArenaPool::new(ArenaConfig::default()).scope(|a| {
        let cx = Cx::new(&g, a, key, bdg, mode);
        for _ in 0..160 {
            low.system.step(dt, &cx).unwrap();
            rsys.step(dt, &cx).unwrap();
            let rl = low.system.state::<ContactState>(lc).unwrap().reaction_n;
            let rr = rsys.state::<ContactState>(cn).unwrap().reaction_n;
            assert!(rl >= 0.0);
            let (wl, wr) = (wc(&low.system, lh), wc(&rsys, hn));
            (r_max, w_min) = (r_max.max(rl), w_min.min(wl));
            (max_dr, max_dw) = (max_dr.max((rl - rr).abs()), max_dw.max((wl - wr).abs()));
            let (ql, qr) = (mq(&low.system, lh), mq(&rsys, hn));
            let dn = |q: &[ModalAcousticState]| {
                q[1].displacement_m_sqrt_kg
                    .hypot(q[2].displacement_m_sqrt_kg)
            };
            assert!((dn(ql) - dn(qr)).abs() < 1e-11);
        }
    });
    let v_end = low.system.state::<Motion>(ls).unwrap().velocity_m_s;
    eprintln!(
        "STRIKE Rmax={r_max:.2}N vend={v_end:.4}m/s wmin={w_min:.3e}m dR={max_dr:.2e} dW={max_dw:.2e}"
    );
    assert!(r_max > 1.0 && v_end > 0.0 && w_min < -1e-6 && max_dr < 1e-9 && max_dw < 1e-11);
    assert!(low.system.ledger().energy_residual_j().abs() < 1e-11);
    for (piston, want) in [(true, "combine piston-face"), (false, "two hertz-contact")] {
        let (mut d, s) = build(cards, xy);
        let mut i = d.interfaces[0].clone();
        let pr = InterfacePair::unordered(i.feature_a, i.feature_b);
        let e = d.catalog.declare(Ed::interface(d.root, "i2", pr)).unwrap();
        (i.name, i.entity) = ("i2".into(), e);
        if piston {
            i.kind = DraftInterfaceKind::PistonFace { state: vec![] };
        }
        d.interfaces.push(i);
        let got = lower(d, &s).err().unwrap();
        assert!(got.contains(want), "got `{got}`, want `{want}`");
    }
}
