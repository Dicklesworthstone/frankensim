//! Two modal heads composed with a real N2 `GasCavity` via `PistonFace` (`.11.9.3`): static
//! balance, in-phase/breathing modes, two-`dt` `ExactPhs` convergence, gas NL, timing & refusals.

use core::f64::consts::{PI, TAU};
use core::num::NonZeroUsize;
use std::{collections::BTreeSet, time::Instant};

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::ContentHash;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_gas::component::ClosedGasParcel;
use fs_matdb::{MaterialCard, NormalizedMaterialCardPack as Pack};
use fs_modal::{SliceOptions, slice_window, time::ModalAcousticState};
use fs_plate::shell::{head::*, profile::ProfileBudget};
use fs_ports::{PortKind, StableId, component::*};
use fs_qty::{Dims, QtyAny};
use fs_scenario::assembly::{draft::*, *};
use fs_scenario::{EntityDeclaration as Ed, *};
use fs_system::{lowering::owners::register_owner_defaults, phs_exact::*, *};

const MYLAR: &[u8] =
    include_bytes!("../../fs-vkplate/packs/head-film-pet-mylar-dupont-matsuzawa1959.fsmatpk");
const N2: &[u8] =
    include_bytes!("../../../examples/assemblies/packs/nitrogen-n2-mcbride2002-lemmon2004.fsmcdpk");
const P: [f64; 6] = [0.1778, 0.00025, 2500.0, 0.3556, 293.15, 90_000.0]; // R, h, T, cav_h, T0, P0

fn with_cx<T>(f: impl FnOnce(&Cx<'_>) -> T) -> T {
    let g = CancelGate::new_clock_free();
    let k = StreamKey {
        seed: 1,
        kernel_id: 1,
        tile: 0,
        iteration: 0,
    };
    ArenaPool::new(ArenaConfig::default()).scope(|a| {
        let cx = Cx::new(&g, a, k, Budget::INFINITE, ExecMode::Deterministic);
        f(&cx)
    })
}

fn drum_parts(ids: [EntityId; 3], mh: ContentHash, gh: ContentHash) -> Vec<DraftPartDecl> {
    let ([h1, h2, gas], [r, h, t_n, ch, t0, p0], q, d) = (ids, P, QtyAny::new, Dims);
    let (len, tmp) = (d([1, 0, 0, 0, 0, 0]), d([0, 0, 0, 1, 0, 0]));
    let prs = d([-1, 1, -2, 0, 0, 0]);
    let mat = |card, k, v: QtyAny| MaterialRef {
        card,
        query_point: vec![(ParameterKey::new(k, v.dims), v)],
        frame: FrameId(0),
        layer: None,
        unknown: vec![],
    };
    let p1 = DraftPartDecl {
        name: "head1".into(),
        entity: h1,
        kind: DraftPartKind::Shell,
        geometry: DraftGeometrySource::TensionedDisk(TensionedDiskGenerator {
            radius: q(r, len),
            thickness: q(h, len),
            tension: q(t_n, Dims([0, 1, -2, 0, 0, 0])),
            radial_intervals: 6,
            azimuths: 24,
        }),
        materials: vec![mat(mh, "T", q(t0, tmp))],
        physics: BTreeSet::from([PhysicsDomain::Mechanics]),
        fidelity: FidelityPolicy {
            rungs: vec![FidelityRung::Modal { mode_count: 1 }],
            error_budget: FidelityErrorBudget {
                relative_tol: 1e-6,
                absolute_tol: QtyAny::dimensionless(1e-12),
                max_work_units: Some(128),
            },
        },
        uncertainty: vec![],
        initial: None,
    };
    let (mut p2, mut pc) = (p1.clone(), p1.clone());
    (p2.name, p2.entity) = ("head2".into(), h2);
    (pc.name, pc.entity, pc.kind) = ("cavity".into(), gas, DraftPartKind::Cavity);
    pc.geometry = DraftGeometrySource::Primitive(PrimitiveShape::Cylinder {
        radius: q(r, len),
        height: q(ch, len),
    });
    pc.materials = vec![mat(gh, "pressure", q(p0, prs))];
    pc.physics = BTreeSet::from([PhysicsDomain::Fluid]);
    pc.fidelity.rungs = vec![FidelityRung::Lumped];
    vec![p1, p2, pc]
}

fn lower_drum(
    cards: (&MaterialCard, &MaterialCard),
    rate: u32,
    (tap, ambient): (f64, bool),
    f1_st: Vec<(ParameterKey, QtyAny)>,
) -> Result<LoweredSystem, LoweringError> {
    let ([r, _, _, _, t0, p0], q) = (P, QtyAny::new);
    let (mut store, mut cat) = (CardStore::new(), EntityCatalog::new());
    let mh = store.insert_material(cards.0.clone());
    let gh = store.insert_material(cards.1.clone());
    let root = cat.declare(Ed::assembly("drum")).unwrap();
    let [h1, h2, gas] =
        ["head1", "head2", "cavity"].map(|n| cat.declare(Ed::part(root, n)).unwrap());
    let mut sf = |p, n| cat.declare(Ed::surface(p, n)).unwrap();
    let (h1i, h1o, h2i, h2o) = (sf(h1, "in"), sf(h1, "out"), sf(h2, "in"), sf(h2, "out"));
    let (gt, gb) = (sf(gas, "top"), sf(gas, "bot"));
    let u = InterfacePair::unordered;
    let f1 = cat.declare(Ed::interface(root, "f1", u(h1i, gt))).unwrap();
    let f2 = cat.declare(Ed::interface(root, "f2", u(h2i, gb))).unwrap();
    let (tim, tmp) = (Dims([0, 0, 1, 0, 0, 0]), Dims([0, 0, 0, 1, 0, 0]));
    let (prs, frc) = (Dims([-1, 1, -2, 0, 0, 0]), Dims([1, 1, -2, 0, 0, 0]));
    let (area, mut sources) = (PI * r * r, Vec::new());
    let mut add = |n: &str, entity, v| {
        let target = SourceTarget::Port {
            entity,
            port: "mechanics".into(),
        };
        sources.push(SourceDecl {
            name: n.into(),
            target,
            history: TimeHistory::Constant(q(v, frc)),
        });
    };
    if ambient {
        add("amb1", h1o, -p0 * area);
        add("amb2", h2o, p0 * area);
    }
    if tap != 0.0 {
        add("tap", h1, tap);
    }
    let s2 = vec![(ParameterKey::new("sign", Dims::NONE), q(-1.0, Dims::NONE))];
    let ifc = |n: &str, entity, a, b, state| DraftInterfaceDecl {
        name: n.into(),
        entity,
        feature_a: a,
        feature_b: b,
        kind: DraftInterfaceKind::PistonFace { state },
    };
    let mut env = DraftEnvironment::earth_lab_air();
    (env.base.ambient_temperature, env.base.ambient_pressure) = (q(t0, tmp), q(p0, prs));
    let dt = q(1.0 / f64::from(rate), tim);
    let tp = TimePlan {
        horizon: q(0.02, tim),
        step_policy: StepPolicy::Fixed { dt },
        default_sample_rate: q(f64::from(rate), Dims([0, 0, -1, 0, 0, 0])),
        scheme: TimeScheme::Explicit,
    };
    let inp = DraftAssemblyInput {
        name: "drum".into(),
        root,
        catalog: cat,
        reference_state: Some(ReferenceState::Assembled),
        environment: env.clone(),
        time_plan: tp,
        parts: drum_parts([h1, h2, gas], mh, gh),
        joints: vec![],
        gear_meshes: vec![],
        interfaces: vec![ifc("f1", f1, h1i, gt, f1_st), ifc("f2", f2, h2i, gb, s2)],
        sources,
        controllers: vec![],
        probes: vec![],
    };
    let mut reg = LoweringRegistry::with_kernel_defaults();
    register_owner_defaults(&mut reg);
    let draft = DraftAssembly::admit(inp, &Scenario::new("drum", 1, env.base)).unwrap();
    reg.lower(&draft, &store)
}

fn exact_phs(w0: f64, z0: f64, kva2: f64, phi_c: f64) -> ExactPortHamiltonianComponent {
    let (w2, c) = (w0 * w0, 2.0 * z0 * w0);
    let (mut h, mut j, mut r) = (vec![0.0; 16], vec![0.0; 16], vec![0.0; 16]);
    (h[0], h[2], h[8], h[10], h[5], h[15]) = (w2 + kva2, -kva2, -kva2, w2 + kva2, 1.0, 1.0);
    (j[1], j[4], j[11], j[14], r[5], r[15]) = (1.0, -1.0, 1.0, -1.0, c, c);
    let qs = Box::new(QuadraticStorage::new(h, 4).unwrap());
    let phs = PortHamiltonian::new(4, 1, j, r, vec![0.0, phi_c, 0.0, 0.0], qs).unwrap();
    let tid = StableId::new("tap").unwrap();
    let pk = PortKind::MechanicalForceVelocity;
    let p = PortSpec::power(tid, pk, PortRole::EffortIn, NonZeroUsize::MIN).unwrap();
    let bdg = LinearSplitBudget {
        max_modes: 8,
        max_omega_dt: 0.45 * PI,
    };
    ExactPortHamiltonianComponent::new(phs, vec![p], bdg).unwrap()
}

fn mode0(low: &LoweredSystem) -> [f64; 4] {
    let ps = &low.receipt.components[0].parameter_sources;
    let q = |nm| ps.iter().find(|p| p.parameter == nm).unwrap().value;
    let bdg = ProfileBudget {
        max_nodes: 10_000,
        max_triangles: 20_000,
        max_feature_evaluations: 100_000,
    };
    let spec = TensionedDiskSpec {
        radius_m: P[0],
        thickness_m: P[1],
        young_pa: q("young_modulus"),
        poisson: q("poisson_ratio"),
        density_kg_m3: q("density"),
        tension_n_m: P[2],
        radial_intervals: 6,
        azimuths: 24,
    };
    let film = TensionedDisk::new(spec, bdg).unwrap();
    let opts = SliceOptions::default();
    let rep = slice_window(&film.model.k, &film.model.m, (0.0, 1e7), &opts).unwrap();
    let m0 = &rep.modes[0];
    let nodes = &film.mesh.nodes;
    let i = nodes
        .iter()
        .position(|&(x, y)| x == 0.0 && y == 0.0)
        .unwrap();
    let a0 = film.modal_area(&m0.phi).unwrap();
    let ph = m0.phi[film.model.dof_map[3 * i].unwrap()];
    [m0.lambda.sqrt(), q("modal_damping_ratio[0]"), a0, ph]
}

fn run_pair(
    cards: (&MaterialCard, &MaterialCard),
    [w0, z0, a0, phi_c]: [f64; 4],
    (rate, steps, tap): (u32, usize, f64),
) -> (Vec<f64>, Vec<f64>, [f64; 7]) {
    let mut low = lower_drum(cards, rate, (tap, true), vec![]).unwrap();
    let m = &low.node_map;
    let (h1, h2, gc) = (m["head1"], m["head2"], m["cavity"]);
    let p0 = low.system.state::<ClosedGasParcel>(gc).unwrap();
    let gs = p0.state().state();
    let gam = gs.specific_heat_cp / (gs.specific_heat_cp - gs.specific_gas_constant);
    let (kb, v0, dt) = (gam * gs.pressure, p0.volume_m3(), 1.0 / f64::from(rate));
    let ex = exact_phs(w0, z0, (kb / v0) * a0 * a0, phi_c);
    let (mut xs, mut xn, mut out) = (vec![0.0; 4], vec![0.0; 4], [0.0; 2]);
    let (mut vp, mut vm) = (Vec::with_capacity(steps), Vec::with_capacity(steps));
    let (mut err, mut sc, mut dvr, mut tc, mut te) = (0.0_f64, 0.0_f64, 0.0_f64, 0.0, 0.0);
    with_cx(|cx| {
        for _ in 0..steps {
            let t = Instant::now();
            low.system.step(dt, cx).unwrap();
            tc += t.elapsed().as_secs_f64();
            let t = Instant::now();
            ex.propose(&xs, dt, &[tap, P[4]], &mut xn, &mut out, cx)
                .unwrap();
            te += t.elapsed().as_secs_f64();
            xs.copy_from_slice(&xn);
            let s1 = low.system.state::<Vec<ModalAcousticState>>(h1).unwrap()[0];
            let s2 = low.system.state::<Vec<ModalAcousticState>>(h2).unwrap()[0];
            vp.push(0.5 * (s1.velocity_m_sqrt_kg_per_s + s2.velocity_m_sqrt_kg_per_s));
            vm.push(0.5 * (s1.velocity_m_sqrt_kg_per_s - s2.velocity_m_sqrt_kg_per_s));
            let (q1, q2) = (s1.displacement_m_sqrt_kg, s2.displacement_m_sqrt_kg);
            err = err.max((q1 - xs[0]).abs()).max((q2 - xs[2]).abs());
            sc = sc.max(xs[0].abs()).max(xs[2].abs());
            let v = low.system.state::<ClosedGasParcel>(gc).unwrap().volume_m3();
            dvr = dvr.max((v - v0).abs() / v0);
        }
    });
    assert!(low.system.ledger().energy_residual_j().abs() < 1e-11);
    let inv = 1e6 / (steps as f64);
    (vp, vm, [err / sc, dvr, kb, v0, gam, tc * inv, te * inv])
}

#[test]
fn composed_two_heads_cavity_oracles_convergence_and_refusals() {
    let mylar = Pack::from_bytes(MYLAR).unwrap();
    let n2 = Pack::from_bytes(N2).unwrap();
    let cards = (mylar.card(), n2.card());
    let mut eq = lower_drum(cards, 48_000, (0.0, true), vec![]).unwrap();
    eq.receipt.verify_completeness().unwrap();
    let p0p = mode0(&eq);
    let mut no_amb = lower_drum(cards, 48_000, (0.0, false), vec![]).unwrap();
    with_cx(|cx| {
        for _ in 0..20 {
            eq.system.step(1.0 / 48_000.0, cx).unwrap();
        }
        no_amb.system.step(1.0 / 48_000.0, cx).unwrap();
    });
    let q_of = |l: &LoweredSystem, k| {
        let id = l.node_map[k];
        l.system.state::<Vec<ModalAcousticState>>(id).unwrap()[0].displacement_m_sqrt_kg
    };
    assert!(q_of(&eq, "head1").abs() < 1e-13 && q_of(&eq, "head2").abs() < 1e-13);
    assert!(q_of(&no_amb, "head1").abs() > 1e-6);
    let (vp, vm, [ec, dvc, kb, v0, gam, uc, xc]) = run_pair(cards, p0p, (48_000, 240, 0.05));
    let (_, _, [ef, dvf, _, _, _, uf, xf]) = run_pair(cards, p0p, (192_000, 960, 0.05));
    let [w0, _, a0, _] = p0p;
    let w_in = ((vp[64] + vp[66]) / (2.0 * vp[65])).acos() * 48_000.0;
    let w_br = ((vm[44] + vm[46]) / (2.0 * vm[45])).acos() * 48_000.0;
    let w_br_r = (w0 * w0 + 2.0 * (kb / v0) * a0 * a0).sqrt();
    let (fi, fi_r, fb, fb_r) = (w_in / TAU, w0 / TAU, w_br / TAU, w_br_r / TAU);
    let (ri, rb) = ((w_in - w0).abs() / w0, (w_br - w_br_r).abs() / w_br_r);
    let nl = 0.5 * (gam + 1.0) * dvf;
    eprintln!(
        "ORACLE in={fi:.4}/{fi_r:.4}Hz({ri:.2e}) br={fb:.4}/{fb_r:.4}Hz({rb:.2e}) K={kb:.1} V={v0:.5} a0={a0:.5} CONV 48k={ec:.3e} 192k={ef:.3e} dv/V={dvc:.3e} nl={nl:.3e} us(comp={uc:.1}/{uf:.1},phs={xc:.1}/{xf:.1})"
    );
    assert!(ri < 1e-4 && rb < 1e-3 && ef < ec && ef < 1e-4 + 10.0 * nl);
    let bad = |n, d, v| {
        let st = vec![(ParameterKey::new(n, d), QtyAny::new(v, d))];
        let e = lower_drum(cards, 48_000, (0.0, true), st).err().unwrap();
        format!("{e}")
    };
    let ea = bad("area", Dims([2, 0, 0, 0, 0, 0]), 0.05);
    let es = bad("sign", Dims::NONE, 2.0);
    assert!(ea.contains("area") && ea.contains("modal shell") && es.contains("sign"));
}
