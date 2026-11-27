//! Independent physics oracles, degenerate-pair gap sweep, and 0-ULP wiring parity
//! for `.11.9.2` (`owners::plate`, `PartLoweringKey(Shell, {Mechanics}, Modal)`).

use core::f64::consts::TAU;
use std::collections::BTreeSet;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{ClaimSet, MaterialCard, NormalizedMaterialCardPack, QueryPoint, SelectionPolicy};
use fs_modal::{ModePair, SliceOptions, slice_window, time::*};
use fs_plate::{
    AssemblyOptions, EdgeSupport, assemble,
    shell::{head::*, profile::ProfileBudget},
};
use fs_qty::semantic::{
    FrequencyConvention::Cyclic, QuantityKind, SemanticType, ValueForm::Static,
};
use fs_qty::{Dims, QtyAny, QuantitySpec};
use fs_scenario::assembly::{MaterialRef, ParameterKey, ReferenceState, draft::*};
use fs_scenario::{EntityCatalog, EntityDeclaration, FrameId, Scenario};
use fs_system::lowering::owners::{plate::DEGENERATE_GAP_TOL, register_owner_defaults};
use fs_system::{CardStore, LoweredSystem, LoweringError, LoweringRegistry, System};

const MYLAR: &[u8] =
    include_bytes!("../../fs-vkplate/packs/head-film-pet-mylar-dupont-matsuzawa1959.fsmatpk");
const J01: f64 = 2.404_825_557_695_773;
const R_M: f64 = 0.1778;
const H_M: f64 = 0.00025;
const T_N_M: f64 = 2500.0;
const RATE_HZ: u32 = 48_000;

fn query_prop(card: &MaterialCard, prop: &str, freq_hz: Option<f64>) -> f64 {
    let qs = |k| QuantitySpec::semantic(SemanticType::new(k, Static));
    let mut pt = QueryPoint::new();
    pt = pt
        .with_quantity("T", qs(QuantityKind::AbsoluteTemperature), 293.15)
        .unwrap();
    if let Some(f) = freq_hz {
        pt = pt
            .with_quantity("frequency", qs(QuantityKind::Frequency(Cyclic)), f)
            .unwrap();
    }
    let pol = SelectionPolicy::SingleClaimOnly;
    card.claims()
        .query(prop, &pt, pol)
        .unwrap()
        .evidence
        .value
        .value
}

fn mylar_props() -> (MaterialCard, [f64; 5]) {
    let pack = NormalizedMaterialCardPack::from_bytes(MYLAR).unwrap();
    let card = pack.card().clone();
    let rho = query_prop(&card, "density", None);
    let e = query_prop(&card, "young_modulus", None);
    let nu = query_prop(&card, "poisson_ratio", None);
    let (d, k0) = (e * H_M.powi(3) / (12.0 * (1.0 - nu * nu)), J01 / R_M);
    let omega = ((T_N_M * k0 * k0 + d * k0.powi(4)) / (rho * H_M)).sqrt();
    let frac = (d * k0 * k0) / (T_N_M + d * k0 * k0);
    (card, [rho, e, nu, omega, frac])
}

fn lower_head(
    card: &MaterialCard,
    n: u32,
    ri: u32,
    az: u32,
    tap: bool,
) -> Result<LoweredSystem, LoweringError> {
    let (mut store, mut cat) = (CardStore::new(), EntityCatalog::new());
    let hash = store.insert_material(card.clone());
    let root = cat.declare(EntityDeclaration::assembly("d")).unwrap();
    let head = cat.declare(EntityDeclaration::part(root, "head")).unwrap();
    let q = QtyAny::new;
    let (len, tim) = (Dims([1, 0, 0, 0, 0, 0]), Dims([0, 0, 1, 0, 0, 0]));
    let (frq, tmp) = (Dims([0, 0, -1, 0, 0, 0]), Dims([0, 0, 0, 1, 0, 0]));
    let target = SourceTarget::Port {
        entity: head,
        port: "mechanics".into(),
    };
    let history = TimeHistory::Constant(q(1.0, Dims([1, 1, -2, 0, 0, 0])));
    let sources = if tap {
        vec![SourceDecl {
            name: "tap".into(),
            target,
            history,
        }]
    } else {
        vec![]
    };
    let env = DraftEnvironment::earth_lab_air();
    let disk = TensionedDiskGenerator {
        radius: q(R_M, len),
        thickness: q(H_M, len),
        tension: q(T_N_M, Dims([0, 1, -2, 0, 0, 0])),
        radial_intervals: ri,
        azimuths: az,
    };
    let qp = vec![(ParameterKey::new("T", tmp), q(293.15, tmp))];
    let mat = MaterialRef {
        card: hash,
        query_point: qp,
        frame: FrameId(0),
        layer: None,
        unknown: vec![],
    };
    let eb = FidelityErrorBudget {
        relative_tol: 1.0e-4,
        absolute_tol: q(1.0e-9, len),
        max_work_units: None,
    };
    let fid = FidelityPolicy {
        rungs: vec![FidelityRung::Modal { mode_count: n }],
        error_budget: eb,
    };
    let part = DraftPartDecl {
        name: "head".into(),
        entity: head,
        kind: DraftPartKind::Shell,
        geometry: DraftGeometrySource::TensionedDisk(disk),
        materials: vec![mat],
        physics: BTreeSet::from([PhysicsDomain::Mechanics]),
        fidelity: fid,
        uncertainty: vec![],
        initial: None,
    };
    let step_policy = StepPolicy::Fixed {
        dt: q(1.0 / f64::from(RATE_HZ), tim),
    };
    let tp = TimePlan {
        horizon: q(0.01, tim),
        step_policy,
        default_sample_rate: q(f64::from(RATE_HZ), frq),
        scheme: TimeScheme::Explicit,
    };
    let inp = DraftAssemblyInput {
        name: "d".into(),
        root,
        catalog: cat,
        reference_state: Some(ReferenceState::Assembled),
        environment: env.clone(),
        time_plan: tp,
        parts: vec![part],
        joints: vec![],
        gear_meshes: vec![],
        interfaces: vec![],
        sources,
        controllers: vec![],
        probes: vec![],
    };
    let draft = DraftAssembly::admit(inp, &Scenario::new("d", 1, env.base)).unwrap();
    let mut reg = LoweringRegistry::with_kernel_defaults();
    register_owner_defaults(&mut reg);
    reg.lower(&draft, &store)
}

fn step_once(sys: &mut System, pool: &ArenaPool, gate: &CancelGate, dt: f64) {
    let key = StreamKey {
        seed: 1,
        kernel_id: 1,
        tile: 0,
        iteration: 0,
    };
    pool.scope(|a| {
        let cx = Cx::new(gate, a, key, Budget::INFINITE, ExecMode::Deterministic);
        sys.step(dt, &cx).unwrap();
    });
}

fn reference_model(
    card: &MaterialCard,
    disk: &TensionedDisk,
    modes: &[ModePair],
) -> (ModalAcousticTimeModel, Vec<f64>) {
    let opts = AssemblyOptions {
        pretension: 0.0,
        support: EdgeSupport::SimplySupported,
    };
    let bend = assemble(&disk.mesh, &disk.section, disk.outer_nodes(), &[], &opts).unwrap();
    let mut prod = vec![0.0; bend.free];
    let nodes = &disk.mesh.nodes;
    let ctr = nodes
        .iter()
        .position(|&(x, y)| x == 0.0 && y == 0.0)
        .and_then(|i| disk.model.dof_map[3 * i])
        .unwrap();
    let mut ref_modes = Vec::with_capacity(modes.len());
    for m in modes {
        let w = m.lambda.sqrt();
        bend.k.spmv(&m.phi, &mut prod);
        let sh = m.phi.iter().zip(&prod).map(|(a, b)| a * b).sum::<f64>() / m.lambda;
        let z = 0.5 * query_prop(card, "loss_tangent_tensile", Some(w / TAU)) * sh;
        let p = fs_math::c64::C64::ZERO;
        ref_modes.push(ModalAcousticMode {
            angular_frequency_rad_s: w,
            damping_ratio: z,
            pressure_per_modal_velocity: p,
        });
    }
    let bdg = ModalAcousticTimeBudget::audible_reference();
    let rmod = ModalAcousticTimeModel::try_new(RATE_HZ, ref_modes, bdg).unwrap();
    (rmod, modes.iter().map(|m| m.phi[ctr]).collect())
}

/// Measure fundamental `(omega, zeta, receipt_share)` from a stepped 1-mode `System`
/// driven by a held 1 N center tap over one full cycle (263 steps at 48 kHz, f0 ~ 182.6 Hz).
fn measure_fundamental(card: &MaterialCard, ri: u32, az: u32) -> (f64, f64, f64) {
    let mut low = lower_head(card, 1, ri, az, true).unwrap();
    low.receipt.verify_completeness().unwrap();
    let ps = &low.receipt.components[0].parameter_sources;
    let rec = ps
        .iter()
        .find(|p| p.parameter == "bending_share[0]")
        .unwrap()
        .value;
    let (node, dt) = (low.node_map["head"], 1.0 / f64::from(RATE_HZ));
    let (pool, gate) = (
        ArenaPool::new(ArenaConfig::default()),
        CancelGate::new_clock_free(),
    );
    let (mut v2, mut v1, mut iv2, mut omega) = (0.0, 0.0, 0.0, 0.0);
    for step in 0..263usize {
        step_once(&mut low.system, &pool, &gate, dt);
        let st = low.system.state::<Vec<ModalAcousticState>>(node).unwrap();
        let v = st[0].velocity_m_sqrt_kg_per_s;
        iv2 += (v1 * v1 + v1 * v + v * v) * (dt / 3.0);
        if step == 66 {
            omega = ((v2 + v) / (2.0 * v1)).acos() / dt;
        }
        (v2, v1) = (v1, v);
    }
    assert!(low.system.ledger().energy_residual_j().abs() < 1.0e-12);
    (
        omega,
        low.system.ledger().sink_heat_j() / (2.0 * omega * iv2),
        rec,
    )
}

#[test]
fn modal_plate_doublet_refusal_and_zero_ulp_parity() {
    let (card, [rho, e_pa, nu, omega_a, _]) = mylar_props();
    let bdg = ProfileBudget {
        max_nodes: 10_000,
        max_triangles: 20_000,
        max_feature_evaluations: 100_000,
    };
    let spec = TensionedDiskSpec {
        radius_m: R_M,
        thickness_m: H_M,
        young_pa: e_pa,
        poisson: nu,
        density_kg_m3: rho,
        tension_n_m: T_N_M,
        radial_intervals: 6,
        azimuths: 24,
    };
    let disk = TensionedDisk::new(spec, bdg).unwrap();
    let opts = SliceOptions::default();
    let sl = |b| {
        slice_window(&disk.model.k, &disk.model.m, (0.0, b), &opts)
            .unwrap()
            .modes
    };
    let modes_7 = sl(20.0 * omega_a * omega_a);
    assert!(modes_7.len() >= 7);
    for i in 1..7 {
        let gap = (modes_7[i].lambda - modes_7[i - 1].lambda) / modes_7[i - 1].lambda;
        eprintln!("mode {i}->{}: rel_gap={gap:.3e}", i + 1);
    }
    assert!((modes_7[2].lambda - modes_7[1].lambda) / modes_7[1].lambda <= DEGENERATE_GAP_TOL);
    let err2 = format!("{}", lower_head(&card, 2, 6, 24, false).err().unwrap());
    assert!(err2.contains("splits a degenerate pair") && err2.contains("choose 1 or 3"));
    let mut low = lower_head(&card, 3, 6, 24, true).unwrap();
    low.receipt.verify_completeness().unwrap();
    let (mut rmod, cphi) = reference_model(&card, &disk, &sl(8.0 * omega_a * omega_a)[..3]);
    let pool = ArenaPool::new(ArenaConfig::default());
    let gate = CancelGate::new_clock_free();
    let (dt, head) = (1.0 / f64::from(RATE_HZ), low.node_map["head"]);
    let bits = |s: &ModalAcousticState| {
        (
            s.displacement_m_sqrt_kg.to_bits(),
            s.velocity_m_sqrt_kg_per_s.to_bits(),
        )
    };
    for step in 0..20 {
        step_once(&mut low.system, &pool, &gate, dt);
        rmod.step(&cphi).unwrap();
        let got = low.system.state::<Vec<ModalAcousticState>>(head).unwrap();
        for (k, (g, r)) in got.iter().zip(rmod.states()).enumerate() {
            assert_eq!(bits(g), bits(r), "m{k}@{step}");
        }
    }
    assert!(low.system.ledger().energy_residual_j().abs() < 1.0e-12);
}

#[test]
fn modal_plate_independent_oracles_and_missing_loss() {
    let (card, [_, _, nu, omega_a, frac_a]) = mylar_props();
    let (w_c, zeta_c, share_c) = measure_fundamental(&card, 6, 24);
    let (w_f, zeta_f, share_f) = measure_fundamental(&card, 12, 48);
    let (err_c, err_f) = (
        (w_c - omega_a).abs() / omega_a,
        (w_f - omega_a).abs() / omega_a,
    );
    let eta = |w| query_prop(&card, "loss_tangent_tensile", Some(w / TAU));
    let (frac_c, frac_f) = (2.0 * zeta_c / eta(w_c), 2.0 * zeta_f / eta(w_f));
    // Simply-supported disk Gaussian curvature scales membrane-limit D k^2/(T + D k^2) by
    // 1 - 2(1 - nu)/J01^2; 2% covers the ~sqrt(D/T)/R ~ 0.97% rim boundary layer + mesh error.
    let frac_kl = frac_a * (1.0 - 2.0 * (1.0 - nu) / (J01 * J01));
    let (b_err_c, b_err_f) = (
        (frac_c - frac_kl).abs() / frac_kl,
        (frac_f - frac_kl).abs() / frac_kl,
    );
    let (fa_hz, fc_hz, ff_hz) = (omega_a / TAU, w_c / TAU, w_f / TAU);
    eprintln!(
        "ORACLE(a): analytic_hz={fa_hz:.4} coarse_hz={fc_hz:.4} (err={err_c:.4e}) fine_hz={ff_hz:.4} (err={err_f:.4e})"
    );
    eprintln!(
        "ORACLE(b): fine={frac_f:.6e} coarse={frac_c:.6e} frac_kl={frac_kl:.6e} membrane={frac_a:.6e} err_kl=({b_err_c:.4e}->{b_err_f:.4e}) rec=({share_c:.6e},{share_f:.6e})"
    );
    assert!(err_c < 0.005 && err_f < 0.0015 && err_f < err_c);
    assert!(b_err_f < 0.02 && b_err_f < b_err_c);
    assert!(
        (frac_f - share_f).abs() / share_f < 1.0e-3 && (frac_c - share_c).abs() / share_c < 1.0e-3
    );
    let (mut no_loss, cl) = (ClaimSet::new(), card.claims());
    for id in cl.observation_ids() {
        let obs = cl.observation(id).unwrap().clone();
        no_loss.register_observation(obs).unwrap();
    }
    for (_, c) in cl.claims_ordered() {
        if c.key.name() != "loss_tangent_tensile" {
            no_loss.insert_claim(c.clone()).unwrap();
        }
    }
    let stripped = MaterialCard::assemble(card.id().clone(), no_loss, vec![]).unwrap();
    let err = lower_head(&stripped, 3, 6, 24, false).err().unwrap();
    assert!(matches!(err, LoweringError::MatDb { .. }));
    assert!(err.to_string().contains("loss_tangent_tensile"));
}
