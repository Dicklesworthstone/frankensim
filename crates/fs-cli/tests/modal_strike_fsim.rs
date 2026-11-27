//! End-to-end `.fsim` drum strike (`.11.9.4`): 2-head + N2 cavity + hickory stick composition,
//! center vs off-center `(1,1)` nodal selectivity, and station/overhang/rim/non-disk refusals.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_dcontact::component::ContactState;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_fsim::{PackCatalog, admit_assembly_source};
use fs_gas::component::ClosedGasParcel;
use fs_nlmodal::component::ModalStructure;
use fs_ports::component::Component;
use fs_system::{
    LoweredSystem, LoweringRegistry, elements::Motion, lowering::owners::register_owner_defaults,
};

const MYLAR: &[u8] =
    include_bytes!("../../fs-vkplate/packs/head-film-pet-mylar-dupont-matsuzawa1959.fsmatpk");
const HICKORY: &[u8] =
    include_bytes!("../../fs-vkplate/packs/stick-hickory-shagbark-fpl-gtr282.fsmatpk");
const N2: &[u8] =
    include_bytes!("../../../examples/assemblies/packs/nitrogen-n2-mcbride2002-lemmon2004.fsmcdpk");
type MState = <ModalStructure as Component>::State;

fn lower_fsim(
    packs: &PackCatalog,
    xy: [f64; 2],
    extra: &str,
    disk_hd: bool,
) -> Result<LoweredSystem, String> {
    let hd_g = if disk_hd {
        ":kind shell :geometry (tensioned-disk :radius 0.1778m :thickness 0.00025m :tension 2500N/m :radial-intervals 6 :azimuths 24) :fidelity (modal 3)"
    } else {
        ":kind rigid-body :geometry (sphere 0.015m) :fidelity lumped"
    };
    let src = format!(
        "(assembly \"drum\" (version 1) (units \"SI\") (seed 0x7) \
         (versions (constellation \"2026-07\") (workspace \"t\")) \
         (budget (wall 60s) (mem 64MiB)) \
         (capability :cores 1cores :mem 64MiB :wall 60s :ops (assembly.*)) \
         (reference-state assembled) (environment :ambient-temp 293.15K :ambient-pressure 90000Pa) \
         (time-plan :horizon 0.01s :dt 0.000020833333333333332s :sample-rate 48000Hz) \
         (part \"stick\" :kind rigid-body :geometry (sphere 0.015m) :material \"carya-ovata\" \
           :query-at (T 293.15K moisture_content 0.12) :physics (mechanics) :fidelity lumped \
           :initial (:position (0m 0m 0m) :orientation (1 0 0 0) \
             :linear-velocity (-2m/s 0m/s 0m/s) :angular-velocity (0rad/s 0rad/s 0rad/s))) \
         (part \"batter\" {hd_g} :material \"polyethylene-terephthalate\" :query-temp 293.15K :physics (mechanics)) \
         (interface \"hit\" :kind hertz-contact :a \"stick\" :b \"batter\" \
           :maximum-force 100000N :maximum-penetration 0.01m :gap 0.0002m \
           :station-x {}m :station-y {}m) {extra})",
        xy[0], xy[1]
    );
    let adm = admit_assembly_source(src.as_bytes(), packs).map_err(|e| e.what)?;
    let mut reg = LoweringRegistry::with_kernel_defaults();
    register_owner_defaults(&mut reg);
    reg.lower(&adm.draft, &adm.card_store)
        .map_err(|e| e.to_string())
}

fn step_peaks(packs: &PackCatalog, xy: [f64; 2]) -> ([f64; 3], [f64; 2]) {
    let mut low = lower_fsim(packs, xy, "", true).unwrap();
    let (hid, dt, mut pk, mut best) = (
        low.node_map["batter"],
        1.0 / 48_000.0,
        [0.0_f64; 3],
        (0.0_f64, [0.0; 2]),
    );
    let (g, k) = (
        CancelGate::new_clock_free(),
        StreamKey {
            seed: 7,
            kernel_id: 4,
            tile: 0,
            iteration: 0,
        },
    );
    ArenaPool::new(ArenaConfig::default()).scope(|a| {
        let cx = Cx::new(&g, a, k, Budget::INFINITE, ExecMode::Deterministic);
        for _ in 0..160 {
            low.system.step(dt, &cx).unwrap();
            let qv = low.system.state::<MState>(hid).unwrap();
            for i in 0..3 {
                pk[i] = pk[i].max(qv[i].displacement_m_sqrt_kg.abs());
            }
            let v = [qv[1].displacement_m_sqrt_kg, qv[2].displacement_m_sqrt_kg];
            if v[0] * v[0] + v[1] * v[1] > best.0 {
                best = (v[0] * v[0] + v[1] * v[1], v);
            }
        }
    });
    (pk, best.1)
}

#[test]
fn fsim_drum_strike_selectivity_and_refusals() {
    let packs = PackCatalog::from_raw_bytes(&[MYLAR, HICKORY, N2]).unwrap();
    let drum_extra = "(part \"resonant\" :kind shell :geometry (tensioned-disk :radius 0.1778m \
        :thickness 0.00025m :tension 2500N/m :radial-intervals 6 :azimuths 24) \
        :material \"polyethylene-terephthalate\" :query-temp 293.15K :physics (mechanics) :fidelity (modal 3)) \
      (part \"cavity\" :kind cavity :geometry (cylinder 0.1778m 0.3556m) :material \"N2\" \
        :query-at (pressure 90000Pa) :physics (fluid) :fidelity lumped) \
      (interface \"f1\" :kind piston-face :a \"batter\" :b \"cavity\") \
      (interface \"f2\" :kind piston-face :a \"resonant\" :b \"cavity\" :sign -1) \
      (source \"amb1\" :part \"batter\" :port \"face\" :history (constant -8938.234549111948N)) \
      (source \"amb2\" :part \"resonant\" :port \"face\" :history (constant 8938.234549111948N))";
    let mut drum = lower_fsim(&packs, [0.06, 0.02], drum_extra, true).unwrap();
    drum.receipt.verify_completeness().unwrap();
    let [ds, dr, dg, dc] =
        ["stick", "resonant", "cavity", "interface.hit"].map(|k| drum.node_map[k]);
    let (dt, mut p_max, mut rq_max, mut r_max) = (1.0 / 48_000.0, 90_000.0_f64, 0.0_f64, 0.0_f64);
    let (g, k) = (
        CancelGate::new_clock_free(),
        StreamKey {
            seed: 7,
            kernel_id: 4,
            tile: 0,
            iteration: 0,
        },
    );
    ArenaPool::new(ArenaConfig::default()).scope(|a| {
        let cx = Cx::new(&g, a, k, Budget::INFINITE, ExecMode::Deterministic);
        for _ in 0..160 {
            drum.system.step(dt, &cx).unwrap();
            let rx = drum.system.state::<ContactState>(dc).unwrap().reaction_n;
            assert!(rx >= 0.0);
            r_max = r_max.max(rx);
            p_max = p_max.max(
                drum.system
                    .state::<ClosedGasParcel>(dg)
                    .unwrap()
                    .state()
                    .state()
                    .pressure,
            );
            rq_max = rq_max.max(
                drum.system.state::<MState>(dr).unwrap()[0]
                    .displacement_m_sqrt_kg
                    .abs(),
            );
        }
    });
    let v_end = drum.system.state::<Motion>(ds).unwrap().velocity_m_s;
    let res = drum.system.ledger().energy_residual_j().abs();
    eprintln!(
        "FSIM DRUM Rmax={r_max:.2}N dp={:.2}Pa res={rq_max:.3e} vend={v_end:.4}m/s ledger={res:.2e}",
        p_max - 90_000.0
    );
    assert!(r_max > 1.0 && p_max > 90_001.0 && v_end > 0.0 && rq_max > 1e-7 && res < 1e-10);
    let ((pc, _), (px, vx), (py, vy)) = (
        step_peaks(&packs, [0.0, 0.0]),
        step_peaks(&packs, [0.0889, 0.0]),
        step_peaks(&packs, [0.0, 0.0889]),
    );
    assert!(pc[0] > 1e-6 && pc[1] < 1e-12 * pc[0] && pc[2] < 1e-12 * pc[0]);
    let (nx, ny) = (vx[0].hypot(vx[1]), vy[0].hypot(vy[1]));
    let cos_xy = (vx[0] * vy[0] + vx[1] * vy[1]).abs() / (nx * ny);
    assert!(nx > 0.1 * px[0] && ny > 0.1 * py[0] && cos_xy < 1e-6);
    eprintln!("FSIM SEL cos_xy={cos_xy:.2e} c1/c0={:.2e}", pc[1] / pc[0]);
    for (xy, disk, want) in [
        ([0.1878, 0.0], true, "outside"),
        ([0.1703, 0.0], true, "overhangs"),
        ([0.1778, 0.0], true, "zero modal participation"),
        ([0.06, 0.02], false, "TensionedDisk"),
    ] {
        let got = lower_fsim(&packs, xy, "", disk).err().unwrap();
        assert!(got.contains(want), "got `{got}`, want `{want}`");
    }
}
