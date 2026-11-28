//! Direct `engine-drum` equivalence oracle (`.11.9.5c-2`, closing `.11.9.2..5`):
//! real `.fsim` 2-head Mylar + N2 `GasCavity` + shagbark-hickory stick + 2 microphones compared
//! against `fs_couple::LinearImpactSystem` + `BaffledPressure::new_multirate` built identically to
//! `apps/instrument-lab/engine-drum/src/{params.rs:219-220, instrument.rs:17, engine.rs:218-303,424-454}`
//! at `1f982f3f8a` (zero dependency edge to `apps/instrument-lab`), plus convergence/scaling
//! attribution of the cavity model-form difference, planted physical-error teeth, and independent
//! physical sign oracles on both batter (`:sign +1`) and resonant (`:sign -1`) heads.

use std::f64::consts::PI;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_couple::modal_acoustic_time::{ModalAcousticState, ModalAcousticTimeBudget};
use fs_couple::pcm_wav::baffled::{BaffledPressure, RayleighMedium, SurfaceSample};
use fs_couple::render::plate::impact::linear::{
    LinearImpactConfig, LinearImpactSystem, VolumeConnection,
};
use fs_couple::render::plate::impact::{BodyPotential, ImpactBody, VolumeSpring};
use fs_couple::render::schedule::force::coupled::ModalCouplingConfig;
use fs_couple::render::schedule::force::coupled::contact::ModalContactConfig;
use fs_couple::render::schedule::force::coupled::contact::multiple::MultiContactConfig;
use fs_dcontact::Obstacle;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_fsim::{PackCatalog, admit_assembly_source};
use fs_gas::component::ClosedGasParcel;
use fs_material::gas::GasState;
use fs_system::acoustic::{AcousticMicrophone, RadiatingModalSurface};
use fs_system::{LoweringRegistry, ParameterProvenance, lowering::owners::register_owner_defaults};

const MYLAR: &[u8] =
    include_bytes!("../../fs-vkplate/packs/head-film-pet-mylar-dupont-matsuzawa1959.fsmatpk");
const HICKORY: &[u8] =
    include_bytes!("../../fs-vkplate/packs/stick-hickory-shagbark-fpl-gtr282.fsmatpk");
const N2: &[u8] =
    include_bytes!("../../../examples/assemblies/packs/nitrogen-n2-mcbride2002-lemmon2004.fsmcdpk");

fn pv(s: &[ParameterProvenance], n: &str) -> f64 {
    s.iter().find(|p| p.parameter == n).unwrap().value
}

fn couple_mic(s: &RadiatingModalSurface, mr: u32, gas: &GasState) -> BaffledPressure {
    // engine.rs:429 downward coords: flip -1 against our +z-up outward_sign on both heads.
    let surf: Vec<_> = s
        .samples
        .iter()
        .map(|p| SurfaceSample {
            position_m: [p.xy_m[0], p.xy_m[1], 0.0],
            area_m2: p.area_m2,
            mode_shape: p.shapes.iter().map(|v| -v).collect(),
        })
        .collect();
    BaffledPressure::new_multirate(
        &surf,
        s.mode_count,
        mr,
        48_000,
        [0.08, 0.05, 0.35], // params.rs:219-220
        RayleighMedium::from(gas),
    )
    .unwrap()
}

fn first_arrival_val(xs: &[f64]) -> (f64, f64) {
    let pk = xs.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
    let &v = xs.iter().find(|v| v.abs() > 1.0e-3 * pk).unwrap();
    (v, pk)
}

fn ed_cfg(mr: u32) -> LinearImpactConfig {
    let mut component = ModalAcousticTimeBudget::audible_reference();
    component.maximum_total_energy_j = 20.0;
    LinearImpactConfig {
        sample_rate_hz: mr,
        max_steps: 40_000,
        maximum_generalized_force: 1e6,
        component,
        coupling: ModalCouplingConfig {
            max_modes: 64,
            max_connections: 8,
            max_setup_terms: 100_000,
            nyquist_guard_fraction: 0.9,
            maximum_total_energy_j: 20.0,
            maximum_abs_pressure_pa: 1e6,
            maximum_abs_connection_force_n: 1e5,
            solve_relative_tolerance: 1e-10,
            energy_absolute_tolerance_j: 1e-9,
            energy_relative_tolerance: 1e-6,
        },
        contact: ModalContactConfig {
            max_iterations: 100,
            maximum_force_n: 1e5,
            maximum_penetration_m: 0.01,
            force_absolute_tolerance_n: 1e-10,
            force_relative_tolerance: 1e-10,
        },
        multiple: MultiContactConfig {
            max_contacts: 32,
            max_sweeps: 100,
            max_setup_terms: 100_000,
        },
    }
}

fn drum_src(ratio: u32, v0: f64, with_cavity: bool, h_scale: f64, t_scale: f64) -> String {
    let (mr, dt) = (48_000 * ratio, 1.0 / f64::from(48_000 * ratio));
    let (ten, h_m, amb) = (
        2500.0 * t_scale,
        0.3556 * h_scale,
        90_000.0 * PI * 0.1778 * 0.1778,
    );
    let cav = if with_cavity {
        format!(
            "(part \"resonant\" :kind shell :geometry (tensioned-disk :radius 0.1778m :thickness 0.00025m \
               :tension 2500N/m :radial-intervals 6 :azimuths 24) :material \"polyethylene-terephthalate\" \
               :query-temp 293.15K :physics (mechanics) :fidelity (modal 3)) \
             (part \"cavity\" :kind cavity :geometry (cylinder 0.1778m {h_m:.12}m) :material \"N2\" \
               :query-at (pressure 90000Pa) :physics (fluid) :fidelity lumped) \
             (interface \"f1\" :kind piston-face :a \"batter\" :b \"cavity\") \
             (interface \"f2\" :kind piston-face :a \"resonant\" :b \"cavity\" :sign -1) \
             (source \"amb1\" :part \"batter\" :port \"face\" :history (constant -{amb:.15}N)) \
             (source \"amb2\" :part \"resonant\" :port \"face\" :history (constant {amb:.15}N)) \
             (probe \"mr\" :kind microphone :target \"resonant\" :position (0.08m 0.05m 0.35m) :sample-rate 48000Hz)"
        )
    } else {
        String::new()
    };
    format!(
        "(assembly \"ed_drum\" (version 1) (units \"SI\") (seed 0x7) \
         (versions (constellation \"2026-07\") (workspace \"t\")) \
         (budget (wall 60s) (mem 64MiB)) \
         (capability :cores 1cores :mem 64MiB :wall 60s :ops (assembly.*)) \
         (reference-state assembled) \
         (environment :ambient-temp 293.15K :ambient-pressure 90000Pa :relative-humidity 0) \
         (time-plan :horizon 0.005s :dt {dt:.17}s :sample-rate {mr}Hz) \
         (part \"stick\" :kind rigid-body :geometry (sphere 0.015m) :material \"carya-ovata\" \
           :query-at (T 293.15K moisture_content 0.12) :physics (mechanics) :fidelity lumped \
           :initial (:position (0m 0m 0m) :orientation (1 0 0 0) \
             :linear-velocity (-{v0:.12}m/s 0m/s 0m/s) :angular-velocity (0rad/s 0rad/s 0rad/s))) \
         (part \"batter\" :kind shell :geometry (tensioned-disk :radius 0.1778m :thickness 0.00025m \
           :tension {ten:.12}N/m :radial-intervals 6 :azimuths 24) :material \"polyethylene-terephthalate\" \
           :query-temp 293.15K :physics (mechanics) :fidelity (modal 3)) \
         (interface \"hit\" :kind hertz-contact :a \"stick\" :b \"batter\" \
           :maximum-force 100000N :maximum-penetration 0.01m :gap 0.0002m :station-x 0.0889m :station-y 0m) \
         (probe \"mb\" :kind microphone :target \"batter\" :position (0.08m 0.05m 0.35m) :sample-rate 48000Hz) \
         {cav})"
    )
}

struct Case {
    rel_b: f64,
    rel_r: f64,
    pk_b: f64,
    pk_r: f64,
    bat_v: f64,
    res_v: f64,
    max_dv: f64,
    gamma: f64,
    e_res: f64,
}

#[allow(clippy::too_many_lines)]
fn run_case(
    packs: &PackCatalog,
    gas: &GasState,
    ratio: u32,
    v0: f64,
    with_cavity: bool,
    (h_scale, t_scale): (f64, f64),
) -> Case {
    let (mr, dt) = (48_000 * ratio, 1.0 / f64::from(48_000 * ratio));
    let lower_s = |s: String| {
        let a = admit_assembly_source(s.as_bytes(), packs).unwrap();
        let mut r = LoweringRegistry::with_kernel_defaults();
        register_owner_defaults(&mut r);
        let l = r.lower(&a.draft, &a.card_store).unwrap();
        (a.draft, l)
    };
    let (draft, mut low) = lower_s(drum_src(ratio, v0, with_cavity, h_scale, t_scale));
    low.receipt.verify_completeness().unwrap();
    let nom = (h_scale - 1.0).abs() == 0.0 && (t_scale - 1.0).abs() == 0.0;
    let ref_pair = (!nom).then(|| lower_s(drum_src(ratio, v0, with_cavity, 1.0, 1.0)));
    let (rd, rlow) = ref_pair.as_ref().map_or((&draft, &low), |(d, l)| (d, l));
    let sb = &rlow.radiating_surfaces[&rd.input().parts[1].entity];
    let mk_mic = |i| {
        AcousticMicrophone::from_draft_and_lowered(&draft.input().probes[i], &draft, &low, dt)
            .unwrap()
    };
    let (mut mb, mut mr_mic) = (mk_mic(0), with_cavity.then(|| mk_mic(1)));
    let m_st = pv(&rlow.receipt.components[0].parameter_sources, "mass");
    let e_st = pv(
        &rlow.receipt.interfaces[0].parameter_sources,
        "effective_modulus",
    );
    let r_st = pv(
        &rlow.receipt.interfaces[0].parameter_sources,
        "effective_radius",
    );
    let damp = |ci: usize| -> Vec<f64> {
        (0..3)
            .map(|k| {
                let z = pv(
                    &rlow.receipt.components[ci].parameter_sources,
                    &format!("modal_damping_ratio[{k}]"),
                );
                2.0 * z * sb.mode_omegas_rad_s[k]
            })
            .collect()
    };
    let area = |s: &RadiatingModalSurface| -> Vec<f64> {
        (0..3)
            .map(|k| {
                s.samples
                    .iter()
                    .map(|p| p.area_m2 * s.outward_sign * p.shapes[k])
                    .sum()
            })
            .collect()
    };
    let hit = sb
        .samples
        .iter()
        .find(|p| (p.xy_m[0] - 0.0889).abs() < 1e-9 && p.xy_m[1].abs() < 1e-9)
        .unwrap();
    let ndof = if with_cavity { 7 } else { 4 };
    let mut crow = vec![0.0; ndof];
    crow[0] = 1.0 / m_st.sqrt();
    for k in 0..3 {
        crow[1 + k] = -sb.outward_sign * hit.shapes[k];
    }
    let obs = Obstacle::new(
        crow,
        1,
        ndof,
        vec![2e-4],
        vec![1.0],
        (4.0 / 3.0) * e_st * r_st.sqrt(),
        1.5,
        "ed".into(),
    )
    .unwrap();
    let s0 = ModalAcousticState {
        displacement_m_sqrt_kg: 0.0,
        velocity_m_sqrt_kg_per_s: v0 * m_st.sqrt(),
    };
    let mut bodies = vec![
        ImpactBody {
            potential: BodyPotential::Linear(vec![0.0]),
            initial: vec![s0],
            damping_per_s: vec![0.0],
        },
        ImpactBody {
            potential: BodyPotential::Linear(sb.mode_omegas_rad_s.clone()),
            initial: vec![ModalAcousticState::default(); 3],
            damping_per_s: damp(1),
        },
    ];
    let (mut vols, mut gamma, mut ed_mr) = (Vec::new(), 1.4, None);
    let ref_a = PI * 0.1778 * 0.1778;
    if with_cavity {
        let sr = &rlow.radiating_surfaces[&rd.input().parts[2].entity];
        bodies.push(ImpactBody {
            potential: BodyPotential::Linear(sr.mode_omegas_rad_s.clone()),
            initial: vec![ModalAcousticState::default(); 3],
            damping_per_s: damp(2),
        });
        let g0 = rlow
            .system
            .state::<ClosedGasParcel>(rlow.node_map["cavity"])
            .unwrap()
            .state()
            .state();
        gamma = g0.gamma;
        let mut areas = vec![0.0];
        areas.extend(area(sb));
        areas.extend(area(sr).into_iter().map(|a| -a));
        vols.push(VolumeConnection {
            spring: VolumeSpring {
                bulk_modulus_pa: gamma * g0.pressure,
                volume_m3: ref_a * 0.3556,
                areas,
            },
            reference_area_m2: ref_a,
        });
        ed_mr = Some(couple_mic(sr, mr, gas));
    }
    let mut ed_mb = couple_mic(sb, mr, gas);
    let gate = CancelGate::new_clock_free();
    let mut ed =
        LinearImpactSystem::new_with_dampers(bodies, vec![obs], vols, vec![], ed_cfg(mr), &gate)
            .unwrap();
    let (mut pb, mut pr, mut eb, mut er, mut epb, mut epr, mut max_dv) = (
        Vec::new(),
        Vec::new(),
        0.0_f64,
        0.0_f64,
        0.0_f64,
        0.0_f64,
        0.0_f64,
    );
    let (mut tb, mut tr, ext) = (
        vec![0.0; ratio as usize * 3],
        vec![0.0; ratio as usize * 3],
        vec![0.0; ndof],
    );
    let v0_cav = ref_a * 0.3556 * h_scale;
    let key = StreamKey {
        seed: 7,
        kernel_id: 5,
        tile: 0,
        iteration: 0,
    };
    ArenaPool::new(ArenaConfig::default()).scope(|a| {
        let cx = Cx::new(&gate, a, key, Budget::INFINITE, ExecMode::Deterministic);
        for _ in 0..240 {
            for s in 0..ratio as usize {
                low.system.step(dt, &cx).unwrap();
                if with_cavity {
                    let v = low
                        .system
                        .state::<ClosedGasParcel>(low.node_map["cavity"])
                        .unwrap()
                        .volume_m3();
                    max_dv = max_dv.max(((v - v0_cav) / v0_cav).abs());
                }
                let gb = mb.advance(&low.system).unwrap();
                let gr = mr_mic
                    .as_mut()
                    .and_then(|m| m.advance(&low.system).unwrap());
                ed.step(&ext, &gate).unwrap();
                let st = ed.state();
                for k in 0..3 {
                    tb[s * 3 + k] = st[3 + 2 * k];
                    if with_cavity {
                        tr[s * 3 + k] = st[9 + 2 * k];
                    }
                }
                if let Some(gb) = gb {
                    let qb = ed_mb.step_trace(&tb).unwrap();
                    pb.push(gb);
                    (eb, epb) = (eb.max((gb - qb).abs()), epb.max(qb.abs()));
                }
                if let (Some(gr), Some(em)) = (gr, ed_mr.as_mut()) {
                    let qr = em.step_trace(&tr).unwrap();
                    pr.push(gr);
                    (er, epr) = (er.max((gr - qr).abs()), epr.max(qr.abs()));
                }
            }
        }
    });
    let (bat_v, _) = first_arrival_val(&pb);
    let res_v = if with_cavity {
        first_arrival_val(&pr).0
    } else {
        0.0
    };
    Case {
        rel_b: eb / epb,
        rel_r: if with_cavity { er / epr } else { 0.0 },
        pk_b: epb,
        pk_r: epr,
        bat_v,
        res_v,
        max_dv,
        gamma,
        e_res: low.system.ledger().energy_residual_j().abs(),
    }
}

#[test]
fn fsim_two_head_drum_matches_engine_drum_mechanics_and_both_microphones() {
    let packs = PackCatalog::from_raw_bytes(&[MYLAR, HICKORY, N2]).unwrap();
    let gas = GasState::try_new_moist_air(293.15, 90_000.0, 0.0).unwrap();

    // 1. No-cavity isolation (stick + HertzContact + ModalStructure + AcousticMicrophone):
    //    matches LinearImpactSystem + BaffledPressure to contact solver tolerance (~1e-10).
    let nc = run_case(&packs, &gas, 1, 2.0, false, (1.0, 1.0));
    eprintln!("NO-CAVITY 1x: rel_b={:.2e} pk_b={:.2}Pa", nc.rel_b, nc.pk_b);
    assert!(nc.pk_b > 1.0 && nc.rel_b <= 1.0e-8);

    // 2. Full 2-head N2 drum at 1x, 2x, 4x: residual is flat from 2x to 4x because the dominant
    //    difference is model form: GasCavity's nonlinear NASA-9 adiabatic ideal gas P(V)=P0(V0/V)^gamma
    //    vs engine-drum's linearized VolumeSpring (K=gamma*P0), whose leading 2nd-order relative
    //    pressure term is pred_nl = ((gamma+1)/2)*max|dV/V0|.
    //    Margin is 1.25 * pred_nl (not noise slack): the +0.25*pred_nl covers omitted higher-order
    //    deterministic terms — cubic (dV/V0)^2 expansion, NASA-9 gamma(T) shift across dT, and
    //    step-held port pressure vs LinearImpactSystem's continuous exp(A dt) VolumeSpring.
    let d1 = run_case(&packs, &gas, 1, 2.0, true, (1.0, 1.0));
    let d2 = run_case(&packs, &gas, 2, 2.0, true, (1.0, 1.0));
    let d4 = run_case(&packs, &gas, 4, 2.0, true, (1.0, 1.0));
    let d_small = run_case(&packs, &gas, 1, 0.2, true, (1.0, 1.0));
    let pred = |d: &Case| 0.5 * (d.gamma + 1.0) * d.max_dv;
    let bound = 1.25 * pred(&d1);
    for (tag, d) in [("1x", &d1), ("2x", &d2), ("4x", &d4), ("0.2m/s", &d_small)] {
        eprintln!(
            "DRUM {tag:6}: rel_b={:.3e} rel_r={:.3e} pred_nl={:.3e} (max|dV/V0|={:.3e})",
            d.rel_b,
            d.rel_r,
            pred(d),
            d.max_dv
        );
    }
    eprintln!("DRUM BOUND: 1.25 * pred_nl(1x) = {bound:.3e}");
    assert!(bound < 2.0e-3 && d1.pk_b > 1.0 && d1.pk_r > 0.1);
    for d in [&d1, &d2, &d4] {
        assert!(d.rel_b <= bound && d.rel_r <= bound && d.e_res < 1.0e-10);
    }
    // Scaling proof: 10x gentler strike -> smaller |dV/V0| -> residual drops well below full strike.
    assert!(d_small.max_dv < 0.2 * d1.max_dv && d_small.rel_b < 0.3 * d1.rel_b);
    // Independent physical sign oracles: downward strike pushes batter away (<0), resonant toward (>0).
    assert!(d1.bat_v < 0.0 && d1.res_v > 0.0);

    // 3. Teeth: planted physical errors (cavity volume x1.1, batter tension x1.05) fail by >30x.
    let bad_vol = run_case(&packs, &gas, 1, 2.0, true, (1.1, 1.0));
    let bad_ten = run_case(&packs, &gas, 1, 2.0, true, (1.0, 1.05));
    eprintln!(
        "TEETH vol*1.1: {:.2e}/{:.2e} ({:.0}x/{:.0}x) | ten*1.05: {:.2e}/{:.2e} ({:.0}x/{:.0}x)",
        bad_vol.rel_b,
        bad_vol.rel_r,
        bad_vol.rel_b / bound,
        bad_vol.rel_r / bound,
        bad_ten.rel_b,
        bad_ten.rel_r,
        bad_ten.rel_b / bound,
        bad_ten.rel_r / bound
    );
    assert!(bad_vol.rel_b > 30.0 * bound && bad_vol.rel_r > 30.0 * bound);
    assert!(bad_ten.rel_b > 30.0 * bound && bad_ten.rel_r > 30.0 * bound);
}
