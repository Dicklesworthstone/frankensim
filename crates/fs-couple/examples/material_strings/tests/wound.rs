//! Material-only winding swaps must change the actual reciprocal pressure path.
use super::*;

#[path = "wound/stretching_memory.rs"]
mod stretching_memory;

#[path = "wound/modal_loss.rs"]
mod modal_loss;

fn wound_spec() -> input::Spec {
    input::Spec::parse(include_str!("../wound-pluck.material-strings")
        .replace("48000 4096 0.01 128", "48000 256 0.01 128").as_bytes()).unwrap()
}
fn wound_ready(cx: &Cx<'_>, gate: &CancelGate, spec: input::Spec, density: f64, band: f64) -> Result<Ready> {
    spec.compile(cx, gate, |name| Ok(if name == "wrap.fsmatpk" {
        // No wrap bending viscosity is needed or silently borrowed from core.
        pack(name, density, band, Some(ETA))
    } else { pack(name, 7800.0, 200_000.0, None) }))
}
fn specimen(ready: &Ready) -> &fs_couple::material_assembly::wound::ResolvedWoundStringSpecimen {
    let source = &ready.performance.bank().source().bank().strings()[0].source;
    let MaterialStringSource::Wound(specimen) = source else { panic!("original winding source retained") };
    specimen
}

#[test]
fn winding_material_swap_reaches_mass_modes_gas_loading_and_real_pressure() {
    scoped(|cx, gate| {
        let mut light = wound_ready(cx, gate, wound_spec(), 8900.0, 200_000.0).unwrap();
        let mut dense = wound_ready(cx, gate, wound_spec(), 17800.0, 200_000.0).unwrap();
        assert_eq!(light.spec.bytes, dense.spec.bytes);
        assert_eq!(specimen(&light).section().core(), specimen(&dense).section().core());
        assert_ne!(specimen(&light).section().winding_material().identity(), specimen(&dense).section().winding_material().identity());
        assert_eq!(specimen(&light).string().tension_n, 20.0);
        assert_eq!(specimen(&dense).string().tension_n, 20.0);
        assert_eq!(specimen(&light).string().bending_stiffness_n_m2,
            specimen(&dense).string().bending_stiffness_n_m2);
        assert!(dense.performance.bank().cards()[0].fundamental_hz() < light.performance.bank().cards()[0].fundamental_hz());
        let original = pressure(&mut light, 37);
        let changed = pressure(&mut dense, 1);
        assert!(original.iter().zip(&changed).any(|(a, b)| (a - b).abs() > 1e-12));
        assert!(original.iter().any(|p| p.abs() > 1e-10));
        let frame = dense.performance.bank().system().last_frame().unwrap();
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        let mut replay = wound_ready(cx, gate, wound_spec(), 8900.0, 200_000.0).unwrap();
        assert_eq!(pressure(&mut replay, 128), original);
    });
}

#[test]
fn wound_extension_uses_spring_plus_core_and_gas_radius_is_the_finished_radius() {
    scoped(|cx, gate| {
        let text = String::from_utf8(wound_spec().bytes).unwrap().replace("tension 20", "extension 0.4998 0.001");
        let ready = wound_ready(cx, gate, input::Spec::parse(text.as_bytes()).unwrap(), 8900.0, 200_000.0).unwrap();
        let source = specimen(&ready);
        let strain = (0.5 - 0.4998) / 0.4998;
        assert_eq!(source.string().tension_n, source.section().axial_stiffness_n() * strain);
        assert!(source.section().axial_stiffness_n() > source.section().core().elastic().axial_stiffness_n());
        assert_eq!(source.string().width_m, 2.0 * source.section().outer_radius_m());
        // Source and mechanical state survive the existing cylinder gas owner.
        let gas_text = text.replace("cavity neck", "cylinder-gas-loading 0.001 30 4\ncavity neck");
        let loaded = wound_ready(cx, gate, input::Spec::parse(gas_text.as_bytes()).unwrap(), 8900.0, 200_000.0).unwrap();
        assert_eq!(specimen(&loaded).section(), source.section());
    });
}

#[test]
fn winding_geometry_source_band_and_unmodeled_torsion_refuse() {
    scoped(|cx, gate| {
        assert!(wound_ready(cx, gate, wound_spec(), 8900.0, 1000.0).is_err());
        let text = String::from_utf8(wound_spec().bytes).unwrap();
        for changed in [text.replace("0.00005 0.000101", "0.00005 0.0001"),
            text.replace("20 4 0 0.8", "20 4 2 0.8")] {
            let spec = input::Spec::parse(changed.as_bytes()).unwrap();
            assert!(wound_ready(cx, gate, spec, 8900.0, 200_000.0).is_err());
        }
        assert!(input::Spec::parse(text.replace("0.00005 0.000101", "0.00004 0.000101").as_bytes()).is_err());
    });
}


const CORE_MEMORY: &str = "bending-loss prony 1 relaxing_bending_modulus bending_relaxation_time";

fn memory_spec(additional: bool) -> input::Spec {
    let declaration = if additional {
        "bending-memory 8 1\nstring-memory 0 1 relaxing_bending_modulus bending_relaxation_time"
    } else { CORE_MEMORY };
    input::Spec::parse(String::from_utf8(wound_spec().bytes).unwrap()
        .replace("cavity neck", &format!("{declaration}\ncavity neck")).as_bytes()).unwrap()
}

/// The core's relaxed modulus, bit-equal to `pack`'s young_modulus.
const EQUILIBRIUM_MODULUS: f64 = 200e9;

fn claimed(base: &NormalizedPack, version: &str, source: &str, band: f64,
    values: &[(&str, Dims, f64)]) -> NormalizedPack {
    let mut claims = base.claims().clone();
    for &(key, dims, value) in values {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(key, dims), value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0).with("omega", 0.0, band),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity, observations: vec![],
            provenance: Provenance { source: source.into(), license: "test-only".into(), artifact: None },
        }).unwrap();
    }
    NormalizedPack::new(base.pack_id(), version,
        base.source_artifact(), base.redistribution_terms(), claims, vec![], vec![]).unwrap()
}

fn memory_pack(name: &str, delta: f64, tau: f64, band: f64, eta: Option<f64>) -> NormalizedPack {
    if name != "wire.fsmatpk" {
        // No relaxation modulus/time (or viscosity) is supplied by the wrap.
        return pack(name, if name == "wrap.fsmatpk" { 8900.0 } else { 7800.0 },
            band, (name == "wrap.fsmatpk").then_some(ETA));
    }
    let mut values = vec![
        ("equilibrium_young_modulus", Pressure::DIMS, EQUILIBRIUM_MODULUS),
        ("relaxing_bending_modulus", Pressure::DIMS, delta),
        ("bending_relaxation_time", Dims([0, 0, 1, 0, 0, 0]), tau),
    ];
    if let Some(eta) = eta { values.push((ETA, DynViscosity::DIMS, eta)); }
    claimed(&pack(name, 7800.0, band, Some(ETA)), "synthetic-wound-core-memory-v1",
        "synthetic wound-core memory integration", band, &values)
}

#[test]
fn wound_core_spectrum_uses_helix_mass_and_drives_shared_pressure_without_a_viscosity_claim() {
    scoped(|cx, gate| {
        let build = |delta| memory_spec(false).compile(cx, gate, |name|
            Ok(memory_pack(name, delta, 0.001, 200_000.0, None))).unwrap();
        let mut a = build(80e9);
        let mut b = build(160e9);
        let core = specimen(&a).section().core();
        let fs_material::section::winding::WoundCoreSection::Prony(core) = core else {
            panic!("the original core spectrum must remain attached");
        };
        assert!(core.elastic().layers()[0].material.property(ETA).is_none());
        assert!(specimen(&a).string().kelvin_voigt_bending.is_none());
        assert_eq!(core.properties()[0].branches,
            vec![("relaxing_bending_modulus".into(), "bending_relaxation_time".into())]);
        assert_eq!(specimen(&a).section().winding_material().identity(),
            specimen(&b).section().winding_material().identity());
        let card = &a.performance.bank().cards()[0];
        let branches = a.performance.bank().system().material_relaxation_branches();
        assert_eq!(branches.len(), card.mode_count);
        // Curvature energy integrated independently over the physical core:
        // (1/2) delta_E integral_A(z^2)dA integral_L(y_xx^2)dx.
        // Modal inertia includes the actual helical wrapping wire.
        let radius = core.elastic().outer_radius_m();
        let moment = std::f64::consts::PI * radius.powi(4) / 4.0;
        let mass = specimen(&a).section().line_density_kg_m();
        assert_eq!(card.linear_density_kg_m, mass);
        for (mode, branch) in branches.iter().enumerate() {
            let k = (mode + 1) as f64 * std::f64::consts::PI / card.length_m;
            let expected = 80e9 * moment * k.powi(4) / mass;
            assert!((branch.stiffness / expected - 1.0).abs() < 2e-14);
            assert_eq!(branch.relaxation_time_s, 0.001);
            assert_eq!(branch.projection[2 * mode], 1.0);
            assert_eq!(branch.projection.iter().filter(|v| **v != 0.0).count(), 1);
        }
        assert!(a.performance.bank().system().material_memory().iter().any(|z| *z != 0.0));
        assert_eq!(a.spec.bytes, b.spec.bytes);
        assert_eq!(a.performance.bank().system().components()[0].modes(),
            b.performance.bank().system().components()[0].modes());
        let p = pressure(&mut a, 37);
        let q = pressure(&mut b, 1);
        assert!(p.iter().any(|x| x.abs() > 1e-10));
        assert!(p.iter().zip(&q).any(|(x, y)| (x - y).abs() > 1e-12));
        assert_eq!(p, pressure(&mut build(80e9), 128));
        let expected = fs_couple::pcm_wav::encode_pcm16_wav(&p, a.spec.rate, a.spec.full_scale).unwrap().0;
        let (wav, _) = build(80e9).write(Cursor::new(Vec::<u8>::new()), gate).unwrap();
        assert_eq!(wav.into_inner(), expected);
        for ready in [&a, &b] {
            let frame = ready.performance.bank().system().last_frame().unwrap();
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        }
    });
}

#[test]
fn wound_core_memory_is_equivalent_through_original_or_layered_binding_with_two_planes_and_gas() {
    scoped(|cx, gate| {
        let build = |additional| {
            let mut spec = memory_spec(additional);
            spec.samples = 96;
            spec.strings[0].polarizations = Some(
                fs_couple::bowed_string::runtime::reciprocal::bank::polarization::TransverseDirections::try_new(
                    [1.0, 0.0], [0.8, 0.6]).unwrap());
            spec.cylinder_gas_loading = Some(CylinderFlowResolution {
                minimum_wavenumber_radius: 0.001, maximum_wavenumber_radius: 30.0, intervals: 2,
            });
            spec.compile(cx, gate, |name| Ok(memory_pack(name, 80e9, 0.001,
                200_000.0, additional.then_some(0.0)))).unwrap()
        };
        let mut original = build(false);
        let mut layered = build(true);
        let source = original.performance.bank().source().bank();
        let fluid = source.cylinder_gas_loading().unwrap();
        let card = &original.performance.bank().cards()[0];
        assert!(card.linear_density_kg_m > fluid.strings()[0].solid_card().linear_density_kg_m);
        assert_eq!(original.performance.bank().system().material_memory().len(), 8);
        assert_eq!(layered.performance.bank().system().material_memory().len(), 8);
        // The two APIs encode the same energy channel with different coordinate
        // scaling (k in stiffness versus sqrt(k) in projection).
        let a = original.performance.bank().system().material_relaxation_branches();
        let b = layered.performance.bank().system().material_relaxation_branches();
        for (left, right) in a.iter().zip(b) {
            assert_eq!(left.relaxation_time_s, right.relaxation_time_s);
            for (x, y) in left.projection.iter().zip(&right.projection) {
                let lhs = left.stiffness * x * x;
                let rhs = right.stiffness * y * y;
                assert!((lhs - rhs).abs() <= 2e-14 * lhs.abs().max(rhs.abs()).max(1.0));
            }
        }
        let p = pressure(&mut original, 1);
        let q = pressure(&mut layered, 37);
        let peak = p.iter().map(|x| x.abs()).fold(0.0, f64::max);
        assert!(peak > 1e-12);
        assert!(p.iter().zip(&q).all(|(x, y)| (x - y).abs() <= 2e-11 * peak));
        for ready in [&original, &layered] {
            let frame = ready.performance.bank().system().last_frame().unwrap();
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        }
    });
}

#[test]
fn wound_core_memory_renders_through_the_stiff_relaxation_band() {
    // dt/tau = 2083 (the arm validated with f64 doubling), 2600 and 3600 (the
    // band where its rounding refused the exact step) and 20833 (1 ns).
    scoped(|cx, gate| {
        let spec = memory_spec(false);
        let dt = 1.0 / f64::from(spec.rate);
        let mut refused = Vec::new();
        for tau in [1e-8, dt / 2600.0, dt / 3600.0, 1e-9] {
            let mut ready = spec.clone().compile(cx, gate, |name|
                Ok(memory_pack(name, 80e9, tau, 200_000.0, None))).unwrap();
            let mut p = vec![0.0; ready.spec.samples as usize];
            if let Some(error) = p.chunks_mut(37).find_map(|chunk| ready.performance.block(chunk).err()) {
                refused.push(format!("dt/tau {:.0}: {error:?}", dt / tau));
                continue;
            }
            assert!(p.iter().all(|x| x.is_finite()) && p.iter().any(|x| x.abs() > 1e-10), "dt/tau {}", dt / tau);
            let frame = ready.performance.bank().system().last_frame().unwrap();
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j,
                "dt/tau {}: {} > {}", dt / tau, frame.energy_residual_j, frame.energy_tolerance_j);
        }
        assert!(refused.is_empty(), "{refused:?}");
    });
}

/// A test-only alteration of the comparisons below (never of the engine), which
/// each check must refuse.
#[derive(Clone, Copy, PartialEq)]
enum Plant { Faithful, StorageDropped, ArmDoubled }

fn kelvin_voigt_limit(plant: Plant) {
    // An arm (delta, tau) on bending strain e carries f + tau f' = delta tau e', so
    // f = eta e' - eta tau e'' + ..., eta = delta tau; with its loss rate f^2 / eta and
    // storage H = f^2 / (2 delta), from relaxed memory,
    //   D_arm = D_KV - 2 H_arm - eta tau^2 int e''^2 + O(tau^3),
    // D_KV = eta int e'^2, the loss of a core viscosity eta on the same motion. At
    // fixed eta that is D(r) = D_inf + b/r in r = dt/tau: D_inf = D_KV and
    // b/r = -2 H_arm, measured, not fitted. So A = D_arm + 2 H_arm does not depend on r
    // to that order.
    //
    // Gate, every sample: arms (delta, tau) and (delta/2, 2 tau) share eta, so both
    // remove energy at the same leading rate eta e'^2 and their motions part only
    // through the terms below. Allowance on A_slow / A_fast - 1 at t = k dt,
    // Omega = pi f_s (the Nyquist guard on every admitted frequency; also above
    // sqrt(3)/t, the e''/e' ratio of motion started from rest within a held-force step):
    //  (a) the tau^2 remainder, at most (Omega tau)^2 for one arm and (2 Omega tau)^2
    //      for the other: 5 (Omega tau)^2;
    //  (b) an arm's residual stiffness, delta (w tau)^2 / (1 + (w tau)^2), is at most
    //      q = (delta / E_eq)(Omega tau)^2 of a mode's stiffness for (delta, tau) and 2 q
    //      for (delta/2, 2 tau). It moves a frequency by dw <= w q / 2 (twice that); a
    //      mode's own loss then moves by <= 4 dw / w and a pair's beat by <= max(dw) t:
    //      3 q (2 + Omega t / 2);
    //  (c) rounding in each arm's positive heat factor, doubled s times in f64 outside
    //      the compensated flow: 2 (2^s - 1) gamma_n per channel, as in gramian.rs, with
    //      n at most the exact path's dense block cap. The arms need s <= 17 here;
    //      2^18 covers the rest. Two channels.
    //
    // Scale, sample 1: a ratio of arms cannot see a wrongly scaled arm, so the fast arm
    // meets a core viscosity eta in an otherwise equal card, all channels summed. Every
    // other channel, G (chiefly the string's own gas drag, in the component channel), is
    // real physics in both cards and dwarfs D_KV, so no channel can be asserted zero.
    // Beside (a) and (b) on D_KV, the allowance carries (b) on G; G's response to the
    // arm-versus-core dissipative difference, (Omega tau)^2 of the core's own reach on
    // G, which a core decaying a mode at most at c = (delta / E_eq) tau Omega^2 bounds by
    // 2 c t; and (c) on every channel of both cards. Its beat term grows with t times
    // G / D_KV, so only k = 1 is checked; the gate holds the law from there.
    scoped(|cx, gate| {
        let (delta, tau, band) = (80e9, 1e-9, 200_000.0);
        // A lossless body plate: the arms are then the only material loss (the core
        // card's material channel is asserted zero every sample).
        let body = claimed(&pack("body.fsmatpk", 7800.0, band, Some(ETA)),
            "synthetic-lossless-body-v1", "synthetic lossless body", band, &[(ETA, DynViscosity::DIMS, 0.0)]);
        let card = |spec: input::Spec, delta: f64, tau: f64, eta: Option<f64>| {
            let packs = |name: &str| if name == "body.fsmatpk" { body.clone() }
                else { memory_pack(name, delta, tau, band, eta) };
            spec.compile(cx, gate, |name| Ok(packs(name))).unwrap()
        };
        let mut fast = card(memory_spec(false), delta, tau, None);
        let mut slow = card(memory_spec(false), 0.5 * delta, 2.0 * tau, None);
        let mut viscous = card(wound_spec(), delta, tau, Some(delta * tau));
        let rate = f64::from(fast.spec.rate);
        let (dt, omega) = (1.0 / rate, std::f64::consts::PI * rate);
        for ready in [&fast, &slow, &viscous] {
            assert!(ready.performance.bank().system().linear_frequency_upper_bound_rad_s() <= omega);
        }
        let nu = fs_phs::linear_exact::MAX_LINEAR_EXPONENTIAL_BLOCK_DIM as f64 * (0.5 * f64::EPSILON);
        let rounding = 2.0 * 2.0_f64.powi(18) * (nu / (1.0 - nu));
        let storage = if plant == Plant::StorageDropped { 0.0 } else { 2.0 };
        let arm = if plant == Plant::ArmDoubled { 2.0 } else { 1.0 };
        let (mut fast_loss, mut slow_loss) = (0.0, 0.0);
        for k in 1..=fast.spec.samples {
            for ready in [&mut fast, &mut slow, &mut viscous] {
                ready.performance.block(&mut [0.0]).unwrap();
            }
            let a = fast.performance.bank().system().coupled_frame().unwrap();
            let b = slow.performance.bank().system().coupled_frame().unwrap();
            let v = viscous.performance.bank().system().coupled_frame().unwrap();
            assert_eq!(v.material_dissipation_j, 0.0, "sample {k}");
            fast_loss += arm * a.material_dissipation_j;
            slow_loss += arm * b.material_dissipation_j;
            let fast_total = fast_loss + storage * arm * a.material_energy_j;
            let slow_total = slow_loss + storage * arm * b.material_energy_j;
            let t = k as f64 * dt;
            let remainder = (omega * tau).powi(2);
            let stiffness = delta / EQUILIBRIUM_MODULUS * remainder * (2.0 + 0.5 * omega * t);
            let allowance = 5.0 * remainder + 3.0 * stiffness + 2.0 * rounding;
            let defect = slow_total / fast_total - 1.0;
            assert!(defect.abs() <= allowance,
                "sample {k}: relaxation-time law defect {defect:e} > {allowance:e}");
            if k == 1 {
                let rest = a.component_dissipation_j + a.connection_dissipation_j + a.radiation_dissipation_j;
                let core = v.component_dissipation_j + v.connection_dissipation_j + v.radiation_dissipation_j;
                let reach = 2.0 * delta / EQUILIBRIUM_MODULUS * tau * omega * omega * t;
                let allowance = (remainder + stiffness) * fast_total + (stiffness + remainder * reach) * rest
                    + rounding * (arm * a.material_dissipation_j + rest + core);
                let defect = fast_total + rest - core;
                assert!(defect.abs() <= allowance, "core scale defect {defect:e} J > {allowance:e} J");
            }
            for frame in [a, b, v] {
                assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j, "sample {k}");
            }
        }
    });
}

#[test]
fn wound_core_memory_at_one_nanosecond_meets_its_kelvin_voigt_limit() {
    kelvin_voigt_limit(Plant::Faithful);
}

#[test]
#[should_panic(expected = "relaxation-time law")]
fn kelvin_voigt_limit_refuses_an_arm_without_its_storage() {
    kelvin_voigt_limit(Plant::StorageDropped);
}

#[test]
#[should_panic(expected = "core scale")]
fn kelvin_voigt_limit_refuses_a_doubled_arm() {
    kelvin_voigt_limit(Plant::ArmDoubled);
}

#[test]
fn wound_core_memory_keeps_source_bands_and_refuses_a_memoryless_card_or_duplicate_owner() {
    scoped(|cx, gate| {
        let spec = memory_spec(false);
        let ready = spec.clone().compile(cx, gate, |name|
            Ok(memory_pack(name, 80e9, 0.001, 200_000.0, None))).unwrap();
        assert!(fs_couple::bowed_string::BowedStringCard::from_wound_material_specimen(
            cx, specimen(&ready), &ready.spec.gas, ready.spec.rate).is_err());
        for (delta, tau, band) in [(-1.0, 0.001, 200_000.0),
            (80e9, 0.0, 200_000.0), (80e9, 0.001, 1000.0)] {
            assert!(spec.clone().compile(cx, gate, |name|
                Ok(memory_pack(name, delta, tau, band, None))).is_err(), "accepted {delta} {tau} {band}");
        }
        // A wrap remains part of the instantaneous source-band admission even
        // though no relaxation arm is projected from its material.
        assert!(spec.clone().compile(cx, gate, |name| Ok(memory_pack(name,
            80e9, 0.001, if name == "wrap.fsmatpk" { 1000.0 } else { 200_000.0 }, None))).is_err());
        let mut duplicate = spec;
        let plain = wound_spec();
        duplicate.bending_memory = Some(bending_memory::Input::new(4, vec![
            bending_memory::Input::uniform_selection(0,
                vec![("relaxing_bending_modulus".into(), "bending_relaxation_time".into())],
                &plain.strings).unwrap(),
        ], &plain.strings).unwrap());
        let mut calls = 0;
        assert!(duplicate.compile(cx, gate, |name| {
            calls += 1; Ok(memory_pack(name, 80e9, 0.001, 200_000.0, Some(0.0)))
        }).is_err());
        assert_eq!(calls, 0);
    });
}
