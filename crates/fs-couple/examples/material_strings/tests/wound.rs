//! Material-only winding swaps must change the actual reciprocal pressure path.
use super::*;

#[path = "wound/stretching_memory.rs"]
mod stretching_memory;

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

fn memory_pack(name: &str, delta: f64, tau: f64, band: f64, eta: Option<f64>) -> NormalizedPack {
    if name != "wire.fsmatpk" {
        // No relaxation modulus/time (or viscosity) is supplied by the wrap.
        return pack(name, if name == "wrap.fsmatpk" { 8900.0 } else { 7800.0 },
            band, (name == "wrap.fsmatpk").then_some(ETA));
    }
    let base = pack(name, 7800.0, band, Some(ETA));
    let mut claims = base.claims().clone();
    let mut values = vec![
        ("equilibrium_young_modulus", Pressure::DIMS, 200e9),
        ("relaxing_bending_modulus", Pressure::DIMS, delta),
        ("bending_relaxation_time", Dims([0, 0, 1, 0, 0, 0]), tau),
    ];
    if let Some(eta) = eta { values.push((ETA, DynViscosity::DIMS, eta)); }
    for (key, dims, value) in values {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(key, dims), value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0).with("omega", 0.0, band),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity, observations: vec![],
            provenance: Provenance { source: "synthetic wound-core memory integration".into(),
                license: "test-only".into(), artifact: None },
        }).unwrap();
    }
    NormalizedPack::new(base.pack_id(), "synthetic-wound-core-memory-v1",
        base.source_artifact(), base.redistribution_terms(), claims, vec![], vec![]).unwrap()
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
fn wound_core_memory_keeps_source_bands_and_refuses_a_memoryless_card_or_duplicate_owner() {
    scoped(|cx, gate| {
        let spec = memory_spec(false);
        let ready = spec.clone().compile(cx, gate, |name|
            Ok(memory_pack(name, 80e9, 0.001, 200_000.0, None))).unwrap();
        assert!(fs_couple::bowed_string::BowedStringCard::from_wound_material_specimen(
            cx, specimen(&ready), &ready.spec.gas, ready.spec.rate).is_err());
        for (delta, tau, band) in [(-1.0, 0.001, 200_000.0),
            (80e9, 1e-9, 200_000.0), (80e9, 0.001, 1000.0)] {
            assert!(spec.clone().compile(cx, gate, |name|
                Ok(memory_pack(name, delta, tau, band, None))).is_err());
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
