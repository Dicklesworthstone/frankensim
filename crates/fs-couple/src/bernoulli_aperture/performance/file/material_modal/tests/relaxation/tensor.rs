//! Complete source tensors reach the original native reed and microphone.
use super::*;
use crate::bernoulli_aperture::dg::loaded::performance::LoadedDgPerformance;
use crate::bernoulli_aperture::performance::file::material_modal::dg::MaterialModalDgValvePerformance;

const SOURCE: &str =
    include_str!("../../../../../../../examples/material_modal_valve/paired-tensor-dg.performance");
const NETWORK: &str =
    include_str!("../../../../../../../examples/material_modal_valve/paired.performance");

fn tensor_source() -> String {
    SOURCE.replace("audio 96000 96000 2048 1", "audio 96000 96000 512 1")
}
fn tensor_network() -> String {
    let prefix = NETWORK
        .split_once("valve 0.0004\n")
        .unwrap()
        .0
        .replace("audio 768000 48000 2400 0.005", "audio 96000 96000 512 1");
    format!(
        "{prefix}valve 0.0004\n{}",
        tensor_source().split_once("valve 0.0004\n").unwrap().1
    )
    .replace("\ncontact ", "\ncoupling 1 1e-16 1e-7 1e-16 1e-8\ncontact ")
}

fn tensor_pack(c33: f64, time: f64, missing: &str, upper: f64) -> NormalizedMaterialCardPack {
    let mut claims = ClaimSet::new();
    let mut identity = [c33.to_le_bytes(), time.to_le_bytes(), upper.to_le_bytes()].concat();
    identity.extend_from_slice(missing.as_bytes());
    let artifact = hash_domain("org.frankensim.test.full-tensor-modal-valve", &identity);
    let provenance = Provenance {
        source: "manufactured full-tensor material fixture".into(),
        license: "CC0-1.0".into(),
        artifact: None,
    };
    let observation = claims
        .register_observation(ObservationDataset {
            specimen: "manufactured tensor blade".into(),
            method: "declared physical regression".into(),
            artifact,
            caveats: "not measured cane or an instrument qualification".into(),
            provenance: provenance.clone(),
        })
        .unwrap();
    for (name, value, dims) in [
        ("density", 900.0, Density::DIMS),
        ("young_modulus", 4e9, Pressure::DIMS),
        ("poisson_ratio", 0.3, Dims::NONE),
        ("strain_limit", 0.02, Dims::NONE),
        ("relaxation_time", time, Dims([0, 0, 1, 0, 0, 0])),
        ("branch_11", 1e8, Pressure::DIMS),
        ("branch_22", 2e8, Pressure::DIMS),
        ("branch_33", c33, Pressure::DIMS),
        ("branch_44", 1.5e8, Pressure::DIMS),
        ("branch_55", 1.2e8, Pressure::DIMS),
        ("branch_66", 0.8e8, Pressure::DIMS),
        ("branch_13", 4e7, Pressure::DIMS),
        ("branch_24", 2e7, Pressure::DIMS),
        ("branch_15", -1e7, Pressure::DIMS),
        ("zero_modulus", 0.0, Pressure::DIMS),
    ] {
        if name == missing {
            continue;
        }
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained()
                    .with("T", 280.0, 320.0)
                    .with("omega", 0.0, upper),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                provenance: provenance.clone(),
                observations: vec![observation.clone()],
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic".into(),
            phase: "solid".into(),
            process: "manufactured".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let normalized = NormalizedPack::new(
        "synthetic-tensor-reed",
        "synthetic-v1",
        artifact,
        "CC0-1.0",
        card.claims().clone(),
        vec![],
        vec![],
    )
    .unwrap();
    NormalizedMaterialCardPack::new(card.id().clone(), normalized).unwrap()
}
fn build_dg(text: &str, block: usize, c33: f64, time: f64) -> MaterialModalDgValvePerformance {
    MaterialModalDgValvePerformance::from_bytes_with_sources(
        text.as_bytes(),
        block,
        &CancelGate::new_clock_free(),
        None,
        |path| {
            assert_eq!(path, "tensor-blade.matpack");
            Ok(tensor_pack(c33, time, "", 100000.0))
        },
    )
    .unwrap()
}
fn dg_valve(renderer: &LoadedDgPerformance) -> &ModalAperture {
    let LoadedAperture::Modal(valve) = renderer.system().aperture() else {
        panic!("original material modal owner");
    };
    valve
}
fn render_dg(source: MaterialModalDgValvePerformance) -> (Vec<f64>, LoadedDgPerformance) {
    let samples = source.info().samples as usize;
    let mut renderer = source.into_renderer();
    let mut pressure = vec![0.0; samples];
    for chunk in pressure.chunks_mut(renderer.config().max_block) {
        renderer.block(chunk).unwrap();
    }
    (pressure, renderer)
}
fn build_network(text: &str, block: usize, c33: f64, time: f64) -> MaterialModalValvePerformance {
    MaterialModalValvePerformance::from_bytes_with_sources(
        text.as_bytes(),
        block,
        &CancelGate::new_clock_free(),
        None,
        |path| {
            assert_eq!(path, "tensor-blade.matpack");
            Ok(tensor_pack(c33, time, "", 100000.0))
        },
    )
    .unwrap()
}

#[test]
fn g3_tensor_source_normal_coupling_and_rates_reach_the_actual_dg_microphone() {
    let text = tensor_source();
    let source = build_dg(&text, 37, 3e8, 0.001);
    assert_eq!(
        (
            source.info().blades,
            source.info().modes,
            source.info().nodes,
            source.info().triangles
        ),
        (2, 4, 8, 4)
    );
    assert_eq!(source.material_records().len(), 2);
    assert!(source.renderer().baffled_receiver().is_some());
    assert!(source.renderer().observer_gas().is_some());
    let valve = dg_valve(source.renderer());
    let original_spec = valve.spec().clone();
    let memory = valve.material_relaxation().unwrap();
    assert!(memory.material_spectra().is_empty());
    assert_eq!(memory.tensor_material_points().len(), 12);
    assert_eq!(memory.memory_sqrt_j(), &[vec![0.0; 72]]);
    assert_eq!(
        memory
            .tensor_material_points()
            .iter()
            .map(|p| p.memory_range().len())
            .sum::<usize>(),
        72
    );
    let (pressure, rendered) = render_dg(source);
    assert!(pressure.iter().all(|p| p.is_finite()));
    assert!(pressure.iter().any(|p| p.abs() > 1e-15));
    assert!(
        dg_valve(&rendered)
            .material_relaxation()
            .unwrap()
            .memory_sqrt_j()[0]
            .iter()
            .any(|z| *z != 0.0)
    );
    for (c33, time) in [(8e8, 0.001), (3e8, 0.003)] {
        let changed = build_dg(&text, 37, c33, time);
        assert_eq!(
            dg_valve(changed.renderer()).spec(),
            &original_spec,
            "geometry, equilibrium stiffness, inertia, port maps and gas stay fixed"
        );
        let (other, changed) = render_dg(changed);
        assert!(
            pressure
                .iter()
                .zip(&other)
                .any(|(a, b)| (a - b).abs() > 1e-15)
        );
        assert_ne!(
            rendered.system().load().state(),
            changed.system().load().state(),
            "the original bore must feel the changed material response"
        );
    }
    let (partitioned, other) = render_dg(build_dg(&text, 1, 3e8, 0.001));
    assert_eq!(pressure, partitioned);
    assert_eq!(
        rendered.system().load().state(),
        other.system().load().state()
    );
    assert_eq!(
        dg_valve(&rendered).modal_states(),
        dg_valve(&other).modal_states()
    );
    assert_eq!(
        dg_valve(&rendered)
            .material_relaxation()
            .unwrap()
            .memory_sqrt_j(),
        dg_valve(&other)
            .material_relaxation()
            .unwrap()
            .memory_sqrt_j()
    );
    assert_eq!(
        encode_pcm16_wav(&pressure, 96000, 1.0).unwrap(),
        encode_pcm16_wav(&partitioned, 96000, 1.0).unwrap()
    );
}

#[test]
fn g3_same_tensor_source_uses_original_two_sided_network_and_external_pressure_owner() {
    let text = tensor_network();
    let source = build_network(&text, 37, 3e8, 0.001);
    assert_eq!(
        valve(&source)
            .material_relaxation()
            .unwrap()
            .tensor_material_points()
            .len(),
        12
    );
    let (pressure, rendered) = render(source);
    let (other, changed) = render(build_network(&text, 37, 8e8, 0.001));
    assert!(pressure.iter().all(|p| p.is_finite()));
    assert!(pressure.iter().any(|p| p.abs() > 1e-15));
    assert!(
        pressure
            .iter()
            .zip(&other)
            .any(|(a, b)| (a - b).abs() > 1e-15)
    );
    assert_ne!(
        rendered
            .source()
            .system()
            .node_frame(ValveSide::Downstream, 0),
        changed
            .source()
            .system()
            .node_frame(ValveSide::Downstream, 0)
    );
    let (partitioned, _) = render(build_network(&text, 1, 3e8, 0.001));
    assert_eq!(pressure, partitioned);
    assert_eq!(
        encode_pcm16_wav(&pressure, 96000, 1.0).unwrap(),
        encode_pcm16_wav(&partitioned, 96000, 1.0).unwrap()
    );
}

#[test]
fn g0_tensor_native_sources_keep_complete_72_state_budget_and_refuse_missing_domains() {
    let text = tensor_source();
    // The accepted 72-state boundary is exercised by both render tests above.
    for refused in [
        text.replace("1 12 72 16777216", "1 12 71 16777216"),
        text.replace("1 12 72 16777216", "1 11 72 16777216"),
        text.replace("1 12 72 16777216", "1 12 72 1"),
        text.replace("isotropic 293.15 100", "isotropic 350 100"),
    ] {
        assert!(
            MaterialModalDgValvePerformance::from_bytes_with_sources(
                refused.as_bytes(),
                37,
                &CancelGate::new_clock_free(),
                None,
                |_| Ok(tensor_pack(3e8, 0.001, "", 100000.0))
            )
            .is_err()
        );
    }
    for (missing, upper) in [
        ("branch_33", 100000.0),
        ("zero_modulus", 100000.0),
        ("", 100.0),
    ] {
        assert!(
            MaterialModalDgValvePerformance::from_bytes_with_sources(
                text.as_bytes(),
                37,
                &CancelGate::new_clock_free(),
                None,
                |_| Ok(tensor_pack(3e8, 0.001, missing, upper))
            )
            .is_err()
        );
        assert!(
            MaterialModalValvePerformance::from_bytes_with_sources(
                tensor_network().as_bytes(),
                37,
                &CancelGate::new_clock_free(),
                None,
                |_| Ok(tensor_pack(3e8, 0.001, missing, upper))
            )
            .is_err()
        );
    }
}
