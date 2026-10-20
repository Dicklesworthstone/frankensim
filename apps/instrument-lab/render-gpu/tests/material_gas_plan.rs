//! Host-only admission of a GPU plan: no CUDA feature, driver or device needed.

#[path = "../src/plan.rs"]
#[allow(dead_code)]
mod plan;

use fs_vkplate::{GasState, Instrument, kit, standard_gas};

fn drum(gas: GasState) -> Instrument {
    let mut spec = kit::kick();
    spec.grid *= 8.0;
    kit::drum_with_gas(&spec, &[], gas).unwrap().instrument
}

#[test]
fn shared_plan_refuses_a_batch_member_with_another_atmosphere() {
    let gas = standard_gas();
    let mut first = drum(gas);
    let image = first.image();
    let plan = plan::Plan::new(&image, &first, 1, 64, true).unwrap();
    assert!(plan.state(&image).is_ok());
    assert!(plan.state(&first.clone().image()).is_ok());

    let compressed = GasState::try_new_moist_air(293.15, 2.0 * gas.pressure, 0.0).unwrap();
    // A coherent state can also share rho and c while carrying different
    // transport data. Admission retains the whole material snapshot.
    let mut other_transport = gas;
    other_transport.dynamic_viscosity *= 2.0;
    other_transport.thermal_conductivity *= 2.0;
    other_transport.validate_coherence().unwrap();
    assert_eq!(other_transport.density, gas.density);
    assert_eq!(other_transport.sound_speed, gas.sound_speed);

    for other in [compressed, other_transport] {
        let mut member = drum(other);
        let other_image = member.image();
        assert!(matches!(plan.state(&other_image), Err(error) if error.contains("atmosphere")));
        // Export also refuses a mismatched image/owner pair before flattening.
        assert!(plan::Plan::new(&other_image, &first, 1, 64, true).is_err());
    }
    assert!(
        plan.state(&image).is_ok(),
        "refusal must not alter the shared plan"
    );
}

fn viscous_material(name: &str) -> fs_vkplate::materials::KitMaterials {
    use fs_matdb::{
        ClaimSet, InterpolationPolicy, MaterialStateId, NormalizedMaterialCardPack, NormalizedPack,
        PropertyClaim, PropertyKey, PropertyValue, Provenance, UncertaintyModel,
    };
    use fs_qty::{
        Density, Dims, DynViscosity, Pressure, QuantitySpec,
        semantic::{FrequencyConvention, QuantityKind, SemanticType, ValueForm},
    };
    let mut claims = ClaimSet::new();
    let temperature = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    let frequency = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::Frequency(FrequencyConvention::Cyclic),
        ValueForm::Static,
    ));
    for (key, value, dims) in [
        ("density", 1390.0, Density::DIMS),
        ("young_modulus", 4e9, Pressure::DIMS),
        ("poisson_ratio", 0.3, Dims::NONE),
        ("plane_stress_viscosity_11", 2e5, DynViscosity::DIMS),
        ("plane_stress_viscosity_22", 1e5, DynViscosity::DIMS),
        ("plane_stress_viscosity_66", 5e4, DynViscosity::DIMS),
        ("plane_stress_viscosity_12", 2e4, DynViscosity::DIMS),
        ("plane_stress_viscosity_16", 0.0, DynViscosity::DIMS),
        ("plane_stress_viscosity_26", 0.0, DynViscosity::DIMS),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(key, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: fs_evidence::ValidityDomain::unconstrained()
                    .with_quantity("T", temperature, 273.15, 323.15)
                    .with_quantity("frequency", frequency, 0.0, 64000.0),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: format!("manufactured GPU host {name}/{key}"),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let normalized = NormalizedPack::new(
        name,
        "manufactured-gpu-host",
        fs_blake3::hash_domain("gpu-host-test", name.as_bytes()),
        "CC0-1.0 manufactured test",
        claims,
        vec![],
        vec![],
    )
    .unwrap();
    let pack = NormalizedMaterialCardPack::new(
        MaterialStateId {
            chemistry: name.into(),
            phase: "solid".into(),
            process: "manufactured numerical test".into(),
            revision: 0,
        },
        normalized,
    )
    .unwrap()
    .to_bytes();
    fs_vkplate::materials::KitMaterials::from_packs_with_bending_loss(
        &pack,
        &[],
        &pack,
        "T",
        293.15,
        9e-6,
        0.005,
        0.05,
        "frequency",
        1000.0,
        [0.3, -0.4],
    )
    .unwrap()
}

#[test]
fn g1_gpu_plan_keeps_the_real_viscous_force_and_nonnegative_loss_operators() {
    let build = |name| {
        let mut spec = kit::kick();
        spec.grid *= 8.0;
        kit::drum_with_materials(&spec, &[], standard_gas(), &viscous_material(name))
            .unwrap()
            .instrument
    };
    let mut ins = build("source A");
    let image = ins.image();
    let plan = plan::Plan::new(&image, &ins, 1, 64, true).unwrap();
    let state = plan.state(&image).unwrap();
    assert_eq!(state.viscous_dissipated, vec![0.0; image.bodies.len()]);
    assert_eq!(plan.visc_ptr.len(), plan.nt + 1);
    assert_eq!(plan.vroot_ptr.len(), 3 * plan.nt + 1);
    let x: Vec<f64> = (0..plan.nt)
        .map(|i| ((i as f64) * 0.53).sin() * 1e-5)
        .collect();
    let mut global_power = 0.0;
    let mut native_power = 0.0;
    for (bi, body) in image.bodies.iter().enumerate() {
        let off = plan.body_off[bi];
        if body.viscosity.is_none() {
            continue;
        }
        for i in 0..body.n {
            let actual = (plan.visc_ptr[off + i] as usize..plan.visc_ptr[off + i + 1] as usize)
                .map(|q| plan.visc_val[q] * x[plan.visc_col[q] as usize])
                .sum::<f64>();
            let expected = (body.viscous.ptr[i]..body.viscous.ptr[i + 1])
                .map(|q| body.viscous.val[q] * x[off + body.viscous.col[q]])
                .sum::<f64>();
            assert_eq!(actual, expected);
        }
        for r in 0..3 * body.n {
            let actual = (plan.vroot_ptr[3 * off + r] as usize
                ..plan.vroot_ptr[3 * off + r + 1] as usize)
                .map(|q| plan.vroot_val[q] * x[plan.vroot_col[q] as usize])
                .sum::<f64>();
            let expected = (body.viscous_root.ptr[r]..body.viscous_root.ptr[r + 1])
                .map(|q| body.viscous_root.val[q] * x[off + body.viscous_root.col[q]])
                .sum::<f64>();
            assert_eq!(actual, expected);
            global_power += actual * actual;
            native_power += expected * expected;
        }
    }
    assert!(global_power > 0.0);
    assert_eq!(global_power, native_power);

    // Equal coefficient numbers from another source are still different data.
    let mut other = build("source B");
    let other_image = other.image();
    assert!(plan.state(&other_image).is_err());
    let mut forged = image.clone();
    forged.bodies[0].viscosity = other_image.bodies[0].viscosity.clone();
    assert!(plan.state(&forged).is_err());
    assert!(plan::Plan::new(&forged, &ins, 1, 64, true).is_err());
    forged = image.clone();
    forged.rate *= 4.0;
    assert!(plan.state(&forged).is_err());
    ins.rate *= 4.0;
    assert!(
        plan::Plan::new(&ins.image(), &ins, 1, 64, true).is_err(),
        "GPU admission must check the actual source frequency band"
    );
}
