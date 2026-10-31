//! Original winding shear data reaches the common mean-strain history and sound.
use super::*;

const GEOMETRY: &str = "geometric-stretching 0.04 64 1000000 1e-10 1e-8";
const CORE: &str = "stretching-memory 8 young_modulus 0";
const WRAP: &str = "winding-memory shear_modulus 1 wrap_delta wrap_time";

fn input_text() -> String {
    String::from_utf8(wound_spec().bytes)
        .unwrap()
        .replace("48000 256 0.01 128", "48000 17 0.01 16")
        .replace("surface 0.002", "surface 0.2")
        .replace(
            "cavity neck",
            &format!("{GEOMETRY}\n{CORE}\n{WRAP}\ncavity neck"),
        )
}

fn source_pack(name: &str, delta: f64, tau: f64, high: f64) -> NormalizedPack {
    let base = pack(
        name,
        if name == "wrap.fsmatpk" {
            8900.0
        } else {
            7800.0
        },
        200_000.0,
        (name == "wrap.fsmatpk").then_some(ETA),
    );
    if name != "wrap.fsmatpk" {
        return base;
    }
    let mut claims = base.claims().clone();
    for (key, value, dims) in [
        ("wrap_delta", delta, Pressure::DIMS),
        ("wrap_time", tau, Dims([0, 0, 1, 0, 0, 0])),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(key, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained()
                    .with("T", 280.0, 320.0)
                    .with("omega", 0.0, high),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured winding spectrum, not instrument calibration".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    NormalizedPack::new(
        base.pack_id(),
        "synthetic-winding-axial-memory-v1",
        base.source_artifact(),
        base.redistribution_terms(),
        claims,
        vec![],
        vec![],
    )
    .unwrap()
}

fn build(
    cx: &Cx<'_>,
    gate: &CancelGate,
    input: &str,
    delta: f64,
    tau: f64,
    high: f64,
) -> Result<Ready> {
    input::Spec::parse(input.as_bytes())?
        .compile(cx, gate, |name| Ok(source_pack(name, delta, tau, high)))
}

#[test]
fn original_wrap_spectrum_preserves_equilibrium_and_starts_at_actual_preload() {
    scoped(|cx, gate| {
        let ready = build(cx, gate, &input_text(), 80e9, 0.001, 200_000.0).unwrap();
        let elastic = build(cx, gate, &input_text(), 0.0, 0.001, 200_000.0).unwrap();
        let bank = ready.performance.bank();
        let image = bank
            .source()
            .bank()
            .geometric_stretching()
            .unwrap()
            .material_relaxation()
            .unwrap();
        assert_eq!(image.channels().len(), 1);
        let channel = &image.channels()[0];
        assert!(channel.winding);
        assert_eq!(
            (
                channel.string,
                channel.layer,
                channel.arm,
                channel.history_binding,
                channel.heat_channel
            ),
            (0, 1, 0, 0, 0)
        );
        let source = &image.sources()[0];
        assert!(source.core.spectra()[0].terms.is_empty());
        let wrap = source.winding.as_ref().unwrap();
        assert_eq!(wrap.shear_spectrum().terms, vec![(80e9, 0.001)]);
        assert_eq!(wrap.section(), specimen(&ready).section());
        // Delta G equals original G here, so the added rigidity equals exactly
        // one original spring. It does not alter that spring's equilibrium EA.
        assert!(
            (channel.relaxing_rigidity_n / wrap.section().winding_axial_stiffness_n() - 1.0).abs()
                < 2e-15
        );
        assert_eq!(
            specimen(&ready).section().axial_stiffness_n(),
            specimen(&elastic).section().axial_stiffness_n()
        );
        assert_eq!(
            bank.cards()[0].linear_density_kg_m,
            elastic.performance.bank().cards()[0].linear_density_kg_m
        );
        assert_eq!(
            bank.system().components()[0].states(),
            elastic.performance.bank().system().components()[0].states()
        );
        let strain = bank.stretching_observation(0).unwrap().additional_strain;
        assert!(strain > 0.0);
        assert_eq!(
            bank.system().nonlinear_material_memory(),
            vec![&[strain][..]]
        );
        assert!(
            elastic
                .performance
                .bank()
                .system()
                .nonlinear_material_memory()
                .is_empty()
        );
    });
}

#[test]
fn winding_relaxation_changes_shared_pressure_and_preserves_callback_replay() {
    scoped(|cx, gate| {
        let text = input_text();
        let mut slow = build(cx, gate, &text, 80e9, 0.003, 200_000.0).unwrap();
        let mut fast = build(cx, gate, &text, 80e9, 0.0002, 200_000.0).unwrap();
        assert_eq!(slow.spec.bytes, fast.spec.bytes);
        assert_eq!(
            slow.catalog["wire.fsmatpk"].pack.content_hash(),
            fast.catalog["wire.fsmatpk"].pack.content_hash()
        );
        assert_ne!(
            slow.catalog["wrap.fsmatpk"].pack.content_hash(),
            fast.catalog["wrap.fsmatpk"].pack.content_hash()
        );
        let initial = slow.performance.bank().system().nonlinear_material_memory()[0].to_vec();
        let expected = pressure(&mut slow, 7);
        let changed = pressure(&mut fast, 16);
        assert!(expected.iter().any(|p| p.abs() > 1e-12));
        assert!(
            expected
                .iter()
                .zip(&changed)
                .any(|(a, b)| (a - b).abs() > 1e-20)
        );
        assert_ne!(
            slow.performance.bank().system().nonlinear_material_memory()[0],
            &initial
        );
        let bank = slow.performance.bank();
        let source = specimen(&slow);
        let observed = bank.stretching_observation(0).unwrap();
        let memory = bank.system().nonlinear_material_memory();
        let image = bank
            .source()
            .bank()
            .geometric_stretching()
            .unwrap()
            .material_relaxation()
            .unwrap();
        let channel = &image.channels()[0];
        let elastic_strain = observed.additional_strain - memory[0][0];
        let expected_tension = source.string().tension_n
            + source.section().axial_stiffness_n() * observed.additional_strain
            + channel.relaxing_rigidity_n * elastic_strain;
        let expected_storage = 0.5
            * bank.cards()[0].length_m
            * (source.section().axial_stiffness_n() * observed.additional_strain.powi(2)
                + channel.relaxing_rigidity_n * elastic_strain.powi(2));
        assert!((observed.tension_n / expected_tension - 1.0).abs() < 2e-14);
        assert!((observed.stretching_energy_j / expected_storage - 1.0).abs() < 2e-14);
        let frame = slow.performance.bank().system().last_frame().unwrap();
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        assert!(frame.nonlinear.unwrap().nonlinear_dissipation_j > 0.0);
        assert_eq!(
            expected,
            pressure(
                &mut build(cx, gate, &text, 80e9, 0.003, 200_000.0).unwrap(),
                1
            )
        );
        let zero = text.replace(&format!("{CORE}\n{WRAP}\n"), "");
        assert_eq!(
            pressure(
                &mut build(cx, gate, &text, 0.0, 0.003, 200_000.0).unwrap(),
                1
            ),
            pressure(
                &mut build(cx, gate, &zero, 0.0, 0.003, 200_000.0).unwrap(),
                7
            )
        );
    });
}

#[test]
fn winding_memory_requires_matching_geometry_source_roles_and_complete_band() {
    let text = input_text();
    for changed in [
        text.replace(&format!("{GEOMETRY}\n"), ""),
        text.replace(&format!("{WRAP}\n"), ""),
        text.replace(CORE, "stretching-memory 0 young_modulus 0"),
        text.replace(
            WRAP,
            "winding-memory shear_modulus 2 wrap_delta wrap_time wrap_delta wrap_time",
        ),
    ] {
        assert!(input::Spec::parse(changed.as_bytes()).is_err());
    }
    scoped(|cx, gate| {
        for (delta, tau, high) in [
            (-1.0, 0.001, 200_000.0),
            (80e9, 0.0, 200_000.0),
            (80e9, 0.001, 1000.0),
        ] {
            assert!(build(cx, gate, &text, delta, tau, high).is_err());
        }
        for declaration in [
            "winding-memory shear_modulus 1 absent_delta wrap_time",
            "winding-memory young_modulus 1 wrap_delta wrap_time",
            "winding-memory shear_modulus 1 wrap_time wrap_delta",
        ] {
            assert!(
                build(
                    cx,
                    gate,
                    &text.replace(WRAP, declaration),
                    80e9,
                    0.001,
                    200_000.0
                )
                .is_err()
            );
        }
        let mut programmatic = input::Spec::parse(text.as_bytes()).unwrap();
        programmatic.stretching_memory.as_mut().unwrap().max_states = 0;
        let mut reads = 0;
        assert!(
            programmatic
                .compile(cx, gate, |name| {
                    reads += 1;
                    Ok(source_pack(name, 80e9, 0.001, 200_000.0))
                })
                .is_err()
        );
        assert_eq!(reads, 0);
        let mut cancelled = build(cx, gate, &text, 80e9, 0.001, 200_000.0).unwrap();
        let initial = cancelled
            .performance
            .bank()
            .system()
            .nonlinear_material_memory()[0]
            .to_vec();
        let stopped = CancelGate::new_clock_free();
        stopped.request();
        assert!(cancelled.write(Cursor::new(Vec::new()), &stopped).is_err());
        assert_eq!(
            cancelled
                .performance
                .bank()
                .system()
                .nonlinear_material_memory()[0],
            &initial
        );
        assert_eq!(cancelled.performance.samples_rendered(), 0);
    });
}

#[test]
fn prepared_winding_mean_strain_memory_preserves_physical_playback() {
    scoped(|cx, gate| {
        let text = input_text();
        let ordinary = build(cx, gate, &text, 80e9, 0.003, 200_000.0).unwrap();
        let mut spec = input::Spec::parse(text.as_bytes()).unwrap();
        spec.prepared_linear = Some(MEMORY_PREPARATION);
        let prepared = spec
            .compile(cx, gate, |name| {
                Ok(source_pack(name, 80e9, 0.003, 200_000.0))
            })
            .unwrap();
        assert_eq!(
            prepared.performance.bank().system().nonlinear_material_memory()[0].len(),
            1
        );
        assert_prepared_memory_parity(ordinary, prepared);
    });
}
