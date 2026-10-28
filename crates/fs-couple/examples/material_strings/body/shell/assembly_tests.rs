//! Complete component data -> common material mechanics -> microphone/PCM.
use super::*;

#[path = "assembly_glue_tests.rs"]
mod material_glue;

const ASSEMBLED: &str = include_str!("../../assembled.body");

fn ready(cx: &Cx<'_>, gate: &CancelGate, text: &str, right_loss: f64) -> Ready {
    let profile = Profile::parse(text.as_bytes()).unwrap();
    spec()
        .compile_with_body(cx, gate, Some(profile), |name| {
            Ok(pack(
                name,
                if name == "right-body.fsmatpk" {
                    right_loss
                } else {
                    1.0
                },
                None,
            ))
        })
        .unwrap()
}

#[test]
fn tied_components_share_modes_original_sources_cavity_and_block_independent_pcm() {
    scoped(|cx, gate| {
        let mut loads = BTreeMap::<String, usize>::new();
        let mut a = spec()
            .compile_with_body(
                cx,
                gate,
                Some(Profile::parse(ASSEMBLED.as_bytes()).unwrap()),
                |name| {
                    *loads.entry(name.into()).or_default() += 1;
                    Ok(pack(name, 1.0, None))
                },
            )
            .unwrap();
        assert_eq!(loads["body.fsmatpk"], 1);
        assert_eq!(loads["wire.fsmatpk"], 1);
        let source = a.performance.bank().source().bank();
        let body = source.material_shell_assembly().unwrap();
        let image = source.material_shell_assembly_image().unwrap();
        assert_eq!(body.specimen.components.len(), 2);
        assert_eq!(image.physical_modes().len(), 2);
        assert_eq!(image.force_ports().len(), 1);
        assert_eq!(image.damping().mode_count(), 3);
        assert!(image.tie_receipt().slave_dofs > 0);
        for mode in 0..3 {
            assert_eq!(image.physical_modes()[0][mode].len(), 54);
            assert_eq!(image.physical_modes()[1][mode].len(), 54);
            // The original rigid mortar enforces the shared physical seam.
            for (master, slave) in [(2, 0), (5, 3), (8, 6)] {
                for dof in 0..6 {
                    let x = image.physical_modes()[0][mode][6 * master + dof];
                    let y = image.physical_modes()[1][mode][6 * slave + dof];
                    assert!((x - y).abs() <= 1e-10 * x.abs().max(y.abs()).max(1.0));
                }
            }
        }
        assert!(
            image
                .volume_shapes_m2_per_sqrt_kg()
                .iter()
                .any(|v| *v != 0.0)
        );
        assert_eq!(
            image.volume_shapes_m2_per_sqrt_kg(),
            image.cavity_volume_shapes_m2_per_sqrt_kg()
        );
        assert_eq!(source.material_cavity().unwrap().cavity.gas(), source.gas());
        assert_eq!(a.body_profile.as_ref().unwrap().bytes, ASSEMBLED.as_bytes());
        let p = pressure(&mut a, 1);
        assert!(p.iter().any(|v| v.abs() > 1e-12));
        let frame = a.performance.bank().system().last_frame().unwrap();
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        let expected = fs_couple::pcm_wav::encode_pcm16_wav(&p, a.spec.rate, a.spec.full_scale)
            .unwrap()
            .0;
        let mut b = ready(cx, gate, ASSEMBLED, 1.0);
        assert_eq!(pressure(&mut b, 37), p);
        let mut c = ready(cx, gate, ASSEMBLED, 1.0);
        let (wav, report) = c
            .write(std::io::Cursor::new(Vec::<u8>::new()), gate)
            .unwrap();
        assert_eq!(wav.into_inner(), expected);
        assert_eq!(report.samples, a.spec.samples);
    });
}

#[test]
fn replacing_one_component_material_or_loss_reaches_shared_audio_without_string_retuning() {
    scoped(|cx, gate| {
        let different = ASSEMBLED.replace("region right @body", "region right right-body.fsmatpk");
        let mut original = ready(cx, gate, ASSEMBLED, 1.0);
        let mut replaced = ready(cx, gate, &different, 1.0);
        let mut lossier = ready(cx, gate, &different, 8.0);
        let modes = |r: &Ready| {
            r.performance.bank().system().components()[1]
                .modes()
                .to_vec()
        };
        assert_ne!(modes(&original), modes(&replaced));
        assert_eq!(modes(&replaced), modes(&lossier));
        let matrix = |r: &Ready| {
            r.performance
                .bank()
                .source()
                .bank()
                .material_shell_assembly_image()
                .unwrap()
                .damping()
                .matrix_per_s()
                .to_vec()
        };
        assert_ne!(matrix(&replaced), matrix(&lossier));
        for r in [&replaced, &lossier] {
            assert_eq!(
                original.performance.bank().cards()[0].fundamental_hz(),
                r.performance.bank().cards()[0].fundamental_hz()
            );
            assert_eq!(
                original.catalog["wire.fsmatpk"].pack.content_hash(),
                r.catalog["wire.fsmatpk"].pack.content_hash()
            );
        }
        let p = pressure(&mut original, 37);
        let q = pressure(&mut replaced, 128);
        let r = pressure(&mut lossier, 1);
        assert!(p.iter().zip(&q).any(|(a, b)| (a - b).abs() > 1e-12));
        assert!(q.iter().zip(&r).any(|(a, b)| (a - b).abs() > 1e-12));
        let encode = |p: &[f64]| {
            fs_couple::pcm_wav::encode_pcm16_wav(p, original.spec.rate, original.spec.full_scale)
                .unwrap()
                .0
        };
        assert_ne!(encode(&p), encode(&q));
        assert_ne!(encode(&q), encode(&r));
    });
}

#[test]
fn component_identity_fluid_roles_tie_addresses_and_bridge_extent_refuse_before_pack_io() {
    for (from, to) in [
        ("component right\n", "component left\n"),
        ("edge-tie left right", "edge-tie left unknown"),
        ("master-node 2", "master-node 90"),
        ("slave-node 3", "slave-node 0"),
        ("midsurface full", "compliant full"),
        ("bridge right 4", "bridge missing 4"),
        (
            "tie-budget 256 128 1e-10 1e-8",
            "tie-budget 256 128 NaN 1e-8",
        ),
    ] {
        assert!(
            Profile::parse(ASSEMBLED.replace(from, to).as_bytes()).is_err(),
            "accepted {to}"
        );
    }
    let start = ASSEMBLED.find("pressure-surfaces").unwrap();
    let end = start + ASSEMBLED[start..].find("end-component").unwrap();
    let missing_fluid_role = format!("{}{}", &ASSEMBLED[..start], &ASSEMBLED[end..]);
    assert!(Profile::parse(missing_fluid_role.as_bytes()).is_err());
    scoped(|cx, gate| {
        let too_many_bridges = ASSEMBLED.replace(
            "bridges 1\nbridge right 4 0 0 1",
            "bridges 2\nbridge left 4 0 0 1\nbridge right 4 0 0 1",
        );
        let mut calls = 0;
        assert!(
            spec()
                .compile_with_body(
                    cx,
                    gate,
                    Some(Profile::parse(too_many_bridges.as_bytes()).unwrap()),
                    |name| {
                        calls += 1;
                        Ok(pack(name, 1.0, None))
                    }
                )
                .is_err()
        );
        assert_eq!(calls, 0);
    });
}

/// The window's mode count is only known once the whole window is solved, so a request beyond it refuses naming the count
/// it holds. Cut between the fixture's second and third modes, the window holds two, and asking for three says so.
#[test]
fn a_request_beyond_the_assembled_window_names_the_modes_it_holds() {
    scoped(|cx, gate| {
        let lambda = ready(cx, gate, ASSEMBLED, 1.0)
            .performance
            .bank()
            .source()
            .bank()
            .material_shell_assembly_image()
            .unwrap()
            .eigenvalue_intervals()
            .to_vec();
        assert_eq!(lambda.len(), 3);
        assert!(
            lambda[1].1 < lambda[2].0,
            "the second and third modes overlap: {lambda:?}"
        );
        let cut = (lambda[1].1 * lambda[2].0).sqrt();
        let narrow = ASSEMBLED.replacen("window 1 100000000\n", &format!("window 1 {cut}\n"), 1);
        let refusal = spec()
            .compile_with_body(
                cx,
                gate,
                Some(Profile::parse(narrow.as_bytes()).unwrap()),
                |name| Ok(pack(name, 1.0, None)),
            )
            .err()
            .expect("a two-mode window cannot hold three modes")
            .to_string();
        assert!(refusal.contains("it holds 2, 3 requested"), "{refusal}");
    });
}

/// The body pack with its constant-loss anchor: `constant_loss_reference_omega` [rad/s], the angular frequency at which
/// its Kelvin-Voigt viscosities are its loss factors. Every other claim is the synthetic pack's own.
fn anchored(plain: &NormalizedPack, reference: f64) -> NormalizedPack {
    use fs_couple::shear_shell::{
        CONSTANT_LOSS_REFERENCE_OMEGA, CONSTANT_LOSS_REFERENCE_OMEGA_DIMS,
    };
    use fs_evidence::ValidityDomain;
    use fs_matdb::{
        InterpolationPolicy, ObservationDataset, PropertyClaim, PropertyKey, PropertyValue,
        Provenance, UncertaintyModel,
    };
    let artifact = plain.source_artifact();
    let provenance = Provenance {
        source: "synthetic constant-loss anchor".into(),
        license: "CC0-1.0".into(),
        artifact: Some(artifact),
    };
    let mut claims = plain.claims().clone();
    let observation = claims
        .register_observation(ObservationDataset {
            specimen: "synthetic".into(),
            method: "declared Kelvin-Voigt reference frequency".into(),
            artifact,
            caveats: "not measured".into(),
            provenance: provenance.clone(),
        })
        .unwrap();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(
                CONSTANT_LOSS_REFERENCE_OMEGA,
                CONSTANT_LOSS_REFERENCE_OMEGA_DIMS,
            ),
            value: PropertyValue::Scalar {
                value: reference,
                dims: CONSTANT_LOSS_REFERENCE_OMEGA_DIMS,
            },
            validity: ValidityDomain::unconstrained()
                .with("T", 280.0, 320.0)
                .with("omega", 0.0, 200_000.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![observation],
            provenance,
        })
        .unwrap();
    let pack = NormalizedPack::new(
        plain.pack_id(),
        plain.compiler(),
        artifact,
        plain.redistribution_terms(),
        claims,
        vec![],
        vec![],
    )
    .unwrap();
    NormalizedPack::from_bytes(&pack.to_bytes()).unwrap()
}

/// An anchored pack reads every retained mode's damping at the mode's own frequency (constant loss factors) instead of
/// letting the Kelvin-Voigt loss grow with frequency. Both components take the anchored body pack, so the modal damping
/// is the Kelvin-Voigt matrix congruent to s_a = sqrt(omega_ref G(omega_a) / omega_a), G the fitted five-branch
/// spectrum. The window is the spectrum's own band, 50 Hz to 10 kHz.
#[test]
fn an_anchored_loss_law_reads_every_modes_damping_at_its_own_frequency() {
    use fs_couple::shear_shell::constant_loss::{ConstantLossSpectrum, SPECTRUM_BAND_RAD_S};
    let reference = 2.0 * std::f64::consts::PI * 200.0;
    let (lo, hi) = SPECTRUM_BAND_RAD_S;
    let text = ASSEMBLED.replacen(
        "window 1 100000000\n",
        &format!("window {} {}\n", lo * lo, hi * hi),
        1,
    );
    scoped(|cx, gate| {
        let compile = |anchor: bool| {
            spec()
                .compile_with_body(
                    cx,
                    gate,
                    Some(Profile::parse(text.as_bytes()).unwrap()),
                    |name| {
                        let plain = pack(name, 1.0, None);
                        Ok(if anchor && name == "body.fsmatpk" {
                            anchored(&plain, reference)
                        } else {
                            plain
                        })
                    },
                )
                .unwrap()
        };
        let (mut kelvin_voigt, mut constant) = (compile(false), compile(true));
        let image = |r: &Ready| {
            let image = r
                .performance
                .bank()
                .source()
                .bank()
                .material_shell_assembly_image()
                .unwrap();
            (
                image.eigenvalue_intervals().to_vec(),
                image.damping().matrix_per_s().to_vec(),
            )
        };
        let ((intervals, viscous), (same, rescaled)) = (image(&kelvin_voigt), image(&constant));
        assert_eq!(intervals, same, "the loss law moved the conservative modes");
        let spectrum = ConstantLossSpectrum::fit();
        let squared = |lambda: f64| reference * spectrum.shape(lambda.sqrt()) / lambda.sqrt();
        // Across the flat spectrum s^2 falls with omega, so each mode's eigenvalue interval brackets its ratio. The
        // margin is a few dozen roundings (the spectrum, the square root, the rescale and the two components' sum).
        let rounding = 64.0 * f64::EPSILON;
        let n = intervals.len();
        for (a, &(low, high)) in intervals.iter().enumerate() {
            let ratio = rescaled[a * n + a] / viscous[a * n + a];
            assert!(
                squared(high) * (1.0 - rounding) <= ratio
                    && ratio <= squared(low) * (1.0 + rounding),
                "mode {a}: {ratio} outside [{}, {}]",
                squared(high),
                squared(low)
            );
        }
        let (p, q) = (
            pressure(&mut kelvin_voigt, 128),
            pressure(&mut constant, 128),
        );
        assert!(p.iter().zip(&q).any(|(a, b)| (a - b).abs() > 1e-12));
    });
}
