//! G3/G4: authored pressure roles reach the original coupled mechanics and PCM.
use super::*;

fn surface_profile(exterior: bool, cavity: bool) -> Profile {
    let mut block = format!(
        "pressure-surfaces {} {}\n",
        if exterior { 32 } else { 0 },
        if cavity { 32 } else { 0 }
    );
    for (role, present) in [("exterior", exterior), ("cavity", cavity)] {
        if present {
            for triangle in 0..32 {
                block.push_str(&format!("{role}-facet {triangle} positive\n"));
            }
        }
    }
    Profile::parse(
        shell_text()
            .replace("bridges 1\n", &format!("{block}bridges 1\n"))
            .as_bytes(),
    )
    .unwrap()
}

fn surfaced_ready(cx: &Cx<'_>, gate: &CancelGate, exterior: bool, cavity: bool) -> Ready {
    shell_spec()
        .compile_with_body(cx, gate, Some(surface_profile(exterior, cavity)), |name| {
            Ok(pack_with_shear(name, 1.0, 1.0, None))
        })
        .unwrap()
}

#[test]
fn exterior_selection_preserves_cavity_reaction_materials_and_mechanical_history() {
    scoped(|cx, gate| {
        let mut original = ready(cx, gate, 1.0, 1.0);
        let mut explicit = surfaced_ready(cx, gate, true, true);
        let mut cavity_only = surfaced_ready(cx, gate, false, true);
        let image = |r: &Ready| {
            r.performance
                .bank()
                .source()
                .bank()
                .material_shell_image()
                .unwrap()
                .clone()
        };
        let a = image(&explicit);
        let b = image(&cavity_only);
        assert_eq!(a.physical_modes(), b.physical_modes());
        assert_eq!(a.damping().matrix_per_s(), b.damping().matrix_per_s());
        assert_eq!(
            a.cavity_volume_shapes_m2_per_sqrt_kg(),
            b.cavity_volume_shapes_m2_per_sqrt_kg()
        );
        assert_eq!(
            a.volume_shapes_m2_per_sqrt_kg(),
            a.cavity_volume_shapes_m2_per_sqrt_kg()
        );
        assert!(a.volume_shapes_m2_per_sqrt_kg().iter().any(|v| *v != 0.0));
        assert!(b.volume_shapes_m2_per_sqrt_kg().iter().all(|v| *v == 0.0));
        let declaration = cavity_only
            .performance
            .bank()
            .source()
            .bank()
            .material_shell()
            .unwrap()
            .specimen
            .pressure_surfaces()
            .unwrap();
        assert!(declaration.exterior.is_empty());
        assert_eq!(declaration.cavity.len(), 32);
        assert_eq!(
            explicit
                .performance
                .bank()
                .source()
                .bank()
                .body_connections(),
            cavity_only
                .performance
                .bank()
                .source()
                .bank()
                .body_connections()
        );

        let p = pressure(&mut explicit, 7);
        assert_eq!(
            pressure(&mut original, 31),
            p,
            "explicit original faces must retain legacy pressure"
        );
        let q = pressure(&mut cavity_only, 96);
        assert!(
            q.iter().any(|v| *v != 0.0),
            "the physically driven cavity neck still radiates"
        );
        assert!(p.iter().zip(&q).any(|(a, b)| a != b));
        for (left, right) in explicit
            .performance
            .bank()
            .system()
            .components()
            .iter()
            .zip(cavity_only.performance.bank().system().components())
        {
            assert_eq!(
                left.states(),
                right.states(),
                "receiver selection must not rewrite mechanics"
            );
        }
        let mut replay = surfaced_ready(cx, gate, false, true);
        let expected =
            fs_couple::pcm_wav::encode_pcm16_wav(&q, replay.spec.rate, replay.spec.full_scale)
                .unwrap()
                .0;
        let (wav, _) = replay
            .write(std::io::Cursor::new(Vec::<u8>::new()), gate)
            .unwrap();
        assert_eq!(wav.into_inner(), expected);
    });
}

#[test]
fn explicit_no_fluid_faces_keep_the_complete_mechanical_body() {
    scoped(|cx, gate| {
        let mut quiet = surfaced_ready(cx, gate, false, false);
        let source = quiet.performance.bank().source().bank();
        let image = source.material_shell_image().unwrap();
        assert_eq!(image.physical_modes().len(), 2);
        assert!(source.material_shell().unwrap().specimen.mass_kg() > 0.0);
        assert!(
            image
                .cavity_volume_shapes_m2_per_sqrt_kg()
                .iter()
                .all(|v| *v == 0.0)
        );
        assert!(pressure(&mut quiet, 13).iter().all(|v| *v == 0.0));
        assert!(
            quiet.performance.bank().system().components()[1]
                .states()
                .iter()
                .any(|state| state.velocity_m_sqrt_kg_per_s != 0.0)
        );
        let frame = quiet.performance.bank().system().last_frame().unwrap();
        assert_eq!(frame.external_work_j, 0.0);
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
    });
}

#[test]
fn pressure_roles_are_explicit_bounded_and_independently_oriented() {
    for block in [
        "pressure-surfaces 33 0\n",
        "pressure-surfaces 1 0\nexterior-facet 32 positive\n",
        "pressure-surfaces 1 0\nexterior-facet 0 automatic\n",
        "pressure-surfaces 2 0\nexterior-facet 0 positive\nexterior-facet 0 negative\n",
        "pressure-surfaces 0 2\ncavity-facet 3 positive\ncavity-facet 3 positive\n",
        "pressure-surfaces 0 1\nexterior-facet 0 positive\n",
    ] {
        let text = shell_text().replace("bridges 1\n", &format!("{block}bridges 1\n"));
        assert!(Profile::parse(text.as_bytes()).is_err(), "accepted {block}");
    }
    let block = "pressure-surfaces 1 1\nexterior-facet 10 positive\ncavity-facet 10 negative\n";
    let profile = Profile::parse(
        shell_text()
            .replace("bridges 1\n", &format!("{block}bridges 1\n"))
            .as_bytes(),
    )
    .unwrap();
    profile.preflight(&shell_spec()).unwrap();
}
