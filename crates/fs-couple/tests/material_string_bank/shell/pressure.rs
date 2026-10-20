//! G3: distinct fluid boundaries on one sourced physical shell and mode basis.
use super::*;
use fs_couple::shear_shell::pressure::ShellPressureSurfaces;
use fs_plate::shell::pressure::{OrientedShellFacet, ShellFace};

fn faces(range: core::ops::Range<usize>, side: ShellFace) -> Vec<OrientedShellFacet> {
    range
        .map(|triangle| OrientedShellFacet { triangle, side })
        .collect()
}
fn build(cx: &Cx<'_>, gate: &CancelGate, body: SharedMaterialShellBody) -> ReciprocalStringBank {
    let mut options = options(0);
    options.network.max_connections = 16;
    options.network.max_setup_terms = 131_072;
    let gas = gas(101_325.0);
    let cavity = fs_couple::bowed_string::runtime::reciprocal::plate::MaterialBodyCavity {
        cavity: fs_couple::cavity_modal::HelmholtzModalCavity::new(
            fs_scenario::HelmholtzCavity {
                volume_m3: 0.01,
                neck_radius_m: 0.01,
                neck_length_m: 0.02,
                inner_end: fs_scenario::NeckEnd::Flanged,
                outer_end: fs_scenario::NeckEnd::Flanged,
            },
            gas,
            fs_couple::cavity_modal::CavityNeckLoss::Resistance(100.0),
        )
        .unwrap(),
        modal_budget: ModalAcousticTimeBudget::audible_reference(),
        listener_m: 1.0,
    };
    ReciprocalStringBank::from_material_shear_shell(
        cx,
        vec![member(7800.0).source, layered().source],
        gas,
        body,
        Some(cavity),
        law(),
        options,
        gate,
    )
    .unwrap()
}

#[test]
fn selected_cavity_row_drives_one_air_spring_independently_of_exterior_observation() {
    with_cx(|cx, gate| {
        let body = shell(700.0, 0.02, 0.17, gate);
        let surfaces = ShellPressureSurfaces {
            exterior: faces(0..9, ShellFace::Positive),
            cavity: faces(9..18, ShellFace::Negative),
        };
        let mut exposed = body.clone();
        exposed.specimen = body
            .specimen
            .with_pressure_surfaces(surfaces.clone(), gate)
            .unwrap();
        let mut internal = body;
        internal.specimen = internal
            .specimen
            .with_pressure_surfaces(
                ShellPressureSurfaces {
                    exterior: vec![],
                    cavity: surfaces.cavity.clone(),
                },
                gate,
            )
            .unwrap();
        let mut a = build(cx, gate, exposed);
        let mut b = build(cx, gate, internal);
        let image = a.source().material_shell_image().unwrap();
        let quiet = b.source().material_shell_image().unwrap();
        assert_eq!(image.source().pressure_surfaces(), Some(&surfaces));
        assert_eq!(image.physical_modes(), quiet.physical_modes());
        assert_eq!(image.force_ports(), quiet.force_ports());
        assert_eq!(
            image.damping().matrix_per_s(),
            quiet.damping().matrix_per_s()
        );
        assert_eq!(
            image.cavity_volume_shapes_m2_per_sqrt_kg(),
            quiet.cavity_volume_shapes_m2_per_sqrt_kg()
        );
        assert_ne!(
            image.volume_shapes_m2_per_sqrt_kg(),
            image.cavity_volume_shapes_m2_per_sqrt_kg()
        );
        assert!(
            quiet
                .volume_shapes_m2_per_sqrt_kg()
                .iter()
                .all(|v| *v == 0.0)
        );
        let cavity = a.source().material_cavity().unwrap();
        let expected = cavity
            .cavity
            .connections(0, 1, image.cavity_volume_shapes_m2_per_sqrt_kg())
            .unwrap();
        assert_eq!(&a.source().body_connections()[..expected.len()], &expected);
        assert_eq!(a.source().body_connections(), b.source().body_connections());
        let weights = a.system().acceleration_pressure_weights().unwrap();
        for (i, &volume) in image.volume_shapes_m2_per_sqrt_kg().iter().enumerate() {
            let expected = gas(101_325.0).density * volume / core::f64::consts::TAU;
            assert!((weights[8 + i] - expected).abs() <= 1e-12 * expected.abs().max(1.0));
        }
        let zero = vec![0.0; a.system().mode_count()];
        let mut difference = 0.0_f64;
        for _ in 0..96 {
            let left = *a.step(&zero, drive()).unwrap();
            let right = *b.step(&zero, drive()).unwrap();
            difference =
                difference.max((left.observer_pressure_pa - right.observer_pressure_pa).abs());
            for frame in [left, right] {
                assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
            }
        }
        assert!(difference > 1e-15);
        for (left, right) in a.system().components().iter().zip(b.system().components()) {
            assert_eq!(
                left.states(),
                right.states(),
                "observer selection changed physical motion"
            );
        }
        assert!(
            b.system().components()[3]
                .states()
                .iter()
                .any(|state| state.velocity_m_sqrt_kg_per_s != 0.0),
            "cavity neck stopped responding"
        );
    });
}

#[test]
fn material_rebinding_keeps_authored_fluid_roles_and_physical_geometry() {
    with_cx(|_, gate| {
        let body = shell(700.0, 0.02, 0.17, gate);
        let surfaces = ShellPressureSurfaces {
            exterior: faces(0..5, ShellFace::Negative),
            cavity: faces(5..18, ShellFace::Positive),
        };
        let source = body
            .specimen
            .with_pressure_surfaces(surfaces.clone(), gate)
            .unwrap();
        let dense = source
            .with_region_material("wood", &shell_material(1400.0, 2000.0), gate)
            .unwrap();
        assert_eq!(dense.pressure_surfaces(), Some(&surfaces));
        assert_eq!(dense.mesh(), source.mesh());
        assert!((dense.mass_kg() / source.mass_kg() - 2.0).abs() < 1e-12);
        let mut duplicate = surfaces.clone();
        duplicate.cavity.push(duplicate.cavity[0]);
        assert!(source.with_pressure_surfaces(duplicate, gate).is_err());
        assert_eq!(source.pressure_surfaces(), Some(&surfaces));
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        assert!(source.with_pressure_surfaces(surfaces, &cancelled).is_err());
    });
}
