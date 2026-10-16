use super::*;

fn declaration(
    tol: f64,
    degree: usize,
    plane: ReferencePlane,
    modes: ModeContent,
) -> BakeDeclaration {
    BakeDeclaration::new(
        ContentHash([1; 32]),
        (10.0, 20.0),
        vec![10.0, 20.0],
        tol,
        plane,
        modes,
        degree,
        Medium::air(),
    )
    .unwrap()
}

#[test]
fn declaration_records_every_physical_choice_before_solving() {
    let original = declaration(1e-3, 2, ReferencePlane::WholeSurface, ModeContent::Plane);
    let make = |geometry, band, frequencies, tol, plane, modes, degree, medium| {
        BakeDeclaration::new(
            geometry,
            band,
            frequencies,
            tol,
            plane,
            modes,
            degree,
            medium,
        )
        .unwrap()
        .hash()
    };
    let g = ContentHash([1; 32]);
    let b = (10.0, 20.0);
    let f = vec![10.0, 20.0];
    let p = ReferencePlane::WholeSurface;
    let m = Medium::air();
    for changed in [
        make(
            ContentHash([2; 32]),
            b,
            f.clone(),
            1e-3,
            p,
            ModeContent::Plane,
            2,
            m,
        ),
        make(
            g,
            (10.0, 21.0),
            vec![10.0, 21.0],
            1e-3,
            p,
            ModeContent::Plane,
            2,
            m,
        ),
        make(
            g,
            b,
            vec![10.0, 15.0, 20.0],
            1e-3,
            p,
            ModeContent::Plane,
            2,
            m,
        ),
        make(g, b, f.clone(), 2e-3, p, ModeContent::Plane, 2, m),
        make(
            g,
            b,
            f.clone(),
            1e-3,
            ReferencePlane::DuctMouth {
                origin_m: [0.0; 3],
                outward_normal: [1.0, 0.0, 0.0],
            },
            ModeContent::Plane,
            2,
            m,
        ),
        make(
            g,
            b,
            f.clone(),
            1e-3,
            p,
            ModeContent::Azimuthal {
                modes: vec![(0, 0), (1, 0)],
            },
            2,
            m,
        ),
        make(g, b, f.clone(), 1e-3, p, ModeContent::Plane, 3, m),
        make(
            g,
            b,
            f.clone(),
            1e-3,
            p,
            ModeContent::Plane,
            2,
            Medium {
                density: m.density * 1.01,
                ..m
            },
        ),
        make(
            g,
            b,
            f,
            1e-3,
            p,
            ModeContent::Plane,
            2,
            Medium {
                sound_speed: m.sound_speed * 1.01,
                ..m
            },
        ),
    ] {
        assert_ne!(original.hash(), changed);
    }
}

#[test]
fn declaration_cap_and_invalid_inputs_refuse() {
    let g = ContentHash([1; 32]);
    let p = ReferencePlane::WholeSurface;
    let m = Medium::air();
    assert!(
        BakeDeclaration::new(
            g,
            (1.0, 4096.0),
            (1..=4096).map(f64::from).collect(),
            1e-3,
            p,
            ModeContent::Plane,
            MAX_SH_DEGREE,
            m
        )
        .is_ok()
    );
    assert!(
        BakeDeclaration::new(
            g,
            (1.0, 4097.0),
            (1..=4097).map(f64::from).collect(),
            1e-3,
            p,
            ModeContent::Plane,
            2,
            m
        )
        .is_err()
    );
    assert!(
        BakeDeclaration::new(
            g,
            (10.0, 20.0),
            vec![10.0, 20.0],
            1e-3,
            p,
            ModeContent::Plane,
            MAX_SH_DEGREE + 1,
            m
        )
        .is_err()
    );
    for f in [
        vec![],
        vec![10.0],
        vec![11.0, 20.0],
        vec![10.0, 19.0],
        vec![10.0, 15.0, 15.0, 20.0],
        vec![10.0, f64::NAN, 20.0],
    ] {
        assert!(
            BakeDeclaration::new(g, (10.0, 20.0), f, 1e-3, p, ModeContent::Plane, 2, m).is_err()
        );
    }
    for tol in [0.0, -1.0, 1.0, f64::NAN, f64::INFINITY] {
        assert!(
            BakeDeclaration::new(
                g,
                (10.0, 20.0),
                vec![10.0, 20.0],
                tol,
                p,
                ModeContent::Plane,
                2,
                m
            )
            .is_err()
        );
    }
    assert!(
        BakeDeclaration::new(
            ContentHash([0; 32]),
            (10.0, 20.0),
            vec![10.0, 20.0],
            1e-3,
            p,
            ModeContent::Plane,
            2,
            m
        )
        .is_err()
    );
}

#[test]
fn dense_cap_is_inclusive_without_allocating_a_dense_system() {
    assert_eq!(dense_status(1000.0, MAX_DENSE_PANELS), None);
    assert_eq!(
        dense_status(1000.0, MAX_DENSE_PANELS + 1),
        Some(BakeStatus::Blocked {
            at_hz: 1000.0,
            needed_panels: MAX_DENSE_PANELS + 1,
            route: "M2.5-W5"
        })
    );
}

#[test]
fn complex_refinement_norm_sees_phase_and_quiet_angular_content() {
    let a = [C64::new(1.0, 0.0), C64::new(0.0, 1.0)];
    let b = [C64::new(0.0, 1.0), C64::new(1.0, 0.0)];
    assert!(relative_change(&a, &b) > 1.0);
    assert_eq!(relative_change(&[C64::ZERO], &[C64::ZERO]), 0.0);
    let quiet = [C64::new(1e-200, 0.0)];
    assert_eq!(relative_change(&quiet, &[C64::new(2e-200, 0.0)]), 0.5);
}

#[test]
fn receipt_strings_escape_run_provenance() {
    assert_eq!(
        json_string("host\"\\\n\t\0"),
        "\"host\\\"\\\\\\n\\t\\u0000\""
    );
}

#[test]
fn driven_plane_and_normal_are_checked_before_solving() {
    let surface =
        SpherePanels::from_triangles(vec![[[0.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]])
            .unwrap();
    let mesh = BakeMesh {
        surface,
        driven: vec![true],
    };
    let plane = ReferencePlane::DuctMouth {
        origin_m: [0.0; 3],
        outward_normal: [1.0, 0.0, 0.0],
    };
    assert!(admit_reference_plane(&mesh, plane).is_ok());
    for changed in [
        ReferencePlane::DuctMouth {
            origin_m: [0.01, 0.0, 0.0],
            outward_normal: [1.0, 0.0, 0.0],
        },
        ReferencePlane::DuctMouth {
            origin_m: [0.0; 3],
            outward_normal: [-1.0, 0.0, 0.0],
        },
    ] {
        assert!(admit_reference_plane(&mesh, changed).is_err());
    }
    let wrong = BakeMesh {
        surface: mesh.surface,
        driven: vec![],
    };
    assert!(admit_reference_plane(&wrong, plane).is_err());
}
