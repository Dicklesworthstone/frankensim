//! The shipped closed geometry selects the original shell/air coupling.
use super::*;

#[test]
fn closed_shell_source_retains_geometry_and_refuses_a_planar_substitute_before_pack_io() {
    let spec =
        input::Spec::parse(include_bytes!("../closed-shell-pluck.material-strings")).unwrap();
    let Some(input::BodyRadiation::ShellBem(options)) = spec.body_radiation else {
        panic!("closed shell source selected a different acoustic model");
    };
    assert_eq!(options.maximum_surface_panels, 80);
    let profile = body::Profile::parse(include_bytes!("../closed-shell.body")).unwrap();
    assert!(profile.has_shell_material_basis());
    profile.preflight(&spec).unwrap();
    scoped(|cx, gate| {
        for profile in [
            None,
            Some(body::Profile::parse(include_bytes!("../orthotropic.body")).unwrap()),
        ] {
            let mut reads = 0;
            assert!(
                spec.clone()
                    .compile_with_body::<NormalizedPack>(cx, gate, profile, |_| {
                        reads += 1;
                        Err("unexpected material read".into())
                    })
                    .is_err()
            );
            assert_eq!(reads, 0, "incompatible body must fail before material I/O");
        }
    });
}
