use super::*;
use crate::dielectric::{BeerLambertAbsorption, CauchyIor, DielectricGlass, GlassProvenance};
use crate::polarization::Mueller;
use fs_geom::{Point3, Vec3};

fn frame(direction: Vec3, horizontal: Vec3) -> PolarizationFrame {
    PolarizationFrame::try_new(direction, horizontal).unwrap()
}

fn camera(frame: PolarizationFrame, power: Stokes) -> PolarizedBeam {
    PolarizedBeam::try_new(
        Point3::new(0.0, 0.0, 0.0),
        frame,
        550.0,
        DielectricGlass::new(
            CauchyIor::try_constant(1.33).unwrap(),
            BeerLambertAbsorption::CLEAR,
            GlassProvenance::Custom,
        ),
        power,
    )
    .unwrap()
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-13 * expected.abs().max(1.0),
        "{actual:.17e} != {expected:.17e}"
    );
}

#[test]
fn unpolarized_connections_recover_the_native_phase_without_prefix_or_pdf_factors() {
    let camera_frame = frame(Vec3::new(0.0, 0.0, -1.0), Vec3::new(0.0, 1.0, 0.0));
    for delta in [0.0, 0.03, 0.5, 0.85] {
        let phase = PhaseFunction::rayleigh_depolarized(delta).unwrap();
        for cosine in [-1.0_f64, -0.8, -0.3, 0.0, 0.4, 0.9, 1.0] {
            let source_frame = frame(
                Vec3::new(((1.0 - cosine) * (1.0 + cosine)).sqrt(), 0.0, cosine),
                Vec3::new(0.0, 1.0, 0.0),
            );
            for prefix_intensity in [1e-200, 1.0, 7.0, 1e200] {
                let factor = phase
                    .rayleigh_connection_factor(
                        camera(camera_frame, Stokes::unpolarized(prefix_intensity)),
                        source_frame,
                        Stokes::unpolarized(1.0),
                    )
                    .unwrap();
                close(factor, phase.eval_cos(cosine).unwrap());
            }
        }
    }
}

#[test]
fn reciprocal_connection_uses_negative_u_positive_v_and_physical_backscatter() {
    let horizontal = Vec3::new(1.0, 0.0, 0.0);
    let camera_frame = frame(Vec3::new(0.0, 0.0, -1.0), horizontal);
    let source_forward = frame(Vec3::new(0.0, 0.0, 1.0), horizontal);
    let source_backward = frame(Vec3::new(0.0, 0.0, -1.0), horizontal);
    let phase = PhaseFunction::Rayleigh;
    let bright = 2.0 * phase.eval_cos(1.0).unwrap();
    let source_u = Stokes::try_new(1.0, 0.0, 1.0, 0.0).unwrap();
    let source_v = Stokes::try_new(1.0, 0.0, 0.0, 1.0).unwrap();
    let source_q = Stokes::try_new(1.0, 1.0, 0.0, 0.0).unwrap();
    let cases = [
        (
            source_forward,
            source_q,
            Stokes {
                i: 1.0,
                q: -1.0,
                u: 0.0,
                v: 0.0,
            },
            0.0,
        ),
        (source_forward, source_u, source_u, 0.0),
        (
            source_forward,
            source_u,
            Stokes {
                i: 1.0,
                q: 0.0,
                u: -1.0,
                v: 0.0,
            },
            bright,
        ),
        (source_forward, source_v, source_v, bright),
        (
            source_forward,
            source_v,
            Stokes {
                i: 1.0,
                q: 0.0,
                u: 0.0,
                v: -1.0,
            },
            0.0,
        ),
        (source_backward, source_v, source_v, 0.0),
        (
            source_backward,
            source_v,
            Stokes {
                i: 1.0,
                q: 0.0,
                u: 0.0,
                v: -1.0,
            },
            bright,
        ),
    ];
    for (source_frame, source, prefix, expected) in cases {
        let factor = phase
            .rayleigh_connection_factor(camera(camera_frame, prefix), source_frame, source)
            .unwrap();
        close(factor, expected);
        if expected == 0.0 {
            assert_eq!(factor, 0.0);
        }
    }
}

#[test]
fn reciprocal_retarder_chain_agrees_with_forward_intensity_for_circular_input() {
    // Physically reverse the order of these reciprocal optical elements and
    // reverse each orientation in the right-handed reverse-travel frame.
    // No transpose or adjoint helper constructs this independent reference.
    let forward = Mueller::linear_polarizer(0.23)
        .compose(Mueller::linear_retarder(-0.61, 1.8))
        .compose(Mueller::linear_retarder(0.37, 0.9));
    let reverse = Mueller::linear_retarder(-0.37, 0.9)
        .compose(Mueller::linear_retarder(0.61, 1.8))
        .compose(Mueller::linear_polarizer(-0.23));
    let prefix = reverse.apply(Stokes::unpolarized(1.0));
    assert!(prefix.u.abs() > 0.01 && prefix.v.abs() > 0.01);
    let camera = camera(
        frame(Vec3::new(0.0, 0.0, -1.0), Vec3::new(1.0, 0.0, 0.0)),
        prefix,
    );
    let source_frame = frame(Vec3::new(0.0, 0.0, 1.0), Vec3::new(1.0, 0.0, 0.0));
    let phase = PhaseFunction::Rayleigh;
    for source in [
        Stokes::unpolarized(1.0),
        Stokes::try_new(1.0, 0.2, -0.3, 0.4).unwrap(),
        Stokes::try_new(1.0, 0.0, 0.0, 1.0).unwrap(),
        Stokes::try_new(1.0, 0.0, 0.0, -1.0).unwrap(),
    ] {
        let factor = phase
            .rayleigh_connection_factor(camera, source_frame, source)
            .unwrap();
        close(
            prefix.i * factor,
            phase.eval_cos(1.0).unwrap() * forward.apply(source).i,
        );
    }
}

#[test]
fn connection_is_frame_covariant_and_reciprocal_for_nonplanar_mixed_states() {
    let prefix = Stokes::try_new(3.0, 1.2, -0.9, 1.8).unwrap();
    let camera = camera(
        frame(Vec3::new(-0.3, 0.4, -0.8), Vec3::new(1.0, 2.0, 0.7)),
        prefix,
    );
    let source_frame = frame(Vec3::new(0.7, -0.2, 0.5), Vec3::new(-0.4, 0.6, 1.0));
    let source = Stokes::try_new(1.0, -0.3, 0.4, -0.5).unwrap();
    let changed_source_frame = frame(source_frame.direction(), Vec3::new(0.7, 1.0, -0.6));
    let changed_source = source_frame.express_in(source, changed_source_frame);
    let changed_camera = camera
        .with_reference_axis(Vec3::new(0.8, -0.6, 0.3))
        .unwrap();
    for delta in [0.0, 0.03, 0.5, 0.85] {
        let phase = PhaseFunction::rayleigh_depolarized(delta).unwrap();
        let expected = phase
            .rayleigh_connection_factor(camera, source_frame, source)
            .unwrap();
        close(
            phase
                .rayleigh_connection_factor(changed_camera, changed_source_frame, changed_source)
                .unwrap(),
            expected,
        );
        let reverse_camera = PolarizedBeam::try_new(
            camera.ray().origin,
            source_frame,
            camera.wavelength_nm(),
            camera.medium(),
            source,
        )
        .unwrap();
        let reverse_source = Stokes {
            i: 1.0,
            q: prefix.q / prefix.i,
            u: prefix.u / prefix.i,
            v: prefix.v / prefix.i,
        };
        close(
            phase
                .rayleigh_connection_factor(reverse_camera, camera.frame(), reverse_source)
                .unwrap(),
            expected,
        );
    }
}

#[test]
fn near_dipole_null_retains_positive_power_and_invalid_optics_refuse() {
    let camera = camera(
        frame(Vec3::new(0.0, 0.0, -1.0), Vec3::new(0.0, 1.0, 0.0)),
        Stokes::unpolarized(1.0),
    );
    let cosine = 1e-9;
    let source_frame = frame(Vec3::new(1.0, 0.0, cosine), Vec3::new(0.0, 1.0, 0.0));
    let source = Stokes::try_new(1.0, -1.0, 0.0, 0.0).unwrap();
    let actual = PhaseFunction::Rayleigh
        .rayleigh_connection_factor(camera, source_frame, source)
        .unwrap();
    let expected = 3.0 / (8.0 * core::f64::consts::PI) * cosine * cosine;
    assert!(actual > 0.0);
    assert!((actual / expected - 1.0).abs() < 3e-14);
    assert_eq!(
        PhaseFunction::Isotropic.rayleigh_connection_factor(camera, source_frame, source),
        Err(PolarizedScatteringError::UnsupportedPhase),
    );
    for invalid in [
        Stokes::unpolarized(0.0),
        Stokes::unpolarized(2.0),
        Stokes {
            i: 1.0,
            q: 1.01,
            u: 0.0,
            v: 0.0,
        },
        Stokes {
            i: 1.0,
            q: 0.0,
            u: 0.0,
            v: f64::NAN,
        },
    ] {
        assert!(
            PhaseFunction::Rayleigh
                .rayleigh_connection_factor(camera, source_frame, invalid)
                .is_err()
        );
    }
    let dark = PolarizedBeam::try_new(
        camera.ray().origin,
        camera.frame(),
        550.0,
        camera.medium(),
        Stokes::unpolarized(0.0),
    )
    .unwrap();
    assert!(
        PhaseFunction::Rayleigh
            .rayleigh_connection_factor(dark, source_frame, source)
            .is_err()
    );
}
