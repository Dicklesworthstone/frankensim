use super::*;
use fs_em::plane_wave::angular_spectrum::AplanaticPupil;

#[test]
fn high_na_pupil_uses_original_medium_and_six_component_focus() {
    // G1/G3: equal refractive index does not imply equal impedance. These
    // original epsilon/mu give k=1 and admittance=4 at omega=1.
    let medium = material(4.0, 0.25);
    let stack =
        PlanarDielectricLayers::try_new(&medium, DielectricResponse::Lossless(&medium), 1.0, &[])
            .unwrap();
    let grid = AngularSpectrumGrid {
        width: 17,
        height: 19,
        spacing_m: [2.0; 2],
    };
    let focus = [4.0, 6.0, 0.8];
    let pupil = AplanaticPupil::try_new(0.95, focus, 0.34).unwrap();
    let mut calls = 0;
    let beam = PlanarAngularSpectrum::prepare_aplanatic_pupil(
        &stack,
        grid,
        pupil,
        |point| {
            calls += 1;
            let [x, y, _] = point.direction;
            [C64::new(x, 0.0), C64::new(y, 0.0)]
        },
        limits(),
        || false,
    )
    .unwrap();
    assert!(calls > 20 && calls < grid.width * grid.height);
    assert!(std::ptr::eq(beam.source(), &stack));
    relative(beam.power().incident_w, 0.34);
    relative(beam.power().transmitted_w, 0.34);
    let field = beam.field(Side::Transmitted, focus[2], || false).unwrap();
    relative(field.normal_power_w(), 0.34);
    let e = field.electric_v_per_m()[3 * grid.width + 2];
    assert!(e[2].abs() > 0.03);
    complex_near(e[0], C64::ZERO, 2e-14);
    complex_near(e[1], C64::ZERO, 2e-14);
}

#[test]
fn high_na_film_diffraction_keeps_power_balance_and_preflight_limits() {
    let incident = material(1.0, 1.0);
    let film = material(2.25, 1.0);
    let stack = PlanarDielectricLayers::try_new(
        &incident,
        DielectricResponse::Lossless(&incident),
        1.0,
        &[DielectricLayer {
            response: DielectricResponse::Lossless(&film),
            thickness_m: 0.73,
        }],
    )
    .unwrap();
    let grid = AngularSpectrumGrid {
        width: 11,
        height: 13,
        spacing_m: [2.0; 2],
    };
    let pupil = AplanaticPupil::try_new(0.9, [0.0, 0.0, 2.0], 0.7).unwrap();
    let beam = PlanarAngularSpectrum::prepare_aplanatic_pupil(
        &stack,
        grid,
        pupil,
        |p| {
            let x = p.pupil_coordinate[0];
            [C64::new(1.0, 0.2 * x), C64::new(0.0, 0.4)]
        },
        limits(),
        || false,
    )
    .unwrap();
    relative(beam.power().incident_w, 0.7);
    assert!(beam.power().reflected_w > 0.0);
    near(beam.power().balance_defect_w(), 0.0, 2e-12);
    let field = beam.field(Side::Transmitted, 0.3, || false).unwrap();
    relative(field.normal_power_w(), beam.power().transmitted_w);

    let mut called = false;
    let mut too_small = limits();
    // The ordinary beam's query allowance does not include this constructor's
    // internally owned N complex-pair pupil samples. Refuse before callback.
    too_small.maximum_workspace_bytes = beam.workspace_bytes();
    assert!(
        PlanarAngularSpectrum::prepare_aplanatic_pupil(
            &stack,
            grid,
            pupil,
            |_| {
                called = true;
                [C64::ONE, C64::ZERO]
            },
            too_small,
            || false,
        )
        .is_err()
    );
    assert!(!called);
    assert!(matches!(
        PlanarAngularSpectrum::prepare_aplanatic_pupil(
            &stack,
            grid,
            pupil,
            |_| [C64::ONE, C64::ZERO],
            limits(),
            || true,
        ),
        Err(AngularSpectrumError::Cancelled)
    ));
}
