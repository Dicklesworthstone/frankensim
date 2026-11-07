use super::tests::{fixture, limits, source};
use super::*;
use fs_material::visco::temperature::TemperatureShiftModel;

fn shift() -> TemperatureShift {
    TemperatureShift::new(
        TemperatureShiftModel::WilliamsLandelFerry {
            c1: 2.0,
            c2_k: 20.0,
        },
        295.0,
        [280.0, 330.0],
    )
    .unwrap()
}

fn near(a: f64, b: f64, scale: f64) {
    assert!(
        (a - b).abs() <= 3e-11 * scale.max(1e-25),
        "{a:e} != {b:e}; scale={scale:e}"
    );
}

#[test]
fn fixed_temperature_scales_original_rotated_resultants_without_compounding() {
    let law = source([4.0, 9.0, 16.0, 1.0, 0.5, -0.75], (10.0, 1e6));
    let original = ShellMaterialViscosity::new(law, 0.005, 0.37).unwrap();
    assert!(original.temperature_shift_state().is_none());
    let (a, b, d) = original.section().matrices();
    for temperature in [285.0, 315.0] {
        let changed = original.at_temperature(&shift(), "T", temperature).unwrap();
        let state = changed.temperature_shift_state().unwrap();
        let factor = state.factor();
        assert_eq!(state.temperature_k(), temperature);
        assert!(std::ptr::eq(changed.source(), original.source()));
        assert_eq!(changed.source().source().state_coordinate("T"), Some(295.0));
        assert_eq!(changed.omega_band_rad_s(), (10.0 / factor, 1e6 / factor));
        assert_eq!(changed.thickness_m(), original.thickness_m());
        assert_eq!(changed.material_angle_rad(), original.material_angle_rad());
        let (aa, bb, dd) = changed.section().matrices();
        assert_eq!(bb, b);
        for i in 0..9 {
            near(
                aa[i],
                factor * a[i],
                factor * a.iter().map(|x| x.abs()).fold(0.0, f64::max),
            );
            near(
                dd[i],
                factor * d[i],
                factor * d.iter().map(|x| x.abs()).fold(0.0, f64::max),
            );
        }
        let repeated = changed.at_temperature(&shift(), "T", temperature).unwrap();
        assert_eq!(repeated.section().matrices(), changed.section().matrices());
        assert_eq!(repeated.omega_band_rad_s(), changed.omega_band_rad_s());
        let reference = changed.at_temperature(&shift(), "T", 295.0).unwrap();
        assert_eq!(reference.section().matrices(), original.section().matrices());
        assert_eq!(reference.omega_band_rad_s(), original.omega_band_rad_s());
    }
    let hot = original.at_temperature(&shift(), "T", 315.0).unwrap();
    assert_eq!(
        hot.at_temperature(&shift(), "T", 285.0)
            .unwrap()
            .section()
            .matrices(),
        original
            .at_temperature(&shift(), "T", 285.0)
            .unwrap()
            .section()
            .matrices()
    );
    assert!(original.at_temperature(&shift(), "missing", 315.0).is_err());
}

#[test]
fn complete_modal_damping_and_physical_channels_scale_at_fixed_temperature() {
    let (mesh, model, modes) = fixture();
    let facet = ShellMaterialViscosity::new(
        source([400.0, 900.0, 1600.0, 100.0, 50.0, -75.0], (0.0, 1e9)),
        0.005,
        0.31,
    )
    .unwrap();
    let gate = CancelGate::new_clock_free();
    let reference = ModalShellMaterialViscosity::project(
        &mesh,
        &[facet.clone()],
        &model,
        &modes,
        limits(),
        &gate,
    )
    .unwrap();
    for temperature in [285.0, 315.0] {
        let changed = facet.at_temperature(&shift(), "T", temperature).unwrap();
        let factor = changed.temperature_shift_state().unwrap().factor();
        let image = ModalShellMaterialViscosity::project(
            &mesh,
            &[changed],
            &model,
            &modes,
            limits(),
            &gate,
        )
        .unwrap();
        assert_eq!(
            image.angular_frequencies_rad_s(),
            reference.angular_frequencies_rad_s()
        );
        assert!(std::ptr::eq(image.facets()[0].source(), facet.source()));
        let expected = reference.modal().matrix_per_s();
        let scale = factor * expected.iter().map(|x| x.abs()).fold(0.0, f64::max);
        for (a, b) in image.modal().matrix_per_s().iter().zip(expected) {
            near(*a, factor * b, scale);
        }
        let velocity = [0.2, -0.1, 0.07, -0.03, 0.04, -0.02];
        let power = |law: &ModalPlateViscosity| {
            law.channel_shapes()
                .iter()
                .map(|row| {
                    let rate = row.iter().zip(velocity).map(|(c, v)| c * v).sum::<f64>();
                    law.reference_resistance_n_s_m() * rate * rate
                })
                .sum::<f64>()
        };
        near(
            power(image.modal()),
            factor * power(reference.modal()),
            factor * power(reference.modal()),
        );
    }
}

#[test]
fn actual_modal_support_uses_shifted_band_even_for_zero_tensor() {
    let (mesh, model, modes) = fixture();
    let gate = CancelGate::new_clock_free();
    let high = modes
        .iter()
        .map(|m| fs_math::det::sqrt(m.lambda))
        .fold(0.0, f64::max);
    for values in [[4.0, 9.0, 16.0, 1.0, 0.5, -0.75], [0.0; 6]] {
        // Warmth moves a limited reference band high enough to admit every
        // actual mode. Inspecting the unshifted source band would reject it.
        let narrow =
            ShellMaterialViscosity::new(source(values, (0.0, 0.5 * high)), 0.005, 0.2).unwrap();
        assert!(
            ModalShellMaterialViscosity::project(
                &mesh,
                &[narrow.clone()],
                &model,
                &modes,
                limits(),
                &gate,
            )
            .is_err()
        );
        let hot = narrow.at_temperature(&shift(), "T", 315.0).unwrap();
        let factor = hot.temperature_shift_state().unwrap().factor();
        let image =
            ModalShellMaterialViscosity::project(&mesh, &[hot], &model, &modes, limits(), &gate)
                .unwrap();
        assert_eq!(image.omega_band_rad_s(), (0.0, 0.5 * high / factor));
        // Cooling narrows a reference band that originally admitted the same
        // basis. The shifted snapshot must refuse those unsupported modes.
        let wide =
            ShellMaterialViscosity::new(source(values, (0.0, 1.1 * high)), 0.005, 0.2).unwrap();
        ModalShellMaterialViscosity::project(
            &mesh,
            &[wide.clone()],
            &model,
            &modes,
            limits(),
            &gate,
        )
        .unwrap();
        let cold = wide.at_temperature(&shift(), "T", 285.0).unwrap();
        assert!(
            ModalShellMaterialViscosity::project(&mesh, &[cold], &model, &modes, limits(), &gate)
                .is_err()
        );
    }
}

#[test]
fn temperature_preserves_weak_physical_shell_work_before_coefficient_rounding() {
    use fs_plate::PlateSection;
    use fs_plate::shell::viscosity::{ShellViscosityBudget, project_viscosity};
    use fs_plate::shell::{ShellSupport, assemble_shell};
    let mesh = ShellMesh::new(
        vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        vec![[0, 1, 2]],
    )
    .unwrap();
    let thickness = 2.0_f64.powi(-8);
    let elastic = PlateSection::isotropic(1e7, 0.25, thickness, 800.0).unwrap();
    let model = assemble_shell(&mesh, &elastic, &[], ShellSupport::Free).unwrap();
    let mode: Vec<_> = mesh
        .nodes
        .iter()
        .flat_map(|p| [-p[0], p[1], 0.0, 0.0, 0.0, 0.0])
        .collect();
    let a = 3.7;
    let b = a + a * 2.0_f64.powi(-48);
    let reference =
        ShellMaterialViscosity::new(source([a, b, 0.0, a, 0.0, 0.0], (0.0, 1e6)), thickness, 0.0)
            .unwrap();
    assert_eq!(
        reference.section().scaled(1.0).unwrap().matrices(),
        reference.section().matrices()
    );
    for factor in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::from_bits(1)] {
        assert!(reference.section().scaled(factor).is_err());
    }
    let power = |facet: &ShellMaterialViscosity| {
        project_viscosity(
            &mesh,
            &[facet.section().clone()],
            &model,
            &[mode.clone()],
            ShellViscosityBudget {
                max_elements: 1,
                max_full_dofs: 18,
                max_dofs: 18,
                max_modes: 1,
                max_work: 1_000_000,
            },
            || false,
        )
        .unwrap()[0]
    };
    for temperature in [315.0, 310.0, 285.0] {
        let shifted = reference
            .at_temperature(&shift(), "T", temperature)
            .unwrap();
        let factor = shifted.temperature_shift_state().unwrap().factor();
        let expected = factor * power(&reference);
        let actual = power(&shifted);
        assert!(
            (actual / expected - 1.0).abs() < 3e-12,
            "weak physical shell shift T={temperature}, aT={factor:e}, actual={actual:e}, expected={expected:e}, relative={:e}",
            actual / expected - 1.0
        );
    }
}
