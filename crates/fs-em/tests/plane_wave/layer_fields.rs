use super::*;
use fs_em::plane_wave::{
    angular_spectrum::{
        AngularSpectrumError, AngularSpectrumGrid, AngularSpectrumLayerLimits,
        AngularSpectrumLimits, AngularSpectrumSide, PlanarAngularSpectrum,
    },
    layers::{DielectricLayer, PlanarDielectricLayers},
};
use fs_math::det;

fn limits() -> AngularSpectrumLimits {
    AngularSpectrumLimits {
        maximum_grid_points: 4096,
        maximum_workspace_bytes: 4 << 20,
        maximum_fft_work: 1 << 20,
    }
}
fn layer_limits() -> AngularSpectrumLayerLimits {
    AngularSpectrumLayerLimits {
        maximum_layer_mode_pairs: 1 << 16,
        maximum_workspace_bytes: 8 << 20,
        maximum_fft_work: 1 << 20,
    }
}
fn relative(a: f64, b: f64, tolerance: f64) {
    near(a, b, tolerance * a.abs().max(b.abs()).max(1e-30));
}

#[test]
fn internal_film_traces_match_both_original_external_fields_and_lossless_power() {
    let incident = material(1.0, 1.0);
    let film = material(3.0, 1.4);
    let stack = PlanarDielectricLayers::try_new(
        &incident,
        DielectricResponse::Lossless(&incident),
        1.0,
        &[DielectricLayer {
            response: DielectricResponse::Lossless(&film),
            thickness_m: 0.8,
        }],
    )
    .unwrap();
    let grid = AngularSpectrumGrid {
        width: 8,
        height: 4,
        spacing_m: [8.0, 9.0],
    };
    let mut input = vec![[C64::ZERO; 2]; grid.width * grid.height];
    input[3] = [C64::new(0.9, 0.2), C64::new(-0.3, 0.8)];
    input[19] = [C64::new(0.2, -0.7), C64::new(0.4, 0.1)];
    let beam = PlanarAngularSpectrum::prepare(&stack, grid, &input, limits(), || false).unwrap();
    let layer = beam.prepare_layer(0, layer_limits(), || false).unwrap();
    assert!(std::ptr::eq(layer.beam(), &beam));
    assert_eq!(layer.absorbed_power_w(), 0.0);
    for (z, side) in [
        (0.0, AngularSpectrumSide::TotalIncidentSide),
        (0.8, AngularSpectrumSide::Transmitted),
    ] {
        let internal = layer.field(z, || false).unwrap();
        let external = beam.field(side, 0.0, || false).unwrap();
        for i in 0..input.len() {
            let ei = internal.field().electric_v_per_m()[i];
            let hi = internal.field().magnetic_a_per_m()[i];
            let ee = external.electric_v_per_m()[i];
            let he = external.magnetic_a_per_m()[i];
            for c in 0..2 {
                complex_near(ei[c], ee[c], 4e-13);
                complex_near(hi[c], he[c], 4e-13);
            }
            complex_near(ei[2].scale(3.0), ee[2], 4e-13);
            complex_near(hi[2].scale(1.4), he[2], 4e-13);
        }
        relative(
            internal.field().normal_power_w(),
            external.normal_power_w(),
            5e-12,
        );
    }
    for z in [0.0, 0.17, 0.53, 0.8] {
        let field = layer.field(z, || false).unwrap();
        assert!(
            field
                .absorbed_power_density_w_m3()
                .iter()
                .all(|v| *v == 0.0)
        );
        relative(
            field.field().normal_power_w(),
            beam.power().transmitted_w,
            5e-12,
        );
    }
}

#[test]
fn diffracting_beam_volume_heat_agrees_with_original_positive_film_loss_and_poynting() {
    let incident = material(E, M);
    let lossy = finite_layers::debye(2.0 * E, 1.2 * E, 1.0, 1.0);
    let thickness = 0.7;
    let stack = PlanarDielectricLayers::try_new(
        &incident,
        DielectricResponse::Lossless(&incident),
        1e9,
        &[DielectricLayer {
            response: DielectricResponse::Debye(&lossy),
            thickness_m: thickness,
        }],
    )
    .unwrap();
    let grid = AngularSpectrumGrid {
        width: 8,
        height: 4,
        spacing_m: [2.0, 3.0],
    };
    let mut input = vec![[C64::ZERO; 2]; grid.width * grid.height];
    input[2] = [C64::ONE, C64::new(0.2, -0.3)];
    input[21] = [C64::new(0.1, 0.8), C64::new(-0.4, 0.2)];
    let beam = PlanarAngularSpectrum::prepare(&stack, grid, &input, limits(), || false).unwrap();
    let layer = beam.prepare_layer(0, layer_limits(), || false).unwrap();
    relative(
        layer.absorbed_power_w(),
        beam.power().film_absorbed_w,
        3e-12,
    );
    let left = layer.field(0.0, || false).unwrap();
    let right = layer.field(thickness, || false).unwrap();
    relative(
        left.field().normal_power_w() - right.field().normal_power_w(),
        layer.absorbed_power_w(),
        4e-12,
    );
    let intervals = 128;
    let dz = thickness / intervals as f64;
    let mut integrated = left.absorbed_power_per_depth_w_m() + right.absorbed_power_per_depth_w_m();
    for i in 1..intervals {
        let field = layer.field(i as f64 * dz, || false).unwrap();
        let coefficient = if i % 2 == 0 { 2.0 } else { 4.0 };
        integrated += coefficient * field.absorbed_power_per_depth_w_m();
    }
    relative(integrated * dz / 3.0, layer.absorbed_power_w(), 2e-7);
    let depth = 0.31;
    let plane = layer.field(depth, || false).unwrap();
    let epsilon = stack.layer_permittivity_f_m(0).unwrap();
    for (e, &heat) in plane
        .field()
        .electric_v_per_m()
        .iter()
        .zip(plane.absorbed_power_density_w_m3())
    {
        relative(
            heat,
            -0.5e9 * epsilon.im * e.iter().map(|v| v.norm_sq()).sum::<f64>(),
            2e-15,
        );
    }
    let delta = 1e-5;
    let before = layer
        .field(depth - delta, || false)
        .unwrap()
        .field()
        .normal_power_w();
    let after = layer
        .field(depth + delta, || false)
        .unwrap()
        .field()
        .normal_power_w();
    relative(
        (before - after) / (2.0 * delta),
        plane.absorbed_power_per_depth_w_m(),
        2e-8,
    );
}

#[test]
fn coherent_fourier_modes_make_spatial_heating_fringes_instead_of_an_intensity_sum() {
    let incident = material(E, M);
    let lossy = finite_layers::debye(2.0 * E, E, 1.0, 1.0);
    let stack = PlanarDielectricLayers::try_new(
        &incident,
        DielectricResponse::Lossless(&incident),
        1e9,
        &[DielectricLayer {
            response: DielectricResponse::Debye(&lossy),
            thickness_m: 0.4,
        }],
    )
    .unwrap();
    let grid = AngularSpectrumGrid {
        width: 8,
        height: 1,
        spacing_m: [2.0, 3.0],
    };
    let mut input = Vec::new();
    for x in 0..grid.width {
        let angle = -core::f64::consts::TAU * x as f64 / grid.width as f64;
        input.push([
            C64::ZERO,
            C64::ONE + C64::new(det::cos(angle), det::sin(angle)),
        ]);
    }
    let beam = PlanarAngularSpectrum::prepare(&stack, grid, &input, limits(), || false).unwrap();
    let layer = beam.prepare_layer(0, layer_limits(), || false).unwrap();
    let depth = 0.19;
    let plane = layer.field(depth, || false).unwrap();
    let kx = core::f64::consts::TAU / (grid.width as f64 * grid.spacing_m[0]);
    let k = 1e9 * det::sqrt(E) * det::sqrt(M);
    let cosine = det::sqrt(1.0 - (kx / k) * (kx / k));
    let zero = stack
        .prepare_layer_field(0, 1.0, Polarization::Te)
        .unwrap()
        .field(depth)
        .unwrap();
    let one = stack
        .prepare_layer_field(0, cosine, Polarization::Te)
        .unwrap()
        .field(depth)
        .unwrap();
    let loss = -0.5e9 * stack.layer_permittivity_f_m(0).unwrap().im;
    let incoherent_heat = loss * (zero.electric_v_per_m.norm_sq() + one.electric_v_per_m.norm_sq());
    let mut largest_interference: f64 = 0.0;
    for x in 0..grid.width {
        let angle = -kx * x as f64 * grid.spacing_m[0];
        let expected = zero.electric_v_per_m
            + one.electric_v_per_m * C64::new(det::cos(angle), det::sin(angle));
        complex_near(plane.field().electric_v_per_m()[x][1], expected, 2e-12);
        relative(
            plane.absorbed_power_density_w_m3()[x],
            loss * expected.norm_sq(),
            3e-11,
        );
        largest_interference =
            largest_interference.max((loss * expected.norm_sq() - incoherent_heat).abs());
    }
    assert!(largest_interference > 0.9 * incoherent_heat);
}

#[test]
fn selected_fields_retain_source_critical_limits_opaque_limits_and_atomic_refusals() {
    let incident = material(1.0, 1.0);
    // Exact q² = 3/4 - (1 - (1/2)²) = 0 for the selected film.
    let critical = material(0.75, 1.0);
    let stack = PlanarDielectricLayers::try_new(
        &incident,
        DielectricResponse::Lossless(&incident),
        1.0,
        &[DielectricLayer {
            response: DielectricResponse::Lossless(&critical),
            thickness_m: 2.0,
        }],
    )
    .unwrap();
    for pol in [Polarization::Te, Polarization::Tm] {
        let prepared = stack.prepare_layer_field(0, 0.5, pol).unwrap();
        assert!(std::ptr::eq(prepared.source(), &stack));
        let left = prepared.field(0.0).unwrap();
        let middle = prepared.field(1.0).unwrap();
        let right = prepared.field(2.0).unwrap();
        complex_near(
            middle.electric_v_per_m.scale(2.0),
            left.electric_v_per_m + right.electric_v_per_m,
            3e-14,
        );
        complex_near(
            middle.conjugate_magnetic_a_per_m.scale(2.0),
            left.conjugate_magnetic_a_per_m + right.conjugate_magnetic_a_per_m,
            3e-14,
        );
    }
    let grid = AngularSpectrumGrid {
        width: 4,
        height: 4,
        spacing_m: [8.0, 8.0],
    };
    let beam =
        PlanarAngularSpectrum::prepare(&stack, grid, &[[C64::ONE; 2]; 16], limits(), || false)
            .unwrap();
    let prepared = beam.prepare_layer(0, layer_limits(), || false).unwrap();
    let reference = prepared.field(0.7, || false).unwrap();
    let mut budget = layer_limits();
    budget.maximum_layer_mode_pairs = prepared.preparation_layer_mode_pairs() - 1;
    assert!(beam.prepare_layer(0, budget, || false).is_err());
    let mut storage = layer_limits();
    storage.maximum_workspace_bytes = prepared.workspace_bytes() - 1;
    assert!(beam.prepare_layer(0, storage, || false).is_err());
    let mut fft = layer_limits();
    fft.maximum_fft_work = beam.field_fft_work() - 1;
    assert!(beam.prepare_layer(0, fft, || false).is_err());
    assert!(beam.prepare_layer(1, layer_limits(), || false).is_err());
    assert!(stack.prepare_layer_field(0, 0.0, Polarization::Te).is_err());
    let mut calls = 0;
    assert!(matches!(
        beam.prepare_layer(0, layer_limits(), || {
            calls += 1;
            calls == 8
        }),
        Err(AngularSpectrumError::Cancelled)
    ));
    let mut calls = 0;
    assert!(matches!(
        prepared.field(0.7, || {
            calls += 1;
            calls == 8
        }),
        Err(AngularSpectrumError::Cancelled)
    ));
    for depth in [-1.0, 2.1, f64::NAN] {
        assert!(prepared.field(depth, || false).is_err());
    }
    assert_eq!(
        reference.field().electric_v_per_m(),
        prepared
            .field(0.7, || false)
            .unwrap()
            .field()
            .electric_v_per_m()
    );
    let lossy = finite_layers::debye(2.0 * E, E, 1.0, 1.0);
    let physical = material(E, M);
    let opaque = PlanarDielectricLayers::try_new(
        &physical,
        DielectricResponse::Lossless(&physical),
        1e9,
        &[DielectricLayer {
            response: DielectricResponse::Debye(&lossy),
            thickness_m: 1e6,
        }],
    )
    .unwrap();
    let field = opaque
        .prepare_layer_field(0, 1.0, Polarization::Te)
        .unwrap();
    assert!(field.absorbed_power_w_m2() > 0.0);
    assert_eq!(field.field(5e5).unwrap().electric_v_per_m, C64::ZERO);
}
