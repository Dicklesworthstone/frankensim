//! G1/G3 spatial coherence: independent Maxwell fringes and Airy heat, with
//! original prepared diffraction as the consumer rather than a scalar sum.
use super::*;
use fs_em::plane_wave::{
    angular_spectrum::{
        AngularSpectrumCoherence, AngularSpectrumCoherenceField, AngularSpectrumCoherenceLimits,
        AngularSpectrumError, AngularSpectrumGrid, AngularSpectrumLayerLimits,
        AngularSpectrumLimits, AngularSpectrumSide as Side, PlanarAngularSpectrum,
    },
    layers::{DielectricLayer, PlanarDielectricLayers},
};
use fs_math::det;

fn limits() -> AngularSpectrumLimits {
    AngularSpectrumLimits {
        maximum_grid_points: 128,
        maximum_workspace_bytes: 1 << 20,
        maximum_fft_work: 1 << 20,
    }
}
fn coherence_limits() -> AngularSpectrumCoherenceLimits {
    AngularSpectrumCoherenceLimits {
        maximum_modes: 8,
        maximum_point_mode_pairs: 1024,
        maximum_workspace_bytes: 2 << 20,
        maximum_fft_work: 1 << 20,
    }
}
fn layer_limits() -> AngularSpectrumLayerLimits {
    AngularSpectrumLayerLimits {
        maximum_layer_mode_pairs: 4096,
        maximum_workspace_bytes: 4 << 20,
        maximum_fft_work: 1 << 20,
    }
}
fn grid() -> AngularSpectrumGrid {
    AngularSpectrumGrid {
        width: 8,
        height: 1,
        spacing_m: [8.0, 3.0],
    }
}
fn phase(x: f64) -> C64 {
    C64::new(det::cos(x), det::sin(x))
}
fn prepare<'s, 'm>(
    source: &'s PlanarDielectricLayers<'m>,
    samples: &[[C64; 2]],
) -> PlanarAngularSpectrum<'s, 'm> {
    PlanarAngularSpectrum::prepare(source, grid(), samples, limits(), || false).unwrap()
}
fn close(a: f64, b: f64) {
    near(a, b, 2e-10 * a.abs().max(b.abs()).max(1e-18));
}

#[test]
fn g1_spatial_coherence_sets_physical_fringe_visibility_without_changing_total_power() {
    let medium = material(1.0, 1.0);
    let stack =
        PlanarDielectricLayers::try_new(&medium, DielectricResponse::Lossless(&medium), 1.0, &[])
            .unwrap();
    let g = grid();
    let kx = core::f64::consts::TAU / (g.width as f64 * g.spacing_m[0]);
    let q = det::sqrt(1.0 - kx * kx);
    for rho in [0.0, 0.6, 1.0] {
        let mut a = Vec::new();
        let mut b = Vec::new();
        for x in 0..g.width {
            let p = phase(-kx * x as f64 * g.spacing_m[0]);
            a.push([C64::ZERO, C64::ONE + p.scale(rho)]);
            b.push([C64::ZERO, p.scale(det::sqrt(1.0 - rho * rho))]);
        }
        let factors = [prepare(&stack, &a), prepare(&stack, &b)];
        let bank =
            AngularSpectrumCoherence::try_new(&factors[0], &factors, coherence_limits(), || false)
                .unwrap();
        let field = bank.field(Side::Transmitted, 0.0, || false).unwrap();
        let mean_flux = 0.5 * (1.0 + q);
        for (x, point) in field.points().iter().enumerate() {
            let fringe = rho * det::cos(kx * x as f64 * g.spacing_m[0]);
            near(point.electric_coherency_v2_m2[1][1].re, 1.0 + fringe, 2e-13);
            near(
                point.complex_poynting_w_m2[2].re,
                mean_flux * (1.0 + fringe),
                3e-13,
            );
            assert_eq!(point.electric_coherency_v2_m2[0][0], C64::ZERO);
        }
        close(
            field.normal_power_w(),
            mean_flux * g.width as f64 * g.spacing_m[0] * g.spacing_m[1],
        );
        close(field.normal_power_w(), bank.power().transmitted_w);
        assert_eq!(bank.power().film_absorbed_w, 0.0);
    }
}

fn same_field(a: &AngularSpectrumCoherenceField, b: &AngularSpectrumCoherenceField) {
    close(a.normal_power_w(), b.normal_power_w());
    for (a, b) in a.points().iter().zip(b.points()) {
        for i in 0..3 {
            close(
                a.magnetic_rms_squared_a2_m2[i],
                b.magnetic_rms_squared_a2_m2[i],
            );
            complex_near(
                a.complex_poynting_w_m2[i],
                b.complex_poynting_w_m2[i],
                2e-12,
            );
            for j in 0..3 {
                complex_near(
                    a.electric_coherency_v2_m2[i][j],
                    b.electric_coherency_v2_m2[i][j],
                    2e-12,
                );
            }
        }
    }
}

#[test]
fn g3_complex_factor_rotation_preserves_polarization_flux_and_internal_heat() {
    let medium = material(E, M);
    let film = finite_layers::debye(2.0 * E, E, 1.0, 1.0);
    let stack = PlanarDielectricLayers::try_new(
        &medium,
        DielectricResponse::Lossless(&medium),
        1e9,
        &[DielectricLayer {
            response: DielectricResponse::Debye(&film),
            thickness_m: 0.4,
        }],
    )
    .unwrap();
    let mut a = Vec::new();
    let mut b = Vec::new();
    let mut ua = Vec::new();
    let mut ub = Vec::new();
    let j = C64::new(0.0, 1.0);
    for x in 0..grid().width {
        let p = phase(-core::f64::consts::TAU * x as f64 / grid().width as f64);
        let av = [C64::new(0.8, 0.2) + p.scale(0.3), j];
        let bv = [p * C64::new(0.2, -0.6), p.scale(0.4) + C64::new(0.2, 0.1)];
        // Apply a complex unitary matrix to SOURCE fields, before any FFT.
        ua.push(std::array::from_fn(|c| {
            (av[c] + j * bv[c]).scale(core::f64::consts::FRAC_1_SQRT_2)
        }));
        ub.push(std::array::from_fn(|c| {
            (j * av[c] + bv[c]).scale(core::f64::consts::FRAC_1_SQRT_2)
        }));
        a.push(av);
        b.push(bv);
    }
    let factors = [prepare(&stack, &a), prepare(&stack, &b)];
    let rotated = [prepare(&stack, &ua), prepare(&stack, &ub)];
    let bank =
        AngularSpectrumCoherence::try_new(&factors[0], &factors, coherence_limits(), || false)
            .unwrap();
    let other =
        AngularSpectrumCoherence::try_new(&rotated[0], &rotated, coherence_limits(), || false)
            .unwrap();
    for side in [
        Side::Incident,
        Side::Reflected,
        Side::TotalIncidentSide,
        Side::Transmitted,
    ] {
        same_field(
            &bank.field(side, 0.31, || false).unwrap(),
            &other.field(side, 0.31, || false).unwrap(),
        );
    }
    let layer = bank.prepare_layer(0, layer_limits(), || false).unwrap();
    let changed = other.prepare_layer(0, layer_limits(), || false).unwrap();
    assert_eq!(layer.preparation_layer_mode_pairs(), 2 * grid().width);
    close(layer.absorbed_power_w(), bank.power().film_absorbed_w);
    close(layer.absorbed_power_w(), changed.absorbed_power_w());
    let left = layer.field(0.19, || false).unwrap();
    let right = changed.field(0.19, || false).unwrap();
    same_field(left.field(), right.field());
    for (&q, &r) in left
        .absorbed_power_density_w_m3()
        .iter()
        .zip(right.absorbed_power_density_w_m3())
    {
        close(q, r);
    }
    let loss = -1e9 * stack.layer_permittivity_f_m(0).unwrap().im;
    for (point, &q) in left
        .field()
        .points()
        .iter()
        .zip(left.absorbed_power_density_w_m3())
    {
        close(
            q,
            loss * (0..3)
                .map(|c| point.electric_coherency_v2_m2[c][c].re)
                .sum::<f64>(),
        );
        assert!(q > 0.0);
    }
    let incident = bank.field(Side::Incident, 0.0, || false).unwrap();
    let ex_ey = (a[0][0] * a[0][1].conj() + b[0][0] * b[0][1].conj()).scale(0.5);
    complex_near(
        incident.points()[0].electric_coherency_v2_m2[0][1],
        ex_ey,
        2e-13,
    );
    assert!(ex_ey.im.abs() > 0.1);
}

// Independent one-film Airy amplitudes, retaining each internal round trip.
fn airy_te(epsilon: C64, kx: f64, depth: f64, thickness: f64) -> C64 {
    let omega = 1e9;
    let q0 = det::sqrt(omega * omega * E * M - kx * kx);
    let q = (epsilon.scale(omega * omega * M) - C64::from_re(kx * kx)).sqrt();
    assert!(q.re > 0.0 && q.im < 0.0);
    let travel = |z: f64| phase(-q.re * z).scale(det::exp(q.im * z));
    let r01 = (C64::from_re(q0) - q) / (C64::from_re(q0) + q);
    let r12 = -r01;
    let amplitude = C64::from_re(2.0 * q0)
        / (C64::from_re(q0) + q)
        / (C64::ONE + r01 * r12 * travel(2.0 * thickness));
    amplitude * (travel(depth) + r12 * travel(2.0 * thickness - depth))
}

#[test]
fn g1_partial_coherence_heat_matches_independent_airy_and_depth_integrated_loss() {
    let medium = material(E, M);
    let film = finite_layers::debye(2.0 * E, E, 1.0, 1.0);
    let d = 0.4;
    let stack = PlanarDielectricLayers::try_new(
        &medium,
        DielectricResponse::Lossless(&medium),
        1e9,
        &[DielectricLayer {
            response: DielectricResponse::Debye(&film),
            thickness_m: d,
        }],
    )
    .unwrap();
    let rho = 0.65;
    let g = grid();
    let kx = core::f64::consts::TAU / (g.width as f64 * g.spacing_m[0]);
    let a: Vec<_> = (0..g.width)
        .map(|x| {
            [
                C64::ZERO,
                C64::ONE + phase(-kx * x as f64 * g.spacing_m[0]).scale(rho),
            ]
        })
        .collect();
    let b: Vec<_> = (0..g.width)
        .map(|x| {
            [
                C64::ZERO,
                phase(-kx * x as f64 * g.spacing_m[0]).scale(det::sqrt(1.0 - rho * rho)),
            ]
        })
        .collect();
    let factors = [prepare(&stack, &a), prepare(&stack, &b)];
    let bank =
        AngularSpectrumCoherence::try_new(&factors[0], &factors, coherence_limits(), || false)
            .unwrap();
    let layer = bank.prepare_layer(0, layer_limits(), || false).unwrap();
    let epsilon = stack.layer_permittivity_f_m(0).unwrap();
    let depth = 0.19;
    let u = airy_te(epsilon, 0.0, depth, d);
    let v = airy_te(epsilon, kx, depth, d);
    let plane = layer.field(depth, || false).unwrap();
    for x in 0..g.width {
        let p = v * phase(-kx * x as f64 * g.spacing_m[0]);
        let expected =
            -0.5e9 * epsilon.im * ((u + p.scale(rho)).norm_sq() + (1.0 - rho * rho) * p.norm_sq());
        close(plane.absorbed_power_density_w_m3()[x], expected);
    }
    let intervals = 64;
    let mut integral = 0.0;
    for i in 0..=intervals {
        let weight = if i == 0 || i == intervals {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        integral += weight
            * layer
                .field(d * i as f64 / intervals as f64, || false)
                .unwrap()
                .absorbed_power_per_depth_w_m();
    }
    near(
        integral * d / (3.0 * intervals as f64),
        layer.absorbed_power_w(),
        3e-7 * layer.absorbed_power_w(),
    );
    close(layer.absorbed_power_w(), bank.power().film_absorbed_w);
}

#[test]
fn g0_dark_rank_zero_source_binding_budgets_and_cancellation() {
    let medium = material(1.0, 1.0);
    let make_stack = || {
        PlanarDielectricLayers::try_new(
            &medium,
            DielectricResponse::Lossless(&medium),
            1.0,
            &[DielectricLayer {
                response: DielectricResponse::Lossless(&medium),
                thickness_m: 0.4,
            }],
        )
        .unwrap()
    };
    let source = make_stack();
    let distinct_source = make_stack();
    let input = vec![[C64::ONE, C64::ZERO]; grid().width];
    let anchor = prepare(&source, &input);
    let foreign = [prepare(&distinct_source, &input)];
    assert!(
        AngularSpectrumCoherence::try_new(&anchor, &foreign, coherence_limits(), || false).is_err()
    );
    let dark = AngularSpectrumCoherence::try_new(
        &anchor,
        &[],
        AngularSpectrumCoherenceLimits {
            maximum_modes: 0,
            maximum_point_mode_pairs: 0,
            maximum_fft_work: 0,
            ..coherence_limits()
        },
        || false,
    )
    .unwrap();
    assert_eq!(dark.power().incident_w, 0.0);
    let zero = dark.field(Side::TotalIncidentSide, 0.2, || false).unwrap();
    assert_eq!(zero.normal_power_w(), 0.0);
    assert!(
        zero.points()
            .iter()
            .all(|p| p.electric_coherency_v2_m2 == [[C64::ZERO; 3]; 3]
                && p.complex_poynting_w_m2 == [C64::ZERO; 3])
    );
    let layer = dark.prepare_layer(0, layer_limits(), || false).unwrap();
    assert_eq!(layer.absorbed_power_w(), 0.0);
    assert!(
        layer
            .field(0.2, || false)
            .unwrap()
            .absorbed_power_density_w_m3()
            .iter()
            .all(|&q| q == 0.0)
    );
    assert!(dark.field(Side::Transmitted, -1.0, || false).is_err());
    let factors = [prepare(&source, &input), prepare(&source, &input)];
    for limits in [
        AngularSpectrumCoherenceLimits {
            maximum_modes: 1,
            ..coherence_limits()
        },
        AngularSpectrumCoherenceLimits {
            maximum_point_mode_pairs: grid().width,
            ..coherence_limits()
        },
        AngularSpectrumCoherenceLimits {
            maximum_fft_work: factors[0].field_fft_work(),
            ..coherence_limits()
        },
        AngularSpectrumCoherenceLimits {
            maximum_workspace_bytes: 1,
            ..coherence_limits()
        },
    ] {
        assert!(AngularSpectrumCoherence::try_new(&anchor, &factors, limits, || false).is_err());
    }
    let bank =
        AngularSpectrumCoherence::try_new(&anchor, &factors, coherence_limits(), || false).unwrap();
    assert!(matches!(
        bank.field(Side::Transmitted, 0.0, || true),
        Err(AngularSpectrumError::Cancelled)
    ));
    let before = bank.power();
    let mut polls = 0;
    assert!(matches!(
        bank.field(Side::Transmitted, 0.0, || {
            polls += 1;
            polls > 10
        }),
        Err(AngularSpectrumError::Cancelled)
    ));
    assert_eq!(bank.power(), before);
    close(
        bank.field(Side::Transmitted, 0.0, || false)
            .unwrap()
            .normal_power_w(),
        before.transmitted_w,
    );
    assert!(
        bank.prepare_layer(
            0,
            AngularSpectrumLayerLimits {
                maximum_fft_work: bank.field_fft_work() - 1,
                ..layer_limits()
            },
            || false
        )
        .is_err()
    );
    let layer = bank.prepare_layer(0, layer_limits(), || false).unwrap();
    assert_eq!(layer.absorbed_power_w(), 0.0);
    assert!(
        layer
            .field(0.2, || false)
            .unwrap()
            .absorbed_power_density_w_m3()
            .iter()
            .all(|&q| q == 0.0)
    );
}
