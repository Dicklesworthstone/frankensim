use super::*;

mod beam;
mod heterogeneous;

fn slab(depths: &[f64]) -> AbsorbingSlab {
    AbsorbingSlab::new(
        depths
            .iter()
            .enumerate()
            .map(|(i, &depth)| {
                AbsorbingLayer::new(1.0, depth, [100.0, 3000.0], format!("test layer {i}")).unwrap()
            })
            .collect(),
        SlabQuadrature::default(),
        || false,
    )
    .unwrap()
}

fn refractive(depths: &[f64], index: f64) -> AbsorbingSlab {
    AbsorbingSlab::new_refractive(
        depths
            .iter()
            .map(|&depth| {
                AbsorbingLayer::new(1.0, depth, [100.0, 3000.0], "authored absorption").unwrap()
            })
            .collect(),
        index,
        "authored constant real dielectric index",
        SlabQuadrature::default(),
        || false,
    )
    .unwrap()
}

#[test]
fn g0_refractive_transparent_and_opaque_limits_close_energy() {
    let transparent = refractive(&[0.0, 0.0], 1.5);
    // Independent Simpson integration of the lossless incoherent two-interface
    // reflectance 2R/(1+R), evaluated separately for s and p polarization.
    let panels = 16_384;
    let mut expected = 0.0;
    for k in 1..=panels {
        let mu = k as f64 / panels as f64;
        let ct = (1.0 - (1.0 - mu * mu) / 2.25).sqrt();
        let rs = ((mu - 1.5 * ct) / (mu + 1.5 * ct)).powi(2);
        let rp = ((1.5 * mu - ct) / (1.5 * mu + ct)).powi(2);
        let weight = if k == panels {
            1.0
        } else if k % 2 == 0 {
            2.0
        } else {
            4.0
        };
        expected += 2.0 * mu * weight / (3.0 * panels as f64) * (rs / (1.0 + rs) + rp / (1.0 + rp));
    }
    close(transparent.reflectance(), expected, 2e-11);
    close(transparent.transmittance(), 1.0 - expected, 2e-11);
    assert!(transparent.reflectance() > 0.1);
    assert_eq!(transparent.face_coupling(), &[[0.0; 2]; 2]);
    let mut heat = [1.0; 2];
    transparent
        .layer_net_fluxes([900.0, 300.0], &[500.0, 400.0], &mut heat)
        .unwrap();
    assert_eq!(heat, [0.0; 2]);

    let opaque = refractive(&[1000.0, 1000.0], 1.5);
    assert_eq!(opaque.transmittance(), 0.0);
    // Thermal exchange between adjacent opaque layers includes ALL internal
    // modes. Omitting trapped rays gives 1 instead of the required n²=2.25.
    close(opaque.layer_coupling()[1], 2.25, 2e-13);
    close(
        opaque.face_coupling()[0][0],
        1.0 - opaque.reflectance(),
        2e-14,
    );
    assert_eq!(opaque.face_coupling()[0][1], 0.0);
    assert!(
        opaque
            .refractive_index_source()
            .unwrap()
            .contains("dielectric")
    );
}

// Independent two-direction sweeps iterate actual boundary-return intensities
// to a fixed point, then integrate with Simpson nodes. No production escape or
// pair-exchange coefficient formulas are used.
fn refractive_ray_oracle(
    depths: &[f64],
    index: f64,
    b: &[f64],
    h: [f64; 2],
) -> ([f64; 2], Vec<f64>) {
    let n = depths.len();
    let panels = 8192;
    let critical = (1.0 - 1.0 / (index * index)).sqrt();
    let mut outgoing = [0.0; 2];
    let mut heat = vec![0.0; n];
    for trapped in [false, true] {
        let limit = if trapped { critical } else { 1.0 };
        for k in 1..=panels {
            let u = limit * k as f64 / panels as f64;
            let mu = if trapped {
                u
            } else {
                (1.0 - (1.0 - u * u) / (index * index)).sqrt()
            };
            let weight = if k == panels {
                1.0
            } else if k % 2 == 0 {
                2.0
            } else {
                4.0
            };
            let q = u * limit * weight / (3.0 * panels as f64)
                * if trapped { index * index } else { 1.0 };
            let rs = if trapped {
                1.0
            } else {
                ((u - index * mu) / (u + index * mu)).powi(2)
            };
            let rp = if trapped {
                1.0
            } else {
                ((index * u - mu) / (index * u + mu)).powi(2)
            };
            for r in [rs, rp] {
                let mut incoming = [0.0; 2];
                let mut escaped = [0.0; 2];
                for _ in 0..100 {
                    let mut right = incoming[0];
                    let mut left = incoming[1];
                    for i in 0..n {
                        let t = (-depths[i] / mu).exp();
                        right = t * right + (1.0 - t) * b[i];
                    }
                    for i in (0..n).rev() {
                        let t = (-depths[i] / mu).exp();
                        left = t * left + (1.0 - t) * b[i];
                    }
                    escaped = [left, right];
                    let next = if trapped {
                        [left, right]
                    } else {
                        [(1.0 - r) * h[0] + r * left, (1.0 - r) * h[1] + r * right]
                    };
                    let converged = (next[0] - incoming[0])
                        .abs()
                        .max((next[1] - incoming[1]).abs())
                        < 1e-11;
                    incoming = next;
                    if converged {
                        break;
                    }
                }
                let mut forward = incoming[0];
                let mut backward = incoming[1];
                for i in 0..n {
                    let t = (-depths[i] / mu).exp();
                    heat[i] += q * (1.0 - t) * (forward - b[i]);
                    forward = t * forward + (1.0 - t) * b[i];
                }
                for i in (0..n).rev() {
                    let t = (-depths[i] / mu).exp();
                    heat[i] += q * (1.0 - t) * (backward - b[i]);
                    backward = t * backward + (1.0 - t) * b[i];
                }
                if !trapped {
                    for side in 0..2 {
                        outgoing[side] += q * (r * h[side] + (1.0 - r) * escaped[side]);
                    }
                }
            }
        }
    }
    (outgoing, heat)
}

#[test]
fn g1_refractive_nonisothermal_heat_matches_independent_reflection_sweeps() {
    let depths = [0.1, 0.7, 1.3, 0.0];
    let b = [9000.0, 3000.0, 700.0, 12000.0];
    let h = [1500.0, 400.0];
    let material = refractive(&depths, 1.5);
    let (outgoing, expected) = refractive_ray_oracle(&depths, 1.5, &b, h);
    let emission = material.emission_to_faces(&b).unwrap();
    let mut actual = [0.0; 4];
    material.layer_net_fluxes(h, &b, &mut actual).unwrap();
    for i in 0..4 {
        close(actual[i], expected[i], 5e-6);
    }
    let computed: [f64; 2] = std::array::from_fn(|i| {
        emission[i] + material.reflectance() * h[i] + material.transmittance() * h[1 - i]
    });
    for side in 0..2 {
        close(computed[side], outgoing[side], 5e-6);
    }
    close(
        actual.iter().sum(),
        h.iter().sum::<f64>() - computed.iter().sum::<f64>(),
        2e-10,
    );
    assert_eq!(actual[3], 0.0);
}

#[test]
fn g3_refractive_subdivision_reversal_and_equilibrium_preserve_physics() {
    let whole = refractive(&[1.0], 1.5);
    let split = refractive(&[0.2, 0.3, 0.5], 1.5);
    let reverse = refractive(&[0.5, 0.3, 0.2], 1.5);
    close(split.transmittance(), whole.transmittance(), 2e-11);
    close(split.reflectance(), whole.reflectance(), 2e-11);
    let mut split_heat = [0.0; 3];
    let mut reverse_heat = [0.0; 3];
    let mut whole_heat = [0.0];
    whole
        .layer_net_fluxes([1700.0, 800.0], &[450.0], &mut whole_heat)
        .unwrap();
    split
        .layer_net_fluxes([1700.0, 800.0], &[450.0; 3], &mut split_heat)
        .unwrap();
    close(split_heat.iter().sum(), whole_heat[0], 1e-7);
    split
        .layer_net_fluxes([1700.0, 800.0], &[900.0, 450.0, 100.0], &mut split_heat)
        .unwrap();
    reverse
        .layer_net_fluxes([800.0, 1700.0], &[100.0, 450.0, 900.0], &mut reverse_heat)
        .unwrap();
    for i in 0..3 {
        close(split_heat[i], reverse_heat[2 - i], 1e-8);
    }
    split
        .layer_net_fluxes([1e100; 2], &[1e100; 3], &mut split_heat)
        .unwrap();
    assert_eq!(split_heat, [0.0; 3]);
    let matched = refractive(&[0.2, 0.3], 1.0);
    assert_eq!(matched.transmittance(), slab(&[0.2, 0.3]).transmittance());
    assert_eq!(matched.reflectance(), 0.0);
}

#[test]
fn g0_refractive_inputs_and_work_limits_refuse_unsupported_models() {
    let layer = AbsorbingLayer::new(1.0, 0.3, [100.0, 3000.0], "test").unwrap();
    for (index, source) in [
        (0.9, "test"),
        (f64::NAN, "test"),
        (f64::INFINITY, "test"),
        (1e200, "test"),
        (1.5, " "),
    ] {
        assert!(
            AbsorbingSlab::new_refractive(
                vec![layer.clone()],
                index,
                source,
                SlabQuadrature::default(),
                || false
            )
            .is_err()
        );
    }
    assert!(
        AbsorbingSlab::new_refractive(
            vec![layer.clone()],
            1.5,
            "test",
            SlabQuadrature {
                maximum_evaluations: 0,
                ..SlabQuadrature::default()
            },
            || false
        )
        .is_err()
    );
    assert!(matches!(
        AbsorbingSlab::new_refractive(vec![layer], 1.5, "test", SlabQuadrature::default(), || true),
        Err(ConductionError::Cancelled { .. })
    ));
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.17e}, expected={expected:.17e}, tolerance={tolerance:.3e}"
    );
}

// Independent direct ray sweeps with composite Simpson angular integration.
// This oracle never forms escape/exchange coefficients or uses production
// angular nodes. Each ray crosses every physical layer before being averaged.
fn ray_oracle(depths: &[f64], b: &[f64], h: [f64; 2]) -> ([f64; 2], Vec<f64>) {
    let n = depths.len();
    let panels = 16_384;
    let mut outgoing = [0.0; 2];
    let mut heat = vec![0.0; n];
    for node in 1..=panels {
        let mu = node as f64 / panels as f64;
        let weight = if node == panels {
            1.0
        } else if node % 2 == 0 {
            2.0
        } else {
            4.0
        };
        let q = 2.0 * mu * weight / (3.0 * panels as f64);
        let mut forward = h[0];
        for i in 0..n {
            let depth = depths[i] / mu;
            let absorption = -(-depth).exp_m1();
            heat[i] += q * absorption * (forward - b[i]);
            forward = (-depth).exp() * forward + absorption * b[i];
        }
        outgoing[1] += q * forward;
        let mut backward = h[1];
        for i in (0..n).rev() {
            let depth = depths[i] / mu;
            let absorption = -(-depth).exp_m1();
            heat[i] += q * absorption * (backward - b[i]);
            backward = (-depth).exp() * backward + absorption * b[i];
        }
        outgoing[0] += q * backward;
    }
    (outgoing, heat)
}

#[test]
fn homogeneous_diffuse_transmission_has_the_exponential_integral_limit() {
    let material = slab(&[1.0]);
    // 2 E_3(1) = E_1(1); NIST DLMF 8.19.3 and its recurrence.
    close(material.transmittance(), 0.21938393439552027, 2e-11);
    assert!((material.transmittance() - (-1.0_f64).exp()).abs() > 0.1);
    for side in 0..2 {
        close(
            material.face_coupling()[0][side],
            1.0 - material.transmittance(),
            2e-15,
        );
    }
    assert_eq!(material.layer_coupling(), &[0.0]);
    assert!(material.quadrature().angular_evaluations > 0);
    assert!(material.quadrature().maximum_partition_residual < 2e-15);
}

#[test]
fn subdividing_isothermal_material_preserves_the_full_angular_path() {
    let whole = slab(&[1.0]);
    let split = slab(&[0.2, 0.3, 0.5]);
    close(split.transmittance(), whole.transmittance(), 2e-11);
    let mixed =
        slab(&[0.2]).transmittance() * slab(&[0.3]).transmittance() * slab(&[0.5]).transmittance();
    assert!((mixed - split.transmittance()).abs() > 0.01);
    let mut whole_heat = [0.0];
    let mut split_heat = [0.0; 3];
    whole
        .layer_net_fluxes([1700.0, 800.0], &[450.0], &mut whole_heat)
        .unwrap();
    split
        .layer_net_fluxes([1700.0, 800.0], &[450.0; 3], &mut split_heat)
        .unwrap();
    close(split_heat.iter().sum(), whole_heat[0], 1e-7);
    let whole_emission = whole.emission_to_faces(&[450.0]).unwrap();
    let split_emission = split.emission_to_faces(&[450.0; 3]).unwrap();
    for side in 0..2 {
        close(split_emission[side], whole_emission[side], 1e-8);
    }
}

#[test]
fn nonisothermal_internal_heat_matches_independent_two_direction_ray_sweeps() {
    let depths = [0.08, 0.7, 2.0, 0.0];
    let b = [9000.0, 3000.0, 700.0, 12000.0];
    let h = [1500.0, 400.0];
    let material = slab(&depths);
    let (outgoing, expected) = ray_oracle(&depths, &b, h);
    let emission = material.emission_to_faces(&b).unwrap();
    let actual_outgoing = [
        material.transmittance() * h[1] + emission[0],
        material.transmittance() * h[0] + emission[1],
    ];
    let mut actual = [0.0; 4];
    material.layer_net_fluxes(h, &b, &mut actual).unwrap();
    for i in 0..4 {
        close(actual[i], expected[i], 2e-7);
    }
    for side in 0..2 {
        close(actual_outgoing[side], outgoing[side], 2e-7);
    }
    close(
        actual.iter().sum(),
        h[0] + h[1] - actual_outgoing[0] - actual_outgoing[1],
        2e-10,
    );
    assert_eq!(actual[3], 0.0);
    assert!(actual[0] < 0.0);
    assert!(actual[2] > 0.0);
}

#[test]
fn reversal_swaps_ports_and_original_layer_heat_without_changing_transport() {
    let forward = slab(&[0.1, 1.2, 0.4]);
    let reverse = slab(&[0.4, 1.2, 0.1]);
    close(forward.transmittance(), reverse.transmittance(), 1e-14);
    let mut a = [0.0; 3];
    let mut b = [0.0; 3];
    forward
        .layer_net_fluxes([2000.0, 600.0], &[1000.0, 400.0, 700.0], &mut a)
        .unwrap();
    reverse
        .layer_net_fluxes([600.0, 2000.0], &[700.0, 400.0, 1000.0], &mut b)
        .unwrap();
    for i in 0..3 {
        close(a[i], b[2 - i], 2e-10);
        for side in 0..2 {
            close(
                forward.face_coupling()[i][side],
                reverse.face_coupling()[2 - i][1 - side],
                2e-14,
            );
        }
        for j in 0..3 {
            assert_eq!(
                forward.layer_coupling()[i * 3 + j],
                forward.layer_coupling()[j * 3 + i]
            );
            assert!(forward.layer_coupling()[i * 3 + j] >= 0.0);
        }
    }
}

#[test]
fn vacuum_thin_and_opaque_limits_keep_absorption_and_internal_exchange() {
    let vacuum = AbsorbingSlab::new(
        vec![AbsorbingLayer::new(0.01, 0.0, [100.0, 3000.0], "vacuum").unwrap()],
        SlabQuadrature {
            maximum_evaluations: 0,
            ..SlabQuadrature::default()
        },
        || false,
    )
    .unwrap();
    assert_eq!(vacuum.transmittance(), 1.0);
    assert_eq!(vacuum.quadrature().angular_evaluations, 0);
    let mut heat = [123.0];
    vacuum
        .layer_net_fluxes([1e100, -1e100], &[1e100], &mut heat)
        .unwrap();
    assert_eq!(heat, [0.0]);

    let thin = slab(&[1e-12]);
    for side in 0..2 {
        close(thin.face_coupling()[0][side] / 2e-12, 1.0, 5e-10);
    }
    // Finite properties may multiply to infinite optical depth. The limiting
    // response is opaque; authored infinite properties themselves are refused.
    let opaque_layer = AbsorbingLayer::new(1e200, 1e200, [100.0, 3000.0], "opaque limit").unwrap();
    let opaque =
        AbsorbingSlab::new(vec![opaque_layer; 3], SlabQuadrature::default(), || false).unwrap();
    assert_eq!(opaque.transmittance(), 0.0);
    close(opaque.layer_coupling()[1], 1.0, 5e-15);
    close(opaque.layer_coupling()[5], 1.0, 5e-15);
    assert_eq!(opaque.layer_coupling()[2], 0.0);
    let mut internal = [0.0; 3];
    opaque
        .layer_net_fluxes([100.0, 100.0], &[900.0, 400.0, 100.0], &mut internal)
        .unwrap();
    for (actual, expected) in internal.into_iter().zip([-1300.0, 200.0, 300.0]) {
        close(actual, expected, 1e-11);
    }
}

#[test]
fn common_background_form_retains_equilibrium_and_small_signed_heat() {
    let material = slab(&[0.4, 0.7]);
    let mut heat = [1.0; 2];
    material
        .layer_net_fluxes([1e100; 2], &[1e100; 2], &mut heat)
        .unwrap();
    assert_eq!(heat, [0.0; 2]);
    let mut shifted = [0.0; 2];
    let b = [40.0, -30.0];
    material
        .layer_net_fluxes([3.0, -5.0], &b, &mut shifted)
        .unwrap();
    material
        .layer_net_fluxes([1003.0, 995.0], &[1040.0, 970.0], &mut heat)
        .unwrap();
    assert_eq!(heat, shifted);
    material
        .layer_net_fluxes([0.0; 2], &[1e-10, 0.0], &mut heat)
        .unwrap();
    assert!(heat[0] < 0.0 && heat[1] > 0.0);
}

#[test]
fn material_temperature_axes_and_nonfinite_results_are_admitted_atomically() {
    for (d, k, support, source) in [
        (0.0, 1.0, [100.0, 1000.0], "source"),
        (1.0, -1.0, [100.0, 1000.0], "source"),
        (1.0, f64::INFINITY, [100.0, 1000.0], "source"),
        (1.0, 1.0, [1000.0, 100.0], "source"),
        (1.0, 1.0, [100.0, 1000.0], " "),
    ] {
        assert!(AbsorbingLayer::new(d, k, support, source).is_err());
    }
    let material = slab(&[0.2, 0.5]);
    material.validate_temperatures(&[100.0, 3000.0]).unwrap();
    assert!(material.validate_temperatures(&[99.0, 3000.0]).is_err());
    assert!(material.validate_temperatures(&[300.0]).is_err());
    assert!(material.emission_to_faces(&[f64::NAN, 400.0]).is_err());
    let mut heat = [123.0, 456.0];
    assert!(
        material
            .layer_net_fluxes([400.0, f64::INFINITY], &[300.0; 2], &mut heat)
            .is_err()
    );
    assert_eq!(heat, [123.0, 456.0]);
    assert!(
        material
            .layer_net_fluxes([400.0; 2], &[300.0], &mut heat)
            .is_err()
    );
    assert_eq!(heat, [123.0, 456.0]);
}

#[test]
fn angular_preparation_refuses_exhausted_budgets_and_honors_cancellation() {
    let layer = AbsorbingLayer::new(1.0, 1.0, [100.0, 3000.0], "test").unwrap();
    assert!(
        AbsorbingSlab::new(
            vec![layer.clone()],
            SlabQuadrature {
                maximum_evaluations: 23,
                ..SlabQuadrature::default()
            },
            || false,
        )
        .is_err()
    );
    assert!(
        AbsorbingSlab::new(
            vec![layer.clone()],
            SlabQuadrature {
                absolute_coefficient_error: 0.0,
                relative_coefficient_error: 0.0,
                ..SlabQuadrature::default()
            },
            || false,
        )
        .is_err()
    );
    let mut calls = 0;
    let cancelled = AbsorbingSlab::new(vec![layer], SlabQuadrature::default(), || {
        calls += 1;
        calls > 8
    });
    assert!(matches!(
        cancelled,
        Err(ConductionError::Cancelled {
            stage: "slab-angular-prepare",
            ..
        })
    ));
    assert!(AbsorbingSlab::new(Vec::new(), SlabQuadrature::default(), || false).is_err());
}
