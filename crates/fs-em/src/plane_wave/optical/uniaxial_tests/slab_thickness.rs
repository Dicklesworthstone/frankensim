//! Original source laws, physical metre scale and both incidence APIs.
use super::*;
use crate::plane_wave::uniaxial::oriented::{
    OrientedUniaxialInterface, slab::OrientedUniaxialSlab,
};

#[test]
fn original_crystal_ports_retain_sourced_thickness_jets_at_spectral_endpoints() {
    let (ordinary, extraordinary) = axes(
        [1.4, 1.6],
        [1.8, 2.0],
        [0.0; 2],
        Some(MU),
        [400.0, 800.0],
        0,
    );
    let crystal = UniaxialOpticalResponse::try_new(&ordinary, &extraordinary, SPAN).unwrap();
    let outside = incident(1.0, MU);
    let exit = incident(1.2, 0.8 * MU);
    let entrance = OrientedUniaxialInterface::try_new(
        DielectricResponse::Lossless(&outside),
        crystal,
        [1.0, 2.0, 3.0],
    )
    .unwrap();
    let depth = 0.5e-6;
    let step = 1e-10;
    let slab = OrientedUniaxialSlab::try_new(&entrance, DielectricResponse::Lossless(&exit), depth)
        .unwrap();
    let electric = [C64::new(0.7, 0.2), C64::new(-0.3, 0.5)];
    // Thickness differentiation never asks the material for an out-of-band
    // neighboring frequency, so both original inclusive endpoints remain valid.
    for wavelength in [400.0, 550.0, 800.0] {
        let omega = frequency(wavelength);
        for cosine_api in [false, true] {
            let actual = if cosine_api {
                slab.thickness_derivatives(omega, 0.8, electric).unwrap()
            } else {
                slab.thickness_derivatives_at_tangential_index(omega, -0.4, electric)
                    .unwrap()
            };
            let samples: [_; 5] = core::array::from_fn(|j| {
                let variant = OrientedUniaxialSlab::try_new(
                    &entrance,
                    DielectricResponse::Lossless(&exit),
                    depth + (j as f64 - 2.0) * step,
                )
                .unwrap();
                if cosine_api {
                    variant.response(omega, 0.8, electric).unwrap()
                } else {
                    variant
                        .response_at_tangential_index(omega, -0.4, electric)
                        .unwrap()
                }
            });
            assert_eq!(actual.response, samples[2]);
            for component in 0..3 {
                let [mm, m, c, p, pp] = samples.map(|s| s.transmitted_electric_v_m[component]);
                let first = (mm - m.scale(8.0) + p.scale(8.0) - pp).scale(1.0 / (12.0 * step));
                let second = (-pp + p.scale(16.0) - c.scale(30.0) + m.scale(16.0) - mm)
                    .scale(1.0 / (12.0 * step * step));
                assert!((first - actual.transmitted.electric_first_v_m2[component]).abs() < 0.03);
                assert!((second - actual.transmitted.electric_second_v_m3[component]).abs() < 2e6);
            }
            assert!(
                (actual.reflected_power_first_per_m + actual.transmitted_power_first_per_m).abs()
                    < 1e-6
            );
            assert!(
                (actual.reflected_power_second_per_m2 + actual.transmitted_power_second_per_m2)
                    .abs()
                    < 100.0
            );
        }
    }
    assert!(
        slab.thickness_derivatives(frequency(900.0), 0.8, electric)
            .is_err()
    );
    assert!(
        slab.thickness_derivatives(frequency(550.0), 0.0, electric)
            .is_err()
    );
    assert!(
        slab.thickness_derivatives_at_tangential_index(frequency(550.0), 2.0, electric)
            .is_err()
    );
}
