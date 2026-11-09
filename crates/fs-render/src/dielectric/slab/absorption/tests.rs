//! G1/G3: physical depth localization, independent bounces and additive heat.

use super::*;
use crate::dielectric::{BeerLambertAbsorption, CauchyIor, DielectricGlass, GlassProvenance};

fn glass(alpha: f64) -> DielectricGlass {
    DielectricGlass::new(
        CauchyIor::try_constant(1.5).unwrap(),
        BeerLambertAbsorption::try_constant(alpha).unwrap(),
        GlassProvenance::Custom,
    )
}

fn relative(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(f64::MIN_POSITIVE),
        "{actual:e} != {expected:e}"
    );
}

fn reflections(ni: f64, nt: f64, ci: f64, ct: f64) -> [f64; 2] {
    [
        ((ni * ci - nt * ct) / (ni * ci + nt * ct)).powi(2),
        ((nt * ci - ni * ct) / (nt * ci + ni * ct)).powi(2),
    ]
}

#[test]
fn g1_depth_cells_match_independent_forward_and_backward_bounce_integrals() {
    let alpha = 42.0;
    let thickness = 0.018;
    for (incident, exit, cosine) in [
        (1.0_f64, 1.0_f64, 0.3_f64),
        (1.0, 1.33, 0.5),
        (1.3, 1.0, 0.7),
        (1.3, 1.0, 0.4), // Rear total reflection: returning loss is essential.
    ] {
        let slab = glass(alpha)
            .prepare_incoherent_slab_absorption_between(550.0, thickness, cosine, incident, exit)
            .unwrap();
        let transverse_squared = incident.powi(2) * (1.0 - cosine.powi(2));
        let internal_cosine = (1.0 - transverse_squared / 2.25).sqrt();
        let exit_cosine_squared = 1.0 - transverse_squared / exit.powi(2);
        let front = reflections(incident, 1.5, cosine, internal_cosine);
        let rear = if exit_cosine_squared <= 0.0 {
            [1.0; 2]
        } else {
            reflections(1.5, exit, internal_cosine, exit_cosine_squared.sqrt())
        };
        let attenuation = (-alpha * thickness / internal_cosine).exp();
        let edges = [0.0, 0.002, 0.009, 0.017, thickness];
        let mut totals = [0.0; 2];
        for interval in edges.windows(2) {
            let [z0, z1] = [interval[0], interval[1]];
            let actual = slab.absorption_between(z0, z1).unwrap();
            let mut expected = [0.0; 2];
            for pol in 0..2 {
                let mut incoming = 1.0 - front[pol];
                for pass in 0..200 {
                    let (begin, end, reflection) = if pass % 2 == 0 {
                        (z0, z1, rear[pol])
                    } else {
                        (thickness - z1, thickness - z0, front[pol])
                    };
                    expected[pol] += incoming
                        * ((-alpha * begin / internal_cosine).exp()
                            - (-alpha * end / internal_cosine).exp());
                    incoming *= attenuation * reflection;
                }
                relative(actual[pol], expected[pol], 2e-14);
                assert!(actual[pol] >= 0.0);
                totals[pol] += actual[pol];
            }
        }
        relative(totals[0], slab.optical().perpendicular.absorptance, 3e-15);
        relative(totals[1], slab.optical().parallel.absorptance, 3e-15);
        assert_eq!(slab.thickness_m(), thickness);
    }
}

#[test]
fn g3_depth_partition_keeps_weak_loss_and_physical_length_scaling() {
    for alpha in [1e-24, 1e-10, 10.0, 1e4] {
        let slab = glass(alpha)
            .prepare_incoherent_slab_absorption_between(550.0, 0.02, 0.7, 1.0, 1.33)
            .unwrap();
        let scaled = glass(alpha / 1e6)
            .prepare_incoherent_slab_absorption_between(550.0, 2e4, 0.7, 1.0, 1.33)
            .unwrap();
        let full = slab.absorption_between(0.0, 0.02).unwrap();
        for count in [1, 7, 128] {
            let mut sum = [0.0; 2];
            for cell in 0..count {
                let start = 0.02 * cell as f64 / count as f64;
                let end = 0.02 * (cell + 1) as f64 / count as f64;
                let value = slab.absorption_between(start, end).unwrap();
                let rescaled = scaled.absorption_between(start * 1e6, end * 1e6).unwrap();
                for pol in 0..2 {
                    relative(value[pol], rescaled[pol], 2e-13);
                    sum[pol] += value[pol];
                }
            }
            relative(sum[0], full[0], 2e-14);
            relative(sum[1], full[1], 2e-14);
        }
        relative(full[0], slab.optical().perpendicular.absorptance, 4e-15);
        relative(full[1], slab.optical().parallel.absorptance, 4e-15);
        if alpha == 1e-24 {
            assert_eq!((-alpha * 0.02).exp(), 1.0);
            assert!(full.iter().all(|x| *x > 0.0));
        }
    }
}

#[test]
fn g1_opaque_bulk_retains_finite_skin_absorption_without_nan_or_uniform_heat() {
    let slab = glass(1e306)
        .prepare_incoherent_slab_absorption_between(550.0, 1e300, 1.0, 1.5, 1.5)
        .unwrap();
    assert_eq!(slab.optical().optical_depth, f64::INFINITY);
    let skin = slab.absorption_between(0.0, 1e-307).unwrap();
    for value in skin {
        relative(value, -(-0.1_f64).exp_m1(), 3e-15);
    }
    assert_eq!(slab.absorption_between(1.0, 1e300).unwrap(), [0.0; 2]);
    assert_eq!(slab.absorption_between(0.0, 1e300).unwrap(), [1.0; 2]);
}

#[test]
fn g3_reverse_coordinates_preserve_a_thin_interval_inside_a_much_thicker_slab() {
    for (alpha, thickness, width) in [(1e-20, 1e10, 1e-10), (1e-300, 1e300, 1e-10)] {
        let slab = glass(alpha)
            .prepare_incoherent_slab_absorption_between(550.0, thickness, 1.0, 1.5, 1.5)
            .unwrap();
        assert_eq!(thickness - width, thickness);
        assert_eq!(
            slab.absorption_between(thickness - width, thickness)
                .unwrap(),
            [0.0; 2],
            "subtracting reversed endpoints loses a physical cell"
        );
        let expected = (-alpha * thickness).exp() * -(-alpha * width).exp_m1();
        assert!(expected > 0.0);
        for value in slab.absorption_between_from_back(0.0, width).unwrap() {
            assert!(value > 0.0);
            relative(value, expected, 3e-15);
        }
        assert_eq!(
            slab.absorption_between_from_back(0.0, thickness).unwrap(),
            slab.absorption_between(0.0, thickness).unwrap()
        );
    }
}

#[test]
fn g0_transparent_tir_empty_and_invalid_intervals_have_exact_semantics() {
    let clear = glass(0.0)
        .prepare_incoherent_slab_absorption_between(550.0, 0.02, 0.7, 1.0, 1.0)
        .unwrap();
    let tir = glass(20.0)
        .prepare_incoherent_slab_absorption_between(550.0, 0.02, 0.1, 2.0, 1.0)
        .unwrap();
    let absorbing = glass(20.0)
        .prepare_incoherent_slab_absorption_between(550.0, 0.02, 0.7, 1.0, 1.0)
        .unwrap();
    for slab in [clear, tir, absorbing] {
        assert_eq!(slab.absorption_between(0.01, 0.01).unwrap(), [0.0; 2]);
        assert_eq!(
            slab.absorption_between_from_back(0.01, 0.01).unwrap(),
            [0.0; 2]
        );
        for (start, end) in [
            (-1e-9, 0.01),
            (0.01, 0.009),
            (0.0, 0.021),
            (0.03, 0.03),
            (f64::NAN, 0.02),
            (0.0, f64::INFINITY),
        ] {
            assert_eq!(
                slab.absorption_between(start, end),
                Err(DielectricError::InvalidAbsorption)
            );
            assert_eq!(
                slab.absorption_between_from_back(start, end),
                Err(DielectricError::InvalidAbsorption)
            );
        }
    }
    assert_eq!(clear.absorption_between(0.0, 0.02).unwrap(), [0.0; 2]);
    assert_eq!(tir.absorption_between(0.0, 0.02).unwrap(), [0.0; 2]);
    assert_eq!(
        clear.absorption_between_from_back(0.0, 0.02).unwrap(),
        [0.0; 2]
    );
    assert_eq!(
        tir.absorption_between_from_back(0.0, 0.02).unwrap(),
        [0.0; 2]
    );
    assert_eq!(
        glass(20.0).prepare_incoherent_slab_absorption_between(200.0, 0.02, 0.7, 1.0, 1.0),
        Err(DielectricError::InvalidWavelength)
    );
}
