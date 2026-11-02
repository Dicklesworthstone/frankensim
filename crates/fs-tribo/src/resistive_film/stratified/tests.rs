use super::*;
use crate::resistive_film::{FilmCell, FilmChannel, FilmChannelReport, FilmLimits, GapPort};

fn profile(values: &[(f64, f64)]) -> GapViscosityProfile {
    GapViscosityProfile::new(
        values
            .iter()
            .map(|&(upper_fraction, dynamic_viscosity_pa_s)| ViscosityLayer {
                upper_fraction,
                dynamic_viscosity_pa_s,
            })
            .collect(),
    )
    .unwrap()
}

fn near(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-300),
        "{actual:e} != {expected:e}"
    );
}

fn geometry() -> ResistiveFilm {
    let gap = GapPort {
        reference_m: 0.0002,
        closure: vec![1.0, -1.0],
    };
    ResistiveFilm::new(
        vec![FilmCell {
            area_m2: 0.0001,
            gap: gap.clone(),
        }],
        [0.02, 0.04]
            .into_iter()
            .map(|length_m| FilmChannel {
                from: 0,
                to: None,
                width_m: 0.01,
                length_m,
                gap: gap.clone(),
            })
            .collect(),
        2,
        0.01,
        FilmLimits {
            minimum_cell_gap_m: 0.00001,
            maximum_gap_m: 0.0005,
            maximum_pressure_pa: 10000.0,
        },
    )
    .unwrap()
}

#[test]
fn uniform_layers_recover_native_flow_tangent_and_work_exactly() {
    // G3: arbitrary subdivision cannot change a uniform Newtonian film.
    let uniform = profile(&[(0.125, 0.01), (0.55, 0.01), (1.0, 0.01)]);
    assert_eq!(uniform.hydraulic_equivalent_viscosity_pa_s(), 0.01);
    assert_eq!(uniform.neutral_plane_fraction(), 0.5);
    near(
        uniform.layer_pressure_work_fractions().iter().sum(),
        1.0,
        2e-15,
    );
    for alpha in [0.0, 2e-4] {
        let native = geometry().with_barus_viscosity(alpha).unwrap();
        let binding = StratifiedFilm::new(&native, vec![uniform.clone(), uniform.clone()]).unwrap();
        assert_eq!(binding.profiles()[0].layers(), uniform.layers());
        assert_eq!(binding.film().channel_viscosities_pa_s(), &[0.01, 0.01]);
        for speed in [-0.001, 0.001] {
            let q = [0.00001, 0.0];
            let v = [speed, 0.0];
            let mut force = [[0.0; 2]; 2];
            let mut pressure = [[0.0]; 2];
            let mut channels = [[FilmChannelReport::default(); 2]; 2];
            let mut tangent = [[0.0; 2]; 2];
            for (i, film) in [&native, binding.film()].into_iter().enumerate() {
                film.evaluate_with_channels_into(
                    &q,
                    &v,
                    &mut force[i],
                    &mut pressure[i],
                    &mut channels[i],
                )
                .unwrap();
                film.tangent_into(&q, &v, &[0.000001, 0.0], &[0.0001, 0.0], &mut tangent[i])
                    .unwrap();
            }
            assert_eq!(force[0], force[1]);
            assert_eq!(pressure[0], pressure[1]);
            assert_eq!(channels[0], channels[1]);
            assert_eq!(tangent[0], tangent[1]);
        }
        let (compiled, retained) = binding.into_parts();
        assert_eq!(retained, vec![uniform.clone(), uniform.clone()]);
        assert_eq!(compiled.pressure_viscosity_coefficient_per_pa(), alpha);
    }
}

#[test]
fn two_equal_thickness_layers_match_closed_form_shear_plane_and_pressure_work() {
    // G2: direct two-layer no-slip solution with mu_lower=1, mu_upper=4.
    // For general a,b, mu_h = 8*a*b*(a+b)/(a*a+14*a*b+b*b).
    let p = profile(&[(0.5, 1.0), (1.0, 4.0)]);
    near(p.hydraulic_equivalent_viscosity_pa_s(), 160.0 / 73.0, 2e-15);
    near(p.neutral_plane_fraction(), 7.0 / 20.0, 2e-15);
    near(p.layer_pressure_work_fractions()[0], 148.0 / 365.0, 2e-15);
    near(p.layer_pressure_work_fractions()[1], 217.0 / 365.0, 2e-15);
    // Neither volume weighting nor arithmetic/harmonic bulk viscosity gives
    // the exact pressure-flow resistance of these perfectly bonded layers.
    assert_ne!(p.hydraulic_equivalent_viscosity_pa_s(), 2.5);
    assert_ne!(p.hydraulic_equivalent_viscosity_pa_s(), 1.6);
}

#[test]
fn symmetric_wall_layers_have_the_exact_cubic_mobility_weights() {
    // G2: for a symmetric core of thickness c, 1/mu_h =
    // (1-c^3)/mu_wall + c^3/mu_core. Wall location matters to shear work.
    let p = profile(&[(0.1, 0.001), (0.9, 0.1), (1.0, 0.001)]);
    let core = 0.8_f64.powi(3) / 0.1;
    let walls = (1.0 - 0.8_f64.powi(3)) / 0.001;
    near(
        p.hydraulic_equivalent_viscosity_pa_s(),
        1.0 / (walls + core),
        3e-15,
    );
    near(p.neutral_plane_fraction(), 0.5, 2e-15);
    near(
        p.layer_pressure_work_fractions()[0],
        0.5 * walls / (walls + core),
        3e-15,
    );
    near(
        p.layer_pressure_work_fractions()[1],
        core / (walls + core),
        3e-15,
    );
    near(
        p.layer_pressure_work_fractions()[2],
        0.5 * walls / (walls + core),
        3e-15,
    );
    assert!(p.hydraulic_equivalent_viscosity_pa_s() < (0.2 * 0.001 + 0.8 * 0.1) / 20.0);
}

#[test]
fn independent_piecewise_velocity_quadrature_recovers_flux_and_each_layers_work() {
    // G2: integrate the continuous quadratic velocity separately, then use
    // two-point Gaussian quadrature. No central-variance implementation is reused.
    let p = profile(&[(0.13, 0.003), (0.4, 0.09), (0.91, 0.012), (1.0, 0.2)]);
    let mut lo = 0.0;
    let mut a = 0.0;
    let mut b = 0.0;
    for layer in p.layers() {
        a += (layer.upper_fraction - lo) / layer.dynamic_viscosity_pa_s;
        b += (layer.upper_fraction.powi(2) - lo * lo) / (2.0 * layer.dynamic_viscosity_pa_s);
        lo = layer.upper_fraction;
    }
    let z0 = b / a;
    let mut velocity_at_lower = 0.0;
    let mut flux = 0.0;
    let mut work = Vec::new();
    lo = 0.0;
    for layer in p.layers() {
        let hi = layer.upper_fraction;
        let width = hi - lo;
        let half = 0.5 * width;
        let points = [
            lo + half - half / 3.0_f64.sqrt(),
            lo + half + half / 3.0_f64.sqrt(),
        ];
        let velocity = |z: f64| {
            velocity_at_lower
                + (z0 * (z - lo) - 0.5 * (z * z - lo * lo)) / layer.dynamic_viscosity_pa_s
        };
        flux += half * (velocity(points[0]) + velocity(points[1]));
        work.push(
            half * ((points[0] - z0).powi(2) + (points[1] - z0).powi(2))
                / layer.dynamic_viscosity_pa_s,
        );
        velocity_at_lower +=
            (z0 * width - 0.5 * (hi * hi - lo * lo)) / layer.dynamic_viscosity_pa_s;
        lo = hi;
    }
    assert!(velocity_at_lower.abs() < 1e-14);
    let dissipated = work.iter().sum::<f64>();
    near(dissipated, flux, 2e-14);
    near(
        p.hydraulic_equivalent_viscosity_pa_s(),
        1.0 / (12.0 * flux),
        2e-14,
    );
    near(p.neutral_plane_fraction(), z0, 2e-15);
    for (actual, expected) in p.layer_pressure_work_fractions().iter().zip(work) {
        near(*actual, expected / dissipated, 2e-14);
    }
}

#[test]
fn reflection_and_uniform_viscosity_scaling_preserve_the_physical_profile() {
    // G3: swap the walls; uniform viscosity-unit scaling changes only resistance.
    let original = profile(&[(0.125, 0.003), (0.375, 0.09), (0.875, 0.012), (1.0, 0.2)]);
    let layers = original.layers();
    let reflected = GapViscosityProfile::new(
        (0..layers.len())
            .rev()
            .map(|i| ViscosityLayer {
                upper_fraction: 1.0
                    - if i == 0 {
                        0.0
                    } else {
                        layers[i - 1].upper_fraction
                    },
                dynamic_viscosity_pa_s: layers[i].dynamic_viscosity_pa_s,
            })
            .collect(),
    )
    .unwrap();
    near(
        reflected.hydraulic_equivalent_viscosity_pa_s(),
        original.hydraulic_equivalent_viscosity_pa_s(),
        3e-15,
    );
    near(
        reflected.neutral_plane_fraction(),
        1.0 - original.neutral_plane_fraction(),
        3e-15,
    );
    for (a, b) in reflected
        .layer_pressure_work_fractions()
        .iter()
        .zip(original.layer_pressure_work_fractions().iter().rev())
    {
        near(*a, *b, 3e-15);
    }
    for scale in [1e-180, 1e180] {
        let scaled = GapViscosityProfile::new(
            layers
                .iter()
                .map(|layer| ViscosityLayer {
                    dynamic_viscosity_pa_s: layer.dynamic_viscosity_pa_s * scale,
                    ..*layer
                })
                .collect(),
        )
        .unwrap();
        near(
            scaled.hydraulic_equivalent_viscosity_pa_s(),
            original.hydraulic_equivalent_viscosity_pa_s() * scale,
            3e-14,
        );
        near(
            scaled.neutral_plane_fraction(),
            original.neutral_plane_fraction(),
            3e-15,
        );
        for (a, b) in scaled
            .layer_pressure_work_fractions()
            .iter()
            .zip(original.layer_pressure_work_fractions())
        {
            near(*a, *b, 3e-14);
        }
    }
}

#[test]
fn narrow_mobile_core_at_extreme_contrast_keeps_its_finite_pressure_mobility() {
    // Raw C-B^2/A loses the O(thickness^3) mobility beneath O(thickness) moments.
    // The exact symmetric-core oracle remains representable at this contrast.
    let thickness = 2.0_f64.powi(-30);
    let p = profile(&[
        (0.5 - 0.5 * thickness, 1e40),
        (0.5 + 0.5 * thickness, 1.0),
        (1.0, 1e40),
    ]);
    let core = thickness.powi(3);
    let walls = (1.0 - core) / 1e40;
    near(
        p.hydraulic_equivalent_viscosity_pa_s(),
        1.0 / (core + walls),
        3e-15,
    );
    assert_eq!(p.neutral_plane_fraction(), 0.5);
    near(
        p.layer_pressure_work_fractions()[1],
        core / (core + walls),
        3e-15,
    );
    near(
        p.layer_pressure_work_fractions()[0],
        0.5 * walls / (core + walls),
        3e-15,
    );
    near(
        p.layer_pressure_work_fractions()[2],
        0.5 * walls / (core + walls),
        3e-15,
    );
}

#[test]
fn heterogeneous_profiles_drive_native_flow_and_barus_pressure_work() {
    let profiles = vec![
        profile(&[(0.5, 0.005), (1.0, 0.02)]),
        profile(&[(0.25, 0.002), (0.75, 0.02), (1.0, 0.002)]),
    ];
    let equivalents = [0.005 * 160.0 / 73.0, 1.0 / (0.875 / 0.002 + 0.125 / 0.02)];
    for alpha in [0.0, 2e-4] {
        let native = geometry().with_barus_viscosity(alpha).unwrap();
        let binding = StratifiedFilm::new(&native, profiles.clone()).unwrap();
        for speed in [-0.001, 0.001] {
            let q = [0.00001, 0.0];
            let mut force = [0.0; 2];
            let mut pressure = [0.0];
            let mut channels = [FilmChannelReport::default(); 2];
            let report = binding
                .film()
                .evaluate_with_channels_into(
                    &q,
                    &[speed, 0.0],
                    &mut force,
                    &mut pressure,
                    &mut channels,
                )
                .unwrap();
            let conductance =
                [0.02, 0.04].map(|length| 0.01 * 0.00019_f64.powi(3) / (12.0 * length));
            let g = [
                conductance[0] / equivalents[0],
                conductance[1] / equivalents[1],
            ];
            let psi = 0.0001 * speed / (g[0] + g[1]);
            let expected_pressure = if alpha == 0.0 {
                psi
            } else {
                -(-alpha * psi).ln_1p() / alpha
            };
            near(pressure[0], expected_pressure, 4e-15);
            near(force[0], 0.0001 * expected_pressure, 4e-15);
            assert_eq!(force[1], -force[0]);
            let mut layer_power = 0.0;
            for i in 0..2 {
                near(channels[i].volume_flow_m3_s, g[i] * psi, 4e-15);
                near(
                    channels[i].dissipated_power_w,
                    g[i] * psi * expected_pressure,
                    6e-15,
                );
                assert!(channels[i].dissipated_power_w > 0.0);
                layer_power += binding.profiles()[i]
                    .layer_pressure_work_fractions()
                    .iter()
                    .map(|fraction| fraction * channels[i].dissipated_power_w)
                    .sum::<f64>();
            }
            near(report.dissipated_power_w, layer_power, 4e-15);
            near(
                report.dissipated_power_w,
                0.0001 * speed * expected_pressure,
                6e-15,
            );
        }
        // Compiling equivalents never edits the original material coefficients.
        assert_eq!(native.channel_viscosities_pa_s(), &[0.01, 0.01]);
    }
}

#[test]
fn incomplete_profiles_and_unrepresentable_layers_refuse_without_floors() {
    let valid = ViscosityLayer {
        upper_fraction: 1.0,
        dynamic_viscosity_pa_s: 0.01,
    };
    assert!(GapViscosityProfile::new(vec![]).is_err());
    assert!(GapViscosityProfile::new(vec![valid; MAX_LAYERS + 1]).is_err());
    for bound in [0.0, -0.1, 0.9, 1.1, f64::NAN, f64::INFINITY] {
        assert!(
            GapViscosityProfile::new(vec![ViscosityLayer {
                upper_fraction: bound,
                ..valid
            }])
            .is_err()
        );
    }
    for viscosity in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            GapViscosityProfile::new(vec![ViscosityLayer {
                dynamic_viscosity_pa_s: viscosity,
                ..valid
            }])
            .is_err()
        );
    }
    for middle in [0.49, 0.5] {
        assert!(
            GapViscosityProfile::new(vec![
                ViscosityLayer {
                    upper_fraction: 0.5,
                    ..valid
                },
                ViscosityLayer {
                    upper_fraction: middle,
                    ..valid
                },
                valid,
            ])
            .is_err()
        );
    }
    assert!(
        GapViscosityProfile::new(vec![
            ViscosityLayer {
                upper_fraction: 0.5,
                dynamic_viscosity_pa_s: f64::MIN_POSITIVE
            },
            ViscosityLayer {
                upper_fraction: 1.0,
                dynamic_viscosity_pa_s: f64::MAX
            },
        ])
        .is_err()
    );
    let maximal = GapViscosityProfile::new(
        (1..=MAX_LAYERS)
            .map(|i| ViscosityLayer {
                upper_fraction: i as f64 / MAX_LAYERS as f64,
                ..valid
            })
            .collect(),
    )
    .unwrap();
    assert_eq!(maximal.layers().len(), MAX_LAYERS);
    assert_eq!(maximal.hydraulic_equivalent_viscosity_pa_s(), 0.01);
    near(
        maximal.layer_pressure_work_fractions().iter().sum(),
        1.0,
        3e-15,
    );
    let native = geometry();
    for count in [0, 1, 3] {
        assert!(StratifiedFilm::new(&native, vec![maximal.clone(); count]).is_err());
    }
    assert_eq!(native.channel_viscosities_pa_s(), &[0.01, 0.01]);
}
