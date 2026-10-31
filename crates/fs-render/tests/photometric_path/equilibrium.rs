//! Isothermal spectral equilibrium through walls and an emitting, scattering gas.
//! This manufactured reference checks the represented visible-band transport;
//! it supplies no measured-material calibration or production error certificate.

use super::*;
use fs_render::photometry::IlluminanceProbe;
use fs_render::photometry::area::{EmissionSidedness, PhotometricAreaEmitter, PhotometricEmitters};
use fs_render::photometry::path::{IlluminanceSampling, PathIlluminanceJob};
use fs_render::spectral::SpectralReflectance;

use super::spectral_reference::cie_reference_xyz;

const TEMPERATURE_K: f64 = 3000.0;
const DEPTH: u32 = 24;
const SAMPLES: u32 = 32_768;

fn planck_reference(wavelength_nm: f64) -> f64 {
    // Original SI constants, independent of the production Planck evaluator.
    // The final factor converts W/(m² sr m) to W/(m² sr nm).
    let h = 6.626_070_15e-34;
    let c = 299_792_458.0;
    let k = 1.380_649e-23;
    let wavelength_m = wavelength_nm * 1.0e-9;
    2.0 * h * c * c / (wavelength_m.powi(5) * (h * c / (wavelength_m * k * TEMPERATURE_K)).exp_m1())
        * 1.0e-9
}

fn equilibrium_match(actual: f64, expected: f64, estimated_se: f64, tail_fraction: f64) {
    assert!(estimated_se.is_finite() && estimated_se >= 0.0);
    // A fixed IID regression rule, not a confidence interval. The small tail
    // allowance belongs only to this manufactured transport operator; it is
    // separate from estimated sampling variation and arithmetic tolerance.
    let tolerance = 6.0 * estimated_se + (tail_fraction + 2.0e-10) * expected;
    assert!(
        actual.is_finite() && (actual - expected).abs() <= tolerance,
        "equilibrium: observed={actual}, expected={expected}, estimated_se={estimated_se}, \
         model_tail_fraction={tail_fraction}, tolerance={tolerance}"
    );
    // An inflated variance estimate must not hide a substantial energy error.
    close(
        actual,
        expected,
        0.03 + tail_fraction,
        "isothermal equilibrium",
    );
}

#[test]
fn g1_isothermal_spectral_enclosure_preserves_lte_equilibrium() {
    with_cx(|cx, _| {
        assert_eq!(SPECTRAL_BINS, 80);
        let mut expected_xyz = [0.0; 3];
        let mut wall_bins = [0.0; SPECTRAL_BINS];
        let mut wall_luminance = 0.0;
        for (bin, wall) in wall_bins.iter_mut().enumerate() {
            // Explicit 80-bin midpoint rule: qualify this declared spectral
            // representation, without using production wavelength/observer
            // helpers or claiming exact integration of the continuous band.
            let wavelength_nm = 382.5 + 5.0 * bin as f64;
            let radiance = planck_reference(wavelength_nm);
            let observer = cie_reference_xyz(wavelength_nm);
            let rho = 0.2 + 0.4 * (wavelength_nm - 380.0) / 400.0;
            *wall = (1.0 - rho) * radiance;
            wall_luminance += 683.0 * observer[1] * *wall * 5.0;
            for channel in 0..3 {
                expected_xyz[channel] += 683.0 * observer[channel] * radiance * 5.0;
            }
        }
        let mut scene = empty_scene([0.13, -0.17, 0.09], [0.0, 0.0, 1.0]);
        let mut enclosure = super::area::diffuse_emitting_cube(0.0);
        enclosure.material = Material::SpectralDiffuse {
            reflectance: SpectralReflectance::table(&[(380.0, 0.2), (780.0, 0.6)]).unwrap(),
        };
        scene.primitives.push(enclosure);
        let surfaces = [PhotometricAreaEmitter::try_new(
            0,
            // Lambertian exitance is pi*Le; six inward 2x2 m faces emit.
            24.0 * PI * wall_luminance,
            1.0,
            SourceSpectrum::try_from_bins(wall_bins).unwrap(),
            EmissionSidedness::Front,
        )
        .unwrap()];
        let emitters = PhotometricEmitters {
            points: &[],
            surfaces: &surfaces,
        };
        let temperatures = [TEMPERATURE_K; 2];
        let probes = [
            IlluminanceProbe::try_new(Point3::new(0.0, 0.0, -0.4), Vec3::new(0.0, 0.0, 1.0))
                .unwrap(),
            IlluminanceProbe::try_new(Point3::new(0.45, -0.3, 0.35), Vec3::new(0.6, 0.0, 0.8))
                .unwrap(),
        ];

        // L_lambda=B_lambda is an exact all-order solution: Le+rho*B=B at
        // every wall, absorption is balanced by kappa_a*B emission, and a
        // normalized reciprocal phase preserves isotropic B under scattering.
        // Both rho and kappa_s/(kappa_a+kappa_s) are at most q=0.6. For any
        // flight transmittance T, the continuation operator satisfies
        // K B <= [q*(1-T)+q*T]*B = q*B, irrespective of the density layers.
        // Thus omitted ideal transport orders are <=q^D*B. Use D-1 here to
        // conservatively cover camera/probe terminal-source conventions; MIS
        // can include additional terminal source terms. This constructed
        // model bound is NOT a bound on floating-point or general-scene error.
        let tail_fraction = 0.6_f64.powi(DEPTH as i32 - 1);
        assert!(tail_fraction < 8.0e-6);
        for (density, strategy) in [
            ([0.5, 0.5], DirectStrategy::BsdfOnly),
            ([0.25, 2.0], DirectStrategy::Mis),
        ] {
            let medium = ScatteringGrid::try_new(
                cx,
                AbsorbingGrid::try_new(
                    cx,
                    VolumeGrid::new([1, 1, 2], &density, [-1.0; 3], [2.0, 2.0, 1.0]),
                    MassAbsorptionSpectrum::table(&[(380.0, 0.5), (780.0, 1.2)]).unwrap(),
                    Some(&temperatures),
                    medium_budget(),
                )
                .unwrap(),
                MassScatteringSpectrum::table(&[(380.0, 0.75), (780.0, 0.3)]).unwrap(),
                PhaseFunction::tabulated(&[
                    (-1.0, 0.15),
                    (-0.6, 0.22),
                    (0.1, 0.55),
                    (0.55, 1.8),
                    (1.0, 4.0),
                ])
                .unwrap(),
            )
            .unwrap();
            let mut image = PathLuminanceJob::new_with_emitters(
                &scene,
                emitters,
                PathMedium::Scattering(&medium),
                TransportOptions::default(),
                Settings {
                    sampler: Sampler::Iid,
                    strategy,
                    ..settings(DEPTH, SAMPLES)
                },
                cx,
            )
            .unwrap();
            image.advance(cx, usize::MAX).unwrap();
            assert!(image.is_complete());
            let value = image.pixel_estimate(0).unwrap();
            assert_eq!(value.samples, SAMPLES);
            for (actual, expected) in value.xyz.into_iter().zip(expected_xyz) {
                close(
                    actual,
                    expected,
                    0.03 + tail_fraction,
                    "Planck equilibrium in independent original-CIE XYZ",
                );
            }
            equilibrium_match(
                value.luminance_cd_m2,
                expected_xyz[1],
                value.estimated_standard_error_cd_m2.unwrap(),
                tail_fraction,
            );
            assert_eq!(value.numerical.kind, fs_evidence::NumericalKind::NoClaim);
            assert!(value.total_error_bound_cd_m2.is_none());

            let mut illuminance = PathIlluminanceJob::new_with_emitters(
                &scene,
                emitters,
                PathMedium::Scattering(&medium),
                &probes,
                TransportOptions::default(),
                IlluminanceSampling {
                    samples_per_probe: SAMPLES,
                    max_depth: DEPTH,
                    sampler: Sampler::Iid,
                    seed: 173,
                },
                cx,
            )
            .unwrap();
            illuminance.advance(cx, usize::MAX).unwrap();
            assert!(illuminance.is_complete());
            for value in illuminance.estimates() {
                let value = value.unwrap();
                assert_eq!(value.samples, SAMPLES);
                assert_eq!(value.direct_lux, 0.0);
                // A complete receiving hemisphere of constant radiance has
                // E=pi*L, independent of probe position and inclination.
                equilibrium_match(
                    value.illuminance_lux,
                    PI * expected_xyz[1],
                    value.estimated_standard_error_lux.unwrap(),
                    tail_fraction,
                );
                assert_eq!(value.numerical.kind, fs_evidence::NumericalKind::NoClaim);
                assert!(value.total_error_bound_lux.is_none());
            }
        }
    });
}
