use super::*;
use crate::thermal_surface::planck_radiance_w_m2_sr_nm;
use crate::volumes::VolumeGrid;
use crate::volumes::transport::AbsorbingGrid;
use crate::volumes::transport::scattering::ScatteringGrid;
use crate::volumes::transport::scattering::tests::{budget, close, with_cx};
use fs_blake3::ContentHash;
use fs_geom::{Point3, Vec3};

fn mixture_budget() -> OpticalMixtureBudget {
    OpticalMixtureBudget { max_components: 8 }
}

#[test]
fn g1_linear_absorption_sum_retains_original_breakpoints_and_sources() {
    with_cx(|cx, _| {
        let a = MassAbsorptionSpectrum::source_table(
            &[(380.0, 1.0), (470.0, 5.0), (780.0, 0.0)],
            ContentHash([1; 32]),
        )
        .unwrap();
        let b = MassAbsorptionSpectrum::source_table(
            &[(380.0, 2.0), (610.0, 0.0), (780.0, 3.0)],
            ContentHash([2; 32]),
        )
        .unwrap();
        let sum = MassAbsorptionSpectrum::weighted_sum(cx, &[(0.2, a), (0.8, b)], mixture_budget())
            .unwrap();
        assert_eq!(sum.count, 4);
        assert_eq!(sum.knots[1].0, 470.0);
        assert_eq!(sum.knots[2].0, 610.0);
        for i in 0..=80 {
            let wavelength = 380.0 + 5.0 * i as f64;
            let first = if wavelength <= 470.0 {
                1.0 + 4.0 * (wavelength - 380.0) / 90.0
            } else {
                5.0 * (780.0 - wavelength) / 310.0
            };
            let second = if wavelength <= 610.0 {
                2.0 * (610.0 - wavelength) / 230.0
            } else {
                3.0 * (wavelength - 610.0) / 170.0
            };
            close(
                sum.at(wavelength).unwrap(),
                0.2 * first + 0.8 * second,
                3e-14,
            );
        }
        assert_eq!(
            MassAbsorptionSpectrum::weighted_sum(cx, &[(1.0, a)], mixture_budget()).unwrap(),
            a
        );
        let mut other_source = a;
        other_source.source_identity = Some(ContentHash([3; 32]));
        let changed = MassAbsorptionSpectrum::weighted_sum(
            cx,
            &[(0.2, other_source), (0.8, b)],
            mixture_budget(),
        )
        .unwrap();
        assert_ne!(sum.source_identity(), changed.source_identity());
        assert_eq!(sum.at(520.0), changed.at(520.0));
    });
}

#[test]
fn g1_linear_plus_molecular_scattering_stays_analytic_between_knots() {
    with_cx(|cx, _| {
        let gas = MassScatteringSpectrum::rayleigh(550.0, 0.5).unwrap();
        let second_gas = MassScatteringSpectrum::rayleigh(650.0, 0.2).unwrap();
        let particles =
            MassScatteringSpectrum::table(&[(380.0, 1.0), (550.0, 4.0), (780.0, 2.0)]).unwrap();
        let mixed = MassScatteringSpectrum::weighted_sum(
            cx,
            &[(0.3, gas), (0.05, second_gas), (0.7, particles)],
            mixture_budget(),
        )
        .unwrap();
        assert_eq!(mixed.spectrum.count, 3);
        assert!(mixed.rayleigh_addend.is_some());
        for wavelength in [380.0_f64, 419.0, 487.0, 550.0, 703.0, 780.0] {
            let particle = if wavelength <= 550.0 {
                1.0 + 3.0 * (wavelength - 380.0) / 170.0
            } else {
                4.0 - 2.0 * (wavelength - 550.0) / 230.0
            };
            let expected = 0.15 * (550.0 / wavelength).powi(4)
                + 0.01 * (650.0 / wavelength).powi(4)
                + 0.7 * particle;
            close(mixed.at(wavelength).unwrap(), expected, 4e-14);
        }
        assert_eq!(
            MassScatteringSpectrum::weighted_sum(cx, &[(1.0, gas)], mixture_budget()).unwrap(),
            gas
        );
        let nested = MassScatteringSpectrum::weighted_sum(
            cx,
            &[(0.4, mixed), (0.6, particles)],
            mixture_budget(),
        )
        .unwrap();
        close(
            nested.at(431.0).unwrap(),
            0.4 * mixed.at(431.0).unwrap() + 0.6 * particles.at(431.0).unwrap(),
            4e-14,
        );
    });
}

#[test]
fn g1_analytic_molecular_sum_reaches_native_flight_emission_and_heat_paths() {
    with_cx(|cx, _| {
        let components = [
            ScatteringMixtureComponent {
                mass_weight: 0.25,
                scattering: MassScatteringSpectrum::rayleigh(550.0, 0.4).unwrap(),
                phase: PhaseFunction::Rayleigh,
            },
            ScatteringMixtureComponent {
                mass_weight: 0.75,
                scattering: MassScatteringSpectrum::rayleigh(600.0, 1.2).unwrap(),
                phase: PhaseFunction::rayleigh_depolarized(0.15).unwrap(),
            },
        ];
        let optics = ScatteringOptics::weighted_sum(cx, &components, mixture_budget()).unwrap();
        // The retained LINEAR part is zero: admission must see the analytic addend.
        assert!(
            optics
                .scattering
                .knots()
                .iter()
                .all(|(_, value)| *value == 0.0)
        );
        let density = [3.0];
        let temperature = [1600.0];
        let absorption = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1; 3], &density, [0.0; 3], [1.0, 1.0, 0.6]),
            MassAbsorptionSpectrum::constant(0.2).unwrap(),
            Some(&temperature),
            budget(8),
        )
        .unwrap();
        let medium =
            ScatteringGrid::try_new(cx, absorption, optics.scattering, optics.phase).unwrap();
        assert!(medium.has_scattering());
        let wavelengths = [405.0_f64, 550.0, 750.0];
        let origin = Point3::new(0.5, 0.5, -0.2);
        let direction = Vec3::new(0.0, 0.0, 1.0);
        let transport = medium
            .transport(cx, origin, direction, 1.0, &wavelengths)
            .unwrap();
        let powers = [2.0, 3.0, 4.0];
        let heat = medium
            .deposit_beam_power(cx, origin, direction, 1.0, &wavelengths, &powers, 1)
            .unwrap();
        let scattering = wavelengths.map(|wavelength| {
            0.1 * (550.0 / wavelength).powi(4) + 0.9 * (600.0 / wavelength).powi(4)
        });
        for lane in 0..3 {
            let total = 0.2 + scattering[lane];
            let transmittance = (-1.8 * total).exp();
            let absorbed = (1.0 - transmittance) * 0.2 / total;
            close(transport.transmittance[lane], transmittance, 2e-13);
            let radiance = planck_radiance_w_m2_sr_nm(1600.0, wavelengths[lane]).unwrap();
            close(
                transport.emission_w_m2_sr_nm[lane],
                radiance * absorbed,
                3e-13,
            );
            close(heat.absorbed_power_w[lane], powers[lane] * absorbed, 3e-13);
            close(
                heat.transmitted_power_w[lane]
                    + heat.absorbed_power_w[lane]
                    + heat.scattered_power_w[lane],
                powers[lane],
                3e-13,
            );
        }
        let flight = medium
            .sample_segment(cx, origin, direction, 1.0, &wavelengths, 0.1, 0.2)
            .unwrap();
        let event = flight.event.unwrap();
        assert_eq!(event.cell_index, 0);
        let distance = -0.2_f64.ln() / (3.0 * (0.2 + scattering[0]));
        close(event.distance_m, 0.2 + distance, 3e-13);
        let transmission = scattering.map(|s| (-3.0 * distance * (0.2 + s)).exp());
        let proposal: f64 = (0..3)
            .map(|lane| transmission[lane] * (0.2 + scattering[lane]) / 3.0)
            .sum();
        for lane in 0..3 {
            close(
                flight.weight[lane],
                transmission[lane] * scattering[lane] / proposal,
                3e-13,
            );
        }
    });
}

#[test]
fn g1_non_scattering_mixture_keeps_incident_power_separate_from_lte_source() {
    with_cx(|cx, _| {
        let component = ScatteringMixtureComponent {
            mass_weight: 0.5,
            scattering: MassScatteringSpectrum::constant(0.0).unwrap(),
            phase: PhaseFunction::Rayleigh,
        };
        let optics = ScatteringOptics::weighted_sum(cx, &[component], mixture_budget()).unwrap();
        let density = [1.0];
        let temperature = [1800.0];
        let absorption = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1; 3], &density, [0.0; 3], [1.0; 3]),
            MassAbsorptionSpectrum::constant(0.4).unwrap(),
            Some(&temperature),
            budget(8),
        )
        .unwrap();
        let medium =
            ScatteringGrid::try_new(cx, absorption, optics.scattering, optics.phase).unwrap();
        assert!(!medium.has_scattering());
        let wavelengths = [550.0];
        let coefficients = medium.coefficients(&wavelengths).unwrap();
        let sample = |include_emission| {
            medium
                .sample_segment_prepared(
                    cx,
                    Point3::new(0.5, 0.5, 0.0),
                    Vec3::new(0.0, 0.0, 1.0),
                    1.0,
                    &wavelengths,
                    0.2,
                    0.3,
                    (&coefficients.0, &coefficients.1, &coefficients.2),
                    include_emission,
                )
                .unwrap()
        };
        let incident = sample(false);
        let emitting = sample(true);
        assert_eq!(incident.weight, emitting.weight);
        assert_eq!(incident.emission_w_m2_sr_nm, [0.0]);
        assert!(emitting.emission_w_m2_sr_nm[0] > 0.0);
    });
}

#[test]
fn g0_exact_knot_capacity_weights_domains_and_cancellation_refuse_explicitly() {
    with_cx(|cx, gate| {
        let table: Vec<_> = (0..16)
            .map(|i| (380.0 + 400.0 * i as f64 / 15.0, 1.0 + i as f64))
            .collect();
        let full = MassAbsorptionSpectrum::table(&table).unwrap();
        let extra =
            MassAbsorptionSpectrum::table(&[(380.0, 0.0), (393.0, 1.0), (780.0, 0.0)]).unwrap();
        assert!(
            MassAbsorptionSpectrum::weighted_sum(
                cx,
                &[(1.0, full), (0.0, extra)],
                mixture_budget()
            )
            .is_ok()
        );
        assert_eq!(
            MassAbsorptionSpectrum::weighted_sum(
                cx,
                &[(1.0, full), (1.0, extra)],
                mixture_budget()
            ),
            Err(VolumeTransportError::Budget)
        );
        assert_eq!(
            MassAbsorptionSpectrum::weighted_sum(cx, &[], mixture_budget()),
            Err(VolumeTransportError::Mixture)
        );
        for weight in [-1.0, f64::NAN, f64::INFINITY] {
            assert_eq!(
                MassAbsorptionSpectrum::weighted_sum(cx, &[(weight, full)], mixture_budget()),
                Err(VolumeTransportError::Mixture)
            );
        }
        let infrared = MassAbsorptionSpectrum::source_table(
            &[(500.0, 0.0), (2500.0, 1.0)],
            ContentHash([9; 32]),
        )
        .unwrap();
        assert_eq!(
            MassAbsorptionSpectrum::weighted_sum(
                cx,
                &[(1.0, full), (0.0, infrared)],
                mixture_budget()
            ),
            Err(VolumeTransportError::Wavelength)
        );
        let zero =
            MassAbsorptionSpectrum::weighted_sum(cx, &[(0.0, full)], mixture_budget()).unwrap();
        assert_eq!(zero.at(500.0), Ok(0.0));
        gate.request();
        assert_eq!(
            MassScatteringSpectrum::weighted_sum(
                cx,
                &[(1.0, MassScatteringSpectrum::constant(1.0).unwrap())],
                mixture_budget()
            ),
            Err(VolumeTransportError::Cancelled)
        );
    });
}
