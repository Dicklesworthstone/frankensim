use super::super::{
    ParticleSizeDistribution,
    tests::{angular_integrals, close, sphere_for_x, with_cx},
};
use super::*;

#[test]
fn g0_single_bin_reduces_to_the_original_sphere_including_native_spectra() {
    with_cx(|cx, _| {
        let sphere = sphere_for_x(1.3, 1.5, 0.01);
        let population = MiePopulation::try_new(&[(sphere, 17.0)]).unwrap();
        let single = sphere
            .at_wavelength(cx, 550.0, MieBudget::default())
            .unwrap();
        let mixed = population
            .at_wavelength(cx, 550.0, MiePopulationBudget::default())
            .unwrap();
        assert_eq!(mixed.cross_sections(), single.cross_sections());
        assert_eq!(
            mixed.mean_particle_mass_kg(),
            sphere.particle_mass_kg().unwrap()
        );
        assert_eq!(mixed.multipoles(), single.multipoles());
        assert_eq!(mixed.total_multipoles(), single.multipoles());
        for mu in [-1.0, -0.3, 0.0, 0.8, 1.0] {
            assert_eq!(mixed.eval_cos(mu), single.eval_cos(mu));
        }
        assert_eq!(
            mixed.phase_table(cx, 0.003, 4096),
            single.phase_table(cx, 0.003, 4096)
        );
        let wavelengths = [380.0, 510.0, 650.0, 780.0];
        assert_eq!(
            population.visible_spectrum(
                cx,
                &wavelengths,
                MiePopulationBudget::default(),
                0.003,
                4096
            ),
            sphere.visible_spectrum(cx, &wavelengths, MieBudget::default(), 0.003, 4096),
        );
    });
}

#[test]
fn g1_number_mass_and_scattered_power_use_their_distinct_physical_weights() {
    with_cx(|cx, _| {
        let small = sphere_for_x(0.3, 1.5, 0.04);
        let large = MieSphere::try_new(7.0 * small.radius_m(), 1.8, 0.2, 1.0, 4300.0).unwrap();
        let population = MiePopulation::try_new(&[(small, 8.0), (large, 2.0)]).unwrap();
        let solution = population
            .at_wavelength(cx, 550.0, MiePopulationBudget::default())
            .unwrap();
        let a = small
            .at_wavelength(cx, 550.0, MieBudget::default())
            .unwrap();
        let b = large
            .at_wavelength(cx, 550.0, MieBudget::default())
            .unwrap();
        let (ca, cb) = (a.cross_sections(), b.cross_sections());
        let mass =
            0.8 * small.particle_mass_kg().unwrap() + 0.2 * large.particle_mass_kg().unwrap();
        let scattering = 0.8 * ca.scattering_m2 + 0.2 * cb.scattering_m2;
        let absorption = 0.8 * ca.absorption_m2 + 0.2 * cb.absorption_m2;
        let actual = solution.cross_sections();
        close(actual.scattering_m2 / scattering, 1.0, 3e-15);
        close(actual.absorption_m2 / absorption, 1.0, 3e-15);
        close(
            actual.mass_scattering_m2_kg / (scattering / mass),
            1.0,
            3e-15,
        );
        close(
            actual.mass_absorption_m2_kg / (absorption / mass),
            1.0,
            3e-15,
        );
        assert!(
            (actual.mass_scattering_m2_kg
                / (0.8 * ca.mass_scattering_m2_kg + 0.2 * cb.mass_scattering_m2_kg)
                - 1.0)
                .abs()
                > 0.1
        );
        close(
            actual.mean_cosine,
            (0.8 * ca.scattering_m2 * ca.mean_cosine + 0.2 * cb.scattering_m2 * cb.mean_cosine)
                / scattering,
            3e-15,
        );
        for mu in [-1.0, -0.37, 0.0, 0.6, 1.0] {
            let expected = (0.8 * ca.scattering_m2 * a.eval_cos(mu).unwrap()
                + 0.2 * cb.scattering_m2 * b.eval_cos(mu).unwrap())
                / scattering;
            close(solution.eval_cos(mu).unwrap() / expected, 1.0, 3e-15);
        }
        // Changing concentration preserves coefficients per particle/mass.
        let rescaled = MiePopulation::try_new(&[(small, 8e300), (large, 2e300)]).unwrap();
        assert_eq!(rescaled, population);
        let (normalization, moment) = angular_integrals(&solution);
        close(normalization, 1.0, 2e-12);
        close(moment, actual.mean_cosine, 2e-12);
        let table = solution.phase_table(cx, 0.002, 4096).unwrap();
        assert!(table.max_sampled_relative_error <= 0.002);
        close(table.phase.mean_cosine(), actual.mean_cosine, 0.002);
    });
}

#[test]
fn g0_zero_scattering_dilutes_mass_without_creating_an_angular_law() {
    with_cx(|cx, _| {
        let matched = sphere_for_x(2.0, 1.0, 0.0);
        let clear = MiePopulation::try_new(&[(matched, 1.0), (matched, 2.0)]).unwrap();
        let solution = clear
            .at_wavelength(cx, 550.0, MiePopulationBudget::default())
            .unwrap();
        assert_eq!(solution.cross_sections().extinction_m2(), 0.0);
        assert_eq!(solution.cross_sections().mean_cosine, 0.0);
        assert_eq!(solution.eval_cos(0.0), Err(MieError::NoScattering));
        assert_eq!(
            solution.phase_table(cx, 0.001, 4096),
            Err(MieError::NoScattering)
        );
        let visible = clear
            .visible_spectrum(
                cx,
                &[380.0, 780.0],
                MiePopulationBudget::default(),
                0.002,
                4096,
            )
            .unwrap();
        assert_eq!(visible.scattering.at(550.0).unwrap(), 0.0);
        assert_eq!(visible.absorption.at(550.0).unwrap(), 0.0);
        assert!(visible.angular_reports.is_empty());
        let active = sphere_for_x(1.0, 1.5, 0.0);
        let mixture = MiePopulation::try_new(&[(active, 1.0), (matched, 1.0)])
            .unwrap()
            .at_wavelength(cx, 550.0, MiePopulationBudget::default())
            .unwrap();
        let native = active
            .at_wavelength(cx, 550.0, MieBudget::default())
            .unwrap();
        close(
            mixture.cross_sections().scattering_m2 / native.cross_sections().scattering_m2,
            0.5,
            1e-15,
        );
        close(
            mixture.cross_sections().mass_scattering_m2_kg
                / native.cross_sections().mass_scattering_m2_kg,
            1.0 / 9.0,
            1e-15,
        );
        assert_eq!(mixture.eval_cos(0.3), native.eval_cos(0.3));
    });
}

fn lognormal_population(
    cx: &Cx<'_>,
    bins: usize,
) -> (MiePopulation, super::super::LognormalRadiusReport) {
    let (sizes, report) =
        ParticleSizeDistribution::lognormal(cx, 1e-9, 1.6, [1e-10, 1e-8], bins).unwrap();
    let particles: Vec<_> = sizes
        .bins()
        .iter()
        .map(|&(radius, weight)| {
            (
                MieSphere::try_new(radius, 1.5, 0.0, 1.0, 1800.0).unwrap(),
                weight,
            )
        })
        .collect();
    (MiePopulation::try_new(&particles).unwrap(), report)
}

#[test]
fn g1_lognormal_radius_refinement_matches_independent_rayleigh_moments() {
    with_cx(|cx, _| {
        let (coarse, report) = lognormal_population(cx, 16);
        let (fine, fine_report) = lognormal_population(cx, 128);
        let (finest, _) = lognormal_population(cx, 256);
        assert_eq!(
            report.omitted_number_fraction,
            fine_report.omitted_number_fraction
        );
        assert!(report.omitted_number_fraction < 1.0e-6);
        assert!(report.omitted_volume_fraction > 100.0 * report.omitted_number_fraction);
        close(
            report.retained_number_fraction + report.omitted_number_fraction,
            1.0,
            2e-15,
        );
        // Independent closed-form truncated lognormal moments, completing the
        // square in exp(q*y) times the normal density. This reference uses
        // no production population reduction or phase table.
        let sigma = 1.6_f64.ln();
        let moment = |q: f64| {
            let z_hi = (10.0_f64.ln() - q * sigma * sigma) / sigma;
            let z_lo = (0.1_f64.ln() - q * sigma * sigma) / sigma;
            let probability = 0.5
                * (fs_math::det::erf(z_hi / 2.0_f64.sqrt())
                    - fs_math::det::erf(z_lo / 2.0_f64.sqrt()));
            1e-9_f64.powf(q) * (0.5 * q * q * sigma * sigma).exp() * probability
                / report.retained_number_fraction
        };
        let contrast = (1.5_f64.powi(2) - 1.0) / (1.5_f64.powi(2) + 2.0);
        let expected =
            2.0 * (core::f64::consts::TAU / 550e-9).powi(4) * contrast.powi(2) * moment(6.0)
                / (1800.0 * moment(3.0));
        let coefficient = |population: &MiePopulation| {
            population
                .at_wavelength(cx, 550.0, MiePopulationBudget::default())
                .unwrap()
                .cross_sections()
                .mass_scattering_m2_kg
        };
        let (a, b, c) = (
            coefficient(&coarse),
            coefficient(&fine),
            coefficient(&finest),
        );
        assert!((b / expected - 1.0).abs() < (a / expected - 1.0).abs() / 10.0);
        assert!((c / expected - 1.0).abs() < 0.002);
        assert!((b / c - 1.0).abs() < 0.003);
        // Spectral preparation retains the full population at original knots.
        let spectrum = coarse
            .visible_spectrum(
                cx,
                &[380.0, 550.0, 780.0],
                MiePopulationBudget::default(),
                0.003,
                4096,
            )
            .unwrap();
        assert_eq!(spectrum.scattering.at(550.0).unwrap(), a);
        assert_eq!(spectrum.absorption.at(550.0).unwrap(), 0.0);
        assert_eq!(spectrum.angular_reports.len(), 3);
        assert!(spectrum.max_relative_spectral_midpoint_error > 0.0);
    });
}

#[test]
fn g0_population_admission_work_budgets_and_cancellation_are_explicit() {
    with_cx(|cx, gate| {
        let sphere = sphere_for_x(1.0, 1.5, 0.0);
        let other_host = MieSphere::try_new(sphere.radius_m(), 1.5, 0.0, 1.33, 1000.0).unwrap();
        for bins in [
            vec![],
            vec![(sphere, 0.0)],
            vec![(sphere, -1.0)],
            vec![(sphere, f64::NAN)],
            vec![(sphere, 1.0), (other_host, 1.0)],
        ] {
            assert_eq!(MiePopulation::try_new(&bins), Err(MieError::Parameters));
        }
        assert_eq!(
            MiePopulation::try_new(&vec![(sphere, 1.0); 257]),
            Err(MieError::Parameters)
        );
        assert_eq!(
            MiePopulation::try_new(&[(sphere, f64::MAX), (sphere, f64::from_bits(1))]),
            Err(MieError::Numerical)
        );
        assert!(ParticleSizeDistribution::try_new(&[(0.0, 1.0)]).is_err());
        let point = ParticleSizeDistribution::lognormal(cx, 1e-7, 1.0, [1e-8, 1e-6], 12).unwrap();
        assert_eq!(point.0.bins(), &[(1e-7, 1.0)]);
        let population = MiePopulation::try_new(&[(sphere, 1.0), (sphere, 1.0)]).unwrap();
        for budget in [
            MiePopulationBudget {
                max_bins: 1,
                ..MiePopulationBudget::default()
            },
            MiePopulationBudget {
                max_total_multipoles: 8,
                ..MiePopulationBudget::default()
            },
        ] {
            assert_eq!(
                population.at_wavelength(cx, 550.0, budget),
                Err(MieError::Budget)
            );
        }
        let prepared = population
            .at_wavelength(cx, 550.0, MiePopulationBudget::default())
            .unwrap();
        gate.request();
        assert_eq!(
            population.at_wavelength(cx, 550.0, MiePopulationBudget::default()),
            Err(MieError::Execution)
        );
        assert_eq!(
            prepared.phase_table(cx, 0.003, 4096),
            Err(MieError::Execution)
        );
        assert_eq!(
            ParticleSizeDistribution::lognormal(cx, 1e-7, 1.5, [1e-8, 1e-6], 16),
            Err(MieError::Execution)
        );
    });
}
