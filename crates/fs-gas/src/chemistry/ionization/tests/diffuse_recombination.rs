use super::collisional_radiative::fixture;
use super::{EV, KB, ME, NA};
mod evolution;
mod illumination;
use crate::chemistry::diffuse_recombination::*;
use fs_matdb::Provenance;
use fs_qty::{Amount, Area, Energy, Volume};
use fs_thermochem::equilibrium::composition::FrozenIonizedGasStateV1;
use fs_thermochem::kinetics::atomic_ionization::photoionization::{
    AtomicPhotoionizationV1, PhotoionizationCrossSectionPointV1,
};
use fs_thermochem::kinetics::atomic_ionization::radiative_recombination::AtomicRadiativeRecombinationV1;

const SIGMA: f64 = 1e-22;

fn parcel(
    x: f64,
    te: f64,
    scale: f64,
    volume: f64,
) -> (FrozenIonizedGasStateV1, AtomicRadiativeRecombinationV1) {
    let (original, source) = fixture(x, te, 0.0, 0.0, 0.0);
    let amounts: Vec<_> = original
        .amounts()
        .iter()
        .map(|x| Amount::new(x.value() * scale))
        .collect();
    let state = FrozenIonizedGasStateV1::new(
        original.thermodynamics(),
        &amounts,
        Volume::new(volume),
        original.heavy_temperature(),
        original.electron_temperature(),
        1e-12,
        || false,
    )
    .unwrap();
    (state, source)
}
fn controls(tolerance: f64) -> DiffuseRecombinationSolve {
    DiffuseRecombinationSolve {
        relative_tolerance: tolerance,
        ..DiffuseRecombinationSolve::default()
    }
}
fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-300),
        "{actual} != {expected}, relative {}",
        (actual / expected - 1.0).abs()
    );
}
fn analytic(state: &FrozenIonizedGasStateV1, source: &AtomicRadiativeRecombinationV1) -> [f64; 3] {
    // Independent complete constant-sigma Milne integrals. Source support in
    // these fixtures extends hundreds of thermal energies beyond the peak.
    let theta = KB * state.electron_temperature().value();
    let chi = source.source().channel().ionization_energy_j_mol() / NA;
    let factor =
        2.0 * (2.0 / std::f64::consts::PI).sqrt() / (ME.powf(1.5) * 299_792_458.0_f64.powi(2));
    let pairs = (NA * state.amounts()[1].value()) * (NA * state.amounts()[2].value())
        / state.volume().value();
    let rate = factor * SIGMA / theta.sqrt()
        * (chi * chi + 2.0 * chi * theta + 2.0 * theta * theta)
        * pairs;
    let kinetic = factor
        * SIGMA
        * theta.sqrt()
        * (chi * chi + 4.0 * chi * theta + 6.0 * theta * theta)
        * pairs;
    [rate, chi * rate + kinetic, kinetic]
}

#[test]
fn g1_uniform_isotropic_slab_recovers_independent_thin_and_thick_escape() {
    let (state, source) = parcel(0.1, 10_000.0, 1.0, 0.2);
    let exact = analytic(&state, &source);
    // Independent scipy.integrate.quad of int_0^1 phi(tau/mu) dmu;
    // equivalently (1/2-E_3(tau))/tau, evaluated directly in the thin case
    // to avoid cancellation in the special-function subtraction.
    for (tau, escaped_fraction) in [
        (1e-6, 0.999_992_630_852_386_5),
        (0.1, 0.837_085_420_917_212_5),
        (1.0, 0.390_308_032_802_239_9),
        (20.0, 0.024_999_999_995_495_445),
        (1e6, 5e-7),
        (1e20, 5e-21),
    ] {
        let area = Area::new(SIGMA * NA * state.amounts()[0].value() / tau);
        let result = evaluate_diffuse_recombination(
            &[DiffuseRecombinationCell {
                state: &state,
                source: &source,
            }],
            area,
            controls(1e-7),
            || false,
        )
        .unwrap();
        let b = result.balance();
        let cell = result.cells()[0].balance();
        close(b.emitted_photons_per_s, exact[0], 2e-11);
        close(b.emitted_power_w, exact[1], 2e-11);
        close(cell.captured_electron_power_w, exact[2], 2e-11);
        close(b.escaped_photons_per_s, exact[0] * escaped_fraction, 3e-7);
        close(
            b.absorbed_photons_per_s,
            exact[0] * (1.0 - escaped_fraction),
            3e-7,
        );
        close(b.escaped_power_w, exact[1] * escaped_fraction, 3e-7);
        close(
            cell.photoelectron_power_w,
            exact[2] * (1.0 - escaped_fraction),
            3e-7,
        );
        close(
            result.interfaces()[0].left_photons_per_s,
            result.interfaces()[1].right_photons_per_s,
            2e-14,
        );
        assert_eq!(result.interfaces()[0].right_photons_per_s, 0.0);
        assert_eq!(result.interfaces()[1].left_photons_per_s, 0.0);
        assert_eq!(cell.nonlocal_absorbed_photons_per_s, 0.0);
        assert!(cell.net_ionization_rate_mol_per_s < 0.0);
        assert!(cell.net_gas_power_w < 0.0);
        close(
            -cell.net_ionization_rate_mol_per_s * NA,
            b.escaped_photons_per_s,
            3e-9,
        );
        close(
            -cell.net_electron_power_w,
            exact[2] * escaped_fraction,
            3e-7,
        );
        close(-cell.net_gas_power_w, b.escaped_power_w, 3e-9);
        assert!(cell.net_photon_flux_residual_per_s.abs() < 3e-13 * b.escaped_photons_per_s);
        assert!(cell.net_power_flux_residual_w.abs() < 3e-13 * b.escaped_power_w);
        assert!(b.photon_residual_per_s.abs() < 3e-14 * exact[0]);
        assert!(b.energy_residual_w.abs() < 3e-14 * exact[1]);
        assert!(cell.gas_energy_residual_w.abs() < 3e-14 * exact[1]);
        assert!(b.maximum_normalized_spectral_estimate <= 0.25);
        assert!(b.maximum_normalized_angular_estimate <= 0.25);
        assert!(b.maximum_normalized_tail_bound <= 0.25);
        assert_eq!(result.cells()[0].state(), &state);
        assert_eq!(result.cells()[0].source(), &source);
    }
}

#[test]
fn g1_bidirectional_nonlocal_absorption_matches_independent_two_slab_integrals() {
    let (a, sa) = parcel(0.2, 10_000.0, 1.0, 0.2); // tau=0.8
    let (b, sb) = parcel(0.15, 30_000.0, 2.0, 0.4); // tau=1.7
    let area = Area::new(SIGMA * NA);
    let cells = [
        DiffuseRecombinationCell {
            state: &a,
            source: &sa,
        },
        DiffuseRecombinationCell {
            state: &b,
            source: &sb,
        },
    ];
    let result = evaluate_diffuse_recombination(&cells, area, controls(1e-7), || false).unwrap();
    let ea = analytic(&a, &sa);
    let eb = analytic(&b, &sb);
    // Independent scalar angular integrals, not a second production sweep:
    // P(i->j)=1/2 int_0^1 phi(tau_i/mu)*(1-exp(-tau_j/mu)) dmu.
    close(
        result.cells()[1].balance().nonlocal_absorbed_photons_per_s,
        ea[0] * 0.205_021_775_055_637_12,
        3e-7,
    );
    close(
        result.cells()[0].balance().nonlocal_absorbed_photons_per_s,
        eb[0] * 0.096_480_835_320_299_82,
        3e-7,
    );
    for c in result.cells() {
        let r = c.balance();
        close(
            r.chemical_storage_rate_w + r.net_electron_power_w,
            r.net_gas_power_w,
            1e-12,
        );
        assert!(r.photon_residual_per_s.abs() < 1e-13 * result.balance().emitted_photons_per_s);
    }
    let reversed =
        evaluate_diffuse_recombination(&[cells[1], cells[0]], area, controls(1e-7), || false)
            .unwrap();
    close(
        result.balance().escaped_power_w,
        reversed.balance().escaped_power_w,
        2e-13,
    );
    close(
        result.interfaces()[0].left_power_w,
        reversed.interfaces()[2].right_power_w,
        2e-13,
    );
    close(
        result.cells()[0].balance().absorbed_power_w,
        reversed.cells()[1].balance().absorbed_power_w,
        2e-13,
    );
}

fn shaped_source(
    original: &AtomicRadiativeRecombinationV1,
    points: &[(f64, f64)],
) -> AtomicRadiativeRecombinationV1 {
    let gap = original.source().channel().ionization_energy_j_mol() / NA;
    let provenance = Provenance {
        source:
            "manufactured frequency-dependent bound-free law for independent transfer verification"
                .into(),
        license: "analytic fixture".into(),
        artifact: None,
    };
    let photo = AtomicPhotoionizationV1::new(
        original.source().channel().clone(),
        provenance.clone(),
        points
            .iter()
            .map(|&(excess_ev, sigma)| PhotoionizationCrossSectionPointV1 {
                photon_energy: Energy::new(gap + excess_ev * EV),
                cross_section: Area::new(sigma),
            })
            .collect(),
    )
    .unwrap();
    AtomicRadiativeRecombinationV1::new(photo, Area::new(SIGMA), provenance).unwrap()
}

#[test]
fn g1_full_spectrum_heats_a_cold_absorber_without_losing_its_high_energy_opacity() {
    let (cold, original) = parcel(0.1, 1000.0, 1.0, 0.2);
    let source = shaped_source(
        &original,
        &[
            (0.0, SIGMA),
            (2.0, SIGMA),
            (4.0, 0.02 * SIGMA),
            (1000.0, 0.02 * SIGMA),
        ],
    );
    let (hot, hot_source) = parcel(0.2, 50_000.0, 1.0, 0.2);
    let result = evaluate_diffuse_recombination(
        &[
            DiffuseRecombinationCell {
                state: &cold,
                source: &source,
            },
            DiffuseRecombinationCell {
                state: &hot,
                source: &hot_source,
            },
        ],
        Area::new(SIGMA * NA),
        controls(1e-7),
        || false,
    )
    .unwrap();
    // Independent scipy QUADPACK integration in photon excess-energy and mu,
    // using the explicit Milne formula and piecewise-linear cross sections.
    // No production moment, panel or transfer evaluator supplies this oracle.
    // Cold opacity falls by 50x over 2..4 eV above threshold; replacing the hot
    // spectrum by its mean photon would miss the strong near-edge absorption.
    for (cell, expected) in result.cells().iter().zip([
        [
            1.515_053_788_311_732_7e27,
            3_322_044_620.014_419_6,
            21_182_663.458_216_477,
            9.879_604_942_551_846e26,
            2_201_932_568.286_31,
            49_453_735.372_160_04,
        ],
        [
            1.552_436_779_383_027e27,
            5_058_789_351.240_075,
            1_676_480_719.983_945_6,
            1.092_779_822_146_368_1e27,
            3_315_198_460.482_351_3,
            934_348_770.759_644_2,
        ],
    ]) {
        let b = cell.balance();
        for (actual, reference) in [
            b.emitted_photons_per_s,
            b.emitted_power_w,
            b.captured_electron_power_w,
            b.self_absorbed_photons_per_s + b.nonlocal_absorbed_photons_per_s,
            b.absorbed_power_w,
            b.photoelectron_power_w,
        ]
        .into_iter()
        .zip(expected)
        {
            close(actual, reference, 3e-7);
        }
    }
    assert!(result.cells()[0].balance().net_electron_power_w > 0.0);
    assert!(result.cells()[1].balance().net_electron_power_w < 0.0);
    let gap = source.source().channel().ionization_energy_j_mol() / NA;
    assert!(
        result.cells()[1].emission_cutoff_j() - gap
            > 20.0 * (result.cells()[0].emission_cutoff_j() - gap)
    );
}

#[test]
fn g3_subdividing_one_uniform_slab_preserves_its_field_and_total_gas_source() {
    let (whole, source) = parcel(0.2, 15_000.0, 1.0, 0.4);
    let (quarter, _) = parcel(0.2, 15_000.0, 0.25, 0.1);
    let area = Area::new(SIGMA * NA / 3.0);
    let coarse = evaluate_diffuse_recombination(
        &[DiffuseRecombinationCell {
            state: &whole,
            source: &source,
        }],
        area,
        controls(1e-7),
        || false,
    )
    .unwrap();
    let quarters = [DiffuseRecombinationCell {
        state: &quarter,
        source: &source,
    }; 4];
    let refined =
        evaluate_diffuse_recombination(&quarters, area, controls(1e-7), || false).unwrap();
    close(
        coarse.balance().escaped_photons_per_s,
        refined.balance().escaped_photons_per_s,
        1e-11,
    );
    close(
        coarse.balance().escaped_power_w,
        refined.balance().escaped_power_w,
        1e-11,
    );
    close(
        coarse.balance().absorbed_power_w,
        refined.balance().absorbed_power_w,
        1e-11,
    );
    let gas: f64 = refined
        .cells()
        .iter()
        .map(|c| c.balance().net_gas_power_w)
        .sum();
    close(coarse.cells()[0].balance().net_gas_power_w, gas, 1e-11);
    assert!(refined.cells()[1].balance().nonlocal_absorbed_photons_per_s > 0.0);
}

#[test]
fn g0_dark_empty_source_domain_work_and_cancellation_are_explicit() {
    let area = Area::new(SIGMA * NA);
    let control = controls(1e-7);
    let empty = evaluate_diffuse_recombination(&[], area, control, || false).unwrap();
    assert!(empty.cells().is_empty());
    assert_eq!(empty.interfaces(), &[DiffuseRecombinationFlux::default()]);
    assert_eq!(empty.balance().emitted_photons_per_s, 0.0);
    let (neutral, source) = parcel(0.0, 1000.0, 1.0, 0.2);
    let dark = evaluate_diffuse_recombination(
        &[DiffuseRecombinationCell {
            state: &neutral,
            source: &source,
        }],
        area,
        control,
        || false,
    )
    .unwrap();
    assert_eq!(
        dark.cells()[0].balance(),
        DiffuseRecombinationCellBalance::default()
    );
    let (active, active_source) = parcel(0.1, 10_000.0, 1.0, 0.2);
    let cells = [DiffuseRecombinationCell {
        state: &active,
        source: &active_source,
    }];
    assert!(matches!(
        evaluate_diffuse_recombination(&cells, area, control, || true),
        Err(DiffuseRecombinationError::Cancelled)
    ));
    let mut polls = 0;
    assert!(matches!(
        evaluate_diffuse_recombination(&cells, area, control, || {
            polls += 1;
            polls >= 500
        }),
        Err(DiffuseRecombinationError::Cancelled)
    ));
    let mut small = control;
    small.maximum_work = 1;
    assert!(matches!(
        evaluate_diffuse_recombination(&cells, area, small, || false),
        Err(DiffuseRecombinationError::WorkBudget)
    ));
    small = control;
    small.maximum_spectral_panels = 1;
    assert!(matches!(
        evaluate_diffuse_recombination(&cells, area, small, || false),
        Err(DiffuseRecombinationError::SpectralBudget)
    ));
    small = control;
    small.maximum_angular_subdivisions = 1;
    small.relative_tolerance = 1e-10;
    small.source_quadrature.relative_tolerance = 1e-12;
    assert!(matches!(
        evaluate_diffuse_recombination(&cells, Area::new(area.value() * 9.0), small, || false),
        Err(DiffuseRecombinationError::AngularBudget)
    ));
    let short = shaped_source(&source, &[(0.0, SIGMA), (0.5, SIGMA)]);
    assert!(matches!(
        evaluate_diffuse_recombination(
            &[
                DiffuseRecombinationCell {
                    state: &neutral,
                    source: &short
                },
                cells[0]
            ],
            area,
            control,
            || false
        ),
        Err(DiffuseRecombinationError::UnresolvedTail { .. })
    ));
    let (_, foreign) = fixture(0.1, 10_000.0, 0.0, 1e8, 0.0);
    assert!(matches!(
        evaluate_diffuse_recombination(
            &[DiffuseRecombinationCell {
                state: &active,
                source: &foreign
            }],
            area,
            control,
            || false
        ),
        Err(DiffuseRecombinationError::Invalid(_))
    ));
    assert_eq!(active, parcel(0.1, 10_000.0, 1.0, 0.2).0);
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(evaluate_diffuse_recombination(&cells, Area::new(bad), control, || false).is_err());
    }
}
