//! G0/G1 adjacent-stage checks on manufactured original sources, not a fitted
//! cross section or experimental qualification of a particular element.

use super::*;
use crate::electron_collision::{
    ElectronCollisionKindV1, ElectronCrossSectionPointV1, ElectronCrossSectionV1,
};
use crate::kinetics::atomic_ionization::photoionization::{
    AtomicPhotoionizationV1, PhotoionizationCrossSectionPointV1, PhotonNumberFluxV1,
};
use crate::kinetics::atomic_ionization::radiative_recombination::AtomicRadiativeRecombinationV1;
use fs_qty::{Area, Energy};

const GAP: f64 = 2.0 * CHI;
const UPPER_MASS: f64 = 4.0 * MP;
const SIGMA: f64 = 1e-22;

fn provenance() -> Provenance {
    Provenance {
        source: "manufactured adjacent atomic stage and explicit all-energy envelope".into(),
        license: "analytic test fixture".into(),
        artifact: None,
    }
}

fn charged_channel(
    lower_charge: i128,
    upper_charge: i128,
) -> Result<AtomicIonizationV1, AtomicIonizationErrorV1> {
    let models = vec![
        source(
            "X1",
            UPPER_MASS + ELECTRON_MASS_KG,
            3.0,
            CHI,
            2.5,
            200_000.0,
        ),
        source("X2", UPPER_MASS, 2.0, CHI + GAP, 2.5, 200_000.0),
        source("e", ELECTRON_MASS_KG, 2.0, 0.0, 2.5, 200_000.0),
    ];
    let ids = models
        .iter()
        .map(|m| m.species().clone())
        .collect::<Vec<_>>();
    let bank = IdealGasEquilibriumV1::new(
        ElementalMatrix::new(
            vec![ElementId::new("X").unwrap()],
            ids.clone(),
            vec![vec![1, 1, 0]],
        )
        .unwrap(),
        ChargeVector::new(ids, vec![lower_charge, upper_charge, -1]).unwrap(),
        models,
    )
    .unwrap();
    let threshold = GAP / AVOGADRO_PER_MOL;
    let cross_section = ElectronCrossSectionV1::tabulated(
        id("X1"),
        ElectronCollisionKindV1::Ionization,
        provenance(),
        [
            (threshold, 0.0),
            (2.0 * threshold, SIGMA),
            (1000.0 * threshold, SIGMA),
        ]
        .map(|(e, s)| ElectronCrossSectionPointV1 {
            energy: Energy::new(e),
            cross_section: Area::new(s),
        })
        .to_vec(),
        Area::new(SIGMA),
    )
    .unwrap();
    AtomicIonizationV1::new(
        bank,
        ReactionId::new("adjacent-stage-ionization").unwrap(),
        id("X2"),
        ElectronImpactRateV1::new(
            id("e"),
            cross_section,
            MaxwellianRateControlsV1 {
                relative_tolerance: 1e-10,
                maximum_subdivisions: 4096,
            },
            0.1,
        )
        .unwrap(),
        AtomicStatisticalWeightsV1 { neutral: 3, ion: 2 },
        1e-12,
    )
}

#[test]
fn adjacent_ionic_stages_match_independent_saha_impact_and_common_temperature_rates() {
    let te = 50_000.0;
    let theta = BOLTZMANN_J_PER_K * te;
    let gap = GAP / AVOGADRO_PER_MOL;
    let x = gap / theta;
    // Exact Maxwellian integral of the ramp-then-constant impact law, with
    // negligible bounded tail beyond the authored 1000-threshold coverage.
    let impact = (8.0 / (core::f64::consts::PI * ELECTRON_MASS_KG)).sqrt()
        * SIGMA
        * theta.sqrt()
        * ((-x).exp() * (1.0 + 2.0 / x) - (-2.0 * x).exp() * (2.0 + 2.0 / x));
    let kc = (4.0 / 3.0) / AVOGADRO_PER_MOL
        * (2.0 * core::f64::consts::PI * ELECTRON_MASS_KG * theta / PLANCK_J_S.powi(2)).powf(1.5)
        * (UPPER_MASS / (UPPER_MASS + ELECTRON_MASS_KG)).powf(1.5)
        * (-x).exp();
    for z in [0, 1, 7] {
        let channel = charged_channel(z, z + 1).unwrap();
        assert_eq!(channel.lower_charge(), z);
        assert_eq!(channel.lower_index(), channel.neutral_index());
        assert_eq!(channel.upper_index(), channel.ion_index());
        close(channel.ionization_energy_j_mol(), GAP, 2e-15, 0.0);
        close(
            channel
                .log_equilibrium_constant(Temperature::new(te))
                .unwrap()
                .exp(),
            kc,
            2e-13,
            0.0,
        );
        let c = [0.8, 0.2, 0.5];
        let value = evaluate(&channel, 4000.0, te, c);
        close(
            value.forward_rate_mol_m3_s,
            AVOGADRO_PER_MOL * impact * c[0] * c[2],
            2e-10,
            0.0,
        );
        close(
            value.reverse_rate_mol_m3_s,
            AVOGADRO_PER_MOL * impact / kc * c[1] * c[2].powi(2),
            2e-10,
            0.0,
        );
        assert!(value.entropy_production_w_m3_k > 0.0);
        let equilibrium = evaluate(&channel, 4000.0, te, [1.0, kc / 0.5, 0.5]);
        assert!(
            equilibrium.progress_rate_mol_m3_s.abs() <= 3e-13 * equilibrium.forward_rate_mol_m3_s
        );
        let mechanism = channel.common_temperature_mechanism();
        let mut work = KineticWorkspaceV1::new(mechanism).unwrap();
        mechanism
            .evaluate(Temperature::new(te), &c, &mut work, || false)
            .unwrap();
        close(
            work.reaction_rates()[0],
            value.progress_rate_mol_m3_s,
            2e-13,
            0.0,
        );
        close(
            work.temperature_derivatives()[0],
            value.electron_temperature_derivative_mol_m3_s_k,
            2e-13,
            0.0,
        );
        for (a, b) in work
            .rate_jacobian()
            .iter()
            .zip(value.concentration_derivatives_per_s)
        {
            close(*a, b, 2e-13, 0.0);
        }
        let rates = work.species_rates();
        assert_eq!(rates[0] + rates[1], 0.0);
        close(
            z as f64 * rates[0] + (z + 1) as f64 * rates[1] - rates[2],
            0.0,
            0.0,
            4e-15 * rates[1].abs(),
        );
    }
}

#[test]
fn ionic_photoionization_and_milne_capture_keep_original_gap_and_level_weights() {
    let channel = charged_channel(1, 2).unwrap();
    let gap = GAP / AVOGADRO_PER_MOL;
    let te = 30_000.0;
    let theta = BOLTZMANN_J_PER_K * te;
    assert!(AtomicPhotoionizationV1::verner_hydrogen(channel.clone()).is_err());
    let photo = AtomicPhotoionizationV1::new(
        channel,
        provenance(),
        [gap, gap + 100.0 * theta]
            .map(|e| PhotoionizationCrossSectionPointV1 {
                photon_energy: Energy::new(e),
                cross_section: Area::new(SIGMA),
            })
            .to_vec(),
    )
    .unwrap();
    let flux = 2e22;
    let p = photo
        .evaluate(
            &[PhotonNumberFluxV1 {
                photon_energy: Energy::new(1.5 * gap),
                number_flux_per_m2_s: flux,
            }],
            || false,
        )
        .unwrap();
    close(p.rate_per_neutral_per_s, SIGMA * flux, 2e-15, 0.0);
    close(
        p.electron_heating_per_neutral_w,
        0.5 * gap * SIGMA * flux,
        2e-15,
        0.0,
    );
    let rr = AtomicRadiativeRecombinationV1::new(photo, Area::new(SIGMA), provenance()).unwrap();
    let value = rr
        .evaluate(
            Temperature::new(te),
            MaxwellianRateControlsV1 {
                relative_tolerance: 1e-10,
                maximum_subdivisions: 4096,
            },
            || false,
        )
        .unwrap();
    // Independent complete polynomial Maxwellian moments. Electron spin is
    // already accounted for in Milne reciprocity; there is no extra factor 2.
    let factor = (3.0 / 2.0) * (2.0 / core::f64::consts::PI).sqrt()
        / (ELECTRON_MASS_KG.powf(1.5) * 299_792_458.0_f64.powi(2))
        * SIGMA;
    let alpha = factor / theta.sqrt() * (gap * gap + 2.0 * gap * theta + 2.0 * theta * theta);
    let kinetic = factor * theta.sqrt() * (gap * gap + 4.0 * gap * theta + 6.0 * theta * theta);
    close(value.coefficient_m3_s, alpha, 2e-10, 0.0);
    close(value.captured_kinetic_energy_j_m3_s, kinetic, 2e-10, 0.0);
    close(
        value.emitted_photon_energy_j_m3_s,
        gap * alpha + kinetic,
        2e-10,
        0.0,
    );
}

#[test]
fn negative_targets_and_nonadjacent_stages_refuse() {
    for (lower, upper) in [(-1, 0), (1, 3), (2, 1)] {
        assert!(
            charged_channel(lower, upper).is_err(),
            "unexpected admission: {lower}->{upper}"
        );
    }
}
