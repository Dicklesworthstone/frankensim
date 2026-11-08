use super::{channel, close};
use crate::electron_collision::{
    AVOGADRO_PER_MOL, ELECTRON_MASS_KG, ELECTRONVOLT_J, ElectronCollisionKindV1,
    ElectronCrossSectionPointV1, ElectronCrossSectionV1, MaxwellianRateControlsV1,
};
use crate::kinetics::ElectronImpactRateV1;
use crate::kinetics::atomic_ionization::photoionization::*;
use crate::kinetics::atomic_ionization::radiative_recombination::{
    AtomicRadiativeRecombinationV1, RadiativeRecombinationErrorV1,
};
use fs_matdb::Provenance;
use fs_qty::{Area, Energy, Temperature};

fn source() -> AtomicPhotoionizationV1 {
    let atomic = channel();
    let gap = atomic.ionization_energy_j_mol() / AVOGADRO_PER_MOL;
    AtomicPhotoionizationV1::new(
        atomic,
        Provenance {
            source: "synthetic linear event table for independent algebra tests".into(),
            license: "analytic test fixture".into(),
            artifact: None,
        },
        vec![
            PhotoionizationCrossSectionPointV1 {
                photon_energy: Energy::new(gap),
                cross_section: Area::new(2e-22),
            },
            PhotoionizationCrossSectionPointV1 {
                photon_energy: Energy::new(2.0 * gap),
                cross_section: Area::new(6e-22),
            },
        ],
    ).unwrap()
}

#[test]
fn photoedge_is_nonzero_below_edge_is_zero_and_upper_support_refuses() {
    let source = source();
    let gap = source.threshold().value();
    assert_eq!(source.channel(), &channel());
    assert_eq!(source.cross_section(Energy::new(0.5 * gap)).unwrap().value(), 0.0);
    assert_eq!(source.cross_section(source.threshold()).unwrap().value(), 2e-22);
    close(source.cross_section(Energy::new(1.5 * gap)).unwrap().value(),
        4e-22, 2e-15, 0.0);
    assert_eq!(source.cross_section(source.maximum_energy()).unwrap().value(), 6e-22);
    assert_eq!(source.cross_section(Energy::new(2.01 * gap)),
        Err(PhotoionizationErrorV1::OutsideEnergyDomain));
}

#[test]
fn discrete_number_flux_moments_close_original_chemical_and_electron_energy() {
    let source = source();
    let gap = source.threshold().value();
    let result = source.evaluate(&[
        PhotonNumberFluxV1 {
            photon_energy: Energy::new(gap), number_flux_per_m2_s: 3e20,
        },
        PhotonNumberFluxV1 {
            photon_energy: Energy::new(1.5 * gap), number_flux_per_m2_s: 2e20,
        },
    ], || false).unwrap();
    // Independent event algebra: 0.06 + 0.08 events/(atom*s).
    close(result.rate_per_neutral_per_s, 0.14, 3e-15, 0.0);
    close(result.absorbed_power_per_neutral_w, 0.18 * gap, 3e-15, 0.0);
    close(result.chemical_power_per_neutral_w, 0.14 * gap, 3e-15, 0.0);
    close(result.electron_heating_per_neutral_w, 0.04 * gap, 3e-15, 0.0);
    close(result.absorbed_power_per_neutral_w,
        result.chemical_power_per_neutral_w + result.electron_heating_per_neutral_w,
        3e-15, 0.0);
    assert_eq!(source.evaluate(&[], || false).unwrap().rate_per_neutral_per_s, 0.0);
}

#[test]
fn source_gap_order_and_local_field_admission_refuse() {
    let source = source();
    let mut points = source.points().unwrap().to_vec();
    points[0].photon_energy = Energy::new(0.9 * source.threshold().value());
    assert!(AtomicPhotoionizationV1::new(
        channel(), source.provenance().clone(), points,
    ).is_err());
    let mut points = source.points().unwrap().to_vec();
    points.swap(0, 1);
    assert!(AtomicPhotoionizationV1::new(
        channel(), source.provenance().clone(), points,
    ).is_err());
    for flux in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(source.evaluate(&[PhotonNumberFluxV1 {
            photon_energy: source.threshold(), number_flux_per_m2_s: flux,
        }], || false).is_err());
    }
    assert_eq!(source.evaluate(&[PhotonNumberFluxV1 {
        photon_energy: Energy::new(source.maximum_energy().value() * 2.0),
        number_flux_per_m2_s: 0.0,
    }], || false), Err(PhotoionizationErrorV1::OutsideEnergyDomain));
}

#[test]
fn cancellation_is_observed_inside_group_work_without_a_partial_result() {
    let source = source();
    let groups = [PhotonNumberFluxV1 {
        photon_energy: source.threshold(), number_flux_per_m2_s: 1e20,
    }; 65];
    let mut polls = 0;
    let result = source.evaluate(&groups, || { polls += 1; polls == 2 });
    assert_eq!(result, Err(PhotoionizationErrorV1::Cancelled));
    assert_eq!(polls, 2);
}

#[test]
fn verner_hydrogen_matches_independent_published_parameter_anchors() {
    let source = AtomicPhotoionizationV1::verner_hydrogen(channel()).unwrap();
    // Direct-power evaluation of the literal Table 1 H I row, independent of
    // the production logarithmic operation tree; these are numerical anchors,
    // not measured-data uncertainties or a bound on the fit's physical error.
    for (energy_ev, expected_m2) in [
        (13.6, 6.346_296_358_990_503e-22),
        (20.0, 2.211_102_984_039_216_3e-22),
        (100.0, 1.939_792_233_567_902_5e-24),
        (1000.0, 1.140_867_331_110_972e-27),
        (50_000.0, 1.715_120_863_209_725_6e-33),
    ] {
        close(
            source.cross_section(Energy::new(energy_ev * ELECTRONVOLT_J)).unwrap().value(),
            expected_m2, 1e-13, 0.0,
        );
    }
}

#[test]
fn verner_hydrogen_keeps_published_edge_and_upper_energy_domain() {
    let source = AtomicPhotoionizationV1::verner_hydrogen(channel()).unwrap();
    assert!(source.points().is_none());
    let edge = source.threshold().value();
    let maximum = source.maximum_energy().value();
    assert_eq!(edge, 13.60 * ELECTRONVOLT_J);
    assert_eq!(maximum, 50_000.0 * ELECTRONVOLT_J);
    assert_eq!(
        source.cross_section(Energy::new(f64::from_bits(edge.to_bits() - 1))).unwrap().value(),
        0.0,
    );
    assert!(source.cross_section(Energy::new(edge)).unwrap().value() > 0.0);
    assert!(source.cross_section(Energy::new(maximum)).unwrap().value() > 0.0);
    assert_eq!(
        source.cross_section(Energy::new(f64::from_bits(maximum.to_bits() + 1))),
        Err(PhotoionizationErrorV1::OutsideEnergyDomain),
    );
}

#[test]
fn verner_hydrogen_refuses_wrong_element_mass_gap_and_statistics() {
    use super::{
        AtomicIonizationV1, AtomicStatisticalWeightsV1, ChargeVector, ElementId,
        ElementalMatrix, IdealGasEquilibriumV1, ReactionId,
    };
    let models = super::models();
    let ids = models.iter().map(|m| m.species().clone()).collect::<Vec<_>>();
    // Species still named H/H+/e, but their declared nuclear row is helium.
    let wrong_element = IdealGasEquilibriumV1::new(
        ElementalMatrix::new(
            vec![ElementId::new("He").unwrap()], ids.clone(), vec![vec![1, 1, 0]],
        ).unwrap(),
        ChargeVector::new(ids, vec![0, 1, -1]).unwrap(),
        models,
    ).unwrap();
    let atomic = AtomicIonizationV1::new(
        wrong_element, ReactionId::new("wrong-element").unwrap(), super::id("H+"),
        super::impact(), super::WEIGHTS, 1e-12,
    ).unwrap();
    assert!(AtomicPhotoionizationV1::verner_hydrogen(atomic).is_err());

    let mut models = super::models();
    models[0] = super::source(
        "H", 2.0 * super::MP + super::ELECTRON_MASS_KG, 2.0, 0.0, 2.5, 200_000.0,
    );
    models[1] = super::source("H+", 2.0 * super::MP, 1.0, super::CHI, 2.5, 200_000.0);
    assert!(AtomicPhotoionizationV1::verner_hydrogen(super::bind(models).unwrap()).is_err());

    let mut models = super::models();
    models[1] = super::source("H+", super::MP, 1.0, 0.95 * super::CHI, 2.5, 200_000.0);
    assert!(AtomicPhotoionizationV1::verner_hydrogen(super::bind(models).unwrap()).is_err());

    let mut models = super::models();
    models[0] = super::source("H", super::MA, 4.0, 0.0, 2.5, 200_000.0);
    let atomic = AtomicIonizationV1::new(
        super::bank(models, 1), ReactionId::new("wrong-statistics").unwrap(), super::id("H+"),
        super::impact(), AtomicStatisticalWeightsV1 { neutral: 4, ion: 1 }, 1e-12,
    ).unwrap();
    assert!(AtomicPhotoionizationV1::verner_hydrogen(atomic).is_err());
}

#[test]
fn verner_heating_uses_original_gap_instead_of_rounded_fit_threshold() {
    let atomic = channel();
    let gap = atomic.ionization_energy_j_mol() / AVOGADRO_PER_MOL;
    let source = AtomicPhotoionizationV1::verner_hydrogen(atomic.clone()).unwrap();
    assert_eq!(source.channel(), &atomic);
    let energy = source.threshold().value();
    let result = source.evaluate(&[PhotonNumberFluxV1 {
        photon_energy: Energy::new(energy), number_flux_per_m2_s: 1e20,
    }], || false).unwrap();
    assert!(energy > gap);
    assert!(result.electron_heating_per_neutral_w > 0.0);
    close(result.chemical_power_per_neutral_w,
        result.rate_per_neutral_per_s * gap, 2e-15, 0.0);
    close(result.electron_heating_per_neutral_w,
        result.rate_per_neutral_per_s * (energy - gap), 2e-15, 0.0);
    close(result.absorbed_power_per_neutral_w,
        result.chemical_power_per_neutral_w + result.electron_heating_per_neutral_w,
        2e-15, 0.0);
}

// Independent ground-state statistical sources with arbitrary display names.
// NIST ASD He I/He II energies and CODATA-2022 alpha mass. The impact table is
// manufactured and is not a helium cross-section calibration.
const HE_GAPS_EV: [f64; 2] = [24.587_389_011, 54.417_765_528_2];
const HE_ALPHA_MASS: f64 = 6.644_657_345_0e-27;

fn helium_channel(
    stage: usize,
    element: &str,
    nucleus_mass: f64,
    weights: [u32; 2],
    gap_ev: f64,
) -> super::AtomicIonizationV1 {
    use super::{AtomicIonizationV1, AtomicStatisticalWeightsV1, ChargeVector, ElementId,
        ElementalMatrix, IdealGasEquilibriumV1, ReactionId, id};
    let gap = gap_ev * ELECTRONVOLT_J;
    let upper_mass = nucleus_mass + if stage == 0 { ELECTRON_MASS_KG } else { 0.0 };
    let models = vec![
        super::source("A", upper_mass + ELECTRON_MASS_KG, f64::from(weights[0]),
            0.0, 2.5, 200_000.0),
        super::source("B", upper_mass, f64::from(weights[1]), gap * AVOGADRO_PER_MOL,
            2.5, 200_000.0),
        super::source("e", ELECTRON_MASS_KG, 2.0, 0.0, 2.5, 200_000.0),
    ];
    let ids = models.iter().map(|m| m.species().clone()).collect::<Vec<_>>();
    let bank = IdealGasEquilibriumV1::new(
        ElementalMatrix::new(vec![ElementId::new(element).unwrap()], ids.clone(),
            vec![vec![1, 1, 0]]).unwrap(),
        ChargeVector::new(ids, vec![stage as i128, stage as i128 + 1, -1]).unwrap(),
        models,
    ).unwrap();
    let impact = ElectronCrossSectionV1::tabulated(
        id("A"), ElectronCollisionKindV1::Ionization,
        Provenance { source: "manufactured helium impact fixture".into(),
            license: "analytic fixture".into(), artifact: None },
        [(gap, 0.0), (10_000.0 * ELECTRONVOLT_J, 1e-22)]
            .map(|(energy, cross_section)| ElectronCrossSectionPointV1 {
                energy: Energy::new(energy), cross_section: Area::new(cross_section),
            }).to_vec(),
        Area::new(1e-22),
    ).unwrap();
    AtomicIonizationV1::new(
        bank, ReactionId::new("helium-photo-fixture").unwrap(), id("B"),
        ElectronImpactRateV1::new(id("e"), impact, MaxwellianRateControlsV1 {
            relative_tolerance: 1e-10, maximum_subdivisions: 4096,
        }, 0.1).unwrap(),
        AtomicStatisticalWeightsV1 { neutral: weights[0], ion: weights[1] }, 1e-12,
    ).unwrap()
}

fn helium_source(stage: usize) -> AtomicPhotoionizationV1 {
    let weights = if stage == 0 { [1, 2] } else { [2, 1] };
    let channel = helium_channel(stage, "He", HE_ALPHA_MASS, weights, HE_GAPS_EV[stage]);
    if stage == 0 {
        AtomicPhotoionizationV1::verner_helium(channel).unwrap()
    } else {
        AtomicPhotoionizationV1::verner_helium_ion(channel).unwrap()
    }
}

#[test]
fn g1_verner_helium_matches_independent_table_one_anchors() {
    // Independently evaluated literal Table 1 rows with 60-digit decimal
    // arithmetic; these check the fit equation, not its physical uncertainty.
    for (stage, anchors) in [
        (0, [(24.59, 7.434_698_694_110_645e-22), (30.0, 5.361_198_799_094_051e-22),
            (40.0, 3.158_999_711_315_489e-22), (54.42, 1.690_780_687_052_980_3e-22),
            (65.39, 1.125_306_418_211_808_6e-22)]),
        (1, [(54.42, 1.587_280_257_538_646_3e-22), (60.0, 1.219_642_943_356_784_3e-22),
            (100.0, 2.960_810_498_707_256e-23), (1000.0, 2.704_322_341_191_720_6e-26),
            (50_000.0, 5.237_174_244_234_447e-32)]),
    ] {
        let source = helium_source(stage);
        for (energy, expected) in anchors {
            close(source.cross_section(Energy::new(energy * ELECTRONVOLT_J)).unwrap().value(),
                expected, 1e-13, 0.0);
        }
    }
}

#[test]
fn g0_helium_support_and_proven_onset_envelope_reach_milne_admission() {
    for stage in [0, 1] {
        let source = helium_source(stage);
        let edge = source.threshold().value();
        let maximum = source.maximum_energy().value();
        let bound = source.cross_section(source.threshold()).unwrap();
        assert_eq!(edge, [24.59, 54.42][stage] * ELECTRONVOLT_J);
        assert_eq!(maximum, [65.39, 50_000.0][stage] * ELECTRONVOLT_J);
        assert!(source.points().is_none());
        assert_eq!(source.cross_section(Energy::new(f64::from_bits(edge.to_bits()-1)))
            .unwrap().value(), 0.0);
        assert_eq!(source.cross_section(Energy::new(f64::from_bits(maximum.to_bits()+1))),
            Err(PhotoionizationErrorV1::OutsideEnergyDomain));
        let mut previous = bound.value();
        for index in 1..=64 {
            let e = if index == 64 { maximum }
                else { edge * (maximum / edge).powf(f64::from(index) / 64.0) };
            let value = source.cross_section(Energy::new(e)).unwrap().value();
            assert!(value > 0.0 && value < previous);
            previous = value;
        }
        assert!(AtomicRadiativeRecombinationV1::new(
            source.clone(), bound, source.provenance().clone()).is_ok());
        assert!(AtomicRadiativeRecombinationV1::new(
            source.clone(), Area::new(0.99 * bound.value()), source.provenance().clone()).is_err());
    }
}

#[test]
fn g0_helium_fit_checks_nucleus_charge_mass_gap_and_ground_state_weights() {
    let constructors = [AtomicPhotoionizationV1::verner_helium,
        AtomicPhotoionizationV1::verner_helium_ion];
    for stage in [0, 1] {
        let weights = if stage == 0 { [1, 2] } else { [2, 1] };
        let gap = HE_GAPS_EV[stage];
        for (element, mass, weights, gap) in [
            ("H", HE_ALPHA_MASS, weights, gap),
            ("He", 0.75 * HE_ALPHA_MASS, weights, gap),
            ("He", HE_ALPHA_MASS, [3, 1], gap),
            ("He", HE_ALPHA_MASS, weights, gap - 0.01),
            ("He", HE_ALPHA_MASS, weights, gap + 0.01),
        ] {
            assert!(constructors[stage](helium_channel(stage, element, mass, weights, gap)).is_err());
        }
        assert!(constructors[1-stage](helium_channel(stage, "He", HE_ALPHA_MASS, weights, gap)).is_err());
        assert!(constructors[stage](channel()).is_err());
    }
}

#[test]
fn g1_helium_photon_heat_and_milne_capture_use_original_level_energy() {
    // Independent adaptive Gauss-Kronrod integration of literal Eq.(1) and
    // the Milne relation, retaining the rounded onsets and finite supports.
    // The He I 10 kK tail is resolved; it is not declared absent above 65.39 eV.
    for (stage, temperature, alpha, kinetic) in [
        (0, 10_000.0, 1.603_818_819_206_269_8e-19, 2.249_660_417_024_480_6e-38),
        (1, 50_000.0, 2.817_838_236_803_637e-19, 1.853_918_442_397_581_6e-37),
    ] {
        let source = helium_source(stage);
        let gap = source.channel().ionization_energy_j_mol() / AVOGADRO_PER_MOL;
        let photon_energy = source.threshold().value();
        let local = source.evaluate(&[PhotonNumberFluxV1 {
            photon_energy: Energy::new(photon_energy), number_flux_per_m2_s: 1e22,
        }], || false).unwrap();
        assert!(photon_energy > gap);
        close(gap, HE_GAPS_EV[stage] * ELECTRONVOLT_J, 2e-15, 0.0);
        close(local.electron_heating_per_neutral_w,
            local.rate_per_neutral_per_s * (photon_energy-gap), 2e-13, 0.0);
        close(local.absorbed_power_per_neutral_w,
            local.chemical_power_per_neutral_w + local.electron_heating_per_neutral_w,
            3e-15, 0.0);
        let source = AtomicRadiativeRecombinationV1::new(
            source.clone(), source.cross_section(source.threshold()).unwrap(),
            source.provenance().clone(),
        ).unwrap();
        let controls = MaxwellianRateControlsV1 { relative_tolerance: 1e-9,
            maximum_subdivisions: 4096 };
        let captured = source.evaluate(Temperature::new(temperature), controls, || false).unwrap();
        close(captured.coefficient_m3_s, alpha, 2e-10, 0.0);
        close(captured.captured_kinetic_energy_j_m3_s, kinetic, 2e-10, 0.0);
        close(captured.emitted_photon_energy_j_m3_s, gap * alpha + kinetic, 2e-10, 0.0);
        if stage == 0 {
            assert!(matches!(source.evaluate(Temperature::new(200_000.0), controls, || false),
                Err(RadiativeRecombinationErrorV1::UnresolvedTail { .. })));
        }
    }
}
