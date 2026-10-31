//! Independent ground-state statistical-mechanics oracles. These are exact
//! ideal-source fixtures, not measured plasma or total recombination validation.
//! The Saha constants below were evaluated separately with 70-digit decimal
//! arithmetic from h, kB, source masses and NIST H ionization energy:
//! <https://physics.nist.gov/cgi-bin/ASD/ie.pl?spectra=H+I>.
//! Maxwellian inverse-channel convention: NIST FLYCHK manual, pp.29–30,
//! <https://nlte.nist.gov/FLY/Doc/Manual_FLYCHK_Nov08.pdf>.

use super::*;
use crate::electron_collision::{MaxwellianRateControlsV1, nist};
use crate::kinetics::KineticWorkspaceV1;
use crate::{
    ChargeVector, ElementId, ElementalMatrix, ElementalReferenceIdV1, InitialStatePolicy,
    MolarEnergyQuantityV1, MolarThermalQuantityV1, NASA9_LAW_ID_V1, NASA9_LAW_VERSION_V1,
    NASA9_STATE_SCHEMA_VERSION_V1, Nasa9RegionV1, Nasa9StandardStateModelV1,
    ReferenceEquationOfStateV1, StandardStateConventionV1, StandardStatePhaseV1,
    required_nasa9_parameter_dims,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{ConstitutiveModelCard, LawId, LawParameter, Provenance};
use fs_qty::{MolarMass, Pressure};
use std::collections::BTreeMap;

const MP: f64 = 1.672_621_925_95e-27;
const MA: f64 = MP + ELECTRON_MASS_KG;
const CHI: f64 = 13.598_434_599_702 * 1.602_176_634e-19 * AVOGADRO_PER_MOL;
const P0: f64 = 100_000.0;
const WEIGHTS: AtomicStatisticalWeightsV1 = AtomicStatisticalWeightsV1 { neutral: 2, ion: 1 };

fn id(name: &str) -> SpeciesId {
    SpeciesId::new(name).unwrap()
}
fn close(actual: f64, expected: f64, relative: f64, absolute: f64) {
    assert!(
        actual.is_finite()
            && expected.is_finite()
            && (actual - expected).abs() <= absolute + relative * expected.abs(),
        "actual={actual:.17e}, expected={expected:.17e}, relative={relative:e}, absolute={absolute:e}"
    );
}

fn source(
    name: &str,
    mass: f64,
    degeneracy: f64,
    formation: f64,
    cp_r: f64,
    maximum_t: f64,
) -> Nasa9StandardStateModelV1 {
    // The statistical weight multiplies the translational partition function;
    // it is not inside its 3/2 power. Electron spin is counted once.
    let entropy_constant = 2.5
        + (degeneracy * BOLTZMANN_J_PER_K / P0
            * (2.0 * core::f64::consts::PI * mass * BOLTZMANN_J_PER_K / PLANCK_J_S.powi(2))
                .powf(1.5))
        .ln();
    let coefficients = [
        0.0,
        0.0,
        cp_r,
        0.0,
        0.0,
        0.0,
        0.0,
        formation / R,
        entropy_constant,
    ];
    let mut parameters = BTreeMap::new();
    for ((name, dims), value) in required_nasa9_parameter_dims()
        .into_iter()
        .take(9)
        .zip(coefficients)
    {
        parameters.insert(name.into(), LawParameter { value, dims });
    }
    parameters.insert(
        "reference_pressure".into(),
        LawParameter {
            value: P0,
            dims: Pressure::DIMS,
        },
    );
    let card = ConstitutiveModelCard {
        law: LawId(NASA9_LAW_ID_V1.into()),
        law_version: NASA9_LAW_VERSION_V1,
        parameters,
        state_schema_version: NASA9_STATE_SCHEMA_VERSION_V1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 1000.0, maximum_t),
        sources: vec![],
        provenance: Provenance {
            source: "Sackur-Tetrode ground-state analytic fixture; NIST H gap".into(),
            license: "analytic test fixture".into(),
            artifact: None,
        },
    };
    Nasa9StandardStateModelV1::new(
        id(name),
        MolarMass::new(mass * AVOGADRO_PER_MOL),
        StandardStateConventionV1::new(
            StandardStatePhaseV1::Gas,
            ReferenceEquationOfStateV1::IdealGas,
            Pressure::new(P0),
            ElementalReferenceIdV1::new("atomic-source-fixture").unwrap(),
        )
        .unwrap(),
        vec![Nasa9RegionV1::from_card(card).unwrap()],
    )
    .unwrap()
}
fn models() -> Vec<Nasa9StandardStateModelV1> {
    vec![
        source("H", MA, 2.0, 0.0, 2.5, 200_000.0),
        source("H+", MP, 1.0, CHI, 2.5, 200_000.0),
        source("e", ELECTRON_MASS_KG, 2.0, 0.0, 2.5, 200_000.0),
    ]
}
fn bank(models: Vec<Nasa9StandardStateModelV1>, nuclear_count: i128) -> IdealGasEquilibriumV1 {
    let ids = models
        .iter()
        .map(|model| model.species().clone())
        .collect::<Vec<_>>();
    IdealGasEquilibriumV1::new(
        ElementalMatrix::new(
            vec![ElementId::new("H").unwrap()],
            ids.clone(),
            vec![vec![nuclear_count, nuclear_count, 0]],
        )
        .unwrap(),
        ChargeVector::new(ids, vec![0, 1, -1]).unwrap(),
        models,
    )
    .unwrap()
}
fn impact() -> ElectronImpactRateV1 {
    ElectronImpactRateV1::new(
        id("e"),
        nist::hydrogen(id("H")),
        MaxwellianRateControlsV1 {
            relative_tolerance: 1e-10,
            maximum_subdivisions: 4096,
        },
        0.1,
    )
    .unwrap()
}
fn bind(
    models: Vec<Nasa9StandardStateModelV1>,
) -> Result<AtomicIonizationV1, AtomicIonizationErrorV1> {
    AtomicIonizationV1::new(
        bank(models, 1),
        ReactionId::new("ground-H-ionization").unwrap(),
        id("H+"),
        impact(),
        WEIGHTS,
        1e-12,
    )
}
fn channel() -> AtomicIonizationV1 {
    bind(models()).unwrap()
}
fn evaluate(
    channel: &AtomicIonizationV1,
    th: f64,
    te: f64,
    c: [f64; 3],
) -> AtomicIonizationEvaluationV1 {
    channel
        .evaluate(Temperature::new(th), Temperature::new(te), c, || false)
        .unwrap()
}

#[test]
fn source_factor_matches_independent_saha_and_original_common_temperature_owner() {
    let channel = channel();
    assert_eq!(channel.thermodynamics().models(), models());
    assert_eq!(channel.statistical_weights(), WEIGHTS);
    assert!(channel.statistical_entropy_discrepancy().abs() < 1e-12);
    close(channel.ionization_energy_j_mol(), CHI, 2e-15, 0.0);
    for (t, expected) in [
        (10_000.0, 0.000_561_623_323_886_333_2),
        (20_000.0, 4.242_724_496_057_948),
    ] {
        close(
            channel
                .log_equilibrium_constant(Temperature::new(t))
                .unwrap()
                .exp(),
            expected,
            3e-13,
            0.0,
        );
    }
    let mechanism = channel.common_temperature_mechanism();
    let mut work = KineticWorkspaceV1::new(mechanism).unwrap();
    let concentrations = [0.95, 0.05, 0.05];
    for t in [7000.0, 12_000.0, 20_000.0] {
        let atomic = evaluate(&channel, t, t, concentrations);
        mechanism
            .evaluate(Temperature::new(t), &concentrations, &mut work, || false)
            .unwrap();
        close(
            atomic.log_equilibrium_constant,
            work.log_equilibrium_constants()[0],
            3e-13,
            3e-13,
        );
        close(
            atomic.progress_rate_mol_m3_s,
            work.reaction_rates()[0],
            3e-12,
            1e-14,
        );
        close(
            atomic.electron_temperature_derivative_mol_m3_s_k,
            work.temperature_derivatives()[0],
            3e-12,
            1e-14,
        );
        for (actual, expected) in atomic
            .concentration_derivatives_per_s
            .into_iter()
            .zip(work.rate_jacobian())
        {
            close(actual, *expected, 3e-12, 1e-14);
        }
        assert_eq!(work.species_rates()[0], -work.species_rates()[1]);
        assert_eq!(work.species_rates()[1], work.species_rates()[2]);
        let mass_rate: f64 = work
            .species_rates()
            .iter()
            .zip(channel.thermodynamics().models())
            .map(|(rate, model)| rate * model.molar_mass().value())
            .sum();
        assert!(mass_rate.abs() <= 2e-15 * work.reaction_rates()[0].abs() * MA * AVOGADRO_PER_MOL);
    }
}

#[test]
fn separate_temperatures_obey_partial_balance_without_heavy_source_extrapolation() {
    let mut sources = models();
    sources[0] = source("H", MA, 2.0, 0.0, 2.5, 10_000.0);
    sources[1] = source("H+", MP, 1.0, CHI, 2.5, 10_000.0);
    let channel = bind(sources).unwrap();
    let te = 30_000.0;
    let kc = channel
        .log_equilibrium_constant(Temperature::new(te))
        .unwrap()
        .exp();
    for th in [1000.0, 4000.0, 10_000.0] {
        let evaluated = evaluate(&channel, th, te, [1.0, kc.sqrt(), kc.sqrt()]);
        assert!(evaluated.progress_rate_mol_m3_s.abs() <= 1e-12 * evaluated.forward_rate_mol_m3_s);
        assert!(evaluated.entropy_production_w_m3_k >= 0.0);
    }
    let cold = evaluate(&channel, 1000.0, te, [0.95, 0.05, 0.05]);
    let warm = evaluate(&channel, 10_000.0, te, [0.95, 0.05, 0.05]);
    assert_eq!(cold, warm);
    let mechanism = channel.common_temperature_mechanism();
    let mut work = KineticWorkspaceV1::new(mechanism).unwrap();
    assert!(
        mechanism
            .evaluate(Temperature::new(te), &[0.95, 0.05, 0.05], &mut work, || {
                false
            })
            .is_err()
    );
    assert!(
        channel
            .log_equilibrium_constant(Temperature::new(999.0))
            .is_err()
    );
}

#[test]
fn rate_derivatives_match_independent_centered_differences() {
    let channel = channel();
    let c = [0.95, 0.05, 0.05];
    let (th, te) = (3000.0, 12_000.0);
    let evaluated = evaluate(&channel, th, te, c);
    for i in 0..3 {
        let delta = c[i] * 1e-5;
        let mut minus = c;
        let mut plus = c;
        minus[i] -= delta;
        plus[i] += delta;
        let numerical = (evaluate(&channel, th, te, plus).progress_rate_mol_m3_s
            - evaluate(&channel, th, te, minus).progress_rate_mol_m3_s)
            / (2.0 * delta);
        close(
            evaluated.concentration_derivatives_per_s[i],
            numerical,
            3e-8,
            1e-12,
        );
    }
    let delta_t = 0.1;
    let numerical = (evaluate(&channel, th, te + delta_t, c).progress_rate_mol_m3_s
        - evaluate(&channel, th, te - delta_t, c).progress_rate_mol_m3_s)
        / (2.0 * delta_t);
    close(
        evaluated.electron_temperature_derivative_mol_m3_s_k,
        numerical,
        2e-8,
        1e-12,
    );
}

fn source_energy_entropy(
    channel: &AtomicIonizationV1,
    th: f64,
    te: f64,
    c: [f64; 3],
) -> (f64, f64) {
    let mut energy = 0.0;
    let mut entropy = 0.0;
    for (slot, index) in [
        channel.neutral_index(),
        channel.ion_index(),
        channel.electron_index(),
    ]
    .into_iter()
    .enumerate()
    {
        let t = if slot == 2 { te } else { th };
        let properties = channel.thermodynamics().models()[index]
            .evaluate_properties(Temperature::new(t))
            .unwrap();
        energy += c[slot] * properties.u().value();
        if c[slot] > 0.0 {
            entropy += c[slot] * (properties.s().value() - R * (c[slot] * R * t / P0).ln());
        }
    }
    (energy, entropy)
}

#[test]
fn independent_source_entropy_and_energy_include_population_dilution_and_formation_gauge() {
    // Shift the electron formation reference and opposite ionic charge reference.
    // The physical gap, rates and total entropy must remain unchanged, while the
    // full heavy-pool formation gain is correctly distinguished from chi.
    let electron_reference = 20_000.0;
    let mut sources = models();
    sources[1] = source("H+", MP, 1.0, CHI - electron_reference, 2.5, 200_000.0);
    sources[2] = source(
        "e",
        ELECTRON_MASS_KG,
        2.0,
        electron_reference,
        2.5,
        200_000.0,
    );
    let shifted = bind(sources).unwrap();
    close(shifted.ionization_energy_j_mol(), CHI, 2e-15, 0.0);
    close(
        shifted.heavy_energy_gap_j_mol(),
        CHI - electron_reference,
        2e-15,
        0.0,
    );
    close(
        shifted.formation_energies_j_mol()[2],
        electron_reference,
        2e-15,
        0.0,
    );
    let (c, th, te) = ([0.95, 0.05, 0.05], 4000.0, 30_000.0);
    let evaluated = evaluate(&shifted, th, te, c);
    close(
        evaluated.progress_rate_mol_m3_s,
        evaluate(&channel(), th, te, c).progress_rate_mol_m3_s,
        1e-13,
        0.0,
    );
    let initial_energy = source_energy_entropy(&shifted, th, te, c).0;
    let delta = 1e-7;
    let trial = |extent: f64| {
        let changed = [c[0] - extent, c[1] + extent, c[2] + extent];
        let next_te = (1.5 * R * c[2] * te - CHI * extent) / (1.5 * R * changed[2]);
        source_energy_entropy(&shifted, th, next_te, changed)
    };
    let (energy_plus, entropy_plus) = trial(delta);
    let (energy_minus, entropy_minus) = trial(-delta);
    close(energy_plus, initial_energy, 3e-15, 1e-9);
    close(energy_minus, initial_energy, 3e-15, 1e-9);
    let entropy_per_extent = (entropy_plus - entropy_minus) / (2.0 * delta);
    close(
        evaluated.entropy_production_w_m3_k,
        entropy_per_extent * evaluated.progress_rate_mol_m3_s,
        2e-7,
        1e-12,
    );
    assert!(evaluated.entropy_production_w_m3_k > 0.0);
}

#[test]
fn exact_absence_has_no_artificial_seed_and_boundary_entropy_is_explicit() {
    let channel = channel();
    let seedless = evaluate(&channel, 4000.0, 20_000.0, [1.0, 0.0, 0.0]);
    assert_eq!(seedless.progress_rate_mol_m3_s, 0.0);
    assert_eq!(seedless.forward_rate_mol_m3_s, 0.0);
    assert_eq!(seedless.reverse_rate_mol_m3_s, 0.0);
    assert_eq!(seedless.electron_temperature_derivative_mol_m3_s_k, 0.0);
    assert_eq!(seedless.entropy_production_w_m3_k, 0.0);
    assert!(seedless.concentration_derivatives_per_s[2] > 0.0);
    let reverse = evaluate(&channel, 4000.0, 20_000.0, [0.0, 1.0, 1.0]);
    assert!(reverse.progress_rate_mol_m3_s < 0.0);
    assert_eq!(reverse.dimensionless_affinity, Some(f64::NEG_INFINITY));
    assert_eq!(reverse.entropy_production_w_m3_k, f64::INFINITY);
    let absent = evaluate(&channel, 4000.0, 20_000.0, [0.0, 0.0, 1.0]);
    assert_eq!(absent.progress_rate_mol_m3_s, 0.0);
    assert_eq!(absent.dimensionless_affinity, None);
    assert_eq!(absent.entropy_production_w_m3_k, 0.0);
}

#[test]
fn invalid_states_and_cancellation_refuse_without_result() {
    let channel = channel();
    for c in [
        [-1.0, 1.0, 1.0],
        [1.0, f64::NAN, 1.0],
        [1.0, 1.0, f64::INFINITY],
    ] {
        assert!(
            channel
                .evaluate(
                    Temperature::new(4000.0),
                    Temperature::new(20_000.0),
                    c,
                    || false
                )
                .is_err()
        );
    }
    for (th, te) in [
        (0.0, 20_000.0),
        (4000.0, f64::NAN),
        (4000.0, 200_001.0),
        (200_000.0, 1000.0),
    ] {
        assert!(
            channel
                .evaluate(
                    Temperature::new(th),
                    Temperature::new(te),
                    [0.9, 0.1, 0.1],
                    || false
                )
                .is_err()
        );
    }
    assert_eq!(
        channel.evaluate(
            Temperature::new(4000.0),
            Temperature::new(20_000.0),
            [0.9, 0.1, 0.1],
            || true
        ),
        Err(AtomicIonizationErrorV1::Cancelled)
    );
    let mut polls = 0;
    let cancelled = channel.evaluate(
        Temperature::new(4000.0),
        Temperature::new(20_000.0),
        [0.9, 0.1, 0.1],
        || {
            polls += 1;
            polls >= 3
        },
    );
    assert_eq!(cancelled, Err(AtomicIonizationErrorV1::Cancelled));
    assert_eq!(polls, 3);
}

#[test]
fn unsupported_internal_modes_statistics_and_below_gap_cross_sections_refuse() {
    let mut sources = models();
    sources[0] = source("H", MA, 2.0, 0.0, 3.5, 200_000.0);
    assert!(matches!(
        bind(sources),
        Err(AtomicIonizationErrorV1::Invalid(_))
    ));
    let mut sources = models();
    sources[0] = source("H", MA, 3.0, 0.0, 2.5, 200_000.0);
    assert!(matches!(
        bind(sources),
        Err(AtomicIonizationErrorV1::StatisticalEntropyMismatch { .. })
    ));
    let mut sources = models();
    sources[1] = source(
        "H+",
        MP,
        1.0,
        14.0 * 1.602_176_634e-19 * AVOGADRO_PER_MOL,
        2.5,
        200_000.0,
    );
    assert!(matches!(
        bind(sources),
        Err(AtomicIonizationErrorV1::Invalid(_))
    ));
    let mut sources = models();
    sources[2] = source("e", 2.0 * ELECTRON_MASS_KG, 2.0, 0.0, 2.5, 200_000.0);
    assert!(matches!(
        bind(sources),
        Err(AtomicIonizationErrorV1::Kinetics(_))
    ));
    assert!(matches!(
        AtomicIonizationV1::new(
            bank(models(), 2),
            ReactionId::new("not-atomic").unwrap(),
            id("H+"),
            impact(),
            WEIGHTS,
            1e-12
        ),
        Err(AtomicIonizationErrorV1::Invalid(_))
    ));
    assert!(matches!(
        AtomicIonizationV1::new(
            bank(models(), 1),
            ReactionId::new("no-states").unwrap(),
            id("H+"),
            impact(),
            AtomicStatisticalWeightsV1 { neutral: 0, ion: 1 },
            1e-12
        ),
        Err(AtomicIonizationErrorV1::Invalid(_))
    ));
}

#[test]
fn transformed_caloric_references_cannot_silently_reinterpret_raw_level_constants() {
    let mut sources = models();
    sources[0] = sources[0]
        .with_caloric_continuity(
            Temperature::new(4000.0),
            MolarEnergyQuantityV1::new(0.0),
            MolarThermalQuantityV1::new(0.0),
        )
        .unwrap();
    assert!(sources[0].caloric_reference().is_some());
    assert!(matches!(
        bind(sources),
        Err(AtomicIonizationErrorV1::Invalid(_))
    ));
}
