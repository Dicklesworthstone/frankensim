//! Independent SI moment oracles and conservation/domain checks. The numerical
//! values were reconstructed with 70-digit Decimal arithmetic from the physical
//! constants and Maxwellian Landau moments, not from production evaluation.

use super::*;
use crate::equilibrium::composition::IdealGasEquilibriumV1;
use crate::{
    ChargeVector, ElementId, ElementalMatrix, ElementalReferenceIdV1, Nasa9RegionV1,
    Nasa9StandardStateModelV1, ReferenceEquationOfStateV1, StandardStateConventionV1,
    StandardStatePhaseV1, required_nasa9_parameter_dims,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_qty::{Amount, MolarMass, Pressure, Volume};
use std::collections::BTreeMap;

const MP: f64 = 1.672_621_925_95e-27;
fn id(name: &str) -> SpeciesId {
    SpeciesId::new(name).unwrap()
}
fn source(name: &str, mass: f64) -> Nasa9StandardStateModelV1 {
    let mut parameters = BTreeMap::new();
    for ((name, dims), value) in required_nasa9_parameter_dims()
        .into_iter()
        .take(9)
        .zip([0.0, 0.0, 2.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    {
        parameters.insert(name.into(), LawParameter { value, dims });
    }
    parameters.insert(
        "reference_pressure".into(),
        LawParameter {
            value: 100_000.0,
            dims: Pressure::DIMS,
        },
    );
    Nasa9StandardStateModelV1::new(
        id(name),
        MolarMass::new(mass * AVOGADRO_PER_MOL),
        StandardStateConventionV1::new(
            StandardStatePhaseV1::Gas,
            ReferenceEquationOfStateV1::IdealGas,
            Pressure::new(100_000.0),
            ElementalReferenceIdV1::new("classical-Coulomb-fixture").unwrap(),
        )
        .unwrap(),
        vec![
            Nasa9RegionV1::from_card(ConstitutiveModelCard {
                law: LawId("nasa9-standard-state".into()),
                law_version: 1,
                parameters,
                state_schema_version: 0,
                initial_state: InitialStatePolicy::ZeroInternalState,
                validity: ValidityDomain::unconstrained().with("T", 100.0, 1e10),
                sources: vec![],
                provenance: Provenance {
                    source: "constant-Cv ideal source for independent Coulomb moment tests".into(),
                    license: "analytic test fixture".into(),
                    artifact: None,
                },
            })
            .unwrap(),
        ],
    )
    .unwrap()
}
fn limits() -> CoulombLimitsV1 {
    CoulombLimitsV1 {
        minimum_coulomb_logarithm: 3.0,
        minimum_particles_per_debye_sphere: 10.0,
        maximum_coulomb_coupling: 0.1,
        maximum_thermal_phase_space_density: 0.1,
        maximum_thermal_energy_rest_energy_ratio: 0.01,
        maximum_quantum_to_classical_ratio: 0.5,
        maximum_relative_drift_thermal_ratio: 0.1,
    }
}
fn state(density: f64, th: f64, te: f64) -> FrozenIonizedGasStateV1 {
    let ids = vec![id("H"), id("H+"), id("e")];
    let bank = IdealGasEquilibriumV1::new(
        ElementalMatrix::new(
            vec![ElementId::new("H").unwrap()],
            ids.clone(),
            vec![vec![1, 1, 0]],
        )
        .unwrap(),
        ChargeVector::new(ids, vec![0, 1, -1]).unwrap(),
        vec![
            source("H", MP + ELECTRON_MASS_KG),
            source("H+", MP),
            source("e", ELECTRON_MASS_KG),
        ],
    )
    .unwrap();
    FrozenIonizedGasStateV1::new(
        &bank,
        &[
            Amount::new(1e21 / AVOGADRO_PER_MOL),
            Amount::new(density / AVOGADRO_PER_MOL),
            Amount::new(density / AVOGADRO_PER_MOL),
        ],
        Volume::new(1.0),
        Temperature::new(th),
        Temperature::new(te),
        1e-12,
        || false,
    )
    .unwrap()
}
fn bind(state: &FrozenIonizedGasStateV1) -> CoulombModelV1 {
    CoulombModelV1::new(
        state,
        CoulombScreeningV1::StaticAllChargedDebye,
        limits(),
        || false,
    )
    .unwrap()
}
fn pair(model: &CoulombModelV1) -> (CoulombSummaryV1, CoulombPairV1) {
    let mut pairs = Vec::new();
    let summary = model
        .evaluate_with_pairs(
            model.initial_state().heavy_temperature(),
            model.initial_state().electron_temperature(),
            |p| pairs.push(p),
            || false,
        )
        .unwrap();
    assert_eq!(pairs.len(), 1);
    (summary, pairs[0])
}
fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        actual.is_finite()
            && expected.is_finite()
            && (actual - expected).abs() <= relative * expected.abs(),
        "actual={actual:.17e}, expected={expected:.17e}, relative={relative:e}"
    );
}

#[test]
fn coulomb_moments_match_independent_finite_mass_si_oracles_for_both_screening_models() {
    let state = state(1e20, 10_000.0, 20_000.0);
    for (screening, debye, logarithm, k, g) in [
        (
            CoulombScreeningV1::StaticAllChargedDebye,
            5.634_559_679_554_459e-7,
            6.513_536_257_806_465,
            0.076_226_076_032_512_2,
            1886.574_047_487_622,
        ),
        (
            CoulombScreeningV1::ElectronDebye,
            9.759_343_643_267_335e-7,
            7.062_842_402_140_52,
            0.082_654_450_768_731_85,
            2045.674_523_021_799_7,
        ),
    ] {
        let model = CoulombModelV1::new(&state, screening, limits(), || false).unwrap();
        let (summary, pair) = pair(&model);
        close(summary.debye_length_m().unwrap(), debye, 4e-14);
        close(pair.coulomb_logarithm(), logarithm, 4e-14);
        close(pair.friction_coefficient_kg_m3_s(), k, 1e-13);
        close(pair.thermal_conductance_w_m3_k(), g, 1e-13);
        close(
            summary.electron_heavy_thermal_conductance_w_m3_k(),
            g,
            1e-13,
        );
        close(
            pair.relative_thermal_variance_m2_s2(),
            303_209_307_203.899_4,
            3e-14,
        );
        close(
            summary.maximum_quantum_to_classical_ratio(),
            0.125_850_584_289_518_15,
            4e-14,
        );
        assert_eq!((pair.first_index(), pair.second_index()), (1, 2));
        assert_eq!(summary.pair_count(), 1);
        // Both directions use ONE reciprocal K, not independently fitted nus.
        let nue = k / (1e20 * ELECTRON_MASS_KG);
        let nui = k / (1e20 * MP);
        close(nue / nui, MP / ELECTRON_MASS_KG, 3e-15);
        close(
            g / (1.5 * BOLTZMANN_J_PER_K * 1e20),
            2.0 * ELECTRON_MASS_KG * nue / (MP + ELECTRON_MASS_KG),
            1e-13,
        );
    }
}

#[test]
fn coulomb_equal_temperature_drift_and_unequal_temperature_signed_correction_conserve_energy() {
    let equal = bind(&state(1e20, 20_000.0, 20_000.0));
    let (_, same) = pair(&equal);
    let w2 = 0.25 * same.maximum_relative_speed_squared_m2_s2();
    let heat = same.drift_heating(w2).unwrap();
    assert_eq!(heat.thermal_exchange_correction_to_first_w_m3(), 0.0);
    close(
        heat.first_power_density_w_m3() / heat.second_power_density_w_m3(),
        ELECTRON_MASS_KG / MP,
        3e-15,
    );
    close(
        heat.first_power_density_w_m3() + heat.second_power_density_w_m3(),
        heat.total_power_density_w_m3(),
        3e-15,
    );
    // At Te/Th=50, drift suppresses some electron-to-ion thermal transfer.
    // The ion drift correction is signed, despite positive total dissipation.
    let hot = bind(&state(1e18, 1000.0, 50_000.0));
    let (_, different) = pair(&hot);
    let w2 = 0.25 * different.maximum_relative_speed_squared_m2_s2();
    let heat = different.drift_heating(w2).unwrap();
    assert!(heat.first_friction_power_density_w_m3() > 0.0);
    assert!(heat.second_friction_power_density_w_m3() > 0.0);
    assert!(heat.thermal_exchange_correction_to_first_w_m3() < 0.0);
    assert!(heat.first_power_density_w_m3() < 0.0);
    assert!(heat.second_power_density_w_m3() > heat.total_power_density_w_m3());
    close(
        heat.first_power_density_w_m3() + heat.second_power_density_w_m3(),
        heat.total_power_density_w_m3(),
        3e-15,
    );
    let thermal_to_ion = different.thermal_conductance_w_m3_k() * 49_000.0;
    close(
        heat.thermal_exchange_correction_to_first_w_m3(),
        -thermal_to_ion * w2 / (2.0 * different.relative_thermal_variance_m2_s2()),
        1e-13,
    );
    assert_eq!(
        different
            .drift_heating(0.0)
            .unwrap()
            .total_power_density_w_m3(),
        0.0
    );
    assert!(matches!(
        different.drift_heating(1.01 * different.maximum_relative_speed_squared_m2_s2()),
        Err(CoulombErrorV1::Domain(_))
    ));
    assert!(matches!(
        different.drift_heating(f64::NAN),
        Err(CoulombErrorV1::Invalid(_))
    ));
}

#[test]
fn coulomb_all_charged_pairs_include_ion_ion_drag_but_only_electron_heavy_thermal_sum() {
    let ids = vec![id("H"), id("H+"), id("He++"), id("e")];
    let bank = IdealGasEquilibriumV1::new(
        ElementalMatrix::new(
            vec![ElementId::new("H").unwrap(), ElementId::new("He").unwrap()],
            ids.clone(),
            vec![vec![1, 1, 0, 0], vec![0, 0, 1, 0]],
        )
        .unwrap(),
        ChargeVector::new(ids, vec![0, 1, 2, -1]).unwrap(),
        vec![
            source("H", MP + ELECTRON_MASS_KG),
            source("H+", MP),
            source("He++", 4.0 * MP),
            source("e", ELECTRON_MASS_KG),
        ],
    )
    .unwrap();
    let amounts = [1e21, 1e18, 2e18, 5e18].map(|n| Amount::new(n / AVOGADRO_PER_MOL));
    let state = FrozenIonizedGasStateV1::new(
        &bank,
        &amounts,
        Volume::new(1.0),
        Temperature::new(10_000.0),
        Temperature::new(20_000.0),
        1e-12,
        || false,
    )
    .unwrap();
    let model = bind(&state);
    let mut pairs = Vec::new();
    let summary = model
        .evaluate_with_pairs(
            state.heavy_temperature(),
            state.electron_temperature(),
            |p| pairs.push(p),
            || false,
        )
        .unwrap();
    assert_eq!(
        pairs
            .iter()
            .map(|p| (p.first_index(), p.second_index()))
            .collect::<Vec<_>>(),
        vec![(1, 2), (1, 3), (2, 3)]
    );
    close(
        summary.electron_heavy_thermal_conductance_w_m3_k(),
        pairs[1].thermal_conductance_w_m3_k() + pairs[2].thermal_conductance_w_m3_k(),
        3e-15,
    );
    assert!(pairs[0].friction_coefficient_kg_m3_s() > 0.0);
    assert_eq!(
        pairs[0]
            .drift_heating(1.0)
            .unwrap()
            .thermal_exchange_correction_to_first_w_m3(),
        0.0
    );
}

#[test]
fn coulomb_population_binding_allows_new_temperatures_but_refuses_inventory_reuse() {
    let initial = state(1e20, 10_000.0, 20_000.0);
    let model = bind(&initial);
    let changed_t = FrozenIonizedGasStateV1::new(
        initial.thermodynamics(),
        initial.amounts(),
        initial.volume(),
        Temperature::new(12_000.0),
        Temperature::new(24_000.0),
        1e-12,
        || false,
    )
    .unwrap();
    assert!(model.is_bound_to_populations(&changed_t));
    let refreshed = bind(&changed_t);
    assert_eq!(
        model
            .evaluate(
                changed_t.heavy_temperature(),
                changed_t.electron_temperature(),
                || false
            )
            .unwrap(),
        refreshed
            .evaluate(
                changed_t.heavy_temperature(),
                changed_t.electron_temperature(),
                || false
            )
            .unwrap()
    );
    assert!(!model.is_bound_to_populations(&state(2e20, 10_000.0, 20_000.0)));
    let changed_v = FrozenIonizedGasStateV1::new(
        initial.thermodynamics(),
        initial.amounts(),
        Volume::new(2.0),
        initial.heavy_temperature(),
        initial.electron_temperature(),
        1e-12,
        || false,
    )
    .unwrap();
    assert!(!model.is_bound_to_populations(&changed_v));
    assert!(matches!(
        model.evaluate(
            Temperature::new(99.0),
            initial.electron_temperature(),
            || false
        ),
        Err(CoulombErrorV1::Species { .. })
    ));
}

#[test]
fn coulomb_fresh_state_refreshes_absent_carriers_density_volume_and_temperature() {
    let neutral = state(0.0, 10_000.0, 20_000.0);
    for screening in [
        CoulombScreeningV1::StaticAllChargedDebye,
        CoulombScreeningV1::ElectronDebye,
    ] {
        let original = CoulombModelV1::new(&neutral, screening, limits(), || false).unwrap();
        let retained = original.clone();
        for (density, volume, th, te) in [
            (0.0, 2.0, 12_000.0, 24_000.0),
            (5e18, 0.25, 12_000.0, 18_000.0),
            (1e20, 2.0, 12_000.0, 24_000.0),
        ] {
            let base = state(density, th, te);
            let actual = FrozenIonizedGasStateV1::new(
                base.thermodynamics(),
                base.amounts(),
                Volume::new(volume),
                base.heavy_temperature(),
                base.electron_temperature(),
                1e-12,
                || false,
            )
            .unwrap();
            assert!(!original.is_bound_to_populations(&actual));
            let refreshed = original.evaluate_at_state(&actual, || false).unwrap();
            let rebound = CoulombModelV1::new(&actual, screening, limits(), || false).unwrap();
            let expected = rebound
                .evaluate(
                    actual.heavy_temperature(),
                    actual.electron_temperature(),
                    || false,
                )
                .unwrap();
            assert_eq!(refreshed, expected);
            assert_eq!(refreshed.pair_count(), usize::from(density > 0.0));
            assert_eq!(
                refreshed.electron_heavy_thermal_conductance_w_m3_k() > 0.0,
                density > 0.0
            );
        }
        // A successful refresh never relaxes the original material binding or
        // changes the coefficients subsequently evaluated at its old inventory.
        assert_eq!(original, retained);
        assert!(original.is_bound_to_populations(&neutral));
        assert_eq!(
            original
                .evaluate(
                    neutral.heavy_temperature(),
                    neutral.electron_temperature(),
                    || false
                )
                .unwrap()
                .pair_count(),
            0
        );
    }
    let charged = state(1e20, 10_000.0, 20_000.0);
    let model = bind(&charged);
    let vanished = model.evaluate_at_state(&neutral, || false).unwrap();
    assert_eq!(vanished.pair_count(), 0);
    assert_eq!(vanished.debye_length_m(), None);
}

#[test]
fn coulomb_fresh_state_refuses_foreign_bank_invalid_new_domain_and_cancellation() {
    let neutral = state(0.0, 10_000.0, 20_000.0);
    let model = bind(&neutral);
    let charged = state(1e20, 10_000.0, 20_000.0);
    let bank = charged.thermodynamics();
    let mut models = bank.models().to_vec();
    // Even a changed neutral source invalidates the complete original bank.
    models[0] = source("H", 2.0 * (MP + ELECTRON_MASS_KG));
    let foreign_bank = IdealGasEquilibriumV1::new(
        bank.elemental_matrix().clone(),
        bank.charge_vector().clone(),
        models,
    )
    .unwrap();
    let foreign = FrozenIonizedGasStateV1::new(
        &foreign_bank,
        charged.amounts(),
        charged.volume(),
        charged.heavy_temperature(),
        charged.electron_temperature(),
        1e-12,
        || false,
    )
    .unwrap();
    assert!(matches!(
        model.evaluate_at_state(&foreign, || false),
        Err(CoulombErrorV1::Invalid(_))
    ));
    assert!(matches!(
        model.evaluate_at_state(&state(1e40, 10_000.0, 20_000.0), || false),
        Err(CoulombErrorV1::Domain(_))
    ));
    assert_eq!(
        model.evaluate_at_state(&charged, || true),
        Err(CoulombErrorV1::Cancelled)
    );
    let mut calls = 0;
    assert_eq!(
        model.evaluate_at_state(&charged, || {
            calls += 1;
            calls >= 7
        }),
        Err(CoulombErrorV1::Cancelled)
    );
    assert_eq!(calls, 7);
}

#[test]
fn coulomb_no_carriers_is_exactly_empty_and_cancellation_returns_no_summary() {
    let neutral = state(0.0, 10_000.0, 20_000.0);
    for screening in [
        CoulombScreeningV1::StaticAllChargedDebye,
        CoulombScreeningV1::ElectronDebye,
    ] {
        let model = CoulombModelV1::new(&neutral, screening, limits(), || false).unwrap();
        let empty = model
            .evaluate_with_pairs(
                neutral.heavy_temperature(),
                neutral.electron_temperature(),
                |_| panic!("no charged pair"),
                || false,
            )
            .unwrap();
        assert_eq!(empty.pair_count(), 0);
        assert_eq!(empty.electron_heavy_thermal_conductance_w_m3_k(), 0.0);
        assert_eq!(empty.debye_length_m(), None);
        assert_eq!(empty.minimum_coulomb_logarithm(), None);
    }
    let initial = state(1e20, 10_000.0, 20_000.0);
    assert_eq!(
        CoulombModelV1::new(
            &initial,
            CoulombScreeningV1::StaticAllChargedDebye,
            limits(),
            || true
        ),
        Err(CoulombErrorV1::Cancelled)
    );
    let model = bind(&initial);
    let mut calls = 0;
    let cancelled = model.evaluate(
        initial.heavy_temperature(),
        initial.electron_temperature(),
        || {
            calls += 1;
            calls >= 7
        },
    );
    assert_eq!(cancelled, Err(CoulombErrorV1::Cancelled));
    assert_eq!(calls, 7);
}

#[test]
fn coulomb_classical_weak_coupling_and_quantum_bounds_are_refusals_not_clamps() {
    let initial = state(1e20, 10_000.0, 20_000.0);
    let mut unrepresentable_drift = limits();
    unrepresentable_drift.maximum_relative_drift_thermal_ratio = 1e-300;
    assert!(matches!(
        CoulombModelV1::new(
            &initial,
            CoulombScreeningV1::StaticAllChargedDebye,
            unrepresentable_drift,
            || false
        ),
        Err(CoulombErrorV1::Arithmetic(_))
    ));
    let mut no_quantum = limits();
    no_quantum.maximum_quantum_to_classical_ratio = 0.1;
    assert!(matches!(
        CoulombModelV1::new(
            &initial,
            CoulombScreeningV1::StaticAllChargedDebye,
            no_quantum,
            || false
        ),
        Err(CoulombErrorV1::Domain(_))
    ));
    let mut many = limits();
    many.minimum_particles_per_debye_sphere = 1000.0;
    assert!(matches!(
        CoulombModelV1::new(
            &initial,
            CoulombScreeningV1::StaticAllChargedDebye,
            many,
            || false
        ),
        Err(CoulombErrorV1::Domain(_))
    ));
    let mut large_log = limits();
    large_log.minimum_coulomb_logarithm = 10.0;
    assert!(matches!(
        CoulombModelV1::new(
            &initial,
            CoulombScreeningV1::StaticAllChargedDebye,
            large_log,
            || false
        ),
        Err(CoulombErrorV1::Domain(_))
    ));
    let mut weak = limits();
    weak.maximum_coulomb_coupling = 1e-4;
    assert!(matches!(
        CoulombModelV1::new(
            &initial,
            CoulombScreeningV1::StaticAllChargedDebye,
            weak,
            || false
        ),
        Err(CoulombErrorV1::Domain(_))
    ));
    assert!(matches!(
        CoulombModelV1::new(
            &state(1e40, 10_000.0, 20_000.0),
            CoulombScreeningV1::StaticAllChargedDebye,
            limits(),
            || false
        ),
        Err(CoulombErrorV1::Domain(_))
    ));
    assert!(matches!(
        CoulombModelV1::new(
            &state(1e20, 10_000.0, 1e9),
            CoulombScreeningV1::StaticAllChargedDebye,
            limits(),
            || false
        ),
        Err(CoulombErrorV1::Domain(_))
    ));
    let mut invalid = limits();
    invalid.minimum_coulomb_logarithm = 1.0;
    assert!(matches!(
        CoulombModelV1::new(
            &initial,
            CoulombScreeningV1::StaticAllChargedDebye,
            invalid,
            || false
        ),
        Err(CoulombErrorV1::Invalid(_))
    ));
}

fn quantum_limits() -> CoulombLimitsV1 {
    // 10^8 K electrons have kT/(m_e c²)=0.01686. This explicit bound admits
    // that nonrelativistic leading-order example, not arbitrary hot plasma.
    CoulombLimitsV1 {
        maximum_thermal_energy_rest_energy_ratio: 0.03,
        ..limits()
    }
}

#[test]
fn quantum_cutoff_admits_hot_electrons_and_matches_independent_si_moments() {
    let cold = state(1e20, 10_000.0, 20_000.0);
    let classical = bind(&cold);
    assert_eq!(
        classical.cutoff(),
        CoulombShortRangeCutoffV1::ClassicalClosestApproach
    );
    let explicit_classical = CoulombModelV1::new_with_cutoff(
        &cold,
        classical.screening(),
        CoulombShortRangeCutoffV1::ClassicalClosestApproach,
        limits(),
        || false,
    )
    .unwrap();
    assert_eq!(classical, explicit_classical);
    let cold_quantum = CoulombModelV1::new_with_cutoff(
        &cold,
        classical.screening(),
        CoulombShortRangeCutoffV1::QuantumDiffraction,
        limits(),
        || false,
    )
    .unwrap();
    // Below the diffraction crossover, the exact retained pair math agrees.
    assert_eq!(pair(&classical), pair(&cold_quantum));

    let hot = state(1e20, 8e7, 1e8);
    for (screening, debye, logarithm, friction, thermal) in [
        (
            CoulombScreeningV1::StaticAllChargedDebye,
            4.600_598_713_389_44e-5,
            17.247_312_063_759_096,
            5.707_510_844_511_341e-7,
            0.014_125_929_597_132_359,
        ),
        (
            CoulombScreeningV1::ElectronDebye,
            6.900_898_070_084_16e-5,
            17.652_777_171_867_26,
            5.841_688_071_260_836e-7,
            0.014_458_014_477_956_465,
        ),
    ] {
        assert!(matches!(
            CoulombModelV1::new(&hot, screening, quantum_limits(), || false),
            Err(CoulombErrorV1::Domain(
                "quantum diffraction exceeds classical pair-scattering limit"
            ))
        ));
        let model = CoulombModelV1::new_with_cutoff(
            &hot,
            screening,
            CoulombShortRangeCutoffV1::QuantumDiffraction,
            quantum_limits(),
            || false,
        )
        .unwrap();
        assert_eq!(
            model.cutoff(),
            CoulombShortRangeCutoffV1::QuantumDiffraction
        );
        let (summary, actual_pair) = pair(&model);
        // 70-digit Decimal reconstruction: b90=1.671191400623426862e-13 m,
        // bq=1.487311365533546238e-12 m. Use bq, not the classical b90.
        close(summary.debye_length_m().unwrap(), debye, 5e-14);
        close(actual_pair.coulomb_logarithm(), logarithm, 5e-14);
        close(actual_pair.friction_coefficient_kg_m3_s(), friction, 9e-14);
        close(actual_pair.thermal_conductance_w_m3_k(), thermal, 9e-14);
        close(
            summary.maximum_quantum_to_classical_ratio(),
            8.899_706_909_565_922,
            5e-14,
        );
        assert_eq!(
            summary,
            model
                .evaluate_point(
                    hot.amounts(),
                    hot.volume(),
                    hot.heavy_temperature(),
                    hot.electron_temperature(),
                    || false,
                )
                .unwrap()
        );
        let w2 = 0.25 * actual_pair.maximum_relative_speed_squared_m2_s2();
        let heat = actual_pair.drift_heating(w2).unwrap();
        close(
            heat.first_power_density_w_m3() + heat.second_power_density_w_m3(),
            friction * w2,
            9e-14,
        );
    }
}

#[test]
fn quantum_cutoff_preserves_all_other_physical_and_original_source_gates() {
    let hot = state(1e20, 8e7, 1e8);
    let construct = |bounds| {
        CoulombModelV1::new_with_cutoff(
            &hot,
            CoulombScreeningV1::StaticAllChargedDebye,
            CoulombShortRangeCutoffV1::QuantumDiffraction,
            bounds,
            || false,
        )
    };
    let model = construct(quantum_limits()).unwrap();
    let summary = pair(&model).0;
    let mut bounds = [quantum_limits(); 5];
    bounds[0].minimum_coulomb_logarithm = summary.minimum_coulomb_logarithm().unwrap() + 1.0;
    bounds[1].minimum_particles_per_debye_sphere =
        2.0 * summary.particles_per_debye_sphere().unwrap();
    bounds[2].maximum_coulomb_coupling = 0.5 * summary.maximum_coulomb_coupling();
    bounds[3].maximum_thermal_phase_space_density =
        0.5 * summary.maximum_thermal_phase_space_density();
    bounds[4].maximum_thermal_energy_rest_energy_ratio =
        0.5 * summary.maximum_thermal_energy_rest_energy_ratio();
    for restrictive in bounds {
        assert!(matches!(
            construct(restrictive),
            Err(CoulombErrorV1::Domain(_))
        ));
    }
    assert!(matches!(
        model.evaluate(Temperature::new(99.0), hot.electron_temperature(), || false),
        Err(CoulombErrorV1::Species { .. })
    ));
    assert!(matches!(
        model.evaluate(hot.heavy_temperature(), Temperature::new(1e9), || false),
        Err(CoulombErrorV1::Domain(_))
    ));
    let mut invalid = quantum_limits();
    invalid.maximum_quantum_to_classical_ratio = f64::NAN;
    assert!(matches!(
        construct(invalid),
        Err(CoulombErrorV1::Invalid(_))
    ));
    assert_eq!(
        model.evaluate(hot.heavy_temperature(), hot.electron_temperature(), || true),
        Err(CoulombErrorV1::Cancelled)
    );
}

#[test]
fn coulomb_trial_refreshes_fusion_ash_with_exact_original_masses_and_charges() {
    let ids = vec![id("D+"), id("He++"), id("T+"), id("e")];
    let bank = IdealGasEquilibriumV1::new(
        ElementalMatrix::new(
            vec![ElementId::new("H").unwrap(), ElementId::new("He").unwrap()],
            ids.clone(),
            vec![vec![1, 0, 1, 0], vec![0, 1, 0, 0]],
        )
        .unwrap(),
        ChargeVector::new(ids, vec![1, 2, 1, -1]).unwrap(),
        vec![
            source("D+", 2.0 * MP),
            source("He++", 4.0 * MP),
            source("T+", 3.0 * MP),
            source("e", ELECTRON_MASS_KG),
        ],
    )
    .unwrap();
    let initial_amounts = [1e20, 0.0, 1e20, 2e20].map(|n| Amount::new(n / AVOGADRO_PER_MOL));
    let final_amounts = [8e19, 2e19, 8e19, 2e20].map(|n| Amount::new(n / AVOGADRO_PER_MOL));
    let initial = FrozenIonizedGasStateV1::new(
        &bank,
        &initial_amounts,
        Volume::new(1.0),
        Temperature::new(8e7),
        Temperature::new(1e8),
        1e-12,
        || false,
    )
    .unwrap();
    let actual = FrozenIonizedGasStateV1::new(
        &bank,
        &final_amounts,
        Volume::new(2.5),
        Temperature::new(9e7),
        Temperature::new(8e7),
        1e-12,
        || false,
    )
    .unwrap();
    for screening in [
        CoulombScreeningV1::StaticAllChargedDebye,
        CoulombScreeningV1::ElectronDebye,
    ] {
        let model = CoulombModelV1::new_with_cutoff(
            &initial,
            screening,
            CoulombShortRangeCutoffV1::QuantumDiffraction,
            quantum_limits(),
            || false,
        )
        .unwrap();
        let original = model.clone();
        let trial = model
            .evaluate_point(
                actual.amounts(),
                actual.volume(),
                actual.heavy_temperature(),
                actual.electron_temperature(),
                || false,
            )
            .unwrap();
        assert_eq!(trial.pair_count(), 6);
        assert_eq!(trial, model.evaluate_at_state(&actual, || false).unwrap());
        let rebound = CoulombModelV1::new_with_cutoff(
            &actual,
            screening,
            CoulombShortRangeCutoffV1::QuantumDiffraction,
            quantum_limits(),
            || false,
        )
        .unwrap();
        assert_eq!(
            trial,
            rebound
                .evaluate(
                    actual.heavy_temperature(),
                    actual.electron_temperature(),
                    || false
                )
                .unwrap()
        );
        assert_eq!(model, original);
        assert!(model.is_bound_to_populations(&initial));
        assert!(!model.is_bound_to_populations(&actual));
        assert_eq!(
            model
                .evaluate_point(
                    initial.amounts(),
                    initial.volume(),
                    initial.heavy_temperature(),
                    initial.electron_temperature(),
                    || false,
                )
                .unwrap()
                .pair_count(),
            3
        );
    }
}

#[test]
fn coulomb_trial_rejects_invalid_axes_populations_volume_temperatures_and_cancellation() {
    let initial = state(1e20, 10_000.0, 20_000.0);
    let model = bind(&initial);
    let eval = |amounts: &[Amount], volume, th, te| {
        model.evaluate_point(amounts, volume, th, te, || false)
    };
    let (v, th, te) = (
        initial.volume(),
        initial.heavy_temperature(),
        initial.electron_temperature(),
    );
    assert!(matches!(
        eval(&initial.amounts()[1..], v, th, te),
        Err(CoulombErrorV1::Invalid(_))
    ));
    for bad in [-1.0, f64::NAN, f64::INFINITY] {
        let mut amounts = initial.amounts().to_vec();
        amounts[0] = Amount::new(bad); // Neutral rows are validated as well.
        assert!(matches!(
            eval(&amounts, v, th, te),
            Err(CoulombErrorV1::Invalid(_))
        ));
    }
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            eval(initial.amounts(), Volume::new(bad), th, te),
            Err(CoulombErrorV1::Invalid(_))
        ));
        assert!(matches!(
            eval(initial.amounts(), v, Temperature::new(bad), te),
            Err(CoulombErrorV1::Invalid(_))
        ));
        assert!(matches!(
            eval(initial.amounts(), v, th, Temperature::new(bad)),
            Err(CoulombErrorV1::Invalid(_))
        ));
    }
    assert!(matches!(
        eval(initial.amounts(), v, Temperature::new(99.0), te),
        Err(CoulombErrorV1::Species { .. })
    ));
    let mut calls = 0;
    assert_eq!(
        model.evaluate_point(initial.amounts(), v, th, te, || {
            calls += 1;
            calls == 7
        }),
        Err(CoulombErrorV1::Cancelled)
    );
    assert_eq!(calls, 7);
}
