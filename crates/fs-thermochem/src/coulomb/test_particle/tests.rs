use super::*;
use crate::equilibrium::composition::IdealGasEquilibriumV1;
use crate::{
    ChargeVector, ElementId, ElementalMatrix, ElementalReferenceIdV1, Nasa9RegionV1,
    Nasa9StandardStateModelV1, ReferenceEquationOfStateV1, StandardStateConventionV1,
    StandardStatePhaseV1, required_nasa9_parameter_dims,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_qty::{MolarMass, Pressure};
use std::collections::BTreeMap;

const ALPHA: f64 = 6.644_657_345e-27;
const DEUTERON: f64 = 3.343_583_776_8e-27;
const TRITON: f64 = 5.007_356_751_2e-27;
const KEV: f64 = 1.602_176_634e-16;

fn id(name: &str) -> SpeciesId {
    SpeciesId::new(name).unwrap()
}
fn index(bank: &IdealGasEquilibriumV1, name: &str) -> usize {
    bank.models()
        .iter()
        .position(|source| source.species() == &id(name))
        .unwrap()
}
fn amounts(bank: &IdealGasEquilibriumV1, particle_counts: &[(&str, f64)]) -> [Amount; 4] {
    let mut result = [Amount::new(0.0); 4];
    for &(name, count) in particle_counts {
        result[index(bank, name)] = Amount::new(count / AVOGADRO_PER_MOL);
    }
    result
}
fn source(name: &str, mass: f64, maximum_temperature: f64) -> Nasa9StandardStateModelV1 {
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
            ElementalReferenceIdV1::new("test-ion-moments").unwrap(),
        )
        .unwrap(),
        vec![
            Nasa9RegionV1::from_card(ConstitutiveModelCard {
                law: LawId("nasa9-standard-state".into()),
                law_version: 1,
                parameters,
                state_schema_version: 0,
                initial_state: InitialStatePolicy::ZeroInternalState,
                validity: ValidityDomain::unconstrained().with("T", 100.0, maximum_temperature),
                sources: vec![],
                provenance: Provenance {
                    source: "constant-Cv independent collision fixture".into(),
                    license: "analytic test fixture".into(),
                    artifact: None,
                },
            })
            .unwrap(),
        ],
    )
    .unwrap()
}
fn model(temperatures: [f64; 2], projectile_tmax: f64) -> TestParticleCoulombV1 {
    // Deliberately nonphysical names: only source masses and charge/element axes
    // identify the product and the structural electron.
    let ids = vec![
        id("product-x"),
        id("target-y"),
        id("target-z"),
        id("carrier-w"),
    ];
    let bank = IdealGasEquilibriumV1::new(
        ElementalMatrix::new(
            vec![ElementId::new("H").unwrap(), ElementId::new("He").unwrap()],
            ids.clone(),
            vec![vec![0, 1, 1, 0], vec![1, 0, 0, 0]],
        )
        .unwrap(),
        ChargeVector::new(ids, vec![2, 1, 1, -1]).unwrap(),
        vec![
            source("product-x", ALPHA, projectile_tmax),
            source("target-y", DEUTERON, 1e10),
            source("target-z", TRITON, 1e10),
            source("carrier-w", ELECTRON_MASS_KG, 1e10),
        ],
    )
    .unwrap();
    let state = FrozenIonizedGasStateV1::new(
        &bank,
        &amounts(
            &bank,
            &[("target-y", 5e19), ("target-z", 5e19), ("carrier-w", 1e20)],
        ),
        Volume::new(1.0),
        Temperature::new(temperatures[0]),
        Temperature::new(temperatures[1]),
        1e-12,
        || false,
    )
    .unwrap();
    let background = CoulombModelV1::new_with_cutoff(
        &state,
        CoulombScreeningV1::StaticAllChargedDebye,
        CoulombShortRangeCutoffV1::QuantumDiffraction,
        CoulombLimitsV1 {
            minimum_coulomb_logarithm: 3.0,
            minimum_particles_per_debye_sphere: 10.0,
            maximum_coulomb_coupling: 0.1,
            maximum_thermal_phase_space_density: 0.1,
            maximum_thermal_energy_rest_energy_ratio: 0.05,
            maximum_quantum_to_classical_ratio: 0.5,
            maximum_relative_drift_thermal_ratio: 0.1,
        },
        || false,
    )
    .unwrap();
    TestParticleCoulombV1::new(background, index(&bank, "product-x"), 0.02).unwrap()
}
fn prepared(model: &TestParticleCoulombV1) -> PreparedTestParticleCoulombV1 {
    let state = model.background().initial_state();
    model
        .prepare(
            state.heavy_temperature(),
            state.electron_temperature(),
            || false,
        )
        .unwrap()
}
fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        actual.is_finite() && (actual - expected).abs() <= relative * expected.abs(),
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

#[test]
fn alpha_moments_match_independent_si_and_incomplete_gamma_oracle() {
    // Independent Python/SciPy calculation: direct SI prefactors and
    // scipy.special.gammainc(1.5,x), not the production erf/series formulation.
    // Static all-charged lambda_D=4.291988413234599e-5 m; frozen logs
    // [D,T,e]=[18.819665195484653,18.819665195484653,16.906033519787382].
    let model = model(
        [
            10.0 * KEV / BOLTZMANN_J_PER_K,
            5.0 * KEV / BOLTZMANN_J_PER_K,
        ],
        1e10,
    );
    let point = prepared(&model);
    assert_eq!(point.target_count(), 3);
    close(
        point.minimum_coulomb_logarithm().unwrap(),
        16.906_033_519_787_382,
        3e-14,
    );
    let moments = point.evaluate(Energy::new(3500.0 * KEV)).unwrap();
    for (pool, (drift, diffusion)) in [
        (-9.665_467_498_790_603e-14, 1.548_578_618_324_872_4e-28),
        (-8.054_701_901_933_429e-12, 6.465_857_000_194_071e-27),
    ]
    .into_iter()
    .enumerate()
    {
        close(moments.mean_energy_drift_j_s[pool], drift, 3e-13);
        close(moments.energy_diffusion_j2_s[pool], diffusion, 3e-13);
    }
    let zero = point.evaluate(Energy::new(0.0)).unwrap();
    for (pool, expected) in [9.711_436_683_183_054e-13, 1.831_503_437_942_622e-14]
        .into_iter()
        .enumerate()
    {
        close(zero.mean_energy_drift_j_s[pool], expected, 3e-13);
        assert!(zero.velocity_drag_frequency_s[pool] > 0.0);
    }
    assert_eq!(zero.energy_diffusion_j2_s, [0.0, 0.0]);
}

#[test]
fn cold_high_speed_limit_has_finite_mass_energy_loss_and_einstein_diffusion() {
    let temperature = 10_000.0;
    let point = prepared(&model([temperature; 2], 1e10));
    let kt = BOLTZMANN_J_PER_K * temperature;
    let e = 3500.0 * KEV;
    let speed = (2.0 * e / ALPHA).sqrt();
    let c =
        2.0 * ELEMENTARY_CHARGE_C.powi(2) / (4.0 * core::f64::consts::PI * VACUUM_PERMITTIVITY_F_M);
    let lambda = (VACUUM_PERMITTIVITY_F_M * kt / (2e20 * ELEMENTARY_CHARGE_C.powi(2))).sqrt();
    let prefactor = 4.0 * core::f64::consts::PI * c.powi(2) * (lambda / (c / kt)).ln() / speed;
    let moments = point.evaluate(Energy::new(e)).unwrap();
    for (pool, density_over_mass) in [5e19 / DEUTERON + 5e19 / TRITON, 1e20 / ELECTRON_MASS_KG]
        .into_iter()
        .enumerate()
    {
        close(
            -moments.mean_energy_drift_j_s[pool],
            prefactor * density_over_mass,
            4e-13,
        );
        close(
            moments.energy_diffusion_j2_s[pool],
            -kt * moments.mean_energy_drift_j_s[pool],
            4e-14,
        );
    }
    close(
        -moments.mean_energy_drift_j_s[1] / (2.0 * e * moments.velocity_drag_frequency_s[1]),
        ALPHA / (ALPHA + ELECTRON_MASS_KG),
        4e-14,
    );
}

#[test]
fn detailed_balance_holds_separately_for_each_pool_across_energy_regimes() {
    let point = prepared(&model(
        [
            10.0 * KEV / BOLTZMANN_J_PER_K,
            5.0 * KEV / BOLTZMANN_J_PER_K,
        ],
        1e10,
    ));
    for energy_kev in [1e-9, 0.05, 0.5, 1.0, 1.5, 5.0, 20.0, 100.0, 3500.0] {
        let e = energy_kev * KEV;
        let h = e * 2e-5;
        let moments = point.evaluate(Energy::new(e)).unwrap();
        let above = point.evaluate(Energy::new(e + h)).unwrap();
        let below = point.evaluate(Energy::new(e - h)).unwrap();
        for pool in 0..2 {
            let derivative =
                (above.energy_diffusion_j2_s[pool] - below.energy_diffusion_j2_s[pool]) / (2.0 * h);
            let density_gradient = moments.energy_diffusion_j2_s[pool]
                * (0.5 / e - 1.0 / (BOLTZMANN_J_PER_K * point.temperatures()[pool].value()));
            let discrepancy =
                (moments.mean_energy_drift_j_s[pool] - derivative - density_gradient).abs();
            assert!(discrepancy <= 3e-8 * (derivative.abs() + density_gradient.abs()));
        }
    }
    let zero = point.evaluate(Energy::new(0.0)).unwrap();
    let tiny_e = 1e-12 * KEV;
    let tiny = point.evaluate(Energy::new(tiny_e)).unwrap();
    for pool in 0..2 {
        close(
            tiny.mean_energy_drift_j_s[pool],
            zero.mean_energy_drift_j_s[pool],
            5e-13,
        );
        close(
            1.5 * tiny.energy_diffusion_j2_s[pool] / tiny_e,
            zero.mean_energy_drift_j_s[pool],
            5e-13,
        );
    }
}

#[test]
fn thermal_tail_changes_sign_and_self_species_targets_are_not_omitted() {
    let model = model([5.0 * KEV / BOLTZMANN_J_PER_K; 2], 1e10);
    let point = prepared(&model);
    for pool in 0..2 {
        assert!(
            point
                .evaluate(Energy::new(0.1 * 5.0 * KEV))
                .unwrap()
                .mean_energy_drift_j_s[pool]
                > 0.0
        );
        assert!(
            point
                .evaluate(Energy::new(10.0 * 5.0 * KEV))
                .unwrap()
                .mean_energy_drift_j_s[pool]
                < 0.0
        );
    }
    let state = model.background().initial_state();
    let amounts = amounts(
        state.thermodynamics(),
        &[("product-x", 1e19), ("carrier-w", 2e19)],
    );
    let self_point = model
        .prepare_point(
            &amounts,
            state.volume(),
            state.heavy_temperature(),
            state.electron_temperature(),
            || false,
        )
        .unwrap();
    assert_eq!(self_point.target_count(), 2);
    assert!(
        self_point
            .evaluate(Energy::new(3500.0 * KEV))
            .unwrap()
            .energy_diffusion_j2_s[0]
            > 0.0
    );
}

#[test]
fn raw_point_refresh_retains_sources_and_accepts_explicit_empty_background() {
    let model = model([1e6, 2e6], 1e10);
    let original = model.clone();
    let state = model.background().initial_state();
    let doubled = state
        .amounts()
        .iter()
        .map(|n| Amount::new(2.0 * n.value()))
        .collect::<Vec<_>>();
    let refreshed = model
        .prepare_point(
            &doubled,
            Volume::new(2.0),
            state.heavy_temperature(),
            state.electron_temperature(),
            || false,
        )
        .unwrap();
    assert_eq!(prepared(&model), refreshed);
    let mut changed = doubled;
    changed[index(state.thermodynamics(), "target-y")] = Amount::new(0.0); // Raw, without charge closure.
    let changed_point = model
        .prepare_point(
            &changed,
            Volume::new(2.0),
            Temperature::new(3e6),
            Temperature::new(4e6),
            || false,
        )
        .unwrap();
    assert_eq!(changed_point.target_count(), 2);
    assert_eq!(
        changed_point.temperatures(),
        [Temperature::new(3e6), Temperature::new(4e6)]
    );
    assert_ne!(
        refreshed.evaluate(Energy::new(KEV)),
        changed_point.evaluate(Energy::new(KEV))
    );
    assert_eq!(model, original);
    let empty = model
        .prepare_point(
            &[Amount::new(0.0); 4],
            Volume::new(1.0),
            state.heavy_temperature(),
            state.electron_temperature(),
            || false,
        )
        .unwrap();
    assert_eq!(empty.target_count(), 0);
    assert_eq!(empty.minimum_coulomb_logarithm(), None);
    let zero = empty.evaluate(Energy::new(KEV)).unwrap();
    assert_eq!(zero.mean_energy_drift_j_s, [0.0; 2]);
    assert_eq!(zero.energy_diffusion_j2_s, [0.0; 2]);
}

#[test]
fn source_domain_cancellation_and_nonrelativistic_bounds_are_explicit() {
    let model = model([1e6, 2e6], 1e10);
    let state = model.background().initial_state();
    for index in [
        state.electron_index().unwrap(),
        state.amounts().len(),
        usize::MAX,
    ] {
        assert!(TestParticleCoulombV1::new(model.background().clone(), index, 0.02).is_err());
    }
    for limit in [0.0, -1.0, 1.0, f64::NAN, f64::INFINITY] {
        assert!(
            TestParticleCoulombV1::new(model.background().clone(), model.projectile_index(), limit)
                .is_err()
        );
    }
    assert_eq!(
        model.prepare(
            state.heavy_temperature(),
            state.electron_temperature(),
            || true
        ),
        Err(CoulombErrorV1::Cancelled)
    );
    assert!(
        model
            .prepare_point(
                &[],
                state.volume(),
                state.heavy_temperature(),
                state.electron_temperature(),
                || false
            )
            .is_err()
    );
    let point = prepared(&model);
    for e in [
        -1.0,
        f64::NAN,
        f64::INFINITY,
        point.maximum_energy().value() * 1.01,
    ] {
        assert!(point.evaluate(Energy::new(e)).is_err());
    }
    let absent_source = self::model([1e6, 2e6], 5e5);
    assert!(
        matches!(absent_source.prepare(Temperature::new(1e6), Temperature::new(2e6), || false),
        Err(CoulombErrorV1::Species { species, .. }) if species == id("product-x"))
    );
}
