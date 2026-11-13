//! G1 source-bound primary creep, reciprocal heat and retained shell history.
//! The independent uniaxial oracle uses authored coefficients, not measured
//! creep or service-life calibration.

use super::*;
use fs_material::SmallStrainLaw;
use fs_material::creep::GAS_CONSTANT_J_MOL_K;
use fs_material::creep::recovery::CreepRecoveryParameters;
use fs_material::creep::recovery::source::{
    CREEP_RECOVERY_ACTIVATION_ENERGY_PROPERTY, CREEP_RECOVERY_REFERENCE_RATE_PROPERTY,
    CREEP_RECOVERY_REFERENCE_TEMPERATURE_PROPERTY, RecoveringNortonSource,
};
use fs_material::creep::source::*;
use fs_material::creep::thermal::FeedbackError as CreepFeedbackError;
use fs_material::creep::thermal::expansion::{ExpandingCreepSource, ThermalCreepState};
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::section::thermal::source::{
    CaloricCreepSource, ThermalSectionSourceError, ThermalSectionSourceFamily,
};
use fs_material::state_point::{
    INVERSE_TEMPERATURE_DIMS, LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point,
};
use fs_qty::QuantitySpec;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

const CREEP_CP: f64 = 2.0;
const STRESS: f64 = 2.0e6;
const RATE: f64 = 0.15;
const PRIMARY_STRAIN: f64 = 0.004;
const ACTIVATION_K: f64 = 6000.0;
const EXPANSION: f64 = 1.0e-6;
const DT: f64 = 0.01;

fn temperature_spec() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}

fn temperature_point(temperature_k: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", temperature_spec(), temperature_k)
        .unwrap()
}

fn fixture() -> Fixture {
    fixture_with_recovery(None)
}

fn fixture_with_recovery(recovery: Option<CreepRecoveryParameters>) -> Fixture {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, E),
        (POISSON_RATIO_PROPERTY, Dims::NONE, 0.2),
        (DENSITY_PROPERTY, Density::DIMS, RHO),
        (CP, CP_DIMS, CREEP_CP),
        (CREEP_REFERENCE_STRESS_PROPERTY, Pressure::DIMS, STRESS),
        (CREEP_REFERENCE_RATE_PROPERTY, CREEP_RATE_DIMS, RATE),
        (CREEP_STRESS_EXPONENT_PROPERTY, Dims::NONE, 2.0),
        (
            CREEP_ACTIVATION_ENERGY_PROPERTY,
            CREEP_ACTIVATION_ENERGY_DIMS,
            ACTIVATION_K * GAS_CONSTANT_J_MOL_K,
        ),
        (
            CREEP_REFERENCE_TEMPERATURE_PROPERTY,
            Temperature::DIMS,
            300.0,
        ),
        (
            CREEP_HARDENING_REFERENCE_STRAIN_PROPERTY,
            Dims::NONE,
            PRIMARY_STRAIN,
        ),
        (CREEP_HARDENING_EXPONENT_PROPERTY, Dims::NONE, 1.0),
        (CREEP_RESIDUAL_RATE_FRACTION_PROPERTY, Dims::NONE, 0.0),
        (
            LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
            INVERSE_TEMPERATURE_DIMS,
            EXPANSION,
        ),
    ] {
        let mut claim = scalar(name, dims, value);
        claim.provenance.source = "manufactured primary-creep shell".into();
        claim.validity =
            ValidityDomain::unconstrained().with_quantity("T", temperature_spec(), 300.0, 320.0);
        claims.insert_claim(claim).unwrap();
    }
    if let Some(parameters) = recovery {
        for (name, dims, value) in [
            (
                CREEP_RECOVERY_REFERENCE_RATE_PROPERTY,
                CREEP_RATE_DIMS,
                parameters.reference_rate_s,
            ),
            (
                CREEP_RECOVERY_ACTIVATION_ENERGY_PROPERTY,
                CREEP_ACTIVATION_ENERGY_DIMS,
                parameters.activation_energy_j_mol,
            ),
            (
                CREEP_RECOVERY_REFERENCE_TEMPERATURE_PROPERTY,
                Temperature::DIMS,
                parameters.reference_temperature_k,
            ),
        ] {
            let mut claim = scalar(name, dims, value);
            claim.provenance.source = "manufactured primary-creep recovery shell".into();
            claim.validity = ValidityDomain::unconstrained().with_quantity(
                "T",
                temperature_spec(),
                300.0,
                320.0,
            );
            claims.insert_claim(claim).unwrap();
        }
    }
    // This pure-creep card deliberately supplies no yield or plastic-hardening
    // property. Preparing a caloric section must not fabricate a J2 material.
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured creeping solid".into(),
            phase: "solid".into(),
            process: "test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let caloric = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new(CP, CP_DIMS),
            density: &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
            temperature_axis: "T",
            lower: &temperature_point(300.0),
            upper: &temperature_point(320.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-5,
    )
    .unwrap();
    Fixture {
        card,
        caloric,
        expansion: true,
    }
}

fn expanding_source(fixture: &Fixture) -> ExpandingCreepSource<NortonCreepSource<'_>> {
    let mechanical = NortonCreepSource::new_with_strain_hardening(
        &fixture.card,
        NortonSourceQuery {
            lower: &temperature_point(300.0),
            upper: &temperature_point(320.0),
            temperature_axis: "T",
            selection: MaterialPropertySelection::SingleClaimOnly,
            strain_limit: 0.04,
        },
        || false,
    )
    .unwrap();
    ExpandingCreepSource::new(
        &fixture.card,
        mechanical,
        &expansion_reference(fixture),
        || false,
    )
    .unwrap()
}

fn expansion_reference(fixture: &Fixture) -> fs_material::state_point::ResolvedMaterialStatePoint {
    resolve_material_state_point(
        &fixture.card,
        &temperature_point(300.0),
        &[ScalarPropertyRequirement::try_new(
            LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
            INVERSE_TEMPERATURE_DIMS,
            ScalarAdmissibility::Finite,
        )
        .unwrap()],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

// At fixed uniaxial stress and endpoint T the primary law is exactly
// dc*(a+old_alpha+dc)=a*dt*rate_ref*(sigma/s_ref)^2*Arrhenius(T).
// This stable quadratic root is independent of the tensor return map.
fn creep_increment(old_alpha: f64, temperature_k: f64, stress_pa: f64) -> f64 {
    let rate = RATE
        * (stress_pa / STRESS).powi(2)
        * (ACTIVATION_K * (1.0 / 300.0 - 1.0 / temperature_k)).exp();
    let rhs = PRIMARY_STRAIN * DT * rate;
    let shifted = PRIMARY_STRAIN + old_alpha;
    2.0 * rhs / (shifted + (shifted * shifted + 4.0 * rhs).sqrt())
}

fn endpoint_temperature(
    old_alpha: f64,
    old_temperature_k: f64,
    old_stress_pa: f64,
    stress_pa: f64,
) -> f64 {
    // Constant isotropic expansion gives s_mechanical=alpha*trace(sigma).
    // Therefore sigma:Delta(epsilon_th)-Delta(T*s) is exactly
    // -alpha*T_old*(sigma_new-sigma_old) on this free-sided axial path.
    let reversible = -EXPANSION * old_temperature_k * (stress_pa - old_stress_pa);
    let residual = |temperature| {
        RHO * CREEP_CP * (temperature - old_temperature_k)
            - stress_pa * creep_increment(old_alpha, temperature, stress_pa)
            - reversible
    };
    let (mut lower, mut upper) = (old_temperature_k - 2.0, old_temperature_k + 4.0);
    assert!(residual(lower) < 0.0 && residual(upper) > 0.0);
    for _ in 0..80 {
        let mid = 0.5 * (lower + upper);
        if residual(mid) > 0.0 {
            upper = mid;
        } else {
            lower = mid;
        }
    }
    0.5 * (lower + upper)
}

#[test]
fn g0_caloric_creep_section_inspects_full_history_and_rejects_a_different_thermal_state() {
    let fixture = fixture();
    let source = CaloricCreepSource::new(
        &fixture.card,
        expanding_source(&fixture),
        &fixture.caloric,
        || false,
    )
    .unwrap();
    let tolerances = EnergyTolerances {
        energy_j_m3: 1e-5,
        yield_pa: 1e-3,
    };
    let virgin = source
        .resolve_section_step(0.0, 0.0, DT)
        .unwrap()
        .law
        .initial_state();
    let strain = [0.003, -0.0006, -0.0006, 0.0, 0.0, 0.0];
    let h = CREEP_CP; // 301 K in this original caloric chart.
    let accepted = source
        .section_increment(0.0, h, &[0.0; 6], &virgin, &strain, tolerances, DT)
        .unwrap()
        .state;
    let retained = accepted.clone();
    assert!(accepted.mechanical.equivalent_creep_strain > 0.0);
    source
        .validate_section_history(h, &strain, &accepted, 1e-3)
        .unwrap();
    source
        .validate_section_history(h, &strain, &accepted, 1e-3)
        .unwrap();
    assert_eq!(accepted, retained);
    // The same strained point really would creep during another interval;
    // accepted-history inspection is not a unit-duration return-map probe.
    let held = source
        .section_increment(h, h, &strain, &accepted, &strain, tolerances, DT)
        .unwrap();
    assert!(
        held.state.mechanical.equivalent_creep_strain > accepted.mechanical.equivalent_creep_strain
    );
    let mismatched = ThermalCreepState::from_parts(accepted.mechanical.clone(), 302.0);
    assert!(matches!(
        source.validate_section_history(h, &strain, &mismatched, 1e-3),
        Err(ThermalSectionSourceError::Creep(
            CreepFeedbackError::Source(CreepSourceError::Invalid(_))
        ))
    ));
    assert!(
        source
            .section_increment(h, h, &strain, &mismatched, &strain, tolerances, DT)
            .is_err()
    );
    assert_eq!(accepted, retained);
    let other_card = Fixture::new(CREEP_CP);
    assert!(
        CaloricCreepSource::new(
            &other_card.card,
            expanding_source(&fixture),
            &fixture.caloric,
            || false,
        )
        .is_err()
    );
}

#[test]
fn g1_primary_creep_shell_self_heats_and_retains_history_through_hold_unload_reload() {
    context(|cx| {
        let fixture = fixture();
        let mesh = mesh();
        let mut thermal = thermal(&fixture, &[0.0], 0.0);
        let source = CaloricCreepSource::new(
            &fixture.card,
            expanding_source(&fixture),
            &fixture.caloric,
            || false,
        )
        .unwrap();
        let layers = vec![vec![ShellThermoplasticLayer {
            cell: CellId(100),
            layer: PlaneStressLayer {
                bottom_z_m: -THICKNESS / 2.0,
                top_z_m: THICKNESS / 2.0,
                subdivisions: 2,
                law: source,
            },
        }]];
        let mut shell =
            ThermoplasticShell::try_new(&thermal, &mesh, layers, discretization(), cx).unwrap();
        let supports = membrane_supports();
        let mut old_alpha = 0.0;
        let mut old_temperature = 300.0;
        let mut old_stress = 0.0;
        let mut previous_loaded_increment = f64::INFINITY;
        for (interval, stress) in [STRESS, STRESS, 0.0, STRESS].into_iter().enumerate() {
            let target = (interval + 1) as f64 * DT;
            let old_solution = shell.solution().clone();
            let old_cells = thermal.cells().clone();
            let expected_temperature =
                endpoint_temperature(old_alpha, old_temperature, old_stress, stress);
            let increment = creep_increment(old_alpha, expected_temperature, stress);
            let next_alpha = old_alpha + increment;
            let reversible = -EXPANSION * old_temperature * (stress - old_stress);
            let step = shell
                .advance(
                    &thermal,
                    &(),
                    target,
                    &supports,
                    &[(6, stress * VOLUME / L)],
                    IncrementSettings::default(),
                    thermal_config(),
                    config(),
                    cx,
                )
                .unwrap();
            let actual_temperature = step.thermal().temperature_k(CellId(100)).unwrap();
            close(actual_temperature, expected_temperature, 5e-8);
            assert!(actual_temperature > old_temperature);
            let free = EXPANSION * (expected_temperature - 300.0);
            let axial = stress / E + next_alpha + free;
            let transverse = -0.2 * stress / E - 0.5 * next_alpha + free;
            close(step.shell().solution().coordinates[1][0], L * axial, 3e-11);
            close(
                step.shell().solution().coordinates[2][1],
                L * transverse,
                3e-11,
            );
            for section in &step.shell().solution().section_states {
                for point in &section.points {
                    let history = &point.history.mechanical;
                    close(history.equivalent_creep_strain, next_alpha, 3e-9);
                    close(history.creep_strain[0], next_alpha, 3e-9);
                    close(history.creep_strain[1], -0.5 * next_alpha, 3e-9);
                    close(history.creep_strain[2], -0.5 * next_alpha, 3e-9);
                    close(point.strain[2], transverse, 3e-9);
                    close(point.history.temperature_k(), expected_temperature, 5e-8);
                }
            }
            let energy = &step.report().mechanical;
            close(energy.creep_heat_j, VOLUME * stress * increment, 3e-10);
            close(energy.physical_heat_j, energy.creep_heat_j, 1e-12);
            close(energy.reversible_heat_j, VOLUME * reversible, 3e-10);
            let cell_heat = step.report().cell_heat_j[&CellId(100)];
            close(cell_heat.creep_j, energy.creep_heat_j, 1e-12);
            close(
                RHO * VOLUME * CREEP_CP * (actual_temperature - old_temperature),
                cell_heat.total_j(),
                3e-10,
            );
            assert!(step.report().combined_residual_j.abs() <= config().energy_tolerance_j);
            close(step.report().external_heat_j, 0.0, 1e-12);
            assert_eq!(step.shell().time_s(), target);
            assert_eq!(step.thermal().time_s(), target);
            assert_eq!(shell.solution(), &old_solution);
            assert_eq!(thermal.cells(), &old_cells);
            if stress > 0.0 {
                // Heating accelerates this exact retained-history step, while
                // primary resistance still slows the repeated loaded holds.
                assert!(increment > creep_increment(old_alpha, old_temperature, stress));
                assert!(increment < previous_loaded_increment);
                previous_loaded_increment = increment;
            } else {
                close(energy.physical_heat_j, 0.0, 1e-12);
                assert!(energy.reversible_heat_j > 0.0);
            }
            if interval == 0 {
                assert!(energy.material_numerical_loss_j > 0.0);
                assert!(energy.reversible_heat_j < 0.0);
            }
            let (next_shell, next_thermal, _, _, _) = step.into_parts();
            shell = next_shell;
            thermal = next_thermal;
            old_alpha = next_alpha;
            old_temperature = expected_temperature;
            old_stress = stress;
        }
        // A late transport observer refusal cannot age any material point or
        // publish the temperature that would have accelerated the next hold.
        let old_solution = shell.solution().clone();
        let old_cells = thermal.cells().clone();
        assert!(
            shell
                .advance(
                    &thermal,
                    &RefuseAfterHeat,
                    5.0 * DT,
                    &supports,
                    &[(6, STRESS * VOLUME / L)],
                    IncrementSettings::default(),
                    thermal_config(),
                    config(),
                    cx,
                )
                .is_err()
        );
        assert_eq!(shell.solution(), &old_solution);
        assert_eq!(thermal.cells(), &old_cells);
    });
}

#[test]
fn g1_shell_recovers_primary_resistance_through_an_unloaded_thermal_dwell_without_healing_creep() {
    let recovery = CreepRecoveryParameters {
        reference_rate_s: 100.0,
        activation_energy_j_mol: 3000.0 * GAS_CONSTANT_J_MOL_K,
        reference_temperature_k: 300.0,
    };
    let fixture = fixture_with_recovery(Some(recovery));
    context(|cx| {
        let mesh = mesh();
        let mut thermal = thermal(&fixture, &[0.0], 0.0);
        let mechanical = RecoveringNortonSource::new(
            &fixture.card,
            NortonSourceQuery {
                lower: &temperature_point(300.0),
                upper: &temperature_point(320.0),
                temperature_axis: "T",
                selection: MaterialPropertySelection::SingleClaimOnly,
                strain_limit: 0.04,
            },
            || false,
        )
        .unwrap();
        let source = CaloricCreepSource::new(
            &fixture.card,
            ExpandingCreepSource::new(
                &fixture.card,
                mechanical,
                &expansion_reference(&fixture),
                || false,
            )
            .unwrap(),
            &fixture.caloric,
            || false,
        )
        .unwrap();
        let layers = vec![vec![ShellThermoplasticLayer {
            cell: CellId(100),
            layer: PlaneStressLayer {
                bottom_z_m: -THICKNESS / 2.0,
                top_z_m: THICKNESS / 2.0,
                subdivisions: 2,
                law: source,
            },
        }]];
        let mut shell =
            ThermoplasticShell::try_new(&thermal, &mesh, layers, discretization(), cx).unwrap();
        let supports = membrane_supports();
        let (mut alpha, mut resistance, mut temperature, mut old_stress, mut time) =
            (0.0, 0.0, 300.0, 0.0, 0.0);
        for (stress, dt) in [
            (STRESS, DT),
            (STRESS, DT),
            (0.0, 4.0 * DT),
            (0.0, 2.0 * DT),
            (STRESS, DT),
        ] {
            let recovery_at = |t: f64| {
                recovery.reference_rate_s
                    * (recovery.activation_energy_j_mol / GAS_CONSTANT_J_MOL_K
                        * (1.0 / recovery.reference_temperature_k - 1.0 / t))
                        .exp()
            };
            let increment_at = |t: f64| {
                // Eliminating r1=(r0+dc)/(1+k dt) gives this independent
                // quadratic at prescribed stress and endpoint temperature.
                let scaled_strain = PRIMARY_STRAIN * (1.0 + dt * recovery_at(t));
                let rate = RATE
                    * (stress / STRESS).powi(2)
                    * (ACTIVATION_K * (1.0 / 300.0 - 1.0 / t)).exp();
                let rhs = scaled_strain * dt * rate;
                let shifted = scaled_strain + resistance;
                2.0 * rhs / (shifted + (shifted * shifted + 4.0 * rhs).sqrt())
            };
            let reversible = -EXPANSION * temperature * (stress - old_stress);
            let residual =
                |t| RHO * CREEP_CP * (t - temperature) - stress * increment_at(t) - reversible;
            let (mut lower, mut upper) = (temperature - 2.0, temperature + 4.0);
            assert!(residual(lower) < 0.0 && residual(upper) > 0.0);
            for _ in 0..80 {
                let mid = 0.5 * (lower + upper);
                if residual(mid) > 0.0 {
                    upper = mid;
                } else {
                    lower = mid;
                }
            }
            let next_temperature = 0.5 * (lower + upper);
            let dc = increment_at(next_temperature);
            let next_resistance = (resistance + dc) / (1.0 + dt * recovery_at(next_temperature));
            let old_solution = shell.solution().clone();
            let old_cells = thermal.cells().clone();
            let step = shell
                .advance(
                    &thermal,
                    &(),
                    time + dt,
                    &supports,
                    &[(6, stress * VOLUME / L)],
                    IncrementSettings::default(),
                    thermal_config(),
                    config(),
                    cx,
                )
                .unwrap();
            close(
                step.thermal().temperature_k(CellId(100)).unwrap(),
                next_temperature,
                5e-8,
            );
            let free = EXPANSION * (next_temperature - 300.0);
            close(
                step.shell().solution().coordinates[1][0],
                L * (stress / E + alpha + dc + free),
                3e-11,
            );
            close(
                step.shell().solution().coordinates[2][1],
                L * (-0.2 * stress / E - 0.5 * (alpha + dc) + free),
                3e-11,
            );
            for section in &step.shell().solution().section_states {
                for point in &section.points {
                    let history = &point.history.mechanical;
                    close(history.native.equivalent_creep_strain, alpha + dc, 3e-9);
                    close(history.native.creep_strain[0], alpha + dc, 3e-9);
                    close(history.resistance_strain, next_resistance, 3e-9);
                    close(point.history.temperature_k(), next_temperature, 5e-8);
                    assert!(history.resistance_strain < history.native.equivalent_creep_strain);
                }
            }
            let energy = &step.report().mechanical;
            close(energy.creep_heat_j, VOLUME * stress * dc, 3e-10);
            close(energy.physical_heat_j, energy.creep_heat_j, 1e-12);
            close(energy.reversible_heat_j, VOLUME * reversible, 3e-10);
            close(
                RHO * VOLUME * CREEP_CP * (next_temperature - temperature),
                step.report().cell_heat_j[&CellId(100)].total_j(),
                3e-10,
            );
            close(
                step.report().combined_residual_j,
                0.0,
                config().energy_tolerance_j,
            );
            if stress == 0.0 {
                assert!(next_resistance < resistance);
                close(energy.physical_heat_j, 0.0, 1e-12);
                // Unloading has reciprocal thermoelastic heat; the next dwell
                // recovers kinetic resistance without any extra caloric source.
                if old_stress == 0.0 {
                    close(energy.reversible_heat_j, 0.0, 1e-12);
                    close(next_temperature, temperature, 1e-12);
                } else {
                    assert!(energy.material_numerical_loss_j > 0.0);
                }
            }
            assert_eq!(shell.solution(), &old_solution);
            assert_eq!(thermal.cells(), &old_cells);
            assert_eq!(shell.time_s(), time);
            let (next_shell, next_thermal, _, _, _) = step.into_parts();
            shell = next_shell;
            thermal = next_thermal;
            alpha += dc;
            resistance = next_resistance;
            temperature = next_temperature;
            old_stress = stress;
            time += dt;
        }
        let old_solution = shell.solution().clone();
        let old_cells = thermal.cells().clone();
        assert!(
            shell
                .advance(
                    &thermal,
                    &RefuseAfterHeat,
                    time + DT,
                    &supports,
                    &[(6, STRESS * VOLUME / L)],
                    IncrementSettings::default(),
                    thermal_config(),
                    config(),
                    cx,
                )
                .is_err()
        );
        assert_eq!(shell.solution(), &old_solution);
        assert_eq!(thermal.cells(), &old_cells);
        assert_eq!(shell.time_s(), time);
        assert_eq!(thermal.time_s(), time);
    });
}
