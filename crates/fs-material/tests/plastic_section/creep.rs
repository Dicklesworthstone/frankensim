use super::*;
use fs_material::creep::plastic::NortonJ2;
use fs_material::creep::recovery::plastic::RecoveringNortonJ2;
use fs_material::creep::recovery::{CreepRecoveryParameters, RecoveringNortonCreep};
use fs_material::creep::{NortonCreep, NortonCreepParameters, StrainHardeningCreep};
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::section::energy::NativeSectionEnergyLaw;
use fs_material::tensor::{contract, von_mises};

const DT: f64 = 0.25;
const STRESS: f64 = 3.0;
const RATE: f64 = 0.004;
const PRIMARY_STRAIN: f64 = 0.002;

fn norton(primary: bool) -> NortonCreep {
    let law = NortonCreep::new(
        elastic(1200.0, 0.2),
        NortonCreepParameters {
            reference_stress_pa: STRESS,
            reference_rate_s: RATE,
            stress_exponent: 1.0,
            activation_energy_j_mol: 0.0,
            reference_temperature_k: 300.0,
            temperature_range_k: [300.0, 500.0],
        },
        300.0,
        DT,
    )
    .unwrap();
    if primary {
        law.with_strain_hardening(StrainHardeningCreep {
            reference_strain: PRIMARY_STRAIN,
            exponent: 1.0,
            residual_rate_fraction: 0.0,
        })
        .unwrap()
    } else {
        law
    }
}

fn gates() -> EnergyTolerances {
    EnergyTolerances {
        energy_j_m3: 1e-9,
        yield_pa: 1e-8,
    }
}

#[test]
fn constant_stress_sections_keep_primary_history_and_its_original_physical_interval() {
    for primary in [false, true] {
        let law = norton(primary);
        assert_eq!(law.physical_duration_s(), Some(DT));
        assert_eq!((&law).physical_duration_s(), Some(DT));
        let section = section(law, 2);
        let mut old = section.initial_state();
        let mut alpha = 0.0;
        let mut previous_increment = f64::INFINITY;
        for _ in 0..5 {
            // Independently solve dc = dt * rate * a/(a + alpha_old + dc).
            // The stable quadratic root also verifies the endpoint history
            // convention; a previous-history multiplier gives a different dc.
            let dc = if primary {
                let b = PRIMARY_STRAIN + alpha;
                let c = PRIMARY_STRAIN * DT * RATE;
                2.0 * c / (b + (b * b + 4.0 * c).sqrt())
            } else {
                DT * RATE
            };
            alpha += dc;
            let strain = [
                STRESS / 1200.0 + alpha,
                -0.2 * STRESS / 1200.0 - 0.5 * alpha,
                0.0,
                0.0,
                0.0,
                0.0,
            ];
            let saved = old.clone();
            let next = section.evaluate(strain, &old).unwrap();
            let energy = section
                .audit_energy(&old, &next.state, gates(), 1e-10)
                .unwrap();
            for point in &next.state.points {
                near(point.history.equivalent_creep_strain, alpha, 2e-14);
                near(point.history.creep_strain[0], alpha, 2e-14);
                near(point.strain[2], strain[1], 2e-14);
            }
            near(next.resultants[0], 0.2 * STRESS, 2e-11);
            for resultant in &next.resultants[1..] {
                near(*resultant, 0.0, 2e-11);
            }
            near(energy.physical_dissipation_j_m2, 0.2 * STRESS * dc, 1e-12);
            near(
                energy.final_state.elastic_j_m2,
                0.2 * 0.5 * STRESS * STRESS / 1200.0,
                1e-12,
            );
            assert_eq!(energy.final_state.hardening_j_m2, 0.0);
            assert_eq!(energy.final_state.kinematic_hardening_j_m2, 0.0);
            near(energy.balance_residual_j_m2, 0.0, 1e-10);
            assert_eq!(old, saved);
            if primary {
                assert!(dc < previous_increment);
                previous_increment = dc;
            }
            old = next.state;
        }
        // An accepted loaded state must still creep over the next dt. Reusing
        // the old state as its own endpoint cannot erase that physical step.
        assert!(section.audit_energy(&old, &old, gates(), 1e-10).is_err());
    }
}

#[test]
fn combined_creep_sections_preserve_the_native_elastic_hardening_and_heat_partition() {
    for primary in [false, true] {
        let creep = norton(primary);
        let plastic = J2Plasticity::new(creep.elastic().clone(), 3.0, 40.0).unwrap();
        let law = NortonJ2::new(creep, plastic).unwrap();
        assert_eq!(law.physical_duration_s(), Some(DT));
        let section = section(law.clone(), 2);
        let loaded = [0.008, -0.001, 0.002, 0.06, -0.01, 0.02];
        let mut old = section.initial_state();
        for step in 0..2 {
            let next = section.evaluate(loaded, &old).unwrap();
            let energy = section
                .audit_energy(&old, &next.state, gates(), 1e-10)
                .unwrap();
            let mut elastic = 0.0;
            let mut hardening = 0.0;
            let mut creep_heat = 0.0;
            let mut plastic_heat = 0.0;
            let mut old_elastic = 0.0;
            let mut old_hardening = 0.0;
            for ((point, before), after) in section
                .points()
                .iter()
                .zip(&old.points)
                .zip(&next.state.points)
            {
                let e = core::array::from_fn(|i| {
                    after.strain[i]
                        - after.history.creep.creep_strain[i]
                        - after.history.plastic.plastic_strain[i]
                });
                let old_e = core::array::from_fn(|i| {
                    before.strain[i]
                        - before.history.creep.creep_strain[i]
                        - before.history.plastic.plastic_strain[i]
                });
                // Storage is reconstructed from the two distinct signed
                // histories, and physical heat from the two accumulated paths.
                let stress = law.creep().elastic().stress(&e, &());
                let old_stress = law.creep().elastic().stress(&old_e, &());
                elastic += point.weight_m * 0.5 * contract(&e, &stress);
                old_elastic += point.weight_m * 0.5 * contract(&old_e, &old_stress);
                hardening += point.weight_m * 0.5 * 40.0 * after.history.plastic.alpha.powi(2);
                old_hardening += point.weight_m * 0.5 * 40.0 * before.history.plastic.alpha.powi(2);
                creep_heat += point.weight_m
                    * von_mises(&stress)
                    * (after.history.creep.equivalent_creep_strain
                        - before.history.creep.equivalent_creep_strain);
                plastic_heat += point.weight_m
                    * 3.0
                    * (after.history.plastic.alpha - before.history.plastic.alpha);
            }
            assert!(creep_heat > 0.0);
            assert!(hardening > 0.0);
            if step == 0 {
                assert!(plastic_heat > 0.0);
            } else {
                near(energy.generalized_work_j_m2, 0.0, 0.0);
                assert!(energy.stored_change_j_m2() < 0.0);
            }
            near(energy.initial.elastic_j_m2, old_elastic, 1e-12);
            near(energy.initial.hardening_j_m2, old_hardening, 1e-12);
            near(energy.final_state.elastic_j_m2, elastic, 1e-12);
            near(energy.final_state.hardening_j_m2, hardening, 1e-12);
            assert_eq!(energy.final_state.kinematic_hardening_j_m2, 0.0);
            near(
                energy.physical_dissipation_j_m2,
                creep_heat + plastic_heat,
                1e-12,
            );
            near(energy.balance_residual_j_m2, 0.0, 1e-10);
            old = next.state;
        }
    }
}

#[test]
fn recovering_sections_keep_resistance_separate_from_permanent_strain_and_heat() {
    let recovery_rate = 2.0;
    let law = RecoveringNortonCreep::new(
        norton(true),
        CreepRecoveryParameters {
            reference_rate_s: recovery_rate,
            activation_energy_j_mol: 0.0,
            reference_temperature_k: 300.0,
        },
    )
    .unwrap();
    assert_eq!(law.physical_duration_s(), Some(DT));
    let section = section(law, 2);
    let mut old = section.initial_state();
    let mut alpha = 0.0;
    let mut resistance = 0.0;
    let g = 1.0 + recovery_rate * DT;
    for _ in 0..5 {
        // At constant axial stress: dc*(a*g + r_old + dc) = dt*rate*a*g.
        // This independent quadratic checks both resistance and permanent path.
        let b: f64 = PRIMARY_STRAIN * g + resistance;
        let c = DT * RATE * PRIMARY_STRAIN * g;
        let dc = 2.0 * c / (b + (b * b + 4.0 * c).sqrt());
        alpha += dc;
        resistance = (resistance + dc) / g;
        let strain = [
            STRESS / 1200.0 + alpha,
            -0.2 * STRESS / 1200.0 - 0.5 * alpha,
            0.0,
            0.0,
            0.0,
            0.0,
        ];
        let next = section.evaluate(strain, &old).unwrap();
        let energy = section
            .audit_energy(&old, &next.state, gates(), 1e-10)
            .unwrap();
        for point in &next.state.points {
            near(point.history.native.equivalent_creep_strain, alpha, 2e-14);
            near(point.history.native.creep_strain[0], alpha, 2e-14);
            near(point.history.resistance_strain, resistance, 2e-14);
            assert!(point.history.resistance_strain < point.history.native.equivalent_creep_strain);
        }
        near(next.resultants[0], 0.2 * STRESS, 2e-11);
        near(energy.physical_dissipation_j_m2, 0.2 * STRESS * dc, 1e-12);
        assert_eq!(energy.final_state.hardening_j_m2, 0.0);
        near(energy.balance_residual_j_m2, 0.0, 1e-10);
        old = next.state;
    }
    let unloaded = section
        .evaluate([alpha, -0.5 * alpha, 0.0, 0.0, 0.0, 0.0], &old)
        .unwrap();
    let energy = section
        .audit_energy(&old, &unloaded.state, gates(), 1e-10)
        .unwrap();
    for point in &unloaded.state.points {
        near(point.history.native.equivalent_creep_strain, alpha, 2e-14);
        near(point.history.native.creep_strain[0], alpha, 2e-14);
        near(point.history.resistance_strain, resistance / g, 2e-14);
    }
    near(energy.physical_dissipation_j_m2, 0.0, 1e-12);
    assert!(energy.numerical_dissipation_j_m2 > 0.0);
    near(energy.balance_residual_j_m2, 0.0, 1e-10);
}

#[test]
fn recovering_joint_section_retains_both_permanent_histories_through_dwell_and_reload() {
    let creep = norton(true);
    let plastic = J2Plasticity::new(creep.elastic().clone(), 3.0, 40.0).unwrap();
    let law = RecoveringNortonJ2::new(
        NortonJ2::new(creep, plastic).unwrap(),
        CreepRecoveryParameters {
            reference_rate_s: 2.0,
            activation_energy_j_mol: 0.0,
            reference_temperature_k: 300.0,
        },
    )
    .unwrap();
    assert_eq!(law.physical_duration_s(), Some(DT));
    assert_eq!((&law).physical_duration_s(), Some(DT));
    let section = section(law, 2);
    let mut old = section.initial_state();
    let (mut alpha_c, mut alpha_p, mut resistance) = (0.0, 0.0_f64, 0.0);
    let g = 1.0 + 2.0 * DT;
    for (step, stress) in [3.2_f64, 3.2, 0.0, 0.0, 3.4].into_iter().enumerate() {
        // Independent uniaxial endpoint equations: the primary quadratic uses
        // recovering resistance, while yield uses permanent plastic hardening.
        let b: f64 = PRIMARY_STRAIN * g + resistance;
        let c = DT * RATE * (stress / STRESS) * PRIMARY_STRAIN * g;
        let dc = 2.0 * c / (b + (b * b + 4.0 * c).sqrt());
        let plastic = alpha_p.max((stress - 3.0) / 40.0);
        let dp = plastic - alpha_p;
        let next_resistance = (resistance + dc) / g;
        let permanent = alpha_c + dc + plastic;
        let strain = [
            stress / 1200.0 + permanent,
            -0.2 * stress / 1200.0 - 0.5 * permanent,
            0.0,
            0.0,
            0.0,
            0.0,
        ];
        let saved = old.clone();
        let next = section.evaluate(strain, &old).unwrap();
        let energy = section
            .audit_energy(&old, &next.state, gates(), 1e-10)
            .unwrap();
        assert_eq!(old, saved);
        for (before, after) in old.points.iter().zip(&next.state.points) {
            let history = &after.history;
            near(
                history.native.creep.equivalent_creep_strain,
                alpha_c + dc,
                2e-14,
            );
            near(history.native.creep.creep_strain[0], alpha_c + dc, 2e-14);
            near(history.native.plastic.alpha, plastic, 2e-14);
            near(history.native.plastic.plastic_strain[0], plastic, 2e-14);
            near(history.resistance_strain, next_resistance, 2e-14);
            near(after.strain[2], strain[1], 2e-14);
            if stress == 0.0 {
                for i in 0..6 {
                    near(
                        history.native.creep.creep_strain[i],
                        before.history.native.creep.creep_strain[i],
                        2e-14,
                    );
                    near(
                        history.native.plastic.plastic_strain[i],
                        before.history.native.plastic.plastic_strain[i],
                        2e-14,
                    );
                }
                assert!(history.resistance_strain < before.history.resistance_strain);
            }
        }
        near(next.resultants[0], 0.2 * stress, 2e-11);
        for resultant in &next.resultants[1..] {
            near(*resultant, 0.0, 2e-11);
        }
        near(
            energy.final_state.elastic_j_m2,
            0.2 * 0.5 * stress * stress / 1200.0,
            1e-12,
        );
        near(
            energy.initial.hardening_j_m2,
            0.2 * 0.5 * 40.0 * alpha_p.powi(2),
            1e-12,
        );
        near(
            energy.final_state.hardening_j_m2,
            0.2 * 0.5 * 40.0 * plastic.powi(2),
            1e-12,
        );
        near(
            energy.physical_dissipation_j_m2,
            0.2 * (stress * dc + 3.0 * dp),
            1e-12,
        );
        assert_eq!(energy.final_state.kinematic_hardening_j_m2, 0.0);
        near(energy.balance_residual_j_m2, 0.0, 1e-10);
        if step == 0 || step == 4 {
            assert!(dc > 0.0 && dp > 0.0);
        } else if stress == 0.0 {
            near(energy.physical_dissipation_j_m2, 0.0, 1e-12);
            assert!(energy.final_state.hardening_j_m2 > 0.0);
            if step == 2 {
                // Elastic unloading loses BE energy without caloric heat;
                // recovery continues during the following zero-work dwell.
                assert!(energy.numerical_dissipation_j_m2 > 0.0);
            } else {
                near(energy.generalized_work_j_m2, 0.0, 1e-12);
            }
        }
        alpha_c += dc;
        alpha_p = plastic;
        resistance = next_resistance;
        old = next.state;
    }
}
