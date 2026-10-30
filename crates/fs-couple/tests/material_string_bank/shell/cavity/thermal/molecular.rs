//! Original source caloric reservoirs on the native shell/string/FE owner.
//! Synthetic collision rates test integration, not empirical gas calibration.
use super::*;
use fs_couple::bowed_string::runtime::reciprocal::bank::cavity::{
    ShellCavityError, ShellCavityMolecularOptions,
};
use fs_material::relaxation::{
    InternalReservoir, RelaxationAbsorption, RelaxationProcess, RelaxationRateLaw, RelaxationState,
    RelaxingMedium, SpeciesRelaxation,
};
use fs_thermochem::SpeciesId;

fn source(gas: GasState, tau: Option<f64>) -> RelaxationState {
    let species = SpeciesId::new("synthetic-cavity-reservoir").unwrap();
    let processes = tau
        .map(|tau| RelaxationProcess {
            reservoir: InternalReservoir::Rotational {
                degrees_of_freedom: 2,
            },
            rate: RelaxationRateLaw::PartnerIndependent {
                p_tau_ref_pa_s: 101325.0 * tau,
                t_ref_k: 293.15,
                exponent: 0.0,
                t_min_k: 200.0,
                t_max_k: 400.0,
            },
        })
        .into_iter()
        .collect();
    let laws = vec![SpeciesRelaxation {
        species: species.clone(),
        processes,
    }];
    RelaxationAbsorption::new(laws)
        .unwrap()
        .evaluate(&RelaxingMedium::new(gas, vec![(species, 1.0)]).unwrap())
        .unwrap()
}
fn molecular_options() -> ShellCavityMolecularOptions {
    ShellCavityMolecularOptions {
        maximum_caloric_coordinates: 8,
        maximum_memory_states: 256,
    }
}
fn molecular_bank(
    cx: &Cx<'_>,
    gate: &CancelGate,
    fe: &HelmholtzFe<'_>,
    gas: GasState,
    source: RelaxationState,
) -> ReciprocalStringBank {
    bank(cx, gate, fe, gas, 700.0)
        .with_material_shell_cavity_thermal_molecular(
            fe,
            gas,
            pressure_resolution(),
            Vec::new(),
            heat_resolution(),
            source,
            molecular_options(),
            gate,
        )
        .unwrap()
}
fn close(a: C64, b: C64) {
    assert!((a - b).abs() <= 3e-8 * b.abs().max(1.0), "{a:?} != {b:?}");
}

// Direct local gas/internal-temperature elimination, BEFORE normalizing any
// capacity, diagonalizing a caloric matrix, or constructing Maxwell branches.
fn direct_bulk(source: &RelaxationState, lambda: Option<f64>, omega: f64) -> C64 {
    let gas = source.gas();
    let s = C64::new(0.0, -omega);
    let cv = source
        .processes()
        .iter()
        .fold(C64::from_re(source.frozen_cv_j_kg_k()), |sum, p| {
            sum + C64::from_re(p.heat_capacity_j_kg_k) / (C64::ONE + s.scale(p.relaxation_time_s))
        });
    match lambda {
        Some(lambda) => {
            C64::from_re(gas.pressure)
                + s.scale(gas.pressure * gas.pressure / gas.temperature)
                    / (s * cv.scale(gas.density) + C64::from_re(gas.thermal_conductivity * lambda))
        }
        None => {
            C64::from_re(gas.pressure) + C64::from_re(gas.pressure * gas.specific_gas_constant) / cv
        }
    }
}

#[test]
fn molecular_fe_bulk_matches_direct_heat_equations_including_unresolved_space() {
    with_cx(|cx, gate| {
        for pressure in [101325.0, 202650.0] {
            with_volume(2, gas(pressure), |fe| {
                for tau in [0.0025, 0.025] {
                    let source = source(gas(pressure), Some(tau));
                    let bank = molecular_bank(cx, gate, fe, gas(pressure), source.clone());
                    let law = bank
                        .source()
                        .spatial_cavity()
                        .unwrap()
                        .thermal_loss()
                        .unwrap();
                    assert_eq!(law.molecular_relaxation(), Some(&source));
                    assert!(law.molecular_complement_memory_states() > 0);
                    assert!(
                        law.memory_mode_indices()
                            .windows(2)
                            .any(|pair| pair[0] == pair[1])
                    );
                    assert_eq!(bank.system().material_memory().len(), law.memory_states());
                    let n = bank
                        .source()
                        .spatial_cavity()
                        .unwrap()
                        .coupling()
                        .cavity_modes();
                    for omega in [0.01, 0.2, 5.0, 100.0, 10000.0] {
                        let got = law.normalized_bulk_response(omega).unwrap();
                        let adiabatic = direct_bulk(&source, None, omega);
                        for i in 0..n {
                            for j in 0..n {
                                let mut gram = 0.0;
                                let mut expected = C64::ZERO;
                                for (row, lambda) in law
                                    .normalized_pressure_overlaps()
                                    .iter()
                                    .zip(&law.diffusion_modes().eigenvalues_per_m2)
                                {
                                    let weight = row[i] * row[j];
                                    gram += weight;
                                    expected = expected
                                        + direct_bulk(&source, Some(*lambda), omega).scale(weight);
                                }
                                expected = expected
                                    + adiabatic.scale(if i == j { 1.0 - gram } else { -gram });
                                close(got[i * n + j], expected);
                            }
                        }
                    }
                    assert!(law.normalized_bulk_response(f64::NAN).is_err());
                }
            });
        }
    });
}

#[test]
fn zero_reservoir_source_recovers_original_fe_thermal_trajectory_bitwise() {
    with_cx(|cx, gate| {
        with_volume(2, gas(101325.0), |fe| {
            let gas = gas(101325.0);
            let mut original = bank(cx, gate, fe, gas, 700.0)
                .with_material_shell_cavity_thermal(
                    fe,
                    pressure_resolution(),
                    Vec::new(),
                    heat_resolution(),
                    gate,
                )
                .unwrap();
            let mut joint = molecular_bank(cx, gate, fe, gas, source(gas, None));
            assert_eq!(
                original.source().body_connections(),
                joint.source().body_connections()
            );
            let a = original
                .source()
                .spatial_cavity()
                .unwrap()
                .thermal_loss()
                .unwrap();
            let b = joint
                .source()
                .spatial_cavity()
                .unwrap()
                .thermal_loss()
                .unwrap();
            assert_eq!(a.decay_rates_per_s(), b.decay_rates_per_s());
            assert_eq!(
                a.normalized_bulk_response(100.0).unwrap(),
                b.normalized_bulk_response(100.0).unwrap()
            );
            assert_eq!(b.molecular_complement_memory_states(), 0);
            let mut force = vec![0.0; joint.system().mode_count()];
            force[4] = 0.02;
            let released = FrictionDrive {
                speed_m_s: 0.0,
                normal_force_n: 0.0,
            };
            let mut scratch = vec![
                0.0;
                joint
                    .source()
                    .spatial_cavity()
                    .unwrap()
                    .observation_state_len()
            ];
            let mut pa = [0.0];
            let mut pb = [0.0];
            for _ in 0..48 {
                assert_eq!(
                    original.step(&force, released).unwrap(),
                    joint.step(&force, released).unwrap()
                );
                original
                    .spatial_cavity_pressures_into(&mut scratch, &mut pa)
                    .unwrap();
                joint
                    .spatial_cavity_pressures_into(&mut scratch, &mut pb)
                    .unwrap();
                assert_eq!(pa, pb);
            }
            assert_eq!(
                original.system().material_memory(),
                joint.system().material_memory()
            );
        })
    });
}

#[test]
fn source_reservoirs_change_native_pressure_sound_and_share_one_energy_history() {
    with_cx(|cx, gate| {
        with_volume(2, gas(101325.0), |fe| {
            let gas = gas(101325.0);
            let mut first = molecular_bank(cx, gate, fe, gas, source(gas, Some(0.025)));
            let mut second = molecular_bank(cx, gate, fe, gas, source(gas, Some(0.0025)));
            let mut force = vec![0.0; first.system().mode_count()];
            force[4] = 0.02;
            let released = FrictionDrive {
                speed_m_s: 0.0,
                normal_force_n: 0.0,
            };
            let mut scratch = vec![
                0.0;
                first
                    .source()
                    .spatial_cavity()
                    .unwrap()
                    .observation_state_len()
            ];
            let mut a = [0.0];
            let mut b = [0.0];
            let mut pressure_difference = 0.0;
            let mut sound_difference = 0.0;
            for _ in 0..96 {
                let af = *first.step(&force, released).unwrap();
                let bf = *second.step(&force, released).unwrap();
                assert!(af.energy_residual_j.abs() <= af.energy_tolerance_j);
                assert!(bf.energy_residual_j.abs() <= bf.energy_tolerance_j);
                first
                    .spatial_cavity_pressures_into(&mut scratch, &mut a)
                    .unwrap();
                second
                    .spatial_cavity_pressures_into(&mut scratch, &mut b)
                    .unwrap();
                pressure_difference += (a[0] - b[0]).abs();
                sound_difference += (af.observer_pressure_pa - bf.observer_pressure_pa).abs();
            }
            assert!(pressure_difference > 0.0);
            assert!(sound_difference > 0.0);
            assert!(first.system().material_memory().iter().any(|z| *z != 0.0));
            let before = first.system().material_memory().to_vec();
            let tick = first.system().samples_rendered();
            let cancelled = CancelGate::new_clock_free();
            cancelled.request();
            assert!(first.step_under_gate(&force, released, &cancelled).is_err());
            assert_eq!(first.system().material_memory(), before);
            assert_eq!(first.system().samples_rendered(), tick);
        })
    });
}

#[test]
fn molecular_fe_admission_refuses_stale_source_and_incomplete_work() {
    with_cx(|cx, gate| {
        with_volume(2, gas(101325.0), |fe| {
            let gas = gas(101325.0);
            for (source, options) in [
                (
                    source(super::super::gas(202650.0), Some(0.025)),
                    molecular_options(),
                ),
                (
                    source(gas, Some(0.025)),
                    ShellCavityMolecularOptions {
                        maximum_caloric_coordinates: 1,
                        ..molecular_options()
                    },
                ),
                (
                    source(gas, Some(0.025)),
                    ShellCavityMolecularOptions {
                        maximum_memory_states: 1,
                        ..molecular_options()
                    },
                ),
            ] {
                assert!(
                    bank(cx, gate, fe, gas, 700.0)
                        .with_material_shell_cavity_thermal_molecular(
                            fe,
                            gas,
                            pressure_resolution(),
                            Vec::new(),
                            heat_resolution(),
                            source,
                            options,
                            gate
                        )
                        .is_err()
                );
            }
            let original = bank(cx, gate, fe, gas, 700.0);
            let cancelled = CancelGate::new_clock_free();
            cancelled.request();
            assert!(matches!(
                original.with_material_shell_cavity_thermal_molecular(
                    fe,
                    gas,
                    pressure_resolution(),
                    Vec::new(),
                    heat_resolution(),
                    source(gas, Some(0.025)),
                    molecular_options(),
                    &cancelled
                ),
                Err(ShellCavityError::Cancelled)
            ));
        })
    });
}
