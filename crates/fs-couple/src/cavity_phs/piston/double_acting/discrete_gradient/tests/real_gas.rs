//! Two real-gas stores keep their original entropy and share one wall effort.

use super::*;
use crate::cavity_phs::piston::tests::real_gas::{real_controls, real_parcel};
use fs_material::thermochemical_gas::nonideal::energy::parcel::IsochoricRealGasParcel;

fn config() -> DoublePistonSolve {
    DoublePistonSolve {
        piston: real_controls(),
        opposing_gas: real_controls().gas,
        volume_tolerance_m3: 1e-18,
    }
}
fn drive(heat_j: [f64; 2]) -> DoublePistonInput {
    DoublePistonInput {
        duration_s: 0.00002,
        heat_j,
        external_force_n: -100.0,
    }
}
fn pair(velocity: f64) -> DoubleActingGasPiston<IsochoricRealGasParcel> {
    DoubleActingGasPiston::new(
        GasPiston::new(real_parcel(300.0, 5e6), 2.0, 0.01, velocity).unwrap(),
        real_parcel(340.0, 4.5e6),
        0.008,
    )
    .unwrap()
}

#[test]
fn both_dense_chambers_preserve_entropy_while_one_wall_exchanges_their_work() {
    for velocity in [-0.4, 0.4] {
        let mut piston = pair(velocity);
        let initial = piston.clone();
        let mut external_work = 0.0;
        for _ in 0..16 {
            let report = piston
                .step_discrete_gradient(drive([0.0; 2]), config(), || false)
                .unwrap();
            external_work += report.external_force_work_j;
            assert!(report.momentum_residual_n_s.abs() <= config().piston.momentum_tolerance_n_s);
            assert!(report.total_energy_residual_j.abs() <= config().piston.energy_tolerance_j);
            assert!(report.volume_exchange_residual_m3.abs() <= config().volume_tolerance_m3);
            let effort = piston.primary().area_m2() * report.gas[0].work_pressure_pa
                - piston.opposing_area_m2() * report.gas[1].work_pressure_pa
                + drive([0.0; 2]).external_force_n;
            assert!(
                (piston.primary().mass_kg() * (report.velocity_m_s[1] - report.velocity_m_s[0])
                    - drive([0.0; 2]).duration_s * effort)
                    .abs()
                    <= config().piston.momentum_tolerance_n_s
            );
        }
        let mut energy = 0.0;
        for (old, new) in [
            (initial.primary().gas(), piston.primary().gas()),
            (initial.opposing(), piston.opposing()),
        ] {
            let entropy = old.mass_kg()
                * old
                    .model()
                    .specific_entropy_change_between_densities(
                        old.density_kg_m3(),
                        old.state().temperature_k(),
                        new.density_kg_m3(),
                        new.state().temperature_k(),
                    )
                    .unwrap();
            let bound = 8.0 * 16.0 * old.mass_kg() * real_controls().gas.energy_tolerance_j_kg
                / real_controls().gas.minimum_temperature_k;
            assert!(entropy.abs() <= bound);
            energy += old.mass_kg()
                * old
                    .model()
                    .specific_internal_energy_change_between_densities(
                        old.density_kg_m3(),
                        old.state().temperature_k(),
                        new.density_kg_m3(),
                        new.state().temperature_k(),
                    )
                    .unwrap();
            assert_eq!(old.mass_kg(), new.mass_kg());
        }
        let kinetic = 0.5
            * piston.primary().mass_kg()
            * (piston.primary().velocity_m_s().powi(2) - initial.primary().velocity_m_s().powi(2));
        assert!(
            (energy + kinetic - external_work).abs() < 16.0 * config().piston.energy_tolerance_j
        );
        assert!(piston.audit_last_window().is_ok());
    }
}

#[test]
fn nonideal_two_face_tangent_and_independent_heat_entropy_are_consistent() {
    let piston = pair(0.2);
    let input = drive([100.0, -30.0]);
    let volume = 0.001001;
    let h = 1e-8;
    let trial = piston
        .discrete_gradient_trial(volume, input, config(), &mut || false)
        .unwrap();
    let plus = piston
        .discrete_gradient_trial(volume + h, input, config(), &mut || false)
        .unwrap();
    let minus = piston
        .discrete_gradient_trial(volume - h, input, config(), &mut || false)
        .unwrap();
    assert!((trial.tangent / ((plus.residual - minus.residual) / (2.0 * h)) - 1.0).abs() < 2e-6);
    let mut next = piston.clone();
    let report = next
        .step_discrete_gradient(input, config(), || false)
        .unwrap();
    for (i, old, new) in [
        (0, piston.primary().gas(), next.primary().gas()),
        (1, piston.opposing(), next.opposing()),
    ] {
        let entropy = old.mass_kg()
            * old
                .model()
                .specific_entropy_change_between_densities(
                    old.density_kg_m3(),
                    old.state().temperature_k(),
                    new.density_kg_m3(),
                    new.state().temperature_k(),
                )
                .unwrap();
        assert!(
            (entropy - input.heat_j[i] / report.thermodynamic_mean[i].temperature_k()).abs() < 1e-9
        );
    }
    assert!(report.total_energy_residual_j.abs() <= config().piston.energy_tolerance_j);
    assert!(next.audit_last_window().is_ok());
}
