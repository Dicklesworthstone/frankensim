//! G1/G3: two original source gases, one momentum, first and second laws.

use super::*;
use crate::cavity_phs::piston::{
    GasPiston,
    tests::{controls as single_controls, drive, gas},
};
use fs_ports::PortedModule;

fn controls() -> DoublePistonSolve {
    DoublePistonSolve {
        piston: single_controls(),
        opposing_gas: single_controls().gas,
        volume_tolerance_m3: 1e-18,
    }
}
/// For equal faces. The solve tries both ends of the primary bracket, and at
/// the far end of the one-gas bracket (1.5 litres) the opposing litre would
/// be halved in one step, which its gas owner refuses for this mixture
/// ("discrete gas compression exceeds positive-tangent branch"). Both ends
/// of this bracket are admitted by both gases.
fn equal_face_controls() -> DoublePistonSolve {
    DoublePistonSolve {
        piston: super::super::super::PistonSolve {
            volume_bracket_m3: [0.0008, 0.0012],
            ..single_controls()
        },
        ..controls()
    }
}
fn input(heat_j: [f64; 2]) -> DoublePistonInput {
    DoublePistonInput {
        duration_s: drive(0.0).duration_s,
        heat_j,
        external_force_n: 0.0,
    }
}
fn pair(velocity: f64, slopes: [f64; 2], opposing_area: f64) -> DoubleActingGasPiston {
    DoubleActingGasPiston::new(
        GasPiston::new(gas(0.75, slopes[0]), 1.0, 0.01, velocity).unwrap(),
        gas(0.25, slopes[1]),
        opposing_area,
    )
    .unwrap()
}
fn entropy(initial: &ClosedGasParcel, final_gas: &ClosedGasParcel) -> f64 {
    initial.mass_kg()
        * initial
            .model()
            .specific_entropy_change_j_kg_k(initial.point(), final_gas.point())
            .unwrap()
}
fn check(step: &DiscreteGradientDoublePistonStep) {
    assert!(step.momentum_residual_n_s.abs() <= controls().piston.momentum_tolerance_n_s);
    assert!(step.mechanical_energy_residual_j.abs() <= controls().piston.energy_tolerance_j);
    assert!(step.total_energy_residual_j.abs() <= controls().piston.energy_tolerance_j);
    assert!(step.volume_exchange_residual_m3.abs() <= controls().volume_tolerance_m3);
    assert!(step.material_trials <= 2 * (controls().piston.maximum_iterations + 3));
}

#[test]
fn adiabatic_motion_preserves_both_nonlinear_source_entropies() {
    for velocity in [-0.2, 0.2] {
        let mut piston = pair(velocity, [0.001, 0.002], 0.01);
        let initial = piston.clone();
        let mut work = 0.0;
        for i in 0..48 {
            let step = piston
                .step_discrete_gradient(
                    DoublePistonInput {
                        external_force_n: if i < 20 { -0.5 } else { 0.0 },
                        ..input([0.0, 0.0])
                    },
                    equal_face_controls(),
                    || false,
                )
                .unwrap();
            work += step.external_force_work_j;
            check(&step);
        }
        let mut stored = 0.0;
        for (old, new) in [
            (initial.primary().gas(), piston.primary().gas()),
            (initial.opposing(), piston.opposing()),
        ] {
            let tolerance =
                8.0 * 48.0 * old.mass_kg() * controls().piston.gas.energy_tolerance_j_kg
                    / controls().piston.gas.minimum_temperature_k;
            assert!(entropy(old, new).abs() <= tolerance);
            stored += old.mass_kg()
                * old
                    .model()
                    .thermodynamics()
                    .specific_internal_energy_change(
                        fs_qty::Temperature::new(old.state().state().temperature),
                        fs_qty::Temperature::new(new.state().state().temperature),
                    )
                    .unwrap()
                    .value();
            assert_eq!(old.mass_kg(), new.mass_kg());
            assert_eq!(new.point().axes().get("specimen-lot"), Some(&42.0));
        }
        let kinetic = 0.5
            * initial.primary().mass_kg()
            * (piston.primary().velocity_m_s().powi(2) - initial.primary().velocity_m_s().powi(2));
        assert!((stored + kinetic - work).abs() < 48.0 * controls().piston.energy_tolerance_j);
        assert!(
            (piston.primary().gas().volume_m3() + piston.opposing().volume_m3() - 0.002).abs()
                < 1e-16
        );
    }
}

#[test]
fn both_constant_capacity_chambers_match_their_own_analytic_isentropes() {
    for velocity in [-8.0, 8.0] {
        let mut piston = pair(velocity, [0.0; 2], 0.006);
        let initial = piston.clone();
        let step = piston
            .step_discrete_gradient(
                DoublePistonInput {
                    external_force_n: -400.0,
                    ..input([0.0, 0.0])
                },
                controls(),
                || false,
            )
            .unwrap();
        assert_ne!(step.displacement_m, 0.0);
        for (old, new) in [
            (initial.primary().gas(), piston.primary().gas()),
            (initial.opposing(), piston.opposing()),
        ] {
            let properties = old.state().thermodynamics().properties();
            let exponent = properties.mass_specific_gas_constant().value()
                / properties.mass_specific_cv().value();
            let exact = old.state().state().temperature
                * (old.volume_m3() / new.volume_m3()).powf(exponent);
            assert!((new.state().state().temperature / exact - 1.0).abs() < 3e-12);
        }
        check(&step);
    }
}

#[test]
fn independent_heat_ports_use_each_materials_entropy_conjugate_temperature() {
    let mut piston = pair(0.1, [0.001, 0.002], 0.006);
    let initial = piston.clone();
    let drive = DoublePistonInput {
        external_force_n: -400.0,
        ..input([20.0, -5.0])
    };
    let step = piston
        .step_discrete_gradient(drive, controls(), || false)
        .unwrap();
    for (i, old, new) in [
        (0, initial.primary().gas(), piston.primary().gas()),
        (1, initial.opposing(), piston.opposing()),
    ] {
        let transfer = drive.heat_j[i] / step.thermodynamic_mean[i].state().temperature;
        assert!((entropy(old, new) - transfer).abs() < 2e-11);
        assert_eq!(step.gas[i].entropy_transfer_j_k, transfer);
    }
    assert!(piston.audit_last_window().is_ok());
    assert!(piston.primary().audit_last_window().is_ok());
    check(&step);
}

#[test]
fn shared_impulse_tangent_differentiates_both_discrete_gibbs_pressure_laws() {
    let piston = pair(-0.1, [0.001, 0.002], 0.006);
    let drive = input([20.0, 7.0]);
    let volume = 0.00102;
    let delta = 1e-8;
    let center = piston
        .discrete_gradient_trial(volume, drive, controls(), &mut || false)
        .unwrap();
    let plus = piston
        .discrete_gradient_trial(volume + delta, drive, controls(), &mut || false)
        .unwrap();
    let minus = piston
        .discrete_gradient_trial(volume - delta, drive, controls(), &mut || false)
        .unwrap();
    let difference = (plus.residual - minus.residual) / (2.0 * delta);
    assert!((center.tangent / difference - 1.0).abs() < 2e-6);
}

#[test]
fn equilibrium_and_failed_proposals_preserve_both_gases_momentum_and_port_windows() {
    let mut piston = pair(0.0, [0.001, 0.002], 0.01);
    let initial = piston.clone();
    let stationary = piston
        .step_discrete_gradient(input([0.0; 2]), controls(), || false)
        .unwrap();
    assert_eq!(piston, initial);
    assert_eq!(stationary.total_energy_residual_j, 0.0);
    assert_eq!(stationary.material_trials, 2);
    let before = piston.clone();
    let transfers = piston.last_window_transfers().map(<[_]>::to_vec);
    let primary_transfers = piston.primary().last_window_transfers().map(<[_]>::to_vec);
    for drive in [
        input([0.0, 1e9]),
        DoublePistonInput {
            duration_s: 0.0,
            ..input([0.0; 2])
        },
    ] {
        assert!(
            piston
                .step_discrete_gradient(drive, controls(), || false)
                .is_err()
        );
        assert_eq!(piston, before);
    }
    let collapsed = DoublePistonSolve {
        piston: super::super::super::PistonSolve {
            volume_bracket_m3: [0.0007, 0.003],
            ..single_controls()
        },
        ..controls()
    };
    assert!(
        piston
            .step_discrete_gradient(input([0.0; 2]), collapsed, || false)
            .is_err()
    );
    let mut probe = before.clone();
    let mut polls = 0;
    probe
        .step_discrete_gradient(input([20.0, 5.0]), equal_face_controls(), || {
            polls += 1;
            false
        })
        .unwrap();
    for cancel_at in [1, polls / 2, polls] {
        let mut current = 0;
        assert!(matches!(
            piston.step_discrete_gradient(input([20.0, 5.0]), equal_face_controls(), || {
                current += 1;
                current == cancel_at
            }),
            Err(PistonError::Cancelled)
        ));
        assert_eq!(piston, before);
        assert_eq!(piston.last_window_transfers().map(<[_]>::to_vec), transfers);
        assert_eq!(
            piston.primary().last_window_transfers().map(<[_]>::to_vec),
            primary_transfers
        );
    }
}

mod real_gas;
