//! Prescribed thermal reference changes retain source datum and viscous memory.
mod plane;

use super::*;
use fs_couple::finite_solid::thermal::visco::{
    ExpansionViscoelastic, ExpansionViscoelasticError, ThermalReferenceOrder,
};
use fs_material::hyper::reference::ReferenceDeformation;
use fs_material::hyper::{Hyperelastic, HyperelasticModel};
use fs_material::visco::temperature::{TemperatureShift, TemperatureShiftModel};
use fs_material::visco_hyper::{BranchEnergy, FiniteStrainBranch, ViscoHyperelastic};

const TOTAL: [f64; 9] = [1.18, 0.06, 0.0, 0.0, 0.94, 0.01, 0.02, 0.0, 1.03];

fn original_law() -> ViscoHyperelastic {
    ViscoHyperelastic::new(
        Hyperelastic::new(
            HyperelasticModel::MooneyRivlin {
                c10: 6.0,
                c01: 2.0,
                kappa: 80.0,
            },
            5.0,
        )
        .unwrap(),
        vec![FiniteStrainBranch {
            energy: BranchEnergy::NeoHookean { mu_pa: 8.0 },
            tau_s: 0.3,
        }],
    )
    .unwrap()
}

fn shift() -> TemperatureShift {
    TemperatureShift::new(
        TemperatureShiftModel::Arrhenius {
            activation_energy_j_mol: 18_000.0,
        },
        300.0,
        [290.0, 310.0],
    )
    .unwrap()
}

fn endpoint(source: &MaterialCard, kinetics: Option<TemperatureShift>) -> ExpansionViscoelastic {
    let thermal = ExpansionThermalReference::try_new(
        source,
        path(source, false, "T", point(300.0)),
        300.0,
        frame(Q),
    )
    .unwrap();
    ExpansionViscoelastic::try_new(original_law(), thermal, kinetics).unwrap()
}

#[test]
fn changing_temperature_retains_memory_and_separates_reference_work_from_heat() {
    let source = card(STRETCHES, "T", false);
    let initial = endpoint(&source, Some(shift()));
    let placed = initial
        .law()
        .step_response(&initial.law().initial_state(), &TOTAL, 0.0)
        .unwrap();
    let history = initial
        .law()
        .step_response(&placed.state, &TOTAL, 0.12)
        .unwrap()
        .state;
    let rem_placed = initial
        .law()
        .step_remainder_response(&initial.law().initial_remainder_state(), &TOTAL, 0.0)
        .unwrap();
    let remainder = initial
        .law()
        .step_remainder_response(&rem_placed.state, &TOTAL, 0.12)
        .unwrap()
        .state;
    assert!(history.dissipated > 0.0);
    assert_ne!(
        history.inverse_viscous_cauchy_green[0],
        [1.0, 1.0, 1.0, 0.0, 0.0, 0.0]
    );
    let original = history.clone();
    let hot = initial.at_temperature(&source, 310.0).unwrap();
    let full = hot.change_reference_from(&initial, &history).unwrap();
    let rem = hot
        .change_remainder_reference_from(&initial, &remainder)
        .unwrap();
    assert_eq!(history, original);
    assert_eq!(full.state.deformation_gradient, TOTAL);
    assert_eq!(
        full.state.inverse_viscous_cauchy_green,
        history.inverse_viscous_cauchy_green
    );
    assert_eq!(full.state.work, history.work);
    assert_eq!(full.state.dissipated, history.dissipated);
    assert_eq!(rem.state.work, remainder.work);
    assert_eq!(rem.state.dissipated, remainder.dissipated);
    near(
        full.reference_work_pa,
        full.state.stored - history.stored,
        2e-12,
    );
    near(
        rem.reference_work_pa,
        rem.state.stored - remainder.stored,
        2e-12,
    );
    near(
        full.state.work + full.state.reference_work - full.state.stored - full.state.dissipated,
        0.0,
        3e-11,
    );
    near(
        rem.state.work + rem.state.reference_work - rem.state.stored - rem.state.dissipated,
        0.0,
        3e-11,
    );
    // Reversing the prescribed temperature at unchanged F and Cv is an
    // energy exchange, not a second irreversible relaxation event.
    let returned = hot.at_temperature(&source, 300.0).unwrap();
    let back = returned.change_reference_from(&hot, &full.state).unwrap();
    near(back.reference_work_pa, -full.reference_work_pa, 3e-11);
    near(back.state.stored, history.stored, 3e-11);
    near(back.state.reference_work, 0.0, 3e-11);
    assert_eq!(back.state.dissipated, history.dissipated);
    assert_eq!(
        back.state.inverse_viscous_cauchy_green,
        history.inverse_viscous_cauchy_green
    );
    // Physical loss occurs only when the retained branches subsequently age.
    let relaxed = hot.law().step_response(&full.state, &TOTAL, 0.08).unwrap();
    assert!(relaxed.dissipated_pa > 0.0);
    assert_eq!(relaxed.state.reference_work, full.state.reference_work);
    near(relaxed.loading_work_pa, 0.0, 3e-11);
}

#[test]
fn thermal_endpoints_always_use_original_length_datum_and_kinetic_calibration() {
    let source = card(STRETCHES, "T", false);
    let original = original_law();
    let calibration = shift();
    let mut current = endpoint(&source, Some(calibration));
    let source_identity = current.thermal_reference().source().identity();
    for temperature in [305.0, 310.0, 296.0, 307.0, 300.0] {
        current = current.at_temperature(&source, temperature).unwrap();
        assert_eq!(current.temperature_k(), temperature);
        assert_eq!(
            current.thermal_reference().source().identity(),
            source_identity
        );
        assert_eq!(current.thermal_reference().frame(), frame(Q));
        let fraction = (temperature - 300.0) / 10.0;
        for (actual, endpoint_stretch) in current
            .thermal_reference()
            .geometry()
            .principal_stretches()
            .into_iter()
            .zip(STRETCHES)
        {
            near(actual, 1.0 + fraction * (endpoint_stretch - 1.0), 2e-14);
        }
        let expected = calibration
            .shift_hyperelastic(&original, temperature)
            .unwrap();
        assert_eq!(current.law().branches, expected.law().branches);
    }
    assert_eq!(
        current.law().reference_deformation(),
        &ReferenceDeformation::identity()
    );
    assert_eq!(current.law().branches, original.branches);
    let fixed = endpoint(&source, None)
        .at_temperature(&source, 310.0)
        .unwrap();
    assert_eq!(fixed.law().branches, original.branches);
}

#[test]
fn thermal_history_refuses_changed_source_frame_datum_and_unsupported_temperature() {
    let source = card(STRETCHES, "T", false);
    let initial = endpoint(&source, Some(shift()));
    let state = initial.law().initial_state();
    let before = state.clone();
    assert!(initial.at_temperature(&source, 310.01).is_err());
    let other = card([1.11, 0.96, 1.04], "T", false);
    assert!(initial.at_temperature(&other, 305.0).is_err());
    let unrotated = ExpansionThermalReference::try_new(
        &source,
        path(&source, false, "T", point(300.0)),
        305.0,
        frame([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
    )
    .unwrap();
    let changed_frame =
        ExpansionViscoelastic::try_new(original_law(), unrotated, Some(shift())).unwrap();
    assert!(matches!(
        changed_frame.change_reference_from(&initial, &state),
        Err(ExpansionViscoelasticError::Incompatible(_))
    ));
    let changed_datum = ExpansionThermalReference::try_new(
        &source,
        path(&source, false, "T", point(295.0)),
        305.0,
        frame(Q),
    )
    .unwrap();
    let changed_datum =
        ExpansionViscoelastic::try_new(original_law(), changed_datum, Some(shift())).unwrap();
    assert!(matches!(
        changed_datum.change_reference_from(&initial, &state),
        Err(ExpansionViscoelasticError::Incompatible(_))
    ));
    let tight_shift = TemperatureShift::new(
        TemperatureShiftModel::Arrhenius {
            activation_energy_j_mol: 18_000.0,
        },
        300.0,
        [298.0, 308.0],
    )
    .unwrap();
    assert!(matches!(
        endpoint(&source, Some(tight_shift)).at_temperature(&source, 309.0),
        Err(ExpansionViscoelasticError::Kinetics(_))
    ));
    let prestrained = original_law().with_reference_deformation(
        ReferenceDeformation::try_new([1.1, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]).unwrap(),
    );
    assert!(
        ExpansionViscoelastic::try_new(prestrained, initial.thermal_reference().clone(), None)
            .is_err()
    );
    assert_eq!(state, before);
}

fn manufactured() -> ReferenceDeformation {
    ReferenceDeformation::try_new([1.1, 0.17, 0.03, 0.0, 0.94, -0.04, 0.0, 0.0, 1.02]).unwrap()
}

#[test]
fn manufactured_prestrain_survives_thermal_cycles_in_both_explicit_orders() {
    let source = card(STRETCHES, "T", false);
    for order in [
        ThermalReferenceOrder::ManufacturedThenThermal,
        ThermalReferenceOrder::ThermalThenManufactured,
    ] {
        let original = original_law().with_reference_deformation(manufactured());
        let thermal = endpoint(&source, None).thermal_reference().clone();
        let mut current = ExpansionViscoelastic::try_new_with_prestrain(
            original.clone(),
            thermal,
            Some(shift()),
            order,
        )
        .unwrap();
        assert_eq!(current.law().reference_deformation(), &manufactured());
        let placed = current
            .law()
            .step_response(&current.law().initial_state(), &TOTAL, 0.0)
            .unwrap();
        let mut history = current
            .law()
            .step_response(&placed.state, &TOTAL, 0.12)
            .unwrap()
            .state;
        let initial = history.clone();
        assert!(initial.dissipated > 0.0);
        for temperature in [305.0, 310.0, 296.0, 307.0, 300.0] {
            let next = current.at_temperature(&source, temperature).unwrap();
            let thermal = *next.thermal_reference().reference_deformation();
            let (left, right) = match order {
                ThermalReferenceOrder::ManufacturedThenThermal => (thermal, manufactured()),
                ThermalReferenceOrder::ThermalThenManufactured => (manufactured(), thermal),
            };
            for i in 0..3 {
                for j in 0..3 {
                    let expected: f64 = (0..3)
                        .map(|k| {
                            left.deformation_gradient()[3 * i + k]
                                * right.deformation_gradient()[3 * k + j]
                        })
                        .sum();
                    near(
                        next.law().reference_deformation().deformation_gradient()[3 * i + j],
                        expected,
                        3e-14,
                    );
                }
            }
            assert_eq!(next.manufactured_reference(), &manufactured());
            assert_eq!(next.reference_order(), order);
            let shifted = shift().shift_hyperelastic(&original, temperature).unwrap();
            assert_eq!(next.law().branches, shifted.law().branches);
            let change = next.change_reference_from(&current, &history).unwrap();
            near(
                change.reference_work_pa,
                change.state.stored - history.stored,
                3e-12,
            );
            assert_eq!(
                change.state.inverse_viscous_cauchy_green,
                initial.inverse_viscous_cauchy_green
            );
            assert_eq!(change.state.work, initial.work);
            assert_eq!(change.state.dissipated, initial.dissipated);
            history = change.state;
            current = next;
        }
        assert_eq!(current.law().reference_deformation(), &manufactured());
        near(history.stored, initial.stored, 3e-11);
        near(history.reference_work, initial.reference_work, 3e-11);
    }
}

#[test]
fn thermal_history_retains_manufacturing_order_even_at_a_commuting_endpoint() {
    let source = card(STRETCHES, "T", false);
    let thermal = endpoint(&source, None).thermal_reference().clone();
    let original = original_law().with_reference_deformation(manufactured());
    let after = ExpansionViscoelastic::try_new_with_prestrain(
        original.clone(),
        thermal.clone(),
        None,
        ThermalReferenceOrder::ManufacturedThenThermal,
    )
    .unwrap();
    let before = ExpansionViscoelastic::try_new_with_prestrain(
        original,
        thermal.clone(),
        None,
        ThermalReferenceOrder::ThermalThenManufactured,
    )
    .unwrap();
    let accepted = after.law().initial_state();
    let preserved = accepted.clone();
    // At the datum Gt = I, both endpoints have identical G. A change in the
    // declared physical composition still changes the future thermal model.
    assert_eq!(
        after.law().reference_deformation(),
        before.law().reference_deformation()
    );
    assert!(matches!(
        before.change_reference_from(&after, &accepted),
        Err(ExpansionViscoelasticError::Incompatible(_))
    ));
    let after_hot = after.at_temperature(&source, 310.0).unwrap();
    let before_hot = before.at_temperature(&source, 310.0).unwrap();
    assert_ne!(
        after_hot.law().reference_deformation(),
        before_hot.law().reference_deformation()
    );
    let changed_manufacture = ExpansionViscoelastic::try_new_with_prestrain(
        original_law(),
        thermal,
        None,
        ThermalReferenceOrder::ManufacturedThenThermal,
    )
    .unwrap();
    assert!(matches!(
        changed_manufacture.change_reference_from(&after, &accepted),
        Err(ExpansionViscoelasticError::Incompatible(_))
    ));
    assert_eq!(accepted, preserved);
}
