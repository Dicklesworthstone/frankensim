//! Finite phase tanks exchange integrated H while both inventories store U.
#[path = "common/nonideal_flash.rs"]
#[allow(dead_code)]
mod common;

use common::{METHANE, MIXTURE_MASS, model, source};
use fs_material::thermochemical_gas::nonideal::{
    PengRobinsonSpeciesBinding, ThermochemicalPengRobinsonModel as Model,
    flash::{
        equilibrium::{
            PhaseEnthalpyControls,
            entropy::{PhaseEntropyControls, volume::PhaseEntropyVolumeControls},
            volume::PhaseEnergyVolumeControls,
        },
        parcel::{
            IsochoricMixtureEquilibriumParcel as Parcel,
            discharge::MixtureEquilibriumDischargeControl,
            transfer::{
                MixtureEquilibriumInflowControl,
                finite::{
                    MixtureEquilibriumTransferControl as Control,
                    MixtureEquilibriumTransferError as E,
                    transfer_adiabatic_mixture_equilibrium as transfer,
                },
            },
        },
    },
};
use fs_qty::{Pressure, Temperature};
use fs_thermochem::{Composition, CompositionBasis, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}; difference {}, allowance {tolerance}",
        actual - expected
    );
}
fn controls(dp: [f64; 2], dt: [f64; 2], rp: [f64; 2], rt: [f64; 2]) -> Control {
    Control::new(
        MixtureEquilibriumDischargeControl::new(PhaseEntropyVolumeControls::new(
            dp.map(Pressure::new),
            PhaseEntropyControls::new(dt.map(Temperature::new)),
        )),
        MixtureEquilibriumInflowControl::new(PhaseEnergyVolumeControls::new(
            rp.map(Pressure::new),
            PhaseEnthalpyControls::new(rt.map(Temperature::new)),
        )),
    )
}
fn dense_controls() -> Control {
    controls(
        [1.8e6, 3.2e6],
        [195.0, 235.0],
        [0.7e6, 1.4e6],
        [195.0, 270.0],
    )
}
fn dilute_controls() -> Control {
    let mut c = controls([6.0, 10.5], [290.0, 360.0], [1.1, 1.5], [290.0, 380.0]);
    c.donor.equilibrium.specific_volume_tolerance_m3_kg = 1e-6;
    c.donor.equilibrium.pressure_tolerance_pa = 1e-10;
    c.receiver.equilibrium.specific_volume_tolerance_m3_kg = 1e-6;
    c.receiver.equilibrium.pressure_tolerance_pa = 1e-10;
    c
}
fn tank(model: Model, mass: f64, temperature: f64, pressure: f64) -> Parcel {
    Parcel::from_pressure_temperature(
        model,
        mass,
        Temperature::new(temperature),
        Pressure::new(pressure),
        dense_controls().receiver.equilibrium.phase.phase_options(),
        || false,
    )
    .unwrap()
}
fn dense_pair(model: &Model) -> [Parcel; 2] {
    [
        tank(model.clone(), 1.0, 220.0, 3e6),
        tank(model.clone(), 10.0, 230.0, 1e6),
    ]
}
fn composition(base: &Model, methane: f64) -> Model {
    let source = base
        .source()
        .with_composition(
            Composition::new(CompositionBasis::MoleFraction, vec![1.0 - methane, methane]).unwrap(),
        )
        .unwrap();
    let bindings = base
        .source()
        .thermodynamics()
        .components()
        .iter()
        .zip(base.mechanics().components())
        .map(|(s, c)| PengRobinsonSpeciesBinding {
            species: s.species().clone(),
            critical_temperature: c.critical_temperature,
            critical_pressure: c.critical_pressure,
            acentric_factor: c.acentric_factor,
        })
        .collect();
    Model::new(
        source,
        bindings,
        base.mechanics().binary_interactions().to_vec(),
    )
    .unwrap()
}

#[test]
fn two_phase_donor_exports_integrated_h_and_both_finite_tanks_conserve_u() {
    let model = source();
    let [mut donor, mut receiver] = dense_pair(&model);
    let before = [donor.clone(), receiver.clone()];
    // Independent 50-digit PR cubic/fugacity/entropy endpoint already used by
    // the native donor tests: final rho=78.02103439554112 kg/m3 at 2 MPa.
    let requested = donor.mass_kg() - donor.volume_m3() * 78.021_034_395_541_12;
    let c = dense_controls();
    let step = transfer(&mut donor, &mut receiver, requested, None, c, || false).unwrap();
    close(donor.state().temperature_k(), 208.202_765_109_369_35, 2e-6);
    close(donor.state().pressure_pa(), 2e6, 0.2);
    close(
        step.discharge.removed_mass_kg,
        0.403_719_163_379_537_4,
        2e-9,
    );
    close(step.transported_enthalpy_j, -8277.049_580_656_221, 3e-4);
    close(
        step.discharge.exported_enthalpy_adjustment_j,
        -1748.002_274_066_815_5,
        2e-4,
    );
    assert!(donor.state().phase().two_phase().is_some());
    assert!(
        donor.state().vapor_mass_fraction().unwrap()
            > before[0].state().vapor_mass_fraction().unwrap()
    );
    assert_eq!(donor.volume_m3(), before[0].volume_m3());
    assert_eq!(receiver.volume_m3(), before[1].volume_m3());
    assert_eq!(donor.model(), before[0].model());
    assert_eq!(
        step.transported_enthalpy_j,
        step.discharge.exported_enthalpy_j
    );
    close(
        receiver.mass_kg() - before[1].mass_kg(),
        step.discharge.removed_mass_kg,
        2e-15,
    );
    // This independent moderate-datum audit subtracts actual state U, never H.
    let changes = [
        donor.internal_energy_j() - before[0].internal_energy_j(),
        receiver.internal_energy_j() - before[1].internal_energy_j(),
    ];
    close(changes[0] + step.transported_enthalpy_j, 0.0, 2e-5);
    close(
        changes[1] - step.transported_enthalpy_j,
        0.0,
        c.energy_tolerance_j,
    );
    close(
        changes[0] + changes[1],
        step.combined_energy_residual_j,
        c.energy_tolerance_j,
    );
    close(
        step.receiver_internal_energy_change_j,
        changes[1],
        c.energy_tolerance_j,
    );
    let entropy = donor.mass_kg() * donor.state().specific_entropy_j_kg_k()
        - before[0].mass_kg() * before[0].state().specific_entropy_j_kg_k()
        + receiver.mass_kg() * receiver.state().specific_entropy_j_kg_k()
        - before[1].mass_kg() * before[1].state().specific_entropy_j_kg_k();
    close(step.entropy_production_j_k, entropy, 2e-8);
    assert!(step.entropy_production_j_k > 0.0);
    assert!(step.mass_residual_kg.abs() < 2e-14);
    assert!(step.species.iter().all(|s| s.residual_kg.abs() < 2e-14));
    assert!(step.combined_energy_residual_j.abs() <= c.energy_tolerance_j);
    assert!(donor.state().pressure_pa() > receiver.state().pressure_pa());
    assert!(
        step.pt_evaluations
            <= c.donor.equilibrium.maximum_pt_evaluations
                + c.receiver.equilibrium.maximum_pt_evaluations
    );
    assert!(
        step.eos_evaluation_upper_bound
            <= c.donor.equilibrium.maximum_eos_evaluation_bound().unwrap()
                + c.receiver
                    .equilibrium
                    .maximum_eos_evaluation_bound()
                    .unwrap()
    );
}

#[test]
fn dilute_finite_transfer_matches_integrated_blowdown_and_distinct_species_filling() {
    let base = source();
    for methane in [0.5, 0.8] {
        let mut donor = tank(composition(&base, methane), 1.0, 350.0, 10.0);
        let mut receiver = tank(base.clone(), 1.0, 300.0, 1.0);
        let volume = receiver.volume_m3();
        let step = transfer(
            &mut donor,
            &mut receiver,
            0.2,
            None,
            dilute_controls(),
            || false,
        )
        .unwrap();
        let ratio = 0.8_f64;
        let donor_mass = methane * 0.016_042_46 + (1.0 - methane) * 0.030_069_04;
        let cp = (3.5 * methane + 4.5 * (1.0 - methane)) * R / donor_mass;
        let cv = cp - R / donor_mass;
        let gamma = cp / cv;
        let exported = cv * 350.0 * (1.0 - ratio.powf(gamma));
        let receiver_cv = 3.0 * R / MIXTURE_MASS;
        let temperature = (receiver_cv * 300.0 + exported) / (receiver_cv + 0.2 * cv);
        let frozen_temperature =
            (receiver_cv * 300.0 + 0.2 * cp * 350.0) / (receiver_cv + 0.2 * cv);
        let n0 = 1.0 / MIXTURE_MASS;
        let incoming = 0.2 / donor_mass;
        close(
            donor.state().temperature_k(),
            350.0 * ratio.powf(gamma - 1.0),
            8e-4,
        );
        close(donor.state().pressure_pa(), 10.0 * ratio.powf(gamma), 2e-5);
        close(step.transported_enthalpy_j, exported, 0.5);
        close(receiver.state().temperature_k(), temperature, 8e-4);
        close(
            receiver.state().pressure_pa(),
            (n0 + incoming) * R * temperature / volume,
            5e-6,
        );
        close(
            receiver.model().mechanics().mole_fractions()[METHANE],
            (0.5 * n0 + methane * incoming) / (n0 + incoming),
            3e-15,
        );
        assert!(
            frozen_temperature - receiver.state().temperature_k() > 1.0,
            "freezing initial donor h over the finite portion overheats the receiving tank"
        );
        assert!(step.discharge.exported_enthalpy_adjustment_j < -1000.0);
        assert!(step.species.iter().all(|s| s.residual_kg.abs() < 1e-14));
        assert!(step.entropy_production_j_k > 0.0);
    }
}

#[test]
fn subdividing_the_same_finite_bulk_portion_preserves_both_endpoints() {
    let base = source();
    let initial = [
        tank(base.clone(), 1.0, 350.0, 10.0),
        tank(base, 1.0, 300.0, 1.0),
    ];
    let [mut whole_donor, mut whole_receiver] = initial.clone();
    let [mut split_donor, mut split_receiver] = initial;
    let c = dilute_controls();
    let whole = transfer(&mut whole_donor, &mut whole_receiver, 0.2, None, c, || {
        false
    })
    .unwrap();
    let first = transfer(&mut split_donor, &mut split_receiver, 0.1, None, c, || {
        false
    })
    .unwrap();
    let second = transfer(&mut split_donor, &mut split_receiver, 0.1, None, c, || {
        false
    })
    .unwrap();
    for (a, b) in [
        (&whole_donor, &split_donor),
        (&whole_receiver, &split_receiver),
    ] {
        close(a.mass_kg(), b.mass_kg(), 4e-16);
        assert_eq!(a.volume_m3(), b.volume_m3());
        close(a.state().temperature_k(), b.state().temperature_k(), 2e-6);
        close(a.state().pressure_pa(), b.state().pressure_pa(), 2e-7);
    }
    close(
        whole.transported_enthalpy_j,
        first.transported_enthalpy_j + second.transported_enthalpy_j,
        5e-4,
    );
    close(
        whole.entropy_production_j_k,
        first.entropy_production_j_k + second.entropy_production_j_k,
        2e-6,
    );
}

#[test]
fn enormous_original_species_datums_leave_paired_physics_and_residuals_unchanged() {
    let [mut a_donor, mut a_receiver] = dense_pair(&source());
    let [mut b_donor, mut b_receiver] =
        dense_pair(&model([1e22, -3e21], [2e18, -7e17], [100.0, 450.0]));
    let amount = a_donor.mass_kg() - a_donor.volume_m3() * 78.021_034_395_541_12;
    let a = transfer(
        &mut a_donor,
        &mut a_receiver,
        amount,
        None,
        dense_controls(),
        || false,
    )
    .unwrap();
    let b = transfer(
        &mut b_donor,
        &mut b_receiver,
        amount,
        None,
        dense_controls(),
        || false,
    )
    .unwrap();
    for (plain, shifted) in [(&a_donor, &b_donor), (&a_receiver, &b_receiver)] {
        assert_eq!(
            plain.state().temperature_k(),
            shifted.state().temperature_k()
        );
        assert_eq!(plain.state().pressure_pa(), shifted.state().pressure_pa());
        assert_eq!(plain.state().phase(), shifted.state().phase());
        assert_eq!(plain.mass_kg(), shifted.mass_kg());
        assert_eq!(plain.volume_m3(), shifted.volume_m3());
    }
    assert_eq!(
        a.discharge.exported_enthalpy_adjustment_j,
        b.discharge.exported_enthalpy_adjustment_j
    );
    assert_eq!(a.donor_energy_residual_j, b.donor_energy_residual_j);
    assert_eq!(a.receiver_energy_residual_j, b.receiver_energy_residual_j);
    assert_eq!(a.combined_energy_residual_j, b.combined_energy_residual_j);
    assert_eq!(a.entropy_production_j_k, b.entropy_production_j_k);
    assert_eq!(a.species, b.species);
    assert_eq!(a.pt_evaluations, b.pt_evaluations);
    assert!(b.transported_enthalpy_j.abs() > 1e19);
}

#[test]
fn idle_late_pressure_refusal_source_mismatch_and_cancellation_are_atomic() {
    let base = source();
    let original = [
        tank(base.clone(), 1.0, 350.0, 10.0),
        tank(base.clone(), 1.0, 300.0, 1.0),
    ];
    let [mut donor, mut receiver] = original.clone();
    let c = dilute_controls();
    let idle = transfer(&mut donor, &mut receiver, 0.0, None, c, || false).unwrap();
    assert_eq!([donor.clone(), receiver.clone()], original);
    assert_eq!(idle.pt_evaluations, 0);
    assert_eq!(idle.eos_evaluation_upper_bound, 0);
    assert_eq!(idle.transported_enthalpy_j, 0.0);
    let mut polls = 0;
    transfer(&mut donor, &mut receiver, 0.2, None, c, || {
        polls += 1;
        false
    })
    .unwrap();
    for stop in [1, polls] {
        let [mut donor, mut receiver] = original.clone();
        let mut calls = 0;
        assert!(matches!(
            transfer(&mut donor, &mut receiver, 0.2, None, c, || {
                calls += 1;
                calls == stop
            }),
            Err(E::Cancelled)
        ));
        assert_eq!([donor, receiver], original);
    }
    let [mut donor, _] = original.clone();
    let mut receiver = tank(base, 1.0, 300.0, 5.0);
    let before = [donor.clone(), receiver.clone()];
    let mut crossing = c;
    crossing.donor.equilibrium.pressure_bracket = [Pressure::new(5.0), Pressure::new(10.5)];
    crossing.donor.equilibrium.phase.temperature_bracket =
        [Temperature::new(280.0), Temperature::new(360.0)];
    crossing.receiver.equilibrium.pressure_bracket = [Pressure::new(6.0), Pressure::new(9.0)];
    assert!(matches!(
        transfer(&mut donor, &mut receiver, 0.3, None, crossing, || false),
        Err(E::PressureHead { .. })
    ));
    assert_eq!([donor, receiver], before);
    let [mut donor, _] = original.clone();
    let mut receiver = tank(model([1.0, 0.0], [0.0; 2], [100.0, 450.0]), 1.0, 300.0, 1.0);
    let before = [donor.clone(), receiver.clone()];
    assert!(matches!(
        transfer(&mut donor, &mut receiver, 0.2, None, c, || false),
        Err(E::Receiver(_))
    ));
    assert_eq!([donor, receiver], before);
}
