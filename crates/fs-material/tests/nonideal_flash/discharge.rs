//! Original equilibrium blowdown: integrated outgoing H, fixed V and changing mass.
use super::{close, common};
use common::{MIXTURE_MASS, model, source};
use fs_material::thermochemical_gas::nonideal::{
    ThermochemicalPengRobinsonModel,
    flash::{
        ThermochemicalFlashError as E,
        equilibrium::entropy::{PhaseEntropyControls, volume::PhaseEntropyVolumeControls},
        parcel::{
            IsochoricMixtureEquilibriumParcel, discharge::MixtureEquilibriumDischargeControl,
        },
    },
};
use fs_qty::{Pressure, Temperature};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

fn controls() -> MixtureEquilibriumDischargeControl {
    MixtureEquilibriumDischargeControl::new(PhaseEntropyVolumeControls::new(
        [Pressure::new(1.8e6), Pressure::new(3.2e6)],
        PhaseEntropyControls::new([Temperature::new(195.0), Temperature::new(235.0)]),
    ))
}

fn parcel(source: ThermochemicalPengRobinsonModel) -> IsochoricMixtureEquilibriumParcel {
    IsochoricMixtureEquilibriumParcel::from_pressure_temperature(
        source,
        1.0,
        Temperature::new(220.0),
        Pressure::new(3e6),
        controls().equilibrium.phase.phase_options(),
        || false,
    )
    .unwrap()
}

fn portion_to_independent_endpoint(vessel: &IsochoricMixtureEquilibriumParcel) -> f64 {
    vessel.mass_kg() - vessel.volume_m3() * 78.021_034_395_541_12
}

#[test]
fn g1_discharge_matches_independent_phase_endpoint_and_integrated_enthalpy() {
    let mut vessel = parcel(source());
    let original = vessel.clone();
    let portion = portion_to_independent_endpoint(&original);
    let report = vessel
        .discharge_adiabatic(portion, Pressure::new(1.5e6), controls(), || false)
        .unwrap();
    // Independent cubic/fugacity/entropy oracle from the original PS/SV tests.
    // Prescribe only outgoing mass and fixed geometry; recover final T and p.
    close(vessel.state().temperature_k(), 208.202_765_109_369_35, 2e-6);
    close(vessel.state().pressure_pa(), 2e6, 0.2);
    close(vessel.state().density_kg_m3(), 78.021_034_395_541_12, 2e-6);
    close(vessel.mass_kg(), 0.596_280_836_620_462_6, 2e-9);
    close(report.removed_mass_kg, 0.403_719_163_379_537_4, 2e-9);
    assert_eq!(vessel.volume_m3().to_bits(), original.volume_m3().to_bits());
    assert_eq!(vessel.model(), original.model());
    assert!(vessel.state().phase().two_phase().is_some());
    assert!(
        vessel.state().vapor_mass_fraction().unwrap()
            > original.state().vapor_mass_fraction().unwrap()
    );
    close(
        report.remaining_specific_internal_energy_change_j_kg,
        -12_591.965_499_954_124,
        1e-4,
    );
    close(
        report.remaining_specific_enthalpy_change_j_kg,
        -9885.548_034_785_53,
        1e-4,
    );
    // h may be negative on the original formation convention. Its conserved
    // transport and the gauge-independent finite-portion adjustment still matter.
    close(report.exported_enthalpy_j, -8277.049_580_656_221, 3e-4);
    close(
        report.mean_exported_enthalpy_j_kg.unwrap(),
        -20_501.998_248_903_896,
        8e-4,
    );
    close(
        report.exported_enthalpy_adjustment_j,
        -1748.002_274_066_815_5,
        2e-4,
    );
    close(report.exported_entropy_j_k, 2223.583_934_597_337, 2e-5);
    close(
        vessel.internal_energy_j() - original.internal_energy_j() + report.exported_enthalpy_j,
        0.0,
        2e-5,
    );
    close(
        vessel.mass_kg() * vessel.state().specific_entropy_j_kg_k()
            - original.mass_kg() * original.state().specific_entropy_j_kg_k()
            + report.exported_entropy_j_k,
        report.entropy_balance_residual_j_k,
        1e-8,
    );
    assert_eq!(report.stored_energy_change_j, -report.exported_enthalpy_j);
    assert!(
        report.entropy_balance_residual_j_k.abs()
            <= vessel.mass_kg() * controls().equilibrium.phase.entropy_tolerance_j_kg_k
    );
    assert!(
        report.volume_residual_m3.abs()
            <= vessel.mass_kg() * controls().equilibrium.specific_volume_tolerance_m3_kg
    );
    assert!(report.mass_residual_kg.abs() <= controls().mass_tolerance_kg);
    assert!(report.pt_evaluations <= controls().equilibrium.maximum_pt_evaluations);
    assert!(report.pressure_evaluations <= controls().equilibrium.maximum_pressure_evaluations);
}

#[test]
fn g1_dilute_discharge_recovers_analytic_blowdown_and_mean_enthalpy() {
    let mut control = MixtureEquilibriumDischargeControl::new(PhaseEntropyVolumeControls::new(
        [Pressure::new(0.6), Pressure::new(1.05)],
        PhaseEntropyControls::new([Temperature::new(290.0), Temperature::new(360.0)]),
    ));
    control.equilibrium.specific_volume_tolerance_m3_kg = 1e-5;
    control.equilibrium.pressure_tolerance_pa = 1e-12;
    let mut vessel = IsochoricMixtureEquilibriumParcel::from_pressure_temperature(
        source(),
        1.0,
        Temperature::new(350.0),
        Pressure::new(1.0),
        control.equilibrium.phase.phase_options(),
        || false,
    )
    .unwrap();
    let original = vessel.clone();
    let report = vessel
        .discharge_adiabatic(0.2, Pressure::new(0.5), control, || false)
        .unwrap();
    let ratio = 0.8_f64;
    let gamma = 4.0 / 3.0;
    let cv = 3.0 * R / MIXTURE_MASS;
    let expected_temperature = 350.0 * ratio.powf(gamma - 1.0);
    let expected_pressure = ratio.powf(gamma);
    let expected_mean_h = cv * 350.0 * (1.0 - ratio.powf(gamma)) / (1.0 - ratio);
    close(vessel.state().temperature_k(), expected_temperature, 1e-4);
    close(vessel.state().pressure_pa(), expected_pressure, 3e-7);
    close(
        report.mean_exported_enthalpy_j_kg.unwrap(),
        expected_mean_h,
        0.1,
    );
    assert!(
        report.mean_exported_enthalpy_j_kg.unwrap()
            < original.state().specific_enthalpy_j_kg() - 10_000.0
    );
    close(
        report.remaining_specific_internal_energy_change_j_kg,
        cv * (expected_temperature - 350.0),
        0.1,
    );
    assert_eq!(vessel.volume_m3(), original.volume_m3());
    assert_eq!(vessel.mass_kg(), 0.8);
}

#[test]
fn g3_discharge_cooling_and_transport_adjustment_ignore_species_datums() {
    let mut plain = parcel(source());
    let mut shifted = parcel(model([1e22, -3e21], [2e18, -7e17], [100.0, 450.0]));
    let portion = portion_to_independent_endpoint(&plain);
    let a = plain
        .discharge_adiabatic(portion, Pressure::new(1.5e6), controls(), || false)
        .unwrap();
    let b = shifted
        .discharge_adiabatic(portion, Pressure::new(1.5e6), controls(), || false)
        .unwrap();
    assert_eq!(
        plain.state().temperature_k(),
        shifted.state().temperature_k()
    );
    assert_eq!(plain.state().pressure_pa(), shifted.state().pressure_pa());
    assert_eq!(plain.state().phase(), shifted.state().phase());
    assert_eq!(plain.mass_kg(), shifted.mass_kg());
    assert_eq!(plain.volume_m3(), shifted.volume_m3());
    assert_eq!(
        a.remaining_specific_internal_energy_change_j_kg,
        b.remaining_specific_internal_energy_change_j_kg
    );
    assert_eq!(
        a.remaining_specific_enthalpy_change_j_kg,
        b.remaining_specific_enthalpy_change_j_kg
    );
    assert_eq!(
        a.exported_enthalpy_adjustment_j,
        b.exported_enthalpy_adjustment_j
    );
    assert_eq!(
        a.mean_exported_enthalpy_change_j_kg,
        b.mean_exported_enthalpy_change_j_kg
    );
    assert_eq!(
        a.entropy_balance_residual_j_k,
        b.entropy_balance_residual_j_k
    );
    assert_eq!(a.pt_evaluations, b.pt_evaluations);
}

#[test]
fn g1_subdivided_discharge_retains_endpoint_and_total_exported_enthalpy() {
    let mut single = parcel(source());
    let mut subdivided = single.clone();
    let portion = portion_to_independent_endpoint(&single);
    let whole = single
        .discharge_adiabatic(portion, Pressure::new(1.5e6), controls(), || false)
        .unwrap();
    let first = subdivided
        .discharge_adiabatic(0.5 * portion, Pressure::new(1.5e6), controls(), || false)
        .unwrap();
    let second = subdivided
        .discharge_adiabatic(0.5 * portion, Pressure::new(1.5e6), controls(), || false)
        .unwrap();
    close(
        subdivided.state().temperature_k(),
        single.state().temperature_k(),
        3e-6,
    );
    close(
        subdivided.state().pressure_pa(),
        single.state().pressure_pa(),
        0.3,
    );
    close(subdivided.mass_kg(), single.mass_kg(), 3e-16);
    assert_eq!(subdivided.volume_m3(), single.volume_m3());
    close(
        first.exported_enthalpy_j + second.exported_enthalpy_j,
        whole.exported_enthalpy_j,
        3e-4,
    );
    close(
        first.exported_entropy_j_k + second.exported_entropy_j_k,
        whole.exported_entropy_j_k,
        3e-8,
    );
}

#[test]
fn g0_discharge_rest_limits_and_pressure_refusals_are_atomic() {
    let mut vessel = parcel(source());
    let original = vessel.clone();
    let portion = portion_to_independent_endpoint(&original);
    let mut capped = controls();
    capped.equilibrium.maximum_pt_evaluations = 1;
    capped.equilibrium.maximum_pressure_evaluations = 1;
    capped.equilibrium.phase.flash_options.log_k_tolerance *= 2.0;
    let rest = vessel
        .discharge_adiabatic(0.0, Pressure::new(4e6), capped, || false)
        .unwrap();
    assert_eq!(vessel, original);
    assert_eq!(rest.pt_evaluations, 0);
    assert_eq!(rest.pressure_evaluations, 0);
    assert_eq!(rest.exported_enthalpy_j, 0.0);
    assert_eq!(rest.mean_exported_enthalpy_j_kg, None);
    for invalid in [
        -1.0,
        f64::NAN,
        f64::INFINITY,
        1.0,
        1.1,
        f64::EPSILON * 0.125,
    ] {
        assert!(
            vessel
                .discharge_adiabatic(invalid, Pressure::new(1.5e6), controls(), || false)
                .is_err()
        );
        assert_eq!(vessel, original);
    }
    assert!(matches!(
        vessel.discharge_adiabatic(portion, Pressure::new(1.5e6), capped, || false),
        Err(E::Invalid(_))
    ));
    assert_eq!(vessel, original);
    for pressure in [0.0, f64::NAN, 3e6, 4e6, 2.5e6] {
        assert!(
            vessel
                .discharge_adiabatic(portion, Pressure::new(pressure), controls(), || false)
                .is_err()
        );
        assert_eq!(vessel, original);
    }
    capped = controls();
    capped.equilibrium.maximum_pt_evaluations = 1;
    assert!(matches!(
        vessel.discharge_adiabatic(portion, Pressure::new(1.5e6), capped, || false),
        Err(E::WorkLimit {
            flash_evaluations: 1
        })
    ));
    assert_eq!(vessel, original);
    capped = controls();
    capped.equilibrium.maximum_pressure_evaluations = 1;
    assert!(matches!(
        vessel.discharge_adiabatic(portion, Pressure::new(1.5e6), capped, || false),
        Err(E::PressureWorkLimit {
            pressure_evaluations: 1
        })
    ));
    assert_eq!(vessel, original);
}

#[test]
fn g4_discharge_final_cancellation_keeps_mass_geometry_and_phase_state() {
    let mut vessel = parcel(source());
    let original = vessel.clone();
    let portion = portion_to_independent_endpoint(&original);
    let mut accepted = original.clone();
    let mut polls = 0;
    let report = accepted
        .discharge_adiabatic(portion, Pressure::new(1.5e6), controls(), || {
            polls += 1;
            false
        })
        .unwrap();
    let mut seen = 0;
    assert_eq!(
        vessel.discharge_adiabatic(portion, Pressure::new(1.5e6), controls(), || {
            seen += 1;
            seen == polls
        }),
        Err(E::Cancelled),
    );
    assert_eq!(seen, polls);
    assert_eq!(vessel, original);
    assert_eq!(
        vessel
            .discharge_adiabatic(portion, Pressure::new(1.5e6), controls(), || false)
            .unwrap(),
        report,
    );
    assert_eq!(vessel, accepted);
}
