//! Closed binary-vessel heating through the actual phase U,V owner.
use super::{close, common};
use common::{MIXTURE_MASS, model, source};
use fs_material::thermochemical_gas::nonideal::ThermochemicalPengRobinsonModel;
use fs_material::thermochemical_gas::nonideal::flash::{
    ThermochemicalFlashError as E,
    equilibrium::{PhaseEnthalpyControls, volume::PhaseEnergyVolumeControls},
    parcel::IsochoricMixtureEquilibriumParcel,
};
use fs_qty::{Pressure, Temperature};

fn controls() -> PhaseEnergyVolumeControls {
    PhaseEnergyVolumeControls::new(
        [Pressure::new(2e6), Pressure::new(8e6)],
        PhaseEnthalpyControls::new([Temperature::new(175.0), Temperature::new(325.0)]),
    )
}

fn parcel(source: ThermochemicalPengRobinsonModel, mass: f64) -> IsochoricMixtureEquilibriumParcel {
    IsochoricMixtureEquilibriumParcel::from_pressure_temperature(
        source,
        mass,
        Temperature::new(220.0),
        Pressure::new(3e6),
        controls().phase.phase_options(),
        || false,
    )
    .unwrap()
}

#[test]
fn g1_closed_binary_vessel_heats_through_phase_change_and_recovers_on_cooling() {
    // Independent coupled cubic/fugacity roots, with phase amounts from the
    // binary material balance and volume held fixed. NASA Cp/R is3.5/4.5.
    // At220K,3MPa: Vbar=0.00017620507091992501m3/mol, u=-39099.94554512886J/kg.
    // At280K on that isochore: p=6938572.708280131Pa, u=178943.68291362084J/kg.
    const DELTA_U: f64 = 218_043.628_458_749_7;
    let mut vessel = parcel(source(), 2.0 * MIXTURE_MASS);
    let initial = vessel.clone();
    assert!(initial.state().phase().two_phase().is_some());
    close(
        initial.volume_m3(),
        2.0 * 0.000_176_205_070_919_925_01,
        1e-14,
    );
    let heat = vessel.mass_kg() * DELTA_U;
    let report = vessel.add_heat(heat, controls(), || false).unwrap();
    close(vessel.state().temperature_k(), 280.0, 2e-6);
    close(vessel.state().pressure_pa(), 6_938_572.708_280_131, 0.2);
    assert!(vessel.state().phase().two_phase().is_none());
    assert!(report.entropy_change_j_k > 0.0);
    close(report.stored_energy_change_j, heat, vessel.mass_kg() * 1e-5);
    assert!(report.energy_residual_j.abs() <= vessel.mass_kg() * 1e-5);
    assert!(report.volume_residual_m3.abs() <= vessel.mass_kg() * 1e-12);
    assert!(report.pt_evaluations >= report.pressure_evaluations);
    assert!(report.eos_evaluation_upper_bound > report.pt_evaluations);
    close(
        vessel.internal_energy_j() - initial.internal_energy_j(),
        heat,
        2e-6,
    );
    assert_eq!(vessel.mass_kg(), initial.mass_kg());
    assert_eq!(vessel.volume_m3(), initial.volume_m3());
    assert_eq!(
        vessel.density_kg_m3().to_bits(),
        initial.density_kg_m3().to_bits()
    );
    assert_eq!(
        vessel.state().constitutive_law(),
        initial.state().constitutive_law()
    );

    let cooled = vessel.add_heat(-heat, controls(), || false).unwrap();
    close(vessel.state().temperature_k(), 220.0, 2e-6);
    close(vessel.state().pressure_pa(), 3e6, 0.2);
    assert!(vessel.state().phase().two_phase().is_some());
    close(
        vessel.state().vapor_mass_fraction().unwrap(),
        initial.state().vapor_mass_fraction().unwrap(),
        1e-7,
    );
    assert!(cooled.entropy_change_j_k < 0.0);
    close(
        report.entropy_change_j_k + cooled.entropy_change_j_k,
        0.0,
        2e-7,
    );
}

#[test]
fn g3_partitioned_heat_and_extensive_scaling_preserve_the_original_inventory() {
    const DELTA_U: f64 = 37_282.242_783_721_6;
    let mut whole = parcel(source(), 2.0 * MIXTURE_MASS);
    let initial = whole.clone();
    let mut partitioned = initial.clone();
    let mut scaled = parcel(source(), 6.0 * MIXTURE_MASS);
    let heat = whole.mass_kg() * DELTA_U;
    whole.add_heat(heat, controls(), || false).unwrap();
    let mut stored = 0.0;
    for _ in 0..4 {
        stored += partitioned
            .add_heat(heat / 4.0, controls(), || false)
            .unwrap()
            .stored_energy_change_j;
    }
    scaled.add_heat(3.0 * heat, controls(), || false).unwrap();
    close(whole.state().temperature_k(), 230.0, 2e-6);
    close(whole.state().pressure_pa(), 3_512_821.298_107_051_3, 0.2);
    close(stored, heat, 4.0 * whole.mass_kg() * 1e-5);
    for state in [&partitioned, &scaled] {
        close(
            state.state().temperature_k(),
            whole.state().temperature_k(),
            2e-6,
        );
        close(
            state.state().pressure_pa(),
            whole.state().pressure_pa(),
            0.2,
        );
        close(
            state.state().vapor_mass_fraction().unwrap(),
            whole.state().vapor_mass_fraction().unwrap(),
            1e-7,
        );
    }
    assert_eq!(partitioned.mass_kg(), initial.mass_kg());
    assert_eq!(partitioned.volume_m3(), initial.volume_m3());
    assert_eq!(partitioned.density_kg_m3(), initial.density_kg_m3());
    close(scaled.volume_m3(), 3.0 * initial.volume_m3(), 1e-18);
}

#[test]
fn g3_small_continuing_heat_is_invariant_to_independent_species_energy_datums() {
    let mut plain = parcel(source(), 2.0 * MIXTURE_MASS);
    let mut shifted = parcel(
        model([1e22, -3e21], [2e18, -7e17], [100.0, 450.0]),
        2.0 * MIXTURE_MASS,
    );
    for _ in 0..3 {
        let t0 = plain.state().temperature_k();
        let a = plain.add_heat(0.125, controls(), || false).unwrap();
        let b = shifted.add_heat(0.125, controls(), || false).unwrap();
        assert!(plain.state().temperature_k() > t0);
        assert_eq!(a, b);
        assert_eq!(
            plain.state().temperature_k(),
            shifted.state().temperature_k()
        );
        assert_eq!(plain.state().pressure_pa(), shifted.state().pressure_pa());
        assert_eq!(plain.state().phase(), shifted.state().phase());
        close(a.stored_energy_change_j, 0.125, plain.mass_kg() * 1e-5);
    }
}

#[test]
fn g0_conserved_initialization_and_unrepresentable_inventory_refusal() {
    let original = parcel(source(), 2.0 * MIXTURE_MASS);
    let restored = IsochoricMixtureEquilibriumParcel::from_internal_energy(
        source(),
        original.mass_kg(),
        original.volume_m3(),
        original.state().specific_internal_energy_j_kg(),
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(original.mass_kg(), restored.mass_kg());
    assert_eq!(original.volume_m3(), restored.volume_m3());
    close(restored.state().temperature_k(), 220.0, 2e-6);
    close(restored.state().pressure_pa(), 3e6, 0.2);
    for (mass, volume) in [
        (0.0, 1.0),
        (1.0, 0.0),
        (-1.0, 1.0),
        (f64::NAN, 1.0),
        (1e-300, 1e300),
    ] {
        assert!(
            IsochoricMixtureEquilibriumParcel::from_internal_energy(
                source(),
                mass,
                volume,
                0.0,
                controls(),
                || false,
            )
            .is_err()
        );
    }
    assert!(
        IsochoricMixtureEquilibriumParcel::from_pressure_temperature(
            source(),
            1e308,
            Temperature::new(220.0),
            Pressure::new(3e6),
            controls().phase.phase_options(),
            || false,
        )
        .is_err()
    );
}

#[test]
fn g4_failed_work_and_final_cancellation_preserve_the_complete_parcel() {
    let mut vessel = parcel(source(), 2.0 * MIXTURE_MASS);
    let before = vessel.clone();
    let mut capped = controls();
    capped.maximum_pt_evaluations = 1;
    assert!(vessel.add_heat(1.0, capped, || false).is_err());
    assert_eq!(vessel, before);
    for heat in [f64::NAN, f64::INFINITY, -f64::INFINITY, 1e-300] {
        assert!(vessel.add_heat(heat, controls(), || false).is_err());
        assert_eq!(vessel, before);
    }
    let zero = vessel.add_heat(0.0, capped, || false).unwrap();
    assert_eq!(zero.stored_energy_change_j, 0.0);
    assert_eq!(zero.pt_evaluations, 0);
    assert_eq!(zero.pressure_evaluations, 0);
    assert_eq!(zero.eos_evaluation_upper_bound, 0);
    assert_eq!(vessel, before);
    let mut counted = before.clone();
    let mut polls = 0;
    counted
        .add_heat(1.0, controls(), || {
            polls += 1;
            false
        })
        .unwrap();
    let mut repeated = 0;
    assert_eq!(
        vessel.add_heat(1.0, controls(), || {
            repeated += 1;
            repeated == polls
        }),
        Err(E::Cancelled),
    );
    assert_eq!(repeated, polls);
    assert_eq!(vessel, before);
    let retried = vessel.add_heat(1.0, controls(), || false).unwrap();
    assert_eq!(vessel, counted);
    assert!(retried.stored_energy_change_j > 0.0);
}
