//! G0/G1/G3/G4 first-law-coupled heat/work/volume equilibrium.
use super::{close, common};
use common::{MIXTURE_MASS, model, source};
use fs_material::thermochemical_gas::nonideal::flash::{
    equilibrium::{
        PhaseEnthalpyControls,
        volume::PhaseEnergyVolumeControls,
    },
    parcel::driven::DrivenMixtureEquilibriumParcel,
};
use fs_qty::{Pressure, Temperature};

fn controls() -> PhaseEnergyVolumeControls {
    PhaseEnergyVolumeControls::new(
        [Pressure::new(1e6), Pressure::new(5e6)],
        PhaseEnthalpyControls::new([Temperature::new(160.0), Temperature::new(340.0)]),
    )
}

#[test]
fn g1_heat_and_boundary_work_recover_known_phase_equilibrium_endpoint() {
    let owner = source();
    let mass = 0.2;
    let mut parcel = DrivenMixtureEquilibriumParcel::from_pressure_temperature(
        owner,
        mass,
        Temperature::new(220.0),
        Pressure::new(3e6),
        controls().phase.phase_options(),
        || false,
    )
    .unwrap();
    let initial_volume = parcel.volume_m3();
    // Independent high-precision PT references already pinned by the parent
    // nonideal-flash suite: U0=-39099.945545128794 J/kg and
    // U1=264458.43920078234 J/kg, rho1=36.471975375322785 kg/m3.
    let delta_u_j = mass * 303_558.384_745_911_1;
    let final_volume = mass / 36.471_975_375_322_785;
    let heat = 80_000.0;
    let work_by = heat - delta_u_j;
    let step = parcel
        .advance(final_volume, heat, work_by, controls(), || false)
        .unwrap();

    close(parcel.state().temperature_k(), 280.0, 5e-7);
    close(parcel.state().pressure_pa(), 3e6, 0.03);
    close(parcel.volume_m3(), final_volume, 0.0);
    close(step.volume_change_m3, final_volume - initial_volume, 0.0);
    close(step.requested_internal_energy_change_j, delta_u_j, 2e-11);
    close(step.stored_energy_change_j, delta_u_j, 3e-6);
    assert!(step.first_law_residual_j.abs() <= mass * controls().internal_energy_tolerance_j_kg);
    assert!(step.volume_residual_m3.abs() <= mass * controls().specific_volume_tolerance_m3_kg);
    assert!(parcel.state().phase().two_phase().is_none());
}

#[test]
fn g3_same_first_law_target_is_invariant_to_large_species_energy_datums() {
    let plain = source();
    let shifted = model([1e22, -3e21], [2e18, -7e17], [100.0, 450.0]);
    let mass = 0.15;
    let mut a = DrivenMixtureEquilibriumParcel::from_pressure_temperature(
        plain,
        mass,
        Temperature::new(220.0),
        Pressure::new(3e6),
        controls().phase.phase_options(),
        || false,
    )
    .unwrap();
    let mut b = DrivenMixtureEquilibriumParcel::from_pressure_temperature(
        shifted,
        mass,
        Temperature::new(220.0),
        Pressure::new(3e6),
        controls().phase.phase_options(),
        || false,
    )
    .unwrap();
    let final_volume = mass / 36.471_975_375_322_785;
    let delta_u = mass * 303_558.384_745_911_1;
    let heat = 60_000.0;
    let work = heat - delta_u;
    let ra = a.advance(final_volume, heat, work, controls(), || false).unwrap();
    let rb = b.advance(final_volume, heat, work, controls(), || false).unwrap();

    assert_eq!(a.state().temperature_k().to_bits(), b.state().temperature_k().to_bits());
    assert_eq!(a.state().pressure_pa().to_bits(), b.state().pressure_pa().to_bits());
    assert_eq!(a.state().density_kg_m3().to_bits(), b.state().density_kg_m3().to_bits());
    assert_eq!(ra.stored_energy_change_j.to_bits(), rb.stored_energy_change_j.to_bits());
    assert_eq!(ra.entropy_change_j_k.to_bits(), rb.entropy_change_j_k.to_bits());
}

#[test]
fn g0_exact_rest_is_bitwise_and_skips_all_phase_work() {
    let mut parcel = DrivenMixtureEquilibriumParcel::from_pressure_temperature(
        source(),
        0.1,
        Temperature::new(220.0),
        Pressure::new(3e6),
        controls().phase.phase_options(),
        || false,
    )
    .unwrap();
    let before = parcel.clone();
    let step = parcel
        .advance(parcel.volume_m3(), 0.0, 0.0, controls(), || false)
        .unwrap();
    assert_eq!(parcel, before);
    assert_eq!(step.stored_energy_change_j, 0.0);
    assert_eq!(step.first_law_residual_j, 0.0);
    assert_eq!(step.pressure_evaluations, 0);
    assert_eq!(step.pt_evaluations, 0);
    assert_eq!(step.eos_evaluation_upper_bound, 0);
}

#[test]
fn g4_refusal_and_final_cancellation_leave_complete_state_unchanged() {
    let make = || {
        DrivenMixtureEquilibriumParcel::from_pressure_temperature(
            source(),
            0.1,
            Temperature::new(220.0),
            Pressure::new(3e6),
            controls().phase.phase_options(),
            || false,
        )
        .unwrap()
    };
    let mut parcel = make();
    let before = parcel.clone();
    assert!(parcel
        .advance(parcel.volume_m3(), f64::NAN, 0.0, controls(), || false)
        .is_err());
    assert_eq!(parcel, before);

    let final_volume = 0.1 / 36.471_975_375_322_785;
    let delta_u = 0.1 * 303_558.384_745_911_1;
    let heat = 40_000.0;
    let work = heat - delta_u;
    let mut control = make();
    let mut polls = 0usize;
    control
        .advance(final_volume, heat, work, controls(), || {
            polls += 1;
            false
        })
        .unwrap();
    let mut replay = 0usize;
    assert!(parcel
        .advance(final_volume, heat, work, controls(), || {
            replay += 1;
            replay == polls
        })
        .is_err());
    assert_eq!(replay, polls);
    assert_eq!(parcel, before);

    parcel
        .advance(final_volume, heat, work, controls(), || false)
        .unwrap();
    assert_eq!(parcel, control);
}
