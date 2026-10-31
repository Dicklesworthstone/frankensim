//! Closed-vessel conservation through the actual production equilibrium owner.
#[path = "common/nonideal_equilibrium.rs"]
mod common;
use common::{MOLAR_MASS, controls, model};
use fs_material::thermochemical_gas::nonideal::equilibrium::{
    PureFluidEquilibriumError as E, PureFluidEquilibriumPhase as Phase,
    parcel::IsochoricEquilibriumParcel,
};
use fs_qty::Temperature;

fn parcel(datum: f64, molar_volume: f64) -> IsochoricEquilibriumParcel {
    IsochoricEquilibriumParcel::new(
        model(datum, [100.0, 250.0]),
        2.0 * MOLAR_MASS,
        2.0 * molar_volume,
        Temperature::new(140.0),
        controls(),
        || false,
    )
    .unwrap()
}

#[test]
fn g1_closed_vessel_heat_cycles_cross_both_phase_boundaries_and_conserve_joules() {
    // Independent original-PR Maxwell-area/energy oracle: heat from140K
    // coexistence to170K vapor or liquid at the same two conserved molar volumes.
    for (molar_volume, heat_j_kg, vapor) in [
        (0.0005, 280_090.004_666_026_16, true),
        (0.00005, 119_577.435_965_805_97, false),
    ] {
        let mut whole = parcel(0.0, molar_volume);
        let initial = whole.clone();
        let mut partitioned = initial.clone();
        assert!(matches!(initial.state().phase(), Phase::Coexistence(_)));
        let heat = whole.mass_kg() * heat_j_kg;
        let tolerance = whole.mass_kg() * controls().energy_tolerance_j_kg;
        let report = whole.add_heat(heat, controls(), || false).unwrap();
        assert_eq!(report.supplied_heat_j, heat);
        assert!((report.stored_energy_change_j - heat).abs() <= tolerance);
        assert!(report.residual_j.abs() <= tolerance);
        assert!(report.eos_evaluations > report.temperature_evaluations);
        assert!((whole.state().temperature_k() - 170.0).abs() < 3e-8);
        assert!(matches!(whole.state().phase(), Phase::Vapor(_)) == vapor);
        assert!(matches!(whole.state().phase(), Phase::Liquid(_)) != vapor);
        let mut stored = 0.0;
        for _ in 0..8 {
            stored += partitioned
                .add_heat(heat / 8.0, controls(), || false)
                .unwrap()
                .stored_energy_change_j;
        }
        assert!((stored - heat).abs() <= 8.0 * tolerance);
        assert!((partitioned.state().temperature_k() - whole.state().temperature_k()).abs() < 1e-7);
        assert!((partitioned.state().pressure_pa() - whole.state().pressure_pa()).abs() < 0.02);
        assert!(
            (whole.internal_energy_j() - initial.internal_energy_j() - heat).abs()
                < 2.0 * tolerance
        );
        for state in [&whole, &partitioned] {
            assert_eq!(state.mass_kg(), initial.mass_kg());
            assert_eq!(state.volume_m3(), initial.volume_m3());
            assert_eq!(
                state.density_kg_m3().to_bits(),
                initial.density_kg_m3().to_bits()
            );
            assert_eq!(
                state.state().density_kg_m3().to_bits(),
                initial.density_kg_m3().to_bits()
            );
        }
        whole.add_heat(-heat, controls(), || false).unwrap();
        assert!((whole.state().temperature_k() - 140.0).abs() < 5e-8);
        assert!(matches!(whole.state().phase(), Phase::Coexistence(_)));
        assert!(
            (whole.state().vapor_mass_fraction() - initial.state().vapor_mass_fraction()).abs()
                < 1e-8
        );
    }
}

#[test]
fn g3_continuing_heat_retains_small_increments_under_large_formation_references() {
    let mut base = parcel(0.0, 0.0005);
    let mut shifted = [parcel(-1e20, 0.0005), parcel(1e20, 0.0005)];
    for _ in 0..16 {
        let before = base.state().temperature_k();
        let accepted = base.add_heat(0.125, controls(), || false).unwrap();
        assert!(base.state().temperature_k() > before);
        assert!(
            (accepted.stored_energy_change_j - 0.125).abs()
                <= base.mass_kg() * controls().energy_tolerance_j_kg
        );
        for parcel in &mut shifted {
            let receipt = parcel.add_heat(0.125, controls(), || false).unwrap();
            assert_eq!(receipt, accepted);
            assert_eq!(parcel.state().temperature_k(), base.state().temperature_k());
            assert_eq!(parcel.state().pressure_pa(), base.state().pressure_pa());
            assert_eq!(parcel.state().phase(), base.state().phase());
        }
    }
}

#[test]
fn g4_phase_work_failure_and_final_cancellation_leave_the_complete_parcel_unchanged() {
    let mut vessel = parcel(0.0, 0.0005);
    let before = vessel.clone();
    let mut limited = controls();
    limited.maximum_eos_evaluations = 4;
    assert!(matches!(
        vessel.add_heat(1.0, limited, || false),
        Err(E::WorkLimit { .. })
    ));
    assert_eq!(vessel, before);
    for heat in [f64::NAN, f64::INFINITY, -f64::INFINITY, 1e-300] {
        assert!(vessel.add_heat(heat, controls(), || false).is_err());
        assert_eq!(vessel, before);
    }
    let zero = vessel.add_heat(0.0, controls(), || false).unwrap();
    assert_eq!(zero.stored_energy_change_j, 0.0);
    assert_eq!(zero.temperature_evaluations, 0);
    assert_eq!(zero.eos_evaluations, 0);
    assert_eq!(vessel, before);
    let mut count_run = before.clone();
    let mut polls = 0;
    count_run
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
        Err(E::Cancelled)
    );
    assert_eq!(repeated, polls);
    assert_eq!(vessel, before);
}
