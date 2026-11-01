//! G0/G1/G3/G4: prescribed-volume work over the original mixture equilibrium.
use super::{close, common};
use common::{METHANE, MIXTURE_MASS, model, source};
use fs_material::thermochemical_gas::nonideal::{
    ThermochemicalPengRobinsonModel,
    flash::{
        ThermochemicalFlashError as E,
        equilibrium::entropy::{PhaseEntropyControls, volume::PhaseEntropyVolumeControls},
        parcel::isentropic::IsentropicMixtureEquilibriumParcel,
    },
};
use fs_qty::{Pressure, Temperature};

const EXPANDED_DENSITY: f64 = 78.021_034_395_541_12;
const EXPANDED_TEMPERATURE: f64 = 208.202_765_109_369_35;
const WORK_PER_MASS: f64 = 12_591.965_499_954_124;

fn controls() -> PhaseEntropyVolumeControls {
    PhaseEntropyVolumeControls::new(
        [Pressure::new(1.8e6), Pressure::new(3.2e6)],
        PhaseEntropyControls::new([Temperature::new(195.0), Temperature::new(235.0)]),
    )
}

fn parcel(
    model: ThermochemicalPengRobinsonModel,
    mass_kg: f64,
) -> IsentropicMixtureEquilibriumParcel {
    IsentropicMixtureEquilibriumParcel::from_pressure_temperature(
        model,
        mass_kg,
        Temperature::new(220.0),
        Pressure::new(3e6),
        controls().phase.phase_options(),
        || false,
    )
    .unwrap()
}

#[test]
fn g1_volume_driven_expansion_and_compression_match_independent_phase_work() {
    let mut vessel = parcel(source(), 2.0 * MIXTURE_MASS);
    let initial = vessel.clone();
    // Independent simultaneous cubic/fugacity/entropy solution from the PS
    // oracle: 220 K, 3 MPa -> 208.20276510936935 K, 2 MPa. The caller prescribes
    // that independently computed volume, not the final pressure or temperature.
    let final_volume = vessel.mass_kg() / EXPANDED_DENSITY;
    let report = vessel
        .change_volume(final_volume, controls(), || false)
        .unwrap();
    close(vessel.state().temperature_k(), EXPANDED_TEMPERATURE, 2e-6);
    close(vessel.state().pressure_pa(), 2e6, 0.2);
    close(vessel.state().density_kg_m3(), EXPANDED_DENSITY, 2e-6);
    let split = vessel.state().phase().two_phase().unwrap();
    close(split.vapor_fraction_mol(), 0.372_659_298_790_543_3, 1e-8);
    close(
        split.liquid_mole_fractions()[METHANE],
        0.301_426_580_035_008_8,
        1e-8,
    );
    close(
        split.vapor_mole_fractions()[METHANE],
        0.834_281_712_348_777_5,
        1e-8,
    );
    close(
        report.boundary_work_j,
        vessel.mass_kg() * WORK_PER_MASS,
        2e-5,
    );
    close(
        report.enthalpy_change_j,
        vessel.mass_kg() * -9885.548_034_785_53,
        2e-5,
    );
    assert_eq!(report.boundary_work_j, -report.stored_energy_change_j);
    close(report.entropy_change_j_k, 0.0, vessel.mass_kg() * 1e-8);
    close(
        report.volume_change_m3,
        final_volume - initial.volume_m3(),
        0.0,
    );
    assert!(
        report.volume_residual_m3.abs()
            <= vessel.mass_kg() * controls().specific_volume_tolerance_m3_kg
    );
    assert_eq!(vessel.volume_m3().to_bits(), final_volume.to_bits());
    assert_eq!(vessel.mass_kg(), initial.mass_kg());
    assert_eq!(vessel.model(), initial.model());
    assert_eq!(
        vessel.state().constitutive_law(),
        initial.state().constitutive_law()
    );
    // Independent state-function identity on the actually published volumes.
    close(
        report.enthalpy_change_j - report.stored_energy_change_j,
        vessel.state().pressure_pa() * vessel.volume_m3()
            - initial.state().pressure_pa() * initial.volume_m3(),
        2e-6,
    );
    let reverse = vessel
        .change_volume(initial.volume_m3(), controls(), || false)
        .unwrap();
    assert!(reverse.boundary_work_j < 0.0);
    close(vessel.state().temperature_k(), 220.0, 2e-6);
    close(vessel.state().pressure_pa(), 3e6, 0.2);
    close(
        vessel.state().vapor_mass_fraction().unwrap(),
        initial.state().vapor_mass_fraction().unwrap(),
        1e-8,
    );
    close(report.boundary_work_j + reverse.boundary_work_j, 0.0, 4e-5);
    close(
        report.enthalpy_change_j + reverse.enthalpy_change_j,
        0.0,
        4e-5,
    );
    assert_eq!(vessel.volume_m3(), initial.volume_m3());
}

#[test]
fn g1_pressure_volume_quadrature_converges_to_source_work_and_scales_extensively() {
    let initial = parcel(source(), 2.0 * MIXTURE_MASS);
    let final_volume = initial.mass_kg() / EXPANDED_DENSITY;
    let expected_work = initial.mass_kg() * WORK_PER_MASS;
    let mut errors = Vec::new();
    for steps in [1, 2, 4] {
        let mut vessel = initial.clone();
        let mut integrated_pressure_work = 0.0;
        let mut source_work = 0.0;
        for step in 1..=steps {
            let next = if step == steps {
                final_volume
            } else {
                initial.volume_m3()
                    + (final_volume - initial.volume_m3()) * step as f64 / steps as f64
            };
            let old_volume = vessel.volume_m3();
            let old_pressure = vessel.state().pressure_pa();
            let report = vessel.change_volume(next, controls(), || false).unwrap();
            integrated_pressure_work +=
                0.5 * (old_pressure + vessel.state().pressure_pa()) * (next - old_volume);
            source_work += report.boundary_work_j;
        }
        close(source_work, expected_work, 8e-5);
        close(vessel.state().temperature_k(), EXPANDED_TEMPERATURE, 2e-6);
        close(vessel.state().pressure_pa(), 2e6, 0.2);
        assert_eq!(vessel.mass_kg(), initial.mass_kg());
        assert_eq!(vessel.model(), initial.model());
        errors.push((integrated_pressure_work - expected_work).abs());
    }
    // Independent trapezoidal p*dV integration has a live, second-order mesh.
    assert!(errors[1] < 0.5 * errors[0], "quadrature errors {errors:?}");
    assert!(errors[2] < 0.5 * errors[1], "quadrature errors {errors:?}");
    assert!(
        errors[2] < 0.005 * expected_work,
        "quadrature errors {errors:?}"
    );
    let mut scaled = parcel(source(), 6.0 * MIXTURE_MASS);
    let report = scaled
        .change_volume(3.0 * final_volume, controls(), || false)
        .unwrap();
    close(report.boundary_work_j, 3.0 * expected_work, 6e-5);
    close(scaled.state().temperature_k(), EXPANDED_TEMPERATURE, 2e-6);
    close(scaled.state().pressure_pa(), 2e6, 0.2);
    assert_eq!(scaled.mass_kg(), 3.0 * initial.mass_kg());
    assert_eq!(scaled.volume_m3(), 3.0 * final_volume);
}

#[test]
fn g3_small_motion_survives_coarse_controls_and_independent_species_datums() {
    let mut plain = parcel(source(), 2.0 * MIXTURE_MASS);
    let mut shifted = parcel(
        model([1e22, -3e21], [2e18, -7e17], [100.0, 450.0]),
        2.0 * MIXTURE_MASS,
    );
    let initial_volume = plain.volume_m3();
    let mut coarse = controls();
    coarse.phase.entropy_tolerance_j_kg_k = 1.0;
    coarse.phase.temperature_tolerance_k = 1.0;
    coarse.specific_volume_tolerance_m3_kg = 1e-3;
    coarse.pressure_tolerance_pa = 1e5;
    for factor in [1.000_001, 1.0] {
        let next = initial_volume * factor;
        let pressure = plain.state().pressure_pa();
        let volume_change = next - plain.volume_m3();
        let a = plain.change_volume(next, coarse, || false).unwrap();
        let b = shifted.change_volume(next, coarse, || false).unwrap();
        assert_eq!(a, b);
        assert_ne!(a.boundary_work_j, 0.0);
        assert_eq!(a.boundary_work_j.signum(), volume_change.signum());
        // Differential first-law reference for this tiny motion; residual
        // safeguards reserve resolution despite the deliberately coarse inputs.
        close(
            a.boundary_work_j,
            pressure * volume_change,
            0.25 * (pressure * volume_change).abs(),
        );
        assert_eq!(
            plain.state().temperature_k(),
            shifted.state().temperature_k()
        );
        assert_eq!(plain.state().pressure_pa(), shifted.state().pressure_pa());
        assert_eq!(plain.state().phase(), shifted.state().phase());
        assert_eq!(plain.volume_m3(), next);
        assert_eq!(plain.mass_kg(), 2.0 * MIXTURE_MASS);
        assert_eq!(shifted.mass_kg(), plain.mass_kg());
    }
}

#[test]
fn g0_entropy_initialization_and_exact_rest_retain_the_full_inventory() {
    let original = parcel(source(), 2.0 * MIXTURE_MASS);
    let mut restored = IsentropicMixtureEquilibriumParcel::from_entropy(
        source(),
        original.mass_kg(),
        original.volume_m3(),
        original.state().specific_entropy_j_kg_k(),
        controls(),
        || false,
    )
    .unwrap();
    close(restored.state().temperature_k(), 220.0, 2e-6);
    close(restored.state().pressure_pa(), 3e6, 0.2);
    assert_eq!(restored.mass_kg(), original.mass_kg());
    assert_eq!(restored.volume_m3(), original.volume_m3());
    close(
        restored.internal_energy_j(),
        original.internal_energy_j(),
        2e-5,
    );
    close(restored.enthalpy_j(), original.enthalpy_j(), 2e-5);
    close(restored.entropy_j_k(), original.entropy_j_k(), 2e-8);
    let before = restored.clone();
    let mut capped = controls();
    capped.maximum_pressure_evaluations = 1;
    capped.maximum_pt_evaluations = 1;
    let report = restored
        .change_volume(restored.volume_m3(), capped, || false)
        .unwrap();
    assert_eq!(restored, before);
    assert_eq!(report.stored_energy_change_j, 0.0);
    assert_eq!(report.enthalpy_change_j, 0.0);
    assert_eq!(report.boundary_work_j, 0.0);
    assert_eq!(report.entropy_change_j_k, 0.0);
    assert_eq!(report.volume_change_m3, 0.0);
    assert_eq!(report.pressure_evaluations, 0);
    assert_eq!(report.pt_evaluations, 0);
    assert_eq!(report.eos_evaluation_upper_bound, 0);
}

#[test]
fn g0_invalid_or_unresolvable_inventory_and_motion_refuse_atomically() {
    let mut vessel = parcel(source(), 2.0 * MIXTURE_MASS);
    let original = vessel.clone();
    for volume in [
        0.0,
        -1.0,
        f64::NAN,
        f64::INFINITY,
        1e-320,
        vessel.volume_m3().next_up(),
    ] {
        assert!(vessel.change_volume(volume, controls(), || false).is_err());
        assert_eq!(vessel, original);
    }
    let mut changed_phase_policy = controls();
    changed_phase_policy.phase.flash_options.log_k_tolerance *= 2.0;
    assert_eq!(
        vessel.change_volume(
            vessel.mass_kg() / EXPANDED_DENSITY,
            changed_phase_policy,
            || false,
        ),
        Err(E::Invalid(
            "continuing isentropic mixture requires its accepted original phase policy",
        )),
    );
    assert_eq!(vessel, original);
    // Rest does not reevaluate the source or install the changed phase policy.
    let rest = vessel
        .change_volume(vessel.volume_m3(), changed_phase_policy, || false)
        .unwrap();
    assert_eq!(rest.pt_evaluations, 0);
    assert_eq!(vessel, original);
    for (mass, volume) in [
        (0.0, 1.0),
        (-1.0, 1.0),
        (f64::NAN, 1.0),
        (1.0, 0.0),
        (1e-300, 1e300),
    ] {
        assert!(
            IsentropicMixtureEquilibriumParcel::from_entropy(
                source(),
                mass,
                volume,
                0.0,
                controls(),
                || false
            )
            .is_err()
        );
    }
    assert!(
        IsentropicMixtureEquilibriumParcel::from_pressure_temperature(
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
fn g4_work_limits_and_final_cancellation_preserve_all_accepted_state() {
    let mut vessel = parcel(source(), 2.0 * MIXTURE_MASS);
    let original = vessel.clone();
    let final_volume = vessel.mass_kg() / EXPANDED_DENSITY;
    let mut capped = controls();
    capped.maximum_pt_evaluations = 1;
    assert!(matches!(
        vessel.change_volume(final_volume, capped, || false),
        Err(E::WorkLimit {
            flash_evaluations: 1
        })
    ));
    assert_eq!(vessel, original);
    let mut pressure_capped = controls();
    pressure_capped.maximum_pressure_evaluations = 1;
    assert!(matches!(
        vessel.change_volume(final_volume, pressure_capped, || false),
        Err(E::PressureWorkLimit {
            pressure_evaluations: 1
        })
    ));
    assert_eq!(vessel, original);
    let mut counted = original.clone();
    let mut polls = 0;
    let expected_report = counted
        .change_volume(final_volume, controls(), || {
            polls += 1;
            false
        })
        .unwrap();
    let mut replay_polls = 0;
    assert_eq!(
        vessel.change_volume(final_volume, controls(), || {
            replay_polls += 1;
            replay_polls == polls
        }),
        Err(E::Cancelled)
    );
    assert_eq!(replay_polls, polls);
    assert_eq!(vessel, original);
    let retried = vessel
        .change_volume(final_volume, controls(), || false)
        .unwrap();
    assert_eq!(retried, expected_report);
    assert_eq!(vessel, counted);
}
