//! Constant-pressure phase change must export the actual expansion work.
use super::{close, common};
use common::{MIXTURE_MASS, model, source};
use fs_material::thermochemical_gas::nonideal::ThermochemicalPengRobinsonModel;
use fs_material::thermochemical_gas::nonideal::flash::{
    ThermochemicalFlashError as E, equilibrium::PhaseEnthalpyControls,
    parcel::isobaric::IsobaricMixtureEquilibriumParcel,
};
use fs_qty::{Pressure, Temperature};

fn controls() -> PhaseEnthalpyControls {
    PhaseEnthalpyControls::new([Temperature::new(175.0), Temperature::new(285.0)])
}

fn parcel(
    model: ThermochemicalPengRobinsonModel,
    mass: f64,
    temperature: f64,
) -> IsobaricMixtureEquilibriumParcel {
    IsobaricMixtureEquilibriumParcel::from_pressure_temperature(
        model,
        mass,
        Temperature::new(temperature),
        Pressure::new(3e6),
        controls().phase_options(),
        || false,
    )
    .unwrap()
}

#[test]
fn g1_isobaric_phase_change_closes_heat_internal_energy_and_expansion_work() {
    let mut body = parcel(source(), 2.0 * MIXTURE_MASS, 180.0);
    let initial = body.clone();
    // Independent 50-digit simultaneous PR cubic/fugacity references from the
    // existing PT battery. The bubble/dew roots are205.7464/250.6892K at3MPa.
    let stages = [
        (
            220.0,
            -16_172.250_165_027_297,
            -39_099.945_545_128_86,
            130.846_120_827_461_9,
            5507.749_287_855_677,
            true,
        ),
        (
            280.0,
            346_713.375_192_810_4,
            264_458.439_200_782_34,
            36.471_975_375_322_785,
            7011.990_134_658_653,
            false,
        ),
    ];
    let mut prior_h = -203_472.419_617_408_38;
    let mut prior_u = -209_511.548_518_424_48;
    let mut prior_s = 4589.992_780_575_417;
    let mut work = 0.0;
    for (t, h, u, rho, s, split) in stages {
        let old_volume = body.volume_m3();
        let heat = body.mass_kg() * (h - prior_h);
        let report = body.add_heat(heat, controls(), || false).unwrap();
        close(body.state().temperature_k(), t, 3e-7);
        close(body.volume_m3(), body.mass_kg() / rho, 3e-12);
        assert_eq!(body.state().phase().two_phase().is_some(), split);
        close(
            report.stored_energy_change_j,
            body.mass_kg() * (u - prior_u),
            2e-6,
        );
        close(
            report.entropy_change_j_k,
            body.mass_kg() * (s - prior_s),
            2e-7,
        );
        assert_eq!(report.volume_change_m3, body.volume_m3() - old_volume);
        assert_eq!(report.boundary_work_j, 3e6 * report.volume_change_m3);
        assert!(report.boundary_work_j > 0.0);
        assert!(report.energy_residual_j.abs() <= body.mass_kg() * 1e-5);
        assert!(report.enthalpy_residual_j.abs() <= body.mass_kg() * 1e-5);
        close(
            report.stored_energy_change_j + report.boundary_work_j,
            heat,
            body.mass_kg() * 1e-5,
        );
        assert_eq!(body.pressure_pa(), 3e6);
        assert_eq!(body.mass_kg(), initial.mass_kg());
        assert_eq!(body.model(), initial.model());
        assert!(report.pt_evaluations > 1);
        assert!(report.eos_evaluation_upper_bound > report.pt_evaluations);
        assert!(
            report.eos_evaluation_upper_bound <= controls().maximum_eos_evaluation_bound().unwrap()
        );
        work += report.boundary_work_j;
        prior_h = h;
        prior_u = u;
        prior_s = s;
    }
    let cooling = body
        .add_heat(-body.mass_kg() * 550_185.794_810_218_8, controls(), || {
            false
        })
        .unwrap();
    close(body.state().temperature_k(), 180.0, 3e-7);
    close(body.volume_m3(), initial.volume_m3(), 3e-12);
    close(work + cooling.boundary_work_j, 0.0, 2e-6);
    assert!(cooling.boundary_work_j < 0.0);
    assert!(cooling.entropy_change_j_k < 0.0);
}

#[test]
fn g3_partitioned_isobaric_heat_and_extensive_scaling_keep_the_same_endpoint() {
    let mut whole = parcel(source(), 2.0 * MIXTURE_MASS, 180.0);
    let mut partitioned = whole.clone();
    let mut scaled = parcel(source(), 6.0 * MIXTURE_MASS, 180.0);
    let heat = whole.mass_kg() * 550_185.794_810_218_8;
    let once = whole.add_heat(heat, controls(), || false).unwrap();
    let mut work = 0.0;
    let mut energy = 0.0;
    for _ in 0..4 {
        let step = partitioned
            .add_heat(heat / 4.0, controls(), || false)
            .unwrap();
        work += step.boundary_work_j;
        energy += step.stored_energy_change_j;
    }
    let big = scaled.add_heat(3.0 * heat, controls(), || false).unwrap();
    close(
        partitioned.state().temperature_k(),
        whole.state().temperature_k(),
        1e-6,
    );
    close(partitioned.volume_m3(), whole.volume_m3(), 3e-12);
    close(
        scaled.state().temperature_k(),
        whole.state().temperature_k(),
        1e-6,
    );
    close(scaled.volume_m3(), 3.0 * whole.volume_m3(), 9e-12);
    close(work, once.boundary_work_j, 5e-6);
    close(energy, once.stored_energy_change_j, 5e-6);
    close(big.boundary_work_j, 3.0 * once.boundary_work_j, 5e-6);
    close(
        big.stored_energy_change_j,
        3.0 * once.stored_energy_change_j,
        5e-6,
    );
}

#[test]
fn g3_small_heat_under_large_species_datums_resolves_original_work_and_energy() {
    let mut plain = parcel(source(), 2.0 * MIXTURE_MASS, 220.0);
    let mut shifted = parcel(
        model([1e22, -3e21], [2e18, -7e17], [100.0, 450.0]),
        2.0 * MIXTURE_MASS,
        220.0,
    );
    let mut coarse = controls();
    coarse.enthalpy.enthalpy_tolerance_j_kg = 1.0;
    for _ in 0..3 {
        let t0 = plain.state().temperature_k();
        let a = plain.add_heat(0.001, coarse, || false).unwrap();
        let b = shifted.add_heat(0.001, coarse, || false).unwrap();
        assert_eq!(a, b);
        assert_eq!(plain.state().phase(), shifted.state().phase());
        assert_eq!(plain.volume_m3(), shifted.volume_m3());
        assert!(plain.state().temperature_k() > t0);
        assert!(a.boundary_work_j > 0.0);
        close(a.enthalpy_change_j, 0.001, 0.00025);
        close(a.stored_energy_change_j + a.boundary_work_j, 0.001, 0.00025);
    }
}

#[test]
fn g0_isobaric_initialization_uses_original_enthalpy_and_refuses_bad_inventory() {
    let original = parcel(source(), 2.0 * MIXTURE_MASS, 220.0);
    let restored = IsobaricMixtureEquilibriumParcel::from_enthalpy(
        source(),
        original.mass_kg(),
        Pressure::new(3e6),
        original.state().specific_enthalpy_j_kg(),
        controls(),
        || false,
    )
    .unwrap();
    close(restored.state().temperature_k(), 220.0, 3e-7);
    close(restored.volume_m3(), original.volume_m3(), 3e-12);
    close(restored.enthalpy_j(), original.enthalpy_j(), 2e-6);
    close(
        restored.internal_energy_j(),
        original.internal_energy_j(),
        2e-6,
    );
    assert_eq!(restored.mass_kg(), original.mass_kg());
    for mass in [0.0, -1.0, f64::NAN, f64::INFINITY, 1e308] {
        assert!(
            IsobaricMixtureEquilibriumParcel::from_pressure_temperature(
                source(),
                mass,
                Temperature::new(220.0),
                Pressure::new(3e6),
                controls().phase_options(),
                || false
            )
            .is_err()
        );
    }
    assert!(
        IsobaricMixtureEquilibriumParcel::from_enthalpy(
            model([1e22, -3e21], [0.0; 2], [100.0, 450.0]),
            original.mass_kg(),
            Pressure::new(3e6),
            1e22,
            controls(),
            || false
        )
        .is_err()
    );
}

#[test]
fn g0_isobaric_rest_preserves_all_bits_and_invalid_heat_cannot_publish() {
    let mut body = parcel(source(), 2.0 * MIXTURE_MASS, 220.0);
    let before = body.clone();
    let mut capped = controls();
    capped.enthalpy.maximum_flash_evaluations = 1;
    let rest = body.add_heat(0.0, capped, || false).unwrap();
    assert_eq!(body, before);
    assert_eq!(rest.pt_evaluations, 0);
    assert_eq!(rest.eos_evaluation_upper_bound, 0);
    assert_eq!(rest.boundary_work_j, 0.0);
    assert_eq!(rest.energy_residual_j, 0.0);
    for heat in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, 1e308] {
        assert!(body.add_heat(heat, controls(), || false).is_err());
        assert_eq!(body, before);
    }
    assert!(body.add_heat(1.0, capped, || false).is_err());
    assert_eq!(body, before);
    assert_eq!(body.add_heat(0.0, controls(), || true), Err(E::Cancelled));
    assert_eq!(body, before);
}

#[test]
fn g4_isobaric_final_cancellation_preserves_pressure_volume_and_original_phase() {
    let initial = parcel(source(), 2.0 * MIXTURE_MASS, 220.0);
    let mut accepted = initial.clone();
    let mut polls = 0;
    let expected = accepted
        .add_heat(1.0, controls(), || {
            polls += 1;
            false
        })
        .unwrap();
    let mut cancelled = initial.clone();
    let mut seen = 0;
    assert_eq!(
        cancelled.add_heat(1.0, controls(), || {
            seen += 1;
            seen == polls
        }),
        Err(E::Cancelled)
    );
    assert_eq!(seen, polls);
    assert_eq!(cancelled, initial);
    assert_eq!(
        cancelled.add_heat(1.0, controls(), || false).unwrap(),
        expected
    );
    assert_eq!(cancelled, accepted);
}
