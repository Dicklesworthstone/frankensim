use super::*;
use crate::fusion::AVOGADRO;
use fs_qty::Amount;
use fs_thermochem::equilibrium::composition::IdealGasEquilibriumV1;
use fs_thermochem::fusion::{KELVIN_PER_KEV, NuclearSpeciesV1};

#[path = "../../../examples/dt_fusion/source.rs"]
mod dt_source;
#[path = "../network/tests/fixture.rs"]
mod network_source;

fn limits() -> DtFusionPlasmaLimitsV1 {
    DtFusionPlasmaLimitsV1 {
        maximum_coulomb_coupling: 0.01,
        maximum_thermal_phase_space_density: 0.01,
        maximum_thermal_energy_rest_energy_ratio: 0.05,
    }
}

fn close(actual: f64, expected: f64, relative: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * expected.abs(),
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}

#[test]
fn g1_both_fusion_banks_compress_and_expand_both_pools_adiabatically() {
    let (dt, dt_state) = dt_source::plasma([3e19, 7e19, 1e18], 0.5, [8.0, 2.0]);
    let (network, network_state) =
        network_source::plasma([2e19, 3e20, 7e19, 4e19, 1e19], 0.5, [8.0, 2.0]);
    for (plasma, initial) in [
        (FusionPlasmaRef::Dt(&dt), &dt_state),
        (FusionPlasmaRef::LightIon(&network), &network_state),
    ] {
        for volume_ratio in [0.125_f64, 8.0] {
            let result = advance_adiabatic_fusion_volume(
                initial,
                plasma,
                Volume::new(initial.volume().value() * volume_ratio),
                limits(),
                || false,
            )
            .unwrap();
            let factor = volume_ratio.powf(-2.0 / 3.0);
            close(
                result.state().heavy_temperature().value(),
                initial.heavy_temperature().value() * factor,
                3e-14,
                0.0,
            );
            close(
                result.state().electron_temperature().value(),
                initial.electron_temperature().value() * factor,
                3e-14,
                0.0,
            );
            close(
                result.state().pressure().value(),
                initial.pressure().value() * factor / volume_ratio,
                3e-14,
                0.0,
            );
            for (work, before, after) in [
                (
                    result.boundary_work_on_heavy(),
                    initial.heavy_internal_energy(),
                    result.state().heavy_internal_energy(),
                ),
                (
                    result.boundary_work_on_electrons(),
                    initial.electron_internal_energy(),
                    result.state().electron_internal_energy(),
                ),
            ] {
                close(work.value(), before.value() * (factor - 1.0), 3e-14, 0.0);
                close(after.value() - before.value(), work.value(), 3e-14, 0.0);
            }
            assert_eq!(result.state().amounts(), initial.amounts());
            close(
                result.energy_residual_j(),
                0.0,
                0.0,
                3e-13 * initial.internal_energy().value(),
            );
        }
    }
}

#[test]
fn g1_adjacent_volumes_retain_nonzero_pressure_work() {
    let (plasma, initial) = dt_source::plasma([3e19, 7e19, 0.0], 1.0, [8.0, 2.0]);
    let final_volume = f64::from_bits(1.0_f64.to_bits() + 1);
    let result = advance_adiabatic_fusion_volume(
        &initial,
        (&plasma).into(),
        Volume::new(final_volume),
        limits(),
        || false,
    )
    .unwrap();
    let fraction = ((-2.0 / 3.0) * (final_volume - 1.0).ln_1p()).exp_m1();
    assert!(result.boundary_work_on_heavy().value() < 0.0);
    assert!(result.boundary_work_on_electrons().value() < 0.0);
    close(
        result.boundary_work_on_heavy().value(),
        initial.heavy_internal_energy().value() * fraction,
        3e-15,
        0.0,
    );
    close(
        result.boundary_work_on_electrons().value(),
        initial.electron_internal_energy().value() * fraction,
        3e-15,
        0.0,
    );
    // A finite pair need not have a finite quotient; the logarithm still does.
    let tiny = Volume::new(f64::from_bits(1));
    let large = Volume::new(f64::MAX);
    close(
        log_volume_ratio(tiny, large).unwrap(),
        f64::MAX.ln() - f64::from_bits(1).ln(),
        2e-15,
        0.0,
    );
    close(
        log_volume_ratio(large, tiny).unwrap(),
        -log_volume_ratio(tiny, large).unwrap(),
        2e-15,
        0.0,
    );
}

#[test]
fn g0_volume_reversal_and_noop_preserve_original_sources_and_inventories() {
    let (plasma, initial) = network_source::plasma([2e19, 3e20, 7e19, 4e19, 1e19], 1.0, [8.0, 2.0]);
    let forward = advance_adiabatic_fusion_volume(
        &initial,
        (&plasma).into(),
        Volume::new(0.37),
        limits(),
        || false,
    )
    .unwrap();
    let reverse = advance_adiabatic_fusion_volume(
        forward.state(),
        (&plasma).into(),
        initial.volume(),
        limits(),
        || false,
    )
    .unwrap();
    assert_eq!(forward.initial_state(), &initial);
    assert_eq!(reverse.state().thermodynamics(), initial.thermodynamics());
    assert_eq!(reverse.state().amounts(), initial.amounts());
    for (actual, expected) in [
        (
            reverse.state().heavy_temperature(),
            initial.heavy_temperature(),
        ),
        (
            reverse.state().electron_temperature(),
            initial.electron_temperature(),
        ),
    ] {
        close(actual.value(), expected.value(), 3e-14, 0.0);
    }
    close(
        forward.total_boundary_work().value() + reverse.total_boundary_work().value(),
        0.0,
        0.0,
        3e-13 * initial.internal_energy().value(),
    );
    let noop = advance_adiabatic_fusion_volume(
        &initial,
        (&plasma).into(),
        initial.volume(),
        limits(),
        || false,
    )
    .unwrap();
    assert_eq!(noop.state(), &initial);
    assert_eq!(noop.state().amounts().as_ptr(), initial.amounts().as_ptr());
    assert_eq!(noop.total_boundary_work().value(), 0.0);
    assert_eq!(noop.energy_residual_j(), 0.0);
}

#[test]
fn g3_volume_motion_refuses_source_physical_domain_and_late_cancellation() {
    let (network, initial) = network_source::plasma_with_helium3_limit(
        [0.0, 3e20, 2e20, 0.0, 0.0],
        1.0,
        [10.0, 3.0],
        12.0 * KELVIN_PER_KEV,
    );
    assert!(matches!(
        advance_adiabatic_fusion_volume(
            &initial,
            (&network).into(),
            Volume::new(0.5),
            limits(),
            || false,
        ),
        Err(DtFusionError::Plasma(DtFusionPlasmaErrorV1::Source { .. }))
    ));
    let (plasma, initial) = dt_source::plasma([1e20, 1e20, 0.0], 1.0, [10.0, 5.0]);
    let original = initial.clone();
    let restrictive = DtFusionPlasmaLimitsV1 {
        maximum_thermal_energy_rest_energy_ratio: 0.02,
        ..limits()
    };
    assert!(matches!(
        advance_adiabatic_fusion_volume(
            &initial,
            (&plasma).into(),
            Volume::new(0.125),
            restrictive,
            || false,
        ),
        Err(DtFusionError::Plasma(DtFusionPlasmaErrorV1::Domain(_)))
    ));
    assert!(matches!(
        advance_adiabatic_fusion_volume(
            &initial,
            (&plasma).into(),
            Volume::new(0.0),
            limits(),
            || false,
        ),
        Err(DtFusionError::Invalid(_))
    ));
    let mut calls = 0;
    advance_adiabatic_fusion_volume(
        &initial,
        (&plasma).into(),
        Volume::new(0.5),
        limits(),
        || {
            calls += 1;
            false
        },
    )
    .unwrap();
    let mut cancelled_calls = 0;
    assert_eq!(
        advance_adiabatic_fusion_volume(
            &initial,
            (&plasma).into(),
            Volume::new(0.5),
            limits(),
            || {
                cancelled_calls += 1;
                cancelled_calls == calls
            },
        ),
        Err(DtFusionError::Cancelled)
    );
    assert_eq!(initial, original);
}

#[test]
fn g3_volume_motion_cannot_cross_a_gap_between_admitted_source_endpoints() {
    use fs_evidence::ValidityDomain;
    use fs_thermochem::{Nasa9RegionV1, Nasa9StandardStateModelV1};

    let (plasma, original) = dt_source::plasma([3e19, 7e19, 0.0], 1.0, [8.0, 2.0]);
    let bank = plasma.thermodynamics();
    let ids = plasma
        .indices()
        .map(|index| bank.models()[index].species().clone());
    let mut models = bank.models().to_vec();
    let source = &models[plasma.indices()[0]];
    let regions =
        [(1e4, 9.0 * KELVIN_PER_KEV), (11.0 * KELVIN_PER_KEV, 2e9)].map(|(lower, upper)| {
            let mut card = source.regions()[0].card().clone();
            card.validity = ValidityDomain::unconstrained().with("T", lower, upper);
            Nasa9RegionV1::from_card(card).unwrap()
        });
    models[plasma.indices()[0]] = Nasa9StandardStateModelV1::new(
        source.species().clone(),
        source.molar_mass(),
        source.convention().clone(),
        regions.into(),
    )
    .unwrap();
    let gapped = DtFusionPlasmaV1::new(
        IdealGasEquilibriumV1::new(
            bank.elemental_matrix().clone(),
            bank.charge_vector().clone(),
            models,
        )
        .unwrap(),
        ids[0].clone(),
        ids[1].clone(),
        ids[2].clone(),
    )
    .unwrap();
    let initial = FrozenIonizedGasStateV1::new(
        gapped.thermodynamics(),
        original.amounts(),
        original.volume(),
        original.heavy_temperature(),
        original.electron_temperature(),
        64.0 * f64::EPSILON,
        || false,
    )
    .unwrap();
    assert!(matches!(
        advance_adiabatic_fusion_volume(
            &initial, (&gapped).into(), Volume::new(0.125), limits(), || false,
        ),
        Err(DtFusionError::Plasma(DtFusionPlasmaErrorV1::Invalid(message)))
            if message.contains("gap")
    ));
}
