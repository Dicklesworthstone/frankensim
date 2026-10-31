//! A source-bound atomic ionization substep changes the original FEEC current.
//! NIST BEB supplies the event rate; the momentum cross section is explicitly
//! manufactured. The two powers calculated from the initial field are HELD
//! over the chemistry interval. The final fixed-voltage field is a subsequent
//! feedback evaluation, not the time integral of an evolving voltage source.
//! Elastic exchange is a separate owner; this test isolates chemical storage.

mod common;
#[path = "../../fs-material/tests/common/ionized_gas.rs"]
mod fixture;

use common::{config, ends, layered, limits};
use fs_em::conduction::{ConductionControls, PlanarConduction};
use fs_gas::chemistry::ionization::{AtomicIonizationSolve, advance_atomic_ionization};
use fs_material::ionized_gas::IonizedGasConductivity;
use fs_qty::{Power, Temperature, Time};
use fs_thermochem::electron_collision::{MaxwellianRateControlsV1, nist};
use fs_thermochem::equilibrium::composition::state::FrozenIonizedGasStateV1;
use fs_thermochem::kinetics::ElectronImpactRateV1;
use fs_thermochem::kinetics::atomic_ionization::{AtomicIonizationV1, AtomicStatisticalWeightsV1};
use fs_thermochem::{ReactionId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};

fn material(state: &FrozenIonizedGasStateV1) -> IonizedGasConductivity {
    let electron = fixture::electron_relaxation(state, 1e-22);
    let ion = fixture::source_two_temperature(
        "H+",
        state.heavy_temperature().value(),
        state.electron_temperature().value(),
        state.pressure().value(),
        1e9,
    );
    IonizedGasConductivity::from_state(
        state.clone(),
        vec![electron, ion],
        fixture::limits(),
        || false,
    )
    .unwrap()
}

fn field_controls() -> ConductionControls {
    ConductionControls {
        linear: config(),
        current_tolerance_a: 1e-6,
        power_tolerance_w: 1e-6,
    }
}

#[test]
fn g1_two_temperature_atomic_ionization_changes_current_and_closes_electrical_work() {
    let initial = fixture::frozen_state(1e-4, 6000.0, 16000.0, 0.2);
    let forward = ElectronImpactRateV1::new(
        fixture::id("e-"),
        nist::hydrogen(fixture::id("H")),
        MaxwellianRateControlsV1 {
            relative_tolerance: 1e-10,
            maximum_subdivisions: 4096,
        },
        0.1,
    )
    .unwrap();
    let channel = AtomicIonizationV1::new(
        initial.thermodynamics().clone(),
        ReactionId::new("source-bound-hydrogen-impact").unwrap(),
        fixture::id("H+"),
        forward,
        AtomicStatisticalWeightsV1 { neutral: 2, ion: 1 },
        1e-12,
    )
    .unwrap();
    let before = material(&initial);
    let sigma_before = before.conductivity_s_m();
    assert!(before.admits_electric_field([2.0, 0.0]));
    let original = PlanarConduction::new(
        layered(0.5, 1.0, 0.2),
        vec![before; 4],
        limits(),
        &mut || false,
    )
    .unwrap();
    let field = original
        .solve(&ends(2.0), field_controls(), &mut || false)
        .unwrap();
    let resistance_before = original
        .resistance(&[2, 5], &[0, 3], field_controls(), &mut || false)
        .unwrap()
        .resistance_ohm();
    let mut powers = [0.0, 0.0];
    let mut volume = 0.0;
    for (cell, e) in field.electric_field_v_m().iter().enumerate() {
        let work = original.conductors()[cell].joule_heating_w_m3(*e).unwrap();
        let cell_volume = original.mesh().face_measure(cell).unwrap();
        powers[0] += cell_volume * work.heavy_power_density_w_m3();
        powers[1] += cell_volume * work.electron_power_density_w_m3();
        volume += cell_volume;
        fixture::relative(
            cell_volume * work.total_power_density_w_m3(),
            field.cell_power_w()[cell],
            1e-13,
        );
    }
    fixture::relative(volume, initial.volume().value(), 1e-14);
    fixture::relative(powers[0] + powers[1], field.power_w(), 1e-13);
    let mut controls =
        AtomicIonizationSolve::new(Temperature::new(2000.0), Temperature::new(20000.0));
    controls.relative_tolerance = 1e-7;
    controls.energy_tolerance_j = 1e-8;
    controls.initial_step_s = 1e-10;
    controls.maximum_step_s = 1e-9;
    let duration = 1e-8;
    let result = advance_atomic_ionization(
        &initial,
        &channel,
        Power::new(powers[1]),
        Power::new(powers[0]),
        Time::new(duration),
        controls,
        || false,
    )
    .unwrap();
    let next = result.state();
    let n0 = initial.amounts();
    let n1 = next.amounts();
    let extent = result.reaction_extent().value();
    assert!(
        extent > 1e-8,
        "the source must produce observable electron-impact growth"
    );
    fixture::relative(n1[1].value() - n0[1].value(), extent, 1e-10);
    fixture::relative(n1[2].value() - n0[2].value(), extent, 1e-10);
    fixture::relative(n1[0].value() + n1[1].value(), 1.0, 1e-12);
    fixture::relative(n1[1].value(), n1[2].value(), 1e-12);
    fixture::relative(next.mass().value(), initial.mass().value(), 1e-12);
    assert_eq!(next.volume(), initial.volume());
    assert!(next.electron_temperature().value() < initial.electron_temperature().value());
    assert_ne!(next.electron_temperature(), next.heavy_temperature());

    // Independent monatomic energy account: electron number changes as well as
    // Te. The source ionization gap enters chemical storage once, without an
    // additional subtraction of the NIST cross-section threshold.
    let electron_sensible = 1.5
        * R
        * (n1[2].value() * next.electron_temperature().value()
            - n0[2].value() * initial.electron_temperature().value());
    let heavy_sensible =
        1.5 * R * (next.heavy_temperature().value() - initial.heavy_temperature().value());
    let chemical_storage = extent * channel.ionization_energy_j_mol();
    let supplied = (powers[0] + powers[1]) * duration;
    assert!((electron_sensible + chemical_storage - powers[1] * duration).abs() < 2e-8);
    assert!((heavy_sensible - powers[0] * duration).abs() < 2e-8);
    assert!((electron_sensible + heavy_sensible + chemical_storage - supplied).abs() < 3e-8);
    let balance = result.balance();
    assert!((balance.electron_deposited_energy_j - powers[1] * duration).abs() < 1e-15);
    assert!((balance.heavy_deposited_energy_j - powers[0] * duration).abs() < 1e-15);
    fixture::relative(balance.chemical_storage_change_j, chemical_storage, 1e-12);
    assert!(balance.energy_residual_j.abs() <= controls.energy_tolerance_j);

    let after = material(next);
    assert!(after.equilibrium().is_none());
    assert_eq!(after.state().amounts(), n1);
    assert!(after.conductivity_s_m() > sigma_before);
    assert!(after.admits_electric_field([2.0, 0.0]));
    let sigma_after = after.conductivity_s_m();
    let updated = PlanarConduction::new(
        layered(0.5, 1.0, 0.2),
        vec![after; 4],
        limits(),
        &mut || false,
    )
    .unwrap();
    let new_field = updated
        .solve(&ends(2.0), field_controls(), &mut || false)
        .unwrap();
    fixture::relative(
        new_field.current_density_a_m2()[0][0],
        -2.0 * sigma_after,
        1e-12,
    );
    fixture::relative(new_field.power_w(), 0.8 * sigma_after, 1e-12);
    assert!(new_field.power_w() > field.power_w());
    let resistance_after = updated
        .resistance(&[2, 5], &[0, 3], field_controls(), &mut || false)
        .unwrap()
        .resistance_ohm();
    assert!(resistance_after < resistance_before);
}
