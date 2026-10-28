//! Independent inventories/energy and retained state on changing atmospheres.
//! These are synthetic thermodynamic/solid fixtures, not instrument validation.
use super::*;
use crate::acoustic_realize::finite_gas::flow::FiniteGasStringFlowDirection;
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_material::thermochemical_gas::energy::parcel::transfer::GasMassTransferSolve;
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::{
    GasNozzle, GasNozzleStepSolve,
};

fn nozzle(area_m2: f64) -> GasNozzle {
    GasNozzle {
        area_m2,
        discharge_coefficient: 0.8,
    }
}
fn flow_controls() -> GasNozzleStepSolve {
    GasNozzleStepSolve {
        nozzle: limits().gas.gas,
        transfer: GasMassTransferSolve {
            donor: limits().gas.gas,
            receiver: limits().gas.gas,
            mass_tolerance_kg: 1e-14,
            energy_tolerance_j: 1e-8,
        },
        maximum_donor_fraction: 0.1,
    }
}
fn exchange(
    _: &fs_material::thermochemical_gas::ResolvedThermochemicalGasState,
    _: fs_material::phase::EquilibriumPhaseState,
) -> Result<ConductanceSample<()>, &'static str> {
    Ok(ConductanceSample {
        conductance_w_k: 0.1,
        source: (),
    })
}
fn inventory(parcel: &ClosedGasParcel) -> BTreeMap<fs_thermochem::SpeciesId, f64> {
    let model = parcel.model().thermodynamics();
    let fractions = model.mole_composition().fractions();
    let molar_mass: f64 = model
        .components()
        .iter()
        .zip(fractions)
        .map(|(species, x)| x * species.molar_mass().value())
        .sum();
    model
        .components()
        .iter()
        .zip(fractions)
        .map(|(species, x)| {
            (
                species.species().clone(),
                parcel.mass_kg() * x * species.molar_mass().value() / molar_mass,
            )
        })
        .collect()
}
fn storage(parcel: &ClosedGasParcel) -> f64 {
    parcel.mass_kg() * parcel.state().specific_internal_energy_j_kg()
}

#[test]
fn gas_injection_conserves_each_species_and_all_storage_without_resetting_layer_or_observer_history()
 {
    with_cx(|cx, _| {
        let sources = Sources::new(4e-5);
        let mut gas = parcel(0.25, 100000.0, 330.0);
        let (mut live, mut thermal) = runtime(cx, &sources, gas.state().state());
        for _ in 0..8 {
            live.step_with_isochoric_gas_exchange(
                cx,
                &mut thermal,
                &mut gas,
                limits(),
                &[0.0; 2],
                &[0.0; 2],
                0.0,
                None,
                contact(),
                no_heat,
                exchange,
            )
            .unwrap();
        }
        let mut reservoir = ClosedGasParcel::from_isochoric(parcel(0.8, 200000.0, 340.0));
        let original = live.clone();
        let old_local = ClosedGasParcel::from_isochoric(gas.clone());
        let old_reservoir = reservoir.clone();
        let old_thermal = thermal.clone();
        let previous = original
            .member
            .observer
            .previous_acceleration
            .clone()
            .unwrap();
        let mut heat_only = live.clone();
        let mut heat_only_cells = thermal.clone();
        let mut heat_only_gas = gas.clone();
        let gas_heat_j = 0.002;
        let heated = heat_only
            .step_with_isochoric_gas_exchange(
                cx,
                &mut heat_only_cells,
                &mut heat_only_gas,
                limits(),
                &[0.002, 0.006],
                &[0.0; 2],
                gas_heat_j,
                Some(OhmicDrive::Voltage(2.0)),
                contact(),
                supplied_heat,
                exchange,
            )
            .unwrap();
        let (expected, native_flow) = reservoir
            .through_nozzle_to(
                &ClosedGasParcel::from_isochoric(heat_only_gas.clone()),
                nozzle(1e-5),
                1.0 / f64::from(RATE),
                flow_controls(),
                || false,
            )
            .unwrap();
        let (expected_reservoir, expected_gas, expected_transfer) = expected.into_parts();
        let frame = live
            .step_with_isochoric_gas_nozzle(
                cx,
                &mut thermal,
                &mut gas,
                &mut reservoir,
                limits(),
                &[0.002, 0.006],
                &[0.0; 2],
                gas_heat_j,
                Some(OhmicDrive::Voltage(2.0)),
                contact(),
                nozzle(1e-5),
                flow_controls(),
                supplied_heat,
                exchange,
            )
            .unwrap();
        assert_eq!(
            frame.direction,
            FiniteGasStringFlowDirection::IntoSurroundings
        );
        assert_eq!(frame.flow, native_flow);
        assert_eq!(frame.transfer, expected_transfer);
        assert_eq!(reservoir, expected_reservoir);
        assert_eq!(gas, expected_gas.into_isochoric());
        assert_eq!(frame.thermal.string.source.next_gas, heat_only_gas);
        assert_eq!(frame.thermal.string.source.gas, heated.string.source.gas);
        assert_eq!(
            frame.thermal.solid_to_gas_heat_j,
            heated.solid_to_gas_heat_j
        );
        assert_eq!(thermal.cells(), heat_only_cells.cells());
        assert_eq!(thermal.time_s(), heat_only_cells.time_s());
        assert_eq!(live.states(), heat_only.states());
        assert_eq!(live.source(), heat_only.source());
        assert_eq!(live.accepted_samples(), 9);
        assert_eq!(live.gas(), gas.state().state());
        assert_ne!(live.modes(), heat_only.modes());
        assert!(frame.transfer.transferred_mass_kg > 0.0);
        assert_eq!(gas.point().axes().get("material-lot"), Some(&7.0));

        let final_local = ClosedGasParcel::from_isochoric(gas.clone());
        let old_a = inventory(&old_local);
        let old_b = inventory(&old_reservoir);
        let new_a = inventory(&final_local);
        let new_b = inventory(&reservoir);
        for species in old_a.keys() {
            assert!(
                (new_a[species] + new_b[species] - old_a[species] - old_b[species]).abs()
                    <= flow_controls().transfer.mass_tolerance_kg
            );
        }
        let local_change = storage(&final_local) - storage(&old_local);
        let reservoir_change = storage(&reservoir) - storage(&old_reservoir);
        let loss = &frame.thermal.string.vibration.dissipation;
        let supplied = frame.thermal.string.vibration.acoustic.input_work_j
            + frame.thermal.external_solid_heat_j
            + gas_heat_j;
        let defect = supplied
            - (independent_internal(&live, 4e-5) - independent_internal(&original, 4e-5)
                + thermal.energy_j()
                - old_thermal.energy_j()
                + local_change
                + reservoir_change
                + loss.authored_loss_j
                + loss.roundoff_residual_j);
        assert!(defect.abs() <= frame.energy_tolerance_j);
        assert!(frame.energy_balance_residual_j.abs() <= frame.energy_tolerance_j);
        assert!((frame.surrounding_energy_change_j - local_change).abs() < 1e-10);
        assert!((frame.reservoir_energy_change_j - reservoir_change).abs() < 1e-10);

        // The native compact dipole at a 1 m observer, independently assembled
        // from the two modal accelerations and the PREVIOUS accepted sample.
        let a = linear_accelerations(&live.member.model, &[0.0; 2]);
        let medium = live.gas();
        let string = live.string();
        let pi = core::f64::consts::PI;
        let scale = medium.density * string.width_m
            / (4.0 * pi * (string.lin_density_kg_m * string.length_m / 2.0).sqrt());
        let pressure = scale * 2.0 * string.length_m / pi * a[0]
            - scale * string.length_m.powi(2) / (2.0 * pi * medium.sound_speed)
                * (a[1] - previous[1])
                * f64::from(RATE);
        near(
            frame.thermal.string.vibration.acoustic.observer_pressure_pa,
            pressure,
            2e-10,
        );
        assert_ne!(
            frame.thermal.string.vibration.acoustic.observer_pressure_pa,
            heated.string.vibration.acoustic.observer_pressure_pa
        );
    });
}

#[test]
fn injected_composition_alters_later_motion_and_pressure_with_the_same_layer_data() {
    with_cx(|cx, _| {
        let sources = Sources::new(4e-5);
        let mut gas_a = parcel(0.25, 100000.0, 300.0);
        let mut gas_b = gas_a.clone();
        let (mut a, mut thermal_a) = runtime(cx, &sources, gas_a.state().state());
        let (mut b, mut thermal_b) = (a.clone(), thermal_a.clone());
        let mut reservoir_a = ClosedGasParcel::from_isochoric(parcel(0.8, 200000.0, 300.0));
        let mut reservoir_b = ClosedGasParcel::from_isochoric(parcel(0.4, 200000.0, 300.0));
        let mut pressure_a = Vec::new();
        let mut pressure_b = Vec::new();
        for sample in 0..16 {
            let left = a
                .step_with_isochoric_gas_nozzle(
                    cx,
                    &mut thermal_a,
                    &mut gas_a,
                    &mut reservoir_a,
                    limits(),
                    &[0.0; 2],
                    &[0.0; 2],
                    0.0,
                    None,
                    contact(),
                    nozzle(1e-5),
                    flow_controls(),
                    no_heat,
                    exchange,
                )
                .unwrap();
            let right = b
                .step_with_isochoric_gas_nozzle(
                    cx,
                    &mut thermal_b,
                    &mut gas_b,
                    &mut reservoir_b,
                    limits(),
                    &[0.0; 2],
                    &[0.0; 2],
                    0.0,
                    None,
                    contact(),
                    nozzle(1e-5),
                    flow_controls(),
                    no_heat,
                    exchange,
                )
                .unwrap();
            if sample == 0 {
                assert_eq!(a.states(), b.states());
                assert_eq!(thermal_a.cells(), thermal_b.cells());
                assert_eq!(a.source(), b.source());
            }
            assert_eq!(a.accepted_samples(), sample + 1);
            assert_eq!(b.accepted_samples(), sample + 1);
            assert_eq!(a.gas(), gas_a.state().state());
            assert_eq!(b.gas(), gas_b.state().state());
            assert!(left.energy_balance_residual_j.abs() <= left.energy_tolerance_j);
            assert!(right.energy_balance_residual_j.abs() <= right.energy_tolerance_j);
            pressure_a.push(left.thermal.string.vibration.acoustic.observer_pressure_pa);
            pressure_b.push(right.thermal.string.vibration.acoustic.observer_pressure_pa);
        }
        assert_ne!(a.gas().sound_speed, b.gas().sound_speed);
        assert_ne!(a.modes(), b.modes());
        assert_ne!(a.states(), b.states());
        assert_ne!(pressure_a, pressure_b);
    });
}

#[test]
fn closed_nozzle_is_the_original_heat_step_and_outflow_retains_direction_and_inventory() {
    with_cx(|cx, _| {
        let sources = Sources::new(4e-5);
        for area in [0.0, 1e-5] {
            let mut gas = parcel(0.25, 200000.0, 330.0);
            let (mut live, mut thermal) = runtime(cx, &sources, gas.state().state());
            let mut reservoir = ClosedGasParcel::from_isochoric(parcel(0.8, 100000.0, 330.0));
            let old_reservoir = reservoir.clone();
            let old_mass = gas.mass_kg();
            let (mut control, mut cells, mut control_gas) =
                (live.clone(), thermal.clone(), gas.clone());
            let expected = control
                .step_with_isochoric_gas_exchange(
                    cx,
                    &mut cells,
                    &mut control_gas,
                    limits(),
                    &[0.0; 2],
                    &[0.0; 2],
                    0.0,
                    None,
                    contact(),
                    no_heat,
                    exchange,
                )
                .unwrap();
            let frame = live
                .step_with_isochoric_gas_nozzle(
                    cx,
                    &mut thermal,
                    &mut gas,
                    &mut reservoir,
                    limits(),
                    &[0.0; 2],
                    &[0.0; 2],
                    0.0,
                    None,
                    contact(),
                    nozzle(area),
                    flow_controls(),
                    no_heat,
                    exchange,
                )
                .unwrap();
            assert_eq!(
                frame.direction,
                FiniteGasStringFlowDirection::OutOfSurroundings
            );
            assert_eq!(live.states(), control.states());
            assert_eq!(thermal.cells(), cells.cells());
            if area == 0.0 {
                assert_eq!(gas, control_gas);
                assert_eq!(reservoir, old_reservoir);
                assert_eq!(frame.transfer.transferred_mass_kg, 0.0);
                assert_eq!(frame.thermal.string.vibration, expected.string.vibration);
                assert_eq!(
                    live.member.observer.previous_acceleration,
                    control.member.observer.previous_acceleration
                );
            } else {
                assert!(gas.mass_kg() < old_mass);
                assert!(reservoir.mass_kg() > old_reservoir.mass_kg());
                assert!(frame.transfer.transferred_mass_kg > 0.0);
                assert!(frame.energy_balance_residual_j.abs() <= frame.energy_tolerance_j);
            }
        }
    });
}

#[test]
fn late_flow_refusals_and_cancellation_preserve_all_four_owners_and_exact_retry() {
    with_cx(|cx, gate| {
        let sources = Sources::new(4e-5);
        let gas = parcel(0.25, 100000.0, 300.0);
        let (original, original_thermal) = runtime(cx, &sources, gas.state().state());
        let reservoir = ClosedGasParcel::from_isochoric(parcel(0.8, 200000.0, 340.0));
        for case in 0..3 {
            let (mut live, mut thermal, mut local, mut other) = (
                original.clone(),
                original_thermal.clone(),
                gas.clone(),
                reservoir.clone(),
            );
            let mut controls = limits();
            let mut flow = flow_controls();
            if case == 0 {
                flow.maximum_donor_fraction = 1e-14;
            }
            if case == 1 {
                controls.gas.maximum_gas_temperature_increment_k = 1e-3;
            }
            let evaluated = Cell::new(false);
            let refused = live.step_with_isochoric_gas_nozzle(
                cx,
                &mut thermal,
                &mut local,
                &mut other,
                controls,
                &[0.002, 0.006],
                &[0.0; 2],
                0.0,
                Some(OhmicDrive::Voltage(2.0)),
                contact(),
                nozzle(1e-5),
                flow,
                supplied_heat,
                |a, b| {
                    evaluated.set(true);
                    if case == 2 {
                        gate.request();
                    }
                    exchange(a, b)
                },
            );
            assert!(evaluated.get() && refused.is_err());
            assert_eq!(local, gas);
            assert_eq!(other, reservoir);
            assert_eq!(thermal.cells(), original_thermal.cells());
            assert_eq!(thermal.time_s(), original_thermal.time_s());
            assert_eq!(live.states(), original.states());
            assert_eq!(live.source(), original.source());
            assert_eq!(live.gas(), original.gas());
            assert_eq!(live.accepted_samples(), original.accepted_samples());
            with_cx(|retry_cx, _| {
                let (mut baseline, mut baseline_cells, mut baseline_gas, mut baseline_reservoir) = (
                    original.clone(),
                    original_thermal.clone(),
                    gas.clone(),
                    reservoir.clone(),
                );
                let actual = live
                    .step_with_isochoric_gas_nozzle(
                        retry_cx,
                        &mut thermal,
                        &mut local,
                        &mut other,
                        limits(),
                        &[0.002, 0.006],
                        &[0.0; 2],
                        0.0,
                        Some(OhmicDrive::Voltage(2.0)),
                        contact(),
                        nozzle(1e-5),
                        flow_controls(),
                        supplied_heat,
                        exchange,
                    )
                    .unwrap();
                let expected = baseline
                    .step_with_isochoric_gas_nozzle(
                        retry_cx,
                        &mut baseline_cells,
                        &mut baseline_gas,
                        &mut baseline_reservoir,
                        limits(),
                        &[0.002, 0.006],
                        &[0.0; 2],
                        0.0,
                        Some(OhmicDrive::Voltage(2.0)),
                        contact(),
                        nozzle(1e-5),
                        flow_controls(),
                        supplied_heat,
                        exchange,
                    )
                    .unwrap();
                assert_eq!(
                    actual.thermal.string.vibration,
                    expected.thermal.string.vibration
                );
                assert_eq!(actual.transfer, expected.transfer);
                assert_eq!(
                    actual.energy_balance_residual_j,
                    expected.energy_balance_residual_j
                );
                assert_eq!(local, baseline_gas);
                assert_eq!(other, baseline_reservoir);
                assert_eq!(thermal.cells(), baseline_cells.cells());
                assert_eq!(live.states(), baseline.states());
            });
        }
    });
}
