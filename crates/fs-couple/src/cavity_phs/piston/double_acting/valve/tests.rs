//! Native parcel and shared-wall integration with independent physical accounts.
use super::*;
use crate::cavity_phs::piston::tests::{controls as single_controls, drive, gas};
use crate::cavity_phs::piston::{GasPiston, PistonSolve};
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_material::thermochemical_gas::energy::parcel::transfer::GasMassTransferSolve;
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

pub(in crate::cavity_phs::piston::double_acting) fn controls() -> DoublePistonSolve {
    DoublePistonSolve {
        piston: single_controls(),
        opposing_gas: single_controls().gas,
        volume_tolerance_m3: 1e-18,
    }
}
pub(in crate::cavity_phs::piston::double_acting) fn exchange() -> GasNozzleStepSolve {
    GasNozzleStepSolve {
        nozzle: single_controls().gas,
        transfer: GasMassTransferSolve {
            donor: single_controls().gas,
            receiver: single_controls().gas,
            mass_tolerance_kg: 1e-14,
            energy_tolerance_j: 1e-8,
        },
        maximum_donor_fraction: 0.1,
    }
}
pub(in crate::cavity_phs::piston::double_acting) fn input() -> DoublePistonInput {
    DoublePistonInput {
        duration_s: drive(0.0).duration_s,
        heat_j: [0.0; 2],
        external_force_n: 0.0,
    }
}
pub(in crate::cavity_phs::piston::double_acting) fn geometry() -> RectangularMeteringValve {
    RectangularMeteringValve {
        reference_primary_volume_m3: 0.001,
        edge_height_at_reference_m: 0.001,
        width_m: 0.01,
        window_height_m: 0.002,
        opens_on_primary_expansion: true,
        discharge_coefficient: 0.8,
    }
}
pub(in crate::cavity_phs::piston::double_acting) fn pair(heat: [f64; 2], slopes: [f64; 2], velocity: f64) -> DoubleActingGasPiston {
    let mut parcels = [gas(0.75, slopes[0]), gas(0.25, slopes[1])];
    for (parcel, heat) in parcels.iter_mut().zip(heat) {
        if heat != 0.0 {
            parcel
                .add_heat(heat, single_controls().gas, || false)
                .unwrap();
        }
    }
    let [primary, opposing] = parcels;
    DoubleActingGasPiston::new(
        GasPiston::new(primary, 1.0, 0.01, velocity).unwrap(),
        opposing,
        0.01,
    )
    .unwrap()
}
pub(in crate::cavity_phs::piston::double_acting) fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.16e} != {expected:.16e}"
    );
}
pub(in crate::cavity_phs::piston::double_acting) fn assert_unchanged(actual: &DoubleActingGasPiston, expected: &DoubleActingGasPiston) {
    assert_eq!(actual, expected);
    // The owner intentionally excludes accounting history from physical PartialEq.
    assert_eq!(actual.ports, expected.ports);
    assert_eq!(actual.primary().ports, expected.primary().ports);
}
pub(in crate::cavity_phs::piston::double_acting) fn species_masses(parcel: &ClosedGasParcel) -> [f64; 2] {
    let x = parcel
        .model()
        .thermodynamics()
        .mole_composition()
        .fractions();
    let molar_mass = x[0] * 0.028 + x[1] * 0.004;
    [
        parcel.mass_kg() * x[0] * 0.028 / molar_mass,
        parcel.mass_kg() * x[1] * 0.004 / molar_mass,
    ]
}
// Integrate the fixture's Cp/R polynomials directly, including both species.
// No reported energy change or production caloric inversion enters this reference.
pub(in crate::cavity_phs::piston::double_acting) fn internal_energy(parcel: &ClosedGasParcel, slope: f64) -> f64 {
    let x = parcel
        .model()
        .thermodynamics()
        .mole_composition()
        .fractions();
    let t = parcel.state().state().temperature;
    let molar_mass = x[0] * 0.028 + x[1] * 0.004;
    parcel.mass_kg() / molar_mass * R * (x[0] * (2.5 * t + 0.5 * slope * t * t) + x[1] * 1.5 * t)
}

#[test]
fn g1_both_upwind_directions_preserve_species_and_close_independent_enclosing_energy() {
    for (heat, donor, direction) in [
        ([100.0, 0.0], 0, ValveDirection::PrimaryToOpposing),
        ([0.0, 100.0], 1, ValveDirection::OpposingToPrimary),
    ] {
        let original = pair(heat, [0.001; 2], 0.03);
        let forcing = DoublePistonInput {
            heat_j: [0.2, 0.1],
            external_force_n: 3.0,
            ..input()
        };
        let proposed = original
            .valve_candidate(geometry(), forcing, exchange(), controls(), 1e-7, || false)
            .unwrap();
        let mut live = original.clone();
        let report = live
            .step_with_valve(geometry(), forcing, exchange(), controls(), 1e-7, || false)
            .unwrap();
        assert_eq!(report, *proposed.step());
        assert_unchanged(&live, proposed.piston());
        assert_eq!(report.direction, Some(direction));
        let before = [original.primary().gas(), original.opposing()];
        let after = [live.primary().gas(), live.opposing()];
        let transferred = report.transfer.transferred_mass_kg;
        assert!(transferred > 0.0);
        near(
            transferred,
            report.flow.mass_flow_kg_s * forcing.duration_s,
            1e-18,
        );
        near(
            after[donor].mass_kg() - before[donor].mass_kg(),
            -transferred,
            1e-14,
        );
        near(
            after[1 - donor].mass_kg() - before[1 - donor].mass_kg(),
            transferred,
            1e-14,
        );
        let donor_species = species_masses(before[donor]);
        for species in 0..2 {
            let expected = transferred * donor_species[species] / before[donor].mass_kg();
            near(report.transfer.species[species].mass_kg, expected, 1e-14);
            near(
                species_masses(after[donor])[species] - donor_species[species],
                -expected,
                1e-14,
            );
            near(
                species_masses(after[1 - donor])[species]
                    - species_masses(before[1 - donor])[species],
                expected,
                1e-14,
            );
        }
        let mut change = 0.0;
        for i in 0..2 {
            let delta = internal_energy(after[i], 0.001) - internal_energy(before[i], 0.001);
            near(report.gas_energy_change_j[i], delta, 1e-8);
            change += delta;
            assert_eq!(after[i].point().axes().get("specimen-lot"), Some(&42.0));
            assert_eq!(after[i].model().thermodynamics().components().len(), 2);
        }
        let kinetic = 0.5
            * live.primary().mass_kg()
            * (live.primary().velocity_m_s().powi(2) - original.primary().velocity_m_s().powi(2));
        let residual = change + kinetic
            - forcing.heat_j.iter().sum::<f64>()
            - forcing.external_force_n * report.piston.displacement_m;
        near(residual, 0.0, 2e-7);
        near(report.total_energy_residual_j, residual, 1e-8);
        near(report.total_mass_residual_kg, 0.0, 1e-14);
        assert_ne!(
            after[1 - donor].state().state().specific_gas_constant,
            before[1 - donor].state().state().specific_gas_constant
        );
        assert_eq!(report.piston.gas[0].supplied_heat_j, forcing.heat_j[0]);
        assert_eq!(report.piston.gas[1].supplied_heat_j, forcing.heat_j[1]);
        near(
            live.primary().gas().volume_m3() + live.opposing().volume_m3(),
            0.002,
            1e-18,
        );
    }
}

#[test]
fn g3_actual_wall_geometry_opens_and_closes_the_bypass_and_zero_flow_is_the_original_step() {
    for closed in [
        RectangularMeteringValve {
            edge_height_at_reference_m: -0.01,
            ..geometry()
        },
        RectangularMeteringValve {
            discharge_coefficient: 0.0,
            ..geometry()
        },
    ] {
        let mut live = pair([100.0, 0.0], [0.001; 2], 0.0);
        let mut reference = live.clone();
        let native = reference.step(input(), controls(), || false).unwrap();
        let report = live
            .step_with_valve(closed, input(), exchange(), controls(), 1e-7, || false)
            .unwrap();
        assert_unchanged(&live, &reference);
        assert_eq!(report.piston, native);
        assert_eq!(report.direction, None);
        assert_eq!(report.transfer.transferred_mass_kg, 0.0);
    }
    let mut live = pair([0.0; 2], [0.001; 2], 0.0);
    let opening = RectangularMeteringValve {
        edge_height_at_reference_m: 0.0,
        ..geometry()
    };
    let forcing = DoublePistonInput {
        external_force_n: 1000.0,
        ..input()
    };
    let first = live
        .step_with_valve(opening, forcing, exchange(), controls(), 1e-7, || false)
        .unwrap();
    assert_eq!(first.opening.area_m2, 0.0);
    assert_eq!(first.transfer.transferred_mass_kg, 0.0);
    assert!(first.endpoint_opening.open_height_m > 0.0);
    near(
        first.endpoint_opening.open_height_m,
        first.piston.displacement_m,
        1e-15,
    );
    near(
        first.endpoint_opening.area_m2,
        opening.width_m * first.piston.displacement_m,
        1e-17,
    );
    let second = live
        .step_with_valve(opening, forcing, exchange(), controls(), 1e-7, || false)
        .unwrap();
    assert_eq!(second.opening, first.endpoint_opening);
    assert_eq!(second.direction, Some(ValveDirection::OpposingToPrimary));
    assert!(second.transfer.transferred_mass_kg > 0.0);

    let mut live = pair([100.0, 0.0], [0.001; 2], 0.2);
    let closing = RectangularMeteringValve {
        edge_height_at_reference_m: 5e-6,
        window_height_m: 8e-6,
        opens_on_primary_expansion: false,
        ..geometry()
    };
    let first = live
        .step_with_valve(closing, input(), exchange(), controls(), 1e-7, || false)
        .unwrap();
    assert!(first.transfer.transferred_mass_kg > 0.0);
    assert!(first.endpoint_opening.edge_height_m < 0.0);
    assert_eq!(first.endpoint_opening.area_m2, 0.0);
    let second = live
        .step_with_valve(closing, input(), exchange(), controls(), 1e-7, || false)
        .unwrap();
    assert_eq!(second.direction, None);
    assert_eq!(second.transfer.transferred_mass_kg, 0.0);
}

#[test]
fn g1_constant_caloric_nozzle_matches_analytic_choking_and_changes_shared_wall_acceleration() {
    for heat in [100.0, 400.0] {
        let original = pair([heat, 0.0], [0.0; 2], 0.0);
        let upstream = original.primary().gas().state().state();
        let downstream_pressure = original.opposing().state().state().pressure;
        let gamma = 3.25_f64 / 2.25;
        let gas_constant = R / (0.75 * 0.028 + 0.25 * 0.004);
        let critical = (2.0 / (gamma + 1.0)).powf(gamma / (gamma - 1.0));
        let ratio = (downstream_pressure / upstream.pressure).max(critical);
        let expected_rate = geometry().discharge_coefficient
            * (geometry().width_m * geometry().edge_height_at_reference_m)
            * upstream.pressure
            / (gas_constant * upstream.temperature).sqrt()
            * (2.0 * gamma / (gamma - 1.0)
                * (ratio.powf(2.0 / gamma) - ratio.powf((gamma + 1.0) / gamma)))
            .sqrt();
        let mut live = original.clone();
        let report = live
            .step_with_valve(geometry(), input(), exchange(), controls(), 1e-7, || false)
            .unwrap();
        near(
            report.flow.mass_flow_kg_s,
            expected_rate,
            5e-9 * expected_rate,
        );
        assert_eq!(
            report.flow.choked,
            downstream_pressure / upstream.pressure <= critical
        );
        let mut closed = original;
        closed.step(input(), controls(), || false).unwrap();
        assert!(live.primary().velocity_m_s() > 0.0);
        assert!(closed.primary().velocity_m_s() - live.primary().velocity_m_s() > 1e-7);
        assert!(
            live.primary().gas().state().state().pressure
                < closed.primary().gas().state().state().pressure
        );
        assert!(
            live.opposing().state().state().pressure > closed.opposing().state().state().pressure
        );
    }
}

#[test]
fn g4_source_inventory_late_pressure_and_final_cancellation_refusals_preserve_exact_retry() {
    let mut source_conflict = pair([100.0, 0.0], [0.001, 0.0], 0.0);
    let original_conflict = source_conflict.clone();
    assert!(
        source_conflict
            .step_with_valve(geometry(), input(), exchange(), controls(), 1e-7, || false,)
            .is_err()
    );
    assert_unchanged(&source_conflict, &original_conflict);

    let mut original = pair([100.0, 0.0], [0.001; 2], 0.03);
    original.step(input(), controls(), || false).unwrap();
    let mut reference = original.clone();
    let mut polls = 0;
    let expected = reference
        .step_with_valve(geometry(), input(), exchange(), controls(), 1e-7, || {
            polls += 1;
            false
        })
        .unwrap();
    for refusal in 0..4 {
        let mut live = original.clone();
        let mut flow = exchange();
        let mut numerical = controls();
        let mut forcing = input();
        if refusal == 0 {
            flow.maximum_donor_fraction = 1e-12;
        }
        if refusal == 1 {
            let volume = original.primary().gas().volume_m3();
            numerical.piston = PistonSolve {
                volume_bracket_m3: [volume + 1e-5, volume + 2e-5],
                ..numerical.piston
            };
        }
        if refusal == 2 {
            forcing.heat_j[1] = 300.0;
            let mut mechanical = original.clone();
            mechanical.step(forcing, controls(), || false).unwrap();
            assert!(
                mechanical.opposing().state().state().pressure
                    > mechanical.primary().gas().state().state().pressure
            );
        }
        let mut count = 0;
        let result = live.step_with_valve(geometry(), forcing, flow, numerical, 1e-7, || {
            count += 1;
            refusal == 3 && count == polls
        });
        assert!(result.is_err());
        if refusal == 3 {
            assert_eq!(count, polls);
            assert_eq!(result, Err(PistonError::Cancelled));
        }
        assert_unchanged(&live, &original);
        let retry = live
            .step_with_valve(geometry(), input(), exchange(), controls(), 1e-7, || false)
            .unwrap();
        assert_eq!(retry, expected);
        assert_unchanged(&live, &reference);
        let mut control = reference.clone();
        for _ in 0..3 {
            assert_eq!(
                live.step_with_valve(geometry(), input(), exchange(), controls(), 1e-7, || false)
                    .unwrap(),
                control
                    .step_with_valve(geometry(), input(), exchange(), controls(), 1e-7, || false)
                    .unwrap()
            );
        }
        assert_unchanged(&live, &control);
    }
}
