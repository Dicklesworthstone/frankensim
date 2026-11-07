//! Coupled-state and transactional checks through the original chamber owner.
use super::*;
use super::super::{
    GasSurfaceError, GasSupportedChamberSolve, advance_resistive_chamber_with_films,
    growing::{
        GasGrowingCoatingSolve, GasIntegratedGrowingCoating,
        advance_growing_integrated_chamber_resolved,
    },
};
use crate::sublimation::gas::{
    GasSurfaceBinding, chamber::GasChamberSolve, tests as fixture,
};
use crate::thermal::conduction::{ConductivityCurve, LayerConductionConfig};
use fs_thermal::sublimation::{
    SolidSurfaceState, SolverConfig, SublimationLaw,
    coating::resistance::{ResistiveSurfaceState, geometry::ConductiveCoatingGeometry},
};
use std::cell::Cell;

fn controls() -> GasResistiveAdvanceSolve {
    GasResistiveAdvanceSolve {
        step: GasSupportedChamberSolve {
            chamber: GasChamberSolve {
                surface: fixture::controls(),
                maximum_surfaces: 8,
            },
            storage: SolverConfig {
                relative_tolerance: 1e-12,
                absolute_energy_tolerance_j: 1e-11,
                max_iterations: 128,
            },
        },
        absolute_mass_tolerance_kg: 1.0,
        absolute_temperature_tolerance_k: 1000.0,
        absolute_specific_enthalpy_tolerance_j_kg: 1e7,
        absolute_heat_tolerance_j: 1e4,
        relative_tolerance: 1e-3,
        minimum_dt_s: 1e-10,
        maximum_dt_s: 0.02,
        maximum_trials: 20000,
    }
}
fn support(h: f64) -> Result<f64, SublimationError> {
    Ok(400.0 + h / 500.0)
}
fn row<'a>(
    law: &'a dyn SublimationLaw,
    binding: &'a GasSurfaceBinding,
    deposit: &'a dyn Fn(f64) -> Result<f64, SublimationError>,
    chart: (f64, f64),
) -> GasResistiveSurface<'a> {
    GasResistiveSurface {
        id: 7,
        law,
        binding,
        state: ResistiveSurfaceState {
            substrate_mass_kg: 0.02,
            substrate_specific_enthalpy_j_kg: -25000.0,
            deposit: SolidSurfaceState {
                mass_kg: 2e-6,
                area_m2: 1e-4,
                specific_enthalpy_j_kg: chart.0 - 50.0 * chart.1,
            },
        },
        conductance_w_k: 0.1,
        substrate_temperature: &support,
        deposit_temperature: deposit,
    }
}

#[test]
fn live_films_follow_the_finite_gas_half_steps_and_match_manual_rebinding() {
    let catalog = fixture::catalog(0.0);
    let gas = fixture::gas(&catalog, 400.0, 3000.0);
    let original_gas = gas.clone();
    let law = fixture::phase();
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let deposit = |h| fixture::temperature(h, chart);
    let initial = row(&law, &binding, &deposit, chart);
    let vapor0 = gas.species_mass_kg(&fixture::species()).unwrap();
    let coefficient = |parcel: &ClosedGasParcel| {
        0.25 * (parcel.state().state().temperature / 400.0)
            * (parcel.species_mass_kg(&fixture::species()).unwrap() / vapor0)
    };
    let mut visits = Vec::new();
    let resolved = advance_resistive_chamber_adaptive_resolved(
        &[initial],
        None,
        |rows: &[GasResistiveSurface<'_>],
         parcel: &ClosedGasParcel,
         output: &mut BTreeMap<u64, f64>,
         _: &mut dyn FnMut() -> bool| {
            visits.push((
                rows[0].state.deposit.mass_kg,
                parcel.species_mass_kg(&fixture::species()).unwrap(),
            ));
            output.insert(7, coefficient(parcel));
            Ok::<(), &'static str>(())
        },
        &BTreeMap::new(),
        &gas,
        &catalog,
        0.02,
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(resolved.audit().chamber.accepted_intervals, 1);
    assert_eq!(resolved.audit().chamber.kernel_evaluations, 3);
    assert_eq!(visits.len(), 3);
    assert_eq!(visits[0], visits[1]);
    assert!(visits[2].0 > visits[1].0);
    assert!(visits[2].1 < visits[1].1);

    let first = advance_resistive_chamber_with_films(
        &[initial],
        None,
        &BTreeMap::from([(7, coefficient(&gas))]),
        &gas,
        &catalog,
        0.01,
        controls().step,
        || false,
    )
    .unwrap();
    let middle = GasResistiveSurface {
        state: first.surfaces()[0].state,
        ..initial
    };
    let expected = advance_resistive_chamber_with_films(
        &[middle],
        None,
        &BTreeMap::from([(7, coefficient(first.gas()))]),
        first.gas(),
        &catalog,
        0.01,
        controls().step,
        || false,
    )
    .unwrap();
    assert_eq!(resolved.endpoint().surfaces()[0].state, expected.surfaces()[0].state);
    assert_eq!(resolved.gas(), expected.gas());
    fixture::near(resolved.audit().chamber.mass_residual_kg, 0.0, 1e-14);
    fixture::near(resolved.audit().chamber.energy_residual_j, 0.0, 1e-7);
    assert_eq!(gas, original_gas);
}

#[test]
fn coefficients_must_be_refreshed_and_source_errors_retain_their_type() {
    let catalog = fixture::catalog(0.0);
    let gas = fixture::gas(&catalog, 400.0, 3000.0);
    let original_gas = gas.clone();
    let law = fixture::phase();
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let deposit = |h| fixture::temperature(h, chart);
    let initial = row(&law, &binding, &deposit, chart);
    let mut calls = 0;
    let stale = advance_resistive_chamber_adaptive_resolved(
        &[initial],
        None,
        |_: &[GasResistiveSurface<'_>],
         _: &ClosedGasParcel,
         output: &mut BTreeMap<u64, f64>,
         _: &mut dyn FnMut() -> bool| {
            calls += 1;
            if calls == 1 {
                output.insert(7, 0.0);
            }
            Ok::<(), u32>(())
        },
        &BTreeMap::new(),
        &gas,
        &catalog,
        0.02,
        controls(),
        || false,
    );
    assert!(matches!(stale, Err(GasResistiveFilmError::Advance(_))));
    assert_eq!(calls, 2);

    let source = advance_resistive_chamber_adaptive_resolved(
        &[initial],
        None,
        |_: &[GasResistiveSurface<'_>],
         _: &ClosedGasParcel,
         _: &mut BTreeMap<u64, f64>,
         _: &mut dyn FnMut() -> bool| Err::<(), u32>(73),
        &BTreeMap::new(),
        &gas,
        &catalog,
        0.02,
        controls(),
        || false,
    );
    assert!(matches!(source, Err(GasResistiveFilmError::Source(73))));
    assert_eq!(gas, original_gas);
}

#[test]
fn zero_time_admission_validates_source_identity_and_cancellation() {
    let catalog = fixture::catalog(0.0);
    let gas = fixture::gas(&catalog, 400.0, 3000.0);
    let law = fixture::phase();
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let deposit = |h| fixture::temperature(h, chart);
    let initial = row(&law, &binding, &deposit, chart);
    let unknown = advance_resistive_chamber_adaptive_resolved(
        &[initial],
        None,
        |_: &[GasResistiveSurface<'_>],
         _: &ClosedGasParcel,
         output: &mut BTreeMap<u64, f64>,
         _: &mut dyn FnMut() -> bool| {
            output.insert(99, 0.1);
            Ok::<(), u32>(())
        },
        &BTreeMap::new(),
        &gas,
        &catalog,
        0.0,
        controls(),
        || false,
    );
    assert!(matches!(unknown, Err(GasResistiveFilmError::Advance(_))));
    let cancelled = Cell::new(false);
    let result = advance_resistive_chamber_adaptive_resolved(
        &[initial],
        None,
        |_: &[GasResistiveSurface<'_>],
         _: &ClosedGasParcel,
         output: &mut BTreeMap<u64, f64>,
         _: &mut dyn FnMut() -> bool| {
            output.insert(7, 0.1);
            cancelled.set(true);
            Ok::<(), u32>(())
        },
        &BTreeMap::new(),
        &gas,
        &catalog,
        0.0,
        controls(),
        || cancelled.get(),
    );
    assert!(matches!(
        result,
        Err(GasResistiveFilmError::Advance(
            GasResistiveAdvanceError::Surface(GasSurfaceError::Surface(
                SublimationError::Cancelled
            ))
        ))
    ));
}

#[test]
fn growing_resolver_observes_the_actual_mass_owned_area() {
    let catalog = fixture::catalog(0.0);
    let gas = fixture::gas(&catalog, 400.0, 3000.0);
    let law = fixture::phase();
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let deposit = |h| fixture::temperature(h, chart);
    let mut initial = row(&law, &binding, &deposit, chart);
    initial.state.deposit.area_m2 = 7.0;
    let curve = ConductivityCurve::new(vec![[200.0, 1e-3], [900.0, 1e-3]]).unwrap();
    let layers = BTreeMap::from([(
        7,
        GasIntegratedGrowingCoating {
            geometry: ConductiveCoatingGeometry::planar(917.0, 1e-4).unwrap(),
            conductivity: &curve,
            contact_resistance_m2_k_w: 2e-4,
            interface_config: LayerConductionConfig {
                max_iterations: 96,
                absolute_heat_tolerance_w: 1e-11,
                relative_heat_tolerance: 1e-12,
            },
        },
    )]);
    let solve = GasGrowingCoatingSolve {
        advance: controls(),
        maximum_deposit_volume_fraction: 0.01,
        absolute_extent_tolerance_m: 1.0,
        absolute_conductance_tolerance_w_k: 1e4,
    };
    let mut calls = 0;
    let result = advance_growing_integrated_chamber_resolved(
        &[initial],
        &layers,
        None,
        |rows: &[GasResistiveSurface<'_>],
         _: &ClosedGasParcel,
         output: &mut BTreeMap<u64, f64>,
         _: &mut dyn FnMut() -> bool| {
            calls += 1;
            assert_eq!(rows[0].state.deposit.area_m2, 1e-4);
            assert_eq!(rows[0].state.deposit.mass_kg, initial.state.deposit.mass_kg);
            output.insert(7, 0.0);
            Ok::<(), u32>(())
        },
        None,
        &gas,
        &catalog,
        0.0,
        solve,
        || false,
    )
    .unwrap();
    assert_eq!(calls, 1);
    assert_eq!(
        result.history().endpoint().surfaces()[0].state.deposit.area_m2,
        1e-4
    );
    assert_eq!(initial.state.deposit.area_m2, 7.0);
    assert_eq!(result.gas(), &gas);
}
