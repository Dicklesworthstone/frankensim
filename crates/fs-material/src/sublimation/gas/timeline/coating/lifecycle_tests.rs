//! G0/G1 lifecycle wiring on the original heat timeline and finite gas owners.
//! Manufactured compatible calorics establish accounting, not empirical fidelity.

use super::adaptive::advance_coating_lifecycle_adaptive_with_gas;
use super::*;
use crate::sublimation::gas::tests as fixture;
use crate::sublimation::gas::timeline::adaptive::GasTimelineAdaptiveSolve;
use fs_thermal::enthalpy::sublimation::coating::CoatingSurface;
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, DepositedCell, DepositionBudget, DepositionSimulation, MaterialId,
    ThermalLaw, ThermalLink, TimedEvent,
};
use fs_thermal::sublimation::{SolverConfig, SublimationLaw, SurfaceKinetics};
use fs_thermochem::SpeciesId;

struct Substrate;
impl ThermalLaw for Substrate {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let t = 400.0 + h / 800.0;
        if (200.0..=900.0).contains(&t) {
            Ok(t)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / 800.0
    }
}

fn boundary(power_w: f64) -> Boundary {
    Boundary {
        ambient_k: 400.0,
        htc_w_m2_k: 0.0,
        source_w: power_w,
        reservoir: None,
    }
}
fn cell(id: u64, power_w: f64) -> DepositedCell {
    DepositedCell {
        id: CellId(id),
        segment_id: id,
        body_id: id,
        material_id: MaterialId(1),
        mass_kg: 0.01,
        specific_h_j_kg: 0.0,
        geometry: CellGeometry {
            start_m: [id as f64, 0.0, 0.0],
            end_m: [id as f64, 0.01, 0.0],
            volume_m3: 0.01 / 2700.0,
            surface_area_m2: 1e-4,
        },
        boundary: boundary(power_w),
    }
}
fn birth(time_s: f64, id: u64, power_w: f64) -> TimedEvent {
    TimedEvent {
        time_s,
        event: Event::Birth(cell(id, power_w)),
    }
}
fn initial<'a>(
    substrate: &'a Substrate,
    law: &'a dyn SublimationLaw,
    deposit_h: &'a dyn Fn(f64) -> Result<f64, SublimationError>,
    binding: &GasSurfaceBinding,
    events: Vec<TimedEvent>,
    initial_mass_kg: f64,
) -> CoatingSimulation<'a> {
    let thermal = DepositionSimulation::new(
        events,
        &[(MaterialId(1), substrate)],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 8,
            max_substeps: 256,
        },
    )
    .unwrap();
    CoatingSimulation::new(
        thermal,
        &[CoatingSurface {
            cell_id: CellId(1),
            law,
            deposit_enthalpy: deposit_h,
            kinetics: SurfaceKinetics {
                accommodation_coefficient: binding.accommodation_coefficient,
                mass_transfer_coefficient_m_s: None,
            },
            initial_mass_kg,
            reservoir: VaporReservoir::Prescribed {
                partial_pressure_pa: 1000.0,
                temperature_k: 400.0,
            },
        }],
        || Ok(()),
    )
    .unwrap()
}
fn selected(binding: &GasSurfaceBinding) -> GasCoatingTimelineSurface<'_> {
    GasCoatingTimelineSurface {
        cell_id: CellId(1),
        binding,
    }
}
fn controls() -> GasSupportedLifecycleSolve {
    let original = fixture::controls();
    GasSupportedLifecycleSolve {
        gas: original.gas,
        storage: SolverConfig {
            relative_tolerance: 1e-12,
            absolute_energy_tolerance_j: 1e-9,
            max_iterations: 128,
        },
        latent_compatibility_tolerance_j_kg: original.latent_compatibility_tolerance_j_kg,
        maximum_deposit_volume_fraction: 0.01,
    }
}
fn heat(dt: f64) -> StepConfig {
    StepConfig {
        max_dt_s: dt,
        energy_tolerance_j: 1e-7,
    }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:e} != {b:e}");
}
fn energy(coating: &CoatingSimulation<'_>, gas: &ClosedGasParcel) -> f64 {
    // All fixtures use zero binding offset. Read actual accepted storage rather
    // than reproducing the phase transfer's reported molecular-energy formula.
    coating.energy_j() + gas.mass_kg() * gas.state().specific_internal_energy_j_kg()
}
fn check_balance(
    original: &CoatingSimulation<'_>,
    gas: &ClosedGasParcel,
    accepted: &GasCoatingTimelineCandidate<'_, '_>,
) {
    let report = accepted.audit();
    close(report.mass_residual_kg, 0.0, 1e-14);
    close(report.energy_residual_j, 0.0, 1e-7);
    close(
        energy(accepted.coating(), accepted.gas())
            - energy(original, gas)
            - report.incoming_energy_j
            - report.thermal.external_heat_j,
        0.0,
        1e-7,
    );
    let carrier = SpeciesId::new("A").unwrap();
    close(
        accepted.gas().species_mass_kg(&carrier).unwrap(),
        gas.species_mass_kg(&carrier).unwrap(),
        1e-14,
    );
    assert_eq!(accepted.coating().cells()[&CellId(1)].mass_kg, 0.01);
}

#[test]
fn g0_undersaturated_bare_support_has_turnover_without_spurious_film() {
    let catalog = fixture::catalog(0.0);
    let gas = fixture::gas(&catalog, 400.0, 10.0);
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let deposit_h = |t| Ok(chart.0 + chart.1 * (t - 400.0));
    let substrate = Substrate;
    let law = fixture::phase();
    let original = initial(
        &substrate,
        &law,
        &deposit_h,
        &binding,
        vec![birth(0.0, 1, 0.0)],
        0.0,
    );
    let accepted = advance_coating_lifecycle_with_gas(
        &original,
        selected(&binding),
        &gas,
        &catalog,
        0.01,
        heat(0.01),
        controls(),
        || Ok(()),
    )
    .unwrap();
    assert_eq!(accepted.coating().deposit_mass_kg(), 0.0);
    assert!(accepted.audit().incoming_mass_kg > 0.0);
    assert_eq!(
        accepted.audit().outgoing_mass_kg,
        accepted.audit().incoming_mass_kg
    );
    assert_eq!(accepted.audit().transferred_mass_kg, 0.0);
    assert_eq!(accepted.coating().temperature_k(CellId(1)).unwrap(), 400.0);
    assert_eq!(accepted.gas(), &gas);
    check_balance(&original, &gas, &accepted);
    assert!(original.cells().is_empty() && original.coatings().is_empty());
}

#[test]
fn g1_birth_then_heat_driven_extinction_retains_contact_and_external_energy() {
    let catalog = fixture::catalog(0.0);
    let gas = fixture::gas(&catalog, 400.0, 1100.0);
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let deposit_h = |t| Ok(chart.0 + chart.1 * (t - 400.0));
    let substrate = Substrate;
    let law = fixture::phase();
    let original = initial(
        &substrate,
        &law,
        &deposit_h,
        &binding,
        vec![
            birth(0.0, 1, 0.0),
            birth(0.0, 2, 0.0),
            TimedEvent {
                time_s: 0.0,
                event: Event::Contact(ThermalLink {
                    id: 7,
                    a: CellId(1),
                    b: CellId(2),
                    area_m2: 1e-5,
                    conductance_w_k: 2.0,
                }),
            },
            TimedEvent {
                time_s: 0.001,
                event: Event::Boundary(CellId(1), boundary(1.6e6)),
            },
            TimedEvent {
                time_s: 0.002,
                event: Event::Boundary(CellId(1), boundary(0.0)),
            },
        ],
        0.0,
    );
    let grown = advance_coating_lifecycle_with_gas(
        &original,
        selected(&binding),
        &gas,
        &catalog,
        0.001,
        heat(0.001),
        controls(),
        || Ok(()),
    )
    .unwrap();
    assert!(grown.coating().deposit_mass_kg() > 0.0);
    assert!(grown.audit().outgoing_mass_kg > 0.0);
    assert_eq!(grown.audit().thermal.external_heat_j, 0.0);
    check_balance(&original, &gas, &grown);

    let bare = advance_coating_lifecycle_with_gas(
        grown.coating(),
        selected(&binding),
        grown.gas(),
        &catalog,
        0.004,
        heat(0.001),
        controls(),
        || Ok(()),
    )
    .unwrap();
    assert_eq!(bare.coating().deposit_mass_kg(), 0.0);
    close(
        bare.audit().transferred_mass_kg,
        grown.coating().deposit_mass_kg(),
        1e-14,
    );
    close(bare.audit().thermal.external_heat_j, 1600.0, 1e-9);
    assert!(bare.coating().temperature_k(CellId(2)).unwrap() > 400.0);
    assert_eq!(bare.coating().cells()[&CellId(2)].mass_kg, 0.01);
    assert_eq!(bare.coating().links().len(), 1);
    check_balance(grown.coating(), grown.gas(), &bare);

    let continued = advance_coating_lifecycle_with_gas(
        bare.coating(),
        selected(&binding),
        bare.gas(),
        &catalog,
        0.006,
        heat(0.001),
        controls(),
        || Ok(()),
    )
    .unwrap();
    assert_eq!(continued.coating().deposit_mass_kg(), 0.0);
    assert_eq!(continued.audit().transferred_mass_kg, 0.0);
    assert!(continued.audit().incoming_mass_kg > 0.0);
    assert!(continued.audit().energy_into_gas_j > 0.0);
    assert!(
        continued.coating().temperature_k(CellId(2)).unwrap()
            > bare.coating().temperature_k(CellId(2)).unwrap()
    );
    check_balance(bare.coating(), bare.gas(), &continued);
}

#[test]
fn g0_adaptive_lifecycle_accepts_the_original_two_half_interval_states() {
    let catalog = fixture::catalog(0.0);
    let gas = fixture::gas(&catalog, 400.0, 2000.0);
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let deposit_h = |t| Ok(chart.0 + chart.1 * (t - 400.0));
    let substrate = Substrate;
    let law = fixture::phase();
    let original = initial(
        &substrate,
        &law,
        &deposit_h,
        &binding,
        vec![birth(0.0, 1, 2.0)],
        0.0,
    );
    let first = advance_coating_lifecycle_with_gas(
        &original,
        selected(&binding),
        &gas,
        &catalog,
        0.001,
        heat(0.001),
        controls(),
        || Ok(()),
    )
    .unwrap();
    let second = advance_coating_lifecycle_with_gas(
        first.coating(),
        selected(&binding),
        first.gas(),
        &catalog,
        0.002,
        heat(0.001),
        controls(),
        || Ok(()),
    )
    .unwrap();
    let adaptive = advance_coating_lifecycle_adaptive_with_gas(
        &original,
        selected(&binding),
        &gas,
        &catalog,
        0.002,
        heat(0.002),
        controls(),
        GasTimelineAdaptiveSolve {
            absolute_mass_tolerance_kg: 1e-5,
            absolute_temperature_tolerance_k: 1.0,
            absolute_specific_enthalpy_tolerance_j_kg: 1000.0,
            relative_tolerance: 0.0,
            minimum_dt_s: 1e-9,
            maximum_trials: 8,
        },
        || Ok(()),
    )
    .unwrap();
    assert_eq!(adaptive.audit().accepted_intervals, 1);
    assert_eq!(adaptive.audit().rejected_trials, 0);
    assert_eq!(
        adaptive.endpoint().coating().cells(),
        second.coating().cells()
    );
    assert_eq!(
        adaptive.endpoint().coating().coatings(),
        second.coating().coatings()
    );
    assert_eq!(adaptive.endpoint().gas(), second.gas());
    assert!(adaptive.endpoint().coating().deposit_mass_kg() > 0.0);
    close(
        adaptive.endpoint().audit().thermal.external_heat_j,
        0.004,
        1e-12,
    );
    close(
        adaptive.endpoint().audit().incoming_mass_kg,
        first.audit().incoming_mass_kg + second.audit().incoming_mass_kg,
        1e-20,
    );
    check_balance(&original, &gas, adaptive.endpoint());
}

#[test]
fn g0_delayed_substrate_birth_receives_only_its_postbirth_exposure() {
    let catalog = fixture::catalog(0.0);
    let gas = fixture::gas(&catalog, 400.0, 2000.0);
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let deposit_h = |t| Ok(chart.0 + chart.1 * (t - 400.0));
    let substrate = Substrate;
    let law = fixture::phase();
    let delayed = initial(
        &substrate,
        &law,
        &deposit_h,
        &binding,
        vec![birth(0.01, 1, 2.0)],
        0.0,
    );
    let born = advance_coating_lifecycle_with_gas(
        &delayed,
        selected(&binding),
        &gas,
        &catalog,
        0.01,
        heat(0.01),
        controls(),
        || Ok(()),
    )
    .unwrap();
    assert_eq!(born.coating().deposit_mass_kg(), 0.0);
    assert_eq!(born.audit().incoming_mass_kg, 0.0);
    assert_eq!(born.audit().outgoing_mass_kg, 0.0);
    assert_eq!(born.audit().thermal.external_heat_j, 0.0);
    assert_eq!(born.gas(), &gas);
    let after_birth = advance_coating_lifecycle_with_gas(
        born.coating(),
        selected(&binding),
        born.gas(),
        &catalog,
        0.02,
        heat(0.01),
        controls(),
        || Ok(()),
    )
    .unwrap();
    let immediate = initial(
        &substrate,
        &law,
        &deposit_h,
        &binding,
        vec![birth(0.0, 1, 2.0)],
        0.0,
    );
    let exposed = advance_coating_lifecycle_with_gas(
        &immediate,
        selected(&binding),
        &gas,
        &catalog,
        0.01,
        heat(0.01),
        controls(),
        || Ok(()),
    )
    .unwrap();
    assert_eq!(after_birth.coating().cells(), exposed.coating().cells());
    assert_eq!(
        after_birth.coating().coatings(),
        exposed.coating().coatings()
    );
    assert_eq!(after_birth.gas(), exposed.gas());
    assert!(after_birth.coating().deposit_mass_kg() > 0.0);
    close(after_birth.audit().thermal.external_heat_j, 0.02, 1e-12);
    check_balance(born.coating(), born.gas(), &after_birth);
}

#[test]
fn g0_zero_time_admits_original_sources_and_late_cancellation_is_transactional() {
    let catalog = fixture::catalog(0.0);
    let gas = fixture::gas(&catalog, 400.0, 2000.0);
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let deposit_h = |t| Ok(chart.0 + chart.1 * (t - 400.0));
    let substrate = Substrate;
    let law = fixture::phase();
    let original = initial(
        &substrate,
        &law,
        &deposit_h,
        &binding,
        vec![birth(0.0, 1, 2.0)],
        0.0,
    );
    let zero = advance_coating_lifecycle_with_gas(
        &original,
        selected(&binding),
        &gas,
        &catalog,
        0.0,
        heat(0.001),
        controls(),
        || Ok(()),
    )
    .unwrap();
    assert_eq!(zero.gas(), &gas);
    assert_eq!(zero.coating().deposit_mass_kg(), 0.0);
    assert_eq!(zero.audit().phase_steps, 0);
    let changed_source = fixture::catalog(1000.0);
    assert!(
        advance_coating_lifecycle_with_gas(
            zero.coating(),
            selected(&binding),
            &gas,
            &changed_source,
            0.0,
            heat(0.001),
            controls(),
            || Ok(()),
        )
        .is_err()
    );

    let mut polls = 0;
    let accepted = advance_coating_lifecycle_with_gas(
        &original,
        selected(&binding),
        &gas,
        &catalog,
        0.002,
        heat(0.001),
        controls(),
        || {
            polls += 1;
            Ok(())
        },
    )
    .unwrap();
    check_balance(&original, &gas, &accepted);
    assert!(polls > 2);
    let mut repeated = 0;
    let refused = advance_coating_lifecycle_with_gas(
        &original,
        selected(&binding),
        &gas,
        &catalog,
        0.002,
        heat(0.001),
        controls(),
        || {
            repeated += 1;
            if repeated == polls - 1 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        },
    );
    assert!(matches!(
        refused,
        Err(GasCoatingTimelineError::Thermal(
            ThermalSublimationError::Thermal(DepositionError::Cancelled)
        ))
    ));
    assert!(original.cells().is_empty() && original.coatings().is_empty());
    assert_eq!(original.pending_events().len(), 1);
    assert_eq!(original.time_s(), 0.0);
    assert_eq!(gas, fixture::gas(&catalog, 400.0, 2000.0));
}
