//! G0/G1/G3: the original phase density reaches actual heat-transfer geometry.
use fs_blake3::ContentHash;
use fs_material::phase::vaporization::IsobaricWaterVaporization;
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
use fs_qty::Temperature;
use fs_thermal::enthalpy::conduction::{
    ConductionObserver, ConductiveLink, ConductivityLaw, PlanarDensityLayer,
};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, StepConfig, ThermalLaw, ThermalLink, TimedEvent,
};
use fs_thermochem::water::SaturatedWaterState;

struct Conductivity;
impl ConductivityLaw for Conductivity {
    fn conductivity_w_m_k(&self, _: f64) -> Result<f64, DepositionError> {
        Ok(5.0)
    }
}

fn curve() -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([0x38; 32]),
        [
            (0.0, 300.0, 0.0, 1000.0),
            (100_000.0, 400.0, 0.0, 900.0),
            (200_000.0, 400.0, 1.0, 600.0),
            (300_000.0, 500.0, 1.0, 500.0),
        ]
        .into_iter()
        .map(|(h, t, fraction, density)| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: h,
            temperature_k: t,
            liquid_mass_fraction: fraction,
            bulk_density_kg_m3: density,
        })
        .collect(),
    )
    .unwrap()
}

fn config() -> StepConfig {
    StepConfig {
        max_dt_s: 0.1,
        energy_tolerance_j: 1e-7,
    }
}

#[test]
fn latent_density_refreshes_both_link_sides_without_changing_cell_inventory() {
    let source = curve();
    let area = 0.01;
    let cells: Vec<_> = [100_000.0, 300_000.0]
        .into_iter()
        .enumerate()
        .map(|(i, h)| {
            let density = ThermalLaw::bulk_density_kg_m3(&source, h).unwrap().unwrap();
            DepositedCell {
                id: CellId(i as u64),
                segment_id: i as u64,
                body_id: 0,
                material_id: MaterialId(0),
                mass_kg: 1.0,
                specific_h_j_kg: h,
                geometry: CellGeometry {
                    start_m: [i as f64, 0.0, 0.0],
                    end_m: [i as f64, 0.1, 0.0],
                    volume_m3: 1.0 / density,
                    surface_area_m2: 2.0 * area,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: if i == 0 { 50_000.0 } else { 0.0 },
                },
            }
        })
        .collect();
    let mut events: Vec<_> = cells
        .iter()
        .cloned()
        .map(|cell| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(cell),
        })
        .collect();
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: 1,
            a: CellId(0),
            b: CellId(1),
            area_m2: area,
            conductance_w_k: 0.0,
        }),
    });
    let initial = DepositionSimulation::new(
        events,
        &[(MaterialId(0), &source)],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
            max_substeps: 100,
        },
    )
    .unwrap()
    .advance_to(0.0, config(), || Ok(()))
    .unwrap()
    .0;
    let observer = ConductionObserver::new(
        &initial,
        &[(MaterialId(0), &Conductivity)],
        &[ConductiveLink {
            link_id: 1,
            half_distances_m: [0.5 / (900.0 * area), 0.5 / (500.0 * area)],
            contact_resistance_m2_k_w: 0.0,
        }],
        &[],
        || Ok(()),
    )
    .unwrap()
    .with_material_density_layers(
        &initial,
        &[
            PlanarDensityLayer {
                cell_id: CellId(0),
                footprint_area_m2: area,
            },
            PlanarDensityLayer {
                cell_id: CellId(1),
                footprint_area_m2: area,
            },
        ],
        || Ok(()),
    )
    .unwrap();

    // Isolate a source-driven latent change on the same thermal owner. Its
    // default link is insulating; the prepared observer is still reusable.
    let midway = initial.advance_to(1.0, config(), || Ok(())).unwrap().0;
    assert_eq!(midway.temperature_k(CellId(0)).unwrap(), 400.0);
    let h = midway.cells()[&CellId(0)].specific_h_j_kg;
    let rho = source
        .values_at_specific_enthalpy(h)
        .unwrap()
        .bulk_density_kg_m3;
    let initial_power = -initial
        .power_balance_observed(&observer, || Ok(()))
        .unwrap()
        .net_power_w[&CellId(1)];
    let changed_power = -midway
        .power_balance_observed(&observer, || Ok(()))
        .unwrap()
        .net_power_w[&CellId(1)];
    let expected = 100.0 * area / (0.5 / (rho * area * 5.0) + 0.5 / (500.0 * area * 5.0));
    assert!((changed_power - expected).abs() < 1e-12);
    assert!(changed_power < 0.93 * initial_power);

    // Evolve the actual conduction observer through the same phase interval.
    let (heated, _, audit) = initial
        .advance_observed(&observer, 0.5, config(), || Ok(()))
        .unwrap();
    assert_eq!(heated.temperature_k(CellId(0)).unwrap(), 400.0);
    assert!(heated.cells()[&CellId(1)].specific_h_j_kg < 300_000.0);
    assert!((heated.energy_j() - initial.energy_j() - 25_000.0).abs() < 1e-7);
    assert!(audit.residual_j.abs() < 1e-7);
    for (i, original) in cells.iter().enumerate() {
        let current = &heated.cells()[&CellId(i as u64)];
        assert_eq!(current.mass_kg, original.mass_kg);
        assert_eq!(current.geometry, original.geometry);
    }
}

#[test]
fn original_water_vapor_inventory_is_not_averaged_as_density_or_temperature() {
    let saturation = SaturatedWaterState::at_temperature(Temperature::new(373.15)).unwrap();
    let law = IsobaricWaterVaporization::from_saturation(saturation).unwrap();
    let liquid_h = saturation.liquid_specific_enthalpy().value();
    let vapor_h = saturation.vapor_specific_enthalpy().value();
    let h = liquid_h + 0.01 * (vapor_h - liquid_h);
    let actual = ThermalLaw::bulk_density_kg_m3(&law, h).unwrap().unwrap();
    let expected = (0.99 / saturation.liquid_density().value()
        + 0.01 / saturation.vapor_density().value())
    .recip();
    assert!((actual / expected - 1.0).abs() < 2e-14);
    assert_eq!(law.temperature_k(liquid_h), law.temperature_k(h));
    assert!(actual < 0.1 * saturation.liquid_density().value());
    assert!(ThermalLaw::bulk_density_kg_m3(&law, vapor_h + 1.0).is_err());
}
