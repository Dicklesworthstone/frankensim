//! Actual timeline integration using compatible manufactured NASA-9 calorics.
use super::*;
use crate::sublimation::gas::GasSurfaceSolve;
use crate::thermochemical_gas::energy::{
    GasEnergySolve, parcel::transfer::species::boundary::GasSpeciesBoundarySolve,
};
use crate::thermochemical_gas::tests::{point, species_with_energy_reference, transport};
use fs_qty::Temperature;
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, DepositionBudget, MaterialId, Reservoir, ThermalLaw, ThermalLink,
    TimedEvent,
};
use fs_thermal::sublimation::SublimationProperties;
use fs_thermochem::{
    Composition, CompositionBasis, SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K,
    mixture::FrozenIdealGasMixtureModelV1,
};

const MOLAR_MASS: f64 = 0.018015;
const LATENT: f64 = 200000.0;
struct Phase;
impl SublimationLaw for Phase {
    fn temperature_range_k(&self) -> [f64; 2] {
        [200.0, 900.0]
    }
    fn properties(&self, t: f64) -> Result<SublimationProperties, SublimationError> {
        if !(200.0..=900.0).contains(&t) {
            return Err(SublimationError::TemperatureOutOfRange {
                temperature_k: t,
                minimum_k: 200.0,
                maximum_k: 900.0,
            });
        }
        let r = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / MOLAR_MASS;
        let p = 1000.0 * fs_math::det::exp(LATENT / r * (1.0 / 400.0 - 1.0 / t));
        Ok(SublimationProperties {
            saturation_pressure_pa: p,
            saturation_pressure_derivative_pa_k: p * LATENT / (r * t * t),
            latent_heat_j_kg: LATENT,
            latent_heat_derivative_j_kg_k: 0.0,
            solid_density_kg_m3: 917.0,
            molar_mass_kg_mol: MOLAR_MASS,
        })
    }
}
#[derive(Clone, Copy)]
struct Caloric {
    h400: f64,
    cp: f64,
}
impl Caloric {
    fn h(self, t: f64) -> f64 {
        self.h400 + self.cp * (t - 400.0)
    }
}
impl ThermalLaw for Caloric {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let t = 400.0 + (h - self.h400) / self.cp;
        if (200.0..=900.0).contains(&t) {
            Ok(t)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / self.cp
    }
}
fn species() -> SpeciesId {
    SpeciesId::new("B").unwrap()
}
fn catalog(reference: f64) -> ThermochemicalGasModel {
    ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![
                species_with_energy_reference("A", 0.02897, 3.5, 0.0, 0.0),
                species_with_energy_reference("B", MOLAR_MASS, 4.0, 0.0, reference),
            ],
            Composition::new(CompositionBasis::MoleFraction, vec![0.99, 0.01]).unwrap(),
        )
        .unwrap(),
        vec![transport("A", 1.7e-5, 0.025), transport("B", 2.0e-5, 0.03)],
        Some(&species()),
    )
    .unwrap()
}
fn gas(catalog: &ThermochemicalGasModel) -> ClosedGasParcel {
    let model = catalog
        .with_composition(
            Composition::new(CompositionBasis::MoleFraction, vec![0.992, 0.008]).unwrap(),
        )
        .unwrap();
    ClosedGasParcel::new(
        model,
        point(350.0, 100000.0).with("local-lot", 9.0).unwrap(),
        0.01,
    )
    .unwrap()
}
fn binding(offset: f64) -> GasSurfaceBinding {
    GasSurfaceBinding {
        species: species(),
        accommodation_coefficient: 0.1,
        solid_enthalpy_offset_j_kg: offset,
    }
}
fn chart(catalog: &ThermochemicalGasModel, binding: &GasSurfaceBinding) -> Caloric {
    let component = catalog
        .thermodynamics()
        .components()
        .iter()
        .find(|component| component.species() == &binding.species)
        .unwrap();
    let evaluated = component.evaluate(Temperature::new(400.0)).unwrap();
    Caloric {
        cp: evaluated.properties().cp().value() / component.molar_mass().value(),
        h400: catalog
            .species_specific_enthalpy_j_kg(&binding.species, 400.0)
            .unwrap()
            - LATENT
            - binding.solid_enthalpy_offset_j_kg,
    }
}
fn controls() -> GasChamberSolve {
    GasChamberSolve {
        surface: GasSurfaceSolve {
            gas: GasSpeciesBoundarySolve {
                gas: GasEnergySolve {
                    minimum_temperature_k: 200.0,
                    maximum_temperature_k: 900.0,
                    energy_tolerance_j_kg: 1e-7,
                    maximum_iterations: 128,
                },
                mass_tolerance_kg: 1e-14,
                energy_tolerance_j: 1e-7,
            },
            maximum_gross_inventory_fraction: 0.25,
            latent_compatibility_tolerance_j_kg: 1e-5,
        },
        maximum_surfaces: 2,
    }
}
fn heat(dt: f64) -> StepConfig {
    StepConfig {
        max_dt_s: dt,
        energy_tolerance_j: 1e-7,
    }
}
fn selection<'a>(
    law: &'a Phase,
    bindings: &'a [GasSurfaceBinding; 2],
    charts: &[Caloric; 2],
) -> [GasTimelineSurface<'a>; 2] {
    std::array::from_fn(|i| GasTimelineSurface {
        cell_id: CellId(i as u64 + 1),
        law,
        binding: &bindings[i],
        pure_solid_enthalpy_range_j_kg: [charts[i].h(200.0), charts[i].h(900.0)],
    })
}
fn boundary(power: f64) -> Boundary {
    Boundary {
        ambient_k: 350.0,
        htc_w_m2_k: 0.0,
        source_w: power,
        reservoir: Some(Reservoir {
            area_m2: 1e-4,
            conductance_w_k: 0.0,
            temperature_k: 350.0,
        }),
    }
}
fn cell(id: u64, t: f64, chart: Caloric, power: f64) -> DepositedCell {
    DepositedCell {
        id: CellId(id),
        segment_id: id,
        body_id: id,
        material_id: MaterialId(id),
        mass_kg: 0.01,
        specific_h_j_kg: chart.h(t),
        geometry: CellGeometry {
            start_m: [id as f64, 0.0, 0.0],
            end_m: [id as f64, 0.01, 0.0],
            volume_m3: 0.01 / 917.0,
            surface_area_m2: 2e-4,
        },
        boundary: boundary(power),
    }
}
fn birth(time_s: f64, cell: DepositedCell) -> TimedEvent {
    TimedEvent {
        time_s,
        event: Event::Birth(cell),
    }
}
fn schedule<'a>(
    charts: &'a [Caloric; 2],
    events: Vec<TimedEvent>,
    steps: usize,
) -> DepositionSimulation<'a> {
    DepositionSimulation::new(
        events,
        &[(MaterialId(1), &charts[0]), (MaterialId(2), &charts[1])],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 8,
            max_substeps: steps,
        },
    )
    .unwrap()
}
fn initial<'a>(charts: &'a [Caloric; 2], power: f64) -> DepositionSimulation<'a> {
    schedule(
        charts,
        vec![
            birth(0.0, cell(1, 400.0, charts[0], power)),
            birth(0.0, cell(2, 300.0, charts[1], power)),
        ],
        1000,
    )
    .advance_to(0.0, heat(1.0), || Ok(()))
    .unwrap()
    .0
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:e} != {b:e}");
}

#[test]
fn original_heat_and_two_real_cell_inventories_share_one_evolving_gas() {
    let catalog = catalog(0.0);
    let bindings = [binding(0.0), binding(0.0)];
    let charts = bindings.each_ref().map(|b| chart(&catalog, b));
    let law = Phase;
    let surfaces = selection(&law, &bindings, &charts);
    let thermal = schedule(
        &charts,
        vec![
            birth(0.0, cell(1, 400.0, charts[0], 1.5)),
            birth(0.0, cell(2, 300.0, charts[1], 1.5)),
            TimedEvent {
                time_s: 0.0,
                event: Event::Contact(ThermalLink {
                    id: 7,
                    a: CellId(1),
                    b: CellId(2),
                    area_m2: 1e-5,
                    conductance_w_k: 0.05,
                }),
            },
        ],
        1000,
    )
    .advance_to(0.0, heat(1.0), || Ok(()))
    .unwrap()
    .0;
    let gas = gas(&catalog);
    let original_carrier = gas.species_mass_kg(&SpeciesId::new("A").unwrap()).unwrap();
    let candidate = advance_deposition_with_gas(
        &thermal,
        &surfaces,
        &gas,
        &catalog,
        0.1,
        heat(0.025),
        controls(),
        || Ok(()),
    )
    .unwrap();
    let audit = candidate.audit();
    assert!(audit.cells[&CellId(1)].transferred_mass_kg > 0.0);
    assert!(audit.cells[&CellId(2)].transferred_mass_kg < 0.0);
    assert!(candidate.thermal().temperature_k(CellId(1)).unwrap() < 400.0);
    assert!(candidate.thermal().temperature_k(CellId(2)).unwrap() > 300.0);
    assert_ne!(candidate.gas().state().state().temperature, 350.0);
    near(audit.thermal.external_heat_j, 0.3, 1e-12);
    near(
        candidate.thermal().mass_kg() + candidate.gas().mass_kg(),
        thermal.mass_kg() + gas.mass_kg(),
        1e-14,
    );
    near(audit.energy_residual_j, 0.0, 1e-7);
    near(audit.mass_residual_kg, 0.0, 1e-14);
    near(
        candidate
            .gas()
            .species_mass_kg(&SpeciesId::new("A").unwrap())
            .unwrap(),
        original_carrier,
        1e-14,
    );
    assert_eq!(candidate.gas().point().axes().get("local-lot"), Some(&9.0));
    assert_eq!(candidate.thermal().links(), thermal.links());
    assert_eq!(thermal.time_s(), 0.0);
}

#[test]
fn caller_order_does_not_change_the_joint_timeline() {
    let catalog = catalog(0.0);
    let bindings = [binding(0.0), binding(0.0)];
    let charts = bindings.each_ref().map(|b| chart(&catalog, b));
    let law = Phase;
    let surfaces = selection(&law, &bindings, &charts);
    let reverse = [surfaces[1], surfaces[0]];
    let thermal = initial(&charts, 0.0);
    let gas = gas(&catalog);
    let a = advance_deposition_with_gas(
        &thermal,
        &surfaces,
        &gas,
        &catalog,
        0.1,
        heat(0.025),
        controls(),
        || Ok(()),
    )
    .unwrap();
    let b = advance_deposition_with_gas(
        &thermal,
        &reverse,
        &gas,
        &catalog,
        0.1,
        heat(0.025),
        controls(),
        || Ok(()),
    )
    .unwrap();
    assert_eq!(a.thermal().cells(), b.thermal().cells());
    assert_eq!(a.gas(), b.gas());
    assert_eq!(a.audit(), b.audit());
}

#[test]
fn endpoint_birth_starts_exchange_after_birth_and_carries_its_datum_once() {
    let catalog = catalog(1.0e7);
    let bindings = [binding(1.0e7), binding(1.0e7 + 100_000.0)];
    let charts = bindings.each_ref().map(|b| chart(&catalog, b));
    let law = Phase;
    let surfaces = selection(&law, &bindings, &charts);
    let thermal = schedule(
        &charts,
        vec![
            birth(0.0, cell(1, 400.0, charts[0], 0.0)),
            birth(0.05, cell(2, 300.0, charts[1], 0.0)),
        ],
        100,
    );
    let candidate = advance_deposition_with_gas(
        &thermal,
        &surfaces,
        &gas(&catalog),
        &catalog,
        0.05,
        heat(0.05),
        controls(),
        || Ok(()),
    )
    .unwrap();
    assert!(candidate.audit().cells[&CellId(1)].transferred_mass_kg > 0.0);
    assert_eq!(candidate.audit().cells[&CellId(2)], GasCellAudit::default());
    assert_eq!(candidate.thermal().cells()[&CellId(2)].mass_kg, 0.01);
    assert_eq!(
        candidate.thermal().cells()[&CellId(2)].specific_h_j_kg,
        charts[1].h(300.0)
    );
    near(
        candidate.audit().incoming_energy_j,
        0.01 * (charts[0].h(400.0) + bindings[0].solid_enthalpy_offset_j_kg)
            + 0.01 * (charts[1].h(300.0) + bindings[1].solid_enthalpy_offset_j_kg),
        1e-7,
    );
    near(candidate.audit().energy_residual_j, 0.0, 1e-7);
}

#[test]
fn endpoint_cover_preserves_start_exposure_then_stops_further_exchange() {
    let catalog = catalog(0.0);
    let bindings = [binding(0.0), binding(0.0)];
    let charts = bindings.each_ref().map(|b| chart(&catalog, b));
    let law = Phase;
    let surfaces = selection(&law, &bindings, &charts);
    let mut covered = boundary(0.0);
    covered.reservoir.as_mut().unwrap().area_m2 = 2e-4;
    let thermal = schedule(
        &charts,
        vec![
            birth(0.0, cell(1, 400.0, charts[0], 0.0)),
            birth(0.0, cell(2, 300.0, charts[1], 0.0)),
            TimedEvent {
                time_s: 0.05,
                event: Event::Boundary(CellId(1), covered),
            },
            TimedEvent {
                time_s: 0.05,
                event: Event::Boundary(CellId(2), covered),
            },
        ],
        100,
    )
    .advance_to(0.0, heat(1.0), || Ok(()))
    .unwrap()
    .0;
    let gas = gas(&catalog);
    let expected = advance_deposition_with_gas(
        &initial(&charts, 0.0),
        &surfaces,
        &gas,
        &catalog,
        0.05,
        heat(0.1),
        controls(),
        || Ok(()),
    )
    .unwrap();
    let actual = advance_deposition_with_gas(
        &thermal,
        &surfaces,
        &gas,
        &catalog,
        0.1,
        heat(0.1),
        controls(),
        || Ok(()),
    )
    .unwrap();
    for id in [CellId(1), CellId(2)] {
        assert_eq!(
            actual.thermal().cells()[&id].mass_kg,
            expected.thermal().cells()[&id].mass_kg
        );
        assert_eq!(
            actual.thermal().cells()[&id].specific_h_j_kg,
            expected.thermal().cells()[&id].specific_h_j_kg
        );
        assert_eq!(actual.thermal().exposed_area_m2(id), Ok(0.0));
    }
    assert_eq!(actual.gas(), expected.gas());
}

#[test]
fn late_cancellation_discards_the_complete_heat_phase_gas_candidate() {
    let catalog = catalog(0.0);
    let bindings = [binding(0.0), binding(0.0)];
    let charts = bindings.each_ref().map(|b| chart(&catalog, b));
    let law = Phase;
    let surfaces = selection(&law, &bindings, &charts);
    let thermal = initial(&charts, 1.5);
    let before = thermal.cells().clone();
    let gas = gas(&catalog);
    let saved_gas = gas.clone();
    let mut calls = 0;
    advance_deposition_with_gas(
        &thermal,
        &surfaces,
        &gas,
        &catalog,
        0.1,
        heat(0.025),
        controls(),
        || {
            calls += 1;
            Ok(())
        },
    )
    .unwrap();
    let mut count = 0;
    let result = advance_deposition_with_gas(
        &thermal,
        &surfaces,
        &gas,
        &catalog,
        0.1,
        heat(0.025),
        controls(),
        || {
            count += 1;
            if count == calls {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        },
    );
    assert!(matches!(
        result,
        Err(GasTimelineError::Thermal(DepositionError::Cancelled))
    ));
    assert_eq!(thermal.cells(), &before);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(gas, saved_gas);
}

#[test]
fn splitting_cannot_renew_the_original_thermal_budget() {
    let catalog = catalog(0.0);
    let bindings = [binding(0.0), binding(0.0)];
    let charts = bindings.each_ref().map(|b| chart(&catalog, b));
    let law = Phase;
    let surfaces = selection(&law, &bindings, &charts);
    let thermal = schedule(
        &charts,
        vec![
            birth(0.0, cell(1, 400.0, charts[0], 0.0)),
            birth(0.0, cell(2, 300.0, charts[1], 0.0)),
        ],
        1,
    );
    let result = advance_deposition_with_gas(
        &thermal,
        &surfaces,
        &gas(&catalog),
        &catalog,
        0.1,
        heat(0.025),
        controls(),
        || Ok(()),
    );
    assert!(matches!(
        result,
        Err(GasTimelineError::Thermal(DepositionError::Budget(_)))
    ));
    assert!(thermal.cells().is_empty());
    assert_eq!(thermal.pending_events().len(), 2);
}

#[test]
fn heat_cannot_leave_the_declared_solid_branch_before_phase_cools_it() {
    let catalog = catalog(0.0);
    let bindings = [binding(0.0), binding(0.0)];
    let charts = bindings.each_ref().map(|b| chart(&catalog, b));
    let law = Phase;
    let mut surfaces = selection(&law, &bindings, &charts);
    surfaces[0].pure_solid_enthalpy_range_j_kg[1] = charts[0].h(400.0) + 1.0;
    let thermal = initial(&charts, 100.0);
    let before = thermal.cells().clone();
    let result = advance_deposition_with_gas(
        &thermal,
        &surfaces,
        &gas(&catalog),
        &catalog,
        0.1,
        heat(0.1),
        controls(),
        || Ok(()),
    );
    assert!(matches!(
        result,
        Err(GasTimelineError::Thermal(DepositionError::MaterialDomain))
    ));
    assert_eq!(thermal.cells(), &before);
}

#[test]
fn fixed_gas_formation_reference_does_not_change_the_coupled_physics() {
    let base = catalog(0.0);
    let shifted = catalog(1.0e7);
    let base_bindings = [binding(0.0), binding(0.0)];
    let shifted_bindings = [binding(1.0e7), binding(1.0e7)];
    let base_charts = base_bindings.each_ref().map(|b| chart(&base, b));
    let shifted_charts = shifted_bindings.each_ref().map(|b| chart(&shifted, b));
    let law = Phase;
    let a = advance_deposition_with_gas(
        &initial(&base_charts, 1.5),
        &selection(&law, &base_bindings, &base_charts),
        &gas(&base),
        &base,
        0.1,
        heat(0.025),
        controls(),
        || Ok(()),
    )
    .unwrap();
    let b = advance_deposition_with_gas(
        &initial(&shifted_charts, 1.5),
        &selection(&law, &shifted_bindings, &shifted_charts),
        &gas(&shifted),
        &shifted,
        0.1,
        heat(0.025),
        controls(),
        || Ok(()),
    )
    .unwrap();
    for id in [CellId(1), CellId(2)] {
        near(
            a.thermal().temperature_k(id).unwrap(),
            b.thermal().temperature_k(id).unwrap(),
            1e-7,
        );
        near(
            a.thermal().cells()[&id].mass_kg,
            b.thermal().cells()[&id].mass_kg,
            1e-14,
        );
    }
    near(
        a.gas().state().state().temperature,
        b.gas().state().state().temperature,
        1e-7,
    );
    near(a.audit().energy_residual_j, 0.0, 1e-7);
    near(b.audit().energy_residual_j, 0.0, 1e-7);
}

#[test]
fn unknown_duplicate_and_empty_selections_leave_both_owners_unchanged() {
    let catalog = catalog(0.0);
    let bindings = [binding(0.0), binding(0.0)];
    let charts = bindings.each_ref().map(|b| chart(&catalog, b));
    let law = Phase;
    let surfaces = selection(&law, &bindings, &charts);
    let thermal = initial(&charts, 0.0);
    let before = thermal.cells().clone();
    let gas = gas(&catalog);
    let saved_gas = gas.clone();
    let mut unknown = surfaces[0];
    unknown.cell_id = CellId(99);
    for invalid in [vec![surfaces[0], surfaces[0]], vec![unknown], vec![]] {
        assert!(
            advance_deposition_with_gas(
                &thermal,
                &invalid,
                &gas,
                &catalog,
                0.1,
                heat(0.025),
                controls(),
                || Ok(())
            )
            .is_err()
        );
    }
    assert_eq!(thermal.cells(), &before);
    assert_eq!(gas, saved_gas);
}

#[test]
fn empty_prebirth_time_jumps_without_spending_unbounded_phase_intervals() {
    let catalog = catalog(0.0);
    let bindings = [binding(0.0), binding(0.0)];
    let charts = bindings.each_ref().map(|b| chart(&catalog, b));
    let law = Phase;
    let surfaces = selection(&law, &bindings, &charts);
    let thermal = schedule(
        &charts,
        vec![
            birth(1000.0, cell(1, 400.0, charts[0], 0.0)),
            birth(1000.0, cell(2, 300.0, charts[1], 0.0)),
        ],
        0,
    );
    let gas = gas(&catalog);
    let mut calls = 0;
    let candidate = advance_deposition_with_gas(
        &thermal,
        &surfaces,
        &gas,
        &catalog,
        1000.0,
        heat(1e-12),
        controls(),
        || {
            calls += 1;
            if calls > 200 {
                Err(DepositionError::Budget("unbounded prebirth phase loop"))
            } else {
                Ok(())
            }
        },
    )
    .unwrap();
    assert_eq!(candidate.thermal().time_s(), 1000.0);
    assert_eq!(candidate.audit().thermal.substeps, 0);
    assert_eq!(candidate.audit().phase_steps, 0);
    assert_eq!(candidate.gas(), &gas);
    assert_eq!(candidate.thermal().cells()[&CellId(1)].mass_kg, 0.01);
    assert_eq!(candidate.thermal().cells()[&CellId(2)].mass_kg, 0.01);
}
