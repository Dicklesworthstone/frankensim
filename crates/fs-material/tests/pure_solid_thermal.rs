//! Whole-solid histories use the original phase chart and thermal/phase solvers.
use fs_blake3::ContentHash;
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase};
use fs_material::thermal::PureSolidThermalLaw;
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionError,
    DepositionSimulation, Event, MaterialId, StepConfig, ThermalLaw, TimedEvent,
    sublimation::{
        ThermalSublimationError,
        chamber::ChamberSurface,
        external::SolidCellExchange,
    },
};
use fs_thermal::sublimation::{
    SolidSurfaceState, SolverConfig, SublimationError, SublimationLaw, SublimationProperties,
    SurfaceKinetics, VaporReservoir, chamber::ChamberSolverConfig,
};

fn curve(shift: f64) -> EquilibriumEnthalpyPhaseCurve {
    let data = [(-140_000.0, 180.0, 0.0), (0.0, 250.0, 0.0),
        (46_000.0, 273.0, 0.0), (100_000.0, 273.0, 1.0), (100_002.0, 275.0, 1.0)];
    EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([0x91; 32]), data.into_iter()
        .map(|(h, t, f)| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: h + shift, temperature_k: t,
            liquid_mass_fraction: f, bulk_density_kg_m3: 917.0,
        }).collect()).unwrap()
}
fn cell(id: u64, h: f64, power: f64) -> DepositedCell {
    DepositedCell {
        id: CellId(id), segment_id: id, body_id: 1, material_id: MaterialId(1),
        mass_kg: 0.02, specific_h_j_kg: h,
        geometry: CellGeometry {
            start_m: [id as f64, 0.0, 0.0], end_m: [id as f64, 0.01, 0.0],
            volume_m3: 0.02 / 917.0, surface_area_m2: 0.001,
        },
        boundary: Boundary {
            ambient_k: 250.0, htc_w_m2_k: 0.0, reservoir: None, source_w: power,
        },
    }
}
fn config() -> StepConfig { StepConfig { max_dt_s: 0.1, energy_tolerance_j: 1e-7 } }
fn initial<'a>(law: &'a dyn ThermalLaw, cells: Vec<DepositedCell>) -> DepositionSimulation<'a> {
    DepositionSimulation::new(cells.into_iter().map(|cell| TimedEvent {
        time_s: 0.0, event: Event::Birth(cell),
    }).collect(), &[(MaterialId(1), law)], DepositionBudget {
        max_cells: 4, max_links: 4, max_events: 8, max_substeps: 100,
    }).unwrap().advance_to(0.0, config(), || Ok(())).unwrap().0
}
struct Ice;
impl SublimationLaw for Ice {
    fn temperature_range_k(&self) -> [f64; 2] { [180.0, 275.0] }
    fn properties(&self, t: f64) -> Result<SublimationProperties, SublimationError> {
        if !(180.0..=275.0).contains(&t) {
            return Err(SublimationError::Invalid("synthetic vapor domain"));
        }
        let p = 100.0 * fs_math::det::exp(6000.0 * (1.0 / 250.0 - 1.0 / t));
        Ok(SublimationProperties {
            saturation_pressure_pa: p, saturation_pressure_derivative_pa_k: p * 6000.0 / (t * t),
            latent_heat_j_kg: 8.314_462_618_153_24 / 0.018 * 6000.0,
            latent_heat_derivative_j_kg_k: 0.0, solid_density_kg_m3: 917.0,
            molar_mass_kg_mol: 0.018,
        })
    }
}
fn kinetics() -> SurfaceKinetics {
    SurfaceKinetics { accommodation_coefficient: 0.01, mass_transfer_coefficient_m_s: None }
}
fn phase_config() -> SolverConfig {
    SolverConfig { relative_tolerance: 1e-11, absolute_energy_tolerance_j: 1e-9, max_iterations: 96 }
}
fn chamber_config() -> ChamberSolverConfig {
    ChamberSolverConfig { surface_solver: phase_config(), relative_mass_tolerance: 1e-10,
        absolute_mass_tolerance_kg: 1e-16, max_iterations: 96, max_surfaces: 4 }
}
fn chamber() -> VaporReservoir {
    VaporReservoir::Finite { vapor_mass_kg: 0.0, volume_m3: 0.01, temperature_k: 250.0 }
}
fn vacuum() -> VaporReservoir {
    VaporReservoir::Prescribed { partial_pressure_pa: 0.0, temperature_k: 250.0 }
}

#[test]
fn solid_domain_retains_original_calorimetry_not_the_melting_temperature_test() {
    for shift in [0.0, 1_048_576.0] {
        let source = curve(shift);
        let solid = source.pure_solid_thermal_law().unwrap();
        assert!(std::ptr::eq(solid.curve(), &source));
        assert_eq!(solid.enthalpy_range_j_kg(), [-140_000.0 + shift, 46_000.0 + shift]);
        assert_eq!(solid.maximum_temperature_slope(), 0.0005);
        assert_eq!(ThermalLaw::maximum_temperature_slope(&source), 1.0);
        for h in [-140_000.0, -20_000.0, 0.0, 45_000.0, 46_000.0] {
            assert_eq!(solid.temperature_k(h + shift), ThermalLaw::temperature_k(&source, h + shift));
        }
        let onset = 46_000.0 + shift;
        for h in [onset.next_up(), 70_000.0 + shift, 100_000.0 + shift] {
            assert_eq!(ThermalLaw::temperature_k(&source, h), Ok(273.0));
            assert_eq!(solid.temperature_k(h), Err(DepositionError::MaterialDomain));
            assert_eq!(solid.specific_entropy_change_j_kg_k(shift, h), Err(DepositionError::MaterialDomain));
        }
        for (a, b) in [(-20_000.0, 20_000.0), (20_000.0, -20_000.0), (46_000.0, 46_000.0)] {
            assert_eq!(solid.specific_entropy_change_j_kg_k(a + shift, b + shift),
                ThermalLaw::specific_entropy_change_j_kg_k(&source, a + shift, b + shift));
        }
        for h in [f64::NAN, f64::INFINITY, (-140_000.0 + shift).next_down()] {
            assert_eq!(solid.temperature_k(h), Err(DepositionError::MaterialDomain));
        }
    }
}

#[test]
fn original_single_and_shared_surface_solves_keep_the_same_solid_response() {
    let source = curve(0.0);
    let solid = PureSolidThermalLaw::new(&source).unwrap();
    let mut bare = initial(&source, vec![cell(1, 0.0, 0.0), cell(2, -20_000.0, 0.0)]);
    let mut guarded = initial(&solid, vec![cell(1, 0.0, 0.0), cell(2, -20_000.0, 0.0)]);
    let (mut a, mut b) = (vacuum(), vacuum());
    let expected = bare.sublimate_cell(CellId(1), &Ice, kinetics(), &mut a,
        0.1, phase_config(), || Ok(())).unwrap();
    let actual = guarded.sublimate_cell(CellId(1), &Ice, kinetics(), &mut b,
        0.1, phase_config(), || Ok(())).unwrap();
    assert_eq!(actual, expected);
    assert_eq!(guarded.cells(), bare.cells());
    assert!(actual.transferred_mass_kg > 0.0);
    assert!((guarded.energy_j() + 400.0 + actual.vapor_energy_j).abs() < 1e-7);
    let surfaces = [1, 2].map(|id| ChamberSurface { cell_id: CellId(id), kinetics: kinetics() });
    let (mut a, mut b) = (chamber(), chamber());
    let expected = bare.sublimate_chamber(&surfaces, &Ice, &mut a,
        0.1, chamber_config(), || Ok(())).unwrap();
    let actual = guarded.sublimate_chamber(&surfaces, &Ice, &mut b,
        0.1, chamber_config(), || Ok(())).unwrap();
    assert_eq!(actual.cell_ids, expected.cell_ids);
    for (got, want) in actual.exchange.surfaces.iter().zip(&expected.exchange.surfaces) {
        assert!((got.transferred_mass_kg - want.transferred_mass_kg).abs() < 1e-14);
        assert!((got.solid.specific_enthalpy_j_kg - want.solid.specific_enthalpy_j_kg).abs() < 1e-6);
    }
    let mass = |gas| match gas {
        VaporReservoir::Finite { vapor_mass_kg, .. } => vapor_mass_kg,
        _ => unreachable!(),
    };
    assert!((mass(a) - mass(b)).abs() < 1e-14);
    for (id, got) in guarded.cells() {
        let want = &bare.cells()[id];
        assert!((got.mass_kg - want.mass_kg).abs() < 1e-14);
        assert!((got.specific_h_j_kg - want.specific_h_j_kg).abs() < 1e-6);
    }
}

#[test]
fn heating_into_a_melting_plateau_refuses_the_whole_heat_phase_advance() {
    let source = curve(0.0);
    let solid = source.pure_solid_thermal_law().unwrap();
    let cells = vec![cell(1, 45_000.0, 400.0)];
    let bare = initial(&source, cells.clone());
    let guarded = initial(&solid, cells);
    let melted = bare.advance_to(0.1, config(), || Ok(())).unwrap().0;
    assert_eq!(melted.temperature_k(CellId(1)), Ok(273.0));
    assert!(melted.cells()[&CellId(1)].specific_h_j_kg > 46_000.0);
    let before = guarded.cells().clone();
    assert!(matches!(guarded.advance_sublimating_to(CellId(1), &Ice, kinetics(), vacuum(),
        0.1, config(), phase_config(), || Ok(())),
        Err(ThermalSublimationError::Thermal(DepositionError::MaterialDomain))));
    let surfaces = [ChamberSurface { cell_id: CellId(1), kinetics: kinetics() }];
    assert!(matches!(guarded.advance_chamber_sublimating_to(&surfaces, &Ice, chamber(),
        0.1, config(), chamber_config(), || Ok(())),
        Err(ThermalSublimationError::Thermal(DepositionError::MaterialDomain))));
    assert_eq!(guarded.cells(), &before);
    assert_eq!(guarded.time_s(), 0.0);
}

#[test]
fn external_exchange_cannot_override_the_material_phase_with_a_wide_caller_interval() {
    let source = curve(0.0);
    let solid = source.pure_solid_thermal_law().unwrap();
    let mut thermal = initial(&solid, vec![cell(1, 0.0, 0.0), cell(2, 0.0, 0.0)]);
    let before = thermal.cells().clone();
    let changes: Vec<_> = [(1, -1000.0), (2, 70_000.0)].into_iter().map(|(id, h)| SolidCellExchange {
        expected: before[&CellId(id)].clone(), interval_exposed_area_m2: 0.001,
        replacement: SolidSurfaceState { mass_kg: 0.019, specific_enthalpy_j_kg: h, area_m2: 0.001 },
        pure_solid_enthalpy_range_j_kg: [-140_000.0, 100_002.0],
    }).collect();
    assert_eq!(thermal.apply_solid_exchanges(&changes, || Ok(())), Err(DepositionError::MaterialDomain));
    assert_eq!(thermal.cells(), &before);
    thermal.apply_solid_exchanges(&changes[..1], || Ok(())).unwrap();
    assert_eq!(thermal.cells()[&CellId(1)].specific_h_j_kg, -1000.0);
    assert_eq!(thermal.cells()[&CellId(2)], before[&CellId(2)]);
}

#[test]
fn final_phase_cancellation_keeps_original_cells_and_vapor_reusable() {
    let source = curve(0.0);
    let solid = source.pure_solid_thermal_law().unwrap();
    let original = initial(&solid, vec![cell(1, 0.0, 0.0)]);
    let mut successful = original.clone();
    let mut gas = chamber();
    let mut polls = 0;
    let expected = successful.sublimate_cell(CellId(1), &Ice, kinetics(), &mut gas,
        0.1, phase_config(), || { polls += 1; Ok(()) }).unwrap();
    let mut attempted = original.clone();
    let mut original_gas = chamber();
    let mut calls = 0;
    assert!(matches!(attempted.sublimate_cell(CellId(1), &Ice, kinetics(), &mut original_gas,
        0.1, phase_config(), || {
            calls += 1;
            if calls == polls { Err(DepositionError::Cancelled) } else { Ok(()) }
        }), Err(ThermalSublimationError::Thermal(DepositionError::Cancelled))));
    assert_eq!(attempted.cells(), original.cells());
    assert_eq!(original_gas, chamber());
    let retried = attempted.sublimate_cell(CellId(1), &Ice, kinetics(), &mut original_gas,
        0.1, phase_config(), || Ok(())).unwrap();
    assert_eq!(retried, expected);
    assert_eq!(attempted.cells(), successful.cells());
}

#[test]
fn liquid_only_and_isolated_solid_endpoint_do_not_supply_a_solid_transport_domain() {
    let liquid = EquilibriumEnthalpyPhaseCurve::try_single_phase(ContentHash([0x92; 32]),
        SolidLiquidPhase::Liquid, curve(0.0).knots()[3..].to_vec()).unwrap();
    assert!(liquid.pure_solid_thermal_law().is_err());
    let mixed = EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([0x93; 32]),
        curve(0.0).knots()[2..].to_vec()).unwrap();
    assert!(mixed.pure_solid_thermal_law().is_err());
}
