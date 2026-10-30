use super::*;
use crate::enthalpy::pyrolysis::{
    Arrhenius, CondensedCaloric, PyrolysisMaterial, VolatileEnthalpy,
};
use crate::enthalpy::{
    CellGeometry, DepositedCell, DepositionBudget, DepositionSimulation, Event, MaterialId,
    StepConfig, ThermalLaw, ThermalLink, TimedEvent,
};

struct Caloric {
    cp: f64,
    offset: f64,
    quadratic: f64,
}
impl ThermalLaw for Caloric {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let relative = h - self.offset;
        let t = if self.quadratic == 0.0 {
            300.0 + relative / self.cp
        } else {
            let b = self.cp + 300.0 * self.quadratic;
            // Positive root in temperature change, avoiding subtraction of
            // two large nearly equal roots around the reference temperature.
            300.0 + 2.0 * relative / (b + (b * b + 2.0 * self.quadratic * relative).sqrt())
        };
        if (200.0..=1200.0).contains(&t) {
            Ok(t)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / (self.cp + 200.0 * self.quadratic)
    }
}
impl CondensedCaloric for Caloric {
    fn specific_enthalpy_j_kg(&self, t: f64) -> Result<f64, DepositionError> {
        if (200.0..=1200.0).contains(&t) {
            Ok(self.offset + (t - 300.0) * (self.cp + 0.5 * self.quadratic * (t + 300.0)))
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
}
static CALORIC: Caloric = Caloric {
    cp: 2000.0,
    offset: 0.0,
    quadratic: 0.0,
};
const AREA: f64 = 0.01;
const DENSITY: f64 = 800.0;
const LENGTH: f64 = 0.008;
const CONDUCTIVITY: f64 = 0.2;

fn boundary(power: f64) -> Boundary {
    Boundary {
        ambient_k: 300.0,
        htc_w_m2_k: 0.0,
        reservoir: None,
        source_w: power,
    }
}
fn fixture<'a>(
    law: &'a dyn CondensedCaloric,
    gas: &'a dyn VolatileEnthalpy,
    temperatures: &[f64],
    weights: &[f64],
    rate: f64,
) -> RecedingSlab<'a> {
    fixture_with_components(law, gas, temperatures, weights, rate, None)
}

fn fixture_with_components<'a>(
    law: &'a dyn CondensedCaloric,
    gas: &'a dyn VolatileEnthalpy,
    temperatures: &[f64],
    weights: &[f64],
    rate: f64,
    char: Option<(&'a dyn CondensedCaloric, f64, CharProperties<'a>)>,
) -> RecedingSlab<'a> {
    let n = temperatures.len();
    assert_eq!(n, weights.len());
    let sum: f64 = weights.iter().sum();
    let mut events = Vec::new();
    let mut widths = Vec::new();
    let mut x = 0.0;
    for (i, (&t, &weight)) in temperatures.iter().zip(weights).enumerate() {
        let end = if i + 1 == n {
            LENGTH
        } else {
            x + LENGTH * weight / sum
        };
        let dx = end - x;
        widths.push(dx);
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(i as u64),
                segment_id: 1,
                body_id: 1,
                material_id: MaterialId(1),
                mass_kg: DENSITY * (AREA * dx),
                specific_h_j_kg: law.specific_enthalpy_j_kg(t).unwrap(),
                geometry: CellGeometry {
                    start_m: [x, 0.0, 0.0],
                    end_m: [end, 0.0, 0.0],
                    volume_m3: AREA * dx,
                    surface_area_m2: 2.0 * AREA,
                },
                boundary: boundary(0.0),
            }),
        });
        x = end;
    }
    for i in 1..n {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: i as u64,
                a: CellId(i as u64 - 1),
                b: CellId(i as u64),
                area_m2: AREA,
                conductance_w_k: series_conductance_w_k(
                    AREA,
                    [0.5 * widths[i - 1], 0.5 * widths[i]],
                    [CONDUCTIVITY; 2],
                )
                .unwrap(),
            }),
        });
    }
    let thermal = DepositionSimulation::new(
        events,
        &[(MaterialId(1), law)],
        DepositionBudget {
            max_cells: 128,
            max_links: 127,
            max_events: 255,
            max_substeps: 10000,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 0.1,
            energy_tolerance_j: 1e-8,
        },
        || Ok(()),
    )
    .unwrap()
    .0;
    let solid = PyrolysisSimulation::new(
        thermal,
        &[PyrolysisMaterial {
            material_id: MaterialId(1),
            virgin: law,
            char_caloric: char.map(|c| c.0),
            volatile: gas,
            char_mass_yield: char.map_or(0.0, |c| c.1),
            kinetics: Arrhenius {
                prefactor_s_inv: rate,
                activation_energy_j_mol: 0.0,
            },
            temperature_range_k: [200.0, 1200.0],
            source: "independent manufactured shrinking slab",
        }],
        || Ok(()),
    )
    .unwrap();
    let geometry = SlabGeometry {
        area_m2: AREA,
        density_kg_m3: DENSITY,
        conductivity_w_m_k: CONDUCTIVITY,
        source: "constant incompressible manufactured condensed solid",
    };
    match char {
        None => RecedingSlab::new(solid, geometry, || Ok(())),
        Some(c) => RecedingSlab::with_char(solid, geometry, c.2, || Ok(())),
    }
    .unwrap()
}

fn approx(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.17e} != {b:.17e} (tol {tolerance:e})"
    );
}

#[test]
fn isothermal_recession_matches_independent_inventory_and_fixed_back_geometry() {
    let gas = |t| CALORIC.specific_enthalpy_j_kg(t);
    let initial = fixture(&CALORIC, &gas, &[600.0; 8], &[1.0; 8], 2.0);
    let mut work = initial.state().implicit_workspace().unwrap();
    let controls = PyrolysisControls {
        maximum_dt_s: 0.1,
        ..Default::default()
    };
    let (next, report) = initial
        .advance_implicit_to(
            0.5,
            controls,
            ImplicitSolve::default(),
            SlabSolve::default(),
            &mut work,
            || Ok(()),
        )
        .unwrap();
    let expected_mass = initial.mass_kg() / 1.2_f64.powi(5);
    approx(next.mass_kg(), expected_mass, 1e-16);
    approx(next.thickness_m(), expected_mass / (DENSITY * AREA), 1e-16);
    assert_eq!(next.back_m().to_bits(), initial.back_m().to_bits());
    approx(next.recession_m(), LENGTH - next.thickness_m(), 1e-16);
    approx(
        report.pyrolysis.volatile_mass_kg,
        initial.mass_kg() - expected_mass,
        1e-16,
    );
    approx(
        next.energy_j() + report.pyrolysis.volatile_enthalpy_j,
        initial.energy_j(),
        1e-8,
    );
    assert_eq!(report.geometry_steps, 5);
    assert_eq!(report.rejected_geometry_attempts, 0);
    assert_eq!(report.remaps, 5);
    for cell in next.state().cells().values() {
        approx(next.state().temperature_k(cell.id).unwrap(), 600.0, 2e-10);
        approx(cell.mass_kg, expected_mass / 8.0, 1e-16);
        approx(cell.geometry.volume_m3 * DENSITY, cell.mass_kg, 1e-17);
    }
    let dx = next.thickness_m() / 8.0;
    approx(
        next.front_half_cell_conductance_w_k().unwrap(),
        2.0 * CONDUCTIVITY * AREA / dx,
        1e-11,
    );
    approx(report.remap_enthalpy_redistribution_j, 0.0, 1e-9);
}

#[test]
fn overlap_remap_preserves_extensive_energy_and_inverts_nonlinear_calorics() {
    let law = Caloric {
        cp: 2000.0,
        offset: -4e6,
        quadratic: 2.0,
    };
    let gas = |t| law.specific_enthalpy_j_kg(t);
    let mut slab = fixture(&law, &gas, &[600.0, 300.0], &[1.0, 1.0], 0.0);
    // An admissible, separately manufactured nonuniform retained inventory.
    // With two donors the reconstruction is constant, giving an exact overlap
    // oracle independent of the integrator and reconstruction code.
    slab.solid
        .thermal
        .cells
        .get_mut(&CellId(0))
        .unwrap()
        .mass_kg *= 0.5;
    slab.solid
        .inventories
        .get_mut(&CellId(0))
        .unwrap()
        .virgin_mass_kg *= 0.5;
    let m = slab.mass_kg();
    let energy = slab.energy_j();
    let hot_h = law.specific_enthalpy_j_kg(600.0).unwrap();
    let cold_h = law.specific_enthalpy_j_kg(300.0).unwrap();
    let expected_h = cold_h + (2.0 / 3.0) * (hot_h - cold_h);
    let expected_t = law.temperature_k(expected_h).unwrap();
    let audit = slab.rezone(SlabSolve::default(), &mut || Ok(())).unwrap();
    approx(slab.mass_kg(), m, 1e-16);
    approx(slab.energy_j(), energy, 1e-8);
    approx(
        slab.state().cells()[&CellId(0)].specific_h_j_kg,
        expected_h,
        1e-9,
    );
    approx(
        slab.state().temperature_k(CellId(0)).unwrap(),
        expected_t,
        1e-10,
    );
    assert!(
        (expected_t - 500.0).abs() > 5.0,
        "must not average temperatures"
    );
    approx(slab.state().temperature_k(CellId(1)).unwrap(), 300.0, 1e-10);
    approx(slab.front_m(), 0.002, 1e-16);
    assert!(audit.redistribution_j > 0.0);
    assert!(audit.mass_kg.abs() < 1e-16 && audit.energy_j.abs() < 1e-8);
}

#[test]
fn monotone_mass_coordinate_reconstruction_refines_and_preserves_energy_datum() {
    let evaluate = |n: usize, offset: f64| {
        let law = Caloric {
            cp: 2000.0,
            offset,
            quadratic: 0.0,
        };
        let gas = |t| law.specific_enthalpy_j_kg(t);
        let weights: Vec<_> = (0..n).map(|i| if i % 2 == 0 { 1.0 } else { 2.0 }).collect();
        let total: f64 = weights.iter().sum();
        let mut x = 0.0;
        let temperatures: Vec<_> = weights
            .iter()
            .map(|w| {
                let end = x + w / total;
                let t = 300.0 + 100.0 * (x * x + x * end + end * end) / 3.0;
                x = end;
                t
            })
            .collect();
        let mut slab = fixture(&law, &gas, &temperatures, &weights, 0.0);
        let energy = slab.energy_j();
        let mass = slab.mass_kg();
        let audit = slab.rezone(SlabSolve::default(), &mut || Ok(())).unwrap();
        let mut error = 0.0;
        let mut result = Vec::new();
        for (i, &id) in slab.ids.iter().enumerate() {
            let a = i as f64 / n as f64;
            let b = (i + 1) as f64 / n as f64;
            let exact_t = 300.0 + 100.0 * (a * a + a * b + b * b) / 3.0;
            let actual_t = slab.state().temperature_k(id).unwrap();
            assert!(actual_t >= temperatures[0] - 1e-9 && actual_t <= temperatures[n - 1] + 1e-9);
            error += (actual_t - exact_t).abs() / n as f64;
            result.push(actual_t);
        }
        assert!(result.windows(2).all(|w| w[0] <= w[1]));
        approx(slab.mass_kg(), mass, 1e-16);
        approx(slab.energy_j(), energy, 2e-8);
        assert!(audit.energy_j.abs() < 2e-8);
        (error, result)
    };
    let errors = [8, 16, 32].map(|n| evaluate(n, 0.0).0);
    assert!(
        errors[0] / errors[1] > 3.0 && errors[1] / errors[2] > 3.0,
        "{errors:?}"
    );
    let base = evaluate(16, 0.0).1;
    let shifted = evaluate(16, 1e8).1;
    for (a, b) in base.into_iter().zip(shifted) {
        approx(a, b, 1e-8);
    }
}

#[test]
fn reaction_disabled_retains_exact_original_implicit_state_and_accounts() {
    let gas = |t| CALORIC.specific_enthalpy_j_kg(t);
    let mut slab = fixture(
        &CALORIC,
        &gas,
        &[650.0, 500.0, 400.0, 300.0],
        &[1.0; 4],
        0.0,
    );
    slab.set_boundary(slab.front_cell_id(), boundary(20.0), || Ok(()))
        .unwrap();
    let controls = PyrolysisControls {
        maximum_dt_s: 0.2,
        ..Default::default()
    };
    let mut work = slab.state().implicit_workspace().unwrap();
    let (original, original_report) = slab
        .state()
        .advance_implicit_to(
            0.2,
            controls,
            ImplicitSolve::default(),
            &mut work,
            || Ok(()),
        )
        .unwrap();
    let mut work = slab.state().implicit_workspace().unwrap();
    let (next, report) = slab
        .advance_implicit_to(
            0.2,
            controls,
            ImplicitSolve::default(),
            SlabSolve::default(),
            &mut work,
            || Ok(()),
        )
        .unwrap();
    assert_eq!(next.state().cells(), original.cells());
    assert_eq!(next.state().inventories(), original.inventories());
    assert_eq!(report.pyrolysis, original_report.pyrolysis);
    assert_eq!(report.remaps, 0);
    assert_eq!(next.recession_m(), 0.0);
}

#[test]
fn complete_geometry_retries_do_not_accumulate_rejected_mass_or_heat() {
    let gas = |t| CALORIC.specific_enthalpy_j_kg(t);
    let slab = fixture(&CALORIC, &gas, &[600.0; 2], &[1.0; 2], 9.0);
    let mut work = slab.state().implicit_workspace().unwrap();
    let controls = PyrolysisControls {
        maximum_dt_s: 0.5,
        ..Default::default()
    };
    let (next, report) = slab
        .advance_implicit_to(
            0.5,
            controls,
            ImplicitSolve::default(),
            SlabSolve::default(),
            &mut work,
            || Ok(()),
        )
        .unwrap();
    assert!(report.rejected_geometry_attempts > 0 && report.geometry_steps > 1);
    assert!(report.maximum_cell_mass_loss_fraction <= 0.25);
    approx(
        next.mass_kg() + report.pyrolysis.volatile_mass_kg,
        slab.mass_kg(),
        1e-15,
    );
    approx(
        next.energy_j() + report.pyrolysis.volatile_enthalpy_j,
        slab.energy_j(),
        1e-7,
    );
    assert_eq!(next.time_s(), 0.5);
    assert_eq!(slab.time_s(), 0.0);
    assert!(next.thickness_m() < 0.0002);
    let refusal = slab.advance_implicit_to(
        0.5,
        controls,
        ImplicitSolve::default(),
        SlabSolve {
            maximum_halvings: 0,
            ..Default::default()
        },
        &mut work,
        || Ok(()),
    );
    assert!(matches!(refusal, Err(SlabError::GeometrySteps)));
    assert_eq!(slab.time_s(), 0.0);
}

#[test]
fn cancellation_after_geometry_preserves_state_and_deterministic_workspace_reuse() {
    let gas = |t| CALORIC.specific_enthalpy_j_kg(t);
    let slab = fixture(&CALORIC, &gas, &[650.0, 500.0, 300.0], &[1.0; 3], 0.4);
    let controls = PyrolysisControls {
        maximum_dt_s: 0.02,
        ..Default::default()
    };
    let mut work = slab.state().implicit_workspace().unwrap();
    let mut calls = 0usize;
    let reference = slab
        .advance_implicit_to(
            0.02,
            controls,
            ImplicitSolve::default(),
            SlabSolve::default(),
            &mut work,
            || {
                calls += 1;
                Ok(())
            },
        )
        .unwrap();
    let snapshot = slab.state().cells().clone();
    let mut count = 0usize;
    let result = slab.advance_implicit_to(
        0.02,
        controls,
        ImplicitSolve::default(),
        SlabSolve::default(),
        &mut work,
        || {
            count += 1;
            if count >= calls - 1 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        },
    );
    assert!(matches!(result, Err(SlabError::Cancelled)));
    assert_eq!(slab.state().cells(), &snapshot);
    assert_eq!(slab.time_s(), 0.0);
    let replay = slab
        .advance_implicit_to(
            0.02,
            controls,
            ImplicitSolve::default(),
            SlabSolve::default(),
            &mut work,
            || Ok(()),
        )
        .unwrap();
    assert_eq!(reference.0.state().cells(), replay.0.state().cells());
    assert_eq!(reference.1, replay.1);
}

#[test]
fn source_geometry_char_and_work_limits_refuse_without_partial_publication() {
    let gas = |t| CALORIC.specific_enthalpy_j_kg(t);
    let slab = fixture(&CALORIC, &gas, &[600.0; 3], &[1.0; 3], 0.1);
    let mut invalid = slab.solid.clone();
    invalid
        .materials
        .get_mut(&MaterialId(1))
        .unwrap()
        .char_mass_yield = 0.1;
    assert!(RecedingSlab::new(invalid, slab.geometry(), || Ok(())).is_err());
    let invalid = slab.solid.clone();
    assert!(
        RecedingSlab::new(
            invalid,
            SlabGeometry {
                density_kg_m3: 2.0 * DENSITY,
                ..slab.geometry()
            },
            || Ok(())
        )
        .is_err()
    );
    let mut invalid = slab.solid.clone();
    invalid
        .thermal
        .cells
        .get_mut(&CellId(1))
        .unwrap()
        .boundary
        .source_w = 1.0;
    assert!(RecedingSlab::new(invalid, slab.geometry(), || Ok(())).is_err());
    let mut work = slab.state().implicit_workspace().unwrap();
    let result = slab.advance_implicit_to(
        1.0,
        PyrolysisControls::default(),
        ImplicitSolve::default(),
        SlabSolve {
            maximum_geometry_steps: 1,
            ..Default::default()
        },
        &mut work,
        || Ok(()),
    );
    assert!(matches!(result, Err(SlabError::GeometrySteps)));
    assert_eq!(slab.time_s(), 0.0);
    let result = slab.advance_implicit_to(
        0.0,
        PyrolysisControls {
            cell_energy_tolerance_j: f64::NAN,
            ..Default::default()
        },
        ImplicitSolve::default(),
        SlabSolve::default(),
        &mut work,
        || Ok(()),
    );
    assert!(matches!(
        result,
        Err(SlabError::Thermal(DepositionError::InvalidInput(_)))
    ));
}

#[path = "charring_tests.rs"]
mod charring;
