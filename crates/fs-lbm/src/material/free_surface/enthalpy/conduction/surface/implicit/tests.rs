//! Analytic stiff capacities use the same production BE/Newton adapter.

use super::*;

struct Sensible;
impl ThermalLaw for Sensible {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if !(-100.0..=1e6).contains(&h) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(300.0 + h)
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0
    }
}
static LAW: Sensible = Sensible;

fn simulation(
    masses: &[f64],
    enthalpies: &[f64],
    link: Option<f64>,
) -> DepositionSimulation<'static> {
    let mut events: Vec<_> = masses
        .iter()
        .zip(enthalpies)
        .enumerate()
        .map(|(i, (&mass, &h))| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(i as u64),
                segment_id: i as u64,
                body_id: 0,
                material_id: MaterialId(0),
                mass_kg: mass,
                specific_h_j_kg: h,
                geometry: CellGeometry {
                    start_m: [0.0; 3],
                    end_m: [1.0, 0.0, 0.0],
                    volume_m3: mass,
                    surface_area_m2: 6.0,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        })
        .collect();
    if let Some(g) = link {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: 0,
                a: CellId(0),
                b: CellId(1),
                area_m2: 1.0,
                conductance_w_k: g,
            }),
        });
    }
    DepositionSimulation::new(
        events,
        &[(MaterialId(0), &LAW)],
        DepositionBudget {
            max_cells: 8,
            max_links: 8,
            max_events: 16,
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 1e-9,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}

fn controls() -> FreeSurfaceImplicitSolve {
    FreeSurfaceImplicitSolve {
        maximum_unknowns: 8,
        maximum_krylov_restart: 4,
        maximum_newton_iterations: 8,
        maximum_residual_evaluations: 256,
        difference_step_j_kg: 0.001,
        cell_energy_tolerance_j: 1e-8,
        newton: NewtonKrylovConfig {
            linear_restart: 4,
            absolute_tolerance: 1e-9,
            relative_tolerance: 1e-12,
            ..NewtonKrylovConfig::default()
        },
    }
}

fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:e} != {b:e}");
}

#[test]
fn g1_implicit_surface_tiny_capacities_use_actual_mass_without_explicit_cfl() {
    let masses = [1e-12, 1e-6, 1.0];
    let initial = simulation(&masses, &[10.0; 3], None);
    let dt = 0.2;
    let g = [1.0, 2.0, 3.0];
    let scales: Vec<_> = masses.iter().zip(g).map(|(&m, g)| m + dt * g).collect();
    let evaluations = Cell::new(0);
    let (next, temperatures, audit, nonlinear, maximum) = solve_interval(
        &initial,
        dt,
        &scales,
        controls(),
        1e-8,
        &evaluations,
        |endpoint| {
            let mut temperatures = Vec::new();
            let mut sources_w = BTreeMap::new();
            for (i, g) in g.into_iter().enumerate() {
                let t = endpoint.temperature_k(CellId(i as u64))?;
                temperatures.push(t);
                sources_w.insert(CellId(i as u64), g * (400.0 - t));
            }
            Ok(StepCoefficients {
                sources_w,
                conductances_w_k: BTreeMap::new(),
                payload: temperatures,
            })
        },
        || false,
    )
    .unwrap();
    assert!(nonlinear.converged && evaluations.get() < 24);
    assert!(maximum < 1e-8);
    for i in 0..3 {
        let expected = 310.0 + dt * g[i] * 90.0 / (masses[i] + dt * g[i]);
        near(temperatures[i], expected, 2e-9);
        near(
            next.temperature_k(CellId(i as u64)).unwrap(),
            expected,
            2e-9,
        );
        assert_eq!(next.cells()[&CellId(i as u64)].mass_kg, masses[i]);
        assert_eq!(initial.cells()[&CellId(i as u64)].specific_h_j_kg, 10.0);
    }
    near(audit.stored_energy_change_j, audit.external_heat_j, 1e-8);
    assert_eq!(audit.substeps, 1);
    assert!(dt > 1e11 * 0.5 * masses[0] / g[0]);
}

#[test]
fn g1_implicit_surface_internal_conduction_conserves_and_refines_first_order() {
    let mut errors = Vec::new();
    for count in [16, 32, 64] {
        let mut thermal = simulation(&[1.0, 2.0], &[10.0, 0.0], Some(2.0));
        let energy = thermal.energy_j();
        let dt = 0.5 / count as f64;
        let work = Cell::new(0);
        let mut config = controls();
        config.maximum_residual_evaluations = 4096;
        for step in 1..=count {
            let (next, (), audit, _, _) = solve_interval(
                &thermal,
                step as f64 * dt,
                &[1.0 + 2.0 * dt, 2.0 + 2.0 * dt],
                config,
                1e-8,
                &work,
                |_| {
                    Ok(StepCoefficients {
                        sources_w: BTreeMap::new(),
                        conductances_w_k: BTreeMap::new(),
                        payload: (),
                    })
                },
                || false,
            )
            .unwrap();
            assert_eq!(audit.external_heat_j, 0.0);
            thermal = next;
        }
        near(thermal.energy_j(), energy, 2e-8);
        let exact = 300.0 + 10.0 / 3.0 + 20.0 / 3.0 * (-1.5_f64).exp();
        errors.push((thermal.temperature_k(CellId(0)).unwrap() - exact).abs());
    }
    for ratio in errors.windows(2).map(|e| e[0] / e[1]) {
        assert!((1.8..2.1).contains(&ratio));
    }
}

#[test]
fn g4_implicit_surface_source_budget_cancel_and_unscaled_energy_gate_refuse() {
    let initial = simulation(&[1e-12], &[0.0], None);
    let resolve = |endpoint: &DepositionSimulation<'_>| -> Result<StepCoefficients<()>, Error> {
        Ok(StepCoefficients {
            sources_w: [(CellId(0), 400.0 - endpoint.temperature_k(CellId(0))?)].into(),
            conductances_w_k: BTreeMap::new(),
            payload: (),
        })
    };
    let mut budget = controls();
    budget.maximum_residual_evaluations = 1;
    let error = solve_interval(
        &initial,
        1.0,
        &[1.0 + 1e-12],
        budget,
        1e-9,
        &Cell::new(0),
        resolve,
        || false,
    )
    .err()
    .unwrap();
    assert!(matches!(
        error,
        Error::ImplicitSurface(FreeSurfaceImplicitFailure::EvaluationBudget)
    ));

    let error = solve_interval(
        &initial,
        1.0,
        &[1.0],
        controls(),
        1e-9,
        &Cell::new(0),
        |_| Err::<StepCoefficients<()>, _>(DepositionError::MaterialDomain.into()),
        || false,
    )
    .err()
    .unwrap();
    assert!(matches!(
        error,
        Error::Thermal(DepositionError::MaterialDomain)
    ));
    let polls = Cell::new(0);
    assert!(
        solve_interval(
            &initial,
            1.0,
            &[1.0],
            controls(),
            1e-9,
            &Cell::new(0),
            resolve,
            || {
                polls.set(polls.get() + 1);
                polls.get() > 8
            }
        )
        .is_err()
    );

    let mut loose = controls();
    loose.newton.absolute_tolerance = 1e9;
    let error = solve_interval(
        &initial,
        1.0,
        &[1.0],
        loose,
        1e-9,
        &Cell::new(0),
        resolve,
        || false,
    )
    .err()
    .unwrap();
    assert!(matches!(
        error,
        Error::Thermal(DepositionError::EnergyBalance(_))
    ));
    let mut oversized = controls();
    oversized.newton.linear_restart = 5;
    assert!(validate(oversized, 1).is_err());
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(initial.cells()[&CellId(0)].mass_kg, 1e-12);
    assert_eq!(initial.cells()[&CellId(0)].specific_h_j_kg, 0.0);
}
