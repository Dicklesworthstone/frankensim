//! Original native caloric/momentum coupling: independent energy and two-cell
//! thermal/ODE references. Synthetic laws are mathematical fixtures only.

use super::*;
use fs_conduction::enthalpy::{
    Boundary, CellGeometry, DepositedCell, DepositionBudget, Event, MaterialId, StepConfig,
    ThermalLaw, ThermalLink, TimedEvent,
};
use fs_material::fluid_transport::rheology::{
    CarreauYasuda, ThermalCarreauYasuda, ViscosityShiftMode,
};
use fs_material::visco::temperature::{TemperatureShift, TemperatureShiftModel};
use std::collections::BTreeMap;

#[path = "tests/adaptive.rs"]
mod adaptive;

#[path = "tests/feedback.rs"]
mod feedback;

#[path = "tests/slip.rs"]
mod slip;

const IDS: [CellId; 2] = [CellId(20), CellId(3)];
const SPEC: TransientCarreauFilmSpec = TransientCarreauFilmSpec {
    density_kg_m3: 1.0,
    thickness_m: 1.0,
    maximum_laminar_shear_rate_per_s: 20.0,
};
const LIMITS: FilmStepLimits = FilmStepLimits {
    maximum_speed_m_s: 20.0,
    maximum_velocity_residual_m_s: 1e-12,
};

struct Caloric {
    cp: f64,
    datum: f64,
}
impl ThermalLaw for Caloric {
    fn caloric_identity(&self) -> Option<[u8; 32]> {
        let mut id = [43; 32];
        id[..8].copy_from_slice(&self.cp.to_bits().to_le_bytes());
        id[8..16].copy_from_slice(&self.datum.to_bits().to_le_bytes());
        Some(id)
    }
    fn temperature_k(&self, energy: f64) -> Result<f64, DepositionError> {
        let temperature = 300.0 + (energy - self.datum) / self.cp;
        if temperature.is_finite() && (280.0..=360.0).contains(&temperature) {
            Ok(temperature)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / self.cp
    }
}

fn source() -> ThermalCarreauYasuda {
    ThermalCarreauYasuda::new(
        CarreauYasuda::new(2.0, 0.25, 2.0, 0.45, 2.0, [0.0, 1e5]).unwrap(),
        TemperatureShift::new(
            TemperatureShiftModel::Arrhenius {
                activation_energy_j_mol: 12_000.0,
            },
            300.0,
            [280.0, 360.0],
        )
        .unwrap(),
        ViscosityShiftMode::VerticalAndRate,
    )
}
fn traction() -> f64 {
    0.5 * source()
        .evaluate(300.0, 0.5)
        .unwrap()
        .dynamic_viscosity_pa_s()
}
fn drive() -> FilmDrive {
    FilmDrive {
        surface_shear_pa: traction(),
        ..FilmDrive::default()
    }
}
fn controls() -> ThermalFilmControls {
    ThermalFilmControls {
        maximum_thermal_cells: 4,
        maximum_endpoint_evaluations: 512,
        maximum_momentum_material_evaluations: 1_000_000,
        energy_scale_j: 0.01,
        finite_difference_step: 1e-5,
        cell_energy_tolerance_j: 2e-10,
        total_energy_tolerance_j: 5e-10,
        combined_energy_tolerance_j: 5e-10,
        solver: NewtonKrylovConfig {
            // Normalized residuals: stop below the independent joule gates,
            // while respecting representable changes to a nonzero h datum.
            absolute_tolerance: 1e-9,
            relative_tolerance: 1e-12,
            ..ThermalFilmControls::default().solver
        },
        ..ThermalFilmControls::default()
    }
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, tolerance {tolerance:e}"
    );
}

// Deliberately nonmonotone CellIds distinguish film order from the native
// thermal unknown order. Explicit areas/masses, not viscosity, define storage.
fn inventory<'a>(
    law: &'a Caloric,
    temperatures: [f64; 2],
    area: f64,
    contact_w_k: f64,
    bath_w_k: f64,
) -> DepositionSimulation<'a> {
    let exposed = area * if contact_w_k > 0.0 { 1.0 } else { 2.0 };
    let mut events: Vec<_> = (0..2)
        .map(|i| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: IDS[i],
                segment_id: i as u64,
                body_id: 1,
                material_id: MaterialId(7),
                mass_kg: 0.5 * area,
                specific_h_j_kg: law.datum + law.cp * (temperatures[i] - 300.0),
                geometry: CellGeometry {
                    start_m: [0.0, 0.5 * i as f64, 0.0],
                    end_m: [0.0, 0.5 * (i + 1) as f64, 0.0],
                    volume_m3: 0.5 * area,
                    surface_area_m2: 2.0 * area,
                },
                boundary: Boundary {
                    ambient_k: 280.0,
                    htc_w_m2_k: if i == 0 { bath_w_k / exposed } else { 0.0 },
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        })
        .collect();
    if contact_w_k > 0.0 {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: 91,
                a: IDS[0],
                b: IDS[1],
                area_m2: area,
                conductance_w_k: contact_w_k,
            }),
        });
    }
    DepositionSimulation::new(
        events,
        &[(MaterialId(7), law as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 4,
            max_links: 4,
            max_events: 8,
            max_substeps: 100,
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
fn coupled(law: &Caloric, area: f64) -> ThermalCarreauFilm<'_> {
    let flow =
        TransientCarreauFilm::from_profile(source(), 300.0, SPEC, vec![0.125, 0.375], 0.0).unwrap();
    ThermalCarreauFilm::new(
        flow,
        inventory(law, [300.0; 2], area, 0.0, 0.0),
        &IDS,
        area,
        [300.0; 2],
        || false,
    )
    .unwrap()
}
fn insulated(
    thermal: &DepositionSimulation<'_>,
    cells: &[CellId],
) -> Result<ThermalFilmExchange<()>, DepositionError> {
    Ok(ThermalFilmExchange {
        // Explicit zero-normal-gradient fluid trace closure at both boundaries.
        boundary_temperatures_k: [
            thermal.temperature_k(cells[0])?,
            thermal.temperature_k(cells[cells.len() - 1])?,
        ],
        heat: StepCoefficients {
            sources_w: BTreeMap::new(),
            conductances_w_k: BTreeMap::new(),
            payload: (),
        }
        .into(),
    })
}

#[test]
fn original_caloric_storage_and_momentum_close_one_isolated_heating_transaction() {
    let law = Caloric {
        cp: 0.05,
        datum: 1000.0,
    };
    let area = 2.5;
    let mut model = coupled(&law, area);
    let before = model.thermal().cells().clone();
    let velocity = model.flow().velocities_m_s().to_vec();
    let dt = 0.08;
    let report = model
        .advance(
            dt,
            drive(),
            LIMITS,
            TransientCarreauControls::default(),
            controls(),
            &[0.0; 2],
            insulated,
            || false,
        )
        .unwrap();
    let mut stored = 0.0;
    for &id in &IDS {
        let cell = &model.thermal().cells()[&id];
        let old = &before[&id];
        assert_eq!(cell.geometry, old.geometry);
        assert_eq!(cell.mass_kg, old.mass_kg);
        assert_eq!(cell.material_id, old.material_id);
        stored += cell.mass_kg * (cell.specific_h_j_kg - old.specific_h_j_kg);
    }
    let kinetic: f64 = velocity
        .iter()
        .zip(model.flow().velocities_m_s())
        .map(|(&old, &new)| 0.25 * area * (new - old) * (new + old))
        .sum();
    near(stored, report.deposited_viscous_heat_j, 5e-10);
    near(kinetic, report.kinetic_energy_change_j, 5e-13);
    near(stored + kinetic, report.mechanical_work_j, 5e-10);
    near(report.other_heat_j, 0.0, 0.0);
    near(report.heat_partition_residual_j, 0.0, 5e-13);
    near(report.combined_energy_residual_j, 0.0, 5e-10);
    assert_eq!(model.thermal().time_s(), dt);
    assert_eq!(model.flow().accepted_time_s(), dt);
    assert_eq!(model.film_cells(), &IDS);
    assert_eq!(model.wetted_area_m2(), area);
    assert_eq!(
        model
            .thermal()
            .material_caloric_identity(MaterialId(7))
            .unwrap(),
        law.caloric_identity(),
    );
    let t = IDS.map(|id| model.thermal().temperature_k(id).unwrap());
    assert!(t.iter().all(|&value| value > 300.0));
    for (face, expected) in [t[0], 0.5 * (t[0] + t[1]), t[1]].into_iter().enumerate() {
        let material = model.flow().material_at_face(face).unwrap();
        near(material.temperature_k(), expected, 2e-12);
        assert_eq!(material.law(), &source());
    }
    let mut frozen =
        TransientCarreauFilm::from_profile(source(), 300.0, SPEC, velocity, 0.0).unwrap();
    frozen
        .step(
            dt,
            drive(),
            LIMITS,
            TransientCarreauControls::default(),
            || false,
        )
        .unwrap();
    assert!(
        model
            .flow()
            .velocities_m_s()
            .iter()
            .zip(frozen.velocities_m_s())
            .any(|(&actual, &fixed)| (actual - fixed).abs() > 1e-5)
    );
}

#[test]
fn original_contact_and_bath_match_an_independent_two_by_two_backward_euler_step() {
    let law = Caloric {
        cp: 2.0,
        datum: 0.0,
    };
    let thermal = inventory(&law, [330.0, 310.0], 1.0, 2.0, 3.0);
    let flow = TransientCarreauFilm::new(source(), 300.0, SPEC, 2).unwrap();
    let mut model =
        ThermalCarreauFilm::new(flow, thermal, &IDS, 1.0, [280.0, 310.0], || false).unwrap();
    let dt = 0.1;
    let report = model
        .advance(
            dt,
            FilmDrive::default(),
            LIMITS,
            TransientCarreauControls::default(),
            controls(),
            &[0.0; 2],
            |candidate, cells| {
                let mut exchange = insulated(candidate, cells)?;
                // The original lower-cell conductance is to a held fluid-side
                // wall trace at 280 K; the free surface is insulated.
                exchange.boundary_temperatures_k[0] = 280.0;
                Ok::<_, DepositionError>(exchange)
            },
            || false,
        )
        .unwrap();
    // Each cell has heat capacity 1 J/K, internal conductance 2 W/K, and
    // cell zero has G=3 W/K to 280 K. Solve the original BE matrix directly.
    let (a, b, d) = (1.0 / dt + 2.0 + 3.0, 2.0, 1.0 / dt + 2.0);
    let rhs = [330.0 / dt + 3.0 * 280.0, 310.0 / dt];
    let expected = [
        (d * rhs[0] + b * rhs[1]) / (a * d - b * b),
        (b * rhs[0] + a * rhs[1]) / (a * d - b * b),
    ];
    for i in 0..2 {
        near(
            model.thermal().temperature_k(IDS[i]).unwrap(),
            expected[i],
            3e-9,
        );
    }
    near(report.other_heat_j, dt * 3.0 * (280.0 - expected[0]), 5e-9);
    near(report.deposited_viscous_heat_j, 0.0, 0.0);
    near(report.mechanical_work_j, 0.0, 0.0);
    near(report.combined_energy_residual_j, 0.0, 5e-10);
    assert!(model.flow().velocities_m_s().iter().all(|&u| u == 0.0));
}

// Independent continuous-in-time two-cell equations. Constitutive inversion
// uses plain bisection of the complete source evaluator, never PreparedCarreau
// or either production integrator. This is a temporal, not spatial, reference.
fn surface_rate(source: &ThermalCarreauYasuda, temperature: f64, stress: f64) -> f64 {
    let (mut low, mut high) = (0.0, 20.0);
    for _ in 0..80 {
        let midpoint = 0.5 * (low + high);
        if midpoint == low || midpoint == high {
            break;
        }
        let value = source
            .evaluate(temperature, midpoint)
            .unwrap()
            .dynamic_viscosity_pa_s()
            * midpoint;
        if value < stress {
            low = midpoint;
        } else {
            high = midpoint;
        }
    }
    0.5 * (low + high)
}
fn ode(source: &ThermalCarreauYasuda, state: [f64; 4]) -> [f64; 4] {
    let [u, v, t0, t1] = state;
    let gamma = [4.0 * u, 2.0 * (v - u)];
    let temperatures = [t0, 0.5 * (t0 + t1)];
    let stress: [f64; 2] = std::array::from_fn(|i| {
        source
            .evaluate(temperatures[i], gamma[i].abs())
            .unwrap()
            .dynamic_viscosity_pa_s()
            * gamma[i]
    });
    let wall_power = stress[0] * u;
    let inner_power = stress[1] * (v - u);
    let surface_power = traction() * 0.25 * surface_rate(source, t1, traction());
    [
        2.0 * (stress[1] - stress[0]),
        2.0 * (traction() - stress[1]),
        (wall_power + 0.5 * inner_power) / 0.025,
        (0.5 * inner_power + surface_power) / 0.025,
    ]
}
fn rk4_reference(intervals: usize) -> [f64; 4] {
    let source = source();
    let mut state = [0.125, 0.375, 300.0, 300.0];
    let h = 0.16 / intervals as f64;
    let shifted =
        |a: [f64; 4], b: [f64; 4], scale: f64| std::array::from_fn(|i| a[i] + scale * b[i]);
    for _ in 0..intervals {
        let a = ode(&source, state);
        let b = ode(&source, shifted(state, a, 0.5 * h));
        let c = ode(&source, shifted(state, b, 0.5 * h));
        let d = ode(&source, shifted(state, c, h));
        state =
            std::array::from_fn(|i| state[i] + h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0);
    }
    state
}

#[test]
fn endpoint_temperature_composition_has_first_order_coupled_time_convergence() {
    let reference = rk4_reference(2048);
    let comparison = rk4_reference(1024);
    for i in 0..4 {
        near(reference[i], comparison[i], 2e-9);
    }
    let law = Caloric {
        cp: 0.05,
        datum: 0.0,
    };
    let mut prior = f64::INFINITY;
    for intervals in [4, 8, 16] {
        let mut model = coupled(&law, 1.0);
        for step in 1..=intervals {
            model
                .advance(
                    0.16 * step as f64 / intervals as f64,
                    drive(),
                    LIMITS,
                    TransientCarreauControls::default(),
                    controls(),
                    &[0.0; 2],
                    insulated,
                    || false,
                )
                .unwrap();
        }
        let error = (0..2)
            .map(|i| (model.thermal().temperature_k(IDS[i]).unwrap() - reference[i + 2]).abs())
            .fold(0.0_f64, f64::max);
        assert!(
            error < 0.57 * prior,
            "intervals={intervals}, error={error:e}, previous={prior:e}"
        );
        prior = error;
    }
    assert!(prior < 4e-4, "{prior:e}");
}

#[test]
fn mapping_source_budgets_and_late_cancellation_preserve_both_original_states() {
    let law = Caloric {
        cp: 0.05,
        datum: 0.0,
    };
    let flow = || TransientCarreauFilm::new(source(), 300.0, SPEC, 2).unwrap();
    assert!(
        ThermalCarreauFilm::new(
            flow(),
            inventory(&law, [300.0; 2], 1.0, 0.0, 0.0),
            &IDS,
            2.0,
            [300.0; 2],
            || false,
        )
        .is_err()
    );
    assert!(
        ThermalCarreauFilm::new(
            flow(),
            inventory(&law, [300.0; 2], 1.0, 0.0, 0.0),
            &[IDS[1], IDS[0]],
            1.0,
            [300.0; 2],
            || false,
        )
        .is_err()
    );
    let mut victim = coupled(&law, 1.0);
    let before = victim.thermal().cells().clone();
    let velocity = victim.flow().velocities_m_s().to_vec();
    let refused = victim.advance(
        0.04,
        drive(),
        LIMITS,
        TransientCarreauControls::default(),
        controls(),
        &[0.0; 2],
        |_, _| Err::<ThermalFilmExchange<()>, _>(DepositionError::MaterialDomain),
        || false,
    );
    assert!(matches!(
        refused,
        Err(ThermalFilmError::Source(DepositionError::MaterialDomain))
    ));
    let queries = Cell::new(0);
    assert!(
        victim
            .advance(
                0.04,
                drive(),
                LIMITS,
                TransientCarreauControls::default(),
                controls(),
                &[0.0; 2],
                |thermal, cells| {
                    queries.set(queries.get() + 1);
                    let mut exchange = insulated(thermal, cells)?;
                    if queries.get() > 1 {
                        exchange.boundary_temperatures_k[1] = 361.0;
                    }
                    Ok::<_, DepositionError>(exchange)
                },
                || false,
            )
            .is_err()
    );
    assert!(queries.get() > 1);
    assert!(
        victim
            .advance(
                0.04,
                drive(),
                LIMITS,
                TransientCarreauControls::default(),
                ThermalFilmControls {
                    maximum_endpoint_evaluations: 1,
                    ..controls()
                },
                &[0.0; 2],
                insulated,
                || false,
            )
            .is_err()
    );
    assert_eq!(victim.thermal().cells(), &before);
    assert_eq!(victim.flow().velocities_m_s(), velocity);
    assert_eq!(victim.thermal().time_s(), 0.0);
    assert_eq!(victim.flow().accepted_time_s(), 0.0);

    let mut reference = coupled(&law, 1.0);
    let mut polls = 0;
    let expected = reference
        .advance(
            0.04,
            drive(),
            LIMITS,
            TransientCarreauControls::default(),
            controls(),
            &[0.0; 2],
            insulated,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut cancelled = coupled(&law, 1.0);
    let mut calls = 0;
    assert!(
        cancelled
            .advance(
                0.04,
                drive(),
                LIMITS,
                TransientCarreauControls::default(),
                controls(),
                &[0.0; 2],
                insulated,
                || {
                    calls += 1;
                    calls == polls
                },
            )
            .is_err()
    );
    assert_eq!(calls, polls);
    assert_eq!(cancelled.thermal().cells(), &before);
    assert_eq!(cancelled.flow().velocities_m_s(), velocity);
    assert_eq!(cancelled.thermal().time_s(), 0.0);
    assert_eq!(cancelled.flow().accepted_time_s(), 0.0);
    let accepted = cancelled
        .advance(
            0.04,
            drive(),
            LIMITS,
            TransientCarreauControls::default(),
            controls(),
            &[0.0; 2],
            insulated,
            || false,
        )
        .unwrap();
    assert_eq!(accepted.flow, expected.flow);
    assert_eq!(
        accepted.deposited_viscous_heat_j,
        expected.deposited_viscous_heat_j
    );
    assert_eq!(cancelled.thermal().cells(), reference.thermal().cells());
    assert_eq!(
        cancelled.flow().velocities_m_s(),
        reference.flow().velocities_m_s()
    );
}

