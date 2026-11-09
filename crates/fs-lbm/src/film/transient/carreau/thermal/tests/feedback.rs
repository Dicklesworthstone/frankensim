//! G1/G3/G4: original caloric feedback, independent one-cell momentum/heat
//! balances, unchanged prescribed sampling and transactional source refusal.

use super::super::adaptive::{
    ThermalFilmAdaptiveControls, ThermalFilmAdaptiveError, ThermalFilmQuery,
};
use super::*;
use std::cell::{Cell, RefCell};
use std::sync::Arc;

const ID: CellId = CellId(17);
const CAPACITY: f64 = 0.02;
const COOLING: f64 = 0.01;

#[derive(Debug)]
struct Actuator;
impl Actuator {
    fn at(&self, temperature: f64, sign: f64) -> FilmDrive {
        let theta = temperature - 300.0;
        FilmDrive {
            acceleration_m_s2: sign * (0.2 + 0.1 * theta),
            wall_velocity_m_s: sign * (0.1 - 0.02 * theta),
            surface_shear_pa: sign * (1.0 - 0.2 * theta),
        }
    }
}

// Deliberately not Clone: the actual accepted source receipt must be moved
// through the native endpoint, not reconstructed from the last callback.
#[derive(Debug)]
struct Receipt {
    original: Arc<Actuator>,
    temperature: f64,
}

fn viscosity() -> ThermalCarreauYasuda {
    // The exact Newtonian-in-shear limit of the original thermal Carreau owner,
    // with rho=h=mu(300 K)=1 and its original positive Arrhenius shift.
    ThermalCarreauYasuda::new(
        CarreauYasuda::new(1.0, 1.0, 2.0, 0.45, 2.0, [0.0, 1e5]).unwrap(),
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
fn model(law: &Caloric) -> ThermalCarreauFilm<'_> {
    let cell = DepositedCell {
        id: ID,
        segment_id: 17,
        body_id: 1,
        material_id: MaterialId(7),
        mass_kg: 1.0,
        specific_h_j_kg: law.datum,
        geometry: CellGeometry {
            start_m: [0.0; 3],
            end_m: [0.0, 1.0, 0.0],
            volume_m3: 1.0,
            surface_area_m2: 2.0,
        },
        boundary: Boundary {
            ambient_k: 300.0,
            htc_w_m2_k: COOLING / 2.0,
            reservoir: None,
            source_w: 0.0,
        },
    };
    let thermal = DepositionSimulation::new(
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Birth(cell),
        }],
        &[(MaterialId(7), law as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 1,
            max_links: 1,
            max_events: 1,
            max_substeps: 4096,
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
    .0;
    let flow =
        TransientCarreauFilm::from_profile(viscosity(), 300.0, SPEC, vec![0.0], 0.0).unwrap();
    ThermalCarreauFilm::new(flow, thermal, &[ID], 1.0, [300.0; 2], || false).unwrap()
}
fn settings() -> ThermalFilmControls {
    ThermalFilmControls {
        maximum_endpoint_evaluations: 100_000,
        maximum_drive_evaluations: 100_000,
        maximum_momentum_material_evaluations: 50_000_000,
        ..controls()
    }
}
fn accuracy(error: f64) -> ThermalFilmAdaptiveControls {
    ThermalFilmAdaptiveControls {
        initial_dt_s: 0.08,
        maximum_dt_s: 0.08,
        minimum_dt_s: 1e-6,
        energy_absolute_tolerance_j: error,
        temperature_absolute_tolerance_k: error / CAPACITY,
        velocity_absolute_tolerance_m_s: error * 5.0,
        transport_absolute_tolerance_m2: error * 5.0,
        relative_tolerance: 0.0,
        window_energy_tolerance_j: 2e-8,
        ..ThermalFilmAdaptiveControls::default()
    }
}
fn sample(
    original: &Arc<Actuator>,
    q: ThermalFilmQuery,
    thermal: &DepositionSimulation<'_>,
    cells: &[CellId],
) -> Result<ThermalFilmDrive<Receipt>, DepositionError> {
    assert_eq!(cells, &[ID]);
    assert_eq!(thermal.time_s(), q.endpoint_s);
    let temperature = thermal.temperature_k(ID)?;
    let sign = if q.segment % 2 == 0 { 1.0 } else { -1.0 };
    Ok(ThermalFilmDrive {
        drive: original.at(temperature, sign),
        payload: Receipt {
            original: original.clone(),
            temperature,
        },
    })
}
fn heat(
    original: &Arc<Actuator>,
    thermal: &DepositionSimulation<'_>,
    cells: &[CellId],
    drive: &ThermalFilmDriveSample<Receipt>,
) -> Result<ThermalFilmExchange<()>, DepositionError> {
    assert!(Arc::ptr_eq(original, &drive.payload.original));
    assert_eq!(
        drive.payload.temperature.to_bits(),
        thermal.temperature_k(ID)?.to_bits()
    );
    assert_eq!(drive.thermal_time_s, Some(thermal.time_s()));
    insulated(thermal, cells)
}

// Exact one-cell midpoint momentum with endpoint-held source values. The
// two half-cell shear losses are physical heat, not the kinetic increment.
fn midpoint(old: f64, dt: f64, temperature: f64, drive: FilmDrive) -> (f64, f64, [f64; 3]) {
    let mu = viscosity()
        .evaluate(temperature, 0.0)
        .unwrap()
        .dynamic_viscosity_pa_s();
    let next = ((1.0 - dt * mu) * old
        + dt * (drive.acceleration_m_s2
            + drive.surface_shear_pa
            + 2.0 * mu * drive.wall_velocity_m_s))
        / (1.0 + dt * mu);
    let mean = 0.5 * old + 0.5 * next;
    let wall_stress = 2.0 * mu * (mean - drive.wall_velocity_m_s);
    let surface = mean + 0.5 * drive.surface_shear_pa / mu;
    let loss = dt
        * (2.0 * mu * (mean - drive.wall_velocity_m_s).powi(2)
            + 0.5 * drive.surface_shear_pa.powi(2) / mu);
    (
        next,
        loss,
        [
            dt * drive.acceleration_m_s2 * mean,
            -dt * wall_stress * drive.wall_velocity_m_s,
            dt * drive.surface_shear_pa * surface,
        ],
    )
}

#[test]
fn g1_actual_endpoint_feedback_closes_the_independent_heat_root_and_retains_source() {
    let law = Caloric {
        cp: CAPACITY,
        datum: 1000.0,
    };
    let original = Arc::new(Actuator);
    let mut coupled = model(&law);
    let calls = Cell::new(0);
    let observed = RefCell::new(Vec::new());
    let dt = 0.1;
    let report = coupled
        .advance_with_thermal_drive(
            dt,
            LIMITS,
            TransientCarreauControls::default(),
            settings(),
            &[0.0],
            |q, t, ids| {
                calls.set(calls.get() + 1);
                observed.borrow_mut().push(t.temperature_k(ID)?);
                sample(&original, q, t, ids)
            },
            |t, ids, drive| heat(&original, t, ids, drive),
            || false,
        )
        .unwrap();
    let residual = |theta: f64| {
        (CAPACITY + dt * COOLING) * theta
            - midpoint(0.0, dt, 300.0 + theta, original.at(300.0 + theta, 1.0)).1
    };
    let (mut lo, mut hi) = (0.0, 5.0);
    assert!(residual(lo) < 0.0 && residual(hi) > 0.0);
    for _ in 0..80 {
        let middle = 0.5 * lo + 0.5 * hi;
        if residual(middle) > 0.0 {
            hi = middle;
        } else {
            lo = middle;
        }
    }
    let temperature = 300.0 + 0.5 * lo + 0.5 * hi;
    let forcing = original.at(temperature, 1.0);
    let (velocity, viscous_heat, work) = midpoint(0.0, dt, temperature, forcing);
    near(
        coupled.thermal().temperature_k(ID).unwrap(),
        temperature,
        2e-8,
    );
    near(coupled.flow().velocities_m_s()[0], velocity, 1e-10);
    near(report.deposited_viscous_heat_j, viscous_heat, 1e-10);
    for (actual, expected) in [
        report.flow.flow.body_work_j_m2,
        report.flow.flow.wall_work_j_m2,
        report.flow.flow.surface_work_j_m2,
    ]
    .into_iter()
    .zip(work)
    {
        near(actual, expected, 1e-10);
    }
    near(
        report.other_heat_j,
        -dt * COOLING * (temperature - 300.0),
        1e-10,
    );
    assert!(report.combined_energy_residual_j.abs() < settings().combined_energy_tolerance_j);
    assert_eq!(report.drive_evaluations, calls.get());
    assert_eq!(report.drive_evaluations, report.endpoint_evaluations);
    assert!(calls.get() > 2);
    assert!(
        observed
            .borrow()
            .iter()
            .any(|&t| (t - temperature).abs() > 0.1)
    );
    assert!(Arc::ptr_eq(&original, &report.drive.payload.original));
    assert_eq!(
        report.drive.payload.temperature.to_bits(),
        coupled.thermal().temperature_k(ID).unwrap().to_bits()
    );
    assert_eq!(report.drive.query.start_s, 0.0);
    assert_eq!(report.drive.query.midpoint_s, 0.05);
    assert_eq!(report.drive.thermal_time_s, Some(dt));
    near(
        report.drive.drive.surface_shear_pa,
        forcing.surface_shear_pa,
        1e-10,
    );
    let mut lagged = model(&law);
    let lagged_report = lagged
        .advance(
            dt,
            original.at(300.0, 1.0),
            LIMITS,
            TransientCarreauControls::default(),
            settings(),
            &[0.0],
            insulated,
            || false,
        )
        .unwrap();
    assert!(lagged.thermal().temperature_k(ID).unwrap() - temperature > 0.5);
    assert!(lagged_report.deposited_viscous_heat_j - viscous_heat > 0.01);
    assert_eq!(lagged_report.drive_evaluations, 0);
    assert_eq!(lagged_report.drive.thermal_time_s, None);
}

// Independent continuous one-cell ODE, integrated only as a test reference.
// Four stages use the physical source/heat equations above, never production
// momentum or thermal solves. The authored sign jump is aligned explicitly.
fn reference(original: &Actuator, end: f64) -> [f64; 2] {
    let steps = 8192;
    let dt = end / steps as f64;
    let mut state = [0.0, 0.0];
    for i in 0..steps {
        let sign = if i < steps / 2 { 1.0 } else { -1.0 };
        let rates = |x: [f64; 2]| {
            let drive = original.at(300.0 + x[1], sign);
            let mu = viscosity()
                .evaluate(300.0 + x[1], 0.0)
                .unwrap()
                .dynamic_viscosity_pa_s();
            [
                drive.acceleration_m_s2 + drive.surface_shear_pa
                    - 2.0 * mu * (x[0] - drive.wall_velocity_m_s),
                (2.0 * mu * (x[0] - drive.wall_velocity_m_s).powi(2)
                    + 0.5 * drive.surface_shear_pa.powi(2) / mu
                    - COOLING * x[1])
                    / CAPACITY,
            ]
        };
        let a = rates(state);
        let b = rates(std::array::from_fn(|j| state[j] + 0.5 * dt * a[j]));
        let c = rates(std::array::from_fn(|j| state[j] + 0.5 * dt * b[j]));
        let d = rates(std::array::from_fn(|j| state[j] + dt * c[j]));
        state =
            std::array::from_fn(|j| state[j] + dt * (a[j] + 2.0 * b[j] + 2.0 * c[j] + d[j]) / 6.0);
    }
    state
}

#[test]
fn g1_adaptive_feedback_refines_original_transients_across_an_exact_drive_event() {
    let law = Caloric {
        cp: CAPACITY,
        datum: 1000.0,
    };
    let original = Arc::new(Actuator);
    let run = |error| {
        let mut model = model(&law);
        let calls = Cell::new(0);
        let report = model
            .advance_adaptive_with_thermal_drive(
                0.16,
                &[0.08],
                LIMITS,
                TransientCarreauControls::default(),
                settings(),
                accuracy(error),
                |q, t, ids| {
                    calls.set(calls.get() + 1);
                    sample(&original, q, t, ids)
                },
                |q, t, ids, drive| {
                    assert_eq!(q, drive.query);
                    heat(&original, t, ids, drive)
                },
                || false,
            )
            .unwrap();
        assert_eq!(report.drive_evaluations, calls.get());
        assert_eq!(report.drive_evaluations, report.endpoint_evaluations);
        assert_eq!(report.drive_samples.len(), report.accepted_steps);
        assert!(report.drive_evaluations > report.accepted_steps);
        let mut end = 0.0;
        for accepted in &report.drive_samples {
            assert_eq!(accepted.query.start_s, end);
            assert_eq!(accepted.thermal_time_s, Some(accepted.query.endpoint_s));
            assert!(Arc::ptr_eq(&accepted.payload.original, &original));
            assert!(accepted.query.start_s < 0.08 || accepted.query.segment == 1);
            assert!(accepted.query.endpoint_s <= 0.08 || accepted.query.segment == 1);
            let sign = if accepted.query.segment == 0 {
                1.0
            } else {
                -1.0
            };
            assert_eq!(
                accepted.drive,
                original.at(accepted.payload.temperature, sign)
            );
            end = accepted.query.endpoint_s;
        }
        assert_eq!(end, 0.16);
        assert!(
            report.combined_energy_residual_j.abs() < accuracy(error).window_energy_tolerance_j
        );
        (model, report)
    };
    let (coarse, coarse_report) = run(2e-5);
    let (fine, fine_report) = run(2e-7);
    let expected = reference(&original, 0.16);
    let difference = |model: &ThermalCarreauFilm<'_>| {
        (model.flow().velocities_m_s()[0] - expected[0]).abs()
            + CAPACITY * (model.thermal().temperature_k(ID).unwrap() - 300.0 - expected[1]).abs()
    };
    assert!(difference(&fine) < difference(&coarse) / 3.0);
    assert!(fine_report.accepted_steps > coarse_report.accepted_steps);
    near(fine.flow().velocities_m_s()[0], expected[0], 5e-5);
    near(
        fine.thermal().temperature_k(ID).unwrap(),
        300.0 + expected[1],
        0.003,
    );
}

#[test]
fn g3_prescribed_callbacks_stay_once_per_interval_and_share_the_original_solver() {
    let law = Caloric {
        cp: CAPACITY,
        datum: 1000.0,
    };
    let forcing = Actuator.at(300.0, 1.0);
    let mut prescribed = model(&law);
    let mut feedback = model(&law);
    let calls = Cell::new(0);
    let held = prescribed
        .advance_adaptive(
            0.04,
            &[],
            LIMITS,
            TransientCarreauControls::default(),
            settings(),
            accuracy(2e-5),
            |_| {
                calls.set(calls.get() + 1);
                Ok::<_, DepositionError>(forcing)
            },
            |_, t, ids| insulated(t, ids),
            || false,
        )
        .unwrap();
    let resolved = feedback
        .advance_adaptive_with_thermal_drive(
            0.04,
            &[],
            LIMITS,
            TransientCarreauControls::default(),
            settings(),
            accuracy(2e-5),
            |_, _, _| {
                Ok::<_, DepositionError>(ThermalFilmDrive {
                    drive: forcing,
                    payload: (),
                })
            },
            |_, t, ids, _| insulated(t, ids),
            || false,
        )
        .unwrap();
    assert_eq!(held.drive_evaluations, calls.get());
    assert_eq!(calls.get(), 3 * held.comparisons);
    assert!(held.endpoint_evaluations > calls.get());
    assert_eq!(resolved.drive_evaluations, resolved.endpoint_evaluations);
    assert_eq!(
        prescribed.flow().velocities_m_s(),
        feedback.flow().velocities_m_s()
    );
    assert_eq!(prescribed.thermal().cells(), feedback.thermal().cells());
    assert_eq!(held.deposited_cell_heat_j, resolved.deposited_cell_heat_j);
    assert!(
        held.drive_samples
            .iter()
            .all(|sample| sample.thermal_time_s.is_none())
    );
    assert!(
        resolved
            .drive_samples
            .iter()
            .all(|sample| sample.thermal_time_s == Some(sample.query.endpoint_s))
    );
}

#[test]
fn g4_feedback_query_budget_late_source_refusal_and_cancellation_publish_nothing() {
    let law = Caloric {
        cp: CAPACITY,
        datum: 1000.0,
    };
    let original = Arc::new(Actuator);
    let mut admitted = model(&law);
    let mut polls = 0;
    let good = admitted
        .advance_adaptive_with_thermal_drive(
            0.04,
            &[0.02],
            LIMITS,
            TransientCarreauControls::default(),
            settings(),
            accuracy(2e-5),
            |q, t, ids| sample(&original, q, t, ids),
            |_, t, ids, drive| heat(&original, t, ids, drive),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    let mut failed = model(&law);
    let before = failed.thermal().cells().clone();
    let count = Cell::new(0);
    let error = failed.advance_adaptive_with_thermal_drive(
        0.04,
        &[0.02],
        LIMITS,
        TransientCarreauControls::default(),
        ThermalFilmControls {
            maximum_drive_evaluations: good.drive_evaluations - 1,
            ..settings()
        },
        accuracy(2e-5),
        |q, t, ids| {
            count.set(count.get() + 1);
            sample(&original, q, t, ids)
        },
        |_, t, ids, drive| heat(&original, t, ids, drive),
        || false,
    );
    assert!(matches!(
        error,
        Err(ThermalFilmAdaptiveError::Film(
            ThermalFilmError::DriveBudget
        ))
    ));
    assert_eq!(count.get(), good.drive_evaluations - 1);
    assert_eq!(failed.thermal().cells(), &before);
    assert_eq!(failed.flow().velocities_m_s(), &[0.0]);
    let error = failed.advance_adaptive_with_thermal_drive(
        0.04,
        &[0.02],
        LIMITS,
        TransientCarreauControls::default(),
        settings(),
        accuracy(2e-5),
        |q, t, ids| {
            if q.start_s >= 0.02 {
                Err(DepositionError::InvalidInput("original actuator refuses"))
            } else {
                sample(&original, q, t, ids)
            }
        },
        |_, t, ids, drive| heat(&original, t, ids, drive),
        || false,
    );
    assert!(matches!(
        error,
        Err(ThermalFilmAdaptiveError::Film(ThermalFilmError::Drive(_)))
    ));
    assert_eq!(failed.thermal().cells(), &before);
    let mut cancelled_polls = 0;
    let error = failed.advance_adaptive_with_thermal_drive(
        0.04,
        &[0.02],
        LIMITS,
        TransientCarreauControls::default(),
        settings(),
        accuracy(2e-5),
        |q, t, ids| sample(&original, q, t, ids),
        |_, t, ids, drive| heat(&original, t, ids, drive),
        || {
            cancelled_polls += 1;
            cancelled_polls == polls
        },
    );
    assert!(matches!(
        error,
        Err(ThermalFilmAdaptiveError::Film(ThermalFilmError::Thermal(
            DepositionError::Cancelled
        )))
    ));
    assert_eq!(failed.thermal().cells(), &before);
    assert_eq!(failed.flow().accepted_time_s(), 0.0);
    let retry = failed
        .advance_adaptive_with_thermal_drive(
            0.04,
            &[0.02],
            LIMITS,
            TransientCarreauControls::default(),
            settings(),
            accuracy(2e-5),
            |q, t, ids| sample(&original, q, t, ids),
            |_, t, ids, drive| heat(&original, t, ids, drive),
            || false,
        )
        .unwrap();
    assert_eq!(failed.thermal().cells(), admitted.thermal().cells());
    assert_eq!(
        failed.flow().velocities_m_s(),
        admitted.flow().velocities_m_s()
    );
    assert_eq!(retry.drive_evaluations, good.drive_evaluations);
}
