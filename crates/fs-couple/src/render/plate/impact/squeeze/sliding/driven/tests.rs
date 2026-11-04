use super::*;
use crate::render::plate::impact::squeeze::sliding::{ShearedFilm, WallVelocityPorts};
use crate::render::plate::impact::squeeze::{
    FilmCell, FilmChannel, FilmLimits, GapPort, ResistiveFilm,
};
use crate::render::plate::impact::{ImpactBody, ImpactConfig};
use fs_exec::CancelGate;

const MASS: f64 = 0.5;
const GAP: f64 = 0.001;
const WIDTH: f64 = 0.02;
const LENGTH: f64 = 0.03;
const MU: f64 = 0.2;
const PRESSURE: f64 = 4000.0;
const DT: f64 = 0.01;
const LAUNCH: f64 = 0.03;

fn bare(speed: f64) -> ImpactSystem {
    let (body, _) = ImpactBody::free_mass(MASS, 0.0, speed).unwrap();
    ImpactSystem::new(
        vec![body],
        vec![],
        vec![],
        vec![],
        ImpactConfig {
            dt_s: DT,
            max_steps: 128,
            maximum_energy_j: 1.0,
            energy_absolute_tolerance_j: 1e-12,
            energy_relative_tolerance: 1e-8,
            maximum_generalized_force: 1.0,
        },
    )
    .unwrap()
}

fn graph() -> ResistiveFilm {
    let gap = GapPort {
        reference_m: GAP,
        closure: vec![0.0],
    };
    ResistiveFilm::new(
        vec![FilmCell {
            area_m2: WIDTH * LENGTH,
            gap: gap.clone(),
        }],
        (0..2)
            .map(|_| FilmChannel {
                from: 0,
                to: None,
                width_m: WIDTH,
                length_m: LENGTH,
                gap: gap.clone(),
            })
            .collect(),
        1,
        MU,
        FilmLimits {
            minimum_cell_gap_m: 1e-6,
            maximum_gap_m: 0.01,
            maximum_pressure_pa: 1e6,
        },
    )
    .unwrap()
}

fn sheared() -> ShearedFilm {
    // The high-pressure edge is oriented against positive physical motion.
    // The low-pressure edge is oriented with it, so entrainment balances at
    // the central cell while the pressure drop drives the moving wall.
    let walls = [-1.0, 1.0]
        .into_iter()
        .map(|sign| WallVelocityPorts {
            lower: vec![sign / MASS.sqrt()],
            upper: vec![0.0],
        })
        .collect();
    ShearedFilm::new(graph(), walls).unwrap()
}

fn driven(pressure: f64) -> PressureDrivenShearedFilm {
    PressureDrivenShearedFilm::new(sheared(), &[Some(pressure), Some(0.0)]).unwrap()
}

fn system(speed: f64, pressure: f64) -> ImpactSystem {
    bare(speed)
        .with_driven_sliding_film(driven(pressure))
        .unwrap()
}

fn speed(system: &ImpactSystem) -> f64 {
    system.state()[1] / MASS.sqrt()
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}; tolerance {tolerance}"
    );
}

fn oracle(before: f64, pressure: f64) -> (f64, f64, f64) {
    // Independent two-channel Couette/Poiseuille solution, with p_cell=P/2.
    // Quadratic free-mass storage makes the original DG step implicit midpoint.
    let conductance = WIDTH * GAP.powi(3) / (12.0 * MU * LENGTH);
    let drag = 2.0 * MU * WIDTH * LENGTH / GAP;
    let drive = WIDTH * GAP * pressure / 2.0;
    let after = ((MASS - 0.5 * DT * drag) * before + DT * drive) / (MASS + 0.5 * DT * drag);
    let midpoint = 0.5 * (before + after);
    let heat = DT * (0.5 * conductance * pressure * pressure + drag * midpoint * midpoint);
    let supply = DT * pressure * (0.5 * conductance * pressure + 0.5 * WIDTH * GAP * midpoint);
    (after, heat, supply)
}

#[test]
fn supplied_sliding_film_advances_native_motion_heat_and_reservoir_work() {
    let mut model = system(LAUNCH, PRESSURE);
    assert_eq!(model.driven_sliding_film_count(), 1);
    assert_eq!(model.sliding_film_count(), 0);
    assert_eq!(model.squeeze_film_count(), 0);
    assert_eq!(model.driven_sliding_film(0).unwrap().port_count(), 1);
    let initial = model.stored_energy_j();
    let mut expected = LAUNCH;
    let mut heat = 0.0;
    let mut supplied = 0.0;
    for tick in 1..=48 {
        let (after, loss, work) = oracle(expected, PRESSURE);
        let frame = model.step(&[0.0], &CancelGate::new_clock_free()).unwrap();
        expected = after;
        near(speed(&model), expected, 2e-10);
        near(frame.dissipated_energy_j, loss, 2e-11);
        near(frame.supplied_work_j, work, 2e-11);
        assert_eq!(frame.sample, tick);
        near(frame.time_s, tick as f64 * DT, 1e-15);
        heat += frame.dissipated_energy_j;
        supplied += frame.supplied_work_j;
    }
    assert!(speed(&model) > LAUNCH);
    near(model.stored_energy_j() + heat, initial + supplied, 2e-10);
}

#[test]
fn prepared_supplied_tangent_matches_reference_and_preserves_failed_steps() {
    for analytic in [false, true] {
        let mut reference = system(LAUNCH, PRESSURE);
        let model = system(LAUNCH, PRESSURE);
        let mut prepared = if analytic {
            model.prepare_analytic().unwrap()
        } else {
            model.prepare().unwrap()
        };
        assert_eq!(prepared.driven_sliding_film_count(), 1);
        for tick in 0..32 {
            if tick == 8 {
                let state = prepared.state().to_vec();
                let samples = prepared.samples();
                let stop = CancelGate::new_clock_free();
                stop.request();
                assert!(prepared.step(&[0.0], &stop).is_err());
                assert!(
                    prepared
                        .step(&[2.0], &CancelGate::new_clock_free())
                        .is_err()
                );
                assert_eq!(prepared.state(), state);
                assert_eq!(prepared.samples(), samples);
            }
            let expected = reference
                .step(&[0.0], &CancelGate::new_clock_free())
                .unwrap();
            let actual = prepared
                .step(&[0.0], &CancelGate::new_clock_free())
                .unwrap();
            for (&a, &b) in prepared.state().iter().zip(reference.state()) {
                near(a, b, 2e-10);
            }
            near(
                actual.dissipated_energy_j,
                expected.dissipated_energy_j,
                2e-11,
            );
            near(actual.supplied_work_j, expected.supplied_work_j, 2e-11);
        }
    }
}

#[test]
fn reverse_wall_motion_can_return_pressure_work_without_negative_heat() {
    let initial_speed = -0.2;
    let mut model = system(initial_speed, PRESSURE);
    let initial_energy = model.stored_energy_j();
    let (after, heat, supplied) = oracle(initial_speed, PRESSURE);
    let frame = model.step(&[0.0], &CancelGate::new_clock_free()).unwrap();
    assert!(frame.supplied_work_j < 0.0);
    assert!(frame.dissipated_energy_j > 0.0);
    near(speed(&model), after, 2e-10);
    near(frame.dissipated_energy_j, heat, 2e-11);
    near(frame.supplied_work_j, supplied, 2e-11);
    near(
        frame.stored_energy_j + heat,
        initial_energy + supplied,
        2e-11,
    );
}

#[test]
fn zero_supply_matches_passive_film_and_all_images_share_the_island_ceiling() {
    let mut passive = bare(LAUNCH).with_sliding_film(sheared()).unwrap();
    let mut supplied = system(LAUNCH, 0.0);
    for _ in 0..32 {
        let a = passive.step(&[0.0], &CancelGate::new_clock_free()).unwrap();
        let b = supplied
            .step(&[0.0], &CancelGate::new_clock_free())
            .unwrap();
        for (&actual, &expected) in supplied.state().iter().zip(passive.state()) {
            near(actual, expected, 2e-10);
        }
        near(b.dissipated_energy_j, a.dissipated_energy_j, 2e-11);
        assert_eq!(b.supplied_work_j, 0.0);
    }
    let mut full = system(LAUNCH, PRESSURE);
    for _ in 0..3 {
        full = full.with_sliding_film(sheared()).unwrap();
    }
    assert_eq!(full.driven_sliding_film_count(), 1);
    assert_eq!(full.sliding_film_count(), 3);
    assert!(full.with_driven_sliding_film(driven(PRESSURE)).is_err());
}

#[test]
fn supplied_films_cannot_silently_disappear_from_equilibrium_analysis() {
    use crate::stationary_equilibrium::EquilibriumSolveOptions;
    use fs_phs::harmonic::HarmonicBudget;
    use fs_phs::harmonic::equilibrium::stationary::StationaryCoordinates;
    use fs_phs::harmonic::equilibrium::{EquilibriumPoint, LinearizationBudget};
    use fs_solver::NewtonKrylovConfig;

    let model = system(LAUNCH, PRESSURE);
    let original = model.state().to_vec();
    let input = [0.0];
    let scales = [1.0; 2];
    let budget = || LinearizationBudget {
        harmonic: HarmonicBudget {
            max_dimension: 2,
            max_dense_bytes: 1_000_000,
            max_dense_work: 1_000_000,
            relative_residual_tolerance: 1e-10,
        },
        max_hessian_actions: 2,
    };
    let point = EquilibriumPoint {
        state: model.state(),
        input: &input,
        balance_scales: &scales,
        relative_tolerance: 1e-10,
    };
    let refusal = "pressure-fed films require supplied-power equilibrium analysis";
    assert!(matches!(
        model.linearize_equilibrium(point, budget(), 2, || false),
        Err(ImpactError::Invalid(message)) if message == refusal
    ));
    let options = EquilibriumSolveOptions {
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-11,
            relative_tolerance: 0.0,
            linear_restart: 2,
            max_linear_cycles: 8,
            forcing_minimum: 1e-12,
            forcing_maximum: 1e-6,
            ..Default::default()
        },
        relative_balance_tolerance: 1e-10,
        max_iterations: 40,
        max_residual_calls: 256,
        max_jacobian_calls: 2048,
        max_workspace_bytes: 1_000_000,
    };
    let coordinates = StationaryCoordinates {
        input: &input,
        state_scales: &scales,
        balance_scales: &scales,
    };
    assert!(matches!(
        model.solve_equilibrium(model.state(), coordinates, budget(), 2, options, || false),
        Err(ImpactError::Invalid(message)) if message == refusal
    ));
    assert_eq!(model.state(), original);
    assert_eq!(model.samples(), 0);
}
