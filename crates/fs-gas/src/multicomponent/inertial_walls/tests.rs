use super::*;
use crate::multicomponent::tests::{gas, source};
use fs_qty::{Length, Temperature};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

const AREA: f64 = 0.7;
const CLOSED: Boundaries = Boundaries::Ends {
    left: Boundary::Reflecting,
    right: Boundary::Reflecting,
};

fn field(references: [f64; 3]) -> MixtureEuler1d {
    // Manufactured, unequal molar masses and temperature-dependent calorics.
    // The third source is present but its species inventory is exactly zero.
    let gas = gas(vec![
        source("A", R, 3.5, 0.4, references[0]),
        source("B", 2.0 * R, 4.5, 0.2, references[1]),
        source("absent", 0.5 * R, 3.0, 0.1, references[2]),
    ]);
    let count = 24;
    let states = (0..count)
        .map(|i| {
            let first = if i < count / 2 { 0.8 } else { 0.2 };
            gas.conserved(&[first, 1.0 - first, 0.0], 0.0, 1.0).unwrap()
        })
        .collect();
    MixtureEuler1d::new(
        gas,
        Length::new(1.0 / count as f64),
        AREA,
        states,
        CLOSED,
        || false,
    )
    .unwrap()
}

fn controls() -> InertialWallsControl {
    InertialWallsControl {
        midpoint_velocity_bracket_m_s: [[-1.0, 1.0]; 2],
        momentum_tolerance_n_s: 1e-12,
        energy_tolerance_j: 1e-11,
        maximum_iterations: 32,
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}"
    );
}

fn momentum(flow: &MixtureEuler1d, walls: [Piston; 2]) -> f64 {
    flow.totals().momentum
        + walls
            .iter()
            .map(|w| w.mass_kg() * w.velocity_m_s())
            .sum::<f64>()
}
fn energy(flow: &MixtureEuler1d, walls: [Piston; 2]) -> f64 {
    flow.totals().energy
        + walls
            .iter()
            .map(|w| 0.5 * w.mass_kg() * w.velocity_m_s().powi(2))
            .sum::<f64>()
}

#[test]
fn heterogeneous_two_wall_motion_preserves_each_species_and_independent_energy_gauges() {
    let mut flow = field([0.0; 3]);
    let mut shifted = field([1e18, -7e17, 3e17]);
    let mut walls = [
        Piston::new(0.25, 0.0).unwrap(),
        Piston::new(1.1, 0.0).unwrap(),
    ];
    let mut other = walls;
    let loads = [
        InertialWallLoad {
            exterior_pressure_pa: 0.2,
            external_force_n: -0.03,
        },
        InertialWallLoad {
            exterior_pressure_pa: 0.3,
            external_force_n: 0.06,
        },
    ];
    let before = flow.totals();
    let initial_momentum = momentum(&flow, walls);
    let initial_energy = energy(&flow, walls);
    let mut external_impulse = 0.0;
    let mut external_work = 0.0;
    for _ in 0..16 {
        let old_faces = [flow.mesh().positions()[0], flow.mesh().positions()[24]];
        let report = flow
            .step_with_inertial_walls(
                &mut walls,
                Time::new(0.0025),
                StepControl::default(),
                loads,
                controls(),
                || false,
            )
            .unwrap();
        let reference = shifted
            .step_with_inertial_walls(
                &mut other,
                Time::new(0.0025),
                StepControl::default(),
                loads,
                controls(),
                || false,
            )
            .unwrap();
        assert_eq!(report, reference);
        assert_eq!(walls, other);
        assert_eq!(flow.states, shifted.states);
        assert_eq!(flow.mesh(), shifted.mesh());
        assert!(report.admitted(controls()));
        assert_eq!(
            report.gas.after,
            flow.totals(),
            "the report must describe the published candidate"
        );
        assert_eq!(
            report.velocity_after_m_s,
            walls.map(|wall| wall.velocity_m_s())
        );
        let faces = [flow.mesh().positions()[0], flow.mesh().positions()[24]];
        for i in 0..2 {
            let force =
                loads[i].external_force_n + [1.0, -1.0][i] * AREA * loads[i].exterior_pressure_pa;
            external_impulse += force * report.gas.elapsed.value();
            external_work += force * (faces[i] - old_faces[i]);
            assert_eq!(report.displacement_m[i], faces[i] - old_faces[i]);
        }
        for change in [
            report.gas.left_boundary_change,
            report.gas.right_boundary_change,
        ] {
            assert_eq!(change.species_count, 3);
            assert!(change.species_masses.iter().all(|mass| *mass == 0.0));
        }
        for k in 0..3 {
            close(
                flow.totals().species_masses[k],
                before.species_masses[k],
                2e-12,
            );
            close(
                report.gas.conservation_residual.species_masses[k],
                0.0,
                3e-13,
            );
        }
        for cell in flow.states.chunks_exact(5) {
            assert_eq!(
                cell[2], 0.0,
                "pressure and motion must never create an absent species"
            );
        }
        close(
            momentum(&flow, walls) - initial_momentum,
            external_impulse,
            5e-11,
        );
        close(energy(&flow, walls) - initial_energy, external_work, 3e-10);
    }
    assert!(walls[0].velocity_m_s() < -0.02);
    assert!(walls[1].velocity_m_s() > 0.003);
    assert!(
        flow.totals().momentum.abs() > 1e-4,
        "resolved gas momentum cannot be discarded"
    );
    let left = flow
        .gas()
        .primitive(flow.cell(0).unwrap(), || false)
        .unwrap();
    let right = flow
        .gas()
        .primitive(flow.cell(23).unwrap(), || false)
        .unwrap();
    assert!(
        (left.pressure - right.pressure).abs() > 0.05,
        "the local source compositions must not collapse into one cavity pressure"
    );
}

#[test]
fn uniform_two_species_composition_reduces_to_original_frozen_two_wall_dynamics() {
    use crate::dynamics::{CaloricDomain, Euler1d, FrozenGas};
    use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
    use fs_thermochem::{Composition, CompositionBasis};

    let models = vec![
        source("A", R, 3.5, 0.4, 0.0),
        source("B", 2.0 * R, 4.5, 0.2, 0.0),
    ];
    let mixed_gas = gas(models.clone());
    let frozen_gas = FrozenGas::new(
        FrozenIdealGasMixtureModelV1::new(
            models,
            Composition::new(CompositionBasis::MassFraction, vec![0.3, 0.7]).unwrap(),
        )
        .unwrap(),
        CaloricDomain::new(Temperature::new(0.001), Temperature::new(50.0)),
        || false,
    )
    .unwrap();
    let count = 20;
    let mut mixed_states = Vec::new();
    let mut frozen_states = Vec::new();
    for i in 0..count {
        let temperature = 0.9 + 0.2 * (i as f64 + 0.5) / count as f64;
        mixed_states.push(mixed_gas.conserved(&[0.3, 0.7], 0.04, temperature).unwrap());
        frozen_states.push(frozen_gas.conserved(1.0, 0.04, temperature).unwrap());
    }
    let mut mixed = MixtureEuler1d::new(
        mixed_gas,
        Length::new(0.05),
        AREA,
        mixed_states,
        CLOSED,
        || false,
    )
    .unwrap();
    let mut frozen = Euler1d::new(
        frozen_gas,
        Length::new(0.05),
        AREA,
        frozen_states,
        CLOSED,
        || false,
    )
    .unwrap();
    let mut walls = [
        Piston::new(0.3, 0.04).unwrap(),
        Piston::new(0.7, 0.04).unwrap(),
    ];
    let mut other = walls;
    let loads = [
        InertialWallLoad {
            exterior_pressure_pa: 0.2,
            external_force_n: 0.0,
        },
        InertialWallLoad {
            exterior_pressure_pa: 0.25,
            external_force_n: 0.0,
        },
    ];
    // This comparison isolates the source calorics, ALE flux and shared wall
    // mechanics from the distinct species-compatible MUSCL slope limiters.
    let gas_controls = StepControl {
        reconstruction: Reconstruction::FirstOrder,
        ..StepControl::default()
    };
    for _ in 0..12 {
        let a = mixed
            .step_with_inertial_walls(
                &mut walls,
                Time::new(0.0005),
                gas_controls,
                loads,
                controls(),
                || false,
            )
            .unwrap();
        let b = frozen
            .step_with_inertial_walls(
                &mut other,
                Time::new(0.0005),
                gas_controls,
                loads,
                controls(),
                || false,
            )
            .unwrap();
        assert_eq!(a.gas.elapsed, b.gas.elapsed);
        for i in 0..2 {
            close(walls[i].velocity_m_s(), other[i].velocity_m_s(), 3e-10);
            close(a.gas_impulse_n_s[i], b.gas_impulse_n_s[i], 3e-11);
        }
        close(
            a.gas.boundary_change.energy,
            b.gas.boundary_change.energy,
            3e-11,
        );
        for (i, q) in frozen.states().iter().enumerate() {
            let cell = mixed.cell(i).unwrap();
            close(cell.density(), q.density, 3e-10);
            close(cell.partial_densities[0], 0.3 * cell.density(), 3e-12);
            close(cell.momentum, q.momentum, 3e-10);
            close(cell.energy, q.energy, 5e-10);
        }
        for (&x, &y) in mixed
            .mesh()
            .positions()
            .iter()
            .zip(frozen.mesh().positions())
        {
            close(x, y, 3e-11);
        }
    }
}

#[test]
fn cancellation_and_failed_mixture_solve_publish_no_species_wall_geometry_or_time() {
    let duration = Time::new(0.0025);
    let original = [
        Piston::new(0.25, 0.0).unwrap(),
        Piston::new(1.1, 0.0).unwrap(),
    ];
    let mut baseline = field([0.0; 3]);
    let mut walls = original;
    let mut polls = 0;
    baseline
        .step_with_inertial_walls(
            &mut walls,
            duration,
            StepControl::default(),
            [InertialWallLoad::default(); 2],
            controls(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    for cutoff in [1, polls / 2, polls] {
        let mut flow = field([0.0; 3]);
        let states = flow.states.clone();
        let mesh = flow.mesh().clone();
        let mut walls = original;
        let mut seen = 0;
        let result = flow.step_with_inertial_walls(
            &mut walls,
            duration,
            StepControl::default(),
            [InertialWallLoad::default(); 2],
            controls(),
            || {
                seen += 1;
                seen == cutoff
            },
        );
        assert_eq!(result, Err(PistonError::Cancelled));
        assert_eq!(
            seen, cutoff,
            "mixture cancellation must never become a gas retry"
        );
        assert_eq!(flow.states, states);
        assert_eq!(flow.mesh(), &mesh);
        assert_eq!(flow.time(), Time::new(0.0));
        assert_eq!(walls, original);
    }
    let mut flow = field([0.0; 3]);
    let states = flow.states.clone();
    let mesh = flow.mesh().clone();
    let mut walls = original;
    let result = flow.step_with_inertial_walls(
        &mut walls,
        duration,
        StepControl {
            maximum_halvings: 0,
            ..StepControl::default()
        },
        [InertialWallLoad::default(); 2],
        InertialWallsControl {
            midpoint_velocity_bracket_m_s: [[-1e-9, 1e-9]; 2],
            ..controls()
        },
        || false,
    );
    assert!(
        matches!(result, Err(PistonError::StepRetries { .. })),
        "{result:?}"
    );
    assert_eq!(flow.states, states);
    assert_eq!(flow.mesh(), &mesh);
    assert_eq!(flow.time(), Time::new(0.0));
    assert_eq!(walls, original);
    flow.boundaries = Boundaries::Periodic;
    assert!(matches!(
        flow.step_with_inertial_walls(
            &mut walls,
            duration,
            StepControl::default(),
            [InertialWallLoad::default(); 2],
            controls(),
            || false
        ),
        Err(PistonError::Invalid(_))
    ));
    assert_eq!(flow.states, states);
    assert_eq!(walls, original);
}
