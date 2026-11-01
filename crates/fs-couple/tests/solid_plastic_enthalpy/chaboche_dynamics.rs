//! G1/G3 cyclic material heat and consistent inertia on one accepted trajectory.
//! All source coefficients are synthetic numerical fixtures.

use super::*;
use fs_couple::solid_enthalpy::plastic::dynamics::{
    ThermoplasticDynamics, ThermoplasticDynamicsReport,
};
use fs_material::plastic::chaboche::ChabocheState;
use fs_material::plastic::thermal::chaboche::ThermalChabocheSource;

type CyclicMotion<'a> = ThermoplasticDynamics<'a, ThermalChabocheSource<'a, 2>>;

fn motion<'a>(
    fixture: &'a Fixture,
    thermal: &DepositionSimulation<'_>,
    velocity: Vec<[f64; 3]>,
    cx: &Cx<'_>,
) -> CyclicMotion<'a> {
    let selected = [Some(MaterialPropertySelection::SingleClaimOnly)];
    let solid = bind(thermal, &[fixture], &NODES, &TETS, cx)
        .with_voce_hardening(&selected, cx)
        .unwrap()
        .with_power_law_overstress(&selected, cx)
        .unwrap()
        .with_chaboche_branches(
            &[(
                MaterialPropertySelection::SingleClaimOnly,
                chaboche::schema(),
            )],
            cx,
        )
        .unwrap();
    ThermoplasticDynamics::try_new(solid, velocity, cx).unwrap()
}

// Exact integral of the square of an affine velocity over a reference tet:
// integral(N_a N_b) = V(1+delta_ab)/20. This does not call the mass assembler.
fn kinetic(velocity: &[[f64; 3]]) -> f64 {
    let mut square = 0.0;
    for component in 0..3 {
        let sum: f64 = velocity.iter().map(|v| v[component]).sum();
        square += sum * sum
            + velocity
                .iter()
                .map(|v| v[component] * v[component])
                .sum::<f64>();
    }
    RHO * V * square / 40.0
}

fn check_energy(
    report: &ThermoplasticDynamicsReport,
    old_velocity: &[[f64; 3]],
    new_velocity: &[[f64; 3]],
) {
    let delta: Vec<_> = old_velocity
        .iter()
        .zip(new_velocity)
        .map(|(old, new)| core::array::from_fn(|i| new[i] - old[i]))
        .collect();
    close(report.inertia.reference_mass_kg, RHO * V, 1e-18);
    close(
        report.inertia.kinetic_before_j,
        kinetic(old_velocity),
        1e-12,
    );
    close(report.inertia.kinetic_after_j, kinetic(new_velocity), 1e-12);
    close(report.inertia.numerical_loss_j, kinetic(&delta), 1e-12);
    close(
        report.inertia.kinetic_change_j,
        kinetic(new_velocity) - kinetic(old_velocity),
        1e-12,
    );
    close(report.inertia.balance_error_j, 0.0, 1e-12);
    let material = &report.material;
    close(
        material.mechanical_residual_j,
        0.0,
        config().energy_tolerance_j,
    );
    close(
        material.thermal_residual_j,
        0.0,
        config().energy_tolerance_j,
    );
    close(
        material.combined_residual_j,
        0.0,
        config().energy_tolerance_j,
    );
    close(material.external_heat_j, 0.0, 1e-10);
    close(
        material.thermal_audit.stored_energy_change_j,
        material.physical_heat_j + material.reversible_heat_j,
        config().energy_tolerance_j,
    );
}

fn check_physical_heat(
    report: &ThermoplasticDynamicsReport,
    before: &ChabocheState<2>,
    after: &ChabocheState<2>,
    temperature: f64,
    duration: f64,
    viscosity: f64,
    exponent: f64,
) {
    let delta = after.plastic.alpha - before.plastic.alpha;
    assert!(delta > 0.0);
    let overstress = 1e6 * (viscosity * delta / (duration * 1e6)).powf(1.0 / exponent);
    let moduli = [
        4e8 - 8e6 * (temperature - 300.0),
        1e8 + 2e6 * (temperature - 300.0),
    ];
    let recovery = [120.0, 15.0];
    let mut recovery_heat = 0.0;
    for branch in 0..2 {
        let backstress = after.backstress_pa[branch];
        let norm_squared: f64 = backstress[..3].iter().map(|x| x * x).sum::<f64>()
            + 2.0 * backstress[3..].iter().map(|x| x * x).sum::<f64>();
        recovery_heat += 1.5 * recovery[branch] * delta * norm_squared / moduli[branch];
    }
    assert!(recovery_heat > 0.0);
    close(
        report.material.physical_heat_j,
        V * ((YIELD + overstress) * delta + recovery_heat),
        1e-10,
    );
}

#[test]
fn g1_free_cyclic_material_impulse_retains_inertia_and_heats_the_original_cell_once() {
    context(|cx| {
        for exponent in [1.0, 3.0] {
            let viscosity = 1000.0;
            let fixture = chaboche::fixture_with_viscosity(Some(exponent), true, viscosity);
            let thermal = initial(&[&fixture], &[], &[0.0]);
            let mut velocity = vec![[0.0; 3]; 4];
            velocity[1][0] = 30.0;
            let dynamics = motion(&fixture, &thermal, velocity, cx);
            let old = dynamics.solid().solution().clone();
            let duration = 2e-6;
            let first = dynamics
                .advance(
                    &thermal,
                    &(),
                    duration,
                    &ROLLERS,
                    &[],
                    IncrementSettings::default(),
                    thermal_config(),
                    config(),
                    cx,
                )
                .unwrap();
            let repeat = dynamics
                .advance(
                    &thermal,
                    &(),
                    duration,
                    &ROLLERS,
                    &[],
                    IncrementSettings::default(),
                    thermal_config(),
                    config(),
                    cx,
                )
                .unwrap();
            assert_eq!(
                first.dynamics().solid().solution(),
                repeat.dynamics().solid().solution()
            );
            assert_eq!(
                first.dynamics().velocity_mps(),
                repeat.dynamics().velocity_mps()
            );
            assert_eq!(dynamics.solid().solution(), &old);
            check_energy(
                first.report(),
                dynamics.velocity_mps(),
                first.dynamics().velocity_mps(),
            );
            let state = first.dynamics().solid().solution().element_states[0].clone();
            let temperature = first.thermal().temperature_k(CellId(100)).unwrap();
            check_physical_heat(
                first.report(),
                &old.element_states[0],
                &state,
                temperature,
                duration,
                viscosity,
                exponent,
            );
            assert!(temperature > 300.0);
            assert!(first.report().inertia.kinetic_change_j < 0.0);
            assert!(first.report().inertia.numerical_loss_j > 0.0);
            assert!(first.report().material.numerical_loss_j > 0.0);
            assert!(state.backstress_pa.iter().all(|branch| branch[0] > 0.0));
            close(first.report().material.mechanical_work_j, 0.0, 1e-12);
            let (warm, moving, observer, _) = first.into_parts();
            let retained = moving.solid().solution().clone();
            let second = moving
                .advance(
                    &warm,
                    &observer,
                    2.0 * duration,
                    &ROLLERS,
                    &[],
                    IncrementSettings::default(),
                    thermal_config(),
                    config(),
                    cx,
                )
                .unwrap();
            check_energy(
                second.report(),
                moving.velocity_mps(),
                second.dynamics().velocity_mps(),
            );
            assert_eq!(moving.solid().solution(), &retained);
            assert_ne!(
                second.dynamics().solid().solution().displacement_m,
                retained.displacement_m
            );
            assert!(
                second.dynamics().solid().solution().element_states[0]
                    .plastic
                    .alpha
                    >= retained.element_states[0].plastic.alpha
            );
            println!(
                "cyclic impulse n={exponent}: T={temperature:.9} K, alpha={:.9}, heat={:.9} J, K={:.9} J",
                state.plastic.alpha,
                second.report().material.physical_heat_j,
                second.report().inertia.kinetic_after_j
            );
        }
    });
}

#[test]
fn g1_cyclic_drive_hold_and_reversal_close_actual_support_and_kinetic_work() {
    context(|cx| {
        let exponent = 3.0;
        let viscosity = 1000.0;
        let fixture = chaboche::fixture_with_viscosity(Some(exponent), true, viscosity);
        let mut thermal = initial(&[&fixture], &[], &[0.0]);
        let mut dynamics = motion(&fixture, &thermal, vec![[0.0; 3]; 4], cx);
        let mut before_hold_stress = 0.0;
        for (index, (strain, duration)) in [(0.012, 2e-5), (0.012, 4e-5), (-0.012, 3e-5)]
            .into_iter()
            .enumerate()
        {
            let target = dynamics.time_s() + duration;
            let old_solution = dynamics.solid().solution().clone();
            // A prescribed isochoric affine pulse has an exact P1 velocity.
            let prescribed: Vec<_> = NODES
                .iter()
                .enumerate()
                .flat_map(|(node, p)| {
                    [
                        (3 * node, strain * p[0]),
                        (3 * node + 1, -0.5 * strain * p[1]),
                        (3 * node + 2, -0.5 * strain * p[2]),
                    ]
                })
                .collect();
            let step = dynamics
                .advance(
                    &thermal,
                    &(),
                    target,
                    &prescribed,
                    &[],
                    IncrementSettings::default(),
                    thermal_config(),
                    config(),
                    cx,
                )
                .unwrap();
            check_energy(
                step.report(),
                dynamics.velocity_mps(),
                step.dynamics().velocity_mps(),
            );
            let solution = step.dynamics().solid().solution();
            let mut support_work = 0.0;
            for (entry, reaction) in prescribed.iter().zip(&solution.increments[0].reactions_n) {
                let (dof, value) = *entry;
                let displacement = value - old_solution.displacement_m[dof / 3][dof % 3];
                close(
                    step.dynamics().velocity_mps()[dof / 3][dof % 3],
                    displacement / duration,
                    1e-12,
                );
                support_work += reaction * displacement;
            }
            close(
                step.report().material.mechanical_work_j,
                support_work,
                1e-12,
            );
            let temperature = step.thermal().temperature_k(CellId(100)).unwrap();
            check_physical_heat(
                step.report(),
                &old_solution.element_states[0],
                &solution.element_states[0],
                temperature,
                duration,
                viscosity,
                exponent,
            );
            if index == 0 {
                before_hold_stress = solution.element_stress[0][0];
            } else if index == 1 {
                assert!(solution.element_stress[0][0] < before_hold_stress);
                close(kinetic(step.dynamics().velocity_mps()), 0.0, 1e-15);
                close(support_work, 0.0, 1e-15);
            } else {
                assert!(solution.element_stress[0][0] < 0.0);
                assert!(
                    solution.element_states[0]
                        .backstress_pa
                        .iter()
                        .all(|x| x[0] < 0.0)
                );
            }
            assert_eq!(dynamics.solid().solution(), &old_solution);
            println!(
                "cyclic inertial drive {index}: T={temperature:.9} K, stress={:.6} Pa, alpha={:.9}, K={:.9} J",
                solution.element_stress[0][0],
                solution.element_states[0].plastic.alpha,
                step.report().inertia.kinetic_after_j
            );
            let (next_thermal, next_dynamics, _, _) = step.into_parts();
            thermal = next_thermal;
            dynamics = next_dynamics;
        }
    });
}
