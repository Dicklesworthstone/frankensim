//! Native directional yielding, laminate coupling and retained Hill shell motion.
//! Synthetic coefficients exercise implementation; they are not a sheet calibration.

use super::*;
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::hill::oriented::OrientedHill48;
use fs_material::plastic::hill::{Hill48, Hill48Coefficients, HillState};
use fs_solid::plastic_shell::dynamics::PlasticShellDensity;
use fs_solid::plastic_shell::energy::{PlasticShellEnergy, ShellEnergyTolerances};

const YIELD: f64 = 1.0;
const HARDENING: f64 = 100.0;

fn material(angle: f64) -> OrientedHill48 {
    let (s, c) = angle.sin_cos();
    OrientedHill48::new(
        Hill48::new(
            IsotropicElastic::new(E, NU, 0.05).unwrap(),
            YIELD,
            HARDENING,
            Hill48Coefficients {
                f: 0.2,
                g: 1.2,
                h: 0.6,
                l: 1.8,
                m: 1.6,
                n: 1.1,
            },
        )
        .unwrap(),
        [[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]],
    )
    .unwrap()
}

fn gates() -> ShellEnergyTolerances {
    ShellEnergyTolerances {
        material: EnergyTolerances {
            energy_j_m3: 1e-9,
            yield_pa: 1e-8,
        },
        section_work_j_m2: 1e-10,
        shell_work_j: 1e-9,
    }
}

fn triangle() -> ShellMesh {
    ShellMesh::new(
        vec![[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        vec![[0, 1, 2]],
    )
    .unwrap()
}

fn membrane_supports(axial: Option<f64>) -> Vec<(usize, f64)> {
    let mut fixed: Vec<_> = (0..3)
        .flat_map(|node| [2, 3, 4].map(move |component| (6 * node + component, 0.0)))
        .collect();
    fixed.extend([(0, 0.0), (1, 0.0), (12, 0.0)]);
    if let Some(strain) = axial {
        fixed.push((6, strain));
    }
    // Width, in-plane shear and relative drilling remain free. On release,
    // the axial length also follows equilibrium rather than a prescribed set.
    fixed
}

fn bending_endpoint(mesh: &ShellMesh, curvature: f64) -> Vec<(usize, f64)> {
    mesh.nodes
        .iter()
        .flat_map(|&[x, _, _]| [0.0, 0.0, 0.5 * curvature * x * x, 0.0, -curvature * x, 0.0])
        .enumerate()
        .collect()
}

fn check_yield_heat(
    section: &LayeredPlaneStressSection<OrientedHill48>,
    old: Option<&SmallStrainShellSolution<HillState>>,
    new: &SmallStrainShellSolution<HillState>,
    energy: &PlasticShellEnergy,
) {
    let mut expected = 0.0;
    for (q, state) in new.section_states.iter().enumerate() {
        for (i, point) in section.points().iter().enumerate() {
            let before = old.map_or(0.0, |old| old.section_states[q].points[i].history.alpha);
            let after = state.points[i].history.alpha;
            assert!(after >= before);
            expected += (0.5 / 3.0) * point.weight_m * YIELD * (after - before);
        }
    }
    near(energy.physical_dissipation_j, expected, 2e-12);
    near(energy.balance_residual_j, 0.0, 1e-9);
}

#[test]
fn g1_directional_membranes_solve_width_shear_and_permanent_extension_on_release() {
    with_cx(|_, cx| {
        let mesh = triangle();
        // q=a|sigma_x|. The associated spatial flow is sign(sigma_x)*n/a.
        // At 45 degrees a²=(F+G+2N)/4 and n_xy=(G-F)/4.
        let directions = [
            (0.0, 1.8_f64, [1.8, -0.6, -1.2, 0.0, 0.0, 0.0]),
            (
                std::f64::consts::FRAC_PI_2,
                0.8,
                [0.8, -0.6, -0.2, 0.0, 0.0, 0.0],
            ),
            (
                std::f64::consts::FRAC_PI_4,
                0.9,
                [0.9, -0.2, -0.7, 0.25, 0.0, 0.0],
            ),
        ];
        let mut peak_stress = [0.0; 3];
        for (direction, (angle, a2, normal)) in directions.into_iter().enumerate() {
            let law = material(angle);
            let section = section(law.clone());
            let mut accepted: Option<SmallStrainShellSolution<HillState>> = None;
            let mut plastic = [0.0; 6];
            let mut alpha = 0.0;
            for axial in [0.0005, 0.0015, 0.004] {
                // Uniaxial Hooke/associated-flow equations give this scalar
                // return independently of the native six-component return.
                let a = a2.sqrt();
                let trial = E * (axial - plastic[0]);
                let gamma = (a * trial - YIELD - HARDENING * alpha).max(0.0) / (E * a2 + HARDENING);
                alpha += gamma;
                for i in 0..6 {
                    plastic[i] += gamma * normal[i] / a;
                }
                let stress = E * (axial - plastic[0]);
                let fixed = membrane_supports(Some(axial));
                let step = problem(&mesh, &section, &fixed, &[], &[1.0], 0.37)
                    .solve_with_energy(accepted.as_ref(), gates(), cx)
                    .unwrap();
                near(step.solution.coordinates[1][0], axial, 1e-17);
                near(
                    step.solution.coordinates[2][1],
                    plastic[1] - NU * stress / E,
                    2e-10,
                );
                near(step.solution.coordinates[1][1], 2.0 * plastic[3], 2e-10);
                for (state, resultants) in step
                    .solution
                    .section_states
                    .iter()
                    .zip(&step.solution.section_resultants)
                {
                    near(resultants[0], H * stress, 2e-9);
                    for value in &resultants[1..] {
                        near(*value, 0.0, 2e-9);
                    }
                    for point in &state.points {
                        near(point.history.alpha, alpha, 2e-10);
                        near(point.strain[2], plastic[2] - NU * stress / E, 2e-10);
                        for (actual, expected) in law
                            .plastic_strain_spatial(&point.history)
                            .iter()
                            .zip(plastic)
                        {
                            near(*actual, expected, 2e-10);
                        }
                    }
                }
                check_yield_heat(&section, accepted.as_ref(), &step.solution, &step.energy);
                near(
                    step.energy.final_state.elastic_j,
                    0.5 * H * stress * stress / (2.0 * E),
                    2e-11,
                );
                near(
                    step.energy.final_state.hardening_j,
                    0.5 * H * 0.5 * HARDENING * alpha * alpha,
                    2e-11,
                );
                peak_stress[direction] = stress;
                accepted = Some(step.solution);
            }
            let loaded = accepted.unwrap();
            let saved = loaded.clone();
            let fixed = membrane_supports(None);
            let released = problem(&mesh, &section, &fixed, &[], &[1.0], 0.37)
                .solve_with_energy(Some(&loaded), gates(), cx)
                .unwrap();
            near(released.solution.coordinates[1][0], plastic[0], 2e-10);
            near(released.solution.coordinates[2][1], plastic[1], 2e-10);
            near(released.solution.coordinates[1][1], 2.0 * plastic[3], 2e-10);
            assert!(released.solution.coordinates[1][0] > 0.001);
            near(released.energy.final_state.elastic_j, 0.0, 2e-12);
            near(released.energy.physical_dissipation_j, 0.0, 2e-12);
            check_yield_heat(
                &section,
                Some(&loaded),
                &released.solution,
                &released.energy,
            );
            for state in &released.solution.section_states {
                for point in &state.points {
                    near(point.history.alpha, alpha, 2e-10);
                }
            }
            assert_eq!(loaded, saved);
        }
        assert!(peak_stress[0] < peak_stress[2] && peak_stress[2] < peak_stress[1]);
    });
}

fn laminate(reverse: bool) -> LayeredPlaneStressSection<OrientedHill48> {
    let angles = if reverse {
        [std::f64::consts::FRAC_PI_2, 0.0]
    } else {
        [0.0, std::f64::consts::FRAC_PI_2]
    };
    LayeredPlaneStressSection::new(
        angles
            .into_iter()
            .enumerate()
            .map(|(i, angle)| PlaneStressLayer {
                bottom_z_m: (i as f64 - 1.0) * H * 0.5,
                top_z_m: i as f64 * H * 0.5,
                subdivisions: 2,
                law: material(angle),
            })
            .collect(),
        ThicknessQuadrature::Gauss3,
        ReductionConfig {
            max_iterations: 48,
            stress_tolerance_pa: 1e-12,
            strain_limit: 0.05,
        },
        SectionBudget {
            maximum_layers: 2,
            maximum_quadrature_points: 12,
        },
    )
    .unwrap()
}

#[test]
fn g3_reversing_the_rolled_layer_stack_reverses_plastic_membrane_bending_coupling() {
    with_cx(|_, cx| {
        let mesh = triangle();
        let sections = [laminate(false), laminate(true)];
        let mut accepted: [Option<SmallStrainShellSolution<HillState>>; 2] = [None, None];
        for curvature in [0.01, 0.04] {
            let fixed = bending_endpoint(&mesh, curvature);
            let mut heat = [0.0; 2];
            for i in 0..2 {
                let step = problem(&mesh, &sections[i], &fixed, &[], &[1.0], 0.37)
                    .solve_with_energy(accepted[i].as_ref(), gates(), cx)
                    .unwrap();
                check_yield_heat(
                    &sections[i],
                    accepted[i].as_ref(),
                    &step.solution,
                    &step.energy,
                );
                heat[i] = step.energy.physical_dissipation_j;
                accepted[i] = Some(step.solution);
            }
            near(heat[0], heat[1], 2e-12);
            let a = accepted[0].as_ref().unwrap();
            let b = accepted[1].as_ref().unwrap();
            for q in 0..3 {
                for i in 0..6 {
                    let parity = if i < 3 { -1.0 } else { 1.0 };
                    near(
                        a.section_resultants[q][i],
                        parity * b.section_resultants[q][i],
                        2e-11,
                    );
                }
                for (left, right) in a.section_states[q]
                    .points
                    .iter()
                    .zip(b.section_states[q].points.iter().rev())
                {
                    near(left.history.alpha, right.history.alpha, 2e-11);
                    for i in 0..6 {
                        near(
                            left.history.plastic_strain[i],
                            -right.history.plastic_strain[i],
                            2e-11,
                        );
                    }
                }
            }
            if curvature == 0.01 {
                near(heat[0], 0.0, 0.0);
                near(a.section_resultants[0][0], 0.0, 2e-12);
            } else {
                assert!(heat[0] > 0.0);
                assert!(a.section_resultants[0][0].abs() > 1e-5);
            }
        }
    });
}

#[test]
fn g1_oblique_hill_bending_keeps_reference_inertia_and_separate_plastic_heat() {
    with_cx(|_, cx| {
        let mesh = triangle();
        let section = section(material(std::f64::consts::FRAC_PI_4));
        let owner = problem(&mesh, &section, &[], &[], &[1.0], 0.37);
        let density = 2.0;
        let dt = 0.25;
        let reference = owner
            .prepare_reference_inertia(PlasticShellDensity::Uniform(density), cx)
            .unwrap();
        let mut accepted: Option<SmallStrainShellSolution<HillState>> = None;
        let mut velocity = vec![[0.0; 6]; 3];
        let mut previous_curvature = 0.0;
        let mut previous_rate = 0.0;
        let kinetic = |rate: f64| {
            // Only node (1,0) moves: integral N1² dA = A/6. The symmetric
            // reference section has no translation/director cross moment.
            0.5 * (density * H * 0.5 / 6.0) * ((0.5 * rate).powi(2) + H * H / 12.0 * rate.powi(2))
        };
        let mut heat = 0.0;
        for curvature in [0.04, 0.02, 0.0] {
            let fixed = bending_endpoint(&mesh, curvature);
            let step = problem(&mesh, &section, &fixed, &[], &[1.0], 0.37)
                .advance_backward_euler_with_energy(
                    accepted.as_ref(),
                    &velocity,
                    &reference,
                    dt,
                    gates(),
                    cx,
                )
                .unwrap();
            let rate = (curvature - previous_curvature) / dt;
            near(
                step.inertia.physical.kinetic_before_j,
                kinetic(previous_rate),
                2e-15,
            );
            near(step.inertia.physical.kinetic_after_j, kinetic(rate), 2e-15);
            near(
                step.inertia.physical.numerical_loss_j,
                kinetic(rate - previous_rate),
                2e-15,
            );
            near(step.inertia.drilling_gauge.kinetic_after_j, 0.0, 1e-25);
            near(step.balance_residual_j(), 0.0, 1e-9);
            check_yield_heat(&section, accepted.as_ref(), &step.solution, &step.energy);
            heat += step.energy.physical_dissipation_j;
            previous_curvature = curvature;
            previous_rate = rate;
            velocity = step.velocity;
            accepted = Some(step.solution);
        }
        assert!(heat > 0.0);
        assert!(
            accepted
                .unwrap()
                .section_states
                .iter()
                .any(|state| state.points.iter().any(|point| point.history.alpha > 0.0))
        );
    });
}
