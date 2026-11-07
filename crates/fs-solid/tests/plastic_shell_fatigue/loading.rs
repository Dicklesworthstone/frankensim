//! G1/G3: physical, independently phased boundary fields reach native fatigue.
//! Synthetic life coefficients are implementation references, not measured fits.

use super::shear::{fixture, fs_gauge, shear_config};
use super::*;
use fs_solid::plastic_shell::fatigue::ShellFatigueLoading;

const A: f64 = 1e-4;
const G: f64 = E / (2.0 * (1.0 + NU));

/// Trace-free in-plane strain: [a, -a, 0, b, 0, 0], with zero rigid spin.
fn affine_values(mesh: &ShellMesh, a: f64, b: f64) -> Vec<f64> {
    let frame = mesh.facet(0).unwrap().frame;
    let origin = mesh.nodes[0];
    let mut values = Vec::new();
    for point in &mesh.nodes {
        let x: f64 = (0..3).map(|c| frame[0][c] * (point[c] - origin[c])).sum();
        let y: f64 = (0..3).map(|c| frame[1][c] * (point[c] - origin[c])).sum();
        for c in 0..6 {
            values.push(if c < 3 {
                frame[0][c] * (a * x + b * y) + frame[1][c] * (b * x - a * y)
            } else {
                0.0
            });
        }
    }
    values
}

fn unchanged<L: SourceBoundFatigueLaw>(
    a: &ShellFatiguePath<'_, '_, L>,
    b: &ShellFatiguePath<'_, '_, L>,
) {
    assert_eq!(a.solution(), b.solution());
    assert_eq!(a.velocity(), b.velocity());
    assert_eq!(a.elapsed_time_s(), b.elapsed_time_s());
    assert_eq!(a.next_step(), b.next_step());
    assert_eq!(a.sample_count(), b.sample_count());
    assert_eq!(a.samples(), b.samples());
    assert_eq!(a.plane_summaries(0), b.plane_summaries(0));
    assert_eq!(a.results(), b.results());
    assert_eq!(a.fatemi_socie_samples(), b.fatemi_socie_samples());
    assert_eq!(a.fatemi_socie_results(), b.fatemi_socie_results());
}

#[test]
fn phased_native_shear_triangle_uses_enclosing_circle_and_whole_cycle_normal_stress() {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let f = fixture(0);
        let section = section(f.law.clone());
        let mesh = mesh();
        let fixed: Vec<_> = (0..18).map(|d| (d, 0.0)).collect();
        let factors = [1.0; 4];
        let p = problem(&mesh, &section, &fixed, &[], &factors);
        let mut shear = fs_gauge(&f.law, &f.calibration);
        shear.normals_in_facet_frame = vec![[1.0, 1.0, 1.0]];
        let reference = p
            .prepare_reference_inertia(PlasticShellDensity::Uniform(1.0), cx)
            .unwrap();
        for dynamic in [false, true] {
            let mut run = ShellFatiguePath::new_with_fatemi_socie(
                &p,
                vec![gauge(&f.law, 1, StrainLifeMeanStress::Uncorrected)],
                vec![shear.clone()],
                config(),
                shear_config(),
                cx,
            )
            .unwrap();
            // a=A*cos(theta), b=sqrt(3)*A*sin(theta), theta=0,120,240,360.
            // On n=(1,1,1)/sqrt(3), the engineering-shear path is an equilateral
            // triangle of circumradius sqrt(8/3)*A. Its longest half-chord is only
            // sqrt(3)/2 of that radius; using a scalar range loses real path geometry.
            let endpoints = [
                (A, 0.0),
                (-0.5 * A, 1.5 * A),
                (-0.5 * A, -1.5 * A),
                (A, 0.0),
            ];
            for (index, (a, b)) in endpoints.into_iter().enumerate() {
                let values = affine_values(&mesh, a, b);
                let loading = ShellFatigueLoading {
                    prescribed_values: &values,
                    nodal_loads: &[],
                };
                if index == 2 {
                    let prefix = run.clone();
                    // A real candidate endpoint is solved before full-tensor
                    // nonclosure refuses. It must consume neither loads nor state.
                    if dynamic {
                        assert!(
                            run.advance_next_backward_euler_cycle_with_loading(
                                &reference, 0.03, loading, cx
                            )
                            .is_err()
                        );
                    } else {
                        assert!(run.advance_next_cycle_with_loading(loading, cx).is_err());
                    }
                    unchanged(&run, &prefix);
                }
                if dynamic && index == 3 {
                    assert!(
                        run.advance_next_backward_euler_cycle_with_loading(
                            &reference, 0.03, loading, cx
                        )
                        .unwrap()
                        .is_some()
                    );
                } else if dynamic {
                    assert!(
                        run.advance_next_backward_euler_with_loading(&reference, 0.03, loading, cx)
                            .unwrap()
                            .is_some()
                    );
                } else if index == 3 {
                    assert!(run.advance_next_cycle_with_loading(loading, cx).unwrap());
                } else {
                    assert!(run.advance_next_with_loading(loading, cx).unwrap());
                }
                let sample = run.fatemi_socie_samples()[0];
                let expected = [a, -a, 0.0, b, 0.0, 0.0];
                for component in 0..6 {
                    close(
                        sample.total_mechanical_strain[component],
                        expected[component],
                        1e-13,
                    );
                    close(
                        sample.stress_pa[component],
                        2.0 * G * expected[component],
                        1e-10,
                    );
                }
                assert_eq!(
                    run.solution().section_states[0].points[1].history.alpha,
                    0.0
                );
                if index == 0 {
                    let prefix = run.clone();
                    run.restart_fatemi_socie_cycle(cx).unwrap();
                    assert_eq!(run.solution(), prefix.solution());
                    assert_eq!(run.plane_summaries(0), prefix.plane_summaries(0));
                    assert_eq!(run.sample_count(), prefix.sample_count());
                }
            }
            let result = run.fatemi_socie_results()[0].unwrap();
            let radius = (8.0_f64 / 3.0).sqrt() * A;
            let maximum = 2.0 * G * A;
            let correction = 1.0 + 0.2 * maximum;
            close(result.engineering_shear_strain_amplitude, radius, 1e-13);
            close(result.maximum_normal_stress_pa, maximum, 1e-10);
            close(result.minimum_normal_stress_pa, -maximum, 1e-10);
            close(result.fatemi_socie_parameter, radius * correction, 1e-13);
            close(
                result.cycles_to_failure.unwrap(),
                0.5 * (0.2 / (radius * correction)).powi(2),
                1e-4,
            );
            assert_eq!(result.samples, 4);
            assert_eq!(run.sample_count(), 5); // Virgin state is still in normal Miner history.
            assert!(run.finish(cx).unwrap()[0].summary.miner_damage > 0.0);
            close(
                run.elapsed_time_s(),
                if dynamic { 0.12 } else { 0.0 },
                1e-15,
            );
        }
    });
}

#[test]
fn independent_forces_and_support_motion_use_original_dynamic_scaling_and_inertia() {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let f = fixture(0);
        let section = section(f.law.clone());
        let mesh = mesh();
        let fixed = membrane_prescribed();
        let factors = [0.5, 1.0, -0.75];
        let p = problem(&mesh, &section, &fixed, &[], &factors);
        let reference = p
            .prepare_reference_inertia(PlasticShellDensity::Uniform(1.0), cx)
            .unwrap();
        let mut run = ShellFatiguePath::new(
            &p,
            vec![gauge(&f.law, 1, StrainLifeMeanStress::Uncorrected)],
            config(),
            cx,
        )
        .unwrap();
        let mut native = problem(&mesh, &section, &fixed, &[], &[0.0])
            .solve(cx)
            .unwrap();
        let mut velocity = vec![[0.0; 6]; 3];
        let dt = 0.03;
        for (index, factor) in factors.into_iter().enumerate() {
            let axial = [A, -0.8 * A, 1.2 * A][index];
            let force = [0.003, -0.005, 0.001][index];
            let values: Vec<_> = fixed
                .iter()
                .map(|&(d, _)| if d == 6 { axial } else { 0.0 })
                .collect();
            // Two additive forces on the free y DOF retain native load semantics.
            let forces = [(13, force), (13, 0.25 * force)];
            let prescribed: Vec<_> = fixed
                .iter()
                .zip(&values)
                .map(|(&(d, _), &v)| (d, factor * v))
                .collect();
            let scaled_forces: Vec<_> = forces.iter().map(|&(d, v)| (d, factor * v)).collect();
            let endpoint = problem(&mesh, &section, &prescribed, &scaled_forces, &[1.0]);
            let mut step = endpoint
                .advance_backward_euler(Some(&native), &velocity, &reference, dt, cx)
                .unwrap();
            let energy = run
                .advance_next_backward_euler_with_loading(
                    &reference,
                    dt,
                    ShellFatigueLoading {
                        prescribed_values: &values,
                        nodal_loads: &forces,
                    },
                    cx,
                )
                .unwrap()
                .unwrap();
            step.solution.increments[0].load_factor = factor;
            assert_eq!(run.solution(), &step.solution);
            assert_eq!(run.velocity(), step.velocity);
            assert_eq!(energy, step.inertia);
            let point = &step.solution.section_states[0].points[1];
            assert_eq!(
                run.samples()[0],
                f.law
                    .accepted_fatigue_sample(&point.strain, &point.history)
                    .unwrap()
            );
            native = step.solution;
            velocity = step.velocity;
        }
        close(run.elapsed_time_s(), 3.0 * dt, 1e-15);
        assert!(run.finish(cx).unwrap()[0].summary.miner_damage > 0.0);
    });
}

#[test]
fn explicit_loading_refuses_support_count_and_invalid_force_before_state_change() {
    with_cx(&CancelGate::new_clock_free(), |cx| {
        let f = fixture(0);
        let section = section(f.law.clone());
        let mesh = mesh();
        let fixed: Vec<_> = (0..18).map(|d| (d, 0.0)).collect();
        let p = problem(&mesh, &section, &fixed, &[], &[1.0, 1.0]);
        let mut run = ShellFatiguePath::new(
            &p,
            vec![gauge(&f.law, 1, StrainLifeMeanStress::Uncorrected)],
            config(),
            cx,
        )
        .unwrap();
        let values = affine_values(&mesh, A, 0.0);
        run.advance_next_with_loading(
            ShellFatigueLoading {
                prescribed_values: &values,
                nodal_loads: &[],
            },
            cx,
        )
        .unwrap();
        let prefix = run.clone();
        assert!(
            run.advance_next_with_loading(
                ShellFatigueLoading {
                    prescribed_values: &values[..17],
                    nodal_loads: &[]
                },
                cx
            )
            .is_err()
        );
        unchanged(&run, &prefix);
        for invalid in [[(18, 1.0)], [(6, f64::NAN)]] {
            assert!(
                run.advance_next_with_loading(
                    ShellFatigueLoading {
                        prescribed_values: &values,
                        nodal_loads: &invalid
                    },
                    cx
                )
                .is_err()
            );
            unchanged(&run, &prefix);
        }
        let mut invalid = values.clone();
        invalid[6] = f64::INFINITY;
        assert!(
            run.advance_next_with_loading(
                ShellFatigueLoading {
                    prescribed_values: &invalid,
                    nodal_loads: &[]
                },
                cx
            )
            .is_err()
        );
        unchanged(&run, &prefix);
        let opposite = affine_values(&mesh, -A, 0.0);
        assert!(
            run.advance_next_with_loading(
                ShellFatigueLoading {
                    prescribed_values: &opposite,
                    nodal_loads: &[]
                },
                cx
            )
            .unwrap()
        );
        assert!(run.is_complete());
    });
}
