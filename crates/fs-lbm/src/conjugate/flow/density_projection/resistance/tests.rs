use super::*;
use crate::conjugate::domain::Voxel;

fn controls(duration_s: f64) -> (DensityProjectionControls, DensityForceControls) {
    // The stiff-step cases exchange up to a million joules. Keep the absolute
    // admission budgets finite and tight at that physical scale.
    let budget = if duration_s > 1.0 { 1e-8 } else { 1e-10 };
    (
        DensityProjectionControls {
            duration_s,
            relative_tolerance: 1e-12,
            maximum_volume_defect_m3_s: 1e-11,
            maximum_energy_defect_j: budget,
        },
        DensityForceControls {
            max_cells: 100,
            maximum_impulse_defect_kg_m_s: budget,
            maximum_energy_defect_j: budget,
        },
    )
}

fn input<'a>(
    density: &'a [f64],
    force: &'a [[f64; 3]],
    acceleration: [f64; 3],
) -> DensityForceInput<'a> {
    DensityForceInput {
        density_kg_m3: density,
        acceleration_m_s2: acceleration,
        force_density_n_m3: force,
    }
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 3e-11 * expected.abs().max(1.0),
        "{actual:.17e} != {expected:.17e}"
    );
}

fn column(axis: usize, axial_widths: Vec<f64>) -> VoxelDomain {
    let mut widths = std::array::from_fn(|_| vec![1.0]);
    widths[axis] = axial_widths;
    widths[(axis + 1) % 3] = vec![2.0];
    widths[(axis + 2) % 3] = vec![3.0];
    VoxelDomain::graded(widths).unwrap()
}

#[test]
#[allow(clippy::too_many_lines)]
fn actual_density_controls_transient_response_and_stiff_darcy_limit_on_every_axis() {
    let gate = CancelGate::new_clock_free();
    for axis in 0..3 {
        let domain = column(axis, vec![3.0]);
        let mut faces = [FlowFace::Wall; 6];
        faces[2 * axis] = FlowFace::Free;
        faces[2 * axis + 1] = FlowFace::Free;
        let initial =
            FlowField::from_face_velocity(
                &domain,
                faces,
                |a, _| {
                    if a == axis { 0.75 } else { 0.0 }
                },
            )
            .unwrap();
        let mut force = [[0.0; 3]];
        force[0][axis] = 2.0;
        for density in [2.0, 128.0] {
            for dt in [0.25, 4096.0] {
                let (projection_controls, force_controls) = controls(dt);
                for difference in [-18.0, 18.0] {
                    let mut baseline: Option<ResistiveDensityProjection> = None;
                    for datum in [0.0, 1e12] {
                        let mut pressure = [None; 6];
                        pressure[2 * axis] = Some(datum + difference);
                        pressure[2 * axis + 1] = Some(datum);
                        let result = initial
                            .force_and_project_density_with_resistance_and_reservoir_pressure(
                                &domain,
                                faces,
                                pressure,
                                input(&[density], &force, [0.0; 3]),
                                &[16.0],
                                &[0.0],
                                projection_controls,
                                force_controls,
                                &gate,
                            )
                            .unwrap();
                        // Integrating the one-dimensional momentum equation
                        // gives rho*du/dt = force + head/L - beta*u. The native
                        // two boundary half-duals together have volume 18 m^3.
                        let drive = 2.0 + difference / 3.0;
                        let speed = (density * 0.75 + dt * drive) / (density + dt * 16.0);
                        let change = speed - 0.75;
                        let mass = density * 18.0;
                        let heat = dt * 16.0 * 18.0 * speed * speed;
                        close(
                            result.projection.field.outward(&domain, 0, 2 * axis),
                            -6.0 * speed,
                        );
                        close(
                            result.projection.field.outward(&domain, 0, 2 * axis + 1),
                            6.0 * speed,
                        );
                        close(result.projection.field.net_outflow(&domain, 0), 0.0);
                        close(
                            result.report.kinetic_energy_before_j,
                            0.5 * mass * 0.75 * 0.75,
                        );
                        close(
                            result.report.kinetic_energy_after_j,
                            0.5 * mass * speed * speed,
                        );
                        close(
                            result.report.kinetic_energy_change_j,
                            mass * (0.75 + 0.5 * change) * change,
                        );
                        close(
                            result.report.backward_euler_loss_j,
                            0.5 * mass * change * change,
                        );
                        close(result.report.resistance_heat_j, heat);
                        close(result.report.cell_resistance_heat_j[0], heat);
                        close(result.report.body_endpoint_work_j, dt * 18.0 * 2.0 * speed);
                        close(result.report.constraint_endpoint_work_j, 0.0);
                        close(result.report.applied_impulse_kg_m_s[axis], dt * 18.0 * 2.0);
                        close(
                            result.report.pressure_impulse_kg_m_s[axis],
                            dt * 6.0 * difference,
                        );
                        close(
                            result.report.resistance_impulse_kg_m_s[axis],
                            -dt * 16.0 * 18.0 * speed,
                        );
                        close(result.report.momentum_change_kg_m_s[axis], mass * change);
                        close(
                            result.projection.report.open_boundary_pressure_work_j,
                            -dt * 6.0 * difference * speed,
                        );
                        assert!(
                            result.report.absolute_energy_defect_j
                                <= force_controls.maximum_energy_defect_j
                        );
                        assert!(
                            result.report.absolute_impulse_defect_kg_m_s
                                <= force_controls.maximum_impulse_defect_kg_m_s
                        );
                        // The steady Darcy velocity is density independent;
                        // its finite-step error retains actual physical inertia.
                        close(
                            speed - drive / 16.0,
                            density / (density + dt * 16.0) * (0.75 - drive / 16.0),
                        );
                        if dt > 1.0 {
                            assert!((speed - drive / 16.0).abs() < 0.002);
                        }
                        if let Some(original) = &baseline {
                            assert_eq!(result.projection.field, original.projection.field);
                            assert_eq!(result.projection.report, original.projection.report);
                            assert_eq!(result.report, original.report);
                        } else {
                            baseline = Some(result);
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn graded_series_resistance_preserves_volume_sources_actual_inertia_and_local_heat() {
    let gate = CancelGate::new_clock_free();
    let (projection_controls, force_controls) = controls(0.4);
    for axis in 0..3 {
        let domain = column(axis, vec![1.0, 2.0]);
        let initial = FlowField::quiescent(&domain);
        let mut faces = [FlowFace::Wall; 6];
        faces[2 * axis] = FlowFace::Free;
        faces[2 * axis + 1] = FlowFace::Free;
        let mut baseline: Option<ResistiveDensityProjection> = None;
        for datum in [0.0, 1e12] {
            let mut pressure = [None; 6];
            pressure[2 * axis] = Some(datum + 8.0);
            pressure[2 * axis + 1] = Some(datum);
            let result = initial
                .force_and_project_density_with_resistance_and_reservoir_pressure(
                    &domain,
                    faces,
                    pressure,
                    input(&[2.0, 5.0], &[[0.0; 3]; 2], [0.0; 3]),
                    &[3.0, 7.0],
                    &[0.8, 0.2],
                    projection_controls,
                    force_controls,
                    &gate,
                )
                .unwrap();
            // Independent series equations: A=6, M=[6,36,30] kg and
            // B=[9,51,42] kg/s. The three pressure drops are W*q, where
            // W=(M+dt*B)/(dt*A^2)=[2/3,47/12,13/4]. Volume sources impose
            // q_middle=q_left+4/5, q_right=q_left+1, and sum(W*q)=8 Pa.
            let q = [97.0 / 470.0, 473.0 / 470.0, 567.0 / 470.0];
            let u = q.map(|q| q / 6.0);
            close(result.projection.field.outward(&domain, 0, 2 * axis), -q[0]);
            close(
                result.projection.field.outward(&domain, 0, 2 * axis + 1),
                q[1],
            );
            close(
                result.projection.field.outward(&domain, 1, 2 * axis + 1),
                q[2],
            );
            close(result.projection.field.net_outflow(&domain, 0), 0.8);
            close(result.projection.field.net_outflow(&domain, 1), 0.2);
            if datum == 0.0 {
                close(result.projection.pressure_correction_pa[0], 5543.0 / 705.0);
                close(result.projection.pressure_correction_pa[1], 7371.0 / 1880.0);
            }
            let kinetic = 0.5 * (6.0 * u[0] * u[0] + 36.0 * u[1] * u[1] + 30.0 * u[2] * u[2]);
            let heat = [
                0.4 * (3.0 * 6.0 / 2.0) * (u[0] * u[0] + u[1] * u[1]),
                0.4 * (7.0 * 12.0 / 2.0) * (u[1] * u[1] + u[2] * u[2]),
            ];
            close(result.report.kinetic_energy_before_j, 0.0);
            close(result.report.kinetic_energy_after_j, kinetic);
            close(result.report.kinetic_energy_change_j, kinetic);
            close(result.report.backward_euler_loss_j, kinetic);
            close(result.report.cell_resistance_heat_j[0], heat[0]);
            close(result.report.cell_resistance_heat_j[1], heat[1]);
            close(result.report.resistance_heat_j, heat[0] + heat[1]);
            close(result.report.pressure_impulse_kg_m_s[axis], 0.4 * 6.0 * 8.0);
            close(
                result.report.resistance_impulse_kg_m_s[axis],
                -0.4 * (9.0 * u[0] + 51.0 * u[1] + 42.0 * u[2]),
            );
            close(
                result.report.momentum_change_kg_m_s[axis],
                6.0 * u[0] + 36.0 * u[1] + 30.0 * u[2],
            );
            close(result.report.body_endpoint_work_j, 0.0);
            close(result.report.constraint_endpoint_work_j, 0.0);
            close(result.projection.report.projection_loss_j, kinetic);
            close(
                result.projection.report.resistance_correction_work_j,
                heat[0] + heat[1],
            );
            close(
                result.projection.report.pressure_source_work_j
                    - result.projection.report.fixed_boundary_pressure_work_j
                    - result.projection.report.open_boundary_pressure_work_j,
                2.0 * kinetic + heat[0] + heat[1],
            );
            close(result.report.energy_defect_j, 0.0);
            if let Some(original) = &baseline {
                // Preserve small hydraulic heads before restoring the datum.
                assert_eq!(result.projection.field, original.projection.field);
                assert_eq!(result.projection.report, original.projection.report);
                assert_eq!(result.report, original.report);
            } else {
                baseline = Some(result);
            }
        }
    }
}

#[test]
#[allow(clippy::too_many_lines)]
fn zero_resistance_preserves_original_projection_and_refusals_preserve_retry() {
    let domain = VoxelDomain::graded([vec![0.5, 1.0, 2.0], vec![0.75, 1.25], vec![1.0]]).unwrap();
    let faces = [
        FlowFace::Fixed,
        FlowFace::Free,
        FlowFace::Wall,
        FlowFace::Wall,
        FlowFace::Wall,
        FlowFace::Wall,
    ];
    let initial = FlowField::from_face_velocity(&domain, faces, |axis, p| match axis {
        0 => 0.2 + 0.1 * p[1],
        1 => 0.03 * p[0],
        _ => 0.0,
    })
    .unwrap();
    let saved = initial.clone();
    let force = [[0.25, -0.125, 0.0]; 6];
    let (projection_controls, force_controls) = controls(0.4);
    let gate = CancelGate::new_clock_free();
    for density in [[2.0; 6], [1.0, 4.0, 8.0, 2.0, 16.0, 3.0]] {
        let forcing = input(&density, &force, [0.125, -0.25, 0.0]);
        let original = initial
            .force_and_project_density(
                &domain,
                faces,
                forcing,
                &[0.0; 6],
                projection_controls,
                force_controls,
                &gate,
            )
            .unwrap();
        let run = |resistance: &[f64], gate: &CancelGate| {
            initial.force_and_project_density_with_resistance(
                &domain,
                faces,
                forcing,
                resistance,
                &[0.0; 6],
                projection_controls,
                force_controls,
                gate,
            )
        };
        let without_drag = run(&[0.0; 6], &gate).unwrap();
        assert_eq!(without_drag.projection, original.projection);
        assert_eq!(without_drag.report.resistance_heat_j, 0.0);
        assert_eq!(without_drag.report.cell_resistance_heat_j, vec![0.0; 6]);
        assert_eq!(without_drag.report.resistance_impulse_kg_m_s, [0.0; 3]);
        assert_eq!(
            without_drag.projection.report.resistance_correction_work_j,
            0.0
        );
        let reference = run(&[1.0, 2.0, 3.0, 4.0, 5.0, 6.0], &gate).unwrap();
        assert!(run(&[1.0; 5], &gate).is_err());
        for invalid in [-1.0, f64::NAN, f64::INFINITY] {
            let mut resistance = [1.0; 6];
            resistance[2] = invalid;
            assert!(run(&resistance, &gate).is_err());
        }
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        assert_eq!(run(&[1.0; 6], &cancelled), Err(ChtError::Cancelled));
        assert_eq!(initial, saved);
        assert_eq!(
            run(&[1.0, 2.0, 3.0, 4.0, 5.0, 6.0], &gate).unwrap(),
            reference
        );
    }
    // A finite coefficient is not an implicit permission to apply fluid drag
    // in an occupied solid, nor a substitute for the exact occupancy mask.
    let occupied = VoxelDomain::graded_from_fn([vec![1.0, 1.0], vec![1.0], vec![1.0]], |p| {
        if p[0] < 1.0 {
            Voxel::Fluid
        } else {
            Voxel::Solid(0)
        }
    })
    .unwrap();
    let still = FlowField::quiescent(&occupied);
    let run = |resistance: &[f64]| {
        still.force_and_project_density_with_resistance(
            &occupied,
            [FlowFace::Wall; 6],
            input(&[2.0, 0.0], &[[0.0; 3]; 2], [0.0; 3]),
            resistance,
            &[0.0; 2],
            projection_controls,
            force_controls,
            &gate,
        )
    };
    let accepted = run(&[1.0, 0.0]).unwrap();
    assert!(run(&[1.0, 1.0]).is_err());
    assert_eq!(run(&[1.0, 0.0]).unwrap(), accepted);
    assert_eq!(still, FlowField::quiescent(&occupied));
}

#[test]
#[allow(clippy::too_many_lines)]
fn hydrostatic_rest_and_held_throughflow_keep_physical_drag_heat_and_support_work() {
    let gate = CancelGate::new_clock_free();
    let (projection_controls, force_controls) = controls(0.4);
    for axis in 0..3 {
        let domain = column(axis, vec![1.0, 2.0]);
        let initial = FlowField::quiescent(&domain);
        let mut faces = [FlowFace::Wall; 6];
        faces[2 * axis + 1] = FlowFace::Free;
        let mut acceleration = [0.0; 3];
        acceleration[axis] = -0.5;
        for datum in [0.0, 1e12] {
            let mut pressure = [None; 6];
            pressure[2 * axis + 1] = Some(datum);
            let result = initial
                .force_and_project_density_with_resistance_and_reservoir_pressure(
                    &domain,
                    faces,
                    pressure,
                    input(&[2.0, 5.0], &[[0.0; 3]; 2], acceleration),
                    &[1e4, 1e-3],
                    &[0.0; 2],
                    projection_controls,
                    force_controls,
                    &gate,
                )
                .unwrap();
            // Integrate actual rho*g over the boundary and interior half-cells.
            // Resistance changes the predictor/metric, but cannot change rest.
            close(result.projection.pressure_correction_pa[0] - datum, 5.5);
            close(result.projection.pressure_correction_pa[1] - datum, 2.5);
            for fluxes in result.projection.field.face_arrays() {
                for &q in fluxes {
                    close(q, 0.0);
                }
            }
            close(result.report.applied_impulse_kg_m_s[axis], -14.4);
            close(result.report.pressure_impulse_kg_m_s[axis], 13.2);
            close(result.report.constraint_impulse_kg_m_s[axis], 1.2);
            close(result.report.momentum_change_kg_m_s[axis], 0.0);
            close(result.report.resistance_heat_j, 0.0);
            close(result.report.backward_euler_loss_j, 0.0);
            close(result.report.body_endpoint_work_j, 0.0);
            close(result.report.constraint_endpoint_work_j, 0.0);
            close(result.report.energy_defect_j, 0.0);
        }
    }
    let domain = VoxelDomain::new(1, 1, 1, 1.0).unwrap();
    let faces = [
        FlowFace::Fixed,
        FlowFace::Free,
        FlowFace::Wall,
        FlowFace::Wall,
        FlowFace::Wall,
        FlowFace::Wall,
    ];
    let initial =
        FlowField::from_face_velocity(&domain, faces, |axis, _| if axis == 0 { 0.75 } else { 0.0 })
            .unwrap();
    let (projection_controls, force_controls) = controls(0.25);
    for datum in [0.0, 1e12] {
        let result = initial
            .force_and_project_density_with_resistance_and_reservoir_pressure(
                &domain,
                faces,
                [None, Some(datum), None, None, None, None],
                input(&[4.0], &[[0.0; 3]], [0.0; 3]),
                &[8.0],
                &[0.0],
                projection_controls,
                force_controls,
                &gate,
            )
            .unwrap();
        close(result.projection.field.flux_x(0, 0, 0), 0.75);
        close(result.projection.field.flux_x(1, 0, 0), 0.75);
        close(result.projection.pressure_correction_pa[0] - datum, 3.0);
        close(result.report.kinetic_energy_before_j, 1.125);
        close(result.report.kinetic_energy_after_j, 1.125);
        close(result.report.kinetic_energy_change_j, 0.0);
        close(result.report.backward_euler_loss_j, 0.0);
        close(result.report.resistance_impulse_kg_m_s[0], -1.5);
        close(result.report.pressure_impulse_kg_m_s[0], 0.75);
        close(result.report.constraint_impulse_kg_m_s[0], 0.75);
        close(result.report.constraint_endpoint_work_j, 0.5625);
        close(result.report.resistance_heat_j, 1.125);
        close(result.report.cell_resistance_heat_j[0], 1.125);
        close(
            result.projection.report.fixed_boundary_pressure_work_j,
            -0.5625,
        );
        close(result.projection.report.open_boundary_pressure_work_j, 0.0);
        close(result.projection.report.pressure_source_work_j, 0.0);
        close(
            result.projection.report.resistance_correction_work_j,
            0.1875,
        );
        close(result.report.energy_defect_j, 0.0);
    }
}
