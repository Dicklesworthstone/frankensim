//! G1/G3: thermal fracture retains native compression, contact forces and source work.

use super::approach::{balanced_force_intervals, minimum_gap};
use super::*;
use fs_contact::sustained::solve::{SolveControls, SolveError};
use fs_solid::brittle::solid::dynamics::fracture::contact::sustained::{
    BrittleFaceSustainedControls, BrittleFaceSustainedError,
};
use fs_solver::NewtonKrylovConfig;

#[path = "sustained/modes.rs"]
mod modes;

pub(super) fn sustained_controls() -> BrittleFaceSustainedControls {
    BrittleFaceSustainedControls {
        solve: SolveControls {
            newton: NewtonKrylovConfig {
                absolute_tolerance: 1e-11,
                relative_tolerance: 0.0,
                ..NewtonKrylovConfig::default()
            },
            maximum_response_evaluations: 512,
            difference_step_n: 1e-3,
            force_residual_tolerance_n: 1e-9,
            maximum_contacts: 3,
            ..SolveControls::default()
        },
        work: midpoint_controls(),
        normal_scale_n_per_m: 4096.0,
        activation_tolerance_m: 1e-10,
        tangential_tolerance_m: 1e-10,
        velocity_tolerance_mps: 1e-8,
        force_tolerance_n: 1e-8,
        complementarity_tolerance_j: 1e-10,
        contact_work_tolerance_j: 1e-10,
        force_moment_tolerance_n_m: 1e-10,
    }
}

/// Bind matter at 400 K, cool to failure, optionally unload the normal strain,
/// and stop the accepted motion before cutting. All three steps are native;
/// no post-cooling initializer can replace rho=6 with the current query's rho=3.
pub(super) fn with_closed_thermal_motion(
    unload_normal: bool,
    maximum_steps: usize,
    run: impl FnOnce(&ThermalBrittleTetMotion<'_>, &[ResolvedThermalBrittle], &Cx<'_>),
) {
    with_cx(|cx| {
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, budget(), 8, cx).unwrap();
        let card = source_card();
        let reference = reference_sources(&card, &geometry);
        let cold = at_temperature(&reference, &card, 300.0);
        let warm = at_temperature(&reference, &card, 350.0);
        let solid = ThermalBrittleTetSolid::try_new(&geometry, &reference, budget(), cx).unwrap();
        let initial =
            ThermalBrittleTetMotion::try_from_solid(&solid, &[[0.0; 3]; 5], 3.0, maximum_steps, cx)
                .unwrap();
        let cooled = initial
            .advance(
                &cold,
                INTERVAL,
                &prescribed(&[[0.0; 3]; 5]),
                &[],
                settings(),
                backward_controls(),
                cx,
            )
            .unwrap()
            .motion;
        let normal_strain = if unload_normal {
            cold[0].free_linear_strain()
        } else {
            0.0
        };
        let displacement = NODES.map(|x| [0.0, 0.0, normal_strain * x[2]]);
        let unloaded = cooled
            .advance(
                &cold,
                INTERVAL,
                &prescribed(&displacement),
                &[],
                settings(),
                backward_controls(),
                cx,
            )
            .unwrap()
            .motion;
        let stopped = unloaded
            .advance(
                &cold,
                INTERVAL,
                &prescribed(&displacement),
                &[],
                settings(),
                backward_controls(),
                cx,
            )
            .unwrap()
            .motion;
        let cut = topology(&geometry, cx);
        let target =
            TetLoadGeometry3::new(cut.nodes_m(), cut.tetrahedra(), budget(), 8, cx).unwrap();
        let split = stopped
            .split_fracture(&cut, &target, fracture_controls(), 1e-12, cx)
            .unwrap()
            .motion;
        assert_eq!(split.accepted_steps(), 3);
        assert_eq!(split.time_s(), 3.0 + 3.0 * INTERVAL);
        assert_eq!(
            split.velocity_mps(),
            &vec![[0.0; 3]; target.nodes_m().len()]
        );
        assert_eq!(split.reference_density_kg_m3(), &[DENSITY; 2]);
        assert!(core::ptr::eq(
            split.reference_density_kg_m3().as_ptr(),
            initial.reference_density_kg_m3().as_ptr()
        ));
        assert_eq!(
            cold[0].band().reference_density_kg_m3().unwrap(),
            0.5 * DENSITY
        );
        close(split.reference_mass_kg(), 2.0, 1e-15);
        close(split.kinetic_energy_j(), 0.0, 1e-15);
        close(split.solid().stored_energy_j(), 0.04, 1e-13);
        close(split.solid().fracture_dissipation_j(), GC, 1e-14);
        close(minimum_gap(&split), 0.0, 1e-15);
        for state in &split.solid().solution().element_states {
            assert!(state.failed());
            close(state.strain()[0], 0.01, 1e-15);
            close(state.strain()[1], 0.01, 1e-15);
            close(
                state.strain()[2],
                if unload_normal { 0.0 } else { 0.01 },
                1e-15,
            );
        }
        run(&split, &warm, cx);
    });
}

fn body_forces(motion: &ThermalBrittleTetMotion<'_>, inward: f64) -> Vec<(usize, f64)> {
    balanced_force_intervals(motion, [0.0; 2], inward)
        .into_iter()
        .map(|(dof, pair)| {
            assert_eq!(pair[0], pair[1]);
            (dof, pair[0])
        })
        .collect()
}

fn retained_matter_and_failure(
    before: &ThermalBrittleTetMotion<'_>,
    after: &ThermalBrittleTetMotion<'_>,
) {
    assert!(core::ptr::eq(
        before.reference_density_kg_m3().as_ptr(),
        after.reference_density_kg_m3().as_ptr()
    ));
    assert_eq!(after.reference_density_kg_m3(), &[DENSITY; 2]);
    close(after.reference_mass_kg(), 2.0, 1e-15);
    close(after.solid().fracture_dissipation_j(), GC, 1e-13);
    assert_eq!(
        before.solid().released_faces(),
        after.solid().released_faces()
    );
    for (old, next) in before
        .solid()
        .solution()
        .element_states
        .iter()
        .zip(&after.solid().solution().element_states)
    {
        assert!(next.failed());
        assert_eq!(next.damage(), old.damage());
        assert_eq!(next.dissipated_j_m3(), old.dissipated_j_m3());
    }
}

#[test]
fn closed_thermal_fragments_recover_compression_with_original_consistent_mass() {
    for compression in [false, true] {
        with_closed_thermal_motion(compression, 6, |motion, _warm, cx| {
            let h = INTERVAL;
            let acceleration = 1.0;
            let result = motion
                .advance_released_face_sustained(
                    h,
                    &[],
                    &body_forces(motion, acceleration),
                    settings(),
                    sustained_controls(),
                    cx,
                )
                .unwrap();
            // One kg per tetrahedron: M=(I+ones)/20. Closed symmetric
            // faces stay at zero normal displacement. When normal mechanical
            // strain starts at zero, native mean compression adds E*d/2 and
            // (1/(5*h^2)+E/12)*d=a/4 at each free apex.
            let young = if compression { 1200.0 } else { 0.0 };
            let denominator = 1.0 + 5.0 * young * h * h / 12.0;
            let d = 1.25 * acceleration * h * h / denominator;
            let geometry = motion.solid().geometry();
            assert_eq!(result.rows.len(), 3);
            assert!(result.solve.nonlinear.converged);
            for (row, force) in result.rows.iter().zip(&result.solve.forces_n) {
                let expected = if geometry.nodes_m()[row.nodes[0]] == [0.0; 3] {
                    acceleration / 2.0 - 3.0 * acceleration / (8.0 * denominator)
                } else {
                    acceleration / 4.0 - acceleration / (8.0 * denominator)
                };
                close(force[0], expected, 1e-8);
                assert_eq!(force[1], 0.0);
                assert_eq!(force[2], 0.0);
            }
            for (element, tet) in geometry.tetrahedra().iter().enumerate() {
                let sign = if element == 0 { 1.0 } else { -1.0 };
                for &node in tet {
                    let increment = if geometry.nodes_m()[node][2] == 0.0 {
                        0.0
                    } else {
                        sign * d
                    };
                    close(
                        result.motion.solid().solution().displacement_m[node][2]
                            - motion.solid().solution().displacement_m[node][2],
                        increment,
                        1e-10,
                    );
                    close(
                        result.motion.velocity_mps()[node][2],
                        2.0 * increment / h,
                        1e-8,
                    );
                }
                close(
                    result.motion.solid().solution().element_stress[element][2],
                    -young * d,
                    1e-8,
                );
                assert!(core::ptr::eq(
                    result.motion.solid().solution().element_states[element]
                        .cohesive_point()
                        .law(),
                    motion.solid().solution().element_states[element]
                        .cohesive_point()
                        .law()
                ));
            }
            let storage = young * d * d / 6.0;
            let kinetic = 0.625 * (acceleration * h / denominator).powi(2);
            let work = 0.5 * acceleration * d;
            close(result.energy.stored_after_j, 0.04 + storage, 1e-10);
            close(result.energy.inertia.kinetic_change_j, kinetic, 1e-10);
            close(inventory(&result.motion).2, kinetic, 1e-10);
            close(result.energy.dead_load_mean_work_j, work, 1e-10);
            close(
                result.energy.external_mean_work_j(),
                storage + kinetic,
                1e-10,
            );
            close(result.energy.parameter_work_j, 0.0, 1e-12);
            close(result.energy.eigenstrain_work_j, 0.0, 1e-12);
            close(result.energy.fracture_dissipation_j, 0.0, 1e-12);
            close(result.normal_contact_work_j, 0.0, 1e-10);
            assert!(result.absolute_contact_work_j <= 1e-10);
            close(result.contact_work_balance_error_j, 0.0, 1e-10);
            close(result.absolute_contact_couple_n_m, 0.0, 1e-10);
            close(result.energy.total_balance_error_j, 0.0, 1e-10);
            close(
                result.native_energy.dead_load_mean_work_j,
                result.energy.dead_load_mean_work_j + result.normal_contact_work_j,
                1e-10,
            );
            assert_eq!(result.energy.inertia.numerical_loss_j, 0.0);
            retained_matter_and_failure(motion, &result.motion);
            assert!(core::ptr::eq(
                motion.solid().materials().as_ptr(),
                result.motion.solid().materials().as_ptr()
            ));
            assert!(minimum_gap(&result.motion) >= 0.0);
            assert!(minimum_gap(&result.motion) <= sustained_controls().activation_tolerance_m);
            assert_eq!(result.motion.time_s(), motion.time_s() + h);
            assert_eq!(result.motion.accepted_steps(), 4);
            assert_eq!(result.motion.remaining_steps(), 2);
            assert_eq!(motion.accepted_steps(), 3);
        });
    }
}

#[test]
fn explicit_temperature_change_generates_contact_pressure_without_healing_or_changing_matter() {
    with_closed_thermal_motion(true, 6, |motion, warm, cx| {
        let geometry = motion.solid().geometry();
        // XY coordinates and both outer apex Z coordinates stay fixed.
        // The six face-normal coordinates are determined by native contact.
        let mut supports = Vec::new();
        for (node, reference) in geometry.nodes_m().iter().enumerate() {
            for axis in 0..3 {
                if axis < 2 || reference[2] != 0.0 {
                    supports.push((
                        3 * node + axis,
                        motion.solid().solution().displacement_m[node][axis],
                    ));
                }
            }
        }
        let heated = motion
            .advance_released_face_sustained_at_materials(
                warm,
                INTERVAL,
                &supports,
                &[],
                settings(),
                sustained_controls(),
                cx,
            )
            .unwrap();
        assert!(core::ptr::eq(
            heated.motion.solid().materials().as_ptr(),
            warm.as_ptr()
        ));
        close(
            heated.motion.solid().materials()[0]
                .band()
                .reference_density_kg_m3()
                .unwrap(),
            4.5,
            1e-14,
        );
        // Holding total strain while free strain rises by .005 gives
        // sigma_zz=-6 Pa and mean sigma_zz=-3 Pa. V*B^T*sigma loads the
        // origin/apex Z DOFs; x/y face corners have zero Z gradient.
        // Hence contact is .5 N at the origin row, zero at the other rows.
        for (row, force) in heated.rows.iter().zip(&heated.solve.forces_n) {
            let expected = if geometry.nodes_m()[row.nodes[0]] == [0.0; 3] {
                0.5
            } else {
                0.0
            };
            close(force[0], expected, 1e-8);
            assert_eq!(force[1], 0.0);
            assert_eq!(force[2], 0.0);
        }
        for state in &heated.motion.solid().solution().element_states {
            close(state.strain()[0], 0.005, 1e-11);
            close(state.strain()[1], 0.005, 1e-11);
            close(state.strain()[2], -0.005, 1e-11);
        }
        for stress in &heated.motion.solid().solution().element_stress {
            close(stress[0], 6.0, 1e-8);
            close(stress[1], 6.0, 1e-8);
            close(stress[2], -6.0, 1e-8);
        }
        // Total V=1/3: transverse U falls .04->.01 J, while compression
        // stores .005 J. Eigenstrain work supplies the signed -.025 J change.
        close(heated.energy.stored_before_j, 0.04, 1e-12);
        close(heated.energy.stored_after_j, 0.015, 1e-10);
        close(heated.energy.eigenstrain_work_j, -0.025, 1e-10);
        close(heated.energy.constitutive_strain_work_j, -0.025, 1e-10);
        close(heated.energy.mechanical_work_j, 0.0, 1e-10);
        close(heated.energy.parameter_work_j, 0.0, 1e-12);
        close(heated.energy.fracture_dissipation_j, 0.0, 1e-12);
        close(heated.energy.external_mean_work_j(), 0.0, 1e-10);
        close(heated.energy.inertia.kinetic_change_j, 0.0, 1e-10);
        close(heated.energy.total_balance_error_j, 0.0, 1e-10);
        close(heated.normal_contact_work_j, 0.0, 1e-10);
        assert!(minimum_gap(&heated.motion) >= 0.0);
        retained_matter_and_failure(motion, &heated.motion);
        assert_eq!(heated.motion.accepted_steps(), 4);

        let cooled = heated
            .motion
            .advance_released_face_sustained_at_materials(
                motion.solid().materials(),
                INTERVAL,
                &supports,
                &[],
                settings(),
                sustained_controls(),
                cx,
            )
            .unwrap();
        assert!(core::ptr::eq(
            cooled.motion.solid().materials().as_ptr(),
            motion.solid().materials().as_ptr()
        ));
        close(cooled.energy.stored_after_j, 0.04, 1e-10);
        close(cooled.energy.eigenstrain_work_j, 0.025, 1e-10);
        close(
            heated.energy.eigenstrain_work_j + cooled.energy.eigenstrain_work_j,
            0.0,
            1e-10,
        );
        close(cooled.energy.external_mean_work_j(), 0.0, 1e-10);
        close(cooled.energy.total_balance_error_j, 0.0, 1e-10);
        close(cooled.energy.fracture_dissipation_j, 0.0, 1e-12);
        retained_matter_and_failure(motion, &cooled.motion);
        assert_eq!(cooled.motion.accepted_steps(), 5);
        assert_eq!(cooled.motion.remaining_steps(), 1);
        assert_eq!(motion.accepted_steps(), 3);
        close(motion.solid().stored_energy_j(), 0.04, 1e-13);
    });
}

#[test]
fn thermal_contact_refusals_preserve_the_original_and_outward_loading_releases_the_faces() {
    with_closed_thermal_motion(false, 7, |motion, _warm, cx| {
        let h = INTERVAL;
        let forces = body_forces(motion, 1.0);
        let snapshot = motion.solid().solution().clone();
        let velocity = motion.velocity_mps().to_vec();
        // Opening admission is deferred only inside unpublished signed-force
        // responses. The public thermal midpoint must still refuse overlap.
        assert!(matches!(
            motion.advance_midpoint(
                motion.solid().materials(),
                h,
                &[],
                &forces,
                settings(),
                midpoint_controls(),
                cx,
            ),
            Err(BrittleTetDynamicError::Solid(
                BrittleTetError::ReleasedFaceClosure { .. }
            ))
        ));
        let mut limited = sustained_controls();
        limited.solve.maximum_response_evaluations = 1;
        let error = motion
            .advance_released_face_sustained(h, &[], &forces, settings(), limited, cx)
            .unwrap_err();
        let BrittleFaceSustainedError::Solve(error) = error else {
            panic!("expected original force-response budget");
        };
        assert!(matches!(
            *error,
            SolveError::ResponseBudget { evaluations: 1, .. }
        ));
        let gate = CancelGate::new_clock_free();
        gate.request();
        ArenaPool::new(ArenaConfig::default()).scope(|arena| {
            let cancelled = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 600,
                    kernel_id: 600,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            assert!(
                motion
                    .advance_released_face_sustained(
                        h,
                        &[],
                        &forces,
                        settings(),
                        sustained_controls(),
                        &cancelled,
                    )
                    .is_err()
            );
        });
        assert_eq!(
            motion.solid().solution().displacement_m,
            snapshot.displacement_m
        );
        assert_eq!(
            motion.solid().solution().element_states,
            snapshot.element_states
        );
        assert_eq!(motion.velocity_mps(), velocity);
        assert_eq!(motion.accepted_steps(), 3);

        let loaded = motion
            .advance_released_face_sustained(h, &[], &forces, settings(), sustained_controls(), cx)
            .unwrap();
        assert!(loaded.solve.forces_n.iter().all(|f| f[0] > 0.0));
        let released = loaded
            .motion
            .advance_released_face_sustained(
                h,
                &[],
                &body_forces(&loaded.motion, -1.0),
                settings(),
                sustained_controls(),
                cx,
            )
            .unwrap();
        assert!(released.solve.forces_n.iter().all(|f| *f == [0.0; 3]));
        close(minimum_gap(&released.motion), h * h, 1e-10);
        for (element, tet) in motion.solid().geometry().tetrahedra().iter().enumerate() {
            let sign = if element == 0 { 1.0 } else { -1.0 };
            for &node in tet {
                let expected = if motion.solid().geometry().nodes_m()[node][2] == 0.0 {
                    -sign * h
                } else {
                    sign * 1.5 * h
                };
                close(released.motion.velocity_mps()[node][2], expected, 1e-8);
            }
        }
        let work = -0.25 * h * h;
        close(released.energy.dead_load_mean_work_j, work, 1e-10);
        close(released.energy.inertia.kinetic_change_j, work, 1e-10);
        close(
            inventory(&released.motion).2 - inventory(&loaded.motion).2,
            work,
            1e-10,
        );
        close(released.energy.stored_after_j, 0.04, 1e-10);
        close(released.energy.fracture_dissipation_j, 0.0, 1e-12);
        close(released.energy.parameter_work_j, 0.0, 1e-12);
        close(released.energy.eigenstrain_work_j, 0.0, 1e-12);
        close(released.energy.total_balance_error_j, 0.0, 1e-10);
        assert_eq!(released.energy.inertia.numerical_loss_j, 0.0);
        assert_eq!(released.normal_contact_work_j, 0.0);
        assert_eq!(released.absolute_contact_work_j, 0.0);
        retained_matter_and_failure(motion, &released.motion);
        assert_eq!(released.motion.accepted_steps(), 5);
        assert_eq!(released.motion.remaining_steps(), 2);
        assert_eq!(released.motion.time_s(), motion.time_s() + 2.0 * h);
    });
}
