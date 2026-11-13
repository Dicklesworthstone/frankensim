//! Several film faces retain one original material inertia and accepted motion.
use super::*;
use fs_couple::material_immersed::dry::{
    MaterialFilmCarrierSi, MaterialFilmPatchSi, ResolvedImmersedSharedSurfaces,
};

type Shared<'a> = ResolvedImmersedSharedSurfaces<'a, ResolvedLiquidState, 1>;

fn carrier_source(material: &ResolvedBody) -> MaterialFilmCarrierSi<'_> {
    let MaterialFilmSurfaceSi::Finite { state, .. } = finite(material) else {
        unreachable!()
    };
    MaterialFilmCarrierSi { material, state }
}

fn patches() -> [MaterialFilmPatchSi; 2] {
    [
        (POINT_BODY, [-2.0, 0.0, 0.0]),
        ([-0.0101, 0.0051, 0.001], [0.0, -3.0, 0.0]),
    ]
    .map(
        |(point_body_m, normal_body)| MaterialFilmPatchSi::Attached {
            carrier: 0,
            point_body_m,
            normal_body,
            cutoff_gap_m: 0.003,
            roughness_gap_m: 1e-5,
            tangential_gaps_m: None,
        },
    )
}

fn corner<'a>(
    liquid: &'a ResolvedLiquidState,
    solid: &'a ResolvedBody,
    material: &'a ResolvedBody,
) -> Shared<'a> {
    let mut source = carrier_source(material);
    source.state.velocity_m_s = [-8e-5, -4e-5, 0.0];
    sphere(liquid, solid)
        .with_shared_surfaces(&[source], &patches(), &CancelGate::new())
        .unwrap()
}

fn kinetic(model: &Shared<'_>) -> f64 {
    let body = model.carriers_lattice()[0];
    model
        .assembly()
        .energy_j(
            0.5 * body.mass * vector(body.velocity).dot(vector(body.velocity))
                + 0.5
                    * body
                        .angular_momentum_body
                        .dot(body.inertia.velocity(body.angular_momentum_body).unwrap()),
        )
        .unwrap()
}

#[test]
fn two_faces_share_material_mobility_pose_and_one_complete_wrench() {
    let liquid = fluid(1000.0, 0.11);
    let solid = body(800.0);
    let material = carrier();
    let mut model = corner(&liquid, &solid, &material);
    let origin = [0.011, -0.013, 0.009];
    let load = [1e-7, -1e-7, 0.0, 1e-11, -2e-11, 0.0];
    assert_eq!(model.carriers_lattice().len(), 1);
    assert_eq!(model.source_patches_si(), patches());
    assert!(core::ptr::eq(
        model.source_carriers_si()[0].material,
        &material
    ));
    let native = model.carriers_lattice()[0];
    close(
        model.assembly().mass_kg(native.mass).unwrap(),
        material.inertia().mass(),
        1e-16,
    );
    for entry in 0..9 {
        close(
            model
                .assembly()
                .inertia_kg_m2(native.inertia.matrix().m[entry])
                .unwrap(),
            material.inertia().inertia_com().m[entry],
            1e-20,
        );
    }
    for step in 1..=3 {
        let before = [
            body_momentum_si(
                model.assembly().state_si(0).unwrap(),
                solid.inertia().mass(),
                origin,
            ),
            body_momentum_si(
                model.carrier_state_si(0).unwrap(),
                material.inertia().mass(),
                origin,
            ),
        ];
        let fluid_before = full_fluid_momentum(model.assembly(), origin);
        let energy_before = kinetic(&model);
        let report = model
            .step_forces_si(
                [[0.0; 6]],
                &[load],
                origin,
                DryBodyCouplingConfig3::default(),
                &CancelGate::new(),
            )
            .unwrap();
        assert_eq!(report.carriers.len(), 1);
        assert_eq!(report.surfaces.len(), 2);
        let shared = &report.carriers[0];
        let carrier_state = model.carrier_state_si(0).unwrap();
        let after = [
            body_momentum_si(
                model.assembly().state_si(0).unwrap(),
                solid.inertia().mass(),
                origin,
            ),
            body_momentum_si(carrier_state, material.inertia().mass(), origin),
        ];
        let fluid_after = full_fluid_momentum(model.assembly(), origin);
        let center = vector(shared.motion.midpoint_center).scale(DX);
        let linear = vector([load[0], load[1], load[2]]).scale(DT);
        let angular =
            vector([load[3], load[4], load[5]]).scale(DT) + (center - vector(origin)).cross(linear);
        let coupling = &report.motion.native.coupling;
        for axis in 0..6 {
            let unit = if axis < 3 {
                model.assembly().impulse_ns(1.0).unwrap()
            } else {
                model.assembly().angular_impulse_nms(1.0).unwrap()
            };
            let external = if axis < 3 {
                array(linear)[axis]
            } else {
                array(angular)[axis - 3]
            } / unit;
            let received: f64 = report
                .surfaces
                .iter()
                .map(|p| {
                    if axis < 3 {
                        p.impulse[axis]
                    } else {
                        p.angular_impulse[axis - 3]
                    }
                })
                .sum();
            let unique = if axis < 3 {
                shared.film.impulse[axis]
            } else {
                shared.film.angular_impulse[axis - 3]
            };
            close(unique, received, 2e-12);
            close(
                (after[1][axis] - before[1][axis]) / unit,
                external + received,
                3e-12,
            );
            let boundary = if axis < 3 {
                coupling.walls.fluid_population_impulse[axis]
            } else {
                coupling.walls.fluid_population_angular_impulse[axis - 3]
                    + coupling.periodic.angular_momentum_change[axis - 3]
            };
            // Count actual carrier momentum once, regardless of face count.
            close(
                fluid_after[axis]
                    - fluid_before[axis]
                    - (report.motion.interior[0][1].momentum[axis]
                        - report.motion.interior[0][0].momentum[axis])
                    + (after[0][axis] - before[0][axis] + after[1][axis] - before[1][axis]) / unit
                    - external
                    - report.motion.physical_external_impulses[0][axis]
                    - report.motion.represented_fluid_forcing_impulse[axis]
                    - boundary,
                0.0,
                4e-10,
            );
        }
        for (index, patch) in patches().into_iter().enumerate() {
            let MaterialFilmPatchSi::Attached {
                point_body_m,
                normal_body,
                cutoff_gap_m,
                roughness_gap_m,
                ..
            } = patch
            else {
                unreachable!()
            };
            let exchange = report.surfaces[index];
            assert_eq!(exchange.carrier, Some(0));
            assert_eq!(exchange.active_spheres, 1);
            let normal_body = vector(normal_body).scale(1.0 / vector(normal_body).norm());
            let midpoint = shared.motion.midpoint_orientation;
            let point = center + midpoint.rotate(vector(point_body_m)).unwrap();
            let normal = midpoint.rotate(normal_body).unwrap();
            let sphere_center = vector(report.motion.native.bodies[0].midpoint_center).scale(DX);
            let distance = (sphere_center - point).dot(normal);
            let foot = sphere_center - normal.scale(distance);
            let response = &shared.film.motion;
            let surface_velocity = vector(response.velocity_midpoint).scale(DX / DT)
                + vector(response.angular_velocity_midpoint)
                    .scale(1.0 / DT)
                    .cross(foot - center);
            let slip = (vector(coupling.bodies[0].velocity_midpoint).scale(DX / DT)
                - surface_velocity)
                .dot(normal);
            let resistance = 6.0
                * core::f64::consts::PI
                * 0.11
                * RADIUS.powi(2)
                * (1.0 / (distance - RADIUS).max(roughness_gap_m) - 1.0 / cutoff_gap_m);
            let expected = array(normal.scale(resistance * slip * DT));
            let (actual_point, actual_normal) = model.surface_geometry_si(index).unwrap();
            let final_point = array(
                vector(carrier_state.center_m)
                    + carrier_state
                        .orientation
                        .rotate(vector(point_body_m))
                        .unwrap(),
            );
            let final_normal = array(carrier_state.orientation.rotate(normal_body).unwrap());
            for axis in 0..3 {
                close(
                    model.assembly().impulse_ns(exchange.impulse[axis]).unwrap(),
                    expected[axis],
                    5e-18,
                );
                close(actual_point[axis], final_point[axis], 1e-16);
                close(actual_normal[axis], final_normal[axis], 1e-14);
            }
        }
        let patch_work: f64 = report.surfaces.iter().map(|p| p.work).sum();
        close(shared.film.work, patch_work, 2e-12);
        close(
            coupling.lubrication.body_work[0] + patch_work + coupling.lubrication.dissipated_energy,
            0.0,
            2e-12,
        );
        close(
            kinetic(&model) - energy_before,
            model
                .assembly()
                .energy_j(shared.motion.kinetic_energy_change)
                .unwrap(),
            1e-20,
        );
        close(shared.motion.work_residual, 0.0, 2e-12);
        let probe = [0.035, 0.044, 0.041];
        assert_eq!(
            model.surface_velocity_at_si(0, probe).unwrap(),
            model.surface_velocity_at_si(1, probe).unwrap()
        );
        assert_eq!(model.assembly().accepted_steps(), step);
    }
    assert_ne!(
        model.carrier_state_si(0).unwrap().orientation,
        So3::identity()
    );
    assert_eq!(model.source_patches_si(), patches());
}

#[test]
fn one_patch_per_carrier_and_empty_registry_reproduce_existing_owners() {
    let liquid = fluid(1000.0, 0.11);
    let solid = body(800.0);
    let material = carrier();
    let gate = CancelGate::new();
    let controls = DryBodyCouplingConfig3::default();
    let mut original = mixed(&liquid, &solid, &material);
    let mut shared = sphere(&liquid, &solid)
        .with_shared_surfaces(
            &[carrier_source(&material)],
            &[
                MaterialFilmPatchSi::Prescribed {
                    plane: prescribed(),
                    axial_spin_rad_s: 0.0,
                },
                patches()[0],
            ],
            &gate,
        )
        .unwrap();
    let load = [1e-7, 2e-8, 0.0, 0.0, -1e-11, 0.0];
    for _ in 0..2 {
        let a = original
            .step_forces_si([[0.0; 6]], &[[0.0; 6], load], [0.0; 3], controls, &gate)
            .unwrap();
        let b = shared
            .step_forces_si([[0.0; 6]], &[load], [0.0; 3], controls, &gate)
            .unwrap();
        assert_eq!(a.motion, b.motion);
        assert_eq!(a.surfaces[1].finite.unwrap(), b.carriers[0].motion);
        assert_eq!(
            a.surfaces[1].film.finite_body.unwrap(),
            b.carriers[0].film.motion
        );
        for i in 0..2 {
            assert_eq!(a.surfaces[i].film.impulse, b.surfaces[i].impulse);
            assert_eq!(
                a.surfaces[i].film.angular_impulse,
                b.surfaces[i].angular_impulse
            );
            assert_eq!(a.surfaces[i].film.work, b.surfaces[i].work);
            assert_eq!(
                original.surface_geometry_si(i).unwrap(),
                shared.surface_geometry_si(i).unwrap()
            );
        }
        assert_eq!(
            original.surface_body_state_si(1).unwrap().unwrap(),
            shared.carrier_state_si(0).unwrap()
        );
        assert_eq!(fields(original.assembly()), fields(shared.assembly()));
        assert_eq!(
            original.assembly().bodies_lattice(),
            shared.assembly().bodies_lattice()
        );
    }
    let mut plain = sphere(&liquid, &solid);
    let mut empty = sphere(&liquid, &solid)
        .with_shared_surfaces(&[], &[], &gate)
        .unwrap();
    let a = plain
        .step_forces_si([[0.0; 6]], [0.0; 3], controls, &gate)
        .unwrap();
    let b = empty
        .step_forces_si([[0.0; 6]], &[], [0.0; 3], controls, &gate)
        .unwrap();
    assert_eq!(a, b.motion);
    assert!(b.carriers.is_empty());
    assert!(b.surfaces.is_empty());
    assert_eq!(fields(&plain), fields(empty.assembly()));
}

#[test]
fn shared_carrier_bindings_rebind_and_failed_trials_keep_exact_retry() {
    let first = fluid(1000.0, 0.11);
    let second = fluid(1000.0, 0.22);
    let solid = body(800.0);
    let material = carrier();
    let gate = CancelGate::new();
    let source = carrier_source(&material);
    let mut invalid = patches();
    let MaterialFilmPatchSi::Attached { carrier, .. } = &mut invalid[1] else {
        unreachable!()
    };
    *carrier = 1;
    assert!(
        sphere(&first, &solid)
            .with_shared_surfaces(&[source], &invalid, &gate)
            .is_err()
    );
    assert!(
        sphere(&first, &solid)
            .with_shared_surfaces(&[source], &[], &gate)
            .is_err()
    );
    assert!(
        sphere(&first, &solid)
            .with_shared_surfaces(&[source; 2], &patches(), &gate)
            .is_err()
    );
    // Reusing a material card is still two bodies when two carrier slots are
    // explicitly declared, with one patch attached to each slot.
    let separate = sphere(&first, &solid)
        .with_shared_surfaces(&[source; 2], &invalid, &gate)
        .unwrap();
    assert_eq!(separate.carriers_lattice().len(), 2);
    let mut actual = corner(&first, &solid, &material);
    let mut reference = corner(&first, &solid, &material);
    let controls = DryBodyCouplingConfig3::default();
    let load = [[1e-7, 2e-8, 0.0, 0.0, -1e-11, 0.0]];
    assert_eq!(
        actual
            .step_forces_si([[0.0; 6]], &load, [0.0; 3], controls, &gate)
            .unwrap(),
        reference
            .step_forces_si([[0.0; 6]], &load, [0.0; 3], controls, &gate)
            .unwrap()
    );
    let populations = fields(actual.assembly());
    let bodies = *actual.assembly().bodies_lattice();
    let carriers = actual.carriers_lattice().to_vec();
    let surfaces = actual.surfaces_lattice().to_vec();
    actual.rebind_fluid(&second, &gate).unwrap();
    reference.rebind_fluid(&second, &gate).unwrap();
    assert!(core::ptr::eq(
        actual.assembly().scaling().material(),
        &second
    ));
    let limited = DryBodyCouplingConfig3 {
        maximum_evaluations: 1,
        ..controls
    };
    assert!(
        actual
            .step_forces_si([[0.0; 6]], &load, [0.0; 3], limited, &gate)
            .is_err()
    );
    let stopped = CancelGate::new();
    stopped.request();
    assert!(
        actual
            .step_forces_si([[0.0; 6]], &load, [0.0; 3], controls, &stopped)
            .is_err()
    );
    for bad in [&[][..], &[[f64::NAN; 6]][..], &[[0.0; 6]; 2][..]] {
        assert!(
            actual
                .step_forces_si([[0.0; 6]], bad, [0.0; 3], controls, &gate)
                .is_err()
        );
    }
    assert_eq!(fields(actual.assembly()), populations);
    assert_eq!(actual.assembly().bodies_lattice(), &bodies);
    assert_eq!(actual.carriers_lattice(), carriers);
    assert_eq!(actual.surfaces_lattice(), surfaces);
    assert_eq!(actual.source_patches_si(), patches());
    assert_eq!(actual.assembly().accepted_steps(), 1);
    assert_eq!(
        actual
            .step_forces_si([[0.0; 6]], &load, [0.0; 3], controls, &gate)
            .unwrap(),
        reference
            .step_forces_si([[0.0; 6]], &load, [0.0; 3], controls, &gate)
            .unwrap()
    );
    assert_eq!(actual.carriers_lattice(), reference.carriers_lattice());
    assert_eq!(actual.surfaces_lattice(), reference.surfaces_lattice());
    assert_eq!(fields(actual.assembly()), fields(reference.assembly()));
    assert_eq!(actual.assembly().accepted_steps(), 2);
}

#[test]
fn live_open_targets_preserve_shared_state_and_enter_the_common_momentum_account() {
    use fs_couple::material_immersed::FaceBoundarySi;
    use fs_lbm::d3q19::{Face3, FaceBoundary3};

    let liquid = fluid(1000.0, 0.11);
    let solid = body(800.0);
    let material = carrier();
    let gate = CancelGate::new();
    let controls = DryBodyCouplingConfig3::default();
    let faces = |speed, pressure| {
        let mut faces = [FaceBoundarySi::Periodic; 6];
        faces[0] = FaceBoundarySi::Velocity {
            velocity_m_s: [speed, 0.0, 0.0],
        };
        faces[1] = FaceBoundarySi::Pressure {
            perturbation_pa: pressure,
        };
        faces
    };
    let make = || {
        let mut state = initial()[0];
        state.velocity_m_s = [0.0; 3];
        let mut geometry = film();
        geometry.maximum_pairs = 2;
        let mut source = carrier_source(&material);
        source.state.velocity_m_s = [-8e-5, -4e-5, 0.0];
        Model::spheres(
            &liquid,
            domain().with_boundaries(faces(0.0, 0.0)),
            [SubmergedSphere::new(&solid, RADIUS).unwrap()],
            [state],
            [0.0; 3],
            geometry,
            &gate,
        )
        .unwrap()
        .with_shared_surfaces(&[source], &patches(), &gate)
        .unwrap()
    };
    let snapshot = |model: &Shared<'_>| {
        (
            fields(model.assembly()),
            *model.assembly().bodies_lattice(),
            model.carriers_lattice().to_vec(),
            model.surfaces_lattice().to_vec(),
            model.assembly().accepted_steps(),
        )
    };
    let origin = [0.011, -0.013, 0.009];
    let mut actual = make();
    let mut reference = make();
    for model in [&mut actual, &mut reference] {
        model
            .step_forces_si([[0.0; 6]], &[[0.0; 6]], origin, controls, &gate)
            .unwrap();
    }
    let accepted = snapshot(&actual);
    assert_eq!(accepted.4, 1);
    let targets = faces(2e-4, 2e-5);
    actual.rebind_open_boundaries_si(targets, &gate).unwrap();
    reference.rebind_open_boundaries_si(targets, &gate).unwrap();
    assert_eq!(snapshot(&actual), accepted);
    let native_targets = actual.assembly().grid().boundaries();
    let FaceBoundary3::Velocity { velocity } = native_targets.face(Face3::XMin) else {
        panic!("velocity reservoir lost");
    };
    close(velocity[0], 0.002, 1e-18);
    let FaceBoundary3::Pressure { density } = native_targets.face(Face3::XMax) else {
        panic!("pressure reservoir lost");
    };
    close(
        density,
        1.0 + 3.0 * 2e-5 * DT * DT / (1000.0 * DX * DX),
        1e-15,
    );

    let stopped = CancelGate::new();
    stopped.request();
    assert!(
        actual
            .rebind_open_boundaries_si(faces(0.0, 0.0), &stopped)
            .is_err()
    );
    assert!(
        actual
            .rebind_open_boundaries_si(faces(0.0, f64::NAN), &gate)
            .is_err()
    );
    assert!(
        actual
            .rebind_open_boundaries_si([FaceBoundarySi::Periodic; 6], &gate)
            .is_err()
    );
    assert!(
        actual
            .step_forces_si(
                [[0.0; 6]],
                &[[0.0; 6]],
                origin,
                DryBodyCouplingConfig3 {
                    maximum_evaluations: 1,
                    ..controls
                },
                &gate,
            )
            .is_err()
    );
    assert_eq!(snapshot(&actual), accepted);
    assert_eq!(actual.assembly().grid().boundaries(), native_targets);

    let before = [
        body_momentum_si(
            actual.assembly().state_si(0).unwrap(),
            solid.inertia().mass(),
            origin,
        ),
        body_momentum_si(
            actual.carrier_state_si(0).unwrap(),
            material.inertia().mass(),
            origin,
        ),
    ];
    let fluid_before = full_fluid_momentum(actual.assembly(), origin);
    let report = actual
        .step_forces_si([[0.0; 6]], &[[0.0; 6]], origin, controls, &gate)
        .unwrap();
    assert_eq!(
        report,
        reference
            .step_forces_si([[0.0; 6]], &[[0.0; 6]], origin, controls, &gate)
            .unwrap()
    );
    assert_eq!(snapshot(&actual), snapshot(&reference));
    assert_eq!(actual.assembly().accepted_steps(), 2);
    let coupling = &report.motion.native.coupling;
    assert!(coupling.open.reconstructed_cells > 0);
    assert!(coupling.open.impulse[0].abs() > 1e-6);
    assert_eq!(report.carriers.len(), 1);
    assert_eq!(report.surfaces.len(), 2);
    let after = [
        body_momentum_si(
            actual.assembly().state_si(0).unwrap(),
            solid.inertia().mass(),
            origin,
        ),
        body_momentum_si(
            actual.carrier_state_si(0).unwrap(),
            material.inertia().mass(),
            origin,
        ),
    ];
    let fluid_after = full_fluid_momentum(actual.assembly(), origin);
    for axis in 0..6 {
        let unit = if axis < 3 {
            actual.assembly().impulse_ns(1.0).unwrap()
        } else {
            actual.assembly().angular_impulse_nms(1.0).unwrap()
        };
        let boundary = if axis < 3 {
            coupling.open.impulse[axis] + coupling.walls.fluid_population_impulse[axis]
        } else {
            coupling.open.angular_impulse[axis - 3]
                + coupling.walls.fluid_population_angular_impulse[axis - 3]
                + coupling.periodic.angular_momentum_change[axis - 3]
        };
        // Actual exterior fluid plus dry sphere plus the ONE carrier receives
        // reservoir momentum once; patch reactions are internal to this sum.
        close(
            fluid_after[axis]
                - fluid_before[axis]
                - (report.motion.interior[0][1].momentum[axis]
                    - report.motion.interior[0][0].momentum[axis])
                + (after[0][axis] - before[0][axis] + after[1][axis] - before[1][axis]) / unit
                - boundary,
            0.0,
            4e-10,
        );
    }
}
