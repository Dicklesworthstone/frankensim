//! Material inertia and actual carrier motion share the dry-sphere film solve.
use super::*;
use fs_couple::material_immersed::dry::{
    MaterialFilmSurfaceSi, PlaneFilmSi, ResolvedImmersedSurfaces,
};
use fs_lbm::d3q19::{DryLubricatedSurfacesExchange3, MovingLubricationSurface3};

type Surfaces<'a> = ResolvedImmersedSurfaces<'a, ResolvedLiquidState, 1>;
const POINT_BODY: [f64; 3] = [-0.003, -0.002, 0.001];

fn carrier() -> ResolvedBody {
    let rotation = So3::exp(So3Tangent::new(Vec3::new(0.35, -0.2, 0.25))).unwrap();
    let geometry = VolumeMoments::cuboid(Vec3::new(0.006, 0.012, 0.008))
        .unwrap()
        .transformed(Se3::from_parts(rotation, Vec3::default()).unwrap())
        .unwrap();
    resolve_body(
        &[RegionBinding {
            key: "film carrier",
            geometry,
            material: &card(
                "film carrier",
                &[("density", Density::DIMS, 1800.0)],
                "solid",
            ),
            point: &point(),
            selection: MaterialPropertySelection::SingleClaimOnly,
        }],
        1,
        || false,
    )
    .unwrap()
}

fn prescribed() -> PlaneFilmSi {
    PlaneFilmSi {
        point_m: [0.0229, 0.04, 0.04],
        normal: [2.0, 0.0, 0.0],
        velocity_m_s: [0.0001, 0.0, 0.0],
        cutoff_gap_m: 0.003,
        roughness_gap_m: 1e-5,
        tangential_gaps_m: None,
    }
}

fn finite(material: &ResolvedBody) -> MaterialFilmSurfaceSi<'_> {
    MaterialFilmSurfaceSi::Finite {
        material,
        state: SphereStateSi {
            center_m: [0.0401, 0.042, 0.039],
            orientation: So3::identity(),
            velocity_m_s: [0.0; 3],
            angular_momentum_body_kg_m2_s: [0.0; 3],
        },
        point_body_m: POINT_BODY,
        normal_body: [-2.0, 0.0, 0.0],
        cutoff_gap_m: 0.003,
        roughness_gap_m: 1e-5,
        tangential_gaps_m: None,
    }
}

fn sphere<'a>(liquid: &'a ResolvedLiquidState, solid: &'a ResolvedBody) -> Model<'a, 1> {
    let mut state = initial()[0];
    state.velocity_m_s = [0.0; 3];
    let mut geometry = film();
    geometry.maximum_pairs = 2;
    Model::spheres(
        liquid,
        domain(),
        [SubmergedSphere::new(solid, RADIUS).unwrap()],
        [state],
        [0.0; 3],
        geometry,
        &CancelGate::new(),
    )
    .unwrap()
}

fn mixed<'a>(
    liquid: &'a ResolvedLiquidState,
    solid: &'a ResolvedBody,
    material: &'a ResolvedBody,
) -> Surfaces<'a> {
    sphere(liquid, solid)
        .with_moving_surfaces(
            &[
                MaterialFilmSurfaceSi::Prescribed {
                    plane: prescribed(),
                    axial_spin_rad_s: 0.0,
                },
                finite(material),
            ],
            &CancelGate::new(),
        )
        .unwrap()
}

fn carrier_energy(model: &Surfaces<'_>, index: usize) -> f64 {
    let MovingLubricationSurface3::Finite(surface) = model.surfaces_lattice()[index] else {
        panic!("expected a finite carrier")
    };
    let body = surface.body();
    let omega = body.inertia.velocity(body.angular_momentum_body).unwrap();
    model
        .assembly()
        .energy_j(
            0.5 * body.mass * vector(body.velocity).dot(vector(body.velocity))
                + 0.5 * body.angular_momentum_body.dot(omega),
        )
        .unwrap()
}

fn mixed_balance(
    model: &Surfaces<'_>,
    before: [SphereStateSi; 2],
    fluid_before: [f64; 6],
    carrier_mass: f64,
    origin: [f64; 3],
    report: &DryLubricatedSurfacesExchange3<1>,
) {
    let assembly = model.assembly();
    let motion = &report.motion;
    let fluid_after = full_fluid_momentum(assembly, origin);
    let sphere_mass = assembly.materials()[0].inertia().mass();
    let old = [
        body_momentum_si(before[0], sphere_mass, origin),
        body_momentum_si(before[1], carrier_mass, origin),
    ];
    let new = [
        body_momentum_si(assembly.state_si(0).unwrap(), sphere_mass, origin),
        body_momentum_si(
            model.surface_body_state_si(1).unwrap().unwrap(),
            carrier_mass,
            origin,
        ),
    ];
    let finite = report.surfaces[1].finite.unwrap();
    let coupling = &motion.native.coupling;
    for axis in 0..6 {
        let unit = if axis < 3 {
            assembly.impulse_ns(1.0).unwrap()
        } else {
            assembly.angular_impulse_nms(1.0).unwrap()
        };
        let prescribed = if axis < 3 {
            report.surfaces[0].film.impulse[axis]
        } else {
            report.surfaces[0].film.angular_impulse[axis - 3]
        };
        let boundary = if axis < 3 {
            coupling.walls.fluid_population_impulse[axis]
        } else {
            coupling.walls.fluid_population_angular_impulse[axis - 3]
                + coupling.periodic.angular_momentum_change[axis - 3]
        };
        // Actual exterior-fluid + dry-sphere + carrier momentum: only the
        // prescribed reaction belongs to an external mechanical reservoir.
        let balance = fluid_after[axis]
            - fluid_before[axis]
            - (motion.interior[0][1].momentum[axis] - motion.interior[0][0].momentum[axis])
            + (new[0][axis] - old[0][axis] + new[1][axis] - old[1][axis]) / unit
            + prescribed
            - motion.physical_external_impulses[0][axis]
            - finite.external_impulse[axis]
            - motion.represented_fluid_forcing_impulse[axis]
            - boundary;
        close(balance, 0.0, 4e-10);
        close(
            (new[1][axis] - old[1][axis]) / unit,
            finite.momentum_change[axis],
            2e-12,
        );
        close(finite.momentum_residual[axis], 0.0, 2e-12);
    }
}

#[test]
fn material_carrier_preserves_full_inertia_current_clock_and_si_wrench() {
    let liquid = fluid(1000.0, 0.11);
    let solid = body(800.0);
    let material = carrier();
    let gate = CancelGate::new();
    let controls = DryBodyCouplingConfig3::default();
    let mut assembly = sphere(&liquid, &solid);
    assembly
        .step_forces_si([[0.0; 6]], [0.0; 3], controls, &gate)
        .unwrap();
    let populations = fields(&assembly);
    let bodies = *assembly.bodies_lattice();
    let mut source = finite(&material);
    let MaterialFilmSurfaceSi::Finite { state, .. } = &mut source else {
        unreachable!()
    };
    state.center_m[0] = 0.06; // no film: isolate the complete authored carrier wrench
    state.orientation = So3::exp(So3Tangent::new(Vec3::new(0.1, 0.2, -0.15))).unwrap();
    state.velocity_m_s = [1e-5, -2e-5, 1e-5];
    state.angular_momentum_body_kg_m2_s = [1e-12, -2e-12, 3e-12];
    let initial = *state;
    let mut model = assembly.with_moving_surfaces(&[source], &gate).unwrap();
    assert_eq!(model.assembly().accepted_steps(), 1);
    assert_eq!(fields(model.assembly()), populations);
    assert_eq!(model.assembly().bodies_lattice(), &bodies);
    let MaterialFilmSurfaceSi::Finite {
        material: retained, ..
    } = model.source_surfaces_si()[0]
    else {
        unreachable!()
    };
    assert!(core::ptr::eq(retained, &material));
    let MovingLubricationSurface3::Finite(surface) = model.surfaces_lattice()[0] else {
        unreachable!()
    };
    let native = surface.body();
    close(
        model.assembly().mass_kg(native.mass).unwrap(),
        material.inertia().mass(),
        1e-16,
    );
    assert!(material.inertia().inertia_com().m[1].abs() > 1e-10);
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
    let returned = model.surface_body_state_si(0).unwrap().unwrap();
    assert_eq!(returned.orientation, initial.orientation);
    for axis in 0..3 {
        close(returned.center_m[axis], initial.center_m[axis], 1e-16);
        close(
            returned.velocity_m_s[axis],
            initial.velocity_m_s[axis],
            1e-18,
        );
        close(
            returned.angular_momentum_body_kg_m2_s[axis],
            initial.angular_momentum_body_kg_m2_s[axis],
            1e-25,
        );
    }
    let probe = [0.055, 0.039, 0.044];
    let omega_body = fs_time::lie::gyrostat::InertiaTensor3::new(material.inertia().inertia_com())
        .unwrap()
        .velocity(vector(initial.angular_momentum_body_kg_m2_s))
        .unwrap();
    let omega_world = initial.orientation.rotate(omega_body).unwrap();
    let expected_velocity = vector(initial.velocity_m_s)
        + omega_world.cross(vector(probe) - vector(initial.center_m));
    for (actual, expected) in model
        .surface_velocity_at_si(0, probe)
        .unwrap()
        .into_iter()
        .zip(array(expected_velocity))
    {
        close(actual, expected, 1e-17);
    }
    let energy_before = carrier_energy(&model, 0);
    let origin = [0.011, -0.013, 0.009];
    let loads = [2e-7, -3e-7, 1e-7, 1e-11, -2e-11, 3e-11];
    let report = model
        .step_forces_si([[0.0; 6]], &[loads], origin, controls, &gate)
        .unwrap();
    let result = report.surfaces[0].finite.unwrap();
    assert_eq!(report.surfaces[0].film.active_spheres, 0);
    let old = body_momentum_si(initial, material.inertia().mass(), origin);
    let new = body_momentum_si(
        model.surface_body_state_si(0).unwrap().unwrap(),
        material.inertia().mass(),
        origin,
    );
    let impulse = vector([loads[0], loads[1], loads[2]]).scale(DT);
    let arm = vector(result.midpoint_center).scale(DX) - vector(origin);
    let angular = vector([loads[3], loads[4], loads[5]]).scale(DT) + arm.cross(impulse);
    for axis in 0..6 {
        let expected = if axis < 3 {
            array(impulse)[axis]
        } else {
            array(angular)[axis - 3]
        };
        close(new[axis] - old[axis], expected, 1e-20);
    }
    close(
        carrier_energy(&model, 0) - energy_before,
        model
            .assembly()
            .energy_j(result.kinetic_energy_change)
            .unwrap(),
        1e-21,
    );
    close(result.work_residual, 0.0, 2e-12);
    assert_eq!(model.assembly().accepted_steps(), 2);
}

#[test]
fn mixed_surfaces_recoil_rotate_and_close_actual_momentum_and_work_over_repeated_steps() {
    let liquid = fluid(1000.0, 0.11);
    let solid = body(800.0);
    let material = carrier();
    let mut model = mixed(&liquid, &solid, &material);
    let origin = [0.011, -0.013, 0.009];
    let mut observed_carrier_work = 0.0;
    for step in 1..=4 {
        let before = [
            model.assembly().state_si(0).unwrap(),
            model.surface_body_state_si(1).unwrap().unwrap(),
        ];
        let fluid_before = full_fluid_momentum(model.assembly(), origin);
        let energy_before = carrier_energy(&model, 1);
        let report = model
            .step_forces_si(
                [[0.0; 6]],
                &[[0.0; 6]; 2],
                origin,
                DryBodyCouplingConfig3::default(),
                &CancelGate::new(),
            )
            .unwrap();
        mixed_balance(
            &model,
            before,
            fluid_before,
            material.inertia().mass(),
            origin,
            &report,
        );
        assert!(report.surfaces[0].finite.is_none());
        assert!(model.surface_body_state_si(0).unwrap().is_none());
        assert_eq!(report.surfaces[0].film.active_spheres, 1);
        assert_eq!(report.surfaces[1].film.active_spheres, 1);
        let film = &report.motion.native.coupling.lubrication;
        assert!(film.dissipated_energy > 0.0);
        close(
            film.body_work[0]
                + report.surfaces.iter().map(|s| s.film.work).sum::<f64>()
                + film.dissipated_energy,
            0.0,
            2e-12,
        );
        let finite = report.surfaces[1].finite.unwrap();
        close(
            carrier_energy(&model, 1) - energy_before,
            model
                .assembly()
                .energy_j(finite.kinetic_energy_change)
                .unwrap(),
            1e-20,
        );
        close(
            finite.kinetic_energy_change - finite.free_rotation_energy_change,
            report.surfaces[1].film.work,
            2e-12,
        );
        close(finite.work_residual, 0.0, 2e-12);
        observed_carrier_work += report.surfaces[1].film.work;
        let state = model.surface_body_state_si(1).unwrap().unwrap();
        let (point, normal) = model.surface_geometry_si(1).unwrap();
        let expected_point =
            vector(state.center_m) + state.orientation.rotate(vector(POINT_BODY)).unwrap();
        let expected_normal = state
            .orientation
            .rotate(Vec3::new(-1.0, 0.0, 0.0))
            .unwrap();
        for axis in 0..3 {
            close(point[axis], array(expected_point)[axis], 1e-16);
            close(normal[axis], array(expected_normal)[axis], 1e-14);
        }
        close(
            model.surface_geometry_si(0).unwrap().0[0],
            prescribed().point_m[0] + step as f64 * DT * prescribed().velocity_m_s[0],
            1e-16,
        );
        assert_eq!(model.assembly().accepted_steps(), step);
    }
    let final_state = model.surface_body_state_si(1).unwrap().unwrap();
    assert!(final_state.velocity_m_s[0] > 0.0);
    assert!(final_state.center_m[0] > 0.0401);
    assert!(vector(final_state.angular_momentum_body_kg_m2_s).norm() > 0.0);
    assert_ne!(model.surface_geometry_si(1).unwrap().1, [-1.0, 0.0, 0.0]);
    assert!(observed_carrier_work > 0.0);
}

#[test]
fn empty_and_prescribed_only_surfaces_keep_existing_trajectories() {
    let liquid = fluid(1000.0, 0.11);
    let solid = body(800.0);
    let gate = CancelGate::new();
    let controls = DryBodyCouplingConfig3::default();
    let mut plain = sphere(&liquid, &solid);
    let mut empty = sphere(&liquid, &solid);
    let loads = [[2e-7, 0.0, 0.0, 0.0, 0.0, 1e-11]];
    assert_eq!(
        plain
            .step_forces_si(loads, [0.0; 3], controls, &gate)
            .unwrap(),
        empty
            .step_forces_si(loads, [0.0; 3], controls, &gate)
            .unwrap()
    );
    let mut empty = empty.with_moving_surfaces(&[], &gate).unwrap();
    let result = empty
        .step_forces_si(loads, &[], [0.0; 3], controls, &gate)
        .unwrap();
    assert_eq!(
        result.motion,
        plain
            .step_forces_si(loads, [0.0; 3], controls, &gate)
            .unwrap()
    );
    assert!(result.surfaces.is_empty());
    assert_eq!(fields(empty.assembly()), fields(&plain));
    assert_eq!(empty.assembly().bodies_lattice(), plain.bodies_lattice());
    assert_eq!(empty.assembly().accepted_steps(), 2);

    let mut original = sphere(&liquid, &solid)
        .with_prescribed_planes(&[prescribed()], &gate)
        .unwrap();
    let mut attached = sphere(&liquid, &solid)
        .with_moving_surfaces(
            &[MaterialFilmSurfaceSi::Prescribed {
                plane: prescribed(),
                axial_spin_rad_s: 0.0,
            }],
            &gate,
        )
        .unwrap();
    for _ in 0..2 {
        let a = original
            .step_forces_si(loads, [0.0; 3], controls, &gate)
            .unwrap();
        let b = attached
            .step_forces_si(loads, &[[0.0; 6]], [0.0; 3], controls, &gate)
            .unwrap();
        assert_eq!(a.motion, b.motion);
        assert_eq!(a.planes[0], b.surfaces[0].film);
        assert_eq!(fields(original.assembly()), fields(attached.assembly()));
        assert_eq!(
            original.assembly().bodies_lattice(),
            attached.assembly().bodies_lattice()
        );
        assert_eq!(
            original.plane_si(0).unwrap().point_m,
            attached.surface_geometry_si(0).unwrap().0
        );
    }
}

#[test]
fn refresh_and_failed_candidates_preserve_all_surface_history_and_exact_retry() {
    let first = fluid(1000.0, 0.11);
    let second = fluid(1000.0, 0.22);
    let wrong_density = fluid(1200.0, 0.22);
    let solid = body(800.0);
    let material = carrier();
    let mut actual = mixed(&first, &solid, &material);
    let mut reference = mixed(&first, &solid, &material);
    let gate = CancelGate::new();
    let controls = DryBodyCouplingConfig3::default();
    let loads = [[2e-7, 0.0, 0.0, 0.0, 0.0, 1e-11]];
    let surface_loads = [[0.0; 6], [-1e-7, 2e-8, 0.0, 0.0, -1e-11, 0.0]];
    assert_eq!(
        actual
            .step_forces_si(loads, &surface_loads, [0.0; 3], controls, &gate)
            .unwrap(),
        reference
            .step_forces_si(loads, &surface_loads, [0.0; 3], controls, &gate)
            .unwrap()
    );
    let populations = fields(actual.assembly());
    let bodies = *actual.assembly().bodies_lattice();
    let surfaces = actual.surfaces_lattice().to_vec();
    actual.rebind_fluid(&second, &gate).unwrap();
    reference.rebind_fluid(&second, &gate).unwrap();
    assert!(core::ptr::eq(
        actual.assembly().scaling().material(),
        &second
    ));
    assert_eq!(actual.surfaces_lattice(), surfaces);
    assert_eq!(fields(actual.assembly()), populations);
    assert_eq!(actual.assembly().bodies_lattice(), &bodies);
    assert!(actual.rebind_fluid(&wrong_density, &gate).is_err());
    let limited = DryBodyCouplingConfig3 {
        maximum_evaluations: 1,
        ..controls
    };
    assert!(
        actual
            .step_forces_si(loads, &surface_loads, [0.0; 3], limited, &gate)
            .is_err()
    );
    let stopped = CancelGate::new();
    stopped.request();
    assert!(
        actual
            .step_forces_si(loads, &surface_loads, [0.0; 3], controls, &stopped)
            .is_err()
    );
    for bad in [
        &[[0.0; 6]][..],
        &[[0.0; 6], [f64::NAN; 6]][..],
        &[[1e-7; 6], [0.0; 6]][..],
    ] {
        assert!(
            actual
                .step_forces_si(loads, bad, [0.0; 3], controls, &gate)
                .is_err()
        );
    }
    assert_eq!(actual.surfaces_lattice(), surfaces);
    assert_eq!(fields(actual.assembly()), populations);
    assert_eq!(actual.assembly().bodies_lattice(), &bodies);
    assert_eq!(actual.assembly().accepted_steps(), 1);
    let accepted = actual
        .step_forces_si(loads, &surface_loads, [0.0; 3], controls, &gate)
        .unwrap();
    assert_eq!(
        accepted,
        reference
            .step_forces_si(loads, &surface_loads, [0.0; 3], controls, &gate)
            .unwrap()
    );
    let MovingLubricationSurface3::Prescribed(before) = surfaces[0] else {
        unreachable!()
    };
    let gap = accepted.motion.native.bodies[0].midpoint_center[0] * DX
        - (before.point()[0] * DX + 0.5 * DT * prescribed().velocity_m_s[0])
        - RADIUS;
    let relative = accepted.motion.native.coupling.bodies[0].velocity_midpoint[0] * DX / DT
        - prescribed().velocity_m_s[0];
    let resistance = 6.0 * core::f64::consts::PI * 0.22 * RADIUS.powi(2)
        * (1.0 / gap.max(prescribed().roughness_gap_m) - 1.0 / prescribed().cutoff_gap_m);
    close(
        actual
            .assembly()
            .impulse_ns(accepted.surfaces[0].film.impulse[0])
            .unwrap(),
        DT * resistance * relative,
        5e-18,
    );
    assert_eq!(actual.surfaces_lattice(), reference.surfaces_lattice());
    assert_eq!(fields(actual.assembly()), fields(reference.assembly()));
    assert_eq!(
        actual.assembly().bodies_lattice(),
        reference.assembly().bodies_lattice()
    );
    assert_eq!(actual.assembly().accepted_steps(), 2);
}
