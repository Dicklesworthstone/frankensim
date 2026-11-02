//! Prescribed film reservoirs on the material assembly's actual accepted clock.
use super::*;
use fs_couple::material_immersed::dry::{PlaneFilmSi, ResolvedImmersedPlanes};
use fs_lbm::d3q19::DryLubricatedPlanesExchange3;

fn plane() -> PlaneFilmSi {
    PlaneFilmSi {
        point_m: [0.0229, 0.0, 0.0],
        normal: [2.0, 0.0, 0.0],
        velocity_m_s: [0.0001, 0.0, 0.0],
        cutoff_gap_m: 0.003,
        roughness_gap_m: 1e-5,
        tangential_gaps_m: None,
    }
}

fn model<'a>(
    liquid: &'a ResolvedLiquidState,
    solid: &'a ResolvedBody,
) -> ResolvedImmersedPlanes<'a, ResolvedLiquidState, 1> {
    let mut state = initial()[0];
    state.velocity_m_s = [0.0; 3];
    Model::spheres(
        liquid,
        domain(),
        [SubmergedSphere::new(solid, RADIUS).unwrap()],
        [state],
        [0.0; 3],
        film(),
        &CancelGate::new(),
    )
    .unwrap()
    .with_prescribed_planes(&[plane()], &CancelGate::new())
    .unwrap()
}

fn check_plane_balance(
    model: &ResolvedImmersedPlanes<'_, ResolvedLiquidState, 1>,
    before: SphereStateSi,
    fluid_before: [f64; 6],
    origin: [f64; 3],
    report: &DryLubricatedPlanesExchange3<1>,
) {
    let assembly = model.assembly();
    let motion = &report.motion;
    let fluid_after = full_fluid_momentum(assembly, origin);
    let mass = assembly.materials()[0].inertia().mass();
    let old = body_momentum_si(before, mass, origin);
    let new = body_momentum_si(assembly.state_si(0).unwrap(), mass, origin);
    let coupling = &motion.native.coupling;
    for axis in 0..6 {
        let unit = if axis < 3 {
            assembly.impulse_ns(1.0).unwrap()
        } else {
            assembly.angular_impulse_nms(1.0).unwrap()
        };
        let plane = if axis < 3 {
            report.planes[0].impulse[axis]
        } else {
            report.planes[0].angular_impulse[axis - 3]
        };
        let boundary = if axis < 3 {
            coupling.walls.fluid_population_impulse[axis]
        } else {
            coupling.walls.fluid_population_angular_impulse[axis - 3]
                + coupling.periodic.angular_momentum_change[axis - 3]
        };
        // The external plane RECEIVES this impulse, so it is added to actual
        // fluid-exterior plus dry-body change to close the physical system.
        let balance = fluid_after[axis]
            - fluid_before[axis]
            - (motion.interior[0][1].momentum[axis] - motion.interior[0][0].momentum[axis])
            + (new[axis] - old[axis]) / unit
            + plane
            - motion.physical_external_impulses[0][axis]
            - motion.represented_fluid_forcing_impulse[axis]
            - boundary;
        close(balance, 0.0, 3e-10);
        close(motion.combined_momentum_residual[axis], 0.0, 3e-10);
        close(motion.body_momentum_residuals[0][axis], 0.0, 3e-10);
    }
}

#[test]
fn prescribed_plane_advances_once_and_receives_the_physical_film_reaction() {
    let liquid = fluid(1000.0, 0.11);
    let solid = body(800.0);
    let mut model = model(&liquid, &solid);
    let origin = [0.011, -0.013, 0.009];
    let before = model.assembly().state_si(0).unwrap();
    let fluid_before = full_fluid_momentum(model.assembly(), origin);
    let wall = model.plane_si(0).unwrap();
    assert_eq!(model.source_planes_si(), &[plane()]);
    assert_eq!(wall.normal, [1.0, 0.0, 0.0]);
    let report = model
        .step_forces_si(
            [[0.0; 6]],
            origin,
            DryBodyCouplingConfig3::default(),
            &CancelGate::new(),
        )
        .unwrap();
    check_plane_balance(&model, before, fluid_before, origin, &report);
    let assembly = model.assembly();
    let center = report.motion.native.bodies[0].midpoint_center[0] * DX;
    let gap = center - (wall.point_m[0] + 0.5 * DT * wall.velocity_m_s[0]) - RADIUS;
    let relative = report.motion.native.coupling.bodies[0].velocity_midpoint[0] * DX / DT
        - wall.velocity_m_s[0];
    let resistance = 6.0
        * core::f64::consts::PI
        * 0.11
        * RADIUS.powi(2)
        * (1.0 / gap.max(wall.roughness_gap_m) - 1.0 / wall.cutoff_gap_m);
    close(
        assembly.impulse_ns(report.planes[0].impulse[0]).unwrap(),
        DT * resistance * relative,
        5e-18,
    );
    close(
        assembly.energy_j(report.planes[0].work).unwrap(),
        assembly.impulse_ns(report.planes[0].impulse[0]).unwrap() * wall.velocity_m_s[0],
        1e-20,
    );
    assert!(report.planes[0].work < 0.0);
    assert!(assembly.state_si(0).unwrap().velocity_m_s[0] > 0.0);
    let film = &report.motion.native.coupling.lubrication;
    assert!(film.dissipated_energy > 0.0);
    close(
        film.body_work[0] + report.planes[0].work + film.dissipated_energy,
        0.0,
        2e-12,
    );
    close(
        model.plane_si(0).unwrap().point_m[0],
        wall.point_m[0] + DT * wall.velocity_m_s[0],
        1e-16,
    );
    assert_eq!(assembly.accepted_steps(), 1);
    assert_eq!(model.source_planes_si(), &[plane()]);
}

#[test]
fn empty_planes_preserve_existing_history_and_every_sphere_only_result() {
    let liquid = fluid(1000.0, 0.11);
    let solids = [body(800.0), body(1800.0)];
    let mut plain = pair(&liquid, [&solids[0], &solids[1]]);
    let mut attached = pair(&liquid, [&solids[0], &solids[1]]);
    let gate = CancelGate::new();
    let controls = DryBodyCouplingConfig3::default();
    let loads = [[2e-7, 0.0, 0.0, 0.0, 0.0, 1e-11], [0.0; 6]];
    assert_eq!(
        plain
            .step_forces_si(loads, [0.0; 3], controls, &gate)
            .unwrap(),
        attached
            .step_forces_si(loads, [0.0; 3], controls, &gate)
            .unwrap()
    );
    let mut attached = attached.with_prescribed_planes(&[], &gate).unwrap();
    assert_eq!(attached.assembly().accepted_steps(), 1);
    let report = attached
        .step_forces_si(loads, [0.0; 3], controls, &gate)
        .unwrap();
    assert_eq!(
        report.motion,
        plain
            .step_forces_si(loads, [0.0; 3], controls, &gate)
            .unwrap()
    );
    assert!(report.planes.is_empty());
    assert_eq!(fields(attached.assembly()), fields(&plain));
    assert_eq!(attached.assembly().bodies_lattice(), plain.bodies_lattice());
    assert_eq!(attached.assembly().accepted_steps(), 2);
}

#[test]
fn plane_clock_rebind_and_failed_native_probes_preserve_exact_retry() {
    let first = fluid(1000.0, 0.11);
    let second = fluid(1000.0, 0.22);
    let wrong_density = fluid(1200.0, 0.22);
    let solid = body(800.0);
    let mut actual = model(&first, &solid);
    let mut reference = model(&first, &solid);
    let gate = CancelGate::new();
    let controls = DryBodyCouplingConfig3::default();
    let loads = [[2e-7, 0.0, 0.0, 0.0, 0.0, 1e-11]];
    assert_eq!(
        actual
            .step_forces_si(loads, [0.0; 3], controls, &gate)
            .unwrap(),
        reference
            .step_forces_si(loads, [0.0; 3], controls, &gate)
            .unwrap()
    );
    let populations = fields(actual.assembly());
    let bodies = *actual.assembly().bodies_lattice();
    let planes = actual.planes_lattice().to_vec();
    actual.rebind_fluid(&second, &gate).unwrap();
    reference.rebind_fluid(&second, &gate).unwrap();
    assert_eq!(actual.planes_lattice(), planes);
    assert_eq!(fields(actual.assembly()), populations);
    assert_eq!(actual.assembly().bodies_lattice(), &bodies);
    assert!(core::ptr::eq(
        actual.assembly().scaling().material(),
        &second
    ));
    assert!(actual.rebind_fluid(&wrong_density, &gate).is_err());
    let limited = DryBodyCouplingConfig3 {
        maximum_evaluations: 1,
        ..controls
    };
    assert!(
        actual
            .step_forces_si(loads, [0.0; 3], limited, &gate)
            .is_err()
    );
    let stopped = CancelGate::new();
    stopped.request();
    assert!(
        actual
            .step_forces_si(loads, [0.0; 3], controls, &stopped)
            .is_err()
    );
    assert!(
        actual
            .step_forces_si([[f64::NAN; 6]], [0.0; 3], controls, &gate)
            .is_err()
    );
    assert_eq!(actual.planes_lattice(), planes);
    assert_eq!(fields(actual.assembly()), populations);
    assert_eq!(actual.assembly().bodies_lattice(), &bodies);
    assert_eq!(actual.assembly().accepted_steps(), 1);
    assert_eq!(
        actual
            .step_forces_si(loads, [0.0; 3], controls, &gate)
            .unwrap(),
        reference
            .step_forces_si(loads, [0.0; 3], controls, &gate)
            .unwrap()
    );
    assert_eq!(actual.planes_lattice(), reference.planes_lattice());
    assert_eq!(fields(actual.assembly()), fields(reference.assembly()));
    assert_eq!(
        actual.assembly().bodies_lattice(),
        reference.assembly().bodies_lattice()
    );
    assert_eq!(actual.assembly().accepted_steps(), 2);
}
