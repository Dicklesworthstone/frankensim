//! Material-defined physical dry assemblies on the original simultaneous owner.
use super::*;
use fs_couple::material_immersed::dry::ResolvedImmersedBodies;
use fs_couple::material_immersed::settling::{SubmergedEllipsoid, SubmergedSphere};
use fs_ga::{Se3, So3Tangent, Vec3};
use fs_lbm::d3q19::{DryBodiesExchange3, DryBodyCouplingConfig3};

type Model<'a, const N: usize> = ResolvedImmersedBodies<'a, ResolvedLiquidState, N>;

fn vector(v: [f64; 3]) -> Vec3 {
    Vec3::new(v[0], v[1], v[2])
}
fn array(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
fn fields<const N: usize>(model: &Model<'_, N>) -> Vec<[f64; Q3]> {
    (0..8)
        .flat_map(|z| {
            (0..8).flat_map(move |y| (0..8).map(move |x| model.grid().populations(x, y, z)))
        })
        .collect()
}
fn states<const N: usize>(model: &Model<'_, N>) -> [SphereStateSi; N] {
    core::array::from_fn(|i| model.state_si(i).unwrap())
}
fn full_fluid_momentum<const N: usize>(model: &Model<'_, N>, origin_m: [f64; 3]) -> [f64; 6] {
    let mut linear = Vec3::default();
    let mut angular = Vec3::default();
    let origin = vector(origin_m).scale(1.0 / DX);
    // There is no imposed fluid force in this reduced-pressure adapter, so
    // actual population moments are physical Guo moments without an F/2 term.
    for z in 0..8 {
        for y in 0..8 {
            for x in 0..8 {
                let position = Vec3::new(x as f64 + 0.5, y as f64 + 0.5, z as f64 + 0.5);
                for (q, value) in model.grid().populations(x, y, z).into_iter().enumerate() {
                    let impulse =
                        Vec3::new(f64::from(E3[q].0), f64::from(E3[q].1), f64::from(E3[q].2))
                            .scale(value);
                    linear = linear + impulse;
                    angular = angular + (position - origin).cross(impulse);
                }
            }
        }
    }
    [
        linear.x, linear.y, linear.z, angular.x, angular.y, angular.z,
    ]
}
fn body_momentum_si(state: SphereStateSi, mass: f64, origin: [f64; 3]) -> [f64; 6] {
    let linear = vector(state.velocity_m_s).scale(mass);
    let angular = state
        .orientation
        .rotate(vector(state.angular_momentum_body_kg_m2_s))
        .unwrap()
        + (vector(state.center_m) - vector(origin)).cross(linear);
    [
        linear.x, linear.y, linear.z, angular.x, angular.y, angular.z,
    ]
}
fn check_balances<const N: usize>(
    model: &Model<'_, N>,
    before: [SphereStateSi; N],
    fluid_before: [f64; 6],
    origin: [f64; 3],
    report: &DryBodiesExchange3<N>,
) {
    let fluid_after = full_fluid_momentum(model, origin);
    let coupling = &report.native.coupling;
    let mut exterior_plus_dry = [0.0; 6];
    for axis in 0..6 {
        close(fluid_before[axis], report.fluid[0].momentum[axis], 2e-10);
        close(fluid_after[axis], report.fluid[1].momentum[axis], 2e-10);
        let boundary = if axis < 3 {
            coupling.walls.fluid_population_impulse[axis]
        } else {
            coupling.walls.fluid_population_angular_impulse[axis - 3]
                + coupling.periodic.angular_momentum_change[axis - 3]
        };
        exterior_plus_dry[axis] = fluid_after[axis]
            - fluid_before[axis]
            - boundary
            - report.represented_fluid_forcing_impulse[axis];
        close(report.represented_fluid_forcing_impulse[axis], 0.0, 0.0);
    }
    for i in 0..N {
        let mass = model.materials()[i].inertia().mass();
        let old = body_momentum_si(before[i], mass, origin);
        let new = body_momentum_si(model.state_si(i).unwrap(), mass, origin);
        for axis in 0..6 {
            let unit = if axis < 3 {
                model.impulse_ns(1.0).unwrap()
            } else {
                model.angular_impulse_nms(1.0).unwrap()
            };
            let observed = (new[axis] - old[axis]) / unit;
            close(observed, report.body_momentum_changes[i][axis], 2e-10);
            let film = if axis < 3 {
                coupling.lubrication.body_impulses[i][axis]
            } else {
                coupling.lubrication.body_angular_impulses[i][axis - 3]
            };
            close(
                observed
                    - report.physical_external_impulses[i][axis]
                    - report.corrected_hydrodynamic_impulses[i][axis]
                    - film,
                0.0,
                3e-10,
            );
            close(report.body_momentum_residuals[i][axis], 0.0, 3e-10);
            // Subtract each measured fictitious interior ONCE, retaining the
            // original dry body; summing all fluid plus dry mass fails here.
            exterior_plus_dry[axis] += observed
                - (report.interior[i][1].momentum[axis] - report.interior[i][0].momentum[axis])
                - report.physical_external_impulses[i][axis];
        }
    }
    for axis in 0..6 {
        close(exterior_plus_dry[axis], 0.0, 3e-10 * N as f64);
        close(
            report.combined_momentum_residual[axis],
            0.0,
            3e-10 * N as f64,
        );
    }
}
fn pair<'a>(fluid: &'a ResolvedLiquidState, solids: [&'a ResolvedBody; 2]) -> Model<'a, 2> {
    Model::spheres(
        fluid,
        domain(),
        solids.map(|s| SubmergedSphere::new(s, RADIUS).unwrap()),
        initial(),
        [0.0; 3],
        film(),
        &CancelGate::new(),
    )
    .unwrap()
}
fn check_normal_film(model: &Model<'_, 2>, report: &DryBodiesExchange3<2>, viscosity: f64) {
    let separation = (vector(report.native.bodies[1].midpoint_center)
        - vector(report.native.bodies[0].midpoint_center))
    .scale(DX);
    let normal = separation.scale(1.0 / separation.norm());
    let gap = separation.norm() - 2.0 * RADIUS;
    let relative = (vector(report.native.coupling.bodies[0].velocity_midpoint)
        - vector(report.native.coupling.bodies[1].velocity_midpoint))
    .scale(DX / DT)
    .dot(normal);
    let resistance = 6.0
        * core::f64::consts::PI
        * viscosity
        * (RADIUS / 2.0).powi(2)
        * (1.0 / gap.max(film().roughness_gap_m) - 1.0 / film().cutoff_gap_m);
    let film = &report.native.coupling.lubrication;
    assert_eq!(film.active_pairs, 1);
    assert!(film.dissipated_energy > 0.0);
    let expected = array(normal.scale(-resistance * relative * DT));
    for axis in 0..3 {
        close(
            model.impulse_ns(film.body_impulses[0][axis]).unwrap(),
            expected[axis],
            5e-18,
        );
        close(
            film.body_impulses[0][axis] + film.body_impulses[1][axis],
            0.0,
            2e-12,
        );
        close(
            film.body_angular_impulses[0][axis] + film.body_angular_impulses[1][axis],
            0.0,
            2e-12,
        );
    }
    close(
        film.body_work.iter().sum::<f64>() + film.dissipated_energy,
        0.0,
        2e-12,
    );
}

#[test]
fn mixed_material_densities_keep_original_dry_inertia_and_reduce_gravity_once() {
    let liquid = fluid(1000.0, 0.11);
    let solids = [body(500.0), body(1000.0), body(2000.0)];
    let descriptors = core::array::from_fn(|i| SubmergedSphere::new(&solids[i], RADIUS).unwrap());
    let initial = [0.02, 0.04, 0.06].map(|x| SphereStateSi {
        center_m: [x, 0.04, 0.04],
        orientation: So3::identity(),
        velocity_m_s: [0.0; 3],
        angular_momentum_body_kg_m2_s: [0.0; 3],
    });
    let gravity = [0.0, 0.0, -1e-4];
    let mut model = Model::spheres(
        &liquid,
        domain(),
        descriptors,
        initial,
        gravity,
        film(),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(model.gravity_m_s2(), gravity);
    assert!(core::ptr::eq(model.scaling().material(), &liquid));
    for i in 0..3 {
        assert!(core::ptr::eq(model.materials()[i], &solids[i]));
        let original = solids[i].inertia();
        close(
            model.mass_kg(model.bodies_lattice()[i].mass).unwrap(),
            original.mass(),
            1e-16,
        );
        for entry in 0..9 {
            close(
                model
                    .inertia_kg_m2(model.bodies_lattice()[i].inertia.matrix().m[entry])
                    .unwrap(),
                original.inertia_com().m[entry],
                1e-20,
            );
        }
        let volume = solids[i].regions()[0].geometry().volume_m3();
        assert_eq!(model.volumes_m3()[i], volume);
        close(model.displaced_masses_kg()[i], 1000.0 * volume, 1e-16);
        close(
            model.reduced_gravity_forces_n()[i][2],
            (original.mass() - 1000.0 * volume) * gravity[2],
            1e-20,
        );
    }
    assert_eq!(model.reduced_gravity_forces_n()[1], [0.0; 3]);
    let origin = [0.011, -0.013, 0.009];
    let before = states(&model);
    let fluid_before = full_fluid_momentum(&model, origin);
    let report = model
        .step_forces_si(
            [[0.0; 6]; 3],
            origin,
            DryBodyCouplingConfig3::default(),
            &CancelGate::new(),
        )
        .unwrap();
    check_balances(&model, before, fluid_before, origin, &report);
    for i in 0..3 {
        close(
            model
                .impulse_ns(report.physical_external_impulses[i][2])
                .unwrap(),
            model.reduced_gravity_forces_n()[i][2] * DT,
            1e-18,
        );
        close(
            model.mass_kg(model.bodies_lattice()[i].mass).unwrap(),
            solids[i].inertia().mass(),
            1e-16,
        );
    }
    assert!(model.state_si(0).unwrap().velocity_m_s[2] > 0.0);
    assert!(model.state_si(2).unwrap().velocity_m_s[2] < 0.0);
    assert_eq!(report.native.coupling.lubrication.active_pairs, 0);
    assert_eq!(model.accepted_steps(), 1);
}

#[test]
fn active_material_film_and_authored_world_wrenches_share_the_native_balance() {
    let liquid = fluid(1000.0, 0.11);
    let solids = [body(800.0), body(1800.0)];
    let mut model = pair(&liquid, [&solids[0], &solids[1]]);
    let origin = [0.011, -0.013, 0.009];
    let loads = [
        [2e-7, 1e-7, 0.0, 0.0, 0.0, 2e-11],
        [0.0, 0.0, 0.0, -1e-11, 0.0, 0.0],
    ];
    let before = states(&model);
    let fluid_before = full_fluid_momentum(&model, origin);
    let report = model
        .step_forces_si(
            loads,
            origin,
            DryBodyCouplingConfig3::default(),
            &CancelGate::new(),
        )
        .unwrap();
    check_balances(&model, before, fluid_before, origin, &report);
    check_normal_film(&model, &report, 0.11);
    for i in 0..2 {
        let impulse = vector([loads[i][0], loads[i][1], loads[i][2]]).scale(DT);
        let arm = vector(report.native.bodies[i].midpoint_center).scale(DX) - vector(origin);
        let angular =
            arm.cross(impulse) + vector([loads[i][3], loads[i][4], loads[i][5]]).scale(DT);
        for axis in 0..3 {
            close(
                model
                    .impulse_ns(report.physical_external_impulses[i][axis])
                    .unwrap(),
                array(impulse)[axis],
                1e-18,
            );
            close(
                model
                    .angular_impulse_nms(report.physical_external_impulses[i][axis + 3])
                    .unwrap(),
                array(angular)[axis],
                1e-20,
            );
        }
    }
}

#[test]
fn density_preserving_rebind_keeps_history_and_failed_assembly_steps_retry_exactly() {
    let first = fluid(1000.0, 0.11);
    let second = fluid(1000.0, 0.22);
    let changed_density = fluid(1200.0, 0.22);
    let solids = [body(800.0), body(1800.0)];
    let mut model = pair(&first, [&solids[0], &solids[1]]);
    let mut reference = pair(&first, [&solids[0], &solids[1]]);
    let gate = CancelGate::new();
    let loads = [[2e-7, 0.0, 0.0, 0.0, 0.0, 1e-11], [0.0; 6]];
    let controls = DryBodyCouplingConfig3::default();
    assert_eq!(
        model
            .step_forces_si(loads, [0.0; 3], controls, &gate)
            .unwrap(),
        reference
            .step_forces_si(loads, [0.0; 3], controls, &gate)
            .unwrap()
    );
    let populations = fields(&model);
    let bodies = *model.bodies_lattice();
    let reduced_weights = *model.reduced_gravity_forces_n();
    let mass_unit = model.mass_kg(1.0).unwrap();
    model.rebind_fluid(&second, &gate).unwrap();
    reference.rebind_fluid(&second, &gate).unwrap();
    assert!(core::ptr::eq(model.scaling().material(), &second));
    assert_eq!(fields(&model), populations);
    assert_eq!(model.bodies_lattice(), &bodies);
    assert_eq!(*model.reduced_gravity_forces_n(), reduced_weights);
    assert_eq!(model.mass_kg(1.0).unwrap(), mass_unit);
    assert_eq!(model.accepted_steps(), 1);
    close(
        model.grid().viscosity(),
        0.22 / 1000.0 * DT / DX.powi(2),
        1e-16,
    );
    let stopped = CancelGate::new();
    stopped.request();
    assert!(matches!(
        model.rebind_fluid(&changed_density, &gate),
        Err(MaterialImmersedError::Invalid(_))
    ));
    assert!(matches!(
        model.rebind_fluid(&first, &stopped),
        Err(MaterialImmersedError::Native(
            PartialSaturationError3::Cancelled
        ))
    ));
    assert!(
        model
            .step_forces_si(loads, [0.0; 3], controls, &stopped)
            .is_err()
    );
    assert!(
        model
            .step_forces_si([[f64::NAN; 6]; 2], [0.0; 3], controls, &gate)
            .is_err()
    );
    let limited = DryBodyCouplingConfig3 {
        maximum_evaluations: 1,
        ..controls
    };
    assert!(
        matches!(model.step_forces_si(loads, [0.0; 3], limited, &gate), Err(MaterialImmersedError::Native(PartialSaturationError3::Invalid { reason })) if reason.contains("budget"))
    );
    assert_eq!(fields(&model), populations);
    assert_eq!(model.bodies_lattice(), &bodies);
    assert_eq!(model.accepted_steps(), 1);
    assert!(core::ptr::eq(model.scaling().material(), &second));
    let accepted = model
        .step_forces_si(loads, [0.0; 3], controls, &gate)
        .unwrap();
    assert_eq!(
        accepted,
        reference
            .step_forces_si(loads, [0.0; 3], controls, &gate)
            .unwrap()
    );
    assert_eq!(fields(&model), fields(&reference));
    assert_eq!(model.bodies_lattice(), reference.bodies_lattice());
    assert_eq!(model.accepted_steps(), 2);
    check_normal_film(&model, &accepted, 0.22);
    let shape = SubmergedSphere::new(&solids[0], RADIUS).unwrap();
    assert!(Model::spheres(&first, domain(), [], [], [0.0; 3], film(), &gate).is_err());
    assert!(
        Model::spheres(
            &first,
            domain(),
            [shape; 9],
            [initial()[0]; 9],
            [0.0; 3],
            film(),
            &gate
        )
        .is_err()
    );
}

#[test]
fn ellipsoid_assembly_retains_the_rotated_original_tensor_and_physical_spin() {
    let liquid = fluid(1000.0, 0.11);
    let axes = [0.006, 0.0045, 0.003];
    let body_from_shape = So3::exp(So3Tangent::new(Vec3::new(0.35, -0.2, 0.25))).unwrap();
    let geometry = VolumeMoments::ellipsoid(vector(axes))
        .unwrap()
        .transformed(Se3::from_parts(body_from_shape, Vec3::default()).unwrap())
        .unwrap();
    let solids = [750.0, 1500.0].map(|density| {
        resolve_body(
            &[RegionBinding {
                key: "ellipsoid",
                geometry,
                material: &card(
                    "dry ellipsoid",
                    &[("density", Density::DIMS, density)],
                    "solid",
                ),
                point: &point(),
                selection: MaterialPropertySelection::SingleClaimOnly,
            }],
            1,
            || false,
        )
        .unwrap()
    });
    let sources = core::array::from_fn(|i| {
        SubmergedEllipsoid::new(&solids[i], axes, body_from_shape).unwrap()
    });
    let initial = [0.0307, 0.0503].map(|x| SphereStateSi {
        center_m: [x, 0.0391, 0.0413],
        orientation: So3::exp(So3Tangent::new(Vec3::new(0.1, 0.2, -0.15))).unwrap(),
        velocity_m_s: [1e-5, -2e-5, 1e-5],
        angular_momentum_body_kg_m2_s: [1e-12, -2e-12, 3e-12],
    });
    let mut model = Model::ellipsoids(
        &liquid,
        domain(),
        sources,
        initial,
        [0.0; 3],
        &CancelGate::new(),
    )
    .unwrap();
    for i in 0..2 {
        assert!(solids[i].inertia().inertia_com().m[1].abs() > 1e-11);
        for entry in 0..9 {
            close(
                model
                    .inertia_kg_m2(model.bodies_lattice()[i].inertia.matrix().m[entry])
                    .unwrap(),
                solids[i].inertia().inertia_com().m[entry],
                1e-20,
            );
        }
        for axis in 0..3 {
            close(
                model.state_si(i).unwrap().angular_momentum_body_kg_m2_s[axis],
                initial[i].angular_momentum_body_kg_m2_s[axis],
                1e-26,
            );
        }
    }
    let origin = [0.011, -0.013, 0.009];
    let before = states(&model);
    let fluid_before = full_fluid_momentum(&model, origin);
    let report = model
        .step_forces_si(
            [[1e-8, 0.0, 0.0, 0.0, 1e-12, -2e-12]; 2],
            origin,
            DryBodyCouplingConfig3::default(),
            &CancelGate::new(),
        )
        .unwrap();
    check_balances(&model, before, fluid_before, origin, &report);
    assert_eq!(report.native.coupling.lubrication.active_pairs, 0);
    assert!(
        model
            .state_si(0)
            .unwrap()
            .angular_momentum_body_kg_m2_s
            .iter()
            .any(|value| *value != 0.0)
    );
    assert_eq!(model.accepted_steps(), 1);
}
