//! Physical units and atomic admission of two walls on one liquid clock.

use super::*;
use fs_lbm::d3q19::freesurface3::PlanarRigidWallMotion3;
use fs_lbm::material::free_surface::MaterialPlanarRigidWallMotion3;

fn tank_domain() -> MaterialFreeSurfaceDomain3<'static> {
    MaterialFreeSurfaceDomain3 {
        cells: [8; 3],
        max_cells: 512,
        ..domain()
    }
}

fn motions(axis: usize) -> [MaterialPlanarRigidWallMotion3; 2] {
    let faces = [
        [PlanarWall3::XMin, PlanarWall3::XMax],
        [PlanarWall3::YMin, PlanarWall3::YMax],
        [PlanarWall3::ZMin, PlanarWall3::ZMax],
    ][axis];
    core::array::from_fn(|end| {
        let mut origin_m = [0.004; 3];
        origin_m[axis] = 0.008 * end as f64;
        let mut velocity_at_origin_m_s = [0.0; 3];
        velocity_at_origin_m_s[(axis + 1) % 3] = if end == 0 { 0.025 } else { -0.02 };
        velocity_at_origin_m_s[(axis + 2) % 3] = if end == 0 { -0.01 } else { 0.012 };
        MaterialPlanarRigidWallMotion3 {
            wall: faces[end],
            origin_m,
            velocity_at_origin_m_s,
            angular_speed_rad_s: if end == 0 { 1.5 } else { -0.8 },
        }
    })
}

fn native_motion(
    motion: MaterialPlanarRigidWallMotion3,
    setup: &MaterialFreeSurfaceScaling<'_>,
) -> PlanarRigidWallMotion3 {
    PlanarRigidWallMotion3 {
        wall: motion.wall,
        origin: motion.origin_m.map(|v| v / setup.cell_size_m()),
        velocity_at_origin: motion
            .velocity_at_origin_m_s
            .map(|v| v / setup.velocity_si(1.0).unwrap()),
        angular_speed: motion.angular_speed_rad_s * setup.timestep_s(),
    }
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[test]
fn g1_pair_si_is_one_native_tick_with_independent_impedance_and_work_on_every_axis() {
    let source = liquid(0, 300.0, None);
    let setup = plan(&source).with_pressure_reference(100_000.0).unwrap();
    let gate = CancelGate::new_clock_free();
    let dx = setup.cell_size_m();
    let cell_mass = source.density_kg_m3() * dx.powi(3);
    let impulse_unit = setup.mass_unit_kg() * setup.velocity_si(1.0).unwrap();
    let work_unit = impulse_unit * setup.velocity_si(1.0).unwrap();
    let mut spin_response = 0.0;
    for j in 0..8 {
        for i in 0..8 {
            let ra = (i as f64 + 0.5 - 4.0) * dx;
            let rb = (j as f64 + 0.5 - 4.0) * dx;
            spin_response += cell_mass * (ra * ra + rb * rb) / 3.0;
        }
    }
    for axis in 0..3 {
        let fluid = setup.free_surface3(tank_domain(), |_, _, _| true).unwrap();
        let before = fluid.lattice().grid().state_bits();
        let motion = motions(axis);
        let (candidate, report) = fluid
            .planar_rigid_wall_pair_candidate_si(motion, 1e-12, &gate)
            .unwrap();
        let (native, raw) = fluid
            .lattice()
            .planar_rigid_wall_pair_candidate(
                motion.map(|m| native_motion(m, &setup)),
                1e-12,
                &gate,
            )
            .unwrap();
        assert_eq!(
            candidate.lattice().grid().state_bits(),
            native.grid().state_bits()
        );
        assert_eq!(fluid.lattice().grid().state_bits(), before);
        assert_eq!(candidate.scaling().material(), &source);
        close(report.mass_before_kg, 512.0 * cell_mass, 2e-18);
        close(report.mass_after_kg, report.mass_before_kg, 2e-18);
        let mut independent_work = 0.0;
        for wall in report.walls {
            let a = (axis + 1) % 3;
            let b = (axis + 2) % 3;
            for d in [a, b] {
                close(
                    wall.fluid_impulse_kg_m_s[d],
                    cell_mass * 64.0 * wall.motion.velocity_at_origin_m_s[d] / 3.0,
                    2e-18,
                );
            }
            close(wall.fluid_response[0][0], cell_mass * 64.0 / 3.0, 2e-18);
            close(wall.fluid_response[1][1], cell_mass * 64.0 / 3.0, 2e-18);
            close(wall.fluid_response[2][2], spin_response, 2e-21);
            close(
                wall.fluid_angular_impulse_kg_m2_s[axis],
                spin_response * wall.motion.angular_speed_rad_s,
                2e-21,
            );
            let work = wall.motion.velocity_at_origin_m_s[a] * wall.fluid_impulse_kg_m_s[a]
                + wall.motion.velocity_at_origin_m_s[b] * wall.fluid_impulse_kg_m_s[b]
                + wall.motion.angular_speed_rad_s * wall.fluid_angular_impulse_kg_m2_s[axis];
            close(wall.fluid_work_j, work, 2e-20);
            independent_work += work;
            assert_eq!(
                wall.wall_impulse_kg_m_s,
                wall.fluid_impulse_kg_m_s.map(|v| -v)
            );
            assert_eq!(
                wall.wall_angular_impulse_kg_m2_s,
                wall.fluid_angular_impulse_kg_m2_s.map(|v| -v)
            );
            assert_eq!(wall.wall_work_j, -wall.fluid_work_j);
            assert_eq!(wall.mass_before_kg, report.mass_before_kg);
            assert_eq!(wall.mass_after_kg, report.mass_after_kg);
        }
        close(report.fluid_work_j, independent_work, 2e-20);
        close(report.fluid_work_j, raw.fluid_work * work_unit, 2e-20);
        for d in 0..3 {
            close(
                report.fluid_impulse_kg_m_s[d],
                raw.fluid_impulse[d] * impulse_unit,
                2e-18,
            );
        }
        assert_eq!(
            report.wall_impulse_kg_m_s,
            report.fluid_impulse_kg_m_s.map(|v| -v)
        );
        assert_eq!(report.wall_work_j, -report.fluid_work_j);
    }
}

#[test]
fn g2_pair_moments_use_a_common_origin_and_pressure_datum_changes_only_support_loads() {
    let source = liquid(0, 300.0, None);
    let setup = plan(&source).with_pressure_reference(100_000.0).unwrap();
    let second_setup = plan(&source).with_pressure_reference(101_000.0).unwrap();
    let first = setup.free_surface3(tank_domain(), |_, _, _| true).unwrap();
    let second = second_setup
        .free_surface3(tank_domain(), |_, _, _| true)
        .unwrap();
    let gate = CancelGate::new_clock_free();
    let mut motion = motions(2);
    motion[0].origin_m = [-0.002, 0.006, 0.001];
    motion[1].origin_m = [0.005, -0.001, 0.010];
    let (candidate, report) = first
        .planar_rigid_wall_pair_candidate_si(motion, 1e-12, &gate)
        .unwrap();
    let (shifted, pressure) = second
        .planar_rigid_wall_pair_candidate_si(motion, 1e-12, &gate)
        .unwrap();
    assert_eq!(
        candidate.lattice().grid().state_bits(),
        shifted.lattice().grid().state_bits()
    );
    assert_eq!(report.fluid_impulse_kg_m_s, pressure.fluid_impulse_kg_m_s);
    assert_eq!(report.fluid_work_j, pressure.fluid_work_j);
    for end in 0..2 {
        let old = report.walls[end];
        let new = pressure.walls[end];
        let mut force_change = [0.0; 3];
        force_change[2] = (if end == 0 { -1.0 } else { 1.0 }) * 1000.0 * old.wetted_area_m2;
        let torque_change = cross(
            core::array::from_fn(|d| old.wetted_centroid_m[d] - old.motion.origin_m[d]),
            force_change,
        );
        for d in 0..3 {
            close(
                new.absolute_wall_force_n.unwrap()[d] - old.absolute_wall_force_n.unwrap()[d],
                force_change[d],
                2e-15,
            );
            close(
                new.absolute_wall_torque_nm.unwrap()[d] - old.absolute_wall_torque_nm.unwrap()[d],
                torque_change[d],
                2e-17,
            );
        }
        assert_eq!(
            old.pressure_perturbation_wall_force_n,
            new.pressure_perturbation_wall_force_n
        );
        assert_eq!(
            old.pressure_perturbation_wall_torque_nm,
            new.pressure_perturbation_wall_torque_nm
        );
    }
    let origin = [0.002, -0.003, 0.012];
    let other_origin = [-0.004, 0.005, -0.001];
    let angular = report.fluid_angular_impulse_at_si(origin).unwrap();
    let shifted_angular = report.fluid_angular_impulse_at_si(other_origin).unwrap();
    let orbital_change = cross(
        core::array::from_fn(|d| origin[d] - other_origin[d]),
        report.fluid_impulse_kg_m_s,
    );
    let (_, raw) = first
        .lattice()
        .planar_rigid_wall_pair_candidate(motion.map(|m| native_motion(m, &setup)), 1e-12, &gate)
        .unwrap();
    let raw_angular = raw
        .fluid_angular_impulse_at(origin.map(|v| v / setup.cell_size_m()))
        .unwrap();
    let angular_unit = setup.mass_unit_kg() * setup.velocity_si(1.0).unwrap() * setup.cell_size_m();
    for d in 0..3 {
        close(shifted_angular[d], angular[d] + orbital_change[d], 2e-20);
        close(angular[d], raw_angular[d] * angular_unit, 2e-20);
    }
    assert_eq!(
        report.wall_angular_impulse_at_si(origin).unwrap(),
        angular.map(|v| -v)
    );
    assert!(report.fluid_angular_impulse_at_si([f64::NAN; 3]).is_err());
}

#[test]
fn g3_caller_order_preserves_the_same_liquid_and_swaps_only_wall_accounts() {
    let source = liquid(0, 300.0, None);
    let fluid = plan(&source)
        .free_surface3(tank_domain(), |_, _, _| true)
        .unwrap();
    let gate = CancelGate::new_clock_free();
    let motion = motions(1);
    let (first, a) = fluid
        .planar_rigid_wall_pair_candidate_si(motion, 1e-12, &gate)
        .unwrap();
    let (second, b) = fluid
        .planar_rigid_wall_pair_candidate_si([motion[1], motion[0]], 1e-12, &gate)
        .unwrap();
    assert_eq!(
        first.lattice().grid().state_bits(),
        second.lattice().grid().state_bits()
    );
    assert_eq!(a.walls, [b.walls[1], b.walls[0]]);
    assert_eq!(a.mass_after_kg, b.mass_after_kg);
    assert_eq!(a.fluid_impulse_kg_m_s, b.fluid_impulse_kg_m_s);
    assert_eq!(a.fluid_work_j, b.fluid_work_j);
}

#[test]
fn g4_stationary_other_wall_recovers_the_original_single_wall_trajectory() {
    let source = liquid(0, 300.0, None);
    let fluid = plan(&source)
        .free_surface3(tank_domain(), |_, _, _| true)
        .unwrap();
    let gate = CancelGate::new_clock_free();
    let mut motion = motions(0);
    motion[1].velocity_at_origin_m_s = [0.0; 3];
    motion[1].angular_speed_rad_s = 0.0;
    let (pair, pair_report) = fluid
        .planar_rigid_wall_pair_candidate_si(motion, 1e-12, &gate)
        .unwrap();
    let (single, single_report) = fluid
        .planar_rigid_wall_motion_candidate_si(motion[0], 1e-12, &gate)
        .unwrap();
    assert_eq!(
        pair.lattice().grid().state_bits(),
        single.lattice().grid().state_bits()
    );
    assert_eq!(pair_report.walls[0], single_report);
    assert_eq!(pair_report.walls[1].fluid_work_j, 0.0);
}

#[test]
fn g5_bad_second_wall_cancellation_and_late_scenario_speed_refuse_atomically() {
    let source = liquid(0, 300.0, None);
    let mut fluid = plan(&source)
        .free_surface3(tank_domain(), |_, _, _| true)
        .unwrap();
    let before = fluid.lattice().grid().state_bits();
    let mass = fluid.ledger_mass_kg().unwrap();
    let counters = fluid.lattice().conversions();
    let motion = motions(2);
    let mut invalid = [motion; 4];
    invalid[0][1].wall = PlanarWall3::ZMin;
    invalid[1][1].velocity_at_origin_m_s[2] = 0.001;
    invalid[2][1].angular_speed_rad_s = 100.0;
    invalid[3][1].origin_m[0] = f64::NAN;
    for pair in invalid {
        assert!(
            fluid
                .step_with_planar_rigid_wall_pair_si(pair, 1e-12, &CancelGate::new_clock_free())
                .is_err()
        );
        assert_eq!(fluid.lattice().grid().state_bits(), before);
        assert_eq!(fluid.ledger_mass_kg().unwrap(), mass);
        assert_eq!(fluid.lattice().conversions(), counters);
    }
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert_eq!(
        fluid.step_with_planar_rigid_wall_pair_si(motion, 1e-12, &cancelled),
        Err(MaterialLbmError::FreeSurface(FreeSurfaceError3::Cancelled))
    );
    assert_eq!(fluid.lattice().grid().state_bits(), before);
    let (expected, receipt) = fluid
        .planar_rigid_wall_pair_candidate_si(motion, 1e-12, &CancelGate::new_clock_free())
        .unwrap();
    assert_eq!(
        fluid
            .step_with_planar_rigid_wall_pair_si(motion, 1e-12, &CancelGate::new_clock_free())
            .unwrap(),
        receipt
    );
    assert_eq!(
        fluid.lattice().grid().state_bits(),
        expected.lattice().grid().state_bits()
    );

    // The constructor starts at physical rest. This held acceleration produces
    // about 0.03 m/s in one tick, decisively above the 0.02 m/s source scenario
    // while remaining well within the native Mach and mass domains.
    let small = plan_free_surface_from_material(&source, 1e-3, 2e-4, 0.02).unwrap();
    let mut forced = small
        .free_surface3(
            MaterialFreeSurfaceDomain3 {
                gravity_m_s2: [0.0, 0.0, -150.0],
                ..tank_domain()
            },
            |_, _, _| true,
        )
        .unwrap();
    let original = forced.lattice().grid().state_bits();
    let gentle = motion.map(|m| MaterialPlanarRigidWallMotion3 {
        velocity_at_origin_m_s: [0.005, 0.0, 0.0],
        angular_speed_rad_s: 0.0,
        ..m
    });
    assert!(matches!(
        forced.step_with_planar_rigid_wall_pair_si(gentle, 1e-12, &CancelGate::new_clock_free()),
        Err(MaterialLbmError::Invalid(reason)) if reason.contains("liquid exceeds")
    ));
    assert_eq!(forced.lattice().grid().state_bits(), original);
}
