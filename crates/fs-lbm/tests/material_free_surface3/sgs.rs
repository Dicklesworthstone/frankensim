//! G1/G3 native wall/SGS composition and physical units; G4 atomic refusals.

use super::*;
use fs_lbm::d3q19::freesurface3::{
    FreeSurfaceKinematicsWorkspace3, PlanarWallFieldControls3, PlanarWallVelocityField3,
};
use fs_lbm::d3q19::freesurface3::sgs::{FreeSurfaceSgsError3, FreeSurfaceSgsWorkspace3};
use fs_lbm::material::free_surface::{
    MaterialFreeSurface3, MaterialFreeSurfaceSgsConfig3, MaterialFreeSurfaceSgsError3,
    MaterialPlanarWallVelocityField3,
};
use fs_lbm::sgs3::{SgsConfig, SgsModel, SgsUse};

fn controls() -> PlanarWallFieldControls3 {
    PlanarWallFieldControls3 {
        maximum_domain_cells: 16 * 16 * 12,
        maximum_wall_faces: 256,
        maximum_relative_mass_defect: 1e-12,
        maximum_relative_face_mass_defect: 1e-12,
    }
}

fn velocities() -> Vec<[f64; 2]> {
    (0..256).map(|i| [
        0.01 + 0.001 * (i % 16) as f64 - 0.0002 * (i / 16) as f64,
        -0.004 + 0.0003 * (i / 16) as f64,
    ]).collect()
}

fn motion(values: &[[f64; 2]]) -> MaterialPlanarWallVelocityField3<'_> {
    MaterialPlanarWallVelocityField3 {
        wall: PlanarWall3::ZMin,
        moment_origin_m: [0.0; 3],
        tangential_velocity_m_s: values,
    }
}

fn config(model: SgsModel, dx: f64) -> MaterialFreeSurfaceSgsConfig3 {
    MaterialFreeSurfaceSgsConfig3 {
        model,
        filter_width_m: 1.7 * dx,
        maximum_eddy_to_molecular: 100.0,
        usage: SgsUse::Exploration,
    }
}

fn native_config(model: SgsModel) -> SgsConfig {
    SgsConfig {
        model,
        filter_width_lu: 1.7,
        maximum_eddy_to_molecular: 100.0,
        usage: SgsUse::Exploration,
    }
}

fn cache(flow: &FreeSurface3) -> FreeSurfaceKinematicsWorkspace3 {
    FreeSurfaceKinematicsWorkspace3::new(flow.grid().active_tiles(), usize::MAX).unwrap()
}

fn near(actual: f64, expected: f64) {
    let tolerance = 2e-12 * expected.abs().max(1e-30);
    assert!((actual - expected).abs() <= tolerance, "{actual} != {expected}");
}

fn primed(source: &ResolvedLiquidSurfaceState) -> MaterialFreeSurface3<'_> {
    let mut flow = plan(source).with_pressure_reference(100_000.0).unwrap()
        .free_surface3(domain(), |_, _, z| z < 6).unwrap();
    flow.step_with_planar_wall_velocity_field_si(
        motion(&velocities()), controls(), &CancelGate::new_clock_free(),
    ).unwrap();
    flow
}

fn same(a: &FreeSurface3, b: &FreeSurface3) {
    assert_eq!(a.grid().state_bits(), b.grid().state_bits());
    assert_eq!(a.ledger_mass().to_bits(), b.ledger_mass().to_bits());
    assert_eq!(a.conversions(), b.conversions());
}

#[test]
fn g3_disabled_wall_sgs_preserves_bgk_trt_populations_and_complete_receipts() {
    let source = liquid(0, 300.0, None);
    let original = primed(&source);
    let values = velocities();
    let physical_motion = motion(&values);
    let gate = CancelGate::new_clock_free();
    let velocity_unit = original.scaling().velocity_si(1.0).unwrap();
    let native_values: Vec<_> = values.iter().map(|v| v.map(|u| u / velocity_unit)).collect();
    let native_motion = PlanarWallVelocityField3 {
        wall: physical_motion.wall,
        moment_origin: [0.0; 3],
        tangential_velocity: &native_values,
    };
    for trt in [false, true] {
        let mut input = original.lattice().clone();
        if trt { input.set_trt_magic(0.25).unwrap(); }
        let mut expected = input.clone();
        let wall = expected.step_with_planar_wall_velocity_field(
            native_motion, controls(), &gate,
        ).unwrap();
        for model in [SgsModel::Wale { cw: 0.0 }, SgsModel::Vreman { cv: 0.0 }] {
            let mut actual = input.clone();
            let mut kinematics = cache(&input);
            let mut sgs = FreeSurfaceSgsWorkspace3::default();
            let (exchange, report) = actual.step_with_planar_wall_velocity_field_and_gradient_sgs(
                native_motion, controls(), &gate, native_config(model), &mut kinematics, &mut sgs,
            ).unwrap();
            same(&actual, &expected);
            assert_eq!(exchange, wall);
            assert_eq!(report.sgs_dissipation_lu, 0.0);
            assert_eq!(report.maximum_eddy_viscosity_lu, 0.0);
        }
    }
    let mut expected = original.clone();
    let wall = expected.step_with_planar_wall_velocity_field_si(
        physical_motion, controls(), &gate,
    ).unwrap();
    for model in [SgsModel::Wale { cw: 0.0 }, SgsModel::Vreman { cv: 0.0 }] {
        let mut actual = original.clone();
        let mut kinematics = cache(original.lattice());
        let mut sgs = FreeSurfaceSgsWorkspace3::default();
        let (exchange, report) = actual.step_with_planar_wall_velocity_field_and_gradient_sgs_si(
            physical_motion, controls(), &gate, config(model, original.scaling().cell_size_m()),
            &mut kinematics, &mut sgs,
        ).unwrap();
        same(actual.lattice(), expected.lattice());
        assert_eq!(exchange, wall);
        assert_eq!(report.sgs_dissipation_w, 0.0);
        assert_eq!(actual.scaling().material(), &source);
    }
}

#[test]
fn g1_enabled_wall_sgs_changes_actual_flow_and_keeps_native_wrench_and_si_power() {
    let source = liquid(0, 300.0, None);
    let original = primed(&source);
    let before = original.lattice().grid().state_bits();
    let values = velocities();
    let physical_motion = motion(&values);
    let gate = CancelGate::new_clock_free();
    let dx = original.scaling().cell_size_m();
    let dt = original.scaling().timestep_s();
    let velocity_unit = original.scaling().velocity_si(1.0).unwrap();
    let native_values: Vec<_> = values.iter().map(|v| v.map(|u| u / velocity_unit)).collect();
    let native_motion = PlanarWallVelocityField3 {
        wall: physical_motion.wall, moment_origin: [0.0; 3], tangential_velocity: &native_values,
    };
    let (molecular, _) = original.planar_wall_velocity_field_candidate_si(
        physical_motion, controls(), &gate,
    ).unwrap();
    for model in [SgsModel::Wale { cw: 0.5 }, SgsModel::Vreman { cv: 0.1 }] {
        let mut kinematics = cache(original.lattice());
        let mut sgs = FreeSurfaceSgsWorkspace3::default();
        let (native, raw_wall, raw_sgs) =
            original.lattice().planar_wall_velocity_field_gradient_sgs_candidate(
                native_motion, controls(), &gate, native_config(model), &mut kinematics, &mut sgs,
            ).unwrap();
        let (physical, wall, report) = original.planar_wall_velocity_field_gradient_sgs_candidate_si(
            physical_motion, controls(), &gate, config(model, dx), &mut kinematics, &mut sgs,
        ).unwrap();
        same(physical.lattice(), &native);
        assert_eq!(original.lattice().grid().state_bits(), before);
        assert_eq!(physical.scaling().material(), &source);
        assert_eq!(report.wet_cells, raw_sgs.wet_cells);
        near(report.liquid_volume_m3, raw_sgs.liquid_volume_lu * dx.powi(3));
        near(report.maximum_eddy_viscosity_m2_s, raw_sgs.maximum_eddy_viscosity_lu * dx.powi(2) / dt);
        let power_unit = source.density_kg_m3() * dx.powi(5) / dt.powi(3);
        near(report.resolved_dissipation_w, raw_sgs.resolved_dissipation_lu * power_unit);
        near(report.sgs_dissipation_w, raw_sgs.sgs_dissipation_lu * power_unit);
        assert!(report.sgs_dissipation_w > 0.0);
        near(report.maximum_eddy_to_molecular, raw_sgs.maximum_eddy_to_molecular);
        let impulse_unit = source.density_kg_m3() * dx.powi(4) / dt;
        for axis in 0..3 {
            near(wall.exchange.fluid_impulse_kg_m_s[axis], raw_wall.fluid_impulse[axis] * impulse_unit);
            near(wall.exchange.fluid_angular_impulse_kg_m2_s[axis], raw_wall.fluid_angular_impulse[axis] * impulse_unit * dx);
            assert_eq!(wall.exchange.wall_impulse_kg_m_s[axis], -wall.exchange.fluid_impulse_kg_m_s[axis]);
            assert_eq!(wall.exchange.wall_angular_impulse_kg_m2_s[axis], -wall.exchange.fluid_angular_impulse_kg_m2_s[axis]);
        }
        near(wall.exchange.fluid_work_j, raw_wall.fluid_work * source.density_kg_m3() * dx.powi(5) / dt.powi(2));
        assert_eq!(wall.exchange.wall_work_j, -wall.exchange.fluid_work_j);
        let mut change: f64 = 0.0;
        for z in 0..12 { for y in 0..16 { for x in 0..16 {
            if let (Some(a), Some(b)) = (native.sample_at(x, y, z), molecular.lattice().sample_at(x, y, z)) {
                for axis in 0..3 { change = change.max((a.velocity[axis] - b.velocity[axis]).abs()); }
            }
        }}}
        assert!(change > 1e-12, "SGS must change physical wall-driven evolution: {change}");
    }
}

#[test]
fn g4_wall_sgs_stencil_viscosity_storage_and_cancel_refusals_are_atomic() {
    let source = liquid(0, 300.0, None);
    let mut flow = primed(&source);
    let before = flow.clone();
    let values = velocities();
    let motion = motion(&values);
    let gate = CancelGate::new_clock_free();
    let config = config(SgsModel::Vreman { cv: 0.1 }, flow.scaling().cell_size_m());
    let mut kinematics = cache(flow.lattice());
    let mut sgs = FreeSurfaceSgsWorkspace3::default();
    let mut short = FreeSurfaceKinematicsWorkspace3::new(0, 0).unwrap();
    assert!(flow.step_with_planar_wall_velocity_field_and_gradient_sgs_si(
        motion, controls(), &gate, config, &mut short, &mut sgs,
    ).is_err());
    assert!(flow.step_with_planar_wall_velocity_field_and_gradient_sgs_si(
        motion, PlanarWallFieldControls3 { maximum_wall_faces: 1, ..controls() },
        &gate, config, &mut kinematics, &mut sgs,
    ).is_err());
    assert!(flow.step_with_planar_wall_velocity_field_and_gradient_sgs_si(
        motion, controls(), &gate,
        MaterialFreeSurfaceSgsConfig3 { maximum_eddy_to_molecular: 0.0, ..config },
        &mut kinematics, &mut sgs,
    ).is_err());
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(flow.step_with_planar_wall_velocity_field_and_gradient_sgs_si(
        motion, controls(), &cancelled, config, &mut kinematics, &mut sgs,
    ), Err(MaterialFreeSurfaceSgsError3::Material(MaterialLbmError::FreeSurface(FreeSurfaceError3::Cancelled)))));
    same(flow.lattice(), before.lattice());
    let mut thin = plan(&source).free_surface3(domain(), |_, _, z| z < 2).unwrap();
    let thin_before = thin.clone();
    assert!(matches!(thin.step_with_planar_wall_velocity_field_and_gradient_sgs_si(
        motion, controls(), &gate, config, &mut kinematics, &mut sgs,
    ), Err(MaterialFreeSurfaceSgsError3::Sgs(FreeSurfaceSgsError3::Cell { .. }))));
    same(thin.lattice(), thin_before.lattice());
    flow.step_with_planar_wall_velocity_field_and_gradient_sgs_si(
        motion, controls(), &gate, config, &mut kinematics, &mut sgs,
    ).unwrap();
}

#[test]
fn g4_late_unrepresentable_si_wall_receipt_does_not_publish_native_candidate() {
    let source = liquid(0, 300.0, None);
    // Deliberately extreme but admitted numerical units. Native angular impulse
    // is finite, while multiplication into physical angular impulse overflows.
    let setup = plan_free_surface_from_material(&source, 1e50, 2e102, 1e-54).unwrap();
    let mut flow = setup.free_surface3(
        MaterialFreeSurfaceDomain3 { cells: [4; 3], max_cells: 64, ..domain() },
        |_, _, _| true,
    ).unwrap();
    let before = flow.clone();
    let values = vec![[1e-55, 0.0]; 16];
    let motion = MaterialPlanarWallVelocityField3 {
        wall: PlanarWall3::ZMin,
        moment_origin_m: [1e208, 0.0, 0.0],
        tangential_velocity_m_s: &values,
    };
    let native_values: Vec<_> = values.iter().map(|v| {
        v.map(|u| u / setup.velocity_si(1.0).unwrap())
    }).collect();
    let gate = CancelGate::new_clock_free();
    let model = SgsModel::Wale { cw: 0.5 };
    let mut kinematics = cache(flow.lattice());
    let mut sgs = FreeSurfaceSgsWorkspace3::default();
    flow.lattice().planar_wall_velocity_field_gradient_sgs_candidate(
        PlanarWallVelocityField3 {
            wall: motion.wall,
            moment_origin: motion.moment_origin_m.map(|v| v / setup.cell_size_m()),
            tangential_velocity: &native_values,
        },
        controls(), &gate, native_config(model), &mut kinematics, &mut sgs,
    ).unwrap();
    assert!(matches!(flow.step_with_planar_wall_velocity_field_and_gradient_sgs_si(
        motion, controls(), &gate, config(model, setup.cell_size_m()), &mut kinematics, &mut sgs,
    ), Err(MaterialFreeSurfaceSgsError3::Material(MaterialLbmError::Invalid(_)))));
    same(flow.lattice(), before.lattice());
}
