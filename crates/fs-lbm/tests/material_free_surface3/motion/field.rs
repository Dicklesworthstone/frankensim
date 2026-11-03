use super::*;
use fs_lbm::d3q19::freesurface3::{PlanarWallFieldControls3, PlanarWallVelocityField3};
use fs_lbm::material::free_surface::MaterialPlanarWallVelocityField3;

fn controls() -> PlanarWallFieldControls3 {
    PlanarWallFieldControls3 {
        maximum_domain_cells: 16 * 16 * 12,
        maximum_wall_faces: 16 * 16,
        maximum_relative_mass_defect: 1e-12,
        maximum_relative_face_mass_defect: 1e-12,
    }
}

fn field() -> Vec<[f64; 2]> {
    (0..256)
        .map(|i| {
            [
                0.01 + 0.001 * (i % 16) as f64 - 0.0002 * (i / 16) as f64,
                -0.004 + 0.0003 * (i / 16) as f64,
            ]
        })
        .collect()
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

#[test]
fn g1_material_face_field_matches_independent_si_impulse_couple_work_and_pressure_loads() {
    let source = liquid(0, 300.0, None);
    let setup = plan(&source).with_pressure_reference(100_000.0).unwrap();
    let original = setup.free_surface3(domain(), |_, _, z| z < 6).unwrap();
    let velocities = field();
    let motion = MaterialPlanarWallVelocityField3 {
        wall: PlanarWall3::ZMin,
        moment_origin_m: [-0.003, 0.009, 0.002],
        tangential_velocity_m_s: &velocities,
    };
    let gate = CancelGate::new_clock_free();
    let (mut candidate, report) = original
        .planar_wall_velocity_field_candidate_si(motion, controls(), &gate)
        .unwrap();
    let velocity_unit = setup.velocity_si(1.0).unwrap();
    let native_velocities: Vec<_> = velocities
        .iter()
        .map(|v| v.map(|u| u / velocity_unit))
        .collect();
    let (native, raw) = original
        .lattice()
        .planar_wall_velocity_field_candidate(
            PlanarWallVelocityField3 {
                wall: motion.wall,
                moment_origin: motion.moment_origin_m.map(|x| x / setup.cell_size_m()),
                tangential_velocity: &native_velocities,
            },
            controls(),
            &gate,
        )
        .unwrap();
    assert_eq!(
        candidate.lattice().grid().state_bits(),
        native.grid().state_bits()
    );
    assert_eq!(candidate.scaling().material(), &source);
    let dx = setup.cell_size_m();
    let cell_mass = source.density_kg_m3() * dx.powi(3);
    let mut total_work = 0.0;
    let mut total_impulse = [0.0; 3];
    for (index, face) in report.faces.iter().enumerate() {
        let u = velocities[index];
        for d in 0..2 {
            close(
                face.exchange.fluid_impulse_kg_m_s[d],
                cell_mass * u[d] / 3.0,
                2e-18,
            );
            close(face.fluid_tangential_response_kg[d], cell_mass / 3.0, 2e-18);
            close(face.velocity_m_s[d], u[d], 2e-17);
        }
        let work = cell_mass * (u[0] * u[0] + u[1] * u[1]) / 3.0;
        close(face.exchange.fluid_work_j, work, 2e-20);
        assert_eq!(face.exchange.wall_work_j, -face.exchange.fluid_work_j);
        total_work += work;
        let force_moment = cross(
            core::array::from_fn(|d| face.center_m[d] - motion.moment_origin_m[d]),
            face.exchange.fluid_impulse_kg_m_s,
        );
        for d in 0..3 {
            close(
                force_moment[d] + face.fluid_couple_at_face_kg_m2_s[d],
                face.exchange.fluid_angular_impulse_kg_m2_s[d],
                2e-18,
            );
            assert_eq!(
                face.wall_couple_at_face_kg_m2_s[d],
                -face.fluid_couple_at_face_kg_m2_s[d]
            );
            total_impulse[d] += face.exchange.fluid_impulse_kg_m_s[d];
        }
        close(
            face.exchange.absolute_wall_force_n.unwrap()[2],
            -100_000.0 * dx * dx,
            2e-15,
        );
        assert!(face.reflected_mass_defect_kg.abs() < 1e-19);
    }
    close(report.exchange.fluid_work_j, total_work, 2e-18);
    for d in 0..3 {
        close(
            report.exchange.fluid_impulse_kg_m_s[d],
            total_impulse[d],
            2e-16,
        );
    }
    close(report.wetted_area_m2, 256.0 * dx * dx, 2e-18);
    close(report.mass_after_kg, report.mass_before_kg, 2e-16);
    close(
        report.maximum_relative_face_mass_defect,
        raw.maximum_relative_face_mass_defect,
        1e-20,
    );

    // The pressure datum changes neither native populations nor mechanical
    // work. Its actual normal load and moment remain an explicit reservoir.
    let no_reference = plan(&source)
        .free_surface3(domain(), |_, _, z| z < 6)
        .unwrap();
    let (without, unreferenced) = no_reference
        .planar_wall_velocity_field_candidate_si(motion, controls(), &gate)
        .unwrap();
    assert_eq!(
        without.lattice().grid().state_bits(),
        candidate.lattice().grid().state_bits()
    );
    assert_eq!(
        unreferenced.exchange.fluid_work_j,
        report.exchange.fluid_work_j
    );
    assert!(unreferenced.exchange.absolute_wall_force_n.is_none());
    assert!(
        unreferenced
            .faces
            .iter()
            .all(|f| f.exchange.absolute_wall_torque_nm.is_none())
    );
    let reference_force = [0.0, 0.0, -100_000.0 * report.wetted_area_m2];
    let centroid = [8.0 * dx, 8.0 * dx, 0.0];
    let reference_torque = cross(
        core::array::from_fn(|d| centroid[d] - motion.moment_origin_m[d]),
        reference_force,
    );
    for d in 0..3 {
        close(
            report.exchange.absolute_wall_torque_nm.unwrap()[d]
                - report.exchange.pressure_perturbation_wall_torque_nm[d],
            reference_torque[d],
            2e-15,
        );
    }
    let first = candidate.lattice().grid().state_bits();
    candidate
        .step_with_planar_wall_velocity_field_si(motion, controls(), &gate)
        .unwrap();
    assert_ne!(first, candidate.lattice().grid().state_bits());
    assert_eq!(candidate.scaling().material(), &source);
}

#[test]
fn g4_local_source_speed_and_late_liquid_speed_refusals_retain_material_history() {
    let source = liquid(0, 300.0, None);
    let mut original = plan(&source)
        .free_surface3(domain(), |_, _, z| z < 6)
        .unwrap();
    let gate = CancelGate::new_clock_free();
    let mut velocities = field();
    velocities[117] = [0.11, 0.0];
    let motion = MaterialPlanarWallVelocityField3 {
        wall: PlanarWall3::ZMin,
        moment_origin_m: [0.0; 3],
        tangential_velocity_m_s: &velocities,
    };
    let before = original.lattice().grid().state_bits();
    let mass = original.ledger_mass_kg().unwrap();
    assert!(
        matches!(original.step_with_planar_wall_velocity_field_si(motion,controls(),&gate),Err(MaterialLbmError::Invalid(reason)) if reason.contains("local wall field"))
    );
    assert_eq!(original.lattice().grid().state_bits(), before);
    assert_eq!(original.ledger_mass_kg().unwrap(), mass);

    let setup = plan_free_surface_from_material(&source, 1e-3, 2e-4, 0.02).unwrap();
    let mut accelerated = setup
        .free_surface3(
            MaterialFreeSurfaceDomain3 {
                gravity_m_s2: [0.0, 0.0, -100.0],
                ..domain()
            },
            |_, _, z| z < 6,
        )
        .unwrap();
    let velocities = vec![[0.005, 0.0]; 256];
    let motion = MaterialPlanarWallVelocityField3 {
        tangential_velocity_m_s: &velocities,
        ..motion
    };
    let bits = accelerated.lattice().grid().state_bits();
    let mass = accelerated.ledger_mass_kg().unwrap();
    let counters = accelerated.lattice().conversions();
    // Initial Guo velocity is admitted. The actual forced native endpoint
    // exceeds the original 0.02 m/s source scenario after the fluid solve.
    assert!(
        matches!(accelerated.step_with_planar_wall_velocity_field_si(motion,controls(),&gate),Err(MaterialLbmError::Invalid(reason)) if reason.contains("liquid exceeds"))
    );
    assert_eq!(accelerated.lattice().grid().state_bits(), bits);
    assert_eq!(accelerated.ledger_mass_kg().unwrap(), mass);
    assert_eq!(accelerated.lattice().conversions(), counters);
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(matches!(
        original.step_with_planar_wall_velocity_field_si(motion, controls(), &cancelled),
        Err(MaterialLbmError::FreeSurface(FreeSurfaceError3::Cancelled))
    ));
    let mut direct = plan(&source)
        .free_surface3(domain(), |_, _, z| z < 6)
        .unwrap();
    let a = original
        .step_with_planar_wall_velocity_field_si(motion, controls(), &gate)
        .unwrap();
    let b = direct
        .step_with_planar_wall_velocity_field_si(motion, controls(), &gate)
        .unwrap();
    assert_eq!(a, b);
    assert_eq!(
        original.lattice().grid().state_bits(),
        direct.lattice().grid().state_bits()
    );
}
