use super::*;
use fs_couple::bernoulli_aperture::modal::pressure::MaterialPressureShape;
use fs_ports::{ModuleWindowAudit, PortedModule};

fn shape(valve: &ModalAperture, blade: usize, uniform: bool) -> MaterialPressureShape {
    let plate = &valve.material_plates()[blade];
    let values = plate.source().chart().mesh.nodes.iter()
        .map(|&(x, _)| if uniform { 1.0 } else { 2.0*x/0.025-1.0 }).collect();
    MaterialPressureShape::from_plate(plate, values, "manufactured signed pressure gradient".into())
}
fn bind(valve: ModalAperture, blade: usize, uniform: bool) -> ModalAperture {
    let trace = shape(&valve, blade, uniform);
    valve.with_material_pressure_shape(blade, trace, 100000, &CancelGate::new_clock_free()).unwrap()
}

#[test]
fn affine_pressure_uses_consistent_face_quadrature_and_original_modes() {
    let original = bound(vec![plate(900., 2e4), plate(1100., 2e4)], 1e6);
    let modes = original.spec().modes.clone();
    let valve = bind(bind(original, 0, false), 1, true);
    assert_eq!(valve.spec().modes, modes); // Uniform acoustic/swept projections are untouched.
    for blade in 0..2 {
        let plate = &valve.material_plates()[blade];
        let port = valve.material_pressure_ports()[blade].as_ref().unwrap();
        let mesh = &plate.source().chart().mesh;
        assert_eq!(port.shape().specimen_identity(), plate.source().specimen_identity());
        let mut net_area = 0.0;
        for (local, mode) in plate.mode_range().enumerate() {
            let mut reference = 0.0;
            for tri in &mesh.tris {
                let [a, b, c] = tri.map(|n| mesh.nodes[n]);
                let area = 0.5*((b.0-a.0)*(c.1-a.1)-(b.1-a.1)*(c.0-a.0));
                // Independent degree-two exact triangle quadrature. Both traces
                // are affine, so their product is integrated exactly by this rule.
                for weights in [[2.0/3.0,1.0/6.0,1.0/6.0],
                    [1.0/6.0,2.0/3.0,1.0/6.0], [1.0/6.0,1.0/6.0,2.0/3.0]] {
                    let phi: f64 = tri.iter().zip(weights)
                        .map(|(&n,w)| w*plate.displacement_shapes()[n][local]).sum();
                    let psi: f64 = tri.iter().zip(weights)
                        .map(|(&n,w)| w*port.shape().nodal_shape()[n]).sum();
                    reference += area/3.0*phi*psi;
                    if local == 0 { net_area += area/3.0*psi; }
                }
            }
            assert!((reference-port.modal_area_m2()[mode]).abs() < 5e-18);
            if blade == 1 {
                assert!((reference-modes[mode].pressure_area_m2).abs() < 3e-16);
            }
        }
        if blade == 0 {
            assert!(net_area.abs() < 1e-18, "manufactured field has zero net force");
            assert!(port.modal_area_m2().iter().any(|a| a.abs() > 1e-8),
                "a balanced pressure gradient must still bend the blade");
        }
    }
}

#[test]
fn spatial_pressure_uses_same_motion_and_closes_work_without_extra_duct_flow() {
    let mut valve = bind(bind(bound(vec![plate(900., 2e4), plate(1100., 2e4)], 1e6), 0, false), 1, true);
    let dt = valve.spec().time_step_s;
    for step in 0..64 {
        let before = valve.modal_states().to_vec();
        let pressures = if step < 16 { [8.0, -2.0] } else { [0.0; 2] };
        let frame = valve.step_with_pressures(drive(0.0), [0.0; 2], pressures).unwrap();
        let mut work = 0.0;
        for blade in 0..2 {
            let port = valve.material_pressure_ports()[blade].as_ref().unwrap();
            let flow: f64 = port.modal_area_m2().iter().zip(before.iter().zip(&frame.modes))
                .map(|(a,(old,new))| a*(new.displacement_m-old.displacement_m)/dt).sum();
            let actual = frame.spatial_pressure[blade];
            assert!((actual.flow_m3_s-flow).abs() < 1e-20);
            assert!((actual.work_j-dt*pressures[blade]*flow).abs() < 1e-24);
            work += actual.work_j;
        }
        assert!((frame.aperture.mechanical_work_j-work).abs() < 1e-24);
        assert!((frame.aperture.bore_flow_m3_s-frame.aperture.jet_flow_m3_s
            -frame.aperture.swept_flow_m3_s).abs() < 1e-16);
        let scale = frame.aperture.stored_energy_j + frame.aperture.dissipated_energy_j.abs()
            + frame.aperture.mechanical_work_j.abs() + frame.aperture.pressure_work_j.abs();
        assert!(frame.aperture.balance_residual_j().abs() < 1e-9*scale.max(1e-30));
        if step == 0 {
            assert_eq!(valve.port_schemas().len(), 6);
            assert_eq!(valve.last_window_transfers().unwrap().len(), 6);
            let ModuleWindowAudit::Available(audit) = valve.audit_last_window().unwrap() else {
                panic!("complete pressure/force port audit");
            };
            assert!(audit.violations().is_empty());
        }
    }
    assert!(valve.modal_states().iter().any(|s| s.displacement_m != 0.0));
    let before = valve.modal_states().to_vec();
    let steps = valve.accepted_steps();
    let energy = valve.stored_energy_j();
    let transfers = valve.last_window_transfers().unwrap().to_vec();
    assert!(valve.step_with_pressures(drive(0.0), [0.0; 2], [f64::NAN, 0.0]).is_err());
    assert_eq!(valve.modal_states(), before);
    assert_eq!(valve.accepted_steps(), steps);
    assert_eq!(valve.stored_energy_j(), energy);
    assert_eq!(valve.last_window_transfers().unwrap(), transfers);
}

#[test]
fn released_spatial_source_preserves_original_pressure_trajectory_bitwise() {
    let mut original = bound(vec![plate(900., 2e4)], 1e6);
    let mut attached = bind(bound(vec![plate(900., 2e4)], 1e6), 0, false);
    for _ in 0..64 {
        let a = original.step(drive(5.0), [0.0; 2]).unwrap();
        let b = attached.step_with_pressures(drive(5.0), [0.0; 2], [0.0; 2]).unwrap();
        assert_eq!(a.aperture, b.aperture);
        assert_eq!(a.modes, b.modes);
    }
}

#[test]
fn mismatched_incomplete_duplicate_and_late_pressure_sources_refuse() {
    let gate = CancelGate::new_clock_free();
    let new = || bound(vec![plate(900., 2e4)], 1e6);
    let other = bound(vec![plate(1100., 2e4)], 1e6);
    assert!(new().with_material_pressure_shape(0, shape(&other, 0, true), 100000, &gate).is_err());
    for values in [vec![], vec![0.0; 21], vec![f64::NAN; 21], vec![1.0; 20]] {
        let v = new();
        let input = MaterialPressureShape::from_plate(&v.material_plates()[0], values, "test".into());
        assert!(v.with_material_pressure_shape(0, input, 100000, &gate).is_err());
    }
    let v = new();
    let input = shape(&v, 0, true);
    assert!(v.with_material_pressure_shape(0, input, 1, &gate).is_err());
    let v = bind(new(), 0, true);
    let input = shape(&v, 0, true);
    assert!(v.with_material_pressure_shape(0, input, 100000, &gate).is_err());
    let mut v = new();
    assert!(v.step_with_pressures(drive(0.0), [0.0; 2], [1.0, 0.0]).is_err());
    v.step(drive(0.0), [0.0; 2]).unwrap();
    let input = shape(&v, 0, true);
    assert!(v.with_material_pressure_shape(0, input, 100000, &gate).is_err());
    let v = new();
    let input = shape(&v, 0, true);
    gate.request();
    assert!(v.with_material_pressure_shape(0, input, 100000, &gate).is_err());
}
