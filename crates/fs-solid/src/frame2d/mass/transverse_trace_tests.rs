use super::*;
use crate::frame2d::{Geometry, SectionModel};

fn frame(length: f64, angle: f64) -> Frame2d {
    let mut frame = Frame2d::with_mass_formulation(Geometry::Linear, FrameMassFormulation::ConsistentReference);
    frame.add_node(0.0, 0.0);
    frame.add_node(length*angle.cos(), length*angle.sin());
    frame.add_element(0, 1, 3, 3.4, &|| SectionModel::Elastic { ea: 200.0, ei: 8.0 }).unwrap();
    frame
}
fn near(a: f64, b: f64) {
    assert!((a-b).abs() < 2e-12*(1.0+a.abs()+b.abs()), "{a:e} != {b:e}");
}

#[test]
fn trace_endpoints_and_rotation_weights_use_original_dofs() {
    let frame = frame(2.0, 0.0);
    let (dofs, left, length) = frame.reference_transverse_trace(0, 0.0).unwrap();
    let (_, right, _) = frame.reference_transverse_trace(0, 1.0).unwrap();
    let (_, mid, _) = frame.reference_transverse_trace(0, 0.5).unwrap();
    assert_eq!(dofs, [0,1,2,3,4,5]);
    assert_eq!(left, [0.0,1.0,0.0,0.0,0.0,0.0]);
    assert_eq!(right, [0.0,0.0,0.0,0.0,1.0,0.0]);
    assert_eq!(mid, [0.0,0.5,0.25,0.0,0.5,-0.25]);
    assert_eq!(length, 2.0);
}

#[test]
fn rigid_rotation_is_reproduced_on_rotated_members() {
    for angle in [0.0_f64, 0.37, 1.2] {
        let length = 2.7;
        let frame = frame(length, angle);
        let theta = 0.012;
        let velocity = [0.0, 0.0, theta, -angle.sin()*theta*length, angle.cos()*theta*length, theta];
        for xi in [0.0, 0.13, 0.5, 0.91, 1.0] {
            let (_, weights, _) = frame.reference_transverse_trace(0, xi).unwrap();
            near(weights.iter().zip(velocity).map(|(b,v)| b*v).sum(), theta*length*xi);
        }
    }
}

#[test]
fn four_point_trace_quadrature_recovers_original_consistent_transverse_inertia() {
    let points = [0.06943184420297371, 0.33000947820757187, 0.6699905217924281, 0.9305681557970262];
    let weights = [0.17392742256872693, 0.32607257743127307, 0.32607257743127307, 0.17392742256872693];
    for angle in [0.0_f64, 0.37, 1.2] {
        for length in [0.3, 1.0, 2.7] {
            let frame = frame(length, angle);
            // Pure transverse translations plus rotations: no axial kinetic energy.
            let velocity = [-angle.sin()*0.3, angle.cos()*0.3, -0.2,
                angle.sin()*0.4, -angle.cos()*0.4, 0.17];
            let mut energy = 0.0;
            for (xi, weight) in points.into_iter().zip(weights) {
                let (_, b, l) = frame.reference_transverse_trace(0, xi).unwrap();
                let speed: f64 = b.iter().zip(velocity).map(|(b,v)| b*v).sum();
                energy += 0.5*3.4*l*weight*speed*speed;
            }
            near(energy, frame.reference_mass().unwrap().kinetic_energy(&velocity));
        }
    }
}

#[test]
fn one_physical_trace_preserves_transverse_force_and_nodal_moment_work() {
    let frame = frame(2.0, 0.7);
    let (_, b, _) = frame.reference_transverse_trace(0, 0.3).unwrap();
    let velocity = [0.1, -0.3, 0.2, 0.7, -0.1, 0.4];
    let force = -2.3;
    let point_power = force*b.iter().zip(velocity).map(|(b,v)| b*v).sum::<f64>();
    let nodal_power: f64 = b.iter().zip(velocity).map(|(b,v)| (force*b)*v).sum();
    near(point_power, nodal_power);
    assert_ne!(b[2], 0.0);
    assert_ne!(b[5], 0.0);
}

#[test]
fn unsupported_trace_requests_refuse_without_changing_the_frame() {
    let mut frame = frame(1.0, 0.0);
    let saved = frame.displacements().to_vec();
    assert!(frame.reference_transverse_trace(1, 0.5).is_err());
    for xi in [-0.1, 1.1, f64::NAN, f64::INFINITY] {
        assert!(frame.reference_transverse_trace(0, xi).is_err());
    }
    assert_eq!(frame.displacements(), saved);
    frame.geometry = Geometry::Corotational;
    assert!(frame.reference_transverse_trace(0, 0.5).is_err());
}

#[test]
fn segment_query_keeps_original_orientation_and_never_moves_geometry() {
    let mut frame=frame(2.0,0.0);
    assert_eq!(frame.reference_member_segment(0).unwrap(),[[0.0,0.0],[2.0,0.0]]);
    frame.u[3]=0.3;
    assert_eq!(frame.reference_member_segment(0).unwrap(),[[0.0,0.0],[2.0,0.0]]);
    let reverse=frame.add_element(1,0,3,1.0,&||SectionModel::Elastic {ea:2.0,ei:1.0}).unwrap();
    assert_eq!(frame.reference_member_segment(reverse).unwrap(),[[2.0,0.0],[0.0,0.0]]);
    assert!(frame.reference_member_segment(99).is_err());
    frame.nodes[0][0]=f64::NAN;
    assert!(frame.reference_member_segment(0).is_err());
    assert_eq!(frame.u[3],0.3);
}
