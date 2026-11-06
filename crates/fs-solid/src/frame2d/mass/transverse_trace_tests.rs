use super::*;
use crate::frame2d::{Geometry, SectionModel};
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer};
use fs_material::state_point::{MaterialPropertySelection, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point};
use fs_qty::{Density, Pressure};

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

fn shear_frame(length: f64, angle: f64, shear_modulus: f64) -> Frame2d {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, value, dims) in [("density", 1000.0, Density::DIMS),
        ("young_modulus", 2e6, Pressure::DIMS), ("shear_modulus", shear_modulus, Pressure::DIMS)]
    {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained(),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance { source: "manufactured frame trace".into(),
                license: "test-only".into(), artifact: None },
        }).unwrap();
        requirements.push(ScalarPropertyRequirement::try_new(
            name, dims, ScalarAdmissibility::StrictlyPositive,
        ).unwrap());
    }
    let card = MaterialCard::assemble(MaterialStateId {
        chemistry: "synthetic shear beam".into(), phase: "solid".into(),
        process: "test".into(), revision: 0,
    }, claims, vec![]).unwrap();
    let material = resolve_material_state_point(
        &card, &QueryPoint::new(), &requirements, MaterialPropertySelection::SingleClaimOnly,
    ).unwrap();
    let section = ConcentricElasticSection::try_new(0.0, vec![ConcentricMaterialLayer {
        outer_radius_m: 0.02, material,
    }]).unwrap().with_shear_modulus("shear_modulus", 0.8).unwrap();
    let mut frame = Frame2d::with_mass_formulation(
        Geometry::Linear, FrameMassFormulation::ConsistentReference,
    );
    frame.add_node(0.0, 0.0);
    frame.add_node(length * angle.cos(), length * angle.sin());
    frame.add_material_shear_element(0, 1, 3, section).unwrap();
    frame
}

#[test]
fn timoshenko_trace_preserves_native_interpolation_rigid_motion_and_point_work() {
    for modulus in [200.0, 2000.0, 2e6] {
        for length in [0.3, 2.7] {
            for angle in [0.0_f64, 0.37, 1.2] {
                let frame = shear_frame(length, angle, modulus);
                let section = frame.material_shear_section(0).unwrap();
                let phi = 12.0 * section.elastic().bending_stiffness_n_m2()
                    / (section.shear_stiffness_n() * length * length);
                let beta = 1.0 / (1.0 + phi);
                let fraction = phi / (1.0 + phi);
                let local = [0.3, -0.2, -0.4, 0.17];
                let velocity = [-angle.sin() * local[0], angle.cos() * local[0], local[1],
                    -angle.sin() * local[2], angle.cos() * local[2], local[3]];
                let delta = local[2] - local[0] - 0.5 * length * (local[1] + local[3]);
                // Independent polynomial solution of constant-section
                // Timoshenko kinematics with the four authored endpoint DOFs.
                let coefficients = [local[0], length * local[1] + fraction * delta,
                    0.5 * length * (local[3] - local[1]) + 3.0 * beta * delta,
                    -2.0 * beta * delta];
                for xi in [0.0, 0.13, 0.5, 0.91, 1.0] {
                    let (dofs, weights, actual_length) = frame.reference_transverse_trace(0, xi).unwrap();
                    assert_eq!(dofs, [0, 1, 2, 3, 4, 5]);
                    near(actual_length, length);
                    let expected = coefficients.iter().rev().fold(0.0, |sum, coefficient| sum * xi + coefficient);
                    let observed = weights.iter().zip(velocity).map(|(b, v)| b * v).sum::<f64>();
                    near(observed, expected);
                    let force = -2.3;
                    near(weights.iter().zip(velocity).map(|(b, v)| (force * b) * v).sum(), force * expected);
                    let rotation = 0.012;
                    let rigid = [0.0, 0.0, rotation, -angle.sin() * rotation * length,
                        angle.cos() * rotation * length, rotation];
                    near(weights.iter().zip(rigid).map(|(b, v)| b * v).sum(), rotation * length * xi);
                }
            }
        }
    }
}

#[test]
fn timoshenko_trace_and_independent_section_rotation_recover_native_inertia() {
    let points = [0.06943184420297371, 0.33000947820757187, 0.6699905217924281, 0.9305681557970262];
    let weights = [0.17392742256872693, 0.32607257743127307, 0.32607257743127307, 0.17392742256872693];
    for modulus in [200.0, 2000.0, 2e6] {
        for length in [0.3, 1.0, 2.7] {
            for angle in [0.0_f64, 0.37, 1.2] {
                let frame = shear_frame(length, angle, modulus);
                let section = frame.material_shear_section(0).unwrap();
                let phi = 12.0 * section.elastic().bending_stiffness_n_m2()
                    / (section.shear_stiffness_n() * length * length);
                let beta = 1.0 / (1.0 + phi);
                let line_mass = section.elastic().line_density_kg_m();
                let spin = section.elastic().layer_integrals().iter()
                    .map(|layer| layer.density_kg_m3() * layer.second_moment_m4()).sum::<f64>();
                let local = [0.3, -0.2, -0.4, 0.17];
                let velocity = [-angle.sin() * local[0], angle.cos() * local[0], local[1],
                    -angle.sin() * local[2], angle.cos() * local[2], local[3]];
                let delta = local[2] - local[0] - 0.5 * length * (local[1] + local[3]);
                let mut energy = 0.0;
                for (xi, weight) in points.into_iter().zip(weights) {
                    let (_, b, _) = frame.reference_transverse_trace(0, xi).unwrap();
                    let speed = b.iter().zip(velocity).map(|(b, v)| b * v).sum::<f64>();
                    // Section spin is not the centerline slope: shear contributes
                    // to translation while this independent rotation stays quadratic.
                    let turning = local[1] + (local[3] - local[1]) * xi
                        + (6.0 * beta * delta / length) * xi * (1.0 - xi);
                    energy += 0.5 * length * weight * (line_mass * speed * speed + spin * turning * turning);
                }
                near(energy, frame.reference_mass().unwrap().kinetic_energy(&velocity));
            }
        }
    }
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
