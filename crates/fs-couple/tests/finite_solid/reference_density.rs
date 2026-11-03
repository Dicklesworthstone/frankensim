//! G1/G3: source-bound graded matter reaches actual motion and point impulses.
use super::*;
use fs_couple::finite_solid::impulse::FiniteSolidPoint;
use fs_ga::{Twist, Vec3};
use fs_mbd_material::density_field::{DensityVertexBinding, TetrahedralDensityBody};
use fs_solid::hyper3::pressure::ClosedTetPressureSurface;
use std::sync::Arc;

const RHO: [f64; 4] = [3.0, 15.0, 3.0, 3.0];

fn array(value: Vec3) -> [f64; 3] {
    [value.x, value.y, value.z]
}
fn vector(value: [f64; 3]) -> Vec3 {
    Vec3::new(value[0], value[1], value[2])
}
fn density_field(rho: [f64; 4], nodes: &[[f64; 3]; 4]) -> Arc<TetrahedralDensityBody> {
    let cards = rho.map(|density| source_card(density, 1.0, false, false));
    let context = point(300.0);
    let bindings: Vec<_> = (0..4)
        .map(|node| DensityVertexBinding {
            position_m: vector(nodes[node]),
            material: &cards[node],
            point: &context,
            selection: MaterialPropertySelection::SingleClaimOnly,
        })
        .collect();
    Arc::new(TetrahedralDensityBody::resolve(&bindings, &TETS, 4, || false).unwrap())
}

#[test]
fn graded_motion_matches_reference_and_deformed_material_inertia() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let density = density_field(RHO, &NODES);
        // Finite affine shear/stretch; density must transport as rho_0/J.
        let u: [[f64; 3]; 4] = NODES.map(|x| [0.2 * x[0] + 0.1 * x[1], -0.1 * x[1], 0.3 * x[2]]);
        let current: Vec<_> = (0..4)
            .map(|a| vector(std::array::from_fn(|c| NODES[a][c] + u[a][c])))
            .collect();
        let twist = Twist::new(Vec3::new(0.3, -0.4, 0.5), Vec3::new(0.2, 0.1, -0.3));
        let velocity: Vec<_> = current
            .iter()
            .map(|&x| array(twist.linear + twist.angular.cross(x)))
            .collect();
        let body = FiniteSolidBody::new_with_reference_density_field(
            problem(&material, &[], &[]),
            Arc::clone(&density),
            &u,
            &velocity,
            0.0,
            Default::default(),
            cx,
        )
        .unwrap();
        let transported = density.deform(&current, || false).unwrap();
        let oracle = transported.inertia();
        near(body.reference_mass_kg(), 1.0, 2e-14);
        near(body.reference_mass_kg(), oracle.mass(), 2e-14);
        near(
            body.kinetic_energy_j().unwrap(),
            oracle.kinetic_energy(twist).unwrap(),
            3e-14,
        );
        let momentum = body.momentum_about([0.0; 3], cx).unwrap();
        let expected = oracle.momentum(twist).unwrap();
        for c in 0..3 {
            near(momentum.linear_kg_m_s[c], array(expected.force)[c], 3e-14);
            near(
                momentum.angular_kg_m2_s[c],
                array(expected.torque)[c],
                3e-14,
            );
        }
        assert!(std::ptr::eq(
            body.reference_density_field().unwrap(),
            density.as_ref()
        ));
        for (vertex, &rho) in body
            .reference_density_field()
            .unwrap()
            .vertices()
            .iter()
            .zip(&RHO)
        {
            assert_eq!(vertex.density_kg_m3(), rho);
            assert_eq!(vertex.material().properties().len(), 1);
        }
        // Same total mass with uniform density has different rotational inertia.
        let averaged = FiniteSolidBody::new(
            problem(&material, &[], &[]),
            &u,
            &velocity,
            0.0,
            Default::default(),
            cx,
        )
        .unwrap();
        assert!(
            (averaged.kinetic_energy_j().unwrap() - body.kinetic_energy_j().unwrap()).abs() > 1e-4
        );
    });
}

#[test]
fn graded_inertia_drives_initial_acceleration_time_motion_and_contact_mobility() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let density = density_field(RHO, &NODES);
        let fixed: Vec<_> = (0..12)
            .filter(|&dof| dof != 3)
            .map(|dof| (dof, 0.0))
            .collect();
        let force = [(3, 0.3)];
        let body = FiniteSolidBody::new_with_reference_density_field(
            problem(&material, &fixed, &force),
            Arc::clone(&density),
            &[[0.0; 3]; 4],
            &[[0.0; 3]; 4],
            0.0,
            Default::default(),
            cx,
        )
        .unwrap();
        // Integral rho*N_1^2 = V*(sum(rho)+2*rho_1)/60 = 0.15 kg.
        let mass = 0.15;
        near(body.native_state().a[0], force[0].1 / mass, 2e-12);
        let gravity = body
            .reference_acceleration_load_n(&[[0.0, -10.0, 0.0]; 4], cx)
            .unwrap();
        assert_eq!(gravity.len(), 4);
        for (node, share) in [0.225, 0.325, 0.225, 0.225].into_iter().enumerate() {
            // These are all restrained Y rows; eliminating them would lose
            // the body's weight and its support-moment accounting.
            assert_eq!(gravity[node].0, 3 * node + 1);
            near(gravity[node].1, -10.0 * share, 2e-14);
        }
        let material_point = FiniteSolidPoint::new(0, [0.0, 1.0, 0.0, 0.0]).unwrap();
        near(
            body.point_directional_inverse_mass_per_kg(
                material_point,
                [1.0, 0.0, 0.0],
                Default::default(),
                cx,
            )
            .unwrap(),
            1.0 / mass,
            2e-12,
        );
        let (impulsed, report) = body
            .point_impulse_candidate(
                material_point,
                [0.03, 0.0, 0.0],
                [0.0; 3],
                Default::default(),
                cx,
            )
            .unwrap();
        near(impulsed.native_state().v[0], 0.2, 2e-12);
        near(report.kinetic_change_j, 0.5 * mass * 0.2_f64.powi(2), 2e-14);
        near(report.impulse_work_j, report.kinetic_change_j, 2e-14);
        assert!(std::ptr::eq(
            impulsed.reference_density_field().unwrap(),
            density.as_ref()
        ));
        let expected = midpoint_oracle(1.0, 0.0, 9.0, 1.0, 0.01, force[0].1);
        let (next, _) = body
            .step_candidate(0.01, &force, [0.0; 3], controls(1e-6), cx)
            .unwrap();
        near(next.native_state().q[0], expected.0 - 1.0, 2e-10);
        near(next.native_state().v[0], expected.1, 2e-9);
        assert!(std::ptr::eq(
            next.reference_density_field().unwrap(),
            density.as_ref()
        ));
        assert_eq!(body.native_state().q, vec![0.0]);
        assert_eq!(body.native_state().v, vec![0.0]);
    });
}

#[test]
fn graded_pressure_initialization_and_uniform_field_parity() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let fixed: Vec<_> = (0..12)
            .filter(|&dof| dof != 3)
            .map(|dof| (dof, 0.0))
            .collect();
        let surface = Arc::new(
            ClosedTetPressureSurface::new(&NODES, &TETS, TetAssemblyBudget::standard(), 4, cx)
                .unwrap(),
        );
        let graded = FiniteSolidBody::new_with_reference_density_field_and_pressure(
            problem(&material, &fixed, &[]),
            density_field(RHO, &NODES),
            Arc::clone(&surface),
            0.9,
            &[[0.0; 3]; 4],
            &[[0.0; 3]; 4],
            0.0,
            Default::default(),
            cx,
        )
        .unwrap();
        near(graded.native_state().a[0], -0.9 * VOLUME / 0.15, 2e-12);
        let uniform = FiniteSolidBody::new_with_reference_density_field_and_pressure(
            problem(&material, &fixed, &[]),
            density_field([6.0; 4], &NODES),
            Arc::clone(&surface),
            0.9,
            &[[0.0; 3]; 4],
            &[[0.0; 3]; 4],
            0.0,
            Default::default(),
            cx,
        )
        .unwrap();
        let original = FiniteSolidBody::new_with_pressure(
            problem(&material, &fixed, &[]),
            surface,
            0.9,
            &[[0.0; 3]; 4],
            &[[0.0; 3]; 4],
            0.0,
            Default::default(),
            cx,
        )
        .unwrap();
        assert_eq!(uniform.native_state(), original.native_state());
        assert_eq!(
            uniform.reference_mass_kg().to_bits(),
            original.reference_mass_kg().to_bits()
        );
        let (next, _) = graded
            .step_candidate_with_pressure(0.001, &[], [0.0; 3], 0.9, controls(1e-5), cx)
            .unwrap();
        assert_eq!(
            next.reference_density_field(),
            graded.reference_density_field()
        );
        assert_eq!(next.reference_mass_kg(), graded.reference_mass_kg());
    });
}

#[test]
fn density_binding_checks_reference_identity_sources_and_cancellation() {
    with_cx(|gate, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let mut shifted = NODES;
        shifted[0][0] = -0.1;
        let make = |density| {
            FiniteSolidBody::new_with_reference_density_field(
                problem(&material, &[], &[]),
                density,
                &[[0.0; 3]; 4],
                &[[0.0; 3]; 4],
                0.0,
                Default::default(),
                cx,
            )
        };
        assert!(make(density_field(RHO, &shifted)).is_err());
        // The density owner canonicalizes connectivity; native FE ordering survives.
        let density = density_field(RHO, &NODES);
        let reordered = [[3, 1, 0, 2]];
        let mut p = problem(&material, &[], &[]);
        p.tetrahedra = &reordered;
        let body = FiniteSolidBody::new_with_reference_density_field(
            p,
            Arc::clone(&density),
            &[[0.0; 3]; 4],
            &[[0.0; 3]; 4],
            0.0,
            Default::default(),
            cx,
        )
        .unwrap();
        assert_eq!(body.problem().tetrahedra, &reordered);
        assert!(body.reference_acceleration_load_n(&[], cx).is_err());
        assert!(
            body.reference_acceleration_load_n(&[[f64::NAN; 3]; 4], cx)
                .is_err()
        );
        let bad_card = source_card(6.0, 1.0, true, false);
        let context = point(300.0)
            .with_quantity("stress", QuantitySpec::dimensional(Pressure::DIMS), 0.0)
            .unwrap();
        let bindings: Vec<_> = NODES
            .iter()
            .map(|&x| DensityVertexBinding {
                position_m: vector(x),
                material: &bad_card,
                point: &context,
                selection: MaterialPropertySelection::SingleClaimOnly,
            })
            .collect();
        let evolving =
            Arc::new(TetrahedralDensityBody::resolve(&bindings, &TETS, 4, || false).unwrap());
        assert!(matches!(make(evolving), Err(FiniteSolidError::Source(_))));
        gate.request();
        assert!(make(density).is_err());
        assert!(
            body.reference_acceleration_load_n(&[[0.0; 3]; 4], cx)
                .is_err()
        );
    });
}
