//! Independent consistent-mass oracles for simultaneous deformable-body impact.
use super::*;
use fs_couple::finite_solid::contact::solids::{
    FiniteSolidContactAssemblyControls, FiniteSolidContactAssemblyError,
    FiniteSolidContactAssemblyReport, FiniteSolidPairContact, SolidContactPoint, contact_candidate,
};

const LOWER: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 0.0, -1.0],
];
fn make_body<'a>(
    material: &'a FiniteSolidMaterial,
    nodes: &'a [[f64; 3]],
    tets: &'a [[usize; 4]],
    fixed: &'a [(usize, f64)],
    velocity: &[[f64; 3]],
    time: f64,
    cx: &Cx<'_>,
) -> FiniteSolidBody<'a> {
    let mut p = problem(material, fixed, &[]);
    p.nodes_m = nodes;
    p.tetrahedra = tets;
    FiniteSolidBody::new(
        p,
        &vec![[0.0; 3]; nodes.len()],
        velocity,
        time,
        Default::default(),
        cx,
    )
    .unwrap()
}
fn point(body: usize, tet: usize, weights: [f64; 4]) -> SolidContactPoint {
    SolidContactPoint {
        body,
        point: FiniteSolidPoint::new(tet, weights).unwrap(),
    }
}
fn pair(id: u64, first: SolidContactPoint, second: SolidContactPoint) -> FiniteSolidPairContact {
    FiniteSolidPairContact {
        id,
        first,
        second,
        normal_world: Vec3::new(0.0, 0.0, 1.0),
        tangent_world: Vec3::new(1.0, 0.0, 0.0),
        activation_tolerance_m: 0.0,
        friction_coefficient: 0.0,
        normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 0.0 },
    }
}
fn assembled_audit(r: &FiniteSolidContactAssemblyReport) {
    same(r.linear_momentum_residual_kg_m_s, Vec3::ZERO, 2e-10);
    same(r.angular_momentum_residual_kg_m2_s, Vec3::ZERO, 2e-10);
    near(r.work_energy_residual_j, 0.0, 2e-10);
    near(r.dissipation_residual_j, 0.0, 2e-10);
    near(r.endpoint_velocity_residual_m_s, 0.0, 2e-10);
    assert!(r.absolute_linear_residual_kg_m_s < 1e-9);
    assert!(r.absolute_angular_residual_kg_m2_s < 1e-9);
    assert!(r.absolute_contact_couple_kg_m2_s < 1e-9);
    assert!(r.absolute_work_residual_j < 1e-8);
}

#[test]
fn different_source_densities_exchange_elastic_and_frictional_impulses_under_a_boost() {
    with_cx(|_, cx| {
        let a_material = bind(&source_card(6.0, 1.0, false, false));
        let b_material = bind(&source_card(12.0, 1.0, false, false));
        for (mu, e) in [(0.0, 1.0), (0.3, 0.0)] {
            let bodies = vec![
                make_body(
                    &a_material,
                    &NODES,
                    &TETS,
                    &[],
                    &[[9.0, -4.0, 15.0]; 4],
                    0.0,
                    cx,
                ),
                make_body(
                    &b_material,
                    &LOWER,
                    &TETS,
                    &[],
                    &[[8.0, -4.0, 16.0]; 4],
                    0.0,
                    cx,
                ),
            ];
            let before: Vec<_> = bodies.iter().map(|b| b.native_state().clone()).collect();
            let mut c = pair(
                7,
                point(0, 0, [0.5, 0.25, 0.25, 0.0]),
                point(1, 0, [0.5, 0.25, 0.25, 0.0]),
            );
            c.friction_coefficient = mu;
            c.normal_law = NormalVelocityLaw::NewtonRestitution { coefficient: e };
            let (next, r) =
                contact_candidate(&bodies, &[c], [0.2, -0.3, 0.4], Default::default(), cx).unwrap();
            // Unit-tet masses are 1 and 2 kg. M^-1=(20 I-4 11^T)/mass;
            // weights (.5,.25,.25,0) give point mobility 3.5/mass per axis.
            let mobility = 5.25;
            let jn = (1.0 + e) / mobility;
            let jt = -mu * jn;
            near(r.contact_mobility_per_kg[0], mobility, 2e-12);
            near(r.contact.impulses[0][0], jn, 2e-11);
            near(r.contact.impulses[0][1], jt, 2e-11);
            near(r.contact.impulses[0][2], 0.0, 2e-12);
            let expected_work =
                jn * (-1.0 + 0.5 * mobility * jn) + jt * (1.0 + 0.5 * mobility * jt);
            near(r.kinetic_change_j, expected_work, 2e-10);
            near(r.contact.dissipated_energy, -expected_work, 2e-11);
            same(
                vec(r.solid_impulses[0].applied.linear_kg_m_s),
                Vec3::new(jt, 0.0, jn),
                2e-12,
            );
            same(
                vec(r.solid_impulses[1].applied.linear_kg_m_s),
                Vec3::new(-jt, 0.0, -jn),
                2e-12,
            );
            assembled_audit(&r);
            for i in 0..2 {
                assert_eq!(bodies[i].native_state(), &before[i]);
                assert_eq!(next[i].native_state().q, before[i].q);
                assert_eq!(next[i].native_state().a, before[i].a);
                assert_eq!(next[i].native_state().t, before[i].t);
                assert_eq!(next[i].native_state().history, before[i].history);
            }
        }
    });
}

#[test]
fn a_shared_middle_solid_keeps_signed_cross_contact_mobility_and_input_order_independence() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let upper = NODES.map(|x| [x[0], x[1], x[2] + 1.0]);
        let bodies = vec![
            make_body(
                &material,
                &upper,
                &TETS,
                &[],
                &[[0.0, 0.0, -1.0]; 4],
                0.0,
                cx,
            ),
            make_body(&material, &NODES, &TETS, &[], &[[0.0; 3]; 4], 0.0, cx),
            make_body(
                &material,
                &LOWER,
                &TETS,
                &[],
                &[[0.0, 0.0, 1.0]; 4],
                0.0,
                cx,
            ),
        ];
        let a = pair(
            20,
            point(0, 0, [1.0, 0.0, 0.0, 0.0]),
            point(1, 0, [0.0, 0.0, 0.0, 1.0]),
        );
        let b = pair(
            10,
            point(1, 0, [1.0, 0.0, 0.0, 0.0]),
            point(2, 0, [1.0, 0.0, 0.0, 0.0]),
        );
        let (next, r) =
            contact_candidate(&bodies, &[a, b], [0.0; 3], Default::default(), cx).unwrap();
        assert_eq!(r.contact_ids, vec![10, 20]);
        // Each vertex inverse mass is16/kg. Different vertices on the middle
        // body have -4/kg mobility, and opposite incidence signs make W12=+4.
        near(r.contact_mobility_per_kg[0], 32.0, 2e-12);
        near(r.contact_mobility_per_kg[3], 4.0, 2e-12);
        near(r.contact_mobility_per_kg[18], 4.0, 2e-12);
        near(r.contact_mobility_per_kg[21], 32.0, 2e-12);
        for impulse in &r.contact.impulses {
            near(impulse[0], 1.0 / 36.0, 2e-11);
        }
        near(r.contact.dissipated_energy, 1.0 / 36.0, 2e-11);
        assembled_audit(&r);
        let permuted =
            contact_candidate(&bodies, &[b, a], [0.0; 3], Default::default(), cx).unwrap();
        assert_eq!(r.contact.impulses, permuted.1.contact.impulses);
        for i in 0..3 {
            assert_eq!(next[i].native_state(), permuted.0[i].native_state());
        }
    });
}

#[test]
fn same_body_touching_arms_retain_both_signed_trace_cross_terms() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        // Two tetrahedral arms share node0. Their distinct tip nodes1 and4
        // coincide at the contact; both reference tetrahedra have volume1/6.
        let nodes = [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [1.0, 0.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.0, 0.0, -1.0],
        ];
        let tets = [[0, 1, 2, 3], [0, 4, 5, 6]];
        let mut velocity = [[0.0; 3]; 7];
        velocity[1][2] = -0.5;
        velocity[4][2] = 0.5;
        let bodies = vec![make_body(&material, &nodes, &tets, &[], &velocity, 0.0, cx)];
        let c = pair(
            5,
            point(0, 0, [0.0, 1.0, 0.0, 0.0]),
            point(0, 1, [0.0, 1.0, 0.0, 0.0]),
        );
        let (_, r) = contact_candidate(&bodies, &[c], [0.0; 3], Default::default(), cx).unwrap();
        // Eliminating the shared-node block gives tip diagonals15.5 and cross
        // mobility+.5: the relative response is15.5+15.5-2*.5=30/kg.
        near(r.contact_mobility_per_kg[0], 30.0, 3e-11);
        near(r.contact.impulses[0][0], 1.0 / 30.0, 2e-11);
        near(r.kinetic_change_j, -1.0 / 60.0, 2e-11);
        assert_eq!(r.solid_mass_solves[0].len(), 6);
        same(
            vec(r.solid_impulses[0].applied.linear_kg_m_s),
            Vec3::ZERO,
            2e-12,
        );
        same(
            vec(r.solid_impulses[0].applied.angular_kg_m2_s),
            Vec3::ZERO,
            2e-12,
        );
        assembled_audit(&r);
    });
}

#[test]
fn constrained_counterpart_preserves_the_native_external_support_impulse() {
    with_cx(|_, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let fixed: Vec<_> = (0..9).map(|dof| (dof, 0.0)).collect();
        let bodies = vec![
            make_body(
                &material,
                &NODES,
                &TETS,
                &[],
                &[[0.0, 0.0, -1.0]; 4],
                0.0,
                cx,
            ),
            make_body(&material, &LOWER, &TETS, &fixed, &[[0.0; 3]; 4], 0.0, cx),
        ];
        let c = pair(
            3,
            point(0, 0, [0.5, 0.25, 0.25, 0.0]),
            point(1, 0, [0.5, 0.25, 0.25, 0.0]),
        );
        let (_, r) = contact_candidate(&bodies, &[c], [0.0; 3], Default::default(), cx).unwrap();
        let j = 2.0 / 7.0;
        near(r.contact.impulses[0][0], j, 2e-11);
        same(
            vec(r.support_impulse.linear_kg_m_s),
            Vec3::new(0.0, 0.0, j),
            2e-11,
        );
        same(
            vec(r.support_impulse.angular_kg_m2_s),
            Vec3::new(0.25 * j, -0.25 * j, 0.0),
            2e-11,
        );
        near(r.solid_impulses[1].kinetic_change_j, 0.0, 2e-13);
        near(r.solid_impulses[1].impulse_work_j, 0.0, 2e-13);
        assembled_audit(&r);
    });
}

#[test]
fn finite_gap_friction_couple_and_late_refusals_never_publish_body_changes() {
    with_cx(|gate, cx| {
        let material = bind(&source_card(6.0, 1.0, false, false));
        let lower = LOWER.map(|x| [x[0], x[1], x[2] - 1e-4]);
        let bodies = vec![
            make_body(
                &material,
                &NODES,
                &TETS,
                &[],
                &[[1.0, 0.0, -1.0]; 4],
                0.0,
                cx,
            ),
            make_body(&material, &lower, &TETS, &[], &[[0.0; 3]; 4], 0.0, cx),
        ];
        let before: Vec<_> = bodies.iter().map(|b| b.native_state().clone()).collect();
        let mut c = pair(
            11,
            point(0, 0, [0.5, 0.25, 0.25, 0.0]),
            point(1, 0, [0.5, 0.25, 0.25, 0.0]),
        );
        c.activation_tolerance_m = 2e-4;
        c.friction_coefficient = 0.3;
        let Err(FiniteSolidContactAssemblyError::Audit(r)) =
            contact_candidate(&bodies, &[c], [0.0; 3], Default::default(), cx)
        else {
            panic!("finite-gap friction couple must refuse")
        };
        assert!(r.angular_momentum_residual_kg_m2_s.y.abs() > 1e-6);
        near(
            r.angular_momentum_residual_kg_m2_s.y,
            1e-4 * r.contact.impulses[0][1],
            2e-11,
        );
        // Opposite sliding directions create two noncentral contact couples
        // whose global moment cancels. Each pair still violates the internal
        // angular balance and must refuse independently of body solve errors.
        let shifted_upper = NODES.map(|x| [x[0] + 2.0, x[1], x[2]]);
        let shifted_lower = lower.map(|x| [x[0] + 2.0, x[1], x[2]]);
        let mut opposing_bodies = bodies.clone();
        opposing_bodies.push(make_body(
            &material,
            &shifted_upper,
            &TETS,
            &[],
            &[[-1.0, 0.0, -1.0]; 4],
            0.0,
            cx,
        ));
        opposing_bodies.push(make_body(
            &material,
            &shifted_lower,
            &TETS,
            &[],
            &[[0.0; 3]; 4],
            0.0,
            cx,
        ));
        let opposing_before: Vec<_> = opposing_bodies
            .iter()
            .map(|body| body.native_state().clone())
            .collect();
        let mut opposite = c;
        opposite.id = 12;
        opposite.first.body = 2;
        opposite.second.body = 3;
        let Err(FiniteSolidContactAssemblyError::Audit(r)) = contact_candidate(
            &opposing_bodies,
            &[c, opposite],
            [0.0; 3],
            Default::default(),
            cx,
        ) else {
            panic!("cancelling pair couples must still refuse")
        };
        same(r.angular_momentum_residual_kg_m2_s, Vec3::ZERO, 2e-10);
        assert!(r.absolute_angular_residual_kg_m2_s < 1e-9);
        assert!(r.absolute_contact_couple_kg_m2_s > 1e-6);
        near(
            r.absolute_contact_couple_kg_m2_s,
            1e-4 * (r.contact.impulses[0][1].abs() + r.contact.impulses[1][1].abs()),
            2e-11,
        );
        for (body, old) in opposing_bodies.iter().zip(&opposing_before) {
            assert_eq!(body.native_state(), old);
        }
        let mut limited = FiniteSolidContactAssemblyControls::default();
        limited.maximum_solids = 1;
        assert!(contact_candidate(&bodies, &[c], [0.0; 3], limited, cx).is_err());
        gate.request();
        assert!(matches!(
            contact_candidate(&bodies, &[c], [0.0; 3], Default::default(), cx),
            Err(FiniteSolidContactAssemblyError::Cancelled)
        ));
        for i in 0..2 {
            assert_eq!(bodies[i].native_state(), &before[i]);
        }
    });
}
