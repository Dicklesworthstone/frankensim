//! Original current geometry and retained hot histories through detected impact.
use super::*;
use fs_couple::inelastic_impact::surface::*;
use fs_material::plastic::finite::thermal::{
    AffineThermalFiniteJ2, FiniteJ2TemperatureSlopes, FiniteJ2ThermalParameters,
};
use fs_solid::hyper3::plastic::thermal::{
    FiniteThermoplasticControls, FiniteThermoplasticTetSolid,
};

const FACE_POINT: TetImpulsePoint = TetImpulsePoint {
    tetrahedron: 0,
    weights: [0.5, 0.25, 0.25, 0.0],
};

fn sphere_pair(point: [f64; 3], velocity: [f64; 3]) -> (RigidBodyState, FiniteSolidSphere) {
    // Offset sphere beneath the z face, moving toward it at relative speed 1.
    let body = RigidBodyState::new(
        Pose::new(
            Vec3::new(point[0] - 0.125, point[1], point[2] - 0.25),
            UnitQuaternion::IDENTITY,
        )
        .unwrap(),
        Vec3::new(
            2.0 * velocity[0],
            2.0 * velocity[1],
            2.0 * (velocity[2] + 1.0),
        ),
        Vec3::ZERO,
    )
    .unwrap();
    let sphere = FiniteSolidSphere {
        id: 7,
        rigid_body: 0,
        center_body_m: Vec3::new(0.125, 0.0, 0.0),
        radius_m: 0.25,
        friction_coefficient: 0.0,
        normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 0.5 },
    };
    (body, sphere)
}

fn face_query(query: &FiniteSolidSurfaceQuery) {
    assert_eq!(query.contacts.len(), 1);
    near(query.contacts[0].gap_m, 0.0, 1e-13);
    let mapped = query.admitted_inelastic_contacts().unwrap();
    let FiniteSolidContactCounterpart::Solid(point) = mapped[0].counterpart else {
        panic!("sphere exterior contact lost its solid counterpart")
    };
    assert_eq!(point.tetrahedron, FACE_POINT.tetrahedron);
    for (actual, expected) in point.weights.into_iter().zip(FACE_POINT.weights) {
        near(actual, expected, 1e-13);
    }
    assert_eq!(mapped[0].id, query.contacts[0].contact.id);
}

pub(super) fn check_plastic_geometry(motion: &FinitePlasticTetMotion<'_>, cx: &Cx<'_>) {
    let trace = motion.point_kinematics(FACE_POINT, cx).unwrap();
    let (body, sphere) = sphere_pair(trace.position_m, trace.velocity_mps);
    let surface = FiniteSolidContactSurface::for_plastic(motion, 4, cx).unwrap();
    face_query(
        &surface
            .plastic_sphere_separations(motion, &[body], &[sphere], Default::default(), cx)
            .unwrap(),
    );
    let event = surface
        .plastic_sphere_contact_candidate(
            motion,
            &[],
            &[body],
            &[properties()],
            &[sphere],
            ORIGIN,
            Default::default(),
            mixed_controls(),
            cx,
        )
        .unwrap();
    face_query(&event.geometry);
    mixed_audit(&event.mechanical.report);
    assert_eq!(
        event.mechanical.solid.solid().element_states(),
        motion.solid().element_states()
    );
    assert_eq!(
        event.mechanical.solid.solid().displacement_m(),
        motion.solid().displacement_m()
    );
    assert_eq!(
        event.mechanical.solid.reference_nodal_density_kg_m3(),
        Some(RHO.as_slice())
    );
}

pub(super) fn check_visco_geometry(motion: &FiniteViscoTetMotion<'_>, cx: &Cx<'_>) {
    let trace = motion.point_kinematics(FACE_POINT, cx).unwrap();
    let (body, sphere) = sphere_pair(trace.position_m, trace.velocity_mps);
    let surface = FiniteSolidContactSurface::for_visco(motion, 4, cx).unwrap();
    face_query(
        &surface
            .visco_sphere_separations(motion, &[body], &[sphere], Default::default(), cx)
            .unwrap(),
    );
    let event = surface
        .visco_sphere_contact_candidate(
            motion,
            &[],
            &[body],
            &[properties()],
            &[sphere],
            ORIGIN,
            Default::default(),
            mixed_controls(),
            cx,
        )
        .unwrap();
    face_query(&event.geometry);
    mixed_audit(&event.mechanical.report);
    assert_eq!(
        event.mechanical.solid.solid().element_states(),
        motion.solid().element_states()
    );
    assert_eq!(
        event.mechanical.solid.solid().time_s(),
        motion.solid().time_s()
    );
    assert_eq!(
        event.mechanical.solid.reference_nodal_density_kg_m3(),
        Some(RHO.as_slice())
    );
}

fn family() -> AffineThermalFiniteJ2 {
    AffineThermalFiniteJ2::new(
        FiniteJ2ThermalParameters {
            reference_temperature_k: 300.0,
            temperature_range_k: [250.0, 650.0],
            reference: FiniteJ2Parameters {
                shear_modulus_pa: 8000.0,
                bulk_modulus_pa: 16000.0,
                yield_stress_pa: 200.0,
                hardening_modulus_pa: 800.0,
            },
            slopes: FiniteJ2TemperatureSlopes {
                shear_modulus_pa_k: -6.0,
                bulk_modulus_pa_k: -8.0,
                yield_stress_pa_k: -0.4,
                hardening_modulus_pa_k: -1.0,
            },
            logarithmic_expansion_per_k: 3e-4,
            logarithmic_expansion_slope_per_k2: 1e-7,
        },
        FiniteJ2Limits {
            minimum_jacobian: 0.25,
            maximum_elastic_log_strain: 0.7,
            maximum_plastic_log_strain: 2.0,
            maximum_accumulated_plastic_strain: 3.0,
        },
    )
    .unwrap()
}

#[test]
fn cached_exterior_detects_hot_deformed_face_and_preserves_actual_momentum_history() {
    with_cx(|gate, cx| {
        let virgin = FiniteThermoplasticTetSolid::try_new(
            &NODES,
            &TETS,
            vec![family()],
            &[300.0],
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        let surface = FiniteSolidContactSurface::for_thermoplastic_tet(&virgin, 4, cx).unwrap();
        // Current x = (1.25 X, .75 Y, Z) + (.25,-.125,.5).
        // The old reference surface is far from the new contact face.
        let stretch = [1.25, 0.75, 1.0];
        let translation = [0.25, -0.125, 0.5];
        let fixed: Vec<_> = (0..12)
            .map(|dof| {
                (
                    dof,
                    (stretch[dof % 3] - 1.0) * NODES[dof / 3][dof % 3] + translation[dof % 3],
                )
            })
            .collect();
        let hot = virgin
            .advance(
                &[450.0],
                &fixed,
                &[],
                [0.0; 3],
                0.1,
                settings(),
                FiniteThermoplasticControls {
                    material: FinitePlasticStepControls {
                        equilibrium_work_tolerance_j: 2e-6,
                        maximum_endpoint_work_defect_j: 800.0,
                        yield_tolerance_pa: 1e-5,
                    },
                    maximum_numerical_kinetic_loss_j: 1.0,
                },
                cx,
            )
            .unwrap()
            .solid
            .with_inertia(&[6.0], &[[0.0; 3]; 4], cx)
            .unwrap();
        assert!(hot.element_states()[0].alpha > 0.1);
        let physical_point = [0.5625, 0.0625, 0.5];
        let (body, sphere) = sphere_pair(physical_point, [0.0; 3]);
        let event = surface
            .thermoplastic_tet_sphere_contact_candidate(
                &hot,
                &[],
                &[body],
                &[properties()],
                &[sphere],
                ORIGIN,
                Default::default(),
                mixed_controls(),
                cx,
            )
            .unwrap();
        face_query(&event.geometry);
        near(
            norm(
                event.geometry.contacts[0]
                    .solid_point_world_m
                    .sub(vector(physical_point)),
            ),
            0.0,
            1e-13,
        );
        let report = &event.mechanical.report;
        // M^-1=20I-4*11' for the unit-total-mass Tet. Its face mobility is 3.5.
        // Offset rigid arm gives (.125)^2*(I^-1)_yy=25/768 with the full tensor.
        let mobility = 3.5 + 0.5 + 25.0 / 768.0;
        let impulse = 1.5 / mobility;
        near(report.contact_mobility_per_kg[0], mobility, 1e-11);
        near(
            event.mechanical.rigid_bodies[0].linear_momentum_world().z,
            2.0 - impulse,
            1e-10,
        );
        near(report.contact.dissipated_energy, 0.375 / mobility, 1e-10);
        near(
            report.solid_impulse.applied.angular_kg_m2_s[0],
            (physical_point[1] - ORIGIN[1]) * impulse,
            1e-10,
        );
        near(
            report.solid_impulse.applied.angular_kg_m2_s[1],
            -(physical_point[0] - ORIGIN[0]) * impulse,
            1e-10,
        );
        mixed_audit(report);
        let after = &event.mechanical.solid;
        assert_eq!(after.element_states(), hot.element_states());
        assert_eq!(after.element_first_piola_pa(), hot.element_first_piola_pa());
        assert_eq!(after.displacement_m(), hot.displacement_m());
        assert_eq!(after.temperatures_k(), &[450.0]);
        assert_eq!(after.free_energy_j(), hot.free_energy_j());
        assert_eq!(after.internal_energy_j(), hot.internal_energy_j());
        assert_eq!(after.reference_mass_kg(), Some(1.0));

        let (penetrating, _) = sphere_pair(
            [
                physical_point[0],
                physical_point[1],
                physical_point[2] + 0.0625,
            ],
            [0.0; 3],
        );
        let raw = surface
            .thermoplastic_tet_sphere_separations(
                &hot,
                &[penetrating],
                &[sphere],
                Default::default(),
                cx,
            )
            .unwrap();
        near(raw.closest[0].gap_m, -0.0625, 1e-13);
        assert!(matches!(
            raw.admitted_inelastic_contacts(),
            Err(FiniteSolidSurfaceError::Penetration { .. })
        ));
        let mut other_nodes = NODES;
        other_nodes[3][2] = 2.0;
        let other = FiniteThermoplasticTetSolid::try_new(
            &other_nodes,
            &TETS,
            vec![family()],
            &[300.0],
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        assert!(matches!(
            surface.thermoplastic_tet_sphere_separations(&other, &[], &[], Default::default(), cx,),
            Err(FiniteSolidSurfaceError::ReferenceMismatch)
        ));
        assert!(
            surface
                .thermoplastic_tet_sphere_separations(
                    &hot,
                    &[body],
                    &[sphere],
                    FiniteSolidSurfaceControls {
                        maximum_query_visits: 0,
                        ..Default::default()
                    },
                    cx,
                )
                .is_err()
        );
        assert_eq!(hot.velocity_mps().unwrap(), [[0.0; 3]; 4]);
        gate.request();
        assert!(matches!(
            surface.thermoplastic_tet_sphere_separations(
                &hot,
                &[body],
                &[sphere],
                Default::default(),
                cx,
            ),
            Err(FiniteSolidSurfaceError::Cancelled)
        ));
    });
}
