//! Independent finite-stretch motion, work, material substitution and objectivity.
//! These are authored numerical fixtures, not measured material validation.
#[path = "finite_solid/thermal.rs"]
mod thermal;
#[path = "finite_solid/expansion.rs"]
mod expansion;
#[path = "finite_solid/pressure.rs"]
mod pressure;
#[path = "finite_solid/impulse.rs"]
mod impulse;
#[path = "finite_solid/contact.rs"]
mod contact;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_couple::finite_solid::*;
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, ConstitutiveModelCard, InitialStatePolicy, InterpolationPolicy, LawId,
    LawParameter, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::state_point::MaterialPropertySelection;
use fs_qty::{Density, Dims, Pressure, QuantitySpec, Temperature};
use fs_solid::hyper3::{HyperTetProblem, HyperTetSettings, HyperelasticLawField};
use fs_solid::linear3::TetAssemblyBudget;

const NODES: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0],
];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];
const VOLUME: f64 = 1.0 / 6.0;

fn with_cx<R>(f: impl FnOnce(&CancelGate, &Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&gate, &Cx::new(
            &gate, arena, StreamKey { seed: 3, kernel_id: 17, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic,
        ))
    })
}

fn point(temperature_k: f64) -> QueryPoint {
    QueryPoint::new().with_quantity(
        "T", QuantitySpec::dimensional(Temperature::DIMS), temperature_k,
    ).unwrap()
}

fn source_card(rho: f64, scale: f64, density_stress_axis: bool, model_stress_axis: bool)
    -> MaterialCard
{
    let artifact = fs_blake3::hash_bytes(
        format!("synthetic finite-solid source: density={rho}; modulus scale={scale}").as_bytes(),
    );
    let provenance = Provenance {
        source: "authored compressible Ogden alpha=2 fixture; not measured data".into(),
        license: "CC0-1.0".into(), artifact: Some(artifact),
    };
    let domain = ValidityDomain::unconstrained().with_quantity(
        "T", QuantitySpec::dimensional(Temperature::DIMS), 280.0, 320.0,
    );
    let stress_domain = |domain: ValidityDomain| domain.with_quantity(
        "stress", QuantitySpec::dimensional(Pressure::DIMS), 0.0, 100.0,
    );
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new("density", Density::DIMS),
        value: PropertyValue::Scalar { value: rho, dims: Density::DIMS },
        validity: if density_stress_axis { stress_domain(domain.clone()) } else { domain.clone() },
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![], provenance: provenance.clone(),
    }).unwrap();
    let model = ConstitutiveModelCard {
        law: LawId("ogden-decoupled".into()), law_version: 1,
        parameters: [
            ("stretch_limit", 3.0, Dims::NONE),
            ("volumetric_kind", 0.0, Dims::NONE),
            ("term_count", 1.0, Dims::NONE),
            ("mu_0", 15.0 * scale, Pressure::DIMS),
            ("alpha_0", 2.0, Dims::NONE),
            ("kappa", 40.0 * scale, Pressure::DIMS),
        ].into_iter().map(|(name, value, dims)|
            (name.into(), LawParameter { value, dims })).collect(),
        state_schema_version: 1, initial_state: InitialStatePolicy::ZeroInternalState,
        validity: if model_stress_axis { stress_domain(domain) } else { domain },
        sources: vec![artifact], provenance,
    };
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: format!("synthetic rho={rho} scale={scale}"), phase: "solid".into(),
            process: "authored analytic fixture".into(), revision: 0,
        }, claims, vec![model],
    ).unwrap()
}

fn bind(card: &MaterialCard) -> FiniteSolidMaterial {
    FiniteSolidMaterial::resolve(
        card, 0, &point(300.0), MaterialPropertySelection::SingleClaimOnly,
    ).unwrap()
}

fn problem<'a>(
    material: &'a FiniteSolidMaterial,
    prescribed: &'a [(usize, f64)],
    force: &'a [(usize, f64)],
) -> HyperTetProblem<'a, FiniteSolidMaterial> {
    HyperTetProblem {
        nodes_m: &NODES, tetrahedra: &TETS,
        material: HyperelasticLawField::Uniform(material),
        prescribed_m: prescribed, nodal_forces_n: force, body_force_n_m3: [0.0; 3],
        budget: TetAssemblyBudget::standard(), settings: HyperTetSettings::default(),
    }
}

fn controls(tolerance_j: f64) -> FiniteSolidControls {
    let mut controls = FiniteSolidControls {
        energy_tolerance_j: tolerance_j, ..FiniteSolidControls::default()
    };
    controls.solve.max_newton_iterations = 32;
    controls.solve.newton.absolute_tolerance = 1.0e-12;
    controls.solve.newton.relative_tolerance = 1.0e-12;
    controls
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!((actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, tolerance {tolerance:e}");
}

// F = diag(lambda, 1, 1). Independent closed forms differentiated by hand,
// deliberately not calling the production energy, Piola, tangent or integrator.
fn nominal_stress(lambda: f64, scale: f64) -> f64 {
    scale * (10.0 * (lambda.powf(1.0 / 3.0) - lambda.powf(-5.0 / 3.0))
        + 40.0 * (lambda - 1.0))
}
fn strain_energy(lambda: f64, scale: f64) -> f64 {
    VOLUME * scale * (
        7.5 * (lambda.powf(4.0 / 3.0) + 2.0 * lambda.powf(-2.0 / 3.0) - 3.0)
        + 20.0 * (lambda - 1.0).powi(2)
    )
}
fn midpoint_oracle(lambda: f64, velocity: f64, rho: f64, scale: f64, h: f64, force: f64)
    -> (f64, f64)
{
    let mass = rho * VOLUME / 10.0;
    let residual = |next: f64| {
        2.0 * mass * (next - lambda - h * velocity) / (h * h)
            + VOLUME * nominal_stress(0.5 * (lambda + next), scale) - force
    };
    let (mut low, mut high) = (0.1, 2.8);
    assert!(residual(low) < 0.0 && residual(high) > 0.0);
    for _ in 0..100 {
        let middle = 0.5 * (low + high);
        if residual(middle) < 0.0 { low = middle; } else { high = middle; }
    }
    let next = 0.5 * (low + high);
    (next, 2.0 * (next - lambda) / h - velocity)
}

#[test]
fn g1_actual_finite_stretch_motion_and_work_match_independent_scalar_equations() {
    with_cx(|_, cx| {
        let card = source_card(6.0, 1.0, false, false);
        let material = bind(&card);
        let fixed: Vec<_> = (0..12).filter(|dof| *dof != 3).map(|dof| (dof, 0.0)).collect();
        let force = [(3, 0.2)];
        let p = problem(&material, &fixed, &force);
        let mut u = [[0.0; 3]; 4]; u[1][0] = 0.3;
        let mut v = [[0.0; 3]; 4]; v[1][0] = 0.1;
        let initial = FiniteSolidBody::new(p, &u, &v, 0.0, Default::default(), cx).unwrap();
        assert!(initial.initialization_report().converged);
        assert_eq!(initial.free_dofs(), &[3]);
        let mut body = initial.clone();
        let (mut lambda, mut velocity) = (1.3, 0.1);
        let mass = 6.0 * VOLUME / 10.0;
        for step in 0..8 {
            let before = body.native_state().clone();
            let expected = midpoint_oracle(lambda, velocity, 6.0, 1.0, 0.01, 0.2);
            let (next, report) = body.step_candidate(
                0.01, &force, [0.0; 3], controls(1.0e-5), cx,
            ).unwrap();
            assert_eq!(body.native_state(), &before);
            near(next.native_state().q[0], expected.0 - 1.0, 2.0e-9);
            near(next.native_state().v[0], expected.1, 2.0e-8);
            near(report.strain_before_j, strain_energy(lambda, 1.0), 1.0e-9);
            near(report.strain_after_j, strain_energy(expected.0, 1.0), 1.0e-9);
            near(report.kinetic_before_j, 0.5 * mass * velocity * velocity, 1.0e-9);
            near(report.kinetic_after_j, 0.5 * mass * expected.1 * expected.1, 1.0e-9);
            near(report.source_work_j, 0.2 * (expected.0 - lambda), 1.0e-10);
            let expected_defect = strain_energy(expected.0, 1.0) - strain_energy(lambda, 1.0)
                + 0.5 * mass * (expected.1 * expected.1 - velocity * velocity)
                - 0.2 * (expected.0 - lambda);
            near(report.integration_defect_j, expected_defect, 2.0e-9);
            near(next.reference_mass_kg(), 1.0, 2.0e-15);
            assert_eq!(report.support_reactions_n.len(), 11);
            assert_eq!(next.native_state().steps, step + 1);
            assert_eq!(next.native_state().history.len(), step + 1);
            (lambda, velocity) = expected;
            body = next;
        }
        assert_eq!(initial.native_state().steps, 0);
        near(initial.native_state().q[0], 0.3, 0.0);
        assert_eq!(material.law().model_card(), &card.models()[0]);
        assert_eq!(material.density().card_identity(), card.content_hash());
        assert_eq!(material.density().material(), card.id());
        assert_eq!(material.law().context(), &material.density().query());
        assert_eq!(material.model_index(), 0);
    });
}

#[test]
fn g3_material_density_and_stiffness_change_motion_on_the_same_reference_mesh() {
    with_cx(|_, cx| {
        let fixed: Vec<_> = (0..12).filter(|dof| *dof != 3).map(|dof| (dof, 0.0)).collect();
        let mut u = [[0.0; 3]; 4]; u[1][0] = 0.25;
        let mut result = Vec::new();
        for (rho, scale, h) in [
            (6.0, 1.0, 0.02), (24.0, 1.0, 0.02),
            (24.0, 1.0, 0.04), (6.0, 2.0, 0.02),
        ] {
            let card = source_card(rho, scale, false, false);
            let material = bind(&card);
            let body = FiniteSolidBody::new(
                problem(&material, &fixed, &[]), &u, &[[0.0; 3]; 4], 0.0,
                Default::default(), cx,
            ).unwrap();
            let (next, _) = body.step_candidate(h, &[], [0.0; 3], controls(1.0e-3), cx).unwrap();
            let expected = midpoint_oracle(1.25, 0.0, rho, scale, h, 0.0);
            near(next.native_state().q[0], expected.0 - 1.0, 1.0e-9);
            near(next.native_state().v[0], expected.1, 1.0e-8);
            near(next.reference_mass_kg(), rho * VOLUME, 1.0e-14);
            result.push(next.native_state().q[0]);
        }
        assert!(result[1] > result[0]); // four times the mass responds more slowly
        near(result[2], result[0], 1.0e-9); // h/sqrt(rho) scaling
        assert!(result[3] < result[0]); // stiffer source restores faster
    });
}

#[test]
fn g1_rigid_rotation_and_translating_acceleration_preserve_mass_and_work() {
    with_cx(|_, cx| {
        let card = source_card(6.0, 1.0, false, false);
        let material = bind(&card);
        let velocity = [0.2, -0.1, 0.3];
        let gravity = [0.4, -0.2, 0.1];
        let body_force = gravity.map(|value| 6.0 * value);
        // Exact ninety-degree rotation around z, then rigid translation.
        let u = NODES.map(|x| [-x[1] - x[0] + 0.125, x[0] - x[1] + 0.25, -0.5]);
        let v = [velocity; 4];
        let mut p = problem(&material, &[], &[]);
        p.body_force_n_m3 = body_force;
        let body = FiniteSolidBody::new(p, &u, &v, 0.0, Default::default(), cx).unwrap();
        let h = 0.02;
        let (next, report) = body.step_candidate(h, &[], body_force, controls(1.0e-8), cx).unwrap();
        let q1 = next.displacement_m();
        let v1 = next.velocity_m_s();
        for node in 0..4 {
            for axis in 0..3 {
                near(q1[node][axis],
                    u[node][axis] + h * velocity[axis] + 0.5 * h * h * gravity[axis], 1.0e-9);
                near(v1[node][axis], velocity[axis] + h * gravity[axis], 1.0e-8);
            }
        }
        let work: f64 = (0..3).map(|axis|
            gravity[axis] * (h * velocity[axis] + 0.5 * h * h * gravity[axis])).sum();
        near(report.source_work_j, work, 1.0e-10);
        near(report.kinetic_after_j - report.kinetic_before_j, work, 1.0e-10);
        near(report.strain_before_j, 0.0, 1.0e-10);
        near(report.strain_after_j, 0.0, 1.0e-10);
        near(report.integration_defect_j, 0.0, 1.0e-10);
        near(next.reference_mass_kg(), 1.0, 2.0e-15);
        assert!(report.support_reactions_n.is_empty());
    });
}

#[test]
fn g1_heterogeneous_regions_use_their_own_reference_mass_under_one_body_force() {
    with_cx(|_, cx| {
        let cards = [
            source_card(6.0, 1.0, false, false),
            source_card(12.0, 2.0, false, false),
        ];
        let materials = cards.each_ref().map(bind);
        let mut nodes = NODES.to_vec();
        nodes.extend(NODES.map(|x| [x[0] + 3.0, x[1], x[2]]));
        let tets = [[0, 1, 2, 3], [4, 5, 6, 7]];
        let body_force = [3.0, 0.6, -0.9];
        let p = HyperTetProblem {
            nodes_m: &nodes, tetrahedra: &tets,
            material: HyperelasticLawField::PerElement(&materials),
            prescribed_m: &[], nodal_forces_n: &[], body_force_n_m3: body_force,
            budget: TetAssemblyBudget::standard(), settings: HyperTetSettings::default(),
        };
        let body = FiniteSolidBody::new(
            p, &[[0.0; 3]; 8], &[[0.0; 3]; 8], 0.0, Default::default(), cx,
        ).unwrap();
        let h = 0.02;
        let (next, report) = body.step_candidate(h, &[], body_force, controls(1.0e-8), cx).unwrap();
        let u = next.displacement_m();
        let v = next.velocity_m_s();
        let mut expected_work = 0.0;
        for (region, density) in [6.0, 12.0].into_iter().enumerate() {
            for node in 4 * region..4 * region + 4 {
                for axis in 0..3 {
                    let acceleration = body_force[axis] / density;
                    near(u[node][axis], 0.5 * h * h * acceleration, 1.0e-9);
                    near(v[node][axis], h * acceleration, 1.0e-8);
                }
            }
            expected_work += body_force.iter()
                .map(|force| 0.5 * VOLUME * h * h * force * force / density).sum::<f64>();
            assert_eq!(materials[region].density().card_identity(), cards[region].content_hash());
        }
        near(next.reference_mass_kg(), 3.0, 1.0e-14);
        near(report.source_work_j, expected_work, 1.0e-10);
        near(report.kinetic_after_j, expected_work, 1.0e-10);
        near(report.strain_after_j, 0.0, 1.0e-10);
        near(report.integration_defect_j, 0.0, 1.0e-10);
    });
}

#[test]
fn g4_source_domain_massless_nodes_and_late_refusal_preserve_accepted_motion() {
    let card = source_card(6.0, 1.0, false, false);
    assert!(FiniteSolidMaterial::resolve(
        &card, 1, &point(300.0), MaterialPropertySelection::SingleClaimOnly,
    ).is_err());
    assert!(FiniteSolidMaterial::resolve(
        &card, 0, &point(330.0), MaterialPropertySelection::SingleClaimOnly,
    ).is_err());
    let context = point(300.0).with_quantity(
        "stress", QuantitySpec::dimensional(Pressure::DIMS), 20.0,
    ).unwrap();
    for source in [
        source_card(6.0, 1.0, true, false),
        source_card(6.0, 1.0, false, true),
    ] {
        // A matching numerical query does not authorize freezing a mechanical axis.
        assert!(FiniteSolidMaterial::resolve(
            &source, 0, &context, MaterialPropertySelection::SingleClaimOnly,
        ).is_err());
    }
    with_cx(|gate, cx| {
        let material = bind(&card);
        let fixed: Vec<_> = (0..12).filter(|dof| *dof != 3).map(|dof| (dof, 0.0)).collect();
        let p = problem(&material, &fixed, &[]);
        let mut u = [[0.0; 3]; 4]; u[1][0] = 0.3;
        let zero = [[0.0; 3]; 4];
        let body = FiniteSolidBody::new(p, &u, &zero, 0.0, Default::default(), cx).unwrap();
        let before = body.native_state().clone();
        let error = body.step_candidate(0.1, &[], [0.0; 3], controls(1.0e-18), cx).unwrap_err();
        assert!(matches!(error, FiniteSolidError::Energy { .. }), "{error}");
        assert_eq!(body.native_state(), &before);
        let (next, report) = body.step_candidate(0.1, &[], [0.0; 3], controls(1.0), cx).unwrap();
        let (stretch, velocity) = midpoint_oracle(1.3, 0.0, 6.0, 1.0, 0.1, 0.0);
        let defect = strain_energy(stretch, 1.0) - strain_energy(1.3, 1.0)
            + 0.5 * (6.0 * VOLUME / 10.0) * velocity * velocity;
        assert!(defect.abs() > 1.0e-10);
        near(report.integration_defect_j, defect, 1.0e-8);
        near(next.native_state().q[0], stretch - 1.0, 1.0e-9);
        assert_eq!(body.native_state(), &before);

        let mut inverted = u; inverted[1][0] = -1.1;
        assert!(FiniteSolidBody::new(p, &inverted, &zero, 0.0, Default::default(), cx).is_err());
        let mut too_stretched = u; too_stretched[1][0] = 2.1;
        assert!(FiniteSolidBody::new(p, &too_stretched, &zero, 0.0, Default::default(), cx).is_err());
        let mut moving_support = zero; moving_support[0][0] = 0.1;
        assert!(FiniteSolidBody::new(
            p, &u, &moving_support, 0.0, Default::default(), cx,
        ).is_err());
        let mut nodes = NODES.to_vec(); nodes.push([2.0, 2.0, 2.0]);
        let mut massless = problem(&material, &[], &[]);
        massless.nodes_m = &nodes;
        assert!(matches!(FiniteSolidBody::new(
            massless, &[[0.0; 3]; 5], &[[0.0; 3]; 5], 0.0, Default::default(), cx,
        ), Err(FiniteSolidError::Invalid(_))));
        assert!(FiniteSolidBody::new(
            p, &u, &zero, 0.0, FiniteSolidInitialization {
                maximum_workspace_components: 1, ..Default::default()
            }, cx,
        ).is_err());

        let mut no_steps = controls(1.0); no_steps.maximum_steps = 0;
        assert!(body.step_candidate(0.01, &[], [0.0; 3], no_steps, cx).is_err());
        gate.request();
        assert!(body.step_candidate(0.01, &[], [0.0; 3], controls(1.0), cx).is_err());
        assert_eq!(body.native_state(), &before);
    });
}
