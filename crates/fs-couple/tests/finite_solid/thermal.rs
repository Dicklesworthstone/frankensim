//! Source-declared thermal length ratios enter native force, tangent and motion.
use super::*;
use fs_couple::finite_solid::thermal::{MeanThermalReference, ThermalReference};
use fs_material::hyper::cards::HyperelasticMaterial;
use fs_material::state_point::INVERSE_TEMPERATURE_DIMS;

const MEAN: &str = "mean_engineering_thermal_expansion_coefficient";

pub(super) fn thermal_card(alpha: f64, scale: f64, stress_axis: bool) -> MaterialCard {
    let base = source_card(6.0, scale, false, false);
    let provenance = base.models()[0].provenance.clone();
    let mut claims = ClaimSet::new();
    let temperature = QuantitySpec::dimensional(Temperature::DIMS);
    let endpoint_stretch = 1.0 + 10.0 * alpha;
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new("density", Density::DIMS),
        value: PropertyValue::Curve {
            abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
            knots: vec![(300.0, 6.0), (310.0, 6.0 / endpoint_stretch.powi(3))],
            dims: Density::DIMS,
        },
        validity: ValidityDomain::unconstrained().with_quantity("T", temperature, 300.0, 310.0),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::LinearInside,
        observations: vec![], provenance: provenance.clone(),
    }).unwrap();
    let mut domain = ValidityDomain::unconstrained()
        .with_quantity("T_reference", temperature, 300.0, 310.0)
        .with_quantity("T_endpoint", temperature, 300.0, 310.0);
    if stress_axis {
        domain = domain.with_quantity("stress", QuantitySpec::dimensional(Pressure::DIMS), 0.0, 100.0);
    }
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(MEAN, INVERSE_TEMPERATURE_DIMS),
        value: PropertyValue::Scalar { value: alpha, dims: INVERSE_TEMPERATURE_DIMS },
        validity: domain, uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![], provenance,
    }).unwrap();
    MaterialCard::assemble(base.id().clone(), claims, base.models().to_vec()).unwrap()
}

fn thermal_query(reference: f64, endpoint: f64) -> QueryPoint {
    let temperature = QuantitySpec::dimensional(Temperature::DIMS);
    QueryPoint::new().with_quantity("T_reference", temperature, reference).unwrap()
        .with_quantity("T_endpoint", temperature, endpoint).unwrap()
}

fn reference(card: &MaterialCard) -> MeanThermalReference {
    MeanThermalReference::resolve(
        card, &PropertyKey::new(MEAN, INVERSE_TEMPERATURE_DIMS),
        &thermal_query(300.0, 310.0), "T_reference", "T_endpoint",
        MaterialPropertySelection::SingleClaimOnly,
    ).unwrap()
}

pub(super) fn prepared(card: &MaterialCard) -> FiniteSolidMaterial {
    FiniteSolidMaterial::resolve_mean_thermal(
        card, 0, reference(card), MaterialPropertySelection::SingleClaimOnly,
    ).unwrap()
}

// Independent derivatives of W0=s^3 W(diag(lambda,1,1)/s), for the authored
// Ogden alpha=2 law. This distinguishes volume, stress and tangent factors.
fn thermal_energy(lambda: f64, s: f64, scale: f64) -> f64 {
    let jth = s.powi(3);
    scale * jth * (
        7.5 * (lambda.powf(4.0 / 3.0) + 2.0 * lambda.powf(-2.0 / 3.0) - 3.0)
            + 20.0 * (lambda / jth - 1.0).powi(2)
    )
}
fn thermal_stress(lambda: f64, s: f64, scale: f64) -> f64 {
    scale * (10.0 * s.powi(3) * (lambda.powf(1.0 / 3.0) - lambda.powf(-5.0 / 3.0))
        + 40.0 * (lambda / s.powi(3) - 1.0))
}
fn thermal_tangent(lambda: f64, s: f64, scale: f64) -> f64 {
    scale * ((10.0 / 3.0) * s.powi(3)
        * (lambda.powf(-2.0 / 3.0) + 5.0 * lambda.powf(-8.0 / 3.0))
        + 40.0 / s.powi(3))
}
fn thermal_midpoint(lambda: f64, velocity: f64, s: f64, scale: f64, h: f64)
    -> (f64, f64)
{
    let mass = 6.0 * VOLUME / 10.0;
    let residual = |next: f64| 2.0 * mass * (next - lambda - h * velocity) / (h * h)
        + VOLUME * thermal_stress(0.5 * (lambda + next), s, scale);
    let (mut low, mut high) = (0.5, 1.7);
    assert!(residual(low) < 0.0 && residual(high) > 0.0);
    for _ in 0..100 {
        let middle = 0.5 * (low + high);
        if residual(middle) < 0.0 { low = middle; } else { high = middle; }
    }
    let next = 0.5 * (low + high);
    (next, 2.0 * (next - lambda) / h - velocity)
}

#[test]
fn g0_thermal_reference_volume_stress_tangent_and_native_assembly_match_closed_forms() {
    with_cx(|_, cx| {
        let card = thermal_card(0.01, 1.0, false);
        let material = prepared(&card);
        let s = 1.1;
        let lambda = 1.25;
        let f = [lambda, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0];
        let p = material.first_piola(&f).unwrap();
        let a = material.consistent_tangent(&f).unwrap();
        near(material.energy_density(&f), thermal_energy(lambda, s, 1.0), 2.0e-12);
        near(p[0], thermal_stress(lambda, s, 1.0), 2.0e-11);
        near(a[0][0], thermal_tangent(lambda, s, 1.0), 2.0e-10);

        let mut u = [0.0; 12]; u[3] = lambda - 1.0;
        let problem = problem(&material, &[], &[]);
        near(problem.stored_energy_j(&u, cx).unwrap(), VOLUME * thermal_energy(lambda, s, 1.0), 2.0e-12);
        near(problem.internal_force(&u, cx).unwrap()[3], VOLUME * thermal_stress(lambda, s, 1.0), 2.0e-11);
        near(problem.tangent_matrix(&u, cx).unwrap().get(3, 3), VOLUME * thermal_tangent(lambda, s, 1.0), 2.0e-10);
        let epsilon = 1.0e-6;
        for column in 0..9 {
            let mut plus = f; plus[column] += epsilon;
            let mut minus = f; minus[column] -= epsilon;
            let pp = material.first_piola(&plus).unwrap();
            let pm = material.first_piola(&minus).unwrap();
            near((material.energy_density(&plus) - material.energy_density(&minus)) / (2.0 * epsilon),
                p[column], 3.0e-8);
            for row in 0..9 {
                near((pp[row] - pm[row]) / (2.0 * epsilon), a[row][column], 3.0e-7);
            }
        }
    });
}

#[test]
fn g1_free_thermal_expansion_and_large_rotation_are_stress_free_with_original_mass() {
    with_cx(|_, cx| {
        let card = thermal_card(0.01, 1.0, false);
        let material = prepared(&card);
        let s = 1.1;
        // F=sR with R a ninety-degree rotation. Fe=R, not a strained body.
        let u = NODES.map(|x| [-s * x[1] - x[0], s * x[0] - x[1], (s - 1.0) * x[2]]);
        let body = FiniteSolidBody::new(
            problem(&material, &[], &[]), &u, &[[0.0; 3]; 4], 0.0, Default::default(), cx,
        ).unwrap();
        let (next, receipt) = body.step_candidate(
            0.01, &[], [0.0; 3], controls(1.0e-8), cx,
        ).unwrap();
        for (actual, expected) in next.displacement_m().iter().flatten().zip(u.iter().flatten()) {
            near(*actual, *expected, 1.0e-10);
        }
        near(receipt.strain_before_j, 0.0, 1.0e-12);
        near(receipt.strain_after_j, 0.0, 1.0e-12);
        near(receipt.kinetic_after_j, 0.0, 1.0e-12);
        near(next.reference_mass_kg(), 1.0, 2.0e-15);
        near(material.reference_density_kg_m3(), 6.0, 0.0);
        let hot_density = FiniteSolidMaterial::resolve(
            &card, 0, &point(310.0), MaterialPropertySelection::SingleClaimOnly,
        ).unwrap().reference_density_kg_m3();
        near(hot_density, 6.0 / s.powi(3), 1.0e-14);
        let thermal = material.thermal_reference().unwrap();
        assert_eq!(material.density().query(), *thermal.reference_context());
        assert_eq!(material.law().context(), thermal.current_context());
        let ThermalReference::Mean(thermal) = thermal else { panic!("mean source required") };
        assert_eq!(thermal.source().resolved().card_identity(), card.content_hash());
        assert_ne!(material.density().query(), *material.law().context());
    });
}

#[test]
fn g1_source_mean_expansion_drives_finite_motion_without_changing_reference_mass() {
    with_cx(|_, cx| {
        let fixed: Vec<_> = (0..12).filter(|dof| *dof != 3).map(|dof| (dof, 0.0)).collect();
        let mut positions = Vec::new();
        for (alpha, scale) in [(0.0, 1.0), (0.005, 1.0), (0.01, 1.0), (-0.005, 1.0), (0.01, 2.0)] {
            let card = thermal_card(alpha, scale, false);
            let material = prepared(&card);
            let mut body = FiniteSolidBody::new(
                problem(&material, &fixed, &[]), &[[0.0; 3]; 4], &[[0.0; 3]; 4],
                0.0, Default::default(), cx,
            ).unwrap();
            let s = 1.0 + 10.0 * alpha;
            let (mut lambda, mut velocity) = (1.0, 0.0);
            for _ in 0..4 {
                let expected = thermal_midpoint(lambda, velocity, s, scale, 0.01);
                let before = body.native_state().clone();
                let (next, receipt) = body.step_candidate(
                    0.01, &[], [0.0; 3], controls(1.0e-5), cx,
                ).unwrap();
                assert_eq!(body.native_state(), &before);
                near(next.native_state().q[0], expected.0 - 1.0, 2.0e-9);
                near(next.native_state().v[0], expected.1, 2.0e-8);
                near(next.reference_mass_kg(), 1.0, 2.0e-15);
                let mass = 6.0 * VOLUME / 10.0;
                let defect = VOLUME * (thermal_energy(expected.0, s, scale) - thermal_energy(lambda, s, scale))
                    + 0.5 * mass * (expected.1.powi(2) - velocity * velocity);
                near(receipt.integration_defect_j, defect, 3.0e-9);
                near(receipt.source_work_j, 0.0, 0.0);
                (lambda, velocity) = expected;
                body = next;
            }
            positions.push(body.native_state().q[0]);
        }
        near(positions[0], 0.0, 1.0e-12);
        assert!(positions[2] > positions[1] && positions[1] > 0.0);
        assert!(positions[3] < 0.0);
        assert!(positions[4] > positions[2]);
    });
}

#[test]
fn g1_each_element_retains_its_own_thermal_reference_on_one_native_motion_owner() {
    with_cx(|_, cx| {
        let card = thermal_card(0.01, 1.0, false);
        let cooler_reference = MeanThermalReference::resolve(
            &card, &PropertyKey::new(MEAN, INVERSE_TEMPERATURE_DIMS),
            &thermal_query(300.0, 305.0), "T_reference", "T_endpoint",
            MaterialPropertySelection::SingleClaimOnly,
        ).unwrap();
        let materials = [
            prepared(&card),
            FiniteSolidMaterial::resolve_mean_thermal(
                &card, 0, cooler_reference, MaterialPropertySelection::SingleClaimOnly,
            ).unwrap(),
        ];
        let mut nodes = NODES.to_vec();
        nodes.extend(NODES.map(|x| [x[0] + 3.0, x[1], x[2]]));
        let tets = [[0, 1, 2, 3], [4, 5, 6, 7]];
        let fixed: Vec<_> = (0..24).filter(|dof| *dof != 3 && *dof != 15)
            .map(|dof| (dof, 0.0)).collect();
        let problem = HyperTetProblem {
            nodes_m: &nodes, tetrahedra: &tets,
            material: HyperelasticLawField::PerElement(&materials),
            prescribed_m: &fixed, nodal_forces_n: &[], body_force_n_m3: [0.0; 3],
            budget: TetAssemblyBudget::standard(), settings: HyperTetSettings::default(),
        };
        let body = FiniteSolidBody::new(
            problem, &[[0.0; 3]; 8], &[[0.0; 3]; 8], 0.0, Default::default(), cx,
        ).unwrap();
        let (next, _) = body.step_candidate(0.01, &[], [0.0; 3], controls(1.0e-5), cx).unwrap();
        for (region, stretch) in [1.1, 1.05].into_iter().enumerate() {
            let expected = thermal_midpoint(1.0, 0.0, stretch, 1.0, 0.01);
            near(next.native_state().q[region], expected.0 - 1.0, 2.0e-9);
            near(next.native_state().v[region], expected.1, 2.0e-8);
            assert_eq!(materials[region].density().card_identity(), card.content_hash());
        }
        assert!(next.native_state().q[0] > next.native_state().q[1]);
        near(next.reference_mass_kg(), 2.0, 4.0e-15);
    });
}

#[test]
fn g4_thermal_reference_refuses_foreign_sources_dynamic_axes_and_wrong_temperature_roles() {
    let card = thermal_card(0.01, 1.0, false);
    let other = thermal_card(0.02, 1.0, false);
    assert!(FiniteSolidMaterial::resolve_mean_thermal(
        &other, 0, reference(&card), MaterialPropertySelection::SingleClaimOnly,
    ).is_err());
    let key = PropertyKey::new(MEAN, INVERSE_TEMPERATURE_DIMS);
    let resolve = |card: &MaterialCard, query: &QueryPoint| MeanThermalReference::resolve(
        card, &key, query, "T_reference", "T_endpoint", MaterialPropertySelection::SingleClaimOnly,
    );
    assert!(resolve(&card, &thermal_query(300.0, 320.0)).is_err());
    let wrong = thermal_query(300.0, 310.0)
        .with_quantity("T_endpoint", QuantitySpec::dimensional(Pressure::DIMS), 310.0).unwrap();
    assert!(resolve(&card, &wrong).is_err());
    let ambiguous = thermal_query(300.0, 310.0).with_quantity(
        "T", QuantitySpec::dimensional(Temperature::DIMS), 305.0,
    ).unwrap();
    assert!(resolve(&card, &ambiguous).is_err());
    let dynamic = thermal_card(0.01, 1.0, true);
    let stress = thermal_query(300.0, 310.0)
        .with_quantity("stress", QuantitySpec::dimensional(Pressure::DIMS), 20.0).unwrap();
    assert!(resolve(&dynamic, &stress).is_err());
    // A curve's independent stress axis need not appear in validity.bounds.
    // The source owner can evaluate it, but this held-temperature mechanical
    // adapter cannot freeze that evolving dependence.
    let mut curve_claims = ClaimSet::new();
    for (_, claim) in card.claims().claims_ordered() {
        let mut claim = (*claim).clone();
        if claim.key.name() == MEAN {
            assert!(claim.validity.bound("stress").is_none());
            claim.value = PropertyValue::Curve {
                abscissa: "stress".into(), abscissa_dims: Pressure::DIMS,
                knots: vec![(0.0, 0.005), (100.0, 0.015)],
                dims: INVERSE_TEMPERATURE_DIMS,
            };
            claim.interpolation = InterpolationPolicy::LinearInside;
        }
        curve_claims.insert_claim(claim).unwrap();
    }
    let curve_card = MaterialCard::assemble(
        card.id().clone(), curve_claims, card.models().to_vec(),
    ).unwrap();
    assert!(fs_material::state_point::resolve_mean_thermal_expansion_strain(
        &curve_card, &key, &stress, "T_reference", "T_endpoint",
        MaterialPropertySelection::SingleClaimOnly,
    ).is_ok());
    assert!(resolve(&curve_card, &stress).is_err());
    let impossible = thermal_card(-0.11, 1.0, false);
    assert!(resolve(&impossible, &thermal_query(300.0, 310.0)).is_err());
    assert!(resolve(&source_card(6.0, 1.0, false, false), &thermal_query(300.0, 310.0)).is_err());
    with_cx(|_, cx| {
        // Fe stays admitted and the very soft response is representable, but
        // the actual total-volume Jacobian overflows and must still refuse.
        let huge = thermal_card(2.0e101, 1.0e-305, false);
        let material = prepared(&huge);
        let stretch = 3.0 * material.thermal_reference().unwrap().principal_stretches()[0];
        let u: Vec<f64> = NODES.iter().flatten().map(|x| (stretch - 1.0) * x).collect();
        let p = problem(&material, &[], &[]);
        assert!(p.stored_energy_j(&u, cx).is_err());
        assert!(p.tangent_matrix(&u, cx).is_err());
    });
}
