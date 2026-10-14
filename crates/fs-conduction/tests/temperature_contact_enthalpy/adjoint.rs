use super::*;
use fs_conduction::LinearConfig;

fn scaled_card(scale: f64) -> InterfaceSystemCard {
    let dims = fs_conduction::AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS;
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(fs_conduction::AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY, dims),
        value: PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: fs_conduction::TEMPERATURE_DIMS,
            knots: vec![(280.0, 0.1 * scale), (380.0, 0.7 * scale)],
            dims,
        },
        validity: ValidityDomain::unconstrained()
            .with("T", 280.0, 380.0)
            .with("normal_pressure", 1e6, 1e6),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::LinearInside,
        observations: vec![],
        provenance: Provenance {
            source: format!("synthetic scaled contact law {scale}, not measurements"),
            license: "internal-test-use".into(),
            artifact: None,
        },
    }).unwrap();
    let side = |name: &str| SurfaceSpec {
        material: MaterialStateId {
            chemistry: name.into(), phase: "solid".into(),
            process: "fixture".into(), revision: 0,
        },
        texture_frame: "fixed reference normal".into(),
    };
    InterfaceSystemCard::assemble(
        side("A"), side("B"),
        SystemContext {
            medium: "dry".into(), third_body: None,
            environment: "vacuum".into(), history: "unaged".into(),
        },
        claims, vec![],
    ).unwrap()
}

fn linear() -> LinearConfig {
    LinearConfig { tolerance: 1e-10, max_iterations: 512, restart: 8 }
}

fn final_temperature(cx: &Cx<'_>, scale: f64) -> f64 {
    let mesh = mesh();
    let bc = ThermalBoundaryBuilder::new(&mesh).adiabatic_remainder().finish().unwrap();
    let card = scaled_card(scale);
    let contacts = TemperatureContacts::new(
        cx, &mesh, &bc, vec![], vec![law(&mesh, &card)], 1,
    ).unwrap();
    let cold = chart(false);
    let hot = chart(true);
    let engine = HeterogeneousEnthalpyBackwardEuler::new(
        cx, &mesh,
        &[
            ReferenceEnthalpyMaterial { curve: &cold, reference_density_kg_m3: 1.0 },
            ReferenceEnthalpyMaterial { curve: &hot, reference_density_kg_m3: 2.0 },
        ],
        &[0, 1],
        EnthalpyBudget { max_vertices: 8, max_elements: 2 },
    ).unwrap();
    let material = ConductivityModel::isotropic(
        ConductivityTable::declared_curve(vec![(280.0, 1.6), (380.0, 3.6)]).unwrap(),
    );
    let source = ScalarField::Nodal(vec![0., 0., 0., 0., 60., 60., 60., 60.]);
    let problem = ConductionProblem {
        mesh: &mesh, boundary: &bc, material: &material,
        element_materials: None, source: &source,
    };
    let mut old = vec![80., 80., 80., 80., 450., 450., 450., 450.];
    for dt in [0.15, 0.35, 0.5] {
        let next = contacts.advance_enthalpy(
            cx, &engine, problem, &old, dt, config(), coupling(),
        ).unwrap();
        old = next.conduction.specific_enthalpy_j_kg;
    }
    engine.phase_curve_at(0).unwrap()
        .state_at_specific_enthalpy(old[0]).unwrap().temperature_k()
}

fn trajectory_gradient(cx: &Cx<'_>, complete_contact_feedback: bool) -> (f64, Vec<f64>) {
    let mesh = mesh();
    let bc = ThermalBoundaryBuilder::new(&mesh).adiabatic_remainder().finish().unwrap();
    let card = scaled_card(1.0);
    let contacts = TemperatureContacts::new(
        cx, &mesh, &bc, vec![], vec![law(&mesh, &card)], 1,
    ).unwrap();
    let cold = chart(false);
    let hot = chart(true);
    let engine = HeterogeneousEnthalpyBackwardEuler::new(
        cx, &mesh,
        &[
            ReferenceEnthalpyMaterial { curve: &cold, reference_density_kg_m3: 1.0 },
            ReferenceEnthalpyMaterial { curve: &hot, reference_density_kg_m3: 2.0 },
        ],
        &[0, 1],
        EnthalpyBudget { max_vertices: 8, max_elements: 2 },
    ).unwrap();
    let material = ConductivityModel::isotropic(
        ConductivityTable::declared_curve(vec![(280.0, 1.6), (380.0, 3.6)]).unwrap(),
    );
    let source = ScalarField::Nodal(vec![0., 0., 0., 0., 60., 60., 60., 60.]);
    let problem = ConductionProblem {
        mesh: &mesh, boundary: &bc, material: &material,
        element_materials: None, source: &source,
    };
    let dts = [0.15, 0.35, 0.5];
    let mut histories = vec![vec![80., 80., 80., 80., 450., 450., 450., 450.]];
    let mut accepted = Vec::new();
    for &dt in &dts {
        let next = contacts.advance_enthalpy(
            cx, &engine, problem, histories.last().unwrap(), dt, config(), coupling(),
        ).unwrap();
        histories.push(next.conduction.specific_enthalpy_j_kg.clone());
        accepted.push(next);
    }

    let mut carry: Option<Vec<f64>> = None;
    let mut resistance = 0.0;
    for index in (0..dts.len()).rev() {
        let endpoint = accepted.pop().unwrap();
        if complete_contact_feedback {
            let step = contacts.linearize_accepted_enthalpy(
                cx, &engine, problem, &histories[index], dts[index], config(),
                coupling().residual_tolerance_j, endpoint, 1,
            ).unwrap();
            let seed = carry.take().unwrap_or_else(|| {
                let mut temperature_seed = vec![0.0; 8];
                temperature_seed[0] = 1.0;
                step.temperature_pullback(cx, &temperature_seed).unwrap()
            });
            let gradient = step.pullback_with_resistance_controls(
                cx, &seed, linear(), 1,
            ).unwrap();
            assert_eq!(gradient.resistance_multipliers.len(), 1);
            assert_eq!(gradient.resistance_multipliers[0].interface, "joint");
            resistance += gradient.resistance_multipliers[0].derivative;
            carry = Some(gradient.storage.transport.previous_specific_enthalpy);
        } else {
            let contact = contacts.linearize_temperature(
                cx, problem.boundary, &endpoint.conduction.temperature, 1,
            ).unwrap();
            let step = engine.linearize_accepted(
                cx, problem, Some(contact.physical_contacts().frozen_interfaces()),
                &histories[index], dts[index], config(), endpoint.conduction,
            ).unwrap();
            let seed = carry.take().unwrap_or_else(|| {
                let mut temperature_seed = vec![0.0; 8];
                temperature_seed[0] = 1.0;
                step.transport().temperature_pullback(cx, &temperature_seed).unwrap()
            });
            let gradient = step.pullback(cx, &seed, linear()).unwrap();
            let mut load_dual = gradient.transport.adjoint.clone();
            for value in &mut load_dual { *value *= dts[index]; }
            resistance += contact.physical_contacts().frozen_interfaces()
                .matching_resistance_scale_pullback(
                    cx, contact.temperature(), &load_dual, 1,
                ).unwrap()[0].derivative;
            carry = Some(gradient.transport.previous_specific_enthalpy);
        }
    }
    (resistance, carry.unwrap())
}

#[test]
fn g3_enthalpy_contact_adjoint_matches_fresh_resistance_laws_through_latent_history() {
    context(&CancelGate::new_clock_free(), |cx| {
        let (actual, initial_bar) = trajectory_gradient(cx, true);
        assert!(actual.abs() > 1e-5);
        assert!(initial_bar[..4].iter().any(|value| value.abs() > 1e-8),
            "later temperature must retain energetic dependence through the earlier latent plateau");

        for epsilon in [2e-4, 1e-4] {
            let expected = (
                final_temperature(cx, 1.0 + epsilon)
                    - final_temperature(cx, 1.0 - epsilon)
            ) / (2.0 * epsilon);
            close(actual, expected, 5e-5 * expected.abs().max(1.0));
        }

        let (frozen, _) = trajectory_gradient(cx, false);
        assert!((frozen - actual).abs() > 0.01 * actual.abs().max(1e-4),
            "freezing R(T) must not masquerade as the complete trajectory derivative");
    });
}

#[test]
fn g0_enthalpy_contact_tangent_has_true_transpose_and_zero_seed_controls() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh();
        let bc = ThermalBoundaryBuilder::new(&mesh).adiabatic_remainder().finish().unwrap();
        let card = scaled_card(1.0);
        let contacts = TemperatureContacts::new(
            cx, &mesh, &bc, vec![], vec![law(&mesh, &card)], 1,
        ).unwrap();
        let cold = chart(false);
        let hot = chart(true);
        let engine = HeterogeneousEnthalpyBackwardEuler::new(
            cx, &mesh,
            &[
                ReferenceEnthalpyMaterial { curve: &cold, reference_density_kg_m3: 1.0 },
                ReferenceEnthalpyMaterial { curve: &hot, reference_density_kg_m3: 2.0 },
            ],
            &[0, 1],
            EnthalpyBudget { max_vertices: 8, max_elements: 2 },
        ).unwrap();
        let material = ConductivityModel::isotropic(
            ConductivityTable::declared_curve(vec![(280.0, 1.6), (380.0, 3.6)]).unwrap(),
        );
        let source = ScalarField::Nodal(vec![0., 0., 0., 0., 60., 60., 60., 60.]);
        let problem = ConductionProblem {
            mesh: &mesh, boundary: &bc, material: &material,
            element_materials: None, source: &source,
        };
        let old = [80., 80., 80., 80., 450., 450., 450., 450.];
        let accepted = contacts.advance_enthalpy(
            cx, &engine, problem, &old, 0.15, config(), coupling(),
        ).unwrap();
        let step = contacts.linearize_accepted_enthalpy(
            cx, &engine, problem, &old, 0.15, config(),
            coupling().residual_tolerance_j, accepted, 1,
        ).unwrap();
        let x = [0.3, -0.4, 0.7, 0.2, -0.6, 0.1, 0.5, -0.2];
        let z = [-0.2, 0.9, 0.1, -0.5, 0.4, -0.3, 0.8, 0.6];
        let jx = step.transport().transport().apply_jacobian(cx, &x).unwrap();
        let jtz = step.transport().transport().apply_jacobian_transpose(cx, &z).unwrap();
        let lhs: f64 = jx.iter().zip(z).map(|(a, b)| a * b).sum();
        let rhs: f64 = x.iter().zip(jtz).map(|(a, b)| a * b).sum();
        close(lhs, rhs, 2e-12 * lhs.abs().max(rhs.abs()).max(1.0));

        let zero = step.pullback_with_resistance_controls(
            cx, &[0.0; 8], linear(), 1,
        ).unwrap();
        assert!(zero.storage.transport.previous_specific_enthalpy.iter().all(|&x| x == 0.0));
        assert!(zero.storage.transport.source_density.iter().all(|&x| x == 0.0));
        assert_eq!(zero.resistance_multipliers[0].derivative, 0.0);
    });
}
