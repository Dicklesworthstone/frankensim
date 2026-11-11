use super::*;
use fs_conduction::{LinearConfig, ThermalBc};

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
            let original = endpoint.conduction.clone();
            let step = engine.linearize_temperature_contact_endpoint(
                cx, problem, &contacts, &histories[index], dts[index], config(),
                endpoint.conduction, coupling().residual_tolerance_j, 1,
            ).unwrap();
            assert_eq!(step.primal().specific_enthalpy_j_kg, original.specific_enthalpy_j_kg);
            assert_eq!(step.primal().temperature, original.temperature);
            assert_eq!(step.primal().liquid_mass_fraction, original.liquid_mass_fraction);
            assert_eq!(step.primal().newton.iterations, original.newton.iterations);
            assert_eq!(step.primal().newton.history, original.newton.history);
            let seed = carry.take().unwrap_or_else(|| {
                let mut temperature_seed = vec![0.0; 8];
                temperature_seed[0] = 1.0;
                step.transport().temperature_pullback(cx, &temperature_seed).unwrap()
            });
            let gradient = step.pullback(cx, &seed, linear()).unwrap();
            let load = step.transport().nodal_load_pullback(cx, &gradient.transport.adjoint).unwrap();
            let controls = endpoint.physical_contacts.frozen_interfaces().resistance_scale_pullback(
                cx, &step.primal().temperature, &load, 1,
            ).unwrap();
            assert_eq!(controls.len(), 1);
            assert_eq!(controls[0].interface, "joint");
            resistance += controls[0].derivative;
            carry = Some(gradient.transport.previous_specific_enthalpy);
        } else {
            let contact = contacts.linearize_temperature(
                cx, problem.boundary, &endpoint.conduction.temperature, 1,
            ).unwrap();
            let step = engine.linearize_accepted(
                cx, problem, Some(endpoint.applied_contacts.frozen_interfaces()),
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
            resistance += endpoint.applied_contacts.frozen_interfaces()
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

fn absolute_kelvin() -> fs_qty::QuantitySpec {
    use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
    fs_qty::QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}

fn two_sided_scaled_card(scale: f64) -> InterfaceSystemCard {
    let dims = fs_conduction::AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS;
    let mut claims = ClaimSet::new();
    let samples = [[256., 256.], [384., 256.], [256., 384.], [384., 384.]]
        .into_iter().map(|p| {
            (p, scale * (0.5 + (p[0] - 320.) / 1024. + 3. * (p[1] - 320.) / 1024.))
        }).collect();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(fs_conduction::AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY, dims),
        value: PropertyValue::Surface2D {
            axes: ["TA".into(), "TB".into()],
            axis_dims: [fs_conduction::TEMPERATURE_DIMS; 2],
            samples, triangles: vec![[0, 1, 2], [1, 3, 2]], dims,
        },
        validity: ValidityDomain::unconstrained().with_quantity("TA", absolute_kelvin(), 256., 384.)
            .with_quantity("TB", absolute_kelvin(), 256., 384.).with("normal_pressure", 1e6, 1e6),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::LinearInside,
        observations: vec![],
        provenance: Provenance {
            source: format!("manufactured ordered contact law, scale {scale}; not measurements"),
            license: "internal-test-use".into(), artifact: None,
        },
    }).unwrap();
    let side = |name: &str| SurfaceSpec {
        material: MaterialStateId {
            chemistry: name.into(), phase: "solid".into(),
            process: "fixture".into(), revision: 0,
        },
        texture_frame: "fixed reference normal".into(),
    };
    InterfaceSystemCard::assemble(side("A"), side("B"), SystemContext {
        medium: "dry".into(), third_body: None,
        environment: "vacuum".into(), history: "unaged".into(),
    }, claims, vec![]).unwrap()
}

fn two_sided_law<'c>(mesh: &ConductionMesh, card: &'c InterfaceSystemCard) -> TemperatureContact<'c> {
    let kelvin = absolute_kelvin();
    let point = QueryPoint::new().with_quantity("TA", kelvin, 320.).unwrap()
        .with_quantity("TB", kelvin, 320.).unwrap().with("normal_pressure", 1e6).unwrap();
    let mut pair = ThermalInterfaces::coincident_face_pairs(mesh).unwrap()[0];
    if mesh.boundary()[pair.side_a].vertices[0] >= 4 {
        std::mem::swap(&mut pair.side_a, &mut pair.side_b);
    }
    TemperatureContact::new_two_sided("joint", vec![pair], card, &point,
        ["TA", "TB"], SelectionPolicy::SingleClaimOnly).unwrap()
}

#[test]
fn g3_two_sided_enthalpy_adjoint_matches_endpoint_resistance_history_and_source_directions() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh();
        let bc = ThermalBoundaryBuilder::new(&mesh).adiabatic_remainder().finish().unwrap();
        let card = two_sided_scaled_card(1.0);
        let contacts = TemperatureContacts::new(cx, &mesh, &bc, vec![],
            vec![two_sided_law(&mesh, &card)], 1).unwrap();
        let cold = chart(false);
        let hot = chart(true);
        let engine = HeterogeneousEnthalpyBackwardEuler::new(cx, &mesh, &[
            ReferenceEnthalpyMaterial { curve: &cold, reference_density_kg_m3: 1.0 },
            ReferenceEnthalpyMaterial { curve: &hot, reference_density_kg_m3: 2.0 },
        ], &[0, 1], EnthalpyBudget { max_vertices: 8, max_elements: 2 }).unwrap();
        let material = ConductivityModel::isotropic(
            ConductivityTable::declared_curve(vec![(280.0, 1.6), (380.0, 3.6)]).unwrap(),
        );
        let source_values = [0., 0., 0., 0., 60., 60., 60., 60.];
        let source = ScalarField::Nodal(source_values.to_vec());
        let problem = ConductionProblem {
            mesh: &mesh, boundary: &bc, material: &material,
            element_materials: None, source: &source,
        };
        let old = [230., 240., 260., 150., 430., 450., 470., 450.];
        let dt = 0.35;
        let target = ContactIterationConfig { residual_tolerance_j: 1e-10, ..coupling() };
        let accepted = contacts.advance_enthalpy(
            cx, &engine, problem, &old, dt, config(), target,
        ).unwrap();
        // Both contact traces have nonzero dT/dh, while the cold apex retains
        // latent energy that cannot be reconstructed from its temperature.
        assert!(accepted.conduction.specific_enthalpy_j_kg[3] > 100.
            && accepted.conduction.specific_enthalpy_j_kg[3] < 200.);
        close(accepted.conduction.temperature[3], 310., 1e-12);
        assert!(accepted.physical_contacts.states()[0].side_a_temperature_k > 310.);
        assert!(accepted.physical_contacts.states()[0].side_b_temperature_k < 350.);
        close(accepted.physical.stored_energy_change_j, 10. * dt, 1e-8);
        let step = engine.linearize_temperature_contact_endpoint(
            cx, problem, &contacts, &old, dt, config(), accepted.conduction,
            target.residual_tolerance_j, 1,
        ).unwrap();
        let temperature_seed = [0.7, -0.2, 0.4, 0.5, -0.3, 0.1, 0.2, -0.1];
        let seed = step.transport().temperature_pullback(cx, &temperature_seed).unwrap();
        for (i, &value) in seed.iter().enumerate() {
            close(value, if i == 3 { 0.0 } else { 0.1 * temperature_seed[i] }, 1e-14);
        }
        let gradient = step.pullback(cx, &seed, linear()).unwrap();
        let load = step.transport().nodal_load_pullback(cx, &gradient.transport.adjoint).unwrap();
        let controls = accepted.physical_contacts.frozen_interfaces().resistance_scale_pullback(
            cx, &step.primal().temperature, &load, 1,
        ).unwrap();
        assert_eq!(controls.len(), 1);
        assert_eq!(controls[0].interface, "joint");
        let history_direction = [7., -3., 4., 11., -5., 8., -6., 2.];
        let source_direction = [5., -2., 8., 11., -7., 4., -3., 6.];
        let derivatives = [
            controls[0].derivative,
            gradient.transport.previous_specific_enthalpy.iter()
                .zip(history_direction).map(|(bar, direction)| bar * direction).sum(),
            gradient.transport.source_density.iter()
                .zip(source_direction).map(|(bar, direction)| bar * direction).sum(),
        ];
        assert!(derivatives.iter().all(|value| value.abs() > 0.01));

        // Every difference re-solves the actual nonlinear endpoint with fresh
        // source data or physical inputs; no frozen contact operator is used.
        let objective = |scale: f64, history: &[f64], density: &[f64]| {
            let card = two_sided_scaled_card(scale);
            let contacts = TemperatureContacts::new(cx, &mesh, &bc, vec![],
                vec![two_sided_law(&mesh, &card)], 1).unwrap();
            let source = ScalarField::Nodal(density.to_vec());
            let endpoint = contacts.advance_enthalpy(cx, &engine,
                ConductionProblem { source: &source, ..problem },
                history, dt, config(), target).unwrap();
            endpoint.conduction.temperature.iter().zip(temperature_seed)
                .map(|(temperature, weight)| temperature * weight).sum::<f64>()
        };
        for epsilon in [2e-4, 1e-4] {
            let plus_h: [f64; 8] = std::array::from_fn(|i| old[i] + epsilon * history_direction[i]);
            let minus_h: [f64; 8] = std::array::from_fn(|i| old[i] - epsilon * history_direction[i]);
            let plus_q: [f64; 8] = std::array::from_fn(|i| source_values[i] + epsilon * source_direction[i]);
            let minus_q: [f64; 8] = std::array::from_fn(|i| source_values[i] - epsilon * source_direction[i]);
            let expected = [
                (objective(1. + epsilon, &old, &source_values)
                    - objective(1. - epsilon, &old, &source_values)) / (2. * epsilon),
                (objective(1., &plus_h, &source_values)
                    - objective(1., &minus_h, &source_values)) / (2. * epsilon),
                (objective(1., &old, &plus_q)
                    - objective(1., &old, &minus_q)) / (2. * epsilon),
            ];
            for (actual, reference) in derivatives.into_iter().zip(expected) {
                close(actual, reference, 2e-5 * reference.abs().max(1.0));
            }
        }
    });
}

#[test]
fn g3_live_enthalpy_endpoint_prescribed_columns_and_direct_goal_match_fresh_bath_solves() {
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh();
        let cold = chart(false);
        let hot = chart(true);
        let engine = HeterogeneousEnthalpyBackwardEuler::new(cx, &mesh, &[
            ReferenceEnthalpyMaterial { curve: &cold, reference_density_kg_m3: 1.0 },
            ReferenceEnthalpyMaterial { curve: &hot, reference_density_kg_m3: 2.0 },
        ], &[0, 1], EnthalpyBudget { max_vertices: 8, max_elements: 2 }).unwrap();
        let material = ConductivityModel::isotropic(
            ConductivityTable::declared_curve(vec![(280.0, 1.6), (380.0, 3.6)]).unwrap(),
        );
        let card = two_sided_scaled_card(1.0);
        let source = ScalarField::Nodal(vec![0., 0., 0., 0., 60., 60., 60., 60.]);
        let old = [230., 240., 260., 150., 430., 450., 470., 450.];
        let dt = 0.35;
        for (side, bath_k) in [(0_usize, 300.0), (1, 340.0)] {
            let boundary = |temperature| ThermalBoundaryBuilder::new(&mesh).region(
                "bath", |face| {
                    let mut vertices = face.vertices;
                    vertices.sort_unstable();
                    vertices == [4 * side, 4 * side + 1, 4 * side + 3].map(|v| v as u32)
                }, ThermalBc::dirichlet(temperature).unwrap(),
            ).unwrap().adiabatic_remainder().finish().unwrap();
            let bc = boundary(bath_k);
            let contacts = TemperatureContacts::new(cx, &mesh, &bc, vec![],
                vec![two_sided_law(&mesh, &card)], 1).unwrap();
            let problem = ConductionProblem { mesh: &mesh, boundary: &bc,
                material: &material, element_materials: None, source: &source };
            let accepted = contacts.advance_enthalpy_newton(
                cx, &engine, problem, &old, dt, coupled_newton(),
            ).unwrap();
            let step = engine.linearize_temperature_contact_endpoint(
                cx, problem, &contacts, &old, dt, config(), accepted.conduction,
                coupling().residual_tolerance_j, 1,
            ).unwrap();
            let mut weights = [0.0; 8];
            weights[4 * (1 - side) + 2] = 1.0;
            weights[4 * side] = 0.4;
            let seed = step.transport().temperature_pullback(cx, &weights).unwrap();
            let gradient = step.pullback(cx, &seed, linear()).unwrap();
            let actual: f64 = gradient.transport.prescribed_temperature.iter().sum();
            assert!(actual > 0.4 + 1e-4, "the bath must also influence the opposite body");
            let load = step.transport().nodal_load_pullback(cx, &gradient.transport.adjoint).unwrap();
            for &(vertex, _) in bc.dirichlet() {
                assert_eq!(gradient.transport.previous_specific_enthalpy[vertex], 0.0);
                assert_eq!(load[vertex], 0.0);
            }
            let objective = |temperature| {
                let bc = boundary(temperature);
                let contacts = TemperatureContacts::new(cx, &mesh, &bc, vec![],
                    vec![two_sided_law(&mesh, &card)], 1).unwrap();
                let solved = contacts.advance_enthalpy_newton(cx, &engine,
                    ConductionProblem { boundary: &bc, ..problem }, &old, dt,
                    coupled_newton()).unwrap();
                solved.conduction.temperature.iter().zip(weights)
                    .map(|(temperature, weight)| temperature * weight).sum::<f64>()
            };
            let epsilon = 1e-3;
            let expected = (objective(bath_k + epsilon) - objective(bath_k - epsilon))
                / (2.0 * epsilon);
            close(actual, expected, 2e-6 * expected.abs().max(1.0));
        }
    });
}

#[test]
fn g4_live_enthalpy_endpoint_keeps_original_gates_and_strict_source_chart_admission() {
    use fs_conduction::transient::enthalpy::adjoint::EnthalpyAdjointError;
    context(&CancelGate::new_clock_free(), |cx| {
        let mesh = mesh();
        let bc = ThermalBoundaryBuilder::new(&mesh).adiabatic_remainder().finish().unwrap();
        let curve = chart(false);
        let engine = HeterogeneousEnthalpyBackwardEuler::new(cx, &mesh, &[
            ReferenceEnthalpyMaterial { curve: &curve, reference_density_kg_m3: 1.0 },
            ReferenceEnthalpyMaterial { curve: &curve, reference_density_kg_m3: 2.0 },
        ], &[0, 1], EnthalpyBudget { max_vertices: 8, max_elements: 2 }).unwrap();
        let material = ConductivityModel::isotropic_declared(2.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let problem = ConductionProblem { mesh: &mesh, boundary: &bc,
            material: &material, element_materials: None, source: &source };
        let original_card = card();
        let contacts = TemperatureContacts::new(cx, &mesh, &bc, vec![],
            vec![law(&mesh, &original_card)], 1).unwrap();
        let old = [50.0; 8];
        let accepted = contacts.advance_enthalpy_newton(
            cx, &engine, problem, &old, 0.2, coupled_newton(),
        ).unwrap().conduction;
        let mut changed = accepted.clone();
        changed.specific_enthalpy_j_kg[0] += 1.0;
        // A loose inner Newton configuration cannot replace the original gate.
        let mut loose = config();
        loose.newton.absolute_tolerance = 10.0;
        assert!(matches!(engine.linearize_temperature_contact_endpoint(
            cx, problem, &contacts, &old, 0.2, loose, changed.clone(), 1e-10, 1,
        ), Err(EnthalpyAdjointError::PrimalResidual { .. })));
        assert!(matches!(engine.linearize_temperature_contact_endpoint(
            cx, problem, &contacts, &old, 0.2, loose, changed, 10.0, 1,
        ), Err(EnthalpyAdjointError::Enthalpy(EnthalpyError::EnergyBalance { .. }))));
        assert!(engine.linearize_temperature_contact_endpoint(
            cx, problem, &contacts, &old, 0.2, config(), accepted.clone(), 1e-10, 0,
        ).is_err());
        let gate = CancelGate::new_clock_free();
        gate.request();
        context(&gate, |cancelled| assert!(matches!(engine.linearize_temperature_contact_endpoint(
            cancelled, problem, &contacts, &old, 0.2, config(), accepted.clone(), 1e-10, 1,
        ), Err(EnthalpyAdjointError::Enthalpy(EnthalpyError::Conduction(
            fs_conduction::ConductionError::Cancelled { .. },
        ))))));

        // Both are valid stationary forward endpoints. Neither a source kink
        // nor a chart kink may acquire a classical adjoint, even for zero heat.
        for (source_card, old) in [
            (card_with_curve(vec![(280., 0.1), (accepted.temperature[0], 0.3), (380., 0.7)]), [50.; 8]),
            (card(), [100.; 8]),
        ] {
            let contacts = TemperatureContacts::new(cx, &mesh, &bc, vec![],
                vec![law(&mesh, &source_card)], 1).unwrap();
            let accepted = contacts.advance_enthalpy_newton(
                cx, &engine, problem, &old, 0.2, coupled_newton(),
            ).unwrap().conduction;
            assert!(engine.linearize_temperature_contact_endpoint(
                cx, problem, &contacts, &old, 0.2, config(), accepted, 1e-10, 1,
            ).is_err());
        }
    });
}
