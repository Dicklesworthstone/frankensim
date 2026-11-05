//! G1/G3: actual spatial Robin reservoirs, without a mean-bath substitution.
use super::*;

fn face_boundary(
    mesh: &ConductionMesh,
    coefficients: [f64; 3],
    baths: [f64; 3],
) -> ThermalBoundary {
    let selected = mesh.boundary()[0].vertices;
    let mut htc = vec![1.0; 4];
    let mut reference = vec![300.0; 4];
    for i in 0..3 {
        htc[selected[i] as usize] = coefficients[i];
        reference[selected[i] as usize] = baths[i];
    }
    ThermalBoundaryBuilder::new(mesh)
        .region("affine reservoir", |face| face.vertices == selected,
            ThermalBc::Robin {
                htc: ScalarField::Nodal(htc),
                t_ref: ScalarField::Nodal(reference),
            })
        .unwrap().adiabatic_remainder().finish().unwrap()
}

fn face_temperatures(mesh: &ConductionMesh, values: [f64; 3]) -> [f64; 4] {
    let mut temperatures = [350.0; 4];
    for (i, &vertex) in mesh.boundary()[0].vertices.iter().enumerate() {
        temperatures[vertex as usize] = values[i];
    }
    temperatures
}

// Independent positive Duffy/Gauss integration of the original rational
// integrand. No isotherm subdivision or logarithmic moment formula is reused.
fn face_reference(
    area: f64,
    htc: [f64; 3],
    temperature: [f64; 3],
    bath: [f64; 3],
) -> (f64, f64, f64) {
    let (nodes, weights) = fs_feec::gauss_legendre(32);
    let delta = std::array::from_fn::<_, 3, _>(|i| temperature[i] - bath[i]);
    let mut bath_entropy = 0.0;
    let mut production = 0.0;
    let mut heat = 0.0;
    for (&x, &wx) in nodes.iter().zip(&weights) {
        let u = 0.5 * (x + 1.0);
        for (&y, &wy) in nodes.iter().zip(&weights) {
            let v = 0.5 * (y + 1.0);
            let n = [u, (1.0 - u) * v, (1.0 - u) * (1.0 - v)];
            let measure = 0.5 * area * wx * wy * (1.0 - u);
            let h: f64 = (0..3).map(|i| n[i] * htc[i]).sum();
            let d: f64 = (0..3).map(|i| n[i] * delta[i]).sum();
            let b: f64 = (0..3).map(|i| n[i] * bath[i]).sum();
            let q = measure * h * d;
            bath_entropy -= q / b;
            heat += q;
            if temperature[0] == temperature[1] && temperature[1] == temperature[2] {
                // The continuous and receiving-node reciprocal T coincide.
                // This form independently retains near-equilibrium production.
                production += q * (d / b) / temperature[0];
            } else {
                let reciprocal: f64 = (0..3).map(|i| n[i] / temperature[i]).sum();
                production += q * (1.0 / b - reciprocal);
            }
        }
    }
    (bath_entropy, production, heat)
}

#[test]
fn spatial_bath_entropy_matches_independent_integration_and_vertex_permutations() {
    let mesh = tetrahedron();
    let material = ConductivityModel::isotropic_declared(2.0).unwrap();
    let source = ScalarField::Uniform(0.0);
    let chart = sensible_curve();
    let htc = [2.0, 3.0, 7.0];
    let t = [320.0, 450.0, 500.0];
    let dt = 0.3;
    with_cx(|cx| {
        let model = EnthalpyBackwardEuler::uniform(
            cx, &mesh, &chart, REFERENCE_DENSITY, budget(),
        ).unwrap();
        for baths in [
            [300.0, 400.0, 600.0],
            [300.0, 300.0, 600.0],
            [300.0, 600.0, 600.0],
            [300.0, 300.0 + 1e-10, 300.0 + 2e-10],
        ] {
            let (expected_bath, expected_production, expected_heat) =
                face_reference(mesh.boundary()[0].area, htc, t, baths);
            for order in [[0, 1, 2], [0, 2, 1], [1, 0, 2], [1, 2, 0], [2, 0, 1], [2, 1, 0]] {
                let coefficients = order.map(|i| htc[i]);
                let references = order.map(|i| baths[i]);
                let temperatures = face_temperatures(&mesh, order.map(|i| t[i]));
                let boundary = face_boundary(&mesh, coefficients, references);
                let h = sensible_h(&temperatures);
                let saved = bits(&h);
                let (endpoint, entropy) = model.evaluate_endpoint_with_entropy(
                    cx, problem(&mesh, &boundary, &material, &source),
                    None, &h, &h, dt, config(),
                ).unwrap();
                near(entropy.robin_entropy_input_j_k, dt * expected_bath, 2e-12);
                near(entropy.robin_generation_j_k, dt * expected_production, 2e-12);
                near(endpoint.convective_out_w, expected_heat, 2e-10);
                near(entropy.balance_residual_j_k, 0.0, 2e-12);
                assert_eq!(entropy.temporal_defect_j_k, 0.0);
                assert!(!entropy.has_unclassified_inputs);
                assert_eq!(bits(&h), saved);
            }
        }
    });
}

#[test]
fn almost_constant_baths_retain_small_positive_entropy_and_exact_local_equilibrium() {
    let mesh = tetrahedron();
    let material = ConductivityModel::isotropic_declared(2.0).unwrap();
    let source = ScalarField::Uniform(0.0);
    let chart = sensible_curve();
    let coefficients = [2.0, 3.0, 7.0];
    let epsilon = 2.0_f64.powi(-20);
    let baths = [300.0 - epsilon, 300.0, 300.0 + 0.5 * epsilon];
    let boundary = face_boundary(&mesh, coefficients, baths);
    with_cx(|cx| {
        let model = EnthalpyBackwardEuler::uniform(
            cx, &mesh, &chart, REFERENCE_DENSITY, budget(),
        ).unwrap();
        let h = sensible_h(&[300.0; 4]);
        let (_, entropy) = model.evaluate_endpoint_with_entropy(
            cx, problem(&mesh, &boundary, &material, &source),
            None, &h, &h, 1.0, config(),
        ).unwrap();
        let (expected_bath, production, _) =
            face_reference(mesh.boundary()[0].area, coefficients, [300.0; 3], baths);
        assert!(production > 0.0 && production < 1e-16);
        near(entropy.robin_generation_j_k, production, 1e-11 * production);
        near(entropy.robin_entropy_input_j_k, expected_bath, 1e-20);
        // The separate closure also sees ordinary sparse-operator roundoff.
        near(entropy.balance_residual_j_k, 0.0, 1e-12);

        // T(x)=B(x) pointwise gives zero Robin heat, even on a varying bath.
        let temperatures = face_temperatures(&mesh, baths);
        let h = sensible_h(&temperatures);
        let actual = mesh.boundary()[0].vertices.map(|v| {
            chart.values_at_specific_enthalpy(h[v as usize]).unwrap().temperature_k
        });
        let equilibrium_boundary = face_boundary(&mesh, coefficients, actual);
        let (endpoint, equilibrium) = model.evaluate_endpoint_with_entropy(
            cx, problem(&mesh, &equilibrium_boundary, &material, &source),
            None, &h, &h, 1.0, config(),
        ).unwrap();
        near(endpoint.convective_out_w, 0.0, 1e-12);
        assert_eq!(equilibrium.robin_entropy_input_j_k, 0.0);
        assert_eq!(equilibrium.robin_generation_j_k, 0.0);
    });
}

#[test]
fn large_bath_contrast_uses_a_logarithmic_integral_without_refining_the_face() {
    let mesh = tetrahedron();
    let area = mesh.boundary()[0].area;
    let material = ConductivityModel::isotropic_declared(2.0).unwrap();
    let source = ScalarField::Uniform(0.0);
    let chart = sensible_curve();
    let low = 1.0;
    let high = 1e12;
    let body = 400.0;
    let htc = 3.0;
    let boundary = face_boundary(&mesh, [htc; 3], [low, low, high]);
    // Independent single-coordinate integral for B=low+(high-low)*N_2:
    // integral 1/B dA = 2A*(ln(high/low)-(1-low/high))/(high*(1-low/high)^2).
    let d = 1.0 - low / high;
    let reciprocal = 2.0 * area / high * ((high / low).ln() - d) / (d * d);
    let bath = htc * (area - body * reciprocal);
    let generation = htc * (body * reciprocal - 2.0 * area
        + area * (2.0 * low + high) / (3.0 * body));
    with_cx(|cx| {
        let model = EnthalpyBackwardEuler::uniform(
            cx, &mesh, &chart, REFERENCE_DENSITY, budget(),
        ).unwrap();
        let h = sensible_h(&[body; 4]);
        let (_, entropy) = model.evaluate_endpoint_with_entropy(
            cx, problem(&mesh, &boundary, &material, &source),
            None, &h, &h, 0.25, config(),
        ).unwrap();
        near(entropy.robin_entropy_input_j_k, 0.25 * bath, 2e-13);
        near(entropy.robin_generation_j_k, 0.25 * generation, 2e-13 * generation);
        let roundoff_scale = entropy.robin_generation_j_k.abs()
            + entropy.algebraic_residual_j_k.abs();
        near(entropy.balance_residual_j_k, 0.0, 128.0 * f64::EPSILON * roundoff_scale);
    });
}

#[test]
fn native_be_and_sdirk_spatial_baths_keep_energy_and_actual_entropy_quadrature() {
    use fs_conduction::transient::enthalpy::heterogeneous::EnthalpyRadiation;
    use fs_conduction::transient::enthalpy::sdirk::{EnthalpySdirk2Input, SDIRK2_GAMMA};
    let mesh = tetrahedron();
    let coefficients = [2.0, 4.0, 7.0];
    let baths = [300.0, 350.0, 400.0];
    let boundary = face_boundary(&mesh, coefficients, baths);
    let material = ConductivityModel::isotropic_declared(2.0).unwrap();
    let source = ScalarField::Uniform(0.0);
    let chart = sensible_curve();
    let old = sensible_h(&[360.0; 4]);
    let saved = bits(&old);
    let dt = 0.01;
    with_cx(|cx| {
        let model = EnthalpyBackwardEuler::uniform(
            cx, &mesh, &chart, REFERENCE_DENSITY, budget(),
        ).unwrap();
        let p = problem(&mesh, &boundary, &material, &source);
        let original = model.advance(cx, p, None, &old, dt, config()).unwrap();
        let (step, entropy) = model.advance_with_entropy(
            cx, p, None, &old, dt, config(), 1e-10,
        ).unwrap();
        assert_eq!(bits(&step.specific_enthalpy_j_kg), bits(&original.specific_enthalpy_j_kg));
        assert_eq!(step.newton, original.newton);
        let temperatures = mesh.boundary()[0].vertices.map(|v| step.temperature[v as usize]);
        let (bath, generation, heat) =
            face_reference(mesh.boundary()[0].area, coefficients, temperatures, baths);
        near(entropy.robin_entropy_input_j_k, dt * bath, 2e-12);
        near(entropy.robin_generation_j_k, dt * generation, 2e-12);
        near(step.stored_energy_change_j, -dt * heat, 1e-9);
        assert!(entropy.robin_generation_j_k > 0.0);
        near(entropy.algebraic_residual_j_k, 0.0, 1e-10);
        near(entropy.balance_residual_j_k, 0.0, 1e-12);

        let input = EnthalpySdirk2Input {
            problem: p, interfaces: None, radiation: EnthalpyRadiation::None, contacts: None,
        };
        let original = model.advance_sdirk2(cx, [input; 2], &old, dt, config()).unwrap();
        let (second_order, entropy) = model.advance_sdirk2_with_entropy(
            cx, [input; 2], &old, dt, config(), 1e-10,
        ).unwrap();
        assert_eq!(bits(&second_order.first_stage.specific_enthalpy_j_kg),
            bits(&original.first_stage.specific_enthalpy_j_kg));
        assert_eq!(bits(&second_order.second_stage.specific_enthalpy_j_kg),
            bits(&original.second_stage.specific_enthalpy_j_kg));
        let mut expected_bath = 0.0;
        let mut expected_production = 0.0;
        let mut expected_heat = 0.0;
        for (stage, weight) in [
            (&second_order.first_stage, (1.0 - SDIRK2_GAMMA) * dt),
            (&second_order.second_stage, SDIRK2_GAMMA * dt),
        ] {
            let temperatures = mesh.boundary()[0].vertices.map(|v| stage.temperature[v as usize]);
            let (bath, production, heat) =
                face_reference(mesh.boundary()[0].area, coefficients, temperatures, baths);
            expected_bath += weight * bath;
            expected_production += weight * production;
            expected_heat += weight * heat;
        }
        near(entropy.robin_entropy_input_j_k, expected_bath, 2e-12);
        near(entropy.robin_generation_j_k, expected_production, 2e-12);
        near(second_order.stored_energy_change_j, -expected_heat, 1e-9);
        near(entropy.algebraic_residual_j_k, 0.0, 1e-10);
        near(entropy.balance_residual_j_k, 0.0, 1e-12);
        with_cancelled_cx(|cancelled| {
            assert!(matches!(
                model.advance_with_entropy(cancelled, p, None, &old, dt, config(), 1e-10),
                Err(EnthalpyEntropyError::Physical(EnthalpyError::Conduction(
                    ConductionError::Cancelled { .. }
                )))
            ));
        });
        assert_eq!(bits(&old), saved);
    });
}
