//! Original source controls against independent gray equations and latent storage.
use super::*;
use fs_conduction::radiation::pointwise::material::adjoint::MaterialEmissivityGradient;
use fs_conduction::radiation::pointwise::material::MaterialRadiationSmoothnessError;
use fs_conduction::transient::enthalpy::EnthalpyError;
use fs_conduction::{ConductionError, LinearConfig};
use fs_qty::Temperature;

const DT: f64 = 0.4;
const SEED: [f64; 8] = [0.7, -0.2, 0.4, 0.3, -0.6, 0.9, 0.2, -0.1];

fn linear() -> LinearConfig {
    LinearConfig { tolerance: 1e-12, max_iterations: 128, restart: 8 }
}

fn source_card(name: &str, axis: Option<&str>, knots: &[(f64, f64)]) -> MaterialCard {
    let (value, interpolation) = if let Some(axis) = axis {
        (PropertyValue::Curve { abscissa: axis.into(),
            abscissa_dims: if axis == "T" { Temperature::DIMS } else { Dims::NONE },
            knots: knots.to_vec(), dims: EMISSIVITY_DIMS }, InterpolationPolicy::LinearInside)
    } else {
        (PropertyValue::Scalar { value: knots[0].1, dims: EMISSIVITY_DIMS },
            InterpolationPolicy::ConstantWithinValidity)
    };
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS), value,
        validity: ValidityDomain::unconstrained().with("T", 300.0, 500.0)
            .with("melt", 0.0, 1.0).with("lot", 7.0, 7.0),
        uncertainty: UncertaintyModel::Unstated, interpolation, observations: Vec::new(),
        provenance: Provenance { source: "independent material enclosure control oracle".into(),
            license: "internal-test-use".into(), artifact: None },
    }).unwrap();
    MaterialCard::assemble(MaterialStateId { chemistry: name.into(), phase: "chart".into(),
        process: "radiation control fixture".into(), revision: 0 }, claims, Vec::new()).unwrap()
}

fn open(cx: &asupersync::Cx<'_>, mesh: &ConductionMesh, cards: &[MaterialCard; 2],
    context: &QueryPoint) -> OpenDiffuseEnclosure
{
    OpenDiffuseEnclosure::new(cx, surfaces(mesh, cards, context),
        vec![vec![0.0, 0.4], vec![0.4, 0.0]],
        ViewFactorEvidence::Analytic { geometry: "reciprocal control oracle".into() },
        ViewFactorTolerance::default()).unwrap()
}

// Direct two-patch radiosity algebra; no production factorization or quadrature.
fn heat(epsilon: [f64; 2], mutual: f64, bath: f64) -> [f64; 2] {
    let black = [SIGMA * 400.0_f64.powi(4), SIGMA * 350.0_f64.powi(4)];
    let bath = SIGMA * bath.powi(4);
    let rhs: [f64; 2] = std::array::from_fn(|i|
        epsilon[i] * black[i] + (1.0 - epsilon[i]) * (1.0 - mutual) * bath);
    let a = (1.0 - epsilon[0]) * mutual;
    let b = (1.0 - epsilon[1]) * mutual;
    let j = [(rhs[0] + a * rhs[1]) / (1.0 - a * b),
        (rhs[1] + b * rhs[0]) / (1.0 - a * b)];
    std::array::from_fn(|i| AREA * epsilon[i] * (black[i] - mutual * j[1 - i] - (1.0 - mutual) * bath))
}

// Controls: four source ordinates, bath K, two histories, densities and W/m³.
// On each latent trace three identical ρ/24 kg nodes share Q/3. The fourth
// node sees only volume heat. This fixed point is contractive for DT=0.4 s.
fn endpoint(controls: &[f64; 11], mutual: f64) -> [f64; 8] {
    let old = [controls[5], controls[6]];
    let rho = [controls[7], controls[8]];
    let source = [controls[9], controls[10]];
    let mut face = old;
    for _ in 0..128 {
        let epsilon = std::array::from_fn(|i| controls[2 * i]
            + (controls[2 * i + 1] - controls[2 * i]) * face[i] / LATENT);
        let q = heat(epsilon, mutual, controls[4]);
        face = std::array::from_fn(|i| old[i] + DT * (source[i] - 8.0 * q[i]) / rho[i]);
    }
    std::array::from_fn(|i| if i % 4 == 3 { old[i / 4] + DT * source[i / 4] / rho[i / 4] }
        else { face[i / 4] })
}

fn five_point(mut objective: impl FnMut(f64) -> f64, delta: f64) -> f64 {
    (objective(-2.0 * delta) - 8.0 * objective(-delta)
        + 8.0 * objective(delta) - objective(2.0 * delta)) / (12.0 * delta)
}

#[test]
fn reflected_latent_adjoint_matches_source_history_density_and_heat_oracles() {
    let mesh = mesh(); let boundary = boundary(&mesh); let context = context();
    let controls = [0.2, 0.8, 0.25, 0.7, 300.0, 520.0, 470.0, 24.0, 30.0, 12.0, 7.0];
    let cards = [source_card("left", Some("melt"), &[(0.0, controls[0]), (1.0, controls[1])]),
        source_card("right", Some("melt"), &[(0.0, controls[2]), (1.0, controls[3])])];
    let material = ConductivityModel::isotropic_declared(2.0).unwrap();
    let source = ScalarField::Nodal((0..8).map(|i| controls[9 + i / 4]).collect());
    let problem = ConductionProblem { mesh: &mesh, boundary: &boundary, material: &material,
        source: &source, element_materials: None };
    let enclosure = closed(&mesh, &cards, &context);
    for datum in [0.0, 4096.0] { for mutual in [1.0, 0.4] {
        let curves = [phase_curve(&cards[0], 400.0, datum), phase_curve(&cards[1], 350.0, datum)];
        let declarations = bindings(&cards, &curves, &context);
        let old: [f64; 8] = std::array::from_fn(|i| datum + controls[5 + i / 4]);
        with_cx(|cx| {
            let open = open(cx, &mesh, &cards, &context);
            let base = if mutual == 1.0 {
                PointwiseEnclosureRadiation::new(cx, &mesh, &boundary, &enclosure,
                    PointwiseEnclosureTolerance::default()).unwrap()
            } else { PointwiseEnclosureRadiation::new_open(cx, &mesh, &boundary, &open, controls[4],
                PointwiseEnclosureTolerance::default()).unwrap() };
            let radiation = MaterialPointwiseEnclosureRadiation::new(cx, &base, &declarations).unwrap();
            let engine = HeterogeneousEnthalpyBackwardEuler::new(cx, &mesh,
                &[ReferenceEnthalpyMaterial { curve: &curves[0], reference_density_kg_m3: controls[7] },
                    ReferenceEnthalpyMaterial { curve: &curves[1], reference_density_kg_m3: controls[8] }],
                &[0, 1], budget()).unwrap();
            let accepted = engine.advance_with_material_enclosure_radiation(cx, problem, None,
                &old, DT, config(), &radiation).unwrap();
            for (actual, expected) in accepted.specific_enthalpy_j_kg.iter().zip(endpoint(&controls, mutual)) {
                close(actual - datum, expected, 2e-8);
            }
            let entries = radiation.adjoint_workspace_entries().unwrap();
            assert!(engine.linearize_accepted_with_material_enclosure_radiation(cx, problem, None,
                &old, DT, config(), &radiation, accepted.clone(), 1e-8, entries - 1).is_err());
            let mut changed = accepted.clone();
            changed.pointwise_enclosure_radiation.as_mut().unwrap().patches[0].material_mean_emissivity = Some(0.9);
            assert!(engine.linearize_accepted_with_material_enclosure_radiation(cx, problem, None,
                &old, DT, config(), &radiation, changed, 1e-8, entries).is_err());
            let bound = engine.linearize_accepted_with_material_enclosure_radiation(cx, problem, None,
                &old, DT, config(), &radiation, accepted, 1e-8, entries).unwrap();
            assert_eq!(bound.transport().temperature_pullback(cx, &SEED).unwrap(), vec![0.0; 8]);
            assert!(bound.transport().pointwise_parameter_pullback(cx, &SEED).is_err());
            let gradient = bound.pullback_with_material_enclosure_controls(cx, &SEED, linear()).unwrap();
            let mut derivatives = Vec::new();
            for (i, patch) in gradient.radiation.patches.iter().enumerate() {
                assert_eq!(patch.material_identity, cards[i].content_hash());
                let MaterialEmissivityGradient::Curve { axis, knot_abscissae, knot_value_derivatives } = &patch.emissivity
                    else { panic!("original curve controls required") };
                assert_eq!(axis, "melt"); assert_eq!(knot_abscissae, &[0.0, 1.0]);
                derivatives.extend(knot_value_derivatives);
            }
            assert_eq!(gradient.radiation.reservoir_temperature.is_none(), mutual == 1.0);
            derivatives.push(gradient.radiation.reservoir_temperature.unwrap_or(0.0));
            for i in 0..2 { derivatives.push(gradient.thermal.transport.previous_specific_enthalpy[4 * i..4 * i + 4].iter().sum()); }
            derivatives.extend(&gradient.thermal.material_reference_density);
            for i in 0..2 { derivatives.push(gradient.thermal.transport.source_density[4 * i..4 * i + 4].iter().sum()); }
            for (control, derivative) in derivatives.into_iter().enumerate() {
                let delta = if control < 4 { 1e-4 } else { 0.01 };
                let expected = five_point(|offset| {
                    let mut perturbed = controls; perturbed[control] += offset;
                    endpoint(&perturbed, mutual).iter().zip(SEED).map(|(h, seed)| h * seed).sum()
                }, delta);
                close(derivative, expected, 3e-7 * expected.abs().max(1.0));
            }
            if mutual == 1.0 {
                let closed = bound.transport().material_enclosure_parameter_pullback(cx, &[1.0; 8]).unwrap();
                for patch in closed.patches {
                    let MaterialEmissivityGradient::Curve { knot_value_derivatives, .. } = patch.emissivity else { unreachable!() };
                    for value in knot_value_derivatives { close(value, 0.0, 2e-10); }
                }
            }
            with_cancelled_cx(|cancelled| assert!(bound
                .pullback_with_material_enclosure_controls(cancelled, &SEED, linear()).is_err()));
        });
    } }
}

#[test]
fn spatial_scalar_and_temperature_ordinate_controls_use_original_nonuniform_reflection() {
    let mesh = mesh(); let boundary = boundary(&mesh); let context = context();
    let cards = [source_card("left", None, &[(0.0, 0.55)]),
        source_card("right", Some("T"), &[(300.0, 0.2), (350.0, 0.5), (500.0, 0.75)])];
    let curves = [phase_curve(&cards[0], 400.0, 0.0), phase_curve(&cards[1], 350.0, 0.0)];
    let declarations = bindings(&cards, &curves, &context);
    let material = ConductivityModel::isotropic_declared(2.0).unwrap(); let source = ScalarField::Uniform(0.0);
    let problem = ConductionProblem { mesh: &mesh, boundary: &boundary, material: &material,
        source: &source, element_materials: None };
    let old = [-30.0, 450.0, 1100.0, 550.0, -50.0, 300.0, 1080.0, 500.0];
    with_cx(|cx| {
        let open = open(cx, &mesh, &cards, &context);
        let base = PointwiseEnclosureRadiation::new_open(cx, &mesh, &boundary, &open, 300.0,
            PointwiseEnclosureTolerance::default()).unwrap();
        let radiation = MaterialPointwiseEnclosureRadiation::new(cx, &base, &declarations).unwrap();
        let engine = HeterogeneousEnthalpyBackwardEuler::new(cx, &mesh,
            &[ReferenceEnthalpyMaterial { curve: &curves[0], reference_density_kg_m3: RHO },
                ReferenceEnthalpyMaterial { curve: &curves[1], reference_density_kg_m3: RHO }],
            &[0, 1], budget()).unwrap();
        let accepted = engine.advance_with_material_enclosure_radiation(cx, problem, None,
            &old, 0.01, config(), &radiation).unwrap();
        let t: [[f64; 3]; 2] = std::array::from_fn(|body|
            std::array::from_fn(|node| accepted.temperature[4 * body + node]));
        // The original temperature source has a corner at 350 K, but its
        // composed h derivative is smooth at this interior latent node.
        assert_eq!(t[1][1], 350.0);
        let bound = engine.linearize_accepted_with_material_enclosure_radiation(cx, problem, None,
            &old, 0.01, config(), &radiation, accepted, 1e-8,
            radiation.adjoint_workspace_entries().unwrap()).unwrap();
        let gradient = bound.transport().material_enclosure_parameter_pullback(cx, &SEED).unwrap();
        let MaterialEmissivityGradient::Scalar { derivative } = &gradient.patches[0].emissivity else { panic!("scalar source") };
        let mut actual = vec![*derivative];
        let MaterialEmissivityGradient::Curve { knot_value_derivatives, .. } = &gradient.patches[1].emissivity else { panic!("temperature source") };
        actual.extend(knot_value_derivatives); actual.push(gradient.reservoir_temperature.unwrap());
        let objective = |values: [f64; 5]| {
            let epsilon = [ [values[0]; 3], t[1].map(|temperature| {
                let (lo, hi, a, b) = if temperature <= 350.0 {
                    (300.0, 350.0, values[1], values[2])
                } else { (350.0, 500.0, values[2], values[3]) };
                a + (b - a) * (temperature - lo) / (hi - lo)
            }) ];
            let emitted: [[f64; 3]; 2] = std::array::from_fn(|i| emitted_load(t[i], epsilon[i]));
            let means = epsilon.map(|e| e.iter().sum::<f64>() / 3.0);
            let bath = SIGMA * values[4].powi(4);
            let rhs: [f64; 2] = std::array::from_fn(|i| emitted[i].iter().sum::<f64>() / AREA
                + (1.0 - means[i]) * 0.6 * bath);
            let a = (1.0 - means[0]) * 0.4; let b = (1.0 - means[1]) * 0.4;
            let j = [(rhs[0] + a * rhs[1]) / (1.0 - a * b), (rhs[1] + b * rhs[0]) / (1.0 - a * b)];
            let mut total = 0.0;
            for i in 0..2 { for node in 0..3 {
                let weight = AREA / 12.0 * (epsilon[i].iter().sum::<f64>() + epsilon[i][node]);
                total += SEED[4 * i + node] * (emitted[i][node] - weight * (0.4 * j[1 - i] + 0.6 * bath));
            } }
            total
        };
        for (index, actual) in actual.into_iter().enumerate() {
            let expected = five_point(|offset| {
                let mut values = [0.55, 0.2, 0.5, 0.75, 300.0]; values[index] += offset;
                objective(values)
            }, if index < 4 { 1e-4 } else { 0.02 });
            close(actual, expected, 2e-7 * expected.abs().max(1.0));
        }
    });
}

#[test]
fn material_enclosure_chart_and_moving_source_corners_refuse_typed_adjoints() {
    let mesh = mesh(); let boundary = boundary(&mesh); let context = context();
    let cards = [source_card("left", Some("melt"), &[(0.0, 0.2), (0.5, 0.6), (1.0, 0.8)]),
        source_card("right", Some("melt"), &[(0.0, 0.2), (0.5, 0.6), (1.0, 0.8)])];
    let curves = [phase_curve(&cards[0], 400.0, 0.0), phase_curve(&cards[1], 400.0, 0.0)];
    let declarations = bindings(&cards, &curves, &context); let enclosure = closed(&mesh, &cards, &context);
    let material = ConductivityModel::isotropic_declared(2.0).unwrap(); let source = ScalarField::Uniform(0.0);
    let problem = ConductionProblem { mesh: &mesh, boundary: &boundary, material: &material,
        source: &source, element_materials: None };
    for h in [0.0, 500.0] { with_cx(|cx| {
        let base = PointwiseEnclosureRadiation::new(cx, &mesh, &boundary, &enclosure,
            PointwiseEnclosureTolerance::default()).unwrap();
        let radiation = MaterialPointwiseEnclosureRadiation::new(cx, &base, &declarations).unwrap();
        let engine = HeterogeneousEnthalpyBackwardEuler::new(cx, &mesh,
            &[ReferenceEnthalpyMaterial { curve: &curves[0], reference_density_kg_m3: RHO },
                ReferenceEnthalpyMaterial { curve: &curves[1], reference_density_kg_m3: RHO }],
            &[0, 1], budget()).unwrap();
        let old = [h; 8];
        let accepted = engine.advance_with_material_enclosure_radiation(cx, problem, None,
            &old, DT, config(), &radiation).unwrap();
        assert_eq!(accepted.specific_enthalpy_j_kg, old);
        let error = engine.linearize_accepted_with_material_enclosure_radiation(cx, problem, None,
            &old, DT, config(), &radiation, accepted, 1e-8,
            radiation.adjoint_workspace_entries().unwrap()).unwrap_err();
        let EnthalpyAdjointError::Enthalpy(EnthalpyError::Conduction(
            ConductionError::RadiationMaterialSmoothness(error))) = error
            else { panic!("typed source/chart error required: {error:?}") };
        if h == 0.0 { assert!(matches!(error, MaterialRadiationSmoothnessError::ChartKink { .. })); }
        else { assert!(matches!(error, MaterialRadiationSmoothnessError::SourceKink { .. })); }
    }); }
}
