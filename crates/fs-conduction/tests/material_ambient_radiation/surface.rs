//! G1/G3: joint temperature/fraction sources and original sample controls.
use super::*;
use fs_conduction::radiation::pointwise::material::MaterialRadiationSmoothnessError;
use fs_conduction::radiation::pointwise::material::adjoint::MaterialEmissivityGradient;
use fs_conduction::transient::enthalpy::EnthalpyError;
use fs_conduction::{ConductionError, LinearConfig};

const COORDINATES: [[f64; 2]; 4] = [[300.0, 0.0], [500.0, 0.0], [300.0, 1.0], [500.0, 1.0]];
const TRIANGLES: [[u32; 3]; 2] = [[1, 3, 2], [0, 1, 2]];
const VALUES: [f64; 4] = [0.1875, 0.3125, 0.6875, 0.8125];

fn surface_card(reverse: bool, samples: Vec<([f64; 2], f64)>, triangles: Vec<[u32; 3]>) -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS),
        value: PropertyValue::Surface2D {
            axes: if reverse { ["melt".into(), "T".into()] } else { ["T".into(), "melt".into()] },
            axis_dims: if reverse { [Dims::NONE, Temperature::DIMS] } else { [Temperature::DIMS, Dims::NONE] },
            samples: samples.into_iter().map(|(point, value)|
                (if reverse { [point[1], point[0]] } else { point }, value)).collect(),
            triangles, dims: EMISSIVITY_DIMS,
        },
        validity: ValidityDomain::unconstrained().with("T", 300.0, 500.0)
            .with("melt", 0.0, 1.0).with("lot", 7.0, 7.0),
        interpolation: InterpolationPolicy::LinearInside, uncertainty: UncertaintyModel::Unstated,
        observations: Vec::new(), provenance: Provenance { source: "manufactured joint thermal emissivity".into(),
            license: "internal-test-use".into(), artifact: None },
    }).unwrap();
    MaterialCard::assemble(MaterialStateId { chemistry: "synthetic".into(), phase: "chart".into(),
        process: "joint thermal source".into(), revision: 0 }, claims, Vec::new()).unwrap()
}

pub(super) fn rectangle(reverse: bool) -> MaterialCard {
    surface_card(reverse, COORDINATES.into_iter().zip(VALUES).collect(), TRIANGLES.to_vec())
}

fn inverse(rhs: [f64; 3], beta: f64) -> [f64; 3] {
    let a = 1.0 + beta;
    let correction = beta * rhs.iter().sum::<f64>() / (a * (a + 3.0 * beta));
    rhs.map(|value| value / a - correction)
}

#[test]
fn latent_energy_full_action_and_original_surface_samples_match_independent_moments() {
    let mesh = mesh(); let boundary = boundary(&mesh); let context = context();
    let material = ConductivityModel::isotropic_declared(2.0).unwrap();
    let source = ScalarField::Uniform(0.0);
    let problem = ConductionProblem { mesh: &mesh, boundary: &boundary, material: &material,
        source: &source, element_materials: None };
    let old = [250.0, 500.0, 750.0, 600.0];
    let seed = [0.7, -0.3, 0.2, 0.1]; let direction = [0.3, -0.7, 0.9, -0.2]; let dt = 0.2;
    for reverse in [false, true] { for bath in [300.0_f64, 450.0] {
        let card = rectangle(reverse); let curve = curve(&card, 0.0);
        // At T=400 this source is epsilon=1/4+f/2. Each face node has mass1.
        let c = STEFAN_BOLTZMANN_W_M2_K4 * (400.0_f64.powi(4) - bath.powi(4));
        let beta = dt * c * 0.5 / LATENT * AREA / 12.0;
        let expected = inverse(std::array::from_fn(|i| old[i] - dt * c * 0.25 * AREA / 3.0), beta);
        let lambda = inverse([seed[0], seed[1], seed[2]], beta);
        with_cx(|cx| {
            let radiation = MaterialPointwiseAmbientRadiation::new(cx, &mesh, &boundary,
                &[patch(&card, &curve, &context, bath)]).unwrap();
            let engine = HeterogeneousEnthalpyBackwardEuler::new(cx, &mesh,
                &[ReferenceEnthalpyMaterial { curve: &curve, reference_density_kg_m3: RHO }],
                &[0], budget()).unwrap();
            let accepted = engine.advance_with_material_ambient_radiation(cx, problem, None,
                &old, dt, config(), &radiation).unwrap();
            for i in 0..3 { close(accepted.specific_enthalpy_j_kg[i], expected[i], 3e-9); }
            assert_eq!(accepted.temperature, vec![400.0; 4]);
            close(accepted.energy_residual_j, 0.0, config().energy_tolerance_j);
            let entries = radiation.adjoint_workspace_entries().unwrap();
            assert!(engine.linearize_accepted_with_material_ambient_radiation(cx, problem, None,
                &old, dt, config(), &radiation, accepted.clone(), 1e-8, entries - 1).is_err());
            let bound = engine.linearize_accepted_with_material_ambient_radiation(cx, problem, None,
                &old, dt, config(), &radiation, accepted, 1e-8, entries).unwrap();
            assert_eq!(bound.transport().temperature_pullback(cx, &seed).unwrap(), vec![0.0; 4]);
            let action = bound.transport().apply_jacobian(cx, &direction).unwrap();
            for i in 0..3 { close(action[i], direction[i]
                + beta * (direction[i] + direction[..3].iter().sum::<f64>()), 4e-13); }
            let gradient = bound.pullback_with_material_ambient_controls(cx, &seed,
                LinearConfig { tolerance: 1e-12, max_iterations: 64, restart: 4 }).unwrap();
            assert_eq!(gradient.radiation.patches[0].material_identity, card.content_hash());
            let MaterialEmissivityGradient::Surface { axes, sample_coordinates, triangles,
                sample_value_derivatives } = &gradient.radiation.patches[0].emissivity
                else { panic!("original surface controls required") };
            assert_eq!(axes, &if reverse { ["melt".to_string(), "T".to_string()] }
                else { ["T".to_string(), "melt".to_string()] });
            assert_eq!(sample_coordinates, &COORDINATES.map(|p| if reverse { [p[1], p[0]] } else { p }).to_vec());
            assert_eq!(triangles, &TRIANGLES.to_vec());
            let mut independent = [0.0; 4];
            for i in 0..3 { for j in 0..3 {
                let f = expected[j] / LATENT;
                let basis = if f <= 0.5 { [0.5 - f, 0.5, f, 0.0] }
                    else { [0.0, 1.0 - f, 0.5, f - 0.5] };
                let moment = AREA / 12.0 * if i == j { 2.0 } else { 1.0 };
                for p in 0..4 { independent[p] -= dt * c * lambda[i] * moment * basis[p]; }
            } }
            for (actual, expected) in sample_value_derivatives.iter().zip(independent) {
                close(*actual, expected, 3e-11);
            }
        });
    } }
}

#[test]
fn joint_temperature_and_fraction_tangent_matches_independent_face_loads() {
    let mesh = mesh(); let boundary = boundary(&mesh); let context = context();
    for reverse in [false, true] {
        let card = rectangle(reverse);
        let curve = EquilibriumEnthalpyPhaseCurve::try_new(card.content_hash(), [
            (-100.0, 300.0, 0.0), (0.0, 360.0, 0.0), (1000.0, 440.0, 1.0), (1200.0, 500.0, 1.0),
        ].map(|(h, t, f)| EnthalpyPhaseKnot { specific_enthalpy_j_kg: h, temperature_k: t,
            liquid_mass_fraction: f, bulk_density_kg_m3: RHO }).to_vec()).unwrap();
        let h = [200.0, 400.0, 800.0, 500.0]; let direction = [0.7, -0.9, 1.3, 2.0];
        let exact = |h: [f64; 4]| {
            let t = [0, 1, 2].map(|i| 360.0 + 0.08 * h[i]);
            let epsilon = [0, 1, 2].map(|i| 0.1875 + 0.000625 * (t[i] - 300.0) + 0.5 * h[i] / LATENT);
            exact_face_load(t, epsilon, 300.0)
        };
        with_cx(|cx| {
            let radiation = MaterialPointwiseAmbientRadiation::new(cx, &mesh, &boundary,
                &[patch(&card, &curve, &context, 300.0)]).unwrap();
            let assembly = radiation.assemble(cx, &h).unwrap();
            for (actual, expected) in assembly.outward_load_w[..3].iter().zip(exact(h)) { close(*actual, expected, 2e-10); }
            let mut action = [0.0; 4];
            assembly.enthalpy_tangent_apply(cx, &direction, &mut action).unwrap();
            let plus = exact(std::array::from_fn(|i| h[i] + 1e-3 * direction[i]));
            let minus = exact(std::array::from_fn(|i| h[i] - 1e-3 * direction[i]));
            for i in 0..3 { close(action[i], (plus[i] - minus[i]) / 2e-3, 2e-7); }
        });
    }
}

#[test]
fn strict_surface_neighborhoods_follow_only_the_locally_active_chart_axes() {
    let mesh = mesh(); let boundary = boundary(&mesh); let context = context();
    let material = ConductivityModel::isotropic_declared(2.0).unwrap(); let source = ScalarField::Uniform(0.0);
    let problem = ConductionProblem { mesh: &mesh, boundary: &boundary, material: &material,
        source: &source, element_materials: None };
    let held_temperature_crease = surface_card(false,
        vec![([400.0, 0.0], 0.25), ([400.0, 1.0], 0.75), ([300.0, 0.5], 0.375), ([500.0, 0.5], 0.75)],
        vec![[0, 2, 1], [0, 1, 3]]);
    let moving_fraction_crease = surface_card(false,
        vec![([300.0, 0.5], 0.25), ([500.0, 0.5], 0.5), ([400.0, 0.0], 0.25), ([400.0, 1.0], 0.75)],
        vec![[0, 2, 1], [0, 1, 3]]);
    for (card, h, bath, admitted) in [
        (held_temperature_crease, 500.0, 300.0, true),
        (moving_fraction_crease, 500.0, 400.0, false),
        (rectangle(false), -50.0, 350.0, true),
        (rectangle(true), 1100.0, 450.0, true),
    ] {
        let curve = curve(&card, 0.0);
        with_cx(|cx| {
            let radiation = MaterialPointwiseAmbientRadiation::new(cx, &mesh, &boundary,
                &[patch(&card, &curve, &context, bath)]).unwrap();
            let engine = EnthalpyBackwardEuler::uniform(cx, &mesh, &curve, RHO, budget()).unwrap();
            let old = [h; 4];
            let accepted = engine.advance_with_material_ambient_radiation(cx, problem, None,
                &old, 0.1, config(), &radiation).unwrap();
            let bound = engine.linearize_accepted_with_material_ambient_radiation(cx, problem, None,
                &old, 0.1, config(), &radiation, accepted, 1e-8, radiation.adjoint_workspace_entries().unwrap());
            if admitted { assert!(bound.is_ok(), "{bound:?}"); }
            else { assert!(matches!(bound, Err(EnthalpyAdjointError::Enthalpy(EnthalpyError::Conduction(
                ConductionError::RadiationMaterialSmoothness(MaterialRadiationSmoothnessError::SourceSurface { .. })
            ))))); }
        });
    }
}
