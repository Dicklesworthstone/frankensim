//! Full source tensors drive the original nonlinear heat and stress owners.
use super::*;
use fs_couple::solid_enthalpy::spatial::SpatialConductivitySource;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

const KEYS: [&str; 3] = ["conductivity_longitudinal", "conductivity_radial", "conductivity_tangential"];

fn pressure_quantity() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static))
}
fn context() -> QueryPoint {
    point(300.).with_quantity("P", pressure_quantity(), 101325.).unwrap()
}
fn rotation() -> [[f64; 3]; 3] {
    let q = 0.5_f64.sqrt();
    [[q, -q, 0.], [q, q, 0.], [0., 0., 1.]]
}
fn binding(r: [[f64; 3]; 3]) -> SpatialConductivitySource {
    SpatialConductivitySource::OrientedOrthotropic {
        principal_properties: KEYS.map(|key| PropertyKey::new(key, CONDUCTIVITY_DIMS)),
        principal_selections: [ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly); 3],
        principal_to_mesh: r,
        principal_frame: ContentHash([73; 32]),
        mesh_frame: solid_config().mesh_frame,
    }
}
fn bind_explicit<'a, 'm: 'a, 'c: 'a>(
    native: &'a HeterogeneousEnthalpyBackwardEuler<'m, 'c>,
    source: SpatialMaterialSource<'a>,
    binding: SpatialConductivitySource,
    cx: &Cx<'_>,
) -> Result<SpatialSolidEnthalpy<'a, 'm, 'c>, SpatialSolidEnthalpyError> {
    SpatialSolidEnthalpy::try_new_with_conductivity(
        native, vec![source], vec![binding], TetAssemblyBudget::standard(),
        solid_config().mesh_frame, 0.1, cx,
    )
}

// Manufactured principal laws at 300, 400, 500 K. No scalar k exists on this
// card. The exact fixed pressure condition must survive every source query.
fn card(values: [[f64; 3]; 3]) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let (base, _) = fixture(1., 1000., 1e-5, 10., 500., 750.);
    let mut claims = ClaimSet::new();
    for (_, claim) in base.claims().claims_ordered() {
        if claim.key.name() != "thermal_conductivity" {
            claims.insert_claim(claim.clone()).unwrap();
        }
    }
    for (name, values) in KEYS.into_iter().zip(values) {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, CONDUCTIVITY_DIMS),
            value: PropertyValue::Curve {
                abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                knots: [300., 400., 500.].into_iter().zip(values).collect(),
                dims: CONDUCTIVITY_DIMS,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), 300., 500.)
                .with_quantity("P", pressure_quantity(), 100000., 110000.),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured orthotropic P1 oracle".into(),
                license: "CC0-1.0".into(), artifact: None,
            },
        }).unwrap();
    }
    let card = MaterialCard::assemble(MaterialStateId {
        chemistry: "manufactured directional conductor".into(),
        phase: "solid".into(), process: "source substitution oracle".into(), revision: 0,
    }, claims, vec![]).unwrap();
    let caloric = ResolvedHeatCapacityCurve::try_new(&card, HeatCapacitySourceQuery {
        heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
        density: &PropertyKey::new("density", Density::DIMS), temperature_axis: "T",
        lower: &point(300.), upper: &point(500.), selection: SelectionPolicy::SingleClaimOnly,
    }, SolidLiquidPhase::Solid, 0., 1e-8).unwrap();
    (card, caloric)
}

#[test]
fn oriented_source_tensor_and_derivative_drive_nonlinear_heat_energy_and_stress() {
    with_cx(|cx| {
        let (card, caloric) = card([[2., 3., 4.], [6., 9., 12.], [1., 1.5, 2.]]);
        let mesh = mesh();
        let native = native(&mesh, caloric.curve(), 6., cx);
        let context = context();
        let model = bind_explicit(&native, source(&card, caloric.curve(), &context), binding(rotation()), cx).unwrap();
        let k = model.conductivity().model_for(0).unwrap();
        for (actual, expected) in [
            (k.tensor_at(400.).unwrap(), [[6., -3., 0.], [-3., 6., 0.], [0., 0., 1.5]]),
            (k.tensor_derivative_at(400.).unwrap(), [[0.02, -0.01, 0.], [-0.01, 0.02, 0.], [0., 0., 0.005]]),
        ] {
            for (a, b) in actual.iter().flatten().zip(expected.iter().flatten()) {
                near(*a, *b, 2e-14);
            }
        }
        let SpatialConductivitySource::OrientedOrthotropic { principal_selections, .. } = model.conductivity_source(0).unwrap() else { panic!("directional binding lost") };
        for (key, selection) in KEYS.iter().zip(principal_selections) {
            assert_eq!(*selection, ClaimSelection::Pinned(card.claims().claims_for(key)[0].0));
        }
        assert!(k.receipts().iter().all(|r| r.query_point.iter().any(|(axis, v)| axis == "P" && *v == 101325.)));
        let geometry = model.geometry(cx).unwrap();
        let loads = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        // Independent P1 manufacture: V=1/6, m_i cp=2.5, mean(T*)=405,
        // grad(T*)=(20,40,20), K(T*) grad(T*)=(0,183,30.5).
        // Each nodal source is 60*V/4=2.5 W. This also changes the mean by
        // one kelvin, so the nonlinear solve must update the source tensor.
        let target = [385., 405., 425., 405.];
        let outflow = [-213.5 / 6., 0., 183. / 6., 30.5 / 6.];
        let old: [f64; 4] = core::array::from_fn(|i| 10. * (target[i] + (outflow[i] - 2.5) / 2.5 - 300.));
        let bits = old.map(f64::to_bits);
        let result = model.advance(&old, 1., &boundary(&mesh), None, &ScalarField::Uniform(60.),
            thermal_config(), &loads, &(0..12).collect::<Vec<_>>(), solid_config(), cx).unwrap();
        for (t, expected) in result.thermal.temperature.iter().zip(target) { near(*t, expected, 2e-8); }
        near(result.thermal.stored_energy_change_j, 10., 1e-8);
        near(result.thermal.energy_residual_j, 0., 1e-8);
        let expected_stress = -737.5 * 1e-5 * (105. + 0.005 * 105. * 105.) / 0.5;
        for &stress in &result.equilibrium.solution.stress.elements()[0].stress_mandel_pa()[..3] {
            near(stress, expected_stress, 1e-7);
        }
        assert_eq!(old.map(f64::to_bits), bits);
        assert_eq!(native.reference_nodal_masses_kg(), &[0.25; 4]);
    });
}

#[test]
fn changing_only_material_orientation_changes_diffusion_and_actual_constrained_stress() {
    with_cx(|cx| {
        let (card, caloric) = card([[11., 1., 11.], [18., 8., 18.], [12., 2., 12.]]);
        let mut nodes = NODES.to_vec();
        nodes.push([2. / 3.; 3]);
        let mesh = ConductionMesh::new(TetComplex::from_tets(5, vec![[0, 1, 2, 3], [4, 1, 2, 3]]), nodes).unwrap();
        let native = HeterogeneousEnthalpyBackwardEuler::new(cx, &mesh,
            &[ReferenceEnthalpyMaterial { curve: caloric.curve(), reference_density_kg_m3: 6. }],
            &[0, 0], EnthalpyBudget { max_vertices: 5, max_elements: 2 },
        ).unwrap();
        let context = context();
        let n = [1. / 3.0_f64.sqrt(); 3];
        let u = [1. / 2.0_f64.sqrt(), -1. / 2.0_f64.sqrt(), 0.];
        let v = [1. / 6.0_f64.sqrt(), 1. / 6.0_f64.sqrt(), -2. / 6.0_f64.sqrt()];
        let r = core::array::from_fn(|i| [n[i], u[i], v[i]]);
        let turned = core::array::from_fn(|i| [u[i], -n[i], v[i]]);
        let mut stresses = Vec::new();
        for (rotation, normal_k_at_400) in [(r, 1.), (turned, 8.)] {
            let model = bind_explicit(&native, source(&card, caloric.curve(), &context), binding(rotation), cx).unwrap();
            let geometry = model.geometry(cx).unwrap();
            let loads = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
            let result = model.advance(&[800., 1000., 1000., 1000., 1200.], 1.,
                &boundary(&mesh), None, &ScalarField::Uniform(0.), thermal_config(),
                &loads, &(0..15).collect::<Vec<_>>(), solid_config(), cx).unwrap();
            // Mirror symmetry preserves the odd apex mode. Element means
            // are 400 +/- 5a, so k_n = k0 + 0.5a; native BE gives
            // 0.1 a^2 + (1 + k0/5)a - 1 = 0, with positive root below.
            let linear = 1. + normal_k_at_400 / 5.;
            let a = 2. / (linear + (linear * linear + 0.4_f64).sqrt());
            near(result.thermal.temperature[0], 400. - 20. * a, 2e-8);
            near(result.thermal.temperature[4], 400. + 20. * a, 2e-8);
            for &t in &result.thermal.temperature[1..4] { near(t, 400., 2e-8); }
            for (element, sign) in [-1., 1.].into_iter().enumerate() {
                let rise = 100. + sign * 5. * a;
                near(result.state.samples()[element].phase.temperature_k(), 300. + rise, 2e-8);
                let stress = -1000. * (1. - 0.0025 * rise) * 1e-5 * (rise + 0.005 * rise * rise) / 0.5;
                for &actual in &result.equilibrium.solution.stress.elements()[element].stress_mandel_pa()[..3] { near(actual, stress, 1e-7); }
            }
            near(result.thermal.stored_energy_change_j, 0., 1e-8);
            near(result.thermal.energy_residual_j, 0., 1e-8);
            stresses.push(result.equilibrium.solution.stress.elements()[0].stress_mandel_pa()[0]);
        }
        assert!((stresses[0] - stresses[1]).abs() > 1e-3);
    });
}

#[test]
fn original_scalar_constructor_and_equal_principal_limit_agree() {
    with_cx(|cx| {
        let (card, caloric) = fixture(2., 1000., 1e-5, 10., 500., 750.);
        let mesh = mesh();
        let native = native(&mesh, caloric.curve(), 6., cx);
        let context = point(300.);
        let original = bind(&native, source(&card, caloric.curve(), &context), cx).unwrap();
        let geometry = original.geometry(cx).unwrap();
        let loads = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let run = |model: &SpatialSolidEnthalpy<'_, '_, '_>| model.advance(
            &[1000., 1200., 800., 1000.], 1., &boundary(&mesh), None,
            &ScalarField::Uniform(0.), thermal_config(), &loads, &FREE, solid_config(), cx,
        ).unwrap().thermal.temperature;
        let expected = run(&original);
        let mut equal = binding(rotation());
        if let SpatialConductivitySource::OrientedOrthotropic { principal_properties, .. } = &mut equal {
            *principal_properties = core::array::from_fn(|_| PropertyKey::new("thermal_conductivity", CONDUCTIVITY_DIMS));
        }
        for (index, selection) in [SpatialConductivitySource::Isotropic, equal].into_iter().enumerate() {
            let selected = bind_explicit(&native, source(&card, caloric.curve(), &context), selection, cx).unwrap();
            let actual = run(&selected);
            if index == 0 { assert_eq!(actual, expected); }
            for (a, b) in actual.iter().zip(&expected) { near(*a, *b, 2e-9); }
        }
    });
}

#[test]
fn invalid_frames_missing_principal_laws_and_lost_conditions_refuse() {
    with_cx(|cx| {
        let (card, caloric) = card([[2.; 3], [6.; 3], [1.; 3]]);
        let mesh = mesh();
        let native = native(&mesh, caloric.curve(), 6., cx);
        let context = context();
        for defect in 0..5 {
            let mut bad = binding(rotation());
            if let SpatialConductivitySource::OrientedOrthotropic { principal_properties, principal_to_mesh, principal_frame, mesh_frame, .. } = &mut bad {
                match defect {
                    0 => *mesh_frame = ContentHash([52; 32]),
                    1 => principal_to_mesh[0][0] = f64::NAN,
                    2 => { for row in principal_to_mesh { row[0] = -row[0]; } }
                    3 => *principal_frame = ContentHash([0; 32]),
                    _ => principal_properties[1] = PropertyKey::new("missing_principal_k", CONDUCTIVITY_DIMS),
                }
            }
            assert!(bind_explicit(&native, source(&card, caloric.curve(), &context), bad, cx).is_err());
        }
        assert!(bind_explicit(&native, source(&card, caloric.curve(), &point(300.)), binding(rotation()), cx).is_err());
    });
}
