//! Focused synthetic-source coupling regressions, not empirical membrane/gas qualification.
mod vibration;
mod motion;

use super::*;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, ConstitutiveModelCard, InitialStatePolicy, InterpolationPolicy, LawId, LawParameter,
    MaterialCard, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue, Provenance,
    QueryPoint, UncertaintyModel,
};
use fs_material::state_point::{MaterialPropertySelection, THERMAL_CONDUCTIVITY_DIMS};
use fs_material::thermochemical_gas::{SpeciesTransportMaterial, ThermochemicalGasModel};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{DynViscosity, MolarMass, Pressure, QuantitySpec, Temperature};
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
use fs_thermochem::{
    Composition, CompositionBasis, ElementalReferenceIdV1, NASA9_LAW_ID_V1, NASA9_LAW_VERSION_V1,
    NASA9_STATE_SCHEMA_VERSION_V1, Nasa9RegionV1, Nasa9StandardStateModelV1,
    ReferenceEquationOfStateV1, SpeciesId, StandardStateConventionV1, StandardStatePhaseV1,
    required_nasa9_parameter_dims,
};

fn quantity(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn species(name: &str, mass: f64, cp_r: f64, slope: f64) -> Nasa9StandardStateModelV1 {
    let mut coefficients = [0.0; 9];
    coefficients[2] = cp_r;
    coefficients[3] = slope;
    let mut parameters = std::collections::BTreeMap::new();
    for ((key, dims), value) in required_nasa9_parameter_dims()
        .into_iter()
        .take(9)
        .zip(coefficients)
    {
        parameters.insert(key.to_owned(), LawParameter { value, dims });
    }
    parameters.insert(
        "reference_pressure".into(),
        LawParameter {
            value: 100000.0,
            dims: Pressure::DIMS,
        },
    );
    let card = ConstitutiveModelCard {
        law: LawId(NASA9_LAW_ID_V1.into()),
        law_version: NASA9_LAW_VERSION_V1,
        parameters,
        state_schema_version: NASA9_STATE_SCHEMA_VERSION_V1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained().with("T", 200.0, 2000.0),
        sources: vec![],
        provenance: Provenance {
            source: "synthetic caloric gas".into(),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    };
    Nasa9StandardStateModelV1::new(
        SpeciesId::new(name).unwrap(),
        MolarMass::new(mass),
        StandardStateConventionV1::new(
            StandardStatePhaseV1::Gas,
            ReferenceEquationOfStateV1::IdealGas,
            Pressure::new(100000.0),
            ElementalReferenceIdV1::new("synthetic").unwrap(),
        )
        .unwrap(),
        vec![Nasa9RegionV1::from_card(card).unwrap()],
    )
    .unwrap()
}
fn transport(name: &str) -> SpeciesTransportMaterial {
    let mut claims = ClaimSet::new();
    for (key, dims, value) in [
        ("dynamic_viscosity", DynViscosity::DIMS, 2e-5),
        ("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS, 0.03),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(key, dims),
                value: PropertyValue::Curve {
                    abscissa: "temperature".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(200.0, value), (2000.0, 2.0 * value)],
                    dims,
                },
                validity: ValidityDomain::unconstrained()
                    .with_quantity(
                        "temperature",
                        quantity(QuantityKind::AbsoluteTemperature),
                        200.0,
                        2000.0,
                    )
                    .with_quantity("pressure", quantity(QuantityKind::Pressure), 1000.0, 1e7),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic transport".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    SpeciesTransportMaterial {
        species: SpeciesId::new(name).unwrap(),
        material: MaterialCard::assemble(
            MaterialStateId {
                chemistry: name.into(),
                phase: "gas".into(),
                process: "synthetic".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap(),
        selection: MaterialPropertySelection::SingleClaimOnly,
    }
}
fn source_gas(x: f64, slope: f64) -> ClosedGasParcel {
    let model = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![
                species("A", 0.028, 3.5, slope),
                species("B", 0.004, 2.5, 0.0),
            ],
            Composition::new(CompositionBasis::MoleFraction, vec![x, 1.0 - x]).unwrap(),
        )
        .unwrap(),
        vec![transport("A"), transport("B")],
        None,
    )
    .unwrap();
    let point = QueryPoint::new()
        .with_quantity(
            "temperature",
            quantity(QuantityKind::AbsoluteTemperature),
            400.0,
        )
        .unwrap()
        .with_quantity("pressure", quantity(QuantityKind::Pressure), 100000.0)
        .unwrap()
        .with("specimen-lot", 42.0)
        .unwrap();
    ClosedGasParcel::new(model, point, 0.001).unwrap()
}

use fs_material::sheet::{MEMBRANE_DIMS, SheetCapabilities, resolve_sheet};
use fs_material::thermochemical_gas::nonideal::{
    PengRobinsonSpeciesBinding, ThermochemicalPengRobinsonModel,
};
use fs_qty::Dims;
use fs_solid::membrane::{MembraneBudget, MembraneDomain, MembraneSolveSettings};
use fs_solid::mesh2::Mesh2;
use fs_solver::NewtonKrylovConfig;
use fs_thermochem::eos::peng_robinson::PengRobinsonPhaseV1;

const V0: f64 = 0.001;
const P0: f64 = 100000.0;
const GAMMA: f64 = 3.25 / 2.25; // Independent constant-Cp oracle for the 75:25 fixture.
const STRETCH: f64 = 1.02;
const SIDE: f64 = 0.1;
const STIFFNESS: f64 = 8000.0;

fn near(actual: f64, expected: f64, atol: f64, rtol: f64) {
    assert!((actual - expected).abs() <= atol + rtol * actual.abs().max(expected.abs()),
        "{actual:e} != {expected:e}");
}
fn gas_controls() -> GasEntropySolve {
    GasEntropySolve { minimum_temperature_k: 200.0, maximum_temperature_k: 1800.0,
        entropy_tolerance_j_kg_k: 1e-10, maximum_iterations: 128 }
}
fn controls() -> GasMembraneSolve {
    GasMembraneSolve { gas: gas_controls(), membrane: MembraneSolveSettings {
        newton: NewtonKrylovConfig { absolute_tolerance: 1e-8,
            relative_tolerance: 1e-10, linear_restart: 8, ..NewtonKrylovConfig::default() },
        max_iterations: 40, max_evaluations: 1000, max_krylov_restart: 8,
    } }
}
fn chart(mesh: Mesh2, free: &[usize]) -> (MembranePatch, Vec<[f64; 3]>, Vec<Option<f64>>) {
    let mut claims = ClaimSet::new();
    for (key, value, dims) in [
        ("areal_density", 0.2, Dims([-2, 1, 0, 0, 0, 0])),
        ("thickness", 0.001, Dims([1, 0, 0, 0, 0, 0])),
        ("warp_stiffness", STIFFNESS, MEMBRANE_DIMS),
        ("weft_stiffness", STIFFNESS, MEMBRANE_DIMS),
        ("poisson_warp_weft", 0.0, Dims::NONE),
        ("shear_stiffness", STIFFNESS / 2.0, MEMBRANE_DIMS),
    ] {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(key, dims), value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained().with("T", 290.0, 310.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity, observations: vec![],
            provenance: Provenance { source: "synthetic membrane; not measured fabric".into(),
                license: "CC0-1.0".into(), artifact: None },
        }).unwrap();
    }
    let card = MaterialCard::assemble(MaterialStateId { chemistry: "synthetic sheet".into(),
        phase: "solid".into(), process: "test".into(), revision: 0 }, claims, vec![]).unwrap();
    let sheet = resolve_sheet(&card, &QueryPoint::new().with("T", 300.0).unwrap(),
        SheetCapabilities { membrane: true, bending: false },
        MaterialPropertySelection::SingleClaimOnly).unwrap();
    let reference: Vec<_> = mesh.nodes.iter()
        .map(|x| [(STRETCH - 1.0) * x[0], (STRETCH - 1.0) * x[1], 0.0]).collect();
    let mut fixed: Vec<_> = reference.iter().flat_map(|u| u.iter().copied()).map(Some).collect();
    for &dof in free { fixed[dof] = None; }
    let count = mesh.elems.len();
    let patch = MembranePatch::try_new(mesh, vec![sheet; count], vec![0.0; count],
        MembraneDomain { minimum_area_ratio: 0.01, maximum_abs_strain: 0.2 },
        MembraneBudget { max_nodes: 32, max_triangles: 32 }, || false).unwrap();
    (patch, reference, fixed)
}
fn pyramid() -> (MembranePatch, Vec<[f64; 3]>, Vec<Option<f64>>) {
    chart(Mesh2 {
        nodes: vec![[0.0, 0.0], [SIDE, 0.0], [0.0, SIDE], [-SIDE, 0.0], [0.0, -SIDE]],
        elems: vec![vec![0, 1, 2], vec![0, 2, 3], vec![0, 3, 4], vec![0, 4, 1]],
        patches: vec![],
    }, &[2])
}
fn volume_gradient() -> f64 { 2.0 * (SIDE * STRETCH).powi(2) / 3.0 }
fn restoring(height: f64) -> f64 {
    2.0 * STIFFNESS * (STRETCH.powi(2) - 1.0) * height
        + 4.0 * STIFFNESS * height.powi(3) / SIDE.powi(2)
}

#[test]
fn g1_sealed_pyramid_equilibrium_and_work_match_independent_isentrope() {
    let (patch, reference, fixed) = pyramid();
    let source = source_gas(0.75, 0.0);
    let height = 0.006;
    let volume = V0 + volume_gradient() * height;
    let pressure = P0 * (V0 / volume).powf(GAMMA);
    let exterior = pressure - restoring(height) / volume_gradient();
    let model = GasMembrane::new(&patch, source.clone(), &reference, &fixed,
        exterior, [0.0007, 0.0015], || false).unwrap();
    let result = model.solve(&[[0.0; 3]; 5], &reference, controls(), || false).unwrap();
    let state = &result.equilibrium;
    near(state.displacement_m[0][2], height, 2e-10, 2e-10);
    near(state.gas.volume_m3(), volume, 2e-12, 0.0);
    near(state.pressure.pressure_pa, pressure, 2e-5, 2e-10);
    near(state.residual_n[0][2], 0.0, 1e-8, 0.0);
    let work = P0 * V0 / (GAMMA - 1.0) * (1.0 - (V0 / volume).powf(GAMMA - 1.0));
    near(state.gas_volume.work_by_gas_j, work, 2e-8, 2e-9);
    let elastic = STIFFNESS * (STRETCH.powi(2) - 1.0) * height.powi(2)
        + STIFFNESS * height.powi(4) / SIDE.powi(2);
    near(state.potential_energy_change_j, elastic - work + exterior * (volume - V0),
        2e-8, 2e-8);
    assert!((work - pressure * (volume - V0)).abs() > 0.02);
    near(state.residual_n.iter().map(|f| f[2]).sum(),
        -(pressure - exterior) * 3.0 * volume_gradient(), 1e-6, 2e-8);
    assert_eq!(state.gas.mass_kg(), source.mass_kg());
    assert_eq!(model.reference_gas(), &MembraneGas::Ideal(source));
    assert!(result.gas_trials < result.evaluations);
}

#[test]
fn g1_inadmissible_newton_trial_backtracks_without_poisoning_the_gas() {
    let (patch, reference, fixed) = pyramid();
    let height = -0.01;
    let volume = V0 + volume_gradient() * height;
    let exterior = P0 * (V0 / volume).powf(GAMMA) - restoring(height) / volume_gradient();
    let lower = V0 + volume_gradient() * (1.05 * height);
    let initial_stiffness = 2.0 * STIFFNESS * (STRETCH.powi(2) - 1.0)
        + GAMMA * P0 / V0 * volume_gradient().powi(2);
    let first_step = (P0 - exterior) * volume_gradient() / initial_stiffness;
    assert!(V0 + volume_gradient() * first_step < lower);
    let model = GasMembrane::new(&patch, source_gas(0.75, 0.0), &reference,
        &fixed, exterior, [lower, 0.0015], || false).unwrap();
    let result = model.solve(&[[0.0; 3]; 5], &reference, controls(), || false).unwrap();
    near(result.equilibrium.displacement_m[0][2], height, 2e-10, 2e-10);
    assert!(result.newton.converged);
}

#[test]
fn g1_zero_pressure_jump_retains_the_finite_gas_stiffness() {
    let (patch, reference, fixed) = pyramid();
    let model = GasMembrane::new(&patch, source_gas(0.75, 0.0), &reference,
        &fixed, P0, [0.0007, 0.0015], || false).unwrap();
    let state = model.evaluate(&reference, &[[0.0; 3]; 5], gas_controls(), || false).unwrap();
    near(state.residual_n[0][2], 0.0, 1e-12, 0.0);
    assert_eq!(state.gas_volume.work_by_gas_j, 0.0);
    let tangent = model.tangent_apply(&reference, &[1.0], gas_controls(), || false).unwrap();
    let structural = 2.0 * STIFFNESS * (STRETCH.powi(2) - 1.0);
    let gas_stiffness = GAMMA * P0 / V0 * volume_gradient().powi(2);
    near(tangent[0], structural + gas_stiffness, 1e-8, 2e-12);
    assert!(gas_stiffness > structural);
}

#[test]
fn g1_six_coordinate_tangent_and_potential_include_geometric_pressure() {
    let (patch, reference, fixed) = chart(Mesh2::triangles(0.3, 0.2, 3, 2),
        &[15, 16, 17, 18, 19, 20]);
    let exterior = 80000.0;
    let model = GasMembrane::new(&patch, source_gas(0.75, 0.0), &reference,
        &fixed, exterior, [0.0007, 0.0015], || false).unwrap();
    let mut u = reference.clone();
    u[5][0] += 0.0004; u[5][1] -= 0.0003; u[5][2] = 0.006;
    u[6][0] -= 0.0006; u[6][1] += 0.0002; u[6][2] = 0.009;
    let mut loads = vec![[0.0; 3]; u.len()];
    loads[5] = [0.3, -0.2, 0.1]; loads[6] = [-0.1, 0.2, -0.3];
    let value = model.evaluate(&u, &loads, gas_controls(), || false).unwrap();
    let h = 1e-6;
    let mut geometric_action = 0.0_f64;
    for (column, &dof) in model.free_dofs().iter().enumerate() {
        let mut lo = u.clone(); let mut hi = u.clone();
        lo[dof / 3][dof % 3] -= h; hi[dof / 3][dof % 3] += h;
        let a = model.evaluate(&lo, &loads, gas_controls(), || false).unwrap();
        let b = model.evaluate(&hi, &loads, gas_controls(), || false).unwrap();
        near((b.potential_energy_change_j - a.potential_energy_change_j) / (2.0 * h),
            value.residual_n[dof / 3][dof % 3], 1e-4, 1e-5);
        let mut direction = vec![0.0; 6]; direction[column] = 1.0;
        let tangent = model.tangent_apply(&u, &direction, gas_controls(), || false).unwrap();
        for (&row, &actual) in model.free_dofs().iter().zip(&tangent) {
            near((b.residual_n[row / 3][row % 3] - a.residual_n[row / 3][row % 3]) / (2.0 * h),
                actual, 4e-4, 4e-6);
        }
        let mut full = vec![[0.0; 3]; u.len()]; full[dof / 3][dof % 3] = 1.0;
        let geometric = patch.pressure_tangent_apply(&u, &full,
            value.pressure.pressure_pa - exterior, || false).unwrap();
        for &row in model.free_dofs() {
            geometric_action = geometric_action.max(geometric[row / 3][row % 3].abs());
        }
    }
    assert!(geometric_action > 1.0, "the geometric-pressure term must be exercised");
}

#[test]
fn g1_temperature_dependent_calorics_follow_original_entropy_and_work() {
    let (patch, reference, fixed) = pyramid();
    let model = GasMembrane::new(&patch, source_gas(0.75, 0.002), &reference,
        &fixed, P0, [0.0007, 0.0015], || false).unwrap();
    let mut u = reference.clone(); u[0][2] = 0.02;
    let value = model.evaluate(&u, &[[0.0; 3]; 5], gas_controls(), || false).unwrap();
    let volume = V0 + volume_gradient() * u[0][2];
    // Independent analytic entropy for Cv/R = 2.25 + 0.0015*T.
    let (mut lo, mut hi) = (200.0_f64, 400.0_f64);
    for _ in 0..80 {
        let t = (lo + hi) / 2.0;
        let entropy = 2.25 * (t / 400.0).ln() + 0.0015 * (t - 400.0) + (volume / V0).ln();
        if entropy > 0.0 { hi = t; } else { lo = t; }
    }
    let t = (lo + hi) / 2.0;
    let expected_pressure = P0 * (t / 400.0) * (V0 / volume);
    let energy = P0 * V0 / 400.0
        * (2.25 * (t - 400.0) + 0.5 * 0.0015 * (t * t - 400.0_f64.powi(2)));
    near(value.pressure.pressure_pa, expected_pressure, 1e-5, 2e-10);
    near(value.gas_volume.stored_energy_change_j, energy, 1e-8, 2e-9);
    let frozen_gamma = 1.0 + 1.0 / (2.25 + 0.0015 * 400.0);
    assert!((expected_pressure - P0 * (V0 / volume).powf(frozen_gamma)).abs() > 1.0);
    let MembraneGas::Ideal(accepted) = &value.gas else { panic!("original gas owner lost") };
    assert_eq!(accepted.point().axes().get("specimen-lot"), Some(&42.0));
}

#[test]
fn g1_peng_robinson_preserves_real_pressure_bulk_modulus_and_work() {
    let ideal = source_gas(0.75, 0.001);
    let source = ThermochemicalPengRobinsonModel::new(ideal.model().clone(), vec![
        PengRobinsonSpeciesBinding { species: SpeciesId::new("A").unwrap(),
            critical_temperature: Temperature::new(190.564),
            critical_pressure: Pressure::new(4.5992e6), acentric_factor: 0.01142 },
        PengRobinsonSpeciesBinding { species: SpeciesId::new("B").unwrap(),
            critical_temperature: Temperature::new(150.687),
            critical_pressure: Pressure::new(4.863e6), acentric_factor: -0.00219 },
    ], vec![0.0, 0.03, 0.03, 0.0]).unwrap();
    let point = QueryPoint::new()
        .with_quantity("temperature", quantity(QuantityKind::AbsoluteTemperature), 300.0).unwrap()
        .with_quantity("pressure", quantity(QuantityKind::Pressure), 5e6).unwrap()
        .with("specimen-lot", 42.0).unwrap();
    let initial = ClosedRealGasParcel::new(source, point, V0, PengRobinsonPhaseV1::Vapor).unwrap();
    let entropy = GasEntropySolve { minimum_temperature_k: 280.0,
        maximum_temperature_k: 450.0, ..gas_controls() };
    let (patch, reference, fixed) = pyramid();
    let model = GasMembrane::new(&patch, initial.clone(), &reference, &fixed,
        5e6, [0.0007, 0.0015], || false).unwrap();
    let mut u = reference.clone(); u[0][2] = 0.006;
    let value = model.evaluate(&u, &[[0.0; 3]; 5], entropy, || false).unwrap();
    let mut direct = initial.clone();
    let work = direct.set_volume_isentropic(value.gas.volume_m3(), entropy, || false).unwrap();
    let MembraneGas::PengRobinson(accepted) = &value.gas else { panic!("PR owner lost") };
    assert_eq!(accepted, &direct);
    assert_eq!(value.gas_volume, work);
    assert_eq!(value.pressure, direct.pressure_response().unwrap());
    assert_eq!(accepted.model(), initial.model());
    assert_eq!(accepted.mass_kg(), initial.mass_kg());
    let mut ideal = ClosedGasParcel::new(initial.model().source().clone(),
        initial.point().clone(), V0).unwrap();
    ideal.set_volume_isentropic(value.gas.volume_m3(), entropy, || false).unwrap();
    assert!((value.pressure.pressure_pa - ideal.state().state().pressure).abs() > 1.0);
    let mut direction = vec![[0.0; 3]; 5]; direction[0][2] = 1.0;
    let structural = patch.tangent_apply(&u, &direction, || false).unwrap()[0][2];
    let expected = structural + value.pressure.bulk_modulus_pa / value.gas.volume_m3()
        * volume_gradient().powi(2);
    near(model.tangent_apply(&u, &[1.0], entropy, || false).unwrap()[0],
        expected, 1e-6, 2e-11);
    assert_eq!(model.reference_gas(), &MembraneGas::PengRobinson(initial));
}

#[test]
fn g4_refused_geometry_budget_and_cancellation_preserve_original_inventory() {
    let (patch, reference, fixed) = pyramid();
    let initial = source_gas(0.75, 0.001);
    let model = GasMembrane::new(&patch, initial.clone(), &reference,
        &fixed, 90000.0, [0.00099, 0.00101], || false).unwrap();
    let mut moved = reference.clone(); moved[0][2] = 0.006;
    assert!(matches!(model.evaluate(&moved, &[[0.0; 3]; 5], gas_controls(), || false),
        Err(GasMembraneError::Invalid(_))));
    moved = reference.clone(); moved[1][2] = 0.0001;
    assert!(matches!(model.evaluate(&moved, &[[0.0; 3]; 5], gas_controls(), || false),
        Err(GasMembraneError::Invalid(_))));
    let mut budget = controls(); budget.membrane.max_evaluations = 1;
    assert!(matches!(model.solve(&[[0.0; 3]; 5], &reference, budget, || false),
        Err(GasMembraneError::EvaluationBudget)));
    assert!(matches!(model.solve(&[[0.0; 3]; 5], &reference, controls(), || true),
        Err(GasMembraneError::Cancelled)));
    let mut open = fixed.clone(); open[5] = None;
    assert!(matches!(GasMembrane::new(&patch, initial.clone(), &reference,
        &open, P0, [0.0007, 0.0015], || false), Err(GasMembraneError::Membrane(_))));
    assert_eq!(model.reference_gas(), &MembraneGas::Ideal(initial));
}
