//! Source-paired variable-Cp transport reaches the original spatial enthalpy solve.
mod support;

use fs_conduction::transient::enthalpy::material::{
    MaterialEnthalpySource, MaterialEnthalpyTransport,
};
use fs_conduction::transient::enthalpy::{EnthalpyBudget, EnthalpyStepConfig};
use fs_conduction::{
    ConductionMesh, ConductionProblem, ConductivityModel, ConductivityTable, ScalarField,
    ThermalBoundaryBuilder,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::phase::{HeatCapacitySourceQuery, SolidLiquidPhase};
use fs_qty::{
    Density, Dims, QuantitySpec, Temperature,
    semantic::{MoistureBasis, QuantityKind, SemanticType, ValueForm},
};
use fs_rep_mesh::TetComplex;
use fs_solver::NewtonKrylovConfig;
use support::{with_cancelled_cx, with_cx};

const CP_DIMS: Dims = Dims([2, 0, -2, -1, 0, 0]);
const K_DIMS: Dims = Dims([1, 1, -3, -1, 0, 0]);

fn moisture() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::Moisture(MoistureBasis::DryMass),
        ValueForm::Static,
    ))
}
fn point(t: f64, water: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "temperature",
            QuantitySpec::dimensional(Temperature::DIMS),
            t,
        )
        .unwrap()
        .with_quantity("moisture", moisture(), water)
        .unwrap()
}
fn card(k_scale: f64, temperature_support: bool) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, lo, hi) in [
        ("rho", Density::DIMS, 10.0, 5.0),
        ("cp", CP_DIMS, 1000.0, 2000.0),
        ("k", K_DIMS, 2.0 * k_scale, 4.0 * k_scale),
        ("k_longitudinal", K_DIMS, 20.0 * k_scale, 40.0 * k_scale),
        ("k_radial", K_DIMS, 2.0 * k_scale, 4.0 * k_scale),
        ("k_tangential", K_DIMS, k_scale, 2.0 * k_scale),
    ] {
        let mut validity =
            ValidityDomain::unconstrained().with_quantity("moisture", moisture(), 0.10, 0.20);
        if temperature_support {
            validity = validity.with("temperature", 300.0, 400.0);
        }
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "temperature".into(),
                    abscissa_dims: Temperature::DIMS,
                    knots: vec![(300.0, lo), (400.0, hi)],
                    dims,
                },
                validity,
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: Vec::new(),
                provenance: Provenance {
                    source: "synthetic independent thermal reference; not measured material data"
                        .into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic thermal fixture".into(),
            phase: "solid".into(),
            process: "numerical reference".into(),
            revision: 0,
        },
        claims,
        Vec::new(),
    )
    .unwrap()
}
fn resolve(
    card: &MaterialCard,
    reference: &QueryPoint,
    conductivity: &PropertyKey,
) -> Result<
    MaterialEnthalpyTransport,
    fs_conduction::transient::enthalpy::material::MaterialEnthalpyError,
> {
    MaterialEnthalpyTransport::from_isotropic_card(
        card,
        MaterialEnthalpySource {
            caloric: HeatCapacitySourceQuery {
                heat_capacity: &PropertyKey::new("cp", CP_DIMS),
                density: &PropertyKey::new("rho", Density::DIMS),
                temperature_axis: "temperature",
                lower: &point(300.0, 0.12),
                upper: &point(400.0, 0.12),
                selection: SelectionPolicy::SingleClaimOnly,
            },
            conductivity,
            reference,
            phase: SolidLiquidPhase::Solid,
            lower_specific_enthalpy_j_kg: 0.0,
            maximum_interpolation_error_j_kg: 0.001,
        },
    )
}
fn material(card: &MaterialCard) -> MaterialEnthalpyTransport {
    resolve(card, &point(325.0, 0.12), &PropertyKey::new("k", K_DIMS)).unwrap()
}
fn resolve_orthotropic(
    card: &MaterialCard,
    principal_axes: [[f64; 3]; 3],
    conductivity: [&PropertyKey; 3],
) -> Result<
    MaterialEnthalpyTransport,
    fs_conduction::transient::enthalpy::material::MaterialEnthalpyError,
> {
    MaterialEnthalpyTransport::from_orthotropic_card(
        card,
        MaterialEnthalpySource {
            caloric: HeatCapacitySourceQuery {
                heat_capacity: &PropertyKey::new("cp", CP_DIMS),
                density: &PropertyKey::new("rho", Density::DIMS),
                temperature_axis: "temperature",
                lower: &point(300.0, 0.12),
                upper: &point(400.0, 0.12),
                selection: SelectionPolicy::SingleClaimOnly,
            },
            conductivity,
            reference: &point(325.0, 0.12),
            phase: SolidLiquidPhase::Solid,
            lower_specific_enthalpy_j_kg: 0.0,
            maximum_interpolation_error_j_kg: 0.001,
        },
        principal_axes,
    )
}
fn orthotropic(card: &MaterialCard, axes: [[f64; 3]; 3]) -> MaterialEnthalpyTransport {
    resolve_orthotropic(
        card,
        axes,
        [
            &PropertyKey::new("k_longitudinal", K_DIMS),
            &PropertyKey::new("k_radial", K_DIMS),
            &PropertyKey::new("k_tangential", K_DIMS),
        ],
    )
    .unwrap()
}

const AXES: [[f64; 3]; 3] = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
const TURNED_AXES: [[f64; 3]; 3] = [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0]];
fn mesh() -> ConductionMesh {
    ConductionMesh::new(
        TetComplex::from_tets(4, vec![[0, 1, 2, 3]]),
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ],
    )
    .unwrap()
}
fn budget() -> EnthalpyBudget {
    EnthalpyBudget {
        max_vertices: 4,
        max_elements: 1,
    }
}
fn config() -> EnthalpyStepConfig {
    EnthalpyStepConfig {
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-9,
            relative_tolerance: 1e-12,
            linear_restart: 4,
            max_linear_cycles: 8,
            forcing_minimum: 1e-12,
            forcing_maximum: 1e-4,
            ..NewtonKrylovConfig::default()
        },
        max_newton_iterations: 16,
        energy_tolerance_j: 1e-7,
    }
}

#[test]
fn integrated_capacity_and_reference_density_keep_the_original_source_context() {
    let card = card(1.0, true);
    let source = material(&card);
    assert_eq!(source.card_identity(), card.content_hash());
    assert_eq!(source.reference_query(), &point(325.0, 0.12));
    assert_eq!(source.reference_density_kg_m3(), 8.75);
    assert_eq!(source.reference_material().reference_density_kg_m3, 8.75);
    for receipt in source.caloric().receipts() {
        card.claims().verify_receipt(receipt).unwrap();
    }
    for receipt in source.conductivity().receipts() {
        card.claims().verify_receipt(receipt).unwrap();
    }
    card.claims()
        .verify_receipt(&source.reference_density().receipt)
        .unwrap();
    // Independent integral: h = 1000*(T-300) + 5*(T-300)^2.
    let state = source
        .caloric()
        .curve()
        .state_at_specific_enthalpy(62_500.0)
        .unwrap();
    assert!((state.temperature_k() - 350.0).abs() < 1e-6);
    assert!((state.bulk_density_kg_m3() - 7.5).abs() < 1e-6);
    assert_eq!(source.reference_density_kg_m3(), 8.75);
    assert_eq!(source.conductivity().tensor_at(350.0).unwrap()[0][0], 3.0);
    assert_eq!(source.identity(), material(&card).identity());
    assert_ne!(
        source.identity(),
        material(&self::card(2.0, true)).identity()
    );
    let other_reference =
        resolve(&card, &point(350.0, 0.12), &PropertyKey::new("k", K_DIMS)).unwrap();
    assert_ne!(source.identity(), other_reference.identity());
}

#[test]
fn source_pairing_refuses_changed_conditions_missing_support_and_wrong_quantities() {
    let card = card(1.0, true);
    for reference in [point(325.0, 0.13), point(401.0, 0.12), point(299.0, 0.12)] {
        assert!(resolve(&card, &reference, &PropertyKey::new("k", K_DIMS)).is_err());
    }
    assert!(
        resolve(
            &card,
            &point(325.0, 0.12),
            &PropertyKey::new("rho", Density::DIMS)
        )
        .is_err()
    );
    assert!(
        resolve(
            &card,
            &point(325.0, 0.12),
            &PropertyKey::new("absent", K_DIMS)
        )
        .is_err()
    );
    assert!(
        resolve(
            &self::card(1.0, false),
            &point(325.0, 0.12),
            &PropertyKey::new("k", K_DIMS)
        )
        .is_err()
    );
    let source = material(&card);
    assert!(source.conductivity().tensor_at(400.01).is_err());
    assert!(
        source
            .caloric()
            .curve()
            .state_at_specific_enthalpy(150_001.0)
            .is_err()
    );
}

#[test]
fn actual_spatial_heating_uses_frozen_reference_mass_and_integrated_cp() {
    let source = material(&card(1.0, true));
    let mesh = mesh();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let heat = ScalarField::Uniform(875.0);
    with_cx(|cx| {
        let stepper = source.prepare_uniform(cx, &mesh, budget()).unwrap();
        let old = vec![0.0; 4];
        let solved = stepper
            .advance(
                cx,
                ConductionProblem {
                    mesh: &mesh,
                    boundary: &boundary,
                    material: source.conductivity(),
                    element_materials: None,
                    source: &heat,
                },
                None,
                &old,
                10.0,
                config(),
            )
            .unwrap();
        // Independent balance: h_gain = (875 W/m3)*10s/(8.75 kg/m3).
        for &h in &solved.specific_enthalpy_j_kg {
            assert!((h - 1000.0).abs() < 1e-7);
        }
        let expected_t = 300.0 + 2000.0 / (1000.0 + (1_020_000.0_f64).sqrt());
        for &t in &solved.temperature {
            assert!((t - expected_t).abs() < 2e-6);
        }
        assert!((solved.stored_energy_change_j - 875.0 * 10.0 / 6.0).abs() < 1e-7);
        assert!(solved.energy_residual_j.abs() < 1e-7);
        assert_eq!(old, vec![0.0; 4]);
    });
    with_cancelled_cx(|cx| {
        assert!(source.prepare_uniform(cx, &mesh, budget()).is_err());
    });
}

#[test]
fn source_conductivity_controls_spatial_equilibration_without_creating_heat() {
    let mesh = mesh();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let heat = ScalarField::Uniform(0.0);
    let old = vec![10_000.0, 20_000.0, 10_000.0, 10_000.0];
    with_cx(|cx| {
        let mut spreads = Vec::new();
        for k_scale in [1.0, 2.0] {
            let source = material(&card(k_scale, true));
            let stepper = source.prepare_uniform(cx, &mesh, budget()).unwrap();
            let solved = stepper
                .advance(
                    cx,
                    ConductionProblem {
                        mesh: &mesh,
                        boundary: &boundary,
                        material: source.conductivity(),
                        element_materials: None,
                        source: &heat,
                    },
                    None,
                    &old,
                    1000.0,
                    config(),
                )
                .unwrap();
            // Independent total energy: each lump has rho0 * V / 4 kg.
            let mass_per_node = 8.75 / 24.0;
            let energy_change: f64 = solved
                .specific_enthalpy_j_kg
                .iter()
                .zip(&old)
                .map(|(now, before)| mass_per_node * (now - before))
                .sum();
            assert!(energy_change.abs() < 1e-7);
            assert!(solved.stored_energy_change_j.abs() < 1e-7);
            assert!(solved.energy_residual_j.abs() < 1e-7);
            assert_eq!(solved.source_w, 0.0);
            assert_eq!(solved.neumann_out_w, 0.0);
            assert_eq!(solved.robin_out_w, 0.0);
            assert_eq!(solved.dirichlet_in_w, 0.0);
            let minimum = solved
                .temperature
                .iter()
                .copied()
                .fold(f64::INFINITY, f64::min);
            let maximum = solved
                .temperature
                .iter()
                .copied()
                .fold(f64::NEG_INFINITY, f64::max);
            assert!(minimum > 300.0 && maximum < 400.0);
            assert!(solved.specific_enthalpy_j_kg[1] < old[1]);
            spreads.push(maximum - minimum);
        }
        // Identical storage and initial state; only the source k(T) doubled.
        assert!(spreads[1] < spreads[0]);
        assert!(spreads[0] - spreads[1] > 1e-3);
        assert_eq!(old, vec![10_000.0, 20_000.0, 10_000.0, 10_000.0]);
    });
}

#[test]
fn g0_directional_sources_rotate_both_conductivity_and_its_nonlinear_tangent() {
    let card = card(1.0, true);
    let s = core::f64::consts::FRAC_1_SQRT_2;
    let axes = [[s, s, 0.0], [-s, s, 0.0], [0.0, 0.0, 1.0]];
    let source = orthotropic(&card, axes);
    let tensor = source.conductivity().tensor_at(350.0).unwrap();
    let tangent = source.conductivity().tensor_derivative_at(350.0).unwrap();
    // Independent rotation of principal k=(30,3,1.5) and dk/dT=(.2,.02,.01).
    let expected = [[16.5, 13.5, 0.0], [13.5, 16.5, 0.0], [0.0, 0.0, 1.5]];
    let derivative = [[0.11, 0.09, 0.0], [0.09, 0.11, 0.0], [0.0, 0.0, 0.01]];
    for i in 0..3 {
        for j in 0..3 {
            assert!((tensor[i][j] - expected[i][j]).abs() < 1e-13);
            assert!((tangent[i][j] - derivative[i][j]).abs() < 1e-15);
        }
    }
    for receipt in source.conductivity().receipts() {
        card.claims().verify_receipt(receipt).unwrap();
    }
    assert_eq!(
        source.reference_density_kg_m3(),
        material(&card).reference_density_kg_m3()
    );
    assert_eq!(
        source.caloric().curve().identity(),
        material(&card).caloric().curve().identity()
    );
    assert_eq!(source.identity(), orthotropic(&card, axes).identity());
    assert_ne!(source.identity(), orthotropic(&card, AXES).identity());
    assert_ne!(
        orthotropic(&card, AXES).identity(),
        orthotropic(&card, TURNED_AXES).identity()
    );
    assert!(source.conductivity().tensor_at(400.01).is_err());
}

#[test]
fn g0_directional_pairing_refuses_missing_sources_and_nonfinite_frames() {
    let card = card(1.0, true);
    let k = PropertyKey::new("k", K_DIMS);
    for missing_or_wrong in [
        PropertyKey::new("absent", K_DIMS),
        PropertyKey::new("rho", Density::DIMS),
    ] {
        assert!(resolve_orthotropic(&card, AXES, [&k, &k, &missing_or_wrong]).is_err());
    }
    assert!(resolve_orthotropic(&self::card(1.0, false), AXES, [&k, &k, &k]).is_err());
    for invalid in [
        f64::NAN,
        f64::INFINITY,
        f64::NEG_INFINITY,
        0.0,
        2.0,
        f64::MAX,
    ] {
        let mut axes = AXES;
        axes[0][0] = invalid;
        assert!(resolve_orthotropic(&card, axes, [&k, &k, &k]).is_err());
        // The original operator constructor shares this admission; NaN must
        // never enter a tensor through a floating-point comparison loophole.
        assert!(
            ConductivityModel::orthotropic(
                axes,
                [
                    ConductivityTable::declared(1.0).unwrap(),
                    ConductivityTable::declared(2.0).unwrap(),
                    ConductivityTable::declared(3.0).unwrap(),
                ],
            )
            .is_err()
        );
    }
}

#[test]
fn g1_g3_data_only_orientation_changes_actual_heat_spreading_with_conserved_energy() {
    let card = card(1.0, true);
    let mesh = mesh();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap();
    let heat = ScalarField::Uniform(0.0);
    let old = vec![10_000.0, 20_000.0, 10_000.0, 10_000.0];
    with_cx(|cx| {
        let run = |axes, initial: &[f64]| {
            let source = orthotropic(&card, axes);
            let stepper = source.prepare_uniform(cx, &mesh, budget()).unwrap();
            let solved = stepper
                .advance(
                    cx,
                    ConductionProblem {
                        mesh: &mesh,
                        boundary: &boundary,
                        material: source.conductivity(),
                        element_materials: None,
                        source: &heat,
                    },
                    None,
                    initial,
                    1000.0,
                    config(),
                )
                .unwrap();
            // Independent mass balance; changing the material frame changes
            // transport but cannot reinterpret the reference mass or Cp.
            let energy: f64 = solved
                .specific_enthalpy_j_kg
                .iter()
                .zip(initial)
                .map(|(now, before)| (8.75 / 24.0) * (now - before))
                .sum();
            assert!(energy.abs() < 1e-7);
            assert!(solved.energy_residual_j.abs() < 1e-7);
            assert_eq!(solved.source_w, 0.0);
            assert_eq!(solved.neumann_out_w, 0.0);
            assert_eq!(solved.robin_out_w, 0.0);
            assert_eq!(solved.dirichlet_in_w, 0.0);
            solved
        };
        let longitudinal_x = run(AXES, &old);
        let longitudinal_y = run(TURNED_AXES, &old);
        assert!(longitudinal_x.temperature[1] < longitudinal_y.temperature[1] - 1e-3);
        // The tetrahedron and lumped storage are symmetric in x/y. Turning
        // both the material and the initial field must turn the solved field.
        let turned = run(TURNED_AXES, &[old[0], old[2], old[1], old[3]]);
        for (i, j) in [(0, 0), (1, 2), (2, 1), (3, 3)] {
            assert!((longitudinal_x.temperature[i] - turned.temperature[j]).abs() < 1e-9);
        }
        assert_eq!(old, vec![10_000.0, 20_000.0, 10_000.0, 10_000.0]);
    });
}
