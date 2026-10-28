//! Source-paired variable-Cp transport reaches the original spatial enthalpy solve.
mod support;

use fs_conduction::transient::enthalpy::material::{
    MaterialEnthalpySource, MaterialEnthalpyTransport,
};
use fs_conduction::transient::enthalpy::{EnthalpyBudget, EnthalpyStepConfig};
use fs_conduction::{ConductionMesh, ConductionProblem, ScalarField, ThermalBoundaryBuilder};
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
