//! Original nodal calorics and native tetrahedral reciprocal thermal work.

use super::*;
use fs_couple::solid_enthalpy::spatial::thermoelastic::{
    ReciprocalSpatialSolid, SpatialCaloricReference, SpatialThermoelasticError,
};
use fs_couple::solid_enthalpy::thermoelastic::ThermoelasticSolve;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

const T0: f64 = 300.0;
const LOW: f64 = 280.0;
const HIGH: f64 = 330.0;
const E0: f64 = 60.0;
const ET: f64 = -0.6;
const NU: f64 = 0.25;
const RHO: f64 = 6.0;

fn absolute() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(QuantityKind::AbsoluteTemperature, ValueForm::Static))
}
fn typed_point(t: f64) -> QueryPoint {
    QueryPoint::new().with_quantity("T", absolute(), t).unwrap()
}
fn physical_controls() -> ThermoelasticSolve {
    ThermoelasticSolve {
        maximum_iterations: 96,
        heat_relaxation: 0.8,
        heat_tolerance_j: 2e-9,
        energy_tolerance_j: 5e-8,
        stress_tolerance_pa: 2e-9,
    }
}
fn chart(alpha: f64, cp: f64) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let mut claims = ClaimSet::new();
    for (name, dims, value, slope) in [
        ("density", Density::DIMS, RHO, 0.0),
        ("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS, cp, 0.0),
        ("thermal_conductivity", CONDUCTIVITY_DIMS, 2.0, 0.0),
        ("young_modulus", Pressure::DIMS, E0, ET),
        ("poisson_ratio", Dims::NONE, NU, 0.0),
        ("linear_thermal_expansion_coefficient", INVERSE_TEMPERATURE_DIMS, alpha, alpha / 100.0),
    ] {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Curve {
                abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                knots: [LOW, HIGH].map(|t| (t, value + slope * (t - T0))).to_vec(), dims,
            },
            validity: ValidityDomain::unconstrained().with_quantity("T", absolute(), LOW, HIGH),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "analytical reciprocal P1 thermoelastic fixture".into(),
                license: "test-only".into(), artifact: None,
            },
        }).unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured nodal solid".into(), phase: "solid".into(),
            process: "analytical source".into(), revision: 0,
        },
        claims, vec![],
    ).unwrap();
    let caloric = ResolvedHeatCapacityCurve::try_new_exact(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
            density: &PropertyKey::new("density", Density::DIMS),
            temperature_axis: "T", lower: &typed_point(LOW), upper: &typed_point(HIGH),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid, 0.0,
    ).unwrap();
    (card, caloric)
}
fn original_source<'a>(
    card: &'a MaterialCard, curve: &'a EquilibriumEnthalpyPhaseCurve,
    context: &'a QueryPoint, cp: f64,
) -> SpatialMaterialSource<'a> {
    let mut result = super::source(card, curve, context);
    result.reference_specific_enthalpy_j_kg = cp * (T0 - LOW);
    result.temperature_interval_k = [LOW, HIGH];
    result
}
fn free_strain(alpha: f64, theta: f64) -> f64 {
    alpha * (theta + theta * theta / 200.0)
}

// Independent isotropic Helmholtz/entropy equations in physical strain, using
// the original engineering constants rather than the material exchange under test.
fn energy(theta: f64, alpha: f64, strain: [f64; 3], volume: f64) -> f64 {
    let e = E0 + ET * theta;
    let bulk = e / (3.0 * (1.0 - 2.0 * NU));
    let shear = e / (2.0 * (1.0 + NU));
    let bulk_t = ET / (3.0 * (1.0 - 2.0 * NU));
    let shear_t = ET / (2.0 * (1.0 + NU));
    let trace = strain.iter().sum::<f64>();
    let dev2 = strain.iter().map(|x| (x - trace / 3.0).powi(2)).sum::<f64>();
    let dilation = trace - 3.0 * free_strain(alpha, theta);
    let psi = volume * (0.5 * bulk * dilation * dilation + shear * dev2);
    let entropy = volume * (
        3.0 * bulk * alpha * (1.0 + theta / 100.0) * dilation
        - 0.5 * bulk_t * dilation * dilation - shear_t * dev2
    );
    psi + (T0 + theta) * entropy
}
fn root(mut lo: f64, mut hi: f64, residual: impl Fn(f64) -> f64) -> f64 {
    assert!(residual(lo) <= 0.0 && residual(hi) >= 0.0);
    for _ in 0..64 {
        let mid = 0.5 * (lo + hi);
        if residual(mid) > 0.0 { hi = mid; } else { lo = mid; }
    }
    0.5 * (lo + hi)
}

#[test]
fn finite_nodal_heat_distinguishes_free_expansion_and_restrained_reciprocal_heat() {
    with_cx(|cx| {
        for alpha in [1e-3, -1e-3] {
            let cp = 10.0;
            let (card, caloric) = chart(alpha, cp);
            let mesh = mesh();
            let native = native(&mesh, caloric.curve(), RHO, cx);
            let context = typed_point(T0);
            let binding = bind(&native, original_source(&card, caloric.curve(), &context, cp), cx).unwrap();
            let geometry = binding.geometry(cx).unwrap();
            let loads = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
            let boundary = boundary(&mesh);
            let power = ScalarField::Uniform(600.0);
            let volume = 1.0 / 6.0;
            let heat = 100.0;
            let old = [cp * (T0 - LOW); 4];
            for restrained in [false, true] {
                let fixed = if restrained { (0..12).collect::<Vec<_>>() } else { FREE.to_vec() };
                let model = ReciprocalSpatialSolid::try_new(
                    &binding, &old, 0.0, &loads, &fixed, solid_config(), physical_controls(),
                    SpatialCaloricReference::StressFreeSolid, cx,
                ).unwrap();
                assert!(model.minimum_nodal_heat_capacities_j_k().iter().all(|c| *c > 0.0));
                let result = model.advance(
                    1.0, &boundary, None, &power, &loads, thermal_config(), physical_controls(), cx,
                ).unwrap();
                let capacity = RHO * volume * cp;
                let theta = if restrained {
                    root(0.0, 20.0, |t| capacity * t + energy(t, alpha, [0.0; 3], volume) - heat)
                } else { heat / capacity };
                for &t in &result.thermal.endpoint.temperature { near(t, T0 + theta, 3e-8); }
                let u = result.solid.equilibrium().solution.displacement.displacement_m();
                near(u[1][0], if restrained { 0.0 } else { free_strain(alpha, theta) }, 3e-10);
                let stress = result.solid.equilibrium().solution.stress.elements()[0].stress_mandel_pa();
                let expected = if restrained {
                    -(E0 + ET * theta) / (1.0 - 2.0 * NU) * free_strain(alpha, theta)
                } else { 0.0 };
                for &sigma in &stress[..3] { near(sigma, expected, 2e-8); }
                near(result.report.external_thermal_heat_j, heat, 1e-9);
                near(result.report.external_work_j, 0.0, 1e-12);
                near(result.report.reversible_heat_j, capacity * theta - heat, 5e-8);
                near(result.report.total_energy_residual_j, 0.0, 5e-8);
                assert!(result.thermal.endpoint.maximum_free_residual_j <= thermal_config().energy_tolerance_j);
                for q in &result.report.nodal_reversible_heat_j {
                    near(q.to_f64(), result.report.reversible_heat_j / 4.0, 1e-12);
                }
                assert_eq!(model.specific_enthalpy_j_kg(), &old);
                assert_eq!(model.time_s(), 0.0);
                assert_eq!(result.solid.time_s(), 1.0);
            }
        }
    });
}

#[test]
fn spatial_closed_heat_preserves_native_graded_inventory_and_local_element_exchange() {
    with_cx(|cx| {
        let cp = 10.0;
        let (card, caloric) = chart(1e-3, cp);
        let mut positions = NODES.to_vec();
        positions.push([1.0, 1.0, 1.0]);
        let mesh = ConductionMesh::new(
            TetComplex::from_tets(5, vec![[0,1,2,3], [1,2,3,4]]), positions,
        ).unwrap();
        let native = HeterogeneousEnthalpyBackwardEuler::new(
            cx, &mesh,
            &[ReferenceEnthalpyMaterial { curve: caloric.curve(), reference_density_kg_m3: RHO }],
            &[0,0], EnthalpyBudget { max_vertices: 5, max_elements: 2 },
        ).unwrap().with_nodal_reference_density(cx, &[6.0, 9.0, 12.0, 15.0, 18.0]).unwrap();
        let masses = native.reference_nodal_masses_kg().to_vec();
        let context = typed_point(T0);
        let binding = bind(&native, original_source(&card, caloric.curve(), &context, cp), cx).unwrap();
        let geometry = binding.geometry(cx).unwrap();
        let loads = geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let boundary = boundary(&mesh);
        let old = [315.0, 305.0, 295.0, 290.0, 300.0].map(|t| cp * (t - LOW));
        let fixed: Vec<_> = (0..15).collect();
        let model = ReciprocalSpatialSolid::try_new(
            &binding, &old, 0.0, &loads, &fixed, solid_config(), physical_controls(),
            SpatialCaloricReference::StressFreeSolid, cx,
        ).unwrap();
        let result = model.advance(
            0.5, &boundary, None, &ScalarField::Uniform(0.0), &loads,
            thermal_config(), physical_controls(), cx,
        ).unwrap();
        let enthalpy = result.solid.specific_enthalpy_j_kg();
        assert_ne!(enthalpy, &old);
        let caloric_change: f64 = masses.iter().zip(enthalpy).zip(old)
            .map(|((&m, &h), old)| m * (h - old)).sum();
        near(caloric_change + result.report.internal_energy_change_j, 0.0, 5e-8);
        near(result.report.external_thermal_heat_j, 0.0, 1e-12);
        near(result.report.external_work_j, 0.0, 1e-12);
        near(result.report.total_energy_residual_j, 0.0, 5e-8);
        let q0 = result.report.elements[0].reversible_heat_j / 4.0;
        let q1 = result.report.elements[1].reversible_heat_j / 4.0;
        assert!((q0 - q1).abs() > 1e-8);
        let expected = [q0, q0+q1, q0+q1, q0+q1, q1];
        for (actual, expected) in result.report.nodal_reversible_heat_j.iter().zip(expected) {
            near(actual.to_f64(), expected, 1e-12);
        }
        assert!(result.thermal.endpoint.residual_j.iter()
            .all(|r| r.abs() <= thermal_config().energy_tolerance_j));
        assert_eq!(masses.as_slice(), native.reference_nodal_masses_kg());
        assert_eq!(model.specific_enthalpy_j_kg(), &old);
    });
}

#[test]
fn loaded_softening_matches_independent_work_and_failed_proposals_are_atomic() {
    with_cx(|cx| {
        let (alpha, cp, volume) = (1e-3, 10.0, 1.0/6.0);
        let (card, caloric) = chart(alpha, cp);
        let mesh = mesh();
        let native = native(&mesh, caloric.curve(), RHO, cx);
        let context = typed_point(T0);
        let binding = bind(&native, original_source(&card, caloric.curve(), &context, cp), cx).unwrap();
        let geometry = binding.geometry(cx).unwrap();
        let before = geometry.assemble(&[], VolumeLoad3::None, &[(3, 0.1)], cx).unwrap();
        let after = geometry.assemble(&[], VolumeLoad3::None, &[(3, 0.15)], cx).unwrap();
        let fixed: Vec<_> = (0..12).filter(|&d| d != 3).collect();
        let old = [cp * (T0 - LOW); 4];
        let model = ReciprocalSpatialSolid::try_new(
            &binding, &old, 0.0, &before, &fixed, solid_config(), physical_controls(),
            SpatialCaloricReference::StressFreeSolid, cx,
        ).unwrap();
        let u0 = model.equilibrium().solution.displacement.displacement_m().to_vec();
        let boundary = boundary(&mesh);
        let power = ScalarField::Uniform(120.0);
        let mut bounded = physical_controls();
        bounded.maximum_iterations = 1;
        assert!(model.advance(
            1.0, &boundary, None, &power, &after, thermal_config(), bounded, cx,
        ).is_err());
        assert!(model.advance(
            100.0, &boundary, None, &power, &after, thermal_config(), physical_controls(), cx,
        ).is_err());
        let result = model.advance(
            1.0, &boundary, None, &power, &after, thermal_config(), physical_controls(), cx,
        ).unwrap();
        let repeat = model.advance(
            1.0, &boundary, None, &power, &after, thermal_config(), physical_controls(), cx,
        ).unwrap();
        assert_eq!(result.thermal.specific_enthalpy_j_kg, repeat.thermal.specific_enthalpy_j_kg);
        let strain = |theta: f64, force: f64| {
            let e = E0 + ET * theta;
            let constrained = e * (1.0 - NU) / ((1.0 + NU) * (1.0 - 2.0 * NU));
            (force / volume + e / (1.0 - 2.0 * NU) * free_strain(alpha, theta)) / constrained
        };
        let start = strain(0.0, 0.1);
        let start_energy = energy(0.0, alpha, [start, 0.0, 0.0], volume);
        let theta = root(0.0, 10.0, |theta| {
            let end = strain(theta, 0.15);
            RHO * volume * cp * theta + energy(theta, alpha, [end, 0.0, 0.0], volume)
                - start_energy - 0.125 * (end - start) - 20.0
        });
        let end = strain(theta, 0.15);
        for &t in &result.thermal.endpoint.temperature { near(t, T0 + theta, 3e-8); }
        near(result.solid.equilibrium().solution.displacement.displacement_m()[1][0], end, 3e-10);
        near(result.report.external_work_j, 0.125 * (end - start), 2e-10);
        near(result.report.internal_energy_change_j,
            energy(theta, alpha, [end, 0.0, 0.0], volume) - start_energy, 2e-8);
        near(result.report.total_energy_residual_j, 0.0, 5e-8);
        assert_eq!(model.specific_enthalpy_j_kg(), &old);
        assert_eq!(model.equilibrium().solution.displacement.displacement_m(), u0.as_slice());
        let gate = CancelGate::new();
        gate.request();
        ArenaPool::new(ArenaConfig::default()).scope(|arena| {
            let cancelled = Cx::new(
                &gate, arena,
                StreamKey { seed: 3, kernel_id: 17, tile: 0, iteration: 0 },
                Budget::INFINITE, ExecMode::Deterministic,
            );
            assert!(matches!(model.advance(
                1.0, &boundary, None, &power, &after, thermal_config(), physical_controls(), &cancelled,
            ), Err(SpatialThermoelasticError::Cancelled)));
        });
        let (weak_card, weak_caloric) = chart(alpha, 1e-4);
        let weak_native = super::native(&mesh, weak_caloric.curve(), RHO, cx);
        let weak_binding = bind(
            &weak_native, original_source(&weak_card, weak_caloric.curve(), &context, 1e-4), cx,
        ).unwrap();
        let weak_geometry = weak_binding.geometry(cx).unwrap();
        let weak_load = weak_geometry.assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        assert!(ReciprocalSpatialSolid::try_new(
            &weak_binding, &[1e-4 * (T0 - LOW); 4], 0.0, &weak_load, &fixed,
            solid_config(), physical_controls(), SpatialCaloricReference::StressFreeSolid, cx,
        ).is_err());
    });
}

#[path = "thermoelastic/dynamic.rs"]
mod dynamic;

#[path = "thermoelastic/contact.rs"]
mod contact;
