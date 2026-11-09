//! Original thermal Hex8 histories drive pointwise H(T, alpha), without new heat.
//! Authored numerical fixtures, not measured material or lifetime qualification.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, HardnessLoadStep, HardnessTestContext, InterpolationPolicy, MaterialCard,
    MaterialStateId, ObservationDataset, PropertyClaim, PropertyKey, PropertyValue, Provenance,
    QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::SmallStrainLaw;
use fs_material::hardness::INDENTATION_HARDNESS_QUANTITY;
use fs_material::hardness::plastic_history::ACCUMULATED_EQUIVALENT_PLASTIC_STRAIN_AXIS as ALPHA;
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::plastic::chaboche::source::BackstressProperties;
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::source::{
    HardeningPropertyBasis, ISOTROPIC_HARDENING_MODULUS_PROPERTY,
    KINEMATIC_HARDENING_MODULUS_PROPERTY,
};
use fs_material::plastic::thermal::{SourceQuery, ThermoplasticSource};
use fs_material::plastic::thermal::family::ThermoplasticSourceFamily;
use fs_material::plastic::voce::source::{
    VOCE_SATURATION_INCREMENT_PROPERTY, VOCE_SATURATION_RATE_PROPERTY,
};
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, POISSON_RATIO_PROPERTY, YIELD_STRESS_PROPERTY,
    YOUNG_MODULUS_PROPERTY,
};
use fs_qty::{Density, Dims, Pressure, QtyAny, QuantitySpec, Temperature};
use fs_solid::plastic3::hex::HexAssemblyBudget;
use fs_solid::plastic3::hex::thermal::{
    ThermalHexHistory, ThermalHexProblem, ThermalHexTolerances,
};
use fs_solid::plastic3::hex::thermal::hardness::ThermalHexHardnessError;
use fs_solid::{Formulation, IncrementSettings};

const SELECT: MaterialPropertySelection = MaterialPropertySelection::SingleClaimOnly;
const CP: f64 = 500.0;
const CP_DIMS: Dims = Dims([2, 0, -2, -1, 0, 0]);
const YIELD_TOL: f64 = 1e-7;
const NODES: [[f64; 3]; 8] = [
    [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 1.0, 0.0], [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0], [1.0, 0.0, 1.0], [1.0, 1.0, 1.0], [0.0, 1.0, 1.0],
];
const CELLS: [[usize; 8]; 1] = [[0, 1, 2, 3, 4, 5, 6, 7]];

fn near(a: f64, b: f64) {
    assert!(a.is_finite() && (a - b).abs() <= 1e-9, "{a:e} != {b:e}");
}
fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(&gate, arena, StreamKey {
            seed: 17, kernel_id: 91, tile: 0, iteration: 0,
        }, Budget::INFINITE, ExecMode::Deterministic), &gate);
    });
}
fn expected_hardness(temperature: f64, alpha: f64) -> f64 {
    500.0 + 4000.0 * alpha - 3.0 * (temperature - 300.0)
}
struct Fixture {
    card: MaterialCard,
    caloric: ResolvedHeatCapacityCurve,
    keys: [PropertyKey; 2],
}
impl Fixture {
    fn new(maximum_hardness_temperature: f64) -> Self {
        let provenance = Provenance {
            source: "synthetic thermal Hex8 hardness".into(), license: "CC0-1.0".into(), artifact: None,
        };
        let mut claims = ClaimSet::new();
        for (name, dims, value) in [
            (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, 1200.0),
            (POISSON_RATIO_PROPERTY, Dims::NONE, 0.3),
            (YIELD_STRESS_PROPERTY, Pressure::DIMS, 1.0),
            (ISOTROPIC_HARDENING_MODULUS_PROPERTY, Pressure::DIMS, 30.0),
            (KINEMATIC_HARDENING_MODULUS_PROPERTY, Pressure::DIMS, 20.0),
            (VOCE_SATURATION_INCREMENT_PROPERTY, Pressure::DIMS, 0.6),
            (VOCE_SATURATION_RATE_PROPERTY, Dims::NONE, 80.0),
            ("fast_C", Pressure::DIMS, 300.0), ("fast_gamma", Dims::NONE, 250.0),
            ("slow_C", Pressure::DIMS, 80.0), ("slow_gamma", Dims::NONE, 30.0),
            ("Cp", CP_DIMS, CP), (DENSITY_PROPERTY, Density::DIMS, 1000.0),
        ] {
            claims.insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims), value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![], provenance: provenance.clone(),
            }).unwrap();
        }
        let observation = claims.register_observation(ObservationDataset {
            specimen: "synthetic thermal Hex8 material".into(),
            method: "synthetic force/projected-area indentation".into(),
            artifact: fs_blake3::hash_domain("fixture", b"thermal Hex8 hardness"),
            caveats: "numerical regression only".into(), provenance: provenance.clone(),
        }).unwrap();
        let keys = [10.0, 20.0].map(|load| {
            PropertyKey::with_quantity("H", INDENTATION_HARDNESS_QUANTITY)
                .with_hardness_test(HardnessTestContext::new(
                    "synthetic spherical indenter, radius 1 mm",
                    vec![HardnessLoadStep::new(
                        QtyAny::new(load, Dims([1, 1, -2, 0, 0, 0])),
                        QtyAny::new(5.0, Dims([0, 0, 1, 0, 0, 0])),
                    ).unwrap()],
                    "synthetic projected-area protocol", observation,
                ).unwrap()).unwrap()
        });
        for (index, key) in keys.iter().enumerate() {
            let samples = [[0.0, 280.0], [0.08, 280.0], [0.08, maximum_hardness_temperature],
                [0.0, maximum_hardness_temperature]].map(|x| {
                    (x, (1.0 + index as f64) * expected_hardness(x[1], x[0]))
                }).to_vec();
            claims.insert_claim(PropertyClaim {
                key: key.clone(), value: PropertyValue::Surface2D {
                    axes: [ALPHA.into(), "T".into()], axis_dims: [Dims::NONE, Temperature::DIMS],
                    samples, triangles: vec![[0, 1, 2], [0, 2, 3]], dims: Pressure::DIMS,
                },
                validity: ValidityDomain::unconstrained().with("T", 280.0, maximum_hardness_temperature)
                    .with_quantity(ALPHA, QuantitySpec::dimensional(Dims::NONE), 0.0, 0.08),
                uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::LinearInside,
                observations: vec![observation], provenance: provenance.clone(),
            }).unwrap();
        }
        let card = MaterialCard::assemble(MaterialStateId {
            chemistry: "synthetic thermoplastic solid".into(), phase: "solid".into(),
            process: "thermal Hex8 hardness fixture".into(), revision: 0,
        }, claims, vec![]).unwrap();
        let caloric = ResolvedHeatCapacityCurve::try_new(&card, HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("Cp", CP_DIMS),
            density: &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS), temperature_axis: "T",
            lower: &QueryPoint::new().with("T", 280.0).unwrap(),
            upper: &QueryPoint::new().with("T", 320.0).unwrap(), selection: SelectionPolicy::SingleClaimOnly,
        }, SolidLiquidPhase::Solid, -20.0 * CP, 1e-8).unwrap();
        Self { card, caloric, keys }
    }
    fn source(&self) -> ThermoplasticSource<'_> {
        ThermoplasticSource::try_new(SourceQuery {
            card: &self.card, caloric: &self.caloric,
            context: QueryPoint::new().with("T", 300.0).unwrap(), temperature_axis: "T".into(),
            selection: SELECT, hardening_basis: HardeningPropertyBasis::MixedPlasticStrainModuli,
            thermal_expansion: None, strain_limit: 0.04,
        }, 0.0).unwrap()
    }
}
fn prescribed(scale: f64) -> Vec<(usize, f64)> {
    NODES.iter().flat_map(|x| [scale * x[0] * (1.0 + 0.4 * x[1]), 0.0, 0.0])
        .enumerate().collect()
}
fn problem<'p, F>(families: &'p [F], fixed: &'p [(usize, f64)], formulation: Formulation)
    -> ThermalHexProblem<'p, F>
{
    ThermalHexProblem {
        nodes_m: &NODES, hexahedra: &CELLS, families, formulation, strain_limit: 0.04,
        rotation_limit_rad: 0.04, prescribed_m: fixed, nodal_forces_n: &[],
        budget: HexAssemblyBudget::default(), settings: IncrementSettings::default(),
    }
}
fn tolerances() -> ThermalHexTolerances {
    ThermalHexTolerances {
        material: EnergyTolerances { energy_j_m3: 1e-8, yield_pa: YIELD_TOL }, energy_j: 1e-8,
    }
}

#[test]
fn g1_actual_point_temperature_and_plastic_history_control_hardness_through_reversal() {
    with_cx(|cx, _| {
        let fixture = Fixture::new(320.0);
        let families = [fixture.source()];
        for formulation in [Formulation::Standard, Formulation::BBar] {
            let fixed = prescribed(0.006);
            let p = problem(&families, &fixed, formulation);
            let temperatures = std::array::from_fn::<_, 8, _>(|q| 295.0 + 2.0 * q as f64);
            let h = temperatures.map(|t| CP * (t - 300.0));
            let initial = vec![families[0].reference().initial_state(); 8];
            let mut step = p.solve_from(ThermalHexHistory {
                displacement_m: &[0.0; 24], quadrature_states: &initial,
                specific_enthalpy_j_kg: &[0.0; 8],
            }, &h, 0.1, tolerances(), cx).unwrap();
            let binding = p.bind_hardness_field(&step, &[(fixture.keys[0].clone(), SELECT)], YIELD_TOL, cx).unwrap();
            let other = p.bind_hardness_field(&step, &[(fixture.keys[1].clone(), SELECT)], YIELD_TOL, cx).unwrap();
            let first = binding.at_accepted(&p, &step, YIELD_TOL, cx).unwrap();
            let second = other.at_accepted(&p, &step, YIELD_TOL, cx).unwrap();
            for q in 0..8 {
                near(second.points()[q].hardness_pa(), 2.0 * first.points()[q].hardness_pa());
                assert_eq!(first.points()[q].key(), &fixture.keys[0]);
            }
            let mut old_field = first;
            for (scale, shift) in [(0.006, 5.0), (-0.006, -5.0)] {
                let fixed = prescribed(scale);
                let p = problem(&families, &fixed, formulation);
                let h = temperatures.map(|t| CP * (t + shift - 300.0));
                let candidate = p.solve_from(step.history(), &h, 0.1, tolerances(), cx).unwrap();
                let before = candidate.clone();
                let field = binding.at_accepted(&p, &candidate, YIELD_TOL, cx).unwrap();
                for q in 0..8 {
                    let alpha = candidate.solution.quadrature_states[q].alpha;
                    near(field.points()[q].hardness_pa(), expected_hardness(temperatures[q] + shift, alpha));
                    assert_eq!(field.points()[q].resolved().state_coordinate(ALPHA), Some(alpha));
                    assert_eq!(field.points()[q].resolved().state_coordinate("T"), Some(temperatures[q] + shift));
                    if shift > 0.0 {
                        near(alpha, step.solution.quadrature_states[q].alpha);
                        near(field.points()[q].hardness_pa(), old_field.points()[q].hardness_pa() - 15.0);
                    } else {
                        assert!(alpha > step.solution.quadrature_states[q].alpha);
                        assert!(field.points()[q].hardness_pa() > old_field.points()[q].hardness_pa());
                    }
                }
                assert_eq!(candidate.solution, before.solution);
                assert_eq!(candidate.energy, before.energy);
                assert_eq!(candidate.quadrature_energy, before.quadrature_energy);
                assert!(field.hardness_pa().last().unwrap() != field.hardness_pa().next().unwrap());
                step = candidate;
                old_field = field;
            }
        }
    });
}

#[test]
fn g0_late_hardness_coverage_and_changed_owner_refuse_without_partial_native_changes() {
    with_cx(|cx, gate| {
        let fixture = Fixture::new(310.0);
        let families = [fixture.source()];
        let fixed = prescribed(0.006);
        let mut p = problem(&families, &fixed, Formulation::BBar);
        let initial = vec![families[0].reference().initial_state(); 8];
        let step = p.solve_from(ThermalHexHistory {
            displacement_m: &[0.0; 24], quadrature_states: &initial,
            specific_enthalpy_j_kg: &[0.0; 8],
        }, &[0.0; 8], 0.1, tolerances(), cx).unwrap();
        let binding = p.bind_hardness_field(&step, &[(fixture.keys[0].clone(), SELECT)], YIELD_TOL, cx).unwrap();
        let accepted = binding.at_accepted(&p, &step, YIELD_TOL, cx).unwrap();
        let mut hot_h = [0.0; 8];
        hot_h[7] = 15.0 * CP;
        let candidate = p.solve_from(step.history(), &hot_h, 0.1, tolerances(), cx).unwrap();
        let before = candidate.clone();
        assert!(matches!(binding.at_accepted(&p, &candidate, YIELD_TOL, cx),
            Err(ThermalHexHardnessError::Source { point: 7, .. })));
        assert_eq!(candidate.solution, before.solution);
        assert_eq!(candidate.energy, before.energy);
        assert_eq!(binding.at_accepted(&p, &step, YIELD_TOL, cx).unwrap(), accepted);
        let other_families = [fixture.source()];
        let foreign = problem(&other_families, &fixed, Formulation::BBar);
        assert!(matches!(binding.at_accepted(&foreign, &step, YIELD_TOL, cx),
            Err(ThermalHexHardnessError::Invalid(_))));
        let mut wrong_temperature = step.clone();
        wrong_temperature.quadrature_temperature_k[7] += 1.0;
        assert!(matches!(binding.at_accepted(&p, &wrong_temperature, YIELD_TOL, cx),
            Err(ThermalHexHardnessError::Invalid(_))));
        p.budget.maximum_quadrature_points = 7;
        assert!(matches!(binding.at_accepted(&p, &step, YIELD_TOL, cx),
            Err(ThermalHexHardnessError::Invalid(_))));
        p.budget.maximum_quadrature_points = 8;
        gate.request();
        assert!(matches!(binding.at_accepted(&p, &step, YIELD_TOL, cx),
            Err(ThermalHexHardnessError::Cancelled)));
    });
}

#[test]
fn g3_thermal_hex_hardness_keeps_complete_chaboche_histories_and_native_energy() {
    with_cx(|cx, _| {
        let fixture = Fixture::new(320.0);
        let source = fixture.source().with_voce_hardening(SELECT).unwrap()
            .with_chaboche_branches(SELECT, [BackstressProperties::new("fast_C", "fast_gamma"),
                BackstressProperties::new("slow_C", "slow_gamma")]).unwrap();
        let families = [source];
        let initial = families[0].resolve_step(0.0, 0.0, 0.1).unwrap().law.initial_state();
        let fixed = prescribed(0.006);
        let p = problem(&families, &fixed, Formulation::BBar);
        let step = p.solve_from(ThermalHexHistory {
            displacement_m: &[0.0; 24], quadrature_states: &vec![initial; 8],
            specific_enthalpy_j_kg: &[0.0; 8],
        }, &[5.0 * CP; 8], 0.1, tolerances(), cx).unwrap();
        let binding = p.bind_hardness_field(&step, &[(fixture.keys[0].clone(), SELECT)], YIELD_TOL, cx).unwrap();
        let field = binding.at_accepted(&p, &step, YIELD_TOL, cx).unwrap();
        for q in 0..8 {
            let state = &step.solution.quadrature_states[q];
            assert_ne!(state.backstress_pa[0], state.backstress_pa[1]);
            near(field.points()[q].hardness_pa(), expected_hardness(305.0, state.plastic.alpha));
        }
        let mut corrupt = step.clone();
        corrupt.solution.quadrature_states[7].backstress_pa[1][0] += 1.0;
        assert!(matches!(binding.at_accepted(&p, &corrupt, YIELD_TOL, cx),
            Err(ThermalHexHardnessError::Native { point: 7, .. })));
        let fixed = prescribed(-0.006);
        let p = problem(&families, &fixed, Formulation::BBar);
        let reversed = p.solve_from(step.history(), &[-5.0 * CP; 8], 0.1, tolerances(), cx).unwrap();
        let before = reversed.clone();
        let values = binding.at_accepted(&p, &reversed, YIELD_TOL, cx).unwrap();
        for q in 0..8 {
            assert!(reversed.solution.quadrature_states[q].plastic.alpha > step.solution.quadrature_states[q].plastic.alpha);
            near(values.points()[q].hardness_pa(), expected_hardness(295.0, reversed.solution.quadrature_states[q].plastic.alpha));
        }
        assert_eq!(reversed.solution, before.solution);
        assert_eq!(reversed.energy, before.energy);
    });
}
