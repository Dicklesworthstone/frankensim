//! G1/G3: native thermal equilibrium is part of material-update acceptance.
//! Authored source curves exercise mechanics and rollback, not calibration.

use std::cell::Cell;
use std::ops::ControlFlow;
use std::panic::{AssertUnwindSafe, catch_unwind};

use fs_blake3::ContentHash;
use fs_cutfem::elastic3::adaptive::AdaptiveElasticity3;
use fs_cutfem::elastic3::{
    CutElasticity3, ElasticMaterialField3, ElasticityError3, ElasticityFormulation3,
    ElasticityOptions3, MaterialRefreshLimits3,
};
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_evidence::ValidityDomain;
use fs_ivl::Interval;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::elastic::thermal::{
    AffineIsotropicThermoelastic, IsothermalThermoelasticInput, IsothermalThermoelasticState,
};
use fs_material::state_point::{
    LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, MaterialPropertySelection,
    POISSON_RATIO_PROPERTY, ScalarAdmissibility, ScalarPropertyRequirement, YOUNG_MODULUS_PROPERTY,
    resolve_material_state_point,
};
use fs_qty::{
    Dims, Pressure, QuantitySpec, Temperature,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_solver::op::LinearOp;

const POINTS: usize = 100_000;
const TOTAL_STRAIN: f64 = 0.1;

fn poll() -> ControlFlow<()> {
    ControlFlow::Continue(())
}

struct Top;
impl CutSdf3 for Top {
    fn value(&self, p: [f64; 3]) -> f64 {
        p[2] - 0.73
    }

    fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval {
        Interval::new(lo[2] - 0.73, hi[2] - 0.73)
    }

    fn derivative_enclose(&self, _: [f64; 3], _: [f64; 3], axis: HeightAxis) -> Interval {
        let value = if axis == HeightAxis::Z { 1.0 } else { 0.0 };
        Interval::new(value, value)
    }
}

fn control<T>(action: impl FnOnce(&mut QuadratureControl3<'_>) -> T) -> T {
    let mut callback = |_| poll();
    let mut control = QuadratureControl3::new(
        QuadratureOptions3 { depth: 0, ..Default::default() },
        &mut callback,
    )
    .unwrap();
    action(&mut control)
}

fn domain() -> HexCell {
    HexCell::try_new([0.0; 3], [1.0; 3]).unwrap()
}

fn options(formulation: ElasticityFormulation3) -> ElasticityOptions3 {
    ElasticityOptions3 { formulation, ghost_gamma: 2.7, ..Default::default() }
}

fn cartesian(
    field: ElasticMaterialField3<'_>,
    formulation: ElasticityFormulation3,
) -> CutElasticity3 {
    control(|control| {
        CutElasticity3::build_with_material_field_and_embedded_dirichlet(
            domain(), [2; 3], &Top, field, &|p| p[0] == 0.0,
            &|_, normal| normal[2] > 0.5, options(formulation),
            Default::default(), Default::default(), control,
        )
        .unwrap()
    })
}

fn tree() -> Octree3 {
    let tree = Octree3::uniform(1, 4, 4096).unwrap();
    tree.refined(&[*tree.leaves().iter().next().unwrap()], poll).unwrap()
}

fn adaptive(
    field: ElasticMaterialField3<'_>,
    tree: &Octree3,
    formulation: ElasticityFormulation3,
) -> AdaptiveElasticity3 {
    control(|control| {
        AdaptiveElasticity3::build_with_material_field_and_embedded_dirichlet(
            domain(), tree, &Top, field, &|p| p[0] == 0.0,
            &|_, normal| normal[2] > 0.5, options(formulation),
            Default::default(), Default::default(), control,
        )
        .unwrap()
    })
}

// Same source chart and reference inventory as cutfem_thermoelastic: changing
// temperature changes C(T) and integrated expansion, while density stays fixed.
fn law() -> AffineIsotropicThermoelastic {
    let quantity = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature, ValueForm::Static,
    ));
    let properties = [
        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, 1000.0, -2.0),
        (POISSON_RATIO_PROPERTY, Dims::NONE, 0.2, 0.0008),
        (
            LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
            Dims([0, 0, 0, -1, 0, 0]), 0.0001, 0.000002,
        ),
    ];
    let mut claims = ClaimSet::new();
    for (name, dims, value, slope) in properties {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims),
            value: PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![
                    (280.0, value - 20.0 * slope),
                    (300.0, value),
                    (360.0, value + 60.0 * slope),
                ],
                dims,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity("T", quantity, 280.0, 360.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "authored CutFEM material-candidate fixture".into(),
                license: "test-only".into(),
                artifact: None,
            },
        }).unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic thermoelastic source".into(),
            phase: "solid".into(),
            process: "native material acceptance test".into(),
            revision: 0,
        },
        claims,
        vec![],
    ).unwrap();
    let query = QueryPoint::new().with_quantity("T", quantity, 300.0).unwrap();
    let requirements: Vec<_> = properties.into_iter().map(|(name, dims, _, _)| {
        ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::Finite).unwrap()
    }).collect();
    let elastic = resolve_material_state_point(
        &card, &query, &requirements[..2], MaterialPropertySelection::SingleClaimOnly,
    ).unwrap();
    let expansion = resolve_material_state_point(
        &card, &query, &requirements[2..], MaterialPropertySelection::SingleClaimOnly,
    ).unwrap();
    AffineIsotropicThermoelastic::new(
        &card, &elastic, &expansion, [280.0, 360.0], 0.05,
    ).unwrap()
}

fn frozen(law: &AffineIsotropicThermoelastic, temperature: f64) -> IsothermalThermoelasticState {
    law.at_temperature(temperature).unwrap()
        .isothermal_elasticity(IsothermalThermoelasticInput {
            reference_density_kg_m3: 3.0,
            mechanical_material_identity: ContentHash([93; 32]),
        }).unwrap()
}

fn temperature(bounds: HexCell) -> f64 {
    let [x, y, z] = bounds.lo();
    316.0 + 3.0 * x + 2.0 * y + z
}

fn probe(nodes: &[[f64; 3]]) -> Vec<f64> {
    nodes.iter().flat_map(|&[x, y, z]| [
        0.002 * x * (1.0 + y), -0.001 * x * z, 0.0007 * x * x,
    ]).collect()
}

fn action(op: &dyn LinearOp, field: &[f64]) -> Vec<f64> {
    let mut result = vec![0.0; op.n()];
    op.apply(field, &mut result);
    result
}

fn close_vectors(actual: &[f64], expected: &[f64]) {
    assert_eq!(actual.len(), expected.len());
    for (&a, &b) in actual.iter().zip(expected) {
        assert!((a - b).abs() <= 3e-10 + 3e-8 * a.abs().max(b.abs()), "{a:e} != {b:e}");
    }
}

// Expand the same physical experiment against the two existing operator APIs;
// no test-only operator abstraction or substitute numerical solve is introduced.
macro_rules! solve_native {
    ($operator:expr, $external:expr) => {{
        let candidate = $operator;
        let prepared = candidate.prepare_thermoelastic(POINTS, poll).unwrap();
        let mut force = candidate.eigenstrain_load(&prepared, poll).unwrap();
        for (value, external) in force.iter_mut().zip(($external).iter()) {
            *value += *external;
        }
        let solution = candidate.solve_controlled(&force, 1e-11, 20_000, 16, |_| poll()).unwrap();
        assert!(solution.iterations() > 0);
        close_vectors(&action(candidate, solution.coefficients()), &force);
        // Algebraic convergence alone is not acceptance of a physical state.
        let response = candidate.thermoelastic_response(
            solution.coefficients(), TOTAL_STRAIN, POINTS, poll,
        ).unwrap();
        assert!(!response.is_empty());
        (solution, prepared, response)
    }};
}

macro_rules! accept_candidate {
    ($op:ident, $fresh:ident, $hot:ident) => {{
        assert_eq!($op.nodes(), $fresh.nodes());
        let scales: Vec<_> = (0..$op.cells())
            .map(|i| 0.6 + 0.25 * i as f64 / $op.cells() as f64).collect();
        $op.set_scales(&scales).unwrap();
        $fresh.set_scales(&scales).unwrap();
        let nodes = $op.nodes().to_vec();
        let fixed = $op.fixed().to_vec();
        let inventory = $op.reference_cell_masses_kg().unwrap().to_vec();
        let occupancy = vec![1.0; $op.cells()];
        let old_source = $op.prepare_thermoelastic(POINTS, poll).unwrap();
        let old_mass = $op.prepare_material_mass(POINTS, poll).unwrap();
        let old_diagonal = old_mass.diagonal(&occupancy, poll).unwrap();
        let field = probe($op.nodes());
        let before = action(&$op, &field);
        let external = $op.body_load(&|_| [0.02, -0.01, 0.03], poll).unwrap();
        let (reference, reference_source, reference_points) = solve_native!(&$fresh, &external);
        let (solution, source, mass, points) = $op.refresh_material_field_with(
            ElasticMaterialField3::ThermoelasticByCell(&$hot),
            MaterialRefreshLimits3::default(),
            |candidate| {
                assert!(candidate.validate_eigenstrain(&old_source).is_err());
                assert!(candidate.validate_material_mass(&old_mass).is_err());
                let (solution, source, points) = solve_native!(candidate, &external);
                close_vectors(solution.coefficients(), reference.coefficients());
                let mass = candidate.prepare_material_mass(POINTS, poll)?;
                assert_eq!(mass.diagonal(&occupancy, poll)?, old_diagonal);
                Ok((solution, source, mass, points))
            },
            poll,
        ).unwrap();
        $op.validate_eigenstrain(&source).unwrap();
        $op.validate_material_mass(&mass).unwrap();
        assert!($op.validate_eigenstrain(&old_source).is_err());
        assert!($op.validate_material_mass(&old_mass).is_err());
        assert_eq!($op.nodes(), nodes.as_slice());
        assert_eq!($op.fixed(), fixed.as_slice());
        assert_eq!($op.scales(), scales.as_slice());
        assert_eq!($op.reference_cell_masses_kg().unwrap(), inventory.as_slice());
        assert_eq!(mass.diagonal(&occupancy, poll).unwrap(), old_diagonal);
        close_vectors(&action(&$op, &field), &action(&$fresh, &field));
        assert!(action(&$op, &field).iter().zip(&before).any(|(a, b)| (a - b).abs() > 1e-4));
        close_vectors(
            &$op.eigenstrain_load(&source, poll).unwrap(),
            &$fresh.eigenstrain_load(&reference_source, poll).unwrap(),
        );
        assert_eq!(points.len(), reference_points.len());
        assert!(points.iter().any(|p| p.temperature_k != points[0].temperature_k));
        for (actual, expected) in points.iter().zip(&reference_points) {
            assert_eq!(actual.temperature_k, expected.temperature_k);
            close_vectors(&actual.mechanical.stress, &expected.mechanical.stress);
        }
        assert_eq!(
            $op.thermoelastic_response(solution.coefficients(), TOTAL_STRAIN, POINTS, poll).unwrap(),
            points,
        );
    }};
}

#[derive(Clone, Copy, Debug)]
enum Failure {
    Rejected,
    Cancelled,
    AttemptPanic,
    CheckpointPanic,
}

fn failures() -> Vec<Failure> {
    let mut cases = vec![Failure::Rejected, Failure::Cancelled];
    if cfg!(panic = "unwind") {
        cases.extend([Failure::AttemptPanic, Failure::CheckpointPanic]);
    }
    cases
}

macro_rules! reject_candidates {
    ($op:ident, $hot:ident) => {{
        let scales: Vec<_> = (0..$op.cells())
            .map(|i| 0.6 + 0.25 * i as f64 / $op.cells() as f64).collect();
        $op.set_scales(&scales).unwrap();
        let old_source = $op.prepare_thermoelastic(POINTS, poll).unwrap();
        let old_mass = $op.prepare_material_mass(POINTS, poll).unwrap();
        let occupancy = vec![1.0; $op.cells()];
        let old_diagonal = old_mass.diagonal(&occupancy, poll).unwrap();
        let old_load = $op.eigenstrain_load(&old_source, poll).unwrap();
        let old_state = $op.cell_thermoelastic_state(0).unwrap().clone();
        let inventory = $op.reference_cell_masses_kg().unwrap().to_vec();
        let field = probe($op.nodes());
        let before = action(&$op, &field);
        let external = $op.body_load(&|_| [0.02, -0.01, 0.03], poll).unwrap();
        for failure in failures() {
            let completed = Cell::new(false);
            // Caller side effects are intentionally observable: candidate-owned
            // preparations can escape, but must be stale after rollback.
            let mut escaped = None;
            let result = catch_unwind(AssertUnwindSafe(|| {
                $op.refresh_material_field_with::<()>(
                    ElasticMaterialField3::ThermoelasticByCell(&$hot),
                    MaterialRefreshLimits3::default(),
                    |candidate| {
                        assert!(candidate.validate_eigenstrain(&old_source).is_err());
                        assert!(candidate.validate_material_mass(&old_mass).is_err());
                        assert_ne!(action(candidate, &field), before);
                        let (_, source, points) = solve_native!(candidate, &external);
                        assert!(points.iter().all(|p| p.temperature_k > 300.0));
                        let mass = candidate.prepare_material_mass(POINTS, poll)?;
                        escaped = Some((source, mass));
                        completed.set(true);
                        match failure {
                            Failure::Rejected => Err(ElasticityError3::Invalid("candidate objective refused")),
                            Failure::AttemptPanic => panic!("candidate physical check panicked"),
                            _ => Ok(()),
                        }
                    },
                    || {
                        if completed.get() {
                            match failure {
                                Failure::Cancelled => return ControlFlow::Break(()),
                                Failure::CheckpointPanic => panic!("final material acceptance checkpoint panicked"),
                                _ => {}
                            }
                        }
                        poll()
                    },
                )
            }));
            assert!(completed.get(), "failure must follow a completed native solve: {failure:?}");
            match failure {
                Failure::Rejected => assert!(matches!(
                    result, Ok(Err(ElasticityError3::Invalid("candidate objective refused")))
                )),
                Failure::Cancelled => assert!(matches!(
                    result, Ok(Err(ElasticityError3::Cancelled))
                )),
                Failure::AttemptPanic | Failure::CheckpointPanic => assert!(result.is_err()),
            }
            $op.validate_eigenstrain(&old_source).unwrap();
            $op.validate_material_mass(&old_mass).unwrap();
            let (source, mass) = escaped.expect("candidate preparation escaped before rejection");
            assert!($op.validate_eigenstrain(&source).is_err());
            assert!($op.validate_material_mass(&mass).is_err());
            assert_eq!(action(&$op, &field), before);
            assert_eq!($op.eigenstrain_load(&old_source, poll).unwrap(), old_load);
            assert_eq!(old_mass.diagonal(&occupancy, poll).unwrap(), old_diagonal);
            assert_eq!($op.scales(), scales.as_slice());
            assert_eq!($op.reference_cell_masses_kg().unwrap(), inventory.as_slice());
            for cell in 0..$op.cells() {
                assert_eq!($op.cell_thermoelastic_state(cell).unwrap(), &old_state);
            }
        }
    }};
}

#[test]
fn accepted_native_temperature_solve_returns_preparations_for_the_committed_operator() {
    let law = law();
    let cold = frozen(&law, 300.0);
    let hot = |bounds| Ok(frozen(&law, temperature(bounds)));
    for formulation in [ElasticityFormulation3::Displacement, ElasticityFormulation3::BBar] {
        let mut op = cartesian(ElasticMaterialField3::ThermoelasticUniform(&cold), formulation);
        let mut fresh = cartesian(ElasticMaterialField3::ThermoelasticByCell(&hot), formulation);
        accept_candidate!(op, fresh, hot);
    }
}

#[test]
fn rejected_cancelled_and_unwinding_native_solves_restore_the_original_thermal_owner() {
    let law = law();
    let cold = frozen(&law, 300.0);
    let hot = |bounds| Ok(frozen(&law, temperature(bounds)));
    let mut op = cartesian(
        ElasticMaterialField3::ThermoelasticUniform(&cold), ElasticityFormulation3::BBar,
    );
    reject_candidates!(op, hot);
}

#[test]
fn adaptive_temperature_acceptance_keeps_hanging_reduction_and_reference_mass() {
    let law = law();
    let cold = frozen(&law, 300.0);
    let hot = |bounds| Ok(frozen(&law, temperature(bounds)));
    let tree = tree();
    for formulation in [ElasticityFormulation3::Displacement, ElasticityFormulation3::BBar] {
        let mut op = adaptive(ElasticMaterialField3::ThermoelasticUniform(&cold), &tree, formulation);
        let mut fresh = adaptive(ElasticMaterialField3::ThermoelasticByCell(&hot), &tree, formulation);
        assert!(op.physical_nodes().len() > op.nodes().len());
        let physical_nodes = op.physical_nodes().to_vec();
        let leaves = op.leaves().to_vec();
        let field = probe(op.nodes());
        let reconstructed = op.physical_displacements(&field).unwrap();
        accept_candidate!(op, fresh, hot);
        assert_eq!(op.physical_nodes(), physical_nodes.as_slice());
        assert_eq!(op.leaves(), leaves.as_slice());
        assert_eq!(op.physical_displacements(&field).unwrap(), reconstructed);
    }
}

#[test]
fn adaptive_rejection_and_unwind_restore_native_bindings_after_reduced_solves() {
    let law = law();
    let cold = frozen(&law, 300.0);
    let hot = |bounds| Ok(frozen(&law, temperature(bounds)));
    let tree = tree();
    let mut op = adaptive(
        ElasticMaterialField3::ThermoelasticUniform(&cold), &tree, ElasticityFormulation3::BBar,
    );
    assert!(op.physical_nodes().len() > op.nodes().len());
    let physical_nodes = op.physical_nodes().to_vec();
    let nodes = op.nodes().to_vec();
    let fixed = op.fixed().to_vec();
    let field = probe(op.nodes());
    let reconstructed = op.physical_displacements(&field).unwrap();
    reject_candidates!(op, hot);
    assert_eq!(op.nodes(), nodes.as_slice());
    assert_eq!(op.physical_nodes(), physical_nodes.as_slice());
    assert_eq!(op.fixed(), fixed.as_slice());
    assert_eq!(op.physical_displacements(&field).unwrap(), reconstructed);
}

#[path = "cutfem_material_candidate/decision.rs"]
mod decision;

#[path = "cutfem_material_candidate/state.rs"]
mod state;
