//! Valid rejected thermoelastic proposals remain owned coupling scratch.
use super::*;

use fs_cutfem::elastic3::MaterialRefreshDecision3;
use fs_cutfem::elastic3::surface::ReferenceLoad3;
use fs_cutfem::elastic3::thermal::ThermoelasticEquilibriumSettings3;

#[derive(Debug)]
enum TrialError {
    Native(ElasticityError3),
    Coupling,
}

impl From<ElasticityError3> for TrialError {
    fn from(error: ElasticityError3) -> Self {
        Self::Native(error)
    }
}

fn settings(physical_dofs: usize) -> ThermoelasticEquilibriumSettings3 {
    ThermoelasticEquilibriumSettings3 {
        refresh: MaterialRefreshLimits3::default(),
        maximum_native_dofs: physical_dofs,
        relative_residual_tolerance: 1e-11,
        maximum_iterations: 20_000,
        poll_iterations: 16,
        maximum_free_residual_n: 2e-8,
        maximum_kinematic_strain: TOTAL_STRAIN,
        maximum_points: POINTS,
        maximum_mass_terms: POINTS,
        maximum_prestress_entries: Some(POINTS),
        moment_origin_m: [0.0; 3],
    }
}

macro_rules! decision_roundtrip {
    ($op:ident, $fresh:ident, $hot:ident, $physical_dofs:expr) => {{
        let settings = settings($physical_dofs);
        let scales: Vec<_> = (0..$op.cells())
            .map(|i| 0.7 + 0.02 * (i % 4) as f64).collect();
        $op.set_scales(&scales).unwrap();
        $fresh.set_scales(&scales).unwrap();
        let field = probe($op.nodes());
        let before = action(&$op, &field);
        let protected = $op.prepare_thermoelastic(POINTS, poll).unwrap();
        let body = |_| [0.02, -0.01, 0.03];
        let original = $op.solve_thermoelastic_equilibrium(
            ReferenceLoad3::body(&body), None, None, settings, |_| poll(),
        ).unwrap();
        // A solve on the current source neither rotates its token nor refreshes C.
        $op.validate_eigenstrain(&protected).unwrap();
        assert_eq!(action(&$op, &field), before);
        let occupancy = vec![1.0; $op.cells()];
        let original_mass = original.mass.diagonal(&occupancy, poll).unwrap();
        let external = $fresh.body_load(&body, poll).unwrap();
        let (reference, _, _) = solve_native!(&$fresh, &external);

        let rejected = $op.refresh_material_field_with_decision::<_, TrialError>(
            ElasticMaterialField3::ThermoelasticByCell(&$hot),
            MaterialRefreshLimits3::default(),
            |candidate| {
                assert!(candidate.validate_eigenstrain(&protected).is_err());
                let result = candidate.solve_thermoelastic_equilibrium(
                    ReferenceLoad3::body(&body), None, None, settings, |_| poll(),
                )?;
                close_vectors(result.solution.coefficients(), reference.coefficients());
                Ok(MaterialRefreshDecision3::Reject(result))
            },
            poll,
        ).unwrap();
        let MaterialRefreshDecision3::Reject(rejected) = rejected else {
            panic!("valid rejection must retain its owned equilibrium");
        };
        assert!(rejected.equilibrium.free_residual_norm <= settings.maximum_free_residual_n);
        $op.validate_eigenstrain(&protected).unwrap();
        $op.validate_material_mass(&original.mass).unwrap();
        $op.validate_prestress(original.prestress.as_ref().unwrap()).unwrap();
        assert!($op.validate_eigenstrain(&rejected.eigenstrain).is_err());
        assert!($op.validate_material_mass(&rejected.mass).is_err());
        assert!($op.validate_prestress(rejected.prestress.as_ref().unwrap()).is_err());
        assert_eq!(action(&$op, &field), before);

        let accepted = $op.refresh_material_field_with_decision::<_, TrialError>(
            ElasticMaterialField3::ThermoelasticByCell(&$hot),
            MaterialRefreshLimits3::default(),
            |candidate| {
                let result = candidate.solve_thermoelastic_equilibrium(
                    ReferenceLoad3::body(&body), None, None, settings, |_| poll(),
                )?;
                Ok(MaterialRefreshDecision3::Accept(result))
            },
            poll,
        ).unwrap();
        let MaterialRefreshDecision3::Accept(accepted) = accepted else {
            panic!("valid acceptance must keep its candidate owner");
        };
        close_vectors(accepted.solution.coefficients(), rejected.solution.coefficients());
        close_vectors(&action(&$op, &field), &action(&$fresh, &field));
        assert_ne!(action(&$op, &field), before);
        $op.validate_eigenstrain(&accepted.eigenstrain).unwrap();
        $op.validate_material_mass(&accepted.mass).unwrap();
        $op.validate_prestress(accepted.prestress.as_ref().unwrap()).unwrap();
        assert!($op.validate_eigenstrain(&protected).is_err());
        assert!($op.validate_material_mass(&rejected.mass).is_err());
        assert_eq!(accepted.mass.diagonal(&occupancy, poll).unwrap(), original_mass);
        assert_eq!($op.scales(), scales.as_slice());
    }};
}

#[test]
fn valid_reject_returns_candidate_equilibrium_then_accepts_identical_native_retry() {
    let law = law();
    let cold = frozen(&law, 300.0);
    let hot = |bounds| Ok(frozen(&law, temperature(bounds)));
    let formulation = ElasticityFormulation3::BBar;
    let mut op = cartesian(ElasticMaterialField3::ThermoelasticUniform(&cold), formulation);
    let mut fresh = cartesian(ElasticMaterialField3::ThermoelasticByCell(&hot), formulation);
    decision_roundtrip!(op, fresh, hot, op.n());
}

#[test]
fn adaptive_decisions_keep_original_hanging_space_and_reference_mass() {
    let law = law();
    let cold = frozen(&law, 300.0);
    let hot = |bounds| Ok(frozen(&law, temperature(bounds)));
    let tree = tree();
    let formulation = ElasticityFormulation3::BBar;
    let mut op = adaptive(ElasticMaterialField3::ThermoelasticUniform(&cold), &tree, formulation);
    let mut fresh = adaptive(ElasticMaterialField3::ThermoelasticByCell(&hot), &tree, formulation);
    assert!(op.physical_nodes().len() > op.nodes().len());
    let physical_nodes = op.physical_nodes().to_vec();
    let field = probe(op.nodes());
    let physical = op.physical_displacements(&field).unwrap();
    decision_roundtrip!(op, fresh, hot, 3 * op.physical_nodes().len());
    assert_eq!(op.physical_nodes(), physical_nodes.as_slice());
    assert_eq!(op.physical_displacements(&field).unwrap(), physical);
}

#[test]
fn decision_errors_keep_coupling_type_and_convert_final_cancellation() {
    let law = law();
    let cold = frozen(&law, 300.0);
    let hot = |bounds| Ok(frozen(&law, temperature(bounds)));
    let mut op = cartesian(
        ElasticMaterialField3::ThermoelasticUniform(&cold), ElasticityFormulation3::Displacement,
    );
    let protected = op.prepare_thermoelastic(POINTS, poll).unwrap();
    let field = probe(op.nodes());
    let before = action(&op, &field);
    let failed = op.refresh_material_field_with_decision::<usize, TrialError>(
        ElasticMaterialField3::ThermoelasticByCell(&hot),
        MaterialRefreshLimits3::default(),
        |candidate| {
            assert_ne!(action(candidate, &field), before);
            Err(TrialError::Coupling)
        },
        poll,
    );
    assert!(matches!(failed, Err(TrialError::Coupling)));
    op.validate_eigenstrain(&protected).unwrap();
    assert_eq!(action(&op, &field), before);
    for accept in [false, true] {
        let completed = Cell::new(false);
        let stopped = op.refresh_material_field_with_decision::<usize, TrialError>(
            ElasticMaterialField3::ThermoelasticByCell(&hot),
            MaterialRefreshLimits3::default(),
            |candidate| {
                assert_ne!(action(candidate, &field), before);
                completed.set(true);
                Ok(if accept {
                    MaterialRefreshDecision3::Accept(17)
                } else {
                    MaterialRefreshDecision3::Reject(17)
                })
            },
            || if completed.get() { ControlFlow::Break(()) } else { poll() },
        );
        assert!(matches!(stopped, Err(TrialError::Native(ElasticityError3::Cancelled))));
        op.validate_eigenstrain(&protected).unwrap();
        assert_eq!(action(&op, &field), before);
    }
}
