//! Supplied-state audits for outer enclosure/convection coupling.
use super::*;

#[test]
fn accepted_enthalpy_enclosure_rechecks_actual_boundary_without_resolving_history() {
    let fixture = Fixture::new(false);
    let source = ScalarField::Uniform(0.0);
    let old = [2000.0, 2000.0, 2000.0, 2000.0, 1000.0, 1000.0, 1000.0, 1000.0];
    let history = old;
    with_cx(|cx| {
        let engine = fixture.engine(cx, false);
        let accepted = engine.advance_with_enclosure_radiation(cx, fixture.problem(&source),
            None, &old, DT, step_config(), &fixture.enclosure, radiation_config()).unwrap();
        let candidate = accepted.conduction.specific_enthalpy_j_kg.clone();
        let exact = engine.evaluate_endpoint_with_enclosure(cx, fixture.problem(&source),
            None, &old, &candidate, DT, step_config(), &fixture.enclosure, radiation_config()).unwrap();
        assert!(exact.residual_norm_j <= accepted.physical_residual_tolerance_j);
        close(exact.energy_residual_j, accepted.physical_energy_residual_j, 1e-10);
        close(exact.convective_out_w, accepted.convective_out_w, 1e-11);
        let changed = fixture.boundary.with_uniform_robin_replacements(&[
            (0, 1.2 * HTC, 310.0), (1, 0.8 * HTC, 290.0),
        ]).unwrap();
        let physical = engine.evaluate_endpoint_with_enclosure(cx,
            ConductionProblem { boundary: &changed, ..fixture.problem(&source) },
            None, &old, &candidate, DT, step_config(), &fixture.enclosure, radiation_config()).unwrap();
        assert!(physical.residual_norm_j > 1e-3,
            "an unchanged solved state must expose its new boundary defect");
        assert_eq!(candidate, accepted.conduction.specific_enthalpy_j_kg);
        assert_eq!(old, history);
        assert_eq!(physical.temperature, exact.temperature);
        let area = 3.0_f64.sqrt();
        let left = physical.temperature[..4].iter().sum::<f64>() / 4.0;
        let right = physical.temperature[4..].iter().sum::<f64>() / 4.0;
        let convection = area * (1.2 * HTC * (left - 310.0) + 0.8 * HTC * (right - 290.0));
        let nodal_mass = RHO * 2.0_f64.sqrt() / 48.0;
        let stored = nodal_mass * candidate.iter().zip(old).map(|(&h, previous)| h - previous).sum::<f64>();
        close(physical.convective_out_w, convection, 1e-10);
        close(physical.stored_energy_change_j, stored, 1e-11);
        close(physical.radiation_out_w, 0.0, 1e-8);
        close(physical.neumann_out_w, 0.0, 0.0);
        close(physical.dirichlet_in_w, 0.0, 0.0);
        close(physical.energy_residual_j, stored + DT * convection, 1e-10);
        for (row, (h, reference)) in physical.convective_robin_fluxes.iter()
            .zip([(1.2 * HTC, 310.0), (0.8 * HTC, 290.0)]) {
            close(row.mean_htc_w_per_m2_k, h, 0.0);
            close(row.mean_reference_temperature_k, reference, 0.0);
            close(row.heat_rate_w, h * area * (row.mean_wall_temperature_k - reference), 1e-11);
        }
    });
}

#[test]
fn enclosure_endpoint_audit_observes_latent_storage_despite_identical_temperature() {
    let fixture = Fixture::new(true);
    let source = ScalarField::Uniform(0.0);
    let old = [2000.0; 8];
    let candidate = [2200.0; 8];
    with_cx(|cx| {
        let engine = fixture.engine(cx, true);
        let physical = engine.evaluate_endpoint_with_enclosure(cx, fixture.problem(&source),
            None, &old, &candidate, DT, step_config(), &fixture.enclosure, radiation_config()).unwrap();
        for temperature in &physical.temperature { close(*temperature, 350.0, 0.0); }
        let stored = 2.0 * RHO * 2.0_f64.sqrt() / 12.0 * 200.0;
        let convection = 2.0 * HTC * 3.0_f64.sqrt() * 50.0;
        close(physical.stored_energy_change_j, stored, 1e-12);
        close(physical.convective_out_w, convection, 1e-12);
        close(physical.radiation_out_w, 0.0, 0.0);
        close(physical.energy_residual_j, stored + DT * convection, 1e-12);
        assert!(physical.residual_norm_j > 1.0,
            "a same-temperature state is not a converged enthalpy state");
        assert_eq!(old, [2000.0; 8]);
        assert_eq!(candidate, [2200.0; 8]);
    });
}
