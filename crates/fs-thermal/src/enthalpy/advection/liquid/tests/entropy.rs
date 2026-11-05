//! G0/G1/G3: moving-liquid entropy and stationary latent storage remain distinct.
use super::*;
use crate::enthalpy::DepositionObserver;
use crate::enthalpy::advection::entropy::AdvectionEntropyAudit;
use std::collections::BTreeMap;

fn balance(a: AdvectionEntropyAudit) {
    close(a.thermal_storage_change_j_k,
        a.advective_boundary_entropy_input_j_k + a.upwind_mixing_j_k
            + a.temporal_defect_j_k + a.other_energy_at_cell_temperature_j_k
            + a.advective_balance_residual_j_k, 2e-13);
}

#[test]
fn isothermal_circulation_and_flushing_do_not_carry_stationary_phase_entropy() {
    for datum in [0.0, 1_099_511_627_776.0, -1_099_511_627_776.0] {
        let law = PhaseLaw { datum, identity: 7, melting_range_k: 0.0 };
        let initial = prepare(&law, &[20.0, 50.0, 80.0], 2.0);
        for rate in [-0.5, 0.5] {
            let original = Advected::new((), flows(3, rate), 3).unwrap();
            for selected in [
                original.clone().with_entropy_audit().with_liquid_carrier(carrier(&law)).unwrap(),
                original.clone().with_liquid_carrier(carrier(&law)).unwrap().with_entropy_audit(),
            ] {
                let (next, observed, energy) = initial.advance_observed(
                    &selected, 1.0, config(1.0), || Ok(()),
                ).unwrap();
                assert_eq!(next.cells(), initial.cells());
                assert_eq!(energy.external_heat_j, 0.0);
                assert_eq!(observed.entropy_audit(), Some(AdvectionEntropyAudit::default()));
                assert_eq!(selected.audit().substeps, 0);
            }
        }
        let single = prepare(&law, &[20.0], 2.0);
        let flushing = Advected::new((), vec![], 0).unwrap()
            .with_boundaries(open::boundaries(CellId(0), CellId(0), 0.5, datum + 100.0), 2).unwrap()
            .with_liquid_carrier(carrier(&law)).unwrap().with_entropy_audit();
        let (same, observed, _) = single.advance_observed(&flushing, 1.0, config(1.0), || Ok(())).unwrap();
        assert_eq!(same.cells(), single.cells());
        assert_eq!(observed.entropy_audit(), Some(AdvectionEntropyAudit::default()));
    }
}

#[test]
fn superheat_melts_solid_with_physical_carrier_entropy_and_separate_time_error() {
    let law = PhaseLaw { datum: 0.0, identity: 7, melting_range_k: 0.0 };
    let initial = prepare(&law, &[150.0, 25.0], 2.0);
    let original = Advected::new((), flows(2, 0.5), 2).unwrap()
        .with_liquid_carrier(carrier(&law)).unwrap();
    let (ordinary, _, original_energy) = initial.advance_observed(&original, 0.125, config(1.0), || Ok(())).unwrap();
    let (next, observed, energy) = initial.advance_observed(
        &original.with_entropy_audit(), 0.125, config(1.0), || Ok(()),
    ).unwrap();
    assert_eq!(next.cells(), ordinary.cells());
    assert_eq!(energy, original_energy);
    assert_eq!(observed.audit().substeps, 1);
    let transferred_mass = 0.5 * 0.125;
    let heat = transferred_mass * 50.0;
    let hot_temperature = 525.0 - 0.5 * heat;
    close(next.cells()[&CellId(0)].specific_h_j_kg, 150.0 - heat, 0.0);
    close(next.cells()[&CellId(1)].specific_h_j_kg, 25.0 + heat, 0.0);
    assert_eq!(next.temperature_k(CellId(1)).unwrap(), 500.0);
    let audit = observed.entropy_audit().unwrap();
    let hot_entropy_change = 2.0 * (hot_temperature / 525.0_f64).ln();
    close(audit.thermal_storage_change_j_k, hot_entropy_change + heat / 500.0, 1e-15);
    close(audit.upwind_mixing_j_k, heat * (1.0 / 500.0 - 1.0 / 525.0), 1e-15);
    close(audit.temporal_defect_j_k, hot_entropy_change + heat / 525.0, 1e-15);
    close(audit.other_energy_at_cell_temperature_j_k, 0.0, 1e-15);
    close(audit.advective_balance_residual_j_k, 0.0, 1e-15);
    assert_eq!(audit.advective_boundary_entropy_input_j_k, 0.0);
    assert!(audit.upwind_mixing_j_k > 0.0 && audit.temporal_defect_j_k < 0.0);
    balance(audit);
}

#[test]
fn liquid_inlet_entropy_closes_latent_storage_and_is_caloric_datum_invariant() {
    let mut audits = Vec::new();
    for datum in [0.0, 1_099_511_627_776.0, -1_099_511_627_776.0] {
        let law = PhaseLaw { datum, identity: 7, melting_range_k: 0.0 };
        let initial = prepare(&law, &[20.0], 2.0);
        let original = Advected::new((), vec![], 0).unwrap()
            .with_liquid_carrier(carrier(&law)).unwrap()
            .with_boundaries(open::boundaries(CellId(0), CellId(0), 0.5, datum + 200.0), 2).unwrap();
        let baseline = initial.advance_observed(&original, 0.25, config(1.0), || Ok(())).unwrap();
        let (next, observed, energy) = initial.advance_observed(
            &original.with_entropy_audit(), 0.25, config(1.0), || Ok(()),
        ).unwrap();
        assert_eq!(next.cells(), baseline.0.cells());
        assert_eq!(energy, baseline.2);
        assert_eq!(observed.audit().substeps, 1);
        assert_eq!(next.temperature_k(CellId(0)).unwrap(), 500.0);
        assert_eq!(energy.external_heat_j, 12.5);
        let audit = observed.entropy_audit().unwrap();
        // A 550 K liquid stream discharges at 500 K: Cp=2 and dm=1/8.
        // Stationary solid contributes storage L*df/Tm, never stream entropy.
        let net_carried_entropy = 0.125 * 2.0 * 1.1_f64.ln();
        close(audit.thermal_storage_change_j_k, 12.5 / 500.0, 1e-15);
        close(audit.advective_boundary_entropy_input_j_k, net_carried_entropy, 1e-15);
        close(audit.upwind_mixing_j_k, 12.5 / 500.0 - net_carried_entropy, 1e-15);
        close(audit.temporal_defect_j_k, 0.0, 1e-15);
        close(audit.other_energy_at_cell_temperature_j_k, 0.0, 1e-15);
        close(audit.advective_balance_residual_j_k, 0.0, 1e-15);
        balance(audit);
        audits.push(audit);
    }
    assert_eq!(audits[0], audits[1]);
    assert_eq!(audits[0], audits[2]);
}

#[derive(Clone)]
struct Cooling;
impl DepositionObserver for Cooling {
    type Error = DepositionError;
    fn sources(&mut self, thermal: &DepositionSimulation<'_>,
        _: &mut dyn FnMut() -> Result<(), DepositionError>)
        -> Result<BTreeMap<CellId, f64>, Self::Error> {
        Ok(thermal.cells().keys().map(|&id| (id, -80.0)).collect())
    }
    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>) -> Result<(), Self::Error> { Ok(()) }
    fn step(&mut self, _: f64, _: f64, _: &BTreeMap<CellId, f64>,
        _: &DepositionSimulation<'_>, _: &mut dyn FnMut() -> Result<(), DepositionError>)
        -> Result<(), Self::Error> { Ok(()) }
}

#[test]
fn freezing_uses_starting_liquid_entropy_and_retains_separate_source_cooling() {
    let law = PhaseLaw { datum: 0.0, identity: 7, melting_range_k: 0.0 };
    let initial = prepare(&law, &[20.0; 2], 2.0);
    let flow = Advected::new(Cooling, flows(2, 0.5), 2).unwrap()
        .with_liquid_carrier(carrier(&law)).unwrap().with_entropy_audit();
    let (frozen, observed, energy) = initial.advance_observed(&flow, 0.25, config(1.0), || Ok(())).unwrap();
    assert!(frozen.cells().values().all(|cell| cell.specific_h_j_kg == 0.0));
    assert_eq!(observed.audit().substeps, 1);
    assert_eq!(energy.external_heat_j, -40.0);
    let audit = observed.entropy_audit().unwrap();
    close(audit.thermal_storage_change_j_k, -40.0 / 500.0, 1e-15);
    close(audit.other_energy_at_cell_temperature_j_k, -40.0 / 500.0, 1e-15);
    assert_eq!(audit.upwind_mixing_j_k, 0.0);
    assert_eq!(audit.advective_boundary_entropy_input_j_k, 0.0);
    assert_eq!(audit.advective_balance_residual_j_k, 0.0);
    close(audit.temporal_defect_j_k, 0.0, 1e-15);
    balance(audit);
    assert!(matches!(frozen.advance_observed(&observed, 0.5, config(1.0), || Ok(())),
        Err(DepositionError::InvalidInput("nonzero liquid mass flow crosses a fully solid cell"))));
    assert_eq!(observed.entropy_audit(), Some(audit));
    assert_eq!(initial.time_s(), 0.0);

    // No flow through a solid needs no mobile carrier at all; storage can cool.
    let stationary = Advected::new(Cooling, flows(2, 0.0), 2).unwrap()
        .with_liquid_carrier(carrier(&law)).unwrap().with_entropy_audit();
    let (_, observed, _) = frozen.advance_observed(&stationary, 0.5, config(1.0), || Ok(())).unwrap();
    balance(observed.entropy_audit().unwrap());
}

#[test]
fn fully_liquid_entropy_matches_bulk_audit_and_cancelled_advance_is_atomic() {
    let law = PhaseLaw { datum: 0.0, identity: 7, melting_range_k: 0.0 };
    let initial = prepare(&law, &[130.0, 170.0, 250.0], 2.0);
    let bulk = Advected::new((), flows(3, 0.5), 3).unwrap().with_entropy_audit();
    let liquid = bulk.clone().with_liquid_carrier(carrier(&law)).unwrap();
    let (reference, reference_observer, reference_energy) = initial.advance_observed(
        &bulk, 0.25, config(1.0), || Ok(()),
    ).unwrap();
    let (next, observed, energy) = initial.advance_observed(&liquid, 0.25, config(1.0), || Ok(())).unwrap();
    assert_eq!(next.cells(), reference.cells());
    assert_eq!(observed.entropy_audit(), reference_observer.entropy_audit());
    assert_eq!(energy, reference_energy);
    let mut polls = 0;
    assert!(matches!(initial.advance_observed(&liquid, 1.0, config(0.01), || {
        polls += 1;
        if polls > 60 { Err(DepositionError::Cancelled) } else { Ok(()) }
    }), Err(DepositionError::Cancelled)));
    assert_eq!(initial.time_s(), 0.0);
    assert_eq!(liquid.audit().substeps, 0);
    assert_eq!(liquid.entropy_audit(), Some(AdvectionEntropyAudit::default()));
}
