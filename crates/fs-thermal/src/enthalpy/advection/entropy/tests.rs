//! G1/G3 balances on the original upwind trajectory, not a second integrator.
use super::*;
use crate::enthalpy::*;

struct Law {
    shift: f64,
}
impl ThermalLaw for Law {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let t = h - self.shift;
        if t.is_finite() && t > 0.0 {
            Ok(t)
        } else {
            Err(DepositionError::MaterialDomain)
        }
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0
    }
    fn specific_entropy_change_j_kg_k(
        &self,
        a: f64,
        b: f64,
    ) -> Result<Option<f64>, DepositionError> {
        let ta = self.temperature_k(a)?;
        self.temperature_k(b)?;
        Ok(Some(((b - a) / ta).ln_1p()))
    }
}
static LAW: Law = Law { shift: 0.0 };
fn cfg(dt: f64) -> StepConfig {
    StepConfig {
        max_dt_s: dt,
        energy_tolerance_j: 1e-8,
    }
}
fn cell(id: u64, h: f64) -> DepositedCell {
    DepositedCell {
        id: CellId(id),
        segment_id: id,
        body_id: 1,
        material_id: MaterialId(1),
        mass_kg: 1.0,
        specific_h_j_kg: h,
        geometry: CellGeometry {
            start_m: [id as f64, 0.0, 0.0],
            end_m: [id as f64 + 1.0, 0.0, 0.0],
            volume_m3: 0.001,
            surface_area_m2: 1.0,
        },
        boundary: Boundary {
            ambient_k: 300.0,
            htc_w_m2_k: 0.0,
            reservoir: None,
            source_w: 0.0,
        },
    }
}
fn prepared<'a>(law: &'a dyn ThermalLaw, cells: Vec<DepositedCell>) -> DepositionSimulation<'a> {
    let events = cells
        .into_iter()
        .map(|c| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(c),
        })
        .collect();
    DepositionSimulation::new(
        events,
        &[(MaterialId(1), law)],
        DepositionBudget {
            max_cells: 8,
            max_links: 8,
            max_events: 16,
            max_substeps: 10_000,
        },
    )
    .unwrap()
    .advance_to(0.0, cfg(1.0), || Ok(()))
    .unwrap()
    .0
}
fn ring(rate: f64) -> Vec<MassFlow> {
    (0..3)
        .map(|id| MassFlow {
            id,
            a: CellId(id),
            b: CellId((id + 1) % 3),
            kg_s: rate,
        })
        .collect()
}
fn open(inlet_h: f64) -> Advected<()> {
    Advected::new((), vec![], 0)
        .unwrap()
        .with_boundaries(
            vec![
                MassBoundary::Inlet {
                    id: 10,
                    cell: CellId(0),
                    kg_s: 1.0,
                    material_id: MaterialId(1),
                    specific_h_j_kg: inlet_h,
                },
                MassBoundary::Outlet {
                    id: 11,
                    cell: CellId(0),
                    kg_s: 1.0,
                },
            ],
            2,
        )
        .unwrap()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a} != {b}");
}
fn balance(a: AdvectionEntropyAudit) {
    close(
        a.thermal_storage_change_j_k,
        a.advective_boundary_entropy_input_j_k
            + a.upwind_mixing_j_k
            + a.temporal_defect_j_k
            + a.other_energy_at_cell_temperature_j_k
            + a.advective_balance_residual_j_k,
        2e-12,
    );
}

#[test]
fn closed_ring_retains_original_trajectory_and_separates_upwind_from_time_error() {
    let initial = [300.0, 600.0, 450.0];
    let s = prepared(&LAW, (0..3).map(|i| cell(i as u64, initial[i])).collect());
    let plain = Advected::new((), ring(0.25), 3).unwrap();
    let requested = plain.clone().with_entropy_audit();
    let (ordinary, _, energy) = s
        .advance_observed(&plain, 0.125, cfg(1.0), || Ok(()))
        .unwrap();
    let (next, observed, audited_energy) = s
        .advance_observed(&requested, 0.125, cfg(1.0), || Ok(()))
        .unwrap();
    assert_eq!(next.cells(), ordinary.cells());
    assert_eq!(energy, audited_energy);
    let a = observed.entropy_audit().unwrap();
    let exact: f64 = (0..3)
        .map(|i| (next.temperature_k(CellId(i as u64)).unwrap() / initial[i]).ln())
        .sum();
    let mixing: f64 = (0..3)
        .map(|i| {
            let ratio = initial[(i + 2) % 3] / initial[i];
            0.125 * 0.25 * (ratio - 1.0 - ratio.ln())
        })
        .sum();
    close(a.thermal_storage_change_j_k, exact, 2e-12);
    close(a.upwind_mixing_j_k, mixing, 2e-13);
    assert!(a.upwind_mixing_j_k > 0.0 && a.temporal_defect_j_k < 0.0);
    assert_eq!(a.advective_boundary_entropy_input_j_k, 0.0);
    close(a.other_energy_at_cell_temperature_j_k, 0.0, 1e-12);
    close(a.advective_balance_residual_j_k, 0.0, 1e-12);
    balance(a);
    assert_eq!(
        requested.entropy_audit(),
        Some(AdvectionEntropyAudit::default())
    );
    assert!(plain.entropy_audit().is_none());
    let reversed: Vec<_> = ring(0.25)
        .into_iter()
        .rev()
        .map(|f| MassFlow {
            a: f.b,
            b: f.a,
            kg_s: -f.kg_s,
            ..f
        })
        .collect();
    let reversed = Advected::new((), reversed, 3).unwrap().with_entropy_audit();
    let (same, other, _) = s
        .advance_observed(&reversed, 0.125, cfg(1.0), || Ok(()))
        .unwrap();
    assert_eq!(next.cells(), same.cells());
    assert_eq!(observed.entropy_audit(), other.entropy_audit());
}

#[test]
fn steady_open_flow_uses_donor_entropy_and_keeps_source_heat_unclassified() {
    let mut reports = Vec::new();
    for shift in [0.0, 1_048_576.0, 2.0_f64.powi(52), -2.0_f64.powi(52)] {
        let law = Law { shift };
        let mut c = cell(0, 400.0 + shift);
        c.boundary.source_w = -200.0; // Holds the well-mixed cell at 400 K.
        let s = prepared(&law, vec![c]);
        let flow = open(600.0 + shift).with_entropy_audit();
        let (next, observed, energy, unified) = s
            .advance_advected_with_entropy(&flow, 0.25, cfg(1.0), || Ok(()))
            .unwrap();
        assert_eq!(next.cells(), s.cells());
        close(energy.external_heat_j, 0.0, 0.0);
        let a = observed.entropy_audit().unwrap();
        close(
            a.advective_boundary_entropy_input_j_k,
            0.25 * 1.5_f64.ln(),
            1e-14,
        );
        close(a.upwind_mixing_j_k, 0.25 * (0.5 - 1.5_f64.ln()), 1e-14);
        close(a.other_energy_at_cell_temperature_j_k, -0.125, 0.0);
        assert_eq!(a.thermal_storage_change_j_k, 0.0);
        assert_eq!(a.temporal_defect_j_k, 0.0);
        let advective = unified.advection.unwrap();
        assert_eq!(advective.boundary_entropy_input_j_k, a.advective_boundary_entropy_input_j_k);
        assert_eq!(advective.upwind_mixing_j_k, a.upwind_mixing_j_k);
        assert_eq!(advective.balance_residual_j_k, a.advective_balance_residual_j_k);
        assert_eq!(unified.unclassified_source_at_cell_j_k, -0.125);
        assert!(unified.has_unclassified_sources);
        close(unified.temporal_defect_j_k, 0.0, 1e-14);
        balance(a);
        reports.push((a, unified));
    }
    assert!(reports.iter().all(|report| report == &reports[0]));
}

#[test]
fn unified_entropy_combines_original_streams_and_actual_baths_without_changing_heat() {
    let initial = [300.0, 600.0, 450.0];
    let mut cells: Vec<_> = (0..3).map(|i| cell(i as u64, initial[i])).collect();
    cells[0].boundary.ambient_k = 280.0;
    cells[0].boundary.htc_w_m2_k = 2.0;
    cells[2].boundary.reservoir = Some(Reservoir {
        area_m2: 0.5, conductance_w_k: 1.0, temperature_k: 500.0,
    });
    let s = prepared(&LAW, cells);
    let flow = Advected::new((), ring(0.25), 3).unwrap();
    let (plain, plain_flow, plain_energy) = s
        .advance_observed(&flow, 0.125, cfg(1.0), || Ok(())).unwrap();
    let (next, observed, energy, unified) = s
        .advance_advected_with_entropy(&flow, 0.125, cfg(1.0), || Ok(())).unwrap();
    assert_eq!(next.cells(), plain.cells());
    assert_eq!(energy, plain_energy);
    assert_eq!(observed.audit(), plain_flow.audit());
    assert_eq!(energy.external_heat_j, 1.25);
    assert_eq!(unified.unclassified_source_at_cell_j_k, 0.0);
    assert!(!unified.has_unclassified_sources);
    let advective = unified.advection.unwrap();
    assert_eq!(advective.boundary_entropy_input_j_k, 0.0);
    let mixing: f64 = (0..3).map(|i| {
        let ratio = initial[(i + 2) % 3] / initial[i];
        0.125 * 0.25 * (ratio - 1.0 - ratio.ln())
    }).sum();
    close(advective.upwind_mixing_j_k, mixing, 1e-14);
    close(unified.reservoir_entropy_input_j_k, 0.125 * (-40.0 / 280.0 + 50.0 / 500.0), 1e-14);
    close(unified.conduction_generation_j_k,
        0.125 * ((-40.0 / 300.0) * (-20.0 / 280.0) + (50.0 / 450.0) * (50.0 / 500.0)), 1e-14);
    let exact: f64 = (0..3).map(|i|
        (next.temperature_k(CellId(i as u64)).unwrap() / initial[i]).ln()).sum();
    close(unified.thermal_storage_change_j_k, exact, 1e-14);
    close(unified.temporal_defect_j_k, observed.entropy_audit().unwrap().temporal_defect_j_k, 1e-14);
    assert!(unified.temporal_defect_j_k < 0.0);
    assert!(flow.entropy_audit().is_none());
}

#[test]
fn unified_entropy_preserves_inner_history_events_and_cancelling_source_flags() {
    #[derive(Clone, Debug, PartialEq)]
    struct History(Vec<(u8, f64)>);
    impl DepositionObserver for History {
        type Error = DepositionError;
        fn begin(&mut self) -> Result<(), Self::Error> {
            self.0.push((0, 0.0));
            Ok(())
        }
        fn sources(&mut self, thermal: &DepositionSimulation<'_>,
            _: &mut dyn FnMut() -> Result<(), DepositionError>,
        ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
            self.0.push((1, thermal.time_s()));
            Ok(BTreeMap::from([(CellId(0), 4.0)]))
        }
        fn requires_step_interval(&self) -> bool { true }
        fn step_sources(&mut self, _: &DepositionSimulation<'_>, _: f64, _: f64, end: f64,
            _: &mut dyn FnMut() -> Result<(), DepositionError>,
        ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
            self.0.push((2, end));
            Ok(BTreeMap::from([(CellId(0), -4.0)]))
        }
        fn event(&mut self, _: &Event, thermal: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
            self.0.push((3, thermal.time_s()));
            Ok(())
        }
        fn step(&mut self, _: f64, _: f64, _: &BTreeMap<CellId, f64>,
            thermal: &DepositionSimulation<'_>, _: &mut dyn FnMut() -> Result<(), DepositionError>,
        ) -> Result<(), Self::Error> {
            self.0.push((4, thermal.time_s()));
            Ok(())
        }
    }
    let c = cell(0, 400.0);
    let boundary = Boundary { htc_w_m2_k: 1.0, ..c.boundary };
    let s = DepositionSimulation::new(vec![
        TimedEvent { time_s: 0.0, event: Event::Birth(c) },
        TimedEvent { time_s: 0.125, event: Event::Boundary(CellId(0), boundary) },
    ], &[(MaterialId(1), &LAW as &dyn ThermalLaw)], DepositionBudget {
        max_cells: 1, max_links: 2, max_events: 2, max_substeps: 100,
    }).unwrap();
    let flow = open(600.0).map_inner(|()| History(vec![(9, -1.0)]));
    let (plain, plain_flow, plain_energy) = s
        .advance_observed(&flow, 0.25, cfg(0.2), || Ok(())).unwrap();
    let (next, observed, energy, entropy) = s
        .advance_advected_with_entropy(&flow, 0.25, cfg(0.2), || Ok(())).unwrap();
    assert_eq!(next.cells(), plain.cells());
    assert_eq!(energy, plain_energy);
    assert_eq!(observed.inner(), plain_flow.inner());
    assert_eq!(observed.audit(), plain_flow.audit());
    assert_eq!(entropy.unclassified_source_at_cell_j_k, 0.0);
    assert!(entropy.has_unclassified_sources);
    assert_eq!(energy.events, 2);
    let (continued, retained, _, _) = next
        .advance_advected_with_entropy(&observed, 0.5, cfg(0.2), || Ok(())).unwrap();
    let (expected, expected_history, _) = plain
        .advance_observed(&plain_flow, 0.5, cfg(0.2), || Ok(())).unwrap();
    assert_eq!(continued.cells(), expected.cells());
    assert_eq!(retained.inner(), expected_history.inner());
    assert_eq!(flow.inner().0, vec![(9, -1.0)]);
}

#[test]
fn entropy_uses_the_same_datum_safe_calorimeter_flux_as_energy() {
    let mut reports = Vec::new();
    for shift in [0.0, 2.0_f64.powi(52), -2.0_f64.powi(52)] {
        let law = Law { shift };
        let mut c = cell(0, 300.0 + shift);
        c.boundary.source_w = 4.0 / 3.0;
        let s = prepared(&law, vec![c]);
        let flow = Advected::new((), vec![], 0).unwrap().with_boundaries(vec![
            MassBoundary::Inlet { id: 10, cell: CellId(0), kg_s: 1.0 / 3.0,
                material_id: MaterialId(1), specific_h_j_kg: 400.0 + shift },
            MassBoundary::Outlet { id: 11, cell: CellId(0), kg_s: 1.0 / 3.0 },
        ], 2).unwrap().with_entropy_audit();
        let (next, observed, energy) = s.advance_observed(&flow, 0.75, cfg(2.0), || Ok(())).unwrap();
        assert_eq!(next.temperature_k(CellId(0)).unwrap(), 326.0);
        close(energy.external_heat_j, 26.0, 1e-13);
        let a = observed.entropy_audit().unwrap();
        close(a.thermal_storage_change_j_k, (326.0_f64 / 300.0).ln(), 1e-14);
        close(a.advective_boundary_entropy_input_j_k, 0.25 * (4.0_f64 / 3.0).ln(), 1e-14);
        close(a.other_energy_at_cell_temperature_j_k, 1.0 / 300.0, 1e-14);
        close(a.advective_balance_residual_j_k, 0.0, 1e-14);
        balance(a);
        reports.push(a);
    }
    assert_eq!(reports[0], reports[1]);
    assert_eq!(reports[0], reports[2]);
}

#[test]
fn closed_loop_time_refinement_reduces_temporal_defect_without_a_spatial_claim() {
    let initial = [300.0, 600.0, 450.0];
    let s = prepared(&LAW, (0..3).map(|i| cell(i as u64, initial[i])).collect());
    let observer = Advected::new((), ring(0.5), 3)
        .unwrap()
        .with_entropy_audit();
    let mut errors = Vec::new();
    for dt in [0.125, 0.0625, 0.03125] {
        let (next, observed, _) = s
            .advance_observed(&observer, 2.0, cfg(dt), || Ok(()))
            .unwrap();
        let a = observed.entropy_audit().unwrap();
        let exact: f64 = (0..3)
            .map(|i| (next.temperature_k(CellId(i as u64)).unwrap() / initial[i]).ln())
            .sum();
        close(a.thermal_storage_change_j_k, exact, 2e-12);
        close(next.energy_j(), s.energy_j(), 1e-9);
        assert!(a.upwind_mixing_j_k > 0.0);
        assert!(a.temporal_defect_j_k < 0.0);
        errors.push(-a.temporal_defect_j_k);
        balance(a);
    }
    for pair in errors.windows(2) {
        assert!((1.8..2.2).contains(&(pair[0] / pair[1])), "{errors:?}");
    }
}

#[test]
fn latent_enthalpy_carries_entropy_once_without_temperature_mixing() {
    struct Plateau;
    impl ThermalLaw for Plateau {
        fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
            if h.is_finite() && (0.0..=1e6).contains(&h) {
                Ok(350.0)
            } else {
                Err(DepositionError::MaterialDomain)
            }
        }
        fn maximum_temperature_slope(&self) -> f64 {
            0.0
        }
        fn specific_entropy_change_j_kg_k(
            &self,
            a: f64,
            b: f64,
        ) -> Result<Option<f64>, DepositionError> {
            self.temperature_k(a)?;
            self.temperature_k(b)?;
            Ok(Some((b - a) / 350.0))
        }
    }
    let s = prepared(&Plateau, vec![cell(0, 100_000.0)]);
    let flow = open(200_000.0).with_entropy_audit();
    let (next, observed, energy, unified) = s
        .advance_advected_with_entropy(&flow, 0.125, cfg(1.0), || Ok(()))
        .unwrap();
    let a = observed.entropy_audit().unwrap();
    close(next.temperature_k(CellId(0)).unwrap(), 350.0, 0.0);
    close(energy.external_heat_j, 12_500.0, 0.0);
    close(a.thermal_storage_change_j_k, 12_500.0 / 350.0, 1e-12);
    close(
        a.advective_boundary_entropy_input_j_k,
        12_500.0 / 350.0,
        1e-12,
    );
    assert_eq!(a.upwind_mixing_j_k, 0.0);
    assert_eq!(a.temporal_defect_j_k, 0.0);
    close(a.other_energy_at_cell_temperature_j_k, 0.0, 1e-12);
    close(unified.thermal_storage_change_j_k, a.thermal_storage_change_j_k, 1e-12);
    assert_eq!(unified.advection.unwrap().upwind_mixing_j_k, 0.0);
    assert_eq!(unified.unclassified_source_at_cell_j_k, 0.0);
    assert!(!unified.has_unclassified_sources);
    close(unified.temporal_defect_j_k, 0.0, 1e-12);
    balance(a);
}

#[test]
fn missing_entropy_and_late_cancellation_leave_original_observers_unchanged() {
    struct Missing;
    impl ThermalLaw for Missing {
        fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
            LAW.temperature_k(h)
        }
        fn maximum_temperature_slope(&self) -> f64 {
            1.0
        }
    }
    let cells = || vec![cell(0, 300.0), cell(1, 600.0), cell(2, 450.0)];
    let missing = prepared(&Missing, cells());
    let plain = Advected::new((), ring(0.5), 3).unwrap();
    assert!(
        missing
            .advance_observed(&plain, 0.125, cfg(1.0), || Ok(()))
            .is_ok()
    );
    let flow = plain.with_entropy_audit();
    assert!(matches!(
        missing.advance_advected_with_entropy(&flow, 0.125, cfg(1.0), || Ok(())),
        Err(DepositionError::InvalidInput(
            "thermal law does not provide advective entropy change"
        ))
    ));
    let s = prepared(&LAW, cells());
    let mut polls = 0;
    assert!(matches!(
        s.advance_advected_with_entropy(&flow, 1.0, cfg(0.01), || {
            polls += 1;
            if polls > 60 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(DepositionError::Cancelled)
    ));
    assert_eq!(s.time_s(), 0.0);
    assert_eq!(flow.audit().substeps, 0);
    assert_eq!(flow.entropy_audit(), Some(AdvectionEntropyAudit::default()));
}

#[test]
fn a_shared_enthalpy_reference_does_not_admit_stationary_species_entropy() {
    struct Bound {
        names: Vec<String>,
        fractions: [f64; 1],
    }
    impl ThermalLaw for Bound {
        fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
            LAW.temperature_k(h)
        }
        fn maximum_temperature_slope(&self) -> f64 {
            1.0
        }
        fn specific_entropy_change_j_kg_k(
            &self,
            a: f64,
            b: f64,
        ) -> Result<Option<f64>, DepositionError> {
            LAW.specific_entropy_change_j_kg_k(a, b)
        }
        fn advected_composition(&self) -> Option<super::super::species::CaloricComposition<'_>> {
            Some(super::super::species::CaloricComposition {
                species: &self.names,
                mass_fractions: &self.fractions,
                enthalpy_reference: "synthetic-energy-reference",
            })
        }
    }
    let law = Bound {
        names: vec!["A".into()],
        fractions: [1.0],
    };
    let s = prepared(&law, vec![cell(0, 300.0)]);
    let flow = open(600.0)
        .with_stationary_species(
            vec!["A".into()],
            vec![super::super::species::InletSpecies {
                boundary_id: 10,
                mass_fractions: vec![1.0],
            }],
            1,
            1,
        )
        .unwrap();
    assert!(
        s.advance_observed(&flow, 0.125, cfg(1.0), || Ok(()))
            .is_ok()
    );
    let audited = flow.with_entropy_audit();
    assert!(matches!(
        s.advance_advected_with_entropy(&audited, 0.125, cfg(1.0), || Ok(())),
        Err(DepositionError::InvalidInput(
            "stationary species entropy requires original source entropy anchors"
        ))
    ));
    assert_eq!(s.time_s(), 0.0);
}
