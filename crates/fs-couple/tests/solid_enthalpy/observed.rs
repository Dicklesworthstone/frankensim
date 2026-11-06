//! Physical heat-owner history reaches the existing material and solid solve.
use super::*;
use fs_thermal::observers::ObserverCouplingError;
use std::collections::BTreeMap;

/// A finite heat store with a one-second decay time, disconnected at stop_s.
/// Only the final native interval determines its heat; step accepts its debit.
#[derive(Clone, Debug, PartialEq)]
struct FiniteHeater {
    remaining_j: f64,
    released_j: f64,
    stop_s: f64,
    pending: Option<(f64, f64, f64, f64)>,
    intervals: Vec<(f64, f64)>,
}
impl FiniteHeater {
    fn new(energy_j: f64, stop_s: f64) -> Self {
        Self { remaining_j: energy_j, released_j: 0.0, stop_s,
            pending: None, intervals: Vec::new() }
    }
}
impl DepositionObserver for FiniteHeater {
    type Error = DepositionError;

    fn requires_step_interval(&self) -> bool { true }

    fn next_transition_time_s(
        &mut self, thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Option<f64>, Self::Error> {
        checkpoint()?;
        Ok((thermal.time_s() < self.stop_s).then_some(self.stop_s))
    }

    fn step_sources(
        &mut self, thermal: &DepositionSimulation<'_>, start_s: f64,
        dt_s: f64, end_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        checkpoint()?;
        if self.pending.is_some() || thermal.time_s() != start_s
            || (start_s < self.stop_s && end_s > self.stop_s) {
            return Err(DepositionError::InvalidInput("heater requires one transition-clipped interval"));
        }
        let heat_j = if start_s < self.stop_s {
            self.remaining_j * -(-dt_s).exp_m1()
        } else { 0.0 };
        self.pending = Some((start_s, dt_s, end_s, heat_j));
        Ok(BTreeMap::from([(CellId(30), heat_j / dt_s)]))
    }

    fn event(&mut self, _: &Event, _: &DepositionSimulation<'_>)
        -> Result<(), Self::Error> { Ok(()) }

    fn step(
        &mut self, start_s: f64, duration_s: f64, _: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        checkpoint()?;
        let (prepared_start, prepared_dt, end_s, heat_j) = self.pending.take()
            .ok_or(DepositionError::InvalidInput("heater interval was not prepared"))?;
        if prepared_start != start_s || prepared_dt != duration_s || end_s != thermal.time_s() {
            return Err(DepositionError::InvalidInput("heater and thermal intervals differ"));
        }
        self.remaining_j -= heat_j;
        self.released_j += heat_j;
        self.intervals.push((start_s, end_s));
        Ok(())
    }
}

#[test]
fn interval_heat_and_retained_source_history_drive_native_solid_equilibrium() {
    with_cx(|cx| {
        let (card, caloric) = material(1000.0, 1e-5, 10.0, 0.25, false);
        let thermal = initial(caloric.curve(), &[(30, 0.0)], 1.0);
        let context = point(300.0);
        let model = SolidEnthalpy::try_new(&thermal, &NODES, &TETS,
            vec![source(30, &card, caloric.curve(), &context)],
            TetAssemblyBudget::standard(), 0.1, config().mesh_frame, cx).unwrap();
        let loads = model.geometry().assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let reference_identity = model.resolve(&thermal, cx).unwrap()
            .materials()[0].material_state_identity;
        let heater = FiniteHeater::new(1000.0, 0.3);
        let original = heater.clone();
        let (warm, warm_history) = model.advance_observed_to(&thermal, &heater,
            0.2, thermal_config(), &loads, &FREE, config(), cx).unwrap();
        let accepted_history = warm_history.clone();
        let (next, history) = model.advance_observed_to(&warm.thermal, &warm_history,
            1.0, thermal_config(), &loads, &FREE, config(), cx).unwrap();

        let heat_j = 1000.0 * -(-0.3_f64).exp_m1();
        let rise = heat_j / 10.0;
        let strain = 1e-5 * (rise + 0.005 * rise * rise);
        near(next.thermal.temperature_k(CellId(30)).unwrap(), 300.0 + rise, 1e-10);
        near(warm.thermal_audit.external_heat_j + next.thermal_audit.external_heat_j,
            heat_j, 1e-9);
        near(warm.thermal_audit.stored_energy_change_j + next.thermal_audit.stored_energy_change_j,
            heat_j, 1e-9);
        near(history.released_j, heat_j, 1e-10);
        near(history.remaining_j + history.released_j, 1000.0, 1e-10);
        assert_eq!(history.pending, None);
        assert!(history.requires_step_interval());
        assert_eq!(history.intervals.len(),
            warm.thermal_audit.substeps + next.thermal_audit.substeps);
        assert!(history.intervals.iter().any(|&(_, end)| end == 0.3));
        assert!(history.intervals.iter().all(|&(start, end)| !(start < 0.3 && end > 0.3)));
        assert_eq!(history.intervals.last().unwrap().1, next.thermal.time_s());
        assert_eq!(next.state.time_s(), next.thermal.time_s());
        assert_eq!(next.thermal.time_s(), 1.0);
        assert_eq!(next.thermal.cells()[&CellId(30)].boundary.source_w, 0.0);
        assert_ne!(next.state.materials()[0].material_state_identity, reference_identity);
        assert_eq!(next.state.samples()[0].elastic.card_identity(), card.content_hash());
        near(next.state.materials()[0].density_kg_m3, 6.0, 1e-12);
        for (node, position) in NODES.iter().enumerate() {
            for (axis, &coordinate) in position.iter().enumerate() {
                near(next.equilibrium.solution.displacement.displacement_m()[node][axis],
                    strain * coordinate, 1e-10);
            }
        }
        for &stress in next.equilibrium.solution.stress.elements()[0].stress_mandel_pa() {
            near(stress, 0.0, 1e-8);
        }

        // Even a failure after the entire heat advance discards the source debit.
        let small = SolidEquilibriumConfig { maximum_nodal_displacement_m: 1e-6, ..config() };
        assert!(matches!(model.advance_observed_to(&thermal, &heater,
            1.0, thermal_config(), &loads, &FREE, small, cx),
            Err(ObserverCouplingError::First(SolidEnthalpyError::Solid(_)))));
        assert_eq!(heater, original);
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(thermal.cells()[&CellId(30)].specific_h_j_kg, 0.0);
        assert_eq!(warm.thermal.time_s(), 0.2);
        assert_eq!(warm_history, accepted_history);
        near(accepted_history.released_j, 1000.0 * -(-0.2_f64).exp_m1(), 1e-10);
    });
}

#[test]
fn observed_source_support_exit_preserves_accepted_heat_and_source_history() {
    with_cx(|cx| {
        let (card, caloric) = material(1000.0, 1e-5, 10.0, 0.25, true);
        let thermal = initial(caloric.curve(), &[(30, 0.0)], 1.0);
        let context = point(300.0);
        let model = SolidEnthalpy::try_new(&thermal, &NODES, &TETS,
            vec![source(30, &card, caloric.curve(), &context)],
            TetAssemblyBudget::standard(), 0.1, config().mesh_frame, cx).unwrap();
        let loads = model.geometry().assemble(&[], VolumeLoad3::None, &[], cx).unwrap();
        let heater = FiniteHeater::new(2000.0, 0.3);
        let (accepted, history) = model.advance_observed_to(&thermal, &heater,
            0.2, thermal_config(), &loads, &FREE, config(), cx).unwrap();
        let before_history = history.clone();
        let before_h = accepted.thermal.cells()[&CellId(30)].specific_h_j_kg;
        let claim = accepted.state.samples()[0].elastic.property("young_modulus")
            .unwrap().answer().receipt.selected;

        // The original selected source ends at 349 K. A different same-name
        // source above 350 K must not be substituted as the heater continues.
        assert!(matches!(model.advance_observed_to(&accepted.thermal, &history,
            1.0, thermal_config(), &loads, &FREE, config(), cx),
            Err(ObserverCouplingError::First(SolidEnthalpyError::Material(_)))));
        assert_eq!(history, before_history);
        assert_eq!(accepted.thermal.time_s(), 0.2);
        assert_eq!(accepted.thermal.cells()[&CellId(30)].specific_h_j_kg, before_h);

        // Retrying a supported endpoint continues the actual accepted source.
        let (retry, history) = model.advance_observed_to(&accepted.thermal, &history,
            0.25, thermal_config(), &loads, &FREE, config(), cx).unwrap();
        near(history.released_j, 2000.0 * -(-0.25_f64).exp_m1(), 1e-10);
        near(retry.thermal.temperature_k(CellId(30)).unwrap(),
            300.0 + history.released_j / 10.0, 1e-10);
        assert_eq!(retry.state.samples()[0].elastic.property("young_modulus")
            .unwrap().answer().receipt.selected, claim);
        assert_eq!(retry.thermal.time_s(), 0.25);
        assert_eq!(heater, FiniteHeater::new(2000.0, 0.3));
        assert_eq!(thermal.time_s(), 0.0);
    });
}
