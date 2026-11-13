//! Original complete material heat observer on the actual inertial history.
use crate::circuit_enthalpy::contacts::heat::PreparedHeatNetwork;
use super::*;

impl DynamicCutThermoelasticSolid3<'_> {
    /// Advance source conduction/contact/reservoir/radiation, real native
    /// inertia and their reversible exchange as the same accepted interval.
    /// The independent observer is retained once; enclosure net roundoff is
    /// excluded from external supply before any accepted owner is replaced.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_heat_network<'t, 'p, 'n, 'a, 'scope, O: DepositionObserver>(
        &mut self, mut native: CutElasticityMut3<'_>,
        thermal: &mut DepositionSimulation<'t>, prepared: &'p PreparedHeatNetwork<'n, 'a>,
        observer: &mut O, target_s: f64, external_n: &[f64],
        thermal_config: StepConfig, config: CutThermoelasticConfig3, cx: &'p Cx<'scope>,
    ) -> Result<DynamicCutThermoelasticReport3, CutThermoelasticHeatError3<O::Error>> {
        super::super::heat::admit(&self.bindings, thermal, prepared, cx)?;
        let seed = prepared.observe(cx, observer.clone(), config.energy_tolerance_j)
            .map_err(CutThermoelasticHeatError3::Heat)?;
        let accepted = self.attempt(&mut native, thermal, &seed, target_s, external_n,
            thermal_config, config, |actual| actual.first.closed_radiation_net_j(), cx)
            .map_err(CutThermoelasticHeatError3::Coupled)?;
        let old_history = core::mem::replace(&mut self.history, accepted.history);
        let old_thermal = core::mem::replace(thermal, accepted.thermal);
        let old_observer = core::mem::replace(observer, accepted.observer.second);
        drop((old_history, old_thermal, old_observer, seed, accepted.observer.first));
        Ok(accepted.report)
    }
}
