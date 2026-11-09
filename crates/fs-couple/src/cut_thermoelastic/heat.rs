//! Complete existing material heat transport inside the same native acceptance.
use crate::circuit_enthalpy::contacts::heat::{HeatNetworkError, PreparedHeatNetwork};
use fs_thermal::observers::ObserverCouplingError;
use super::*;

#[derive(Debug)]
pub enum CutThermoelasticHeatError3<E = Infallible> {
    Heat(HeatNetworkError),
    Coupled(CutThermoelasticError3<ObserverCouplingError<HeatNetworkError, E>>),
    Invalid(&'static str),
    Cancelled,
}
impl<E: core::fmt::Debug> core::fmt::Display for CutThermoelasticHeatError3<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "material heat/native cut solid: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for CutThermoelasticHeatError3<E> {}

impl CutThermoelasticSolid3<'_> {
    /// Advance the SAME finite caloric cells with the original source-backed
    /// conduction/contact/reservoir/radiation observer and compatible CutFEM
    /// equilibrium. Electrical networks refuse because this owner does not
    /// represent circuit work.
    ///
    /// Internal enclosure radiation's integrated net roundoff is removed from
    /// external supply INSIDE the native candidate acceptance. Both heat and
    /// combined accounts use the same original joule tolerance. No second heat
    /// step, source debit or observer composition occurs after acceptance.
    /// Pass the independent observer again on continuation; this method retains
    /// its accepted history and reconstructs only the native heat wrapper.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_heat_network<'t, 'p, 'n, 'a, 'scope, O: DepositionObserver>(
        &mut self, mut native: CutElasticityMut3<'_>,
        thermal: &mut DepositionSimulation<'t>, prepared: &'p PreparedHeatNetwork<'n, 'a>,
        observer: &mut O, target_s: f64, nodal_forces_n: &[f64],
        thermal_config: StepConfig, config: CutThermoelasticConfig3, cx: &'p Cx<'scope>,
    ) -> Result<CutThermoelasticReport3, CutThermoelasticHeatError3<O::Error>> {
        admit(&self.bindings, thermal, prepared, cx)?;
        let seed = prepared.observe(cx, observer.clone(), config.energy_tolerance_j)
            .map_err(CutThermoelasticHeatError3::Heat)?;
        let accepted = self.attempt(&mut native, thermal, &seed, target_s, nodal_forces_n,
            thermal_config, config, |actual| actual.first.closed_radiation_net_j(), cx)
            .map_err(CutThermoelasticHeatError3::Coupled)?;
        // Keep ALL observer-owning values alive until all accepted owners move.
        // In particular seed.second is a user observer with a possible Drop.
        let old_history = core::mem::replace(&mut self.history, accepted.history);
        let old_thermal = core::mem::replace(thermal, accepted.thermal);
        let old_observer = core::mem::replace(observer, accepted.observer.second);
        drop((old_history, old_thermal, old_observer, seed, accepted.observer.first));
        Ok(accepted.report)
    }

}

pub(super) fn admit<E>(
        bindings: &Bindings<'_>, thermal: &DepositionSimulation<'_>, prepared: &PreparedHeatNetwork<'_, '_>,
        cx: &Cx<'_>,
    ) -> Result<(), CutThermoelasticHeatError3<E>> {
        cx.checkpoint().map_err(|_| CutThermoelasticHeatError3::Cancelled)?;
        let network = prepared.thermal_material_network().map_err(CutThermoelasticHeatError3::Heat)?;
        for expected in &bindings.cells {
            cx.checkpoint().map_err(|_| CutThermoelasticHeatError3::Cancelled)?;
            let actual = network.cells.iter().find(|cell| cell.cell == expected.cell)
                .ok_or(CutThermoelasticHeatError3::Invalid("native cut-solid cell lacks its original heat-network binding"))?;
            if expected.reference.material_id != actual.material_id
                || !std::ptr::eq(expected.curve, actual.curve)
                || expected.point.material.content_hash() != actual.point.material.content_hash()
                || expected.point.temperature_axis != actual.point.temperature_axis
                || expected.point.liquid_fraction_axis != actual.point.liquid_fraction_axis
            {
                return Err(CutThermoelasticHeatError3::Invalid("heat and native cut solid require identical material and caloric owners"));
            }
            let cell = thermal.cells().get(&expected.cell)
                .ok_or(CutThermoelasticHeatError3::Invalid("inactive native cut-solid thermal cell"))?;
            let phase = expected.curve.state_at_specific_enthalpy(cell.specific_h_j_kg)
                .map_err(|error| CutThermoelasticHeatError3::Coupled(CutThermoelasticError3::Phase(error)))?;
            let mechanical = expected.point.at(phase)
                .map_err(|error| CutThermoelasticHeatError3::Coupled(CutThermoelasticError3::State(error)))?;
            let transport = actual.point.at(phase)
                .map_err(|error| CutThermoelasticHeatError3::Coupled(CutThermoelasticError3::State(error)))?;
            if mechanical != transport {
                return Err(CutThermoelasticHeatError3::Invalid("heat and native cut solid material conditions disagree"));
            }
        }
        Ok(())
    }
