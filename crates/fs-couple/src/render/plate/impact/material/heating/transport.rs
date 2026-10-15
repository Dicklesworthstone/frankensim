//! Contact heat is an internal source in the existing uniform thermal carrier.
//! A transport owner proposes all stores before either physical owner commits.
use super::*;
use fs_material::phase::{UniformEnthalpyStep,UniformEnthalpyStepInput};

/// Accepted heat-transfer records, retained with their concrete solver reports.
#[derive(Clone, Debug)]
pub struct TransportedImpactFrame<R> {
    /// Combined accounting; its residual also subtracts the boundary heat below.
    pub impact: HeatedImpactFrame,
    /// Actual signed boundary heat supplied to all thermal stores this step [J].
    pub boundary_heat_j: f64,
    /// Original thermal proposals/reports in reservoir order, without erasure.
    pub transfers: Vec<UniformEnthalpyStep<R>>,
}
impl<'a> HeatedImpactSystem<'a> {
    /// Couple deposited contact work to an existing thermal transport owner.
    ///
    /// Every store is visited, even during separation or elastic contact: cooling
    /// does not stop when contact stops. The request uses the current chart state,
    /// invariant mass/geometry, actual dt and only that store's contact heat.
    /// Pass e.g. `|_, input| environment.advance(cx,input).map_err(...)` for
    /// fs-conduction, or another owner of UniformEnthalpyStep. The callback must
    /// propose, not debit live reservoirs: irreversible external side effects
    /// cannot be rolled back. Original reports publish only with complete success.
    ///
    /// A returned state must belong to the SAME chart. Per-step and cumulative
    /// physical energy defects must fit the existing absolute heat allowance;
    /// a callback cannot enlarge it by returning a larger residual tolerance.
    /// No correction rescales heat or repairs a solver state. Adiabatic stepping
    /// can resume from the actual returned state without replaying boundary heat.
    /// This reporting path allocates bounded proposal storage (at most 64 stores).
    pub fn step_with_transport<R>(&mut self, forces:&[f64], gate:&CancelGate,
        mut transport:impl FnMut(&str,UniformEnthalpyStepInput<'_>)
            -> Result<UniformEnthalpyStep<R>,ImpactError>,
    ) -> Result<TransportedImpactFrame<R>,ImpactError> {
        let Self {mechanics,stores,routes,heat,heat_tolerance_j,..}=self;
        let mut proposals=Vec::with_capacity(stores.len());
        let mut inventories=Vec::with_capacity(stores.len());
        let (mut deposited,mut change,mut boundary)=(0.0,0.0,0.0);
        let mechanical=mechanics.step_checked(forces,gate,|system,next,velocity,record| {
            deposited=stage_heat(system,next,velocity,record,routes,heat)?;
            for (i,s) in stores.iter().enumerate() {
                if gate.is_requested() {return Err(ImpactError::Cancelled);}
                let p=transport(&s.key,UniformEnthalpyStepInput {
                    curve:s.curve,initial:s.state,mass_kg:s.mass_kg,volume_m3:s.volume_m3,
                    surface_area_m2:s.surface_area_m2,internal_heat_j:heat[i],duration_s:system.config.dt_s,
                })?;
                let h=p.state.specific_enthalpy_j_kg();
                if !p.external_heat_j.is_finite() || !p.energy_residual_tolerance_j.is_finite()
                    || p.energy_residual_tolerance_j<0.0 || p.energy_residual_tolerance_j>*heat_tolerance_j
                    || s.curve.state_at_specific_enthalpy(h).map_err(|e|ImpactError::Owner(e.to_string()))?!=p.state {
                    return Err(invalid("thermal transport returned foreign material state or invalid source/budget"));
                }
                let delta=s.mass_kg*(h-s.state.specific_enthalpy_j_kg());
                let (received,received_correction)=accumulate(s.received_heat_j,s.received_correction_j,heat[i]);
                let (external,external_correction)=accumulate(s.boundary_heat_j,s.boundary_correction_j,p.external_heat_j);
                let inventory=s.mass_kg*(h-s.initial.specific_enthalpy_j_kg());
                let local=delta-heat[i]-p.external_heat_j;
                let cumulative=inventory-received-external;
                if ![delta,received,received_correction,external,external_correction,inventory,local,cumulative].into_iter().all(f64::is_finite)
                    || local.abs()>*heat_tolerance_j || cumulative.abs()>*heat_tolerance_j {
                    return Err(invalid("thermal transport violates local or accumulated energy balance"));
                }
                change+=delta;boundary+=p.external_heat_j;
                proposals.push(p);inventories.push((received,received_correction,external,external_correction));
            }
            let defect=change-deposited-boundary;
            if ![change,boundary,defect].into_iter().all(f64::is_finite) || defect.abs()>*heat_tolerance_j {
                return Err(invalid("combined thermal transport energy balance exceeds its allowance"));
            }
            Ok(())
        })?;
        // The unchanged mechanical gate has now accepted; publication cannot fail.
        for ((s,p),&(received,rc,external,ec)) in stores.iter_mut().zip(&proposals).zip(&inventories) {
            s.state=p.state;s.received_heat_j=received;s.received_correction_j=rc;
            s.boundary_heat_j=external;s.boundary_correction_j=ec;
            // Preserve the actual solver endpoint, including its disclosed residual.
            s.adiabatic_base=p.state;s.adiabatic_heat_j=received;
        }
        Ok(TransportedImpactFrame {boundary_heat_j:boundary,transfers:proposals,
            impact:HeatedImpactFrame {
                deposited_heat_j:deposited,thermal_energy_change_j:change,
                remaining_mechanical_loss_j:mechanical.dissipated_energy_j-deposited,
                combined_energy_residual_j:mechanical.balance_residual_j+change-deposited-boundary,
                mechanical,
            }})
    }
}
fn accumulate(sum:f64,correction:f64,value:f64)->(f64,f64) {
    let y=value-correction;let next=sum+y;(next,(next-sum)-y)
}
