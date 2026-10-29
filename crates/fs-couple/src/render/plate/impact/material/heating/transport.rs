//! Selected physical loss is an internal source in the existing thermal carrier.
//! A transport owner proposes all stores before either physical owner commits.
use super::*;
use fs_contact::normal_patch::PointNormalPatchReceipt;
use fs_material::phase::{UniformEnthalpyStep,UniformEnthalpyStepInput};

/// One accepted reservoir and its actual contact work for a joint thermal solve.
/// Keys and order are stable; the canonical input owns all physical coordinates.
#[derive(Clone, Copy, Debug)]
pub struct ReservoirThermalInput<'a> {
    pub key: &'a str,
    pub thermal: UniformEnthalpyStepInput<'a>,
}

/// Read-only, unpublished collision endpoint and simultaneous thermal requests.
/// The mechanics are frozen during transport: this is not temperature feedback
/// into stiffness. Observations are instantaneous endpoints, not averaged contact
/// area, continuous collision detection, or an extra source of contact heat.
pub struct ImpactThermalTrial<'a> {
    system: &'a ImpactSystem,
    next: &'a [f64],
    inputs: &'a [ReservoirThermalInput<'a>],
}
impl<'a> ImpactThermalTrial<'a> {
    /// All stores at the SAME old thermal time, in original reservoir order.
    pub fn reservoirs(&self) -> &[ReservoirThermalInput<'a>] { self.inputs }
    /// Exact mechanical duration shared by every reservoir request.
    pub fn duration_s(&self) -> f64 { self.system.config.dt_s }
    /// One-based proposed mechanical sample; this has not yet been published.
    pub fn sample(&self) -> u64 { self.system.sample+1 }
    /// Original immutable contact source, including both ordered bulk states
    /// and the complete interface card. Observing neither mutates nor advances it.
    pub fn contact_source(&self,key:&str)
        ->Result<&fs_contact::normal_patch::source::SourcePointContact,ImpactError> {
        self.system.material_contacts().iter().find(|c|c.key()==key)
            .map(|c|c.source()).ok_or_else(||invalid("unknown thermal-trial contact"))
    }
    /// Physical area of an admitted analytic circular/elliptic contact [m2].
    /// Separation returns exactly zero. Elliptic axes, not an equivalent radius,
    /// supply the area. A finite-gap pressure-equivalent footprint is NOT admitted
    /// as a physical transfer area. This remains an endpoint, not a time average.
    pub fn contact_patch_area_m2(&self,key:&str)->Result<f64,ImpactError> {
        let p=self.contact_observation(key)?;
        if p.finite_gap_response_identity.is_some() {
            return Err(invalid("pressure-equivalent footprint is not a physical thermal contact area"));
        }
        if p.approach_m<=0.0 {return Ok(0.0);}
        let product=match p.elliptic_patch_axes {
            Some(axes)=>axes.semi_major_axis_m*axes.semi_minor_axis_m,
            None=>p.patch_radius_m*p.patch_radius_m,
        };
        let area=std::f64::consts::PI*product;
        if !positive(area) {return Err(invalid("thermal contact area is not representable"));}
        Ok(area)
    }
    /// Original material owner at the proposed mechanical endpoint. Use its
    /// footprint/pressure to resolve a thermal interface, not to redeposit loss.
    /// Source envelopes and geometry are checked by the original contact owner.
    pub fn contact_observation(&self, key:&str) -> Result<PointNormalPatchReceipt,ImpactError> {
        let (contact,c)=self.system.material_contacts().iter().enumerate()
            .find(|(_,c)|c.key()==key).ok_or_else(||invalid("unknown thermal-trial contact"))?;
        let (d,v)=c.kinematics(self.next,self.next)
            .map_err(|source|ImpactError::MaterialContact {contact,source})?;
        c.source().evaluate(&format!("{key}/thermal-trial/{}",self.sample()),d,v,self.duration_s())
            .map_err(|source|ImpactError::MaterialContact {contact,source})
    }
}

/// Joint thermal proposal, including independently declared outside heat.
/// For an isolated collection this is zero even when its stores exchange heat.
#[derive(Clone, Debug)]
pub struct CoupledThermalStep<R> {
    pub transfers: Vec<UniformEnthalpyStep<R>>,
    /// Heat from outside the ENTIRE collection [J], excluding all internal links.
    pub boundary_heat_j: f64,
}

/// Accepted heat-transfer records, retained with their concrete solver reports.
#[derive(Clone, Debug)]
pub struct TransportedImpactFrame<R> {
    /// Combined accounting; its residual also subtracts the boundary heat below.
    pub impact: HeatedImpactFrame,
    /// Actual signed boundary heat supplied to all thermal stores this step [J].
    /// Internal exchange cancels in this sum; individual stores retain both signs.
    pub boundary_heat_j: f64,
    /// Original thermal proposals/reports in reservoir order, without erasure.
    pub transfers: Vec<UniformEnthalpyStep<R>>,
}
impl<'a> HeatedImpactSystem<'a> {
    /// Couple deposited contact work to independent thermal transport owners.
    /// Every store is visited, including after separation. For interacting stores
    /// use `step_with_coupled_transport`: sequential callbacks must not debit a
    /// live partner or separately replay its heat. This convenience path shares
    /// all chart, energy, cancellation and publication gates with that joint path.
    pub fn step_with_transport<R>(&mut self, forces:&[f64], gate:&CancelGate,
        mut transport:impl FnMut(&str,UniformEnthalpyStepInput<'_>)
            -> Result<UniformEnthalpyStep<R>,ImpactError>,
    ) -> Result<TransportedImpactFrame<R>,ImpactError> {
        self.step_with_coupled_transport(forces,gate,|trial| {
            let mut proposals=Vec::with_capacity(trial.inputs.len());
            for input in trial.inputs {
                if gate.is_requested() {return Err(ImpactError::Cancelled);}
                proposals.push(transport(input.key,input.thermal)?);
            }
            let boundary_heat_j=proposals.iter().map(|p|p.external_heat_j).sum();
            Ok(CoupledThermalStep {transfers:proposals,boundary_heat_j})
        })
    }

    /// Propose interacting finite reservoirs together on one collision step.
    ///
    /// The callback is called once with ALL old thermal states, each store's
    /// actual deposited contact work, and the unpublished mechanical endpoint.
    /// Return exactly one canonical proposal per input, in the same order. A
    /// finite exchange Q is -Q in its donor's `external_heat_j` and +Q in its
    /// recipient's; it is NOT counted again as internal contact dissipation.
    /// The signed sum must match the independently declared collection boundary
    /// heat. In an isolated pair that declaration is zero, so a missing donor
    /// debit refuses even when the recipient's local balance would pass.
    ///
    /// Reuse a thermal network/contact solver in this callback. It must propose,
    /// not mutate live reservoirs: arbitrary callback side effects cannot be
    /// undone. Concrete reports publish only with complete mechanical acceptance.
    /// Every state must belong to its exact original chart. Local, cumulative
    /// and total energy defects must fit the existing absolute heat allowance;
    /// no returned allowance can enlarge it and no state/heat is repaired.
    ///
    /// The adapter does not infer topology or a heat-transfer law. Endpoint
    /// contact gating is a time discretization, not event-resolved conduction.
    /// Bounded request/proposal storage allocates for at most 64 reservoirs.
    pub fn step_with_coupled_transport<R>(&mut self, forces:&[f64], gate:&CancelGate,
        mut transport:impl FnMut(ImpactThermalTrial<'_>)
            -> Result<CoupledThermalStep<R>,ImpactError>,
    ) -> Result<TransportedImpactFrame<R>,ImpactError> {
        let Self {mechanics,stores,routes,film_routes,memory_routes,heat,heat_tolerance_j,..}=self;
        let mut proposals=Vec::new();
        let mut inventories=Vec::with_capacity(stores.len());
        let (mut deposited,mut change,mut boundary)=(0.0,0.0,0.0);
        let mechanical=mechanics.step_checked(forces,gate,|system,next,velocity,effort,record| {
            deposited=stage_heat(system,next,velocity,effort,record,routes,film_routes,memory_routes,heat)?;
            if gate.is_requested() {return Err(ImpactError::Cancelled);}
            let inputs:Vec<_>=stores.iter().enumerate().map(|(i,s)|ReservoirThermalInput {
                key:&s.key,thermal:UniformEnthalpyStepInput {
                    curve:s.curve,initial:s.state,mass_kg:s.mass_kg,volume_m3:s.volume_m3,
                    surface_area_m2:s.surface_area_m2,internal_heat_j:heat[i],duration_s:system.config.dt_s,
                },
            }).collect();
            let proposed=transport(ImpactThermalTrial {system,next,inputs:&inputs})?;
            let declared_boundary=proposed.boundary_heat_j;
            proposals=proposed.transfers;
            if proposals.len()!=stores.len() {
                return Err(invalid("coupled thermal transport must return exactly one proposal per reservoir"));
            }
            for ((i,s),p) in stores.iter().enumerate().zip(&proposals) {
                if gate.is_requested() {return Err(ImpactError::Cancelled);}
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
                inventories.push((received,received_correction,external,external_correction));
            }
            if !declared_boundary.is_finite() || (boundary-declared_boundary).abs()>*heat_tolerance_j {
                return Err(invalid("thermal exchange does not close the declared collection boundary"));
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

