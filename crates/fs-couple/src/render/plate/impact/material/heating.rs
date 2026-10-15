//! Contact dissipation -> finite material enthalpy, on one accepted time step.
//!
//! Contact force/storage and stepping remain with fs-dcontact/fs-phs. Thermal
//! charts remain with fs-material. Only the selected contact losses are routed;
//! modal damping, felt crush, radiation and fluid losses are NOT relabelled heat.
//! The mechanical constitutive state is frozen: this is one-way thermalization,
//! not temperature-dependent contact, expansion, melting mechanics or flash T.
/// Existing thermal transport owners applied on the same accepted collision step.
pub mod transport;
/// Explicit unloaded source refresh from thermal state, on the same transaction.
pub mod feedback;
use fs_exec::CancelGate;
use fs_material::phase::{EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState};
use fs_phs::PreparedStepRecord;
use super::super::{ImpactError, ImpactFrame, ImpactSystem, PreparedImpactSystem,
    MAX_IMPACT_MODES, invalid};

/// One finite, uniform thermal destination, with its original material chart.
/// Mass and geometry are explicit specimen data, not inferred from a surface.
#[derive(Clone, Debug, PartialEq)]
pub struct CaloricReservoir<'a> {
    key: String,
    curve: &'a EquilibriumEnthalpyPhaseCurve,
    mass_kg: f64,
    volume_m3: f64,
    surface_area_m2: f64,
    initial: EquilibriumPhaseState,
    state: EquilibriumPhaseState,
    received_heat_j: f64,
    received_correction_j: f64,
    boundary_heat_j: f64,
    boundary_correction_j: f64,
    adiabatic_base: EquilibriumPhaseState,
    adiabatic_heat_j: f64,
}
impl<'a> CaloricReservoir<'a> {
    /// Attach an admitted initial state and invariant mass to a uniform store.
    /// The caller owns lumping, density/volume consistency and negligible pV work.
    pub fn new(key: String, curve: &'a EquilibriumEnthalpyPhaseCurve,
        initial: EquilibriumPhaseState, mass_kg: f64, volume_m3: f64,
        surface_area_m2: f64) -> Result<Self, ImpactError> {
        if !key_valid(&key) || ![mass_kg, volume_m3, surface_area_m2].into_iter().all(positive)
            || curve.state_at_specific_enthalpy(initial.specific_enthalpy_j_kg())
                .map_err(|e|ImpactError::Owner(e.to_string()))? != initial {
            return Err(invalid("contact heat store needs a key, physical specimen and its own chart state"));
        }
        Ok(Self { key,curve,mass_kg,volume_m3,surface_area_m2,initial,state:initial,received_heat_j:0.0,received_correction_j:0.0,
            boundary_heat_j:0.0,boundary_correction_j:0.0,adiabatic_base:initial,adiabatic_heat_j:0.0 })
    }
    /// Stable destination name, independent of contact and mode indexing.
    pub fn key(&self) -> &str { &self.key }
    /// Borrow of the original source-bound equilibrium material chart.
    pub const fn curve(&self) -> &'a EquilibriumEnthalpyPhaseCurve { self.curve }
    /// Conserved thermal mass [kg].
    pub const fn mass_kg(&self) -> f64 { self.mass_kg }
    /// Prescribed specimen volume [m3], not updated by thermal expansion.
    pub const fn volume_m3(&self) -> f64 { self.volume_m3 }
    /// Prescribed exposed surface area [m2].
    pub const fn surface_area_m2(&self) -> f64 { self.surface_area_m2 }
    /// Accepted temperature, liquid mass fraction, density and specific enthalpy.
    pub const fn state(&self) -> EquilibriumPhaseState { self.state }
    /// Cumulative deposited contact heat since attachment [J].
    pub const fn received_heat_j(&self) -> f64 { self.received_heat_j }
    /// Cumulative signed boundary transfer into this store [J].
    pub const fn boundary_heat_j(&self) -> f64 { self.boundary_heat_j }
    /// Actual chart storage change since attachment [J]; invariant to h datum.
    pub fn storage_change_j(&self) -> f64 {
        self.mass_kg*(self.state.specific_enthalpy_j_kg()-self.initial.specific_enthalpy_j_kg())
    }
}

/// Explicit partition of ONE contact's irreversible work between two stores.
/// No implicit equal split or material-name-based partition is supplied. A
/// store may receive several contacts; both names may designate the same store.
#[derive(Clone, Debug, PartialEq)]
pub struct ContactHeatRoute {
    /// Existing immutable material contact key, not a mutable array position.
    pub contact: String,
    /// First destination's stable key.
    pub to_a: String,
    /// Second destination's stable key.
    pub to_b: String,
    /// Fraction of contact work delivered to A, in [0,1]; B gets the remainder.
    pub fraction_to_a: f64,
}
#[derive(Clone, Copy)]
struct Route { contact: usize, a: usize, b: usize, fraction: f64 }

/// Combined accepted mechanics and thermal accounting. Mechanical dissipation
/// is still reported by its original owner; deposited heat must NOT be counted
/// as an additional loss from the combined system.
#[derive(Clone, Debug)]
pub struct HeatedImpactFrame {
    /// Unmodified accepted mechanical result, including every loss channel.
    pub mechanical: ImpactFrame,
    /// Selected contact loss deposited into thermal destinations this step [J].
    pub deposited_heat_j: f64,
    /// Actual sum of mass times chart enthalpy changes this step [J].
    pub thermal_energy_change_j: f64,
    /// Mechanical losses not routed into these stores [J]. Signed roundoff in
    /// this subtraction is disclosed, not silently labelled physical cooling.
    pub remaining_mechanical_loss_j: f64,
    /// Mechanical balance plus thermal change minus contact and boundary heat [J].
    pub combined_energy_residual_j: f64,
}

/// Existing prepared mechanics with an atomic, source-bound thermal destination.
/// Successful prepared stepping remains the original equation. This wrapper
/// allocates no hot-path buffers; chart evaluation is outside the Newton solve.
/// No hard-real-time guarantee, spatial heat solve or thermal feedback is claimed.
/// Internal adaptive substepping is not exposed by this fixed-step wrapper.
pub struct HeatedImpactSystem<'a> {
    mechanics: PreparedImpactSystem,
    stores: Vec<CaloricReservoir<'a>>,
    routes: Vec<Route>,
    candidate: Vec<EquilibriumPhaseState>,
    next_received: Vec<f64>,
    next_correction: Vec<f64>,
    heat: Vec<f64>,
    heat_tolerance_j: f64,
}
impl<'a> HeatedImpactSystem<'a> {
    /// Attach explicit stores and routes at the mechanics' CURRENT accepted time.
    /// Already dissipated work is not replayed. Unknown/duplicate contacts or
    /// stores and invalid partitions refuse. Unrouted losses remain external.
    pub fn new(mechanics: PreparedImpactSystem, stores: Vec<CaloricReservoir<'a>>,
        routes: Vec<ContactHeatRoute>, heat_tolerance_j: f64) -> Result<Self, ImpactError> {
        if stores.is_empty() || stores.len()>64 || routes.is_empty() || routes.len()>32
            || !positive(heat_tolerance_j) {
            return Err(invalid("contact heating needs 1..64 stores, 1..32 routes and a positive joule allowance"));
        }
        for (i,s) in stores.iter().enumerate() {
            if stores[..i].iter().any(|old|old.key==s.key) {
                return Err(invalid("duplicate contact heat store"));
            }
        }
        let mut compiled: Vec<Route> = Vec::with_capacity(routes.len());
        for route in routes {
            if !route.fraction_to_a.is_finite() || !(0.0..=1.0).contains(&route.fraction_to_a) {
                return Err(invalid("heat partition must be finite and in [0,1]"));
            }
            let contact=mechanics.material_contacts().iter().position(|c|c.key()==route.contact)
                .ok_or_else(||invalid("heat route refers to an unknown source contact"))?;
            if compiled.iter().any(|r|r.contact==contact) {return Err(invalid("contact heat must be routed at most once"));}
            let a=stores.iter().position(|s|s.key==route.to_a).ok_or_else(||invalid("unknown first heat store"))?;
            let b=stores.iter().position(|s|s.key==route.to_b).ok_or_else(||invalid("unknown second heat store"))?;
            compiled.push(Route {contact,a,b,fraction:route.fraction_to_a});
        }
        // Canonical contact order makes input route ordering immaterial.
        compiled.sort_by_key(|r|r.contact);
        let candidate=stores.iter().map(|s|s.state).collect();
        let next_received=stores.iter().map(|s|s.received_heat_j).collect();
        let next_correction=stores.iter().map(|s|s.received_correction_j).collect();
        let heat=vec![0.0;stores.len()];
        Ok(Self {mechanics,stores,routes:compiled,candidate,next_received,next_correction,heat,heat_tolerance_j})
    }
    /// Read-only mechanics: callers cannot advance it behind the thermal clock.
    pub const fn mechanics(&self) -> &PreparedImpactSystem { &self.mechanics }
    /// Accepted stores in their original input order, unchanged on refusal.
    pub fn reservoirs(&self) -> &[CaloricReservoir<'a>] { &self.stores }
    /// Select the existing analytic/finite-difference Newton image, without reset.
    pub fn set_analytic_newton(&mut self, enabled: bool) {self.mechanics.set_analytic_newton(enabled);}
    /// Set the existing numerical work ceiling; no thermal/model change.
    pub fn set_iteration_limit(&mut self, limit: usize) -> Result<(), ImpactError> {
        self.mechanics.set_iteration_limit(limit)
    }
    /// Release both accepted physical owners together without resetting either.
    pub fn into_parts(self) -> (PreparedImpactSystem, Vec<CaloricReservoir<'a>>) {
        (self.mechanics,self.stores)
    }
    /// Advance motion and thermal storage exactly once, or neither. Heat uses
    /// the original contact loss at the solver midpoint and its ACTUAL port
    /// velocity (Gonzalez effort), not endpoint power times dt. Thermal limits
    /// and rounding are checked before mechanical history/state publication.
    pub fn step(&mut self, forces: &[f64], gate: &CancelGate) -> Result<HeatedImpactFrame, ImpactError> {
        let Self {mechanics,stores,routes,candidate,next_received,next_correction,heat,heat_tolerance_j}=self;
        let mut deposited=0.0; let mut thermal_change=0.0;
        let mechanical=mechanics.step_checked(forces,gate,|system,next,velocity,record| {
            deposited=stage_heat(system,next,velocity,record,routes,heat)?;
            thermal_change=0.0;
            for (i,store) in stores.iter().enumerate() {
                if gate.is_requested() {return Err(ImpactError::Cancelled);}
                if heat[i]==0.0 {
                    candidate[i]=store.state;next_received[i]=store.received_heat_j;
                    next_correction[i]=store.received_correction_j;continue;
                }
                let increment=heat[i]-store.received_correction_j;
                let total=store.received_heat_j+increment;
                next_correction[i]=(total-store.received_heat_j)-increment;
                let h=store.adiabatic_base.specific_enthalpy_j_kg()
                    +(total-store.adiabatic_heat_j)/store.mass_kg;
                if !total.is_finite() || !h.is_finite() {return Err(invalid("contact thermal inventory overflow"));}
                candidate[i]=store.curve.state_at_specific_enthalpy(h)
                    .map_err(|e|ImpactError::Owner(e.to_string()))?;
                next_received[i]=total;
                let actual_h=candidate[i].specific_enthalpy_j_kg();
                thermal_change+=store.mass_kg*(actual_h-store.state.specific_enthalpy_j_kg());
                // Check cumulative loss of small increments as well as this step.
                let inventory=store.mass_kg*(actual_h-store.initial.specific_enthalpy_j_kg());
                if !inventory.is_finite() || (inventory-total-store.boundary_heat_j).abs()>*heat_tolerance_j {
                    return Err(invalid("contact thermal inventory exceeds its joule allowance"));
                }
            }
            if !thermal_change.is_finite() || (thermal_change-deposited).abs()>*heat_tolerance_j {
                return Err(invalid("contact heat transfer exceeds its joule allowance"));
            }
            Ok(())
        })?;
        // All fallible work (including final cancellation and mechanical gates)
        // precedes publication. No callbacks or allocations occur in this block.
        for (i,store) in stores.iter_mut().enumerate() {
            store.state=candidate[i];store.received_heat_j=next_received[i];
            store.received_correction_j=next_correction[i];
        }
        Ok(HeatedImpactFrame {deposited_heat_j:deposited,thermal_energy_change_j:thermal_change,
            remaining_mechanical_loss_j:mechanical.dissipated_energy_j-deposited,
            combined_energy_residual_j:mechanical.balance_residual_j+thermal_change-deposited,
            mechanical})
    }
}
fn positive(v:f64)->bool {v.is_finite() && v>0.0}
fn key_valid(key:&str)->bool {!key.trim().is_empty() && key.len()<=128}

fn stage_heat(system:&ImpactSystem, next:&[f64], velocity:&[f64], record:PreparedStepRecord,
    routes:&[Route], heat:&mut[f64]) -> Result<f64,ImpactError> {
    let n=system.modes;
    let mut midpoint=[0.0;2*MAX_IMPACT_MODES];let mut power=[0.0;32];
    for i in 0..2*n {midpoint[i]=0.5*(system.x[i]+next[i]);}
    let count=system.material_contacts().len();
    if !system.contact.dissipated_power_by_obstacle_into(&midpoint[..2*n],velocity,&mut power[..count]) {
        return Err(invalid("original contact owner refused step heat evaluation"));
    }
    heat.fill(0.0);let mut deposited=0.0;
    for r in routes {
        let q=system.config.dt_s*power[r.contact];
        let a=q*r.fraction; let b=q-a;
        if !q.is_finite() || q<0.0 || (power[r.contact]>0.0 && q==0.0) {
            return Err(invalid("contact heat is unrepresentable"));
        }
        heat[r.a]+=a;heat[r.b]+=b;deposited+=q;
    }
    let roundoff=256.0*f64::EPSILON*(record.dissipated.abs()+deposited.abs());
    if heat.iter().any(|q|!q.is_finite() || *q<0.0) || !deposited.is_finite()
        || deposited>record.dissipated+roundoff {
        return Err(invalid("routed contact heat exceeds the time owner's dissipation"));
    }
    Ok(deposited)
}
