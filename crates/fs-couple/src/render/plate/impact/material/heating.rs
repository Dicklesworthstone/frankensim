//! Selected physical dissipation -> finite material enthalpy on one time step.
//!
//! Contact force/storage and stepping remain with fs-dcontact/fs-phs. Thermal
//! charts remain with fs-material. Selected contact, incompressible film and
//! Maxwell-arm losses use their original owners at the actual discrete-gradient effort.
//! Modal damping, felt crush, radiation and compressible-film free-energy loss
//! are not relabelled heat. No dissipation is charged to mechanics twice.
//! The mechanical constitutive state is frozen: this is one-way thermalization,
//! not temperature-dependent contact, expansion, melting mechanics or flash T.
/// Existing thermal transport owners applied on the same accepted collision step.
pub mod transport;
/// Explicit unloaded source refresh from thermal state, on the same transaction.
pub mod feedback;
/// Joint contact heating through the existing material enthalpy network.
pub mod network;
/// Live sourced liquid viscosity, joint film work and the original heat owners.
pub mod liquid_film;
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
    /// Cumulative selected mechanical heat since attachment [J].
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

/// Explicit partition of one incompressible pressure island's viscous work.
/// The index is the existing immutable film insertion order. Destinations may
/// be the lubricant's own caloric store or explicitly modeled thermal sinks;
/// no wall/fluid heat partition is inferred from geometry or material names.
#[derive(Clone, Debug, PartialEq)]
pub struct SqueezeFilmHeatRoute {
    /// Existing pressure island, in `ImpactSystem::squeeze_film_count` order.
    pub film: usize,
    /// First finite destination's stable key.
    pub to_a: String,
    /// Second finite destination's stable key; may equal the first.
    pub to_b: String,
    /// Fraction delivered to A, in [0,1]; B receives the remainder.
    pub fraction_to_a: f64,
}
#[derive(Clone, Copy)]
struct FilmRoute { film: usize, a: usize, b: usize, fraction: f64 }

/// Destination of one original Maxwell arm's irreversible material work.
/// Branch indices follow `ImpactSystem::with_relaxation_branches` insertion
/// order, including arms compiled from a material shell spectrum. A branch's
/// heat belongs to one explicitly modeled caloric specimen; neither a material
/// name nor a mechanical mode determines that destination. Several arms may
/// heat the same store, and unselected arms remain outside the thermal boundary.
#[derive(Clone, Debug, PartialEq)]
pub struct RelaxationHeatRoute {
    /// Existing scalar material-memory branch, not a mechanical mode index.
    pub branch: usize,
    /// Stable key of the finite thermal destination.
    pub to: String,
}
#[derive(Clone, Copy)]
struct MemoryRoute { branch: usize, to: usize }

/// Combined accepted mechanics and thermal accounting. Mechanical dissipation
/// is still reported by its original owner; deposited heat must NOT be counted
/// as an additional loss from the combined system.
#[derive(Clone, Debug)]
pub struct HeatedImpactFrame {
    /// Unmodified accepted mechanical result, including every loss channel.
    pub mechanical: ImpactFrame,
    /// Selected contact, viscous-film and material-memory loss this step [J].
    pub deposited_heat_j: f64,
    /// Actual sum of mass times chart enthalpy changes this step [J].
    pub thermal_energy_change_j: f64,
    /// Mechanical losses not routed into these stores [J]. Signed roundoff in
    /// this subtraction is disclosed, not silently labelled physical cooling.
    pub remaining_mechanical_loss_j: f64,
    /// Mechanical balance plus thermal change minus deposited and boundary heat [J].
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
    film_routes: Vec<FilmRoute>,
    memory_routes: Vec<MemoryRoute>,
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
        Self::new_with_films(mechanics, stores, routes, Vec::new(), heat_tolerance_j)
    }

    /// Attach contact and Reynolds-film heat to the SAME material stores and
    /// thermal network. A film-only system is admitted. Each selected channel
    /// is routed once, at the pressure solve's midpoint configuration and actual
    /// work-conjugate velocity, then deposited before mechanical publication.
    ///
    /// The incompressible film has no stored energy, so its nonnegative port
    /// work is physical viscous heat. Compressible-film free-energy dissipation
    /// is deliberately excluded: it is not automatically finite internal heat.
    /// Fluid viscosity is frozen over the sample; temperature feedback requires
    /// an explicit constitutive refresh, not a second thermal or pressure step.
    pub fn new_with_films(mechanics: PreparedImpactSystem,
        stores: Vec<CaloricReservoir<'a>>, routes: Vec<ContactHeatRoute>,
        film_routes: Vec<SqueezeFilmHeatRoute>, heat_tolerance_j: f64,
    ) -> Result<Self, ImpactError> {
        Self::new_with_losses(mechanics, stores, routes, film_routes, Vec::new(), heat_tolerance_j)
    }

    /// Route selected physical loss channels through the same accepted step.
    ///
    /// Material-memory heat uses `dt * (1/tau) * discrete_effort_z^2` from the
    /// original Maxwell resistance and the successful time owner's COMPLETE
    /// discrete gradient. It includes the nonlinear Gonzalez correction; it is
    /// neither endpoint power nor a mechanical energy-deficit estimate. The
    /// arm's stored energy stays in mechanics and is never deposited as heat.
    ///
    /// A memory-only system is supported. Each arm is routed at most once.
    /// This thermalizes an admitted constitutive law without changing its
    /// coefficients: temperature-dependent modulus/relaxation feedback and
    /// spatial heat allocation require their own explicit physical binding.
    pub fn new_with_losses(mechanics: PreparedImpactSystem,
        stores: Vec<CaloricReservoir<'a>>, routes: Vec<ContactHeatRoute>,
        film_routes: Vec<SqueezeFilmHeatRoute>, memory_routes: Vec<RelaxationHeatRoute>,
        heat_tolerance_j: f64,
    ) -> Result<Self, ImpactError> {
        if stores.is_empty() || stores.len()>64
            || (routes.is_empty() && film_routes.is_empty() && memory_routes.is_empty())
            || routes.len()>32 || film_routes.len()>4 || memory_routes.len()>256
            || !positive(heat_tolerance_j) {
            return Err(invalid("mechanical heating needs 1..64 stores, up to 32 contacts, four films and 256 material arms, at least one route and a positive joule allowance"));
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
        let mut compiled_films: Vec<FilmRoute> = Vec::with_capacity(film_routes.len());
        for route in film_routes {
            if route.film >= mechanics.squeeze_film_count()
                || compiled_films.iter().any(|r| r.film == route.film) {
                return Err(invalid("film heat must name an existing pressure island exactly once"));
            }
            if !route.fraction_to_a.is_finite() || !(0.0..=1.0).contains(&route.fraction_to_a) {
                return Err(invalid("film heat partition must be finite and in [0,1]"));
            }
            let a=stores.iter().position(|s|s.key==route.to_a).ok_or_else(||invalid("unknown first film heat store"))?;
            let b=stores.iter().position(|s|s.key==route.to_b).ok_or_else(||invalid("unknown second film heat store"))?;
            compiled_films.push(FilmRoute {film:route.film,a,b,fraction:route.fraction_to_a});
        }
        compiled_films.sort_by_key(|r|r.film);
        let mut compiled_memory: Vec<MemoryRoute> = Vec::with_capacity(memory_routes.len());
        for route in memory_routes {
            if route.branch >= mechanics.relaxation_branch_count()
                || compiled_memory.iter().any(|r|r.branch == route.branch) {
                return Err(invalid("material heat must name an existing relaxation arm exactly once"));
            }
            let to=stores.iter().position(|s|s.key==route.to)
                .ok_or_else(||invalid("unknown material-memory heat store"))?;
            compiled_memory.push(MemoryRoute {branch:route.branch,to});
        }
        compiled_memory.sort_by_key(|r|r.branch);
        let candidate=stores.iter().map(|s|s.state).collect();
        let next_received=stores.iter().map(|s|s.received_heat_j).collect();
        let next_correction=stores.iter().map(|s|s.received_correction_j).collect();
        let heat=vec![0.0;stores.len()];
        Ok(Self {mechanics,stores,routes:compiled,film_routes:compiled_films,memory_routes:compiled_memory,
            candidate,next_received,next_correction,heat,heat_tolerance_j})
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
    /// the original selected loss at the solver midpoint and its ACTUAL
    /// discrete effort, not endpoint power times dt. Thermal limits
    /// and rounding are checked before mechanical history/state publication.
    pub fn step(&mut self, forces: &[f64], gate: &CancelGate) -> Result<HeatedImpactFrame, ImpactError> {
        let Self {mechanics,stores,routes,film_routes,memory_routes,candidate,next_received,next_correction,heat,heat_tolerance_j}=self;
        let mut deposited=0.0; let mut thermal_change=0.0;
        let mechanical=mechanics.step_checked(forces,gate,|system,next,velocity,effort,record| {
            deposited=stage_heat(system,next,velocity,effort,record,routes,film_routes,memory_routes,heat)?;
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

fn stage_heat(system:&ImpactSystem, next:&[f64], velocity:&[f64], effort:&[f64], record:PreparedStepRecord,
    routes:&[Route], film_routes:&[FilmRoute], memory_routes:&[MemoryRoute], heat:&mut[f64]) -> Result<f64,ImpactError> {
    let n=system.modes;
    let mut midpoint=[0.0;2*MAX_IMPACT_MODES];let mut power=[0.0;32];
    for i in 0..2*n {midpoint[i]=0.5*(system.x[i]+next[i]);}
    let count=system.material_contacts().len();
    if !routes.is_empty()
        && !system.contact.dissipated_power_by_obstacle_into(&midpoint[..2*n],velocity,&mut power[..count]) {
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
    if !film_routes.is_empty() {
        let mut configuration=[0.0;MAX_IMPACT_MODES];
        let mut force=[0.0;MAX_IMPACT_MODES];
        let mut pressure=[0.0;fs_tribo::resistive_film::MAX_CELLS];
        for i in 0..n {configuration[i]=midpoint[2*i];}
        for r in film_routes {
            let film=&system.squeeze_films[r.film];
            let report=film.evaluate_into(&configuration[..n],velocity,&mut force[..n],
                &mut pressure[..film.cell_count()]).map_err(|e|invalid(e.0))?;
            let q=system.config.dt_s*report.dissipated_power_w;
            if !q.is_finite() || q<0.0 || (report.dissipated_power_w>0.0 && q==0.0) {
                return Err(invalid("viscous film heat is unrepresentable"));
            }
            let a=q*r.fraction;
            heat[r.a]+=a;heat[r.b]+=q-a;deposited+=q;
        }
    }
    if !memory_routes.is_empty() {
        let memory=system.relaxation.as_ref()
            .ok_or_else(||invalid("material heat requires the original relaxation owner"))?;
        for r in memory_routes {
            let power=memory.branch_dissipated_power(effort,r.branch)?;
            let q=system.config.dt_s*power;
            if !q.is_finite() || q<0.0 || (power>0.0 && q==0.0) {
                return Err(invalid("material relaxation heat is unrepresentable"));
            }
            heat[r.to]+=q;deposited+=q;
        }
    }
    let roundoff=256.0*f64::EPSILON*(record.dissipated.abs()+deposited.abs());
    if heat.iter().any(|q|!q.is_finite() || *q<0.0) || !deposited.is_finite()
        || deposited>record.dissipated+roundoff {
        return Err(invalid("routed physical heat exceeds the time owner's dissipation"));
    }
    Ok(deposited)
}
