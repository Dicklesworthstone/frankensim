//! Reciprocal elastic/finite-enthalpy feedback on the original P1 load owner.
//!
//! Each trial advances the ORIGINAL heat timeline from accepted history, then
//! solves native compatible equilibrium at its actual local temperatures.
//! Endpoint equilibrium gives W=sum V*(sigma0+sigma1):delta(epsilon)/2
//! = (f0+f1).delta(u)/2 for unchanged homogeneous supports. This is a declared
//! trapezoidal load/work quadrature, not exact integration of an arbitrary
//! changing-load/temperature path. Local reversible heat is W_e-delta(U_e).
//! Source-derived entropy distinguishes this exchange from simply converting
//! elastic coefficient work to heat. No solver residual is deposited as heat.
//!
//! The bounded initial law is isotropic affine E(T), fixed nu and constant
//! instantaneous alpha, on a fully solid constant-reference-density interval
//! of each original caloric chart. Later phase transitions remain in that chart.
//! Spatial temperatures, native stress/deformation, source pins, reference mass,
//! and the original thermal observer are retained atomically. The endpoint
//! coupling is generally first order; it does not establish exact reversible
//! entropy conservation or add inertia, plasticity, phase flow, moving supports,
//! follower loads or deformation-dependent heat-transfer geometry.

/// Reciprocal finite heat coupled to the original inertial Newmark trajectory.
pub mod dynamic;

/// Local isentropic wave/modal tangent from the same accepted thermal history.
pub mod acoustic;

/// Complete original material heat transport on the same reciprocal trajectory.
pub mod heat;

use super::*;
use std::convert::Infallible;
use fs_material::elastic::thermal::{
    AffineIsotropicThermoelastic, ThermoelasticError, ThermoelasticState,
};
use fs_thermal::enthalpy::HeatTransfer;

pub const MAX_THERMOELASTIC_ITERATIONS:usize=256;

#[derive(Clone,Copy,Debug)]
pub struct ThermoelasticSolve {
    pub maximum_iterations:usize,
    /// Fixed-point relaxation of the per-element reversible joules, in (0,1].
    pub heat_relaxation:f64,
    /// Sum of absolute per-element applied/source heat mismatches [J].
    pub heat_tolerance_j:f64,
    /// Independent local/global mechanical and total energy gates [J].
    pub energy_tolerance_j:f64,
    /// Native stress versus original material potential admission [Pa].
    pub stress_tolerance_pa:f64,
}
#[derive(Debug)]
pub enum ThermoelasticFeedbackError<E=Infallible> {
    Invalid(&'static str),
    Source(ThermoelasticError),
    Solid(SolidEnthalpyError),
    Thermal(DepositionError),
    Observer(E),
    NotConverged {iterations:usize,heat_defect_j:f64,combined_residual_j:f64},
    Cancelled,
}
impl<E> From<DepositionError> for ThermoelasticFeedbackError<E> {
    fn from(error:DepositionError)->Self {
        if error==DepositionError::Cancelled {Self::Cancelled} else {Self::Thermal(error)}
    }
}
impl<E> From<SolidEnthalpyError> for ThermoelasticFeedbackError<E> {
    fn from(error:SolidEnthalpyError)->Self {
        match error {
            SolidEnthalpyError::Cancelled|SolidEnthalpyError::Thermal(DepositionError::Cancelled)=>
                Self::Cancelled,
            other=>Self::Solid(other),
        }
    }
}
impl<E> From<ThermoelasticError> for ThermoelasticFeedbackError<E> {
    fn from(error:ThermoelasticError)->Self {Self::Source(error)}
}
impl<E:core::fmt::Debug> core::fmt::Display for ThermoelasticFeedbackError<E> {
    fn fmt(&self,f:&mut core::fmt::Formatter<'_>)->core::fmt::Result {
        write!(f,"reciprocal spatial thermoelasticity: {self:?}")
    }
}
impl<E:std::error::Error+'static> std::error::Error for ThermoelasticFeedbackError<E> {}
fn check<E>(cx:&Cx<'_>)->Result<(),ThermoelasticFeedbackError<E>> {
    cx.checkpoint().map_err(|_|ThermoelasticFeedbackError::Cancelled)
}
fn finite<E>(value:f64)->Result<f64,ThermoelasticFeedbackError<E>> {
    if value.is_finite() {Ok(value)} else {
        Err(ThermoelasticFeedbackError::Invalid("thermoelastic work or heat is unrepresentable"))
    }
}
fn controls<E>(config:ThermoelasticSolve)->Result<(),ThermoelasticFeedbackError<E>> {
    if config.maximum_iterations==0 || config.maximum_iterations>MAX_THERMOELASTIC_ITERATIONS
        || !config.heat_relaxation.is_finite() || config.heat_relaxation<=0.0 || config.heat_relaxation>1.0
        || [config.heat_tolerance_j,config.energy_tolerance_j,config.stress_tolerance_pa]
            .iter().any(|v|!v.is_finite()||*v<=0.0) {
        return Err(ThermoelasticFeedbackError::Invalid("bounded iterations, relaxation and positive finite tolerances required"));
    }
    Ok(())
}

/// Retained native equilibrium and its original finite thermal/material history.
pub struct ThermoelasticSolid<'s,'m> {
    source:&'s SolidEnthalpy<'m>,
    laws:Vec<AffineIsotropicThermoelastic>,
    state:SolidEnthalpyState<'m>,
    equilibrium:SolidEquilibrium,
    potentials:Vec<ThermoelasticState>,
    fixed_dofs:Vec<usize>,
    nodal_forces_n:Vec<(usize,f64)>,
    solve:SolidEquilibriumConfig,
}
impl<'s,'m> ThermoelasticSolid<'s,'m> {
    /// Start from an actual native load-bearing equilibrium. One source interval
    /// per tetrahedron must cover the original stress-free reference and the
    /// current state, and fit inside a fully solid interval of the SAME caloric
    /// chart. A later phase transition elsewhere on that chart is retained.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(source:&'s SolidEnthalpy<'m>,thermal:&DepositionSimulation<'_>,
        loads:&TetLoadCase3<'_,'_>,fixed_dofs:&[usize],temperature_ranges_k:&[[f64;2]],
        solve:SolidEquilibriumConfig,config:ThermoelasticSolve,cx:&Cx<'_>)
        ->Result<Self,ThermoelasticFeedbackError> {
        check(cx)?; controls(config)?;
        if temperature_ranges_k.len()!=source.cells.len() {
            return Err(ThermoelasticFeedbackError::Invalid("one thermal source interval per original tetrahedron required"));
        }
        let state=source.resolve(thermal,cx)?;
        let mut laws=Vec::with_capacity(source.cells.len());
        for (index,(bound,range)) in source.cells.iter().zip(temperature_ranges_k).enumerate() {
            check(cx)?;
            if !matches!(bound.source.stiffness_model,ThermalStiffnessModel::Isotropic)
                || !matches!(bound.source.expansion_model,ThermalExpansionModel::Isotropic) {
                return Err(ThermoelasticFeedbackError::Invalid("reciprocal potential requires the declared isotropic stiffness and expansion models"));
            }
            let sample=&state.samples[index];
            let law=AffineIsotropicThermoelastic::new(bound.source.point.material,&sample.elastic,
                sample.expansion.reference(),*range,source.maximum_absolute_mandel_strain)?;
            let knots=bound.source.curve.knots();
            if range[0]<knots[0].temperature_k || range[1]>knots[knots.len()-1].temperature_k {
                return Err(ThermoelasticFeedbackError::Invalid("caloric chart must cover the whole thermoelastic source interval"));
            }
            let density=bound.reference_cell.mass_kg/bound.reference_volume_m3;
            let mut minimum=f64::INFINITY;
            for pair in knots.windows(2) {
                check(cx)?;
                let dt=pair[1].temperature_k-pair[0].temperature_k;
                if dt<=0.0 {
                    if pair[0].temperature_k>=range[0] && pair[0].temperature_k<=range[1] {
                        return Err(ThermoelasticFeedbackError::Invalid("the admitted thermoelastic interval cannot intersect an isothermal phase segment"));
                    }
                    continue;
                }
                if pair[0].temperature_k>=range[1] || pair[1].temperature_k<=range[0] {continue;}
                if pair.iter().any(|knot|knot.liquid_mass_fraction!=0.0
                    || !same_measure(knot.bulk_density_kg_m3,density)) {
                    return Err(ThermoelasticFeedbackError::Invalid("the admitted caloric interval must be fully solid at constant reference density"));
                }
                minimum=minimum.min(finite(bound.reference_cell.mass_kg
                    *(pair[1].specific_enthalpy_j_kg-pair[0].specific_enthalpy_j_kg)/dt)?);
            }
            if !minimum.is_finite()
                || finite(minimum+law.minimum_heat_capacity_j_k(bound.reference_volume_m3)?)?<=0.0 {
                return Err(ThermoelasticFeedbackError::Invalid("caloric plus elastic-field capacity must be positive over the full admitted strain/temperature envelope"));
            }
            laws.push(law);
        }
        let equilibrium=state.solve(loads,fixed_dofs,solve,cx)?;
        let potentials=potentials(&state,&equilibrium,&laws,config,cx)?;
        check(cx)?;
        Ok(Self {source,laws,state,equilibrium,potentials,fixed_dofs:fixed_dofs.to_vec(),
            nodal_forces_n:loads.nodal_forces_n().to_vec(),solve})
    }
    pub fn state(&self)->&SolidEnthalpyState<'m> {&self.state}
    pub fn equilibrium(&self)->&SolidEquilibrium {&self.equilibrium}
    pub fn potentials(&self)->&[ThermoelasticState] {&self.potentials}
    pub fn time_s(&self)->f64 {self.state.time_s()}

    fn admit_thermal<E>(&self,thermal:&DepositionSimulation<'_>,accepted:bool,cx:&Cx<'_>)
        ->Result<SolidEnthalpyState<'m>,ThermoelasticFeedbackError<E>> {
        check(cx)?;
        let state=self.source.resolve(thermal,cx)?;
        if accepted && (state.time_s.to_bits()!=self.state.time_s.to_bits()
            || state.materials!=self.state.materials || state.strains!=self.state.strains
            || state.samples.iter().zip(&self.state.samples).any(|(a,b)|a.phase!=b.phase)) {
            return Err(ThermoelasticFeedbackError::Invalid("accepted mechanical and caloric histories must remain exactly paired"));
        }
        for (law,sample) in self.laws.iter().zip(&state.samples) {
            check(cx)?;
            law.at_materials(&sample.elastic,sample.expansion.current())?;
        }
        Ok(state)
    }

    fn mechanics<E>(&self,state:SolidEnthalpyState<'m>,loads:&TetLoadCase3<'_,'_>,
        config:ThermoelasticSolve,cx:&Cx<'_>)->Result<MechanicalProposal<'m>,ThermoelasticFeedbackError<E>> {
        let equilibrium=state.solve(loads,&self.fixed_dofs,self.solve,cx)?;
        let potentials=potentials(&state,&equilibrium,&self.laws,config,cx)?;
        let mut cells=Vec::with_capacity(self.source.cells.len());
        let mut work=0.0; let mut internal_change=0.0; let mut heat=0.0;
        for (i,bound) in self.source.cells.iter().enumerate() {
            check(cx)?;
            let old=&self.equilibrium.solution.stress.elements()[i];
            let next=&equilibrium.solution.stress.elements()[i];
            let a=&self.potentials[i]; let b=&potentials[i];
            let free0=self.state.strains[i].free_strain_mandel();
            let free1=state.strains[i].free_strain_mandel();
            let sigma0=old.stress_mandel_pa(); let sigma1=next.stress_mandel_pa();
            let e0=old.elastic_strain_mandel(); let e1=next.elastic_strain_mandel();
            let volume=bound.reference_volume_m3;
            let mechanical_work_j=finite(0.5*volume*(0..6).map(|k|
                (sigma0[k]+sigma1[k])*(e1[k]+free1[k]-e0[k]-free0[k])).sum::<f64>())?;
            let expansion_work_j=finite(-0.5*volume*(0..6).map(|k|
                (sigma0[k]+sigma1[k])*(free1[k]-free0[k])).sum::<f64>())?;
            let c0=self.state.materials[i].stiffness_mandel_pa();
            let c1=state.materials[i].stiffness_mandel_pa();
            let stiffness_work_j=finite(0.5*volume*(0..6).map(|row|
                e0[row]*(0..6).map(|col|(c1[row][col]-c0[row][col])*e1[col]).sum::<f64>()).sum::<f64>())?;
            let free_energy_residual_j=finite(b.free_energy_j-a.free_energy_j
                -mechanical_work_j-expansion_work_j-stiffness_work_j)?;
            let internal_energy_change_j=finite(b.internal_energy_j-a.internal_energy_j)?;
            let reversible_heat_j=finite(mechanical_work_j-internal_energy_change_j)?;
            let temperature_k=[self.state.samples[i].phase.temperature_k(),state.samples[i].phase.temperature_k()];
            let potential_residual_j=finite(reversible_heat_j+expansion_work_j+stiffness_work_j
                +temperature_k[1]*b.entropy_j_k-temperature_k[0]*a.entropy_j_k)?;
            if free_energy_residual_j.abs()>config.energy_tolerance_j
                || potential_residual_j.abs()>config.energy_tolerance_j {
                return Err(ThermoelasticFeedbackError::Invalid("native stress/material work does not close the source thermodynamic potential"));
            }
            work=finite(work+mechanical_work_j)?;
            heat=finite(heat+reversible_heat_j)?;
            internal_change=finite(internal_change+internal_energy_change_j)?;
            cells.push(ThermoelasticCellStep {cell:bound.source.cell,temperature_k,
                potentials:[*a,*b],mechanical_work_j,expansion_work_j,stiffness_work_j,
                reversible_heat_j,free_energy_residual_j});
        }
        let before=self.equilibrium.solution.displacement.displacement_m();
        let after=equilibrium.solution.displacement.displacement_m();
        let mut external_work_j=0.0;
        for &(dof,force) in self.nodal_forces_n.iter().chain(loads.nodal_forces_n()) {
            check(cx)?;
            external_work_j=finite(external_work_j+0.5*force*(after[dof/3][dof%3]-before[dof/3][dof%3]))?;
        }
        let mechanical_residual_j=finite(work-external_work_j)?;
        if mechanical_residual_j.abs()>config.energy_tolerance_j {
            return Err(ThermoelasticFeedbackError::Invalid("actual averaged nodal work disagrees with native element stress work"));
        }
        Ok(MechanicalProposal {state,equilibrium,potentials,cells,external_work_j,
            internal_energy_change_j:internal_change,reversible_heat_j:heat,mechanical_residual_j})
    }

    /// Advance native equilibrium and the same spatial thermal/observer history
    /// together. All iterations restart from the accepted heat state; only the
    /// balanced final pair and its original observer are returned. Supports and
    /// reference geometry stay fixed. Loads may change; their actual endpoint
    /// nodal forces define the declared trapezoidal work convention.
    #[allow(clippy::too_many_arguments)]
    pub fn advance<'t,O:DepositionObserver>(&self,thermal:&DepositionSimulation<'t>,
        observer:&O,target_s:f64,loads:&TetLoadCase3<'_,'_>,thermal_config:StepConfig,
        config:ThermoelasticSolve,cx:&Cx<'_>)
        ->Result<ThermoelasticCandidate<'s,'m,'t,O>,ThermoelasticFeedbackError<O::Error>> {
        check(cx)?; controls(config)?;
        let dt=target_s-self.time_s();
        if !dt.is_finite() || dt<=0.0 || !target_s.is_finite()
            || !thermal_config.energy_tolerance_j.is_finite()
            || thermal_config.energy_tolerance_j<0.0
            || thermal_config.energy_tolerance_j>config.energy_tolerance_j {
            return Err(ThermoelasticFeedbackError::Invalid("one positive physical interval and compatible thermal energy tolerance required"));
        }
        self.admit_thermal(thermal,true,cx)?;
        let mut applied=vec![0.0;self.source.cells.len()];
        let mut defects=[f64::INFINITY;2];
        for iteration in 1..=config.maximum_iterations {
            check(cx)?;
            let mut watts=BTreeMap::new();
            for (bound,&heat) in self.source.cells.iter().zip(&applied) {
                let power=finite(heat/dt)?;
                if heat!=0.0 && power==0.0 {
                    return Err(ThermoelasticFeedbackError::Invalid("reversible heat underflows its power conversion"));
                }
                if power!=0.0 {watts.insert(bound.source.cell,power);}
            }
            let coupled=ReversibleHeat {solid:self,inner:observer.clone(),cx,watts};
            let (thermal_next,coupled,audit)=thermal.advance_observed(&coupled,target_s,thermal_config,
                ||cx.checkpoint().map_err(|_|DepositionError::Cancelled))?;
            let state=self.admit_thermal(&thermal_next,false,cx)?;
            let actual=self.mechanics(state,loads,config,cx)?;
            let applied_heat_j=finite(applied.iter().sum())?;
            let heat_defect_j=finite(applied.iter().zip(&actual.cells)
                .map(|(used,cell)|(used-cell.reversible_heat_j).abs()).sum())?;
            let external_heat_j=finite(audit.external_heat_j-applied_heat_j)?;
            let thermal_residual_j=finite(audit.stored_energy_change_j-audit.incoming_energy_j
                -external_heat_j-actual.reversible_heat_j)?;
            let combined_residual_j=finite(actual.internal_energy_change_j+audit.stored_energy_change_j
                -audit.incoming_energy_j-actual.external_work_j-external_heat_j)?;
            defects=[heat_defect_j,combined_residual_j];
            if heat_defect_j<=config.heat_tolerance_j
                && thermal_residual_j.abs()<=config.energy_tolerance_j
                && combined_residual_j.abs()<=config.energy_tolerance_j {
                let report=ThermoelasticReport {iterations:iteration,thermal_audit:audit,
                    cells:actual.cells,mechanical_work_j:actual.external_work_j,
                    internal_energy_change_j:actual.internal_energy_change_j,
                    reversible_heat_j:actual.reversible_heat_j,applied_heat_j,external_heat_j,
                    mechanical_residual_j:actual.mechanical_residual_j,thermal_residual_j,
                    combined_residual_j,heat_defect_j};
                let solid=Self {source:self.source,laws:self.laws.clone(),state:actual.state,
                    equilibrium:actual.equilibrium,potentials:actual.potentials,
                    fixed_dofs:self.fixed_dofs.clone(),nodal_forces_n:loads.nodal_forces_n().to_vec(),
                    solve:self.solve};
                check(cx)?;
                return Ok(ThermoelasticCandidate {thermal:thermal_next,solid,observer:coupled.inner,report});
            }
            for (value,cell) in applied.iter_mut().zip(&actual.cells) {
                *value=finite(config.heat_relaxation.mul_add(cell.reversible_heat_j-*value,*value))?;
            }
        }
        check(cx)?;
        Err(ThermoelasticFeedbackError::NotConverged {iterations:config.maximum_iterations,
            heat_defect_j:defects[0],combined_residual_j:defects[1]})
    }
}

fn potentials<E>(state:&SolidEnthalpyState<'_>,equilibrium:&SolidEquilibrium,
    laws:&[AffineIsotropicThermoelastic],config:ThermoelasticSolve,cx:&Cx<'_>)
    ->Result<Vec<ThermoelasticState>,ThermoelasticFeedbackError<E>> {
    let mut potentials=Vec::with_capacity(laws.len());
    for (i,(law,sample)) in laws.iter().zip(&state.samples).enumerate() {
        check(cx)?;
        let source=law.at_materials(&sample.elastic,sample.expansion.current())?;
        let native=&equilibrium.solution.stress.elements()[i];
        let free=state.strains[i].free_strain_mandel();
        let expected=source.free_strain()?;
        for k in 0..6 {
            let target=if k<3 {expected}else{0.0};
            if (free[k]-target).abs()>4096.0*f64::EPSILON
                *law.maximum_absolute_mandel_strain() {
                return Err(ThermoelasticFeedbackError::Invalid("native integrated expansion differs from the admitted potential"));
            }
        }
        let strain=core::array::from_fn(|k|native.elastic_strain_mandel()[k]+free[k]);
        let response=source.state(strain,sample.reference_volume_m3)?;
        if (response.free_energy_j-native.elastic_energy_j()).abs()>config.energy_tolerance_j
            || response.stress_mandel_pa.iter().zip(native.stress_mandel_pa())
                .any(|(a,b)|(a-b).abs()>config.stress_tolerance_pa) {
            return Err(ThermoelasticFeedbackError::Invalid("native stress/energy differs from the same original thermal material potential"));
        }
        potentials.push(response);
    }
    Ok(potentials)
}

struct MechanicalProposal<'m> {
    state:SolidEnthalpyState<'m>,equilibrium:SolidEquilibrium,potentials:Vec<ThermoelasticState>,
    cells:Vec<ThermoelasticCellStep>,external_work_j:f64,internal_energy_change_j:f64,
    reversible_heat_j:f64,mechanical_residual_j:f64,
}
#[derive(Clone,Debug,PartialEq)]
pub struct ThermoelasticCellStep {
    pub cell:CellId,
    pub temperature_k:[f64;2],
    pub potentials:[ThermoelasticState;2],
    pub mechanical_work_j:f64,
    pub expansion_work_j:f64,
    pub stiffness_work_j:f64,
    pub reversible_heat_j:f64,
    pub free_energy_residual_j:f64,
}
#[derive(Debug)]
pub struct ThermoelasticReport {
    pub iterations:usize,
    /// Original thermal audit includes the applied reversible source once.
    pub thermal_audit:EnergyAudit,
    pub cells:Vec<ThermoelasticCellStep>,
    pub mechanical_work_j:f64,
    pub internal_energy_change_j:f64,
    pub reversible_heat_j:f64,
    pub applied_heat_j:f64,
    pub external_heat_j:f64,
    pub mechanical_residual_j:f64,
    pub thermal_residual_j:f64,
    pub combined_residual_j:f64,
    pub heat_defect_j:f64,
}
pub struct ThermoelasticCandidate<'s,'m,'t,O> {
    pub thermal:DepositionSimulation<'t>,
    pub solid:ThermoelasticSolid<'s,'m>,
    pub observer:O,
    pub report:ThermoelasticReport,
}
impl<'s,'m,'t,O> ThermoelasticCandidate<'s,'m,'t,O> {
    pub fn into_parts(self)->(DepositionSimulation<'t>,ThermoelasticSolid<'s,'m>,O,ThermoelasticReport) {
        (self.thermal,self.solid,self.observer,self.report)
    }
}

#[derive(Clone)]
struct ReversibleHeat<'a,'s,'m,'cx,O> {
    solid:&'a ThermoelasticSolid<'s,'m>,inner:O,cx:&'a Cx<'cx>,watts:BTreeMap<CellId,f64>,
}
impl<O:DepositionObserver> DepositionObserver for ReversibleHeat<'_,'_,'_,'_,O> {
    type Error=ThermoelasticFeedbackError<O::Error>;
    fn begin(&mut self)->Result<(),Self::Error> {
        self.inner.begin().map_err(ThermoelasticFeedbackError::Observer)
    }
    fn requires_step_interval(&self)->bool {self.inner.requires_step_interval()}
    fn step_sources(&mut self,thermal:&DepositionSimulation<'_>,start:f64,dt:f64,end:f64,
        check:&mut dyn FnMut()->Result<(),DepositionError>)->Result<BTreeMap<CellId,f64>,Self::Error> {
        self.inner.step_sources(thermal,start,dt,end,check).map_err(ThermoelasticFeedbackError::Observer)
    }
    fn sources(&mut self,thermal:&DepositionSimulation<'_>,
        checkpoint:&mut dyn FnMut()->Result<(),DepositionError>)->Result<BTreeMap<CellId,f64>,Self::Error> {
        self.solid.admit_thermal(thermal,false,self.cx)?;
        let mut rates=self.inner.sources(thermal,checkpoint).map_err(ThermoelasticFeedbackError::Observer)?;
        for (&cell,&power) in &self.watts {
            checkpoint()?; let rate=rates.entry(cell).or_insert(0.0);
            *rate=finite(*rate+power)?;
        }
        Ok(rates)
    }
    fn conductances(&mut self,thermal:&DepositionSimulation<'_>,
        check:&mut dyn FnMut()->Result<(),DepositionError>)->Result<BTreeMap<u64,f64>,Self::Error> {
        self.inner.conductances(thermal,check).map_err(ThermoelasticFeedbackError::Observer)
    }
    fn heat_transfers(&mut self,thermal:&DepositionSimulation<'_>,
        check:&mut dyn FnMut()->Result<(),DepositionError>)->Result<Vec<HeatTransfer>,Self::Error> {
        self.inner.heat_transfers(thermal,check).map_err(ThermoelasticFeedbackError::Observer)
    }
    fn limit_timestep_s(&mut self,thermal:&DepositionSimulation<'_>,proposed:f64,
        check:&mut dyn FnMut()->Result<(),DepositionError>)->Result<f64,Self::Error> {
        self.inner.limit_timestep_s(thermal,proposed,check).map_err(ThermoelasticFeedbackError::Observer)
    }
    fn event(&mut self,event:&Event,thermal:&DepositionSimulation<'_>)->Result<(),Self::Error> {
        self.inner.event(event,thermal).map_err(ThermoelasticFeedbackError::Observer)?;
        self.solid.admit_thermal(thermal,false,self.cx).map(|_|())
    }
    fn step(&mut self,start:f64,dt:f64,old:&BTreeMap<CellId,f64>,thermal:&DepositionSimulation<'_>,
        check:&mut dyn FnMut()->Result<(),DepositionError>)->Result<(),Self::Error> {
        self.inner.step(start,dt,old,thermal,check).map_err(ThermoelasticFeedbackError::Observer)?;
        self.solid.admit_thermal(thermal,false,self.cx).map(|_|())
    }
}
