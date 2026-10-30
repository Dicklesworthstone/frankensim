//! Restrained uniaxial thermoviscoelastic members in series with one cohesive joint.
//!
//! Temperature and total end-to-end extension are inputs; joint opening and force
//! are solved outputs. For each member, the existing backward-Euler material owner
//! gives sigma_i=C_i epsilon_i+b_i. Eliminating the common axial force gives
//! S F + g = D - free, S=sum L_i/(A_i C_i), free=-sum L_i b_i/C_i.
//! Thus the existing cohesive equilibrium owner solves g/S + A_b t(g)=(D-free)/S.
//! No additional constitutive law, finite-difference tangent or root solver is used.
//!
//! This is a small-strain, quasistatic axial coupon/restraint, not a 3D peel or
//! growing-mesh solver. Member temperatures are uniform; all members and the bond
//! already exist at initialization. Frames, reference areas and lengths stay fixed.
//! Thermal transport and contact updates belong to the caller. No fracture heat,
//! healing, fitted adhesion, dynamic snap-back or equilibrium uniqueness is claimed.

/// Original bulk and ordered interface sources on the same joint history.
pub mod source;

use std::collections::BTreeMap;
use fs_exec::Cx;
use fs_material::cohesive::{BilinearCohesive, CohesiveEnergy, CohesivePoint};
use fs_material::thermovisco::{Energy, Point, ThermoMaxwell};
use fs_solver::{LinearOp, NewtonReport};
use super::{CohesivePatch, CohesiveState, CohesiveSystem, JumpTerm};

/// One uniform axial member with its own immutable material and birth reference.
#[derive(Debug, Clone, PartialEq)]
pub struct Member<'a> {
    id: u64,
    length_m: f64,
    area_m2: f64,
    point: Point<'a>,
}
impl<'a> Member<'a> {
    /// Stress-free birth at zero total strain; lengths/areas are reference values.
    pub fn new(id: u64, length_m: f64, area_m2: f64, law: &'a ThermoMaxwell,
        birth_temperature_k: f64) -> Result<Self, Error> {
        positive(length_m)?; positive(area_m2)?; positive(length_m*area_m2)?;
        let point = law.birth(birth_temperature_k,0.0)
            .map_err(|error| Error::Material { member: id,error })?;
        Ok(Self { id,length_m,area_m2,point })
    }
    pub fn id(&self) -> u64 { self.id }
    pub fn reference_length_m(&self) -> f64 { self.length_m }
    pub fn reference_area_m2(&self) -> f64 { self.area_m2 }
    pub fn point(&self) -> &Point<'a> { &self.point }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Budget {
    pub max_members: usize,
    pub max_branch_states: usize,
    pub max_krylov_basis_entries: usize,
}

#[derive(Debug, Clone, Copy)]
pub struct Config {
    /// Existing cohesive force, local energy, iteration and displacement controls.
    pub equilibrium: super::Config,
    /// Absolute end-to-end extension compatibility gate, m.
    pub extension_tolerance_m: f64,
    /// Absolute combined discrete work/storage gate, J; separate from accuracy.
    pub energy_tolerance_j: f64,
}

#[derive(Debug, Clone)]
pub enum Error {
    Invalid(&'static str),
    Budget,
    Cancelled,
    Material { member: u64, error: fs_material::thermovisco::Error },
    Equilibrium(super::Error),
    Compatibility(f64),
    MemberForce { member: u64, residual_n: f64 },
    Energy(f64),
}
impl From<super::Error> for Error {
    fn from(error: super::Error) -> Self {
        match error { super::Error::Cancelled => Self::Cancelled, other => Self::Equilibrium(other) }
    }
}
impl core::fmt::Display for Error {
    fn fmt(&self,f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f,"thermal cohesive joint refused: {self:?}")
    }
}
impl std::error::Error for Error {}

/// Accepted physical history. Clone is an in-memory checkpoint, not a codec.
#[derive(Debug, Clone, PartialEq)]
pub struct Joint<'a> {
    time_s: f64,
    extension_m: f64,
    members: BTreeMap<u64,Member<'a>>,
    interface: u64,
    area_m2: f64,
    bond: CohesivePoint<'a>,
    budget: Budget,
}

/// Actual extensive endpoint accounts. All energy quantities are joules.
#[derive(Debug, Clone)]
pub struct Report {
    pub force_n: f64,
    pub opening_m: f64,
    pub extension_residual_m: f64,
    pub maximum_force_residual_n: f64,
    /// Constitutive densities integrated over each original member volume.
    pub bulk_energy_j: Energy,
    pub cohesive_energy_j: CohesiveEnergy,
    /// Physical viscosity of the surviving bond, separate from intrinsic fracture
    /// and bulk Maxwell loss. No heat partition or deposition is inferred.
    pub bond_viscous_dissipated_j: f64,
    /// Endpoint-force quadrature F_new*(D_new-D_old), not exact actuator work.
    pub endpoint_work_j: f64,
    /// F_new*(g_new-g_old) minus exactly integrated cohesive path work. Signed
    /// quadrature correction, NOT physical dissipation or a passivity claim.
    pub cohesive_work_correction_j: f64,
    pub energy_residual_j: f64,
    pub nonlinear: NewtonReport,
    pub operator_evaluations: usize,
    pub material_point_updates: usize,
}

struct SeriesStiffness(f64);
impl LinearOp for SeriesStiffness {
    fn n(&self) -> usize { 1 }
    fn apply(&self,x: &[f64],y: &mut [f64]) { y[0]=self.0*x[0]; }
}

impl<'a> Joint<'a> {
    /// Initialize a bonded, stress-free assembly at a nonnegative physical time.
    /// Birth temperatures may differ; each member retains its own free-strain
    /// reference. No bond is inferred from a thermal contact or material name.
    pub fn new(time_s: f64, members: Vec<Member<'a>>, interface: u64, area_m2: f64,
        law: &'a BilinearCohesive, budget: Budget) -> Result<Self, Error> {
        if !time_s.is_finite() || time_s < 0.0 { return Err(Error::Invalid("invalid start time")); }
        positive(area_m2)?;
        if members.is_empty() || members.len() > budget.max_members { return Err(Error::Budget); }
        let mut ordered = BTreeMap::new();
        let mut branches = 0usize;
        for member in members {
            branches = branches.checked_add(member.point.law().branch_count()).ok_or(Error::Budget)?;
            if branches > budget.max_branch_states { return Err(Error::Budget); }
            if ordered.insert(member.id,member).is_some() { return Err(Error::Invalid("duplicate member ID")); }
        }
        Ok(Self { time_s,extension_m: 0.0,members: ordered,interface,area_m2,
            bond: law.bonded_point(),budget })
    }
    pub fn time_s(&self) -> f64 { self.time_s }
    pub fn extension_m(&self) -> f64 { self.extension_m }
    pub fn interface_id(&self) -> u64 { self.interface }
    pub fn bond(&self) -> CohesivePoint<'a> { self.bond }
    pub fn members(&self) -> &BTreeMap<u64,Member<'a>> { &self.members }
    pub fn stored_energy_j(&self) -> f64 {
        self.members.values().map(|m| m.length_m*m.area_m2*m.point.stored_energy_j_m3()).sum::<f64>()
            +self.area_m2*self.bond.stored_energy_j_m2()
    }

    /// Endpoint temperature is used by the existing implicit material update.
    /// Subdivide thermal/load paths for accuracy and at reversals. A changing
    /// stiffness never resets branch strain or cohesive damage. The temporary
    /// reduced operator is immutable during each native equilibrium solve.
    /// All member and bond histories publish together only after force, length,
    /// constitutive and combined-energy checks. Cancellation refuses the whole
    /// interval. Only bounded constructor/clone operations lack internal polls.
    pub fn advance_controlled(&self, target_s: f64, extension_m: f64,
        temperatures_k: &BTreeMap<u64,f64>, config: Config,
        cancelled: &mut impl FnMut() -> bool) -> Result<(Self,Report),Error> {
        poll(cancelled)?;
        let dt=finite(target_s-self.time_s)?; positive(dt)?; finite(extension_m)?;
        positive(config.extension_tolerance_m)?; positive(config.energy_tolerance_j)?;
        if temperatures_k.len()!=self.members.len() { return Err(Error::Invalid("one temperature per member required")); }
        let mut affine=Vec::with_capacity(self.members.len());
        let mut compliance=0.0;
        let mut free=0.0;
        for (&id,member) in &self.members {
            poll(cancelled)?;
            let temperature=*temperatures_k.get(&id).ok_or(Error::Invalid("missing member temperature"))?;
            let law=member.point.linearize(temperature,dt).map_err(|error| Error::Material { member: id,error })?;
            let s=member.length_m/member.area_m2/law.tangent_pa; positive(s)?;
            compliance=finite(compliance+s)?;
            free=finite(free-member.length_m*(law.offset_pa/law.tangent_pa))?;
            affine.push((id,temperature,law));
        }
        let stiffness=SeriesStiffness(1.0/compliance); positive(stiffness.0)?;
        let drive=finite(stiffness.0*finite(extension_m-free)?)?;
        let patch=CohesivePatch::new(self.interface,self.area_m2,self.bond.law(),
            vec![JumpTerm { dof: 0,local: [1.0,0.0,0.0] }])?;
        let system=CohesiveSystem::new(&stiffness,vec![patch],super::Budget {
            max_dofs: 1,max_patches: 1,max_jump_terms: 1,
            max_krylov_basis_entries: self.budget.max_krylov_basis_entries })?;
        // Bind the SAME physical point to the new algorithmic series operator.
        // This private adapter cannot replace laws, damage or geometric maps.
        let old=CohesiveState { system: &system,displacements: vec![self.bond.jump_m()[0]],
            points: vec![self.bond] };
        let solved=old.solve_over_time_controlled(&[drive],&[[0.0;3]],dt,config.equilibrium,cancelled)?;
        let gap=solved.state.displacements_m()[0];
        // Preserve the original law lifetime; no second equilibrium is solved.
        let bond=self.bond.trial_over_time([gap,0.0,0.0],dt,config.equilibrium.local_energy_tolerance_j_m2)
            .map_err(|error| super::Error::Material { patch: self.interface,error })?;
        if bond.equilibrium.point!=solved.state.points()[0] { return Err(Error::Invalid("inconsistent native cohesive endpoint")); }
        let force=finite(self.area_m2*bond.traction_pa[0])?;
        let mut next=self.clone();
        let mut energy=Energy::default();
        let mut actual_extension=gap;
        let mut maximum_force_residual_n=solved.maximum_force_residual_n;
        for (id,temperature,law) in affine {
            poll(cancelled)?;
            let member=&self.members[&id];
            let strain=finite((force/member.area_m2-law.offset_pa)/law.tangent_pa)?;
            let volume=member.length_m*member.area_m2;
            let local_gate=config.energy_tolerance_j/self.members.len() as f64/volume;
            let step=member.point.advance(temperature,strain,dt,local_gate)
                .map_err(|error| Error::Material { member: id,error })?;
            let force_residual=finite(member.area_m2*step.stress_pa-force)?;
            maximum_force_residual_n=maximum_force_residual_n.max(force_residual.abs());
            if force_residual.abs()>config.equilibrium.force_tolerance_n {
                return Err(Error::MemberForce { member: id,residual_n: force_residual });
            }
            actual_extension=finite(actual_extension+member.length_m*strain)?;
            accumulate(&mut energy,step.energy,volume)?;
            next.members.get_mut(&id).expect("retained member").point=step.point;
        }
        let extension_residual=finite(actual_extension-extension_m)?;
        if extension_residual.abs()>config.extension_tolerance_m { return Err(Error::Compatibility(extension_residual)); }
        next.time_s=target_s; next.extension_m=extension_m; next.bond=bond.equilibrium.point;
        let cohesive=solved.cohesive_energy_j;
        let viscous=solved.viscous_dissipated_j;
        let cohesive_work=finite(cohesive.work+viscous)?;
        let work=finite(force*(extension_m-self.extension_m))?;
        let correction=finite(force*(gap-self.bond.jump_m()[0])-cohesive_work)?;
        let storage=finite(next.stored_energy_j()-self.stored_energy_j())?;
        let residual=finite(storage+energy.dissipated+energy.numerical_dissipation
            +cohesive.dissipated+viscous+correction-work-energy.thermal_strain_work-energy.modulus_work)?;
        let conjugacy=finite(work-correction-energy.mechanical_work-cohesive_work)?;
        let storage_check=finite(storage-energy.stored_change-cohesive.stored_change)?;
        for defect in [residual,conjugacy,storage_check] {
            if defect.abs()>config.energy_tolerance_j { return Err(Error::Energy(defect)); }
        }
        poll(cancelled)?;
        Ok((next,Report { force_n: force,opening_m: gap,extension_residual_m: extension_residual,
            maximum_force_residual_n,bulk_energy_j: energy,cohesive_energy_j: cohesive,bond_viscous_dissipated_j: viscous,
            endpoint_work_j: work,cohesive_work_correction_j: correction,energy_residual_j: residual,
            nonlinear: solved.nonlinear,operator_evaluations: solved.operator_evaluations,
            material_point_updates: self.members.len() }))
    }

    /// Same operation within the existing executor cancellation scope.
    pub fn advance_to(&self,cx: &Cx<'_>,target_s: f64,extension_m: f64,
        temperatures_k: &BTreeMap<u64,f64>,config: Config) -> Result<(Self,Report),Error> {
        self.advance_controlled(target_s,extension_m,temperatures_k,config,
            &mut || cx.checkpoint().is_err())
    }
}

fn accumulate(sum: &mut Energy,part: Energy,volume: f64) -> Result<(),Error> {
    for (target,value) in [(&mut sum.mechanical_work,part.mechanical_work),
        (&mut sum.thermal_strain_work,part.thermal_strain_work),(&mut sum.modulus_work,part.modulus_work),
        (&mut sum.stored_change,part.stored_change),(&mut sum.dissipated,part.dissipated),
        (&mut sum.numerical_dissipation,part.numerical_dissipation),(&mut sum.residual,part.residual)] {
        *target=finite(*target+volume*value)?;
    }
    Ok(())
}
fn finite(value: f64) -> Result<f64,Error> {
    if value.is_finite() { Ok(value) } else { Err(Error::Invalid("nonfinite joint arithmetic")) }
}
fn positive(value: f64) -> Result<(),Error> {
    if value.is_finite() && value>0.0 { Ok(()) } else { Err(Error::Invalid("positive finite data required")) }
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(),Error> {
    if cancelled() { Err(Error::Cancelled) } else { Ok(()) }
}

#[cfg(test)]
mod tests;

#[cfg(test)]
#[path = "thermal_joint/source_tests.rs"]
mod source_tests;
