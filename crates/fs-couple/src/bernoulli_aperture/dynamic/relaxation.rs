//! Regional hereditary bending projected through the actual retained plate mode.
//! Equilibrium stiffness remains in the plate; only additional Maxwell arms enter
//! this state. Prestress, pressure area and the lay law are never relaxed here.
use super::{DynamicAperture, invalid};
use crate::acoustic_realize::AcousticRealizeError;
use fs_material::visco::GeneralizedMaxwell;
use fs_phs::RelaxationBranch;
use fs_plate::{PlateSection, PlateSectionField, dkt_stiffness};

/// One isotropic equilibrium section and its explicitly supplied relaxation law.
#[derive(Clone, Debug)]
pub struct PlateRelaxationRegion {
    /// Unique original triangle indices; regions must partition the whole plate.
    pub triangles: Vec<usize>,
    /// E_inf must reproduce the source section's actual equilibrium D matrix.
    /// E_j and tau are supplied material data, not inferred from a damping ratio.
    pub material: GeneralizedMaxwell,
    /// Constant Poisson ratio of this proportional isotropic bending law.
    pub poisson_ratio: f64,
    /// Declared material-use band [Hz]. Checks cover the retained equilibrium and
    /// instantaneous frequencies, NOT the spectrum of arbitrary transient forcing.
    pub band_hz: (f64, f64),
    /// Source/condition or an explicit authored-model label.
    pub provenance: String,
}

/// Complete material map and explicit temporal accuracy admission.
#[derive(Clone, Debug)]
pub struct PlateRelaxationSpec {
    /// In source order; every triangle appears exactly once.
    pub regions: Vec<PlateRelaxationRegion>,
    /// At most 64 supplied regional terms. Zero-modulus terms carry no state;
    /// positive arms are never truncated.
    pub max_branches: usize,
    /// dt/tau ceiling in (0,2]. Two is the nonoscillatory midpoint-memory limit;
    /// a smaller value is a caller-selected resolution allowance, not a certificate.
    pub max_dt_over_tau: f64,
    /// dt*omega_instantaneous ceiling in (0,1], chosen by the caller.
    pub max_angular_step: f64,
}

/// Initial viscous deformation is part of the physical state, not a solver seed.
#[derive(Clone, Debug)]
pub enum InitialApertureMemory {
    /// Each arm is relaxed at the supplied current plate displacement.
    Relaxed,
    /// Zero prior viscous displacement, which can store energy in a bent plate.
    Unrelaxed,
    /// One equivalent opening-coordinate viscous displacement [m] per positive
    /// arm, in region/term order. Uniform history within each regional arm.
    ViscousDisplacementM(Vec<f64>),
}

/// Retained source map and accepted material history. No mutable state escapes.
#[derive(Clone, Debug)]
pub struct ApertureRelaxation {
    spec: PlateRelaxationSpec,
    branches: Vec<RelaxationBranch>,
    region_stiffness_n_m: Vec<f64>,
    z: Vec<f64>,
    energy_j: f64,
    dissipated_j: f64,
}
impl ApertureRelaxation {
    /// Exact supplied laws, regions, source labels and numerical allowances.
    #[must_use]
    pub const fn spec(&self) -> &PlateRelaxationSpec { &self.spec }
    /// Discrete regional phi^T K_bending phi in opening coordinates [N/m].
    /// This excludes membrane prestress and does not modify equilibrium storage.
    #[must_use]
    pub fn region_bending_stiffness_n_m(&self) -> &[f64] { &self.region_stiffness_n_m }
    /// Positive projected relaxing arms; projection=[1] in opening metres.
    #[must_use]
    pub fn branches(&self) -> &[RelaxationBranch] { &self.branches }
    /// Accepted energy-normalized viscous coordinates z [sqrt(J)].
    #[must_use]
    pub fn memory_sqrt_j(&self) -> &[f64] { &self.z }
    /// Current arm storage [J], already included in total aperture/tube energy.
    #[must_use]
    pub const fn stored_energy_j(&self) -> f64 { self.energy_j }
    /// Accumulated material loss since construction [J], not contact or jet loss.
    #[must_use]
    pub const fn dissipated_energy_j(&self) -> f64 { self.dissipated_j }

    pub(super) fn restoring_force(&self, old: f64, next: f64, dt: f64)
        -> Result<f64, AcousticRealizeError>
    {
        let mut force = 0.0;
        for (arm, &z) in self.branches.iter().zip(&self.z) {
            force += arm.scalar_midpoint_step(old,next,z,dt).map_err(map)?.1;
        }
        if !force.is_finite() { return Err(invalid("aperture material restoring force overflowed")); }
        Ok(force)
    }
    pub(super) fn preview(&self, old: f64, next: f64, dt: f64)
        -> Result<MemoryTrial, AcousticRealizeError>
    {
        let mut z_next=Vec::with_capacity(self.z.len());
        let (mut energy,mut loss)=(0.0,0.0);
        for (arm,&z) in self.branches.iter().zip(&self.z) {
            let (new,_,stored,dissipated)=arm.scalar_midpoint_step(old,next,z,dt).map_err(map)?;
            z_next.push(new);energy+=stored;loss+=dissipated;
        }
        let cumulative=self.dissipated_j+loss;
        if ![energy,loss,cumulative].iter().all(|x| x.is_finite()) {
            return Err(invalid("aperture material storage or loss overflowed"));
        }
        Ok(MemoryTrial {z_next,energy,loss,cumulative})
    }
    pub(super) fn accept(&mut self, trial: MemoryTrial) {
        self.z=trial.z_next;self.energy_j=trial.energy;self.dissipated_j=trial.cumulative;
    }
}
pub(super) struct MemoryTrial {
    z_next: Vec<f64>,
    pub(super) energy: f64,
    pub(super) loss: f64,
    cumulative: f64,
}
fn map(e: fs_phs::PhsError) -> AcousticRealizeError { AcousticRealizeError::Nonlinear(e.to_string()) }

impl DynamicAperture {
    /// Attach a supplied regional Maxwell spectrum to a cold plate-bound valve.
    /// The same momentum/contact/pressure solve includes its midpoint restoring
    /// force; material history publishes only with the complete coupled step.
    /// Compatible with spatial closure, ordinary tubes, branched networks and
    /// their existing passive boundary/cavity/wall loads.
    ///
    /// Requires zero separate modal damping to avoid counting material loss twice.
    /// The actual DKT matrices project each region; no stiffness, area, mass or
    /// resonance is retuned. Proportional isotropic bending and one structural
    /// mode only; neither prestress nor contact stiffness receives Maxwell arms.
    /// This is supplied constitutive data, not automatic wet-cane identification.
    ///
    /// # Errors
    /// Missing plate, repeated/late attachment, unmatched equilibrium sections,
    /// incomplete regions, invalid material/history, or exceeded state/time/band
    /// budgets. No positive material arm is dropped to admit a specimen.
    pub fn with_plate_relaxation(mut self, spec: PlateRelaxationSpec, initial: InitialApertureMemory)
        -> Result<Self, AcousticRealizeError>
    {
        if self.accepted_steps!=0 || self.relaxation.is_some() || self.spec.damping_ratio!=0.0 {
            return Err(invalid("plate relaxation requires one cold attachment and zero separate modal damping"));
        }
        let plate=self.plate.as_ref().ok_or_else(||invalid("material bending requires the retained physical plate"))?;
        let chart=plate.chart();let count=chart.mesh.tris.len();
        if spec.regions.is_empty() || spec.regions.len()>count || spec.max_branches>64
            || !spec.max_dt_over_tau.is_finite() || spec.max_dt_over_tau<=0.0 || spec.max_dt_over_tau>2.0
            || !spec.max_angular_step.is_finite() || spec.max_angular_step<=0.0 || spec.max_angular_step>1.0
        {return Err(invalid("plate relaxation requires complete regions and explicit bounded time/state allowances"));}
        let total=spec.regions.iter().try_fold(0usize,|sum,r|sum.checked_add(r.material.terms.len()))
            .ok_or_else(||invalid("material branch count overflow"))?;
        if total>spec.max_branches {return Err(invalid("plate material exceeds its branch allowance"));}
        let mut covered=vec![false;count];let mut region_stiffness=Vec::new();let mut branches=Vec::new();
        for region in &spec.regions {
            if region.provenance.trim().is_empty() || region.triangles.is_empty() || region.triangles.len()>count
                || !region.band_hz.0.is_finite() || !region.band_hz.1.is_finite()
                || region.band_hz.0<0.0 || region.band_hz.1<=region.band_hz.0 {
                return Err(invalid("material region requires explicit source, finite band and nonempty triangle set"));
            }
            // The model has public fields: re-admit them before use.
            GeneralizedMaxwell::new(region.material.e_inf,region.material.terms.clone())
                .map_err(|e|AcousticRealizeError::Nonlinear(e.to_string()))?;
            let mut stiffness=0.0;
            for &i in &region.triangles {
                if i>=count || covered[i] {return Err(invalid("material regions must cover each source triangle exactly once"));}
                covered[i]=true;
                let section=match chart.section_field() {
                    PlateSectionField::Uniform(s)=>s,PlateSectionField::PerElement(s)=>&s[i],
                };
                let equilibrium=PlateSection::isotropic(region.material.e_inf,region.poisson_ratio,
                    section.thickness,section.density).map_err(AcousticRealizeError::Plate)?;
                if equilibrium.d!=section.d {return Err(invalid("Maxwell equilibrium modulus and Poisson ratio must match the actual isotropic plate section"));}
                let tri=chart.mesh.tris[i];
                let x=tri.map(|n|chart.mesh.nodes[n].0);let y=tri.map(|n|chart.mesh.nodes[n].1);
                let (ke,_)=dkt_stiffness(&x,&y,&section.d,i).map_err(AcousticRealizeError::Plate)?;
                let mut q=[0.0;9];
                for (local,&node) in tri.iter().enumerate() {q[3*local..3*local+3].copy_from_slice(&plate.shape_per_opening()[node]);}
                let mut element=0.0;
                for a in 0..9 {for b in 0..9 {element+=q[a]*ke[9*a+b]*q[b];}}
                if !element.is_finite() || element<0.0 {return Err(invalid("regional projected bending energy is not finite nonnegative"));}
                stiffness+=element;
            }
            if !stiffness.is_finite() || stiffness<=0.0 {return Err(invalid("material region has no resolved positive bending energy"));}
            region_stiffness.push(stiffness);
            for &(modulus,tau) in &region.material.terms {
                if modulus==0.0 {continue;} // exactly zero is no physical arm
                if self.spec.time_step_s/tau>spec.max_dt_over_tau {return Err(invalid("material relaxation is unresolved at the admitted mechanical clock"));}
                let arm=RelaxationBranch {projection:vec![1.0],stiffness:stiffness*(modulus/region.material.e_inf),relaxation_time_s:tau};
                arm.scalar_relaxed_memory(0.0).map_err(map)?;branches.push(arm);
            }
        }
        if covered.iter().any(|x|!*x) {return Err(invalid("material map omits plate triangles"));}
        let instantaneous=self.spec.stiffness_n_m+branches.iter().map(|a|a.stiffness).sum::<f64>();
        let omega=(instantaneous/self.spec.mass_kg).sqrt();
        if !omega.is_finite() || omega*self.spec.time_step_s>spec.max_angular_step {
            return Err(invalid("instantaneous material stiffness exceeds the declared time resolution"));
        }
        let equilibrium_hz=(self.spec.stiffness_n_m/self.spec.mass_kg).sqrt()/core::f64::consts::TAU;
        let instant_hz=omega/core::f64::consts::TAU;
        if spec.regions.iter().any(|r|equilibrium_hz<r.band_hz.0 || instant_hz>r.band_hz.1) {
            return Err(invalid("retained equilibrium or instantaneous mode exceeds a supplied material band"));
        }
        let displacement=self.state.opening_m-self.spec.aperture.rest_opening_m;
        if let InitialApertureMemory::ViscousDisplacementM(values)=&initial {
            if values.len()!=branches.len() {return Err(invalid("initial material history must name every positive regional arm"));}
        }
        let mut z=Vec::with_capacity(branches.len());let mut energy=0.0;
        for (i,arm) in branches.iter().enumerate() {
            let viscous=match &initial {
                InitialApertureMemory::Relaxed=>displacement,InitialApertureMemory::Unrelaxed=>0.0,
                InitialApertureMemory::ViscousDisplacementM(values)=>values[i],
            };
            let value=arm.scalar_relaxed_memory(viscous).map_err(map)?;
            energy+=arm.scalar_stored_energy(displacement,value).map_err(map)?;z.push(value);
        }
        if !energy.is_finite() || !(energy+self.energy_at(self.state)).is_finite() {
            return Err(invalid("initial material energy overflow"));
        }
        self.relaxation=Some(ApertureRelaxation {spec,branches,region_stiffness_n_m:region_stiffness,
            z,energy_j:energy,dissipated_j:0.0});
        Ok(self)
    }

    /// Accepted material state and source map, absent on the original elastic path.
    #[must_use]
    pub fn relaxation(&self) -> Option<&ApertureRelaxation> {self.relaxation.as_ref()}
}
