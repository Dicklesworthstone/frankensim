//! One supplied chamber, its actual head, distributed pressure and felt striker.
//! Only composition: fs-plate/fs-modal own reduction, fs-material owns felt,
//! and fs-phs owns the prepared nonlinear time step.
use super::{spec::Spec, Error};
use fs_couple::{modal_acoustic_time::ModalAcousticState,
    render::plate::impact::{BodyPotential, ImpactBody, ImpactConfig, SubsteppedImpactSystem,
        cavity::{CavityCoupling, meridian::{Budget, MeridianCavity}},
        compliant::{CompliantJaw, PadSide, PadSite, striker::FeltStriker},
        felt::KelvinBranch, membrane::MembranePotential},
    vibroacoustic::{AcousticMedium, CavityModes}};
use fs_exec::CancelGate;
use fs_material::fiber::WoolFelt;
use fs_plate::{ModePair, shell::{head::TensionedDisk, profile::ProfileBudget,
    head::nonlinear::MembraneReductionBudget}};
use std::{collections::BTreeMap, f64::consts::{PI, TAU}};

pub struct Model {
    pub mechanics: SubsteppedImpactSystem,
    pub cavity: CavityCoupling,
    pub air: MeridianCavity,
    pub film: TensionedDisk,
    pub modes: Vec<ModePair>,
    pub tip_weight: f64,
    pub pickup: Vec<f64>,
    pub force: Vec<f64>,
    tensioner: Option<super::tuning::Tensioner>,
    edge_port: Option<(usize,f64)>,
}
pub fn rim(film: &TensionedDisk) -> Vec<(usize, usize)> {
    let mut edges = BTreeMap::<(usize, usize), Vec<(usize, usize)>>::new();
    for t in &film.mesh.tris { for k in 0..3 {
        let (a,b) = (t[k],t[(k+1)%3]);
        edges.entry((a.min(b),a.max(b))).or_default().push((a,b));
    }}
    edges.values().filter(|e| e.len() == 1).map(|e| e[0]).collect()
}
/// Whole-disk containment in the actual convex polygon, not just four samples.
pub fn footprint(film: &TensionedDisk, centre: [f64;2], radius: f64) -> Result<Vec<[f64;2]>, Error> {
    if centre.iter().any(|v| !v.is_finite()) || !radius.is_finite() || radius <= 0.0 {
        return Err("invalid mallet footprint".into());
    }
    let edges = rim(film);
    if edges.len() < 3 { return Err("head has no closed rim".into()); }
    for (a,b) in edges {
        let (x,y) = film.mesh.nodes[a]; let (u,v) = film.mesh.nodes[b];
        let distance = ((u-x)*(centre[1]-y) - (v-y)*(centre[0]-x))/(u-x).hypot(v-y);
        if !distance.is_finite() || distance <= radius { return Err("entire mallet face must lie strictly inside the actual head rim".into()); }
    }
    let r = radius / 2.0_f64.sqrt();
    Ok(vec![[centre[0]+r,centre[1]], [centre[0],centre[1]+r],
        [centre[0]-r,centre[1]], [centre[0],centre[1]-r]])
}
/// Signed OUTWARD head displacement; the mechanical coordinate is downward.
/// Air shapes, norms and overlaps all come from the SAME meridional basis.
pub fn interface(film: &TensionedDisk, modes: &[ModePair], air: &MeridianCavity,
    gate: &CancelGate) -> Result<(CavityModes, Vec<f64>), Error> {
    use fs_couple::render::plate::impact::cavity::meridian::end_disk::{
        DiskPort, EndDisk, InterfaceBudget};
    let interface = air.end_disk_interface(film,modes,DiskPort {
        end: EndDisk::Upper, first_mode: 1, structural_modes: modes.len()+1, positive_z: -1.0 },
        InterfaceBudget { maximum_points: 12_288, maximum_terms: 1_000_000 },gate)?;
    Ok((interface.cavity,interface.overlaps))
}
impl Model {
    pub fn new(spec: &Spec, steps: u64, gate: &CancelGate) -> Result<Self, Error> {
        Self::with_exterior(spec,steps,gate,|_,_| Ok(None))
    }
    /// Prepare the acoustic scene before playing, then bind its SAME modal
    /// addresses to the shared reciprocal radiation-storage owner.
    pub fn with_exterior<F>(spec: &Spec, steps: u64, gate: &CancelGate, bake: F) -> Result<Self, Error>
    where F: FnOnce(&TensionedDisk, &[ModePair])
        -> Result<Option<fs_couple::render::plate::impact::radiation::Model>, Error>
    {
        Self::with_tensioner(spec,steps,gate,None,bake)
    }
    pub fn with_tensioner<F>(spec:&Spec,steps:u64,gate:&CancelGate,
        tensioner:Option<super::tuning::Tensioner>,bake:F)->Result<Self,Error>
    where F:FnOnce(&TensionedDisk,&[ModePair])
        ->Result<Option<fs_couple::render::plate::impact::radiation::Model>,Error>
    {
        Self::with_execution(spec,steps,gate,tensioner,super::execution::Settings::default(),bake)
    }
    pub fn with_execution<F>(spec:&Spec,steps:u64,gate:&CancelGate,
        tensioner:Option<super::tuning::Tensioner>,execution:super::execution::Settings,bake:F)->Result<Self,Error>
    where F:FnOnce(&TensionedDisk,&[ModePair])
        ->Result<Option<fs_couple::render::plate::impact::radiation::Model>,Error>
    {
        execution.validate()?;
        spec.validate()?;
        if let Some(t)=&tensioner {t.admit(steps)?;}
        if steps == 0 || gate.is_requested() { return Err("positive step budget and uncancelled preparation required".into()); }
        let film = TensionedDisk::new(spec.head, ProfileBudget {
            max_nodes: 2048, max_triangles: 4096, max_feature_evaluations: 100_000 })?;
        let modes = fs_modal::slice_window(&film.model.k, &film.model.m,
            ((TAU*spec.head_band_hz[0]).powi(2),(TAU*spec.head_band_hz[1]).powi(2)),
            &fs_plate::SliceOptions::default())?.modes;
        if modes.is_empty() || modes.len() > 32 { return Err("head window must contain 1..32 complete modes; none are truncated".into()); }
        let head_end=modes.len()+1;
        let n=head_end+usize::from(tensioner.is_some());
        let air = MeridianCavity::new(spec.profile.clone(), AcousticMedium { rho0: spec.air[0], c0: spec.air[1] },
            Budget { radial_intervals: spec.air_mesh[0], maximum_azimuthal_order: spec.air_mesh[1],
                maximum_frequency_hz: spec.air_maximum_hz, maximum_modes: 256, maximum_nodes: 512,
                maximum_matrix_entries: 2_000_000, relative_tolerance: 1e-7 }, gate)?;
        let (basis,mut overlaps) = interface(&film,&modes,&air,gate)?;
        // Fixed-reference radial edge motion sweeps no cavity volume directly.
        // Preserve all old head rows; append a zero actuator row BEFORE gas inertia.
        overlaps.resize(n*basis.omegas.len(),0.0);
        let drag: Vec<_> = basis.omegas.iter().map(|w| if *w == 0.0 {0.0} else {spec.air[2]}).collect();
        let coupling = CavityCoupling::new(&basis,n,&overlaps,&drag)?;
        let mut fixed: Vec<_> = rim(&film).iter().flat_map(|&(a,b)| [a,b]).collect();
        fixed.sort_unstable(); fixed.dedup();
        let budget=MembraneReductionBudget {max_modes:32,max_nodes:2048,max_facet_pairs:2_000_000,
            max_solve_entries:16_000_000,relative_tolerance:1e-6};
        let mut damping:Vec<_>=modes.iter().map(|m|2.0*spec.head_damping_ratio*m.lambda.sqrt()).collect();
        let potential=match &tensioner {
            Some(t)=>{
                let law=fs_couple::render::plate::impact::dilation::DilationPotential::from_disk(&film,&modes,t.config,budget)?;
                damping.push(t.drag_n_s_m/t.config.effective_mass_kg);
                BodyPotential::Dilation(law)
            },
            None=>BodyPotential::Membrane(MembranePotential::from_pencil(&film.mesh,&film.section,
                &film.model,&modes,&fixed,budget,spec.maximum_slope)?),
        };
        let head=ImpactBody {potential,initial:vec![ModalAcousticState::default();n-1],damping_per_s:damping};
        let m = spec.mallet; let f = spec.felt;
        let points = footprint(&film,[m[5],m[6]],m[1])?;
        let face = fs_couple::render::plate::impact::linear::wire::film_shapes(&film,&modes,&points)?;
        let sites:Vec<_>=face.into_iter().map(|shapes| {
            let mut weights:Vec<_>=std::iter::once(0.0).chain(shapes).collect();weights.resize(n,0.0);
            PadSite {weights,area_m2:PI*m[1]*m[1]/4.0}
        }).collect();
        let jaw = CompliantJaw { side: PadSide::Negative, mass_kg: m[0], drag_n_s_m: 0.0,
            initial_gap_m: m[3], initial_velocity_m_s: m[4], thickness_m: m[2],
            law: WoolFelt::new(f[0],f[1],f[2],f[3],f[4],f[5])?, prior_maximum_strain: f[6],
            creep: spec.creep.iter().map(|v| KelvinBranch { stiffness_n_m: v[0], viscosity_n_s_m: v[1] }).collect() };
        let mallet = FeltStriker::new(n,0,&sites,&jaw)?;
        let tip_weight = mallet.port.inverse_sqrt_mass;
        let (reference,cavity) = coupling.build(vec![mallet.body,head],vec![],mallet.pads,
            ImpactConfig { dt_s: 1.0/f64::from(spec.rate_hz), max_steps: steps, maximum_energy_j: 20.0,
                energy_absolute_tolerance_j: 1e-9, energy_relative_tolerance: 1e-6,
                maximum_generalized_force: 1e6 },gate)?;
        let reference = match bake(&film,&modes)? {
            Some(load) => cavity.attach_exterior_radiation(reference,&load,&(1..head_end).collect::<Vec<_>>(),
                fs_couple::render::plate::impact::radiation::MAX_JOINT_STATES/2)?,
            None => reference,
        };
        if gate.is_requested() { return Err("kettle preparation cancelled before publication".into()); }
        let edge_port=reference.dilation_force_port(1);
        if tensioner.is_some()!=edge_port.is_some() {return Err("lost edge force coordinate during composition".into());}
        let force = vec![0.0;cavity.total_modes()];
        let pickup = fs_couple::render::plate::impact::linear::wire::film_shapes(&film,&modes,&[[m[5],m[6]]])?.remove(0);
        Ok(Self { mechanics: execution.prepare(reference)?, cavity, air, film, modes, tip_weight, pickup, force, tensioner, edge_port })
    }
    /// Exact-sample force assignment enters the SAME accepted mechanical tick.
    /// Force scratch is restored on both success and refusal; sample/history
    /// publication remains owned by the shared impact solver, never by a scheduler.
    pub fn step(&mut self,gate:&CancelGate)->Result<fs_couple::render::plate::impact::ImpactFrame,Error> {
        if let (Some(t),Some((index,weight)))=(&self.tensioner,self.edge_port) {
            let previous=self.force[index];
            self.force[index]=previous+t.force_n(self.mechanics.samples())*weight;
            let result=self.mechanics.step(&self.force,gate);self.force[index]=previous;
            Ok(result?)
        } else {Ok(self.mechanics.step(&self.force,gate)?)}
    }
    pub fn head_displacement(&self) -> f64 {
        self.pickup.iter().enumerate().map(|(j,b)| b*self.mechanics.state()[2*(j+1)]).sum()
    }
}
