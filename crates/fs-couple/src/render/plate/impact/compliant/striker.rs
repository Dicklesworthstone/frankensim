//! Finite-area felt strikers in an existing physical coordinate layout.
//!
//! This is a cold permutation of MovingPads, not another material/contact law.
//! Every footprint site has independent felt and Kelvin history but shares ONE
//! inertial striker. The caller replaces its old hard-tip body/contact during
//! construction; installing both would double-count the same physical contact.
use super::{CompliantJaw, JawPort, MovingPads, PadSite};
use super::super::{ImpactBody, ImpactError, MAX_IMPACT_MODES, felt::FeltPad, invalid};

/// One replacement striker and its finite-area compression-only contacts.
/// Other mechanical coordinates keep their exact addresses. No state is changed
/// in an already-running instrument, and no extra acoustic source is introduced.
#[derive(Debug)]
pub struct FeltStriker {
    pub body: ImpactBody,
    pub pads: Vec<FeltPad>,
    pub port: JawPort,
}
impl FeltStriker {
    /// A flat (or otherwise authored) face above a nonconforming target.
    /// `clearances_m[i]` is the nonnegative EXTRA reference gap at site i,
    /// in addition to the common `jaw.initial_gap_m`. The caller derives it
    /// from geometry; it is not a material thickness or a contact-weight scale.
    ///
    /// All sites still share one inertia and independently retain the original
    /// WoolFelt/Kelvin histories. Only the affine contact offset changes. Zero
    /// clearances preserve `new` exactly, including its coordinate origin.
    ///
    /// # Errors
    /// Invalid shape, negative/nonfinite or overflowing gap, or any refusal
    /// from the existing finite-area striker and pad owners.
    pub fn new_with_clearances(total: usize, coordinate: usize, sites: &[PadSite],
        jaw: &CompliantJaw, clearances_m: &[f64]) -> Result<Self, ImpactError>
    {
        if clearances_m.len() != sites.len() || clearances_m.iter().any(|d|
            !d.is_finite() || *d < 0.0 || !(jaw.initial_gap_m + d).is_finite())
        {
            return Err(invalid("felt striker requires one finite nonnegative clearance per site"));
        }
        let mut result = Self::new(total, coordinate, sites, jaw)?;
        for (pad, &clearance) in result.pads.iter_mut().zip(clearances_m) {
            // Keep the common inertial origin at -initial_gap. Adding a site
            // offset MUST NOT move that origin, or the gap is counted twice.
            if clearance != 0.0 { pad.precompression_m = -clearance; }
            pad.validate(total)?;
        }
        Ok(result)
    }

    /// `sites` use the complete target layout, including a reserved, ZERO
    /// `coordinate` entry. The supplied jaw owns effective mass, launch speed,
    /// felt properties and whole-footprint Kelvin elements. A negative-side jaw
    /// strikes a positive-downward head with positive force. Positive-side jaws
    /// use inward-positive motion too; their physical world-axis sign reverses.
    ///
    /// Striker q is measured from the uncompressed contact plane: initial tip
    /// displacement is -gap. MovingPads instead uses zero initial travel and an
    /// explicit gap, so this adapter translates the coordinate origin along
    /// with the pad offset. Compression, kinetic energy and force work agree.
    /// The original four-site and total-mode budgets still apply.
    pub fn new(total: usize, coordinate: usize, sites: &[PadSite], jaw: &CompliantJaw)
        -> Result<Self, ImpactError>
    {
        if total < 2 || total > MAX_IMPACT_MODES || coordinate >= total {
            return Err(invalid("felt striker requires a reserved coordinate and a bounded surface basis"));
        }
        if sites.len() > super::MAX_PAD_SITES || sites.iter().any(|s|
            s.weights.len() != total || s.weights[coordinate] != 0.0)
        {
            return Err(invalid("felt footprint must not couple the surface to its reserved striker"));
        }
        let surface: Vec<_> = sites.iter().map(|s| {
            let mut weights=s.weights.clone();weights.remove(coordinate);
            PadSite { weights, area_m2:s.area_m2 }
        }).collect();
        // Reuse all existing physical admission, area partitioning and material
        // lowering. The temporary order is [all non-striker modes, striker].
        let mut compiled=MovingPads::new(total-1,&surface,std::slice::from_ref(jaw))?;
        let mut body=compiled.bodies.remove(0);
        let mut port=compiled.ports[0];port.coordinate=coordinate;
        let position=-jaw.initial_gap_m/port.inverse_sqrt_mass;
        if !position.is_finite() || (jaw.initial_gap_m>0.0 && position==0.0) {
            return Err(invalid("felt striker initial gap is not representable"));
        }
        body.initial[0].displacement_m_sqrt_kg=position;
        for pad in &mut compiled.pads {
            // Exact permutation; no re-normalization of a head's mode shape.
            pad.weights[coordinate..].rotate_right(1);
            pad.precompression_m=0.0;
        }
        Ok(Self { body, pads:compiled.pads, port })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::PadSide;
    use super::super::super::{BodyPotential, ImpactConfig, ImpactSystem};
    use fs_exec::CancelGate;
    use fs_material::fiber::WoolFelt;
    use crate::modal_acoustic_time::ModalAcousticState;
    fn jaw() -> CompliantJaw {
        CompliantJaw { side:PadSide::Negative,mass_kg:0.02,drag_n_s_m:0.0,
            initial_gap_m:0.00002,initial_velocity_m_s:0.5,thickness_m:0.006,
            law:WoolFelt::new(1e6,0.2,2.2,3.0,0.15,0.7).unwrap(),
            prior_maximum_strain:0.0,creep:vec![] }
    }
    fn config() -> ImpactConfig {
        ImpactConfig { dt_s:2e-6,max_steps:2000,maximum_energy_j:1.0,
            energy_absolute_tolerance_j:1e-10,energy_relative_tolerance:1e-7,
            maximum_generalized_force:1e4 }
    }
    #[test]
    fn relocation_preserves_every_surface_coefficient_and_physical_gap() {
        for coordinate in 0..4 {
            let mut row=vec![2.0,-3.0,4.0,5.0];row[coordinate]=0.0;
            let sites=[PadSite {weights:row.clone(),area_m2:0.0001}];
            let j=jaw();let compiled=FeltStriker::new(4,coordinate,&sites,&j).unwrap();
            let weight=compiled.port.inverse_sqrt_mass;
            assert_eq!(compiled.port.coordinate,coordinate);
            for (i,&b) in compiled.pads[0].weights.iter().enumerate() {
                assert_eq!(b,if i==coordinate {weight}else{-row[i]});
            }
            let q=compiled.body.initial[0].displacement_m_sqrt_kg;
            assert!((q*weight+j.initial_gap_m).abs()<1e-19);
            let v=compiled.body.initial[0].velocity_m_sqrt_kg_per_s;
            assert!((0.5*v*v-0.5*j.mass_kg*j.initial_velocity_m_s.powi(2)).abs()<1e-16);
            assert_eq!(compiled.pads[0].precompression_m,0.0);
        }
    }
    #[test]
    fn finite_footprint_does_not_average_opposite_surface_motion_before_contact() {
        let mut j=jaw();j.initial_gap_m=0.0;j.initial_velocity_m_s=0.0;
        let a=FeltStriker::new(3,0,&[
            PadSite {weights:vec![0.0,1.0,1.0],area_m2:0.0002},
            PadSite {weights:vec![0.0,1.0,-1.0],area_m2:0.0002}],&j).unwrap();
        let surface=ImpactBody {potential:BodyPotential::Linear(vec![0.0;2]),
            initial:vec![ModalAcousticState::default();2],damping_per_s:vec![0.0;2]};
        let s=ImpactSystem::new(vec![a.body,surface],vec![],a.pads,vec![],config()).unwrap();
        let mut x=s.state().to_vec();x[4]=0.0002;
        let mut g=vec![0.0;x.len()];use fs_phs::Storage;
        s.mechanical.gradient(&x,&mut g);
        assert!(s.mechanical.hamiltonian(&x)>0.0 && g[4]>0.0 && g[0]>0.0);
        // A single mean row [0,1,0] would produce exactly zero contact here.
    }
    #[test]
    fn felt_strike_exchanges_impulse_and_retains_history_across_rejected_steps() {
        let make=|| {
            let (target,weight)=ImpactBody::free_mass(0.1,0.0,0.0).unwrap();
            let a=FeltStriker::new(2,0,&vec![PadSite {weights:vec![0.0,weight],area_m2:0.0001};4],&jaw()).unwrap();
            ImpactSystem::new(vec![a.body,target],vec![],a.pads,vec![],config()).unwrap()
                .prepare_analytic().unwrap()
        };
        let mut a=make();let mut clean=make();let gate=CancelGate::new_clock_free();
        let initial=a.stored_energy_j();let mut loss=0.0;
        for tick in 0..1200 {
            if tick==300 {
                let state=a.state().to_vec();let history=a.felt_history(0).unwrap();
                let cancel=CancelGate::new_clock_free();cancel.request();
                assert!(a.step(&[0.0;2],&cancel).is_err());
                assert!(a.step(&[1e9,0.0],&gate).is_err());
                assert_eq!(a.state(),state);assert_eq!(a.felt_history(0).unwrap(),history);
            }
            let frame=a.step(&[0.0;2],&gate).unwrap();clean.step(&[0.0;2],&gate).unwrap();
            assert_eq!(a.state(),clean.state());loss+=frame.dissipated_energy_j;
            assert!((a.stored_energy_j()+loss-initial).abs()<1e-7);
            let momentum=0.02_f64.sqrt()*a.state()[1]+0.1_f64.sqrt()*a.state()[3];
            assert!((momentum-0.01).abs()<1e-8);
        }
        assert!(a.state()[2]>0.0 && loss>0.0 && a.felt_history(0).unwrap().eps_max>0.0);
    }
    #[test]
    fn missing_surface_self_contact_and_excess_footprint_refuse() {
        let j=jaw();let site=PadSite {weights:vec![0.0,2.0],area_m2:0.001};
        assert!(FeltStriker::new(1,0,&[],&j).is_err());
        assert!(FeltStriker::new(2,2,&[site.clone()],&j).is_err());
        assert!(FeltStriker::new(2,0,&vec![site.clone();5],&j).is_err());
        for row in [vec![1.0,2.0],vec![0.0],vec![0.0,0.0],vec![0.0,f64::NAN]] {
            assert!(FeltStriker::new(2,0,&[PadSite {weights:row,..site.clone()}],&j).is_err());
        }
    }
}

#[cfg(test)]
mod profiled_tests {
    use super::*;
    use super::super::{PadSide, super::{ImpactConfig, ImpactSystem, felt::KelvinBranch}};
    use fs_exec::CancelGate;
    use fs_material::fiber::WoolFelt;
    use fs_phs::Storage;

    fn jaw() -> CompliantJaw {
        CompliantJaw { side:PadSide::Positive, mass_kg:0.02, drag_n_s_m:0.0,
            initial_gap_m:0.00002, initial_velocity_m_s:0.5, thickness_m:0.006,
            law:WoolFelt::new(1e6,0.2,2.2,3.0,0.15,0.7).unwrap(),
            prior_maximum_strain:0.0, creep:vec![] }
    }
    fn cfg() -> ImpactConfig {
        ImpactConfig { dt_s:2e-6,max_steps:2000,maximum_energy_j:1.0,
            energy_absolute_tolerance_j:1e-10,energy_relative_tolerance:1e-7,
            maximum_generalized_force:1e4 }
    }
    fn sites() -> Vec<PadSite> {
        vec![PadSite { weights:vec![0.0,2.0],area_m2:0.0001 },
            PadSite { weights:vec![0.0,2.0],area_m2:0.0003 }]
    }
    #[test]
    fn local_clearances_preserve_one_mass_area_and_material_time_constants() {
        let mut j=jaw();j.creep.push(KelvinBranch { stiffness_n_m:1000.,viscosity_n_s_m:2. });
        let base=FeltStriker::new(2,0,&sites(),&j).unwrap();
        let same=FeltStriker::new_with_clearances(2,0,&sites(),&j,&[0.,0.]).unwrap();
        let profiled=FeltStriker::new_with_clearances(2,0,&sites(),&j,&[0.,0.001]).unwrap();
        assert_eq!(base.body.initial[0].displacement_m_sqrt_kg,profiled.body.initial[0].displacement_m_sqrt_kg);
        assert_eq!(base.body.initial[0].velocity_m_sqrt_kg_per_s,profiled.body.initial[0].velocity_m_sqrt_kg_per_s);
        for ((a,b),c) in base.pads.iter().zip(&same.pads).zip(&profiled.pads) {
            assert_eq!(a.precompression_m.to_bits(),b.precompression_m.to_bits());
            assert_eq!(a.weights,c.weights);assert_eq!(a.area_m2,c.area_m2);
            assert_eq!(a.creep[0].stiffness_n_m,c.creep[0].stiffness_n_m);
            assert_eq!(a.creep[0].viscosity_n_s_m,c.creep[0].viscosity_n_s_m);
        }
        assert_eq!(profiled.pads[1].precompression_m,-0.001);
        for invalid in [vec![],vec![0.,-1.],vec![0.,f64::NAN],vec![f64::INFINITY,0.]] {
            assert!(FeltStriker::new_with_clearances(2,0,&sites(),&j,&invalid).is_err());
        }
    }
    #[test]
    fn a_recessed_site_does_not_contact_early_or_duplicate_the_high_site_force() {
        let mut j=jaw();j.initial_velocity_m_s=0.;
        let make=|sites:&[PadSite],gaps:&[f64]| {
            let a=FeltStriker::new_with_clearances(2,0,sites,&j,gaps).unwrap();
            let target=ImpactBody::free_mass(0.25,0.,0.).unwrap().0;
            ImpactSystem::new(vec![a.body,target],vec![],a.pads,vec![],cfg()).unwrap()
        };
        let a=make(&sites(),&[0.,0.001]);let single=make(&sites()[..1],&[0.]);
        let mut x=a.state().to_vec();x[0]=0.0002*0.02_f64.sqrt();
        let mut g=vec![0.;x.len()];let mut one=g.clone();
        a.mechanical.gradient(&x,&mut g);single.mechanical.gradient(&x,&mut one);
        assert_eq!(g,one);assert_eq!(a.mechanical.hamiltonian(&x),single.mechanical.hamiltonian(&x));
        assert!(g[0]>0. && g[2]>0.);
        assert!((g[0]*0.02_f64.sqrt()-g[2]*0.25_f64.sqrt()).abs()<1e-12);
    }
    #[test]
    fn profiled_felt_collision_has_reciprocal_impulse_loss_and_exact_retry() {
        let make=|| {
            let a=FeltStriker::new_with_clearances(2,0,&sites(),&jaw(),&[0.,0.0001]).unwrap();
            let target=ImpactBody::free_mass(0.25,0.,0.).unwrap().0;
            ImpactSystem::new(vec![a.body,target],vec![],a.pads,vec![],cfg()).unwrap().prepare_analytic().unwrap()
        };
        let mut a=make();let mut clean=make();let initial=a.stored_energy_j();
        let gate=CancelGate::new_clock_free();let mut loss=0.;
        for tick in 0..1000 {
            if tick==400 {
                let old=a.state().to_vec();let h=a.felt_history(0).unwrap();
                let cancelled=CancelGate::new_clock_free();cancelled.request();
                assert!(a.step(&[0.;2],&cancelled).is_err());
                assert_eq!(a.state(),old);assert_eq!(a.felt_history(0).unwrap(),h);
            }
            let f=a.step(&[0.;2],&gate).unwrap();clean.step(&[0.;2],&gate).unwrap();
            assert_eq!(a.state(),clean.state());loss+=f.dissipated_energy_j;
            assert!((f.stored_energy_j+loss-initial).abs()<1e-7);
            // Jaw axis is down and target axis up: physical signed momentum.
            assert!((a.state()[1]*0.02_f64.sqrt()-a.state()[3]*0.25_f64.sqrt()-0.01).abs()<1e-8);
        }
        assert!(loss>0. && a.state()[2]<0.);
        assert!(a.felt_history(0).unwrap().eps_max>0.);
    }
}
