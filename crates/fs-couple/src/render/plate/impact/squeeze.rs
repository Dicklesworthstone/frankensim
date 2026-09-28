//! Reynolds-film resistance in the EXISTING joint mechanical time equation.
//! fs-tribo owns the pressure law; fs-phs owns stepping and dissipation receipts.
//! Contact, film pressure and its state-dependent tangent are evaluated at the
//! same trial configuration and discrete-gradient effort, never one tick late.
pub use fs_tribo::resistive_film::{FilmCell, FilmChannel, FilmError, FilmLimits,
    FilmReport, GapPort, ResistiveFilm};
use super::{ImpactError, ImpactSystem, MAX_IMPACT_MODES, invalid};
use fs_tribo::resistive_film::MAX_CELLS;

impl ImpactSystem {
    /// Add one immutable fluid pressure island without replacing contacts,
    /// resetting body/history state, or creating a second time integrator.
    /// Construction is cold. The film's ports must cover every mechanical mode,
    /// including exact zeros for unrelated bodies. At most four islands.
    pub fn with_squeeze_film(mut self, film: ResistiveFilm) -> Result<Self, ImpactError> {
        if self.gas_film.is_some() || film.port_count() != self.modes || self.squeeze_films.len() >= 4 {
            return Err(invalid("squeeze film mode count, duplicate gas image or four-island ceiling failed"));
        }
        let mut q = [0.0; MAX_IMPACT_MODES];
        for i in 0..self.modes { q[i] = self.x[2*i]; }
        let mut force = [0.0; MAX_IMPACT_MODES]; let mut pressure = [0.0; MAX_CELLS];
        film.evaluate_into(&q[..self.modes], &[0.0; MAX_IMPACT_MODES][..self.modes],
            &mut force[..self.modes], &mut pressure[..film.cell_count()])
            .map_err(|e| invalid(e.0))?;
        self.squeeze_films.push(film);
        Ok(self)
    }

    pub(super) fn has_nonlinear_dissipation(&self) -> bool {
        self.contact_loss || !self.squeeze_films.is_empty() || self.gas_film.is_some()
    }

    // Endpoint admission is necessary too: an admissible midpoint must not
    // publish collapsed/trapped fluid at the end of an otherwise accepted tick.
    pub(super) fn validate_squeeze_configuration(&self, x: &[f64]) -> Result<(), ImpactError> {
        if let Some(gas)=&self.gas_film {gas.observe(x)?;}
        if self.squeeze_films.is_empty() { return Ok(()); }
        if x.len() < 2*self.modes { return Err(invalid("invalid squeeze state prefix")); }
        let mut q = [0.0; MAX_IMPACT_MODES]; let zero = [0.0; MAX_IMPACT_MODES];
        let mut force = [0.0; MAX_IMPACT_MODES]; let mut p = [0.0; MAX_CELLS];
        for i in 0..self.modes { q[i] = x[2*i]; }
        for film in &self.squeeze_films {
            film.evaluate_into(&q[..self.modes], &zero[..self.modes],
                &mut force[..self.modes], &mut p[..film.cell_count()])
                .map_err(|e| invalid(e.0))?;
        }
        Ok(())
    }

    pub(super) fn dissipative_flow_into(&self, x: &[f64], e: &[f64], out: &mut [f64]) -> bool {
        // The contact owner validates shapes and fills all entries, preserving
        // its existing Hunt--Crossley reaction and zeros on memory coordinates.
        if !self.contact.dissipative_flow_into(x,e,out) { return false; }
        if self.gas_film.as_ref().is_some_and(|gas|!gas.add_flow(x,e,out)) {return false;}
        if self.squeeze_films.is_empty() { return true; }
        let mut q=[0.0; MAX_IMPACT_MODES]; let mut v=[0.0; MAX_IMPACT_MODES];
        let mut f=[0.0; MAX_IMPACT_MODES]; let mut p=[0.0; MAX_CELLS];
        for i in 0..self.modes { q[i]=x[2*i]; v[i]=e[2*i+1]; }
        for film in &self.squeeze_films {
            if film.evaluate_into(&q[..self.modes], &v[..self.modes],
                &mut f[..self.modes], &mut p[..film.cell_count()]).is_err() { return false; }
            for i in 0..self.modes { out[2*i+1] += f[i]; }
        }
        out.iter().all(|x| x.is_finite())
    }

    pub(super) fn dissipative_flow_tangent_into(&self, x: &[f64], e: &[f64],
        dx: &[f64], de: &[f64], out: &mut [f64]) -> bool
    {
        if !self.contact.dissipative_flow_tangent_into(x,e,dx,de,out) { return false; }
        if self.gas_film.as_ref().is_some_and(|gas|!gas.add_flow_tangent(x,e,dx,de,out)) {return false;}
        if self.squeeze_films.is_empty() { return true; }
        let mut q=[0.0; MAX_IMPACT_MODES]; let mut v=[0.0; MAX_IMPACT_MODES];
        let mut dq=[0.0; MAX_IMPACT_MODES]; let mut dv=[0.0; MAX_IMPACT_MODES];
        let mut f=[0.0; MAX_IMPACT_MODES];
        for i in 0..self.modes { q[i]=x[2*i]; v[i]=e[2*i+1]; dq[i]=dx[2*i]; dv[i]=de[2*i+1]; }
        for film in &self.squeeze_films {
            if film.tangent_into(&q[..self.modes], &v[..self.modes], &dq[..self.modes],
                &dv[..self.modes], &mut f[..self.modes]).is_err() { return false; }
            for i in 0..self.modes { out[2*i+1] += f[i]; }
        }
        out.iter().all(|x| x.is_finite())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::{ImpactBody, ImpactConfig};
    use fs_exec::CancelGate;

    fn pair() -> ImpactSystem {
        let (a,wa)=ImpactBody::free_mass(0.2,0.0,0.01).unwrap();
        let (b,wb)=ImpactBody::free_mass(0.4,0.0,0.0).unwrap();
        let g=GapPort { reference_m:0.001, closure:vec![wa,-wb] };
        let film=ResistiveFilm::new(vec![FilmCell { area_m2:0.002,gap:g.clone() }],
            vec![FilmChannel { from:0,to:None,width_m:0.05,length_m:0.01,gap:g }],2,1.8e-5,
            FilmLimits { minimum_cell_gap_m:1e-6,maximum_gap_m:0.01,maximum_pressure_pa:1000.0 }).unwrap();
        ImpactSystem::new(vec![a,b],vec![],vec![],vec![],ImpactConfig {
            dt_s:1e-5,max_steps:256,maximum_energy_j:1.0,energy_absolute_tolerance_j:1e-12,
            energy_relative_tolerance:1e-7,maximum_generalized_force:1.0 }).unwrap()
            .with_squeeze_film(film).unwrap()
    }
    fn momentum(x: &[f64]) -> f64 { 0.2_f64.sqrt()*x[1]+0.4_f64.sqrt()*x[3] }

    #[test]
    fn air_drives_resting_body_and_reference_energy_accounts_for_the_loss() {
        let mut s=pair();let initial=s.stored_energy_j();let p0=momentum(s.state());let mut loss=0.0;
        for _ in 0..128 {let f=s.step(&[0.0;2],&CancelGate::new_clock_free()).unwrap();
            assert!(f.dissipated_energy_j>=0.0);loss+=f.dissipated_energy_j;}
        assert!(s.state()[3]>0.0);assert!(s.stored_energy_j()<initial);
        assert!((momentum(s.state())-p0).abs()<1e-10);
        assert!((s.stored_energy_j()+loss-initial).abs()<1e-10);
    }

    #[test]
    fn prepared_images_keep_pressure_in_the_joint_equation_and_retry_atomically() {
        let mut reference=pair();
        for _ in 0..128 {reference.step(&[0.0;2],&CancelGate::new_clock_free()).unwrap();}
        for analytic in [false,true] {
            let mut p=if analytic {pair().prepare_analytic().unwrap()}else{pair().prepare().unwrap()};
            let initial=p.stored_energy_j();let mut loss=0.0;
            for tick in 0..128 {
                if tick==32 {
                    let before=p.state().to_vec();let samples=p.samples();
                    assert!(p.step(&[2.0,0.0],&CancelGate::new_clock_free()).is_err());
                    assert_eq!(p.state(),before);assert_eq!(p.samples(),samples);
                }
                loss+=p.step(&[0.0;2],&CancelGate::new_clock_free()).unwrap().dissipated_energy_j;
            }
            for (a,b) in p.state().iter().zip(reference.state()) {assert!((a-b).abs()<1e-8);}
            assert!((p.stored_energy_j()+loss-initial).abs()<1e-10);
        }
    }

    #[test]
    fn mechanical_prefix_only_and_endpoint_gap_refusal() {
        let s=pair();let mut x=s.state().to_vec();x.extend([0.3,-0.2]);
        let e=vec![0.0,0.01,0.0,0.0,0.7,-0.8];let mut out=vec![99.0;6];
        assert!(s.dissipative_flow_into(&x,&e,&mut out));
        assert_eq!(out[0],0.0);assert_eq!(out[2],0.0);assert_eq!(&out[4..],&[0.0;2]);
        x[0]=0.01;
        assert!(s.validate_squeeze_configuration(&x).is_err());
        assert_eq!(s.samples(),0);
    }
}
