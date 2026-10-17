//! Same-tick hammer/damper contact and geometric extension, with one physical commit.
//! The extra conservative force is resolved against the original Schur image.
//! Hammer/jack prediction and acoustic/damper splitting occur once per tick,
//! not once per trial; all felt/Prony history stays frozen during these solves.
use super::{Instrument,Error,felt};
use super::super::linear::string_stretching;

/// Set FS_PIANO_CONTACT_TRACE=1 before the first contact solve. The cached opt-in
/// retains at most one worst normalized force residual per original sweep. No
/// output on success, no heap buffer, and no changes to any acceptance decision.
fn contact_trace_enabled()->bool {
    static ENABLED:std::sync::OnceLock<bool>=std::sync::OnceLock::new();
    *ENABLED.get_or_init(||std::env::var_os("FS_PIANO_CONTACT_TRACE")
        .is_some_and(|value|value=="1"))
}

#[derive(Clone,Copy)]
struct ContactSweepTrace {
    contact:usize,midi:u8,ratio:f64,force_n:f64,expected_n:f64,tolerance_n:f64,
    slope_n_m:f64,start_m:f64,end_m:f64,free_m:f64,compliance_m_n:f64,
    gap_m:f64,reconstructed_gap_m:f64,primitives:[f64;4],
    felt_parameters:[f64;6],felt_state:[f64;2],
    quadrature:(&'static str,bool),
    thickness_m:f64,area_m2:f64,rate_n_s_m_p:f64,dt_s:f64,
}
impl ContactSweepTrace {
    fn emit(&self,sweep:usize,out:&mut impl std::io::Write) {
        let Self {contact,midi,ratio,force_n,expected_n,tolerance_n,slope_n_m,
            start_m,end_m,free_m,compliance_m_n,gap_m,reconstructed_gap_m,
            primitives:[a,b,pa,pb],quadrature:(quadrature,legacy_midpoint_threshold),
            felt_parameters:[f_ref_pa,eps_ref,p,q,crush_fraction,eps_densify],
            felt_state:[eps_max,sig_max_pa],thickness_m,area_m2,rate_n_s_m_p,dt_s}=*self;
        let residual_n=(force_n-expected_n).abs();
        let delta_strain=b-a;let delta_primitive_pa=pb-pa;
        let midpoint=quadrature=="midpoint";
        let _=writeln!(out,"piano-contact-trace sweep={sweep} contact={contact} midi={midi} \
            ratio={ratio:.17e} force_n={force_n:.17e} expected_n={expected_n:.17e} \
            residual_n={residual_n:.17e} tolerance_n={tolerance_n:.17e} slope_n_m={slope_n_m:.17e} \
            start_m={start_m:.17e} end_m={end_m:.17e} free_m={free_m:.17e} \
            compliance_m_n={compliance_m_n:.17e} gap_m={gap_m:.17e} \
            reconstructed_gap_m={reconstructed_gap_m:.17e} strain_start={a:.17e} strain_end={b:.17e} \
            primitive_start_pa={pa:.17e} primitive_end_pa={pb:.17e} \
            delta_strain={delta_strain:.17e} delta_primitive_pa={delta_primitive_pa:.17e} midpoint={midpoint} \
            quadrature={quadrature} legacy_midpoint_threshold={legacy_midpoint_threshold} \
            f_ref_pa={f_ref_pa:.17e} eps_ref={eps_ref:.17e} p={p:.17e} q={q:.17e} \
            crush_fraction={crush_fraction:.17e} eps_densify={eps_densify:.17e} \
            eps_max={eps_max:.17e} sig_max_pa={sig_max_pa:.17e} \
            thickness_m={thickness_m:.17e} area_m2={area_m2:.17e} rate_n_s_m_p={rate_n_s_m_p:.17e} dt_s={dt_s:.17e}");
    }
}

impl Instrument {
    /// Attach supplied axial rigidities in this instrument's original scale.
    /// Preparation is cold, before excitation. No running stiffness changes,
    /// inferred winding rigidity, modal replacement or hidden retuning.
    pub fn configure_string_stretching(&mut self,spec:&string_stretching::Specification)->Result<(),String> {
        if self.accounting.input_work_j!=0. || self.hammers.iter().any(|h|h.active||h.held)
            || self.bank.q.iter().chain(&self.bank.v).any(|x|*x!=0.) {
            return Err("string stretching must be prepared before piano excitation".into());
        }
        self.bank.configure_string_stretching(&self.courses,spec)
    }

    // Called after the ONE hammer/jack prediction. Only numerical arrays change
    // here; the existing energy/material gate and output-sample rollback own
    // publication. Even a rejected final trial consumes no mechanical time.
    pub(super) fn solve_string_contact_step(&mut self)->Result<(),Error> {
        for _ in 0..string_stretching::MAX_TRIALS {
            self.solve_contacts_for_trial()?;
            let strings_converged=self.bank.correct_string_stretching_step().map_err(Error::Contact)?;
            let dampers_converged=if let Some(dampers)=&mut self.felt_dampers {
                let converged=dampers.correct(&self.bank.next_q).map_err(Error::Contact)?;
                if !converged {
                    self.bank.set_external_force(dampers.generalized_force()).map_err(Error::Contact)?;
                }
                converged
            }else{true};
            if strings_converged&&dampers_converged {return Ok(());}
            self.bank.predict();
        }
        Err(Error::Contact("coupled string/hammer/damper solve exhausted 64 trials"))
    }

    // Original contact equation, global 32-sweep limit, per-hammer blocks and
    // force tolerance. Re-evaluate against each candidate extension force rather
    // than accepting contacts against yesterday's tension or a split kick.
    fn solve_contacts_for_trial(&mut self)->Result<(),Error> {
        let nc=self.contacts.len();self.active.clear();
        self.force.fill(0.0);
        for c in 0..nc {
            let ci=self.bank.strings[self.bank.contact_strings[c]].course;
            self.gap[c]=self.contact_hammer_projection[c]*self.hammer_free[ci]-self.bank.free_contact[c]-self.contact_recession_m[c];
            if self.hammers[ci].active&&self.contacts[c].enabled {
                self.active.push(c);self.force[c]=self.contacts[c].force;
            }
        }
        for &i in &self.active {for &j in &self.active{self.gap[i]-=self.contact_h[i*nc+j]*self.force[j];}}
        let mut converged=self.active.is_empty();
        // Only the opt-in initializes snapshots or performs diagnostic arithmetic.
        let mut trace=contact_trace_enabled().then(||[None;32]);
        for sweep in 0..32 {
            let mut worst:Option<(usize,f64,f64,f64,f64)>=None;
            if converged{break;}
            self.contact_sweep()?;
            converged=true;
            for &i in &self.active {
                let ci=self.bank.strings[self.bank.contact_strings[i]].course;
                let material=&self.creep[i];let old=&self.contacts[i];
                let start=old.overlap-material.deformation(&old.memory);
                let end=self.gap[i]-material.free_deformation(&old.memory)-material.compliance()*self.force[i];
                let (expected,slope)=felt::average_with_rate(&self.laws[ci],&old.state,start,end,
                    self.contact_thickness_m[i],self.contact_areas[i],self.source_rate_n_s_m_p[i],
                    1.0/f64::from(self.bank.rate));
                // Scalar sites converge to what the felt law can resolve: 1e-9 N, or the
                // force that felt::solve_with_rate's rounding floor moves through the felt's
                // tangent (stiff with R_H). The string-stretching iteration around this solve
                // tests 1e-10 relative, so the old 1e-5 N gate stalled it. Sites solved as one
                // finite-footprint block keep that block solver's own 1e-5 N gate.
                let tolerance=if self.contact_block_selected(ci) {1e-5+1e-8*expected.abs()} else {
                    let diagonal=self.contact_h[i*nc+i];
                    let free=self.gap[i]+diagonal*self.force[i]-material.free_deformation(&old.memory);
                    let floor=16.0*f64::EPSILON*(end.abs()+((diagonal+material.compliance())*self.force[i]).abs()+free.abs());
                    1e-9+1e-11*expected.abs()+floor*slope.abs()
                };
                if !expected.is_finite()||(self.force[i]-expected).abs()>tolerance{converged=false;}
                if trace.is_some() {
                    let ratio=if expected.is_finite()&&self.force[i].is_finite()
                        &&tolerance.is_finite()&&tolerance>0.0 {
                        (self.force[i]-expected).abs()/tolerance
                    }else{f64::INFINITY};
                    if worst.is_none_or(|(_,r,_,_,_)|ratio>r) {
                        worst=Some((i,ratio,expected,slope,tolerance));
                    }
                }
            }
            if let (Some(trace),Some((i,ratio,expected_n,slope_n_m,tolerance_n)))=(&mut trace,worst) {
                let ci=self.bank.strings[self.bank.contact_strings[i]].course;
                let material=&self.creep[i];let old=&self.contacts[i];
                let start_m=old.overlap-material.deformation(&old.memory);
                let end_m=self.gap[i]-material.free_deformation(&old.memory)-material.compliance()*self.force[i];
                let mut reconstructed_gap_m=self.contact_hammer_projection[i]*self.hammer_free[ci]
                    -self.bank.free_contact[i]-self.contact_recession_m[i];
                for &j in &self.active {reconstructed_gap_m-=self.contact_h[i*nc+j]*self.force[j];}
                // Copy fixed scalar data only; retain the committed branch so a
                // failed quadrature can be reproduced without changing the solve.
                let law=&self.laws[ci];
                trace[sweep]=Some(ContactSweepTrace {contact:i,midi:self.courses[ci].midi,ratio,
                    force_n:self.force[i],expected_n,tolerance_n,slope_n_m,start_m,end_m,
                    free_m:self.gap[i]+self.contact_h[i*nc+i]*self.force[i]-material.free_deformation(&old.memory),
                    compliance_m_n:self.contact_h[i*nc+i]+material.compliance(),gap_m:self.gap[i],
                    reconstructed_gap_m,primitives:felt::diagnostic_primitives(law,&old.state,
                        start_m,end_m,self.contact_thickness_m[i]),
                    quadrature:felt::diagnostic_quadrature(law,&old.state,
                        start_m,end_m,self.contact_thickness_m[i]),
                    felt_parameters:[law.f_ref,law.eps_ref,law.p,law.q,law.crush_fraction,law.eps_densify],
                    felt_state:[old.state.eps_max,old.state.sig_max],
                    thickness_m:self.contact_thickness_m[i],area_m2:self.contact_areas[i],
                    rate_n_s_m_p:self.source_rate_n_s_m_p[i],dt_s:1.0/f64::from(self.bank.rate)});
            }
        }
        if !converged{
            if let Some(trace)=trace {
                use std::io::Write;
                let mut out=std::io::stderr().lock();
                let _=writeln!(out,"piano-contact-trace refusal sweeps=32 active_contacts={} mechanics_rate_hz={}",
                    self.active.len(),self.bank.rate);
                for (sweep,entry) in trace.iter().enumerate() {
                    if let Some(entry)=entry {entry.emit(sweep+1,&mut out);}
                }
            }
            return Err(Error::NoConvergence);
        }
        self.bank.finish(&self.force);
        Ok(())
    }
}

#[cfg(test)]
#[path="string_engine_tests.rs"]
mod tests;
