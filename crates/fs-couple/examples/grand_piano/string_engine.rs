//! Same-tick contact and geometric extension, with one physical commit.
//! The extra conservative force is resolved against the original Schur image.
//! Hammer/jack prediction and acoustic/damper splitting occur once per tick,
//! not once per trial; all felt/Prony history stays frozen during these solves.
use super::{Instrument,Error,felt};
use super::super::linear::string_stretching;

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
            if self.bank.correct_string_stretching_step().map_err(Error::Contact)? {return Ok(());}
            self.bank.predict();
        }
        Err(Error::Contact("geometric string/contact solve exhausted 64 trials"))
    }

    // Original contact equation, global 32-sweep limit, per-hammer blocks and
    // force tolerance. Re-evaluate against each candidate extension force rather
    // than accepting contacts against yesterday's tension or a split kick.
    fn solve_contacts_for_trial(&mut self)->Result<(),Error> {
        let nc=self.contacts.len();self.active.clear();
        self.force.fill(0.0);
        for c in 0..nc {
            let ci=self.bank.strings[self.bank.contact_strings[c]].course;
            self.gap[c]=self.hammer_free[ci]-self.bank.free_contact[c];
            if self.hammers[ci].active&&self.contacts[c].enabled {
                self.active.push(c);self.force[c]=self.contacts[c].force;
            }
        }
        for &i in &self.active {for &j in &self.active{self.gap[i]-=self.contact_h[i*nc+j]*self.force[j];}}
        let mut converged=self.active.is_empty();
        for _ in 0..32 {
            if converged{break;}
            self.contact_sweep()?;
            converged=true;
            for &i in &self.active {
                let ci=self.bank.strings[self.bank.contact_strings[i]].course;let c=self.courses[ci];
                let material=&self.creep[i];let old=&self.contacts[i];
                let start=old.overlap-material.deformation(&old.memory);
                let end=self.gap[i]-material.free_deformation(&old.memory)-material.compliance()*self.force[i];
                let expected=felt::average(&self.laws[ci],&old.state,start,end,
                    c.felt_thickness_m,self.contact_areas[i]).0;
                if !expected.is_finite()||(self.force[i]-expected).abs()>1e-5+1e-8*expected.abs(){converged=false;}
            }
        }
        if !converged{return Err(Error::NoConvergence);}
        self.bank.finish(&self.force);
        Ok(())
    }
}

#[cfg(test)]
#[path="string_engine_tests.rs"]
mod tests;
