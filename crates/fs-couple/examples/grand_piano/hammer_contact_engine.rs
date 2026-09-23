//! Per-hammer blocks inside the existing global contact iteration.
//! Cross-key soundboard reactions still use the complete original compliance.
use super::{Instrument, Error, Course, Bank, felt};
use felt::block::{Site, Workspace, MAX_SITES};

pub(super) struct Prepared { workspace: Workspace, finite: Vec<bool> }
impl Prepared {
    pub(super) fn new(bank:&Bank,courses:&[Course])->Result<Option<Self>,String> {
        let mut count=vec![0usize;courses.len()];
        for &si in &bank.contact_strings {count[bank.strings[si].course]+=1;}
        if count.iter().any(|n|*n>MAX_SITES) {return Err("per-hammer simultaneous contact budget exceeded".into());}
        let finite:Vec<_>=count.iter().zip(courses).map(|(n,c)|*n>usize::from(c.unison)).collect();
        if !finite.iter().any(|v|*v) {return Ok(None);}
        Ok(Some(Self {workspace:Workspace::new()?,finite}))
    }
}

impl Instrument {
    pub(super) fn contact_sweep(&mut self)->Result<(),Error> {
        let mut first=0;
        while first<self.active.len() {
            let ci=self.bank.strings[self.bank.contact_strings[self.active[first]]].course;
            let mut end=first+1;
            while end<self.active.len()
                && self.bank.strings[self.bank.contact_strings[self.active[end]]].course==ci {end+=1;}
            if end-first>1 && self.contact_solver.as_ref().is_some_and(|s|s.finite[ci]) {
                self.contact_block(first,end,ci)?;
            } else {
                // Keep the original point/unselected arithmetic and ordering.
                for index in first..end {self.scalar_contact(index,ci)?;}
            }
            first=end;
        }
        Ok(())
    }

    fn scalar_contact(&mut self,index:usize,ci:usize)->Result<(),Error> {
        let nc=self.contacts.len();let i=self.active[index];
        let c=self.courses[ci];let diagonal=self.contact_h[i*nc+i];
        let material=&self.creep[i];let old=&self.contacts[i];
        let start=old.overlap-material.deformation(&old.memory);
        let free=self.gap[i]+diagonal*self.force[i]-material.free_deformation(&old.memory);
        let next=felt::solve(&self.laws[ci],&old.state,start,free,
            diagonal+material.compliance(),c.felt_thickness_m,self.contact_areas[i]).map_err(Error::Contact)?;
        let change=next-self.force[i];self.force[i]=next;
        for &j in &self.active {self.gap[j]-=self.contact_h[j*nc+i]*change;}
        Ok(())
    }

    fn contact_block(&mut self,first:usize,end:usize,ci:usize)->Result<(),Error> {
        let n=end-first;let nc=self.contacts.len();
        if n>MAX_SITES {return Err(Error::Contact("simultaneous hammer site budget exceeded"));}
        let c=self.courses[ci];let initial=self.active[first];
        let mut sites=[Site {law:&self.laws[ci],history:&self.contacts[initial].state,
            start_m:0.,thickness_m:c.felt_thickness_m,area_m2:self.contact_areas[initial]};MAX_SITES];
        let mut a=[0.;MAX_SITES*MAX_SITES];let mut free=[0.;MAX_SITES];
        let mut warm=[0.;MAX_SITES];let mut forces=[0.;MAX_SITES];
        for row in 0..n {
            let i=self.active[first+row];let material=&self.creep[i];let old=&self.contacts[i];
            sites[row]=Site {law:&self.laws[ci],history:&old.state,
                start_m:old.overlap-material.deformation(&old.memory),
                thickness_m:c.felt_thickness_m,area_m2:self.contact_areas[i]};
            free[row]=self.gap[i]-material.free_deformation(&old.memory);warm[row]=self.force[i];
            for col in 0..n {
                let j=self.active[first+col];let h=self.contact_h[i*nc+j];
                free[row]+=h*self.force[j];a[row*n+col]=h;
            }
            // Prony strain is local. Its free history is above, and only its
            // exact held-force response belongs on this diagonal, once.
            a[row*n+row]+=material.compliance();
        }
        self.contact_solver.as_mut().ok_or(Error::Contact("missing prepared hammer block"))?
            .workspace.solve(&sites[..n],&a[..n*n],&free[..n],&warm[..n],&mut forces[..n])
            .map_err(Error::Contact)?;
        // Publish the block into NUMERICAL scratch only after it converges.
        // Every other hammer feels these same signed board reactions. Physical
        // motion/history/air still commit at the original complete-sample gate.
        for col in 0..n {
            let i=self.active[first+col];let change=forces[col]-self.force[i];self.force[i]=forces[col];
            for &j in &self.active {self.gap[j]-=self.contact_h[j*nc+i]*change;}
        }
        Ok(())
    }
}

#[cfg(test)]
#[path = "hammer_contact_engine_tests.rs"]
mod tests;
