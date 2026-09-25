//! Parallel impedance loading through the existing lossless port interconnection.
use crate::{PhsError, PortHamiltonian, Storage, interconnect};

impl PortHamiltonian {
    /// Cold-append storage coordinates without moving any original state or
    /// external port. Their structure rows begin at zero; a caller-supplied
    /// dissipative port may evolve them inside the SAME discrete-gradient step.
    /// The wrapper receives the owned original storage and its exact extent.
    /// It must preserve that prefix's law and supply the complete new gradient.
    /// This is a composition seam, not a certificate for the caller's storage.
    ///
    /// Matrices are copied only during construction. Callers own state/work caps;
    /// dimensions and the original skew/PSD structure are re-admitted unchanged.
    pub fn with_appended_storage(self, additional_states: usize,
        wrap: impl FnOnce(Box<dyn Storage>, usize) -> Box<dyn Storage>)
        -> Result<Self, PhsError>
    {
        let error=||PhsError::Dimension{what:"appended storage dimensions"};
        if additional_states==0 {return Err(error());}
        let n=self.n.checked_add(additional_states).ok_or_else(error)?;
        let square=n.checked_mul(n).ok_or_else(error)?;
        let ports=n.checked_mul(self.m).ok_or_else(error)?;
        let mut j=vec![0.0;square];let mut r=vec![0.0;square];let mut g=vec![0.0;ports];
        for i in 0..self.n {
            j[i*n..i*n+self.n].copy_from_slice(&self.j[i*self.n..(i+1)*self.n]);
            r[i*n..i*n+self.n].copy_from_slice(&self.r[i*self.n..(i+1)*self.n]);
        }
        g[..self.g.len()].copy_from_slice(&self.g);
        Self::new(n,self.m,j,r,g,wrap(self.storage,self.n))
    }

    /// Attach an admitted load while retaining every original external drive port.
    ///
    /// `load_port_to_drive[j]` selects the original port whose flow drives load
    /// port j. The reaction is subtracted from that port's external effort:
    /// `u_base = u_external - S^T y_load`, `u_load = S y_base`.
    /// Thus internal port power cancels, including signed cross-port loads.
    /// Original states stay first; all load states follow. The external port
    /// count, ordering, and meaning do not change. Repeated selections represent
    /// multiple loads on one flow and sum their reactions, not duplicate energy.
    ///
    /// This is only composition. Storage, damping and time discretization remain
    /// with their existing owners; a positive-real approximation still requires
    /// separate evidence that it represents the intended physical impedance.
    ///
    /// # Errors
    /// Incomplete/out-of-range port maps, unrepresentable dimensions, or failure
    /// of the existing composite skew/PSD admission. Callers own resource caps.
    pub fn with_parallel_load(self, load: Self, load_port_to_drive: &[usize])
        -> Result<Self, PhsError>
    {
        if load.m == 0 || load_port_to_drive.len() != load.m
            || load_port_to_drive.iter().any(|&p| p >= self.m)
        { return Err(PhsError::BadPortPairing); }
        let error = || PhsError::Dimension { what: "parallel load dimensions" };
        let dim = self.n.checked_add(load.n).ok_or_else(error)?;
        dim.checked_mul(dim).ok_or_else(error)?;
        let external = self.m;
        let expanded = external.checked_add(load.m).ok_or_else(error)?;
        let extent = self.n.checked_mul(expanded).ok_or_else(error)?;
        let mut g = vec![0.0; extent];
        for row in 0..self.n {
            g[row*expanded..row*expanded+external]
                .copy_from_slice(&self.g[row*external..(row+1)*external]);
            for (j,&p) in load_port_to_drive.iter().enumerate() {
                g[row*expanded+external+j] = self.g[row*external+p];
            }
        }
        // Duplicate only PORT columns, never storage or state. Consuming the
        // duplicate ports with interconnect leaves all original drives external.
        let host = Self::new(self.n, expanded, self.j, self.r, g, self.storage)?;
        let pairs: Vec<_> = (0..load.m).map(|j| (external+j,j)).collect();
        interconnect(host, load, &pairs)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{QuadraticStorage, modal_bank_ports, step};
    fn host() -> PortHamiltonian {
        // Two free unit inertias, with arbitrary positions and velocities.
        PortHamiltonian::new(4,2,
            vec![0.,1.,0.,0., -1.,0.,0.,0., 0.,0.,0.,1., 0.,0.,-1.,0.],
            vec![0.;16],vec![0.,0., 1.,0., 0.,0., 0.,1.],
            Box::new(QuadraticStorage::new(vec![0.,0.,0.,0., 0.,1.,0.,0.,
                0.,0.,0.,0., 0.,0.,0.,1.],4).unwrap())).unwrap()
    }
    fn load(z:f64)->PortHamiltonian {
        modal_bank_ports(&[30.],&[z],&[&[3.],&[-2.]]).unwrap()
    }
    #[test]
    fn signed_multiport_load_preserves_external_drives_and_reacts_reciprocally() {
        let s=host().with_parallel_load(load(0.2),&[1,0]).unwrap();
        assert_eq!(s.state_dim(),6);assert_eq!(s.port_dim(),2);
        let (j,_,g)=s.structure();
        assert_eq!(j[1*6+5],2.);assert_eq!(j[3*6+5],-3.);
        assert_eq!(j[5*6+1],-2.);assert_eq!(j[5*6+3],3.);
        assert_eq!(&g[..8],&[0.,0.,1.,0.,0.,0.,0.,1.]);
        assert!(g[8..].iter().all(|x|*x==0.));
        let x=[0.,0.2,0.,-0.1,0.,0.3];
        assert_eq!(s.output(&x),vec![0.2,-0.1]);
        let effort=s.effort(&x);
        let exchange=(0..6).map(|i|effort[i]*(0..6).map(|k|j[i*6+k]*effort[k]).sum::<f64>()).sum::<f64>();
        assert!(exchange.abs()<1e-14);
    }
    #[test]
    fn loaded_drive_work_and_acoustic_storage_share_the_same_step_ledger() {
        for z in [0.,0.25] {
            let s=host().with_parallel_load(load(z),&[0,1]).unwrap();
            let mut x=vec![0.,0.3,0.,-0.2,0.,0.];let initial=s.hamiltonian(&x);
            let (mut work,mut loss)=(0.,0.);
            for k in 0..300 {
                let f=if k<60 {[0.2,-0.1]}else{[0.,0.]};
                let r=step(&s,&x,&f,0.001).unwrap();work+=r.supplied;loss+=r.dissipated;x=r.x;
            }
            assert!((s.hamiltonian(&x)+loss-initial-work).abs()<1e-9);
            assert!(x[4].abs()+x[5].abs()>1e-5);
            if z==0. {assert_eq!(loss,0.);}else{assert!(loss>0.);}
        }
    }
    #[test]
    fn repeated_ports_sum_reactions_and_invalid_maps_refuse() {
        let s=host().with_parallel_load(load(0.2),&[1,1]).unwrap();
        let (j,_,_)=s.structure();assert_eq!(j[3*6+5],-1.);assert_eq!(j[5*6+3],1.);
        for map in [vec![],vec![0],vec![0,2],vec![usize::MAX,0]] {
            assert!(host().with_parallel_load(load(0.2),&map).is_err());
        }
    }
}
