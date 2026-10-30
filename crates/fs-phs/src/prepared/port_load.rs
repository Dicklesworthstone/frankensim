//! Parallel impedance loading through the existing lossless port interconnection.
use crate::{PhsError, PortHamiltonian, Storage, interconnect};

// Constitutive callbacks are immutable. All evolving memory stays in x.
struct SharedStorage(std::rc::Rc<dyn Storage>);
impl Storage for SharedStorage {
    fn shared_clone(&self) -> Option<Box<dyn Storage>> {
        Some(Box::new(Self(self.0.clone())))
    }
    fn hamiltonian(&self, x: &[f64]) -> f64 { self.0.hamiltonian(x) }
    fn gradient(&self, x: &[f64], out: &mut [f64]) { self.0.gradient(x, out); }
    fn quadratic_part(&self) -> Option<crate::linear_exact::QuadraticPart> {
        self.0.quadratic_part()
    }
    fn exact_quadratic_matrix(&self) -> Option<&[f64]> { self.0.exact_quadratic_matrix() }
    fn midpoint_energy_remainder(&self, a: &[f64], b: &[f64], mid: &[f64]) -> Option<f64> {
        self.0.midpoint_energy_remainder(a, b, mid)
    }
}

impl PortHamiltonian {
    /// Expose additional work-conjugate ports without changing the physical law.
    ///
    /// `columns` is row-major state-by-additional-port. Original external ports
    /// remain first, followed by these columns. Every state, structure entry and
    /// constitutive callback is retained unchanged. This cold composition seam
    /// neither infers a force/volume map nor creates new storage or dynamics.
    /// Callers supply the physical units and own subsequent port interconnection.
    pub fn with_additional_ports(self, additional_ports: usize, columns: &[f64])
        -> Result<Self, PhsError>
    {
        let error = || PhsError::Dimension { what: "additional physical port columns" };
        if additional_ports == 0 || self.n.checked_mul(additional_ports) != Some(columns.len())
            || columns.iter().any(|v| !v.is_finite())
        { return Err(error()); }
        let m = self.m.checked_add(additional_ports).ok_or_else(error)?;
        let mut g = vec![0.0; self.n.checked_mul(m).ok_or_else(error)?];
        for i in 0..self.n {
            g[i*m..i*m+self.m].copy_from_slice(&self.g[i*self.m..(i+1)*self.m]);
            g[i*m+self.m..(i+1)*m]
                .copy_from_slice(&columns[i*additional_ports..(i+1)*additional_ports]);
        }
        let mut system = Self::new(self.n, m, self.j, self.r, g, self.storage)?;
        system.resistance_witness = self.resistance_witness;
        Ok(system)
    }

    /// Read-only constitutive law. Physical histories belong to the explicit x.
    #[must_use]
    pub fn storage(&self) -> &dyn Storage { self.storage.as_ref() }

    /// Cold-copy J/R/G while sharing the original immutable storage callbacks.
    ///
    /// This permits transactional preparation of a joint owner without moving
    /// or resetting its source's accepted state. Every callback, including the
    /// stable nonlinear energy remainder, is delegated unchanged. There is no
    /// approximation or inferred quadratic law. Allocations occur here only;
    /// shared callbacks retain the original owner's allocation behavior.
    ///
    /// Storage must obey its existing read-only trial contract. The caller owns
    /// state synchronization and a single accepted clock; this shares no x.
    #[must_use]
    pub fn share_storage(&mut self) -> Self {
        let j = self.j.clone();
        let r = self.r.clone();
        let g = self.g.clone();
        if let Some(storage)=self.storage.shared_clone() {
            return Self { n:self.n, m:self.m, j, r, g, storage,
                resistance_witness:self.resistance_witness.clone() };
        }
        let original = core::mem::replace(&mut self.storage,
            Box::new(crate::QuadraticStorage { n: 0, q: Vec::new() }));
        let shared: std::rc::Rc<dyn Storage> = std::rc::Rc::from(original);
        self.storage = Box::new(SharedStorage(shared.clone()));
        Self { n:self.n, m:self.m, j, r, g, storage:Box::new(SharedStorage(shared)),
            resistance_witness:self.resistance_witness.clone() }
    }

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
        let witness=crate::compose_resistance_witness(&[&self],n,&r)?;
        let mut system=Self::new(n,self.m,j,r,g,wrap(self.storage,self.n))?;
        system.resistance_witness=witness;
        Ok(system)
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
        let extent = self.n.checked_mul(load.m).ok_or(PhsError::Dimension {
            what: "parallel load column dimensions",
        })?;
        let mut columns = vec![0.0; extent];
        for row in 0..self.n {
            for (column, &drive) in load_port_to_drive.iter().enumerate() {
                columns[row*load.m+column] = self.g[row*self.m+drive];
            }
        }
        self.with_parallel_load_columns(load, &columns)
    }

    /// Attach a load at physical flow ports distinct from the external drives.
    ///
    /// `load_columns` is the finite row-major host-state-by-load-port matrix B.
    /// The load receives `u_load = B^T grad H_host`; its reaction on the host is
    /// `-B y_load`. Both maps use the SAME columns, so the internal power is
    /// exactly zero even for signed, nonorthogonal or repeated projections.
    /// For canonical mechanics, entries on momentum rows select velocities;
    /// physical area weights select a cavity's volume flow. They need not be
    /// proportional to an external point-force projection.
    ///
    /// Every original external drive, state and constitutive law stays in its
    /// original order. All load states follow, and no load port remains external.
    /// The existing interconnect owner assembles and admits the full skew
    /// coupling. This supplies no time integrator, reduction or impedance fit.
    /// Callers own state/work capacities and the physical units of each column.
    ///
    /// # Errors
    /// An incomplete/nonfinite column matrix, zero load ports, unrepresentable
    /// dimensions, or the existing composite skew/PSD admission refusal.
    pub fn with_parallel_load_columns(self, load: Self, load_columns: &[f64])
        -> Result<Self, PhsError>
    {
        let error = || PhsError::Dimension { what: "parallel load dimensions" };
        if load.m == 0 || self.n.checked_mul(load.m) != Some(load_columns.len())
            || load_columns.iter().any(|x| !x.is_finite())
        { return Err(PhsError::BadPortPairing); }
        let dim = self.n.checked_add(load.n).ok_or_else(error)?;
        dim.checked_mul(dim).ok_or_else(error)?;
        let external = self.m;
        let expanded = external.checked_add(load.m).ok_or_else(error)?;
        let extent = self.n.checked_mul(expanded).ok_or_else(error)?;
        let mut g = vec![0.0; extent];
        for row in 0..self.n {
            g[row*expanded..row*expanded+external]
                .copy_from_slice(&self.g[row*external..(row+1)*external]);
            g[row*expanded+external..(row+1)*expanded]
                .copy_from_slice(&load_columns[row*load.m..(row+1)*load.m]);
        }
        // Extra columns expose physical ports, not new states or storage.
        // Consuming them leaves every original external drive untouched.
        let host = Self::new(self.n, expanded, self.j, self.r, g, self.storage)?;
        let pairs: Vec<_> = (0..load.m).map(|j| (external+j,j)).collect();
        interconnect(host, load, &pairs)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{QuadraticStorage, modal_bank_ports, step};
    #[test]
    fn additional_physical_ports_keep_original_memory_and_reciprocal_work() {
        let make=||crate::mass_spring_damper(2.0,7.0,0.3).unwrap()
            .with_relaxation_branches(vec![crate::RelaxationBranch {
                projection:vec![1.0,0.0],stiffness:3.0,relaxation_time_s:0.2,
            }]).unwrap();
        let x=[0.25,-0.1,0.03];
        let source=make();
        let extended=make().with_additional_ports(2,&[0.0,0.0,-0.5,2.0,0.0,0.0]).unwrap();
        assert_eq!(extended.state_dim(),source.state_dim());
        assert_eq!(extended.port_dim(),3);
        assert_eq!(extended.hamiltonian(&x),source.hamiltonian(&x));
        assert_eq!(extended.effort(&x),source.effort(&x));
        let velocity=x[1]/2.0;
        assert_eq!(extended.output(&x),vec![velocity,-0.5*velocity,2.0*velocity]);
        // Three physical forces at different work-conjugate leverages give
        // exactly the same total force 1 - 0.5*2 + 2*3 = 6 on the original mass.
        let a=step(&extended,&x,&[1.0,2.0,3.0],0.001).unwrap();
        let b=step(&source,&x,&[6.0],0.001).unwrap();
        assert_eq!(a.x,b.x);
        assert_eq!(a.delta_h,b.delta_h);
        assert_eq!(a.dissipated,b.dissipated);
        assert!((a.supplied-b.supplied).abs()<1e-15);
        assert!(make().with_additional_ports(0,&[]).is_err());
        assert!(make().with_additional_ports(1,&[0.0,1.0]).is_err());
        assert!(make().with_additional_ports(1,&[0.0,f64::NAN,0.0]).is_err());
    }
    #[test]
    fn repeated_cold_sharing_preserves_the_complete_storage_without_wrapper_growth() {
        let mut source=crate::mass_spring_damper(2.0,7.0,0.3).unwrap()
            .with_relaxation_branches(vec![crate::RelaxationBranch {
                projection:vec![1.0,0.0],stiffness:3.0,relaxation_time_s:0.2,
            }]).unwrap();
        let x=[0.25,-0.1,0.03];let next=[0.24,-0.09,0.031];
        let mid:Vec<_>=x.iter().zip(next).map(|(a,b)|f64::midpoint(*a,b)).collect();
        let energy=source.hamiltonian(&x);let effort=source.effort(&x);
        let remainder=source.storage().midpoint_energy_remainder(&x,&next,&mid);
        let quadratic=source.storage().quadratic_part().unwrap();
        let first=source.share_storage();
        let identity=core::ptr::from_ref(source.storage()).cast::<()>().addr();
        for _ in 0..1024 {
            let copy=source.share_storage();
            assert_eq!(core::ptr::from_ref(source.storage()).cast::<()>().addr(),identity);
            assert_eq!(copy.structure(),source.structure());
            assert_eq!(copy.hamiltonian(&x),energy);assert_eq!(copy.effort(&x),effort);
            assert_eq!(copy.storage().midpoint_energy_remainder(&x,&next,&mid),remainder);
            let q=copy.storage().quadratic_part().unwrap();
            assert_eq!(q.exact,quadratic.exact);assert_eq!(q.matrix,quadratic.matrix);
        }
        let a=step(&source,&x,&[0.1],0.001).unwrap();
        let b=step(&first,&x,&[0.1],0.001).unwrap();
        assert_eq!(a.x,b.x);assert_eq!(a.dissipated,b.dissipated);assert_eq!(a.supplied,b.supplied);
    }
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
