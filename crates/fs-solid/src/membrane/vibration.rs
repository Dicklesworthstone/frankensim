//! Small vibrations about the SAME source-bound, dead-load equilibrium.
//!
//! The membrane owner supplies its exact material + geometric tangent and
//! conserved reference mass; `fs-modal` supplies the generalized eigensolve.
//! This is a bounded dense lowering, not a new eigensolver or an unstressed
//! plate approximation. All free coordinates are retained. A zero/negative or
//! unresolved eigenvalue refuses rather than being clipped into an audible mode.
//! Pressure/volume ports below are linearized about the accepted geometry;
//! they do not add follower-pressure prestress or a radiation/damping model.

use super::{MembraneError, MembranePatch, MembraneSolution, MembraneSolveSettings,
    finite_vectors, poll};
use fs_modal::{ModePair, ModalError, eigh_gen_dense};

/// Hard limit on the existing non-interruptible dense eigensolver leaf.
/// Cancellation is checked before/after that leaf and at each tangent column.
pub const MAX_DENSE_FREE_DOFS: usize = 128;

/// Explicit work and eigenpair acceptance bounds.
#[derive(Clone, Copy, Debug)]
pub struct VibrationSettings {
    /// Caller allocation/work cap, additionally limited by MAX_DENSE_FREE_DOFS.
    pub max_free_dofs: usize,
    /// Maximum M-inverse residual / absolute eigenvalue, strictly in (0,1).
    pub relative_residual_tolerance: f64,
}

/// Physical equilibrium, spectrum or projection refusal.
#[derive(Debug)]
pub enum VibrationError {
    /// Original material/geometry/nonlinear failure, preserved without fallback.
    Membrane(MembraneError),
    /// Original generalized eigensolver failure.
    Modal(ModalError),
    /// A zero, negative, nonfinite or insufficiently resolved structural mode.
    UnresolvedMode { index: usize, lambda: f64, residual: f64 },
}
impl From<MembraneError> for VibrationError {
    fn from(error: MembraneError) -> Self { Self::Membrane(error) }
}
impl core::fmt::Display for VibrationError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result { write!(f,"membrane vibration: {self:?}") }
}
impl std::error::Error for VibrationError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self { Self::Membrane(e) => Some(e), Self::Modal(e) => Some(e), _ => None }
    }
}

/// Immutable equilibrium, complete modal basis and its original material owner.
/// q [m sqrt(kg)], qdot [m sqrt(kg)/s], phi [1/sqrt(kg)]. No guessed damping,
/// target frequencies, independent mass input or mode-selection heuristic.
#[derive(Debug)]
pub struct MembraneVibration<'a> {
    patch: &'a MembranePatch,
    equilibrium: MembraneSolution,
    free_dofs: Vec<usize>,
    modes: Vec<ModePair>,
    pressure_load_per_pa: Vec<[f64;3]>,
    volume_coupling: Vec<f64>,
}

impl MembranePatch {
    /// Solve the original equilibrium and lower its exact tangent using the
    /// SAME prescribed DOFs and material/reference mass. Keeping these two
    /// operations together prevents a foreign/stale solution or a changed
    /// constraint list from silently becoming a supposedly prestressed pencil.
    /// The full positive discrete spectrum is checked, not just an audible
    /// window which could conceal an unstable mode. No continuum stability,
    /// mesh-convergence or finite-amplitude vibration claim is made.
    pub fn solve_vibration(
        &self, dead_load_n: &[[f64;3]], prescribed: &[(usize,f64)],
        initial_displacement_m: &[[f64;3]], solve: MembraneSolveSettings,
        settings: VibrationSettings, mut cancelled: impl FnMut() -> bool,
    ) -> Result<MembraneVibration<'_>, VibrationError> {
        poll(&mut cancelled)?;
        let nfull = 3*self.mesh().nodes.len();
        if !settings.relative_residual_tolerance.is_finite()
            || !(0.0..1.0).contains(&settings.relative_residual_tolerance)
            || settings.relative_residual_tolerance == 0.0 {
            return Err(MembraneError::Invalid("relative modal residual tolerance must be in (0,1)").into());
        }
        let mut fixed = vec![false;nfull];
        for &(dof,value) in prescribed {
            poll(&mut cancelled)?;
            if dof>=nfull || fixed[dof] || !value.is_finite() {
                return Err(MembraneError::Invalid("unique finite prescribed displacements required").into());
            }
            fixed[dof]=true;
        }
        let free_dofs:Vec<_>=(0..nfull).filter(|&i|!fixed[i]).collect();
        let n=free_dofs.len();
        if n==0 || n>settings.max_free_dofs || n>MAX_DENSE_FREE_DOFS {
            return Err(MembraneError::Budget("positive bounded free-coordinate count for dense vibration").into());
        }
        let equilibrium=self.solve(dead_load_n,prescribed,initial_displacement_m,solve,&mut cancelled)?;
        // A converged kink or source-domain endpoint need not admit a unique
        // two-sided vibration tangent. Keep the original material refusal.
        self.admit_linearization(&equilibrium.displacement_m,&mut cancelled)?;
        let mut k=vec![0.0;n*n];
        let mut m=vec![0.0;n*n];
        let mut direction=vec![[0.0;3];self.mesh().nodes.len()];
        for (col,&dof) in free_dofs.iter().enumerate() {
            poll(&mut cancelled)?;
            m[col*n+col]=self.nodal_mass_kg()[dof/3];
            direction[dof/3][dof%3]=1.0;
            let action=self.tangent_apply(&equilibrium.displacement_m,&direction,&mut cancelled)?;
            direction[dof/3][dof%3]=0.0;
            for (row,&r) in free_dofs.iter().enumerate() { k[row*n+col]=action[r/3][r%3]; }
        }
        // Check the exact assembled action BEFORE the shared dense owner does
        // its roundoff symmetrization. Never symmetrize a nonsymmetric physics law.
        let scale=k.iter().fold(0.0_f64,|a,b|a.max(b.abs()));
        let allowance=64.0*f64::EPSILON*n as f64*scale;
        for i in 0..n { for j in 0..i {
            if (k[i*n+j]-k[j*n+i]).abs()>allowance {
                return Err(MembraneError::Invalid("membrane tangent is not symmetric to roundoff").into());
            }
        } }
        poll(&mut cancelled)?;
        let modes=eigh_gen_dense(&k,&m,n).map_err(VibrationError::Modal)?;
        poll(&mut cancelled)?;
        if modes.len()!=n {
            return Err(MembraneError::Invalid("dense modal owner returned an incomplete basis").into());
        }
        for (index,mode) in modes.iter().enumerate() {
            if mode.phi.len()!=n || mode.phi.iter().any(|v|!v.is_finite())
                || !mode.lambda.is_finite() || !mode.residual.is_finite() || mode.residual<0.0
                || !mode.interval.0.is_finite() || !mode.interval.1.is_finite()
                || mode.interval.0<=0.0
                || mode.residual>settings.relative_residual_tolerance*mode.lambda {
                return Err(VibrationError::UnresolvedMode {index,lambda:mode.lambda,residual:mode.residual});
            }
        }
        // Integral N_a n dA = oriented area vector / 3, on the CURRENT
        // equilibrium triangle. Relative edges avoid global-translation loss.
        let mut pressure_load_per_pa=vec![[0.0;3];self.mesh().nodes.len()];
        for triangle in &self.mesh().elems {
            poll(&mut cancelled)?;
            let edge=|node:usize| -> [f64;3] {
                std::array::from_fn(|axis| {
                    let reference=if axis<2 {self.mesh().nodes[node][axis]-self.mesh().nodes[triangle[0]][axis]} else {0.0};
                    reference+(equilibrium.displacement_m[node][axis]-equilibrium.displacement_m[triangle[0]][axis])
                })
            };
            let a=edge(triangle[1]);let b=edge(triangle[2]);
            let cross=[a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]];
            for &node in triangle { for axis in 0..3 {pressure_load_per_pa[node][axis]+=cross[axis]/6.0;} }
        }
        if !finite_vectors(&pressure_load_per_pa,self.mesh().nodes.len()) {
            return Err(MembraneError::Invalid("unrepresentable equilibrium pressure port").into());
        }
        let mut result=MembraneVibration {patch:self,equilibrium,free_dofs,modes,pressure_load_per_pa,volume_coupling:vec![]};
        result.volume_coupling=result.project_nodal_force(&result.pressure_load_per_pa)?;
        poll(&mut cancelled)?;
        Ok(result)
    }
}

impl MembraneVibration<'_> {
    /// Original immutable geometry, grain, source bundles and reference masses.
    pub fn patch(&self) -> &MembranePatch {self.patch}
    /// Original accepted nonlinear equilibrium, including support reactions.
    pub fn equilibrium(&self) -> &MembraneSolution {&self.equilibrium}
    /// Reduced-to-full DOF map, with prescribed perturbations identically zero.
    pub fn free_dofs(&self) -> &[usize] {&self.free_dofs}
    /// ALL mass-normalized positive modes from the existing modal owner.
    pub fn modes(&self) -> &[ModePair] {&self.modes}
    /// Per-node force / uniform signed pressure jump [m2], on the accepted
    /// geometry. Positive pressure acts along the reference-oriented normal.
    pub fn pressure_load_per_pa(&self) -> &[[f64;3]] {&self.pressure_load_per_pa}
    /// B [m2/sqrt(kg)]: Q=B.qdot, f_modal=B*p. This is an oriented swept-volume
    /// port, not by itself a sealed volume, a radiation model or fluid loading.
    pub fn volume_coupling(&self) -> &[f64] {&self.volume_coupling}
    /// Project nodal dead forces [N] to generalized forces [N/sqrt(kg)].
    /// The same phi reconstructs velocities, so virtual work is preserved.
    pub fn project_nodal_force(&self, force_n:&[[f64;3]]) -> Result<Vec<f64>,MembraneError> {
        if !finite_vectors(force_n,self.patch.mesh().nodes.len()) {
            return Err(MembraneError::Invalid("one finite nodal force per reference node required"));
        }
        let result:Vec<_>=self.modes.iter().map(|mode|mode.phi.iter().zip(&self.free_dofs)
            .map(|(p,&dof)|p*force_n[dof/3][dof%3]).sum::<f64>()).collect();
        if result.iter().any(|x|!x.is_finite()) {return Err(MembraneError::Invalid("unrepresentable modal force"));}
        Ok(result)
    }
    /// Reconstruct a small displacement increment [m] (or velocity [m/s] when
    /// supplied generalized velocities). This does not advance or mutate the
    /// accepted equilibrium or validate a finite-amplitude nonlinear trajectory.
    pub fn reconstruct(&self, coordinates:&[f64]) -> Result<Vec<[f64;3]>,MembraneError> {
        if coordinates.len()!=self.modes.len() || coordinates.iter().any(|x|!x.is_finite()) {
            return Err(MembraneError::Invalid("one finite value per modal coordinate required"));
        }
        let mut result=vec![[0.0;3];self.patch.mesh().nodes.len()];
        for (mode,q) in self.modes.iter().zip(coordinates) {
            for (p,&dof) in mode.phi.iter().zip(&self.free_dofs) {result[dof/3][dof%3]+=p*q;}
        }
        if !finite_vectors(&result,self.patch.mesh().nodes.len()) {return Err(MembraneError::Invalid("unrepresentable modal reconstruction"));}
        Ok(result)
    }
}
