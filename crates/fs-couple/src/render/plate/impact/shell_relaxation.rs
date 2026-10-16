//! Supplied hereditary FLEXURAL response in the actual curved-shell basis.
//! Nonlinear membrane metric and numerical drilling stay in equilibrium storage.
use super::{ImpactError,ImpactSystem,InitialMemory,RelaxationBranch,invalid};
use super::super::BodyPotential;
use fs_material::visco::GeneralizedMaxwell;

/// One shell's additional Maxwell bending spectrum. Ratios refer to the
/// equilibrium DKT stiffness already in its geometric reduction. This does not
/// replace a measured total modal decay or infer bronze properties from a name.
#[derive(Clone,Debug)]
pub struct ShellBendingSpectrum {
    /// Original body index, not a modal/state address.
    pub body: usize,
    /// (excess flexural modulus / equilibrium modulus, relaxation time [s]).
    /// Zero excess contributes no storage. No positive branch is truncated.
    pub branches: Vec<(f64,f64)>,
    /// Declared material validity window [Hz], including instantaneous stiffness.
    /// The zero-frequency translation is excluded from the lower-band check.
    pub band_hz: [f64;2],
}

impl ImpactSystem {
    /// Attach all selected shells in one cold material-memory operation.
    ///
    /// Factor the actual physical bending form K_b=L L^T in its nonzero-energy
    /// subspace, retaining every cross-mode term. Each supplied branch is
    /// H=(sqrt(ratio)*L^T*q-z)^2/2 with R_zz=1/tau, realized by the EXISTING
    /// fs-phs relaxation owner. Only exactly zero rows are absent from material
    /// storage (e.g. rigid translation); singular nonzero bending subspaces
    /// refuse, never use a pseudo-inverse, diagonal approximation or rank cutoff.
    ///
    /// All original mechanical coordinates, forces, felt/contact histories and
    /// acoustic addresses stay fixed. Callers must remove the selected shell's
    /// separate intrinsic modal damping at construction, not its stand/mufflers
    /// or radiation load. Membrane stretching and drilling do not relax here.
    ///
    /// # Errors
    /// Invalid/mismatched bodies or spectra, unresolvable factors, exceeded
    /// material/state budget or instantaneous frequency/relaxation-time guards.
    /// The existing one-attachment/before-first-step rule remains unchanged.
    pub fn with_shell_bending_relaxation(self, spectra: &[ShellBendingSpectrum],
        initial: InitialMemory, maximum_branches: usize) -> Result<Self, ImpactError>
    {
        if spectra.is_empty() {return Ok(self);}
        if spectra.len()>self.mechanical.bodies.len() || maximum_branches>256 {
            return Err(invalid("shell bending memory exceeds body or branch budget"));
        }
        let mut arms=Vec::new();
        for (selection,spec) in spectra.iter().enumerate() {
            if spectra[..selection].iter().any(|s|s.body==spec.body)
                || spec.branches.len()>8 || spec.band_hz.iter().any(|v|!v.is_finite())
                || spec.band_hz[0]<0. || spec.band_hz[1]<=spec.band_hz[0] {
                return Err(invalid("shell bending memory needs distinct bodies, <=8 branches and a finite ordered band"));
            }
            // Reuse material admission; do not implement another Maxwell law.
            let law=GeneralizedMaxwell::new(1.,spec.branches.clone())
                .map_err(|e|ImpactError::Owner(e.to_string()))?;
            let Some(BodyPotential::Shell(shell))=self.mechanical.bodies.get(spec.body) else {
                return Err(invalid("bending memory body must be the original curved shell"));
            };
            if shell.omegas().iter().any(|w|*w>0. && (*w<std::f64::consts::TAU*spec.band_hz[0]
                || *w>std::f64::consts::TAU*spec.band_hz[1])) {
                return Err(invalid("shell eigenfrequency is outside its supplied material band"));
            }
            let positive=law.terms.iter().filter(|t|t.0>0.).count();
            if positive==0 {continue;}
            let n=shell.mode_count();
            // A cold operator projection, never a per-step mesh operation.
            let matrix=shell.bending_stiffness(&vec![1.;shell.facet_count()],2_000_000)
                .map_err(|e|ImpactError::Owner(e.to_string()))?;
            let active:Vec<_>=(0..n).filter(|&i|matrix[i*n..(i+1)*n].iter().any(|v|*v!=0.)).collect();
            if active.is_empty() {continue;}
            if arms.len().checked_add(active.len()*positive).is_none_or(|m|m>maximum_branches) {
                return Err(invalid("complete shell bending spectrum exceeds the material memory budget"));
            }
            let ratio: f64=law.terms.iter().map(|t|t.0).sum();
            let mut tangent=vec![0.;n*n];let zero=vec![0.;n];let mut direction=zero.clone();let mut column=zero.clone();
            for j in 0..n {
                direction[j]=1.;shell.hessian_vector(&zero,&direction,&mut column);direction[j]=0.;
                for i in 0..n {tangent[i*n+j]=column[i]+ratio*matrix[i*n+j];}
            }
            let bound=(0..n).map(|i|tangent[i*n..(i+1)*n].iter().map(|v|v.abs()).sum::<f64>())
                .fold(0.0_f64,f64::max).sqrt();
            if !ratio.is_finite() || tangent.iter().any(|v|!v.is_finite()) || !bound.is_finite()
                || bound*self.config.dt_s>=0.9*std::f64::consts::PI || bound>std::f64::consts::TAU*spec.band_hz[1]
                || law.terms.iter().any(|&(r,t)|r>0. && self.config.dt_s/t>0.25) {
                return Err(invalid("instantaneous shell bending or material relaxation is unresolved by the supplied band/clock"));
            }
            let size=active.len();let mut dense=vec![0.;size*size];
            for (i,&a) in active.iter().enumerate() {for (j,&b) in active.iter().enumerate() {dense[i*size+j]=matrix[a*n+b];}}
            let factor=fs_la::factor::cholesky(&dense,size).map_err(|_|invalid("shell bending form needs a resolved positive subspace; no rank truncation"))?;
            for i in 0..size {for j in 0..size {
                let actual=(0..size).map(|k|factor.l(i,k)*factor.l(j,k)).sum::<f64>();
                let scale=(dense[i*size+i].sqrt()*dense[j*size+j].sqrt()).max(f64::MIN_POSITIVE);
                if !actual.is_finite() || (actual-dense[i*size+j]).abs()>1e-10*scale {
                    return Err(invalid("shell bending factor does not reconstruct the physical energy"));
                }
            }}
            let start=self.mechanical.bodies[..spec.body].iter().map(BodyPotential::count).sum::<usize>();
            for &(ratio,tau) in &law.terms {
                if ratio==0. {continue;}
                for i in 0..size {
                    let mut projection=vec![0.;self.x.len()];
                    for j in i..size {projection[2*(start+active[j])]=factor.l(j,i);}
                    arms.push(RelaxationBranch {projection,stiffness:ratio,relaxation_time_s:tau});
                }
            }
        }
        self.with_relaxation_branches(arms,initial,maximum_branches)
    }
}

#[cfg(test)]
#[path="shell_relaxation_tests.rs"]
mod tests;
