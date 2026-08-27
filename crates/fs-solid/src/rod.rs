//! Geometrically exact Cosserat rods (bead tfz.14): Lie-group nodal
//! state — positions in R³ plus canonical `fs-ga::So3` nodal frames updated
//! MULTIPLICATIVELY on the group (never additive quaternion arithmetic) —
//! with the full strain set: axial/shear
//! `Γ = Rᵀ r′ − e₁` and bending/torsion `κ = Log(Rᵢ⁻¹Rᵢ₊₁)/L₀`,
//! both built from RELATIVE quantities so rigid motions produce
//! exactly zero strain (the objectivity battery checks the energy is
//! invariant, not just small).
//!
//! Statics: total-energy formulation with finite-difference residual
//! (body/right rotational perturbations — `R <- R Exp(delta)`, never
//! component nudging) and FD
//! tangent, dense-LU Newton with load stepping. Fixture-scale by
//! design (≤ a few hundred DOFs); analytic tangents and SE(3)
//! DYNAMICS under fs-time's symplectic integrators are the recorded
//! successor scope.

use crate::SolidError;
use fs_ga::{GaError, So3, Vec3};

/// Diagonal section stiffness of a Cosserat rod.
#[derive(Debug, Clone, Copy)]
pub struct RodSection {
    /// Axial stiffness EA.
    pub ea: f64,
    /// Shear stiffness GA (both transverse directions).
    pub ga: f64,
    /// Torsional stiffness GJ.
    pub gj: f64,
    /// Bending stiffness EI (both bending directions).
    pub ei: f64,
}

/// A discrete rod: reference = straight along +x with uniform segment
/// length; state = nodal positions plus canonical `SO(3)` frames (body→world).
#[derive(Debug, Clone)]
pub struct Rod {
    /// Nodal positions.
    pub positions: Vec<[f64; 3]>,
    /// Validated, canonical nodal frames mapping body vectors into the world.
    pub frames: Vec<So3>,
    /// Reference segment length.
    pub l0: f64,
    /// Section stiffness.
    pub section: RodSection,
}

/// Dead end-loading at the rod tip.
#[derive(Debug, Clone, Copy, Default)]
pub struct TipLoad {
    /// Dead force on the last node (world frame).
    pub force: [f64; 3],
    /// Dead moment on the last node (body frame of the tip).
    pub moment: [f64; 3],
}

fn group_error(error: GaError) -> SolidError {
    SolidError::InternalInvariant {
        what: format!("canonical rod-frame operation refused: {error}"),
    }
}

impl Rod {
    /// A straight reference rod along +x with `segments` segments of
    /// total length `length`.
    #[must_use]
    #[allow(clippy::cast_precision_loss)]
    pub fn straight(length: f64, segments: usize, section: RodSection) -> Rod {
        let l0 = length / segments as f64;
        Rod {
            positions: (0..=segments).map(|i| [i as f64 * l0, 0.0, 0.0]).collect(),
            frames: vec![So3::identity(); segments + 1],
            l0,
            section,
        }
    }

    /// Segment strains: (Γ − e₁-relative, κ), both in the segment
    /// frame.
    /// # Errors
    /// Returns an internal-invariant refusal if a canonical group operation
    /// detects accumulated non-unit drift.
    pub fn strains(&self, seg: usize) -> Result<([f64; 3], [f64; 3]), SolidError> {
        let (a, b) = (seg, seg + 1);
        let relative = self.frames[b]
            .body_minus(self.frames[a])
            .map_err(group_error)?;
        let midpoint = self.frames[a]
            .body_plus(relative.scale(0.5))
            .map_err(group_error)?;
        let dr = Vec3::new(
            (self.positions[b][0] - self.positions[a][0]) / self.l0,
            (self.positions[b][1] - self.positions[a][1]) / self.l0,
            (self.positions[b][2] - self.positions[a][2]) / self.l0,
        );
        // Γ = Rᵀ r′ − e₁ (rotate world tangent into the mid frame).
        let local = midpoint.inverse().rotate(dr).map_err(group_error)?;
        let gamma = [local.x - 1.0, local.y, local.z];
        let angular = relative.angular;
        let kappa = [
            angular.x / self.l0,
            angular.y / self.l0,
            angular.z / self.l0,
        ];
        Ok((gamma, kappa))
    }

    /// Total internal strain energy.
    /// # Errors
    /// Propagates canonical frame-operation refusals from [`Rod::strains`].
    pub fn energy(&self) -> Result<f64, SolidError> {
        let s = &self.section;
        let mut e = 0.0;
        for seg in 0..self.positions.len() - 1 {
            let (g, k) = self.strains(seg)?;
            e += 0.5
                * self.l0
                * (s.ea * g[0] * g[0]
                    + s.ga * (g[1] * g[1] + g[2] * g[2])
                    + s.gj * k[0] * k[0]
                    + s.ei * (k[1] * k[1] + k[2] * k[2]));
        }
        Ok(e)
    }

    /// Free DOFs: everything except node 0 (clamped position + frame);
    /// layout per free node: [dx, dy, dz, θx, θy, θz].
    fn ndof(&self) -> usize {
        6 * (self.positions.len() - 1)
    }

    fn apply_increment(&mut self, delta: &[f64], scale: f64) -> Result<(), SolidError> {
        for node in 1..self.positions.len() {
            let k = 6 * (node - 1);
            for c in 0..3 {
                self.positions[node][c] += scale * delta[k + c];
            }
            let increment = Vec3::new(
                scale * delta[k + 3],
                scale * delta[k + 4],
                scale * delta[k + 5],
            );
            self.frames[node] = self.frames[node]
                .body_plus(fs_ga::So3Tangent::new(increment))
                .map_err(group_error)?;
        }
        Ok(())
    }

    /// Potential Π = E_int − F·r_tip (dead force; the dead tip moment
    /// enters the residual directly — it has no global potential under
    /// multiplicative updates).
    fn potential(&self, load: &TipLoad, factor: f64) -> Result<f64, SolidError> {
        let tip = self.positions[self.positions.len() - 1];
        Ok(self.energy()?
            - factor * (load.force[0] * tip[0] + load.force[1] * tip[1] + load.force[2] * tip[2]))
    }

    /// FD residual in body/right perturbation coordinates: ∂Π/∂dof,
    /// minus the body-frame tip moment on the last node's rotational DOFs.
    fn residual(&self, load: &TipLoad, factor: f64) -> Result<Vec<f64>, SolidError> {
        let n = self.ndof();
        let mut r = vec![0.0f64; n];
        // FD scales: the tangent is FD-of-FD — the residual step must
        // sit well above energy roundoff or the nested difference is
        // noise (measured: a junk Newton direction that creeps).
        let eps = 3e-6;
        let mut probe = self.clone();
        for k in 0..n {
            let mut d = vec![0.0f64; n];
            d[k] = eps;
            probe.clone_from(self);
            probe.apply_increment(&d, 1.0)?;
            let ep = probe.potential(load, factor)?;
            probe.clone_from(self);
            probe.apply_increment(&d, -1.0)?;
            let em = probe.potential(load, factor)?;
            r[k] = (ep - em) / (2.0 * eps);
        }
        let tipk = n - 3;
        r[tipk] -= factor * load.moment[0];
        r[tipk + 1] -= factor * load.moment[1];
        r[tipk + 2] -= factor * load.moment[2];
        Ok(r)
    }

    /// Newton statics under load stepping; returns residual norms per
    /// step (evidence).
    ///
    /// # Errors
    /// [`SolidError::NewtonStalled`] with the history on failure.
    pub fn solve_static(
        &mut self,
        load: &TipLoad,
        steps: usize,
        tol: f64,
    ) -> Result<Vec<f64>, SolidError> {
        let n = self.ndof();
        let mut finals = Vec::new();
        for step in 1..=steps {
            #[allow(clippy::cast_precision_loss)]
            let factor = step as f64 / steps as f64;
            let mut history = Vec::new();
            let mut converged = false;
            for _ in 0..40 {
                let r = self.residual(load, factor)?;
                let rn = r.iter().map(|x| x * x).sum::<f64>().sqrt();
                history.push(rn);
                if rn < tol {
                    converged = true;
                    break;
                }
                // FD tangent (fixture-scale dense).
                let eps = 3e-4;
                let mut kmat = vec![0.0f64; n * n];
                let mut probe = self.clone();
                for col in 0..n {
                    let mut d = vec![0.0f64; n];
                    d[col] = eps;
                    probe.clone_from(self);
                    probe.apply_increment(&d, 1.0)?;
                    let rp = probe.residual(load, factor)?;
                    probe.clone_from(self);
                    probe.apply_increment(&d, -1.0)?;
                    let rm = probe.residual(load, factor)?;
                    for row in 0..n {
                        kmat[row * n + col] = (rp[row] - rm[row]) / (2.0 * eps);
                    }
                }
                let f = fs_la::factor::lu(&kmat, n).map_err(|_| SolidError::SolveFailed {
                    iters: 0,
                    rel_residual: f64::INFINITY,
                })?;
                let mut d: Vec<f64> = r.iter().map(|x| -x).collect();
                f.solve(&mut d);
                // Backtracking on the residual norm (dead moments make
                // Π alone an incomplete merit function).
                let mut alpha = 1.0f64;
                let mut accepted = false;
                for _ in 0..20 {
                    let mut trial = self.clone();
                    trial.apply_increment(&d, alpha)?;
                    let rt = trial.residual(load, factor)?;
                    let rtn = rt.iter().map(|x| x * x).sum::<f64>().sqrt();
                    if rtn < rn {
                        *self = trial;
                        accepted = true;
                        break;
                    }
                    alpha *= 0.5;
                }
                if !accepted {
                    return Err(SolidError::NewtonStalled { history });
                }
            }
            if !converged {
                return Err(SolidError::NewtonStalled { history });
            }
            finals.push(*history.last().expect("nonempty"));
        }
        Ok(finals)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn linear_tip_deflection_matches_beam_theory() {
        let section = RodSection {
            ea: 1e4,
            ga: 1e4,
            gj: 1.0,
            ei: 1.0,
        };
        let mut rod = Rod::straight(1.0, 8, section);
        let p = 0.01; // PL²/EI = 0.01: linear regime
        let load = TipLoad {
            force: [0.0, p, 0.0],
            moment: [0.0; 3],
        };
        let hist = rod.solve_static(&load, 1, 1e-9).expect("linear solve");
        let want = p / (3.0 * section.ei); // PL³/3EI
        let got = rod.positions[8][1];
        assert!(
            (got - want).abs() < 0.02 * want,
            "tip {got} vs beam theory {want}; residual history {hist:?}"
        );
    }
}
