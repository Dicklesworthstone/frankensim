//! Exact zero-pressure traces on the original hierarchical acoustic space.
//!
//! A pressure-release face constrains every nonzero trace basis function,
//! including edge and face modes. Eliminating these original coordinates
//! preserves the symmetric stiffness/mass pencil and its physical inertia;
//! no large admittance, penalty stiffness or artificial boundary modes enter.

use super::{AcousticFeError, HelmholtzFe, HelmholtzModes, MAX_MODES, tri_quadrature};
use fs_modal::{SliceOptions, SliceReport, slice_window};
use fs_sparse::{Coo, Csr};
use std::collections::BTreeSet;

/// Mixed rigid / ideal pressure-release cavity eigenpairs.
///
/// The wrapper retains the boundary condition independently of the shared
/// eigenpair representation. These modes cannot be passed directly into
/// [`super::ReducedHelmholtz`]: its static enrichment vectors are unconstrained
/// and do not preserve a pressure-release trace. The full-space harmonic
/// pressure-release solve is the driven counterpart.
///
/// Ideal `p = 0` is an explicitly prescribed boundary. It does not approximate
/// exterior radiation, an opening's end correction or a finite impedance.
#[derive(Clone, Debug)]
pub struct PressureReleaseModes {
    /// Full original-coordinate shapes, with exactly zero constrained entries.
    /// Residuals and inertia counts belong to the free-coordinate pencil.
    pub eigenpairs: HelmholtzModes,
    /// Canonically sorted, unique original boundary triangles.
    pub pressure_release_faces: Vec<[u32; 3]>,
    /// Original pressure degrees of freedom retained after exact elimination.
    pub free_dofs: usize,
}

/// Shared by modal, harmonic and sparse transient solves on the same geometry.
pub(super) struct PressureReleaseConstraints {
    pub(super) free_dofs: Vec<usize>,
    pub(super) constrained: Vec<bool>,
    pub(super) original_to_free: Vec<usize>,
    pub(super) faces: BTreeSet<[u32; 3]>,
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), AcousticFeError> {
    if cancelled() {
        Err(AcousticFeError::Cancelled)
    } else {
        Ok(())
    }
}

/// What an exact pressure-release elimination does when it leaves no free pressure degree of freedom.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum EmptyFreePressure {
    /// A fluid-only solve has nothing left to solve: refuse.
    Refuse,
    /// The caller carries wall states of its own, so its complete pencil is not empty: admit the empty set. The coupled
    /// pencil refuses only when the pressure and the wall states together are 0.
    AdmitWithWallStates,
}

impl HelmholtzFe<'_> {
    /// The free and the constrained pressure coordinates of the released faces; see [`EmptyFreePressure`] for the case in
    /// which no free coordinate remains.
    pub(super) fn pressure_release_constraints(
        &self,
        faces: &[[u32; 3]],
        empty: EmptyFreePressure,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PressureReleaseConstraints, AcousticFeError> {
        poll(cancelled)?;
        let mut unique = BTreeSet::new();
        let mut constrained = vec![false; self.ndof()];
        let quadrature = tri_quadrature(self.order + 1);
        for &face in faces {
            poll(cancelled)?;
            let mut key = face;
            key.sort_unstable();
            if !unique.insert(key) {
                return Err(AcousticFeError::InvalidInput {
                    what: "pressure-release boundary faces must be unique",
                });
            }
            // The same original trace enumeration used by boundary loads and
            // diffusion. Vertex-only constraints would leave high-order face
            // pressure live and change the physical boundary condition.
            for dof in self.face_rule(face, &quadrature)?.dofs {
                constrained[dof] = true;
            }
        }
        let mut free_dofs = Vec::with_capacity(self.ndof());
        let mut original_to_free = vec![usize::MAX; self.ndof()];
        for (global, &fixed) in constrained.iter().enumerate() {
            poll(cancelled)?;
            if !fixed {
                original_to_free[global] = free_dofs.len();
                free_dofs.push(global);
            }
        }
        if free_dofs.is_empty() && empty == EmptyFreePressure::Refuse {
            return Err(AcousticFeError::InvalidInput {
                what: "pressure-release boundary leaves no free pressure degrees of freedom",
            });
        }
        poll(cancelled)?;
        Ok(PressureReleaseConstraints {
            free_dofs,
            constrained,
            original_to_free,
            faces: unique,
        })
    }

    /// Cavity eigenpairs with exact `p = 0` on selected boundary faces and
    /// rigid natural conditions on all remaining faces.
    ///
    /// The original sparse eigenvalue owner solves the restricted `(K, M)`
    /// pencil and supplies its inertia count and mass-norm residuals. Shapes
    /// are lifted to the original basis with exact zeros on every constrained
    /// vertex, edge and face coordinate. Empty `pressure_release_faces` calls
    /// the existing rigid-wall solver unchanged, including its constant mode.
    ///
    /// Cancellation is checked during restriction and around `slice_window`;
    /// the original eigensolver has no internal cancellation callback.
    ///
    /// # Errors
    /// Invalid, repeated or internal boundary faces; no free pressure DOFs;
    /// cancellation; and the original modal-window/cap refusals.
    pub fn modes_with_pressure_release(
        &self,
        window_hz: (f64, f64),
        pressure_release_faces: &[[u32; 3]],
        options: &SliceOptions,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PressureReleaseModes, AcousticFeError> {
        let constraints = self.pressure_release_constraints(
            pressure_release_faces,
            EmptyFreePressure::Refuse,
            cancelled,
        )?;
        let eigenpairs = if constraints.faces.is_empty() {
            self.modes(window_hz, options)?
        } else {
            let c = self.medium.sound_speed_m_s;
            let to_lambda = |frequency: f64| {
                let k = core::f64::consts::TAU * frequency / c;
                k * k.abs()
            };
            let report = constraints.modes(
                &self.stiffness,
                &self.mass,
                (to_lambda(window_hz.0), to_lambda(window_hz.1)),
                options,
                cancelled,
            )?;
            HelmholtzModes {
                lambda: report.modes.iter().map(|mode| mode.lambda).collect(),
                frequency_hz: report
                    .modes
                    .iter()
                    .map(|mode| c * mode.lambda.max(0.0).sqrt() / core::f64::consts::TAU)
                    .collect(),
                residual: report.modes.iter().map(|mode| mode.residual).collect(),
                shapes: report.modes.into_iter().map(|mode| mode.phi).collect(),
                below_low: report.below_low,
                below_high: report.below_high,
            }
        };
        poll(cancelled)?;
        Ok(PressureReleaseModes {
            eigenpairs,
            free_dofs: constraints.free_dofs.len(),
            pressure_release_faces: constraints.faces.into_iter().collect(),
        })
    }
}

impl PressureReleaseConstraints {
    pub(super) fn restrict(
        &self,
        matrix: &Csr,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Csr, AcousticFeError> {
        if matrix.nrows() != self.constrained.len() || matrix.ncols() != self.constrained.len() {
            return Err(AcousticFeError::InvalidInput {
                what: "pressure-release restriction requires the original square operator",
            });
        }
        let mut result = Coo::new(self.free_dofs.len(), self.free_dofs.len());
        for (local, &global) in self.free_dofs.iter().enumerate() {
            poll(cancelled)?;
            let (columns, values) = matrix.row(global);
            for (&column, &value) in columns.iter().zip(values) {
                let reduced = self.original_to_free[column];
                if reduced != usize::MAX {
                    result.push(local, reduced, value);
                }
            }
        }
        let result = result.assemble();
        poll(cancelled)?;
        Ok(result)
    }

    pub(super) fn modes(
        &self,
        stiffness: &Csr,
        mass: &Csr,
        window: (f64, f64),
        options: &SliceOptions,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SliceReport, AcousticFeError> {
        let stiffness = self.restrict(stiffness, cancelled)?;
        let mass = self.restrict(mass, cancelled)?;
        poll(cancelled)?;
        let mut report = slice_window(&stiffness, &mass, window, options)?;
        poll(cancelled)?;
        if report.modes.len() > MAX_MODES {
            return Err(AcousticFeError::TooManyModes {
                found: report.modes.len(),
                cap: MAX_MODES,
            });
        }
        for mode in &mut report.modes {
            poll(cancelled)?;
            if mode.phi.len() != self.free_dofs.len() {
                return Err(AcousticFeError::InvalidInput {
                    what: "pressure-release eigensolver lost the original free basis",
                });
            }
            let mut shape = vec![0.0; self.constrained.len()];
            for (&global, &value) in self.free_dofs.iter().zip(&mode.phi) {
                shape[global] = value;
            }
            mode.phi = shape;
        }
        Ok(report)
    }
}

#[cfg(test)]
mod tests;
