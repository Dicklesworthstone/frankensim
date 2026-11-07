//! Optional accuracy control of the native force/energy endpoint iteration.

use super::{
    COUPLING_PASSES, NsfError, NsfGrid3, Q, THETA0, body_force,
    material::NsfCellMaterial,
    open_boundary::{BoundaryFlux, OpenFaces},
};
use crate::lattice3::Lattice3;

/// Explicit endpoint tolerances; the default grid retains its fixed three passes.
///
/// All four tolerances are finite, positive and dimensionless. Reference scales
/// are the ORIGINAL gas reference, never the magnitude of a vanishing force,
/// velocity or the stored energy datum. With `c_ref = sqrt(R T_ref)` and
/// `cv_ref = R / (gamma_ref - 1)`, the physical scales per unit mass are
/// `c_ref²/dx` for force, `c_ref` for velocity and `cv_ref T_ref` for both energy
/// defects. Their lattice forms are `rho*THETA0`, `sqrt(THETA0)` and
/// `rho*THETA0/(gamma_ref-1)`, respectively.
///
/// This controls the nonlinear endpoint defect. It does not estimate spatial
/// truncation error, change the Mach/material/stability envelope, or qualify
/// additional flow regimes. Source caloric-inversion tolerances can impose a
/// residual floor; exhausting the budget refuses rather than relaxing a limit.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NsfCouplingControl {
    /// Total outer force/velocity/material endpoint evaluations; at least four.
    /// The first three are the unchanged original correctors. Every later
    /// pass probes that candidate and, only if needed and budget remains,
    /// performs another original corrector. Every probe adds one actual-state
    /// stress/energy-flux evaluation. No second collision/stream occurs.
    pub maximum_endpoint_evaluations: usize,
    /// Maximum componentwise force defect divided by `rho*THETA0`.
    pub force_relative_tolerance: f64,
    /// Maximum physical velocity defect divided by `sqrt(THETA0)`.
    pub velocity_relative_tolerance: f64,
    /// Caloric closure defect divided by the fixed reference thermal energy.
    pub caloric_relative_tolerance: f64,
    /// Local conservative energy-equation defect on that same reference scale.
    pub energy_relative_tolerance: f64,
}

impl NsfCouplingControl {
    pub fn validate(self) -> Result<(), NsfError> {
        if self.maximum_endpoint_evaluations <= COUPLING_PASSES
            || !self
                .tolerances()
                .into_iter()
                .all(|value| value.is_finite() && value > 0.0)
        {
            return Err(NsfError::Input(
                "NSF endpoint convergence needs at least four endpoint evaluations and finite positive tolerances",
            ));
        }
        Ok(())
    }

    fn tolerances(self) -> [f64; 4] {
        [
            self.force_relative_tolerance,
            self.velocity_relative_tolerance,
            self.caloric_relative_tolerance,
            self.energy_relative_tolerance,
        ]
    }

    pub(super) fn accepts(self, report: NsfCouplingReport) -> bool {
        report
            .residuals()
            .into_iter()
            .zip(self.tolerances())
            .all(|(residual, tolerance)| residual <= tolerance)
    }
}

/// Endpoint equation with the largest residual-to-requested-tolerance ratio.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NsfCouplingComponent {
    Force,
    Velocity,
    Caloric,
    Energy,
}

/// Measured residuals of the candidate actually accepted by the native step.
/// On [`NsfError::CouplingNotConverged`], these instead describe its unaccepted
/// final candidate. All norms are infinity norms over fluid cells/components.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NsfCouplingReport {
    /// Actual native conservative correctors completed for this candidate.
    pub corrector_passes: usize,
    /// Outer force/velocity/material evaluations, including the final probe.
    /// Exactly `corrector_passes + 1`, on both success and budget refusal.
    pub endpoint_evaluations: usize,
    /// Actual-state stress/energy-flux probes. Total native flux evaluations
    /// are `corrector_passes + residual_probes`; none occur in fixed-pass mode.
    pub residual_probes: usize,
    /// Stored force versus the thermodynamic/source force at the candidate.
    pub force_relative_residual: f64,
    /// Stored physical velocity versus the newly evaluated half-force velocity.
    pub velocity_relative_residual: f64,
    /// Stored E versus source internal energy plus ACTUAL candidate kinetic energy.
    pub caloric_relative_residual: f64,
    /// Stored E versus the native energy update re-evaluated at that same E,
    /// temperature, physical velocity, Guo stress and source/boundary data.
    pub energy_relative_residual: f64,
    /// Equation with the largest requested-tolerance ratio (first on ties).
    pub worst_component: NsfCouplingComponent,
    /// Largest residual divided by its explicitly requested tolerance.
    /// May be infinity for a representably tiny tolerance; acceptance compares
    /// residuals directly, without rounding a tolerance ratio down to one.
    pub maximum_tolerance_ratio: f64,
}

impl NsfCouplingReport {
    fn residuals(self) -> [f64; 4] {
        [
            self.force_relative_residual,
            self.velocity_relative_residual,
            self.caloric_relative_residual,
            self.energy_relative_residual,
        ]
    }

    fn finish(&mut self, control: NsfCouplingControl) {
        let components = [
            NsfCouplingComponent::Force,
            NsfCouplingComponent::Velocity,
            NsfCouplingComponent::Caloric,
            NsfCouplingComponent::Energy,
        ];
        for ((residual, tolerance), component) in self
            .residuals()
            .into_iter()
            .zip(control.tolerances())
            .zip(components)
        {
            let ratio = residual / tolerance;
            if ratio > self.maximum_tolerance_ratio {
                self.maximum_tolerance_ratio = ratio;
                self.worst_component = component;
            }
        }
    }
}

// The original predictor tuple, named so probes borrow the retained histories
// instead of overwriting them. These are the SAME native work vectors.
pub(super) struct Endpoint {
    pub(super) force: Vec<[f64; 3]>,
    pub(super) velocity: Vec<[f64; 3]>,
    pub(super) energy: Vec<f64>,
    pub(super) divergence: Vec<f64>,
    pub(super) wall_flux: f64,
    pub(super) wall_heat_flux: Vec<f64>,
    pub(super) boundary_flux: BoundaryFlux,
    pub(super) faces: Option<OpenFaces>,
    pub(super) material: Vec<Option<NsfCellMaterial>>,
}

impl<L: Lattice3<27>> NsfGrid3<L> {
    /// Select bounded endpoint accuracy control on this same native solver.
    /// Every accepted step reports all four actual endpoint residuals through
    /// [`Self::last_coupling_report`]. No tolerance or pass count is defaulted.
    pub fn with_coupling_control(mut self, control: NsfCouplingControl) -> Result<Self, NsfError> {
        self.set_coupling_control(Some(control))?;
        Ok(self)
    }

    /// Change the accuracy control for subsequent steps; `None` selects the
    /// original three-pass arithmetic. Changing it does not rebind physical
    /// fields or restart accepted flux histories. The last completed step's
    /// report remains available until another step succeeds.
    /// Invalid controls leave both the selection and the accepted state intact.
    pub fn set_coupling_control(
        &mut self,
        control: Option<NsfCouplingControl>,
    ) -> Result<(), NsfError> {
        if let Some(control) = control {
            control.validate()?;
        }
        self.coupling_control = control;
        Ok(())
    }

    /// Selected endpoint controls, or `None` for the original three passes.
    #[must_use]
    pub const fn coupling_control(&self) -> Option<NsfCouplingControl> {
        self.coupling_control
    }

    /// Last successful step's actual endpoint residuals, or `None` before an
    /// accuracy-controlled step or after a successful fixed-pass step.
    /// Refusal/cancellation retains the prior successful report.
    #[must_use]
    pub const fn last_coupling_report(&self) -> Option<NsfCouplingReport> {
        self.last_coupling_report
    }

    // Borrow the already evaluated native predictor. Only the actual
    // candidate stress/flux must be recomputed, using its HELD force and
    // physical velocity. Source callbacks are never invoked by this probe.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn probe_coupling_endpoint(
        &self,
        control: NsfCouplingControl,
        passes: usize,
        candidate: &Endpoint,
        theta: &[f64],
        populations: &[[f64; Q]],
        rho: &[f64],
        momentum: &[[f64; 3]],
        predictor_force: &[[f64; 3]],
        predictor_velocity: &[[f64; 3]],
        predictor_faces: Option<&OpenFaces>,
        material: &[Option<NsfCellMaterial>],
        body_acceleration: &[body_force::BodyAcceleration],
        source_next: &[(f64, [f64; 3], f64)],
        expected_energy: impl Fn(usize, f64, [f64; 3]) -> f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<NsfCouplingReport, NsfError> {
        body_force::checkpoint(cancelled)?;
        // The SGS/enthalpy-contact predictor internally corrects velocity.
        // Its probe must instead resolve faces at the actual retained velocity.
        // Other predictors already used this same candidate as u_guess.
        let contact = if self.coupled_sgs_contact() {
            let faces = self.resolve_open_faces(
                self.open.as_ref(),
                self.characteristic_pressure_offset_pa,
                self.characteristic_inflow_temperature_k,
                self.characteristic_impedance_pa_s_m,
                self.characteristic_patch,
                rho,
                &candidate.velocity,
                theta,
                material,
            )?;
            let mut force = self.force_field_with_body_acceleration(
                rho,
                theta,
                faces.as_ref(),
                body_acceleration,
            );
            if self.sources.is_some() {
                for i in 0..self.n {
                    if !self.solid[i] {
                        for a in 0..3 {
                            force[i][a] += source_next[i].1[a];
                        }
                    }
                }
            }
            let velocity = self.velocity_field(rho, momentum, &force);
            Some((faces, force, velocity))
        } else {
            None
        };
        let (faces, required_force, required_velocity) = match &contact {
            Some((faces, force, velocity)) => {
                (faces.as_ref(), force.as_slice(), velocity.as_slice())
            }
            None => (predictor_faces, predictor_force, predictor_velocity),
        };
        body_force::checkpoint(cancelled)?;
        let transport = self.sgs_field(rho, &candidate.velocity, theta, material, faces)?;
        let stress = self.stress_field(
            populations,
            &candidate.force,
            rho,
            theta,
            material,
            &transport,
        );
        let (divergence, _, _) = self.flux_divergence(
            &candidate.energy,
            rho,
            &candidate.velocity,
            theta,
            &stress,
            &transport,
            faces,
        )?;
        body_force::checkpoint(cancelled)?;
        let mut report = NsfCouplingReport {
            corrector_passes: passes - 1,
            endpoint_evaluations: passes,
            residual_probes: passes - COUPLING_PASSES,
            force_relative_residual: 0.0,
            velocity_relative_residual: 0.0,
            caloric_relative_residual: 0.0,
            energy_relative_residual: 0.0,
            worst_component: NsfCouplingComponent::Force,
            maximum_tolerance_ratio: 0.0,
        };
        let reference_cv = 1.0 / (self.law.gamma - 1.0);
        let velocity_scale = THETA0.sqrt();
        for i in 0..self.n {
            if i % 64 == 0 {
                body_force::checkpoint(cancelled)?;
            }
            if self.solid[i] {
                continue;
            }
            let force_scale = rho[i] * THETA0;
            let energy_scale = force_scale * reference_cv;
            if !(force_scale.is_finite()
                && force_scale > 0.0
                && energy_scale.is_finite()
                && energy_scale > 0.0)
            {
                return Err(NsfError::Input(
                    "NSF coupling reference scale is not representable",
                ));
            }
            for a in 0..3 {
                report.force_relative_residual =
                    report.force_relative_residual.max(relative_difference(
                        candidate.force[i][a],
                        required_force[i][a],
                        force_scale,
                        i,
                    )?);
                report.velocity_relative_residual =
                    report.velocity_relative_residual.max(relative_difference(
                        candidate.velocity[i][a],
                        required_velocity[i][a],
                        velocity_scale,
                        i,
                    )?);
            }
            let u = candidate.velocity[i];
            let internal = material
                .get(i)
                .and_then(Option::as_ref)
                .map_or(reference_cv * theta[i], |value| {
                    value.sensible_energy_lattice
                });
            let caloric = rho[i] * (internal + 0.5 * (u[0] * u[0] + u[1] * u[1] + u[2] * u[2]));
            report.caloric_relative_residual = report.caloric_relative_residual.max(
                relative_difference(candidate.energy[i], caloric, energy_scale, i)?,
            );
            let energy = expected_energy(i, divergence[i], u);
            report.energy_relative_residual = report.energy_relative_residual.max(
                relative_difference(candidate.energy[i], energy, energy_scale, i)?,
            );
        }
        report.finish(control);
        Ok(report)
    }
}

fn relative_difference(a: f64, b: f64, scale: f64, cell: usize) -> Result<f64, NsfError> {
    if !(a.is_finite() && b.is_finite()) {
        return Err(NsfError::NonFinite { cell });
    }
    let difference = a - b;
    let residual = if difference.is_finite() {
        difference.abs() / scale
    } else {
        // Opposite finite values may overflow their unscaled difference.
        (a / scale - b / scale).abs()
    };
    if residual.is_finite() {
        Ok(residual)
    } else {
        Err(NsfError::NonFinite { cell })
    }
}

#[cfg(test)]
mod tests;
