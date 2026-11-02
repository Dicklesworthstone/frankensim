//! One adaptive receiver bank shared by every original molecular source mode.
use super::*;

/// Absolute propagated quadrature discrepancies for one physical receiver.
/// These are estimates for integrating fixed solved traces. They exclude
/// boundary discretization, material/source error and cancellation-conditioned
/// relative error. They are not acoustic amplitudes or an extra loss channel.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct MolecularReceiverError {
    pub pressure_pa: f64,
    pub particle_velocity_m_s: [f64; 3],
}

/// Actual shared integration work and source-column discrepancy estimates.
#[derive(Debug, Clone)]
pub struct MolecularReceiverQuadrature {
    /// Integrated once per receiver/frequency, independent of source count.
    pub kernel_evaluations: usize,
    /// Source-major, receiver-minor order, matching the harmonic transfer bank.
    /// Values correspond to unit source rates in each declared source's units.
    pub source_fields: Vec<MolecularReceiverError>,
}

impl MolecularBemObservation<'_> {
    /// Select adaptive retained-triangle receivers for all harmonic and FIR
    /// samples of this original source. The boundary solve/formulation, load
    /// projections, gas law and physical receiver positions are unchanged.
    /// Quadrature is prepared once per frequency and reused for every source.
    ///
    /// The additional row-storage and minimum quadrature work are admitted
    /// here, before any boundary solve. prepare_response retains each source's
    /// propagated pressure/velocity estimates and actual integration count.
    /// This choice is part of the operator hash, preventing reuse of centroid
    /// or differently controlled receiver bakes. It does not certify BEM mesh
    /// convergence or the relative accuracy of a nearly cancelling field.
    pub fn with_near_field(
        mut self,
        options: MolecularNearFieldOptions,
    ) -> Result<Self, MolecularObservationError> {
        options
            .admit(
                self.boundary.surface.areas().len(),
                self.boundary.points.len(),
            )
            .map_err(MolecularError::from)?;
        self.near_field_options = Some(options);
        Ok(self)
    }

    pub fn near_field_options(&self) -> Option<MolecularNearFieldOptions> {
        self.near_field_options
    }

    pub(super) fn observation_samples(
        &self,
        solutions: &[MolecularRadiation],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Vec<C64>, Option<MolecularReceiverQuadrature>), MolecularObservationError> {
        if solutions.is_empty() || solutions.len() != self.boundary.sources.len() {
            return Err(invalid(
                "molecular receivers require the complete original source bank",
            ));
        }
        let receivers = self.boundary.points.len();
        let count = solutions
            .len()
            .checked_mul(receivers)
            .ok_or_else(|| invalid("molecular receiver bank extent overflow"))?;
        let entries = count
            .checked_mul(4)
            .ok_or_else(|| invalid("molecular receiver component extent overflow"))?;
        let mut samples = Vec::with_capacity(entries);
        let Some(options) = self.near_field_options else {
            for solution in solutions {
                BemObservation::flatten(
                    solution.observe(
                        self.boundary.surface,
                        &self.boundary.points,
                        self.boundary.minimum_clearance_m,
                        &mut *cancelled,
                    )?,
                    &mut samples,
                );
            }
            return Ok((samples, None));
        };
        let prepared = solutions[0].prepare_near_field(
            self.boundary.surface,
            &self.boundary.points,
            self.boundary.minimum_clearance_m,
            options,
            &mut *cancelled,
        )?;
        let mut quadrature = MolecularReceiverQuadrature {
            kernel_evaluations: prepared.kernel_evaluations(),
            source_fields: Vec::with_capacity(count),
        };
        let momentum = self.source().gas().density * solutions[0].omega_rad_s();
        for solution in solutions {
            let fields = prepared.evaluate_fields(solution, cancelled)?;
            for ((pressure, gradient), (pressure_error, gradient_error)) in fields
                .pressure_pa
                .into_iter()
                .zip(fields.gradient_pa_per_m)
                .zip(
                    fields
                        .quadrature_error_estimate_pa
                        .into_iter()
                        .zip(fields.quadrature_error_estimate_pa_per_m),
                )
            {
                let velocity_error = gradient_error.map(|error| error / momentum);
                if velocity_error.iter().any(|error| !error.is_finite()) {
                    return Err(invalid("molecular receiver velocity discrepancy overflow"));
                }
                samples.push(pressure);
                samples.extend(gradient);
                quadrature.source_fields.push(MolecularReceiverError {
                    pressure_pa: pressure_error,
                    particle_velocity_m_s: velocity_error,
                });
            }
        }
        Ok((samples, Some(quadrature)))
    }
}
