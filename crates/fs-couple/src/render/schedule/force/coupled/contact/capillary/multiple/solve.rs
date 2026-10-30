//! Joint original-source residual, analytic Newton tangent, and error propagation.
use super::*;
use fs_contact::capillary_bridge::contact_line::ContactLineError;

impl MultipleCircularWettingModalSystem<'_> {
    pub(super) fn gap_tolerance(&self, i: usize, gap: f64) -> Result<f64, Error> {
        Ok(finite(
            self.config.gap_absolute_tolerance_m
                + self.config.gap_relative_tolerance
                    * self.points[i]
                        .gap
                        .abs()
                        .max(gap.abs())
                        .max(self.free_gaps[i].abs()),
        )?)
    }
    pub(super) fn force_tolerance(&self, scale: f64) -> Result<f64, Error> {
        Ok(finite(
            self.config.force_absolute_tolerance_n + self.config.force_relative_tolerance * scale,
        )?)
    }

    fn evaluate(
        points: &[Point<'_>],
        gaps: &[f64],
        steps: &mut [Option<AsymmetricCircularContactLineStep>],
        dt: f64,
        maximum_delta: f64,
        evaluations: &mut usize,
        ceiling: usize,
        gate: &CancelGate,
    ) -> Result<(), Error> {
        for (i, point) in points.iter().enumerate() {
            poll(Some(gate))?;
            if !gaps[i].is_finite() || (gaps[i] - point.gap).abs() > maximum_delta {
                return Err(Error::Bridge {
                    index: i,
                    error: ContactLineError::Domain("joint mechanical gap increment").into(),
                });
            }
            if *evaluations == ceiling {
                return Err(Error::Budget);
            }
            *evaluations += 1;
            steps[i] = Some(
                point
                    .source
                    .advance(point.state, gaps[i], dt, gate)
                    .map_err(|e| Error::Bridge {
                        index: i,
                        error: e.into(),
                    })?,
            );
        }
        Ok(())
    }

    fn residual_norm(
        &self,
        gaps: &[f64],
        steps: &[Option<AsymmetricCircularContactLineStep>],
    ) -> Result<f64, Error> {
        let n = self.points.len();
        let mut norm = 0.0_f64;
        for i in 0..n {
            let mut r = finite(gaps[i] - self.free_gaps[i])?;
            for (j, step) in steps.iter().enumerate() {
                r = finite(
                    r - self.compliance[i * n + j]
                        * step
                            .as_ref()
                            .expect("complete original probe")
                            .mean_normal_force_n,
                )?;
            }
            norm = norm.max(finite(r.abs() / self.gap_tolerance(i, gaps[i])?)?);
        }
        Ok(norm)
    }

    fn factor(&mut self, gate: &CancelGate) -> Result<(), Error> {
        let n = self.points.len();
        for i in 0..n {
            self.residual[i] = finite(self.free_gaps[i] - self.gaps[i])?;
            for j in 0..n {
                let step = self.steps[j].as_ref().expect("complete original probe");
                self.residual[i] = finite(
                    self.residual[i] + self.compliance[i * n + j] * step.mean_normal_force_n,
                )?;
                self.jacobian[i * n + j] = finite(
                    (if i == j { 1.0 } else { 0.0 })
                        - self.compliance[i * n + j] * step.mean_force_gap_derivative_n_per_m,
                )?;
            }
        }
        poll(Some(gate))?;
        self.lu.factor(&self.jacobian)?;
        poll(Some(gate))?;
        Ok(())
    }

    /// Absolute inverse times the original force-to-gap uncertainty, retaining
    /// the full cross-bridge response and the actual local conditioning.
    fn accuracy(&mut self, gate: &CancelGate) -> Result<(bool, f64), Error> {
        let n = self.points.len();
        // Preserve the residual in candidate_gaps while reusing the RHS buffer.
        self.candidate_gaps.copy_from_slice(&self.residual);
        for j in 0..n {
            poll(Some(gate))?;
            self.residual.fill(0.0);
            self.residual[j] = 1.0;
            self.lu
                .solve_factored_into(&self.residual, &mut self.correction)?;
            for i in 0..n {
                self.inverse[i * n + j] = self.correction[i].abs();
            }
        }
        self.residual.copy_from_slice(&self.candidate_gaps);
        let mut gain = 0.0_f64;
        for i in 0..n {
            let mut uncertainty = 0.0;
            for j in 0..n {
                uncertainty = finite(
                    uncertainty
                        + self.compliance[i * n + j].abs()
                            * self.steps[j]
                                .as_ref()
                                .expect("complete original probe")
                                .estimated_force_error_n,
                )?;
            }
            self.correction[i] = uncertainty;
        }
        let mut accepted = true;
        for i in 0..n {
            gain = gain.max(finite(self.inverse[i * n..(i + 1) * n].iter().sum())?);
            self.gap_errors[i] = dot(&self.inverse[i * n..(i + 1) * n], &self.correction)?;
            let tolerance = self.gap_tolerance(i, self.gaps[i])?;
            if self.gap_errors[i] > tolerance {
                return Err(Error::Bridge {
                    index: i,
                    error: CapillaryModalError::GapAccuracy {
                        residual_m: self.residual[i],
                        estimated_error_m: self.gap_errors[i],
                        tolerance_m: tolerance,
                    },
                });
            }
            accepted &= finite(self.residual[i].abs() + self.gap_errors[i])? <= tolerance;
        }
        if gain > self.config.maximum_compatibility_gain {
            return Err(Error::CompatibilityGain {
                gain,
                maximum: self.config.maximum_compatibility_gain,
            });
        }
        Ok((accepted, gain))
    }

    pub(super) fn solve(&mut self, gate: &CancelGate) -> Result<(usize, usize, f64), Error> {
        let n = self.points.len();
        // Reserve one final original-source audit for EACH bridge.
        let ceiling = self.config.maximum_evaluations - n;
        let mut evaluations = 0;
        Self::evaluate(
            &self.points,
            &self.gaps,
            &mut self.steps,
            self.network.dt,
            self.config.maximum_gap_increment_m,
            &mut evaluations,
            ceiling,
            gate,
        )?;
        for iteration in 0..=self.config.maximum_iterations {
            poll(Some(gate))?;
            let merit = self.residual_norm(&self.gaps, &self.steps)?;
            self.factor(gate)?;
            if merit <= 1.0 {
                let (accepted, gain) = self.accuracy(gate)?;
                if accepted {
                    return Ok((iteration, evaluations, gain));
                }
            }
            if iteration == self.config.maximum_iterations {
                return Err(Error::NotConverged);
            }
            self.lu
                .solve_factored_into(&self.residual, &mut self.correction)?;
            let mut fraction = 1.0;
            let mut accepted = false;
            for _ in 0..self.config.maximum_backtracks {
                poll(Some(gate))?;
                for i in 0..n {
                    self.candidate_gaps[i] = finite(self.gaps[i] + fraction * self.correction[i])?;
                }
                match Self::evaluate(
                    &self.points,
                    &self.candidate_gaps,
                    &mut self.candidate_steps,
                    self.network.dt,
                    self.config.maximum_gap_increment_m,
                    &mut evaluations,
                    ceiling,
                    gate,
                ) {
                    Ok(()) => {
                        let next =
                            self.residual_norm(&self.candidate_gaps, &self.candidate_steps)?;
                        if next < merit || next <= 0.125 {
                            core::mem::swap(&mut self.gaps, &mut self.candidate_gaps);
                            core::mem::swap(&mut self.steps, &mut self.candidate_steps);
                            accepted = true;
                            break;
                        }
                    }
                    Err(Error::Bridge {
                        error: CapillaryModalError::ContactLine(ContactLineError::Domain(_)),
                        ..
                    }) => {}
                    Err(e) => return Err(e),
                }
                fraction *= 0.5;
            }
            if !accepted {
                return Err(Error::NotConverged);
            }
        }
        Err(Error::NotConverged)
    }
}
