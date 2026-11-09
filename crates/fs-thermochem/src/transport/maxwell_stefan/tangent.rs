//! Analytic directional response of the original Maxwell–Stefan equations.
//!
//! A prime denotes a derivative with respect to one dimensionless parameter.
//! Species and molecular masses are fixed; the caller owns the differentiable
//! thermodynamic and pair-coefficient sources and their validity domains.

use super::{
    AmountConcentration, Driving, Frame, MAX_MAXWELL_STEFAN_SPECIES_V1, MaxwellStefanErrorV1,
    MaxwellStefanEvaluationV1, MaxwellStefanMolarEvaluationV1, MaxwellStefanMolarSolveV1,
    MaxwellStefanSolveV1, MaxwellStefanV1, MaxwellStefanWorkspaceV1, Solve, filled, poll,
    scaled_product_ratio, sum,
};

/// Operation-tree version of the analytic directional evaluation.
pub const MAXWELL_STEFAN_TANGENT_EVALUATOR_VERSION_V1: u32 = 1;

/// Source derivatives in the owner's canonical species order.
///
/// Derivatives are with respect to one dimensionless parameter, not time unless
/// the caller supplies the appropriate scaling. The composition derivative must
/// lie in the simplex tangent plane; driving-force derivatives must sum to zero.
/// Pair derivatives may have either sign but must be symmetric with zero diagonal.
/// At exact-zero fractions this is a local algebraic derivative. A physical
/// one-sided path and source differentiability remain the caller's responsibility.
#[derive(Clone, Copy, Debug)]
pub struct MaxwellStefanDirectionV1<'a> {
    /// dx_i/ds, dimensionless.
    pub mole_fraction_derivatives: &'a [f64],
    /// dd_i/ds [1/m], including every authored thermodynamic-force derivative.
    pub driving_force_derivatives_per_m: &'a [f64],
    /// dc/ds [mol/m³].
    pub concentration_derivative_mol_m3: f64,
    /// dD_ij/ds [m²/s], row-major; no diagonal or hidden mixing rule.
    pub binary_diffusivity_derivatives_m2_s: &'a [f64],
}

/// The original reusable solver plus two retained primal-flux arrays.
/// Ordinary scalar workspace size and evaluation behavior are unchanged.
#[derive(Clone, Debug)]
pub struct MaxwellStefanTangentWorkspaceV1 {
    solver: MaxwellStefanWorkspaceV1,
    value_mass: Vec<f64>,
    value_molar: Vec<f64>,
}
impl MaxwellStefanTangentWorkspaceV1 {
    /// Allocated numerical capacities, including the original factorization.
    /// Fixed bounded stack scratch and fixed-size owner values are excluded,
    /// following the ordinary workspace convention.
    pub fn workspace_bytes(&self) -> usize {
        self.solver.workspace_bytes()
            + (self.value_mass.capacity() + self.value_molar.capacity())
                * core::mem::size_of::<f64>()
    }
}

/// Value and derivative relative to the mass-averaged reference velocity.
#[derive(Clone, Copy, Debug)]
pub struct MaxwellStefanTangentEvaluationV1<'a> {
    value: MaxwellStefanEvaluationV1<'a>,
    mass_derivative: &'a [f64],
    molar_derivative: &'a [f64],
    gradient_derivative_residual: f64,
    frame_derivative_residual: f64,
}
impl<'a> MaxwellStefanTangentEvaluationV1<'a> {
    /// Original physical value, residuals and pair-friction entropy diagnostic.
    pub const fn value(&self) -> MaxwellStefanEvaluationV1<'a> {
        self.value
    }
    /// dj_i/ds [kg/(m² s)] in the same mass frame.
    pub const fn mass_flux_derivatives(&self) -> &'a [f64] {
        self.mass_derivative
    }
    /// dz_i/ds [mol/(m² s)] in that mass frame.
    pub const fn molar_flux_derivatives(&self) -> &'a [f64] {
        self.molar_derivative
    }
    /// Largest differentiated original-row residual after division by c [1/m].
    pub const fn maximum_gradient_derivative_residual_per_m(&self) -> f64 {
        self.gradient_derivative_residual
    }
    /// Signed sum of the published mass-flux derivatives [kg/(m² s)].
    pub const fn mass_frame_derivative_residual_kg_m2_s(&self) -> f64 {
        self.frame_derivative_residual
    }
}

/// Value and derivative relative to the distinct molar-averaged velocity.
#[derive(Clone, Copy, Debug)]
pub struct MaxwellStefanMolarTangentEvaluationV1<'a> {
    inner: MaxwellStefanTangentEvaluationV1<'a>,
}
impl<'a> MaxwellStefanMolarTangentEvaluationV1<'a> {
    /// Original physical value and diagnostics in the molar frame.
    pub const fn value(&self) -> MaxwellStefanMolarEvaluationV1<'a> {
        MaxwellStefanMolarEvaluationV1 {
            evaluation: self.inner.value,
        }
    }
    /// dz_i/ds [mol/(m² s)], whose sum closes the molar derivative budget.
    pub const fn molar_flux_derivatives(&self) -> &'a [f64] {
        self.inner.molar_derivative
    }
    /// dj_i/ds [kg/(m² s)]; their sum need not vanish in the molar frame.
    pub const fn mass_flux_derivatives(&self) -> &'a [f64] {
        self.inner.mass_derivative
    }
    /// Largest differentiated original-row residual after division by c [1/m].
    pub const fn maximum_gradient_derivative_residual_per_m(&self) -> f64 {
        self.inner.gradient_derivative_residual
    }
    /// Signed sum of the published molar-flux derivatives [mol/(m² s)].
    pub const fn molar_frame_derivative_residual_mol_m2_s(&self) -> f64 {
        self.inner.frame_derivative_residual
    }
}

impl MaxwellStefanV1 {
    /// Allocate the original solver and two O(N) retained value arrays once.
    pub fn tangent_workspace(
        &self,
    ) -> Result<MaxwellStefanTangentWorkspaceV1, MaxwellStefanErrorV1> {
        Ok(MaxwellStefanTangentWorkspaceV1 {
            solver: self.workspace()?,
            value_mass: filled(self.species.len(), 0.0)?,
            value_molar: filled(self.species.len(), 0.0)?,
        })
    }

    /// Evaluate source-owned fluxes and their analytic directional derivative.
    ///
    /// Differentiate A(x,D) z = -c d and the SAME constant-mass frame:
    /// A z' = -c d' - c' d - A' z. The original solver evaluates both right-hand
    /// sides with an identical matrix, reusing exact factors when applicable.
    /// Every source is freshly admitted; exact-absent and trace species remain.
    /// The supplied controls apply independently to the original physical value
    /// and the differentiated equations. A residual refusal may concern either.
    ///
    /// This is a Jacobian-vector operation with respect to the supplied source
    /// data, not an inferred EOS/activity derivative or a derivative of entropy.
    /// No finite difference or extrapolated physical state is evaluated.
    /// Caller inputs stay unchanged; cancellation publishes no partial result.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_driving_forces_tangent<'a>(
        &self,
        mole_fractions: &[f64],
        driving_forces_per_m: &[f64],
        concentration: AmountConcentration,
        binary_diffusivities_m2_s: &[f64],
        direction: MaxwellStefanDirectionV1<'_>,
        controls: MaxwellStefanSolveV1,
        workspace: &'a mut MaxwellStefanTangentWorkspaceV1,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MaxwellStefanTangentEvaluationV1<'a>, MaxwellStefanErrorV1> {
        self.tangent_in_frame(
            mole_fractions,
            driving_forces_per_m,
            concentration,
            binary_diffusivities_m2_s,
            direction,
            Solve {
                relative_tolerance: controls.relative_tolerance,
                absolute_gradient_tolerance_per_m: controls.absolute_gradient_tolerance_per_m,
                absolute_frame_tolerance: controls.absolute_mass_flux_tolerance_kg_m2_s,
                frame: Frame::MassAveraged,
            },
            workspace,
            cancelled,
        )
    }

    /// The same analytic operation with sum(z_i') = 0 in the molar frame.
    /// Source, derivative-domain and error-control responsibilities are identical
    /// to evaluate_driving_forces_tangent; no barycentric conversion is inferred.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_molar_driving_forces_tangent<'a>(
        &self,
        mole_fractions: &[f64],
        driving_forces_per_m: &[f64],
        concentration: AmountConcentration,
        binary_diffusivities_m2_s: &[f64],
        direction: MaxwellStefanDirectionV1<'_>,
        controls: MaxwellStefanMolarSolveV1,
        workspace: &'a mut MaxwellStefanTangentWorkspaceV1,
        cancelled: impl FnMut() -> bool,
    ) -> Result<MaxwellStefanMolarTangentEvaluationV1<'a>, MaxwellStefanErrorV1> {
        self.tangent_in_frame(
            mole_fractions,
            driving_forces_per_m,
            concentration,
            binary_diffusivities_m2_s,
            direction,
            Solve {
                relative_tolerance: controls.relative_tolerance,
                absolute_gradient_tolerance_per_m: controls.absolute_gradient_tolerance_per_m,
                absolute_frame_tolerance: controls.absolute_molar_flux_tolerance_mol_m2_s,
                frame: Frame::MolarAveraged,
            },
            workspace,
            cancelled,
        )
        .map(|inner| MaxwellStefanMolarTangentEvaluationV1 { inner })
    }

    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn tangent_in_frame<'a>(
        &self,
        x: &[f64],
        driving: &[f64],
        concentration: AmountConcentration,
        diffusion: &[f64],
        direction: MaxwellStefanDirectionV1<'_>,
        controls: Solve,
        work: &'a mut MaxwellStefanTangentWorkspaceV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<MaxwellStefanTangentEvaluationV1<'a>, MaxwellStefanErrorV1> {
        poll(&mut cancelled)?;
        let n = self.species.len();
        admit_direction(n, direction, controls, &mut cancelled)?;
        let value = self.evaluate_in_frame(
            x,
            Driving::Thermodynamic {
                forces: driving,
                concentration,
            },
            diffusion,
            controls,
            &mut work.solver,
            &mut cancelled,
        )?;
        work.value_mass.copy_from_slice(value.mass_fluxes);
        work.value_molar.copy_from_slice(value.molar_fluxes);
        let diagnostics = (
            value.maximum_gradient_residual_per_m,
            value.frame_residual,
            value.drag_entropy_production_w_m3_k,
        );
        let c = concentration.value();
        let mut effective = [0.0; MAX_MAXWELL_STEFAN_SPECIES_V1];
        for i in 0..n {
            poll(&mut cancelled)?;
            let mut terms = [0.0; 4 * MAX_MAXWELL_STEFAN_SPECIES_V1 + 2];
            terms[0] = direction.driving_force_derivatives_per_m[i];
            terms[1] = scaled_product_ratio(
                &[direction.concentration_derivative_mol_m3, driving[i]],
                &[c],
            )?;
            for j in 0..n {
                if i != j {
                    let pair =
                        source_terms(i, j, n, x, diffusion, direction, &work.value_molar, c)?;
                    terms[2 + 4 * j..6 + 4 * j].copy_from_slice(&pair);
                }
            }
            effective[i] = sum(terms[..4 * n + 2].iter().copied());
            if !effective[i].is_finite() {
                return Err(MaxwellStefanErrorV1::Arithmetic(
                    "linearized diffusion driving force",
                ));
            }
        }
        let tangent = self.evaluate_in_frame(
            x,
            Driving::Thermodynamic {
                forces: &effective[..n],
                concentration,
            },
            diffusion,
            controls,
            &mut work.solver,
            &mut cancelled,
        )?;
        // Reconstruct the differentiated ORIGINAL equations from published
        // fluxes. A check of only the assembled effective driving vector could
        // hide a source-term cancellation or representation defect.
        let mut maximum = 0.0_f64;
        let mut rows_close = true;
        for i in 0..n {
            poll(&mut cancelled)?;
            let mut terms = [0.0; 6 * MAX_MAXWELL_STEFAN_SPECIES_V1 + 2];
            terms[0] = direction.driving_force_derivatives_per_m[i];
            terms[1] = scaled_product_ratio(
                &[direction.concentration_derivative_mol_m3, driving[i]],
                &[c],
            )?;
            for j in 0..n {
                if i == j {
                    continue;
                }
                terms[2 + 6 * j..6 + 6 * j].copy_from_slice(&source_terms(
                    i,
                    j,
                    n,
                    x,
                    diffusion,
                    direction,
                    &work.value_molar,
                    c,
                )?);
                terms[6 + 6 * j] = scaled_product_ratio(
                    &[x[j], tangent.molar_fluxes[i]],
                    &[c, diffusion[i * n + j]],
                )?;
                terms[7 + 6 * j] = -scaled_product_ratio(
                    &[x[i], tangent.molar_fluxes[j]],
                    &[c, diffusion[i * n + j]],
                )?;
            }
            let residual = sum(terms[..6 * n + 2].iter().copied());
            let scale = sum(terms[..6 * n + 2].iter().map(|term| term.abs()));
            if !residual.is_finite() || !scale.is_finite() {
                return Err(MaxwellStefanErrorV1::Arithmetic(
                    "differentiated diffusion residual",
                ));
            }
            maximum = maximum.max(residual.abs());
            rows_close &= residual.abs()
                <= controls.absolute_gradient_tolerance_per_m + controls.relative_tolerance * scale;
        }
        if !rows_close {
            return Err(match controls.frame {
                Frame::MassAveraged => MaxwellStefanErrorV1::Residual {
                    maximum_gradient_residual_per_m: maximum,
                    mass_frame_residual_kg_m2_s: tangent.frame_residual,
                },
                Frame::MolarAveraged => MaxwellStefanErrorV1::MolarResidual {
                    maximum_gradient_residual_per_m: maximum,
                    molar_frame_residual_mol_m2_s: tangent.frame_residual,
                },
            });
        }
        poll(&mut cancelled)?;
        Ok(MaxwellStefanTangentEvaluationV1 {
            value: MaxwellStefanEvaluationV1 {
                mass_fluxes: &work.value_mass,
                molar_fluxes: &work.value_molar,
                maximum_gradient_residual_per_m: diagnostics.0,
                frame_residual: diagnostics.1,
                drag_entropy_production_w_m3_k: diagnostics.2,
            },
            mass_derivative: tangent.mass_fluxes,
            molar_derivative: tangent.molar_fluxes,
            gradient_derivative_residual: maximum,
            frame_derivative_residual: tangent.frame_residual,
        })
    }
}

fn admit_direction(
    n: usize,
    d: MaxwellStefanDirectionV1<'_>,
    controls: Solve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), MaxwellStefanErrorV1> {
    if d.mole_fraction_derivatives.len() != n
        || d.driving_force_derivatives_per_m.len() != n
        || d.binary_diffusivity_derivatives_m2_s.len() != n * n
        || !d.concentration_derivative_mol_m3.is_finite()
        || d.mole_fraction_derivatives
            .iter()
            .chain(d.driving_force_derivatives_per_m)
            .any(|v| !v.is_finite())
    {
        return Err(MaxwellStefanErrorV1::Invalid(
            "diffusion derivative dimensions or values",
        ));
    }
    let dx_sum = sum(d.mole_fraction_derivatives.iter().copied());
    let dx_scale = d
        .mole_fraction_derivatives
        .iter()
        .map(|v| v.abs())
        .fold(1.0, f64::max);
    if !dx_sum.is_finite() || dx_sum.abs() > 64.0 * f64::EPSILON * n as f64 * dx_scale {
        return Err(MaxwellStefanErrorV1::Invalid(
            "composition derivative leaves the simplex plane",
        ));
    }
    let force_sum = sum(d.driving_force_derivatives_per_m.iter().copied());
    let force_scale = sum(d.driving_force_derivatives_per_m.iter().map(|v| v.abs()));
    if !force_sum.is_finite()
        || !force_scale.is_finite()
        || force_sum.abs()
            > controls.absolute_gradient_tolerance_per_m + controls.relative_tolerance * force_scale
    {
        return Err(MaxwellStefanErrorV1::Invalid(
            "driving derivatives do not sum to zero",
        ));
    }
    for i in 0..n {
        poll(cancelled)?;
        if d.binary_diffusivity_derivatives_m2_s[i * n + i] != 0.0 {
            return Err(MaxwellStefanErrorV1::Invalid(
                "nonzero binary-diffusion derivative diagonal",
            ));
        }
        for j in i + 1..n {
            let value = d.binary_diffusivity_derivatives_m2_s[i * n + j];
            if !value.is_finite() || value != d.binary_diffusivity_derivatives_m2_s[j * n + i] {
                return Err(MaxwellStefanErrorV1::Invalid(
                    "binary-diffusion derivatives are not finite symmetric",
                ));
            }
        }
    }
    Ok(())
}

// Four separately represented terms of (A' z)_i/c for one original pair.
// Scaling all factors before rounding avoids forming D'/D or 1/(c*D²).
#[allow(clippy::too_many_arguments)]
fn source_terms(
    i: usize,
    j: usize,
    n: usize,
    x: &[f64],
    diffusion: &[f64],
    d: MaxwellStefanDirectionV1<'_>,
    z: &[f64],
    concentration: f64,
) -> Result<[f64; 4], MaxwellStefanErrorV1> {
    let pair = diffusion[i * n + j];
    let derivative = d.binary_diffusivity_derivatives_m2_s[i * n + j];
    Ok([
        scaled_product_ratio(
            &[d.mole_fraction_derivatives[j], z[i]],
            &[concentration, pair],
        )?,
        -scaled_product_ratio(
            &[d.mole_fraction_derivatives[i], z[j]],
            &[concentration, pair],
        )?,
        -scaled_product_ratio(&[x[j], z[i], derivative], &[concentration, pair, pair])?,
        scaled_product_ratio(&[x[i], z[j], derivative], &[concentration, pair, pair])?,
    ])
}

#[cfg(test)]
mod tests;
