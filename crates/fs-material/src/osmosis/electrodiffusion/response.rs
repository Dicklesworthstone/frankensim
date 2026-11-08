//! Analytic response of the original open-circuit membrane law.
//!
//! At fixed voltage, `dJ_i=kL_i dcL_i-kR_i dcR_i+h_i [B(dcL-dcR)]_i`,
//! where `h_i=dJ_i/dA_i>=0`. Eliminating the same floating electric field gives
//! `dpsi=-sum(z_i dJ_i)/sum(z_i² h_i)` and `dJ_i+=z_i h_i dpsi`.
//! Neither reciprocal concentrations nor numerical perturbations are required.

use super::*;
use crate::osmosis::nonideal::VirialSolutionModel;

/// Signed inventory/load derivatives with respect to one common parameter.
/// These are tangent variations, NOT new physical endpoint values. Temperature,
/// membrane coefficients, fixed charge and separately trapped amounts stay fixed.
#[derive(Clone, Copy, Debug)]
pub struct ElectrodiffusionDirection<'a> {
    pub left_amounts_mol: &'a [f64],
    pub right_amounts_mol: &'a [f64],
    pub volumes_m3: [f64; 2],
    pub pressure_drop_pa: f64,
}

/// Directional derivatives per unit of the caller's common parameter.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ElectrodiffusionDirectionalReport {
    pub volume_flux_derivative_m_s: f64,
    pub osmotic_pressure_drop_derivative_pa: f64,
    pub driving_pressure_derivative_pa: f64,
    /// None on an infinite/indeterminate branch restricted to fixed conducting
    /// support. No finite derivative of the singular voltage is inferred.
    pub dimensionless_potential_derivative: Option<f64>,
    /// Actual compensated sum of z_i times the returned ionic flux derivatives.
    pub charge_flux_derivative_mol_m2_s: f64,
}

/// Compact immutable response retaining the ORIGINAL model and one electric solve.
/// No dense Jacobian or heap allocation is retained. Sparse concentration
/// directions cost O(N + N*nnz); a dense direction costs O(N²).
#[derive(Debug)]
pub struct ElectrodiffusionResponse<'a> {
    law: &'a ElectrodiffusiveMembrane,
    model: Option<&'a VirialSolutionModel>,
    face: Face,
    volumes: [f64; 2],
    trapped: [f64; 2],
    forward: [f64; MAX_OSMOTIC_SOLUTES],
    reverse: [f64; MAX_OSMOTIC_SOLUTES],
    slope: [f64; MAX_OSMOTIC_SOLUTES],
    unavailable_absent: [bool; MAX_OSMOTIC_SOLUTES],
    pressure_gradient_rt: [[f64; MAX_OSMOTIC_SOLUTES]; 2],
    charge_scale: f64,
    slope_scale: f64,
    susceptibility: f64,
    report: ElectrodiffusionReport,
}

impl ElectrodiffusiveMembrane {
    /// Evaluate flux and prepare its analytic response using one electrical solve.
    /// Failure or cancellation leaves the caller's complete output unchanged.
    pub fn flux_response_into<'a>(
        &'a self,
        problem: ElectrodiffusionProblem<'a>,
        controls: ElectrodiffusionControls,
        output: &mut [f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<ElectrodiffusionResponse<'a>> {
        self.flux_response_into_counted(problem, controls, output, &mut 0, cancelled)
    }

    /// Rejected electrical work remains in the shared counter. Repeated tangent
    /// applications do not solve or change the potential and consume no further
    /// potential evaluations. All original flux/thermodynamic admission applies.
    pub fn flux_response_into_counted<'a>(
        &'a self,
        problem: ElectrodiffusionProblem<'a>,
        controls: ElectrodiffusionControls,
        output: &mut [f64],
        potential_evaluations: &mut usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ElectrodiffusionResponse<'a>> {
        let n = self.valences.len();
        shape(n, output.len())?;
        let mut flux = [0.0; MAX_OSMOTIC_SOLUTES];
        let report = self.flux_into_counted(
            problem,
            controls,
            &mut flux[..n],
            potential_evaluations,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        let mut response = ElectrodiffusionResponse {
            law: self,
            model: problem.left.nonideal_model(),
            face: self.prepare(problem)?,
            volumes: [
                problem.left.solvent_volume_m3(),
                problem.right.solvent_volume_m3(),
            ],
            trapped: problem.trapped_particles_mol,
            forward: [0.0; MAX_OSMOTIC_SOLUTES],
            reverse: [0.0; MAX_OSMOTIC_SOLUTES],
            slope: [0.0; MAX_OSMOTIC_SOLUTES],
            unavailable_absent: [false; MAX_OSMOTIC_SOLUTES],
            pressure_gradient_rt: [[1.0; MAX_OSMOTIC_SOLUTES]; 2],
            charge_scale: 1.0,
            slope_scale: 0.0,
            susceptibility: 0.0,
            report,
        };
        let psi = match report.potential {
            DonnanPotential::Finite {
                dimensionless_left_minus_right,
                ..
            } => Some(dimensionless_left_minus_right),
            _ => None,
        };
        for i in 0..n {
            poll(&mut cancelled)?;
            let z = self.valences[i];
            let p = self.permeances[i].value();
            if p > 0.0 && (psi.is_some() || z == 0) {
                let (left, right) = (response.face.left[i], response.face.right[i]);
                let channel = (|| -> Result<_> {
                    let a = finite(
                        response.face.excess_drop[i] + f64::from(z) * psi.unwrap_or(0.0),
                        "channel response drift is not representable",
                    )?;
                    let (forward, reverse, weight) = channel_coefficients(p, a)?;
                    let slope = finite(
                        p * if a >= 0.0 {
                            (1.0 - weight) * left + weight * right
                        } else {
                            weight * left + (1.0 - weight) * right
                        },
                        "channel response susceptibility is not representable",
                    )?;
                    Ok((forward, reverse, slope))
                })();
                match channel {
                    Ok((forward, reverse, slope)) => {
                        response.forward[i] = forward;
                        response.reverse[i] = reverse;
                        response.slope[i] = slope;
                    }
                    // Exact absent pools are harmless on fixed-support paths.
                    // Their unavailable donor derivative is never fabricated.
                    Err(ElectrodiffusionError::Material(OsmosisError::Arithmetic(_)))
                        if left == 0.0 && right == 0.0 =>
                    {
                        response.unavailable_absent[i] = true
                    }
                    Err(error) => return Err(error),
                }
            }
            if response.slope[i] > 0.0 && z != 0 {
                response.charge_scale = response.charge_scale.max(f64::from(z).abs());
                response.slope_scale = response.slope_scale.max(response.slope[i]);
            }
            if let Some(model) = response.model {
                let (mut left, mut right) = (Sum::default(), Sum::default());
                left.add(1.0);
                right.add(1.0);
                for j in 0..n {
                    let b = model.interaction_matrix_m3_mol()[i * n + j];
                    left.add(b * response.face.left[j]);
                    right.add(b * response.face.right[j]);
                }
                response.pressure_gradient_rt[0][i] = left.finish()?;
                response.pressure_gradient_rt[1][i] = right.finish()?;
            }
        }
        if response.slope_scale > 0.0 {
            let mut susceptibility = Sum::default();
            for i in 0..n {
                if self.valences[i] == 0 {
                    continue;
                }
                let z = f64::from(self.valences[i]) / response.charge_scale;
                susceptibility.add(z * z * (response.slope[i] / response.slope_scale));
            }
            response.susceptibility = susceptibility.finish()?;
        }
        poll(&mut cancelled)?;
        output.copy_from_slice(&flux[..n]);
        Ok(response)
    }
}

impl ElectrodiffusionResponse<'_> {
    pub const fn report(&self) -> &ElectrodiffusionReport {
        &self.report
    }
    pub fn species_count(&self) -> usize {
        self.law.valences.len()
    }

    /// Differentiate all original species, water and the floating field together.
    /// Finite-potential zero endpoints are supported without log derivatives.
    /// At an infinite/indeterminate voltage, a direction must keep every absent
    /// conducting charged endpoint absent; otherwise SingularDerivative is
    /// returned. Blocked ions may vary and still affect original water/activity
    /// forces. Neutral-only and water-only fixed-support directions remain valid.
    /// An absent channel with unrepresentable donor coefficients likewise admits
    /// only directions keeping that pool exactly absent.
    /// Output is atomic on refusal; every call is allocation-free.
    pub fn directional_derivative_into(
        &self,
        direction: ElectrodiffusionDirection<'_>,
        output: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ElectrodiffusionDirectionalReport> {
        let n = self.species_count();
        for actual in [
            direction.left_amounts_mol.len(),
            direction.right_amounts_mol.len(),
            output.len(),
        ] {
            shape(n, actual)?;
        }
        if direction
            .left_amounts_mol
            .iter()
            .chain(direction.right_amounts_mol)
            .any(|d| !d.is_finite())
            || direction.volumes_m3.iter().any(|d| !d.is_finite())
            || !direction.pressure_drop_pa.is_finite()
        {
            return Err(ElectrodiffusionError::Invalid(
                "finite original inventory/load direction required",
            ));
        }
        poll(&mut cancelled)?;
        let regular = matches!(self.report.potential, DonnanPotential::Finite { .. });
        let mut dc = [[0.0; MAX_OSMOTIC_SOLUTES]; 2];
        let mut difference = [0.0; MAX_OSMOTIC_SOLUTES];
        let mut pressure = Sum::default();
        for (side, amounts) in [direction.left_amounts_mol, direction.right_amounts_mol]
            .into_iter()
            .enumerate()
        {
            let concentrations = if side == 0 {
                &self.face.left
            } else {
                &self.face.right
            };
            for i in 0..n {
                if self.unavailable_absent[i] && amounts[i] != 0.0 {
                    return Err(ElectrodiffusionError::SingularDerivative);
                }
                if !regular
                    && self.law.valences[i] != 0
                    && self.law.permeances[i].value() > 0.0
                    && concentrations[i] == 0.0
                    && amounts[i] != 0.0
                {
                    return Err(ElectrodiffusionError::SingularDerivative);
                }
                dc[side][i] = finite(
                    (-concentrations[i]).mul_add(direction.volumes_m3[side], amounts[i])
                        / self.volumes[side],
                    "inventory concentration tangent is not representable",
                )?;
                pressure.add(
                    (if side == 0 { 1.0 } else { -1.0 })
                        * self.pressure_gradient_rt[side][i]
                        * dc[side][i],
                );
            }
            pressure.add(
                (if side == 0 { -1.0 } else { 1.0 })
                    * (self.trapped[side] / self.volumes[side])
                    * (direction.volumes_m3[side] / self.volumes[side]),
            );
        }
        for i in 0..n {
            difference[i] = finite(
                dc[0][i] - dc[1][i],
                "concentration-difference tangent is not representable",
            )?;
        }
        let mut excess: [Sum; MAX_OSMOTIC_SOLUTES] = core::array::from_fn(|_| Sum::default());
        if let Some(model) = self.model {
            // A charge-preserving Newton column changes only a free ion and a
            // pivot. Skip zero columns so assembling all such columns is O(N²).
            for j in 0..n {
                if difference[j] == 0.0 {
                    continue;
                }
                poll(&mut cancelled)?;
                for (i, sum) in excess.iter_mut().enumerate().take(n) {
                    if self.slope[i] != 0.0 {
                        sum.add(model.interaction_matrix_m3_mol()[i * n + j] * difference[j]);
                    }
                }
            }
        }
        let mut candidate = [0.0; MAX_OSMOTIC_SOLUTES];
        let mut scale = 0.0_f64;
        for (i, excess) in excess.into_iter().enumerate().take(n) {
            let mut flux = Sum::default();
            flux.add(self.forward[i] * dc[0][i]);
            flux.add(-self.reverse[i] * dc[1][i]);
            flux.add(self.slope[i] * excess.finish()?);
            candidate[i] = flux.finish()?;
            if self.law.valences[i] != 0 {
                scale = scale.max(candidate[i].abs());
            }
        }
        let potential = if regular {
            if !positive(self.slope_scale) || !positive(self.susceptibility) {
                return Err(ElectrodiffusionError::SingularDerivative);
            }
            let mut charge = Sum::default();
            if scale > 0.0 {
                for (i, value) in candidate.iter().enumerate().take(n) {
                    if self.law.valences[i] == 0 {
                        continue;
                    }
                    let term =
                        (f64::from(self.law.valences[i]) / self.charge_scale) * (*value / scale);
                    if *value != 0.0 && term == 0.0 {
                        return Err(ElectrodiffusionError::Material(OsmosisError::Arithmetic(
                            "scaled ionic charge tangent is not representable",
                        )));
                    }
                    charge.add(term);
                }
            }
            let normalized = -charge.finish()?;
            let dpsi = voltage_tangent(
                scale,
                self.slope_scale,
                normalized,
                self.charge_scale * self.susceptibility,
            )?;
            for (i, value) in candidate.iter_mut().enumerate().take(n) {
                *value = finite(
                    (f64::from(self.law.valences[i]) * self.slope[i]).mul_add(dpsi, *value),
                    "open-circuit ionic flux tangent is not representable",
                )?;
            }
            Some(dpsi)
        } else {
            None
        };
        let mut current = Sum::default();
        for (i, value) in candidate.iter().enumerate().take(n) {
            current.add(f64::from(self.law.valences[i]) * *value);
        }
        let osmotic = finite(
            self.face.rt * pressure.finish()?,
            "osmotic pressure tangent is not representable",
        )?;
        let driving = finite(
            direction.pressure_drop_pa - osmotic,
            "water driving tangent is not representable",
        )?;
        let report = ElectrodiffusionDirectionalReport {
            volume_flux_derivative_m_s: finite(
                self.law.hydraulic.value() * driving,
                "water flux tangent is not representable",
            )?,
            osmotic_pressure_drop_derivative_pa: osmotic,
            driving_pressure_derivative_pa: driving,
            dimensionless_potential_derivative: potential,
            charge_flux_derivative_mol_m2_s: current.finish()?,
        };
        poll(&mut cancelled)?;
        output.copy_from_slice(&candidate[..n]);
        Ok(report)
    }
}

// Scale ratios before multiplication. The logarithmic fallback recovers finite
// tangents when an intermediate scale ratio over/underflows; it never supplies
// a concentration floor or silently replaces an unrepresentable nonzero field.
fn voltage_tangent(
    numerator_scale: f64,
    denominator_scale: f64,
    normalized_charge: f64,
    normalized_conductance: f64,
) -> Result<f64> {
    if numerator_scale == 0.0 || normalized_charge == 0.0 {
        return Ok(0.0);
    }
    let ratio = numerator_scale / denominator_scale;
    let factor = normalized_charge / normalized_conductance;
    let direct = ratio * factor;
    if positive(ratio) && direct.is_finite() && direct != 0.0 {
        return Ok(direct);
    }
    let magnitude = det::exp(
        det::ln(numerator_scale) - det::ln(denominator_scale) + det::ln(normalized_charge.abs())
            - det::ln(normalized_conductance),
    );
    if !positive(magnitude) {
        return Err(ElectrodiffusionError::Material(OsmosisError::Arithmetic(
            "floating voltage tangent is not representable",
        )));
    }
    Ok(magnitude.copysign(normalized_charge))
}
