//! Causal, positive real spatial admittance on the original surface traces.
//!
//! Each authored real pole has `tau q_dot + q = p` and outward velocity
//! `a(x) q`, where a >= 0 and tau > 0. Together with Y0(x) >= 0 this gives
//! `Y(x,s) = Y0(x) + sum a_j(x)/(1+s tau_j)`. It is a parallel family of
//! resistive/inertive branches; it does not realize arbitrary measured wall
//! data, complex poles, elastic resonance or lateral wall coupling.
//! The physical memory stores `tau integral(a q²)/2` and dissipates
//! `integral(a q²)`. A real harmonic sample alone cannot supply this law.
//!
//! Auxiliary differential boundary states and passive rational admittances:
//! Bilbao and Hamilton, ICA 2016, "Passive time-domain numerical designs for
//! room acoustics", https://www.pure.ed.ac.uk/ws/portalfiles/portal/28188217/ICA_2016_bilbaohamilton.pdf.
//! This owner admits the explicit positive-residue, real-pole subset only.

use super::super::{HeterogeneousError, HeterogeneousHelmholtzFe, invalid, poll};
use super::{SpatialBoundaryAdmittance, SpatialResistiveBoundary};
use fs_sparse::{Coo, Csr};

/// Maximum explicitly authored real relaxation poles on one physical patch.
pub const MAX_ADMITTANCE_RELAXATION_TERMS: usize = 64;

/// A causal wall law with original spatial residue forms and real stable poles.
///
/// Every component has identical original faces, geometry/material ownership,
/// and declared quadrature order. Coefficients need not be proportional, so
/// their spatial pressure correlations remain independent. Preparation retains
/// the sampled forms; neither harmonic reuse nor transient stepping resamples.
pub struct SpatialRelaxingBoundary<'m, 'f, 'c> {
    pub(in crate::acoustics::helmholtz::heterogeneous) instantaneous:
        SpatialResistiveBoundary<'m, 'f, 'c>,
    pub(in crate::acoustics::helmholtz::heterogeneous) terms: Vec<AdmittanceRelaxation>,
}

pub(in crate::acoustics::helmholtz::heterogeneous) struct AdmittanceRelaxation {
    pub conductance: Csr,
    pub relaxation_time_s: f64,
}

impl<'f, 'c> HeterogeneousHelmholtzFe<'f, 'c> {
    /// Compose an explicitly authored causal spatial admittance.
    ///
    /// `instantaneous` is Y0(x), or zero when absent. Each term supplies the
    /// nonnegative amplitude a_j(x) [m/(Pa s)] and its positive constant tau_j
    /// [s]. The sampled coefficient fields must belong to this exact model and
    /// have identical faces and quadrature orders. No fitting, pole inference,
    /// frequency-local extrapolation or material callback occurs here.
    /// A spatially varying relaxation time must be resolved into separately
    /// authored patches/spectra; it is not averaged by this owner.
    pub fn spatial_relaxing_boundary<'m>(
        &'m self,
        instantaneous: Option<&SpatialResistiveBoundary<'_, '_, '_>>,
        terms: &[(&SpatialResistiveBoundary<'_, '_, '_>, f64)],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialRelaxingBoundary<'m, 'f, 'c>, HeterogeneousError> {
        poll(cancelled)?;
        if terms.is_empty() || terms.len() > MAX_ADMITTANCE_RELAXATION_TERMS {
            return Err(invalid("one to 64 explicit spatial admittance poles required").into());
        }
        let reference = terms[0].0;
        reference.validate(self)?;
        for (field, tau) in terms {
            poll(cancelled)?;
            field.validate(self)?;
            if field.faces != reference.faces
                || field.quadrature_order != reference.quadrature_order
            {
                return Err(invalid(
                    "spatial admittance components require identical original faces and quadrature",
                )
                .into());
            }
            if !tau.is_finite() || *tau <= 0.0 {
                return Err(
                    invalid("positive finite spatial admittance relaxation time required").into(),
                );
            }
        }
        if let Some(field) = instantaneous {
            field.validate(self)?;
            if field.faces != reference.faces
                || field.quadrature_order != reference.quadrature_order
            {
                return Err(invalid(
                    "instantaneous admittance must use the same original patch and quadrature",
                )
                .into());
            }
        }
        let instantaneous = SpatialResistiveBoundary {
            model: self,
            quadrature_order: reference.quadrature_order,
            faces: reference.faces.clone(),
            conductance: instantaneous.map_or_else(
                || Coo::new(self.ndof(), self.ndof()).assemble(),
                |field| field.conductance.clone(),
            ),
        };
        let mut retained = Vec::with_capacity(terms.len());
        for &(field, relaxation_time_s) in terms {
            poll(cancelled)?;
            retained.push(AdmittanceRelaxation {
                conductance: field.conductance.clone(),
                relaxation_time_s,
            });
        }
        poll(cancelled)?;
        Ok(SpatialRelaxingBoundary {
            instantaneous,
            terms: retained,
        })
    }
}

impl<'m, 'f, 'c> SpatialRelaxingBoundary<'m, 'f, 'c> {
    /// Original canonical boundary triangles.
    pub fn faces(&self) -> &[[u32; 3]] {
        self.instantaneous.faces()
    }

    /// Original common surface quadrature order.
    pub fn quadrature_order(&self) -> usize {
        self.instantaneous.quadrature_order()
    }

    /// Number of explicitly authored memory branches.
    pub fn relaxation_terms(&self) -> usize {
        self.terms.len()
    }

    /// Original constant relaxation time of a branch [s].
    pub fn relaxation_time_s(&self, term: usize) -> Option<f64> {
        self.terms.get(term).map(|arm| arm.relaxation_time_s)
    }

    /// Reuse the SAME causal law in exp(+i omega t) harmonic acoustics.
    /// Both absorption and reactance are retained. No fitted time law is
    /// reconstructed from these frequency samples.
    pub fn at_frequency(
        &self,
        frequency_hz: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialBoundaryAdmittance<'m, 'f, 'c>, HeterogeneousError> {
        poll(cancelled)?;
        let mut result = self.instantaneous.at_frequency(frequency_hz)?;
        let omega = core::f64::consts::TAU * frequency_hz;
        let n = self.instantaneous.model.ndof();
        let mut real = Coo::new(n, n);
        let mut imaginary = Coo::new(n, n);
        add_form(
            &mut real,
            &self.instantaneous.conductance,
            n,
            1.0,
            cancelled,
        )?;
        for term in &self.terms {
            poll(cancelled)?;
            let x = omega * term.relaxation_time_s;
            // Avoid overflowing x²; the inverse product also handles a large
            // but finite omega/tau pair whose product itself overflows.
            let (re, im) = if x <= 1.0 {
                let denominator = 1.0 + x * x;
                (1.0 / denominator, -x / denominator)
            } else {
                let inverse = (1.0 / omega) / term.relaxation_time_s;
                let denominator = 1.0 + inverse * inverse;
                (inverse * inverse / denominator, -inverse / denominator)
            };
            add_form(&mut real, &term.conductance, n, re, cancelled)?;
            add_form(&mut imaginary, &term.conductance, n, im, cancelled)?;
        }
        result.conductance = finish_form(real, n, cancelled)?;
        result.susceptance = finish_form(imaginary, n, cancelled)?;
        poll(cancelled)?;
        Ok(result)
    }

    pub(in crate::acoustics::helmholtz::heterogeneous) fn midpoint_condensate(
        &self,
        duration_s: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialResistiveBoundary<'m, 'f, 'c>, HeterogeneousError> {
        let n = self.instantaneous.model.ndof();
        let mut combined = Coo::new(n, n);
        add_form(
            &mut combined,
            &self.instantaneous.conductance,
            n,
            1.0,
            cancelled,
        )?;
        for term in &self.terms {
            let (_, pressure_weight) = midpoint_weights(duration_s, term.relaxation_time_s)?;
            add_form(
                &mut combined,
                &term.conductance,
                n,
                pressure_weight,
                cancelled,
            )?;
        }
        Ok(SpatialResistiveBoundary {
            model: self.instantaneous.model,
            quadrature_order: self.quadrature_order(),
            faces: self.faces().to_vec(),
            conductance: finish_form(combined, n, cancelled)?,
        })
    }
}

pub(in crate::acoustics::helmholtz::heterogeneous) fn midpoint_weights(
    duration_s: f64,
    tau: f64,
) -> Result<(f64, f64), HeterogeneousError> {
    let half = 0.5 * duration_s;
    if !half.is_finite() || half <= 0.0 {
        return Err(invalid("positive representable admittance time step required").into());
    }
    let scale = half.max(tau);
    let a = tau / scale;
    let b = half / scale;
    let history = a / (a + b);
    let pressure = b / (a + b);
    if ![history, pressure]
        .into_iter()
        .all(|x| x.is_finite() && x > 0.0)
    {
        return Err(invalid("admittance memory time scales are unrepresentable").into());
    }
    Ok((history, pressure))
}

fn add_form(
    output: &mut Coo,
    matrix: &Csr,
    n: usize,
    factor: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), HeterogeneousError> {
    for row in 0..n {
        poll(cancelled)?;
        let (cols, entries) = matrix.row(row);
        for (&col, &entry) in cols.iter().zip(entries) {
            let value = factor * entry;
            if !value.is_finite() {
                return Err(invalid("spatial causal admittance form overflow").into());
            }
            if value != 0.0 {
                output.push(row, col, value);
            }
        }
    }
    Ok(())
}

fn finish_form(
    input: Coo,
    n: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Csr, HeterogeneousError> {
    let result = input.assemble();
    for row in 0..n {
        poll(cancelled)?;
        if result.row(row).1.iter().any(|x| !x.is_finite()) {
            return Err(invalid("spatial causal admittance sum overflow").into());
        }
    }
    Ok(result)
}
