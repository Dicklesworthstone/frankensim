//! Symmetric coordinate increments of the original nonlocal shear operator.
//!
//! A finite temperature change has one contribution per ORIGINAL temperature
//! coordinate. A forward coordinate sweep and a reverse coordinate sweep each
//! telescope to H(T1)-H(T0). Averaging their coordinate contributions gives a
//! time-reversal-symmetric discrete gradient without assigning the stack's free
//! or internal energy to arbitrary layers. This is a declared finite-step path
//! approximation: it is consistent under refinement, but is not an exact
//! straight-path integral or a continuum entropy-production certificate.

use super::{
    LaminateShearThermoelasticError, ResolvedLaminateShearThermoelastic, check, evaluate, finite,
    invalid,
};

/// A finite operator increment resolved by original temperature coordinate.
///
/// For endpoint shear strains `g0`, `g1` and reference area A, the parameter
/// work belonging to coordinate i is `A g0^T delta_H_i g1 / 2`. Together with
/// the endpoint entropy conjugates it gives the caloric transfer
/// `Q_i = -A g0^T delta_H_i g1 / 2 - delta(T_i S_i)`.
/// The sum is mechanical shear work minus the change in the one stack internal
/// energy. This carrier never distributes a solver or floating-point residual
/// into those physical transfers.
#[derive(Clone, Debug, PartialEq)]
pub struct LaminateShearTemperatureIncrement {
    coordinate_increments: Vec<[f64; 4]>,
    residual: [f64; 4],
}

impl LaminateShearTemperatureIncrement {
    /// Symmetric coordinate increments [Pa m], in original spatial layer order.
    /// Each matrix uses the section xz/yz axes, row major. An unchanged actual
    /// temperature has an exactly zero coordinate contribution.
    pub fn stiffness_coordinate_increments_pa_m(&self) -> &[[f64; 4]] {
        &self.coordinate_increments
    }

    /// Floating-point `sum_i delta_H_i - (H1-H0)` [Pa m]. This is bounded on
    /// construction and is diagnostic only; it is not allocated to any layer.
    pub const fn stiffness_increment_residual_pa_m(&self) -> &[f64; 4] {
        &self.residual
    }
}

impl ResolvedLaminateShearThermoelastic {
    /// Resolve a finite shear-operator change by symmetric coordinate sweeps.
    ///
    /// Endpoints must share the exact original source owner, geometry,
    /// reference receipts, admitted temperature intervals and strain envelopes.
    /// Thus every hybrid corner uses an already admitted original source chart;
    /// independently rebound materials or different affine intervals refuse.
    /// Each sweep changes a layer's cached exact coefficient record, then
    /// evaluates the original selected closure with its retained physical axes.
    /// It needs O(n²) scalar work, O(n) storage, no source re-query and no new
    /// derivative preparation.
    /// Prepare once per endpoint pair and reuse across spatial quadrature.
    pub fn temperature_increment(
        &self,
        end: &Self,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<LaminateShearTemperatureIncrement, LaminateShearThermoelasticError> {
        check(&mut cancelled)?;
        if !core::ptr::eq(
            self.current.source().layers(),
            end.current.source().layers(),
        ) || !core::ptr::eq(
            self.current.source().reference_section(),
            end.current.source().reference_section(),
        ) || self.layers.len() != end.layers.len()
        {
            return Err(invalid(
                "shear temperature increments require the same original source owner and geometry",
            ));
        }
        let mut changed = Vec::with_capacity(self.layers.len());
        for (i, (before, after)) in self.layers.iter().zip(&end.layers).enumerate() {
            check(&mut cancelled)?;
            if before.temperature_range.map(f64::to_bits)
                != after.temperature_range.map(f64::to_bits)
                || before.strain_limit.to_bits() != after.strain_limit.to_bits()
            {
                return Err(invalid(
                    "shear increment endpoints need identical original affine intervals and strain envelopes",
                ));
            }
            if self.temperatures[i] != end.temperatures[i] {
                changed.push(i);
            }
        }
        let mut coordinate_increments = vec![[0.0; 4]; self.layers.len()];
        if changed.is_empty() {
            if self.stiffness != end.stiffness {
                return Err(invalid(
                    "unchanged source temperatures cannot change the native shear operator",
                ));
            }
            check(&mut cancelled)?;
            return Ok(LaminateShearTemperatureIncrement {
                coordinate_increments,
                residual: [0.0; 4],
            });
        }

        let mut roundoff_scale = core::array::from_fn::<_, 4, _>(|k| {
            self.stiffness[k].abs().max(end.stiffness[k].abs())
        });
        for reverse in [false, true] {
            let mut layers = self.layers.clone();
            let mut previous = self.stiffness;
            for step in 0..changed.len() {
                check(&mut cancelled)?;
                let coordinate = changed[if reverse {
                    changed.len() - 1 - step
                } else {
                    step
                }];
                layers[coordinate] = end.layers[coordinate].clone();
                // The native source endpoints are authoritative. Intermediate
                // corners use the same smooth expression whose agreement with
                // those native operators was checked at endpoint preparation.
                let next = if step + 1 == changed.len() {
                    end.stiffness
                } else {
                    evaluate(
                        &layers,
                        self.current.prepared().section().shear_model(),
                        &|_, coefficients| coefficients[0],
                        &mut cancelled,
                    )?
                };
                for k in 0..4 {
                    let delta = finite(next[k] - previous[k])?;
                    coordinate_increments[coordinate][k] =
                        finite(coordinate_increments[coordinate][k] + 0.5 * delta)?;
                    roundoff_scale[k] = finite(roundoff_scale[k] + delta.abs())?;
                }
                previous = next;
            }
        }

        let mut residual = [0.0; 4];
        for k in 0..4 {
            let mut sum = 0.0;
            for increment in &coordinate_increments {
                check(&mut cancelled)?;
                sum = finite(sum + increment[k])?;
            }
            residual[k] = finite(sum - finite(end.stiffness[k] - self.stiffness[k])?)?;
            let tolerance = finite(
                128.0 * f64::EPSILON * (self.layers.len() as f64 + 1.0) * roundoff_scale[k],
            )?;
            if residual[k].abs() > tolerance {
                return Err(invalid(
                    "symmetric shear coordinate increments fail the operator telescoping identity",
                ));
            }
        }
        check(&mut cancelled)?;
        Ok(LaminateShearTemperatureIncrement {
            coordinate_increments,
            residual,
        })
    }
}
