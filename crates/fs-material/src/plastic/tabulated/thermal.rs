//! Shared affine-temperature arithmetic for an explicitly interpreted yield table.
//!
//! The consuming constitutive law owns the stress measure, reference temperature
//! and admitted thermal interval. This immutable preparation owns only the fixed
//! plastic-strain grid, affine yield ordinates and exact hardening-storage
//! derivative. It never converts Cauchy and Kirchhoff stress or selects sources.

use super::{HardeningPoint, TabulatedError, TabulatedHardening};
use std::sync::Arc;

/// `Y_i(T) = Y_i(T_ref) + s_i (T-T_ref)` on the original fixed strain knots.
/// The associated energetic storage is `integral_0^alpha (Y(p,T)-Y(0,T)) dp`.
#[derive(Clone, Debug, PartialEq)]
pub struct AffineTabulatedHardening {
    prepared: Arc<Prepared>,
}

#[derive(Debug, PartialEq)]
struct Prepared {
    reference: TabulatedHardening,
    yield_slopes_pa_k: Vec<f64>,
    slope_gradients_pa_k: Vec<f64>,
    storage_slopes_j_m3_k: Vec<f64>,
    maximum_absolute_storage_slope_j_m3_k: f64,
}

impl AffineTabulatedHardening {
    /// Prepare one finite derivative for each original absolute yield ordinate.
    /// The consumer must also admit both thermal endpoint curves and its full
    /// native history bounds before using this family in constitutive updates.
    pub fn new(
        reference: TabulatedHardening,
        yield_stress_slopes_pa_k: Vec<f64>,
    ) -> Result<Self, TabulatedError> {
        let points = reference.points();
        if yield_stress_slopes_pa_k.len() != points.len()
            || yield_stress_slopes_pa_k
                .iter()
                .any(|slope| !slope.is_finite())
        {
            return Err(TabulatedError::Parameters(
                "each fixed plastic-strain knot requires one finite absolute-yield temperature derivative",
            ));
        }
        let origin = yield_stress_slopes_pa_k[0];
        let mut gradients = Vec::with_capacity(points.len() - 1);
        let mut storage = Vec::with_capacity(points.len());
        let mut maximum = 0.0_f64;
        storage.push(0.0);
        for index in 0..points.len() - 1 {
            let width = points[index + 1].plastic_strain - points[index].plastic_strain;
            let lower = finite(yield_stress_slopes_pa_k[index] - origin)?;
            let upper = finite(yield_stress_slopes_pa_k[index + 1] - origin)?;
            let rise =
                finite(yield_stress_slopes_pa_k[index + 1] - yield_stress_slopes_pa_k[index])?;
            let gradient = finite(rise / width)?;
            if rise != 0.0 && gradient == 0.0 {
                return Err(TabulatedError::Numerical(
                    "tabulated temperature derivative is below representable slope resolution",
                ));
            }
            let start = storage[index];
            let end = finite(start + finite((0.5 * lower + 0.5 * upper) * width)?)?;
            maximum = maximum.max(end.abs());
            // A signed primitive can reach its largest magnitude inside a
            // segment even when the stored derivative vanishes at both knots.
            if (lower < 0.0 && upper > 0.0) || (lower > 0.0 && upper < 0.0) {
                let distance = -lower / gradient;
                let extremum = finite(start + finite((0.5 * lower) * distance)?)?;
                maximum = maximum.max(extremum.abs());
            }
            gradients.push(gradient);
            storage.push(end);
        }
        Ok(Self {
            prepared: Arc::new(Prepared {
                reference,
                yield_slopes_pa_k: yield_stress_slopes_pa_k,
                slope_gradients_pa_k: gradients,
                storage_slopes_j_m3_k: storage,
                maximum_absolute_storage_slope_j_m3_k: maximum,
            }),
        })
    }

    /// Original curve at the consuming family's explicit reference temperature.
    pub fn reference(&self) -> &TabulatedHardening {
        &self.prepared.reference
    }

    /// Explicit derivatives of the original absolute yield ordinates [Pa/K].
    pub fn yield_stress_slopes_pa_k(&self) -> &[f64] {
        &self.prepared.yield_slopes_pa_k
    }

    /// Largest magnitude of the exact signed storage derivative on the full
    /// original strain interval, including extrema between strain knots.
    pub fn maximum_absolute_storage_slope_j_m3_k(&self) -> f64 {
        self.prepared.maximum_absolute_storage_slope_j_m3_k
    }

    /// Fixed-history derivative of recoverable hardening energy [J/(m³ K)].
    /// Its negative contributes to entropy. This bounded query allocates nothing.
    pub fn hardening_energy_temperature_derivative_j_m3_k(
        &self,
        alpha: f64,
    ) -> Result<f64, TabulatedError> {
        let data = &self.prepared;
        let points = data.reference.points();
        let maximum = data.reference.maximum_plastic_strain();
        if !alpha.is_finite() || alpha < 0.0 || alpha > maximum {
            return Err(TabulatedError::Domain {
                plastic_strain: alpha,
                maximum,
            });
        }
        if alpha == maximum {
            return Ok(data.storage_slopes_j_m3_k[points.len() - 1]);
        }
        let index = points.partition_point(|point| point.plastic_strain <= alpha) - 1;
        let width = alpha - points[index].plastic_strain;
        finite(
            data.storage_slopes_j_m3_k[index]
                + ((data.yield_slopes_pa_k[index] - data.yield_slopes_pa_k[0])
                    + 0.5 * data.slope_gradients_pa_k[index] * width)
                    * width,
        )
    }

    /// Prepare an endpoint outside material-point iterations. The caller must
    /// first admit this temperature change against its original source interval.
    /// Reference preparation is shared; a changed endpoint checks every ordinate
    /// and segment through the original bounded hardening constructor.
    pub(crate) fn at_temperature_change(
        &self,
        temperature_change_k: f64,
    ) -> Result<TabulatedHardening, TabulatedError> {
        finite(temperature_change_k)?;
        if temperature_change_k == 0.0 {
            return Ok(self.reference().clone());
        }
        let points = self
            .reference()
            .points()
            .iter()
            .zip(self.yield_stress_slopes_pa_k())
            .map(|(point, slope)| {
                Ok(HardeningPoint {
                    plastic_strain: point.plastic_strain,
                    yield_stress_pa: finite(
                        slope.mul_add(temperature_change_k, point.yield_stress_pa),
                    )?,
                })
            })
            .collect::<Result<Vec<_>, TabulatedError>>()?;
        TabulatedHardening::new(points)
    }
}

fn finite(value: f64) -> Result<f64, TabulatedError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(TabulatedError::Numerical(
            "nonfinite affine hardening arithmetic",
        ))
    }
}
