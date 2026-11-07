//! Affine temperature evolution of the original bounded Kirchhoff yield table.
//!
//! Plastic-strain knots stay fixed. Each absolute yield ordinate evolves as
//! `Y_i(T)=Y_i(T_ref)+s_i*(T-T_ref)`. Admitting the complete table at both ends
//! of the family's temperature interval proves positivity and nondecreasing
//! resistance throughout it. No temperature or plastic-strain extrapolation,
//! interpolation between incompatible strain coordinates, or heat fraction is
//! inferred. Thermal softening changes resistance without erasing history.
//!
//! The energetic convention is the original integral of `Y(alpha,T)-Y(0,T)`.
//! Its temperature derivative is prepared exactly on the original strain grid.
//! Sharing that immutable preparation keeps frozen endpoints cheap and retains
//! entropy even where the entire endpoint curve becomes a yield plateau.
//! Fixed-history free energy is affine in temperature, so this addition has no
//! mechanical heat-capacity term. The existing elastic/expansion bound remains.

use super::{AffineThermalFiniteJ2, FiniteJ2ThermalError, finite};
use crate::plastic::finite::FiniteJ2Error;
use crate::plastic::tabulated::{HardeningPoint, TabulatedError, TabulatedHardening};
use std::sync::Arc;

/// One reference curve and an explicitly supplied derivative for each ordinate.
/// All stresses and derivatives use the native original-volume Kirchhoff basis.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteJ2ThermalTabulatedHardening {
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

impl FiniteJ2ThermalTabulatedHardening {
    /// Prepare the signed derivative of hardening storage. The reference curve
    /// supplies fixed strain knots; every absolute yield ordinate needs its own
    /// finite derivative. The receiving thermal family admits the temperature
    /// span, original initial-yield coefficients and native history bounds.
    pub fn new(
        reference: TabulatedHardening,
        yield_stress_slopes_pa_k: Vec<f64>,
    ) -> Result<Self, FiniteJ2ThermalError> {
        let points = reference.points();
        if yield_stress_slopes_pa_k.len() != points.len()
            || yield_stress_slopes_pa_k.iter().any(|s| !s.is_finite())
        {
            return Err(FiniteJ2ThermalError::Parameters(
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
                return Err(FiniteJ2ThermalError::Numerical(
                    "tabulated temperature derivative is below representable slope resolution",
                ));
            }
            let start = storage[index];
            let end = finite(start + finite((0.5 * lower + 0.5 * upper) * width)?)?;
            maximum = maximum.max(end.abs());
            // Signed primitives may have an interior extremum even when their
            // knot values are small. Admit that extremum, not only endpoints.
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

    /// Exact immutable curve at the receiving family's reference temperature.
    pub fn reference(&self) -> &TabulatedHardening {
        &self.prepared.reference
    }

    /// Explicit derivatives of the original absolute yield ordinates [Pa/K].
    pub fn yield_stress_slopes_pa_k(&self) -> &[f64] {
        &self.prepared.yield_slopes_pa_k
    }

    /// `partial_T integral_0^alpha (Y(p,T)-Y(0,T)) dp` at fixed history.
    /// The derivative is signed and independent of temperature. Its negative
    /// contributes to mechanical entropy; no material-point allocation occurs.
    pub fn hardening_energy_temperature_derivative_j_m3_k(
        &self,
        alpha: f64,
    ) -> Result<f64, FiniteJ2ThermalError> {
        let data = &self.prepared;
        let points = data.reference.points();
        let maximum = data.reference.maximum_plastic_strain();
        if !alpha.is_finite() || alpha < 0.0 || alpha > maximum {
            return Err(FiniteJ2Error::Tabulated(TabulatedError::Domain {
                plastic_strain: alpha,
                maximum,
            })
            .into());
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

    /// Prepare one endpoint outside constitutive iterations. At the reference
    /// temperature, preserve the original curve's prepared arithmetic exactly.
    pub(super) fn frozen_table(
        &self,
        temperature_change_k: f64,
    ) -> Result<TabulatedHardening, FiniteJ2ThermalError> {
        if temperature_change_k == 0.0 {
            return Ok(self.reference().clone());
        }
        let points: Result<Vec<HardeningPoint>, FiniteJ2ThermalError> = self
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
            .collect();
        TabulatedHardening::new(points?)
            .map_err(|error| FiniteJ2ThermalError::Material(FiniteJ2Error::Tabulated(error)))
    }
}

impl AffineThermalFiniteJ2 {
    /// Select one bounded table whose yield ordinates evolve affinely with
    /// temperature. The table replaces the complete isotropic yield radius:
    /// original H and dH/dT must both be zero, and initial yield/value derivative
    /// must match its first knot. The complete supplied strain/temperature
    /// rectangle is admitted. Different selected curves cannot replace history.
    ///
    /// Perzyna and simultaneous Norton/primary creep retain this selection.
    /// Pure Norton creep has no yield hardening and refuses; Voce and tabulated
    /// hardening are exclusive. This does not infer annealing or texture change.
    pub fn with_tabulated_hardening(
        mut self,
        hardening: FiniteJ2ThermalTabulatedHardening,
    ) -> Result<Self, FiniteJ2ThermalError> {
        let p = self.parameters;
        if self.norton_creep
            || self.voce_hardening.is_some()
            || p.reference.hardening_modulus_pa != 0.0
            || p.slopes.hardening_modulus_pa_k != 0.0
            || p.reference.yield_stress_pa != hardening.reference().initial_yield_stress_pa()
            || p.slopes.yield_stress_pa_k != hardening.yield_stress_slopes_pa_k()[0]
            || self.limits.maximum_accumulated_plastic_strain
                > hardening.reference().maximum_plastic_strain()
            || self
                .tabulated_hardening
                .as_ref()
                .is_some_and(|old| old != &hardening)
        {
            return Err(FiniteJ2ThermalError::Parameters(
                "thermal table needs matching initial yield and derivative, zero linear hardening, full history support, a plastic family and no competing or replacement selection",
            ));
        }
        finite(
            p.temperature_range_k[1] * hardening.prepared.maximum_absolute_storage_slope_j_m3_k,
        )?;
        self.tabulated_hardening = Some(hardening);
        // Affine ordinate values and differences achieve extrema at these
        // endpoints. Native admission also checks every return stiffness.
        for temperature in p.temperature_range_k {
            self.at_temperature_for_inspection(temperature)?;
        }
        Ok(self)
    }

    /// Exact original curve and temperature derivatives, shared by endpoints.
    pub fn tabulated_hardening(&self) -> Option<&FiniteJ2ThermalTabulatedHardening> {
        self.tabulated_hardening.as_ref()
    }
}
