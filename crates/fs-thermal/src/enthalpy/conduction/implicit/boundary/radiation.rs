//! A resolved exterior surface on the original enthalpy boundary channel.
//!
//! The large isothermal surroundings fill the radiative view. Outward surface
//! flux is `h_c (T_s-T_b) + epsilon sigma (T_s^4-T_b^4)`, following
//! <https://mooseframework.inl.gov/docs/site/source/components/HSBoundaryRadiation.html>.
//! The original finite-volume half-cell supplies `k(h) (T_c-T_s)/d`, with an
//! optional series coating/contact resistance. No surface storage is added.

use super::super::ImplicitConductionLaw;
use super::surface::BoundGraySurface;
use super::{
    ImplicitBoundaryGeometry, ImplicitBoundaryLaw, ImplicitBoundaryState, ImplicitGraySurfaceLaw,
    ImplicitGraySurfaceState,
};
use crate::enthalpy::{DepositionError, ThermalLaw, finite, nonnegative, positive};

// CODATA value, rounded consistently with the existing fs-conduction owner.
// https://physics.nist.gov/cgi-bin/cuu/Value?sigma|search_for=stefan+boltzmann
const SIGMA: f64 = 5.670_374_419e-8;

/// Bounded scalar surface-temperature work, separate from thermal time error.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GrayAmbientControls {
    /// At most this many safeguarded Newton surface evaluations, 1..=128.
    pub maximum_iterations: usize,
    /// Absolute allowance for the original surface power balance [W].
    pub absolute_power_tolerance_w: f64,
    /// Relative allowance on the actual surface powers, in [0,1).
    pub relative_power_tolerance: f64,
}
impl Default for GrayAmbientControls {
    fn default() -> Self {
        Self {
            maximum_iterations: 64,
            absolute_power_tolerance_w: 1e-10,
            relative_power_tolerance: 1e-12,
        }
    }
}

/// Complete surface solution at one original caloric endpoint.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GrayAmbientResponse {
    /// Secant and complete dQ_out/dh consumed by the native Newton operator.
    pub boundary: ImplicitBoundaryState,
    /// Exterior temperature after the half-cell and coating/contact drop [K].
    pub surface_temperature_k: f64,
    /// Exactly the power applied by the native secant, positive outward [W].
    pub outward_power_w: f64,
    /// Surface convection contribution, positive outward [W].
    pub convection_power_w: f64,
    /// Surface gray-radiation contribution, positive outward [W].
    pub radiation_power_w: f64,
    /// Applied power minus the sum of the two surface contributions [W].
    pub surface_residual_w: f64,
    /// Actual safeguarded scalar evaluations; zero for insulation/zero area.
    pub iterations: usize,
}

/// Nonlinear convection and radiation on one existing ambient/reservoir face.
///
/// Attach using `with_boundary_laws` on the original chain or graph. Both
/// mechanisms use the supplied bath temperature; distinct gas and radiative
/// bath temperatures require a different boundary model. Emissivity and htc
/// are explicit constants or the original `with_surface_law` response, never
/// inferred from a material name or the cell's authored htc.
/// This law REPLACES the selected channel and must not also enter as a source.
///
/// The original k(h) is cellwise uniform in each trial, as in the native
/// finite-volume conduction law. Its full derivative is retained, including
/// phase dependence on a latent plateau. This does not integrate a varying
/// k(T) profile inside a half-cell or resolve a finite radiation enclosure.
/// The geometry's resistance denotes an additional coating/contact BEFORE
/// the exterior surface; convection htc is supplied separately and counted
/// once. No moving area or surface heat capacity is introduced. The face
/// kernel allocates nothing; a source law retains its own evaluation costs.
/// The geometry-aware trait query includes original density-driven thickness
/// in this same surface solve, including on an isothermal phase plateau.
pub struct GrayAmbientBoundary<'a> {
    material: &'a dyn ImplicitConductionLaw,
    emissivity: f64,
    convection_w_m2_k: f64,
    surface_law: Option<BoundGraySurface<'a>>,
    controls: GrayAmbientControls,
}

impl<'a> GrayAmbientBoundary<'a> {
    /// Bind the original caloric/conductivity law and explicit surface data.
    pub fn new(
        material: &'a dyn ImplicitConductionLaw,
        emissivity: f64,
        convection_w_m2_k: f64,
        controls: GrayAmbientControls,
    ) -> Result<Self, DepositionError> {
        if !emissivity.is_finite()
            || !(0.0..=1.0).contains(&emissivity)
            || !nonnegative(convection_w_m2_k)
            || !(1..=128).contains(&controls.maximum_iterations)
            || !nonnegative(controls.absolute_power_tolerance_w)
            || !nonnegative(controls.relative_power_tolerance)
            || controls.relative_power_tolerance >= 1.0
        {
            return Err(DepositionError::InvalidInput(
                "gray boundary needs emissivity in [0,1], nonnegative htc and bounded surface controls",
            ));
        }
        Ok(Self {
            material,
            emissivity,
            convection_w_m2_k,
            surface_law: None,
            controls,
        })
    }

    /// Replace zero fixed emissivity/htc with one original enthalpy surface law.
    /// The same caloric object and finite source bounds are required. Actual
    /// source support is checked once at each trial's original bulk enthalpy,
    /// before insulation shortcuts. The source owns both exterior coefficients;
    /// a second attachment or nonzero authored coefficients refuse.
    pub fn with_surface_law(
        mut self,
        law: &'a dyn ImplicitGraySurfaceLaw,
    ) -> Result<Self, DepositionError> {
        if self.surface_law.is_some() || self.emissivity != 0.0 || self.convection_w_m2_k != 0.0 {
            return Err(DepositionError::InvalidInput(
                "gray surface source must uniquely replace zero fixed emissivity and htc",
            ));
        }
        self.surface_law = Some(BoundGraySurface::new(self.material.caloric_law(), law)?);
        Ok(self)
    }

    /// Resolve surface temperature, original heat balance and analytic tangent.
    /// A failed source, surface solve or checkpoint publishes no response.
    pub fn evaluate(
        &self,
        specific_h_j_kg: f64,
        bath_temperature_k: f64,
        geometry: ImplicitBoundaryGeometry,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<GrayAmbientResponse, DepositionError> {
        self.evaluate_with_geometry_tangent(
            specific_h_j_kg,
            bath_temperature_k,
            geometry,
            0.0,
            checkpoint,
        )
    }

    fn evaluate_with_geometry_tangent(
        &self,
        specific_h_j_kg: f64,
        bath_temperature_k: f64,
        geometry: ImplicitBoundaryGeometry,
        relative_distance_slope_kg_j: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<GrayAmbientResponse, DepositionError> {
        checkpoint()?;
        if !positive(bath_temperature_k)
            || !nonnegative(geometry.area_m2)
            || !positive(geometry.half_distance_m)
            || !nonnegative(geometry.surface_resistance_m2_k_w)
            || !relative_distance_slope_kg_j.is_finite()
        {
            return Err(DepositionError::InvalidInput(
                "invalid gray boundary temperature or geometry",
            ));
        }
        let state = self.material.state_at_specific_enthalpy(specific_h_j_kg)?;
        let t = state.temperature_k;
        let k = state.conductivity_w_m_k;
        if !positive(t)
            || !nonnegative(state.temperature_slope_kg_k_j)
            || !nonnegative(k)
            || !state.conductivity_slope_kg_w_m_k_j.is_finite()
            || (k == 0.0 && state.conductivity_slope_kg_w_m_k_j != 0.0)
        {
            return Err(DepositionError::MaterialDomain);
        }
        let coefficients = match &self.surface_law {
            Some(law) => law.evaluate(specific_h_j_kg, checkpoint)?,
            None => ImplicitGraySurfaceState::fixed(self.emissivity, self.convection_w_m2_k),
        };
        if geometry.area_m2 == 0.0
            || k == 0.0
            || (coefficients.emissivity == 0.0
                && coefficients.convection_w_m2_k == 0.0
                && !coefficients.has_enthalpy_response())
        {
            checkpoint()?;
            return Ok(GrayAmbientResponse {
                boundary: ImplicitBoundaryState {
                    conductance_w_k: 0.0,
                    power_enthalpy_derivative_w_kg_j: 0.0,
                },
                surface_temperature_k: if k == 0.0 { bath_temperature_k } else { t },
                outward_power_w: 0.0,
                convection_power_w: 0.0,
                radiation_power_w: 0.0,
                surface_residual_w: 0.0,
                iterations: 0,
            });
        }
        let bulk_resistance = finite(geometry.half_distance_m / k)?;
        let resistance = finite(bulk_resistance + geometry.surface_resistance_m2_k_w)?;
        let bulk = finite(resistance.recip())?;
        if !positive(bulk_resistance) || !positive(bulk) {
            return Err(DepositionError::Arithmetic);
        }
        let delta = t - bath_temperature_k;
        // x=(Ts-Tb)/(Tc-Tb) preserves a small thermal difference without
        // subtracting rounded surface temperatures or fourth powers.
        let mut x = fraction(
            bulk,
            self.surface_coefficients(t, bath_temperature_k, coefficients)?
                .0,
        );
        let (mut lo, mut hi) = (0.0, 1.0);
        for iterations in 1..=self.controls.maximum_iterations {
            checkpoint()?;
            let surface = delta.mul_add(x, bath_temperature_k);
            let (surface_g, surface_slope, radiation_g) =
                self.surface_coefficients(surface, bath_temperature_k, coefficients)?;
            let g = finite(geometry.area_m2 * series(bulk, surface_g))?;
            let power = finite(g * delta)?;
            let surface_delta = finite(x * delta)?;
            let convection =
                finite(geometry.area_m2 * coefficients.convection_w_m2_k * surface_delta)?;
            let radiation = finite(geometry.area_m2 * radiation_g * surface_delta)?;
            let surface_power = finite(convection + radiation)?;
            let conducted = finite(geometry.area_m2 * bulk * ((1.0 - x) * delta))?;
            let residual = finite(power - surface_power)?;
            let allowance = finite(
                self.controls.absolute_power_tolerance_w
                    + self.controls.relative_power_tolerance
                        * power.abs().max(surface_power.abs()).max(conducted.abs()),
            )?;
            if residual.abs().max((conducted - surface_power).abs()) <= allowance {
                // dQ/dh = A F'(Ts)/(a+F'(Ts)) * [a T'(h)+a'(h)(Tc-Ts)].
                // a'=a*(bulk R / total R)*(k'/k-d'/d). The original
                // density can change thickness even on a latent plateau.
                let bulk_derivative = finite(
                    bulk * (bulk_resistance / resistance)
                        * (state.conductivity_slope_kg_w_m_k_j / k
                            - relative_distance_slope_kg_j),
                )?;
                let mut derivative = finite(
                    geometry.area_m2
                        * fraction(surface_slope, bulk)
                        * finite(
                            bulk * state.temperature_slope_kg_k_j
                                + bulk_derivative * ((1.0 - x) * delta),
                        )?,
                )?;
                if coefficients.has_enthalpy_response() {
                    // F(Ts,h) now varies directly with original h. From
                    // a(T-Ts)=F, Qh=A*[Fs*(a*Th+ah*(T-Ts))+a*Fh]/(a+Fs).
                    // Keep Fh at fixed Ts, including emissivity onset at zero.
                    let source_derivative = coefficients.flux_enthalpy_derivative(
                        surface,
                        bath_temperature_k,
                        surface_delta,
                        surface_delta,
                    )?;
                    derivative = finite(
                        derivative
                            + geometry.area_m2 * fraction(bulk, surface_slope) * source_derivative,
                    )?;
                }
                checkpoint()?;
                return Ok(GrayAmbientResponse {
                    boundary: ImplicitBoundaryState {
                        conductance_w_k: g,
                        power_enthalpy_derivative_w_kg_j: derivative,
                    },
                    surface_temperature_k: surface,
                    outward_power_w: power,
                    convection_power_w: convection,
                    radiation_power_w: radiation,
                    surface_residual_w: residual,
                    iterations,
                });
            }
            let defect = finite(bulk * (1.0 - x) - surface_g * x)?;
            if defect > 0.0 {
                lo = x;
            } else {
                hi = x;
            }
            let proposed = x + defect / finite(bulk + surface_slope)?;
            let next = if proposed > lo && proposed < hi && proposed != x {
                proposed
            } else {
                lo + 0.5 * (hi - lo)
            };
            if next == x || next <= lo || next >= hi {
                break;
            }
            x = next;
        }
        Err(DepositionError::InvalidInput(
            "gray boundary surface balance did not converge",
        ))
    }

    // Factoring T^4-Tb^4 retains near-equilibrium radiative flux exactly once.
    fn surface_coefficients(
        &self,
        t: f64,
        bath: f64,
        coefficients: ImplicitGraySurfaceState,
    ) -> Result<(f64, f64, f64), DepositionError> {
        let (radiation, slope) = if coefficients.emissivity == 0.0 {
            (0.0, 0.0)
        } else {
            let scale = coefficients.emissivity * SIGMA;
            (
                finite(scale * (t + bath) * t.mul_add(t, bath * bath))?,
                finite(((4.0 * scale * t) * t) * t)?,
            )
        };
        Ok((
            finite(coefficients.convection_w_m2_k + radiation)?,
            finite(coefficients.convection_w_m2_k + slope)?,
            radiation,
        ))
    }
}

impl ImplicitBoundaryLaw for GrayAmbientBoundary<'_> {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self.material.caloric_law()
    }

    fn state_at_specific_enthalpy(
        &self,
        specific_h_j_kg: f64,
        bath_temperature_k: f64,
        geometry: ImplicitBoundaryGeometry,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ImplicitBoundaryState, DepositionError> {
        self.evaluate(specific_h_j_kg, bath_temperature_k, geometry, checkpoint)
            .map(|response| response.boundary)
    }

    fn state_at_specific_enthalpy_with_geometry_tangent(
        &self,
        specific_h_j_kg: f64,
        bath_temperature_k: f64,
        geometry: ImplicitBoundaryGeometry,
        half_distance_relative_enthalpy_derivative_kg_j: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Option<ImplicitBoundaryState>, DepositionError> {
        self.evaluate_with_geometry_tangent(
            specific_h_j_kg,
            bath_temperature_k,
            geometry,
            half_distance_relative_enthalpy_derivative_kg_j,
            checkpoint,
        )
        .map(|response| Some(response.boundary))
    }
}

// Positive series coefficients, avoiding a*b and an overflowing a+b.
fn series(a: f64, b: f64) -> f64 {
    if a <= b {
        a / (1.0 + a / b)
    } else {
        b / (1.0 + b / a)
    }
}
fn fraction(a: f64, b: f64) -> f64 {
    if a <= b {
        let ratio = a / b;
        ratio / (1.0 + ratio)
    } else {
        1.0 / (1.0 + b / a)
    }
}
