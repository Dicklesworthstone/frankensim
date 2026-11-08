//! Convection and gray radiation to distinct baths through one exterior surface.
//!
//! The two flux laws follow MOOSE's `HSBoundaryAmbientConvection` and
//! `HSBoundaryRadiation` formulations. Their prescribed temperatures need not
//! agree; the original finite half-cell is shared, never duplicated.

use super::super::ImplicitConductionLaw;
use super::surface::BoundGraySurface;
use super::{
    BathTemperatureSchedule, GrayAmbientBoundary, GrayAmbientControls, ImplicitBoundaryExchangeLaw,
    ImplicitBoundaryExchangeState, ImplicitBoundaryGeometry, ImplicitGraySurfaceLaw,
    ImplicitGraySurfaceState,
};
use crate::enthalpy::backward_euler::PrescribedBathHeat;
use crate::enthalpy::{DepositionError, ThermalLaw, finite, nonnegative, positive};

// The same CODATA constant as the original gray ambient boundary.
const SIGMA: f64 = 5.670_374_419e-8;

/// A resolved shared surface and the two exchanges actually used by the owner.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GrayTwoBathResponse {
    /// Bath 0 is convection; bath 1 is gray radiation. Both carry actual Kelvin.
    pub boundary: ImplicitBoundaryExchangeState,
    pub surface_temperature_k: f64,
    /// Sum of the two signed powers, positive outward from the original cell.
    pub outward_power_w: f64,
    pub convection_power_w: f64,
    pub radiation_power_w: f64,
    /// Conducted power minus the sum of the two exterior powers [W].
    pub surface_residual_w: f64,
    pub iterations: usize,
}

/// One half-cell and surface exchanging heat with distinct gas and radiative baths.
///
/// The original selected channel supplies the gas temperature, including its
/// scheduled boundary events. The large black radiative surrounding has the
/// explicitly fixed temperature supplied here, or the continuous history
/// attached with `with_radiative_temperature_history`. Attach through the original
/// chain/graph `with_boundary_exchanges`; each exchange then retains its actual
/// reservoir entropy rather than using an inferred effective temperature.
///
/// Emissivity and convection htc are fixed or supplied by the original
/// `with_surface_law` bulk-enthalpy response. The geometry's resistance denotes
/// coating/contact BEFORE this common exterior surface; it is not a second
/// convection film. The original k(h) is cellwise uniform in a trial, with its
/// full derivative retained even on a latent plateau. The surface has no
/// storage, so the original caloric object remains the only body-energy owner.
/// This is a gray full-view boundary, not a finite radiosity enclosure or a
/// spectral photon-entropy model. The surface kernel allocates no memory; a
/// source law retains its own evaluation costs. Gas boundary events obey native
/// time integrator's admission rules. The geometry-aware trait query admits
/// original density-driven half-thickness while retaining both actual bath
/// powers, their distinct temperatures and the same shared surface balance.
pub struct GrayTwoBathBoundary<'a> {
    material: &'a dyn ImplicitConductionLaw,
    emissivity: f64,
    convection_w_m2_k: f64,
    surface_law: Option<BoundGraySurface<'a>>,
    radiative_temperature_k: f64,
    radiative_history: Option<&'a BathTemperatureSchedule>,
    controls: GrayAmbientControls,
}

impl<'a> GrayTwoBathBoundary<'a> {
    pub fn new(
        material: &'a dyn ImplicitConductionLaw,
        emissivity: f64,
        convection_w_m2_k: f64,
        radiative_temperature_k: f64,
        controls: GrayAmbientControls,
    ) -> Result<Self, DepositionError> {
        if !positive(radiative_temperature_k) {
            return Err(DepositionError::InvalidInput(
                "gray radiative bath temperature must be positive and finite",
            ));
        }
        // Keep the existing emissivity, htc and numerical-control admission.
        GrayAmbientBoundary::new(material, emissivity, convection_w_m2_k, controls)?;
        Ok(Self {
            material,
            emissivity,
            convection_w_m2_k,
            surface_law: None,
            radiative_temperature_k,
            radiative_history: None,
            controls,
        })
    }

    /// Replace zero fixed emissivity/htc with one original enthalpy surface law.
    /// Its caloric identity, bounds and source coordinates have the same
    /// meaning as `GrayAmbientBoundary::with_surface_law`. Each trial resolves
    /// the source once, then both baths exchange heat at the same solved surface
    /// and retain their actual temperatures for entropy accounting.
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

    /// Replace the fixed radiative bath with one prescribed continuous history.
    ///
    /// The original native time owner samples this history at each physical
    /// endpoint and SDIRK stage, independently of the gas bath. Its intervals
    /// stop at every knot. Both bath powers retain their actual temperatures
    /// for heat and entropy accounting, with one shared half-cell resistance.
    /// The closed source time range never extrapolates, including when the
    /// emissivity is zero. This borrows the original immutable schedule and
    /// introduces neither surface storage nor another heat source.
    pub fn with_radiative_temperature_history(
        mut self,
        history: &'a BathTemperatureSchedule,
    ) -> Result<Self, DepositionError> {
        if self.radiative_history.is_some() {
            return Err(DepositionError::InvalidInput(
                "radiative bath already has a prescribed temperature history",
            ));
        }
        self.radiative_history = Some(history);
        Ok(self)
    }

    /// Resolve one shared surface, its signed bath powers and exact total tangent.
    /// Failed material queries, surface convergence or cancellation return no
    /// response. Numerical acceptance controls the NET surface balance, even
    /// when large opposing bath powers nearly cancel.
    pub fn evaluate(
        &self,
        specific_h_j_kg: f64,
        gas_temperature_k: f64,
        geometry: ImplicitBoundaryGeometry,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<GrayTwoBathResponse, DepositionError> {
        self.evaluate_with_geometry_tangent(
            specific_h_j_kg,
            gas_temperature_k,
            geometry,
            0.0,
            None,
            checkpoint,
        )
    }

    /// Resolve both actual bath powers at a physical time without advancing it.
    /// A history-owned boundary requires this query (or the native time-aware
    /// trait query); the timeless query refuses instead of using a stale bath.
    pub fn evaluate_at_time(
        &self,
        specific_h_j_kg: f64,
        gas_temperature_k: f64,
        geometry: ImplicitBoundaryGeometry,
        time_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<GrayTwoBathResponse, DepositionError> {
        self.evaluate_with_geometry_tangent(
            specific_h_j_kg,
            gas_temperature_k,
            geometry,
            0.0,
            Some(time_s),
            checkpoint,
        )
    }

    fn evaluate_with_geometry_tangent(
        &self,
        specific_h_j_kg: f64,
        gas_temperature_k: f64,
        geometry: ImplicitBoundaryGeometry,
        relative_distance_slope_kg_j: f64,
        time_s: Option<f64>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<GrayTwoBathResponse, DepositionError> {
        checkpoint()?;
        if time_s.is_some_and(|time| !time.is_finite() || time < 0.0) {
            return Err(DepositionError::InvalidInput(
                "invalid boundary exchange time",
            ));
        }
        let radiative_temperature_k = match (self.radiative_history, time_s) {
            (Some(history), Some(time)) => history.temperature_at(time)?,
            (Some(_), None) => {
                return Err(DepositionError::InvalidInput(
                    "scheduled radiative bath requires physical time",
                ));
            }
            (None, _) => self.radiative_temperature_k,
        };
        if !positive(gas_temperature_k)
            || !nonnegative(geometry.area_m2)
            || !positive(geometry.half_distance_m)
            || !nonnegative(geometry.surface_resistance_m2_k_w)
            || !relative_distance_slope_kg_j.is_finite()
        {
            return Err(DepositionError::InvalidInput(
                "invalid two-bath gray boundary temperature or geometry",
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
            || (coefficients.emissivity == 0.0
                && coefficients.convection_w_m2_k == 0.0
                && (k == 0.0 || !coefficients.has_enthalpy_response()))
        {
            checkpoint()?;
            return Ok(self.response(
                t, gas_temperature_k, radiative_temperature_k, [0.0; 2], 0.0, 0.0, 0,
            )?);
        }
        let (bulk, conductivity_weight) = if k == 0.0 {
            // The massless exterior can still exchange heat between its two
            // baths, while its disconnected bulk receives no resolved heat.
            (0.0, 0.0)
        } else {
            let bulk_resistance = finite(geometry.half_distance_m / k)?;
            let resistance = finite(bulk_resistance + geometry.surface_resistance_m2_k_w)?;
            let bulk = finite(resistance.recip())?;
            if !positive(bulk_resistance) || !positive(bulk) {
                return Err(DepositionError::Arithmetic);
            }
            (bulk, bulk_resistance / resistance)
        };
        // Work in Ts-Tgas rather than subtracting rounded surface Kelvin. This
        // retains small near-equilibrium fluxes, including the one-bath limit.
        let cell_delta = t - gas_temperature_k;
        let radiative_delta = radiative_temperature_k - gas_temperature_k;
        let mut lo = 0.0_f64.min(cell_delta).min(radiative_delta);
        let mut hi = 0.0_f64.max(cell_delta).max(radiative_delta);
        let initial_radiation =
            Self::radiation_coefficients(t, radiative_temperature_k, coefficients.emissivity)?.0;
        let initial_total = finite(bulk + coefficients.convection_w_m2_k + initial_radiation)?;
        if !positive(initial_total) {
            return Err(DepositionError::Arithmetic);
        }
        let guess = finite(
            (bulk / initial_total) * cell_delta
                + (initial_radiation / initial_total) * radiative_delta,
        )?;
        let mut offset = if guess >= lo && guess <= hi {
            guess
        } else {
            lo + 0.5 * (hi - lo)
        };
        for iterations in 1..=self.controls.maximum_iterations {
            checkpoint()?;
            let surface = finite(gas_temperature_k + offset)?;
            if !positive(surface) {
                return Err(DepositionError::Arithmetic);
            }
            let (radiation_g, radiation_slope) =
                Self::radiation_coefficients(surface, radiative_temperature_k, coefficients.emissivity)?;
            let surface_slope = finite(coefficients.convection_w_m2_k + radiation_slope)?;
            let convection_flux = finite(coefficients.convection_w_m2_k * offset)?;
            let radiation_flux = finite(radiation_g * (offset - radiative_delta))?;
            let conducted_flux = finite(bulk * (cell_delta - offset))?;
            let convection = finite(geometry.area_m2 * convection_flux)?;
            let radiation = finite(geometry.area_m2 * radiation_flux)?;
            let outward = finite(convection + radiation)?;
            let conducted = finite(geometry.area_m2 * conducted_flux)?;
            let residual = finite(conducted - outward)?;
            let allowance = finite(
                self.controls.absolute_power_tolerance_w
                    + self.controls.relative_power_tolerance * conducted.abs().max(outward.abs()),
            )?;
            if residual.abs() <= allowance {
                let bulk_derivative = if k == 0.0 {
                    0.0
                } else {
                    finite(bulk * conductivity_weight
                        * (state.conductivity_slope_kg_w_m_k_j / k
                            - relative_distance_slope_kg_j))?
                };
                // Qh=A*f'(Ts)/(a+f'(Ts))*[a*Th+a_h*(Tc-Ts)]. Both bath
                // temperatures are held fixed; a_h includes the original
                // k'/k-d'/d response without an apparent heat capacity.
                let mut derivative = finite(
                    geometry.area_m2
                        * fraction(surface_slope, bulk)
                        * finite(
                            bulk * state.temperature_slope_kg_k_j
                                + bulk_derivative * (cell_delta - offset),
                        )?,
                )?;
                if bulk != 0.0 && coefficients.has_enthalpy_response() {
                    // Direct coefficient response at fixed Ts contributes
                    // A*a/(a+Fs)*Fh after the SAME surface balance is resolved.
                    // A disconnected bulk has exactly zero total derivative,
                    // while its two baths may still exchange finite heat.
                    let source_derivative = coefficients.flux_enthalpy_derivative(
                        surface,
                        radiative_temperature_k,
                        offset,
                        offset - radiative_delta,
                    )?;
                    derivative = finite(
                        derivative
                            + geometry.area_m2 * fraction(bulk, surface_slope) * source_derivative,
                    )?;
                }
                checkpoint()?;
                return self.response(
                    surface,
                    gas_temperature_k,
                    radiative_temperature_k,
                    [convection, radiation],
                    derivative,
                    residual,
                    iterations,
                );
            }
            let defect = finite(conducted_flux - finite(convection_flux + radiation_flux)?)?;
            if defect > 0.0 {
                lo = offset;
            } else {
                hi = offset;
            }
            let proposed = offset + defect / finite(bulk + surface_slope)?;
            let next = if proposed > lo && proposed < hi && proposed != offset {
                proposed
            } else {
                lo + 0.5 * (hi - lo)
            };
            if next == offset || next <= lo || next >= hi {
                break;
            }
            offset = next;
        }
        Err(DepositionError::InvalidInput(
            "two-bath gray boundary surface balance did not converge",
        ))
    }

    fn radiation_coefficients(
        t: f64,
        bath: f64,
        emissivity: f64,
    ) -> Result<(f64, f64), DepositionError> {
        if emissivity == 0.0 {
            return Ok((0.0, 0.0));
        }
        let scale = emissivity * SIGMA;
        Ok((
            finite(scale * (t + bath) * t.mul_add(t, bath * bath))?,
            finite(((4.0 * scale * t) * t) * t)?,
        ))
    }

    fn response(
        &self,
        surface_temperature_k: f64,
        gas_temperature_k: f64,
        radiative_temperature_k: f64,
        powers: [f64; 2],
        derivative: f64,
        surface_residual_w: f64,
        iterations: usize,
    ) -> Result<GrayTwoBathResponse, DepositionError> {
        let boundary = ImplicitBoundaryExchangeState {
            baths: [
                PrescribedBathHeat {
                    outward_power_w: powers[0],
                    temperature_k: gas_temperature_k,
                },
                PrescribedBathHeat {
                    outward_power_w: powers[1],
                    temperature_k: radiative_temperature_k,
                },
            ],
            power_enthalpy_derivative_w_kg_j: derivative,
        }
        .validate()?;
        Ok(GrayTwoBathResponse {
            boundary,
            surface_temperature_k,
            outward_power_w: boundary.outward_power_w()?,
            convection_power_w: powers[0],
            radiation_power_w: powers[1],
            surface_residual_w,
            iterations,
        })
    }
}

impl ImplicitBoundaryExchangeLaw for GrayTwoBathBoundary<'_> {
    fn supports_current_geometry(&self) -> bool {
        true
    }

    fn caloric_law(&self) -> &dyn ThermalLaw {
        self.material.caloric_law()
    }

    fn state_at_specific_enthalpy(
        &self,
        specific_h_j_kg: f64,
        bath_temperature_k: f64,
        geometry: ImplicitBoundaryGeometry,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ImplicitBoundaryExchangeState, DepositionError> {
        self.evaluate(specific_h_j_kg, bath_temperature_k, geometry, checkpoint)
            .map(|response| response.boundary)
    }

    fn state_at_specific_enthalpy_and_time(
        &self,
        specific_h_j_kg: f64,
        bath_temperature_k: f64,
        geometry: ImplicitBoundaryGeometry,
        time_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ImplicitBoundaryExchangeState, DepositionError> {
        self.evaluate_at_time(
            specific_h_j_kg,
            bath_temperature_k,
            geometry,
            time_s,
            checkpoint,
        )
        .map(|response| response.boundary)
    }

    fn state_at_specific_enthalpy_with_geometry_tangent(
        &self,
        specific_h_j_kg: f64,
        bath_temperature_k: f64,
        geometry: ImplicitBoundaryGeometry,
        half_distance_relative_enthalpy_derivative_kg_j: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Option<ImplicitBoundaryExchangeState>, DepositionError> {
        self.evaluate_with_geometry_tangent(
            specific_h_j_kg,
            bath_temperature_k,
            geometry,
            half_distance_relative_enthalpy_derivative_kg_j,
            None,
            checkpoint,
        )
        .map(|response| Some(response.boundary))
    }

    fn state_at_specific_enthalpy_with_geometry_tangent_and_time(
        &self,
        specific_h_j_kg: f64,
        bath_temperature_k: f64,
        geometry: ImplicitBoundaryGeometry,
        half_distance_relative_enthalpy_derivative_kg_j: f64,
        time_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Option<ImplicitBoundaryExchangeState>, DepositionError> {
        self.evaluate_with_geometry_tangent(
            specific_h_j_kg,
            bath_temperature_k,
            geometry,
            half_distance_relative_enthalpy_derivative_kg_j,
            Some(time_s),
            checkpoint,
        )
        .map(|response| Some(response.boundary))
    }

    fn next_transition_time_s(&self, after_s: f64) -> Option<f64> {
        self.radiative_history
            .and_then(|history| history.next_transition_time_s(after_s))
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
