//! Source-integrated Cp/Cv storage exchanging with two original exterior baths.
//!
//! The existing material-surface kernel closes the exact conductivity integral
//! against convection and gray radiation. This adapter supplies the complete
//! original-energy tangent and keeps the two actual bath powers for the native
//! heat and entropy owners. It introduces neither surface storage nor a clock.

use super::CaloricConductionLayer;
use crate::thermal::conduction::{
    ConvectiveLayerError, LayerConductionConfig, ReservoirLayerGeometry,
};
use fs_thermal::enthalpy::backward_euler::PrescribedBathHeat;
use fs_thermal::enthalpy::conduction::implicit::boundary::{
    GrayTwoBathResponse, ImplicitBoundaryExchangeLaw, ImplicitBoundaryExchangeState,
    ImplicitBoundaryGeometry, ImplicitGraySurfaceLaw, ImplicitGraySurfaceState,
};
use fs_thermal::enthalpy::{DepositionError, ThermalLaw};

// The same rounded CODATA value as the original gray boundary owners.
const SIGMA: f64 = 5.670_374_419e-8;

struct BoundSurface<'a> {
    law: &'a dyn ImplicitGraySurfaceLaw,
    bounds_j_kg: [f64; 2],
}

#[derive(Clone, Copy)]
struct SurfaceTangent {
    temperature_derivative_w_k: f64,
    energy_derivative_w_kg_j: f64,
}

/// One source-integrated material half-cell with convection and gray radiation.
///
/// Attach with the original implicit chain/graph's `with_boundary_exchanges`.
/// The selected channel supplies the gas temperature, including its scheduled
/// changes. The separate large black radiative surrounding has the explicitly
/// fixed temperature supplied here. Each bath retains its own signed power and
/// actual temperature even when the original material receives no net heat.
///
/// The layer is planar, homogeneous, source-free and quasi-steady. Its original
/// Cp or Cv law remains the sole finite-body heat inventory. Both bath laws act
/// on the SAME solved surface; material k(T) is integrated to that surface.
/// Only the material path requires conductivity support, not the bath states.
///
/// This is a full-view gray boundary without coating/contact resistance,
/// participating radiation or a finite enclosure. A nonzero authored surface
/// resistance refuses rather than becoming an inferred second convective film.
/// The kernel allocates no memory; an optional original source retains its own
/// evaluation costs. Spatial and temporal accuracy remain the native solver's
/// responsibility rather than a consequence of this constitutive closure.
pub struct IntegratedCaloricGrayBoundary<'l, 'a, 'm> {
    layer: &'l CaloricConductionLayer<'a, 'm>,
    geometry: ImplicitBoundaryGeometry,
    flux_geometry: ReservoirLayerGeometry,
    material_multiplier_m: f64,
    fixed: ImplicitGraySurfaceState,
    surface: Option<BoundSurface<'l>>,
    radiative_temperature_k: f64,
    config: LayerConductionConfig,
}

impl<'l, 'a, 'm> IntegratedCaloricGrayBoundary<'l, 'a, 'm> {
    /// Bind the original positive area and centre-to-face distance, with zero
    /// coating/contact resistance. HTC owns only the exterior convective film.
    /// Zero emissivity or HTC explicitly disables that exchange. To attach an
    /// original surface source, initially declare BOTH fixed coefficients zero.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        layer: &'l CaloricConductionLayer<'a, 'm>,
        geometry: ImplicitBoundaryGeometry,
        emissivity: f64,
        convection_w_m2_k: f64,
        radiative_temperature_k: f64,
        config: LayerConductionConfig,
    ) -> Result<Self, DepositionError> {
        let flux_geometry =
            ReservoirLayerGeometry::planar(geometry.area_m2, geometry.half_distance_m)?;
        let fixed = ImplicitGraySurfaceState {
            emissivity,
            emissivity_enthalpy_derivative_kg_j: 0.0,
            convection_w_m2_k,
            convection_enthalpy_derivative_kg_w_m2_k_j: 0.0,
        };
        if fixed.admit().is_err()
            || !radiative_temperature_k.is_finite()
            || radiative_temperature_k <= 0.0
            || geometry.surface_resistance_m2_k_w != 0.0
            || config.max_iterations == 0
            || !config.absolute_heat_tolerance_w.is_finite()
            || config.absolute_heat_tolerance_w < 0.0
            || !config.relative_heat_tolerance.is_finite()
            || config.relative_heat_tolerance < 0.0
        {
            return Err(DepositionError::InvalidInput(
                "integrated gray boundary needs passive coefficients, a positive bath, zero coating resistance and bounded solve controls",
            ));
        }
        Ok(Self {
            layer,
            geometry,
            flux_geometry,
            material_multiplier_m: geometry.area_m2 / geometry.half_distance_m,
            fixed,
            surface: None,
            radiative_temperature_k,
            config,
        })
    }

    /// Uniquely replace zero fixed coefficients with one original bulk-energy
    /// surface source. Its supported h/u bounds must lie inside the original
    /// caloric/conductivity support and retain the exact same caloric object.
    ///
    /// This preserves `ImplicitGraySurfaceLaw`'s BULK-state closure: the source
    /// is evaluated once per energy trial, then its coefficients are held fixed
    /// during the surface solve. Its complete original energy derivatives are
    /// included afterward. Sources requiring the solved exterior temperature,
    /// evolving surface phases or changing flow states need their own closure;
    /// the adapter does not substitute those into the source's coordinates.
    pub fn with_surface_law(
        mut self,
        law: &'l dyn ImplicitGraySurfaceLaw,
    ) -> Result<Self, DepositionError> {
        let [lo, hi] = law.enthalpy_bounds_j_kg();
        let [source_lo, source_hi] = self.layer.energy_bounds_j_kg();
        if self.surface.is_some()
            || self.fixed.emissivity != 0.0
            || self.fixed.convection_w_m2_k != 0.0
            || !lo.is_finite()
            || !hi.is_finite()
            || lo >= hi
            || lo < source_lo
            || hi > source_hi
            || !core::ptr::addr_eq(
                self.layer.caloric_law() as *const dyn ThermalLaw,
                law.caloric_law() as *const dyn ThermalLaw,
            )
        {
            return Err(DepositionError::InvalidInput(
                "integrated gray surface source must uniquely replace zero fixed coefficients on its original supported caloric object",
            ));
        }
        self.surface = Some(BoundSurface {
            law,
            bounds_j_kg: [lo, hi],
        });
        Ok(self)
    }

    /// Resolve one material surface, both bath powers and their total tangent
    /// with respect to the original h or u. Neither physical history nor source
    /// state changes. The returned surface residual follows the existing
    /// integral kernel's declared power tolerance, including opposing drives.
    pub fn evaluate(
        &self,
        specific_energy_j_kg: f64,
        gas_temperature_k: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<GrayTwoBathResponse, DepositionError> {
        checkpoint()?;
        let material = self.layer.state_at_specific_energy(specific_energy_j_kg)?;
        let coefficients = if let Some(surface) = &self.surface {
            let [lo, hi] = surface.bounds_j_kg;
            if !(lo..=hi).contains(&specific_energy_j_kg) {
                return Err(DepositionError::MaterialDomain);
            }
            checkpoint()?;
            let state = surface
                .law
                .state_at_specific_enthalpy(specific_energy_j_kg)?;
            state.admit()?;
            checkpoint()?;
            state
        } else {
            self.fixed
        };
        let curve = self.layer.conductivity_curve();
        let response = self
            .flux_geometry
            .evaluate_material_to_two_baths(
                curve,
                material.temperature_k,
                [gas_temperature_k, self.radiative_temperature_k],
                curve.temperature_range_k(),
                self.config,
                |surface| self.exterior_exchange(surface, gas_temperature_k, coefficients),
                &mut *checkpoint,
            )
            .map_err(|error| match error {
                ConvectiveLayerError::Thermal(error) | ConvectiveLayerError::Film(error) => error,
                ConvectiveLayerError::SurfaceBalance { .. } => DepositionError::InvalidInput(
                    "integrated caloric gray surface balance did not converge",
                ),
            })?;
        let cell_slope = finite(self.material_multiplier_m * material.conductivity_w_m_k)?;
        let surface_slope = finite(
            self.material_multiplier_m
                * curve.conductivity_w_m_k(response.surface_temperature_k)?,
        )?;
        if cell_slope <= 0.0 || surface_slope <= 0.0 {
            return Err(DepositionError::Arithmetic);
        }
        let exterior_slope = response.film.temperature_derivative_w_k;
        // Differentiate M*(K(Tcell)-K(Ts)) = F(Ts,e). The integral's two
        // endpoint conductivities already include the entire k(T) response;
        // adding a bulk dk/de term would double count it. Source F_e is taken
        // at fixed Ts and may be signed, including zero-coefficient onset.
        let storage = (fraction(exterior_slope, surface_slope) * cell_slope)
            * material.temperature_slope_kg_k_j;
        let derivative = finite(
            fraction(surface_slope, exterior_slope)
                .mul_add(response.film.energy_derivative_w_kg_j, storage),
        )?;
        let [convection, radiation] = response.powers_to_baths_w;
        let outward = finite(convection + radiation)?;
        let boundary = ImplicitBoundaryExchangeState {
            baths: [
                PrescribedBathHeat {
                    outward_power_w: convection,
                    temperature_k: gas_temperature_k,
                },
                PrescribedBathHeat {
                    outward_power_w: radiation,
                    temperature_k: self.radiative_temperature_k,
                },
            ],
            power_enthalpy_derivative_w_kg_j: derivative,
        };
        checkpoint()?;
        Ok(GrayTwoBathResponse {
            boundary,
            surface_temperature_k: response.surface_temperature_k,
            outward_power_w: outward,
            convection_power_w: convection,
            radiation_power_w: radiation,
            surface_residual_w: response.heat_residual_w,
            iterations: response.iterations,
        })
    }

    fn exterior_exchange(
        &self,
        surface: f64,
        gas: f64,
        source: ImplicitGraySurfaceState,
    ) -> Result<([f64; 2], [f64; 2], SurfaceTangent), DepositionError> {
        let area = self.geometry.area_m2;
        let gas_delta = surface - gas;
        let radiation_delta = surface - self.radiative_temperature_k;
        let convection_g = finite(area * source.convection_w_m2_k)?;
        let mut radiation_g = 0.0;
        let mut radiation_slope = 0.0;
        let mut radiation_energy = 0.0;
        if source.emissivity != 0.0 || source.emissivity_enthalpy_derivative_kg_j != 0.0 {
            let sum = finite(surface + self.radiative_temperature_k)?;
            let squares = finite(surface.mul_add(
                surface,
                self.radiative_temperature_k * self.radiative_temperature_k,
            ))?;
            if source.emissivity != 0.0 {
                let scale = finite((area * SIGMA) * source.emissivity)?;
                radiation_g = finite((scale * sum) * squares)?;
                radiation_slope = finite(((4.0 * scale * surface) * surface) * surface)?;
                if radiation_g <= 0.0 || radiation_slope <= 0.0 {
                    return Err(DepositionError::Arithmetic);
                }
            }
            if source.emissivity_enthalpy_derivative_kg_j != 0.0 {
                // No division by emissivity: a supported source branch can
                // start at epsilon=0 with a nonzero original energy tangent.
                let scale = finite((area * SIGMA) * source.emissivity_enthalpy_derivative_kg_j)?;
                radiation_energy = finite(((scale * sum) * squares) * radiation_delta)?;
            }
        }
        if source.convection_w_m2_k > 0.0 && convection_g == 0.0 {
            return Err(DepositionError::Arithmetic);
        }
        let convection = finite(convection_g * gas_delta)?;
        let radiation = finite(radiation_g * radiation_delta)?;
        let convection_energy =
            finite((area * source.convection_enthalpy_derivative_kg_w_m2_k_j) * gas_delta)?;
        Ok((
            [convection, radiation],
            [convection_g, radiation_g],
            SurfaceTangent {
                temperature_derivative_w_k: finite(convection_g + radiation_slope)?,
                energy_derivative_w_kg_j: finite(convection_energy + radiation_energy)?,
            },
        ))
    }
}

impl ImplicitBoundaryExchangeLaw for IntegratedCaloricGrayBoundary<'_, '_, '_> {
    fn caloric_law(&self) -> &dyn ThermalLaw {
        self.layer.caloric_law()
    }

    fn state_at_specific_enthalpy(
        &self,
        specific_h_j_kg: f64,
        bath_temperature_k: f64,
        geometry: ImplicitBoundaryGeometry,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<ImplicitBoundaryExchangeState, DepositionError> {
        checkpoint()?;
        if geometry != self.geometry {
            return Err(DepositionError::InvalidInput(
                "integrated caloric gray boundary requires its original area, distance and zero coating resistance",
            ));
        }
        self.evaluate(specific_h_j_kg, bath_temperature_k, checkpoint)
            .map(|response| response.boundary)
    }
}

// a/(a+b), with nonnegative a,b and at least one positive. Avoid overflowing
// their sum, retaining the exact zero-exchange limit and its source tangent.
fn fraction(a: f64, b: f64) -> f64 {
    if a <= b {
        let ratio = a / b;
        ratio / (1.0 + ratio)
    } else {
        1.0 / (1.0 + b / a)
    }
}

fn finite(value: f64) -> Result<f64, DepositionError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(DepositionError::Arithmetic)
    }
}
