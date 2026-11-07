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
    ImplicitBoundaryGeometry, ImplicitGrayExteriorFluxState, ImplicitGrayExteriorLaw,
    ImplicitGraySurfaceLaw, ImplicitGraySurfaceState,
};
use fs_thermal::enthalpy::{DepositionError, ThermalLaw};

// The same rounded CODATA value as the original gray boundary owners.
const SIGMA: f64 = 5.670_374_419e-8;

struct BoundSurface<'a> {
    law: &'a dyn ImplicitGraySurfaceLaw,
    bounds_j_kg: [f64; 2],
}

struct BoundExterior<'a> {
    law: &'a dyn ImplicitGrayExteriorLaw,
    energy_bounds_j_kg: [f64; 2],
    temperature_bounds_k: [f64; 2],
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
/// on the SAME exterior surface. Material k(T) is integrated to its supported
/// material face; an optional coating/contact resistance separates that face
/// from the exterior surface. Only the material path requires conductivity
/// support, not the exterior coating face or bath states.
///
/// This is a full-view gray boundary with an optional fixed, zero-storage
/// coating/contact resistance before its shared exterior surface. That
/// resistance excludes the convective film and radiation. Participating
/// radiation and finite enclosure storage require their own physical owners.
/// The kernel allocates no memory; an optional original source retains its own
/// evaluation costs. Spatial and temporal accuracy remain the native solver's
/// responsibility rather than a consequence of this constitutive closure.
/// The geometry-aware native query includes density-driven normal thickness
/// through the same coating/exterior balance, at fixed area and bath states.
pub struct IntegratedCaloricGrayBoundary<'l, 'a, 'm> {
    layer: &'l CaloricConductionLayer<'a, 'm>,
    geometry: ImplicitBoundaryGeometry,
    flux_geometry: ReservoirLayerGeometry,
    material_multiplier_m: f64,
    coating_resistance_k_w: f64,
    fixed: ImplicitGraySurfaceState,
    surface: Option<BoundSurface<'l>>,
    exterior: Option<BoundExterior<'l>>,
    radiative_temperature_k: f64,
    config: LayerConductionConfig,
}

impl<'l, 'a, 'm> IntegratedCaloricGrayBoundary<'l, 'a, 'm> {
    /// Bind the original positive area and centre-to-face distance, with a
    /// nonnegative areal coating/contact resistance in m² K/W. It is converted
    /// to total K/W once; HTC owns only the exterior convective film.
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
            || !geometry.surface_resistance_m2_k_w.is_finite()
            || geometry.surface_resistance_m2_k_w < 0.0
            || config.max_iterations == 0
            || !config.absolute_heat_tolerance_w.is_finite()
            || config.absolute_heat_tolerance_w < 0.0
            || !config.relative_heat_tolerance.is_finite()
            || config.relative_heat_tolerance < 0.0
        {
            return Err(DepositionError::InvalidInput(
                "integrated gray boundary needs passive coefficients, a positive bath, nonnegative coating resistance and bounded solve controls",
            ));
        }
        let coating_resistance_k_w = geometry.surface_resistance_m2_k_w / geometry.area_m2;
        if !coating_resistance_k_w.is_finite()
            || (geometry.surface_resistance_m2_k_w > 0.0 && coating_resistance_k_w == 0.0)
        {
            return Err(DepositionError::Arithmetic);
        }
        Ok(Self {
            layer,
            geometry,
            flux_geometry,
            material_multiplier_m: geometry.area_m2 / geometry.half_distance_m,
            coating_resistance_k_w,
            fixed,
            surface: None,
            exterior: None,
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
    /// included afterward. A source requiring the solved exterior temperature
    /// uses `with_exterior_law`; this method never substitutes that temperature
    /// into a source declared on the bulk coordinate.
    pub fn with_surface_law(
        mut self,
        law: &'l dyn ImplicitGraySurfaceLaw,
    ) -> Result<Self, DepositionError> {
        let [lo, hi] = law.enthalpy_bounds_j_kg();
        let [source_lo, source_hi] = self.layer.energy_bounds_j_kg();
        if self.surface.is_some()
            || self.exterior.is_some()
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

    /// Uniquely replace zero fixed coefficients with an actual exterior law.
    ///
    /// Energy bounds lie inside the SAME original caloric/conductivity chart.
    /// Exterior temperature support is independent of that chart, the bulk
    /// temperature and either bath. Each surface trial queries its original
    /// source using `(e_bulk, Ts, Tgas)`; no query uses a substituted bulk state.
    ///
    /// Both source partials enter the exact implicit tangent, including signed
    /// exterior-temperature feedback and coefficient onset at zero. A final
    /// branch needs positive combined material/exterior differential stiffness
    /// to eliminate this zero-storage surface stably; a singular or locally
    /// unstable branch refuses. This local check and the existing bracket do
    /// not claim global uniqueness for an arbitrary nonlinear source.
    pub fn with_exterior_law(
        mut self,
        law: &'l dyn ImplicitGrayExteriorLaw,
    ) -> Result<Self, DepositionError> {
        let [lo, hi] = law.energy_bounds_j_kg();
        let [source_lo, source_hi] = self.layer.energy_bounds_j_kg();
        let [surface_lo, surface_hi] = law.surface_temperature_bounds_k();
        if self.surface.is_some()
            || self.exterior.is_some()
            || self.fixed.emissivity != 0.0
            || self.fixed.convection_w_m2_k != 0.0
            || !lo.is_finite()
            || !hi.is_finite()
            || lo >= hi
            || lo < source_lo
            || hi > source_hi
            || !surface_lo.is_finite()
            || surface_lo <= 0.0
            || !surface_hi.is_finite()
            || surface_lo >= surface_hi
            || !core::ptr::addr_eq(
                self.layer.caloric_law() as *const dyn ThermalLaw,
                law.caloric_law() as *const dyn ThermalLaw,
            )
        {
            return Err(DepositionError::InvalidInput(
                "integrated exterior source needs sole ownership, its original supported caloric object and positive surface bounds",
            ));
        }
        self.exterior = Some(BoundExterior {
            law,
            energy_bounds_j_kg: [lo, hi],
            temperature_bounds_k: [surface_lo, surface_hi],
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
        self.evaluate_with_geometry_tangent(
            specific_energy_j_kg,
            gas_temperature_k,
            self.geometry,
            0.0,
            checkpoint,
        )
    }

    fn evaluate_with_geometry_tangent(
        &self,
        specific_energy_j_kg: f64,
        gas_temperature_k: f64,
        geometry: ImplicitBoundaryGeometry,
        relative_distance_slope_kg_j: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<GrayTwoBathResponse, DepositionError> {
        checkpoint()?;
        let (flux_geometry, material_multiplier_m) = if geometry == self.geometry {
            (self.flux_geometry, self.material_multiplier_m)
        } else {
            let current =
                ReservoirLayerGeometry::planar(geometry.area_m2, geometry.half_distance_m)?;
            (current, current.conductance_multiplier_m())
        };
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
        // An exterior source owns its actual Ts support independently of the
        // bulk chart and bath states. Fixed/bulk coefficients retain the exact
        // original support and arithmetic, including the zero-coating path.
        let surface_span = if let Some(exterior) = &self.exterior {
            let [lo, hi] = exterior.energy_bounds_j_kg;
            if !(lo..=hi).contains(&specific_energy_j_kg) {
                return Err(DepositionError::MaterialDomain);
            }
            exterior.temperature_bounds_k
        } else if self.coating_resistance_k_w == 0.0 {
            curve.temperature_range_k()
        } else {
            [
                material
                    .temperature_k
                    .min(gas_temperature_k)
                    .min(self.radiative_temperature_k),
                material
                    .temperature_k
                    .max(gas_temperature_k)
                    .max(self.radiative_temperature_k),
            ]
        };
        let response = flux_geometry
            .evaluate_material_to_two_baths_with_resistance(
                curve,
                material.temperature_k,
                [gas_temperature_k, self.radiative_temperature_k],
                surface_span,
                self.coating_resistance_k_w,
                self.config,
                |surface| {
                    if let Some(exterior) = &self.exterior {
                        // The original kernel brackets Ts inside the source
                        // interval and polls around every callback. Its saved
                        // payload carries the final partials without a second
                        // source query after convergence.
                        let source = exterior.law.flux_state_at_surface(
                            specific_energy_j_kg,
                            surface,
                            gas_temperature_k,
                        )?;
                        source.admit()?;
                        self.exterior_source_exchange(surface, gas_temperature_k, source)
                    } else {
                        self.exterior_exchange(
                            surface,
                            gas_temperature_k,
                            coefficients,
                            coefficients.convection_w_m2_k,
                        )
                    }
                },
                &mut *checkpoint,
            )
            .map_err(|error| match error {
                ConvectiveLayerError::Thermal(error) | ConvectiveLayerError::Film(error) => error,
                ConvectiveLayerError::SurfaceBalance { .. } => DepositionError::InvalidInput(
                    "integrated caloric gray surface balance did not converge",
                ),
            })?;
        let mut cell_slope = finite(material_multiplier_m * material.conductivity_w_m_k)?;
        let mut surface_slope = finite(
            material_multiplier_m
                * curve.conductivity_w_m_k(response.material_face_temperature_k)?,
        )?;
        if cell_slope <= 0.0 || surface_slope <= 0.0 {
            return Err(DepositionError::Arithmetic);
        }
        let mut coating_response = 1.0;
        if self.coating_resistance_k_w > 0.0 {
            // Eliminate Tmaterial = Texterior + R*Q analytically. The effective
            // material derivatives are Gcell/(1+R*Gface) and Gface/(1+R*Gface).
            // Split the arithmetic to retain tiny R without 1/R overflow and
            // large R without R*Gface overflow in the denominator.
            let ratio = self.coating_resistance_k_w * surface_slope;
            if ratio <= 1.0 {
                let factor = 1.0 / (1.0 + ratio);
                coating_response = factor;
                cell_slope *= factor;
                surface_slope *= factor;
            } else {
                let coating_slope = self.coating_resistance_k_w.recip();
                let combined = coating_slope / (1.0 + coating_slope / surface_slope);
                coating_response = combined / surface_slope;
                cell_slope = finite(combined * (cell_slope / surface_slope))?;
                surface_slope = combined;
            }
            if cell_slope <= 0.0 || surface_slope <= 0.0 {
                return Err(DepositionError::Arithmetic);
            }
        }
        // Q=C*integral(Tmaterial..Tcell) k dT, C=A/d. Eliminating
        // Tmaterial=Texterior+Rcoat*Q leaves the geometric forcing
        // -Q*(d'/d)/(1+Rcoat*C*k(Tmaterial)). Use the retained integral
        // power and SAME accepted source payload, without a second solve.
        let geometry_forcing = if relative_distance_slope_kg_j == 0.0 {
            0.0
        } else {
            finite(
                (-response.solid_power_to_surface_w * coating_response)
                    * relative_distance_slope_kg_j,
            )?
        };
        let exterior_slope = response.film.temperature_derivative_w_k;
        // Differentiate the material/coat/exterior balance with the eliminated
        // leg slopes. The integral's two endpoint conductivities already
        // include the entire k(T) response; adding a bulk dk/de term would
        // double count it. Source F_e is taken at fixed exterior Ts and may
        // be signed, including zero-coefficient onset.
        let derivative = if exterior_slope >= 0.0 {
            // Preserve the original fixed/bulk coefficient arithmetic exactly.
            let exterior_weight = fraction(exterior_slope, surface_slope);
            let mut storage = (exterior_weight * cell_slope) * material.temperature_slope_kg_k_j;
            if relative_distance_slope_kg_j != 0.0 {
                storage = finite(exterior_weight.mul_add(geometry_forcing, storage))?;
            }
            finite(
                fraction(surface_slope, exterior_slope)
                    .mul_add(response.film.energy_derivative_w_kg_j, storage),
            )?
        } else {
            let denominator = finite(surface_slope + exterior_slope)?;
            if denominator <= 0.0 {
                return Err(DepositionError::InvalidInput(
                    "integrated exterior source has a singular or locally unstable surface balance",
                ));
            }
            // A passive source can have negative differential feedback. Keep
            // that signed response; normalized positive-only fractions would
            // silently turn it into a different Newton/adjoint tangent.
            let mut body = finite(cell_slope * material.temperature_slope_kg_k_j)?;
            if relative_distance_slope_kg_j != 0.0 {
                body = finite(body + geometry_forcing)?;
            }
            finite((exterior_slope / denominator).mul_add(
                body,
                (surface_slope / denominator) * response.film.energy_derivative_w_kg_j,
            ))?
        };
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
        convection_temperature_tangent_w_m2_k: f64,
    ) -> Result<([f64; 2], [f64; 2], SurfaceTangent), DepositionError> {
        let area = self.geometry.area_m2;
        let gas_delta = surface - gas;
        let radiation_delta = surface - self.radiative_temperature_k;
        let convection_g = finite(area * source.convection_w_m2_k)?;
        let convection_slope = finite(area * convection_temperature_tangent_w_m2_k)?;
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
                temperature_derivative_w_k: finite(convection_slope + radiation_slope)?,
                energy_derivative_w_kg_j: finite(convection_energy + radiation_energy)?,
            },
        ))
    }

    fn exterior_source_exchange(
        &self,
        surface: f64,
        gas: f64,
        source: ImplicitGrayExteriorFluxState,
    ) -> Result<([f64; 2], [f64; 2], SurfaceTangent), DepositionError> {
        // Form the convective flux tangent directly. Reconstructing it from
        // dh/dTs would reject valid equilibrium limits or subtract large,
        // nearly cancelling coefficient contributions.
        let (powers, conductances, mut tangent) = self.exterior_exchange(
            surface,
            gas,
            source.coefficients,
            source.convection_temperature_tangent_w_m2_k,
        )?;
        let area = self.geometry.area_m2;
        if source.emissivity_temperature_derivative_per_k != 0.0 {
            // Factored fourth-power difference also retains the source tangent
            // at epsilon=0. No quotient by emissivity or radiative power occurs.
            let scale = finite((area * SIGMA) * source.emissivity_temperature_derivative_per_k)?;
            let sum = finite(surface + self.radiative_temperature_k)?;
            let squares = finite(surface.mul_add(
                surface,
                self.radiative_temperature_k * self.radiative_temperature_k,
            ))?;
            tangent.temperature_derivative_w_k = finite(finite((scale * sum) * squares)?.mul_add(
                surface - self.radiative_temperature_k,
                tangent.temperature_derivative_w_k,
            ))?;
        }
        Ok((powers, conductances, tangent))
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
                "integrated caloric gray boundary requires its original area, distance and coating resistance",
            ));
        }
        self.evaluate(specific_h_j_kg, bath_temperature_k, checkpoint)
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
        checkpoint()?;
        if geometry.area_m2 != self.geometry.area_m2
            || geometry.surface_resistance_m2_k_w != self.geometry.surface_resistance_m2_k_w
            || !half_distance_relative_enthalpy_derivative_kg_j.is_finite()
        {
            return Err(DepositionError::InvalidInput(
                "density-driven caloric gray boundary requires its original area, coating and finite length tangent",
            ));
        }
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
