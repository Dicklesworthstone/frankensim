//! One source-neutral correlation chain rule and native surface heat response.
//!
//! Fluid owners provide their original transport, viscosity and expansion
//! tangents. The same correlation evidence, thermal direction, equilibrium
//! heat limits and conducting-surface interface apply to liquids and gases.

use super::*;

/// Source-specific wall inputs, resolved by the original fluid owner. The
/// shared evaluator applies one correction to the same physical tube section.
pub(super) enum TurbulentWallCorrection {
    Liquid(TransportSample),
    Air,
}

pub(super) trait ExteriorFluid {
    fn caloric(&self) -> &dyn ThermalLaw;
    fn energy_bounds(&self) -> [f64; 2];
    fn settings(&self) -> ExteriorConvectionSettings;
    /// Explicit fluid-specific correction on the existing turbulent section.
    /// Unselected sources retain the original constant-property behavior.
    fn turbulent_wall_correction(
        &self,
        _temperature: f64,
    ) -> Result<Option<TurbulentWallCorrection>, ExteriorConvectionError> {
        Ok(None)
    }
    fn sample_transport(
        &self,
        temperature: f64,
        derivatives: bool,
    ) -> Result<TransportSample, ExteriorConvectionError>;
    fn wall_viscosity(&self, temperature: f64)
    -> Result<(f64, f64, bool), ExteriorConvectionError>;
    fn expansion(
        &self,
        temperature: f64,
        sample: &TransportSample,
        derivative: bool,
    ) -> Result<(f64, f64), ExteriorConvectionError>;
}

pub(super) fn admit(
    caloric: &dyn ThermalLaw,
    energy_bounds_j_kg: [f64; 2],
    settings: ExteriorConvectionSettings,
) -> Result<(), ExteriorConvectionError> {
    use ExteriorConvectionError::Invalid;
    let [lo, hi] = energy_bounds_j_kg;
    if !lo.is_finite()
        || !hi.is_finite()
        || lo >= hi
        || !positive_interval(settings.fluid_temperature_bounds_k)
        || !positive_interval(settings.surface_temperature_bounds_k)
        || !settings.characteristic_length.value().is_finite()
        || settings.characteristic_length.value() <= 0.0
        || !settings.emissivity.is_finite()
        || !(0.0..=1.0).contains(&settings.emissivity)
    {
        return Err(Invalid(
            "finite energy and positive temperature/length bounds and passive emissivity are required",
        ));
    }
    for energy in energy_bounds_j_kg {
        let temperature = caloric
            .temperature_k(energy)
            .map_err(ExteriorConvectionError::Thermal)?;
        if !temperature.is_finite() || temperature <= 0.0 {
            return Err(Invalid("the original caloric endpoints must be physical"));
        }
    }
    if matches!(settings.law, SourceCellConvectionLaw::ForcedTubeSectionLiquidWall { .. }) {
        return Err(Invalid(
            "implicit liquid wall-Pr correction uses ForcedTubeSection with with_gnielinski_wall_correction",
        ));
    }
    if let SourceCellConvectionLaw::AidingVerticalPlate { settings } = settings.law {
        settings.validate()?;
    }
    if matches!(
        settings.law,
        SourceCellConvectionLaw::ForcedRectangularDuctSection { .. }
    ) {
        return Err(Invalid(
            "rectangular developing convection requires a whole four-wall area owner; a pointwise exterior law cannot admit that geometry",
        ));
    }
    settings
        .law
        .admit_geometry(settings.characteristic_length)?;
    settings
        .law
        .admit_property_temperature(settings.property_temperature)?;
    Ok(())
}

pub(super) fn evaluate(
    source: &impl ExteriorFluid,
    surface_temperature_k: f64,
    bath_temperature_k: f64,
) -> Result<ExteriorConvectionResponse, ExteriorConvectionError> {
    use ExteriorConvectionError::Invalid;
    let settings = source.settings();
    if !surface_temperature_k.is_finite()
        || !(settings.surface_temperature_bounds_k[0]..=settings.surface_temperature_bounds_k[1])
            .contains(&surface_temperature_k)
        || !bath_temperature_k.is_finite()
        || bath_temperature_k <= 0.0
    {
        return Err(Invalid(
            "exterior or bath temperature leaves its physical support",
        ));
    }
    let (property_temperature_k, film_factor) = match settings.property_temperature {
        ConvectionPropertyTemperature::Reservoir => (bath_temperature_k, 0.0),
        ConvectionPropertyTemperature::FilmMean => {
            (0.5 * surface_temperature_k + 0.5 * bath_temperature_k, 0.5)
        }
    };
    let sample = source.sample_transport(property_temperature_k, film_factor != 0.0)?;
    let state = sample.state;
    let length = settings.characteristic_length;
    let rho = state.density_kg_m3();
    let cp = state.specific_heat_j_kg_k();
    let k = state.thermal_conductivity_w_m_k();
    let mu = state.dynamic_viscosity_pa_s();
    let [drho, dcp, dk, dmu] = sample.derivatives.map(|slope| slope * film_factor);
    let dpr = finite(state.prandtl() * (dmu / mu + dcp / cp - dk / k))?;
    let varying_transport =
        film_factor != 0.0 && sample.derivatives.iter().any(|slope| *slope != 0.0);
    if let SourceCellConvectionLaw::ForcedPlateSection {
        correlation,
        speed_m_s,
        start,
    } = settings.law
    {
        let re = state.reynolds(speed_m_s, length.value())?;
        let dre = finite(re * (drho / rho - dmu / mu))?;
        let section = evaluate_flat_plate_section(
            correlation,
            re,
            state.prandtl(),
            start.value() / length.value(),
        )?;
        let coefficient = section
            .heat_transfer_coefficient_with_sensitivities(ThermalConductivity::new(k), length)?;
        // Require each potentially moving source group even at a zero
        // first-order chain factor; a fixed group needs no invented slope
        // at its empirical boundary. Use the SECTION's partials: the
        // turbulent leading-edge correction cancels at positive starts.
        let derivative = finite(
            coefficient_partial(&coefficient, "Re", dre, drho != 0.0 || dmu != 0.0)?
                + coefficient_partial(
                    &coefficient,
                    "Pr",
                    dpr,
                    dmu != 0.0 || dcp != 0.0 || dk != 0.0,
                )?
                + coefficient_partial(&coefficient, "k", dk, dk != 0.0)?,
        )?;
        let (_, upstream, nusselt, section_nusselt) = section.into_parts();
        return Ok(ExteriorConvectionResponse {
            property_temperature_k,
            transport: state,
            surface_transport: None,
            turbulent_property_correction: None,
            surface_dynamic_viscosity_pa_s: None,
            nusselt,
            plate_section: Some(PlateSectionDetails {
                start,
                end: length,
                upstream,
                nusselt: section_nusselt,
            }),
            tube_section: None,
            coefficient,
            coefficient_surface_derivative_w_m2_k2: derivative,
        });
    }
    if let SourceCellConvectionLaw::ForcedTubeSection {
        speed_m_s,
        start,
        end,
    } = settings.law
    {
        let re = state.reynolds(speed_m_s, length.value())?;
        let mut section = evaluate_turbulent_circular_duct_section(
            re,
            state.prandtl(),
            start.value() / length.value(),
            end.value() / length.value(),
        )?;
        let wall = source.turbulent_wall_correction(surface_temperature_k)?;
        let correction = match &wall {
            Some(TurbulentWallCorrection::Liquid(sample)) => {
                Some(TurbulentTubePropertyCorrection::Liquid {
                    wall_prandtl: sample.state.prandtl(),
                })
            }
            Some(TurbulentWallCorrection::Air) => Some(TurbulentTubePropertyCorrection::Air {
                bulk_temperature: Temperature::new(bath_temperature_k),
                wall_temperature: Temperature::new(surface_temperature_k),
            }),
            None => None,
        };
        let coefficient = if let Some(correction) = correction {
            section = section.with_property_correction(correction)?;
            section
                .heat_transfer_coefficient_with_sensitivities(ThermalConductivity::new(k), length)?
        } else {
            section.heat_transfer_coefficient(ThermalConductivity::new(k), length)?
        };
        let derivative = match &wall {
            Some(TurbulentWallCorrection::Liquid(sample)) => {
                let [_, dcp, dk, dmu] = sample.derivatives;
                let dpr = finite(
                    sample.state.prandtl()
                        * (dmu / sample.state.dynamic_viscosity_pa_s()
                            + dcp / sample.state.specific_heat_j_kg_k()
                            - dk / sample.state.thermal_conductivity_w_m_k()),
                )?;
                coefficient_partial(
                    &coefficient,
                    "Pr_wall",
                    dpr,
                    dcp != 0.0 || dk != 0.0 || dmu != 0.0,
                )?
            }
            Some(TurbulentWallCorrection::Air) => {
                coefficient_partial(&coefficient, "T_wall", 1.0, true)?
            }
            None => 0.0,
        };
        let (_, _, upstream, nusselt, section_nusselt) = section.into_parts();
        // Bulk Re/Pr/k stay fixed during surface elimination. Only the
        // explicitly selected liquid-property or air-temperature chain moves.
        return Ok(ExteriorConvectionResponse {
            property_temperature_k,
            transport: state,
            surface_transport: match wall {
                Some(TurbulentWallCorrection::Liquid(sample)) => Some(sample.state),
                _ => None,
            },
            turbulent_property_correction: correction,
            surface_dynamic_viscosity_pa_s: None,
            nusselt,
            plate_section: None,
            tube_section: Some(TurbulentTubeSectionDetails {
                diameter: length,
                start,
                end,
                upstream,
                nusselt: section_nusselt,
            }),
            coefficient,
            coefficient_surface_derivative_w_m2_k2: derivative,
        });
    }
    if let SourceCellConvectionLaw::ForcedLaminarTubeSection {
        speed_m_s,
        start,
        end,
    } = settings.law
    {
        let re = state.reynolds(speed_m_s, length.value())?;
        let section = evaluate_circular_duct_section(
            re,
            state.prandtl(),
            start.value() / length.value(),
            end.value() / length.value(),
        )?;
        let coefficient = section.heat_transfer_coefficient(ThermalConductivity::new(k), length)?;
        let (_, _, upstream, nusselt, section_nusselt) = section.into_parts();
        // Bulk-held Hausen has no wall-property correction. Its actual
        // section coefficient therefore has exactly zero wall derivative.
        return Ok(ExteriorConvectionResponse {
            property_temperature_k,
            transport: state,
            surface_transport: None,
            turbulent_property_correction: None,
            surface_dynamic_viscosity_pa_s: None,
            nusselt,
            plate_section: None,
            tube_section: Some(TubeSectionDetails {
                diameter: length,
                start,
                end,
                upstream,
                nusselt: section_nusselt,
            }),
            coefficient,
            coefficient_surface_derivative_w_m2_k2: 0.0,
        });
    }
    let mut wall_viscosity = None;
    let mut varying_wall_viscosity = false;
    let (nusselt, dnu) = match settings.law {
        SourceCellConvectionLaw::Forced {
            correlation,
            speed_m_s,
            settings,
        } => {
            if correlation == CorrelationId::DittusBoelter
                && ((surface_temperature_k < bath_temperature_k
                    && settings.direction == ThermalDirection::HeatingFluid)
                    || (surface_temperature_k > bath_temperature_k
                        && settings.direction == ThermalDirection::CoolingFluid))
            {
                return Err(Invalid(
                    "Dittus-Boelter direction disagrees with actual surface/bath heat flow",
                ));
            }
            let re = state.reynolds(speed_m_s, length.value())?;
            let dre = finite(re * (drho / rho - dmu / mu))?;
            let mut inputs =
                CorrelationInputs::forced(re, state.prandtl()).with_direction(settings.direction);
            if let Some(ratio) = settings.length_over_hydraulic_diameter {
                inputs = inputs.with_length_ratio(ratio);
            }
            if let Some(ratio) = settings.aspect_ratio {
                inputs = inputs.with_aspect_ratio(ratio);
            }
            let dratio = if matches!(correlation,
                CorrelationId::SiederTateTurbulent | CorrelationId::SiederTateLaminar)
            {
                let (wall, dwall, varies) = source.wall_viscosity(surface_temperature_k)?;
                varying_wall_viscosity = varies;
                wall_viscosity = Some(wall);
                let ratio = finite(mu / wall)?;
                inputs = inputs.with_viscosity_ratio(ratio);
                finite(-ratio * dwall / wall)?
            } else {
                0.0
            };
            let nu = inputs.evaluate_with_sensitivities(correlation)?;
            let derivative = partial(&nu, "Re", dre)?
                + partial(&nu, "Pr", dpr)?
                + partial(&nu, "mu_bulk_over_surface", dratio)?;
            (nu, finite(derivative)?)
        }
        SourceCellConvectionLaw::ForcedSphere { speed_m_s } => {
            let re = state.reynolds(speed_m_s, length.value())?;
            let (wall, dwall, varies) = source.wall_viscosity(surface_temperature_k)?;
            varying_wall_viscosity = varies;
            wall_viscosity = Some(wall);
            let ratio = finite(mu / wall)?;
            let nu = CorrelationInputs::forced(re, state.prandtl())
                .with_viscosity_ratio(ratio)
                .evaluate_with_sensitivities(CorrelationId::WhitakerSphere)?;
            let derivative = partial(&nu, "mu_bulk_over_surface", finite(-ratio * dwall / wall)?)?;
            (nu, derivative)
        }
        SourceCellConvectionLaw::AidingVerticalPlate { settings } => {
            let (beta, dbeta_dt) =
                source.expansion(property_temperature_k, &sample, film_factor != 0.0)?;
            let delta = surface_temperature_k - bath_temperature_k;
            let groups = settings.groups(
                state,
                beta,
                Temperature::new(surface_temperature_k),
                Temperature::new(bath_temperature_k),
                length,
            )?;
            let constant_zero_drive =
                settings.gravity_m_s2 == 0.0 || (beta == 0.0 && dbeta_dt == 0.0);
            if groups.rayleigh == 0.0 && !constant_zero_drive {
                return Err(Invalid(
                    "zero-Rayleigh mixed convection needs the native heat-tangent interface",
                ));
            }
            let dre = finite(groups.reynolds * (drho / rho - dmu / mu))?;
            let factor = finite(
                settings.gravity_m_s2 * length.value().powi(3)
                    / (state.kinematic_viscosity_m2_s() * state.thermal_diffusivity_m2_s()),
            )?;
            // Differentiate the signed drive before taking its magnitude;
            // neither beta nor Delta T is used as a divisor.
            let dra = finite(
                factor
                    * (beta
                        + delta
                            * (film_factor * dbeta_dt
                                + beta * (2.0 * drho / rho + dcp / cp - dmu / mu - dk / k)))
                    * beta.signum()
                    * delta.signum(),
            )?;
            let nu = CorrelationInputs::mixed(groups.reynolds, groups.rayleigh, groups.prandtl)
                .evaluate_with_sensitivities(CorrelationId::ManningQureshiAidingVerticalPlate)?;
            require_partial(&nu, "Ra", groups.rayleigh != 0.0)?;
            let derivative =
                partial(&nu, "Re", dre)? + partial(&nu, "Ra", dra)? + partial(&nu, "Pr", dpr)?;
            (nu, finite(derivative)?)
        }
        SourceCellConvectionLaw::VerticalPlate { gravity_m_s2 }
        | SourceCellConvectionLaw::Natural { gravity_m_s2, .. } => {
            let geometry = match settings.law {
                SourceCellConvectionLaw::Natural { geometry, .. } => geometry,
                _ => NaturalConvectionGeometry::VerticalPlate,
            };
            let (beta, dbeta_dt) =
                source.expansion(property_temperature_k, &sample, film_factor != 0.0)?;
            let delta = surface_temperature_k - bath_temperature_k;
            let ra = state.rayleigh(gravity_m_s2, beta, delta, length.value())?;
            let constant_zero_drive = geometry == NaturalConvectionGeometry::Sphere
                && (gravity_m_s2 == 0.0
                    || (beta == 0.0 && (film_factor == 0.0 || dbeta_dt == 0.0)));
            if ra == 0.0 && !constant_zero_drive {
                return Err(Invalid(
                    "zero-Rayleigh natural convection has no finite coefficient tangent on this interface",
                ));
            }
            // No division by beta or Delta T: signed zero/negative source
            // expansion must not turn into a manufactured ideal-gas law.
            let factor = finite(
                gravity_m_s2 * length.value().powi(3)
                    / (state.kinematic_viscosity_m2_s() * state.thermal_diffusivity_m2_s()),
            )?;
            let dra = finite(
                factor
                    * (beta
                        + delta
                            * (film_factor * dbeta_dt
                                + beta * (2.0 * drho / rho + dcp / cp - dmu / mu - dk / k))),
            )?;
            let nu = geometry.evaluate_with_sensitivities(ra, state.prandtl())?;
            // A zero first-order chain factor does not justify selecting a
            // derivative at a discontinuous natural-regime boundary.
            require_partial(&nu, "Ra", ra != 0.0)?;
            let derivative = partial(&nu, "Ra", dra * ra.signum())? + partial(&nu, "Pr", dpr)?;
            (nu, finite(derivative)?)
        }
        SourceCellConvectionLaw::ForcedPlateSection { .. }
        | SourceCellConvectionLaw::ForcedTubeSection { .. }
        | SourceCellConvectionLaw::ForcedLaminarTubeSection { .. } => {
            unreachable!("finite sections return their own coefficient above")
        }
        SourceCellConvectionLaw::ForcedRectangularDuctSection { .. } => {
            unreachable!("pointwise rectangular geometry is refused at construction")
        }
        SourceCellConvectionLaw::ForcedTubeSectionLiquidWall { .. } => {
            return Err(Invalid(
                "implicit liquid wall-Pr correction requires its original prepared liquid builder",
            ));
        }
    };
    // Moving source properties can have a zero aggregate group derivative
    // at a local extremum. That does not establish differentiability of a
    // missing empirical partial at a regime boundary.
    if varying_transport {
        for axis in ["Re", "Pr"] {
            // A stopped external stream stays at Re=0 regardless of the
            // fluid-property response; it needs no Reynolds derivative.
            let varies = nusselt
                .groups()
                .get(axis)
                .is_some_and(|&value| axis != "Re" || value != 0.0);
            require_partial(&nusselt, axis, varies)?;
        }
    }
    require_partial(&nusselt, "mu_bulk_over_surface", varying_wall_viscosity)?;
    let coefficient = nusselt.heat_transfer_coefficient(ThermalConductivity::new(k), length)?;
    let derivative = finite((dnu * k + nusselt.evidence().value * dk) / length.value())?;
    Ok(ExteriorConvectionResponse {
        property_temperature_k,
        transport: state,
        surface_transport: None,
        turbulent_property_correction: None,
        surface_dynamic_viscosity_pa_s: wall_viscosity,
        nusselt,
        plate_section: None,
        tube_section: None,
        coefficient,
        coefficient_surface_derivative_w_m2_k2: derivative,
    })
}

pub(super) fn flux_state_at_surface(
    source: &impl ExteriorFluid,
    energy: f64,
    surface: f64,
    bath: f64,
) -> Result<ImplicitGrayExteriorFluxState, DepositionError> {
    let settings = source.settings();
    let energy_bounds = source.energy_bounds();
    // At exact sphere equilibrium Nu-2 = O(|Delta T|^(1/4)); mixed flow
    // similarly approaches its pure-forced coefficient (zero at rest).
    // The buoyant contribution to signed heat is o(Delta T), even when
    // dh/dTs diverges. Return the heat tangent, not a fictitious finite
    // coefficient derivative. Plate/cylinder cards retain their minima.
    if surface == bath
        && matches!(
            settings.law,
            SourceCellConvectionLaw::Natural {
                geometry: NaturalConvectionGeometry::Sphere,
                ..
            } | SourceCellConvectionLaw::AidingVerticalPlate { .. }
        )
    {
        if !energy.is_finite()
            || !(energy_bounds[0]..=energy_bounds[1]).contains(&energy)
            || !surface.is_finite()
            || !(settings.surface_temperature_bounds_k[0]
                ..=settings.surface_temperature_bounds_k[1])
                .contains(&surface)
        {
            return Err(DepositionError::MaterialDomain);
        }
        source.caloric().temperature_k(energy)?;
        let resolve = || -> Result<f64, ExteriorConvectionError> {
            let sample = source.sample_transport(surface, false)?;
            let (beta, _) = source.expansion(surface, &sample, false)?;
            let length = settings.characteristic_length;
            let nu = match settings.law {
                SourceCellConvectionLaw::Natural { gravity_m_s2, .. } => {
                    let ra = sample
                        .state
                        .rayleigh(gravity_m_s2, beta, 0.0, length.value())?;
                    NaturalConvectionGeometry::Sphere.evaluate(ra, sample.state.prandtl())?
                }
                SourceCellConvectionLaw::AidingVerticalPlate { settings } => {
                    let groups = settings.groups(
                        sample.state,
                        beta,
                        Temperature::new(surface),
                        Temperature::new(bath),
                        length,
                    )?;
                    match fs_convection::evaluate(
                        CorrelationId::ManningQureshiAidingVerticalPlate,
                        CorrelationInputs::mixed(groups.reynolds, groups.rayleigh, groups.prandtl),
                    ) {
                        Ok(nu) => nu,
                        // The original card validates all groups, including
                        // Pr, before its dedicated zero-driver refusal.
                        // Only at EXACT thermal equilibrium do heat and its
                        // tangent have the zero resting-flow limit. There
                        // is still no Nusselt evidence for this state.
                        Err(CorrelationError::InvalidGroup {
                            axis: "max(Re,Ra)", ..
                        }) if groups.reynolds == 0.0 && groups.rayleigh == 0.0 => return Ok(0.0),
                        Err(error) => return Err(error.into()),
                    }
                }
                _ => unreachable!("equilibrium branch only admits sphere and mixed plate"),
            };
            Ok(nu
                .heat_transfer_coefficient(
                    ThermalConductivity::new(sample.state.thermal_conductivity_w_m_k()),
                    length,
                )?
                .value
                .value())
        };
        let htc = resolve().map_err(|_| DepositionError::MaterialDomain)?;
        return Ok(ImplicitGrayExteriorFluxState {
            coefficients: ImplicitGraySurfaceState {
                emissivity: settings.emissivity,
                emissivity_enthalpy_derivative_kg_j: 0.0,
                convection_w_m2_k: htc,
                convection_enthalpy_derivative_kg_w_m2_k_j: 0.0,
            },
            emissivity_temperature_derivative_per_k: 0.0,
            convection_temperature_tangent_w_m2_k: htc,
        });
    }
    let state = state_at_surface(source, energy, surface, bath)?;
    let tangent = state.coefficients.convection_w_m2_k
        + (surface - bath) * state.convection_temperature_derivative_w_m2_k2;
    if !tangent.is_finite() {
        return Err(DepositionError::Arithmetic);
    }
    Ok(ImplicitGrayExteriorFluxState {
        coefficients: state.coefficients,
        emissivity_temperature_derivative_per_k: state.emissivity_temperature_derivative_per_k,
        convection_temperature_tangent_w_m2_k: tangent,
    })
}

pub(super) fn state_at_surface(
    source: &impl ExteriorFluid,
    energy: f64,
    surface: f64,
    bath: f64,
) -> Result<ImplicitGrayExteriorState, DepositionError> {
    let settings = source.settings();
    let energy_bounds = source.energy_bounds();
    if !energy.is_finite() || !(energy_bounds[0]..=energy_bounds[1]).contains(&energy) {
        return Err(DepositionError::MaterialDomain);
    }
    source.caloric().temperature_k(energy)?;
    let response = evaluate(source, surface, bath).map_err(|_| DepositionError::MaterialDomain)?;
    Ok(ImplicitGrayExteriorState {
        coefficients: ImplicitGraySurfaceState {
            emissivity: settings.emissivity,
            emissivity_enthalpy_derivative_kg_j: 0.0,
            convection_w_m2_k: response.coefficient.value.value(),
            convection_enthalpy_derivative_kg_w_m2_k_j: 0.0,
        },
        emissivity_temperature_derivative_per_k: 0.0,
        convection_temperature_derivative_w_m2_k2: response.coefficient_surface_derivative_w_m2_k2,
    })
}

fn require_partial(
    evaluation: &NusseltEvaluation,
    axis: &'static str,
    varies: bool,
) -> Result<(), ExteriorConvectionError> {
    if varies && !evaluation.evidence().sensitivity.d_qoi.contains_key(axis) {
        Err(ExteriorConvectionError::MissingDerivative { axis })
    } else {
        Ok(())
    }
}
fn partial(
    evaluation: &NusseltEvaluation,
    axis: &'static str,
    input_derivative: f64,
) -> Result<f64, ExteriorConvectionError> {
    if input_derivative == 0.0 {
        return Ok(0.0);
    }
    let derivative = evaluation
        .evidence()
        .sensitivity
        .d_qoi
        .get(axis)
        .ok_or(ExteriorConvectionError::MissingDerivative { axis })?;
    finite(derivative * input_derivative)
}
fn coefficient_partial(
    coefficient: &Evidence<HeatTransferCoefficient>,
    axis: &'static str,
    input_derivative: f64,
    source_varies: bool,
) -> Result<f64, ExteriorConvectionError> {
    if !source_varies {
        return Ok(0.0);
    }
    let derivative = coefficient
        .sensitivity
        .d_qoi
        .get(axis)
        .ok_or(ExteriorConvectionError::MissingDerivative { axis })?;
    finite(derivative * input_derivative)
}
