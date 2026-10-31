//! Film-temperature properties on the original finite gas/body heat residual.
//! Film evaluation has no independent energy or mass; bulk gas keeps both.
use super::super::evaluate_natural_from_material;
use super::*;
use fs_material::thermochemical_gas::ThermochemicalGasModel;
use fs_qty::Temperature;

/// The existing correlation/geometry, with an explicit mean-film convention.
#[derive(Clone, Copy, Debug)]
pub enum FilmGasHeatLaw {
    Forced(ForcedGasSurface),
    Natural(NaturalGasSurface),
}

/// Source bank and constitutive boundary law for one finite gas/body patch.
/// The bank may contain more species than the current gas, or different initial
/// fractions. Each query reconstructs the ACTUAL active gas composition and
/// preserves its original source claims before evaluating film properties.
#[derive(Clone, Copy)]
pub struct FilmGasSurface<'a> {
    pub model: &'a ThermochemicalGasModel,
    pub law: FilmGasHeatLaw,
}

/// The actual bulk gas/body states and the separate property evaluation.
/// Bulk gas remains the thermodynamic/acoustic state; film never replaces it.
#[derive(Debug)]
pub struct FilmGasConvectionSample {
    pub interface: FiniteGasConvectionSample,
    /// None only for an explicitly isolated zero-area patch.
    pub property_gas: Option<ResolvedThermochemicalGasState>,
}

impl FilmGasHeatLaw {
    fn area(self) -> Result<f64, MaterialConvectionError> {
        match self {
            Self::Forced(surface) => {
                surface.validate()?;
                Ok(surface.area_m2)
            }
            Self::Natural(surface) => {
                surface.validate()?;
                Ok(surface.area_m2)
            }
        }
    }
}

impl FilmGasSurface<'_> {
    /// Re-query source-pinned properties at (T_wall + T_bulk)/2 and actual bulk
    /// pressure. Both signed heat flow and Rayleigh forcing retain the FULL
    /// wall-to-bulk temperature contrast. The film is not a replacement gas
    /// endpoint, a second parcel, or an invented finite-capacity reservoir.
    ///
    /// The original correlation must support this property convention and the
    /// declared geometry. The model remains dilute ideal, frozen composition
    /// during this evaluation, with no boundary-layer or chemical heat solver.
    /// The callback works with existing powered/reactive/multi-boundary owners;
    /// their complete exposed-area checks and source budgets remain mandatory.
    pub fn resolve(
        self,
        gas: &ResolvedThermochemicalGasState,
        body: EquilibriumPhaseState,
    ) -> Result<ConductanceSample<FilmGasConvectionSample>, MaterialConvectionError> {
        let area_m2 = self.law.area()?;
        let wall = body.temperature_k();
        let bulk = gas.state().temperature;
        if !wall.is_finite() || wall <= 0.0 {
            return Err(invalid("film convection needs a positive body temperature"));
        }
        if area_m2 == 0.0 {
            return Ok(ConductanceSample {
                conductance_w_k: 0.0,
                source: FilmGasConvectionSample {
                    interface: FiniteGasConvectionSample {
                        gas: gas.clone(),
                        body,
                        area_m2,
                        evaluation: None,
                    },
                    property_gas: None,
                },
            });
        }
        let film = self
            .model
            .evaluate_temperature_from_snapshot(gas, Temperature::new(0.5 * wall + 0.5 * bulk))?;
        let evaluated = match self.law {
            FilmGasHeatLaw::Forced(surface) => {
                let mut settings = surface.settings;
                if wall > bulk {
                    settings.direction = ThermalDirection::HeatingFluid;
                } else if wall < bulk {
                    settings.direction = ThermalDirection::CoolingFluid;
                }
                Some(evaluate_forced_from_material(
                    surface.correlation,
                    &film,
                    surface.speed_m_s,
                    surface.characteristic_length,
                    settings,
                )?)
            }
            FilmGasHeatLaw::Natural(surface) => {
                if wall == bulk {
                    None
                } else {
                    Some(evaluate_natural_from_material(
                        surface.geometry,
                        &film,
                        surface.gravity_m_s2,
                        1.0 / film.state().temperature,
                        Temperature::new(wall),
                        Temperature::new(bulk),
                        surface.characteristic_length,
                    )?)
                }
            }
        };
        let evaluation = evaluated.map(|value| {
            let (transport, nusselt, coefficient) = value.into_parts();
            GasConvectionEvaluation {
                transport,
                nusselt,
                coefficient,
            }
        });
        let conductance_w_k = match &evaluation {
            None => 0.0,
            Some(value) => {
                let conductance = area_m2 * value.coefficient.value.value();
                if !conductance.is_finite() || conductance <= 0.0 {
                    return Err(invalid(
                        "film convective conductance overflowed or vanished",
                    ));
                }
                conductance
            }
        };
        Ok(ConductanceSample {
            conductance_w_k,
            source: FilmGasConvectionSample {
                interface: FiniteGasConvectionSample {
                    gas: gas.clone(),
                    body,
                    area_m2,
                    evaluation,
                },
                property_gas: Some(film),
            },
        })
    }

    /// Advance heat, source-backed gas energy, finite body enthalpy and piston
    /// momentum together through the existing heat/motion residual. The same
    /// proposed bulk gas determines pressure work, storage and acoustic state;
    /// only the transfer coefficient uses its source-pinned film properties.
    /// Both original physical endpoints are returned together or neither.
    #[allow(clippy::too_many_arguments)]
    pub fn finite_body_candidate(
        self,
        piston: &GasPiston,
        input: PistonInput,
        body: UniformEnthalpyStepInput<'_>,
        controls: PistonSolve,
        exchange: HeatExchangeSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<
        ThermalPistonCandidate<FilmGasConvectionSample>,
        ThermalPistonError<MaterialConvectionError>,
    > {
        if cancelled() {
            return Err(ThermalPistonError::Cancelled);
        }
        let area = self.law.area().map_err(ThermalPistonError::Interface)?;
        if !body.surface_area_m2.is_finite()
            || body.surface_area_m2 <= 0.0
            || area > body.surface_area_m2
        {
            return Err(ThermalPistonError::Invalid(
                "film convective patch exceeds the finite body's exposed area",
            ));
        }
        piston.thermal_candidate(
            input,
            body,
            controls,
            exchange,
            |gas, wall| self.resolve(gas, wall),
            cancelled,
        )
    }
}
