//! Film-temperature properties on the original finite gas/body heat residual.
//! Film evaluation has no independent energy or mass; bulk gas keeps both.
mod source;
pub use source::{FilmGasSource, PengRobinsonFilmSource};

use super::super::{
    AidingVerticalPlateSettings, evaluate_aiding_vertical_plate_from_material,
    evaluate_natural_from_material, mixed,
};
use super::*;
use fs_material::fluid_transport::{FluidTransport, ThermalFluidTransport};
use fs_material::thermochemical_gas::ThermochemicalGasModel;
use fs_qty::Temperature;

/// The existing correlation/geometry, with an explicit mean-film convention.
#[derive(Clone, Copy, Debug)]
pub enum FilmGasHeatLaw {
    Forced(ForcedGasSurface),
    Natural(NaturalGasSurface),
    /// Laminar aiding mixed convection with explicitly bounded density departure.
    AidingVerticalPlate(AidingVerticalPlateGasSurface),
}

/// One isothermal vertical plate patch in a prescribed aiding free stream.
/// Its mean-film gas is a property query only; the bulk gas/body keep all storage.
/// The complete heat-solve bracket must remain in the aiding and density-valid
/// domain, including trial endpoints that the existing solver evaluates first.
#[derive(Clone, Copy, Debug)]
pub struct AidingVerticalPlateGasSurface {
    pub area_m2: f64,
    /// Plate height along the vertical free stream [m].
    pub characteristic_length: Length,
    pub settings: AidingVerticalPlateSettings,
}

impl AidingVerticalPlateGasSurface {
    fn validate(self) -> Result<(), MaterialConvectionError> {
        self.settings.validate()?;
        if !self.area_m2.is_finite() || self.area_m2 < 0.0
            || !self.characteristic_length.value().is_finite()
            || self.characteristic_length.value() <= 0.0
        {
            return Err(invalid("aiding mixed gas convection needs nonnegative area and positive plate height"));
        }
        Ok(())
    }
}

/// Source bank and constitutive boundary law for one finite gas/body patch.
/// The bank may contain more species than the current gas, or different initial
/// fractions. Each query reconstructs the ACTUAL active gas composition and
/// preserves its original source claims before evaluating film properties.
pub struct FilmGasSurface<'a, M: FilmGasSource = ThermochemicalGasModel> {
    pub model: &'a M,
    pub law: FilmGasHeatLaw,
}
impl<M: FilmGasSource> Copy for FilmGasSurface<'_, M> {}
impl<M: FilmGasSource> Clone for FilmGasSurface<'_, M> {
    fn clone(&self) -> Self { *self }
}

/// The actual bulk gas/body states and the separate property evaluation.
/// Bulk gas remains the thermodynamic/acoustic state; film never replaces it.
#[derive(Debug)]
pub struct FilmGasConvectionSample<S = ResolvedThermochemicalGasState> {
    pub interface: FiniteGasConvectionSample<S>,
    /// None only for an explicitly isolated zero-area patch.
    pub property_gas: Option<S>,
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
            Self::AidingVerticalPlate(surface) => {
                surface.validate()?;
                Ok(surface.area_m2)
            }
        }
    }
}

impl<M: FilmGasSource> FilmGasSurface<'_, M> {
    /// Re-query source-pinned properties at (T_wall + T_bulk)/2 and actual bulk
    /// pressure. Both signed heat flow and Rayleigh forcing retain the FULL
    /// wall-to-bulk temperature contrast. The film is not a replacement gas
    /// endpoint, a second parcel, or an invented finite-capacity reservoir.
    ///
    /// The original correlation must support this property convention and the
    /// declared geometry. The original ideal or explicit-branch real-gas owner
    /// retains frozen composition and its single-phase source domain. This is
    /// a property query with no boundary-layer or chemical heat solver.
    /// The callback works with existing powered/reactive/multi-boundary owners;
    /// their complete exposed-area checks and source budgets remain mandatory.
    pub fn resolve(
        self,
        gas: &M::State,
        body: EquilibriumPhaseState,
    ) -> Result<ConductanceSample<FilmGasConvectionSample<M::State>>, MaterialConvectionError> {
        let area_m2 = self.law.area()?;
        let wall = body.temperature_k();
        let bulk = gas.temperature_k();
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
            .evaluate_property_temperature(gas, Temperature::new(0.5 * wall + 0.5 * bulk))?;
        let expansion = || -> Result<f64, MaterialConvectionError> {
            let thermal = film.thermal_transport_state()?;
            if thermal.transport() != film.transport_state()? || thermal.temperature_k() != film.temperature_k() {
                return Err(invalid("film thermal and momentum views must describe the same material state"));
            }
            thermal.volumetric_expansion_per_k()
                .ok_or_else(|| invalid("film buoyancy requires source-owned volumetric expansion"))
        };
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
                        expansion()?,
                        Temperature::new(wall),
                        Temperature::new(bulk),
                        surface.characteristic_length,
                    )?)
                }
            }
            FilmGasHeatLaw::AidingVerticalPlate(surface) => {
                if wall == bulk && surface.settings.speed_up_m_s == 0.0 {
                    mixed::inputs(film.transport_state()?, surface.settings,
                        expansion()?, Temperature::new(wall),
                        Temperature::new(bulk), surface.characteristic_length)?;
                    None
                } else {
                    Some(evaluate_aiding_vertical_plate_from_material(
                        &film, surface.settings, expansion()?,
                        Temperature::new(wall), Temperature::new(bulk),
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
    pub fn finite_body_candidate<G: PistonGas<State = M::State>>(
        self,
        piston: &GasPiston<G>,
        input: PistonInput,
        body: UniformEnthalpyStepInput<'_>,
        controls: PistonSolve,
        exchange: HeatExchangeSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<
        ThermalPistonCandidate<FilmGasConvectionSample<M::State>, G>,
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
