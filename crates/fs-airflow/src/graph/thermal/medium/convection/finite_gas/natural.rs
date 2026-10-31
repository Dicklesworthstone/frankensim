//! Natural convection on the original finite gas/body heat residual.
//!
//! The shared geometry-specific correlation reads one actual source-resolved gas
//! state and both current temperatures. This adapter owns no energy store,
//! flow field, gas law or timestepper.

use super::super::{NaturalConvectionGeometry, evaluate_vertical_plate_from_thermal_material};
use super::{
    ConductanceSample, EquilibriumPhaseState, FiniteGasConvectionSample, GasConvectionEvaluation,
    GasPiston, HeatExchangeSolve, Length, MaterialConvectionError, PistonGas, PistonGasState,
    PistonInput, PistonSolve, ThermalPistonCandidate, ThermalPistonError, UniformEnthalpyStepInput,
    invalid,
};
use fs_material::fluid_transport::ThermalFluidTransport;
use fs_qty::Temperature;

/// Exposed isothermal natural-convection surface with declared geometry/gravity.
///
/// The original Churchill-Chu card consumes transport at the supplied BULK gas
/// state. This uniform-reservoir approximation requires the caller to admit
/// its temperature contrast; there is no independently resolved film state,
/// enclosure circulation, neighbouring-wall correction or measured-fidelity
/// claim. That bulk-property shortcut remains limited to vertical plates;
/// horizontal plates/cylinders and spheres require `FilmGasSurface` and its
/// source bank.
/// No flow speed, independent beta or material-name switch is introduced.
#[derive(Clone, Copy, Debug)]
pub struct NaturalGasSurface {
    /// Geometry and horizontal face orientation. Select the length below
    /// according to this geometry; exposed patch area is not a perimeter.
    pub geometry: NaturalConvectionGeometry,
    /// Actual exposed patch area [m2]; zero explicitly isolates this patch.
    pub area_m2: f64,
    /// Height, cylinder/sphere diameter or horizontal plate area/perimeter [m].
    pub characteristic_length: Length,
    /// Magnitude of gravitational acceleration [m/s2]. The sphere admits zero
    /// gravity through its sourced Nu=2 conductive limit; other geometries
    /// require strictly positive gravity and their original Rayleigh ranges.
    pub gravity_m_s2: f64,
}

impl NaturalGasSurface {
    pub(super) fn validate(self) -> Result<(), MaterialConvectionError> {
        if !self.area_m2.is_finite()
            || self.area_m2 < 0.0
            || !self.characteristic_length.value().is_finite()
            || self.characteristic_length.value() <= 0.0
            || !self.gravity_m_s2.is_finite()
            || self.gravity_m_s2 < 0.0
            || (self.gravity_m_s2 == 0.0 && self.geometry != NaturalConvectionGeometry::Sphere)
        {
            return Err(invalid(
                "finite gas natural convection requires nonnegative area, positive characteristic length and positive gravity (a sphere also admits zero gravity)",
            ));
        }
        Ok(())
    }

    /// Resolve natural convection at the ORIGINAL coupling owner's gas/body
    /// query states. Pressure, composition, density, Cp, viscosity and k remain
    /// bound to that same gas evaluation. Heating and cooling use the vertical
    /// plate's admitted symmetry; the signed heat remains with the outer owner.
    ///
    /// The source-owned thermal view supplies the local isobaric expansion
    /// `beta = -(1/rho) d(rho)/dT` at that same state. Ideal and real gases retain
    /// their own EOS derivative alongside their density and Cp. This local
    /// derivative does not establish a finite-temperature validity interval
    /// for the bulk-property approximation.
    ///
    /// Exact equal temperatures have exactly zero transfer. Return G=0 and no
    /// correlation evaluation for that state, retaining both original sources;
    /// this does not claim the Ra=0 heat-transfer coefficient. Any nonzero
    /// contrast must pass the original correlation's Rayleigh-domain admission.
    /// Zero area likewise returns explicit isolation. Nonzero subnormal heat
    /// response is not silently clipped to either of these exact limits.
    pub fn resolve<S: PistonGasState + ThermalFluidTransport>(
        self,
        gas: &S,
        body: EquilibriumPhaseState,
    ) -> Result<ConductanceSample<FiniteGasConvectionSample<S>>, MaterialConvectionError> {
        self.validate()?;
        let wall_temperature_k = body.temperature_k();
        let gas_temperature_k = gas.temperature_k();
        if !wall_temperature_k.is_finite() || wall_temperature_k <= 0.0 {
            return Err(invalid(
                "finite gas natural convection needs a positive body temperature",
            ));
        }
        if self.area_m2 == 0.0 || wall_temperature_k == gas_temperature_k {
            return Ok(ConductanceSample {
                conductance_w_k: 0.0,
                source: FiniteGasConvectionSample {
                    gas: gas.clone(),
                    body,
                    area_m2: self.area_m2,
                    evaluation: None,
                },
            });
        }
        if self.geometry != NaturalConvectionGeometry::VerticalPlate {
            return Err(invalid(
                "horizontal and spherical natural convection require source-resolved FilmGasSurface properties",
            ));
        }
        let evaluated = evaluate_vertical_plate_from_thermal_material(
            gas,
            self.gravity_m_s2,
            Temperature::new(wall_temperature_k),
            self.characteristic_length,
        )?;
        let (transport, nusselt, coefficient) = evaluated.into_parts();
        let conductance = self.area_m2 * coefficient.value.value();
        if !conductance.is_finite() || conductance <= 0.0 {
            return Err(invalid(
                "finite gas natural conductance overflowed or vanished",
            ));
        }
        Ok(ConductanceSample {
            conductance_w_k: conductance,
            source: FiniteGasConvectionSample {
                gas: gas.clone(),
                body,
                area_m2: self.area_m2,
                evaluation: Some(GasConvectionEvaluation {
                    transport,
                    nusselt,
                    coefficient,
                }),
            },
        })
    }

    /// Propose buoyant heat exchange, finite body enthalpy and gas/piston motion
    /// together through the original thermal residual. Neither input mutates;
    /// publish both returned endpoints only after the enclosing step accepts.
    ///
    /// The patch must fit the declared body boundary and have one heat owner.
    /// The original uniform-body and negligible body pV-work assumptions remain;
    /// the interface does not infer body conductivity or a Biot certificate.
    /// Gravity drives the correlation, not a resolved momentum field. There is
    /// no added gravitational storage, viscous heat or radiation path.
    #[allow(clippy::too_many_arguments)]
    pub fn finite_body_candidate<G: PistonGas>(
        self,
        piston: &GasPiston<G>,
        input: PistonInput,
        body: UniformEnthalpyStepInput<'_>,
        controls: PistonSolve,
        exchange: HeatExchangeSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<
        ThermalPistonCandidate<FiniteGasConvectionSample<G::State>, G>,
        ThermalPistonError<MaterialConvectionError>,
    >
    where
        G::State: ThermalFluidTransport,
    {
        if cancelled() {
            return Err(ThermalPistonError::Cancelled);
        }
        self.validate().map_err(ThermalPistonError::Interface)?;
        if !body.surface_area_m2.is_finite()
            || body.surface_area_m2 <= 0.0
            || self.area_m2 > body.surface_area_m2
        {
            return Err(ThermalPistonError::Invalid(
                "natural convective patch exceeds the finite body's exposed area",
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
