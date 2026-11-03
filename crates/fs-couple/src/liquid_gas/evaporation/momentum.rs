//! Mean-drift phase transport between native liquid and one finite gas parcel.
//!
//! Condensing mass arrives with the retained gas mean velocity. Evaporating
//! mass leaves with the actual native liquid velocity. Perfectly inelastic
//! receiver mixing deposits only the corresponding relative kinetic loss in
//! the original receiver's caloric owner. This is a declared mean-drift
//! closure, not a resolved vapor jet, molecular half-range flux or recoil law.
//! Pressure traction remains the separate original native boundary operation.
//!
//! Gas EOS storage stays internal energy. The explicit gas bulk momentum below
//! is retained alongside it; kinetic energy is derived, never added to EOS
//! internal energy or counted again in vapor enthalpy.

use super::*;
use fs_lbm::material::free_surface::enthalpy::FreeSurfaceMassInflow;
use fs_material::thermochemical_gas::energy::{GasEnergyError, GasHeatStep};

mod kinematics;
use kinematics::kinetic_energy;

/// Original finite gas plus an explicitly initialized bulk momentum.
///
/// This is a homogeneous mean-drift source reservoir, not a resolved gas-flow
/// solver. Its fixed connected volume, source composition and internal caloric
/// state retain exactly the meaning of ClosedGasParcel. Scenario authors must
/// supply the initial velocity in the native liquid's Cartesian frame.
#[derive(Clone, Debug, PartialEq)]
pub struct MovingGasParcel {
    parcel: ClosedGasParcel,
    momentum_kg_m_s: [f64; 3],
}

impl MovingGasParcel {
    pub fn new(
        parcel: ClosedGasParcel,
        velocity_m_s: [f64; 3],
    ) -> Result<Self, EvaporatingLiquidGasError> {
        let momentum = velocity_m_s.map(|u| parcel.mass_kg() * u);
        if velocity_m_s.iter().any(|u| !u.is_finite())
            || (0..3).any(|d| velocity_m_s[d] != 0.0 && momentum[d] == 0.0)
        {
            return Err(EvaporatingLiquidGasError::Invalid(
                "finite representable initial gas drift required",
            ));
        }
        Self::from_momentum(parcel, momentum)
    }

    fn from_momentum(
        parcel: ClosedGasParcel,
        momentum_kg_m_s: [f64; 3],
    ) -> Result<Self, EvaporatingLiquidGasError> {
        let velocity = velocity(parcel.mass_kg(), momentum_kg_m_s)?;
        if !kinetic_energy(parcel.mass_kg(), velocity).is_finite() {
            return Err(EvaporatingLiquidGasError::Invalid(
                "gas bulk kinetic energy is unrepresentable",
            ));
        }
        Ok(Self {
            parcel,
            momentum_kg_m_s,
        })
    }

    pub const fn parcel(&self) -> &ClosedGasParcel {
        &self.parcel
    }
    pub const fn momentum_kg_m_s(&self) -> [f64; 3] {
        self.momentum_kg_m_s
    }
    pub fn velocity_m_s(&self) -> [f64; 3] {
        self.momentum_kg_m_s.map(|p| p / self.parcel.mass_kg())
    }
    pub fn kinetic_energy_j(&self) -> f64 {
        kinetic_energy(self.parcel.mass_kg(), self.velocity_m_s())
    }
    pub fn into_parts(self) -> (ClosedGasParcel, [f64; 3]) {
        (self.parcel, self.momentum_kg_m_s)
    }
}

/// Closed mechanical account of the phase SOURCE before native flow.
///
/// The following native flow step can exchange momentum with pressure and
/// stationary walls and can dissipate resolved kinetic energy numerically.
/// Those effects are not phase mixing heat and are not included here.
#[derive(Clone, Debug)]
pub struct PhaseMomentumReport {
    /// The one physical phase/flow interval [s], not an extra source timestep.
    pub duration_s: f64,
    /// Mean velocity of the actual donor. Zero transfer uses the liquid value.
    pub donor_velocity_m_s: [f64; 3],
    pub gas_momentum_before_kg_m_s: [f64; 3],
    pub gas_momentum_after_kg_m_s: [f64; 3],
    pub gas_velocity_before_m_s: [f64; 3],
    pub gas_velocity_after_m_s: [f64; 3],
    /// Actual original native population change, including material refresh.
    pub liquid_kinetic_energy_change_j: f64,
    /// Actual gas kinetic change computed from increments, not large totals.
    pub gas_kinetic_energy_change_j: f64,
    /// Physical receiver loss already stored in original native liquid h [J].
    pub liquid_mixing_heat_j: f64,
    /// Physical receiver loss requested from original gas calorics [J].
    pub gas_mixing_heat_j: f64,
    /// Actual original source-caloric update; None when no gas mixing occurs.
    pub gas_heat: Option<GasHeatStep>,
    /// Actual liquid plus gas momentum change [kg m/s], component by component.
    pub momentum_residual_kg_m_s: [f64; 3],
    /// Sum of absolute components, gated by controls.native momentum tolerance.
    pub absolute_momentum_residual_kg_m_s: f64,
    /// Actual liquid plus gas kinetic change plus physical liquid mixing heat
    /// and actual gas heat storage [J]. Caloric representation in the liquid
    /// is accounted separately by native.native_phase_caloric_defect_j.
    pub energy_residual_j: f64,
}

#[derive(Clone, Debug)]
pub struct MovingGasEvaporationReport {
    /// Original source/flow reports. energy_residual_j additionally includes
    /// the closed source mechanical account below. Original phase endpoints
    /// precede receiver mixing; actual final states are in the candidate.
    pub native: EvaporatingLiquidGasReport,
    pub momentum: PhaseMomentumReport,
}

pub struct MovingGasEvaporationCandidate<'s, 'm, 'catalog> {
    liquid: MaterialFreeSurfaceEnthalpy3<'s, 'm>,
    gas: MovingGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    report: MovingGasEvaporationReport,
}
impl<'s, 'm, 'catalog> MovingGasEvaporationCandidate<'s, 'm, 'catalog> {
    pub const fn liquid(&self) -> &MaterialFreeSurfaceEnthalpy3<'s, 'm> {
        &self.liquid
    }
    pub const fn gas(&self) -> &MovingGasParcel {
        &self.gas
    }
    pub const fn catalog(&self) -> &'catalog ThermochemicalGasModel {
        self.catalog
    }
    pub const fn report(&self) -> &MovingGasEvaporationReport {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        MaterialFreeSurfaceEnthalpy3<'s, 'm>,
        MovingGasParcel,
        &'catalog ThermochemicalGasModel,
        MovingGasEvaporationReport,
    ) {
        (self.liquid, self.gas, self.catalog, self.report)
    }
}

/// Advance native phase, reciprocal donor momentum and one original flow step.
///
/// Optional sensible film heat uses the original coefficient/solve pair;
/// optional pressure uses the original maximum density-excursion bound. Both
/// preserve the original phase entry's clock, source, geometry, capacity and
/// connected fixed-volume assumptions. The final mixing-heated gas EOS supplies
/// optional pressure loading. Pressure forces are not an additional phase
/// impulse, and this entry does not claim a closed whole-flow FSI balance.
///
/// The original Stefan-film and latent/enthalpy source solve chooses the
/// represented transfer. The native local-velocity exchange handles
/// evaporation, including zero-flux heat-only intervals; original native
/// complete mixing handles condensation at the retained upstream gas drift.
/// Relative receiver losses heat the original liquid or gas once. No velocity
/// is inferred from a saturation pressure, mass flux or species identity.
///
/// Mixing and phase are a first-order split over the SAME retained interval;
/// refining that original interval resolves this splitting error. A source,
/// numerical, native-speed, conservation or cancellation refusal preserves
/// both accepted owners and their clock.
#[allow(clippy::too_many_arguments)]
pub fn advance_liquid_with_moving_gas<'s, 'm, 'catalog>(
    liquid: &MaterialFreeSurfaceEnthalpy3<'s, 'm>,
    gas: &MovingGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    surface: LiquidGasEvaporationSurface<'_>,
    controls: EvaporatingLiquidGasSolve,
    heat_controls: Option<(f64, GasFilmHeatSolve)>,
    maximum_density_excursion: Option<f64>,
    gate: &CancelGate,
) -> Result<MovingGasEvaporationCandidate<'s, 'm, 'catalog>, EvaporatingLiquidGasError> {
    let (candidate, motion) = super::advance_native_with_momentum(
        liquid,
        gas.parcel(),
        catalog,
        surface,
        controls,
        heat_controls,
        maximum_density_excursion,
        Some(gas.momentum_kg_m_s()),
        gate,
    )?;
    let momentum = motion.ok_or(EvaporatingLiquidGasError::Invalid(
        "moving gas phase source did not retain its momentum account",
    ))?;
    let (liquid, parcel, catalog, native) = candidate.into_parts();
    let gas = MovingGasParcel::from_momentum(parcel, momentum.gas_momentum_after_kg_m_s)?;
    poll(gate)?;
    Ok(MovingGasEvaporationCandidate {
        liquid,
        gas,
        catalog,
        report: MovingGasEvaporationReport { native, momentum },
    })
}

fn velocity(mass: f64, momentum: [f64; 3]) -> Result<[f64; 3], EvaporatingLiquidGasError> {
    kinematics::velocity(mass, momentum).ok_or(EvaporatingLiquidGasError::Invalid(
        "gas mass or drift is unrepresentable",
    ))
}

fn gas_error(error: GasEnergyError) -> EvaporatingLiquidGasError {
    if error == GasEnergyError::Cancelled {
        EvaporatingLiquidGasError::Cancelled
    } else {
        GasEvaporationError::from(error).into()
    }
}

pub(super) fn admit_mixed_endpoints(
    liquid: &MaterialFreeSurfaceEnthalpy3<'_, '_>,
    gas: &ClosedGasParcel,
    catalog: &ThermochemicalGasModel,
    surface: LiquidGasEvaporationSurface<'_>,
    controls: EvaporatingLiquidGasSolve,
    phase: &GasEvaporationStep,
    liquid_h: f64,
) -> Result<(), EvaporatingLiquidGasError> {
    let temperature = liquid.liquid_temperature_at_specific_enthalpy_k(liquid_h)?;
    let properties = fs_thermal::evaporation::properties(surface.law, temperature)
        .map_err(GasEvaporationError::from)?;
    if properties.saturation_pressure_pa >= gas.state().state().pressure {
        return Err(GasEvaporationError::Surface(EvaporationError::Boiling {
            saturation_pressure_pa: properties.saturation_pressure_pa,
            total_pressure_pa: gas.state().state().pressure,
        })
        .into());
    }
    let residual = catalog
        .species_specific_enthalpy_j_kg(&surface.binding.species, temperature)
        .map_err(GasEvaporationError::from)?
        - (liquid_h + surface.binding.liquid_enthalpy_offset_j_kg)
        - properties.latent_heat_j_kg;
    if !residual.is_finite() || residual.abs() > controls.phase.latent_compatibility_tolerance_j_kg
    {
        return Err(GasEvaporationError::CaloricMismatch {
            temperature_k: temperature,
            residual_j_kg: residual,
            tolerance_j_kg: controls.phase.latent_compatibility_tolerance_j_kg,
        }
        .into());
    }
    if (temperature - phase.surface_temperature_k[0]).abs()
        > controls.phase.maximum_temperature_increment_k
        || (gas.state().state().temperature - phase.gas_temperature_k[0]).abs()
            > controls.phase.maximum_temperature_increment_k
    {
        return Err(GasEvaporationError::StepSizeLimit {
            maximum_duration_s: 0.5 * phase.duration_s,
            reason: "phase plus receiver mixing exceeds the original temperature increment bound",
        }
        .into());
    }
    Ok(())
}

pub(super) fn exchange(
    liquid: &mut MaterialFreeSurfaceEnthalpy3<'_, '_>,
    original_gas: &ClosedGasParcel,
    gas: &mut ClosedGasParcel,
    original_momentum: [f64; 3],
    transfer: FreeSurfaceMassExchange,
    controls: EvaporatingLiquidGasSolve,
    gate: &CancelGate,
) -> Result<(FreeSurfaceMassExchangeReport, PhaseMomentumReport), EvaporatingLiquidGasError> {
    type Error = EvaporatingLiquidGasError;
    poll(gate)?;
    let gas_velocity = velocity(original_gas.mass_kg(), original_momentum)?;
    let [x, y, z] = transfer.coordinates;
    let liquid_velocity = liquid
        .flow()
        .sample_at_si(x, y, z)
        .map_err(FreeSurfaceEnthalpyError::from)?
        .ok_or(Error::Invalid("phase donor has no native liquid velocity"))?
        .velocity_m_s;
    let (native, liquid_kinetic, liquid_heat, donor) = if transfer.mass_kg > 0.0 {
        let report = liquid.mix_interface_inflow(
            &[FreeSurfaceMassInflow {
                exchange: transfer,
                velocity_m_s: gas_velocity,
            }],
            controls.native,
            gate,
        )?;
        (
            report.exchange,
            report.kinetic_energy_change_j,
            report.mixing_heat_j,
            gas_velocity,
        )
    } else {
        let report =
            liquid.exchange_interface_mass_with_kinetic(&[transfer], controls.native, gate)?;
        (
            report.exchange,
            report.kinetic_energy_change_j,
            0.0,
            liquid_velocity,
        )
    };
    let gas_dm = -native.mass_change_kg;
    let motion = kinematics::advance(
        original_gas.mass_kg(),
        gas.mass_kg(),
        original_momentum,
        gas_dm,
        donor,
    )
    .ok_or(Error::Invalid(
        "gas phase momentum or physical mixing heat is unrepresentable",
    ))?;
    let next_momentum = motion.momentum;
    let next_velocity = motion.velocity;
    let gas_heat_j = motion.mixing_heat_j;
    let gas_heat = if gas_heat_j != 0.0 {
        Some(
            gas.add_heat(gas_heat_j, controls.phase.gas.gas, || gate.is_requested())
                .map_err(gas_error)?,
        )
    } else {
        None
    };
    let gas_kinetic = motion.kinetic_change_j;
    let residual = core::array::from_fn(|d| {
        native.momentum_change_kg_m_s[d] + (next_momentum[d] - original_momentum[d])
    });
    let absolute_momentum = residual.iter().map(|p| p.abs()).sum::<f64>();
    let energy = liquid_kinetic
        + gas_kinetic
        + liquid_heat
        + gas_heat
            .as_ref()
            .map_or(0.0, |heat| heat.stored_energy_change_j);
    if ![
        liquid_kinetic,
        gas_kinetic,
        liquid_heat,
        absolute_momentum,
        energy,
    ]
    .iter()
    .all(|v| v.is_finite())
        || absolute_momentum > controls.native.momentum_tolerance_kg_m_s
    {
        return Err(Error::Invalid(
            "combined phase-source momentum exceeded its native budget",
        ));
    }
    if energy.abs() > controls.energy_tolerance_j {
        return Err(Error::EnergyBalance {
            residual_j: energy,
            tolerance_j: controls.energy_tolerance_j,
        });
    }
    poll(gate)?;
    Ok((
        native,
        PhaseMomentumReport {
            duration_s: liquid.flow().scaling().timestep_s(),
            donor_velocity_m_s: donor,
            gas_momentum_before_kg_m_s: original_momentum,
            gas_momentum_after_kg_m_s: next_momentum,
            gas_velocity_before_m_s: gas_velocity,
            gas_velocity_after_m_s: next_velocity,
            liquid_kinetic_energy_change_j: liquid_kinetic,
            gas_kinetic_energy_change_j: gas_kinetic,
            liquid_mixing_heat_j: liquid_heat,
            gas_mixing_heat_j: gas_heat_j,
            gas_heat,
            momentum_residual_kg_m_s: residual,
            absolute_momentum_residual_kg_m_s: absolute_momentum,
            energy_residual_j: energy,
        },
    ))
}
