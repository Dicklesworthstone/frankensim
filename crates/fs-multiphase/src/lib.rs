#![forbid(unsafe_code)]
#![warn(missing_docs)]

//! L3 (`FLUX`) multiphase wall boiling and filmwise condensation constitutive laws.
//!
//! Provides:
//! - [`boiling`]: Full five-regime Nukiyama pool-boiling curve with Rohsenow
//!   nucleate boiling, exact implicit Schluender/Thome-Shakir binary mixture
//!   suppression driven by `fs-thermochem`'s NRTL `Delta T_bp = T_dew - T_bub`,
//!   Lienhard-Dhir/Zuber CHF (`C_CHF = 0.149`), `C^1`-smooth monotone cubic
//!   Hermite transition boiling, Berenson/Bromley + Henry film boiling /
//!   Leidenfrost minimum, exact analytical `dq''/dT_w`, and explicit refusal of
//!   Nukiyama boiling on pure triglyceride oils.
//! - [`condensation`]: Nusselt laminar, wavy-laminar, and turbulent filmwise
//!   condensation (vertical plate, inclined plate, horizontal cylinder) with
//!   Colburn-Hougen non-condensable air boundary-layer degradation.
//! - [`BoilingCondensationInterfaceV1`]: Unified single-source surface
//!   phase-change evaluator dispatching strictly on `T_w > T_sat` (boiling) or
//!   `T_w < T_dew` (condensation) with zero double-counting against sub-boiling
//!   convection.
//! - [`vessel`]: Stirred multiphase liquid/vapor/oil/emulsion vessel
//!   [`fs_ports::component::Component`] (`MultiphaseVesselComponent`) unifying
//!   stockpots, saucepans, skillets, pressure cookers, and commercial deep fryers
//!   with dryout regularization, NRTL wine deglazing & fond dissolution, Carreau
//!   sauce rheology, Bikerman boil-over foam, butter water-flash to *beurre
//!   noisette* browning, and deep-fryer cold-zone/hydrolysis/TPM aging.
//! - [`cycle`]: Closed-loop Péneloux-translated Peng--Robinson vapor-compression
//!   refrigeration cycle components (`RefrigerationCompressorComponent`,
//!   `CondenserEvaporatorHXComponent`, `ExpansionThrottleComponent`, and
//!   `IceHarvesterComponent`) across `R-290`, `R-134a`, `R-404A`, `R-448A`, and
//!   subcritical/transcritical `R-744`.

/// Full five-regime Nukiyama pool boiling, Schluender-Thome mixture suppression,
/// Lienhard-Dhir CHF, `C^1` transition boiling, and Bromley/Leidenfrost film boiling.
pub mod boiling;

/// Nusselt filmwise condensation and Colburn-Hougen non-condensable gas degradation.
pub mod condensation;

/// Closed-loop Peng--Robinson vapor-compression cycle components and air-cooled ice harvester.
pub mod cycle;

/// Stirred multiphase liquid/vapor/oil/emulsion vessel [`fs_ports::component::Component`].
pub mod vessel;

pub use boiling::{
    BoilingCurveV1, BoilingErrorV1, BoilingEvaluationV1, BoilingFluidKindV1,
    BoilingFluidPropertiesV1, BoilingRegimeBoundariesV1, BoilingRegimeV1, BoilingStepAdviceV1,
    BoilingSurfaceCardV1, LIENHARD_DHIR_FLAT_PLATE_C_CHF_V1, LeidenfrostModeV1,
    ROHSENOW_SURFACE_LAW_ID_V1, STANDARD_GRAVITY_M_S2, ZUBER_INFINITE_PLATE_C_CHF_V1,
};
pub use condensation::{
    CondensationErrorV1, CondensationEvaluationV1, CondensationFluidPropertiesV1,
    CondensationGeometryV1, CondensationModelV1, CondensationRegimeV1, NonCondensableGasV1,
};
pub use cycle::{
    COMPRESSOR_PROBE_COUNT, ClosedCycleSummaryV1, CompressionResultV1, CompressorStateV1,
    CondenserEvaporatorHXComponent, CondenserEvaporatorHXConfig, CondenserEvaporatorHXStateV1,
    CycleErrorV1, CycleOperatingPointV1, ExpansionThrottleComponent, ExpansionThrottleConfig,
    ExpansionThrottleStateV1, HX_PROBE_COUNT, HxModeV1, ICE_HARVESTER_PROBE_COUNT,
    IIR_REFERENCE_TEMPERATURE_K, IIR_SAT_LIQUID_ENTHALPY_J_KG, IIR_SAT_LIQUID_ENTROPY_J_KG_K,
    IceHarvesterComponent, IceHarvesterConfig, IceHarvesterStateV1, PR_REFRIGERANT_LAW_ID_V1,
    PR_REFRIGERANT_LAW_VERSION_V1, REFRIGERANT_STREAM_WIDTH, RefrigerantModelV1,
    RefrigerantStatePointV1, RefrigerationCompressorComponent, RefrigerationCompressorConfig,
    SaturationEnvelopeV1, THROTTLE_PROBE_COUNT, ThrottleResultV1,
};
pub use fs_thermochem::caloric::proximate::JANAF_FORMATION_DATUM_ID_V1;
pub use vessel::{
    COOKING_OIL_LAW_ID_V1, COOKING_OIL_LAW_VERSION_V1, LIQUID_STREAM_WIDTH, MOLAR_MASS_ETOH_KG_MOL,
    MOLAR_MASS_FFA_KG_MOL, MOLAR_MASS_GLYCEROL_KG_MOL, MOLAR_MASS_H2O_KG_MOL,
    MOLAR_MASS_NACL_KG_MOL, MOLAR_MASS_SUCROSE_KG_MOL, MOLAR_MASS_TRIGLYCERIDE_KG_MOL,
    MultiphaseVesselComponent, MultiphaseVesselConfig, MultiphaseVesselState, OilCardV1,
    OpenVentSpecV1, PROBE_FFA_PCT, PROBE_FOAM_H_M, PROBE_FOND_ADHERED_KG, PROBE_M_LIQ_KG,
    PROBE_SAFETY_OR_BOILOVER_FLAG, PROBE_T_COLDZONE_K, PROBE_T_LIQ_K, PROBE_T_SMOKE_K,
    PROBE_TPM_PCT, PROBE_VISC_PA_S, PROBE_X_ETOH, PROBE_X_SOLIDS, VAPOR_STREAM_WIDTH,
    VESSEL_PROBE_COUNT, VesselErrorV1, VesselInitialStateV1, VesselSafetyFlagsV1,
};

use fs_qty::Temperature;

/// Active wall phase-change mode selected by [`BoilingCondensationInterfaceV1`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SurfacePhaseChangeModeV1 {
    /// Wall temperature exceeds local saturation/bubble temperature (`T_w > T_sat`).
    Boiling(BoilingRegimeV1),
    /// Wall temperature lies below local vapor dew-point temperature (`T_w < T_dew`).
    Condensation(CondensationRegimeV1),
    /// Wall temperature lies between `T_dew` and `T_sat` (or equal to boundary):
    /// no wall phase change; sensible heat transfer is owned by `fs-convection`.
    NoPhaseChange,
}

/// Unified single-source wall phase-change evaluation (signed heat flux positive
/// FROM wall INTO fluid).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SurfacePhaseChangeEvaluationV1 {
    /// Active phase-change mode.
    pub mode: SurfacePhaseChangeModeV1,
    /// Signed heat flux from wall into fluid [W/m^2] (`> 0` for boiling, `< 0` for condensation).
    pub signed_heat_flux_w_m2: f64,
    /// Analytical tangent `d(signed_heat_flux)/dT_w` [W/(m^2 K)].
    pub dq_dtw_w_m2_k: f64,
    /// Signed phase-change mass flux [kg/(m^2 s)] (`> 0` evaporation, `< 0` condensation).
    pub signed_vaporization_flux_kg_m2_s: f64,
    /// Vapor-phase specific enthalpy on `janaf-298.15K-formation` datum [J/kg].
    pub vapor_enthalpy_janaf_j_kg: f64,
    /// Liquid-phase specific enthalpy on `janaf-298.15K-formation` datum [J/kg].
    pub liquid_enthalpy_janaf_j_kg: f64,
}

/// Unified surface two-phase heat- and mass-transfer interface pairing a prepared
/// boiling curve (`T_w > T_sat`) with a condensation model (`T_w < T_dew`).
///
/// Guarantees single-source ownership: boiling and condensation never overlap
/// because `T_dew <= T_sat(P)` for every vapor-air mixture at total pressure `P`.
#[derive(Clone, Debug, PartialEq)]
pub struct BoilingCondensationInterfaceV1 {
    boiling: BoilingCurveV1,
    condensation: CondensationModelV1,
}

impl BoilingCondensationInterfaceV1 {
    /// Construct a unified boiling/condensation interface from prepared models.
    ///
    /// # Errors
    /// Returns [`BoilingErrorV1::InvalidParameter`] if condensation dew temperature
    /// exceeds boiling saturation temperature (which would violate single-source
    /// non-overlapping ownership).
    pub fn new(
        boiling: BoilingCurveV1,
        condensation: CondensationModelV1,
    ) -> Result<Self, BoilingErrorV1> {
        if condensation.dew_temperature().value() > boiling.saturation_temperature().value() + 1e-6
        {
            return Err(BoilingErrorV1::InvalidParameter(
                "condensation dew temperature cannot exceed boiling saturation temperature",
            ));
        }
        Ok(Self {
            boiling,
            condensation,
        })
    }

    /// Borrow the underlying prepared Nukiyama boiling curve.
    #[must_use]
    pub const fn boiling(&self) -> &BoilingCurveV1 {
        &self.boiling
    }

    /// Borrow the underlying prepared Nusselt/Colburn-Hougen condensation model.
    #[must_use]
    pub const fn condensation(&self) -> &CondensationModelV1 {
        &self.condensation
    }

    /// Evaluate signed wall phase-change heat flux (`> 0` wall-to-fluid boiling,
    /// `< 0` fluid-to-wall condensation) and exact tangent `dq''/dT_w` at `wall_temperature`.
    pub fn evaluate(
        &self,
        wall_temperature: Temperature,
    ) -> Result<SurfacePhaseChangeEvaluationV1, BoilingErrorV1> {
        let tw = wall_temperature.value();
        if !tw.is_finite() || tw <= 0.0 {
            return Err(BoilingErrorV1::InvalidTemperature { temperature_k: tw });
        }
        let t_sat = self.boiling.saturation_temperature().value();
        let t_dew = self.condensation.dew_temperature().value();

        if tw > t_sat {
            let b = self.boiling.evaluate(wall_temperature)?;
            return Ok(SurfacePhaseChangeEvaluationV1 {
                mode: SurfacePhaseChangeModeV1::Boiling(b.regime),
                signed_heat_flux_w_m2: b.heat_flux_w_m2,
                dq_dtw_w_m2_k: b.dq_dtw_w_m2_k,
                signed_vaporization_flux_kg_m2_s: b.vaporization_mass_flux_kg_m2_s,
                vapor_enthalpy_janaf_j_kg: b.vapor_enthalpy_janaf_j_kg,
                liquid_enthalpy_janaf_j_kg: b.liquid_enthalpy_janaf_j_kg,
            });
        }

        if tw < t_dew {
            let c = self.condensation.evaluate(wall_temperature).map_err(|e| {
                BoilingErrorV1::Unrepresentable(match e {
                    CondensationErrorV1::InvalidTemperature { .. } => {
                        "condensation wall temperature invalid"
                    }
                    CondensationErrorV1::ConvergenceFailed { .. } => {
                        "Colburn-Hougen interface balance did not converge"
                    }
                    _ => "condensation evaluation failed",
                })
            })?;
            // Condensation delivers heat FROM fluid INTO wall, so wall-to-fluid flux is negative.
            // d(-q_cond)/dT_w = -d(q_cond)/dT_w > 0 (raising T_w reduces condensation flux).
            return Ok(SurfacePhaseChangeEvaluationV1 {
                mode: SurfacePhaseChangeModeV1::Condensation(c.regime),
                signed_heat_flux_w_m2: -c.heat_flux_w_m2,
                dq_dtw_w_m2_k: -c.dq_dtw_w_m2_k,
                signed_vaporization_flux_kg_m2_s: -c.condensation_mass_flux_kg_m2_s,
                vapor_enthalpy_janaf_j_kg: c.vapor_enthalpy_janaf_j_kg,
                liquid_enthalpy_janaf_j_kg: c.liquid_enthalpy_janaf_j_kg,
            });
        }

        let props = self.boiling.fluid_properties();
        Ok(SurfacePhaseChangeEvaluationV1 {
            mode: SurfacePhaseChangeModeV1::NoPhaseChange,
            signed_heat_flux_w_m2: 0.0,
            dq_dtw_w_m2_k: 0.0,
            signed_vaporization_flux_kg_m2_s: 0.0,
            vapor_enthalpy_janaf_j_kg: props.vapor_enthalpy_janaf_j_kg,
            liquid_enthalpy_janaf_j_kg: props.liquid_enthalpy_janaf_j_kg,
        })
    }
}
