//! Nusselt laminar, wavy-laminar, and turbulent filmwise condensation with
//! Colburn-Hougen non-condensable air boundary-layer resistance.
//!
//! Supports vertical plates, inclined plates, and horizontal cylinders, parameterized
//! by explicit ambient pressure `P`, vapor/air mixture composition, and gas-side
//! convective mass-transfer conductance. Owns wall condensation strictly when
//! `T_w < T_dew(P, w_air)` so there is zero double counting with sub-boiling
//! single-phase convection.
//!
//! All transcendentals route strictly through `fs_math::det`.

use crate::boiling::STANDARD_GRAVITY_M_S2;
use core::fmt;
use fs_math::det::{ln, pow, sin};
use fs_qty::{Pressure, Temperature};
use fs_thermochem::caloric::proximate::water_state_enthalpy_janaf;
use fs_thermochem::water::if97::regions::If97FluidState;
use fs_thermochem::water::if97::{
    If97Region, If97SaturationPair, If97State, saturation_pressure, saturation_temperature,
};
use fs_thermochem::water::transport::If97TransportState;

/// Molar mass of ordinary water vapor [kg/mol].
pub const WATER_MOLAR_MASS_KG_MOL_V1: f64 = 0.018_015_28;
/// Molar mass of standard dry air [kg/mol].
pub const DRY_AIR_MOLAR_MASS_KG_MOL_V1: f64 = 0.028_964_7;
/// Default steam-air Lewis number `Le = alpha / D_AB` (dimensionless).
pub const STEAM_AIR_LEWIS_NUMBER_V1: f64 = 0.85;
/// Laminar-to-wavy-laminar condensate film Reynolds threshold (Incropera Ch. 10.7).
pub const NUSSELT_WAVY_RE_THRESHOLD_V1: f64 = 30.0;
/// Wavy-laminar-to-turbulent condensate film Reynolds threshold (Incropera Ch. 10.7).
pub const NUSSELT_TURBULENT_RE_THRESHOLD_V1: f64 = 1_800.0;

/// Dedicated typed refusal for filmwise condensation setup and evaluation.
#[derive(Clone, Debug, PartialEq)]
pub enum CondensationErrorV1 {
    /// Total pressure is non-finite, non-positive, or outside vapor saturation domain.
    InvalidPressure {
        /// Rejected pressure [Pa].
        pressure_pa: f64,
    },
    /// Wall or bulk temperature is non-finite or non-positive.
    InvalidTemperature {
        /// Rejected temperature [K].
        temperature_k: f64,
    },
    /// Condensation surface geometry (height, length, diameter, or tilt angle) is invalid.
    InvalidGeometry(&'static str),
    /// Non-condensable air mass fraction is non-finite or outside `[0, 1]`.
    InvalidNonCondensableFraction {
        /// Rejected mass fraction.
        air_mass_fraction: f64,
    },
    /// Invalid fluid property or gas-side mass-transfer parameter.
    InvalidParameter(&'static str),
    /// Underlying `fs-thermochem` saturation or transport query failed.
    Thermochem(&'static str),
    /// Bounded Colburn-Hougen interface balance did not converge.
    ConvergenceFailed {
        /// Iterations executed.
        iterations: usize,
        /// Final heat-flux balance residual [W/m^2].
        residual_w_m2: f64,
    },
    /// Intermediate or derived quantity became non-finite or unrepresentable.
    Unrepresentable(&'static str),
}

impl fmt::Display for CondensationErrorV1 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidPressure { pressure_pa } => {
                write!(f, "condensation pressure {pressure_pa} Pa is invalid")
            }
            Self::InvalidTemperature { temperature_k } => {
                write!(f, "condensation temperature {temperature_k} K is invalid")
            }
            Self::InvalidGeometry(msg) => write!(f, "invalid condensation geometry: {msg}"),
            Self::InvalidNonCondensableFraction { air_mass_fraction } => write!(
                f,
                "non-condensable air mass fraction {air_mass_fraction} must be finite in [0, 1]"
            ),
            Self::InvalidParameter(msg) => write!(f, "invalid condensation parameter: {msg}"),
            Self::Thermochem(msg) => write!(f, "condensation thermochem lookup failed: {msg}"),
            Self::ConvergenceFailed {
                iterations,
                residual_w_m2,
            } => write!(
                f,
                "Colburn-Hougen balance failed after {iterations} iterations (residual {residual_w_m2:e} W/m^2)"
            ),
            Self::Unrepresentable(msg) => write!(f, "unrepresentable condensation quantity: {msg}"),
        }
    }
}

impl core::error::Error for CondensationErrorV1 {}

/// Condensing surface geometry for Nusselt filmwise condensation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum CondensationGeometryV1 {
    /// Vertical flat plate or vertical wall of height `height_m > 0` [m].
    VerticalPlate {
        /// Vertical condensate drainage height `L` [m].
        height_m: f64,
    },
    /// Inclined flat plate of length `length_m > 0` [m] tilted at
    /// `tilt_from_horizontal_rad` in `(0, pi/2]` rad (`pi/2` = vertical).
    InclinedPlate {
        /// Plate length along drainage direction `L` [m].
        length_m: f64,
        /// Inclination angle from horizontal [rad] in `(0, pi/2]`.
        tilt_from_horizontal_rad: f64,
    },
    /// Single horizontal circular cylinder / tube of outer diameter `diameter_m > 0` [m].
    HorizontalCylinder {
        /// Outer cylinder diameter `D` [m].
        diameter_m: f64,
    },
}

impl CondensationGeometryV1 {
    /// Validate geometry dimensions and angles fail-closed.
    pub fn validate(self) -> Result<(), CondensationErrorV1> {
        match self {
            Self::VerticalPlate { height_m } => {
                if !height_m.is_finite() || height_m <= 0.0 {
                    return Err(CondensationErrorV1::InvalidGeometry(
                        "vertical plate height_m must be finite and strictly positive",
                    ));
                }
            }
            Self::InclinedPlate {
                length_m,
                tilt_from_horizontal_rad,
            } => {
                if !length_m.is_finite() || length_m <= 0.0 {
                    return Err(CondensationErrorV1::InvalidGeometry(
                        "inclined plate length_m must be finite and strictly positive",
                    ));
                }
                if !tilt_from_horizontal_rad.is_finite()
                    || tilt_from_horizontal_rad <= 0.0
                    || tilt_from_horizontal_rad > core::f64::consts::FRAC_PI_2 + 1e-12
                {
                    return Err(CondensationErrorV1::InvalidGeometry(
                        "inclined plate angle must lie in (0, pi/2] rad",
                    ));
                }
            }
            Self::HorizontalCylinder { diameter_m } => {
                if !diameter_m.is_finite() || diameter_m <= 0.0 {
                    return Err(CondensationErrorV1::InvalidGeometry(
                        "horizontal cylinder diameter_m must be finite and strictly positive",
                    ));
                }
            }
        }
        Ok(())
    }
}

/// Active filmwise condensation regime.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CondensationRegimeV1 {
    /// Wall temperature `T_w >= T_dew` or pure dry air (`w_air == 1`): zero condensation.
    NoCondensation,
    /// Smooth laminar Nusselt condensate film (`Re_delta <= 30`, or horizontal cylinder).
    LaminarFilm,
    /// Wavy-laminar condensate film (`30 < Re_delta <= 1800`, Kutateladze correlation).
    WavyLaminarFilm,
    /// Turbulent condensate film (`Re_delta > 1800`, Labuntsov correlation).
    TurbulentFilm,
}

/// Explicit condensate/vapor properties at the condensing condition (for analytical
/// constant-property checks or non-water condensables).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CondensationFluidPropertiesV1 {
    /// Condensate liquid density `rho_l` [kg/m^3].
    pub liquid_density_kg_m3: f64,
    /// Vapor density `rho_v` [kg/m^3].
    pub vapor_density_kg_m3: f64,
    /// Condensate dynamic viscosity `mu_l` [Pa s].
    pub liquid_viscosity_pa_s: f64,
    /// Condensate thermal conductivity `k_l` [W/(m K)].
    pub liquid_conductivity_w_m_k: f64,
    /// Condensate isobaric specific heat `c_{p,l}` [J/(kg K)].
    pub liquid_cp_j_kg_k: f64,
    /// Latent heat of condensation `h_fg` [J/kg].
    pub latent_heat_j_kg: f64,
    /// Vapor specific enthalpy on `janaf-298.15K-formation` datum [J/kg].
    pub vapor_enthalpy_janaf_j_kg: f64,
    /// Condensate liquid specific enthalpy on `janaf-298.15K-formation` datum [J/kg].
    pub liquid_enthalpy_janaf_j_kg: f64,
}

impl CondensationFluidPropertiesV1 {
    /// Validate positivity and finiteness of all condensation properties.
    pub fn validate(&self) -> Result<(), CondensationErrorV1> {
        let positives = [
            self.liquid_density_kg_m3,
            self.vapor_density_kg_m3,
            self.liquid_viscosity_pa_s,
            self.liquid_conductivity_w_m_k,
            self.liquid_cp_j_kg_k,
            self.latent_heat_j_kg,
        ];
        if !positives.into_iter().all(|v| v.is_finite() && v > 0.0) {
            return Err(CondensationErrorV1::InvalidParameter(
                "condensation fluid properties must be finite and strictly positive",
            ));
        }
        if self.liquid_density_kg_m3 <= self.vapor_density_kg_m3 {
            return Err(CondensationErrorV1::InvalidParameter(
                "condensate density must strictly exceed vapor density",
            ));
        }
        if !self.vapor_enthalpy_janaf_j_kg.is_finite()
            || !self.liquid_enthalpy_janaf_j_kg.is_finite()
        {
            return Err(CondensationErrorV1::InvalidParameter(
                "JANAF condensation enthalpies must be finite",
            ));
        }
        Ok(())
    }
}

/// Bulk non-condensable air specification for Colburn-Hougen boundary-layer degradation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NonCondensableGasV1 {
    /// Bulk non-condensable dry-air mass fraction `w_air` in `[0, 1]` (`0` = 100% pure steam).
    pub air_mass_fraction: f64,
    /// Bulk gas-mixture temperature `T_bulk` [K] (`>= T_dew`).
    pub bulk_temperature_k: f64,
    /// Gas-side sensible convective heat-transfer coefficient `h_g` [W/(m^2 K)]
    /// (e.g., `25.0` for forced combi-oven fan, `8.0` for insulated food carrier lid headspace).
    pub gas_convective_htc_w_m2_k: f64,
    /// Bulk gas-mixture isobaric specific heat `c_{p,g}` [J/(kg K)] (default `~1050-1900`).
    pub gas_cp_j_kg_k: f64,
    /// Steam-air Lewis number `Le` (default `0.85`).
    pub lewis_number: f64,
}

impl NonCondensableGasV1 {
    /// Pure saturated steam (`w_air = 0.0`), zero non-condensable diffusion resistance.
    #[must_use]
    pub const fn pure_steam(saturation_temperature_k: f64) -> Self {
        Self {
            air_mass_fraction: 0.0,
            bulk_temperature_k: saturation_temperature_k,
            gas_convective_htc_w_m2_k: 25.0,
            gas_cp_j_kg_k: 2_050.0,
            lewis_number: STEAM_AIR_LEWIS_NUMBER_V1,
        }
    }

    /// Steam-air mixture typical of a combi-oven cavity or hot-holding insulated food carrier headspace.
    pub fn steam_air(
        air_mass_fraction: f64,
        bulk_temperature: Temperature,
        gas_convective_htc_w_m2_k: f64,
    ) -> Result<Self, CondensationErrorV1> {
        if !air_mass_fraction.is_finite() || !(0.0..=1.0).contains(&air_mass_fraction) {
            return Err(CondensationErrorV1::InvalidNonCondensableFraction { air_mass_fraction });
        }
        let t_b = bulk_temperature.value();
        if !t_b.is_finite() || t_b <= 0.0 {
            return Err(CondensationErrorV1::InvalidTemperature { temperature_k: t_b });
        }
        if !gas_convective_htc_w_m2_k.is_finite() || gas_convective_htc_w_m2_k <= 0.0 {
            return Err(CondensationErrorV1::InvalidParameter(
                "gas_convective_htc_w_m2_k must be finite and positive",
            ));
        }
        // Mass-weighted dry-air (~1010 J/kg-K) and water-vapor (~2000 J/kg-K) specific heat.
        let cp_g = air_mass_fraction * 1_010.0 + (1.0 - air_mass_fraction) * 2_000.0;
        Ok(Self {
            air_mass_fraction,
            bulk_temperature_k: t_b,
            gas_convective_htc_w_m2_k,
            gas_cp_j_kg_k: cp_g,
            lewis_number: STEAM_AIR_LEWIS_NUMBER_V1,
        })
    }

    /// Validate all non-condensable gas fields fail-closed.
    pub fn validate(&self) -> Result<(), CondensationErrorV1> {
        if !self.air_mass_fraction.is_finite() || !(0.0..=1.0).contains(&self.air_mass_fraction) {
            return Err(CondensationErrorV1::InvalidNonCondensableFraction {
                air_mass_fraction: self.air_mass_fraction,
            });
        }
        if !self.bulk_temperature_k.is_finite() || self.bulk_temperature_k <= 0.0 {
            return Err(CondensationErrorV1::InvalidTemperature {
                temperature_k: self.bulk_temperature_k,
            });
        }
        if !self.gas_convective_htc_w_m2_k.is_finite()
            || self.gas_convective_htc_w_m2_k <= 0.0
            || !self.gas_cp_j_kg_k.is_finite()
            || self.gas_cp_j_kg_k <= 0.0
            || !self.lewis_number.is_finite()
            || self.lewis_number <= 0.0
        {
            return Err(CondensationErrorV1::InvalidParameter(
                "gas convective HTC, cp_g, and Lewis number must be finite positive",
            ));
        }
        Ok(())
    }
}

/// Evaluation of Nusselt / Colburn-Hougen filmwise condensation at one wall temperature `T_w`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CondensationEvaluationV1 {
    /// Evaluated wall temperature `T_w` [K].
    pub wall_temperature: Temperature,
    /// Bulk steam-air dew-point temperature `T_dew` [K].
    pub dew_temperature: Temperature,
    /// Liquid-vapor interface temperature `T_i` [K] (`T_i == T_dew` for pure steam;
    /// `T_w < T_i < T_dew` when non-condensable air is present).
    pub interface_temperature: Temperature,
    /// Active film condensation regime.
    pub regime: CondensationRegimeV1,
    /// Total fluid-to-wall condensation heat flux `q'' >= 0` [W/m^2].
    pub heat_flux_w_m2: f64,
    /// Effective overall condensation HTC `h_eff = q'' / (T_dew - T_w)` [W/(m^2 K)].
    pub effective_htc_w_m2_k: f64,
    /// Pure condensate film HTC `h_film = q'' / (T_i - T_w)` [W/(m^2 K)].
    pub film_htc_w_m2_k: f64,
    /// Analytical (or interface-consistent) derivative `dq''/dT_w <= 0` [W/(m^2 K)].
    pub dq_dtw_w_m2_k: f64,
    /// Condensate film Reynolds number `Re_delta = 4 * Gamma / mu_l`.
    pub film_reynolds: f64,
    /// Condensate mass flux per unit wall area `m_dot_cond'' >= 0` [kg/(m^2 s)].
    pub condensation_mass_flux_kg_m2_s: f64,
    /// Vapor specific enthalpy on `janaf-298.15K-formation` datum [J/kg].
    pub vapor_enthalpy_janaf_j_kg: f64,
    /// Condensate liquid specific enthalpy on `janaf-298.15K-formation` datum [J/kg].
    pub liquid_enthalpy_janaf_j_kg: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
enum CondensationFluidSourceV1 {
    IapwsSteam,
    Fixed(CondensationFluidPropertiesV1),
}

/// Prepared Nusselt filmwise condensation + Colburn-Hougen non-condensable gas model.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CondensationModelV1 {
    pressure_pa: f64,
    vapor_partial_pressure_pa: f64,
    dew_temperature_k: f64,
    geometry: CondensationGeometryV1,
    gas: NonCondensableGasV1,
    source: CondensationFluidSourceV1,
    allow_wavy_turbulent: bool,
}

impl CondensationModelV1 {
    /// Prepare a steam/steam-air filmwise condensation model backed by IAPWS-IF97 +
    /// R12-08 + R15-11 at caller-specified total pressure `pressure`.
    pub fn steam_air(
        pressure: Pressure,
        geometry: CondensationGeometryV1,
        gas: NonCondensableGasV1,
    ) -> Result<Self, CondensationErrorV1> {
        let p = pressure.value();
        if !p.is_finite() || p <= 0.0 {
            return Err(CondensationErrorV1::InvalidPressure { pressure_pa: p });
        }
        geometry.validate()?;
        gas.validate()?;

        let w_air = gas.air_mass_fraction;
        let w_v = 1.0 - w_air;
        let (p_v_bulk, t_dew) = if w_v <= 0.0 {
            // Pure dry air: no condensable steam; set dew temperature below any liquid state.
            (0.0, 0.0)
        } else {
            let n_v = w_v / WATER_MOLAR_MASS_KG_MOL_V1;
            let n_air = w_air / DRY_AIR_MOLAR_MASS_KG_MOL_V1;
            let y_v = n_v / (n_v + n_air);
            let pv = y_v * p;
            let td = saturation_temperature(Pressure::new(pv))
                .map_err(|_| CondensationErrorV1::InvalidPressure { pressure_pa: pv })?
                .value();
            (pv, td)
        };

        Ok(Self {
            pressure_pa: p,
            vapor_partial_pressure_pa: p_v_bulk,
            dew_temperature_k: t_dew,
            geometry,
            gas,
            source: CondensationFluidSourceV1::IapwsSteam,
            allow_wavy_turbulent: true,
        })
    }

    /// Prepare a pure-steam condensation model (`w_air = 0`) at caller-specified `pressure`.
    pub fn pure_steam(
        pressure: Pressure,
        geometry: CondensationGeometryV1,
    ) -> Result<Self, CondensationErrorV1> {
        let p = pressure.value();
        if !p.is_finite() || p <= 0.0 {
            return Err(CondensationErrorV1::InvalidPressure { pressure_pa: p });
        }
        let t_sat = saturation_temperature(pressure)
            .map_err(|_| CondensationErrorV1::InvalidPressure { pressure_pa: p })?
            .value();
        Self::steam_air(pressure, geometry, NonCondensableGasV1::pure_steam(t_sat))
    }

    /// Prepare a constant-property condensation model (ideal for exact analytical
    /// Nusselt verification or non-water condensables).
    pub fn fixed_properties(
        pressure: Pressure,
        dew_temperature: Temperature,
        geometry: CondensationGeometryV1,
        properties: CondensationFluidPropertiesV1,
        allow_wavy_turbulent: bool,
    ) -> Result<Self, CondensationErrorV1> {
        let p = pressure.value();
        let td = dew_temperature.value();
        if !p.is_finite() || p <= 0.0 {
            return Err(CondensationErrorV1::InvalidPressure { pressure_pa: p });
        }
        if !td.is_finite() || td <= 0.0 {
            return Err(CondensationErrorV1::InvalidTemperature { temperature_k: td });
        }
        geometry.validate()?;
        properties.validate()?;
        Ok(Self {
            pressure_pa: p,
            vapor_partial_pressure_pa: p,
            dew_temperature_k: td,
            geometry,
            gas: NonCondensableGasV1::pure_steam(td),
            source: CondensationFluidSourceV1::Fixed(properties),
            allow_wavy_turbulent,
        })
    }

    /// Restrict plate correlations to pure laminar Nusselt (disabling wavy/turbulent
    /// film transitions) or re-enable wavy-laminar and turbulent film regimes.
    #[must_use]
    pub const fn with_wavy_turbulent(mut self, allow: bool) -> Self {
        self.allow_wavy_turbulent = allow;
        self
    }

    /// Bulk dew-point temperature `T_dew(P, w_air)` [K] (`0.0` for pure dry air).
    #[must_use]
    pub const fn dew_temperature(&self) -> Temperature {
        Temperature::new(self.dew_temperature_k)
    }

    /// Bulk vapor partial pressure `p_{v,bulk}` [Pa].
    #[must_use]
    pub const fn vapor_partial_pressure(&self) -> Pressure {
        Pressure::new(self.vapor_partial_pressure_pa)
    }

    /// Non-condensable gas specification.
    #[must_use]
    pub const fn non_condensable_gas(&self) -> &NonCondensableGasV1 {
        &self.gas
    }

    /// Evaluate filmwise condensation heat flux, HTC, condensate mass flux, and
    /// derivative `dq''/dT_w` at `wall_temperature`.
    pub fn evaluate(
        &self,
        wall_temperature: Temperature,
    ) -> Result<CondensationEvaluationV1, CondensationErrorV1> {
        let tw = wall_temperature.value();
        if !tw.is_finite() || tw <= 0.0 {
            return Err(CondensationErrorV1::InvalidTemperature { temperature_k: tw });
        }

        // Single-source handoff guard: zero condensation whenever T_w >= T_dew or w_air == 1.
        if self.gas.air_mass_fraction >= 1.0 || tw >= self.dew_temperature_k {
            let (h_v_j, h_l_j) = self.default_janaf_enthalpies();
            return Ok(CondensationEvaluationV1 {
                wall_temperature,
                dew_temperature: Temperature::new(self.dew_temperature_k),
                interface_temperature: Temperature::new(self.dew_temperature_k),
                regime: CondensationRegimeV1::NoCondensation,
                heat_flux_w_m2: 0.0,
                effective_htc_w_m2_k: 0.0,
                film_htc_w_m2_k: 0.0,
                dq_dtw_w_m2_k: 0.0,
                film_reynolds: 0.0,
                condensation_mass_flux_kg_m2_s: 0.0,
                vapor_enthalpy_janaf_j_kg: h_v_j,
                liquid_enthalpy_janaf_j_kg: h_l_j,
            });
        }

        // Pure steam (no non-condensable air): interface sits at T_i = T_dew.
        if self.gas.air_mass_fraction == 0.0 {
            return self.evaluate_pure_film(tw, self.dew_temperature_k);
        }

        // Non-condensable air present: solve Colburn-Hougen interface balance on (T_w, T_dew).
        self.evaluate_colburn_hougen(tw)
    }

    fn default_janaf_enthalpies(&self) -> (f64, f64) {
        match self.source {
            CondensationFluidSourceV1::Fixed(p) => {
                (p.vapor_enthalpy_janaf_j_kg, p.liquid_enthalpy_janaf_j_kg)
            }
            CondensationFluidSourceV1::IapwsSteam => {
                if self.dew_temperature_k > 273.16
                    && let Ok(pair) =
                        If97SaturationPair::at_temperature(Temperature::new(self.dew_temperature_k))
                    && let (Ok(hv), Ok(hl)) = (
                        water_state_enthalpy_janaf(pair.vapor()),
                        water_state_enthalpy_janaf(pair.liquid()),
                    )
                {
                    return (hv.value(), hl.value());
                }
                (0.0, 0.0)
            }
        }
    }

    fn film_properties_at(
        &self,
        t_interface_k: f64,
        t_wall_k: f64,
    ) -> Result<CondensationFluidPropertiesV1, CondensationErrorV1> {
        match self.source {
            CondensationFluidSourceV1::Fixed(p) => Ok(p),
            CondensationFluidSourceV1::IapwsSteam => {
                let pair = If97SaturationPair::at_temperature(Temperature::new(t_interface_k))
                    .map_err(|_| {
                        CondensationErrorV1::Thermochem("IAPWS saturation pair at T_i failed")
                    })?;
                let t_film = (0.5 * (t_interface_k + t_wall_k)).max(273.20);
                // Evaluate subcooled condensate liquid at P_total (or slightly above P_sat(T_film))
                // so Region 1 stays strictly in the liquid phase.
                let p_sat_film = saturation_pressure(Temperature::new(t_film))
                    .map_err(|_| CondensationErrorV1::Thermochem("IAPWS P_sat(T_film) failed"))?
                    .value();
                let p_liq = self.pressure_pa.max(p_sat_film + 100.0);
                let liq_state = If97State::at_pressure_temperature(
                    Pressure::new(p_liq),
                    Temperature::new(t_film),
                    If97Region::Region1,
                )
                .unwrap_or(*pair.liquid());
                let liq_tr = If97TransportState::from_state(If97FluidState::Gibbs(liq_state))
                    .map_err(|_| {
                        CondensationErrorV1::Thermochem("IAPWS condensate liquid transport failed")
                    })?;
                let h_v_janaf = water_state_enthalpy_janaf(pair.vapor())
                    .map_err(|_| CondensationErrorV1::Thermochem("JANAF vapor enthalpy failed"))?
                    .value();
                let h_l_janaf = water_state_enthalpy_janaf(&liq_state)
                    .map_err(|_| CondensationErrorV1::Thermochem("JANAF liquid enthalpy failed"))?
                    .value();

                let props = CondensationFluidPropertiesV1 {
                    liquid_density_kg_m3: liq_state.density().value(),
                    vapor_density_kg_m3: pair.vapor().density().value(),
                    liquid_viscosity_pa_s: liq_tr.dynamic_viscosity().value(),
                    liquid_conductivity_w_m_k: liq_tr.thermal_conductivity().value(),
                    liquid_cp_j_kg_k: liq_state.cp().value(),
                    latent_heat_j_kg: pair.latent_enthalpy().value(),
                    vapor_enthalpy_janaf_j_kg: h_v_janaf,
                    liquid_enthalpy_janaf_j_kg: h_l_janaf,
                };
                props.validate()?;
                Ok(props)
            }
        }
    }

    fn evaluate_pure_film(
        &self,
        tw: f64,
        t_interface_k: f64,
    ) -> Result<CondensationEvaluationV1, CondensationErrorV1> {
        let props = self.film_properties_at(t_interface_k, tw)?;
        let (regime, h_film, q_film, dq_dtw_film, re_delta, m_dot) = nusselt_film_eval(
            self.geometry,
            props,
            t_interface_k - tw,
            self.allow_wavy_turbulent,
        )?;
        let dt_overall = self.dew_temperature_k - tw;
        let h_eff = if dt_overall > 0.0 {
            q_film / dt_overall
        } else {
            h_film
        };
        Ok(CondensationEvaluationV1 {
            wall_temperature: Temperature::new(tw),
            dew_temperature: Temperature::new(self.dew_temperature_k),
            interface_temperature: Temperature::new(t_interface_k),
            regime,
            heat_flux_w_m2: q_film,
            effective_htc_w_m2_k: h_eff,
            film_htc_w_m2_k: h_film,
            dq_dtw_w_m2_k: dq_dtw_film,
            film_reynolds: re_delta,
            condensation_mass_flux_kg_m2_s: m_dot,
            vapor_enthalpy_janaf_j_kg: props.vapor_enthalpy_janaf_j_kg,
            liquid_enthalpy_janaf_j_kg: props.liquid_enthalpy_janaf_j_kg,
        })
    }

    #[allow(clippy::too_many_lines)]
    fn evaluate_colburn_hougen(
        &self,
        tw: f64,
    ) -> Result<CondensationEvaluationV1, CondensationErrorV1> {
        let t_dew = self.dew_temperature_k;
        let p_tot = self.pressure_pa;
        let p_v_bulk = self.vapor_partial_pressure_pa;
        let p_nc_bulk = p_tot - p_v_bulk;
        if !p_nc_bulk.is_finite() || p_nc_bulk <= 0.0 {
            return self.evaluate_pure_film(tw, t_dew);
        }

        // Chilton-Colburn gas-side mass-transfer conductance g_m [kg/(m^2 s)]:
        // g_m = h_g / (c_{p,g} * Le^(2/3))
        let g_m = self.gas.gas_convective_htc_w_m2_k
            / (self.gas.gas_cp_j_kg_k * pow(self.gas.lewis_number, 2.0 / 3.0));
        let h_g = self.gas.gas_convective_htc_w_m2_k;
        let t_bulk = self.gas.bulk_temperature_k.max(t_dew);

        let balance = |t_i: f64| -> Result<
            (
                f64,
                f64,
                f64,
                f64,
                f64,
                CondensationRegimeV1,
                CondensationFluidPropertiesV1,
            ),
            CondensationErrorV1,
        > {
            let props = self.film_properties_at(t_i, tw)?;
            let (regime, h_film, q_film, _, re_delta, _) =
                nusselt_film_eval(self.geometry, props, t_i - tw, self.allow_wavy_turbulent)?;
            let p_v_i = saturation_pressure(Temperature::new(t_i))
                .map_err(|_| CondensationErrorV1::Thermochem("IAPWS P_sat(T_i) failed"))?
                .value();
            let p_nc_i = (p_tot - p_v_i).max(1e-6);
            // Stefan logarithmic driving force: ln((P - p_{v,i}) / (P - p_{v,bulk})) >= 0 for T_i <= T_dew.
            let m_dot_diff = (g_m * ln(p_nc_i / p_nc_bulk)).max(0.0);
            let q_sens = h_g * (t_bulk - t_i);
            let q_gas = q_sens + m_dot_diff * props.latent_heat_j_kg;
            let residual = q_film - q_gas;
            Ok((
                residual, q_film, h_film, m_dot_diff, re_delta, regime, props,
            ))
        };

        // Bracket on [tw + eps, t_dew].
        let span = t_dew - tw;
        let mut lo = tw + (1e-6 * span).min(1e-5);
        let mut hi = t_dew;

        let (r_lo, _, _, _, _, _, _) = balance(lo)?;
        let (r_hi, _, _, _, _, _, _) = balance(hi)?;

        // If superheated bulk sensible heating exceeds film conduction even at T_i = T_dew,
        // the surface stays dry (zero condensation).
        if r_hi <= 0.0 {
            let (h_v_j, h_l_j) = self.default_janaf_enthalpies();
            return Ok(CondensationEvaluationV1 {
                wall_temperature: Temperature::new(tw),
                dew_temperature: Temperature::new(t_dew),
                interface_temperature: Temperature::new(t_dew),
                regime: CondensationRegimeV1::NoCondensation,
                heat_flux_w_m2: 0.0,
                effective_htc_w_m2_k: 0.0,
                film_htc_w_m2_k: 0.0,
                dq_dtw_w_m2_k: 0.0,
                film_reynolds: 0.0,
                condensation_mass_flux_kg_m2_s: 0.0,
                vapor_enthalpy_janaf_j_kg: h_v_j,
                liquid_enthalpy_janaf_j_kg: h_l_j,
            });
        }
        if r_lo >= 0.0 {
            // Very tiny air fraction / extreme film resistance: interface essentially at lo.
            hi = lo;
        }

        let mut t_i = 0.5 * (lo + hi);
        let mut last_res = f64::INFINITY;
        for _ in 0..80 {
            t_i = 0.5 * (lo + hi);
            let (res, q_film, _, _, _, _, _) = balance(t_i)?;
            last_res = res;
            if res.abs() <= 1e-10 * q_film.abs().max(1.0) || (hi - lo) <= 1e-12 {
                break;
            }
            if res < 0.0 {
                lo = t_i;
            } else {
                hi = t_i;
            }
        }
        let (res, q_film, h_film, m_dot_diff, re_delta, regime, props) = balance(t_i)?;
        if res.abs() > 1e-3 * q_film.abs().max(1.0) && (hi - lo) > 1e-9 {
            return Err(CondensationErrorV1::ConvergenceFailed {
                iterations: 80,
                residual_w_m2: last_res,
            });
        }

        // Compute tangent dq''/dT_w via central difference on the resolved Colburn-Hougen curve.
        let dt_overall = t_dew - tw;
        let h_eff = q_film / dt_overall;
        let eps = (1e-4 * dt_overall).clamp(1e-6, 1e-3);
        let dq_dtw = if tw - eps > 273.2 && tw + eps < t_dew {
            let qp = self.solve_colburn_hougen_flux(tw + eps)?;
            let qm = self.solve_colburn_hougen_flux(tw - eps)?;
            (qp - qm) / (2.0 * eps)
        } else {
            -h_eff
        };

        Ok(CondensationEvaluationV1 {
            wall_temperature: Temperature::new(tw),
            dew_temperature: Temperature::new(t_dew),
            interface_temperature: Temperature::new(t_i),
            regime,
            heat_flux_w_m2: q_film,
            effective_htc_w_m2_k: h_eff,
            film_htc_w_m2_k: h_film,
            dq_dtw_w_m2_k: dq_dtw,
            film_reynolds: re_delta,
            condensation_mass_flux_kg_m2_s: m_dot_diff,
            vapor_enthalpy_janaf_j_kg: props.vapor_enthalpy_janaf_j_kg,
            liquid_enthalpy_janaf_j_kg: props.liquid_enthalpy_janaf_j_kg,
        })
    }

    fn solve_colburn_hougen_flux(&self, tw: f64) -> Result<f64, CondensationErrorV1> {
        let t_dew = self.dew_temperature_k;
        let p_tot = self.pressure_pa;
        let p_nc_bulk = p_tot - self.vapor_partial_pressure_pa;
        let g_m = self.gas.gas_convective_htc_w_m2_k
            / (self.gas.gas_cp_j_kg_k * pow(self.gas.lewis_number, 2.0 / 3.0));
        let h_g = self.gas.gas_convective_htc_w_m2_k;
        let t_bulk = self.gas.bulk_temperature_k.max(t_dew);

        let mut lo = tw + (1e-6 * (t_dew - tw)).min(1e-5);
        let mut hi = t_dew;

        // Match `evaluate_colburn_hougen`'s dry-surface check at `T_i = T_dew`.
        let props_hi = self.film_properties_at(hi, tw)?;
        let (_, _, q_film_hi, _, _, _) =
            nusselt_film_eval(self.geometry, props_hi, hi - tw, self.allow_wavy_turbulent)?;
        if q_film_hi <= h_g * (t_bulk - hi) {
            return Ok(0.0);
        }

        for _ in 0..72 {
            let t_i = 0.5 * (lo + hi);
            let props = self.film_properties_at(t_i, tw)?;
            let (_, _, q_film, _, _, _) =
                nusselt_film_eval(self.geometry, props, t_i - tw, self.allow_wavy_turbulent)?;
            let p_v_i = saturation_pressure(Temperature::new(t_i))
                .map_err(|_| CondensationErrorV1::Thermochem("IAPWS P_sat(T_i) failed"))?
                .value();
            let p_nc_i = (p_tot - p_v_i).max(1e-6);
            let m_dot = (g_m * ln(p_nc_i / p_nc_bulk)).max(0.0);
            let q_gas = h_g * (t_bulk - t_i) + m_dot * props.latent_heat_j_kg;
            if (q_film - q_gas).abs() <= 1e-11 * q_film.abs().max(1.0) || (hi - lo) <= 1e-12 {
                return Ok(q_film);
            }
            if q_film < q_gas {
                lo = t_i;
            } else {
                hi = t_i;
            }
        }
        let t_i = 0.5 * (lo + hi);
        let props = self.film_properties_at(t_i, tw)?;
        let (_, _, q_film, _, _, _) =
            nusselt_film_eval(self.geometry, props, t_i - tw, self.allow_wavy_turbulent)?;
        Ok(q_film)
    }
}

/// Core Nusselt (1916) / Kutateladze / Labuntsov filmwise condensation evaluation.
///
/// Returns `(regime, h_film, q_film, dq_dtw, re_delta, m_dot_cond)`.
#[allow(clippy::too_many_lines)]
fn nusselt_film_eval(
    geometry: CondensationGeometryV1,
    props: CondensationFluidPropertiesV1,
    dt: f64,
    allow_wavy_turbulent: bool,
) -> Result<(CondensationRegimeV1, f64, f64, f64, f64, f64), CondensationErrorV1> {
    if !dt.is_finite() || dt <= 0.0 {
        return Err(CondensationErrorV1::InvalidTemperature { temperature_k: dt });
    }
    let rho_l = props.liquid_density_kg_m3;
    let rho_v = props.vapor_density_kg_m3;
    let mu_l = props.liquid_viscosity_pa_s;
    let k_l = props.liquid_conductivity_w_m_k;
    let cp_l = props.liquid_cp_j_kg_k;
    let h_fg = props.latent_heat_j_kg;
    // Rohsenow modified latent heat (Incropera Eq. 10.24).
    let h_fg_prime = h_fg + 0.68 * cp_l * dt;

    let (coeff, g_eff, length_m, is_plate) = match geometry {
        CondensationGeometryV1::VerticalPlate { height_m } => {
            (0.943, STANDARD_GRAVITY_M_S2, height_m, true)
        }
        CondensationGeometryV1::InclinedPlate {
            length_m,
            tilt_from_horizontal_rad,
        } => (
            0.943,
            STANDARD_GRAVITY_M_S2 * sin(tilt_from_horizontal_rad),
            length_m,
            true,
        ),
        CondensationGeometryV1::HorizontalCylinder { diameter_m } => {
            (0.729, STANDARD_GRAVITY_M_S2, diameter_m, false)
        }
    };

    // Laminar Nusselt correlation (Incropera Eqs. 10.30 / 10.40):
    // h_lam = C * [ g_eff * rho_l * (rho_l - rho_v) * k_l^3 * h'_fg / (mu_l * dT * L) ]^(1/4)
    let group =
        (g_eff * rho_l * (rho_l - rho_v) * k_l * k_l * k_l * h_fg_prime) / (mu_l * dt * length_m);
    let h_lam = coeff * pow(group, 0.25);
    let q_lam = h_lam * dt;
    let m_dot_lam = q_lam / h_fg_prime;
    let re_lam = (4.0 * m_dot_lam * length_m) / mu_l;

    // Exact derivative of q_lam(dT) = C * [... / (mu_l * L)]^(1/4) * dT^(3/4) * (h_fg + 0.68*cp_l*dT)^(1/4)
    // w.r.t. T_w (where d(dT)/dT_w = -1):
    let dq_lam_dtw = -q_lam * (0.75 / dt + (0.17 * cp_l) / h_fg_prime);

    if !is_plate || !allow_wavy_turbulent || re_lam <= NUSSELT_WAVY_RE_THRESHOLD_V1 {
        return Ok((
            CondensationRegimeV1::LaminarFilm,
            h_lam,
            q_lam,
            dq_lam_dtw,
            re_lam,
            m_dot_lam,
        ));
    }

    // Dimensionless film group P = k_l * L * dT / (mu_l * h'_fg * (nu_l^2 / g_eff)^(1/3)).
    let nu_l = mu_l / rho_l;
    let len_scale = pow((nu_l * nu_l) / g_eff, 1.0 / 3.0);
    let p_dim = (k_l * length_m * dt) / (mu_l * h_fg_prime * len_scale);

    // Wavy-laminar (Kutateladze / Incropera Eqs. 10.37-10.38).
    let re_wavy = pow(3.70 * p_dim + 4.8, 0.82);
    if re_wavy <= NUSSELT_TURBULENT_RE_THRESHOLD_V1 {
        let denom = 1.08 * pow(re_wavy, 1.22) - 5.2;
        if !denom.is_finite() || denom <= 0.0 {
            return Err(CondensationErrorV1::Unrepresentable(
                "wavy-laminar Nusselt denominator non-positive",
            ));
        }
        let h_wavy = (k_l / len_scale) * (re_wavy / denom);
        let q_wavy = h_wavy * dt;
        let m_dot = q_wavy / h_fg_prime;
        // Analytical derivative via d(p_dim)/dT_w = -p_dim * (1/dT - 0.68*cp_l/h'_fg):
        let dp_ddt = p_dim * (1.0 / dt - (0.68 * cp_l) / h_fg_prime);
        let dre_ddt = 0.82 * (3.70 * dp_ddt) * pow(3.70 * p_dim + 4.8, -0.18);
        let ddenom_dre = 1.08 * 1.22 * pow(re_wavy, 0.22);
        let dh_dre = (k_l / len_scale) * (denom - re_wavy * ddenom_dre) / (denom * denom);
        let dq_dtw = -(h_wavy + dt * dh_dre * dre_ddt);
        return Ok((
            CondensationRegimeV1::WavyLaminarFilm,
            h_wavy,
            q_wavy,
            dq_dtw,
            re_wavy,
            m_dot,
        ));
    }

    // Turbulent film (Labuntsov / Incropera Eq. 10.39).
    let pr_l = (mu_l * cp_l) / k_l;
    let sqrt_pr = pow(pr_l, 0.5);
    let base = 0.0690 * sqrt_pr * p_dim - 151.0 * sqrt_pr + 253.0;
    if !base.is_finite() || base <= 0.0 {
        return Err(CondensationErrorV1::Unrepresentable(
            "turbulent Nusselt film Reynolds base non-positive",
        ));
    }
    let re_turb = pow(base, 4.0 / 3.0);
    let denom = 8_750.0 + (58.0 / sqrt_pr) * (pow(re_turb, 0.75) - 253.0);
    if !denom.is_finite() || denom <= 0.0 {
        return Err(CondensationErrorV1::Unrepresentable(
            "turbulent Nusselt denominator non-positive",
        ));
    }
    let h_turb = (k_l / len_scale) * (re_turb / denom);
    let q_turb = h_turb * dt;
    let m_dot = q_turb / h_fg_prime;

    let dp_ddt = p_dim * (1.0 / dt - (0.68 * cp_l) / h_fg_prime);
    let dre_ddt = (4.0 / 3.0) * (0.0690 * sqrt_pr * dp_ddt) * pow(base, 1.0 / 3.0);
    let ddenom_dre = (58.0 / sqrt_pr) * 0.75 * pow(re_turb, -0.25);
    let dh_dre = (k_l / len_scale) * (denom - re_turb * ddenom_dre) / (denom * denom);
    let dq_dtw = -(h_turb + dt * dh_dre * dre_ddt);

    Ok((
        CondensationRegimeV1::TurbulentFilm,
        h_turb,
        q_turb,
        dq_dtw,
        re_turb,
        m_dot,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pure_laminar_vertical_plate_matches_analytical_nusselt_and_fd_tangent() {
        let props = CondensationFluidPropertiesV1 {
            liquid_density_kg_m3: 960.0,
            vapor_density_kg_m3: 0.60,
            liquid_viscosity_pa_s: 2.8e-4,
            liquid_conductivity_w_m_k: 0.68,
            liquid_cp_j_kg_k: 4_216.0,
            latent_heat_j_kg: 2.257e6,
            vapor_enthalpy_janaf_j_kg: -13.29e6,
            liquid_enthalpy_janaf_j_kg: -15.55e6,
        };
        let model = CondensationModelV1::fixed_properties(
            Pressure::new(101_325.0),
            Temperature::new(373.15),
            CondensationGeometryV1::VerticalPlate { height_m: 0.05 },
            props,
            false,
        )
        .unwrap();
        let tw = 368.15; // dT = 5 K -> laminar Re_delta < 30
        let ev = model.evaluate(Temperature::new(tw)).unwrap();
        assert_eq!(ev.regime, CondensationRegimeV1::LaminarFilm);

        let dt = 5.0;
        let h_fg_p = props.latent_heat_j_kg + 0.68 * props.liquid_cp_j_kg_k * dt;
        let expected_h = 0.943
            * pow(
                (STANDARD_GRAVITY_M_S2
                    * props.liquid_density_kg_m3
                    * (props.liquid_density_kg_m3 - props.vapor_density_kg_m3)
                    * props.liquid_conductivity_w_m_k
                    * props.liquid_conductivity_w_m_k
                    * props.liquid_conductivity_w_m_k
                    * h_fg_p)
                    / (props.liquid_viscosity_pa_s * dt * 0.05),
                0.25,
            );
        assert!((ev.film_htc_w_m2_k - expected_h).abs() / expected_h < 1e-12);

        let eps = 1e-4;
        let qp = model
            .evaluate(Temperature::new(tw + eps))
            .unwrap()
            .heat_flux_w_m2;
        let qm = model
            .evaluate(Temperature::new(tw - eps))
            .unwrap()
            .heat_flux_w_m2;
        let fd = (qp - qm) / (2.0 * eps);
        assert!((ev.dq_dtw_w_m2_k - fd).abs() / fd.abs() < 1e-9);
    }

    #[test]
    fn colburn_hougen_degrades_htc_monotonically_with_air_mass_fraction() {
        let p = Pressure::new(101_325.0);
        let geom = CondensationGeometryV1::VerticalPlate { height_m: 0.10 };
        let tw = Temperature::new(330.0);

        let mut prev_h = f64::INFINITY;
        for w_air in [0.0, 0.05, 0.25, 0.50, 0.70] {
            let gas =
                NonCondensableGasV1::steam_air(w_air, Temperature::new(373.15), 25.0).unwrap();
            let model = CondensationModelV1::steam_air(p, geom, gas).unwrap();
            let ev = model.evaluate(tw).unwrap();
            assert!(
                ev.effective_htc_w_m2_k < prev_h,
                "air mass fraction {w_air} must strictly reduce effective HTC"
            );
            prev_h = ev.effective_htc_w_m2_k;
        }
    }
}
