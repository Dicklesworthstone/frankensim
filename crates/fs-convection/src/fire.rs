//! Nominal, source-declared compartment-fire correlations in coherent SI.
//!
//! Heskestad geometry uses total heat release, while plume buoyancy uses its
//! explicitly declared convective part and actual ambient rho, Cp, T and g.
//! The source's linear flaming-region entrainment is distinguished from its
//! above-flame plume. Centerline similarity expressions refuse flame-region
//! queries; neither temperatures nor velocities are silently clipped.
//!
//! These inexpensive estimates do not evolve a compartment, burn species,
//! resolve ventilation, or certify fire safety. Every result retains a static
//! model card with qualitative applicability and **unquantified** discrepancy.
//! No finite uncertainty band or numerical calibration box is invented. This
//! distinct nominal wrapper does not manufacture a generic certified Evidence.
//! All evaluators are allocation-free bounded scalar algebra using deterministic
//! elementary functions; there is no unbounded work requiring cancellation.

use crate::{HeatTransferCoefficient, SourceProvenance};
use fs_evidence::ValidityDomain;
use fs_math::det;
use fs_qty::{
    Acceleration, Area, Density, Length, MassFlowRate, Power, Qty, Temperature, Velocity,
};

/// Mass-specific heat capacity in J/(kg K).
pub type SpecificHeatCapacityV1 = Qty<2, 0, -2, -1, 0, 0>;

/// Quantity-specific prediction ratios and uncertainty-aware validation summaries.
pub mod validation;
/// Mass-specific energy in J/kg.
pub type SpecificEnergyV1 = Qty<2, 0, -2, 0, 0, 0>;

/// Stable implemented empirical relation, separate from convection's Nu cards.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FireCorrelationIdV1 {
    /// Mean visible flame height and virtual origin from total heat release.
    HeskestadFlame,
    /// Unobstructed plume and linear flaming-region entrainment.
    HeskestadEntrainment,
    /// Above-flame centerline temperature and velocity similarity expressions.
    HeskestadCenterline,
    /// McCaffrey--Quintiere--Harkleroad naturally ventilated hot layer.
    MqhHotLayer,
    /// Thomas's minimum heat-release indicator for flashover.
    ThomasFlashover,
    /// Babrauskas's ventilation-factor flashover indicator.
    BabrauskasFlashover,
    /// Kawagoe-type post-flashover ventilation-limited burning capacity.
    KawagoeVentilation,
}
impl FireCorrelationIdV1 {
    /// Complete, stable source-card order; no preference between models implied.
    pub const ALL: [Self; 7] = [
        Self::HeskestadFlame,
        Self::HeskestadEntrainment,
        Self::HeskestadCenterline,
        Self::MqhHotLayer,
        Self::ThomasFlashover,
        Self::BabrauskasFlashover,
        Self::KawagoeVentilation,
    ];
    /// Original source, applicability and unquantified model-form uncertainty.
    pub const fn card(self) -> &'static FireModelCardV1 {
        match self {
            Self::HeskestadFlame => &HESKESTAD_FLAME,
            Self::HeskestadEntrainment => &HESKESTAD_ENTRAINMENT,
            Self::HeskestadCenterline => &HESKESTAD_CENTERLINE,
            Self::MqhHotLayer => &MQH,
            Self::ThomasFlashover => &THOMAS,
            Self::BabrauskasFlashover => &BABRAUSKAS,
            Self::KawagoeVentilation => &KAWAGOE,
        }
    }
}

/// Model discrepancy is explicitly unknown, not zero or a guessed percentage.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FireDiscrepancyV1 {
    /// No retained calibration artifact quantifies error for this application.
    Unquantified,
}
/// Empirical arithmetic is a nominal prediction, not a validation promotion.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FirePredictionStatusV1 {
    /// Source formula evaluated under its declared assumptions.
    Nominal,
}
/// Source card for a nominal fire relation with no fabricated numerical band.
///
/// Unlike the generic ModelCard's mandatory numeric discrepancy, this card can
/// express unquantified model error without publishing a nonfinite solver value.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FireModelCardV1 {
    /// Stable implemented relation.
    pub correlation: FireCorrelationIdV1,
    /// Stable model name.
    pub name: &'static str,
    /// Evaluator/formula version.
    pub version: u32,
    /// Bibliographic authority, including equation or section.
    pub source: SourceProvenance,
    /// Primary document location.
    pub source_url: &'static str,
    /// Nominal evaluation status; no certified-result conversion is provided.
    pub status: FirePredictionStatusV1,
    /// Explicit unknown discrepancy.
    pub discrepancy: FireDiscrepancyV1,
    /// Source/application conditions; positivity alone does not establish them.
    pub domain_assumptions: &'static [&'static str],
    /// No quantitative calibration box is retained for these nominal models.
    pub quantitative_validity_domain: Option<&'static ValidityDomain>,
    /// Physical effects and interpretation beyond the implemented relation.
    pub known_failures: &'static [&'static str],
}

/// Typed nominal result retaining its exact static source card.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FireEstimateV1<T> {
    value: T,
    correlation: FireCorrelationIdV1,
}
impl<T: Copy> FireEstimateV1<T> {
    /// Coherent-SI nominal value; a detached value has no certification authority.
    pub const fn value(&self) -> T {
        self.value
    }
}
impl<T> FireEstimateV1<T> {
    /// Source, nominal status, applicability and unquantified discrepancy.
    pub const fn card(&self) -> &'static FireModelCardV1 {
        self.correlation.card()
    }
}

/// Refusal preserves invalid-input bits or names the physical regime boundary.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum FireCorrelationErrorV1 {
    /// A dimensional input or heat-release partition is inadmissible.
    InvalidInput {
        /// Input name, including units where useful.
        field: &'static str,
        /// Exact offered IEEE-754 bits.
        bits: u64,
    },
    /// The selected empirical branch does not represent this physical regime.
    Regime {
        /// Selected source relation.
        correlation: FireCorrelationIdV1,
        /// Concrete reason no nominal prediction can be returned.
        reason: &'static str,
    },
    /// Required finite physical arithmetic overflowed or underflowed.
    Arithmetic(&'static str),
}
impl core::fmt::Display for FireCorrelationErrorV1 {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "nominal fire correlation: {self:?}")
    }
}
impl std::error::Error for FireCorrelationErrorV1 {}

/// Actual ambient properties for the dimensional Heskestad equations.
///
/// Cp is evaluated at the ambient state supplied by the caller. This local
/// constant-property plume approximation does not replace the owning gas EOS.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FireAmbientV1 {
    temperature: Temperature,
    density: Density,
    specific_heat: SpecificHeatCapacityV1,
    gravity: Acceleration,
}
impl FireAmbientV1 {
    /// Admit positive finite absolute T, density, Cp and gravity, all in SI.
    pub fn new(
        temperature: Temperature,
        density: Density,
        specific_heat: SpecificHeatCapacityV1,
        gravity: Acceleration,
    ) -> Result<Self, FireCorrelationErrorV1> {
        for (name, value) in [
            ("ambient temperature K", temperature.value()),
            ("ambient density kg/m3", density.value()),
            ("ambient Cp J/(kg K)", specific_heat.value()),
            ("gravity m/s2", gravity.value()),
        ] {
            positive_input(value, name)?;
        }
        Ok(Self {
            temperature,
            density,
            specific_heat,
            gravity,
        })
    }
    /// Actual absolute ambient temperature.
    pub const fn temperature(&self) -> Temperature {
        self.temperature
    }
    /// Actual ambient mass density.
    pub const fn density(&self) -> Density {
        self.density
    }
    /// Caller-evaluated ambient mass-specific Cp.
    pub const fn specific_heat(&self) -> SpecificHeatCapacityV1 {
        self.specific_heat
    }
    /// Declared gravitational acceleration magnitude.
    pub const fn gravity(&self) -> Acceleration {
        self.gravity
    }
}

/// Branch used by the source's piecewise entrainment prescription.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HeskestadRegionV1 {
    /// Exactly zero total and convective heat release; no plume exists.
    Dormant,
    /// Below the mean flame height: m(z)=m(L) z/L, CFAST Eq. (3.17).
    Flaming,
    /// At or above the mean flame height: CFAST Eq. (3.15).
    AboveFlame,
}
/// Entrainment across a horizontal plane above the fuel source.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HeskestadEntrainmentV1 {
    /// Nominal gas mass transported upwards [kg/s]; not a fuel production rate.
    pub mass_rate: MassFlowRate,
    /// Explicit executed source branch.
    pub region: HeskestadRegionV1,
}
/// Unclipped above-flame centerline similarity prediction.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HeskestadCenterlineV1 {
    /// Temperature difference above ambient [K], not an absolute temperature.
    pub temperature_rise: Temperature,
    /// Absolute centerline temperature [K].
    pub temperature: Temperature,
    /// Upward centerline velocity [m/s].
    pub velocity: Velocity,
}

/// One immutable source and ambient state for Heskestad's unobstructed plume.
///
/// The CFAST TN 1889v1 dimensional equations are evaluated directly in SI.
/// Q*=Q/(rho Cp T sqrt(gD) D^2); geometry uses Q, while entrainment and
/// centerline properties use Qc. No default radiative fraction is assumed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HeskestadFireV1 {
    total_heat_release: Power,
    convective_heat_release: Power,
    diameter: Length,
    ambient: FireAmbientV1,
    dimensionless_heat_release: f64,
    flame_height: Length,
    virtual_origin: Option<Length>,
    mass_scale: f64,
    mass_offset: f64,
    temperature_scale: f64,
    velocity_scale: f64,
}
impl HeskestadFireV1 {
    /// Admit an explicit 0<=Qc<=Q partition and positive source diameter.
    /// Exactly Q=Qc=0 is dormant. A positive source whose mean-flame fit is
    /// nonpositive refuses rather than inventing or clipping a flame height.
    pub fn new(
        total_heat_release: Power,
        convective_heat_release: Power,
        diameter: Length,
        ambient: FireAmbientV1,
    ) -> Result<Self, FireCorrelationErrorV1> {
        nonnegative_input(total_heat_release.value(), "total HRR W")?;
        nonnegative_input(convective_heat_release.value(), "convective HRR W")?;
        positive_input(diameter.value(), "source diameter m")?;
        if convective_heat_release.value() > total_heat_release.value() {
            return Err(invalid(
                "convective HRR exceeds total HRR W",
                convective_heat_release.value(),
            ));
        }
        let mut fire = Self {
            total_heat_release,
            convective_heat_release,
            diameter,
            ambient,
            dimensionless_heat_release: 0.0,
            flame_height: Length::new(0.0),
            virtual_origin: None,
            mass_scale: 0.0,
            mass_offset: 0.0,
            temperature_scale: 0.0,
            velocity_scale: 0.0,
        };
        if total_heat_release.value() == 0.0 {
            return Ok(fire);
        }
        let (rho, cp, t, g, d) = (
            ambient.density.value(),
            ambient.specific_heat.value(),
            ambient.temperature.value(),
            ambient.gravity.value(),
            diameter.value(),
        );
        let sqrt_g = det::sqrt(g);
        let heat_scale = positive_result(rho * cp * t * sqrt_g, "ambient buoyancy scale")?;
        let diameter_scale = positive_result(d * d * det::sqrt(d), "source diameter scale")?;
        let q_star = positive_result(
            (total_heat_release.value() / heat_scale) / diameter_scale,
            "dimensionless total heat release",
        )?;
        let power = pow_positive(q_star, 2.0 / 5.0);
        let flame_height = finite_result(d * (-1.02 + 3.7 * power), "mean flame height")?;
        if flame_height <= 0.0 {
            return Err(FireCorrelationErrorV1::Regime {
                correlation: FireCorrelationIdV1::HeskestadFlame,
                reason: "positive fire has no positive mean flame height under this source fit",
            });
        }
        let origin = finite_result(d * (-1.02 + 1.4 * power), "virtual origin")?;
        positive_result(flame_height - origin, "flame height above virtual origin")?;
        fire.dimensionless_heat_release = q_star;
        fire.flame_height = Length::new(flame_height);
        fire.virtual_origin = Some(Length::new(origin));
        if convective_heat_release.value() == 0.0 {
            return Ok(fire);
        }
        // S=Qc/(rho Cp T sqrt(g)) has units m^(5/2). Rewriting the source
        // equations through S avoids redundant powers and unit conversion.
        let s = positive_result(
            convective_heat_release.value() / heat_scale,
            "convective buoyancy scale",
        )?;
        let s_one_third = positive_result(pow_positive(s, 1.0 / 3.0), "convective cube root")?;
        fire.mass_scale = positive_result(0.196 * rho * sqrt_g * s_one_third, "entrainment scale")?;
        // Algebraic expansion of Eq3.15: the correction is 0.196*2.9*Qc/(Cp*T).
        fire.mass_offset = positive_result(
            (0.196 * 2.9) * (convective_heat_release.value() / cp / t),
            "entrainment correction",
        )?;
        fire.temperature_scale = positive_result(
            9.1 * t * s_one_third * s_one_third,
            "centerline temperature scale",
        )?;
        fire.velocity_scale =
            positive_result(3.4 * sqrt_g * s_one_third, "centerline velocity scale")?;
        Ok(fire)
    }
    /// Original total heat release; used for flame geometry.
    pub const fn total_heat_release(&self) -> Power {
        self.total_heat_release
    }
    /// Original convective heat release; used for plume buoyancy.
    pub const fn convective_heat_release(&self) -> Power {
        self.convective_heat_release
    }
    /// Original source diameter.
    pub const fn diameter(&self) -> Length {
        self.diameter
    }
    /// Original actual ambient property declaration.
    pub const fn ambient(&self) -> FireAmbientV1 {
        self.ambient
    }
    /// Dimensionless Q*, exactly zero for a dormant source.
    pub const fn dimensionless_heat_release(&self) -> f64 {
        self.dimensionless_heat_release
    }
    /// Nominal mean flame height [m]; exactly zero for a dormant source.
    pub const fn mean_flame_height(&self) -> Length {
        self.flame_height
    }
    /// Virtual origin relative to the fuel surface; absent when no fire exists.
    pub const fn virtual_origin(&self) -> Option<Length> {
        self.virtual_origin
    }
    /// Source card for the retained flame geometry.
    pub const fn card(&self) -> &'static FireModelCardV1 {
        FireCorrelationIdV1::HeskestadFlame.card()
    }

    /// Source-prescribed entrainment at a nonnegative height above the fuel.
    /// The flaming branch is explicit and continuous at L. This is an
    /// unobstructed uniform-ambient plume, not an upper-layer penetration cap.
    pub fn entrainment(
        &self,
        height: Length,
    ) -> Result<FireEstimateV1<HeskestadEntrainmentV1>, FireCorrelationErrorV1> {
        nonnegative_input(height.value(), "height above fuel m")?;
        let (mass, region) = if self.total_heat_release.value() == 0.0 {
            (0.0, HeskestadRegionV1::Dormant)
        } else if height.value() < self.flame_height.value() {
            (
                self.plume_mass(self.flame_height.value())?
                    * (height.value() / self.flame_height.value()),
                HeskestadRegionV1::Flaming,
            )
        } else {
            (
                self.plume_mass(height.value())?,
                HeskestadRegionV1::AboveFlame,
            )
        };
        nonnegative_result(mass, "entrainment mass rate")?;
        Ok(estimate(
            FireCorrelationIdV1::HeskestadEntrainment,
            HeskestadEntrainmentV1 {
                mass_rate: MassFlowRate::new(mass),
                region,
            },
        ))
    }

    /// Unclipped centerline temperature and velocity above the mean flame.
    /// No 900K temperature plateau, velocity maximum, hot-layer correction,
    /// radial Gaussian profile or near-flame interpolation is inferred.
    pub fn centerline(
        &self,
        height: Length,
    ) -> Result<FireEstimateV1<HeskestadCenterlineV1>, FireCorrelationErrorV1> {
        nonnegative_input(height.value(), "centerline height above fuel m")?;
        let mut rise = 0.0;
        let mut velocity = 0.0;
        if self.total_heat_release.value() > 0.0 {
            if height.value() < self.flame_height.value() {
                return Err(FireCorrelationErrorV1::Regime {
                    correlation: FireCorrelationIdV1::HeskestadCenterline,
                    reason: "centerline similarity requires height at or above the mean flame",
                });
            }
            if self.convective_heat_release.value() > 0.0 {
                let z = self.relative_height(height.value())?;
                let z_one_third =
                    positive_result(pow_positive(z, 1.0 / 3.0), "centerline height cube root")?;
                let z_five_thirds =
                    positive_result(z * z_one_third * z_one_third, "centerline height power")?;
                rise = positive_result(
                    self.temperature_scale / z_five_thirds,
                    "centerline temperature rise",
                )?;
                velocity =
                    positive_result(self.velocity_scale / z_one_third, "centerline velocity")?;
            }
        }
        let temperature = positive_result(
            self.ambient.temperature.value() + rise,
            "absolute centerline temperature",
        )?;
        Ok(estimate(
            FireCorrelationIdV1::HeskestadCenterline,
            HeskestadCenterlineV1 {
                temperature_rise: Temperature::new(rise),
                temperature: Temperature::new(temperature),
                velocity: Velocity::new(velocity),
            },
        ))
    }
    fn relative_height(&self, height: f64) -> Result<f64, FireCorrelationErrorV1> {
        let Some(origin) = self.virtual_origin else {
            return Err(FireCorrelationErrorV1::Arithmetic(
                "dormant plume has no virtual origin",
            ));
        };
        positive_result(height - origin.value(), "height above plume virtual origin")
    }
    fn plume_mass(&self, height: f64) -> Result<f64, FireCorrelationErrorV1> {
        if self.convective_heat_release.value() == 0.0 {
            return Ok(0.0);
        }
        let z = self.relative_height(height)?;
        let root = positive_result(pow_positive(z, 1.0 / 3.0), "plume height cube root")?;
        positive_result(
            self.mass_scale * z * root * root + self.mass_offset,
            "above-flame entrainment",
        )
    }
}

/// A single vertical opening or an explicitly justified equivalent opening.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FireOpeningV1 {
    area: Area,
    height: Length,
    ventilation_factor: f64,
}
impl FireOpeningV1 {
    /// Declare area [m²] and height [m]. Zero area denotes a closed opening;
    /// only the ventilation-capacity relation admits that zero-flow limit.
    pub fn new(area: Area, height: Length) -> Result<Self, FireCorrelationErrorV1> {
        nonnegative_input(area.value(), "opening area m2")?;
        positive_input(height.value(), "opening height m")?;
        let ventilation_factor = if area.value() == 0.0 {
            0.0
        } else {
            positive_result(
                area.value() * det::sqrt(height.value()),
                "opening ventilation factor",
            )?
        };
        Ok(Self {
            area,
            height,
            ventilation_factor,
        })
    }
    /// Declared opening area, possibly exactly zero.
    pub const fn area(&self) -> Area {
        self.area
    }
    /// Declared positive opening height.
    pub const fn height(&self) -> Length {
        self.height
    }
    fn require_open(self, correlation: FireCorrelationIdV1) -> Result<(), FireCorrelationErrorV1> {
        if self.area.value() == 0.0 {
            Err(FireCorrelationErrorV1::Regime {
                correlation,
                reason: "this naturally ventilated correlation requires a nonzero opening",
            })
        } else {
            Ok(())
        }
    }
}

/// MQH hot-layer absolute temperature and its rise above the declared ambient.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MqhHotLayerV1 {
    /// Nominal layer temperature rise [K].
    pub temperature_rise: Temperature,
    /// Nominal absolute layer temperature [K].
    pub temperature: Temperature,
}
/// MQH natural-ventilation estimate with an explicit effective wall coefficient.
///
/// NIST SP1169 Eq2.1 uses Q in kW and h in kW/(m² K); the typed API accepts W
/// and W/(m² K), converting both explicitly. Enclosing area includes the floor,
/// ceiling and walls, excluding openings. This pre-flashover nominal relation
/// does not derive wall thermal penetration, infer h, or cap the temperature.
pub fn mqh_hot_layer_temperature(
    total_heat_release: Power,
    ambient: Temperature,
    opening: FireOpeningV1,
    enclosing_area: Area,
    effective_wall_coefficient: HeatTransferCoefficient,
) -> Result<FireEstimateV1<MqhHotLayerV1>, FireCorrelationErrorV1> {
    let id = FireCorrelationIdV1::MqhHotLayer;
    opening.require_open(id)?;
    nonnegative_input(total_heat_release.value(), "total HRR W")?;
    positive_input(ambient.value(), "ambient temperature K")?;
    positive_input(
        enclosing_area.value(),
        "enclosing area excluding openings m2",
    )?;
    positive_input(
        effective_wall_coefficient.value(),
        "effective wall coefficient W/(m2 K)",
    )?;
    let rise = if total_heat_release.value() == 0.0 {
        0.0
    } else {
        let q_kw = positive_result(total_heat_release.value() * 1e-3, "MQH HRR kW")?;
        let h_kw = positive_result(
            effective_wall_coefficient.value() * 1e-3,
            "MQH wall coefficient kW/(m2 K)",
        )?;
        let denominator = positive_result(
            opening.ventilation_factor * h_kw * enclosing_area.value(),
            "MQH denominator",
        )?;
        let group = positive_result((q_kw / denominator) * q_kw, "MQH heat-release group")?;
        positive_result(
            6.85 * pow_positive(group, 1.0 / 3.0),
            "MQH temperature rise",
        )?
    };
    let temperature = positive_result(ambient.value() + rise, "MQH absolute temperature")?;
    Ok(estimate(
        id,
        MqhHotLayerV1 {
            temperature_rise: Temperature::new(rise),
            temperature: Temperature::new(temperature),
        },
    ))
}

/// Thomas's nominal flashover HRR indicator: 7.8 AT + 378 Ao sqrt(Ho) [kW].
/// AT is the entire enclosing surface area excluding openings. The returned
/// threshold is not a prediction of ignition, event time, or actual flashover.
pub fn thomas_flashover_heat_release(
    enclosing_area: Area,
    opening: FireOpeningV1,
) -> Result<FireEstimateV1<Power>, FireCorrelationErrorV1> {
    let id = FireCorrelationIdV1::ThomasFlashover;
    opening.require_open(id)?;
    positive_input(
        enclosing_area.value(),
        "enclosing area excluding openings m2",
    )?;
    let watts = positive_result(
        7800.0 * enclosing_area.value() + 378000.0 * opening.ventilation_factor,
        "Thomas flashover HRR",
    )?;
    Ok(estimate(id, Power::new(watts)))
}
/// Babrauskas's nominal flashover HRR indicator: 750 Ao sqrt(Ho) [kW].
/// Ventilation geometry alone cannot determine whether or when flashover occurs.
pub fn babrauskas_flashover_heat_release(
    opening: FireOpeningV1,
) -> Result<FireEstimateV1<Power>, FireCorrelationErrorV1> {
    let id = FireCorrelationIdV1::BabrauskasFlashover;
    opening.require_open(id)?;
    let watts = positive_result(
        750000.0 * opening.ventilation_factor,
        "Babrauskas flashover HRR",
    )?;
    Ok(estimate(id, Power::new(watts)))
}

/// Nominal ventilation capacity; neither actual oxygen supply nor combustion.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct KawagoeCapacityV1 {
    /// Source approximation for incoming air, 0.5 Ao sqrt(Ho) [kg/s].
    pub air_mass_rate: MassFlowRate,
    /// Capacity for heat release inside the compartment [W].
    pub heat_release_capacity: Power,
}
/// Post-flashover, naturally ventilated capacity from the Kawagoe approximation.
///
/// Qcap=efficiency*heat_per_kg_air*(0.5 Ao sqrt(Ho)). Heat per mass of reacting
/// air is explicit, not an inferred 3 MJ/kg constant. Efficiency is in [0,1].
/// This does not replace actual vent flow, oxygen species balances or a
/// fuel-controlled growth model; external combustion is not included.
pub fn kawagoe_ventilation_capacity(
    opening: FireOpeningV1,
    heat_per_mass_air: SpecificEnergyV1,
    combustion_efficiency: f64,
) -> Result<FireEstimateV1<KawagoeCapacityV1>, FireCorrelationErrorV1> {
    positive_input(
        heat_per_mass_air.value(),
        "heat per mass of reacting air J/kg",
    )?;
    if !combustion_efficiency.is_finite() || !(0.0..=1.0).contains(&combustion_efficiency) {
        return Err(invalid("combustion efficiency", combustion_efficiency));
    }
    let mass = nonnegative_result(
        0.5 * opening.ventilation_factor,
        "Kawagoe incoming air rate",
    )?;
    let capacity = nonnegative_result(
        (mass * combustion_efficiency) * heat_per_mass_air.value(),
        "ventilation HRR capacity",
    )?;
    Ok(estimate(
        FireCorrelationIdV1::KawagoeVentilation,
        KawagoeCapacityV1 {
            air_mass_rate: MassFlowRate::new(mass),
            heat_release_capacity: Power::new(capacity),
        },
    ))
}

fn estimate<T>(correlation: FireCorrelationIdV1, value: T) -> FireEstimateV1<T> {
    FireEstimateV1 { value, correlation }
}
fn invalid(field: &'static str, value: f64) -> FireCorrelationErrorV1 {
    FireCorrelationErrorV1::InvalidInput {
        field,
        bits: value.to_bits(),
    }
}
fn positive_input(value: f64, field: &'static str) -> Result<(), FireCorrelationErrorV1> {
    if value.is_finite() && value > 0.0 {
        Ok(())
    } else {
        Err(invalid(field, value))
    }
}
fn nonnegative_input(value: f64, field: &'static str) -> Result<(), FireCorrelationErrorV1> {
    if value.is_finite() && value >= 0.0 {
        Ok(())
    } else {
        Err(invalid(field, value))
    }
}
fn finite_result(value: f64, stage: &'static str) -> Result<f64, FireCorrelationErrorV1> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(FireCorrelationErrorV1::Arithmetic(stage))
    }
}
fn positive_result(value: f64, stage: &'static str) -> Result<f64, FireCorrelationErrorV1> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(FireCorrelationErrorV1::Arithmetic(stage))
    }
}
fn nonnegative_result(value: f64, stage: &'static str) -> Result<f64, FireCorrelationErrorV1> {
    if value.is_finite() && value >= 0.0 {
        Ok(value)
    } else {
        Err(FireCorrelationErrorV1::Arithmetic(stage))
    }
}
fn pow_positive(value: f64, power: f64) -> f64 {
    det::exp(power * det::ln(value))
}

const CFAST_URL: &str = "https://nvlpubs.nist.gov/nistpubs/TechnicalNotes/NIST.TN.1889v1.pdf";
const MQH_URL: &str = "https://nvlpubs.nist.gov/nistpubs/SpecialPublications/NIST.SP.1169.pdf";
const FAST_URL: &str = "https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=910113";
const FLASHOVER_URL: &str = "https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=912911";

const fn card(
    correlation: FireCorrelationIdV1,
    name: &'static str,
    source: SourceProvenance,
    url: &'static str,
    domain: &'static [&'static str],
    failures: &'static [&'static str],
) -> FireModelCardV1 {
    FireModelCardV1 {
        correlation,
        name,
        version: 1,
        source,
        source_url: url,
        status: FirePredictionStatusV1::Nominal,
        discrepancy: FireDiscrepancyV1::Unquantified,
        domain_assumptions: domain,
        quantitative_validity_domain: None,
        known_failures: failures,
    }
}
const HESKESTAD_FLAME: FireModelCardV1 = card(
    FireCorrelationIdV1::HeskestadFlame,
    "fire.heskestad-flame",
    SourceProvenance {
        citation: "Peacock et al., CFAST Technical Reference Guide, NIST TN 1889v1 (2015), Eqs. 3.16-3.17; Heskestad",
        identifier: "doi:10.6028/NIST.TN.1889v1",
    },
    CFAST_URL,
    &[
        "Buoyancy-dominated unobstructed axisymmetric fire; total HRR defines Q*",
        "Actual homogeneous ambient rho, Cp, absolute T and gravity are explicit",
        "A positive fire must produce a positive fitted mean flame height; zero HRR is explicitly dormant",
    ],
    &[
        "No wind, wall/corner attachment, oxygen limitation or resolved flame chemistry",
        "Mean visible flame height is an empirical intermittency statistic, not a solid flame boundary",
    ],
);
const HESKESTAD_ENTRAINMENT: FireModelCardV1 = card(
    FireCorrelationIdV1::HeskestadEntrainment,
    "fire.heskestad-entrainment",
    SourceProvenance {
        citation: "Peacock et al., CFAST Technical Reference Guide, NIST TN 1889v1 (2015), Eqs. 3.15-3.17; Heskestad",
        identifier: "doi:10.6028/NIST.TN.1889v1",
    },
    CFAST_URL,
    &[
        "Unobstructed axisymmetric plume in a homogeneous ambient",
        "Convective HRR supplies buoyancy; total HRR supplies virtual origin",
        "Below L use the source's m(L) z/L branch, never extrapolate the plume formula through its virtual origin",
    ],
    &[
        "No wall/corner modification or hot-layer penetration correction",
        "No fuel mass source, species chemistry or upper-layer energy balance is implied",
    ],
);
const HESKESTAD_CENTERLINE: FireModelCardV1 = card(
    FireCorrelationIdV1::HeskestadCenterline,
    "fire.heskestad-centerline",
    SourceProvenance {
        citation: "Heskestad centerline similarity expressions in NIST SP1169 Eq.4.1 and NIST TN1889v1 Eqs.3.20-3.21 (unclipped plume branches)",
        identifier: "doi:10.6028/NIST.TN.1889v1",
    },
    CFAST_URL,
    &[
        "Unconfined homogeneous-ambient centerline at or above the mean flame",
        "Uses actual ambient rho, Cp and absolute kelvin temperature",
    ],
    &[
        "Not a flaming/intermittent-region model, hot-upper-layer correction or radial field",
        "No source-specific 900K plateau or velocity cap is applied; near-flame values remain nominal",
    ],
);
const MQH: FireModelCardV1 = card(
    FireCorrelationIdV1::MqhHotLayer,
    "fire.mqh-hot-layer",
    SourceProvenance {
        citation: "Overholt, NIST SP1169 (2014), Eq.2.1; McCaffrey, Quintiere and Harkleroad natural-ventilation correlation",
        identifier: "doi:10.6028/NIST.SP.1169",
    },
    MQH_URL,
    &[
        "Pre-flashover naturally ventilated compartment with a representative hot upper layer",
        "Enclosing area excludes openings; effective wall heat-transfer coefficient is explicitly supplied",
        "Empirical coefficients use kW and kW/(m2 K); public inputs are coherent SI",
    ],
    &[
        "Not forced ventilation, a sealed room, a transient wall solver or a post-flashover model",
        "No unreported calibration domain or confidence interval is inferred from the source's validation examples",
    ],
);
const THOMAS: FireModelCardV1 = card(
    FireCorrelationIdV1::ThomasFlashover,
    "fire.thomas-flashover",
    SourceProvenance {
        citation: "Peacock et al., NIST SP921 (2000), section5.13, Thomas's flashover correlation",
        identifier: "NIST-SP-921-2000-sec5.13",
    },
    FAST_URL,
    &[
        "Natural ventilation through one vertical opening or an explicitly justified equivalent opening",
        "All enclosing surfaces excluding openings contribute to AT",
    ],
    &[
        "A nominal HRR indicator is not a binary physical flashover decision",
        "Fuel arrangement, linings, exposure time and geometry can change actual flashover substantially",
    ],
);
const BABRAUSKAS: FireModelCardV1 = card(
    FireCorrelationIdV1::BabrauskasFlashover,
    "fire.babrauskas-flashover",
    SourceProvenance {
        citation: "Peacock, Reneke, Bukowski and Babrauskas, Defining Flashover for Fire Hazard Calculations, Fire Safety Journal32 (1999), Babrauskas correlation",
        identifier: "NIST-publication-912911",
    },
    FLASHOVER_URL,
    &[
        "Ventilation-factor estimate for compartment fires with a vertical opening",
        "Q_fo=750 Ao sqrt(Ho) in kW is a threshold estimate, not a state evolution law",
    ],
    &[
        "No heat-loss material properties or exposure-time dependence",
        "No actual transition event, ignition mechanism or experimental validation of this implementation",
    ],
);
const KAWAGOE: FireModelCardV1 = card(
    FireCorrelationIdV1::KawagoeVentilation,
    "fire.kawagoe-ventilation-capacity",
    SourceProvenance {
        citation: "Peacock et al., NIST SP921 (2000), section5.14, ventilation limit based on Kawagoe",
        identifier: "NIST-SP-921-2000-sec5.14",
    },
    FAST_URL,
    &[
        "Post-flashover natural ventilation with the source's nominal air-inflow coefficient 0.5",
        "Heat per kilogram of reacting air and combustion efficiency are explicitly authored",
        "One vertical vent or a justified equivalent vent at comparable elevations",
    ],
    &[
        "Not an actual vent-flow or oxygen-species solve and not an early fuel-controlled burning model",
        "Mechanical ventilation, wind, materially different ambient conditions and external burning require separate models",
    ],
);

#[cfg(test)]
mod tests;
