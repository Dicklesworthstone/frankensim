//! Selected average net-ligament tensile limit for a finite-width centre crack.
//!
//! Force equilibrium gives sigma_net = sigma_gross / (1 - 2a/W), where W is
//! full plate width and a is physical half crack length. The original typed
//! secant geometry and tensile-positive gross nominal-stress convention are
//! required. The caller selects a sourced tensile limit at the same fixed
//! external conditions as the Paris/Walker fit; no allowable is inferred.
//!
//! This is an additional model limit, not a crack-tip stress, toughness,
//! residual-strength certification, necking or plastic-collapse model. Reaching
//! it refuses growth without declaring rupture or changing mechanical stiffness.
//! AFGROW DTD Handbook §2.1.2 distinguishes net-section and fracture criteria:
//! <https://www.afgrow.net/applications/DTDHandbook/sections/page2_1_2.aspx>.

use std::sync::Arc;

use super::{
    CrackGeometry, CrackGrowthError, CrackGrowthForecast, CrackGrowthStep, ParisCrackGrowth,
    ParisLaw, StressCycle, TabulatedCrackGeometry,
};
use crate::strength::IsotropicTensileLimit;

/// Immutable tensile source and exact original geometry/calibration binding.
#[derive(Clone, Debug, PartialEq)]
pub struct CentreCrackTensileLimit {
    limit: Arc<IsotropicTensileLimit>,
    law: Arc<ParisLaw>,
    geometry: Arc<TabulatedCrackGeometry>,
    plate_width_m: f64,
}

impl CentreCrackTensileLimit {
    /// Select the average net-section criterion for this sourced growth model.
    /// Constant Y, arbitrary tables and caller-declared Paris fits refuse.
    pub fn try_new(
        growth: &ParisCrackGrowth,
        limit: IsotropicTensileLimit,
    ) -> Result<Self, CentreCrackTensileError> {
        let source = growth
            .law()
            .source()
            .ok_or(CentreCrackTensileError::Invalid(
                "a net-section guard requires the original source-backed Paris/Walker fit",
            ))?;
        let resolved = limit.resolved();
        if !source.matches_material_state(
            resolved.card_identity(),
            resolved.material(),
            &resolved.query(),
        ) || resolved.query_point().len()
            != source
                .query()
                .axes()
                .keys()
                .filter(|axis| {
                    *axis != &source.request().stress_ratio_axis
                        && *axis != &source.request().delta_k_squared_axis
                })
                .count()
        {
            return Err(CentreCrackTensileError::Invalid(
                "tensile and growth sources must share their original card, material and typed operating query",
            ));
        }
        let CrackGeometry::Tabulated(geometry) = growth.geometry() else {
            return Err(CentreCrackTensileError::Invalid(
                "typed centre-crack geometry is required",
            ));
        };
        let centre =
            geometry
                .finite_width_centre_crack()
                .ok_or(CentreCrackTensileError::Invalid(
                    "an arbitrary Y(a) table does not establish net ligament width",
                ))?;
        Ok(Self {
            limit: Arc::new(limit),
            law: growth.law.clone(),
            geometry: geometry.clone(),
            plate_width_m: centre.plate_width_m(),
        })
    }

    pub fn limit(&self) -> &IsotropicTensileLimit {
        &self.limit
    }
    pub fn geometry(&self) -> &TabulatedCrackGeometry {
        &self.geometry
    }
    pub fn law(&self) -> &ParisLaw {
        &self.law
    }
    pub fn plate_width_m(&self) -> f64 {
        self.plate_width_m
    }

    /// Allocation-free average ligament assessment within the retained table.
    /// Compression is outside this tensile-only crack model, not clipped away.
    pub fn assess(
        &self,
        half_length_m: f64,
        gross_stress_pa: f64,
    ) -> Result<CentreCrackTensileAssessment, CentreCrackTensileError> {
        let (minimum, maximum) = self.geometry.length_bounds();
        if !half_length_m.is_finite() || half_length_m < minimum || half_length_m > maximum {
            return Err(CrackGrowthError::CrackLengthOutsideDomain {
                length_m: half_length_m,
                minimum,
                maximum,
            }
            .into());
        }
        if !gross_stress_pa.is_finite() || gross_stress_pa < 0.0 {
            return Err(CrackGrowthError::StressOutsideDomain {
                minimum_pa: gross_stress_pa,
                maximum_pa: gross_stress_pa,
            }
            .into());
        }
        // Divide first: 2a can overflow even though a/W is within support.
        let cracked_fraction = 2.0 * (half_length_m / self.plate_width_m);
        let ligament_fraction = 1.0 - cracked_fraction;
        let net_stress_pa = gross_stress_pa / ligament_fraction;
        let utilization = net_stress_pa / self.limit.limit_pa();
        if cracked_fraction <= 0.0
            || ligament_fraction <= 0.0
            || !net_stress_pa.is_finite()
            || !utilization.is_finite()
            || (gross_stress_pa > 0.0 && utilization == 0.0)
        {
            return Err(CentreCrackTensileError::Arithmetic);
        }
        Ok(CentreCrackTensileAssessment {
            half_length_m,
            gross_stress_pa,
            ligament_fraction,
            net_stress_pa,
            utilization,
        })
    }

    /// Refuse equality or exceedance without changing any growth state.
    pub fn check(
        &self,
        half_length_m: f64,
        gross_stress_pa: f64,
    ) -> Result<CentreCrackTensileAssessment, CentreCrackTensileError> {
        let assessment = self.assess(half_length_m, gross_stress_pa)?;
        if assessment.reaches_limit() {
            Err(CentreCrackTensileError::LimitReached(assessment))
        } else {
            Ok(assessment)
        }
    }

    /// Integrate one block, publishing only if its peak is below the selected
    /// limit at both ends. Nondecreasing a makes the endpoint the maximum net
    /// demand throughout a constant-amplitude block. Existing fit/geometry,
    /// exposure and cancellation limits still apply.
    pub fn advance_with_cancel(
        &self,
        growth: &mut ParisCrackGrowth,
        cycle: StressCycle,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<CrackGrowthStep, CentreCrackTensileError> {
        if cancelled() {
            return Err(CrackGrowthError::Cancelled.into());
        }
        self.check_binding(growth)?;
        self.check(growth.state().crack_length_m, cycle.maximum_pa)?;
        let step = growth.preview(cycle, cancelled)?;
        self.check(step.after.crack_length_m, cycle.maximum_pa)?;
        if cancelled() {
            return Err(CrackGrowthError::Cancelled.into());
        }
        growth.state = step.after;
        Ok(step)
    }

    pub fn advance(
        &self,
        growth: &mut ParisCrackGrowth,
        cycle: StressCycle,
    ) -> Result<CrackGrowthStep, CentreCrackTensileError> {
        self.advance_with_cancel(growth, cycle, &mut || false)
    }

    /// Conditional cycles to a = (W/2)(1 - sigma_max/limit), using the existing
    /// analytical inverse. A target beyond geometry or calibration support
    /// refuses; no extrapolated lifetime, toughness or actual rupture is inferred.
    /// The forecast is read-only and requires a strictly positive stress range.
    pub fn cycles_to_limit_with_cancel(
        &self,
        growth: &ParisCrackGrowth,
        minimum_pa: f64,
        maximum_pa: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<CrackGrowthForecast, CentreCrackTensileError> {
        if cancelled() {
            return Err(CrackGrowthError::Cancelled.into());
        }
        self.check_binding(growth)?;
        self.check(growth.state().crack_length_m, maximum_pa)?;
        let ratio = maximum_pa / self.limit.limit_pa();
        let target = (0.5 * self.plate_width_m) * (1.0 - ratio);
        if !target.is_finite() || target <= 0.0 || (maximum_pa > 0.0 && ratio == 0.0) {
            return Err(CentreCrackTensileError::Arithmetic);
        }
        Ok(growth.cycles_to_length_with_cancel(minimum_pa, maximum_pa, target, cancelled)?)
    }

    pub fn cycles_to_limit(
        &self,
        growth: &ParisCrackGrowth,
        minimum_pa: f64,
        maximum_pa: f64,
    ) -> Result<CrackGrowthForecast, CentreCrackTensileError> {
        self.cycles_to_limit_with_cancel(growth, minimum_pa, maximum_pa, &mut || false)
    }

    fn check_binding(&self, growth: &ParisCrackGrowth) -> Result<(), CentreCrackTensileError> {
        let same_geometry = matches!(growth.geometry(), CrackGeometry::Tabulated(table)
            if Arc::ptr_eq(table, &self.geometry) || table == &self.geometry);
        if !same_geometry || !(Arc::ptr_eq(&growth.law, &self.law) || growth.law() == self.law()) {
            return Err(CentreCrackTensileError::Invalid(
                "net-section guard belongs to a different original geometry or growth calibration",
            ));
        }
        Ok(())
    }
}

/// Average nominal-section observation; no local crack-tip stress is claimed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CentreCrackTensileAssessment {
    half_length_m: f64,
    gross_stress_pa: f64,
    ligament_fraction: f64,
    net_stress_pa: f64,
    utilization: f64,
}
impl CentreCrackTensileAssessment {
    pub fn half_length_m(&self) -> f64 {
        self.half_length_m
    }
    pub fn gross_stress_pa(&self) -> f64 {
        self.gross_stress_pa
    }
    pub fn ligament_fraction(&self) -> f64 {
        self.ligament_fraction
    }
    pub fn net_stress_pa(&self) -> f64 {
        self.net_stress_pa
    }
    pub fn utilization(&self) -> f64 {
        self.utilization
    }
    pub fn reaches_limit(&self) -> bool {
        self.utilization >= 1.0
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum CentreCrackTensileError {
    Invalid(&'static str),
    Growth(CrackGrowthError),
    LimitReached(CentreCrackTensileAssessment),
    Arithmetic,
}
impl From<CrackGrowthError> for CentreCrackTensileError {
    fn from(error: CrackGrowthError) -> Self {
        Self::Growth(error)
    }
}
impl core::fmt::Display for CentreCrackTensileError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            f,
            "centre-crack net-section tensile assessment refused: {self:?}"
        )
    }
}
impl std::error::Error for CentreCrackTensileError {}
