//! Certified chamber-volume functions over a prescribed angle parameter.
//!
//! A chamber family builds one exact-distance chart for the complete closed
//! chamber at each angle.  `fs-query::geometric_moments` then performs the
//! certified spatial quadrature.  Named boundaries and a proven closure
//! statement are mandatory; a mnemonic engine formula cannot mint volume
//! authority for missing seals, ports, rotor flanks, or as-built geometry.

use crate::{MotionError, ProofState};
use fs_exec::Cx;
use fs_geom::{Aabb, Chart};
use fs_ivl::Interval;
use fs_query::geometric_moments;

/// G1 closed-form oracle for one declared ideal nominal Wankel chamber.
///
/// This is an independent comparison value, not a chamber chart and not a
/// certificate for finite seals, rotor-flank construction, ports, recesses,
/// deformation, or as-built geometry.  Its convention is
///
/// `V(alpha) = V_min + sqrt(3)/2 * e * (R1 + R2) * b
///             * (1 - sin(2 alpha / 3 + phase))`,
///
/// with `R1 = generating_radius + housing_parallel_transfer` and
/// `R2 = generating_radius + rotor_parallel_transfer`.  This convention is
/// Peden et al., SAE 2018-01-1452, Eq. (2); callers supply `V_min` so the
/// declared recess/seal convention remains explicit.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct IdealWankelVolumeOracle {
    /// Declared minimum chamber volume, including any declared ideal recess.
    pub minimum_volume_m3: f64,
    /// Shaft eccentricity.
    pub eccentricity_m: f64,
    /// Nominal rotor generating radius.
    pub generating_radius_m: f64,
    /// Parallel transfer applied to the ideal housing.
    pub housing_parallel_transfer_m: f64,
    /// Parallel transfer applied to the ideal rotor.
    pub rotor_parallel_transfer_m: f64,
    /// Axial housing depth.
    pub housing_depth_m: f64,
    /// Phase inside the volume sinusoid; the conventional reference is π/6.
    pub volume_phase_radians: f64,
}

impl IdealWankelVolumeOracle {
    fn validate(self) -> Result<Self, MotionError> {
        for (value, what) in [
            (self.minimum_volume_m3, "ideal Wankel minimum volume"),
            (self.eccentricity_m, "ideal Wankel eccentricity"),
            (self.generating_radius_m, "ideal Wankel generating radius"),
            (
                self.housing_parallel_transfer_m,
                "ideal Wankel housing parallel transfer",
            ),
            (
                self.rotor_parallel_transfer_m,
                "ideal Wankel rotor parallel transfer",
            ),
            (self.housing_depth_m, "ideal Wankel housing depth"),
            (self.volume_phase_radians, "ideal Wankel volume phase"),
        ] {
            if !value.is_finite() {
                return Err(MotionError::NonFiniteInput { what });
            }
        }
        let housing_radius = self.generating_radius_m + self.housing_parallel_transfer_m;
        let rotor_radius = self.generating_radius_m + self.rotor_parallel_transfer_m;
        if self.minimum_volume_m3 < 0.0
            || self.eccentricity_m <= 0.0
            || self.generating_radius_m <= 0.0
            || self.housing_depth_m <= 0.0
            || !housing_radius.is_finite()
            || housing_radius <= 0.0
            || !rotor_radius.is_finite()
            || rotor_radius <= 0.0
        {
            return Err(MotionError::InvalidGeometry {
                what: "ideal Wankel oracle needs nonnegative V_min and positive finite e, R1, R2, and depth",
            });
        }
        Ok(self)
    }

    /// Evaluate the G1 comparison formula at a finite shaft angle.
    pub fn volume_at(self, shaft_angle_radians: f64) -> Result<f64, MotionError> {
        let valid = self.validate()?;
        if !shaft_angle_radians.is_finite() {
            return Err(MotionError::NonFiniteInput {
                what: "ideal Wankel shaft angle",
            });
        }
        let housing_radius = valid.generating_radius_m + valid.housing_parallel_transfer_m;
        let rotor_radius = valid.generating_radius_m + valid.rotor_parallel_transfer_m;
        let amplitude = 0.5
            * fs_math::det::sqrt(3.0)
            * valid.eccentricity_m
            * (housing_radius + rotor_radius)
            * valid.housing_depth_m;
        let phase = (2.0 / 3.0) * shaft_angle_radians + valid.volume_phase_radians;
        let volume = valid.minimum_volume_m3 + amplitude * (1.0 - fs_math::det::sin(phase));
        if !volume.is_finite() || volume < 0.0 {
            return Err(MotionError::InvalidGeometry {
                what: "ideal Wankel formula produced a non-finite or negative volume",
            });
        }
        Ok(volume)
    }
}

/// Content of a named closed chamber model.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ChamberDefinition {
    name: String,
    boundary_names: Vec<String>,
    seal_convention: String,
    closure: ProofState,
}

impl ChamberDefinition {
    /// Construct a chamber definition.  At least two distinct nonempty named
    /// boundaries and a nonempty seal/closure convention are required.
    pub fn new(
        name: impl Into<String>,
        boundary_names: Vec<String>,
        seal_convention: impl Into<String>,
        closure: ProofState,
    ) -> Result<Self, MotionError> {
        let name = name.into();
        let seal_convention = seal_convention.into();
        if name.trim().is_empty() || seal_convention.trim().is_empty() {
            return Err(MotionError::InvalidGeometry {
                what: "chamber name and seal convention must be nonempty",
            });
        }
        if boundary_names.len() < 2 || boundary_names.iter().any(|name| name.trim().is_empty()) {
            return Err(MotionError::InvalidGeometry {
                what: "a chamber needs at least two nonempty named boundaries",
            });
        }
        for (index, boundary) in boundary_names.iter().enumerate() {
            if boundary_names[..index].contains(boundary) {
                return Err(MotionError::InvalidGeometry {
                    what: "chamber boundary names must be unique",
                });
            }
        }
        Ok(Self {
            name,
            boundary_names,
            seal_convention,
            closure,
        })
    }

    /// Stable chamber name.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Names of every chart contributing to the closed boundary.
    #[must_use]
    pub fn boundary_names(&self) -> &[String] {
        &self.boundary_names
    }

    /// Declared ideal/finite seal and clearance convention.
    #[must_use]
    pub fn seal_convention(&self) -> &str {
        &self.seal_convention
    }

    /// Proof state of complete boundary closure under that convention.
    #[must_use]
    pub fn closure(&self) -> ProofState {
        self.closure
    }
}

/// Additive abstract-volume uncertainty beyond the spatial quadrature band,
/// in cubic metres.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct ChamberVolumeErrors {
    /// Chart/representation conversion error.
    pub chart_conversion_m3: f64,
    /// Intended-motion versus constructed-tube error.
    pub motion_tube_m3: f64,
    /// Seal/closure-model error under the declared boundary convention.
    pub boundary_closure_m3: f64,
    /// Other explicitly declared model-form error.
    pub model_form_m3: f64,
}

fn add_nonnegative_upper(left: f64, right: f64) -> f64 {
    let sum = left + right;
    if sum.total_cmp(&0.0).is_eq() {
        0.0
    } else {
        sum.next_up()
    }
}

impl ChamberVolumeErrors {
    fn validate(self) -> Result<Self, MotionError> {
        for (value, what) in [
            (
                self.chart_conversion_m3,
                "chamber chart-conversion volume error",
            ),
            (self.motion_tube_m3, "chamber motion-tube volume error"),
            (
                self.boundary_closure_m3,
                "chamber boundary-closure volume error",
            ),
            (self.model_form_m3, "chamber model-form volume error"),
        ] {
            if !value.is_finite() || value < 0.0 {
                return Err(MotionError::InvalidEvidence { what });
            }
        }
        Ok(self)
    }

    /// Outward-rounded total error radius in cubic metres.
    pub fn total_upper(self) -> Result<f64, MotionError> {
        let valid = self.validate()?;
        let total = add_nonnegative_upper(
            add_nonnegative_upper(valid.chart_conversion_m3, valid.motion_tube_m3),
            add_nonnegative_upper(valid.boundary_closure_m3, valid.model_form_m3),
        );
        if !total.is_finite() {
            return Err(MotionError::InvalidEvidence {
                what: "chamber total volume error must be finite",
            });
        }
        Ok(total)
    }
}

/// Builds the exact-distance chart for a complete chamber at one angle.
///
/// The returned chart must represent the region bounded by every name in the
/// supplied definition under its seal convention.  `geometric_moments`
/// independently refuses any chart that does not advertise and evidence the
/// global `ExactDistance` theorem.
pub trait ChamberChartFamily: Send + Sync {
    /// Concrete exact-distance chamber chart.
    type Chamber: Chart;

    /// Build the chamber at `angle_radians`.
    fn chart_at(
        &self,
        definition: &ChamberDefinition,
        angle_radians: f64,
        cx: &Cx<'_>,
    ) -> Result<Self::Chamber, MotionError>;
}

/// Certified `V(theta)` result.
#[derive(Debug, Clone, PartialEq)]
pub struct ChamberVolumeReceipt {
    /// Named chamber definition used to construct the chart.
    pub definition: ChamberDefinition,
    /// Angle in radians.
    pub angle_radians: f64,
    /// Final volume enclosure after all additive error inflation.
    pub volume_m3: Interval,
    /// Raw `fs-query` quadrature enclosure.
    pub spatial_quadrature_m3: Interval,
    /// Width of the spatial quadrature band (logged, not double-counted).
    pub spatial_quadrature_band_m3: f64,
    /// Additional named error components.
    pub errors: ChamberVolumeErrors,
    /// Requested maximum spatial cell spacing in metres.
    pub h_m: f64,
    /// Cells certified wholly inside.
    pub sure_cells: u64,
    /// Boundary-straddling cells retained in the upper bound.
    pub band_cells: u64,
}

/// Evaluate a certified chamber-volume function at one angle.
pub fn chamber_volume_at<F: ChamberChartFamily>(
    family: &F,
    definition: &ChamberDefinition,
    angle_radians: f64,
    integration_domain: &Aabb,
    h_m: f64,
    errors: ChamberVolumeErrors,
    cx: &Cx<'_>,
) -> Result<ChamberVolumeReceipt, MotionError> {
    if !angle_radians.is_finite() {
        return Err(MotionError::NonFiniteInput {
            what: "chamber angle",
        });
    }
    if definition.closure() != ProofState::Proven {
        return Err(MotionError::InvalidEvidence {
            what: "chamber boundary closure must be Proven before volume authority",
        });
    }
    let errors = errors.validate()?;
    cx.checkpoint().map_err(|_| MotionError::Cancelled)?;
    let chart = family.chart_at(definition, angle_radians, cx)?;
    let moments = geometric_moments(&chart, integration_domain, h_m, cx)?;
    let raw = Interval::new(moments.volume.lo, moments.volume.hi);
    let extra = errors.total_upper()?;
    let lower = (raw.lo() - extra).next_down().max(0.0);
    let upper = (raw.hi() + extra).next_up();
    if !(lower.is_finite() && upper.is_finite() && lower <= upper) {
        return Err(MotionError::InvalidEvidence {
            what: "inflated chamber volume must be a finite ordered enclosure",
        });
    }
    Ok(ChamberVolumeReceipt {
        definition: definition.clone(),
        angle_radians,
        volume_m3: Interval::new(lower, upper),
        spatial_quadrature_m3: raw,
        spatial_quadrature_band_m3: raw.width(),
        errors,
        h_m: moments.h,
        sure_cells: moments.sure_cells,
        band_cells: moments.band_cells,
    })
}

/// Reusable certified chamber-volume function with fixed geometry,
/// quadrature, and error policy.
#[derive(Debug, Clone)]
pub struct ChamberVolumeFunction<F> {
    family: F,
    definition: ChamberDefinition,
    integration_domain: Aabb,
    h_m: f64,
    errors: ChamberVolumeErrors,
}

impl<F: ChamberChartFamily> ChamberVolumeFunction<F> {
    /// Construct a reusable `V(theta)` evaluator.
    pub fn new(
        family: F,
        definition: ChamberDefinition,
        integration_domain: Aabb,
        h_m: f64,
        errors: ChamberVolumeErrors,
    ) -> Result<Self, MotionError> {
        if !(integration_domain.is_finite()
            && integration_domain.max.x > integration_domain.min.x
            && integration_domain.max.y > integration_domain.min.y
            && integration_domain.max.z > integration_domain.min.z)
        {
            return Err(MotionError::InvalidGeometry {
                what: "chamber integration domain must be finite with positive 3-D extent",
            });
        }
        if !h_m.is_finite() || h_m <= 0.0 {
            return Err(MotionError::InvalidConfiguration {
                what: "chamber quadrature spacing h_m must be finite and positive",
            });
        }
        let errors = errors.validate()?;
        Ok(Self {
            family,
            definition,
            integration_domain,
            h_m,
            errors,
        })
    }

    /// Evaluate `V(theta)` with the fixed certified policy.
    pub fn at(&self, angle_radians: f64, cx: &Cx<'_>) -> Result<ChamberVolumeReceipt, MotionError> {
        chamber_volume_at(
            &self.family,
            &self.definition,
            angle_radians,
            &self.integration_domain,
            self.h_m,
            self.errors,
            cx,
        )
    }

    /// The immutable named boundary definition.
    #[must_use]
    pub fn definition(&self) -> &ChamberDefinition {
        &self.definition
    }
}
