//! Bounded nominal compression references from an original stress/strain table.
//!
//! The caller may explicitly author linear interpolation between discrete
//! observations. That choice never changes the source claim's interpolation
//! policy or receipts. No zero-strain point is invented: a coupon starts at the
//! first admitted observation and reports loading work relative to that state.
//! Positive shortening, engineering compression strain and nominal stress use
//! original length/area. Unloading, holds, recovery, creep, failure, lateral
//! contraction, dissipated heat and recoverable energy are not supplied.

use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::{
    ClaimSelection, InterpolationPolicy, MatDbError, MaterialAnswer, MaterialCard, PropertyKey,
    PropertyValue, QueryPoint,
};
use fs_qty::{Density, Dims, Pressure, QuantitySpec};

const SPEED: Dims = Dims([1, 0, -1, 0, 0, 0]);
const RATE: Dims = Dims([0, 0, -1, 0, 0, 0]);

/// Source interpolation and an authored interpolation of discrete observations
/// have different authority, even when they produce identical numerical values.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CompressionReferenceInterpolation {
    SourceInterpolationOnly,
    PiecewiseLinearBetweenObservations,
}

#[derive(Clone, Debug, PartialEq)]
pub struct CompressionReferenceRequest {
    pub density: PropertyKey,
    pub stress: PropertyKey,
    /// Exact dimensionless engineering-compression coordinate, positive in compression.
    pub strain_axis: String,
    /// Exact source m/s coordinate, enforced against each shortening/time ramp.
    pub crosshead_speed_axis: Option<String>,
    pub interpolation: CompressionReferenceInterpolation,
}

#[derive(Clone, Debug, PartialEq)]
pub enum CompressionReferenceError {
    Source(MatDbError),
    Invalid(&'static str),
}
impl core::fmt::Display for CompressionReferenceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Source(error) => error.fmt(f),
            Self::Invalid(message) => write!(f, "compression reference: {message}"),
        }
    }
}
impl std::error::Error for CompressionReferenceError {}
impl From<MatDbError> for CompressionReferenceError {
    fn from(error: MatDbError) -> Self {
        Self::Source(error)
    }
}
type Result<T> = core::result::Result<T, CompressionReferenceError>;
fn require(condition: bool, message: &'static str) -> Result<()> {
    if condition {
        Ok(())
    } else {
        Err(CompressionReferenceError::Invalid(message))
    }
}
fn positive(value: f64, message: &'static str) -> Result<f64> {
    require(value.is_finite() && value > 0.0, message)?;
    Ok(value)
}
fn scalar_key(key: &PropertyKey, dims: Dims) -> bool {
    key.quantity() == QuantitySpec::dimensional(dims)
        && key.hardness_test().is_none()
        && key.elastic_component().is_none()
        && key.strain_component().is_none()
        && key.stress_component().is_none()
}
fn query(
    card: &MaterialCard,
    key: &PropertyKey,
    point: &QueryPoint,
    selection: ClaimSelection,
) -> Result<MaterialAnswer> {
    let answer = match selection {
        ClaimSelection::Policy(policy) => card.claims().query_typed(key, point, policy)?,
        ClaimSelection::Pinned(pin) => card.claims().query_pinned_typed(key, point, pin)?,
    };
    card.claims().verify_receipt(&answer.receipt)?;
    Ok(answer)
}

/// An immutable, source-bound reference, with every original knot answer retained.
#[derive(Clone, Debug)]
pub struct ResolvedCompressionReference {
    request: CompressionReferenceRequest,
    card_identity: ContentHash,
    density: MaterialAnswer,
    observations: Vec<(f64, MaterialAnswer)>,
    slopes: Vec<f64>,
    crosshead_speed_m_s: Option<f64>,
    identity: ContentHash,
}
impl ResolvedCompressionReference {
    /// Resolve a single nominal compression curve and strain-independent reference
    /// density. The supplied point must itself be source-supported. Every original
    /// knot is queried under the same selection and frozen external context.
    /// Scalar strengths, moduli, true stresses, surfaces and rate-dependent curves
    /// cannot substitute. Discrete observations require explicit model authorship.
    pub fn resolve(
        card: &MaterialCard,
        point: &QueryPoint,
        request: &CompressionReferenceRequest,
        density_selection: ClaimSelection,
        stress_selection: ClaimSelection,
    ) -> Result<Self> {
        require(
            scalar_key(&request.density, Density::DIMS)
                && scalar_key(&request.stress, Pressure::DIMS),
            "density and nominal compression stress need complete scalar dimension-only keys",
        )?;
        require(
            matches!(
                request.stress.name(),
                "compressive-stress"
                    | "compressive_stress"
                    | "compression-stress"
                    | "compression_stress"
            ),
            "select nominal compression stress, not modulus, strength or true stress",
        )?;
        require(
            point.axis_quantities().get(&request.strain_axis)
                == Some(&QuantitySpec::dimensional(Dims::NONE)),
            "compression strain must be an explicit dimensionless query coordinate",
        )?;
        let density = query(card, &request.density, point, density_selection)?;
        let stress = query(card, &request.stress, point, stress_selection)?;
        let density_claim = card.claims().claim(density.receipt.selected).ok_or(
            CompressionReferenceError::Invalid("missing selected density"),
        )?;
        require(
            matches!(density_claim.value, PropertyValue::Scalar { .. })
                && density_claim.interpolation == InterpolationPolicy::ConstantWithinValidity
                && !density_claim
                    .validity
                    .bounds()
                    .contains_key(&request.strain_axis),
            "reference density must be a strain-independent validity-wide scalar",
        )?;
        positive(
            density.evidence.value.value,
            "positive finite reference density required",
        )?;
        let claim = card.claims().claim(stress.receipt.selected).ok_or(
            CompressionReferenceError::Invalid("missing selected compression curve"),
        )?;
        let knots = match &claim.value {
            PropertyValue::Curve {
                abscissa,
                abscissa_dims,
                knots,
                dims,
            } if abscissa == &request.strain_axis
                && *abscissa_dims == Dims::NONE
                && *dims == Pressure::DIMS =>
            {
                knots
            }
            _ => {
                return Err(CompressionReferenceError::Invalid(
                    "an actual engineering compression-strain table is required",
                ));
            }
        };
        require(
            claim.interpolation == InterpolationPolicy::LinearInside
                || (claim.interpolation == InterpolationPolicy::TabulatedOnly
                    && request.interpolation
                        == CompressionReferenceInterpolation::PiecewiseLinearBetweenObservations),
            "discrete source points require explicit piecewise-linear reference authorship",
        )?;
        require(
            (2..=128).contains(&knots.len()),
            "compression reference needs 2..128 source knots",
        )?;

        // A fixed context does not silently own a frequency/strain-rate protocol.
        // Crosshead speed is the one admitted rate and is checked during loading.
        for selected in [density_claim, claim] {
            for (axis, quantity) in selected.validity.axis_quantities() {
                require(
                    quantity.dims() != RATE,
                    "frequency or strain-rate applicability needs a separate physical owner",
                )?;
                if quantity.dims() == SPEED {
                    require(
                        request.crosshead_speed_axis.as_deref() == Some(axis.as_str()),
                        "every source speed constraint needs its explicit crosshead-speed owner",
                    )?;
                }
            }
        }
        let crosshead_speed_m_s = if let Some(axis) = &request.crosshead_speed_axis {
            require(
                claim.validity.axis_quantities().get(axis)
                    == Some(&QuantitySpec::dimensional(SPEED))
                    && claim.validity.bounds().contains_key(axis),
                "crosshead-speed owner must name the selected stress claim's typed m/s coordinate",
            )?;
            Some(positive(
                *point
                    .axes()
                    .get(axis)
                    .ok_or(CompressionReferenceError::Invalid("missing source speed"))?,
                "positive finite source speed required",
            )?)
        } else {
            None
        };

        let mut observations = Vec::with_capacity(knots.len());
        for &(strain, _) in knots {
            require(
                strain.is_finite() && (0.0..1.0).contains(&strain),
                "engineering compression strains must be in [0,1)",
            )?;
            let at = point.clone().with(&request.strain_axis, strain)?;
            // Preserve all original knots. A narrowed validity domain is not
            // silently repaired by dropping inconvenient source observations.
            let answer = query(card, &request.stress, &at, stress_selection)?;
            require(
                answer.receipt.selected == stress.receipt.selected,
                "the selected source claim must remain unchanged over its knots",
            )?;
            require(
                answer.evidence.value.value.is_finite() && answer.evidence.value.value >= 0.0,
                "nominal compression stress must be finite and nonnegative",
            )?;
            observations.push((strain, answer));
        }
        let mut slopes = Vec::with_capacity(knots.len() - 1);
        for pair in observations.windows(2) {
            let dx = positive(
                pair[1].0 - pair[0].0,
                "source strains must strictly increase",
            )?;
            let slope = (pair[1].1.evidence.value.value - pair[0].1.evidence.value.value) / dx;
            require(
                slope.is_finite() && slope >= 0.0,
                "compression reference requires nondecreasing stress with finite segment tangents",
            )?;
            slopes.push(slope);
        }
        let mut h = DomainHasher::new("org.frankensim.fs-material.compression-reference.v1");
        h.update(card.content_hash().as_bytes());
        h.update(density.receipt.content_hash().as_bytes());
        h.update(stress.receipt.content_hash().as_bytes());
        h.update(&[match request.interpolation {
            CompressionReferenceInterpolation::SourceInterpolationOnly => 0,
            CompressionReferenceInterpolation::PiecewiseLinearBetweenObservations => 1,
        }]);
        for (strain, answer) in &observations {
            h.update(&strain.to_bits().to_le_bytes());
            h.update(answer.receipt.content_hash().as_bytes());
        }
        if let Some(axis) = &request.crosshead_speed_axis {
            h.update(&(axis.len() as u64).to_le_bytes());
            h.update(axis.as_bytes());
        }
        Ok(Self {
            request: request.clone(),
            card_identity: card.content_hash(),
            density,
            observations,
            slopes,
            crosshead_speed_m_s,
            identity: h.finalize(),
        })
    }
    pub const fn request(&self) -> &CompressionReferenceRequest {
        &self.request
    }
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
    pub const fn card_identity(&self) -> ContentHash {
        self.card_identity
    }
    pub const fn density_answer(&self) -> &MaterialAnswer {
        &self.density
    }
    pub fn observations(&self) -> &[(f64, MaterialAnswer)] {
        &self.observations
    }
    pub fn strain_domain(&self) -> [f64; 2] {
        [
            self.observations[0].0,
            self.observations[self.observations.len() - 1].0,
        ]
    }
    fn sample(&self, strain: f64) -> Result<(f64, f64)> {
        let [lo, hi] = self.strain_domain();
        require(
            strain.is_finite() && strain >= lo && strain <= hi,
            "compression strain outside original source span; no extrapolation or invented origin",
        )?;
        let i = self
            .observations
            .partition_point(|p| p.0 <= strain)
            .saturating_sub(1)
            .min(self.slopes.len() - 1);
        let stress = self.observations[i].1.evidence.value.value
            + self.slopes[i] * (strain - self.observations[i].0);
        require(
            stress.is_finite() && stress >= 0.0,
            "unrepresentable model stress",
        )?;
        Ok((stress, self.slopes[i]))
    }
    fn incremental_work_density(&self, from: f64, to: f64) -> Result<f64> {
        // Integrate each crossed segment directly, avoiding subtraction of large
        // cumulative work values for a small late increment.
        let mut work = 0.0;
        for pair in self.observations.windows(2) {
            let a = from.max(pair[0].0);
            let b = to.min(pair[1].0);
            if b > a {
                let ya = self.sample(a)?.0;
                let yb = self.sample(b)?.0;
                work += (b - a) * (0.5 * ya + 0.5 * yb);
            }
        }
        require(
            work.is_finite() && work >= 0.0,
            "unrepresentable model loading work",
        )?;
        Ok(work)
    }
    /// Begin at the first observed compression, with zero subsequent model work.
    /// Geometry is authored original geometry, not an inferred test specimen.
    pub fn coupon(&self, length_m: f64, area_m2: f64) -> Result<CompressionReferenceCoupon> {
        positive(length_m, "positive finite original length required")?;
        positive(area_m2, "positive finite original area required")?;
        let volume_m3 = positive(length_m * area_m2, "unrepresentable reference volume")?;
        let mass_kg = positive(
            volume_m3 * self.density.evidence.value.value,
            "unrepresentable reference mass",
        )?;
        let initial_shortening_m = self.strain_domain()[0] * length_m;
        require(
            initial_shortening_m.is_finite(),
            "unrepresentable initial shortening",
        )?;
        let mut h = DomainHasher::new("org.frankensim.fs-material.compression-reference-coupon.v1");
        h.update(self.identity.as_bytes());
        h.update(&length_m.to_bits().to_le_bytes());
        h.update(&area_m2.to_bits().to_le_bytes());
        Ok(CompressionReferenceCoupon {
            material: self.clone(),
            length_m,
            area_m2,
            volume_m3,
            mass_kg,
            initial_shortening_m,
            strain: self.strain_domain()[0],
            elapsed_s: 0.0,
            work_j: 0.0,
            identity: h.finalize(),
            last: None,
        })
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct CompressionReferenceResponse {
    pub elapsed_s: f64,
    pub shortening_m: f64,
    pub compressive_strain: f64,
    pub nominal_stress_pa: f64,
    pub force_n: f64,
    /// One-sided segment stiffness; right side at interior knots, left at the end.
    pub tangent_n_m: f64,
    /// Model loading work SINCE the first observed compression, not from zero.
    pub model_loading_work_j: f64,
    pub incremental_model_loading_work_j: f64,
    pub identity: ContentHash,
}

#[derive(Clone, Debug)]
pub struct CompressionReferenceCoupon {
    material: ResolvedCompressionReference,
    length_m: f64,
    area_m2: f64,
    volume_m3: f64,
    mass_kg: f64,
    initial_shortening_m: f64,
    strain: f64,
    elapsed_s: f64,
    work_j: f64,
    identity: ContentHash,
    last: Option<CompressionReferenceResponse>,
}
impl CompressionReferenceCoupon {
    pub const fn initial_shortening_m(&self) -> f64 {
        self.initial_shortening_m
    }
    pub const fn mass_kg(&self) -> f64 {
        self.mass_kg
    }
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
    pub fn last_response(&self) -> Option<&CompressionReferenceResponse> {
        self.last.as_ref()
    }

    /// Advance a strictly increasing shortening ramp, preserving source speed
    /// when stated. Rejected steps leave every history field unchanged.
    pub fn advance_shortening(
        &mut self,
        shortening_m: f64,
        duration_s: f64,
    ) -> Result<CompressionReferenceResponse> {
        positive(duration_s, "finite positive ramp duration required")?;
        let strain = shortening_m / self.length_m;
        require(
            shortening_m.is_finite() && strain.is_finite() && strain > self.strain,
            "compression reference forbids holds, unloading and reversal",
        )?;
        let (stress, tangent) = self.material.sample(strain)?;
        if let Some(speed) = self.material.crosshead_speed_m_s {
            let actual = (strain - self.strain) * self.length_m / duration_s;
            require(
                actual.is_finite() && (actual / speed - 1.0).abs() <= 1e-10,
                "imposed crosshead speed differs from the selected source protocol",
            )?;
        }
        let elapsed_s = self.elapsed_s + duration_s;
        require(
            elapsed_s.is_finite() && elapsed_s > self.elapsed_s,
            "unrepresentable elapsed time",
        )?;
        let force_n = stress * self.area_m2;
        let tangent_n_m = tangent * self.area_m2 / self.length_m;
        let increment = self
            .material
            .incremental_work_density(self.strain, strain)?
            * self.volume_m3;
        let work = self.work_j + increment;
        require(
            [force_n, tangent_n_m, increment, work]
                .iter()
                .all(|x| x.is_finite() && *x >= 0.0)
                && (increment == 0.0 || work > self.work_j),
            "unrepresentable force, stiffness or loading work",
        )?;
        let mut h = DomainHasher::new("org.frankensim.fs-material.compression-reference-step.v1");
        h.update(self.identity.as_bytes());
        for value in [
            shortening_m,
            duration_s,
            elapsed_s,
            strain,
            stress,
            force_n,
            tangent_n_m,
            increment,
            work,
        ] {
            h.update(&value.to_bits().to_le_bytes());
        }
        let response = CompressionReferenceResponse {
            elapsed_s,
            shortening_m,
            compressive_strain: strain,
            nominal_stress_pa: stress,
            force_n,
            tangent_n_m,
            model_loading_work_j: work,
            incremental_model_loading_work_j: increment,
            identity: h.finalize(),
        };
        self.strain = strain;
        self.elapsed_s = elapsed_s;
        self.work_j = work;
        self.identity = response.identity;
        self.last = Some(response.clone());
        Ok(response)
    }
}
