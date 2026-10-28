//! Source-bound, small incremental uniaxial Prony reference specimens.
//!
//! Coefficients are explicit scalar claims, never inferred from hardness,
//! density, a static compression curve, or a material name. The source query
//! and its original receipts remain attached to the reference. Authored time,
//! strain and rate limits bound model use; they do not establish experimental
//! validation, particularly when a spectrum was fitted in the frequency domain.
//!
//! Strain and force controls share one branch history. Positive compression
//! uses original length, area and volume. The initially relaxed incremental
//! state is a deviation about the source operating point, not a reconstruction
//! of unknown preload stress or energy. No densification, permanent set, damage,
//! lateral contraction, temperature evolution or nonlinear foam law is inferred.

use super::stress_control::StressControlledMaxwell;
use super::{GeneralizedMaxwell, PronyState, ViscoError};
use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::{
    InterpolationPolicy, MatDbError, MaterialAnswer, MaterialCard, PropertyKey, PropertyValue,
    QueryPoint, SelectionPolicy,
};
use fs_qty::{Density, Dims, Pressure, QuantitySpec};

const TIME: Dims = Dims([0, 0, 1, 0, 0, 0]);
const LENGTH: Dims = Dims([1, 0, 0, 0, 0, 0]);
const AREA: Dims = Dims([2, 0, 0, 0, 0, 0]);

/// Explicit model-use limits, not inferred source-validity claims.
#[derive(Debug, Clone, PartialEq)]
pub struct PronyReferenceLimits {
    pub minimum_incremental_strain: f64,
    pub maximum_incremental_strain: f64,
    pub maximum_absolute_strain_rate_s: f64,
    pub maximum_elapsed_s: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct PronyReferenceRequest {
    pub limits: PronyReferenceLimits,
    /// Bind the specimen's original length to this selected source coordinate.
    pub original_length_axis: Option<String>,
    /// Bind original area, not diameter or current compressed area.
    pub original_area_axis: Option<String>,
}

#[derive(Debug, Clone, PartialEq)]
pub enum PronyReferenceError {
    Source(MatDbError),
    Law(ViscoError),
    Invalid(&'static str),
}
impl core::fmt::Display for PronyReferenceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Source(e) => e.fmt(f),
            Self::Law(e) => e.fmt(f),
            Self::Invalid(message) => write!(f, "Prony reference: {message}"),
        }
    }
}
impl std::error::Error for PronyReferenceError {}
impl From<MatDbError> for PronyReferenceError {
    fn from(e: MatDbError) -> Self {
        Self::Source(e)
    }
}
impl From<ViscoError> for PronyReferenceError {
    fn from(e: ViscoError) -> Self {
        Self::Law(e)
    }
}
type Result<T> = core::result::Result<T, PronyReferenceError>;
fn require(ok: bool, message: &'static str) -> Result<()> {
    if ok {
        Ok(())
    } else {
        Err(PronyReferenceError::Invalid(message))
    }
}
fn positive(x: f64, message: &'static str) -> Result<f64> {
    require(x.is_finite() && x > 0.0, message)?;
    Ok(x)
}
fn parameter(
    card: &MaterialCard,
    point: &QueryPoint,
    name: &str,
    dims: Dims,
    selection: SelectionPolicy,
    properties: &mut Vec<(String, MaterialAnswer)>,
) -> Result<f64> {
    let answer = card
        .claims()
        .query_typed(&PropertyKey::new(name, dims), point, selection)?;
    card.claims().verify_receipt(&answer.receipt)?;
    let claim =
        card.claims()
            .claim(answer.receipt.selected)
            .ok_or(PronyReferenceError::Invalid(
                "selected coefficient claim is missing",
            ))?;
    require(
        matches!(claim.value, PropertyValue::Scalar { .. })
            && claim.interpolation == InterpolationPolicy::ConstantWithinValidity,
        "each Prony coefficient must be a validity-wide scalar; sampled curves are unsupported",
    )?;
    let value = answer.evidence.value.value;
    require(value.is_finite(), "nonfinite source coefficient")?;
    properties.push((name.to_owned(), answer));
    Ok(value)
}
fn hash_number(h: &mut DomainHasher, x: f64) {
    h.update(&x.to_bits().to_le_bytes());
}
fn hash_axis(h: &mut DomainHasher, axis: &Option<String>) {
    match axis {
        None => {
            h.update(&[0]);
        }
        Some(name) => {
            h.update(&[1]);
            h.update(&(name.len() as u64).to_le_bytes());
            h.update(name.as_bytes());
        }
    }
}

/// Immutable coefficients plus their exact original source answers.
#[derive(Debug, Clone)]
pub struct ResolvedPronyReference {
    request: PronyReferenceRequest,
    card_identity: ContentHash,
    query: QueryPoint,
    properties: Vec<(String, MaterialAnswer)>,
    stress_law: StressControlledMaxwell,
    density_kg_m3: Option<f64>,
    original_length_m: Option<f64>,
    original_area_m2: Option<f64>,
    identity: ContentHash,
}
impl ResolvedPronyReference {
    /// Resolve `prony-branch-count`, `prony-equilibrium-modulus`, and indexed
    /// `prony-modulus-N`/`prony-relaxation-time-N` claims in Pa and seconds.
    /// Density is optional for this quasistatic reference; absence is retained.
    /// At most 64 strictly positive arms are admitted; zero arms are elastic.
    pub fn resolve(
        card: &MaterialCard,
        point: &QueryPoint,
        request: &PronyReferenceRequest,
        selection: SelectionPolicy,
    ) -> Result<Self> {
        let limits = &request.limits;
        require(
            limits.minimum_incremental_strain.is_finite()
                && limits.maximum_incremental_strain.is_finite()
                && limits.minimum_incremental_strain <= 0.0
                && limits.maximum_incremental_strain >= 0.0
                && limits.minimum_incremental_strain < limits.maximum_incremental_strain,
            "finite ordered strain limits must contain the initial zero increment",
        )?;
        positive(
            limits.maximum_absolute_strain_rate_s,
            "positive finite maximum strain rate required",
        )?;
        positive(
            limits.maximum_elapsed_s,
            "positive finite maximum elapsed time required",
        )?;
        let mut properties = Vec::new();
        let count = parameter(
            card,
            point,
            "prony-branch-count",
            Dims::NONE,
            selection,
            &mut properties,
        )?;
        require(
            (0.0..=64.0).contains(&count) && count.fract() == 0.0,
            "Prony branch count must be an integer in 0..64",
        )?;
        let e_inf = positive(
            parameter(
                card,
                point,
                "prony-equilibrium-modulus",
                Pressure::DIMS,
                selection,
                &mut properties,
            )?,
            "positive finite equilibrium modulus required",
        )?;
        let mut terms = Vec::with_capacity(count as usize);
        for branch in 1..=count as usize {
            let e = positive(
                parameter(
                    card,
                    point,
                    &format!("prony-modulus-{branch}"),
                    Pressure::DIMS,
                    selection,
                    &mut properties,
                )?,
                "positive finite branch modulus required",
            )?;
            let tau = positive(
                parameter(
                    card,
                    point,
                    &format!("prony-relaxation-time-{branch}"),
                    TIME,
                    selection,
                    &mut properties,
                )?,
                "positive finite relaxation time required",
            )?;
            terms.push((e, tau));
        }
        let density_kg_m3 = if card.claims().claims_for("density").is_empty() {
            None
        } else {
            Some(positive(
                parameter(
                    card,
                    point,
                    "density",
                    Density::DIMS,
                    selection,
                    &mut properties,
                )?,
                "positive finite source density required",
            )?)
        };
        // Unused coordinates must not masquerade as constraints on this model.
        for (axis, _) in point.axes() {
            require(
                properties.iter().any(|(_, answer)| {
                    card.claims()
                        .claim(answer.receipt.selected)
                        .is_some_and(|claim| {
                            claim.validity.bounds().contains_key(axis)
                                && claim.validity.axis_quantities().get(axis)
                                    == point.axis_quantities().get(axis)
                        })
                }),
                "query coordinate is not owned by the selected coefficient claims",
            )?;
        }
        let geometry = |axis: &Option<String>, dims: Dims| -> Result<Option<f64>> {
            let Some(axis) = axis else { return Ok(None) };
            require(
                !axis.is_empty()
                    && point.axis_quantities().get(axis) == Some(&QuantitySpec::dimensional(dims)),
                "source geometry needs a selected, explicitly dimensioned coordinate",
            )?;
            Ok(Some(positive(
                *point.axes().get(axis).ok_or(PronyReferenceError::Invalid(
                    "missing original geometry coordinate",
                ))?,
                "positive finite original source geometry required",
            )?))
        };
        let original_length_m = geometry(&request.original_length_axis, LENGTH)?;
        let original_area_m2 = geometry(&request.original_area_axis, AREA)?;
        let stress_law = StressControlledMaxwell::new(&GeneralizedMaxwell::new(e_inf, terms)?)?;
        let mut h = DomainHasher::new("org.frankensim.fs-material.prony-reference.v1");
        h.update(card.content_hash().as_bytes());
        for (_, answer) in &properties {
            h.update(answer.receipt.content_hash().as_bytes());
        }
        for x in [
            limits.minimum_incremental_strain,
            limits.maximum_incremental_strain,
            limits.maximum_absolute_strain_rate_s,
            limits.maximum_elapsed_s,
        ] {
            hash_number(&mut h, x);
        }
        hash_axis(&mut h, &request.original_length_axis);
        hash_axis(&mut h, &request.original_area_axis);
        Ok(Self {
            request: request.clone(),
            card_identity: card.content_hash(),
            query: point.clone(),
            properties,
            stress_law,
            density_kg_m3,
            original_length_m,
            original_area_m2,
            identity: h.finalize(),
        })
    }
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
    pub const fn card_identity(&self) -> ContentHash {
        self.card_identity
    }
    pub const fn query(&self) -> &QueryPoint {
        &self.query
    }
    pub fn properties(&self) -> &[(String, MaterialAnswer)] {
        &self.properties
    }
    pub fn equilibrium_modulus_pa(&self) -> f64 {
        self.stress_law.law().e_inf
    }
    pub fn branches(&self) -> &[(f64, f64)] {
        &self.stress_law.law().terms
    }
    pub const fn request(&self) -> &PronyReferenceRequest {
        &self.request
    }

    /// Start a relaxed incremental specimen using fixed original geometry.
    pub fn coupon(&self, length_m: f64, area_m2: f64) -> Result<PronyReferenceCoupon> {
        positive(
            length_m,
            "positive finite original specimen length required",
        )?;
        positive(area_m2, "positive finite original specimen area required")?;
        for (actual, source) in [
            (length_m, self.original_length_m),
            (area_m2, self.original_area_m2),
        ] {
            if let Some(source) = source {
                require(
                    (actual - source).abs() <= 8.0 * f64::EPSILON * source.abs(),
                    "specimen geometry differs from the bound original source geometry",
                )?;
            }
        }
        let volume_m3 = positive(length_m * area_m2, "unrepresentable original volume")?;
        let mass_kg = self
            .density_kg_m3
            .map(|rho| positive(rho * volume_m3, "unrepresentable source mass"))
            .transpose()?;
        let mut h = DomainHasher::new("org.frankensim.fs-material.prony-specimen.v1");
        h.update(self.identity.as_bytes());
        hash_number(&mut h, length_m);
        hash_number(&mut h, area_m2);
        Ok(PronyReferenceCoupon {
            reference: self.clone(),
            state: self.stress_law.law().state(),
            length_m,
            area_m2,
            volume_m3,
            mass_kg,
            elapsed_s: 0.0,
            last: None,
            identity: h.finalize(),
        })
    }
}

/// Every mechanical quantity is incremental about the fixed operating point.
/// Work may decrease during unloading; integrated physical dissipation may not.
#[derive(Debug, Clone, PartialEq)]
pub struct PronyReferenceResponse {
    pub elapsed_s: f64,
    pub incremental_strain: f64,
    pub incremental_displacement_m: f64,
    pub incremental_stress_pa: f64,
    pub incremental_force_n: f64,
    pub stored_energy_j: f64,
    pub model_work_j: f64,
    pub dissipated_energy_j: f64,
    pub incremental_model_work_j: f64,
    pub incremental_dissipated_energy_j: f64,
    pub energy_balance_residual_j: f64,
    pub identity: ContentHash,
}

#[derive(Debug, Clone)]
pub struct PronyReferenceCoupon {
    reference: ResolvedPronyReference,
    state: PronyState,
    length_m: f64,
    area_m2: f64,
    volume_m3: f64,
    mass_kg: Option<f64>,
    elapsed_s: f64,
    last: Option<PronyReferenceResponse>,
    identity: ContentHash,
}
impl PronyReferenceCoupon {
    pub const fn mass_kg(&self) -> Option<f64> {
        self.mass_kg
    }
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
    pub const fn last_response(&self) -> Option<&PronyReferenceResponse> {
        self.last.as_ref()
    }

    fn admit_time(&self, dt: f64) -> Result<f64> {
        positive(
            dt,
            "positive finite step duration required; instantaneous jumps are unsupported",
        )?;
        let elapsed = self.elapsed_s + dt;
        require(
            elapsed.is_finite()
                && elapsed > self.elapsed_s
                && elapsed <= self.reference.request.limits.maximum_elapsed_s,
            "step exceeds authored time limit or cannot advance the time coordinate",
        )?;
        Ok(elapsed)
    }
    fn admit_path(&self, strain: [f64; 2], rate: [f64; 2]) -> Result<()> {
        let limits = &self.reference.request.limits;
        require(
            strain.iter().chain(rate.iter()).all(|x| x.is_finite()),
            "unrepresentable history bounds",
        )?;
        require(
            strain[0] >= limits.minimum_incremental_strain
                && strain[1] <= limits.maximum_incremental_strain,
            "whole-ramp strain bounds exceed authored limits; shorten the step or revise the model scope",
        )?;
        require(
            rate[0] >= -limits.maximum_absolute_strain_rate_s
                && rate[1] <= limits.maximum_absolute_strain_rate_s,
            "whole-ramp strain-rate bounds exceed the authored limit",
        )?;
        Ok(())
    }
    /// Exact linear strain ramp. Equal endpoint strain performs stress relaxation.
    pub fn advance_strain(
        &mut self,
        strain: f64,
        duration_s: f64,
    ) -> Result<PronyReferenceResponse> {
        let elapsed = self.admit_time(duration_s)?;
        require(strain.is_finite(), "finite target strain required")?;
        let rate = (strain - self.state.strain) / duration_s;
        self.admit_path(
            [strain.min(self.state.strain), strain.max(self.state.strain)],
            [rate, rate],
        )?;
        let mut next = self.state.clone();
        // Integrate each increment directly; do not subtract large cumulative
        // ledgers to recover a potentially tiny new work or heat contribution.
        next.work = 0.0;
        next.dissipated = 0.0;
        self.reference
            .stress_law
            .law()
            .step(&mut next, strain, duration_s);
        let (work, heat) = (next.work, next.dissipated);
        next.work = self.state.work + work;
        next.dissipated = self.state.dissipated + heat;
        self.accept(next, work, heat, elapsed, 0, strain, duration_s)
    }
    /// Exact linear nominal-force ramp. Constant force creeps; zero force recovers.
    /// Conservative whole-step modal bounds can refuse a large step even if its
    /// endpoint is in range. Subdivide that history to obtain tighter bounds.
    pub fn advance_force(
        &mut self,
        force_n: f64,
        duration_s: f64,
    ) -> Result<PronyReferenceResponse> {
        let elapsed = self.admit_time(duration_s)?;
        require(
            force_n.is_finite(),
            "finite target incremental force required",
        )?;
        let mut next = self.state.clone();
        let response =
            self.reference
                .stress_law
                .step(&mut next, force_n / self.area_m2, duration_s)?;
        self.admit_path(response.strain_bounds, response.strain_rate_bounds_s)?;
        self.accept(
            next,
            response.incremental_work_j_m3,
            response.incremental_dissipated_j_m3,
            elapsed,
            1,
            force_n,
            duration_s,
        )
    }
    fn accept(
        &mut self,
        next: PronyState,
        work: f64,
        heat: f64,
        elapsed: f64,
        control: u8,
        target: f64,
        dt: f64,
    ) -> Result<PronyReferenceResponse> {
        let law = self.reference.stress_law.law();
        let mut storage = 0.5 * (law.e_inf * next.strain) * next.strain;
        for (&(e, _), &q) in law.terms.iter().zip(&next.q) {
            storage += 0.5 * (e * q) * q;
        }
        require(
            next.q.iter().all(|q| q.is_finite())
                && heat.is_finite()
                && heat >= 0.0
                && storage.is_finite()
                && storage >= 0.0
                && next.dissipated >= self.state.dissipated,
            "nonfinite state, storage or physical dissipation",
        )?;
        let values = [
            elapsed,
            next.strain,
            next.strain * self.length_m,
            next.stress,
            next.stress * self.area_m2,
            storage * self.volume_m3,
            next.work * self.volume_m3,
            next.dissipated * self.volume_m3,
            work * self.volume_m3,
            heat * self.volume_m3,
            (next.work - storage - next.dissipated) * self.volume_m3,
        ];
        require(
            values.iter().all(|v| v.is_finite()),
            "unrepresentable specimen response",
        )?;
        let mut h = DomainHasher::new("org.frankensim.fs-material.prony-history-step.v1");
        h.update(self.identity.as_bytes());
        h.update(&[control]);
        hash_number(&mut h, target);
        hash_number(&mut h, dt);
        for value in values {
            hash_number(&mut h, value);
        }
        let response = PronyReferenceResponse {
            elapsed_s: values[0],
            incremental_strain: values[1],
            incremental_displacement_m: values[2],
            incremental_stress_pa: values[3],
            incremental_force_n: values[4],
            stored_energy_j: values[5],
            model_work_j: values[6],
            dissipated_energy_j: values[7],
            incremental_model_work_j: values[8],
            incremental_dissipated_energy_j: values[9],
            energy_balance_residual_j: values[10],
            identity: h.finalize(),
        };
        self.state = next;
        self.elapsed_s = elapsed;
        self.identity = response.identity;
        self.last = Some(response.clone());
        Ok(response)
    }
}
