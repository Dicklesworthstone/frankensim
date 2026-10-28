//! Source-bound axial comparisons from reported tensile or compression moduli.
//!
//! This is an explicitly authored, zero-origin linear MODULUS REFERENCE.
//! A reported compression modulus may be fitted over a finite strain interval;
//! it is not promoted to a measured initial tangent or a complete stress curve.
//! The caller's strain ceiling is a model restriction, not source calibration.
//!
//! The existing generalized-Maxwell owner's equilibrium-spring limit supplies
//! stress and ramp work. Its use here does not assert a measured equilibrium
//! modulus, absent real dissipation, recoverable energy, isotropy, or Poisson
//! ratio. Only strictly monotonic loading in the selected direction is exposed:
//! unloading, holds/relaxation, cyclic response and failure are unsupported.
//! Work is labelled model loading work, never recovered material energy.

use crate::state_point::{
    MaterialPropertySelection, MaterialStatePointError, ResolvedMaterialStatePoint,
    ScalarPropertyRequirement, is_fixed_external_coordinate, resolve_material_state_point,
};
use crate::visco::{GeneralizedMaxwell, PronyState};
use fs_blake3::{ContentHash, DomainHasher};
use fs_matdb::{MaterialCard, QueryPoint};
use fs_qty::{Density, Dims, Pressure, QuantitySpec};

const SPEED: Dims = Dims([1, 0, -1, 0, 0, 0]);

/// One reported axial test mode; tension and compression are not interchangeable.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AxialLoadingMode {
    Tension,
    Compression,
}
impl AxialLoadingMode {
    fn sign(self) -> f64 {
        match self {
            Self::Tension => 1.0,
            Self::Compression => -1.0,
        }
    }
    fn accepts_name(self, name: &str) -> bool {
        // Explicit property-vocabulary admission, not material-name inference.
        // Unknown aliases require a deliberate source-manifest role mapping.
        match self {
            Self::Tension => matches!(name, "tensile-modulus" | "tensile_modulus"),
            Self::Compression => matches!(
                name,
                "compression-modulus"
                    | "compression_modulus"
                    | "compressive-modulus"
                    | "compressive_modulus"
            ),
        }
    }
}

/// Complete caller-selected source roles and authored model restrictions.
#[derive(Clone, Debug, PartialEq)]
pub struct AxialReferenceRequest {
    pub density: ScalarPropertyRequirement,
    pub modulus: ScalarPropertyRequirement,
    pub mode: AxialLoadingMode,
    /// Positive magnitude below one; not a measured linearity/failure bound.
    pub maximum_absolute_strain: f64,
    /// Exact source crosshead-speed axis, when the modulus declares one.
    /// Its typed m/s coordinate is enforced against displacement/time steps.
    /// None declares no source-supplied crosshead-speed constraint.
    pub crosshead_speed_axis: Option<String>,
}

#[derive(Clone, Debug, PartialEq)]
pub enum AxialReferenceError {
    Source(MaterialStatePointError),
    Invalid(&'static str),
}
impl core::fmt::Display for AxialReferenceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Source(error) => error.fmt(f),
            Self::Invalid(what) => write!(f, "axial modulus reference: {what}"),
        }
    }
}
impl std::error::Error for AxialReferenceError {}
impl From<MaterialStatePointError> for AxialReferenceError {
    fn from(error: MaterialStatePointError) -> Self {
        Self::Source(error)
    }
}
type Result<T> = core::result::Result<T, AxialReferenceError>;
fn require(condition: bool, what: &'static str) -> Result<()> {
    if condition {
        Ok(())
    } else {
        Err(AxialReferenceError::Invalid(what))
    }
}
fn positive(value: f64, what: &'static str) -> Result<f64> {
    require(value.is_finite() && value > 0.0, what)?;
    Ok(value)
}
fn scalar_requirement(requirement: &ScalarPropertyRequirement, dims: Dims) -> bool {
    requirement.quantity() == QuantitySpec::dimensional(dims)
        && requirement.hardness_test().is_none()
        && requirement.elastic_component().is_none()
        && requirement.strain_component().is_none()
        && requirement.stress_component().is_none()
}

/// Frozen reference retaining the original property bundle and its uncertainty.
/// No synthetic modulus, Poisson, loss or strength claim is inserted in the card.
#[derive(Clone, Debug)]
pub struct ResolvedAxialModulusReference {
    source: ResolvedMaterialStatePoint,
    request: AxialReferenceRequest,
    density_kg_m3: f64,
    modulus_pa: f64,
    crosshead_speed_m_s: Option<f64>,
    law: GeneralizedMaxwell,
    identity: ContentHash,
}
impl ResolvedAxialModulusReference {
    /// Resolve density and the explicitly named tensile/compression modulus at
    /// one complete source point. Selection remains the original matdb policy.
    /// Flexural/strength/secant-point stress keys and wrong loading modes refuse.
    pub fn resolve(
        card: &MaterialCard,
        point: &QueryPoint,
        request: &AxialReferenceRequest,
        selection: MaterialPropertySelection,
    ) -> Result<Self> {
        require(
            scalar_requirement(&request.density, Density::DIMS)
                && scalar_requirement(&request.modulus, Pressure::DIMS),
            "density and modulus require complete scalar dimension-only keys",
        )?;
        require(
            request.mode.accepts_name(request.modulus.name()),
            "select the matching reported tensile/compression modulus; no flexural, strength or point-stress substitution",
        )?;
        require(
            request.maximum_absolute_strain.is_finite()
                && request.maximum_absolute_strain > 0.0
                && request.maximum_absolute_strain < 1.0,
            "authored strain ceiling must be finite and in (0,1)",
        )?;
        let source = resolve_material_state_point(
            card,
            point,
            &[request.density.clone(), request.modulus.clone()],
            selection,
        )?;
        for property in source.properties() {
            for (axis, _) in property.source_axes() {
                require(
                    is_fixed_external_coordinate(
                        axis,
                        point.axis_quantities().get(axis).copied(),
                        point.axes().get(axis).copied(),
                    ),
                    "a sampled strain/rate/frequency/history law cannot become a fixed axial modulus",
                )?;
            }
        }
        let density = source
            .property(request.density.name())
            .ok_or(AxialReferenceError::Invalid("resolved density missing"))?;
        let modulus =
            source
                .property(request.modulus.name())
                .ok_or(AxialReferenceError::Invalid(
                    "resolved axial modulus missing",
                ))?;
        let rho = positive(density.value_si(), "positive finite density required")?;
        let value = positive(
            modulus.value_si(),
            "positive finite reported modulus required",
        )?;
        let claim = card
            .claims()
            .claim(modulus.answer().receipt.selected)
            .ok_or(AxialReferenceError::Invalid(
                "selected modulus claim missing",
            ))?;
        let speed = if let Some(axis) = &request.crosshead_speed_axis {
            require(
                claim.validity.axis_quantities().get(axis)
                    == Some(&QuantitySpec::dimensional(SPEED))
                    && claim.validity.bounds().contains_key(axis),
                "crosshead-speed owner must name the selected modulus's typed m/s axis",
            )?;
            Some(positive(
                *point.axes().get(axis).ok_or(AxialReferenceError::Invalid(
                    "missing source crosshead speed",
                ))?,
                "source crosshead speed must be positive and finite",
            )?)
        } else {
            None
        };
        for (axis, quantity) in claim.validity.axis_quantities() {
            if quantity.dims() == SPEED {
                require(
                    request.crosshead_speed_axis.as_deref() == Some(axis.as_str()),
                    "source crosshead-speed applicability needs its explicit displacement/time owner",
                )?;
            }
        }
        let law = GeneralizedMaxwell::new(value, Vec::new()).map_err(|_| {
            AxialReferenceError::Invalid("native equilibrium spring refused modulus")
        })?;
        let mut h = DomainHasher::new("org.frankensim.fs-material.axial-modulus-reference.v1");
        h.update(source.identity().as_bytes());
        h.update(&[match request.mode {
            AxialLoadingMode::Tension => 0,
            AxialLoadingMode::Compression => 1,
        }]);
        h.update(&request.maximum_absolute_strain.to_bits().to_le_bytes());
        if let Some(axis) = &request.crosshead_speed_axis {
            h.update(&[1]);
            h.update(&(axis.len() as u64).to_le_bytes());
            h.update(axis.as_bytes());
        } else {
            h.update(&[0]);
        }
        Ok(Self {
            source,
            request: request.clone(),
            density_kg_m3: rho,
            modulus_pa: value,
            crosshead_speed_m_s: speed,
            law,
            identity: h.finalize(),
        })
    }

    pub const fn source(&self) -> &ResolvedMaterialStatePoint {
        &self.source
    }
    pub const fn request(&self) -> &AxialReferenceRequest {
        &self.request
    }
    pub const fn reported_modulus_pa(&self) -> f64 {
        self.modulus_pa
    }
    pub const fn density_kg_m3(&self) -> f64 {
        self.density_kg_m3
    }
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Begin a fresh reference specimen of uniform original area and length.
    /// The source test's geometry is not inferred from the requested specimen.
    /// Inertia uses original rho*A*L; no lateral contraction or updated density
    /// is invented. The reduced stiffness is exactly E*A/L.
    pub fn coupon(&self, length_m: f64, area_m2: f64) -> Result<AxialReferenceCoupon> {
        positive(length_m, "positive finite reference length required")?;
        positive(area_m2, "positive finite original area required")?;
        let volume = positive(length_m * area_m2, "reference volume is unrepresentable")?;
        let mass = positive(
            self.density_kg_m3 * volume,
            "reference mass is unrepresentable",
        )?;
        let stiffness = positive(
            self.modulus_pa * area_m2 / length_m,
            "axial stiffness is unrepresentable",
        )?;
        let mut h = DomainHasher::new("org.frankensim.fs-material.axial-reference-coupon.v1");
        h.update(self.identity.as_bytes());
        h.update(&length_m.to_bits().to_le_bytes());
        h.update(&area_m2.to_bits().to_le_bytes());
        Ok(AxialReferenceCoupon {
            material: self.clone(),
            length_m,
            area_m2,
            volume_m3: volume,
            mass_kg: mass,
            stiffness_n_m: stiffness,
            state: self.law.state(),
            elapsed_s: 0.0,
            displacement_m: 0.0,
            identity: h.finalize(),
            last: None,
        })
    }
}

/// Accepted nominal reference response. Work is the integral of this model's
/// loading force, not an assertion about recoverable material energy or heat.
#[derive(Clone, Debug, PartialEq)]
pub struct AxialReferenceResponse {
    pub elapsed_s: f64,
    pub displacement_m: f64,
    pub engineering_strain: f64,
    pub reference_stress_pa: f64,
    pub reference_force_n: f64,
    pub model_loading_work_j: f64,
    pub incremental_model_loading_work_j: f64,
    pub mass_kg: f64,
    pub identity: ContentHash,
}

/// Monotonic axial reference driven by displacement. Source state stays fixed.
#[derive(Clone, Debug)]
pub struct AxialReferenceCoupon {
    material: ResolvedAxialModulusReference,
    length_m: f64,
    area_m2: f64,
    volume_m3: f64,
    mass_kg: f64,
    stiffness_n_m: f64,
    state: PronyState,
    elapsed_s: f64,
    displacement_m: f64,
    identity: ContentHash,
    last: Option<AxialReferenceResponse>,
}
impl AxialReferenceCoupon {
    pub const fn material(&self) -> &ResolvedAxialModulusReference {
        &self.material
    }
    pub const fn stiffness_n_m(&self) -> f64 {
        self.stiffness_n_m
    }
    pub const fn mass_kg(&self) -> f64 {
        self.mass_kg
    }
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
    pub fn last_response(&self) -> Option<&AxialReferenceResponse> {
        self.last.as_ref()
    }

    /// Advance one strictly monotonic linear displacement ramp. Holds,
    /// unloading and reversal refuse instead of manufacturing a relaxation or
    /// recovery law. Known source crosshead speed must match the imposed ramp.
    /// All validation and representability checks finish before state changes.
    pub fn advance_displacement(
        &mut self,
        displacement_m: f64,
        duration_s: f64,
    ) -> Result<AxialReferenceResponse> {
        positive(duration_s, "finite positive loading duration required")?;
        require(displacement_m.is_finite(), "finite displacement required")?;
        let strain = displacement_m / self.length_m;
        let signed_increment =
            self.material.request.mode.sign() * (displacement_m - self.displacement_m);
        require(
            signed_increment.is_finite()
                && signed_increment > 0.0
                && strain.is_finite()
                && self.material.request.mode.sign() * strain > 0.0,
            "loading mode forbids a hold, unloading or reversal",
        )?;
        require(
            strain.abs() <= self.material.request.maximum_absolute_strain,
            "authored reference strain ceiling exceeded",
        )?;
        if let Some(required_speed) = self.material.crosshead_speed_m_s {
            let imposed_speed = signed_increment / duration_s;
            // Arithmetic tolerance only, not a widened source rate domain.
            require(
                imposed_speed.is_finite()
                    && (imposed_speed / required_speed - 1.0).abs() <= 1.0e-10,
                "imposed crosshead speed differs from the source protocol",
            )?;
        }
        let elapsed = positive(self.elapsed_s + duration_s, "elapsed duration overflow")?;
        require(
            elapsed > self.elapsed_s,
            "elapsed duration increment is unrepresentable",
        )?;
        let mut candidate = self.state.clone();
        let stress = self.material.law.step(&mut candidate, strain, duration_s);
        let force = stress * self.area_m2;
        let work = candidate.work * self.volume_m3;
        let increment = (candidate.work - self.state.work) * self.volume_m3;
        require(
            stress.is_finite()
                && force.is_finite()
                && self.material.request.mode.sign() * force > 0.0,
            "reference stress or force is unrepresentable",
        )?;
        positive(work, "model loading work is unrepresentable")?;
        positive(
            increment,
            "incremental model loading work is unrepresentable",
        )?;
        let mut h = DomainHasher::new("org.frankensim.fs-material.axial-reference-step.v1");
        h.update(self.identity.as_bytes());
        for x in [
            displacement_m,
            duration_s,
            elapsed,
            strain,
            stress,
            force,
            work,
            increment,
        ] {
            h.update(&x.to_bits().to_le_bytes());
        }
        let response = AxialReferenceResponse {
            elapsed_s: elapsed,
            displacement_m,
            engineering_strain: strain,
            reference_stress_pa: stress,
            reference_force_n: force,
            model_loading_work_j: work,
            incremental_model_loading_work_j: increment,
            mass_kg: self.mass_kg,
            identity: h.finalize(),
        };
        self.state = candidate;
        self.elapsed_s = elapsed;
        self.displacement_m = displacement_m;
        self.identity = response.identity;
        self.last = Some(response.clone());
        Ok(response)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::state_point::ScalarAdmissibility;
    use fs_evidence::ValidityDomain;
    use fs_matdb::{
        ClaimSet, InterpolationPolicy, MaterialStateId, PropertyClaim, PropertyKey, PropertyValue,
        Provenance, UncertaintyModel,
    };

    fn fixture(curved_tension: bool) -> MaterialCard {
        let mut claims = ClaimSet::new();
        for (name, value, dims) in [
            ("density", 1000.0, Density::DIMS),
            ("tensile-modulus", 2.0e9, Pressure::DIMS),
            ("compression-modulus", 1.0e9, Pressure::DIMS),
            ("flexural-modulus", 3.0e9, Pressure::DIMS),
            ("compressive-stress", 20.0e6, Pressure::DIMS),
        ] {
            let mut validity = ValidityDomain::unconstrained().with_quantity(
                "source-grade",
                QuantitySpec::dimensional(Dims::NONE),
                1.0,
                1.0,
            );
            if name != "density" {
                validity = validity.with_quantity(
                    "test-speed",
                    QuantitySpec::dimensional(SPEED),
                    0.001,
                    0.001,
                );
            }
            let curved = curved_tension && name == "tensile-modulus";
            if curved {
                validity = validity.with_quantity(
                    "strain",
                    QuantitySpec::dimensional(Dims::NONE),
                    0.0,
                    0.01,
                );
            }
            claims
                .insert_claim(PropertyClaim {
                    key: PropertyKey::new(name, dims),
                    value: if curved {
                        PropertyValue::Curve {
                            abscissa: "strain".into(),
                            abscissa_dims: Dims::NONE,
                            knots: vec![(0.0, value), (0.01, value * 0.8)],
                            dims,
                        }
                    } else {
                        PropertyValue::Scalar { value, dims }
                    },
                    validity,
                    uncertainty: UncertaintyModel::Unstated,
                    interpolation: if curved {
                        InterpolationPolicy::LinearInside
                    } else {
                        InterpolationPolicy::ConstantWithinValidity
                    },
                    observations: vec![],
                    provenance: Provenance {
                        source: "synthetic axial-reference numerical fixture".into(),
                        license: "test-only".into(),
                        artifact: None,
                    },
                })
                .unwrap();
        }
        MaterialCard::assemble(
            MaterialStateId {
                chemistry: "synthetic axial reference".into(),
                phase: "solid".into(),
                process: "test fixture".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap()
    }
    fn point(speed: f64) -> QueryPoint {
        QueryPoint::new()
            .with_quantity("source-grade", QuantitySpec::dimensional(Dims::NONE), 1.0)
            .unwrap()
            .with_quantity("test-speed", QuantitySpec::dimensional(SPEED), speed)
            .unwrap()
    }
    fn request(mode: AxialLoadingMode) -> AxialReferenceRequest {
        AxialReferenceRequest {
            density: ScalarPropertyRequirement::try_new(
                "density",
                Density::DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )
            .unwrap(),
            modulus: ScalarPropertyRequirement::try_new(
                match mode {
                    AxialLoadingMode::Tension => "tensile-modulus",
                    AxialLoadingMode::Compression => "compression-modulus",
                },
                Pressure::DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )
            .unwrap(),
            mode,
            maximum_absolute_strain: 0.002,
            crosshead_speed_axis: Some("test-speed".into()),
        }
    }
    fn resolve(mode: AxialLoadingMode) -> ResolvedAxialModulusReference {
        ResolvedAxialModulusReference::resolve(
            &fixture(false),
            &point(0.001),
            &request(mode),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap()
    }

    #[test]
    fn g1_distinct_axial_modes_have_independent_ea_over_l_force_mass_and_work() {
        for mode in [AxialLoadingMode::Tension, AxialLoadingMode::Compression] {
            let material = resolve(mode);
            let mut coupon = material.coupon(0.1, 1.0e-4).unwrap();
            let e = if mode == AxialLoadingMode::Tension {
                2.0e9
            } else {
                1.0e9
            };
            assert_eq!(coupon.stiffness_n_m(), e * 1.0e-4 / 0.1);
            let first = coupon
                .advance_displacement(mode.sign() * 0.0001, 0.1)
                .unwrap();
            let force = mode.sign() * e * 1.0e-4 * 0.001;
            assert!((first.reference_force_n / force - 1.0).abs() < 1.0e-14);
            assert!((first.mass_kg / 0.01 - 1.0).abs() < 1.0e-14);
            assert!(
                (first.model_loading_work_j / (0.5 * force * first.displacement_m) - 1.0).abs()
                    < 1.0e-14
            );
            let second = coupon
                .advance_displacement(mode.sign() * 0.0002, 0.1)
                .unwrap();
            assert!((second.reference_force_n / (2.0 * force) - 1.0).abs() < 1.0e-14);
            assert!(
                (second.model_loading_work_j / first.model_loading_work_j - 4.0).abs() < 1.0e-14
            );
            let mut single = material.coupon(0.1, 1.0e-4).unwrap();
            let one = single
                .advance_displacement(mode.sign() * 0.0002, 0.2)
                .unwrap();
            assert!((one.model_loading_work_j / second.model_loading_work_j - 1.0).abs() < 1.0e-14);
            assert_ne!(one.identity, second.identity);
            assert_eq!(material.source().properties().len(), 2);
            assert!(
                material.source().properties().iter().all(|p| p
                    .answer()
                    .evidence
                    .value
                    .uncertainty
                    == UncertaintyModel::Unstated)
            );
        }
        assert_ne!(
            resolve(AxialLoadingMode::Tension).identity(),
            resolve(AxialLoadingMode::Compression).identity()
        );
    }

    #[test]
    fn g0_source_roles_rate_and_nonlinear_tables_refuse() {
        let card = fixture(false);
        let mut req = request(AxialLoadingMode::Tension);
        assert!(
            ResolvedAxialModulusReference::resolve(
                &card,
                &point(0.002),
                &req,
                MaterialPropertySelection::SingleClaimOnly
            )
            .is_err()
        );
        assert!(
            ResolvedAxialModulusReference::resolve(
                &card,
                &QueryPoint::new(),
                &req,
                MaterialPropertySelection::SingleClaimOnly
            )
            .is_err()
        );
        req.crosshead_speed_axis = None;
        assert!(
            ResolvedAxialModulusReference::resolve(
                &card,
                &point(0.001),
                &req,
                MaterialPropertySelection::SingleClaimOnly
            )
            .is_err()
        );
        for name in [
            "compression-modulus",
            "flexural-modulus",
            "compressive-stress",
        ] {
            let mut req = request(AxialLoadingMode::Tension);
            req.modulus = ScalarPropertyRequirement::try_new(
                name,
                Pressure::DIMS,
                ScalarAdmissibility::StrictlyPositive,
            )
            .unwrap();
            assert!(
                ResolvedAxialModulusReference::resolve(
                    &card,
                    &point(0.001),
                    &req,
                    MaterialPropertySelection::SingleClaimOnly
                )
                .is_err()
            );
        }
        let p = point(0.001)
            .with_quantity("strain", QuantitySpec::dimensional(Dims::NONE), 0.001)
            .unwrap();
        assert!(
            ResolvedAxialModulusReference::resolve(
                &fixture(true),
                &p,
                &request(AxialLoadingMode::Tension),
                MaterialPropertySelection::SingleClaimOnly
            )
            .is_err()
        );
    }

    #[test]
    fn g0_reversal_hold_overstrain_and_protocol_errors_leave_history_unchanged() {
        let material = resolve(AxialLoadingMode::Tension);
        assert!(material.coupon(0.0, 1.0e-4).is_err());
        assert!(material.coupon(0.1, f64::MAX).is_err());
        let mut coupon = material.coupon(0.1, 1.0e-4).unwrap();
        let accepted = coupon.advance_displacement(0.0001, 0.1).unwrap();
        for (displacement, duration) in [
            (-0.0001, 0.1),
            (0.0001, 0.1),
            (0.00005, 0.1),
            (0.0003, 0.2),
            (0.0002, 0.2),
            (0.0002, 0.0),
            (f64::NAN, 0.1),
        ] {
            assert!(coupon.advance_displacement(displacement, duration).is_err());
            assert_eq!(coupon.last_response(), Some(&accepted));
            assert_eq!(coupon.identity(), accepted.identity);
        }
        let next = coupon.advance_displacement(0.0002, 0.1).unwrap();
        assert!((next.model_loading_work_j / 0.04 - 1.0).abs() < 1.0e-14);
    }
}
