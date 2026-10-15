//! Immutable material/interface cards lowered to the shared point-contact law.
//!
//! The same admitted Hertz coefficient can feed a potential-based solver and
//! the original finite-patch evaluator. No material-name dispatch, fitted
//! restitution, or reimplementation of elliptic Hertz geometry is used.
use crate::interface_binding::{BoundNormalContactLaw, BoundNormalContactModel};
use super::{ApplicabilityInput, InputUncertainty, NormalPatchError, NormalPatchGeometry,
    NormalPatchIdentity, NormalPatchLaw, NormalPatchReceipt, NormalPatchRequest,
    PointNormalPatchReceipt};

/// Geometry of the RELATIVE local gap; not an inferred material property.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum PointContactGeometry {
    /// Axisymmetric quadratic gap of radius R [m].
    SpherePlane { effective_radius_m: f64 },
    /// Actual two-positive-principal-curvature gap [1/m], maximum first.
    EllipticParaboloid {
        maximum_curvature_m_inverse: f64,
        minimum_curvature_m_inverse: f64,
    },
}

/// Caller-selected operating box, checked against the original model card.
/// This does not extend its evidence or qualify the geometry approximation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PointContactEnvelope {
    /// Inclusive indentation ceiling [m], strictly positive.
    pub maximum_indentation_m: f64,
    /// Inclusive magnitude ceiling on closure/separation rate [m/s].
    pub maximum_absolute_rate_m_per_s: f64,
    /// Maximum requested sample duration [s], strictly positive.
    pub maximum_step_s: f64,
    /// Finite half-space depth from the geometry owner [m].
    pub half_space_depth_m: f64,
    /// Smallest relevant layer thickness from the geometry owner [m].
    pub layer_thickness_m: f64,
    /// Original caller-declared uncertainty, never promoted to certification.
    pub uncertainty: InputUncertainty,
}

/// Frozen, sourced point contact suitable for any mass/potential/port solver.
///
/// The selected law has F_elastic = K delta^(3/2), U = 2 K delta^(5/2)/5,
/// and F_loss = chi F_elastic delta_dot. K is obtained from the EXISTING
/// finite-patch evaluator, including its true elliptic shape calculation.
/// Consumers must retain this owner and enforce its envelope on accepted
/// configurations and on the actual rate used by their dissipative port.
#[derive(Clone, Debug, PartialEq)]
pub struct SourcePointContact {
    model: BoundNormalContactModel,
    geometry: PointContactGeometry,
    envelope: PointContactEnvelope,
    template: NormalPatchRequest,
    stiffness_n_per_m_3_2: f64,
    dissipation_s_per_m: f64,
    limiting_response: PointNormalPatchReceipt,
}

impl SourcePointContact {
    /// Compile one explicit contact from both bulk cards and the interface card.
    ///
    /// Checks both rate extremes at maximum indentation. Patch size, strain,
    /// depth/layer ratios and peak pressure are monotone on this nonadhesive
    /// Hertz operating box; the negative-rate corner also checks that the
    /// Hunt--Crossley total force cannot turn tensile. Nothing is clamped.
    /// Material state, geometry and constitutive choice stay frozen afterwards.
    pub fn new(model: BoundNormalContactModel, geometry: PointContactGeometry,
        envelope: PointContactEnvelope) -> Result<Self, NormalPatchError> {
        for (field, value) in [
            ("maximum_indentation_m", envelope.maximum_indentation_m),
            ("maximum_step_s", envelope.maximum_step_s),
            ("half_space_depth_m", envelope.half_space_depth_m),
            ("layer_thickness_m", envelope.layer_thickness_m),
        ] {
            if !value.is_finite() || value <= 0.0 {
                return Err(NormalPatchError::InvalidInput { field });
            }
        }
        if !envelope.maximum_absolute_rate_m_per_s.is_finite()
            || envelope.maximum_absolute_rate_m_per_s < 0.0 {
            return Err(NormalPatchError::InvalidInput { field: "maximum_absolute_rate_m_per_s" });
        }
        let normal = model.normal_state();
        let elastic = normal.elastic();
        let temperature_k = elastic.state_coordinate("T")
            .ok_or(NormalPatchError::InvalidInput { field: "source temperature coordinate T" })?;
        let modulus = elastic.reduced_modulus_pa();
        let chi = match model.law() {
            BoundNormalContactLaw::ElasticHertz => 0.0,
            BoundNormalContactLaw::HuntCrossleyPoint { dissipation_s_per_m } => dissipation_s_per_m,
        };
        let (patch_geometry, law) = match (geometry, model.law()) {
            (PointContactGeometry::SpherePlane { effective_radius_m }, BoundNormalContactLaw::ElasticHertz) =>
                (NormalPatchGeometry::SpherePlane, NormalPatchLaw::HertzSpherePlane {
                    effective_radius_m, reduced_modulus_pa: modulus }),
            (PointContactGeometry::SpherePlane { effective_radius_m }, BoundNormalContactLaw::HuntCrossleyPoint { .. }) =>
                (NormalPatchGeometry::SpherePlane, NormalPatchLaw::HuntCrossleySphere {
                    effective_radius_m, reduced_modulus_pa: modulus, dissipation_s_per_m: chi }),
            (PointContactGeometry::EllipticParaboloid { maximum_curvature_m_inverse, minimum_curvature_m_inverse },
                BoundNormalContactLaw::ElasticHertz) =>
                (NormalPatchGeometry::EllipticParaboloid, NormalPatchLaw::HertzEllipticParaboloid {
                    maximum_principal_curvature_m_inverse: maximum_curvature_m_inverse,
                    minimum_principal_curvature_m_inverse: minimum_curvature_m_inverse,
                    reduced_modulus_pa: modulus }),
            (PointContactGeometry::EllipticParaboloid { maximum_curvature_m_inverse, minimum_curvature_m_inverse },
                BoundNormalContactLaw::HuntCrossleyPoint { .. }) =>
                (NormalPatchGeometry::EllipticParaboloid, NormalPatchLaw::HuntCrossleyEllipticParaboloid {
                    maximum_principal_curvature_m_inverse: maximum_curvature_m_inverse,
                    minimum_principal_curvature_m_inverse: minimum_curvature_m_inverse,
                    reduced_modulus_pa: modulus, dissipation_s_per_m: chi }),
        };
        let card = model.model_card();
        let template = NormalPatchRequest {
            identity: NormalPatchIdentity {
                model_id: format!("{}:v{}:{}", card.law.0, card.law_version, card.content_hash()),
                source_id: format!("fs-contact/normal-model-binding/{}", model.identity()),
                state_id: "source-point-contact/elastic-reference".into(),
            },
            interface: elastic.interface().interface().clone(), law, geometry: patch_geometry,
            indentation_m: envelope.maximum_indentation_m, indentation_rate_m_per_s: 0.0,
            step_s: envelope.maximum_step_s, line_load_n_per_m: 0.0,
            applicability: ApplicabilityInput {
                half_space_depth_m: envelope.half_space_depth_m,
                layer_thickness_m: envelope.layer_thickness_m,
                yield_strength_pa: elastic.limiting_yield_stress_pa(),
                characteristic_rate_m_per_s: model.characteristic_rate_m_per_s(),
                temperature_k, adhesion_energy_j_per_m2: normal.adhesion_energy_j_per_m2(),
            },
            limits: model.limits(), uncertainty: envelope.uncertainty,
        };
        let reference = point(&template)?;
        let d = envelope.maximum_indentation_m;
        let denominator = d * d.sqrt();
        let stiffness = reference.normal_force_n / denominator;
        if !denominator.is_finite() || denominator <= 0.0 || !stiffness.is_finite() || stiffness <= 0.0 {
            return Err(NormalPatchError::Overflow { field: "source point Hertz coefficient" });
        }
        let mut corner = template.clone();
        corner.indentation_rate_m_per_s = -envelope.maximum_absolute_rate_m_per_s;
        corner.identity.state_id = "source-point-contact/separation-limit".into();
        point(&corner)?;
        corner.indentation_rate_m_per_s = envelope.maximum_absolute_rate_m_per_s;
        corner.identity.state_id = "source-point-contact/closure-limit".into();
        let limiting_response = point(&corner)?;
        Ok(Self { model, geometry, envelope, template,
            stiffness_n_per_m_3_2: stiffness, dissipation_s_per_m: chi, limiting_response })
    }

    /// Original immutable bulk, interface, model and source-property receipts.
    pub const fn model(&self) -> &BoundNormalContactModel { &self.model }
    /// Geometry used to obtain the force coefficient.
    pub const fn geometry(&self) -> PointContactGeometry { self.geometry }
    /// Operating box already checked against the selected card's limits.
    pub const fn envelope(&self) -> PointContactEnvelope { self.envelope }
    /// K [N/m^(3/2)] for the original Hertz potential, not a fitted spring.
    pub const fn stiffness_n_per_m_3_2(&self) -> f64 { self.stiffness_n_per_m_3_2 }
    /// Card-selected chi [s/m]. Zero is the elastic law.
    pub const fn dissipation_s_per_m(&self) -> f64 { self.dissipation_s_per_m }
    /// Worst closure corner's original finite-patch response and authority.
    pub const fn limiting_response(&self) -> &PointNormalPatchReceipt { &self.limiting_response }

    /// Allocation-free operating-box check for a numerical consumer.
    /// A negative geometric gap must be converted to positive indentation by
    /// the geometry/solver owner; this method never silently clips its inputs.
    pub fn validate_kinematics(&self, indentation_m: f64, rate_m_per_s: f64,
        step_s: f64) -> Result<(), NormalPatchError> {
        if !indentation_m.is_finite() || indentation_m < 0.0 {
            return Err(NormalPatchError::InvalidInput { field: "source point indentation" });
        }
        if !rate_m_per_s.is_finite() {
            return Err(NormalPatchError::InvalidInput { field: "source point rate" });
        }
        if !step_s.is_finite() || step_s <= 0.0 {
            return Err(NormalPatchError::InvalidInput { field: "source point step" });
        }
        for (ratio, value, limit) in [
            ("source point indentation", indentation_m, self.envelope.maximum_indentation_m),
            ("source point absolute rate", rate_m_per_s.abs(), self.envelope.maximum_absolute_rate_m_per_s),
            ("source point sample duration", step_s, self.envelope.maximum_step_s),
        ] {
            if value > limit { return Err(NormalPatchError::OutsideApplicability { ratio, value, limit }); }
        }
        Ok(())
    }

    /// Evaluate pressure, footprint, force, tangent and loss through the original
    /// law owner. This cold observation allocates receipt identities; the
    /// numerical path can use the admitted coefficient and allocation-free guard.
    /// Work here is instantaneous power times step, NOT a time-integrated loss.
    pub fn evaluate(&self, state_id: &str, indentation_m: f64, rate_m_per_s: f64,
        step_s: f64) -> Result<PointNormalPatchReceipt, NormalPatchError> {
        self.validate_kinematics(indentation_m, rate_m_per_s, step_s)?;
        let mut request = self.template.clone();
        request.identity.state_id = state_id.into();
        request.indentation_m = indentation_m;
        request.indentation_rate_m_per_s = rate_m_per_s;
        request.step_s = step_s;
        point(&request)
    }
}

fn point(request: &NormalPatchRequest) -> Result<PointNormalPatchReceipt, NormalPatchError> {
    match request.evaluate()? {
        NormalPatchReceipt::Point(receipt) => Ok(receipt),
        NormalPatchReceipt::Line(_) => Err(NormalPatchError::GeometryLawMismatch { geometry: request.geometry }),
    }
}
