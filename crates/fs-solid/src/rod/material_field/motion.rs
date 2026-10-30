//! Frozen source-bound Cosserat configurations for native dynamic consumers.
//!
//! The mechanical owner is shared by centered and spatial-inertia views. Both
//! use ORIGINAL reference inventory, physical section rotary inertia and the
//! same analytical force/source/support admission. Neither is a time integrator.
//! The centered entry point retains its original API and refuses eccentric mass;
//! [`SpatialMaterialRod::into_spatial_motion_model`] retains that mass's centroid
//! and complete central tensor instead of discarding translation/rotation coupling.

use std::sync::Arc;
use fs_exec::CancelGate;
use fs_ga::So3;
use fs_material::section::ConcentricShearSection;
use fs_material::state_point::{DENSITY_PROPERTY, YOUNG_MODULUS_PROPERTY};
use crate::rod::{Rod, derivatives::RodEnergy};
use super::{SpatialMaterialRod, SpatialRodError, checkpoint};

/// Conservative ambient body loads and fluid-state refresh on the same model.
pub mod environment;
pub use environment::RodBodyLoads;

/// Temperature-only source updates retaining the original dynamic inventory.
pub mod thermal;

pub mod spatial {
//! Reference centroid and full rotary tensor for eccentric rod attachments.
//!
//! A strand contributes its declared half-mass at each rotating attachment.
//! Successive central-inertia merges avoid subtracting two large parallel-axis
//! tensors. This is still endpoint lumping, not distributed sagging-cable mass.
//! General rigid-inertia validation and principal-axis lowering remain with the
//! existing multibody consumers; no eigensolver or motion law is added here.

use fs_exec::CancelGate;
use super::{PreparedRodMotion, RodNodalInertia, RodRootBoundary};
use super::super::{SpatialMaterialRod, SpatialRodError, checkpoint};

pub type PreparedSpatialRodMotion = PreparedRodMotion<RodNodalSpatialInertia>;

/// Material-frame centroid offset and rotary inertia ABOUT that centroid.
/// Products of inertia are retained. Mass never follows the deformed length.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RodNodalSpatialInertia {
    mass_kg: f64,
    center_of_mass_body_m: [f64; 3],
    inertia_com_body_kg_m2: [[f64; 3]; 3],
}
impl RodNodalSpatialInertia {
    pub const fn mass_kg(&self) -> f64 { self.mass_kg }
    /// Offset from the material node to its lumped centroid, material axes [m].
    pub const fn center_of_mass_body_m(&self) -> [f64; 3] { self.center_of_mass_body_m }
    /// Symmetric central tensor in MATERIAL axes [kg m2], not principal moments.
    pub const fn inertia_com_body_kg_m2(&self) -> [[f64; 3]; 3] { self.inertia_com_body_kg_m2 }

    fn centered(inertia: RodNodalInertia) -> Self {
        let [x, y, z] = inertia.principal_kg_m2();
        Self { mass_kg: inertia.mass_kg(), center_of_mass_body_m: [0.0; 3],
            inertia_com_body_kg_m2: [[x, 0.0, 0.0], [0.0, y, 0.0], [0.0, 0.0, z]] }
    }

    // I_AB,c = I_A,c + (m_A*m_B/(m_A+m_B)) [d]x [d]x^T for point B.
    // The reduced mass and scaled offsets avoid m_A*m_B and m*r*r overflow.
    fn add_point_mass(&mut self, mass_kg: f64, offset_body_m: [f64; 3])
        -> Result<(), SpatialRodError> {
        let total = self.mass_kg + mass_kg;
        if !mass_kg.is_finite() || mass_kg <= 0.0 || !total.is_finite()
            || offset_body_m.iter().any(|v| !v.is_finite()) {
            return Err(SpatialRodError::Invalid("finite positive attachment inventory required"));
        }
        let weight = mass_kg / total;
        let reduced = if mass_kg <= self.mass_kg {
            mass_kg * (self.mass_kg / total)
        } else { self.mass_kg * weight };
        let delta: [f64; 3] = std::array::from_fn(|i| offset_body_m[i] - self.center_of_mass_body_m[i]);
        let factor = fs_math::det::sqrt(reduced);
        let q = delta.map(|d| factor * d);
        let center: [f64; 3] = std::array::from_fn(|i| weight.mul_add(delta[i], self.center_of_mass_body_m[i]));
        let mut tensor = self.inertia_com_body_kg_m2;
        for i in 0..3 {
            let j = (i + 1) % 3;
            let k = (i + 2) % 3;
            tensor[i][i] += q[j] * q[j] + q[k] * q[k];
            for j in i + 1..3 {
                tensor[i][j] -= q[i] * q[j];
                tensor[j][i] = tensor[i][j];
            }
        }
        if !reduced.is_finite() || reduced <= 0.0
            || delta.iter().chain(&q).chain(&center).any(|v| !v.is_finite())
            || tensor.iter().flatten().any(|v| !v.is_finite()) {
            return Err(SpatialRodError::Invalid("eccentric nodal inertia is unrepresentable"));
        }
        self.mass_kg = total;
        self.center_of_mass_body_m = center;
        self.inertia_com_body_kg_m2 = tensor;
        Ok(())
    }
}

impl SpatialMaterialRod {
    /// Freeze the same mechanics with full eccentric attachment inventory.
    /// Each link keeps its original half-mass at each body-frame offset. The
    /// original centered API remains strict rather than returning false moments.
    /// Source validity, support geometry and all reference receipts are shared.
    /// A time consumer must separately admit its supported constraint class.
    pub fn into_spatial_motion_model(self, root: RodRootBoundary, gate: &CancelGate)
        -> Result<PreparedSpatialRodMotion, SpatialRodError> {
        self.prepare_motion(root, gate, |source, reference| {
            let mut inertia: Vec<_> = reference.into_iter().map(RodNodalSpatialInertia::centered).collect();
            for link in &source.tension_links {
                checkpoint(gate)?;
                for attachment in link.attachments() {
                    inertia[attachment.node].add_point_mass(0.5 * link.reference_mass_kg(), attachment.offset_body_m)?;
                }
            }
            checkpoint(gate)?;
            Ok(inertia)
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn central_merge_retains_products_of_inertia_and_direct_kinetic_energy() {
        let base = RodNodalInertia { mass_kg: 2.0, principal_kg_m2: [0.4, 0.3, 0.2] };
        let mut combined = RodNodalSpatialInertia::centered(base);
        let points = [(0.5, [0.3, -0.2, 0.1]), (1.3, [-0.1, 0.4, 0.2])];
        for (m, r) in points { combined.add_point_mass(m, r).unwrap(); }
        let center = combined.center_of_mass_body_m();
        assert_eq!(combined.mass_kg(), 3.8);
        assert!(combined.inertia_com_body_kg_m2()[0][1] != 0.0);
        for case in 0..32 {
            let v = [0.03 * case as f64, -0.2, 0.1];
            let w = [0.2, 0.01 * case as f64, -0.3];
            let cross = |r: [f64; 3]| [w[1]*r[2]-w[2]*r[1], w[2]*r[0]-w[0]*r[2], w[0]*r[1]-w[1]*r[0]];
            let speed2 = |r| cross(r).into_iter().zip(v).map(|(a,b)| (a+b)*(a+b)).sum::<f64>();
            let direct = 0.5 * base.mass_kg() * speed2([0.0; 3])
                + 0.5 * (0..3).map(|i| base.principal_kg_m2()[i]*w[i]*w[i]).sum::<f64>()
                + points.into_iter().map(|(m,r)| 0.5*m*speed2(r)).sum::<f64>();
            let tensor = combined.inertia_com_body_kg_m2();
            let reduced = 0.5 * combined.mass_kg() * speed2(center)
                + 0.5 * (0..3).map(|i| (0..3).map(|j| w[i]*tensor[i][j]*w[j]).sum::<f64>()).sum::<f64>();
            assert!((direct-reduced).abs() <= 2e-14*direct.max(1.0));
        }
    }

    #[test]
    fn symmetric_eccentric_mass_keeps_its_rotary_inertia_even_with_zero_centroid() {
        let mut i = RodNodalSpatialInertia::centered(RodNodalInertia {
            mass_kg: 2.0, principal_kg_m2: [0.4, 0.3, 0.2] });
        i.add_point_mass(1.0, [0.0, 0.3, 0.0]).unwrap();
        i.add_point_mass(1.0, [0.0, -0.3, 0.0]).unwrap();
        assert!(i.center_of_mass_body_m().iter().all(|c| c.abs() < 1e-16));
        assert!((i.inertia_com_body_kg_m2()[0][0] - 0.58).abs() < 1e-15);
        assert!((i.inertia_com_body_kg_m2()[2][2] - 0.38).abs() < 1e-15);
    }
}

}
pub use spatial::{PreparedSpatialRodMotion, RodNodalSpatialInertia};

/// Explicit root policy; releasing the static root never happens implicitly.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RodRootBoundary { Clamped, Free }

/// Caller-owned candidate; every evaluation validates dimensions and supports.
#[derive(Clone, Debug, PartialEq)]
pub struct RodMotionPose {
    pub positions_m: Vec<[f64; 3]>,
    pub frames: Vec<So3>,
}

/// Centered row-lumped mass and section principal inertia in material axes.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RodNodalInertia {
    mass_kg: f64,
    principal_kg_m2: [f64; 3],
}
impl RodNodalInertia {
    pub const fn mass_kg(&self) -> f64 { self.mass_kg }
    pub const fn principal_kg_m2(&self) -> [f64; 3] { self.principal_kg_m2 }
}

/// Complete native force, including fixed nodes, for impulse/reaction accounting.
#[derive(Clone, Debug, PartialEq)]
pub struct RodMotionResponse {
    /// Elastic/link storage plus the original gravity and hydrostatic potentials [J].
    pub potential_energy_j: f64,
    /// Force ON each material node, world frame [N].
    pub force_n: Vec<[f64; 3]>,
    /// Couple ABOUT the material node, expressed in its body frame [N m].
    /// A consumer moving about a different centroid must shift this wrench.
    pub torque_body_nm: Vec<[f64; 3]>,
}

/// Immutable sources, supports and reference inertia shared across a trajectory.
///
/// The default inertia is the original centered model. The spatial alias changes
/// only the inventory representation; force evaluation and material coordinates
/// are identical. No arbitrary externally supplied inertia can be installed.
#[derive(Clone, Debug)]
pub struct PreparedRodMotion<I = RodNodalInertia> {
    source: Arc<SpatialMaterialRod>,
    inertia: Vec<I>,
    fixed_position: Vec<[bool; 3]>,
    fixed_frame: Vec<bool>,
}

impl SpatialMaterialRod {
    /// Freeze the original centered-inertia lane. Eccentric links still refuse;
    /// use `into_spatial_motion_model` to preserve their full mass properties.
    /// Finite source-frequency bands cannot authorize arbitrary transients.
    pub fn into_motion_model(self, root: RodRootBoundary, gate: &CancelGate)
        -> Result<PreparedRodMotion, SpatialRodError> {
        self.prepare_motion(root, gate, |source, mut inertia| {
            for link in &source.tension_links {
                checkpoint(gate)?;
                if link.attachments().iter().flat_map(|a| a.offset_body_m).any(|v| v != 0.0) {
                    return Err(SpatialRodError::Invalid(
                        "dynamic offset links require the spatial-inertia motion model"));
                }
                for a in link.attachments() {
                    inertia[a.node].mass_kg += 0.5 * link.reference_mass_kg();
                }
            }
            if inertia.iter().any(|i| !i.mass_kg.is_finite() || i.mass_kg <= 0.0) {
                return Err(SpatialRodError::Invalid("positive finite nodal reference mass required"));
            }
            Ok(inertia)
        })
    }

    // One cold-path source/reference/support owner for both inertia views.
    fn prepare_motion<I>(self, root: RodRootBoundary, gate: &CancelGate,
        inventory: impl FnOnce(&Self, Vec<RodNodalInertia>) -> Result<Vec<I>, SpatialRodError>)
        -> Result<PreparedRodMotion<I>, SpatialRodError> {
        checkpoint(gate)?;
        for section in &self.sections { admit_section(section)?; checkpoint(gate)?; }
        let reference = self.reference_mass();
        let mut inertia: Vec<_> = reference.nodal_mass_kg().iter().map(|&mass_kg|
            RodNodalInertia { mass_kg, principal_kg_m2: [0.0; 3] }).collect();
        let half_length = 0.5 * reference.segment_length_m();
        for (segment, section) in reference.sections().iter().enumerate() {
            checkpoint(gate)?;
            admit_section(section)?;
            let mut diameter = 0.0;
            for (layer, geometry) in section.elastic().layers().iter()
                .zip(section.elastic().layer_integrals()) {
                let rho = layer.material.property(DENSITY_PROPERTY)
                    .ok_or(SpatialRodError::Invalid("reference inertia lost its original density"))?;
                let value = half_length * (rho.value_si() * geometry.second_moment_m4());
                if !value.is_finite() || value <= 0.0 {
                    return Err(SpatialRodError::Invalid("reference rotary inertia is unrepresentable"));
                }
                diameter += value;
            }
            for node in [segment, segment + 1] {
                for (sum, value) in inertia[node].principal_kg_m2.iter_mut()
                    .zip([2.0 * diameter, diameter, diameter]) { *sum += value; }
            }
        }
        for link in &self.tension_links {
            checkpoint(gate)?;
            for layer in link.section().layers() {
                for name in [DENSITY_PROPERTY, YOUNG_MODULUS_PROPERTY] {
                    admit_property(&layer.material, name)?;
                }
            }
        }
        if inertia.iter().any(|i| !i.mass_kg.is_finite() || i.mass_kg <= 0.0
            || i.principal_kg_m2.iter().any(|v| !v.is_finite() || *v <= 0.0)) {
            return Err(SpatialRodError::Invalid("positive finite nodal mass and rotary inertia required"));
        }
        let inertia = inventory(&self, inertia)?;
        let mut fixed_position = vec![[false; 3]; inertia.len()];
        let mut fixed_frame = vec![false; inertia.len()];
        if root == RodRootBoundary::Clamped { fixed_position[0] = [true; 3]; fixed_frame[0] = true; }
        for support in &self.supports {
            fixed_position[support.node] = support.position_m.map(|v| v.is_some());
            fixed_frame[support.node] = support.frame.is_some();
        }
        let model = PreparedRodMotion { source: Arc::new(self), inertia, fixed_position, fixed_frame };
        model.evaluate(&model.initial_pose(), gate)?;
        checkpoint(gate)?;
        Ok(model)
    }
}

fn admit_section(section: &ConcentricShearSection) -> Result<(), SpatialRodError> {
    for layer in section.elastic().layers() {
        for name in [DENSITY_PROPERTY, YOUNG_MODULUS_PROPERTY, section.shear_property()] {
            admit_property(&layer.material, name)?;
        }
    }
    Ok(())
}
fn admit_property(material: &fs_material::state_point::ResolvedMaterialStatePoint, name: &str)
    -> Result<(), SpatialRodError> {
    let property = material.property(name)
        .ok_or(SpatialRodError::Invalid("dynamic rod lost an original section property"))?;
    if !property.is_constant_at_fixed_temperature()
        || property.angular_frequency_band_rad_s()?.is_some_and(|(lo, hi)| lo > 0.0 || hi < f64::MAX) {
        return Err(SpatialRodError::Invalid(
            "transient rod needs frozen elastic coefficients without a finite frequency validity band"));
    }
    Ok(())
}

impl<I> PreparedRodMotion<I> {
    pub fn source(&self) -> &SpatialMaterialRod { &self.source }
    pub fn inertia(&self) -> &[I] { &self.inertia }
    pub fn fixed_positions(&self) -> &[[bool; 3]] { &self.fixed_position }
    pub fn fixed_frames(&self) -> &[bool] { &self.fixed_frame }
    pub fn initial_pose(&self) -> RodMotionPose {
        RodMotionPose { positions_m: self.source.rod.positions.clone(), frames: self.source.rod.frames.clone() }
    }

    /// Evaluate native forces ABOUT material nodes, never mass centroids.
    /// No tangent is allocated. Every source, support and branch-cut check is
    /// shared by both inertia views; cancellation retains no candidate state.
    pub fn evaluate(&self, pose: &RodMotionPose, gate: &CancelGate)
        -> Result<RodMotionResponse, SpatialRodError> {
        let rod = self.rod_at_pose(pose, gate)?;
        let energy = super::derivatives::SpatialEnergy::new(&self.source);
        let potential_energy_j = energy.value(&rod)?;
        let response = energy.linearize(&rod, &[], false, &|| checkpoint(gate))?
            .ok_or(SpatialRodError::Invalid("native dynamic rod requires analytical forces"))?;
        let force_n = response.gradient.iter().map(|g| [-g[0], -g[1], -g[2]]).collect();
        let torque_body_nm = response.gradient.iter().map(|g| [-g[3], -g[4], -g[5]]).collect();
        checkpoint(gate)?;
        Ok(RodMotionResponse { potential_energy_j, force_n, torque_body_nm })
    }

    // All material/body-load observations use the same current pose and support
    // admission. A source's retained static pose must not substitute for motion.
    fn rod_at_pose(&self, pose: &RodMotionPose, gate: &CancelGate)
        -> Result<Rod, SpatialRodError> {
        checkpoint(gate)?;
        let count = self.inertia.len();
        if pose.positions_m.len() != count || pose.frames.len() != count
            || pose.positions_m.iter().flatten().any(|v| !v.is_finite()) {
            return Err(SpatialRodError::Invalid("dynamic rod pose must match the finite reference mesh"));
        }
        for node in 0..count {
            checkpoint(gate)?;
            for c in 0..3 {
                if self.fixed_position[node][c]
                    && pose.positions_m[node][c] != self.source.rod.positions[node][c] {
                    return Err(SpatialRodError::Invalid("dynamic rod moved a stationary position support"));
                }
            }
            if self.fixed_frame[node] && pose.frames[node] != self.source.rod.frames[node] {
                return Err(SpatialRodError::Invalid("dynamic rod moved a stationary frame support"));
            }
        }
        checkpoint(gate)?;
        Ok(Rod { positions: pose.positions_m.clone(), frames: pose.frames.clone(),
            l0: self.source.rod.l0, section: self.source.rod.section })
    }
}

#[cfg(test)]
mod source_tests {
//! G0/G3/G4: actual source-bound eccentric inventory and native force ownership.
use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer, ConcentricShearSection};
use fs_material::state_point::{MaterialPropertySelection, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point};
use fs_qty::{Density, Pressure};
use crate::rod::material_field::{SpatialMaterialRod, SpatialRodError};
use crate::rod::material_field::motion::RodRootBoundary;
use crate::rod::material_field::tension::{MaterialTensionLink, TensionLinkAttachment};

fn section() -> ConcentricShearSection {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name,value,dims) in [("density",1200.0,Density::DIMS),
        ("young_modulus",1e6,Pressure::DIMS),("shear_modulus",4e5,Pressure::DIMS)] {
        claims.insert_claim(PropertyClaim { key:PropertyKey::new(name,dims),
            value:PropertyValue::Scalar { value,dims }, validity:ValidityDomain::unconstrained(),
            interpolation:InterpolationPolicy::ConstantWithinValidity, uncertainty:UncertaintyModel::Unstated,
            observations:vec![],provenance:Provenance { source:"manufactured eccentric rod".into(),
                license:"test-only".into(),artifact:None } }).unwrap();
        requirements.push(ScalarPropertyRequirement::try_new(name,dims,ScalarAdmissibility::StrictlyPositive).unwrap());
    }
    let card=MaterialCard::assemble(MaterialStateId { chemistry:"manufactured".into(),
        phase:"solid".into(),process:"test".into(),revision:0 },claims,vec![]).unwrap();
    let material=resolve_material_state_point(&card,&QueryPoint::new(),&requirements,
        MaterialPropertySelection::SingleClaimOnly).unwrap();
    ConcentricElasticSection::try_new(0.0,vec![ConcentricMaterialLayer { outer_radius_m:0.02,material }])
        .unwrap().with_shear_modulus("shear_modulus",0.8).unwrap()
}
fn specimen(offset:[f64;3]) -> SpatialMaterialRod {
    let section=section();
    let link=MaterialTensionLink::try_new(section.elastic().clone(),[
        TensionLinkAttachment { node:0,offset_body_m:offset },
        TensionLinkAttachment { node:2,offset_body_m:offset },
    ],1.0,1,0.2).unwrap();
    SpatialMaterialRod::try_straight(1.0,vec![section;2]).unwrap()
        .with_tension_links(vec![link],&CancelGate::new()).unwrap()
}
fn near(a:f64,b:f64) { assert!((a-b).abs()<=4e-14*a.abs().max(b.abs()).max(1.0),"{a:e} != {b:e}"); }

#[test]
fn g0_offset_inventory_matches_direct_centroid_and_parallel_axis_tensor() {
    let gate=CancelGate::new();
    let offset=[0.03,-0.04,0.02];
    let rod=specimen(offset);
    let link_half=0.5*rod.tension_links()[0].reference_mass_kg();
    let bare=SpatialMaterialRod::try_straight(1.0,vec![section();2]).unwrap()
        .into_motion_model(RodRootBoundary::Free,&gate).unwrap();
    let expected_total=rod.total_reference_mass_kg();
    let model=rod.into_spatial_motion_model(RodRootBoundary::Free,&gate).unwrap();
    near(model.inertia().iter().map(|i|i.mass_kg()).sum(),expected_total);
    for node in [0,2] {
        let original=bare.inertia()[node];
        let added=model.inertia()[node];
        let reduced=original.mass_kg()*link_half/(original.mass_kg()+link_half);
        let square=offset.iter().map(|r|r*r).sum::<f64>();
        for i in 0..3 {
            near(added.center_of_mass_body_m()[i],link_half*offset[i]/added.mass_kg());
            for j in 0..3 {
                let expected=if i==j { original.principal_kg_m2()[i]+reduced*(square-offset[i]*offset[i]) }
                    else { -reduced*offset[i]*offset[j] };
                near(added.inertia_com_body_kg_m2()[i][j],expected);
            }
        }
    }
    assert_eq!(model.inertia()[1].center_of_mass_body_m(),[0.0;3]);
}

#[test]
fn g3_centered_view_is_exactly_compatible_and_uses_identical_native_forces() {
    let gate=CancelGate::new();
    let rod=specimen([0.0;3]);
    let a=rod.clone().into_motion_model(RodRootBoundary::Free,&gate).unwrap();
    let b=rod.into_spatial_motion_model(RodRootBoundary::Free,&gate).unwrap();
    for (x,y) in a.inertia().iter().zip(b.inertia()) {
        assert_eq!(x.mass_kg(),y.mass_kg());
        assert_eq!(y.center_of_mass_body_m(),[0.0;3]);
        for i in 0..3 { for j in 0..3 {
            assert_eq!(y.inertia_com_body_kg_m2()[i][j],if i==j {x.principal_kg_m2()[i]} else {0.0});
        }}
    }
    let mut pose=a.initial_pose();
    pose.positions_m[2][0]+=0.005;
    assert_eq!(a.evaluate(&pose,&gate).unwrap(),b.evaluate(&pose,&gate).unwrap());
}

#[test]
fn g0_gravity_force_and_offset_couple_use_the_same_inventory_once() {
    let gate=CancelGate::new();
    let offset=[0.03,-0.04,0.02];
    let g=[0.0,-9.0,1.0];
    let rod=specimen(offset).with_gravity(g,&gate).unwrap();
    let model=rod.into_spatial_motion_model(RodRootBoundary::Free,&gate).unwrap();
    let response=model.evaluate(&model.initial_pose(),&gate).unwrap();
    for (node,i) in model.inertia().iter().enumerate() {
        let force=g.map(|v|i.mass_kg()*v);
        let c=i.center_of_mass_body_m();
        let moment=[c[1]*force[2]-c[2]*force[1],c[2]*force[0]-c[0]*force[2],c[0]*force[1]-c[1]*force[0]];
        for axis in 0..3 { near(response.force_n[node][axis],force[axis]); near(response.torque_body_nm[node][axis],moment[axis]); }
    }
}

#[test]
fn g4_centered_refusal_support_admission_and_cancellation_remain_explicit() {
    let gate=CancelGate::new();
    assert!(specimen([0.0,0.01,0.0]).into_motion_model(RodRootBoundary::Free,&gate).is_err());
    let model=specimen([0.0,0.01,0.0]).into_spatial_motion_model(RodRootBoundary::Clamped,&gate).unwrap();
    let mut pose=model.initial_pose();
    pose.positions_m[0][0]+=0.01;
    assert!(model.evaluate(&pose,&gate).is_err());
    gate.request();
    assert_eq!(model.evaluate(&model.initial_pose(),&gate),Err(SpatialRodError::Cancelled));
    assert!(specimen([0.0,0.01,0.0]).into_spatial_motion_model(RodRootBoundary::Free,&gate).is_err());
}

}
