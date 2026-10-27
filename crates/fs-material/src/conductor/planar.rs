//! Planar Ohmic response from scalar or fixed-principal-axis material data.
//!
//! Both in-plane principal coefficients are independent source requirements.
//! Their fixed frame maps into mesh x/y by `sigma = Q diag(sigma1,sigma2) Q^T`.
//! Resistivity sources are evaluated before inversion in their principal frame;
//! componentwise reciprocal of a rotated matrix is NOT tensor inversion.
//! This is a symmetric positive DC law, without Hall response, electrical
//! memory, rotating temperature-dependent eigenvectors or inferred normal
//! transport. The physical plane must be decoupled from the omitted direction.

use super::*;
use fs_blake3::{ContentHash, DomainHasher};

/// Material-axis order 1, 2; both independent principal values are required.
pub const PRINCIPAL_CONDUCTIVITY_PROPERTIES: [&str; 2] =
    ["electrical_conductivity_1", "electrical_conductivity_2"];
/// Material-axis order 1, 2, in ohm metres.
pub const PRINCIPAL_RESISTIVITY_PROPERTIES: [&str; 2] =
    ["electrical_resistivity_1", "electrical_resistivity_2"];

/// Explicit, fixed, proper in-plane material-to-geometry rotation.
/// A finite angle generates its orthonormal basis through the common tensor
/// rotation owner; arbitrary unvalidated matrix entries are never admitted.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlanarElectricalFrame {
    material_frame: ContentHash,
    geometry_frame: ContentHash,
    angle_rad: f64,
}
impl PlanarElectricalFrame {
    /// Material axis 1 is counterclockwise from geometry x, in radians.
    /// Both frame identities must be nonzero; equal identities require zero
    /// rotation. A mesh consumer owns the correspondence to its coordinates.
    pub fn new(
        material_frame: ContentHash,
        geometry_frame: ContentHash,
        angle_rad: f64,
    ) -> Result<Self, ConductorError> {
        if material_frame == ContentHash([0; 32])
            || geometry_frame == ContentHash([0; 32])
            || !angle_rad.is_finite()
            || (material_frame == geometry_frame && angle_rad != 0.0)
        {
            return Err(ConductorError::InvalidInput {
                quantity: "named planar electrical frames and finite proper rotation",
            });
        }
        Ok(Self {
            material_frame,
            geometry_frame,
            angle_rad,
        })
    }
    pub const fn material_frame(&self) -> ContentHash {
        self.material_frame
    }
    pub const fn geometry_frame(&self) -> ContentHash {
        self.geometry_frame
    }
    pub const fn angle_rad(&self) -> f64 {
        self.angle_rad
    }
}

/// Native conductivity coefficient, in S/m, retaining scalar arithmetic when
/// the input law is scalar. Tensor symmetry does not imply isotropy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum PlanarConductivity {
    Scalar(f64),
    Tensor([[f64; 2]; 2]),
}
impl PlanarConductivity {
    pub fn matrix_s_m(self) -> [[f64; 2]; 2] {
        match self {
            Self::Scalar(s) => [[s, 0.0], [0.0, s]],
            Self::Tensor(s) => s,
        }
    }
    /// Physical current density J=sigma E. Field owners retain their numerical
    /// representability checks and independently integrate E.J and source work.
    pub fn current_density_a_m2(self, electric_field_v_m: [f64; 2]) -> [f64; 2] {
        let [x, y] = electric_field_v_m;
        match self {
            Self::Scalar(s) => [s * x, s * y],
            Self::Tensor(s) => [
                s[0][0].mul_add(x, s[0][1] * y),
                s[1][0].mul_add(x, s[1][1] * y),
            ],
        }
    }
}

mod sealed {
    pub trait Sealed {}
}

/// Shared source-bearing response for the existing planar field and heat
/// owners. Implementations are the canonical scalar, orthotropic and mixed
/// source objects below; a bare numerical tensor does not acquire authority.
pub trait PlanarConductorMaterial: sealed::Sealed + Clone {
    fn material(&self) -> &ResolvedMaterialStatePoint;
    fn source_basis(&self) -> ElectricalPropertyBasis;
    fn frame(&self) -> Option<PlanarElectricalFrame>;
    fn source_property_names(&self) -> &'static [&'static str];
    fn conductivity(&self) -> PlanarConductivity;
    /// Rebuild the SAME law/basis/orientation around an atomically refreshed
    /// original material bundle. Never choose alternate properties by name.
    fn with_material(&self, material: ResolvedMaterialStatePoint) -> Result<Self, ConductorError>;
    fn source_properties(&self) -> impl Iterator<Item = &ResolvedScalarProperty> {
        self.source_property_names().iter().map(|name| {
            self.material()
                .property(name)
                .expect("canonical conductor retains every admitted electrical property")
        })
    }
}

impl sealed::Sealed for ResolvedConductorState {}
impl PlanarConductorMaterial for ResolvedConductorState {
    fn material(&self) -> &ResolvedMaterialStatePoint {
        self.material()
    }
    fn source_basis(&self) -> ElectricalPropertyBasis {
        self.source_basis()
    }
    fn frame(&self) -> Option<PlanarElectricalFrame> {
        None
    }
    fn source_property_names(&self) -> &'static [&'static str] {
        match self.source_basis() {
            ElectricalPropertyBasis::Resistivity => &[ELECTRICAL_RESISTIVITY_PROPERTY],
            ElectricalPropertyBasis::Conductivity => &[ELECTRICAL_CONDUCTIVITY_PROPERTY],
        }
    }
    fn conductivity(&self) -> PlanarConductivity {
        PlanarConductivity::Scalar(self.resistivity_ohm_m().recip())
    }
    fn with_material(&self, material: ResolvedMaterialStatePoint) -> Result<Self, ConductorError> {
        Self::try_from_material_with_basis(material, self.source_basis())
    }
}

/// Immutable source principal values plus their explicit geometry-frame law.
/// Material names do not select anisotropy and missing axes never imply zero.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedOrthotropicConductor2 {
    material: ResolvedMaterialStatePoint,
    basis: ElectricalPropertyBasis,
    frame: PlanarElectricalFrame,
    principal_conductivity_s_m: [f64; 2],
    conductivity_s_m: [[f64; 2]; 2],
    identity: ContentHash,
}
impl ResolvedOrthotropicConductor2 {
    /// Adopt a complete multiphysics source bundle without dropping other roles.
    /// The frame maps both independently measured principal coefficients into
    /// the same geometry coordinates used by E and J. Source axes must describe
    /// fixed external conditions; frequency/field/strain-dependent laws refuse.
    pub fn try_from_material(
        material: ResolvedMaterialStatePoint,
        basis: ElectricalPropertyBasis,
        frame: PlanarElectricalFrame,
    ) -> Result<Self, ConductorError> {
        let names = principal_names(basis);
        let mut principal = [0.0; 2];
        for (axis, name) in names.iter().enumerate() {
            let p = material
                .property(name)
                .ok_or(ConductorError::InvalidInput { quantity: name })?;
            if p.requirement().quantity() != QuantitySpec::dimensional(basis.dimensions())
                || !p.is_static_at_fixed_external_conditions()
                || !p.value_si().is_finite()
                || p.value_si() <= 0.0
            {
                return Err(ConductorError::InvalidInput {
                    quantity: "positive static principal electrical coefficients",
                });
            }
            principal[axis] = match basis {
                ElectricalPropertyBasis::Conductivity => p.value_si(),
                ElectricalPropertyBasis::Resistivity => p.value_si().recip(),
            };
            if !principal[axis].is_finite() || principal[axis] <= 0.0 {
                return Err(ConductorError::InvalidDerived {
                    quantity: "principal conductivity from source",
                });
            }
        }
        // Reuse second-order tensor rotation, with physical tensor shear. The
        // zero normal entry is a computational embedding of a declared 2-D
        // law, not an assertion that a 3-D material is a perfect insulator.
        let q = crate::tensor::rotation([0.0, 0.0, 1.0], frame.angle_rad);
        let rotated = crate::tensor::rotate(&[principal[0], principal[1], 0.0, 0.0, 0.0, 0.0], &q);
        let matrix = [[rotated[0], rotated[3]], [rotated[3], rotated[1]]];
        let scale = principal[0].max(principal[1]);
        let determinant = (matrix[0][0] / scale).mul_add(
            matrix[1][1] / scale,
            -(matrix[0][1] / scale) * (matrix[1][0] / scale),
        );
        if matrix.iter().flatten().any(|v| !v.is_finite())
            || !determinant.is_finite()
            || determinant <= 0.0
        {
            return Err(ConductorError::InvalidDerived {
                quantity: "representably positive rotated electrical tensor",
            });
        }
        let mut hash =
            DomainHasher::new("org.frankensim.fs-material.planar-orthotropic-conductor.v1");
        hash.update(material.identity().as_bytes());
        hash.update(frame.material_frame.as_bytes());
        hash.update(frame.geometry_frame.as_bytes());
        hash.update(&frame.angle_rad.to_bits().to_le_bytes());
        hash.update(&[match basis {
            ElectricalPropertyBasis::Resistivity => 0,
            ElectricalPropertyBasis::Conductivity => 1,
        }]);
        Ok(Self {
            material,
            basis,
            frame,
            principal_conductivity_s_m: principal,
            conductivity_s_m: matrix,
            identity: hash.finalize(),
        })
    }
    pub const fn material(&self) -> &ResolvedMaterialStatePoint {
        &self.material
    }
    pub const fn source_basis(&self) -> ElectricalPropertyBasis {
        self.basis
    }
    pub const fn frame(&self) -> PlanarElectricalFrame {
        self.frame
    }
    pub const fn principal_conductivity_s_m(&self) -> &[f64; 2] {
        &self.principal_conductivity_s_m
    }
    pub const fn conductivity_s_m(&self) -> &[[f64; 2]; 2] {
        &self.conductivity_s_m
    }
    /// Binds source state, independent basis, frame identities and orientation.
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }
}
impl sealed::Sealed for ResolvedOrthotropicConductor2 {}
impl PlanarConductorMaterial for ResolvedOrthotropicConductor2 {
    fn material(&self) -> &ResolvedMaterialStatePoint {
        self.material()
    }
    fn source_basis(&self) -> ElectricalPropertyBasis {
        self.basis
    }
    fn frame(&self) -> Option<PlanarElectricalFrame> {
        Some(self.frame)
    }
    fn source_property_names(&self) -> &'static [&'static str] {
        principal_names(self.basis)
    }
    fn conductivity(&self) -> PlanarConductivity {
        PlanarConductivity::Tensor(self.conductivity_s_m)
    }
    fn with_material(&self, material: ResolvedMaterialStatePoint) -> Result<Self, ConductorError> {
        Self::try_from_material(material, self.basis, self.frame)
    }
}

/// Heterogeneous fields may combine ordinary scalar and oriented conductors
/// without replacing their source receipts or creating separate field solvers.
#[derive(Clone, Debug, PartialEq)]
pub enum ResolvedPlanarConductor {
    Isotropic(ResolvedConductorState),
    Orthotropic(ResolvedOrthotropicConductor2),
}
impl From<ResolvedConductorState> for ResolvedPlanarConductor {
    fn from(value: ResolvedConductorState) -> Self {
        Self::Isotropic(value)
    }
}
impl From<ResolvedOrthotropicConductor2> for ResolvedPlanarConductor {
    fn from(value: ResolvedOrthotropicConductor2) -> Self {
        Self::Orthotropic(value)
    }
}
impl sealed::Sealed for ResolvedPlanarConductor {}
impl PlanarConductorMaterial for ResolvedPlanarConductor {
    fn material(&self) -> &ResolvedMaterialStatePoint {
        match self {
            Self::Isotropic(s) => s.material(),
            Self::Orthotropic(s) => s.material(),
        }
    }
    fn source_basis(&self) -> ElectricalPropertyBasis {
        match self {
            Self::Isotropic(s) => s.source_basis(),
            Self::Orthotropic(s) => s.source_basis(),
        }
    }
    fn frame(&self) -> Option<PlanarElectricalFrame> {
        match self {
            Self::Isotropic(_) => None,
            Self::Orthotropic(s) => Some(s.frame()),
        }
    }
    fn source_property_names(&self) -> &'static [&'static str] {
        match self {
            Self::Isotropic(s) => s.source_property_names(),
            Self::Orthotropic(s) => s.source_property_names(),
        }
    }
    fn conductivity(&self) -> PlanarConductivity {
        match self {
            Self::Isotropic(s) => s.conductivity(),
            Self::Orthotropic(s) => s.conductivity(),
        }
    }
    fn with_material(&self, material: ResolvedMaterialStatePoint) -> Result<Self, ConductorError> {
        match self {
            Self::Isotropic(s) => s.with_material(material).map(Self::Isotropic),
            Self::Orthotropic(s) => s.with_material(material).map(Self::Orthotropic),
        }
    }
}

fn principal_names(basis: ElectricalPropertyBasis) -> &'static [&'static str; 2] {
    match basis {
        ElectricalPropertyBasis::Conductivity => &PRINCIPAL_CONDUCTIVITY_PROPERTIES,
        ElectricalPropertyBasis::Resistivity => &PRINCIPAL_RESISTIVITY_PROPERTIES,
    }
}

/// Resolve both principal axes at one complete material state using the usual
/// source selection/pins; no scalar fallback or inferred normal response.
pub fn resolve_orthotropic_conductor2(
    card: &MaterialCard,
    point: &QueryPoint,
    selection: MaterialPropertySelection,
    basis: ElectricalPropertyBasis,
    frame: PlanarElectricalFrame,
) -> Result<ResolvedOrthotropicConductor2, ConductorError> {
    let requirements = principal_names(basis)
        .iter()
        .map(|name| {
            ScalarPropertyRequirement::try_new(
                *name,
                basis.dimensions(),
                ScalarAdmissibility::StrictlyPositive,
            )
        })
        .collect::<Result<Vec<_>, _>>()?;
    ResolvedOrthotropicConductor2::try_from_material(
        resolve_material_state_point(card, point, &requirements, selection)?,
        basis,
        frame,
    )
}
