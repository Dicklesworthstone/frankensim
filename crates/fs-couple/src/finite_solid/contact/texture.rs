//! A declared material line supplies the physical tangent direction of contact.
//!
//! Exactly one participating surface carries the selected line: the native
//! rigid orientation maps a body direction by R, or the original P1 solid
//! geometry maps a reference direction by its TOTAL deformation gradient F.
//! Its projection into the actual contact plane defines a right-handed unit
//! tangent frame. A vanishing/unresolved projection refuses; no arbitrary world
//! axis replaces the material direction. This query moves no witness or state.
//!
//! This is a single-line orientation model. Stretch and shear may rotate that
//! line, but do not infer new friction coefficients or push forward an ellipse
//! metric. The perpendicular tangent is the current orthogonal complement.
//! Choosing which surface supplies an interface texture is an explicit material
//! model; textures of two surfaces are never averaged. Keep that choice fixed
//! for a retained source identity, or explicitly rebind its history.
//!
//! Physical line evolution and coordinate transport of constitutive history
//! are distinct. `fs_tribo::tangential::transport` owns isometric history changes
//! with a declared common material spin and resolved normal path. Passing zero
//! physical spin merely because a new tangent gauge was selected would retain
//! an old world direction, rather than convecting its material texture. This
//! geometry query deliberately performs no history transport or source update.

use super::mixed::{FiniteSolidContactCounterpart, FiniteSolidMixedContact};
use super::{array, contact_frame, magnitude, unit};
use crate::finite_solid::FiniteSolidBody;
use fs_exec::Cx;
use fs_mbd::{RigidBodyState, Vec3};
use fs_solid::hyper3::HyperTetError;
use fs_tribo::partial_slip::{PartialSlipError, TangentFrame};

/// One explicitly selected participating material surface, with a fixed line.
/// Magnitude is irrelevant; every direction must be finite and nonzero.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum FiniteSolidContactTexture {
    /// The first rigid point's body frame carries this material direction.
    FirstRigid { direction_body: Vec3 },
    /// The second rigid point's body frame carries this material direction.
    /// A solid counterpart refuses this choice.
    SecondRigid { direction_body: Vec3 },
    /// The solid point's ORIGINAL tetrahedron carries this reference direction.
    /// Coordinates are those of the original reference mesh, not an elastic
    /// thermal reference. A rigid counterpart refuses this choice.
    Solid { direction_reference: Vec3 },
}

/// Current physical line and the original tribology owner's contact frame.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteSolidContactTextureFrame {
    /// Normalized carried R a0 or F a0, before projection into the plane.
    pub material_direction_world: Vec3,
    /// Unit normal, tangent and bitangent; longitudinal is the material tangent.
    pub frame: TangentFrame,
    /// Norm of the projected unit line; values near zero are ill-conditioned.
    pub tangential_fraction: f64,
}

impl FiniteSolidContactTextureFrame {
    /// Bind this physical material tangent to an existing mixed-contact geometry.
    /// The original normal, witnesses, source coefficients and law are retained.
    #[must_use]
    pub fn bind<F>(self, mut contact: FiniteSolidMixedContact<F>) -> FiniteSolidMixedContact<F> {
        let t = self.frame.longitudinal();
        contact.tangent_world = Vec3::new(t[0], t[1], t[2]);
        contact
    }
}

/// Geometry, attachment or cancellation refusal; no body/history is modified.
#[derive(Debug)]
pub enum FiniteSolidContactTextureError {
    Invalid(&'static str),
    Solid(HyperTetError),
    Frame(PartialSlipError),
    Cancelled,
}
impl From<HyperTetError> for FiniteSolidContactTextureError {
    fn from(error: HyperTetError) -> Self {
        Self::Solid(error)
    }
}
impl From<PartialSlipError> for FiniteSolidContactTextureError {
    fn from(error: PartialSlipError) -> Self {
        Self::Frame(error)
    }
}
impl core::fmt::Display for FiniteSolidContactTextureError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "material-attached contact texture: {self:?}")
    }
}
impl std::error::Error for FiniteSolidContactTextureError {}

impl FiniteSolidBody<'_> {
    /// Query a material-attached tangent at an original mixed contact.
    ///
    /// The native cached reference gradients evaluate the selected P1 element
    /// from only its four original nodes, including prescribed displacements.
    /// No full-mesh displacement buffer, gradient preparation, constitutive
    /// evaluation or mechanical solve is performed. P1 F is constant inside the
    /// element, so the original material point's weights need no interpolation.
    ///
    /// The existing mechanical owner must still admit positions, gaps, masses,
    /// current source laws and impulse balances. This query admits only the
    /// selected attachment, native total F and contact-frame conditioning.
    pub fn contact_texture_frame<F>(
        &self,
        bodies: &[RigidBodyState],
        contact: &FiniteSolidMixedContact<F>,
        texture: FiniteSolidContactTexture,
        cx: &Cx<'_>,
    ) -> Result<FiniteSolidContactTextureFrame, FiniteSolidContactTextureError> {
        use FiniteSolidContactTextureError as E;
        let poll = || cx.checkpoint().map_err(|_| E::Cancelled);
        poll()?;
        let rigid_direction = |body: usize, direction: Vec3| {
            let state = bodies
                .get(body)
                .ok_or(E::Invalid("texture body is outside native body array"))?;
            let direction = unit(direction).map_err(E::Invalid)?;
            Ok::<_, E>(state.pose().orientation().rotate_body_to_world(direction))
        };
        let carried = match texture {
            FiniteSolidContactTexture::FirstRigid { direction_body } => {
                rigid_direction(contact.point.body, direction_body)?
            }
            FiniteSolidContactTexture::SecondRigid { direction_body } => {
                let FiniteSolidContactCounterpart::Rigid(point) = contact.counterpart else {
                    return Err(E::Invalid(
                        "second rigid texture requires a rigid counterpart",
                    ));
                };
                rigid_direction(point.body, direction_body)?
            }
            FiniteSolidContactTexture::Solid {
                direction_reference,
            } => {
                let FiniteSolidContactCounterpart::Solid(point) = contact.counterpart else {
                    return Err(E::Invalid("solid texture requires a solid counterpart"));
                };
                let a = array(unit(direction_reference).map_err(E::Invalid)?);
                let f = self.prepared.deformation_gradient_with(
                    point.tetrahedron(),
                    |node| {
                        std::array::from_fn(|axis| {
                            let dof = 3 * node + axis;
                            self.free
                                .binary_search(&dof)
                                .map_or(self.fixed[dof], |row| self.state.q[row])
                        })
                    },
                    cx,
                )?;
                let direction: [f64; 3] = std::array::from_fn(|row| {
                    f[3 * row].mul_add(a[0], f[3 * row + 1].mul_add(a[1], f[3 * row + 2] * a[2]))
                });
                Vec3::new(direction[0], direction[1], direction[2])
            }
        };
        let material_direction_world = unit(carried).map_err(E::Invalid)?;
        // Reuse the mechanical owner's robust cross-product construction and
        // its 128-epsilon angular-resolution refusal, rather than normal-part
        // subtraction that can invent a direction for a nearly parallel line.
        let axes =
            contact_frame(contact.normal_world, material_direction_world).map_err(E::Invalid)?;
        let frame = TangentFrame::new(array(axes[0]), array(axes[1]))?;
        let tangential_fraction = magnitude(axes[0].cross(material_direction_world));
        poll()?;
        Ok(FiniteSolidContactTextureFrame {
            material_direction_world,
            frame,
            tangential_fraction,
        })
    }
}
