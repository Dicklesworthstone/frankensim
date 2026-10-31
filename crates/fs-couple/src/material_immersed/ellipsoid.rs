//! Ellipsoidal source geometry on the shared retained-fluid carrier owner.

use fs_ga::{Mat3, Se3};
use fs_lbm::d3q19::BodyFrameEllipsoid3;

use super::*;

/// One uniform ellipsoidal material region, centered at the carrier COM.
/// The shape may be rotated relative to the authored material/body frame;
/// both its sampled surface and complete inertia use this SAME rotation.
/// Matching integral moments is a consistency check, not a geometry proof.
#[derive(Clone, Copy, Debug)]
pub struct AdditionalMaterialEllipsoid<'a> {
    material: &'a ResolvedBody,
    semiaxes_m: [f64; 3],
    body_from_shape: So3,
}

impl<'a> AdditionalMaterialEllipsoid<'a> {
    /// Bind the original material inertia to matching centered ellipsoid
    /// geometry. `body_from_shape` rotates its principal shape axes into the
    /// material/body coordinates used by [`BodyStateSi`].
    pub fn new(
        material: &'a ResolvedBody,
        semiaxes_m: [f64; 3],
        body_from_shape: So3,
    ) -> Result<Self, Error> {
        let pose =
            Se3::from_parts(body_from_shape, Vec3::new(0.0, 0.0, 0.0)).map_err(|source| {
                Error::Native(PartialSaturationError3::Motion {
                    source: source.into(),
                })
            })?;
        let geometry = VolumeMoments::ellipsoid(vector(semiaxes_m))?.transformed(pose)?;
        if material.regions().len() != 1 || material.regions()[0].geometry() != geometry {
            return Err(Error::Invalid(
                "additional material must be one matching centered uniform ellipsoidal region",
            ));
        }
        Ok(Self {
            material,
            semiaxes_m,
            body_from_shape,
        })
    }

    /// Original density/query receipts and complete material-frame inertia.
    pub const fn material(self) -> &'a ResolvedBody {
        self.material
    }

    /// Principal semiaxes shared by the material volume and sampled surface [m].
    pub const fn semiaxes_m(self) -> [f64; 3] {
        self.semiaxes_m
    }

    /// Fixed proper rotation from principal shape axes to material/body axes.
    pub const fn body_from_shape(self) -> So3 {
        self.body_from_shape
    }
}

impl AdditionalMaterialEllipsoid<'_> {
    pub(super) fn lattice_shape(self, cell_size_m: f64) -> Result<NativeEllipsoid3, Error> {
        Ok(NativeEllipsoid3 {
            shape: BodyFrameEllipsoid3::new(vector(convert(self.semiaxes_m, 1.0 / cell_size_m)?))?,
            shape_from_body: self.body_from_shape.inverse().matrix(),
        })
    }
}

/// Material-derived ellipsoids with resolved native fluid torque and motion.
///
/// This uses retained interior fluid PLUS the additional sourced carrier mass
/// and full tensor. Native midpoint geometry, reciprocal fluid/body exchange,
/// and free rigid rotation share one atomic step. Separated enclosing spheres
/// are required by the native motion owner; near-contact lubrication and a
/// displaced-fluid rigid-solid closure are outside this model.
pub type RetainedFluidEllipsoids<'a, F, const N: usize> =
    RetainedFluidBodies<'a, F, N, AdditionalMaterialEllipsoid<'a>>;

impl<'a, F: FluidTransport + ?Sized, const N: usize> RetainedFluidEllipsoids<'a, F, N> {
    /// Construct quiescent fluid and source-derived, freely rotating ellipsoids.
    /// The canonical material density fixes mass and every inertia component;
    /// neither a representative radius nor a diagonal-only tensor is used.
    pub fn new(
        fluid: &'a F,
        domain: ImmersedDomain,
        sources: [AdditionalMaterialEllipsoid<'a>; N],
        initial: [BodyStateSi; N],
        gate: &CancelGate,
    ) -> Result<Self, Error> {
        Self::new_with_surfaces(
            fluid,
            domain,
            sources,
            sources.map(|s| s.material),
            initial,
            |scaling| {
                let mut shapes = Vec::with_capacity(N);
                for source in sources {
                    poll(gate)?;
                    shapes.push(source.lattice_shape(scaling.cell_size_m())?);
                }
                Ok(CarrierSurfaces::Ellipsoids(
                    shapes
                        .try_into()
                        .map_err(|_| Error::Invalid("ellipsoid count"))?,
                ))
            },
            gate,
        )
    }
}

/// Precompute the fixed frame transform once; membership retains the native
/// ellipsoid predicate without allocation or group validation per sample.
#[derive(Clone, Copy)]
pub(super) struct NativeEllipsoid3 {
    shape: BodyFrameEllipsoid3,
    shape_from_body: Mat3,
}

impl BodyFrameShape3 for NativeEllipsoid3 {
    fn bounding_radius(&self) -> f64 {
        self.shape.bounding_radius()
    }

    fn contains(&self, point: Vec3) -> bool {
        self.shape.contains(self.shape_from_body.apply(point))
    }
}
