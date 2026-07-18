//! G0 coverage for the typed manifold-layout authority (7tv.22.3).
//!
//! This tranche separates point storage, optimization coordinates, and
//! intrinsic tangent dimension at the public type boundary. It does not yet
//! migrate retraction, projection, transport, solver, or packing consumers.

#![deny(unsafe_code)]

use fs_opt::{
    MANIFOLD_LAYOUT_SCHEMA_VERSION, Manifold, ManifoldLayoutError, ParamDim, PointDim, TangentDim,
};

fn assert_layout(manifold: Manifold, point: u32, parameter: u32, tangent: u32) {
    let layout = manifold.layout().expect("valid descriptor has a layout");
    assert_eq!(layout.schema_version(), MANIFOLD_LAYOUT_SCHEMA_VERSION);
    assert_eq!(layout.manifold(), manifold);

    let point_dim: PointDim = layout.point_dim();
    let param_dim: ParamDim = layout.param_dim();
    let tangent_dim: TangentDim = layout.tangent_dim();
    assert_eq!(point_dim.get(), point);
    assert_eq!(param_dim.get(), parameter);
    assert_eq!(tangent_dim.get(), tangent);

    assert_eq!(manifold.point_dim(), Some(point));
    assert_eq!(manifold.param_dim(), Some(parameter));
    assert_eq!(manifold.tangent_dim(), Some(tangent));
    assert_eq!(
        manifold.layout(),
        Ok(layout),
        "layout recomputation must be exact"
    );
}

/// G0: each descriptor-domain-valid family exposes the declared storage,
/// retraction-coordinate, and intrinsic tangent dimensions without conflating
/// their Rust types.
#[test]
fn dimension_families_have_one_exact_typed_table() {
    assert_layout(Manifold::Rn { dim: 3 }, 3, 3, 3);
    assert_layout(Manifold::Sphere { ambient: 4 }, 4, 4, 3);
    assert_layout(Manifold::So3, 4, 3, 3);
    assert_layout(Manifold::Stiefel { n: 4, p: 2 }, 8, 8, 5);
    assert_layout(Manifold::Stiefel { n: 1, p: 1 }, 1, 1, 0);
}

/// G0: exact representability boundaries are independent from later
/// deployment caps, and the first overflowing Stiefel storage formula refuses.
#[test]
fn checked_layout_formulas_hold_at_the_u32_boundary() {
    assert_layout(Manifold::Rn { dim: u32::MAX }, u32::MAX, u32::MAX, u32::MAX);
    assert_layout(
        Manifold::Sphere { ambient: u32::MAX },
        u32::MAX,
        u32::MAX,
        u32::MAX - 1,
    );
    assert_layout(
        Manifold::Stiefel {
            n: u32::from(u16::MAX),
            p: u32::from(u16::MAX),
        },
        4_294_836_225,
        4_294_836_225,
        2_147_385_345,
    );

    let n = 65_536;
    let p = 65_536;
    assert_eq!(
        Manifold::Stiefel { n, p }.layout(),
        Err(ManifoldLayoutError::StiefelPointDimensionOverflow { n, p })
    );
}

/// G0: domain-invalid descriptors fail through stable structured variants.
/// The legacy raw formula accessors remain unchanged compatibility projections;
/// only the typed layout confers descriptor-domain validity, while deployment
/// admission remains separately governed by `AdmissionCaps`.
#[test]
fn invalid_raw_descriptors_never_produce_typed_layouts() {
    assert_eq!(
        Manifold::Rn { dim: 0 }.layout(),
        Err(ManifoldLayoutError::ZeroEuclideanDimension)
    );
    for ambient in [0, 1] {
        assert_eq!(
            Manifold::Sphere { ambient }.layout(),
            Err(ManifoldLayoutError::DegenerateSphere { ambient })
        );
    }
    for (n, p) in [(0, 0), (4, 0), (2, 3)] {
        assert_eq!(
            Manifold::Stiefel { n, p }.layout(),
            Err(ManifoldLayoutError::InvalidStiefelFrame { n, p })
        );
    }

    assert_eq!(Manifold::Rn { dim: 0 }.point_dim(), Some(0));
    assert_eq!(Manifold::Sphere { ambient: 1 }.tangent_dim(), Some(0));
    assert_eq!(Manifold::Stiefel { n: 2, p: 3 }.point_dim(), Some(6));
}
