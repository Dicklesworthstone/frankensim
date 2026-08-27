//! G0/G3 coverage for typed factor and ordered-product manifold layouts.
//!
//! This tranche separates point storage, optimization coordinates, and
//! intrinsic tangent dimension at the public type boundary. Product operations
//! compose the live factor-local validation and retraction contracts; they do
//! not add a robot-specific manifold variant or wire-schema representation.

#![deny(unsafe_code)]

use fs_opt::{
    MANIFOLD_LAYOUT_SCHEMA_VERSION, Manifold, ManifoldLayoutError, OptError,
    PRODUCT_MANIFOLD_LAYOUT_SCHEMA_VERSION, ParamDim, PointDim, ProductCoordinate, ProductFactor,
    ProductFactorId, ProductManifold, ProductManifoldError, TangentDim,
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
        Manifold::Stiefel { n: u32::MAX, p: 1 },
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

fn factor(id: u32, manifold: Manifold) -> ProductFactor {
    ProductFactor::new(ProductFactorId::new(id), manifold)
}

fn bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|value| value.to_bits()).collect()
}

/// G0: a heterogeneous configuration retains declaration order, stable factor
/// identity, checked cumulative dimensions, and three non-interchangeable
/// offset spaces.
#[test]
fn heterogeneous_product_layout_has_exact_typed_offsets() {
    let product = ProductManifold::new(vec![
        factor(10, Manifold::Rn { dim: 3 }),
        factor(20, Manifold::So3),
        factor(30, Manifold::Sphere { ambient: 4 }),
        factor(40, Manifold::Stiefel { n: 3, p: 2 }),
    ])
    .expect("valid heterogeneous configuration");
    let layout = product.layout();

    assert_eq!(
        layout.schema_version(),
        PRODUCT_MANIFOLD_LAYOUT_SCHEMA_VERSION
    );
    assert_eq!(layout.point_dim().get(), 17);
    assert_eq!(layout.param_dim().get(), 16);
    assert_eq!(layout.tangent_dim().get(), 12);

    let expected = [
        (10, Manifold::Rn { dim: 3 }, 0, 0, 0),
        (20, Manifold::So3, 3, 3, 3),
        (30, Manifold::Sphere { ambient: 4 }, 7, 6, 6),
        (40, Manifold::Stiefel { n: 3, p: 2 }, 11, 10, 9),
    ];
    assert_eq!(layout.factors().len(), expected.len());
    for (index, (block, (id, manifold, point, parameter, tangent))) in
        layout.factors().iter().zip(expected).enumerate()
    {
        assert_eq!(block.index(), index as u32);
        assert_eq!(block.factor().id().get(), id);
        assert_eq!(block.factor().manifold(), manifold);
        assert_eq!(block.manifold_layout(), manifold.layout().unwrap());
        assert_eq!(block.point_offset().get(), point);
        assert_eq!(block.param_offset().get(), parameter);
        assert_eq!(block.tangent_offset().get(), tangent);
    }
}

/// G0: block slicing resolves by identity while still requiring one complete
/// payload in the correct coordinate space.
#[test]
fn product_blocks_slice_points_parameters_and_tangents_without_conflation() {
    let translation = ProductFactorId::new(101);
    let orientation = ProductFactorId::new(202);
    let product = ProductManifold::new(vec![
        ProductFactor::new(translation, Manifold::Rn { dim: 3 }),
        ProductFactor::new(orientation, Manifold::So3),
    ])
    .unwrap();
    let point = [10.0, 11.0, 12.0, 20.0, 21.0, 22.0, 23.0];
    let parameter = [30.0, 31.0, 32.0, 40.0, 41.0, 42.0];
    let tangent = [50.0, 51.0, 52.0, 60.0, 61.0, 62.0];

    assert_eq!(
        product.layout().point_block(translation, &point).unwrap(),
        &point[..3]
    );
    assert_eq!(
        product.layout().point_block(orientation, &point).unwrap(),
        &point[3..]
    );
    assert_eq!(
        product
            .layout()
            .parameter_block(orientation, &parameter)
            .unwrap(),
        &parameter[3..]
    );
    assert_eq!(
        product
            .layout()
            .tangent_block(translation, &tangent)
            .unwrap(),
        &tangent[..3]
    );
    assert_eq!(
        product.layout().point_block(translation, &point[..6]),
        Err(ProductManifoldError::PayloadLength {
            coordinate: ProductCoordinate::Point,
            expected: 7,
            got: 6,
        })
    );
    assert_eq!(
        product
            .layout()
            .tangent_block(ProductFactorId::new(999), &tangent),
        Err(ProductManifoldError::UnknownFactor {
            id: ProductFactorId::new(999),
        })
    );
}

/// G0/G3: R3 x SO(3) retraction is exactly the concatenation of the existing
/// factor operations, including the SO(3) three-parameter/four-point split.
#[test]
fn r3_so3_retraction_is_factor_operation_equal_and_replayable() {
    let translation = Manifold::Rn { dim: 3 };
    let orientation = Manifold::So3;
    let product = ProductManifold::new(vec![factor(1, translation), factor(2, orientation)])
        .expect("R3 x SO(3)");
    let point = [1.0, -2.0, 0.5, 1.0, 0.0, 0.0, 0.0];
    let parameter = [0.25, 0.5, -1.0, 0.0, 0.0, std::f64::consts::FRAC_PI_2];

    product.validate_point(&point).expect("valid product point");
    product
        .validate_parameter(&parameter)
        .expect("valid product parameter");
    let translated = translation.retract(&point[..3], &parameter[..3]).unwrap();
    let rotated = orientation.retract(&point[3..], &parameter[3..]).unwrap();
    let mut expected = translated;
    expected.extend_from_slice(&rotated);

    let landed = product.retract(&point, &parameter).unwrap();
    let replay = product.retract(&point, &parameter).unwrap();
    assert_eq!(bits(&landed), bits(&expected));
    assert_eq!(bits(&replay), bits(&landed));
    product
        .validate_point(&landed)
        .expect("landed blockwise point remains on every factor");
}

/// G0/G3: heterogeneous joint-like factors use the same generic blockwise
/// operation path, with exact equality to every live factor implementation.
#[test]
fn heterogeneous_retraction_is_the_exact_ordered_factor_concatenation() {
    let factors = [
        factor(10, Manifold::Rn { dim: 1 }),
        factor(20, Manifold::Sphere { ambient: 3 }),
        factor(30, Manifold::So3),
        factor(40, Manifold::Stiefel { n: 3, p: 2 }),
    ];
    let product = ProductManifold::new(factors.to_vec()).unwrap();
    let point = [
        2.0, // R1
        1.0, 0.0, 0.0, // S2
        1.0, 0.0, 0.0, 0.0, // SO(3)
        1.0, 0.0, 0.0, 0.0, 1.0, 0.0, // Stiefel(3, 2), column-major
    ];
    let parameter = [
        0.5, // R1
        0.0, 0.25, 0.0, // S2 ambient parameter
        0.1, -0.2, 0.3, // SO(3) body parameter
        0.0, 0.0, 0.1, 0.0, 0.0, 0.2, // Stiefel ambient parameter
    ];

    let mut expected = Vec::new();
    let mut point_offset = 0;
    let mut param_offset = 0;
    for factor in factors {
        let manifold = factor.manifold();
        let point_end = point_offset + manifold.point_dim().unwrap() as usize;
        let param_end = param_offset + manifold.param_dim().unwrap() as usize;
        expected.extend_from_slice(
            &manifold
                .retract(
                    &point[point_offset..point_end],
                    &parameter[param_offset..param_end],
                )
                .unwrap(),
        );
        point_offset = point_end;
        param_offset = param_end;
    }

    let landed = product.retract(&point, &parameter).unwrap();
    assert_eq!(bits(&landed), bits(&expected));
    product.validate_point(&landed).unwrap();
}

/// G0: construction refuses ambiguous, invalid, or unrepresentable product
/// layouts before any block can be sliced or operated on.
#[test]
fn product_layout_refuses_empty_duplicate_invalid_and_overflowing_factors() {
    assert_eq!(
        ProductManifold::new(Vec::new()),
        Err(ProductManifoldError::EmptyProduct)
    );
    assert_eq!(
        ProductManifold::new(vec![
            factor(7, Manifold::Rn { dim: 1 }),
            factor(7, Manifold::So3),
        ]),
        Err(ProductManifoldError::DuplicateFactorId {
            id: ProductFactorId::new(7),
            first_index: 0,
            duplicate_index: 1,
        })
    );
    assert_eq!(
        ProductManifold::new(vec![factor(8, Manifold::Rn { dim: 0 })]),
        Err(ProductManifoldError::FactorLayout {
            id: ProductFactorId::new(8),
            index: 0,
            source: ManifoldLayoutError::ZeroEuclideanDimension,
        })
    );
    assert_eq!(
        ProductManifold::new(vec![
            factor(9, Manifold::Rn { dim: u32::MAX }),
            factor(10, Manifold::Rn { dim: 1 }),
        ]),
        Err(ProductManifoldError::DimensionOverflow {
            coordinate: ProductCoordinate::Point,
            id: ProductFactorId::new(10),
            index: 1,
            offset: u32::MAX,
            factor_dim: 1,
        })
    );
}

/// G0: aggregate shape checks and factor-local finite/domain refusals retain
/// deterministic product attribution. The first declared bad block wins.
#[test]
fn product_validation_and_retraction_fail_closed_with_factor_identity() {
    let product = ProductManifold::new(vec![
        factor(11, Manifold::Rn { dim: 3 }),
        factor(22, Manifold::So3),
    ])
    .unwrap();
    let identity_point = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0];
    let zero_parameter = [0.0; 6];

    assert_eq!(
        product.retract(&identity_point[..6], &zero_parameter),
        Err(ProductManifoldError::PayloadLength {
            coordinate: ProductCoordinate::Point,
            expected: 7,
            got: 6,
        })
    );
    assert_eq!(
        product.retract(&identity_point, &zero_parameter[..5]),
        Err(ProductManifoldError::PayloadLength {
            coordinate: ProductCoordinate::Parameter,
            expected: 6,
            got: 5,
        })
    );

    let quiet_nan = f64::from_bits(0x7ff8_0000_0000_1234);
    let bad_first_point = [quiet_nan, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0];
    assert!(matches!(
        product.validate_point(&bad_first_point),
        Err(ProductManifoldError::FactorOperation {
            id,
            index: 0,
            source: OptError::RetractionNonFinite {
                input: "manifold operation point",
                component: 0,
                bits,
            },
        }) if id == ProductFactorId::new(11) && bits == quiet_nan.to_bits()
    ));

    let bad_rotation = [0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0];
    assert!(matches!(
        product.validate_point(&bad_rotation),
        Err(ProductManifoldError::FactorOperation {
            id,
            index: 1,
            source: OptError::RetractionDomain {
                manifold: "SO(3)",
                what: "point must have unit norm before retraction",
                location: None,
                ..
            },
        }) if id == ProductFactorId::new(22)
    ));

    let bad_parameter = [0.0, 0.0, 0.0, 0.0, quiet_nan, 0.0];
    assert!(matches!(
        product.validate_parameter(&bad_parameter),
        Err(ProductManifoldError::FactorOperation {
            id,
            index: 1,
            source: OptError::RetractionNonFinite {
                input: "product manifold parameter",
                component: 1,
                bits,
            },
        }) if id == ProductFactorId::new(22) && bits == quiet_nan.to_bits()
    ));
    assert!(matches!(
        product.retract(&identity_point, &bad_parameter),
        Err(ProductManifoldError::FactorOperation {
            id,
            index: 1,
            source: OptError::RetractionNonFinite {
                input: "retraction step",
                component: 1,
                bits,
            },
        }) if id == ProductFactorId::new(22) && bits == quiet_nan.to_bits()
    ));
}

/// G5: rebuilding the same declared configuration reproduces the complete
/// layout, factor table, offsets, and dimensions exactly.
#[test]
fn product_layout_rebuild_replays_exactly() {
    let factors = vec![
        factor(3, Manifold::Rn { dim: 1 }),
        factor(5, Manifold::Sphere { ambient: 3 }),
        factor(8, Manifold::So3),
        factor(13, Manifold::Stiefel { n: 4, p: 2 }),
    ];
    let first = ProductManifold::new(factors.clone()).unwrap();
    let second = ProductManifold::new(factors).unwrap();
    assert_eq!(first.layout(), second.layout());
}
