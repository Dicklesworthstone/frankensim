//! Product geometry is the concatenation of the actual single-factor geometry.
use fs_opt::{
    Manifold, ProductCoordinate, ProductDifferentialError, ProductFactor,
    ProductFactorId, ProductManifold, ProductManifoldError,
};

fn fixture() -> (ProductManifold, Vec<f64>, Vec<f64>) {
    let factors = vec![
        ProductFactor::new(ProductFactorId::new(99), Manifold::Rn { dim: 2 }),
        ProductFactor::new(ProductFactorId::new(7), Manifold::Sphere { ambient: 3 }),
        ProductFactor::new(ProductFactorId::new(42), Manifold::So3),
        ProductFactor::new(ProductFactorId::new(3), Manifold::Stiefel { n: 3, p: 2 }),
    ];
    let product = ProductManifold::new(factors).unwrap();
    let x = vec![0.2, -0.4, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
        1.0, 0.0, 0.0, 0.0, 1.0, 0.0];
    let ambient: Vec<f64> = (0..15).map(|i| (i as f64 - 6.0) * 0.07).collect();
    (product, x, ambient)
}

fn bits(x: &[f64]) -> Vec<u64> { x.iter().map(|x| x.to_bits()).collect() }
fn dot(x: &[f64], y: &[f64]) -> f64 { x.iter().zip(y).map(|(x, y)| x*y).sum() }

#[test]
fn differential_blocks_match_factor_authority_bit_for_bit() {
    let (product, x, ambient) = fixture();
    let gradient = product.parameter_gradient(&x, &ambient).unwrap();
    assert_eq!(x.len(), 15);
    assert_eq!(gradient.len(), 14);
    product.validate_parameter_tangent(&x, &gradient).unwrap();
    let curve = product.retract_curve(&x, &gradient, 0.2).unwrap();
    let step: Vec<f64> = gradient.iter().map(|g| 0.2*g).collect();
    assert_eq!(bits(&curve.point), bits(&product.retract(&x, &step).unwrap()));
    for block in product.layout().factors() {
        let factor = block.factor();
        let id = factor.id();
        let p = product.layout().point_block(id, &x).unwrap();
        let a = product.layout().point_block(id, &ambient).unwrap();
        let g = product.layout().parameter_block(id, &gradient).unwrap();
        assert_eq!(bits(g), bits(&factor.manifold().parameter_gradient(p, a).unwrap()));
        let local = factor.manifold().retract_curve(p, g, 0.2).unwrap();
        assert_eq!(bits(&local.point), bits(product.layout().point_block(id, &curve.point).unwrap()));
        assert_eq!(bits(&local.velocity), bits(product.layout().parameter_block(id, &curve.velocity).unwrap()));
    }
}

#[test]
fn curve_gradient_pairing_matches_scalar_finite_differences() {
    let (product, x, ambient) = fixture();
    let d = product.parameter_gradient(&x, &ambient).unwrap();
    let a = 0.3;
    let curve = product.retract_curve(&x, &d, a).unwrap();
    let g = product.parameter_gradient(&curve.point, &ambient).unwrap();
    let h = 1e-6;
    let plus = product.retract_curve(&x, &d, a+h).unwrap();
    let minus = product.retract_curve(&x, &d, a-h).unwrap();
    let fd = (dot(&ambient, &plus.point)-dot(&ambient, &minus.point))/(2.0*h);
    assert!((fd-dot(&g, &curve.velocity)).abs() < 2e-9);
}

#[test]
fn transport_matches_each_factor_and_lands_in_product_tangent() {
    let (product, x, ambient) = fixture();
    let v = product.parameter_gradient(&x, &ambient).unwrap();
    let step: Vec<f64> = v.iter().map(|v| 0.3*v).collect();
    let to = product.retract(&x, &step).unwrap();
    let transported = product.transport_parameter(&x, &step, &to, &v).unwrap();
    product.validate_parameter_tangent(&to, &transported).unwrap();
    for block in product.layout().factors() {
        let f = block.factor();
        let layout = product.layout();
        let local = f.manifold().transport_parameter(
            layout.point_block(f.id(), &x).unwrap(),
            layout.parameter_block(f.id(), &step).unwrap(),
            layout.point_block(f.id(), &to).unwrap(),
            layout.parameter_block(f.id(), &v).unwrap(),
        ).unwrap();
        assert_eq!(bits(&local), bits(layout.parameter_block(f.id(), &transported).unwrap()));
    }
}

#[test]
fn aggregate_shape_is_checked_before_any_factor_math() {
    let (product, mut x, ambient) = fixture();
    x[2] = 5.0; // Also invalid, but complete envelope failure has priority.
    let bad = &ambient[..14];
    assert!(matches!(product.parameter_gradient(&x, bad),
        Err(ProductDifferentialError::Geometry(ProductManifoldError::PayloadLength {
            coordinate: ProductCoordinate::Point, expected: 15, got: 14,
        }))));
    assert!(matches!(product.retract_curve(&x, &[0.0; 13], 0.1),
        Err(ProductDifferentialError::Geometry(ProductManifoldError::PayloadLength {
            coordinate: ProductCoordinate::Parameter, expected: 14, got: 13,
        }))));
}

#[test]
fn failures_retain_nonsequential_factor_identity_and_declaration_index() {
    let (product, mut x, ambient) = fixture();
    x[2] = 2.0;
    assert!(matches!(product.parameter_gradient(&x, &ambient),
        Err(ProductDifferentialError::Geometry(ProductManifoldError::FactorOperation {
            id, index: 1, ..
        })) if id == ProductFactorId::new(7)));
    let (product, x, ambient) = fixture();
    let v = product.parameter_gradient(&x, &ambient).unwrap();
    let step: Vec<f64> = v.iter().map(|v| 0.1*v).collect();
    let mut to = product.retract(&x, &step).unwrap();
    to[9] += 1e-4;
    let before = bits(&to);
    assert!(matches!(product.transport_parameter(&x, &step, &to, &v),
        Err(ProductDifferentialError::Geometry(ProductManifoldError::FactorOperation {
            id, index: 3, ..
        })) if id == ProductFactorId::new(3)));
    assert_eq!(bits(&to), before);
}
