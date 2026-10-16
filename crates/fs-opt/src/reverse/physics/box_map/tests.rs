use super::*;
use std::cell::Cell;
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey, TilePool};

#[derive(Debug, Default)]
struct Polynomial { samples: Cell<usize>, products: Cell<usize> }
impl PhysicsModel for Polynomial {
    fn signature(&self) -> PhysicsSignature<'_> {
        PhysicsSignature { study: "polynomial", parameter_count: 2,
            parameter_dims: Dims::NONE, residual_dims: Dims::NONE }
    }
    fn value_gradient(&self, y: &[f64], _: Option<&Cx<'_>>) -> Result<PhysicsSample, PhysicsError> {
        self.samples.set(self.samples.get()+1);
        Ok(PhysicsSample { value: y[0]*y[0]*y[0]+0.5*y[0]*y[1]*y[1]+2.0*y[1],
            gradient: vec![3.0*y[0]*y[0]+0.5*y[1]*y[1], y[0]*y[1]+2.0] })
    }
    fn second_order(&self) -> Option<&dyn PhysicsHessian> { Some(self) }
}
impl PhysicsHessian for Polynomial {
    fn hessian_vector_product(&self, y: &[f64], v: &[f64], _: Option<&Cx<'_>>) -> Result<Vec<f64>, PhysicsError> {
        self.products.set(self.products.get()+1);
        Ok(vec![6.0*y[0]*v[0]+y[1]*v[1], y[1]*v[0]+y[0]*v[1]])
    }
}
fn close(actual: &[f64], expected: &[f64], tolerance: f64) {
    assert_eq!(actual.len(), expected.len());
    for (&a, &e) in actual.iter().zip(expected) {
        assert!((a-e).abs() <= tolerance*(1.0+e.abs()), "{a:e} != {e:e}");
    }
}
fn bounds() -> LogisticBox {
    LogisticBox::new(&[0.1, -0.5], &[2.0, 3.0], Dims::NONE).unwrap()
}

#[test]
fn complete_coordinate_curvature_matches_closed_form_and_gradient_differences() {
    let bounds = bounds(); let native = Polynomial::default(); let model = bounds.bind(&native).unwrap();
    let x = [0.4, -0.7]; let v = [0.3, -0.4]; let w = [-0.2, 0.6];
    let y = bounds.decode(&x, None).unwrap();
    let slope: Vec<_> = y.iter().enumerate().map(|(i,y)|
        (y-bounds.lower[i])*(bounds.upper[i]-y)/bounds.widths[i]).collect();
    let second: Vec<_> = y.iter().enumerate().map(|(i,y)|
        slope[i]*(1.0-2.0*(y-bounds.lower[i])/bounds.widths[i])).collect();
    let physical = native.value_gradient(&y, None).unwrap();
    let scaled: Vec<_> = slope.iter().zip(v).map(|(a,v)| a*v).collect();
    let hp = native.hessian_vector_product(&y, &scaled, None).unwrap();
    let expected: Vec<_> = (0..2).map(|i| slope[i]*hp[i]+second[i]*physical.gradient[i]*v[i]).collect();
    let hv = model.hessian_vector_product(&x, &v, None).unwrap();
    close(&hv, &expected, 1e-13);
    let hw = model.hessian_vector_product(&x, &w, None).unwrap();
    close(&[w[0]*hv[0]+w[1]*hv[1]], &[v[0]*hw[0]+v[1]*hw[1]], 1e-13);
    let eps = 1e-5;
    let plus: Vec<_> = x.iter().zip(v).map(|(x,v)| x+eps*v).collect();
    let minus: Vec<_> = x.iter().zip(v).map(|(x,v)| x-eps*v).collect();
    let gp = model.value_gradient(&plus, None).unwrap().gradient;
    let gm = model.value_gradient(&minus, None).unwrap().gradient;
    close(&hv, &gp.iter().zip(gm).map(|(a,b)| (a-b)/(2.0*eps)).collect::<Vec<_>>(), 1e-8);
    let naive: Vec<_> = slope.iter().zip(hp).map(|(a,h)| a*h).collect();
    assert!(hv.iter().zip(naive).any(|(a,b)| (a-b).abs() > 0.01));
}

#[test]
fn inverse_and_saturation_gates_run_before_the_physical_provider() {
    let bounds = bounds(); let native = Polynomial::default(); let model = bounds.bind(&native).unwrap();
    for x in [[0.0, 0.0], [-20.0, 20.0], [1.3, -2.1]] {
        let y = bounds.decode(&x, None).unwrap();
        close(&bounds.encode(&y, None).unwrap(), &x, 1e-8);
    }
    for bad in [f64::NAN, f64::INFINITY, -1000.0, 1000.0] {
        assert!(model.value_gradient(&[bad, 0.0], None).is_err());
        assert!(model.hessian_vector_product(&[0.0, 0.0], &[bad, f64::NAN], None).is_err());
    }
    assert!(bounds.encode(&[bounds.lower[0], 0.0], None).is_err());
    assert!(model.hessian_vector_product(&[0.0, 0.0], &[1.0], None).is_err());
    assert_eq!(native.samples.get(), 0); assert_eq!(native.products.get(), 0);
    assert_eq!(model.hessian_vector_product(&[0.0, 0.0], &[0.0, 0.0], None).unwrap(), vec![0.0;2]);
    assert_eq!(native.samples.get(), 1); assert_eq!(native.products.get(), 0);
}

#[derive(Debug)]
struct FirstOrder;
impl PhysicsModel for FirstOrder {
    fn signature(&self) -> PhysicsSignature<'_> {
        PhysicsSignature { study: "first", parameter_count: 1,
            parameter_dims: Dims([1,0,0,0,0,0]), residual_dims: Dims([2,0,0,0,0,0]) }
    }
    fn value_gradient(&self, y: &[f64], _: Option<&Cx<'_>>) -> Result<PhysicsSample, PhysicsError> {
        Ok(PhysicsSample { value: y[0]*y[0], gradient: vec![2.0*y[0]] })
    }
}
#[test]
fn dimensional_admission_and_first_order_capability_are_preserved() {
    let physical = FirstOrder;
    let wrong = LogisticBox::new(&[0.1], &[0.9], Dims::NONE).unwrap();
    assert!(wrong.bind(&physical).is_err());
    let bounds = LogisticBox::new(&[0.1], &[0.9], physical.signature().parameter_dims).unwrap();
    let model = bounds.bind(&physical).unwrap();
    assert_eq!(model.signature().parameter_dims, Dims::NONE);
    assert_eq!(model.signature().residual_dims, physical.signature().residual_dims);
    assert_eq!(model.signature().study, "first");
    assert!(model.second_order().is_none());
    assert!(model.hessian_vector_product(&[0.0], &[1.0], None).is_err());
    close(&model.value_gradient(&[0.0], None).unwrap().gradient, &[0.2], 1e-14);
    for (lo, hi) in [(1.0,1.0), (2.0,1.0), (-f64::MAX,f64::MAX), (1.0,f64::from_bits(1.0_f64.to_bits()+1))] {
        assert!(LogisticBox::new(&[lo], &[hi], Dims::NONE).is_err());
    }
}

#[derive(Debug)]
struct Malformed(u8);
impl PhysicsModel for Malformed {
    fn signature(&self) -> PhysicsSignature<'_> {
        PhysicsSignature { study: "bad", parameter_count: 2, parameter_dims: Dims::NONE, residual_dims: Dims::NONE }
    }
    fn value_gradient(&self, _: &[f64], _: Option<&Cx<'_>>) -> Result<PhysicsSample, PhysicsError> {
        if self.0 == 3 { return Err(PhysicsError::Domain("native domain refusal")); }
        Ok(PhysicsSample { value: 1.0, gradient: match self.0 {
            0 => vec![1.0], 1 => vec![f64::NAN, 1.0], _ => vec![1.0;2],
        } })
    }
    fn second_order(&self) -> Option<&dyn PhysicsHessian> { Some(self) }
}
impl PhysicsHessian for Malformed {
    fn hessian_vector_product(&self, _: &[f64], _: &[f64], _: Option<&Cx<'_>>) -> Result<Vec<f64>, PhysicsError> {
        Ok(if self.0 == 2 { vec![1.0] } else { vec![1.0, f64::NAN] })
    }
}
#[test]
fn malformed_native_derivatives_and_zero_direction_domain_errors_are_not_hidden() {
    let bounds = bounds();
    for mode in 0..5 {
        let native = Malformed(mode); let model = bounds.bind(&native).unwrap();
        assert!(model.hessian_vector_product(&[0.0;2], &[1.0;2], None).is_err());
    }
    let native = Malformed(3); let model = bounds.bind(&native).unwrap();
    assert_eq!(model.hessian_vector_product(&[0.0;2], &[0.0;2], None),
        Err(PhysicsError::Domain("native domain refusal")));
}

#[test]
fn explicit_cancellation_precedes_decoding_and_native_work() {
    let pool = TilePool::for_host(1, 73); let gate = CancelGate::new_clock_free();
    pool.arena_pool().scope(|arena| {
        let cx = Cx::new(&gate, arena, StreamKey { seed: 73, kernel_id: 92, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic);
        let bounds = bounds(); let native = Polynomial::default(); let model = bounds.bind(&native).unwrap();
        gate.request();
        assert_eq!(model.value_gradient(&[0.0;2], Some(&cx)), Err(PhysicsError::Cancelled));
        assert_eq!(model.hessian_vector_product(&[0.0;2], &[1.0;2], Some(&cx)), Err(PhysicsError::Cancelled));
        assert_eq!(bounds.encode(&[0.5;2], Some(&cx)), Err(PhysicsError::Cancelled));
        assert_eq!(native.samples.get(), 0);
    });
}
