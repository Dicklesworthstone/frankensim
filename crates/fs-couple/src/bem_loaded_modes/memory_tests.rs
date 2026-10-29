use super::*;
use super::super::tests::{oscillator, options};

fn branch(projection: Vec<f64>, stiffness: f64, tau: f64) -> RelaxationBranch {
    RelaxationBranch { projection, stiffness, relaxation_time_s: tau }
}
fn fixture() -> BemLoadedModes {
    let mut model = oscillator(0.0);
    // The dark mechanical coordinate has characteristic polynomial
    // (s+200)(s^2+100s+90000). No radiation parameter is fitted.
    model.stiffness[3] = 60_000.0;
    model.viscosity[3] = 0.0;
    model.settings.max_dimension = 12;
    model.settings.max_operator_evaluations = 300;
    model.with_relaxation_branches(&[branch(vec![0.0,1.0],50_000.0,1.0/300.0)],1,||false).unwrap()
}
fn local_options() -> PoleOptions {
    let mut o = options(); o.max_dimension = 12; o
}

#[test]
fn g1_retained_memory_recovers_oscillatory_and_nonoscillatory_exact_poles() {
    let model = fixture();
    let expected = [C64::new(87_500.0_f64.sqrt(),-50.0),C64::new(0.0,-200.0)];
    for exact in expected {
        let pole = model.resonance_pole(exact+C64::new(0.0,3.0),
            &[C64::ZERO,C64::ONE],&local_options(),||false).unwrap();
        assert!((pole.angular_frequency_rad_s-exact).abs()<1e-8,"{pole:?}");
        assert!(pole.mechanical_residual<1e-10 && pole.boundary_residual<1e-10 && pole.memory_residual<1e-10);
        assert_eq!(pole.elastic_memory_sqrt_j.len(),1);
        assert_eq!(pole.surface_pressure_pa.len(),8);
        let s = C64::new(pole.angular_frequency_rad_s.im,-pole.angular_frequency_rad_s.re);
        let e = pole.elastic_memory_sqrt_j[0];
        let c = 50_000.0_f64.sqrt();
        assert!(((C64::ONE+s.scale(1.0/300.0))*e
            -s.scale(c/300.0)*pole.modal_displacement[1]).abs()<1e-7);
    }
}

#[test]
fn g1_contour_crosses_bare_material_rate_without_meromorphic_count_error() {
    let model = fixture();
    let o = BeynOptions {
        nodes:128,probe_columns:2,rank_relative_tolerance:1e-8,
        count_tolerance:1e-6,quadrature_tolerance:1e-5,
        contour_min_reciprocal_condition:1e-12,linear_residual_tolerance:1e-11,
        svd_tolerance:1e-13,svd_max_sweeps:64,eigen_max_sweeps:128,
        max_dimension:12,max_dense_bytes:1_000_000,max_work:100_000_000,
        refinement:local_options(),
    };
    // Contains BOTH the actual -200i pole and the bare -300i material rate.
    // The latter is not a pole of this coupled system. Eliminating memory
    // would put a meromorphic singularity inside the contour and corrupt count.
    let result = model.resonance_spectrum(Ellipse {
        center:C64::new(0.0,-200.0),real_radius:80.0,imaginary_radius:130.0,
    },&o,||false).unwrap();
    assert_eq!(result.argument_count,1);
    assert_eq!(result.poles.len(),1);
    assert!((result.poles[0].angular_frequency_rad_s-C64::new(0.0,-200.0)).abs()<1e-8);
    let bare = model.original_pencil(C64::new(0.0,-300.0),||false).unwrap();
    assert!(bare.matrix.iter().chain(&bare.derivative).all(|&z| finite(z)));
}

#[test]
fn g0_full_signed_memory_keeps_cross_terms_and_equal_time_hidden_histories() {
    let mut base = oscillator(0.2);
    base.settings.max_dimension = 12;
    let model = base.with_relaxation_branches(&[
        branch(vec![1.0,-2.0],9.0,0.01),branch(vec![2.0,-4.0],9.0,0.01),
    ],2,||false).unwrap();
    assert_eq!(model.memory_count(),2); // no equal-time state merging
    let w = C64::new(0.0,-100.0);
    let a = model.original_pencil(w,||false).unwrap();
    let n = model.dimension();
    // x=p=0, e=(2,-1): zero net mechanical force, but genuine material
    // storage and relaxation. Both states must remain even though it is dark.
    for i in 0..n {
        let residual = a.matrix[i*n+10].scale(2.0)-a.matrix[i*n+11];
        assert!(residual.abs()<1e-12,"row={i}: {residual:?}");
    }
    // Every cross term of c c^T remains in the independent Schur image.
    let w = C64::new(350.0,-7.0);
    let a = model.original_pencil(w,||false).unwrap();
    let s = C64::new(w.im,-w.re);
    let factor = s.scale(0.01)/(C64::ONE+s.scale(0.01));
    for i in 0..2 { for j in 0..2 {
        let correction = -(a.matrix[i*n+10]*a.matrix[10*n+j]/a.matrix[10*n+10]
            +a.matrix[i*n+11]*a.matrix[11*n+j]/a.matrix[11*n+11]);
        let c = [3.0,-6.0];
        assert!((correction-factor.scale(5.0*c[i]*c[j])).abs()<1e-10);
    }}
}

#[test]
fn g1_memory_analytic_derivative_includes_each_original_internal_equation() {
    let model = fixture();
    for w in [C64::new(400.0,-30.0),C64::new(0.0,-300.0)] {
        let a = model.original_pencil(w,||false).unwrap();
        for direction in [C64::ONE,C64::new(0.0,1.0)] {
            let h = direction.scale(0.001);
            let plus = model.original_pencil(w+h,||false).unwrap();
            let minus = model.original_pencil(w-h,||false).unwrap();
            for i in 0..a.matrix.len() {
                let fd = (plus.matrix[i]-minus.matrix[i])/h.scale(2.0);
                assert!((fd-a.derivative[i]).abs()<3e-7*(1.0+a.derivative[i].abs()));
            }
        }
    }
}

#[test]
fn g0_memory_admission_counts_all_states_and_empty_attachment_preserves_bits() {
    let model = oscillator(0.2);
    let w = C64::new(350.0,-3.0);
    let before = model.original_pencil(w,||false).unwrap();
    let empty = model.with_relaxation_branches(&[],0,||false).unwrap();
    let after = empty.original_pencil(w,||false).unwrap();
    assert_eq!(before.matrix,after.matrix);
    assert_eq!(before.derivative,after.derivative);
    // Existing exact ten-state cap refuses an eleventh material coordinate.
    assert!(oscillator(0.2).with_relaxation_branches(&[branch(vec![1.0,0.0],1.0,0.1)],1,||false).is_err());
    let model = fixture();
    assert!(model.with_relaxation_branches(&[branch(vec![1.0,0.0],1.0,0.1)],1,||false).is_err());
    let result = oscillator(0.2).with_relaxation_branches(&[],0,||true);
    assert!(matches!(result,Err(BemModesError::Cancelled)));
    for arm in [branch(vec![1.0],1.0,0.1),branch(vec![0.0,0.0],1.0,0.1),
        branch(vec![1.0,0.0],-1.0,0.1),branch(vec![1.0,0.0],1.0,0.0)] {
        let mut model = oscillator(0.0);model.settings.max_dimension=12;
        assert!(model.with_relaxation_branches(&[arm],1,||false).is_err());
    }
    let model=fixture();
    assert!(matches!(model.resonance_pole(C64::new(0.0,-300.0),
        &[C64::ZERO,C64::ONE],&local_options(),||false),Err(BemModesError::SeedMemorySingular{branch:0})));
}
