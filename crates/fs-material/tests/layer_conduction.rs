//! Independent slab balances and derivatives, not iteration telemetry oracles.
use fs_material::thermal::conduction::{
    ConductivityCurve, LayerConductionConfig, TwoLayerConduction, MAX_CONDUCTIVITY_KNOTS,
};
use fs_thermal::enthalpy::{DepositionError, series_contact_conductance_w_k};

fn config() -> LayerConductionConfig {
    LayerConductionConfig { max_iterations: 96, absolute_heat_tolerance_w: 1e-9,
        relative_heat_tolerance: 1e-12 }
}
fn constant(k: f64) -> ConductivityCurve {
    ConductivityCurve::new(vec![[200.0, k], [1000.0, k]]).unwrap()
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a-b).abs() <= tolerance, "{a:.16e} != {b:.16e}");
}

#[test]
fn constant_layers_recover_series_bulk_and_contact_resistance_in_both_directions() {
    let a = constant(8.0);
    let b = constant(2.0);
    for r in [0.0, 0.01, 0.1] {
        let layer = TwoLayerConduction::new([&a, &b], [0.02, 0.03], 0.004, r).unwrap();
        let g = series_contact_conductance_w_k(0.004, [0.02, 0.03], [8.0, 2.0], r).unwrap();
        for t in [[600.0, 300.0], [300.0, 600.0], [400.0; 2]] {
            let result = layer.evaluate(t, config(), || Ok(())).unwrap();
            close(result.power_a_to_b_w, g * (t[0]-t[1]), 1e-9);
            close(result.conductance_w_k, g, 1e-11);
            close(result.endpoint_derivative_w_k[0], g, 1e-12);
            close(result.endpoint_derivative_w_k[1], -g, 1e-12);
            let j = result.power_a_to_b_w / 0.004;
            close(result.interface_temperature_k[0], t[0] - j * 0.02/8.0, 1e-9);
            close(result.interface_temperature_k[1], t[1] + j * 0.03/2.0, 1e-9);
            close(result.interface_temperature_k[0] - result.interface_temperature_k[1], j*r, 1e-9);
        }
    }
}

#[test]
fn large_temperature_gradient_integrates_k_instead_of_freezing_cell_endpoints() {
    // k(T)=2+0.01*T. The exact integrated heat for these lengths is 156 W.
    let k = ConductivityCurve::new(vec![[200.0, 4.0], [1000.0, 12.0]]).unwrap();
    let model = TwoLayerConduction::new([&k, &k], [0.02, 0.03], 0.004, 0.0).unwrap();
    let result = model.evaluate([600.0, 300.0], config(), || Ok(())).unwrap();
    let exact = 0.004/0.05 * (2.0 * 300.0 + 0.005 * (600.0*600.0-300.0*300.0));
    close(result.power_a_to_b_w, exact, 1e-9);
    let interface = ((48.4_f64).sqrt()-2.0)/0.01;
    close(result.interface_temperature_k[0], interface, 1e-9);
    close(result.interface_temperature_k[1], interface, 1e-9);
    close(result.endpoint_derivative_w_k[0], 0.004/0.05 * 8.0, 1e-10);
    close(result.endpoint_derivative_w_k[1], -0.004/0.05 * 5.0, 1e-10);
    let frozen = series_contact_conductance_w_k(0.004, [0.02, 0.03], [8.0, 5.0], 0.0).unwrap()*300.0;
    assert!((frozen-exact).abs() > 14.0, "endpoint-only k is a different constitutive approximation");
}

#[test]
fn piecewise_integrals_retain_breakpoints_small_differences_and_exact_reversal() {
    let curve = ConductivityCurve::new(vec![[200.0, 1.0], [400.0, 5.0], [600.0, 3.0], [900.0, 9.0]]).unwrap();
    close(curve.integral_w_m(200.0, 900.0).unwrap(), 600.0+800.0+1800.0, 1e-12);
    let seam = 400.0_f64;
    let a = seam.next_down();
    let b = seam.next_up();
    let tiny = curve.integral_w_m(a,b).unwrap();
    assert!(tiny > 0.0);
    close(tiny/(5.0*(b-a)), 1.0, 2e-15);
    assert_eq!(curve.integral_w_m(b,a), Ok(-tiny));
    assert_eq!(curve.integral_w_m(seam,seam), Ok(0.0));
    close(curve.integral_w_m(270.0, 850.0).unwrap(),
        curve.integral_w_m(270.0, 450.0).unwrap()+curve.integral_w_m(450.0,850.0).unwrap(), 1e-12);
}

#[test]
fn unequal_nonlinear_layers_have_reciprocal_heat_and_complete_endpoint_derivatives() {
    let a = ConductivityCurve::new(vec![[200.0, 2.0], [450.0, 8.0], [700.0, 5.0], [1000.0, 10.0]]).unwrap();
    let b = ConductivityCurve::new(vec![[200.0, 12.0], [500.0, 3.0], [1000.0, 7.0]]).unwrap();
    for r in [0.0, 0.02] {
        for t in [[650.0, 310.0], [320.0, 680.0]] {
            let model = TwoLayerConduction::new([&a,&b], [0.02,0.035], 0.004,r).unwrap();
            let result = model.evaluate(t,config(),||Ok(())).unwrap();
            let qi = [0.004/0.02 * a.integral_w_m(result.interface_temperature_k[0],t[0]).unwrap(),
                0.004/0.035 * b.integral_w_m(t[1],result.interface_temperature_k[1]).unwrap()];
            for q in qi { close(q, result.power_a_to_b_w, 2e-9); }
            assert!(result.power_a_to_b_w*(1.0/t[1]-1.0/t[0]) > 0.0);
            let swapped = TwoLayerConduction::new([&b,&a], [0.035,0.02], 0.004,r).unwrap()
                .evaluate([t[1],t[0]],config(),||Ok(())).unwrap();
            assert_eq!(result.power_a_to_b_w,-swapped.power_a_to_b_w);
            assert_eq!(result.interface_temperature_k, [swapped.interface_temperature_k[1],swapped.interface_temperature_k[0]]);
            for i in 0..2 {
                let mut plus=t; let mut minus=t;
                plus[i]+=0.001; minus[i]-=0.001;
                let fd=(model.evaluate(plus,config(),||Ok(())).unwrap().power_a_to_b_w
                    -model.evaluate(minus,config(),||Ok(())).unwrap().power_a_to_b_w)/0.002;
                close(result.endpoint_derivative_w_k[i],fd,1e-6);
            }
            let scaled = TwoLayerConduction::new([&a,&b], [0.02,0.035], 0.008,r).unwrap()
                .evaluate(t,config(),||Ok(())).unwrap();
            close(scaled.power_a_to_b_w,2.0*result.power_a_to_b_w,3e-9);
        }
    }
}

#[test]
fn unsupported_sources_work_exhaustion_and_final_cancellation_return_no_response() {
    for knots in [vec![], vec![[200.0,1.0]], vec![[200.0,1.0],[200.0,2.0]],
        vec![[200.0,1.0],[400.0,0.0]], vec![[200.0,1.0],[400.0,f64::NAN]],
        vec![[200.0,1.0];MAX_CONDUCTIVITY_KNOTS+1]] {
        assert!(ConductivityCurve::new(knots).is_err());
    }
    let curve = ConductivityCurve::new(vec![[200.0,4.0],[1000.0,12.0]]).unwrap();
    for bad in [199.0,1001.0,f64::NAN,f64::INFINITY] {
        assert!(curve.conductivity_w_m_k(bad).is_err());
        assert!(curve.integral_w_m(bad,bad).is_err());
    }
    let model = TwoLayerConduction::new([&curve;2],[0.02,0.03],0.004,0.0).unwrap();
    assert!(matches!(model.evaluate([600.0,300.0], LayerConductionConfig {
        max_iterations:1,..config()
    },||Ok(())),Err(DepositionError::Budget(_))));
    let mut calls=0;
    let expected=model.evaluate([600.0,300.0],config(),||{calls+=1;Ok(())}).unwrap();
    let mut retry_calls=0;
    assert_eq!(model.evaluate([600.0,300.0],config(),||{
        retry_calls+=1;
        if retry_calls==calls {Err(DepositionError::Cancelled)} else {Ok(())}
    }),Err(DepositionError::Cancelled));
    assert_eq!(model.evaluate([600.0,300.0],config(),||Ok(())).unwrap(),expected);
}

#[test]
fn near_perfect_contact_does_not_amplify_kelvin_roundoff_into_false_energy_loss() {
    let k = ConductivityCurve::new(vec![[200.0,4.0],[1000.0,12.0]]).unwrap();
    let perfect = TwoLayerConduction::new([&k;2],[0.02,0.03],0.004,0.0).unwrap()
        .evaluate([600.0,300.0],config(),||Ok(())).unwrap();
    for r in [1e-16,1e-14,1e-12,1e-10,1e-8] {
        let got = TwoLayerConduction::new([&k;2],[0.02,0.03],0.004,r).unwrap()
            .evaluate([600.0,300.0],config(),||Ok(())).unwrap();
        assert!(got.layer_heat_residual_w.abs() <= 1e-9+1e-12*got.power_a_to_b_w.abs());
        assert!(got.contact_temperature_residual_k.abs() <= got.contact_roundoff_allowance_k
            +(1e-9+1e-12*got.power_a_to_b_w.abs())/0.004*r);
        assert!(got.power_a_to_b_w <= perfect.power_a_to_b_w+1e-9);
        close(got.power_a_to_b_w,perfect.power_a_to_b_w,30000.0*r+1e-9);
    }
}
