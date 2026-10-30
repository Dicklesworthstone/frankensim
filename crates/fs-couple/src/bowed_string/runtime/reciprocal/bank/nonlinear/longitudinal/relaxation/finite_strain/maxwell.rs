//! One Maxwell strain coordinate, eliminated inside the existing joint trial.
//! Returns relaxed strain, work-conjugate midpoint stress coefficient, and heat.
pub(super) fn step(before:f64,after:f64,memory:f64,rate:f64,stiffness:f64)
    -> Option<(f64,f64,f64)> {
    if [before,after,memory,rate,stiffness].iter().any(|v|!v.is_finite())
        || rate<=0.0 || rate>2.0 || stiffness<=0.0 {return None;}
    let residual=(f64::midpoint(before,after)-memory)/(1.0+0.5*rate);
    let next=memory+rate*residual;
    let force=stiffness*residual;
    let loss=stiffness*rate*residual*residual;
    [next,force,loss].iter().all(|v|v.is_finite()).then_some((next,force,loss))
}

#[cfg(test)]
mod tests {
    use super::step;
    #[test]
    fn g0_step_closes_storage_work_and_positive_heat() {
        for rate in [1e-8,0.001,0.2,1.0,2.0] {for sign in [-1.0,1.0] {
            let (a,b,z,c)=(sign*0.001,sign*0.003,sign*0.0007,2500.0);
            let (next,stress,heat)=step(a,b,z,rate,c).unwrap();
            let initial=0.5*c*(a-z).powi(2);let final_energy=0.5*c*(b-next).powi(2);
            assert!(heat>=0.0);
            assert!((final_energy-initial+heat-stress*(b-a)).abs()<1e-14);
        }}
    }
    #[test]
    fn g1_fixed_strain_relaxation_has_the_original_time_constant() {
        let (epsilon,z0,c,total,tau)=(0.002,0.0003,1700.0,0.002,0.001);
        let mut errors=Vec::new();
        for steps in [16,32,64] {
            let mut z=z0;let rate=total/(steps as f64*tau);let mut heat=0.0;
            for _ in 0..steps {let next=step(epsilon,epsilon,z,rate,c).unwrap();z=next.0;heat+=next.2;}
            let expected=epsilon+(z0-epsilon)*(-total/tau).exp();
            errors.push((z-expected).abs());
            assert!((0.5*c*(epsilon-z).powi(2)+heat-0.5*c*(epsilon-z0).powi(2)).abs()<1e-14);
        }
        for pair in errors.windows(2) {assert!((pair[0]/pair[1]-4.0).abs()<0.02);}
    }
    #[test]
    fn g1_endpoint_strain_tangent_is_the_implicit_maxwell_modulus() {
        let (a,b,z,rate,c)=(0.001,0.003,0.0007,0.25,2500.0);let h=1e-7;
        let plus=step(a,b+h,z,rate,c).unwrap();let minus=step(a,b-h,z,rate,c).unwrap();
        let gain=1.0/(2.0+rate);
        assert!(((plus.1-minus.1)/(2.0*h)-c*gain).abs()<1e-8);
        assert!(((plus.0-minus.0)/(2.0*h)-rate*gain).abs()<1e-10);
    }
    #[test]
    fn g0_relaxed_state_is_stationary_and_invalid_steps_refuse() {
        assert_eq!(step(0.002,0.002,0.002,0.2,1500.0),Some((0.002,0.0,0.0)));
        for rate in [0.0,-1.0,2.01,f64::INFINITY] {assert!(step(0.0,0.0,0.0,rate,1.0).is_none());}
        assert!(step(0.0,0.0,0.0,0.2,-1.0).is_none());
        assert!(step(f64::NAN,0.0,0.0,0.2,1.0).is_none());
    }
}
