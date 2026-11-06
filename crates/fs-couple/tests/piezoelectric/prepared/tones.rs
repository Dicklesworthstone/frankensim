//! G1/G4 continuous source work, charge, resonance and phase ownership.
//! These authored material fixtures are analytic/numerical references.
use super::{controls, near, rod};
use fs_couple::piezoelectric::{PiezoelectricBody, PiezoelectricCouplingError};
use fs_couple::piezoelectric::prepared::{
    PiezoelectricTone, PreparedPiezoelectricTones,
};
use fs_phs::linear_exact::{DenseExponentialBudget, LinearSplitBudget};

fn budgets() -> (LinearSplitBudget, DenseExponentialBudget) {
    (
        LinearSplitBudget { max_modes: 64, max_omega_dt: 3.0 },
        DenseExponentialBudget {
            max_dimension: 128, max_dense_bytes: 32_000_000,
            max_multiply_terms: 1_000_000_000,
        },
    )
}
fn runtime(body: PiezoelectricBody, tones: &[PiezoelectricTone], h: f64)
    -> PreparedPiezoelectricTones {
    let (split,dense) = budgets();
    PreparedPiezoelectricTones::new(body,tones,h,split,dense,|| false).unwrap()
}
fn dc(force: f64, current: f64) -> PiezoelectricTone {
    PiezoelectricTone {
        angular_frequency_rad_s: 0.0,
        force_cosine_n: force, force_sine_n: 0.0,
        current_cosine_a: current, current_sine_a: 0.0,
    }
}

// Closed resonant solution in the original physical coordinates. The source
// is F=fc*cos(wt)+fs*sin(wt), I=ic*cos(wt), with zero initial x,v,Q.
// The work formulas integrate F*v and I*V separately, without an energy remainder.
fn resonant(t: f64, d: f64) -> [f64; 5] {
    let (fc,fs,ic) = (0.004,0.0015,0.002);
    let c = 0.25*(2.0-10.0*d*d);
    let theta = 2.5*d;
    let w = fs_math::det::sqrt(2.5+theta*theta/c); // m = 1
    let z = w*t;
    let sn = fs_math::det::sin(z);
    let cs = fs_math::det::cos(z);
    let sin2 = fs_math::det::sin(2.0*z);
    let cos2 = fs_math::det::cos(2.0*z);
    let se = fs+theta*ic/(c*w);
    let x = fc*t*sn/(2.0*w)+se*(sn-z*cs)/(2.0*w*w);
    let v = fc*(sn+z*cs)/(2.0*w)+se*t*sn/2.0;
    let q = ic*sn/w;
    let a = -t*cos2/(4.0*w)+sin2/(8.0*w*w); // integral t*sin*cos
    let b = sn*sn/(2.0*w);                  // integral sin*cos
    let dd = t*t/4.0+t*sin2/(4.0*w)+(cos2-1.0)/(8.0*w*w);
    let ee = t/2.0-sin2/(4.0*w);
    let ff = t*t/4.0-t*sin2/(4.0*w)-(cos2-1.0)/(8.0*w*w);
    let mechanical = fc*fc*b/(2.0*w)+fc*fc*dd/2.0+fc*se*a/2.0
        + fs*fc*ee/(2.0*w)+fs*fc*a/2.0+fs*se*ff/2.0;
    let electrical = q*q/(2.0*c)-theta*ic/c
        *(fc*a/(2.0*w)+se*(b-w*dd)/(2.0*w*w));
    [x,v,q,mechanical,electrical]
}

#[test]
fn g1_exact_resonance_retains_zero_and_weak_coupling_and_both_source_work_integrals() {
    for d in [0.0,1e-8,0.001,-0.001] {
        let c = 0.25*(2.0-10.0*d*d);
        let theta = 2.5*d;
        let frequency = fs_math::det::sqrt(2.5+theta*theta/c);
        let tones = [PiezoelectricTone {
            angular_frequency_rad_s: frequency,
            force_cosine_n: 0.004, force_sine_n: 0.0015,
            current_cosine_a: 0.002, current_sine_a: 0.0,
        }];
        let body = PiezoelectricBody::new(rod(d,3.0),0.0,0.0,0.0).unwrap();
        let source = body.rod().material().clone();
        let mut prepared = runtime(body,&tones,0.05);
        let mut previous = resonant(0.0,d);
        for interval in 1..=64 {
            let time = interval as f64*0.05;
            let expected = resonant(time,d);
            let report = prepared.step(controls(),|| false).unwrap();
            near(prepared.body().extension_m(),expected[0],2e-12);
            near(prepared.body().velocity_m_s(),expected[1],2e-12);
            near(prepared.body().charge_c(),expected[2],2e-12);
            near(report.mechanical_supply_j,expected[3]-previous[3],2e-13);
            near(report.electrical_supply_j,expected[4]-previous[4],2e-13);
            near(report.integrated_current_c,expected[2]-previous[2],2e-13);
            near(report.charge_residual_c,0.0,2e-13);
            near(report.energy_residual_j,0.0,2e-13);
            near(prepared.phase_quadratures()[0],fs_math::det::cos(frequency*time),2e-12);
            near(prepared.phase_quadratures()[1],fs_math::det::sin(frequency*time),2e-12);
            assert_eq!(prepared.accepted_intervals(),interval);
            assert_eq!(prepared.elapsed_time_s(),time);
            assert_eq!(prepared.body().rod().material(),&source);
            previous = expected;
        }
    }
}

fn sources(tones: &[PiezoelectricTone], time: f64) -> [f64; 2] {
    let mut total = [0.0; 2];
    for tone in tones {
        let angle = tone.angular_frequency_rad_s*time;
        let (c,s) = (fs_math::det::cos(angle),fs_math::det::sin(angle));
        total[0] += tone.force_cosine_n*c+tone.force_sine_n*s;
        total[1] += tone.current_cosine_a*c+tone.current_sine_a*s;
    }
    total
}

// Independent x,v,Q equations plus the original physical integrands.
// No native matrix, Gramian, source carrier coordinate or storage callback.
fn derivative(state: [f64; 7], time: f64, d: f64, tones: &[PiezoelectricTone]) -> [f64; 7] {
    let theta = 2.5*d;
    let capacitance = 0.25*(2.0-10.0*d*d);
    let voltage = (state[2]-theta*state[0])/capacitance;
    let [force,current] = sources(tones,time);
    [state[1],force-2.5*state[0]+theta*voltage,current,
        force*state[1],current*voltage,current,voltage]
}
fn reference(state: &mut [f64; 7], start: f64, h: f64, d: f64, tones: &[PiezoelectricTone]) {
    state[3..].fill(0.0);
    let dt = h/512.0;
    for sub in 0..512 {
        let t = start+sub as f64*dt;
        let a = derivative(*state,t,d,tones);
        let b = derivative(core::array::from_fn(|i| state[i]+0.5*dt*a[i]),t+0.5*dt,d,tones);
        let c = derivative(core::array::from_fn(|i| state[i]+0.5*dt*b[i]),t+0.5*dt,d,tones);
        let e = derivative(core::array::from_fn(|i| state[i]+dt*c[i]),t+dt,d,tones);
        for i in 0..7 { state[i] += dt*(a[i]+2.0*b[i]+2.0*c[i]+e[i])/6.0; }
    }
}
fn spectrum() -> [PiezoelectricTone; 3] {
    [
        dc(0.001,0.0001),
        PiezoelectricTone {
            angular_frequency_rad_s: 1.25,
            force_cosine_n: 0.003, force_sine_n: -0.001,
            current_cosine_a: 0.0008, current_sine_a: 0.0003,
        },
        PiezoelectricTone {
            angular_frequency_rad_s: 3.1,
            force_cosine_n: -0.002, force_sine_n: 0.0007,
            current_cosine_a: -0.0004, current_sine_a: 0.0009,
        },
    ]
}

#[test]
fn g1_multitone_source_interference_matches_independent_motion_work_charge_and_voltage() {
    let tones = spectrum();
    for d in [-0.2,0.2] {
        let body = PiezoelectricBody::new(rod(d,3.0),0.002,0.001,0.0004).unwrap();
        let mut prepared = runtime(body,&tones,0.07);
        let mut expected = [0.002,0.001,0.0004,0.0,0.0,0.0,0.0];
        let mut differs_from_held_force = false;
        let mut differs_from_held_current = false;
        for interval in 0..24 {
            let start = interval as f64*0.07;
            let old_x = prepared.body().extension_m();
            reference(&mut expected,start,0.07,d,&tones);
            let report = prepared.step(controls(),|| false).unwrap();
            near(prepared.body().extension_m(),expected[0],2e-12);
            near(prepared.body().velocity_m_s(),expected[1],2e-12);
            near(prepared.body().charge_c(),expected[2],2e-12);
            near(report.mechanical_supply_j,expected[3],2e-13);
            near(report.electrical_supply_j,expected[4],2e-13);
            near(report.integrated_current_c,expected[5],2e-13);
            near(report.mean_voltage_v,expected[6]/0.07,2e-12);
            near(report.mean_velocity_m_s,(expected[0]-old_x)/0.07,2e-12);
            near(report.energy_residual_j,0.0,2e-13);
            near(report.charge_residual_c,0.0,2e-13);
            let [force,current] = sources(&tones,start);
            differs_from_held_force |= (report.mechanical_supply_j
                -force*(prepared.body().extension_m()-old_x)).abs() > 1e-10;
            differs_from_held_current |= (report.electrical_supply_j
                -current*expected[6]).abs() > 1e-10;
            for (tone,phase) in tones.iter().zip(prepared.phase_quadratures().chunks_exact(2)) {
                let angle = tone.angular_frequency_rad_s*prepared.elapsed_time_s();
                near(phase[0],fs_math::det::cos(angle),3e-12);
                near(phase[1],fs_math::det::sin(angle),3e-12);
            }
        }
        assert!(differs_from_held_force);
        assert!(differs_from_held_current);
        assert_eq!(prepared.tones(),tones.as_slice());
    }
}

#[test]
fn g1_dc_carrier_matches_the_original_held_current_owner_without_resetting_charge() {
    let body = PiezoelectricBody::new(rod(0.2,3.0),0.002,-0.001,0.0004).unwrap();
    let tones = [dc(0.004,0.0003)];
    let mut continuous = runtime(body.clone(),&tones,0.08);
    let mut held = super::prepared(body,0.08);
    for _ in 0..16 {
        let a = continuous.step(controls(),|| false).unwrap();
        let b = held.step(0.004,0.0003,controls(),|| false).unwrap();
        for (left,right) in continuous.body().state().into_iter().zip(held.body().state()) {
            near(left,right,2e-12);
        }
        near(a.mechanical_supply_j,b.mechanical_supply_j,2e-13);
        near(a.electrical_supply_j,b.electrical_supply_j,2e-13);
        near(a.integrated_current_c,0.0003*0.08,2e-14);
        near(a.mean_voltage_v,b.mean_voltage_v,2e-12);
        near(continuous.phase_quadratures()[0],1.0,2e-12);
        near(continuous.phase_quadratures()[1],0.0,2e-12);
    }
}

#[test]
fn g4_late_cancellation_preserves_body_phase_clock_and_report_then_retries_identically() {
    let body = PiezoelectricBody::new(rod(0.2,3.0),0.002,0.001,0.0004).unwrap();
    let tones = spectrum();
    let mut probe = runtime(body.clone(),&tones,0.07);
    let mut polls = 0;
    let expected = probe.step(controls(),|| { polls += 1; false }).unwrap();
    assert!(polls > tones.len());
    let mut prepared = runtime(body.clone(),&tones,0.07);
    let phase = prepared.phase_quadratures().to_vec();
    let mut seen = 0;
    assert!(matches!(prepared.step(controls(),|| {
        seen += 1; seen == polls
    }),Err(PiezoelectricCouplingError::Cancelled)));
    assert_eq!(prepared.body(),&body);
    assert_eq!(prepared.phase_quadratures(),phase.as_slice());
    assert_eq!(prepared.accepted_intervals(),0);
    assert_eq!(prepared.elapsed_time_s(),0.0);
    assert!(prepared.last_step().is_none());
    let retry = prepared.step(controls(),|| false).unwrap();
    assert_eq!(retry,expected);
    assert_eq!(prepared.body(),probe.body());
    assert_eq!(prepared.phase_quadratures(),probe.phase_quadratures());
    assert_eq!(prepared.accepted_intervals(),1);
    assert_eq!(prepared.elapsed_time_s(),0.07);
    let accepted_body = prepared.body().clone();
    let accepted_phase = prepared.phase_quadratures().to_vec();
    let mut invalid = controls();
    invalid.charge_tolerance_c = f64::NAN;
    assert!(prepared.step(invalid,|| false).is_err());
    assert_eq!(prepared.body(),&accepted_body);
    assert_eq!(prepared.phase_quadratures(),accepted_phase.as_slice());
    assert_eq!(prepared.accepted_intervals(),1);
    assert_eq!(prepared.last_step(),Some(&retry));
}

#[test]
fn g4_complete_tone_window_refuses_hidden_strain_excursion_before_phase_advances() {
    let frequency = fs_math::det::sqrt(2.5);
    let duration = core::f64::consts::TAU/frequency;
    let tones = [PiezoelectricTone {
        angular_frequency_rad_s: frequency,
        force_cosine_n: 0.12, force_sine_n: 0.0,
        current_cosine_a: 0.0, current_sine_a: 0.0,
    }];
    let body = PiezoelectricBody::new(rod(0.0,3.0),0.0,0.0,0.0).unwrap();
    let mut prepared = runtime(body.clone(),&tones,duration);
    // Exact resonant x is zero again at the endpoint, but at 3/4 period
    // |x|=3*pi*F/(4*k)=0.113... m, beyond the original 0.1 m strain limit.
    assert!(3.0*core::f64::consts::PI*0.12/(4.0*2.5) > 0.1);
    assert!(prepared.step(controls(),|| false).is_err());
    assert_eq!(prepared.body(),&body);
    assert_eq!(prepared.phase_quadratures(),&[1.0,0.0]);
    assert_eq!(prepared.accepted_intervals(),0);
    assert!(prepared.last_step().is_none());
}

#[test]
fn g4_tone_preparation_admits_full_source_extent_and_refuses_invalid_programs() {
    let body = PiezoelectricBody::new(rod(0.2,3.0),0.0,0.0,0.0).unwrap();
    let (split,dense) = budgets();
    let tone = spectrum()[1];
    assert!(PreparedPiezoelectricTones::new(body.clone(),&[tone],0.1,
        LinearSplitBudget { max_modes: 4, ..split },dense,|| false).is_err());
    assert!(PreparedPiezoelectricTones::new(body.clone(),&[tone],0.1,split,
        DenseExponentialBudget { max_dimension: 7, ..dense },|| false).is_err());
    assert!(PreparedPiezoelectricTones::new(body.clone(),&[tone,tone],0.1,split,dense,|| false).is_err());
    for invalid in [
        PiezoelectricTone { angular_frequency_rad_s: f64::NAN, ..tone },
        PiezoelectricTone { current_cosine_a: f64::INFINITY, ..tone },
        PiezoelectricTone { angular_frequency_rad_s: -1.0, ..tone },
        PiezoelectricTone { angular_frequency_rad_s: 0.0, ..tone },
    ] {
        assert!(PreparedPiezoelectricTones::new(body.clone(),&[invalid],0.1,
            split,dense,|| false).is_err());
    }
    assert!(matches!(PreparedPiezoelectricTones::new(body,&[tone],0.1,
        split,dense,|| true),Err(PiezoelectricCouplingError::Cancelled)));
}
