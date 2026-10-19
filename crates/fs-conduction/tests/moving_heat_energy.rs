//! Analytic schedule-energy integrals and their authored-control derivatives.
mod support;
use fs_conduction::bc::history::march::moving::{MovingHeatBudget,MovingHeatSource};
use fs_conduction::{COMPONENT_POWER_DIMS,ConductionError};
use fs_qty::{Dims,QtyAny};
use fs_scenario::signal::{Interp,TimeSignal};
use support::{with_cx,with_cancelled_cx};

fn source(power: TimeSignal) -> MovingHeatSource {
    MovingHeatSource::new(std::array::from_fn(|_| TimeSignal::Constant(
        QtyAny::new(0.,Dims([1,0,0,0,0,0])))),[1.;3],power,0.,MovingHeatBudget {
        max_vertices:8,max_elements:6,min_active_vertices:1,power_relative_tolerance:1e-10,
    }).unwrap()
}
fn table(interp: Interp) -> MovingHeatSource {
    source(TimeSignal::Table {times:vec![1.,2.,4.],values:vec![2.,4.,8.],dims:COMPONENT_POWER_DIMS,interp})
}

#[test]
fn g1_short_pulses_linear_segments_and_clamped_tails_have_analytic_energy() {
    with_cx(|cx| {
        let pulse=source(TimeSignal::Table {times:vec![0.,1.,1.125,2.],values:vec![0.,120.,0.,0.],
            dims:COMPONENT_POWER_DIMS,interp:Interp::Hold});
        let e=pulse.integrate_absorbed_power(cx,0.,2.,4).unwrap();
        assert_eq!(e.energy_j,15.);
        assert_eq!(e.power_value_weights_s,vec![1.,0.125,0.875,0.]);
        let ramp=source(TimeSignal::Ramp {t_start:2.,t_end:5.,
            from:QtyAny::new(2.,COMPONENT_POWER_DIMS),to:QtyAny::new(8.,COMPONENT_POWER_DIMS)});
        assert_eq!(ramp.integrate_absorbed_power(cx,0.,7.,2).unwrap().energy_j,35.);
        assert_eq!(table(Interp::Linear).integrate_absorbed_power(cx,0.,5.,3).unwrap().power_value_weights_s,vec![1.5,1.5,2.]);
        assert_eq!(table(Interp::Linear).integrate_absorbed_power(cx,0.,5.,3).unwrap().energy_j,25.);
        assert_eq!(table(Interp::Hold).integrate_absorbed_power(cx,0.,5.,3).unwrap().energy_j,20.);
        let one=source(TimeSignal::Table {times:vec![1.],values:vec![3.],dims:COMPONENT_POWER_DIMS,interp:Interp::Hold});
        assert_eq!(one.integrate_absorbed_power(cx,-1.,2.,1).unwrap().energy_j,9.);
    });
}

#[test]
fn g3_integral_derivatives_match_original_schedule_perturbations_and_partitioning() {
    with_cx(|cx| {
        for interp in [Interp::Hold,Interp::Linear] {
            let source=table(interp);
            let all=source.integrate_absorbed_power(cx,0.25,5.,3).unwrap();
            let first=source.integrate_absorbed_power(cx,0.25,1.75,3).unwrap();
            let last=source.integrate_absorbed_power(cx,1.75,5.,3).unwrap();
            assert!((all.energy_j-first.energy_j-last.energy_j).abs()<1e-13);
            assert!((all.power_value_weights_s.iter().sum::<f64>()-4.75).abs()<1e-14);
            for i in 0..3 {
                let(mut plus,mut minus)=(vec![2.,4.,8.],vec![2.,4.,8.]);
                plus[i]+=1e-4;minus[i]-=1e-4;
                let p=source.with_absorbed_power_values(&plus,0.,3).unwrap().integrate_absorbed_power(cx,0.25,5.,3).unwrap();
                let m=source.with_absorbed_power_values(&minus,0.,3).unwrap().integrate_absorbed_power(cx,0.25,5.,3).unwrap();
                assert!(((p.energy_j-m.energy_j)/2e-4-all.power_value_weights_s[i]).abs()<1e-9);
            }
        }
    });
}

#[test]
fn g0_integral_refusals_and_zero_duration_never_change_the_source() {
    let original=table(Interp::Hold);let saved=original.clone();
    with_cancelled_cx(|cx| assert!(matches!(original.integrate_absorbed_power(cx,0.,5.,3),Err(ConductionError::Cancelled {..}))));
    with_cx(|cx| {
        for (a,b) in [(2.,1.),(f64::NAN,1.),(0.,f64::INFINITY),(-f64::MAX,f64::MAX)] {
            assert!(original.integrate_absorbed_power(cx,a,b,3).is_err());
        }
        assert!(original.integrate_absorbed_power(cx,0.,5.,2).is_err());
        let zero=original.integrate_absorbed_power(cx,2.,2.,3).unwrap();
        assert_eq!(zero.energy_j,0.);assert_eq!(zero.power_value_weights_s,vec![0.;3]);
        let huge=source(TimeSignal::Constant(QtyAny::new(f64::MAX,COMPONENT_POWER_DIMS)));
        assert!(huge.integrate_absorbed_power(cx,0.,2.,1).is_err());
        assert_eq!(original,saved);
    });
}
