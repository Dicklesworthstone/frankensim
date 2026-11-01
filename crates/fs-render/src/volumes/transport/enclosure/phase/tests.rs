use super::*;
use crate::volumes::transport::{OpticalMixtureBudget};
use crate::volumes::transport::scattering::{MassScatteringSpectrum, ScatteringMixtureComponent, SpectralPhaseMixture};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use core::f64::consts::PI;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>, &CancelGate)->R)->R {
    let gate=CancelGate::new_clock_free();
    let pool=fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena|f(&Cx::new(&gate,arena,StreamKey{seed:173,kernel_id:0x5350_454e,tile:0,iteration:0},
        Budget::INFINITE,ExecMode::Deterministic),&gate))
}
fn contains(x:Interval,y:f64) {
    assert!(x.lo()<=y && y<=x.hi(),"{x:?} does not enclose {y}");
}
#[test]
fn spectral_interpolation_and_gray_limits_have_independent_phase_formulas() {
    with_cx(|cx,_| {
        let phase=PhaseFunction::spectral_table(&[(380.0,PhaseFunction::Isotropic),(780.0,PhaseFunction::Rayleigh)]).unwrap();
        for lambda in [380.0,400.0,580.0,700.0,780.0] {
            for mu in [-1.0,-0.3,0.0,0.8,1.0] {
                let t=(lambda-380.0)/400.0;
                let truth=((1.0-t)+t*0.75*(1.0+mu*mu))/(4.0*PI);
                let bound=phase.enclosure_at(cx,lambda,Interval::point(mu)).unwrap();
                contains(bound,truth);
                assert!(bound.width()<1e-13);
            }
        }
        let a=PhaseFunction::Isotropic.enclosure_cosine(cx,Interval::new(-1.0,1.0)).unwrap();
        assert_eq!(phase.enclosure_at(cx,380.0,Interval::new(-1.0,1.0)).unwrap(),a);
        assert!(matches!(phase.enclosure_cosine(cx,Interval::point(0.0)),Err(VolumeTransportError::Wavelength)));
        assert!(matches!(phase.enclosure_at(cx,379.0,Interval::point(0.0)),Err(VolumeTransportError::Wavelength)));
    });
}
#[test]
fn cross_section_weighted_mixtures_are_not_equal_weight_phase_averages() {
    with_cx(|cx,_| {
        let components=[
            ScatteringMixtureComponent{mass_weight:0.25,scattering:MassScatteringSpectrum::constant(0.5).unwrap(),phase:PhaseFunction::Isotropic},
            ScatteringMixtureComponent{mass_weight:0.75,scattering:MassScatteringSpectrum::rayleigh(400.0,2.0).unwrap(),phase:PhaseFunction::Rayleigh},
        ];
        let phase=PhaseFunction::SpectralMixture(SpectralPhaseMixture::try_new(cx,&components,OpticalMixtureBudget::default()).unwrap());
        for lambda in [400.0,550.0,700.0] {
            for mu in [-1.0,0.0,0.5,1.0] {
                let r=400.0/lambda;let a=0.125;let b=1.5*r*r*r*r;
                let truth=(a+b*0.75*(1.0+mu*mu))/(a+b)/(4.0*PI);
                let bound=phase.enclosure_at(cx,lambda,Interval::point(mu)).unwrap();
                contains(bound,truth);assert!(bound.width()<1e-12);
            }
        }
        let zeros=[ScatteringMixtureComponent{mass_weight:0.0,..components[0].clone()}];
        let empty=PhaseFunction::SpectralMixture(SpectralPhaseMixture::try_new(cx,&zeros,OpticalMixtureBudget::default()).unwrap());
        assert_eq!(empty.enclosure_at(cx,500.0,Interval::point(1.0)).unwrap(),
            PhaseFunction::Isotropic.enclosure_cosine(cx,Interval::point(1.0)).unwrap());
    });
}
#[test]
fn angular_interior_peaks_and_invalid_laws_are_not_hidden() {
    with_cx(|cx,gate| {
        let peak=PhaseFunction::tabulated(&[(-1.0,0.0),(-0.01,0.0),(0.0,1.0),(0.01,0.0),(1.0,0.0)]).unwrap();
        let bound=peak.enclosure_at(cx,500.0,Interval::new(-0.2,0.2)).unwrap();
        assert_eq!(bound.lo(),0.0);contains(bound,peak.eval_cos(0.0).unwrap());
        for g in [1.0,f64::NAN,f64::INFINITY] {
            assert!(matches!(PhaseFunction::HenyeyGreenstein{g}.enclosure_at(cx,500.0,Interval::point(0.0)),Err(VolumeTransportError::Phase)));
        }
        assert!(PhaseFunction::Isotropic.enclosure_cosine(cx,Interval::new(-1.01,0.0)).is_err());
        gate.request();
        assert!(matches!(peak.enclosure_at(cx,500.0,Interval::point(0.0)),Err(VolumeTransportError::Cancelled)));
    });
}
