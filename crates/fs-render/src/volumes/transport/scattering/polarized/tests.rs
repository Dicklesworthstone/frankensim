use super::*;
use crate::dielectric::{BeerLambertAbsorption, CauchyIor, DielectricGlass, GlassProvenance};
use crate::volumes::{VolumeGrid, transport::{AbsorbingGrid, AbsorbingGridBudget, MassAbsorptionSpectrum}};
use crate::volumes::transport::scattering::MassScatteringSpectrum;
use fs_geom::Point3;

fn beam(power: Stokes) -> PolarizedBeam {
    PolarizedBeam::try_new(
        Point3::new(0.5,0.5,-0.5),
        PolarizationFrame::try_new(Vec3::new(0.0,0.0,1.0),Vec3::new(1.0,0.0,0.0)).unwrap(),
        550.0,
        DielectricGlass::new(CauchyIor::try_constant(1.0).unwrap(),BeerLambertAbsorption::CLEAR,GlassProvenance::Custom),
        power,
    ).unwrap()
}
fn near(a:f64,b:f64) { assert!((a-b).abs() <= 3e-12*a.abs().max(b.abs()).max(1.0),"{a:.17e} != {b:.17e}"); }
fn same(a:Stokes,b:Stokes) {
    for (a,b) in [a.i,a.q,a.u,a.v].into_iter().zip([b.i,b.q,b.u,b.v]) { near(a,b); }
}

#[test]
fn rayleigh_matrix_matches_scalar_and_retains_physical_states() {
    for delta in [0.0,0.03,0.5,0.85] {
        let phase=PhaseFunction::rayleigh_depolarized(delta).unwrap();
        for j in 0..=64 {
            let c=-1.0+2.0*j as f64/64.0;
            let matrix=phase.rayleigh_mueller(c).unwrap();
            assert_eq!(matrix.0[0][0],phase.eval_cos(c).unwrap());
            for (q,u,v) in [(0.0,0.0,0.0),(1.0,0.0,0.0),(-1.0,0.0,0.0),
                (0.0,1.0,0.0),(0.0,0.0,1.0),(0.3,-0.4,0.5)] {
                let out=matrix.apply(Stokes::try_new(1.0,q,u,v).unwrap());
                assert!(out.i>=0.0);
                assert!(out.q.hypot(out.u).hypot(out.v)<=out.i+3e-15);
            }
        }
        let out=phase.rayleigh_mueller(0.0).unwrap().apply(Stokes::unpolarized(1.0));
        near((out.i-out.q)/(out.i+out.q),delta);
        let forward=phase.rayleigh_mueller(1.0).unwrap();
        near(forward.0[2][2]/forward.0[0][0],1.0-delta);
        near(forward.0[3][3]/forward.0[0][0],1.0-2.0*delta);
    }
}

#[test]
fn right_angle_dipole_extinction_and_polarized_sample_weights() {
    let phase=PhaseFunction::Rayleigh;
    let unpolarized=phase.sample_polarized(beam(Stokes::unpolarized(1.0)),[0.5,0.0]).unwrap();
    same(unpolarized.beam.power(),Stokes { i:1.0,q:1.0,u:0.0,v:0.0 });
    near(unpolarized.pdf,3.0/(16.0*core::f64::consts::PI));
    // The native azimuth-zero sample goes along +x: an x dipole cannot emit
    // along itself, while a y dipole has twice the unpolarized intensity.
    let dark=phase.sample_polarized(beam(Stokes::try_new(1.0,1.0,0.0,0.0).unwrap()),[0.5,0.0]).unwrap();
    near(dark.beam.power().i,0.0);
    let bright=phase.sample_polarized(beam(Stokes::try_new(1.0,-1.0,0.0,0.0).unwrap()),[0.5,0.0]).unwrap();
    near(bright.beam.power().i,2.0);
    assert!(bright.beam.power().i>1.0);
}

#[test]
fn phase_frames_are_covariant_and_backscatter_reverses_handedness() {
    let initial=beam(Stokes::try_new(1.0,0.2,-0.3,0.4).unwrap());
    let changed=initial.with_reference_axis(Vec3::new(1.0,2.0,3.0)).unwrap();
    for delta in [0.0,0.03,0.5] {
        let phase=PhaseFunction::rayleigh_depolarized(delta).unwrap();
        for sample in [[0.0,0.0],[0.5,0.3],[0.93,0.7]] {
            let a=phase.sample_polarized(initial,sample).unwrap();
            let b=phase.sample_polarized(changed,sample).unwrap();
            let b=b.beam.with_reference_axis(a.beam.frame().horizontal()).unwrap();
            same(a.beam.power(),b.power());
            assert_eq!(a.beam.medium(),initial.medium());
            assert_eq!(a.beam.ray().origin,initial.ray().origin);
        }
    }
    let circular=beam(Stokes::try_new(1.0,0.0,0.0,1.0).unwrap());
    same(PhaseFunction::Rayleigh.sample_polarized(circular,[0.0,0.0]).unwrap().beam.power(),
        Stokes { i:1.0,q:0.0,u:0.0,v:-1.0 });
}

#[test]
fn angular_estimator_conserves_arbitrary_incident_polarization() {
    for delta in [0.0,0.03,0.85] {
        let phase=PhaseFunction::rayleigh_depolarized(delta).unwrap();
        for state in [Stokes::unpolarized(1.0),Stokes::try_new(1.0,1.0,0.0,0.0).unwrap(),
            Stokes::try_new(1.0,0.0,1.0,0.0).unwrap(),Stokes::try_new(1.0,0.0,0.0,1.0).unwrap()] {
            let mut sum=0.0;
            for j in 0..32 { for k in 0..64 {
                sum+=phase.sample_polarized(beam(state),[(j as f64+0.5)/32.0,(k as f64+0.5)/64.0])
                    .unwrap().beam.power().i/2048.0;
            } }
            near(sum,1.0);
        }
    }
}

#[test]
fn original_free_flight_owns_cell_positions_albedo_and_escape_probability() {
    super::super::tests::with_cx(|cx,_| {
        let density=[1.0,0.0,3.0];
        let absorption=AbsorbingGrid::try_new(cx,
            VolumeGrid::new([1,1,3],&density,[0.0;3],[1.0;3]),
            MassAbsorptionSpectrum::constant(0.25).unwrap(),None,
            AbsorbingGridBudget { max_grid_cells:3,max_steps_per_segment:16 }).unwrap();
        let grid=ScatteringGrid::try_new(cx,absorption,MassScatteringSpectrum::constant(0.75).unwrap(),PhaseFunction::Rayleigh).unwrap();
        let polarized=RayleighGrid::try_new(cx,&grid).unwrap();
        let input=beam(Stokes::unpolarized(1.0));
        for u in [(-0.5_f64).exp(),(-2.0_f64).exp(),(-10.0_f64).exp()] {
            let native=grid.sample_segment(cx,input.ray().origin,input.ray().dir,4.0,&[550.0],0.5,u).unwrap();
            let result=polarized.sample_beam_segment(cx,input,4.0,[u,0.5,0.0]).unwrap();
            assert_eq!(result.event,native.event);
            assert_eq!(result.escape_probability,native.escape_probability);
            near(result.beam.power().i,native.weight[0]);
            if let Some(event)=result.event {
                assert_ne!(event.cell_index,1);
                assert_eq!(result.beam.ray().origin,event.position);
                near(result.beam.power().i,0.75);
                near(result.beam.power().q,0.75);
                assert!(result.phase_pdf.is_some());
            } else {
                assert_eq!(result.beam.ray().origin,input.ray().at(4.0));
                same(result.beam.power(),input.power());
                near(result.escape_probability,(-4.0_f64).exp());
                assert_eq!(result.phase_pdf,None);
            }
        }
    });
}

#[test]
fn absorption_only_limit_and_unsupported_models_refuse_without_double_loss() {
    super::super::tests::with_cx(|cx,_| {
        let density=[1.0];
        let make=|coefficient| ScatteringGrid::try_new(cx,AbsorbingGrid::try_new(cx,
            VolumeGrid::new([1,1,1],&density,[0.0;3],[1.0;3]),
            MassAbsorptionSpectrum::constant(0.4).unwrap(),None,
            AbsorbingGridBudget { max_grid_cells:1,max_steps_per_segment:8 }).unwrap(),
            MassScatteringSpectrum::constant(coefficient).unwrap(),PhaseFunction::Isotropic).unwrap();
        let grid=make(0.0);
        let polarized=RayleighGrid::try_new(cx,&grid).unwrap();
        let initial=beam(Stokes::try_new(1.0,0.0,0.6,0.8).unwrap());
        let out=polarized.sample_beam_segment(cx,initial,2.0,[0.5,0.5,0.5]).unwrap();
        same(out.beam.power(),scale_power(initial.power(),(-0.4_f64).exp()).unwrap());
        assert_eq!(out.event,None);
        assert!(matches!(RayleighGrid::try_new(cx,&make(1.0)),Err(PolarizedScatteringError::UnsupportedPhase)));
        for (index,absorption) in [(1.33,0.0),(1.0,0.1)] {
            let medium=DielectricGlass::new(CauchyIor::try_constant(index).unwrap(),
                BeerLambertAbsorption::try_constant(absorption).unwrap(),GlassProvenance::Custom);
            let invalid=PolarizedBeam::try_new(initial.ray().origin,initial.frame(),550.0,medium,initial.power()).unwrap();
            assert_eq!(polarized.sample_beam_segment(cx,invalid,2.0,[0.5;3]),Err(PolarizedScatteringError::CarrierMedium));
        }
        assert!(polarized.sample_beam_segment(cx,initial,f64::INFINITY,[0.5;3]).is_err());
        assert!(polarized.sample_beam_segment(cx,initial,2.0,[1.0,0.5,0.5]).is_err());
    });
}
