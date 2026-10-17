use super::*;
use super::super::{geometry::Course, linear::BoardMode, steinway_scale};

fn medium() -> Medium { Medium { density: 1.2, sound_speed: 343.0 } }

#[test]
fn actual_piano_target_sweep_preserves_cancelling_power_and_half_spaces() {
    let cells: Vec<_> = [-1.0, 0.0, 1.0].into_iter().zip([1.0, -2.0, 1.0])
        .map(|(x, shape)| AreaSample { position_m: [0.02*x, 0.0], area_m2: 1e-4,
            mode_shape: vec![shape, 1.0] }).collect();
    let source = PreparedSurface::new(&cells).unwrap();
    let omega: Vec<_> = [100.0, 400.0, 1000.0].into_iter().map(|hz|TAU*hz).collect();
    // Exercise the production target path BEFORE either full or projected fit.
    let one = sample_impedances(&source, &omega, medium(), 1).unwrap();
    let two = sample_impedances(&source, &omega, medium(), 2).unwrap();
    for (i, &w) in omega.iter().enumerate() {
        let legacy = source.impedance(w,medium()).unwrap();
        assert!(legacy[0].re < 0.0 && one[i][0].re > 0.0);
        assert_eq!(one[i][1],one[i][2]);
        assert_eq!(one[i],source.power_consistent_impedance(w,medium()).unwrap());
        for j in 0..4 {
            assert_eq!(two[i][j],one[i][j].scale(2.0));
            assert_eq!(one[i][j].im,legacy[j].im);
        }
    }
    // A bad late frequency cannot return a partially prepared target bank.
    assert!(sample_impedances(&source,&[omega[0],f64::NAN],medium(),1).is_err());
}

#[test]
fn loaded_physical_shapes_reach_the_unchanged_passive_fitter_and_projection() {
    let course = steinway_scale::courses().unwrap().into_iter().find(|c|c.midi==69).unwrap();
    let course = Course { unison:1, duplex_length_m:0.0, ..course };
    let mut bridge = [0.0;88]; bridge[69-21]=2.0;
    let board = [BoardMode { frequency_hz:250.0, damping_ratio:0.0, bridge, volume:0.001 }];
    let bank = Bank::new(&[course],&board,192_000,8_000.0,20,false).unwrap();
    let surface = [SurfaceSample { position_m:[0.0,0.0,0.0],area_m2:0.001,mode_shape:vec![1.0] }];
    let options = Options::parse("1,40,500,65").unwrap();
    let loaded_shape = bank.project_board_shape(&[1.0]).unwrap();
    assert!((loaded_shape[0]-1.0).abs()>1e-5,"must use the actual string-loaded basis");
    let source = PreparedSurface::new(&[AreaSample { position_m:[0.0,0.0],area_m2:0.001,
        mode_shape:loaded_shape }]).unwrap();
    // A complete load is sampled linearly and fitted with the travel-time dictionary;
    // one cell has no travel time, so one spacing spans the band.
    let last=(options.samples-1) as f64;
    let mut omega:Vec<_>=(0..options.samples).map(|i|
        TAU*(options.min_hz+(options.max_hz-options.min_hz)*i as f64/last)).collect();
    omega[0]=TAU*options.min_hz;omega[options.samples-1]=TAU*options.max_hz;
    // Independent assembly of the precise fitting inputs, not the adapter's
    // sampling helper. The full adapter must produce the identical realization.
    let data:Vec<_>=omega.iter().map(|&w|source.power_consistent_impedance(w,medium()).unwrap()).collect();
    let expected=fit::fit_resolved(&omega,&data,1,omega[options.samples-1]-omega[0]).unwrap();
    let full=prepare(&surface,&bank,&options,medium(),48_000,None).unwrap();
    assert!(full.report.to_string().contains("Power-consistent area-point resistance"));
    for hz in [40.0,193.0,500.0] {
        assert_eq!(full.model.impedance(TAU*hz).unwrap(),expected.model.impedance(TAU*hz).unwrap());
    }
    // The explicit projection keeps the logarithmic grid and the shared fixed
    // dictionary; a lossless one-port coordinate map leaves both fits on the load.
    let projected=prepare(&surface,&bank,&options,medium(),48_000,Some(1)).unwrap();
    assert!(projected.basis.is_some());
    assert_eq!(projected.report.board_modes,1);
    assert_eq!(projected.report.acoustic_ports,1);
    assert!(projected.report.surface_rms_error.unwrap()<1e-12);
    for hz in [40.0,193.0,500.0] {
        let load=source.power_consistent_impedance(TAU*hz,medium()).unwrap()[0];
        for model in [&full.model,&projected.model] {
            let z=model.impedance(TAU*hz).unwrap()[0];
            assert!((z-load).abs()<0.05*load.abs(),"fit left the sampled load at {hz} Hz");
        }
    }
}
