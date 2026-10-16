//! Keep the real board/lid geometry and the original radiation-power gates.
//! The regression crosses the former scene-radius switch without changing
//! the source geometry, material cards, receiver locations or frequency grid.
use super::*;
use super::super::{exterior_geometry,playback};

#[test]
fn native_skin_and_lid_use_the_same_source_equation_for_pressure_and_reaction() {
    let (board,courses,_,_)=super::super::tests::small_source_inputs();
    let text=exterior_geometry::tests::specification()
        .replace("moving,skin","moving,soundboard_skin")
        .replace("band-hz,40,400,17","band-hz,40,300,41")
        .replace("board-band-hz,400","board-band-hz,300");
    let spec=Specification::read(&format!("{text}receiver-m,0.05,0.05,1\n")).unwrap();
    let controls=playback::Controls::from_texts(&courses,None,None,None).unwrap();
    let scene=super::super::prepare_controlled_body(&board,courses,None,spec,
        &playback::Options::default(),controls,true).unwrap();
    let skin_count=scene.boundary.surface.areas().len();
    let skin_weights=scene.boundary.weights.clone();
    let asset=exterior_geometry::tests::box_obj("Lid",[20.,20.,100.],[50.,50.,4.]);
    let manifest="frankensim-piano-rigid-assembly-v1\nsource,estimated,existing rigid-scene regression\npart,lid,lid.obj,Lid,0.001,0,0,0\npose,lid,0.02,0.02,0.1,1,0,0,30,0,0,0\n";
    let assembly=exterior_geometry::rigid::Assembly::from_text(manifest,|_|Ok(asset.clone())).unwrap();
    let boundary=assembly.attach(scene.boundary).unwrap();
    let old_switch=0.5*scene.spec.medium.sound_speed/(TAU*boundary.radius);
    assert!(old_switch>280. && old_switch<287.);
    assert_eq!(boundary.surface.areas().len(),28);
    assert_eq!(boundary.components,2);
    for (row,original) in boundary.weights.iter().zip(&skin_weights) {
        assert_eq!(&row[..skin_count],original);
        assert!(row[skin_count..].iter().all(|x|*x==0.));
    }
    let policy=GeometryPolicy::new(&boundary.surface).unwrap();
    let pressure=boundary.sample(&scene.spec).unwrap();
    // Use the entire old failing grid; the two existing full playback/CLI
    // regressions remain unchanged and cover the downstream fitted paths.
    for (f,&w) in pressure.omega.iter().enumerate() {
        let k=w/scene.spec.medium.sound_speed;
        assert_eq!(policy.formulation(k).unwrap(),Formulation::PlainCbie);
        let loading=sample(&boundary,&scene.spec,w).unwrap();
        for (channel,row) in loading.receiver_transfer.iter().enumerate() {
            for (input,&h) in row.iter().enumerate() {
                let acceleration_to_velocity=pressure.values[channel][input][f]*C64::new(0.,-w);
                assert!((acceleration_to_velocity-h).abs()<1e-9*(1.+h.abs()));
            }
        }
        if [0,20,37,38,40].contains(&f) {
            let fields:Vec<Vec<C64>>=boundary.weights.iter().map(|row|
                row.iter().map(|v|C64::new(*v,0.)).collect()).collect();
            let refs:Vec<_>=fields.iter().map(Vec::as_slice).collect();
            let direct=helmholtz::solve_radiation_batch(&boundary.surface,k,scene.spec.medium,
                &refs,Formulation::PlainCbie).unwrap();
            let r=boundary.weights.len();
            for (j,solution) in direct.iter().enumerate() {
                assert!(solution.radiated_power_roundoff_interval.0>0.);
                for (i,shape) in boundary.weights.iter().enumerate() {
                    let impedance=shape.iter().zip(boundary.surface.areas()).zip(&solution.pressure)
                        .fold(C64::ZERO,|z,((g,a),p)|z+p.scale(g*a));
                    assert_eq!(loading.impedance[i*r+j],impedance);
                }
            }
        }
    }
    assert!(scene.piano.bank.q.iter().chain(&scene.piano.bank.v).all(|x|*x==0.));
    assert_eq!(scene.piano.accounting.input_work_j,0.);
}
