//! Geometry import -> source strings/felt/shanks -> the existing physical audio.
use super::{import, crowned_board, tests::{fixture, SPEC, PAIRS}};
use super::super::super::super::{audio::AudioStream, engine::{Instrument, ShankGeometry},
    performance::Performance, steinway_scale, steinway_d, mesh_import};
use std::fmt::Write;

#[test]
fn source_piano_pressure_changes_with_imported_thickness_and_preserves_block_determinism() {
    let run = |factor, split| {
        let imported = import(&fixture(0.015,factor),SPEC,PAIRS).unwrap();
        let board = crowned_board::CrownedBoard::read(&imported.fsb).unwrap().prepare(&[69],400.).unwrap();
        let mut courses = steinway_scale::courses().unwrap(); courses.retain(|c|c.midi==69);
        let materials = courses.iter().map(steinway_scale::hammer_material).collect::<Result<Vec<_>,_>>().unwrap();
        let piano = Instrument::new_with_course_shanks(courses,&board.modes,48_000,4,24,true,
            materials,ShankGeometry::published()).unwrap();
        let score = Performance::read("sample,event,key,value\n0,note_on,69,2\n1024,note_off,69,0\n",&[69],2400).unwrap();
        let mut stream = AudioStream::new(piano,score,Some(&board.surface),[0.675,1.,1.],
            fs_bem::helmholtz::Medium::air(),1.).unwrap();
        let mut pressure = vec![0.;2400];
        if split { for block in pressure.chunks_mut(127) { stream.render_block(block).unwrap(); } }
        else { stream.render_block(&mut pressure).unwrap(); }
        assert!(pressure.iter().all(|p|p.is_finite()));
        assert!(pressure.iter().any(|p|p.abs()>1e-10));
        let piano = stream.instrument();
        assert!(piano.accounting.felt_loss_j>0.);
        assert!(piano.accounting.shank_loss_j>0.);
        assert!((piano.accounting.input_work_j-piano.energy_j()-piano.accounting.dissipated_j()).abs()<1e-7);
        let (wav,_) = fs_couple::pcm_wav::encode_pcm16_wav(&pressure,48_000,2.).unwrap();
        assert_eq!(&wav[..4],b"RIFF");
        (pressure,wav)
    };
    let a = run(1.0,false); let repeat = run(1.0,true); let b = run(1.4,false);
    assert_eq!(a,repeat);
    assert_ne!(a.0,b.0); assert_ne!(a.1,b.1);
}

#[test]
fn the_source_model_d_retains_all_bridges_supports_and_bonded_rib_paths() {
    let native = steinway_d::build(4).unwrap().geometry;
    let (source,spec) = mesh_import::export(&native).unwrap();
    let doc = fs_io::obj::read_obj_document(&source).unwrap();
    let n = doc.soup.positions.len();
    let mut obj = String::new();
    // Constant-thickness test skins around the source planform, NOT measured
    // Model D surfaces or a replacement for the source tapered section cards.
    for z in [0.004,-0.004] { for p in &doc.soup.positions {
        writeln!(obj,"v {:.17e} {:.17e} {z}",p.x,p.y).unwrap();
    }}
    obj.push_str("o soundboard\n");
    for r in &doc.regions {
        writeln!(obj,"usemtl {}",r.material.as_ref().unwrap()).unwrap();
        for face in r.triangles.clone() {
            let t = doc.soup.triangles[face];
            writeln!(obj,"f {} {} {}",t[0]+1,t[1]+1,t[2]+1).unwrap();
        }
    }
    obj.push_str("o underside\n");
    for t in &doc.soup.triangles {
        writeln!(obj,"f {} {} {}",t[0] as usize+n+1,t[2] as usize+n+1,t[1] as usize+n+1).unwrap();
    }
    let skins = "frankensim-board-skins-v1\nlower,underside\nthickness,geometry\nthickness-range,0.007,0.009\npairing,projected,1e-8\n";
    let imported = import(&obj,&spec,skins).unwrap();
    assert_eq!(imported.source_vertices.len(),n);
    for prefix in ["fixed,","stiffener,","bridge,"] {
        assert_eq!(native.lines().filter(|l|l.starts_with(prefix)).count(),
            imported.fsb.lines().filter(|l|l.starts_with(prefix)).count());
    }
    assert_eq!(imported.fsb.lines().filter(|l|l.starts_with("bridge,")).count(),88);
    let beams = |s: &str| s.lines().filter(|l|l.starts_with("stiffener,"))
        .map(|l|l.split(',').skip(1).map(|x|x.parse::<f64>().unwrap()).collect::<Vec<_>>()).collect::<Vec<_>>();
    assert_eq!(beams(&native),beams(&imported.fsb));
    let board = crowned_board::CrownedBoard::read(&imported.fsb).unwrap();
    assert_eq!(board.max_height_m,0.);
    assert!(board.mass_kg.is_finite() && board.mass_kg>0.);
}
