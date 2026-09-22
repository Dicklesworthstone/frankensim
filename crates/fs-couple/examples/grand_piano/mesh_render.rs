//! Imported structural geometry with the existing source-derived Model D
//! strings, per-key wool/Prony cards, shanks, bridge coupling and air observer.
//! This is an offline composition, not another piano solver or PCM encoder.
use super::{audio::AudioStream, board_geometry::{BoardGeometry, PreparedBoard},
    engine::{Instrument, ShankGeometry}, geometry::Course, performance::{self, Performance},
    steinway_scale};

pub const RATE: u32 = 48_000;
const BOARD_BAND_HZ: f64 = 400.0;
const MICROPHONE_M: [f64;3] = [0.675,1.0,1.0];

pub fn frames(seconds: &str) -> Result<usize,String> {
    let seconds:f64 = seconds.parse().map_err(|_|"invalid render duration")?;
    if !seconds.is_finite() || !(0.05..=60.0).contains(&seconds) {
        return Err("render duration must be finite and in 0.05..60 seconds".into());
    }
    Ok((seconds*f64::from(RATE)).round() as usize)
}
fn prepare(text:&str, scale:Vec<Course>) -> Result<(Instrument,PreparedBoard),String> {
    let keys:Vec<_> = scale.iter().map(|c|c.midi).collect();
    let board=BoardGeometry::read(text)?.prepare(&keys,BOARD_BAND_HZ)?;
    let materials=scale.iter().map(steinway_scale::hammer_material).collect::<Result<Vec<_>,_>>()?;
    let piano=Instrument::new_with_course_shanks(scale,&board.modes,RATE,4,24,true,
        materials,ShankGeometry::published())?;
    Ok((piano,board))
}

pub struct Rendered {
    pub wav:Vec<u8>,
    pub peak_pa:f64,
    pub report:String,
}
fn encode(piano:Instrument, board:&PreparedBoard, score:Performance, count:usize) -> Result<Rendered,String> {
    let mut stream=AudioStream::new(piano,score,Some(&board.surface),MICROPHONE_M,
        fs_bem::helmholtz::Medium::air(),1.0)?;
    let mut pressure=vec![0.0;count];
    for block in pressure.chunks_mut(256) { stream.render_block(block).map_err(|e|e.to_string())?; }
    let peak_pa=pressure.iter().fold(0.0_f64,|p,&x|p.max(x.abs()));
    let (wav,clips)=fs_couple::pcm_wav::encode_pcm16_wav(&pressure,RATE,2.0).map_err(|e|e.to_string())?;
    let piano=stream.instrument();
    let losses=piano.accounting.dissipated_j();
    let report=format!("{} frames at {} Hz; peak {:.9e} Pa, {} clips at 2 Pa PCM full scale; no normalization.\n\
        Panel/ribs/bridges: {:.9} m^2, {:.9} kg; {} modes through {} Hz.\n\
        Input {:.12e} J; stored {:.12e} J; component loss {:.12e} J; closure {:.12e} J.\n\
        Source: {}. Flat board, source/estimated materials, raw source string tensions; no measured digital-twin claim.\n\
        Receiver {:?} m, infinite-baffle Rayleigh pressure; no room, lid scattering or radiation backreaction.",
        count,RATE,peak_pa,clips,board.area_m2,board.mass_kg,board.modes.len(),BOARD_BAND_HZ,
        piano.accounting.input_work_j,piano.energy_j(),losses,
        piano.accounting.input_work_j-piano.energy_j()-losses,board.provenance,MICROPHONE_M);
    Ok(Rendered{wav,peak_pa,report})
}

/// All 88 source courses remain in the coupled system, including silent keys.
/// MIDI is a gesture schedule only; its pitch numbers do not retune strings.
/// Default MIDI mapping: channel 1, velocity 127 -> 4.5 m/s, switch sustain.
pub fn render(text:&str, count:usize, midi:Option<&str>) -> Result<Rendered,String> {
    if !(2_400..=2_880_000).contains(&count) { return Err("render frame budget exceeded".into()); }
    let scale=steinway_scale::courses()?;
    let keys:Vec<_> = scale.iter().map(|c|c.midi).collect();
    // Admit the schedule before doing the expensive modal preparation.
    let score=match midi {
        Some(path)=>{
            let parsed=performance::midi::load(path,&keys,RATE,count as u64,
                performance::midi::Mapping::default())?;
            println!("MIDI channel 1: {} launches, {} other-channel messages ignored; uncalibrated velocity 127 -> 4.5 m/s",
                parsed.report.selected_note_ons,parsed.report.other_channel_messages);
            Performance::from_events(parsed.events,&keys,count as u64)?
        }
        None=>Performance::demonstration(&keys,RATE,count as u64,Some(69),Some(2.0))?,
    };
    let (piano,board)=prepare(text,scale)?;
    encode(piano,&board,score,count)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn panel()->String {
        let mut text=String::from("frankensim-board-geometry-si-v1\nsource,estimated,synthetic mesh-to-source-instrument regression\n");
        for (i,p) in [[0.,0.],[1.,0.],[1.,1.],[0.,1.],[0.5,0.5]].iter().enumerate() {
            text.push_str(&format!("node,{i},{},{}\n",p[0],p[1]));
        }
        for (i,t) in [[0,1,4],[1,2,4],[2,3,4],[3,0,4]].iter().enumerate() {
            text.push_str(&format!("triangle,{i},{},{},{},0.008,450,1e10,8e8,0.3,6e8,0\n",t[0],t[1],t[2]));
        }
        text.push_str("support,simply_supported\nfixed,0\nfixed,1\nfixed,2\nfixed,3\ndamping,0.01\npretension,0\nbridge,69,0,0,0,1\n");
        text
    }
    #[test]
    fn duration_bounds_precede_files_and_allocations() {
        for seconds in ["NaN","inf","-1","0","0.049","61","nonsense"] { assert!(frames(seconds).is_err()); }
        assert_eq!(frames("0.05").unwrap(),2400); assert_eq!(frames("60").unwrap(),2_880_000);
        assert!(render("not a board",1,None).is_err());
    }
    #[test]
    fn imported_panel_drives_source_felt_shanks_and_physical_pressure_deterministically() {
        // One real source course keeps this regression cheap; the public render
        // retains every source course, and refuses this deliberately incomplete board.
        let (obj,spec)=super::super::mesh_import::export(&panel()).unwrap();
        let imported=super::super::mesh_import::import(&obj,&spec).unwrap();
        let run=|| {
            let mut scale=steinway_scale::courses().unwrap(); scale.retain(|c|c.midi==69);
            let (piano,board)=prepare(&imported.fsb,scale).unwrap();
            assert!((board.mass_kg-3.6).abs()<1e-10);
            let score=Performance::read("sample,event,key,value\n0,note_on,69,2\n1024,note_off,69,0\n",&[69],2400).unwrap();
            encode(piano,&board,score,2400).unwrap()
        };
        let a=run(); let b=run();
        assert!(a.peak_pa.is_finite() && a.peak_pa>0.0);
        assert_eq!(&a.wav[..4],b"RIFF"); assert_eq!(a.wav,b.wav); assert_eq!(a.report,b.report);
        assert!(prepare(&imported.fsb,steinway_scale::courses().unwrap()).is_err());
    }
}
