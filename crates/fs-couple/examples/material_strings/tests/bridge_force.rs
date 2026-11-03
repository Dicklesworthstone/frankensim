//! `--bridge-force`: the same WAV bytes, with each accepted sample's connection reactions beside them.
//! Gauntlet tier G2/G3: an executable fixture, and a metamorphic identity (the read-out cannot change the solve).
use super::*;

const SAMPLES:usize=16_384;
const RATE_HZ:f64=48_000.0;

fn bowed()->input::Spec {
    // pluck-and-bow: the bow lands at 10 ms and lifts at 40 ms, so the rest of the run is the string ringing down.
    let text=include_str!("../pluck-and-bow.material-strings").replace("48000 4096 0.01 128",&format!("48000 {SAMPLES} 0.01 128"));
    input::Spec::parse(text.as_bytes()).unwrap()
}

/// Power of `x` at `f_hz`, by direct summation over the whole window (Hann-weighted); the windows here are a few
/// thousand samples and a few hundred frequencies, so no FFT is needed.
fn line_power(x:&[f64],f_hz:f64)->f64 {
    let n=x.len() as f64;
    let (mut re,mut im)=(0.0,0.0);
    for (k,v) in x.iter().enumerate() {
        let w=0.5-0.5*(2.0*std::f64::consts::PI*k as f64/(n-1.0)).cos();
        let phase=2.0*std::f64::consts::PI*f_hz*k as f64/RATE_HZ;
        re+=w*v*phase.cos();im-=w*v*phase.sin();
    }
    re*re+im*im
}

#[test]
fn a_bridge_force_render_writes_the_same_wav_and_the_strings_own_force_at_its_f0() {
    scoped(|cx,gate| {
        let mut plain=compile(cx,gate,bowed(),7800.0);
        let (wav,summary)=plain.write(Cursor::new(Vec::<u8>::new()),gate).unwrap();
        let mut forced=compile(cx,gate,bowed(),7800.0);
        let mut forces=Vec::<u8>::new();
        let (forced_wav,forced_summary)=forced.write_bridge_forces(Cursor::new(Vec::<u8>::new()),&mut forces,false,gate).unwrap();
        // The read-out cannot change the solve: same PCM bytes, same summary, same final mechanical state.
        assert_eq!(forced_wav.into_inner(),wav.into_inner());
        assert_eq!(forced_summary,summary);
        assert_eq!(forced.performance.bank().system().last_frame(),plain.performance.bank().system().last_frame());
        let text=String::from_utf8(forces).unwrap();
        let mut lines=text.lines();
        let header:Vec<&str>=lines.next().unwrap().split('\t').collect();
        assert_eq!(header[0],"sample");
        assert!(header.len()>=2,"one string means at least one connection column");
        for (k,name) in header[1..].iter().enumerate() { assert_eq!(*name,format!("connection_{k}_n")); }
        let rows:Vec<Vec<f64>>=lines.map(|line|line.split('\t').map(|v|v.parse::<f64>().unwrap()).collect()).collect();
        assert_eq!(rows.len(),SAMPLES);
        for (k,row) in rows.iter().enumerate() {
            assert_eq!(row.len(),header.len());
            assert_eq!(row[0],k as f64+1.0);
            assert!(row.iter().all(|v|v.is_finite()));
        }
        // This card has one string, so connection 0 is that string's own spring/dashpot bridge.
        let bridge:Vec<f64>=rows.iter().map(|row|row[1]).collect();
        assert!(bridge.iter().any(|f|f.abs()>1e-9),"the bridge force must not be identically zero");
        // After the bow lifts (sample ~1920) the string rings freely. Its fundamental from the card: the vibrating length is the
        // SPEAKING length, the string row's span times its bridge station (MEASURED_STRINGS.md): 0.5 m x 0.8 = 0.4 m; 0.3 mm
        // radius, 20 N, 7800 kg/m^3 (the density `compile` gives the wire): f0 = sqrt(T/(rho pi r^2))/(2 L) = 119.0 Hz.
        let f0=(20.0/(7800.0*std::f64::consts::PI*0.0003_f64.powi(2))).sqrt()/(2.0*0.5*0.8);
        let ring=&bridge[4096..];
        let grid:Vec<f64>=(0..=800).map(|k|40.0+0.5*k as f64).collect();
        let power:Vec<f64>=grid.iter().map(|&f|line_power(ring,f)).collect();
        let (best,_)=grid.iter().zip(&power).max_by(|a,b|a.1.partial_cmp(b.1).unwrap()).unwrap();
        eprintln!("bridge-force strongest line {best} Hz; string f0 from the card {f0} Hz");
        assert!((best-f0).abs()<0.05*f0,"strongest bridge-force line {best} Hz is not the string's f0 {f0} Hz");
    });
}
