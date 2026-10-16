//! Board modal input is a physical port table, not a bank of assigned pitches.
//! A measured or fs-plate-derived table can replace the explicitly authored
//! demonstration without changing the runtime or string material/geometry.
use super::linear::{BoardMode, MAX_BOARD_MODES};
use fs_math::det;
use std::f64::consts::PI;

pub const HEADER: &str = "mode,frequency_hz,damping_ratio,volume_m2_per_sqrt_kg,midi,bridge_per_sqrt_kg";

/// Temporary illustrative board, NOT a measured/ribbed/crowned Steinway board.
/// Mode shapes carry kg normalization. The FE geometry producer is a separate
/// cold front door; an arbitrary 1 kg board is never silently assumed.
pub fn demonstration() -> Vec<BoardMode> {
    [(95.0,1.0,1.0),(183.0,2.0,1.0),(297.0,1.0,2.0),(416.0,3.0,1.0)]
        .iter().map(|&(hz,m,n)| {
            let mass_root=det::sqrt(8.0/4.0); // authored 8 kg rectangular surrogate
            let mut bridge=[0.0;88];
            for (k,g) in bridge.iter_mut().enumerate() {
                let x=0.12+0.75*k as f64/87.0;
                let y=0.62-0.28*k as f64/87.0;
                *g=det::sin(m*PI*x)*det::sin(n*PI*y)/mass_root;
            }
            let volume=2.4*(1.0-det::cos(m*PI))*(1.0-det::cos(n*PI))/(m*n*PI*PI*mass_root);
            BoardMode {frequency_hz:hz,damping_ratio:0.015,bridge,volume}
        }).collect()
}

/// Every retained mode must specify the bridge for every admitted key. Missing
/// values are refused, not filled with a convenient invented spatial shape.
pub fn read(text:&str,keys:&[u8])->Result<Vec<BoardMode>,String> {
    if keys.is_empty()||keys.iter().any(|k|!(21..=108).contains(k)){return Err("invalid admitted board key set".into());}
    let mut modes:Vec<BoardMode>=Vec::new();
    let mut seen:Vec<[bool;88]>=Vec::new();
    let mut header=false;
    for (line,raw) in text.lines().enumerate() {
        let raw=raw.trim();if raw.is_empty()||raw.starts_with('#'){continue;}
        if !header {if raw!=HEADER{return Err("expected unit-explicit board header".into());}header=true;continue;}
        let f:Vec<&str>=raw.split(',').map(str::trim).collect();
        if f.len()!=6{return Err(format!("board line {}: expected six columns",line+1));}
        let mode=f[0].parse::<usize>().map_err(|_|"invalid board mode index")?;
        let key=f[4].parse::<u8>().map_err(|_|"invalid board MIDI key")?;
        if mode>=MAX_BOARD_MODES || !(21..=108).contains(&key){return Err("board mode/key out of range".into());}
        let parse=|i:usize|f[i].parse::<f64>().map_err(|_|format!("board line {}: invalid scalar",line+1));
        let (hz,zeta,volume,g)=(parse(1)?,parse(2)?,parse(3)?,parse(5)?);
        if [hz,zeta,volume,g].iter().any(|x|!x.is_finite())||hz<=0.0||zeta<0.0{return Err("invalid board scalar".into());}
        if mode>modes.len(){return Err("board modes must appear in contiguous index order".into());}
        if mode==modes.len(){modes.push(BoardMode{frequency_hz:hz,damping_ratio:zeta,volume,bridge:[0.0;88]});seen.push([false;88]);}
        let b=&mut modes[mode];
        if b.frequency_hz!=hz||b.damping_ratio!=zeta||b.volume!=volume{return Err("inconsistent repeated board mode properties".into());}
        let k=usize::from(key-21);
        if seen[mode][k]{return Err("duplicate board mode/key row".into());}
        b.bridge[k]=g;seen[mode][k]=true;
    }
    if modes.is_empty(){return Err("empty soundboard".into());}
    for (m,rows) in seen.iter().enumerate(){for &key in keys{
        if !rows[usize::from(key-21)]{return Err(format!("board mode {m} has no bridge shape for key {key}"));}
    }}
    Ok(modes)
}

pub fn write(modes:&[BoardMode])->String {
    let mut text=format!("{HEADER}\n");
    for (i,b) in modes.iter().enumerate(){for (k,g) in b.bridge.iter().enumerate(){
        text.push_str(&format!("{i},{:.17e},{:.17e},{:.17e},{},{:.17e}\n",b.frequency_hz,b.damping_ratio,b.volume,k+21,g));
    }}
    text
}

#[cfg(test)]
mod tests {
    use super::*;
    fn table(count:usize)->String {
        let mut text=format!("{HEADER}\n");
        for mode in 0..count {
            text.push_str(&format!("{mode},{},0.01,0.1,69,{}\n",90.0+mode as f64*10.0,mode as f64*0.001));
        }
        text
    }
    #[test]
    fn wider_board_tables_keep_every_supplied_mode_and_bridge_coefficient() {
        for count in [33,64,MAX_BOARD_MODES] {
            let modes=read(&table(count),&[69]).unwrap();
            assert_eq!(modes.len(),count);
            for (i,m) in modes.iter().enumerate() {
                assert_eq!(m.frequency_hz,90.0+i as f64*10.0);
                assert_eq!(m.bridge[48],i as f64*0.001);
            }
            let roundtrip=read(&write(&modes),&[69]).unwrap();
            assert_eq!(roundtrip.len(),count);
            assert_eq!(roundtrip[count-1].bridge,modes[count-1].bridge);
        }
    }
    #[test]
    fn larger_capacity_does_not_admit_missing_measurements_or_truncate_overbudget_input() {
        assert!(read(&table(MAX_BOARD_MODES),&[60,69]).is_err());
        assert!(read(&table(MAX_BOARD_MODES+1),&[69]).is_err());
    }
}
