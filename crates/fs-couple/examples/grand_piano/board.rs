//! Board modal input is a physical port table, not a bank of assigned pitches.
//! A measured or fs-plate-derived table can replace the explicitly authored
//! demonstration without changing the runtime or string material/geometry.
use super::linear::{BoardMode, MAX_BOARD_MODES};
use fs_math::det;
use std::f64::consts::PI;

pub const HEADER: &str = "mode,frequency_hz,damping_ratio,volume_m2_per_sqrt_kg,midi,bridge_per_sqrt_kg";
/// `HEADER` plus the gradient (dw/dx, dw/dy) of the mass-normalized mode shape at the key's bridge station.
pub const HEADER_WITH_SLOPES: &str = "mode,frequency_hz,damping_ratio,volume_m2_per_sqrt_kg,midi,bridge_per_sqrt_kg,slope_x_per_m_sqrt_kg,slope_y_per_m_sqrt_kg";

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
    // The six-column table, or the same table with the two bridge-station slope columns that
    // `write_with_slopes` adds. The slopes are validated and then dropped: no runtime path uses them.
    let mut columns=0;
    for (line,raw) in text.lines().enumerate() {
        let raw=raw.trim();if raw.is_empty()||raw.starts_with('#'){continue;}
        if columns==0 {
            columns=match raw{HEADER=>6,HEADER_WITH_SLOPES=>8,_=>return Err("expected unit-explicit board header".into())};
            continue;
        }
        let f:Vec<&str>=raw.split(',').map(str::trim).collect();
        if f.len()!=columns{return Err(format!("board line {}: expected {columns} columns",line+1));}
        let mode=f[0].parse::<usize>().map_err(|_|"invalid board mode index")?;
        let key=f[4].parse::<u8>().map_err(|_|"invalid board MIDI key")?;
        if mode>=MAX_BOARD_MODES || !(21..=108).contains(&key){return Err("board mode/key out of range".into());}
        let parse=|i:usize|f[i].parse::<f64>().map_err(|_|format!("board line {}: invalid scalar",line+1));
        let (hz,zeta,volume,g)=(parse(1)?,parse(2)?,parse(3)?,parse(5)?);
        if [hz,zeta,volume,g].iter().any(|x|!x.is_finite())||hz<=0.0||zeta<0.0{return Err("invalid board scalar".into());}
        if columns==8&&!(parse(6)?.is_finite()&&parse(7)?.is_finite()){return Err("invalid board slope scalar".into());}
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

/// `write`, with `slopes[mode][key-21] = [dw/dx, dw/dy]` appended to every row (see `HEADER_WITH_SLOPES`).
pub fn write_with_slopes(modes:&[BoardMode],slopes:&[[[f64;2];88]])->Result<String,String> {
    if slopes.len()!=modes.len(){return Err("one slope table per board mode is required".into());}
    let mut text=format!("{HEADER_WITH_SLOPES}\n");
    for (i,b) in modes.iter().enumerate(){for (k,g) in b.bridge.iter().enumerate(){
        let [sx,sy]=slopes[i][k];
        text.push_str(&format!("{i},{:.17e},{:.17e},{:.17e},{},{:.17e},{:.17e},{:.17e}\n",b.frequency_hz,b.damping_ratio,b.volume,k+21,g,sx,sy));
    }}
    Ok(text)
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
    fn slope_table(count:usize)->(Vec<BoardMode>,Vec<[[f64;2];88]>) {
        let modes=read(&table(count),&[69]).unwrap();
        let slopes=(0..count).map(|m|{let mut s=[[0.0;2];88];s[48]=[0.25*m as f64+0.5,-1.5*m as f64];s}).collect();
        (modes,slopes)
    }
    #[test]
    fn slope_columns_round_trip_and_leave_the_runtime_table_unchanged() {
        let (modes,slopes)=slope_table(5);
        let text=write_with_slopes(&modes,&slopes).unwrap();
        assert!(text.starts_with(HEADER_WITH_SLOPES));
        let row=text.lines().find(|l|l.split(',').nth(4)==Some("69")&&l.starts_with("3,")).unwrap();
        let f:Vec<f64>=row.split(',').map(|x|x.parse().unwrap()).collect();
        assert_eq!((f[6],f[7]),(1.25,-4.5));
        let back=read(&text,&[69]).unwrap();
        assert_eq!(back.len(),modes.len());
        for (a,b) in modes.iter().zip(&back){assert_eq!(a.bridge,b.bridge);assert_eq!(a.frequency_hz,b.frequency_hz);}
        // The six-column table still reads as before.
        assert_eq!(read(&write(&modes),&[69]).unwrap().len(),modes.len());
    }
    #[test]
    fn slope_columns_are_checked_not_trusted() {
        let (modes,slopes)=slope_table(2);
        let good=write_with_slopes(&modes,&slopes).unwrap();
        // A six-column row under the slope header, and a slope header over six columns, both refuse.
        let short:String=good.lines().map(|l|{let f:Vec<&str>=l.split(',').collect();if f[0]=="mode"{l.to_owned()}else{f[..6].join(",")}}).collect::<Vec<_>>().join("\n");
        assert!(read(&short,&[69]).is_err());
        let wide=write(&modes).replacen(HEADER,HEADER_WITH_SLOPES,1);
        assert!(read(&wide,&[69]).is_err());
        // Poison the x slope of the one row that carries a nonzero slope (mode 1, key 69) and prove it was edited.
        let mut edited=0;
        let poisoned:String=good.lines().map(|l|{
            let mut f:Vec<String>=l.split(',').map(str::to_owned).collect();
            if f[0]=="1"&&f[4]=="69"{f[6]="NaN".into();edited+=1;}
            f.join(",")
        }).collect::<Vec<_>>().join("\n");
        assert_eq!(edited,1);
        assert!(read(&good,&[69]).is_ok());
        assert!(read(&poisoned,&[69]).is_err());
        assert!(write_with_slopes(&modes,&slopes[..1]).is_err());
    }
    #[test]
    fn larger_capacity_does_not_admit_missing_measurements_or_truncate_overbudget_input() {
        assert!(read(&table(MAX_BOARD_MODES),&[60,69]).is_err());
        assert!(read(&table(MAX_BOARD_MODES+1),&[69]).is_err());
    }
}
