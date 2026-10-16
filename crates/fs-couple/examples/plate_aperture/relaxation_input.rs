//! A supplied homogeneous Maxwell spectrum for the physical example.
//! The library supports complete regional maps; this CLI covers the whole plate.
use fs_couple::bernoulli_aperture::dynamic::relaxation::{PlateRelaxationSpec,PlateRelaxationRegion,InitialApertureMemory};
use fs_material::visco::GeneralizedMaxwell;
use std::io::Read;

pub fn load(path:&str, triangles:usize)->Result<(PlateRelaxationSpec,InitialApertureMemory),String> {
    let mut bytes=Vec::new();std::fs::File::open(path).map_err(|e|e.to_string())?
        .take(65537).read_to_end(&mut bytes).map_err(|e|e.to_string())?;
    let text=std::str::from_utf8(&bytes).map_err(|e|e.to_string())?;
    decode(text,triangles,path)
}
fn decode(text:&str,triangles:usize,source:&str)->Result<(PlateRelaxationSpec,InitialApertureMemory),String> {
    if text.len()>65536 {return Err("relaxation input exceeds 64 KiB".into());}
    let mut rows=text.lines().map(str::trim).filter(|s|!s.is_empty() && !s.starts_with('#'));
    if rows.next()!=Some("frankensim-plate-relaxation-v1") {return Err("expected frankensim-plate-relaxation-v1".into());}
    let mut read=|name:&str,count:usize|->Result<Vec<&str>,String> {
        let fields:Vec<_>=rows.next().ok_or_else(||format!("missing {name}"))?.split_whitespace().collect();
        if fields.len()!=count+1 || fields[0]!=name {return Err(format!("expected {name} with {count} fields"));}
        Ok(fields[1..].to_vec())
    };
    let scalar=|s:&str|->Result<f64,String> {
        let x:f64=s.parse().map_err(|_|"invalid scalar".to_string())?;
        if !x.is_finite() {return Err("material input must be finite".into());}Ok(x)
    };
    let e=scalar(read("equilibrium_pa",1)?[0])?;
    let nu=scalar(read("poisson_ratio",1)?[0])?;
    let b=read("band_hz",2)?;let band=(scalar(b[0])?,scalar(b[1])?);
    let max_dt=scalar(read("max_dt_over_tau",1)?[0])?;
    let max_angle=scalar(read("max_angular_step",1)?[0])?;
    let initial=match read("initial",1)?[0] {
        "relaxed"=>InitialApertureMemory::Relaxed,"unrelaxed"=>InitialApertureMemory::Unrelaxed,
        _=>return Err("initial must be explicitly relaxed or unrelaxed".into()),
    };
    let count:usize=read("branches",1)?[0].parse().map_err(|_|"invalid branch count".to_string())?;
    if count>64 {return Err("at most 64 supplied branches".into());}
    let mut terms=Vec::with_capacity(count);
    for _ in 0..count {let p=read("branch",2)?;terms.push((scalar(p[0])?,scalar(p[1])?));}
    if rows.next().is_some() {return Err("unexpected trailing material record".into());}
    let material=GeneralizedMaxwell::new(e,terms).map_err(|e|e.to_string())?;
    Ok((PlateRelaxationSpec {regions:vec![PlateRelaxationRegion {
        triangles:(0..triangles).collect(),material,poisson_ratio:nu,band_hz:band,
        provenance:format!("supplied spectrum file {source}; no inferred material identification"),
    }],max_branches:64,max_dt_over_tau:max_dt,max_angular_step:max_angle},initial))
}
#[cfg(test)]
mod tests {
    use super::*;
    const INPUT:&str=include_str!("illustrative-relaxation.txt");
    #[test]
    fn supplied_branches_and_initial_history_are_not_inferred() {
        let (s,initial)=decode(INPUT,24,"test").unwrap();
        assert_eq!(s.regions[0].material.e_inf,4e9);assert_eq!(s.regions[0].triangles.len(),24);
        assert_eq!(s.regions[0].material.terms,vec![(2e9,0.001),(1e9,0.01)]);
        assert!(matches!(initial,InitialApertureMemory::Relaxed));
        for input in [INPUT.replace("branches 2","branches 18446744073709551615"),
            INPUT.replace("initial relaxed","initial guessed"),INPUT.replace("4000000000","NaN"),
            INPUT.replace("branch 1000000000 0.01",""),format!("{INPUT}\nignored")]
        {assert!(decode(&input,24,"test").is_err());}
    }
}
