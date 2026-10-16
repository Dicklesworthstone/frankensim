//! Example-only electrical-interface data; the production consumer owns no fit.
use super::{midpoint_options, simulate_interfaces};
use fs_lbm::material::electrothermal::magnetic::transport::midpoint::interfaces::{
    InterfaceCircuit, InterfaceMidpointOptions, InterfaceResponse,
};
use fs_material::MaterialError;
use std::{error::Error, io::Read};

pub(super) struct Table { knots: Vec<[f64;3]> }
#[derive(Debug, Clone, PartialEq)]
pub(super) struct FilmSample { temperature_k:f64, bracket:[usize;2] }
impl Table {
    fn parse(text:&str)->Result<Self,Box<dyn Error>> {
        if text.len()>65_536 {return Err("interface table exceeds 64 KiB".into());}
        let mut knots:Vec<[f64;3]>=Vec::new();
        for line in text.lines() {
            let line=line.split('#').next().unwrap_or("").trim();
            if line.is_empty() {continue;}
            if knots.len()==256 {return Err("interface table exceeds 256 knots".into());}
            let row=line.split_whitespace().map(str::parse).collect::<Result<Vec<f64>,_>>()?;
            if row.len()!=3 || !row.iter().all(|v|v.is_finite()) || row[0]<=0.0 || row[1]<0.0 || !(0.0..=1.0).contains(&row[2]) {
                return Err("expected Kelvin, nonnegative contact ohm-m2 and heat fraction in [0,1]".into());
            }
            if knots.last().is_some_and(|p|p[0]>=row[0]) {return Err("interface temperatures must strictly increase".into());}
            knots.push([row[0],row[1],row[2]]);
        }
        if knots.len()<2 {return Err("interface table requires 2..256 knots".into());}
        Ok(Self {knots})
    }
    fn read(path:&str)->Result<Self,Box<dyn Error>> {
        let mut text=String::new();
        std::fs::File::open(path)?.take(65_537).read_to_string(&mut text)?;
        Self::parse(&text)
    }
    pub(super) fn at(&self,temperatures:[f64;2])->Result<InterfaceResponse<FilmSample>,MaterialError> {
        let lo=self.knots[0][0];let hi=self.knots[self.knots.len()-1][0];
        if !temperatures.iter().all(|t|t.is_finite() && *t>=lo && *t<=hi) {
            return Err(MaterialError::State {what:format!("surface temperatures {temperatures:?} K outside interface table support")});
        }
        // This EXAMPLE explicitly declares a film-temperature interpolation.
        // The production API does not silently average the supplied surfaces.
        let t=0.5*temperatures[0]+0.5*temperatures[1];
        let upper=self.knots.partition_point(|p|p[0]<t).max(1).min(self.knots.len()-1);
        let a=self.knots[upper-1];let b=self.knots[upper];let f=(t-a[0])/(b[0]-a[0]);
        Ok(InterfaceResponse {material:FilmSample {temperature_k:t,bracket:[upper-1,upper]},
            resistance_ohm_m2:(1.0-f)*a[1]+f*b[1], heat_fraction_to_first:(1.0-f)*a[2]+f*b[2]})
    }
}
pub(super) fn options()->InterfaceMidpointOptions {
    InterfaceMidpointOptions {material:midpoint_options(), resistance_relative_tolerance:1e-8,
        partition_absolute_tolerance:1e-8,circuit:InterfaceCircuit::Wired}
}
pub(super) fn run(bulk:&super::Table,steps:usize,voltage:f64,b:f64,path:&str)->Result<(),Box<dyn Error>> {
    let interface=Table::read(path)?;
    println!("# synthetic bulk and interface tables; declared film-temperature law; electrical-only fixed faces; coupled midpoint");
    println!("time_s,mean_K,mean_vx_m_s,current_A,joule_W,lorentz_W,min_contact_ohm_m2,max_contact_ohm_m2,interface_relative,coupling_defect_J,nonlinear_trials");
    let totals=simulate_interfaces(bulk,steps,voltage,b,true,Some(&interface),|r| {
        println!("{:.6},{:.9},{:.9},{:.12},{:.9},{:.9},{:.9},{:.9},{:.6e},{:.6e},{}",
            r.time_s,r.mean_k,r.mean_vx,r.current_a,r.joule_w,r.mechanical_w,r.min_contact,r.max_contact,
            r.interface_relative,r.defect_j,r.iterations);
    })?;
    eprintln!("electrode_J={:.9} joule_J={:.9} midpoint_mechanical_J={:.9} stored_J={:.9} wall_J={:.9} coupling_defect_J={:.6e}",
        totals.electrode_j,totals.joule_j,totals.mechanical_j,totals.stored_j,totals.wall_j,totals.defect());
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn declared_interpolation_checks_both_surface_domains_and_all_columns() {
        let table=Table::parse("300 0.01 0.25\n320 0.03 0.75\n").unwrap();
        let r=table.at([305.0,315.0]).unwrap();
        assert!((r.resistance_ohm_m2-0.02).abs()<1e-15);assert_eq!(r.heat_fraction_to_first,0.5);
        assert_eq!(r.material.temperature_k,310.0);assert_eq!(r.material.bracket,[0,1]);
        // An in-range mean does not admit an unsupported individual surface.
        assert!(table.at([290.0,330.0]).is_err());
        for t in ["300 0.01 0.5", "300 -1 0.5\n320 0.03 0.5", "300 0.01 2\n320 0.03 0.5",
            "300 0.01 0.5\n300 0.03 0.5", "300 NaN 0.5\n320 0.03 0.5"] { assert!(Table::parse(t).is_err()); }
        let zeros=Table::parse("300 0 0.5\n320 0 0.5").unwrap();assert_eq!(zeros.at([305.0,315.0]).unwrap().resistance_ohm_m2,0.0);
    }
    #[test]
    fn interface_data_alone_changes_actual_current_heat_and_motion() {
        let bulk=super::super::Table::parse(include_str!("electrical_constant.txt")).unwrap();
        let a=Table::parse(include_str!("interface_constant.txt")).unwrap();
        let b=Table::parse(include_str!("interface_temperature_dependent.txt")).unwrap();
        let mut ra=Vec::new();let mut rb=Vec::new();
        let x=simulate_interfaces(&bulk,40,10.0,0.5,true,Some(&a),|r|ra.push(r)).unwrap();
        let y=simulate_interfaces(&bulk,40,10.0,0.5,true,Some(&b),|r|rb.push(r)).unwrap();
        assert!((x.joule_j-y.joule_j).abs()>1e-6);
        assert!((x.last.mean_vx-y.last.mean_vx).abs()>1e-8);
        assert!((x.last.mean_k-y.last.mean_k).abs()>1e-6);
        assert!(ra[0].current_a>rb[0].current_a);
        for r in ra.iter().chain(&rb) {assert!(r.interface_relative<=options().resistance_relative_tolerance);}
        assert!(ra[20..].iter().any(|r|r.joule_w>0.0 && r.mechanical_w<0.0));
        assert!(x.defect().abs()<1e-7 && y.defect().abs()<1e-7);
    }
    #[test]
    fn zero_resistance_interfaces_match_the_existing_midpoint_electrical_path() {
        let bulk=super::super::Table::parse(include_str!("electrical_constant.txt")).unwrap();
        let zero=Table::parse("290 0 0.75\n330 0 0.75").unwrap();
        let a=simulate_interfaces(&bulk,4,2.0,0.0,true,Some(&zero),|_|{}).unwrap();
        let b=super::super::simulate_mode(&bulk,4,2.0,0.0,true,|_|{}).unwrap();
        assert_eq!(a.last.current_a,b.last.current_a);assert_eq!(a.last.mean_k,b.last.mean_k);
        assert_eq!(a.last.mean_vx,b.last.mean_vx);assert_eq!(a.joule_j,b.joule_j);
    }
}
