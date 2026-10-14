//! Example-only ingress. All head interpolation and source admission stay in fs-airflow.
use std::collections::BTreeMap;
use fs_airflow::SourceProvenance;
use fs_airflow::graph::thermal::medium::series::enthalpy::coupled::feedback::pressure::parallel::pump::{PumpHeadCurve, PumpHeadPoint};
use fs_qty::{Length, VolumetricFlowRate};
use super::Result;

pub struct PumpTable { pub curve: PumpHeadCurve, pub license: String }
impl PumpTable {
    pub fn parse(text: &str) -> Result<Self> {
        if text.len()>128*1024 { return Err("pump table exceeds 128 KiB".into()); }
        let mut headers=BTreeMap::new();let mut points=Vec::new();
        for (i,line) in text.lines().enumerate() {
            let line=line.split('#').next().unwrap_or("").trim();
            if line.is_empty() { continue; }
            if let Some((key,value))=line.split_once('=') {
                let (key,value)=(key.trim(),value.trim());
                if !matches!(key,"source"|"identifier"|"license"|"gravity_m_s2"|"minimum_nu_m2_s"|"maximum_nu_m2_s")
                    || value.is_empty() || headers.insert(key,value).is_some()
                { return Err(format!("invalid or duplicate pump header on line {}",i+1).into()); }
            } else {
                if points.len()==1024 { return Err("pump table exceeds 1024 samples".into()); }
                let row=line.split_whitespace().map(str::parse::<f64>)
                    .collect::<std::result::Result<Vec<_>,_>>()?;
                if row.len()!=2 { return Err(format!("pump line {} needs Q_m3_s H_m",i+1).into()); }
                points.push(PumpHeadPoint { flow: VolumetricFlowRate::new(row[0]),head: Length::new(row[1]) });
            }
        }
        let get=|key| headers.get(key).copied().ok_or_else(||format!("missing pump {key}= header"));
        let curve=PumpHeadCurve::try_new(points,get("gravity_m_s2")?.parse()?,
            (get("minimum_nu_m2_s")?.parse()?,get("maximum_nu_m2_s")?.parse()?),
            SourceProvenance::new(get("source")?,get("identifier")?))?;
        Ok(Self { curve,license:get("license")?.to_owned() })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn source_curve_ingress_retains_license_and_refuses_bad_units_rows_or_domains() {
        let text=include_str!("pump.txt");
        let source=PumpTable::parse(text).unwrap();
        assert_eq!(source.license,"CC0-1.0");
        assert_eq!(source.curve.points().len(),2);
        for bad in [text.replace("gravity_m_s2=10","gravity_m_s2=NaN"),
            text.replace("minimum_nu_m2_s=0.000001","minimum_nu_m2_s=0.001"),
            text.replace("0.00002 0","0 0"),text.replace("license=CC0-1.0", ""),
            format!("{text}\nunknown=1"),format!("{text}\n0.1 1 2")]
        { assert!(PumpTable::parse(&bad).is_err()); }
    }
}
