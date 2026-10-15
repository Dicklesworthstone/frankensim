//! Input only: all property, enthalpy and expansion mathematics stay shared.
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::state_point::*;
use fs_qty::{Density, Dims, Pressure, QuantitySpec, Temperature};
use super::{Error, invalid};

pub const MAX_BYTES: usize = 65_536;
pub struct Source { pub card: MaterialCard, pub caloric: ResolvedHeatCapacityCurve, pub reference_k: f64 }
pub fn point(t: f64) -> Result<QueryPoint, Error> {
    Ok(QueryPoint::new().with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), t)?
        .with_quantity("pressure", QuantitySpec::dimensional(Pressure::DIMS), 100_000.)?)
}
impl Source {
    pub fn parse(text: &str, label: &str) -> Result<Self, Error> {
        if text.len() > MAX_BYTES { return Err(invalid("material input exceeds 64 KiB")); }
        let mut lines = text.lines().enumerate().filter_map(|(i,line)| {
            let line = line.split('#').next().unwrap_or("").trim();
            (!line.is_empty()).then_some((i+1,line))
        });
        if lines.next().map(|(_,s)|s) != Some("frankensim-thermal-vibration-v1") {
            return Err(invalid("missing thermal vibration material v1 header"));
        }
        let mut rows: Vec<[f64; 6]> = Vec::new();
        for (number,line) in lines {
            if rows.len() == 256 { return Err(invalid("material table exceeds 256 rows")); }
            let fields: Vec<_> = line.split_whitespace().collect();
            if fields.len()!=6 { return Err(invalid(&format!("line {number}: expected T density Cp E nu alpha"))); }
            let mut row = [0.0_f64;6];
            for (v,s) in row.iter_mut().zip(fields) { *v=s.parse()?; }
            if row.iter().any(|v|!v.is_finite()) || [0,1,2,3].iter().any(|&i|row[i]<=0.)
                || row[4]<=-1. || row[4]>=0.5 || rows.last().is_some_and(|r|row[0]<=r[0]) {
                return Err(invalid(&format!("line {number}: finite data, positive T/rho/Cp/E, -1<nu<0.5 and increasing T required")));
            }
            rows.push(row);
        }
        if rows.len()<2 { return Err(invalid("at least two temperature rows required")); }
        let (lower,upper)=(rows[0][0],rows[rows.len()-1][0]);
        let properties=[(DENSITY_PROPERTY,Density::DIMS),
            (SPECIFIC_HEAT_CAPACITY_PROPERTY,SPECIFIC_HEAT_CAPACITY_DIMS),
            (YOUNG_MODULUS_PROPERTY,Pressure::DIMS),(POISSON_RATIO_PROPERTY,Dims([0;6])),
            (LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,INVERSE_TEMPERATURE_DIMS)];
        let mut claims=ClaimSet::new();
        for (column,(key,dims)) in properties.into_iter().enumerate() {
            claims.insert_claim(PropertyClaim { key:PropertyKey::new(key,dims),
                value:PropertyValue::Curve {abscissa:"T".into(),abscissa_dims:Temperature::DIMS,
                    knots:rows.iter().map(|r|(r[0],r[column+1])).collect(),dims},
                validity:ValidityDomain::unconstrained()
                    .with_quantity("T",QuantitySpec::dimensional(Temperature::DIMS),lower,upper)
                    .with_quantity("pressure",QuantitySpec::dimensional(Pressure::DIMS),100_000.,100_000.),
                interpolation:InterpolationPolicy::LinearInside,uncertainty:UncertaintyModel::Unstated,
                observations:vec![],provenance:Provenance {source:format!("caller table {label}; not independently verified"),
                    license:"caller-supplied".into(),artifact:None}})?;
        }
        let card=MaterialCard::assemble(MaterialStateId {chemistry:label.into(),phase:"solid".into(),
            process:"caller-supplied isotropic thermal vibration data".into(),revision:0},claims,vec![])?;
        let caloric=ResolvedHeatCapacityCurve::try_new(&card,HeatCapacitySourceQuery {
            heat_capacity:&PropertyKey::new(SPECIFIC_HEAT_CAPACITY_PROPERTY,SPECIFIC_HEAT_CAPACITY_DIMS),
            density:&PropertyKey::new(DENSITY_PROPERTY,Density::DIMS),temperature_axis:"T",
            lower:&point(lower)?,upper:&point(upper)?,selection:SelectionPolicy::SingleClaimOnly,
        },SolidLiquidPhase::Solid,0.,1e-8)?;
        Ok(Self {card,caloric,reference_k:lower})
    }
}
