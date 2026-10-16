//! Bounded data ingress. The material card and existing source owners do the
//! interpolation, constitutive admission and caloric construction.
use std::io::Read;
use fs_conduction::lumped::LumpedThermalTransport;
use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
    NormalizedMaterialCardPack, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::state_point::{DENSITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS,
    SPECIFIC_HEAT_CAPACITY_PROPERTY, THERMAL_CONDUCTIVITY_DIMS, THERMAL_CONDUCTIVITY_PROPERTY};
use fs_material::thermovisco::Limits;
use fs_material::thermovisco::source::{ResolvedThermoMaxwell, ThermoMaxwellSourceQuery};
use fs_qty::{Density, Dims, Pressure, Temperature, Time};
use super::{Result, VOLUME_M3, surface_area};
const ONE: Dims = Dims([0;6]);
const MAX_BYTES: u64 = 16*1024*1024;

pub(super) struct Source {
    pub card: MaterialCard,
    pub mechanical: ResolvedThermoMaxwell,
    pub caloric: ResolvedHeatCapacityCurve,
    pub transport: LumpedThermalTransport,
}
impl Source {
    pub fn from_card(card: MaterialCard, branches: usize, low: f64, high: f64,
        gate: &CancelGate) -> Result<Self> {
        if branches > 30 { return Err("at most 30 explicitly declared branches are supported".into()); }
        let lower=QueryPoint::new().with("T",low)?;
        let upper=QueryPoint::new().with("T",high)?;
        let branch_keys:Vec<_>=(0..branches).map(|i| (
            PropertyKey::new(format!("relaxing_young_modulus_{i}"),Pressure::DIMS),
            PropertyKey::new(format!("relaxation_time_{i}"),Time::DIMS))).collect();
        let mechanical=ResolvedThermoMaxwell::try_new(&card,ThermoMaxwellSourceQuery {
            equilibrium_modulus:&PropertyKey::new("equilibrium_young_modulus",Pressure::DIMS),
            branches:&branch_keys,modulus_scale:&PropertyKey::new("modulus_scale",ONE),
            relaxation_scale:&PropertyKey::new("relaxation_scale",ONE),
            free_strain:&PropertyKey::new("free_thermal_strain",ONE),temperature_axis:"T",
            lower:&lower,upper:&upper,selection:SelectionPolicy::SingleClaimOnly,
        },Limits {max_strain:0.04,max_strain_rate_s:10.0},|| gate.is_requested())?;
        let caloric=ResolvedHeatCapacityCurve::try_new(&card,HeatCapacitySourceQuery {
            heat_capacity:&PropertyKey::new(SPECIFIC_HEAT_CAPACITY_PROPERTY,SPECIFIC_HEAT_CAPACITY_DIMS),
            density:&PropertyKey::new(DENSITY_PROPERTY,Density::DIMS),temperature_axis:"T",
            lower:&lower,upper:&upper,selection:SelectionPolicy::SingleClaimOnly,
        },SolidLiquidPhase::Solid,0.0,1e-7)?;
        let transport=LumpedThermalTransport::from_material_card_without_radiation(&card,
            &PropertyKey::new(THERMAL_CONDUCTIVITY_PROPERTY,THERMAL_CONDUCTIVITY_DIMS),"T",
            &[lower,upper],SelectionPolicy::SingleClaimOnly)?;
        // Initial feedback admission checks the common card/context and the
        // existing fixed-density, constant-energy, solid-only physical limits.
        mechanical.caloric_point(&caloric,300.0,0.0,VOLUME_M3,surface_area())?;
        if gate.is_requested() { return Err("cancelled material ingress".into()); }
        Ok(Self {card,mechanical,caloric,transport})
    }
    pub fn pack(path:&str, branches:usize, gate:&CancelGate) -> Result<Self> {
        let bytes=read(path)?;
        let pack=NormalizedMaterialCardPack::from_bytes(&bytes)?;
        Self::from_card(pack.card().clone(),branches,280.0,320.0,gate)
    }
    pub fn table(path:&str, gate:&CancelGate) -> Result<Self> {
        Self::table_text(std::str::from_utf8(&read(path)?)?,gate)
    }
    pub fn table_text(text:&str, gate:&CancelGate) -> Result<Self> {
        if text.len() as u64 > MAX_BYTES { return Err("material table exceeds 16 MiB".into()); }
        let mut rows=text.lines().filter_map(|line| {
            let line=line.split('#').next().unwrap_or("").trim();
            (!line.is_empty()).then_some(line)
        });
        if rows.next()!=Some("frankensim.thermovisco.v1") { return Err("unsupported material-table header".into()); }
        let mut equilibrium=None; let mut branches=Vec::new(); let mut thermal=Vec::<[f64;7]>::new();
        for line in rows {
            if gate.is_requested() { return Err("cancelled material table".into()); }
            let mut words=line.split_whitespace();
            let kind=words.next().ok_or("empty material row")?;
            let mut values=Vec::new();
            for word in words {
                if values.len()==7 { return Err("too many material row values".into()); }
                let value:f64=word.parse()?;
                if !value.is_finite() { return Err("material values must be finite".into()); }
                values.push(value);
            }
            match (kind,values.as_slice()) {
                ("equilibrium",[e]) if equilibrium.is_none() => equilibrium=Some(*e),
                ("branch",[e,tau]) if branches.len()<30 => branches.push((*e,*tau)),
                ("state",[t,scale,shift,free,cp,rho,k]) if thermal.len()<1024 => {
                    if *t<=0.0 || thermal.last().is_some_and(|old| old[0]>=*t) {
                        return Err("thermal rows need strictly increasing positive kelvin".into());
                    }
                    thermal.push([*t,*scale,*shift,*free,*cp,*rho,*k]);
                }
                _ => return Err("unknown, repeated, malformed or over-budget material row".into()),
            }
        }
        let equilibrium=equilibrium.ok_or("equilibrium modulus is required")?;
        if thermal.len()<2 { return Err("at least two thermal state rows are required".into()); }
        let low=thermal[0][0]; let high=thermal[thermal.len()-1][0];
        let mut identity=fs_blake3::DomainHasher::new("org.frankensim.example.thermovisco-table.v1");
        identity.update(text.as_bytes());
        let provenance=Provenance { source:format!("declared thermovisco input {} (not independent measurements)",identity.finalize().to_hex()),
            license:"input terms unspecified".into(),artifact:None };
        let mut claims=ClaimSet::new();
        let mut constant=|name:String,dims: Dims,value:f64| -> Result<()> {
            claims.insert_claim(PropertyClaim {key:PropertyKey::new(name,dims),
                value:PropertyValue::Scalar {value,dims},validity:ValidityDomain::unconstrained().with("T",low,high),
                uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::ConstantWithinValidity,
                observations:vec![],provenance:provenance.clone()})?;
            Ok(())
        };
        constant("equilibrium_young_modulus".into(),Pressure::DIMS,equilibrium)?;
        for (i,&(e,tau)) in branches.iter().enumerate() {
            constant(format!("relaxing_young_modulus_{i}"),Pressure::DIMS,e)?;
            constant(format!("relaxation_time_{i}"),Time::DIMS,tau)?;
        }
        for (column,name,dims) in [(1,"modulus_scale",ONE),(2,"relaxation_scale",ONE),
            (3,"free_thermal_strain",ONE),(4,SPECIFIC_HEAT_CAPACITY_PROPERTY,SPECIFIC_HEAT_CAPACITY_DIMS),
            (5,DENSITY_PROPERTY,Density::DIMS),(6,THERMAL_CONDUCTIVITY_PROPERTY,THERMAL_CONDUCTIVITY_DIMS)] {
            claims.insert_claim(PropertyClaim {key:PropertyKey::new(name,dims),
                value:PropertyValue::Curve {abscissa:"T".into(),abscissa_dims:Temperature::DIMS,
                    knots:thermal.iter().map(|row|(row[0],row[column])).collect(),dims},
                validity:ValidityDomain::unconstrained().with("T",low,high),uncertainty:UncertaintyModel::Unstated,
                interpolation:InterpolationPolicy::LinearInside,observations:vec![],provenance:provenance.clone()})?;
        }
        let card=MaterialCard::assemble(MaterialStateId {chemistry:"declared thermovisco table".into(),
            phase:"solid model assumption".into(),process:"user-supplied coefficients; unqualified".into(),revision:0},claims,vec![])?;
        Self::from_card(card,branches.len(),low,high,gate)
    }
}
fn read(path:&str) -> Result<Vec<u8>> {
    let mut bytes=Vec::new(); std::fs::File::open(path)?.take(MAX_BYTES+1).read_to_end(&mut bytes)?;
    if bytes.len() as u64 > MAX_BYTES { return Err("material file exceeds 16 MiB".into()); }
    Ok(bytes)
}
