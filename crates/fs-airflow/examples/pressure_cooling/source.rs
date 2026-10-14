//! Small example-only ingress: authored SI rows become ordinary material cards.
use std::collections::BTreeMap;

use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, QuantitySpec, Temperature};

use super::Result;

pub struct ThermalTable {
    pub card: MaterialCard,
    pub lower_k: f64,
    pub upper_k: f64,
    pub pressure_pa: f64,
}

// Existing liquid callers retain their entry point and identical card data.
pub type LiquidTable = ThermalTable;

impl ThermalTable {
    pub fn parse(text: &str) -> Result<Self> { Self::parse_phase(text, false) }

    /// Solid T/rho/Cp/k data uses the same bounded parser and material-card owner.
    /// No viscosity is invented or required for a solid thermal calculation.
    pub fn parse_solid(text: &str) -> Result<Self> { Self::parse_phase(text, true) }

    fn parse_phase(text: &str, solid: bool) -> Result<Self> {
        if text.len() > 128 * 1024 {
            return Err("material table exceeds 128 KiB".into());
        }
        let mut headers = BTreeMap::new();
        let mut rows: Vec<Vec<f64>> = Vec::new();
        let columns = if solid { 4 } else { 5 };
        for (index, line) in text.lines().enumerate() {
            let line = line.split('#').next().unwrap_or("").trim();
            if line.is_empty() { continue; }
            if let Some((key, value)) = line.split_once('=') {
                let (key, value) = (key.trim(), value.trim());
                if !matches!(key, "material" | "source" | "license" | "pressure_pa")
                    || value.is_empty() || headers.insert(key, value).is_some()
                {
                    return Err(format!("invalid or duplicate header on line {}", index + 1).into());
                }
            } else {
                if rows.len() == 256 { return Err("material table exceeds 256 rows".into()); }
                let values = line.split_whitespace().map(str::parse::<f64>)
                    .collect::<std::result::Result<Vec<_>, _>>()?;
                if values.len() != columns {
                    let expected = if solid { "T rho Cp k" } else { "T rho Cp k mu" };
                    return Err(format!("line {} needs {expected}", index + 1).into());
                }
                let row = values;
                if row.iter().any(|v| !v.is_finite() || *v <= 0.0)
                    || rows.last().is_some_and(|prior| prior[0] >= row[0])
                {
                    return Err(format!("line {} needs positive SI values and increasing T", index + 1).into());
                }
                rows.push(row);
            }
        }
        let required = |key| headers.get(key).copied()
            .ok_or_else(|| format!("missing {key}= header"));
        let name = required("material")?;
        let source = required("source")?;
        let license = required("license")?;
        let pressure_pa: f64 = required("pressure_pa")?.parse()?;
        if rows.len() < 2 || !pressure_pa.is_finite() || pressure_pa <= 0.0 {
            return Err("at least two rows and a positive absolute source pressure are required".into());
        }
        let lower_k = rows[0][0];
        let upper_k = rows[rows.len() - 1][0];
        let mut claims = ClaimSet::new();
        for (column, property, dims) in [
            (1, "density", Dims([-3, 1, 0, 0, 0, 0])),
            (2, "specific_heat_capacity", Dims([2, 0, -2, -1, 0, 0])),
            (3, "thermal_conductivity", Dims([1, 1, -3, -1, 0, 0])),
            (4, "dynamic_viscosity", Dims([-1, 1, -1, 0, 0, 0])),
        ] {
            if column >= columns { continue; }
            let quantity = if property == "thermal_conductivity" {
                QuantitySpec::semantic(SemanticType::new(QuantityKind::ThermalConductivity, ValueForm::Static))
            } else { QuantitySpec::dimensional(dims) };
            claims.insert_claim(PropertyClaim {
                key: PropertyKey::with_quantity(property, quantity),
                value: PropertyValue::Curve {
                    abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                    knots: rows.iter().map(|row| (row[0], row[column])).collect(), dims,
                },
                validity: ValidityDomain::unconstrained().with("T", lower_k, upper_k)
                    .with("pressure", pressure_pa, pressure_pa),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance { source: source.into(), license: license.into(), artifact: None },
            })?;
        }
        let card = MaterialCard::assemble(MaterialStateId {
            chemistry: name.into(), phase: if solid { "solid" } else { "liquid" }.into(),
            process: if solid { "authored stored-solid example table" }
                else { "authored pressure-cooling example table" }.into(), revision: 0,
        }, claims, vec![])?;
        Ok(Self { card, lower_k, upper_k, pressure_pa })
    }

    pub fn point(&self, temperature_k: f64) -> Result<QueryPoint> {
        Ok(QueryPoint::new().with("T", temperature_k)?.with("pressure", self.pressure_pa)?)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn solid_and_liquid_ingress_share_thermal_claims_without_inventing_viscosity() {
        let liquid = include_str!("constant.txt");
        let solid = liquid.replace(" 0.002", "");
        let a = ThermalTable::parse(liquid).unwrap();
        let b = ThermalTable::parse_solid(&solid).unwrap();
        assert_ne!(a.card.content_hash(), b.card.content_hash());
        for (name, dims) in [
            ("density", Dims([-3,1,0,0,0,0])),
            ("specific_heat_capacity", Dims([2,0,-2,-1,0,0])),
            ("thermal_conductivity", Dims([1,1,-3,-1,0,0])),
        ] {
            let quantity = if name == "thermal_conductivity" {
                QuantitySpec::semantic(SemanticType::new(QuantityKind::ThermalConductivity, ValueForm::Static))
            } else { QuantitySpec::dimensional(dims) };
            let key = PropertyKey::with_quantity(name, quantity);
            let at = a.point(300.0).unwrap();
            let qa = a.card.claims().query_typed(&key, &at, fs_matdb::SelectionPolicy::SingleClaimOnly).unwrap();
            let qb = b.card.claims().query_typed(&key, &at, fs_matdb::SelectionPolicy::SingleClaimOnly).unwrap();
            assert_eq!(qa.evidence.value.value, qb.evidence.value.value);
        }
        let mu = PropertyKey::with_quantity("dynamic_viscosity", QuantitySpec::dimensional(Dims([-1,1,-1,0,0,0])));
        assert!(b.card.claims().query_typed(&mu, &b.point(300.0).unwrap(), fs_matdb::SelectionPolicy::SingleClaimOnly).is_err());
        assert!(ThermalTable::parse_solid(liquid).is_err());
        assert!(ThermalTable::parse(&solid).is_err());
        assert!(ThermalTable::parse_solid(&solid.replace("pressure_pa=100000", "")).is_err());
    }
}
