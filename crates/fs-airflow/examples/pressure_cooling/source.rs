//! Small example-only ingress: authored SI rows become ordinary material cards.
use std::collections::BTreeMap;

use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, QuantitySpec, Temperature};

use super::Result;

pub struct LiquidTable {
    pub card: MaterialCard,
    pub lower_k: f64,
    pub upper_k: f64,
    pub pressure_pa: f64,
}

impl LiquidTable {
    pub fn parse(text: &str) -> Result<Self> {
        if text.len() > 128 * 1024 {
            return Err("material table exceeds 128 KiB".into());
        }
        let mut headers = BTreeMap::new();
        let mut rows: Vec<[f64; 5]> = Vec::new();
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
                let row: [f64; 5] = values.try_into()
                    .map_err(|_| format!("line {} needs T rho Cp k mu", index + 1))?;
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
            chemistry: name.into(), phase: "liquid".into(),
            process: "authored pressure-cooling example table".into(), revision: 0,
        }, claims, vec![])?;
        Ok(Self { card, lower_k, upper_k, pressure_pa })
    }

    pub fn point(&self, temperature_k: f64) -> Result<QueryPoint> {
        Ok(QueryPoint::new().with("T", temperature_k)?.with("pressure", self.pressure_pa)?)
    }
}
