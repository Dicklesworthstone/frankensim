//! Cold table input. All interpolation, caloric integration and source selection
//! remain with fs-matdb/fs-material; this module does not implement material laws.
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_material::conductor::{ELECTRICAL_RESISTIVITY_PROPERTY, ELECTRICAL_RESISTIVITY_DIMS};
use fs_material::state_point::*;
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_qty::{Density, Pressure, QuantitySpec, Temperature};
use super::{Error, invalid};

pub const HEADER: &str = "frankensim-electrothermal-material-v1";
pub const MAX_BYTES: usize = 65536;

pub struct Source {
    pub card: MaterialCard,
    pub caloric: ResolvedHeatCapacityCurve,
    pub reference_k: f64,
}
pub fn point(t: f64) -> Result<QueryPoint, Error> {
    Ok(QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), t)?
        .with_quantity("pressure", QuantitySpec::dimensional(Pressure::DIMS), 100_000.)?)
}

impl Source {
    /// SI columns: T, reference density, Cp, alpha, electrical rho, thermal k.
    /// The density at the first row initializes mass, which then stays fixed.
    pub fn parse(text: &str, label: &str) -> Result<Self, Error> {
        if text.len() > MAX_BYTES { return Err(invalid("material input exceeds 64 KiB")); }
        let mut lines = text.lines().enumerate().filter_map(|(i, line)| {
            let line = line.split('#').next().unwrap_or("").trim();
            (!line.is_empty()).then_some((i + 1, line))
        });
        if lines.next().map(|(_,s)| s) != Some(HEADER) {
            return Err(invalid("missing electrothermal material v1 header"));
        }
        let mut rows: Vec<[f64; 6]> = Vec::new();
        for (line, text) in lines {
            if rows.len() == 256 { return Err(invalid("material table exceeds 256 rows")); }
            let columns = text.split_whitespace().collect::<Vec<_>>();
            if columns.len() != 6 { return Err(invalid(&format!("line {line}: expected six SI columns"))); }
            let mut row = [0.0; 6];
            for (out, value) in row.iter_mut().zip(columns) {
                *out = value.parse().map_err(|_| invalid(&format!("line {line}: invalid number")))?;
            }
            if row.iter().any(|x| !x.is_finite())
                || [0,1,2,4,5].into_iter().any(|i| row[i] <= 0.)
                || rows.last().is_some_and(|previous| row[0] <= previous[0]) {
                return Err(invalid(&format!("line {line}: finite data, positive T/rho/Cp/resistivity/k and increasing T required")));
            }
            rows.push(row);
        }
        if rows.len() < 2 { return Err(invalid("at least two material rows required")); }
        let lower = rows[0][0];
        let upper = rows[rows.len()-1][0];
        let properties = [
            (DENSITY_PROPERTY, Density::DIMS),
            (SPECIFIC_HEAT_CAPACITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS),
            (LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, INVERSE_TEMPERATURE_DIMS),
            (ELECTRICAL_RESISTIVITY_PROPERTY, ELECTRICAL_RESISTIVITY_DIMS),
            (THERMAL_CONDUCTIVITY_PROPERTY, THERMAL_CONDUCTIVITY_DIMS),
        ];
        let mut claims = ClaimSet::new();
        for (column, (name, dims)) in properties.into_iter().enumerate() {
            claims.insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve { abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                    knots: rows.iter().map(|r| (r[0], r[column+1])).collect(), dims },
                validity: ValidityDomain::unconstrained()
                    .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), lower, upper)
                    .with_quantity("pressure", QuantitySpec::dimensional(Pressure::DIMS), 100_000., 100_000.),
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated, observations: vec![],
                provenance: Provenance { source: format!("user table {label}; not independently verified"),
                    license: "caller-supplied".into(), artifact: None },
            })?;
        }
        let card = MaterialCard::assemble(MaterialStateId { chemistry: label.into(), phase: "solid".into(),
            process: "user-supplied electrothermal table".into(), revision: 0 }, claims, vec![])?;
        let caloric = ResolvedHeatCapacityCurve::try_new(&card, HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new(SPECIFIC_HEAT_CAPACITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS),
            density: &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
            temperature_axis: "T", lower: &point(lower)?, upper: &point(upper)?,
            selection: fs_matdb::SelectionPolicy::SingleClaimOnly,
        }, SolidLiquidPhase::Solid, 0., 1e-8)?;
        Ok(Self { card, caloric, reference_k: lower })
    }

    pub fn transport(&self, t: f64) -> Result<fs_lbm::conjugate::enthalpy::transport::EnthalpyTransport<ResolvedMaterialStatePoint>,
        fs_lbm::conjugate::enthalpy::JouleEnthalpyError> {
        use fs_lbm::conjugate::{ChtError, enthalpy::transport::EnthalpyTransport};
        let source = || -> Result<_, Error> {
            let required = [
                ScalarPropertyRequirement::try_new(ELECTRICAL_RESISTIVITY_PROPERTY,
                    ELECTRICAL_RESISTIVITY_DIMS, ScalarAdmissibility::StrictlyPositive)?,
                ScalarPropertyRequirement::try_new(THERMAL_CONDUCTIVITY_PROPERTY,
                    THERMAL_CONDUCTIVITY_DIMS, ScalarAdmissibility::StrictlyPositive)?,
            ];
            Ok(resolve_material_state_point(&self.card, &point(t)?, &required,
                MaterialPropertySelection::SingleClaimOnly)?)
        };
        let material = source().map_err(|e| ChtError::InvalidInput {
            field: "electrothermal table", reason: e.to_string(),
        })?;
        let resistivity = material.property(ELECTRICAL_RESISTIVITY_PROPERTY)
            .ok_or(ChtError::InvalidInput { field: "resistivity", reason: "required property absent".into() })?.value_si();
        let k = material.property(THERMAL_CONDUCTIVITY_PROPERTY)
            .ok_or(ChtError::InvalidInput { field: "conductivity", reason: "required property absent".into() })?.value_si();
        EnthalpyTransport::try_new(material, Some(resistivity), k)
    }
}
