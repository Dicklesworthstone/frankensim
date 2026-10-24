//! Example input is explicitly synthetic; the shared material owners evaluate it.
use std::error::Error;

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::liquid::{DYNAMIC_VISCOSITY_DIMS, DYNAMIC_VISCOSITY_PROPERTY};
use fs_material::state_point::{
    DENSITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS, SPECIFIC_HEAT_CAPACITY_PROPERTY,
    THERMAL_CONDUCTIVITY_DIMS, THERMAL_CONDUCTIVITY_PROPERTY,
};
use fs_qty::{Density, QuantitySpec, Temperature};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

pub struct MaterialInput {
    pub card: MaterialCard,
    pub cp: f64,
    pub reference_density: f64,
    pub temperatures: [f64; 2],
}

pub fn absolute_temperature() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature, ValueForm::Static,
    ))
}

pub fn point(temperature: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", absolute_temperature(), temperature).expect("finite admitted temperature")
        .with("pressure", 100_000.0).expect("finite scenario pressure")
}

impl MaterialInput {
    /// Cp followed by 2..64 rows of T, density, dynamic viscosity, conductivity.
    /// No interpolation, extrapolation or gas/liquid selection happens here.
    pub fn parse(text: &str) -> Result<Self, Box<dyn Error>> {
        if text.len() > 65_536 { return Err("material input exceeds 64 KiB".into()); }
        let rows: Vec<Vec<f64>> = text.lines().filter_map(|line| {
            let line = line.split('#').next().unwrap_or("").trim();
            (!line.is_empty()).then_some(line)
        }).map(|line| line.split_whitespace().map(str::parse).collect())
            .collect::<Result<_, _>>()?;
        if !(3..=65).contains(&rows.len()) || rows[0].len() != 1 {
            return Err("expected Cp then 2..64 rows of T rho mu k, all in SI units".into());
        }
        let cp = rows[0][0];
        if !cp.is_finite() || cp <= 0.0 { return Err("Cp must be finite and positive".into()); }
        for (i, row) in rows[1..].iter().enumerate() {
            if row.len() != 4 || row.iter().any(|x| !x.is_finite() || *x <= 0.0)
                || (i > 0 && row[0] <= rows[i][0]) {
                return Err("material rows need positive finite T/rho/mu/k and increasing T".into());
            }
        }
        let temperatures = [rows[1][0], rows.last().expect("admitted rows")[0]];
        let validity = ValidityDomain::unconstrained()
            .with_quantity("T", absolute_temperature(), temperatures[0], temperatures[1])
            .with("pressure", 100_000.0, 100_000.0);
        let mut claims = ClaimSet::new();
        for (name, dims, column) in [
            (DENSITY_PROPERTY, Density::DIMS, 1),
            (DYNAMIC_VISCOSITY_PROPERTY, DYNAMIC_VISCOSITY_DIMS, 2),
            (THERMAL_CONDUCTIVITY_PROPERTY, THERMAL_CONDUCTIVITY_DIMS, 3),
            (SPECIFIC_HEAT_CAPACITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS, 0),
        ] {
            let quantity = if name == THERMAL_CONDUCTIVITY_PROPERTY {
                QuantitySpec::semantic(SemanticType::new(
                    QuantityKind::ThermalConductivity, ValueForm::Static,
                ))
            } else { QuantitySpec::dimensional(dims) };
            let value = if column == 0 { PropertyValue::Scalar { value: cp, dims } }
                else { PropertyValue::Curve {
                    abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                    knots: rows[1..].iter().map(|row| (row[0], row[column])).collect(), dims,
                } };
            claims.insert_claim(PropertyClaim {
                key: PropertyKey::with_quantity(name, quantity), value,
                validity: validity.clone(), uncertainty: UncertaintyModel::Unstated,
                interpolation: if column == 0 { InterpolationPolicy::ConstantWithinValidity }
                    else { InterpolationPolicy::LinearInside },
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic material_buoyant_flow input; not measured substance data".into(),
                    license: "example input".into(), artifact: None,
                },
            })?;
        }
        let card = MaterialCard::assemble(MaterialStateId {
            chemistry: "synthetic Newtonian liquid".into(), phase: "liquid".into(),
            process: "unqualified example table".into(), revision: 0,
        }, claims, vec![])?;
        Ok(Self { card, cp, reference_density: rows[1][1], temperatures })
    }
}
