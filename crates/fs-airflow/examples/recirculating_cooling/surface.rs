//! Example-only surface/environment ingress into the existing radiation owner.
//! A scalar emissivity is deliberately constant within its stated T domain.
//! Convection coefficient, air temperature and radiative reservoir are separate
//! environmental declarations, not properties silently inferred from a coating.
use std::collections::BTreeMap;

use fs_conduction::{AmbientRadiationPatch, SurfaceEmissivity, ThermalBc, SURFACE_EMISSIVITY_PROPERTY};
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, SelectionPolicy, UncertaintyModel};
use fs_qty::{Dims, QuantitySpec};
use super::{INITIAL_K, Result};

pub const EXPOSED_REGION: &str = "exposed-bottom";

pub struct SurfaceTable {
    pub patch: AmbientRadiationPatch,
    pub convection_w_m2_k: f64,
    pub air_temperature_k: f64,
    pub license: String,
}
impl SurfaceTable {
    pub fn parse(text: &str) -> Result<Self> {
        if text.len() > 128*1024 { return Err("surface description exceeds 128 KiB".into()); }
        let mut values = BTreeMap::new();
        for (index, line) in text.lines().enumerate() {
            let line = line.split('#').next().unwrap_or("").trim();
            if line.is_empty() { continue; }
            let (key, value) = line.split_once('=')
                .ok_or_else(|| format!("surface line {} needs key=value", index+1))?;
            let (key, value) = (key.trim(), value.trim());
            if !matches!(key, "material" | "source" | "license" | "emissivity"
                | "minimum_surface_k" | "maximum_surface_k" | "air_temperature_k"
                | "convection_w_m2_k" | "radiation_temperature_k")
                || value.is_empty() || values.insert(key, value).is_some()
            { return Err(format!("unknown, empty or duplicate surface field on line {}", index+1).into()); }
        }
        let get = |key| values.get(key).copied().ok_or_else(|| format!("missing surface {key}= field"));
        let number = |key| -> Result<f64> { Ok(get(key)?.parse()?) };
        let (low, high) = (number("minimum_surface_k")?, number("maximum_surface_k")?);
        let epsilon = number("emissivity")?;
        let air_temperature_k = number("air_temperature_k")?;
        let convection_w_m2_k = number("convection_w_m2_k")?;
        let ambient = number("radiation_temperature_k")?;
        if [low, high, air_temperature_k, ambient].iter().any(|v| !v.is_finite() || *v<=0.0)
            || low>high || INITIAL_K<low || INITIAL_K>high
        { return Err("finite positive environment temperatures and a surface span containing 300 K required".into()); }
        // No tiny artificial h or inferred air temperature. Pure-radiation
        // boundaries are outside the existing ambient-patch solver's model.
        ThermalBc::robin(convection_w_m2_k, air_temperature_k)?;
        let mut claims = ClaimSet::new();
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::with_quantity(SURFACE_EMISSIVITY_PROPERTY, QuantitySpec::dimensional(Dims::NONE)),
            value: PropertyValue::Scalar { value: epsilon, dims: Dims::NONE },
            validity: ValidityDomain::unconstrained().with("T", low, high),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance { source: get("source")?.into(), license: get("license")?.into(), artifact: None },
        })?;
        let card = MaterialCard::assemble(MaterialStateId {
            chemistry: get("material")?.into(), phase: "solid".into(),
            process: "authored constant surface finish".into(), revision: 0,
        }, claims, vec![])?;
        let emissivity = SurfaceEmissivity::from_card(EXPOSED_REGION, &card,
            INITIAL_K, SelectionPolicy::SingleClaimOnly)?;
        Ok(Self { patch: AmbientRadiationPatch::new(EXPOSED_REGION, emissivity, ambient)?,
            convection_w_m2_k, air_temperature_k, license: get("license")?.to_owned() })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    const LOW: &str = include_str!("surface_low_emissivity.txt");
    const HIGH: &str = include_str!("surface_high_emissivity.txt");

    #[test]
    fn surface_data_retains_source_and_separate_environment_without_defaults() {
        let a = SurfaceTable::parse(LOW).unwrap();
        let b = SurfaceTable::parse(HIGH).unwrap();
        assert_eq!(a.patch.emissivity().value(), 0.1);
        assert_eq!(b.patch.emissivity().value(), 0.9);
        assert_eq!(a.air_temperature_k, 300.0);
        assert_eq!(a.patch.ambient_temperature_k(), 250.0);
        assert_eq!(a.license, "CC0-1.0");
        assert_ne!(a.patch.emissivity().card_identity(), b.patch.emissivity().card_identity());
        for bad in [LOW.replace("emissivity=0.1", "emissivity=NaN"),
            LOW.replace("emissivity=0.1", "emissivity=1.1"),
            LOW.replace("convection_w_m2_k=5", "convection_w_m2_k=0"),
            LOW.replace("air_temperature_k=300", ""),
            LOW.replace("minimum_surface_k=250", "minimum_surface_k=350"),
            LOW.replace("license=CC0-1.0", ""), format!("{LOW}\nemissivity=0.2")]
        { assert!(SurfaceTable::parse(&bad).is_err()); }
    }

    #[test]
    fn radiation_option_requires_storage_and_keeps_existing_numeric_arguments() {
        let base: Vec<String> = ["fluid.txt", "pump.txt", "0.1", "1", "2", "20", "0", "299", "302", "20", "80"]
            .into_iter().map(str::to_owned).collect();
        let (_, _, _, legacy) = super::super::parse_args(&base).unwrap();
        let mut args = vec!["--radiation".into(), "surface.txt".into()];
        args.extend(base.clone());
        assert!(super::super::parse_radiation_args(&args).is_err());
        args.splice(2..2, ["--solid".into(), "solid.txt".into()]);
        let (surface, solid, fluid, pump, parsed) = super::super::parse_radiation_args(&args).unwrap();
        assert_eq!((surface, solid, fluid, pump), (Some("surface.txt"), Some("solid.txt"), "fluid.txt", "pump.txt"));
        assert_eq!(parsed, legacy);
        assert!(super::super::parse_radiation_args(&["--radiation".into()]).is_err());
        assert!(super::super::parse_radiation_args(&base).unwrap().0.is_none());
    }

    #[test]
    fn surface_only_changes_reach_actual_fem_and_separated_joint_energy() {
        use super::super::{LiquidTable, PumpTable, Settings, ThermalTable, run_with_surface};
        let fluid = LiquidTable::parse(include_str!("../pressure_cooling/constant.txt")).unwrap();
        let pump = PumpTable::parse(include_str!("../pressure_cooling/pump.txt")).unwrap();
        let solid = ThermalTable::parse_solid(include_str!("solid_low_capacity.txt")).unwrap();
        let low = SurfaceTable::parse(LOW).unwrap();
        let high = SurfaceTable::parse(HIGH).unwrap();
        let settings = Settings { mass_kg: 0.1, dt_s: 1.0, steps: 2, body_w: 20.0,
            reservoir_w: 0.0, inlet_bracket_k: (299.0,302.0), pressure_bracket_pa: (20.0,80.0) };
        let a = run_with_surface(&fluid,&pump,Some(&solid),Some(&low),settings).unwrap();
        let b = run_with_surface(&fluid,&pump,Some(&solid),Some(&high),settings).unwrap();
        assert!(b[1].reservoir_k < a[1].reservoir_k);
        assert!(b[1].walls_k.iter().zip(&a[1].walls_k).all(|(b,a)|b<a));
        for rows in [&a, &b] {
            for row in rows {
                let storage = row.storage.as_ref().unwrap();
                let radiation = row.radiation.as_ref().unwrap();
                assert!(radiation.outward_j>0.0);
                assert!((storage.external_j + radiation.outward_j + radiation.external_convection_j - 20.0).abs()<1e-7);
                assert!(storage.combined_residual_j.abs()<=5e-4);
                assert!(storage.interface_residual_j.abs()<=5e-4);
                assert!((row.solid_heat_w-row.return_j).abs()<=5e-4);
                assert!(radiation.enthalpy_solves>=row.solid_calls);
            }
        }
        let hot = SurfaceTable::parse(&HIGH.replace("radiation_temperature_k=250", "radiation_temperature_k=400")).unwrap();
        let rows = run_with_surface(&fluid,&pump,Some(&solid),Some(&hot),
            Settings { body_w: 0.0, steps: 1, ..settings }).unwrap();
        assert!(rows[0].radiation.as_ref().unwrap().outward_j<0.0);
        assert!(rows[0].reservoir_k>INITIAL_K);
        assert!(run_with_surface(&fluid,&pump,None,Some(&high),settings).is_err());
    }
}
