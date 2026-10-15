//! Explicit example heat drives, evaluated on the shared midpoint material.
//! The conductive option is a lumped gas layer k*A/d coupled to an infinite
//! reservoir, not resolved convection, a Nusselt correlation or a wall PDE.
use super::Result;
use fs_material::thermochemical_gas::ResolvedThermochemicalGasState;

#[derive(Clone, Copy, Debug)]
pub enum HeatDrive {
    PrescribedPower(f64),
    /// Existing driver API supplies whole-step energy without W/J roundtripping.
    PrescribedEnergy(f64),
    ConductiveReservoir { temperature_k: f64, area_m2: f64, distance_m: f64 },
}
impl HeatDrive {
    pub fn parse(args: &[String]) -> Result<Self> {
        let drive = match args {
            [] => Self::PrescribedPower(0.0),
            [power] => Self::PrescribedPower(power.parse()?),
            [mode, temperature, area, distance] if mode == "wall" => Self::ConductiveReservoir {
                temperature_k: temperature.parse()?, area_m2: area.parse()?, distance_m: distance.parse()?,
            },
            _ => return Err("heat drive: [heat_w] or wall temperature_K area_m2 distance_m required".into()),
        };
        drive.validate()?;
        Ok(drive)
    }

    pub fn validate(self) -> Result<()> {
        let valid = match self {
            Self::PrescribedPower(watts) => watts.is_finite(),
            Self::PrescribedEnergy(joules) => joules.is_finite(),
            Self::ConductiveReservoir { temperature_k, area_m2, distance_m } => {
                temperature_k.is_finite() && temperature_k > 0.0
                    && area_m2.is_finite() && area_m2 >= 0.0
                    && distance_m.is_finite() && distance_m > 0.0
            }
        };
        if valid { Ok(()) } else { Err("invalid heat power or reservoir/film geometry".into()) }
    }

    /// Whole-step heat from the ACTUAL midpoint, with no frozen conductivity.
    /// Area is the heat-transfer area, distinct from the mechanical piston area.
    pub fn energy_j(self, dt: f64, midpoint: &ResolvedThermochemicalGasState) -> Result<f64> {
        self.validate()?;
        if !dt.is_finite() || dt <= 0.0 { return Err("positive finite heat interval required".into()); }
        let power = match self {
            Self::PrescribedPower(watts) => watts,
            Self::PrescribedEnergy(joules) => return Ok(joules),
            Self::ConductiveReservoir { temperature_k, area_m2, distance_m } => {
                let gas = midpoint.state();
                let conductance = gas.thermal_conductivity*(area_m2/distance_m);
                let delta = temperature_k-gas.temperature;
                let watts = conductance*delta;
                if !conductance.is_finite() || conductance < 0.0
                    || (area_m2 != 0.0 && conductance == 0.0)
                    || !watts.is_finite() || (conductance != 0.0 && delta != 0.0 && watts == 0.0)
                {
                    return Err("material gas-layer heat flow is unrepresentable".into());
                }
                watts
            }
        };
        let heat = power*dt;
        if !heat.is_finite() || (power != 0.0 && heat == 0.0) {
            return Err("integrated heat is not representable in joules".into());
        }
        Ok(heat)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::source;
    const A: &str = include_str!("mixture-a.tsv");

    #[test]
    fn only_material_conductivity_changes_the_conductive_heat_rate() {
        let model = source::load(A, "synthetic original").unwrap();
        let changed = A.replace("0.025", "0.050").replace("0.150", "0.300");
        let changed = source::load(&changed, "synthetic doubled conductivity").unwrap();
        let point = source::point(400.0, 100000.0).unwrap();
        let original = model.evaluate(&point).unwrap();
        let doubled = changed.evaluate(&point).unwrap();
        let wall = HeatDrive::ConductiveReservoir { temperature_k: 500.0, area_m2: 0.02, distance_m: 0.001 };
        let q = wall.energy_j(0.0005, &original).unwrap();
        assert!(q > 0.0);
        assert!((wall.energy_j(0.0005, &doubled).unwrap()/q-2.0).abs() < 1e-12);
        assert_eq!(original.state().specific_heat_cp, doubled.state().specific_heat_cp);
        assert_eq!(original.state().sound_speed, doubled.state().sound_speed);
        let cold = HeatDrive::ConductiveReservoir { temperature_k: 300.0, area_m2: 0.02, distance_m: 0.001 };
        assert!((cold.energy_j(0.0005, &original).unwrap()+q).abs() < 1e-12);
        let isolated = HeatDrive::ConductiveReservoir { temperature_k: 500.0, area_m2: 0.0, distance_m: 0.001 };
        assert_eq!(isolated.energy_j(0.0005, &original).unwrap(), 0.0);
    }

    #[test]
    fn heat_drive_inputs_are_explicit_and_invalid_data_refuses() {
        let parse = |args: &[&str]| HeatDrive::parse(&args.iter().map(|s| (*s).to_owned()).collect::<Vec<_>>());
        assert!(parse(&[]).is_ok());
        assert!(parse(&["-100"]).is_ok());
        assert!(parse(&["wall", "500", "0.02", "0.001"]).is_ok());
        for args in [vec!["NaN"], vec!["wall"], vec!["unknown", "500", "0.02", "0.001"],
            vec!["wall", "0", "0.02", "0.001"], vec!["wall", "500", "-1", "0.001"],
            vec!["wall", "500", "0.02", "0"], vec!["wall", "inf", "0.02", "0.001"]] {
            assert!(parse(&args).is_err());
        }
    }
}
