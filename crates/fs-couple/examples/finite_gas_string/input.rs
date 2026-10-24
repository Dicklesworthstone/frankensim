//! Explicit bounded native input; all physical values are SI.

use fs_couple::cavity_phs::piston::thermal::HeatExchangeSolve;
use fs_couple::acoustic_realize::thermoelastic::ThermoelasticStringSolve;
use fs_material::conductor::OhmicDrive;
use fs_material::thermochemical_gas::energy::GasEnergySolve;

pub const MAX_INPUT_BYTES: usize = 8192;
pub const MAX_SAMPLES: usize = 262_144;

#[derive(Clone, Copy, Debug)]
pub enum Prestress {
    FixedTension(f64),
    FixedSupports,
}

#[derive(Clone, Copy, Debug)]
pub struct Thermoelastic {
    pub reference_length_m: f64,
    pub reference_temperature_k: f64,
    pub strain_limit: f64,
    pub solve: ThermoelasticStringSolve,
}

#[derive(Clone, Debug)]
pub struct Spec {
    pub sample_rate_hz: u32,
    pub samples: usize,
    pub full_scale_pa: f64,
    pub length_m: f64,
    pub radius_m: f64,
    pub prestress: Prestress,
    pub modes: usize,
    pub station: f64,
    pub force_n: f64,
    pub release_sample: usize,
    pub listener_m: f64,
    pub solid_temperature_k: f64,
    pub solid_minimum_k: f64,
    pub solid_maximum_k: f64,
    pub query_omega_rad_s: f64,
    pub enthalpy_interpolation_j_kg: f64,
    pub maximum_solid_increment_k: f64,
    pub gas_volume_m3: f64,
    pub gas_solve: GasEnergySolve,
    pub maximum_gas_increment_k: f64,
    pub transfer_area_m2: f64,
    pub gas_distance_m: f64,
    pub solid_distance_m: f64,
    pub exchange: HeatExchangeSolve,
    pub drive: Option<OhmicDrive>,
    pub solid_power_w: f64,
    pub gas_power_w: f64,
    pub thermoelastic: Option<Thermoelastic>,
}

struct Fields<'a>(std::str::SplitAsciiWhitespace<'a>);
impl Fields<'_> {
    fn keyword(&mut self, expected: &str) -> super::Result<()> {
        if self.0.next() != Some(expected) {
            return Err(format!("expected {expected}").into());
        }
        Ok(())
    }
    fn value<T>(&mut self, label: &str) -> super::Result<T>
    where
        T: std::str::FromStr,
        T::Err: std::error::Error + 'static,
    {
        Ok(self
            .0
            .next()
            .ok_or_else(|| format!("missing {label}"))?
            .parse()?)
    }
}

pub fn parse(text: &str) -> super::Result<Spec> {
    if text.len() > MAX_INPUT_BYTES {
        return Err("string scenario exceeds 8192 bytes".into());
    }
    let mut f = Fields(text.split_ascii_whitespace());
    f.keyword("frankensim-finite-gas-string-v1")?;
    f.keyword("audio")?;
    let (sample_rate_hz, samples, full_scale_pa) = (
        f.value("sample rate")?,
        f.value("samples")?,
        f.value("full scale")?,
    );
    f.keyword("string")?;
    let length_m = f.value("length")?;
    let radius_m = f.value("radius")?;
    let prestress = match f.0.next().ok_or("missing string tension or fixed_support")? {
        "fixed_support" => Prestress::FixedSupports,
        value => Prestress::FixedTension(value.parse()?),
    };
    let modes = f.value("modes")?;
    f.keyword("release")?;
    let (station, force_n, release_sample) = (
        f.value("station")?,
        f.value("point force")?,
        f.value("release sample")?,
    );
    f.keyword("listener")?;
    let listener_m = f.value("listener distance")?;
    f.keyword("solid")?;
    let (
        solid_temperature_k,
        solid_minimum_k,
        solid_maximum_k,
        query_omega_rad_s,
        enthalpy_interpolation_j_kg,
        maximum_solid_increment_k,
    ) = (
        f.value("solid T")?,
        f.value("solid minimum T")?,
        f.value("solid maximum T")?,
        f.value("query omega")?,
        f.value("enthalpy interpolation allowance")?,
        f.value("solid temperature step limit")?,
    );
    f.keyword("gas")?;
    let gas_volume_m3 = f.value("gas volume")?;
    let gas_solve = GasEnergySolve {
        minimum_temperature_k: f.value("gas minimum T")?,
        maximum_temperature_k: f.value("gas maximum T")?,
        energy_tolerance_j_kg: f.value("gas energy tolerance")?,
        maximum_iterations: f.value("gas iterations")?,
    };
    let maximum_gas_increment_k = f.value("gas temperature step limit")?;
    f.keyword("conduction")?;
    let (transfer_area_m2, gas_distance_m, solid_distance_m) = (
        f.value("transfer area")?,
        f.value("gas layer distance")?,
        f.value("solid layer distance")?,
    );
    f.keyword("exchange")?;
    let maximum_heat: f64 = f.value("exchange heat bracket")?;
    let exchange = HeatExchangeSolve {
        heat_bracket_j: [-maximum_heat, maximum_heat],
        source_tolerance_j: f.value("exchange source tolerance")?,
        energy_tolerance_j: f.value("exchange energy tolerance")?,
        maximum_iterations: f.value("exchange iterations")?,
    };
    f.keyword("drive")?;
    let drive = match f.0.next() {
        Some("none") => None,
        Some("current") => Some(OhmicDrive::Current(f.value("current")?)),
        Some("voltage") => Some(OhmicDrive::Voltage(f.value("voltage")?)),
        _ => return Err("drive requires none, current A or voltage V".into()),
    };
    f.keyword("heat")?;
    let (solid_power_w, gas_power_w) = (f.value("solid power")?, f.value("gas power")?);
    let thermoelastic = match f.0.next() {
        None => None,
        Some("thermoelastic") => {
            let reference_length_m = f.value("reference stress-free length")?;
            let reference_temperature_k = f.value("reference temperature")?;
            let strain_limit = f.value("linear strain limit")?;
            let bound: f64 = f.value("reversible heat bracket")?;
            Some(Thermoelastic {
                reference_length_m,
                reference_temperature_k,
                strain_limit,
                solve: ThermoelasticStringSolve {
                    reversible_heat_bracket_j: [-bound, bound],
                    energy_tolerance_j: f.value("reversible heat tolerance")?,
                    maximum_iterations: f.value("reversible heat iterations")?,
                },
            })
        }
        _ => return Err("only a thermoelastic record may follow heat".into()),
    };
    if f.0.next().is_some() {
        return Err("string scenario contains trailing fields".into());
    }
    let spec = Spec {
        sample_rate_hz,
        samples,
        full_scale_pa,
        length_m,
        radius_m,
        prestress,
        modes,
        station,
        force_n,
        release_sample,
        listener_m,
        solid_temperature_k,
        solid_minimum_k,
        solid_maximum_k,
        query_omega_rad_s,
        enthalpy_interpolation_j_kg,
        maximum_solid_increment_k,
        gas_volume_m3,
        gas_solve,
        maximum_gas_increment_k,
        transfer_area_m2,
        gas_distance_m,
        solid_distance_m,
        exchange,
        drive,
        solid_power_w,
        gas_power_w,
        thermoelastic,
    };
    spec.validate()?;
    Ok(spec)
}

impl Spec {
    pub fn validate(&self) -> super::Result<()> {
        self.gas_solve.validate()?;
        if !(8000..=192_000).contains(&self.sample_rate_hz)
            || !(1..=MAX_SAMPLES).contains(&self.samples)
            || !(1..=32).contains(&self.modes)
            || self.release_sample > self.samples
        {
            return Err("audio requires rate 8000..192000, 1..262144 samples, 1..32 modes and a release inside the run".into());
        }
        if [
            self.full_scale_pa,
            self.length_m,
            self.radius_m,
            self.listener_m,
            self.solid_temperature_k,
            self.solid_minimum_k,
            self.solid_maximum_k,
            self.query_omega_rad_s,
            self.enthalpy_interpolation_j_kg,
            self.maximum_solid_increment_k,
            self.gas_volume_m3,
            self.maximum_gas_increment_k,
            self.gas_distance_m,
            self.solid_distance_m,
            self.exchange.heat_bracket_j[1],
            self.exchange.source_tolerance_j,
            self.exchange.energy_tolerance_j,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x <= 0.0)
        {
            return Err(
                "physical magnitudes, brackets and tolerances must be finite and positive".into(),
            );
        }
        match (self.prestress, self.thermoelastic) {
            (Prestress::FixedTension(tension), None) if tension.is_finite() && tension > 0.0 => {}
            (Prestress::FixedSupports, Some(thermal)) => {
                if [thermal.reference_length_m, thermal.reference_temperature_k,
                    thermal.strain_limit, thermal.solve.reversible_heat_bracket_j[1],
                    thermal.solve.energy_tolerance_j]
                    .iter().any(|value| !value.is_finite() || *value <= 0.0)
                    || !(self.solid_minimum_k..=self.solid_maximum_k)
                        .contains(&thermal.reference_temperature_k)
                    || thermal.solve.reversible_heat_bracket_j[0]
                        != -thermal.solve.reversible_heat_bracket_j[1]
                    || !(1..=128).contains(&thermal.solve.maximum_iterations)
                {
                    return Err("invalid thermoelastic reference, strain bound or heat controls".into());
                }
            }
            _ => return Err("use positive fixed tension without thermoelastic, or fixed_support with its thermoelastic record".into()),
        }
        if !self.station.is_finite()
            || !(0.0..=1.0).contains(&self.station)
            || !self.force_n.is_finite()
            || !self.solid_power_w.is_finite()
            || !self.gas_power_w.is_finite()
            || self.drive.is_some_and(|drive| match drive {
                OhmicDrive::Current(i) | OhmicDrive::Voltage(i) => !i.is_finite(),
            })
        {
            return Err("station must be in [0,1]; force, drive and powers must be finite".into());
        }
        if self.solid_minimum_k >= self.solid_maximum_k
            || !(self.solid_minimum_k..=self.solid_maximum_k).contains(&self.solid_temperature_k)
        {
            return Err(
                "initial solid temperature must lie in a nonempty declared solid interval".into(),
            );
        }
        let surface = 2.0 * std::f64::consts::PI * self.radius_m * (self.length_m + self.radius_m);
        if !surface.is_finite()
            || !self.transfer_area_m2.is_finite()
            || self.transfer_area_m2 < 0.0
            || self.transfer_area_m2 > surface
            || self.exchange.heat_bracket_j[0] != -self.exchange.heat_bracket_j[1]
            || !(1..=128).contains(&self.exchange.maximum_iterations)
        {
            return Err(
                "transfer area exceeds the cylinder surface or exchange controls are invalid"
                    .into(),
            );
        }
        Ok(())
    }
}
