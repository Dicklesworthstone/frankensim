//! Explicit bounded native input; all physical values are SI.

use fs_couple::cavity_phs::piston::thermal::HeatExchangeSolve;
use fs_couple::acoustic_realize::thermoelastic::ThermoelasticStringSolve;
use fs_material::conductor::OhmicDrive;
use fs_material::thermochemical_gas::energy::GasEnergySolve;
use fs_couple::acoustic_realize::concentric::ConcentricStringThermalSolve;
use fs_solver::NewtonKrylovConfig;
use fs_material::thermochemical_gas::energy::parcel::transfer::{GasMassTransferSolve,
    nozzle::{GasNozzle, GasNozzleStepSolve}};

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
pub struct Sleeve {
    pub pack: String,
    pub outer_radius_m: f64,
    pub interface_area_m2: f64,
    pub resistance_geometry_per_m: [f64; 2],
}

#[derive(Clone, Debug)]
pub struct Bonded {
    pub core_radius_m: f64,
    pub reference_length_m: f64,
    pub reference_temperature_k: f64,
    pub strain_limit: f64,
    pub sleeves: Vec<Sleeve>,
    pub solve: ConcentricStringThermalSolve,
}

/// Authored restriction and a second finite inventory, never an infinite bath.
#[derive(Clone, Debug)]
pub struct Nozzle {
    pub reservoir_record: String,
    pub reservoir_volume_m3: f64,
    pub geometry: GasNozzle,
    pub mass_tolerance_kg: f64,
    pub energy_tolerance_j: f64,
    pub maximum_donor_fraction: f64,
}
impl Nozzle {
    pub fn controls(&self, gas: GasEnergySolve) -> GasNozzleStepSolve {
        GasNozzleStepSolve { nozzle: gas, transfer: GasMassTransferSolve {
            donor: gas, receiver: gas, mass_tolerance_kg: self.mass_tolerance_kg,
            energy_tolerance_j: self.energy_tolerance_j,
        }, maximum_donor_fraction: self.maximum_donor_fraction }
    }
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
    pub nozzle: Option<Nozzle>,
    pub transfer_area_m2: f64,
    pub gas_distance_m: f64,
    pub solid_distance_m: f64,
    pub exchange: HeatExchangeSolve,
    pub drive: Option<OhmicDrive>,
    pub solid_power_w: f64,
    pub gas_power_w: f64,
    pub thermoelastic: Option<Thermoelastic>,
    pub bonded: Option<Bonded>,
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
    let nozzle = if f.0.clone().next() == Some("nozzle") {
        f.keyword("nozzle")?;
        Some(Nozzle {
            reservoir_record: f.0.next().ok_or("missing reservoir gas record")?.to_owned(),
            reservoir_volume_m3: f.value("reservoir volume")?,
            geometry: GasNozzle { area_m2: f.value("nozzle area")?,
                discharge_coefficient: f.value("nozzle discharge coefficient")? },
            mass_tolerance_kg: f.value("gas transfer mass tolerance")?,
            energy_tolerance_j: f.value("gas transfer energy tolerance")?,
            maximum_donor_fraction: f.value("maximum donor fraction")?,
        })
    } else { None };
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
    let (thermoelastic, bonded) = match f.0.next() {
        None => (None, None),
        Some("thermoelastic") => {
            let reference_length_m = f.value("reference stress-free length")?;
            let reference_temperature_k = f.value("reference temperature")?;
            let strain_limit = f.value("linear strain limit")?;
            let bound: f64 = f.value("reversible heat bracket")?;
            (Some(Thermoelastic {
                reference_length_m,
                reference_temperature_k,
                strain_limit,
                solve: ThermoelasticStringSolve {
                    reversible_heat_bracket_j: [-bound, bound],
                    energy_tolerance_j: f.value("reversible heat tolerance")?,
                    maximum_iterations: f.value("reversible heat iterations")?,
                },
            }), None)
        }
        Some("bonded") => {
            let core_radius_m = f.value("core radius")?;
            let count: usize = f.value("sleeve count")?;
            if !(1..=15).contains(&count) { return Err("bonded strings require 1..15 sleeves".into()); }
            let reference_length_m = f.value("reference length")?;
            let reference_temperature_k = f.value("reference temperature")?;
            let strain_limit = f.value("strain limit")?;
            let mut sleeves = Vec::with_capacity(count);
            for _ in 0..count {
                f.keyword("layer")?;
                let pack = f.0.next().ok_or("missing sleeve pack path")?.to_owned();
                sleeves.push(Sleeve { pack, outer_radius_m: f.value("layer radius")?,
                    interface_area_m2: f.value("interface area")?,
                    resistance_geometry_per_m: [f.value("inner half resistance geometry")?,
                        f.value("outer half resistance geometry")?] });
            }
            f.keyword("layered_solve")?;
            let solve = ConcentricStringThermalSolve {
                maximum_cells: count + 1,
                maximum_thermal_trials: f.value("thermal trial budget")?,
                maximum_newton_iterations: f.value("thermal iteration budget")?,
                difference_step_j: f.value("finite difference energy step")?,
                maximum_temperature_increment_k: maximum_solid_increment_k,
                cell_energy_tolerance_j: f.value("cell energy tolerance")?,
                thermal_energy_tolerance_j: f.value("thermal energy tolerance")?,
                total_energy_tolerance_j: f.value("combined energy tolerance")?,
                newton: NewtonKrylovConfig { absolute_tolerance: f.value("nonlinear residual tolerance")?,
                    relative_tolerance: 0.0, linear_restart: count + 1,
                    ..NewtonKrylovConfig::default() },
            };
            (None, Some(Bonded { core_radius_m, reference_length_m,
                reference_temperature_k, strain_limit, sleeves, solve }))
        }
        _ => return Err("only thermoelastic or bonded may follow heat".into()),
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
        nozzle,
        transfer_area_m2,
        gas_distance_m,
        solid_distance_m,
        exchange,
        drive,
        solid_power_w,
        gas_power_w,
        thermoelastic,
        bonded,
    };
    spec.validate()?;
    Ok(spec)
}

impl Spec {
    pub fn validate(&self) -> super::Result<()> {
        self.gas_solve.validate()?;
        if let Some(nozzle) = &self.nozzle {
            if nozzle.reservoir_record.is_empty() || nozzle.reservoir_record.len() > 512
                || !std::path::Path::new(&nozzle.reservoir_record).components()
                    .all(|part| matches!(part, std::path::Component::Normal(_)))
                || [nozzle.reservoir_volume_m3, nozzle.mass_tolerance_kg,
                    nozzle.energy_tolerance_j, nozzle.maximum_donor_fraction]
                    .iter().any(|v| !v.is_finite() || *v <= 0.0)
                || nozzle.maximum_donor_fraction >= 1.0
                || !nozzle.geometry.area_m2.is_finite() || nozzle.geometry.area_m2 < 0.0
                || !nozzle.geometry.discharge_coefficient.is_finite()
                || !(0.0..=1.0).contains(&nozzle.geometry.discharge_coefficient)
            { return Err("nozzle requires a relative gas record, finite volume/area, coefficient in [0,1], positive transfer tolerances and donor fraction in (0,1)".into()); }
        }
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
        match (self.prestress, self.thermoelastic, &self.bonded) {
            (Prestress::FixedTension(tension), None, None) if tension.is_finite() && tension > 0.0 => {}
            (Prestress::FixedSupports, Some(thermal), None) => {
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
            (Prestress::FixedSupports, None, Some(bonded)) => {
                bonded.solve.newton.validate()?;
                if bonded.sleeves.is_empty() || bonded.sleeves.len() > 15
                    || bonded.solve.maximum_cells != bonded.sleeves.len() + 1
                    || !(1..=4096).contains(&bonded.solve.maximum_thermal_trials)
                    || !(1..=128).contains(&bonded.solve.maximum_newton_iterations)
                    || [bonded.core_radius_m, bonded.reference_length_m, bonded.reference_temperature_k,
                        bonded.strain_limit, bonded.solve.difference_step_j,
                        bonded.solve.cell_energy_tolerance_j, bonded.solve.thermal_energy_tolerance_j,
                        bonded.solve.total_energy_tolerance_j, bonded.solve.newton.absolute_tolerance]
                        .iter().any(|v| !v.is_finite() || *v <= 0.0)
                    || bonded.strain_limit >= 1.0
                    || !(self.solid_minimum_k..=self.solid_maximum_k).contains(&bonded.reference_temperature_k)
                    || bonded.solve.maximum_temperature_increment_k != self.maximum_solid_increment_k
                { return Err("invalid bonded reference, layer count or thermal solve controls".into()); }
                let mut inner = bonded.core_radius_m;
                for sleeve in &bonded.sleeves {
                    if sleeve.pack.is_empty() || sleeve.pack.len() > 512
                        || !std::path::Path::new(&sleeve.pack).components()
                            .all(|part| matches!(part, std::path::Component::Normal(_)))
                        || !sleeve.outer_radius_m.is_finite() || sleeve.outer_radius_m <= inner
                        || [sleeve.interface_area_m2, sleeve.resistance_geometry_per_m[0],
                            sleeve.resistance_geometry_per_m[1]].iter().any(|v| !v.is_finite() || *v <= 0.0)
                        || sleeve.interface_area_m2 > 2.0 * std::f64::consts::PI * inner * self.length_m
                    { return Err("bonded layers need relative packs, increasing radii and valid interface geometry".into()); }
                    inner = sleeve.outer_radius_m;
                }
                if inner.to_bits() != self.radius_m.to_bits() {
                    return Err("outermost bonded layer must equal the declared string radius".into());
                }
                if self.transfer_area_m2 > 2.0 * std::f64::consts::PI * self.radius_m * self.length_m {
                    return Err("bonded gas contact must lie on the outer cylindrical surface".into());
                }
            }
            _ => return Err("use fixed tension without thermal support, or fixed_support with one thermoelastic or bonded record".into()),
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
