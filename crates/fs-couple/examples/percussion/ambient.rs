//! The gas a percussion run compresses and radiates into. It is DATA, never a
//! constant in code:
//! - an imported description (drum, shell or hi-hat file) declares its own
//!   `ambient` record and refuses without one;
//! - a built-in reference instrument loads the committed record
//!   `examples/dry-air-20c.gas` through the gas owner;
//! - `--gas-record FILE` overrides either. The command line sets it once per
//!   process, before any instrument is built; builders call [`effective`].
//! Every receipt names which of the three supplied the air.
use std::io::Read;
use std::sync::OnceLock;
use fs_bem::helmholtz::Medium;
use fs_couple::vibroacoustic::AcousticMedium;
use fs_material::gas::GasState;
use fs_material::gas_species::input::read_gas_record;
use super::Error;

/// The committed reference record shared by every built-in instrument.
const REFERENCE_RECORD: &str = include_str!("../dry-air-20c.gas");
const MAX_RECORD_BYTES: u64 = 8192;

/// Which declaration supplied a run's air.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Source {
    /// The imported description's own `ambient` row.
    Description,
    /// The committed reference record of a built-in instrument.
    Reference,
    /// `--gas-record FILE`.
    CommandLine,
}

/// A resolved surrounding gas and the declaration it came from.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Ambient {
    pub gas: GasState,
    pub source: Source,
}

impl Ambient {
    /// The committed reference record, read through the gas owner.
    pub fn reference() -> Self {
        let gas = parse(REFERENCE_RECORD)
            .expect("the committed examples/dry-air-20c.gas record resolves");
        Self { gas, source: Source::Reference }
    }

    /// A description's `ambient` row: one complete shared gas declaration.
    pub fn described(record: &str) -> Result<Self, Error> {
        Ok(Self { gas: parse(record)?, source: Source::Description })
    }

    /// Exterior BEM medium: this gas's density and sound speed.
    pub fn medium(self) -> Medium {
        Medium { density: self.gas.density, sound_speed: self.gas.sound_speed }
    }

    /// Interior cavity medium: the same density and sound speed.
    pub fn acoustic_medium(self) -> AcousticMedium {
        AcousticMedium { rho0: self.gas.density, c0: self.gas.sound_speed }
    }

    /// Adiabatic bulk modulus rho c^2 (gamma p for an ideal gas) [Pa].
    pub fn bulk_modulus_pa(self) -> f64 {
        self.gas.density * self.gas.sound_speed * self.gas.sound_speed
    }

    /// One-line receipt naming the declaration and the resolved state.
    pub fn receipt(self) -> String {
        let source = match self.source {
            Source::Description => "description ambient row",
            Source::Reference => "built-in reference record examples/dry-air-20c.gas",
            Source::CommandLine => "--gas-record",
        };
        format!("ambient: {source}; T={} K, p={} Pa, rho={} kg/m3, c={} m/s, water mole fraction {}",
            self.gas.temperature, self.gas.pressure, self.gas.density, self.gas.sound_speed,
            self.gas.water_mole_fraction)
    }
}

/// Exactly one shared gas declaration and nothing else.
pub fn parse(record: &str) -> Result<GasState, Error> {
    let mut words = record.split_ascii_whitespace();
    let gas = read_gas_record(&mut words, None, None).map_err(|e| format!("ambient gas refused: {e}"))?;
    if words.next().is_some() {
        return Err("ambient gas record has trailing fields".into());
    }
    Ok(*gas.state())
}

/// Remove `--gas-record FILE`, returning its resolved gas (if supplied).
pub fn option(args: &mut Vec<String>) -> Result<Option<Ambient>, Error> {
    let mut positions = args.iter().enumerate().filter(|(_, s)| s.as_str() == "--gas-record");
    let Some((start, _)) = positions.next() else { return Ok(None) };
    if positions.next().is_some() { return Err("--gas-record may be supplied only once".into()); }
    let path = args.get(start + 1).filter(|s| !s.is_empty() && !s.starts_with("--"))
        .ok_or("--gas-record requires a gas declaration file")?.clone();
    let mut text = String::new();
    std::fs::File::open(&path).map_err(|e| format!("--gas-record {path}: {e}"))?
        .take(MAX_RECORD_BYTES + 1).read_to_string(&mut text)
        .map_err(|e| format!("--gas-record {path}: {e}"))?;
    if text.len() as u64 > MAX_RECORD_BYTES { return Err("--gas-record exceeds 8 KiB".into()); }
    let gas = parse(&text).map_err(|e| format!("--gas-record {path}: {e}"))?;
    args.drain(start..start + 2);
    Ok(Some(Ambient { gas, source: Source::CommandLine }))
}

/// `--gas-record` replaces a description's or the reference's air.
pub fn resolve(declared: Ambient, cli: Option<Ambient>) -> Ambient { cli.unwrap_or(declared) }

static COMMAND_LINE: OnceLock<Ambient> = OnceLock::new();

/// Record this process's `--gas-record` override; refuses a second one.
pub fn set_command_line(cli: Ambient) -> Result<(), Error> {
    COMMAND_LINE.set(cli).map_err(|_| "--gas-record was already declared for this run".into())
}

/// The air an instrument uses: its declared air unless the command line overrode it.
pub fn effective(declared: Ambient) -> Ambient { resolve(declared, COMMAND_LINE.get().copied()) }

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_reference_is_the_committed_record_through_the_gas_owner() {
        let reference = Ambient::reference();
        assert_eq!(reference.gas, GasState::try_new_moist_air(293.15, 101_325.0, 0.0).unwrap());
        assert_eq!(reference.source, Source::Reference);
        assert!((reference.bulk_modulus_pa() / (reference.gas.gamma * reference.gas.pressure) - 1.0).abs() < 1e-12);
        assert!(reference.receipt().contains("built-in reference record"));
    }

    #[test]
    fn records_that_do_not_resolve_refuse_and_never_become_air() {
        for bad in ["", "293.15 101325", "293.15 101325 0 7", "mixture 293.15 101325 mole 1 Xe 1",
            "thermochemical 293.15 101325 mole 1 N2 1", "293.15 0 0"] {
            assert!(parse(bad).is_err(), "{bad}");
            assert!(Ambient::described(bad).is_err(), "{bad}");
        }
        let co2 = Ambient::described("mixture 293.15 101325 mole 1 CO2 1").unwrap();
        assert!(co2.gas.sound_speed < 0.8 * Ambient::reference().gas.sound_speed);
        assert_eq!(co2.source, Source::Description);
    }

    #[test]
    fn the_command_line_record_overrides_and_is_removed_from_the_arguments() {
        let dir = std::env::temp_dir().join(format!("fs-percussion-ambient-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("thin.gas");
        std::fs::write(&path, "293.15 70927.5 0\n").unwrap();
        let mut args: Vec<String> = ["drum", "64", "--gas-record", path.to_str().unwrap()]
            .iter().map(|s| (*s).to_string()).collect();
        let cli = option(&mut args).unwrap().unwrap();
        assert_eq!(args, ["drum", "64"]);
        assert_eq!(cli.source, Source::CommandLine);
        assert_eq!(resolve(Ambient::reference(), Some(cli)), cli);
        assert_eq!(resolve(Ambient::reference(), None), Ambient::reference());
        let mut twice: Vec<String> = ["--gas-record", "a", "--gas-record", "b"].iter().map(|s| (*s).to_string()).collect();
        assert!(option(&mut twice).is_err());
        std::fs::write(&path, "293.15 101325\n").unwrap();
        let mut bad: Vec<String> = ["--gas-record", path.to_str().unwrap()].iter().map(|s| (*s).to_string()).collect();
        assert!(option(&mut bad).is_err());
        assert_eq!(bad.len(), 2, "a refused record leaves the arguments intact");
    }
}
