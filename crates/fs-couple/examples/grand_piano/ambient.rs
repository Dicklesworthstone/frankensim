//! The air the piano radiates into. It is DATA, never a constant in code. The
//! built-in instruments (the demonstration and preset pianos, and the Model D
//! strings around an imported board) load the committed record
//! `examples/moist-air-20c-65rh.gas` (20 C, 101325 Pa, 65 % RH; the water mole
//! fraction follows from the gas owner's Buck-1996 saturation law) through the
//! gas owner. The Model D soundboard reads the SAME record for its wood's
//! moisture content (`steinway_d::BoardEnvironment`). `--gas-record FILE`
//! overrides it; a `thermochemical ...` record resolves only through supplied
//! `--gas-component` source packs. Every report names which of the two
//! supplied the air.
use fs_bem::helmholtz::Medium;
use fs_material::gas::GasState;
use fs_material::gas_species::input::{ThermochemicalSources, read_gas_record};
use std::io::Read;

/// The committed reference record shared by every built-in instrument.
const REFERENCE_RECORD: &str = include_str!("../moist-air-20c-65rh.gas");
const MAX_RECORD_BYTES: u64 = 8192;

/// Which declaration supplied the air.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Source {
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
            .expect("the committed examples/moist-air-20c-65rh.gas record resolves");
        Self { gas, source: Source::Reference }
    }

    /// `--gas-record FILE`: exactly one shared gas declaration.
    pub fn load(path: &str) -> Result<Self, String> { Self::load_with(path, None) }

    /// `--gas-record FILE`, resolving a thermochemical record through `sources`.
    pub fn load_with(path: &str, sources: Option<&ThermochemicalSources<'_>>) -> Result<Self, String> {
        let mut text = String::new();
        std::fs::File::open(path).map_err(|e| format!("--gas-record {path}: {e}"))?
            .take(MAX_RECORD_BYTES + 1).read_to_string(&mut text)
            .map_err(|e| format!("--gas-record {path}: {e}"))?;
        if text.len() as u64 > MAX_RECORD_BYTES { return Err("--gas-record exceeds 8 KiB".into()); }
        let gas = parse_with(&text, sources).map_err(|e| format!("--gas-record {path}: {e}"))?;
        Ok(Self { gas, source: Source::CommandLine })
    }

    /// `--gas-record FILE` when given, otherwise the reference record.
    pub fn selected(path: Option<&str>) -> Result<Self, String> { Self::selected_with(path, None) }

    /// As `selected`, resolving a thermochemical record through `sources`.
    pub fn selected_with(path: Option<&str>, sources: Option<&ThermochemicalSources<'_>>)
        -> Result<Self, String> {
        path.map_or_else(|| Ok(Self::reference()), |path| Self::load_with(path, sources))
    }

    /// Exterior medium: this gas's density and sound speed.
    pub fn medium(self) -> Medium {
        Medium { density: self.gas.density, sound_speed: self.gas.sound_speed }
    }

    /// One-line report naming the declaration and the resolved state.
    pub fn receipt(self) -> String {
        let source = match self.source {
            Source::Reference => "built-in reference record examples/moist-air-20c-65rh.gas",
            Source::CommandLine => "--gas-record",
        };
        format!("Ambient air: {source}; T={} K, p={} Pa, rho={} kg/m3, c={} m/s, water mole fraction {}.",
            self.gas.temperature, self.gas.pressure, self.gas.density, self.gas.sound_speed,
            self.gas.water_mole_fraction)
    }
}

/// Exactly one shared gas declaration and nothing else; never defaults to air.
pub fn parse(record: &str) -> Result<GasState, String> { parse_with(record, None) }

/// As `parse`; a thermochemical record resolves only through `sources`.
pub fn parse_with(record: &str, sources: Option<&ThermochemicalSources<'_>>) -> Result<GasState, String> {
    let mut words = record.split_ascii_whitespace();
    let gas = read_gas_record(&mut words, None, sources).map_err(|e| format!("ambient gas refused: {e}"))?;
    if words.next().is_some() { return Err("ambient gas record has trailing fields".into()); }
    Ok(*gas.state())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_reference_is_the_committed_humid_record_and_dry_air_is_the_former_room_air() {
        let reference = Ambient::reference();
        assert_eq!(reference.gas, GasState::try_new_moist_air(293.15, 101_325.0, 0.65).unwrap());
        assert!(reference.receipt().contains("built-in reference record examples/moist-air-20c-65rh.gas"));
        // 65 % RH air is lighter and faster than fs-bem's dry room-air medium...
        assert!(reference.medium().density < Medium::air().density);
        assert!(reference.medium().sound_speed > Medium::air().sound_speed);
        // ...which the explicit dry record still reproduces, bit for bit.
        let dry = Ambient { gas: parse("293.15 101325 0").unwrap(), source: Source::CommandLine };
        assert_eq!(dry.medium().density.to_bits(), Medium::air().density.to_bits());
        assert_eq!(dry.medium().sound_speed.to_bits(), Medium::air().sound_speed.to_bits());
    }

    #[test]
    fn a_record_that_does_not_resolve_refuses() {
        for bad in ["", "293.15 101325", "293.15 101325 0 7", "mixture 293.15 101325 mole 1 Xe 1", "293.15 0 0"] {
            assert!(parse(bad).is_err(), "{bad}");
        }
        assert!(Ambient::load("/nonexistent/air.gas").is_err());
        assert_eq!(Ambient::selected(None).unwrap(), Ambient::reference());
    }
}
