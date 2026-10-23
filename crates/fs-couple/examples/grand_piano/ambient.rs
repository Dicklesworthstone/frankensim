//! The air the piano radiates into. It is DATA, never a constant in code. The
//! built-in instruments (the demonstration and preset pianos, and the Model D
//! strings around an imported board) load the committed record
//! `examples/dry-air-20c.gas` through the gas owner. `--gas-record FILE`
//! overrides it. Every report names which of the two supplied the air.
use fs_bem::helmholtz::Medium;
use fs_material::gas::GasState;
use fs_material::gas_species::input::read_gas_record;
use std::io::Read;

/// The committed reference record shared by every built-in instrument.
const REFERENCE_RECORD: &str = include_str!("../dry-air-20c.gas");
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
            .expect("the committed examples/dry-air-20c.gas record resolves");
        Self { gas, source: Source::Reference }
    }

    /// `--gas-record FILE`: exactly one shared gas declaration.
    pub fn load(path: &str) -> Result<Self, String> {
        let mut text = String::new();
        std::fs::File::open(path).map_err(|e| format!("--gas-record {path}: {e}"))?
            .take(MAX_RECORD_BYTES + 1).read_to_string(&mut text)
            .map_err(|e| format!("--gas-record {path}: {e}"))?;
        if text.len() as u64 > MAX_RECORD_BYTES { return Err("--gas-record exceeds 8 KiB".into()); }
        let gas = parse(&text).map_err(|e| format!("--gas-record {path}: {e}"))?;
        Ok(Self { gas, source: Source::CommandLine })
    }

    /// `--gas-record FILE` when given, otherwise the reference record.
    pub fn selected(path: Option<&str>) -> Result<Self, String> {
        path.map_or_else(|| Ok(Self::reference()), Self::load)
    }

    /// Exterior medium: this gas's density and sound speed.
    pub fn medium(self) -> Medium {
        Medium { density: self.gas.density, sound_speed: self.gas.sound_speed }
    }

    /// One-line report naming the declaration and the resolved state.
    pub fn receipt(self) -> String {
        let source = match self.source {
            Source::Reference => "built-in reference record examples/dry-air-20c.gas",
            Source::CommandLine => "--gas-record",
        };
        format!("Ambient air: {source}; T={} K, p={} Pa, rho={} kg/m3, c={} m/s, water mole fraction {}.",
            self.gas.temperature, self.gas.pressure, self.gas.density, self.gas.sound_speed,
            self.gas.water_mole_fraction)
    }
}

/// Exactly one shared gas declaration and nothing else; never defaults to air.
pub fn parse(record: &str) -> Result<GasState, String> {
    let mut words = record.split_ascii_whitespace();
    let gas = read_gas_record(&mut words, None, None).map_err(|e| format!("ambient gas refused: {e}"))?;
    if words.next().is_some() { return Err("ambient gas record has trailing fields".into()); }
    Ok(*gas.state())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_reference_is_the_committed_record_and_matches_the_former_room_air() {
        let reference = Ambient::reference();
        assert_eq!(reference.gas, GasState::try_new_moist_air(293.15, 101_325.0, 0.0).unwrap());
        // fs-bem's room-air medium is the same dry-air state, bit for bit.
        assert_eq!(reference.medium().density.to_bits(), Medium::air().density.to_bits());
        assert_eq!(reference.medium().sound_speed.to_bits(), Medium::air().sound_speed.to_bits());
        assert!(reference.receipt().contains("built-in reference record"));
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
