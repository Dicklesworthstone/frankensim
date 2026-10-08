//! Cold export of the existing microphone coordinate map. No mechanics step,
//! observer, eigenanalysis or alternate pressure projection is introduced.
use super::linear::Bank;
use std::io::Write;

const DOMAIN: &str = "fs.piano.loaded-board-pressure-basis.v1";

pub fn write(out: &mut impl Write, bank: &Bank) -> Result<(), String> {
    let (basis, omega2) = bank.board_pressure_basis()?;
    let count = bank.board_count;
    let dimension = u64::try_from(count)
        .map_err(|e| e.to_string())?
        .to_le_bytes();
    // Protocol: rows and columns as little-endian u64, row-major IEEE-754 f64
    // map entries, then the loaded diagonal omega-squared entries, all LE.
    let mut hash = fs_blake3::DomainHasher::new(DOMAIN);
    hash.update(&dimension);
    hash.update(&dimension);
    for value in basis.iter().chain(omega2) {
        hash.update(&value.to_le_bytes());
    }
    writeln!(
        out,
        "{{\n  \"format\": \"{DOMAIN}\",\n  \"basis_id\": \"{}\",",
        hash.finalize()
    )
    .map_err(|e| e.to_string())?;
    writeln!(out, "  \"bare_mode_count\": {count},\n  \"loaded_mode_count\": {count},\n  \"mechanical_sample_rate_hz\": {},",
        bank.rate).map_err(|e| e.to_string())?;
    writeln!(out, "  \"map_order\": \"row-major: bare row, loaded column\",\n  \"pressure_mode_order\": \"loaded column index, zero-based\",\n  \"frequency_role\": \"split-step diagonal reference; not full coupled instrument poles\",")
        .map_err(|e| e.to_string())?;
    write!(out, "  \"bare_from_loaded\": [").map_err(|e| e.to_string())?;
    numbers(out, basis)?;
    write!(out, "],\n  \"loaded_reference_omega_squared_rad2_s2\": [")
        .map_err(|e| e.to_string())?;
    numbers(out, omega2)?;
    write!(out, "],\n  \"loaded_reference_frequency_hz\": [").map_err(|e| e.to_string())?;
    for (i, value) in omega2.iter().enumerate() {
        if i != 0 {
            write!(out, ",").map_err(|e| e.to_string())?;
        }
        write!(
            out,
            "{:.17e}",
            fs_math::det::sqrt(*value) / std::f64::consts::TAU
        )
        .map_err(|e| e.to_string())?;
    }
    writeln!(out, "]\n}}").map_err(|e| e.to_string())?;
    out.flush().map_err(|e| e.to_string())
}

fn numbers(out: &mut impl Write, values: &[f64]) -> Result<(), String> {
    for (i, value) in values.iter().enumerate() {
        if i != 0 {
            write!(out, ",").map_err(|e| e.to_string())?;
        }
        write!(out, "{value:.17e}").map_err(|e| e.to_string())?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn repeated_export_is_byte_stable_and_preserves_the_prepared_bank() {
        let scale = super::super::geometry::demonstration_scale().unwrap();
        let bank = Bank::new(
            &[scale[48]],
            &super::super::board::demonstration(),
            192_000,
            21_600.0,
            12,
            true,
        )
        .unwrap();
        let before = (
            bank.q.clone(),
            bank.v.clone(),
            bank.energy(),
            bank.volume_velocity(),
        );
        let mut first = Vec::new();
        let mut second = Vec::new();
        write(&mut first, &bank).unwrap();
        write(&mut second, &bank).unwrap();
        assert_eq!(first, second);
        assert_eq!(
            before,
            (
                bank.q.clone(),
                bank.v.clone(),
                bank.energy(),
                bank.volume_velocity()
            )
        );
        let text = String::from_utf8(first).unwrap();
        assert!(text.contains("\"loaded_mode_count\": 4"));
        assert!(text.contains("not full coupled instrument poles"));
        assert!(text.contains("\"loaded_reference_omega_squared_rad2_s2\""));
        let array = |field: &str| -> Vec<f64> {
            text.split(&format!("\"{field}\": ["))
                .nth(1)
                .unwrap()
                .split(']')
                .next()
                .unwrap()
                .split(',')
                .map(|value| value.parse().unwrap())
                .collect()
        };
        let map = array("bare_from_loaded");
        let omega2 = array("loaded_reference_omega_squared_rad2_s2");
        let (expected_map, expected_omega2) = bank.board_pressure_basis().unwrap();
        assert_eq!(map, expected_map);
        assert_eq!(omega2, expected_omega2);
        let frequencies = array("loaded_reference_frequency_hz");
        assert_eq!(
            frequencies,
            omega2
                .iter()
                .map(|value| fs_math::det::sqrt(*value) / std::f64::consts::TAU)
                .collect::<Vec<_>>()
        );
        // Rebuild the documented identity preimage from the emitted values,
        // using the existing one-shot owner rather than the writer's streaming path.
        let dimension = u64::try_from(bank.board_count).unwrap().to_le_bytes();
        let mut payload = [dimension, dimension].concat();
        for value in map.iter().chain(&omega2) {
            payload.extend(value.to_le_bytes());
        }
        let expected_id = fs_blake3::hash_domain(DOMAIN, &payload);
        assert!(text.contains(&format!("\"basis_id\": \"{expected_id}\"")));
    }
}
