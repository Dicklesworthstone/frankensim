//! Caller-supplied per-key hammer materials for the existing contact solver.
//!
//! Format (SI):
//! ```text
//! frankensim-hammer-materials-v1
//! felt,69,400000,0.2,2.5,3.2,0.25,0.8,2500000
//! branch,69,2000000,0.0002
//! branch,69,500000,0.004
//! ```
//! felt: key, reference stress [Pa], reference strain, loading exponent,
//! unloading exponent, crush fraction, densification strain, equilibrium
//! Prony modulus [Pa]. branch: key, modulus [Pa], relaxation time [s].
//! Branches follow their felt record; an explicit felt with no branches is
//! legal. Every admitted key needs exactly one card, including unplayed keys.
//! The original creep realization replaces its instantaneous spring with the
//! unilateral felt envelope; this is NOT a parallel stress sum or an output EQ.
//! Importing declared parameters does not certify a coupon fit, bandwidth,
//! source licence, calibrated instrument or agreement with measurements.

use std::collections::{BTreeMap, BTreeSet};
use fs_material::{WoolFelt, visco::GeneralizedMaxwell};

pub type Material = (WoolFelt, GeneralizedMaxwell);
pub const HEADER: &str = "frankensim-hammer-materials-v1";

/// Cold input only; constitutive validation remains with fs-material.
/// Return cards in COURSE order, not numeric-key or file order.
pub fn read(text: &str, keys: &[u8]) -> Result<Vec<Material>, String> {
    let wanted: BTreeSet<_> = keys.iter().copied().collect();
    if keys.is_empty() || keys.len() > 88 || wanted.len() != keys.len()
        || keys.iter().any(|k| !(21..=108).contains(k)) {
        return Err("hammer materials need distinct admitted keys in 21..=108".into());
    }
    let mut cards: BTreeMap<u8, Material> = BTreeMap::new();
    let mut header = false;
    for (index, raw) in text.lines().enumerate() {
        let line = index + 1;
        let row = raw.split('#').next().unwrap_or("").trim();
        if row.is_empty() { continue; }
        let error = |why: &str| format!("hammer material line {line}: {why}");
        if !header {
            if row != HEADER { return Err(error("expected frankensim-hammer-materials-v1")); }
            header = true;
            continue;
        }
        let fields: Vec<_> = row.split(',').map(str::trim).collect();
        let expected = match fields[0] {
            "felt" => 9,
            "branch" => 4,
            _ => return Err(error("unknown record; expected felt or branch")),
        };
        if fields.len() != expected { return Err(error("wrong number of fields")); }
        let key: u8 = fields[1].parse().map_err(|_| error("invalid key"))?;
        if !wanted.contains(&key) { return Err(error("key is absent from the supplied string scale")); }
        let scalar = |column: usize| -> Result<f64, String> {
            let x: f64 = fields[column].parse().map_err(|_| error("invalid scalar"))?;
            if !x.is_finite() { return Err(error("nonfinite material parameter")); }
            Ok(x)
        };
        if fields[0] == "felt" {
            if cards.contains_key(&key) { return Err(error("duplicate felt card")); }
            let law = WoolFelt::new(scalar(2)?, scalar(3)?, scalar(4)?, scalar(5)?,
                scalar(6)?, scalar(7)?).map_err(|e| error(&e.to_string()))?;
            let e_inf = scalar(8)?;
            if e_inf <= 0.0 { return Err(error("creep solid needs positive equilibrium modulus")); }
            let prony = GeneralizedMaxwell::new(e_inf, Vec::new())
                .map_err(|e| error(&e.to_string()))?;
            cards.insert(key, (law, prony));
        } else {
            let modulus = scalar(2)?;
            let tau = scalar(3)?;
            let (_, prony) = cards.get_mut(&key)
                .ok_or_else(|| error("branch must follow its key's felt card"))?;
            if prony.terms.len() >= 8 { return Err(error("at most eight Prony branches per key")); }
            let mut terms = prony.terms.clone();
            terms.push((modulus, tau));
            let next = GeneralizedMaxwell::new(prony.e_inf, terms)
                .map_err(|e| error(&e.to_string()))?;
            if !(next.e_inf + next.terms.iter().map(|t| t.0).sum::<f64>()).is_finite() {
                return Err(error("instantaneous modulus overflow"));
            }
            *prony = next;
        }
    }
    if !header { return Err("missing hammer material header".into()); }
    keys.iter().map(|key| cards.remove(key)
        .ok_or_else(|| format!("missing hammer material for key {key}; no fallback is permitted")))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    fn card(key: u8, stress: f64) -> String {
        format!("felt,{key},{stress},0.2,2.5,3.2,0.25,0.8,2500000\n")
    }
    #[test]
    fn material_identity_and_course_order_survive_import() {
        let text = format!("# declared SI parameters\n{HEADER}\n{}branch,60,2000000,0.0002\n{}",
            card(60, 300000.0), card(69, 400000.0));
        let cards = read(&text, &[69, 60]).unwrap();
        assert_eq!(cards[0].0.f_ref, 400000.0);
        assert!(cards[0].1.terms.is_empty());
        assert_eq!(cards[1].0.f_ref, 300000.0);
        assert_eq!(cards[1].1.e_inf, 2500000.0);
        assert_eq!(cards[1].1.terms, vec![(2000000.0, 0.0002)]);
    }
    #[test]
    fn incomplete_ambiguous_or_invalid_materials_never_select_defaults() {
        let valid = format!("{HEADER}\n{}", card(69, 400000.0));
        assert!(read(&valid, &[60, 69]).is_err());
        assert!(read(&valid, &[69, 69]).is_err());
        assert!(read(&valid, &[60]).is_err());
        assert!(read(&valid, &[]).is_err());
        for text in [String::new(), card(69, 400000.0),
            format!("{valid}{}", card(69, 500000.0)),
            format!("{valid}branch,69,-1,0.001\n"),
            format!("{valid}branch,69,1,0\n"),
            format!("{valid}branch,69,NaN,0.001\n"),
            format!("{valid}branch,69,1,inf\n"),
            format!("{valid}unknown,69,1,2\n"),
            format!("{HEADER}\nbranch,69,1,0.001\n{}", card(69, 400000.0)),
            format!("{valid}{}", "branch,69,1,0.001\n".repeat(9)),
            valid.replace("2500000", "0"), valid.replace("400000", "NaN")] {
            assert!(read(&text, &[69]).is_err(), "accepted {text}");
        }
    }
    #[test]
    fn imported_materials_change_contact_mechanics_not_the_string_scale() {
        use crate::{board, engine, geometry};
        let course = geometry::demonstration_scale().unwrap()[48];
        let build = |stress| {
            let text = format!("{HEADER}\n{}branch,69,2000000,0.0002\nbranch,69,500000,0.004\n", card(69, stress));
            engine::Instrument::new_with_course_felts(vec![course], &board::demonstration(),
                48_000, 4, 12, true, read(&text, &[69]).unwrap()).unwrap()
        };
        let mut soft = build(300000.0);
        let mut hard = build(600000.0);
        assert_eq!(soft.bank.q, hard.bank.q);
        assert_eq!(soft.bank.v, hard.bank.v);
        assert_eq!(soft.bank.modes.iter().map(|m| m.omega).collect::<Vec<_>>(),
            hard.bank.modes.iter().map(|m| m.omega).collect::<Vec<_>>());
        soft.note_on(69, 2.0).unwrap(); hard.note_on(69, 2.0).unwrap();
        let mut changed = 0.0_f64;
        for _ in 0..1500 {
            let a = soft.step().unwrap(); let b = hard.step().unwrap();
            assert!(a.is_finite() && b.is_finite());
            changed = changed.max((a - b).abs());
        }
        assert!(changed > 1e-12, "material must change the physical board response");
        for piano in [&soft, &hard] {
            assert!(piano.accounting.felt_loss_j > 0.0);
            assert!(piano.accounting.felt_relaxation_loss_j > 0.0);
            assert!((piano.accounting.input_work_j - piano.energy_j()
                - piano.accounting.dissipated_j()).abs() < 1e-7);
        }
    }
}
