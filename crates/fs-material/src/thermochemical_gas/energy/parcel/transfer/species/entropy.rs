//! Entropy differences for closed, fixed-volume selective gas transport.
//!
//! A constant entropy datum cancels for each conserved source-qualified species.
//! Compute that cancellation before rounding. Species inventory residuals remain
//! in the transport owner's mass accounts; their arbitrary reference entropy is
//! deliberately not interpreted as physical entropy generation. No reactions,
//! external mass sources or moving boundaries are admitted by this helper.

use super::super::snapshot::MAX_SNAPSHOT_PARCELS;
use super::{ClosedGasParcel, GasEnergyError, SpeciesId, source_error};
use fs_qty::Temperature;
use fs_thermochem::{Nasa9StandardStateModelV1, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};
use std::collections::BTreeMap;

#[derive(Debug)]
pub(super) enum Error {
    Invalid(&'static str),
    Parcel { id: u64, error: GasEnergyError },
    Cancelled,
}
fn check(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
#[derive(Default)]
struct Accumulator {
    value: f64,
    correction: f64,
}
impl Accumulator {
    fn add(&mut self, term: f64) -> Result<(), Error> {
        let next = self.value + term;
        let lost = if self.value.abs() >= term.abs() {
            (self.value - next) + term
        } else {
            (term - next) + self.value
        };
        self.correction += lost;
        self.value = next;
        if ![term, self.value, self.correction]
            .iter()
            .all(|x| x.is_finite())
        {
            return Err(Error::Invalid(
                "closed gas entropy arithmetic is unrepresentable",
            ));
        }
        Ok(())
    }
    fn total(self) -> Result<f64, Error> {
        let value = self.value + self.correction;
        if value.is_finite() {
            Ok(value)
        } else {
            Err(Error::Invalid("closed gas entropy sum is unrepresentable"))
        }
    }
}
struct Reference<'a> {
    model: &'a Nasa9StandardStateModelV1,
    temperature_k: f64,
    volume_m3: f64,
    amount_mol: f64,
    mass_change: Accumulator,
}
// Source cards are immutable content identities; no material-law approximation.
fn same_source(a: &Nasa9StandardStateModelV1, b: &Nasa9StandardStateModelV1) -> bool {
    a.species() == b.species()
        && a.molar_mass().value().to_bits() == b.molar_mass().value().to_bits()
        && a.convention() == b.convention()
        && a.regions()
            .iter()
            .map(|r| r.card_identity())
            .eq(b.regions().iter().map(|r| r.card_identity()))
}
fn amount(parcel: &ClosedGasParcel, index: usize) -> Result<f64, Error> {
    let mixture = parcel.model().thermodynamics();
    let molar_mass = parcel
        .state()
        .thermodynamics()
        .properties()
        .molar_mass()
        .value();
    let fraction = mixture.mole_composition().fractions()[index];
    let value = (parcel.mass_kg() / molar_mass) * fraction;
    if value.is_finite() && value >= 0.0 && (fraction == 0.0 || value > 0.0) {
        Ok(value)
    } else {
        Err(Error::Invalid(
            "closed gas species amount vanished or overflowed",
        ))
    }
}
// ln(a/b), including nearly equal positive values and overflowing raw ratios.
fn log_ratio(a: f64, b: f64) -> f64 {
    if a == b {
        return 0.0;
    }
    let (lo, hi, sign) = if a > b { (b, a, 1.0) } else { (a, b, -1.0) };
    let delta = (hi - lo) / lo;
    let rounded = 1.0 + delta;
    sign * if rounded == 1.0 {
        delta
    } else if rounded.is_finite() {
        fs_math::det::ln(rounded) * (delta / (rounded - 1.0))
    } else {
        fs_math::det::ln(hi) - fs_math::det::ln(lo)
    }
}
// b*ln(b/a)-(b-a). The small-increment series avoids cancellation of the
// linear term. With |z|<=1/8, 22 terms put truncation below binary64 roundoff.
fn mixing_remainder(a: f64, b: f64) -> f64 {
    if a == b {
        return 0.0;
    }
    let z = (b - a) / a;
    if z.abs() <= 0.125 {
        let mut polynomial = 1.0 / (22.0 * 21.0);
        for k in (2..22).rev() {
            polynomial = 1.0 / f64::from(k * (k - 1)) - z * polynomial;
        }
        // Multiply in this order so z*z need not be representable by itself.
        ((a * z) * z) * polynomial
    } else {
        b * log_ratio(b, a) - (b - a)
    }
}
fn thermal_change(
    model: &Nasa9StandardStateModelV1,
    a: f64,
    b: f64,
    id: u64,
) -> Result<f64, Error> {
    model
        .molar_entropy_change(Temperature::new(a), Temperature::new(b))
        .map(|s| s.value())
        .map_err(|e| Error::Parcel {
            id,
            error: source_error(e),
        })
}
// Relative partial entropy, with its constant reference canceled algebraically.
// Partial-pressure ratios use n*T/V, avoiding an underflowing x*p product.
fn potential(reference: &Reference<'_>, n: f64, t: f64, v: f64, id: u64) -> Result<f64, Error> {
    let logarithm = log_ratio(n, reference.amount_mol) + log_ratio(t, reference.temperature_k)
        - log_ratio(v, reference.volume_m3);
    Ok(thermal_change(reference.model, reference.temperature_k, t, id)? - R * logarithm)
}
fn contributions(
    reference: &Reference<'_>,
    n0: f64,
    n1: f64,
    t0: f64,
    t1: f64,
    volume: f64,
    id: u64,
) -> Result<[f64; 3], Error> {
    if n0 == n1 && t0 == t1 {
        return Ok([0.0; 3]);
    }
    if n0 == 0.0 {
        return Ok([
            n1 * (potential(reference, n1, t1, volume, id)? + R),
            0.0,
            0.0,
        ]);
    }
    if n1 == 0.0 {
        return Ok([
            -n0 * (potential(reference, n0, t0, volume, id)? + R),
            0.0,
            0.0,
        ]);
    }
    let dn = n1 - n0;
    let caloric = thermal_change(reference.model, t0, t1, id)? - R * log_ratio(t1, t0);
    let carried = if dn == 0.0 {
        0.0
    } else {
        dn * potential(reference, n0, t0, volume, id)?
    };
    Ok([n1 * caloric, carried, -R * mixing_remainder(n0, n1)])
}

/// Same closed species inventory, same fixed volumes, original source models.
/// The report cancels the entropy datum of every conserved species. The bounded
/// species-mass residual is checked separately and is not multiplied by that
/// arbitrary datum. This is endpoint entropy, not an unconditional time theorem.
pub(super) fn closed_entropy_change<'a>(
    states: impl IntoIterator<Item = (u64, &'a ClosedGasParcel, &'a ClosedGasParcel)>,
    mass_tolerance_kg: f64,
    mut cancelled: impl FnMut() -> bool,
) -> Result<f64, Error> {
    check(&mut cancelled)?;
    if !mass_tolerance_kg.is_finite() || mass_tolerance_kg <= 0.0 {
        return Err(Error::Invalid(
            "positive closed-species mass tolerance required",
        ));
    }
    let mut rows = Vec::new();
    for row in states {
        check(&mut cancelled)?;
        if rows.len() == MAX_SNAPSHOT_PARCELS || row.1.volume_m3() != row.2.volume_m3() {
            return Err(Error::Invalid(
                "bounded entropy audit requires unchanged parcel volumes",
            ));
        }
        rows.push(row);
    }
    let mut references: BTreeMap<SpeciesId, Vec<Reference<'a>>> = BTreeMap::new();
    for &(_, old, _) in &rows {
        for (index, model) in old.model().thermodynamics().components().iter().enumerate() {
            check(&mut cancelled)?;
            let n = amount(old, index)?;
            if n == 0.0 {
                continue;
            }
            let group = references.entry(model.species().clone()).or_default();
            if !group.iter().any(|r| same_source(r.model, model)) {
                group.push(Reference {
                    model,
                    temperature_k: old.state().state().temperature,
                    volume_m3: old.volume_m3(),
                    amount_mol: n,
                    mass_change: Accumulator::default(),
                });
            }
        }
    }
    let mut entropy = Accumulator::default();
    for &(id, old, next) in &rows {
        let before = old.model().thermodynamics().components();
        let after = next.model().thermodynamics().components();
        let (mut i, mut j) = (0, 0);
        while i < before.len() || j < after.len() {
            check(&mut cancelled)?;
            let ordering = match (before.get(i), after.get(j)) {
                (Some(a), Some(b)) => a.species().cmp(b.species()),
                (Some(_), None) => std::cmp::Ordering::Less,
                (None, Some(_)) => std::cmp::Ordering::Greater,
                (None, None) => unreachable!(),
            };
            let (model, n0, n1) = match ordering {
                std::cmp::Ordering::Less => {
                    let result = (&before[i], amount(old, i)?, 0.0);
                    i += 1;
                    result
                }
                std::cmp::Ordering::Greater => {
                    let result = (&after[j], 0.0, amount(next, j)?);
                    j += 1;
                    result
                }
                std::cmp::Ordering::Equal => {
                    if !same_source(&before[i], &after[j]) {
                        return Err(Error::Invalid("closed species entropy source changed"));
                    }
                    let result = (&before[i], amount(old, i)?, amount(next, j)?);
                    i += 1;
                    j += 1;
                    result
                }
            };
            if n0 == 0.0 && n1 == 0.0 {
                continue;
            }
            let reference = references
                .get_mut(model.species())
                .and_then(|group| group.iter_mut().find(|r| same_source(r.model, model)))
                .ok_or(Error::Invalid(
                    "closed entropy endpoint introduced a new species source",
                ))?;
            for term in contributions(
                reference,
                n0,
                n1,
                old.state().state().temperature,
                next.state().state().temperature,
                old.volume_m3(),
                id,
            )? {
                entropy.add(term)?;
            }
            reference
                .mass_change
                .add((n1 - n0) * model.molar_mass().value())?;
        }
    }
    for group in references.into_values() {
        for reference in group {
            check(&mut cancelled)?;
            if reference.mass_change.total()?.abs() > mass_tolerance_kg {
                return Err(Error::Invalid(
                    "closed entropy audit did not conserve each source-qualified species",
                ));
            }
        }
    }
    check(&mut cancelled)?;
    entropy.total()
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn g0_mixing_remainder_keeps_the_positive_quadratic_signal() {
        let z = 2.0_f64.powi(-40);
        for scale in [1e-200, 1.0, 1e200] {
            let a = scale;
            let b = scale * (1.0 + z);
            let actual_z = (b - a) / a;
            let expected = ((a * actual_z) * actual_z) * 0.5;
            let remainder = mixing_remainder(a, b);
            assert!(remainder > 0.0);
            assert!((remainder / expected - 1.0).abs() < actual_z);
        }
        assert_eq!(mixing_remainder(1.0, 1.0), 0.0);
    }
}
