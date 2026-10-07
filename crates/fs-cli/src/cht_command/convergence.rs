//! `cooling-cht` grid-convergence verification. A scene with a
//! `grid_convergence` block (`{"splits": [a, b, c]}`, three increasing
//! whole numbers) is solved three times, every cell of its declared grid
//! split into `a`, `b` and `c` equal parts per axis, and every result
//! quantity gets the procedure of Celik et al. (2008, J. Fluids Eng. 130,
//! 078001; ASME V&V 20): the observed order `p` from the three solutions
//! (a fixed-point iteration when the two refinement ratios differ), the
//! Richardson extrapolation, and the fine-grid convergence index
//! `GCI = 1.25 |phi_1 - phi_2| / (r_21^p - 1)` (Roache's three-grid safety
//! factor), reported in the quantity's own units.
//!
//! No-claims: the GCI is an asymptotic error band, not a bound, and means
//! something only when the grids are in the asymptotic range (monotone
//! convergence at an order near the schemes'); oscillatory, divergent and
//! indeterminate quantities get no band (their spread is reported); the
//! order used for the band and the extrapolation is capped at 2, the
//! schemes' formal order (a higher observed order is pre-asymptotic and
//! would shrink the band); geometry not aligned with the declared grid
//! (STL solids, off-grid boxes) changes its staircase between levels,
//! which the study measures but cannot separate from discretization
//! error; Monte Carlo radiation factors re-sample on every level; every
//! level inherits the `cooling-cht` no-claims.

use std::fmt::Write as _;
use std::time::Instant;

use super::json::JsonValue as J;
use super::study::{json_text, parallel_map, quantity, quantity_names, workers};
use super::{Failure, NO_CLAIM, Result, Scene, bad, execute_within_budget, num, quote};

const CONVERGENCE_SCHEMA: &str = "frankensim.cooling-cht.convergence.v1";
/// Roache's safety factor for a three-grid study.
const SAFETY_FACTOR: f64 = 1.25;
/// The schemes' formal order: the cap on the order used for the band.
const FORMAL_ORDER: f64 = 2.0;
/// The largest split per axis (the cell cap still applies).
const MAX_SPLIT: usize = 16;

/// How a quantity behaves over the three grids (Celik's ratio
/// `R = (phi_2 - phi_1) / (phi_3 - phi_2)`, 1 the finest grid).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum Behaviour {
    /// `0 < R < 1`: the band applies.
    Monotone,
    /// `R < 0`: the solutions alternate.
    Oscillatory,
    /// `R >= 1`: the change grows with refinement.
    Divergent,
    /// All three agree to round-off.
    Converged,
    /// The two finer agree to round-off but the coarse does not.
    Indeterminate,
}

impl Behaviour {
    fn name(self) -> &'static str {
        match self {
            Self::Monotone => "monotone",
            Self::Oscillatory => "oscillatory",
            Self::Divergent => "divergent",
            Self::Converged => "converged",
            Self::Indeterminate => "indeterminate",
        }
    }
}

/// The convergence analysis of one quantity.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct Analysis {
    pub behaviour: Behaviour,
    /// Observed order (monotone only).
    pub observed_order: Option<f64>,
    /// The order used: observed, capped at the formal order.
    pub order_used: Option<f64>,
    /// Richardson extrapolation from the two finer grids.
    pub extrapolated: Option<f64>,
    /// Fine-grid convergence index, in the quantity's units.
    pub gci: Option<f64>,
    /// Largest minus smallest of the three values.
    pub spread: f64,
}

/// Celik's observed order for refinement ratios `r21 = h2/h1`,
/// `r32 = h3/h2` and changes `e21 = phi2 - phi1`, `e32 = phi3 - phi2`:
/// `p = |ln|e32/e21| + q(p)| / ln r21`, `q(p) = ln((r21^p - s) /
/// (r32^p - s))`, `s = sign(e32/e21)`, by fixed-point iteration from
/// `q = 0` (exact in one step for equal ratios).
pub(super) fn observed_order(r21: f64, r32: f64, e21: f64, e32: f64) -> Option<f64> {
    let ratio = e32 / e21;
    let s = ratio.signum();
    let mut p = ratio.abs().ln().abs() / r21.ln();
    for _ in 0..200 {
        let q = ((r21.powf(p) - s) / (r32.powf(p) - s)).ln();
        let next = (ratio.abs().ln() + q).abs() / r21.ln();
        if !next.is_finite() {
            return None;
        }
        if (next - p).abs() <= 1e-12 * p.max(1.0) {
            return Some(next);
        }
        p = next;
    }
    None
}

/// Analyse `values = [coarse, medium, fine]` solved with cell splits
/// `splits = [a, b, c]` (the grid spacing scales as `1 / split`).
pub(super) fn analyse(values: [f64; 3], splits: [usize; 3]) -> Analysis {
    let [phi3, phi2, phi1] = values;
    let r21 = splits[2] as f64 / splits[1] as f64;
    let r32 = splits[1] as f64 / splits[0] as f64;
    let (e21, e32) = (phi2 - phi1, phi3 - phi2);
    let spread = values.iter().copied().fold(f64::NEG_INFINITY, f64::max)
        - values.iter().copied().fold(f64::INFINITY, f64::min);
    let scale = values.iter().fold(0.0f64, |m, v| m.max(v.abs()));
    let zero = |e: f64| e.abs() <= 1e-12 * scale;
    let none = |behaviour| Analysis {
        behaviour,
        observed_order: None,
        order_used: None,
        extrapolated: None,
        gci: None,
        spread,
    };
    if zero(e21) && zero(e32) {
        return Analysis {
            extrapolated: Some(phi1),
            gci: Some(0.0),
            ..none(Behaviour::Converged)
        };
    }
    if zero(e21) {
        return none(Behaviour::Indeterminate);
    }
    let ratio = e21 / e32;
    if zero(e32) || ratio >= 1.0 {
        return none(Behaviour::Divergent);
    }
    if ratio < 0.0 {
        return none(Behaviour::Oscillatory);
    }
    let Some(p) = observed_order(r21, r32, e21, e32) else {
        return none(Behaviour::Indeterminate);
    };
    let used = p.min(FORMAL_ORDER);
    let growth = r21.powf(used) - 1.0;
    Analysis {
        behaviour: Behaviour::Monotone,
        observed_order: Some(p),
        order_used: Some(used),
        extrapolated: Some(phi1 + (phi1 - phi2) / growth),
        gci: Some(SAFETY_FACTOR * e21.abs() / growth),
        spread,
    }
}

fn parse_splits(block: &J) -> Result<[usize; 3]> {
    let message = "grid_convergence.splits must be three increasing whole numbers in 1..=16 (cells split per axis on each level, e.g. [1, 2, 3])";
    let items = block
        .get("splits")
        .and_then(J::as_array)
        .filter(|items| items.len() == 3)
        .ok_or_else(|| bad(message))?;
    let mut splits = [0usize; 3];
    for (slot, item) in splits.iter_mut().zip(items) {
        match item.as_f64() {
            Some(v) if v >= 1.0 && v.fract() == 0.0 && v <= MAX_SPLIT as f64 => {
                *slot = v as usize;
            }
            _ => return Err(bad(message)),
        }
    }
    if !(splits[0] < splits[1] && splits[1] < splits[2]) {
        return Err(bad(message));
    }
    Ok(splits)
}

fn opt(value: Option<f64>) -> Result<String> {
    value.map_or_else(|| Ok("null".to_string()), num)
}

/// Run the grid-convergence study declared in `root`.
#[allow(clippy::too_many_lines)] // admission, three solves, one report
pub(super) fn run(root: &J, base: &std::path::Path, json_mode: bool) -> Result<String> {
    let block = root.get("grid_convergence").expect("checked by the caller");
    let splits = parse_splits(block)?;
    let mut template = root.clone();
    if let J::Object(entries) = &mut template {
        entries.retain(|(key, _)| key != "grid_convergence");
    }
    // Field files come from the finest level only (levels run in parallel).
    let mut coarse_template = template.clone();
    if let J::Object(entries) = &mut coarse_template {
        entries.retain(|(key, _)| key != "output");
    }
    // Every level is admitted before any is solved.
    let mut scenes = Vec::with_capacity(3);
    for (level, &split) in splits.iter().enumerate() {
        let source = if level == 2 {
            &template
        } else {
            &coarse_template
        };
        let scene = Scene::from_root_split(source, base, split).map_err(|failure| Failure {
            code: failure.code,
            message: format!("grid_convergence split {split}: {}", failure.message),
        })?;
        scenes.push(scene);
    }
    let cells: Vec<usize> = scenes
        .iter()
        .map(|scene| scene.grid.dims().iter().product())
        .collect();
    let workers = workers(block, "grid_convergence", 3)?;
    let started = Instant::now();
    let outcomes = parallel_map(&scenes, workers, |scene| {
        let begun = Instant::now();
        execute_within_budget(scene, true)
            .and_then(|text| {
                J::parse(text.trim()).map_err(|e| bad(format!("unreadable result: {e:?}")))
            })
            .map(|result| (result, begun.elapsed().as_secs_f64()))
    });
    let mut results = Vec::with_capacity(3);
    for (outcome, split) in outcomes.into_iter().zip(splits) {
        let (result, wall) = outcome.map_err(|failure| Failure {
            code: failure.code,
            message: format!("grid_convergence split {split}: {}", failure.message),
        })?;
        results.push((result, wall));
    }
    let wall_s = started.elapsed().as_secs_f64();
    let fine = &results[2].0;
    let rows: Vec<(String, [f64; 3], Analysis)> = quantity_names(fine)
        .into_iter()
        .filter_map(|name| {
            let values = [0, 1, 2].map(|k| quantity(&results[k].0, &name));
            let values = [values[0]?, values[1]?, values[2]?];
            values
                .iter()
                .all(|v| v.is_finite())
                .then(|| (name, values, analyse(values, splits)))
        })
        .collect();
    if json_mode {
        let mut out = format!(
            "{{\"schema\":{},\"status\":\"completed\",\"splits\":[{},{},{}],\"threads\":{workers},\"levels\":[",
            quote(CONVERGENCE_SCHEMA),
            splits[0],
            splits[1],
            splits[2]
        );
        for (k, ((_, wall), split)) in results.iter().zip(splits).enumerate() {
            if k > 0 {
                out.push(',');
            }
            let _ = write!(
                out,
                "{{\"split\":{split},\"cells\":{},\"wall_s\":{}}}",
                cells[k],
                num(*wall)?
            );
        }
        out.push_str("],\"quantities\":[");
        for (i, (name, values, analysis)) in rows.iter().enumerate() {
            if i > 0 {
                out.push(',');
            }
            let _ = write!(
                out,
                "{{\"name\":{},\"values\":[{},{},{}],\"convergence\":{},\"observed_order\":{},\"order_used\":{},\"extrapolated\":{},\"gci\":{},\"spread\":{}}}",
                quote(name),
                num(values[0])?,
                num(values[1])?,
                num(values[2])?,
                quote(analysis.behaviour.name()),
                opt(analysis.observed_order)?,
                opt(analysis.order_used)?,
                opt(analysis.extrapolated)?,
                opt(analysis.gci)?,
                num(analysis.spread)?
            );
        }
        let _ = writeln!(
            out,
            "],\"fine_result\":{},\"wall_s\":{},\"evidence\":\"Estimated\",\"no_claim\":{}}}",
            json_text(fine),
            num(wall_s)?,
            quote(&format!(
                "three-grid Celik/Roache GCI (safety factor 1.25, order capped at 2): an asymptotic band, not a bound, and only for monotone quantities; off-grid geometry restaircases between levels; each level: {NO_CLAIM}"
            ))
        );
        Ok(out)
    } else {
        let mut out = format!(
            "status=completed\nsplits={},{},{}\n",
            splits[0], splits[1], splits[2]
        );
        for (k, ((_, wall), split)) in results.iter().zip(splits).enumerate() {
            let _ = writeln!(
                out,
                "level split={split} cells={} wall_s={wall:.3}",
                cells[k]
            );
        }
        for (name, values, analysis) in &rows {
            let shown = |v: Option<f64>| v.map_or_else(|| "none".to_string(), |v| format!("{v}"));
            let _ = writeln!(
                out,
                "quantity={name} values={},{},{} convergence={} observed_order={} extrapolated={} gci={}",
                values[0],
                values[1],
                values[2],
                analysis.behaviour.name(),
                shown(analysis.observed_order),
                shown(analysis.extrapolated),
                shown(analysis.gci)
            );
        }
        let _ = writeln!(out, "wall_s={wall_s:.3}\nevidence=Estimated");
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn equal_ratios_recover_order_and_limit_exactly() {
        // phi(h) = 7 + 3 h^2 on h = 1, 1/2, 1/4.
        let phi = |h: f64| 3.0f64.mul_add(h * h, 7.0);
        let a = analyse([phi(1.0), phi(0.5), phi(0.25)], [1, 2, 4]);
        assert_eq!(a.behaviour, Behaviour::Monotone);
        assert!((a.observed_order.unwrap() - 2.0).abs() < 1e-12);
        assert!((a.extrapolated.unwrap() - 7.0).abs() < 1e-12);
        // GCI = 1.25 |phi1 - phi2| / 3 = 1.25 * 3 * (1/4 - 1/16) / 3.
        assert!((a.gci.unwrap() - 1.25 * 0.1875).abs() < 1e-12);
    }

    #[test]
    fn unequal_ratios_iterate_to_the_observed_order() {
        // phi(h) = 1 + 0.4 h^1.5 on splits 1, 2, 3 (ratios 2 and 1.5).
        let phi = |s: f64| 0.4f64.mul_add((1.0 / s).powf(1.5), 1.0);
        let a = analyse([phi(1.0), phi(2.0), phi(3.0)], [1, 2, 3]);
        assert_eq!(a.behaviour, Behaviour::Monotone);
        assert!((a.observed_order.unwrap() - 1.5).abs() < 1e-9, "{a:?}");
        assert!((a.extrapolated.unwrap() - 1.0).abs() < 1e-9, "{a:?}");
    }

    #[test]
    fn super_formal_orders_are_capped_and_others_get_no_band() {
        // Fourth-order data: the band uses order 2 (wider than order 4's).
        let phi = |s: f64| 1.0 + (1.0 / s).powi(4);
        let a = analyse([phi(1.0), phi(2.0), phi(4.0)], [1, 2, 4]);
        assert!((a.observed_order.unwrap() - 4.0).abs() < 1e-9);
        assert_eq!(a.order_used, Some(2.0));
        assert!(a.gci.unwrap() > (phi(2.0) - phi(4.0)) / 15.0);
        let oscillating = analyse([1.0, 0.9, 0.95], [1, 2, 3]);
        assert_eq!(oscillating.behaviour, Behaviour::Oscillatory);
        assert!(oscillating.gci.is_none());
        assert!((oscillating.spread - 0.1).abs() < 1e-15);
        let growing = analyse([1.0, 1.1, 1.3], [1, 2, 3]);
        assert_eq!(growing.behaviour, Behaviour::Divergent);
        assert!(growing.gci.is_none());
        assert_eq!(
            analyse([2.0, 2.0, 2.0], [1, 2, 3]).behaviour,
            Behaviour::Converged
        );
        assert_eq!(
            analyse([2.5, 2.0, 2.0], [1, 2, 3]).behaviour,
            Behaviour::Indeterminate
        );
    }
}
