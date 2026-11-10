//! Material property rows read from the FrankenSim matdb seed packs
//! (`data/matdb/seed-v1/<pack>/properties.tsv`), embedded at compile time.
//!
//! The seed TSV is the provenance-carrying source the `fs-matdb` pack compiler
//! consumes. The watch only needs the `scalar` and `curve` rows, so this module
//! reads those two row kinds directly and converts their units to SI. Every
//! number the watch uses from here can be traced to the cited table row.

/// One embedded seed pack: its matdb id and the raw `properties.tsv` text.
#[derive(Clone, Copy, Debug)]
pub struct SeedPack {
    pub id: &'static str,
    text: &'static str,
}

macro_rules! seed {
    ($name:ident, $id:literal) => {
        pub const $name: SeedPack = SeedPack {
            id: $id,
            text: include_str!(concat!(
                "../../../../data/matdb/seed-v1/",
                $id,
                "/properties.tsv"
            )),
        };
    };
}

seed!(GOLD, "gold-pure-nbs-c447");
seed!(SILVER, "silver-pure-nbs-c447");
seed!(COPPER, "copper-pure-nbs-c447");
seed!(COPPER_FUSION, "copper-fusion-nasa-tp3287");
seed!(COPPER_LIQUID, "copper-liquid-nasa-tp3287");
seed!(BRASS_C36000, "brass-c36000-cda-reference");
seed!(BERYLLIUM_COPPER, "beryllium-copper-c17200-ngk-reference");
seed!(MUSIC_WIRE, "steel-music-wire-string-heptonstall-ligo-p1400229");
seed!(STEEL_1018, "steel-1018-onlinemetals-reference");
seed!(INVAR, "invar-fe36ni-nist-thermal-reference");
seed!(SAPPHIRE, "sapphire-al2o3-crystran-reference");
seed!(N42, "jinshan-n42-pristine-temperature");
seed!(PENNZANE, "pennzane-shf-x-2000-bearing-oil");

/// Convert a value in a matdb unit spelling to SI. Unknown spellings are
/// returned unchanged only when they are already SI.
fn to_si(value: f64, unit: &str) -> Option<f64> {
    Some(match unit {
        "MPa" => value * 1.0e6,
        "GPa" => value * 1.0e9,
        "kA/m" => value * 1.0e3,
        "kJ/m3" => value * 1.0e3,
        "g/mL" => value * 1.0e3,
        "mm^2/s" => value * 1.0e-6,
        "degC" => value + 273.15,
        "Pa" | "1" | "K" | "T" | "K^-1" | "K-1" | "kg/m3" | "kg/m^3" | "J/kg/K" | "W/m/K"
        | "J/kg" | "Ohm*m" => value,
        _ => return None,
    })
}

impl SeedPack {
    fn rows(self, kind: &'static str) -> impl Iterator<Item = Vec<&'static str>> {
        self.text
            .lines()
            .map(|line| line.split('\t').collect::<Vec<_>>())
            .filter(move |cols| cols.first() == Some(&kind))
    }

    /// SI value of the `scalar` row whose property key is `property`.
    #[must_use]
    pub fn scalar(self, property: &str) -> Option<f64> {
        self.rows("scalar")
            .find(|cols| cols.get(1) == Some(&property))
            .and_then(|cols| {
                let value: f64 = cols.get(4)?.trim().parse().ok()?;
                to_si(value, cols.get(5)?.trim())
            })
    }

    /// Like [`Self::scalar`] but panics with the pack id when the row is missing:
    /// the embedded packs are fixed at compile time, so a miss is a code bug.
    #[must_use]
    pub fn require(self, property: &str) -> f64 {
        self.scalar(property)
            .unwrap_or_else(|| panic!("matdb pack {} has no scalar row {property}", self.id))
    }

    /// Knots `(axis, value)` of the `curve` row `property`, value in SI.
    #[must_use]
    pub fn curve(self, property: &str) -> Option<Curve> {
        let cols = self
            .rows("curve")
            .find(|cols| cols.get(1) == Some(&property))?;
        let unit = cols.get(6)?.trim();
        let mut knots = Vec::new();
        for pair in cols.get(7)?.split(',') {
            let (x, y) = pair.split_once(':')?;
            let x: f64 = x.trim().parse().ok()?;
            let y: f64 = y.trim().parse().ok()?;
            knots.push((x, to_si(y, unit)?));
        }
        if knots.len() < 2 {
            return None;
        }
        knots.sort_by(|a, b| a.0.total_cmp(&b.0));
        Some(Curve { knots })
    }
}

/// Piecewise-linear curve through matdb knots, clamped at both ends.
#[derive(Clone, Debug)]
pub struct Curve {
    pub knots: Vec<(f64, f64)>,
}

impl Curve {
    #[must_use]
    pub fn eval(&self, x: f64) -> f64 {
        let k = &self.knots;
        if x <= k[0].0 {
            return k[0].1;
        }
        let last = k[k.len() - 1];
        if x >= last.0 {
            return last.1;
        }
        let i = k.partition_point(|p| p.0 <= x).max(1);
        let (x0, y0) = k[i - 1];
        let (x1, y1) = k[i];
        y0 + (y1 - y0) * (x - x0) / (x1 - x0)
    }

    /// Central-difference slope dy/dx at `x` over the local knot spacing.
    #[must_use]
    pub fn slope(&self, x: f64) -> f64 {
        let h = 2.0;
        (self.eval(x + h) - self.eval(x - h)) / (2.0 * h)
    }

    #[must_use]
    pub fn domain(&self) -> (f64, f64) {
        (self.knots[0].0, self.knots[self.knots.len() - 1].0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn embedded_packs_resolve_the_rows_the_watch_uses() {
        assert_eq!(GOLD.require("density"), 19300.0);
        assert!((GOLD.require("young-modulus") - 71.016e9).abs() < 1.0);
        assert_eq!(GOLD.require("thermal-conductivity"), 296.0);
        assert_eq!(SILVER.require("specific-heat"), 234.304);
        assert_eq!(COPPER.require("melting-point"), 1356.35);
        assert!(COPPER_FUSION.require("latent-heat-fusion") > 2.0e5);
        assert!((MUSIC_WIRE.require("young-modulus") - 212.0e9).abs() < 1.0);
        assert_eq!(MUSIC_WIRE.require("loss-tangent-465um"), 1.4e-4);
        assert!((PENNZANE.require("density-25degC") - 840.0).abs() < 1e-9);
        assert!((PENNZANE.require("kinematic-viscosity-40degC") - 107e-6).abs() < 1e-12);
        assert!((N42.require("remanence-25degC") - 1.275).abs() < 1e-12);
        let e = INVAR.curve("young-modulus").expect("invar E(T)");
        assert!(e.eval(293.15) > 1.3e11 && e.eval(293.15) < 1.7e11);
        let k = INVAR.curve("conductivity").expect("invar k(T)");
        assert!(k.eval(77.0) < k.eval(293.0));
    }
}
