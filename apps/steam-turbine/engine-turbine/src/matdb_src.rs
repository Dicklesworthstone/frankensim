//! Material property rows read from FrankenSim's `data/matdb/seed-v1` packs,
//! embedded at compile time so both native and `wasm32` targets share identical
//! sourced numbers with zero filesystem calls.

/// One embedded `seed-v1` pack: its canonical matdb identifier and raw TSV text.
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

// Rotor, shaft, casing, blading, seals, frame, bearing, lube oil, generator, fluids
seed!(AISI_4140, "aisi-4140-rc33");
seed!(STEEL_4130, "steel-4130-sheet-normalized-mil-hdbk-5j");
seed!(STEEL_1018, "steel-1018-onlinemetals-reference");
seed!(STEEL_A36, "steel-a36-onlinemetals-reference");
seed!(SS_14016, "stainless-14016-thyssenkrupp-reference");
seed!(SS_17_4PH, "stainless-17-4ph-h1025-bar-mil-hdbk-5j");
seed!(SS_316L, "stainless-14404-thyssenkrupp-reference");
seed!(INCONEL_718, "inconel-718-aged-special-metals-2007");
seed!(INCONEL_625, "inconel-625-special-metals-2013");
seed!(CAST_IRON, "gray-cast-iron-s2-s");
seed!(BRONZE_C93200, "bronze-c93200-cda-reference");
seed!(TIN_PURE, "tin-pure-nbs-c447");
seed!(LEAD_PURE, "lead-pure-nbs-c447");
seed!(SHELL_OIL_S2, "shell-heat-transfer-oil-s2-2011");
seed!(COPPER_IACS, "copper-annealed-iacs-nbs-hb100");
seed!(COPPER_PURE, "copper-pure-nbs-c447");
seed!(M19_STEEL, "nasa-tn-d-8184-m19-material-deck");
seed!(M19_STEINMETZ, "torrent-2018-m19-steinmetz-inputs");
seed!(WATER_LIQ, "water-liquid-iapws-sr6-08");
seed!(STEAM_TRANSPORT, "water-vapor-transport-iapws2008-2011");
seed!(DRY_AIR, "air-dry-ussa1976");

fn to_si(value: f64, unit: &str) -> Option<f64> {
    Some(match unit {
        "MPa" => value * 1.0e6,
        "GPa" => value * 1.0e9,
        "kPa" => value * 1.0e3,
        "kJ/kg" => value * 1.0e3,
        "kA/m" => value * 1.0e3,
        "kJ/m3" => value * 1.0e3,
        "g/mL" => value * 1.0e3,
        "mm^2/s" | "mm2/s" => value * 1.0e-6,
        "mm" => value * 1.0e-3,
        "um" | "µm" => value * 1.0e-6,
        "degC" => value + 273.15,
        "10^-6/K" | "um/(m*K)" => value * 1.0e-6,
        "Pa" | "1" | "K" | "T" | "K^-1" | "K-1" | "1/K" | "kg/m3" | "kg/m^3" | "J/kg/K"
        | "J/(kg*K)" | "W/m/K" | "W/(m*K)" | "J/kg" | "Ohm*m" | "Pa*s" | "m^2/s" | "m2/s"
        | "W/kg" | "Hz" | "m/s" | "A/m" | "m" | "N/m" | "kg/mol" | "J/mol/K" => value,
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

    /// Like [`Self::scalar`], panicking with the pack id if a compile-time row is absent.
    #[must_use]
    pub fn require(self, property: &str) -> f64 {
        self.scalar(property)
            .unwrap_or_else(|| panic!("matdb pack {} has no scalar row {property}", self.id))
    }

    /// Knots `(axis, value)` of the `curve` row `property`, value converted to SI.
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

    /// Like [`Self::curve`], panicking if absent.
    #[must_use]
    pub fn require_curve(self, property: &str) -> Curve {
        self.curve(property)
            .unwrap_or_else(|| panic!("matdb pack {} has no curve row {property}", self.id))
    }
}

/// Piecewise-linear curve through matdb knots, clamped at domain ends.
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

    #[must_use]
    pub fn domain(&self) -> (f64, f64) {
        (self.knots[0].0, self.knots[self.knots.len() - 1].0)
    }
}
