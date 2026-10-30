//! The bar card: one struck bar as data (geometry, wood, mallet, gas, receiver, render controls).
//!
//! Format `frankensim.bar-card.v1`: a header line, then tab-separated `key value...` rows. `#` starts a
//! comment. Every physical value is expected to carry its source in a `source` row of the card; the
//! parser does not interpret `source`/`note` rows. Paths in the card are relative to the repository root.
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

const HEADER: &str = "frankensim.bar-card.v1";

/// The bar's thickness along its length. The cut is on the lower face; the top face stays flat.
#[derive(Clone, Debug, PartialEq)]
pub enum Profile {
    Uniform,
    /// A straight section of `straight_length_m` centred on the bar, flanked by circular arcs of the
    /// grindstone `stone_radius_m` that rise from `middle_thickness_m` to the full thickness.
    Arc {
        stone_radius_m: f64,
        middle_thickness_m: f64,
        straight_length_m: f64,
    },
}

/// A Helmholtz resonator under the bar: a closed box of `volume_m3` with one opening of `opening_area_m2`
/// in a plate `neck_length_m` thick (the circular neck of the same area; both ends flanged), over which
/// the bar's lower face sits across the interval `wet_x_m` along its length, over its whole width.
#[derive(Clone, Debug, PartialEq)]
pub struct Resonator {
    pub volume_m3: f64,
    pub opening_area_m2: f64,
    pub neck_length_m: f64,
    pub wet_x_m: [f64; 2],
    /// The printed resonance of this resonator alone, for the report only.
    pub measured_hz: Option<f64>,
}

#[derive(Clone, Debug)]
pub struct Card {
    pub name: String,
    pub text: String,
    pub length_m: f64,
    pub width_m: f64,
    pub thickness_m: f64,
    pub profile: Profile,
    pub density_kg_m3: f64,
    /// (E_L, E_T, E_R) along the mesh axes (x = length, y = width, z = thickness).
    pub youngs_pa: [f64; 3],
    /// (nu_12, nu_13, nu_23) in the mesh axes.
    pub poisson: [f64; 3],
    /// (G_12, G_23, G_31) in the mesh axes.
    pub shear_pa: [f64; 3],
    /// Grid counts of rung 1 along (x, y, z); rung r uses r times these.
    pub mesh_counts: [usize; 3],
    pub mesh_rungs: Vec<usize>,
    pub band_hz: f64,
    pub windows: usize,
    pub strike_fraction: f64,
    pub receiver_m: [f64; 3],
    pub bem_stride: usize,
    pub gas_path: PathBuf,
    pub loss_path: Option<PathBuf>,
    pub mallet_contact_pack: String,
    pub mallet_contact_path: PathBuf,
    pub mallet_properties_path: PathBuf,
    pub mallet_speed_m_s: f64,
    pub duration_s: f64,
    pub substeps: u32,
    pub rate_hz: u32,
    pub full_scale_pa: f64,
    /// Printed modal frequencies of the same bar, for the report only (never used to build the bar).
    pub reference_vertical_hz: Vec<f64>,
    pub resonator: Option<Resonator>,
}

impl Profile {
    /// Thickness at distance `x` from the bar's end.
    pub fn thickness(&self, length: f64, full: f64, x: f64) -> f64 {
        match *self {
            Profile::Uniform => full,
            Profile::Arc {
                stone_radius_m: r,
                middle_thickness_m: h_m,
                straight_length_m: s,
            } => {
                let d = (x - 0.5 * length).abs() - 0.5 * s;
                if d <= 0.0 {
                    return h_m;
                }
                let depth = full - h_m;
                let reach = (2.0 * r * depth - depth * depth).max(0.0).sqrt();
                if d >= reach {
                    full
                } else {
                    (h_m + r - (r * r - d * d).sqrt()).min(full)
                }
            }
        }
    }
}

struct Rows(BTreeMap<String, Vec<String>>);

impl Rows {
    fn get(&self, key: &str) -> Result<&[String], String> {
        self.0
            .get(key)
            .map(Vec::as_slice)
            .ok_or_else(|| format!("bar card is missing the row `{key}`"))
    }

    fn number(&self, key: &str) -> Result<f64, String> {
        let v = self.get(key)?;
        if v.len() != 1 {
            return Err(format!("row `{key}` needs exactly one value"));
        }
        parse_f64(key, &v[0])
    }

    fn numbers<const N: usize>(&self, key: &str) -> Result<[f64; N], String> {
        let v = self.get(key)?;
        if v.len() != N {
            return Err(format!("row `{key}` needs exactly {N} values"));
        }
        let mut out = [0.0; N];
        for (slot, text) in out.iter_mut().zip(v) {
            *slot = parse_f64(key, text)?;
        }
        Ok(out)
    }

    fn count(&self, key: &str) -> Result<usize, String> {
        let x = self.number(key)?;
        if x.fract() != 0.0 || x < 1.0 {
            return Err(format!("row `{key}` must be a positive integer"));
        }
        Ok(x as usize)
    }

    fn text(&self, key: &str) -> Result<&str, String> {
        let v = self.get(key)?;
        if v.len() != 1 {
            return Err(format!("row `{key}` needs exactly one value"));
        }
        Ok(&v[0])
    }
}

fn parse_f64(key: &str, text: &str) -> Result<f64, String> {
    let x: f64 = text
        .parse()
        .map_err(|_| format!("row `{key}`: `{text}` is not a number"))?;
    if x.is_finite() {
        Ok(x)
    } else {
        Err(format!("row `{key}`: `{text}` is not finite"))
    }
}

/// One `scalar` row of a matdb properties file, converted to SI by its unit column. Only units whose factor
/// is exact are accepted; anything else is refused rather than guessed.
pub fn pack_scalar(path: &Path, property: &str) -> Result<f64, String> {
    let text = std::fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
    for line in text.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.first() != Some(&"scalar") || f.get(3) != Some(&property) {
            continue;
        }
        let value: f64 = f
            .get(4)
            .and_then(|v| v.parse().ok())
            .filter(|v: &f64| v.is_finite())
            .ok_or_else(|| format!("{}: `{property}` is not a finite number", path.display()))?;
        let factor = match f.get(5).copied() {
            Some("Pa" | "kg/m3" | "1") => 1.0,
            Some("MPa") => 1e6,
            Some("GPa") => 1e9,
            other => {
                return Err(format!(
                    "{}: `{property}` has unit {other:?}, which the bar card does not convert",
                    path.display()
                ));
            }
        };
        return Ok(value * factor);
    }
    Err(format!("{}: no scalar `{property}`", path.display()))
}

impl Card {
    /// Read and validate a card; `root` resolves the card's relative paths.
    pub fn load(path: &Path, root: &Path) -> Result<Card, String> {
        let text = std::fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
        Card::parse(&text, root)
    }

    pub fn parse(text: &str, root: &Path) -> Result<Card, String> {
        let mut lines = text
            .lines()
            .filter(|l| !l.trim().is_empty() && !l.starts_with('#'));
        if lines.next().map(str::trim) != Some(HEADER) {
            return Err(format!("bar card must start with `{HEADER}`"));
        }
        let mut rows = BTreeMap::new();
        for line in lines {
            let mut fields = line.split('\t').map(str::trim);
            let key = fields.next().unwrap_or_default().to_string();
            if rows
                .insert(key.clone(), fields.map(str::to_string).collect())
                .is_some()
            {
                return Err(format!("row `{key}` appears twice"));
            }
        }
        let rows = Rows(rows);
        let length_m = rows.number("length_m")?;
        let width_m = rows.number("width_m")?;
        let thickness_m = rows.number("thickness_m")?;
        let profile = match rows.get("profile")?.first().map(String::as_str) {
            Some("uniform") => Profile::Uniform,
            Some("arc") => {
                let p = rows.numbers::<3>("profile_arc")?;
                Profile::Arc {
                    stone_radius_m: p[0],
                    middle_thickness_m: p[1],
                    straight_length_m: p[2],
                }
            }
            other => return Err(format!("profile must be `uniform` or `arc`, got {other:?}")),
        };
        // The wood is three cited packs, never numbers typed into the card: the density's pack, the
        // longitudinal modulus' pack and the pack whose Handbook ratios stand in for the unmeasured ones.
        let density_kg_m3 = pack_scalar(&root.join(rows.text("wood_density")?), "density")?;
        let e_l = pack_scalar(
            &root.join(rows.text("wood_youngs_longitudinal")?),
            "young_modulus_longitudinal",
        )?;
        let ratios_pack = root.join(rows.text("wood_ratios")?);
        let ratio = |property: &str| pack_scalar(&ratios_pack, property);
        // Mesh axes: x = length (L), y = width (T), z = thickness (R), as in Rucz et al. Sec. III.A.
        let youngs_pa = [e_l, e_l * ratio("et_over_el")?, e_l * ratio("er_over_el")?];
        let poisson = [ratio("nu_lt")?, ratio("nu_lr")?, ratio("nu_tr")?];
        let shear_pa = [
            e_l * ratio("glt_over_el")?,
            e_l * ratio("grt_over_el")?,
            e_l * ratio("glr_over_el")?,
        ];
        let card = Card {
            name: rows.text("name")?.to_string(),
            text: text.to_string(),
            length_m,
            width_m,
            thickness_m,
            profile,
            density_kg_m3,
            youngs_pa,
            poisson,
            shear_pa,
            mesh_counts: {
                let c = rows.numbers::<3>("mesh_counts")?;
                [c[0] as usize, c[1] as usize, c[2] as usize]
            },
            mesh_rungs: rows
                .get("mesh_rungs")?
                .iter()
                .map(|r| r.parse::<usize>().map_err(|_| format!("bad rung `{r}`")))
                .collect::<Result<_, _>>()?,
            band_hz: rows.number("band_hz")?,
            windows: rows.count("windows")?,
            strike_fraction: rows.number("strike_fraction")?,
            receiver_m: rows.numbers::<3>("receiver_m")?,
            bem_stride: rows.count("bem_stride")?,
            gas_path: root.join(rows.text("gas")?),
            loss_path: match rows.text("loss")? {
                "none" => None,
                p => Some(root.join(p)),
            },
            mallet_contact_pack: rows.text("mallet_contact_pack")?.to_string(),
            mallet_contact_path: root.join(rows.text("mallet_contact")?),
            mallet_properties_path: root.join(rows.text("mallet_properties")?),
            mallet_speed_m_s: rows.number("mallet_speed_m_s")?,
            duration_s: rows.number("duration_s")?,
            substeps: rows.count("substeps")? as u32,
            rate_hz: rows.count("rate_hz")? as u32,
            full_scale_pa: rows.number("full_scale_pa")?,
            reference_vertical_hz: rows
                .get("reference_vertical_hz")
                .map(|v| v.iter().filter_map(|x| x.parse().ok()).collect())
                .unwrap_or_default(),
            resonator: match rows.get("resonator") {
                Err(_) => None,
                Ok(kind) if kind.first().map(String::as_str) == Some("helmholtz") => {
                    Some(Resonator {
                        volume_m3: rows.number("resonator_volume_m3")?,
                        opening_area_m2: rows.number("resonator_opening_area_m2")?,
                        neck_length_m: rows.number("resonator_neck_length_m")?,
                        wet_x_m: rows.numbers::<2>("resonator_wet_x_m")?,
                        measured_hz: rows.number("resonator_measured_hz").ok(),
                    })
                }
                Ok(_) => return Err("resonator must be `helmholtz`".into()),
            },
        };
        card.validate()?;
        Ok(card)
    }

    fn validate(&self) -> Result<(), String> {
        let positive = [
            ("length_m", self.length_m),
            ("width_m", self.width_m),
            ("thickness_m", self.thickness_m),
            ("density_kg_m3", self.density_kg_m3),
            ("band_hz", self.band_hz),
            ("mallet_speed_m_s", self.mallet_speed_m_s),
            ("duration_s", self.duration_s),
            ("full_scale_pa", self.full_scale_pa),
        ];
        for (name, v) in positive {
            if !(v > 0.0) {
                return Err(format!("`{name}` must be positive"));
            }
        }
        let [nx, ny, nz] = self.mesh_counts;
        if nx % 4 != 0 || ny % 2 != 0 || nz == 0 {
            return Err("mesh_counts need nx divisible by 4, ny even, nz >= 1".into());
        }
        if self.mesh_rungs.is_empty() || self.mesh_rungs.iter().any(|&r| r == 0) {
            return Err("mesh_rungs must be positive integers".into());
        }
        let finest = self.mesh_rungs.iter().copied().max().unwrap_or(1);
        if self
            .mesh_counts
            .iter()
            .any(|c| (c * finest) % self.bem_stride != 0)
        {
            return Err("bem_stride must divide the finest mesh counts".into());
        }
        if !(0.0..=1.0).contains(&self.strike_fraction) {
            return Err("strike_fraction must lie in [0, 1]".into());
        }
        if let Some(r) = &self.resonator {
            let [x0, x1] = r.wet_x_m;
            if !(r.volume_m3 > 0.0
                && r.opening_area_m2 > 0.0
                && r.neck_length_m > 0.0
                && 0.0 <= x0
                && x0 < x1
                && x1 <= self.length_m)
            {
                return Err("resonator needs positive volume, opening and neck, and a wet interval inside the bar".into());
            }
        }
        if let Profile::Arc {
            stone_radius_m,
            middle_thickness_m,
            straight_length_m,
        } = self.profile
        {
            let depth = self.thickness_m - middle_thickness_m;
            let reach = (2.0 * stone_radius_m * depth - depth * depth)
                .max(0.0)
                .sqrt();
            if !(stone_radius_m > 0.0
                && middle_thickness_m > 0.0
                && middle_thickness_m < self.thickness_m
                && straight_length_m >= 0.0
                && straight_length_m + 2.0 * reach <= self.length_m)
            {
                return Err("the cut does not fit the bar: check stone radius, middle thickness and straight length".into());
            }
        }
        Ok(())
    }
}
