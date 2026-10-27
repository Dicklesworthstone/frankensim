//! Explicit SI geometry and constitutive inputs. No manufacturer lookup or fitted pitch.
use super::Error;
use fs_couple::render::plate::impact::cavity::meridian::Station;
use fs_plate::shell::head::TensionedDiskSpec;
use std::{collections::BTreeMap, io::Read};

#[derive(Clone, Debug)]
pub struct Spec {
    pub profile: Vec<Station>,
    pub exterior: Vec<Station>,
    pub head: TensionedDiskSpec,
    pub head_band_hz: [f64; 2],
    pub head_damping_ratio: f64,
    pub maximum_slope: f64,
    pub air: [f64; 3], // rho [kg/m3], c [m/s], momentum drag [1/s]
    pub air_mesh: [usize; 2],
    pub air_maximum_hz: f64,
    pub rate_hz: u32,
    // Effective mass, radius, thickness, gap, speed, strike x, strike y.
    pub mallet: [f64; 7],
    // Reference stress/strain, loading/unloading powers, crush, densification, history.
    pub felt: [f64; 7],
    pub creep: Vec<[f64; 2]>,
}
fn scalars<const N: usize>(row: &[&str]) -> Result<[f64; N], Error> {
    if row.len() != N { return Err("wrong number of physical values".into()); }
    let mut out = [0.0_f64; N];
    for (out, text) in out.iter_mut().zip(row) {
        *out = text.parse::<f64>()?;
        if !out.is_finite() { return Err("all physical values must be finite".into()); }
    }
    Ok(out)
}
impl Spec {
    pub fn load(path: &str) -> Result<Self, Error> {
        let mut text = String::new();
        std::fs::File::open(path)?.take(65_537).read_to_string(&mut text)?;
        Self::read(&text)
    }
    pub fn read(text: &str) -> Result<Self, Error> {
        if text.len() > 65_536 { return Err("kettle specification exceeds 64 KiB".into()); }
        let mut lines = text.lines().map(|s| s.split('#').next().unwrap_or("").trim())
            .filter(|s| !s.is_empty());
        if lines.next() != Some("frankensim-kettle-v1") { return Err("expected frankensim-kettle-v1".into()); }
        let mut rows = BTreeMap::new();
        let mut profile = Vec::new();
        let mut exterior = Vec::new();
        let mut creep = Vec::new();
        for line in lines {
            let fields: Vec<_> = line.split(',').map(str::trim).collect();
            match fields[0] {
                "station" | "outer_station" => {
                    let [z_m, radius_m] = scalars(&fields[1..])?;
                    let target = if fields[0] == "station" { &mut profile } else { &mut exterior };
                    target.push(Station { z_m, radius_m });
                    if target.len() > 64 { return Err("at most 64 wall stations".into()); }
                }
                "creep" => {
                    creep.push(scalars(&fields[1..])?);
                    if creep.len() > 4 { return Err("at most four whole-face Kelvin branches".into()); }
                }
                "head" | "head_mesh" | "head_band_hz" | "maximum_slope" | "air"
                | "air_mesh" | "air_maximum_hz" | "rate_hz" | "mallet" | "felt" => {
                    if rows.insert(fields[0], fields[1..].to_vec()).is_some() {
                        return Err(format!("duplicate {} record", fields[0]).into());
                    }
                }
                _ => return Err(format!("unknown kettle record: {}", fields[0]).into()),
            }
        }
        let row = |key: &str| -> Result<&Vec<&str>, Error> {
            rows.get(key).ok_or_else(|| format!("missing {key} record").into())
        };
        let [thickness_m, young_pa, poisson, density_kg_m3, tension_n_m, head_damping_ratio] = scalars(row("head")?)?;
        let counts = |key: &str| -> Result<[usize; 2], Error> {
            let r = row(key)?;
            if r.len() != 2 { return Err("mesh requires two integer counts".into()); }
            Ok([r[0].parse()?, r[1].parse()?])
        };
        let [radial_intervals, azimuths] = counts("head_mesh")?;
        let rate = row("rate_hz")?;
        if rate.len() != 1 { return Err("rate_hz requires one integer".into()); }
        let result = Self {
            head: TensionedDiskSpec { radius_m: profile.last().ok_or("missing chamber profile")?.radius_m,
                thickness_m, young_pa, poisson, density_kg_m3, tension_n_m, radial_intervals, azimuths },
            profile, exterior, head_band_hz: scalars(row("head_band_hz")?)?, head_damping_ratio,
            maximum_slope: scalars::<1>(row("maximum_slope")?)?[0], air: scalars(row("air")?)?,
            air_mesh: counts("air_mesh")?, air_maximum_hz: scalars::<1>(row("air_maximum_hz")?)?[0],
            rate_hz: rate[0].parse()?, mallet: scalars(row("mallet")?)?, felt: scalars(row("felt")?)?, creep,
        };
        result.validate()?;
        Ok(result)
    }
    pub fn validate(&self) -> Result<(), Error> {
        let h = self.head;
        if self.profile.len() < 3 || self.profile.len() > 64
            || self.profile[0].radius_m != 0.0
            || self.profile.iter().any(|s| !s.z_m.is_finite() || !s.radius_m.is_finite() || s.radius_m < 0.0)
            || self.profile[1..].iter().any(|s| s.radius_m <= 0.0)
            || self.profile.windows(2).any(|s| s[1].z_m <= s[0].z_m)
            || h.radius_m != self.profile.last().ok_or("missing rim")?.radius_m
        { return Err("a sealed single-head chamber needs an ordered wall, axis bottom and positive-radius rim".into()); }
        if !h.radius_m.is_finite() || h.radius_m <= 0.0 || !h.tension_n_m.is_finite() || h.tension_n_m <= 0.0
            || !(1..=16).contains(&h.radial_intervals) || !(8..=96).contains(&h.azimuths)
            || !self.head_damping_ratio.is_finite() || self.head_damping_ratio < 0.0
            || !self.maximum_slope.is_finite() || !(0.0 < self.maximum_slope && self.maximum_slope < 1.0)
            || !(8_000..=1_536_000).contains(&self.rate_hz)
            || self.head_band_hz.iter().any(|f| !f.is_finite())
            || self.head_band_hz[0] < 0.0 || self.head_band_hz[1] <= self.head_band_hz[0]
            || self.head_band_hz[1] / f64::from(self.rate_hz) >= 0.45
            || self.air.iter().any(|v| !v.is_finite()) || self.air[0] <= 0.0 || self.air[1] <= 0.0
            || self.air[2] < 0.0 || self.air[2] > 100_000.0
            || !(2..=16).contains(&self.air_mesh[0]) || self.air_mesh[1] > 8
            || !self.air_maximum_hz.is_finite() || self.air_maximum_hz <= 0.0
            || self.air_maximum_hz / f64::from(self.rate_hz) >= 0.45
        { return Err("invalid material, clock, slope, mesh or spectral budget".into()); }
        fs_plate::PlateSection::isotropic(h.young_pa, h.poisson, h.thickness_m, h.density_kg_m3)?;
        if self.mallet.iter().chain(&self.felt).any(|v| !v.is_finite())
            || self.mallet[..3].iter().any(|v| *v <= 0.0) || self.mallet[3] < 0.0
            || !(0.0..=20.0).contains(&self.mallet[4])
            || self.felt[6] < 0.0 || self.felt[6] >= self.felt[5]
            || self.creep.len() > 4 || self.creep.iter().flatten().any(|v| !v.is_finite() || *v <= 0.0)
        { return Err("invalid physical mallet, felt history or Kelvin data".into()); }
        fs_material::fiber::WoolFelt::new(self.felt[0], self.felt[1], self.felt[2], self.felt[3], self.felt[4], self.felt[5])?;
        if !self.exterior.is_empty() { self.validate_exterior()?; }
        Ok(())
    }
    /// Separate physical outer wall; no automatic uniform-thickness offset.
    pub fn validate_exterior(&self) -> Result<(), Error> {
        let wall = &self.exterior;
        if wall.len() < 3 || wall.len() > 64 || self.profile.len() < 3
            || wall[0].radius_m != 0.0 || wall[0].z_m >= self.profile[0].z_m
            || wall.iter().any(|p| !p.z_m.is_finite() || !p.radius_m.is_finite() || p.radius_m < 0.0)
            || wall[1..].iter().any(|p| p.radius_m <= 0.0)
            || wall.windows(2).any(|p| p[1].z_m <= p[0].z_m)
            || wall.last().unwrap().z_m != self.profile.last().unwrap().z_m
            || wall.last().unwrap().radius_m <= self.head.radius_m
        { return Err("audio needs an explicit closed outer wall below the inner bottom and beyond the head rim, at the same top plane".into()); }
        // Difference of two affine radial graphs is affine on each common
        // interval: checking the union of their knots proves enclosure there.
        let radius = |profile: &[Station], z: f64| -> Option<f64> {
            profile.windows(2).find(|p| p[0].z_m <= z && z <= p[1].z_m).map(|p| {
                let t = (z-p[0].z_m)/(p[1].z_m-p[0].z_m);
                (1.0-t)*p[0].radius_m+t*p[1].radius_m
            })
        };
        for p in self.profile.iter().chain(wall) {
            if let Some(inner) = radius(&self.profile,p.z_m) {
                let outer = radius(wall,p.z_m).ok_or("outer wall does not cover inner chamber")?;
                if !outer.is_finite() || outer <= inner { return Err("outer wall intersects or enters the supplied inner chamber".into()); }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod scalar_tests {
    use super::{scalars, Spec};

    #[test]
    fn physical_scalars_keep_double_precision_and_signed_values() {
        let values = scalars::<4>(&["1e-300", "1e300", "-2.5", "-0.0"]).unwrap();
        assert_eq!(values, [1e-300_f64, 1e300_f64, -2.5, -0.0]);
        assert!(values[3].is_sign_negative());
        assert_eq!(scalars::<0>(&[]).unwrap(), [0.0_f64; 0]);
    }

    #[test]
    fn nonfinite_malformed_and_wrong_width_scalars_refuse() {
        for text in ["NaN", "inf", "-inf", "1e999", "", "not-a-number"] {
            assert!(scalars::<1>(&[text]).is_err(), "accepted {text}");
        }
        assert!(scalars::<2>(&["1"]).is_err());
        assert!(scalars::<1>(&["1", "2"]).is_err());
    }

    #[test]
    fn full_instrument_card_rejects_nonfinite_material_and_geometry() {
        let card = include_str!("illustrative-parabolic.fsk");
        let baseline = Spec::read(card).unwrap();
        assert_eq!(baseline.head.radius_m, 0.3302);
        for changed in [
            card.replace("air,1.2,343,0", "air,NaN,343,0"),
            card.replace("station,0,0", "station,inf,0"),
            card.replace("creep,3000,6", "creep,3000,1e999"),
        ] {
            assert_ne!(changed, card);
            assert!(Spec::read(&changed).is_err());
        }
    }
}
