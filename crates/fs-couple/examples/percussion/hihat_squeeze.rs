//! Physical thin-gap air loading; no oscillator, output filter, or second clock.
//! This adapter samples the two actual inner skins and compiles a polar pressure
//! graph for fs-flux. The shared impact owner supplies motion and pressure work.
use super::{Error, Shell, ShellFace, Spec};
use fs_couple::render::plate::impact::squeeze::{
    FilmCell, FilmChannel, FilmLimits, GapPort, ResistiveFilm,
};
use std::io::Read;
use std::path::Path;

const MAX_BYTES: usize = 8192;

pub(super) struct Config {
    inner_m: f64,
    outer_m: f64,
    radial: usize,
    azimuths: usize,
    viscosity_pa_s: f64,
    inner_open: bool,
    outer_open: bool,
    limits: FilmLimits,
}

pub(super) fn option(args: &mut Vec<String>) -> Result<Option<Config>, Error> {
    let positions: Vec<_> = args.iter().enumerate()
        .filter_map(|(i, s)| (s == "--squeeze-film").then_some(i)).collect();
    if positions.len() > 1 { return Err("--squeeze-film may be supplied only once".into()); }
    let Some(&i) = positions.first() else { return Ok(None); };
    let path = args.get(i + 1).filter(|s| !s.is_empty() && !s.starts_with("--"))
        .ok_or("--squeeze-film requires an input file; see SQUEEZE_FILM.md")?;
    let mut text = String::new();
    std::fs::File::open(Path::new(path))?.take((MAX_BYTES + 1) as u64).read_to_string(&mut text)?;
    let config = Config::parse(&text)?;
    args.drain(i..i + 2);
    Ok(Some(config))
}

impl Config {
    fn parse(text: &str) -> Result<Self, Error> {
        if text.len() > MAX_BYTES { return Err("squeeze-film input exceeds 8 KiB".into()); }
        let (mut header, mut annulus, mut viscosity, mut limits, mut boundary) =
            (false, None, None, None, None);
        for (line, raw) in text.lines().enumerate() {
            let row = raw.split('#').next().unwrap_or("").trim();
            if row.is_empty() { continue; }
            let bad = || -> Error { format!("squeeze-film line {}: invalid or duplicate record", line + 1).into() };
            if !header {
                if row != "frankensim-squeeze-film-v1" { return Err(bad()); }
                header = true;
                continue;
            }
            let f: Vec<_> = row.split(',').map(str::trim).collect();
            let number = |i: usize| -> Result<f64, Error> {
                let x: f64 = f[i].parse().map_err(|_| bad())?;
                if !x.is_finite() { return Err(bad()); }
                Ok(x)
            };
            match f[0] {
                "annulus" if f.len() == 5 && annulus.is_none() => {
                    let nr: usize = f[3].parse().map_err(|_| bad())?;
                    let nt: usize = f[4].parse().map_err(|_| bad())?;
                    annulus = Some((number(1)?, number(2)?, nr, nt));
                }
                "viscosity_pa_s" if f.len() == 2 && viscosity.is_none() => viscosity = Some(number(1)?),
                "limits" if f.len() == 4 && limits.is_none() => limits = Some(FilmLimits {
                    minimum_cell_gap_m: number(1)?, maximum_gap_m: number(2)?, maximum_pressure_pa: number(3)?,
                }),
                "boundary" if f.len() == 3 && boundary.is_none() => {
                    let open = |s: &str| match s { "open" => Ok(true), "sealed" => Ok(false), _ => Err(bad()) };
                    boundary = Some((open(f[1])?, open(f[2])?));
                }
                _ => return Err(bad()),
            }
        }
        let (inner_m, outer_m, radial, azimuths) = annulus.ok_or("missing film annulus")?;
        let (inner_open, outer_open) = boundary.ok_or("missing inner/outer film boundary conditions")?;
        let out = Self { inner_m, outer_m, radial, azimuths, inner_open, outer_open,
            viscosity_pa_s: viscosity.ok_or("missing film viscosity")?, limits: limits.ok_or("missing film limits")? };
        out.validate()?;
        Ok(out)
    }

    fn validate(&self) -> Result<(), Error> {
        let cells = self.radial.checked_mul(self.azimuths).ok_or("film grid count overflow")?;
        if !self.inner_m.is_finite() || !self.outer_m.is_finite() || self.inner_m <= 0.0
            || self.outer_m <= self.inner_m || self.radial == 0 || self.azimuths < 3 || cells > 64
            || !self.viscosity_pa_s.is_finite() || self.viscosity_pa_s <= 0.0
            || !self.limits.minimum_cell_gap_m.is_finite() || self.limits.minimum_cell_gap_m <= 0.0
            || !self.limits.maximum_gap_m.is_finite() || self.limits.maximum_gap_m <= self.limits.minimum_cell_gap_m
            || !self.limits.maximum_pressure_pa.is_finite() || self.limits.maximum_pressure_pa <= 0.0
            || !(self.inner_open || self.outer_open)
        { return Err("film needs an ordered annulus, 1..64 cells, >=3 azimuths, positive SI limits and an ambient drain".into()); }
        Ok(())
    }

    pub(super) fn compile(&self, spec: &Spec, upper: &Shell, lower: &Shell,
        upper_start: usize, lower_start: usize, total: usize) -> Result<ResistiveFilm, Error>
    {
        self.compile_with(total, |[x, y]| {
            let a = upper.port([x, y], ShellFace::Negative)?;
            // Proper lower-body rotation, exactly as in inter-cymbal contact.
            let b = lower.port([x, -y], ShellFace::Negative)?;
            let mut closure = vec![0.0; total];
            for (k, w) in a.weights.iter().enumerate() { closure[upper_start + k] = -w; }
            for (k, w) in b.weights.iter().enumerate() { closure[lower_start + k] = -w; }
            Ok(GapPort { reference_m: spec.separation + a.position_m[2] + b.position_m[2], closure })
        })
    }

    // This seam also checks the polar finite-volume geometry against an analytic
    // flat-annulus solution, independently of shell meshing/eigenmode selection.
    fn compile_with(&self, total: usize,
        mut sample: impl FnMut([f64; 2]) -> Result<GapPort, Error>) -> Result<ResistiveFilm, Error>
    {
        self.validate()?;
        let dr = (self.outer_m - self.inner_m) / self.radial as f64;
        let angle = 2.0 * std::f64::consts::PI / self.azimuths as f64;
        let point = |r: f64, theta: f64| [r * theta.cos(), r * theta.sin()];
        let mut cells = Vec::with_capacity(self.radial * self.azimuths);
        let mut channels = Vec::new();
        for ring in 0..self.radial {
            let a = self.inner_m + ring as f64 * dr;
            let b = self.inner_m + (ring + 1) as f64 * dr;
            let r = (a + b) * 0.5;
            for sector in 0..self.azimuths {
                let theta = (sector as f64 + 0.5) * angle;
                let from = ring * self.azimuths + sector;
                cells.push(FilmCell { area_m2: 0.5 * (b*b - a*a) * angle,
                    gap: sample(point(r, theta))? });
                // Every face aperture is sampled on the actual geometry. A rim
                // face can close while adjacent fluid VOLUME cells remain open.
                if ring + 1 < self.radial {
                    channels.push(FilmChannel { from, to: Some(from + self.azimuths),
                        width_m: b * angle, length_m: dr, gap: sample(point(b, theta))? });
                } else if self.outer_open {
                    channels.push(FilmChannel { from, to: None, width_m: b * angle,
                        length_m: dr * 0.5, gap: sample(point(b, theta))? });
                }
                if ring == 0 && self.inner_open {
                    channels.push(FilmChannel { from, to: None, width_m: a * angle,
                        length_m: dr * 0.5, gap: sample(point(a, theta))? });
                }
                channels.push(FilmChannel { from, to: Some(ring * self.azimuths + (sector + 1) % self.azimuths),
                    width_m: dr, length_m: r * angle, gap: sample(point(r, (sector + 1) as f64 * angle))? });
            }
        }
        Ok(ResistiveFilm::new(cells, channels, total, self.viscosity_pa_s, self.limits)?)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    const INPUT: &str = "frankensim-squeeze-film-v1\nannulus,0.02,0.1,16,4\nviscosity_pa_s,0.000018\nlimits,0.000001,0.002,1000\nboundary,open,open\n";
    fn flat(gap: f64) -> GapPort { GapPort { reference_m: gap, closure: vec![0.0,1.0,-1.0,0.0] } }

    #[test]
    fn strict_complete_input_and_unambiguous_cli() {
        assert!(Config::parse(INPUT).is_ok());
        for text in [INPUT.replace("16,4", "17,4"), INPUT.replace("16,4", "1,2"),
            INPUT.replace("0.000018", "NaN"), INPUT.replace("open,open", "sealed,sealed"),
            INPUT.replace("0.02,0.1", "0.1,0.02"), format!("{INPUT}viscosity_pa_s,0.000018"),
            INPUT.replace("boundary,open,open\n", ""), format!("{INPUT}unknown,1")]
        { assert!(Config::parse(&text).is_err(), "{text}"); }
        let mut args = vec!["hihat".into(), "model.fshh".into()];
        let before = args.clone();assert!(option(&mut args).unwrap().is_none());assert_eq!(args,before);
        for mut args in [vec!["--squeeze-film".into()],
            vec!["--squeeze-film".into(), "--analytic-newton".into()],
            vec!["--squeeze-film".into(), "x".into(), "--squeeze-film".into(), "y".into()]] {
            let before = args.clone();assert!(option(&mut args).is_err());assert_eq!(args,before);
        }
    }

    #[test]
    fn annulus_load_converges_and_unrelated_bodies_receive_no_air_force() {
        let mut config = Config::parse(INPUT).unwrap();
        let mu = config.viscosity_pa_s;let a = config.inner_m;let b = config.outer_m;
        let h = 0.001_f64;let v = 0.01;
        let exact = 1.5 * std::f64::consts::PI * mu * v / h.powi(3)
            * (b.powi(4) - a.powi(4) - (b*b - a*a).powi(2) / (b/a).ln());
        let mut previous = f64::INFINITY;
        for rings in [2,4,8,16] {
            config.radial = rings;
            let film = config.compile_with(4, |_| Ok(flat(h))).unwrap();
            let mut force = [0.0;4];let mut pressure = vec![0.0;film.cell_count()];
            film.evaluate_into(&[0.0;4], &[0.0,v,0.0,0.0], &mut force, &mut pressure).unwrap();
            let error = (force[1]/exact - 1.0).abs();assert!(error < previous);previous = error;
            assert_eq!(force[0],0.0);assert_eq!(force[3],0.0);assert_eq!(force[1],-force[2]);
            film.evaluate_into(&[0.0;4], &[7.0,v,v,-4.0], &mut force, &mut pressure).unwrap();
            assert_eq!(force,[0.0;4]);
        }
        assert!(previous < 0.01);
    }

    #[test]
    fn closing_outer_drain_changes_pressure_and_full_trapping_refuses() {
        let config = Config::parse(INPUT).unwrap();
        let open = config.compile_with(4, |_| Ok(flat(0.001))).unwrap();
        let closed = config.compile_with(4, |[x,y]| {
            Ok(flat(if x.hypot(y) > config.outer_m - 1e-10 { 0.0 } else { 0.001 }))
        }).unwrap();
        let mut f_open = [0.0;4];let mut f_closed = [0.0;4];let mut pressure = vec![0.0;open.cell_count()];
        let v = [0.0,0.0001,0.0,0.0];
        open.evaluate_into(&[0.0;4],&v,&mut f_open,&mut pressure).unwrap();
        closed.evaluate_into(&[0.0;4],&v,&mut f_closed,&mut pressure).unwrap();
        assert!(f_closed[1] > f_open[1]);
        let trapped = config.compile_with(4, |[x,y]| {
            let r=x.hypot(y);Ok(flat(if r>config.outer_m-1e-10 || r<config.inner_m+1e-10 {0.0}else{0.001}))
        }).unwrap();
        let before = f_closed;
        assert!(trapped.evaluate_into(&[0.0;4],&v,&mut f_closed,&mut pressure).is_err());
        assert_eq!(f_closed,before);
    }
}
