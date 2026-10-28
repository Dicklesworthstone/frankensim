//! Declared two-head drum geometry/materials, lowered through existing owners.
//! No authored resonances, retuned samples, implicit material defaults or new FEM.
use std::io::Read;
use fs_plate::shell::head::TensionVariation;
use super::{Error, ModePair, TensionedDisk, TensionedDiskSpec, mesh_budget};
use super::ambient::Ambient;
use fs_material::gas_species::input::ThermochemicalSources;
use fs_matdb::NormalizedMaterialCardPack;
use fs_couple::material_source::MaterialRecord;
use std::{collections::BTreeMap, path::Path, sync::Arc};

#[path="drum_material.rs"]
pub mod material;
pub use material::HeadMaterial;

pub const HEADER: &str = "frankensim-drum-spec-v1";
const MAX_BYTES: usize = 65_536;
/// A file's head loss field: no citable film or edge loss, entered as zero drag.
pub const STATED_UNKNOWN_HEAD_LOSS: &str = "stated-unknown";
/// A file's head loss field: the `--head-relaxation` material law supplies it.
pub const RELAXATION_HEAD_LOSS: &str = "relaxation";
const UNSOURCED_HEAD_DAMPING: &str = "an authored head damping ratio has no source; declare stated-unknown (zero modal drag, with a receipt) or relaxation (with --head-relaxation)";

/// Where a head's loss comes from. A file states it with a token and never
/// carries an authored damping ratio.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum HeadLoss {
    /// No citable film or edge loss: zero modal drag, printed as a stated unknown.
    StatedUnknown,
    /// The head's loss is the `--head-relaxation` material law: zero modal drag.
    Relaxation,
    /// `Spec::reference()`'s authored ratio: a library fixture, which no file can declare.
    Fixture,
}
impl HeadLoss {
    /// The receipt a render prints for this head.
    pub fn receipt(self, head: &str, damping_ratio: f64) -> String {
        match self {
            Self::StatedUnknown => format!("head {head} film and edge loss: STATED UNKNOWN, entered as 0; its decays are upper bounds"),
            Self::Relaxation => format!("head {head} loss: the --head-relaxation material law; no modal drag"),
            Self::Fixture => format!("head {head} loss: library fixture damping ratio {damping_ratio} (Spec::reference(), not a file)"),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Head {
    pub thickness_m: f64,
    pub young_pa: f64,
    pub poisson: f64,
    pub density_kg_m3: f64,
    pub tension_n_m: f64,
    /// Viscous modal ratio: q'' + 2*zeta*omega*q' + omega^2*q = force.
    pub damping_ratio: f64,
}
#[derive(Clone, Debug, PartialEq)]
pub struct Spec {
    pub radius_m: f64,
    pub depth_m: f64,
    pub outer_radius_m: f64,
    /// Batter then resonant head; both use the same geometric mesh.
    pub heads: [Head; 2],
    /// Original source for each head, retained beside its actual coefficients.
    pub head_materials: [Option<HeadMaterial>; 2],
    /// Where each head's loss comes from; only a Fixture carries a nonzero ratio.
    pub head_loss: [HeadLoss; 2],
    /// Equilibrated installed tensor variations; zero preserves uniform heads.
    pub tension_variation: [TensionVariation; 2],
    pub radial_intervals: usize,
    pub azimuths: usize,
    pub band_hz: [f64; 2],
    /// The air inside the shell and around it: compressed by the heads and
    /// radiated into. Never a constant; see the `ambient` module.
    pub ambient: Ambient,
}
impl Spec {
    /// Original estimated 14 x 6.5 inch host, including its original loss law.
    pub fn reference() -> Self {
        let batter = Head { thickness_m: 0.000254, young_pa: 4e9, poisson: 0.38,
            density_kg_m3: 1390.0, tension_n_m: 3000.0, damping_ratio: 0.001 };
        Self { radius_m: 0.1778-0.0075, depth_m: 0.1651, outer_radius_m: 0.1778,
            heads: [batter, Head { thickness_m: 0.0000762, tension_n_m: 1500.0, ..batter }],
            head_materials: [None, None],
            head_loss: [HeadLoss::Fixture; 2],
            tension_variation: [TensionVariation::default(); 2],
            radial_intervals: 5, azimuths: 32, band_hz: [80.0, 500.0],
            ambient: Ambient::reference() }
    }
    pub fn validate(&self) -> Result<(), Error> {
        if [self.radius_m, self.depth_m, self.outer_radius_m].iter().any(|v| !v.is_finite() || *v<=0.0)
            || self.outer_radius_m<=self.radius_m || !self.volume_m3().is_finite() || self.volume_m3()<=0.0 {
            return Err("drum geometry needs positive finite clear radius/depth and an outer radius beyond the head".into());
        }
        if !(1..=32).contains(&self.radial_intervals) || !(8..=128).contains(&self.azimuths)
            || self.band_hz.iter().any(|v| !v.is_finite()) || self.band_hz[0]<0.0
            || self.band_hz[1]<=self.band_hz[0]
            || !(std::f64::consts::TAU*self.band_hz[1]).powi(2).is_finite() {
            return Err("drum needs 1..32 radial intervals, 8..128 azimuths and a finite increasing nonnegative frequency window".into());
        }
        for field in self.tension_variation { field.validate()?; }
        for (head, source) in self.heads.iter().zip(&self.head_materials) {
            if let Some(source) = source { source.validate(head)?; }
        }
        for (head,loss) in self.heads.iter().zip(self.head_loss) {
            if loss!=HeadLoss::Fixture && head.damping_ratio!=0.0 {
                return Err("a stated-unknown or relaxation head has zero modal drag".into());
            }
        }
        for head in self.heads {
            if [head.thickness_m,head.young_pa,head.poisson,head.density_kg_m3,
                head.tension_n_m,head.damping_ratio].iter().any(|v| !v.is_finite())
                || head.tension_n_m<=0.0 || head.damping_ratio<0.0 {
                return Err("head parameters must be finite with positive tension and passive damping".into());
            }
            // The existing section owner remains the authority on E, nu, h, rho.
            fs_plate::PlateSection::isotropic(head.young_pa,head.poisson,
                head.thickness_m,head.density_kg_m3)?;
        }
        Ok(())
    }
    /// A `relaxation` head needs the `--head-relaxation` law that supplies its loss.
    pub fn admit_head_loss(&self, relaxation_supplied: bool) -> Result<(), Error> {
        if !relaxation_supplied && self.head_loss.contains(&HeadLoss::Relaxation) {
            return Err("a head declaring loss relaxation needs --head-relaxation FILE, the law that supplies it".into());
        }
        Ok(())
    }
    pub fn volume_m3(&self) -> f64 {
        std::f64::consts::PI*self.radius_m*self.radius_m*self.depth_m
    }
    pub fn admit_clock(&self, dt_s: f64, audio: bool) -> Result<(), Error> {
        self.validate()?;
        if !dt_s.is_finite() || dt_s<=0.0 || self.band_hz[1]*dt_s>=0.45 {
            return Err("declared head window crosses the mechanical Nyquist guard".into());
        }
        // The existing fixed-receiver BEM bake has this declared band. Do not
        // silently extend its authority when admitting a different instrument.
        if audio && (self.band_hz[0]<40.0 || self.band_hz[1]>1640.0) {
            return Err("drum audio requires a head window inside the existing 40..1640 Hz acoustic bake; broader windows remain CSV-only".into());
        }
        Ok(())
    }
    pub fn admit_snare(&self, wires: super::snare::SnareSet) -> Result<(), Error> {
        wires.mode_count()?;
        if (0.5*wires.length_m).hypot(0.5*wires.width_m)>=self.radius_m {
            return Err("the unchanged snare span does not fit the supplied head; wires are not automatically shortened or retensioned".into());
        }
        Ok(())
    }
    pub fn head(&self, index: usize) -> Result<TensionedDisk, Error> {
        self.validate()?;
        let h = self.heads.get(index).ok_or("head index must be batter or resonant")?;
        Ok(TensionedDisk::new_with_tension_variation(TensionedDiskSpec { radius_m: self.radius_m,
            thickness_m: h.thickness_m, young_pa: h.young_pa, poisson: h.poisson,
            density_kg_m3: h.density_kg_m3, tension_n_m: h.tension_n_m,
            radial_intervals: self.radial_intervals, azimuths: self.azimuths },self.tension_variation[index],mesh_budget())?)
    }
    pub fn prepare(&self, dt_s: f64, audio: bool) -> Result<(Vec<TensionedDisk>,Vec<Vec<ModePair>>),Error> {
        self.admit_clock(dt_s,audio)?;
        let mut films = Vec::with_capacity(2);
        let mut modes = Vec::with_capacity(2);
        let pi = std::f64::consts::PI;
        for i in 0..2 {
            let film = self.head(i)?;
            if let Some(source) = &self.head_materials[i] {
                eprintln!("head {} material: {}", ["batter", "resonant"][i], source.receipt());
            }
            if !self.tension_variation[i].is_zero() {
                eprintln!("head {i}: base tension {} N/m, equilibrated tensor variation {:?}; fixed rim, no inferred lug forces or measured-map claim",
                    self.heads[i].tension_n_m,self.tension_variation[i]);
            }
            let pairs = fs_modal::slice_window(&film.model.k,&film.model.m,
                ((2.0*pi*self.band_hz[0]).powi(2),(2.0*pi*self.band_hz[1]).powi(2)),
                &fs_plate::SliceOptions::default())?.modes;
            if pairs.is_empty() { return Err(format!("head {i} frequency window is empty; no fallback modes are inserted").into()); }
            films.push(film); modes.push(pairs);
        }
        if modes.iter().map(Vec::len).sum::<usize>()>63 {
            return Err("two heads exceed the 63-coordinate host/radiation budget; reduce the explicit window, not the returned mode set".into());
        }
        Ok((films,modes))
    }
    pub fn load(path: &str) -> Result<Self, Error> {
        Self::load_with_gas_sources(path,None)
    }
    pub fn load_with_gas_sources(path: &str,sources: Option<&ThermochemicalSources<'_>>) -> Result<Self, Error> {
        let mut text = String::new();
        std::fs::File::open(path)?.take((MAX_BYTES+1) as u64).read_to_string(&mut text)?;
        let parent = Path::new(path).parent().unwrap_or_else(|| Path::new("."));
        Self::read_with_all_sources(&text,sources,super::material_packs::section_pack_loader(parent))
    }
    /// Geometry, both heads, mesh, band and the `ambient,<gas record>` row are
    /// mandatory: an imported drum declares its own air. Optional per-head
    /// tension_variation and tension_harmonic rows add a physical equilibrated tensor, not pitches.
    /// Each head states its loss as stated-unknown or relaxation; an authored damping ratio refuses.
    /// Order is free and comments start with #.
    pub fn read(text: &str) -> Result<Self, Error> {
        Self::read_with_gas_sources(text,None)
    }
    pub fn read_with_gas_sources(text: &str,sources: Option<&ThermochemicalSources<'_>>) -> Result<Self, Error> {
        Self::read_with_all_sources(text,sources,|name| Err(format!("head_material {name} requires an explicit material-card loader")))
    }
    /// Resolve all source records from already-loaded data. Only the native
    /// file loader opens paths; malformed or duplicate declarations still refuse.
    pub fn read_with_all_sources(text: &str,sources: Option<&ThermochemicalSources<'_>>,
        mut load: impl FnMut(&str) -> Result<NormalizedMaterialCardPack, String>) -> Result<Self, Error>
    {
        if text.len()>MAX_BYTES { return Err("drum specification exceeds 64 KiB".into()); }
        let (mut header,mut geometry,mut mesh,mut band,mut ambient) = (false,None,None,None,None);
        let mut heads = [None,None];
        let mut head_materials = [None,None];
        let mut packs: BTreeMap<String, Arc<NormalizedMaterialCardPack>> = BTreeMap::new();
        let mut head_loss = [HeadLoss::StatedUnknown; 2];
        let mut variations = [None,None];
        let mut variation_seen = [false; 2];
        let mut harmonic_seen = [[false; fs_plate::shell::head::tension::TENSION_HARMONIC_COUNT];2];
        for (number,raw) in text.lines().enumerate() {
            let row = raw.split('#').next().unwrap_or("").trim();
            if row.is_empty() { continue; }
            let error = |why: &str| -> Error { format!("drum specification line {}: {why}",number+1).into() };
            if !header {
                if row!=HEADER { return Err(error("expected frankensim-drum-spec-v1")); }
                header=true; continue;
            }
            let fields: Vec<_> = row.split(',').map(str::trim).collect();
            let count = match fields[0] { "geometry"=>4, "head"=>8, "head_material"=>6, "mesh"|"band_hz"=>3, "tension_variation"=>9, "tension_harmonic"=>5,
                "ambient"=>2,
                _=>return Err(error("unknown record")) };
            if fields.len()!=count { return Err(error("wrong number of fields")); }
            let value = |i: usize| -> Result<f64,Error> {
                let v: f64 = fields[i].parse().map_err(|_| error("invalid scalar"))?;
                if !v.is_finite() { return Err(error("nonfinite scalar")); } Ok(v)
            };
            match fields[0] {
                "geometry" => {
                    if geometry.is_some() { return Err(error("duplicate geometry")); }
                    geometry=Some([value(1)?,value(2)?,value(3)?]);
                }
                "head"|"head_material" => {
                    let i=match fields[1] {"batter"=>0,"resonant"=>1,_=>return Err(error("unknown head"))};
                    if heads[i].is_some() { return Err(error("duplicate head")); }
                    let source_backed = fields[0] == "head_material";
                    head_loss[i]=match fields[if source_backed {4} else {7}] {
                        STATED_UNKNOWN_HEAD_LOSS=>HeadLoss::StatedUnknown,
                        RELAXATION_HEAD_LOSS=>HeadLoss::Relaxation,
                        other if other.parse::<f64>().is_ok() => return Err(error(UNSOURCED_HEAD_DAMPING)),
                        _=>return Err(error("unknown head loss; declare stated-unknown or relaxation")),
                    };
                    if source_backed {
                        let record = MaterialRecord::parse(&fields[5].split_ascii_whitespace().collect::<Vec<_>>())
                            .map_err(|e| error(&e.to_string()))?;
                        let pack = match packs.get(&record.pack) {
                            Some(pack) => pack.clone(),
                            None => {
                                let pack = Arc::new(load(&record.pack).map_err(|e| error(&e))?);
                                packs.insert(record.pack.clone(), pack.clone());
                                pack
                            }
                        };
                        let (head, source) = HeadMaterial::resolve(record, pack, value(2)?, value(3)?)
                            .map_err(|e| error(&e.to_string()))?;
                        heads[i] = Some(head); head_materials[i] = Some(source);
                    } else {
                        heads[i]=Some(Head {thickness_m:value(2)?,young_pa:value(3)?,poisson:value(4)?,
                            density_kg_m3:value(5)?,tension_n_m:value(6)?,damping_ratio:0.0});
                    }
                }
                "tension_variation" => {
                    let i=match fields[1] {"batter"=>0,"resonant"=>1,_=>return Err(error("unknown tension-variation head"))};
                    if variation_seen[i] { return Err(error("duplicate head tension variation")); }
                    variation_seen[i] = true;
                    // Affine and harmonic records are independent. Preserve any
                    // harmonics already supplied, including explicit zero rows.
                    let variation = variations[i].get_or_insert_with(TensionVariation::default);
                    variation.constant_n_m = [value(2)?, value(3)?, value(4)?];
                    variation.gradient_n_m2 = [value(5)?, value(6)?, value(7)?, value(8)?];
                }
                "tension_harmonic" => {
                    let i=match fields[1] {"batter"=>0,"resonant"=>1,_=>return Err(error("unknown tension-harmonic head"))};
                    let order=fields[2].parse::<usize>().map_err(|_|error("invalid tension-harmonic order"))?;
                    if !(fs_plate::shell::head::tension::MIN_TENSION_HARMONIC_ORDER
                        ..=fs_plate::shell::head::tension::MAX_TENSION_HARMONIC_ORDER).contains(&order) {
                        return Err(error("tension-harmonic order must be 2..12"));
                    }
                    let slot=order-fs_plate::shell::head::tension::MIN_TENSION_HARMONIC_ORDER;
                    if harmonic_seen[i][slot] {return Err(error("duplicate head tension harmonic"));}
                    harmonic_seen[i][slot]=true;
                    variations[i].get_or_insert_with(TensionVariation::default).harmonics_n_m[slot]=[value(3)?,value(4)?];
                }
                "mesh" => {
                    if mesh.is_some() { return Err(error("duplicate mesh")); }
                    mesh=Some([fields[1].parse::<usize>().map_err(|_| error("invalid radial count"))?,
                        fields[2].parse::<usize>().map_err(|_| error("invalid azimuth count"))?]);
                }
                "band_hz" => {
                    if band.is_some() { return Err(error("duplicate frequency window")); }
                    band=Some([value(1)?,value(2)?]);
                }
                "ambient" => {
                    if ambient.is_some() { return Err(error("duplicate ambient")); }
                    ambient=Some(Ambient::described_with_gas_sources(fields[1],sources).map_err(|e| error(&e.to_string()))?);
                }
                _ => unreachable!(),
            }
        }
        if !header { return Err("missing drum specification header".into()); }
        let [radius_m,depth_m,outer_radius_m]=geometry.ok_or("missing drum geometry")?;
        let [radial_intervals,azimuths]=mesh.ok_or("missing drum mesh")?;
        let result=Self {radius_m,depth_m,outer_radius_m,radial_intervals,azimuths,
            band_hz:band.ok_or("missing head frequency window")?,
            tension_variation:variations.map(Option::unwrap_or_default),
            heads:[heads[0].ok_or("missing batter head")?,heads[1].ok_or("missing resonant head")?],
            head_loss,head_materials,
            ambient:ambient.ok_or("missing drum ambient row: an imported drum must declare its air")?};
        result.validate()?; Ok(result)
    }
}

/// Remove only this option, leaving the shared playing/numerics controls intact.
pub fn option(args: &mut Vec<String>) -> Result<Option<String>,Error> {
    let mut positions=args.iter().enumerate().filter(|(_,s)|s.as_str()=="--drum-spec");
    let Some((start,_))=positions.next() else {return Ok(None);};
    if positions.next().is_some() {return Err("--drum-spec may be supplied only once".into());}
    let path=args.get(start+1).filter(|s|!s.is_empty() && !s.starts_with("--"))
        .ok_or("--drum-spec requires a file path")?.clone();
    args.drain(start..start+2); Ok(Some(path))
}
pub fn admit_command(path: Option<&str>, command: &str) -> Result<(),Error> {
    let base=command.strip_suffix("-wav").or_else(||command.strip_suffix("-mic")).unwrap_or(command);
    if path.is_some() && !matches!(base,"drum"|"drum-stretch"|"drum-modal"|"snare"|"snare-off") {
        return Err("--drum-spec applies to two-head drum/snare commands, not shell-profile cymbals".into());
    }
    Ok(())
}

#[cfg(test)]
#[path = "drum_spec_tests.rs"]
mod tests;

#[cfg(test)]
#[path="drum_material_tests.rs"]
mod material_tests;

#[cfg(test)]
#[path="tension_render_tests.rs"]
mod tension_tests;
