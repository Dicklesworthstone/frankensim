//! Profile-defined enclosed air and rigid exterior, composed with the original
//! head, stick, contact, material-memory, time and pressure owners.
use super::{cavity, drum_spec, Error, ImpactBody, Mechanics, ModePair, Obstacle, TensionedDisk};
use fs_couple::render::plate::impact::{MAX_IMPACT_MODES, cavity::CavityCoupling,
    cavity::meridian::{Budget, MeridianCavity, Station, end_disk::{DiskPort,EndDisk,InterfaceBudget}}, damping::ViscousDamper, felt::FeltPad};
use fs_couple::vibroacoustic::AcousticMedium;
use fs_exec::CancelGate;
use std::io::Read;

const MAX_BYTES: usize = 65_536;
/// Axial stations run DOWN from the batter head. The pressure and outer-skin
/// radii are independent supplied measurements; no uniform thickness is inferred.
#[derive(Clone, Debug)]
pub struct Spec {
    pub stations: Vec<Station>,
    pub outer_radii_m: Vec<f64>,
    pub budget: Budget,
}
impl Spec {
    pub fn load(path: &str) -> Result<Self, Error> {
        let mut text = String::new();
        std::fs::File::open(path)?.take((MAX_BYTES + 1) as u64).read_to_string(&mut text)?;
        Self::read(&text)
    }
    pub fn read(text: &str) -> Result<Self, Error> {
        if text.len() > MAX_BYTES { return Err("chamber profile exceeds 64 KiB".into()); }
        let (mut header, mut mesh, mut band, mut limits) = (false, None, None, None);
        let mut stations = Vec::new(); let mut outer = Vec::new();
        for (line, raw) in text.lines().enumerate() {
            let row = raw.split('#').next().unwrap_or("").trim();
            if row.is_empty() { continue; }
            let error = |why: &str| -> Error { format!("chamber line {}: {why}", line + 1).into() };
            if !header {
                if row != "frankensim-chamber-profile-v1" { return Err(error("wrong profile header")); }
                header = true; continue;
            }
            let f: Vec<_> = row.split(',').map(str::trim).collect();
            let number = |i: usize| -> Result<f64, Error> {
                let v: f64 = f[i].parse().map_err(|_| error("invalid scalar"))?;
                if !v.is_finite() { return Err(error("nonfinite scalar")); } Ok(v)
            };
            let integer = |i: usize| -> Result<usize, Error> {
                f[i].parse().map_err(|_| error("invalid count"))
            };
            match f[0] {
                "station" if f.len() == 4 => {
                    if stations.len() == 128 { return Err(error("at most 128 stations")); }
                    stations.push(Station { z_m: number(1)?, radius_m: number(2)? });
                    outer.push(number(3)?);
                }
                "mesh" if f.len() == 3 && mesh.is_none() => { mesh = Some([integer(1)?, integer(2)?]); }
                "band_hz" if f.len() == 2 && band.is_none() => { band = Some(number(1)?); }
                "budget" if f.len() == 5 && limits.is_none() => {
                    limits = Some((integer(1)?, integer(2)?, integer(3)?, number(4)?));
                }
                _ => return Err(error("unknown, duplicate or malformed record")),
            }
        }
        let mesh = mesh.ok_or("chamber needs mesh,radial_intervals,maximum_azimuthal_order")?;
        let (nodes, entries, modes, tolerance) = limits.ok_or("chamber needs budget,nodes,matrix_entries,modes,residual_tolerance")?;
        let result = Self { stations, outer_radii_m: outer, budget: Budget {
            radial_intervals: mesh[0], maximum_azimuthal_order: mesh[1],
            maximum_frequency_hz: band.ok_or("chamber needs band_hz,maximum")?,
            maximum_nodes: nodes, maximum_matrix_entries: entries, maximum_modes: modes,
            relative_tolerance: tolerance,
        }};
        if !header { return Err("missing chamber header".into()); }
        result.validate()?; Ok(result)
    }
    pub(super) fn validate(&self) -> Result<(), Error> {
        let b = self.budget;
        if !(2..=128).contains(&self.stations.len()) || self.outer_radii_m.len() != self.stations.len()
            || !(2..=32).contains(&b.radial_intervals) || b.maximum_azimuthal_order > 8
            || !(1..=512).contains(&b.maximum_nodes) || !(1..=8).contains(&b.maximum_modes)
            || b.maximum_matrix_entries == 0 || b.maximum_matrix_entries > 2_000_000
            || !b.maximum_frequency_hz.is_finite() || b.maximum_frequency_hz <= 0.0
            || !b.relative_tolerance.is_finite() || b.relative_tolerance <= 0.0 || b.relative_tolerance > 1e-3
            || self.stations[0].z_m != 0.0
            || self.stations.iter().zip(&self.outer_radii_m).any(|(s,r)|
                !s.z_m.is_finite() || !s.radius_m.is_finite() || s.radius_m <= 0.0
                || !r.is_finite() || *r <= s.radius_m)
            || self.stations.windows(2).any(|s| s[1].z_m <= s[0].z_m) {
            return Err("chamber needs ordered SI stations, positive inner/outer radii and bounded acoustic controls".into());
        }
        Ok(())
    }
    pub fn admit(&self, drum: &drum_spec::Spec, dt_s: f64, neck: bool, barrel: bool) -> Result<(), Error> {
        self.validate()?; drum.admit_clock(dt_s, false)?;
        let last = self.stations.len() - 1;
        // Exact common seams: neither air nor BEM snaps a measured lip to a
        // different head. Input decimal roundoff has no hidden geometry repair.
        if self.stations[last].z_m != drum.depth_m
            || [0,last].into_iter().any(|i| self.stations[i].radius_m != drum.radius_m
                || self.outer_radii_m[i] != drum.outer_radius_m) {
            return Err("chamber end radii/depth must equal the supplied head and bearing-edge geometry".into());
        }
        if self.budget.maximum_frequency_hz * dt_s >= 0.45 {
            return Err("chamber pressure window crosses the mechanical Nyquist guard".into());
        }
        if neck || barrel {
            return Err("profiled chamber with a sidewall neck or elastic barrel needs matching noncylindrical wall coupling; cylinder adapters are refused".into());
        }
        Ok(())
    }
    pub fn volume_m3(&self) -> f64 {
        self.stations.windows(2).map(|s| {
            let a=s[0].radius_m; let b=s[1].radius_m;
            core::f64::consts::PI/3.0*(s[1].z_m-s[0].z_m)*(a*a+a*b+b*b)
        }).sum()
    }
    fn air(&self, gate: &CancelGate) -> Result<MeridianCavity, Error> {
        Ok(MeridianCavity::new(self.stations.clone(), AcousticMedium {rho0:1.2,c0:343.0}, self.budget, gate)?)
    }
    /// Compile the whole original solid layout before appending any gas inertia.
    /// Only actual head surfaces overlap pressure: sticks/mutes/shafts are zeros.
    #[allow(clippy::too_many_arguments)]
    pub fn build(&self, films: &[TensionedDisk], modes: &[Vec<ModePair>], bodies: Vec<ImpactBody>,
        contacts: Vec<Obstacle>, pads: Vec<FeltPad>, dampers: Vec<ViscousDamper>,
        steps: u64, dt_s: f64, prepared: bool, snares: bool, drag: f64)
        -> Result<(Mechanics,cavity::InteriorPressure), Error> {
        cavity::validate_drag(drag)?;
        let gate = CancelGate::new_clock_free(); let air = self.air(&gate)?;
        let depth = self.stations.last().ok_or("empty chamber")?.z_m;
        let structural = bodies.iter().try_fold(0usize,|n,b|n.checked_add(b.initial.len()))
            .ok_or("chamber structural address overflow")?;
        let mut configuration = super::mechanics::coupled_config(steps,dt_s,snares)?;
        if drag > 0.0 { configuration.coupling.max_connections=32; configuration.coupling.max_setup_terms=1_000_000; }
        let maximum = if prepared {configuration.coupling.max_modes} else {MAX_IMPACT_MODES};
        if structural > maximum {return Err("chamber structural mode budget exceeded".into());}
        let sampled = air.sample(&[[0.0,0.0,0.0]],8,&gate)?;
        let overlap=head_overlaps(&air,films,modes,structural,&gate)?;
        let damping:Vec<_>=sampled.omegas.iter().map(|w|if *w==0.0 {0.0}else{drag}).collect();
        let coupling=CavityCoupling::new_with_mode_budget(&sampled,structural,&overlap,&damping,maximum)?;
        let first=air.values_at([0.0,0.0,0.25*depth])?;
        let second=air.values_at([0.0,0.0,0.75*depth])?;
        eprintln!("profiled chamber: volume_m3={}, pressure_hz={:?}; supplied piecewise-linear inner/outer walls; rigid sidewall, no equivalent cylinder, no calibration or real-time claim",
            air.volume_m3(),sampled.omegas.iter().map(|w|w/core::f64::consts::TAU).collect::<Vec<_>>());
        let (system,coupling)=if prepared {
            if !pads.is_empty() {return Err("prepared linear chamber cannot discard material-history pads".into());}
            let radius=films.first().ok_or("chamber needs a head")?.spec.radius_m;
            let (s,c)=coupling.build_linear_with_dampers(bodies,contacts,dampers,
                core::f64::consts::PI*radius*radius,configuration,&gate)?;
            (Mechanics::Prepared(s),c)
        } else {
            let (s,c)=coupling.build_with_dampers(bodies,contacts,pads,dampers,super::config(steps,dt_s),&gate)?;
            (Mechanics::Reference(s),c)
        };
        Ok((system,cavity::InteriorPressure {coupling,first,second}))
    }
}
/// Compose the shared end-disk interface; do not copy its surface integration.
fn head_overlaps(air:&MeridianCavity,films:&[TensionedDisk],modes:&[Vec<ModePair>],
    structural:usize,gate:&CancelGate)->Result<Vec<f64>,Error> {
    if films.len()!=2 || modes.len()!=2 || structural==0 || structural>4096 {
        return Err("chamber needs the two original head pencils and a bounded solid layout".into());
    }
    let mut result=vec![0.0;structural*air.modes().len()];let mut first=1;
    for (head,(film,pairs)) in films.iter().zip(modes).enumerate() {
        let port=DiskPort {end:if head==0 {EndDisk::Lower}else{EndDisk::Upper},
            first_mode:first,structural_modes:structural,positive_z:1.0};
        let interface=air.end_disk_interface(film,pairs,port,
            InterfaceBudget {maximum_points:12000,maximum_terms:1_000_000},gate)?;
        for (sum,value) in result.iter_mut().zip(interface.overlaps) {*sum+=value;}
        first+=pairs.len();
    }
    Ok(result)
}
pub fn option(args: &mut Vec<String>) -> Result<Option<String>, Error> {
    let mut matches=args.iter().enumerate().filter(|(_,a)|a.as_str()=="--chamber-profile");
    let Some((at,_))=matches.next() else {return Ok(None)};
    if matches.next().is_some() {return Err("--chamber-profile may appear only once".into());}
    let path=args.get(at+1).filter(|v|!v.is_empty() && !v.starts_with("--"))
        .ok_or("--chamber-profile needs an input path")?.clone();
    args.drain(at..at+2); Ok(Some(path))
}

#[cfg(test)]
mod tests;
