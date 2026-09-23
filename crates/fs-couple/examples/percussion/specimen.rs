//! Supplied physical shell geometry. Parsing is cold; all mesh, section, modal
//! and nonlinear mechanics stay with their existing owners. A file is input,
//! not proof that its numbers were measured or a calibrated manufacturer's CAD.
use super::Error;
use fs_plate::shell::profile::{AnnularRelief, ProfileStation,
    SurfaceIndentation, revolve};
use fs_plate::shell::survey::MeshShell;
use std::io::{Read, Write as _};
use std::fmt::Write as _;

#[path = "shell_mesh.rs"]
mod mesh_input;
use std::path::{Path, PathBuf};

const MAX_BYTES: usize = 1_048_576;

/// One explicit material, meridian and optional stress-free reference details.
/// Azimuth count and frequency window are numerical choices, not timbre controls.
pub struct Specimen {
    pub stations: Vec<ProfileStation>,
    pub rings: Vec<AnnularRelief>,
    pub dents: Vec<SurfaceIndentation>,
    pub azimuths: usize,
    pub young_pa: f64,
    pub poisson: f64,
    pub density_kg_m3: f64,
    pub band_hz: [f64; 2],
    // Profile fields above are used only by the original revolved input.
    sampled: Option<MeshShell>,
    default_strike: Option<[f64;2]>,
}
impl Specimen {
    /// Original estimated splash, unchanged when no external geometry is selected.
    pub fn reference() -> Self {
        Self {
            stations: [(0.00615,0.022,0.0016),(0.012,0.0215,0.0016),(0.020,0.020,0.0015),
                (0.039,0.008,0.0011),(0.050,0.0055,0.0009),(0.070,0.003,0.00075),
                (0.085,0.0015,0.0006),(0.1016,0.0,0.0005)].map(|(radius_m,height_m,thickness_m)|
                    ProfileStation {radius_m,height_m,thickness_m}).to_vec(),
            rings: vec![], dents: vec![], azimuths: 32,
            young_pa: 112.6e9, poisson: 0.342, density_kg_m3: 8607.0,
            band_hz: [50.0,1200.0], sampled: None, default_strike: None,
        }
    }

    pub fn load(path: &Path) -> Result<Self, Error> {
        // A changing file or a pipe cannot exceed the same in-memory input limit.
        let mut text = String::new();
        std::fs::File::open(path)?.take((MAX_BYTES+1) as u64).read_to_string(&mut text)?;
        Self::parse(&text)
    }

    /// Strict, line-oriented numeric records with SI units and no hidden material
    /// defaults. Comments begin with '#'. Unknown fields and duplicate singleton
    /// records refuse, so a misspelled thickness/material cannot go unused.
    pub fn parse(text: &str) -> Result<Self, Error> {
        if text.len()>MAX_BYTES { return Err("shell profile exceeds the 1 MiB input limit".into()); }
        let mut material=None; let mut azimuths=None; let mut band=None;
        let mut stations=Vec::new(); let mut rings=Vec::new(); let mut dents=Vec::new();
        let mut header=false;
        for (line,raw) in text.lines().enumerate() {
            let data=raw.split('#').next().unwrap_or("").trim();
            if data.is_empty() { continue; }
            if !header {
                if data!="frankensim-shell-profile-v1" { return Err("shell profile needs frankensim-shell-profile-v1 header".into()); }
                header=true; continue;
            }
            let fields:Vec<_>=data.split(',').map(str::trim).collect();
            let invalid=||format!("shell profile line {}: unknown, duplicated or malformed record",line+1);
            let numbers=|count:usize|->Result<Vec<f64>,Error> {
                if fields.len()!=count+1 {return Err(invalid().into());}
                fields[1..].iter().map(|s| {
                    let value=s.parse::<f64>().map_err(|_|invalid())?;
                    if !value.is_finite() {return Err(invalid().into());}
                    Ok(value)
                }).collect()
            };
            match fields[0] {
                "material" if material.is_none()=>{let v=numbers(3)?;material=Some([v[0],v[1],v[2]]);}
                "azimuths" if azimuths.is_none() && fields.len()==2=>{
                    let value=fields[1].parse::<usize>().map_err(|_|invalid())?;
                    if !(8..=1024).contains(&value) {return Err("shell azimuths must be in 8..=1024".into());}
                    azimuths=Some(value);
                }
                "band_hz" if band.is_none()=>{let v=numbers(2)?;
                    if v[0]<=0.0 || v[1]<=v[0] {return Err("shell band needs 0 < lower < upper Hz".into());}
                    band=Some([v[0],v[1]]);
                }
                "station" if stations.len()<1024=>{let v=numbers(3)?;
                    stations.push(ProfileStation {radius_m:v[0],height_m:v[1],thickness_m:v[2]});}
                "ring" if rings.len()<1024=>{let v=numbers(4)?;
                    rings.push(AnnularRelief {radius_m:v[0],half_width_m:v[1],height_delta_m:v[2],thickness_delta_m:v[3]});}
                "dent" if dents.len()<1024=>{let v=numbers(5)?;
                    dents.push(SurfaceIndentation {center_m:[v[0],v[1]],radius_m:v[2],height_delta_m:v[3],thickness_delta_m:v[4]});}
                _=>return Err(invalid().into()),
            }
        }
        let [young_pa,poisson,density_kg_m3]=material.ok_or("shell profile requires material E_Pa,nu,rho_kg_m3")?;
        if !header || stations.len()<2 {return Err("shell profile requires at least two meridian stations".into());}
        // Physical section and geometric admission belong to revolve, below.
        Ok(Self {stations,rings,dents,young_pa,poisson,density_kg_m3,
            azimuths:azimuths.ok_or("shell profile requires azimuths")?,
            band_hz:band.ok_or("shell profile requires band_hz")?,sampled:None,default_strike:None})
    }

    /// Material and every detail enter the actual mass/stiffness geometry. Keep
    /// the existing work ceilings, and refuse visibly undersampled imported
    /// features rather than passing an omitted hammer/lathe detail off as solved.
    pub fn build(&self) -> Result<MeshShell, Error> {
        if let Some(shell) = &self.sampled { return Ok(shell.clone()); }
        let shell=revolve(&self.stations,self.azimuths,self.young_pa,self.poisson,self.density_kg_m3,
            &self.rings,&self.dents,super::mesh_budget())?;
        let inner=self.stations[0].radius_m; let outer=self.stations.last().unwrap().radius_m;
        if self.rings.iter().any(|r| r.radius_m+r.half_width_m<=inner || r.radius_m-r.half_width_m>=outer)
            || self.dents.iter().any(|d| {
                let r=d.center_m[0].hypot(d.center_m[1]);
                r+d.radius_m<=inner || r-d.radius_m>=outer
            }) {return Err("shell detail support does not intersect the supplied surface".into());}
        if shell.underresolved_features!=0 {
            return Err(format!("{} shell details underresolved: refine meridian stations/azimuths; max edge={} m",
                shell.underresolved_features,shell.max_edge_m).into());
        }
        Ok(shell.into())
    }

    /// Full supplied 3D geometry with explicit per-node thickness and materials.
    /// No meridian, alloy-name lookup or reference profile is substituted.
    pub fn parse_mesh(text: &str) -> Result<Self, Error> {
        let (shell,band_hz,strike) = mesh_input::parse(text)?;
        Ok(Self { stations:vec![],rings:vec![],dents:vec![],azimuths:0,
            young_pa:0.,poisson:0.,density_kg_m3:0.,band_hz,
            sampled:Some(shell),default_strike:Some(strike) })
    }

    pub fn default_strike_position(&self) -> Result<[f64;2],Error> {
        if let Some(p) = self.default_strike { return Ok(p); }
        let inner = self.stations.first().ok_or("shell has no default strike geometry")?.radius_m;
        let outer = self.stations.last().ok_or("shell has no default strike geometry")?.radius_m;
        Ok([(inner+2.0*outer)/3.0,0.0])
    }

    /// Export a declared revolved profile as explicit SI samples, before any
    /// modal reduction. Round-trip floats retain the original geometry exactly.
    pub fn mesh_text(&self) -> Result<String,Error> {
        if self.sampled.is_some() { return Err("profile-to-mesh export requires a profile input".into()); }
        let shell = self.build()?; let strike = self.default_strike_position()?;
        let mut text = format!("{}\n# Declared profile inputs, not measurement evidence.\n",mesh_input::HEADER);
        writeln!(text,"band_hz,{:.17e},{:.17e}",self.band_hz[0],self.band_hz[1])?;
        writeln!(text,"strike,{:.17e},{:.17e}",strike[0],strike[1])?;
        writeln!(text,"material,0,{:.17e},{:.17e},{:.17e}",self.young_pa,self.poisson,self.density_kg_m3)?;
        for (i,(p,h)) in shell.mesh.nodes.iter().zip(&shell.nodal_thickness_m).enumerate() {
            writeln!(text,"node,{i},{:.17e},{:.17e},{:.17e},{:.17e}",p[0],p[1],p[2],h)?;
        }
        for (i,t) in shell.mesh.tris.iter().enumerate() {
            writeln!(text,"triangle,{i},{},{},{},0",t[0],t[1],t[2])?;
        }
        Ok(text)
    }

    pub fn input_label(&self) -> &'static str {
        if self.sampled.is_some() { "supplied 3D mesh" } else { "supplied profile" }
    }
}

/// Extract geometry selection before the existing playing parser. Reject an
/// unrelated command before file I/O; this option never changes drum physics.
pub fn option(args: &mut Vec<String>) -> Result<Option<PathBuf>,Error> {
    path_option(args,"--shell-profile")
}
pub fn mesh_option(args: &mut Vec<String>) -> Result<Option<PathBuf>,Error> {
    path_option(args,"--shell-mesh")
}
fn path_option(args: &mut Vec<String>, flag: &str) -> Result<Option<PathBuf>,Error> {
    let mut positions=args.iter().enumerate().filter(|(_,a)| a.as_str()==flag);
    let Some((i,_))=positions.next() else {return Ok(None);};
    if positions.next().is_some() || i+1==args.len() || args[i+1].starts_with("--") {
        return Err(format!("{flag} requires exactly one input path").into());
    }
    let path=PathBuf::from(&args[i+1]);args.drain(i..i+2);Ok(Some(path))
}
pub fn admit_command(path: Option<&Path>,command: &str) -> Result<(),Error> {
    if path.is_some() && !matches!(command,"splash"|"splash-wav"|"splash-mic") {
        return Err("--shell-profile applies to the splash shell execution path, not drum or snare commands".into());
    }
    Ok(())
}

/// Reject ambiguous geometry selections before opening either input.
pub fn admit_selection(profile: Option<&Path>, mesh: Option<&Path>, command: &str) -> Result<(),Error> {
    if profile.is_some() && mesh.is_some() { return Err("choose --shell-profile or --shell-mesh, not both".into()); }
    if mesh.is_some() && !matches!(command,"splash"|"splash-wav"|"splash-mic") {
        return Err("--shell-mesh applies only to the physical splash shell path".into());
    }
    admit_command(profile,command)
}

pub fn load_selection(profile: Option<&Path>, mesh: Option<&Path>) -> Result<Option<Specimen>,Error> {
    if profile.is_some() && mesh.is_some() { return Err("choose one explicit shell geometry input".into()); }
    if let Some(path) = mesh {
        let mut text = String::new();
        std::fs::File::open(path)?.take((mesh_input::MAX_BYTES+1) as u64).read_to_string(&mut text)?;
        return Specimen::parse_mesh(&text).map(Some);
    }
    profile.map(Specimen::load).transpose()
}

/// Cold, no-clobber export front door. It performs no eigensolve or playback.
pub fn export_command(args: &[String]) -> Result<bool,Error> {
    if args.first().map(String::as_str) != Some("export-shell-mesh") { return Ok(false); }
    if !(2..=3).contains(&args.len()) {
        return Err("usage: percussion export-shell-mesh OUTPUT.fss [INPUT.profile]".into());
    }
    let input = if args.len() == 3 { Specimen::load(Path::new(&args[2]))? } else { Specimen::reference() };
    let text = input.mesh_text()?;
    let mut file = std::fs::OpenOptions::new().write(true).create_new(true).open(&args[1])?;
    file.write_all(text.as_bytes())?; file.flush()?;
    eprintln!("wrote explicit shell mesh to {}; physical values retain their input/estimated status, not a measured specimen",args[1]);
    Ok(true)
}

#[cfg(test)]
mod tests {
    use super::*;
    const INPUT:&str=include_str!("estimated-splash.profile");

    #[test]
    fn external_reference_reconstructs_the_identical_mechanical_and_acoustic_geometry() {
        let expected=Specimen::reference().build().unwrap();
        let actual=Specimen::parse(INPUT).unwrap().build().unwrap();
        assert_eq!(actual.mesh.nodes,expected.mesh.nodes);assert_eq!(actual.mesh.tris,expected.mesh.tris);
        assert_eq!(actual.nodal_thickness_m,expected.nodal_thickness_m);
        assert_eq!(actual.mass_kg.to_bits(),expected.mass_kg.to_bits());
        for (a,b) in actual.sections.iter().zip(&expected.sections) {assert_eq!(a.d,b.d);}
    }
    #[test]
    fn supplied_thickness_material_and_details_change_physical_mass_and_stiffness() {
        let old=Specimen::parse(INPUT).unwrap().build().unwrap();
        let mut supplied=Specimen::parse(INPUT).unwrap();
        for s in &mut supplied.stations {s.thickness_m*=2.0;}
        let thick=supplied.build().unwrap();
        assert_eq!(thick.mesh.nodes,old.mesh.nodes);
        assert!((thick.mass_kg/old.mass_kg-2.0).abs()<1e-12);
        assert!((thick.sections[0].d[0]/old.sections[0].d[0]-8.0).abs()<1e-12);
        supplied.young_pa*=1.5;supplied.density_kg_m3*=0.5;
        let material=supplied.build().unwrap();
        assert!((material.sections[0].d[0]/thick.sections[0].d[0]-1.5).abs()<1e-12);
        assert!((material.mass_kg/thick.mass_kg-0.5).abs()<1e-12);
        let detailed=Specimen::parse(&format!("{INPUT}ring,0.04,0.1,0.0001,-0.00001\ndent,0.04,0,0.1,-0.0002,0.00002\n"))
            .unwrap().build().unwrap();
        assert_ne!(detailed.mesh.nodes,old.mesh.nodes);assert_ne!(detailed.nodal_thickness_m,old.nodal_thickness_m);
        assert_ne!(detailed.mass_kg,old.mass_kg);
    }
    #[test]
    fn malformed_or_unresolved_data_cannot_silently_fall_back_to_the_reference() {
        for suffix in ["material,1e9,0.3,1000", "azimuths,64", "band_hz,20,2000", "dnet,1,2,3,4,5",
            "station,NaN,0,0.001", "station,0.2,0,", "ring,0.1,1,2"] {
            assert!(Specimen::parse(&format!("{INPUT}{suffix}\n")).is_err());
        }
        for text in [INPUT.replace("material,112.6e9,0.342,8607\n",""),
            INPUT.replace("band_hz,50,1200","band_hz,1200,50"),
            INPUT.replace("azimuths,32","azimuths,1000000")] {assert!(Specimen::parse(&text).is_err());}
        for suffix in ["dent,0.04,0,0.0001,-0.0001,0", "ring,2,0.1,0,0.001", "station,0.1,0,0.001"] {
            assert!(Specimen::parse(&format!("{INPUT}{suffix}\n")).unwrap().build().is_err());
        }
        assert!(Specimen::parse(&" ".repeat(MAX_BYTES+1)).is_err());
    }
    #[test]
    fn selection_preserves_physical_controls_and_rejects_unrelated_commands() {
        let mut args:Vec<String>=["splash-mic","100","--shell-profile","scan.txt","--strike-speed-m-s","4"].map(String::from).to_vec();
        let path=option(&mut args).unwrap().unwrap();assert_eq!(path,Path::new("scan.txt"));
        let (args,stroke)=super::super::playing::parse(args).unwrap();
        assert_eq!(args,["splash-mic","100"]);assert_eq!(stroke.speed_m_s,4.0);
        for command in ["splash","splash-wav","splash-mic"] {admit_command(Some(&path),command).unwrap();}
        for command in ["drum","drum-stretch","snare"] {assert!(admit_command(Some(&path),command).is_err());}
        for text in ["--shell-profile","--shell-profile --strike-speed-m-s 4","--shell-profile a --shell-profile b"] {
            let mut args:Vec<_>=text.split_whitespace().map(String::from).collect();let before=args.clone();
            assert!(option(&mut args).is_err());assert_eq!(args,before);
        }
    }

    #[test]
    fn supplied_shell_reaches_nonlinear_motion_and_rejects_displaced_hardware() {
        let input=Specimen::parse(INPUT).unwrap();
        let mut experiment=super::super::splash_with_specimen(192,2e-6,false,
            super::super::Stroke::default(),Some(input)).unwrap();
        let initial=experiment.system.state().to_vec();
        experiment.system=experiment.system.into_prepared_nonlinear().unwrap();
        let gate=fs_exec::CancelGate::new_clock_free();let mut peak=0.0_f64;
        for _ in 0..192 {
            let frame=experiment.system.step(&experiment.force,&gate).unwrap();
            assert!(frame.balance_residual_j.abs()<1e-7);
            let x=experiment.system.state();
            let motion=experiment.observer_a.iter().enumerate().map(|(i,b)|b*x[2*i]).sum::<f64>();
            peak=peak.max(motion.abs());
        }
        assert!(peak>0.0);assert_ne!(experiment.system.state(),initial);
        // The host's pads are at radius 12 mm; a wider mounting hole cannot
        // silently move those pads onto the first available metal triangle.
        let mut absent=Specimen::reference();absent.stations.drain(..2);
        let shell=absent.build().unwrap();
        assert!(super::super::playing::shell_location(&shell.mesh.nodes,&shell.mesh.tris,[0.012,0.0]).is_err());
    }
}
