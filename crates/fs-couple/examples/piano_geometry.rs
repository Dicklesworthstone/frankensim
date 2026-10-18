//! Export the existing source-derived Model D at explicit mechanical resolution.
//! The production binary is std-only; no eigensolve or alternate physics model.
//! Full builder regressions live in grand_piano, alongside its actual reader.
#[cfg(not(test))]
#[path = "grand_piano/steinway_d.rs"]
mod steinway_d;

use std::collections::BTreeSet;
use std::io::Write;
use std::path::Path;

const USAGE: &str = "piano_geometry [--divisions 4..32]
    [--rib-normal-step-m positive_metres (admits divisions through 64)]
    [--bridge-map estimated|rt0425] [--format fsb|obj]
    [--panel-rotary-mass lumped|consistent-cubic]
    [--stiffener-mass lumped|consistent-hermite|consistent-eccentric] --output fresh_path
Defaults: 8 transverse divisions, original row spacing, rt0425 bridge map, FSB.
The rt0425 map projects the source's published coupling points onto the bridges.
The estimated map (stations spaced evenly along each bridge) puts most agraffes
over the soundboard itself; it only reproduces earlier boards.
The step is the maximum distance perpendicular to the ribs between mesh rows.
All source feature stations are retained; this is not acoustic-only refinement.
Output contains the same source-derived panel/ribs/bridges and estimated material
fields, resampled at the requested mesh. It is not measured or convergence-certified.
FSB feeds grand_piano --board-geometry PATH and the existing exterior loader.
The mass options select existing physical operators in FSB, not material changes.
They require FSB; consistent-cubic also needs --edge-cubic-board-mass at rendering.
No mode cap, source band, string tension or rendering default is changed.";

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Format { Fsb, Obj }

#[derive(Debug)]
struct Options {
    divisions: usize,
    rib_normal_step_m: Option<f64>,
    published_contacts: bool,
    format: Format,
    output: String,
    panel_rotary_mass: Option<&'static str>,
    stiffener_mass: Option<&'static str>,
}
impl Options {
    fn parse(args: &[String]) -> Result<Self, String> {
        let mut options=Self {divisions:8,rib_normal_step_m:None,
            published_contacts:true,format:Format::Fsb,output:String::new(),
            panel_rotary_mass:None,stiffener_mass:None};
        let mut seen=BTreeSet::new();
        let mut args=args.iter();
        while let Some(flag)=args.next() {
            if !seen.insert(flag.as_str()) {return Err(format!("duplicate option {flag}"));}
            let value=args.next().filter(|s|!s.starts_with("--"))
                .ok_or_else(||format!("missing value for {flag}"))?;
            match flag.as_str() {
                "--divisions" => options.divisions=value.parse().map_err(|_|"invalid division count")?,
                "--rib-normal-step-m" => options.rib_normal_step_m=Some(value.parse()
                    .map_err(|_|"invalid rib-normal step")?),
                "--bridge-map" => options.published_contacts=match value.as_str() {
                    "estimated"=>false,"rt0425"=>true,_=>return Err("bridge map must be estimated or rt0425".into()),
                },
                "--format" => options.format=match value.as_str() {
                    "fsb"=>Format::Fsb,"obj"=>Format::Obj,_=>return Err("format must be fsb or obj".into()),
                },
                "--panel-rotary-mass" => options.panel_rotary_mass=Some(match value.as_str() {
                    "lumped"=>"lumped", "consistent-cubic"=>"consistent-cubic",
                    _=>return Err("panel rotary mass must be lumped or consistent-cubic".into()),
                }),
                "--stiffener-mass" => options.stiffener_mass=Some(match value.as_str() {
                    "lumped"=>"lumped", "consistent-hermite"=>"consistent-hermite",
                    "consistent-eccentric"=>"consistent-eccentric",
                    _=>return Err("stiffener mass must be lumped, consistent-hermite or consistent-eccentric".into()),
                }),
                "--output" => options.output=value.clone(),
                _=>return Err(format!("unknown option {flag}")),
            }
        }
        let maximum=if options.rib_normal_step_m.is_some(){64}else{32};
        if !(4..=maximum).contains(&options.divisions) {
            return Err(format!("divisions require 4..{maximum}; counts above 32 need an explicit rib-normal step"));
        }
        if options.rib_normal_step_m.is_some_and(|step|!step.is_finite() || step<=0.0) {
            return Err("rib-normal step must be finite positive metres".into());
        }
        if options.format!=Format::Fsb
            && (options.panel_rotary_mass.is_some() || options.stiffener_mass.is_some()) {
            return Err("physical mass selections require FSB; OBJ cannot retain kinetic operators".into());
        }
        if options.output.is_empty() {return Err("a nonempty --output path is required".into());}
        Ok(options)
    }

    /// Add only explicit selections understood by the shared board reader.
    /// Source material/geometry rows are unchanged. No mass contribution is
    /// assembled here, and an existing selection is never silently overwritten.
    fn apply_mass_selection(&self, text: &mut String) -> Result<(),String> {
        let rows=[("panel-rotary-mass",self.panel_rotary_mass),
            ("stiffener-mass",self.stiffener_mass)];
        if rows.iter().all(|(_,value)|value.is_none()) { return Ok(()); }
        if self.format!=Format::Fsb {
            return Err("physical mass selections require FSB".into());
        }
        for (tag,value) in &rows {
            if value.is_some() && text.lines().any(|line| {
                let line=line.trim();
                !line.starts_with('#') && line.split(',').next().map(str::trim)==Some(*tag)
            }) {
                return Err(format!("source geometry already selects {tag}; no override was applied"));
            }
        }
        // Validate every requested row before changing text, so an error in the
        // second row cannot leave the first selection partially applied.
        if !text.ends_with('\n') {text.push('\n');}
        for (tag,value) in rows {
            if let Some(value)=value {text.push_str(&format!("{tag},{value}\n"));}
        }
        Ok(())
    }
}

fn write_fresh(path: &Path, bytes: &[u8]) -> Result<(),String> {
    let mut output=std::fs::OpenOptions::new().write(true).create_new(true).open(path)
        .map_err(|e|format!("create {}: {e}",path.display()))?;
    output.write_all(bytes).and_then(|()|output.sync_all())
        .map_err(|e|format!("write {}: {e}; partial output may remain",path.display()))
}

#[cfg(not(test))]
fn run(options: &Options) -> Result<(),String> {
    // Complete physical geometry and capacity admission before opening output.
    let preset=match options.rib_normal_step_m {
        Some(step)=>steinway_d::build_resolved(options.divisions,step,options.published_contacts),
        None if options.published_contacts=>steinway_d::build_with_rt0425_contacts(options.divisions),
        None=>steinway_d::build_with_estimated_contacts(options.divisions),
    }?;
    let nodes=preset.geometry.lines().filter(|s|s.starts_with("node,")).count();
    let triangles=preset.geometry.lines().filter(|s|s.starts_with("triangle,")).count();
    let mut text=match options.format {Format::Fsb=>preset.geometry,Format::Obj=>preset.obj};
    options.apply_mass_selection(&mut text)?;
    write_fresh(Path::new(&options.output),text.as_bytes())?;
    println!("Exported {}: {nodes} mechanical nodes, {triangles} panel triangles, 88 bridge stations",
        options.output);
    println!("Transverse divisions {}; rib-normal step {:?} m; bridge map {}",
        options.divisions,options.rib_normal_step_m,
        if options.published_contacts {"RT-0425 projected"} else {"estimated"});
    if options.panel_rotary_mass.is_some() || options.stiffener_mass.is_some() {
        println!("Panel rotary inertia {:?}; stiffener inertia {:?} (None preserves reader defaults)",
            options.panel_rotary_mass,options.stiffener_mass);
    }
    println!("Source-derived geometry with estimated fields; mesh resolution is not an accuracy certificate.");
    Ok(())
}

#[cfg(not(test))]
fn main() {
    let args:Vec<String>=std::env::args().skip(1).collect();
    if args.as_slice()==["--help"] || args.as_slice()==["-h"] {println!("{USAGE}");return;}
    if let Err(error)=Options::parse(&args).and_then(|options|run(&options)) {
        eprintln!("piano_geometry: {error}\n{USAGE}");std::process::exit(1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn parse(args: &[&str])->Result<Options,String> {
        Options::parse(&args.iter().map(|s|s.to_string()).collect::<Vec<_>>())
    }
    #[test]
    fn default_and_explicit_resolution_have_no_hidden_source_switch() {
        let default=parse(&["--output","a.fsb"]).unwrap();
        assert_eq!(default.divisions,8);assert!(default.rib_normal_step_m.is_none());
        assert!(default.published_contacts);assert_eq!(default.format,Format::Fsb);
        let fine=parse(&["--output","a.obj","--divisions","48","--rib-normal-step-m","0.02",
            "--bridge-map","estimated","--format","obj"]).unwrap();
        assert_eq!(fine.divisions,48);assert_eq!(fine.rib_normal_step_m,Some(0.02));
        assert!(!fine.published_contacts);assert_eq!(fine.format,Format::Obj);
        assert!(parse(&["--output","a.fsb","--bridge-map","rt0425"]).unwrap().published_contacts);
        assert!(USAGE.contains("not measured"));
    }
    #[test]
    fn malformed_or_ambiguous_options_refuse() {
        for tail in [vec!["--divisions","33"],vec!["--divisions","3"],
            vec!["--divisions","65","--rib-normal-step-m","0.02"],
            vec!["--rib-normal-step-m","NaN"],vec!["--rib-normal-step-m","inf"],
            vec!["--rib-normal-step-m","0"],vec!["--rib-normal-step-m","-1"],
            vec!["--format","wav"],vec!["--bridge-map","other"],vec!["--what","1"],
            vec!["--divisions","8","--divisions","16"],vec!["--divisions"]] {
            let mut args=vec!["--output","fresh.fsb"];args.extend(tail);
            assert!(parse(&args).is_err(),"{args:?}");
        }
        assert!(parse(&[]).is_err());assert!(parse(&["--output",""]).is_err());
    }
    #[test]
    fn physical_inertia_choices_reach_fsb_without_changing_source_rows() {
        let options=parse(&["--output","a.fsb","--panel-rotary-mass","consistent-cubic",
            "--stiffener-mass","consistent-eccentric"]).unwrap();
        let source="frankensim-board-geometry-si-v1\nsource,mixed,source drawing\n# preserved geometry\n";
        let mut text=source.to_string();
        options.apply_mass_selection(&mut text).unwrap();
        assert_eq!(text,format!("{source}panel-rotary-mass,consistent-cubic\nstiffener-mass,consistent-eccentric\n"));
        let default=parse(&["--output","a.fsb"]).unwrap();
        let mut untouched=source.to_string();default.apply_mass_selection(&mut untouched).unwrap();
        assert_eq!(untouched,source);
        // A future source builder adding either row cannot accidentally double
        // inertia or silently change the requested image.
        for tag in ["panel-rotary-mass","stiffener-mass"] {
            let mut duplicate=format!("{source}  {tag},lumped\n");
            let before=duplicate.clone();
            assert!(options.apply_mass_selection(&mut duplicate).is_err());
            assert_eq!(duplicate,before);
        }
    }
    #[test]
    fn unsupported_or_unrecordable_inertia_choices_refuse() {
        for tail in [vec!["--panel-rotary-mass","other"],vec!["--stiffener-mass","other"],
            vec!["--format","obj","--panel-rotary-mass","lumped"],
            vec!["--format","obj","--stiffener-mass","consistent-hermite"],
            vec!["--panel-rotary-mass","lumped","--panel-rotary-mass","consistent-cubic"]] {
            let mut args=vec!["--output","fresh"];args.extend(tail);
            assert!(parse(&args).is_err(),"{args:?}");
        }
        for panel in ["lumped","consistent-cubic"] {
            for beam in ["lumped","consistent-hermite","consistent-eccentric"] {
                assert!(parse(&["--output","fresh","--panel-rotary-mass",panel,
                    "--stiffener-mass",beam]).is_ok());
            }
        }
    }
    #[test]
    fn export_preserves_bytes_and_never_overwrites_an_existing_file() {
        let root=std::env::temp_dir().join(format!("piano-geometry-{}-{}",std::process::id(),
            std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos()));
        std::fs::create_dir(&root).unwrap();let path=root.join("board.fsb");
        let bytes=b"frankensim-board-geometry-si-v1\n# exact caller bytes\n";
        write_fresh(&path,bytes).unwrap();assert_eq!(std::fs::read(&path).unwrap(),bytes);
        assert!(write_fresh(&path,b"replacement").is_err());
        assert_eq!(std::fs::read(path).unwrap(),bytes);
        assert!(write_fresh(&root,b"not a file").is_err());
    }
}
