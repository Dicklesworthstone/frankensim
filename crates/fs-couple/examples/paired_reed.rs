//! Two independently reduced blades, one real tube, and a physical WAV observer.
//! Usage: paired_reed INPUT.paired NEW_OUTPUT.wav FULL_SCALE_PA
//! Input sections are explicitly authored numerical data, not measured cane.
use fs_couple::bernoulli_aperture::network::{NetworkNode, TubeNetworkSpec, TubeSection};
use fs_couple::bernoulli_aperture::paired::{BladeState, PairedAperture, PairedApertureLimits};
use fs_couple::bernoulli_aperture::paired_network::PairedApertureNetwork;
use fs_couple::bernoulli_aperture::plate::{PlateApertureOptions, PlateApertureReduction};
use fs_couple::bernoulli_aperture::radiation::BaffledRadiationLoad;
use fs_couple::bernoulli_aperture::tube::TubeDrive;
use fs_couple::pcm_wav::{baffled::CircularOutletReceiver, encode_pcm16_wav};
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;
use fs_material::gas::GasState;
use fs_plate::{AssemblyOptions, EdgeSupport, PlateChart, PlateMesh, PlateSection, SliceOptions};
use std::io::{Read, Write};

const MAX_INPUT_BYTES: u64 = 16_384;
const MAX_SAMPLES: u64 = 1_000_000;
const SCHEMA: &str = "frankensim-paired-reed-v1";

#[derive(Clone, Debug)]
struct Input {
    rate: u32,
    samples: u64,
    gas: [f64; 3], // Kelvin, static Pa, relative humidity
    blades: [[f64; 5]; 2], // E, nu, rho, thickness, damping ratio
    pressure_pa: f64,
    attack: u64,
    release: u64,
    release_duration: u64,
}
fn value<T: std::str::FromStr>(text: &str) -> Result<T, String> {
    text.parse().map_err(|_| format!("invalid numeric source field: {text}"))
}
fn row<'a>(lines: &mut impl Iterator<Item = &'a str>, tag: &str, count: usize)
    -> Result<Vec<&'a str>, String>
{
    let fields: Vec<_> = lines.next().ok_or_else(|| format!("missing {tag} record"))?
        .split_ascii_whitespace().collect();
    if fields.len() != count + 1 || fields[0] != tag {
        return Err(format!("expected {tag} with exactly {count} fields"));
    }
    Ok(fields[1..].to_vec())
}
impl Input {
    fn parse(text: &str) -> Result<Self, String> {
        if text.len() as u64 > MAX_INPUT_BYTES { return Err("paired input exceeds byte budget".into()); }
        let mut lines = text.lines().map(|l| l.split('#').next().unwrap_or("").trim())
            .filter(|l| !l.is_empty());
        if lines.next() != Some(SCHEMA) { return Err("unsupported paired input schema".into()); }
        let audio = row(&mut lines, "audio", 2)?;
        let rate: u32 = value(audio[0])?;
        let samples: u64 = value(audio[1])?;
        if !(8000..=192000).contains(&rate) || samples == 0 || samples > MAX_SAMPLES {
            return Err("audio requires 8000..=192000 Hz and 1..=1000000 samples".into());
        }
        let g = row(&mut lines, "gas", 3)?;
        let gas = [value(g[0])?, value(g[1])?, value(g[2])?];
        let mut blades = [[0.0; 5]; 2];
        for blade in &mut blades {
            let r = row(&mut lines, "blade", 5)?;
            for (slot, source) in blade.iter_mut().zip(r) { *slot = value(source)?; }
            if !blade.iter().all(|x: &f64| x.is_finite()) || blade[4] < 0.0 {
                return Err("blade fields must be finite with nonnegative damping ratio".into());
            }
            // Keep the section owner's elastic, thickness and density domain.
            PlateSection::isotropic(blade[0], blade[1], blade[3], blade[2]).map_err(|e|e.to_string())?;
        }
        let p = row(&mut lines, "pressure", 4)?;
        let pressure_pa: f64 = value(p[0])?;
        let attack: u64 = value(p[1])?;
        let release: u64 = value(p[2])?;
        let release_duration: u64 = value(p[3])?;
        if !pressure_pa.is_finite() || pressure_pa < 0.0 || attack > release
            || release >= samples || release.checked_add(release_duration).is_none_or(|n| n >= samples) {
            return Err("pressure needs finite nonnegative Pa and ordered attack/release followed by ringdown".into());
        }
        GasState::try_new_moist_air(gas[0],gas[1],gas[2]).map_err(|e|format!("{e:?}"))?;
        if lines.next().is_some() { return Err("unexpected trailing paired source record".into()); }
        Ok(Self {rate,samples,gas,blades,pressure_pa,attack,release,release_duration})
    }
    // A physical piecewise-linear blowing pressure, sampled at the midpoint of
    // each held mechanical interval. Zero pressure does not clear stored sound.
    fn pressure(&self, sample: u64) -> f64 {
        if sample >= self.release + self.release_duration { 0.0 }
        else if sample >= self.release {
            self.pressure_pa * ((self.release+self.release_duration-sample) as f64-0.5)
                / self.release_duration as f64
        } else if sample < self.attack {
            self.pressure_pa * (sample as f64+0.5) / self.attack as f64
        } else { self.pressure_pa }
    }
}

fn plate(data: [f64; 5], gate: &CancelGate) -> Result<PlateApertureReduction,String> {
    let section=PlateSection::isotropic(data[0],data[1],data[3],data[2]).map_err(|e|e.to_string())?;
    let mesh=PlateMesh::rectangle(0.025,0.01,6,2);
    let root=mesh.nodes.iter().enumerate().filter(|(_,p)|p.0==0.0).map(|(i,_)|i).collect();
    let slit_edges=mesh.boundary_edges().into_iter().filter(|&(a,b)|
        (mesh.nodes[a].0-0.025).abs()<1e-12 && (mesh.nodes[b].0-0.025).abs()<1e-12)
        .map(|(a,b)|[a,b]).collect();
    let chart=PlateChart::with_boundary_and_regions(mesh,section,root,vec![]).map_err(|e|e.to_string())?;
    PlateApertureReduction::from_chart(chart,PlateApertureOptions {
        assembly:AssemblyOptions {pretension:0.0,support:EdgeSupport::Clamped},
        eigenvalue_window:(1.0,(core::f64::consts::TAU*450.0).powi(2)),mode_index:0,
        eigensolver:SliceOptions::default(),max_nodes:100,max_triangles:100,
        slit_edges,rest_opening_m:0.0002,damping_ratio:data[4],
        max_slit_mode_variation:0.25,max_slope:0.1,
    },gate).map_err(|e|e.to_string())
}

fn build(input: &Input, gate: &CancelGate) -> Result<PairedApertureNetwork,String> {
    let plates=[plate(input.blades[0],gate)?,plate(input.blades[1],gate)?];
    let dt=1.0/f64::from(input.rate);
    for p in &plates {
        if dt*(p.stiffness_n_m()/p.mass_kg()).sqrt()>0.3 {
            return Err("plate mode exceeds the example's 0.3-radian time-resolution allowance".into());
        }
    }
    let gas=GasState::try_new_moist_air(input.gas[0],input.gas[1],input.gas[2]).map_err(|e|format!("{e:?}"))?;
    // This is an explicitly baffled circular outlet, not an oboe bell model.
    let radius=0.007;
    let band=1000.0_f64.min(0.05*f64::from(input.rate));
    let radiation=BaffledRadiationLoad::new(radius,gas.density,gas.sound_speed,dt,band,gate)
        .map_err(|e|e.to_string())?;
    let graph=TubeNetworkSpec {nodes:vec![NetworkNode::Inlet,radiation.termination()],
        sections:vec![TubeSection {nodes:[0,1],length_m:0.25,radius_m:radius,max_length_error_m:0.005}],
        sound_speed_m_s:gas.sound_speed,max_wave_memory_bytes:1<<20};
    let z=graph.inlet_impedance_with_gases(&[gas]).map_err(|e|e.to_string())?;
    let contact=Obstacle::new(vec![-1.0],1,1,vec![0.0],vec![1.0],1e8,2.0,
        "authored shared blade contact; not a measured cane interface".into())
        .and_then(|c|c.with_internal_loss(0.5)).map_err(|e|e.to_string())?;
    let pair=PairedAperture::from_plates(plates,gas,z,dt,input.samples,
        PairedApertureLimits {maximum_abs_displacement_m:0.001,maximum_abs_velocity_m_s:10.0,
            maximum_energy_j:1.0,max_evaluations:256,velocity_tolerance_m_s:1e-15,relative_tolerance:1e-13},
        [BladeState::default();2],contact).map_err(|e|e.to_string())?;
    PairedApertureNetwork::new(pair,graph,vec![gas],gate).map_err(|e|e.to_string())
}

fn run() -> Result<(),String> {
    let args: Vec<_>=std::env::args().skip(1).collect();
    let [input_path,output_path,scale]=args.as_slice() else {
        return Err("usage: paired_reed INPUT.paired NEW_OUTPUT.wav FULL_SCALE_PA".into());
    };
    let scale:f64=value(scale)?;
    if !scale.is_finite() || scale<=0.0 {return Err("full scale must be explicit positive finite pascals".into());}
    let mut bytes=Vec::new();
    std::fs::File::open(input_path).map_err(|e|e.to_string())?.take(MAX_INPUT_BYTES+1)
        .read_to_end(&mut bytes).map_err(|e|e.to_string())?;
    let input=Input::parse(std::str::from_utf8(&bytes).map_err(|e|e.to_string())?)?;
    let gate=CancelGate::new_clock_free();
    let mut model=build(&input,&gate)?;
    let mut receiver=model.baffled_receiver(1,input.rate,CircularOutletReceiver {
        position_m:[0.0,0.0,0.25],radial_rings:4,angular_points:16,
        maximum_frequency_hz:1000.0_f64.min(0.05*f64::from(input.rate)),
    }).map_err(|e|e.to_string())?;
    for (i,p) in model.aperture().plates().expect("original independent sources").iter().enumerate() {
        eprintln!("blade {i}: mass={:.17e} kg, stiffness={:.17e} N/m, pressure_area={:.17e} m2, eigenvalue_interval={:?}",
            p.mass_kg(),p.stiffness_n_m(),p.pressure_area_m2(),p.eigenvalue_interval());
    }
    eprintln!("actual gas density={} kg/m3, speed={} m/s; represented tube={} m; output is one-way baffled exterior Pa",
        model.source_gases()[0].density,model.source_gases()[0].sound_speed,model.represented_sections()[0].represented_length_m);
    let mut pressure=Vec::with_capacity(input.samples as usize);
    let mut maximum_balance=0.0_f64;
    for i in 0..input.samples {
        let frame=model.step(TubeDrive {upstream_pressure_pa:input.pressure(i),body_flow_m3_s:0.0})
            .map_err(|e|format!("sample {i}: {e}"))?;
        maximum_balance=maximum_balance.max(frame.coupled.balance_residual_j().abs());
        let flow=model.node_frame(1).expect("admitted physical outlet").net_flow_into_node_m3_s;
        pressure.push(receiver.step(&[flow])?);
    }
    let (wav,clips)=encode_pcm16_wav(&pressure,input.rate,scale).map_err(|e|e.to_string())?;
    // No file is created until physics, receiver and encoding have all succeeded.
    // create_new refuses an existing destination; it never overwrites a run.
    let mut out=std::fs::OpenOptions::new().write(true).create_new(true).open(output_path)
        .map_err(|e|format!("new WAV destination: {e}"))?;
    out.write_all(&wav).and_then(|_|out.flush()).map_err(|e|e.to_string())?;
    eprintln!("samples={}, PCM clipped_samples={clips}, maximum_energy_residual={maximum_balance:.17e} J, final_storage={:.17e} J; no peak normalization",
        input.samples,model.stored_energy_j());
    Ok(())
}
fn main() {
    if let Err(error)=run() {eprintln!("paired reed refused: {error}");std::process::exit(1);}
}

#[cfg(test)]
mod tests {
    use super::*;
    const BASE: &str=include_str!("paired_reed/equal-blades.paired");
    const HEAVIER: &str=include_str!("paired_reed/heavier-second-blade.paired");
    #[test]
    fn source_pressure_releases_without_retuning_for_material() {
        let a=Input::parse(BASE).unwrap();let b=Input::parse(HEAVIER).unwrap();
        assert_eq!(a.blades[0],b.blades[0]);
        assert_eq!(a.blades[1][0],b.blades[1][0]);
        assert_eq!(2.0*a.blades[1][2],b.blades[1][2]);
        for i in 0..a.samples {assert_eq!(a.pressure(i),b.pressure(i));}
        assert_eq!(a.pressure(a.release+a.release_duration),0.0);
        for bad in [BASE.replace("900","NaN"),BASE.replace("101325","-1"),
            BASE.replace("pressure 5 480 4800 480","pressure 5 6000 4800 480"),
            format!("{BASE}\nextra 1\n")] {assert!(Input::parse(&bad).is_err());}
    }
    #[test]
    fn actual_plate_sources_survive_density_only_substitution_and_drive_the_outlet() {
        let gate=CancelGate::new_clock_free();
        let input=Input::parse(BASE).unwrap();let changed=Input::parse(HEAVIER).unwrap();
        let p0=plate(input.blades[0],&gate).unwrap();let p1=plate(changed.blades[1],&gate).unwrap();
        assert!((p1.mass_kg()/p0.mass_kg()-2.0).abs()<1e-7);
        assert!((p1.stiffness_n_m()/p0.stiffness_n_m()-1.0).abs()<1e-7);
        assert!((p1.pressure_area_m2()/p0.pressure_area_m2()-1.0).abs()<1e-7);
        let mut a=build(&input,&gate).unwrap();let mut b=build(&changed,&gate).unwrap();
        assert_eq!(b.aperture().plates().unwrap()[1].chart().mesh.nodes,p1.chart().mesh.nodes);
        let mut difference=0.0;
        for i in 0..256 {
            let drive=TubeDrive {upstream_pressure_pa:input.pressure(i),body_flow_m3_s:0.0};
            a.step(drive).unwrap();b.step(drive).unwrap();
            difference+=(a.node_frame(1).unwrap().net_flow_into_node_m3_s
                -b.node_frame(1).unwrap().net_flow_into_node_m3_s).abs();
        }
        assert!(difference>0.0);
        let source=b.aperture().plates().unwrap();let node=source[1].options().slit_edges[0][0];
        let state=b.aperture().blade_states()[1];
        let actual=b.aperture().nodal_motion(1,node).unwrap();
        assert_eq!(actual,source[1].nodal_motion(node,source[1].options().rest_opening_m+state.displacement_m,state.velocity_m_s).unwrap());
    }
}
