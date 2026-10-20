//! Material packs + physical geometry/loads + canonical gestures -> one coupled
//! string/body/cavity state and the existing streamed PCM writer. No preset bank.
#[path = "material_strings/input.rs"]
mod input;
#[path = "material_strings/construction.rs"]
mod construction;
#[path = "material_strings/body.rs"]
mod body;
#[path = "material_strings/stops.rs"]
mod stops;
use std::{collections::BTreeMap, error::Error, fs::{File,OpenOptions}, io::{Read,Seek,Write}, path::Path};
use fs_couple::bowed_string::runtime::reciprocal::bank::{
    MaterialStringSource,ReciprocalStringBank,StringBankOptions,StringBowLocation};
use fs_couple::bowed_string::runtime::reciprocal::bank::plate::{SharedMaterialPlateBody,SharedPlateBridge};
use fs_couple::bowed_string::runtime::reciprocal::bank::torsion::{StringBodyRotationConnection,StringBridgeLocation};
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::MaterialBowContact;
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::performance::{MultiBowPerformance,MultiBowPerformanceConfig};
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::performance::force::{StringForceCoupling,StringForcePort,StringForceProgram};
use fs_couple::bowed_string::runtime::reciprocal::plate::MaterialBodyCavity;
use fs_couple::cavity_modal::{HelmholtzModalCavity,CavityNeckLoss};
use fs_couple::modal_acoustic_time::ModalAcousticTimeBudget;
use fs_couple::pcm_wav::observation::PressureRenderer;
use fs_couple::pcm_wav::stream::{Pcm16WavStream,WavStreamSummary,render_pressure_pcm16,ScheduledWavProgress};
use fs_couple::render::schedule::force::{ForceInitialization,ForceRenderConfig,ModalForceEvent};
use fs_couple::render::schedule::force::coupled::{ModalCouplingConfig,MAX_MODAL_CONNECTIONS,PreparedModalBudget};
use fs_couple::air_path::cylinder::CylinderFlowResolution;
use fs_couple::render::schedule::force::coupled::contact::friction::ModalFrictionConfig;
use fs_couple::render::schedule::force::coupled::contact::friction::multiple::MultiFrictionConfig;
use fs_couple::string_specimen::{with_uniform_circular_material_state,
    with_uniform_circular_material_and_constraints, StringGeometryConstraint, StringPrestress,
    KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY as ETA};
use fs_couple::thin_plate::{with_uniform_plate_material_state,PlateMaterialModel,PlateThicknessConstraint};
use fs_exec::{CancelGate,Cx};
use fs_matdb::{NormalizedPack,NormalizedInterfacePack,MaterialCard,MaterialStateId,QueryPoint};
use fs_material::gas::GasState;
use fs_material::state_point::{ResolvedMaterialStatePoint,ScalarPropertyRequirement as Requirement,
    ScalarAdmissibility as Admit,MaterialPropertySelection as Pick,resolve_material_state_point};
use fs_qty::{Density,Pressure,DynViscosity,Dims};
use fs_scenario::{PrestressedString,ThinPlate,HelmholtzCavity,gesture::{GestureSchedule,GestureTarget,GestureValue}};
use fs_tribo::FrictionLaw;

type Result<T> = std::result::Result<T,Box<dyn Error>>;
const SHEAR: &str = "shear_modulus";
const SHEAR_ETA: &str = "kelvin_voigt_shear_viscosity";
fn modal_budget() -> ModalAcousticTimeBudget { ModalAcousticTimeBudget::audible_reference() }
fn poll(gate: &CancelGate) -> Result<()> { if gate.is_requested() { Err("material-string construction cancelled".into()) } else { Ok(()) } }
fn options(spec:&input::Spec,bow:StringBowLocation,tensor_body:bool,
    additional_body_connections:usize) -> Result<StringBankOptions> {
    let (max_modes,max_connections,max_setup_terms)=spec
        .network_dimensions_with_body_connections(tensor_body,additional_body_connections)?;
    Ok(StringBankOptions { sample_rate_hz:spec.rate,bow,string_budget:modal_budget(),
        network:ModalCouplingConfig { max_modes,max_connections,max_setup_terms,
            nyquist_guard_fraction:0.9,maximum_total_energy_j:100.0,maximum_abs_pressure_pa:1e6,
            maximum_abs_connection_force_n:1000.0,solve_relative_tolerance:1e-10,
            energy_absolute_tolerance_j:1e-12,energy_relative_tolerance:1e-8 },
        friction:ModalFrictionConfig { max_iterations:96,maximum_force_n:100.0,maximum_slip_speed_m_s:10.0,
            force_absolute_tolerance_n:1e-10,force_relative_tolerance:1e-9,sticking_tolerance_m:1e-12 } })
}
fn template(w:&input::Wire) -> PrestressedString {
    // Extension-derived tension is filled by the original material/support owner.
    let tension = match w.prestress { StringPrestress::FixedTension(t) => t, _ => 0.0 };
    PrestressedString { length_m:w.length,tension_n:tension,lin_density_kg_m:1.0,axial_stiffness_n:1.0,
        width_m:1.0,n_modes:w.bending_modes,damping_ratio:0.0,rayleigh:None,bending_stiffness_n_m2:0.0,
        kelvin_voigt_bending:None,relaxation_bending:None,polarization_detune:0.0,moving_end:false }
}
// The retained pack includes normalization/uncertainty/rights metadata. Its ID
// labels source identity, not a new assertion about chemistry or measured quality.
struct CatalogEntry { pack: NormalizedPack, card: MaterialCard }
type Catalog = BTreeMap<String,CatalogEntry>;
fn resolve(name:&str, catalog:&mut Catalog, loader:&mut impl FnMut(&str)->Result<NormalizedPack>,
    spec:&input::Spec, body:bool, twist:bool) -> Result<ResolvedMaterialStatePoint>
{
    let mut requirements=vec![Requirement::try_new("density",Density::DIMS,Admit::StrictlyPositive)?,
        Requirement::try_new("young_modulus",Pressure::DIMS,Admit::StrictlyPositive)?,
        Requirement::try_new(ETA,DynViscosity::DIMS,Admit::Finite)?];
    if body { requirements.push(Requirement::try_new("poisson_ratio",Dims::NONE,Admit::Finite)?); }
    if twist {
        requirements.push(Requirement::try_new(SHEAR,Pressure::DIMS,Admit::StrictlyPositive)?);
        requirements.push(Requirement::try_new(SHEAR_ETA,DynViscosity::DIMS,Admit::Finite)?);
    }
    resolve_properties(name,catalog,loader,spec,&requirements)
}
fn resolve_properties(name:&str, catalog:&mut Catalog,
    loader:&mut impl FnMut(&str)->Result<NormalizedPack>, spec:&input::Spec,
    requirements:&[Requirement]) -> Result<ResolvedMaterialStatePoint>
{
    if !catalog.contains_key(name) {
        let pack=loader(name)?;
        let card=MaterialCard::assemble(MaterialStateId { chemistry:format!("source-pack:{}",pack.pack_id()),
            phase:"declared solid model".into(),process:"frozen source-pack realization".into(),revision:0 },
            pack.claims().clone(),vec![])?;
        catalog.insert(name.into(),CatalogEntry { pack,card });
    }
    // This compact input profile explicitly declares bare SI T and omega axes.
    // Incompatible typed/context axes refuse in the original resolver; they are
    // never removed, guessed, coerced or supplied from an unrelated source.
    let query=QueryPoint::new().with("T",spec.temperature)?.with("omega",spec.omega)?;
    Ok(resolve_material_state_point(&catalog[name].card,&query,requirements,Pick::SingleClaimOnly)?)
}
struct Ready { spec: input::Spec, catalog: Catalog, performance: MultiBowPerformance,
    body_profile: Option<body::Profile>, rotational_bridge_baselines_m: Vec<f64>,
    stop_profile: Option<stops::Profile> }
impl input::Spec {
    fn compile(self,cx:&Cx<'_>,gate:&CancelGate,loader:impl FnMut(&str)->Result<NormalizedPack>)->Result<Ready> {
        self.compile_with_body(cx,gate,None,loader)
    }
    fn compile_with_body(self,cx:&Cx<'_>,gate:&CancelGate,profile:Option<body::Profile>,
        loader:impl FnMut(&str)->Result<NormalizedPack>)->Result<Ready> {
        self.compile_with_profiles(cx,gate,profile,None,loader)
    }
    fn compile_with_profiles(self,cx:&Cx<'_>,gate:&CancelGate,profile:Option<body::Profile>,
        stop_profile:Option<stops::Profile>,
        mut loader:impl FnMut(&str)->Result<NormalizedPack>)->Result<Ready> {
        poll(gate)?;
        if let Some(body)=&profile { body.preflight(&self)?; }
        if let Some(stops)=&stop_profile { stops.validate_performance(&self,gate)?; }
        // Count every fluid inertia and physical connection before pack I/O.
        // The exact owner still applies its stricter connected-island bound.
        let additional_body_connections = profile.as_ref()
            .map_or(0, body::Profile::additional_body_connections);
        self.network_dimensions_with_body_connections(profile.is_some(), additional_body_connections)?;
        let mut baselines = Vec::with_capacity(self.rotational_bridges.len());
        for joint in &self.rotational_bridges {
            let body = profile.as_ref().ok_or("rotational bridges require --body-profile with distinct physical bridge footprints")?;
            baselines.push(body.bridge_separation_m(&self, joint.positive_bridge_string, joint.negative_bridge_string)?);
        }
        let interface_bytes = profile.as_ref().map(body::Profile::interface_pack_bytes)
            .transpose()?.unwrap_or(0);
        let mut remaining = input::MAX_TOTAL_PACK_BYTES.checked_sub(interface_bytes)
            .ok_or("combined material and interface pack byte budget exceeded")?;
        // Keep the legacy bulk callback, with one combined cold source budget.
        // Catalog caching invokes it once per original bulk path.
        let mut loader = |name: &str| -> Result<NormalizedPack> {
            poll(gate)?;
            let pack = loader(name)?;
            poll(gate)?;
            let size = pack.to_bytes().len();
            if size > input::MAX_PACK_BYTES || size > remaining {
                return Err("combined material and interface pack byte budget exceeded".into());
            }
            remaining -= size;
            Ok(pack)
        };
        let mut catalog=Catalog::new(); let mut sources=Vec::with_capacity(self.strings.len());
        for wire in &self.strings {
            poll(gate)?;
            sources.push(wire.resolve_source(&mut catalog, &mut loader, &self, gate)?);
        }
        let mut contacts=Vec::with_capacity(self.schedule.tracks().len());
        for track in self.schedule.tracks() {
            let GestureTarget::BowStroke { string }=track.target else { return Err("non-bow track".into()); };
            let GestureValue::Bow { station,.. }=track.initial else { return Err("non-bow value".into()); };
            contacts.push(StringBowLocation { string:usize::from(string),station_fraction:station });
        }
        let settings=options(&self,contacts[0],profile.is_some(),additional_body_connections)?;
        let bank=if let Some(body_profile)=&profile {
            body_profile.compile(&self,sources,&mut catalog,&mut loader,settings.clone(),cx,gate)?
        } else {
            let material=resolve(&self.body.pack,&mut catalog,&mut loader,&self,true,false)?;
            let body=SharedMaterialPlateBody {
                specimen:with_uniform_plate_material_state(ThinPlate { length_m:self.body.length,width_m:self.body.width,
                    thickness_m:1.0,density_kg_m3:1.0,e1_pa:1.0,e2_pa:1.0,nu12:0.0,g12_pa:1.0,material_angle_rad:0.0,
                    damping_ratio:0.0,thermoelastic:None,kelvin_voigt_bending:None,n_modes:self.body.modes,
                    geometric_nonlinearity:false,pretension_n_m:0.0,clamped:false },&material,PlateMaterialModel::Isotropic,
                    PlateThicknessConstraint::FixedThickness(self.body.thickness))?.with_kelvin_voigt_bending_loss()?,
                bridges:self.strings.iter().map(|w|w.bridge).collect(),modal_budget:modal_budget(),listener_m:self.body.listener,
            };
            match self.cavity.clone() {
                Some(cavity)=>ReciprocalStringBank::from_material_plate_with_cavity(cx,sources,self.gas,body,cavity,
                    self.law.clone(),settings.clone(),gate)?,
                None=>ReciprocalStringBank::from_material_plate(cx,sources,self.gas,body,self.law.clone(),settings.clone(),gate)?,
            }
        };
        // Select the active prepared owner before either attachment can enlarge
        // its connected linear block. No record means the ordinary owner.
        let bank = if let Some(budget) = self.prepared_linear {
            bank.with_prepared_linear_propagation(budget, gate)?
        } else { bank };
        let bank = if self.rotational_bridges.is_empty() { bank } else {
            let mut joints = Vec::with_capacity(self.rotational_bridges.len());
            for (joint, &separation_m) in self.rotational_bridges.iter().zip(&baselines) {
                poll(gate)?;
                let rotation = bank.source().rotation_between_bridges(
                    StringBridgeLocation { string: joint.positive_bridge_string, connection: 0 },
                    StringBridgeLocation { string: joint.negative_bridge_string, connection: 0 }, separation_m,
                )?;
                joints.push(StringBodyRotationConnection {
                    string: joint.string, station_fraction: joint.station_fraction, body: rotation,
                    stiffness_n_m_per_rad: joint.stiffness_n_m_per_rad,
                    damping_n_m_s_per_rad: joint.damping_n_m_s_per_rad, rest_angle_rad: joint.rest_angle_rad,
                });
            }
            bank.with_rotational_body_connections(joints, gate)?
        };
        let bank = if let Some(resolution) = self.cylinder_gas_loading {
            bank.with_cylinder_gas_loading(resolution, gate)?
        } else { bank };
        let additional=contacts[1..].iter().map(|location|MaterialBowContact { location:*location,
            law:self.law.clone(),friction:settings.friction }).collect();
        let bank=bank.with_additional_bows(additional,MultiFrictionConfig { max_ports:16,max_sweeps:128,
            max_setup_terms:1_048_576 },gate)?;
        let bank=match &stop_profile { Some(stops)=>stops.bind(bank,gate)?, None=>bank };
        let names=self.schedule.tracks().iter().map(|t|t.id.clone()).collect();
        let mut performance=bank.into_scheduled(self.schedule.clone(),names,MultiBowPerformanceConfig {
            samples:self.samples,max_block:self.block,max_work:input::MAX_WORK,max_events:input::MAX_CONTROLS })?;
        if let Some(forces)=&self.forces { performance=performance.with_string_forces(forces.clone(),gate)?; }
        if let Some(stops)=&stop_profile { performance=stops.bind_motion(performance,gate)?; }
        performance.validate_sample_rate(self.rate)?;performance.validate_sample_count(self.samples)?;poll(gate)?;
        Ok(Ready { spec:self,catalog,performance,body_profile:profile,rotational_bridge_baselines_m:baselines,stop_profile })
    }
}
impl Ready {
    fn write<W:Write+Seek>(&mut self,output:W,gate:&CancelGate)->Result<(W,WavStreamSummary)> {
        let mut stream=Pcm16WavStream::new(output,self.spec.rate,self.spec.full_scale,self.spec.block)?;
        let mut scratch=vec![0.0;self.spec.block];
        let progress=render_pressure_pcm16(&mut self.performance,&mut stream,gate,&mut scratch,self.spec.samples)?;
        let result=stream.finish()?;
        if !matches!(progress,ScheduledWavProgress::Completed { .. }) {
            return Err(format!("cancelled: only {} samples finalized",result.1.samples).into());
        }
        Ok(result)
    }
}
fn read_bounded(path:&Path,cap:usize)->Result<Vec<u8>> {
    let mut bytes=Vec::new();File::open(path)?.take(cap as u64+1).read_to_end(&mut bytes)?;
    if bytes.len()>cap {return Err(format!("{} exceeds {cap} bytes",path.display()).into());} Ok(bytes)
}
fn main()->Result<()> {
    let args:Vec<_>=std::env::args_os().collect();
    if args.len()<3 {
        return Err("usage: material_strings PERFORMANCE NEW_OUTPUT.wav [--body-profile PROFILE] [--stop-profile PROFILE]".into());
    }
    let (body_path,stop_path)=stops::profile_paths(&args[3..])?;
    let mut profile=body_path.map(|p|body::Profile::parse(&read_bounded(p,body::MAX_BYTES)?)).transpose()?;
    let stop_profile=stop_path.map(|p|stops::Profile::parse(&read_bounded(p,stops::MAX_BYTES)?)).transpose()?;
    let path=Path::new(&args[1]);let parent=path.parent().unwrap_or(Path::new("."));
    let spec=input::Spec::parse(&read_bounded(path,input::MAX_INPUT_BYTES)?)?;
    let gate=CancelGate::new_clock_free();let pool=fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let mut remaining=input::MAX_TOTAL_PACK_BYTES;
    if let Some(stops) = &stop_profile { stops.validate_performance(&spec,&gate)?; }
    if let Some(body) = &mut profile {
        body.load_interfaces(&spec,&gate,|name| {
            let bytes=read_bounded(&parent.join(name),input::MAX_PACK_BYTES.min(remaining))?;
            remaining=remaining.checked_sub(bytes.len()).ok_or("combined source pack byte budget exceeded")?;
            Ok(NormalizedInterfacePack::from_bytes(&bytes)?)
        })?;
    }
    let mut ready=pool.scope(|arena| {
        let cx=Cx::new(&gate,arena,fs_exec::StreamKey {seed:0,kernel_id:0,tile:0,iteration:0},
            fs_exec::Budget::INFINITE,fs_exec::ExecMode::Deterministic);
        spec.compile_with_profiles(&cx,&gate,profile,stop_profile,|name| {
            let bytes=read_bounded(&parent.join(name),input::MAX_PACK_BYTES.min(remaining))?;
            remaining=remaining.checked_sub(bytes.len()).ok_or("material pack byte budget exceeded")?;
            Ok(NormalizedPack::from_bytes(&bytes)?)
        })
    })?;
    for (path,source) in &ready.catalog {
        eprintln!("source={path:?}; pack={:?}; hash={}; artifact={}",source.pack.pack_id(),
            source.pack.content_hash(),source.pack.source_artifact());
    }
    if let Some(profile)=&ready.body_profile { eprintln!("body_profile={profile:?}"); }
    if let Some(profile)=&ready.stop_profile { eprintln!("stop_profile_bytes={}; stops={}",
        profile.bytes.len(),ready.performance.bank().source().stops().len()); }
    for (joint, baseline) in ready.spec.rotational_bridges.iter().zip(&ready.rotational_bridge_baselines_m) {
        eprintln!("rotational_bridge={joint:?}; body_footprint_separation_m={baseline}");
    }
    if let Some(budget)=ready.performance.bank().source().bank().prepared_linear_budget() {
        eprintln!("prepared_linear={budget:?}");
    }
    if let Some(loading)=ready.performance.bank().source().bank().cylinder_gas_loading() {
        let additional = ready.body_profile.as_ref()
            .map_or(0, body::Profile::additional_body_connections);
        let (modes,links,terms)=ready.spec
            .network_dimensions_with_body_connections(ready.body_profile.is_some(),additional)?;
        eprintln!("cylinder_gas_loading={:?}; complete_modes={modes}; connection_allowance={links}; setup_terms_allowance={terms}",
            loading.resolution());
    }
    eprintln!("input_bytes={}; gas={:?}; preload_J={:?}",ready.spec.bytes.len(),ready.spec.gas,ready.performance.preload_energy_j());
    // Complete source, geometry, frequency, schedule and preload admission is
    // finished before creating output. Never overwrite any existing path.
    let output=OpenOptions::new().write(true).create_new(true).open(&args[2])?;
    let (_,summary)=ready.write(output,&gate).map_err(|e|format!("WAV not completed; keep/discard its partial prefix explicitly: {e}"))?;
    eprintln!("samples={}; rate_Hz={}; clips={}; full_scale_Pa={}",summary.samples,summary.sample_rate_hz,
        summary.clipped_samples,summary.full_scale_pa);
    Ok(())
}

#[cfg(test)]
#[path = "material_strings/tests.rs"]
mod tests;
