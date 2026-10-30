//! Material packs + physical geometry/loads + canonical gestures -> one coupled
//! string/body/cavity state and the existing streamed PCM writer. No preset bank.
#[path = "material_strings/input.rs"]
mod input;
#[path = "material_strings/construction.rs"]
mod construction;
#[path = "material_strings/material_inputs.rs"]
mod material_inputs;
#[path = "material_strings/body.rs"]
mod body;
#[path = "material_strings/body_response.rs"]
mod body_response;
#[path = "material_strings/coating_strain.rs"]
mod coating_strain;
#[path = "material_strings/stops.rs"]
mod stops;
#[path = "material_strings/bow_interfaces.rs"]
mod bow_interfaces;
#[path = "material_strings/shear_layers.rs"]
mod shear_layers;
#[path = "material_strings/contact_thermal.rs"]
mod contact_thermal;
#[path = "material_strings/bridge_interfaces.rs"]
mod bridge_interfaces;
#[path = "material_strings/axial_bearings.rs"]
mod axial_bearings;
#[path = "material_strings/translational_bridges.rs"]
mod translational_bridges;
#[path = "material_strings/bending_memory.rs"]
mod bending_memory;
#[path = "material_strings/torsion_memory.rs"]
mod torsion_memory;
#[path = "material_strings/fibers.rs"]
mod fibers;
// Native entry points share one bounded source loader and phase policy.
#[path = "../src/bin/music_render/thermochemical.rs"]
mod thermochemical;
#[path = "../src/bin/music_render/wind_input/attachments.rs"]
mod attachments;
#[path = "material_strings/native_input.rs"]
mod native_input;
#[path = "../src/bin/music_render/json.rs"]
mod json;
use json::json_string;
use std::{collections::BTreeMap, error::Error, fs::{File,OpenOptions}, io::{Read,Seek,Write}, path::Path};
use fs_couple::bowed_string::runtime::reciprocal::bank::{
    MaterialStringSource,ReciprocalStringBank,StringBankOptions,StringBowLocation,residual};
use fs_couple::bowed_string::runtime::reciprocal::bank::plate::{SharedMaterialPlateBody,SharedPlateBridge};
use fs_couple::bowed_string::runtime::reciprocal::bank::torsion::{StringBodyRotationConnection,StringBridgeLocation,
    SHEAR_VISCOSITY_REFERENCE_OMEGA as TORSION_REFERENCE,SHEAR_VISCOSITY_REFERENCE_OMEGA_DIMS as TORSION_REFERENCE_DIMS};
use fs_couple::bowed_string::runtime::reciprocal::bank::radiation::rayleigh::PlateRayleighReport;
use fs_couple::bowed_string::runtime::reciprocal::bank::radiation::bem::ShellBemReport;
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::{MaterialBowContact,MultiBowStringBank};
use fs_couple::bowed_string::runtime::reciprocal::bank::nonlinear::longitudinal::{LongitudinalStringResolution,LongitudinalMaterialLoss};
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::performance::{MultiBowPerformance,MultiBowPerformanceConfig};
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::performance::bow_stations::BowStationProgram;
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::performance::force::{
    MaterialForcePort, MaterialForceProgram, StringForceCoupling, StringForcePort,
};
use fs_couple::bowed_string::runtime::reciprocal::plate::MaterialBodyCavity;
use fs_couple::cavity_modal::{HelmholtzModalCavity,CavityNeckLoss};
use fs_couple::modal_acoustic_time::ModalAcousticTimeBudget;
use fs_couple::pcm_wav::observation::{DecimatedRenderer,PressureRenderer};
use fs_couple::pcm_wav::stream::{Pcm16WavStream,WavStreamSummary,render_pressure_pcm16,ScheduledWavProgress};
use fs_couple::render::schedule::force::{ForceInitialization,ForceRenderConfig,ModalForceEvent};
use fs_couple::render::schedule::force::coupled::{ModalCouplingConfig,ModalPreparation,MAX_MODAL_CONNECTIONS,PreparedModalBudget};
use fs_couple::air_path::cylinder::CylinderFlowResolution;
use fs_couple::render::schedule::force::coupled::contact::friction::ModalFrictionConfig;
use fs_couple::render::schedule::force::coupled::contact::friction::interface::DryFrictionBinding;
use fs_couple::render::schedule::force::coupled::contact::friction::nonlinear::NonlinearModalConfig;
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
use material_inputs::MaterialPack;

type Result<T> = std::result::Result<T,Box<dyn Error>>;
const SHEAR: &str = "shear_modulus";
const SHEAR_ETA: &str = "kelvin_voigt_shear_viscosity";
fn modal_budget() -> ModalAcousticTimeBudget { ModalAcousticTimeBudget::audible_reference() }
fn poll(gate: &CancelGate) -> Result<()> { if gate.is_requested() { Err("material-string construction cancelled".into()) } else { Ok(()) } }
fn options(spec:&input::Spec,bow:StringBowLocation,tensor_body:bool,
    additional_body_connections:usize,additional_body_modes:usize) -> Result<StringBankOptions> {
    let (max_modes,max_connections,max_setup_terms)=spec
        .network_dimensions_with_body_modes(tensor_body,additional_body_connections,additional_body_modes)?;
    // The N*(L+1)^2 screen leaves no room for a zero-stiffness coordinate (a cavity neck
    // plug) in a static preload: its free-block solve costs free*L^2 + free^2*L + free^3
    // with free <= L, so a full-rank viscous body plus a cavity was refused by construction.
    let max_setup_terms=max_setup_terms
        .saturating_add(max_connections.saturating_pow(3).saturating_mul(3));
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
        kelvin_voigt_bending:None,relaxation_bending:None,polarization_detune:0.0,moving_end:false,
        tension_loss_factor:w.tension_loss.as_ref().map(|row|row.eta_f) }
}
// The retained pack includes normalization/uncertainty/rights metadata. Its ID
// labels source identity, not a new assertion about chemistry or measured quality.
struct CatalogEntry { pack: MaterialPack, card: MaterialCard }
type Catalog = BTreeMap<String,CatalogEntry>;
fn resolve(name:&str, catalog:&mut Catalog, loader:&mut impl FnMut(&str)->Result<MaterialPack>,
    spec:&input::Spec, body:bool, twist:bool, loss:&input::BendingLoss) -> Result<ResolvedMaterialStatePoint>
{
    let mut requirements=vec![Requirement::try_new("density",Density::DIMS,Admit::StrictlyPositive)?,
        Requirement::try_new("young_modulus",Pressure::DIMS,Admit::StrictlyPositive)?];
    match loss {
        input::BendingLoss::KelvinVoigt =>
            requirements.push(Requirement::try_new(ETA,DynViscosity::DIMS,Admit::Finite)?),
        input::BendingLoss::Prony { properties } => {
            if body { return Err("string memory selection cannot replace a body constitutive law".into()); }
            requirements.push(Requirement::try_new(
                fs_couple::string_specimen::EQUILIBRIUM_YOUNG_MODULUS_PROPERTY,
                Pressure::DIMS,Admit::StrictlyPositive)?);
            for pair in properties {
                for (name,dims,admission) in [(&pair.modulus,Pressure::DIMS,Admit::NonNegative),
                    (&pair.relaxation_time,Dims([0,0,1,0,0,0]),Admit::StrictlyPositive)] {
                    if let Some(prior) = requirements.iter().find(|r| r.name() == name) {
                        if prior.dims() != dims {
                            return Err("selected string source property has incompatible physical quantities".into());
                        }
                    } else { requirements.push(Requirement::try_new(name.as_str(),dims,admission)?); }
                }
            }
        },
    }
    if body { requirements.push(Requirement::try_new("poisson_ratio",Dims::NONE,Admit::Finite)?); }
    if twist {
        requirements.push(Requirement::try_new(SHEAR,Pressure::DIMS,Admit::StrictlyPositive)?);
        requirements.push(Requirement::try_new(SHEAR_ETA,DynViscosity::DIMS,Admit::Finite)?);
        // A pack that wrote its shear viscosity at a stated reference frequency has constant torsional loss
        // (fs_couple::bowed_string::runtime::reciprocal::bank::torsion::SHEAR_VISCOSITY_REFERENCE_OMEGA).
        load_pack(name,catalog,loader)?;
        if !catalog[name].pack.claims().claims_for(TORSION_REFERENCE).is_empty() {
            requirements.push(Requirement::try_new(TORSION_REFERENCE,TORSION_REFERENCE_DIMS,Admit::StrictlyPositive)?);
        }
    }
    if !body && let LongitudinalMaterialLoss::KelvinVoigt { viscosity_property } = &spec.longitudinal_loss {
        requirements.push(Requirement::try_new(viscosity_property.as_str(),DynViscosity::DIMS,Admit::Finite)?);
    }
    if !body && let Some(memory)=&spec.bending_memory {
        memory.add_requirements(name,&spec.strings,&mut requirements)?;
    }
    if !body {
        for wire in &spec.strings {
            wire.add_torsion_memory_requirements(name,&mut requirements)?;
        }
    }
    resolve_properties(name,catalog,loader,spec,&requirements)
}
fn load_pack(name:&str, catalog:&mut Catalog, loader:&mut impl FnMut(&str)->Result<MaterialPack>) -> Result<()> {
    if !catalog.contains_key(name) {
        let pack=loader(name)?;
        let card=pack.card()?;
        catalog.insert(name.into(),CatalogEntry { pack,card });
    }
    Ok(())
}
fn resolve_properties(name:&str, catalog:&mut Catalog,
    loader:&mut impl FnMut(&str)->Result<MaterialPack>, spec:&input::Spec,
    requirements:&[Requirement]) -> Result<ResolvedMaterialStatePoint>
{
    load_pack(name,catalog,loader)?;
    let query=material_inputs::query(&catalog[name].card,spec.temperature,spec.omega)?;
    Ok(resolve_material_state_point(&catalog[name].card,&query,requirements,Pick::SingleClaimOnly)?)
}
#[derive(Debug)]
enum BodyRadiationReport {
    Rayleigh(PlateRayleighReport),
    ShellBem(ShellBemReport),
}
impl core::fmt::Display for BodyRadiationReport {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Rayleigh(report) => write!(formatter, "{report:?}"),
            Self::ShellBem(report) => write!(formatter, "{report:?}"),
        }
    }
}
struct Ready { spec: input::Spec, catalog: Catalog, performance: MultiBowPerformance,
    body_profile: Option<body::Profile>, rotational_bridge_baselines_m: Vec<f64>,
    stop_profile: Option<stops::Profile>, body_radiation: Option<BodyRadiationReport> }
/// The complete physical bank before gesture scheduling and microphone setup.
/// Both analysis and playback pass through this one source-binding sequence.
struct ReadyBank { spec: input::Spec, catalog: Catalog, bank: MultiBowStringBank,
    body_profile: Option<body::Profile>, rotational_bridge_baselines_m: Vec<f64>,
    stop_profile: Option<stops::Profile>, body_radiation: Option<BodyRadiationReport>,
    bow_station_program: Option<BowStationProgram> }
impl ReadyBank {
    fn into_performance(self, gate: &CancelGate) -> Result<Ready> {
        if matches!(self.bank.source().bank().preparation(), ModalPreparation::HarmonicOnly(_)) {
            return Err("a harmonic-only bank cannot be scheduled for playback".into());
        }
        self.spec.schedule_bank(self.bank, self.catalog, self.body_profile,
            self.rotational_bridge_baselines_m, self.stop_profile, self.body_radiation,
            self.bow_station_program, gate)
    }
}
impl input::Spec {
    fn compile<P: Into<MaterialPack>>(self,cx:&Cx<'_>,gate:&CancelGate,loader:impl FnMut(&str)->Result<P>)->Result<Ready> {
        self.compile_with_body(cx,gate,None,loader)
    }
    fn compile_with_body<P: Into<MaterialPack>>(self,cx:&Cx<'_>,gate:&CancelGate,profile:Option<body::Profile>,
        loader:impl FnMut(&str)->Result<P>)->Result<Ready> {
        self.compile_with_profiles(cx,gate,profile,None,loader)
    }
    fn compile_with_profiles<P: Into<MaterialPack>>(self,cx:&Cx<'_>,gate:&CancelGate,profile:Option<body::Profile>,
        stop_profile:Option<stops::Profile>,
        loader:impl FnMut(&str)->Result<P>)->Result<Ready> {
        self.compile_bank_with_profiles(cx,gate,profile,stop_profile,false,loader)?
            .into_performance(gate)
    }
    fn compile_bank_with_profiles<P: Into<MaterialPack>>(mut self,cx:&Cx<'_>,gate:&CancelGate,
        profile:Option<body::Profile>,stop_profile:Option<stops::Profile>,harmonic_only:bool,
        mut loader:impl FnMut(&str)->Result<P>)->Result<ReadyBank> {
        poll(gate)?;
        let preparation = if harmonic_only {
            if !profile.as_ref().is_some_and(body::Profile::has_shell_material_basis) {
                return Err("loaded body-response-only requires a curved or assembled shell profile".into());
            }
            ModalPreparation::HarmonicOnly(self.prepared_linear.ok_or(
                "loaded body-response-only requires an explicit prepared-linear resource declaration")?)
        } else { self.prepared_linear.into() };
        if self.longitudinal_waves.is_none()
            && !matches!(&self.longitudinal_loss, LongitudinalMaterialLoss::Elastic) {
            return Err("longitudinal material loss requires explicit axial coordinates".into());
        }
        if self.gas != *self.gas_record.state() {
            return Err("material-string gas coefficients must retain their original source record".into());
        }
        self.validate_cavity_gas()?;
        if self.cavity_gas_record.is_some() && self.cavity.is_none()
            && !profile.as_ref().is_some_and(body::Profile::has_spatial_cavity) {
            return Err("cavity_gas requires an actual lumped or spatial body cavity".into());
        }
        if self.cavity_molecular.is_some()
            && !profile.as_ref().is_some_and(body::Profile::has_thermal_cavity) {
            return Err("cavity_molecular requires a body profile with an original interior FE and air-thermal law".into());
        }
        self.validate_shell_observer()?;
        let bow_station_program = self.bow_station_program(gate)?;
        if let Some(layers)=self.friction.shear_layers() { layers.preflight(&self)?; }
        if let Some(body)=&profile { body.preflight(&self)?; }
        if self.one_way_observer.is_some() && profile.is_none() {
            return Err("one_way_observer requires a body profile that declares its closed acoustic skin".into());
        }
        if let Some(stops)=&stop_profile { stops.validate_performance(&self,gate)?; }
        if let Some(radiation) = self.body_radiation {
            radiation.validate(self.rate)?;
            if matches!(preparation, ModalPreparation::Prepared(_))
                && !matches!(radiation, input::BodyRadiation::Skin(_)) {
                return Err("prepared playback admits computed body radiation only through a skin load".into());
            }
            match radiation {
                input::BodyRadiation::Rayleigh(options) => {
                    if !profile.as_ref().is_some_and(body::Profile::has_planar_material_basis)
                        || (options.sides == 2 && self.cavity.is_some()) {
                        return Err("computed Rayleigh radiation needs a regional flat-plate profile and one exterior side when a cavity is modeled".into());
                    }
                },
                input::BodyRadiation::ShellBem(_) => {
                    if !profile.as_ref().is_some_and(body::Profile::has_shell_material_basis) {
                        return Err("computed BEM radiation requires a curved shell or assembled shell body profile".into());
                    }
                },
                input::BodyRadiation::Skin(_) => {
                    if !profile.as_ref().is_some_and(|body| body.has_shell_material_basis() && body.has_acoustic_skin()) {
                        return Err("a computed skin load requires a shell body profile that declares its closed acoustic skin".into());
                    }
                },
            }
        }
        // Count every fluid inertia and physical connection before pack I/O.
        // The exact owner still applies its stricter connected-island bound.
        let additional_body_connections = profile.as_ref()
            .map_or(0, body::Profile::additional_body_connections);
        let additional_body_modes = profile.as_ref()
            .map_or(0, body::Profile::additional_body_modes);
        self.network_dimensions_with_body_modes(profile.is_some(), additional_body_connections,
            additional_body_modes)?;
        let mut baselines = Vec::with_capacity(self.rotational_bridges.len());
        for joint in &self.rotational_bridges {
            let body = profile.as_ref().ok_or("rotational bridges require --body-profile with distinct physical bridge footprints")?;
            baselines.push(body.bridge_separation_m(&self, joint.positive_bridge_string, joint.negative_bridge_string)?);
        }
        let mut contacts=Vec::with_capacity(self.schedule.tracks().len());
        for track in self.schedule.tracks() {
            let GestureTarget::BowStroke { string }=track.target else { return Err("non-bow track".into()); };
            let GestureValue::Bow { station,.. }=track.initial else { return Err("non-bow value".into()); };
            contacts.push(StringBowLocation { string:usize::from(string),station_fraction:station });
        }
        let mut settings=options(&self,contacts[0],profile.is_some(),additional_body_connections,
            additional_body_modes)?;
        let memory = self.string_memory_states()?
            .checked_add(self.bending_memory.as_ref().map_or(0,bending_memory::Input::max_memory_states))
            .and_then(|n| n.checked_add(profile.as_ref().map_or(0, |body| body.material_memory_states_for(&self))))
            .and_then(|n| n.checked_add(self.bow_fibers.as_ref()
                .and_then(|fibers| fibers.memory.as_ref()).map_or(0, |memory| memory.max_states)))
            .ok_or("complete material memory extent overflow")?;
        // A prepared owner keeps the skin load's poles as its air bank beside the transition, charged to the prepared
        // byte and work budget (the library's own split); only the exact owner holds them in its state.
        let acoustic_states = if self.prepared_linear.is_some() { 0 }
            else { self.body_radiation.map_or(0, |r| 2*r.maximum_poles()) };
        let internal_states = acoustic_states
            .checked_add(memory).ok_or("complete constitutive/acoustic state extent overflow")?;
        if internal_states != 0 {
            let core = settings.network.max_modes.checked_mul(2).ok_or("material setup dimension overflow")?;
            let terms = core.checked_add(internal_states).and_then(|d| d.checked_mul(d))
                .and_then(|d| core.checked_mul(internal_states).and_then(|m| d.checked_add(m)))
                .and_then(|d| settings.network.max_setup_terms.checked_add(d))
                .ok_or("complete material setup work overflow")?;
            if (memory > 256 && self.prepared_linear.is_none()) || terms > input::MAX_WORK as usize {
                return Err("complete material memory exceeds the original state or setup-work ceiling".into());
            }
            settings.network.max_setup_terms = terms;
        }
        let bridge_interface_bytes = self.bridge_interfaces.pack_bytes()?;
        let axial_bearing_bytes = self.axial_bearings.pack_bytes()?;
        let axial_bearings = self.axial_bearings.bindings(self.temperature, self.omega)?.to_vec();
        let normal_source_bytes = stop_profile.as_ref().map(stops::Profile::source_pack_bytes)
            .transpose()?.unwrap_or(0);
        let interface_bytes = profile.as_ref().map(body::Profile::interface_pack_bytes)
            .transpose()?.unwrap_or(0).checked_add(self.friction.pack_bytes()?)
            .and_then(|bytes| bytes.checked_add(bridge_interface_bytes))
            .and_then(|bytes| bytes.checked_add(axial_bearing_bytes))
            .and_then(|bytes| bytes.checked_add(normal_source_bytes))
            .ok_or("combined interface pack byte count overflow")?;
        self.bridge_interfaces.apply(self.temperature, self.omega, &mut self.strings)?;
        let mut remaining = input::MAX_TOTAL_PACK_BYTES.checked_sub(interface_bytes)
            .ok_or("combined material and interface pack byte budget exceeded")?;
        // Every original artifact counts in the same cold source budget,
        // including a canonical card's complete associated model pack.
        // Catalog caching invokes the loader once per original bulk path.
        let mut loader = |name: &str| -> Result<MaterialPack> {
            poll(gate)?;
            let pack: MaterialPack = loader(name)?.into();
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
        self.friction.validate_surfaces(&self, &catalog)?;
        let bank=if let Some(body_profile)=&profile {
            body_profile.compile_with_preparation(&self,sources,&mut catalog,&mut loader,
                settings.clone(),preparation,cx,gate)?
        } else {
            if self.body.modes > input::MAX_PLATE_BODY_MODES {
                return Err(format!("flat plate bodies admit 1..={} retained body modes",
                    input::MAX_PLATE_BODY_MODES).into());
            }
            let material=resolve(&self.body.pack,&mut catalog,&mut loader,&self,true,false,&input::BendingLoss::KelvinVoigt)?;
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
                    self.bow_law(0)?.clone(),settings.clone(),gate)?,
                None=>ReciprocalStringBank::from_material_plate(cx,sources,self.gas,body,self.bow_law(0)?.clone(),settings.clone(),gate)?,
            }
        };
        let bank = if self.bridge_interfaces.bindings().is_empty() { bank } else {
            bank.with_material_bridge_interfaces(self.bridge_interfaces.bindings().to_vec(), gate)?
        };
        // Select the active prepared owner before either attachment can enlarge
        // its connected linear block. No record means the ordinary owner.
        // A shell body bank is already built on this record (it cannot exist without it); binding it again refuses.
        let bank = match self.prepared_linear {
            Some(budget) if bank.source().prepared_linear_budget().is_none() =>
                bank.with_prepared_linear_propagation(budget, gate)?,
            _ => bank,
        };
        // Omitted bending modes in series with every bridge, before any stop binds (pexfp.299).
        let bank = if self.residual_flexibility { bank.with_residual_flexibility(gate)? } else { bank };
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
        let bank = translational_bridges::bind(bank, &self.translational_bridges, gate)?;
        let bank = if profile.as_ref().is_some_and(body::Profile::has_thermal_cavity) {
            bank // Gas inertia precedes the material/thermal memory prefix.
        } else if let Some(resolution) = self.cylinder_gas_loading {
            bank.with_cylinder_gas_loading(resolution, gate)?
        } else { bank };
        let bank = if let Some((slope, config)) = self.geometric_stretching {
            bank.with_geometric_stretching(slope, config, gate)?
        } else { bank };
        let bank = if let Some((resolution, config)) = &self.longitudinal_waves {
            bank.with_longitudinal_waves_and_bearings(vec![resolution.clone(); self.strings.len()],
                self.longitudinal_loss.clone(), axial_bearings, *config, gate)?
        } else { bank };
        let bank = if profile.as_ref().is_some_and(body::Profile::has_thermal_cavity) {
            bank // The source material prefix was bound before thermal air memory.
        } else {
            let bank = match &profile {
                Some(body) => body.bind_material_memory(bank, gate)?,
                None => bank,
            };
            match &self.bending_memory {
                Some(memory) => memory.bind(bank, gate)?,
                None => bank,
            }
        };
        // The original body basis survives fluid, torsion and material-history
        // binding. Compute its air reaction once, before bows or preload run.
        let (bank, body_radiation) = match self.body_radiation {
            Some(input::BodyRadiation::Rayleigh(options)) => {
                let (bank, report) = bank.with_material_plate_rayleigh(options, gate)?;
                (bank, Some(BodyRadiationReport::Rayleigh(report)))
            },
            Some(input::BodyRadiation::ShellBem(options)) => {
                let (bank, report) = bank.with_material_shell_bem_with_angular_refinement(
                    options,self.shell_bem_angular_refinement,gate)?;
                (bank, Some(BodyRadiationReport::ShellBem(report)))
            },
            Some(input::BodyRadiation::Skin(load)) => {
                let skin = profile.as_ref().and_then(|body| body.acoustic_skin(&self))
                    .ok_or("a computed skin load requires a body profile that declares its closed acoustic skin")?;
                let (bank, report) = bank.with_material_skin_bem(load.options, &skin,
                    load.maximum_anchor_distance_m, load.angular_refinement, gate)?;
                (bank, Some(BodyRadiationReport::ShellBem(report)))
            },
            None => (bank, None),
        };
        let joint=MultiFrictionConfig { max_ports:self.bow_contact_count(),max_sweeps:128,
            max_setup_terms:1_048_576 };
        let mut bank=if let Some(fibers)=&self.bow_fibers {
            let ribbon=fibers.compile(&self,&mut catalog,&mut loader,gate)?;
            bank.with_tensioned_fiber_ribbon(ribbon,fibers.station_fraction,joint,gate)?
        } else if let Some((width,cells))=self.bow_ribbon {
            bank.with_uniform_bow_ribbon(width,cells,joint,gate)?
        } else {
            let additional=contacts[1..].iter().enumerate().map(|(index,location)|Ok(MaterialBowContact {
                location:*location,law:self.bow_law(index+1)?.clone(),friction:settings.friction }))
                .collect::<Result<Vec<_>>>()?;
            bank.with_additional_bows(additional,joint,gate)?
        };
        if let Some((_,cells))=self.bow_ribbon {
            if let Some(binding)=self.friction.binding(&self.schedule.tracks()[0].id)? {
                // The source declaration describes the WHOLE ribbon. Its load
                // and area share the same quadrature, preserving nominal pressure.
                let cell=DryFrictionBinding::new(binding.card().clone(),binding.temperature_kelvin(),
                    binding.nominal_area_m2()/cells as f64)?;
                for port in 0..cells { bank=bank.with_bow_interface(port,cell.clone())?; }
            }
        } else {
            for (port,track) in self.schedule.tracks().iter().enumerate() {
                if let Some(binding)=self.friction.binding(&track.id)? {
                    bank=bank.with_bow_interface(port,binding.clone())?;
                }
            }
        }
        let bank=match &stop_profile {
            Some(stops)=>stops.bind(bank,&self,&mut catalog,&mut loader,gate)?, None=>bank };
        let bank=if let Some(layers)=self.friction.shear_layers() {
            layers.bind(&self,bank,&mut catalog,&mut loader,gate)?
        } else { bank };
        poll(gate)?;
        Ok(ReadyBank { spec:self,catalog,bank,body_profile:profile,
            rotational_bridge_baselines_m:baselines,stop_profile,body_radiation,bow_station_program })
    }

    fn schedule_bank(self, bank:MultiBowStringBank, catalog:Catalog, profile:Option<body::Profile>,
        baselines:Vec<f64>, stop_profile:Option<stops::Profile>, body_radiation:Option<BodyRadiationReport>,
        bow_station_program:Option<BowStationProgram>, gate:&CancelGate) -> Result<Ready> {
        let performance_config=MultiBowPerformanceConfig { samples:self.samples,max_block:self.block,
            max_work:input::MAX_WORK,max_events:input::MAX_CONTROLS };
        let mut performance=if self.bow_ribbon.is_some() {
            bank.into_scheduled_ribbon(self.schedule.clone(),self.schedule.tracks()[0].id.clone(),performance_config)?
        } else {
            let names=self.schedule.tracks().iter().map(|t|t.id.clone()).collect();
            bank.into_scheduled(self.schedule.clone(),names,performance_config)?
        };
        if let Some(forces)=&self.forces { performance=performance.with_material_forces(forces.clone(),gate)?; }
        if let Some(stops)=&stop_profile { performance=stops.bind_motion(performance,gate)?; }
        if let Some(program)=bow_station_program {
            performance=performance.with_bow_station_motion(program,gate)?;
        }
        if let Some(limits)=self.bow_fibers.as_ref().and_then(|fibers|fibers.motion) {
            let GestureValue::Bow {velocity_m_per_s,..}=self.schedule.tracks()[0].initial else {
                return Err("fiber motion requires an explicit initial bow speed".into());
            };
            performance=performance.with_fiber_support_motion(velocity_m_per_s,limits,gate)?;
        }
        if let Some(program)=self.bow_fibers.as_ref().and_then(|fiber|fiber.support_program.as_ref()) {
            performance=performance.with_fiber_support_program(program.clone(),gate)?;
        }
        if let Some(program)=self.bow_fibers.as_ref().and_then(|fiber|fiber.normal_program.as_ref()) {
            performance=performance.with_fiber_normal_support_program(program.clone(),gate)?;
        }
        if let Some(observer)=&self.plate_observer {
            performance=performance.with_material_plate_observer(observer.clone(),gate)?;
        }
        if let Some(observer)=&self.shell_observer {
            performance=performance.with_material_shell_observer_with_quadrature(
                observer.clone(),self.shell_receiver_quadrature,gate)?;
        }
        if let Some(observer)=&self.one_way_observer {
            let skin=profile.as_ref().and_then(|body|body.acoustic_skin(&self))
                .ok_or("one_way_observer requires a body profile that declares its closed acoustic skin")?;
            let store=profile.as_ref().and_then(|body|body.observer_store());
            performance=performance.with_one_way_shell_observer(observer.clone(),&skin,store,gate)?;
        }
        performance.validate_sample_rate(self.rate)?;performance.validate_sample_count(self.samples)?;poll(gate)?;
        Ok(Ready { spec:self,catalog,performance,body_profile:profile,rotational_bridge_baselines_m:baselines,stop_profile,body_radiation })
    }
}
impl Ready {
    fn write<W:Write+Seek>(&mut self,output:W,gate:&CancelGate)->Result<(W,WavStreamSummary)> {
        self.write_observed(output,None,false,gate)
    }
    /// `write` one sample per callback, with each accepted sample's bow ports in `trace`: the physical sample, port,
    /// regime, slip velocity [m/s] (the step's slip distance times the rate), traction and normal force [N]. The WAV
    /// bytes are `write`'s: the solve is causal and the PCM writer encodes each sample alone.
    fn write_traced<W:Write+Seek>(&mut self,output:W,trace:&mut impl Write,gate:&CancelGate)->Result<(W,WavStreamSummary)> {
        self.write_observed(output,Some(trace),false,gate)
    }
    fn write_observed<W:Write+Seek>(&mut self,output:W,trace:Option<&mut dyn Write>,
        decimate:bool,gate:&CancelGate)->Result<(W,WavStreamSummary)> {
        match trace {
            Some(trace)=>write_traced_pressure(&mut self.performance,output,trace,&self.spec,
                decimate,gate,|performance|performance.bank().system().last_frame()),
            None=>write_pressure(&mut self.performance,output,&self.spec,decimate,gate),
        }
    }
}

/// Observe the original accepted mechanical steps before any output resampling.
/// The callback may write an offline trace; failure prevents a completed WAV and
/// poisons this observation adapter without pretending to roll back mechanics.
struct ObservedPressure<'a,S,F> {
    source:&'a mut S,
    observe:F,
    failure:Option<Box<dyn Error>>,
}
impl<S:PressureRenderer,F:FnMut(&S)->Result<()>> PressureRenderer for ObservedPressure<'_,S,F> {
    fn samples_rendered(&self)->u64 { self.source.samples_rendered() }
    fn max_block_len(&self)->usize { 1 }
    fn validate_sample_rate(&self,rate:u32)->std::result::Result<(),fs_couple::render::RenderError> {
        if self.failure.is_some() { return Err(fs_couple::render::RenderError::Poisoned); }
        self.source.validate_sample_rate(rate)
    }
    fn validate_sample_count(&self,count:u64)->std::result::Result<(),fs_couple::render::RenderError> {
        if self.failure.is_some() { return Err(fs_couple::render::RenderError::Poisoned); }
        self.source.validate_sample_count(count)
    }
    fn block(&mut self,output:&mut [f64])->std::result::Result<(),fs_couple::render::RenderError> {
        self.validate_sample_count(output.len() as u64)?;
        if output.len()!=1 {
            return Err(fs_couple::render::RenderError::Sizing { what:"accepted-step observation requires one mechanical sample" });
        }
        self.source.block(output)?;
        if let Err(error)=(self.observe)(self.source) {
            self.failure=Some(error);
            return Err(fs_couple::render::RenderError::Control { what:"accepted mechanical-step observation failed" });
        }
        Ok(())
    }
}

/// All native string writers use this same pressure clock, encoder and scaling.
/// Thermal feedback and tracing are source adapters, never output-rate physics.
fn write_pressure<S:PressureRenderer,W:Write+Seek>(source:&mut S,output:W,spec:&input::Spec,
    decimate:bool,gate:&CancelGate)->Result<(W,WavStreamSummary)> {
    let (mut renderer,samples,rate):(Box<dyn PressureRenderer+'_>,u64,u32)=if decimate {
        if spec.rate<=48_000 { return Err("--decimate requires mechanics above 48000 Hz".into()); }
        let renderer=DecimatedRenderer::new(source,spec.rate,48_000,spec.block)?;
        let samples=renderer.output_samples_for(spec.samples)?;
        (Box::new(renderer),samples,48_000)
    } else { (Box::new(source),spec.samples,spec.rate) };
    let mut stream=Pcm16WavStream::new(output,rate,spec.full_scale,spec.block)?;
    let mut scratch=vec![0.0;spec.block.min(renderer.max_block_len())];
    let progress=render_pressure_pcm16(&mut renderer,&mut stream,gate,&mut scratch,samples)?;
    let result=stream.finish()?;
    if !matches!(progress,ScheduledWavProgress::Completed { .. }) {
        return Err(format!("cancelled: only {} output samples finalized",result.1.samples).into());
    }
    Ok(result)
}

fn write_traced_pressure<S:PressureRenderer,W:Write+Seek>(source:&mut S,output:W,
    trace:&mut dyn Write,spec:&input::Spec,decimate:bool,gate:&CancelGate,
    frame:impl for<'s> Fn(&'s S)->Option<&'s fs_couple::render::schedule::force::coupled::contact::friction::multiple::MultiFrictionFrame>)
    ->Result<(W,WavStreamSummary)> {
    writeln!(trace,"sample\tport\tregime\tslip_velocity_m_s\ttraction_n\tnormal_force_n")?;
    let mut observed=ObservedPressure { source,failure:None,observe:|source:&S|->Result<()> {
        let frame=frame(source).ok_or("a rendered sample left no accepted bow frame")?;
        for (port,p) in frame.ports.iter().enumerate() {
            writeln!(trace,"{}\t{port}\t{:?}\t{:e}\t{:e}\t{:e}",frame.sample,p.regime,
                p.slip_distance_m*f64::from(spec.rate),p.traction_n,p.normal_force_n)?;
        }
        Ok(())
    } };
    let result=write_pressure(&mut observed,output,spec,decimate,gate);
    if let Some(error)=observed.failure { return Err(error); }
    result
}
fn read_bounded(path:&Path,cap:usize)->Result<Vec<u8>> {
    let mut bytes=Vec::new();File::open(path)?.take(cap as u64+1).read_to_end(&mut bytes)?;
    if bytes.len()>cap {return Err(format!("{} exceeds {cap} bytes",path.display()).into());} Ok(bytes)
}
/// `material_strings assemble-shell TABLE OPTIONS NEW_PROFILE`: a `frankensim.shell-assembly.v1` table written as a v4
/// body profile with the choices of a `frankensim-shell-assembly-options-v1` file (BODY_PROFILES.md). The profile is
/// admitted before it is written. Printed: the body length and width [m] and reference thickness the performance's body
/// record takes, how far each bridge point is from its node, and the table's cavity row for the cavity record.
fn assemble_shell(args:&[std::ffi::OsString])->Result<()> {
    let [table,options,out]=args else {
        return Err("usage: material_strings assemble-shell TABLE OPTIONS NEW_PROFILE".into());
    };
    let options=body::ShellTableOptions::parse(std::str::from_utf8(&read_bounded(Path::new(options),body::MAX_BYTES)?)?)?;
    let converted=options.convert(std::str::from_utf8(&read_bounded(Path::new(table),body::MAX_SHELL_BYTES)?)?)?;
    body::Profile::parse(converted.profile.as_bytes())?;
    OpenOptions::new().write(true).create_new(true).open(out)?.write_all(converted.profile.as_bytes())?;
    println!("body-length-m {}\nbody-width-m {}\nreference-thickness-m {}",converted.length_m,converted.width_m,
        options.reference_thickness_m());
    for (k,snap) in converted.bridge_snaps_m.iter().enumerate() { println!("bridge-snap-m {k} {snap}"); }
    if let Some([volume,radius,length])=converted.cavity {
        println!("cavity-volume-m3 {volume}\ncavity-neck-radius-m {radius}\ncavity-neck-length-m {length}");
    }
    Ok(())
}
/// Take `NAME PATH` (`--observer-store DIR`, `--bow-trace FILE`) out of the option pairs; the rest pass on.
fn path_option<'a>(args:&'a [std::ffi::OsString],name:&str)->Result<(Option<&'a Path>,Vec<std::ffi::OsString>)> {
    let (mut path,mut rest)=(None,Vec::new());
    for pair in args.chunks(2) {
        if pair[0]!=name { rest.extend_from_slice(pair); continue; }
        if path.is_some() { return Err(format!("{name} cannot be repeated").into()); }
        path=Some(Path::new(pair.get(1).ok_or_else(|| format!("{name} needs a path"))?));
    }
    Ok((path,rest))
}
fn main()->Result<()> {
    let raw=std::env::args_os().collect::<Vec<_>>();
    let (workers,raw)=workers_option(&raw)?;
    fs_bem::helmholtz::with_dense_workers(workers,|| run(&raw))
}
/// Take `--workers N` (threads for the dense BEM air loads; 1 when absent) out
/// of the arguments. A scheduling budget only: every load is bit-identical.
fn workers_option(raw:&[std::ffi::OsString])->Result<(usize,Vec<std::ffi::OsString>)> {
    let (mut workers,mut rest,mut args)=(None,Vec::with_capacity(raw.len()),raw.iter());
    while let Some(arg)=args.next() {
        if arg!="--workers" { rest.push(arg.clone()); continue; }
        if workers.is_some() { return Err("--workers cannot be repeated".into()); }
        workers=Some(args.next().and_then(|n|n.to_str()).and_then(|n|n.parse::<usize>().ok())
            .filter(|n|(1..=256).contains(n)).ok_or("--workers needs a thread count from 1 to 256")?);
    }
    Ok((workers.unwrap_or(1),rest))
}
fn run(raw:&[std::ffi::OsString])->Result<()> {
    if raw.get(1).is_some_and(|a| a=="assemble-shell") { return assemble_shell(&raw[2..]); }
    let (args, gas_options)=thermochemical::Options::extract_os(raw)?;
    let (args, response_options)=body_response::Options::extract(&args)?;
    let (args, coating_options)=coating_strain::Options::extract(&args)?;
    let response_only=response_options.as_ref().is_some_and(body_response::Options::analysis_only);
    let loaded_response_only=response_only
        && response_options.as_ref().is_some_and(body_response::Options::loaded);
    let positional_count=if response_only {2} else {3};
    let (args, source_options)=if response_only { native_input::SourceOptions::extract_after(&args,2)? }
        else { native_input::SourceOptions::extract(&args)? };
    if args.len()<positional_count {
        return Err("usage: material_strings PERFORMANCE NEW_OUTPUT.wav [--body-response|--loaded-body-response NEW.tsv MIN_HZ MAX_HZ POINTS BRIDGE_SELECTION] [options], or material_strings PERFORMANCE --body-response-only|--loaded-body-response-only NEW.tsv MIN_HZ MAX_HZ POINTS BRIDGE_SELECTION --body-profile PROFILE [source options]; loaded selections: BRIDGE_INDEX, all, or coherent with repeated --body-response-force BRIDGE REAL_N IMAG_N".into());
    }
    if gas_options.record.is_some() {
        return Err("material-string gas records belong in the performance file".into());
    }
    let path=Path::new(&args[1]);let parent=path.parent().unwrap_or(Path::new("."));
    let sources=source_options.load(parent)?;
    let (store_dir,options)=path_option(&args[positional_count..],"--observer-store")?;
    if sources.is_restricted() && store_dir.is_some() {
        return Err("--self-contained refuses --observer-store; observer setup must be computed from attached physical sources".into());
    }
    let (trace_path,options)=path_option(&options,"--bow-trace")?;
    let (thermal_path,options)=path_option(&options,"--contact-thermal")?;
    if response_only && (store_dir.is_some() || trace_path.is_some() || thermal_path.is_some() || source_options.decimate) {
        return Err("response-only commands accept body/source options; playback observer stores, traces, contact heat and decimation need a WAV render".into());
    }
    let thermal_profile=thermal_path.map(|p|contact_thermal::Profile::parse(
        &sources.read_option(p,contact_thermal::MAX_BYTES)?)).transpose()?;
    let (body_path,stop_path)=stops::profile_paths(&options)?;
    if response_only && stop_path.is_some() {
        return Err("response-only commands do not admit --stop-profile; no nonlinear contact linearization is selected".into());
    }
    let mut profile=body_path.map(|p|body::Profile::parse(&sources.read_option(p,body::MAX_SHELL_BYTES)?)).transpose()?;
    if response_options.is_some() && !profile.as_ref().is_some_and(body::Profile::has_shell_material_basis) {
        return Err("--body-response requires --body-profile with curved or assembled shell mechanics".into());
    }
    if let Some(options)=coating_options {
        options.validate_context(response_only,
            profile.as_ref().is_some_and(body::Profile::has_shell_material_basis))?;
    }
    if let Some(dir)=store_dir {
        profile.as_mut().ok_or("--observer-store requires --body-profile")?.keep_observer_setups(dir)?;
    }
    let mut stop_profile=stop_path.map(|p|stops::Profile::parse(&sources.read_option(p,stops::MAX_BYTES)?)).transpose()?;
    let gas_sources=if sources.is_restricted() {
        gas_options.load_with(|p|sources.read_option(p,input::MAX_PACK_BYTES).map_err(|e|e.to_string()))?
    } else { gas_options.load()? };
    let gas_bindings=gas_sources.bindings();
    let gas_registry=gas_sources.registry(&gas_bindings);
    let mut spec=input::Spec::parse_with_gas_sources(
        &read_bounded(path,input::MAX_INPUT_BYTES)?,gas_registry.as_ref())?;
    if response_options.as_ref().is_some_and(body_response::Options::loaded)
        && spec.plate_observer.is_some() {
        return Err("loaded shell response requires its shell or acoustic-skin receiver; a planar microphone is incompatible".into());
    }
    if let Some(profile)=&thermal_profile { profile.configure(&mut spec); }
    let gate=std::sync::Arc::new(CancelGate::new_clock_free());
    let pool=fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let mut remaining=input::MAX_TOTAL_PACK_BYTES;
    if let Some(stops) = &stop_profile { stops.validate_performance(&spec,&gate)?; }
    if let Some(body) = &mut profile {
        body.load_interfaces(&spec,&gate,|name| {
            let bytes=sources.read_relative(name,input::MAX_PACK_BYTES.min(remaining))?;
            remaining=remaining.checked_sub(bytes.len()).ok_or("combined source pack byte budget exceeded")?;
            Ok(NormalizedInterfacePack::from_bytes(&bytes)?)
        })?;
    }
    if response_only && !loaded_response_only {
        let mut catalog=Catalog::new();
        let projected=profile.as_ref().ok_or("--body-response-only requires --body-profile")?
            .project_shell(&spec,&mut catalog,&mut |name| {
                let bytes=sources.read_relative(name,input::MAX_PACK_BYTES.min(remaining))?;
                remaining=remaining.checked_sub(bytes.len()).ok_or("material pack byte budget exceeded")?;
                MaterialPack::from_bytes(&bytes)
            },&gate)?;
        for (path,source) in &catalog {
            eprintln!("body_source={path:?}; pack={:?}; hash={}; artifact={}",source.pack.pack_id(),
                source.pack.content_hash(),source.pack.source_artifact());
        }
        return response_options.as_ref().expect("response-only option checked")
            .write(projected.shell(),projected.assembly(),&gate);
    }
    spec.friction.load(&gate,|name| {
        let bytes=sources.read_relative(name,input::MAX_PACK_BYTES.min(remaining))?;
        remaining=remaining.checked_sub(bytes.len()).ok_or("combined source pack byte budget exceeded")?;
        Ok(NormalizedInterfacePack::from_bytes(&bytes)?)
    })?;
    spec.bridge_interfaces.load(spec.temperature, spec.omega, &gate, |name| {
        let bytes=sources.read_relative(name,input::MAX_PACK_BYTES.min(remaining))?;
        remaining=remaining.checked_sub(bytes.len()).ok_or("combined source pack byte budget exceeded")?;
        Ok(NormalizedInterfacePack::from_bytes(&bytes)?)
    })?;
    spec.axial_bearings.load(spec.temperature, spec.omega, &gate, |name| {
        let bytes=sources.read_relative(name,input::MAX_PACK_BYTES.min(remaining))?;
        remaining=remaining.checked_sub(bytes.len()).ok_or("combined source pack byte budget exceeded")?;
        Ok(NormalizedInterfacePack::from_bytes(&bytes)?)
    })?;
    if let Some(stops) = &mut stop_profile {
        stops.load_sources(&gate, |name| {
            let bytes=sources.read_relative(name,input::MAX_PACK_BYTES.min(remaining))?;
            remaining=remaining.checked_sub(bytes.len()).ok_or("combined source pack byte budget exceeded")?;
            Ok(bytes)
        })?;
    }
    let bound=pool.scope(|arena| {
        let cx=Cx::new(&gate,arena,fs_exec::StreamKey {seed:0,kernel_id:0,tile:0,iteration:0},
            fs_exec::Budget::INFINITE,fs_exec::ExecMode::Deterministic);
        spec.compile_bank_with_profiles(&cx,&gate,profile,stop_profile,loaded_response_only,|name| {
            let bytes=sources.read_relative(name,input::MAX_PACK_BYTES.min(remaining))?;
            remaining=remaining.checked_sub(bytes.len()).ok_or("material pack byte budget exceeded")?;
            MaterialPack::from_bytes(&bytes)
        })
    })?;
    if loaded_response_only {
        let work=bound.spec.prepared_linear.ok_or(
            "loaded body-response-only lost its explicit work declaration")?.max_multiply_terms;
        for (path,source) in &bound.catalog {
            eprintln!("coupled_source={path:?}; pack={:?}; hash={}; artifact={}",source.pack.pack_id(),
                source.pack.content_hash(),source.pack.source_artifact());
        }
        eprintln!("loaded_response_preparation=harmonic_only; time_propagator=absent; state_clock=unadvanced");
        return response_options.as_ref().expect("loaded response-only option checked")
            .write_loaded_for_spec(&bound.bank,&bound.spec,bound.body_profile.as_ref(),work,&gate);
    }
    let mut ready=bound.into_performance(&gate)?;
    if let Some(options)=coating_options {
        ready.performance=ready.performance.with_material_shell_coating_monitor(
            options.budget(),std::sync::Arc::clone(&gate))?;
    }
    let contact_heat=thermal_profile.as_ref().map(|profile|profile.prepare(&ready,&gate)).transpose()?;
    if let Some(heat)=&contact_heat { heat.describe(); }
    for (path,source) in &ready.catalog {
        eprintln!("source={path:?}; pack={:?}; hash={}; artifact={}",source.pack.pack_id(),
            source.pack.content_hash(),source.pack.source_artifact());
    }
    let table=input::tension_loss_rows()?;
    for (string,wire) in ready.spec.strings.iter().enumerate() {
        let others=|id:Option<&str>| table.iter().filter(|row|Some(row.id.as_str())!=id)
            .map(|row|format!("{} {}",row.id,row.eta_f)).collect::<Vec<_>>().join(", ");
        match &wire.tension_loss {
            Some(row)=>eprintln!("tension_loss_string={string}; row={}; eta_F={}; source={:?}; other_rows=[{}]",
                row.id,row.eta_f,row.source,others(Some(&row.id))),
            None=>eprintln!("tension_loss_string={string}; eta_F=stated unknown, entered as 0 (no tension loss: air and the pack's bending law only); rows=[{}]",
                others(None)),
        }
    }
    eprintln!("bow_friction={:?}",ready.spec.friction);
    if let Some(forces)=ready.performance.material_force_program() {
        eprintln!("physical_force_ports={:?}; initial_force_N={:?}; initialization={:?}; authored_events={}; preload_energy_J={:?}; shared mechanical clock and coupled solver",
            forces.ports,forces.initial_force_n,forces.initialization,forces.events.len(),
            ready.performance.preload_energy_j());
    }
    for port in 0..ready.spec.bow_contact_count() {
        if let Some(layer)=ready.performance.bank().bow_shear_layer(port) {
            eprintln!("bow_shear_layer_port={port}; source={}; area_m2={}; thickness_m={}; initial_T_K={}; stiffness_N_m={}; yield_N={}; plastic_resistance_N_s_m={}; max_elastic_shear={}",
                layer.resolved().identity(),layer.area_m2(),layer.thickness_m(),layer.temperature_kelvin(),
                layer.stiffness_n_m(),layer.yield_force_n(),layer.plastic_resistance_n_s_m(),layer.maximum_elastic_shear_strain());
        }
    }
    if !ready.spec.bridge_interfaces.bindings().is_empty() {
        eprintln!("bridge_interfaces={:?}",ready.spec.bridge_interfaces);
    }
    if ready.spec.axial_bearings.len() != 0 {
        eprintln!("axial_bearings={:?}",ready.spec.axial_bearings);
    }
    if let Some(ribbon)=ready.performance.bank().source().ribbon() { eprintln!("bow_ribbon={ribbon:?}"); }
    if let Some(fibers)=ready.performance.bank().source().fiber_ribbon() {
        eprintln!("bow_fibers={:?}; resolution={:?}; mass_kg={}; station={}; support_motion={:?}; rigid-end material fibers, not validated thermal rosin",
            fibers.model().geometry(),fibers.model().resolution(),fibers.model().mass_kg(),fibers.station_fraction(),
            ready.performance.bank().fiber_support_motion_limits());
        if let Some(strain)=fibers.model().strain() {
            eprintln!("bow_fiber_strain={:?}; joint_solve={:?}; moderate-slope fixed-end material strain",
                strain.limits(),strain.config());
        }
    }
    if let Some(motion)=ready.performance.fiber_support_program() {
        eprintln!("bow_fiber_support_motion={motion:?}; material-derived frame inertia and laboratory kinetic correction");
    }
    if let Some(motion)=ready.performance.fiber_normal_support_program() {
        eprintln!("bow_fiber_normal_support_motion={motion:?}; normal material inertia and solved contact load on the same physical clock");
    }
    if let Some(profile)=&ready.body_profile { eprintln!("body_profile={profile:?}"); }
    if let Some(memory)=ready.performance.bank().source().bank().material_shell_coating_relaxation() {
        let band=memory.source_band_rad_s();
        eprintln!("coating_memory=prony-continuation; body_modes={}; histories={}; source_fit_band_Hz=[{},{}]; continuation=passive_model_outside_measured_band; initial_history=relaxed; full_modal_cross_terms=true",
            memory.mode_count(),memory.coupling_rows_per_s().len(),
            band.0/core::f64::consts::TAU,band.1/core::f64::consts::TAU);
        for (fit,certificate) in memory.fit_certificates().iter().enumerate() {
            eprintln!("coating_fit={fit}; sampled_storage_relative_error={}; sampled_loss_relative_error={}; verification_points={}; errors=numerical_fit_not_measurement_uncertainty",
                certificate.sup_rel_storage_err,certificate.sup_rel_loss_err,certificate.verification_points);
        }
    }
    if let Some(report)=&ready.body_radiation { eprintln!("computed_body_radiation={report}"); }
    // A skin load's first-order receipt at each loaded coordinate's own frequency,
    // in the one-way receipt's terms: zeta = Re Z/(2 w), shift = Im Z/(2 w).
    if let Some(binding)=ready.performance.bank().source().bank().body_radiation()
        .filter(|binding|binding.material_shell_bem().is_some_and(|r|r.skin_identity.is_some())) {
        let components=ready.performance.bank().system().components();
        let body=ready.spec.strings.len();
        let structural=components[body].modes().len();
        let ports=binding.model().ports;
        for (port,&coordinate) in binding.source_body_coordinates().iter().enumerate() {
            let w=if coordinate<structural { components[body].modes()[coordinate].angular_frequency_rad_s }
                else { components[body+1].modes()[coordinate-structural].angular_frequency_rad_s };
            let z=binding.model().impedance(w)?[port*ports+port];
            eprintln!("two_way_load_port={port}; coordinate={coordinate}; f_Hz={}; radiation_zeta={}; added_mass_shift={}",
                w/core::f64::consts::TAU,z.re/(2.0*w),0.5*z.im/w);
        }
    }
    if let Some(report)=ready.performance.material_plate_observer() {
        eprintln!("plate_observer={:?}; specimen={}; basis={}; body_modes={}; surface_cells={}; maximum_Hz={}; minimum_disc_ppw={}; delay_support_samples={:?}; numerical_delay_samples={}; replay_terms={}; gas_law={:?}; distinct_paths={}; gas_kernel_terms={}; same signed loaded plate surface; finite-distance pressure replaces compact pressure; infinite baffle, no room or closed-shell claim",
            report.options.position_m, report.specimen_identity, report.basis_identity,
            report.body_modes, report.surface_cells, report.options.maximum_frequency_hz,
            report.minimum_equivalent_disc_ppw, report.delay_samples, report.filter_delay_samples,
            report.replay_terms, report.gas.law, report.gas.filters.len(), report.gas.kernel_terms);
    }
    if let Some(report)=ready.performance.material_shell_observer() {
        let transfer=report.transfer();
        eprintln!("shell_observer={:?}; basis={}; loading={}; body_modes={}; surface_panels={}; pass_Hz={}; taper_stop_Hz={}; geometric_delays_s={:?}; reconstruction_latency_samples={}; sweep_terms={}; replay_terms={}; closed-shell finite-distance pressure replaces compact pressure; homogeneous exterior without molecular absorption; declared latency and finite window retained",
            report.options.receiver,report.basis_identity,report.loading_identity,
            report.body_modes,report.surface_panels,transfer.pass_band_upper_hz,transfer.taper_stop_hz,
            transfer.geometric_delays_s,transfer.filter_latency_samples,
            report.sweep_terms,report.replay_terms);
    }
    if let Some(report)=ready.performance.one_way_shell_observer() {
        let transfer=report.transfer();
        let feedback=match report.two_way_load {
            Some(_) => "the skin's own computed air load reacts on the two-way body modes and the mouth; the rest stay one-way",
            None => "no air reaction enters the mechanics",
        };
        eprintln!("one_way_observer={:?}; body_modes={}; apertures={}; surface_panels={}; two_way_load={:?}; two_way_modes={}; largest_anchor_distance_m={}; plain_cbie_limit_rad_per_m={}; pass_Hz={}; taper_stop_Hz={}; geometric_delays_s={:?}; reconstruction_latency_samples={}; boundary_terms={}; replay_terms={}; rule={:?}; setup_key={}; setup_loaded={}; setup_entry={:?}; skin finite-distance pressure replaces compact pressure; {feedback}; homogeneous exterior without molecular absorption",
            report.options.receiver,report.body_modes,report.apertures,report.surface_panels,
            report.two_way_load.map(|key|key.to_hex()),report.two_way_modes,
            report.largest_anchor_distance_m,report.plain_cbie_limit_rad_per_m,transfer.pass_band_upper_hz,
            transfer.taper_stop_hz,transfer.geometric_delays_s,transfer.filter_latency_samples,
            report.boundary_terms,report.replay_terms,report.options.rule,report.setup_key.to_hex(),
            report.setup_loaded,ready.body_profile.as_ref().and_then(|body|body.observer_store_entry(&report.setup_key)));
        for (mode,bound) in report.modes.iter().enumerate() {
            eprintln!("one_way_mode={mode}; f_Hz={}; wood_zeta={}; radiation_zeta={}; added_mass_shift={}; half_power_bandwidth_Hz={}; nearest_spacing_Hz={}; largest_mixing={}; mixing_mode={:?}; one_way={}",
                bound.frequency_hz,bound.wood_zeta,bound.radiation_zeta,bound.added_mass_shift,
                bound.half_power_bandwidth_hz,bound.nearest_spacing_hz,bound.largest_mixing,
                bound.largest_mixing_mode,bound.one_way);
        }
        eprintln!("one_way_fed_back_modes={:?}",report.fed_back_modes());
    }
    if let Some(profile)=&ready.stop_profile { eprintln!("stop_profile_bytes={}; stops={}",
        profile.bytes.len(),ready.performance.bank().source().stops().len()); }
    for (joint, baseline) in ready.spec.rotational_bridges.iter().zip(&ready.rotational_bridge_baselines_m) {
        eprintln!("rotational_bridge={joint:?}; body_footprint_separation_m={baseline}");
    }
    for joint in &ready.spec.translational_bridges {
        eprintln!("translational_bridge={joint:?}; body_map=original_shell_global_displacement");
    }
    if let Some(budget)=ready.performance.bank().source().bank().prepared_linear_budget() {
        eprintln!("prepared_linear={budget:?}");
    }
    if ready.spec.residual_flexibility {
        let bank=ready.performance.bank();
        for (index,(member,card)) in bank.source().bank().strings().iter().zip(bank.cards()).enumerate() {
            for connection in &member.connections {
                let condensed=residual::condense_bridge(&member.source,card,connection)?;
                eprintln!("residual_bridge string={index}; station={}; G_m_N={:e}; k_eff_N_m={:e}; c_eff_N_s_m={:e}",
                    connection.station_fraction,condensed.residual_m_n,condensed.stiffness_n_m,condensed.damping_n_s_m);
            }
        }
        eprintln!("residual_stop_compliance_m_N={:?}; opposed_stop_pairs={:?}",
            bank.system().obstacle_residual_compliance(),bank.system().obstacle_pairs());
    }
    eprintln!("linear_frequency_upper_bound_rad_s={:e}",ready.performance.bank().system().linear_frequency_upper_bound_rad_s());
    if let Some(loading)=ready.performance.bank().source().bank().cylinder_gas_loading() {
        let additional = ready.body_profile.as_ref()
            .map_or(0, body::Profile::additional_body_connections);
        let air_modes = ready.body_profile.as_ref()
            .map_or(0, body::Profile::additional_body_modes);
        let (modes,links,terms)=ready.spec
            .network_dimensions_with_body_modes(ready.body_profile.is_some(),additional,air_modes)?;
        eprintln!("cylinder_gas_loading={:?}; complete_modes={modes}; connection_allowance={links}; setup_terms_allowance={terms}",
            loading.resolution());
    }
    eprintln!("input_bytes={}; gas={:?}; preload_J={:?}",ready.spec.bytes.len(),ready.spec.gas,ready.performance.preload_energy_j());
    let gas_source_json=thermochemical::provenance(core::slice::from_ref(&ready.spec.gas_record));
    if !gas_source_json.is_empty() { eprintln!("gas_source={{\"scope\":\"frozen material-string atmosphere\"{gas_source_json}}}"); }
    if let Some(interior) = &ready.spec.cavity_gas_record {
        let source = thermochemical::provenance(core::slice::from_ref(interior));
        eprintln!("cavity_gas={:?}; cavity_gas_source={{\"scope\":\"frozen cavity interior\"{source}}}", interior.state());
    }
    // Complete source, geometry, frequency, schedule and preload admission is
    // finished before creating output. Never overwrite any existing path.
    if source_options.decimate {
        if ready.spec.rate<=48_000 { return Err("--decimate requires mechanics above 48000 Hz".into()); }
        let renderer=DecimatedRenderer::new(&mut ready.performance,ready.spec.rate,48_000,ready.spec.block)?;
        let info=renderer.info();
        renderer.output_samples_for(ready.spec.samples)?;
        eprintln!("pressure_observation={info:?}; mechanical state and gesture clock unchanged");
    }
    if let Some(response)=response_options {
        if response.loaded() {
            let work=ready.spec.prepared_linear.map_or(input::MAX_WORK,|budget|budget.max_multiply_terms);
            response.write_loaded_for_spec(ready.performance.bank(),&ready.spec,
                ready.body_profile.as_ref(),work,&gate)?;
        } else {
            let source=ready.performance.bank().source().bank();
            response.write(source.material_shell_image(),source.material_shell_assembly_image(),&gate)?;
        }
    }
    let output=OpenOptions::new().write(true).create_new(true).open(&args[2])?;
    let mut trace=trace_path.map(|path|OpenOptions::new().write(true).create_new(true).open(path)
        .map(std::io::BufWriter::new)).transpose()?;
    let trace_writer=trace.as_mut().map(|writer|writer as &mut dyn Write);
    let written=match &contact_heat {
        Some(heat)=>heat.write_observed(&mut ready,output,trace_writer,source_options.decimate,&gate),
        None=>ready.write_observed(output,trace_writer,source_options.decimate,&gate),
    };
    coating_strain::describe(&ready.performance,&mut std::io::stderr().lock())?;
    if let Some(trace)=trace.as_mut() { trace.flush()?; }
    let (_,summary)=written.map_err(|e|format!("WAV not completed; keep/discard its partial prefix explicitly: {e}"))?;
    if let Some(frame)=ready.performance.fiber_support_frame() {
        eprintln!("bow_fiber_support_frame={frame:?}");
    }
    eprintln!("samples={}; rate_Hz={}; clips={}; full_scale_Pa={}",summary.samples,summary.sample_rate_hz,
        summary.clipped_samples,summary.full_scale_pa);
    // Steps whose stop missed its tolerance at the actual endpoint and passed it after a headroom re-solve.
    eprintln!("obstacle_resolves={}",ready.performance.bank().system().obstacle_resolves());
    Ok(())
}

#[cfg(test)]
#[path = "material_strings/tests.rs"]
mod tests;

