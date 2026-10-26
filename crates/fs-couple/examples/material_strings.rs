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
#[path = "material_strings/bow_interfaces.rs"]
mod bow_interfaces;
#[path = "material_strings/shear_layers.rs"]
mod shear_layers;
#[path = "material_strings/bridge_interfaces.rs"]
mod bridge_interfaces;
#[path = "material_strings/bending_memory.rs"]
mod bending_memory;
#[path = "material_strings/torsion_memory.rs"]
mod torsion_memory;
#[path = "material_strings/fibers.rs"]
mod fibers;
// Native entry points share one bounded source loader and phase policy.
#[path = "../src/bin/music_render/thermochemical.rs"]
mod thermochemical;
#[path = "../src/bin/music_render/json.rs"]
mod json;
use json::json_string;
use std::{collections::BTreeMap, error::Error, fs::{File,OpenOptions}, io::{Read,Seek,Write}, path::Path};
use fs_couple::bowed_string::runtime::reciprocal::bank::{
    MaterialStringSource,ReciprocalStringBank,StringBankOptions,StringBowLocation};
use fs_couple::bowed_string::runtime::reciprocal::bank::plate::{SharedMaterialPlateBody,SharedPlateBridge};
use fs_couple::bowed_string::runtime::reciprocal::bank::torsion::{StringBodyRotationConnection,StringBridgeLocation,
    SHEAR_VISCOSITY_REFERENCE_OMEGA as TORSION_REFERENCE,SHEAR_VISCOSITY_REFERENCE_OMEGA_DIMS as TORSION_REFERENCE_DIMS};
use fs_couple::bowed_string::runtime::reciprocal::bank::radiation::rayleigh::PlateRayleighReport;
use fs_couple::bowed_string::runtime::reciprocal::bank::radiation::bem::ShellBemReport;
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::MaterialBowContact;
use fs_couple::bowed_string::runtime::reciprocal::bank::nonlinear::longitudinal::{LongitudinalStringResolution,LongitudinalMaterialLoss};
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
        tension_loss_factor:None }
}
// The retained pack includes normalization/uncertainty/rights metadata. Its ID
// labels source identity, not a new assertion about chemistry or measured quality.
struct CatalogEntry { pack: NormalizedPack, card: MaterialCard }
type Catalog = BTreeMap<String,CatalogEntry>;
fn resolve(name:&str, catalog:&mut Catalog, loader:&mut impl FnMut(&str)->Result<NormalizedPack>,
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
fn load_pack(name:&str, catalog:&mut Catalog, loader:&mut impl FnMut(&str)->Result<NormalizedPack>) -> Result<()> {
    if !catalog.contains_key(name) {
        let pack=loader(name)?;
        let card=MaterialCard::assemble(MaterialStateId { chemistry:format!("source-pack:{}",pack.pack_id()),
            phase:"declared solid model".into(),process:"frozen source-pack realization".into(),revision:0 },
            pack.claims().clone(),vec![])?;
        catalog.insert(name.into(),CatalogEntry { pack,card });
    }
    Ok(())
}
fn resolve_properties(name:&str, catalog:&mut Catalog,
    loader:&mut impl FnMut(&str)->Result<NormalizedPack>, spec:&input::Spec,
    requirements:&[Requirement]) -> Result<ResolvedMaterialStatePoint>
{
    load_pack(name,catalog,loader)?;
    // This compact input profile explicitly declares bare SI T and omega axes.
    // Incompatible typed/context axes refuse in the original resolver; they are
    // never removed, guessed, coerced or supplied from an unrelated source.
    let query=QueryPoint::new().with("T",spec.temperature)?.with("omega",spec.omega)?;
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
impl input::Spec {
    fn compile(self,cx:&Cx<'_>,gate:&CancelGate,loader:impl FnMut(&str)->Result<NormalizedPack>)->Result<Ready> {
        self.compile_with_body(cx,gate,None,loader)
    }
    fn compile_with_body(self,cx:&Cx<'_>,gate:&CancelGate,profile:Option<body::Profile>,
        loader:impl FnMut(&str)->Result<NormalizedPack>)->Result<Ready> {
        self.compile_with_profiles(cx,gate,profile,None,loader)
    }
    fn compile_with_profiles(mut self,cx:&Cx<'_>,gate:&CancelGate,profile:Option<body::Profile>,
        stop_profile:Option<stops::Profile>,
        mut loader:impl FnMut(&str)->Result<NormalizedPack>)->Result<Ready> {
        poll(gate)?;
        if self.longitudinal_waves.is_none()
            && !matches!(&self.longitudinal_loss, LongitudinalMaterialLoss::Elastic) {
            return Err("longitudinal material loss requires explicit axial coordinates".into());
        }
        if self.gas != *self.gas_record.state() {
            return Err("material-string gas coefficients must retain their original source record".into());
        }
        self.validate_shell_observer()?;
        if let Some(layers)=self.friction.shear_layers() { layers.preflight(&self)?; }
        if let Some(body)=&profile { body.preflight(&self)?; }
        if self.one_way_observer.is_some() && profile.is_none() {
            return Err("one_way_observer requires a body profile that declares its closed acoustic skin".into());
        }
        if let Some(stops)=&stop_profile { stops.validate_performance(&self,gate)?; }
        if let Some(radiation) = self.body_radiation {
            radiation.validate(self.rate)?;
            if self.prepared_linear.is_some() && !matches!(radiation, input::BodyRadiation::Skin(_)) {
                return Err("computed body radiation requires the exact owner with acoustic storage".into());
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
            .and_then(|n| n.checked_add(profile.as_ref().map_or(0,body::Profile::material_memory_states)))
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
            if memory > 256 || terms > input::MAX_WORK as usize {
                return Err("complete string/body memory exceeds the original state or setup-work ceiling".into());
            }
            settings.network.max_setup_terms = terms;
        }
        let bridge_interface_bytes = self.bridge_interfaces.pack_bytes()?;
        let interface_bytes = profile.as_ref().map(body::Profile::interface_pack_bytes)
            .transpose()?.unwrap_or(0).checked_add(self.friction.pack_bytes()?)
            .and_then(|bytes| bytes.checked_add(bridge_interface_bytes))
            .ok_or("combined interface pack byte count overflow")?;
        self.bridge_interfaces.apply(self.temperature, self.omega, &mut self.strings)?;
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
        self.friction.validate_surfaces(&self, &catalog)?;
        let bank=if let Some(body_profile)=&profile {
            body_profile.compile(&self,sources,&mut catalog,&mut loader,settings.clone(),cx,gate)?
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
        let bank = if profile.as_ref().is_some_and(body::Profile::has_thermal_cavity) {
            bank // Gas inertia precedes the material/thermal memory prefix.
        } else if let Some(resolution) = self.cylinder_gas_loading {
            bank.with_cylinder_gas_loading(resolution, gate)?
        } else { bank };
        let bank = if let Some((slope, config)) = self.geometric_stretching {
            bank.with_geometric_stretching(slope, config, gate)?
        } else { bank };
        let bank = if let Some((resolution, config)) = &self.longitudinal_waves {
            bank.with_longitudinal_waves_and_loss(vec![resolution.clone(); self.strings.len()],
                self.longitudinal_loss.clone(), *config, gate)?
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
        let bank=match &stop_profile { Some(stops)=>stops.bind(bank,gate)?, None=>bank };
        let bank=if let Some(layers)=self.friction.shear_layers() {
            layers.bind(&self,bank,&mut catalog,&mut loader,gate)?
        } else { bank };
        let performance_config=MultiBowPerformanceConfig { samples:self.samples,max_block:self.block,
            max_work:input::MAX_WORK,max_events:input::MAX_CONTROLS };
        let mut performance=if self.bow_ribbon.is_some() {
            bank.into_scheduled_ribbon(self.schedule.clone(),self.schedule.tracks()[0].id.clone(),performance_config)?
        } else {
            let names=self.schedule.tracks().iter().map(|t|t.id.clone()).collect();
            bank.into_scheduled(self.schedule.clone(),names,performance_config)?
        };
        if let Some(forces)=&self.forces { performance=performance.with_string_forces(forces.clone(),gate)?; }
        if let Some(stops)=&stop_profile { performance=stops.bind_motion(performance,gate)?; }
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
/// Take `--observer-store DIR` out of the option pairs; the rest are profile options.
fn observer_store_option(args:&[std::ffi::OsString])->Result<(Option<&Path>,Vec<std::ffi::OsString>)> {
    let (mut store,mut rest)=(None,Vec::new());
    for pair in args.chunks(2) {
        if pair[0]!="--observer-store" { rest.extend_from_slice(pair); continue; }
        if store.is_some() { return Err("--observer-store cannot be repeated".into()); }
        store=Some(Path::new(pair.get(1).ok_or("--observer-store needs a directory")?));
    }
    Ok((store,rest))
}
fn main()->Result<()> {
    let raw=std::env::args_os().collect::<Vec<_>>();
    if raw.get(1).is_some_and(|a| a=="assemble-shell") { return assemble_shell(&raw[2..]); }
    let (args, gas_options)=thermochemical::Options::extract_os(&raw)?;
    if args.len()<3 {
        return Err("usage: material_strings PERFORMANCE NEW_OUTPUT.wav [--body-profile PROFILE] [--stop-profile PROFILE] [--observer-store DIR] [--gas-component SPECIES NASA9 TRANSPORT --gas-water none|ID --gas-phase frozen|liquid-water]".into());
    }
    if gas_options.record.is_some() {
        return Err("material-string gas records belong in the performance file".into());
    }
    let (store_dir,options)=observer_store_option(&args[3..])?;
    let (body_path,stop_path)=stops::profile_paths(&options)?;
    let mut profile=body_path.map(|p|body::Profile::parse(&read_bounded(p,body::MAX_SHELL_BYTES)?)).transpose()?;
    if let Some(dir)=store_dir {
        profile.as_mut().ok_or("--observer-store requires --body-profile")?.keep_observer_setups(dir)?;
    }
    let stop_profile=stop_path.map(|p|stops::Profile::parse(&read_bounded(p,stops::MAX_BYTES)?)).transpose()?;
    let path=Path::new(&args[1]);let parent=path.parent().unwrap_or(Path::new("."));
    let gas_sources=gas_options.load()?;
    let gas_bindings=gas_sources.bindings();
    let gas_registry=gas_sources.registry(&gas_bindings);
    let mut spec=input::Spec::parse_with_gas_sources(
        &read_bounded(path,input::MAX_INPUT_BYTES)?,gas_registry.as_ref())?;
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
    spec.friction.load(&gate,|name| {
        let bytes=read_bounded(&parent.join(name),input::MAX_PACK_BYTES.min(remaining))?;
        remaining=remaining.checked_sub(bytes.len()).ok_or("combined source pack byte budget exceeded")?;
        Ok(NormalizedInterfacePack::from_bytes(&bytes)?)
    })?;
    spec.bridge_interfaces.load(spec.temperature, spec.omega, &gate, |name| {
        let bytes=read_bounded(&parent.join(name),input::MAX_PACK_BYTES.min(remaining))?;
        remaining=remaining.checked_sub(bytes.len()).ok_or("combined source pack byte budget exceeded")?;
        Ok(NormalizedInterfacePack::from_bytes(&bytes)?)
    })?;
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
    eprintln!("bow_friction={:?}",ready.spec.friction);
    for port in 0..ready.spec.bow_contact_count() {
        if let Some(layer)=ready.performance.bank().bow_shear_layer(port) {
            eprintln!("bow_shear_layer_port={port}; source={}; area_m2={}; thickness_m={}; frozen_T_K={}; stiffness_N_m={}; yield_N={}; plastic_resistance_N_s_m={}; max_elastic_shear={}",
                layer.resolved().identity(),layer.area_m2(),layer.thickness_m(),layer.temperature_kelvin(),
                layer.stiffness_n_m(),layer.yield_force_n(),layer.plastic_resistance_n_s_m(),layer.maximum_elastic_shear_strain());
        }
    }
    if !ready.spec.bridge_interfaces.bindings().is_empty() {
        eprintln!("bridge_interfaces={:?}",ready.spec.bridge_interfaces);
    }
    if let Some(ribbon)=ready.performance.bank().source().ribbon() { eprintln!("bow_ribbon={ribbon:?}"); }
    if let Some(fibers)=ready.performance.bank().source().fiber_ribbon() {
        eprintln!("bow_fibers={:?}; resolution={:?}; mass_kg={}; station={}; frozen-station rigid-end model, not thermal rosin",
            fibers.model().geometry(),fibers.model().resolution(),fibers.model().mass_kg(),fibers.station_fraction());
    }
    if let Some(profile)=&ready.body_profile { eprintln!("body_profile={profile:?}"); }
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
    if let Some(budget)=ready.performance.bank().source().bank().prepared_linear_budget() {
        eprintln!("prepared_linear={budget:?}");
    }
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
    // Complete source, geometry, frequency, schedule and preload admission is
    // finished before creating output. Never overwrite any existing path.
    let output=OpenOptions::new().write(true).create_new(true).open(&args[2])?;
    let (_,summary)=ready.write(output,&gate).map_err(|e|format!("WAV not completed; keep/discard its partial prefix explicitly: {e}"))?;
    eprintln!("samples={}; rate_Hz={}; clips={}; full_scale_Pa={}",summary.samples,summary.sample_rate_hz,
        summary.clipped_samples,summary.full_scale_pa);
    // Steps whose stop missed its tolerance at the actual endpoint and passed it after a headroom re-solve.
    eprintln!("obstacle_resolves={}",ready.performance.bank().system().obstacle_resolves());
    Ok(())
}

#[cfg(test)]
#[path = "material_strings/tests.rs"]
mod tests;

