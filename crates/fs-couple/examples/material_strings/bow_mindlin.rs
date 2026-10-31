//! Explicit circular Hertz-patch data reaches the original source Mindlin law.
//! No fitted stiffness, coefficient, contact radius or loading path is inferred.
use super::*;
use fs_contact::interface_binding::{NormalContactModelSelection, mindlin::SourceMindlinContact};
use fs_contact::normal_patch::{source::{PointContactEnvelope, PointContactGeometry,
    state::PointContactMaterialCards}};
use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::MultiBowStringBank;
use fs_matdb::{InterfaceSystemCard, NormalizedModelPack};
use fs_material::state_point::{append_held_environment, query_with_source_conventions};
use fs_qty::{QuantitySpec, semantic::{FrequencyConvention, QuantityKind, SemanticType, ValueForm}};
use fs_tribo::InputAuthority;

/// One explicitly resolved circular Hertz asperity at one fixed string station.
/// This is not a cylinder/ribbon geometry approximation or a per-cell multiplier.
#[derive(Clone, Debug)]
pub struct Declaration {
    pub track: String,
    pub receiver_pack: String,
    pub interface_pack: String,
    pub normal_model_pack: String,
    pub friction_property: String,
    pub radius_m: f64,
    pub indentation_m: f64,
    pub envelope: PointContactEnvelope,
    pub maximum_history_intervals: usize,
}
#[derive(Clone)]
pub struct Inputs {
    declarations: Vec<Declaration>,
    interfaces: BTreeMap<String, NormalizedInterfacePack>,
    models: BTreeMap<String, NormalizedModelPack>,
}
impl core::fmt::Debug for Inputs {
    fn fmt(&self, formatter: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        let interfaces = self.interfaces.iter().map(|(path, pack)| (path, pack.content_hash())).collect::<Vec<_>>();
        let models = self.models.iter().map(|(path, pack)| (path, pack.content_hash())).collect::<Vec<_>>();
        formatter.debug_struct("SourceMindlinInputs").field("declarations", &self.declarations)
            .field("interfaces", &interfaces).field("models", &models).finish()
    }
}
impl Inputs {
    pub fn new(declarations: Vec<Declaration>) -> Result<Self> {
        if declarations.is_empty() || declarations.len() > 16
            || declarations.iter().enumerate().any(|(i,d)| d.track.is_empty() || d.track.len() > 128
                || declarations[..i].iter().any(|previous| previous.track == d.track)) {
            return Err("Mindlin sources require 1..=16 distinct original track IDs".into());
        }
        for d in &declarations {
            if !(1..=fs_tribo::mindlin::MAX_HISTORY_INTERVALS).contains(&d.maximum_history_intervals)
                || d.indentation_m > d.envelope.maximum_indentation_m
                || d.envelope.maximum_absolute_rate_m_per_s < 0.0
                || [d.envelope.uncertainty.radius_relative, d.envelope.uncertainty.modulus_relative,
                    d.envelope.uncertainty.load_relative].iter().any(|x| *x < 0.0) {
                return Err("bow-mindlin requires bounded positive indentation/history and nonnegative rate/uncertainty".into());
            }
            Requirement::try_new(&d.friction_property, Dims::NONE, Admit::NonNegative)?;
        }
        Ok(Self { declarations, interfaces: BTreeMap::new(), models: BTreeMap::new() })
    }
    pub fn track_count(&self) -> usize { self.declarations.len() }
    pub fn contains_track(&self, track: &str) -> bool { self.declarations.iter().any(|d| d.track == track) }
    pub fn has_source(&self, name: &str) -> bool { self.interfaces.contains_key(name) }
    pub fn preflight(&self, spec: &input::Spec) -> Result<()> {
        spec.friction.validate_tracks(&spec.schedule)?;
        if spec.bow_ribbon.is_some() || spec.bow_fibers.is_some()
            || spec.bow_station_motion.is_some() || spec.bow_residual_flexibility {
            return Err("fixed circular Mindlin patches require fixed point-bow geometry; ribbon/cylinder subdivision, migration and residual-tail history are not inferred".into());
        }
        if spec.contact_thermal || spec.dry_contact_thermal {
            return Err("fixed-source Mindlin playback is isothermal; finite-cell thermal feedback needs explicit constitutive history/energy transfer".into());
        }
        for track in spec.schedule.tracks().iter().filter(|track| self.contains_track(&track.id)) {
            for value in std::iter::once(&track.initial).chain(track.events.iter().map(|event| &event.value)) {
                if let GestureValue::Bow { normal_force_n, .. } = value {
                    if *normal_force_n != 0.0 {
                        return Err("bow-mindlin requires zero additional prescribed normal load in every gesture; retained geometry/materials supply the Hertz preload".into());
                    }
                }
            }
        }
        Ok(())
    }
    pub fn load(&mut self, gate: &CancelGate, mut loader: impl FnMut(&str) -> Result<Vec<u8>>) -> Result<()> {
        let mut interfaces = BTreeMap::new();
        let mut models = BTreeMap::new();
        let mut remaining = input::MAX_TOTAL_PACK_BYTES;
        let mut read = |path: &str| -> Result<Vec<u8>> {
            poll(gate)?;
            let bytes = loader(path)?;
            if bytes.len() > input::MAX_PACK_BYTES || bytes.len() > remaining {
                return Err("Mindlin source artifacts exceed the shared material byte budget".into());
            }
            remaining -= bytes.len();
            Ok(bytes)
        };
        for d in &self.declarations {
            if !interfaces.contains_key(&d.interface_pack) {
                interfaces.insert(d.interface_pack.clone(), NormalizedInterfacePack::from_bytes(&read(&d.interface_pack)?)?);
            }
            if !models.contains_key(&d.normal_model_pack) {
                models.insert(d.normal_model_pack.clone(), NormalizedModelPack::from_bytes(&read(&d.normal_model_pack)?)?);
            }
        }
        poll(gate)?;
        self.interfaces = interfaces; self.models = models;
        Ok(())
    }
    pub fn pack_bytes(&self) -> Result<usize> {
        if self.declarations.iter().any(|d| !self.interfaces.contains_key(&d.interface_pack)
            || !self.models.contains_key(&d.normal_model_pack)) {
            return Err("load every Mindlin interface and normal-model source before construction".into());
        }
        self.interfaces.values().map(|p| p.to_bytes().len()).chain(self.models.values().map(|p| p.to_bytes().len()))
            .try_fold(0usize, |sum, size| sum.checked_add(size).ok_or_else(|| "Mindlin source byte count overflow".into()))
    }
    pub fn bind(&self, spec: &input::Spec, bank: MultiBowStringBank, catalog: &mut Catalog,
        loader: &mut impl FnMut(&str) -> Result<MaterialPack>, gate: &CancelGate) -> Result<MultiBowStringBank>
    {
        self.preflight(spec)?;
        let mut bindings = Vec::new();
        for (port, track) in spec.schedule.tracks().iter().enumerate() {
            let Some(d) = self.declarations.iter().find(|d| d.track == track.id) else { continue; };
            poll(gate)?;
            let member = fs_scenario::gesture::path::legacy_member(&track.target, "bow")
                .ok_or("Mindlin source requires an original material-string bow track")?;
            let wire = spec.strings.get(member).ok_or("Mindlin track names an absent string")?;
            let (outer, depth) = match &wire.section {
                input::WireSection::Measured { .. } => return Err("whole-string measurements do not identify Mindlin contact material".into()),
                input::WireSection::Homogeneous => (&wire.pack, wire.radius),
                input::WireSection::CloseWound { pack, wire_radius_m, .. } => (pack, *wire_radius_m),
                input::WireSection::Bonded { core_radius_m, sleeves } => {
                    let layer = sleeves.last().ok_or("Mindlin contact needs an actual outer constituent")?;
                    let inner = sleeves.iter().rev().nth(1).map_or(*core_radius_m, |s| s.outer_radius_m);
                    (&layer.pack, layer.outer_radius_m - inner)
                }
            };
            if d.envelope.layer_thickness_m > depth || d.envelope.half_space_depth_m > depth {
                return Err("Mindlin half-space/layer envelope exceeds the string's actual outer material depth".into());
            }
            load_pack(&d.receiver_pack, catalog, loader)?;
            let a = &catalog.get(outer).ok_or("Mindlin outer material is absent from the original catalog")?.card;
            let b = &catalog[&d.receiver_pack].card;
            let environment = |name: &str| spec.material_states.get(name).map_or(&[][..], Vec::as_slice);
            let bulk = [material_inputs::query_with_environment(a, spec.temperature, spec.omega, environment(outer))?,
                material_inputs::query_with_environment(b, spec.temperature, spec.omega, environment(&d.receiver_pack))?];
            let interface = self.interfaces.get(&d.interface_pack).ok_or("Mindlin interface pack is not loaded")?;
            let models = self.models.get(&d.normal_model_pack).ok_or("Mindlin model pack is not loaded")?;
            let absolute = QuantitySpec::semantic(SemanticType::new(QuantityKind::AbsoluteTemperature, ValueForm::Static));
            let angular = QuantitySpec::semantic(SemanticType::new(QuantityKind::Frequency(FrequencyConvention::Angular), ValueForm::Static));
            let point = query_with_source_conventions(interface.claims_pack().claims(),
                &[("T", absolute, spec.temperature), ("omega", angular, spec.omega)])?;
            let point = append_held_environment(interface.claims_pack().claims(), point, environment(&d.interface_pack))?;
            let card = InterfaceSystemCard::assemble(interface.card().surface_a().clone(), interface.card().surface_b().clone(),
                interface.card().context().clone(), interface.claims_pack().claims().clone(), models.models().to_vec())?;
            let materials = PointContactMaterialCards { surface_a: a, surface_b: b, interface: &card,
                bulk_selection: [Pick::SingleClaimOnly, Pick::SingleClaimOnly], interface_selection: Pick::SingleClaimOnly,
                model_selection: NormalContactModelSelection::SingleSupported, authority: InputAuthority::CallerDeclared };
            let normal = materials.resolve_at_fixed_environments([&bulk[0], &bulk[1]], &point,
                PointContactGeometry::SpherePlane { effective_radius_m: d.radius_m },
                PointContactEnvelope { maximum_step_s: 1.0 / f64::from(spec.rate), ..d.envelope })?;
            let source = SourceMindlinContact::resolve(normal, d.indentation_m,
                Requirement::try_new(&d.friction_property, Dims::NONE, Admit::NonNegative)?,
                Pick::SingleClaimOnly, d.maximum_history_intervals)?;
            bindings.push((port, source));
        }
        poll(gate)?;
        Ok(bank.with_bow_mindlin_contacts(bindings, gate)?)
    }
}
