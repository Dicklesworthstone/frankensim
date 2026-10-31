//! Explicit source inputs for the original Hertz/Hunt--Crossley point owner.
//! No contact constants, specimen properties or relative curvatures are inferred.
use super::*;
use fs_contact::interface_binding::NormalContactModelSelection;
use fs_contact::normal_patch::InputUncertainty;
use fs_contact::normal_patch::source::{PointContactEnvelope, PointContactGeometry, SourcePointContact};
use fs_contact::normal_patch::source::state::PointContactMaterialCards;
use fs_matdb::{InterfaceSystemCard, NormalizedModelPack};
use fs_material::state_point::{append_held_environment, query_with_source_conventions};
use fs_qty::{QuantitySpec, semantic::{FrequencyConvention, QuantityKind, SemanticType, ValueForm}};
use fs_tribo::InputAuthority;
use fs_couple::bowed_string::runtime::reciprocal::bank::shell::MaterialShellPoint;

#[derive(Clone, Debug)]
enum ReceiverMaterial {
    Pack(String),
    Body,
}

#[derive(Clone, Debug)]
pub(super) struct Declaration {
    pub(super) gap_m: f64,
    receiver: ReceiverMaterial,
    interface: String,
    model: String,
    geometry: PointContactGeometry,
    envelope: PointContactEnvelope,
}
impl Declaration {
    pub(super) fn parse(values: &[&str], gap_m: f64, maximum_indentation_m: f64,
        receiver: Receiver) -> Result<Self> {
        let path = |word: &str| -> Result<String> {
            if word.is_empty() || word.len() > 512 || !Path::new(word).components()
                .all(|p| matches!(p, std::path::Component::Normal(_))) {
                return Err("normal source pack paths must stay relative to the performance directory".into());
            }
            Ok(word.into())
        };
        let maximum_curvature_m_inverse = scalar(values[5])?;
        let minimum_curvature_m_inverse = scalar(values[6])?;
        let envelope = PointContactEnvelope {
            maximum_indentation_m,
            maximum_absolute_rate_m_per_s: scalar(values[7])?,
            // Replaced by the actual mechanical clock before source admission.
            maximum_step_s: 0.0,
            half_space_depth_m: scalar(values[8])?,
            layer_thickness_m: scalar(values[9])?,
            uncertainty: InputUncertainty { radius_relative: scalar(values[10])?,
                modulus_relative: scalar(values[11])?, load_relative: scalar(values[12])? },
        };
        if minimum_curvature_m_inverse <= 0.0
            || maximum_curvature_m_inverse < minimum_curvature_m_inverse
            || envelope.maximum_absolute_rate_m_per_s < 0.0
            || envelope.half_space_depth_m <= 0.0 || envelope.layer_thickness_m <= 0.0
            || [envelope.uncertainty.radius_relative, envelope.uncertainty.modulus_relative,
                envelope.uncertainty.load_relative].iter().any(|u| *u < 0.0) {
            return Err("material-stop requires positive relative curvatures/depths and nonnegative rate/uncertainties".into());
        }
        let receiver = match receiver {
            Receiver::Shell(_) if values[2] == "@receiver" => ReceiverMaterial::Body,
            Receiver::Fixed if values[2] != "@receiver" => ReceiverMaterial::Pack(path(values[2])?),
            _ => return Err("shell material-stop requires @receiver; fixed contact requires an explicit bulk pack".into()),
        };
        Ok(Self { gap_m, receiver, interface: path(values[3])?,
            model: path(values[4])?,
            geometry: PointContactGeometry::EllipticParaboloid {
                maximum_curvature_m_inverse, minimum_curvature_m_inverse }, envelope })
    }
}

#[derive(Clone, Debug, Default)]
pub(super) struct Sources {
    interfaces: BTreeMap<String, NormalizedInterfacePack>,
    models: BTreeMap<String, NormalizedModelPack>,
}
impl Profile {
    /// Whether a held state names an original decoded material-stop interface.
    pub fn has_material_interface(&self, path: &str) -> bool {
        self.sources.interfaces.contains_key(path)
    }

    /// All-or-nothing cold decode through the same bounded attachment/file owner.
    /// Each distinct interface/model path is read once and retained verbatim as a
    /// canonical source pack. Bulk packs continue through the common catalog.
    pub fn load_sources(&mut self, gate: &CancelGate,
        mut loader: impl FnMut(&str) -> Result<Vec<u8>>) -> Result<()> {
        let mut sources = Sources::default();
        let mut remaining = input::MAX_TOTAL_PACK_BYTES;
        let mut read = |path: &str| -> Result<Vec<u8>> {
            poll(gate)?;
            let bytes = loader(path)?;
            if bytes.len() > input::MAX_PACK_BYTES || bytes.len() > remaining {
                return Err("normal source packs exceed the shared material byte budget".into());
            }
            remaining -= bytes.len();
            Ok(bytes)
        };
        for entry in &self.entries {
            if let Law::Material(d) = &entry.law {
                if !sources.interfaces.contains_key(&d.interface) {
                    let pack = NormalizedInterfacePack::from_bytes(&read(&d.interface)?)?;
                    sources.interfaces.insert(d.interface.clone(), pack);
                }
                if !sources.models.contains_key(&d.model) {
                    let pack = NormalizedModelPack::from_bytes(&read(&d.model)?)?;
                    sources.models.insert(d.model.clone(), pack);
                }
            }
        }
        poll(gate)?;
        self.sources = sources;
        Ok(())
    }

    pub fn source_pack_bytes(&self) -> Result<usize> {
        for entry in &self.entries {
            if let Law::Material(d) = &entry.law {
                if !self.sources.interfaces.contains_key(&d.interface) || !self.sources.models.contains_key(&d.model) {
                    return Err("load every material-stop interface and model pack before construction".into());
                }
            }
        }
        self.sources.interfaces.values().map(|p| p.to_bytes().len())
            .chain(self.sources.models.values().map(|p| p.to_bytes().len()))
            .try_fold(0usize, |n, size| n.checked_add(size).ok_or_else(|| "normal source pack byte count overflow".into()))
    }
}

impl Sources {
    pub(super) fn resolve(&self, entries: &[Entry], shell_points: &[MaterialShellPoint<'_>],
        spec: &input::Spec, catalog: &mut Catalog,
        loader: &mut impl FnMut(&str) -> Result<MaterialPack>, gate: &CancelGate)
        -> Result<Vec<(usize, SourcePointContact)>> {
        let mut sources = Vec::new();
        for (index, entry) in entries.iter().enumerate() {
            let Law::Material(d) = &entry.law else { continue; };
            poll(gate)?;
            let wire = spec.strings.get(entry.string).ok_or("material stop names an absent string")?;
            // The contact sees the actual outer constituent, not a core material
            // or a homogenized section modulus. The declared contact envelope
            // cannot exceed that constituent's real radial depth.
            let (outer, depth) = match &wire.section {
                input::WireSection::Measured { .. } => return Err("whole-string measurements do not identify a material contact surface".into()),
                input::WireSection::Homogeneous => (&wire.pack, wire.radius),
                input::WireSection::CloseWound { pack, wire_radius_m, .. } => (pack, *wire_radius_m),
                input::WireSection::Bonded { core_radius_m, sleeves } => {
                    let layer = sleeves.last().ok_or("material stop needs an actual outer layer")?;
                    let inner = sleeves.iter().rev().nth(1).map_or(*core_radius_m, |s| s.outer_radius_m);
                    (&layer.pack, layer.outer_radius_m-inner)
                },
            };
            if d.envelope.layer_thickness_m > depth || d.envelope.half_space_depth_m > depth {
                return Err("material-stop half-space/layer depth exceeds the string's actual outer constituent".into());
            }
            let receiver = match &d.receiver {
                ReceiverMaterial::Pack(path) => { load_pack(path, catalog, loader)?; None },
                ReceiverMaterial::Body => {
                    let Receiver::Shell(part) = entry.receiver else {
                        return Err("source receiver lost its original shell geometry".into());
                    };
                    let point = shell_points.get(part).ok_or("source receiver names an absent shell point")?;
                    if d.envelope.layer_thickness_m > point.thickness_m
                        || d.envelope.half_space_depth_m > point.thickness_m {
                        return Err("material-stop half-space/layer depth exceeds the receiver's actual section thickness".into());
                    }
                    // The current point owner is isotropic. Merely finding a
                    // similarly named isotropic property elsewhere on an
                    // orthotropic card cannot replace its actual body model.
                    if !point.isotropic_elasticity || ["density", "young_modulus", "poisson_ratio"].iter()
                        .any(|name| point.material.property(name).is_none()) {
                        return Err("source point contact requires the receiver's actual isotropic elastic state".into());
                    }
                    Some(point.material)
                },
            };
            let a = &catalog.get(outer).ok_or("outer string material is absent from the catalog")?.card;
            let b = match &d.receiver {
                ReceiverMaterial::Pack(path) => &catalog[path].card,
                ReceiverMaterial::Body => {
                    let actual = receiver.ok_or("source receiver lost its original material")?;
                    &catalog.values().find(|pack| pack.card.content_hash() == actual.card_identity())
                        .ok_or("receiver's original immutable material card is absent from the catalog")?.card
                },
            };
            let environment = |path: &str| spec.material_states.get(path)
                .map_or(&[][..], Vec::as_slice);
            let point = material_inputs::query_with_environment(a, spec.temperature, spec.omega,
                environment(outer))?;
            let receiver_point = match &d.receiver {
                ReceiverMaterial::Pack(path) => material_inputs::query_with_environment(
                    b, spec.temperature, spec.omega, environment(path))?,
                ReceiverMaterial::Body => receiver.ok_or("source receiver lost its original material")?.query(),
            };
            let interface = self.interfaces.get(&d.interface).ok_or("material stop interface pack is not loaded")?;
            let models = self.models.get(&d.model).ok_or("material stop model pack is not loaded")?;
            let absolute = QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature, ValueForm::Static));
            let angular = QuantitySpec::semantic(SemanticType::new(
                QuantityKind::Frequency(FrequencyConvention::Angular), ValueForm::Static));
            let interface_point = query_with_source_conventions(interface.claims_pack().claims(),
                &[("T", absolute, spec.temperature), ("omega", angular, spec.omega)])?;
            let interface_point = append_held_environment(interface.claims_pack().claims(),
                interface_point, environment(&d.interface))?;
            // This explicit record associates two existing normalized artifacts;
            // the original interface/model cards still enforce ordered surfaces,
            // property validity, exact supported law schema and source authority.
            let card = InterfaceSystemCard::assemble(
                interface.card().surface_a().clone(), interface.card().surface_b().clone(),
                interface.card().context().clone(), interface.claims_pack().claims().clone(),
                models.models().to_vec())?;
            let materials = PointContactMaterialCards {
                surface_a: a, surface_b: b, interface: &card,
                bulk_selection: [Pick::SingleClaimOnly, Pick::SingleClaimOnly],
                interface_selection: Pick::SingleClaimOnly,
                model_selection: NormalContactModelSelection::SingleSupported,
                authority: InputAuthority::CallerDeclared,
            };
            let envelope = PointContactEnvelope {
                maximum_step_s: 1.0/f64::from(spec.rate), ..d.envelope };
            let independent_environment = [&point, &receiver_point, &interface_point].iter()
                .any(|p| p.axes().keys().any(|axis| axis != "T" && axis != "omega"));
            let source = if independent_environment {
                materials.resolve_at_fixed_environments([&point, &receiver_point],
                    &interface_point, d.geometry, envelope)
            } else {
                materials.resolve_at_temperatures([&point, &receiver_point],
                    &interface_point, d.geometry, envelope)
            }?;
            if let Some(actual) = receiver {
                let contact = source.model().normal_state().elastic().surface_b_state().resolved();
                if contact.card_identity() != actual.card_identity()
                    || contact.query_point() != actual.query_point()
                    || contact.axis_quantities() != actual.axis_quantities()
                    || ["density", "young_modulus", "poisson_ratio"].iter().any(|name| {
                        match (actual.property(name), contact.property(name)) {
                            (Some(body), Some(patch)) => body.answer() != patch.answer(),
                            _ => true,
                        }
                    }) {
                    return Err("normal contact must retain the receiver body's exact elastic source receipts and state".into());
                }
            }
            sources.push((index, source));
        }
        poll(gate)?;
        Ok(sources)
    }
}
