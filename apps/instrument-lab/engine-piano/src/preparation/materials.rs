//! Canonical material packs reach the original live string/bridge/hammer owner.
//! Geometry is explicit; neither a material name nor a density implies radius.
//! These are perfectly bonded concentric sections at a fixed resolved state,
//! with sourced Kelvin–Voigt bending loss. Helical winding, interlayer slip and
//! evolving temperature require their own physical owners and are not inferred.

use std::{collections::{BTreeMap, BTreeSet}, sync::Arc};
use fs_matdb::{NormalizedMaterialCardPack, QueryPoint};
use fs_material::{section::{ConcentricElasticSection, ConcentricKelvinVoigtSection,
    ConcentricMaterialLayer, MAX_CONCENTRIC_LAYERS}, state_point::{DENSITY_PROPERTY,
    YOUNG_MODULUS_PROPERTY, MaterialPropertySelection, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point}};
use fs_qty::{Density, DynViscosity, Pressure};
use crate::refmods::geometry::Course;

pub const MAX_PACK_BYTES: usize = 16 * 1024 * 1024;
pub const MAX_PACKS: usize = 128;
pub const HEADER: &str = "frankensim-piano-materials-v1";
pub const COLUMNS: &str = "key,outer_radius_m,material_index";

/// Explicit browser operating point in absolute kelvin and cyclic hertz.
/// Axis names select original source coordinates; units and conventions are
/// fixed by this entry and cannot be inferred from an axis spelling.
pub fn operating_point(temperature_axis: &str, temperature_k: f64,
    frequency_axis: &str, frequency_hz: f64) -> Result<QueryPoint, String> {
    use fs_qty::{QuantitySpec, semantic::{FrequencyConvention, QuantityKind,
        SemanticType, ValueForm}};
    if [temperature_axis, frequency_axis].iter().any(|axis|
        axis.trim().is_empty() || axis.trim() != *axis || axis.len() > 128)
        || temperature_axis == frequency_axis
        || !temperature_k.is_finite() || temperature_k <= 0.0
        || !frequency_hz.is_finite() || frequency_hz <= 0.0 {
        return Err("declare distinct temperature/frequency axes and positive finite K/Hz values".into());
    }
    let typed = |kind| QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static));
    QueryPoint::new().with_quantity(temperature_axis, typed(QuantityKind::AbsoluteTemperature), temperature_k)
        .and_then(|point| point.with_quantity(frequency_axis,
            typed(QuantityKind::Frequency(FrequencyConvention::Cyclic)), frequency_hz))
        .map_err(|e| e.to_string())
}

/// Admitted source owners and their actual section integrals. Cloning a live
/// engine retains these immutable packs, including original source receipts.
#[derive(Clone, Debug)]
pub struct StringMaterials {
    sections: BTreeMap<u8, ConcentricKelvinVoigtSection>,
    packs: Vec<Arc<NormalizedMaterialCardPack>>,
}
impl StringMaterials {
    /// Read original normalized packs, then assign them by explicit zero-based
    /// index to the spatially ordered annuli in each key's geometry rows. All
    /// three coefficients resolve from one complete query per original card.
    /// Missing/ambiguous claims, wrong dimensions or unsupported state refuse.
    pub fn from_packs(bytes: &[&[u8]], geometry: &str, point: &QueryPoint,
        viscosity_property: &str) -> Result<Self, String> {
        if bytes.is_empty() || bytes.len() > MAX_PACKS {
            return Err("string materials require 1..128 canonical packs".into());
        }
        let mut size = 0usize;
        for pack in bytes {
            size = size.checked_add(pack.len()).ok_or("material pack byte overflow")?;
            if pack.is_empty() || size > MAX_PACK_BYTES {
                return Err("nonempty material packs must fit the 16 MiB input budget".into());
            }
        }
        if geometry.len() > super::MAX_SCALE_BYTES {
            return Err("string material geometry exceeds the 256 KiB input budget".into());
        }
        if viscosity_property.trim().is_empty() || viscosity_property.len() > 128
            || viscosity_property.trim() != viscosity_property
            || [DENSITY_PROPERTY, YOUNG_MODULUS_PROPERTY].contains(&viscosity_property) {
            return Err("declare a distinct nonempty solid bending-viscosity property".into());
        }
        let requirements = [
            ScalarPropertyRequirement::try_new(DENSITY_PROPERTY, Density::DIMS,
                ScalarAdmissibility::StrictlyPositive),
            ScalarPropertyRequirement::try_new(YOUNG_MODULUS_PROPERTY, Pressure::DIMS,
                ScalarAdmissibility::StrictlyPositive),
            ScalarPropertyRequirement::try_new(viscosity_property, DynViscosity::DIMS,
                ScalarAdmissibility::Finite),
        ].into_iter().collect::<Result<Vec<_>, _>>().map_err(|e| e.to_string())?;
        let mut packs = Vec::with_capacity(bytes.len());
        let mut states = Vec::with_capacity(bytes.len());
        for bytes in bytes {
            let pack = Arc::new(NormalizedMaterialCardPack::from_bytes(bytes)
                .map_err(|e| e.to_string())?);
            let state = resolve_material_state_point(pack.card(), point, &requirements,
                MaterialPropertySelection::SingleClaimOnly).map_err(|e| e.to_string())?;
            packs.push(pack);
            states.push(state);
        }
        let mut lines = geometry.lines().map(|line| line.split('#').next().unwrap_or("").trim())
            .filter(|line| !line.is_empty());
        if lines.next() != Some(HEADER) || lines.next() != Some(COLUMNS) {
            return Err(format!("string material geometry requires {HEADER} and {COLUMNS}"));
        }
        let mut layers = BTreeMap::<u8, Vec<ConcentricMaterialLayer>>::new();
        for line in lines {
            let fields = line.split(',').map(str::trim).collect::<Vec<_>>();
            if fields.len() != 3 { return Err("string material row needs key,radius,index".into()); }
            let midi = fields[0].parse::<u8>().map_err(|_| "invalid string material key")?;
            if !(21..=108).contains(&midi) { return Err("string material key must be 21..108".into()); }
            let radius = fields[1].parse::<f64>().map_err(|_| "invalid string outer radius")?;
            let index = fields[2].parse::<usize>().map_err(|_| "invalid material pack index")?;
            let material = states.get(index).ok_or("string material pack index is out of range")?;
            let course = layers.entry(midi).or_default();
            if course.len() >= MAX_CONCENTRIC_LAYERS {
                return Err("string section exceeds the 128-layer input budget".into());
            }
            course.push(ConcentricMaterialLayer { outer_radius_m: radius, material: material.clone() });
        }
        if layers.is_empty() { return Err("string material geometry has no sections".into()); }
        let sections = layers.into_iter().map(|(key, layers)| {
            let section = ConcentricElasticSection::try_new(0.0, layers)
                .and_then(|section| section.with_kelvin_voigt_bending(viscosity_property))
                .map_err(|e| format!("string material key {key}: {e}"))?;
            Ok((key, section))
        }).collect::<Result<BTreeMap<_, _>, String>>()?;
        Ok(Self { sections, packs })
    }

    /// Bounded WASM framing of the unchanged canonical bytes; this is not an
    /// alternate material serialization. Every byte must belong to one pack.
    pub fn from_packed_bytes(bytes: &[u8], lengths: &[u32], geometry: &str,
        point: &QueryPoint, viscosity_property: &str) -> Result<Self, String> {
        if bytes.len() > MAX_PACK_BYTES || lengths.is_empty() || lengths.len() > MAX_PACKS {
            return Err("string material packet exceeds the pack count or byte budget".into());
        }
        let mut offset = 0usize;
        let mut packs = Vec::with_capacity(lengths.len());
        for &length in lengths {
            let end = offset.checked_add(length as usize).ok_or("material framing overflow")?;
            packs.push(bytes.get(offset..end).ok_or("truncated material pack packet")?);
            offset = end;
        }
        if offset != bytes.len() { return Err("trailing bytes after string material packs".into()); }
        Self::from_packs(&packs, geometry, point, viscosity_property)
    }

    pub fn sections(&self) -> &BTreeMap<u8, ConcentricKelvinVoigtSection> { &self.sections }
    pub fn packs(&self) -> &[Arc<NormalizedMaterialCardPack>] { &self.packs }

    /// Bind the complete selected scale. Installed length/tension, hammer,
    /// unison and duplex geometry survive; μ, EA, EI and ηI use source integrals.
    /// Requiring complete coverage makes the damping policy unambiguous.
    pub fn apply(&self, courses: &[Course]) -> Result<Vec<Course>, String> {
        let keys = courses.iter().map(|course| course.midi).collect::<BTreeSet<_>>();
        if keys.len() != courses.len() || keys != self.sections.keys().copied().collect() {
            return Err("string materials must cover exactly every selected scale key".into());
        }
        courses.iter().map(|course| course.with_bonded_kelvin_voigt_section(
            &self.sections[&course.midi])).collect()
    }

    pub(crate) fn verify_courses(&self, courses: &[Course]) -> Result<(), String> {
        if self.apply(courses)?.as_slice() != courses {
            return Err("live string coefficients differ from their bound material sections".into());
        }
        Ok(())
    }

    /// Actual source/state identities and coefficients installed in playback.
    pub fn json(&self) -> String {
        let sections = self.sections.iter().map(|(key, section)| {
            let elastic = section.elastic();
            let identities = elastic.layers().iter().map(|layer|
                format!("\"{}\"", layer.material.identity().to_hex()))
                .collect::<Vec<_>>().join(",");
            let radii = elastic.layers().iter().map(|layer| layer.outer_radius_m)
                .collect::<Vec<_>>();
            let band = section.omega_band_rad_s();
            format!(concat!("{{\"key\":{},\"stateIdentities\":[{}],\"outerRadiiM\":{:?},",
                "\"lineDensityKgM\":{},\"axialRigidityN\":{},\"flexuralRigidityNm2\":{},",
                "\"viscousRigidityNm2S\":{},\"omegaBandRadS\":[{},{}]}}"), key, identities, radii,
                elastic.line_density_kg_m(), elastic.axial_stiffness_n(),
                elastic.bending_stiffness_n_m2(), section.viscous_stiffness_n_m2_s(), band.0, band.1)
        }).collect::<Vec<_>>().join(",");
        format!("{{\"model\":\"bonded-concentric-kelvin-voigt\",\"importedPacks\":{},\"sections\":[{}]}}",
            self.packs.len(), sections)
    }
}
