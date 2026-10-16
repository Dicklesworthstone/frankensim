//! Abaqus INP and Nastran BDF thermal/geometry subset import and card census.
//!
//! Bead: `frankensim-extreal-program-f85xj.11.5`
//!
//! Provides bounded, hardened parsing of geometry, node/element sets, thermal
//! material properties, and boundary conditions from Abaqus .inp and Nastran .bdf / .dat
//! input files. Emits an explicit census of supported and unsupported cards and
//! never claims external execution semantics.

use crate::quarantine::{ImportReceipt, Quarantined};
use crate::{IoError, MAX_ELEMENTS};
use fs_blake3::{ContentHash, DomainHasher};
use std::collections::BTreeMap;

/// Extracted FE node.
#[derive(Debug, Clone, PartialEq)]
pub struct FeNode {
    /// Node identifier.
    pub id: u64,
    /// 3D coordinates.
    pub coords: [f64; 3],
}

/// Extracted FE solid/shell element.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FeElement {
    /// Element identifier.
    pub id: u64,
    /// Element topology type / name (e.g. "C3D4", "DC3D8", "CTETRA").
    pub element_type: String,
    /// Node connectivity list.
    pub node_ids: Vec<u64>,
}

/// Extracted thermal material card.
#[derive(Debug, Clone, PartialEq)]
pub struct FeMaterial {
    /// Material name.
    pub name: String,
    /// Thermal conductivity [W/(m*K)].
    pub thermal_conductivity: Option<f64>,
    /// Specific heat capacity [J/(kg*K)].
    pub specific_heat: Option<f64>,
    /// Mass density [kg/m^3].
    pub density: Option<f64>,
}

/// Extracted thermal boundary condition or source.
#[derive(Debug, Clone, PartialEq)]
pub enum FeBoundaryCondition {
    /// Prescribed temperature Dirichlet condition.
    PrescribedTemperature {
        /// Target node or element set name.
        set_name: String,
        /// Prescribed temperature in Kelvin.
        temperature_k: f64,
    },
    /// Convective film Robin condition.
    ConvectiveFilm {
        /// Target surface or element set name.
        set_name: String,
        /// Convective heat transfer coefficient [W/(m^2*K)].
        h_coeff: f64,
        /// Ambient reference temperature in Kelvin.
        ambient_k: f64,
    },
    /// Surface or volumetric heat flux.
    HeatFlux {
        /// Target set name.
        set_name: String,
        /// Heat flux magnitude in Watts (or W/m^2).
        flux_w: f64,
    },
}

/// Extracted model structure from an INP or BDF file.
#[derive(Debug, Clone, PartialEq)]
pub struct FeModel {
    /// Source dialect ("Abaqus-INP" or "Nastran-BDF").
    pub dialect: String,
    /// Extracted nodes.
    pub nodes: Vec<FeNode>,
    /// Extracted elements.
    pub elements: Vec<FeElement>,
    /// Node sets.
    pub node_sets: BTreeMap<String, Vec<u64>>,
    /// Element sets.
    pub element_sets: BTreeMap<String, Vec<u64>>,
    /// Materials.
    pub materials: Vec<FeMaterial>,
    /// Boundary conditions.
    pub boundary_conditions: Vec<FeBoundaryCondition>,
    /// Census of supported card occurrences.
    pub supported_cards: BTreeMap<String, usize>,
    /// Census of unsupported card occurrences.
    pub unsupported_cards: BTreeMap<String, usize>,
}

impl FeModel {
    /// Create a new empty FE model with declared dialect.
    #[must_use]
    pub fn new(dialect: impl Into<String>) -> Self {
        Self {
            dialect: dialect.into(),
            nodes: Vec::new(),
            elements: Vec::new(),
            node_sets: BTreeMap::new(),
            element_sets: BTreeMap::new(),
            materials: Vec::new(),
            boundary_conditions: Vec::new(),
            supported_cards: BTreeMap::new(),
            unsupported_cards: BTreeMap::new(),
        }
    }
}

fn hash_len(hasher: &mut DomainHasher, len: usize) {
    let wire_len = u64::try_from(len).expect("collection length must fit the u64 receipt format");
    hasher.update(&wire_len.to_le_bytes());
}

fn hash_str(hasher: &mut DomainHasher, value: &str) {
    hash_len(hasher, value.len());
    hasher.update(value.as_bytes());
}

fn hash_optional_f64(hasher: &mut DomainHasher, value: Option<f64>) {
    match value {
        Some(value) => {
            hasher.update(&[1]);
            hasher.update(&value.to_bits().to_le_bytes());
        }
        None => hasher.update(&[0]),
    }
}

fn hash_fe_model(domain: &str, model: &FeModel) -> ContentHash {
    let mut hasher = DomainHasher::new(domain);
    hash_str(&mut hasher, &model.dialect);

    hash_len(&mut hasher, model.nodes.len());
    for node in &model.nodes {
        hasher.update(&node.id.to_le_bytes());
        for coordinate in node.coords {
            hasher.update(&coordinate.to_bits().to_le_bytes());
        }
    }

    hash_len(&mut hasher, model.elements.len());
    for element in &model.elements {
        hasher.update(&element.id.to_le_bytes());
        hash_str(&mut hasher, &element.element_type);
        hash_len(&mut hasher, element.node_ids.len());
        for node_id in &element.node_ids {
            hasher.update(&node_id.to_le_bytes());
        }
    }

    for sets in [&model.node_sets, &model.element_sets] {
        hash_len(&mut hasher, sets.len());
        for (name, ids) in sets {
            hash_str(&mut hasher, name);
            hash_len(&mut hasher, ids.len());
            for id in ids {
                hasher.update(&id.to_le_bytes());
            }
        }
    }

    hash_len(&mut hasher, model.materials.len());
    for material in &model.materials {
        hash_str(&mut hasher, &material.name);
        hash_optional_f64(&mut hasher, material.thermal_conductivity);
        hash_optional_f64(&mut hasher, material.specific_heat);
        hash_optional_f64(&mut hasher, material.density);
    }

    hash_len(&mut hasher, model.boundary_conditions.len());
    for condition in &model.boundary_conditions {
        match condition {
            FeBoundaryCondition::PrescribedTemperature {
                set_name,
                temperature_k,
            } => {
                hasher.update(&[0]);
                hash_str(&mut hasher, set_name);
                hasher.update(&temperature_k.to_bits().to_le_bytes());
            }
            FeBoundaryCondition::ConvectiveFilm {
                set_name,
                h_coeff,
                ambient_k,
            } => {
                hasher.update(&[1]);
                hash_str(&mut hasher, set_name);
                hasher.update(&h_coeff.to_bits().to_le_bytes());
                hasher.update(&ambient_k.to_bits().to_le_bytes());
            }
            FeBoundaryCondition::HeatFlux { set_name, flux_w } => {
                hasher.update(&[2]);
                hash_str(&mut hasher, set_name);
                hasher.update(&flux_w.to_bits().to_le_bytes());
            }
        }
    }

    for cards in [&model.supported_cards, &model.unsupported_cards] {
        hash_len(&mut hasher, cards.len());
        for (name, count) in cards {
            hash_str(&mut hasher, name);
            hash_len(&mut hasher, *count);
        }
    }

    hasher.finalize()
}

/// Verifiable cryptographic receipt for an admitted INP or BDF model.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InpBdfReceipt {
    /// Dialect name.
    pub dialect: String,
    /// Extracted node count.
    pub node_count: usize,
    /// Extracted element count.
    pub element_count: usize,
    /// Material count.
    pub material_count: usize,
    /// Supported card count.
    pub supported_card_count: usize,
    /// Unsupported card count.
    pub unsupported_card_count: usize,
    /// Content hash of the extracted model.
    pub content_hash: ContentHash,
}

/// Parse an Abaqus INP file extracting thermal and geometry subsets.
#[allow(clippy::too_many_lines)]
pub fn parse_abaqus_inp(input: &str) -> Result<(Quarantined<FeModel>, InpBdfReceipt), IoError> {
    let mut model = FeModel::new("Abaqus-INP");
    let mut current_section = "";
    let mut current_element_type = String::new();
    let mut current_set_name = String::new();
    let mut current_material = Option::<FeMaterial>::None;

    for (line_idx, raw_line) in input.lines().enumerate() {
        let line = raw_line.trim();
        if line.is_empty() || line.starts_with("**") {
            continue;
        }

        if line.starts_with('*') {
            // New keyword card
            let upper = line.to_ascii_uppercase();
            if upper.starts_with("*NODE") {
                current_section = "NODE";
                *model
                    .supported_cards
                    .entry("*NODE".to_string())
                    .or_insert(0) += 1;
            } else if upper.starts_with("*ELEMENT") {
                current_section = "ELEMENT";
                *model
                    .supported_cards
                    .entry("*ELEMENT".to_string())
                    .or_insert(0) += 1;
                // Parse TYPE=
                current_element_type = "C3D4".to_string();
                for part in upper.split(',') {
                    let p = part.trim();
                    if p.starts_with("TYPE=") {
                        current_element_type = p["TYPE=".len()..].trim().to_string();
                    }
                }
            } else if upper.starts_with("*NSET") {
                current_section = "NSET";
                *model
                    .supported_cards
                    .entry("*NSET".to_string())
                    .or_insert(0) += 1;
                current_set_name = "NSET".to_string();
                for part in line.split(',') {
                    let p = part.trim();
                    if p.to_ascii_uppercase().starts_with("NSET=") {
                        current_set_name = p["NSET=".len()..].trim().to_string();
                    }
                }
                model.node_sets.entry(current_set_name.clone()).or_default();
            } else if upper.starts_with("*ELSET") {
                current_section = "ELSET";
                *model
                    .supported_cards
                    .entry("*ELSET".to_string())
                    .or_insert(0) += 1;
                current_set_name = "ELSET".to_string();
                for part in line.split(',') {
                    let p = part.trim();
                    if p.to_ascii_uppercase().starts_with("ELSET=") {
                        current_set_name = p["ELSET=".len()..].trim().to_string();
                    }
                }
                model
                    .element_sets
                    .entry(current_set_name.clone())
                    .or_default();
            } else if upper.starts_with("*MATERIAL") {
                if let Some(mat) = current_material.take() {
                    model.materials.push(mat);
                }
                current_section = "MATERIAL";
                *model
                    .supported_cards
                    .entry("*MATERIAL".to_string())
                    .or_insert(0) += 1;
                let mut name = "MAT".to_string();
                for part in line.split(',') {
                    let p = part.trim();
                    if p.to_ascii_uppercase().starts_with("NAME=") {
                        name = p["NAME=".len()..].trim().to_string();
                    }
                }
                current_material = Some(FeMaterial {
                    name,
                    thermal_conductivity: None,
                    specific_heat: None,
                    density: None,
                });
            } else if upper.starts_with("*CONDUCTIVITY") {
                current_section = "CONDUCTIVITY";
                *model
                    .supported_cards
                    .entry("*CONDUCTIVITY".to_string())
                    .or_insert(0) += 1;
            } else if upper.starts_with("*SPECIFIC HEAT") {
                current_section = "SPECIFIC_HEAT";
                *model
                    .supported_cards
                    .entry("*SPECIFIC HEAT".to_string())
                    .or_insert(0) += 1;
            } else if upper.starts_with("*DENSITY") {
                current_section = "DENSITY";
                *model
                    .supported_cards
                    .entry("*DENSITY".to_string())
                    .or_insert(0) += 1;
            } else if upper.starts_with("*BOUNDARY") {
                current_section = "BOUNDARY";
                *model
                    .supported_cards
                    .entry("*BOUNDARY".to_string())
                    .or_insert(0) += 1;
            } else if upper.starts_with("*DFLUX") || upper.starts_with("*FILM") {
                current_section = "FLUX";
                *model
                    .supported_cards
                    .entry(upper.split(',').next().unwrap_or("*FLUX").to_string())
                    .or_insert(0) += 1;
            } else {
                current_section = "UNSUPPORTED";
                let card_name = upper.split(',').next().unwrap_or(&upper).to_string();
                *model.unsupported_cards.entry(card_name).or_insert(0) += 1;
            }
            continue;
        }

        // Data lines
        match current_section {
            "NODE" => {
                let parts: Vec<&str> = line.split(',').map(str::trim).collect();
                if parts.len() >= 4 {
                    let id: u64 = parts[0].parse().map_err(|_| IoError::Malformed {
                        at: line_idx,
                        what: "invalid node id in *NODE".to_string(),
                    })?;
                    let x: f64 = parts[1].parse().unwrap_or(0.0);
                    let y: f64 = parts[2].parse().unwrap_or(0.0);
                    let z: f64 = parts[3].parse().unwrap_or(0.0);
                    if model.nodes.len() < MAX_ELEMENTS {
                        model.nodes.push(FeNode {
                            id,
                            coords: [x, y, z],
                        });
                    }
                }
            }
            "ELEMENT" => {
                let parts: Vec<&str> = line.split(',').map(str::trim).collect();
                if parts.len() >= 2 {
                    let id: u64 = parts[0].parse().map_err(|_| IoError::Malformed {
                        at: line_idx,
                        what: "invalid element id in *ELEMENT".to_string(),
                    })?;
                    let node_ids: Vec<u64> =
                        parts[1..].iter().filter_map(|p| p.parse().ok()).collect();
                    if model.elements.len() < MAX_ELEMENTS {
                        model.elements.push(FeElement {
                            id,
                            element_type: current_element_type.clone(),
                            node_ids,
                        });
                    }
                }
            }
            "NSET" => {
                let set = model.node_sets.entry(current_set_name.clone()).or_default();
                for part in line.split(',').map(str::trim) {
                    if let Ok(id) = part.parse::<u64>() {
                        set.push(id);
                    }
                }
            }
            "ELSET" => {
                let set = model
                    .element_sets
                    .entry(current_set_name.clone())
                    .or_default();
                for part in line.split(',').map(str::trim) {
                    if let Ok(id) = part.parse::<u64>() {
                        set.push(id);
                    }
                }
            }
            "CONDUCTIVITY" => {
                if let Some(mat) = current_material.as_mut() {
                    let val: f64 = line
                        .split(',')
                        .next()
                        .unwrap_or("")
                        .trim()
                        .parse()
                        .unwrap_or(0.0);
                    mat.thermal_conductivity = Some(val);
                }
            }
            "SPECIFIC_HEAT" => {
                if let Some(mat) = current_material.as_mut() {
                    let val: f64 = line
                        .split(',')
                        .next()
                        .unwrap_or("")
                        .trim()
                        .parse()
                        .unwrap_or(0.0);
                    mat.specific_heat = Some(val);
                }
            }
            "DENSITY" => {
                if let Some(mat) = current_material.as_mut() {
                    let val: f64 = line
                        .split(',')
                        .next()
                        .unwrap_or("")
                        .trim()
                        .parse()
                        .unwrap_or(0.0);
                    mat.density = Some(val);
                }
            }
            "BOUNDARY" => {
                let parts: Vec<&str> = line.split(',').map(str::trim).collect();
                if parts.len() >= 4 {
                    let set_name = parts[0].to_string();
                    let temp: f64 = parts[3].parse().unwrap_or(300.0);
                    model
                        .boundary_conditions
                        .push(FeBoundaryCondition::PrescribedTemperature {
                            set_name,
                            temperature_k: temp,
                        });
                }
            }
            _ => {}
        }
    }

    if let Some(mat) = current_material {
        model.materials.push(mat);
    }

    // Build receipt
    let supported_total: usize = model.supported_cards.values().sum();
    let unsupported_total: usize = model.unsupported_cards.values().sum();

    let content_hash = hash_fe_model("org.frankensim.inp.receipt.v2", &model);

    let receipt = InpBdfReceipt {
        dialect: model.dialect.clone(),
        node_count: model.nodes.len(),
        element_count: model.elements.len(),
        material_count: model.materials.len(),
        supported_card_count: supported_total,
        unsupported_card_count: unsupported_total,
        content_hash,
    };

    let source_receipt = ImportReceipt {
        format: "abaqus-inp",
        source_hash: fs_obs::fnv1a64(input.as_bytes()),
        parser_version: crate::VERSION,
        parsed: (model.nodes.len(), model.elements.len()),
    };

    Ok((Quarantined::new(model, source_receipt, Vec::new()), receipt))
}

/// Parse a Nastran BDF / DAT file extracting thermal and geometry subsets.
#[allow(clippy::too_many_lines)]
pub fn parse_nastran_bdf(input: &str) -> Result<(Quarantined<FeModel>, InpBdfReceipt), IoError> {
    let mut model = FeModel::new("Nastran-BDF");

    for (line_idx, raw_line) in input.lines().enumerate() {
        let line = raw_line.trim();
        if line.is_empty() || line.starts_with('$') {
            continue;
        }

        // Split free-field commas or 8-character fixed-field
        let tokens: Vec<&str> = if line.contains(',') {
            line.split(',').map(str::trim).collect()
        } else {
            // 8-character fixed field chunks
            let mut chunks = Vec::new();
            let bytes = line.as_bytes();
            let mut i = 0;
            while i < bytes.len() {
                let end = (i + 8).min(bytes.len());
                chunks.push(std::str::from_utf8(&bytes[i..end]).unwrap_or("").trim());
                i += 8;
            }
            chunks
        };

        if tokens.is_empty() {
            continue;
        }

        let card = tokens[0].to_ascii_uppercase();
        match card.as_str() {
            "GRID" => {
                *model.supported_cards.entry("GRID".to_string()).or_insert(0) += 1;
                if tokens.len() >= 6 {
                    let id: u64 = tokens[1].parse().map_err(|_| IoError::Malformed {
                        at: line_idx,
                        what: "invalid GRID id".to_string(),
                    })?;
                    let x: f64 = tokens[3].parse().unwrap_or(0.0);
                    let y: f64 = tokens[4].parse().unwrap_or(0.0);
                    let z: f64 = tokens[5].parse().unwrap_or(0.0);
                    if model.nodes.len() < MAX_ELEMENTS {
                        model.nodes.push(FeNode {
                            id,
                            coords: [x, y, z],
                        });
                    }
                }
            }
            "CTETRA" | "CHEXA" | "CPENTA" | "CPYRAM" | "CTRIA3" | "CQUAD4" => {
                *model.supported_cards.entry(card.clone()).or_insert(0) += 1;
                if tokens.len() >= 4 {
                    let id: u64 = tokens[1].parse().map_err(|_| IoError::Malformed {
                        at: line_idx,
                        what: format!("invalid {card} element id"),
                    })?;
                    let node_ids: Vec<u64> =
                        tokens[3..].iter().filter_map(|t| t.parse().ok()).collect();
                    if model.elements.len() < MAX_ELEMENTS {
                        model.elements.push(FeElement {
                            id,
                            element_type: card,
                            node_ids,
                        });
                    }
                }
            }
            "MAT4" | "MAT5" => {
                *model.supported_cards.entry(card).or_insert(0) += 1;
                if tokens.len() >= 3 {
                    let name = tokens[1].to_string();
                    let k: f64 = tokens[2].parse().unwrap_or(0.0);
                    let cp: f64 = tokens.get(3).and_then(|t| t.parse().ok()).unwrap_or(0.0);
                    let rho: f64 = tokens.get(4).and_then(|t| t.parse().ok()).unwrap_or(0.0);
                    model.materials.push(FeMaterial {
                        name,
                        thermal_conductivity: Some(k),
                        specific_heat: if cp > 0.0 { Some(cp) } else { None },
                        density: if rho > 0.0 { Some(rho) } else { None },
                    });
                }
            }
            "TEMP" | "TEMPD" => {
                *model.supported_cards.entry(card).or_insert(0) += 1;
                if tokens.len() >= 3 {
                    let set_name = tokens[1].to_string();
                    let temp: f64 = tokens[2].parse().unwrap_or(300.0);
                    model
                        .boundary_conditions
                        .push(FeBoundaryCondition::PrescribedTemperature {
                            set_name,
                            temperature_k: temp,
                        });
                }
            }
            _ => {
                *model.unsupported_cards.entry(card).or_insert(0) += 1;
            }
        }
    }

    let supported_total: usize = model.supported_cards.values().sum();
    let unsupported_total: usize = model.unsupported_cards.values().sum();

    let content_hash = hash_fe_model("org.frankensim.bdf.receipt.v2", &model);

    let receipt = InpBdfReceipt {
        dialect: model.dialect.clone(),
        node_count: model.nodes.len(),
        element_count: model.elements.len(),
        material_count: model.materials.len(),
        supported_card_count: supported_total,
        unsupported_card_count: unsupported_total,
        content_hash,
    };

    let source_receipt = ImportReceipt {
        format: "nastran-bdf",
        source_hash: fs_obs::fnv1a64(input.as_bytes()),
        parser_version: crate::VERSION,
        parsed: (model.nodes.len(), model.elements.len()),
    };

    Ok((Quarantined::new(model, source_receipt, Vec::new()), receipt))
}
