//! Private frozen sidecar codec. Length-framed UTF-8, little-endian counts,
//! finite shortest-round-trip quantities and exactly six SI exponents.
use super::*;
use crate::entity::{read_assembly_catalog, write_assembly_catalog};

const MAGIC: &[u8; 8] = b"FSASMBL\0";

pub(crate) struct Writer {
    bytes: Vec<u8>,
    budget: AssemblyBudget,
    error: Option<AssemblyError>,
}
impl Writer {
    pub(crate) fn new(budget: AssemblyBudget) -> Self {
        Self {
            bytes: Vec::new(),
            budget,
            error: None,
        }
    }
    pub(crate) fn put(&mut self, value: &[u8]) {
        if self.error.is_some() {
            return;
        }
        if value.len() > self.budget.max_wire_bytes.saturating_sub(self.bytes.len()) {
            self.error = Some(AssemblyError::BudgetExceeded {
                resource: "wire bytes",
                cap: self.budget.max_wire_bytes,
            });
            return;
        }
        if self.bytes.try_reserve(value.len()).is_err() {
            self.error = Some(invalid("allocation", "assembly wire"));
            return;
        }
        self.bytes.extend_from_slice(value);
    }
    pub(crate) fn tag(&mut self, value: u8) {
        self.put(&[value]);
    }
    pub(crate) fn integer(&mut self, value: u64) {
        self.put(&value.to_le_bytes());
    }
    pub(crate) fn count(&mut self, value: usize, limit: usize, resource: &'static str) {
        if let Err(e) = cap(value, limit, resource) {
            self.error = Some(e);
            return;
        }
        self.integer(value as u64);
    }
    pub(crate) fn text(&mut self, value: &str) {
        self.count(value.len(), self.budget.max_text_bytes, "text bytes");
        self.put(value.as_bytes());
    }
    pub(crate) fn hash(&mut self, value: ContentHash) {
        self.put(&value.0);
    }
    pub(crate) fn entity(&mut self, value: EntityId) {
        self.tag(match value.kind() {
            EntityKind::Assembly => 1,
            EntityKind::Part => 2,
            EntityKind::Region => 3,
            EntityKind::Surface => 4,
            EntityKind::Interface => 5,
        });
        self.hash(value.digest());
    }
    pub(crate) fn maybe_entity(&mut self, value: Option<EntityId>) {
        self.tag(u8::from(value.is_some()));
        if let Some(v) = value {
            self.entity(v);
        }
    }
    pub(crate) fn maybe_text(&mut self, value: Option<&str>) {
        self.tag(u8::from(value.is_some()));
        if let Some(v) = value {
            self.text(v);
        }
    }
    pub(crate) fn qty(&mut self, q: QtyAny) {
        if !q.value.is_finite() {
            self.error = Some(invalid("quantity", "non-finite wire value"));
            return;
        }
        let value = if q.value == 0.0 { 0.0 } else { q.value };
        self.text(&value.to_string());
        for exponent in q.dims.0 {
            self.tag(exponent as u8);
        }
    }
    pub(crate) fn key(&mut self, key: &ParameterKey) {
        self.text(&key.name);
        for exponent in key.dims.0 {
            self.tag(exponent as u8);
        }
    }
    pub(crate) fn source(&mut self, source: &ToleranceSource) {
        match source {
            ToleranceSource::Drawing { sheet, note } => {
                self.tag(0);
                self.text(sheet);
                self.text(note);
            }
            ToleranceSource::Standard { clause } => {
                self.tag(1);
                self.text(clause);
            }
            ToleranceSource::Assumed { rationale } => {
                self.tag(2);
                self.text(rationale);
            }
            ToleranceSource::Measured { ledger_row, method } => {
                self.tag(3);
                self.text(ledger_row);
                self.text(method);
            }
            ToleranceSource::Reconstructed { ledger_row, route } => {
                self.tag(4);
                self.text(ledger_row);
                self.text(route);
            }
        }
    }
    fn bound(&mut self, b: &DimensionBound) {
        self.key(&b.property);
        self.qty(b.lower);
        self.qty(b.upper);
        self.source(&b.source);
    }
    fn parameters(&mut self, p: &[(ParameterKey, QtyAny)]) {
        self.count(p.len(), self.budget.max_nested_items, "parameters");
        for (k, q) in p {
            self.key(k);
            self.qty(*q);
        }
    }
    pub(crate) fn finish(self) -> Result<Vec<u8>, AssemblyError> {
        self.error.map_or(Ok(self.bytes), Err)
    }
}

pub(crate) struct Reader<'a> {
    bytes: &'a [u8],
    offset: usize,
    pub(crate) budget: AssemblyBudget,
}
impl<'a> Reader<'a> {
    pub(crate) fn new(bytes: &'a [u8], budget: AssemblyBudget) -> Result<Self, AssemblyError> {
        cap(bytes.len(), budget.max_wire_bytes, "wire bytes")?;
        Ok(Self {
            bytes,
            offset: 0,
            budget,
        })
    }
    pub(crate) fn fail(&self, detail: impl Into<String>) -> AssemblyError {
        AssemblyError::Wire {
            offset: self.offset,
            detail: detail.into(),
        }
    }
    pub(crate) fn take(&mut self, n: usize) -> Result<&'a [u8], AssemblyError> {
        let end = self
            .offset
            .checked_add(n)
            .ok_or_else(|| self.fail("byte offset overflow"))?;
        let result = self
            .bytes
            .get(self.offset..end)
            .ok_or_else(|| self.fail("truncated record"))?;
        self.offset = end;
        Ok(result)
    }
    pub(crate) fn tag(&mut self) -> Result<u8, AssemblyError> {
        Ok(self.take(1)?[0])
    }
    pub(crate) fn boolean(&mut self) -> Result<bool, AssemblyError> {
        match self.tag()? {
            0 => Ok(false),
            1 => Ok(true),
            _ => Err(self.fail("invalid option tag")),
        }
    }
    pub(crate) fn integer(&mut self) -> Result<u64, AssemblyError> {
        Ok(u64::from_le_bytes(
            self.take(8)?
                .try_into()
                .map_err(|_| self.fail("integer width"))?,
        ))
    }
    pub(crate) fn count(
        &mut self,
        limit: usize,
        resource: &'static str,
    ) -> Result<usize, AssemblyError> {
        let value =
            usize::try_from(self.integer()?).map_err(|_| AssemblyError::BudgetExceeded {
                resource,
                cap: limit,
            })?;
        cap(value, limit, resource)?;
        if value > self.bytes.len().saturating_sub(self.offset) {
            return Err(self.fail("count exceeds remaining bytes"));
        }
        Ok(value)
    }
    pub(crate) fn text(&mut self) -> Result<String, AssemblyError> {
        let len = self.count(self.budget.max_text_bytes, "text bytes")?;
        let bytes = self.take(len)?;
        let value = std::str::from_utf8(bytes).map_err(|_| self.fail("invalid UTF-8"))?;
        Ok(value.to_owned())
    }
    pub(crate) fn hash(&mut self) -> Result<ContentHash, AssemblyError> {
        Ok(ContentHash(
            self.take(32)?
                .try_into()
                .map_err(|_| self.fail("hash width"))?,
        ))
    }
    pub(crate) fn entity(&mut self) -> Result<EntityId, AssemblyError> {
        let kind = match self.tag()? {
            1 => EntityKind::Assembly,
            2 => EntityKind::Part,
            3 => EntityKind::Region,
            4 => EntityKind::Surface,
            5 => EntityKind::Interface,
            _ => return Err(self.fail("unknown entity kind")),
        };
        Ok(EntityId::from_wire(kind, self.hash()?))
    }
    pub(crate) fn maybe_entity(&mut self) -> Result<Option<EntityId>, AssemblyError> {
        if self.boolean()? {
            Ok(Some(self.entity()?))
        } else {
            Ok(None)
        }
    }
    pub(crate) fn maybe_text(&mut self) -> Result<Option<String>, AssemblyError> {
        if self.boolean()? {
            Ok(Some(self.text()?))
        } else {
            Ok(None)
        }
    }
    pub(crate) fn dims(&mut self) -> Result<Dims, AssemblyError> {
        let b = self.take(6)?;
        Ok(Dims([
            b[0] as i8, b[1] as i8, b[2] as i8, b[3] as i8, b[4] as i8, b[5] as i8,
        ]))
    }
    pub(crate) fn qty(&mut self) -> Result<QtyAny, AssemblyError> {
        let text = self.text()?;
        let v: f64 = text.parse().map_err(|_| self.fail("invalid quantity"))?;
        if !v.is_finite() || (if v == 0.0 { 0.0 } else { v }).to_string() != text {
            return Err(self.fail("noncanonical quantity"));
        }
        Ok(QtyAny::new(v, self.dims()?))
    }
    pub(crate) fn key(&mut self) -> Result<ParameterKey, AssemblyError> {
        Ok(ParameterKey::new(self.text()?, self.dims()?))
    }
    pub(crate) fn source(&mut self) -> Result<ToleranceSource, AssemblyError> {
        Ok(match self.tag()? {
            0 => ToleranceSource::Drawing {
                sheet: self.text()?,
                note: self.text()?,
            },
            1 => ToleranceSource::Standard {
                clause: self.text()?,
            },
            2 => ToleranceSource::Assumed {
                rationale: self.text()?,
            },
            3 => ToleranceSource::Measured {
                ledger_row: self.text()?,
                method: self.text()?,
            },
            4 => ToleranceSource::Reconstructed {
                ledger_row: self.text()?,
                route: self.text()?,
            },
            _ => return Err(self.fail("unknown tolerance source")),
        })
    }
    fn bound(&mut self) -> Result<DimensionBound, AssemblyError> {
        Ok(DimensionBound {
            property: self.key()?,
            lower: self.qty()?,
            upper: self.qty()?,
            source: self.source()?,
        })
    }
    fn parameters(&mut self) -> Result<Vec<(ParameterKey, QtyAny)>, AssemblyError> {
        let n = self.count(self.budget.max_nested_items, "parameters")?;
        let mut result = Vec::new();
        for _ in 0..n {
            result.push((self.key()?, self.qty()?));
        }
        Ok(result)
    }
}

pub(super) fn encode(
    value: &AssemblyDescription,
    b: AssemblyBudget,
) -> Result<Vec<u8>, AssemblyError> {
    let mut w = Writer::new(b);
    w.put(MAGIC);
    w.put(&ASSEMBLY_DESCRIPTION_SCHEMA_VERSION.to_le_bytes());
    w.hash(value.scenario_ir_hash);
    w.entity(value.input.root);
    write_assembly_catalog(&value.input.catalog, &mut w, b);
    w.tag(
        match value
            .input
            .reference_state
            .ok_or(AssemblyError::MissingReferenceState)?
        {
            ReferenceState::Unloaded => 0,
            ReferenceState::Strung => 1,
            ReferenceState::Assembled => 2,
            ReferenceState::Played => 3,
        },
    );
    w.count(value.input.parts.len(), b.max_parts, "parts");
    for p in &value.input.parts {
        w.entity(p.entity);
        w.tag(match p.kind {
            PartKind::String => 0,
            PartKind::Shell => 1,
            PartKind::Solid => 2,
            PartKind::Duct => 3,
            PartKind::Cavity => 4,
            PartKind::ValvePattern => 5,
            PartKind::DrivenBody => 6,
        });
        match &p.geometry {
            GeometrySource::Parametric { generator, params } => {
                w.tag(0);
                w.text(generator);
                w.parameters(params);
            }
            GeometrySource::Imported {
                fingerprint,
                representation,
                ledger_row,
            } => {
                w.tag(1);
                w.hash(fingerprint.hash());
                w.tag(match representation {
                    Representation::Mesh => 0,
                    Representation::Nurbs => 1,
                    Representation::Sdf => 2,
                    Representation::Voxel => 3,
                });
                w.text(ledger_row);
            }
        }
        w.count(p.materials.len(), b.max_nested_items, "materials");
        for m in &p.materials {
            w.hash(m.card);
            w.parameters(&m.query_point);
            w.integer(u64::from(m.frame.0));
            w.tag(u8::from(m.layer.is_some()));
            if let Some(layer) = m.layer {
                w.integer(u64::from(layer));
            }
            w.count(m.unknown.len(), b.max_nested_items, "unknown properties");
            for UnknownProperty::Unknown { property, bound } in &m.unknown {
                w.key(property);
                w.tag(u8::from(bound.is_some()));
                if let Some(bound) = bound {
                    w.bound(bound);
                }
            }
        }
        w.count(
            p.uncertainty.len(),
            b.max_nested_items,
            "uncertainty bounds",
        );
        for bound in &p.uncertainty {
            w.bound(bound);
        }
    }
    w.count(value.input.interfaces.len(), b.max_interfaces, "interfaces");
    for i in &value.input.interfaces {
        w.entity(i.entity);
        match i.kind {
            InterfaceKind::PowerPort { schema } => {
                w.tag(0);
                w.hash(schema);
            }
            InterfaceKind::Contact { law } => {
                w.tag(1);
                w.hash(law);
            }
            InterfaceKind::Glue { law } => {
                w.tag(2);
                w.hash(law);
            }
        }
    }
    w.count(
        value.input.configurations.len(),
        b.max_configurations,
        "configurations",
    );
    for c in &value.input.configurations {
        w.text(&c.name);
        w.count(
            c.assignments.len(),
            b.max_assignments_per_state,
            "assignments per state",
        );
        for (p, q) in &c.assignments {
            w.entity(p.entity);
            w.key(&p.parameter);
            w.qty(*q);
        }
        match &c.provenance {
            StateProvenance::Declared { rationale } => {
                w.tag(0);
                w.text(rationale);
            }
            StateProvenance::PlayingMeasured { ledger_row } => {
                w.tag(1);
                w.text(ledger_row);
            }
            StateProvenance::SpeechImaging {
                ledger_row,
                subject,
                token,
                playing_envelope,
            } => {
                w.tag(2);
                w.text(ledger_row);
                w.text(subject);
                w.text(token);
                w.tag(u8::from(playing_envelope.is_some()));
                if let Some(h) = playing_envelope {
                    w.hash(*h);
                }
            }
        }
    }
    w.count(value.input.observers.len(), b.max_nested_items, "observers");
    for o in &value.input.observers {
        w.entity(o.entity);
        w.integer(u64::from(o.frame.0));
        w.key(&o.quantity);
        for p in o.position {
            w.qty(p);
        }
    }
    w.finish()
}

fn u32_value(r: &mut Reader<'_>) -> Result<u32, AssemblyError> {
    u32::try_from(r.integer()?).map_err(|_| r.fail("u32 overflow"))
}
pub(super) fn decode(
    bytes: &[u8],
    scenario: &Scenario,
    b: AssemblyBudget,
) -> Result<AssemblyDescription, AssemblyError> {
    let mut r = Reader::new(bytes, b)?;
    if r.take(8)? != MAGIC {
        return Err(r.fail("assembly magic"));
    }
    let version = u32::from_le_bytes(r.take(4)?.try_into().map_err(|_| r.fail("version width"))?);
    if version != ASSEMBLY_DESCRIPTION_SCHEMA_VERSION {
        return Err(AssemblyError::UnknownVersion { found: version });
    }
    let scenario_ir_hash = r.hash()?;
    let root = r.entity()?;
    let catalog = read_assembly_catalog(&mut r)?;
    let reference_state = Some(match r.tag()? {
        0 => ReferenceState::Unloaded,
        1 => ReferenceState::Strung,
        2 => ReferenceState::Assembled,
        3 => ReferenceState::Played,
        _ => return Err(r.fail("reference state")),
    });
    let mut parts = Vec::new();
    let n = r.count(b.max_parts, "parts")?;
    for _ in 0..n {
        let entity = r.entity()?;
        let kind = match r.tag()? {
            0 => PartKind::String,
            1 => PartKind::Shell,
            2 => PartKind::Solid,
            3 => PartKind::Duct,
            4 => PartKind::Cavity,
            5 => PartKind::ValvePattern,
            6 => PartKind::DrivenBody,
            _ => return Err(r.fail("part kind")),
        };
        let geometry = match r.tag()? {
            0 => GeometrySource::Parametric {
                generator: r.text()?,
                params: r.parameters()?,
            },
            1 => {
                let fingerprint = GeometryFingerprint::from_hash(r.hash()?);
                let representation = match r.tag()? {
                    0 => Representation::Mesh,
                    1 => Representation::Nurbs,
                    2 => Representation::Sdf,
                    3 => Representation::Voxel,
                    _ => return Err(r.fail("representation")),
                };
                GeometrySource::Imported {
                    fingerprint,
                    representation,
                    ledger_row: r.text()?,
                }
            }
            _ => return Err(r.fail("geometry source")),
        };
        let mut materials = Vec::new();
        let n = r.count(b.max_nested_items, "materials")?;
        for _ in 0..n {
            let card = r.hash()?;
            let query_point = r.parameters()?;
            let frame = FrameId(u32_value(&mut r)?);
            let layer = if r.boolean()? {
                Some(u32_value(&mut r)?)
            } else {
                None
            };
            let mut unknown = Vec::new();
            let n = r.count(b.max_nested_items, "unknown properties")?;
            for _ in 0..n {
                let property = r.key()?;
                let bound = if r.boolean()? { Some(r.bound()?) } else { None };
                unknown.push(UnknownProperty::Unknown { property, bound });
            }
            materials.push(MaterialRef {
                card,
                query_point,
                frame,
                layer,
                unknown,
            });
        }
        let mut uncertainty = Vec::new();
        let n = r.count(b.max_nested_items, "uncertainty bounds")?;
        for _ in 0..n {
            uncertainty.push(r.bound()?);
        }
        parts.push(PartDecl {
            entity,
            kind,
            geometry,
            materials,
            uncertainty,
        });
    }
    let mut interfaces = Vec::new();
    let n = r.count(b.max_interfaces, "interfaces")?;
    for _ in 0..n {
        let entity = r.entity()?;
        let kind = match r.tag()? {
            0 => InterfaceKind::PowerPort { schema: r.hash()? },
            1 => InterfaceKind::Contact { law: r.hash()? },
            2 => InterfaceKind::Glue { law: r.hash()? },
            _ => return Err(r.fail("interface kind")),
        };
        interfaces.push(InterfaceDecl { entity, kind });
    }
    let mut configurations = Vec::new();
    let n = r.count(b.max_configurations, "configurations")?;
    for _ in 0..n {
        let name = r.text()?;
        let mut assignments = Vec::new();
        let n = r.count(b.max_assignments_per_state, "assignments per state")?;
        for _ in 0..n {
            assignments.push((
                ParameterPath {
                    entity: r.entity()?,
                    parameter: r.key()?,
                },
                r.qty()?,
            ));
        }
        let provenance = match r.tag()? {
            0 => StateProvenance::Declared {
                rationale: r.text()?,
            },
            1 => StateProvenance::PlayingMeasured {
                ledger_row: r.text()?,
            },
            2 => {
                let ledger_row = r.text()?;
                let subject = r.text()?;
                let token = r.text()?;
                let playing_envelope = if r.boolean()? { Some(r.hash()?) } else { None };
                StateProvenance::SpeechImaging {
                    ledger_row,
                    subject,
                    token,
                    playing_envelope,
                }
            }
            _ => return Err(r.fail("state provenance")),
        };
        configurations.push(ConfigurationState {
            name,
            assignments,
            provenance,
        });
    }
    let mut observers = Vec::new();
    let n = r.count(b.max_nested_items, "observers")?;
    for _ in 0..n {
        observers.push(ObserverDecl {
            entity: r.entity()?,
            frame: FrameId(u32_value(&mut r)?),
            quantity: r.key()?,
            position: [r.qty()?, r.qty()?, r.qty()?],
        });
    }
    if r.offset != bytes.len() {
        return Err(AssemblyError::TrailingBytes);
    }
    let result = AssemblyDescription {
        scenario_ir_hash,
        input: AssemblyInput {
            root,
            catalog,
            parts,
            interfaces,
            configurations,
            observers,
            reference_state,
        },
    };
    result.verify_scenario(scenario)?;
    validate(&result.input, scenario, b)?;
    if encode(&result, b)? != bytes {
        return Err(r.fail("noncanonical assembly encoding"));
    }
    Ok(result)
}
