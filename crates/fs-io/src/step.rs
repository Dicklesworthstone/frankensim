//! Bounded ISO 10303-21 clear-text structure parsing and canonical writing.
//!
//! This module deliberately stops at the exchange-file syntax boundary. It
//! validates the Part-21 envelope, mandatory header records, entity-instance
//! graph, parameter syntax, resource limits, duplicate identifiers, and
//! dangling references. It does **not** interpret an EXPRESS schema, tessellate
//! CAD geometry, establish AP conformance, or confer geometric authority.

use crate::IoError;
use std::collections::BTreeSet;
use std::fmt::Write as _;

/// Receipt/parser schema for this bounded syntax representation.
pub const STEP_SYNTAX_VERSION: &str = "part21-syntax-v1";

const HARD_MAX_STEP_NESTING: usize = 256;
const HARD_MAX_STEP_STRING_BYTES: usize = 32_769;
const HARD_MAX_STEP_IDENTIFIER_BYTES: usize = 255;

/// Conservative resource limits for the syntax-only Part-21 kernel.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct StepLimits {
    /// Maximum source or canonical-output byte length.
    pub max_bytes: usize,
    /// Maximum number of lexical tokens.
    pub max_tokens: usize,
    /// Maximum number of DATA-section instances.
    pub max_instances: usize,
    /// Maximum total number of parameter values, including nested values.
    pub max_values: usize,
    /// Maximum aggregate/typed-parameter nesting depth.
    pub max_nesting: usize,
    /// Maximum encoded string-token length, including delimiters.
    pub max_string_bytes: usize,
    /// Maximum decimal number-token length.
    pub max_number_bytes: usize,
    /// Maximum identifier or enumeration length in bytes.
    pub max_identifier_bytes: usize,
    /// Maximum number of schema identifiers in `FILE_SCHEMA`.
    pub max_schema_ids: usize,
    /// Maximum simple-entity components in one complex instance.
    pub max_components_per_instance: usize,
}

impl Default for StepLimits {
    fn default() -> Self {
        Self {
            max_bytes: 256 * 1024 * 1024,
            max_tokens: 16_000_000,
            max_instances: 1_000_000,
            max_values: 8_000_000,
            max_nesting: 128,
            max_string_bytes: 32_769,
            max_number_bytes: 4_096,
            max_identifier_bytes: 255,
            max_schema_ids: 32,
            max_components_per_instance: 1_024,
        }
    }
}

/// A non-authoritative application-protocol hint derived from FILE_SCHEMA.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StepProfileHint {
    /// A schema label names AP203 or CONFIG_CONTROL_DESIGN.
    Ap203,
    /// A schema label names AP214 or AUTOMOTIVE_DESIGN.
    Ap214,
    /// Labels point at both AP203 and AP214 families.
    Ambiguous,
    /// No supported profile label was recognized.
    Other,
}

impl StepProfileHint {
    /// Stable receipt label.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Ap203 => "ap203",
            Self::Ap214 => "ap214",
            Self::Ambiguous => "ambiguous",
            Self::Other => "other",
        }
    }
}

/// One Part-21 parameter value in the syntax-only representation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum StepValue {
    /// `$`: an omitted value.
    Omitted,
    /// `*`: a derived value.
    Derived,
    /// A decimal integer or real token, preserved without binary conversion.
    Number(String),
    /// A STEP string with doubled apostrophes decoded.
    String(String),
    /// A `.NAME.` enumeration token, without the surrounding periods.
    Enumeration(String),
    /// A `#N` entity reference.
    Reference(u64),
    /// A parenthesized aggregate.
    Aggregate(Vec<StepValue>),
    /// A typed parameter such as `LENGTH_MEASURE(1.0)`.
    Typed {
        /// Type name.
        name: String,
        /// Type parameters.
        parameters: Vec<StepValue>,
    },
}

/// One simple entity component.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StepEntity {
    /// Entity/type name. Canonical output uppercases it.
    pub name: String,
    /// Ordered parameter values.
    pub parameters: Vec<StepValue>,
}

/// One DATA-section instance, simple or complex.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StepInstance {
    /// Source-local Part-21 instance identifier (the number in `#N`).
    pub id: u64,
    /// One component for a simple entity; multiple for a complex entity.
    pub components: Vec<StepEntity>,
}

/// Mandatory Part-21 header records.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StepHeader {
    /// `FILE_DESCRIPTION` parameters.
    pub file_description: Vec<StepValue>,
    /// `FILE_NAME` parameters.
    pub file_name: Vec<StepValue>,
    /// `FILE_SCHEMA` parameters.
    pub file_schema: Vec<StepValue>,
}

/// A syntax-only Part-21 document.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StepDocument {
    /// Mandatory header records.
    pub header: StepHeader,
    /// DATA instances. Canonical writing orders them by `id`.
    pub instances: Vec<StepInstance>,
}

impl StepDocument {
    /// Return the declared FILE_SCHEMA strings after structural validation.
    ///
    /// # Errors
    /// [`IoError`] if the header is not in the supported bounded shape.
    pub fn schema_identifiers(&self) -> Result<Vec<String>, IoError> {
        schema_identifiers(&self.header, &accessor_limits())
    }

    /// Return the non-authoritative AP203/AP214 label hint.
    ///
    /// # Errors
    /// [`IoError`] if `FILE_SCHEMA` is malformed.
    pub fn profile_hint(&self) -> Result<StepProfileHint, IoError> {
        Ok(profile_hint(&self.schema_identifiers()?))
    }

    /// Admit only a caller-selected FILE_SCHEMA identifier.
    ///
    /// Matching is ASCII case-insensitive but otherwise exact. This is a
    /// declaration gate, not proof that the instance graph conforms to the
    /// named EXPRESS schema.
    ///
    /// # Errors
    /// [`IoError::Unsupported`] if the expected schema was not declared, or a
    /// structural error if FILE_SCHEMA itself is malformed.
    pub fn require_declared_schema(&self, expected: &str) -> Result<(), IoError> {
        let limits = accessor_limits();
        validate_string(expected, &limits)?;
        if expected.is_empty() {
            return Err(malformed(0, "expected FILE_SCHEMA identifier is empty"));
        }
        let schemas = schema_identifiers(&self.header, &limits)?;
        if schemas
            .iter()
            .any(|schema| schema.eq_ignore_ascii_case(expected))
        {
            return Ok(());
        }
        Err(IoError::Unsupported {
            what: format!(
                "expected FILE_SCHEMA {expected:?}; document declares {:?}",
                schemas
            ),
        })
    }
}

fn accessor_limits() -> StepLimits {
    StepLimits {
        max_schema_ids: usize::MAX,
        ..StepLimits::default()
    }
}

/// Syntax receipt emitted with a parsed Part-21 document.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StepStructureReceipt {
    /// Non-cryptographic FNV-1a fingerprint of the exact source bytes.
    source_fingerprint: u64,
    /// Non-cryptographic FNV-1a fingerprint of canonical-layout bytes.
    canonical_layout_fingerprint: u64,
    /// FILE_SCHEMA strings, preserving their declared spelling.
    schema_identifiers: Vec<String>,
    /// Non-authoritative AP-family hint.
    profile_hint: StepProfileHint,
    /// Number of DATA instances.
    instance_count: usize,
    /// Number of entity references across all instance parameters.
    reference_count: usize,
    /// Largest source-local instance identifier, or zero for an empty DATA section.
    max_instance_id: u64,
    /// Exact limits under which the source and canonical layout were admitted.
    limits: StepLimits,
}

impl StepStructureReceipt {
    /// Exact-source FNV-1a fingerprint. This is not collision-resistant.
    #[must_use]
    pub const fn source_fingerprint(&self) -> u64 {
        self.source_fingerprint
    }

    /// Canonical-layout FNV-1a fingerprint. This is not collision-resistant.
    #[must_use]
    pub const fn canonical_layout_fingerprint(&self) -> u64 {
        self.canonical_layout_fingerprint
    }

    /// Declared FILE_SCHEMA identifiers.
    #[must_use]
    pub fn schema_identifiers(&self) -> &[String] {
        &self.schema_identifiers
    }

    /// Non-authoritative application-protocol label hint.
    #[must_use]
    pub const fn profile_hint(&self) -> StepProfileHint {
        self.profile_hint
    }

    /// Number of DATA instances.
    #[must_use]
    pub const fn instance_count(&self) -> usize {
        self.instance_count
    }

    /// Number of entity references.
    #[must_use]
    pub const fn reference_count(&self) -> usize {
        self.reference_count
    }

    /// Largest source-local instance identifier, or zero for empty DATA.
    #[must_use]
    pub const fn max_instance_id(&self) -> u64 {
        self.max_instance_id
    }

    /// Limits used for admission and canonical-layout serialization.
    #[must_use]
    pub const fn limits(&self) -> StepLimits {
        self.limits
    }

    /// Canonical JSON suitable for a syntax-stage ledger event.
    #[must_use]
    pub fn to_json(&self) -> String {
        let mut out = String::from(
            "{\"kind\":\"step-part21-structure-receipt\",\"authority\":\"syntax-only\",",
        );
        let _ = write!(
            out,
            "\"syntax_version\":\"{}\",\"crate_version\":\"{}\",\
             \"source_fingerprint_fnv1a64\":\"{:016x}\",\
             \"canonical_layout_fingerprint_fnv1a64\":\"{:016x}\",\
             \"profile_hint\":\"{}\",\"schemas\":[",
            STEP_SYNTAX_VERSION,
            crate::VERSION,
            self.source_fingerprint,
            self.canonical_layout_fingerprint,
            self.profile_hint.as_str()
        );
        for (index, schema) in self.schema_identifiers.iter().enumerate() {
            if index > 0 {
                out.push(',');
            }
            push_json_string(&mut out, schema);
        }
        let _ = write!(
            out,
            "],\"instances\":{},\"references\":{},\"max_instance_id\":{},\
             \"limits\":{{\"bytes\":{},\"tokens\":{},\"instances\":{},\"values\":{},\
             \"nesting\":{},\"string_bytes\":{},\"number_bytes\":{},\"identifier_bytes\":{},\
             \"schema_ids\":{},\"components_per_instance\":{}}},\
             \"no_claim\":\"no EXPRESS-schema, geometry, tessellation, B-rep, AP-conformance, or certificate authority\"}}",
            self.instance_count,
            self.reference_count,
            self.max_instance_id,
            self.limits.max_bytes,
            self.limits.max_tokens,
            self.limits.max_instances,
            self.limits.max_values,
            self.limits.max_nesting,
            self.limits.max_string_bytes,
            self.limits.max_number_bytes,
            self.limits.max_identifier_bytes,
            self.limits.max_schema_ids,
            self.limits.max_components_per_instance
        );
        out
    }
}

/// A parsed document paired with its source-bound syntax receipt.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ParsedStep {
    /// Parsed syntax tree.
    document: StepDocument,
    /// Versioned limits, non-cryptographic fingerprints, and graph counts.
    receipt: StepStructureReceipt,
}

impl ParsedStep {
    /// Borrow the admitted document. The paired receipt cannot become stale.
    #[must_use]
    pub const fn document(&self) -> &StepDocument {
        &self.document
    }

    /// Borrow the immutable source-bound syntax receipt.
    #[must_use]
    pub const fn receipt(&self) -> &StepStructureReceipt {
        &self.receipt
    }

    /// Consume the parse result and discard its receipt before mutation.
    #[must_use]
    pub fn into_document(self) -> StepDocument {
        self.document
    }
}

/// Parse the bounded syntax-only Part-21 subset with default limits.
///
/// # Errors
/// [`IoError`] on malformed, unsupported, oversized, duplicate-ID, or
/// dangling-reference input.
pub fn parse_step(bytes: &[u8]) -> Result<ParsedStep, IoError> {
    parse_step_with_limits(bytes, StepLimits::default())
}

/// Parse the bounded syntax-only Part-21 subset with explicit limits.
///
/// # Errors
/// [`IoError`] on malformed, unsupported, oversized, duplicate-ID, or
/// dangling-reference input.
pub fn parse_step_with_limits(bytes: &[u8], limits: StepLimits) -> Result<ParsedStep, IoError> {
    validate_limits(limits)?;
    if bytes.len() > limits.max_bytes {
        return Err(resource(format!(
            "Part-21 input is {} bytes; cap is {}",
            bytes.len(),
            limits.max_bytes
        )));
    }
    if !bytes.is_ascii() {
        let at = bytes.iter().position(|byte| !byte.is_ascii()).unwrap_or(0);
        return Err(IoError::Unsupported {
            what: format!(
                "non-ASCII Part-21 source byte at {at}; encoded character directives are staged"
            ),
        });
    }

    let mut parser = Parser::new(bytes, limits);
    let document = parser.parse_document()?;
    let stats = validate_document(&document, limits)?;
    let canonical = write_step_with_limits(&document, limits)?;
    let schemas = schema_identifiers(&document.header, &limits)?;
    let receipt = StepStructureReceipt {
        source_fingerprint: fs_obs::fnv1a64(bytes),
        canonical_layout_fingerprint: fs_obs::fnv1a64(&canonical),
        profile_hint: profile_hint(&schemas),
        schema_identifiers: schemas,
        instance_count: document.instances.len(),
        reference_count: stats.reference_count,
        max_instance_id: stats.max_instance_id,
        limits,
    };
    Ok(ParsedStep { document, receipt })
}

/// Canonically serialize a syntax-only Part-21 document with default limits.
///
/// DATA instances are ordered by numeric ID. Parameter order, complex-entity
/// component order, and number-token spelling are preserved because changing
/// them without an EXPRESS schema can change semantics.
///
/// # Errors
/// [`IoError`] if a caller-constructed document is malformed, exceeds limits,
/// contains duplicate identifiers, or contains dangling references.
pub fn write_step(document: &StepDocument) -> Result<Vec<u8>, IoError> {
    write_step_with_limits(document, StepLimits::default())
}

/// Canonically serialize with explicit resource limits.
///
/// # Errors
/// [`IoError`] under the same conditions as [`write_step`].
pub fn write_step_with_limits(
    document: &StepDocument,
    limits: StepLimits,
) -> Result<Vec<u8>, IoError> {
    validate_limits(limits)?;
    validate_document(document, limits)?;
    preflight_output_tokens(document, limits)?;

    let mut out = String::new();
    push_bounded(&mut out, "ISO-10303-21;\nHEADER;\n", limits.max_bytes)?;
    write_header_record(
        &mut out,
        "FILE_DESCRIPTION",
        &document.header.file_description,
        limits,
    )?;
    write_header_record(&mut out, "FILE_NAME", &document.header.file_name, limits)?;
    write_header_record(
        &mut out,
        "FILE_SCHEMA",
        &document.header.file_schema,
        limits,
    )?;
    push_bounded(&mut out, "ENDSEC;\nDATA;\n", limits.max_bytes)?;

    let mut instances: Vec<&StepInstance> = Vec::new();
    instances
        .try_reserve(document.instances.len())
        .map_err(|_| resource("allocation refused while ordering Part-21 instances"))?;
    instances.extend(document.instances.iter());
    instances.sort_unstable_by_key(|instance| instance.id);
    for instance in instances {
        write_instance(&mut out, instance, limits)?;
    }
    push_bounded(&mut out, "ENDSEC;\nEND-ISO-10303-21;\n", limits.max_bytes)?;
    Ok(out.into_bytes())
}

fn preflight_output_tokens(document: &StepDocument, limits: StepLimits) -> Result<(), IoError> {
    let mut count = 0usize;
    add_output_tokens(&mut count, 4, limits.max_tokens)?; // magic; HEADER;
    for parameters in [
        document.header.file_description.as_slice(),
        document.header.file_name.as_slice(),
        document.header.file_schema.as_slice(),
    ] {
        add_output_tokens(&mut count, 1, limits.max_tokens)?; // record name
        count_parameter_tokens(parameters, &mut count, limits.max_tokens)?;
        add_output_tokens(&mut count, 1, limits.max_tokens)?; // semicolon
    }
    add_output_tokens(&mut count, 4, limits.max_tokens)?; // ENDSEC; DATA;
    for instance in &document.instances {
        add_output_tokens(&mut count, 2, limits.max_tokens)?; // #id and '='
        if instance.components.len() == 1 {
            count_entity_tokens(&instance.components[0], &mut count, limits.max_tokens)?;
        } else {
            add_output_tokens(&mut count, 1, limits.max_tokens)?; // '('
            for component in &instance.components {
                count_entity_tokens(component, &mut count, limits.max_tokens)?;
            }
            add_output_tokens(&mut count, 1, limits.max_tokens)?; // ')'
        }
        add_output_tokens(&mut count, 1, limits.max_tokens)?; // semicolon
    }
    add_output_tokens(&mut count, 5, limits.max_tokens)?; // ENDSEC; end magic; EOF
    Ok(())
}

fn count_entity_tokens(entity: &StepEntity, count: &mut usize, cap: usize) -> Result<(), IoError> {
    add_output_tokens(count, 1, cap)?;
    count_parameter_tokens(&entity.parameters, count, cap)
}

fn count_parameter_tokens(
    values: &[StepValue],
    count: &mut usize,
    cap: usize,
) -> Result<(), IoError> {
    add_output_tokens(count, 2, cap)?; // parentheses
    add_output_tokens(count, values.len().saturating_sub(1), cap)?; // commas
    for value in values {
        match value {
            StepValue::Aggregate(values) => count_parameter_tokens(values, count, cap)?,
            StepValue::Typed { parameters, .. } => {
                add_output_tokens(count, 1, cap)?; // type name
                count_parameter_tokens(parameters, count, cap)?;
            }
            StepValue::Omitted
            | StepValue::Derived
            | StepValue::Number(_)
            | StepValue::String(_)
            | StepValue::Enumeration(_)
            | StepValue::Reference(_) => add_output_tokens(count, 1, cap)?,
        }
    }
    Ok(())
}

fn add_output_tokens(count: &mut usize, added: usize, cap: usize) -> Result<(), IoError> {
    *count = count
        .checked_add(added)
        .ok_or_else(|| resource("canonical Part-21 token count overflow"))?;
    if *count > cap {
        return Err(resource(format!(
            "canonical Part-21 token count exceeds cap {cap}"
        )));
    }
    Ok(())
}

#[derive(Debug)]
struct DocumentStats {
    reference_count: usize,
    max_instance_id: u64,
}

fn validate_limits(limits: StepLimits) -> Result<(), IoError> {
    if limits.max_bytes == 0
        || limits.max_tokens == 0
        || limits.max_instances == 0
        || limits.max_values == 0
        || limits.max_nesting == 0
        || limits.max_string_bytes == 0
        || limits.max_number_bytes == 0
        || limits.max_identifier_bytes == 0
        || limits.max_schema_ids == 0
        || limits.max_components_per_instance == 0
    {
        return Err(resource("all Part-21 limits must be nonzero"));
    }
    if limits.max_nesting > HARD_MAX_STEP_NESTING {
        return Err(resource(format!(
            "Part-21 nesting cap {} exceeds hard stack-safe ceiling {HARD_MAX_STEP_NESTING}",
            limits.max_nesting
        )));
    }
    if limits.max_string_bytes > HARD_MAX_STEP_STRING_BYTES {
        return Err(resource(format!(
            "Part-21 string-token cap {} exceeds syntax ceiling {HARD_MAX_STEP_STRING_BYTES}",
            limits.max_string_bytes
        )));
    }
    if limits.max_identifier_bytes > HARD_MAX_STEP_IDENTIFIER_BYTES {
        return Err(resource(format!(
            "Part-21 identifier cap {} exceeds syntax ceiling {HARD_MAX_STEP_IDENTIFIER_BYTES}",
            limits.max_identifier_bytes
        )));
    }
    Ok(())
}

fn validate_document(
    document: &StepDocument,
    limits: StepLimits,
) -> Result<DocumentStats, IoError> {
    if document.instances.len() > limits.max_instances {
        return Err(resource(format!(
            "Part-21 instance count {} exceeds cap {}",
            document.instances.len(),
            limits.max_instances
        )));
    }
    validate_header_shape(&document.header, &limits)?;

    let mut value_count = 0usize;
    validate_header_values(
        "FILE_DESCRIPTION",
        &document.header.file_description,
        &limits,
        &mut value_count,
    )?;
    validate_header_values(
        "FILE_NAME",
        &document.header.file_name,
        &limits,
        &mut value_count,
    )?;
    validate_header_values(
        "FILE_SCHEMA",
        &document.header.file_schema,
        &limits,
        &mut value_count,
    )?;

    let mut ids = Vec::new();
    ids.try_reserve(document.instances.len())
        .map_err(|_| resource("allocation refused for Part-21 instance identifiers"))?;
    let mut max_instance_id = 0u64;
    for instance in &document.instances {
        if instance.id == 0 {
            return Err(malformed(0, "Part-21 instance identifier #0 is invalid"));
        }
        ids.push(instance.id);
        max_instance_id = max_instance_id.max(instance.id);
    }
    ids.sort_unstable();
    if let Some(duplicate) = ids.windows(2).find(|pair| pair[0] == pair[1]) {
        return Err(malformed(
            0,
            format!("duplicate Part-21 instance identifier #{}", duplicate[0]),
        ));
    }

    for instance in &document.instances {
        if instance.components.is_empty() {
            return Err(malformed(
                0,
                format!("Part-21 instance #{} has no entity component", instance.id),
            ));
        }
        if instance.components.len() > limits.max_components_per_instance {
            return Err(resource(format!(
                "Part-21 instance #{} has {} components; cap is {}",
                instance.id,
                instance.components.len(),
                limits.max_components_per_instance
            )));
        }
        let mut component_names = BTreeSet::new();
        for component in &instance.components {
            validate_identifier(&component.name, limits, "entity name")?;
            let canonical_name = component.name.to_ascii_uppercase();
            if !component_names.insert(canonical_name.clone()) {
                return Err(malformed(
                    0,
                    format!(
                        "Part-21 instance #{} repeats complex component {canonical_name}",
                        instance.id
                    ),
                ));
            }
            validate_values(&component.parameters, 1, &limits, &mut value_count)?;
        }
    }

    let mut reference_count = 0usize;
    for instance in &document.instances {
        for component in &instance.components {
            validate_resolved_references(&component.parameters, &ids, &mut reference_count)?;
        }
    }
    Ok(DocumentStats {
        reference_count,
        max_instance_id,
    })
}

fn validate_header_shape(header: &StepHeader, limits: &StepLimits) -> Result<(), IoError> {
    let [
        StepValue::Aggregate(descriptions),
        StepValue::String(implementation_level),
    ] = header.file_description.as_slice()
    else {
        return Err(malformed(
            0,
            "FILE_DESCRIPTION must contain a nonempty string aggregate and implementation-level string",
        ));
    };
    validate_nonempty_string_aggregate(descriptions, "FILE_DESCRIPTION descriptions", limits)?;
    let valid_implementation_level =
        implementation_level
            .split_once(';')
            .is_some_and(|(edition, conformance)| {
                !edition.is_empty()
                    && edition.bytes().all(|byte| byte.is_ascii_digit())
                    && !conformance.is_empty()
                    && conformance.bytes().all(|byte| byte.is_ascii_digit())
            });
    if !valid_implementation_level {
        return Err(malformed(
            0,
            "FILE_DESCRIPTION implementation level must be decimal-edition;decimal-conformance",
        ));
    }

    let [
        StepValue::String(_),
        StepValue::String(_),
        StepValue::Aggregate(authors),
        StepValue::Aggregate(organizations),
        StepValue::String(_),
        StepValue::String(_),
        StepValue::String(_),
    ] = header.file_name.as_slice()
    else {
        return Err(malformed(
            0,
            "FILE_NAME must contain name, timestamp, author aggregate, organization aggregate, preprocessor, originating system, and authorization strings",
        ));
    };
    validate_nonempty_string_aggregate(authors, "FILE_NAME authors", limits)?;
    validate_nonempty_string_aggregate(organizations, "FILE_NAME organizations", limits)?;
    let _ = schema_identifiers(header, limits)?;
    Ok(())
}

fn validate_nonempty_string_aggregate(
    values: &[StepValue],
    field: &str,
    limits: &StepLimits,
) -> Result<(), IoError> {
    if values.is_empty() {
        return Err(malformed(0, format!("{field} must not be empty")));
    }
    for value in values {
        let StepValue::String(value) = value else {
            return Err(malformed(0, format!("{field} members must be strings")));
        };
        validate_string(value, limits)?;
    }
    Ok(())
}

fn validate_header_values(
    record: &str,
    values: &[StepValue],
    limits: &StepLimits,
    value_count: &mut usize,
) -> Result<(), IoError> {
    validate_values(values, 1, limits, value_count)?;
    if let Some(reference) = first_reference(values) {
        return Err(malformed(
            0,
            format!("{record} may not reference DATA instance #{reference}"),
        ));
    }
    Ok(())
}

fn validate_values(
    values: &[StepValue],
    depth: usize,
    limits: &StepLimits,
    count: &mut usize,
) -> Result<(), IoError> {
    if depth > limits.max_nesting {
        return Err(resource(format!(
            "Part-21 parameter nesting exceeds cap {}",
            limits.max_nesting
        )));
    }
    for value in values {
        *count = count
            .checked_add(1)
            .ok_or_else(|| resource("Part-21 parameter count overflow"))?;
        if *count > limits.max_values {
            return Err(resource(format!(
                "Part-21 parameter count exceeds cap {}",
                limits.max_values
            )));
        }
        match value {
            StepValue::Omitted | StepValue::Derived | StepValue::Reference(_) => {}
            StepValue::Number(number) => validate_number(number, 0, limits)?,
            StepValue::String(value) => validate_string(value, limits)?,
            StepValue::Enumeration(name) => {
                validate_identifier(name, *limits, "enumeration")?;
            }
            StepValue::Aggregate(nested) => {
                validate_values(nested, depth + 1, limits, count)?;
            }
            StepValue::Typed { name, parameters } => {
                validate_identifier(name, *limits, "typed-parameter name")?;
                if parameters.len() != 1 {
                    return Err(malformed(
                        0,
                        format!("Part-21 typed parameter {name} must contain exactly one value"),
                    ));
                }
                validate_values(parameters, depth + 1, limits, count)?;
            }
        }
    }
    Ok(())
}

fn validate_string(value: &str, limits: &StepLimits) -> Result<(), IoError> {
    let encoded_len = value
        .bytes()
        .try_fold(2usize, |length, byte| {
            length.checked_add(if byte == b'\'' { 2 } else { 1 })
        })
        .ok_or_else(|| resource("Part-21 encoded string length overflow"))?;
    if encoded_len > limits.max_string_bytes {
        return Err(resource(format!(
            "Part-21 encoded string token is {encoded_len} bytes; cap is {}",
            limits.max_string_bytes
        )));
    }
    if !value.is_ascii() {
        return Err(IoError::Unsupported {
            what: "non-ASCII constructed STEP string; encoded character directives are staged"
                .to_string(),
        });
    }
    if value.contains('\\') {
        return Err(IoError::Unsupported {
            what: "Part-21 encoded-character directives are staged".to_string(),
        });
    }
    if let Some((index, _)) = value
        .bytes()
        .enumerate()
        .find(|(_, byte)| *byte < b' ' && !matches!(*byte, b'\t' | b'\n' | b'\r'))
    {
        return Err(malformed(
            index,
            "control byte in constructed Part-21 string",
        ));
    }
    Ok(())
}

fn schema_identifiers(header: &StepHeader, limits: &StepLimits) -> Result<Vec<String>, IoError> {
    let [StepValue::Aggregate(schemas)] = header.file_schema.as_slice() else {
        return Err(malformed(
            0,
            "FILE_SCHEMA must contain exactly one aggregate of schema strings",
        ));
    };
    if schemas.is_empty() {
        return Err(malformed(0, "FILE_SCHEMA schema aggregate is empty"));
    }
    if schemas.len() > limits.max_schema_ids {
        return Err(resource(format!(
            "FILE_SCHEMA declares {} schemas; cap is {}",
            schemas.len(),
            limits.max_schema_ids
        )));
    }
    let mut result = Vec::new();
    result
        .try_reserve(schemas.len())
        .map_err(|_| resource("allocation refused for FILE_SCHEMA identifiers"))?;
    let mut unique = BTreeSet::new();
    for schema in schemas {
        let StepValue::String(schema) = schema else {
            return Err(malformed(
                0,
                "FILE_SCHEMA aggregate members must be strings",
            ));
        };
        validate_string(schema, limits)?;
        if schema.is_empty() {
            return Err(malformed(0, "FILE_SCHEMA contains an empty identifier"));
        }
        let folded = schema.to_ascii_uppercase();
        if !unique.insert(folded) {
            return Err(malformed(0, "FILE_SCHEMA contains a duplicate identifier"));
        }
        result.push(schema.clone());
    }
    Ok(result)
}

fn profile_hint(schemas: &[String]) -> StepProfileHint {
    let mut ap203 = false;
    let mut ap214 = false;
    for schema in schemas {
        let schema = schema.to_ascii_uppercase();
        ap203 |=
            schema == "AP203" || schema.starts_with("AP203_") || schema == "CONFIG_CONTROL_DESIGN";
        ap214 |= schema == "AP214"
            || schema.starts_with("AP214_")
            || schema == "AUTOMOTIVE_DESIGN"
            || schema.starts_with("AUTOMOTIVE_DESIGN_");
    }
    match (ap203, ap214) {
        (true, false) => StepProfileHint::Ap203,
        (false, true) => StepProfileHint::Ap214,
        (true, true) => StepProfileHint::Ambiguous,
        (false, false) => StepProfileHint::Other,
    }
}

fn first_reference(values: &[StepValue]) -> Option<u64> {
    for value in values {
        match value {
            StepValue::Reference(id) => return Some(*id),
            StepValue::Aggregate(values)
            | StepValue::Typed {
                parameters: values, ..
            } => {
                if let Some(reference) = first_reference(values) {
                    return Some(reference);
                }
            }
            StepValue::Omitted
            | StepValue::Derived
            | StepValue::Number(_)
            | StepValue::String(_)
            | StepValue::Enumeration(_) => {}
        }
    }
    None
}

fn validate_resolved_references(
    values: &[StepValue],
    ids: &[u64],
    count: &mut usize,
) -> Result<(), IoError> {
    for value in values {
        match value {
            StepValue::Reference(id) => {
                *count = count
                    .checked_add(1)
                    .ok_or_else(|| resource("Part-21 reference count overflow"))?;
                if ids.binary_search(id).is_err() {
                    return Err(malformed(0, format!("dangling Part-21 reference #{id}")));
                }
            }
            StepValue::Aggregate(values)
            | StepValue::Typed {
                parameters: values, ..
            } => validate_resolved_references(values, ids, count)?,
            StepValue::Omitted
            | StepValue::Derived
            | StepValue::Number(_)
            | StepValue::String(_)
            | StepValue::Enumeration(_) => {}
        }
    }
    Ok(())
}

fn write_header_record(
    out: &mut String,
    name: &str,
    parameters: &[StepValue],
    limits: StepLimits,
) -> Result<(), IoError> {
    push_bounded(out, name, limits.max_bytes)?;
    write_parameters(out, parameters, limits)?;
    push_bounded(out, ";\n", limits.max_bytes)
}

fn write_instance(
    out: &mut String,
    instance: &StepInstance,
    limits: StepLimits,
) -> Result<(), IoError> {
    push_bounded(out, "#", limits.max_bytes)?;
    push_bounded(out, &instance.id.to_string(), limits.max_bytes)?;
    push_bounded(out, "=", limits.max_bytes)?;
    if instance.components.len() == 1 {
        write_entity(out, &instance.components[0], limits)?;
    } else {
        push_bounded(out, "(", limits.max_bytes)?;
        for component in &instance.components {
            write_entity(out, component, limits)?;
        }
        push_bounded(out, ")", limits.max_bytes)?;
    }
    push_bounded(out, ";\n", limits.max_bytes)
}

fn write_entity(out: &mut String, entity: &StepEntity, limits: StepLimits) -> Result<(), IoError> {
    push_bounded(out, &entity.name.to_ascii_uppercase(), limits.max_bytes)?;
    write_parameters(out, &entity.parameters, limits)
}

fn write_parameters(
    out: &mut String,
    values: &[StepValue],
    limits: StepLimits,
) -> Result<(), IoError> {
    push_bounded(out, "(", limits.max_bytes)?;
    for (index, value) in values.iter().enumerate() {
        if index > 0 {
            push_bounded(out, ",", limits.max_bytes)?;
        }
        write_value(out, value, limits)?;
    }
    push_bounded(out, ")", limits.max_bytes)
}

fn write_value(out: &mut String, value: &StepValue, limits: StepLimits) -> Result<(), IoError> {
    match value {
        StepValue::Omitted => push_bounded(out, "$", limits.max_bytes),
        StepValue::Derived => push_bounded(out, "*", limits.max_bytes),
        StepValue::Number(number) => push_bounded(out, number, limits.max_bytes),
        StepValue::String(value) => {
            push_bounded(out, "'", limits.max_bytes)?;
            for part in value.split_inclusive('\'') {
                push_bounded(out, part, limits.max_bytes)?;
                if part.ends_with('\'') {
                    push_bounded(out, "'", limits.max_bytes)?;
                }
            }
            push_bounded(out, "'", limits.max_bytes)
        }
        StepValue::Enumeration(name) => {
            push_bounded(out, ".", limits.max_bytes)?;
            push_bounded(out, &name.to_ascii_uppercase(), limits.max_bytes)?;
            push_bounded(out, ".", limits.max_bytes)
        }
        StepValue::Reference(id) => {
            push_bounded(out, "#", limits.max_bytes)?;
            push_bounded(out, &id.to_string(), limits.max_bytes)
        }
        StepValue::Aggregate(values) => write_parameters(out, values, limits),
        StepValue::Typed { name, parameters } => {
            push_bounded(out, &name.to_ascii_uppercase(), limits.max_bytes)?;
            write_parameters(out, parameters, limits)
        }
    }
}

fn push_bounded(out: &mut String, value: &str, cap: usize) -> Result<(), IoError> {
    let next = out
        .len()
        .checked_add(value.len())
        .ok_or_else(|| resource("Part-21 output length overflow"))?;
    if next > cap {
        return Err(resource(format!(
            "canonical Part-21 output exceeds byte cap {cap}"
        )));
    }
    out.try_reserve(value.len())
        .map_err(|_| resource("allocation refused for canonical Part-21 output"))?;
    out.push_str(value);
    Ok(())
}

fn push_json_string(out: &mut String, value: &str) {
    out.push('"');
    for character in value.chars() {
        match character {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            character if character.is_control() => {
                let _ = write!(out, "\\u{:04x}", u32::from(character));
            }
            character => out.push(character),
        }
    }
    out.push('"');
}

fn validate_identifier(identifier: &str, limits: StepLimits, kind: &str) -> Result<(), IoError> {
    if identifier.is_empty() || identifier.len() > limits.max_identifier_bytes {
        return Err(malformed(
            0,
            format!(
                "Part-21 {kind} length {} is outside 1..={}",
                identifier.len(),
                limits.max_identifier_bytes
            ),
        ));
    }
    let mut bytes = identifier.bytes();
    let Some(first) = bytes.next() else {
        return Err(malformed(0, format!("empty Part-21 {kind}")));
    };
    if !first.is_ascii_alphabetic()
        || !bytes.all(|byte| byte.is_ascii_alphanumeric() || byte == b'_')
    {
        return Err(malformed(
            0,
            format!("invalid Part-21 {kind} {identifier:?}"),
        ));
    }
    Ok(())
}

fn validate_number(number: &str, at: usize, limits: &StepLimits) -> Result<(), IoError> {
    if number.len() > limits.max_number_bytes {
        return Err(resource(format!(
            "Part-21 number token is {} bytes; cap is {}",
            number.len(),
            limits.max_number_bytes
        )));
    }
    let bytes = number.as_bytes();
    if bytes.is_empty() {
        return Err(malformed(at, "empty Part-21 number"));
    }
    let mut index = 0usize;
    if matches!(bytes[index], b'+' | b'-') {
        index += 1;
    }
    let integer_start = index;
    while index < bytes.len() && bytes[index].is_ascii_digit() {
        index += 1;
    }
    let integer_digits = index - integer_start;
    let mut has_decimal_point = false;
    if index < bytes.len() && bytes[index] == b'.' {
        has_decimal_point = true;
        index += 1;
        while index < bytes.len() && bytes[index].is_ascii_digit() {
            index += 1;
        }
    }
    if integer_digits == 0 {
        return Err(malformed(at, format!("invalid Part-21 number {number:?}")));
    }
    if index < bytes.len() && bytes[index] == b'E' {
        if !has_decimal_point {
            return Err(malformed(
                at,
                format!("Part-21 real exponent requires a decimal point in {number:?}"),
            ));
        }
        index += 1;
        if index < bytes.len() && matches!(bytes[index], b'+' | b'-') {
            index += 1;
        }
        let exponent_start = index;
        while index < bytes.len() && bytes[index].is_ascii_digit() {
            index += 1;
        }
        if index == exponent_start {
            return Err(malformed(
                at,
                format!("Part-21 exponent has no digits in {number:?}"),
            ));
        }
    }
    if index != bytes.len() {
        return Err(malformed(at, format!("invalid Part-21 number {number:?}")));
    }
    Ok(())
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum TokenKind {
    StartMagic,
    EndMagic,
    Identifier(String),
    Reference(u64),
    Number(String),
    String(String),
    Enumeration(String),
    LeftParen,
    RightParen,
    Comma,
    Equal,
    Semicolon,
    Omitted,
    Derived,
    Eof,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct Token {
    kind: TokenKind,
    at: usize,
}

struct Lexer<'a> {
    bytes: &'a [u8],
    index: usize,
    tokens: usize,
    limits: StepLimits,
}

impl<'a> Lexer<'a> {
    fn new(bytes: &'a [u8], limits: StepLimits) -> Self {
        Self {
            bytes,
            index: 0,
            tokens: 0,
            limits,
        }
    }

    fn next(&mut self) -> Result<Token, IoError> {
        self.skip_trivia()?;
        let at = self.index;
        if at == self.bytes.len() {
            return self.emit(TokenKind::Eof, at);
        }
        if self.consume_magic(b"ISO-10303-21") {
            return self.emit(TokenKind::StartMagic, at);
        }
        if self.consume_magic(b"END-ISO-10303-21") {
            return self.emit(TokenKind::EndMagic, at);
        }
        let byte = self.bytes[self.index];
        match byte {
            b'(' => self.single(TokenKind::LeftParen, at),
            b')' => self.single(TokenKind::RightParen, at),
            b',' => self.single(TokenKind::Comma, at),
            b'=' => self.single(TokenKind::Equal, at),
            b';' => self.single(TokenKind::Semicolon, at),
            b'$' => self.single(TokenKind::Omitted, at),
            b'*' => self.single(TokenKind::Derived, at),
            b'#' => self.reference(at),
            b'\'' => self.string(at),
            b'.' => self.enumeration(at),
            b'+' | b'-' | b'0'..=b'9' => self.number(at),
            byte if byte.is_ascii_alphabetic() => self.identifier(at),
            b'"' => Err(IoError::Unsupported {
                what: format!("Part-21 binary literal at byte {at} is outside this syntax subset"),
            }),
            _ => Err(malformed(
                at,
                format!("unexpected Part-21 byte 0x{byte:02x}"),
            )),
        }
    }

    fn emit(&mut self, kind: TokenKind, at: usize) -> Result<Token, IoError> {
        self.tokens = self
            .tokens
            .checked_add(1)
            .ok_or_else(|| resource("Part-21 token count overflow"))?;
        if self.tokens > self.limits.max_tokens {
            return Err(resource(format!(
                "Part-21 token count exceeds cap {}",
                self.limits.max_tokens
            )));
        }
        Ok(Token { kind, at })
    }

    fn single(&mut self, kind: TokenKind, at: usize) -> Result<Token, IoError> {
        self.index += 1;
        self.emit(kind, at)
    }

    fn consume_magic(&mut self, magic: &[u8]) -> bool {
        let Some(candidate) = self.bytes.get(self.index..self.index + magic.len()) else {
            return false;
        };
        if candidate != magic {
            return false;
        }
        let boundary = self.bytes.get(self.index + magic.len()).copied();
        if boundary.is_some_and(|byte| byte.is_ascii_alphanumeric() || byte == b'_') {
            return false;
        }
        self.index += magic.len();
        true
    }

    fn skip_trivia(&mut self) -> Result<(), IoError> {
        loop {
            while self
                .bytes
                .get(self.index)
                .is_some_and(u8::is_ascii_whitespace)
            {
                self.index += 1;
            }
            if self.bytes.get(self.index..self.index + 2) != Some(b"/*") {
                return Ok(());
            }
            let comment_at = self.index;
            self.index += 2;
            let Some(end) = self.bytes[self.index..]
                .windows(2)
                .position(|window| window == b"*/")
            else {
                return Err(malformed(comment_at, "unterminated Part-21 comment"));
            };
            self.index += end + 2;
        }
    }

    fn reference(&mut self, at: usize) -> Result<Token, IoError> {
        self.index += 1;
        let start = self.index;
        while self.bytes.get(self.index).is_some_and(u8::is_ascii_digit) {
            self.index += 1;
        }
        if self.index == start {
            return Err(malformed(at, "Part-21 '#' has no instance digits"));
        }
        let id = parse_u64_decimal(&self.bytes[start..self.index], at)?;
        if id == 0 {
            return Err(malformed(at, "Part-21 instance identifier #0 is invalid"));
        }
        self.emit(TokenKind::Reference(id), at)
    }

    fn string(&mut self, at: usize) -> Result<Token, IoError> {
        self.index += 1;
        let mut value = String::new();
        loop {
            let Some(&byte) = self.bytes.get(self.index) else {
                return Err(malformed(at, "unterminated Part-21 string"));
            };
            self.index += 1;
            if byte == b'\'' {
                if self.bytes.get(self.index) == Some(&b'\'') {
                    self.index += 1;
                    push_char_bounded(
                        &mut value,
                        '\'',
                        self.limits.max_string_bytes,
                        "Part-21 string",
                    )?;
                    continue;
                }
                break;
            }
            if byte == b'\\' {
                return Err(IoError::Unsupported {
                    what: format!(
                        "Part-21 encoded-character directive at byte {} is staged",
                        self.index - 1
                    ),
                });
            }
            if byte < b' ' && !matches!(byte, b'\t' | b'\n' | b'\r') {
                return Err(malformed(self.index - 1, "control byte in Part-21 string"));
            }
            push_char_bounded(
                &mut value,
                char::from(byte),
                self.limits.max_string_bytes,
                "Part-21 string",
            )?;
        }
        let encoded_len = self.index - at;
        if encoded_len > self.limits.max_string_bytes {
            return Err(resource(format!(
                "Part-21 encoded string token is {encoded_len} bytes; cap is {}",
                self.limits.max_string_bytes
            )));
        }
        self.emit(TokenKind::String(value), at)
    }

    fn enumeration(&mut self, at: usize) -> Result<Token, IoError> {
        self.index += 1;
        let start = self.index;
        while self
            .bytes
            .get(self.index)
            .is_some_and(|byte| byte.is_ascii_alphanumeric() || *byte == b'_')
        {
            self.index += 1;
        }
        if self.index == start || self.bytes.get(self.index) != Some(&b'.') {
            return Err(malformed(at, "malformed Part-21 enumeration"));
        }
        if self.index - start > self.limits.max_identifier_bytes {
            return Err(resource(format!(
                "Part-21 enumeration exceeds byte cap {}",
                self.limits.max_identifier_bytes
            )));
        }
        if self.bytes[start..self.index]
            .iter()
            .any(u8::is_ascii_lowercase)
        {
            return Err(malformed(
                at,
                "Part-21 enumeration keywords must be uppercase",
            ));
        }
        let name = ascii_string(&self.bytes[start..self.index], "enumeration")?;
        validate_identifier(&name, self.limits, "enumeration")?;
        self.index += 1;
        self.emit(TokenKind::Enumeration(name), at)
    }

    fn identifier(&mut self, at: usize) -> Result<Token, IoError> {
        let start = self.index;
        self.index += 1;
        while self
            .bytes
            .get(self.index)
            .is_some_and(|byte| byte.is_ascii_alphanumeric() || *byte == b'_')
        {
            self.index += 1;
        }
        if self.index - start > self.limits.max_identifier_bytes {
            return Err(resource(format!(
                "Part-21 identifier exceeds byte cap {}",
                self.limits.max_identifier_bytes
            )));
        }
        if self.bytes[start..self.index]
            .iter()
            .any(u8::is_ascii_lowercase)
        {
            return Err(malformed(at, "Part-21 keywords must be uppercase"));
        }
        let value = ascii_string(&self.bytes[start..self.index], "identifier")?;
        self.emit(TokenKind::Identifier(value), at)
    }

    fn number(&mut self, at: usize) -> Result<Token, IoError> {
        let start = self.index;
        if self
            .bytes
            .get(self.index)
            .is_some_and(|byte| matches!(*byte, b'+' | b'-'))
        {
            self.index += 1;
        }
        while self.bytes.get(self.index).is_some_and(u8::is_ascii_digit) {
            self.index += 1;
        }
        if self.bytes.get(self.index) == Some(&b'.') {
            self.index += 1;
            while self.bytes.get(self.index).is_some_and(u8::is_ascii_digit) {
                self.index += 1;
            }
        }
        if self.bytes.get(self.index).is_some_and(|byte| *byte == b'E') {
            self.index += 1;
            if self
                .bytes
                .get(self.index)
                .is_some_and(|byte| matches!(*byte, b'+' | b'-'))
            {
                self.index += 1;
            }
            while self.bytes.get(self.index).is_some_and(u8::is_ascii_digit) {
                self.index += 1;
            }
        }
        let token_len = self.index - start;
        if token_len > self.limits.max_number_bytes {
            return Err(resource(format!(
                "Part-21 number token is {token_len} bytes; cap is {}",
                self.limits.max_number_bytes
            )));
        }
        let value = ascii_string(&self.bytes[start..self.index], "number")?;
        validate_number(&value, at, &self.limits)?;
        self.emit(TokenKind::Number(value), at)
    }
}

fn push_char_bounded(
    out: &mut String,
    character: char,
    cap: usize,
    what: &str,
) -> Result<(), IoError> {
    let next = out
        .len()
        .checked_add(character.len_utf8())
        .ok_or_else(|| resource(format!("{what} length overflow")))?;
    if next > cap {
        return Err(resource(format!("{what} exceeds byte cap {cap}")));
    }
    out.try_reserve(character.len_utf8())
        .map_err(|_| resource(format!("allocation refused for {what}")))?;
    out.push(character);
    Ok(())
}

fn ascii_string(bytes: &[u8], what: &str) -> Result<String, IoError> {
    let value = core::str::from_utf8(bytes)
        .map_err(|_| malformed(0, format!("non-ASCII Part-21 {what}")))?;
    let mut owned = String::new();
    owned
        .try_reserve(value.len())
        .map_err(|_| resource(format!("allocation refused for Part-21 {what}")))?;
    owned.push_str(value);
    Ok(owned)
}

fn parse_u64_decimal(bytes: &[u8], at: usize) -> Result<u64, IoError> {
    let mut value = 0u64;
    for byte in bytes {
        value = value
            .checked_mul(10)
            .and_then(|value| value.checked_add(u64::from(*byte - b'0')))
            .ok_or_else(|| malformed(at, "Part-21 instance identifier overflows u64"))?;
    }
    Ok(value)
}

struct Parser<'a> {
    lexer: Lexer<'a>,
    lookahead: Option<Token>,
    values: usize,
}

impl<'a> Parser<'a> {
    fn new(bytes: &'a [u8], limits: StepLimits) -> Self {
        Self {
            lexer: Lexer::new(bytes, limits),
            lookahead: None,
            values: 0,
        }
    }

    fn parse_document(&mut self) -> Result<StepDocument, IoError> {
        self.expect_kind(TokenKind::StartMagic, "ISO-10303-21")?;
        self.expect_kind(TokenKind::Semicolon, "';' after ISO-10303-21")?;
        self.expect_identifier("HEADER")?;
        self.expect_kind(TokenKind::Semicolon, "';' after HEADER")?;

        let file_description = self.parse_required_header("FILE_DESCRIPTION")?;
        let file_name = self.parse_required_header("FILE_NAME")?;
        let file_schema = self.parse_required_header("FILE_SCHEMA")?;
        self.expect_identifier("ENDSEC")?;
        self.expect_kind(TokenKind::Semicolon, "';' after HEADER ENDSEC")?;
        self.expect_identifier("DATA")?;
        self.expect_kind(TokenKind::Semicolon, "';' after DATA")?;

        let mut instances = Vec::new();
        loop {
            let token = self.peek()?.clone();
            match token.kind {
                TokenKind::Reference(_) => {
                    if instances.len() >= self.lexer.limits.max_instances {
                        return Err(resource(format!(
                            "Part-21 instance count exceeds cap {}",
                            self.lexer.limits.max_instances
                        )));
                    }
                    instances
                        .try_reserve(1)
                        .map_err(|_| resource("allocation refused for Part-21 instances"))?;
                    instances.push(self.parse_instance()?);
                }
                TokenKind::Identifier(name) if name == "ENDSEC" => break,
                _ => {
                    return Err(malformed(token.at, "expected DATA instance or ENDSEC"));
                }
            }
        }
        self.expect_identifier("ENDSEC")?;
        self.expect_kind(TokenKind::Semicolon, "';' after DATA ENDSEC")?;
        self.expect_kind(TokenKind::EndMagic, "END-ISO-10303-21")?;
        self.expect_kind(TokenKind::Semicolon, "';' after END-ISO-10303-21")?;
        self.expect_kind(TokenKind::Eof, "end of Part-21 input")?;

        Ok(StepDocument {
            header: StepHeader {
                file_description,
                file_name,
                file_schema,
            },
            instances,
        })
    }

    fn parse_required_header(&mut self, name: &str) -> Result<Vec<StepValue>, IoError> {
        self.expect_identifier(name)?;
        let parameters = self.parse_parameter_list(1)?;
        self.expect_kind(TokenKind::Semicolon, "';' after header record")?;
        Ok(parameters)
    }

    fn parse_instance(&mut self) -> Result<StepInstance, IoError> {
        let id = self.take_reference("instance identifier")?;
        self.expect_kind(TokenKind::Equal, "'=' after instance identifier")?;
        let components = if matches!(&self.peek()?.kind, TokenKind::LeftParen) {
            self.take()?;
            let mut components = Vec::new();
            while !matches!(&self.peek()?.kind, TokenKind::RightParen) {
                if components.len() >= self.lexer.limits.max_components_per_instance {
                    return Err(resource(format!(
                        "complex Part-21 component count exceeds cap {}",
                        self.lexer.limits.max_components_per_instance
                    )));
                }
                components
                    .try_reserve(1)
                    .map_err(|_| resource("allocation refused for complex entity components"))?;
                components.push(self.parse_entity()?);
            }
            self.take()?;
            if components.len() < 2 {
                return Err(malformed(
                    0,
                    "complex Part-21 instance must contain at least two components",
                ));
            }
            components
        } else {
            vec![self.parse_entity()?]
        };
        self.expect_kind(TokenKind::Semicolon, "';' after DATA instance")?;
        Ok(StepInstance { id, components })
    }

    fn parse_entity(&mut self) -> Result<StepEntity, IoError> {
        let name = self.take_identifier("entity name")?;
        let parameters = self.parse_parameter_list(1)?;
        Ok(StepEntity { name, parameters })
    }

    fn parse_parameter_list(&mut self, depth: usize) -> Result<Vec<StepValue>, IoError> {
        if depth > self.lexer.limits.max_nesting {
            return Err(resource(format!(
                "Part-21 parameter nesting exceeds cap {}",
                self.lexer.limits.max_nesting
            )));
        }
        self.expect_kind(TokenKind::LeftParen, "'('")?;
        let mut values = Vec::new();
        if matches!(&self.peek()?.kind, TokenKind::RightParen) {
            self.take()?;
            return Ok(values);
        }
        loop {
            values
                .try_reserve(1)
                .map_err(|_| resource("allocation refused for Part-21 parameters"))?;
            values.push(self.parse_value(depth)?);
            if matches!(&self.peek()?.kind, TokenKind::Comma) {
                self.take()?;
                if matches!(&self.peek()?.kind, TokenKind::RightParen) {
                    return Err(malformed(
                        self.peek()?.at,
                        "trailing comma in Part-21 parameter list",
                    ));
                }
                continue;
            }
            self.expect_kind(TokenKind::RightParen, "')'")?;
            return Ok(values);
        }
    }

    fn parse_value(&mut self, depth: usize) -> Result<StepValue, IoError> {
        self.values = self
            .values
            .checked_add(1)
            .ok_or_else(|| resource("Part-21 parameter count overflow"))?;
        if self.values > self.lexer.limits.max_values {
            return Err(resource(format!(
                "Part-21 parameter count exceeds cap {}",
                self.lexer.limits.max_values
            )));
        }
        let token = self.take()?;
        match token.kind {
            TokenKind::Omitted => Ok(StepValue::Omitted),
            TokenKind::Derived => Ok(StepValue::Derived),
            TokenKind::Reference(id) => Ok(StepValue::Reference(id)),
            TokenKind::Number(value) => Ok(StepValue::Number(value)),
            TokenKind::String(value) => Ok(StepValue::String(value)),
            TokenKind::Enumeration(value) => Ok(StepValue::Enumeration(value)),
            TokenKind::LeftParen => {
                self.lookahead = Some(Token {
                    kind: TokenKind::LeftParen,
                    at: token.at,
                });
                Ok(StepValue::Aggregate(self.parse_parameter_list(depth + 1)?))
            }
            TokenKind::Identifier(name) => {
                let parameters = self.parse_parameter_list(depth + 1)?;
                if parameters.len() != 1 {
                    return Err(malformed(
                        token.at,
                        format!("Part-21 typed parameter {name} must contain exactly one value"),
                    ));
                }
                Ok(StepValue::Typed { name, parameters })
            }
            _ => Err(malformed(token.at, "expected Part-21 parameter value")),
        }
    }

    fn peek(&mut self) -> Result<&Token, IoError> {
        if self.lookahead.is_none() {
            self.lookahead = Some(self.lexer.next()?);
        }
        self.lookahead
            .as_ref()
            .ok_or_else(|| malformed(self.lexer.index, "missing Part-21 lookahead"))
    }

    fn take(&mut self) -> Result<Token, IoError> {
        if let Some(token) = self.lookahead.take() {
            Ok(token)
        } else {
            self.lexer.next()
        }
    }

    fn expect_kind(&mut self, expected: TokenKind, description: &str) -> Result<(), IoError> {
        let token = self.take()?;
        if core::mem::discriminant(&token.kind) == core::mem::discriminant(&expected) {
            Ok(())
        } else {
            Err(malformed(token.at, format!("expected {description}")))
        }
    }

    fn expect_identifier(&mut self, expected: &str) -> Result<(), IoError> {
        let token = self.take()?;
        match token.kind {
            TokenKind::Identifier(value) if value == expected => Ok(()),
            _ => Err(malformed(token.at, format!("expected {expected}"))),
        }
    }

    fn take_identifier(&mut self, description: &str) -> Result<String, IoError> {
        let token = self.take()?;
        match token.kind {
            TokenKind::Identifier(value) => Ok(value),
            _ => Err(malformed(token.at, format!("expected {description}"))),
        }
    }

    fn take_reference(&mut self, description: &str) -> Result<u64, IoError> {
        let token = self.take()?;
        match token.kind {
            TokenKind::Reference(value) => Ok(value),
            _ => Err(malformed(token.at, format!("expected {description}"))),
        }
    }
}

fn malformed(at: usize, what: impl Into<String>) -> IoError {
    IoError::Malformed {
        at,
        what: what.into(),
    }
}

fn resource(what: impl Into<String>) -> IoError {
    IoError::ResourceBound { what: what.into() }
}
