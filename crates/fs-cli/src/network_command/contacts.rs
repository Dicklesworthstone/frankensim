//! Explicit fixed-resistance contacts between matching, separately owned P1
//! traces. Reuses fs-conduction's contact geometry, assembly and flux reporting.
//! Inline cards encode CALLER DECLARATIONS, not independently measured material
//! authority. No contact heat is counted as an external source or air exchange.
use super::*;
use fs_conduction::{InterfaceFacePair, InterfaceFlux, InterfaceResistance,
    InterfaceSurface, ThermalInterfaces, AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS as RD,
    AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY as RP};
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy,
    SurfaceSpec, SystemContext, UncertaintyModel};

mod sensitivity;
use sensitivity::Trace;

#[derive(Debug)]
struct Declaration {
    name: String,
    source: String,
    side_a_material: String,
    side_b_material: String,
    resistance: f64,
    pair_count: usize,
    traces: Vec<Trace>,
}

#[derive(Debug)]
pub(super) struct Contacts {
    pub interfaces: ThermalInterfaces,
    declarations: Vec<Declaration>,
    vertex_count: usize,
}

impl Contacts {
    /// Check every coincident pair and all external/contact face ownership
    /// before numerical work. An absent declaration never means perfect contact
    /// or an insulated gap between duplicated coincident faces.
    pub fn parse(value: Option<&J>, mesh: &ConductionMesh, surfaces: &[Surface],
        adiabatic: bool) -> Result<Option<Self>> {
        let candidates = ThermalInterfaces::coincident_face_pairs(mesh).map_err(producer)?;
        let Some(value) = value else {
            if !candidates.is_empty() { return Err(bad("coincident solid traces require explicit solid.contacts")); }
            if !adiabatic && surfaces.iter().map(|s| s.faces.len()).sum::<usize>() != mesh.boundary().len() {
                return Err(bad("every non-contact exterior face must be cooled unless adiabatic_remainder is true"));
            }
            return Ok(None);
        };
        let entries = array(value, "solid.contacts", 4096)?;
        if entries.is_empty() { return Err(bad("solid.contacts must be nonempty when supplied")); }
        let slots: BTreeMap<_, _> = mesh.boundary().iter().enumerate().map(|(i, f)| (f.vertices, i)).collect();
        let mut owned: BTreeSet<_> = surfaces.iter().flat_map(|s| s.faces.iter().copied()).collect();
        let mut names: BTreeSet<_> = surfaces.iter().map(|s| s.name.clone()).collect();
        let mut declarations = Vec::new();
        let mut bound = Vec::new();
        let mut pending_traces = Vec::new();
        for entry in entries {
            object(entry, &["name", "source", "side_a_material", "side_b_material",
                "resistance_m2_k_w", "face_pairs"], "contact")?;
            let mut row = Declaration {
                name: string(get(entry, "name")?, "contact.name")?,
                source: string(get(entry, "source")?, "contact.source")?,
                side_a_material: string(get(entry, "side_a_material")?, "contact.side_a_material")?,
                side_b_material: string(get(entry, "side_b_material")?, "contact.side_b_material")?,
                resistance: positive(get(entry, "resistance_m2_k_w")?, "contact.resistance_m2_k_w")?,
                pair_count: 0,
                traces: Vec::new(),
            };
            if !names.insert(row.name.clone()) { return Err(bad("contact and cooling surface names must be distinct")); }
            if !(1.0 / row.resistance).is_finite() { return Err(bad("contact conductance density is not representable")); }
            let mut pairs = Vec::new();
            for pair in array(get(entry, "face_pairs")?, "contact.face_pairs", mesh.boundary().len()/2)? {
                object(pair, &["side_a", "side_b"], "contact face pair")?;
                let mut side_a = indices::<3>(get(pair, "side_a")?, "contact.side_a", mesh.vertex_count())?;
                let mut side_b = indices::<3>(get(pair, "side_b")?, "contact.side_b", mesh.vertex_count())?;
                side_a.sort_unstable(); side_b.sort_unstable();
                let a = *slots.get(&side_a).ok_or_else(|| bad("contact side A is not an exterior trace triangle"))?;
                let b = *slots.get(&side_b).ok_or_else(|| bad("contact side B is not an exterior trace triangle"))?;
                if !owned.insert(side_a) || !owned.insert(side_b) {
                    return Err(bad("a contact face is repeated or also has an external cooling owner"));
                }
                pairs.push(InterfaceFacePair { side_a: a, side_b: b });
            }
            if pairs.is_empty() { return Err(bad("a contact requires at least one face pair")); }
            row.pair_count = pairs.len();
            pending_traces.push(pairs.clone());
            bound.push(InterfaceSurface::new(row.name.clone(), pairs, resistance(&row)?).map_err(producer)?);
            declarations.push(row);
        }
        if !adiabatic && owned.len() != mesh.boundary().len() {
            return Err(bad("every non-contact exterior face must be cooled unless adiabatic_remainder is true"));
        }
        // Only ownership is being bound here, not a physical solve. Coefficients
        // and references change later, but the selected trace faces do not.
        let mut boundary = ThermalBoundaryBuilder::new(mesh);
        for surface in surfaces {
            boundary = boundary.region(&surface.name, |face| surface.faces.contains(&face.vertices),
                ThermalBc::robin(1.0, 300.0).map_err(producer)?).map_err(producer)?;
        }
        let boundary = boundary.adiabatic_remainder().finish().map_err(producer)?;
        let interfaces = ThermalInterfaces::new(mesh, &boundary, bound).map_err(producer)?;
        // Preserve the complete admitted pairing; do not infer an effective
        // scalar surface from aggregate heat or a representative jump.
        for (row, pairs) in declarations.iter_mut().zip(pending_traces) {
            row.traces = Trace::bind(mesh, &pairs)?;
        }
        declarations.sort_by(|a, b| a.name.cmp(&b.name));
        Ok(Some(Self { interfaces, declarations, vertex_count: mesh.vertex_count() }))
    }

    pub fn render(&self, fluxes: &[InterfaceFlux]) -> Result<String> {
        if fluxes.len() != self.declarations.len() { return Err(bad("contact flux/declaration count mismatch")); }
        let mut rows = Vec::new();
        for row in &self.declarations {
            let flux = fluxes.iter().find(|flux| flux.interface == row.name)
                .ok_or_else(|| bad("contact is missing its evaluated flux"))?;
            rows.push(format!("{{\"name\":{},\"source\":{},\"side_a_material\":{},\"side_b_material\":{},\"resistance_m2_k_w\":{},\"face_pairs\":{},\"area_m2\":{},\"conductance_w_k\":{},\"mean_jump_a_minus_b_k\":{},\"heat_a_to_b_w\":{},\"authority\":\"caller-declared constant contact; inline card does not add measured material authority\",\"uncertainty\":null}}",
                quote(&row.name), quote(&row.source), quote(&row.side_a_material), quote(&row.side_b_material),
                num(row.resistance)?, row.pair_count, num(flux.area_m2)?, num(flux.conductance_w_per_k)?,
                num(flux.mean_jump_k)?, num(flux.heat_rate_a_to_b_w)?));
        }
        Ok(format!("[{}]", rows.join(",")))
    }
}

fn resistance(row: &Declaration) -> Result<InterfaceResistance> {
    // A real inline claim whose provenance explicitly says caller declaration.
    // No fabricated experiment, material-data license, validity measurement or
    // uncertainty band is used to get past the card-backed contact API.
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim { key: PropertyKey::new(RP, RD),
        value: PropertyValue::Scalar { value: row.resistance, dims: RD },
        validity: ValidityDomain::unconstrained(), uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity, observations: Vec::new(),
        provenance: Provenance { source: format!("caller-declared cooling-network contact {}: {}", row.name, row.source),
            license: "unspecified".into(), artifact: None } }).map_err(producer)?;
    let side = |label: &str| SurfaceSpec { material: MaterialStateId {
        chemistry: label.to_string(), phase: "caller-unspecified".into(),
        process: "caller-declared contact label".into(), revision: 0 },
        texture_frame: "caller-unspecified".into() };
    let card = InterfaceSystemCard::assemble(side(&row.side_a_material), side(&row.side_b_material),
        SystemContext { medium: "caller-unspecified".into(), third_body: None,
            environment: "caller-unspecified".into(), history: "caller-unspecified".into() },
        claims, Vec::new()).map_err(producer)?;
    InterfaceResistance::from_card(&row.name, &card, &QueryPoint::new(), SelectionPolicy::SingleClaimOnly).map_err(producer)
}

#[cfg(test)]
mod tests;
