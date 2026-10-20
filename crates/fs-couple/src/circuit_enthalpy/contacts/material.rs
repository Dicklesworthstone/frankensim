//! Same ordered interface card and thermal-link geometry for contact current/heat.
mod prepared;
pub use prepared::PreparedElectricalContactNetwork;
use super::*;
use crate::circuit_enthalpy::material::{MaterialNetwork,MaterialNetworkError,MaterialNetworkSample};
use fs_material::conductor::ElectricalPropertyBasis;
use fs_material::thermal::interface::ThermalInterfacePoint;
use fs_material::thermal::source::ThermalMaterialError;
use fs_material::state_point::{MaterialPropertySelection,MaterialStatePointError};
use fs_material::electrical_interface::resolve_electrical_interface_state;
use fs_matdb::PropertyValue;

/// One electrical contact using an existing thermal link's A/B order and area.
pub struct MaterialElectricalContact<'a> {
    /// Stable thermal-link ID, also the electrical contact ID. One binding per link.
    pub link: u64,
    /// Bulk circuit branch containing this series contact; not a resistor-only index.
    pub branch: usize,
    /// Original interface card, both bulk identities and explicit evolving axes.
    pub point: ThermalInterfacePoint<'a>,
    /// Original contact claim policy or exact property pin.
    pub selection: MaterialPropertySelection,
    /// Explicit heat fraction entering the link's A cell in [0,1].
    pub heat_fraction_to_first: f64,
}
/// Data-only assembly; the original material network and all physical owners remain.
pub struct ElectricalContactNetwork<'a> {
    /// Existing caloric, electrical and conductive material bindings.
    pub bulk: &'a MaterialNetwork<'a>,
    /// Complete ordered-interface/branch assignments.
    pub contacts: &'a [MaterialElectricalContact<'a>],
    /// Explicit resistivity/conductivity basis for every original bulk resistor.
    pub electrical_bases: &'a [(usize,ElectricalPropertyBasis)],
}
#[derive(Debug)]
pub enum ElectricalContactNetworkError {
    /// Incomplete binding, identity, unsupported contact law or partition.
    Invalid(&'static str),
    /// Existing bulk material resolver refused.
    Bulk(MaterialNetworkError),
    /// Existing shared two-sided thermal condition resolver refused.
    State(ThermalMaterialError),
    /// Existing electrical interface resolver refused.
    Electrical(MaterialStatePointError),
    /// Invalid contact geometry or circuit reduction.
    Circuit(MaterialCircuitError),
    /// Thermal geometry or cancellation refusal.
    Thermal(DepositionError),
}
impl core::fmt::Display for ElectricalContactNetworkError {
    fn fmt(&self,f:&mut core::fmt::Formatter<'_>)->core::fmt::Result {write!(f,"electrical contact network: {self:?}")}
}
impl std::error::Error for ElectricalContactNetworkError {}

impl ElectricalContactNetwork<'_> {
    fn routes(&self,thermal:&DepositionSimulation<'_>)->Result<Vec<ContactHeatCells>,ElectricalContactNetworkError> {
        if self.contacts.len()>MAX_SERIES_CONTACTS {return Err(ElectricalContactNetworkError::Invalid("contact capacity exceeded"));}
        let mut routes=BTreeMap::new();
        for binding in self.contacts {
            let link=thermal.links().get(&binding.link).ok_or(ElectricalContactNetworkError::Invalid("electrical contact thermal link is inactive"))?;
            if !binding.heat_fraction_to_first.is_finite() || !(0.0..=1.0).contains(&binding.heat_fraction_to_first)
                || routes.insert(binding.link,ContactHeatCells {id:binding.link,cells:[link.a,link.b],
                    heat_fraction_to_first:binding.heat_fraction_to_first}).is_some() {
                return Err(ElectricalContactNetworkError::Invalid("duplicate contact or invalid heat fraction"));
            }
        }
        Ok(routes.into_values().collect())
    }
    /// Complete explicit electrical heat destinations in stable contact-ID order.
    /// No partition is inferred; each original thermal link supplies its A/B cells.
    pub fn heat_destinations(&self,thermal:&DepositionSimulation<'_>)
        ->Result<Vec<ContactHeatCells>,ElectricalContactNetworkError> {
        self.routes(thermal)
    }
    /// Fresh same-trial bulk and interface resolution. The thermal link owns
    /// contact area and source A/B order; missing area never becomes a short.
    /// All phase/card identities are checked by the original state resolver.
    pub fn resolve(&self,thermal:&DepositionSimulation<'_>,mut cancelled:impl FnMut()->bool)
        ->Result<ContactThermalCoefficients<MaterialNetworkSample>,ElectricalContactNetworkError> {
        use ElectricalContactNetworkError as Error;
        if cancelled() {return Err(Error::Thermal(DepositionError::Cancelled));}
        self.routes(thermal)?;
        let bulk=self.bulk.resolve_with_bases(thermal,self.electrical_bases,&mut cancelled).map_err(Error::Bulk)?;
        self.resolve_contacts(thermal,bulk,None,cancelled)
    }

    // Both fresh-policy and prepared callers use this same interface evaluator.
    fn resolve_contacts(&self,thermal:&DepositionSimulation<'_>,bulk:CircuitThermalCoefficients<MaterialNetworkSample>,
        pins:Option<&BTreeMap<u64,MaterialPropertySelection>>,mut cancelled:impl FnMut()->bool)
        ->Result<ContactThermalCoefficients<MaterialNetworkSample>,ElectricalContactNetworkError> {
        use ElectricalContactNetworkError as Error;
        if cancelled() {return Err(Error::Thermal(DepositionError::Cancelled));}
        self.routes(thermal)?;
        let mut contacts=Vec::with_capacity(self.contacts.len());
        for binding in self.contacts {
            if cancelled() {return Err(Error::Thermal(DepositionError::Cancelled));}
            let link=&thermal.links()[&binding.link];
            if !self.bulk.resistors.iter().any(|r|r.destination.branch==binding.branch
                && (r.destination.cell==link.a || r.destination.cell==link.b)) {
                return Err(Error::Invalid("contact branch must belong to one of its two material cells"));
            }
            let phases=[link.a,link.b].map(|id|bulk.source.cells.iter().find(|c|c.cell==id).map(|c|c.phase));
            let [Some(a),Some(b)]=phases else {return Err(Error::Invalid("missing contact-side material state"));};
            let point=binding.point.at([a,b]).map_err(Error::State)?;
            let selection=match pins {
                Some(pins)=>pins.get(&binding.link).ok_or(Error::Invalid("prepared contact claim is absent"))?.clone(),
                None=>binding.selection.clone(),
            };
            let resolved=resolve_electrical_interface_state(binding.point.interface,&point,selection).map_err(Error::Electrical)?;
            let property=resolved.resolved().property(ELECTRICAL_CONTACT_RESISTANCE_PROPERTY)
                .ok_or(Error::Invalid("contact source property absent"))?;
            let claim=binding.point.interface.claims().claim(property.answer().receipt.selected)
                .ok_or(Error::Invalid("contact source claim absent"))?;
            match &claim.value {
                PropertyValue::Scalar {..}=>{},
                PropertyValue::Curve {abscissa,abscissa_dims,..} if binding.point.is_thermal_axis(abscissa,*abscissa_dims)=>{},
                _=>return Err(Error::Invalid("contact law must be static or use declared thermal/phase coordinates")),
            }
            if property.angular_frequency_band_rad_s().map_err(|e|Error::State(ThermalMaterialError::Coordinate(format!("{e:?}"))))?
                .is_some_and(|(lower,_)|lower>0.0) {
                return Err(Error::Invalid("electrical contact source excludes static response"));
            }
            contacts.push(SeriesElectricalContact::new(binding.link,binding.branch,resolved,link.area_m2).map_err(Error::Circuit)?);
        }
        if cancelled() {return Err(Error::Thermal(DepositionError::Cancelled));}
        Ok(ContactThermalCoefficients {bulk,contacts})
    }

    /// Data-only closed-loop contact heating on the shared thermal/circuit step.
    /// Initially selected bulk and contact sources cannot change in a trial.
    /// Contact/radiation heat transport can be composed through the generic
    /// `contacts::step_candidate` callback; this convenience path uses the
    /// original material network's bulk thermal links and authored boundaries.
    #[allow(clippy::too_many_arguments)]
    pub fn step_candidate<'m>(&self,thermal:&DepositionSimulation<'m>,target_s:f64,
        circuit:&MaterialCircuit,ports:&[f64],controls:CircuitEnthalpySolve,guess_j:&[f64],cancelled:impl FnMut()->bool)
        ->Result<CircuitEnthalpyCandidate<'m,ContactThermalSample<MaterialNetworkSample>>,CircuitEnthalpyError<ElectricalContactNetworkError>> {
        let cancel=RefCell::new(cancelled);
        if (*cancel.borrow_mut())() {return Err(CircuitEnthalpyError::Cancelled);}
        let routes=self.routes(thermal).map_err(CircuitEnthalpyError::Source)?;
        let prepared=match self.prepare(thermal,||(*cancel.borrow_mut())()) {
            Err(ElectricalContactNetworkError::Thermal(DepositionError::Cancelled))
            |Err(ElectricalContactNetworkError::Bulk(MaterialNetworkError::Thermal(DepositionError::Cancelled)))
                =>return Err(CircuitEnthalpyError::Cancelled),
            other=>other.map_err(CircuitEnthalpyError::Source)?,
        };
        let destinations=self.bulk.resistors.iter().map(|r|r.destination).collect::<Vec<_>>();
        let result=super::step_candidate(thermal,target_s,CircuitEnthalpyInput {circuit,ports,destinations:&destinations},
            &routes,controls,guess_j,|state|prepared.resolve(state,||(*cancel.borrow_mut())()),
            ||(*cancel.borrow_mut())());
        match result {
            Err(CircuitEnthalpyError::Source(ElectricalContactNetworkError::Thermal(DepositionError::Cancelled)))
            |Err(CircuitEnthalpyError::Source(ElectricalContactNetworkError::Bulk(MaterialNetworkError::Thermal(DepositionError::Cancelled))))
                =>Err(CircuitEnthalpyError::Cancelled),
            other=>other,
        }
    }

    /// Drive an electric or magnetic gap from this SAME bulk/contact/thermal
    /// material network. Actual endpoint winding/interface temperatures alter
    /// resistance, current, force and motion while all heat is deposited once.
    /// Field material and inertia retain the actuator's admitted fixed snapshots.
    #[allow(clippy::too_many_arguments)]
    pub fn field_motion_candidate<'m>(&self,thermal:&DepositionSimulation<'m>,target_s:f64,
        motion:&crate::material_circuit::motion::MaterialFieldMotion,ports:&[f64],external_force_n:f64,
        controls:CircuitEnthalpySolve,guess_j:&[f64],cancelled:impl FnMut()->bool)
        ->Result<crate::circuit_enthalpy::motion::FieldMotionEnthalpyCandidate<'m,ContactThermalSample<MaterialNetworkSample>>,
            CircuitEnthalpyError<ElectricalContactNetworkError>> {
        let cancel=RefCell::new(cancelled);
        if (*cancel.borrow_mut())() {return Err(CircuitEnthalpyError::Cancelled);}
        let routes=self.routes(thermal).map_err(CircuitEnthalpyError::Source)?;
        let prepared=match self.prepare(thermal,||(*cancel.borrow_mut())()) {
            Err(ElectricalContactNetworkError::Thermal(DepositionError::Cancelled))
            |Err(ElectricalContactNetworkError::Bulk(MaterialNetworkError::Thermal(DepositionError::Cancelled)))
                =>return Err(CircuitEnthalpyError::Cancelled),
            other=>other.map_err(CircuitEnthalpyError::Source)?,
        };
        let destinations=self.bulk.resistors.iter().map(|r|r.destination).collect::<Vec<_>>();
        let result=crate::circuit_enthalpy::motion::step_with_contacts_candidate(thermal,target_s,
            crate::circuit_enthalpy::motion::FieldMotionEnthalpyInput {motion,ports,external_force_n,destinations:&destinations},
            &routes,controls,guess_j,|state|prepared.resolve(state,||(*cancel.borrow_mut())()),
            ||(*cancel.borrow_mut())());
        match result {
            Err(CircuitEnthalpyError::Source(ElectricalContactNetworkError::Thermal(DepositionError::Cancelled)))
            |Err(CircuitEnthalpyError::Source(ElectricalContactNetworkError::Bulk(MaterialNetworkError::Thermal(DepositionError::Cancelled))))
                =>Err(CircuitEnthalpyError::Cancelled),
            other=>other,
        }
    }
}
