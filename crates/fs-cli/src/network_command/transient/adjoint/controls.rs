//! Opt-in physical controls of the same accepted trajectory. No extra primal
//! or adjoint solve per component/contact, and no division by baseline watts.
use super::*;
use fs_conduction::transient::backward_euler::StepLinearization;

// Bound both accumulator storage and the number of rows the result can emit.
const MAX_COMPONENT_INTERVAL_ROWS: usize = 65_536;

#[derive(Debug, Clone, Copy, Default)]
pub(super) struct Options {
    components: bool,
    contacts: bool,
}
impl Options {
    pub(super) fn parse(value: &J) -> Result<Self> {
        Ok(Self {
            components: value.get("component_power").map(|v| boolean(v,"adjoint.component_power"))
                .transpose()?.unwrap_or(false),
            contacts: value.get("contact_resistance").map(|v| boolean(v,"adjoint.contact_resistance"))
                .transpose()?.unwrap_or(false),
        })
    }

    fn dimensions(self, request: &Request, schedule: &Schedule) -> Result<(usize,usize)> {
        let components = if self.components {
            let map = request.solid_data.component_map.as_ref()
                .ok_or_else(|| bad("component-power adjoints require solid.component_power footprints"))?;
            let audit = request.solid_data.power.as_ref()
                .ok_or_else(|| bad("component-power adjoints require the admitted PowerMap audit"))?;
            if map.components().len() != audit.rows().len() {
                return Err(bad("component-power footprint/audit arity mismatch"));
            }
            for (component,row) in map.components().iter().zip(audit.rows()) {
                if component.name() != row.name() || !row.bound_volume_m3().is_finite()
                    || row.bound_volume_m3() <= 0.0 {
                    return Err(bad("component-power footprint normalization changed"));
                }
            }
            map.components().len()
        } else { 0 };
        let contacts = if self.contacts {
            request.contacts.as_ref().map(|c| c.interfaces.surface_count())
                .filter(|&n| n>0).ok_or_else(|| bad("contact-resistance adjoints require solid.contacts"))?
        } else { 0 };
        let entries = components.checked_mul(schedule.intervals.len())
            .ok_or_else(|| budget("component-control count overflow"))?;
        if entries > MAX_COMPONENT_INTERVAL_ROWS {
            return Err(budget("component-power adjoints exceed 65536 component-interval rows"));
        }
        Ok((components,contacts))
    }

    /// Retained accumulator payload and vector headers share the existing
    /// max_checkpoint_bytes allowance. Temporary FEM vectors are workspace.
    pub(super) fn storage_bytes(self, request: &Request, schedule: &Schedule) -> Result<usize> {
        let (components,contacts) = self.dimensions(request,schedule)?;
        let entries = components.checked_mul(schedule.intervals.len())
            .and_then(|n| n.checked_add(contacts))
            .ok_or_else(|| budget("trajectory-control count overflow"))?;
        entries.checked_mul(std::mem::size_of::<f64>())
            .and_then(|n| n.checked_add((usize::from(self.components)+usize::from(self.contacts))
                * std::mem::size_of::<Vec<f64>>()))
            .ok_or_else(|| budget("trajectory-control byte count overflow"))
    }
}

pub(super) struct Accumulation {
    components: Option<Vec<f64>>,
    contacts: Option<Vec<f64>>,
    component_count: usize,
}
fn zeros(count: usize) -> Result<Vec<f64>> {
    let mut values = Vec::new();
    values.try_reserve_exact(count).map_err(|_| budget("cannot allocate admitted trajectory controls"))?;
    values.resize(count,0.0);
    Ok(values)
}
impl Accumulation {
    pub(super) fn new(options: Options, request: &Request, schedule: &Schedule) -> Result<Self> {
        let (component_count,contact_count) = options.dimensions(request,schedule)?;
        let component_entries = component_count.checked_mul(schedule.intervals.len())
            .ok_or_else(|| budget("component-control count overflow"))?;
        Ok(Self {
            components: if options.components { Some(zeros(component_entries)?) } else { None },
            contacts: if options.contacts { Some(zeros(contact_count)?) } else { None },
            component_count,
        })
    }

    /// The total nodal-load adjoint has already passed the same coupled
    /// residual gate as history/fan/radiation controls. Its inverse includes
    /// C/dt; neither the source nor contact contraction is multiplied by dt.
    pub(super) fn record(&mut self, request: &Request, cx: &Cx<'_>,
        step: &StepLinearization<'_>, interval: usize, lambda: &[f64]) -> Result<()> {
        poll(cx)?;
        if let Some(values) = self.components.as_mut() {
            let density = step.source_density_pullback(cx,lambda).map_err(producer)?;
            let map = request.solid_data.component_map.as_ref().ok_or_else(|| bad("missing component map"))?;
            let audit = request.solid_data.power.as_ref().ok_or_else(|| bad("missing component audit"))?;
            let start = interval.checked_mul(self.component_count)
                .ok_or_else(|| budget("component interval offset overflow"))?;
            let end = start.checked_add(self.component_count)
                .ok_or_else(|| budget("component interval offset overflow"))?;
            let sums = values.get_mut(start..end).ok_or_else(|| bad("unknown component-control interval"))?;
            for ((sum,component),row) in sums.iter_mut().zip(map.components()).zip(audit.rows()) {
                poll(cx)?;
                let mut numerator = 0.0;
                for (index,&vertex) in component.vertices().iter().enumerate() {
                    if index % 512 == 0 { poll(cx)?; }
                    numerator = finite(numerator + *density.get(vertex)
                        .ok_or_else(|| bad("component footprint left the source-density field"))?)?;
                }
                // This is the original producer's geometric normalization,
                // not a reconstructed unit-power or zero-power ratio.
                *sum = finite(*sum + finite(numerator/row.bound_volume_m3())?)?;
            }
        }
        if let Some(sums) = self.contacts.as_mut() {
            let contacts = request.contacts.as_ref().ok_or_else(|| bad("missing trajectory contacts"))?;
            let values = contacts.trajectory_log_resistance_gradients(cx,&step.primal().temperature,lambda)?;
            if values.len() != sums.len() { return Err(bad("contact-control arity changed")); }
            for (sum,value) in sums.iter_mut().zip(values) { *sum = finite(*sum+value)?; }
        }
        poll(cx)
    }

    /// No fragment at all when neither control was requested: ordinary
    /// adjoints keep their existing result and do no additional integration.
    pub(super) fn report_fragment(&self, request: &Request, cx: &Cx<'_>, schedule: &Schedule)
        -> Result<String> {
        if self.components.is_none() && self.contacts.is_none() { return Ok(String::new()); }
        let components = match &self.components {
            Some(values) => self.component_report(request,cx,schedule,values)?,
            None => "null".into(),
        };
        let contacts = match &self.contacts {
            Some(values) => request.contacts.as_ref().ok_or_else(|| bad("missing trajectory contacts"))?
                .trajectory_sensitivity_json(values)?,
            None => "null".into(),
        };
        Ok(format!(",\"component_power_sensitivities\":{components},\"contact_resistance_sensitivities\":{contacts}"))
    }

    fn component_report(&self, request: &Request, cx: &Cx<'_>, schedule: &Schedule,
        values: &[f64]) -> Result<String> {
        let map = request.solid_data.component_map.as_ref().ok_or_else(|| bad("missing component map"))?;
        let mut base = zeros(self.component_count)?;
        let mut intervals = Vec::new();
        for (ordinal,interval) in schedule.intervals.iter().enumerate() {
            poll(cx)?;
            let mut rows = Vec::new();
            for (index,component) in map.components().iter().enumerate() {
                let derivative = values[ordinal*self.component_count+index];
                let watts = match &interval.workload {
                    Workload::Scale(scale) => {
                        base[index] = finite(base[index]+finite(scale*derivative)?)?;
                        finite(scale*component.watts())?
                    }
                    Workload::Components(powers) => *powers.get(component.name())
                        .ok_or_else(|| bad("workload omits a component control"))?,
                };
                rows.push(format!("{{\"component\":{},\"applied_power_w\":{},\"dtemperature_dpower_w_k_per_w\":{},\"dtemperature_dpower_multiplier_k\":{}}}",
                    quote(component.name()),num(watts)?,num(derivative)?,num(finite(watts*derivative)?)?));
            }
            intervals.push(format!("{{\"interval\":{ordinal},\"rows\":[{}]}}",rows.join(",")));
        }
        let base_rows = map.components().iter().zip(base).map(|(component,derivative)| {
            Ok(format!("{{\"component\":{},\"base_power_w\":{},\"dtemperature_dbase_power_w_k_per_w\":{}}}",
                quote(component.name()),num(component.watts())?,num(derivative)?))
        }).collect::<Result<Vec<_>>>()?.join(",");
        Ok(format!("{{\"method\":\"consistent-P1-source-transpose\",\"intervals\":[{}],\"base_rows\":[{base_rows}],\"scope\":\"interval derivative changes one component's applied watts at every occurrence of that base interval; fixed original P1 footprint, including overlapping and zero-power components; an absolute interval control can be replayed by component_powers_w with every other component unchanged; base_rows change original solid.component_power watts only in power_scale intervals, not explicitly overridden component_powers_w intervals; update declared system total with base watts; no renormalization of other components, extra timestep factor, continuous-time peak bound or unique derivative at maximum ties; at zero watts the admissible direction is one-sided\"}}",intervals.join(",")))
    }
}
