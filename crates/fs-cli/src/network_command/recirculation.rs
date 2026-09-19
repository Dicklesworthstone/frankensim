//! File-driven imposed exhaust return, sharing the actual transport producer.
//! Fresh temperatures remain caller inputs; mixed intake temperatures are
//! solved afresh after every change of flow, convection or solid temperature.
use super::*;
use fs_airflow::graph::thermal::transport::{TransportMarch, recirculation::{
    MAX_RECIRCULATION_LINKS, MAX_RECIRCULATION_SUPPLIES, RecirculationLink,
}};

mod sensitivity;
mod design;

const MODEL: &str = "prescribed-adiabatic-return";
const SCOPE: &str = "imposed receiving-supply fractions; positive fresh makeup; instantaneous adiabatic return with externally imposed pressure reset; no return-duct hydraulics, fan heat, humidity or residence time; not an experimentally validated or uncertainty-certified model";

#[derive(Debug)]
pub(super) struct Policy {
    links: Vec<RecirculationLink>,
    tolerance_k: f64,
    source: String,
    design: Option<design::Design>,
}

impl Policy {
    pub(super) fn parse_request(value: &J, inlets: &[TransportInlet], root: &J) -> Result<Self> {
        let Some(spec) = value.get("design") else { return Self::parse(value, inlets); };
        let fields = value.as_object().ok_or_else(|| bad("recirculation must be an object"))?;
        let plain = J::Object(fields.iter().filter(|(key, _)| key != "design").cloned().collect());
        let mut policy = Self::parse(&plain, inlets)?;
        policy.design = Some(design::Design::parse(spec, root, &policy)?);
        Ok(policy)
    }

    pub(super) fn run_design(&self, cx: &Cx<'_>) -> Option<Result<String>> {
        self.design.as_ref().map(|design| design.solve(cx))
    }

    pub(super) fn parse(value: &J, inlets: &[TransportInlet]) -> Result<Self> {
        object(value, &["model", "source", "temperature_tolerance_k", "links"], "recirculation")?;
        if get(value, "model")?.as_str() != Some(MODEL) {
            return Err(bad("recirculation.model must be prescribed-adiabatic-return"));
        }
        let source = string(get(value, "source")?, "recirculation.source")?;
        let tolerance_k = positive(get(value, "temperature_tolerance_k")?,
            "recirculation.temperature_tolerance_k")?;
        let entries = array(get(value, "links")?, "recirculation.links", MAX_RECIRCULATION_LINKS)?;
        if entries.is_empty() { return Err(bad("recirculation.links must be nonempty")); }
        let mut links = Vec::with_capacity(entries.len());
        for entry in entries {
            object(entry, &["supply_node", "return_node", "fraction"], "recirculation link")?;
            let supply_node = integer_raw(get(entry, "supply_node")?, "recirculation.supply_node")?;
            let return_node = integer_raw(get(entry, "return_node")?, "recirculation.return_node")?;
            let fraction = number(get(entry, "fraction")?, "recirculation.fraction")?;
            if !inlets.iter().any(|inlet| inlet.node == supply_node)
                || supply_node == return_node || !(0.0..1.0).contains(&fraction) {
                return Err(bad("return links require a declared fresh supply, a distinct exhaust node, and a fraction in [0,1)"));
            }
            links.push(RecirculationLink { supply_node, return_node, fraction });
        }
        links.sort_by_key(|link| (link.supply_node, link.return_node));
        let mut fractions = BTreeMap::<usize, f64>::new();
        for (i, link) in links.iter().enumerate() {
            if i > 0 && (links[i - 1].supply_node, links[i - 1].return_node)
                == (link.supply_node, link.return_node) {
                return Err(bad("duplicate recirculation supply/return pair"));
            }
            if link.fraction == 0.0 { continue; }
            let sum = fractions.entry(link.supply_node).or_default();
            *sum += link.fraction;
            if !sum.is_finite() || *sum >= 1.0 {
                return Err(bad("each recycled supply requires a strictly positive fresh-air fraction"));
            }
        }
        if fractions.len() > MAX_RECIRCULATION_SUPPLIES {
            return Err(bad("recirculation supply cap exceeded"));
        }
        Ok(Self { links, tolerance_k, source, design: None })
    }

    pub(super) fn bind<'flow>(&self, cx: &Cx<'_>, network: TransportNetwork<'flow>)
        -> Result<TransportNetwork<'flow>> {
        poll(cx)?;
        if self.design.is_some() {
            return Err(bad("return-fraction designs must construct a resolved candidate before transport"));
        }
        // Recheck ACTUAL supply/exhaust signs and capacities on every candidate.
        network.with_recirculation(cx, self.links.clone(), self.tolerance_k).map_err(producer)
    }

    fn report(&self, march: &TransportMarch) -> Result<String> {
        let mut configured = Vec::with_capacity(self.links.len());
        for link in &self.links {
            configured.push(format!("{{\"supply_node\":{},\"return_node\":{},\"fraction\":{}}}",
                link.supply_node, link.return_node, num(link.fraction)?));
        }
        let prefix = format!("{{\"model\":{},\"source\":{},\"scope\":{},\"temperature_tolerance_k\":{},\"configured_links\":[{}]",
            quote(MODEL), quote(&self.source), quote(SCOPE), num(self.tolerance_k)?, configured.join(","));
        let Some(report) = &march.recirculation else {
            if self.links.iter().any(|link| link.fraction > 0.0) {
                return Err(producer("return-air configuration was lost before result publication"));
            }
            return Ok(format!("{prefix},\"status\":\"once-through-zero-returns\",\"supplies\":[],\"streams\":[]}}"));
        };
        let active: Vec<_> = self.links.iter().filter(|link| link.fraction > 0.0).collect();
        if active.len() != report.streams.len() || active.iter().zip(&report.streams).any(|(link, stream)| {
            link.supply_node != stream.supply_node || link.return_node != stream.return_node
                || link.fraction.to_bits() != stream.fraction.to_bits()
        }) {
            return Err(producer("return-air result does not match its configured topology"));
        }
        let supplies = report.supplies.iter().map(|s| Ok(format!(
            "{{\"node\":{},\"fresh_temperature_k\":{},\"mixed_temperature_k\":{},\"fresh_capacity_w_per_k\":{},\"return_capacity_w_per_k\":{}}}",
            s.node, num(s.fresh_temperature_k)?, num(s.mixed_temperature_k)?,
            num(s.fresh_capacity_w_per_k)?, num(s.return_capacity_w_per_k)?)))
            .collect::<Result<Vec<_>>>()?.join(",");
        let streams = report.streams.iter().map(|s| Ok(format!(
            "{{\"supply_node\":{},\"return_node\":{},\"fraction\":{},\"capacity_w_per_k\":{},\"temperature_k\":{}}}",
            s.supply_node, s.return_node, num(s.fraction)?, num(s.capacity_w_per_k)?, num(s.temperature_k)?)))
            .collect::<Result<Vec<_>>>()?.join(",");
        Ok(format!("{prefix},\"status\":\"solved\",\"supplies\":[{supplies}],\"streams\":[{streams}],\"external_heat_gain_w\":{},\"heat_imbalance_w\":{},\"max_mixing_residual_k\":{},\"inlet_gradient_semantics\":\"derivatives with respect to fresh makeup, including implicit return feedback; return fractions held fixed\"}}",
            num(report.external_heat_gain_w)?, num(report.heat_imbalance_w)?, num(report.max_mixing_residual_k)?))
    }
}

pub(super) fn attach(output: String, request: &Request, evaluated: &Evaluation) -> Result<String> {
    let march = &evaluated.coupled.transport;
    let (report, sensitivities) = match &request.recirculation {
        Some(policy) => (policy.report(march)?, policy.sensitivity_json(request, evaluated)?),
        None if march.recirculation.is_none() => return Ok(output),
        None => return Err(producer("unrequested return-air model in cooling result")),
    };
    let prefix = output.strip_suffix("}\n").ok_or_else(|| bad("invalid cooling JSON framing"))?;
    Ok(format!("{prefix},\"recirculation\":{report},\"recirculation_sensitivity\":{sensitivities}}}\n"))
}

/// Actual outer sensible-heat gain, excluding air circulating internally.
pub(super) fn external_heat_gain(march: &TransportMarch) -> f64 {
    march.recirculation.as_ref().map_or(march.external_heat_gain_w, |r| r.external_heat_gain_w)
}

/// A bounded endpoint observation; do not repeat topology/provenance per step.
pub(super) fn history_field(request: &Request, march: &TransportMarch) -> Result<String> {
    let Some(policy) = &request.recirculation else {
        if march.recirculation.is_some() { return Err(producer("unrequested transient return model")); }
        return Ok(String::new());
    };
    let Some(report) = &march.recirculation else {
        if policy.links.iter().any(|link| link.fraction > 0.0) {
            return Err(producer("return-air feedback missing from an accepted endpoint"));
        }
        return Ok(",\"recirculation\":{\"status\":\"once-through-zero-returns\"}".into());
    };
    let supplies = report.supplies.iter().map(|s| Ok(format!(
        "{{\"node\":{},\"mixed_temperature_k\":{}}}", s.node, num(s.mixed_temperature_k)?)))
        .collect::<Result<Vec<_>>>()?.join(",");
    Ok(format!(",\"recirculation\":{{\"mixed_supplies\":[{supplies}],\"external_heat_gain_w\":{},\"heat_imbalance_w\":{},\"max_mixing_residual_k\":{}}}",
        num(report.external_heat_gain_w)?, num(report.heat_imbalance_w)?, num(report.max_mixing_residual_k)?))
}

#[cfg(test)]
mod tests;
#[cfg(test)]
mod transient_tests;
