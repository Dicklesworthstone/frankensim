//! Rebind the nominal hydraulic point from the SAME retained branch capacities.
//! There is no replacement flow solve and no report-JSON reconstruction.
use super::{Cx, ProjectSpec, SolveRefusal, bad, conduction_error, finite, poll,
    AIR_SPECIFIC_HEAT_J_KG_K, conjugate};
use fs_airflow::{FanArrangement, EnclosureNetwork, LeakageElement, LossNetwork,
    sharp_edged_orifice_loss, conjugate::AirPath};
use fs_airflow::composite::speed::{operating_flow_speed_gradient, SpeedDerivativeError};
use fs_project::fansystem::{FanSystemDecl, FanSystemTopology, lower_fan_system};
use fs_project::spec::dims;
use fs_qty::{Area, Density, VolumetricFlowRate};
use std::collections::BTreeMap;

/// None is a nonsmooth physical control, never a zero hydraulic derivative.
pub(super) fn flow_weights(
    cx: &Cx<'_>, spec: &ProjectSpec, paths: &[AirPath], laws: &[conjugate::AirflowLaw],
) -> Result<Option<Vec<f64>>, SolveRefusal> {
    poll(cx)?;
    let cooling=spec.cooling.as_ref().ok_or_else(|| bad("missing hydraulic cooling declaration"))?;
    let system=cooling.fan_system.as_ref().ok_or_else(|| bad("missing fan system"))?;
    let topology=match &system.topology {
        FanSystemTopology::Series(_) => FanArrangement::Series,
        FanSystemTopology::Parallel(_) => FanArrangement::Parallel,
        FanSystemTopology::Single => return Err(bad("multi-bank derivative requires the declared composition")),
    };
    let memory=spec.budgets.as_ref().map_or(0,|b| b.memory_bytes);
    let max_points=usize::try_from(memory/128).unwrap_or(usize::MAX).min(8192);
    let mut points=0_usize;
    for bank in &system.banks {
        poll(cx)?;
        points=points.checked_add(bank.curve.points.len()).ok_or_else(|| bad("fan curve count overflow"))?;
    }
    if system.banks.len()>64 || points>max_points || cooling.vents.is_empty() || cooling.vents.len()>64 {
        return Err(bad("multi-bank derivative exceeds 64 banks/vents or its declared curve-point budget"));
    }
    // Reuse the project owner's validated lowering, one member at a time.
    // Building a second composite curve would add needless quadratic work.
    let mut banks=Vec::with_capacity(system.banks.len());
    for bank in &system.banks {
        poll(cx)?;
        banks.push(lower_fan_system(&FanSystemDecl { version:system.version,
            banks:vec![bank.clone()],topology:FanSystemTopology::Single })
            .map_err(|e| bad(format!("fan member lowering refused: {}",e.detail)))?.system_bank);
    }
    let envelope=spec.envelope.as_ref().ok_or_else(|| bad("missing hydraulic density envelope"))?;
    // Same arithmetic and physical constants as flow_network_receipt.
    let ambient_mid=(envelope.ambient_lo.value+envelope.ambient_hi.value)/2.0;
    let rho=finite(envelope.pressure.value/(super::super::super::AIR_SPECIFIC_GAS_CONSTANT*ambient_mid))?;
    if rho<=0.0 { return Err(bad("nonpositive hydraulic density")); }
    let density=Density::new(rho);
    let mut resistances=BTreeMap::new();
    let mut branches=Vec::with_capacity(cooling.vents.len());
    for vent in &cooling.vents {
        poll(cx)?;
        if vent.area.dims!=dims::AREA { return Err(bad("vent area units changed")); }
        let element=sharp_edged_orifice_loss(format!("vent:{}",vent.region),Area::new(vent.area.value),density)
            .map_err(|e| bad(format!("native vent law refused: {e}")))?;
        if resistances.insert(vent.region.as_str(),element.resistance.value()).is_some() {
            return Err(bad("ambiguous native vent branch ownership"));
        }
        branches.push(LossNetwork::Element(element));
    }
    let primary=if branches.len()==1 { branches.pop().expect("nonempty vents") }
        else { LossNetwork::parallel(branches).map_err(|e| bad(format!("native vent topology refused: {e}")))? };
    let leakage=cooling.airflow_leakage.as_ref().ok_or_else(|| bad("missing mandatory hydraulic leakage"))?;
    if leakage.area.dims!=dims::AREA { return Err(bad("leakage area units changed")); }
    let leak=sharp_edged_orifice_loss("leakage",Area::new(leakage.area.value),density)
        .map_err(|e| bad(format!("native leakage law refused: {e}")))?;
    let resistance=EnclosureNetwork::new(primary,LeakageElement::new(leak)).equivalent_resistance();
    if !(resistance.value().is_finite() && resistance.value()>0.0) {
        return Err(bad("invalid equivalent hydraulic resistance"));
    }
    let mut start=0_usize;
    let mut total=None::<f64>;
    for path in paths {
        poll(cx)?;
        let end=start.checked_add(path.segments().len()).ok_or_else(|| bad("air path count overflow"))?;
        if path.segments().is_empty() || end>laws.len() { return Err(bad("missing native air-path laws")); }
        let branch=&laws[start].branch;
        if path.segments().iter().zip(&laws[start..end]).any(|(s,l)| s.region()!=l.target.as_str() || l.branch.as_str()!=branch.as_str()) {
            return Err(bad("hydraulic derivative does not own the retained branch ordering"));
        }
        let r=*resistances.get(branch.as_str()).ok_or_else(|| bad("air path names no native vent"))?;
        let branch_flow=finite((path.capacity_rate_w_per_k()/AIR_SPECIFIC_HEAT_J_KG_K)/rho)?;
        // All native vents and leakage are parallel quadratic paths, so
        // q_branch/Q=sqrt(R_total/R_branch). Keep leakage in R_total.
        let q=finite(branch_flow*fs_math::det::sqrt(finite(r/resistance.value())?))?;
        if q<=0.0 { return Err(bad("nonpositive retained hydraulic flow")); }
        if let Some(previous)=total {
            if (q-previous).abs()>1024.0*f64::EPSILON*q.abs().max(previous.abs()) {
                return Err(bad("retained air capacities do not share the declared hydraulic operating point"));
            }
        } else { total=Some(q); }
        start=end;
    }
    if start!=laws.len() { return Err(bad("hydraulic derivative omits declared air segments")); }
    let q=total.ok_or_else(|| bad("no retained hydraulic flow"))?;
    // The owner rechecks the nominal fan/loss equations at this exact field's
    // reconstructed point. This is not differentiation of tolerance corners.
    match operating_flow_speed_gradient(cx,&banks,topology,VolumetricFlowRate::new(q),resistance,1e-8,max_points) {
        Ok(g) => Ok(Some(g.log_flow_per_log_speed)),
        Err(SpeedDerivativeError::NonSmooth { .. }) => Ok(None),
        Err(SpeedDerivativeError::Cancelled) => Err(conduction_error("cli-solve-cancelled",
            "multi-bank hydraulic derivative interrupted","resume the accepted pipeline prefix")),
        Err(e) => Err(bad(format!("multi-bank hydraulic derivative refused: {e}"))),
    }
}
