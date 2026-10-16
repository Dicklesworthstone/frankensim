//! Implicit wall-feedback transpose for the implemented mean-patch radiation.
//! The response can contain steady conduction or C/dt + J(T_new). Neither this
//! loop nor its IQN history is differentiated. Consistent Robin face integrals
//! remain owned by RobinResponse; radiation is never added to the air heat.
use super::*;
use fs_conduction::adjoint::robin::RobinResponse;
use fs_couple::iqn_ils::IqnIls;

/// Exact combined rows and radiation driving point for one accepted solid call.
/// These are solver inputs, not a claim that a finite-tolerance root is exact.
pub(super) struct Point {
    pub htc: BTreeMap<String, f64>,
    pub references: Vec<f64>,
    pub driving: Vec<f64>,
    pub radiative_htc: Vec<f64>,
}

pub(in crate::network_command) struct PatchGradient {
    pub surface: String,
    pub log_emissivity: f64,
    pub ambient_temperature: f64,
}

pub(in crate::network_command) struct Gradient {
    pub thermal: fan_gradient::CoolingGradient,
    pub patches: Vec<PatchGradient>,
    pub threshold: f64,
}

struct Row {
    air_fraction: f64,
    radiation_fraction: f64,
    air_offset: f64,
    ambient_offset: f64,
    wall_log_slope: f64,
    ambient_log_slope: f64,
}
impl Row {
    fn pullback(&self, reference: f64, log_htc: f64) -> Result<[f64; 5]> {
        let radiation = finite(self.radiation_fraction
            * (log_htc + self.ambient_offset * reference), "radiation coefficient pullback")?;
        Ok([
            finite(self.air_fraction * reference, "convective reference pullback")?,
            finite(self.air_fraction * (log_htc + self.air_offset * reference),
                "convective coefficient pullback")?,
            radiation,
            finite(self.radiation_fraction * reference + radiation * self.ambient_log_slope,
                "radiation surroundings pullback")?,
            finite(radiation * self.wall_log_slope, "radiation temperature feedback")?,
        ])
    }
}

/// d ln(h_rad)/d(T_wall,T_ambient), normalized before squaring.
fn coefficient_slopes(wall: f64, ambient: f64) -> Result<(f64, f64)> {
    if !(wall.is_finite() && ambient.is_finite() && wall > 0.0 && ambient > 0.0) {
        return Err(producer("radiation tangent requires positive finite temperatures"));
    }
    let scale = wall.max(ambient);
    let w = wall / scale;
    let a = ambient / scale;
    let common = 1.0 / (w + a);
    let denominator = w*w + a*a;
    Ok((finite((common + 2.0*w/denominator)/scale, "radiation wall slope")?,
        finite((common + 2.0*a/denominator)/scale, "radiation ambient slope")?))
}

/// Temperature functionals only. Heat-output derivatives require their own
/// direct terms and must not be smuggled in as combined-Robin heat weights.
#[allow(clippy::too_many_arguments)]
pub(super) fn pullback(
    policy: &Policy, request: &Request, cx: &Cx<'_>, network: &TransportNetwork<'_>,
    solid: &RobinResponse, point: &Point, references: &[f64], htc: &BTreeMap<String, f64>,
    nodal_weights: &[f64], wall_weights: &[f64], derivations: &[convection::Derived],
) -> Result<Gradient> {
    poll(cx)?;
    let names = network.regions();
    let n = names.len();
    if n == 0 || references.len() != n || point.references.len() != n
        || point.driving.len() != n || point.radiative_htc.len() != n
        || solid.ports().len() != n || wall_weights.len() != n
        || nodal_weights.len() != solid.temperature().len()
    {
        return Err(bad("radiation derivative requires the complete accepted interface and objective"));
    }
    let walls = solid.wall_means(cx, solid.temperature()).map_err(producer)?;
    let air = network.linearize(cx, &walls).map_err(producer)?;
    let mut rows = Vec::with_capacity(n);
    for (i, (&name, port)) in names.iter().zip(solid.ports()).enumerate() {
        poll(cx)?;
        let ha = *htc.get(name).ok_or_else(|| bad("missing radiation air coefficient"))?;
        let total = *point.htc.get(name).ok_or_else(|| bad("missing combined radiation coefficient"))?;
        let hr = point.radiative_htc[i];
        let air_error = finite(air.primal().reference_temperatures_k[i] - references[i],
            "radiation tangent air binding")?;
        if port.name != name || port.htc_w_m2_k.to_bits() != total.to_bits()
            || port.reference_k.to_bits() != point.references[i].to_bits()
            || air_error.abs() > request.limits.temperature
        {
            return Err(producer("radiation derivative is not bound to the accepted solid/air rows"));
        }
        let (ambient, wall_log_slope, ambient_log_slope) = match policy.patches.get(name) {
            Some(patch) => {
                let (w,a) = coefficient_slopes(point.driving[i], patch.ambient_k)?;
                (patch.ambient_k,w,a)
            }
            None => (port.reference_k,0.0,0.0),
        };
        rows.push(Row { air_fraction: ha/total, radiation_fraction: hr/total,
            air_offset: references[i]-port.reference_k, ambient_offset: ambient-port.reference_k,
            wall_log_slope, ambient_log_slope });
    }
    let mut current = vec![0.0; n];
    let zero_heat = vec![0.0; n];
    let mut accelerator = IqnIls::new(n, acceleration::POLICY).map_err(producer)?;
    let mut last_residual = 0.0;
    for iteration in 1..=request.limits.derivative {
        poll(cx)?;
        let weights = current.iter().zip(wall_weights)
            .map(|(a,b)| finite(a+b, "radiation objective seed")).collect::<Result<Vec<_>>>()?;
        let gradient = solid.pullback(cx, nodal_weights, &weights, &zero_heat).map_err(producer)?;
        let chained = rows.iter().enumerate().map(|(i,row)|
            row.pullback(gradient.references[i], gradient.log_htc[i])).collect::<Result<Vec<_>>>()?;
        let mut air_weights = air.zero_objective();
        air_weights.references = chained.iter().map(|row| row[0]).collect();
        let transported = air.pullback(cx, &air_weights).map_err(producer)?;
        let proposal = transported.walls.iter().zip(&chained)
            .map(|(a,row)| finite(a+row[4], "coupled radiation transpose proposal"))
            .collect::<Result<Vec<_>>>()?;
        let mut residual = 0.0_f64;
        let mut scale = 0.0_f64;
        for (&old,&next) in current.iter().zip(&proposal) {
            residual = residual.max(finite(next-old, "radiation transpose residual")?.abs());
            scale = scale.max(old.abs()).max(next.abs());
        }
        let threshold = finite(request.limits.relative + request.limits.relative*scale,
            "radiation transpose threshold")?;
        last_residual = residual;
        poll(cx)?;
        if residual <= threshold {
            let log_htc = transported.log_conductances.iter().zip(&chained)
                .map(|(air,row)| finite(air+row[1], "total radiative cooling coefficient derivative"))
                .collect::<Result<Vec<_>>>()?;
            // For temperature objectives, scaling every air hA and flow
            // capacity equally leaves air temperatures unchanged. Only the
            // air-conductance terms enter this homogeneity identity.
            let flow_scale = transported.log_conductances.iter().try_fold(0.0, |sum,value|
                finite(sum-value, "radiative cooling capacity derivative"))?;
            let thermal = CoupledGradient { inlets: transported.inlets, log_htc,
                nodal_load: gradient.nodal_load, interface_adjoint: air_weights.references,
                interface_residual: residual, iterations: iteration };
            let thermal = fan_gradient::from_flow_response(request, cx, &names, derivations, thermal, flow_scale)?;
            let patches = names.iter().enumerate().filter(|(_,name)| policy.patches.contains_key(**name))
                .map(|(i,&name)| PatchGradient { surface: name.to_string(),
                    log_emissivity: chained[i][2], ambient_temperature: chained[i][3] }).collect();
            return Ok(Gradient { thermal, patches, threshold });
        }
        if iteration < request.limits.derivative {
            current = accelerator.step(&current, &proposal, request.limits.relaxation)
                .map_err(producer)?.values;
        }
    }
    Err(Failure { code: "cooling-network-radiation-budget", message: format!(
        "total radiation adjoint exhausted {} sweeps; last unrelaxed equation residual {last_residual}; no partial gradient published",
        request.limits.derivative) })
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn secant_slopes_keep_both_temperature_paths_and_equal_temperature_limit() {
        for (w,a) in [(280.0,310.0),(500.0,250.0),(300.0,300.0)] {
            let (dw,da) = coefficient_slopes(w,a).unwrap();
            let h = |w:f64,a:f64| (w+a)*(w*w+a*a);
            let e = 0.001;
            assert!((dw-(h(w+e,a).ln()-h(w-e,a).ln())/(2.0*e)).abs()<1e-10);
            assert!((da-(h(w,a+e).ln()-h(w,a-e).ln())/(2.0*e)).abs()<1e-10);
            assert!((w*dw+a*da-3.0).abs()<1e-14);
        }
        let (dw,da)=coefficient_slopes(300.0,300.0).unwrap();
        assert!((dw-0.005).abs()<1e-16);
        assert_eq!(dw,da);
    }
    #[test]
    fn combined_robin_chain_keeps_air_and_radiation_controls_distinct() {
        let row=Row { air_fraction:0.8,radiation_fraction:0.2,air_offset:2.0,
            ambient_offset:-8.0,wall_log_slope:0.005,ambient_log_slope:0.004 };
        for (a,b) in row.pullback(3.0,-7.0).unwrap().into_iter().zip([2.4,-0.8,-6.2,0.5752,-0.031]) {
            assert!((a-b).abs()<1e-14);
        }
        assert!(coefficient_slopes(0.0,300.0).is_err());
        assert!(row.pullback(f64::INFINITY,0.0).is_err());
    }
}
