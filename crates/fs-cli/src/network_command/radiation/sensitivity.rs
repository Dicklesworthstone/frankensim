//! Total derivative of the implemented mean-patch radiation closure.
//!
//! Let w be all wall means, r=R(w,p) the air references, and T=S(r,h_rad(w,p),p)
//! the existing combined-Robin solid solve. The implicit transpose equation is
//! mu = (d inputs/dw)^T S'^T (objective_T + C^T mu). Each sweep uses ONE actual
//! FEM pullback and ONE transport pullback. No dense interface matrix, finite
//! differences, pointwise-T^4 substitution or solver-iteration differentiation.
//!
//! The nominal state is the accepted numerical fixed point. Partial derivatives
//! are evaluated at the exact combined rows that produced its final field; their
//! difference from an exact-root derivative is not an outward-rounded bound.
use super::*;
use fs_airflow::graph::thermal::coupled_transport::sensitivity::CoupledObjective;
use fs_conduction::adjoint::robin::RobinGradient;
use fs_couple::iqn_ils::IqnIls;

/// Exact inputs to the final solid evaluation, retained only on a gradient run.
/// Repeating these inputs must reproduce its field before binding derivatives.
pub(super) struct Binding {
    pub(super) htc: BTreeMap<String, f64>,
    pub(super) references: Vec<f64>,
    pub(super) driving: Vec<f64>,
    pub(super) radiative_htc: Vec<f64>,
    pub(super) config: SolveConfig,
}

/// Chain rule for one actual combined Robin row. Heat is integrated with the
/// producer's consistent face mass, not a product of average temperature jumps.
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

/// d ln(h_rad)/d(T_wall,T_ambient), with normalization before squaring.
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

#[allow(clippy::too_many_arguments)]
pub(super) fn pullback(
    policy: &Policy, request: &Request, cx: &Cx<'_>, network: &TransportNetwork<'_>,
    inner: &Inner, references: &[f64], htc: &BTreeMap<String, f64>,
    objective: &objective::ObjectiveState, derivations: &[convection::Derived],
) -> Result<(fan_gradient::CoolingGradient, String)> {
    poll(cx)?;
    let names = network.regions();
    let n = names.len();
    let binding = inner.binding.as_ref().ok_or_else(|| bad("missing radiation derivative binding"))?;
    if n == 0 || references.len() != n || binding.references.len() != n
        || binding.driving.len() != n || binding.radiative_htc.len() != n {
        return Err(bad("radiation derivative requires the complete final interface"));
    }
    let boundary = request.boundary(&names, &binding.references, &binding.htc)?;
    let material = fs_conduction::ConductivityModel::isotropic_declared(request.conductivity).map_err(producer)?;
    let uniform = ScalarField::Uniform(request.source);
    let source = request.solid_data.nodal_source.as_ref().unwrap_or(&uniform);
    let problem = ConductionProblem { mesh: &request.mesh, boundary: &boundary,
        material: &material, element_materials: request.solid_data.element_materials.as_ref(), source };
    let solid = match &request.contacts {
        Some(contact) => RobinLinearization::new_with_interfaces(cx, problem, &contact.interfaces,
            binding.config.clone(), &names),
        None => RobinLinearization::new(cx, problem, binding.config.clone(), &names),
    }.map_err(producer)?;
    if solid.temperature().len() != inner.conduction.temperature.len()
        || solid.temperature().iter().zip(&inner.conduction.temperature)
            .any(|(a,b)| a.to_bits() != b.to_bits()) {
        return Err(producer("radiation adjoint reconstruction did not reproduce the accepted solid field"));
    }
    let walls = solid.wall_means(cx, solid.temperature()).map_err(producer)?;
    let air = network.linearize(cx, &walls).map_err(producer)?;
    let mut rows = Vec::with_capacity(n);
    for (i, (&name, port)) in names.iter().zip(solid.ports()).enumerate() {
        poll(cx)?;
        if port.name != name || (air.primal().reference_temperatures_k[i] - references[i]).abs()
            > request.limits.temperature {
            return Err(producer("radiation derivative is not bound to the accepted air fixed point"));
        }
        let ha = htc[name];
        let hr = binding.radiative_htc[i];
        let total = port.htc_w_m2_k;
        let (ambient, wall_log_slope, ambient_log_slope) = match policy.patches.get(name) {
            Some(patch) => {
                let (w,a) = coefficient_slopes(binding.driving[i], patch.ambient_k)?;
                (patch.ambient_k,w,a)
            }
            None => (port.reference_k,0.0,0.0),
        };
        rows.push(Row { air_fraction: ha/total, radiation_fraction: hr/total,
            air_offset: references[i]-port.reference_k, ambient_offset: ambient-port.reference_k,
            wall_log_slope, ambient_log_slope });
    }
    let mut seed = CoupledObjective { nodal_temperatures: vec![0.0; solid.temperature().len()],
        wall_temperatures: vec![0.0; n], solid_heat_rates: vec![0.0; n], air: air.zero_objective() };
    objective.seed(&mut seed);
    let mut current = vec![0.0; n];
    let mut accelerator = IqnIls::new(n, acceleration::POLICY).map_err(producer)?;
    let mut last_residual = 0.0;
    for iteration in 1..=request.limits.derivative {
        poll(cx)?;
        let weights = current.iter().zip(&seed.wall_temperatures)
            .map(|(a,b)| finite(a+b, "radiation objective seed")) .collect::<Result<Vec<_>>>()?;
        let RobinGradient { references: reference_gradient, log_htc: coefficient_gradient,
            nodal_load, .. } = solid.pullback(cx, &seed.nodal_temperatures, &weights, &seed.solid_heat_rates)
            .map_err(producer)?;
        let chained = rows.iter().enumerate().map(|(i,row)|
            row.pullback(reference_gradient[i], coefficient_gradient[i])).collect::<Result<Vec<_>>>()?;
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
            // The explicit objective contains temperatures only. Uniformly
            // scaling hA AND every flow capacity leaves air temperatures fixed.
            // Subtract only AIR conductance derivatives, never solid ones.
            let flow_scale = transported.log_conductances.iter().try_fold(0.0, |sum,value|
                finite(sum-value, "radiative cooling capacity derivative"))?;
            let thermal = CoupledGradient { inlets: transported.inlets, log_htc, nodal_load,
                interface_adjoint: air_weights.references, interface_residual: residual, iterations: iteration };
            let gradient = fan_gradient::from_flow_response(request, cx, &names, derivations, thermal, flow_scale)?;
            let mut patches = Vec::new();
            for (i,&name) in names.iter().enumerate() {
                if let Some(patch) = policy.patches.get(name) {
                    patches.push(format!("{{\"surface\":{},\"dobjective_dlog_emissivity_k\":{},\"dobjective_demissivity_k\":{},\"dobjective_dambient_temperature\":{}}}",
                        quote(name),num(chained[i][2])?,num(finite(chained[i][2]/patch.emissivity.value(),
                            "absolute emissivity derivative")?)?,num(chained[i][3])?));
                }
            }
            let report = format!("{{\"method\":\"implicit-wall-feedback-adjoint\",\"iterations\":{iteration},\"equation_residual\":{},\"equation_threshold\":{},\"reconstruction_solid_solves\":1,\"surfaces\":[{}],\"scope\":\"total steady derivative of the admitted mean-patch Robin closure, including k-prime, contacts, radiation-coefficient feedback and mixed-air feedback; full consistent face integrals, not products of mean jumps; one exact-input reconstruction must reproduce the accepted field bits; no iteration-history AD, pointwise T^4 substitution, transient derivative or error certificate; maximum objectives retain selected-active-vertex semantics; emissivity bounds restrict admissible perturbation directions\"}}",
                num(residual)?,num(threshold)?,patches.join(","));
            poll(cx)?;
            return Ok((gradient,report));
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
    fn secant_slopes_keep_both_temperature_paths_and_the_equal_temperature_limit() {
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
    fn combined_robin_chain_distinguishes_air_and_radiation_controls() {
        let row=Row { air_fraction:0.8,radiation_fraction:0.2,air_offset:2.0,
            ambient_offset:-8.0,wall_log_slope:0.005,ambient_log_slope:0.004 };
        let g=row.pullback(3.0,-7.0).unwrap();
        for (a,b) in g.into_iter().zip([2.4,-0.8,-6.2,0.5752,-0.031]) {
            assert!((a-b).abs()<1e-14);
        }
        assert!(coefficient_slopes(0.0,300.0).is_err());
        assert!(row.pullback(f64::INFINITY,0.0).is_err());
    }
}
