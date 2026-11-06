//! One compliant contact restraining the original assembly's single null mode.
//! K q = f - b R and n^T K = 0 imply R = (n^T f)/(n^T b). The existing
//! Schur response solves the compatible deformation; translating along n then
//! imposes the original contact closure. No penalty support or load projection.
use super::*;

pub(in super::super) fn initialize(
    network: &mut CoupledModalSystem, external: &[f64], contact: &ModalContact,
    config: ModalContactConfig, direction: &[f64], gate: &CancelGate,
) -> Result<f64, ModalCouplingError> {
    poll(Some(gate))?;
    network.require_static_initialization(external)?;
    if direction.len() != network.mode_count() {
        return Err(invalid("static support needs a complete mass-normalized null direction"));
    }
    admit_setup(network, 1, MultiContactConfig {
        max_contacts: 1, max_sweeps: 1, max_setup_terms: network.config.max_setup_terms,
    })?;
    let column = contact_column(network, contact, config)?;
    // Admits normalization, EVERY spring and elastic mode, and complete rank.
    // In particular an unrelated free coordinate cannot be quietly fixed.
    let response = StaticResponse::new_with_gauges(network, &[direction.to_vec()], gate)?;
    let storage = ContactStorage::new(Box::new(ZeroStorage), 1, vec![contact.law.clone()])
        .map_err(ModalCouplingError::ContactLaw)?;
    let point = Point { column, storage, config, gap: contact.law.gaps()[0] };
    let roundoff_factor = network.config.solve_relative_tolerance
        .min(64.0*f64::EPSILON*network.mode_count() as f64);
    let (coupling, coupling_scale) = projection(direction, &point.column)?;
    if coupling.abs() <= roundoff_factor*coupling_scale {
        return Err(invalid("contact cannot resolve the declared free direction"));
    }
    let (load, _) = projection(direction, external)?;
    let reaction = finite(load/coupling)?;
    if reaction < 0.0 {
        return Err(invalid("a unilateral contact cannot balance a tensile required reaction"));
    }
    if load != 0.0 && reaction == 0.0 {
        return Err(invalid("static support reaction is below representable range"));
    }
    limit("static contact force", reaction, config.maximum_force_n)?;
    let mut forces = Vec::with_capacity(external.len());
    let mut roundoff = Vec::with_capacity(external.len());
    for (&f, &b) in external.iter().zip(&point.column) {
        poll(Some(gate))?;
        let mapped = finite(b*reaction)?;
        if b != 0.0 && reaction != 0.0 && mapped == 0.0 {
            return Err(invalid("static support force is below representable range"));
        }
        // Retain the original subtraction, including its rounding residual.
        forces.push(finite(f-mapped)?);
        roundoff.push(finite(roundoff_factor*finite(f.abs()+mapped.abs())?)?);
    }
    let mut q = response.solve_with_roundoff(network, &forces, true, true, Some(&roundoff), gate)?;
    let closure = dot(&point.column, &q)?;
    let (target, iterations) = if reaction == 0.0 {
        // There is no unique pose under zero resultant load. Keep the admitted
        // zero-gauge solution if open; otherwise select its touching translate.
        (closure.min(point.gap), 0)
    } else {
        closure_for_reaction(&point, reaction, gate)?
    };
    translate(&mut q, direction, finite((target-closure)/coupling)?)?;
    // A stiff law can expose a rounded attachment contraction. One correction
    // uses the measured closure; it changes neither the law nor its tolerance.
    let actual = dot(&point.column, &q)?;
    let expected = point.force(actual)?;
    if finite(reaction-expected)?.abs() > force_tolerance(reaction, expected, config)? {
        translate(&mut q, direction, finite((target-actual)/coupling)?)?;
    }
    let energy = network.stage_static_equilibrium_with_allowance(&q, &forces, Some(&roundoff), gate)?;
    publish(network, external, std::slice::from_ref(&point), None, &[reaction], energy, iterations, gate)
}

pub(super) fn projection(a: &[f64], b: &[f64]) -> Result<(f64, f64), ModalCouplingError> {
    let mut sum = 0.0;
    let mut scale = 0.0;
    for (&a, &b) in a.iter().zip(b) {
        let term = finite(a*b)?;
        if a != 0.0 && b != 0.0 && term == 0.0 {
            return Err(invalid("static support projection is below representable range"));
        }
        sum = finite(sum+term)?;
        scale = finite(scale+term.abs())?;
    }
    Ok((sum, scale))
}

fn translate(q: &mut [f64], direction: &[f64], distance: f64) -> Result<(), ModalCouplingError> {
    for (q, &a) in q.iter_mut().zip(direction) { *q = finite(*q+a*distance)?; }
    Ok(())
}

// Invert the ACTUAL ContactStorage gradient, not a fitted or duplicated law.
// Contact closure is monotone; both the bracket and scalar work are bounded by
// the original penetration/iteration caps. Positive reactions never select the
// zero-force endpoint merely because an absolute tolerance is larger.
pub(super) fn closure_for_reaction(point: &Point, reaction: f64, gate: &CancelGate)
    -> Result<(f64, usize), ModalCouplingError>
{
    let mut lo = point.gap;
    let mut hi = finite(point.gap+point.config.maximum_penetration_m)?;
    let upper = point.force(hi)?;
    if upper < reaction || hi == lo {
        return Err(invalid("required contact reaction exceeds the representable penetration budget"));
    }
    let mut best = (hi, finite(upper-reaction)?, force_tolerance(reaction, upper, point.config)?);
    if best.1.abs() <= 0.0625*best.2 { return Ok((hi, 0)); }
    let mut iterations = 0;
    for iteration in 1..=point.config.max_iterations {
        poll(Some(gate))?;
        let x = f64::midpoint(lo, hi);
        if x == lo || x == hi { break; }
        iterations = iteration;
        let force = point.force(x)?;
        let residual = finite(force-reaction)?;
        let tolerance = force_tolerance(reaction, force, point.config)?;
        if force > 0.0 && residual.abs()/tolerance < best.1.abs()/best.2 {
            best = (x, residual, tolerance);
        }
        if force > 0.0 && residual.abs() <= 0.0625*tolerance { return Ok((x, iteration)); }
        if residual < 0.0 { lo = x; } else { hi = x; }
    }
    // As in the shared reaction solver, exhausted floating-point refinement
    // may still meet the caller's original tolerance. The final actual-state
    // gate remains mandatory, including after translating the modal states.
    if best.1.abs() <= best.2 { return Ok((best.0, iterations)); }
    Err(ModalCouplingError::ContactSolve {
        residual_n: best.1, tolerance_n: best.2, iterations,
    })
}
