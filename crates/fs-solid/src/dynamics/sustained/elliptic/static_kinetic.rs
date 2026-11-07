//! Coupled static capture and kinetic breakaway with distinct ellipse shapes.

use super::*;

#[allow(clippy::too_many_arguments)]
pub(crate) fn advance<B: Native, P, E>(
    body: &mut B,
    contacts: &[SolidSustainedEllipticContact<'_>],
    kinetic_coefficients_mu: &[[f64; 2]],
    initial_forces_n: &[[f64; 3]],
    external_force_n: &[f64],
    controls: SolidSustainedControls,
    cx: &Cx<'_>,
    admit_path: impl FnOnce(&SolidSustainedEllipticPath<'_, '_, B::Owner>) -> Result<P, E>,
) -> Result<SolidSustainedEllipticStep<P, B::Step>, SolidSustainedError<E, B::Error>> {
    poll(cx)?;
    let m = contacts.len();
    let n = body.state().displacement.len();
    if m == 0
        || m > super::super::MAX_CONTACTS
        || m > controls.solve.maximum_contacts
        || kinetic_coefficients_mu.len() != m
        || initial_forces_n.len() != m
    {
        return Err(SolidSustainedError::Invalid(
            "one bounded complete kinetic pair and warm start per static ellipse required",
        ));
    }
    let bytes = m as u128
        * (core::mem::size_of::<SolidSustainedEllipticContact<'_>>()
            + core::mem::size_of::<usize>()
            + 2 * core::mem::size_of::<[f64; 3]>()
            + core::mem::size_of::<[f64; 2]>()
            + 2 * core::mem::size_of::<bool>()) as u128
        + 256;
    let selection_terms = 8 * n as u128 * m as u128 + 128 * m as u128;
    bound(
        "directional friction regime scratch bytes",
        bytes,
        controls.maximum_workspace_bytes as u128,
    )?;
    bound(
        "directional friction regime admission work",
        selection_terms,
        controls.maximum_response_multiply_terms as u128,
    )?;
    let mut remaining = controls;
    remaining.maximum_workspace_bytes -= bytes as usize;
    remaining.maximum_response_multiply_terms -= selection_terms as u64;
    let mut order: Vec<_> = (0..m).collect();
    order.sort_unstable_by_key(|&index| contacts[index].id);
    if order
        .windows(2)
        .any(|w| contacts[w[0]].id == contacts[w[1]].id)
    {
        return Err(SolidSustainedError::Invalid("duplicate contact identity"));
    }
    let mut selected: Vec<_> = order.iter().map(|&i| contacts[i]).collect();
    let kinetic: Vec<_> = order.iter().map(|&i| kinetic_coefficients_mu[i]).collect();
    let mut initial: Vec<_> = order.iter().map(|&i| initial_forces_n[i]).collect();
    let mut can_switch = vec![false; m];
    for (i, contact) in selected.iter_mut().enumerate() {
        poll(cx)?;
        let static_mu = contact.friction.coefficients_mu;
        if (0..2).any(|axis| {
            !kinetic[i][axis].is_finite()
                || kinetic[i][axis] < 0.0
                || kinetic[i][axis] > static_mu[axis]
        }) || contact
            .velocity_weights
            .iter()
            .any(|row| row.len() != n || row.iter().any(|value| !value.is_finite()))
            || contact
                .prescribed_surface_velocity_mps
                .iter()
                .any(|value| !value.is_finite())
        {
            return Err(SolidSustainedError::Invalid(
                "complete finite traces and kinetic principal strengths within each static bound required",
            ));
        }
        // The same original axis is used for both shapes. Both must be admitted
        // even if the initial velocity would not exercise one of them.
        EllipticForceCoordinates::prepare(contact.friction).map_err(law_error)?;
        let kinetic_friction = EllipticFriction {
            coefficients_mu: kinetic[i],
            ..contact.friction
        };
        EllipticForceCoordinates::prepare(kinetic_friction).map_err(law_error)?;
        let speed = tangent_norm(contact, body.state().velocity, 1.0)?;
        if speed > 0.0 {
            contact.friction = kinetic_friction;
        } else {
            can_switch[i] = kinetic[i] != static_mu;
        }
    }

    let terms_per_response = super::super::response_terms(body, m);
    let (_, staging_terms) = staging_requirements(n, m);
    let mut response_evaluations = 0usize;
    let mut feasibility_refinements = 0usize;
    let mut resolves = 0usize;
    let mut admission = Some(admit_path);
    let mut switch = vec![false; m];
    let mut warm_start = vec![[0.0; 3]; m];
    loop {
        poll(cx)?;
        bound(
            "directional friction regime selection work",
            selection_terms,
            remaining.maximum_response_multiply_terms as u128,
        )?;
        let mut trial_controls = remaining;
        trial_controls.maximum_response_multiply_terms -= selection_terms as u64;
        let mut retry = false;
        let mut trial_evaluations = 0usize;
        let result = advance_selected(
            body,
            &selected,
            &initial,
            external_force_n,
            trial_controls,
            cx,
            |path, solve, coordinates| {
                trial_evaluations = solve.response_evaluations;
                response_evaluations += trial_evaluations;
                feasibility_refinements += solve.feasibility_refinements;
                for (i, contact) in selected.iter().enumerate() {
                    poll(cx)?;
                    if can_switch[i]
                        && (tangent_norm(contact, path.velocity_after_mps, 1.0)?
                            > controls.velocity_tolerance_mps
                            || tangent_norm(contact, path.displacement_change_m, path.duration_s)?
                                > controls.displacement_tolerance_m)
                    {
                        switch[i] = true;
                        retry = true;
                    }
                    warm_start[i] = coordinates[i]
                        .force_to_physical(solve.forces_n[i])
                        .map_err(law_error)?;
                }
                if retry {
                    // This private control return precedes the native audit and
                    // commit. The flag is inspected only with this Path error;
                    // no physical or numerical refusal is treated as a retry.
                    Err(SolidSustainedError::Invalid("directional breakaway trial"))
                } else {
                    Ok(())
                }
            },
            |path| {
                admission
                    .take()
                    .expect("only the final directional regime admits a physical path")(
                    path
                )
                .map_err(SolidSustainedError::Path)
            },
        );
        match result {
            Err(SolidSustainedError::Path(_)) if retry => {
                // Every discarded solve pays for its coordinate staging and
                // reserves its original physical audit. No retry replenishes
                // the native product or shared solver response budget.
                let spent = selection_terms
                    + staging_terms
                    + (trial_evaluations as u128 + 1) * terms_per_response;
                bound(
                    "directional friction regime response work",
                    spent,
                    remaining.maximum_response_multiply_terms as u128,
                )?;
                remaining.maximum_response_multiply_terms -= spent as u64;
                remaining.solve.maximum_response_evaluations = remaining
                    .solve
                    .maximum_response_evaluations
                    .saturating_sub(trial_evaluations);
                if remaining.solve.maximum_response_evaluations < 2 {
                    return Err(SolidSustainedError::Budget {
                        resource: "directional friction regime response evaluations",
                        required: response_evaluations as u128 + 2,
                        maximum: controls.solve.maximum_response_evaluations as u128,
                    });
                }
                for i in 0..m {
                    if switch[i] {
                        selected[i].friction.coefficients_mu = kinetic[i];
                        can_switch[i] = false;
                        switch[i] = false;
                    }
                }
                initial.copy_from_slice(&warm_start);
                resolves += 1;
            }
            Ok(mut step) => {
                // Pure report assembly only after original native acceptance.
                step.transformed.friction_regime_resolves = resolves;
                step.transformed.solve.response_evaluations = response_evaluations;
                step.transformed.solve.feasibility_refinements = feasibility_refinements;
                return Ok(step);
            }
            Err(error) => return Err(flatten(error)),
        }
    }
}

fn tangent_norm<E, N>(
    contact: &SolidSustainedEllipticContact<'_>,
    native: &[f64],
    surface_duration: f64,
) -> Result<f64, SolidSustainedError<E, N>> {
    let mut relative = [0.0; 2];
    for a in 0..2 {
        relative[a] = finite(
            dot(contact.velocity_weights[a + 1], native)
                - Dd::from_f64(surface_duration)
                    * Dd::from_f64(contact.prescribed_surface_velocity_mps[a + 1]),
        )?;
    }
    let norm = relative[0].hypot(relative[1]);
    if norm.is_finite() {
        Ok(norm)
    } else {
        Err(SolidSustainedError::Invalid(
            "unrepresentable original tangent motion",
        ))
    }
}
