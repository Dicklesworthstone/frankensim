//! Passive rank-one factors lowered into the existing complete modal owner.
use super::*;
use crate::render::schedule::force::coupled::ModalColumnConnection;

pub(crate) fn condensed_columns(
    models: &[ModalAcousticTimeModel],
    links: Vec<ModalConnection>,
    pairs: &[Option<BridgePair>],
    config: ModalCouplingConfig,
    gate: &CancelGate,
) -> Result<Vec<ModalColumnConnection>, ReciprocalStringError> {
    let mut offsets = Vec::with_capacity(models.len());
    let mut count = 0usize;
    for model in models {
        poll(gate)?;
        offsets.push(count);
        count = count
            .checked_add(model.modes().len())
            .filter(|n| *n <= config.max_modes)
            .ok_or_else(|| invalid("coupled residual bridges exceed the modal budget"))?;
    }
    let mut expanded = links.len();
    for pair in pairs.iter().flatten().filter(|pair| pair.law.coupled) {
        // One extra pure damper only when its row differs from the first
        // spring. A zero first spring leaves no spurious zero-law column.
        let law = pair.law;
        let separate = law.back != 0.0 && law.first.damping_n_s_m > 0.0;
        if separate {
            expanded = expanded
                .checked_add(usize::from(law.first.stiffness_n_m > 0.0))
                .ok_or_else(|| invalid("coupled residual connection count overflow"))?;
        }
    }
    let extent = expanded
        .checked_add(1)
        .ok_or_else(|| invalid("coupled residual setup overflow"))?;
    let terms = count
        .checked_mul(extent)
        .and_then(|n| n.checked_mul(extent))
        .ok_or_else(|| invalid("coupled residual setup overflow"))?;
    if expanded > config.max_connections || terms > config.max_setup_terms {
        return Err(invalid(
            "coupled residual bridge factors exceed the declared native connection/setup budgets",
        ));
    }
    let mut result = Vec::with_capacity(expanded);
    let mut consumed = vec![false; links.len()];
    for (index, link) in links.iter().enumerate() {
        poll(gate)?;
        if consumed[index] {
            continue;
        }
        let first = column(link, models, &offsets, count)?;
        let pair = pairs
            .get(index)
            .and_then(Option::as_ref)
            .filter(|p| p.law.coupled);
        let Some(pair) = pair else {
            result.push(first);
            continue;
        };
        let second_index = pairs.len() + pair.directed;
        let second_link = links
            .get(second_index)
            .ok_or_else(|| invalid("coupled residual bridge lost its original support"))?;
        let mut second = column(second_link, models, &offsets, count)?;
        let law = pair.law;
        for (b, a) in second.column.iter_mut().zip(&first.column) {
            *b = checked(*b - law.mix * a)?;
        }
        second.stiffness_n_m = law.second.stiffness_n_m;
        second.damping_n_s_m = law.second.damping_n_s_m;
        second.rest_extension_m =
            checked(second.rest_extension_m - law.mix * first.rest_extension_m)?;
        let separate = law.back != 0.0 && law.first.damping_n_s_m > 0.0;
        let mut elastic = first.clone();
        elastic.stiffness_n_m = law.first.stiffness_n_m;
        elastic.damping_n_s_m = if separate {
            0.0
        } else {
            law.first.damping_n_s_m
        };
        if elastic.stiffness_n_m > 0.0 || elastic.damping_n_s_m > 0.0 {
            result.push(elastic);
        }
        if separate {
            let mut damping = first;
            for (a, b) in damping.column.iter_mut().zip(&second.column) {
                *a = checked(*a - law.back * b)?;
            }
            damping.stiffness_n_m = 0.0;
            damping.damping_n_s_m = law.first.damping_n_s_m;
            damping.rest_extension_m = 0.0;
            result.push(damping);
        }
        result.push(second);
        consumed[second_index] = true;
    }
    Ok(result)
}

fn column(
    link: &ModalConnection,
    models: &[ModalAcousticTimeModel],
    offsets: &[usize],
    count: usize,
) -> Result<ModalColumnConnection, ReciprocalStringError> {
    let mut column = vec![0.0; count];
    for (attachment, sign) in [(&link.left, 1.0), (&link.right, -1.0)] {
        let model = models
            .get(attachment.component)
            .ok_or_else(|| invalid("residual bridge attachment names an unknown component"))?;
        if attachment.shapes.len() != model.modes().len() {
            return Err(invalid(
                "residual bridge attachment differs from its original modal basis",
            ));
        }
        for (i, shape) in attachment.shapes.iter().enumerate() {
            let slot = &mut column[offsets[attachment.component] + i];
            *slot = checked(*slot + sign * shape)?;
        }
    }
    Ok(ModalColumnConnection {
        column,
        stiffness_n_m: link.stiffness_n_m,
        damping_n_s_m: link.damping_n_s_m,
        rest_extension_m: link.rest_extension_m,
    })
}
