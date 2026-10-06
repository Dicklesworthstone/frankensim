//! Named contact controls on the SAME accepted native field and complete dual.
use std::collections::BTreeMap;
use super::{Cx, ProjectSpec, RungSolved, SolveRefusal, bad, json_string, lower, poll, row};

pub(super) const OUTPUT: &str = "temperature-max-contact-adjoint";
pub(super) const SCOPE: &str = "This contact-extended request additionally varies each named matching interface's resistance via R''_face -> s R''_face at s=1. Every retained face resistance and both independent P1 temperature traces are used with the same complete thermal dual. A resistance map keeps its spatial pattern; different interfaces remain independent even when they share a card. Contact geometry, pairing, pressure, finish, card selection and temperature dependence stay fixed. Contact rows are local normalized resistance controls, not absolute resistance derivatives, card uncertainty, a new material claim or a nonmatching-contact capability.";

pub(super) fn requested(spec: &ProjectSpec) -> bool {
    spec.outputs.as_deref().unwrap_or(&[]).iter().any(|r| r.name == OUTPUT)
}

pub(super) fn rows(
    cx: &Cx<'_>, spec: &ProjectSpec, solved: &RungSolved, lambda: &[f64], max_faces: usize,
) -> Result<Vec<String>, SolveRefusal> {
    poll(cx)?;
    let declared = spec.interface_cards.as_deref().unwrap_or(&[]);
    let data = solved.adjoint_data.as_ref().ok_or_else(|| bad("contact derivative lost the accepted operator"))?;
    let Some(interfaces) = data.interfaces.as_ref() else {
        if declared.is_empty() { return Ok(Vec::new()); }
        return Err(bad("declared contact controls have no retained contact operator"));
    };
    if interfaces.surface_count() != declared.len() || lambda.len() != solved.mesh.vertex_count() {
        return Err(bad("contact controls do not match the complete native interface/dual binding"));
    }
    for &(v, prescribed) in data.boundary.dirichlet() {
        poll(cx)?;
        if lambda[v] != 0.0 || solved.solution.temperature[v] != prescribed {
            return Err(bad("contact controls require unchanged prescribed temperatures and zero fixed-node dual entries"));
        }
    }
    let gradients = interfaces.matching_resistance_scale_pullback(cx,
        &solved.solution.temperature, lambda, max_faces).map_err(lower)?;
    let mut by_name = BTreeMap::new();
    for gradient in gradients {
        poll(cx)?;
        if by_name.insert(gradient.interface.clone(), gradient).is_some() {
            return Err(bad("contact gradient names are not unique"));
        }
    }
    let mut result = Vec::with_capacity(declared.len());
    for (ordinal, binding) in declared.iter().enumerate() {
        poll(cx)?;
        let gradient = by_name.remove(&binding.interface)
            .ok_or_else(|| bad("missing or repeated declared contact-control identity"))?;
        let card = gradient.card_identity.to_hex();
        if binding.card != card || gradient.face_pairs == 0 {
            return Err(bad("contact gradient differs from its retained material card or has no bound face pairs"));
        }
        let mut entry = row("contact-resistance-multiplier", &binding.interface, ordinal,
            "1", gradient.derivative)?;
        entry.pop();
        entry.push_str(&format!(",\"reference_value\":1,\"interface_card\":{},\"mapped\":{},\"face_pairs\":{},\"parameterization\":\"R''_face -> s R''_face, evaluated at s=1\"}}",
            json_string(&card), gradient.mapped, gradient.face_pairs));
        result.push(entry);
    }
    if !by_name.is_empty() { return Err(bad("a native contact gradient has no declared owner")); }
    poll(cx)?;
    Ok(result)
}
