//! The existing HEX8 projection, retained unchanged by every memory operator.

use super::{
    Element, Hex8ElasticProblem, IntegrationPoint, TetViscoField, ViscoElementError, checkpoint,
};
use crate::linear::Formulation;
use crate::linear3::hex8::tensor_rule;
use crate::linear3::isoparametric::StrainProjection;
use fs_exec::Cx;

pub(super) fn prepare(
    prepared: &mut Vec<Element>,
    problem: &Hex8ElasticProblem<'_>,
    field: TetViscoField<'_>,
    cx: &Cx<'_>,
) -> Result<(), ViscoElementError> {
    let rule = tensor_rule(problem.stiffness_points_per_dir);
    for (element, conn) in problem.hexahedra.iter().enumerate() {
        checkpoint(cx)?;
        let nodes = conn.map(|node| problem.nodes_m[node]);
        let projection = StrainProjection::prepare(
            nodes,
            element,
            problem.formulation,
            &rule,
            |nodes, reference, element, point| problem.point(nodes, reference, element, point),
            cx,
        )?;
        let mut points = Vec::with_capacity(rule.len());
        for (index, &(reference, weight)) in rule.iter().enumerate() {
            checkpoint(cx)?;
            let compatible_point;
            let point = if let Some(prepared) = projection.point(index) {
                prepared
            } else {
                compatible_point = problem.point(nodes, reference, element, index)?;
                &compatible_point
            };
            let b = projection.columns::<24>(point);
            points.push(IntegrationPoint {
                b: (0..24).map(|j| core::array::from_fn(|i| b[i][j])).collect(),
                compatible_gradients: (problem.formulation == Formulation::BBar)
                    .then(|| point.gradients.to_vec()),
                volume_m3: weight * point.determinant,
            });
        }
        prepared.push(Element {
            dofs: conn
                .iter()
                .flat_map(|node| [3 * node, 3 * node + 1, 3 * node + 2])
                .collect(),
            law: field.at(element).law.clone(),
            density_kg_m3: field.at(element).density_kg_m3,
            material_state_identity: field.at(element).material_state_identity,
            points,
        });
    }
    Ok(())
}
