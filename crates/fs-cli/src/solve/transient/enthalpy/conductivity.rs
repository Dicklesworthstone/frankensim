//! Regional phase conductivity lowers to the shared enthalpy transport owner.

use fs_project::{TransientEnthalpy, TransientRegionEnthalpy};

use super::*;

pub(super) fn bind<'m, 'c>(
    cx: &Cx<'_>,
    engine: HeterogeneousEnthalpyBackwardEuler<'m, 'c>,
    storage: &TransientEnthalpy,
    assignment: &[usize],
    deadline: Option<(Instant, f64)>,
) -> Result<HeterogeneousEnthalpyBackwardEuler<'m, 'c>, SolveRefusal> {
    if !storage.materials.iter().any(|row| row.phase_conductivity.is_some()) {
        return Ok(engine);
    }
    let identity = fs_conduction::material::LiquidMassFractionConductivity::declared(1.0, 1.0)
        .map_err(lower)?;
    let mut models = Vec::with_capacity(storage.materials.len());
    for row in &storage.materials {
        poll(cx, deadline)?;
        models.push(match &row.phase_conductivity {
            Some(declaration) => declaration.model()
                .map_err(|v| conduction_error(v.code, v.what, v.fix))?,
            None => identity,
        });
    }
    let mut element_laws = Vec::with_capacity(assignment.len());
    for (element, &material) in assignment.iter().enumerate() {
        if element % 256 == 0 {
            poll(cx, deadline)?;
        }
        element_laws.push(*models.get(material)
            .ok_or_else(|| bad("phase conductivity has an unknown reference-material assignment"))?);
    }
    engine.with_phase_conductivity(cx, element_laws).map_err(failed)
}

pub(super) fn receipt(row: &TransientRegionEnthalpy) -> Result<String, SolveRefusal> {
    let Some(declaration) = &row.phase_conductivity else {
        return Ok(String::new());
    };
    let model = declaration.model()
        .map_err(|v| conduction_error(v.code, v.what, v.fix))?;
    Ok(format!(
        ",\"phase_conductivity\":{{\"law\":\"linear-liquid-mass-fraction\",\"solid_multiplier\":{},\"liquid_multiplier\":{},\"source\":{},\"coordinate\":\"element-mean-nodal-liquid-mass-fraction\",\"authority\":\"Declared\"}}",
        number(model.solid_multiplier())?,
        number(model.liquid_multiplier())?,
        json_string(&declaration.source),
    ))
}
