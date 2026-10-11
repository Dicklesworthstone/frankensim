//! Bind natural walls to the same full solid/air tangent as forced convection.
//! Re-evaluate the physical card on the unchanged field before the core owner
//! checks its complete residual. Independent ambient controls retain BOTH
//! film-temperature partials; their wall feedback is included only once.

use std::collections::BTreeMap;
use fs_conduction::{ScalarField, ThermalBc};
use fs_conduction::adjoint::RobinGradient;
use fs_project::ThermalBoundaryCondition;
use super::{Cx, ProjectSpec, RungSolved, SolveRefusal, bad, finite,
    natural_feedback, poll, zeros};
use super::super::natural;

pub(super) struct Bound {
    pub replacements: Vec<(usize, f64, f64)>,
    pub wall_slopes: Vec<f64>,
    ambient_log_slopes: Vec<Option<f64>>,
}

impl Bound {
    pub fn ambient(
        &self, names: &[&str], gradient: &RobinGradient,
    ) -> Result<BTreeMap<String, f64>, SolveRefusal> {
        if gradient.references.len() != names.len() || gradient.log_htc.len() != names.len()
            || self.ambient_log_slopes.len() != names.len()
        {
            return Err(bad("mixed adjoint lost the natural control ordering"));
        }
        let mut result = BTreeMap::new();
        for (i, &name) in names.iter().enumerate() {
            if let Some(slope) = self.ambient_log_slopes[i] {
                result.insert(name.to_string(), finite(gradient.references[i]
                    + gradient.log_htc[i] * slope)?);
            }
        }
        Ok(result)
    }
}

pub(super) fn bind(
    cx: &Cx<'_>, spec: &ProjectSpec, solved: &RungSolved,
    names: &[&str], regions: &[usize], means: &[f64],
) -> Result<Bound, SolveRefusal> {
    poll(cx)?;
    let m = names.len();
    if m > 64 || regions.len() != m || means.len() != m {
        return Err(bad("mixed adjoint needs matching bounded natural/air trace arrays"));
    }
    let data = solved.adjoint_data.as_ref().ok_or_else(|| bad("mixed adjoint lost its operator"))?;
    let setup = spec.cooling.as_ref().and_then(|c| c.conduction.as_ref())
        .ok_or_else(|| bad("mixed adjoint lost its conduction declaration"))?;
    let mut bound = Bound {
        replacements: Vec::new(), wall_slopes: zeros(m)?, ambient_log_slopes: vec![None; m],
    };
    let laws = natural::natural_laws(setup)?;
    for law in &laws {
        poll(cx)?;
        let i = names.iter().position(|&name| name == law.target)
            .ok_or_else(|| bad("mixed natural law has no retained Robin trace"))?;
        if data.air_paths.iter().flat_map(|path| path.segments())
            .any(|segment| segment.region() == law.target)
        {
            return Err(bad("natural and forced convection must occupy distinct traces"));
        }
        let declared = setup.boundaries.iter().find(|row| row.target == law.target)
            .ok_or_else(|| bad("mixed natural law has no physical declaration"))?;
        let ThermalBoundaryCondition::NaturalConvection {
            characteristic_length, ambient_temperature, correlation,
        } = &declared.condition else { return Err(bad("mixed natural law changed its declaration")); };
        let ThermalBc::Robin { htc: ScalarField::Uniform(_), t_ref: ScalarField::Uniform(ambient) }
            = &data.boundary.conditions()[regions[i]]
        else { return Err(bad("mixed natural adjoint needs its original uniform Robin row")); };
        if correlation != natural_feedback::CARD || *ambient != ambient_temperature.value {
            return Err(bad("mixed natural adjoint requires the original Churchill-Chu card and ambient"));
        }
        let pressure = spec.envelope.as_ref()
            .ok_or_else(|| bad("mixed natural adjoint has no pressure envelope"))?.pressure.value;
        let h = natural::coefficient(law, means[i] - *ambient, pressure)?.htc_w_m2_k;
        let [wall, reservoir] = natural_feedback::partials(
            characteristic_length.value, means[i], *ambient, h)?;
        bound.replacements.push((regions[i], h, *ambient));
        bound.wall_slopes[i] = wall;
        bound.ambient_log_slopes[i] = Some(finite(reservoir / h)?);
    }
    poll(cx)?;
    Ok(bound)
}
