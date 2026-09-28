//! Explicit native constraint selection; volume-only studies do not invent a
//! stress allowance or emit unmeasured stress evidence.
use super::*;

#[derive(Debug, Clone)]
pub(crate) enum ProjectedControls {
    Stress(projected::Controls),
    Volume(projected::volume::Controls),
}

impl ProjectedControls {
    pub(crate) fn canonical(&self, out: &mut String) {
        match self {
            Self::Stress(policy) => policy.canonical(out),
            Self::Volume(policy) => policy.canonical(out),
        }
    }
}

pub(crate) fn parse_controls(fields: &[Node], target: f64) -> Result<Option<ProjectedControls>> {
    let volume = fields.windows(2).any(|pair|
        matches!(&pair[0].kind, NodeKind::Keyword(key) if key == "constraint-mode")
        && matches!(&pair[1].kind, NodeKind::Symbol(value) if value == "projected-volume"));
    if volume {
        projected::volume::Controls::parse(fields, target)
            .map(|policy| Some(ProjectedControls::Volume(policy)))
    } else {
        projected::parse_controls(fields, target)
            .map(|policy| policy.map(ProjectedControls::Stress))
    }
}

pub(super) fn stress_controls(spec: &ElasticitySpec) -> Result<&projected::Controls> {
    match &spec.projected {
        Some(ProjectedControls::Stress(policy)) => Ok(policy),
        _ => Err(malformed("this producer requires an explicit projected-stress policy")),
    }
}
