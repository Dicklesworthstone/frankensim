//! Signed physical response declarations for the native affine load family.
use super::*;
use fs_topopt::sdf3::response::limits::ResponseLimitSense3;
use fs_topopt::sdf3::robust::envelope::affine_load_space::response::AffineResponseOptions3;
use std::collections::BTreeSet;

pub(crate) const MAX_LIMITS: usize = 64;

#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) enum Observation {
    /// Mean of one displacement component over a fixed native cut-volume band.
    MeanDisplacement { axis: usize, x_fraction: [f64; 2] },
    /// Native numerical traction integrated over the declared embedded support.
    SupportForce { axis: usize },
    /// The same traction tested by axis cross (position - about_m).
    SupportMoment { axis: usize, about_m: [f64; 3] },
}

impl Observation {
    pub(crate) fn kind(self) -> &'static str {
        match self {
            Self::MeanDisplacement { .. } => "mean-displacement",
            Self::SupportForce { .. } => "support-force",
            Self::SupportMoment { .. } => "support-moment",
        }
    }

    pub(crate) fn unit(self) -> &'static str {
        match self {
            Self::MeanDisplacement { .. } => "m",
            Self::SupportForce { .. } => "N",
            Self::SupportMoment { .. } => "N*m",
        }
    }

    pub(crate) fn axis(self) -> usize {
        match self {
            Self::MeanDisplacement { axis, .. }
            | Self::SupportForce { axis }
            | Self::SupportMoment { axis, .. } => axis,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Limit {
    pub name: String,
    pub observation: Observation,
    pub sense: ResponseLimitSense3,
    pub bound: f64,
    pub scale: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub(crate) struct Family {
    pub options: AffineResponseOptions3,
    pub limits: Vec<Limit>,
}

impl Family {
    pub(crate) fn validate(&self, spec: &Spec) -> Result<()> {
        if self.options.max_limits != MAX_LIMITS
            || self.limits.is_empty()
            || self.limits.len() > MAX_LIMITS
            || !(1e-6..=1e6).contains(&self.options.sharpness)
        {
            return Err(invalid(
                "affine-responses requires 1..=64 signed limits and sharpness in [1e-6,1e6]",
            ));
        }
        let mut names = BTreeSet::new();
        for limit in &self.limits {
            if limit.name.is_empty()
                || limit.name.len() > 64
                || !limit
                    .name
                    .bytes()
                    .all(|c| c.is_ascii_alphanumeric() || c == b'_' || c == b'-')
                || !names.insert(&limit.name)
            {
                return Err(invalid(
                    "response names must be unique 1..=64 byte ASCII letters, digits, underscores or hyphens",
                ));
            }
            if limit.observation.axis() > 2
                || !limit.bound.is_finite()
                || limit.bound.abs() > 1e15
                || !(1e-12..=1e15).contains(&limit.scale)
            {
                return Err(invalid(
                    "response bounds must be finite with magnitude <=1e15; scales must be in [1e-12,1e15] in the response's physical units",
                ));
            }
            match limit.observation {
                Observation::MeanDisplacement {
                    x_fraction: [a, b], ..
                } => {
                    let cells = f64::from(1_u32 << spec.level);
                    if !a.is_finite()
                        || !b.is_finite()
                        || !(0.0 <= a && a < b && b <= 1.0)
                        || (a * cells).fract() != 0.0
                        || (b * cells).fract() != 0.0
                    {
                        return Err(invalid(
                            "mean-displacement x-fraction must satisfy 0<=a<b<=1 and align to the initial octree planes",
                        ));
                    }
                }
                Observation::SupportForce { .. } | Observation::SupportMoment { .. }
                    if spec.fixed.embedded_penalty().is_none() =>
                {
                    return Err(invalid(
                        "support-force and support-moment require an authored embedded support",
                    ));
                }
                Observation::SupportMoment { about_m, .. }
                    if about_m.iter().any(|x| !x.is_finite() || x.abs() > 1e12) =>
                {
                    return Err(invalid(
                        "support-moment about-m requires finite reference coordinates with magnitude <=1e12 m",
                    ));
                }
                _ => {}
            }
        }
        Ok(())
    }

    /// Two vectors per limit at the common 50,000-DOF admission ceiling.
    /// This is a structural workspace envelope, not measured resident memory.
    pub(crate) fn memory_envelope(&self) -> usize {
        self.limits.len().saturating_mul(16 * 50_000)
    }
}

fn axis(node: &Node) -> Result<usize> {
    match &node.kind {
        NodeKind::Symbol(s) | NodeKind::Str(s) => match s.as_str() {
            "x" => Ok(0),
            "y" => Ok(1),
            "z" => Ok(2),
            _ => Err(invalid("response axis must be x, y or z")),
        },
        _ => Err(invalid("response axis must be x, y or z")),
    }
}

fn sense(node: &Node) -> Result<ResponseLimitSense3> {
    match &node.kind {
        NodeKind::Symbol(s) | NodeKind::Str(s) if s == "at-most" => Ok(ResponseLimitSense3::AtMost),
        NodeKind::Symbol(s) | NodeKind::Str(s) if s == "at-least" => {
            Ok(ResponseLimitSense3::AtLeast)
        }
        _ => Err(invalid("response sense must be at-most or at-least")),
    }
}

pub(super) fn parse(node: &Node) -> Result<(Family, f64)> {
    let values = fields(
        node,
        "affine-responses",
        &["sharpness", "joint-sharpness", "limits"],
    )?;
    let rows = super::super::super::list(values[2], "affine response limits")?;
    if rows.is_empty() || rows.len() > MAX_LIMITS {
        return Err(invalid(
            "affine-responses requires 1..=64 original signed inequalities",
        ));
    }
    let mut limits = Vec::with_capacity(rows.len());
    for row in rows {
        let items = super::super::super::list(row, "affine response limit")?;
        let kind = match items.first().map(|v| &v.kind) {
            Some(NodeKind::Symbol(name)) => name.as_str(),
            _ => return Err(invalid("expected a named physical response observation")),
        };
        let (fields, observation, sense_at) = match kind {
            "mean-displacement" => {
                let fields = fields(
                    row,
                    kind,
                    &["name", "axis", "x-fraction", "sense", "bound-m", "scale-m"],
                )?;
                let observation = Observation::MeanDisplacement {
                    axis: axis(fields[1])?,
                    x_fraction: super::super::super::pair(fields[2], "response x-fraction")?,
                };
                (fields, observation, 3)
            }
            "support-force" => {
                let fields = fields(row, kind, &["name", "axis", "sense", "bound-n", "scale-n"])?;
                let observation = Observation::SupportForce {
                    axis: axis(fields[1])?,
                };
                (fields, observation, 2)
            }
            "support-moment" => {
                let fields = fields(
                    row,
                    kind,
                    &["name", "axis", "about-m", "sense", "bound-n-m", "scale-n-m"],
                )?;
                let observation = Observation::SupportMoment {
                    axis: axis(fields[1])?,
                    about_m: vector(fields[2])?,
                };
                (fields, observation, 3)
            }
            _ => {
                return Err(invalid(
                    "affine responses must be mean-displacement, support-force or support-moment",
                ));
            }
        };
        let name = match &fields[0].kind {
            NodeKind::Str(name) => name.clone(),
            _ => return Err(invalid("each response requires a quoted name")),
        };
        limits.push(Limit {
            name,
            observation,
            sense: sense(fields[sense_at])?,
            bound: scalar(fields[sense_at + 1])?,
            scale: scalar(fields[sense_at + 2])?,
        });
    }
    Ok((
        Family {
            options: AffineResponseOptions3 {
                max_limits: MAX_LIMITS,
                sharpness: scalar(values[0])?,
            },
            limits,
        },
        scalar(values[1])?,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    const SOURCE: &str = include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../examples/marquee/bracket-3d-affine-responses.fsim"
    ));

    #[test]
    fn original_signed_si_limits_and_their_joint_policy_bind_source_identity() {
        let spec = super::super::parse(SOURCE).unwrap();
        let family = spec.affine_responses.as_ref().unwrap();
        assert_eq!(family.limits.len(), 5);
        assert_eq!(family.limits[0].sense, ResponseLimitSense3::AtLeast);
        assert_eq!(family.limits[0].bound, -0.02 / 1024.0);
        assert_eq!(family.limits[0].observation.unit(), "m");
        assert_eq!(family.limits[2].observation.unit(), "N");
        assert_eq!(family.limits[4].observation.unit(), "N*m");
        assert_eq!(
            spec.affine_load_space.unwrap().response_options,
            Some(family.options)
        );
        assert_eq!(super::super::parse(&spec.canonical).unwrap().id, spec.id);
        for source in [
            SOURCE.replace(":bound-m -0.00001953125", ":bound-m 0.0"),
            SOURCE.replace(":about-m (0.0 0.0 0.0)", ":about-m (0.5 0.0 0.0)"),
        ] {
            assert_ne!(super::super::parse(&source).unwrap().id, spec.id);
        }
    }

    #[test]
    fn response_admission_refuses_wrong_units_supports_patches_and_ambiguous_policy() {
        for (from, to) in [
            ("tip-z-upper", "tip-z-lower"),
            ("tip-z-lower", "not a scalar identifier"),
            (":axis z", ":axis diagonal"),
            (":x-fraction (0.5 1.0)", ":x-fraction (0.3 1.0)"),
            (":bound-m -0.00001953125", ":bound-n -0.00001953125"),
            (":scale-m 0.00001953125", ":scale-m 0.0"),
            (":sense at-least", ":sense absolute"),
            (":about-m (0.0 0.0 0.0)", ":about-m (0.0 0.0)"),
            (":sharpness 40.0", ":sharpness 0.0"),
            (
                "(embedded :axis x :fraction (0.0 0.5) :penalty 32.0)",
                "left",
            ),
        ] {
            assert!(
                super::super::parse(&SOURCE.replace(from, to)).is_err(),
                "admitted {to}"
            );
        }
        let mismatched = SOURCE.replace(
            "(affine-responses\n    :sharpness 40.0\n    :joint-sharpness 40.0",
            "(affine-responses\n    :sharpness 40.0\n    :joint-sharpness 50.0",
        );
        assert_ne!(mismatched, SOURCE);
        assert!(super::super::parse(&mismatched).is_err());

        let mut spec = super::super::parse(SOURCE).unwrap();
        let seed = spec.affine_responses.as_ref().unwrap().limits[0].clone();
        for count in [64, 65] {
            let family = spec.affine_responses.as_mut().unwrap();
            family.limits = (0..count)
                .map(|i| {
                    let mut limit = seed.clone();
                    limit.name = format!("signed-{i}");
                    limit
                })
                .collect();
            assert_eq!(family.memory_envelope(), count * 16 * 50_000);
            assert_eq!(
                spec.affine_responses
                    .as_ref()
                    .unwrap()
                    .validate(&spec)
                    .is_ok(),
                count == 64
            );
        }
        spec.affine_responses.as_mut().unwrap().limits.clear();
        assert!(
            spec.affine_responses
                .as_ref()
                .unwrap()
                .validate(&spec)
                .is_err()
        );
    }
}
