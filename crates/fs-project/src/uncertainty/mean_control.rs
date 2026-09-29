//! Deterministic pre-sampling mean controls. Secants retain complete native
//! endpoint solves; nominal adjoints use one native solve at marginal means.

use super::{Node, UniformParameter, Result, error, fields, integer, list, symbol};

/// Chosen before any probability observations; never selected for a favorable result.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MeanControlMethod {
    CoordinateSecant,
    NominalAdjoint,
}
impl MeanControlMethod {
    #[must_use]
    pub const fn tag(self) -> &'static str {
        match self { Self::CoordinateSecant => "coordinate-secant", Self::NominalAdjoint => "nominal-adjoint" }
    }
}

/// Explicit additional native-solve allowance for frozen mean coefficients.
/// The whole study wall budget includes calibration. A failed probe refuses the
/// requested control; it never switches method or enters the probability sample.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MeanControlPolicy {
    /// Original lifetime cap for completed or terminally refused probes, 0..=64.
    pub max_solves: usize,
    pub method: MeanControlMethod,
}
impl MeanControlPolicy {
    /// No calibration is necessary when every centered input is identically zero.
    #[must_use]
    pub fn probe_count(self, parameters: &[UniformParameter]) -> usize {
        match self.method {
            MeanControlMethod::CoordinateSecant => probe_count(parameters),
            MeanControlMethod::NominalAdjoint => usize::from(parameters.iter().any(|p| p.low != p.high)),
        }
    }
    /// Exact predeclared calibration point, never a random draw.
    pub fn probe(self, parameters: &[UniformParameter], ordinal: usize) -> Result<Vec<f64>> {
        match self.method {
            MeanControlMethod::CoordinateSecant => probe(parameters, ordinal),
            MeanControlMethod::NominalAdjoint if ordinal == 0 && self.probe_count(parameters) == 1 => Ok(means(parameters)),
            MeanControlMethod::NominalAdjoint => Err(error("nominal adjoint probe ordinal exceeds its plan")),
        }
    }
}

pub(super) fn parse(node: &Node, parameters: &[UniformParameter]) -> Result<MeanControlPolicy> {
    let nodes = list(node)?;
    let name = nodes.first().ok_or_else(|| error("empty mean-control policy"))?;
    let method = if symbol(name, "coordinate-secant").is_ok() { MeanControlMethod::CoordinateSecant }
        else { symbol(name, "nominal-adjoint")?; MeanControlMethod::NominalAdjoint };
    let f = fields(&nodes[1..], &["max-solves"])?;
    let max_solves = usize::try_from(integer(f["max-solves"])?).map_err(|_| error("probe cap overflow"))?;
    let policy = MeanControlPolicy { max_solves, method };
    if max_solves > 64 || policy.probe_count(parameters) > max_solves {
        return Err(error("mean-control probes must fit max-solves <=64: two per variable for secants, one total for nominal adjoints"));
    }
    for p in parameters {
        if p.low == p.high { continue; }
        match method {
            MeanControlMethod::CoordinateSecant if !(p.high - p.low).is_finite() =>
                return Err(error("mean-control support width is not finite")),
            MeanControlMethod::NominalAdjoint if matches!(p.target,
                super::Target::AirInletTemperature | super::Target::FanSpeedRatio) =>
                return Err(error("nominal adjoint does not yet supply air-inlet or fan-speed derivatives; explicitly select coordinate-secant for those controls")),
            _ => {},
        }
    }
    Ok(policy)
}

/// Number of predeclared whole-model SECANT solves, never random draws.
#[must_use]
pub fn probe_count(parameters: &[UniformParameter]) -> usize {
    parameters.iter().filter(|p| p.low != p.high).count() * 2
}

/// Analytic physical marginal means. Dependence never changes these means.
#[must_use]
pub fn means(parameters: &[UniformParameter]) -> Vec<f64> {
    parameters.iter().map(|p| p.low.midpoint(p.high)).collect()
}

/// Exact next secant probe: low then high for each nonconstant coordinate.
/// Points may be off a singular copula's support; these are deterministic model
/// probes, not probability observations. The project separately admits each.
pub fn probe(parameters: &[UniformParameter], ordinal: usize) -> Result<Vec<f64>> {
    let variable = ordinal / 2;
    let index = parameters.iter().enumerate().filter(|(_, p)| p.low != p.high)
        .nth(variable).map(|(i, _)| i).ok_or_else(|| error("mean-control probe ordinal exceeds its plan"))?;
    let mut point = means(parameters);
    point[index] = if ordinal % 2 == 0 { parameters[index].low } else { parameters[index].high };
    Ok(point)
}

/// Freeze full-support secant slopes after ALL probes succeed. Never fit to,
/// select using, or include these values in the later random observation set.
/// A nonlinear/nonsmooth model is permitted; no local derivative is claimed.
pub fn coefficients(parameters: &[UniformParameter], values: &[f64]) -> Result<Vec<f64>> {
    if values.len() != probe_count(parameters) || !values.iter().all(|v| v.is_finite()) {
        return Err(error("mean control requires its complete finite probe sequence"));
    }
    let mut ordinal = 0;
    let mut gradient = Vec::with_capacity(parameters.len());
    for parameter in parameters {
        let slope = if parameter.low == parameter.high { 0.0 } else {
            let width = parameter.high - parameter.low;
            if !(width.is_finite() && width > 0.0) {
                return Err(error("mean-control support width is invalid"));
            }
            let value = (values[ordinal + 1] - values[ordinal]) / width;
            ordinal += 2;
            value
        };
        if !slope.is_finite() { return Err(error("mean-control secant is not representable")); }
        gradient.push(slope);
    }
    Ok(gradient)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::uncertainty::{Target, UncertaintyStudy};
    const STUDY: &str = r#"(fsim-uncertainty-study :version 3 :project "p.fsim"
      :samples 4 :seed 7 :wall-time 120s :method monte-carlo :correlation independent
      :qoi "temperature-max" :geometry ((mesh :role "solid" :path "s.stl" :unit "m" :max-hole-edges 0))
      :materials () :interfaces () :mean-control (coordinate-secant :max-solves 4)
      :parameters ((uniform :name "power" :target power :entity "solid" :low 2W :high 6W)
                   (uniform :name "ambient" :target convection-temperature :entity "wall" :low 290K :high 310K)))"#;
    fn params() -> Vec<UniformParameter> { UncertaintyStudy::parse(STUDY).unwrap().parameters().to_vec() }
    #[test]
    fn predeclared_points_and_affine_coefficients_are_physical() {
        let p = params();
        assert_eq!(means(&p), [4.0,300.0]);
        let points: Vec<_> = (0..probe_count(&p)).map(|i| probe(&p,i).unwrap()).collect();
        assert_eq!(points, [vec![2.0,300.0],vec![6.0,300.0],vec![4.0,290.0],vec![4.0,310.0]]);
        let values: Vec<_> = points.iter().map(|x| 5.0 + 3.0*x[0] + 2.0*x[1]).collect();
        assert_eq!(coefficients(&p,&values).unwrap(), [3.0,2.0]);
        assert!(probe(&p,4).is_err());
    }
    #[test]
    fn singleton_inputs_spend_no_probe_and_center_exactly() {
        let mut p=params(); p[0].high=p[0].low;
        assert_eq!(probe_count(&p),2); assert_eq!(probe(&p,0).unwrap(),[2.0,290.0]);
        assert_eq!(coefficients(&p,&[10.0,50.0]).unwrap(),[0.0,2.0]);
        p[1].high=p[1].low;
        assert_eq!(probe_count(&p),0); assert_eq!(coefficients(&p,&[]).unwrap(),[0.0,0.0]);
    }
    #[test]
    fn policy_versions_and_exact_probe_allowances_are_admitted() {
        let s=UncertaintyStudy::parse(STUDY).unwrap();
        assert_eq!(s.mean_control(),Some(MeanControlPolicy{max_solves:4,method:MeanControlMethod::CoordinateSecant}));
        assert_eq!(UncertaintyStudy::parse(s.canonical()).unwrap(),s);
        for text in [STUDY.replace(":version 3",":version 1"),
            STUDY.replace(":max-solves 4",":max-solves 3"),
            STUDY.replace(":max-solves 4",":max-solves 65"),
            STUDY.replace("coordinate-secant","adjoint"),
            STUDY.replace(":max-solves 4",":max-solves 4 :other 1"),
            STUDY.replace(" :mean-control (coordinate-secant :max-solves 4)","")] {
            assert!(UncertaintyStudy::parse(&text).is_err(),"{text}");
        }
    }
    #[test]
    fn probe_values_are_not_interpreted_as_probability_observations() {
        let p=params();
        assert!(coefficients(&p,&[1.0]).is_err());
        assert!(coefficients(&p,&[1.0,2.0,3.0,f64::NAN]).is_err());
        let mut p=p; p[0].low=0.0; p[0].high=f64::from_bits(1);
        assert!(coefficients(&p,&[0.0,1.0,2.0,3.0]).is_err());
    }
    #[test]
    fn fixed_count_qmc_and_copula_keep_their_original_laws() {
        let text=STUDY.replace("monte-carlo","quasi-monte-carlo :qmc (owen-scrambled-sobol :replicates 2 :samples-per-replicate 2)")
            .replace("independent","(gaussian-copula :latent-correlation ((1 1) (1 1)))");
        let s=UncertaintyStudy::parse(&text).unwrap();
        assert!(s.qmc().is_some()); assert!(s.compliance().is_none());
        assert_eq!(s.latent_correlation().unwrap(), &[vec![1.0,1.0],vec![1.0,1.0]]);
        assert_eq!(probe(s.parameters(),0).unwrap(),[2.0,300.0]);
    }
    #[test]
    fn large_finite_midpoint_does_not_overflow() {
        let p=vec![UniformParameter{name:"flux".into(),entity:"wall".into(),target:Target::HeatFlux,
            low:f64::MAX/2.0,high:f64::MAX}];
        assert!(means(&p)[0].is_finite());
        assert_eq!(coefficients(&p,&[0.0,0.0]).unwrap(),[0.0]);
    }
    #[test]
    fn nominal_adjoint_uses_one_mean_point_for_all_physical_parameters() {
        let text=STUDY.replace("coordinate-secant :max-solves 4","nominal-adjoint :max-solves 1");
        let s=UncertaintyStudy::parse(&text).unwrap();
        let policy=s.mean_control().unwrap();
        assert_eq!(policy.method,MeanControlMethod::NominalAdjoint);
        assert_eq!(policy.probe_count(s.parameters()),1);
        assert_eq!(policy.probe(s.parameters(),0).unwrap(),[4.0,300.0]);
        assert!(policy.probe(s.parameters(),1).is_err());
        assert!(UncertaintyStudy::parse(&text.replace(":max-solves 1",":max-solves 0")).is_err());
        assert_eq!(UncertaintyStudy::parse(s.canonical()).unwrap(),s);
    }
    #[test]
    fn nominal_singletons_skip_calibration_and_unsupported_variables_refuse() {
        let text=STUDY.replace("coordinate-secant :max-solves 4","nominal-adjoint :max-solves 0")
            .replace(":high 6W",":high 2W").replace(":high 310K",":high 290K");
        let s=UncertaintyStudy::parse(&text).unwrap();
        assert_eq!(s.mean_control().unwrap().probe_count(s.parameters()),0);
        let unsupported=STUDY.replace("coordinate-secant :max-solves 4","nominal-adjoint :max-solves 1")
            .replace("convection-temperature","air-inlet-temperature");
        assert!(UncertaintyStudy::parse(&unsupported).is_err());
        let accepted=unsupported.replace(":high 310K",":high 290K");
        assert!(UncertaintyStudy::parse(&accepted).is_ok());
    }
}
