//! Native Latin-hypercube layout admission; sampling belongs to fs-uq.

use super::{Node, Result, error, fields, integer, list, symbol};

/// Fixed replicated stratification, distinct from a power-of-two Sobol net.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LhsLayout {
    /// Independently randomized hypercubes, at least two.
    pub replicates: usize,
    /// Strata per coordinate per hypercube, at least two; arbitrary integers.
    pub samples_per_replicate: usize,
}

pub(super) fn parse(node: &Node, samples: usize) -> Result<LhsLayout> {
    let nodes = list(node)?;
    symbol(nodes.first().ok_or_else(|| error("empty LHS layout"))?, "replicated-latin-hypercube")?;
    let f = fields(&nodes[1..], &["replicates", "samples-per-replicate"])?;
    let replicates = usize::try_from(integer(f["replicates"])?)
        .map_err(|_| error("LHS replicate count overflow"))?;
    let samples_per_replicate = usize::try_from(integer(f["samples-per-replicate"])?)
        .map_err(|_| error("LHS stratum count overflow"))?;
    if !(2..=256).contains(&replicates) || samples_per_replicate < 2
        || replicates.checked_mul(samples_per_replicate) != Some(samples) {
        return Err(error("LHS requires at least two replicates times a stratum count >=2, exactly matching :samples"));
    }
    Ok(LhsLayout { replicates, samples_per_replicate })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::uncertainty::UncertaintyStudy;

    const LAYOUT: &str = ":lhs (replicated-latin-hypercube :replicates 2 :samples-per-replicate 3)";
    fn source() -> String {
        format!(r#"(fsim-uncertainty-study
            :version 1 :project "base.fsim" :samples 6 :seed 29 :wall-time 120s
            :method latin-hypercube {LAYOUT} :correlation independent :qoi "temperature-max"
            :geometry ((mesh :role "enclosure" :path "plate.stl" :unit "m" :max-hole-edges 0))
            :materials ("aa6061.fsmcdpk") :interfaces ()
            :parameters ((uniform :name "power" :target power :entity "air" :low 4W :high 6W)))"#)
    }

    #[test]
    fn g0_non_power_of_two_layout_round_trips_without_becoming_mc_or_qmc() {
        let study = UncertaintyStudy::parse(&source()).unwrap();
        assert_eq!(study.lhs(), Some(LhsLayout { replicates: 2, samples_per_replicate: 3 }));
        assert!(study.qmc().is_none());
        assert!(!study.sobol_sensitivity());
        assert!(study.compliance().is_none());
        assert!(study.mean_control().is_none());
        assert_eq!(UncertaintyStudy::parse(study.canonical()).unwrap(), study);
    }

    #[test]
    fn g0_layout_is_mandatory_exclusive_and_exactly_budgeted() {
        let original = source();
        for bad in [
            original.replace(LAYOUT, ""),
            original.replace(":method latin-hypercube", ":method monte-carlo"),
            original.replace(":samples 6", ":samples 7"),
            original.replace(":replicates 2", ":replicates 1"),
            original.replace(":samples-per-replicate 3", ":samples-per-replicate 1"),
            original.replace(LAYOUT, &format!("{LAYOUT} {LAYOUT}")),
            original.replace(LAYOUT, &format!("{LAYOUT} :qmc (owen-scrambled-sobol :replicates 2 :samples-per-replicate 3)")),
        ] {
            assert!(UncertaintyStudy::parse(&bad).is_err(), "{bad}");
        }
    }

    #[test]
    fn g0_dependence_and_sequential_policies_cannot_silently_change_the_design() {
        let original = source();
        for bad in [
            original.replace(":correlation independent", ":correlation (gaussian-copula :latent-correlation ((1)))"),
            original.replace(":version 1", ":version 2 :compliance (bernoulli-mixture :required-probability 0.95 :alpha 0.05 :min-samples 2)"),
            original.replace(":version 1", ":version 3 :mean-control ()"),
        ] {
            assert!(UncertaintyStudy::parse(&bad).is_err(), "{bad}");
        }
        let mc = original.replace(":method latin-hypercube", ":method monte-carlo").replace(LAYOUT, "");
        assert!(UncertaintyStudy::parse(&mc).unwrap().lhs().is_none());
    }

    #[test]
    fn g0_frozen_mean_controls_preserve_the_explicit_lhs_layout() {
        for method in ["coordinate-secant :max-solves 2", "nominal-adjoint :max-solves 1"] {
            let text = source().replace(":version 1",
                &format!(":version 3 :mean-control ({method})"));
            let study = UncertaintyStudy::parse(&text).unwrap();
            assert_eq!(study.lhs(), Some(LhsLayout { replicates: 2, samples_per_replicate: 3 }));
            assert!(study.mean_control().is_some());
            assert!(study.qmc().is_none() && study.compliance().is_none());
            assert_eq!(UncertaintyStudy::parse(study.canonical()).unwrap(), study);
            assert!(UncertaintyStudy::parse(&text.replace(":version 3", ":version 1")).is_err());
        }
    }
}
