//! Shared production policy for steady solves and immutable-old-state transient
//! coupling trials. No hidden retry, added physics solve, or enlarged budget.
use super::*;
use fs_airflow::conjugate::IqnIlsConfig;
use fs_airflow::graph::thermal::transport::TransportError;

pub(super) const POLICY: IqnIlsConfig = IqnIlsConfig {
    max_history: 8,
    relative_rank_tolerance: 1e-10,
};

/// Same callback seam for the steady and transient producers. Each call starts
/// fresh history: a different timestep, workload, fan speed or material state
/// must never inherit secants for a different interface map.
pub(super) fn solve_coupled_transport<F>(
    cx: &Cx<'_>, network: &TransportNetwork<'_>, config: &ConjugateConfig, solid: F,
) -> std::result::Result<CoupledTransportSolution, TransportError>
where F: FnMut(&Cx<'_>, &[f64]) -> std::result::Result<Vec<SolidRegionState>, AirflowError>,
{
    fs_airflow::graph::thermal::coupled_transport::solve_coupled_transport_iqn(
        cx, network, config, POLICY, solid,
    )
}

pub(super) fn render(relaxation: f64, has_adjoint: bool) -> Result<String> {
    Ok(format!(
        "{{\"method\":\"iqn-ils\",\"scope\":\"complete-mixed-network-interface\",\"max_history\":{},\"relative_rank_tolerance\":{},\"fallback\":\"declared-relaxation\",\"fallback_omega\":{},\"adjoint_method\":{},\"acceptance\":\"fresh unrelaxed interface residuals and independent branch/solid energy gates; derivative residual checked separately\"}}",
        POLICY.max_history, num(POLICY.relative_rank_tolerance)?, num(relaxation)?,
        if has_adjoint { "\"iqn-ils\"" } else { "null" },
    ))
}
