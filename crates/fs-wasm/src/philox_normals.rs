//! Annus Mirabilis `philox_normals` export: typed envelope over fs-rand.
//!
//! An unexplained empty `Vec<f64>` is not a refusal.

pub use fs_rand::philox_normals::{
    admit_philox_normals, philox_normals, philox_normals_admitted, Refusal, TargetKind,
    KERNEL_VERSION, PHILOX_NORMALS_MAX_COUNT, PhiloxNormalsSpec,
};

/// JSON envelope for tests and for hosts that do not consume `JsValue`.
/// Never `[]`.
pub fn philox_normals_envelope_json(
    seed: u64,
    stream_kernel: u32,
    tile: u32,
    start_index: u64,
    count: usize,
) -> String {
    match philox_normals(seed, stream_kernel, tile, start_index, count) {
        Ok(values) => {
            let mut body = String::from("{\"ok\":{");
            body.push_str(&format!("\"kernel\":\"{KERNEL_VERSION}\","));
            body.push_str("\"export\":\"philox_normals\",");
            body.push_str("\"layout\":{\"length\":");
            body.push_str(&values.len().to_string());
            body.push_str(",\"indexRule\":\"draws\"},");
            body.push_str("\"valuesBits\":[");
            for (i, z) in values.iter().enumerate() {
                if i > 0 {
                    body.push(',');
                }
                body.push_str(&format!("\"{:016x}\"", z.to_bits()));
            }
            body.push_str("]}}");
            body
        }
        Err(r) if r.target_kind == TargetKind::ExecutionOutcome => format!(
            "{{\"execution\":{{\"code\":\"{}\",\"message\":\"{}\",\"ranked_repairs\":[{}],\"details\":{}}}}}",
            r.code,
            escape_json(&r.message),
            r.ranked_repairs
                .iter()
                .map(|s| format!("\"{}\"", escape_json(s)))
                .collect::<Vec<_>>()
                .join(","),
            r.details
        ),
        Err(r) => format!(
            "{{\"refusal\":{{\"code\":\"{}\",\"message\":\"{}\",\"ranked_repairs\":[{}],\"details\":{}}}}}",
            r.code,
            escape_json(&r.message),
            r.ranked_repairs
                .iter()
                .map(|s| format!("\"{}\"", escape_json(s)))
                .collect::<Vec<_>>()
                .join(","),
            r.details
        ),
    }
}

fn escape_json(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for ch in s.chars() {
        match ch {
            '\\' => out.push_str("\\\\"),
            '"' => out.push_str("\\\""),
            '\n' => out.push_str("\\n"),
            c => out.push(c),
        }
    }
    out
}
