//! JSON string escaping shared by native render entry points.
pub(super) fn json_string(value: &str) -> String {
    use std::fmt::Write as _;
    let mut out = String::new();
    out.push('"');
    for ch in value.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            ch if ch <= '\u{001f}' => {
                write!(&mut out, "\\u{:04x}", u32::from(ch)).expect("writing to String");
            }
            ch => out.push(ch),
        }
    }
    out.push('"');
    out
}
