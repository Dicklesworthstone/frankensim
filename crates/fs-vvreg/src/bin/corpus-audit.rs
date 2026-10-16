//! Deterministic audit of the seeded V&V corpus.

use std::io::Write as _;

fn main() {
    let report = fs_vvreg::corpus::corpus().audit();
    let rendered = report.render_table();
    if std::io::stdout()
        .lock()
        .write_all(rendered.as_bytes())
        .is_err()
    {
        std::process::exit(2);
    }
    if !report.is_clean() {
        std::process::exit(1);
    }
}
