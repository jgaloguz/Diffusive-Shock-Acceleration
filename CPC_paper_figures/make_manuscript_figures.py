"""Rebuild the seven manuscript figures from saved numerical outputs."""
from pathlib import Path
import argparse
import hashlib
import json
import shutil
from tempfile import TemporaryDirectory

import make_plots as plots


FIGURES = {
    "one_count_limit": "dsa_one_count_limit_readable",
    "path_intensity_available": "dsa_path_density_all_readable",
    "expected_amplification": "dsa_expected_oversampling_readable",
    "induced_drifts": "dsa_induced_drift_readable",
    "spectrum_comparison": "dsa_spectrum_comparison_ibm",
    "performance_metrics": "dsa_performance_metrics_ibm",
    "spatial_marginal_original_curves": "dsa_spatial_marginals_readable",
}


def main():
    folder = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=folder / "data")
    parser.add_argument("--output", type=Path, default=folder / "figures")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    plots.style()
    with TemporaryDirectory(prefix="dsa-figures-") as temporary:
        project = plots.Project(args.data, Path(temporary))
        plots.performance(project)
        plots.spectrum(project)
        plots.one_count(project)
        plots.likelihood(project)
        plots.path_intensity(project)
        plots.original_marginal_typography(project)

        for source, target in FIGURES.items():
            for extension in ("pdf", "png"):
                shutil.copy2(project.output / f"{source}.{extension}",
                             args.output / f"{target}.{extension}")
        shutil.copy2(project.output / "metrics_audit.csv",
                     args.output / "metrics_audit.csv")
        project.diagnostics["manuscript_filenames"] = FIGURES
        project.diagnostics["input_sha256"] = {
            str(path.relative_to(args.data)):
                hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sorted(args.data.rglob("*")) if path.is_file()
        }
        (args.output / "figure_manifest.json").write_text(
            json.dumps(project.diagnostics, indent=2) + "\n")

    print(f"Wrote the seven manuscript figures to {args.output.resolve()}")
    print("No simulation or backward solve was run.")


if __name__ == "__main__":
    main()
