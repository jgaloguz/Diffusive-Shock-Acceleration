#!/usr/bin/env python3
"""Replot saved DSA outputs; this script never reruns, smooths, or resamples simulations.

Run: python make_plots.py --data data --output figures
The figures use Robert's IBM color-blind-safe line-and-marker palette.
See README.txt and diagnostics.json for data qualifications.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import re

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
from matplotlib.ticker import NullLocator
from pypdf import PdfReader, PdfWriter

IBM = ["#648fff", "#fe6100", "#ffb000", "#dc267f", "#785ef0"]
IBM_ORDER = [IBM[0], IBM[4], IBM[3], IBM[1], IBM[2]]
MARKERS = ["o", "s", "^", "D", "v"]
BLACK = "#202020"
CONV = 1.0 / 2.99792414e4  # Saved SPECTRUM F -> cm^-3 (MeV/c)^-1.
AU_CM = 1.495978707e13
MP_MEV = 938.2720813
P0_MEV = np.sqrt(1.0 * (1.0 + 2.0 * MP_MEV))
TAU_SECONDS = 30 * 86400.0
DENSITY_LABEL = r"Particle density (cm$^{-3}$ (MeV/$c$)$^{-1}$)"
FAMILIES = [
    ("split_0", "Continuous splitting"),
    ("split_1", "Discrete splitting"),
    ("imps_0", "Splitting-matched IS"),
    ("imps_1", "Shock-directed IS"),
]
METHODS = [
    ("base", "Direct integration", IBM_ORDER[0], "o"),
    ("split_0", "Continuous splitting", IBM_ORDER[1], "s"),
    ("split_100", "Discrete splitting", IBM_ORDER[2], "^"),
    ("imps_0", "Splitting-matched IS", IBM_ORDER[3], "D"),
    ("imps_1", "Shock-directed IS", IBM_ORDER[4], "v"),
]


def style():
    plt.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["DejaVu Sans"],
        "font.size": 12.5, "axes.titlesize": 12.5, "axes.labelsize": 12.5,
        "font.weight": "normal", "axes.titleweight": "normal",
        "axes.labelweight": "normal", "axes.grid": False,
        "xtick.labelsize": 11.5, "ytick.labelsize": 11.5, "legend.fontsize": 11.5,
        "axes.linewidth": .65, "axes.spines.top": False,
        "axes.spines.right": False, "axes.axisbelow": True,
        "pdf.fonttype": 42, "ps.fonttype": 42,
        "savefig.facecolor": "white", "savefig.dpi": 220,
        "legend.frameon": False, "legend.handlelength": 2.0,
        "mathtext.fontset": "dejavusans",
    })


def figure_axes(*args, **kwargs):
    """Keep type sizes consistent after figures are fitted to manuscript width."""
    scale = kwargs.get("figsize", (6.7, 4))[0] / 6.7
    plt.rcParams.update({
        "font.size": 12.5*scale, "axes.titlesize": 12.5*scale,
        "axes.labelsize": 12.5*scale, "xtick.labelsize": 11.5*scale,
        "ytick.labelsize": 11.5*scale, "legend.fontsize": 11.5*scale,
        "legend.title_fontsize": 11.5*scale,
    })
    return plt.subplots(*args, **kwargs)


class Project:
    def __init__(self, data, output):
        self.data = data
        self.output = output
        output.mkdir(parents=True, exist_ok=True)
        self.x = self.load("dsa_analytic_pos.dat")
        self.E = self.load("dsa_analytic_enr.dat")
        self.El = self.load("dsa_analytic_enr_left.dat")
        self.p = self.load("dsa_analytic_mom.dat")  # p / p0, NOT physical p.
        self.dp = self.load("dsa_analytic_dmom.dat")  # Delta p / p0.
        self.times = self.load("dsa_analytic_time.dat")[:4]
        self.dz = float(np.diff(self.x)[0])
        self.figures = []
        self.captions = {}
        self.diagnostics = {"units": {"density_factor": CONV,
            "density": "particles cm^-3 (MeV/c)^-1", "horizontal_coordinate":
            "Proton kinetic energy; F remains a density per momentum."}}
        self.metrics = {}
        self.collect_metrics()

    def load(self, name):
        a = np.loadtxt(self.data / name)
        if not np.all(np.isfinite(a)):
            raise ValueError(f"Nonfinite input: {name}")
        return a

    def paths(self, name):
        return self.load(f"dsa_forward_path_dens_pp_{name}.dat").reshape(4, 100, 100)

    def nseed(self, folder, logstem):
        text = (self.data / folder / "slurm_output" / (logstem + ".err")).read_text()
        nr = int(re.search(r"Number of CPUs: (\d+)", text)[1])
        ns = int(re.search(r"Number of trajectories per CPU: (\d+)", text)[1])
        return nr * ns

    def collect_metrics(self):
        rows = []
        for stem, label in FAMILIES:
            folder = f"e2metric_{stem}"
            costs = self.load(f"{folder}/steps_per_particle.dat")
            records = []
            for param, cost, cost_batch_sd in costs:
                if param == 0:
                    file = "dsa_forward_mom_pp_base.dat"
                    logstem = "base"
                else:
                    fmt = f"{param:.1f}" if stem == "imps_1" else f"{param:.2f}"
                    file = f"dsa_forward_mom_pp_{stem}_{fmt}.dat"
                    logstem = ("A0_" if stem == "imps_1" else "alpha_") + fmt
                a = self.load(f"{folder}/{file}")
                N = self.nseed(folder, logstem)
                r = a[:, 5] / a[:, 4]
                # Raw moments of the independent seed contribution Y:
                # E[(Y/H - 1)^2] = e2, E[Y/H] = r.
                second = a[:, 6] + 2 * r - 1
                var = np.maximum(a[:, 6] - (r - 1)**2, 0)
                se = np.sqrt(var / (N - 1))
                ess = N * r**2 / second
                rec = dict(param=float(param), cost=float(cost), a=a,
                    e2=a[:, 6], c=cost*a[:, 6], N=N, ratio=r,
                    mean_se=se, effective_contributors=ess)
                records.append(rec)
            for rec in records:
                rec["gain"] = records[0]["c"] / rec["c"]
                rows.append({"method": label, "parameter": rec["param"],
                    "independent_seeds": rec["N"], "steps_per_seed": rec["cost"],
                    "gain_at_300GeV": np.interp(3e5, self.El, rec["gain"]),
                    "H_over_reference_at_300GeV": np.interp(3e5, self.El, rec["ratio"]),
                    "MC_standard_error_at_300GeV": np.interp(3e5, self.El, rec["mean_se"]),
                    "effective_contributors_at_300GeV": np.interp(3e5, self.El, rec["effective_contributors"]),
                    "effective_contributors_at_904GeV": rec["effective_contributors"][-1]})
            self.metrics[stem] = records
        with (self.output / "metrics_audit.csv").open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=rows[0].keys())
            writer.writeheader()
            writer.writerows(rows)
        self.diagnostics["metrics"] = rows

    def save(self, fig, name, caption):
        fig.savefig(self.output / (name + ".pdf"), bbox_inches="tight", pad_inches=.06)
        fig.savefig(self.output / (name + ".png"), bbox_inches="tight", pad_inches=.06)
        self.figures.append(name)
        self.captions[name] = caption
        plt.close(fig)


def log_energy(ax, lower_edge=False):
    ax.set_xscale("log")
    ax.set_xlim(1, 1e6)
    ax.set_xticks([1, 1e2, 1e4, 1e6])
    ax.xaxis.set_minor_locator(NullLocator())


def param_line(ax, x, y, i, label=None, **kw):
    ax.plot(x, y, color=IBM[i], marker=MARKERS[i], markersize=3.6,
        markeredgewidth=.35, markeredgecolor="white", linewidth=1.3,
        markevery=(i, 8), label=label, **kw)


def parameter_handles(records, shock=False):
    handles = [Line2D([], [], color=IBM[i], marker=MARKERS[i], linewidth=1.4,
        markersize=4, label=(r"$A_0=" + f'{r["param"]:g}' + r"$" if shock
            else r"$\alpha=" + f'{r["param"]:g}' + r"$"))
        for i, r in enumerate(records[1:])]
    return handles


def shock_legend(ax, records):
    ax.legend(handles=parameter_handles(records, shock=True), ncol=2,
        loc="upper left", handlelength=1.3, columnspacing=.9,
        borderaxespad=.5)


def shared_shock_legend(fig, records):
    fig.legend(handles=parameter_handles(records, shock=True), ncol=4,
        loc="outside lower center", handlelength=1.2, columnspacing=.7)


def shared_parameter_legend(fig, records, direct=False):
    handles = parameter_handles(records)
    if direct:
        handles = [Line2D([], [], color=BLACK, linestyle="--", linewidth=1.5,
            label="Direct")] + handles
    fig.legend(handles=handles, loc="outside upper center",
        ncol=3 if len(handles)>5 else len(handles), handlelength=1.3,
        columnspacing=.9)


def performance(p):
    fig, placeholder = figure_axes(figsize=(7.3, 9.0), layout="constrained")
    placeholder.remove()
    left, right, legends = fig.subfigures(1, 3,
        width_ratios=[1,1,.42], wspace=.015)
    axs = np.column_stack([left.subplots(4,1,sharex=True),
        right.subplots(4,1,sharex=True)])
    for row, (stem, title) in enumerate(FAMILIES):
        records = p.metrics[stem]
        linestyle = "-"
        for col, key in enumerate(["e2", "c"]):
            ax = axs[row, col]
            ax.set_yscale("log")
            ax.plot(p.El, records[0][key], color=BLACK, linestyle=linestyle,
                linewidth=1.5, zorder=2)
            for i, rec in enumerate(records[1:]):
                ax.plot(p.El, rec[key], color=IBM_ORDER[i], linestyle=linestyle,
                    marker=MARKERS[i], markersize=3.6, linewidth=1.3,
                    markeredgewidth=.35, markeredgecolor="white", markevery=(i,8))
            log_energy(ax, True)
            if col == 0:
                ax.set_ylim(1, 1e8)
                ax.set_yticks([1, 1e2, 1e4, 1e6, 1e8])
                ax.set_title(title, loc="left")
            else:
                ax.set_ylim(1e5, 1e13)
                ax.set_yticks([1e5, 1e7, 1e9, 1e11, 1e13])
                ax.set_title(" ",loc="left")
            ax.yaxis.set_minor_locator(NullLocator())
            if row == 3:
                ax.set_xlabel("$E$ (MeV)")
    left.supylabel("Mean squared relative error",fontsize=plt.rcParams["axes.labelsize"])
    right.supylabel("Cost-adjusted error (integration steps)",fontsize=plt.rcParams["axes.labelsize"])
    lgs = legends.add_gridspec(4,1)
    alpha_ax = legends.add_subplot(lgs[:3,0])
    shock_ax = legends.add_subplot(lgs[3,0])
    for ax in (alpha_ax,shock_ax):
        ax.set_axis_off()
    alpha_handles = [Line2D([],[],color=BLACK,linewidth=1.5,label=r"$\alpha=0$")]
    shock_handles = [Line2D([],[],color=BLACK,linewidth=1.5,label=r"$A_0=0$")]
    for i,r in enumerate(p.metrics["split_0"][1:]):
        alpha_handles.append(Line2D([],[],color=IBM_ORDER[i],marker=MARKERS[i],
            linewidth=1.4,markersize=4,label=r"$\alpha="+f'{r["param"]:g}'+r"$"))
    for i,r in enumerate(p.metrics["imps_1"][1:]):
        shock_handles.append(Line2D([],[],color=IBM_ORDER[i],marker=MARKERS[i],
            linewidth=1.4,markersize=4,label=r"$A_0="+f'{r["param"]:g}'+r"$"))
    alpha_ax.legend(handles=alpha_handles,loc="upper left",borderaxespad=0,
        handlelength=1.3,handletextpad=.5)
    shock_legend = shock_ax.legend(handles=shock_handles,loc="upper left",borderaxespad=0,
        handlelength=1.3,handletextpad=.5)
    # Align the last legend with the top of the fourth data row, even though
    # its containing subfigure has a separate constrained layout.
    fig.canvas.draw()
    row_top = axs[3,1].get_window_extent().y1
    legend_top = shock_ax.transAxes.inverted().transform((0,row_top))[1]
    shock_legend.set_bbox_to_anchor((0,legend_top),transform=shock_ax.transAxes)
    p.save(fig, "performance_metrics", "Single-seed mean squared relative error (left) and cost-adjusted error (right), using every saved momentum threshold. The latter is mean integration steps per seed times single-seed error, measured in integration steps. E is the proton kinetic energy corresponding to the lower momentum threshold p. All curves are solid; the upper right alpha legend applies to the first three rows, and the lower right A0 legend applies to the last row, with A0 measured in au^-1. Black curves show direct stochastic integration (zero bias); the colored parameters increase in IBM order: blue, purple, magenta, orange, yellow. All values remain unsmoothed; markers are shown on every eighth threshold. Discrete splitting uses 100 thresholds.")

    fig, axs = figure_axes(2, 2, figsize=(7.3, 5.4), sharex=True, sharey=True, layout="constrained")
    for ax, (stem, title) in zip(axs.flat, FAMILIES):
        recs = p.metrics[stem]
        ax.set_title(title, loc="left")
        for i, rec in enumerate(recs[1:]):
            param_line(ax, p.El, rec["gain"], i)
        ax.axhline(1, color=BLACK, linestyle="--", linewidth=1.1)
        ax.set_yscale("log")
        ax.set_ylim(1e-2, 1e5)
        ax.set_yticks([1e-2, 1, 1e2, 1e4])
        ax.yaxis.set_minor_locator(NullLocator())
        log_energy(ax, True)
    for ax in axs[-1]:
        ax.set_xlabel("$E$ (MeV)")
    for ax in axs[:, 0]:
        ax.set_ylabel("Efficiency gain")
    shared_parameter_legend(fig, p.metrics["split_0"])
    shared_shock_legend(fig, p.metrics["imps_1"])
    p.save(fig, "efficiency_gains", "Step-based efficiency gain c_direct/c_method from the same saved data as performance_metrics. Values above the black dashed line at one favor enhanced sampling. E is the proton kinetic energy corresponding to the lower momentum threshold p. The shared alpha legend applies to continuous splitting, discrete splitting, and matched IS; the shock-directed panel has its own A0 legend, with A0 measured in au^-1. All panels use identical scales. This is a compact alternative to the eight-panel metric figure, not an additional experiment or a wall-clock comparison.")


def spectrum(p):
    fig, ax = figure_axes(figsize=(6.7, 4.2), layout="constrained")
    handles, budgets = [], []
    for stem, label, color, marker in METHODS:
        a = p.load(f"spectrum_comparison/dsa_forward_mom_pp_{stem}.dat")
        keep = a[:, 1] > 0
        ax.plot(p.E[keep], CONV*a[keep, 1], linestyle="none", marker=marker,
            color=color, markersize=5.8, markeredgewidth=.45,
            markeredgecolor="white", alpha=.9, zorder=3)
        handles.append(Line2D([], [], color=color, marker=marker, linestyle="none",
            markersize=6.3, label=label))
        log = (p.data / "spectrum_comparison/slurm_output" / (stem+".out")).read_text()
        cost = float(re.search(r"Steps per seed particle = ([0-9.]+)", log)[1])
        N = p.nseed("spectrum_comparison", stem)
        budgets.append({"method": label, "independent_seeds": N,
            "steps_per_seed": cost, "total_steps": N*cost})
    reference = p.load("spectrum_comparison/dsa_forward_mom_pp_base.dat")[:, 2]
    ax.plot(p.E, CONV*reference, color=BLACK, linewidth=1.4, zorder=4)
    handles.append(Line2D([], [], color=BLACK, linewidth=1.4, label="Sharp-shock reference"))
    ax.set_yscale("log")
    ax.set_ylim(1e-21, 8e-9)
    ax.set_yticks([1e-9, 1e-12, 1e-15, 1e-18, 1e-21])
    ax.set_ylabel(DENSITY_LABEL)
    ax.set_xlabel("Proton kinetic energy (MeV)")
    log_energy(ax)
    ax.yaxis.set_minor_locator(NullLocator())
    ax.legend(handles=handles, loc="lower left",
        labelspacing=.35, borderaxespad=.7)
    p.diagnostics["spectrum_budgets"] = budgets
    p.save(fig, "spectrum_comparison", "Saved spectrum calculations at 30 days in the spatial bin [0.1,0.2] au. The curve is the saved pointwise sharp-shock reference at x=0.15 au; no new bin averaging of the analytic solution is used. F is particle density per unit momentum while the horizontal coordinate E is proton kinetic energy. Zero-count bins have no point on the logarithmic scale. CAUTION: the supplied logs imply 9.12e9 steps for matched IS versus 2.42e9--2.86e9 for the other curves; these outputs do not have equal step budgets. All enhanced calculations use alpha=1 or A0=1 au^-1. See diagnostics.json.")


def one_count(p):
    # A known 128,000-seed saved direct run replaces the unavailable original
    # 12,800-seed illustration. No histogram thinning or synthetic points.
    a = p.load("spectrum_comparison/dsa_forward_mom_pp_base.dat")
    N = p.nseed("spectrum_comparison", "base")
    dp_phys = p.dp * P0_MEV
    one = TAU_SECONDS / (N * p.dz * AU_CM * dp_phys)
    fig, ax = figure_axes(figsize=(6.7, 4.2), layout="constrained")
    ax.loglog(p.E, CONV*a[:, 2], color=BLACK, linestyle="--", linewidth=1.9,
        label="Sharp-shock reference")
    keep = a[:, 1] > 0
    ax.loglog(p.E[keep], CONV*a[keep, 1], linestyle="none", marker="o",
        markersize=6.4, color=IBM[0], markeredgecolor="white", markeredgewidth=.5,
        label="Direct integration")
    for mult, color in [(1, IBM[2]), (2, IBM[1]), (3, IBM[3])]:
        ax.loglog(p.E, mult*one, color=color, linestyle="--", linewidth=1.3,
            label=f"{mult} {'count' if mult == 1 else 'counts'}")
    log_energy(ax)
    ax.set_ylabel(DENSITY_LABEL)
    ax.set_xlabel("Proton kinetic energy (MeV)")
    ax.set_ylim(1e-20, 1e-8)
    ax.yaxis.set_minor_locator(NullLocator())
    ax.legend(loc="lower left", ncol=1)
    p.save(fig, "one_count_limit", f"Particle density F versus proton kinetic energy, with physical units on the axes. Blue circles show the supplied {N:,}-seed direct run, and the dashed black curve is the pointwise sharp-shock reference. Yellow, orange, and magenta dashed curves show the one-, two-, and three-count levels, respectively, for 100 momentum bins and a 0.1-au spatial bin. The current root direct output has been overwritten by a much larger run; this figure uses the identifiable spectrum_comparison run. Its N differs from the manuscript's original 12,800-seed illustration, so the caption must be updated if substituted.")


def likelihood(p):
    a = p.load("feynman_kac/dsa_likelihood_space_time_alpha=1.00.dat")
    x, h = a[0], a[1:]
    time_colors = [IBM[0], IBM[4], IBM[3], IBM[1]]
    fig, ax = figure_axes(figsize=(6.7, 3.7), layout="constrained")
    for i, (t, row) in enumerate(zip(p.times, h)):
        ax.plot(x, row, color=time_colors[i], linewidth=1.7,
            label=f"{int(t)} {'day' if t == 1 else 'days'}")
    ax.set_xlim(-12, 1)
    ax.set_ylim(0, 8.7)
    ax.set_xlabel("$x$ (au)")
    ax.set_ylabel("Future-oversampling factor")
    ax.legend(loc="upper left", ncol=1)
    p.save(fig, "expected_amplification", "Expected future amplification hbar(t,x) for alpha=1 and terminal time 30 days, directly from the saved 10,000-point spatial slices. No backward solve was rerun.")

    # Same three-point, second-order nonuniform-grid derivative as Juan's solver.
    lnh = np.log(h)
    grad = np.gradient(lnh, x, axis=1, edge_order=2)
    u = 2.5 - 1.5*np.tanh(4*x/.01)
    kappa = .668459*u*u/16
    drift = 2*kappa[None, :]*grad
    shock_drift = -2*kappa*np.tanh(x/.01)
    fig, ax = figure_axes(figsize=(6.7, 3.7), layout="constrained")
    for i, (t, row) in enumerate(zip(p.times, drift)):
        ax.plot(x, row, color=time_colors[i], linewidth=1.6,
            label=f"{int(t)} {'day' if t == 1 else 'days'}")
    ax.plot(x, shock_drift, color=BLACK, linestyle="--", linewidth=1.3,
        label="Shock-directed")
    ax.set_xlabel("$x$ (au)")
    ax.set_ylim(-1.85, 1.45)
    ax.set_xlim(-12, 1)
    ax.set_ylabel(r"Drift ($10^7$ cm s$^{-1}$)")
    ax.legend(loc="lower left", ncol=1)
    p.save(fig, "induced_drifts", "Artificial drift V from the saved hbar slices, with alpha=1 and terminal time 30 days. The black dashed curve is shock-directed drift with A0=1 au^-1. Spatial derivatives are recomputed by the same nonuniform-grid three-point formula as the original solver; no new PDE solve is used.")


def check_continuous_path_file(p):
    paths = p.paths("split_0")
    mom = p.load("dsa_forward_mom_pp_split_0.dat")
    # The unweighted/weighted ratio in one spatial bin is approximately ell(p).
    # This estimates the configuration; it is NOT a precise alpha metadata field.
    idx = int(np.argmin(np.abs(p.x-.15)))
    ratio = paths[-1, idx] / np.maximum(mom[:, 1], 1e-200)
    keep = (mom[:, 1] > 0) & (paths[-1, idx] > 0) & (p.p > 2) & (p.p < 300)
    exponent = float(np.polyfit(np.log(p.p[keep]), np.log(ratio[keep]), 1)[0])
    mass = np.sum(paths*p.dp[None, None, :], axis=(1, 2))*p.dz
    direct = np.sum(p.paths("base")*p.dp[None, None, :], axis=(1, 2))*p.dz
    p.diagnostics["continuous_path_file"] = {
        "estimated_amplification_exponent": exponent,
        "terminal_intensity_ratio_to_direct": float(mass[-1]/direct[-1]),
        "expected_ratio_for_alpha_1_in_manuscript": "about 5",
        "last_modified": "2026-09-24T12:45:32Z (later than the September 22 figure)",
        "status": "The legacy root file appears overwritten by an alpha=1.5 run and is excluded from alpha=1 comparisons. Appendix B instead uses the separately named alpha=1 Expanse rerun with 12.8 million seeds."}
    return exponent


def path_intensity(p):
    check_continuous_path_file(p)
    rows = [("base", "Direct integration"),
        ("split_10", "Discrete splitting"),
        ("imps_0", "Splitting-matched IS"),
        ("imps_1", "Shock-directed IS")]
    # A warm, perceptually ordered scalar map; IBM colors remain on line plots.
    cmap = plt.get_cmap("inferno").copy()
    cmap.set_bad("white")
    xedges = np.r_[p.x-.5*p.dz, p.x[-1]+.5*p.dz]
    Eedges = np.r_[p.El, 1e6]
    fig, axs = figure_axes(4, 4, figsize=(8.0, 8.0), sharex=True, sharey=True,
        layout="constrained")
    # Reserve row-heading space without allowing a long heading to stretch
    # the first column or create large gaps between the heat-map panels.
    fig.set_layout_engine(None)
    fig.subplots_adjust(left=.13, right=.99, bottom=.17, top=.89,
        wspace=.12, hspace=.6)
    norm = LogNorm(vmin=1e-16, vmax=2e-8)
    for r, (stem, label) in enumerate(rows):
        a = CONV*p.paths(stem)  # Apply conversion exactly once.
        for c, t in enumerate(p.times):
            ax = axs[r, c]
            pcm = ax.pcolormesh(xedges, Eedges,
                np.ma.masked_less_equal(a[c].T, 0), cmap=cmap, norm=norm,
                shading="flat", rasterized=True)
            ax.set_yscale("log")
            ax.set_xlim(-3.1, 3.1)
            ax.set_ylim(1, 1e6)
            ax.set_xticks([-2, 0, 2])
            ax.set_yticks([1, 1e2, 1e4, 1e6])
            ax.yaxis.set_minor_locator(NullLocator())
            if r == 0:
                ax.set_title(f"{int(t)} {'day' if t == 1 else 'days'}", pad=35)
            if c == 0:
                ax.set_ylabel("$E$ (MeV)")
            if r == 3:
                ax.set_xlabel("$x$ (au)")
        axs[r, 0].text(0, 1.035, label, transform=axs[r, 0].transAxes,
            va="bottom")
    cbax = fig.add_axes([.50,.04,.48,.015])
    cb = fig.colorbar(pcm, cax=cbax, orientation="horizontal",
        ticks=[1e-16, 1e-14, 1e-12, 1e-10, 1e-8])
    cb.ax.xaxis.set_minor_locator(NullLocator())
    fig.text(.005,.0475,r"Path intensity (cm$^{-3}$ (MeV/$c$)$^{-1}$)",
        ha="left",va="center",fontsize=plt.rcParams["axes.labelsize"])
    p.save(fig, "path_intensity_available", "Unweighted ancestral path intensities in particles cm^-3 (MeV/c)^-1, counted per terminal lineage and scaled per initial seed as in the original outputs. Columns show 1, 3, 10, and 30 days. Discrete splitting uses 10 thresholds and alpha=1; matched IS uses alpha=1; shock-directed IS uses A0=1 au^-1. Momentum-density units are retained while displaying proton kinetic energy E. All panels share one logarithmic warm inferno color scale; zero bins are white, and histogram cells use flat shading with explicit bin edges. Physical-unit conversion is applied once, correcting the double conversion in the repository's --method all plotting branch. IMPORTANT: this is an alternative four-row comparison: matched IS replaces continuous splitting because the supplied continuous-path file appears to have been overwritten by an alpha=1.5 run. The script does not silently relabel that file. Spatial display is cropped to [-3.1,3.1] au; normalization still uses the original complete output.")


def marginal(p):
    # Retain the two-method diagnostic alongside the full three-way appendix plot.
    sets = [("ltest", "Terminal-weighted direct", BLACK, "o"),
        ("imps_0", "Splitting-matched IS", IBM[4], "D")]
    qs, masses = {}, {}
    for stem, _, _, _ in sets:
        a = np.sum(p.paths(stem)*p.dp[None, None, :], axis=2)
        mass = np.sum(a, axis=1)*p.dz
        qs[stem] = a/mass[:, None]
        masses[stem] = mass.tolist()
        assert np.allclose(np.sum(qs[stem], axis=1)*p.dz, 1)
    fig, axs = figure_axes(2, 2, figsize=(7.3, 4.9), sharex="col",
        sharey="row", gridspec_kw={"height_ratios": [2, 1]}, layout="constrained")
    for c, ti in enumerate([1, 3]):
        ax = axs[0, c]
        for stem, label, color, marker in sets:
            ax.semilogy(p.x, qs[stem][ti], color=color, linewidth=1.5,
                marker=marker, markersize=3.3, markevery=3, label=label)
        ax.set_title(f"{int(p.times[ti])} days", loc="left")
        ax.set_ylim(1e-3, 5)
        ref = qs["ltest"][ti]
        keep = np.logical_and.reduce([qs[stem][ti] > 1e-3 for stem, *_ in sets])
        diff = np.full(ref.shape, np.nan)
        diff[keep] = 100*(qs["imps_0"][ti, keep]/ref[keep]-1)
        axs[1, c].plot(p.x, diff, color=IBM[4], linewidth=1.1, marker="D",
            markersize=3, markevery=3)
        axs[1, c].axhline(0, color=BLACK, linewidth=.8)
        axs[1, c].set_xlabel("$x$ (au)")
        axs[1, c].set_xlim(-1.7, 3)
    axs[0, 0].set_ylabel(r"$q$ (au$^{-1}$)")
    axs[1, 0].set_ylabel("Difference (%)")
    handles, labels = axs[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside upper center", ncol=2)
    p.diagnostics["marginal_integrals_code_units"] = masses
    p.save(fig, "spatial_marginal_available", "Normalized spatial marginals for terminal-weighted direct paths and matched IS at alpha=1. The lower panels show SIGNED fractional differences in percent, so agreement and direction of error are both readable; all bins with either q<=1e-3 au^-1 are masked. The full three-method comparison, including the recovered alpha=1 splitting output, is provided in spatial_marginal_original_curves. Normalization is over the full [-5,5]-au histogram, not over the plotted window. No uncertainty intervals can be recovered for these spatial histograms from the saved aggregate grids alone.")


def original_marginal_typography(p):
    """Draw all three alpha=1 curves from numerical histograms.

    The stable output name is retained for existing manuscript references.
    No original-image pixels, reconstructed paths, or alpha=1.5 outputs are used.
    """
    recovered = "split_0_alpha_1.00_expanse_20261123"
    sets = [("ltest", "Weighted direct", IBM[3]),
            (recovered, "Continuous splitting", IBM[0]),
            ("imps_0", "Splitting-matched IS", IBM[1])]
    qs, masses = {}, {}
    for stem, _, _ in sets:
        spatial = (p.paths(stem)*p.dp[None,None,:]).sum(axis=2)
        masses[stem] = spatial.sum(axis=1)*p.dz
        qs[stem] = spatial/masses[stem][:,None]
        assert np.allclose(qs[stem].sum(axis=1)*p.dz,1)
    fig, axs = figure_axes(2, 2, figsize=(7.3, 5.5), sharex=True,
        sharey="row", layout="constrained")
    fig.get_layout_engine().set(hspace=.06)
    xlimits = [(-1.2,2.7)]*2
    xticks = [[-1,0,1,2]]*2
    density_limits = [(1e-3,6)]*2
    difference_limits = [(.01,45)]*2
    difference_ticks = [[.01,.1,1,10]]*2
    for col, ti in enumerate([1,3]):
        ref = qs["ltest"][ti]
        axs[0,col].semilogy(p.x,np.where(ref>0,ref,np.nan),
            color=sets[0][2],linewidth=1.1)
        # Distinct marker shapes remain visible when values coincide.
        for (stem,_,color),marker in zip(sets[1:],["o","x"]):
            axs[0,col].semilogy(p.x,np.where(qs[stem][ti]>0,qs[stem][ti],np.nan),
                color=color,linestyle="none",marker=marker,markersize=4.4,
                markerfacecolor="none",markeredgewidth=.9)
        axs[0,col].set_title(f"{int(p.times[ti])} days",loc="center")
        keep = np.logical_and.reduce([qs[stem][ti]>1e-3 for stem,_,_ in sets])
        for (stem,_,color),marker in zip(sets[1:],["o","x"]):
            diff = np.full(p.x.shape,np.nan)
            diff[keep] = np.abs(qs[stem][ti,keep]/ref[keep]-1)
            axs[1,col].semilogy(p.x,100*diff,color=color,linestyle="none",
                marker=marker,markersize=4.4,markerfacecolor="none",markeredgewidth=.9)
        axs[1,col].set_xlabel("$x$ (au)")
        for ax in axs[:,col]:
            ax.set_xlim(*xlimits[col])
            ax.set_xticks(xticks[col])
            ax.set_xticklabels([f"{v:g}" for v in xticks[col]])
            ax.yaxis.set_minor_locator(NullLocator())
        axs[0,col].set_ylim(*density_limits[col])
        axs[0,col].set_yticks([1e-3,1e-2,1e-1,1])
        axs[1,col].set_ylim(*difference_limits[col])
        axs[1,col].set_yticks(difference_ticks[col])
        axs[1,col].set_yticklabels([f"{v:g}" for v in difference_ticks[col]])
    axs[0,0].set_ylabel(r"Spatial density (au$^{-1}$)")
    axs[1,0].set_ylabel("Relative difference (%)")
    handles = [Line2D([],[],color=sets[0][2],linewidth=1.1,label=sets[0][1])]
    handles += [Line2D([],[],color=color,marker=marker,linestyle="none",markersize=5,
        markerfacecolor="none",markeredgewidth=.9,label=label)
        for (_,label,color),marker in zip(sets[1:],["o","x"])]
    fig.legend(handles=handles,loc="outside upper center",ncol=3,
        handlelength=1.4,columnspacing=.9)
    p.diagnostics["appendix_recovered_spatial_data"] = {
        "continuous_splitting_source": f"dsa_forward_path_dens_pp_{recovered}.dat",
        "independent_initial_seeds":12800000,"Expanse_job_array":55996956,
        "unnormalized_integrals_code_units":{k:v.tolist() for k,v in masses.items()},
        "all_three_curves_from_numerical_histograms":True,"original_image_used":False,
        "x_limits_au":xlimits,"density_y_limits_au^-1":density_limits,
        "difference_y_limits_percent":difference_limits}
    p.save(fig,"spatial_marginal_original_curves", "Three-curve Appendix B comparison at alpha=1, regenerated entirely from numerical histograms. Continuous splitting uses the 12.8-million-seed Expanse rerun (array 55996956), merged across four independent batches by initial-seed count before normalizing. Weighted direct is a thin magenta line, continuous splitting uses blue open circles, and matched IS uses orange crosses. All panels share x limits, and y limits match within each row; normalization remains over the full [-5,5]-au histogram. Top: normalized spatial density in au^-1. Bottom: absolute fractional differences expressed as percentages (100*d); all three marginals must exceed 1e-3 au^-1. Markers show the saved bin centers without offsets. No original-image pixels, smoothing, tracing, or resampling are used.")


def diagnostics_plot(p):
    selected = [("split_0", 0, "Direct integration", BLACK, "o"),
        ("split_0", 1.5, "Continuous splitting", IBM[0], "s"),
        ("split_1", 1.5, "Discrete splitting", IBM[1], "^"),
        ("imps_0", 1.5, "Splitting-matched IS", IBM[4], "D"),
        ("imps_1", 1, "Shock-directed IS", IBM[3], "v")]
    fig, axs = figure_axes(2, 1, figsize=(7.1, 5.4), sharex=True,
        gridspec_kw={"height_ratios": [1.1, 1.3]}, layout="constrained")
    for stem, param, label, color, marker in selected:
        rec = next(r for r in p.metrics[stem] if r["param"] == param)
        axs[0].plot(p.El, rec["ratio"], color=color, linewidth=1.25)
        idx = np.arange(2, 100, 8)
        axs[0].errorbar(p.El[idx], rec["ratio"][idx], yerr=rec["mean_se"][idx],
            color=color, marker=marker, markersize=3.5, linestyle="none",
            linewidth=.9, capsize=1.5, label=label)
        axs[1].plot(p.El, rec["effective_contributors"], color=color,
            marker=marker, markersize=3.5, markevery=8, linewidth=1.3)
    axs[0].axhline(1, color=BLACK, linewidth=.7, linestyle="--")
    axs[0].set_ylim(.35, 1.55)
    axs[0].set_ylabel(r"$\widehat H/H_{\rm ref}$")
    axs[1].set_yscale("log")
    axs[1].set_ylim(1, 2e7)
    axs[1].yaxis.set_minor_locator(NullLocator())
    axs[1].set_ylabel("Effective seeds")
    axs[1].set_xlabel("$E$ (MeV)")
    for ax in axs:
        log_energy(ax, True)
    handles, labels = axs[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside upper center", ncol=2)
    p.save(fig, "sampling_diagnostics", "Diagnostic only. Continuous splitting, discrete splitting and matched IS use alpha=1.5; shock-directed IS uses A0=1 au^-1. Top: integrated tail estimates relative to the sharp-shock reference, with one-standard-error Monte Carlo bars reconstructed from the saved seed first and second moments and the seed counts in the matching logs. E is the proton kinetic energy corresponding to the lower momentum threshold p. These are uncertainties of the mean H estimate, NOT error bars for e2 or c. Bottom: (sum Y)^2/sum Y^2 for the independent seed contributions Y, including all descendants within each splitting seed. This equals the number of contributing seeds for direct unweighted Bernoulli contributions and is a moment-based effective count for weighted or branching contributions. Substantial reference offsets remain when Monte Carlo error is small; more seeds do not resolve that systematic/reference discrepancy. The saved data do not identify its cause.")


def finish(p):
    writer = PdfWriter()
    order = ["efficiency_gains", "spectrum_comparison", "performance_metrics",
        "path_intensity_available", "expected_amplification", "induced_drifts",
        "spatial_marginal_original_curves", "spatial_marginal_available",
        "one_count_limit", "sampling_diagnostics"]
    for name in order:
        writer.append(PdfReader(p.output/(name+".pdf")))
    with (p.output/"dsa_plots_ibm.pdf").open("wb") as f:
        writer.write(f)
    p.diagnostics["source_sha256"] = {str(f.relative_to(p.data)):
        hashlib.sha256(f.read_bytes()).hexdigest() for f in sorted(p.data.rglob("*")) if f.is_file()}
    (p.output/"diagnostics.json").write_text(json.dumps(p.diagnostics, indent=2)+"\n")
    (p.output/"captions.txt").write_text("\n\n".join(name+"\n"+p.captions[name] for name in order)+"\n")
    print("Created",len(p.figures),"figures and dsa_plots_ibm.pdf")
    print("Continuous path file fitted exponent:", p.diagnostics["continuous_path_file"]["estimated_amplification_exponent"])
    print("Spectrum steps:",[(r["method"],r["total_steps"]) for r in p.diagnostics["spectrum_budgets"]])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data", type=Path, default=Path(__file__).parent/"data")
    ap.add_argument("--output", type=Path, default=Path(__file__).parent/"figures")
    args = ap.parse_args()
    style()
    p = Project(args.data, args.output)
    performance(p)
    spectrum(p)
    one_count(p)
    likelihood(p)
    path_intensity(p)
    marginal(p)
    original_marginal_typography(p)
    diagnostics_plot(p)
    finish(p)


if __name__ == "__main__":
    main()
