"""End-to-end tests of `pipeline batch-effects` on synthetic count files with a
known batch effect."""

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp
import yaml

from rosmap_processing.pipelines.batch_effects import load_spec, run_batch_effects
from rosmap_processing.utils.obs_summary import summarize_obs

N_GENES = 80
BATCH_GENES = range(0, 10)       # up in batch B2
CELLTYPE_GENES = range(10, 30)   # up in astrocytes
DATASET_GENES = range(30, 40)    # up in the second dataset


def _counts(path, rng, donors, batch_of, dataset_shift=False, cells_per_type=60, layer=None):
    """batch_of: donor -> batch (whole donors per batch), or a list of batches
    that each cell is assigned to at random (pooled libraries)."""
    rates = rng.gamma(1.0, 1.0, size=N_GENES)
    blocks, obs = [], []
    for d, donor in enumerate(donors):
        donor_effect = rng.lognormal(0, 0.2, size=N_GENES)
        for celltype in ("Excitatory", "Astrocytes"):
            r = np.tile(rates * donor_effect, (cells_per_type, 1))
            if celltype == "Astrocytes":
                r[:, list(CELLTYPE_GENES)] *= 6
            batch = np.array([batch_of[donor]] * cells_per_type) if isinstance(batch_of, dict) \
                else rng.choice(batch_of, size=cells_per_type)
            r[np.ix_(batch == "B2", list(BATCH_GENES))] *= 5
            if dataset_shift:
                r[:, list(DATASET_GENES)] *= 5
            depth = rng.uniform(0.5, 2.0, size=(cells_per_type, 1))
            blocks.append(rng.poisson(r * depth))
            obs.append(pd.DataFrame({
                "projid": 1000 + d, "individualID": donor, "batch": batch,
                "source_file": f"data/v2/counts/X/{celltype}.h5ad",
                "cogdx": 1 if d % 2 else 4, "msex": d % 3 % 2,
                "age_death": "90+" if d % 4 == 0 else str(75 + d), "pmi": 5.0 + d,
            }, index=[f"{donor}_{celltype}_{i}" for i in range(cells_per_type)]))
    obs = pd.concat(obs)
    var = pd.DataFrame(index=[f"G{i}" for i in range(N_GENES)])
    counts = sp.csr_matrix(np.vstack(blocks).astype(np.float32))
    adata = ad.AnnData(counts, obs=obs, var=var)
    if layer:  # counts in a layer, something else (not counts) in X, as in SeaAD
        adata.layers[layer] = counts
        adata.X = counts * 0.5
    adata.write_h5ad(path)
    return var.index


def _gene_files(tmp_path, names):
    (tmp_path / "with_batch_genes.txt").write_text("\n".join(names[:40]) + "\n")
    (tmp_path / "without_batch_genes.txt").write_text("\n".join(names[10:]) + "\n")
    pd.DataFrame({"var_name": names, "gene_id": [f"ENSG{i:011d}" for i in range(N_GENES)],
                  "ref_name": names, "method": "symbol"}).to_csv(tmp_path / "map.tsv", sep="\t", index=False)
    (tmp_path / "ensembl.txt").write_text("\n".join(f"ENSG{i:011d}.3" for i in range(N_GENES)) + "\n")


CELLTYPES = [
    {"column": "source_file", "match": "Excitatory", "label": "Exc"},
    {"column": "source_file", "match": "Astrocytes", "label": "Ast"},
]
COVARIATES = {
    "diagnosis": {"column": "cogdx", "recode": {"1": "no", "4": "yes"}},
    "sex": {"column": "msex"},
    "age": {"column": "age_death", "numeric": True},
    "pmi": {"column": "pmi", "numeric": True},
}
SETTINGS = dict(seed=1, max_cells_per_donor=80, min_genes=1, min_cells_per_donor=50, n_pcs=10,
                k_neighbors=10, n_permutations=4, min_donors=3, min_cells_per_pseudobulk=5,
                max_cells_silhouette=500, chunk_size=300)


def _write_spec(tmp_path, name, datasets, gene_sets, **extra):
    spec = {"output_dir": str(tmp_path / name), "datasets": datasets, "gene_sets": gene_sets,
            **SETTINGS, **extra}
    path = tmp_path / f"{name}.yaml"
    path.write_text(yaml.safe_dump(spec))
    return path


def _all_output_text(out):
    return "".join(p.read_text() for p in out.iterdir() if p.is_file())


def test_within_dataset_finds_the_batch_effect(tmp_path):
    rng = np.random.default_rng(0)
    donors = [f"R{100 + i}" for i in range(12)]
    batch_of = {d: ["B1", "B2", "B3"][i // 4] for i, d in enumerate(donors)}
    names = _counts(tmp_path / "counts.h5ad", rng, donors, batch_of)
    _gene_files(tmp_path, list(names))
    spec = _write_spec(tmp_path, "within", {"one": {
        "counts": [str(tmp_path / "missing.h5ad"), str(tmp_path / "counts.h5ad")],
        "donor": "projid", "celltype": CELLTYPES, "batch": ["batch"], "covariates": COVARIATES,
    }}, {"with batch genes": str(tmp_path / "with_batch_genes.txt"),
         "without batch genes": str(tmp_path / "without_batch_genes.txt")})

    out = run_batch_effects(spec)
    cells = pd.read_csv(out / "cells.tsv", sep="\t")
    assert set(cells["celltype"]) == {"all", "Exc", "Ast"}
    assert cells.loc[cells["celltype"] == "all", "cells_sampled"].item() == 12 * 80

    m = pd.read_csv(out / "cell_metrics.tsv", sep="\t")
    pcr = m[(m["label"] == "batch:batch") & (m["metric"] == "pc_regression") & (m["celltype"] == "mean")]
    excess = pcr.set_index("gene_set")["excess"]
    assert excess["with batch genes"] > 0.1
    assert excess["with batch genes"] > excess["without batch genes"] + 0.1
    mixing = m[(m["label"] == "batch:batch") & (m["metric"] == "mixing_ratio") & (m["celltype"] == "Exc")]
    mixing = mixing.set_index("gene_set")["value"]
    assert mixing["with batch genes"] < mixing["without batch genes"]
    celltype_pcr = m[(m["label"] == "celltype") & (m["celltype"] == "all") & (m["metric"] == "pc_regression")]
    assert (celltype_pcr["value"] > 0.2).all()

    d = pd.read_csv(out / "donor_metrics.tsv", sep="\t")
    assert {"partial_r2_median", "r2_alone_median", "balanced_accuracy"} <= set(d["metric"])
    conf = pd.read_csv(out / "confounding.tsv", sep="\t")
    assert set(conf["covariate"]) == {"diagnosis", "sex", "age", "pmi"}
    assert (pd.read_csv(out / "summary.tsv", sep="\t")["gene_set"].nunique()) == 2

    # Aggregates only: no donor ID appears in any output
    text = _all_output_text(out)
    assert not any(donor in text for donor in donors)
    assert not any(str(1000 + i) in text.split() for i in range(12))

    with pytest.raises(FileExistsError):
        run_batch_effects(spec)


def test_matched_donors_across_datasets(tmp_path):
    rng = np.random.default_rng(1)
    donors = [f"R{200 + i}" for i in range(8)]
    batch_of = {d: "B1" for d in donors}
    names = _counts(tmp_path / "a.h5ad", rng, donors, batch_of)
    # The second dataset has two extra donors that must be left out, and a shift
    _counts(tmp_path / "b.h5ad", rng, donors + ["R998", "R999"],
            {**batch_of, "R998": "B1", "R999": "B1"}, dataset_shift=True)
    _gene_files(tmp_path, list(names))
    dataset = lambda f: {"counts": str(tmp_path / f), "gene_map": str(tmp_path / "map.tsv"),
                         "donor": "projid", "match_on": "individualID", "celltype": CELLTYPES,
                         "covariates": COVARIATES}
    spec = _write_spec(tmp_path, "matched", {"a": dataset("a.h5ad"), "b": dataset("b.h5ad")},
                       {"ensembl": str(tmp_path / "ensembl.txt")}, match_donors=True)

    out = run_batch_effects(spec)
    cells = pd.read_csv(out / "cells.tsv", sep="\t").set_index(["dataset", "celltype"])
    assert cells.loc[("b", "all"), "donors_sampled"] == 8
    assert cells.loc[("a", "Exc"), "cells_balanced"] == cells.loc[("b", "Exc"), "cells_balanced"] > 0
    assert cells.loc[("a", "Exc"), "donors_balanced"] == cells.loc[("b", "Exc"), "donors_balanced"]
    m = pd.read_csv(out / "cell_metrics.tsv", sep="\t")
    pcr = m[(m["label"] == "dataset") & (m["metric"] == "pc_regression") & (m["celltype"] == "mean")]
    assert pcr["excess"].item() > 0.05
    d = pd.read_csv(out / "donor_metrics.tsv", sep="\t")
    paired = d[(d["label"] == "a vs b") & (d["celltype"] == "Exc")].set_index("metric")["value"]
    assert {"share_dataset", "share_donor", "share_residual", "donor_correlation"} <= set(paired.index)
    assert paired["share_dataset"] > 0
    assert not any(donor in _all_output_text(out) for donor in donors)


def test_different_donors_across_datasets(tmp_path):
    rng = np.random.default_rng(3)
    donors_a = [f"R{400 + i}" for i in range(8)]
    donors_b = [f"H{500 + i}" for i in range(5)]  # fewer: balanced for the cell-level metrics
    names = _counts(tmp_path / "a.h5ad", rng, donors_a, {d: "B1" for d in donors_a})
    _counts(tmp_path / "b.h5ad", rng, donors_b, {d: "B1" for d in donors_b}, dataset_shift=True,
            layer="UMIs")
    _gene_files(tmp_path, list(names))
    dataset = lambda f, **kw: {"counts": str(tmp_path / f), "gene_map": str(tmp_path / "map.tsv"),
                               "donor": "individualID", "celltype": CELLTYPES,
                               "covariates": COVARIATES, **kw}
    spec = _write_spec(tmp_path, "across", {"a": dataset("a.h5ad"), "b": dataset("b.h5ad", layer="UMIs")},
                       {"ensembl": str(tmp_path / "ensembl.txt")})

    out = run_batch_effects(spec)
    m = pd.read_csv(out / "cell_metrics.tsv", sep="\t")
    pcr = m[(m["label"] == "dataset") & (m["metric"] == "pc_regression") & (m["celltype"] == "mean")]
    assert pcr["excess"].item() > 0.05
    cells = pd.read_csv(out / "cells.tsv", sep="\t").set_index(["dataset", "celltype"])
    for t in ("Exc", "Ast"):
        a, b = cells.loc[("a", t)], cells.loc[("b", t)]
        assert a["cells_balanced"] == b["cells_balanced"] > 0
        assert a["donors_balanced"] == b["donors_balanced"] <= 5
        assert a["cells_balanced"] % a["donors_balanced"] == 0
    exc = m[(m["label"] == "dataset") & (m["metric"] == "pc_regression") & (m["celltype"] == "Exc")]
    assert exc["n_cells"].item() == 2 * cells.loc[("a", "Exc"), "cells_balanced"]
    d = pd.read_csv(out / "donor_metrics.tsv", sep="\t")
    accuracy = d[(d["label"] == "dataset") & (d["metric"] == "balanced_accuracy") & (d["celltype"] == "mean")]
    assert accuracy["value"].item() > 0.9
    conf = pd.read_csv(out / "confounding.tsv", sep="\t")
    assert set(conf["grouping"]) == {"dataset"}
    assert not any(donor in _all_output_text(out) for donor in donors_a + donors_b)


def test_counts_layer_must_hold_counts(tmp_path):
    rng = np.random.default_rng(4)
    donors = [f"R{600 + i}" for i in range(4)]
    names = _counts(tmp_path / "c.h5ad", rng, donors, {d: "B1" for d in donors}, layer="UMIs")
    _gene_files(tmp_path, list(names))
    spec = _write_spec(tmp_path, "x_not_counts", {"c": {
        "counts": str(tmp_path / "c.h5ad"), "donor": "projid", "celltype": CELLTYPES}},
        {"g": str(tmp_path / "with_batch_genes.txt")}, min_donors=2)
    with pytest.raises(ValueError, match="not counts"):
        run_batch_effects(spec)  # X holds counts * 0.5, the counts are in layer UMIs


def test_batch_that_splits_donors(tmp_path):
    rng = np.random.default_rng(5)
    donors = [f"R{700 + i}" for i in range(8)]
    names = _counts(tmp_path / "pooled.h5ad", rng, donors, ["B1", "B2", "B3"])
    _gene_files(tmp_path, list(names))
    spec = _write_spec(tmp_path, "pooled", {"one": {
        "counts": str(tmp_path / "pooled.h5ad"), "donor": "projid", "celltype": CELLTYPES,
        "batch": ["batch"], "covariates": COVARIATES,
    }}, {"with batch genes": str(tmp_path / "with_batch_genes.txt"),
         "without batch genes": str(tmp_path / "without_batch_genes.txt")})

    out = run_batch_effects(spec)
    m = pd.read_csv(out / "cell_metrics.tsv", sep="\t")
    pcr = m[(m["label"] == "batch:batch") & (m["metric"] == "pc_regression") & (m["celltype"] == "Exc")]
    assert set(pcr["null_type"]) == {"within"}
    excess = pcr.set_index("gene_set")["excess"]
    assert excess["with batch genes"] > 0.1
    assert abs(excess["without batch genes"]) < 0.05
    # No donor-level value, so no donor-level batch rows
    d = pd.read_csv(out / "donor_metrics.tsv", sep="\t")
    assert "batch:batch" not in set(d["label"])
    assert "batch:batch" not in set(pd.read_csv(out / "confounding.tsv", sep="\t")["grouping"])


def _legacy_h5ad(path, n=40):
    """obs in the anndata < 0.8 layout, as in the SeaAD download: categoricals
    as codes referring to their categories, names with '/' as nested groups."""
    import h5py

    columns = {"Donor ID": (np.repeat([0, 1], n // 2), ["H1", "H2"]),
               "Sex": (np.tile([0, 1], n // 2), ["Female", "Male"]),
               "Hispanic/Latino": (np.r_[np.full(n - 1, 0), -1], ["No", "Yes"])}
    with h5py.File(path, "w") as f:
        obs = f.create_group("obs")
        obs.attrs.update({"_index": "cell", "encoding-type": "dataframe", "encoding-version": "0.1.0",
                          "column-order": np.array([*columns, "Age at Death"], dtype=h5py.string_dtype())})
        obs.create_dataset("cell", data=np.array([f"c{i}" for i in range(n)], dtype=h5py.string_dtype()))
        for name, (codes, categories) in columns.items():
            cats = obs.create_dataset(f"__categories/{name}", data=np.array(categories, dtype=h5py.string_dtype()))
            obs.create_dataset(name, data=codes.astype(np.int8)).attrs["categories"] = cats.ref
        obs.create_dataset("Age at Death", data=np.linspace(70, 95, n))


def test_reads_the_old_anndata_layout(tmp_path):
    import h5py

    from rosmap_processing.utils.obs_summary import obs_column_names, read_obs_column

    _legacy_h5ad(tmp_path / "old.h5ad")
    with h5py.File(tmp_path / "old.h5ad", "r") as f:
        assert obs_column_names(f) == ["Donor ID", "Sex", "Hispanic/Latino", "Age at Death"]
        sex = read_obs_column(f, "Sex")
        assert list(sex.categories) == ["Female", "Male"] and sex[1] == "Male"
        hispanic = pd.Series(read_obs_column(f, "Hispanic/Latino"))
        assert hispanic.isna().sum() == 1
    table = summarize_obs(tmp_path / "old.h5ad", donor_column="Donor ID").set_index("column")
    assert table.loc["Sex", "level"] == "cell" and table.loc["Sex", "n_unique"] == 2
    assert table.loc["Hispanic/Latino", "pct_missing"] == 2.5


def test_spec_rejects_unknown_keys(tmp_path):
    path = tmp_path / "bad.yaml"
    path.write_text(yaml.safe_dump({"output_dir": "x", "datasets": {"a": {}}, "gene_sets": {"g": "x"},
                                    "n_permutation": 3}))
    with pytest.raises(ValueError, match="Unknown key"):
        load_spec(path)


def test_obs_summary_has_no_values(tmp_path):
    rng = np.random.default_rng(2)
    donors = [f"R{300 + i}" for i in range(4)]
    _counts(tmp_path / "c.h5ad", rng, donors, {d: "B1" if i < 2 else "B2" for i, d in enumerate(donors)},
            cells_per_type=5)
    table = summarize_obs(tmp_path / "c.h5ad").set_index("column")
    assert table.loc["batch", "level"] == "donor"
    assert table.loc["source_file", "level"] == "cell"
    assert table.loc["batch", "n_unique"] == 2
    text = table.to_csv()
    assert not any(d in text for d in donors) and "B1" not in text and "Excitatory" not in text
