import pickle
import pandas as pd

from pathlib import Path
from datetime import datetime

from .result import AtlasQueryResult


class K562Atlas:
    """
    PerturbAtlas implementation for the K562 PerturbAtlas-v1.0 release.

    The atlas stores perturbation-associated gene sets, perturbation metadata,
    module anchors, and lightweight metadata describing the reference atlas.
    """

    VERSION = "1.0.0"

    def __init__(
        self,
        strong_up,
        perturbation_metadata,
        module_anchors=None,
        version="v1.0",
        cell_line="K562",
        atlas_name="PerturbAtlas-K562-v1.0",
        metadata=None
    ):
        self.version = version
        self.cell_line = cell_line
        self.atlas_name = atlas_name

        self.strong_up = self._normalize_strong_up(strong_up)
        self.perturbation_metadata = perturbation_metadata.copy()

        if "perturbation" in self.perturbation_metadata.columns:
            self.perturbation_metadata["perturbation"] = (
                self.perturbation_metadata["perturbation"]
                .astype(str)
                .str.strip()
                .str.upper()
            )

        if module_anchors is None:
            module_anchors = {
                "myeloid": ["CEBPA", "CEBPB", "CEBPE", "SPI1"],
                "erythroid": ["KLF1", "IKZF3", "DUSP9", "PRDM1"],
                "developmental": ["HOXA13", "HOXC13", "TP73", "ZBTB10"]
            }

        self.module_anchors = {
            str(module).strip(): [
                str(anchor).strip().upper()
                for anchor in anchors
                if pd.notna(anchor) and str(anchor).strip() != ""
            ]
            for module, anchors in module_anchors.items()
        }

        self.reference_genes = set().union(*self.strong_up.values())

        self.metadata = metadata or {}
        self.metadata.setdefault("atlas_name", self.atlas_name)
        self.metadata.setdefault("version", self.version)
        self.metadata.setdefault("package_version", self.VERSION)
        self.metadata.setdefault("created_at", datetime.now().isoformat())
        self.metadata.setdefault("organism", "Homo sapiens")
        self.metadata.setdefault("cell_line", self.cell_line)
        self.metadata.setdefault("dataset", "K562 Perturb-seq")
        self.metadata.setdefault("atlas_type", "Gene-set perturbation atlas")
        self.metadata.setdefault("n_reference_genes", len(self.reference_genes))
        self.metadata.setdefault("n_perturbations", len(self.strong_up))

    @staticmethod
    def _normalize_strong_up(strong_up):
        """
        Normalize perturbation names and gene symbols to uppercase strings.
        Also removes duplicate and empty gene names.
        """
        if strong_up is None:
            raise ValueError("strong_up cannot be None.")

        if len(strong_up) == 0:
            raise ValueError("strong_up cannot be empty.")

        normalized = {}

        for pert, genes in strong_up.items():
            pert_name = str(pert).strip().upper()

            if pert_name == "":
                continue

            gene_set = {
                str(gene).strip().upper()
                for gene in genes
                if pd.notna(gene) and str(gene).strip() != ""
            }

            normalized[pert_name] = gene_set

        if len(normalized) == 0:
            raise ValueError("strong_up contains no valid perturbation signatures.")

        return normalized

    def _sync_metadata(self):
        """
        Keep legacy attributes and metadata dictionary consistent.

        This is useful after loading older pickle files that may not contain
        all metadata fields introduced during the engineering pass.
        """
        if not hasattr(self, "version"):
            self.version = "v1.0"

        if not hasattr(self, "cell_line"):
            self.cell_line = "K562"

        if not hasattr(self, "atlas_name"):
            self.atlas_name = "PerturbAtlas-K562-v1.0"

        if not hasattr(self, "metadata") or self.metadata is None:
            self.metadata = {}

        self.metadata.setdefault("atlas_name", self.atlas_name)
        self.metadata.setdefault("version", self.version)
        self.metadata.setdefault("package_version", self.VERSION)
        self.metadata.setdefault("created_at", "unknown")
        self.metadata.setdefault("organism", "Homo sapiens")
        self.metadata.setdefault("cell_line", self.cell_line)
        self.metadata.setdefault("dataset", "K562 Perturb-seq")
        self.metadata.setdefault("atlas_type", "Gene-set perturbation atlas")
        self.metadata["n_reference_genes"] = len(self.reference_genes)
        self.metadata["n_perturbations"] = len(self.strong_up)

    def info(self):
        """
        Print a human-readable summary of the atlas.
        """
        self._sync_metadata()

        print("=" * 50)
        print(self.metadata.get("atlas_name", self.atlas_name))
        print("=" * 50)

        print("\nAtlas type:")
        print(f"  {self.metadata.get('atlas_type', 'Gene-set perturbation atlas')}")

        print("\nReference:")
        print(f"  Cell line : {self.metadata.get('cell_line', self.cell_line)}")
        print(f"  Organism  : {self.metadata.get('organism', 'Homo sapiens')}")
        print(f"  Dataset   : {self.metadata.get('dataset', 'K562 Perturb-seq')}")
        print(f"  Version   : {self.metadata.get('version', self.version)}")

        created_at = self.metadata.get("created_at", None)
        if created_at not in (None, "", "unknown"):
            print(f"  Created   : {created_at}")

        print("\nContent:")
        print(f"  Perturbations  : {len(self.strong_up)}")
        print(f"  Reference genes: {len(self.reference_genes)}")
        print(f"  Modules        : {', '.join(self.module_anchors.keys())}")

        print("\nModule anchors:")
        for module, anchors in self.module_anchors.items():
            print(f"  {module:15s}: {', '.join(anchors)}")

    def summary_dict(self):
        """
        Return atlas metadata as a dictionary for programmatic use.
        """
        self._sync_metadata()

        return {
            "atlas_name": self.metadata.get("atlas_name", self.atlas_name),
            "version": self.metadata.get("version", self.version),
            "package_version": self.metadata.get("package_version", self.VERSION),
            "created_at": self.metadata.get("created_at"),
            "organism": self.metadata.get("organism", "Homo sapiens"),
            "cell_line": self.metadata.get("cell_line", self.cell_line),
            "dataset": self.metadata.get("dataset", "K562 Perturb-seq"),
            "n_perturbations": len(self.strong_up),
            "n_reference_genes": len(self.reference_genes),
            "n_modules": len(self.module_anchors),
            "modules": list(self.module_anchors.keys()),
            "metadata_columns": list(self.perturbation_metadata.columns)
        }

    def __repr__(self):
        return (
            f"K562Atlas("
            f"name='{getattr(self, 'atlas_name', 'PerturbAtlas-K562-v1.0')}', "
            f"perturbations={len(self.strong_up)}, "
            f"reference_genes={len(self.reference_genes)}, "
            f"modules={len(self.module_anchors)}"
            f")"
        )

    def __len__(self):
        return len(self.strong_up)

    def perturbations(self):
        """
        Return the available perturbation names in the atlas.
        """
        return sorted(self.strong_up.keys())

    def modules(self):
        """
        Return the available biological modules in the atlas.
        """
        return list(self.module_anchors.keys())

    def save(self, path):
        """
        Save the atlas object to disk using pickle.

        Example:
            atlas.save("PerturbAtlas-K562-v1.0.pkl")
        """
        self._sync_metadata()

        path = Path(path)

        if path.parent != Path("."):
            path.parent.mkdir(parents=True, exist_ok=True)

        with open(path, "wb") as f:
            pickle.dump(self, f)

        return str(path)

    @classmethod
    def load(cls, path):
        """
        Load a saved K562Atlas object from disk.

        Example:
            atlas = K562Atlas.load("PerturbAtlas-K562-v1.0.pkl")
        """
        path = Path(path)

        if not path.exists():
            raise FileNotFoundError(f"Atlas file not found: {path}")

        with open(path, "rb") as f:
            atlas = pickle.load(f)

        if not isinstance(atlas, cls):
            raise TypeError(
                f"Loaded object is not a {cls.__name__}. "
                f"Found: {type(atlas).__name__}"
            )

        # Upgrade older saved objects to the current engineering-pass structure.
        atlas.strong_up = cls._normalize_strong_up(atlas.strong_up)

        if "perturbation" in atlas.perturbation_metadata.columns:
            atlas.perturbation_metadata["perturbation"] = (
                atlas.perturbation_metadata["perturbation"]
                .astype(str)
                .str.strip()
                .str.upper()
            )

        atlas.module_anchors = {
            str(module).strip(): [
                str(anchor).strip().upper()
                for anchor in anchors
                if pd.notna(anchor) and str(anchor).strip() != ""
            ]
            for module, anchors in atlas.module_anchors.items()
        }

        atlas.reference_genes = set().union(*atlas.strong_up.values())
        atlas._sync_metadata()

        return atlas

    def _normalize_query_genes(self, query_genes):
        if query_genes is None:
            raise ValueError("query_genes cannot be None.")

        if isinstance(query_genes, str):
            raise TypeError(
                "query_genes must be a list, set, tuple, or similar iterable "
                "of gene names, not a single string."
            )

        query_genes = list(query_genes)

        if len(query_genes) == 0:
            raise ValueError("query_genes cannot be empty.")

        normalized = {
            str(gene).strip().upper()
            for gene in query_genes
            if pd.notna(gene) and str(gene).strip() != ""
        }

        if len(normalized) == 0:
            raise ValueError(
                "query_genes contains no valid gene names after cleaning."
            )

        return normalized

    def _jaccard(self, genes_a, genes_b):
        genes_a = set(genes_a)
        genes_b = set(genes_b)

        union = genes_a | genes_b

        if len(union) == 0:
            return 0

        return len(genes_a & genes_b) / len(union)

    def score_modules(self, query_genes):
        query_genes = self._normalize_query_genes(query_genes)
        query_genes = query_genes & self.reference_genes

        if len(query_genes) == 0:
            raise ValueError(
                "None of the query genes are present in the atlas after normalization."
            )

        module_scores = {}

        for module, anchors in self.module_anchors.items():
            scores = []

            for anchor in anchors:
                if anchor not in self.strong_up:
                    continue

                scores.append(
                    self._jaccard(query_genes, self.strong_up[anchor])
                )

            module_scores[module] = (
                sum(scores) / len(scores)
                if len(scores) > 0
                else 0
            )

        return module_scores

    def retrieve(
        self,
        query_genes,
        top_n=15,
        min_recall=0.7,
        min_signature_size=20
    ):
        query_genes = self._normalize_query_genes(query_genes)
        query_genes = query_genes & self.reference_genes

        if len(query_genes) == 0:
            raise ValueError(
                "None of the query genes are present in the atlas after normalization."
            )

        rows = []

        for pert, genes in self.strong_up.items():
            genes = set(genes)

            if len(genes) < min_signature_size:
                continue

            overlap = len(query_genes & genes)
            recall = overlap / len(query_genes)
            precision = overlap / len(genes) if len(genes) > 0 else 0

            if recall < min_recall:
                continue

            rows.append({
                "perturbation": pert,
                "query_recall": recall,
                "precision": precision,
                "shared_genes": overlap,
                "signature_size": len(genes)
            })

        result = pd.DataFrame(rows)

        if result.empty:
            return result

        result = result.merge(
            self.perturbation_metadata,
            on="perturbation",
            how="left"
        )

        sort_columns = [
            "query_recall",
            "shared_genes"
        ]

        if "master_score_norm" in result.columns:
            sort_columns.append("master_score_norm")

        result = result.sort_values(
            sort_columns,
            ascending=False
        )

        return result.head(top_n)

    def query(
        self,
        query_genes,
        top_n=15,
        min_recall=0.7,
        min_signature_size=20,
        module_filter=True
    ):
        clean_query_genes = self._normalize_query_genes(query_genes)

        module_scores = self.score_modules(clean_query_genes)

        predicted_module = max(
            module_scores,
            key=module_scores.get
        )

        results = self.retrieve(
            query_genes=clean_query_genes,
            top_n=top_n,
            min_recall=min_recall,
            min_signature_size=min_signature_size
        )

        if module_filter and not results.empty:
            if "predicted_module" not in results.columns:
                raise KeyError(
                    "The perturbation metadata must contain a 'predicted_module' "
                    "column when module_filter=True."
                )

            module_consistent_results = results[
                results["predicted_module"] == predicted_module
            ].copy()
        else:
            module_consistent_results = results

        return AtlasQueryResult(
            predicted_module=predicted_module,
            module_scores=module_scores,
            all_results=results,
            module_consistent_results=module_consistent_results,
            query_genes=clean_query_genes
        )
