class AtlasQueryResult:
    def __init__(
        self,
        predicted_module,
        module_scores,
        all_results,
        module_consistent_results,
        query_genes
    ):
        self.predicted_module = predicted_module
        self.module_scores = module_scores
        self.all_results = all_results
        self.module_consistent_results = module_consistent_results
        self.query_genes = list(query_genes)

    def module_margin(self):
        scores = sorted(
            self.module_scores.values(),
            reverse=True
        )

        if len(scores) < 2:
            return None

        return scores[0] - scores[1]

    def top_module_agreement(self, top_n=5):
        if self.all_results.empty:
            return 0

        if "predicted_module" not in self.all_results.columns:
            return 0

        top = self.all_results.head(top_n)

        return (
            top["predicted_module"] == self.predicted_module
        ).mean()

    def mean_top_recall(self, top_n=5):
        if self.all_results.empty:
            return 0

        if "query_recall" not in self.all_results.columns:
            return 0

        return self.all_results.head(top_n)["query_recall"].mean()

    def evidence(self, top_n=5):
        margin = self.module_margin()
        agreement = self.top_module_agreement(top_n=top_n)
        mean_recall = self.mean_top_recall(top_n=top_n)

        return {
            "module_margin": margin,
            "top_module_agreement": agreement,
            "mean_top_recall": mean_recall,
            "query_size": len(self.query_genes),
            "n_retrieved": len(self.all_results)
        }

    def confidence_score(self, top_n=5):
        ev = self.evidence(top_n=top_n)

        margin = ev["module_margin"] or 0
        agreement = ev["top_module_agreement"]
        mean_recall = ev["mean_top_recall"]

        margin_scaled = min(margin / 0.05, 1)

        score = (
            0.4 * margin_scaled
            + 0.4 * agreement
            + 0.2 * mean_recall
        )

        return score

    def evidence_score(self, top_n=5):
        """
        Alias for confidence_score().

        Validation benchmarks showed that this score behaves as an empirical
        evidence score. This method gives the API a clearer name while keeping
        confidence_score() for backward compatibility.
        """
        return self.confidence_score(top_n=top_n)

    def confidence_label(self, top_n=5):
        score = self.confidence_score(top_n=top_n)

        if score >= 0.85:
            return "Very High"
        elif score >= 0.70:
            return "High"
        elif score >= 0.50:
            return "Moderate"
        else:
            return "Low"

    def confidence_stars(self, top_n=5):
        score = self.confidence_score(top_n=top_n)

        n_stars = round(score * 5)
        n_stars = max(0, min(5, n_stars))

        return "★" * n_stars + "☆" * (5 - n_stars)

    def summary(self, top_n=5):
        print("=" * 50)
        print("PerturbAtlas-K562-v1.0 Query Result")
        print("=" * 50)

        print("\nPredicted module:")
        print(f"  {self.predicted_module}")

        print("\nModule scores:")
        sorted_scores = sorted(
            self.module_scores.items(),
            key=lambda x: x[1],
            reverse=True
        )

        for module, score in sorted_scores:
            print(f"  {module:15s} {score:.4f}")

        print("\nEvidence score:")
        print(
            f"  {self.confidence_stars(top_n)} "
            f"{self.confidence_label(top_n)} "
            f"({self.evidence_score(top_n):.2f})"
        )

        print("\nTop perturbations:")
        if self.module_consistent_results.empty:
            print("  No module-consistent perturbations found.")
        else:
            top = self.module_consistent_results.head(top_n)

            for i, row in enumerate(top.itertuples(), start=1):
                print(
                    f"  {i}. {row.perturbation:10s} "
                    f"| recall={row.query_recall:.2f} "
                    f"| precision={row.precision:.2f} "
                    f"| module={row.predicted_module}"
                )

        print("\nQuery size:")
        print(f"  {len(self.query_genes)} genes")

    def top_perturbations(self, n=10, module_consistent=True):
        if module_consistent:
            return self.module_consistent_results.head(n)
        else:
            return self.all_results.head(n)

    def to_csv(self, path, module_consistent=True):
        if module_consistent:
            self.module_consistent_results.to_csv(path, index=False)
        else:
            self.all_results.to_csv(path, index=False)

    def interpret(self, top_n=5, min_confidence=0.50):
        module_text = {
            "myeloid": "myeloid differentiation / innate immune-like program",
            "erythroid": "erythroid / hemoglobin-associated program",
            "developmental": "developmental or cell-state transition program"
        }

        ev = self.evidence(top_n=top_n)
        score = self.evidence_score(top_n=top_n)
        label = self.confidence_label(top_n=top_n)

        print("=" * 50)
        print("PerturbAtlas-K562-v1.0 Interpretation")
        print("=" * 50)

        print("\nPredicted biological program:")
        print(
            f"  {self.predicted_module} "
            f"({module_text.get(self.predicted_module, 'unannotated program')})"
        )

        print("\nEvidence score:")
        print(
            f"  {self.confidence_stars(top_n)} "
            f"{label} "
            f"({score:.2f})"
        )

        print("\nEvidence:")
        margin = ev["module_margin"]

        if margin is None:
            print("  Module score margin       : NA")
        else:
            print(f"  Module score margin       : {margin:.4f}")

        print(f"  Top-{top_n} module agreement : {ev['top_module_agreement']:.2f}")
        print(f"  Mean top-{top_n} recall      : {ev['mean_top_recall']:.2f}")
        print(f"  Retrieved perturbations   : {ev['n_retrieved']}")
        print(f"  Query genes               : {ev['query_size']}")

        print("\nMost similar perturbations:")
        if self.module_consistent_results.empty:
            print("  No module-consistent perturbations found.")
        else:
            top = self.module_consistent_results.head(top_n)

            for i, row in enumerate(top.itertuples(), start=1):
                print(
                    f"  {i}. {row.perturbation} "
                    f"(recall={row.query_recall:.2f}, "
                    f"precision={row.precision:.2f})"
                )

        print("\nBiological interpretation:")

        if score < min_confidence:
            print(
                "  Insufficient evidence for a confident biological interpretation. "
                "The query produced a low evidence score, weak module separation, "
                "limited retrieval support, or inconsistent perturbation matches. "
                "This result should be treated as exploratory rather than as a "
                "reliable perturbational program assignment."
            )
            return

        if self.predicted_module == "myeloid":
            print(
                "  The query resembles the K562 myeloid-like perturbation "
                "program. Strong matches to CEBP/SPI1-axis perturbations "
                "suggest activation of a myeloid differentiation or innate "
                "immune-like state."
            )

        elif self.predicted_module == "erythroid":
            print(
                "  The query resembles the K562 erythroid-like perturbation "
                "program. Overlap with erythroid-associated perturbations "
                "suggests activation of hemoglobinization, erythroid "
                "maturation, or red-cell-like transcriptional features."
            )

        elif self.predicted_module == "developmental":
            print(
                "  The query resembles a broader developmental or cell-state "
                "transition program. This module is less lineage-specific than "
                "the myeloid and erythroid programs and should be interpreted "
                "as a regulatory/state-transition signal rather than a single "
                "well-defined lineage."
            )

        else:
            print(
                "  No curated biological interpretation is available for this module."
            )
